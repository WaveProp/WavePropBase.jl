@info "including gmshIO.jl from WavePropBase/IO"

using .Gmsh

"""
macro gmsh(ex)

Use a `try` block to initialize `gmsh`, execute `ex`, and finalize `gmsh`
regardless of how `ex` exits.
"""
macro gmsh(ex)
    return quote
        gmsh.initialize()
        try
            set_verbosity(0)
            $(esc(ex))
        finally
            # make sure we finalize gmsh if something goes wrong
            gmsh.finalize()
        end
    end
end

struct GmshEntity <: AbstractEntity
    dim::UInt8
    gmshtag::Int64
    tag::Int64
    boundary::Vector{GmshEntity}
    model::String
    function GmshEntity(d::Integer, gmshtag::Integer, model, boundary=GmshEntity[])
        msg = "an elementary entities in the boundary has the wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d-1 msg
        end
        tag = new_tag(d)
        ent = new(d, gmshtag, tag, boundary, model)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

gmshtag(e::GmshEntity) = e.gmshtag

"""
    import_domain([model;dim=3])

Construct a [`Domain`](@ref) from the `gmsh` `model` with all entities of
dimension `dim`; by defaul the current `gmsh` model is used.

!!! note
    This function assumes that `gmsh` has been initialized, and
    does not handle its finalization.
"""
function gmsh_import_domain(model=gmsh.model.getCurrent();dim=3)
    Ω = Domain() # Create empty domain
    gmsh_import_domain!(Ω,model;dim)
    return Ω
end

"""
    gmsh_import_domain!(Ω::Domain,[model;dim=3])

Like [`gmsh_import_domain`](@ref), but appends entities to `Ω` instead of
creating a new domain.

!!! note
    This function assumes that `gmsh` has been initialized, and does not handle its
    finalization.
"""
function gmsh_import_domain!(Ω::Domain,model=gmsh.model.getCurrent();dim=3)
    old_model = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model)
    dim_tags = gmsh.model.getEntities(dim)
    for (_, tag) in dim_tags
        ent = GmshEntity(dim, tag, model)
        _fill_entity_boundary!(ent,model)
        push!(Ω, ent)
    end
    gmsh.model.setCurrent(old_model)
    return Ω
end

"""
    gmsh_import_mesh(Ω;[dim=3])

Create a `GenericMesh` for the entities in `Ω`. Passing `dim=2` will create a
two-dimensional mesh by projecting the original mesh onto the `x,y` plane.

!!! danger
    This function assumes that `gmsh` has been initialized, and does not handle its
    finalization.
"""
function gmsh_import_mesh(Ω::Domain;dim=3)
    msh = GenericMesh{3,Float64}()
    gmsh_import_mesh!(msh,Ω)
    if dim == 3
        return msh
    elseif dim == 2
        return convert_to_2d(msh)
    else
        error("`dim` value must be `2` or `3`")
    end
end

"""
    meshgen!(msh,Ω)

Similar to [`meshgen`](@ref), but append information to `msh` instead of
creating a new mesh.

!!! danger
    This function assumes that `gmsh` has been initialized, and does not handle its
    finalization.
"""
function gmsh_import_mesh!(msh::GenericMesh,Ω::Domain)
    _, coord, _ = gmsh.model.mesh.getNodes()
    nodes = reinterpret(SVector{3,Float64}, coord) |> collect
    shift = length(msh.nodes) # gmsh node tags need to be shifted
    append!(msh.nodes,nodes)
    els = elements(msh)
    e2t = ent2tags(msh)
    # Recursively populate the dictionaries
    _domain_to_mesh!(els, e2t, Ω, shift)
    return msh
end

"""
    gmsh_read_geo(fname::String;dim=3)

Read a `.geo` file and generate a [`Domain`](@ref) with all entities of
dimension `dim`.

!!! danger
    This function assumes that `gmsh` has been initialized, and does not handle its
    finalization.
"""
function gmsh_read_geo(fname;dim=3)
    Ω = Domain() # Create empty domain
    try
        gmsh.open(fname)
    catch
        @error "could not open $fname"
    end
    gmsh_import_domain!(Ω;dim)
    return Ω
end

"""
    gmsh_read_msh(fname::String;dim=3)

Read `fname` and create a `Domain` and a `GenericMesh` structure with all
entities in `Ω` of dimension `dim`.

!!! danger
    This function assumes that `gmsh` has been initialized, and does not handle its
    finalization.
"""
function gmsh_read_msh(fname;dim=3)
    Ω = Domain()
    try
        gmsh.open(fname)
    catch
        @error "could not open $fname"
    end
    Ω   = gmsh_import_domain(;dim)
    msh = gmsh_import_mesh(Ω;dim)
    return Ω,msh
end

"""
    _fill_entity_boundary!

Use the `gmsh` API to add the boundary of an `ElementaryEntity`.

This is a helper function, and should not be called by itself.
"""
function _fill_entity_boundary!(ent,model)
    combine  = true # FIXME: what should we use here?
    oriented = false
    dim_tags = gmsh.model.getBoundary((geometric_dimension(ent), gmshtag(ent)),combine,oriented)
    for (d, t) in dim_tags
        # if haskey(ENTITIES,(d,t))
        #     bnd = ENTITIES[(d,t)]
        # else
            bnd = GmshEntity(d, t, model)
            _fill_entity_boundary!(bnd,model)
        # end
        push!(ent.boundary, bnd)
    end
    return ent
end

"""
    _domain_to_mesh!(elements, ent2tag, Ω::Domain)

Recursively populating the dictionaries `elements` and `ent2tag`.
"""
function _domain_to_mesh!(elements, ent2tag, Ω::Domain, shift)
    isempty(Ω) && (return elements, ent2tag)
    for ω in Ω
        _ent_to_mesh!(elements, ent2tag, ω, shift)
    end
    Γ = skeleton(Ω)
    _domain_to_mesh!(elements, ent2tag, Γ, shift)
end

"""
    _ent_to_mesh!(elements, ent2tag, ω::ElementaryEntity)

For each element type used to mesh `ω`:
- push into `elements::Dict` the pair `etype=>ntags`;
- push into `ent2tag::Dict` the pair `etype=>etags`;

where:
- `etype::DataType` determines the type of the element (see
    [`_type_tag_to_etype`](@ref));
- `ntags::Matrix{Int}` gives the indices of the nodes defining those
    elements;
- `etags::Vector{Int}` gives the indices of those elements in `elements`.
"""
function _ent_to_mesh!(elements, ent2tag, ω::GmshEntity, shift)
    ω in keys(ent2tag) && (return elements, ent2tag)
    etypes_to_etags = Dict{DataType,Vector{Int}}()
    # Loop on GMSH element types (integer)
    type_tags, _, ntagss = gmsh.model.mesh.getElements(geometric_dimension(ω),gmshtag(ω))
    for (type_tag, ntags) in zip(type_tags, ntagss)
        _, _, _, Np, _ = gmsh.model.mesh.getElementProperties(type_tag)
        ntags = reshape(ntags, Int(Np), :)
        etype = _type_tag_to_etype(type_tag)
        if etype in keys(elements)
            etag = size(elements[etype], 2) .+ collect(1:size(ntags,2))
            ntags = hcat(elements[etype], ntags .+ shift)
        else
            etag = collect(1:size(ntags, 2))
        end
        push!(elements, etype => ntags)
        push!(etypes_to_etags, etype => etag)
    end
    push!(ent2tag, ω => etypes_to_etags)
    return elements, ent2tag
end

function set_meshsize(hmax, hmin=hmax)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", hmin)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", hmax)
end

function set_meshorder(order)
    gmsh.option.setNumber("Mesh.ElementOrder", order)
end

function set_verbosity(i)
    gmsh.option.setNumber("General.Verbosity",i)
end

function summary(model)
    gmsh.model.setCurrent(model)
    @printf("List of entities in model `%s`: \n", model)
    @printf("|%10s|%10s|%10s|\n","name","dimension","tag")
    ents = gmsh.model.getEntities()
    # pgroups = gmsh.model.getPhysicalGroups()
    for ent in ents
        name = gmsh.model.getEntityName(ent...)
        dim, tag = ent
        @printf("|%10s|%10d|%10d|\n", name, dim, tag)
    end
    println()
end

function summary()
    models = gmsh.model.list()
    for model in models
        summary(model)
    end
end

"""
    _type_tag_to_etype(tag)

Mapping of `gmsh` element types, encoded as an integer, to the internal
equivalent of those. This function assumes `gmsh` has been initilized.
"""
function _type_tag_to_etype(tag)
    T = SVector{3,Float64} # point type
    name,dim,order,num_nodes,ref_nodes,num_primary_nodes  = gmsh.model.mesh.getElementProperties(tag)
    num_nodes = Int(num_nodes) #convert to Int64
    if occursin("Point",name)
        etype = LagrangePoint{3,Float64}
    elseif occursin("Line",name)
    etype = LagrangeLine{num_nodes,T}
    elseif occursin("Triangle",name)
        etype = LagrangeTriangle{num_nodes,T}
    elseif occursin("Quadrilateral",name)
        etype = LagrangeSquare{num_nodes,T}
    elseif occursin("Tetrahedron",name)
        etype = LagrangeTetrahedron{num_nodes,T}
    else
        error("unable to parse gmsh element of family $name")
    end
    return etype
end

"""
    _etype_to_type_tag(etype)

Mapping of internal element types, to the integer tag of `gmsh` elements. This
function assumes `gmsh` has been initialized.
"""
function _etype_to_type_tag(el::LagrangeElement)
    etype = typeof(el)
    tag = 1
    while true
        E   = _type_tag_to_etype(tag)
        E === etype && (return tag)
        tag = tag + 1
    end
end
