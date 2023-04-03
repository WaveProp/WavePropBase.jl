"""
    vtk_mesh_file(mesh::GenericMesh[, Ω::Domain], name::String)

Creates a `VTK` file (.vtu) with name `name` containing the mesh information. It
is possible to export only a `Domain` (i.e. only a part of the mesh).

Note that all the heavy lifting is done by the package
[`WriteVTK`](https://juliavtk.github.io/WriteVTK.jl/stable/); we refer to its
documentation for more information.

!!! note
    The output is not an actual file (on disk). To save it, simply write:
    ```julia
        vtk_mesh_file(mesh, "my_mesh") |> vtk_save
    ```
To add some data (scalar or vector values at the mesh nodes or mesh cells) for
visualization purposes, you can do for instance:

```julia
vtkfile = vtk_mesh_file(mesh, name)
vtkfile["my_point_data", VTKPointData()] = pdata
vtkfile["my_cell_data", VTKCellData()] = cdata
vtk_save(vtkfile)
```

It is possible also to export a partition `Ωs::Vector{Domain}` using
`Multiblock` files (.vtm), for instance like so

    vtmfile = vtk_multiblock(name)
    for (Ω, pdata) in zip(Ωs, pdatas)
        vtkfile = vtk_mesh_file(mesh, Ω, name)
        vtkfile["my_point_data", VTKPointData()] = pdata
    end
    vtk_save(vtmfile)

To save a sequence of solutions (time steps, iterations), simply append the
number of the element to the file name. Paraview will recognize the sequence
automatically.
"""
function vtk_mesh_file(mesh::GenericMesh, name::String)
    points = _vtk_points(mesh)
    cells = _vtk_cells(mesh)
    return vtk_grid(name * ".vtu", points, cells)
end
function vtk_mesh_file(mesh::GenericMesh, Ω::Domain, name::String)
    points = _vtk_points(mesh)
    cells = _vtk_cells(mesh, Ω)
    return vtk_grid(name * ".vtu", points, cells)
end

"""
    _vtk_points(mesh::GenericMesh)

Creates the matrix of nodes in the format required by `WriteVTK`.
"""
function _vtk_points(mesh::GenericMesh)
    vtx = zeros(Float64, 3, length(nodes(mesh)))
    for (i, nd) in enumerate(nodes(mesh))
        vtx[1:ambient_dimension(mesh), i] = nd
    end
    return vtx
end

"""
    _vtk_cells(mesh::GenericMesh, E::DataType)
    _vtk_cells(mesh::GenericMesh, Ω::Domain)

Creates the vector of all elements contained in the mesh in the format
required by `WriteVTK` for a particular element type `E<:AbstractElement`
or a `Domain` instance.
"""
function _vtk_cells end

function _vtk_cells(tags, E::Union{Type{<:LagrangeElement},Type{<:SVector}})
    vtk_cell_type, ind = etype_to_vtk_cell_type(E)
    return [MeshCell(vtk_cell_type, tags[ind, i]) for i in axes(tags, 2)]
end
function _vtk_cells(mesh::GenericMesh, Ω::Domain)
    cells = MeshCell[]
    # Loop on `ElementaryEntity`
    for ω in Ω
        # Loop on `AbstractElement`
        for (E, ind) in ent2tagsdict(mesh)[ω]
            # Subset corresponding to the `ElementaryEntity`
            tags = etype2data(mesh)[E][:, ind]
            append!(cells, _vtk_cells(tags, E))
        end
    end
    return cells
end
function _vtk_cells(mesh::GenericMesh)
    cells = MeshCell[]
    # Loop on `AbstractElement`
    for (E, tags) in etype2data(mesh)
        # Export only the cells of the largest geometrical dimension
        if domain_dimension(E) == ambient_dimension(mesh)
            append!(cells, _vtk_cells(tags, E))
        end
    end
    return cells
end

"""
    etype_to_vtk_cell_type(E::DataType)

Map internal element types to a tuple containing:
- the corresponding `WriteVTK` cell types (following the convention chosen by
  `VTK`, see below);
- the indices in the `elements` column that defines the element. This is because
  we want to support the export of more than just the flat elements available in
  the `VTK` specification, hence which may require a conversion of some sort.

See VTK specification [Fig. 2] on
[http://www.vtk.org/VTK/img/file-formats.pdf](http://www.vtk.org/VTK/img/file-formats.pdf)
"""
function etype_to_vtk_cell_type(E::Type{<:LagrangeElement{T,N}}) where {T,N}
    vtk_type =
    if T === ReferenceLine
        VTKCellTypes.VTK_LAGRANGE_CURVE
    elseif T === ReferenceTriangle
        VTKCellTypes.VTK_LAGRANGE_TRIANGLE
    elseif T === ReferenceSquare
        VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
    elseif T === ReferenceTetrahedron
        VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON
    end
    return vtk_type, collect(1:N)
end
etype_to_vtk_cell_type(E::Type{<:SVector}) = (VTKCellTypes.VTK_VERTEX, collect(1:1))
