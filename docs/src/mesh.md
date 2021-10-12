# [Mesh module](@id mesh-section)

## Overview 

## Mesh structures

### `GenericMesh`

### `CartesianMesh`

### `SubMesh`

```@example
using GmshSDK
@gmsh begin    
    tag1 = gmsh.model.occ.addSphere(0,0,0,1)
    tag2 = gmsh.model.occ.addSphere(0,0,0,1)
    gmsh.model.occ.synchronize()
    #Ω    = GmshSDK.Domain()
    #msh  = GmshSDK.meshgen(Ω)
end
```

```@index
Modules = [WavePropBase.Mesh]
```