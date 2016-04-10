# HexGenerator
-----
## Introduction


- Given a surface sliced into several cylinders, on which there exists a **Strebel Differential**. The **Strebel Differential** induces a quad mesh on the surface. Furthermore, it also induces a **Hex mesh**. This software is used to generate the hexahedron mesh from the strebel differential on the cylinders, with tetrahedron mesh and surface as input.

- However, it still can't merge the hex meshes of the cylinders **automatically**. The *VolumeViewerQt* software can merge the hex cylinders manually, with clicking the corresponding vertices.

## Dependencies

- Comprehensive MeshLib, including **HMeshLib** and **TMeshLib**.
- Eigen

## Authors

- Jerome - Hex mesh generation, HMeshLib, TMeshLib
- David Gu - Triangle Mesh MeshLib