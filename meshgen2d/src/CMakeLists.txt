
SET(meshgen2d_SRCS BGGridMesh.cpp BGTriangleMesh.cpp BGVertex.cpp Body.cpp 
   BoundaryElement.cpp BoundaryLayer.cpp Connect.cpp Element.cpp 
   GeometryEdge.cpp main.cpp Mesh.cpp MeshNode.cpp MeshParser.cpp Node.cpp 
   PQ.cpp QuadElement.cpp QuadLayer.cpp SSSFVoronoiSegment.cpp 
   SSMFVoronoiSegment.cpp TriangleElement.cpp TriangleNELayer.cpp 
   Vertex.cpp VoronoiSegment.cpp VoronoiVertex.cpp VSVertex.cpp)

INCLUDE_DIRECTORIES("include")
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/meshgen2d/src/include")
FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/meshgen2d/src/include")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/meshgen2d/src/include")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/meshgen2d/src")

ADD_EXECUTABLE(Mesh2D ${meshgen2d_SRCS})

INSTALL(TARGETS Mesh2D RUNTIME DESTINATION "bin" RENAME ElmerMesh2D)
#INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/Mesh2D 
              #DESTINATION "bin" 
              #RENAME ElmerMesh2D)
