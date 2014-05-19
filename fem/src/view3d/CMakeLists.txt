SET(view3d_SOURCES ViewFactors.c TestModel.c 
    LineUtil.c TriangleUtil.c BiLinearUtil.c 
    BiQuadraticUtil.c BiCubicUtil.c RayTrace.c 
    VectorUtil.c LUDecomp.c 2d_4node.c 16node_quad.c second.c)

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/fem/src/view3d")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}")
ADD_LIBRARY(view3d STATIC ${view3d_SOURCES})
