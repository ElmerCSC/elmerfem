// compilation: just include all occ-headers and link against all occ-libs

#include "occheaders.h"

int main()
{
  // set points:
  gp_Pnt p1(-1, -1, -1);
  gp_Pnt p2(+1, -1, -1);
  gp_Pnt p3(+1, +1, -1);
  gp_Pnt p4(-1, +1, -1);
  
  // make edges from points:
  TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(p1, p2);
  TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(p2, p3);
  TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(p3, p4);
  TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(p4, p1);
  
  // make wire frame from edges:
  TopoDS_Wire wire = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);
  
  // make face from wire frame:
  TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);
  
  // make cube from face (extrude):
  gp_Vec vec(0 , 0 , +2);
  TopoDS_Shape cube = BRepPrimAPI_MakePrism(face , vec);
  
  // fillet edges of cube:
  BRepFilletAPI_MakeFillet fillet(cube);

  TopExp_Explorer explorer(cube, TopAbs_EDGE);  
  while(explorer.More())
    {
      TopoDS_Edge edge = TopoDS::Edge(explorer.Current());
      fillet.Add(0.1, edge);
      explorer.Next();
    }
  
  cube = fillet.Shape();
  
  // write filleted cube to brep file:
  BRepTools::Write(cube, "filleted_cube.brep");

  return 0;
}
