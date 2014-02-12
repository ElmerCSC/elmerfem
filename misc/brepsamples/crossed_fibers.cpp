// compilation: just include all occ-headers and link against all occ-libs

#include "occheaders.h"

int main()
{
  double x = 0.50;
  double y = 1.0 - sqrt(1.0 - x*x);

  // First body:
  //=============

  // make a circle:
  gp_Pnt circleCenterA(0, 0, 0);
  gp_Dir circleNormalA(1, 0, 0);
  gp_Ax1 circleAxisA(circleCenterA, circleNormalA);
  Standard_Real circleRadiusA = 0.05;

  gp_Circ circleA;
  circleA.SetLocation(circleCenterA);
  circleA.SetRadius(circleRadiusA);
  circleA.SetAxis(circleAxisA);
  
  Handle_Geom_Circle geomCircleA = GC_MakeCircle(circleA);
  TopoDS_Edge circleEdgeA = BRepBuilderAPI_MakeEdge(geomCircleA);
  TopoDS_Wire circleEdgeWireA = BRepBuilderAPI_MakeWire(circleEdgeA);
  TopoDS_Face faceA = BRepBuilderAPI_MakeFace(circleEdgeWireA);

  // make a wire path for sweep:
  gp_Pnt p1A(-x, y, 0);
  gp_Pnt p2A( x, y, 0);

  Handle(Geom_TrimmedCurve) tcA = GC_MakeArcOfCircle(p1A, circleCenterA, p2A);
  TopoDS_Edge tcEdgeA = BRepBuilderAPI_MakeEdge(tcA);
  TopoDS_Wire tcEdgeWireA = BRepBuilderAPI_MakeWire(tcEdgeA);

  TopoDS_Shape bodyA = BRepOffsetAPI_MakePipe(tcEdgeWireA, faceA);

  // Second body:
  //==============

  // make a circle:
  gp_Pnt circleCenterB(0, 0.06, 0);
  gp_Dir circleNormalB(0, 0, 1);
  gp_Ax1 circleAxisB(circleCenterB, circleNormalB);
  Standard_Real circleRadiusB = 0.05;

  gp_Circ circleB;
  circleB.SetLocation(circleCenterB);
  circleB.SetRadius(circleRadiusB);
  circleB.SetAxis(circleAxisB);
  
  Handle_Geom_Circle geomCircleB = GC_MakeCircle(circleB);
  TopoDS_Edge circleEdgeB = BRepBuilderAPI_MakeEdge(geomCircleB);
  TopoDS_Wire circleEdgeWireB = BRepBuilderAPI_MakeWire(circleEdgeB);
  TopoDS_Face faceB = BRepBuilderAPI_MakeFace(circleEdgeWireB);

  // make a wire path for sweep:
  gp_Pnt p1B(0, 0.06 - y, -x);
  gp_Pnt p2B(0, 0.06 - y,  x);

  Handle(Geom_TrimmedCurve) tcB = GC_MakeArcOfCircle(p1B, circleCenterB, p2B);
  TopoDS_Edge tcEdgeB = BRepBuilderAPI_MakeEdge(tcB);
  TopoDS_Wire tcEdgeWireB = BRepBuilderAPI_MakeWire(tcEdgeB);

  TopoDS_Shape bodyB = BRepOffsetAPI_MakePipe(tcEdgeWireB, faceB);

  // fuse:
  BRepAlgoAPI_Fuse final(bodyA, bodyB);
  TopoDS_Shape result = final.Shape();
  
  // write body to brep file:
  BRepTools::Write(result, "crossed_fibers.brep");

  return 0;
}
