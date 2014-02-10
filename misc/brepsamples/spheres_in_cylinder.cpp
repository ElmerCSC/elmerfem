// compilation: just include all occ-headers and link against all occ-libs

#include "occheaders.h"

int main()
{
  // make a cylinder:
  //------------------
  gp_Pnt orig(0, 0, 0);
  gp_Dir normal(0, 0, 1);
  Standard_Real radius = 0.25;
  Standard_Real height = 1.0;

  gp_Ax2 axis(orig, normal);
  TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(axis, radius, height);
  
  // make a sphere:
  //---------------
  gp_Pnt sphere_orig(0.15, 0, 0.75);
  Standard_Real sphere_radius = 0.025;

  TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(sphere_orig, sphere_radius);

  // make another sphere:
  //---------------------
  gp_Pnt sphere_orig2(0, 0.15, 0.8);
  Standard_Real sphere_radius2 = 0.025;

  TopoDS_Shape sphere2 = BRepPrimAPI_MakeSphere(sphere_orig2, sphere_radius2);

  // subtract sphere from cylinder:
  //--------------------------------
  BRepAlgoAPI_Cut tmp(cylinder, sphere);
  cylinder = tmp.Shape();

  BRepAlgoAPI_Cut tmp2(cylinder, sphere2);
  cylinder = tmp2.Shape();

  // write result in brep file:
  //---------------------------
  BRepTools::Write(cylinder, "sedimentation.brep");

  return 0;
}
