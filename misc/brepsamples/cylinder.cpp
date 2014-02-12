// compilation: just include all occ-headers and link against all occ-libs

#include "occheaders.h"

int main()
{
  // make a cylinder:
  gp_Pnt orig(0, 0, 0);
  gp_Dir normal(0, 0, 1);
  gp_Ax2 axis(orig, normal);
  Standard_Real radius = 0.2;
  Standard_Real height = 1.0;
  TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(axis, radius, height);
  
  // write in brep file:
  BRepTools::Write(cylinder, "cylinder.brep");

  return 0;
}
