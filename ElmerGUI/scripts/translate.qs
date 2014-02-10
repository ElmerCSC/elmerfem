egp.ReadPostFile("case.ep");
egp.SetSurfaces(true);
egp.SetIsoSurfaces(true);
egp.SetFeatureEdges(true);
egp.SetText(true);
preferences.SetFeatureAngle(45);
surfaces.SetFieldName("Null");
surfaces.SetOpacity(20);
isoSurfaces.SetFieldName("nodes_z");
isoSurfaces.SetColorName("Temperature");
isoSurfaces.SetContours(true);
isoSurfaces.KeepFieldLimits(true);
text.SetMessage("");
text.SetSize(24);
text.SetBold(true);
egp.SetOrientation(45, 45, 0);
egp.ResetCamera();
egp.Render();
var zmin = egp.GetMinZ();
var zmax = egp.GetMaxZ();
for(var i = 0; i < 100; i++)
{
   var z = zmin + i/100.0 * (zmax-zmin);
   isoSurfaces.SetMinFieldVal(z);
   isoSurfaces.SetMaxFieldVal(z);
   text.SetMessage("z = " + z);
   egp.Redraw();
}
