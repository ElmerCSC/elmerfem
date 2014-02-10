egp.ReadPostFile("case.ep")
egp.SetSurfaces(1)
for(var i = 0; i < 180; i++)
{
    egp.RotateY(2.0)
    egp.ResetCamera()
    egp.Render()
    egp.SavePngFile("frame" + i + ".png")
}
