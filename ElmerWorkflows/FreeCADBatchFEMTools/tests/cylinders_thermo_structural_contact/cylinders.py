# A test case for FreeCAD automatic scripting with thermo-structural (Elmer test is a modified fem/tests/ThermalBiMetal2)
# original date: November 2019
# Author: Eelis Takala
# email: eelis.takala@gmail.com
doc = App.newDocument('Cylinders')
import sys
import os
import time
import math
import FreeCAD
import Part

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import FreeCADBatchFEMTools as FBFT

script_directory = os.path.dirname(__file__)

def cylinder(r1, r2, h, symmetry_planes=None, name='cyl', mesh_size=10.):
    if symmetry_planes == None: symmetry_planes = ['yz','zx', 'xy']
    planes_original = list(symmetry_planes)
    #if symmetry_planes == None: symmetry_planes = ['yz']
    #if symmetry_planes == None: symmetry_planes = ['zx']
    #if symmetry_planes == None: symmetry_planes = ['xy']

    cyl = doc.addObject("PartDesign::Body","Cylinder")
    sketch = cyl.newObject('Sketcher::SketchObject',name+' sketch')
    sketch.MapMode = 'FlatFace'
    sketch.addGeometry(Part.Circle(App.Vector(0.000000,0.000000,0),App.Vector(0,0,1),r2),False)
    sketch.addGeometry(Part.Circle(App.Vector(0.000000,0.000000,0),App.Vector(0,0,1),r1),False)
    doc.recompute()

    pad = cyl.newObject("PartDesign::Pad",name + " pad")
    pad.Profile = sketch
    pad.Length = h
    pad.Length2 = h
    pad.Type = 0
    pad.UpToFace = None
    pad.Reversed = 0
    pad.Midplane = 1
    pad.Offset = 0.000000
    doc.recompute()

    cyl = FBFT.reduce_half_symmetry(cyl, name, App, doc, planes=symmetry_planes)
    doc.recompute()

    # create entities dict
    #face_picks = [('alpha0', 0), ('b2', 1), ('b3', 2), ('b4', 3), ('b5', 4), ('b6', 5)]
    #faces = FBFT.pick_faces_from_geometry(cyl, face_picks)
    #faces = []
    #for plane in planes_original:
        #FBFT.add_symmetry_plane_faces_in_entity_list(faces, cyl, plane)

    #if symmetry_planes == None: symmetry_planes = ['yz']
    #if symmetry_planes == None: symmetry_planes = ['zx']
    #if symmetry_planes == None: symmetry_planes = ['xy']

    face_picks = [('outer', 0), ('yz', 1), ('xy', 2), ('zx', 3), ('end', 4), ('inner', 5)]
    faces = FBFT.pick_faces_from_geometry(cyl, face_picks)

    solids = []
    FBFT.add_entity_in_list(solids, name, cyl, {'mesh size': mesh_size})
    entities_dict = FBFT.create_entities_dict(name, faces, solids, cyl)

    return entities_dict
 
def create_geometry():

    mesh_size_max = 100.
    cyl_entities_list = []
    cyl_entities_list.append(cylinder(20, 10, 40., name='cyl1', mesh_size=3.))
    cyl_entities_list.append(cylinder(30, 20, 20., name='cyl2', mesh_size=4.))
    entities_dict = FBFT.merge_entities_dicts(cyl_entities_list, 'All', default_mesh_size=mesh_size_max,
                                                                  add_prefixes={'solids': False, 'faces': True})
    solid_objects = FBFT.get_solids_from_entities_dict(entities_dict)

    mesh_object, compound_filter = FBFT.create_mesh_object_and_compound_filter(solid_objects, mesh_size_max, doc, separate_boundaries=True)

    FBFT.find_boundaries_with_entities_dict(mesh_object, compound_filter, entities_dict, doc, separate_boundaries=True)
    FBFT.find_bodies_with_entities_dict(mesh_object, compound_filter, entities_dict, doc)
    FBFT.define_mesh_sizes(mesh_object, compound_filter, entities_dict, doc, point_search=False, ignore_list=['air'])

    FBFT.fit_view()
    FBFT.create_mesh(mesh_object)
    FBFT.export_unv(os.path.join(script_directory, 'cylinders.unv'), mesh_object)
     
try:
    create_geometry()

except Exception:
    import traceback
    print(str(traceback.format_exc()))

exit()
