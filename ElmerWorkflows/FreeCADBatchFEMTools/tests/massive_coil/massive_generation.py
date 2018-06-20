"""
  Massive coil generation and solving test for 
  FreeCADBatchFEMTools - A library for using FreeCAD for FEM preprocessing in batch mode

  Copyright 1st May 2018 - , Trafotek Oy, Finland
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
 
  You should have received a copy of the GNU Lesser General Public
  License along with this library (in file ../LGPL-2.1); if not, write 
  to the Free Software Foundation, Inc., 51 Franklin Street, 
  Fifth Floor, Boston, MA  02110-1301  USA

  Authors: Eelis Takala, Janne Keranen
  Emails:  eelis.takala@gmail.com, janne.sami.keranen@vtt.fi
  Address: Trafotek Oy
           Kaarinantie 700
           20540 Turku
           Finland

  Original Date: May 2018
"""
doc = App.newDocument("Massive winding")
import os
import sys
import Part
import ObjectsFem

PWD = os.path.dirname(os.path.realpath(__file__))
module_path = PWD + "/../.."
sys.path.insert(0, module_path)

from FreeCADBatchFEMTools import fit_view
from FreeCADBatchFEMTools import pick_faces_from_geometry
from FreeCADBatchFEMTools import add_entity_in_list
from FreeCADBatchFEMTools import create_entities_dict
from FreeCADBatchFEMTools import merge_entities_dicts
from FreeCADBatchFEMTools import get_solids_from_entities_dict
from FreeCADBatchFEMTools import create_mesh_object_and_compound_filter
from FreeCADBatchFEMTools import find_boundaries_with_entities_dict
from FreeCADBatchFEMTools import find_bodies_with_entities_dict
from FreeCADBatchFEMTools import define_mesh_sizes
from FreeCADBatchFEMTools import create_mesh
from FreeCADBatchFEMTools import run_elmergrid

def create_core(h, w, a, b, name, mesh_sizes=None):
    
    core = doc.addObject('PartDesign::Body',name+'_obj')
    sketch = core.newObject('Sketcher::SketchObject',name+' sketch')
    #sketch.Support = (doc.XY_Plane, [''])
    sketch.MapMode = 'FlatFace'

    sketch.addGeometry(Part.LineSegment(App.Vector(0,-h/2, 0),App.Vector(0,h/2, 0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(0,h/2, 0),App.Vector(w/2, h/2, 0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(w/2,h/2,0),App.Vector(w/2,a/2,0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(w/2,a/2,0),App.Vector(b/2,a/2,0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(b/2,a/2,0),App.Vector(b/2,-a/2,0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(b/2,-a/2,0),App.Vector(w/2,-a/2,0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(w/2,-a/2,0),App.Vector(w/2,-h/2,0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(w/2,-h/2,0),App.Vector(0,-h/2,0)),False)

    doc.recompute()
    fit_view()

    rev = core.newObject("PartDesign::Revolution",name + " rev")
    rev.Profile = sketch
    rev.ReferenceAxis = (sketch,['V_Axis'])
    rev.Angle=180
    rev.Reversed = 1
    doc.recompute()
#    #fit_view()

    # Here we define the entities dictionary.
    face_picks=[
               ('infinity_xz0',6),
               ('infinity_xz1',0),
               ('cylinder_lateral_surface',1),
               ('cylinder_lateral_surface',5),
               ('xy0',7),
               ('xy0',8)
               ]
    faces = pick_faces_from_geometry(core, face_picks)

    solids = []
    add_entity_in_list(solids, name, core, mesh_sizes=mesh_sizes)

    entities_out = create_entities_dict(name, faces, solids, main_object=core)
    return entities_out
    
def create_coil (r, R, h, z, name, face_entities=True, mesh_sizes=None):
    coil = doc.addObject('PartDesign::Body',name+'_obj')
    sketch = coil.newObject('Sketcher::SketchObject',name+' sketch')
    #sketch.Support = (doc.XY_Plane, [''])
    sketch.MapMode = 'FlatFace'

    sketch.addGeometry(Part.LineSegment(App.Vector(r,-h/2+z, 0),App.Vector(r,h/2+z, 0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(r,h/2+z, 0),App.Vector(R,h/2+z, 0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(R,h/2+z, 0),App.Vector(R,-h/2+z, 0)),False)
    sketch.addGeometry(Part.LineSegment(App.Vector(R,-h/2+z, 0),App.Vector(r,-h/2+z, 0)),False)
    doc.recompute()
    fit_view()

    rev = coil.newObject("PartDesign::Revolution",name + " rev")
    rev.Profile = sketch
    rev.ReferenceAxis = (sketch,['V_Axis'])
    rev.Angle=180
    rev.Reversed = 1
    doc.recompute()

    coil_faces=coil.Shape.Faces

    # Here we define the entities dictionary. Note that each coil contain the same 
    # names for each geometric entity. That is why the 'name' key is required
    # in the dictionary.

    # Naming of the faces can be avoided by giving face_entities=False as argument
    faces=[]
    if(face_entities):
        face_picks=[('gamma0',4),
                    ('gamma1',5)]

        faces = pick_faces_from_geometry(coil, face_picks)

    solids=[]
    add_entity_in_list(solids, name, coil, mesh_sizes=mesh_sizes)
    
    entities_out = create_entities_dict(name, faces, solids, main_object=coil)

    return entities_out

def create_air(coil_inner_radius, coil_outer_radius, coil_height, airgap, z, coil_entities, name, mesh_sizes=None):
    
    airpluscoil_entities=create_coil(coil_inner_radius-airgap, coil_outer_radius+airgap, coil_height+airgap*2, z, name = "air and coil", face_entities=False)
	
    air = doc.addObject("Part::Cut","air"+'_obj')
    air.Base = airpluscoil_entities['main object']
    air.Tool = coil_entities['main object']
    doc.recompute()

    # Here we define the entities dictionary.
    face_picks=[
               ('xy0',3),
               ('xy0',4),
               ('cylinder_lateral_surface',5)
               ]

    faces = pick_faces_from_geometry(air, face_picks)

    solids=[]
    add_entity_in_list(solids, name, air, mesh_sizes=mesh_sizes)
    entities_out = create_entities_dict(name, faces, solids, main_object=air)

    return entities_out


def create_solids(model_parameters):
    coil_inner_radius = model_parameters['coil_inner_radius']
    coil_outer_radius = model_parameters['coil_outer_radius']
    coil_height = model_parameters['coil_height']
    z = model_parameters['z']
    core_height = model_parameters['core_height']
    core_diameter = model_parameters['core_diameter']
    core_center_length = model_parameters['core_center_length']
    core_center_diameter = model_parameters['core_center_diameter']
    airgap = model_parameters['airgap']
    default_mesh_size = model_parameters['default_mesh_size']
    coil_mesh_sizes = model_parameters['coil_mesh_sizes']
    core_mesh_sizes = model_parameters['core_mesh_sizes']
    air_mesh_sizes = model_parameters['air_mesh_sizes']
    
    coil_entities = create_coil(coil_inner_radius, coil_outer_radius, 
                                coil_height, z, name = "L1", mesh_sizes=coil_mesh_sizes)
    core_entities = create_core(core_height,core_diameter,core_center_length,
                                core_center_diameter, name = 'core', mesh_sizes=core_mesh_sizes)
    air_entities = create_air(coil_inner_radius, coil_outer_radius, coil_height, 
                              airgap, z, coil_entities, name = 'air', mesh_sizes=air_mesh_sizes)
    
    doc.recompute()
    name = 'All'
    solids_entities = merge_entities_dicts([core_entities, coil_entities, air_entities], 
            name, default_mesh_size = default_mesh_size, add_prefixes={'solids':False, 'faces':True})
    return solids_entities

##### Tests #####
def test_all():
    model_parameters = {}
    model_parameters['coil_inner_radius'] = 0.11/2
    model_parameters['coil_outer_radius'] = 0.11/2+0.1
    model_parameters['coil_height'] = 0.1
    model_parameters['z'] =.0
    model_parameters['core_height'] =.31 
    model_parameters['core_diameter'] =.32
    model_parameters['core_center_length'] =.11
    model_parameters['core_center_diameter'] =.1
    model_parameters['airgap'] =.005
    model_parameters['default_mesh_size'] = 0.04  #If default mesh size is smaller than given mesh sizes, fallback value is used.
    model_parameters['coil_mesh_sizes'] = {
            'L1':0.02,
            'alpha0':0.01
            }
    model_parameters['core_mesh_sizes'] = {'core':0.04}
    model_parameters['air_mesh_sizes'] = {'air':0.009}
 
    default_mesh_size = model_parameters['default_mesh_size']

    entities_dict = create_solids(model_parameters)
    solid_objects = get_solids_from_entities_dict(entities_dict)

    mesh_object, compound_filter = create_mesh_object_and_compound_filter(solid_objects,default_mesh_size, doc)
    find_boundaries_with_entities_dict(mesh_object, compound_filter, entities_dict, doc)
    find_bodies_with_entities_dict(mesh_object, compound_filter, entities_dict, doc)
    define_mesh_sizes(mesh_object, compound_filter, entities_dict, doc)
    fit_view()
    create_mesh(mesh_object)
    export_path=PWD+"/circuits_harmonic_massive/1962.unv"
    run_elmergrid(export_path, mesh_object)
#################

test_all() #This test is with all the definitions and examples
if not FreeCAD.GuiUp:
	exit()

