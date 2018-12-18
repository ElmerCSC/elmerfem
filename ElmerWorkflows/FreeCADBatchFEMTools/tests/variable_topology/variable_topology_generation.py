"""
  Geometry generation for with variable topology and solving test for 
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
doc = App.newDocument("Variable_geometry")
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
from FreeCADBatchFEMTools import remove_compare_faces_from_list
from FreeCADBatchFEMTools import faces_with_vertices_in_symmetry_plane
from FreeCADBatchFEMTools import create_xor_object

def reduce_half_symmetry(solid, name, planes=None, reversed_direction = False):
    doc.recompute()
    if planes==None: return solid
    plane = planes.pop()
    doc.recompute()
    reduced_name = name + '_' + plane
    tool_box = doc.addObject("Part::Box","CutBox"+reduced_name)
    x = 10. * solid.Shape.BoundBox.XLength
    y = 10. * solid.Shape.BoundBox.YLength
    z = 10. * solid.Shape.BoundBox.ZLength
    center=solid.Shape.CenterOfMass
    tool_box.Length = x
    tool_box.Width = y 
    tool_box.Height = z
    if plane == 'zx':
        tool_box.Placement = App.Placement(App.Vector(center.x-x/2.,0,center.z-z/2.),App.Rotation(App.Vector(0,0,1),0))
    elif plane == 'xy':
        tool_box.Placement = App.Placement(App.Vector(center.x-x/2.,center.y-y/2.,0),App.Rotation(App.Vector(0,0,1),0))
    elif plane == 'yz':
        tool_box.Placement = App.Placement(App.Vector(0,center.y-y/2.,center.z-z/2.),App.Rotation(App.Vector(0,0,1),0))
    else:
        raise ValueError("Wrong keyword for plane variable, should be: zx, xy or yz!")
    
    if reversed_direction:
        half_symmetry = doc.addObject("Part::MultiCommon",reduced_name)
        half_symmetry.Shapes = [solid, tool_box]
    else:
        half_symmetry = doc.addObject("Part::Cut", reduced_name)
        half_symmetry.Base = solid
        half_symmetry.Tool = tool_box

    if len(planes) > 0:
        return reduce_half_symmetry(half_symmetry, reduced_name, planes, reversed_direction)

    return half_symmetry

def create_bar(box_size, bar_lenght, center, name, symmetry_planes=None, mesh_size=None):
    bar_obj = doc.addObject('Part::Box', name+'_obj')
    bar_obj.Height = box_size
    bar_obj.Width = bar_lenght
    bar_obj.Length = box_size
    bar_obj.Placement.Base = FreeCAD.Vector(center)- FreeCAD.Vector(box_size/2,0,box_size/2)
    doc.recompute()
     
    bar_obj = reduce_half_symmetry(bar_obj, name, planes=symmetry_planes)

    doc.recompute()
    face_picks=[]
    face_picks=[('alpha0',0),
                ('alpha1',5),
                ('beta0',3),
                ('beta1',1),
                ('gamma0',4),
                ('gamma1',2)]
 
    faces = pick_faces_from_geometry(bar_obj, face_picks)

    solids = []
    add_entity_in_list(solids, name, bar_obj, {'mesh size':mesh_size})
    entities_dict = create_entities_dict(name, faces, solids, main_object=bar_obj)

    return entities_dict

def create_air(bars, x, y, z, name, symmetry_planes=None, mesh_size=None):
    
    air_with_bars = doc.addObject('Part::Box', name+ '_whole')
    air_with_bars.Height = z
    air_with_bars.Width = y
    air_with_bars.Length = x
    air_with_bars.Placement.Base = FreeCAD.Vector(-x/2,-y/2,-z/2)
    air_with_bars = reduce_half_symmetry(air_with_bars, name, planes=symmetry_planes)

    bar_objs = [bar['main object'] for bar in bars]
    solid_objects = [air_with_bars]+bar_objs

    air = create_xor_object(solid_objects, doc)

    doc.recompute()

    faces_in_symmetry_plane = faces_with_vertices_in_symmetry_plane(air.Shape.Faces, plane='zx')

    # remove faces that belong to bar
    bar_face_objects = []
    for bar in bars:
        bar_face_objects += [face_entity['geometric object'] for face_entity in bar['faces']]

    remove_compare_faces_from_list(bar_face_objects, faces_in_symmetry_plane)

    faces = []
    add_entity_in_list(faces, 'gamma1', faces_in_symmetry_plane[0])

    solids = []
    add_entity_in_list(solids, name, air, {'mesh size':mesh_size})
    entities_dict = create_entities_dict(name, faces, solids, main_object=air)
    doc.recompute()
    return entities_dict


def create_solids(model_parameters):
    x = model_parameters['x']
    y = model_parameters['y']
    z = model_parameters['z']
    n1 = model_parameters['n1'] 
    n2 = model_parameters['n2']
    box_size = model_parameters['box_size']
    airgap = model_parameters['airgap']
    default_mesh_size = model_parameters['default_mesh_size']
    bar_mesh_sizes = model_parameters['bar_mesh_sizes']
    air_mesh_sizes = model_parameters['air_mesh_sizes']
    bar_centers = []
    for j in range (0, n2):
        for i in range (0, n1):
            bar_centers.append(((airgap+box_size/2)+i*(airgap+box_size)-x/2,-y/2,(airgap+box_size/2)+j*(airgap+box_size)-z/2))
    bars = []
    for i,center in enumerate(bar_centers):
        name = "bar_" + str(i)
        bar = create_bar(box_size, y, center, name = name, symmetry_planes=['zx'], mesh_size=bar_mesh_sizes)
        bars.append(bar)

    air_entities = create_air(bars, x, y, z, name = 'air', symmetry_planes=['zx'], mesh_size=air_mesh_sizes)

    entities_dicts = [air_entities] + bars

    entities_dict = merge_entities_dicts(entities_dicts, 
            'variable_topology', default_mesh_size=default_mesh_size)

    doc.recompute()

    return entities_dict

def create_sif(model_parameters):

    T_BC=model_parameters['T_BC']
    n1=model_parameters['n1']
    n2=model_parameters['n2']

    filepart1='''
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "variable_topology"
  Include variable_topology/mesh.names
  Results Directory ""
End

Simulation
  Max Output Level = 5
  Coordinate System = Cartesian !3D
  Coordinate Mapping(3) = 1 2 3
  
  Steady State Max Iterations = 1
  Output Intervals = 1
  Simulation Type = Steady state !Transient
  Timestepping Method = BDF
  BDF Order = 1
  Post File = variable_topology.vtu

  Timestep Sizes(1) = 0.1
  Timestep Intervals(1) = 20
  Output Intervals(1) = 1

End 


Constants
  Gravity(4) = 0 -1 0 9.82
  Stefan Boltzmann = 5.67e-08
End

Body 1
  Target Bodies(1) = $ air
  Name = "Air"
  Equation = 1
  Material = 1
  Initial Condition = 1
End

Body 2
  Target Bodies(''' 

    filepart1+=str(n1*n2)+ ') = $ '

    for i in range (0, n1*n2-1):
        filepart1=filepart1+'bar_'+str(i)+' \\ \n    '
    filepart1=filepart1+'bar_'+str(i+1)

    filepart2='''
  Name = "Bars"
  Equation = 1
  Material = 2

  Initial Condition = 1
End

Solver 1
  Equation = Heat Equation
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Exec Solver = Always
  Stabilize = True
  Bubbles = False
  Lumped Mass Matrix = False
  Optimize Bandwidth = True
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-7
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Relaxation Factor = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-10
  BiCGstabl polynomial degree = 2
  Linear System Preconditioning = ILU1
  Linear System ILUT Tolerance = 1.0e-3
  Linear System Abort Not Converged = False
  Linear System Residual Output = 20
  Linear System Precondition Recompute = 1
End

Solver 2 
  Exec Solver = Always
  Equation = "save scalars"
  Procedure = "SaveData" "SaveScalars"
!  Variable 1 = Temperature
!  Operator 1 = int mean
  Variable 1 = Temperature
  Operator 1 = max abs
  Variable 2 = Temperature
  Operator 2 = min abs
  Variable 3 = Time
  Filename = scalars.dat
End 

Equation 1
  Name = "HeatEquation"
  Active Solvers(1) = 1
  Convection = "constant"
End

Initial Condition 1
  Name = "InitialState"
  Temperature = 0.0
End


Material 1
  Name = "Air"
  Heat Conductivity = Real 1.0 !0.026
  Density = Real 1.0 !1.1644

  Heat Capacity = Real 1.0 !0.00121
End

Material 2
  Name = "Copper"
  Heat Conductivity = Real 1.0 !401
  Density = Real 1.0 !8940

  Heat Capacity = Real 1.0 !3.45
End


Boundary Condition 1
  Target Boundaries(1) = $ air_gamma1
  Name = "Symmetry_air"

  Heat Flux = Real 0.0
End

    '''

    BCstring1='''
Boundary Condition '''

    BCstring2='''
  Target Boundaries(1) = $ bar_'''

    BCstring3='''_gamma1
  Name = "Bar '''
    BCstring4=''' BC"
  Temperature = Real '''

    BCstring5='''
End

'''

    filepart3=""
    num=0
    for j in range (0, n2):
        for i in range (0, n1):    
            filepart3=filepart3+BCstring1+str(num+2)+BCstring2+str(num)+BCstring3+str(num)+BCstring4+str(T_BC[j][i])+BCstring5
            num=num+1

    export_path=PWD+u"/variable_topology.sif"
    file = open(export_path, "w+") 
    file.write(filepart1+filepart2+filepart3)
    file.close()

def get_test_params(test_number=1):
    if test_number == 1:
        n1=5
        n2=2
        T_BC=[[1.,2.,3.,4., 5.],[2.,3.,4.,5.,6.]]
    elif test_number == 2:
        n1=5
        n2=3
        T_BC=[[1.,2.,3.,4., 5.],[2.,3.,4.,5.,6.], [3.,4.,5.,6.,7.]]
    elif test_number == 3:
        n1=4
        n2=2
        T_BC=[[1.,2.,3.,4.],[2.,3.,4.,5.]]
    elif test_number == 4:
        n1=3
        n2=3
        T_BC=[[1.,2.,3.],[2.,3.,4.], [3.,4.,5]]
    elif test_number == 5:
        n1=10
        n2=1
        T_BC=[[1.,2.,3.,4., 5.,2.,3.,4.,5.,6.]]
    elif test_number == 6:
        n1=2
        n2=1
        T_BC=[[1.,6.]]
    return n1, n2, T_BC

##### Tests #####
def test_all():
    n1, n2, T_BC = get_test_params(6)
    box_size=10.
    airgap=5.
    depth=25.
    model_parameters = {}
    model_parameters['n1']=n1
    model_parameters['n2']=n2
    model_parameters['box_size']=box_size
    model_parameters['airgap']=airgap
    model_parameters['T_BC']=T_BC
    model_parameters['x'] =n1*box_size+(n1+1)*airgap
    model_parameters['y'] =depth
    model_parameters['z'] =n2*box_size+(n2+1)*airgap
    model_parameters['default_mesh_size'] = 10  #If default mesh size is smaller than given mesh sizes, fallback value is used.
    model_parameters['bar_mesh_sizes'] = 2
    model_parameters['air_mesh_sizes'] = 2
 
    default_mesh_size = model_parameters['default_mesh_size']

    entities_dict = create_solids(model_parameters)
    solid_objects = get_solids_from_entities_dict(entities_dict)

    fit_view()

    mesh_object, compound_filter = create_mesh_object_and_compound_filter(solid_objects, default_mesh_size, doc)

    find_boundaries_with_entities_dict(mesh_object, compound_filter, entities_dict, doc)
    find_bodies_with_entities_dict(mesh_object, compound_filter, entities_dict, doc)

    define_mesh_sizes(mesh_object, compound_filter, entities_dict, doc)
    create_mesh(mesh_object)
    export_path=PWD+u"/variable_topology.unv"
    run_elmergrid(export_path, mesh_object)
    create_sif(model_parameters)
#################

test_all() #This test is with all the definitions and examples
if not FreeCAD.GuiUp:
	exit()

