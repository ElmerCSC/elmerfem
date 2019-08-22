"""
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
from __future__ import print_function
import Fem
import FreeCAD
import BOPTools.SplitFeatures
import ObjectsFem
import femmesh.gmshtools
import math
import itertools


def fit_view():
    """
    If GUI is available, fit the view so that the geometry can be seen
    """
    if FreeCAD.GuiUp:
        import FreeCADGui
        FreeCADGui.ActiveDocument.activeView().viewAxonometric()
        FreeCADGui.SendMsgToActiveView("ViewFit")

def isclose(a, b, rel_tol=1e-4, abs_tol=1e-4):
    """
    Returns True if a and b are close to each other (within absolute or relative tolerance).

    :param a: float
    :param b: float
    :param rel_tol: float
    :param abs_tol: float
    :return: bool
    """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def vectors_are_same(vec1,vec2,tol=1e-4):
    """
    Compares vectors vec1 and vec2. If they are same within a tolerance returns
    True if not returns false.

    :param vec1: Vector 1 
    :param vec2: Vector 2 to be compared with Vector 1
    :return: Boolean
    """
    vec3 = vec1.sub(vec2)
    return isclose(vec3.Length, 0., abs_tol=tol)

def faces_with_vertices_in_symmetry_plane(face_object_list, plane=None, abs_tol=1e-4):
    """
    Returns faces from a list of FreeCAD face objects. The returned faces have to 
    be in a defined symmetry plane. The face is in symmetry plane if all of its points
    and the center of mass are in the plane.

    :param face_object_list: list of FreeCAD face objects
    :param plane: symmetry plane ('zx', 'xy' or 'yz')

    :return: list of FreeCAD face objects that are in the given symmetry plane
    """
    if plane is None: return None
    face_object_list_out = []
    for face_object in face_object_list:
        vertices = face_object.Vertexes
        center_of_mass = face_object.CenterOfMass
        if plane=='zx': center_compare_value = center_of_mass.y
        elif plane=='xy': center_compare_value = center_of_mass.z
        elif plane=='yz': center_compare_value = center_of_mass.x
        else: raise ValueError("Wrong keyword for plane variable, should be: zx, xy or yz!")
        for i, vertex in enumerate(vertices):
            if plane=='zx': compare_value = vertex.Y
            elif plane=='xy': compare_value = vertex.Z
            elif plane=='yz': compare_value = vertex.X
            else: raise ValueError("Wrong keyword for plane variable, should be: zx, xy or yz!")
            if not isclose(compare_value, 0., abs_tol=abs_tol): break
        if i==len(vertices)-1 and isclose(center_compare_value, 0., abs_tol=abs_tol): face_object_list_out.append(face_object)
    return face_object_list_out

def faces_same_center_of_masses(face1, face2, tolerance=0.0001):    
    """
    Compare two faces by comparing if they have same centers of mass with the tolerance.

    :param face1: FreeCAD face object
    :param face2: FreeCAD face object
    :param tolerance: float

    """
    return vectors_are_same(face1.CenterOfMass, face2.CenterOfMass, tolerance)

def faces_are_same(face1, face2, tolerance=1e-4):
    """
    Return true if face1 is same as face2. The faces are same if they have
    the same center of masses and same vertices.

    :param face1: FreeCAD face object
    :param face2: FreeCAD face object

    :return: bool
    """
    return faces_same_center_of_masses(face1, face2, tolerance) and faces_have_same_vertices(face1, face2, tolerance)

def is_face_in_list(search_face, face_object_list, tolerance=1e-4):
    """
    Returns true if search_face is in the face_object_list. Compares faces with 
    face_compare method.

    :param search_face: FreeCAD face object
    :param face_object_list: list of FreeCAD face objects
    """
    for face_object in face_object_list:
        if faces_are_same(search_face, face_object): return True
    return False

def remove_compare_face_from_list(cface, face_object_list, tolerance=1e-4):
    """
    Removes the first FreeCAD face object that matches in the FreeCAD face object list. 
    Uses face_compare to determine if the face is to be removed. 

    :param cface: a FreeCAD face object to be compared
    :param face_object_list: the list of FreeCAD face objects where the face 
                             is to be removed in case of a match

    :return: list of FreeCAD face objects that are removed from the original list of 
             FreeCAD face objects.
    """
 
    for i, face_object in enumerate(face_object_list):
        if faces_are_same(cface, face_object): 
            return face_object_list.pop(i)
    return None

def remove_compare_faces_from_list(compare_face_object_list, face_object_list):
    """
    Removes all the face objects in compare_face_object_list that match to the face objects in 
    the face_object_list. Uses face_compare to determine if the face is to be removed. 

    :param compare_face_object_list: list of FreeCAD face objects to be compared
    :param face_object_list: original list of FreeCAD face objects

    :return: list of FreeCAD face objects that are removed from the original list of 
             FreeCAD face objects.
    """
    removed = []
    for face_object in compare_face_object_list:
        removed.append(remove_compare_face_from_list(face_object, face_object_list))
    return removed 

def faces_have_same_vertices(face1, face2, tolerance=0.0001):
    """
    Compare two faces by comparing that they have same number of vertices and 
    the vertices are in identical coordinates with the tolerance. Return 
    truth value to the faces are the same in this regard.

    :param face1: FreeCAD face object
    :param face2: FreeCAD face object
    :return: bool
    """
    face_vertices_found = []
    for vertex in face2.Vertexes:
        for cvertex in face1.Vertexes:
            if vectors_are_same(vertex.Point,cvertex.Point, tolerance):                   
                face_vertices_found.append(1)
    return len(face_vertices_found) == len(face2.Vertexes) and len(face_vertices_found) == len(face1.Vertexes)

def is_point_inside_face(face, vector, tolerance=0.0001):
    """
    Returns True if point is inside face.

    WARNING: This function calls function face.isInside which does NOT respect tolerance
             https://forum.freecadweb.org/viewtopic.php?t=31524

    :param face: FreeCAD face object
    :param vector: Vector
    :param tolerance: float

    :return: bool
    """
    return face.isInside(vector, tolerance, True)

def is_point_inside_solid(solid, vector, tolerance=0.0001, include_faces=True):
    """
    Returns True if point is inside solid.

    :param solid: FreeCAD solid object
    :param vector: Vector
    :param tolerance: float
    :param include_faces: bool

    :return: bool
    """
    return solid.isInside(vector, tolerance, include_faces)

def is_point_inside_solid_with_round(solid, vector, tolerance=0.0001, round_digits=6):
    """
    Returns True if point is inside solid (faces included) with
    additional tolerance (8 points checked).
    Tries upper and lower rounding of coordinates with precision round_digits.

    :param solid: FreeCAD solid object
    :param vector: Vector
    :param tolerance: float
    :param round_digits: integer

    :return: bool
    """
    rounding = 10**round_digits
    x_floor, x_ceil = math.floor(rounding*vector.x)/rounding, math.ceil(rounding*vector.x)/rounding
    y_floor, y_ceil = math.floor(rounding*vector.y)/rounding, math.ceil(rounding*vector.y)/rounding
    z_floor, z_ceil = math.floor(rounding*vector.z)/rounding, math.ceil(rounding*vector.z)/rounding
    for coordinates in itertools.product([x_floor, x_ceil], [y_floor, y_ceil], [z_floor, z_ceil]):
        vector.x = coordinates[0]
        vector.y = coordinates[1]
        vector.z = coordinates[2]
        if is_point_inside_solid(solid, vector, tolerance):
            return True
    return False

def is_edge_in_solid(solid, edge, tolerance=0.0001):
    """
    Returns True if edge inside solid by comparing is edge vertices inside solid.

    :param solid: FreeCAD solid object
    :param edge: FreeCAD edge object
    :param tolerance: float

    :return: bool
    """
    for vertex in edge.Vertexes:
        if not is_point_inside_solid(solid, vertex.Point, tolerance):
            return False
    return True

def get_point_from_solid(solid, tolerance=0.0001):
    """
    Returns point from given solid.

    :param solid: FreeCAD solid object
    :param tolerance: float

    :return: None or FreeCAD vector object
    """
    x_min, y_min, z_min = solid.BoundBox.XMin, solid.BoundBox.YMin, solid.BoundBox.ZMin
    x_len, y_len, z_len = solid.BoundBox.XLength, solid.BoundBox.YLength, solid.BoundBox.ZLength
    for split_count in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
                        197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
                        311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
                        431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541]:
        x_split_len, y_split_len, z_split_len = x_len/split_count, y_len/split_count, z_len/split_count
        for i in range(1, split_count):
            x_test = x_min + i*x_split_len
            for j in range(1, split_count):
                y_test = y_min + j*y_split_len
                for k in range(1, split_count):
                    test_point = FreeCAD.Vector(x_test, y_test, z_min + k*z_split_len)
                    if is_point_inside_solid(solid, test_point, tolerance, include_faces=False):
                        return test_point
    return None

def get_point_from_face(face):
    """
    Returns point from given face.

    :param face: FreeCAD face object.

    :return: None or FreeCAD vector object
    """
    u_min, u_max, v_min, v_max = face.ParameterRange
    u_len, v_len = u_max-u_min, v_max-v_min
    # use primes so same points are not checked multiple times
    for split_count in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
                        197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
                        311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
                        431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541]:
        u_split_len, v_split_len = u_len/float(split_count), v_len/float(split_count)
        for i in range(1, split_count):
            u_test = u_min + i*u_split_len
            for j in range(1, split_count):
                v_test = v_min + j*v_split_len
                if face.isPartOfDomain(u_test, v_test):
                    return face.valueAt(u_test, v_test)
    return None

def is_face_in_face(face1, face2, tolerance=0.0001):
    """
    Returns True if all vertices and point in face1 also belongs to face2

    :param face1: FreeCAD face object
    :param face2: FreeCAD face object
    :param tolerance: float

    :return: bool
    """
    for vertex in face1.Vertexes:
        if not is_point_inside_face(face2, vertex.Point, tolerance):
            return False
    point_in_face1 = get_point_from_face(face1)
    if point_in_face1 is not None:
        if not is_point_inside_face(face2, point_in_face1, tolerance):
            return False
    else:
        raise ValueError('Face point not found')
    return True

def is_face_in_solid(solid, face, tolerance=0.0001, use_round=True):
    """
    Checks if all face vertices and one additional point from face are inside solid.
    If use_round is True calls function meth:`is_point_inside_solid_with_round` to check
    if point inside solid.

    :param solid: FreeCAD solid object
    :param face: FreeCAD face object
    :param tolerance: float
    :param use_round: bool

    :return: bool
    """
    if use_round:
        is_point_in_solid_func = is_point_inside_solid_with_round
    else:
        is_point_in_solid_func = is_point_inside_solid
    for vertex in face.Vertexes:
        if not is_point_in_solid_func(solid, vertex.Point, tolerance):
            return False
    point_in_face = get_point_from_face(face)
    if point_in_face is not None:
        if not is_point_in_solid_func(solid, point_in_face, tolerance):
            return False
    else:
        raise ValueError('Face point not found')
    return True

def is_compound_filter_solid_in_solid(compound_filter_solid, solid, tolerance=0.0001, point_search=True):
    """
    If point_search is True:
        returns True if all faces of compound_filter_solid are inside solid
    else:
        returns compound_filter_solid.common(solid).Volume > 0

    :param compound_filter_solid: FreeCAD solid object
    :param solid: FreeCAD solid object
    :param tolerance: float (only used with point_search)
    :param point_search: bool

    :return: bool
    """
    if point_search:
        point_in_solid = get_point_from_solid(compound_filter_solid, tolerance)
        if point_in_solid is None:
            raise ValueError('Solid point not found')
        return is_point_inside_solid(solid, point_in_solid, tolerance)
    return compound_filter_solid.common(solid).Volume > 0

def solids_are_the_same(solid1, solid2):
    """
    Compare two solids by comparing have they same number of faces and the faces are identical. 
    Return truth value to the solids are the same in this regard.
    
    :param solid1: FreeCAD solid object
    :param solid2: FreeCAD solid object
    :return: bool
    """
    solid_faces_found = []
    for face in solid2.Faces:
        for cface in solid1.Faces:
            if faces_same_center_of_masses(cface, face) and faces_have_same_vertices(cface, face):                   
                solid_faces_found.append(1)
    return len(solid_faces_found) == len(solid2.Faces) and len(solid_faces_found) == len(solid1.Faces)

def create_boolean_compound(solid_objects, doc):
    """
    Creates a FreeCAD boolean compound for the list of FreeCAD solid objects.
    This is needed when mesh is computed for the whole geometry. Note that
    there is also a create_mesh_object_and_compound_filter for meshing purpose.

    :param solid_objects: list of FreeCAD solid geometry objects.
    :param doc: FreeCAD document object.
    :return: FreeCAD compound object.
    """
    doc.recompute()
    comp_obj = BOPTools.SplitFeatures.makeBooleanFragments(name='Compsolid')
    comp_obj.Objects = solid_objects 
    comp_obj.Mode = "CompSolid"
    comp_obj.Proxy.execute(comp_obj)
    comp_obj.purgeTouched()
    return comp_obj

def create_xor_object(solid_objects, doc):
    """
    Creates a FreeCAD xor object for the list of FreeCAD solid objects.

    :param solid_objects: list of FreeCAD solid geometry objects.
    :param doc: FreeCAD document object.
    :return: FreeCAD xor object
    """
    doc.recompute()
    xor_object = BOPTools.SplitFeatures.makeXOR(name='XOR')
    xor_object.Objects = solid_objects
    xor_object.Proxy.execute(xor_object)
    xor_object.purgeTouched()
    return xor_object

def create_compound_filter(compsolid):
    """
    Create a compound filter. This is needed for meshing. Note that
    there is also a create_mesh_object_and_compound_filter for meshing purpose.

    :param compsolid: FreeCAD compound object for example from create_boolean_compound
    :return: FreeCAD compound filter object
    """
    import CompoundTools.CompoundFilter
    compound_filter = CompoundTools.CompoundFilter.makeCompoundFilter(name='CompoundFilter')
    compound_filter.Base = compsolid
    compound_filter.FilterType = 'window-volume' #???
    compound_filter.Proxy.execute(compound_filter) #???
    return compound_filter

def create_mesh_object(compound_filter, CharacteristicLength, doc):
    """
    Creates a mesh object that controls the mesh definitions.

    :param compound_filter: FreeCAD compound filter object
    :param CharacteristicLength: Default mesh size characteristic length
    :param doc: FreeCAD document object.
    :return: FreeCAD mesh object
    """
    # Make a FEM mesh and mesh groups for material bodies and boundary conditions
    mesh_object = ObjectsFem.makeMeshGmsh(doc, 'GMSHMeshObject')
    mesh_object.Part = compound_filter
    mesh_object.CharacteristicLengthMax = CharacteristicLength
    mesh_object.Algorithm3D = 'New Delaunay'
    mesh_object.ElementOrder = u"1st"
    return mesh_object

def create_mesh(mesh_object):
    """
    Create mesh mesh with Gmsh.

    :param mesh_object: FreeCAD mesh object
    """
    gmsh_mesh = femmesh.gmshtools.GmshTools(mesh_object)
    error = gmsh_mesh.create_mesh()
    print(error)

def create_mesh_object_and_compound_filter(solid_objects, CharacteristicLength, doc):
    """
    Creates FreeCAD mesh and compound filter objects. Uses create_boolean_compound and 
    create_compound_filter, create_mesh_object methods.

    :param solid_objects: list of FreeCAD solid geometry objects
    :param CharacteristicLength: Default mesh size characteristic length
    :param doc: FreeCAD document object.
    """
    boolean_compound = create_boolean_compound(solid_objects, doc)
    compound_filter = create_compound_filter(boolean_compound)
    mesh_object = create_mesh_object(compound_filter, CharacteristicLength, doc)
    return mesh_object, compound_filter

def run_elmergrid(export_path, mesh_object, out_dir=None):
    """
    Run ElmerGrid as an external process if it found in the operating system.

    :param export_path: path where the result is written
    :param mesh_object: FreeCAD mesh object that is to be exported
    :param out_dir: directory where to write mesh files (if not given unv file name is used)
    """
    # Export to UNV file for Elmer
    export_objects = [mesh_object]
    Fem.export(export_objects, export_path)

    optional_path_to_elmer = ''  # optional path to ElmerGrid
    elmerGrid_command = optional_path_to_elmer + 'ElmerGrid 8 2 ' + export_path + ' -autoclean -names'
    if out_dir is not None:
        elmerGrid_command += ' -out ' + out_dir

    from PySide import QtCore, QtGui
    try:
            process = QtCore.QProcess()
            process.startDetached(elmerGrid_command)
            FreeCAD.Console.PrintMessage('Running ' + elmerGrid_command + '\n')
            FreeCAD.Console.PrintMessage('Finished ' + elmerGrid_command + '\n')
    except:
            FreeCAD.Console.PrintError('Error')
            QtGui.QMessageBox.critical(None, 'Error', 'Error!!', QtGui.QMessageBox.Abort)

def export_unv(export_path, mesh_object):
    """
    Exports UNV file for Elmer.

    :param export_path: string
    :param mesh_object: Mesh object
    """
    Fem.export([mesh_object], export_path)

def find_compound_filter_boundary(compound_filter, face):
    """
    Find which face in the compound filter object is the face as the one given in second argument
    Returns the name of the face in compound filter.

    :param compound_filter: FreeCAD compound filter
    :param face: FreeCAD face object
    :return: string
    """
    faces = compound_filter.Shape.Faces
    face_found = None
    for num, cface in enumerate(faces):
        if faces_have_same_vertices(cface, face):
            face_found = num
    if face_found is None: return None
    string = "Face" + str(face_found+1)
    return string

def find_compound_filter_solid(compound_filter, solid):
    """
    Find which solid in the compound filter object is the solid as the one given in second argument
    Returns the name of the solid in compound filter.

    :param compound_filter: FreeCAD compound filter
    :param solid: FreeCAD solid object
    :return: string
    """
    solids = compound_filter.Shape.Solids
    solid_found = None
    for num, csolid in enumerate(solids):
        if solids_are_the_same(csolid, solid):
            solid_found = num
    if solid_found is None: return None
    string = "Solid" + str(solid_found+1)
    return string

def find_compound_filter_boundaries(compound_filter, face):
    """
    Finds all faces in the compound filter object which are inside given face.
    Returns a tuple containing all names of the faces in compound filter.

    :param compound_filter: FreeCAD compound filter
    :param face: FreeCAD face object

    :return: tuple
    """
    face_name_list = []
    for num, cface in enumerate(compound_filter.Shape.Faces):
        if is_face_in_face(cface, face):
            face_name_list.append("Face" + str(num+1))
    if len(face_name_list) == 0:
        raise ValueError("Faces not found")
    return tuple(face_name_list)

def find_compound_filter_solids(compound_filter, solid, point_search=True):
    """
    Finds all solids in the compound filter object which are inside given solid.
    Returns a tuple containing all names of the solids in compound filter.

    :param compound_filter: FreeCAD compound filter
    :param solid: FreeCAD solid object
    :param point_search: bool

    :return: tuple
    """
    solid_name_list = []
    for num, csolid in enumerate(compound_filter.Shape.Solids):
        if is_compound_filter_solid_in_solid(csolid, solid, point_search=point_search):
            solid_name_list.append("Solid" + str(num+1))
    if len(solid_name_list) == 0:
        raise ValueError("Solids not found")
    return tuple(solid_name_list)

""" 
There is no topological naming in FreeCAD. Further, the face numbers are changed in 
making boolean compound. Hence, then making the original solids, the important solids 
and faces are named with generating the following entities dictionary: 

Entities dict definition 
    example_entities = {
      'name' : 'This is the name of the example',
      'faces' : [
                 {'name':'face1',
                  'geometric object':face1_geom_object,
                  'mesh size':mesh_size},
                  ..., 
                 {'name':'facen',
                  'geometric object':facen_geom_object,
                  'mesh size':mesh_size}
                ]
      'solids' : [
                  {'name':'solid_1',
                   'geometric object':solid1_geom_object,
                   'mesh size':mesh_size_1}, 
                   ...,
                  {'name':'solid_n',
                   'geometric object':solidn_geom_object,
                   'mesh size':mesh_size_n}
                 ]
      'main object': main_geom_object
     }
In 'faces' and 'solids' have lists that are so called entity lists.
Entity is defined as a dict containing the name, geometric object and 
the mesh size. In principle one could dynamically add more keys in this
entity property dict:

                  {'name':'solid_n',
                   'geometric object':solidn_geom_object,
                   'mesh size':mesh_size_n}

main_geom_object is for storing a main geometry. For example usually 
when a single solid is created, it contains many face and one solid. 
it is handy to store the solid under the 'main object' key.
"""
def add_entity_in_list(entity_list, name, geom_object, mesh_sizes=None):
    """
    Add entity in list of entities. The mesh sizes can be defined by 
    providing the following dictionary:

    mesh_sizes={
                entity_name_1:mesh_size_for_entity_name_1,
                ...
                entity_name_n:mesh_size_for_entity_name_n,
                mesh size:mesh size if name is not in the dict
               }

    :param entity_list: [entity_1, ..., entity_n]
    :param name: string (name of the entity to be added)
    :param geom_object: geometric object of the entity
    :param mesh_sizes: dict
    """
    if mesh_sizes is not None:
        if name in mesh_sizes: mesh_size = mesh_sizes[name]
        elif 'mesh size' in mesh_sizes: mesh_size = mesh_sizes['mesh size']
        else: mesh_size = None
    else:
        mesh_size = None
    entity_list.append({'name': name,
                        'geometric object': geom_object,
                        'mesh size': mesh_size})

def create_entities_dict(name, face_entity_list, solid_entity_list, main_object=None):
    """
    Helper method for creating an entities dictionary.

    :param name: name of the collection of entities 
    :param face_entity_list: [face_entity_1, ..., face_entity_n]
    :param solid_entity_list: [solid_entity_1, ..., solid_entity_n]
    :param main_object: main object (usually a solid when the entity only has one)
    :return: entities_dict
    """
    entities_dict = {
            'name': name,
            'faces': face_entity_list,
            'solids': solid_entity_list,
            'main object': main_object
            }
    return entities_dict

def pick_faces_from_geometry(geom_object, face_picks, mesh_sizes=None):
    """
    Helper function for picking faces from a geometry object.
    The mesh sizes can be defined by providing the following dictionary:
    mesh_sizes={
                entity_name_1:mesh_size_for_entity_name_1,
                ...
                entity_name_n:mesh_size_for_entity_name_n,
                mesh size:mesh size if name is not in the dict
               }

    :param geom_object: FreeCAD geometric object where the faces are picked
    :face_picks: tuple('name of the face', int(face_number))
    :mesh_sizes: dict
    :return:
    """
    faces = []
    face_objects = geom_object.Shape.Faces
    for face_pick in face_picks:
        face_name = face_pick[0]
        face_number = face_pick[1]
        add_entity_in_list(faces, face_name, face_objects[face_number], mesh_sizes)
    return faces

def merge_entities_dicts(entities_dicts, name, default_mesh_size=None, add_prefixes={'solids':False, 'faces':True}):
    """ 
    This method merges all the entities_dicts and optionally prefixes the entity names with the 
    name of the entity. As default the solids are not prefixed but the faces are.

    :param entities_dicts: [entities_dict_1, ..., entities_dict_n]
    :name: string
    :default_mesh_size: float
    :add_prefixes: {'solids':bool, 'faces':bool}
    """
    entities_out = {'name': name}
    faces = []
    solids = []
    for d in entities_dicts:
        for face in d['faces']:
            if face['mesh size'] is None: face['mesh size'] = default_mesh_size
            if add_prefixes['faces']: face_name = d['name'] + '_' + face['name']
            else: face_name = face['name']
            add_entity_in_list(faces, face_name, face['geometric object'], {'mesh size':face['mesh size']})
        for solid in d['solids']:
            if add_prefixes['solids']: solid_name = d['name'] + '_' + solid['name']
            else: solid_name = solid['name']
            if solid['mesh size'] is None: solid['mesh size'] = default_mesh_size
            add_entity_in_list(solids, solid_name, solid['geometric object'], {'mesh size':solid['mesh size']})
    entities_out['faces'] = faces
    entities_out['solids'] = solids
    return entities_out

def get_solids_from_entities_dict(entities_dict):
    """
    Return a list of solids from entities dictionary.

    :param entities_dict: entities dictionary
    :return: [solid object 1, ..., solid object n]
    """
    solid_objects = [solid_dict_list['geometric object'] for solid_dict_list in entities_dict['solids']] 
    return solid_objects

def create_mesh_group_and_set_mesh_size(mesh_object, doc, name, mesh_size):
    """
    Creates mesh group with function ObjectsFem.makeMeshGroup.
    Adds property 'mesh_size' to created group and returns object.

    :param mesh_object: FreeCAD mesh object
    :param doc: FreeCAD document object.
    :param name: string
    :param mesh_size: float

    :return: MeshGroup object
    """
    # The third argument of makeMeshGroup is True, as we want to use labels,
    # not the names which cannot be changed.
    #
    # WARNING: No other object should have same label than this Mesh Group,
    # otherwise FreeCAD adds numbers to the end of the label to make it unique
    obj = ObjectsFem.makeMeshGroup(doc, mesh_object, True, name+'_group')
    obj.Label = name
    obj.addProperty('App::PropertyFloat', 'mesh_size')
    obj.mesh_size = mesh_size
    return obj

def merge_boundaries(mesh_object, compound_filter, doc, face_entity_dict, compound_face_names, face_name_list,
                     surface_objs, surface_objs_by_compound_face_names, surface_object=None):
    """
    If face in compound_faces is already added to surface:
        - renames surface if there was only one face in existing
        - removes face from existing surface and creates a new surface for merged face
    Creates new surface object (MeshGroup) for compound_faces if needed and surface_object is not given.

    :param mesh_object: FreeCAD mesh object
    :param compound_filter: FreeCAD compound filter
    :param doc: FreeCAD document object.
    :param face_entity_dict: dictionary
    :param compound_face_names: tuple containing compound face names in face
    :param face_name_list: list containing already handled face names
    :param surface_objs: list containing created surface objects same order as in face_name_list
    :param surface_objs_by_compound_face_names: dictionary (for checking if face needs to be merged)
    :param surface_object: None or already created surface object

    :return: tuple containing surface object and tuple containing filtered compound names
    """
    filtered_compound_faces = []
    for cface_name in compound_face_names:
        if cface_name in surface_objs_by_compound_face_names:
            surf_obj = surface_objs_by_compound_face_names[cface_name]
            old_face_name = surf_obj.Label
            new_face_name = '{}_{}'.format(old_face_name, face_entity_dict['name'])
            old_found_cface_names = surf_obj.References[0][1]
            filtered_old_found_cface_names = [cfname_i for cfname_i in old_found_cface_names if cfname_i != cface_name]
            if len(filtered_old_found_cface_names) == 0:
                # existing mesh object with new label
                surf_obj.Label = new_face_name
                # update face name in face_name_list
                index_found = face_name_list.index(old_face_name)
                face_name_list[index_found] = new_face_name
            else:
                # update references for existing mesh group
                surf_obj.References = [(compound_filter, tuple(filtered_old_found_cface_names))]
                # handle merged boundary
                if new_face_name in face_name_list:
                    # add merged boundary to existing mesh group
                    surface_index = face_name_list.index(new_face_name)
                    found_cface_names = surface_objs[surface_index].References[0][1]
                    surface_objs[surface_index].References = [(compound_filter, found_cface_names+tuple([cface_name]))]
                else:
                    # create new mesh group for merged boundary
                    surface_objs.append(create_mesh_group_and_set_mesh_size(mesh_object, doc, new_face_name,
                                                                            face_entity_dict['mesh size']))
                    surface_objs[-1].References = [(compound_filter, tuple([cface_name]))]
                    face_name_list.append(new_face_name)
                    surface_objs_by_compound_face_names[cface_name] = surface_objs[-1]
        else:
            filtered_compound_faces.append(cface_name)
            # create new mesh group only once if needed
            if surface_object is None:
                surface_object = create_mesh_group_and_set_mesh_size(mesh_object, doc, face_entity_dict['name'],
                                                                     face_entity_dict['mesh size'])
            surface_objs_by_compound_face_names[cface_name] = surface_object

    return surface_object, tuple(filtered_compound_faces)

def find_boundaries_with_entities_dict(mesh_object, compound_filter, entities_dict, doc):
    """
    For all faces in entities_dict, the same face in compound filter is added to a Mesh Group.
    All faces with same name in entities_dict are merged into one Mesh Group with the original name. 

    :param mesh_object: FreeCAD mesh object
    :param compound_filter: FreeCAD compound filter
    :param entities_dict: entities dictionary
    :param doc: FreeCAD document object.

    :return: list containing MeshGroup objects with mesh size.
    """
    surface_objs = []
    face_name_list = []
    surface_objs_by_cface_names = {}
    for num, face in enumerate(entities_dict['faces']):
        if face['name'] in face_name_list:
            # Old name, do not create new MeshGroup
            index_found = face_name_list.index(face['name'])
            found_cface_names = surface_objs[index_found].References[0][1]
            cface_names = find_compound_filter_boundaries(compound_filter, face['geometric object'])
            surface_obj, filtered_cface_names = merge_boundaries(mesh_object, compound_filter, doc, face,
                                                                 cface_names, face_name_list, surface_objs,
                                                                 surface_objs_by_cface_names,
                                                                 surface_object=surface_objs[index_found])
            if len(filtered_cface_names) > 0:
              surface_obj.References = [(compound_filter, found_cface_names+filtered_cface_names)]
        else:
            # New name, create new MeshGroup
            cface_names = find_compound_filter_boundaries(compound_filter, face['geometric object'])
            surface_obj, filtered_cface_names = merge_boundaries(mesh_object, compound_filter, doc, face,
                                                                 cface_names, face_name_list, surface_objs,
                                                                 surface_objs_by_cface_names, surface_object=None)
            if len(filtered_cface_names) > 0:  # new surface_obj is already created
                surface_obj.References = [(compound_filter, filtered_cface_names)]
                surface_objs.append(surface_obj)
                face_name_list.append(face['name'])
    return surface_objs

def find_bodies_with_entities_dict(mesh_object, compound_filter, entities_dict, doc, point_search=True):
    """
    For all solids in entities_dict, the same solid in compound filter is added to a Mesh Group.
    All solids with same name in entities_dict are merged into one Mesh Group with the original name.

    :param mesh_object: FreeCAD mesh object
    :param compound_filter: FreeCAD compound filter
    :param entities_dict: entities dictionary
    :param doc: FreeCAD document object.
    :param point_search: bool

    :return: list containing MeshGroup objects with mesh size.
    """
    solid_objs = []
    solid_name_list = []
    for num, solid in enumerate(entities_dict['solids']):
        if solid['name'] in solid_name_list:
            # Old name, do not create new MeshGroup
            index_found = solid_name_list.index(solid['name'])
            found_csolid_names = solid_objs[index_found].References[0][1]
            csolid_names = find_compound_filter_solids(compound_filter, solid['geometric object'].Shape, point_search)
            found_csolid_names = found_csolid_names + csolid_names
            solid_objs[index_found].References = [(compound_filter, found_csolid_names)]
        else:
            # New name, create new MeshGroup
            solid_objs.append(create_mesh_group_and_set_mesh_size(mesh_object, doc, solid['name'], solid['mesh size']))
            csolid_names = find_compound_filter_solids(compound_filter, solid['geometric object'].Shape, point_search)
            solid_objs[-1].References = [(compound_filter, csolid_names)]
            solid_name_list.append(solid['name'])
    return solid_objs

def define_mesh_sizes_with_mesh_groups(mesh_object, mesh_group_list, doc, ignore_list=None):
    """
    Meshregions are needed to have regionwise mesh density parameters.
    The mesh element length is the third parameter given in makeMeshRegion.
    Each mesh group in mesh_group_list needs to know its
    mesh size (created with function :meth:`create_mesh_group_and_set_mesh_size`).

    :param mesh_object: FreeCAD mesh object
    :param mesh_group_list: list containing MeshGroups
    :param doc: FreeCAD document object.
    :param ignore_list: None or list containing solid names which mesh size is not defined.
    """
    if ignore_list is None:
        ignore_list = []
    for mesh_group in mesh_group_list:
        if mesh_group.Label not in ignore_list:
            solid_obj = ObjectsFem.makeMeshRegion(doc, mesh_object, mesh_group.mesh_size, mesh_group.Name+'_region')
            solid_obj.References = [(mesh_group.References[0][0], mesh_group.References[0][1])]

def define_mesh_sizes(mesh_object, compound_filter, entities_dict, doc, point_search=True, ignore_list=None):
    """
    Meshregions are needed to have regionwise mesh density parameters. 
    The mesh element length is the third parameter given in makeMeshRegion. 

    :param mesh_object: FreeCAD mesh object
    :param compound_filter: FreeCAD compound filter
    :param entities_dict: entities dictionary
    :param doc: FreeCAD document object.
    :param point_search: bool
    :param ignore_list: None or list containing solid names which mesh size is not defined.
    """
    if ignore_list is None:
        ignore_list = []
    solid_objs = []
    solid_name_list = []
    for num, solid in enumerate(entities_dict['solids']):
        if solid['name'] in ignore_list:
            continue
        if solid['name'] in solid_name_list:
            # Old name, do not create new MeshGroup
            index_found = solid_name_list.index(solid['name'])
            found_csolid_names = solid_objs[index_found].References[0][1]
            csolid_names = find_compound_filter_solids(compound_filter, solid['geometric object'].Shape, point_search)
            solid_objs[index_found].References = [(compound_filter, found_csolid_names+csolid_names)]
        else:
            # New name, create new MeshGroup
            solid_objs.append(ObjectsFem.makeMeshRegion(doc, mesh_object, solid['mesh size'], solid['name']+'_region'))
            csolid_names = find_compound_filter_solids(compound_filter, solid['geometric object'].Shape, point_search)
            solid_objs[-1].References = [(compound_filter, csolid_names)]
            solid_name_list.append(solid['name'])


