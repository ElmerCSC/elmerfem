doc = App.newDocument('single cube geometry')
import sys
import os
import time
import Part
import FreeCAD

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import FreeCADBatchFEMTools


def print_line(line_str, start_time):
    """
    Prints elapsed time and flushes stdout.

    :param line_str: A string.
    :param start_time: A float.
    """
    elapsed_time = round(time.time() - start_time, 1)
    print('Elapsed time: {} seconds'.format(elapsed_time))
    print(line_str)
    sys.stdout.flush()
    return elapsed_time

def import_step_part(doc, path, name):
    """
    Imports part from step file and returns it.

    :param doc: FreeCAD document.
    :param path: Path to step file.
    :param name: A string.

    :return: FreeCAD document object 'Part::Feature'.
    """
    part = doc.addObject('Part::Feature', name + '_obj')
    part.Shape = Part.read(path)
    doc.recompute()
    return part

def get_cube_entity_dict(cube_name, mesh_size, directory, x_shift):
    """
    Import cube and create entity dictionary.

    :param cube_name: A string.
    :param mesh_size: Mesh size of cube.
    :param directory: Path to directory of this script.
    :param x_shift: Numeric.

    :return: A list containing entities dictionaries.
    """
    cube = import_step_part(doc, os.path.join(directory, '..', 'stepfiles', 'cube200.step'), cube_name)
    cube.Placement = App.Placement(App.Vector(x_shift, 0, 0), App.Rotation(App.Vector(0, 0, 1), 0))
    # create entities dict
    face_picks = [('alpha0', 5), ('alpha1', 0), ('beta0', 4), ('beta1', 2), ('gamma0', 3), ('gamma1', 1)]
    faces = FreeCADBatchFEMTools.pick_faces_from_geometry(cube, face_picks)
    solids = []
    FreeCADBatchFEMTools.add_entity_in_list(solids, cube_name, cube, {'mesh size': mesh_size})
    return FreeCADBatchFEMTools.create_entities_dict(cube_name, faces, solids, cube)

def create_geometry(directory):
    """
    Imports cube two times and sets them touching each other. Exports unv file with separate boundaries.

    :param directory: Path to directory of this script.

    :return: A list containing tuples (name, execution time).
    """
    total_start_time = time.time()
    separate_bounds = True
    entities_list = []
    return_time_list = []
    mesh_size = 20
    mesh_size_max = mesh_size*2

    print_line('Creating cube geometry...', total_start_time)
    entities_list.append(get_cube_entity_dict('cube1', mesh_size, directory, 0))
    entities_list.append(get_cube_entity_dict('cube2', mesh_size*2., directory, 200))
    print_line("Merging entities dictionaries...", total_start_time)
    entities_dict = FreeCADBatchFEMTools.merge_entities_dicts(entities_list, 'All', default_mesh_size=mesh_size_max,
                                                              add_prefixes={'solids': False, 'faces': True})
    print_line("Getting solids from entities dictionaries...", total_start_time)
    solid_objects = FreeCADBatchFEMTools.get_solids_from_entities_dict(entities_dict)
    print_line("Creating mesh object and compound filter...", total_start_time)
    [mesh_object,
     compound_filter] = FreeCADBatchFEMTools.create_mesh_object_and_compound_filter(solid_objects, mesh_size_max, doc,
                                                                                    separate_boundaries=separate_bounds)
    print_line("Finding boundaries...", total_start_time)
    start_time = time.time()
    FreeCADBatchFEMTools.find_boundaries_with_entities_dict(mesh_object, compound_filter, entities_dict, doc,
                                                            separate_boundaries=separate_bounds)
    return_time_list.append(('-Find boundaries:', time.time() - start_time))
    print_line("Finding bodies...", total_start_time)
    start_time = time.time()
    body_mesh_groups = FreeCADBatchFEMTools.find_bodies_with_entities_dict(mesh_object, compound_filter,
                                                                           entities_dict, doc)
    return_time_list.append(('-Find solids:', time.time() - start_time))
    print_line("Defining mesh sizes...", total_start_time)
    start_time = time.time()
    FreeCADBatchFEMTools.define_mesh_sizes_with_mesh_groups(mesh_object, body_mesh_groups, doc)
    return_time_list.append(('-Define mesh sizes:', time.time() - start_time))
    FreeCADBatchFEMTools.fit_view()
    start_time = time.time()
    print_line("Creating mesh...", total_start_time)
    FreeCADBatchFEMTools.create_mesh(mesh_object, directory=directory)
    return_time_list.append(('-Create mesh:', time.time() - start_time))
    print_line("Exporting unv...", total_start_time)
    FreeCADBatchFEMTools.export_unv(os.path.join(directory, 'separate_boundaries_doublecubemeshtest.unv'), mesh_object)
    print_line("Geometry done", total_start_time)
    return return_time_list + [('-Total time:', time.time() - total_start_time)]


script_directory = os.path.dirname(__file__)

try:
    elapsed_times = create_geometry(script_directory)
except Exception:
    import traceback
    print(str(traceback.format_exc()))
else:
    print('\nexecution times with 1 cube:')
    for time_tuple in elapsed_times:
        print(time_tuple[0], round(time_tuple[1], 1))
    sys.stdout.flush()

if not FreeCAD.GuiUp:
    exit()
