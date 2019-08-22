doc = App.newDocument('cube geometry')
import sys
import os
import time
import math
import FreeCAD

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import FreeCADBatchFEMTools


def cut_xy_plane(part, cut_len, cut_width, cut_height):
    doc.recompute()
    doc.recompute()
    reduced_name = part.Name + '_xy'
    tool_box = doc.addObject("Part::Box", "CutBox" + reduced_name + '_obj')
    tool_box.Length = cut_len
    tool_box.Width = cut_width
    tool_box.Height = cut_height
    c_x = part.Shape.BoundBox.Center.x
    c_y = part.Shape.BoundBox.Center.y
    tool_box.Placement = App.Placement(App.Vector(c_x-cut_len/2., c_y-cut_width/2., -cut_height),
                                       App.Rotation(App.Vector(0, 0, 1), 0))

    cut_part = doc.addObject("Part::Cut", reduced_name + '_obj')
    cut_part.Base = part
    cut_part.Tool = tool_box
    return cut_part

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

def create_air_geometry(entities_list, length, width, height, shift, mesh_size):
    """
    Creates air around created parts in part_list.

    :param entities_list: A list containing entities dicts.
    :param length: A float.
    :param width: A float.
    :param height: A float.
    :param shift: A float.
    :param mesh_size: A float or None.

    :return: Air part object.
    """
    bounding_box_part = doc.addObject("Part::Box", "air_obj")
    bounding_box_part.Length = length
    bounding_box_part.Width = width
    bounding_box_part.Height = height
    bounding_box_part.Placement = App.Placement(App.Vector(-shift, -shift, 0),
                                                App.Rotation(App.Vector(0, 0, 1), 0))
    doc.recompute()

    solid_objects = [bounding_box_part] + [solid['main object'] for solid in entities_list]
    air_part = FreeCADBatchFEMTools.create_xor_object(solid_objects, doc)
    doc.recompute()
    # remove already created faces
    faces_in_symmetry_plane = FreeCADBatchFEMTools.faces_with_vertices_in_symmetry_plane(air_part.Shape.Faces,
                                                                                         plane='xy')
    used_faces = []
    for part in entities_list:
        used_faces.extend([face_entity['geometric object'] for face_entity in part['faces']])
    FreeCADBatchFEMTools.remove_compare_faces_from_list(used_faces, faces_in_symmetry_plane)
    faces = []
    for face in faces_in_symmetry_plane:
        FreeCADBatchFEMTools.add_entity_in_list(faces, 'beta0', face)
    solids = []
    FreeCADBatchFEMTools.add_entity_in_list(solids, 'air', air_part, {'mesh size': mesh_size})
    entities_dict = FreeCADBatchFEMTools.create_entities_dict('air', faces, solids, main_object=air_part)
    doc.recompute()

    return entities_dict

def create_sphere_geometry(nof_spheres, mesh_size, diameter, distance):
    """
    Create sphere geometry. Spheres are moved down 10% and cut from xy plane
    because caused problems with air otherwise.

    :param nof_spheres: Number of spheres.
    :param mesh_size: Mesh size of spheres.
    :param diameter: Diameter of spheres.
    :param distance: Distance between spheres.

    :return: A list containing entities dictionaries.
    """
    sphere_entities_dict_list = []
    r = diameter/2.0
    shift = diameter/10.0
    n = int(math.ceil(math.sqrt(nof_spheres)))
    bbox_len = n*(diameter+distance)
    counter = 0
    for i in range(n):
        for j in range(n):
            sphere_name = 'sphere{:04d}'.format(counter + 1)
            sphere = doc.addObject('Part::Sphere', sphere_name + '_obj')
            sphere.Radius = r
            sphere.Placement = App.Placement(App.Vector(r + i*(diameter+distance), r + j*(diameter+distance), r-shift),
                                             App.Rotation(App.Vector(0, 0, 1), 0))
            doc.recompute()
            sphere = cut_xy_plane(sphere, bbox_len, bbox_len, diameter)
            doc.recompute()
            # create entities dict
            face_picks = [('alpha1', 0), ('beta0', 1)]
            faces = FreeCADBatchFEMTools.pick_faces_from_geometry(sphere, face_picks)
            solids = []
            FreeCADBatchFEMTools.add_entity_in_list(solids, sphere_name, sphere, {'mesh size': mesh_size})
            entities_dict = FreeCADBatchFEMTools.create_entities_dict(sphere_name, faces, solids, sphere)
            sphere_entities_dict_list.append(entities_dict)
            doc.recompute()
            counter += 1
            if counter == nof_spheres:
                return sphere_entities_dict_list

def create_cube_geometry(nof_cubes, mesh_size, length, width, height, cube_distance):
    """
    Create cube geometry.

    :param nof_cubes: Number of cubes.
    :param mesh_size: Mesh size of cubes.
    :param length: Length of cubes.
    :param width: Width of cubes.
    :param height: Height of cubes.
    :param cube_distance: Distance between cubes.

    :return: A list containing entities dictionaries.
    """
    cube_entities_dict_list = []
    n = int(math.ceil(math.sqrt(nof_cubes)))
    counter = 0
    for i in range(n):
        for j in range(n):
            cube_name = 'cube{:04d}'.format(counter + 1)
            cube = doc.addObject('Part::Box', cube_name + '_obj')
            cube.Length = length
            cube.Width = width
            cube.Height = height
            cube.Placement = App.Placement(App.Vector(i * (length+cube_distance), j * (width+cube_distance), 0),
                                           App.Rotation(App.Vector(0, 0, 1), 0))
            doc.recompute()
            # create entities dict
            face_picks = [('alpha0', 5), ('alpha1', 0), ('beta0', 4), ('beta1', 2), ('gamma0', 3), ('gamma1', 1)]
            faces = FreeCADBatchFEMTools.pick_faces_from_geometry(cube, face_picks)
            solids = []
            FreeCADBatchFEMTools.add_entity_in_list(solids, cube_name, cube, {'mesh size': mesh_size})
            entities_dict = FreeCADBatchFEMTools.create_entities_dict(cube_name, faces, solids, cube)
            cube_entities_dict_list.append(entities_dict)
            doc.recompute()
            counter += 1
            if counter == nof_cubes:
                return cube_entities_dict_list


def create_geometry(nof_cubes, mesh_size, find_boundaries_from_entities_dict, find_solids_from_entities_dict, directory,
                    add_air, point_search, use_spheres, mesh_size_air):
    """
    Creates geometry fof cubes and exports unv file.

    :param nof_cubes: Number of cubes/spheres created.
    :param mesh_size: Mesh size for cubes/spheres.
    :param find_boundaries_from_entities_dict: 'True' or 'False'.
    :param find_solids_from_entities_dict: 'True' or 'False'.
    :param directory: Path to directory of this script.
    :param add_air: A boolean. Is air added to geometry.
    :param point_search: A boolean.
    :param use_spheres: A boolean.
    :param mesh_size_air: A float.

    :return: A list containing tuples (name, execution time).
    """
    total_start_time = time.time()
    entities_list = []
    return_time_list = []
    cube_edge_length, cube_distance = 200, 20
    mesh_size_max = max(mesh_size, mesh_size_air)
    if use_spheres:
        print_line('Creating sphere geometry...', total_start_time)
        entities_list.extend(create_sphere_geometry(nof_cubes, mesh_size, cube_edge_length, cube_distance))
    else:
        print_line('Creating cube geometry...', total_start_time)
        entities_list.extend(create_cube_geometry(nof_cubes, mesh_size, cube_edge_length, cube_edge_length,
                                                  cube_edge_length, cube_distance))
    if add_air:
        n = int(math.ceil(math.sqrt(nof_cubes)))
        air_edge_length = n*(cube_edge_length+cube_distance) + cube_distance
        air_height = cube_edge_length+2*cube_distance
        print_line("Creating air geometry...", total_start_time)
        entities_list.append(create_air_geometry(entities_list, air_edge_length, air_edge_length, air_height,
                                                 cube_distance, mesh_size_air))
    entities_list.reverse()  # give air first to get mesh sizes correctly
    print_line("Merging entities dictionaries...", total_start_time)
    entities_dict = FreeCADBatchFEMTools.merge_entities_dicts(entities_list, 'All', default_mesh_size=mesh_size_max,
                                                              add_prefixes={'solids': False, 'faces': True})
    print_line("Getting solids from entities dictionaries...", total_start_time)
    solid_objects = FreeCADBatchFEMTools.get_solids_from_entities_dict(entities_dict)
    print_line("Creating mesh object and compound filter...", total_start_time)
    mesh_object, compound_filter = FreeCADBatchFEMTools.create_mesh_object_and_compound_filter(solid_objects,
                                                                                               mesh_size_max, doc)
    if find_boundaries_from_entities_dict.lower() == 'true':
        print_line("Finding boundaries...", total_start_time)
        start_time = time.time()
        FreeCADBatchFEMTools.find_boundaries_with_entities_dict(mesh_object, compound_filter,
                                                                entities_dict, doc)
        return_time_list.append(('-Find boundaries:', time.time() - start_time))
    if find_solids_from_entities_dict.lower() == 'true':
        print_line("Finding bodies...", total_start_time)
        start_time = time.time()
        body_mesh_groups = FreeCADBatchFEMTools.find_bodies_with_entities_dict(mesh_object, compound_filter,
                                                                               entities_dict, doc, point_search)
        return_time_list.append(('-Find solids:', time.time() - start_time))
    print_line("Defining mesh sizes...", total_start_time)
    start_time = time.time()
    # ignore air from mesh definition, mesh_object forces air to largest mesh
    if find_solids_from_entities_dict.lower() == 'true':
        FreeCADBatchFEMTools.define_mesh_sizes_with_mesh_groups(mesh_object, body_mesh_groups, doc, ignore_list=['air'])
    else:  # in this case mesh sizes needs to be found
        FreeCADBatchFEMTools.define_mesh_sizes(mesh_object, compound_filter, entities_dict, doc, point_search,
                                               ignore_list=['air'])
    return_time_list.append(('-Define mesh sizes:', time.time() - start_time))
    FreeCADBatchFEMTools.fit_view()
    start_time = time.time()
    print_line("Creating mesh...", total_start_time)
    FreeCADBatchFEMTools.create_mesh(mesh_object)
    return_time_list.append(('-Create mesh:', time.time() - start_time))
    print_line("Exporting unv...", total_start_time)
    FreeCADBatchFEMTools.export_unv(os.path.join(directory, 'cubemeshtest.unv'), mesh_object)
    print_line("Geometry done", total_start_time)
    return return_time_list + [('-Total time:', time.time() - total_start_time)]


script_directory = os.path.dirname(__file__)

with open(os.path.join(script_directory, 'cubemeshtestparameters.txt')) as f:
    [number_of_cubes, cube_mesh_size, find_boundaries, find_solids, create_air, air_mesh_size,
     find_with_points, create_spheres, append_file] = f.read().split()

create_air = create_air.lower() == 'true'
find_with_points = find_with_points.lower() == 'true'
create_spheres = create_spheres.lower() == 'true'
if air_mesh_size == 'None':
    air_mesh_size = cube_mesh_size
try:
    elapsed_times = create_geometry(int(number_of_cubes), float(cube_mesh_size), find_boundaries, find_solids,
                                    script_directory, create_air, find_with_points, create_spheres, float(air_mesh_size))
except Exception:
    import traceback
    print(str(traceback.format_exc()))
else:
    if create_air:
        info_line = '\nexecution times with {} cubes (air: {}, spheres: {}, fb: {}, '.format(number_of_cubes, create_air,
                                                                                             create_spheres, find_boundaries)
        info_line += 'fs: {}, mesh_size: {}, air_mesh_size: {}, point_search: {}):'.format(find_solids, cube_mesh_size,
                                                                                           air_mesh_size, find_with_points)
    else:
        info_line = '\nexecution times with {} cubes (air: {}, spheres: {}, fb: {}, '.format(number_of_cubes, create_air,
                                                                                             create_spheres, find_boundaries)
        info_line += 'fs: {}, mesh_size: {}, point_search: {}):'.format(find_solids, cube_mesh_size, find_with_points)
    print(info_line)
    for time_tuple in elapsed_times:
        print(time_tuple[0], round(time_tuple[1], 1))
    sys.stdout.flush()
    if append_file.lower() == 'true':
        with open(os.path.join(script_directory, 'cubemeshtestexecutiontimes.txt'), 'a') as f:
            f.write(info_line + '\n')
            for time_tuple in elapsed_times:
                f.write('{} {}\n'.format(time_tuple[0], round(time_tuple[1], 1)))
            f.write('\n')
if not FreeCAD.GuiUp:
  exit()
