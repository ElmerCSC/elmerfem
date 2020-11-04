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

  Authors: Eelis Takala, Sami Rannikko
  Emails:  eelis.takala@gmail.com
  Address: Trafotek Oy
           Kaarinantie 700
           20540 Turku
           Finland

  Original Date: April 2020
"""
import os
import shutil


def parse_geo_id_list_string_to_list(geo_list_string):
    """
    Collects ids from geo string.

    :param geo_list_string: A string e.g. '{15, 16};'

    :return: A list e.g ['15', '16'].
    """
    id_list = []
    for id_str in geo_list_string.replace(';', '').replace('{', '').replace('}', '').split(','):
        id_list.append(id_str.strip())
    return id_list


def collect_geometry_ids_from_geo_file(geo_file_path):
    """
    Collects geometry ids from geo file to dictionary.

    :param geo_file_path: Path to geo file.

    :return: A dictionary.
    """
    geom_id_dict = {}
    with open(geo_file_path, 'r') as geo_file:
        for line in geo_file:
            if line.startswith('Physical Volume(') or line.startswith('Physical Surface('):
                left_hand_side, right_hand_side = line.split('=')  # Physical Volume("SolidName") = {2, 3};
                volume_or_surface_name = left_hand_side.split('"')[1]
                geom_id_dict[volume_or_surface_name] = parse_geo_id_list_string_to_list(right_hand_side)
    return geom_id_dict


def _get_transfinite_line_geo_file_line(line_param_dict):
    """
    Returns transfinite line geo file line 'Transfinite Line {59, 60} = 9 Using Progression 1;'.

    :param line_param_dict: A dictionary containing transfinite line params.

    :return: A string.
    """
    if line_param_dict.get('comment', ''):
        comment = '  // {}'.format(line_param_dict['comment'])
    else:
        comment = ''
    return 'Transfinite Line {{{}}} = {} Using Progression {};{}\n'.format(', '.join(line_param_dict['lines']),
                                                                           line_param_dict['points'],
                                                                           line_param_dict['progression'], comment)


def _get_transfinite_surface_geo_file_line(surface_name, transfinite_param_dict, geometry_id_dict,
                                           exact_surface_equality):
    """
    Returns transfinite surface geo file line e.g. 'Transfinite Surface {42, 44, 43, 41, 131, 130} Right;'.

    :param surface_name: A string.
    :param transfinite_param_dict:  A dictionary containing transfinite parameters.
    :param geometry_id_dict: A dictionary (names to geometry ids)
    :param exact_surface_equality: A boolean.

    :return: A string.
    """
    geo_file_surface_name = None
    if exact_surface_equality:
        geo_file_surface_name = surface_name
    else:
        for geo_file_surface_name_i in geometry_id_dict:
            if surface_name in geo_file_surface_name_i:  # 'A1_alpha1' in 'A1_alpha1_A2_alpha0'
                geo_file_surface_name = geo_file_surface_name_i
                break
    if geo_file_surface_name is None or geo_file_surface_name not in geometry_id_dict:
        raise ValueError('Surface {} not found from geo file'.format(surface_name))
    surface_id_str_list = geometry_id_dict[geo_file_surface_name]
    direction_key = '{}_direction'.format(surface_name)
    if transfinite_param_dict.get(direction_key, ''):
        return 'Transfinite Surface {{{}}} {};  // {}\n'.format(', '.join(surface_id_str_list),
                                                                transfinite_param_dict[direction_key],
                                                                geo_file_surface_name)
    else:
        return 'Transfinite Surface {{{}}};  // {}\n'.format(', '.join(surface_id_str_list), geo_file_surface_name)


def add_transfinite_lines_to_geo_file(directory, transfinite_param_list, file_name='shape2mesh.geo',
                                      exact_surface_equality=False):
    """
    Adds transfinite lines to file file_name.

    :param directory: Gmsh temp file location.
    :param transfinite_param_list: A list containing dictionaries {'volume': 'name',
                                                                   'surface_list': [s_name, s_name2]}.
    :param file_name: A string.
    :param exact_surface_equality: A boolean. When given surfaces are compared from geo file name use exact names.
                                   Default values is False because surfaces can be merged (e.g. 'A1_alpha0_A2_alpha1').

    """
    geo_path = os.path.join(directory, file_name)
    geo_path_cp = os.path.join(directory, 'shape2mesh_cp.geo')
    shutil.copy(geo_path, geo_path_cp)
    geom_id_dict = collect_geometry_ids_from_geo_file(geo_path)
    with open(geo_path, 'w') as geo_file:
        with open(geo_path_cp, 'r') as original_geo_file_cp:
            lines_added = False
            for i, line in enumerate(original_geo_file_cp):
                if i == 1:
                    geo_file.write('General.ExpertMode = 1;  // Allow all mesh algorithms for transfinite\n')
                if not lines_added and line.startswith('// Characteristic Length'):
                    lines_added = True
                    for p_dict in transfinite_param_list:
                        geo_file.write('// Transfinite {}\n'.format(p_dict['volume']))
                        for line_p_dict in p_dict.get('line_params', []):
                            geo_file.write(_get_transfinite_line_geo_file_line(line_p_dict))
                        for surface_name in p_dict['surface_list']:
                            geo_file.write(_get_transfinite_surface_geo_file_line(surface_name, p_dict, geom_id_dict,
                                                                                  exact_surface_equality))
                        geo_file.write('Transfinite Volume {{{}}};\n'.format(', '.join(geom_id_dict[p_dict['volume']])))
                    geo_file.write('\n')
                geo_file.write(line)
    os.remove(geo_path_cp)

