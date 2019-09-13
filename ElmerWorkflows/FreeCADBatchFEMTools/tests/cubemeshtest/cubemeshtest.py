import argparse
import os
import subprocess

description = 'Writes parameters to file cubemeshtestparameters.txt and runs cubemeshtest_freecadscript.py with FreeCAD.'
description += ' Edge length of cube is 200 (diameter 200 with spheres). Distance between two cubes (spheres) is 20.'

ps_help = 'Is points used for finding boundaries and solids (section and common search otherwise)'
af_help = 'Appends execution times to file cubemeshtestexecutiontimes.txt if used'

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-n', '--number_of_cubes', type=int, required=True)
parser.add_argument('-m', '--cube_mesh_size', type=float, default=50.0)
parser.add_argument('-fb', '--find_boundaries', action='store_true')
parser.add_argument('-fs', '--find_solids', action='store_true')
parser.add_argument('-ca', '--create_air', action='store_true', help='Is air created around cubes')
parser.add_argument('-am', '--air_mesh_size', type=float, help='If not given cube_mesh_size is used')
parser.add_argument('-ps', '--point_search', action='store_true', help=ps_help)
parser.add_argument('-s', '--spheres', action='store_true', help='Use spheres instead of cubes')
parser.add_argument('-af', '--append_file', action='store_true', help=af_help)
parser.add_argument('-fe', '--freecad-executable', type=str, default="FreeCAD", help='give the path to FreeCAD executable')

args = parser.parse_args()

directory = os.path.dirname(os.path.realpath(__file__))
freecadscript_name = os.path.join(directory, 'cubemeshtest_freecadscript.py')
if not os.path.isfile(freecadscript_name):
    print("cubemeshtest_freecadscript.py does not exist, check that you are in correct directory")
else:
    with open(os.path.join(directory, 'cubemeshtestparameters.txt'), 'w') as f:
        f.write('{} {} {} {} {} {} {} {} {}'.format(args.number_of_cubes, args.cube_mesh_size, args.find_boundaries,
                                                    args.find_solids, args.create_air, args.air_mesh_size,
                                                    args.point_search, args.spheres, args.append_file))
    try:
        p = subprocess.Popen([args.freecad_executable, '-c', freecadscript_name])
        p.communicate()
    except:
        print("Running FreeCAD failed!!! Try to give the correct FreeCAD executable as an argument (--freecad-executable, -fe)")
        exit()
