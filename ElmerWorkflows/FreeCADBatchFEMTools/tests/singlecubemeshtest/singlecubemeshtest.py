import argparse
import os
import subprocess

description = 'Runs script singlecubemeshtest_freecadscript.py with FreeCAD.'

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-fe', '--freecad-executable', type=str, default="FreeCAD", help='give the path to FreeCAD executable')

args = parser.parse_args()

directory = os.path.dirname(os.path.realpath(__file__))
freecadscript_name = os.path.join(directory, 'singlecubemeshtest_freecadscript.py')
if not os.path.isfile(freecadscript_name):
    print("singlecubemeshtest_freecadscript.py does not exist, check that you are in correct directory")
else:
    try:
        p = subprocess.Popen([args.freecad_executable, '-c', freecadscript_name])
        p.communicate()
    except Exception:
        print("Running FreeCAD failed!!! Try to give the correct FreeCAD executable as an argument (--freecad-executable, -fe)")
