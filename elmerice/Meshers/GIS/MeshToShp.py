#!/usr/bin/python
#==========================================
#
# FILE: MeshToShp.py 
# USAGE : python MeshToShp.py  [-h] -d <inputdir>
# DESCRIPTION:  Create shapefiles from elmer mesh
#
# BUGS: ---
#
# AUTHOR:   F. Gillet-Chaulet
# ORGANIZATION: IGE(CNRS-France)
#
# VERSION: V1
# CREATED: 2020-05-02
# MODIFIED:
#
#==========================================
import sys, getopt,os

def main(argv):
   import shapefile
   
   found_d=False
   bc_only=False
   try:
      opts, args = getopt.getopt(argv,"hd:",["bc"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-d"):
         dir_name= arg
         found_d=True
      elif opt in ("--bc"):
         bc_only=True

   if not found_d:
      print('missing mandatory mesh dir name')
      usage()
      sys.exit()

   nodefile=os.path.join(dir_name, 'mesh.nodes') 
   vertices = {}
   with open(nodefile) as fin:
     for  line in fin:
         (j,k,x,y,z) = line.split()
         vertices[int(j)] = [float(x), float(y)]
 
   outputdir="{}_shp".format(dir_name)
   try: 
       # Create target Directory
       os.makedirs(outputdir,exist_ok=True)

   except OSError as error: 
    print("Directory {} can not be created" .format(outputdir)) 
  
   ## bc as polylines
   bcfile=os.path.join(outputdir, 'boundaries')
   shp=shapefile.Writer(bcfile, shapefile.POLYLINE)
   shp.field('enum', 'N')
   shp.field('etype', 'N')
   shp.field('BCId', 'N')

   boundaryfile = os.path.join(dir_name, 'mesh.boundary')
   with open(boundaryfile) as fin:
       for line in fin:
         l = line.split()
         enum,bc,etype=int(l[0]),int(l[1]),int(l[4])
         nv=etype%100
         evertices = [vertices[int(i)] for i in l[5:5+nv]]
         shp.line([evertices])
         shp.record(enum,etype,bc)

   shp.close()

   if not bc_only :
      ## elemnts as polygons
      bcfile=os.path.join(outputdir, 'elements')
      shp=shapefile.Writer(bcfile, shapefile.POLYGON)
      shp.field('enum', 'N')
      shp.field('etype', 'N')
      shp.field('BodyId', 'N')

      Efile = os.path.join(dir_name, 'mesh.elements')
      with open(Efile) as fin:
          for line in fin:
            l = line.split()
            enum,bd,etype=int(l[0]),int(l[1]),int(l[2])
            nv=etype%100
            evertices = [vertices[int(i)] for i in l[3:3+nv]]

            # As its elements so no hole
            # should we check the rotation order?
            #  seems not
            # this would be the solution wiyh shapely
            #polygon = shapely.geometry.Polygon(evertices)
            #if not polygon.exterior.is_ccw:
            #  evertices=evertices[::-1]

            # spyshp autoimatically add lastpt=firstpt
            # to close polygons
            shp.poly([evertices])
            shp.record(enum,etype,bd)

      shp.close()

   print("Shapefiles for Elmer mesh have been created")
   print("You can define the projection with gdal tools")
   print("  gdalsrsinfo  -o wkt \"EPSG:XYZW\" > {}/elements.prj".format(outputdir))
   print("  gdalsrsinfo  -o wkt \"EPSG:XYZW\" > {}/boundaries.prj".format(outputdir))

def usage():
   print('usage: MeshToShp.py  [-h] -d <inputfile>')
   print('options:')
   print('   -h [print help]')
   print('   -d <mesh dir. name>')
   print('   --bc [output only boundary elements]')


if __name__ == "__main__":
   main(sys.argv[1:])
