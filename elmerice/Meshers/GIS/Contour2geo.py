#!/usr/bin/python
#========================================== 
#
# FILE: Contour2geo.py
# USAGE : python Contour2geo.py -r res [-h] [-i <inputfile>] [-o <outputfile>] [--spline] ')
# DESCRIPTION:  Create a geometry file (.geo) for Gmsh from a closed contour
#
# BUGS: ---
#
# AUTHOR:   F. Gillet-Chaulet
# ORGANIZATION: IGE(CNRS-France)
#
# VERSION: V1 
# CREATED: 2020-05-02
# MODIFIED: 
#  * 2020-05-19: prescribe mesh size using uniform background fied
#
#========================================== 
import sys, getopt

def main(argv):
   import numpy  as np

   inputfile = 'Contour.txt'
   outputfile = 'Contour.geo'
   spline = False

   found_r=False
   try:
      opts, args = getopt.getopt(argv,"hi:o:r:",["spline","gmsh3"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-r"):
         el_size= arg
         found_r=True
      elif opt in ("-i"):
         inputfile = arg
      elif opt in ("-o"):
         outputfile = arg
      elif opt in ("--spline"):
         spline = True
   
   if not found_r:
      print('missing mandatory argument -r')
      usage()
      sys.exit(2)
    
   print('Input file is : %s'%(inputfile))
   print('Output file is :  %s'%(outputfile))
   if (spline):
      print('make contour using splines')
   else :
      print('make contour using compound lines')

###############
###############
   if inputfile.lower().endswith('.shp'):
      import shapefile
      sf = shapefile.Reader(inputfile)

      #number of layers
      nlayers=len(sf)
      print('found {0} features in the shapefile'.format(nlayers))

      # index array for BC conditions
      index=np.arange(nlayers, dtype=np.uint8)

      # get shape records
      shapeRecs = sf.shapeRecords()
      c=0
      for r in shapeRecs:
          # BC given by attribute
          if hasattr(r.record, 'BC'):
              index[c]=r.record.BC
          else:
          # or shape order
              index[c]=c+1
          c=c+1

      print('BC ordering by feature: ')
      print(index)

      # sort shapes by BC
      sarray=np.argsort(index)
      
      # start with first shape
      lnum=0
      shape = shapeRecs[sarray[lnum]]

      # it works only is last point is first point of next shape
      # check this
      if ((lnum+1)>(nlayers-1)) :
          nextshape=shapeRecs[sarray[0]]
      else:
          nextshape=shapeRecs[sarray[lnum+1]]
      if (shape.shape.points[-1] != nextshape.shape.points[0]):
          print('sorry only works if last point is first of next feature')
          sys.exit(2)

      # take all points except last
      npp=len(shape.shape.points)-1
      pts=shape.shape.points[0:npp]
      # there will be npp lines for this BC
      lines=[npp]

      for i in range(1,len(sarray)):
        lnum=i
        shape = shapeRecs[sarray[lnum]]
        # it works only is last point is first point of next shape
        # check this
        if ((lnum+1)>(nlayers-1)) :
          nextshape=shapeRecs[sarray[0]]
        else:
          nextshape=shapeRecs[sarray[lnum+1]]
        if (shape.shape.points[-1] != nextshape.shape.points[0]):
          print('sorry only works if last point is first of next feature')
          sys.exit(2)

        npp=len(shape.shape.points)-1
        pts.extend(shape.shape.points[0:npp])
        lines.append(npp)

      Contour = np.array(pts)
      x = Contour[:,0]
      y = Contour[:,1]
      Npt = len(x)

   else:
      # ascii files
      Contour = np.loadtxt(inputfile)
      x = Contour[:,0]
      y = Contour[:,1]

      if x[0]==x[-1] and y[0]==y[-1]:
        #print('Same first and last points in contour file')
        Npt = len(x)-1
      else:
        Npt = len(x)
      # only one closed contour
      lines=[Npt]

   print('found %i points'%Npt)

   # Open the output file
   geo = open(outputfile, 'w')
   geo.write('// This a a geo file created using the python script Contour2geo.py // \n')
   geo.write('Mesh.Algorithm=5; \n')
   geo.write('// To controle the element size, one can directly modify the lc value in the geo file // \n')
   geo.write('lc = {0} ; \n'.format(el_size))
   geo.write('// Mesh size near the boundary from prescribed value  //\n')
   geo.write('Mesh.CharacteristicLengthFromCurvature = 0; \n')
   geo.write('Mesh.CharacteristicLengthFromPoints = 1; \n')
   geo.write('// Give a backgroung field with uniform value for the mesh size  // \n')
   geo.write('Mesh.CharacteristicLengthExtendFromBoundary = 0; \n')
   geo.write('Field[1] = MathEval; \n')
   geo.write('Field[1].F = Sprintf("%g",lc); \n')
   geo.write('Background Field = 1; \n')

   # write the points coordinates (x,y,0,lc)
   np=0
   for j in range(0,Npt):
     np=np+1
     geo.write('Point({0}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(x[j],y[j])+r'}'+'; \n')

   # if spline
   if spline:

      l0=1
      nl=1
      for j in range(0,(len(lines)-1)):
        geo.write('Spline({0}) = '.format(nl)+r'{')
        lastp=lines[j]
        for k in range(0,lastp):
          geo.write('{0},'.format(l0+k))
        geo.write('{0}'.format(l0+lastp) +r'}' + '; \n')
        l0=l0+lastp
        nl=nl+1

      geo.write('Spline({0}) = '.format(nl)+r'{')
      lastp=lines[-1]
      for k in range(0,lastp):
          geo.write('{0},'.format(l0+k))   
      geo.write('1}; \n')

      if len(lines) == 1 :
          geo.write('Line Loop(1) = {1}; \n')
      else:
          geo.write('Line Loop(1) = '+r'{'+'{0}:{1}'.format(1,len(lines)) +r'}' +'; \n')

      geo.write('Plane Surface(1) = {1}; \n')
      geo.write('Physical Surface(1) = {1}; \n')

      for j in range(0,len(lines)):
        geo.write('Physical Line({0}) = '.format(j+1) +r'{' +'{0}'.format(j+1) +r'}'+'; \n')

# else it is lines, as a spline might not work in all cases
   else:

      nl=0
      for j in range(0,Npt-1):
        nl=nl+1
        geo.write('Line({0}) = '.format(nl)+r'{'+'{0},{1}'.format(j+1,j+2)+r'}'+'; \n')
      geo.write('Line({0}) = '.format(nl+1)+r'{'+'{0},{1}'.format(j+2,1)+r'}'+'; \n')

      geo.write('Curve Loop(1) = '+r'{'+'{0}:{1}'.format(1,nl+1)+r'};'+' \n')
      geo.write('Plane Surface(1) = {1}; \n')

      # create physical curves
      l0=1
      phy=0
      for j in range(0,(len(lines)-1)):
        phy=phy+1
        lf=l0+lines[j]-1
        geo.write('Physical Curve({0}) = '.format(phy)+r'{'+'{0}:{1}'.format(l0,lf)+r'};'+' \n')
        l0=lf+1
     
      phy=phy+1
      lf=l0+lines[-1]-1              
      geo.write('Physical Curve({0}) = '.format(phy)+r'{'+'{0}:{1}'.format(l0,lf)+r'};'+' \n')
    
      geo.write('Physical Surface(1) = {1}; \n')

   geo.close()


def usage():
   print('usage: Contour2geo.py -r res [-h] [-i <inputfile>] [-o <outputfile>] [--spline] ')
   print('options:')
   print('   -h [print help]')
   print('   -r res [resolution]')
   print('   -i <inputfile> [default:Contour.txt]')
   print('   -o <outputfile> [default:Contour.geo]')
   print('   --spline [using splines instead of compound lines]')


if __name__ == "__main__":
   main(sys.argv[1:])
