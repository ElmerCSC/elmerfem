#----------------------------------------------------------------------
#                       qmake project file for libng
#----------------------------------------------------------------------
include(../ElmerGUI.pri)

TARGET = ng
TEMPLATE = lib
CONFIG -= qt debug
CONFIG += staticlib release warn_off
DESTDIR = ngcore
OBJECTS_DIR = obj
DEFINES += NO_PARALLEL_THREADS
INCLUDEPATH = libsrc/include

unix: QMAKE_CXXFLAGS += -ffriend-injection

SOURCES = libsrc/opti/linopt.cpp \
          libsrc/opti/bfgs.cpp \
          libsrc/opti/linsearch.cpp \
          libsrc/meshing/global.cpp \
          libsrc/meshing/bisect.cpp \
          libsrc/meshing/meshtool.cpp \
          libsrc/meshing/refine.cpp \
          libsrc/meshing/ruler3.cpp \
          libsrc/meshing/improve3.cpp \
          libsrc/meshing/smoothing3.cpp \
          libsrc/meshing/adfront3.cpp \
          libsrc/meshing/tetrarls.cpp \
          libsrc/meshing/prism2rls.cpp \
          libsrc/meshing/pyramidrls.cpp \
          libsrc/meshing/pyramid2rls.cpp \
          libsrc/meshing/netrule3.cpp \
          libsrc/meshing/ruler2.cpp \
          libsrc/meshing/meshclass.cpp \
          libsrc/meshing/improve2.cpp \
          libsrc/meshing/smoothing2.cpp \
          libsrc/meshing/smoothing2.5.cpp \
          libsrc/meshing/adfront2.cpp \
          libsrc/meshing/netrule2.cpp \
          libsrc/meshing/triarls.cpp \
          libsrc/meshing/geomsearch.cpp \
          libsrc/meshing/secondorder.cpp \
          libsrc/meshing/meshtype.cpp \
          libsrc/meshing/parser3.cpp \
          libsrc/meshing/meshing2.cpp \
          libsrc/meshing/quadrls.cpp \
          libsrc/meshing/specials.cpp \
          libsrc/meshing/parser2.cpp \
          libsrc/meshing/meshing3.cpp \
          libsrc/meshing/meshfunc.cpp \
          libsrc/meshing/localh.cpp \
          libsrc/meshing/improve2gen.cpp \
          libsrc/meshing/delaunay.cpp \
          libsrc/meshing/boundarylayer.cpp \
          libsrc/meshing/msghandler.cpp \
          libsrc/meshing/meshfunc2d.cpp \
          libsrc/meshing/topology.cpp \
          libsrc/meshing/clusters.cpp \
          libsrc/meshing/curvedelems_new.cpp \
          libsrc/meshing/hprefinement.cpp \
          libsrc/meshing/validate.cpp \
          libsrc/interface/nglib.cpp \
          libsrc/gprim/geomtest3d.cpp \
          libsrc/gprim/geom2d.cpp \
          libsrc/gprim/geom3d.cpp \
          libsrc/gprim/adtree.cpp \
          libsrc/gprim/transform3d.cpp \
          libsrc/gprim/geomfuncs.cpp \
          libsrc/linalg/polynomial.cpp \
          libsrc/linalg/densemat.cpp \
          libsrc/linalg/vector.cpp \
          libsrc/csg/algprim.cpp \
          libsrc/csg/brick.cpp \
          libsrc/csg/manifold.cpp \
          libsrc/csg/bspline2d.cpp \
          libsrc/csg/meshsurf.cpp \
          libsrc/csg/csgeom.cpp \
          libsrc/csg/polyhedra.cpp \
          libsrc/csg/curve2d.cpp \
          libsrc/csg/singularref.cpp \
          libsrc/csg/edgeflw.cpp \
          libsrc/csg/solid.cpp \
          libsrc/csg/explicitcurve2d.cpp \
          libsrc/csg/specpoin.cpp \
          libsrc/csg/gencyl.cpp \
          libsrc/csg/revolution.cpp \
          libsrc/csg/genmesh.cpp \
          libsrc/csg/spline3d.cpp \
          libsrc/csg/surface.cpp \
          libsrc/csg/identify.cpp \
          libsrc/csg/triapprox.cpp \
          libsrc/csg/csgparser.cpp \
          libsrc/csg/extrusion.cpp \
          libsrc/geom2d/geom2dmesh.cpp \
          libsrc/geom2d/spline.cpp \
          libsrc/geom2d/splinegeometry.cpp \
          libsrc/geom2d/genmesh2d.cpp \
          libsrc/stlgeom/meshstlsurface.cpp \
          libsrc/stlgeom/stlline.cpp \
          libsrc/stlgeom/stltopology.cpp \
          libsrc/stlgeom/stltool.cpp \
          libsrc/stlgeom/stlgeom.cpp \
          libsrc/stlgeom/stlgeomchart.cpp \
          libsrc/stlgeom/stlgeommesh.cpp \
          libsrc/general/moveablemem.cpp \
          libsrc/general/ngexception.cpp \
          libsrc/general/table.cpp \
          libsrc/general/optmem.cpp \
          libsrc/general/spbita2d.cpp \
          libsrc/general/hashtabl.cpp \
          libsrc/general/sort.cpp \
          libsrc/general/flags.cpp \
          libsrc/general/seti.cpp \
          libsrc/general/bitarray.cpp \
          libsrc/general/symbolta.cpp \
          libsrc/general/mystring.cpp \
          libsrc/general/profiler.cpp
#          libsrc/general/array.cpp
