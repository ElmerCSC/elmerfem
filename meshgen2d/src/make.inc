# Elmer Mesh generator makefile, objects, dependencies and rules
# Needs Gnu make!
# Martti Verho
# 25.03.99

# Updated for Mesh99
# Reino Ruusu
# 24.1.2000

# Added BoundaryLayer
# Reino Ruusu
# 24.1.2003

all: $(TARGET) install

# Result directory
BINDIR  = $(ELMER_PREFIX)$(PATH_SEP)bin

#-Search paths
vpath %.$(CPP_EXT) src
vpath %.$(C_EXT) src
vpath %.h include

OBJECTS=BGGridMesh.$(OBJ_EXT) \
BGTriangleMesh.$(OBJ_EXT) \
BGVertex.$(OBJ_EXT) \
Body.$(OBJ_EXT) \
BoundaryElement.$(OBJ_EXT) \
BoundaryLayer.$(OBJ_EXT) \
Connect.$(OBJ_EXT) \
Element.$(OBJ_EXT) \
GeometryEdge.$(OBJ_EXT) \
main.$(OBJ_EXT) \
Mesh.$(OBJ_EXT) \
MeshNode.$(OBJ_EXT) \
MeshParser.$(OBJ_EXT) \
Node.$(OBJ_EXT) \
PQ.$(OBJ_EXT) \
QuadElement.$(OBJ_EXT) \
QuadLayer.$(OBJ_EXT) \
SSMFVoronoiSegment.$(OBJ_EXT) \
SSSFVoronoiSegment.$(OBJ_EXT) \
TriangleElement.$(OBJ_EXT) \
TriangleNELayer.$(OBJ_EXT) \
Vertex.$(OBJ_EXT) \
VoronoiSegment.$(OBJ_EXT) \
VoronoiVertex.$(OBJ_EXT) \
VSVertex.$(OBJ_EXT)

# Objects dependencies
BGGridMesh.$(OBJ_EXT):
BGMesh.$(OBJ_EXT):
BGVertex.$(OBJ_EXT):
Body.$(OBJ_EXT):
BoundaryElement.$(OBJ_EXT):
Connect.$(OBJ_EXT):
Element.$(OBJ_EXT):
GeometryEdge.$(OBJ_EXT):
main.$(OBJ_EXT):
Mesh.$(OBJ_EXT):
MeshNode.$(OBJ_EXT):
MeshParser.$(OBJ_EXT):
Node.$(OBJ_EXT):
PQ.$(OBJ_EXT):
QuadElement.$(OBJ_EXT):
QuadLayer.$(OBJ_EXT):
SSMFVoronoiSegment.$(OBJ_EXT):
SSSFVoronoiSegment.$(OBJ_EXT):
TriangleElement.$(OBJ_EXT):
TriangleNELayer.$(OBJ_EXT):
Vertex.$(OBJ_EXT):
VoronoiSegment.$(OBJ_EXT):
VoronoiVertex.$(OBJ_EXT):
VSVertex.$(OBJ_EXT):

#-Object rules
%.$(OBJ_EXT) : %.$(CPP_EXT) 
	$(CC) $(CFLAGS) $(INCLUDES) $<	
	
%.$(OBJ_EXT) : %.$(C_EXT) 
	$(C) $(CFLAGS) $(INCLUDES) $<

#-Exe build
$(TARGET):  $(OBJECTS)
	$(LINK) $(LIBRARIES) $(OBJECTS) $(LIBS) $(LINKOUT)
	
install:
	$(CP) $(TARGET) $(BINDIR)

build: $(TARGET)

clean:
	$(RM) *.$(OBJ_EXT) $(TARGET) 

cleanobject:
	$(RM) *.$(OBJ_EXT)

cleanlink:
	$(RM) $(TARGET)
