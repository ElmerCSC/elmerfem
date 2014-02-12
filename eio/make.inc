# Elmer EIO library makefile, objects and dependencies
# Martti Verho
# 20.11.98

all: cleanobject $(TARGET) install
lib: all	

# Result directories
INCLDIR = $(ELMER_PREFIX)$(PATH_SEP)include
LIBDIR  = $(ELMER_PREFIX)$(PATH_SEP)lib

#-Search paths
vpath %.$(CPP_EXT) src
vpath %.$(C_EXT) src
vpath %.h include

#-Objects
OBJECTS = \
	EIOConstraintAgent.$(OBJ_EXT) \
	EIODualMeshAgent.$(OBJ_EXT) \
	EIOGeometryAgent.$(OBJ_EXT) \
	EIOMeshAgent.$(OBJ_EXT) \
	EIOModelDataAgent.$(OBJ_EXT) \
	EIOModelManager.$(OBJ_EXT) \
	EIOPartWriter.$(OBJ_EXT) \
	EIOSolverAgent.$(OBJ_EXT) \
	eio_api.$(OBJ_EXT)

###	EIOPartReader.$(OBJ_EXT)

#-Objects dependencies (missing currently actually!)
EIOConstraintAgent.$(OBJ_EXT):
EIODualMeshAgent.$(OBJ_EXT):
EIOGeometryAgent.$(OBJ_EXT):
EIOMeshAgent.$(OBJ_EXT):
EIOModelDataAgent.$(OBJ_EXT):
EIOModelManager.$(OBJ_EXT):
### EIOPartReader.$(OBJ_EXT): 
EIOPartWriter.$(OBJ_EXT):
EIOSolverAgent.$(OBJ_EXT):
eio_api.$(OBJ_EXT):


#-Object rules
%.$(OBJ_EXT) : %.$(CPP_EXT) 
	$(CC) $(CFLAGS) $(INCLUDES) $<	
	
%.$(OBJ_EXT) : %.$(C_EXT) 
	(C) $(CFLAGS) $(INCLUDES) $<

#-Library (C or F binding)
$(TARGET): $(OBJECTS)
	$(AR) $(LINKOUT) $(OBJECTS)

install:
	$(CP) $(TARGET) $(LIBDIR)
	$(CP) .$(PATH_SEP)include$(PATH_SEP)eio_api.h  $(INCLDIR)

build: cleanobject $(TARGET)
	
clean:
	$(RM) *.$(OBJ_EXT) $(TARGET)

cleanobject:
	$(RM) *.$(OBJ_EXT)

cleanlink:
	$(RM) $(TARGET)
