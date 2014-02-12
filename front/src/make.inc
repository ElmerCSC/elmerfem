# ELMER Front makefile, objects and dependencies
# Needs gnu make!
# Martti Verho
# 17.03.03

all: $(TARGET) install

# Result directory
BINDIR = $(ELMER_PREFIX)$(PATH_SEP)bin

#-Search paths
vpath %.$(CPP_EXT) src
vpath %.$(C_EXT) src
vpath %.h include
vpath %.hpp include
vpath %.h ../FrontLib/include

OBJECTS = \
	ecif_body.$(OBJ_EXT) \
	ecif_body2D.$(OBJ_EXT) \
	ecif_body3D.$(OBJ_EXT) \
	ecif_bodyLayer.$(OBJ_EXT) \
	ecif_bodyElement.$(OBJ_EXT) \
	ecif_bodyElement1D.$(OBJ_EXT) \
	ecif_bodyElement2D.$(OBJ_EXT) \
	ecif_bodyElement3D.$(OBJ_EXT) \
	ecif_bodyElementGroup.$(OBJ_EXT) \
	ecif_bodyElementLoop.$(OBJ_EXT) \
	ecif_bodyForce.$(OBJ_EXT) \
	ecif_bodyParameter.$(OBJ_EXT) \
	ecif_boundaryCondition.$(OBJ_EXT) \
	ecif_boundaryParameter.$(OBJ_EXT) \
	ecif_boundbox.$(OBJ_EXT) \
	ecif_calculator.$(OBJ_EXT) \
	ecif_const.$(OBJ_EXT) \
	ecif_constant.$(OBJ_EXT) \
	ecif_control.$(OBJ_EXT) \
	ecif_coordinate.$(OBJ_EXT) \
	ecif_datafile.$(OBJ_EXT) \
	ecif_def_trx.$(OBJ_EXT) \
	ecif_equation.$(OBJ_EXT) \
	ecif_equationVariables.$(OBJ_EXT) \
	ecif_func.$(OBJ_EXT) \
	ecif_geometry.$(OBJ_EXT) \
	ecif_gridH.$(OBJ_EXT) \
	ecif_gridParameter.$(OBJ_EXT) \
	ecif_initialCondition.$(OBJ_EXT) \
	ecif_input.$(OBJ_EXT) \
	ecif_inputAbaqus.$(OBJ_EXT) \
	ecif_inputFront.$(OBJ_EXT) \
	ecif_inputEmf.$(OBJ_EXT) \
	ecif_inputEgf.$(OBJ_EXT) \
	ecif_inputElmer.$(OBJ_EXT) \
	ecif_inputFidap.$(OBJ_EXT) \
	ecif_inputIdeas.$(OBJ_EXT) \
	ecif_inputIdeasWF.$(OBJ_EXT) \
	ecif_inputIges.$(OBJ_EXT) \
	ecif_inputThetis.$(OBJ_EXT) \
	ecif_main.$(OBJ_EXT) \
	ecif_material.$(OBJ_EXT) \
	ecif_mesh.$(OBJ_EXT) \
	ecif_model.$(OBJ_EXT) \
	ecif_modelMeshManager.$(OBJ_EXT) \
	ecif_modelOutputManager.$(OBJ_EXT) \
	ecif_modelObject.$(OBJ_EXT) \
	ecif_modelParameter.$(OBJ_EXT) \
	ecif_model_aux.$(OBJ_EXT) \
	ecif_nurbs.$(OBJ_EXT) \
	ecif_parameter.$(OBJ_EXT) \
	ecif_parameterField.$(OBJ_EXT) \
	ecif_process.$(OBJ_EXT) \
	ecif_renderer.$(OBJ_EXT) \
	ecif_renderer_OGL.$(OBJ_EXT) \
	ecif_simulationParameter.$(OBJ_EXT) \
	ecif_solver.$(OBJ_EXT) \
	ecif_solverControl.$(OBJ_EXT) \
	ecif_timer.$(OBJ_EXT) \
	ecif_timestep.$(OBJ_EXT) \
	ecif_userinterface_TCL.$(OBJ_EXT) \
	ecif_userSettings.$(OBJ_EXT)

# Object dependencies
ecif_body.$(OBJ_EXT): ecif_def.h ecif_body.h ecif_bodyLayer.h ecif_bodyElement.h ecif_model.h
ecif_body2D.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h
ecif_body2D.$(OBJ_EXT): ecif_def.h ecif_body3D.h ecif_bodyElement.h
ecif_bodyLayer.$(OBJ_EXT): ecif_def.h ecif_body.h ecif_bodyLayer.h
ecif_bodyElement.$(OBJ_EXT):ecif_def.h ecif_bodyElement.h ecif_bodyElementGroup.h ecif_geometry.h ecif_model.h
ecif_bodyElement1D.$(OBJ_EXT): ecif_def.h ecif_bodyElement1D.h ecif_geometry.h ecif_renderer.h
ecif_bodyElement2D.$(OBJ_EXT): ecif_def.h ecif_bodyElement2D.h ecif_geometry.h ecif_renderer.h
ecif_bodyElement3D.$(OBJ_EXT): ecif_def.h ecif_bodyElement3D.h ecif_geometry.h ecif_renderer.h 
ecif_bodyElementGroup.$(OBJ_EXT): ecif_def.h ecif_bodyElement.h ecif_bodyElementGroup.h
ecif_bodyElementLoop.$(OBJ_EXT): ecif_def.h ecif_bodyElementLoop.h
ecif_bodyForce.$(OBJ_EXT): ecif_bodyForce.h ecif_parameter.h
ecif_bodyParameter.$(OBJ_EXT): ecif_bodyParameter.h ecif_parameter.h
ecif_boundaryCondition.$(OBJ_EXT): ecif_boundaryCondition.h ecif_bodyElementGroup.h ecif_parameter.h 
ecif_boundaryParameter.$(OBJ_EXT): ecif_boundaryParameter.h ecif_parameter.h
ecif_boundbox.$(OBJ_EXT): ecif_boundbox.h ecif_def_stl.h
ecif_calculator.$(OBJ_EXT): ecif_calculator.h ecif_parameter.h
ecif_const.$(OBJ_EXT): ecif_const.h
ecif_control.$(OBJ_EXT): ecif_def.h ecif_body.h ecif_bodyElement.h ecif_control.h \
ecif_input.h ecif_inputAbaqus.h ecif_inputEmf.h ecif_inputEgf.h ecif_inputElmer.h \
ecif_inputFidap.h ecif_inputIdeas.h ecif_inputIges.h ecif_inputThetis.h \
ecif_material.h ecif_model.h ecif_renderer_OGL.h ecif_userinterface.h ecif_timer.h 
ecif_coordinate.$(OBJ_EXT): ecif_coordinate.h ecif_parameter.h
ecif_datafile.$(OBJ_EXT): ecif_datafile.h ecif_parameter.h
ecif_def_trx.$(OBJ_EXT): ecif_def.h ecif_def_trx.h
ecif_equation.$(OBJ_EXT):ecif_equation.h ecif_parameter.h
ecif_equationVariables.$(OBJ_EXT):ecif_equationVariables.h ecif_parameter.h
ecif_func.$(OBJ_EXT): ecif_func.h ecif_def.h 
ecif_geometry.$(OBJ_EXT):ecif_def.h ecif_boundbox.h ecif_geometry.h 
ecif_gridH.$(OBJ_EXT): ecif_gridH.h
ecif_gridParameter.$(OBJ_EXT): ecif_gridParameter.h ecif_parameter.h
ecif_initialCondition.$(OBJ_EXT): ecif_initialCondition.h ecif_parameter.h
ecif_input.$(OBJ_EXT):ecif_def.h ecif_input.h ecif_model.h 
ecif_inputAbaqus.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputAbaqus.h ecif_model.h
ecif_inputFront.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputFront.h ecif_model.h
ecif_inputEmf.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputEmf.h ecif_model.h
ecif_inputEgf.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputEgf.h ecif_model.h
ecif_inputElmer.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputElmer.h ecif_model.h
ecif_inputFidap.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputFidap.h ecif_model.h
ecif_inputIdeas.$(OBJ_EXT): ecif_def.h ecif_body.h ecif_bodyElement.h ecif_geometry.h \
ecif_inputIdeas.h ecif_inputIdeasWF.h ecif_model.h
ecif_inputIdeasInputWF.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputIdeasWF.h ecif_model.h
ecif_inputIges.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputIges.h ecif_model.h
ecif_inputThetis.$(OBJ_EXT): ecif_def.h ecif_body2D.h ecif_bodyElement.h ecif_bodyElement2D.h \
ecif_geometry.h ecif_inputThetis.h ecif_model.h
ecif_main.$(OBJ_EXT): ecif_control.h ecif_renderer.h ecif_userinterface.h 
ecif_material.$(OBJ_EXT): ecif_material.h ecif_parameter.h
ecif_mesh.$(OBJ_EXT): ecif_def.h ecif_mesh.h
ecif_model.$(OBJ_EXT): ecif_def.h ecif_model.h
ecif_modelObject.$(OBJ_EXT): ecif_modelObject.h
ecif_modelParameter.$(OBJ_EXT): ecif_modelParameter.h ecif_parameter.h
ecif_model_aux.$(OBJ_EXT): ecif_def.h ecif_model_aux.h
ecif_nurbs.$(OBJ_EXT): ecif_nurbs.h
ecif_parameter.$(OBJ_EXT): ecif_parameter.h
ecif_parameterField.$(OBJ_EXT): ecif_parameterField.h
ecif_renderer.$(OBJ_EXT): ecif_def.h ecif_renderer.h 
ecif_renderer_OGL.$(OBJ_EXT): ecif_def.h ecif_body.h ecif_bodyElement.h ecif_control.h ecif_model.h \
ecif_renderer_OGL.h ecif_renderer_OGL_UNIX.hpp ecif_renderer_OGL_UNIX.hpp 
ecif_simulationParameter.$(OBJ_EXT): ecif_simulationParameter.h ecif_parameter.h
ecif_solver.$(OBJ_EXT): ecif_solver.h ecif_parameter.h
ecif_solverControl.$(OBJ_EXT): ecif_solverControl.h ecif_parameter.h
ecif_stdafx.$(OBJ_EXT): ecif_stdafx.h
ecif_timer.$(OBJ_EXT): ecif_timer.h
ecif_timestep.$(OBJ_EXT): ecif_parameter.h ecif_timestep.h 
ecif_userinterface_TCL.$(OBJ_EXT): ecif_def.h ecif_body.h ecif_bodyElement.h ecif_bodyForce.h \
ecif_boundaryCondition.h ecif_constant.h ecif_control.h ecif_coordinate.h ecif_datafile.h ecif_def.h \
ecif_equation.h ecif_geometry.h ecif_gridParameter.h ecif_initialCondition.h ecif_material.h \
ecif_model.h ecif_solver.h ecif_timer.h ecif_timestep.h ecif_userinterface_TCL.h
ecif_userSettings.$(OBJ_EXT): ecif_parameter.h ecif_userSettings.h 

#-Object rules
%.$(OBJ_EXT) : %.$(CPP_EXT) 
	$(CC) $(CFLAGS) $(INCLUDES) $<

%.$(OBJ_EXT) : %.$(C_EXT) 
	$(C) $(CFLAGS) $(INCLUDES) $<


#-Create executable

build: $(TARGET)

$(TARGET):$(OBJECTS)
	$(LINK) $(LIBRARIES) $(OBJECTS) $(LIBS) $(LINKOUT)

# NOTE: install in Windows (xcopy) with gmake needs
# very special path separator handling!
install:
	$(CP) $(TARGET) $(BINDIR)
	$(CP) $(SCRIPT_SOURCEDIR)$(PATH_SEP)* $(BINDIR)
	$(INSTALL) tcl $(TCL_INSTDIR)
	$(INSTALL) lib $(LIB_INSTDIR)

clean:
	$(RM) *.$(OBJ_EXT) $(TARGET)

cleanobject:
	$(RM) *.$(OBJ_EXT)

cleanlink:
	$(RM) $(TARGET)
