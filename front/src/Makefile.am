INCLUDES = \
        -I$(top_builddir)/include       \
        -I$(top_srcdir)/include         

bin_PROGRAMS=ElmerFront

COMMON_SRCS = \
	ecif_body2D.cpp \
	ecif_body3D.cpp \
	ecif_body.cpp \
	ecif_bodyElement1D.cpp \
	ecif_bodyElement2D.cpp \
	ecif_bodyElement3D.cpp \
	ecif_bodyElement.cpp \
	ecif_bodyElementGroup.cpp \
	ecif_bodyElementLoop.cpp \
	ecif_bodyForce.cpp \
	ecif_bodyLayer.cpp \
	ecif_bodyParameter.cpp \
	ecif_boundaryCondition.cpp \
	ecif_boundaryParameter.cpp \
	ecif_boundbox.cpp \
	ecif_calculator.cpp \
	ecif_constant.cpp \
	ecif_const.cpp \
	ecif_control.cpp \
	ecif_coordinate.cpp \
	ecif_datafile.cpp \
	ecif_def_trx.cpp \
	ecif_equation.cpp \
	ecif_equationVariables.cpp \
	ecif_func.cpp \
	ecif_geometry.cpp \
	ecif_gridH.cpp \
	ecif_gridParameter.cpp \
	ecif_initialCondition.cpp \
	ecif_inputAbaqus.cpp \
	ecif_input.cpp \
	ecif_inputEgf.cpp \
	ecif_inputElmer.cpp \
	ecif_inputEmf.cpp \
	ecif_inputFidap.cpp \
	ecif_inputFront.cpp \
	ecif_inputIdeas.cpp \
	ecif_inputIdeasWF.cpp \
	ecif_inputIges.cpp \
	ecif_inputThetis.cpp \
	ecif_main.cpp \
	ecif_material.cpp \
	ecif_mesh.cpp \
	ecif_model_aux.cpp \
	ecif_model.cpp \
	ecif_modelMeshManager.cpp \
	ecif_modelObject.cpp \
	ecif_modelOutputManager.cpp \
	ecif_modelParameter.cpp \
	ecif_nurbs.cpp \
	ecif_parameter.cpp \
	ecif_parameterField.cpp \
	ecif_process.cpp \
	ecif_renderer.cpp \
	ecif_renderer_OGL.cpp \
	ecif_simulationParameter.cpp \
	ecif_solverControl.cpp \
	ecif_solver.cpp \
	ecif_timer.cpp \
	ecif_timestep.cpp \
	ecif_userinterface_TCL.cpp \
	frontlib.cpp \
	ecif_userSettings.cpp

if HAVE_WINDRES

.rc.o:
	windres $^ -o $@
%.o : %.rc
	windres $^ -o $@

ElmerFront_SOURCES = $(COMMON_SRCS) ElmerFront-res.rc

else

ElmerFront_SOURCES =  $(COMMON_SRCS)

endif

install-data-local:
	@echo "Installing tcl scripts to ${prefix}/share/elmerfront/tcl"

	mkdir -p "$(DESTDIR)${prefix}/share/elmerfront/tcl"
	cp -r ./tcl/* "$(DESTDIR)${prefix}/share/elmerfront/tcl"
	rm -rf `find $(DESTDIR)$(prefix)/share/elmerfront/tcl -name .svn`

	@echo "Installing *.edf *.rgb to ${prefix}/share/elmerfront/lib"
	mkdir -p "$(DESTDIR)${prefix}/share/elmerfront/lib"
	cp -r ./lib/* "$(DESTDIR)${prefix}/share/elmerfront/lib"
	rm -rf `find $(DESTDIR)$(prefix)/share/elmerfront/lib -name .svn`



