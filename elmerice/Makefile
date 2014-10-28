##############################################
# Makefile for Elmer/Ice
# (c) 2012 CSC - IT Center for Science Ltd.
# Author: T. Zwinger
# Original Date: 2012-12-03
#############################################

PLATFORM := $(shell uname -s)
ELMERICE_HOME ?= $(ELMER_HOME)/share/elmersolver
EXECUTABLES := ElmerIceSolvers.$(PLATFORM) ElmerIceUSF.$(PLATFORM)
INSTALLED := $(ELMERICE_HOME)/include/defgrid.mod
MODULES := defgrid.mod
MKDIR_P := mkdir -p
SOLVERSOURCE := \
Solvers/AIFlowSolve_nlD2.F90 \
Solvers/AIFlowSolve_nlS2.F90 \
Solvers/CaffeSolver.F90 \
Solvers/ComputeDevStressNS.F90 \
Solvers/ComputeEigenValues.F90 \
Solvers/ComputeNormal.F90 \
Solvers/ComputeStrainRate.F90 \
Solvers/DeformationalHeat.F90 \
Solvers/EPLSolver.F90 \
Solvers/fAandfB_in.F90 \
Solvers/FabricSolve.F90 \
Solvers/Flowdepth.F90 \
Solvers/ForceToStress.F90 \
Solvers/GetHydrostaticLoads.F90 \
Solvers/GolfLaw.F90 \
Solvers/GroundedSolver.F90 \
Solvers/IntegratedVelocity.F90 \
Solvers/IDSSolver.F90 \
Solvers/PorousSolve.F90 \
Solvers/pointwise.F90 \
Solvers/SIASolver.F90 \
Solvers/SSASolver.F90 \
Solvers/ThicknessSolver.F90 \
Solvers/TemperateIce.F90 \
Solvers/ExportVertically.F90 \
Solvers/AdjointSolver.F90 \
Solvers/DJDBeta_Adjoint.F90 \
Solvers/DJDmu_Adjoint.F90 \
Solvers/CostSolver_Adjoint.F90 \
Solvers/DJDBeta_Robin.F90 \
Solvers/DJDmu_Robin.F90 \
Solvers/CostSolver_Robin.F90 \
Solvers/m1qn3.F \
Solvers/Grid2DInterpolator.F90 \
Solvers/Optimize_m1qn3Parallel.F90 \
Solvers/OutputStrainHeating.F90 \
Solvers/UpdateExport.F90 \
Solvers/IntegrateVertically.F90\
Solvers/EnthalpySolver.F90
USFSOURCE := \
UserFunctions/Buoyancy.F90 \
UserFunctions/USF_Contact.F90 \
UserFunctions/CaffeFlow.F90 \
UserFunctions/USF_LateralFriction.F90 \
UserFunctions/USF_ShapeFactor.F90 \
UserFunctions/USF_Sliding.F90 \
UserFunctions/USF_Zs.F90\
UserFunctions/USF_GetFrictionHeating.F90 \
UserFunctions/USF_WaterTransfer.F90

solverobjects := $(filter %.o,$(patsubst %.F90,%.o,$(SOLVERSOURCE))) $(filter %.o,$(patsubst %.f,%.o,$(SOLVERSOURCE)))
usfobjects := $(patsubst %.F90,%.o,$(USFSOURCE))

all: compile


compile:  $(EXECUTABLES)

install: directories
	@if test '$(PLATFORM)' = 'Darwin'; then \
           echo "cp ElmerIceSolvers.dylib $(ELMERICE_HOME)/lib"; \
           echo "cp ElmerIceUSF.dylib $(ELMERICE_HOME)/lib"; \
	   cp ElmerIceSolvers.dylib $(ELMERICE_HOME)/lib; \
           cp ElmerIceUSF.dylib $(ELMERICE_HOME)/lib; \
        else \
           if test '$(PLATFORM)' = 'Linux'; then \
                echo "cp ElmerIceSolvers.so $(ELMERICE_HOME)/lib"; \
                echo "cp ElmerIceUSF.so $(ELMERICE_HOME)/lib"; \
	   	cp ElmerIceSolvers.so $(ELMERICE_HOME)/lib; \
                cp ElmerIceUSF.so $(ELMERICE_HOME)/lib; \
	   fi \
        fi

	cp $(MODULES) $(ELMERICE_HOME)/include

directories: $(ELMERICE_HOME)/lib $(ELMERICE_HOME)/include

$(ELMERICE_HOME)/lib:
	$(MKDIR_P) $(ELMERICE_HOME)/lib
$(ELMERICE_HOME)/include:
	$(MKDIR_P) $(ELMERICE_HOME)/include

ElmerIceSolvers.$(PLATFORM): $(solverobjects)
	@if test '$(PLATFORM)' = 'Darwin'; then \
	   echo "Dynamic linking solvers to $(PLATFORM)"; \
           elmerf90 $(solverobjects) -o ElmerIceSolvers.dylib; \
        else \
           if test '$(PLATFORM)' = 'Linux'; then \
	   	echo "Dynamic linking solvers to $(PLATFORM)"; \
                elmerf90  $(solverobjects) -o ElmerIceSolvers.so; \
	   fi \
        fi

ElmerIceUSF.$(PLATFORM): $(usfobjects)
	@if test '$(PLATFORM)' = 'Darwin'; then \
	   echo "Dynamic linking user functions to $(PLATFORM)"; \
           elmerf90 $(usfobjects) -o ElmerIceUSF.dylib; \
        else \
           if test '$(PLATFORM)' = 'Linux'; then \
	   	echo "Dynamic linking user functions to $(PLATFORM)"; \
                elmerf90 $(usfobjects) -o ElmerIceUSF.so; \
	   fi \
	fi

clean:
	@if test '$(PLATFORM)' = 'Darwin'; then \
	   echo "rm -f ElmerIceSolvers.dylib ElmerIceUSF.dylib"; \
	   rm -f ElmerIceSolvers.dylib ElmerIceUSF.dylib; \
        else \
           if test '$(PLATFORM)' = 'Linux'; then \
               echo "rm -f ElmerIceSolvers.so ElmerIceUSF.so"; \
	       rm -f ElmerIceSolvers.so ElmerIceUSF.so; \
	   fi \
	fi

	rm -f  Solvers/*.o UserFunctions/*.o $(MODULES)

purge: clean
	@if test '$(PLATFORM)' = 'Darwin'; then \
           echo "rm -f  $(ELMERICE_HOME)/lib/ElmerIceSolvers.dylib $(ELMERICE_HOME)/lib/ElmerIceUSF.dylib"; \
	   rm -f  $(ELMERICE_HOME)/lib/ElmerIceSolvers.dylib $(ELMERICE_HOME)/lib/ElmerIceUSF.dylib; \
        else \
           if test '$(PLATFORM)' = 'Linux'; then \
	       echo "rm -f  $(ELMERICE_HOME)/lib/ElmerIceSolvers.so $(ELMERICE_HOME)/lib/ElmerIceUSF.so"; \
	       rm -f  $(ELMERICE_HOME)/lib/ElmerIceSolvers.so $(ELMERICE_HOME)/lib/ElmerIceUSF.so; \
	   fi \
	fi

	rm -f $(INSTALLED)

%.o: %.F90 
	elmerf90-nosh -c $< -o $@

%.o: %.F
	elmerf90-nosh -c $< -o $@
