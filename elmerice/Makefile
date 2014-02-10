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
Solvers/AIFlowSolve_nlD2.f90 \
Solvers/AIFlowSolve_nlS2.f90 \
Solvers/CaffeSolver.f90 \
Solvers/ComputeDevStressNS.f90 \
Solvers/ComputeEigenValues.f90 \
Solvers/ComputeNormal.f90 \
Solvers/ComputeStrainRate.f90 \
Solvers/DeformationalHeat.f90 \
Solvers/EPLSolver.f90 \
Solvers/fAandfB_in.f90 \
Solvers/FabricSolve.f90 \
Solvers/Flowdepth.f90 \
Solvers/ForceToStress.f90 \
Solvers/GetHydrostaticLoads.f90 \
Solvers/GolfLaw.f90 \
Solvers/GroundedSolver.f90 \
Solvers/IntegratedVelocity.f90 \
Solvers/IDSSolver.f90 \
Solvers/PorousSolve.f90 \
Solvers/pointwise.f90 \
Solvers/SIASolver.f90 \
Solvers/SSASolver.f90 \
Solvers/ThicknessSolver.f90 \
Solvers/TemperateIce.f90 \
Solvers/ExportVertically.f90 \
Solvers/AdjointSolver.f90 \
Solvers/DJDBeta_Adjoint.f90 \
Solvers/DJDmu_Adjoint.f90 \
Solvers/CostSolver_Adjoint.f90 \
Solvers/DJDBeta_Robin.f90 \
Solvers/DJDmu_Robin.f90 \
Solvers/CostSolver_Robin.f90 \
Solvers/m1qn3.f \
Solvers/Grid2DInterpolator.f90 \
Solvers/Optimize_m1qn3Parallel.f90 \
Solvers/OutputStrainHeating.f90 \
Solvers/UpdateExport.f90 \
Solvers/IntegrateVertically.f90
USFSOURCE := \
UserFunctions/Buoyancy.f90 \
UserFunctions/USF_Contact.f90 \
UserFunctions/CaffeFlow.f90 \
UserFunctions/USF_LateralFriction.f90 \
UserFunctions/USF_ShapeFactor.f90 \
UserFunctions/USF_Sliding.f90 \
UserFunctions/USF_Zs.f90\
UserFunctions/USF_GetFrictionHeating.f90 \
UserFunctions/USF_WaterTransfer.f90

solverobjects := $(filter %.o,$(patsubst %.f90,%.o,$(SOLVERSOURCE))) $(filter %.o,$(patsubst %.f,%.o,$(SOLVERSOURCE)))
usfobjects := $(patsubst %.f90,%.o,$(USFSOURCE))

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

%.o: %.f90 
	elmerf90-nosh -c $< -o $@

%.o: %.f
	elmerf90-nosh -c $< -o $@
