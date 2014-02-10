#/*****************************************************************************/
# *
# *  Elmer, A Finite Element Software for Multiphysical Problems
# *
# *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
# * 
# *  This program is free software; you can redistribute it and/or
# *  modify it under the terms of the GNU General Public License
# *  as published by the Free Software Foundation; either version 2
# *  of the License, or (at your option) any later version.
# * 
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program (in file fem/GPL-2); if not, write to the 
# *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
# *  Boston, MA 02110-1301, USA.
# *
# *****************************************************************************/
   
#***********************************************************************
#Program:   ELMER Front 
#Module:    ecif_tk_globParams.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Global variables
#
#************************************************************************

	global UserSetting

  # Solver keyword (SOLVER.KEYWORDS-file) sections
	# and short names for them
	# ========================
	# NOTE: Keyword Type info is stored like:
	# SKWD(BF,heatsource,0) if not given
	# SKWD(BF,heatsource,1) if given
	#
	global SKWD 
	
	# Short names for sections
	#			                            
  set SKWD(BodyForce)          BF
  set SKWD(Body)               BO
  set SKWD(BoundaryCondition)  BC
  set SKWD(Boundary)           BA
  set SKWD(Constant)           CN
  set SKWD(Equation)           EQ
  set SKWD(InitialCondition)   IC
  set SKWD(Material)           MT
  set SKWD(Simulation)         SI
  set SKWD(Solver)             SL


#==========================================#
#                                          #
#             FRONT_INIT                   #
#                                          #
#==========================================#
#
# Initialize global arries
#
#
proc FRONT_INIT {} {

  # All Front global arries
	# =======================
  global Model ObjectTable MatcDefinitions
  global Common Info Status ModelFlags 
  global BodyProperty BodyDisplay LabelDisplay
  global Boundary BoundaryDisplay VertexDisplay MeshDefine MeshSelect PostFileSelect
  global ProcessIds ProcessTable
  global InputFileInfo ModelInfo BodyInfo
  global RendererInfo SystemInfo SourceList TclCommand
  global \
      BodyParameter \
      BoundaryParameter \
      BodyForce \
      BoundaryCondition \
      Coordinate \
      Calculator \
      Constant \
      Datafile \
      Equation \
      EquationVariable \
      GridParameter \
      GridH \
      InitialCondition \
      Material \
      ModelParameter \
      ModelProperty \
      Processor \
      SimulationParameter \
      Solver \
      SolverControl \
      SolverOrder \
      SolverSystem \
      Timestep \
      UserSetting \
      Variables

  # These are arries which are "panel" related
  #
  set Common(allPanelArries) {
     BodyParameter BoundaryParameter ModelParameter SimulationParameter
     BodyForce BoundaryCondition Boundary Calculator Constant Coordinate Datafile
     Equation EquationVariable MeshDefine GridParameter GridH InitialCondition Material
     ModelProperty Processor Solver SolverControl SolverOrder Timestep UserSetting Variables
  }

  # These panels are handled by "StandardPanel" procs. They need special
  # calls for cancel etc. when called "outside"
  #
  set Common(standardPanelArries) {
     BodyParameter BoundaryParameter
     BodyForce BoundaryCondition Equation
     GridParameter GridH InitialCondition Material
     Calculator Solver
  }

  # These parameter arries must be checked when an equation
  # has been changed
  #
  set Common(maskedParameters) {
    BodyForce BoundaryCondition InitialCondition Material
  }

  set Common(currentArray) ""
  set Common(currentStandardPanel) ""

  #---Parameter panel arries menu names
  #
  set BodyForce(menuName)           "Body forces"
  set BodyParameter(menuName)       "Body parameters"
  set BoundaryCondition(menuName)   "Boundary conditions"
  set Boundary(menuName)            "Boundaries"
  set BoundaryParameter(menuName)   "Boundary parameters"
  set Coordinate(menuName)          "Coordinate settings"
  set Calculator(menuName)          "Calculators"
  set Constant(menuName)            "Physical constants"
  set Datafile(menuName)            "Datafiles"
  set Equation(menuName)            "Equations"
  set EquationVariable(menuName)    "Equation variables"
  set GridParameter(menuName)       "Mesh structure"
  set GridH(menuName)               "Mesh density"
  set InitialCondition(menuName)    "Initial conditions"
  set Material(menuName)            "Materials"
  set MeshDefine(menuName)          "Mesh define"
  set ModelParameter(menuName)      "Model parameters"
  set ModelProperty(menuName)       "Model name and directories"
  set SimulationParameter(menuName) "Simulation parameters"
  set Solver(menuName)              "Solver settings"
  set SolverControl(menuName)       "Solver control"
  set Processor(menuName)           "Processor settings"
  set SolverOrder(menuName)         "Equation solving order"
  set Timestep(menuName)            "Timestep settings"
  set UserSetting(menuName)         "Settings"
  set Variables(menuName)           "Variables"

  #---Parameter panel array codes
  set BodyForce(parameterType)           "BF"
  set BodyParameter(parameterType)       "BodyP"
  set BoundaryCondition(parameterType)   "BC"
  set Boundary(parameterType)            "BA"
  set BoundaryParameter(parameterType)   "BndrP"
  set Coordinate(parameterType)          "CO"
  set Calculator(parameterType)          "CL"
  set Constant(parameterType)            "CN"
  set Datafile(parameterType)            "DF"
  set Equation(parameterType)            "EQ"
  set EquationVariable(parameterType)    "EV"
  set GridParameter(parameterType)       "GR"
  set GridH(parameterType)               "GH"
  set InitialCondition(parameterType)    "IC"
  set Material(parameterType)            "MT"
  set MeshDefine(parameterType)          "MD"
  set ModelParameter(parameterType)      "ModelP"
  set ModelProperty(parameterType)       "MP"
  set SimulationParameter(parameterType) "SimuP"
  set Solver(parameterType)              "SL"
  set SolverControl(parameterType)       "SLC"
  set Processor(parameterType)           "PR"
  set SolverOrder(parameterType)         "SO"
  set Timestep(parameterType)            "TS"
  set UserSetting(parameterType)         "US"
  set Variables(parameterType)           "VA"


  #
  #---Panel array names in Common
  set Common(panelArray,BF)     "BodyForce"
  set Common(panelArray,BodyP)  "BodyParameter"
  set Common(panelArray,BC)     "BoundaryCondition"
  set Common(panelArray,BA)     "Boundary"
  set Common(panelArray,BndrP)  "BoundaryParameter"
  set Common(panelArray,CL)     "Calculator"
  set Common(panelArray,CN)     "Constant"
  set Common(panelArray,CO)     "Coordinate"
  set Common(panelArray,DF)     "Datafile"
  set Common(panelArray,EQ)     "Equation"
  set Common(panelArray,EV)     "EquationVariable"
  set Common(panelArray,GR)     "GridParameter"
  set Common(panelArray,GH)     "GridH"
  set Common(panelArray,IC)     "InitialCondition"
  set Common(panelArray,MD)     "MeshDefine"
  set Common(panelArray,MT)     "Material"
  set Common(panelArray,ModelP) "ModelParameter"
  set Common(panelArray,MP)     "ModelProperty"
  set Common(panelArray,PR)     "Processor"
  set Common(panelArray,SimuP)  "SimulationParameter"
  set Common(panelArray,SL)     "Solver"
  set Common(panelArray,SLC)    "SolverControl"
  set Common(panelArray,SO)     "SolverOrder"
  set Common(panelArray,TS)     "Timestep"
  set Common(panelArray,US)     "UserSetting"
  set Common(panelArray,VA)     "Variables"

  #==========================================================================#
  #==========================================================================#

  #---Status fields
  set Status(useAbsoluteBase) 1  ;# Absolute base like Model(nofBodies); relative like Model(nofBodieWithEquation)
  MenuExec::initStatusFields

  #---Model environment, properties panel stuff, initial values
  #
  set Model(CASE_NAME) ""
  set Model(GEOMETRY_DIMENSION)   ""
  set Model(SIMULATION_DIMENSION) ""
  set Model(GEOMETRY_TYPE) "None"
  set Model(cadPath) ""
  set Model(cadFileType) 0        ;# Default: Elmer egf
  set Model(meshPath) ""
  set Model(meshFileType) 2       ;# Default: Elmer 2D
  set Model(CAD_OPEN_DIR) ""
  set Model(MESH_OPEN_DIR) ""
  set Model(MESH_SAVE_DIR) ""
  set Model(EMF_HOME) ""          ;# Default directory for Ecf model files
  set Model(EMF_PATH) ""          ;# Ecf model file full name (dir + filename)
  set Model(EMF_PATH_IN) ""       ;# Ecf model infile full name (dir + filename)
  set Model(EMF_FILE) ""          ;# Ecf model file name (without directory)
  set Model(EMF_SAVE_DIR) ""      ;# Current save directory for Ecf model files
  set Model(EMF_OPEN_DIR) ""      ;# Current read directory for Ecf model files
  set Model(outElmerPostMeshFile) ""
  set Model(outThetisMeshFile) ""

  set Model(meshInputUnit) 1.0    ;# External mesh scaling factor (1.0 <--> m, 0.001 <--> mm ect)
  set Model(meshNames) ""         ;# All meshdir names under MESHDIR
  set Model(meshHs) ""            ;# All meshH values
  set Model(meshFs) ""            ;# All meshF values
  set Model(nofActiveMeshes) 0    ;# Nof different meshes in active solvers

  set Model(status) 0             ;# Model status ok
  set Model(statusMessage) ""

  set Model(hasUserDefinitions) 0
  set Model(hasMatcDefinitions) 0
  set Model(inDefinitionFile) ""
  set Model(inDefinitionFileDir) ""
  set Model(outDefinitionFile) ""

  set Model(allDefinitions) ""
  set Model(newDefinitions) ""
  set Model(removedFieldsList) ""

  set Model(matcInputFile_emf) ""
  set Model(matcInputFile_sif) ""

  set Model(Solver,askNonExistingMesh) 1
  set Model(Solver,solving) 0

  set ModelFlags(GEOMETRY_TYPE_CAD) 0
  set ModelFlags(GEOMETRY_TYPE_MESH) 0
  set ModelFlags(DRAW_SOURCE_CAD) 0
  set ModelFlags(DRAW_SOURCE_MESH) 0
  set ModelFlags(DRAW_TARGET_BODIES) 0
  set ModelFlags(DRAW_TARGET_SURFACES) 0
  set ModelFlags(SELECT_OBJECTS_EXTEND) 0

  set ModelFlags(LABEL_DISPLAY_NODE) 0
  set ModelFlags(LABEL_DISPLAY_ELEMENT) 0
  set ModelFlags(LABEL_DISPLAY_VERTEX) 0
  set ModelFlags(LABEL_DISPLAY_EDGE) 0
  set ModelFlags(LABEL_DISPLAY_FACE) 0
  set ModelFlags(LABEL_DISPLAY_BODY) 0

  # Init model directories
  #set ModelProperty(MODEL_NAME) ""
  #ModelProperty::initCaseDirectoryVariables

  set MeshDefine(MESH_H,default) -1
  set MeshDefine(MESH_F,default) 1.0

  set PostFileSelect(postFile) ""
  set PostFileSelect(addDirectory) ""
  set PostFileSelect(fileList) ""
  set PostFileSelect(serverId) 0

  set RendererInfo(LINE_WIDTH_GRANULARITY) ""
  set RendererInfo(LINE_WIDTH_RANGE) ""

  #---Names for ELMER components starting scripts

  # If process exec messages are shown is message area
  set Info(displayExecMessages) 0

  # If Matc definition should be drop from model file
  set Info(dropModelMatcDefinitions) 0

  #--Info array variable names for the processes
  #
  # NOTE: The variable name (the first part before comma) and the text
  # MUST match, because the text is used as the array variable name!!!
  #
  # This (a bit stupid) indirective naming is used because we can call
  # eg. solver-process stuff more safely like:
  # set process_name $Info(Solver,processName)
  # set exec_name $Info($process,proc,WIN32)
  # and NOT like:
  # set process_name "Solver"
  # where erronous writing could go unnoticed!
  #
  
  #--Internal names for the processes
  #  Do NOT change these values, they are 'keys'!!!
  #
  set Info(Emf2Db,processName)          "Emf2DB"
  set Info(GebhardtFactors,processName) "GebhardtFactors"
  set Info(Mesh2D,processName)          "Mesh2D"
  set Info(Mesh3D,processName)          "Mesh3D"
  set Info(Results,processName)         "Results"
  set Info(Post,processName)            "Post"
  set Info(Solver,processName)          "Solver"
  set Info(Viewfactors,processName)     "Viewfactors"
  set Info(F90,processName)             "F90"

  #--External executable names for the processes
  set Info(Emf2Db,proc,WIN32)           "Emf2Db.exe"
  set Info(GebhardtFactors,proc,WIN32)  "GebhardtFactors.exe"
  set Info(Mesh2D,proc,WIN32)           "ElmerMesh2D.exe"
  set Info(Mesh3D,proc,WIN32)           "ElmerMesh3D.exe"  ;# Not in use:-)
  set Info(Results,proc,WIN32)          "ElmerPost.exe"
  set Info(Post,proc,WIN32)             "ElmerPost.exe"
  set Info(Solver,proc,WIN32)           "ElmerSolver.exe"
  set Info(Viewfactors,proc,WIN32)      "ViewFactors.exe"
  set Info(F90,proc,WIN32)              "f90.exe"

  set Info(Emf2Db,proc,UNIX)          "Emf2Db"
  set Info(GebhardtFactors,proc,UNIX) "GebhardtFactors"
  set Info(Mesh2D,proc,UNIX)          "ElmerMesh2D"
  set Info(Mesh3D,proc,UNIX)          "Mesh3D"
  set Info(Results,proc,UNIX)         "ElmerPost"
  set Info(Post,proc,UNIX)            "ElmerPost"
  set Info(Solver,proc,UNIX)          "ElmerSolver"
  set Info(Viewfactors,proc,UNIX)     "ViewFactors"
  set Info(F90,proc,UNIX)             "f90"

  #--Short display names for the processes
  set Info(Emf2Db,processTag)          "Emf2DB"
  set Info(GebhardtFactors,processTag) "GebFact"
  set Info(Mesh2D,processTag)          "Mesh2D"
  set Info(Mesh3D,processTag)          "Mesh3D"
  set Info(Results,processTag)         "ElmerPost"
  set Info(Post,processTag)            "ElmerPost"
  set Info(Solver,processTag)          "Solver"
  set Info(Viewfactors,processTag)     "Viewfact"
  set Info(F90,processTag)             "F90Comp"

  set Info(Emf2Db,checkAllDone)          0
  set Info(GebhardtFactors,checkAllDone) 1
  set Info(Mesh2D,checkAllDone)          1
  set Info(Mesh3D,checkAllDone)          1
  set Info(Results,checkAllDone)         0
  set Info(Post,checkAllDone)            0
  set Info(Solver,checkAllDone)          1
  set Info(Viewfactors,checkAllDone)     1
  set Info(F90,checkAllDone)             0

  set Info(Results,tclServer)   "ElmerPost-$Info(PID)"
  set Info(Results,tclCommand)  ""

  set Info(defaultLogDirectoryName) "LOGDIR"
  set Info(meshDirectoryName) "MESHDIR"
  set Info(meshInfoFileName) "ELMERFRONT_MESHINFO"

  #---Names for ELMER logfiles
  set Info(Emf2Db,logfile) ElmerEmf2Db
  set Info(GebhardtFactors,logfile) ElmerGebhardtFactors
  set Info(Mesh2D,logfile) ElmerMesh2D
  set Info(Mesh3D,logfile) ElmerMesh3D
  set Info(Results,logfile) ElmerPost
  set Info(Post,logfile) ElmerPost
  set Info(Solver,logfile) ElmerSolver
  set Info(Viewfactors,logfile) ElmerViewfactors
  set Info(F90,logfile) ElmerF90Compiler

  #---User entered procedure libraries settings in the
  #   procedure entry panel are stored here (Note this is
  #   gloabl value, common for all procedure panels during
  #   the session
  set Info(procedureLibraries) ""
  set Info(lastProcedureDir) ""

  #---Background process priority levels (Mesh, Solver, *factors etc.)
  set Info(ElmerRun,lowPriority) "LOW_PRIORITY"
  set Info(ElmerRun,normalPriority) "NORMAL_PRIORITY"
  set Info(ElmerRun,highPriority) "HIGH_PRIORITY"
  set Info(bgPriorityLevel) $Info(ElmerRun,lowPriority)

  #---Process run sequence nubers to identify log files
  set Info(ElmerRun,runnumber) 0

  set Info(browser,id) 0
  set Info(browser,smallDelay)   50
  set Info(browser,mediumDelay) 200
  set Info(browser,largeDelay)  800
  
  set Info(browser,maxCharCount) [expr 4*1000*1000]  ;# Max text size (4MB) for a process browser window

  set Info(processTableMonitorInterval) 50  ;# In milliseconds!
  set Info(processTableRefreshInterval) 50  ;# In milliseconds!
  	
	#--Paths and timestamps for SOLVER.KEYWORDS files in use
	# Elmer /lib; model file's directory
	set Info(file_SKWD_lib) ""
	set Info(mtime_SKWD_lib) 0
	set Info(file_SKWD_model) ""
	set Info(mtime_SKWD_model) 0

  #---User levels
  set Info(noviceUser)   1
  set Info(advancedUser) 2
  set Info(powerUser) 3

  set Info(userLevel) $Info(advancedUser)
  set Info(userLevel,new) $Info(advancedUser)
  set Info(userLevel,prev) $Info(advancedUser)

  #---Initially process list is empty
  set ProcessIds ""

  set Info(nofEquationedPanels) 4 
  set Info(outerType) 2 
  set Info(innerType) 4
  #
	set Info(defaultMatcPrecision) 8

  set Info(NO_INDEX) -1
  set Info(arguments) ""
  set Info(results) ""
  set Info(cmdSeparator) @
  set Info(argSeparator) ^
  set Info(objectSeparator) !
  set Info(fieldSeparator) &

  set Info(matcMarker) $
  
  #
  set Info(menuDots) "..."
  set Info(dispListSeparatorIn) "|"
  set Info(dispListSeparatorExtra1) ""
  set Info(dispListSeparatorExtra2) ""
  set Info(leftFieldSizeSep) "("
  set Info(rightFieldSizeSep) ")"
  set Info(valueFlag) "=="
  set Info(fileFlag)  "=:"
  set Info(procFlag)  "=."
  set Info(tableFlag) "=="
  set Info(matcFlag) "==$"
  set Info(dispListSeparatorOut) "|"  ;#-NOTE: outsep = extra1+insep+extra2
  set Info(groupDataSeparator) ";" 
  set Info(arrayDataSeparator) ";" 
  set Info(dataListSeparator) ";" 
  set Info(includePathSeparator) ";" 

  #-NOTE: We keep inactive data to avoid loosing possibly
  # complicated parameter definitions when an equation is
  # turned off!!!
  #
  set Info(keepInactiveFields) 1

  #
  # NOTE: Keep these character unique!!!
  set Info(maskStart) "¤"
  set Info(variableMaskStart) "%"
  set Info(variableMaskStop) "¥"

  set Info(maskTrueMark) "½"
  set Info(maskFalseMark) "§"
  set Info(maskSepMark) "å"

  set Info(inactiveMark) "-" 
  set Info(selectedMark) "#" 
  set Info(negationMark) "~" 
  set Info(paramMarked) "*" 
  set Info(paramUnmarked) " " 
  set Info(selectionBoxTrueMarker) "\[X\]"
  set Info(selectionBoxFalseMarker) "\[ \]"
  #
  set Info(defaultApplyState) disabled
  set Info(defaultCancelState) normal  ;# NOTE: Cancel is active all the time!
  set Info(messageIcon) info
  set Info(anywayOk) "Click OK to continue anyway"
  set Info(continueOk) "Click OK to continue"
  set Info(noMoreNotifying) [join [list \
        "Click Yes to continue and ignore the message for the session.\n\n" \
        "Click No to continue and ignore the meassage only this time." ]]
  set Info(setModelNameMsg) "You can set the name using menus: Problem/Name and dir..." 
  set Info(fieldMsg) "Data in the following field(s) is incorrect:\n"
  set Info(numberFieldMsg) " (value for a number)"
  set Info(fileFieldMsg) " (filename)"
  set Info(showDataFormatErrorMsg) 1
  set Info(messageWindow,maxSize) 100

  #
  set Info(markMainWindowTitle) 1
  set Info(windowList) ""
  set Info(thisWindow) {}
  set Info(hasRadiation) 0
  set Info(drawBodies) 0
  set Info(drawBoundaries) 0
  set Info(drawMesh) 0
  set Info(cntrlDblSelected) 0
  set Info(controlPressed) 0

  set Info(NEXT_ACTIVE_SELECTION_TOLERANCE) ""
  set Info(nextActiveSelectionTolerance) ""
  set Info(lastAppliedSelectionMethod) ""
  set Info(lastAppliedSelectionMode) ""

  set Info(currentMoveAxis) ""
  set Info(currentRotateAxis) ""

  set Info(informAboutRemovedFields) 0
  
  # Input file info lists
  set Info(solverKeywordFiles) ""
  set Info(matcFiles) ""
  set Info(colorFiles) ""
  set Info(definitionFiles) ""
  
  set Info(browserDir) $Info(workingDirectory)
  set Info(browserPath) ""

  set Info(editorDir) $Info(workingDirectory)
  set Info(editorPath) ""

  # Id counters for creating unique widgets
  #
  set Info(directoryEntryId) 0
  set Info(fileEntryId) 0
  set Info(fileSelectPanelId) 0

  set MatcDefinitions(defaultSaveFile) "./matc.txt"

  #====================================================================================#
  #====================================================================================#
  #====================================================================================#
  #                              FIELD PROPERTIES                                      #  
  #====================================================================================#
  #====================================================================================#
  #====================================================================================#

  #
  #---PROPERTY fields in Common array
  #
  # NOTE for WidgetType: "grouper" are not real field types, they are for grouping:
  #  "data grouper" serves as a common name for all sub-fields
  #  "screen grouper" serves as frame through which sub-fields can be reached
  #
  set Common(WidgetClass)			    01;#-Widget class:
                                      # OutputField
                                      # RadioButtonsContainer MemberOutputField
                                      # ParentField MemberField
                                      # PanelField WorkField
  # NOTE: !!! Remeber to add any new widget type to the widget handling procs:
  # Widget::setFieldStatus
  # Widget::configureField
  #
  set Common(WidgetType)			    02;#-WidgetType:
                                      # -Panel widget types:   
                                      # BrowsableDirectory BrowsableFile IncludeFile
                                      # Button CheckButton
                                      # CheckBox Entry Label Separator ListBox
                                      # OptionMenu RadioButton RadioButtonOptionMenu
                                      # SelectionBox Text
                                      # -Non-panel types:
                                      # GroupData GroupScreen GroupDataAndScreen
  set Common(WidgetPacking)       03;#-Widget packing: Horizontal (-side left), Vertical (-side top)
  set Common(WidgetAnchor)        04;#-Widget anchor in packing: w(est) c(enter) e(ast) n(orth) nw etc
  set Common(WidgetLabel)			    05;#-Widget own label (not used in all cases!)
  set Common(WidgetCommandProc)   06;#-Widget specific command proc (not used in all cases!)
  set Common(WidgetWidth)			    07;#-Widget length (not used in all cases!)
  set Common(WidgetHeight)		    08;#-Widget hight (not used in all cases!)
  set Common(FieldClass)	        09;#-Field class: Normal OutputOnly (ie. not used from model file!)
  set Common(FieldValueType)	    10;#-Value types:
                                     # Integer Real Logical String
                                     # File Directory
                                     # Procedure Library Function Variable
  set Common(FieldDataType)	      11;#-Data type: Scalar Array
  set Common(FieldDataSize) 	    12;#-Data size: 1, N, N 2, N N , DIM, (2*DIM), DIM DIM etc.
  set Common(MinDataSize1) 	      13;#-Data size: possible minimum size1
  set Common(MinDataSize2) 	      14;#-Data size: possible minimum size2
  set Common(MaxDataSize1) 	      15;#-Data size: possible maximum size1
  set Common(MaxDataSize2) 	      16;#-Data size: possible maximum size2
  set Common(FieldDataSep) 	      17;#-Data separator: like space or ';' (default ;).
  set Common(FieldFormat)		      18;#-Field (tcl) format string
  set Common(Label)			          19;#-Field's label on the screen
  set Common(LabelByCoordinate)		20;#-If coordinate label should be added to the label (0/1)
  set Common(LabelWidth)			    21;#-Field's label width on the screen (for aliging the left edges of the widgets!)
  set Common(ScreenPadX)	  	    22;#-Field's padx amount on the screen (default 4)
  set Common(ScreenPadXMembers)	  23;#-Field's padx amount on the screen for group members (default 0)
  set Common(ScreenPadY)	  	    24;#-Field's pady amount on the screen (default 2)
  set Common(ScreenPadYMembers)	  25;#-Field's padx amount on the screen for group members (default 0)
  set Common(ScreenIndent)		    26;#-Field's indent amount on the screen
  set Common(ScreenIndentMembers)	27;#-Group field's memeber indent amount on the screen
  set Common(ScreenMemberArea)	  28;#-1 --> Group members are packed into separate area
  set Common(UnitString)			    29;#-Unit string
  set Common(CoordinateIndex)		  30;#-Coordinate index for the field ): 0=no coordinate index
  set Common(Procedurable)		    31;#-Procedure checkbutton: 0=no, 1=yes
  set Common(Tableable)		        32;#-Procedure checkbutton: 0=no, 1=yes
  set Common(Variabled)		        33;#-If field can have an argument variable (relevan only for Proc/Table)
  set Common(HasAbsCheckBox)		  34;#-Use abs-path checkbox: 0=no, 1=yes
  set Common(HasUseCheckBox)		  35;#-Use field value checkbox: 0=no, 1=yes
  set Common(UseAppendMode)		    36;#-Append mode for directories (like Include path) 0=no, 1=yes
  set Common(Unit)				        37;#-Field's unit (data stored as unit*value)
  set Common(Limits)			        38;#-Limits: Unlimited Interval Set Constant
  set Common(InitialValue)		    39;#-Field inital value
  set Common(AlwaysInitialize)		40;#-Start always with initial value (overwrite "loaded" value)
  set Common(DropInitialValued)		41;#-If field should from from output if it has initial value
  set Common(AlwaysOutput)		    42;#-If field should be always output to sif (even when value is None/False)
  set Common(Display)	            43;#-Display (displayed when active): 1/0
  set Common(InitiallyActive)	    44;#-InitiallyActive: 1/0
  set Common(ActivitySlaves)      45;#-Which field are controlled {list-of-fields}
  set Common(ActivityParents)     46;#-Who is controlling my acitivty {parent-name}
  set Common(ActivityValues)      47;#-What values are controlling my acitivty {list-of-activity-is-true-values}
  set Common(ActivityOnlyByMask)  48;#-Only mask controlled activity: 1/0
  set Common(NeededLevel)				  49;#-If field must be entered when active (0=not needed, 1=needed,but not with a parameter (INCLUDE) file, 2=always needed)
  set Common(PanelPage)			      50;#-Variable's panel page nbr (1 2 ...)
  set Common(SubPanel)			      51;#-Variable's sub  panel nbr (0 1 2 ...) NOTE: 0 is special: not in values area
  set Common(Mask)			          52;#-Field mask to select fields by equation type
                                         # NOTE If two mask:
                                         # 1. mask = when field is displayed (packed) in the panel ("total" mask)
                                         # 2. mask = when field is active in the panel ("case" or "state" mask)
                                         # So ,2. mask specifies more, a field cannot be active if it is not packed!
  set Common(EquationVars)	      53;#-Stores directly sub equation names like Equation(ADVECTION_DIFFUSION_EQUATION_vars) for AD-equation
  set Common(IndexVariable)	      54;#-Indexing variable in EquationVariable-array like ADVECTION_DIFFUSION_EQUATION
  set Common(IndexType)	          55;#-Indexing type: Pre/Post like (Oxygen)DIFFUSIVITY or DIFFUSIVITY(Oxygen)
  set Common(IndexAlsoMask)	      56;#-If indexed variables should create also indexed masks like (¤AD3#1)
  set Common(ChildWidgetClass)		57;#-Widget class for a child field (typically an indexed field)
  set Common(ExcludeFromSif)		  58;#-An otherwise normal output field is NOT ouput to Solver Input File
  set Common(SifName)		          59;#-Specially formed Sif-name (normally not needed!)
                                     #-NOTE: it is defiend however for variables to be used in table/proc entries!
  set Common(OutputSifType)		    60;#-Normally true, but not for ex. Include is without (File) type kwd in Sif
  set Common(TemplateVariable)		61;#-FieldProperty source variable for "indexed" variables like AD_VAR for <Oxygen>AD_VARA
                                     # NOTE: This property is in proc DataField::selectArrayFields
                                     # it is not meant to given explicitely to any field here!
  set Common(FieldTarget)		      62;#-Sif output block for the field. Only accepted value: Solver, used only for Equation panel fields

  #====================================================================================#
  #====================================================================================#
  #====================================================================================#

  # Set default values for each property
  set prop ""
  set vals ""

  lappend prop {WidgetClass WidgetType WidgetPacking WidgetAnchor WidgetLabel WidgetCommandProc WidgetWidth WidgetHeight}
  lappend vals {OutputField Entry Vertical e "" "" 12 1}

  lappend prop {FieldClass FieldValueType FieldDataType}
  lappend vals {Normal Real Scalar}

  lappend prop {FieldDataSize MinDataSize1 MinDataSize2 MaxDataSize1 MaxDataSize2}
  lappend vals {1 "" "" "" ""}

  lappend prop {FieldDataSep FieldFormat}
  lappend vals {";" "% -12.12g"}

  lappend prop {Label LabelByCoordinate LabelWidth}
  lappend vals {"" 0 0}

  lappend prop {ScreenPadX ScreenPadXMembers ScreenPadY ScreenPadYMembers}
  lappend vals {0 0 0 0}

  lappend prop {ScreenIndent ScreenIndentMembers ScreenMemberArea UnitString}
  lappend vals {0 0 1 ""}

  lappend prop {CoordinateIndex Procedurable Tableable Variabled}
  lappend vals {0 1 1 1}

  lappend prop {HasAbsCheckBox HasUseCheckBox UseAppendMode}
  lappend vals {0 0 0}

  lappend prop {Unit Limits InitialValue AlwaysInitialize DropInitialValued AlwaysOutput}
  lappend vals {1 {Unlimited} "" 0 0 0}

  lappend prop {Display InitiallyActive ActivitySlaves ActivityParents ActivityValues ActivityOnlyByMask}
  lappend vals {1 1 "" "" "" 1}

  lappend prop {NeededLevel PanelPage SubPanel Mask}
  lappend vals {0 1 1 ""}

  lappend prop {EquationVars IndexVariable IndexType IndexAlsoMask ChildWidgetClass}
  lappend vals {"" "" "" 1 OutputField}

  lappend prop {ExcludeFromSif SifName OutputSifType FieldTarget}
  lappend vals {0 "" 1 ""}

  # Name of all field properties
  set Common(allFieldProperties) [join $prop]

  # Default values for each property
  set Common(defaultFieldPropertyValues) [join $vals]

  # ====================
  # Set default values #
  # ====================
  #
  DataField::setFieldProperties1 "" default \
     $Common(allFieldProperties) \
     $Common(defaultFieldPropertyValues)

  # =========================================
  # These properties are not by equation!!! #
  # =========================================
  #
  set Common(propertiesNotByEquation) {
    WidgetClass
    WidgetType
    FieldValueType
    FieldDataType
    Mask
  }

  #==============================================================================#
  # Meaning of problem masks which control the state of the entry fields in panels 
  #
  # FL (Navier-Stokes laminar flow)
  # FT1 (Navier-Stokes flow, turbulent_ke)
  # H HT1 HT2 HT3 (Heat, transfer: conduction, constant convection, computed convection)
  # HP1 HP2 HP3(Heat phase change method: Spatial 1, Spatial 2, Temporal)
  # SA SM1 ST1 (Stress analysis: mechanical, thermal)
  # AD AD1 AD2 AD3 (Advection-Diffusion: conduction, constant convection, computed conv.)

  # dimG2 dimG3 (geometry dimension masks)
  # dimS2 dimS3 (simulation dimension masks)
  # nB nP nF nE nV (n=id, for strictly by Body BodyPair Face Edge Vertex )
  # bcD bcN2 bcN3 (Boundary Condition masks: bcD=Dirichlet bcN2=2D Neumann bcN3=3D Neumann) 

  # MG1 (Multigrid solver generated variable)
  #==============================================================================#

  # NOTE: Do not add this to the equation mask!
  # This mask is based on Timestep Transient/Steady State value
  # It is used currently only in the proc: Panel::getCurrentVariables
  #
  set Common(timeMask) "(¤T)"


  # NOTE: first, empty row for index value 0
  # which means: no coordinate dependency
  #
  # None, Cartesian, Axi and Cylindric, Polar
  #
  set Common(CoordinateSymbols) {
    { {}  {}  {}  }
    { (X) (Y) (Z) }
    { (R) (Z) (P) }
    { (R) (T) (P) }
  }

  FRONT_SET_INITIAL_GENERIC_FIELD_PROPERTIES   

  FRONT_SET_INITIAL_PANEL_FIELD_PROPERTIES

  FRONT_UPDATE_FIELD_PROPERTIES

  FRONT_SET_PANEL_SETTINGS

}
### End FRONT_INIT proc




#===================================================#
#                                                   #
#        FRONT SHOW INIT MESSAGES                   #
#                                                   #
#===================================================#
#
# This can be called when UI widget has been build!
#
#
proc FRONT_SHOW_INIT_MESSAGES {} {
  global Info

  # Check succes and display messages
  # =================================

  # If a settings-file was read
  if { $Info(userSettingFilesReadInfo) != "" } {

    foreach info $Info(userSettingFilesReadInfo) {
      Message::showMessage $info $Info(remMsgColor)
    }
  }

  # If something was wrong with the settings-file
  if { $Info(userSettingFilesErrorInfo) != "" } {

    Message::showMessage "WARNING when reading settings file:" $Info(wrnMsgColor)
    foreach info $Info(userSettingFilesErrorInfo) {
      Message::showMessage $info 
    }
    Message::showMessage " "
  }

  # If some command line arguments
  if { $Info(commandLineInfo) != "" } {

    Message::showMessage "Command line arguments:" $Info(remMsgColor)

    foreach info $Info(commandLineInfo) {
      set msg [lindex $info 0]
      set color [lindex $info 1]
      Message::showMessage $msg $color
    }
  }
}
### End FRONT_SHOW_INIT_MESSAGES proc



#===================================================#
#                                                   #
#        READ DEFAULT DEFINITION FILES              #
#                                                   #
#===================================================#
#
# This can be called when UI widget has been build!
#
#
proc READ_DEFAULT_DEFINITION_FILES {} {
  global Info Model

  set Model(defaultDefinitions) ""

  # Read default definition files
  # NOTE: For each file, the result is collected to the global
  # array:  Model(newDefinitions)
  #
  foreach def_file $Model(defaultDefinitionFiles) {
    
    UserDefined::readDefinitionFile $def_file
    
    # Copy added definitions to the deafult definitions
    # log-array
    #
    foreach def_line $Model(newDefinitions) {
      lappend Model(defaultDefinitions) $def_line
    }
  }

}
### End READ_DEFAULT_DEFINITION_FILES proc



#===================================================#
#                                                   #
#             FRONT_POST_INIT                       #
#                                                   #
#===================================================#
#
# This can be called when UI widget has been build!
#
#
proc FRONT_POST_INIT { {read_default_def_file 1 } {read_in_def_file 1} {preset_model_data 1} } {
  global BoundaryCondition Common Equation Model

  # Add user defined equations from default def files
  # =================================================
  if { $read_default_def_file } {
    READ_DEFAULT_DEFINITION_FILES
  }

  # Add loaded user defined equations
  # =================================
  if { $read_in_def_file && $Model(inDefinitionFile) != "" } {
    UserDefined::readDefinitionFile $Model(inDefinitionFile)
  }

  #-Check possible setting files data
  UserSetting::checkData

  # Create "allField"-vars, set some initial field properties
  # =========================================================
  foreach globArray $Common(allPanelArries) {
    upvar #0 $globArray theArray

    #-Work list for field name variables
    # NOTE: (allFields) will be possibly expanded
    # later by indexing variables!
    if { ![info exists theArray(initialFields)] } {
      set theArray(initialFields) ""
    }

    set theArray(allFields) $theArray(initialFields)

    # Initialize some array fields properties
    # NOTE: Indexed fields will get their properties from parent fields
    # which are initialized (if needed) here
    DataField::setDefaultFieldProperties $globArray $theArray(allFields)
  }

  #-Set finally those defaults which have not been defined via arries
  DataField::setDefaultFieldProperties "#" $Common(allFieldNames)

  # BoundaryCondition specific
  # ==========================
  set BoundaryCondition(dirichletVnames) $BoundaryCondition(initialDirichletVnames)
  set BoundaryCondition(neumannVnames) $BoundaryCondition(initialNeumannVnames)


  # Check parameter data
  # ====================
  
  if { [info exist Equation(ids)] &&
       $Equation(ids) != ""
     } {
    #-Create Equation parameter masks
    StdPanelInit::updateEquationDataAndMasks 

    #-Create ObjectTable masks
    StdPanelInit::createObjectTableMasks

    #-Set current problem mask value
    Equation::constructProblemMask

    #-Set problem indices based on body masks
    Equation::updateEquationIndices 

    #-Clean Equation parameters
    set modified1 0
    #Panel::checkParameters Equation modified1

    #-Clean all mask related arries
    set modified2 0
    Panel::checkMaskedParameters modified2
  }

  # Init equation index data
  # ========================
  Equation::setOwnEquationIndices

  if { [info exists Equation(activeIndices)] &&
       $Equation(activeIndices) != ""
     } {
    set modified 0
    Solver::checkSolverData $Equation(activeIndices) modified

  } else {
    Solver::initSolverSystemData
  }

  if { $preset_model_data } {
    MenuExec::presetModelData

  } else {
    
    #-User has possibly added some new constants!
    Panel::initFields Constant
    
    if { $Model(GEOMETRY_DIMENSION) != "" } {
      Constant::panelSave 0
    }

    Panel::initFields Variables
    
    # EquationVariable data needs special handling
    #
    if { $Model(GEOMETRY_DIMENSION) != "" } {
      Panel::initFields EquationVariable
      EquationVariable::updateFieldData
      EquationVariable::formActiveParameter
      Util::cpp_exec equationVariablesPanelOk
      EquationVariable::updateAllVariableNames
    }
  }

	# Load possible SOLVR.KEYWORDS files
	# ==================================
	UserDefined::loadSolverKeywordsFiles

}
### End FRONT_POST_INIT proc




####################################################################################
####################################################################################


#==============================================#
#                                              # 
#  FRONT_SET_INITIAL_GENERIC_FIELD_PROPERTIES  #
#                                              # 
#==============================================#
#
# Generic initial field property definitions #
#
proc FRONT_SET_INITIAL_GENERIC_FIELD_PROPERTIES {} {
  global Common

  #=============#
  # Field names #
  #=============#
  #
  set Common(allFieldNames) {
    ACTIVE
    ACTIVE_IN_MESHING
    ADVECTION_DIFFUSION_EQUATION
    ADVECTION_DIFFUSION_EQUATION_vars
    AUTO_LOAD_MESH
    AUTO_SAVE_MODEL
    AUTO_SAVE_SOLVER_INPUT
    BDF_ORDER
    BM_GEBHARDT_FACTORS_LOG
    BM_GEBHARDT_FACTORS_NONE
    BM_GEBHARDT_FACTORS_SHELL
    BM_MESH_LOG
    BM_MESH_NONE
    BM_MESH_SHELL
    BM_PROCEDURE_COMPILER_LOG
    BM_PROCEDURE_COMPILER_NONE
    BM_PROCEDURE_COMPILER_SHELL
    BM_SOLVER_LOG
    BM_SOLVER_NONE
    BM_SOLVER_SHELL
    BM_VIEW_FACTORS_LOG
    BM_VIEW_FACTORS_NONE
    BM_VIEW_FACTORS_SHELL
    BOUNDARY_LAYER_THICKNESS
    BOUSSINESQ
    BROWSER_COMMAND
    BROWSE_MODE_GEBHARDT_FACTORS
    BROWSE_MODE_MESH
    BROWSE_MODE_PROCEDURE_COMPILER
    BROWSE_MODE_SOLVER
    BROWSE_MODE_VIEW_FACTORS
    BUBBLES
    CHECK_KEYWORDS
    CHECK_LATENT_HEAT_RELEASE
    CLEAR_ENTRY_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_PANELS_BUTTON
    COMPRESSIBILITY_MODEL
    CONVECTION
    CONVECTION_VELOCITY_1
    CONVECTION_VELOCITY_2
    CONVECTION_VELOCITY_3
    CONV_COMPUTED
    CONV_CONSTANT
    CONV_NONE
    COORDINATE
    COORDINATE_1
    COORDINATE_2
    COORDINATE_3
    COORDINATE_MAPPING
    COORDINATE_MAPPING_1
    COORDINATE_MAPPING_2
    COORDINATE_MAPPING_3
    COORDINATE_SYSTEM
    CURRENT_MESH_DIRECTORY
    DEFAULT_AUTO_SAVE_EXTERNAL_MESH
    DEFAULT_CAD_FILES_DIRECTORY
    DEFAULT_INCLUDE_PATH
    DEFAULT_RESULTS_DIRECTORY
    DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY
    DEFAULT_MODEL_DIRECTORY
    DEFAULT_LOG_DIRECTORY
    DEFAULT_USE_MODEL_SETTINGS
    DENSITY
    DIFFUSION
    DIFFUSION_FLUX
    DIFFUSION_SOURCE
    DIFFUSIVITY
    DISPLACEMENT
    DISPLACEMENT_1
    DISPLACEMENT_2
    DISPLACEMENT_3
		ECHO_ON
    EDIT_BUTTON
    EDITOR_COMMAND
    EMISSIVITY
    ENTHALPY
    EQUATION
    EQUATION_LIST
    EXEC_SOLVER
    EXTERNAL_PRESSURE
    EXTERNAL_TEMPERATURE
    FLOW_BODYFORCE_1
    FLOW_BODYFORCE_2
    FLOW_BODYFORCE_3
    FLOW_FORCE_BC
    FONT_SIZES
    FORCE_1
    FORCE_2
    FORCE_3
    FREE_MOVING
    FREE_SURFACE
    FUNCTION_NAME
    GEBHARDT_FACTORS
    GRAVITY
    GRAVITY_1
    GRAVITY_2
    GRAVITY_3
    GRAVITY_ABS
    HEAT_CAPACITY
    HEAT_CONDUCTIVITY
    HEAT_EQUATION
    HEAT_EQUATION_vars
    HEAT_EXPANSION_COEFFICIENT
    HEAT_FLUX
    HEAT_FLUX_BC
    HEAT_SOURCE
    HEAT_TRANSFER_COEFFICIENT
    HYDROSTATIC_PRESSURE
    INCLUDE
    INCLUDE,Browse
    INCLUDE,Use
    INCLUDE_PATH
    KE_C1
    KE_C2
    KE_CLIP
    KE_CMU
    KE_SIGMAE
    KE_SIGMAK
    KE_TURBULENCE
    KE_TURBULENCE_vars
    KINETIC_DISSIPATION
    KINETIC_ENERGY
    LATENT_HEAT
    LIBARY_NAME
    LINEAR_SYSTEM_ABORT_NOT_CONVERGED
    LINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    LINEAR_SYSTEM_DIRECT_METHOD
    LINEAR_SYSTEM_ILUT_TOLERANCE
    LINEAR_SYSTEM_ITERATIVE_METHOD
    LINEAR_SYSTEM_MAX_ITERATIONS
    LINEAR_SYSTEM_METHOD
    LINEAR_SYSTEM_MULTIGRID_METHOD
    LINEAR_SYSTEM_PRECONDITIONING
    LINEAR_SYSTEM_RESIDUAL_OUTPUT
    LINEAR_SYSTEM_SOLVER
    LISTBOX1
    LISTBOX2
    LISTBOX3
    LISTBOX4
    LUMPED_MASS_MATRIX
    MAX_OUTPUT_LEVEL
    MESH
    MESH_BG_MESH
    MESH_BG_MESH_FILE
    MESH_BG_MESH_ACTIVE
    MESH_BG_MESH_CONTROL
    MESH_DENSITY_TYPE
    MESH_DT_H
    MESH_DT_N
    MESH_DT_R
    MESH_DT_NONE
    MESH_ELEMENT_ORDER
    MESH_ELEMENT_TYPE
    MESH_ET_LINE
    MESH_ET_TRIANGLE
    MESH_ET_QUAD
    MESH_ET_TETRA
    MESH_ET_BRICK
    MESH_F
    MESH_H
    MESH_N
    MESH_R
    MESH_INDEX
    MESH_INPUT_FILE
    MESH_NAME
    MESH_NAME_LIST
    MESH_LAYER_TYPE
    MESH_QUADGRID_N1
    MESH_QUADGRID_N2
    MESH_SEED_TYPE
    MESH_ST_EXPLICIT
    MESH_ST_IMPLICIT
    MESH_SEED_EDGE
    MESH_SEED_DIRECTION
    MESH_SEED_VERTEX
    MG_CONVERGENCE_TOLERANCE
    MG_EQUAL_SPLIT
    MG_LEVELS
    MG_MAX_ITERATIONS
    MG_MESH_NAME
    MG_PRECONDITIONING
    MG_PRE_SMOOTHING_ITERATIONS
    MG_POST_SMOOTHING_ITERATIONS
    MG_SMOOTHER
    MIN_OUTPUT_LEVEL
    MODEL_DESCRIPTION
    MODEL_DIRECTORY
    MODEL_NAME
    NAVIER-STOKES
    NAVIER-STOKES_vars
    NEWMARK_BETA
    NOF_PROCESSORS
    NONLINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    NONLINEAR_SYSTEM_MAX_ITERATIONS
    NONLINEAR_SYSTEM_NEWTON_AFTER_ITERATIONS
    NONLINEAR_SYSTEM_NEWTON_AFTER_TOLERANCE
    NONLINEAR_SYSTEM_RELAXATION_FACTOR
    NON_BLOCKING
    NORMAL_FORCE
    NORMAL-TANGENTIAL_DISPLACEMENT
    NORMAL-TANGENTIAL_VELOCITY
    OUTPUT_CALLER
    OUTPUT_FILE
    OUTPUT_INTERVALS
    OUTPUT_LEVEL
    OUTPUT_PREFIX
    PCM_NONE
    PCM_SPATIAL_1
    PCM_SPATIAL_2
    PCM_TEMPORAL
    PERIODIC_BOUNDARIES
    PHASE_CHANGE_INTERVALS
    PHASE_CHANGE_MODEL
    PLANE_STRESS
    POISSON_RATIO
    POST_FILE
    PRESSURE
    PRESSURE_1
    PRESSURE_2
    PRESSURE_3
    PROBLEM_DESCRIPTION
    PROBLEM_NAME
    PROCEDURE
    PROCEDURE_BUTTON
    RADIATION
    RADIATION_BOUNDARIES
    RADIATION_TARGET_BODY
    RADIATION_TARGET_BODIES
    RADIATION_TARGET_BODY_NAME
    RD_DIFFUSE_GRAY
    RD_IDEALIZED
    RD_NONE
    REFERENCE_PRESSURE
    REFERENCE_TEMPERATURE
    RELOAD_INPUT_FILE
    RESTART_FILE
    RESTART_POSITION
    RESULT_MESH_NAME
    RESULTS_DIRECTORY
    SIMULATION_TYPE
    SM_MECHANICAL
    SM_THERMAL
    SOLVER_INPUT_FILE
    SOLVING_ORDER
    SORET_DIFFUSIVITY
    SPECIFIC_HEAT_RATIO
    STABILIZE
    STEADY_STATE_CONVERGENCE_TOLERANCE
    STEADY_STATE_MAX_ITERATIONS
    STEADY_STATE_OUTPUT_INTERVAL
    STEFAN_BOLTZMANN
    STRESS_ANALYSIS
    STRESS_ANALYSIS_vars
    STRESS_BODYFORCE_1
    STRESS_BODYFORCE_2
    STRESS_BODYFORCE_3
    STRESS_MODEL
    SURFACE_ROUGHNESS
    SURFACE_TENSION_COEFFICIENT
    SURFACE_TENSION_EXPANSION_COEFFICIENT
    TABLE_BUTTON
    TEMPERATURE
    LOG_DIRECTORY
    TIME
    TIME_DERIVATIVE_ORDER
    TIMESTEPPING_METHOD
    TIMESTEP_INTERVALS
    TIMESTEP_SIZES
    TM_KE
    TM_NONE
    TURBULENCE_MODEL
    VARIABLE
    VARIABLE_DOFS
    VARIABLE_DOFS_ALL
    VARIABLE_LIST
    VARIABLE_BUTTON
    VELOCITY
    VELOCITY_1
    VELOCITY_2
    VELOCITY_3
    VIEW_FACTORS
    VISCOSITY
    WALL_LAW
    YOUNGS_MODULUS
  }


  # These specific property values are needed for nearly all fields
  #
  set stdProp2    { Label UnitString SubPanel Mask }
  set stdPropLUFM { Label UnitString SubPanel Mask }
  set stdPropLUF  { Label UnitString SubPanel }
  set stdPropLUM  { Label UnitString Mask }
  set stdPropLFM  { Label SubPanel Mask }
  set stdPropLU   { Label UnitString }
  set stdPropLF   { Label SubPanel }
  set stdPropLM   { Label Mask }
  set stdPropL    { Label }
  set stdPropUFM  { UnitString SubPanel Mask }
  set stdPropUM   { UnitString Mask }
  set stdPropFM   { SubPanel Mask }
  set stdPropF    { SubPanel }
  set stdPropM    { Mask }

  namespace import DataField::getFieldPropert* DataField::setFieldPropert*

  #===========================#
  # Field specific properties #
  #===========================#

  set fn ACTIVE
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
      { InitialValue True }
      { ExcludeFromSif 1}
	  }

  set fn ACTIVE_IN_MESHING
  setFieldProperties1 "" $fn $stdPropLF {"Include body in mesh" 0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType CheckBox }
      { WidgetPacking Horizontal }
      { WidgetWidth 0 }
      { Limits {Set {1 0}} }
      { InitialValue 1 }
    }

  set fn ADVECTION_DIFFUSION_EQUATION
  setFieldProperties1 "" $fn $stdPropL {"ADVECTION DIFFUSION EQ."}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue "" }
	  }

  # Stores list of active AddDiff variables in the equation parameter
  # like ADVECTION_DIFFUSION_EQUATION_vars=Oxygen;Nitrogen;...
  #
  set fn ADVECTION_DIFFUSION_EQUATION_vars
  #setFieldProperties1 "" $fn $stdPropM {(¤AD)}
  setFieldProperties2 "" $fn \
    { { WidgetClass OutputField }
      { Display 0 }
      { FieldValueType String }
      { ExcludeFromSif 1 }
    }

  set fn AUTO_LOAD_MESH
  setFieldProperties1 "" $fn $stdPropL {"Auto load mesh"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { Limits {Set {1 0}} }
      { InitialValue 1 }
    }

  set fn AUTO_SAVE_MODEL
  setFieldProperties1 "" $fn $stdPropL {"Auto save model file"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { Limits {Set {1 0}} }
      { ActivityOnlyByMask 0 }
      { InitialValue 1 }
    }

  set fn AUTO_SAVE_SOLVER_INPUT
  setFieldProperties1 "" $fn $stdPropL {"Auto save solver input file"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { Limits {Set {1 0}} }
      { ActivityOnlyByMask 0 }
      { InitialValue 1 }
    }

  set fn BDF_ORDER
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn InitialValue 1

  set fn BOUNDARY_LAYER_THICKNESS
  setFieldProperties1 "" $fn $stdPropLUFM {"Bndr layer thickness" "\[m\]" 1 (¤FT1) }

  set fn BOUSSINESQ
  setFieldProperties1 "" $fn $stdPropLUFM {"Boussinesq" "(Yes/No)" 1 (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
  }

  set fn BROWSER_COMMAND
  setFieldProperties1 "" $fn $stdPropL {"Browser command"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 40 }
      { FieldValueType String }
      { FieldDataSize "N" }
      { FieldFormat "% -39s" }
    }

  set fn BROWSE_MODE_GEBHARDT_FACTORS
  setFieldProperties1 "" $fn $stdPropL {"Gebhardt factors browse:"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { ScreenMemberArea 0}
      { InitialValue "Logfile" }
    }

  set fn BM_GEBHARDT_FACTORS_LOG
  setFieldProperties1 "" $fn $stdPropL {"Log wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Logfile"} }
    }

  set fn BM_GEBHARDT_FACTORS_SHELL
  setFieldProperties1 "" $fn $stdPropL {"Shell wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Shell"} }
    }

  set fn BM_GEBHARDT_FACTORS_NONE
  setFieldProperties1 "" $fn $stdPropL {"None"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "None"} }
    }

  set fn BROWSE_MODE_MESH
  setFieldProperties1 "" $fn $stdPropL {"Mesh generating browse:"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { ScreenMemberArea 0}
      { InitialValue "Logfile" }
    }

  set fn BM_MESH_LOG
  setFieldProperties1 "" $fn $stdPropL {"Log wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Logfile"} }
    }

  set fn BM_MESH_SHELL
  setFieldProperties1 "" $fn $stdPropL {"Shell wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Shell"} }
    }

  set fn BM_MESH_NONE
  setFieldProperties1 "" $fn $stdPropL {"None"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "None"} }
    }

  set fn BROWSE_MODE_PROCEDURE_COMPILER
  setFieldProperties1 "" $fn $stdPropL {"Procedure compiler browse:"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { ScreenMemberArea 0}
      { InitialValue "Logfile" }
    }

  set fn BM_PROCEDURE_COMPILER_LOG
  setFieldProperties1 "" $fn $stdPropL {"Log wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Logfile"} }
    }

  set fn BM_PROCEDURE_COMPILER_SHELL
  setFieldProperties1 "" $fn $stdPropL {"Shell wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Shell"} }
    }

  set fn BM_PROCEDURE_COMPILER_NONE
  setFieldProperties1 "" $fn $stdPropL {"None"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "None"} }
    }

  set fn BROWSE_MODE_SOLVER
  setFieldProperties1 "" $fn $stdPropL {"Solver browse:"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { ScreenMemberArea 0}
      { InitialValue "Logfile" }
    }

  set fn BM_SOLVER_LOG
  setFieldProperties1 "" $fn $stdPropL {"Log wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Logfile"} }
    }

  set fn BM_SOLVER_SHELL
  setFieldProperties1 "" $fn $stdPropL {"Shell wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Shell"} }
    }

  set fn BM_SOLVER_NONE
  setFieldProperties1 "" $fn $stdPropL {"None"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "None"} }
    }

  set fn BROWSE_MODE_VIEW_FACTORS
  setFieldProperties1 "" $fn $stdPropL {"View Factors browse:"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { ScreenMemberArea 0}
      { InitialValue "Logfile" }
    }

  set fn BM_VIEW_FACTORS_LOG
  setFieldProperties1 "" $fn $stdPropL {"Log wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Logfile"} }
    }

  set fn BM_VIEW_FACTORS_SHELL
  setFieldProperties1 "" $fn $stdPropL {"Shell wind."}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Shell"} }
    }

  set fn BM_VIEW_FACTORS_NONE
  setFieldProperties1 "" $fn $stdPropL {"None"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "None"} }
    }

  set fn BUBBLES
  setFieldProperties1 "" $fn $stdPropL {"Bubbles"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn CHECK_KEYWORDS
  setFieldProperties1 "" $fn $stdPropLF {"When unknown sif keyword:" 1 }
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 12 }
      { FieldValueType String }
      { Limits {Set {"Ignore" "Warn" "Abort"}} }
      { InitialValue "Warn" }
      { OutputSifType 0 }
	  }


  set fn CHECK_LATENT_HEAT_RELEASE
  setFieldProperties1 "" $fn $stdPropLFM {"Latent heat release:" 3 {(¤H) (¤HP)} }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
      { ScreenIndent 3 }
	  }

  set fn CLEAR_ENTRY_BUTTON
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType Button }
      { WidgetCommandProc Panel::clearEntryButtonCommandProc }
      { WidgetLabel "Clear entry" }
      { WidgetWidth 0 }
      { WidgetPacking Horizontal }
    }

  set fn CLEAR_PANEL_BUTTON
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType Button }
      { WidgetCommandProc Panel::clearPanelButtonCommandProc }
      { WidgetLabel "Clear panel" }
      { WidgetWidth 0 }
      { WidgetPacking Horizontal }
    }

  set fn CLEAR_PANELS_BUTTON
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType Button }
      { WidgetCommandProc Panel::clearPanelsButtonCommandProc }
      { WidgetLabel "Clear all" }
      { WidgetWidth 0 }
      { WidgetPacking Horizontal }
    }

  set fn CLEAR_UNDO
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType Button }
      { WidgetCommandProc Panel::undoClearButtonCommandProc }
      { WidgetLabel "Undo clear" }
      { WidgetWidth 0 }
      { WidgetPacking Horizontal }
    }

  set fn COMPRESSIBILITY_MODEL
  setFieldProperties1 "" $fn $stdPropLFM {"Compressibility model" 2 (¤FL)|(¤HT2)|(¤HT3)|(¤AD2)|(¤AD3) }
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 18 }
      { FieldValueType String }
      { Limits {Set {"Incompressible" "Perfect Gas Equation 1" } } }
      { InitialValue "Incompressible" }
      { DropInitialValued 1 }
    }

  # This list is not in use!
  #    { Limits {Set {"Incompressible" "User Defined 1"  "User Defined 2"
  #                   "Perfect Gas Equation 1" "Perfect Gas Equation 2" "Perfect Gas Equation 3" } } }

  set fn CONVECTION
  setFieldProperties1 "" $fn $stdPropLFM {"Convection:" 2 (¤H)|(¤AD) }
  setFieldProperties2 "" $fn \
    { { WidgetClass  RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { WidgetPacking Vertical }
      { FieldValueType String }
      { InitialValue "None" }
	  }

  set fn CONV_NONE
  setFieldProperties1 "" $fn $stdPropLFM {"None" 2 (¤H)|(¤AD) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 12 }
      { Limits {Constant "None"} }
      { ScreenIndent 6 }
	  }

  set fn CONV_COMPUTED
  setFieldProperties1 "" $fn $stdPropLFM {"Computed" 2 (¤H)|(¤AD) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 12 }
      { Limits {Constant "Computed"} }
      { ScreenIndent 6 }
	  }

  set fn CONV_CONSTANT
  setFieldProperties1 "" $fn $stdPropLFM {"Constant" 2 (¤H)|(¤AD) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 12 }
      { Limits {Constant "Constant"} }
      { ScreenIndent 6 }
	  }

  set fn CONVECTION_VELOCITY_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUFM {"Convection Velocity" "\[m/s\]" 2 {(¤H)|(¤AD) (¤HT2)|(¤AD2)} }
    setFieldProperty    "" $fn$i SifName "Convection Velocity $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
    setFieldProperty    "" $fn$i InitialValue 0
  }
  setFieldProperty  "" CONVECTION_VELOCITY_2 Mask {(¤H)|(¤AD) ((¤HT2)|(¤AD2))&((¤dimS2)|(¤dimS3))}
  setFieldProperty  "" CONVECTION_VELOCITY_3 Mask {(¤H)|(¤AD) ((¤HT2)|(¤AD2))&(¤dimS3)}

  set fn COORDINATE
  setFieldProperties1 "" $fn $stdPropLUM {"Coordinate" "\[m\]" ""}
	setFieldProperty "" $fn FieldDataSize "DIM"
  setFieldProperty "" $fn FieldDataType Array

  set fn COORDINATE_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUF {"Coordinate" "\[m\]" 0}
    setFieldProperty    "" $fn$i SifName "Coordinate $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" COORDINATE_2 Mask ((¤dimS2)|(¤dimS3))
  setFieldProperty    "" COORDINATE_3 Mask (¤dimS3)


  set fn COORDINATE_MAPPING
  setFieldProperty "" $fn FieldDataType Array
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn InitialValue {1 2 3}

  set fn COORDINATE_SYSTEM
  setFieldProperty "" $fn FieldValueType String
  setFieldProperty "" $fn InitialValue "None"

  set fn CURRENT_MESH_DIRECTORY
  setFieldProperty "" $fn FieldValueType File

  set fn DEFAULT_MODEL_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"Model directory"}
  setFieldProperties2 "" $fn \
    { { FieldClass OutputOnly } 
      { WidgetWidth 40 }
      { WidgetType  BrowsableDirectory }
      { FieldValueType String }
      { FieldFormat "% -39s" }
      { InitialValue "\$Info(ELMER_MODEL_DIRECTORY)" }
    }

  set fn DEFAULT_CAD_FILES_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"Cad files directory"}
  setFieldProperties2 "" $fn \
    { { FieldClass OutputOnly } 
      { WidgetWidth 40 }
      { WidgetType  BrowsableDirectory }
      { FieldValueType String }
      { FieldFormat "% -39s" }
    }

  set fn DEFAULT_INCLUDE_PATH
  setFieldProperties1 "" $fn $stdPropL {"Include path"}
  setFieldProperties2 "" $fn \
    { { FieldClass OutputOnly } 
      { WidgetWidth 40 }
      { WidgetType  BrowsableDirectory }
      { UseAppendMode 1 }
      { FieldValueType String }
      { FieldFormat "% -39s" }
      { InitialValue "\$Info(ELMER_INCLUDE_PATH)" }
    }

  set fn DEFAULT_LOG_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"Log directory"}
  setFieldProperties2 "" $fn \
    { { FieldClass OutputOnly } 
      { WidgetWidth 40 }
      { WidgetType  BrowsableDirectory }
      { FieldValueType String }
      { FieldFormat "% -39s" }
      { InitialValue "\$Info(LOG_DIRECTORY)" }
    }

  set fn DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"External mesh files directory"}
  setFieldProperties2 "" $fn \
    { { FieldClass OutputOnly } 
      { WidgetWidth 40 }
      { WidgetType  BrowsableDirectory }
      { FieldValueType String }
      { FieldFormat "% -39s" }
    }

  set fn DEFAULT_RESULTS_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"Results directory"}
  setFieldProperties2 "" $fn \
    { { FieldClass OutputOnly } 
      { WidgetWidth 40 }
      { WidgetType  BrowsableDirectory }
      { FieldValueType String }
      { FieldFormat "% -39s" }
      { InitialValue "\$Info(ELMER_RESULTS_DIRECTORY)" }
    }

  set fn DEFAULT_USE_MODEL_SETTINGS
  setFieldProperties1 "" $fn $stdPropL {"Apply model file settings"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { Limits {Set {1 0}} }
      { InitialValue 1 }
    }

  set fn DEFAULT_AUTO_SAVE_EXTERNAL_MESH
  setFieldProperties1 "" $fn $stdPropL {"Auto save external mesh in Elmer format"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { Limits {Set {1 0}} }
      { InitialValue 1 }
    }

  set fn DENSITY
  setFieldProperties1 "" $fn $stdPropLUM {"Density" "\[kg/m3]" (¤FL)|(¤H)|(¤S)|(¤AD) }
  setFieldProperty    "" $fn  NeededLevel 1
  	
  # NOTE: SifName="" --> creates names like Oxygen
  set fn DIFFUSION
  setFieldProperties1 "" $fn $stdPropUM {"\[kg/m3\]" (¤AD) }
  setFieldProperty    "" $fn SifName ""
  setFieldProperty    "" $fn IndexType "Pre"
  setFieldProperty    "" $fn IndexVariable ADVECTION_DIFFUSION_EQUATION
  setFieldProperty    "" $fn WidgetClass ParentField

  # NOTE: SifName="Flux" --> creates names like Oxygen Flux
  set fn DIFFUSION_FLUX
  setFieldProperties1 "" $fn $stdPropLUFM {"Flux" "\[kg/m2s\] " 2 (¤AD) }
  setFieldProperty    "" $fn SifName "Flux"
  setFieldProperty    "" $fn IndexType "Pre"
  setFieldProperty    "" $fn IndexVariable ADVECTION_DIFFUSION_EQUATION
  setFieldProperty    "" $fn WidgetClass ParentField

  # NOTE: SifName="Source" --> creates names like Oxygen Source
  set fn DIFFUSION_SOURCE
  setFieldProperties1 "" $fn $stdPropLUM {"Source" "\[kg/m3s\]" (¤AD) }
  setFieldProperty    "" $fn SifName "Source"
  setFieldProperty    "" $fn IndexType "Pre"
  setFieldProperty    "" $fn IndexVariable ADVECTION_DIFFUSION_EQUATION
  setFieldProperty    "" $fn WidgetClass ParentField


  # NOTE: SifName="Diffusivity" --> creates names like Oxygen Diffusivity
  set fn DIFFUSIVITY
  setFieldProperties1 "" $fn $stdPropLUM {"Diffusivity" "\[kg/m\]" (¤AD) }
  setFieldProperty    "" $fn SifName "Diffusivity"
  setFieldProperties2 "" $fn \
    { {WidgetClass ParentField}
      {FieldDataSize "DIM DIM"}
      {IndexType "Pre"}
      {IndexVariable ADVECTION_DIFFUSION_EQUATION }
      {NeededLevel 1}
    }

  set fn DISPLACEMENT
  setFieldProperties1 "" $fn $stdPropLUM {"Displacement" "\[m\]" (¤S)}
	setFieldProperty "" $fn FieldDataSize "DIM"
  setFieldProperty "" $fn FieldDataType Array

  set fn DISPLACEMENT_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUM {"Displacement" "\[m\]" (¤S) }
    setFieldProperty    "" $fn$i SifName "Displacement $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" DISPLACEMENT_2 Mask {(¤S) (¤S)&((¤dimS2)|(¤dimS3))}
  setFieldProperty    "" DISPLACEMENT_3 Mask {(¤S) (¤S)&(¤dimS3)}

  set fn ECHO_ON
  setFieldProperties1 "" $fn $stdPropLF {"Echo solver input file:" 1 }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
      { OutputSifType 0 }
	  }


  set fn EDIT_BUTTON
  setFieldProperties1 "" $fn $stdPropLF {"   " 0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType Button }
      { WidgetCommandProc Panel::editButtonCommandProc }
      { WidgetLabel "Edit" }
      { WidgetWidth 0 }
      { WidgetPacking Horizontal }
    }

  set fn EDITOR_COMMAND
  setFieldProperties1 "" $fn $stdPropL {"Editor command"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 40 }
      { FieldValueType String }
      { FieldFormat "% -39s" }
    }

  set fn EMISSIVITY
  setFieldProperties1 "" $fn $stdPropLUFM {"Emissivity" "(0...1)" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetWidth 10 }
      { Limits {Interval {0.0 1.0} } }
      { NeededLevel 2 }
      { ActivityOnlyByMask 0 }
      { ActivityParents RADIATION }
      { ActivityValues {"Diffuse gray"} }
    }

  set fn ENTHALPY
  setFieldProperties1 "" $fn $stdPropLUM {"Enthalpy" "\[J/m3]" (¤H) }

  set fn EQUATION
  setFieldProperty "" $fn FieldValueType String

  set fn EQUATION_LIST
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType SelectionBox }
      { WidgetLabel "Active equations:" }
      { WidgetHeight 8 }
      { WidgetWidth 24 }
      { WidgetPacking Vertical }
      { Limits {Set {@Equation::updateEquationList}} }
      { ActivityOnlyByMask 0 }
     }

  set fn EXEC_SOLVER
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 15 }
      { FieldValueType String }
      { Limits {Set { "Always" "Before Simulation" "After Simulation"
                      "Before Timestep" "After Timestep" "Never" } }}
    }

  set fn EXTERNAL_PRESSURE
  setFieldProperties1 "" $fn $stdPropLUFM {"External pressure" "\[Pa\]" 2 (¤FL) }

  set fn EXTERNAL_TEMPERATURE
  setFieldProperties1 "" $fn $stdPropLUM {"External temperature" "\[K\]" (¤H) }

  set fn FLOW_BODYFORCE_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUM {"Flow Bodyforce" "\[N/m3]" (¤FL) }
    setFieldProperty    "" $fn$i SifName "Flow Bodyforce $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" FLOW_BODYFORCE_2 Mask {(¤FL) (¤FL)&((¤dimS2)|(¤dimS3))}
  setFieldProperty    "" FLOW_BODYFORCE_3 Mask {(¤FL) (¤FL)&(¤dimS3)}

  set fn FLOW_FORCE_BC
  setFieldProperties1 "" $fn $stdPropFM {0 (¤FL) }
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
	  }

  set fn FONT_SIZES
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn FieldDataType Array
  setFieldProperty "" $fn FieldDataSize 3

  set fn FORCE_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUFM {"Force" "\[N/m2\]" 2 (¤S) }
    setFieldProperty    "" $fn$i SifName "Force $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" FORCE_2 Mask {(¤S) (¤S)&((¤dimS2)|(¤dimS3))}
  setFieldProperty    "" FORCE_3 Mask {(¤S) (¤S)&(¤dimS3)}

  set fn FREE_MOVING
  setFieldProperties1 "" $fn $stdPropLUM {"Free moving" "(Yes/No)" (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn FREE_SURFACE
  setFieldProperties1 "" $fn $stdPropLUM {"Free surface" "(Yes/No)" (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn FUNCTION_NAME
  setFieldProperties1 "" $fn $stdPropL "Function name"
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField }
      { FieldValueType Function }
      { FieldDataSize 1 }
    }

  set fn GRAVITY 
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldDataType Array }
      { InitialValue {0 -1 0 9.82} }
    }

  set fn GRAVITY_1
  setFieldProperties1 "" $fn $stdPropL {"X-direction"}
  setFieldProperty "" $fn WidgetClass PanelField
  setFieldProperty "" $fn Procedurable 0

  set fn GRAVITY_2
  setFieldProperties1 "" $fn $stdPropL {"Y-direction"}
  setFieldProperty "" $fn WidgetClass PanelField
  setFieldProperty "" $fn Procedurable 0

  set fn GRAVITY_3
  setFieldProperties1 "" $fn $stdPropL {"Z-direction"}
  setFieldProperty "" $fn WidgetClass PanelField
  setFieldProperty "" $fn Procedurable 0

  set fn GRAVITY_ABS
  setFieldProperties1 "" $fn $stdPropLU {"Intensity" "\[m/s2\]"}
  setFieldProperty "" $fn WidgetClass PanelField
  setFieldProperty "" $fn Procedurable 0

  set fn GEBHARDT_FACTORS
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue "GebhardtFactors.dat"

  # Note: This is needed also for compressible flow!
  # Activity is controled in Material(COMP...) edit and fillProcs
  set fn HEAT_CAPACITY
  setFieldProperties1 "" $fn $stdPropLUFM {"Heat capacity" "\[J/kgK]" { {#0 2} {#2 1} {#4 2} } (¤FL)|(¤H)|(¤AD2)|(¤AD3) }
  setFieldProperties2 "" $fn \
    { { ActivityOnlyByMask 0 }
    }

  set fn HEAT_CONDUCTIVITY
  setFieldProperties1 "" $fn $stdPropLUM {"Heat conductivity" "\[W/mK]" (¤H) }
  setFieldProperties2 "" $fn \
    { { FieldDataSize "DIM DIM" }
      { NeededLevel 1 }
	  }

  set fn HEAT_EQUATION
  setFieldProperties1 "" $fn $stdPropL {"HEAT EQUATION"}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
	  }

  set fn HEAT_EQUATION_vars
  setFieldProperties1 "" $fn $stdPropM {(¤H)}
  setFieldProperties2 "" $fn \
    { { WidgetClass OutputField }
      { Display 0 }
      { FieldValueType String }
      { ExcludeFromSif 1 }
    }

  set fn HEAT_EXPANSION_COEFFICIENT
  setFieldProperties1 "" $fn $stdPropLUM {"Heat expansion coeff." "\[1/K\]" (¤FL) }

  set fn HEAT_FLUX
  setFieldProperties1 "" $fn $stdPropLUM {"Heat flux" "\[W/m2\]" (¤H) }

  set fn HEAT_FLUX_BC
  setFieldProperties1 "" $fn $stdPropFM {0 (¤H) }
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
	  }

  set fn HEAT_SOURCE
  setFieldProperties1 "" $fn $stdPropLUM {"Heat Source" "\[W/kg]" (¤H) }

  set fn HEAT_TRANSFER_COEFFICIENT
  setFieldProperties1 "" $fn $stdPropLUM {"Heat transfer coeff." "\[W/m2K\]" (¤H) }


  set fn HYDROSTATIC_PRESSURE
  setFieldProperties1 "" $fn $stdPropLFM {"Calc hydrostatic pressure" 2 (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
	  }

  set fn INCLUDE
  setFieldProperties1 "" $fn $stdPropLF {"Include file" 0}
  setFieldProperties2 "" $fn \
    { { WidgetType IncludeFile }
      { WidgetWidth 40 }
      { FieldValueType File }
      { OutputSifType 0 }
      { FieldFormat "% -39s" }
      { ActivityOnlyByMask 0 }
      { InitiallyActive 0 }
      { HasUseCheckBox 1 } 
      { Procedurable 0 } 
      { Tableable 0 } 
      { NeededLevel 3 }
    }

  set fn INCLUDE,Browse
  setFieldProperties2 "" $fn \
    { { WidgetClass WorkField }
      { WidgetType Button }
    }

  set fn INCLUDE,Use
  setFieldProperties2 "" $fn \
    { { WidgetClass WorkField }
      { FieldValueType Logical }
      { InitialValue 0 }
      { AlwaysInitialize 1 }
      { WidgetType CheckBox }
      { ActivityOnlyByMask 0 }
      { ExcludeFromSif 1 }
    }

  set fn INCLUDE_PATH
  setFieldProperties1 "" $fn $stdPropL {"Include path"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 40 }
      { FieldValueType String }
      { FieldFormat "% -39s" }
    }

  set fn KE_TURBULENCE
  setFieldProperties1 "" $fn $stdPropFM {0 (¤FT1) }
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue "" }
	  }

  set fn KE_TURBULENCE_vars
  setFieldProperties1 "" $fn $stdPropM {(¤FT1)}
  setFieldProperties2 "" $fn \
    { { WidgetClass OutputField }
      { Display 0 }
      { FieldValueType String }
      { ExcludeFromSif 1 }
    }

  set fn KE_CLIP
  setFieldProperties1 "" $fn $stdPropLFM {"KE Clip" 2 {(¤FL) (¤FT1)} }
  setFieldProperties2 "" $fn \
    { { InitialValue 1.0e-06 }
      { NeededLevel 1 }
      { ScreenIndent 6 }
	  }

  set fn KE_C1
  setFieldProperties1 "" $fn $stdPropLM {"KE C1" (¤FT1) }
  setFieldProperty    "" $fn InitialValue 1.44

  set fn KE_C2
  setFieldProperties1 "" $fn $stdPropLM {"KE C2" (¤FT1) }
  setFieldProperty    "" $fn InitialValue 1.92

  set fn KE_CMU
  setFieldProperties1 "" $fn $stdPropLM {"KE Cmu" (¤FT1) }
  setFieldProperty    "" $fn InitialValue 0.09

  set fn KE_SIGMAE
  setFieldProperties1 "" $fn $stdPropLM {"KE SigmaE" (¤FT1) }
  setFieldProperty    "" $fn InitialValue 1.30
  setFieldProperty    "" $fn SifName "KE SigmaE"

  set fn KE_SIGMAK
  setFieldProperties1 "" $fn $stdPropLM {"KE SigmaK" (¤FT1) }
  setFieldProperty    "" $fn InitialValue 1.00
  setFieldProperty    "" $fn SifName "KE SigmaK"

  set fn KINETIC_ENERGY
  setFieldProperties1 "" $fn $stdPropLUM { "Kinetic energy" "\[J\]" (¤FT1) }
  setFieldProperty    "" $fn SifName "Kinetic energy"

  set fn KINETIC_DISSIPATION
  setFieldProperties1 "" $fn $stdPropLUM {"Kinetic dissipation" "\[J\]" (¤FT1) }
  setFieldProperty    "" $fn SifName "Kinetic dissipation"

  set fn LATENT_HEAT
  setFieldProperties1 "" $fn $stdPropLUFM {"Latent heat" "\[J/m3\]" 2 (¤HP) }

  set fn LIBRARY_FILE
  setFieldProperties1 "" $fn $stdPropL "Library file"
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField }
      { WidgetType  BrowsableFile }
      { HasAbsCheckBox 1 }
      { FieldValueType Library }
      { FieldDataSize 1 }
    }

  set fn LINEAR_SYSTEM_ABORT_NOT_CONVERGED
  setFieldProperties1 "" $fn $stdPropLF {"Abort when not converged" 2}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
    }

  set fn LINEAR_SYSTEM_CONVERGENCE_TOLERANCE
  setFieldProperties1 "" $fn $stdPropLF {"Tolerance" 2}
  setFieldProperty "" $fn Procedurable 0

  set fn LINEAR_SYSTEM_DIRECT_METHOD
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType String }
      { Limits {Set { "BANDED" "UMFPack" } }}
    }

  set fn LINEAR_SYSTEM_ITERATIVE_METHOD
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType String }
      { Limits {Set { "BiCGStab" "TFQMR" "CGS" "CG" "GMRES" } }}
    }

  set fn LINEAR_SYSTEM_MAX_ITERATIONS
  setFieldProperties1 "" $fn $stdPropLF {"Max iterations" 2}
  setFieldProperty "" $fn FieldValueType Integer

  set fn LINEAR_SYSTEM_METHOD
  setFieldProperties1 "" $fn $stdPropLF {"Solving method" 2}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField }
      { WidgetType OptionMenu }
      { WidgetWidth 12 }
      { FieldValueType String }
    }

  # NOTE: This field is needed to store a different method list
  # for the multigrid solver, it is output as "Iterative Method" in Sif!
  #
  set fn LINEAR_SYSTEM_MULTIGRID_METHOD
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { SifName "Linear System Iterative Method" }
      { FieldValueType String }
      { Limits {Set { "Jacobi" "CG" "BiCGStab" } }}
    }

  set fn LINEAR_SYSTEM_PRECONDITIONING
  setFieldProperties1 "" $fn $stdPropLF {"Preconditioning" 2}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 10 }
      { FieldValueType String }
      { Limits {Set {
        "None" "Diagonal" "Multigrid"
        "ILU0" "ILU1" "ILU2" "ILU3" "ILU4" "ILU5" "ILU6" "ILU7" "ILU8" "ILU9" "ILUT"} }}
    }

  set fn LINEAR_SYSTEM_ILUT_TOLERANCE
  setFieldProperties1 "" $fn $stdPropLF {"ILUT Tolerance" 2}
  setFieldProperty "" $fn InitiallyActive 0
  setFieldProperty "" $fn ActivityOnlyByMask 0
  setFieldProperty "" $fn Procedurable 0

  set fn LINEAR_SYSTEM_RESIDUAL_OUTPUT
  setFieldProperties1 "" $fn $stdPropLF {"Residual output" 2}
  setFieldProperty "" $fn FieldValueType Integer

  set fn LINEAR_SYSTEM_SOLVER
  setFieldProperties1 "" $fn $stdPropLF {"Solver type" 2}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 12 }
      { FieldValueType String }
      { Limits {Set { "Direct" "Iterative" "Multigrid"} }}
    }

  set fn LISTBOX1
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetType ListBox }
      { WidgetHeight 6 }
      { WidgetWidth 25 }
      { InitialValue "" }
    }

  set fn LISTBOX2
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetType ListBox }
      { WidgetHeight 6 }
      { WidgetWidth 25 }
      { InitialValue "" }
    }

  set fn LISTBOX3
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetType ListBox }
      { WidgetHeight 6 }
      { WidgetWidth 25 }
      { InitialValue "" }
    }

  set fn LISTBOX4
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetType ListBox }
      { WidgetHeight 6 }
      { WidgetWidth 25 }
      { InitialValue "" }
    }

  set fn LOG_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"Log directory"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 40 }
      { FieldValueType String }
      { FieldFormat "% -39s" }
    }

  set fn LUMPED_MASS_MATRIX
  setFieldProperties1 "" $fn $stdPropL {"Lumped mass matrix"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
    }

  set fn MAX_OUTPUT_LEVEL
  setFieldProperties1 "" $fn $stdPropLF {"Max info output level" 1}
  setFieldProperties2 "" $fn \
    { { FieldValueType Integer }
      { Limits {Interval {0 31}} }
      { InitialValue 31 }
    }

  set fn MESH
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { FieldValueType File }
    }

  set fn MESH_BG_MESH
  setFieldProperties1 "" $fn $stdPropL {"Background mesh"}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { ScreenPadY 2 }
      { WidgetWidth 16 }
      { FieldValueType String }
      { Limits {Set {"Delaunay" "Grid" "SparseGrid" "DenseGrid" "External"} } }
      { InitialValue "Delaunay" }
    }

  set fn MESH_BG_MESH_FILE
  setFieldProperties1 "" $fn $stdPropL {"Bg mesh file"}
  setFieldProperties2 "" $fn \
    { { WidgetType BrowsableFile }
      { HasAbsCheckBox 1 }
      { ScreenPadY 2 }
      { WidgetWidth 15 }
    }

  set fn MESH_BG_MESH_ACTIVE
  setFieldProperties1 "" $fn $stdPropL {"Use bg mesh file"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { InitialValue 0 }
    }

  set fn MESH_BG_MESH_CONTROL
  setFieldProperties1 "" $fn $stdPropL {"Use as control file"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { InitialValue 0 }
    }

  set fn MESH_DENSITY_TYPE
  setFieldProperties1 "" $fn $stdPropL {"Mesh density value type"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { InitialValue "R" }
    }

  set fn MESH_DT_H
  setFieldProperties1 "" $fn $stdPropL {"Mesh H"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 6 }
      { Limits {Constant "H"} }
    }

  set fn MESH_DT_R
  setFieldProperties1 "" $fn $stdPropL {"Relative"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 6 }
      { Limits {Constant "R"} }
    }

  set fn MESH_DT_N
  setFieldProperties1 "" $fn $stdPropL {"Nof elements"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 9 }
      { Limits {Constant "N"} }
    }

  set fn MESH_DT_NONE
  setFieldProperties1 "" $fn $stdPropL {"None"}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 4 }
      { Limits {Constant "None"} }
    }

  set fn MESH_F
  setFieldProperties1 "" $fn $stdPropL {"Mesh F"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 16 }
      { InitialValue 1.0 }
      { NeededLevel 3 }
      { Limits {Interval {0 Inf}} }
      { Procedurable 0 } 
      { Tableable 0 } 
    }

  set fn MESH_H
  setFieldProperties1 "" $fn $stdPropLU {"Mesh H" "\[m\]"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 16 }
      { InitialValue "\$MeshDefine(MESH_H,default)" }
      { NeededLevel 3 }
		  { FieldDataSize N }
		  { MaxDataSize1 2 }
      { Limits {Interval {0 Inf}} }
      { Procedurable 0 } 
      { Tableable 0 } 
      { ActivityOnlyByMask 0 }
    }

  set fn MESH_N
  setFieldProperties1 "" $fn $stdPropL {"Nof elements"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 16 }
      { InitialValue "" }
      { NeededLevel 3 }
      { Limits {Interval {(0 Inf)}} }
      { FieldFormat "%16d" }
      { Procedurable 0 } 
      { Tableable 0 } 
      { ActivityOnlyByMask 0 }
    }

  set fn MESH_R
  setFieldProperties1 "" $fn $stdPropL {"Relative mesh H"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 16 }
      { InitialValue 1.0 }
      { NeededLevel 3 }
      { FieldDataType Array }
		  { FieldDataSize N 1}
		  { MaxDataSize1 2 }
      { Procedurable 0 } 
      { Tableable 0 } 
      { ActivityOnlyByMask 0 }
    }

  set fn MESH_ELEMENT_ORDER
  setFieldProperties1 "" $fn $stdPropL {"Element order"}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { ScreenPadY 2 }
      { WidgetWidth 8 }
      { FieldValueType String }
      { Limits {Set {"Linear" "Parabolic" } } }
      { InitialValue "Linear" }
    }

  set fn MESH_ELEMENT_TYPE
  setFieldProperties1 "" $fn $stdPropL {"Element type"}
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { InitialValue "Triangle" }
    }

  set fn MESH_ET_LINE
  setFieldProperties1 "" $fn $stdPropLM {"Line" (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 4 }
      { Limits {Constant "Line"} }
    }

  set fn MESH_ET_TRIANGLE
  setFieldProperties1 "" $fn $stdPropLM {"Triangle" (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Triangle"} }
    }

  set fn MESH_ET_QUAD
  setFieldProperties1 "" $fn $stdPropLM {"Quad" (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 4 }
      { Limits {Constant "Quad"} }
    }

  set fn MESH_ET_TETRA
  setFieldProperties1 "" $fn $stdPropLM {"Tetra" (¤dimG3) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 5 }
      { Limits {Constant "Tetra"} }
    }

  set fn MESH_ET_BRICK
  setFieldProperties1 "" $fn $stdPropLM {"Brick" (¤dimG3) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 5 }
      { Limits {Constant "Brick"} }
    }


  set fn MESH_INDEX
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { WidgetClass PanelField } 
      { InitialValue $Info(NO_INDEX) }
    }

  set fn MESH_INPUT_FILE
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue "\@Datafile::getDefaultMeshInputFile"
  setFieldProperty "" $fn ExcludeFromSif 1

  set fn MESH_NAME
  setFieldProperties1 "" $fn $stdPropLF {"Mesh name" 0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetWidth 24 }
      { FieldValueType String }
      { ScreenPadY 3 }
      { InitiallyActive 0 }
      { InitialValue "" }
      { Procedurable 0 } 
      { Tableable 0 } 
    }

  set fn MESH_NAME_LIST
  setFieldProperties1 "" $fn $stdPropLF {"Select mesh" 0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType OptionMenu }
      { WidgetWidth 9 }
      { FieldValueType String }
      { Limits {Set {$Model(meshNames)} } }
      { InitialValue "Mesh1" }
    }

  set fn MESH_LAYER_TYPE
  setFieldProperties1 "" $fn $stdPropL {"Meshing method"}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { ScreenPadY 3 }
      { WidgetWidth 16 }
      { WidgetPacking Vertical }
      { FieldValueType String }
      { Limits {Set { 
                  {"BoundaryMesh"} 
                  {"MovingFront" "VoronoiVertex" "SSSFMovingFront" "SSMFMovingFront"} 
                  {"TriangleNEGrid" "TriangleNWGrid" "TriangleUJNEGrid" "TriangleUJNWGrid" "TriangleFBNEGrid" "TriangleFBNWGrid"} 
                  {"QuadGrid"} 
                  {"TetraOctree"}
                  {"BrickOctree"} }}}
      { InitialValue "MovingFront" }
    }

  set fn MESH_QUADGRID_N1
  setFieldProperties1 "" $fn $stdPropLFM {"Nof elements (1st and 3rd edge)" 2 {(¤dimG2) (¤dimG2)} }
  setFieldProperties2 "" $fn \
    { { WidgetWidth 12 }
      { ScreenPadY 8 }
      { NeededLevel 3 }
      { Limits {Interval {(0 Inf)} } }
      { FieldFormat "%12d" }
      { Procedurable 0 } 
      { Tableable 0 } 
      { ActivityOnlyByMask 0 }
    }

  set fn MESH_QUADGRID_N2
  setFieldProperties1 "" $fn $stdPropLFM {"Nof elements (2nd and 4th edge)" 2 {(¤dimG2) (¤dimG2)} }
  setFieldProperties2 "" $fn \
    { { WidgetWidth 12 }
      { ScreenPadY 8 }
      { NeededLevel 3 }
      { Limits {Interval {(0 Inf)} } }
      { FieldFormat "%12d" }
      { Procedurable 0 } 
      { Tableable 0 } 
      { ActivityOnlyByMask 0 }
    }

  set fn MESH_SEED_EDGE
  setFieldProperties1 "" $fn $stdPropLFM {"Meshing seed edge" 2 (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { FieldValueType Integer }
      { WidgetWidth 20 }
      { ScreenPadY 8 }
      { NeededLevel 3 }
      { Procedurable 0 }
      { InitialValue "" }
      { ActivityOnlyByMask 0 }
    }

  # Not in use
  set fn MESH_SEED_DIRECTION
  setFieldProperties1 "" $fn $stdPropLUFM {"Direction" "(x,y)" 2 (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { FieldDataType Array }
      { FieldDataSize 3 }
      { WidgetWidth 20 }
      { Procedurable 0 }
      { InitialValue "" }
      { ActivityOnlyByMask 0 }
    }

  # Not in use
  set fn MESH_SEED_VERTEX
  setFieldProperties1 "" $fn $stdPropLUFM {"Inner node" "(x,y)" 2 (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { FieldDataType Array }
      { FieldDataSize 3 }
      { ScreenPadY 5 }
      { WidgetWidth 20 }
      { Procedurable 0 }
      { InitialValue "" }
      { ActivityOnlyByMask 0 }
    }

  set fn MESH_SEED_TYPE
  setFieldProperties1 "" $fn $stdPropLFM {"Mesh seed type" 2 (¤dimG2) }
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { InitialValue "Implicit" }
      { ActivityOnlyByMask 0 }
    }

  # Not in use
  set fn MESH_ST_EXPLICIT
  setFieldProperties1 "" $fn $stdPropLF {"Inner node" 2}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 8 }
      { Limits {Constant "Explicit"} }
    }

  # Not in use
  set fn MESH_ST_IMPLICIT
  setFieldProperties1 "" $fn $stdPropLF {"Boundary" 2}
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Implicit"} }
    }

  set fn MG_EQUAL_SPLIT
  setFieldProperties1 "" $fn $stdPropLF {"MG equal split" 2}
  setFieldProperty "" $fn WidgetType CheckBox
  setFieldProperty "" $fn FieldValueType Logical

  set fn MG_LEVELS
  setFieldProperties1 "" $fn $stdPropLF {"MG levels" 2}
  setFieldProperty "" $fn FieldValueType Integer

  set fn MG_MESH_NAME
  setFieldProperties1 "" $fn $stdPropLF {"MG mesh name" 2}
  setFieldProperty "" $fn WidgetWidth 18
  setFieldProperty "" $fn FieldValueType File

  set fn MG_CONVERGENCE_TOLERANCE
  setFieldProperties1 "" $fn $stdPropLF {"MG tolerance" 2}
  setFieldProperty "" $fn Procedurable 0

  set fn MG_MAX_ITERATIONS
  setFieldProperties1 "" $fn $stdPropLF {"MG max iterations" 2}
  setFieldProperty "" $fn FieldValueType Integer

  set fn MG_PRECONDITIONING
  setFieldProperties1 "" $fn $stdPropLF {"MG preconditioning" 2}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 10 }
      { FieldValueType String }
      { Limits {Set {
        "None" 
        "ILU0" "ILU1" "ILU2" "ILU3" "ILU4" "ILU5" "ILU6" "ILU7" "ILU8" "ILU9" "ILUT"} }}
    }

  set fn MG_ILUT_TOLERANCE
  setFieldProperties1 "" $fn $stdPropLF {"MG ILUT Tolerance" 2}
  setFieldProperty "" $fn InitiallyActive 0
  setFieldProperty "" $fn ActivityOnlyByMask 0
  setFieldProperty "" $fn Procedurable 0


  set fn MG_SMOOTHER
  setFieldProperties1 "" $fn $stdPropLF {"MG smoother" 2}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 10 }
      { FieldValueType String }
      { Limits {Set { "Jacobi" "CG" "BiCGStab" } }}
    }
  
  set fn MG_PRE_SMOOTHING_ITERATIONS
  setFieldProperties1 "" $fn $stdPropLF {"MG pre smoothing iters." 2}
  setFieldProperty "" $fn FieldValueType Integer

  set fn MG_POST_SMOOTHING_ITERATIONS
  setFieldProperties1 "" $fn $stdPropLF {"MG post smoothing iters." 2}
  setFieldProperty "" $fn FieldValueType Integer

  set fn MIN_OUTPUT_LEVEL
  setFieldProperties1 "" $fn $stdPropLF {"Min info output level" 1}
  setFieldProperties2 "" $fn \
    { { FieldValueType Integer }
      { Limits {Interval {0 31}} }
      { InitialValue 0 }
    }

  set fn MODEL_DESCRIPTION
  setFieldProperty "" $fn FieldValueType String

  set fn MODEL_DIRECTORY
  setFieldProperty "" $fn FieldValueType File

  set fn MODEL_NAME
  setFieldProperty "" $fn FieldValueType String

  set fn NAVIER-STOKES
  setFieldProperties1 "" $fn $stdPropL {"NAVIER-STOKES EQUATION"}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
    }

  set fn NAVIER-STOKES_vars
  setFieldProperties1 "" $fn $stdPropM {(¤FL1)}
  setFieldProperties2 "" $fn \
    { { WidgetClass OutputField }
      { Display 0 }
      { FieldValueType String }
      { ExcludeFromSif 1 }
    }

  set fn NEWMARK_BETA
  setFieldProperty "" $fn FieldValueType Real
  setFieldProperty "" $fn InitialValue 1.0

  set fn NOF_PROCESSORS
  setFieldProperties1 "" $fn $stdPropL {"Number of processors:"}
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn InitialValue 1

  set fn NON_BLOCKING
  setFieldProperties1 "" $fn $stdPropLUFM {"Non blocking surface" "(Yes/No)" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn NONLINEAR_SYSTEM_CONVERGENCE_TOLERANCE
  setFieldProperty "" $fn Label "Tolerance"
  setFieldProperty "" $fn Procedurable 0

  set fn NONLINEAR_SYSTEM_MAX_ITERATIONS
  setFieldProperty "" $fn Label "Max iterations"
  setFieldProperty "" $fn FieldValueType Integer

  set fn NONLINEAR_SYSTEM_NEWTON_AFTER_ITERATIONS
  setFieldProperty "" $fn Label "Newton after iteration"
  setFieldProperty "" $fn FieldValueType Integer

  set fn NONLINEAR_SYSTEM_NEWTON_AFTER_TOLERANCE
  setFieldProperty "" $fn Label "Newton after tolerance"
  setFieldProperty "" $fn Procedurable 0

  set fn NONLINEAR_SYSTEM_RELAXATION_FACTOR
  setFieldProperty "" $fn Label "Relaxation factor"
  setFieldProperty "" $fn Procedurable 0

  set fn NORMAL_FORCE
  setFieldProperties1 "" $fn $stdPropLUFM {"Normal force" "\[N/m2\]" 2 (¤S) }

  set fn NORMAL-TANGENTIAL_DISPLACEMENT
  setFieldProperties1 "" $fn $stdPropLUM {"Normal tangential displ." "(Yes/No)" (¤S) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn NORMAL-TANGENTIAL_VELOCITY
  setFieldProperties1 "" $fn $stdPropLUM {"Normal tangential veloc." "(Yes/No)" (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn OUTPUT_CALLER
  setFieldProperties1 "" $fn $stdPropLF {"Output calling subroutine info" 1}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue True }
      { AlwaysOutput 1 }
    }

  set fn OUTPUT_FILE
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue "\@Datafile::getDefaultOutputFile"

  set fn OUTPUT_INTERVALS
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn FieldDataType Array

  set fn OUTPUT_LEVEL
  setFieldProperties1 "" $fn $stdPropLF {"Info output level" 1}
  setFieldProperties2 "" $fn \
    { { FieldValueType Integer } 
      { FieldDataType Array }
      { FieldDataSize  "N" }
      { MaxDataSize1  32 }
      { WidgetWidth 20 }
      { InitialValue "" }
      { Limits {Set {0 1}} }
      { Tableable 1 } 
      { Variabled 0 } 
    }

  set fn PCM_NONE
  setFieldProperties1 "" $fn $stdPropLFM {"None" 3 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 7 }
      { Limits {Constant "None"} }
      { ScreenIndent 6 }
	  }

  set fn PCM_SPATIAL_1
  setFieldProperties1 "" $fn $stdPropLFM {"Spatial 1" 3 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 7 }
      { Limits {Constant "Spatial 1"} }
      { ScreenIndent 6 }
	  }

  set fn PCM_SPATIAL_2
  setFieldProperties1 "" $fn $stdPropLFM {"Spatial 2" 3 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 7 }
      { Limits {Constant "Spatial 2"} }
      { ScreenIndent 6 }
	  }

  set fn PCM_TEMPORAL
  setFieldProperties1 "" $fn $stdPropLFM {"Temporal" 3 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 7 }
      { Limits {Constant "Temporal"} }
      { ScreenIndent 6 }
	  }

  set fn PERIODIC_BOUNDARIES
  setFieldProperty "" $fn FieldDataType Array
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn FieldFormat ""

  set fn PHASE_CHANGE_INTERVALS
  setFieldProperties1 "" $fn $stdPropLUM {"Phase change intervals" "\[K\]" (¤HP) }
  setFieldProperties2 "" $fn \
    { { FieldDataType Array }
      { FieldDataSize "2 N" }
	  }

  set fn PHASE_CHANGE_MODEL
  setFieldProperties1 "" $fn $stdPropLFM {"Phase change model:" 3 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { WidgetPacking Vertical }
      { FieldValueType String }
      { InitialValue "None" }
      { ScreenIndent 0 }
	  }

  set fn PLANE_STRESS
  setFieldProperties1 "" $fn $stdPropLFM {"Plane stress" 2 {(¤S) (¤dimS2)&(¤S)} }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
      { ScreenIndent 0 }
	  }

  set fn POISSON_RATIO
  setFieldProperties1 "" $fn $stdPropLM {"Poisson ratio" (¤S) }
  setFieldProperty    "" $fn  NeededLevel 1

  set fn POST_FILE
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue "\@Datafile::getDefaultPostFile"

  set fn PRESSURE
  setFieldProperties1 "" $fn $stdPropLUM {"Pressure" "\[Pa\]" (¤FL) }	
  setFieldProperty "" $fn SifName "Pressure"

  set fn PRESSURE_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUFM {"Pressure" "\[N/m2\]" 2 (¤FL) }
    setFieldProperty    "" $fn$i SifName "Pressure $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" PRESSURE_2 Mask {(¤FL) (¤FL)&((¤dimS2)|(¤dimS3))}
  setFieldProperty    "" PRESSURE_3 Mask {(¤FL) (¤FL)&(¤dimS3)}

  set fn PROBLEM_DESCRIPTION
  setFieldProperty "" $fn FieldValueType String

  set fn PROBLEM_NAME
  setFieldProperty "" $fn FieldValueType String

  set fn PROCEDURE
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType File }
      { FieldDataSize 2 }
    }

  set fn PROCEDURE_BUTTON
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType CheckBox }
      { WidgetPacking Horizontal }
      { WidgetLabel "Proc" }
      { WidgetCommandProc Panel::procedureCheckBoxCommandProc }
      { WidgetWidth 0 }
      { Limits {Set {1 0}} }
      { InitialValue 0 }
      { InitiallyActive 0 }
      { ActivityOnlyByMask 0}
    }

  set fn RADIATION
  setFieldProperties1 "" $fn $stdPropLFM {"Radiation" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { FieldValueType String }
      { InitialValue "" }
    }

  # Not in use!
  set fn RADIATION_BOUNDARIES
  setFieldProperties1 "" $fn $stdPropFM {2 (¤H) }
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Integer }
      { ActivityOnlyByMask 0 }
      { InitialValue "" }
    }

  set fn RADIATION_TARGET_BODY
  setFieldProperties1 "" $fn $stdPropFM {2 (¤H) }
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Integer }
      { ActivityOnlyByMask 0 }
      { ActivityParents RADIATION }
      { ActivityValues {"Diffuse gray"} }
      { InitialValue "" }
    }

  # Experimental, no in use!
  set fn RADIATION_TARGET_BODIES
  setFieldProperties1 "" $fn $stdPropFM {2 (¤H) }
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Integer }
      { ActivityOnlyByMask 0 }
      { InitialValue "" }
    }

  set fn RADIATION_TARGET_BODY_NAME
  setFieldProperties1 "" $fn $stdPropLFM {"Radiation target body" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField }
      { WidgetType OptionMenu }
      { WidgetWidth 9 }
      { FieldValueType String }
      { Limits {Set {""} } }
      { InitialValue "" }
      { ActivityOnlyByMask 0 }
    }

  set fn RD_IDEALIZED
  setFieldProperties1 "" $fn $stdPropLFM {"Idealized" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 7 }
      { Limits {Constant "Idealized"} }
    }

  set fn RD_NONE
  setFieldProperties1 "" $fn $stdPropLFM {"None" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 4 }
      { Limits {Constant ""} }
    }

  set fn RD_DIFFUSE_GRAY
  setFieldProperties1 "" $fn $stdPropLFM {"Diffuse gray" 2 (¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetPacking Horizontal }
      { WidgetWidth 9 }
      { Limits {Constant "Diffuse gray"} }
    }

  set fn REFERENCE_PRESSURE
  setFieldProperties1 "" $fn $stdPropLUFM {"Reference pressure" "\[Pa\]" 2 (¤FL)|(¤HT2)|(¤HT3)|(¤AD2)|(¤AD3) }
  setFieldProperties2 "" $fn \
    { { InitialValue 0.0 }
      { NeededLevel 2 }
      { ActivityOnlyByMask 0 }
      { InitiallyActive 0 }
      { ActivityParents COMPRESSIBILITY_MODEL }
      { ActivityValues {"~Incompressible"} }
    }


  set fn REFERENCE_TEMPERATURE
  setFieldProperties1 "" $fn $stdPropLUM {"Reference temperature" "\[K\]" (¤FL)|(¤HT2)|(¤HT3)|(¤AD2)|(¤AD3)|(¤S) }

  set fn RELOAD_INPUT_FILE
  setFieldProperties1 "" $fn $stdPropL "ReloadInput file"
  setFieldProperties2 "" $fn \
    { { WidgetType  BrowsableFile }
      { HasAbsCheckBox 1 }
      { HasUseCheckBox 1 } 
      { FieldValueType File }
      { FieldDataSize 1 }
      { WidgetWidth 25 }
    }

  set fn RESTART_FILE
  setFieldProperty "" $fn FieldValueType File

  set fn RESTART_POSITION
  setFieldProperties2 "" $fn { 
      { FieldValueType Integer }
      { FieldFormat "%10d" }
      { Procedurable 0 }
      { Tableable 0 } 
    }

  set fn RESULT_MESH_NAME
  setFieldProperties1 "" $fn $stdPropLF {"Result mesh name" 0}
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue ""

  set fn RESULTS_DIRECTORY
  setFieldProperties1 "" $fn $stdPropL {"Results directory"}
  setFieldProperties2 "" $fn \
    { { WidgetWidth 40 }
      { FieldValueType File }
      { FieldFormat "% -39s" }
    }

  set fn SELECT_MESH_WIDGET
  setFieldProperties1 "" $fn $stdPropLF {"Select mesh" 0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField }
      { WidgetType OptionMenu }
      { WidgetWidth 9 }
      { FieldValueType String }
      { Limits {Set {$Model(meshNames)} } }
      { InitialValue "\$Model(currentMeshName)" }
    }

  set fn SIF_MODEL_INCLUDE_FILE
  setFieldProperties1 "" $fn $stdPropL {"General include file"}
  setFieldProperties2 "" $fn \
    { { WidgetType IncludeFile }
      { HasUseCheckBox 1 }
      { WidgetWidth 40 }
      { FieldValueType File }
      { FieldFormat "% -39s" }
    }

  set fn SIF_SIMULATION_INCLUDE_FILE
  setFieldProperties1 "" $fn $stdPropL {"Simulation include file"}
  setFieldProperties2 "" $fn \
    { { WidgetType IncludeFile }
      { HasUseCheckBox 1 }
      { WidgetWidth 40 }
      { FieldValueType File }
      { FieldFormat "% -39s" }
    }

  set fn SIMULATION_TYPE
  setFieldProperty "" $fn FieldValueType String
  setFieldProperty "" $fn InitialValue "Steady State"

  set fn SM_MECHANICAL
  setFieldProperties1 "" $fn $stdPropLFM {"Mechanical" 2 (¤S) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 9 }
      { Limits {Constant "Mechanical"} }
      { ScreenIndent 6 }
    }

  set fn SM_THERMAL
  setFieldProperties1 "" $fn $stdPropLFM {"Thermal" 2 (¤S)&(¤H) }
  setFieldProperties2 "" $fn \
    { { WidgetClass MemberOutputField }
      { WidgetType RadioButton }
      { WidgetWidth 9 }
      { Limits {Constant "Thermal"} }
      { ScreenIndent 6 }
	  }

  set fn SOLVER_INPUT_FILE
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue "\@Datafile::getDefaultSolverInputFile"

  set fn SOLVING_ORDER
  setFieldProperty "" $fn FieldValueType Integer

  set fn SORET_DIFFUSIVITY
  setFieldProperties1 "" $fn $stdPropLUM {" Soret Diffusivity" "\[kg/m\]" (¤AD) }
  setFieldProperty    "" $fn IndexType "Pre"
  setFieldProperty    "" $fn IndexVariable ADVECTION_DIFFUSION_EQUATION
  setFieldProperty    "" $fn WidgetClass ParentField

  set fn SPECIFIC_HEAT_RATIO
  setFieldProperties1 "" $fn $stdPropLFM {"Specific Heat Ratio" 2 (¤FL)|(¤HT2)|(¤HT3)|(¤AD2)|(¤AD3) }
  setFieldProperties2 "" $fn \
    { { InitialValue 1.4 }
      { NeededLevel 2 }
      { ActivityOnlyByMask 0 }
      { InitiallyActive 0 }
      { ActivityParents COMPRESSIBILITY_MODEL }
      { ActivityValues {"~Incompressible"} }
    }

  set fn STABILIZE
  setFieldProperties1 "" $fn $stdPropL {"Stabilize"}
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn STEADY_STATE_CONVERGENCE_TOLERANCE
  setFieldProperties1 "" $fn $stdPropLF {"Steady state tolerance" 2}
  setFieldProperty "" $fn Procedurable 0

  set fn STEADY_STATE_MAX_ITERATIONS
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn InitialValue 20

  set fn STEADY_STATE_OUTPUT_INTERVAL
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn InitialValue 1
  setFieldProperty "" $fn SifName "Output Intervals"

  set fn STEFAN_BOLTZMANN
  setFieldProperties1 "" $fn $stdPropLUF { "Stefan-Boltzmann constant" "\[W/m2K4\]" 2}
  setFieldProperty "" $fn InitialValue 5.67e-08
  setFieldProperty "" $fn Procedurable 0

  set fn STRESS_ANALYSIS
  setFieldProperties1 "" $fn $stdPropL {"STRESS ANALYSIS"}
  setFieldProperties2 "" $fn \
    { { Display 0 }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
      { InitialValue False }
	  }

  set fn STRESS_ANALYSIS_vars
  setFieldProperties1 "" $fn $stdPropM {(¤S)}
  setFieldProperties2 "" $fn \
    { { WidgetClass OutputField }
      { Display 0 }
      { FieldValueType String }
      { ExcludeFromSif 1 }
    }

  set fn STRESS_BODYFORCE_
  foreach i {1 2 3} {
    setFieldProperties1 "" $fn$i $stdPropLUM {"Stress Bodyforce" "\[N/m3]" (¤S) }
    setFieldProperty    "" $fn$i SifName "Stress Bodyforce $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" STRESS_BODYFORCE_2 Mask {(¤S) (¤S)&((¤dimS2)|(¤dimS3))}
  setFieldProperty    "" STRESS_BODYFORCE_3 Mask {(¤S) (¤S)&(¤dimS3)}


  set fn STRESS_MODEL
  setFieldProperties1 "" $fn $stdPropFM {2 (¤S) }
  setFieldProperties2 "" $fn \
    { { WidgetClass RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { WidgetPacking Vertical }
      { FieldValueType String }
	  }

  set fn SURFACE_ROUGHNESS
  setFieldProperties1 "" $fn $stdPropLM  {"Surface roughness" (¤FT1) }
  setFieldProperty "" $fn InitialValue 9

  set fn SURFACE_TENSION_COEFFICIENT
  setFieldProperties1 "" $fn $stdPropLUFM {"Surf. tension coeff." "\[N/m\]" 2 (¤FL) }

  set fn SURFACE_TENSION_EXPANSION_COEFFICIENT
  setFieldProperties1 "" $fn $stdPropLUFM {"Surf. tension exp. coeff" "\[1/K\]" 2 (¤FL) }

  set fn TABLE_BUTTON
  setFieldProperties1 "" $fn $stdPropF {0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType CheckBox }
      { WidgetPacking Horizontal }
      { WidgetLabel "Table" }
      { WidgetCommandProc Panel::tableCheckBoxCommandProc }
      { WidgetWidth 0 }
      { Limits {Set {1 0}} }
      { InitialValue 0 }
      { InitiallyActive 0 }
      { ActivityOnlyByMask 0}
    }

  set fn TEMPERATURE
  setFieldProperties1 "" $fn $stdPropLUM {"Temperature" "\[K\]" (¤H) }
  setFieldProperty "" $fn SifName "Temperature"
  setFieldProperty "" $fn SifTypeGiven 1

  # NOTE: We do not set (¤T) mask, because time-variable
  # used as step counter for Steady-state problems!!!
  set fn TIME
  setFieldProperties1 "" $fn $stdPropLUF {"Time" "\[s\]" 0 }
  setFieldProperty "" $fn SifName "Time"

  set fn TIME_DERIVATIVE_ORDER
  setFieldProperty "" $fn FieldValueType Integer

  set fn TIMESTEPPING_METHOD
  setFieldProperties1 "" $fn $stdPropL {"Time stepping method:"}
  setFieldProperties2 "" $fn \
    { { WidgetType OptionMenu }
      { WidgetWidth 18 }
      { FieldValueType String }
      { Limits {Set {"Newmark" "BDF" } } }
      { InitialValue "Newmark" }
    }

  set fn TIMESTEP_INTERVALS
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn FieldDataType Array
  setFieldProperty "" $fn Tableable 0

  set fn TIMESTEP_SIZES
  setFieldProperty "" $fn FieldDataType Array
  setFieldProperty "" $fn Tableable 0

  set fn TM_NONE
  setFieldProperties1 "" $fn $stdPropLFM {"None" 2 (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType RadioButton }
      { WidgetWidth 7 }
      { Limits {Constant ""} }
      { ScreenIndent 6 }
	  }

  set fn TM_KE
  setFieldProperties1 "" $fn $stdPropLFM {"KE" 2 (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetType RadioButton }
      { WidgetWidth 7 }
      { Limits {Constant "KE"} }
      { ScreenIndent 6 }
	  }

  set fn TURBULENCE_MODEL
  setFieldProperties1 "" $fn $stdPropLFM {"Turbulence model:" 2 (¤FL) }
  setFieldProperties2 "" $fn \
    { { WidgetClass  RadioButtonsContainer }
      { WidgetType GroupDataAndScreen }
      { WidgetPacking Vertical }
      { InitialValue "" }
      { ExcludeFromSif 1}
	  }

  set fn VARIABLE
  setFieldProperty "" $fn FieldDataSize "N N"
  setFieldProperty "" $fn FieldValueType String

  set fn VARIABLE_DOFS
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn Display 0

# NOTE: This is used for Calculators for storing
# dofs-array entered by the user. The actual (possibly dimension
# dependent dof-value is stored in VARIABLE_DOFS)
 set fn VARIABLE_DOFS_ALL
  setFieldProperty "" $fn FieldValueType Integer
  setFieldProperty "" $fn ExcludeFromSif 1


  set fn VARIABLE_BUTTON
  setFieldProperties1 "" $fn $stdPropLF {"   " 0}
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType Button }
      { WidgetCommandProc "EquationVariable::openPanel" }
      { WidgetLabel "Variables" }
      { WidgetWidth 0 }
      { WidgetPacking Horizontal }
      { InitiallyActive 0 }
      { ActivityOnlyByMask 0 }
    }

  set fn VARIABLE_LIST
  setFieldProperties2 "" $fn \
    { { WidgetClass PanelField } 
      { WidgetType SelectionBox }
      { WidgetLabel "Variables:" }
      { WidgetHeight 6 }
      { WidgetWidth 16 }
      { WidgetPacking Vertical }
      { SubPanel 3 }
      { Limits {Set {@Equation::updateVariableList}} }
      { ActivityOnlyByMask 0 }
      { Mask {(¤AD) (¤AD2)} }

     }

  set fn VELOCITY
  setFieldProperties1 "" $fn $stdPropLUM {"Velocity" "\[m/s\]" (¤FL)}
	setFieldProperty "" $fn FieldDataSize "DIM"
  setFieldProperty "" $fn FieldDataType Array

  set fn VELOCITY_
  foreach i { 1 2 3 } {
    setFieldProperties1 "" $fn$i $stdPropLUM {"Velocity" "\[m/s\]" (¤FL) }
    setFieldProperty    "" $fn$i SifName "Velocity $i"
    setFieldProperty    "" $fn$i CoordinateIndex $i
    setFieldProperty    "" $fn$i LabelByCoordinate 1
  }
  setFieldProperty    "" VELOCITY_2 Mask {(¤FL) (¤FL)&((¤dimS2)|(¤dimS3))}
  setFieldProperty    "" VELOCITY_3 Mask {(¤FL) (¤FL)&(¤dimS3)}


  set fn VIEW_FACTORS
  setFieldProperty "" $fn FieldValueType File
  setFieldProperty "" $fn InitialValue "Viewfactors.dat"

  set fn VISCOSITY
  setFieldProperties1 "" $fn $stdPropLUM {"Viscosity" "\[kg/ms\]" (¤FL) }
  setFieldProperty    "" $fn  NeededLevel 1

  set fn WALL_LAW
  setFieldProperties1 "" $fn $stdPropLUM {"Wall law" "(Yes/No)" (¤FT) }
  setFieldProperties2 "" $fn \
    { { WidgetType CheckBox }
      { FieldValueType Logical }
      { Limits {Set {True False}} }
    }

  set fn YOUNGS_MODULUS
  setFieldProperties1 "" $fn $stdPropLM {"Young's modulus" (¤S) }
  setFieldProperty    "" $fn  NeededLevel 1

  namespace forget DataField::getFieldPropert* DataField::setFieldPropert*

}
### End FRONT_SET_INITIAL_GENERIC_FIELD_PROPERITES proc






#===================================================#
#                                                   #
#     FRONT_SET_INITIAL_PANEL_FIELD_PROPERTIES      #
#                                                   #
#===================================================#
#
# Panel specific field property definitions 
#
proc FRONT_SET_INITIAL_PANEL_FIELD_PROPERTIES {} {
  global Common Info Model
  
  foreach arr $Common(allPanelArries) {
    global $arr
  }

  # These specific propertiy values are needed for nearly all fields
  #
  set stdProp2    { Label UnitString SubPanel Mask }
  set stdPropLUFM { Label UnitString SubPanel Mask }
  set stdPropLUF  { Label UnitString SubPanel }
  set stdPropLUM  { Label UnitString Mask }
  set stdPropLFM  { Label SubPanel Mask }
  set stdPropLU   { Label UnitString }
  set stdPropLF   { Label SubPanel }
  set stdPropLM   { Label Mask }
  set stdPropL    { Label }
  set stdPropUFM  { UnitString SubPanel Mask }
  set stdPropUM   { UnitString Mask }
  set stdPropFM   { SubPanel Mask }
  set stdPropF    { SubPanel }
  set stdPropM    { Mask }

  namespace import DataField::getFieldPropert* DataField::setFieldPropert*


#============#
# BodyForce  #
#============#

  #-Panel variables
  set BodyForce(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
    BOUSSINESQ
    FLOW_BODYFORCE_1 FLOW_BODYFORCE_2 FLOW_BODYFORCE_3
    STRESS_BODYFORCE_1 STRESS_BODYFORCE_2 STRESS_BODYFORCE_3
    HEAT_SOURCE
    DIFFUSION_SOURCE
  }


  set bf BodyForce
  set BF $BodyForce(parameterType)
  set fn INCLUDE
  set Common($BF,$fn,editProc) exists

  setFieldProperties1 $bf $fn $stdPropLF {"Parameter file" 0}


#===============#
# BodyParameter #
#===============#

  set BodyParameter(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
  }


#==========#
# Boundary #
#==========#
  set Boundary(initialFields) ""


#====================#
# Boundary Condition #
#====================#

  #-Panel variables
  set BoundaryCondition(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
    VELOCITY_1 VELOCITY_2 VELOCITY_3
    NORMAL-TANGENTIAL_VELOCITY
    DISPLACEMENT_1 DISPLACEMENT_2 DISPLACEMENT_3
    NORMAL-TANGENTIAL_DISPLACEMENT
    FORCE_1 FORCE_2 FORCE_3
    NORMAL_FORCE
    PRESSURE_1 PRESSURE_2 PRESSURE_3
    PRESSURE
    EXTERNAL_PRESSURE
    FREE_MOVING
    FREE_SURFACE
    KINETIC_ENERGY
    KINETIC_DISSIPATION 
    SURFACE_ROUGHNESS
    WALL_LAW
    BOUNDARY_LAYER_THICKNESS
    SURFACE_TENSION_COEFFICIENT
    SURFACE_TENSION_EXPANSION_COEFFICIENT
    TEMPERATURE
    HEAT_FLUX
    HEAT_TRANSFER_COEFFICIENT
    EXTERNAL_TEMPERATURE
    RADIATION
    EMISSIVITY
    RADIATION_TARGET_BODY
    RADIATION_TARGET_BODY_NAME
    FLOW_FORCE_BC
    HEAT_FLUX_BC
    DIFFUSION
    DIFFUSION_FLUX
  }


  set bc BoundaryCondition
  set BC $BoundaryCondition(parameterType)

  set fn INCLUDE
  set Common($BC,$fn,editProc) exists
  setFieldProperties1 $bc $fn $stdPropLF {"Parameter file" 0}

  set fn RADIATION
  set Common($BC,$fn,fillProc) exits
  set Common($BC,$fn,group) {RD_NONE RD_IDEALIZED RD_DIFFUSE_GRAY}


  set BoundaryCondition(initialDirichletVnames) {
    VELOCITY_1
    VELOCITY_2
    VELOCITY_3
    NORMAL-TANGENTIAL_VELOCITY
    PRESSURE
    TEMPERATURE
    DISPLACEMENT_1
    DISPLACEMENT_2
    DISPLACEMENT_3
    NORMAL-TANGENTIAL_DISPLACEMENT
    DIFFUSION
  }

  set BoundaryCondition(initialNeumannVnames) {
    FORCE_1 FORCE_2 FORCE_3
    NORMAL_FORCE
    PRESSURE_1 PRESSURE_2 PRESSURE_3
    EXTERNAL_PRESSURE
    FREE_MOVING
    FREE_SURFACE
    KINETIC_ENERGY
    KINETIC_DISSIPATION 
    SURFACE_ROUGHNESS
    WALL_LAW
    BOUNDARY_LAYER_THICKNESS
    SURFACE_TENSION_COEFFICIENT
    SURFACE_TENSION_EXPANSION_COEFFICIENT
    HEAT_FLUX
    HEAT_TRANSFER_COEFFICIENT
    EXTERNAL_TEMPERATURE
    RADIATION
    RD_NONE
    RD_IDEALIZED
    RD_DIFFUSE_GRAY
    EMISSIVITY
    RADIATION_TARGET_BODY
    RADIATION_TARGET_BODY_NAME
    FLOW_FORCE_BC
    HEAT_FLUX_BC
    DIFFUSION_FLUX
  }

  set BoundaryCondition(subElementFields) $BoundaryCondition(initialDirichletVnames)



#===================#
# BoundaryParameter #
#===================#

  set BoundaryParameter(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
  }


#============#
# Calculator #
#============#

  # Panel variables
  set Calculator(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    ACTIVE
    EXEC_SOLVER
    EQUATION
    LIBRARY_FILE
    FUNCTION_NAME
    PROCEDURE
    VARIABLE
    VARIABLE_DOFS
    VARIABLE_DOFS_ALL
    SOLVING_ORDER
  }


  set cl Calculator
  set CL $Calculator(parameterType)

  set lw 14
  set ww 32

  set fn EQUATION
  setFieldProperty $cl $fn WidgetType Entry
  setFieldProperty $cl $fn FieldType String
  setFieldProperty $cl $fn SubPanel 1
  setFieldProperty $cl $fn Label "Equation"
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn WidgetWidth 40
  setFieldProperty $cl $fn WidgetAnchor w
  setFieldProperty $cl $fn NeededLevel 0

  set fn ACTIVE
  setFieldProperty $cl $fn WidgetType CheckBox
  setFieldProperty $cl $fn SubPanel 1
  setFieldProperty $cl $fn Label "Active"
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn WidgetAnchor w

  set fn LIBRARY_FILE
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn WidgetAnchor w
  setFieldProperty $cl $fn WidgetWidth $ww
  setFieldProperty $cl $fn NeededLevel 2

  set fn FUNCTION_NAME
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn WidgetAnchor w
  setFieldProperty $cl $fn WidgetWidth $ww
  setFieldProperty $cl $fn NeededLevel 2

  set fn VARIABLE
  setFieldProperty $cl $fn Label "Variable name"
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn SubPanel 1
  setFieldProperty $cl $fn WidgetType Entry
  setFieldProperty $cl $fn WidgetWidth $ww
  setFieldProperty $cl $fn WidgetAnchor w
  setFieldProperty $cl $fn FieldType String
  setFieldProperty $cl $fn InitialValue ""
  setFieldProperty $cl $fn NeededLevel 0

  set fn VARIABLE_DOFS_ALL
  setFieldProperty $cl $fn Label "Variable DOFs"
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn SubPanel 1
  setFieldProperty $cl $fn WidgetType Entry
  setFieldProperty $cl $fn WidgetWidth 10
  setFieldProperty $cl $fn WidgetAnchor w
  setFieldProperty $cl $fn Display 1
  setFieldProperty $cl $fn FieldDataType Array
  setFieldProperty $cl $fn FieldDataSize  "N"
  setFieldProperty $cl $fn MaxDataSize1  3
  setFieldProperty $cl $fn NeededLevel 0
  setFieldProperty $cl $fn Tableable 0

  set fn EXEC_SOLVER
  setFieldProperties1 $cl $fn $stdPropLF {"Exec calculator" 1}
  setFieldProperty $cl $fn InitialValue { "Always" }
  setFieldProperty $cl $fn LabelWidth $lw
  setFieldProperty $cl $fn WidgetAnchor w

  set fn SOLVING_ORDER
  setFieldProperty $cl $fn Display 0
  setFieldProperty $cl $fn ExcludeFromSif 1



#==========#
# Constant #
#==========#

  set cn Constant
  set CN $Constant(parameterType)

  set fn GRAVITY_HEADER
  setFieldProperties1 $cn $fn $stdPropL {"Gravity vector:"}
  setFieldProperty $cn $fn WidgetType Label

  set fn PHYSICAL_CONSTANTS_HEADER
  setFieldProperties1 $cn $fn $stdPropLF {"Physical constants:" 2}
  setFieldProperty $cn $fn WidgetType Label

  set Constant(initialFields) { 
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    GRAVITY_HEADER
    GRAVITY_1
    GRAVITY_2
    GRAVITY_3
    GRAVITY_ABS
    GRAVITY
    STEFAN_BOLTZMANN
  }


#============#
# Coordinate #
#============#

  set Coordinate(initialFields) {
    COORDINATE_SYSTEM
    COORDINATE_MAPPING
  }

  set Coordinate(coordinateLabelIndex) 0

  set Coordinate(coordinateLabels) {
    { X Y Z }
    { R Z P }
    { R T P }
  }


#==========#
# Datafile #
#==========#

  set Datafile(initialFields) {
    SOLVER_INPUT_FILE
    MESH_INPUT_FILE
    OUTPUT_FILE
    RESTART_FILE
    RESTART_POSITION
    POST_FILE
    GEBHARDT_FACTORS
    VIEW_FACTORS
  }


#==========#
# Equation #
#==========#

  # Predefined equation indices
  # 0 = Navier-Stokes
  # 1 = KE Turbulence
  # 2 = Heat equation
  # 3 = Stress Analysis
  # 4 = Advection Diffusion Equation
  set Equation(predefinedIndices) {0 1 2 3 4}
  set Equation(indices) {0 1 2 3 4}
  set Equation(equationIndices) ""

  #-Panel variables
  # NOTE: STRESS_MODEL /Mechanical/Thermal)removed 01.10.99 /Mve
  # There is no support for these keywords in Solver!
  set Equation(STRESS_MODEL) ""

  #
  set Equation(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
    VARIABLE_BUTTON
    EQUATION_LIST
    VARIABLE_LIST
    NAVIER-STOKES 
      NAVIER-STOKES_vars 
      HYDROSTATIC_PRESSURE 
      TURBULENCE_MODEL
        KE_TURBULENCE
          KE_TURBULENCE_vars
          KE_CLIP
    HEAT_EQUATION
      HEAT_EQUATION_vars
      CONVECTION
      PHASE_CHANGE_MODEL
      CHECK_LATENT_HEAT_RELEASE
    STRESS_ANALYSIS
      STRESS_ANALYSIS_vars
      PLANE_STRESS
    ADVECTION_DIFFUSION_EQUATION
      ADVECTION_DIFFUSION_EQUATION_vars
  }


  set eq Equation
  set EQ $Equation(parameterType)

  #--------------
  # NAVIER-STOKES
  #--------------
  setFieldProperties2 Equation NAVIER-STOKES \
    { { EquationIndex 0 }
      { EquationMask (¤FL1) }
      { EquationLabel "Navier-Stokes" }
      { EquationField NAVIER-STOKES }
      { EquationVars NAVIER-STOKES_vars }
      { IsMultiVar 0 }
      { Variable "Flow Solution" }
      { VariableComponents {Velocity; Pressure} }
      { VariableDofs { 1 2 3; 1 } }
      { IncludeInEquationPanel 1 }
      { SystemInfoIndices 0 }
      { SifName "Navier-Stokes" }
      { StatusLineName "Flow" }
      { ProcessTableName NS }
	  }


  #--------------
  # KE TURBULENCE
  #--------------
  setFieldProperties2 Equation KE_TURBULENCE \
    { { EquationIndex 1 }
      { EquationMask (¤FT1) }
      { EquationLabel "KE turbulence" }
      { EquationField KE_TURBULENCE}
      { EquationVars KE_TURBULENCE_vars}
      { IsMultiVar 0 }
      { Variable "KE Solution" }
      { VariableComponents {Kinetic Energy; Kinetic Dissipation} }
      { VariableDofs { 1; 1 } }
      { IncludeInEquationPanel 0 }
      { SystemInfoIndices {0 1} }
      { SifName "KE Turbulence" }
      { StatusLineName "KE-Turb" }
      { ProcessTableName {KE KE-D} }
	  }


  # NOTE: We store this variable in the model data, but it is
  # not in the Solver-input file
  # This variable is needed to get value for the variables
  # like KE_TURBULENCE 
  #
  set fn TURBULENCE_MODEL
  set Common($EQ,$fn,group) {TM_NONE TM_KE}
  set Common($EQ,$fn,groupEditProc) exists


  #--------------
  # HEAT EQUATION
  #--------------
  setFieldProperties2 Equation HEAT_EQUATION \
    { { EquationIndex 2 }
      { EquationMask (¤H) }
      { EquationLabel "Heat equation" }
      { EquationField HEAT_EQUATION}
      { EquationVars HEAT_EQUATION_vars}
      { IsMultiVar 0 }
      { Variable "Temperature" }
      { VariableComponents "" }
      { VariableDofs 1 }
      { IncludeInEquationPanel 1 }
      { SystemInfoIndices 0 }
      { SifName "Heat Equation" }
      { StatusLineName "Heat" }
      { ProcessTableName DC }
	  }

  set fn CONVECTION
  set Common($EQ,$fn,group) {CONV_NONE CONV_CONSTANT CONV_COMPUTED}
  set Common($EQ,$fn,groupEditProc) exists

  set fn PHASE_CHANGE_MODEL
  set Common($EQ,$fn,group) {PCM_NONE PCM_SPATIAL_1 PCM_SPATIAL_2 PCM_TEMPORAL}
  set Common($EQ,$fn,groupEditProc) exists


  #----------------
  # STRESS ANALYSIS
  #----------------
  setFieldProperties2 Equation STRESS_ANALYSIS \
    { { EquationIndex 3 }
      { EquationMask (¤S) }
      { EquationLabel "Stress analysis" }
      { EquationField STRESS_ANALYSIS}
      { EquationVars STRESS_ANALYSIS_vars}
      { IsMultiVar 0 }
      { Variable "Displacement" }
      { VariableComponents "" }
      { VariableDofs {1 2 3} }
      { IncludeInEquationPanel 1 }
      { SystemInfoIndices 0 }
      { SifName "Stress Analysis" }
      { StatusLineName "Stress" }
      { ProcessTableName SA }
	  }

  # NOTE: not in use currently
  # MVe 09/99
  set fn STRESS_MODEL
  set Common($EQ,$fn,group)  {SM_MECHANICAL SM_THERMAL}
  set Common($EQ,$fn,groupEditProc) exists


  #-----------------------------
  # ADVECTION DIFFUSION EQUATION
  #-----------------------------
  setFieldProperties2 Equation ADVECTION_DIFFUSION_EQUATION \
    { { EquationIndex 4 }
      { EquationMask (¤AD) }
      { EquationLabel "Advection diffusion" }
      { EquationField ADVECTION_DIFFUSION_EQUATION }
      { EquationVars ADVECTION_DIFFUSION_EQUATION_vars }
      { IsMultiVar 1 }
      { Variable "" }
      { VariableComponents "" }
      { VariableDofs 1 }
      { IncludeInEquationPanel 1 }
      { SystemInfoIndices 0 }
      { SifName "Advection Diffusion Equation" }
      { StatusLineName "AdvDiff" }
      { ProcessTableName AD }
	  }

  set fn ADVECTION_DIFFUSION_EQUATION
  set Common($EQ,$fn,fillProc) exists



#===================#
# EquationVariable  #
#===================#

# NOTE: AlwaysAddInitialValue is needed to prevent a "real" initial
# like Species for AdvDiff to be added always to the list of varaibles
# In conratry, user defined varaibles are always added, otherwise they
# would not been seen at all, if the user would add a differently named
# variable compared to waht is stored in model file.
# Ref. proc EquationVariable::updateFieldData
#
  set EquationVariable(initialFields) {
    NAVIER-STOKES
    KE_TURBULENCE
    HEAT_EQUATION
    STRESS_ANALYSIS
    ADVECTION_DIFFUSION_EQUATION
  }

  set ev EquationVariable

  set fn NAVIER-STOKES
  setFieldProperties2 $ev $fn \
    { { Display 0 }
      { FieldValueType String }
      { InitialValue "" }
      { AlwaysAddInitialValue 0 }
	  }

  set fn KE_TURBULENCE
  setFieldProperties2 $ev $fn \
    { { Display 0 }
      { FieldValueType String }
      { InitialValue "" }
      { AlwaysAddInitialValue 0 }
	  }

  set fn HEAT_EQUATION
  setFieldProperties2 $ev $fn \
    { { Display 0 }
      { FieldValueType String }
      { InitialValue "" }
      { AlwaysAddInitialValue 0 }
	  }

  set fn STRESS_ANALYSIS
  setFieldProperties2 $ev $fn \
    { { Display 0 }
      { FieldValueType String }
      { InitialValue "" }
      { AlwaysAddInitialValue 0 }
	  }
  
  # NOTE: Only AdvDiff has here a value!!!
  #
  set fn ADVECTION_DIFFUSION_EQUATION
  setFieldProperties2 $ev $fn \
    { { Display 0 }
      { FieldValueType String }
      { InitialValue "Species" }
      { AlwaysAddInitialValue 0 }
	  }


#=======#
# GridH #
#=======#

  set gh GridH
  set GH $GridH(parameterType)

  set fn MESH_DENSITY_TYPE
  set Common($GH,$fn,group) {MESH_DT_NONE MESH_DT_H MESH_DT_R MESH_DT_N} 

  #-Panel variables
  set GridH(initialFields) {
    MESH_NAME
    MESH_DENSITY_TYPE
    MESH_H
    MESH_R
    MESH_N
  }


#===============#
# GridParameter #
#===============#

  #-Panel variables
  set GridParameter(initialFields) {
    ACTIVE_IN_MESHING
    MESH_NAME
    MESH_ELEMENT_TYPE
    MESH_ELEMENT_ORDER
    MESH_LAYER_TYPE
    MESH_BG_MESH
    MESH_BG_MESH_FILE
    MESH_DENSITY_TYPE
    MESH_H
    MESH_R
    MESH_QUADGRID_N1
    MESH_QUADGRID_N2
    LISTBOX1
    MESH_SEED_EDGE
    MESH_SEED_TYPE
  }

  # End of list will be like this
  # when vertex/direction is implemented
  #MESH_QUADGRID_N2
  #MESH_SEED_TYPE
  #LISTBOX1
  #MESH_SEED_EDGE
  #MESH_SEED_VERTEX
  #MESH_SEED_DIRECTION

  set GridParameter(lineGroup) {
    MESH_H MESH_R
    MESH_DT_NONE MESH_DT_H MESH_DT_R
  }

  set GridParameter(triangleGroup) {
    MESH_LAYER_TYPE MESH_BG_MESH
    MESH_H MESH_R
    MESH_DT_NONE MESH_DT_H MESH_DT_R
  }

  set GridParameter(seededTriangleGroup)  { 
    MESH_SEED_TYPE MESH_SEED_EDGE 
  }

  # Not in use. This is the "future" group!
  set GridParameter(seededTriangleGroup2)  { 
    MESH_SEED_EDGE MESH_SEED_VERTEX MESH_SEED_DIRECTION
    MESH_SEED_TYPE MESH_ST_IMPLICIT MESH_ST_EXPLICIT
  }

  set GridParameter(quadGridGroup) {
    MESH_QUADGRID_N1 MESH_QUADGRID_N2
  }

  # Not in use!
  set GridParameter(meshMethodIndices,LINE) 0
  set GridParameter(meshMethodIndices,TRIANGLE) 1
  set GridParameter(meshMethodIndices,TRIANGLE_GRID) 2
  set GridParameter(meshMethodIndices,QUAD_GRID) 3
  set GridParameter(meshMethodIndices,TETRA) 4
  set GridParameter(meshMethodIndices,BRICK) 5

  set gr GridParameter
  set GR $GridParameter(parameterType)

  set fn MESH_ELEMENT_TYPE
  set Common($GR,$fn,group) {MESH_ET_LINE MESH_ET_TRIANGLE MESH_ET_QUAD MESH_ET_TETRA MESH_ET_BRICK}

  set fn MESH_DENSITY_TYPE
  set Common($GR,$fn,group) {MESH_DT_NONE MESH_DT_H MESH_DT_R} 

  # NOTE!!!
  # This is currently only output data, value is always the default
  # value, ie. "Implicit"
  set fn MESH_SEED_TYPE
  set Common($GR,$fn,group) {MESH_ST_IMPLICIT MESH_ST_EXPLICIT} 
  setFieldProperty $gr $fn WidgetClass OutputField
  setFieldProperty $gr $fn Display 0

  set fn MESH_SEED_VERTEX
  setFieldProperty $gr $fn Display 0

  set fn MESH_SEED_DIRECTION
  setFieldProperty $gr $fn Display 0

  set fn LISTBOX1
  setFieldProperties1 $gr $fn $stdPropLFM {"Edges" 2 (¤dimG2) }
  setFieldProperty $gr $fn WidgetHeight 5

  set fn MESH_BG_MESH_FILE
  setFieldProperty $gr $fn ActivityOnlyByMask 0

#===================#
# Initial condition #
#===================#

  #-Panel variables
  set InitialCondition(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
    VELOCITY_1 VELOCITY_2 VELOCITY_3
    KINETIC_ENERGY
    KINETIC_DISSIPATION 
    PRESSURE
    TEMPERATURE
    DISPLACEMENT_1 DISPLACEMENT_2 DISPLACEMENT_3
    DIFFUSION
  }


  set ic InitialCondition
  set IC $InitialCondition(parameterType)

  set fn INCLUDE
  set Common($IC,$fn,editProc) exists
  setFieldProperties1 $ic $fn $stdPropLF {"Parameter file" 0}


#====================#
# Material Parameter #
#====================#

  #-Panel variables
  set Material(initialFields) { 
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
    DENSITY
    VISCOSITY
    HEAT_EXPANSION_COEFFICIENT
    REFERENCE_TEMPERATURE 
    KE_CMU
    KE_C1
    KE_C2
    KE_SIGMAK
    KE_SIGMAE
    COMPRESSIBILITY_MODEL
    REFERENCE_PRESSURE 
    SPECIFIC_HEAT_RATIO
    HEAT_CAPACITY
    HEAT_CONDUCTIVITY
    ENTHALPY
    PHASE_CHANGE_INTERVALS
    CONVECTION_VELOCITY_1 CONVECTION_VELOCITY_2 CONVECTION_VELOCITY_3
    YOUNGS_MODULUS
    POISSON_RATIO
    DIFFUSIVITY
  }


  set mt Material
  set MT $Material(parameterType)

  set fn INCLUDE
  set Common($MT,$fn,editProc) exists
  setFieldProperties1 $mt $fn $stdPropLF {"Parameter file" 0}


#============#
# MeshDefine #
#============#
  set MeshDefine(initialFields) {
    MESH_NAME
    MESH_INDEX
    MESH_H
    MESH_F
    MESH_BG_MESH_FILE
    MESH_BG_MESH_ACTIVE
    MESH_BG_MESH_CONTROL
  }

  set md MeshDefine
  set MD $MeshDefine(parameterType)

  # These variable names are used to make it easier to handlee MeshDefine panel
  # updates etc.
  #
  set MeshDefine(fieldNames) { 
    MESH_NAME MESH_H MESH_F
    MESH_BG_MESH_FILE MESH_BG_MESH_ACTIVE MESH_BG_MESH_CONTROL
  }

  set MeshDefine(listNames) { 
    meshNames meshHs meshFs
    meshBgMeshFiles meshBgMeshActives meshBgMeshControls
  }


#================#
# ModelParameter #
#================#

  set ModelParameter(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
  }



#===============#
# ModelProperty #
#===============#

  # Problem directory value types
  # Directories are: INCLIDE_PATH, RESULTS_DIRECTORY, LOG_DIRECTORY
  #
  set ModelProperty(valueTypes) {Session Model Settings}

  set ModelProperty(initialFields) {
    MODEL_NAME
    MODEL_DESCRIPTION
    PROBLEM_NAME
    PROBLEM_DESCRIPTION
    MODEL_DIRECTORY
    CURRENT_MESH_DIRECTORY
    INCLUDE_PATH
    RESULTS_DIRECTORY
    LOG_DIRECTORY
  }


#===========#
# Processor #
#===========#

  set Processor(initialFields) {
    NOF_PROCESSORS
  }

  set Processor(max_nof_processors) 192


#=====================#
# SimulationParameter #
#=====================#

  set SimulationParameter(initialFields) {
    INCLUDE,Use
    INCLUDE,Browse
    INCLUDE
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    PROCEDURE_BUTTON
    TABLE_BUTTON
  }


#========#
# Solver #
#========#
# This was after Variable-Dofs     TIME_DERIVATIVE_ORDER

  # Panel variables
  #
  # NOTE: This is Sif output order, not screen order
  #
  set Solver(initialFields) {
    ACTIVE
    EXEC_SOLVER
    EQUATION
    SOLVING_ORDER
    MESH
    VARIABLE
    VARIABLE_DOFS
    LIBRARY_FILE
    FUNCTION_NAME
    PROCEDURE
    INCLUDE
    INCLUDE,Use
    INCLUDE,Browse
    LINEAR_SYSTEM_LABEL
    LINEAR_SYSTEM_SOLVER
    LINEAR_SYSTEM_METHOD
    LINEAR_SYSTEM_DIRECT_METHOD
    LINEAR_SYSTEM_ITERATIVE_METHOD
    LINEAR_SYSTEM_MULTIGRID_METHOD
    LINEAR_SYSTEM_MAX_ITERATIONS
    LINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    LINEAR_SYSTEM_ABORT_NOT_CONVERGED
    LINEAR_SYSTEM_PRECONDITIONING
    LINEAR_SYSTEM_ILUT_TOLERANCE
    LINEAR_SYSTEM_RESIDUAL_OUTPUT
    STEADY_STATE_CONVERGENCE_TOLERANCE
    MG_LABEL
    MG_LEVELS
    MG_EQUAL_SPLIT
    MG_MESH_NAME
    MG_SMOOTHER
    MG_PRE_SMOOTHING_ITERATIONS
    MG_POST_SMOOTHING_ITERATIONS
    MG_MAX_ITERATIONS
    MG_CONVERGENCE_TOLERANCE
    MG_PRECONDITIONING
    MG_ILUT_TOLERANCE
    STABILIZE
    BUBBLES
    LUMPED_MASS_MATRIX
    NONLINEAR_SYSTEM_LABEL
    NONLINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    NONLINEAR_SYSTEM_MAX_ITERATIONS
    NONLINEAR_SYSTEM_NEWTON_AFTER_ITERATIONS
    NONLINEAR_SYSTEM_NEWTON_AFTER_TOLERANCE
    NONLINEAR_SYSTEM_RELAXATION_FACTOR
  }

  set sl Solver
  set SL $Solver(parameterType)

  #-Default values per system and other Solver specific settings
  #
  # Predefined equation indices
  # 0 = Navier-Stokes
  # 1 = KE Turbulence
  # 2 = Heat equation
  # 3 = Stress Analysis
  # 4 = Advection Diffusion Equation
  
  # NOTE: The first property value (without the #n system index)
  # is a generic value which is applied if nothing is specified for
  # the solver with #n index or through edf-file!!!

  set iv InitialValue
  set ia InitiallyActive
  set lw 12

  set fn EQUATION
  setFieldProperty $sl $fn SubPanel 0
  setFieldProperty $sl $fn $iv { @Solver::getEquationSifName }
  setFieldProperty $sl $fn Display 0

  set fn VARIABLE
  setFieldProperty $sl $fn SubPanel 0
  setFieldProperty $sl $fn $iv { @Solver::getEquationVariable }
  setFieldProperty $sl $fn Display 0

  set fn VARIABLE_DOFS
  setFieldProperty $sl $fn $iv { @Solver::getEquationVariableDofs }
  setFieldProperty $sl $fn Display 0

  set fn LIBRARY_FILE
  setFieldProperty $sl $fn Label "Library file"
  setFieldProperty $sl $fn LabelWidth $lw
  setFieldProperty $sl $fn WidgetAnchor w
  setFieldProperty $sl $fn WidgetWidth 32
  setFieldProperty $sl $fn SubPanel 0
  setFieldProperty $sl $fn NeedeLevel 1
  setFieldProperty $sl $fn InitiallyActive { 1 {#0 0} {#1 0} {#2 0} {#3 0} }

  set fn FUNCTION_NAME
  setFieldProperty $sl $fn Label "Function name"
  setFieldProperty $sl $fn LabelWidth $lw
  setFieldProperty $sl $fn WidgetAnchor w
  setFieldProperty $sl $fn WidgetWidth 32
  setFieldProperty $sl $fn SubPanel 0
  setFieldProperty $sl $fn NeedeLevel 1
  setFieldProperty $sl $fn InitiallyActive { 1 {#0 0} {#1 0} {#2 0} {#3 0} }

  set fn PROCEDURE
  setFieldProperty $sl $fn $iv { {#4 "AdvectionDiffusion;AdvectionDiffusionSolver"} }

  set fn EXEC_SOLVER
  setFieldProperties1 $sl $fn $stdPropLF {"Exec solver" 1}
  setFieldProperty $sl $fn $iv { "Always" }

  set fn SOLVING_ORDER
  setFieldProperty $sl $fn $iv { {#0 2} {#1 3} {#2 1} {#3 4} {#4 5} }
  setFieldProperty $sl $fn ExcludeFromSif 1
  setFieldProperty $sl $fn Display 0

  set fn MESH
  setFieldProperty $sl $fn SubPanel 0
  setFieldProperty $sl $fn $iv {"@Solver::getDefaultMesh"}
  setFieldProperty $sl $fn Display 0

  set fn ACTIVE
  setFieldProperty $sl $fn $iv {False}

  set fn INCLUDE
  setFieldProperty $sl $fn $iv {""}
  setFieldProperty $sl $fn LabelWidth $lw

  set fn LINEAR_SYSTEM_LABEL
  setFieldProperty $sl $fn Label "LINEAR SYSTEM SETTINGS"
  setFieldProperty $sl $fn WidgetType Label
  setFieldProperty $sl $fn PanelPage 1
  setFieldProperty $sl $fn SubPanel 2

  set fn LINEAR_SYSTEM_SOLVER
  setFieldProperty $sl $fn $iv { Iterative }
  setFieldProperty $sl $fn $ia { 1 }

  set fn LINEAR_SYSTEM_DIRECT_METHOD
  setFieldProperty $sl $fn $iv { "BANDED" }
  setFieldProperty $sl $fn $ia { 0 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"Direct"}

  set fn LINEAR_SYSTEM_ITERATIVE_METHOD
  setFieldProperty $sl $fn $iv { BiCGStab }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"Iterative"}

  set fn LINEAR_SYSTEM_METHOD
  setFieldProperty $sl $fn $iv { BiCGStab }

  set fn LINEAR_SYSTEM_MULTIGRID_METHOD
  setFieldProperty $sl $fn $iv { Jacobi }
  setFieldProperty $sl $fn $ia { 0 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"Multigrid"}

  set fn LINEAR_SYSTEM_CONVERGENCE_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-08 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"~Multigrid"}

  set fn LINEAR_SYSTEM_MAX_ITERATIONS
  setFieldProperty $sl $fn $iv { 300 {#0 500} {#1 500} {#2 350} {#3 300} {#4 300} }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"~Multigrid"}

  set fn LINEAR_SYSTEM_ABORT_NOT_CONVERGED
  setFieldProperty $sl $fn $iv { True }
  setFieldProperty $sl $fn AlwaysOutput { 1 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"~Multigrid"}

  set fn LINEAR_SYSTEM_PRECONDITIONING
  setFieldProperty $sl $fn $iv { ILU0 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"~Multigrid"}

  set fn LINEAR_SYSTEM_ILUT_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-03 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"~Multigrid"}

  set fn LINEAR_SYSTEM_RESIDUAL_OUTPUT
  setFieldProperty $sl $fn $iv { 1 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER}
  setFieldProperty $sl $fn ActivityValues {"~Multigrid"}

  set fn STEADY_STATE_CONVERGENCE_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-05 }

  set fn STABILIZE
  setFieldProperty $sl $fn $iv { True }
  setFieldProperty $sl $fn $ia { 1 {#1 0} {#3 0} }
  setFieldProperty $sl $fn Display { 1 {#1 0}  {#3 0} }

  set fn BUBBLES
  setFieldProperty $sl $fn $iv { False }
  setFieldProperty $sl $fn $ia { 1 {#1 0} {#3 0} }
  setFieldProperty $sl $fn Display { 1 {#1 0} {#3 0} }

  set fn LUMPED_MASS_MATRIX
  setFieldProperty $sl $fn $iv { False }

  set fn MG_LABEL
  setFieldProperty $sl $fn Label "MULTIGRID SETTINGS"
  setFieldProperty $sl $fn WidgetType Label
  setFieldProperty $sl $fn PanelPage 2
  setFieldProperty $sl $fn SubPanel 2

  set fn MG_LEVELS
  setFieldProperty $sl $fn $iv { 1 }
  setFieldProperty $sl $fn $ia { 0 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_EQUAL_SPLIT
  setFieldProperty $sl $fn $iv { True }
  setFieldProperty $sl $fn $ia { 0 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_MESH_NAME
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_SMOOTHER
  setFieldProperty $sl $fn $iv { Jacobi }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_PRE_SMOOTHING_ITERATIONS
  setFieldProperty $sl $fn $iv { 5 }
  setFieldProperty $sl $fn $ia { 0 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_POST_SMOOTHING_ITERATIONS
  setFieldProperty $sl $fn $iv { 5 }
  setFieldProperty $sl $fn $ia { 0 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_CONVERGENCE_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-08 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_MAX_ITERATIONS
  setFieldProperty $sl $fn $iv { 300 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_PRECONDITIONING
  setFieldProperty $sl $fn $iv { ILU0 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn MG_ILUT_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-03 }
  setFieldProperty $sl $fn PanelPage { 2 }
  setFieldProperty $sl $fn ActivityParents {LINEAR_SYSTEM_SOLVER LINEAR_SYSTEM_PRECONDITIONING}
  setFieldProperty $sl $fn ActivityValues { {"Multigrid"} {"Multigrid"} }

  set fn NONLINEAR_SYSTEM_LABEL
  setFieldProperty $sl $fn Label "LINEARIZATION SETTINGS"
  setFieldProperty $sl $fn WidgetType Label
  setFieldProperty $sl $fn PanelPage 1
  setFieldProperty $sl $fn SubPanel 1

  set fn NONLINEAR_SYSTEM_CONVERGENCE_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-05 }

  set fn NONLINEAR_SYSTEM_MAX_ITERATIONS
  setFieldProperty $sl $fn $iv { 1 {#0 5} }

  set fn NONLINEAR_SYSTEM_NEWTON_AFTER_ITERATIONS
  setFieldProperty $sl $fn $iv { 3 }

  set fn NONLINEAR_SYSTEM_NEWTON_AFTER_TOLERANCE
  setFieldProperty $sl $fn $iv { 1.0e-02 }

  set fn NONLINEAR_SYSTEM_RELAXATION_FACTOR
  setFieldProperty $sl $fn $iv { 1.0 }

  set fn TIME_DERIVATIVE_ORDER
  setFieldProperty $sl $fn Display 0
  setFieldProperty $sl $fn $iv {{#3 2}}

  #==============================================================================#
  # Meaning of the system indices
  #
  # 0 = Navier-Stokes
  # 1 = KE Turbulence
  # 2 = Heat equation
  # 3 = Stress Equations
  # 3 = Advection Diffusion
  #==============================================================================#
  #set SolverSystem(indices) {0 1 2 3 4}
  #set SolverSystem(labels) {"Navier-Stokes" "KE Turbulence" "Heat equation" "Stress analysis" "Advection Diffusion"}
  #set SolverSystem(names) {NAVIER-STOKES KE_TURBULENCE HEAT_EQUATION STRESS_ANALYSIS ADVECTION_DIFFUSION_EQUATION}
  #set SolverSystem(outNames) {"Navier-Stokes" "KE Turbulence" "Heat Equation" "Stress Analysis" "Advection Diffusion Equation"}

  set Solver(direct_group) {
    LINEAR_SYSTEM_DIRECT_METHOD
  }

  set Solver(iter_group) {
    LINEAR_SYSTEM_ITERATIVE_METHOD
    LINEAR_SYSTEM_MAX_ITERATIONS
    LINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    LINEAR_SYSTEM_ABORT_NOT_CONVERGED
    LINEAR_SYSTEM_PRECONDITIONING
    LINEAR_SYSTEM_ILUT_TOLERANCE
    LINEAR_SYSTEM_RESIDUAL_OUTPUT
  }

  set Solver(multigrid_group) {
    LINEAR_SYSTEM_MULTIGRID_METHOD
    MG_LEVELS
    MG_EQUAL_SPLIT
    MG_SMOOTHER
    MG_PRE_SMOOTHING_ITERATIONS
    MG_POST_SMOOTHING_ITERATIONS
    MG_MESH_NAME
    MG_CONVERGENCE_TOLERANCE
    MG_MAX_ITERATIONS
    MG_PRECONDITIONING
    MG_ILUT_TOLERANCE
  }

  set Solver(commonFields) {
    EXEC_SOLVER
    STABILIZE
    BUBBLES
    LUMPED_MASS_MATRIX
  }

  set Solver(linearSystemFields) {
    LINEAR_SYSTEM_LABEL
    LINEAR_SYSTEM_SOLVER
    LINEAR_SYSTEM_METHOD
    LINEAR_SYSTEM_DIRECT_METHOD
    LINEAR_SYSTEM_ITERATIVE_METHOD
    LINEAR_SYSTEM_MAX_ITERATIONS
    LINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    LINEAR_SYSTEM_ABORT_NOT_CONVERGED
    LINEAR_SYSTEM_PRECONDITIONING
    LINEAR_SYSTEM_ILUT_TOLERANCE
    LINEAR_SYSTEM_RESIDUAL_OUTPUT
    STEADY_STATE_CONVERGENCE_TOLERANCE
  }

  set Solver(nonlinearSystemFields) {
    NONLINEAR_SYSTEM_LABEL
    NONLINEAR_SYSTEM_CONVERGENCE_TOLERANCE
    NONLINEAR_SYSTEM_MAX_ITERATIONS
    NONLINEAR_SYSTEM_NEWTON_AFTER_ITERATIONS
    NONLINEAR_SYSTEM_NEWTON_AFTER_TOLERANCE
    NONLINEAR_SYSTEM_RELAXATION_FACTOR
  }

  set Solver(multigridFields) {
    MG_LABEL
    MG_LEVELS
    MG_EQUAL_SPLIT
    MG_SMOOTHER
    MG_PRE_SMOOTHING_ITERATIONS
    MG_POST_SMOOTHING_ITERATIONS
    MG_MESH_NAME
    MG_CONVERGENCE_TOLERANCE
    MG_MAX_ITERATIONS
    MG_PRECONDITIONING
    MG_ILUT_TOLERANCE
    LINEAR_SYSTEM_MULTIGRID_METHOD
  }

  # Make the common fields visible for all predefined solvers!!!
  #
  set m1 ((¤H)|(¤F)|(¤S)|(¤AD))
  set m2 ((¤H)|(¤F)|(¤S)|(¤AD))

  foreach fn $Solver(commonFields) {
    setFieldProperty Solver $fn Mask [list $m1 $m2]
  }
  foreach fn $Solver(linearSystemFields) {
    setFieldProperty Solver $fn Mask [list $m1 $m2]
  }
  foreach fn $Solver(nonlinearSystemFields) {
    setFieldProperty Solver $fn Mask [list $m1 $m2]
  }
  foreach fn $Solver(multigridFields) {
    setFieldProperty Solver $fn Mask [list $m1 $m2]
  }
  


#================#
# SolverControl  #
#================#

  set sc SolverControl
  set SC $SolverControl(parameterType)

  set fn EXEC_SOLVER
  setFieldProperties1 $sc $fn $stdPropLF {"Reload file" 1}
  setFieldProperty $sc $fn $iv { "Always" }
  setFieldProperty $sc $fn "ExcludeFromSif" { 1 }

  set SolverControl(initialFields) {
    CLEAR_PANELS_BUTTON
    CLEAR_PANEL_BUTTON
    CLEAR_ENTRY_BUTTON
    EDIT_BUTTON
    TABLE_BUTTON
    ECHO_ON
    CHECK_KEYWORDS
    MIN_OUTPUT_LEVEL
    MAX_OUTPUT_LEVEL
    OUTPUT_LEVEL
    OUTPUT_CALLER
    RELOAD_INPUT_FILE
    EXEC_SOLVER
  }


#==========#
# Timestep #
#==========#

  set Timestep(initialFields) {
    ACTIVE
    TIMESTEPPING_METHOD
    NEWMARK_BETA
    BDF_ORDER
    TIMESTEP_INTERVALS
    TIMESTEP_SIZES
    OUTPUT_INTERVALS
    SIMULATION_TYPE
    STEADY_STATE_MAX_ITERATIONS
    STEADY_STATE_OUTPUT_INTERVAL
  }



#===============#
# User settings #
#===============#

  set UserSetting(initialFields) {
    DEFAULT_MODEL_DIRECTORY
    DEFAULT_CAD_FILES_DIRECTORY
    DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY
    DEFAULT_INCLUDE_PATH
    DEFAULT_RESULTS_DIRECTORY
    DEFAULT_LOG_DIRECTORY
    DEFAULT_USE_MODEL_SETTINGS
    DEFAULT_AUTO_SAVE_EXTERNAL_MESH
    AUTO_LOAD_MESH
    AUTO_SAVE_MODEL
    AUTO_SAVE_SOLVER_INPUT
    BROWSER_COMMAND
    EDITOR_COMMAND
    BROWSE_MODE_GEBHARDT_FACTORS
    BROWSE_MODE_MESH
    BROWSE_MODE_SOLVER
    BROWSE_MODE_VIEW_FACTORS
    BROWSE_MODE_PROCEDURE_COMPILER
    FONT_SIZES
  }


  set us UserSetting
  set US $UserSetting(parameterType)

  set fn BROWSE_MODE_GEBHARDT_FACTORS
  set Common($US,$fn,fillProc) exists
  set Common($US,$fn,group) {BM_GEBHARDT_FACTORS_LOG BM_GEBHARDT_FACTORS_SHELL BM_GEBHARDT_FACTORS_NONE }

  set fn BROWSE_MODE_MESH
  set Common($US,$fn,fillProc) exists
  set Common($US,$fn,group) {BM_MESH_LOG BM_MESH_SHELL BM_MESH_NONE }

  set fn BROWSE_MODE_PROCEDURE_COMPILER
  set Common($US,$fn,fillProc) exists
  set Common($US,$fn,group) {BM_PROCEDURE_COMPILER_LOG BM_PROCEDURE_COMPILER_SHELL BM_PROCEDURE_COMPILER_NONE }

  set fn BROWSE_MODE_SOLVER
  set Common($US,$fn,fillProc) exists
  set Common($US,$fn,group) {BM_SOLVER_LOG BM_SOLVER_SHELL BM_SOLVER_NONE }

  set fn BROWSE_MODE_VIEW_FACTORS
  set Common($US,$fn,fillProc) exists
  set Common($US,$fn,group) {BM_VIEW_FACTORS_LOG BM_VIEW_FACTORS_SHELL BM_VIEW_FACTORS_NONE }

  # Not shown on settings screen!
  set fn FONT_SIZES
  setFieldProperty UserSetting $fn Display 0


#===========#
# Variables #
#===========#

  set Variables(initialFields) {
    COORDINATE COORDINATE_1 COORDINATE_2 COORDINATE_3
    DISPLACEMENT DISPLACEMENT_1 DISPLACEMENT_2 DISPLACEMENT_3
    KINETIC_ENERGY KINETIC_DISSIPATION
    PRESSURE
    TEMPERATURE
    TIME
    VELOCITY VELOCITY_1 VELOCITY_2 VELOCITY_3
  }


  # Check that both problem and display masks will be set
  #======================================================
  foreach fn $Common(allFieldNames) {

    set masks [getFieldProperty "" $fn Mask]

    if { 1 == [llength $masks] } {
      lappend masks $masks
      setFieldProperty "" $fn Mask $masks
    }
  }


  # Set BoundaryCondition specific masks
  # ====================================
  # NOTE: This cannot be done until BOTH problem and display masks are set!!
  set m ""
  lappend m {}
  lappend m &(¤bcD)

  DataField::addFieldMaskItem BoundaryCondition \
                              $BoundaryCondition(initialDirichletVnames) \
                              $m

  # NOTE: Neumann type conditions can be set to
  # faces and edges in 3D.
  # The user must know what to do with edges!
  # 
  set m ""
  lappend m {}
  lappend m &(¤bcN)

  DataField::addFieldMaskItem BoundaryCondition \
                              $BoundaryCondition(initialNeumannVnames) \
                              $m

  namespace forget DataField::getFieldPropert* DataField::setFieldPropert*
}
### End FRONT_SET_INITIAL_PANEL_FIELD_PROPERITES proc




#=================================#
#                                 # 
#  FRONT_UPDATE_FIELD_PROPERTIES  #
#                                 # 
#=================================#
#
# Update field property definitions after generic settings #
#
proc FRONT_UPDATE_FIELD_PROPERTIES {} {
  global Common


  #=========================#
  # Field's activity slaves #
  #=========================#
  #
  # NOTE: Only standard arries processed here!
  #
  foreach arr $Common(standardPanelArries) {
    
    upvar #0 $arr theArr
    
    foreach fld $theArr(initialFields) {
      set parents [DataField::getFieldProperty $arr $fld ActivityParents]
      foreach parent $parents {
        set slaves [DataField::getFieldProperty $arr $parent ActivitySlaves]
        lappend slaves $fld
        DataField::setFieldProperty $arr $parent ActivitySlaves $slaves
      }
    }
  }

}
### End FRONT_UPDATE_FIELD_PROPERITES proc




#==========================================#
#                                          #
#        FRONT_SET_PANEL_SETTINGS          # 
#                                          #
#==========================================#
#
# Set panel screen and other control values
#
proc FRONT_SET_PANEL_SETTINGS {} {
  global Common Info Model ObjectTable MatcDefinitions
  global Status ModelFlags 
  global BodyProperty BodyDisplay LabelDisplay
  global Boundary BoundaryDisplay VertexDisplay MeshDefine MeshSelect PostFileSelect
  global ProcessIds ProcessTable
  global InputFileInfo ModelInfo BodyInfo
  global RendererInfo SystemInfo SourceList TclCommand

  # Create some more global arrynames
  foreach arr $Common(allPanelArries) {
    global $arr
  }

  #--Fixed size for values area in the standard panels
  set wid [ expr ceil( (450.0 / 1024) * $Info(maxWinSizeX) )]

  if { $wid < 450 } {
    set wid 450
  }

  set hig [ expr ceil( (250.0 /  768) * $Info(maxWinSizeY) ) ]

  if { $hig < 250 } {
    set hig 250
  }

  set Common(valuesFrameW) $wid
  set Common(valuesFrameH) $hig

  #====================================================================

  # Set deafult values for panel array vars
  # =======================================
  #
  foreach arr $Common(allPanelArries) {
    upvar #0 $arr theArr
		
		set theArr(panelEqNbr) 0
		set theArr(panelPageNbr) 1

    set theArr(labelLength) 17

    set theArr(hasProcButton) 0
    set theArr(hasTableButton) 0

    set theArr(panelsByEquation) 0

    set theArr(hasBodies) 0
    set theArr(hasBoundaries) 0
    set theArr(hasParameters) 0
    set theArr(hasTarget) 0

    set theArr(hasMaskedFields) 0
    set theArr(hasMaskedParameter) 0
    set theArr(hasMaskedTarget) 0
    set theArr(acceptEmptyMask) 0

    set theArr(canMultiSelectObjects) 0
    set theArr(canMultiSelectBoundaries) 0

    set theArr(hasAddProc) 0
    set theArr(hasApplyProc) 0

    set theArr(panelSizeKnown) 0

    set theArr(currentProcedureDir) ""
  }

  #-----------------
  # BODY FORCE panel
  #-----------------
  set BodyForce(winTitle) "Body forces"
  set BodyForce(winName) .bfw
  set BodyForce(winId) -1
  set BodyForce(winGeometry) -20+35

  set BodyForce(hasProcButton) 1
  set BodyForce(hasTableButton) 1
  set BodyForce(parameterName) "body force definition"
  set BodyForce(parameterEntryName) "A body force definition"
  set BodyForce(objectHdr) "Bodies:"
  set BodyForce(valuesHdrBase) "body forces:"
  set BodyForce(paramsHdr) "Body force sets:"

  set BodyForce(panelsByEquation) 1
  set BodyForce(hasBodies) 1
  set BodyForce(hasParameters) 1
  set BodyForce(canMultiSelectObjects) 1
  set BodyForce(hasTarget) 1
  set BodyForce(hasMaskedFields) 1
  set BodyForce(hasMaskedParameter) 1
  set BodyForce(hasMaskedTarget) 1
  set BodyForce(panelOkProc) bodyForcePanelOk
  set BodyForce(parameterId) "-Force" 
  set BodyForce(defaultName) "BodyForce" 


  #---------------------
  # BODY PARAMETER panel
  #---------------------
  set BodyParameter(winTitle) "Body parameters"
  set BodyParameter(winName) .bParamw
  set BodyParameter(winId) -1
  set BodyParameter(winGeometry) -20+35

  set BodyParameter(hasProcButton) 1
  set BodyParameter(hasTableButton) 1
  set BodyParameter(parameterName) "body parameter definition"
  set BodyParameter(parameterEntryName) "A body parameter"
  set BodyParameter(objectHdr) "Bodies:"
  set BodyParameter(valuesHdrBase) "Body parameters:"
  set BodyParameter(paramsHdr) "Body paramerter sets:"

  set BodyParameter(hasBodies) 1
  set BodyParameter(hasParameters) 1
  set BodyParameter(canMultiSelectObjects) 1
  set BodyParameter(hasTarget) 1
  set BodyParameter(panelOkProc) bodyParameterPanelOk
  set BodyParameter(parameterId) "-Parameter" 
  set BodyParameter(defaultName) "BodyParameter" 


  #-------------
  # BODIES panel
  #-------------
  set BodyProperty(winTitle) "Bodies"
  set BodyProperty(winName) .bw
  set BodyProperty(winId) -1
  set BodyProperty(winGeometry) -100+100


  #-------------------
  # BODY DISPLAY panel
  #-------------------
  set BodyDisplay(winTitle) "Display bodies"
  set BodyDisplay(winName) .bodyDispw
  set BodyDisplay(winId) -1
  set BodyDisplay(winGeometry) -100+100


  #----------------
  # BODY INFO panel
  #----------------
  set BodyInfo(winTitle) "Body info"
  set BodyInfo(winName) .biw
  set BodyInfo(winId) -1
  set BodyInfo(winGeometry) +420+120


  #-------------------------
  # BOUNDARY CONDITION panel
  #-------------------------
  set BoundaryCondition(winTitle) "Boundary conditions"
  set BoundaryCondition(winName) .bcw
  set BoundaryCondition(winGeometry) -15+10
  set BoundaryCondition(winId) -1

  set BoundaryCondition(hasProcButton) 1
  set BoundaryCondition(hasTableButton) 1
  set BoundaryCondition(parameterName) "boundary condition"
  set BoundaryCondition(parameterEntryName) "A boundary condition"
  set BoundaryCondition(objectHdr) "Bodies:"
  set BoundaryCondition(boundaryHdr) "Boundaries:"
  set BoundaryCondition(valuesHdrBase) "boundary conditions:"
  set BoundaryCondition(paramsHdr) "Boundary condition sets:"
  set BoundaryCondition(labelLength) 25

  set BoundaryCondition(panelsByEquation) 1
  set BoundaryCondition(hasBodies) 1 
  set BoundaryCondition(hasBoundaries) 1 
  set BoundaryCondition(hasParameters) 1 
  set BoundaryCondition(canMultiSelectBoundaries) 1 
  set BoundaryCondition(hasTarget) 1
  set BoundaryCondition(hasMaskedFields) 1
  set BoundaryCondition(hasMaskedParameter) 1
  set BoundaryCondition(hasMaskedTarget) 1
  set BoundaryCondition(hasParameterTargetCheckProc) 1 ;# For RADIATION check!
  set BoundaryCondition(panelOkProc) boundaryConditionPanelOk
  set BoundaryCondition(parameterId) "-Cond" 
  set BoundaryCondition(defaultName) "Constraint" 


  #-----------------
  # BOUNDARIES panel
  #-----------------
  set Boundary(winTitle) "Boundaries"
  set Boundary(winName) .bndrw
  set Boundary(winId) -1
  set Boundary(winGeometry) -15+10
  #set Boundary(winGeometry) -100+100

  set Boundary(parameterName) ""
  set Boundary(objectHdr) "Bodies:"
  set Boundary(boundaryHdr) "Boundaries:"
  set Boundary(valuesHdrBase) ""
  set Boundary(paramsHdr) ""
  set Boundary(labelLength) 25


  set Boundary(defaultName) "Boundary" 

  set Boundary(hasBodies) 1 
  set Boundary(hasBoundaries) 1 
  set Boundary(canMultiSelectBoundaries) 1 
  set Boundary(hasTarget) 1
  set Boundary(panelOkProc) boundaryPanelOk
  set Boundary(parameterId) "" 


  #-----------------------
  # BOUNDARY DISPLAY panel
  #-----------------------
  set BoundaryDisplay(winTitle) "Display boundaries"
  set BoundaryDisplay(winName) .bndrDispw
  set BoundaryDisplay(winId) -1
  set BoundaryDisplay(winGeometry) -100+100


  #-------------------------
  # BOUNDARY PARAMETER panel
  #-------------------------
  set BoundaryParameter(winTitle) "Boundary parameters"
  set BoundaryParameter(winName) .bndrParamw
  set BoundaryParameter(winGeometry) -15+10
  set BoundaryParameter(winId) -1

  set BoundaryParameter(hasProcButton) 1
  set BoundaryParameter(hasTableButton) 1
  set BoundaryParameter(parameterName) "boundary parameter"
  set BoundaryParameter(parameterEntryName) "A boundary parameter"
  set BoundaryParameter(objectHdr) "Bodies:"
  set BoundaryParameter(boundaryHdr) "Boundaries:"
  set BoundaryParameter(valuesHdrBase) "Boundary parameters:"
  set BoundaryParameter(paramsHdr) "Boundary parameter sets:"
  set BoundaryParameter(labelLength) 25

  set BoundaryParameter(hasBodies) 1 
  set BoundaryParameter(hasBoundaries) 1 
  set BoundaryParameter(hasParameters) 1 
  set BoundaryParameter(canMultiSelectBoundaries) 1 
  set BoundaryParameter(hasTarget) 1
  set BoundaryParameter(panelOkProc) boundaryParameterPanelOk
  set BoundaryParameter(parameterId) "-Param" 
  set BoundaryParameter(defaultName) "BoundaryParameter" 


  #-----------------
  # CALCULATOR panel
  #-----------------
  set Calculator(winTitle) "Calculators"
  set Calculator(winName) .clw
  set Calculator(winId) -1
  set Calculator(winGeometry) +400+100

  set Calculator(parameterName) "calculator"
  set Calculator(parameterEntryName) "A calculator"
  set Calculator(paramsHdr) "Calculators"
  set Calculator(valuesHdr) "Calculator settings:"
  set Calculator(defaultName) "Calculator" 

  set Calculator(hasParameters) 1
  set Calculator(acceptEmptyMask) 1
  set Calculator(panelOkProc) calculatorPanelOk


  #---------------
  # CONSTANT panel
  #---------------
  set Constant(winTitle) "Physical constants"
  set Constant(winName) .cnw
  set Constant(winId) -1
  set Constant(winGeometry) +400+100

  set Constant(valuesHdr) "Physical constants and parameters"
  set Constant(defaultName) "Constant" 

  set Constant(parameterId) 1
  set Constant(acceptEmptyMask) 1
  set Constant(panelOkProc) constantPanelOk


  #-----------------
  # COORDINATE panel
  #-----------------
  set Coordinate(winTitle) "Coordinate settings"
  set Coordinate(winName) .crdw
  set Coordinate(winId) -1
  set Coordinate(winGeometry) +400+100


  #---------------
  # DATAFILE panel
  #---------------
  set Datafile(winTitle) "Datafiles"
  set Datafile(winName) .dfw
  set Datafile(winId) -1
  set Datafile(winGeometry) +400+100

  set Datafile(defaultName) "Datafile" 


  #---------------
  # EQUATION panel
  #---------------
  set Equation(winTitle) "Equations"
  set Equation(winName) .eqw
  set Equation(winId) -1
  set Equation(winGeometry) -20+35

  set Equation(hasProcButton) 1
  set Equation(hasTableButton) 1
  set Equation(parameterName) "body equation definition"
  set Equation(parameterEntryName) "An equation definition"
  set Equation(objectHdr) "Bodies:"
  set Equation(valuesHdrBase) "settings:"
  set Equation(paramsHdr) "Equation sets:"

  set Equation(panelsByEquation) 1
  set Equation(hasBodies) 1
  set Equation(hasParameters) 1
  set Equation(canMultiSelectObjects) 1
  set Equation(hasTarget) 1
  set Equation(hasMaskedFields) 1
  set Equation(hasMaskedTarget) 1
  set Equation(hasAddProc) 1
  set Equation(hasApplyProc) 1
  set Equation(panelOkProc) equationPanelOk
  set Equation(parameterId) "-Eq" 
  set Equation(defaultName) "Equation" 


  #-----------------------
  # EQUATIONVARIABLE panel
  #-----------------------
  set EquationVariable(winTitle) "Equation Variables"
  set EquationVariable(winName) .eqvw
  set EquationVariable(winId) -1
  set EquationVariable(winGeometry) +400+100

  set EquationVariable(parameterId) 1
  set EquationVariable(defaultName) "EquationVariable" 
  set EquationVariable(dataWasChanged) 0


  #-------------
  # GRID H panel
  #-------------
  set GridH(winTitle) "Mesh density"
  set GridH(winName) .ghw
  set GridH(winId) -1
  set GridH(winGeometry) -20+35

  set GridH(hasParameterTargetCheckProc) 1 ;# For MESH_N and vertices!
  set GridH(parameterName) "mesh density parameter"
  set GridH(parameterEntryName) "A mesh density parameter"
  set GridH(objectHdr) "Bodies:"
  set GridH(boundaryHdr) "Boundaries:"
  set GridH(valuesHdrBase) "Values for mesh density parameter:"
  set GridH(paramsHdr) "Mesh density sets:"

  set GridH(hasBodies) 1
  set GridH(hasBoundaries) 1
  set GridH(hasParameters) 1
  set GridH(canMultiSelectBoundaries) 1
  set GridH(hasTarget) 1
  set GridH(hasMaskedFields) 1
  set GridH(hasMaskedParameter) 1
  set GridH(hasMaskedTarget) 1
  set GridH(acceptEmptyMask) 1
  set GridH(hasAddProc) 1
  set GridH(hasApplyProc) 1
  set GridH(panelOkProc) gridHPanelOk
  set GridH(parameterId) "-MeshC" 
  set GridH(defaultName) "MeshDensity" 


  #---------------------
  # GRID PARAMETER panel
  #---------------------
  set GridParameter(winTitle) "Mesh structure"
  set GridParameter(winName) .grw
  set GridParameter(winId) -1
  set GridParameter(winGeometry) -20+35

  set GridParameter(hasParameterTargetCheckProc) 1 ;# For STRUCTURED mesh check!
  set GridParameter(parameterName) "mesh structure parameter"
  set GridParameter(parameterEntryName) "A mesh structure parameter"
  set GridParameter(objectHdr) "Bodies:"
  set GridParameter(valuesHdrBase) "Values for mesh structure parameter:"
  set GridParameter(paramsHdr) "Mesh structure sets:"

  set GridParameter(hasBodies) 1
  set GridParameter(hasParameters) 1
  set GridParameter(hasTarget) 1
  set GridParameter(hasMaskedFields) 1
  set GridParameter(hasMaskedParameter) 1
  set GridParameter(hasMaskedTarget) 1
  set GridParameter(hasAddProc) 1
  set GridParameter(hasApplyProc) 1
  set GridParameter(panelOkProc) gridParameterPanelOk
  set GridParameter(parameterId) "-Mesh" 
  set GridParameter(defaultName) "MeshStructure" 


  #------------------------
  # INITIAL CONDITION panel
  #------------------------
  set InitialCondition(winTitle) "Initial conditions"
  set InitialCondition(winName) .icw
  set InitialCondition(winId) -1
  set InitialCondition(winGeometry) -20+35

  set InitialCondition(parameterName) "initial condition"
  set InitialCondition(parameterEntryName) "An initial condition"
  set InitialCondition(objectHdr) "Bodies:"
  set InitialCondition(valuesHdrBase) "initial conditions:"
  set InitialCondition(paramsHdr) "Initial condition sets:"

  set InitialCondition(hasProcButton) 1
  set InitialCondition(hasTableButton) 1
  set InitialCondition(panelsByEquation) 1
  set InitialCondition(hasBodies) 1
  set InitialCondition(hasParameters) 1
  set InitialCondition(canMultiSelectObjects) 1
  set InitialCondition(hasTarget) 1
  set InitialCondition(hasMaskedFields) 1
  set InitialCondition(hasMaskedParameter) 1
  set InitialCondition(hasMaskedTarget) 1
  set InitialCondition(panelOkProc) initialConditionPanelOk
  set InitialCondition(parameterId) "-Cond" 
  set InitialCondition(defaultName) "InitialCondition" 

  #----------------------
  # INPUT FILE INFO panel
  #----------------------
  set InputFileInfo(winTitle) "Input file info"
  set InputFileInfo(winName) .miw
  set InputFileInfo(winId) -1
  set InputFileInfo(winGeometry) +420+120


  #--------------------
  # LABEL DISPLAY panel
  #--------------------
  set LabelDisplay(winTitle) "Labels"
  set LabelDisplay(winName) .ldw
  set LabelDisplay(winId) -1
  set LabelDisplay(winGeometry) -100+100

  #-----------------------
  # MATC DEFINITIONS panel
  #-----------------------
  set MatcDefinitions(winTitle) "Define matc functions and variables"
  set MatcDefinitions(winName) .matcDefsw
  set MatcDefinitions(winId) -1
  set MatcDefinitions(winGeometry) -100+100


  #--------------------------
  # MATERIAL PROPERTIES panel
  #--------------------------
  set Material(winTitle) "Materials"
  set Material(winName) .mtw
  set Material(winId) -1
  set Material(winGeometry) -20+35

  set Material(parameterName) "material definition"
  set Material(parameterEntryName) "A material definition"
  set Material(objectHdr) "Bodies:"
  set Material(valuesHdrBase) "material properties:"
  set Material(paramsHdr) "Material property sets:"

  set Material(hasProcButton) 1
  set Material(hasTableButton) 1
  set Material(panelsByEquation) 1
  set Material(hasBodies) 1
  set Material(hasParameters) 1
  set Material(canMultiSelectObjects) 1
  set Material(hasTarget) 1
  set Material(hasMaskedFields) 1
  set Material(hasMaskedParameter) 1
  set Material(hasMaskedTarget) 1
  set Material(panelOkProc) materialPanelOk
  set Material(parameterId) "-Mater" 
  set Material(defaultName) "Material" 


  #------------------
  # MESH DEFINE panel
  #------------------
  set MeshDefine(winTitle) "Define mesh"
  set MeshDefine(winName) .mHw
  set MeshDefine(winId) -1
  set MeshDefine(winGeometry) -100+100


  #-----------------
  # MODEL INFO panel
  #-----------------
  set ModelInfo(winTitle) "Model info"
  set ModelInfo(winName) .miw
  set ModelInfo(winId) -1
  set ModelInfo(winGeometry) +420+120


  #----------------------
  # MODEL PARAMETER panel
  #----------------------
  set ModelParameter(winTitle) "Model parameters"
  set ModelParameter(winName) .simupw
  set ModelParameter(winId) -1
  set ModelParameter(winGeometry) +400+100

  set ModelParameter(hasProcButton) 1
  set ModelParameter(hasTableButton) 1
  set ModelParameter(valuesHdr) "Model parameters"
  set ModelParameter(defaultName) "Model parameter" 

  set ModelParameter(parameterId) 1

  set ModelParameter(acceptEmptyMask) 1
  set ModelParameter(panelOkProc) modelParameterPanelOk


  #-----------------------
  # MODEL PROPERTIES panel
  #-----------------------
  set ModelProperty(winTitle) "Model name and directories"
  set ModelProperty(winName) .mpw
  set ModelProperty(winId) -1
  set ModelProperty(winGeometry) -100+100


  #----------------------
  # POSTFILE SELECT panel
  #----------------------
  set PostFileSelect(winTitle) "Select result file"
  set PostFileSelect(winName) .pfsw
  set PostFileSelect(winId) -1
  set PostFileSelect(winGeometry) +420+120


  #----------------
  # PROCESSOR panel
  #----------------
  set Processor(winTitle) "Processor settings"
  set Processor(winName) .prw
  set Processor(winId) -1
  set Processor(winGeometry) +420+120


  #--------------------
  # PROCESS TABLE panel
  #--------------------
  set ProcessTable(winTitle) "Process table"
  set ProcessTable(winName) .prtw
  set ProcessTable(winId) -1
  set ProcessTable(winGeometry) +420+120


  #---------------------------
  # SIMULATION PARAMETER panel
  #---------------------------
  set SimulationParameter(winTitle) "Simulation parameters"
  set SimulationParameter(winName) .simupw
  set SimulationParameter(winId) -1
  set SimulationParameter(winGeometry) +400+100

  set SimulationParameter(hasProcButton) 1
  set SimulationParameter(hasTableButton) 1
  set SimulationParameter(valuesHdr) "Simulation parameters"
  set SimulationParameter(defaultName) "Simulation" 

  set SimulationParameter(parameterId) 1

  set SimulationParameter(acceptEmptyMask) 1
  set SimulationParameter(panelOkProc) simulationParameterPanelOk


  #-------------
  # SOLVER panel
  #-------------
  set Solver(winTitle) "Solver settings"
  set Solver(winName) .slw
  set Solver(winId) -1
  set Solver(winGeometry) +420+120
  set Solver(valuesHdrBase) "solver parameters"

  set Solver(panelsByEquation) 1
  set Solver(targetMask) ""
  set Solver(hasMaskedFields) 1
  set Solver(hasParameters) 1
  set Solver(acceptEmptyMask) 1
  set Solver(hasMaskedTarget) 1
  set Solver(panelOkProc) solverParameterPanelOk
  set Solver(defaultName) "Solver" 

  #---------------------
  # SOLVER CONTROL panel
  #---------------------
  set SolverControl(winTitle) "Solver control"
  set SolverControl(winName) .slcw
  set SolverControl(winId) -1
  set SolverControl(winGeometry) +400+100

  set SolverControl(hasTableButton) 1
  set SolverControl(valuesHdr) "Solver control parameters"
  set SolverControl(defaultName) "Solver control" 

  set SolverControl(parameterId) 1

  set SolverControl(acceptEmptyMask) 1
  set SolverControl(panelOkProc) solverControlPanelOk

  #-------------------
  # SOLVER ORDER panel
  #-------------------
  set SolverOrder(winTitle) "Equation solving order"
  set SolverOrder(winName) .sow
  set SolverOrder(winId) -1
  set SolverOrder(winGeometry) +420+120


  #------------------
  # SYSTEM INFO panel
  #------------------
  set SystemInfo(winTitle) "System info"
  set SystemInfo(winName) .siw
  set SystemInfo(winId) -1
  set SystemInfo(winGeometry) +420+120


  #---------------
  # TIMESTEP panel
  #---------------
  set Timestep(winTitle) "Timestep settings"
  set Timestep(winName) .tsw
  set Timestep(winId) -1
  set Timestep(winGeometry) +420+120

  set Timestep(defaultName) "Timestep" 


  #-------------------
  # USER SETTING panel
  #-------------------
  set UserSetting(winTitle) "Settings"
  set UserSetting(winName) .usw
  set UserSetting(winId) -1
  set UserSetting(winGeometry) -20+25
  set UserSetting(valuesHdr) "Values for user settings:"

  set UserSetting(acceptEmptyMask) 1
  set UserSetting(panelOkProc) userSettingPanelOk

  #-----------------------
  # VERTEX ISPLAY panel
  #-----------------------
  set VertexDisplay(winTitle) "Display vertices"
  set VertexDisplay(winName) .vrtxDispw
  set VertexDisplay(winId) -1
  set VertexDisplay(winGeometry) -100+100

}
### End FRONT_SET_PANEL_SETTINGS proc





###################################################
#                                                 # 
#            Field help stuff                     #
#                                                 # 
###################################################

# Accepted field help panel names (in compressed format!) and
# their gui equivalent formats
#
# 1. group: Menu based names
# 2. group: Sif based names
# 3. group: Alternative (all imaginable) versions
#
# 1. value: Panel name in input (in compressed format)
# 2. value: Gui array name
#
set Info(fieldHelpPanelNameTable) {

    # Menu names
    # ==========
    #--Edit menu
    {Bodies BodyProperty}
    {Boundaries Boundary}
    
    #--Problem menu
    {ModelNameAndDirectories ModelProperty}
    {Datafiles Datafile}
    {ModelParameters ModelParameter}
    {SimulationParameters SimulationParameter}
    {CoordinateSettings Coordinate}
    {TimestepSettings Timestep}
    {PhysicalConstants Constant}
    {Equations Equation}
    {Calculators Calculator}
    {EquationSolvingOrder SolverOrder}

    #--Model menu
    {BodyParameters BodyParameter}
    {BoundaryParameters BoundaryParameter}
    {InitialConditions InitialCondition}
    {BodyForces BodyForce}
    {Materials Material}
    {BoundaryConditions BoundaryCondition}
    
    #--Mesh menu
    {DefineMesh MeshDefine}
    {MeshDensity GridH}
    {MeshStructure GridParameter}

    #--Solver menu
    {SolverSettings Solver}
    {SolverControl SolverControl}

    # Sif names
    # =========
    {Constant Constant}
    {Equation Equation}
    {BodyForce BodyForce}
    {BoundaryCondition BoundaryCondition}
    {InitialCondition InitialCondition}
    {Material Material}
    {Simulation SimulationParameter}
    {Solver Solver}

    # Any other format
    # ================
    {EquationParameter Equation}
    {EquationParameters Equation}
    {EquationSettings Equation}
    {BodyForceParameter BodyForce}
    {BodyForceParameters BodyForce}
    {BodyForceSettings BodyForce}
    {BodyParameter BodyParameter}
    {BodyParameterSettings BodyParameter}
    {BoundaryConditionParameter BoundaryCondition}
    {BoundaryConditionParameters BoundaryCondition}
    {BoundaryConditionSettings BoundaryCondition}
    {BoundaryParameter BoundaryParameter}
    {BoundaryParameterSettings BoundaryParameter}
    {InitialConditionParameter InitialCondition}
    {InitialConditionParameters InitialCondition}
    {InitialConditionSettings InitialCondition}
    {MaterialParameter Material}
    {MaterialParameters Material}
    {MaterialParameterSettings Material}
    {Model ModelParameter}
    {ModelParameter ModelParameter}
    {ModelParameterSettings ModelParameter}
    {PhysicalConstant Constant}
    {PhysicalConstantSettings Constant}
    {SimulationParameter SimulationParameter}
    {SimulationParameterSettings SimulationParameter}
    {Solvers Solver}
    {SolverParameter Solver}
    {SolverParameters Solver}
    {SolverParameterSettings Solver}
    {SolverControlParameter SolverControl}
    {SolverControlParameters SolverControl}
    {SolverControlSettings SolverControl}
}



#===================================================#
#                                                   #
#           FRONT_READ_FIELD_HELP_FILE              #
#                                                   #
#===================================================#
#
# NOTE: No help file delivered with Elmer distribution!!!
#
#
proc FRONT_READ_FIELD_HELP_FILE {filename {warn_no_success 1} } {
  global Info

  # Try to open file: elmer-home/lib/frontFieldHelp.fhf
  #
  set cmsg ""
  if { [ catch {set in [open $filename "r"]} cmsg ] } {
    
    if { $warn_no_success} {
      Message::showMessage ""
      Message::showMessage "WARNING*** Field help file $filename not read:" $Info(wrnMsgColor)
      Message::showMessage "           $cmsg" $Info(wrnMsgColor)
      Message::showMessage ""
    }    

    return 0
  }

  Message::showMessage "Reading field help file:  $filename" $Info(remMsgColor)
  
  set in_panel 0
  set in_field 0
  set txt ""
  
  # Read field help file
  # ====================
  while {![eof $in]} {

    # Reset help text, if no more in a field
    if { !$in_field } {

      # ===============================
      # Store current field's help text
      # ===============================
      if { $txt != "" } {
        set theArray(help,$fld) $txt
      }

      set txt ""
    }

    # No comments or continuation marks inside help text!
    #
    if { $in_field } {        
      set line [Util::readInputFileLine $in 0 0]

    # Other lines may contain comments and continuation marks!
    #
    } else {
      set line [Util::readInputFileLine $in]
    }

    set line [Util::stringTrim $line]

    # Change text delimiters for somethin a bit easier to handle
    # in Tcl!
    set line [string map { "\{" "\<" } $line]
    set line [string map { "\}" "\>" } $line]

    # Keywords do not contain blanks, so this works
    set kwd [string trim [lindex $line 0]]
    set kwd_value [string trim [join [lrange $line 1 end]]]

    # Line can also start with a field name (or End Panel), but check possible
    # text delimiter first!
    #
    set fld_name $line
    set fld_name [string trim [join [split $fld_name "<"]]]

    # Include help file
    # =================
    if { !$in_panel && [string equal -nocase $kwd "include"] } {

      # If include-name is relative, pick directory from the original filename
      #
      if { [string equal -nocase "relative" [file pathtype $kwd_value]] } {
        set kwd_value [file join [file dirname $filename] $kwd_value]
      }

      FRONT_READ_FIELD_HELP_FILE $kwd_value 1
    }      

    # New panel section
    # =================
    if { !$in_field && [string equal -nocase $kwd "panel"] } {
      set globArray ""
      
      set name_tbl $Info(fieldHelpPanelNameTable)

      if { [Util::convertPanelNameToGuiFormat $name_tbl $kwd_value globArray] } {
        upvar #0 $globArray theArray
        set in_panel 1
      }
  
      continue
     }

    # Panel section ended
    # ===================
    if { $in_panel && !$in_field &&
         ([string equal -nocase $fld_name "end panel"] ||
          [string equal -nocase $fld_name "endpanel"]
         )
       } {
      set globArray ""
      set in_panel 0

      continue
    }

    # New field started
    # =================
    if {$in_panel && !$in_field} {

      set fld [DataField::fieldNameSifToGui $fld_name]
      set in_field 1
      set txt ""
 
      set pos [string first "<" $line]

      if { $pos == -1 } {
        set pos 0
      }
      
      # Remove keyword and field name from input-line
      #
      set line [string range $line $pos end]
    }

    # Read field's help text line
    # ===========================
    if { $in_panel && $in_field } {

      # Remove text delimiters
      #
      if { [string index $line 0] == "<" } {
          set line [string range $line 1 end]
      }

      if { [string index $line end] == ">" } {
          set line [string range $line 0 end-1]
          set in_field 0
      }

      # Add the new piece of text
      # =========================
      append txt " "
      append txt $line
      set txt [string trim $txt]

    }
  }
}
### End FRONT_READ_FIELD_HELP_FILE proc
