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
#Module:    ecif_tk_procsInterface.tcl
#Language:  Tcl
#Date:      28.10.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Cpp-->tcl interface procedures
#
#************************************************************************

proc Interface::acceptParameters {} {
  global Info

  set data [split $Info(arguments) $Info(argSeparator)]

  set msg [lindex $data 0]

  set pdatas [lindex $data 1]
  set pdatas [split $pdatas $Info(objectSeparator)]

  set Info(results) ""

  set w .copyParamsW
  set wgeom  -100+100

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $msg
  wm geometry $w $wgeom
  
  set Info(results) [CheckBoxList::create $w $pdatas "" 40 15]

  MenuBuild::configureMenuButtonOption File copy_parameters state normal
}
 

#-Procedure reads all parameter data when called from cpp-environment
#
proc Interface::applyAllParameterDataPre {} {
  global Info Common Model ModelFlags Boundary
  global ModelParameter SimulationParameter
  global Constant Coordinate Datafile
  global BodyParameter BoundaryParameter
  global BodyForce BoundaryCondition Constant Coordinate Datafile
  global Equation EquationVariable GridH GridParameter InitialCondition
  global Material ModelProperty Solver SolverControl Timestep UserSetting

  set Info(informAboutRemovedFields) 1
  set Model(removedFieldsList) ""

  # Init non-standard panels (except Solver). This also sets
  # default values for possible new fields and drops possible 
  # obsolate fields from parameters!
  # NOTE: These are parameters which have only one active
  # parameter id at a time. That is why we initialize also the
  # fields, so that clients can right away call Coordinate(COORDINATE_SYSTEM) 
  # etc.

  # ==============================================
  # Possible pre converse of older version data!!!
  # ==============================================
  set Info(doConversion) 0
  foreach globArray $Common(allPanelArries) {
    set msg ""
    if { [catch { Interface::preConverseInputData $globArray } msg] } {
      set msg [string map {"\"" " "} $msg]
      set aname [string toupper $globArray]
      set msg [list "NOTE: Old version data for $aname was not pre-converted because:\n" $msg]
      Message::dispMessage $msg
    }
  }
  # ==============================================

  set Info(markMainWindowTitle) 0

  GridH::updateObjectParameterId
  GridParameter::updateObjectParameterId

  set mmodified 0 ;# Model modfied
  set pmodified 0 ;# Paramters modified

  # These arries need at least one parameter instance to be saved or
  # at least instantiated !!!
  # 0 <--> do not save
  #
  set arr_defs { Constant Coordinate Datafile SimulationParameter
                 SolverControl Timestep {EquationVariable 0}
               }

  foreach arr_def $arr_defs {

    set arr [lindex $arr_def 0]
    set save 1

    if { 1 < [llength $arr_def] } {
      set save [lindex $arr_def 1]
    }

    upvar #0 $arr theArr

    set theArr(dataChanged) 0
    set theArr(dataModified) 0

    set active_id 1

    foreach id $theArr(ids) {
 
      # Last active will be the active!
      if { "True" == [DataField::getFieldValue $arr $id ACTIVE] } {
        set active_id $id
      }

      set mask [Util::getArrayValue $arr targetMask]

      set pmodified1 0
      set theArr($id,data) [DataField::checkParameterData $arr $id pmodified1 $mask]

      if {$pmodified1} {
        set pmodified 1
      }
    }

    # If no parameter define, form at least one!
    #
    if { $theArr(ids) == "" || $theArr(1,data) == "" } { 
      Panel::initFields $arr
      DataField::formDefaultNonStandardParameter $arr 1
      set theArr(ids) 1
    }

    Panel::resetFields $arr

    Panel::initFields $arr

    # Create data fields from the active paramter!
    set theArr(parameterId) $active_id
    DataField::formDataFields $arr $theArr($theArr(parameterId),data)

    if { $save } {
      execNsProc $arr panelSave 0
    }

  } ;#--For each instantiated


  # What is the idea here??? Is a grid parameter really needed
  # if we have a mesh?
  # Inactivated because it causes an error when parameters are copied
  # to a model which has a mesh, MVe 1.2.2003
  #
  if { $ModelFlags(GEOMETRY_TYPE_CAD) &&
       $Model(Mesh,exists) &&
       !$Model(Mesh,hasParameter)
     } {
    #MeshDefine::panelSave 0
  }

  set Info(markMainWindowTitle) 1

  if { $mmodified || $pmodified } {

    set Model(Front,needsUpdate) 1
    
    if { $mmodified } {
      Util::setMainWindowTitle
    }

  } else {
    MenuExec::resetUpdateFlags
  }

}


proc Interface::applyAllParameterDataPost {} {

  global Info Common Model ModelFlags Boundary
  global ModelParameter SimulationParameter
  global Constant Coordinate Datafile
  global BodyParameter BoundaryParameter
  global BodyForce BoundaryCondition Constant Coordinate Datafile
  global Equation EquationVariable GridH GridParameter InitialCondition
  global Material ModelProperty Solver SolverControl Timestep UserSetting

  set Info(informAboutRemovedFields) 1
  set Model(removedFieldsList) ""

  #-EquationVariable data needs special handling because
  # model-file variables may not contain the new (default) equation
  # variables given in user definitions (this problem results normally
  # if the user has changed the name of an equation variable compared to what
  # is saved in the model file
  #

set times ""
set t ""
set t1 [clock clicks -milliseconds]

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]
#MSG "- EquationVariable::updateFieldData $t"
  EquationVariable::updateFieldData

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- EquationVariable::formActiveParameter $t"
  EquationVariable::formActiveParameter

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- Util::cpp_exec equationVariablesPanelOk $t"
  Util::cpp_exec equationVariablesPanelOk

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- EquationVariable::updateAllVariableNames $t"
  EquationVariable::updateAllVariableNames

  #-Set default coordinate system if not given
  if { $Coordinate(COORDINATE_SYSTEM) == "None" } {
    Coordinate::setDefaultCoordinateSystem
  }

  #-Apply coord-system (eg. Model(SIMULATION_DIMENSION)!)
  Coordinate::applyCoordinateDimension $Coordinate(COORDINATE_SYSTEM)

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- Equation::setOtherMasks $t"

  #-Dimension and bc masks
  Equation::setOtherMasks

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- DataField::setVariables $t"

  #-Update variable components
  DataField::setVariables

  #-Expand fields which depend on an "indexing" variables
  # (eg. a variable defined in EquationVariable)
  DataField::setAllFields $Common(allPanelArries)
#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- StdPanelInit::updateEquationDataAndMasks $t"

  #-Create Equation parameter masks
  StdPanelInit::updateEquationDataAndMasks 

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- StdPanelInit::createObjectTableMasks $t"

  #-Create ObjectTable masks
  StdPanelInit::createObjectTableMasks

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- Equation::constructProblemMask $t"

  #-Set current problem mask value
  Equation::constructProblemMask

  #-Set problem indices based on body masks
  Equation::updateEquationIndices 

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- StdPanelInit::createObjectAndBoundaryLists $t"

  #-Init all mask related panels and check parameters
  foreach arr $Common(maskedParameters) {
    StdPanelInit::createObjectAndBoundaryLists $arr
  }

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- Panel::checkInputParameters $t"

  # Remove all non-input fields from all parameters
  Panel::checkInputParameters $Common(allPanelArries) "allFields"

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- Panel::checkParameters $t"

  set pmodified 0 ;# Parameters modified
  Panel::checkParameters $Common(maskedParameters) pmodified

  # ===============================================
  # Possible post converse of older version data!!!
  # ===============================================
  if { $Info(doConversion) } {
    foreach globArray $Common(allPanelArries) {

      set msg ""
      if { [catch { Interface::postConverseInputData $globArray } msg] } {
        set msg [string map {"\"" " "} $msg]
        set aname [string toupper $globArray]
        set msg [list "NOTE: Old version data for $aname was not post-converted because:\n" $msg]
        Message::dispMessage $msg
      }
    }
  }
  # ==============================================

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- Solver::checkSolverData $t"

  #-Solver data need special handling because there can be many solvers
  # contrary to other non-standard panel objects
  # 
  set smodified 0 ;# Solvers modified
  Solver::checkSolverData $Equation(activeIndices) smodified

  Util::updateFlagNameValues

  # Replace possible old value with Newmark stuff
  Timestep::checkTimesteppingMethodValue

  if { $Timestep(dataModified) } {
    Message::showMessage "NOTE: An obsolete Timestepping keyword has been changed!"
    Message::showMessage "Saving the model file recommended!"
    Timestep::panelSave 1
  }

  #-Other globals
  set Boundary(firstTime) 1
  set Boundary(isUpdated) 0

  set Info(markMainWindowTitle) 1

  if { $pmodified || $smodified } {
    set Model(Front,needsUpdate) 1
    Util::setMainWindowTitle
  } else {
    MenuExec::resetUpdateFlags
  }

  # If some fields were removed from the model!
  # ===========================================
  if { $Model(removedFieldsList) != "" &&
       $Model(hasUserDefinitions)
     } {
    set msg ""
    lappend msg "WARNING: The following panel fields were removed from the model!\n\n"
    lappend msg "(Possibly a definition file should be loaded BEFORE reading the model file)\n\n"
    DataField::informAboutRemovedFields $msg
  }

  set Info(informAboutRemovedFields) 0
  #set Model(removedFieldsList) ""

#set t2 [clock clicks -milliseconds]; set t [expr $t2 - $t1]; lappend times $t
#MSG "- List: $times"
#MSG "- All done $t"
}


#-Procedure reads all body-level data when called from cpp-environment
#
proc Interface::applyBodyData {} {
  global ObjectTable

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }
    
    set cpp_color $ObjectTable($id,clr)

    set ObjectTable($id,clr) "#$cpp_color"

    set ObjectTable($id,dspl) 1
  }
}


#-Procedure to read inner/outer boundary-level data 
#
proc Interface::applyBoundaryData {} {
  global Model ObjectTable
  
  # Edges as subobjects of Faces
  # Vertices as subobjects of Faces if no edges defined
  #
  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    #StdPanelInit::setObjectTableSubIds "F" "E" 0
    #StdPanelInit::setObjectTableSubIds "F" "V" 1

  # Vertices as subobjects of Edges
  #
  } else {
    #StdPanelInit::setObjectTableSubIds "E" "V"
  }

  # Set display mode
  #
  foreach id $ObjectTable(ids) {

    if { $Model(GEOMETRY_DIMENSION) == "3D" } {

      if { $ObjectTable($id,tp) != "F" &&
           $ObjectTable($id,tp) != "V" 
         } {
        continue
      }

    } else {
      if { $ObjectTable($id,tp) != "E" &&
           $ObjectTable($id,tp) != "V" 
         } {
        continue
      }

    }

    set ObjectTable($id,dspl) 1
  }

}


#-Procedure reads all model-level data when called from cpp-environment
#
proc Interface::applyModelData {} {
  global Info MeshDefine Model ModelProperty 

  # Check that input file version is sensible
  set mfn $Info(FRONT_INPUT_VERSION_NBR) ;# Input model file version
  set efn $Info(FRONT_MAIN_VERSION_NBR)  ;# Current Elmer Front main version

  if { $mfn > $efn } {

    set msg "WARNING: Elmer Front version ($efn) is older than the model file version ($mfn) !"
    Message::dispOkMessage $msg"
  }

  #--Pick possible current mesh name
  if { $Model(meshNames) != "" &&
       $Model(currentMeshIndex) != $Info(NO_INDEX)
     } {
    set old_current_mesh [lindex $Model(meshNames) $Model(currentMeshIndex)]

  } else {
    set old_current_mesh ""
  }

  set old_mesh_names $Model(meshNames)

  #--Find active mesh names
  #set Model(meshNames) [ModelProperty::findModelMeshNames]

  #--Find old mesh data by the new mesh name list
  set old2new_indices ""

  if { $old_mesh_names != "" } {

    foreach old_name $old_mesh_names {
      lappend old2new_indices  [lsearch $Model(meshNames) $old_name]
    }
  }

  #--Form mesh-hs, mesh-fs anf bg-mesh var values for
  #  current meshes
  set mesh_hs ""
  set mesh_fs ""

  set mesh_bg ""
  set mesh_bg_act ""
  set mesh_bg_ctrl ""

  # These parameter values means that a mesh
  # was not identified!
  set def_h -1
  set def_f -1
  set def_bg ""
  set def_bg_act 0
  set def_bg_ctrl 0
  
  # "Sparse" Model-array lists
  foreach mn $Model(meshNames) {
    lappend mesh_hs $def_h
    lappend mesh_fs $def_f
    lappend mesh_bg $def_bg
    lappend mesh_bg_act $def_bg_act
    lappend mesh_bg_ctrl $def_bg_ctrl
  }

  foreach idx $old2new_indices \
          mh  $Model(meshHs)   \
          mf  $Model(meshFs) { 

    # If old mesh not found any more
    if { $idx == -1 || $idx == "" } {
      continue
    }
    
    # Use model values if they are "actual"
    if { $mh != "" && $mf != "" } {
      set mesh_hs [lreplace $mesh_hs $idx $idx $mh]
      set mesh_fs [lreplace $mesh_fs $idx $idx $mf]

      # Update "sparse" bg-lists from indexed input lists!!!
      #
      set bg_idx [lsearch $Model(meshBgMeshFileIndices) $idx]

      if { $bg_idx == -1 } {
        continue
      }

      set mesh_bg [ lreplace $mesh_bg $idx $idx [lindex $Model(meshBgMeshFiles) $bg_idx] ]
      set mesh_bg_act [ lreplace $mesh_bg_act $idx $idx [lindex $Model(meshBgMeshActives) $bg_idx] ]
      set mesh_bg_ctrl [ lreplace $mesh_bg_ctrl $idx $idx [lindex $Model(meshBgMeshControls) $bg_idx] ]
    }
  }

  set Model(meshHs) $mesh_hs
  set Model(meshFs) $mesh_fs
  set Model(meshBgMeshFiles) $mesh_bg
  set Model(meshBgMeshActives) $mesh_bg_act
  set Model(meshBgMeshControls) $mesh_bg_ctrl

  # Init current mesh name, if missing
  if { $old_current_mesh == "" } {
    set old_current_mesh [lindex $Model(meshNames) 0]
  }

  #--Set current mesh index and name
  set mesh_index [lsearch $Model(meshNames) $old_current_mesh]

  if { $mesh_index != -1 } {
    set Model(currentMeshIndex) $mesh_index
    set Model(currentMeshName) [lindex $Model(meshNames) $mesh_index]

  } else {
    set Model(currentMeshIndex) $Info(NO_INDEX)
    set Model(currentMeshName) ""
  }

  #--Set directory info
  ModelProperty::setModelCaseName
  ModelProperty::setCurrentMeshDirectory
}


proc Interface::applyModelFlags {} {
  global ModelFlags

  if { $ModelFlags(DRAW_TARGET_BODIES) } {
    set ModelFlags(DRAW_TARGET) "DRAW_TARGET_BODIES"

  } elseif { $ModelFlags(DRAW_TARGET_SURFACES) } {
    set ModelFlags(DRAW_TARGET) "DRAW_TARGET_SURFACES"

  } else {
    set ModelFlags(DRAW_TARGET) "DRAW_TARGET_EDGES"
  }

}


proc Interface::applyModelGeometryDimension {} {
  global Info Model

  if { $Model(GEOMETRY_DIMENSION) == "2D" } {
    set Info(drawEdges) 1
  } else {
    set Info(drawSurfaces) 1
  }

  Message::showMessage "Model dimension: $Model(GEOMETRY_DIMENSION)"
}


# Apply model status
# 0 <--> ok, !0 <--> some error
#
proc Interface::applyModelStatus {} {
  global Model

  if { $Model(status) == 0 } {
    set Model(statusMessage) "" 
  } else {
  Util::cpp_exec putModelStatusMessage
  }
}


proc Interface::applyObjectTableData {} {
}


proc Interface::applyStatsData {} {
  global Model

  if { $Model(minEdgeSize) > 1e+010 } {
    set Model(minEdgeSize) ""
  }

}


#---Generic Debug command
#
proc Interface::debug {} {
  global UserSetting

  Message::showMessage "Models dir = $UserSetting(DEFAULT_MODEL_DIRECTORY)"
}

#---Initalization commands

proc Interface::activateFileMenus {} {
  MenuBuild::activateFileMenus
}


proc Interface::checkmeshInfoTs {} {
  global Info

  set ts $Info(arguments)
}


# Set option value in one button
#
proc Interface::configureButtonOption {} {
  global Info

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set button  [lindex $arg_list 0]
  set option  [lindex $arg_list 1]
  set value   [lindex $arg_list 2]

  MenuBuild::configureButtonOption $button $option $value
}


# Configure one or more buttons state
#
proc Interface::configureButtons {} {
  global Info

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set buttons [lindex $arg_list 0]
  set value   [lindex $arg_list 1]

  MenuBuild::configureButtons $buttons $value
}


# Set state for the menu commands in a group "menu"
#
proc Interface::configureMenuButtons {} {
  global Info

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set menu          [lindex $arg_list 0]
  set buttons       [lindex $arg_list 1]
  set value         [lindex $arg_list 2]

  MenuBuild::configureMenuButtons $menu $buttons $value
}


# Set option value in one menu button (command)
#
proc Interface::configureMenuButton {} {
  global Info

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set menu        [lindex $arg_list 0]
  set button      [lindex $arg_list 1]
  set option      [lindex $arg_list 2]
  set value       [lindex $arg_list 2]

  MenuBuild::configureMenuButton $menu $button $option $value
}


proc Interface::createBodiesMenus {} {
  global gMW
  MenuBuild::createBodyListMenu $gMW(Model,bodyList) ""
}

 
# Receive command to display error-message
#
proc Interface::displayErrorMsg {} {
  global Info __GUIServerPort

  set sep $Info(argSeparator)
  set data [split $Info(arguments) $sep]
  set err_level [lindex $data 0]

  if { $err_level == 0 } {
    set Info(messageIcon) warning
  } else {
    set Info(messageIcon) error
  }

  set msg [lindex $data 1]
  Message::dispOkMessage $msg

  # ERROR LEVEL 99 <===> STOPPING FRONT PROGRAM !!!
  if {$err_level == 99} {
    Util::cpp_exec cpp_exit
    exit
  }
}


proc Interface::displayBodySelectPanel {} {
  Util::execProc BodyDisplay::openPanel
}


proc Interface::displayBoundarySelectPanel {} {
  Util::execProc BoundaryDisplay::openPanel
}


proc Interface::displayLabelSelectPanel {} {
  Util::execProc LabelDisplay::openPanel
}


proc Interface::fieldNameGuiToSif {} {
  global Info
  
  set gui_name $Info(arguments)
  set Info(results) [DataField::fieldNameGuiToSif $gui_name]
}


proc Interface::fieldNameSifToGui {} {
  global Info

  set sif_name $Info(arguments)
  set Info(results) [DataField::fieldNameSifToGui $sif_name]
}


proc Interface::getColorName {} {
  global Info 

}


# Check if an equation panel field (possibly for a specific equation) should
# be output into a solver section
#
proc Interface::getIsSolverTargetField {} {
  global Info 

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set eq_name [lindex $arg_list 0]
  set fld_name [lindex $arg_list 1]

  set Info(results) 0

  set field_target [DataField::getFieldProperty Equation $fld_name FieldTarget]

  # Is not a redirected field
  #
  if { $field_target == "" } {
    set Info(results) 0
    return
  }

  # Ok, is redirected and equation name does not matter
  if { $eq_name == "" } {
    set Info(results) 1
    return
  }

  set eq_field [DataField::fieldNameSifToGui $eq_name]

  set emask [DataField::getFieldProperty Equation $eq_field EquationMask]

  # Is redirected but equation name does not match
  if { ![DataField::fieldDisplayMaskMatches Equation $fld_name $emask] } {
    set Info(results) 0
    return
  }

  set Info(results) 1
}


# Get mif-file name
#
proc Interface::getMeshInputFileName {} {
  global Info

  set Info(results) [MenuExec::getMeshInputFileName]
}


proc Interface::getNextActiveSelectionTolerance {} {
  global Info Model

  set Info(nextActiveSelectionTolerance) $Info(arguments)
  set Info(NEXT_ACTIVE_SELECTION_TOLERANCE) $Info(arguments)
}


proc Interface::getParallelInfo {} {
  global Info Processor

  set arg_list [split $Info(arguments) $Info(fieldSeparator)]
  set index 0

  foreach var $Processor(allFields) {
    set Processor($var) [lindex $arg_list $index]
    incr index
  }
}


proc Interface::getParameterFieldInfo {} {
  global Info

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set globArray [lindex $arg_list 0]
  set fld [lindex $arg_list 1]

  set sif_name [DataField::getFieldProperty $globArray $fld SifName "" 0]
  set vtype [DataField::getFieldProperty $globArray $fld FieldValueType]
  set output_sif_type [DataField::getFieldProperty $globArray $fld OutputSifType]
  set always_output [DataField::getFieldProperty $globArray $fld AlwaysOutput]
  set dtype [DataField::getFieldProperty $globArray $fld FieldDataType]

  # NOTE: Argument order in the list is fixed!!!
  # (ref. UserInterface_TCL::getParameterInfo() in cpp side)
  #
  set Info(results) ""

  #--1. Possible specially set Sif-name
  lappend Info(results) $sif_name

  #--2. Field value type (Integer Real etc.)
  lappend Info(results) $vtype

  #--3. Is sif output field (0/1)
  if { [Panel::isSifField $globArray $fld] } {
    lappend Info(results) 1
  } else {
    lappend Info(results) 0
  }

  #--4. Output field value type to sif (0/1)
  lappend Info(results) $output_sif_type

  #--5. Is field should be always output (even when value is (False/None)) (0/1)
  lappend Info(results) $always_output

  #--6. Is array (0/1)
  if { $dtype == "Scalar" } {
    lappend Info(results) 0
  } else {
    lappend Info(results) 1
  }

  #--7. Is quoted (0/1)
  if { $vtype == "String" } {
    lappend Info(results) 1
  } else {
    lappend Info(results) 0
  }

  #--8. Is file name (0/1)
  if { $vtype == "File" } {
    lappend Info(results) 1
  } else {
    lappend Info(results) 0
  }

  #--9. Is procedure name (0/1)
  if { $vtype == "Procedure" } {
    lappend Info(results) 1
  } else {
    lappend Info(results) 0
  }

}


# Check from array SKWD if value type is given for
# a keyword in a solver section
# Section and keyword are given from cpp-side via the
# Info(arguments) and result is returned using
# the global Info(results) transfer variable
# 
proc Interface::getSolverKeywordTypeGiven {} {
  global Info SKWD

  set arg_list [split $Info(arguments) $Info(argSeparator)]

  set section [lindex $arg_list 0]
  set kwd [lindex $arg_list 1]

	set section [string tolower $section]
	set kwd [string tolower $kwd]

	set sec ""

	switch $section {
 		"simulation" {set sec "SI"}
		"constants" {set sec "CN"}
		"body force" {set sec "BF"}
		"initial condition" {set sec "IC"}
		"boundary condition" {set sec "BC"}
		"equation" {set sec "EQ"}
		"material" {set sec "MT"}
		"solver" {set sec "SL"}
		"body" {set sec "BO"}
		"boundary" {set sec "BA"}
	}

	if { $sec == "" || ![info exists SKWD($sec,$kwd) ] } {
		set Info(results) 0
	} else {
		set Info(results) $SKWD($sec,$kwd)
	}
}


proc Interface::getUseVariableNameInEquationName {} {
  global Info

  set eq_index $Info(arguments)
  set Info(results) [Equation::useVariableNameInEquationName $eq_index]
}


proc Interface::initEquationedPanels {} {

   set reset 0

   StdPanelInit::initPanelData InitialCondition $reset
   StdPanelInit::initPanelData BoundaryCondition $reset
   StdPanelInit::initPanelData BodyForce $reset
   StdPanelInit::initPanelData Material $reset
}


# Inform gui about a Matc-file which was read by the model
#
proc Interface::matcFileWasRead {} {
  global Info

  set fname $Info(arguments);
  
  # Store filename in info-list
  Util::addToFileList Info(matcFiles) $fname

  set msg "***MATC file loaded: $fname"

  Message::showMessage $msg $Info(remMsgColor) 
}


# Inform gui about a color-file which was read by the model
#
proc Interface::colorFileWasRead {} {
  global Info

  set fname $Info(arguments);
  
  # Store filename in info-list
  Util::addToFileList Info(colorFiles) $fname

  set msg "***Color file loaded: $fname"

  Message::showMessage $msg $Info(remMsgColor) 
}


# For debugging!
#
proc Interface::printFlags {} {
  global ModelFlags

  set flags {
    GEOMETRY_TYPE_CAD GEOMETRY_TYPE_MESH
    DRAW_SOURCE_CAD DRAW_SOURCE_MESH
    DRAW_TARGET_BODIES DRAW_TARGET_SURFACES 
  }

  foreach flag $flags {
    set var "ModelFlags($flag"
    append var ")"
    set var [set $var]
    Message::showMessage "$flag = $var"
  }

}


proc Interface::resetObjectTable {} {

  MenuExec::resetObjectTable
}


proc Interface::resetObjectTableByType {} {
  global Info

  set type [string toupper $Info(arguments)]

  MenuExec::resetObjectTableByType $type
}


### Update FLAGS ###

# If there was a change data related to the target
#
proc Interface::setNeedsUpdate {} {
  global Model Info

  set target $Info(arguments)
#print "$target needs update!"
  set Model($target,needsUpdate) 1
}


# If target data was updated
#
proc Interface::setWasUpdated {} {
  global Model Info

  set target $Info(arguments)

#print "$target was updated!"
  set Model($target,needsUpdate) 0
}


# If a mesh has been generated for the model
#
proc Interface::setModelHasElmerMesh {} {
  global Model

  set Model(Mesh,exists) 1
}


# If any (applied) mesh parameter set
#
proc Interface::setModelHasMeshParameter {} {
  global Model

  set Model(Mesh,hasParameter) 1
}


# If model has any geometry related Matc-definitions
#
proc Interface::setModelHasMatcDefinitions {} {
  global Model

  set Model(hasMatcDefinitions) 1

  MenuBuild::configureMenuButtonOption Edit updateCadGeometry state normal
  MenuBuild::configureMenuButtonOption Edit dropMatcDefinitions state normal
}


# If Mesh was edited
#
proc Interface::setMeshEdited {} {
  global Model

  set Model(Mesh,edited) 1
}


# If (internal or external) Mesh exists
#
proc Interface::setMeshExists {} {
  global Model

  set Model(Mesh,exists) 1
}


proc Interface::setMeshInputUnit {} {
  global Info Model

  set unit [lindex $Info(arguments) 0]
  
  if { $unit > 0 } {
    set Model(meshInputUnit) $unit
  } else {
    set Model(meshInputUnit) 1.0
  }
}


##### MISCELLANOUS DATA UPDATE

proc Interface::saveModelPropertyData {} {
  ModelProperty::panelSave
}


proc Interface::selectBody {} {
  global Info BodyDisplay

  set sp $Info(fieldSeparator)

  set bd1_id [lindex $Info(arguments) 0]
  set lr1_id [lindex $Info(arguments) 1]
  set bd2_id [lindex $Info(arguments) 2]
  set lr2_id [lindex $Info(arguments) 3]
  set mode [lindex $Info(arguments) 4]

  ListBox::selectBody $bd1_id $lr1_id $bd2_id $lr2_id $mode

  # NOTE: For bodies we do not actually know the selection mode
  # but the client can toggle when handling this call!
  #
  if { [winfo exist $BodyDisplay(winName)] } {
    BodyDisplay::setSelectionMode $bd1_id
    BodyDisplay::setSelectionMode $bd2_id
  }
}


# Set element selected in the list box
#
proc Interface::selectBoundary {} {
  global Info

  set bndr_id [lindex $Info(arguments) 0]
  set bd1_id [lindex $Info(arguments) 1]
  set lr1_id [lindex $Info(arguments) 2]
  set bd2_id [lindex $Info(arguments) 3]
  set lr2_id [lindex $Info(arguments) 4]
  set extend  [lindex $Info(arguments) 5]

  ListBox::selectBoundary $bndr_id $bd1_id $lr1_id $bd2_id $lr2_id $extend
}


proc Interface::setBoundarySelectionMode {} {
  global Info BoundaryDisplay VertexDisplay

  set bndr_id [lindex $Info(arguments) 0]
  set mode    [lindex $Info(arguments) 1]
  set do_upd  [lindex $Info(arguments) 2]

  Object::setBoundarySelectionMode $bndr_id $mode

  set tp [Object::getType $bndr_id]

  if { ($tp == "F" || $tp == "FG" || $tp == "E" || $tp == "EG") &&
       [winfo exist $BoundaryDisplay(winName)]
    } {
    BoundaryDisplay::setSelectionMode $bndr_id $mode $do_upd

  } elseif { ($tp == "V" || $tp == "VG") &&
             [winfo exist $VertexDisplay(winName)]
    } {
    VertexDisplay::setSelectionMode $bndr_id $mode $do_upd
  }
}


proc Interface::selectBoundary2 {} {
   global Info

#MSG "SelectElements received: $Info(arguments)"
}


proc Interface::selectBoundaries {} {
   global Info

#MSG "SelectBoundaries received: $Info(arguments)"
}


proc Interface::markSelectedBoundaries {} {

  ListBox::markSelectedBoundaries
}


proc Interface::rendererReset {} {

  MenuExec::rendererReset
}


proc Interface::setCurrentTimestamp {} {

  Util::setCurrentTimestamp
}


proc Interface::setCurrentMeshH {} {
  global Info Model

  set mesh_h $Info(arguments)
  set index $Model(curentMeshIndex)

  set $Model(meshHs) [lreplace $Model(meshHs) $index $index $mesh_h]
}


proc Interface::setCurrentMeshF {} {
  global Info Model 

  set mesh_f $Info(arguments)
  set index $Model(curentMeshIndex)

  set $Model(meshFs) [lreplace $Model(meshFs) $index $index $mesh_f]
}


proc Interface::setInitialMeshH {} {
  global Info MeshDefine

  set MeshDefine(MESH_H,default) $Info(arguments)

  set MeshDefine(MESH_F,default) 1.0
}


proc Interface::setInitialState {} {
  # Init timestep data 
  Panel::initFields Timestep
  Timestep::applyTimestepType
}


# EMF_FILE name is set to be the same as
# inModelFile (when it is checked in cpp-environment)
#
proc Interface::setModelFilePath {} {
  global Model

  if {$Model(inModelFile) != ""} {

    # By default in/out names the same (of course :-)
    set Model(EMF_PATH_IN) $Model(inModelFile)

    # Check that model file name extension is .emf
    set mfn [file root $Model(inModelFile)]
    append mfn ".emf"

    set Model(EMF_PATH) $mfn

    set Model(EMF_FILE) [file tail $Model(EMF_PATH)]
  }
}


# Activated form cpp if change in Emissivity
# boundary condition has been noticed
# NOT called currently MVe 05.11.97
#
proc Interface::setParameterFieldValueState {} {
  global Info Model

  set arg_list [split $Info(arguments) $Info(argSeparator)]
  set parameter_id  [lindex $arg_list 0]
  set field_name [lindex $arg_list 1]
  set has_value  [lindex $arg_list 2]
  set value_is_changed [lindex $arg_list 3]

  if { $field_name == "EMISSIVITY" } {
    set Model(Emissivity,hasValue) \
      [$Model(Emissivity,hasValue) || $has_value]

    set Model(Emissivity,valueHasChanged) \
      [$Model(Emissivity,valueHasChanged) || $value_has_changed]

    if { $Model(Emissivity,valueHasChanged) } {
      set Model(GebhardtFactors,needsUpdate,data)
    }
  }

#print "Gui: Emissivity has changed"
}


proc Interface::setRotatePriority {} {
  global Info

  set axis $Info(arguments)

  MenuExec::setRotatePriority $axis
}


# Set body force field in the status area
proc Interface::setStatusBodyForces {} {
  global Model Status

  set fld "bodyForces"
  set count $Model(nofBodiesWithBodyForce)

  if { $Status(useAbsoluteBase) } {
    set base $Model(nofBodies)

  } else {
    set base $Model(nofBodiesWithEquation)
  }

  if { $count == 0} {
    set st "None"
  } else {
    set st "$count / $base"
  }

  set Status($fld) $st

  Widget::configure1 $Status($fld,entryWidget) -fg "Black"
}


# Set inner boundary condtion field in the status area
proc Interface::setStatusInnerBoundaryConditions {} {
  global Info Model ObjectTable Status

  set fld "innerBoundaryConditions"
  set count $Model(nofInnerBoundariesWithCondition)

  set missing_param 0
  set ib_count 0
   
  if { $Model(GEOMETRY_DIMENSION) =="3D" } {
    set btp "F"
  } else {
    set btp "E"
  }

  foreach id $ObjectTable(ids) {
    
    if { $ObjectTable($id,tp) != $btp } {
      continue
    }

    set pr_id $ObjectTable($id,prId)

    # An inner boundary
    if { $ObjectTable($pr_id,tp) == "BP" } {

      set bd1_id $ObjectTable($pr_id,pr1Id)
      set bd2_id $ObjectTable($pr_id,pr2Id)

      # If both bodies have an equation but the boundary
      # does not have a condition
      if { $ObjectTable($bd1_id,eq) != $Info(NO_INDEX) &&
           $ObjectTable($bd2_id,eq) != $Info(NO_INDEX)
         } {
      
        incr ib_count

        if { $ObjectTable($id,bc) == $Info(NO_INDEX) } {
          set missing_param 1
        }
      }
    }
  }

  if { $Status(useAbsoluteBase) } {
    set base $Model(nofInnerBoundaries)

  } else {
    set base $ib_count
  }

  if { $count == 0} {
    set st "None"
  } else {
    set st "$count / $base"
  }

  set Status($fld) $st

  # NOTE: We do not mark missing parameter definitons
  # for inner boundaries, it is a bit hard to know
  # when it really is missing!
  Widget::configure1 $Status($fld,entryWidget) -fg "Black"
}


# Set outer boundary condtion field in the status area
#
proc Interface::setStatusOuterBoundaryConditions {} {
  global Info Model ObjectTable Status

  set fld "outerBoundaryConditions"
  set count $Model(nofOuterBoundariesWithCondition)

  set missing_param 0
  set ob_count 0
   
  if { $Model(GEOMETRY_DIMENSION) =="3D" } {
    set btp "F"
  } else {
    set btp "E"
  }

  foreach id $ObjectTable(ids) {
    
    if { $ObjectTable($id,tp) != $btp } {
      continue
    }

    set pr_id $ObjectTable($id,prId)

    # An outer boundary
    if { $ObjectTable($pr_id,tp) == "B" } {
      
      # If the body has an equation but boundary
      # does not have a condition
      #
      if { $ObjectTable($pr_id,eq) != $Info(NO_INDEX) } {
  
        incr ob_count

        if { $ObjectTable($id,bc) == $Info(NO_INDEX) } {
          set missing_param 1
        }
      }
    }
  }

  if { $Status(useAbsoluteBase) } {
    set base $Model(nofOuterBoundaries)

  } else {
    set base $ob_count
  }

  if { $count == 0} {
    set st "None"
  } else {
    set st "$count / $base"
  }

  set Status($fld) $st

  if {$missing_param} {
    Widget::configure1 $Status($fld,entryWidget) -fg "Red"
  } else {
    Widget::configure1 $Status($fld,entryWidget) -fg "Black"
  }

}


# Set equation field in the status area
#
proc Interface::setStatusEquations {} {
  global Equation Info Model Status

  set pmask $Equation(problemMask)

  set fld "equations"
  set count $Model(nofBodiesWithEquation)

  set base $Model(nofBodies)

  #-No equations attached
  if { $count == 0 } {
    set st "None"

  #-Equations attached
  } else {
    set st "$count / $base"
    append st "    "

    foreach ename $Equation(allFields) {

      if { ![Equation::isEquation $ename] } {
        continue
      }

      set emask [DataField::getFieldProperty Equation $ename EquationMask]

      if { [Util::patternMatchesString $emask $pmask 0] } {
        set sl_nm [DataField::getFieldProperty Equation $ename StatusLineName]
        append st "$sl_nm  "
      }
    }
  }

  set Status($fld) $st

  if { $count <= 0 } {
    Widget::configure1 $Status($fld,entryWidget) -fg "Red"

  } else {
    Widget::configure1 $Status($fld,entryWidget) -fg "Black"
  }

}


# Set initialy condtion field in the status area
#
proc Interface::setStatusInitialConditions {} {
  global Info Model ObjectTable Status Timestep

  set fld initialConditions
  set count $Model(nofBodiesWithInitialCondition)

  set missing_param 0

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }
    
    if { $ObjectTable($id,eq) != $Info(NO_INDEX) &&
         $ObjectTable($id,ic) == $Info(NO_INDEX)
       } {
      set missing_param 1
      break
    }
  }

  if { $Status(useAbsoluteBase) } {
    set base $Model(nofBodies)
  } else {
    set base $Model(nofBodiesWithEquation)
  }

  if { $count == 0} {
    set st "None"
  } else {
    set st "$count / $base"
  }

  set Status($fld) $st

  if { $missing_param &&
       $Timestep(SIMULATION_TYPE) != "Steady State"
     } {
    Widget::configure1 $Status($fld,entryWidget) -fg "Red"
  } else {
    Widget::configure1 $Status($fld,entryWidget) -fg "Black"
  }

}


# Set material field in the status area
proc Interface::setStatusMaterials {} {
  global Info Model ObjectTable Status

  set fld "materials"
  set count $Model(nofBodiesWithMaterial)

  set missing_param 0

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    if { $ObjectTable($id,eq) != $Info(NO_INDEX) &&
         $ObjectTable($id,mt) == $Info(NO_INDEX)
       } {
      set missing_param 1
      break
    }
  }

  if { $Status(useAbsoluteBase) } {
    set base $Model(nofBodies)
  } else {
    set base $Model(nofBodiesWithEquation)
  }

  if { $count == 0} {
    set st "None"
  } else {
    set st "$count / $base"
  }

  set Status($fld) $st

  if {$missing_param} {
    Widget::configure1 $Status($fld,entryWidget) -fg "Red"
  } else {
    Widget::configure1 $Status($fld,entryWidget) -fg "Black"
  }

}


# Set nof active meshes field in the status area
#
proc Interface::setStatusMeshes {} {
  global Model ModelFlags Status

  set fld "meshes"

  set nof_meshes [llength $Model(activeMeshNames)]

  if { $nof_meshes == 0} {
    set st "None"

  } else {
    set st $nof_meshes
  }

  set Status($fld) $st

  if { $Model(nofBodiesWithEquation) > 0 &&
       $ModelFlags(GEOMETRY_TYPE_CAD) && $nof_meshes == 0
     } {
    Widget::configure1 $Status($fld,entryWidget) -fg "Red"
  } else {
    Widget::configure1 $Status($fld,entryWidget) -fg "Black"
  }

}


# Set timestep field in the status area
#
proc Interface::setStatusTimestamps {} {
  global Model Status Timestep

  set fld "timesteps"

  if { $Model(nofTimesteps) == 0 } {
    set Status($fld) "None"

  } elseif { $Timestep(SIMULATION_TYPE) == "Steady State" } {
      set Status($fld) "Max $Model(nofTimesteps)"

  } else {
      set Status($fld) $Model(nofTimesteps)
  }

  Widget::configure1 $Status($fld,entryWidget) -fg "Black"
}


# Receive command to display message in the 
# main message window
#
proc Interface::showMessage {} {
  global Info __GUIServerPort
  
  set sep $Info(argSeparator)
  set data [split $Info(arguments) $sep]

  set msg [lindex $data 0]

  # Try to unfold the if it happens to be a list
  # NOTE: Character like "\{\[" could induce an error
  # when applying the join-command, so we use catch!
  #
  if { [catch { set message [join $msg] } ] } {
    set message $msg
  }

  set message [string trim $message]

  set nof_line_feeds [lindex $data 1]
  set append [lindex $data 2]

  # Try to give some time for gui if we are
  # in a lengthy process!
  Util::updateGui
  Util::updateGui 0

  # Color for error/warning/note lines!
  if { 0 == [string compare -nocase -length 6 "***ERR" $message] } {
    set color $Info(errMsgColor)

  } elseif { 0 == [string compare -nocase -length 13 "***MATC ERROR" $message] } {
    set color $Info(errMsgColor)

  } elseif { 0 == [string compare -nocase -length 10 "MATC ERROR" $message] } {
    set color $Info(errMsgColor)

  } elseif { 0 == [string compare -nocase -length 7 "***WARN" $message] } {
    set color $Info(wrnMsgColor)

  } elseif { 0 == [string compare -nocase -length 3 "***" $message] } {
    set color $Info(remMsgColor)

  } else {
    set color ""
  }

  Message::showMessage $message $color $nof_line_feeds $append 
}


proc Interface::variableNameGuiToSif {} {
  global Info
  
  set gui_name $Info(arguments)
  set Info(results) [DataField::variableNameGuiToSif $gui_name]
}


proc Interface::variableNameSifToGui {} {
  global Info

  set sif_name $Info(arguments)
  set Info(results) [DataField::variableNameSifToGui $sif_name]
}




# ================================
# ================================
# CONVERSION OF OLDER VERSION DATA 
# ================================
# ================================


# ============
# Pre Converse
# ============

# Possible pre converse of older version data
#
proc Interface::preConverseInputData {globArray} {
  global Info
  upvar #0 $globArray theArray

  set theArray(doConversion) 0

  if { $Info(FRONT_INPUT_VERSION_NBR) < 0 } {
    return 
  }

  Interface::preConverseInputDataTo6 $globArray
  Interface::preConverseInputDataTo7 $globArray
  Interface::preConverseInputDataTo9 $globArray
}


# Pre-converse to version 6
#
proc Interface::preConverseInputDataTo6 {globArray} {
  global Info Common
  upvar #0 $globArray theArray

  if { $Info(FRONT_INPUT_VERSION_NBR) >= 6 } {
    return 
  }

  # Only Equation data needs pre conversion
  # =======================================
  if { $globArray != $Common(panelArray,EQ) } {
    return 
  }

  global Equation
  set modified 0

  foreach pid $Equation(ids) {
    
    # Converse DIFFUSION_VARIABLES --> AD_VAR
    #
    set data1 $Equation($pid,data)
    set data2 [string map {">DIFFUSION_VARIABLES=" ">AD_VAR="} $data1]

    set Equation($pid,data) $data2

    if { $data1 != $data2 } {
      set modified 1
      Message::showMessage "Equation parameter $pid modified!"
    }
  }

  if { $modified } {
    set theArray(doConversion) 1
    set Info(doConversion) 1
  }
}


# Pre-converse to version 7
#
proc Interface::preConverseInputDataTo7 {globArray} {
  global Info Common
  upvar #0 $globArray theArray

  if { $Info(FRONT_INPUT_VERSION_NBR) >= 7 } {
    return 
  }

  # Only BoundaryCondition,Equation and Solver data needs pre conversion
  # ==================================================
  if {  $globArray != $Common(panelArray,BC)  &&
        $globArray != $Common(panelArray,EQ)           &&
        $globArray != $Common(panelArray,SL)
     } {
    return 
  }

  if { $globArray == $Common(panelArray,BC) } {
  
    global BoundaryCondition
    set modified 0

    foreach pid $BoundaryCondition(ids) {
    
      # Converse NORMAL_TANGENTIAL_DISPLACEMENT --> NORMAL-TANGENTIAL_DISPLACEMENT
      # Converse NORMAL_TANGENTIAL_VELOCITY --> NORMAL-TANGENTIAL_VELOCITY
      #
      set data1 $BoundaryCondition($pid,data)
      set data2 [string map {"NORMAL_TANGENTIAL_DISPLACEMENT=" "NORMAL-TANGENTIAL_DISPLACEMENT="} $data1]
      set data2 [string map {"NORMAL_TANGENTIAL_VELOCITY=" "NORMAL-TANGENTIAL_VELOCITY="} $data2]

      set BoundaryCondition($pid,data) $data2

      if { $data1 != $data2 } {
        set modified 1
        Message::showMessage "Boundary Condition parameter $pid modified!"
      }
    }
  }

  if { $globArray == "Equation" } {
  
    global Equation
    set modified 0

    foreach pid $Equation(ids) {
    
      # Converse NAVIER_STOKES --> NAVIER-STOKES
      # Converse AD_VARS --> ADVECTION_DIFFUSION_EQUATION_vars
      #
      set data1 $Equation($pid,data)
      set data2 [string map {"NAVIER_STOKES=" "NAVIER-STOKES="} $data1]
      set data2 [string map {"AD_VARS=" "ADVECTION_DIFFUSION_EQUATION_vars="} $data2]

      set Equation($pid,data) $data2

      if { $data1 != $data2 } {
        set modified 1
        Message::showMessage "Equation parameter $pid modified!"
      }
    }
  }

  if { $globArray == $Common(panelArray,SL) } {
  
    global Solver
    set modified 0

    foreach pid $Solver(ids) {
    
      # Converse field name EQUATION_NAME --> EQUATION
      # Converse field name SOLVER_MESH --> MESH
      # Converse field value KE_Turbulence --> KE Turbulence (in EQUATION_NAME)
      # Converse field value Heat_Equation --> Heat Equation (in EQUATION_NAME)
      # Converse field value Stress_Analysis --> Stress Analysis (in EQUATION_NAME)
      # Converse field value Advection_Diffusion_Equation --> Advection Diffusion Equation (in EQUATION_NAME)
      #
      set data1 $Solver($pid,data)
      set data2 [string map {"EQUATION_NAME=" "EQUATION="} $data1]
      set data2 [string map {"SOLVER_MESH=" "MESH="} $data2]
      set data2 [string map {"KE_Turbulence" "KE Turbulence"} $data2]
      set data2 [string map {"Heat_Equation" "Heat Equation"} $data2]
      set data2 [string map {"Stress_Analysis" "Stress Analysis"} $data2]
      set data2 [string map {"Advection_Diffusion_Equation" "Advection Diffusion Equation"} $data2]

      set Solver($pid,data) $data2

      if { $data1 != $data2 } {
        set modified 1
        Message::showMessage "Solver parameter $pid modified!"
      }
    }
  }

  if { $modified } {
    set theArray(doConversion) 1
    set Info(doConversion) 1
  }
}


# Pre-converse to version 9
#
proc Interface::preConverseInputDataTo9 {globArray} {
  global Info Common
  upvar #0 $globArray theArray

  if { $Info(FRONT_INPUT_VERSION_NBR) >= 9 } {
    return 
  }

  # Only Solver data needs pre conversion
  # ==================================================
  if {  $globArray != $Common(panelArray,SL) } {
    return 
  }

  global Solver
  set modified 0

  foreach pid $Solver(ids) {
  
    # Converse field name MULTIGRID_LEVELS --> MG_LEVELS
    #
    set data1 $Solver($pid,data)
    set data2 [string map {"MULTIGRID_LEVELS=" "MG_LEVELS="} $data1]

    set Solver($pid,data) $data2

    if { $data1 != $data2 } {
      set modified 1
      Message::showMessage "Solver parameter $pid modified!"
    }
  }

  if { $modified } {
    set theArray(doConversion) 1
    set Info(doConversion) 1
  }
}



# =============
# Post Converse
# =============

# Possible post converse of older version data
#
proc Interface::postConverseInputData {globArray} {
  global Info
  upvar #0 $globArray theArray

  if { $Info(FRONT_INPUT_VERSION_NBR) < 0 } {
    return 
  }

  if { !$theArray(doConversion) } {
    return
  }

  Interface::postConverseInputDataTo7 $globArray
}


# Post-converse to version 7
#
proc Interface::postConverseInputDataTo7 {globArray} {
  global Info Common
  upvar #0 $globArray theArray

  if { $Info(FRONT_INPUT_VERSION_NBR) >= 7 } {
    return 
  }

  # Only Equation data needs post conversion
  # ========================================
  if { $globArray != $Common(panelArray,EQ) } {
    return 
  }

  global Equation

  foreach pid $Equation(ids) {

    Panel::resetFields Equation

    StdPanelInit::initFieldValues Equation

    # Create the updated equation paramter line
    # -----------------------------------------
    DataField::formDataFields Equation $Equation($pid,data)

    Equation(NAVIER-STOKES,fillProc)
    Equation(HEAT_EQUATION,fillProc)
    Equation(STRESS_ANALYSIS,fillProc)
    Equation(ADVECTION_DIFFUSION_EQUATION,fillProc)

    set Equation($pid,data) [DataField::formParameterLine Equation]

    #Panel::unsetFields $arr
    Panel::unsetAllFieldData $globArray
  }

  # Store updated data into model (Note: this does not save the model file, the user
  # has to save it actively!)
  #
  Util::cpp_exec readConvertedEquationData
}


# end ecif_tk_procsInterface.tcl
# ******************************
