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
#Module:    ecif_tk_procsMenuExec.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Menu handling procedures
#
#************************************************************************


#####################################################################
#
#                   Initialization procs
#
#####################################################################

# Apply user setting for directories and 
# auto OPEN model file if model-name was given!!!
#
proc MenuExec::applySettings {} {
  global Info Model ModelProperty UserSetting

  # Copy possible command line values
  foreach vn {MODEL_DIRECTORY MODEL_NAME PROBLEM_NAME} {
    if { $Info($vn) != "" } {
      set ModelProperty($vn) $Info($vn)
      set Info($vn) ""
    }
  }

  #--Set working directory
  if { $Info(workingDirectory) == "" } {
    set Info(workingDirectory) [MenuExec::getCurrentDirectory]
  }
  
  #--Set other default directoies
  #
  if { $Model(CAD_OPEN_DIR) == "" } { 
    set Model(CAD_OPEN_DIR) $Info(workingDirectory)
  }

  if { $Model(MESH_OPEN_DIR) == "" } { 
    set Model(MESH_OPEN_DIR) $Info(workingDirectory)
  }

  if { $Model(EMF_OPEN_DIR) == "" } { 
    set Model(EMF_OPEN_DIR) $Info(workingDirectory)
  }


  # Model directory
  # ===============
  set rc  [ModelProperty::setDefaultDirectoryValue MODEL_DIRECTORY]

  # If not ok
  if { $rc != 1 } {

    # If invalid value
    if { $rc == -1 } {
      set dir $UserSetting(DEFAULT_MODEL_DIRECTORY)
      Message::showMessage "NOTE: Invalid model directory value: $dir"
    }

    # Put current directory as default
    ModelProperty::setDefaultDirectoryValue MODEL_DIRECTORY $Info(workingDirectory)
  }

  ModelProperty::applyModelDirectoryValue

  # Include path
  # ============
  set dir $UserSetting(DEFAULT_INCLUDE_PATH)
  if { $dir != "" } {
    set rc [ModelProperty::setDefaultDirectoryValue INCLUDE_PATH $dir]
  }

  # Log directory
  # =============
  set dir $Info(defaultLogDirectoryName)

  if { $dir != "" } {
    set rc [ModelProperty::setDefaultDirectoryValue LOG_DIRECTORY $dir]
  }

  # Results directory
  # =================
  set dir $UserSetting(DEFAULT_RESULTS_DIRECTORY)
  if { $dir != ""  } {
    set rc [ModelProperty::setDefaultDirectoryValue RESULTS_DIRECTORY $dir]
  }

  # Cad files directory
  # ===================
  set dir $UserSetting(DEFAULT_CAD_FILES_DIRECTORY)
  if { $Model(CAD_OPEN_DIR) == "" &&  $dir != "" && [file isdirectory $dir] } {
    set Model(CAD_OPEN_DIR) $dir
  }

  # Mesh files directory
  # ====================
  set dir $UserSetting(DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY)
  if { $dir != "" && [file isdirectory $dir] } {

    if { $Model(MESH_OPEN_DIR) == "" } {
      set Model(MESH_OPEN_DIR) $dir
    }

    if { $Model(MESH_SAVE_DIR) == "" } {
      set Model(MESH_SAVE_DIR) $dir
    }
  }

  # Font sizes
  # ==========
  # Target widgets in order: Entry, Table Entry, and Message Area
  #
  if { $UserSetting(FONT_SIZES) != "" } {
    
    foreach sz $UserSetting(FONT_SIZES) \
            trg {fontSizeEntry fontSizeTable fontSizeMsg} {

      if { $sz == "" || $trg == "" } {
        break
      }

      set Info($trg) $sz
    }
  }

  # Other data
  # ==========

  #--Set and check case directories
  if { $Model(CASE_NAME) != "" } {
    set force 0
    set show_msg 1
    MenuExec::applySettingsToCaseDirectories $force $show_msg
  }

  set model_dir $ModelProperty(MODEL_DIRECTORY,absolute)

  #--Default model open and save dirs
  set Model(EMF_OPEN_DIR) $model_dir
  set Model(EMF_SAVE_DIR) $model_dir

  #--Open model file, if name given !!!
  #--------------------------------
  #
  if { $ModelProperty(MODEL_NAME) != "" } {

    set emf_path [file join $model_dir $ModelProperty(MODEL_NAME).emf]
  
    if { [file exists $emf_path] } {
      Message::showMessage "Opening model file: $emf_path"
      #MenuExec::presetModelData
      MenuExec::loadNewModel $emf_path

    } else {
      Message::showMessage "NOTE: cannot open model file: $emf_path"
    }
  }

} ;# End applySettings


proc MenuExec::checkUserLevel {} {
  global Info

  set Info(userLevel,prev) $Info(userLevel)

  if { $Info(userLevel,new) == 0 } {
    set Info(userLevel) $Info(userLevel,prev)
  } else {
    set Info(userLevel) $Info(userLevel,new)
  }

  MenuBuild::applyUserLevel
  UserSetting::applyUserLevel
}


# Reset model data before loading
# a new model
#
proc MenuExec::presetModelData {} {
  global Common Coordinate Info Model ModelProperty ObjectTable
  global Solver Equation

  set Model(CASE_NAME) ""
  set Model(GEOMETRY_DIMENSION)   ""
  set Model(SIMULATION_DIMENSION) ""
  set Coordinate(COORDINATE_MAPPING) {1 2 3}

  set Model(EMF_PATH_IN) ""
  set Model(EMF_PATH) ""
  set Model(EMF_FILE) ""
  set Model(inModelFile) ""
  set Model(cadPath) ""

  set Model(sifModelIncludeFile) ""
  set Model(sifSimulationIncludeFile) ""

  set Model(hasMatcDefinitions) 0

  set Model(meshNames) ""
  set Model(nofActiveMeshes) 0
  set Model(activeMeshIndices) ""
  set Model(activeMeshNames) ""

  set Model(meshBgMeshFiles) ""
  set Model(meshBgMeshActives) ""
  set Model(meshBgMeshControls) ""
  set Model(currentMeshBgMeshFile) ""
  set Model(currentMeshBgMeshControlFile) ""

  MenuExec::resetObjectTable

  # Reset panel arries
  # NOTE: Skip UserSetting!
  # =======================     
  # (Its values are kept during the whole session!)
  foreach arr $Common(allPanelArries) {

    if { $arr == "UserSetting" } {
      continue
    }

    MenuExec::resetArrayData $arr
  }


  # These arries need a bit special handling, because
  # they do not have id-parameters!
  Panel::initFields ModelProperty "" 1
  Panel::initFields Variables

  #Panel::initFields EquationVariable
  #EquationVariable::formActiveParameter
  #EquationVariable::updateAllVariableNames

  set force 1
  set show_msg 0
  ModelProperty::initCaseDirectoryVariables

  #---General model level flags
  set Model(meshNames) ""
  set Model(meshHs) ""
  set Model(meshFs) ""
  set Model(meshBgMeshFiles) ""
  set Model(meshBgMeshActives) ""
  set Model(meshBgMeshControls) ""

  set Model(currentMeshIndex) $Info(NO_INDEX)
  set Model(currentMeshName) ""
  set Model(currentBgMeshFile) ""
  set Model(currentBgMeshControlFile) ""

  set Model(hasDiffuseGrayRadiation) 0
  set Model(nofMeshZeroVelocityElements) 0
  set Model(nofBodiesWithEquation) 0

  set Model(Database,needsUpdate) 0
  set Model(Front,needsUpdate) 0
  set Model(GebhardtFactors,needsUpdate,data) 0
  set Model(GebhardtFactors,needsUpdate,mesh) 0
  set Model(GebhardtFactors,verifyFilename) 0
  set Model(Mesh,exists) 0
  set Model(Mesh,generating) 0
  set Model(Mesh,edited) 0
  set Model(Mesh,hasParameter) 0
  set Model(Mesh,needsUpdate) 0
  set Model(Mesh,added) 0
  set Model(MeshFromDB,needsUpdate) 0
  set Model(Results,activated) 0
  set Model(Solver,inputFileExists) 0
  set Model(Solver,needsUpdate) 0
  set Model(Solver,askZeroVelocityElements) 1
  set Model(Solver,askNonExistingMesh) 1
  set Model(Viewfactors,needsUpdate,data) 0
  set Model(Viewfactors,needsUpdate,mesh) 0
  set Model(Viewfactors,verifyFilename) 0
  set Model(removedFieldsList) ""

  set Equation(totalMask) ""
  set Equation(problemMask) ""
  set Equation(activeIndices) ""
  set Equation(updateValuesArea) 0
  Equation::initMasks
    
  MenuExec::resetModelFlags
  MenuExec::initStatusFields

  #-Apply user setting values
  MenuExec::applySettings

  #-Read possible color tables (via cpp)
  MenuExec::readColorFile "front"
  MenuExec::readColorFile "elmer"
  MenuExec::readColorFile "user"
  MenuExec::readColorFile "model"

  #-Read possible matc variable definitions (via cpp)
  MenuExec::readMatcFile "elmer"
  MenuExec::readMatcFile "user"
  MenuExec::readMatcFile "model"
}


proc MenuExec::resetObjectTable {} {
  global ObjectTable

  set ObjectTable(count) 0

  if { [info exist ObjectTable(ids)] } {
    Util::unsetArrayIdVariables ObjectTable $ObjectTable(ids)
    unset ObjectTable(ids)
  }

  set ObjectTable(ids) ""
}


proc MenuExec::resetObjectTableByType {type} {
  global ObjectTable

  if { ![info exist ObjectTable(ids)] } {
	  return
  }

  set ids_to_keep ""

  set delete_count 0
  set ids_to_delete ""
  
  foreach id $ObjectTable(ids) {
    if { $ObjectTable($id,tp) == $type } {
      lappend ids_to_delete $id
      incr delete_count
    } else {
      lappend ids_to_keep $id
    }
  }

  if { $delete_count > 0 } {
    Util::unsetArrayIdVariables ObjectTable $ids_to_delete
  }

  set ObjectTable(ids) $ids_to_keep
}


proc MenuExec::resetArrayData { arr } {
  upvar #0 $arr theArr
  
  set theArr(panelSizeKnown) 0

  #Panel::unsetFields $arr
  Panel::unsetAllFieldData $arr

  if { [info exist theArr(ids)] } {
    Util::unsetArrayIdVariables $arr  $theArr(ids)
    unset theArr(ids)
  }

  if { [info exist theArr(ids,old)] } {
    unset theArr(ids,old)
  }

  Util::unsetArrayVariables $arr objectIndex
  Util::unsetArrayVariables $arr boundaryIndex
  Util::unsetArrayVariables $arr parameterIndex
}


#
proc MenuExec::resetUpdateFlags {} {
  global Model

  set Model(Front,needsUpdate) 0
  set Model(Database,needsUpdate) 0
  set Model(GebhardtFactors,needsUpdate,data) 0
  set Model(GebhardtFactors,needsUpdate,mesh) 0
  set Model(Viewfactors,needsUpdate,data) 0
  set Model(Viewfactors,needsUpdate,mesh) 0
}


# Update model data after loading
# a new model
#
proc MenuExec::updateModelData {} {
  global Model

  #MenuExec::resetUpdateFlags

  Solver::updateActiveMeshInfo

  set check_factors 0
  BoundaryCondition::updateModelStatus $check_factors

  Interface::setStatusMeshes
}


# Set status fields to initial values
#
proc MenuExec::initStatusFields {} {
  global Status

  set Status(equations) "None"
  set Status(timesteps) "None"
  set Status(materials) "None"
  set Status(initialConditions) "None"
  set Status(bodyForces) "None"
  set Status(innerBoundaryConditions) "None"
  set Status(outerBoundaryConditions) "None"
  set Status(meshes) "None"
}


# =========================
# EQUATION DEFINITONS STUFF
# =========================

# Load user defined equations definition file
proc MenuExec::loadDefinitionFile {} {
  global Common Model

  set types {
    { {All} {*} }
    { {Txt} {.txt} }
    { {Elmer def file} {.edf .Edf .EDF} }
  }

  set title "Load Definition File"

  # Use the latest directory
  set dir $Model(inDefinitionFileDir)

  #--Display dialog and get path
  set path [ MenuExec::openFile $title $types $dir ]

  if {$path == ""} {
    return
  }

  set Model(inDefinitionFile) $path
  set Model(inDefinitionFileDir) [file dirname $path]
  
  # NOTE: We do not read default def files or reset model data here!!!
  set read_default_def_file 0
  set read_in_def_file 1
  set preset_model_data 0

  FRONT_POST_INIT $read_default_def_file $read_in_def_file $preset_model_data

  # Mark panel size unknown (possible new fields!)
  #
  foreach arr $Common(allPanelArries) {
    upvar #0 $arr theArr
    set theArr(panelSizeKnown) 0
  }
} 


# Load user defined equations definition file
proc MenuExec::saveDefinitionFileAs {} {
  global Model

  set types {
    { {Elmer def file} {.edf .Edf .EDF} }
    { {Txt} {.txt} }
    { {All} {*} }
  }

  set title "Save Definition File"

  set default_path $Model(outDefinitionFile)

  #--Display dialog and get saved path
  set path [MenuExec::saveFile $title $types $default_path]

  if {$path == ""} {
    return
  }

  set append 0
  MenuExec::saveDefinitionFile $path $append
} 


# Saves current definitions in a file
#
proc MenuExec::saveDefinitionFile {path {append 0}} {
  global Model

  if { $path == "" } {
    return 0
  }

  if { $Model(allDefinitions) == "" } {
    return 0
  }

  set msg ""

  if { $append } {
    set mode "a+"
  } else {
    set mode "w"
  }

  if { [catch {set def_ch [open $path $mode]} msg] } {
    set msg [string map {"\"" " "} $msg]
    Message::showMessage "NOTE: cannot save definition file: $msg"
    return 0
  }

  foreach line $Model(allDefinitions) {
    puts $def_ch $line
  }

  close $def_ch

  return 1
}


# Clear current definitions
#
proc MenuExec::clearCurrentDefinitions { {verify 1}} {
  global Info Model

  set msg [list \
     "This command removes all definitions loaded by the user!\n\n" \
     "\nPress OK to remove definitions, otherwise press Cancel" \
     ]

  if { $verify &&
       ![Message::verifyAction $Info(powerUser) $msg]
     } {
    return
  }

  set Model(hasUserDefinitions) 0
  set Model(inDefinitionFile) ""
  set Model(allDefinitions) ""

  FRONT_SET_INITIAL_GENERIC_FIELD_PROPERTIES
  FRONT_SET_INITIAL_PANEL_FIELD_PROPERTIES

  # NOTE: We do not reset model data or read any new def file here!!!
  set read_default_def_file 0
  set read_in_def_file 0
  set preset_model_data 0

  FRONT_POST_INIT $read_default_def_file $read_in_def_file $preset_model_data
} 



#####################################################################
#
#                        CHECKING procs
#
#####################################################################

# Procedure checks what to do when
# model file is not saved
#
proc MenuExec::askModelOk {process_name {automatic_save 0} } {
  global Info Model

  if { $Model(Front,needsUpdate) == 0 } {
    return "ok"
  }

  if {$automatic_save} {
    set rc [saveModelFile]
    return $rc
  }

  set msg [list "Model file is not saved!\n\n" \
                 $Info(anywayOk) ]

  set Info(messageIcon) warning
  set result [ Message::dispOkCancelMessage $msg ]

  return $result
}


# Procedure checks what to do when mesh is missing when
# sif file is saved
#
proc MenuExec::askSolverOk { {path ""}} {
  global Model Info

  if { $Model(Mesh,exists) ||
       $Model(Mesh,generating) ||
       !$Model(Solver,askNonExistingMesh)
     } {
    return "ok"
  }
  
  # A warning of a non-existing mesh!
  #
  if { $path != "" } {
    set msg [list "Saving solver input file as:\n\n" \
                  "$path \n\n" \
                  "NOTE: Elmer mesh for the model does not exist\n" \
                  "and the problem cannot be possibly solved!\n\n" ]
  } else {
    set msg [list "NOTE: Elmer mesh for the model does not exist\n" \
                  "and the problem cannot be possibly solved!\n\n" ]
  }

  set Info(messageIcon) warning
  
  #--Solving
  #
  if { $Model(Solver,solving) } {
    set result [ Message::dispOkCancelMessage $msg ]

    # User cancelled sif-file saving
    #
    if { $result == "cancel" } {
      set Model(Solver,solving) 0
      return "cancel"
    }
  
  #--Just saving sif file
  #
  } else {
    Message::dispOkMessage $msg
  }

  # If user proceeds, we will not notify about this any more
  # for this model file!
  #
  set Model(Solver,askNonExistingMesh) 0

  return "ok"
}


# Checks that Cad geometry data is ok 
#
proc MenuExec::checkCadData {process_name} {
  global Info Model ModelFlags

  #-These processes are meaningless without
  # Cad geometry
  if { $process_name == $Info(Mesh2D,processName)  ||
       $process_name == $Info(Mesh3D,processName)
     } {

    if { !$ModelFlags(GEOMETRY_TYPE_CAD) } {
      return "cancel"
    }
  }

  return "ok"
}


# Check that model DB-directory is given
#
proc MenuExec::checkDBDirectory {} {
  global Info ModelProperty 

  if { $ModelProperty(MODEL_DIRECTORY) == "" } {
    set msg [list "Model directory not defined!\n" \
                  "Use command: Save model file as... in File menu!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  return "ok"
}


# Check that model DB-name is given
#
proc MenuExec::checkDBName { {name_needed 1} } {
  global ModelProperty Info

  if {$ModelProperty(MODEL_NAME) == ""} {

    #-Model must have name
    if { $name_needed } {
      set msg [list "NOTE: Model name is not defined, model cannot be saved!\n\n" \
                    $Info(setModelNameMsg) ]

      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return "cancel"

    #-Model name optional
    } else {
      set msg [list "NOTE: Model name is not defined, model cannot be saved!\n\n" \
                    "$Info(setModelNameMsg) \n\n" \
                    $Info(anywayOk) ]

      set Info(messageIcon) warning
      return [ Message::dispOkCancelMessage $msg ]
    }
  }

  return "ok"
}


# Checks if GebhardtFactors should be calculated (when running the Solver)
# and that the file exists (when running theSolver) and notifies if default file
# name is used (when running the GebhardtFactors)
#
proc MenuExec::checkGebhardtFactors {process_name} {
  global Equation Info Datafile Model ModelProperty

  #---If not an interesting process
  if { $process_name != $Info(Solver,processName) && $process_name != $Info(GebhardtFactors,processName) } {
    return "ok"
  }

  #---Proc Solver  checking
  if { $process_name == $Info(Solver,processName) &&
       $Model(hasDiffuseGrayRadiation) 
     } {
    
    # If the user has not turned of the checking for the session
    # or for this case for viewfactors check status against mesh
    #
    if { $Model(Viewfactors,needsUpdate,mesh) == 0 &&
         $Model(GebhardtFactors,needsUpdate,mesh) != -1
       } {

      # Check if mesh is more recent than the GebhardtFactors file
      set eidx [DataField::getFieldProperty Equation HEAT_EQUATION EquationIndex]
      set dc_mh [MenuExec::getSolverMeshHeaderFile $eidx 0]
      set dc_mn [MenuExec::getSolverMeshName $eidx 0]
      set gf_fn [MenuExec::getMeshRelatedFile $dc_mn GEBHARDT_FACTORS]

      if { $gf_fn == "" ||
           [ file mtime $dc_mh] > [file mtime $gf_fn]
         } {

        set Model(GebhardtFactors,needsUpdate,mesh) 1

        set msg [list "Gebhard factors should be updated (MESH is more recent)!\n\n" \
                      $Info(noMoreNotifying) ] 

        set Info(messageIcon) warning
        set result [Message::dispCancelYesNoMessage $msg]

        if { $result == "cancel" } {
          return "cancel"
        } elseif { $result == "yes" } {
          set Model(GebhardtFactors,needsUpdate,mesh) -1
          return "ok"
        }
      }
    }


    # If emissivity has been changed (checked elswhere!)
    # ==============================
    if { $Model(GebhardtFactors,needsUpdate,data) == 1 } {
      set msg [list "Gebhard factors should be updated (EMISSIVITY has been changed)!\n\n" \
                    $Info(noMoreNotifying) ] 
     
      set Info(messageIcon) warning

      set Info(messageIcon) warning
      set result [Message::dispCancelYesNoMessage $msg]

      if { $result == "cancel" } {
        return "cancel"
      } elseif { $result == "yes" } {
        set Model(GebhardtFactors,needsUpdate,data) -1
        return "ok"
      }
    }
  }

  #---Proc GebhardtFactors checking

  if { $Model(GebhardtFactors,verifyFilename) } { 

    set msg [list "Default filename $Datafile(GEBHARDT_FACTORS)\n" \
                  "will be used for Gebhardt factors!\n\n"             \ 
                  $Info(anywayOk) ]
  
    if { ![Message::verifyAction $Info(advancedUser) $msg] } {
      return "cancel"

    } else {
      set Model(GebhardtFactors,verifyFilename) 0
    }
  }

  #-Check that output file is writable
  set reserved_files ""

  set path [MenuExec::getMeshRelatedFile $Model(currentMeshName) GEBHARDT_FACTORS]

  if { ![Util::isWritableFile $path] } {
    append reserved_files "Gebhardt factors file"
  }

  if { $reserved_files != "" } {
    set msg [list "Cannot open the Gebhardt factors output file for writing!\n" \
                  "Enter a new file name in the Problem/Data files panel or wait until\n" \
                  "the process which is using the file, has finished!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  return "ok"
}


# Checks that mesh parameter data is ok 
# for mesh generator
#
proc MenuExec::checkMeshData {process_name} {
  global Info Model ModelProperty Status

  if { $process_name != $Info(Mesh2D,processName) &&
       $process_name != $Info(Mesh3D,processName)
     } {
    return "ok"
  }

  #---Check that mesh directory exists
  if {"ok" != [Util::checkAndCreateDirectory2B \
                  $ModelProperty(CURRENT_MESH_DIRECTORY,absolute) \
                  "Elmer mesh files"]
     } {
    return "cancel"
  }

  return "ok"
}


# Check that model outfile exist
#
proc MenuExec::checkModelFile {} {
  global Info Model

  if {$Model(EMF_FILE) == ""} {
    set msg [list "Model file not saved!\n" \
                  "Save model file: (File/Save model file)" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  return "ok"
}


proc MenuExec::checkProcessEnvironment {process_name_var} {
  global Info Model Datafile ModelProperty UserSetting
  upvar $process_name_var process_name

  #---Preprocess (meshing)
  if { $process_name == $Info(Mesh2D,processName) ||
       $process_name == $Info(Mesh3D,processName)
     } {
      set Model(Mesh,generating) 1 

  #---Solve
  } elseif { $process_name == $Info(Solver,processName) } {
    set Model(Solver,solving) 1

    set rc [MenuExec::askSolverOk]

    if { $rc != "ok" } {
      return "cancel"
    }
  }

  #---Postprocess
  #
  # NOTE: Postprocessor can be started although model
  # file is not available or stored!
  if { $process_name == $Info(Post,processName) } {
    return "ok"
  }

  if { $process_name == $Info(Results,processName) } {
    # If no model opened, start plain Postprocessor!
    #
    if { ![info exists Datafile(POST_FILE)] ||
         $Datafile(POST_FILE) == ""
       } {
      set process_name $Info(Post,processName)
      return "ok"
    }

    if { "ok" != [MenuExec::checkPostFile $process_name] } {
      return "cancel"
    }

    return "ok"
  }
   

  #---Check that DB-dir is available and model file exists
  #
  if { "ok" != [MenuExec::checkDBDirectory] } {
    return "cancel"
  }

  if { "ok" != [MenuExec::checkDBName 1] } {
    return "cancel"
  }

  set automatic_save 0
  #---Check that model file is saved
  #
  if { $process_name == $Info(GebhardtFactors,processName) ||
       $process_name == $Info(Mesh2D,processName)          ||
       $process_name == $Info(Mesh3D,processName)          ||
       $process_name == $Info(Solver,processName)          ||
       $process_name == $Info(Viewfactors,processName) } {
    set automatic_save $UserSetting(AUTO_SAVE_MODEL)

  } elseif {
    $process_name == $Info(Emf2Db,processName) } {
    set automatic_save 1

  } else {
    set automatic_save 0
  }

  if { "ok" != [ MenuExec::askModelOk $process_name $automatic_save ] } {
    return "cancel"
  }

  #---Check that mesh data is uptodate if process needs it
  #
  if { "ok" != [MenuExec::checkMeshData $process_name] } {
    return "cancel"
  }

  #---Check that solver related data is ok
  #
  if { "ok" != [MenuExec::checkSolverData $process_name] } {
    return "cancel"
  }

  if { "ok" != [MenuExec::checkViewfactors $process_name] } {
    return "cancel"
  }

  if { "ok" != [MenuExec::checkGebhardtFactors $process_name] } {
    return "cancel"
  }

  return "ok"
}


# Checks that ElmerPost file exists
#
proc MenuExec::checkPostFile {process_name} {
  global Info Model

  if { $process_name != $Info(Results,processName) } {
    return "ok"
  }

  set path [MenuExec::getMeshRelatedFile $Model(currentMeshName) POST_FILE]

  if { ![file exists $path] || ![file readable $path] } {
    set msg [list "Postprocessor file: $path is not available\n\n" \
                   $Info(continueOk) ]
    
    if { ![Message::verifyAction $Info(noviceUser) $msg cancel] } {
      return "cancel"
    }
  }

  return "ok"
}


# Checks that data is ok for solver
#
proc MenuExec::checkSolverData {process_name} {
  global Info Model ModelProperty Solver
  #--Proc Solver  checking

  if { $process_name != $Info(Solver,processName) } {
    return "ok"
  }

  #---Check that Mesh data exists
  set model_dir $ModelProperty(MODEL_DIRECTORY,absolute)

  if { ![file isdirectory $model_dir] } {

    set msg [list "Directory for the model data:\n" \
                  "$model_dir \n"    \
                  "does not exist, cannot continue!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  #-Check that all equations have an existing mesh
  set msg_info ""
  if { ![MenuExec::checkEquationMeshes msg_info] } {

    set msg [list "NOTE: The following equations do not have a mesh\n" \
                  "in the bodies they should be calculated:\n\n" \
                  "$msg_info \n\n" \
                  "Cannot continue!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  #---Check that Model data is ok

  # No equations!
  if { $Model(nofBodiesWithEquation) <= 0 } {
    set msg [list "No equations defined!\n" \
                  "Nothing to solve!\n\n"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  # If error in model (status != 0)
  Util::cpp_exec checkModelStatus
  if { $Model(status) != 0 } {

    set msg [list "Model is not completely defined!\n" \
                  "$Model(statusMessage)\n"      \
                  "Cannot possibly solve the problem!\n\n" \
                   $Info(anywayOk) ]

    if { ![Message::verifyAction $Info(noviceUser) $msg] } {
      return "cancel"
    }
  }

  #-Check that solver-output files are writable
  set reserved_files ""

  foreach sid $Solver(activeIds) { 

    set ename [DataField::getFieldValue Solver $sid EQUATION]

    set path [MenuExec::getSolverResultFile $sid OUTPUT_FILE]

    if { ![Util::isWritableFile $path] } {
      append reserved_files "Result file for: $ename  "
      break
    }

    set path [MenuExec::getSolverResultFile $sid POST_FILE]

    if { ![Util::isWritableFile $path] } {
      append reserved_files "Postprocessor file for: $ename "
      break
    }
  }

  # If even one of the files is reserved, recommend stop!
  #
  if { $reserved_files != "" } {
    set msg [list "The following Solver output file(s) could be unavailable for writing:\n\n" \
                  "$reserved_files\n\n" \
                  "Enter a new file name in the Problem/Data files panel or wait until\n" \
                  "the process which is using the files, has finished!" \
                   $Info(anywayOk) ]

    if { ![Message::verifyAction $Info(powerUser) $msg] } {
      return "cancel"
    }

  }


  # Check zero-velocity elements
  # ----------------------------
  Util::cpp_exec checkMeshCornerElements

  if { $Model(Solver,askZeroVelocityElements) &&
       $Model(nofMeshZeroVelocityElements) > 0
     } {

    set msg [list "Mesh contains zero-velocity element(s).\nThis may create singularity problems in Solver!\n\n" \
                   "Press YES to correct the mesh by splitting the elements\n\n" \
                   "Press NO to continue without correcting the mesh"]  

    set Info(messageIcon) warning
    set result [Message::dispCancelYesNoMessage $msg]

    if { $result == "cancel" } {  
      return "cancel"
    }

    if { $result == "yes" } {

      # Correct zero-velocity elements
      # ----------------------------
      Util::cpp_exec correctMeshZeroVelocityElements
    
      # Save modified mesh
      # ------------------
      set Model(Mesh,edited) 1
      MenuExec::saveModelFile
    }

  }
  
  #Finally, all ok?!
  return "ok"
}


# Check if a body is exlcuded from the given mesh
#
proc MenuExec::bodyExcludedFromMesh { id mesh_index } {
  global ObjectTable
  
  if { [info exists ObjectTable($id,tp)] ||
       $ObjectTable($id,tp) != "B"
     } {
    return 0
  }

  # If some of the body layers are in the mesh, then the body itself
  # is in the mesh!
  #
  foreach lr_id $ObjectTable($id,lrIds) {
    if { -1 == [lsearch $ObjectTable($lr_id,excldMshIndcs) $mesh_index] } {
      return 0
    }
  }

  return 1
}


proc MenuExec::checkEquationMeshes { {msg_var ""} } {
  global Equation Info ObjectTable

  if { $msg_var != "" } {
    upvar $msg_var msg
  }

  set msg ""
 
  foreach id $ObjectTable(ids) {

    # If not a body object, or the body does not have an equation
    #
    if { ![info exists ObjectTable($id,eq)] ||
         $ObjectTable($id,eq) == $Info(NO_INDEX)
       } {
      continue
    }
    
    set eq_id $ObjectTable($id,eq)
    set sids [Solver::findSolverIdsForEquation $eq_id]

    if { $sids == "" } {
      continue
    }

    foreach sid $sids {

      set mesh_index [Solver::getMeshIndex $sid]

      if { [MenuExec::bodyExcludedFromMesh $id $mesh_index] } {

        append mrow "Equation: "
        set eq_nm [DataField::getFieldValue Solver $sid EQUATION]
        append mrow [string map {"_" " "} $eq_nm]
        append mrow "   "

        set mrow "Body: "
        append mrow $ObjectTable($id,nm)

        append msg $mrow
        append msg "\n"
      }

    } ;# Each solver for equation

  } ; # Each object

  if { $msg != "" } {
    return 0

  } else {
    return 1
  }
}


#-Check that Solver input files exists
# in the model directory
#
proc MenuExec::checkSolverInputFile {} {
  global Datafile Model

  set Model(Solver,inputFileExists) 0

  if { [file exists $Datafile(SOLVER_INPUT_FILE)] } {
    set Model(Solver,inputFileExists) 1
  }
}


# Checks if Viewfactors should be calculated (when running the Solver)
# and that the file exists (when running theSolver) and notifies if default file
# name is used (when running the Viewfactors)
#
proc MenuExec::checkViewfactors {process_name} {
  global Equation Info Model ModelProperty

#print "process=$process_name"

  #---If not an interesting process
  if { $process_name != $Info(Solver,processName) &&
       $process_name != $Info(Viewfactors,processName) } {
    return "ok"
  }

  #---Proc Solver  checking
  if { $process_name == $Info(Solver,processName) &&
       $Model(hasDiffuseGrayRadiation)
     } {

    # Check if mesh is more recent than the Viewfactors file
    # ======================================================
    # Check only if the checking is not turned of by the user
    if { $Model(Viewfactors,needsUpdate,mesh) != -1 } {

      set eidx [DataField::getFieldProperty Equation HEAT_EQUATION EquationIndex]
      set dc_mh [MenuExec::getSolverMeshHeaderFile $eidx 0]
      set dc_mn [MenuExec::getSolverMeshName $eidx 0]
      set vf_fn [MenuExec::getMeshRelatedFile $dc_mn VIEW_FACTORS]

      if { $vf_fn == "" ||
           [ file mtime $dc_mh] > [file mtime $vf_fn ]
         } {

        set Model(Viewfactors,needsUpdate,mesh) 1

        set msg [list "Viewfactors should be updated (MESH is more recent)!\n\n" \
                      $Info(noMoreNotifying) ] 
      
        set Info(messageIcon) warning
        set result [Message::dispCancelYesNoMessage $msg]

        if { $result == "cancel" } {
          return "cancel"
        } elseif { $result == "yes" } {
          set Model(Viewfactors,needsUpdate,mesh) -1
          return "ok"
        }
      }
    }
  }

  #---Proc Viewfactors checking

  #-Check that output file is writable
  set reserved_files ""

  set path [MenuExec::getMeshRelatedFile $Model(currentMeshName) VIEW_FACTORS]

  if { ![Util::isWritableFile $path] } {
    append reserved_files "View factors file "
  }

  if { $reserved_files != "" } {
    set msg [list "Cannot open the View factors output file for writing!\n" \
                  "Enter a new file name in the Problem/Data files panel or wait until\n" \
                  "the process which is using the file, has finished!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return "cancel"
  }

  return "ok"
}


# Proc check if radiation boundary condition is applied
# This is needed when allowing linearization button for
# HEAT-system solver
#
proc MenuExec::updateRadiationStatus  {} {
  global Info BoundaryCondition ObjectTable

  if { ![info exists BoundaryCondition(ids)] } {
    return
  }

  #-Extra fill characters around separators (normally spaces )
  set disp_sep $Info(dispListSeparatorIn)
  set extra1 $Info(dispListSeparatorExtra1)
  set extra2 $Info(dispListSeparatorExtra2)

  set appliedParamIds ""

  foreach id $ObjectTable(ids) {

    if { ![info exist ObjectTable($id,bc)] } {
      continue
    }

    set pid $ObjectTable($id,bc)

    if { $pid == $Info(NO_INDEX) } {
      continue
    }

    if { -1 == [lsearch $appliedParamIds $pid ] } {
      lappend appliedParamIds $pid
    }
  }

  foreach pid $BoundaryCondition(ids) {

    if { -1 == [lsearch $appliedParamIds $pid] } {
      continue
    }

    set param $BoundaryCondition($pid,data)

    # Form values list. Pick first var-type  
    # indicator (like (Sc)) from value and collect them into a list
    # which is used for the current parameter.
    set tmp [split $param $disp_sep]

    set values_list ""

    foreach val $tmp {
      set size [DataField::extractFieldSize cleaned_result $val]
      lappend values_list $cleaned_result
    }
    set variable_id "RADIATION"

    #--Find radiation field value
    set index [lsearch -regexp $values_list ^$extra1$variable_id]

    #-Field not found
    if {$index == -1} {
      continue
    }

    set data [lindex $values_list $index]

    #-Variable-id (like TEMP=) is trimmed away.
    regsub $extra1$variable_id $data {} rest

    #-Trim fill-character(s) extra1 before value-type indicator
    set rest [string trimleft $rest $extra1]

    #-Is value a procedure name?
    set is_proc 0
    if {$Info(procFlag) == [string range $rest 0 1]} {
      set is_proc 1
    }

    set rest [string trimleft  $rest $Info(valueFlag)$Info(procFlag)$extra2]
    set rest [string trimright $rest $extra2]
    set rest [string trimright $rest $Info(groupDataSeparator)] 
    set value $rest
    
    if {$value != "none"} {
      set Info(hasRadiation) 1
      return
    }
  }

  set Info(hasRadiation) 0
}


# Check possible open and modified panel before
# loading a new model
proc MenuExec::checkOpenPanels {} {
  global Common Info

  MenuExec::closeAllWindows

  # Check if some arries are still modified
  #
  set mod_arries ""

  foreach globArr $Common(allPanelArries) {
    upvar #0 $globArr theArr

    if { [info exist theArr(dataChanged)] && $theArr(dataChanged) ||
         [info exist theArr(dataModified)] && $theArr(dataModified)
       } {
      lappend mod_arries $globArr
    }
  }

  # Some still modified!
  if { $mod_arries != "" } {
    set msg [list "NOTE: Data for the following panels is modified\n:" \
                  "$mod_arries\n" \
                  $Info(anywayOk)]

    if { ![Message::verifyAction $Info(powerUser) $msg] } {
      return 0
    }
  }
      
  # Ok
  return 1
}


#####################################################################
#
#                      FILE MENU  procs                            
#
#####################################################################

#------------#
# OPEN procs #
#------------#

# Open Cad input file.
proc MenuExec::openCadFile {} {
  global Info Model ModelFlags ModelProperty

  if { "ok" != [MenuExec::askModelOk "exit"] } {
    return
  }

  if { ![MenuExec::checkOpenPanels] } {
    return
  }

  # Supported Cad types and their default extensions
  # -------------------------------------------------
  set cad_types {"Elmer" "Iges" "Ideas"}
  set rb_labels {"Elmer 2D" "Iges 2D" "Ideas 2D"}

  set cad_extensions {
    { {Elmer} {.egf .Egf .EGF} }
    { {Iges} {.igs .Igs .IGS} }
    { {Ideas} {.unv .Unv .UNV} }
  }

  #--Build extensions list for the file browser

  # Helper variable
  set all_extension  {{All} {*}}

  #---Collect all cad extension into one list to be used for the 
  #   Default-option
  #
  #-First add 'All'-option
  set default_extensions [list $all_extension]

  #-Then add each separate mesh-type as a separate option
  foreach me $cad_extensions {
    lappend default_extensions $me
  }

  #---Now build the final extensions-list to be used for the file browser
  #
  set extensions ""

  #-We have to add 'All' option for each separate cad-type
  foreach me $cad_extensions {
    lappend extensions [list $me $all_extension]
  }

  #--Display FileOpen dialog and get the path
  #
  global cad_path_ cad_type_index_
  set cad_path_ ""
  set cad_type_index_ ""

  incr Info(fileSelectPanelId)
  set wname ".w$Info(fileSelectPanelId)"
  set hdr "Select Cad file type"

  set wtitle "Open Cad File"
  set path $Model(cadPath)
  set dir $Model(CAD_OPEN_DIR)
  set idx $Model(cadFileType)

  Widget::dispFileSelectPanel $wname $wtitle $hdr cad_path_ cad_type_index_ $rb_labels $extensions 0 $dir $path $idx
  tkwait variable cad_path_

  set cad_path $cad_path_
  set cad_type [lindex $cad_types $cad_type_index_]

  #--Nothing selected or cancelled
  #
  if {$cad_path == "" || $cad_path == "!cancel!" } {
    return
  }

  MenuExec::presetModelData

  set Model(cadPath) $cad_path
  set Model(CAD_OPEN_DIR) [file dirname $cad_path]
  set Model(cadFileType) $cad_type_index_

  Util::cpp_exec openCadFile "\"$cad_path\" \"$cad_type\""

  #-Example of a direct call when using a slave interpreter
  #set Info(arguments) $cad_path
  #cpp_openCadFile

  if { $Model(GEOMETRY_DIMENSION) != "" } {
    ModelProperty::openPanel
  }
} 



# Open mesh file.
proc MenuExec::openMeshFile {} {
  global Info MeshDefine Model ModelFlags ModelProperty

  if { "ok" != [MenuExec::askModelOk "exit"] } {
    return
  }

  if { ![MenuExec::checkOpenPanels] } {
    return
  }

  set Model(meshOpenMode) "new_model"

  # If opening mesh to an existing model which is based on an
  # external mesh (ie. mesh geometry)
  #
  if { !$ModelFlags(GEOMETRY_TYPE_CAD) &&
       $ModelProperty(MODEL_NAME) != "" 
     } {
     
    # Note: Model(meshOpenMode) is updated
    set win .meshOpenWin
    MenuExec::openMeshDialog $win Model meshOpenMode

    tkwait window $win

    if { $Model(meshOpenMode) == "cancel" } {
      return
    }
  }
  
  # Supported mesh types and their default extensions
  # -------------------------------------------------
  set mesh_types {
    "Abaqus2D" "Abaqus3D"
    "Elmer2D" "Elmer3D"
    "Fidap2D" "Fidap3D"
    "Ideas2D" "Ideas3D"
  }
  set rb_labels {
    "Abaqus 2D" "Abaqus 3D"
    "Elmer 2D" "Elmer 3D"
    "Fidap neutral file 2D" "Fidap neutral file 3D"
    "Ideas universal file 2D" "Ideas universal file 3D"
  }

  set mesh_extensions {
    { {Abaqus2D} {.inp .Inp .INP} }
    { {Abaqus3D} {.inp .Inp .INP} }
    { {Elmer2D} {mesh.header Mesh.header Mesh.Header MESH.HEADER} }
    { {Elmer3D} {mesh.header Mesh.header Mesh.Header MESH.HEADER} }
    { {Fidap2D} {.fdneut .Fdneut .FDNEUT} }
    { {Fidap3D} {.fdneut .Fdneut .FDNEUT} }
    { {Ideas2D} {.unv .Unv .UNV} }
    { {Ideas3D} {.unv .Unv .UNV} }
  }

  if { $Info(ELMER_FRONT_THETIS_SUPPORT) } {
    lappend mesh_types "Thetis2D Thetis3D"
    lappend rb_labels "Thetis 2D" "Thetis 3D"
    lappend mesh_extensions {{Thetis2D} {.tmf .Tmf .TMF}}
    lappend mesh_extensions {{Thetis3D} {.tmf .Tmf .TMF}}
  }


  #--Build extensions list for the file browser

  # Helper variable
  set all_extension  {{All} {*}}

  #---Collect all mesh extension into one list to be used for the 
  #   Default-option
  #
  #-First add 'All'-option
  set default_extensions [list $all_extension]

  #-Then add each separate mesh-type as a separate option
  foreach me $mesh_extensions {
    lappend default_extensions $me
  }

  #---Now build the final extensions-list to be used for the file browser
  #
  set extensions ""

  #-We have to add 'All' option for each separate mesh-type
  foreach me $mesh_extensions {
    lappend extensions [list $me $all_extension]
  }


  #--Reopen current mesh (eg. a remeshed external mesh!)
  if { $Model(meshOpenMode) == "update_mesh" } {
    set fn [lindex [file split $Model(meshPath)] end]
    set wtitle "Reload Mesh File (current is: $fn)"
    set dir [file dirname $Model(meshPath)]
    set path $Model(meshPath)
    set idx $Model(meshFileType)

  #--Read new mesh
  } else {
    set wtitle "Open Mesh File"
    set dir $Model(MESH_OPEN_DIR)
    set path $Model(meshPath)
    set idx $Model(meshFileType)
  }

  #--Display FileOpen dialog and get the path
  #
  global mesh_path_ mesh_type_index_
  set mesh_path_ ""
  set mesh_type_index_ ""

  incr Info(fileSelectPanelId)
  set wname ".w$Info(fileSelectPanelId)"
  set hdr "Select mesh file type"

  #--Display selection widget
  #
  Widget::dispFileSelectPanel $wname $wtitle $hdr mesh_path_ mesh_type_index_ $rb_labels $extensions 0 $dir $path $idx
  tkwait variable mesh_path_

  set mesh_path $mesh_path_
  set mesh_type [lindex $mesh_types $mesh_type_index_]

  #--Nothing selected or cancelled
  #
  if {$mesh_path == "" || $mesh_path == "!cancel!" } {
    return
  }

  set Model(nofMeshZeroVelocityElements) 0

  #--For a new model, reset all data
  if { $Model(meshOpenMode) == "new_model" } {
    MenuExec::presetModelData
  }

  set directory [file dirname $mesh_path]
  set filename [file tail $mesh_path]
  set extension [file extension $filename]

  set Model(MESH_OPEN_DIR) $directory
  set Model(meshPath) $mesh_path
  set Model(meshFileType) $mesh_type_index_

  Panel::initFields ModelProperty

  #--New model 
  if { $Model(meshOpenMode) == "new_model" } {

    Util::cpp_exec openMeshFile "1 \"$Model(meshPath)\" \"$mesh_type\""

    if { $Model(GEOMETRY_DIMENSION) != "" } {
      ModelProperty::openPanel
      tkwait window $ModelProperty(winName)
    }

  #--New mesh added or updated
  } else {
    Util::cpp_exec openMeshFile "0 \"$Model(meshPath)\" \"$mesh_type\""
  }

  #--For these we need a new mesh name
  #
  if { ($ModelProperty(MODEL_NAME) != "" )     &&
       ($Model(meshOpenMode) == "new_model" ||
        $Model(meshOpenMode) == "add_mesh")
     } {

    set add_new_mesh 1

    set rc [MeshDefine::openPanel $add_new_mesh]

    tkwait window $MeshDefine(winName)

    if { $rc != "cancel" } {

      set Model(Mesh,added) 1

      #--Set mesh name and directory
      set Model(currentMeshName) $MeshDefine(MESH_NAME)
      ModelProperty::setCurrentMeshDirectory
    }
  }
} 


# OLD VERSIONS: Cad/Mesh Open
# ===========================
# Open Cad input file.
proc MenuExec::openCadFile_old {} {
  global Info Model ModelFlags ModelProperty

  if { "ok" != [MenuExec::askModelOk "exit"] } {
    return
  }

  if { ![MenuExec::checkOpenPanels] } {
    return
  }

  set types {
    { {All} {*} }
    { {Elmer} {.egf .Egf .EGF} }
    { {Iges} {.igs .Igs .IGS} }
    { {Ideas} {.unv .Unv .UNV} }
  }
  set title "Open Cad File"

  #--Display dialog and get path
  set path [ MenuExec::openFile $title $types $Model(CAD_OPEN_DIR) ]


  if {$path == ""} {
    return
  }

  MenuExec::presetModelData

  set Model(cadPath) $path
  set Model(CAD_OPEN_DIR) [file dirname $path]

  Util::cpp_exec openCadFile $Model(cadPath)

  #-Example of a direct call when using a slave interpreter
  #set Info(arguments) $path
  #cpp_openCadFile

  if { $Model(GEOMETRY_DIMENSION) != "" } {
    ModelProperty::openPanel
  }
} 


# Open mesh file.
proc MenuExec::openMeshFile_old {} {
  global Info MeshDefine Model ModelFlags ModelProperty

  if { "ok" != [MenuExec::askModelOk "exit"] } {
    return
  }

  if { ![MenuExec::checkOpenPanels] } {
    return
  }

  set Model(meshOpenMode) "new_model"

  # If opening mesh to an existing model which is based on an
  # external mesh (ie. mesh geometry)
  #
  if { !$ModelFlags(GEOMETRY_TYPE_CAD) &&
       $ModelProperty(MODEL_NAME) != "" 
     } {
     
    # Note: Model(meshOpenMode) is updated
    set win .meshOpenWin
    MenuExec::openMeshDialog $win Model meshOpenMode

    tkwait window $win

    if { $Model(meshOpenMode) == "cancel" } {
      return
    }
  }

  set types {
    { {All} {*} }
    { {Abaqus} {.inp .Inp .INP} }
    { {Elmer}  {mesh.header Mesh.header Mesh.Header} }
    { {Fidap}  {.fdneut .Fdneut .FDNEUT} }
    { {Ideas}  {.unv .Unv .UNV} }
  }

  if { $Info(ELMER_FRONT_THETIS_SUPPORT) } {
    lappend types { {Thetis} {.tmf .Tmf .TMF} }
  }

  set path ""
  
  #--Reopen current mesh (eg. a remeshed external mesh!)
  if { $Model(meshOpenMode) == "update_mesh" } {

    set fn [lindex [file split $Model(meshPath)] end]
    set title "Reopen Mesh File (current is: $fn)"
    set dir [file dirname $Model(meshPath)]

  #--Read new mesh
  } else {
    set title "Open Mesh File"
    set dir $Model(MESH_OPEN_DIR)
  }

  #--Display FileOpen dialog and get the path
  set path [ MenuExec::openFile $title $types $dir ]

  if {$path == ""} {
    return
  }

  set Model(nofMeshZeroVelocityElements) 0

  #--For a new model, reset all data
  if { $Model(meshOpenMode) == "new_model" } {
    MenuExec::presetModelData
  }

  set directory [file dirname $path]
  set filename [file tail $path]
  set extension [file extension $filename]

  set Model(MESH_OPEN_DIR) $directory
  set Model(meshPath) $path

  Panel::initFields ModelProperty

  # Mesh type is selected by the extension!
  #
  set mesh_type ""

  #--New model 
  if { $Model(meshOpenMode) == "new_model" } {

    Util::cpp_exec openMeshFile "1 \"$Model(meshPath)\" \"$mesh_type\" "

    if { $Model(GEOMETRY_DIMENSION) != "" } {
      ModelProperty::openPanel
      tkwait window $ModelProperty(winName)
    }

  #--New mesh added or updated
  } else {
    Util::cpp_exec openMeshFile "0 \"$Model(meshPath)\" \"$mesh_type\" "
  }

  #--For these we need a new mesh name
  #
  if { ($ModelProperty(MODEL_NAME) != "" )     &&
       ($Model(meshOpenMode) == "new_model" ||
        $Model(meshOpenMode) == "add_mesh")
     } {

    set add_new_mesh 1

    set rc [MeshDefine::openPanel $add_new_mesh]

    tkwait window $MeshDefine(winName)

    if { $rc != "cancel" } {

      set Model(Mesh,added) 1

      #--Set mesh name and directory
      set Model(currentMeshName) $MeshDefine(MESH_NAME)
      ModelProperty::setCurrentMeshDirectory
    }
  }
} 


proc MenuExec::setRotatePriority {axis} {
  global Info Model

  if { $Model(GEOMETRY_DIMENSION) == "2D" &&
       ($axis == "X" || $axis == "Y")
     } {
   return
  }

  # Toggle current value if needed
  if { $Info(currentRotateAxis) == $axis } {
    set Info(currentRotateAxis) ""

  } else {
    set Info(currentRotateAxis) $axis
  }

  MenuExec::rendererSetRotatePriorities
}


# Note: updates globArray(varname) variable
#
proc MenuExec::openMeshDialog { winname arrayname varname} {
  global Info Model ModelFlags

  set var $arrayname
  append var "($varname)"

  set w $winname

  toplevel $w
  focus $w 

  wm geometry $w +500+100
  wm title $w "Open mesh file"

  set fpx $Info(framePadX2)
  set fpy $Info(framePadY2)

  #--Widgets
  frame $w.fl
  frame $w.fb

  #-Boundaries edited, not possible to update or add new meshes!
  if { $ModelFlags(GEOMETRY_EDITED_BOUNDARIES) } {
    set add_state disabled
    set upd_state disabled

  #-Ok, new meshes accepted
  } else {
    set add_state normal
    set upd_state normal
  }

  label $w.fl.l1 -text "\n\nSelect mesh open mode\n\n"
  pack $w.fl.l1 -expand 1 -side top  

  button $w.fb.b1 -text "New model" \
                  -command "set $var new_model; destroy $w; return"

  button $w.fb.b2 -text "Update mesh" \
                  -command "set $var update_mesh; destroy $w; return" \
                  -state $upd_state

  button $w.fb.b3 -text "Add mesh" \
                  -command "set $var add_mesh; destroy $w; return" \
                  -state $add_state

  button $w.fb.b4 -text "Cancel" \
                  -command "set $var cancel; destroy $w; return"

  pack $w.fb.b1 $w.fb.b2 $w.fb.b3 $w.fb.b4 -side left -padx $fpx  
  
  focus $w.fb.b1

  #--Packing
  pack $w.fb $w.fl -expand 1 -side bottom -fill both -padx $fpx -pady $fpy 
}


# Read model file.
#
proc MenuExec::openModelFile {} {
  global Info Model ModelProperty UserSetting

  if { "ok" != [MenuExec::askModelOk "exit"] } {
    return
  }

  if { ![MenuExec::checkOpenPanels] } {
    return
  }

  set types {
    { {Emf}  {.emf .Emf .EMF} }
    { {All}  {*} }
  }

  set title "Open Model File"

  #--Display dialog and get path
  set path [ MenuExec::openFile $title $types $Model(EMF_OPEN_DIR) ]

  if {$path == ""} {
    return
  }

  #--Load model file
  MenuExec::presetModelData

  if { ![MenuExec::loadNewModel $path] } {
    return
  }

  MenuExec::updateModelData

  # Delete possible old log files
  #
  set msg ""
  set log_files [Util::getFileList $ModelProperty(LOG_DIRECTORY,absolute) Elmer*.log.* ]

  catch {
    foreach lf $log_files {
      set lf [file join $ModelProperty(LOG_DIRECTORY,absolute) $lf]
      file delete $lf
    }
  } msg

  if { "" != $log_files && $msg == "" } {
    Message::showMessage "NOTE: Old logfiles deleted!"
  }

	#--Load possible updated SOLVER.KEYWORD files
	UserDefined::loadSolverKeywordsFiles
} 


# Load model file (emf-file)
#
# Return: 1 <--> ok, 0 <--> error
#
proc MenuExec::loadNewModel { model_file_path } {
  global Info Model ModelProperty
   
  #-Get model file directory and path
  set Model(inModelFile) $model_file_path
  set ModelProperty(MODEL_DIRECTORY) [file dirname $model_file_path]
  set Model(EMF_OPEN_DIR) [file dirname $model_file_path]
  set Model(EMF_SAVE_DIR) [file dirname $model_file_path]

  #-Apply new model directory value
  ModelProperty::applyModelDirectoryValue

  #-NOTE: EMF_PATH_IN/OUT, EMF_FILE are set via cpp, by calling
  # Gui after the model file name is checked in cpp
  Util::cpp_exec openModelFile $model_file_path

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return 0
  }

  #-Apply new problem directory values
  set force_init 0
  set show_msg 0
  MenuExec::applySettingsToCaseDirectories $force_init $show_msg

  #-Set at least initial value
  set force_initial 0
  Panel::initField ModelProperty MESH_NAME $force_initial

  #-Store checked values into cpp-model
  Util::cpp_exec modelPropertyPanelOk

  return 1
}


# Read mesh into model
proc MenuExec::loadMesh { {number ""} } {
  global Model ProcessTable

  if { $number == "" ||
       ![info exists ProcessTable($number,terminated)] ||
       $ProcessTable($number,terminated)
     } {

    Util::setFlagValue DRAW_TARGET DRAW_TARGET_BODIES 1

    set Model(nofMeshZeroVelocityElements) 0

    Util::cpp_exec loadMesh

  } else {
    after 1000 "MenuExec::loadMesh $number"
  }

}


# Delete mesh from model
proc MenuExec::unloadMesh { {msg ""} } {
  Util::cpp_exec unloadMesh $msg
}


proc MenuExec::monitorMesh {number} {
  global ProcessTable MeshDefine Model

  if { [info exists ProcessTable($number,terminated)] &&
       !$ProcessTable($number,terminated)
     } {
    after 500 "MenuExec::monitorMesh $number"

  } else {
    if { [info exist MeshDefine(solveButton)] &&
         [winfo exists $MeshDefine(solveButton)]
       } {
      $MeshDefine(solveButton) configure -state normal
    }
  }
}


proc MenuExec::resetModelFlags {} {
  global Model
  set Model(Front,needsUpdate) 0
  set Model(GebhardtFactors,needsUpdate,data) 0
  set Model(GebhardtFactors,needsUpdate,mesh) 0
  set Model(Solver,needsUpdate) 0
  set Model(Viewfactors,needsUpdate,data) 0
  set Model(Viewfactors,needsUpdate,mesh) 0
}


proc MenuExec::setModelFile {} {
  global Info Model

  set path ""
  set title "Set default name for the model file"
  set label "Model file name: "

  if { "ok" != [Widget::dispDataEntry path $Model(EMF_FILE) $title $label] } {
    return
  }

  set Model(EMF_PATH_IN) $path
  set Model(EMF_PATH) $path
  set Model(EMF_FILE) [file tail $path]
  set Model(EMF_SAVE_DIR) [file dirname $path]
}


# Copy parameters from an emf-file
#
proc MenuExec::copyParameters {} {
  global Info Model UserSetting

  set types {
    { {Emf}  {.emf .Emf .EMF} }
    { {All}  {*} }
  }

  set title "Select Model File To Copy Parameters"

  #--Display dialog and get path
  set path [ MenuExec::openFile $title $types $Model(EMF_OPEN_DIR) ]

  if {$path == ""} {
    return
  }

  MenuBuild::configureMenuButtonOption File copy_parameters state disabled

  Util::cpp_exec copyParameters $path
} 


#------------#
# SAVE procs #
#------------#

# ==========
# MODEL FILE
# ==========

# Save model file with existing names
# Check names and directories
#
proc MenuExec::saveModelFile {} {
  global Info Model ModelFlags ModelProperty

  if {"ok" != [MenuExec::checkDBName 1] } {
    return "cancel"
  }

  #--Check what to do with removed fields if the original model
  #  has definitions
  if { $Model(removedFieldsList) != "" && $Model(hasUserDefinitions) } {

    set msg [list "Note: Some fields have been removed from the model\n\n"         \
                  "(possibly owing to a missing definition file)\n\n"              \
                  "and they will removed permanently, if model file is saved!\n\n\n"  \
                  $Info(anywayOk) ]

    set Info(messageIcon) warning
    set result [ Message::dispOkCancelMessage $msg ]

    if { $result == "cancel" } {
      return "cancel"
    }
  }

  #--Check model directory
  if { $ModelProperty(MODEL_DIRECTORY) == "" } {
    set title "Set model directory"
    set label "Enter directory name"
    set dir [MenuExec::checkDirectory "" $title $label]

    if { [file isdirectory $dir] } {
      set ModelProperty(MODEL_DIRECTORY) $dir
      set ModelProperty(MODEL_DIRECTORY,absolute) [Util::makePathAbsolute $dir]
      
      # Store new values into Model
      Util::cpp_exec modelPropertyPanelOk
      
    } else {
      return "cancel"
    }
  }

  set needed 1
  if { -1 == [ModelProperty::checkModelDirectory $needed] } {
    return "cancel"
  }

  # If no mesh name defined for an external mesh
  if { $ModelFlags(GEOMETRY_TYPE_MESH) && 
       $Model(meshNames) == ""
     } {
    set msg "No mesh name defined!"
    set Info(messageIcon) error
    Message::dispOkMessage  $msg
    return "cancel"
  }

  set rc "ok"

  #--Default is given
  set msg ""
  if {[catch { set rc [MenuExec::saveCheckedModelFile $Model(EMF_PATH)] } msg]} {
    set msg [list "Could not save data because:" $msg]

    set Info(messageIcon) error
    Message::dispOkMessage  $msg
    return "cancel"
  }

  return $rc
}


# Ask user for the name of the model file
# NOTE: Not in use!
#
proc MenuExec::saveModelFileAs {} {
  global Info Model ModelProperty

  #--Check model name
  if {"ok" != [MenuExec::checkDBName 1] } {
    return "cancel"
  }

  #--Check model directory
  if { $ModelProperty(MODEL_DIRECTORY) == "" } {
    set title "Set model directory"
    set label "Enter directory name"
    set dir [MenuExec::checkDirectory "" $title $label]

    if { [file isdirectory $dir] } {
      set ModelProperty(MODEL_DIRECTORY) $dir
    } else {
      return "cancel"
    }
  }

  set needed 1
  if { -1 == [ModelProperty::checkModelDirectory $needed] } {
    return "cancel"
  }

  #--Model file extension
  set types {
    { {Ecif Files}  {.emf} }
    { {All Files}  {*} }
  }
  set title "Save Model File"

  #--Display dialog and get path
  set path [MenuExec::saveFile $title $types $Model(EMF_PATH)]

  if {$path == ""} {
    return "cancel"
  }

  set Model(EMF_PATH) $path
  set Model(EMF_FILE) [file tail $path]
  set Model(EMF_SAVE_DIR) [file dirname $path]

  set rc "ok"

  set msg ""
  if {[catch { set rc [MenuExec::saveCheckedModelFile $Model(EMF_PATH)] } msg]} {
    set msg [list "Could not save data because:" $msg]

    set Info(messageIcon) error
    Message::dispOkMessage  $msg
    return "cancel"
  }

  return $rc
}


# Saves model file
# Names and directories are checked!
#
proc MenuExec::saveCheckedModelFile {path {model_file_only 0} {show_msg 1}} {
  global Info Model UserSetting ModelProperty
  global ModelFlags

  # Set possible user definitions status
  # ====================================
	#
  if { $Model(allDefinitions) != "" } {

    set Model(hasUserDefinitions) 1

    # NOTE: Currently no automatic save!
    #
    #set def_path [file root $path]
    #append def_path ".edf"

    #-We always overwrite any existing def file!    
    #set append 0

    # Def-file not saved any more with the model! MVe 16.02.01
    #if { [MenuExec::saveDefinitionFile $def_path $append] } {
    #  set Model(outDefinitionFile) $def_path
    #}

  }

	# Load possible updated SOLVER.KEYWORD files
	# ==========================================
	#
	UserDefined::loadSolverKeywordsFiles


  # Update model file info into cpp-side
  # ====================================
	#
  set save_time [Util::setCurrentTimestamp]

  #-Model file history
  set Model(modified) "$save_time User=$Info(USER) Host=$Info(HOST)"
  Util::cpp_exec readModelFileModified $Model(modified)

  if { $Model(created) == "" } {
    set Model(created) $Model(modified)
    Util::cpp_exec readModelFileCreated $Model(created)
  }

  Util::cpp_exec readModelHasUserDefinitions $Model(hasUserDefinitions)

  #-Read model status indicator from C++
  Util::cpp_exec checkModelStatus
  
  #-ModelFile save time
  Util::cpp_exec readModelFileTime $save_time

  #-Set Front update timestamp!
  if { $Model(Front,needsUpdate) } {
    set Model(Front,timestamp) $save_time
    Util::cpp_exec readTimestamp [list "Front" $save_time]
  }

  #-Set Database ts already here, if it will be anyway
  # updated!
  # NOTE: Database itself must be updated after model file
  # is saved, because it will need the model file contents!

  if {!$model_file_only} {
    set Model(Database,timestamp) $save_time
    Util::cpp_exec readTimestamp [list "Database" $save_time]
  }

  # Save model file
  # ===============
  # NOTE: This must be done last, because otherwise latest values for
  # status, model-time etc are not saved into the file!
  #

  #-Make a backup before saving!
  set path_tmp ""

  if { [file exist $path] } {
    set path_tmp $path.tmp
    file copy -force $path $path_tmp
  }

  Util::cpp_exec saveModelFile $path

  if {$show_msg} {
    Message::showMessage "Model file saved as $path"
  }

  #-Remove the backup after saving!
  catch {file delete -force $path_tmp}

  set Model(Front,needsUpdate) 0
  set Model(removedFieldsList) ""

  # If updating only model file
  if {$model_file_only} {
    return "ok"
  }

  # Update problem directory state variables
  # ========================================
	#
  ModelProperty::updateCaseDirectoryVariables

  # Mesh saving is conditional
  # ==========================
  # Check if the Elmer format mesh should be stored:
  # if we have an external mesh or
  # user has renamed current mesh
  #
  if { !$Model(Mesh,generating) } {

    if { ( $ModelFlags(GEOMETRY_TYPE_MESH) && !$Model(Mesh,exists) ) ||
         ( $Model(Mesh,exists) && ( $Model(Mesh,added) || $Model(Mesh,edited)) )
       } {
  
      # If no Elmer mesh
      #
      if {!$Model(Mesh,exists)} {
           
        if {$UserSetting(DEFAULT_AUTO_SAVE_EXTERNAL_MESH)} {
          # Ask user if we should save the mesh in Elmer format
          set msg [list "NOTE: Elmer mesh does not exist\n\n" \
                        "The external mesh will be saved in Elmer format\n" \
                        "using the name: $Model(currentMeshName) \n\n" \
                         $Info(continueOk) ]

          if { [Message::verifyAction $Info(powerUser) $msg] } {
            MenuExec::saveElmerMeshFile
          }
        }

      # Elmer mesh exists, but it has been
      # edited (split/combined) or renamed
      } elseif { $Model(Mesh,added) || $Model(Mesh,edited) } {

        MenuExec::saveElmerMeshFile
      }
    }

  } ;# End mesh saving related


  # Solver sif-file update is conditional
  # =====================================

  #-Solver input file saving rc!
  set rc "ok"

  set solver_input_file [MenuExec::getSolverInputFile]

  #-If auto-save flag is off, we save only if Solver
  # input file does not exist
  #
  if { !$UserSetting(AUTO_SAVE_SOLVER_INPUT) } {

    MenuExec::checkSolverInputFile

    if { !$Model(Solver,inputFileExists) } {
      set rc [MenuExec::saveSolverInputFile $solver_input_file]
    }

  #-Otherwise we always save
  #
  } else {
    set rc [MenuExec::saveSolverInputFile $solver_input_file]
  }

  return $rc
}



# =================
# SOLVER INPUT FILE
# =================

proc MenuExec::saveSolverInputFile {path} {
  global Model Info

  set result [MenuExec::askSolverOk $path]

  # User cancelled sif-file saving
  #
  if { $result == "cancel" } {
    return "cancel"
  }

  set dir [file dirname $path]

  set dir_msg "solver input file"

  if { "cancel" == [Util::checkAndCreateDirectory2B $dir $dir_msg] } {
    return "cancel"
  }

  # Make a backup before saving!
  set path_tmp ""
  if { [file exist $path] } {
    set path_tmp $path.tmp
    file copy -force $path $path_tmp
  }

  Util::cpp_exec saveSolverInputFile $path

  # Remove the backup after saving!
  catch {file delete -force $path_tmp}

  set Model(Solver,needsUpdate) 0
  set Model(Solver,inputFileExists) 1

  return "ok"

}




# ===========================================
# ELMER MESH FILE (in model or as "external")
# ===========================================

# Save mesh result file.
proc MenuExec::saveElmerMeshFile { {as_external 0} } {
  global Info Model ModelProperty ModelFlags

  # As "External" mesh, ask the directory
  # =================
  if {$as_external} {
    set title "Set mesh directory"
    set label "Enter a directory name"
    set dir [MenuExec::checkDirectory "" $title $label]

    if { $dir == "!cancel!" } {
      return "cancel"
    }  

    # Check that mesh directory exists
    set ask 1
    if {"ok" != [Util::checkAndCreateDirectory2B $dir "Mesh files" $ask]} {
      return "cancel"
    }

  # As Model mesh
  # =============
  } else {
    ModelProperty::setCurrentMeshDirectory

    #---Check that mesh directory exists
    set ask 0
    if {"ok" != [Util::checkAndCreateDirectory2B \
                   $ModelProperty(CURRENT_MESH_DIRECTORY,absolute) \
                   "Elmer mesh files" $ask]} {
      return "cancel"
    }

    set dir $ModelProperty(CURRENT_MESH_DIRECTORY,absolute)
  }

  set Model(Mesh,exists) 1
  set Model(Mesh,edited) 0
  set Model(Mesh,added) 0

  Util::cpp_exec saveElmerMeshFile $dir
}


# ==================================
# ELEMRPOST MESH (mesh in ep-format)
# ==================================

# Save mesh result file.
proc MenuExec::saveElmerPostMeshFile {} {
  global Info Model

  set types {
    { {Mesh Files} {.ep .EP .Ep .p .P} }
    { {All Files} {*} }
  }
  set title "Save mesh in ElmerPost format"

  set default_path $Model(outElmerPostMeshFile)

  #--Display dialog and get path
  set path [MenuExec::saveFile $title $types $default_path]

  if {$path == ""} {
    return
  }

  set Model(outElmerPostMeshFile) $path

  Util::cpp_exec saveElmerPostMeshFile $path
}



# ===================================
# THETIS MESH (mesh in thetis format)
# ===================================

# Save mesh result file.
proc MenuExec::saveThetisMeshFile {} {
  global Info Model
  set types {
    { {Thetis Mesh Files} {.tmf} }
    { {Text Files} {.txt } }
    { {All Files} {*} }
  }
  set title "Save mesh in Thetis format"

  set default_path $Model(outThetisMeshFile)

  #--Display dialog and get path
  set path [MenuExec::saveFile $title $types $default_path]

  if {$path == ""} {
    return
  }
  
  set Model(outThetisMeshFile) $path

  Util::cpp_exec saveThetisMeshFile $path
}



#-----------------------#
# OTHER File menu procs #
#-----------------------#

proc MenuExec::browseFile { {edit_mode 0} } {
  global Info

  set types {
    { {All} {*} }
  }

  if { $edit_mode } {
    set title "Edit File"
  } else {
    set title "Browse File"
  }

  if { $edit_mode } {
    set open_dir $Info(editorDir)
    set open_fil $Info(editorPath)
  } else {
    set open_dir $Info(browserDir)
    set open_fil $Info(browserPath)
  }

  #--Display dialog and get path
  set path [MenuExec::openFile $title $types $open_dir $open_fil]

  if { $path == "" } return

  if { $edit_mode } {
    set Info(editorDir) [file dirname $path]
    set Info(editorPath) $path
  } else {
    set Info(browserDir) [file dirname $path]
    set Info(browserPath) $path
  }

  if { $edit_mode } {
    FileBrowser::browse .winFileEedit $path 1 1
  } else {
    FileBrowser::browse .winFileBrowse $path 
  }
}


proc MenuExec::browseEgfFile { {edit_mode 0} } {
  global Info Model

  if { $Model(cadPath) != "" } {
    set path $Model(cadPath) 

  } else {
    set types {
      { {Egf} {.egf .Egf .EGF} }
      { {All} {*} }
    }

    if { $edit_mode } {
      set title "Edit Egf File"
    } else {
      set title "Browse Egf File"
    }

    #--Display dialog and get path
    set path [MenuExec::openFile $title $types $Model(CAD_OPEN_DIR)]
  }

  if { $path == "" } return

  if { $edit_mode } {
    FileBrowser::browse .winEgfFilEedit $path 1 1
  } else {
    FileBrowser::browse .winEgfFileBrowse $path 
  }
}


# Browse current definitions
#
proc MenuExec::browseCurrentDefinitions { {edit_mode 0 } } {

  if { $edit_mode } {
    FileBrowser::browse .winDefBrowse {Model allDefinitions} 0 1 1 
  } else {
    FileBrowser::browse .winDefEdit {Model allDefinitions} 0
  }
} 


# NOTE: No edit-mode for model files!
#
proc MenuExec::browseEmfFile {} {
  global Info Model

  # For event testing:  
  # global InitialCondition mainWindow
  # event generate $InitialCondition(winName)  <Configure> -above  $mainWindow  
  # event generate $InitialCondition(winName)  <Circulate> -place PlaceOnTop  
 
  if {$Model(EMF_PATH) != ""} {
    FileBrowser::browse .winModelFileBrowse $Model(EMF_PATH) 

  } else {
    set msg [list "No model file defined!\n" \
                  "Command not valid" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
  }
}


proc MenuExec::browseSifFile { {edit_mode 0} } {
  global Info Datafile
 
  set solver_input_file [MenuExec::getSolverInputFile]

  if { $solver_input_file != "" } {

    if { ![file isfile $solver_input_file] } {
      set msg [list "File not found:\n" \
                    "$solver_input_file\n" \
                    "Save first model (File/Save model file) !" ]
                    
      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return
    }

    if { $edit_mode } {
      FileBrowser::browse .winSolverFileEdit $solver_input_file 1 1
    } else {
      FileBrowser::browse .winSolverFileBrowse $solver_input_file 
    }

  } else {
    set msg [list "No solver input file defined!\n" \
                  "Define solver input file (Problem/Datafiles)\n" \
                  "and save model!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
  }
}


proc MenuExec::browseSolverReloadFile { {edit_mode 0} } {
  global Info Datafile

  set reload_file [MenuExec::getSolverReloadFile]
 
  if { $reload_file != "" } {

    if { ![file isfile $reload_file] } {
      set msg [list "File not found:\n" \
                    "$reload_file !"]
                    

      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return
    }

    if { $edit_mode } {
      FileBrowser::browse .winSolverReloadEdit $reload_file 1 1
    } else {
      FileBrowser::browse .winSolverReloadBrowse $reload_file 
    }

  } else {
    set msg [list "No solver reload file defined!\n" \
                  "Define solver input file (Solver/Solver control)\n" \
                  "and save model!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
  }
}


proc MenuExec::browseLogFile {} {
  global Info
 
  FileBrowser::browse .winLogFile $Info(ElmerRun,logfile)
}




#-Program EXIT
#
proc MenuExec::cifExit {} {
  global Info Model
  global ProcessIds ProcessTable

  set do_exit 1

  catch {
  ###### Start of CATCH #####

  set result [MenuExec::askModelOk "exit"]

  if { $result  != "ok" } {
    set do_exit 0
    return
  }

  #---Check unstopped processes and delete possible logfiles
  #
  set p_info ""
  foreach nbr $ProcessIds {

    # List unstopped process
    if { !$ProcessTable($nbr,terminated) } {
      append p_info "Process: $ProcessTable($nbr,name) "
      append p_info "Nbr: $nbr "
      append p_info "PID: $ProcessTable($nbr,pid) "
      append p_info "Model: $ProcessTable($nbr,model) "

      if { $ProcessTable($nbr,logfile) != "" } {
        append p_info "Logfile: $ProcessTable($nbr,logfile) "
      }
      append p_info "\n"

    # Delete possible logfile
    } else {
      set msg ""
      catch { [file delete $ProcessTable($nbr,logfile)] msg }
    }
  }

  # Confirm exit, if active processes
  if { $p_info != "" } {
    set msg [list "NOTE: The following processes are still active!\n\n"  \
                  "$p_info\n\n" \
                  "$Info(anywayOk)" ]

    if { ![Message::verifyAction $Info(advancedUser) $msg] } {
      set do_exit 0
      return
    }
  }

  }
  ###### End of CATCH #####

  # If exit was cancelled by the user
  #
  if { !$do_exit } {
    return
  }

  #---Exit Elmer Front!!!
  Util::cpp_exec cpp_exit
  exit
}


#-Open file cover procedure
#
proc MenuExec::openFile {title types {default_dir ""} {default_file ""} } {

 #---Tk Open file dialog
  set path [tk_getOpenFile            \
             -title $title            \
             -filetypes $types        \
             -initialdir $default_dir \
             -initialfile $default_file \
           ]

  return $path
}


#-Save file cover procedure
#
proc MenuExec::saveFile {title types default_path} {

  set file [file tail $default_path]
  set dir [file dirname $default_path]

  #---Tk Save file dialog 
  set path [tk_getSaveFile      \
             -title $title      \
             -filetypes $types  \
             -initialdir $dir   \
             -initialfile $file \
           ]

  return $path
}


#####################################################################
#
#                    BUTTON AREA procs
#
#####################################################################

#---Renderer display button procs

#-Move (translate) model
#
proc MenuExec::rendererMove { direction {axis ""} } {
  global Info

  if { $axis == "" } {
    set axis $Info(currentMoveAxis)
  }

  #-Map axis to integer
  if {$axis == "X"} {
    set axis 0
  } elseif {$axis == "Y"} {
    set axis 1
  } else {
    set axis 2
  }

  set data [list $axis $direction]
  Util::cpp_exec rendererTranslateModel $data
}

#-Set rotation priorities
#         
proc MenuExec::rendererSetRotatePriorities { } {
  global Info

  set data {0 0 0}

  switch $Info(currentRotateAxis) {

    "X" {
      set data {1 0 0}
    }

    "Y" {
      set data {0 1 0}
    }

    "Z" {
      set data {0 0 1}
    }
  }

  Util::cpp_exec rendererSetRotatePriorities $data

}


#-Rotate model
#
proc MenuExec::rendererRotate { direction {axis ""} } {
  global Info

  if { $axis == "" } {
    set axis $Info(currentRotateAxis)
  }

  #-Map axis to integer
  #-None checked
  if {$axis == ""} {
    return
  #-X checked
  } elseif {$axis == "X"} {
    set axis 0
  #-Y checked
  } elseif {$axis == "Y"} {
    set axis 1
  #-Z checked
  } else {
    set axis 2
  }
   
  set data [list $axis $direction]
  Util::cpp_exec rendererRotateModel $data
}


#-Scale model
#
proc MenuExec::rendererScale {direction} {
  global Info

  set data $direction
  Util::cpp_exec rendererScaleModel $data
}


#-Reset graphics screen
#
proc MenuExec::rendererReset {} {
  Util::cpp_exec rendererResetModel
}


#-Display (refresh) graphics screen
#
proc MenuExec::rendererDisplay {} {
  global ModelFlags

  # Set at least bodies to be drawn!
  if { $ModelFlags(DRAW_TARGET_BODIES)    == 0  &&
       $ModelFlags(DRAW_TARGET_SURFACES)  == 0  &&
       $ModelFlags(DRAW_TARGET_EDGES)     == 0
     } {
    Util::setFlagValue DRAW_TARGET DRAW_TARGET_BODIES 1
  } 

  Util::cpp_exec rendererDisplayModel
}


proc MenuExec::show_results {} {
  global Model

  # If we have only one mesh
  if { 1 == $Model(nofActiveMeshes) } {

    set postFile [list [MenuExec::getSolverResultFile 0 POST_FILE]]
    MenuExec::prepare_elmer_exec Results $postFile

  } else {
    PostFileSelect::openPanel
  }
}


proc MenuExec::info_ {} {
  ModelInfo::openPanel
}


#-Stop Front task (by Break button)
#
proc MenuExec::doBreak {} {
  global Info

  #Message::showMessage "Break not active, sorry!"
  Util::cpp_exec doBreak
}




#####################################################################
#
#                       EDIT MENU  procs
#
#####################################################################

proc MenuExec::removeCadGeometry {} {
  global ModelFlags Model ModelProperty Info

  if {$ModelFlags(GEOMETRY_TYPE_CAD) == 0} {
    set msg "No Cad geometry available, nothing to remove!"

    set Info(messageIcon) info
    Message::dispOkMessage $msg
    return 
  }

  if {$ModelFlags(GEOMETRY_TYPE_MESH) == 0} {
    set msg "No other (mesh) geometry available, cannot remove Cad geometry!"

    set Info(messageIcon) info
    Message::dispOkMessage $msg
    return 
  }

  set msg [list \
     "This command removes Cad geometry from the model!\n\n" \
     "If you do that NOTE the following:\n\n" \
     "--You CANNOT create a new mesh, because geometry is lost!\n\n" \
     "--You CAN edit boundaries which will be based on mesh elements only!\n\n" \
     "Change model name before saving the changes (if you want to keep the original model)!\n\n" \
     "\nPress OK to remove Cad geometry, otherwise press Cancel" \
     ]

  if { ![Message::verifyAction $Info(powerUser) $msg] } {
    return
  }

  Util::cpp_exec removeCadGeometry

  set Model(Front,needsUpdate) 1

  Message::showMessage "Cad geometry removed!"

  # Update main windo title
  Util::updateMainWindowTitle
}


# This proc removes all inactive parameer field from
# all mask dependent parameters
#
proc MenuExec::removeInactiveParameterData {} {
  global Info Equation

  # NOTE: If turned off
  # ===================
  #set msg "Sorry, not in use currently!"
  #Message::dispOkMessage $msg
  #return
  # ===========================

  set msg [list \
     "This command removes all inactive parameter data from the model!\n\n" \
     "NOTE: you have to save the model after this command!\n\n" \
     "\nPress OK to remove inactive parameter data, otherwise press Cancel" \
     ]

  if { ![Message::verifyAction $Info(powerUser) $msg] } {
    return
  }

  #-Set flag value for skipping inactive data
  set tmp $Info(keepInactiveFields)
  set Info(keepInactiveFields) 0 

  #-Clean Equation parameters
  set modified1 0
  Panel::checkParameters Equation modified1

  #-Clean all mask related arries
  set modified2 0
  Panel::checkMaskedParameters modified2

  #-Clean solver data
  set modified2 0
  Solver::checkSolverData $Equation(activeIndices) modified2

  #-Restore flag
  set Info(keepInactiveFields) $tmp 
 
}


# Set Matc input file name to be stored in emf-file
#
proc MenuExec::setMatcInputFile {type} {
  global Info Model

  # Matc-input file for model (emf) file
  #
  if { $type == "emf" } {
    set title "Matc input for Model file"
    set mcf $Model(matcInputFile_emf)

  # Matc-input file for sif-file
  #
  } elseif { $type == "sif" } {
    set title "Matc input for Sif-file"
    set mcf $Model(matcInputFile_sif)
  }

  set label "Enter file name"

  global mcf_res ;# Reference var must be global in dispFileEntry!

  set mcf_res ""

  incr Info(fileEntryId)
  set wname ".w$Info(fileEntryId)"

  Widget::dispFileEntry $wname mcf_res $mcf $title $label 
  tkwait variable mcf_res

  if { $mcf_res != "!cancel!"  } {
    
    #-If ok value   
    if { [file isfile $mcf_res] } {

      set mcf $mcf_res

    #-Empty file name
    } elseif { $mcf_res == "" } {

      set msg [list "No Matc input file will be in use!\n\n" \
                     $Info(anywayOk) ]

      set Info(messageIcon) warning
      set result [ Message::dispOkCancelMessage $msg ]

      # Set empty name
      if { $result != "cancel" } {
        set mcf ""

      # Close entry
      } else {
        unset mcf_res
        return
      }

    #-Invalid file name
    } else {

      set msg [list "A non existing file name!\n\n$mcf_res\n\n" \
                     $Info(anywayOk) ]

      set Info(messageIcon) error
      set result [ Message::dispOkCancelMessage $msg ]

      # Use invalid name (user knows better!)
      if { $result != "cancel" } {
        set mcf $mcf_res

      # Close entry
      } else {
        unset mcf_res
        return
      }
    }

    # Set matc input file in model
    #
    if { $type == "emf" } {
      set Model(matcInputFile_emf) $mcf
      Util::cpp_exec setMatcInputFileEmf $Model(matcInputFile_emf)
    } elseif { $type == "sif" } {
      set Model(matcInputFile_sif) $mcf
      Util::cpp_exec setMatcInputFileSif $Model(matcInputFile_sif)
    }
  }

  set Model(Front,needsUpdate) 1
  Util::updateMainWindowTitle

  unset mcf_res
}


# Update model geometry (by Matc-definitions)
#
proc MenuExec::verifyUpdateCadGeometry  {} {
  global Info Model

  # Warning if mesh(es) exists
  if { $Model(Mesh,exists) } {
    set msg [list "Matc-dependent model geometry will be updated and" \
                  "\ncurrent meshes may not then match the geometry!\n\n" \
                  $Info(anywayOk) ]
  } else {
    set msg "Matc-dependent model geometry will be updated!"
  }

  set Info(messageIcon) warning
  set result [ Message::dispOkCancelMessage $msg ]
    
  # User cancelled
  if { $result == "cancel" } {
    return
  }

  # Update (Matc-dependent) geometry
  Util::cpp_exec updateCadGeometry
}

proc MenuExec::verifyDropMatcDefs {} {
  global Info Model

  # Warning
  if { $Info(dropModelMatcDefinitions) && $Model(hasMatcDefinitions) } {

    set msg [list "Model has Matc dependent geometry and these definitions will NOT be stored!\n\n" \
                   $Info(anywayOk) ]

    set Info(messageIcon) warning
    set result [ Message::dispOkCancelMessage $msg ]
    
    # User cancelled
    if { $result == "cancel" } {
      set Info(dropModelMatcDefinitions) 0
    }
  }

  # Set menu text
  MenuBuild::setCheckMenuLabel Edit "Drop Matc definitions from model" Info(dropModelMatcDefinitions)
}


proc MenuExec::setMeshInputUnit {} {
  global Model

  global unit
  set unit ""

  set title "Set mesh unit for external mesh input"
  set label "Enter input unit:"

  Widget::dispDataEntry unit $Model(meshInputUnit) $title $label
  tkwait variable unit

  set result $unit
  unset unit

  if { $result == "!cancel!" } {
    return
  }

  if { $result == "" || $result <= 0 || [catch {expr $result}] } {
    Message::dispOkMessage "ERROR: Invalid value: $result (must be a positive number)"
    return
  }

  Util::cpp_exec setMeshInputUnit $result

}

#####################################################################
#
#                   ELMER PROCESS EXEC procs
#
#####################################################################

proc MenuExec::send_WIN32 {target command} {
  return [dde eval $target $command]
}


##########################
### PREPARE ELMER EXEC ###
##########################

# This procedure constructs proper arguments for the calls
# to the external Elmer procedures
#
proc MenuExec::prepare_elmer_exec {process_name { process_args ""}  } {
  global Info Model ModelProperty Datafile Timestep
  global Processor UserSetting
  global ProcessTable

  #--General environment checking
  if { "ok" != [MenuExec::checkProcessEnvironment process_name] } {
    return -1
  }

  #--System specific call name specifier
  if { $Info(platform) == "windows" } {
    set sys WIN32
  } else {
    set sys UNIX
  }

  # Set program name and arguments
  # ==============================
  set exec_name $Info($process_name,proc,$sys)
  set arg ""
  set use_eval 1
  set priority_level $Info(bgPriorityLevel)
  set browser_delay 0
  set process_info ""
  set add_to_process_table 1
  set show_in_process_list 1
  set browse_parser "none"
  set browse_parser_needs_runnumber 1

  # GebhardtFactors
  # ---------------
  if { $process_name == $Info(GebhardtFactors,processName) } {

    set arg [MenuExec::getSolverInputFile]
    set browse_mode $UserSetting(BROWSE_MODE_GEBHARDT_FACTORS)
    set browser_delay $Info(browser,largeDelay)

  # Mesh
  # ----
  } elseif { $process_name == $Info(Mesh2D,processName) ||
             $process_name == $Info(Mesh3D,processName)
           } {
    
    # Create mif-file
    set mif_name [MenuExec::getMeshInputFileName]
    set fname [file join $ModelProperty(CURRENT_MESH_DIRECTORY,absolute) $mif_name]
    Util::cpp_exec saveMeshInputFile $fname

    set arg1 "\"$fname\" \"$ModelProperty(CURRENT_MESH_DIRECTORY,absolute)\""

    #-Control background mesh file in use
    #
    if { $Model(currentMeshBgMeshControlFile) != "" } {
      set arg2 "--bgcontrol=\"$Model(currentMeshBgMeshControlFile)\""

    #-Normal background mesh file in use
    #
    } elseif { $Model(currentMeshBgMeshFile) != "" } {
      set arg2 "--bgmesh=\"$Model(currentMeshBgMeshFile)\""

    #-No background mesh stuff
    } else {
      set arg2 ""
    }
    
    set arg [string trim "$arg1 $arg2"]

    set browse_mode $UserSetting(BROWSE_MODE_MESH)
    set browser_delay $Info(browser,largeDelay)
    set process_info [ProcessTable::formMeshInfo]
    set Model($process_name,ElmerMeshProcessing) 1

  # Results
  # -------
  # NOTE:
  # --Windows: we start first the server by a command line
  # call and then send the initializing commands via DDE in
  # proc MenuExec::reset_elmer_exec. This is because of the
  # limited size of the command line buffer.
  # --Unix: we send init level commands by the command line call
  # and use DDE only for refresh  commands
  #
  } elseif { $process_name == $Info(Results,processName) } {

    set postFile [lindex $process_args 0]
    set server $Info(Results,tclServer)

    # These variables are used to control ElmerPost
    # call mode
    if { ![info exists Info(Results,$server,mode)] } {
      set Info(Results,$server,mode) ""
      set Info(Results,$server,created) 0
    }

    #--Check server state and form proper ElmerPost command
    set rc [MenuExec::prepare_ElmerPost_exec $postFile                \
                                             $server                  \
                                             Info(Results,tclCommand)]

    #--If server exist, use send to deliver update commands
    if { $rc == 1 } {

      set Info(Results,$server,created) 0

      # Try to send update command
      if { ![MenuExec::sendTclCommands $server $Info(Results,tclCommand)  10 msg] } {
        Message::showMessage "$server not updated because: $msg" 

	    } elseif { $Info(displayExecMessages) } {
        Message::showMessage "$server updated: $Info(Results,tclCommand)" 
      }

      return 0

    #--If server doesn't exist, create it using a command line call!
    } elseif { $rc == 0 && $Info(Results,$server,mode) != "" } {

      set Info(Results,$server,created) 1

      set arg "-id"
      lappend arg $Info(Results,tclServer)
      lappend arg "-file"
      lappend arg $postFile
      lappend arg $Info(Results,tclCommand)

    #--Not found or exists, but not accessible!
    } else {
      set Info(Results,$server,created) 0
      return 0
    }

    #-Create a new post session
    set use_eval 0
    set priority_level "NORMAL_PRIORITY"
    set add_to_process_table 0
    set browse_mode "None"

  # Post
  # ----
  # NOTE: This is used when no model file is opened, then
  # "Results" is changed to "Post" command in checkProcessEnvironment!
  } elseif { $process_name == $Info(Post,processName) } {

    set arg ""
    set use_eval 0
    set priority_level "NORMAL_PRIORITY"
    set add_to_process_table 0
    set browse_mode "None"
 
  # Solver
  # ------
  # NOTE: We do not add model-directory to the solver-input-file
  # name, when we write solver startinfo data! MVe 29.09.99
  #
  } elseif { $process_name == $Info(Solver,processName) } {

    # For exe-version
    set startinfo_file [file join $ModelProperty(MODEL_DIRECTORY,absolute) "ELMERSOLVER_STARTINFO"]
    set ch [open $startinfo_file "w+"]
    puts $ch $Datafile(SOLVER_INPUT_FILE)
    puts $ch $Processor(NOF_PROCESSORS) 
    close $ch
    
    set reload_file [DataField::getFieldValue SolverControl 1 RELOAD_INPUT_FILE active]

    if { $active && $reload_file != "" } {
      set rereadinfo_file [file join $ModelProperty(MODEL_DIRECTORY,absolute) "ELMERSOLVER_REREADINFO"]
      set ch [open $rereadinfo_file "w+"]
      puts $ch $reload_file
      puts $ch $Processor(NOF_PROCESSORS) 
      close $ch
    }

    # For batch-version       
    set arg "$Datafile(SOLVER_INPUT_FILE)"
    append arg " $Processor(NOF_PROCESSORS)"
    set browse_mode $UserSetting(BROWSE_MODE_SOLVER)
    set browser_delay $Info(browser,largeDelay)
    set process_info [ProcessTable::formSolverInfo]
    set browse_parser "Process::solverBrowseParser"
    set browse_parser_needs_runnumber 1

  # Viewfactors
  # -----------
  } elseif { $process_name == $Info(Viewfactors,processName) } {

    set arg [MenuExec::getSolverInputFile]
    set browse_mode $UserSetting(BROWSE_MODE_VIEW_FACTORS)
    set browser_delay $Info(browser,largeDelay)
 
  # Error
  # -----
  } else { 
    set msg "Unknown program name: $process_name \n"
    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return -1
  }

  #--Complete name for the executable  
  set home  $Info(ELMER_HOME)
  set exec_dir  [file join $home "bin"] 
  set exec_cmd [file join $exec_dir $exec_name]

  #--Check that executable exists
  if { ![file executable $exec_cmd] } {
    set msg [list "Can't find executable program:\n" \
                  "$exec_cmd" \
                  "\nProgram not run!"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return -1
  }

  #--Other run arguments
  set browse_mode [string tolower $browse_mode]
  set process_tag $Info($process_name,processTag)
  set logfile_base [file join $ModelProperty(LOG_DIRECTORY,absolute) \
                         $Info($process_name,logfile) ]

  #--Check that logfile can be created!
  if { $browse_mode == "logfile" } {
    
    if { ![file isdirectory $ModelProperty(LOG_DIRECTORY,absolute)] } {

      set dir_msg "for the log file"
      set cre_msg "Press Yes to create the directory and continue\n"
      append cre_msg "Press No to continue and use system shell window as log!"

      set rc [Util::checkAndCreateDirectory3B $ModelProperty(LOG_DIRECTORY,absolute) $dir_msg 1 $cre_msg]
      
      # Do not continue
      if { $rc == "cancel" } {
        return -1
      }
      
      # Continue without logfile
      if { $rc == "no" } {
        set logfile_base ""
        set browse_mode "shell"
      }
    }
  }

  #--Changedir to model's directory (store current)
  set Info(currentDirectory) [pwd]
  Util::cpp_exec setCurrentDirectory $ModelProperty(MODEL_DIRECTORY,absolute)

  # Exec the program
  # ================
  set runnumber [MenuExec::do_elmer_exec \
      $process_name \
      $use_eval $exec_cmd $arg \
      $priority_level \
      $add_to_process_table \
      $show_in_process_list \
      $browse_mode \
      $exec_name \
      $process_tag \
      $logfile_base \
      $process_info]

  #--Restore current directoy
  Util::cpp_exec setCurrentDirectory $Info(currentDirectory)
  
  # Succes, reset flags etc.
  # ========================
  if { $runnumber != -1 } {

    if { $add_to_process_table } {
      set ProcessTable($runnumber,browserDelay) $browser_delay
    }

    # If possible logfile browser needs to know the runnumber 
    # (ex. to be able to update process table etc.)
    if { $browse_parser != "none" && $browse_parser_needs_runnumber } {
      append browse_parser " $runnumber"
    }

    MenuExec::reset_elmer_exec $process_name $process_args $runnumber $browse_mode $browse_parser

    return $runnumber

  # No success
  # ==========
  } else {
    return -1
  }
}


# ElmerPost specific
# ==================
#
proc MenuExec::prepare_ElmerPost_exec {post_file server script_argument} {
  global Info
  upvar $script_argument script

  set rc 1
  set script ""

  # Check server state if something to send
  # =======================================
  if { $Info(Results,$server,mode) != "" } {

    set try_count 10

    set src  [MenuExec::existsTclServer $server $try_count]

    #---Server is not responding!
    if { $src == 0 || $src == -1 } {

      set Info(Results,$server,mode) ""
      
      #--Server just created, but not found!
      if { $src == 0 && $Info(Results,$server,created) } {
	      set msg [list "NOTE: Elmer Post server not found!\n\n" \
                      "We cannot send any update commands to Elmer Post!\n\n" \
                      "Use menus in Elmer Post for updating the results!"]

        set Info(messageIcon) error
        Message::dispOkMessage $msg
	      return 0

      #--Server just created, but not accessible!
      } elseif { $src == -1 && $Info(Results,$server,created) } {
        set msg [list "NOTE: Elmer Post server is not accessible!\n\n" \
                      "We cannot send any update commands to Elmer Post!\n\n"\
                      "Use menus in Elmer Post for updating the results!"]

        set Info(messageIcon) error
        Message::dispOkMessage $msg
	      return 0
      
      #---Server must be created (again)
      } else {
        set Info(Results,$server,mode) ""
      }
    }
  }


  # No server exists yet, we will use command line call
  # ===================================================
  if { $Info(Results,$server,mode) == "" } {
    set Info(Results,$server,mode) "start"
    set rc 0
  }

  # Form proper ep command
  # ======================
  set script [MenuExec::formElmerPostCommand $post_file $Info(Results,$server,mode)]

  #--Set next command mode
  if { $Info(Results,$server,mode) == "start" } {

    #-Windows
    if { $Info(platform) == "windows" } {
      set Info(Results,$server,mode) "init"

    #-Unix (directly from start to refresh!)
    } else {
      set Info(Results,$server,mode) "refresh"
    }

  } elseif { $Info(Results,$server,mode) == "init" } {
    set Info(Results,$server,mode) "refresh"
  }
 
  return $rc
}



#####################
### DO ELMER EXEC ###
#####################

# Do the actual Elmer exec call
# =============================
# Arguments (when not self explaining :-):
# process_name: proc key to Info arries, like "Solver"
# exec_cmd: actual call (path, exec-name etc.)
# arg : arguments for the exec_cmd
# exec_name: name for executable file like "mesh.exe"
# process_tag: a short name for the process table like "Mesh2D", "GebFact" etc.
# process_info: process specific ino like db-name, Mesh-H etc.
#
proc MenuExec::do_elmer_exec {
    process_name
    use_eval 
    exec_cmd arg 
    {priority "NORMAL_PRIORITY"} 
    {add_to_process_table 0} 
    {show_in_process_list 0} 
    {browse_mode "none"} 
    {exec_name "Unknown"}
    {process_tag "Unknown"} 
    {logfile_base ""} 
    {process_info ""} 
    } {

  global Info Model ProcessTable errorCode

  #--If nothing to do!
  if {$exec_cmd == "" } {
    return
  }

  # To be sure!
  if { !$add_to_process_table } {
    set show_in_process_list 0
  }

  incr Info(ElmerRun,runnumber)

  set browse_mode [string tolower $browse_mode]

  # Output directed to a logfile
  # ============================
  if { $browse_mode == "logfile" && $logfile_base != "" } {

    # Try to find a free logfile,runnumber combi
    # Increment runnumber, if necessary
    #
    set success 0
    set counter 0

    while { !$success &&
            $counter < 100
          } {
      incr counter

      # Full name for the logfile!!!
      set testfile $logfile_base.log.$Info(ElmerRun,runnumber)

      # Delete possible old version of the log file
      set msg ""
      set catch_rc [ catch [file delete -force $testfile] msg]

      if {!$catch_rc} {
        set success 1
        set logfile $testfile
      } else {
        incr Info(ElmerRun,runnumber)
      }
    }

    if {!$success} {
      set msg [list "Program $exec_cmd was not run properly, because:\n" \
                    "could not open a logfile for it!"]

      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return -1
    }

  # No logfile
  # ==========
  } else {
    set logfile ""
  }
  
  set pid 0
  set runnumber $Info(ElmerRun,runnumber)
  set starttime [clock seconds]
  set use_log 0


  ###############################
  #--Run the program using catch!
  ###############################

  set msg ""
  set code [ catch {

    #-ElmerPost is run always using exec
    ###########
    # 
    if { $process_name == $Info(Results,processName) } {

      set arg1 [lindex $arg 0] ;# -id
      set arg2 [lindex $arg 1] ;# server name
      set arg3 [lindex $arg 2] ;# -file
      set arg4 [lindex $arg 3] ;# filename
      set arg5 [lindex $arg 4] ;# commands

      # Make easier to print (for messages)
      set arg [join $arg]

      set pid [exec $exec_cmd $arg1 $arg2 $arg3 $arg4 $arg5 & ]
      
    #-WIN32 (other than ElmerPost)
    #######
    } elseif { $Info(platform) == "windows" } {
      
      #-Send process output to shell console
      if { $browse_mode == "shell" } {
        set use_log 0

      #-Send process output to Elmer log file
      } elseif { $browse_mode == "logfile" } {
        set use_log 1

      #-Nothing special with the ouput
      } else {
        set use_log 0
      }

      # Set final logfile mode
      if { $use_log } {
        set logf $logfile
      } else {
        set logf $browse_mode
      }
      
      set data ""
      lappend data $exec_cmd       ;# exe-path
      lappend data $arg            ;# arguments
      lappend data $runnumber      ;# id nbr
      lappend data $exec_name      ;# name
      lappend data $priority       ;# priority level
      lappend data $logf           ;# logfile or none/shell

      Util::cpp_exec processStart $data

    #-UNIX  (other than ElmerPost)
    ######
    } else {

       #-Send process output to shell console
      if { $browse_mode == "shell" } {

        if { $use_eval } {
          set pid [eval exec $exec_cmd $arg & ]
        } else {
          set pid [exec $exec_cmd $arg & ]
        }
      
      #-Send process output to Elmer log file
      } elseif { $browse_mode == "logfile" } {

        if { $use_eval } {
          set pid [eval exec $exec_cmd $arg >>& $logfile & ]
        } else {
          set pid [exec $exec_cmd $arg >>& $logfile & ]
        }

      #-"None", no visible ouput
      } else {

        if { $use_eval } {
          set pid [eval exec $exec_cmd $arg > /dev/null & ]
        } else {
          set pid [exec $exec_cmd $arg > /dev/null & ]
        }
      }
    }

  } msg ]
  ############## End catch Run

  # Display process run message in the message area
  # ===============================================
  if { $Info(displayExecMessages) } {
    Message::showMessage "Process $runnumber: $exec_name $arg" 
  }

  # If something went wrong
  # =======================
  if { $code != 0 } {
    set msg [list "Program $exec_cmd was not run properly, because:\n" \
                  "$msg"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return -1
  }

  # If process is added to process table
  # ====================================
  if { $add_to_process_table } {

    set check_all_done $Info($process_name,checkAllDone)

    ProcessTable::addEntry $pid $runnumber $priority \
                           $show_in_process_list $check_all_done \
                           $process_tag $logfile $starttime $process_info

    if { $show_in_process_list } {
      Process::monitor $runnumber
    }
  }

  return $runnumber
}


########################
### RESET ELMER EXEC ###
########################

# This procedure resets Elemer exec related flags if the call
# to the executable was succeful
#
proc MenuExec::reset_elmer_exec {process_name process_args runnumber {browse_mode ""} {browse_parser "none"} } {
  global Equation Info Model ModelProperty ProcessTable UserSetting

  #---Possibly browse the log file
  #
  if { $browse_mode == "logfile" } {
    MenuExec::browseProcessLog $runnumber $browse_parser
  }

  set set_timestamp 0
  set save_model_file 0

  #---Update flags and timestamps

  # Emf2Db
  # ======
  if { $process_name == $Info(Emf2Db,processName) } {
      set Model(Database,needsUpdate) 0

  # GebhardtFactors
  # ===============
  } elseif { $process_name == $Info(GebhardtFactors,processName) } {
    set set_timestamp 1
    set save_model_file 1
    set Model(GebhardtFactors,needsUpdate,data) 0
    set Model(GebhardtFactors,needsUpdate,mesh) 0

  # Mesh
  # ====
  } elseif { $process_name == $Info(Mesh2D,processName) ||
             $process_name == $Info(Mesh3D,processName) 
           } {

    #-Wait until mesh generator at end
    MenuExec::monitorMesh $runnumber

    if { !$ProcessTable($runnumber,terminated) } {
      tkwait variable ProcessTable($runnumber,terminated)
    }

    set Model(Mesh,generating) 0

    #-If success
    if { $ProcessTable($runnumber,state) == "all done" } {
      
      set set_timestamp 1
      set save_model_file 1
      set Model(MeshFromDB,needsUpdate) 1
      MenuBuild::configureButtonOption loadMesh state disabled

      MenuExec::unloadMesh "--Removing old mesh data from the memory"

      #-Load mesh from DB
      if { $UserSetting(AUTO_LOAD_MESH) } {
        MenuExec::loadMesh $runnumber
      } else {
        MenuBuild::configureFileLoadMeshMenus 1
        MenuBuild::configureButtonOption loadMesh state normal
      }

      #-Set Viewfactors and GebhardtFactors update flags if heat equation
      # mesh was updated
      set eidx [DataField::getFieldProperty Equation HEAT_EQUATION EquationIndex]
      if { $Model(hasDiffuseGrayRadiation) &&
           $Model(currentMeshName) == [MenuExec::getSolverMeshName $eidx 0]
         } {
        
        # If checking is not turend of by the user
        if { $Model(Viewfactors,needsUpdate,mesh) != -1 } { 
          set Model(Viewfactors,needsUpdate,mesh) 1
        }

        # If checking is not turend of by the user
        if { $Model(GebhardtFactors,needsUpdate,mesh) != -1 } {
          set Model(GebhardtFactors,needsUpdate,mesh) 1
        }
      }

      #-Update mesh names
      #set Model(meshNames) [ModelProperty::findModelMeshNames]

    } ;# If mesh generating succesful

  # Results
  # =======
  } elseif { $process_name == $Info(Results,processName) } {

    set server $Info(Results,tclServer)

    #--If ElmerPost should be initialized (in "init"-mode)
    if { $Info(Results,$server,mode) == "init" } {
      MenuExec::prepare_elmer_exec $process_name $process_args
    }
    
  # Solver
  # ======
  } elseif { $process_name == $Info(Solver,processName) } {

    # Warning for zero velocity condition corner bulk elements  
    if { $Model(nofMeshZeroVelocityElements) > 0 } {
      set warning "***WARNING: Mesh contains zero-velocity element(s)!"
      Process::addWarning $runnumber $warning     
    }
      
    # NOTE: Order is important here!    
    set ProcessTable($runnumber,solver,postFiles) [MenuExec::getSolverResultFiles POST_FILE]
    set ProcessTable($runnumber,solver,outputFiles) [MenuExec::getSolverResultFiles OUTPUT_FILE]
    PostFileSelect::updateFileList $runnumber

    MenuBuild::configureButtons results 1
    set set_timestamp 1
    set save_model_file 1
    set Model(Solver,solving) 0

  # Viewfactors
  # ===========
  } elseif { $process_name == $Info(Viewfactors,processName) } {

    set set_timestamp 1
    set save_model_file 1
    set Model(Viewfactors,needsUpdate,data) 0
    set Model(Viewfactors,needsUpdate,mesh) 0
  }

  # We have polled enough the process, we do not need
  # it any more!
  Process::markRemoveable $runnumber

  # Save timestamp data in the model
	#
  if { $set_timestamp } {
    set ts [Util::setCurrentTimestamp]
    set Model($process_name,timestamp) $ts
    Util::cpp_exec readTimestamp [list $process_name $ts]
	}    

  # Save model file
	#
	if { $UserSetting(AUTO_SAVE_MODEL) && $save_model_file } {
    set model_file_only 1
    set show_msg 0
    MenuExec::saveCheckedModelFile $Model(EMF_PATH) $model_file_only $show_msg
  }

}


proc MenuExec::browseProcessLog {runnumber { browse_parser "none" } } {
  global Info ProcessTable

  set logfile $ProcessTable($runnumber,logfile)
  set logchannel $ProcessTable($runnumber,channel)

  set must_exist 0
  set can_edit 0
  set can_delete_cancelled 0

  # Open browser
  FileBrowser::browse .winLogFile $logfile $must_exist $can_edit $can_delete_cancelled \
	                    $browse_parser $logchannel $runnumber
}


proc MenuExec::getElmerPostHeader {post_datafile} {

  if { ![file exist $post_datafile] ||
       ![file readable $post_datafile]
     } {
    return ""
  }

  set msg ""
  set channel_id -1

  if { [catch {set channel_id [open $post_datafile] } msg] } {
    return ""
  }

  set ep_header [gets $channel_id]

  close $channel_id

  return $ep_header
}


# Send command(s) to a remote tcl server
# Give 10 * wait-time test time for the server
# before given up
#
proc MenuExec::sendTclCommands {server cmd {max_try_count 10} {message ""} } {

  if { $message != "" } {
    upvar $message msg
  }

  set try_count 0
  set msg ""

  # Try to send commands
  # ====================
  while {$try_count < $max_try_count} {

    incr try_count

    #--If send was ok
    if { ![catch {send $server [join $cmd] } msg] } {
      return 1
    
    #--Otherewise try again
    } else {
      Util::doWait 500
    }
  }

  # No success!
  # ===========
  return 0
}


# Test if a remote tcl server is alive
#
proc MenuExec::existsTclServer {server {max_try_count 10} } {
  global Info

  set server_state ""
  set try_count 0

  while { $try_count <= $max_try_count } {

    set msg ""
    incr try_count

    # Win32
    # =====
    if { $Info(platform) == "windows" } {
      
      #--Problems
      if { [catch {send $server ""} msg] } {
        set server_state 0

      #--Ok, we can stop
      } else {
        return 1
	    }
	    
    # Unix
    # ====
    } else {

      #--Problems
      if { [catch {send $server ""} msg] } {

        #-"no application named"
        if { [regexp (no) $msg] &&
             [regexp (application) $msg] &&
             [regexp (named) $msg]
           } {
          set server_state 0

        #-"server insecure"
        } elseif { [regexp (server) $msg] &&
                  [regexp (insecure) $msg]
                 } {
          set server_state -1
        }

      #--Ok, we can stop
      } else {
        return 1
	    }
    }

    # Wait a while before the next try
    # ================================
    Util::doWait 100
  }

  return $server_state
}


proc MenuExec::formElmerPostCommand { post_file command_mode} {
  global Info Model Timestep

  set ep_cmd ""

  set ep_header [MenuExec::getElmerPostHeader $post_file]

  # Read nof timesteps
  # ==================
  # 4. item is the nof timesteps
  set nof_ts [lindex $ep_header 3]

  # Read first variable type and name
  # =================================
  set tmp [split $ep_header ":"]
  set tmp_len [llength $tmp]

  #-Last item is the 1. block is the type of the first variable
  set ep_var1_type [lindex [lindex $tmp 0] end]

  #-If we have a single variable
  if { $tmp_len == 2 } {
    set ep_var1_name [join [lindex $tmp 1] ]

  #-All but the last items in the 2. block define the name of the first variable
  } else {
    set end_pos [llength [lindex $tmp 1]]
    incr end_pos -2
    set ep_var1_name [join [lrange [lindex $tmp 1] 0 $end_pos]]
  }

  #-For a vector type, abs-value is always defined, we
  # start with it
  if { "vector" == [string tolower $ep_var1_type] } {
    append ep_var1_name "_abs"
  }
  
  # Timestep info
  # =============

  #---Transient problem: start with the first timestep
  if { $Timestep(SIMULATION_TYPE) == "transient" } {
    set step_info "1 1 1"

  #---Steady state problem: start with the last timestep
  } elseif { $nof_ts != "" } {
    set step_info "$nof_ts $nof_ts 1"
  } else {
    set step_info ""
  }

  # File info
  # =========
  set file_info "readfile "
  append file_info \\\"$post_file\\\"
  append file_info " "

  # Write ep-command
  # ================

  #---Start 
  if { $command_mode == "start" } {
  
    #-Windows start
    if { $Info(platform) == "windows" } {
      set ep_cmd " $step_info;"

    #-Unix start
    } else {
      set ep_cmd " $step_info;"
      append ep_cmd "set DisplayStyle(ColorMesh) 1;"
      append ep_cmd "set MeshStyle  1;"
      append ep_cmd "set MeshColor  $ep_var1_name;"
      append ep_cmd "set DisplayStyle(ColorScale) 1;"
      append ep_cmd "set ColorScaleStyle  1;"
      append ep_cmd "set ColorScaleColor  $ep_var1_name;"
      append ep_cmd "set ColorScaleX -0.8;"
      append ep_cmd "set ColorScaleY -0.95;"
      append ep_cmd "set ColorScaleDecimals 4;"
      append ep_cmd "UpdateObject;"
      append ep_cmd "mesh_edit;"      ;# Open ColorMeshEdit-panel
      #append ep_cmd "colscale_edit;" ;# Open ColorScaleEdit-panel
      append ep_cmd "time_set_loop;"  ;# Update time value in TimeStep-loop panel (if possibly open) 
      append ep_cmd "display;"
    }

  #---Init
  } elseif { $command_mode == "init" } {

    #-Windows init
    if { $Info(platform) == "windows" } {
      append ep_cmd "set DisplayStyle(ColorMesh) 1;"
      append ep_cmd "set MeshStyle  1;"
      append ep_cmd "set MeshColor  $ep_var1_name;"
      append ep_cmd "set DisplayStyle(ColorScale) 1;"
      append ep_cmd "set ColorScaleStyle  1;"
      append ep_cmd "set ColorScaleColor  $ep_var1_name;"
      append ep_cmd "set ColorScaleX -0.8;"
      append ep_cmd "set ColorScaleY -0.95;"
      append ep_cmd "set ColorScaleDecimals 4;"
      append ep_cmd "UpdateObject;"
      append ep_cmd "mesh_edit;"      ;# Open ColorMeshEdit-panel
      #append ep_cmd "colscale_edit;" ;# Open ColorScaleEdit-panel
      append ep_cmd "time_set_loop;"  ;# Update time value in TimeStep-loop panel (if possibly open) 
      append ep_cmd "display;"

    #-Unix init (do nothing!)
    } else {
      set ep_cmd ""
    }

  #---Refresh
  } elseif { $command_mode == "refresh" } {
    set ep_cmd $file_info
    append ep_cmd "$step_info;"
    append ep_cmd "time_set_loop;"  ;# Update time value in TimeStep-loop panel (if possibly open) 
    append ep_cmd "display;"
  }

  # Return command
  # ==============
  return $ep_cmd
}


# Elmerpost commands delivered using a file
# (Not in use)
#
proc MenuExec::formElmerPostCommandFile {post_file} {
  global Model Timestep

  set msg ""
  set channel_id -1

  set cmd_file [file join $ModelProperty(LOG_DIRECTORY,absolute) "ELMERPOST_COMMAND_FILE"]

  if { [catch { set channel_id [open $cmd_file "w"] } msg] } {
    return ""
  }

  set ep_header [MenuExec::getElmerPostHeader $post_file]

  #---Read nof timesteps
  # 4. item is the nof timesteps
  set nof_ts [lindex $ep_header 3]

  #---Read first variable type and name
  set tmp [split $ep_header ":"]
  set tmp_len [llength $tmp]

  # last item is the 1. block is the type of the first variable
  set ep_var1_type [lindex [lindex $tmp 0] end]

  # If we have a single variable
  if { $tmp_len == 2 } {
    set ep_var1_name [join [lindex $tmp 1] ]

  # All but the last items in the 2. block define the name of the first variable
  } else {
    set end_pos [llength [lindex $tmp 1]]
    incr end_pos -2
    set ep_var1_name [join [lrange [lindex $tmp 1] 0 $end_pos]]
  }

  # For a vector type, abs-value is always defined, we
  # start with it
  if { "vector" == [string tolower $ep_var1_type] } {
    append ep_var1_name "_abs"
  }

  #---Tansient problem: start with the first timestep
  if { $Timestep(SIMULATION_TYPE) == "transient" } {
    set step_info "1 1 1"

  #---Steady state problem: start with the last timestep
  } elseif { $nof_ts != "" } {
    set step_info "$nof_ts $nof_ts 1"
  } else {
    set step_info ""
  }

  #---Write ep-command file
  puts $channel_id "readfile \"$post_file\" $step_info;"
  puts $channel_id "set DisplayStyle(ColorMesh) 1;"
  puts $channel_id "set MeshStyle  1;"
  puts $channel_id "set MeshColor  $ep_var1_name;"
  puts $channel_id "set DisplayStyle(ColorScale) 1;"
  puts $channel_id "set ColorScaleStyle  1;"
  puts $channel_id "set ColorScaleColor  $ep_var1_name;"
  puts $channel_id "set ColorScaleX -0.8;"
  puts $channel_id "set ColorScaleY -0.95;"
  puts $channel_id "set ColorScaleDecimals 4;"
  puts $channel_id "UpdateObject;"
  puts $channel_id "mesh_edit;"
  puts $channel_id "time-set_loop;"
  #puts $channel_id "colscale_edit;"

  close $channel_id

  return $cmd_file
}



#####################################################################
#
#              MODEL DIRECTORY AND FILE utilities
#
#####################################################################

# Read color table file (via cpp)
#
proc MenuExec::readColorFile { {type "front"} } {
  global Info Model ModelProperty
  
  set path ""
  set fn "rgb.txt"
  set places {ELMER_FRONT_BUILD_LIB ELMER_FRONT_INSTALL_LIB}

  #-Search in Front-lib
  #
  if { $type == "front" } {

    if { [info exists Info(ELMER_FRONT_HOME)] &&
         $Info(ELMER_FRONT_HOME) != ""
       } {
      set path [file join $Info(ELMER_FRONT_HOME) lib]

    }
  } elseif {$type == "elmer"} {
	foreach place $places {

	    if { [info exists Info($place)] && $Info($place) != "" } {
		set path [file join $Info($place) $fn]
		if { [file exist $path] } {
		    set path $Info($place)
		}
	    }
	}
  #-Search in elmer-user's home directory
  #
  } elseif { $type == "user" } {
    set path $Info(ELMER_USER_HOME)

  #-Search in model file directory
  #
  } elseif { $type == "model" } {
    
    if { $ModelProperty(MODEL_DIRECTORY) != "" } {
      set path $ModelProperty(MODEL_DIRECTORY,absolute)
    }
  }

  # Read color table file if file found
  #
  if { $path != "" } {
    set path [file join $path $fn]
    if { [file exist $path] } {
      Util::cpp_exec readColorFile $path
    }
  }

}


# Check if a Matc-file is found in a predefined directory and read the file
# via cpp call
#
# NOTE: File name must match exactly to "matc.txt" !!!
#
proc MenuExec::readMatcFile { {type "front"} } {
  global Info Model ModelProperty
  
  set path ""
  set fn "matc.txt"

  #-Search in elmer home directory
  #
  if { $type == "elmer" } {
    set path $Info(ELMER_FRONT_HOME)

  #-Search in elmer-user's home directory
  #
  } elseif { $type == "user" } {
    set path $Info(ELMER_USER_HOME)

  #-Search in model file directory
  #
  } elseif { $type == "model" } {
    
    if { $ModelProperty(MODEL_DIRECTORY) != "" } {
      set path $ModelProperty(MODEL_DIRECTORY,absolute)
    }
  }

  # Read Matc file if file found
  #
  if { $path != "" } {
    set path [file join $path $fn]

    if { [file exist $path] } {
      Util::cpp_exec readMatcFile $path
    }
  }
}


# Get Mesh2D input file name (mif-file) for the model
#
proc MenuExec::getMeshInputFileName {} {
  global ModelProperty

    set mif_name $ModelProperty(MODEL_NAME)
    append mif_name ".mif"

    return $mif_name
}


# Get filename for a mesh related file (like elmerpost file) for the model
#
proc MenuExec::getMeshRelatedFiles { mesh_names file_var {absolute 1} } {
  global Solver

  set files ""

  foreach mesh_name $mesh_names {
    
    set fname [MenuExec::getMeshRelatedFile $mesh_name $file_var]

    # Absolute path
    if { $absolute } {
      lappend files $fname

    # Meshname + filename (like Mesh1/CircleInBox.ep)
    } else {
      set fn [file tail $fname]
      set mn [lindex [file split [file dirname $fname] ] end]
      lappend files [file join $mn $fn]
    }
  }

  return $files
}


# Absolute path for mesh related files like
# GebhardtFactors and ViewFactors
# NOTE: Datafile is suppoused to have all relevant data
# in variable field, no paramter line i schecked
# (perhaps this should be changed to work via parameter
# ids?)
#
proc MenuExec::getMeshRelatedFile { mesh_name file_var} {
  global Datafile Info ModelProperty

  if { $mesh_name == "" || 
       ![info exists Datafile($file_var)]
     } {
    return ""
  }

  set fname $Datafile($file_var)


  # If fname is not relative, we use it as it is
  if { "relative" != [file pathtype $fname] } {
    return $fname
  }

  # If results' directory given,use it as thebase, add the meshname
  if { $ModelProperty(RESULTS_DIRECTORY) != "" } {
    return [file join $ModelProperty(RESULTS_DIRECTORY,absolute) \
                      $mesh_name \
                      $fname]
  }

  # Default solver output directory, add the meshname
  return [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                    $Info(meshDirectoryName) \
                    $mesh_name \
                    $fname ]
  
}


# Get absolute path name for the solver's mesh file' header
# Solver defined by the system and subsystem indices
#
proc MenuExec::getSolverMeshHeaderFile {system_index subsystem_index} {
  global Info ModelProperty

  set mname [MenuExec::getSolverMeshName $system_index $subsystem_index]

  if { $mname == "" } {
    return ""
  }

  # Add paths
  return [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                    $Info(meshDirectoryName) \
                    $mname \
                    "mesh.header"]

}


# Get mesh file name for the solver
# Solver define by the system and subsystem indices
# NOTE: This is in fact a directory name under MESHDIR, but
# no path information is return!
#
proc MenuExec::getSolverMeshName {system_index subsystem_index} {
  global Info Solver SolverSystem

  set sid [lindex $SolverSystem($system_index,ids) $subsystem_index]

  if { $sid == $Info(NO_INDEX) || $sid == "" } {
    return ""
  }

  return  [DataField::getFieldValue Solver $sid MESH]
}


# Absolute path name for a solver's mesh directory
# Solver identified by solver id
#
proc MenuExec::getSolverMeshDirectory { sid } {
  global Info ModelProperty Solver

  if { $sid == $Info(NO_INDEX) || $solver_id == "" } {
    return ""
  }

  # Get solver mesh name
  set mesh_name [DataField::getFieldValue Solver $sid MESH]

  # Add paths
  return [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                    $Info(meshDirectoryName) \
                    $mesh_name ]

}


proc MenuExec::getSolverMeshDirectories {} {
  global Solver

  set dirs ""

  foreach sid $Solver(activeIds) {
    lappend dirs [MenuExec::getSolverMeshDirectory $id]
  }

  return $dirs
}


proc MenuExec::getSolverResultFileEquations {result_file} {
  global Equation Solver 

  set equations ""

  set mesh [lindex [file split [file dirname $result_file] ] end]

  foreach id $Solver(activeIds) {

    set smname [DataField::getFieldValue Solver $id MESH]

    if { $smname == "" || $smname == $mesh } {

      set eq_var [DataField::fieldNameSifToGui [DataField::getFieldValue Solver $id EQUATION]]

      if { ![info exists Equation($eq_var)] } {
        continue
      }

      set dsp_name [DataField::getFieldProperty Equation $eq_var StatusLineName]

      append equations $dsp_name]
      append equations "  "
    }
  }

  set equations [string trim $equations]

  return [list $equations]
}


proc MenuExec::getSolverResultFiles { file_var {absolute 1} } {
  global Solver

  set files ""

  foreach id $Solver(activeIds) {
    
    set fname [MenuExec::getSolverResultFile $id $file_var]

    # Absolute path
    if { $absolute } {
      lappend files $fname

    # Meshname + filename (like Mesh1/CircleInBox.ep)
    } else {
      set fn [file tail $fname]
      set mn [lindex [file split [file dirname $fname] ] end]
      lappend files [file join $mn $fn]
    }
  }

  return $files
}


# Absolute path name for solver result files for the solver with "id"
# Use solver-id to get the propeer mesh-name
# Used fex. for Output and Post result files
#
proc MenuExec::getSolverResultFile { solver_id file_var } {
  global Info Datafile ModelProperty Solver

  if { ![info exists Datafile($file_var)] } {
    return ""
  }

  set fname $Datafile($file_var)

  #--If "generic" query, return just filename
  #
  if { $solver_id == $Info(NO_INDEX) || $solver_id == "" } {
    return $fname
  }

  #--If "generic single solver" query, return value for the first
  #  active Solver
  if { $solver_id == 0 } {
    set solver_id [lindex $Solver(activeIds) 0]

    return [MenuExec::getSolverResultFile $solver_id $file_var]
  }

  #--If fname is not relative, we use it as it is
  if { "relative" != [file pathtype $fname] } {
    return $fname
  }
  
  #--Get current base directory for the solver and mesh
  set dir_name [MenuExec::getSolverResultDirectory $solver_id]

  #--Check that directory exists
  set rc [Util::checkAndCreateDirectory2B $dir_name "Solver result files directory"]

  if { $rc != "ok" } {
    return ""
  }

  return [file join $dir_name $fname]
}


# Absolute dirdctory path for solver result files for the solver with "id"
# Use solver-id to get the propeer mesh-name
# Used fex. for Output and Post result files
#
proc MenuExec::getSolverResultDirectory { solver_id  } {
  global Info ModelProperty Solver

  #--If "generic single solver" query, return value for the first
  #  active Solver
  if { $solver_id == 0 || $solver_id == $Info(NO_INDEX) || $solver_id == "" } {
    set solver_id [lindex $Solver(activeIds) 0]
    return [MenuExec::getSolverResultDirectory $solver_id]
  }

  # Get solver mesh file name (directory in fact!)
  set mesh_name [DataField::getFieldValue Solver $solver_id MESH]
 
  if { $mesh_name == "" } {
     set mesh_name [Solver::getDefaultMesh]
     }

  #--If result's directory given, use it as the base, add the meshname
  if { $ModelProperty(RESULTS_DIRECTORY) != "" } {
    return [file join $ModelProperty(RESULTS_DIRECTORY,absolute) \
                      $mesh_name]
  }

  # Default solver output directory, add the meshname
  return [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                    $Info(meshDirectoryName) \
                    $mesh_name]

}


# Full name for the solver input file
# NOTE: simple filename, model-directory is added
#
proc MenuExec::getSolverInputFile {} {
  global Datafile ModelProperty

  return [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                    $Datafile(SOLVER_INPUT_FILE)]
}


# Full name for the solver reload file
#
proc MenuExec::getSolverReloadFile {} {
  global SolverControl ModelProperty

  set reload_file [DataField::getFieldValue SolverControl 1 RELOAD_INPUT_FILE]

  if { $reload_file == "" } {
    return ""
  }

  return [file join $ModelProperty(MODEL_DIRECTORY,absolute) $reload_file]
}


# Check that directory exists, if not, open
# a display for entry
#
proc MenuExec::checkDirectory {dir title label} {
  global Info

  if { [file isdirectory $dir] } {
    return $dir
  }
  
  global tmp_dir
  set tmp_dir ""

  if { $dir != "" } {
    set try_dir [file dirname $dir]
  } else {
    set try_dir $dir
  }
  
  incr Info(directoryEntryId)
  set wname ".w$Info(directoryEntryId)"

  Widget::dispDirectoryEntry $wname tmp_dir $try_dir $title $label
  tkwait variable tmp_dir
  
  set res_dir $tmp_dir

  unset tmp_dir

  return $res_dir
}


proc MenuExec::getCurrentDirectory {} {
  return [pwd]
}


proc MenuExec::setWorkingDirectory { {initial_value ""} } {
  global Info Model

  if { $initial_value == "" } {
    set initial_value $Info(workingDirectory)
  }

  set title "Set working directory"
  set label "Enter directory name"

  global wkd ;# Reference var must be global in dispDirectoryEntry!

  set wkd ""

  incr Info(directoryEntryId)
  set wname ".w$Info(directoryEntryId)"
  
  Widget::dispDirectoryEntry $wname wkd $Info(workingDirectory) $title $label 
  tkwait variable wkd

  if { $wkd != "!cancel!"  } {
    
    #-If ok value   
    if { [file isdirectory $wkd] } {

      set Info(workingDirectory) $wkd
      Message::showMessage "Working directory: $Info(workingDirectory)"

    #-Invalid directory name
    } else {

      set msg [list "A non existing directory name!\n\n$wkd\n\n" \
                     $Info(anywayOk) ]

      set Info(messageIcon) error
      set result [ Message::dispOkCancelMessage $msg ]

      # Try again
      if { $result != "cancel" } {
        MenuExec::setWorkingDirectory $wkd

      # Close entry
      } else {
        unset wkd
        return
      }
    }

    # Put current working directory as default model directory (if not defined)
    ModelProperty::setDefaultDirectoryValue MODEL_DIRECTORY $Info(workingDirectory)

    set Model(CAD_OPEN_DIR) $Info(workingDirectory)
    set Model(MESH_OPEN_DIR) $Info(workingDirectory)
    set Model(EMF_OPEN_DIR) $Info(workingDirectory)

    # Set pwd
    Util::cpp_exec setCurrentDirectory $Info(workingDirectory)
  }

  unset wkd
}


# Apply user settings values to problem directories
#
proc MenuExec::applySettingsToCaseDirectories { {force_to_settings 0 } {show_usage_message 1} } {
  global Info Model ModelProperty UserSetting

  # Model directory
  # ===============
  if { $force_to_settings ||
       $ModelProperty(MODEL_NAME) == "" ||
       $ModelProperty(MODEL_DIRECTORY) == ""
     } {
    set ModelProperty(MODEL_DIRECTORY) $UserSetting(DEFAULT_MODEL_DIRECTORY)
  }
  
  set ModelProperty(MODEL_DIRECTORY,absolute) [Util::makePathAbsolute $ModelProperty(MODEL_DIRECTORY)]

  # Default for log directory (is always needed!)
  # =========================
  set def_val $Info(defaultLogDirectoryName)
  ModelProperty::setDefaultDirectoryValue LOG_DIRECTORY $def_val

  # Other problem directories
  # =========================
  set dir_names { INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY }
  foreach vn $dir_names {
    
    # Store settings values to settings-variable
    set ModelProperty($vn,settings) $UserSetting(DEFAULT_$vn)
    
    # If value is missing or should force anyway
    # 
    if { $force_to_settings || 
         ( $ModelProperty($vn) == ""  &&
           $UserSetting(DEFAULT_$vn) != ""
         )
       } {
      set ModelProperty($vn) $UserSetting(DEFAULT_$vn)
      set ModelProperty($vn,valueType) [ModelProperty::valueTypeMenuVar2Key settings]
    }
  }

  # =================
  #-Check new values
  ModelProperty::checkCaseDirectoryValues

  # Show info
  # =========
  if { $show_usage_message } {
    MenuExec::showCaseDirectoriesMessage
  }
}


# Show current values for problem directories
# NOTE: If model name is not yet defined, we show
# only models-directory value
#
proc MenuExec::showCaseDirectoriesMessage {} {
  global Info ModelProperty Model

  set clr $Info(remMsgColor)

  Message::showMessage "Model dir = $ModelProperty(MODEL_DIRECTORY,absolute)" $clr
  
  if { $Model(CASE_NAME) != "" } {
    Message::showMessage "Inc path = $ModelProperty(INCLUDE_PATH,absolute)" $clr
    Message::showMessage "Res dir  = $ModelProperty(RESULTS_DIRECTORY,absolute)" $clr
    Message::showMessage "Log dir  = $ModelProperty(LOG_DIRECTORY,absolute)" $clr
    Message::showMessage ""
  }
}



#####################################################################
#
#                    WINDOW MENU procs
#
#####################################################################

proc MenuExec::closeAllWindows  {} {
  global Info

  foreach row $Info(windowList) {

    set arr [lindex $row 0]
    set win [winfo atomname [lindex $row 1]]

    set modified 0

    upvar #0 $arr theArray

    if { ( [info exists theArray(dataChanged)] && $theArray(dataChanged) )  ||
         ( [info exists theArray(dataModified)] && $theArray(dataModified) )
       } {
      set modified 1
    }

    catch {MenuExec::closePanel $arr $win $modified}
  }
}


proc MenuExec::closeUnmodifiedWindows  {} {
  global Info

  foreach row $Info(windowList) {

    set arr [lindex $row 0]
    set win [winfo atomname [lindex $row 1]]

    upvar #0 $arr theArray

    set modified 0

    if { ( [info exist theArray(dataChanged)] && $theArray(dataChanged) )  ||
         ( [info exist theArray(dataModified)] && $theArray(dataModified) )
       } {
      set modified 1
    }

    if { !$modified} {
      catch {MenuExec::closePanel $arr $win $modified}
    }
  }
}


proc MenuExec::closePanel {globArray win is_modified} {
  global Common

  set is_standard_panel 0
    
  #-Check if panel is one of the standard panels, they
  # need special handling for the cancel
  foreach arr $Common(standardPanelArries) {

    if { $arr == $globArray } {
      set is_standard_panel 1
      break
    }
  }

  #-Standard panel cancel
  if {$is_standard_panel} {

    set tmp $Common(currentStandardPanel)
    set Common(currentStandardPanel) $globArray

    if {$is_modified} {
      StdPanelExec::panelCancel $win $globArray

    } else {
      StdPanelExec::panelOk $globArray
    }

    set Common(currentStandardPanel) $tmp

  #-Other cancel
  } else {

    if {$is_modified} {
      execNsProc $globArray panelCancel $win

    } else {
      execNsProc $globArray panelOk $win
    }

  }

}


#####################################################################
#
#                    HELP MENU procs
#
#####################################################################

proc MenuExec::helpIndex {} {

  Message::dispOkMessage "Sorry, not yet implemented!"
}


proc MenuExec::helpContents {} {

  Message::dispOkMessage "Sorry, not yet implemented!"
}


# ecif_tk_procsMenu.tcl
# ********************
