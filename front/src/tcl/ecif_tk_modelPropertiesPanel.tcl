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
#Module:    ecif_tk_modelProperties.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting general model properties.
#
#************************************************************************

#--Edit model file paths
#
proc ModelProperty::openPanel {} {
  global Model ModelProperty Info

  set w $ModelProperty(winName)
  set wgeom $ModelProperty(winGeometry)

  set id [winfo atom $w]
  set ModelProperty(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow ModelProperty $id $ModelProperty(winTitle) $wgeom] } {
    return
  }  

  # Update problem directory state variables
  ModelProperty::updateCaseDirectoryVariables

  # Initialize field states
  ModelProperty::resetFieldStates
	
	set ModelProperty(dataChanged) 0
  set ModelProperty(dataModified) 0

  toplevel $w 
  focus $w
  set this $w

  #--Window properties
  wm title $w $ModelProperty(winTitle)
  wm geometry $w $wgeom

  #-Init vars if they don't have values
  Panel::initFields ModelProperty

  Panel::backupFields ModelProperty

  #-Init extra panel vars
  set flds {
    MODEL_DIRECTORY,absolute
    INCLUDE_PATH,valueType
    INCLUDE_PATH,model,save
    RESULTS_DIRECTORY,valueType
    RESULTS_DIRECTORY,model,save
    LOG_DIRECTORY,valueType
    LOG_DIRECTORY,model,save
  }
  Panel::initFields ModelProperty $flds

  #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)

  #---Outer frames
  set f1 [frame $w.f1] ;#File path entries
  set f11 [frame $w.f1.ef1] ;# Model name
  set f12 [frame $w.f1.ef2] ;# Problem name
  set f13 [frame $w.f1.ef3] ;# Result mesh name
  set f14 [frame $w.f1.ef4] ;# Model description
  set f15 [frame $w.f1.ef5] ;# Problem description
  set f16 [frame $w.f1.ef6] ;# Model directory
  set f17 [frame $w.f1.ef7] ;# Include path
  set f18 [frame $w.f1.ef8] ;# Results directory
  set f19 [frame $w.f1.ef9] ;# Log files directory
  set f2 [frame $w.f2] ;#Apply+Ok+cancel buttons frame
  #
  set lbl_wid  17
  set desc_wid 53
  set dir_wid  45
  set name_wid 32

  #---Model name
  set fld MODEL_NAME
  set f $f11

  label $f.l  -text "Model name:" -width $lbl_wid -anchor w
  set ModelProperty(allWidgets,label,$fld) $f.l

  entry $f.e  -textvariable ModelProperty($fld) -width $name_wid \
              -font $Info(entryFont)
  focus $f.e 

  #---Problem name
  set fld PROBLEM_NAME
  set f $f12

  label $f.l  -text "Problem name:" -width $lbl_wid -anchor w
  set ModelProperty(allWidgets,label,$fld) $f.l

  entry $f.e  -textvariable ModelProperty($fld) -width $name_wid \
              -font $Info(entryFont)
 
  # If model name given, give focus to problem name
  if { $ModelProperty(MODEL_NAME) != "" } {
    focus $f.e 
  }


  #---Result mesh name
  set fld RESULT_MESH_NAME
  set f $f13

  label $f.l  -text "Result mesh name:" -width $lbl_wid -anchor w
  set ModelProperty(allWidgets,label,$fld) $f.l

  entry $f.e  -textvariable ModelProperty($fld) -width $name_wid \
              -font $Info(entryFont)
  

  #---Model description
  set fld MODEL_DESCRIPTION
  set f $f14
  set ModelProperty(allWidgets,label,$fld) $f.l

  label $f.l  -text "Model description:" -width $lbl_wid -anchor w

  text $f.e   -width $desc_wid -height 2
  set ModelProperty($fld,widget) $f.e
  $ModelProperty($fld,widget) insert end $ModelProperty($fld)

  #---Problem description
  set fld PROBLEM_DESCRIPTION
  set f $f15

  label $f.l  -text "Problem description:" -width $lbl_wid -anchor w
  set ModelProperty(allWidgets,label,$fld) $f.l

  text $f.e   -width $desc_wid -height 2 
  set ModelProperty($fld,widget) $f.e
  $ModelProperty($fld,widget) insert end $ModelProperty($fld)

  #---Model directory
  set fld MODEL_DIRECTORY
  set f $f16

  label $f.l  -text "Model directory:" -width $lbl_wid -anchor w
  set ModelProperty(allWidgets,label,$fld) $f.l

  entry $f.e  -textvariable ModelProperty($fld) -width $dir_wid \
              -font $Info(entryFont)
  button $f.b -text ... -width 1 \
                -command "Util::fillDirectoryEntry $f.e -panel ModelProperty" 


  #---Problem directories
  # Include Path
  # Results Directory
  # Log Directory
  set flds  { INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY }
  set lbls  { "Include path:" "Results directory:" "Log directory:"}
  set frms  { f17 f18 f19 }
  set apnds { 1 0 0 }

  foreach frm $frms fld $flds lbl $lbls apnd $apnds {
    set f [set $frm]

    #-Label
    label $f.l  -text $lbl -width $lbl_wid -anchor w
    set ModelProperty(allWidgets,label,$fld) $f.l

    #-Directory entry
    entry $f.e  -textvariable ModelProperty($fld) -width $dir_wid \
                -font $Info(entryFont)
    
    #-Browse button
    button $f.b -text ... -width 1 \
                -command "Util::fillDirectoryEntry $f.e -panel ModelProperty -append $apnd"

    #-"Save in session" checkButton
    #set wdg [checkbutton $f.cb_session_save -variable ModelProperty($fld,session,save) -text "Save in\nsession"]

    #-"Save in model" checkb
    set wdg [checkbutton $f.cb_model_save -variable ModelProperty($fld,model,save) -text "Save in\nmodel"]
    set ModelProperty(allWidgets,$fld,model,save) $wdg
      
    #-Value type optionMenu
    set wdg $f.vtype
    set ModelProperty(allWidgets,$fld,vtype) $wdg
    set om [Panel::createOptionMenuWidget $wdg \
                   ModelProperty($fld,valueType) $ModelProperty(valueTypes)]

    $wdg configure -indicatoron 0 -width 6 -activebackground $Info(optionMenuBg)

    set check_proc "\"ModelProperty::applyDirectoryValueType $fld\""
    bind $om <<MenuSelect>> "PanelCheck::checkOptionMenu $om $wdg ModelProperty $fld,valueType $check_proc {%s}"
  }
 
  #---Buttons

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $f2.ok -text OK -command "ModelProperty::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "ModelProperty::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f2.apply -text Apply -command "ModelProperty::panelApply" \
                               -state $ap]

  focus $ok_btn
  set ModelProperty(applyButton)  $ap_btn
  set ModelProperty(cancelButton) $cn_btn
  
  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1
  
  # Pack and set bindings for fields
  set ids {
    11 12 
    14 15
    16
    17 18 19
  }

  set vars {
    MODEL_NAME PROBLEM_NAME
    MODEL_DESCRIPTION  PROBLEM_DESCRIPTION
    MODEL_DIRECTORY
    INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY
  }

  foreach i $ids var $vars {

    #-Frame
    set f [set f$i]

    #- Label and entry
    pack $f.l $f.e -side left -anchor w

    #-Browse button
    if {[winfo exists $f.b]} {
      pack $f.b -side left -anchor w -padx $fpx1
    }

    #-Value type option menu
    set fld_hn $var.VALUE_TYPE
    set wdg $f.vtype

    if {[winfo exists $wdg]} {
      pack $wdg -side left -anchor w  -padx $fpx1
      bind $wdg <ButtonPress-3> "Widget::genericButton-3 ModelProperty $fld_hn $wdg"
    } 

    #-Save in session checkbox
    set fld $var,session,save
    set fld_hn $var,SAVE_IN_SESSION
    set wdg $f.cb_session_save

    if {[winfo exists $wdg]} {
      pack $wdg -side left -anchor w
      bind  $wdg <ButtonPress-1> "Panel::panelDataChanged 1 ModelProperty $wdg {%A %K}"
      Widget::setCheckBoxBindings non_standard ModelProperty $fld $wdg
      bind $wdg <ButtonPress-3> "Widget::genericButton-3 ModelProperty $fld_hn $wdg"

    } 

    #-Save in model checkbox
    set fld $var,model,save
    set fld_hn $var,SAVE_IN_MODEL
    set wdg $f.cb_model_save

    if {[winfo exists $wdg]} {
      pack $wdg -side left -anchor w  -padx $fpx1
      bind  $wdg <ButtonPress-1> "Panel::panelDataChanged 1 ModelProperty $wdg {%A %K}"
      Widget::setCheckBoxBindings non_standard ModelProperty $fld $wdg
      bind $wdg <ButtonPress-3> "Widget::genericButton-3 ModelProperty $fld_hn $wdg"
    } 

    pack $f -side top -pady $fpy2 -anchor w

    # Entry bindings
    set wdg $f.e
    set ModelProperty(allWidgets,$var) $wdg
    bind  $wdg <KeyRelease> "+Panel::panelDataChanged 1 ModelProperty $wdg {%A %K}"
    bind $wdg <KeyRelease> "+Widget::entryKeyRelease ModelProperty $var $wdg {%A %K}"
    bind $wdg <KeyPress-Return> "+ModelProperty::checkEntry $var"
    Widget::setEntryBindings non_standard ModelProperty $var $wdg

    set ModelProperty($var,err) 0
    set ModelProperty($var,mod) 0
    set ModelProperty($var,prev) $ModelProperty($var)
  }

  pack $f1 -side top -expand 1 -anchor w -padx $fpx2 -pady $fpy2
  pack $f2 -side top -expand 1 -padx $fpx2 -pady $fpy2

  # Set field label bindings for right-button help
  Widget::setLabelBindings ModelProperty

	# Check model name data
	# =====================
	#
	if { $Model(EMF_FILE) == "" &&
	     $ModelProperty(MODEL_NAME) != ""
		 } {
		ModelProperty::checkModelNames
		ModelProperty::panelCheck
	}

}


# Set initial values
#
proc ModelProperty::initCaseDirectoryVariables {} {
  global ModelProperty

  set ModelProperty(MODEL_DIRECTORY) ""
  set ModelProperty(MODEL_DIRECTORY,absolute) ""

  set ModelProperty(CURRENT_MESH_DIRECTORY) ""
  set ModelProperty(CURRENT_MESH_DIRECTORY,absolute) ""

  #set ModelProperty(INCLUDE_PATH) ""
  set ModelProperty(INCLUDE_PATH,valueType) [ModelProperty::valueTypeMenuVar2Key session]
  set ModelProperty(INCLUDE_PATH,absolute) ""
  set ModelProperty(INCLUDE_PATH,settings) ""
  set ModelProperty(INCLUDE_PATH,session) ""
  set ModelProperty(INCLUDE_PATH,model) ""
  set ModelProperty(INCLUDE_PATH,session,save) 0
  set ModelProperty(INCLUDE_PATH,model,save) 0

  #set ModelProperty(RESULTS_DIRECTORY)  ""
  set ModelProperty(RESULTS_DIRECTORY,valueType) [ModelProperty::valueTypeMenuVar2Key session]
  set ModelProperty(RESULTS_DIRECTORY,absolute)  ""
  set ModelProperty(RESULTS_DIRECTORY,settings) ""
  set ModelProperty(RESULTS_DIRECTORY,session) ""
  set ModelProperty(RESULTS_DIRECTORY,model) ""
  set ModelProperty(RESULTS_DIRECTORY,session,save) 0
  set ModelProperty(RESULTS_DIRECTORY,model,save) 0

  #set ModelProperty(LOG_DIRECTORY)  ""
  set ModelProperty(LOG_DIRECTORY,valueType) [ModelProperty::valueTypeMenuVar2Key session]
  set ModelProperty(LOG_DIRECTORY,absolute)  ""
  set ModelProperty(LOG_DIRECTORY,settings) ""
  set ModelProperty(LOG_DIRECTORY,session) ""
  set ModelProperty(LOG_DIRECTORY,model) ""
  set ModelProperty(LOG_DIRECTORY,session,save) 0
  set ModelProperty(LOG_DIRECTORY,model,save) 0
}


# Update case directory values
#
proc ModelProperty::updateCaseDirectoryVariables {} {
  global ModelProperty

  set dirvars {INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY}

  foreach fld $dirvars {

    if { $ModelProperty($fld,model,save) } {
      set ModelProperty($fld,valueType) Model
      set ModelProperty($fld,model) $ModelProperty($fld)
    }
  }
}


# Set case name (= model name + problem name)
#
proc ModelProperty::setModelCaseName {} {
  global Model ModelProperty

  #-Add problem name to model-file name, if defined
  if { $ModelProperty(PROBLEM_NAME) != "" } {

    set Model(CASE_NAME) $ModelProperty(MODEL_NAME)
    append Model(CASE_NAME) "."
    append Model(CASE_NAME) $ModelProperty(PROBLEM_NAME)

  } else {
    set Model(CASE_NAME) $ModelProperty(MODEL_NAME)
  }
}


# Reset panel field states
#
proc ModelProperty::resetFieldStates {} {
  global ModelProperty

  # Reset field states
  foreach fld $ModelProperty(allFields) {

    if { !$ModelProperty($fld,trg) } {
      continue
    }

    set ModelProperty($fld,prev) $ModelProperty($fld)
    set ModelProperty($fld,err) 0
    set ModelProperty($fld,mod) 0
    Widget::setFieldStatus ModelProperty $fld
  }

  # Dir "save in model" checkbuttons and value type option menu
  foreach fld {INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY} {

    #set ModelProperty($fld,session,save,prev) $ModelProperty($fld,session,save)
    #set ModelProperty($fld,session,save,err) 0
    #set ModelProperty($fld,session,save,mod) 0

    set ModelProperty($fld,model,save,prev) $ModelProperty($fld,model,save)
    set ModelProperty($fld,model,save,err) 0
    set ModelProperty($fld,model,save,mod) 0
    Widget::setFieldStatus ModelProperty $fld,model,save CheckBox

    set ModelProperty($fld,valueType,prev) $ModelProperty($fld,valueType) 
    set ModelProperty($fld,vtype,err) 0
    set ModelProperty($fld,vtype,mod) 0
    Widget::setFieldStatus ModelProperty $fld,vtype OptionMenu

  }
}


# Store panel data in cpp-side
#
proc ModelProperty::panelSave {} {
  global Model ModelProperty Datafile

  # Store old values
  Panel::backupFields ModelProperty

  if { $ModelProperty(MODEL_DIRECTORY,mod)    ||
       $ModelProperty(INCLUDE_PATH,mod)       ||
       $ModelProperty(RESULTS_DIRECTORY,mod)  ||
       $ModelProperty(LOG_DIRECTORY,mod)
     } {
    MenuExec::showCaseDirectoriesMessage
  }

  # Reset field states
  ModelProperty::resetFieldStates

  # Model file needs update
  set Model(Front,needsUpdate) 1

  # Store new values into Model
  Util::cpp_exec modelPropertyPanelOk

  Interface::applyModelData

  if { $Model(activeMeshIndices) == "" } {
    Solver::updateActiveMeshInfo
    Util::cpp_exec meshDefinePanelOk
  }


  Panel::panelDataChanged 0 ModelProperty 
  Panel::panelDataModified 0 ModelProperty 

  #-Datfile must be updated!
  set Datafile(SOLVER_INPUT_FILE) [Datafile::getDefaultSolverInputFile]
  set Datafile(MESH_INPUT_FILE) [Datafile::getDefaultMeshInputFile]
  set Datafile(OUTPUT_FILE) [Datafile::getDefaultOutputFile]
  set Datafile(POST_FILE) [Datafile::getDefaultPostFile]
  Datafile::panelSave

  # Update window title
  Util::updateMainWindowTitle

}


proc ModelProperty::panelOk {w} {
  global ModelProperty

  #---No changes
  if { !$ModelProperty(dataChanged) } {
    Panel::cancel $w; return
  }

  # Error in data, dont close the panel
  if { -1 == [ModelProperty::panelCheck] }  {
    return
  }

  ModelProperty::panelSave

  Panel::cancel $w
}


proc ModelProperty::panelApply {} {
  global ModelProperty

  #---No changes
  if { !$ModelProperty(dataChanged) } {
    return
  }

  # Error in data, dont close the panel
  if { -1 == [ModelProperty::panelCheck] }  {
    return
  }

  ModelProperty::panelSave
}


proc ModelProperty::panelCancel {w} {

  if { ![Panel::verifyCancel ModelProperty] } {
    return
  }

  #-Restore old values
  Panel::restoreFields ModelProperty

  Panel::cancel $w
}


# ModelProperty panel check
#
proc ModelProperty::panelCheck {} {
  global Info Model ModelProperty Datafile

  # Check model and problem name
  # ============================
  # NOTE: A name must be ok, if not
  # it is forced to empty string to disable
  # storing with improper (path)name

  if { -1 == [ModelProperty::checkModelNames] } {
    return -1
  }

  # Model directory
  # ===============
  ModelProperty::checkModelDirectoryEntry

  ModelProperty::applyModelDirectoryValue

  set needed 0
  if { -1 == [ModelProperty::checkModelDirectory $needed] } {
    return -1
  }

  # Check problem directory values
  # ==============================
  if { -1 == [ModelProperty::checkCaseDirectoryValues] } {
    return -1
  }

  #-Update problem directory value types
  foreach fld { INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY } {
    if { $ModelProperty($fld,mod) } {
      set ModelProperty($fld,valueType) Session
      set ModelProperty($fld,session) $ModelProperty($fld)
    }
  }

  # Trim descriptions
  # =================
  #set ModelProperty(DESCRIPTION) [$ModelProperty(description_widget) get 1.0 end ]

  foreach fld { MODEL_DESCRIPTION PROBLEM_DESCRIPTION } {

    set str [$ModelProperty($fld,widget) get 1.0 end ]
    set str [string trimright $str]

    $ModelProperty($fld,widget) delete 1.0 end
    $ModelProperty($fld,widget) insert 1.0 $str

    set ModelProperty($fld) $str
  }

}


proc ModelProperty::checkModelNames {} {
  global Info ModelProperty Model

  # Check model and problem name
  # ============================
  # NOTE: A name must be ok, if not
  # it is forced to empty string to disable
  # storing with improper (path)name

  foreach vn {MODEL_NAME PROBLEM_NAME RESULT_MESH_NAME} \
          ns {model problem mesh} {

    if { $ModelProperty($vn) == "" } {
      continue
    }

    set ModelProperty($vn) [string trim $ModelProperty($vn)]
    set data $ModelProperty($vn)

    # Check that name is proper for a directory
    #
    if { ![ModelProperty::checkNameForFile $data $ns] } {
      return -1
    }

    # Name ok
    set ModelProperty($vn) $data
  }

  #-No model name given
  if { $ModelProperty(MODEL_NAME) == "" } {

    set Info(messageIcon) warning
    set msg [list "No model name given, model cannot be saved!" ]

    if { "cancel" == [Message::dispOkCancelMessage $msg] } {
      return -1
    } else {
      return 1
    }
  }

  # Update other MODEL_NAME dependent names
  # =======================================
	
  ModelProperty::setModelCaseName

  set Model(EMF_FILE) "$Model(CASE_NAME).emf"

}


#-Check that name is simple file or directory name
# (so that it can be also a filename for the model file!)
#
proc ModelProperty::checkNameForFile {name name_type} {
  global Info

  set data [string trim $name " " ]


  if { "." != [file dirname $data] } {
    set msg [list "Invalid $name_type name!\n" \
                  "Name should not contain any path separators!" ]

    set Info(messageIcon) error
    Message::dispOkMessage $msg 
    return 0
  }

  if { [regexp -nocase "\[^A-Z^0-9^\n^\ ^\_^\-\]" $data] } {
    set msg [ list "Sorry, only letters, digits and characters  (-, ,_) \n" \
             "are legal in $name_type name!" ]
    set Info(messageIcon) error
    Message::dispOkMessage $msg 
    return 0
  }

  return 1
}


proc ModelProperty::checkModelDirectory { {needed 0} } {
  global Info Model ModelProperty 

  # Full model path (model directory + model name)
  # ===============
  # Create if needed (and if user wants!)
  # Set model case name

  set dir $ModelProperty(MODEL_DIRECTORY,absolute)

  if {$needed} {

    if { "ok" != [Util::checkAndCreateDirectory2B $dir "for the Model Directory" 1 "" 1 ""] } {
      return -1
    }

  } else {

    set dir_msg "for the model files"
    set cre_msg "NOTE: Cannot save the model file if the directory does not exist!\n\n"
    append cre_msg "Press  Yes to create the directory and continue\n"
    append cre_msg "Press  No  to continue without creating the directory"
    
    set rc [Util::checkAndCreateDirectory3B $dir $dir_msg 1 $cre_msg]

    if { $rc == "cancel" } {
      return -1
    }
    
    if { $rc == "no" } {
      return 1
    }
  }
}


# Check values set for problem directories
#
proc ModelProperty::checkCaseDirectoryValues {} {
  global Info ModelProperty Model

  # We use model-directory-absolute as the base
  set base $ModelProperty(MODEL_DIRECTORY,absolute)

  # Include path
  # ------------
  set fld INCLUDE_PATH
  set inc_sep $Info(includePathSeparator)

  set session $ModelProperty($fld,session,save)
  set model $ModelProperty($fld,model,save)

  set dirs $ModelProperty($fld)
  set dirs [string trim [string trim $dirs] $inc_sep]
  set Modelproperty($fld) $dirs

  set dirs [split $ModelProperty(INCLUDE_PATH) $inc_sep]
  set inc_path_abs ""

  foreach dir $dirs {

    set dir [string trim $dir]
    set dir_abs [Util::makePathAbsolute $dir $base] 

    if { $dir_abs != "" && ![file isdirectory $dir_abs] } {
      Message::showMessage "WARNING! Cannot find Include Path: $dir_abs !" $Info(wrnMsgColor)
    }

    append inc_path_abs $dir_abs 
    append inc_path_abs "; "
  }

  set ModelProperty($fld,absolute) [string trim [string trim $inc_path_abs] $inc_sep]


  # Results directory
  # -----------------
  set fld RESULTS_DIRECTORY
  set session $ModelProperty($fld,session,save)
  set model $ModelProperty($fld,model,save)
  set dir $ModelProperty($fld)

  PanelCheck::checkFields ModelProperty $fld

  set ModelProperty($fld,absolute) [Util::makePathAbsolute $dir $base]
  set dir_abs $ModelProperty($fld,absolute)

  if { $dir_abs != ""  } {
    
    set dir_msg "for the result files"
    set cre_msg "WARNING: Cannot run Solver if results directory does not exist!"
    
    if { "cancel" == [Util::checkAndCreateDirectory3B $dir_abs $dir_msg 1 $cre_msg] } {
      return -1
    }
  }


  # Log files directory
  # -------------------------
  set fld LOG_DIRECTORY
  set session $ModelProperty($fld,session,save)
  set model $ModelProperty($fld,model,save)
  set dir $ModelProperty($fld)
    
  PanelCheck::checkFields ModelProperty $fld

  set ModelProperty($fld,absolute) [Util::makePathAbsolute $dir $base]
  set dir_abs $ModelProperty($fld,absolute)
    
  if { $dir_abs != "" && ![file isdirectory $dir_abs] } {
    #Message::showMessage "WARNING! Cannot find Log Directory: $dir_abs !" $Info(wrnMsgColor)
  }

  # Data ok
  # =======
  return 1
}


# Apply selection done in a problem directory entry's
# value type option menu
#
proc ModelProperty::applyDirectoryValueType {fld} {
  global ModelProperty

  if { $ModelProperty($fld,valueType,prev) == 
       $ModelProperty($fld,valueType)
     } {

    return
  }

  Panel::panelDataChanged 1 ModelProperty 

  # Map menu keyword to variable code
  set key $ModelProperty($fld,valueType)
  set code [ModelProperty::valueTypeMenuKey2Var $key]

  # Set selected value type as active
  set ModelProperty($fld) $ModelProperty($fld,$code)
}


# A Cover for calling proper entry check proc based
# on argument field
#
proc ModelProperty::checkEntry {fld} {
  global ModelProperty

  switch $fld {

    MODEL_DIRECTORY {
      ModelProperty::checkModelDirectoryEntry
    }
  }
}


# Check model directory entry value
#
proc ModelProperty::checkModelDirectoryEntry {} {
  global Info ModelProperty

  PanelCheck::checkFields ModelProperty MODEL_DIRECTORY
}


# Apply model directory entry value
#
proc ModelProperty::applyModelDirectoryValue {} {
  global Info Model ModelProperty

  # Form absolute version of the model-directory
  # If given as relative, add working directory
  # 
  set dir [Util::formFullPath $ModelProperty(MODEL_DIRECTORY)]
 
  set base $Info(workingDirectory)

  set ModelProperty(MODEL_DIRECTORY,absolute) [Util::makePathAbsolute $dir $base]

  #---Update other MODEL_DIRECTORY dependent datafile names
  set Model(EMF_OPEN_DIR) $dir
  set Model(EMF_SAVE_DIR) $dir

  set Model(EMF_PATH) [file join $dir $Model(EMF_FILE)]

  # Update mesh names
  #set Model(meshNames) [ModelProperty::findModelMeshNames]

	# Load possible updated or new SOLVER.KEYWORD files
	#
  UserDefined::loadSolverKeywordsFiles
}


# Set default directory value
# NOTE: This is meant for ModelProperty - UserSetting dir pairs!!!
# Returns:
#  0: Not defined
# -1: Invalid value
#  1: Ok
#
proc ModelProperty::setDefaultDirectoryValue { dirvar {default_value ""} {abs_base ""} } {
  global Info ModelProperty UserSetting

  #--Check that variable is valid and that we have
  #  default value available
  if { ![info exists ModelProperty($dirvar)] ||
       ( $default_value == "" && 
        ![info exists UserSetting(DEFAULT_$dirvar)]
       )
     } {
    return 0 ; # NOT DEFINED
  }

  #--Set default user setting value if not yet defined
  if { [info exists UserSetting(DEFAULT_$dirvar)] &&
       $UserSetting(DEFAULT_$dirvar) == ""
     } {
    set UserSetting(DEFAULT_$dirvar) $default_value
  }

  #--Pick default value
  if { $default_value == "" } {
    set default_value $UserSetting(DEFAULT_$dirvar)
  }

  #--Set directory if corresponding ModelProperty value
  #  is not set earlier
  if { $ModelProperty($dirvar) == "" &&
       $default_value != ""
     } {
    set ModelProperty($dirvar) $UserSetting(DEFAULT_$dirvar)
  }

  #--Set absolute directory value

  set ModelProperty($dirvar,absolute) [Util::makePathAbsolute \
                                              $ModelProperty($dirvar) $abs_base]

  #--Check if directory is valid (using the absolute path!)
  if { $ModelProperty($dirvar) != "" &&
       ![file isdirectory $ModelProperty($dirvar,absolute)]
     } {
    return -1 ;# INVALID VALUE!
  }

  # Ok
  return 1
}


# Set value for mesh files directory
#
proc ModelProperty::setCurrentMeshDirectory {} {
  global Info ModelProperty Model

  set ModelProperty(CURRENT_MESH_DIRECTORY) \
        [ModelProperty::getMeshDirectory $Model(currentMeshName) 0]

  set ModelProperty(CURRENT_MESH_DIRECTORY,absolute) \
        [ModelProperty::getMeshDirectory $Model(currentMeshName) 1]
}


# Get value for the mesh files directory for "mesh_name"
#
proc ModelProperty::getMeshDirectory {mesh_name {absolute 1} } {
  global Info ModelProperty Model

  set dir [file join $Info(meshDirectoryName) $mesh_name]

  set base $ModelProperty(MODEL_DIRECTORY,absolute)

  if { $absolute } {
    set result [Util::makePathAbsolute $dir $base]

  } else {
    set result $dir
  }

  return $result
}
  

# Map visible menu keyword to an internal variable code
# for the problem directory entry option menu
#
proc ModelProperty::valueTypeMenuKey2Var {key} {
  global ModelProperty

  set var ""

  # Code menu keyword
  switch $key {

    Settings {
      set var settings
    }

    Model {
      set var model
    }

    Session {
      set var session
    }
  }

  return $var
}


# Map internal variable code to a visible menu keyword
# for the problem directory entry option menu
#
proc ModelProperty::valueTypeMenuVar2Key {var} {
  global ModelProperty

  set key ""

  # Find menu keyword
  switch $var {

    settings {
      set key Settings
    }

    model {
      set key Model
    }

    session {
      set key Session
    }
  }

  return $key
}


# Check if user really accepts the change to a new mesh name
# Ok: return 1
# No: return 0
# 
proc  ModelProperty::verifyMeshNameChange {old_name new_name} {
  global Info

  if { $new_name != $old_name &&
       $old_name != ""
     } {

    set msg [list "NOTE: Mesh name was changed and it will be saved in the directory: \n\n" \
                  "...MESHDIR/$new_name\n\n" \
                  "when model file is saved!\n\n" \
                  "Press Ok to continue!" ]


    if { ![Message::verifyAction $Info(advancedUser) $msg ok warning] } {
        return 0
    } else {
        return 1
    }
  }

  return 1
}


proc ModelProperty::findModelMeshNames {} {
  global Info Model ModelProperty

  if { $ModelProperty(MODEL_DIRECTORY,absolute) == "" } {
    return ""
  } 

  set mdir [file join $ModelProperty(MODEL_DIRECTORY,absolute) $Info(meshDirectoryName)] 

  # These meshes are in the mesh directory
  set dir_names [Util::getSubdirectoryList $mdir]


  if { [info exists Model(meshNames)] } {
    set old_names $Model(meshNames)
  } else {
    set old_names ""
  }

  set new_names ""

  # Pick first all old model mesh names which match the dir meshes
  foreach mn $old_names {
    
    set pos [lsearch $dir_names $mn]
    
    # Append each matched mesh to the new list and remove
    # it from the dir-meshes list
    if { $pos != -1 } {
      lappend new_names $mn
      set dir_names [lreplace $dir_names $pos $pos]
    }
  }

  # Append what is left in the dir-meshes list to the list of
  # new model meshes
  foreach mn $dir_names {
    lappend new_names $mn
  }

        
  return $new_names
}


# Not in use? MVe 25.09.00
#
proc ModelProperty::updateMeshNamesAndIndices {} {
  global Info Model

  set old_mesh_names $Model(meshNames)
  set new_mesh_names [ ModelProperty::findModelMeshNames]

  set old_indices ""

  foreach old_nm $old_mesh_names {

    set idx  [lsearch $new_mesh_names $old_nm]

    if { -1 != $idx } {
      lappend old_indices $idx
    }
  }

  ModelProperty::updateGridVariable "gr" $Info(NO_INDEX) $old_indices
  ModelProperty::updateGridVariable "H" -1.0 $old_indices
}


# Update specified mesh parameter related variable.
# NOTE: Nof values must be the same as nof-meshes!
# Arg: var_vame (like "gr" or "H")
#      default_value (like $Info(NO_INDEX) or -1.0)
#      old_indices (indices for old value in the new list)
#
proc ModelProperty::updateGridVariable { var_name default_value old_indices } {
  global Info Model ObjectTable

  set new_values ""

  # Init new values (for each mesh-name!)
  foreach nm $Model(meshNames) {
    lappend new_values $default_value
  }

  # Loop all objects
  # ================
  foreach id $ObjectTable(ids) {

    # Skip if no grid variable defined
    if { ![info exists ObjectTable($id,$var_name)] } {
      continue
    }

    # Insert old values into proper (new) places
    # in the new values list
    #
    set old_values $ObjectTable($id,$var_name)
    
    foreach old_value $old_values idx $old_indices {

      # If old mesh is not any more in use
      # ----------------------------------
      if { $idx == $Info(NO_INDEX) } {
        continue
      }
    
      # Insert valid value
      # ------------------
      set new_values [lreplace $new_value $idx $idx $old_value]
    } 

    # Store new values
    # ----------------
    set ObjectTable($id,$var_name) $new_values

  }

}


#End ecif_modelPropertiesPanel.tcl

