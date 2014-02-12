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
#Module:    ecif_tk_userSettings.tcl
#Language:  Tcl
#Date:      07.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for user settings for Elmer Front
#
#************************************************************************
 
#------User settings panel proc------
#
# This procedure displays the model statistic-info screen
#
proc UserSetting::openPanel { } {
  # Global variables
  global Info UserSetting ModelFlags

  set w $UserSetting(winName)
  set wgeom $UserSetting(winGeometry)

  set id [winfo atom $w]
  set UserSetting(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow UserSetting $id $UserSetting(winTitle) $wgeom] } {
    return
  }  

  set UserSetting(dataChanged) 0
  set UserSetting(dataModified) 0
 
  toplevel $w
  focus $w 
  set this $w

  wm title $w $UserSetting(winTitle)
  wm geometry $w $wgeom 

  Panel::initFields UserSetting
  Panel::backupFields UserSetting

  StdPanelCreate::setNofValuesAreaFrames UserSetting

  # Initial values for save flags
  set UserSetting(saveForSession) 1
  set UserSetting(saveInModel) 0
  set UserSetting(saveInFile) 0
   

  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #-We set each row of data into a separate frame.
  frame $w.f1      ;# Values area
  frame $w.f2      ;# Buttons area outer
  frame $w.fSave   ;# Save buttons area 
  frame $w.fOk     ;# Ok,cancel buttons area inner


  #-Values 
  StdPanelCreate::createValuesArea $w.f1 UserSetting
  StdPanelCreate::packValuesArea $w.f1 UserSetting

  UserSetting::applyUserLevel

  pack $w.f1 $w.f2 -side top -expand 1 -fill both -padx $fpx3 -pady $fpy2

  # For buttons which needs model
  set mdl_state normal
  if { !($ModelFlags(GEOMETRY_TYPE_CAD) || $ModelFlags(GEOMETRY_TYPE_MESH)) } {
    set mdl_state disabled
  }

  pack $w.fSave $w.fOk -in $w.f2 -side top -pady $fpy2

  #-Save, Cancel buttons
  checkbutton $w.save_for_session -text "Save for session" \
                                  -variable UserSetting(saveForSession)
  checkbutton $w.save_in_model   -text "Save in model" \
                                 -variable UserSetting(saveInModel) -state $mdl_state
  checkbutton $w.save_in_file    -text "Save in file" \
                                 -variable gloUserSetting(saveInFile)

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.ok     -text OK     -command "UserSetting::panelOk $this"]
  set cn_btn [button $w.cancel -text Cancel -command "UserSetting::panelCancel $this" -state $ca]
  set ap_btn [button $w.apply  -text Apply  -command "UserSetting::panelApply" -state $ap]

  focus $ok_btn
  set UserSetting(applyButton)  $ap_btn
  set UserSetting(cancelButton) $cn_btn

  pack $w.save_for_session $w.save_in_model $w.save_in_file -in $w.fSave -side left -padx $fpx3

  pack $ok_btn $cn_btn $ap_btn -in $w.fOk -side left -padx $fpx3
}


proc UserSetting::panelSave {} {
  global UserSetting

  if { $UserSetting(saveForSession) } {
    UserSetting::saveForSession
  }

  if { $UserSetting(saveInModel) } {
    UserSetting::saveInModel
  }

  if { $UserSetting(saveInFile) } {
    UserSetting::saveInSettingsFile
  }

  Panel::panelDataChanged 0 UserSetting
  Panel::panelDataModified 0 UserSetting
  StdPanelExec::setValuesAreaStatus UserSetting 0

}


# Save data into model
#
proc UserSetting::saveForSession {} {
  global Info Model UserSetting

  set panel $UserSetting(parameterType)

  #--Update problem directories
  #
  # If some of the problem directories has been modified, show also message
  set dir_vars { MODEL_DIRECTORY INCLUDE_PATH RESULTS_DIRECTORY LOG_DIRECTORY }

  set modified 0
  foreach vn $dir_vars {
    if {$UserSetting(DEFAULT_$vn,mod)} {
      set modified 1
    }
  }

  set force 0
  set show_msg $modified
  MenuExec::applySettingsToCaseDirectories $force $show_msg

  ModelProperty::applyModelDirectoryValue

  #--Other directories (CAD-files, MESH-files)
  if { $UserSetting(DEFAULT_CAD_FILES_DIRECTORY,mod) } {
    set Model(CAD_OPEN_DIR) $UserSetting(DEFAULT_CAD_FILES_DIRECTORY)
  }

  if { $UserSetting(DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY,mod) } {
    set Model(MESH_OPEN_DIR) $UserSetting(DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY)
  }

  #--Store old values
  Panel::backupFields UserSetting

  StdPanelExec::setValuesAreaStatus UserSetting 0
}


# Save data also into model
#
proc UserSetting::saveInModel {} {
  global Info UserSetting Model

  #--Save first values into session
  #UserSetting::saveForSession

  #--Form parameter data
  set UserSetting(ids) 1
  DataField::formNonStandardParameter UserSetting 1 "UserSetting1"

  set Model(Front,needsUpdate) 1

  #--Write data into model
  Util::cpp_exec userSettingPanelOk
}


# Save setting in the (defaults) settings file
#
proc UserSetting::saveInSettingsFile {} {
  global Info UserSetting

  #--Save first values into session
  #UserSetting::saveForSession

  set path [MenuExec::saveFile "User settings default file" "" $UserSetting(filePath) ]

  if {$path == ""} {
    return
  }

  set UserSetting(filePath) $path

  Util::cpp_exec saveUserSettingFile $UserSetting(filePath)
}


# OK proc 
proc UserSetting::panelOk {w} {
  global UserSetting

  #---No changes
  if { !$UserSetting(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![UserSetting::checkPanelData] } {
    return
  }

  UserSetting::panelSave

  #-Reset oldvalues, if not saved for session
  if { !$UserSetting(saveForSession) } {
    UserSetting::panelCancel $w 1

  #-Otherwise just close the window
  } else {
    Panel::cancel $w
  }

}


proc UserSetting::panelApply {} {
  global UserSetting

  #---No changes
  if { !$UserSetting(dataChanged) } {
    Panel::cancel $w; return
  }

  set code [UserSetting::checkPanelData]

  #---Error
  if { ![UserSetting::checkPanelData] } {
    return
  }

  UserSetting::panelSave
}


proc UserSetting::panelCancel { w {force_cancel 0} } {
  global Info UserSetting 

  if { ![Panel::verifyCancel UserSetting] } {
    return
  }

  #if { !$force_cancel && [Panel::panelVarsChanged UserSetting] } {
  #  set msg [list "NOTE: Panel data was changed, but not saved!\n" \
  #                $Info(anywayOk) ]
  #
  #  set Info(messageIcon) warning
  #
  #  if { "cancel" == [ Message::dispOkCancelMessage $msg ] } {
  #    return 
  #  }
  #}

  #---Reset into old values
  Panel::restoreFields UserSetting

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc UserSetting::checkPanelData {} {
  global UserSetting ModelProperty Model

  return 1
}


# NOTE: This is used to check automatically loaded
# settings files data
# (Panel data it checked when it is entered!)
#
proc UserSetting::checkData {} {
  global Info UserSetting

  set dvars {
    DEFAULT_MODEL_DIRECTORY
    DEFAULT_CAD_FILES_DIRECTORY
    DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY
    DEFAULT_RESULTS_DIRECTORY
    DEFAULT_INCLUDE_PATH
    DEFAULT_LOG_DIRECTORY
  }

  foreach dv $dvars {
    set UserSetting($dv) [Util::makeTclFilename $UserSetting($dv)]
  }
}


# Apply (changed) UserLevel effects for the panel
#
proc UserSetting::applyUserLevel {} {
  global Info UserSetting Model

  if { ![winfo exists $UserSetting(winName)] } {
    return
  }

  if { $Info(userLevel) < $Info(powerUser) } {
    set state disabled
  } else {
    set state normal
  }

  # These buttons are UserLevel dependent
  foreach cb { AUTO_SAVE_MODEL AUTO_SAVE_SOLVER_INPUT} {
    $UserSetting(allWidgets,$cb) configure -state $state
  }
}


# end ecif_tk_userSettingsPanel.tcl
# ********************
