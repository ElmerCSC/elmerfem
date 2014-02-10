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
#Module:    ecif_tk_solverControlPanel.tcl
#Language:  Tcl
#Date:      13.02.01
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the solver (ouput) control parameters
#
#************************************************************************


# This procedure displays the solver control panel
#
#------SolverControl panel proc------
#
proc SolverControl::openPanel {} {
  global Info SolverControl Model

  set w $SolverControl(winName)
  set wgeom $SolverControl(winGeometry)

  set Info(thisWindow) $w
  set this $w

  #--Store windows-id
  set id [winfo atom $w]
  set SolverControl(winId) $id

  if { 1 == [Util::checkPanelWindow SolverControl $id $SolverControl(winTitle) $wgeom] } {
    return
  }  

  set SolverControl(dataChanged) 0
  set SolverControl(dataModified) 0

  Panel::resetFields SolverControl
  Panel::initFields SolverControl
  set id $SolverControl(parameterId)
  if { [info exists SolverControl($id,data)] } {
    DataField::formDataFields SolverControl $SolverControl($id,data)
  }
  Panel::backupFields SolverControl

  toplevel $w
  focus $w

  wm title $w $SolverControl(winTitle)
  wm geometry $w $wgeom 

  #-----WIDGET CREATION
  frame $w.f1 ;#--Fields
  frame $w.fB ;#--Buttons

  StdPanelCreate::setNofValuesAreaFrames SolverControl
  StdPanelCreate::createValuesArea $w.f1 SolverControl
  PanelCheck::execPanelFillProcs SolverControl
  StdPanelExec::setValuesAreaActivity SolverControl ""
  StdPanelCreate::packValuesArea $w.f1 SolverControl

  set SolverControl(dataChanged) 0
  set SolverControl(dataModified) 0

  #---WIDGET PACKING
  set fpx $Info(framePadX1)
  set fpy $Info(framePadY1)

  #-----Fields
  pack $w.f1 -side top  -anchor nw -fill x -padx $fpx -pady $fpy

  #-----Buttons packing widgets packing
  pack $w.fB -side top  -padx $fpx -pady $fpy

  #-----Apply, Ok and cancel buttons creating and packing

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "SolverControl::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "SolverControl::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command SolverControl::panelApply \
                                 -state $ap]

  focus $ok_btn
  set SolverControl(applyButton)  $ap_btn
  set SolverControl(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx 
  
  #-----Initialization
  #-Nothing so far

  # Set field label bindings for right-button help
  Widget::setLabelBindings SolverControl
}


proc SolverControl::panelSave { {inform_front 1} } {
  global Info SolverControl Model

  #--Store old values
  Panel::backupFields SolverControl

  #--Form parameter data
  set SolverControl(ids) 1
  DataField::formNonStandardParameter SolverControl 1 "SolverControl1"

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 SolverControl 
  Panel::panelDataModified 0 SolverControl 
  StdPanelExec::setValuesAreaStatus SolverControl 0

  Util::cpp_exec solverControlPanelOk
}


proc SolverControl::panelOk {w} {
  global SolverControl

  #---No changes
  if { !$SolverControl(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![SolverControl::checkPanelData] } {
    return
  }

  #---Save data
  SolverControl::panelSave
  Panel::cancel $w
} 


proc SolverControl::panelApply {} {
  global SolverControl

  #---No changes
  if { !$SolverControl(dataChanged) } {
    return
  }

  #---Error in data
  if { ![SolverControl::checkPanelData] } {
    return
  }

  SolverControl::panelSave
}


proc SolverControl::panelCancel {w} {
  global Info SolverControl

  if { ![Panel::verifyCancel SolverControl] } {
    return
  }

  #---Reset into old values
  Panel::restoreFields SolverControl

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc SolverControl::checkPanelData {} {
  global SolverControl

  #-Check fields
  if { ![SolverControl::checkData] } {
    StdPanelExec::setValuesAreaStatus SolverControl
    return 0
  }

  return 1
}


# Return 1 = ok, 0 = error
#
proc SolverControl::checkData {} {
  global Info SolverControl Model

  set totmsg ""

  #-Check numeric variables
  foreach fld $SolverControl(allFields) {

    # Check only current target fields
    if { !$SolverControl($fld,act) } {
      continue
    }
    
    if { ![info exist SolverControl($fld)] } {
      continue
    }

    set msg [PanelCheck::checkFieldValue SolverControl $fld]

    #-Error messages are collected
    if {$msg != ""} {

      set SolverControl($fld,err) 1

      set label [DataField::getFieldProperty SolverControl $fld Label]
      #lappend totmsg [string toupper $label]
      lappend totmsg [string toupper $fld]
      lappend totmsg "\n\n$msg"

    #-Field Ok
    } else {
      set SolverControl($fld,err) 0
    }
  }

  #-Ok (no msg-variable was created!)
  if { $totmsg == "" } {
    return 1

  #-Error message is displayed
  } else {
    set totmsg [linsert $totmsg 0 $Info(fieldMsg)]

    set Info(messageIcon) error
    Message::dispOkMessage $totmsg  "$Info(FRONT_NAME) message!" $SolverControl(winName)
    return 0
  }  

  return 1
}


# end ecif_tk_solverControlPanel.tcl
# ********************
