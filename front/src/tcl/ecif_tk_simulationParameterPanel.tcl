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
#Module:    ecif_tk_simulationParameterPanel.tcl
#Language:  Tcl
#Date:      13.02.01
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the user defined simulation section parameters
#
#************************************************************************


# This procedure displays the user defined simulation parameters
#
#------SimulationParameter definitions  proc------
#
proc SimulationParameter::openPanel {} {
  global Info SimulationParameter Model

  set w $SimulationParameter(winName)
  set wgeom $SimulationParameter(winGeometry)

  set Info(thisWindow) $w
  set this $w

  #--Store windows-id
  set id [winfo atom $w]
  set SimulationParameter(winId) $id

  if { 1 == [Util::checkPanelWindow SimulationParameter $id $SimulationParameter(winTitle) $wgeom] } {
    return
  }  

  set SimulationParameter(dataChanged) 0
  set SimulationParameter(dataModified) 0

  toplevel $w
  focus $w

  wm title $w $SimulationParameter(winTitle)
  wm geometry $w $wgeom 

  Panel::resetFields SimulationParameter

  Panel::initFields SimulationParameter

  set id $SimulationParameter(parameterId)
  if { [info exists SimulationParameter($id,data)] } {
    DataField::formDataFields SimulationParameter $SimulationParameter($id,data)
  }

  Panel::backupFields SimulationParameter

  #-----WIDGET CREATION
  frame $w.f1 ;#--Fields
  frame $w.fB ;#--Buttons

  StdPanelCreate::setNofValuesAreaFrames SimulationParameter
  StdPanelCreate::createValuesArea $w.f1 SimulationParameter
  PanelCheck::execPanelFillProcs SimulationParameter
  StdPanelExec::setValuesAreaActivity SimulationParameter ""
  StdPanelCreate::packValuesArea $w.f1 SimulationParameter

  set SimulationParameter(dataChanged) 0
  set SimulationParameter(dataModified) 0

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

  set ok_btn [button $w.fB.ok -text OK -command "SimulationParameter::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "SimulationParameter::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command SimulationParameter::panelApply \
                                 -state $ap]

  focus $ok_btn
  set SimulationParameter(applyButton)  $ap_btn
  set SimulationParameter(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx 
  
  #-----Initialization
  #-Nothing so far

  # Set field label bindings for right-button help
  Widget::setLabelBindings SimulationParameter

}


proc SimulationParameter::panelSave { {inform_front 1} } {
  global Info SimulationParameter Model

  #--Store old values
  Panel::backupFields SimulationParameter

  #--Form parameter data
  set SimulationParameter(ids) 1
  DataField::formNonStandardParameter SimulationParameter 1 "Simulation1"

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 SimulationParameter 
  Panel::panelDataModified 0 SimulationParameter 
  StdPanelExec::setValuesAreaStatus SimulationParameter 0

  Util::cpp_exec simulationParameterPanelOk
}


proc SimulationParameter::panelOk {w} {
  global SimulationParameter

  #---No changes
  if { !$SimulationParameter(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![SimulationParameter::checkPanelData] } {
    return
  }

  #---Save data
  SimulationParameter::panelSave
  Panel::cancel $w
} 


proc SimulationParameter::panelApply {} {
  global SimulationParameter

  #---No changes
  if { !$SimulationParameter(dataChanged) } {
    return
  }

  #---Error in data
  if { ![SimulationParameter::checkPanelData] } {
    return
  }

  SimulationParameter::panelSave
}


proc SimulationParameter::panelCancel {w} {
  global Info SimulationParameter

  if { ![Panel::verifyCancel SimulationParameter] } {
    return
  }

  #---Reset into old values
  Panel::restoreFields SimulationParameter

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc SimulationParameter::checkPanelData {} {
  global SimulationParameter

  #-Check fields
  if { ![SimulationParameter::checkData] } {
    StdPanelExec::setValuesAreaStatus SimulationParameter
    return 0
  }

  return 1
}


# Return 1 = ok, 0 = error
#
proc SimulationParameter::checkData {} {
  global Info SimulationParameter Model

  set totmsg ""

  #-Check numeric variables
  foreach fld $SimulationParameter(allFields) {

    # Check only current target fields
    if { !$SimulationParameter($fld,act) } {
      continue
    }
    
    if { ![info exist SimulationParameter($fld)] } {
      continue
    }

    set msg [PanelCheck::checkFieldValue SimulationParameter $fld]

    #-Error messages are collected
    if {$msg != ""} {

      set SimulationParameter($fld,err) 1

      set label [DataField::getFieldProperty SimulationParameter $fld Label]
      #lappend totmsg [string toupper $label]
      lappend totmsg [string toupper $fld]
      lappend totmsg "\n\n$msg"

    #-Field Ok
    } else {
      set SimulationParameter($fld,err) 0
    }
  }

  #-Ok (no msg-variable was created!)
  if { $totmsg == "" } {
    return 1

  #-Error message is displayed
  } else {
    set totmsg [linsert $totmsg 0 $Info(fieldMsg)]

    set Info(messageIcon) error
    Message::dispOkMessage $totmsg  "$Info(FRONT_NAME) message!" $SimulationParameter(winName)
    return 0
  }  

  return 1
}


# end ecif_tk_simulationParameterPanel.tcl
# ********************
