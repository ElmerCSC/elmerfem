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
#Module:    ecif_tk_modelParameterPanel.tcl
#Language:  Tcl
#Date:      13.02.01
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the user defined model (generic section) parameters
#
#************************************************************************


# This procedure displays the user defined model parameters
#
#------ModelParameter definitions  proc------
#
proc ModelParameter::openPanel {} {
  global Info ModelParameter Model

  set w $ModelParameter(winName)
  set wgeom $ModelParameter(winGeometry)

  set Info(thisWindow) $w
  set this $w

  #--Store windows-id
  set id [winfo atom $w]
  set ModelParameter(winId) $id

  if { 1 == [Util::checkPanelWindow ModelParameter $id $ModelParameter(winTitle) $wgeom] } {
    return
  }  

  set ModelParameter(dataChanged) 0
  set ModelParameter(dataModified) 0

  toplevel $w
  focus $w

  wm title $w $ModelParameter(winTitle)
  wm geometry $w $wgeom 


  Panel::resetFields ModelParameter

  Panel::initFields ModelParameter

  set id $ModelParameter(parameterId)
  if { [info exists ModelParameter($id,data)] } {
    DataField::formDataFields ModelParameter $ModelParameter($id,data)
  }

  Panel::backupFields ModelParameter

  #-----WIDGET CREATION
  frame $w.f1 ;#--Fields
  frame $w.fB ;#--Buttons


  StdPanelCreate::setNofValuesAreaFrames ModelParameter
  StdPanelCreate::createValuesArea $w.f1 ModelParameter
  PanelCheck::execPanelFillProcs ModelParameter
  StdPanelExec::setValuesAreaActivity ModelParameter ""
  StdPanelCreate::packValuesArea $w.f1 ModelParameter

  set ModelParameter(dataChanged) 0
  set ModelParameter(dataModified) 0

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

  set ok_btn [button $w.fB.ok -text OK -command "ModelParameter::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "ModelParameter::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command ModelParameter::panelApply \
                                 -state $ap]

  focus $ok_btn
  set ModelParameter(applyButton)  $ap_btn
  set ModelParameter(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx 
  
  #-----Initialization
  #-Nothing so far

  # Set field label bindings for right-button help
  Widget::setLabelBindings ModelParameter
}


proc ModelParameter::panelSave { {inform_front 1} } {
  global Info ModelParameter Model

  #--Store old values
  Panel::backupFields ModelParameter

  #--Form parameter data
  set ModelParameter(ids) 1
  DataField::formNonStandardParameter ModelParameter 1 "Simulation1"

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 ModelParameter 
  Panel::panelDataModified 0 ModelParameter 

  Util::cpp_exec modelParameterPanelOk
}


proc ModelParameter::panelOk {w} {
  global ModelParameter

  #---No changes
  if { !$ModelParameter(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![ModelParameter::checkPanelData] } {
    return
  }

  #---Save data
  ModelParameter::panelSave
  Panel::cancel $w
} 


proc ModelParameter::panelApply {} {
  global ModelParameter

  #---No changes
  if { !$ModelParameter(dataChanged) } {
    return
  }

  #---Error in data
  if { ![ModelParameter::checkPanelData] } {
    return
  }

  ModelParameter::panelSave
}


proc ModelParameter::panelCancel {w} {
  global Info ModelParameter

  if { ![Panel::verifyCancel ModelParameter] } {
    return
  }

  #---Reset into old values
  Panel::restoreFields ModelParameter

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc ModelParameter::checkPanelData {} {
  global Info ModelParameter Model

   # Ok
  return 1
}


# end ecif_tk_modelParameterPanel.tcl
# ********************
