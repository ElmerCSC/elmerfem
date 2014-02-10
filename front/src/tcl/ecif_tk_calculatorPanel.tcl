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
#Module:    ecif_tk_calculatorPanel.tcl
#Language:  Tcl
#Date:      24.01.01
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for calculator solver settings
#
#************************************************************************


#------Calculator definitions  proc------
# This is a "standard" panel
#
proc Calculator::openPanel { } {
  global Info Calculator Model

  set w $Calculator(winName)
  set wgeom $Calculator(winGeometry)

  set Info(thisWindow) $w
  set this $w

  #--Store windows-id
  set id [winfo atom $w]
  set Calculator(winId) $id

  if { 1 == [Util::checkPanelWindow Calculator $id $Calculator(winTitle) $wgeom] } {
    return
  }  

  toplevel $w
  focus $w

  wm title $w $Calculator(winTitle)
  wm geometry $w $wgeom 

  Panel::initFields Calculator

  StdPanelCreate::setNofValuesAreaFrames Calculator

  #-----WIDGET CREATION
  frame $w.f1 ;#--Parameter box
  frame $w.f1.fParamBox
  frame $w.f1.fParamBox.fParams

  frame $w.f2 ;#--Buttons1
  frame $w.f2.fButtons1 ;#(Add etc.)

  frame $w.f3 ;#--Values area
  frame $w.f3.fValues -relief groove -bd 2

  frame $w.fB ;#--Buttons2 ;#(ok etc.)
  frame $w.fB.fButtons2 

  set fpx0 0
  set fpy0 0
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)
  set expnd 1
  set fill both
  set anch w

  #-Parameters
  StdPanelCreate::createParamBoxArea $w.f1.fParamBox Calculator "" 50
  Calculator::createButtons1Area $w.f2.fButtons1 Calculator
  StdPanelCreate::createValuesArea $w.f3.fValues Calculator
  StdPanelCreate::createButtons2Area $w.fB.fButtons2 Calculator $w

  StdPanelExec::setValuesAreaActivity Calculator ""


  #---WIDGET PACKING
  set fpx $Info(framePadX1)
  set fpy $Info(framePadY1)

  #-----Parameter box
  pack $w.f1 -side top  -anchor nw -fill x -padx $fpx -pady $fpy

  pack $w.f1.fParamBox -side left -padx 2m -expand $expnd -fill y -anchor c \
                       -padx $fpx0 -pady $fpy0

  #-----Buttons-1 frame (Attach etc, Body etc. names)
  pack $w.f2 -side top  -anchor nw -fill x -padx $fpx -pady $fpy

  pack $w.f2.fButtons1 -expand $expnd -fill $fill \
                       -padx $fpx0 -pady $fpy0

  #-----Values area
  pack $w.f3 -side top  -anchor nw -fill x -padx $fpx -pady $fpy

  pack $w.f3.fValues -expand $expnd -fill $fill -anchor $anch \
                     -padx $fpx1 -pady $fpy1

  #-----Buttons packing widgets packing
  pack $w.fB -side top  -padx $fpx -pady $fpy

  pack $w.fB.fButtons2 -expand $expnd  -padx $fpx0 -pady $fpy0

  StdPanelCreate::packParamBoxArea $w.f1.fParamBox Calculator
  Calculator::packButtons1Area $w.f2.fButtons1 Calculator
  StdPanelCreate::packValuesArea $w.f3.fValues Calculator
  StdPanelCreate::packButtons2Area $w.fB.fButtons2 Calculator

  #-----Initialization
  StdPanelInit::createWidgetBindings Calculator
  StdPanelInit::initPanelData Calculator

  # Set field label bindings for right-button help
  Widget::setLabelBindings Calculator
}


# Special version for Calculators (only Add etc. buttons)
proc Calculator::createButtons1Area {frame globArray} {
  global Common Info ModelFlags
  upvar #0 $globArray theArray

  frame $frame.sf1
  frame $frame.sf2

  # Add , Update and Delete buttons
  button $frame.sf1.add -text Add -command "StdPanelExec::addParameter $globArray"
  button $frame.sf1.update -text Update -command "StdPanelExec::updateParameter $globArray"
  button $frame.sf1.delete -text Delete -command "StdPanelExec::deleteParameter $globArray"

  set wdg $frame.sf1.add;     bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  set theArray(panelAddButton) $wdg

  set wdg $frame.sf1.update;  bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  $wdg configure -state disabled
  set theArray(panelUpdateButton) $wdg

  set wdg $frame.sf1.delete;  bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"

  if { ![info exists theArray(ids)] ||
       0 == [llength $theArray(ids)]
     } {
    $wdg configure -state disabled
  }

  set theArray(panelDeleteButton) $wdg

  # These updates make data "unmodified"
  set wdg $frame.sf1.add;     bind $wdg <ButtonRelease-1> "+Panel::panelDataModified 0 $globArray $wdg {%A %K}"
  set wdg $frame.sf1.update;  bind $wdg <ButtonRelease-1> "+Panel::panelDataModified 0 $globArray $wdg {%A %K}"

  #-Store buttons to panel specific array
  set m $frame.sf1
  set theArray(buttons1) [list $m.add $m.update $m.delete]

  # Parameter name is editable
  # ==========================
  label $frame.sf2.pname_label -text "Name: "

  if { $theArray(hasBoundaries) } {
    set wid 18
  } else {
    set wid 25
  }

  # Parameter name variable
  set pn_var $globArray
  append pn_var (parameterName)

  entry $frame.sf2.pname -textvariable $pn_var -font $Info(entryFont)  -width $wid

  set wdg $frame.sf2.pname
  set theArray(parameterNameWidget) $wdg
  set theArray(parameterName,err) 0
  set theArray(parameterName,mod) 0

  # Set bindings
  set pn "parameterName"
  bind $wdg <KeyRelease> "+StdPanelExec::updateCurrentParameterName $globArray $wdg"
  bind $wdg <KeyRelease> "+Panel::panelDataModified 1 $globArray $wdg"
  bind $wdg <KeyRelease> "+Widget::entryKeyRelease $globArray $pn $wdg {%A %K}"
  bind $wdg <FocusIn> "+Panel::setProcAndTableButtonStates $globArray disabled disabled"
  bind $wdg <KeyPress-Escape> "+Widget::entryKeyPress-Escape $globArray $pn $wdg {%A %K}"

}


# Calculator sepcific
proc Calculator::packButtons1Area {frame globArray} {
  global Common Info
  upvar #0 $globArray theArray

  set px 5; set py 4

  # Add and Delete buttons
  pack $frame.sf1.delete -side right -pady $py -padx $px -anchor e
  pack $frame.sf1.update -side right -pady $py -padx $px -anchor e
  pack $frame.sf1.add  -side right -pady $py -padx $px -anchor e

  # Calculator name
  pack $frame.sf2.pname -side right -pady 1 -padx 2 -anchor w
  pack $frame.sf2.pname_label -side right -pady 1 -padx 2 -anchor w 
  
  pack $frame.sf1 $frame.sf2 -side left -fill y -padx 3
}


# Update variable dofs (after coordinate changes etc)
#
proc Calculator::updateVariableDofs { {is_modified ""} } {
  global Calculator

  if { $is_modified != "" } {
    upvar $is_modified modified
  }

  set modified 0

  foreach id $Calculator(ids) {

    set all_dofs [DataField::getFieldValue Calculator $id VARIABLE_DOFS_ALL]

    set old_dofs [DataField::getFieldValue Calculator $id VARIABLE_DOFS]
    set new_dofs [Calculator::getVariableDofs $all_dofs]

    if { $old_dofs != $new_dofs } {
      DataField::setFieldValue Calculator $id VARIABLE_DOFS $new_dofs
      set modified 1
    }
  }
}


proc Calculator::getVariableDofs {all_dofs} {
  global Calculator Model

  set values $all_dofs
  
  if { 2 == [llength $values] } {
    set values [linsert $values 0 1]
  }

  #-No dimension dependence
  #
  if { 1 == [llength $values] } {
    set dofs $values

  #-Depends on simulation dimension
  #
  } elseif { $Model(SIMULATION_DIMENSION) == "1D" } {
    set dofs [lindex $values 0]

  } elseif { $Model(SIMULATION_DIMENSION) == "2D" } {
    set dofs [lindex $values 1]

  } else {
    set dofs [lindex $values end]
  }

  return $dofs
}


# end ecif_tk_calculatorPanel.tcl
# ********************
