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
#Module:    ecif_tk_constantPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the physical constants
#
#************************************************************************


# This procedure displays the physical constant definition panel
#
#------Constant definitions  proc------
#
proc Constant::openPanel { } {
  global Info Constant Model

  set w $Constant(winName)
  set wgeom $Constant(winGeometry)

  set Info(thisWindow) $w
  set this $w

  #--Store windows-id
  set id [winfo atom $w]
  set Constant(winId) $id

  if { 1 == [Util::checkPanelWindow Constant $id $Constant(winTitle) $wgeom] } {
    return
  }  

  set Constant(dataChanged) 0
  set Constant(dataModified) 0

  toplevel $w
  focus $w

  wm title $w $Constant(winTitle)
  wm geometry $w $wgeom 

  # 
  set id $Constant(parameterId)

  Panel::resetFields Constant

  Panel::initFields Constant
  
  if { [info exists Constant($id,data)] } {
    DataField::formDataFields Constant $Constant($id,data)
  }

  Panel::backupFields Constant

  StdPanelCreate::setNofValuesAreaFrames Constant

  # Form helper variables for screen handling
  Constant::formGravityComponents
    
  #-----WIDGET CREATION
  frame $w.f1 ;#--Gravity and other physical constants
  frame $w.fB ;#--Buttons

  StdPanelCreate::createValuesArea $w.f1 Constant
  StdPanelExec::setValuesAreaActivity Constant ""
  StdPanelCreate::packValuesArea $w.f1 Constant

  if { $Model(SIMULATION_DIMENSION) == "2D" } {
    Widget::configureEntry $Constant(allWidgets,GRAVITY_3) disabled
  }

  #---WIDGET PACKING
  set fpx $Info(framePadX1)
  set fpy $Info(framePadY1)

  #-----Gravity vector and physical constants frame
  pack $w.f1 -side top  -anchor nw -fill x -padx $fpx -pady $fpy

  #-----Buttons packing widgets packing
  pack $w.fB -side top  -padx $fpx -pady $fpy

  #-----Apply, Ok and cancel buttons creating and packing

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "Constant::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "Constant::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command Constant::panelApply \
                                 -state $ap]

  focus $ok_btn
  set Constant(applyButton)  $ap_btn
  set Constant(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx 
  
  #-----Initialization
  #-Nothing so far

  # Set field label bindings for right-button help
  Widget::setLabelBindings Constant
}


proc Constant::panelSave { {inform_front 1} } {
  global Info Constant Model

  Constant::formGravityComponents

  # Form final gravity variable
  set gravity ""
  append gravity $Constant(GRAVITY_1)
  append gravity " "
  append gravity $Constant(GRAVITY_2)
  append gravity " "
  append gravity $Constant(GRAVITY_3)
  append gravity " "
  append gravity $Constant(GRAVITY_ABS)

  set Constant(GRAVITY) $gravity

  set panel $Constant(parameterType)

  #--Store old values
  Panel::backupFields Constant

  set Constant(GRAVITY,dataSize) 4

  #--Form parameter data
  set Constant(ids) 1
  DataField::formNonStandardParameter Constant 1 "Constant1"

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 Constant 
  Panel::panelDataModified 0 Constant 

  Util::cpp_exec constantPanelOk
}


proc Constant::panelOk {w} {
  global Constant

  #---No changes
  if { !$Constant(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![Constant::checkPanelData] } {
    return
  }

  #---Save data
  Constant::panelSave
  Panel::cancel $w
} 


proc Constant::panelApply {} {
  global Constant

  #---No changes
  if { !$Constant(dataChanged) } {
    return
  }

  #---Error in data
  if { ![Constant::checkPanelData] } {
    return
  }

  Constant::panelSave
}


proc Constant::panelCancel {w} {
  global Info Constant

  if { ![Panel::verifyCancel Constant] } {
    return
  }

  #---Reset into old values
  Panel::restoreFields Constant

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc Constant::checkPanelData {} {
  global Info Constant Model

  #--Check gravity field
  set gravity_error 0

  set gravityVars {GRAVITY_1 GRAVITY_2 GRAVITY_3 GRAVITY_ABS}

  foreach var $gravityVars {
    
    set msg [DataField::checkValue res_value res_size $Constant($var) "%10.4f" "nbr" 1]

    if {$msg != ""} {
      set gravity_error 1
      
      set lbl $var

      set Info(messageIcon) error
      Message::dispOkMessage  [list  "$lbl\n\n" "Incorrect value!\n"]

    } else {
      set Constant($var) $res_value
    }
  }

  # If data ok, normalize vector
  if { !$gravity_error } {

    set vector [Util::normalizeVector $Constant(GRAVITY_1) \
                                      $Constant(GRAVITY_2) \
                                      $Constant(GRAVITY_3)]
 
    lappend vector $Constant(GRAVITY_ABS)

    set Constant(GRAVITY) ""

    foreach i {0 1 2 3} var $gravityVars {
      
      set Constant($var,prev) $Constant($var)
      set Constant($var,err,prev) $Constant($var,err)
      set Constant($var,mod,prev) $Constant($var,mod)
      
      set Constant($var) [lindex $vector $i]
      set Constant($var,err) 0
      set Constant($var,mod) 0 ;# All fields are recalculated!
      
      lappend Constant(GRAVITY) $Constant($var)

      if { $i == 2 && $Model(SIMULATION_DIMENSION) == "2D" } {
        continue
      }

      Widget::setFieldStatus Constant $var
    }
  }

  #--Check other constants
  set constant_error 0

  set constantVars STEFAN_BOLTZMANN

  foreach var $constantVars {

    set msg [DataField::checkValue res_value res_size $Constant($var) "%10.4e" "nbr" 1]

    if {$msg != ""} {
      set constant_error 1

      set lbl $var

      set Info(messageIcon) error
      Message::dispOkMessage "Incorrect value in field: $lbl" \
                     "$Info(FRONT_NAME) error message" \
                     $Constant(winName) 
    } else {
      set Constant($var) $res_value
    }
  }

  if {$gravity_error || $constant_error} {
    return 0
  }
 
  # Ok
  return 1
}


# Helper proc for handling gravity vector
# NOTE: These named compnents are used in screen handling
# and also when gravity is updated outside via
# Constant::panelSave proc!
#
proc Constant::formGravityComponents {} {
  global Constant Model

  # Helper variables
  set Constant(GRAVITY_1)   [lindex $Constant(GRAVITY) 0]
  set Constant(GRAVITY_2)   [lindex $Constant(GRAVITY) 1]
  set Constant(GRAVITY_3)   [lindex $Constant(GRAVITY) 2]
  set Constant(GRAVITY_ABS) [lindex $Constant(GRAVITY) 3]

  if { $Model(SIMULATION_DIMENSION) == "2D" } {
    set Constant(GRAVITY_3) 0
  }

  set vector [Util::normalizeVector $Constant(GRAVITY_1) \
                                    $Constant(GRAVITY_2) \
                                    $Constant(GRAVITY_3)]

  set Constant(GRAVITY_1)   [lindex $vector 0]
  set Constant(GRAVITY_2)   [lindex $vector 1]
  set Constant(GRAVITY_3)   [lindex $vector 2]
 
}

# end ecif_tk_constantDefPanel.tcl
# ********************
