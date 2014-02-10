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
#Module:    ecif_tk_equationVariablesPanel.tcl
#Language:  Tcl
#Date:      17.05.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining variables for specific equations
#
#************************************************************************


proc EquationVariable::openPanel { callingPanel } {
  ## This procedure displays the input and output files definition panel
  ## Global variables
  global Info Common Model EquationVariable Equation
  global platform
  
  if { $Equation(equationIndex) == "" } {
    return
  }

  set w $EquationVariable(winName)
  set wgeom $EquationVariable(winGeometry)

  set id [winfo atom $w]
  set EquationVariable(winId) $id
  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow EquationVariable $id $EquationVariable(winTitle) $wgeom] } {
    return
  }  

  set EquationVariable(dataChanged) 0
  set EquationVariable(dataModified) 0
  set EquationVariable(dataWasChanged) 0

  set this $w
  toplevel $w
  focus $w

  wm title $w $EquationVariable(winTitle)
  wm geometry $w $wgeom 

  Panel::resetFields EquationVariable

  # We init here in the case that some new fields
  # have been added by user definition files!
  #
  Panel::initFields EquationVariable

  set id $EquationVariable(parameterId)
  DataField::formDataFields EquationVariable $EquationVariable($id,data)
  

  Panel::backupFields EquationVariable

  set Info(equationVariablesPanelApplied) 0

 #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)

  #---Outer frames
  set f1 [frame $w.f1] ;#Title
  set f2 [frame $w.f2] ;#Help text
  set f3 [frame $w.f3] ;#Text entry
  set f4 [frame $w.f4] ;#Apply+Ok+cancel buttons frame
  #
  set ent_wid 30
  set ent_hig 16

  set eq_var [Equation::getFieldPropertyByIndex $Equation(equationIndex) EquationField]

  if { $eq_var == ""} {
    return
  }

  set EquationVariable(currentField) $eq_var

  set eq_lbl [DataField::getFieldProperty Equation $eq_var Label]

  # Label
  #
  label $f1.l  -text "$eq_lbl variables"
  pack $f1.l -side top

  set table_hdr    "\n\nEnter one variable name per line:"

  label $f2.l -text $table_hdr

  pack $f2.l -side top
  
  # Text widget
  #
  set wdg [text $f3.e  -width $ent_wid -height $ent_hig -font $Info(entryFont)]

  bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 EquationVariable $wdg {%A %K}"
  bind  $wdg <KeyPress-Return> "EquationVariable::panelCheck
  "
  pack $wdg -side top -padx $fpx2 -pady $fpy2 
  
  set EquationVariable(allWidgets,variableNames) $wdg

  # Buttons
  #

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $f4.ok -text OK -command "EquationVariable::panelOk $this"]
  set cn_btn [button $f4.cancel -text Cancel -command "EquationVariable::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f4.apply -text Apply -command "EquationVariable::panelApply" \
                               -state $ap]

  focus $ok_btn
  set EquationVariable(applyButton)  $ap_btn
  set EquationVariable(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1

  pack $f1 $f2 $f3 $f4 -side top -expand 1 -padx $fpx2 -pady $fpy2

  # Init panel
  DataField::data2TextWidget EquationVariable $EquationVariable(currentField) $wdg
}


proc EquationVariable::panelSave { {inform_front 1} {update_equation_panel 1} } {
  global Common EquationVariable Equation Model

  set panel $EquationVariable(parameterType)

  #--Update variables and components
  # NOTE: This must be done before backup!!!
  #
  set old_names [EquationVariable::getOldVariableNames]
  set added_names [EquationVariable::getAddedVariableNames]
  set removed_names [EquationVariable::getRemovedVariableNames]

  set ename $EquationVariable(currentField)
  set emask [DataField::getFieldProperty Equation $ename EquationMask]
  set edofs [DataField::getFieldProperty Equation $ename VariableDofs "" 0]

  foreach rn $removed_names {
    DataField::unsetVariableComponents $rn $edofs
  }

  foreach an $added_names {
    DataField::setVariableComponents $ename $an 
  }
  
  #Panel::initFields Variables

  #--Store old values 
  Panel::backupFields EquationVariable

  #--Form active parameter
  EquationVariable::formActiveParameter

  #--Update list of all variables EquationVariable(allVariableNames) 
  EquationVariable::updateAllVariableNames

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 EquationVariable
  Panel::panelDataModified 0 EquationVariable

  if { $update_equation_panel } {
    Equation::updateEquationVariables $old_names $added_names $removed_names
  }

  #--Set panel states
  Panel::panelDataModified 1 Equation 
  Panel::panelDataChanged 1 Equation 
  set EquationVariable(dataWasChanged) 1

  Util::cpp_exec equationVariablesPanelOk
}


proc EquationVariable::panelOk {w} {
  global EquationVariable
  
  #---No changes
  if { !$EquationVariable(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![EquationVariable::panelCheck 1] } {
    return
  }

  EquationVariable::panelSave
  Panel::cancel $w
}


proc EquationVariable::panelApply {} {
  global EquationVariable

  #---No changes
  if { !$EquationVariable(dataChanged) } {
    return
  }

  #---Error in data
  if { ![EquationVariable::panelCheck 1] } {
    return
  }

  EquationVariable::panelSave
}


proc EquationVariable::panelCancel {w} {
  global EquationVariable

  if { ![Panel::verifyCancel EquationVariable] } {
    return
  }

  #-Restore old values
  Panel::restoreFields EquationVariable

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc EquationVariable::panelCheck { {trim_empty_lines 0} } {
  global EquationVariable Model Info

  set is_ok 1

  #--Check widget data
  set ew $EquationVariable(allWidgets,variableNames)
  set data [$ew get 0.0 end]

  set fld $EquationVariable(currentField)

  #--Pick widget data
  DataField::textWidget2Data EquationVariable $fld $ew $trim_empty_lines

  #--Form data list
  set data [split $EquationVariable($fld) $Info(dataListSeparator)]

  if { $trim_empty_lines } {
    DataField::data2TextWidget EquationVariable $fld $ew
  }

  #--Check that no spaces or duplicate names
  # NOTE: duplicate checking is non-case sensitive because
  # Fortran does not differentiate case in variable names
  #
  set index 1
  foreach d1 $data {

    set d1 [Util::stringTrim $d1]

    if { -1 != [string first " "  $d1] ||
         -1 != [string first "\t" $d1]
       } {
      #set msg "Space is not allowed in a variable name!  ($d1)"

      #set Info(messageIcon) error
      #Message::dispOkMessage $msg
      #set is_ok 0
      #break
    }

    set data2 [lrange $data $index end]
    incr index

    foreach d2 $data2 {

      set d2 [Util::stringTrim $d2]

      if { [string tolower $d1] == [string tolower $d2] } {
        set msg "Duplicate variable name: $d1 $d2"
        set Info(messageIcon) error
        Message::dispOkMessage $msg
        set is_ok 0
      }
    }
  }

  #--Result
  return $is_ok
}


# ============
# Helper procs
# ============

# We have currently only one EquationVariable parameter!!!
#
proc EquationVariable::formActiveParameter {} {
  global EquationVariable

  set EquationVariable(ids) 1
  DataField::formNonStandardParameter EquationVariable 1 "EquationVariable1"
}


# Find old variable names
#
proc EquationVariable::getOldVariableNames {} {
  global EquationVariable

  set ename $EquationVariable(currentField)
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]

  set old_names [split $EquationVariable($ename,old) $sep]

  return $old_names
}


# Find variable names that were added
#
proc EquationVariable::getAddedVariableNames {} {
  global EquationVariable

  set new_names ""

  set ename $EquationVariable(currentField)
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]

  set old_names [split $EquationVariable($ename,old) $sep]
  set cur_names [split $EquationVariable($ename) $sep]

  foreach cur_name $cur_names {

    if { "" == [Equation::getMatchingVariableName $cur_name $old_names] } {
      lappend new_names $cur_name
    }
  }

  return $new_names
}


# Find variable names that were removed
#
proc EquationVariable::getRemovedVariableNames {} {
  global EquationVariable

  set removed_names ""

  set ename $EquationVariable(currentField)
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]

  set old_names [split $EquationVariable($ename,old) $sep]
  set cur_names [split $EquationVariable($ename) $sep]

  foreach old_name $old_names {

    if { "" == [Equation::getMatchingVariableName $old_name $cur_names] } {
      lappend removed_names $old_name
    }
  }

  return $removed_names
}


# Check that all equation variable names are store in the
# equation-name fields
# NOTE: This problem could occur if a user definition changes
# the name in an existing multi-var equation (all user defined
# equations are by multi-var default and has an initial value for the
# equation variable!)
#
proc EquationVariable::updateFieldData {} {
  global Equation EquationVariable

  foreach eq_name $Equation(allFields) {

    if { ![Equation::isEquation $eq_name] } {
      continue
    }

    set vn [DataField::getFieldProperty EquationVariable $eq_name InitialValue]

    if { $vn == "" } {
      continue
    }
    
    # Add to the existing field (if not there)
    #
    if { [info exists EquationVariable($eq_name)] } {

      set sep [DataField::getFieldProperty EquationVariable $eq_name FieldDataSep]
      set old_list [split $EquationVariable($eq_name) $sep]
      
      # This belongs always to the list!
      set EquationVariable($eq_name) $vn

      # Add each of the olds which is different from '$vn', to the list
      #
      foreach old $old_list {

        if { ![string equal -nocase $vn $old] } {
          append EquationVariable($eq_name) "$sep$old"
        }
      }
   
    }
  }
}


# Update the list of all variable names
#
proc EquationVariable::updateAllVariableNames {} {
  global EquationVariable

  set EquationVariable(allVariableNames) ""

  foreach fld $EquationVariable(allFields) {

    set sep [DataField::getFieldProperty EquationVariable $fld FieldDataSep]

    set names [split $EquationVariable($fld) $sep]

    foreach name $names {
      lappend EquationVariable(allVariableNames) $name
    }
  }
}


#End ecif_equationVariables.tcl

