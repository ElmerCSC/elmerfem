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
#Module:    ecif_tk_standardPanelCheck.tcl
#Language:  Tcl
#Date:      05.10.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Standard panels checking procedures.
#
#************************************************************************

# ========================== #
# Panel specific check procs #
# ========================== #

# *fillProcs:
# ==========
# are run when field data is filled after selecting a parameter etc.

# *editProcs:
# ==========
# are run when field's value is modified on the screen

# *checkProcs :
# ===========
# are run when field's value is being checked # before add or update.
# These fields are normally not on the screen, so there is no editProc for them!
#
# NOTE: checkProcs are called from procedure PanelCheck::evalFieldCheckProc
# in procedure StdPanelCheck::checkFieldValue
#
# NOTE: These field should located after the fields they are using for checking the
# value in the array-name(allFields) list! An example is the Calculator PROCEDURE
# field which is built from LIBRARY_FILE and FUNCTION_NAME fields
#
# NOTE: These procs should return "" when ok, otherwise some error message!!!


################################
### BOUNDARY CONDITION panel ###
################################

# Check if argument parameter has one common parent body in all its atachemnts
# This is needed eg. when deciding if DiffuseGray button should be active
# for the paramter
#
proc BoundaryCondition::parameterHasCommonParent {pid} {
  global BoundaryCondition Info ObjectTable

  set counter 0
  set parent_ids ""
  set attach_count ""

  foreach oid $ObjectTable(ids) {

    if { ![info exists ObjectTable($oid,bc)] ||
         $ObjectTable($oid,bc) != $pid
       } {
      continue
    }

    set pr_id $ObjectTable($oid,prId)

    # No parent: vertex, edge in 3D
    if { $pr_id == $Info(NO_INDEX) } {
      continue
    }

    set pr1_id $ObjectTable($pr_id,pr1Id)
    set pr2_id $ObjectTable($pr_id,pr2Id)

    incr counter

    if { -1 == [lsearch $parent_ids $pr1_id] } {
      lappend parent_ids $pr1_id
    }

    if { -1 == [lsearch $parent_ids $pr2_id] } {
      lappend parent_ids $pr2_id
    }

    if { ![info exist attach_count$pr1_id] } {
      set attach_count$pr1_id 1
    } else {
      incr attach_count$pr1_id
    }

    if { ![info exist attach_count$pr2_id] } {
      set attach_count$pr2_id 1
    } else {
      incr attach_count$pr2_id
    }
  }

  if { [llength $parent_ids] == 0 } {
    return 1
  }

  foreach id $parent_ids {

    set value [set attach_count$id]


    if { $value >= $counter } {
      return 1
    }
  }

  return 0
}


# Check if a given boundary condition accepts current target
# Problem here is the radiation target body matching!
# NOTE: Target body id is an absolute id in the sif-file. That's why
# we must be careful with targets!
#
proc BoundaryCondition::parameterTargetCheckProc {paramter_id {accept_var ""} } {
  global BoundaryCondition Info ObjectTable

  if { $accept_var != "" } {
    upvar #1 $accept_var accept
  }

  set accept 1

  set pid $paramter_id

  set mode [DataField::getFieldValue BoundaryCondition $pid RADIATION]

  if { $mode == "" || ![string equal -nocase $mode "Diffuse gray"] } {
    return $accept
  }

  # Diffusegray radiation defined --> target bodies must match!
  #
  #set oid $BoundaryCondition($pid,oid)
  #set parent_id $ObjectTable($oid,prId)

  set target_body [DataField::getFieldValue BoundaryCondition $pid RADIATION_TARGET_BODY]

  set oid $BoundaryCondition(objectId)
  set bd1_id $ObjectTable($oid,pr1Id)
  set bd2_id $ObjectTable($oid,pr2Id)

  set bd1_tag $ObjectTable($bd1_id,tg)
  set bd2_tag $ObjectTable($bd2_id,tg)

  # Body2 is "Outside"
  if { $bd1_tag == $bd2_tag } {
    set bd2_tag $Info(NO_INDEX)
  }

  #if { $parent_id != $BoundaryCondition(objectId) }

  if { $target_body == $bd1_tag ||
       $target_body == $bd2_tag
     } {
    
    return $accept
  }

  set accept 0


  return $accept
}


proc BoundaryCondition(RADIATION,fillProc) {} {

  BoundaryCondition(RADIATION,groupEditProc)
}


proc BoundaryCondition(RADIATION,groupEditProc) { {fld ""} } {
  global BoundaryCondition Model ObjectTable

  set BoundaryCondition(EMISSIVITY,act) 0
  set BoundaryCondition(RADIATION_TARGET_BODY_NAME,act) 0
  set BoundaryCondition(RADIATION_TARGET_BODY,act) 0

  set value $BoundaryCondition(RADIATION)

  if { $value != "None" && $value != "" } {
    set BoundaryCondition(EMISSIVITY,act) 1
  }

  # Diffuse gray radiation
  # ======================
  if { $value != "" && [string equal -nocase $value "Diffuse gray"] } {

    set BoundaryCondition(RADIATION_TARGET_BODY_NAME,act) 1
    set BoundaryCondition(RADIATION_TARGET_BODY,act) 1
  }

  if { ![Widget::arrayWidgetExists BoundaryCondition RD_DIFFUSE_GRAY] } {
    return
  }

  if { ![BoundaryCondition::parameterHasCommonParent $BoundaryCondition(parameterId)] } {
    set BoundaryCondition(RD_DIFFUSE_GRAY,act) 0
  }

  Widget::configureField BoundaryCondition RD_DIFFUSE_GRAY
  Widget::configureField BoundaryCondition EMISSIVITY

  # NOTE: Name-list widget state is handled in its own check-proc
  # because widget must be rebuilt if name list changes!
  BoundaryCondition(RADIATION_TARGET_BODY,fillProc)
  BoundaryCondition(RADIATION_TARGET_BODY_NAME,fillProc)
}


# This proc puts the radiation target body name-list in the OptionMenu
# widget according to the selected boundary
#
proc BoundaryCondition(RADIATION_TARGET_BODY_NAME,fillProc) {} {
  global BoundaryCondition Info ObjectTable

  if { ![Widget::arrayWidgetExists BoundaryCondition RADIATION_TARGET_BODY_NAME] } {
    return
  }

  set fld RADIATION_TARGET_BODY_NAME

  set wdg $BoundaryCondition(allWidgets,$fld)

  set bd1 $BoundaryCondition(body1Id)
  set bd2 $BoundaryCondition(body2Id)

  if {$bd1 == "" || $bd1 == $Info(NO_INDEX) } {
    return
  }

  set body1 "($ObjectTable($bd1,tg)) "
  append body1 $ObjectTable($bd1,nm)

  if { $bd2 != $Info(NO_INDEX) } {
    set body2 "($ObjectTable($bd2,tg)) "
    append body2 $ObjectTable($bd2,nm)

  } else {
    set body2 "Outside"
  }

  set names [list $body1 $body2]
  
  if { $BoundaryCondition(RADIATION) != "" &&   
       [string equal -nocase $BoundaryCondition(RADIATION) "Diffuse gray"]
     } {
    set wdg_state normal
  } else {
    set wdg_state disabled
  }

  # Pick current MenuSelect proc
  set ms_proc [ bind $BoundaryCondition(allWidgets,$fld,om) <<MenuSelect>>]

  destroy $wdg

  # Build new option-menu with current parent names
  set om [Panel::createOptionMenuWidget $wdg BoundaryCondition($fld) $names]
  $wdg configure -indicatoron 0 -width 9 
  #$wdg configure -activeforeground black
  #$wdg configure -disabledforeground darkGray
  Widget::configureS $wdg $wdg_state

  BoundaryCondition(RADIATION_TARGET_BODY_NAME,checkWidgetState)

  bind $om <<MenuSelect>> $ms_proc

  # Save and pack rebuilt widget
  set BoundaryCondition(allWidgets,$fld) $wdg
  pack $wdg

}


# This proc picks the user selected name value from the OptionMenu
# widget body-name variable and then sets corresponding value 
# to the target-body variable
#
proc BoundaryCondition(RADIATION_TARGET_BODY_NAME,editProc) {} {
  global BoundaryCondition Info ObjectTable

  set name $BoundaryCondition(RADIATION_TARGET_BODY_NAME)
  set name_prev $BoundaryCondition(RADIATION_TARGET_BODY_NAME,prev)
  set BoundaryCondition(RADIATION_TARGET_BODY_NAME,prev) $name

  if { $name != $name_prev } {
    Panel::panelDataModified 1 BoundaryCondition
  }    
  
  #-Outside
  if { $name == "Outside" } {
    set bd_tag $Info(NO_INDEX)
  
  #-Pick body tag
  } else {
    set bd_tag [string trim [lindex [split $name ")"] 0] "("]
  }

  set BoundaryCondition(RADIATION_TARGET_BODY) $bd_tag
}


# This proc checks if the radiation target body widget must be disabled to prevent
# changing the target body
#
proc BoundaryCondition(RADIATION_TARGET_BODY_NAME,checkWidgetState) { {force_state 0} {value 0} } {
  global BoundaryCondition ObjectTable

  set wdg $BoundaryCondition(allWidgets,RADIATION_TARGET_BODY_NAME)
  
  if { ![winfo exists $wdg] } {
    return
  }

  #-Calculate nof parent bodies for the parameter attachments
  #
  set parent_ids ""
  set pid $BoundaryCondition(parameterId)

  foreach oid $ObjectTable(ids) {

    if { [info exists ObjectTable($oid,bc)] &&
         $ObjectTable($oid,bc) == $pid
       } {
      
      if { -1 == [lsearch $parent_ids $ObjectTable($oid,prId)] } {
        lappend parent_ids $ObjectTable($oid,prId)
      }
    }
  }

  # If no attachments or attachments only for the boundaries of the
  # current body, we can allow changing the target bodies, ie. activate
  # the widget
  #   
  if {  0 == [llength $parent_ids] ||
       (1 == [llength $parent_ids] &&
        $BoundaryCondition(objectId) == [lindex $parent_ids 0]
       ) } {
    return
  }

  Widget::configureS $wdg "disabled"
}


# This proc puts the correct body-name value to be displayed in the
# OptionMenu widget, based on target-body value, when a boundary
# or boundary condition is selected
# NOTE: Current version with only one target body in the result list!
#
proc BoundaryCondition(RADIATION_TARGET_BODY,fillProc) {} {
  global BoundaryCondition Info ObjectTable

  if { ![Widget::arrayWidgetExists BoundaryCondition RADIATION_TARGET_BODY_NAME] } {
    return
  }


  set oid $BoundaryCondition(objectId)

  set bd1_id $ObjectTable($oid,pr1Id)
  set bd2_id $ObjectTable($oid,pr2Id)

  set bd1_tag $ObjectTable($bd1_id,tg)
  set bd2_tag $ObjectTable($bd2_id,tg)

  # Second 'body' is 'Outside'
  if { $bd2_tag == $bd1_tag } {
    set bd2_id $Info(NO_INDEX)
    set bd2_tag $Info(NO_INDEX)
  }


  # Select correct body tag for the target body
  #
  set bd_tag $BoundaryCondition(RADIATION_TARGET_BODY)

  # First of the body pair is ok
  if { $bd_tag == $bd1_tag } {
    set bd_id $bd1_id

  # Second of the body pair is ok
  } elseif { $bd_tag == $bd2_tag } {
    set bd_id $bd2_id

  # We are adding: must be one of the bodies!
  } elseif { $bd_tag == "" || $bd_tag == $Info(NO_INDEX) || $BoundaryCondition(mode) == "ADD" } {
    set bd_id $bd1_id
    set bd_tag $bd1_tag

  # An existing paramter was selected, we use its target body tag
  } else {
      set bd_id [Object::findIdByTag "B" $bd_tag]
  }

  set BoundaryCondition(RADIATION_TARGET_BODY) $bd_tag

  # Set target body name 
  # =====================
  #-Outside
  if { $bd_tag == $Info(NO_INDEX) } {  
    set name "Outside"

  #-Body tag and name
  } else {
    set name "($bd_tag) "
    append name $ObjectTable($bd_id,nm)
  }

  set BoundaryCondition(RADIATION_TARGET_BODY_NAME) $name
  set BoundaryCondition(RADIATION_TARGET_BODY_NAME,prev) $name
}


# This proc puts the correct body-name value to be displayed in the
# OptionMenu widget, based on target-body value, when a boundary
# or boundary condition is selected
# NOT IN USE, MVe 01.08.00!
#
proc BoundaryCondition_experimental(RADIATION_TARGET_BODIES,fillProc) {} {
  global BoundaryCondition Info ObjectTable

  if { ![Widget::arrayWidgetExists BoundaryCondition RADIATION_TARGET_BODY_NAME] } {
    return
  }

  set oid $BoundaryCondition(objectId)

  set bd1_id $ObjectTable($oid,pr1Id)
  set bd2_id $ObjectTable($oid,pr2Id)

  set bd1_tag $ObjectTable($bd1_id,tg)
  set bd2_tag $ObjectTable($bd2_id,tg)

  set bndr_tag $ObjectTable($BoundaryCondition(boundaryId),tg)

  # Pre version 5 file conversion!
  if { $BoundaryCondition(RADIATION_TARGET_BODIES) == "" &&
       [info exists BoundaryCondition(RADIATION_TARGET_BODY)] &&
       $BoundaryCondition(RADIATION_TARGET_BODY) != ""
     } {

    set bd_tag $BoundaryCondition(RADIATION_TARGET_BODY)

    lappend BoundaryCondition(RADIATION_BOUNDARIES) $bndr_tag
    lappend BoundaryCondition(RADIATION_TARGET_BODIES) $bd_tag
  }

  set idx [lsearch $BoundaryCondition(RADIATION_BOUNDARIES) $bndr_tag]

  #--If not yet stored
  if { $idx == -1 } {

    set bd_tag $BoundaryCondition(RADIATION_TARGET_BODY)

    if { $bd_tag == $Info(NO_INDEX) } {
      
      if { $ObjectTable($oid,tp) == "B" } {
        set bd_id $Info(NO_INDEX)

      } else {
        set bd_id $bd2_id
      }

    } elseif { $bd_tag == $bd1_tag } {
      set bd_id $bd1_id

    } elseif { $bd_tag == $bd2_tag } {
      set bd_id $bd2_id

    } else {
      set bd_tag $bd1_tag
      set bd_id $bd1_id
    }

  # Check parent body by boundary's index
  } else {

    set bd_tag [lindex $BoundaryCondition(RADIATION_TARGET_BODIES) $idx]

    if { $bd_tag != $Info(NO_INDEX) } {
      set bd_id [Object::findIdByTag "B" $bd_tag]
    }
  }

  set BoundaryCondition(RADIATION_TARGET_BODY) $bd_tag

  # Set target body name 
  # =====================
  #-Outside
  if { $bd_tag == $Info(NO_INDEX) } {  
    set name "Outside"

  #-Body tag and name
  } else {
    set name "($bd_tag) "
    append name $ObjectTable($bd_id,nm)
  }

  set BoundaryCondition(RADIATION_TARGET_BODY_NAME) $name
  set BoundaryCondition(RADIATION_TARGET_BODY_NAME,prev) $name
}



# Checks and marks if there are some force related flow boundary conditions
#
proc BoundaryCondition(FLOW_FORCE_BC,fillProc) {} {
  global BoundaryCondition

  set value ""

  if { $BoundaryCondition(EXTERNAL_PRESSURE) != ""
       ||
       ( $BoundaryCondition(PRESSURE_1) != ""  ||
         $BoundaryCondition(PRESSURE_2) != ""  ||
         $BoundaryCondition(PRESSURE_3) != ""  
       )
       ||
       $BoundaryCondition(SURFACE_TENSION_COEFFICIENT) != "" 
       ||
       ( [info exists BoundaryCondition(WALL_LAW)] &&
         $BoundaryCondition(WALL_LAW) == "True"
       )
     } {
    set value "True"
  }

  set BoundaryCondition(FLOW_FORCE_BC) $value
}


# Checks and marks if there are some heat flux related boundary conditions
#
proc BoundaryCondition(HEAT_FLUX_BC,fillProc) {} {
  global BoundaryCondition

  set value ""

  if { ($BoundaryCondition(RADIATION) == "Idealized" ||
        $BoundaryCondition(RADIATION) == "Diffuse gray"
       ) 
       ||
       $BoundaryCondition(HEAT_FLUX) != ""
       ||
       ($BoundaryCondition(HEAT_TRANSFER_COEFFICIENT) != "" &&
        $BoundaryCondition(EXTERNAL_TEMPERATURE) != ""
       )
     } {
    set value "True"
  }

  set BoundaryCondition(HEAT_FLUX_BC) $value
}


proc BoundaryCondition(NORMAL-TANGENTIAL_DISPLACEMENT,fillProc) {} {

  BoundaryCondition(NORMAL-TANGENTIAL_DISPLACEMENT,editProc)
}

proc BoundaryCondition(NORMAL-TANGENTIAL_DISPLACEMENT,editProc) {} {
  global BoundaryCondition

  if { ![Widget::arrayWidgetExists BoundaryCondition NORMAL-TANGENTIAL_DISPLACEMENT] } {
    return
  }

  StdPanelExec::setNormalTangentialLabels \
    BoundaryCondition \
    DISPLACEMENT \
    $BoundaryCondition(NORMAL-TANGENTIAL_DISPLACEMENT)
}


proc BoundaryCondition(NORMAL-TANGENTIAL_VELOCITY,fillProc) {} {

  BoundaryCondition(NORMAL-TANGENTIAL_VELOCITY,editProc)
}


proc BoundaryCondition(NORMAL-TANGENTIAL_VELOCITY,editProc) {} {
  global BoundaryCondition

  if { ![Widget::arrayWidgetExists BoundaryCondition NORMAL-TANGENTIAL_VELOCITY] } {
    return
  }

  StdPanelExec::setNormalTangentialLabels \
    BoundaryCondition \
    VELOCITY \
    $BoundaryCondition(NORMAL-TANGENTIAL_VELOCITY)
}


#-Procedure checks velocity-group field data
# Not in use!
#
proc BoundaryCondition(VELOCITY,checkProc) { } {
  return [Common(VectorGroup,checkProc) BoundaryCondition VELOCITY]
}



########################
### CALCULATOR panel ###
########################

proc Calculator(ACTIVE,fillProc) {} {
  Calculator(ACTIVE,editProc)
}

proc Calculator(ACTIVE,editProc) {} {
  global Calculator

  if { $Calculator(ACTIVE) == "True" } {
    set active 1
  } else {
    set active 0
  }

  set flds { INCLUDE,Use
             EQUATION LIBRARY_FILE FUNCTION_NAME
             VARIABLE VARIABLE_DOFS
           }

  foreach fld $flds {
    set Calculator($fld,act) $active
    Widget::configureField Calculator $fld
  }

  Widget::configureField Calculator INCLUDE

}


proc Calculator(EQUATION,fillProc) {} {
  global Calculator

  #set Calculator(parameterName) [string trim $Calculator(EQUATION)]
}


proc Calculator(PROCEDURE,fillProc) {} {
  global Calculator

  set vals [split $Calculator(PROCEDURE) ";"]

  set Calculator(LIBRARY_FILE) [lindex $vals 0]
  set Calculator(FUNCTION_NAME) [lindex $vals 1]
}

proc Calculator(PROCEDURE,checkProc) {} {
  global Calculator

  set lb $Calculator(LIBRARY_FILE)
  set fn $Calculator(FUNCTION_NAME)

  set Calculator(PROCEDURE) [join [list $lb $fn] ";"]

  return ""
}


proc Calculator(VARIABLE_DOFS_ALL,checkProc) {} {
  global Calculator Model

  set values $Calculator(VARIABLE_DOFS_ALL)
  
  if { 2 == [llength $values] } {
    set values [linsert $values 0 1]
  }

  #-No dimension dependence
  #
  if { 1 == [llength $values] } {
    set dof $values

  #-Depends on simulation dimension
  #
  } elseif { $Model(SIMULATION_DIMENSION) == "1D" } {
    set dof [lindex $values 0]

  } elseif { $Model(SIMULATION_DIMENSION) == "2D" } {
    set dof [lindex $values 1]

  } else {
    set dof [lindex $values end]
  }

  set Calculator(VARIABLE_DOFS) [Calculator::getVariableDofs $values]

  return ""
}



######################
### EQUATION panel ###
######################

proc Equation(EQUATION_LIST,fillProc) {} {
  global Equation EquationVariable Model

  Equation::updateEquationList 
  
  if { $Model(GEOMETRY_DIMENSION) == "" ||
       ![info exists Equation(allWidgets,EQUATION_LIST)] ||
       ![winfo exists $Equation(allWidgets,EQUATION_LIST)]
     } {
    return
  }

  set sindex [lsearch $Equation(ownIndices) $Equation(equationIndex)]
  $Equation(allWidgets,EQUATION_LIST) selection set $sindex

  Equation::setVariableWidgetsActivity $Equation(equationIndex)
  Equation::updateVariableList $Equation(equationIndex)
}


proc Equation(EQUATION_LIST,editProc) { arguments } {
  global Equation

  #-Pick arguments
  set event_mode [lindex $arguments 0] ;# Button and keys pressed
  set lb_x [lindex $arguments 1] ;# Local lb x-coordinate for selection
  set lb_y [lindex $arguments 2] ;# Local lb y-coordinate for selection 

  set current_row ""
  set picked_row ""
  set sindex ""

  set current_row [$Equation(allWidgets,EQUATION_LIST) curselection]

  #-Selected ("clicked") equation row
  #
  if { $lb_x != "" && $lb_y != "" } {
    set picked_row [$Equation(allWidgets,EQUATION_LIST) index @$lb_x,$lb_y]
  } else {
    set picked_row $current_row
  }


  if { $current_row == "" && $picked_row == "" } {
    return
  }

  if { $event_mode == "ButtonRelease-1" &&
       ($lb_x != "" && $lb_x < 25)
     } {
    set event_mode "ButtonRelease-3"
  }

  #-Btn-1: Equation change only
  # ==================================
  if { $event_mode == "ButtonRelease-1" } {
    set select_row $picked_row
    set toggle_row ""
    set sindex $picked_row
    
  #-Btn-3: Equation change+toggle
  # ==================================
  } elseif { $event_mode == "ButtonRelease-3" } {
    set select_row $picked_row
    set toggle_row $picked_row
    set sindex $picked_row


  #-Ctrl-Btn-1: Toggle current row
  # ==================================
  } elseif { $event_mode == "Control-ButtonRelease-1" } {
    set select_row ""
    set toggle_row $current_row
    set sindex $current_row

  #-Ctrl-Btn-3: Toggle picked row
  # ==================================
  } elseif { $event_mode == "Control-ButtonRelease-3" } {
    set select_row ""
    set toggle_row $picked_row
    set sindex $current_row
  }
  
  # Select equation
  #
  if { $select_row != "" } {

    set ename [lindex $Equation(ownNames) $select_row]
    set eindex [lindex $Equation(ownIndices) $select_row]

    set Equation(equationLabel) [Equation::getFieldPropertyByIndex $eindex EquationLabel]
    StdPanelExec::equationMenuSelectionProc Equation $eindex
  }

  # Toggle equation on/off value
  #
  if { $toggle_row != "" } {

    set ename [lindex $Equation(ownNames) $toggle_row]
    set eindex [lindex $Equation(ownIndices) $toggle_row]

    # Toggle value
    #
    if { $Equation($ename) == "True" } {
      set Equation($ename) "False"
    } else {
      set Equation($ename) "True"
    }

    #-Pick current top row index in the listbox
    set  ix [ $Equation(allWidgets,EQUATION_LIST) nearest 0]

    #-Update equation list itself (X-markers!)
    Equation::updateEquationList

    # Keep same rows in the listbox
    $Equation(allWidgets,EQUATION_LIST) see [expr 7 + $ix]

    #-Update variable list stuff
    if { $Equation(equationIndex) == $eindex } {
      Equation::setVariableWidgetsActivity $eindex
    }

    Equation::updateVariableList $eindex

    Panel::panelDataModified 1 Equation
  }


  #-Some further checking
  #
  if { "" != [info commands Equation($ename,editProc) ] } {
    Equation($ename,editProc)

  } else {
    if { 1 == [DataField::getFieldProperty Equation $ename IsUserDefinedField "" 0] } {
      Equation::setUserDefinedEquationMask $ename
    }
  }

  # Update Equation panel field activity
  StdPanelExec::setValuesAreaActivity Equation $Equation(targetMask)
  
  if { $sindex != "" } {
    $Equation(allWidgets,EQUATION_LIST) selection set $sindex
    $Equation(allWidgets,EQUATION_LIST) activate $sindex
  }
}


proc Equation(VARIABLE_LIST,editProc) {arguments} {
  global Equation EquationVariable

  if { !$Equation(VARIABLE_LIST,act) } {
    return
  }

  set event_mode [lindex $arguments 0] ;# Button an d keys pressed
  set lb_x [lindex $arguments 1] ;# Local lb x-coordinate for selection
  set lb_y [lindex $arguments 2] ;# Local lb y-coordinate for selection 

  set sindex [$Equation(allWidgets,VARIABLE_LIST) index @$lb_x,$lb_y]

  if { $sindex == "" || $sindex < 0 } {
    return
  } 

  set eindex $Equation(equationIndex)
  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  #-All possible variable names

  set all_names [Equation::getAllVariableNames $eindex]

  #-Find currently active variable names
  set act_names [Equation::getActiveVariableNames $eindex]

  #-Find selected variable and "toggle" its state
  set sel_name [lindex $all_names $sindex]

  set idx [List::nocaseSearch $act_names $sel_name]

  #-Turn on (= add to the actives' list)
  if { $idx == -1 } {
    lappend act_names $sel_name

  #-Turn off (= delete from the actives' list)
  } else {
    set act_names [lreplace $act_names $idx $idx]
  }

  #-Sort names in actives-list by the names order in all-list
  #
  set sorted_list ""

  foreach name $all_names {

    if { -1 != [List::nocaseSearch $act_names $name] } {
      lappend sorted_list $name
    }
  }

  #-Update actives-list variable
  set act_names_var [DataField::getFieldProperty Equation $ename EquationVars]
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]
  set act_names [join $sorted_list $sep]

  set Equation($act_names_var) $act_names

  #-Update variable-list widget
  Equation::updateVariableList $eindex
}


# NOTE: Be careful with the (...,groupEditProc)-procs
# because they are also called in the "fillGroupFromDispList"
# (..procsDataField.tcl) proc when an equation parameter is
# selected in the equation list box!!!

#======================================#
#     Generic UserDefinedField         #
#======================================#

proc Equation(UserDefinedField,fillProc) { fld } {
 
  if { "" != [DataField::getFieldProperty Equation $fld EquationIndex "" 0] } {
    Equation::setUserDefinedEquationMask $fld
  }
}


proc Equation(UserDefinedField,editProc) { fld } {

  if { "" != [DataField::getFieldProperty Equation $fld EquationIndex "" 0] } {
    Equation::setUserDefinedEquationMask $fld
  }
}


#===========================================#
#                 FLOW equation             #
#===========================================#

# NOTE: This is called always when Navier-Stokes panel is activated
#
proc Equation(NAVIER-STOKES,fillProc) {} {
  global Equation

  #-We can not update values area now (loop!), so
  # store old value befor setting ther global var
  #
  set flag $Equation(updateValuesArea)
  set Equation(updateValuesArea) 0

  Equation(NAVIER-STOKES,editProc)
  Equation(TURBULENCE_MODEL,groupEditProc)
  
  set Equation(updateValuesArea) $flag
}


# NOTE: This is called when Navier-Stokes box is clicked
#
proc Equation(NAVIER-STOKES,editProc) {} {
  global Equation

  #-If no flow 
  if {$Equation(NAVIER-STOKES) != "True"} {
    set Equation(flow_type) ""
    set Equation(flow_laminarType) ""
    set Equation(flow_turbulentType) ""
    set Equation(TURBULENCE_MODEL) ""
    set Equation(HYDROSTATIC_PRESSURE) "False"

  #-When flow is checked on, default is laminar (NS) flow
  } else {
    set Equation(flow_type) "了L"
    set Equation(flow_laminarType) "了L1"
  }

  #-Broadcast type change to other
  Equation::updateTargetMask
}


#-Group (container) widget Turbulence Model widgets
# NOTE: This is grouper variable, it is not meant to be output
# into Solver-input file, but we need it to set the value
# for variables like like KE_TURBULENCE which are not on screen
#
proc Equation(TURBULENCE_MODEL,groupEditProc) { {fld ""} } {
  global Equation

  # Reset all group's target output variables
  set Equation(flow_turbulentType) ""
  set Equation(KE_TURBULENCE) "False"

  # Is some flow
  if { $Equation(flow_type) != "" } {

    switch $Equation(TURBULENCE_MODEL) {
      "" {
        set Equation(flow_turbulentType) ""
      }
      "KE" {
        set Equation(flow_turbulentType) "了T1"
        set Equation(KE_TURBULENCE) "True"
      }
    }
  }

  #-Broadcast type change
  Equation::updateTargetMask
}



#=======================================#
#               HEAT equation           #
#=======================================#

# NOTE: This is called when Heat equation panel is activated
#
proc Equation(HEAT_EQUATION,fillProc) {} {
  global Equation

  #-We can not update values area now (loop!), so
  # store old value befor setting ther global var
  #
  set flag $Equation(updateValuesArea)
  set Equation(updateValuesArea) 0

  Equation(HEAT_EQUATION,editProc)
  Equation(CONVECTION,groupEditProc)
  Equation(PHASE_CHANGE_MODEL,groupEditProc)

  set Equation(updateValuesArea) $flag
}


# NOTE: This is called when Heat Equation check box is clicked
#
proc Equation(HEAT_EQUATION,editProc) {} {
  global Equation

  #-If no heat
  if {$Equation(HEAT_EQUATION) != "True"} {

    set Equation(heat_type) ""
    set Equation(heat_transferType) ""
    set Equation(heat_phaseType) ""

    set Equation(PHASE_CHANGE_MODEL) ""
    set Equation(CHECK_LATENT_HEAT_RELEASE) ""

    # If no adv-diff, disable convetion fields
    # ========================================
    if { $Equation(ADVECTION_DIFFUSION_EQUATION) != "True" } {
      set Equation(CONVECTION) "None"
    }

  #-Heat equation was selected, we have at least:
  } else {
    set Equation(heat_type) "人"
    set Equation(heat_transferType) "人T1"
  }

  #-Broadcast type change
  Equation::updateTargetMask
}


#-Convection group (container) widget
#
proc Equation(CONVECTION,groupEditProc) { {fld ""} } {
  Equation(CONVECTION,heatCheckProc)  
  Equation(CONVECTION,diffusionCheckProc)  
}


#-Convection heat tranfer part check proc
#
proc Equation(CONVECTION,heatCheckProc) { {fld ""} } {
  global Equation

  if { $Equation(heat_transferType) == "" } {
    return
  }

  if { ![info exists Equation(CONVECTION)] } {
    return
  }

  switch $Equation(CONVECTION) {
    "None" {
      set Equation(heat_transferType) "人T1"
    }
    "Constant" {
      set Equation(heat_transferType) "人T1 人T2"
    }
    "Computed" {
      set Equation(heat_transferType) "人T1 人T3"
    }
  }
    
  #-Broadcast type change
  Equation::updateTargetMask
}


proc Equation(PHASE_CHANGE_MODEL,groupEditProc) { {fld ""} } {
  global Equation

  switch $Equation(PHASE_CHANGE_MODEL) {
    "None" {
      set Equation(heat_phaseType) ""
    }
    "Spatial 1" {
      set Equation(heat_phaseType) "人P1"
    }
    "Spatial 2" {
      set Equation(heat_phaseType) "人P2"
    }
    "Temporal" {
      set Equation(heat_phaseType) "人P3"
    }
  }

  #-Broadcast change
  Equation::updateTargetMask
}



#=========================================#
#                STRESS analysis          #
#=========================================#

# NOTE: This is called when Stress Analysis panel is activated
#
proc Equation(STRESS_ANALYSIS,fillProc) {} {
  global Equation

  #-We can not update values area now (loop!), so
  # store old value befor setting ther global var
  #
  set flag $Equation(updateValuesArea)
  set Equation(updateValuesArea) 0

  Equation(STRESS_ANALYSIS,editProc)
  Equation(STRESS_MODEL,groupEditProc)

  set Equation(updateValuesArea) $flag
}


# NOTE: This is called when Stress Analysis check box is clicked
#
proc Equation(STRESS_ANALYSIS,editProc) {} {
  global Equation Model

  #-No stress analysis selected
  if {$Equation(STRESS_ANALYSIS) != "True"} {
    set Equation(stress_type) ""
    set Equation(stress_mechanicalType) ""
    set Equation(stress_thermalType) ""
    set Equation(STRESS_MODEL) ""

  } else {

    set Equation(stress_type) "又"
  }

  Equation::setUserDefinedFieldsActivity STRESS_ANALYSIS

  #-Broadcast type change to other
  Equation::updateTargetMask
}


proc Equation(STRESS_MODEL,groupEditProc) { {fld ""} } {
  global Equation

  if { $Equation(stress_type) == "" } {
    return
  }

  set Equation(stress_mechanicalType) ""
  set Equation(stress_thermalType) ""

  switch $Equation(STRESS_MODEL) {

    "" {
      if { $Equation(heat_type) != "" } {
        set Equation(stress_mechanicalType) "又T1"
        #set Equation(STRESS_MODEL) "Thermal"
      } else {
        set Equation(stress_thermalType) "又M1"
        #set Equation(STRESS_MODEL) "Mechanical"
      }
    }

    "Mechanical" {
      set Equation(stress_mechanicalType) "又M1"
      set Equation(stress_thermalType) ""
    }

    "Thermal" {
      set Equation(stress_thermalType) "又T1"
      set Equation(stress_mechanicalType) ""
    }
  }

  #-Broadcast type change
  Equation::updateTargetMask
}


#===========================#
#   ADVECTION-DIFFUSION     #
#===========================#


# NOTE: This is called when Advection-Diffusion panel is activated
#
proc Equation(ADVECTION_DIFFUSION_EQUATION,fillProc) {} {
  global Equation EquationVariable

  #-We can not update values area now (loop!), so
  # store old value befor setting ther global var
  #
  set flag $Equation(updateValuesArea)
  set Equation(updateValuesArea) 0

  set Equation(diffusion_type) "" 
  set Equation(diffusion_transferType) "" 

  # No advection-diffusion
  # ======================
  if { $Equation(ADVECTION_DIFFUSION_EQUATION) != "True" } {

    set Equation(updateValuesArea) $flag

    #-Broadcast type change
    Equation::updateTargetMask

    return
  }

  set Equation(diffusion_type) "乙D" 

  Equation(CONVECTION,diffusionCheckProc)  
  Equation(ADVECTION_DIFFUSION_EQUATION,editProc)

  #-Restore global flag var
  set Equation(updateValuesArea) $flag
}


#-Convection advection-diffusion field check proc (done when panel opened
# or in general when data contenst is searched
#
proc Equation(CONVECTION,diffusionCheckProc) { {fld ""} } {
  global Equation

  if { $Equation(diffusion_type) == "" } {
    return
  }

  if { ![info exists Equation(CONVECTION)] } {
    return
  }

  switch $Equation(CONVECTION) {
    "None"     { set Equation(diffusion_transferType) "乙D1"}
    "Constant" { set Equation(diffusion_transferType) "乙D1 乙D2"}
    "Computed" { set Equation(diffusion_transferType) "乙D1 乙D3"}
  }
    
  #-Broadcast type change
  Equation::updateTargetMask
}



# NOTE: This is called when Advection-Diffusion check box is clicked
#
proc Equation(ADVECTION_DIFFUSION_EQUATION,editProc) {} {
  global Equation EquationVariable Info

  # No diffusion selected
  # =====================
  if { $Equation(ADVECTION_DIFFUSION_EQUATION) != "True" } {
    
    # If no heat, disable convection fields
    # =====================================
    if { $Equation(HEAT_EQUATION) != "True" } {
      #set Equation(CONVECTION) "None"
    }

    set Equation(diffusion_type) ""
    set Equation(diffusion_transferType) ""

  # Diffusion selected
  # ==================
  } else {
    set Equation(diffusion_type) "乙D"
  }

  #-Broadcast type change to other
  Equation::updateTargetMask
}



####################
### GRID H panel ###
####################

# Procedure checks if a parameter (pid) is applied to
# vertices. In this case the parameter cannot be
# modified to a nof-elements type density value!
#
proc GridH::parameterHasVertexParents {pid} {
  global GridH ObjectTable

  foreach id $ObjectTable(ids) {

    # Not a vertex
    # We should add 3D edges here also?
    #
    if { "V" != [Object::getType $id] } {
      continue
    }


    # If parameter attached to this boundary, we can stop here!
    #
    if { [info exists ObjectTable($id,gh)] &&
         $ObjectTable($id,gh) == $pid
       } {
      return 1
    }
  }

  return 0
}


# Nof elements-type density parameter cannot be attached to a vertex!
#
proc GridH::parameterTargetCheckProc {paramter_id {accept_var ""} } {
  global GridH Info ObjectTable

  if { $accept_var != "" } {
    upvar #1 $accept_var accept
  }

  set accept 1

  set pid $paramter_id

  set mode [DataField::getFieldValue GridH $pid MESH_DENSITY_TYPE]
  
  if { $mode == "N" &&
       "V" == [Object::getType $GridH(targetId)]
     } {
    set accept 0
  }

  return $accept
}


proc GridH::checkPanel {} {
  global GridH Model ObjectTable

  set id $GridH(targetId)

  if { "V" == [Object::getType $GridH(targetId)] } {
    set GridH(MESH_DT_N,act) 0
  } else {
    set GridH(MESH_DT_N,act) 1
  }

  Widget::configureField GridH MESH_DT_N
}


# Set current mesh index based on selection in the mesh
# option menu
#
proc GridH(MESH_NAME,fillProc) {} {
  global GridH Model

  set GridH(MESH_NAME) [lindex $Model(meshNames) $GridH(MESH_INDEX)]
}


# Set mesh density type field states when parameter is selected
#
proc GridH(MESH_DENSITY_TYPE,fillProc) {} {
  global GridH

  GridH(MESH_DENSITY_TYPE,groupEditProc)
}


# Set mesh density type field states when radiobutton is selected
#
proc GridH(MESH_DENSITY_TYPE,groupEditProc) { {fld ""} } {
  global GridH ObjectTable

  set accept_dt_n 1

  if { "V" == [Object::getType $GridH(targetId)] } {
    set GridH(MESH_DT_N,act) 0
    set accept_dt_n 0
  }

  set GridH(MESH_H,act) 0
  set GridH(MESH_N,act) 0
  set GridH(MESH_R,act) 0

  switch $GridH(MESH_DENSITY_TYPE) {

    "H" {
      set GridH(MESH_H,act) 1
    }

    "N" {
      if { $accept_dt_n } {
        set GridH(MESH_N,act) 1
      }
    }
    
    "R" {
      set GridH(MESH_R,act) 1
    }
  }

  Widget::configureField GridH MESH_H
  Widget::configureField GridH MESH_R
  Widget::configureField GridH MESH_N

  GridH(MESH_DT_N,editProc)
}


# If parameter is attached also to vertices and the mesh density
# type "Number of elements" is selected, user cannot update this
# paramter
#
proc GridH(MESH_DT_N,editProc) {} {
  global GridH Model ObjectTable

  set pid $GridH(parameterId)

  if { [GridH::parameterHasVertexParents $pid] } {
    
    if { "N" == $GridH(MESH_DENSITY_TYPE) } {
      set GridH(updateAllowed) 0

    } else {
      set GridH(updateAllowed) 1
      #set GridH(MESH_N) ""
    }
  }

  Panel::panelDataModified $GridH(dataModified) GridH ""
}



# Update ObjectTable gh-id from ghIds when panel is opend or
# mesh has been changed in the panel
#
proc GridH::updateObjectParameterId {} {
  global GridH Info MeshDefine Model ObjectTable

  if { [info exists MeshDefine(edited,GridH)] &&
       $MeshDefine(edited,GridH)
     } {
    return
  }

  if { ![info exists GridH(MESH_NAME)] } {
    set GridH(MESH_NAME) $Model(currentMeshName)
    set GridH(MESH_INDEX) $Model(currentMeshIndex)
  }

  foreach id $ObjectTable(ids) {
 
    if { [info exists ObjectTable($id,gh)] } {

      # Default value
      set ObjectTable($id,gh) $Info(NO_INDEX)

      # If ids list given, try to find one for the current
      # mesh
      if { $GridH(MESH_INDEX) != $Info(NO_INDEX) &&
           [info exists ObjectTable($id,ghIds)] &&
           $ObjectTable($id,ghIds) != ""
         } { 

        set idx [lsearch $ObjectTable($id,ghMshIndcs) $GridH(MESH_INDEX)]
        # If something for the current mesh!
        if { $idx != -1 } {

          set ObjectTable($id,gh) [lindex $ObjectTable($id,ghIds) $idx]
        }
      }
    }

  }
}


# Update ObjectTable gh-id to ghIds (and ghMshIndcs)
# Used after ok/apply and when mesh has been changed
# in the panel
#
proc GridH::updateObjectParameterIds {} {
  global Info GridH ObjectTable

  if { ![info exist GridH(MESH_INDEX)] ||
       $GridH(MESH_INDEX) == ""        ||
       $GridH(MESH_INDEX) == $Info(NO_INDEX)
     } {
    return
  }

  set msh_idx $GridH(MESH_INDEX)

  foreach id $ObjectTable(ids) {

    if { ![info exists ObjectTable($id,gh)] } {
      continue
    }

    set ids $ObjectTable($id,ghIds)
    set indices $ObjectTable($id,ghMshIndcs)
    set gid $ObjectTable($id,gh)

    if { $gid == $Info(NO_INDEX) } {
      set result [List::removeFromIndexedIdList $indices $ids $msh_idx $gid]
    } else {
      set result [List::addToIndexedIdList $indices $ids $msh_idx $gid]
    }

    set ObjectTable($id,ghMshIndcs) [lindex $result 0]
    set ObjectTable($id,ghIds) [lindex $result 1]

  } ;# Foreach object

  Panel::panelDataChanged 0 GridH
}



############################
### GRID PARAMETER panel ###
############################

# Procedure checks if a parameter (pid) is applied to
# non-structured bodies. In this case the parameter cannot be
# modified to a structured (quadgrid mesh) parameter!
#
proc GridParameter::parameterHasNonStructuredParents {pid} {
  global GridParameter ObjectTable

  foreach id $ObjectTable(ids) {

    # Not a body or boyd-layer
    #
    if { "BL" != [Object::getType $id] &&
         "BL" != [Object::getType $id]
       } {
      continue
    }

    # If structure body, no problems
    #
    if { [info exists ObjectTable($id,accptStrMsh)] &&
         $ObjectTable($id,accptStrMsh)
       } {
      continue
    }

    # If parameter attached to this body, we can stop here!
    #
    if { [info exists ObjectTable($id,gr)] &&
         $ObjectTable($id,gr) == $pid
       } {
      return 1
    }
  }

  return 0
}


# Structured mesh parameter cannot be attached to a non-structurerd body!
#
proc GridParameter::parameterTargetCheckProc {paramter_id {accept_var ""} } {
  global GridParameter Info ObjectTable

  if { $accept_var != "" } {
    upvar #1 $accept_var accept
  }

  set accept 1


  set pid $paramter_id

  set mode [DataField::getFieldValue GridParameter $pid MESH_ELEMENT_TYPE]
  
  if { $mode == "Quad" && !$ObjectTable($GridParameter(objectId),accptStrMsh) } {
    set accept 0
  }

  return $accept
}


proc GridParameter::checkPanel {} {
  global GridParameter Model ObjectTable

  set id $GridParameter(objectId)

  if { $ObjectTable($id,excldMsh) } {
    set GridParameter(ACTIVE_IN_MESHING) 0
  } else {
    set GridParameter(ACTIVE_IN_MESHING) 1
  }

  # For open bodies line is the only option
  #
  if { "OPN" == [Object::getClass $id] ||
       $GridParameter(MESH_ELEMENT_TYPE) == "Line"  
     } {
    set GridParameter(MESH_ET_LINE,act) 1
    set GridParameter(MESH_ET_TRIANGLE,act) 0
    set GridParameter(MESH_ET_QUAD,act) 0
    set GridParameter(MESH_ELEMENT_TYPE) "Line"
    set GridParameter(MESH_LAYER_TYPE) "BoundaryMesh"
    GridParameter(MESH_ET_LINE,editProc)

  } else {
    set GridParameter(MESH_ET_LINE,act) 0
    set GridParameter(MESH_ET_TRIANGLE,act) 1
    set GridParameter(MESH_ET_QUAD,act) 1
  }

  if { !$ObjectTable($id,accptStrMsh) } {
    set GridParameter(MESH_ET_QUAD,act) 0
  } else {
    set GridParameter(MESH_ET_QUAD,act) 1
  }

  Widget::configureField GridParameter MESH_ET_LINE
  Widget::configureField GridParameter MESH_ET_TRIANGLE
  Widget::configureField GridParameter MESH_ET_QUAD
}


proc GridParameter::updateObjectBoundariesWdg {} {
  global GridParameter Model ObjectTable

  set list ""

  foreach id $ObjectTable($GridParameter(objectId),sbIds) {
    
    set row ""
    append row "("
    append row $ObjectTable($id,tg)
    append row ") "
    append row $ObjectTable($id,nm)

    lappend list $row
  }

  ListBox::fill $GridParameter(allWidgets,LISTBOX1) $list
}


# Set state for group of widgets
#
proc GridParameter::setGroupState {group activity} {
  global GridParameter

  foreach fld $group {
    set GridParameter($fld,act) $activity
    Widget::configureField GridParameter $fld
  }
}


proc GridParameter(ACTIVE_IN_MESHING,editProc) {} {
  global GridParameter Model ObjectTable

  if { $Model(nofBodies) == 1 } {
    set GridParameter(ACTIVE_IN_MESHING) 1
  }

  set id $GridParameter(objectId)

  if { !$GridParameter(ACTIVE_IN_MESHING) } {
    set ObjectTable($id,excldMsh) 1
  } else {
    set ObjectTable($id,excldMsh) 0
  }
}


proc GridParameter(LISTBOX1,editProc) {key} {
  global GridParameter ObjectTable

  set ids $ObjectTable($GridParameter(objectId),sbIds)

  set lb $GridParameter(allWidgets,LISTBOX1) 
 
  set sindex [$lb curselection]

  if { $sindex == "" } {
    return
  }

  # Normal selection, highlight corresponding entry!
  #
  if { $key == "ButtonRelease-1" } {
  
    set n1 $GridParameter(MESH_QUADGRID_N1,prev)
    set n2 $GridParameter(MESH_QUADGRID_N2,prev)
    set wdg1 $GridParameter(allWidgets,MESH_QUADGRID_N1)
    set wdg2 $GridParameter(allWidgets,MESH_QUADGRID_N2)

    #wdg1 configure -bg white
    #wdg2 configure -bg white

    if { $sindex == 0 || $sindex == 2 } {
      #wdg1 configure -bg green
      focus $wdg1
      $wdg1 delete 0 end; $wdg1 insert 0 $n1
      $wdg1 icursor end
    } else {
      #wdg2 configure -bg green
      focus $wdg2 
      $wdg2 delete 0 end; $wdg2 insert 0 $n2
      $wdg2 icursor end
    }
    return
  }

  set id [lindex $ids $sindex]

  set prntId $ObjectTable($id,prId)

  set bndr $id
  set bd1 $ObjectTable($prntId,pr1Id)
  set bd2 $ObjectTable($prntId,pr2Id)
  set lr1 -1
  set lr2 -1
  set accept_body_change 1
  set update_gui 1

  #--NOTE: Apply boundary selection via cpp !!!
  #
  set GridParameter(boundaryIds) $ids
  set GridParameter(boundaryLB) $lb
  set GridParameter(boundaryIndex) $sindex

  set data "$bndr $bd1 $lr1 $bd2 $lr2 $accept_body_change $update_gui"

  Util::cpp_exec boundarySelected $data
}


# Set proper widget states when panel data is filled
# (after selecting a paramter etc.)
#
proc GridParameter(MESH_ELEMENT_TYPE,fillProc) {} {
  global GridParameter ObjectTable

  if { $GridParameter(MESH_ELEMENT_TYPE) == "Line" } {
    GridParameter(MESH_ET_LINE,editProc)

  } elseif { $GridParameter(MESH_ELEMENT_TYPE) == "Triangle" } {
    GridParameter(MESH_ET_TRIANGLE,editProc)

  } elseif { $GridParameter(MESH_ELEMENT_TYPE) == "Quad" } {

    if {!$ObjectTable($GridParameter(targetId),accptStrMsh)} {
      GridParameter::setGroupState $GridParameter(quadGridGroup) 0

    } else {      
      set GridParameter(MESH_LAYER_TYPE) "QuadGrid"
      GridParameter(MESH_ET_QUAD,editProc)
    }
  }

  if { ![Util::panelExists GridParameter] } {
    return
  }
}


proc GridParameter(MESH_ET_LINE,editProc) {} {
  global GridParameter Model ObjectTable

  GridParameter::setGroupState $GridParameter(triangleGroup) 0
  GridParameter::setGroupState $GridParameter(seededTriangleGroup) 0

  GridParameter::setGroupState $GridParameter(quadGridGroup) 0

  GridParameter::setGroupState $GridParameter(lineGroup) 1

  set GridParameter(MESH_QUADGRID_N1) ""
  set GridParameter(MESH_QUADGRID_N2) ""

  GridParameter::updateMeshLayerMenu

  GridParameter(MESH_LAYER_TYPE,fillProc)

  GridParameter(MESH_DENSITY_TYPE,fillProc)

  if { ![Util::panelExists GridParameter] } {
    return
  }

}


proc GridParameter(MESH_ET_TRIANGLE,editProc) {} {
  global GridParameter Model ObjectTable

  GridParameter::setGroupState $GridParameter(lineGroup) 0
  GridParameter::setGroupState $GridParameter(quadGridGroup) 0
  GridParameter::setGroupState $GridParameter(triangleGroup) 1

  if {!$ObjectTable($GridParameter(targetId),accptStrMsh)} {
    set GridParameter(MESH_QUADGRID_N1) ""
    set GridParameter(MESH_QUADGRID_N2) ""
  }

  if { [string match "SS*" $GridParameter(MESH_LAYER_TYPE)] } {
    GridParameter::setGroupState $GridParameter(seededTriangleGroup) 1
  } else {
    GridParameter::setGroupState $GridParameter(seededTriangleGroup) 0
  }

  GridParameter::updateMeshLayerMenu

  GridParameter(MESH_LAYER_TYPE,fillProc)

  GridParameter(MESH_DENSITY_TYPE,fillProc)

  if { ![Util::panelExists GridParameter] } {
    return
  }

}


proc GridParameter(MESH_ET_QUAD,editProc) {} {
  global GridParameter Model ObjectTable

  set pid $GridParameter(parameterId)

  GridParameter::setGroupState $GridParameter(lineGroup) 0
  GridParameter::setGroupState $GridParameter(triangleGroup) 0
  GridParameter::setGroupState $GridParameter(seededTriangleGroup) 0
  GridParameter::setGroupState $GridParameter(quadGridGroup) 1

  GridParameter::updateMeshLayerMenu

  GridParameter(MESH_LAYER_TYPE,fillProc)

  GridParameter(MESH_DENSITY_TYPE,fillProc)

  if { ![Util::panelExists GridParameter] } {
    return
  }

}


proc GridParameter::updateMeshLayerMenu {} {
  global GridParameter

  if { ![Util::panelExists GridParameter] } {
    return
  }

  if { $GridParameter(MESH_ELEMENT_TYPE) == "Line" } {
    set idx $GridParameter(meshMethodIndices,LINE)

  } elseif { $GridParameter(MESH_ELEMENT_TYPE) == "Triangle" } {
    set idx $GridParameter(meshMethodIndices,TRIANGLE)

  } elseif { $GridParameter(MESH_ELEMENT_TYPE) == "Quad" } {
    set idx $GridParameter(meshMethodIndices,QUAD_GRID)

  } else {
    return
  }

  set fld MESH_LAYER_TYPE

  # Pick current MenuSelect proc
  #
  set wdg $GridParameter(allWidgets,$fld)
  set ms_proc [ bind $GridParameter(allWidgets,$fld,om) <<MenuSelect>>]

  # Get new menu options
  #
  set tmp [lindex [DataField::getFieldProperty GridParameter $fld Limits] 1]
  set names [lindex $tmp $idx]

  destroy $wdg

  # Build new option-menu
  #
  set om [Panel::createOptionMenuWidget $wdg GridParameter($fld) $names]
  $wdg configure -indicatoron 0 -width 16 
  Widget::configureS $wdg active

  bind $om <<MenuSelect>> $ms_proc

  # Save and pack rebuilt widget
  #
  set GridParameter(allWidgets,$fld) $wdg
  pack $wdg
}


proc GridParameter(MESH_BG_MESH,fillProc) {} {

  GridParameter(MESH_BG_MESH,editProc)
}


proc GridParameter(MESH_BG_MESH,editProc) {} {
  global GridParameter MeshDefine

  # If mesh level bg-mesh is active and it is not
  # a control type-mesh, then body level bg-meshes
  # are not active
  #
  # NOTE: This test not is use currently!!!, MVe 15.02.01
  #
  if { 0 && 
       $MeshDefine(MESH_BG_MESH_ACTIVE) &&
       $MeshDefine(MESH_BG_MESH_FILE) != "" &&
       !$MeshDefine(MESH_BG_MESH_CONTROL)
     } {
    set GridParameter(MESH_BG_MESH_FILE,act) 0

  # Otherwise External bg-mesh type activates
  # the file name entry
  #
  } else {

    if { $GridParameter(MESH_BG_MESH) == "External" } {
      set GridParameter(MESH_BG_MESH_FILE,act) 1
    } else {
      set GridParameter(MESH_BG_MESH_FILE,act) 0
    }
  }

  Widget::configureField GridParameter MESH_BG_MESH_FILE 
}


# Set mesh density type field states when parameter is selected
#
proc GridParameter(MESH_DENSITY_TYPE,fillProc) {} {
  global GridParameter

  GridParameter(MESH_DENSITY_TYPE,groupEditProc)
}


# Set mesh density type field states when radiobutton is selected
#
proc GridParameter(MESH_DENSITY_TYPE,groupEditProc) { {fld ""} } {
  global GridParameter

  set GridParameter(MESH_H,act) 0
  set GridParameter(MESH_R,act) 0

  switch $GridParameter(MESH_DENSITY_TYPE) {

    "H" {
      set GridParameter(MESH_H,act) 1
    }
    
    "R" {
      set GridParameter(MESH_R,act) 1
    }
  }

  if { ![Util::panelExists GridParameter] } {
    return
  }

  Widget::configureField GridParameter MESH_H
  Widget::configureField GridParameter MESH_R
}


# Set mesh seed type states by option menu value when
# parameter is selected
#
proc GridParameter(MESH_LAYER_TYPE,fillProc) {} {
  global GridParameter

  if { [string match "SS*" $GridParameter(MESH_LAYER_TYPE)] } {
    GridParameter::setGroupState $GridParameter(seededTriangleGroup) 1
  } else {
    GridParameter::setGroupState $GridParameter(seededTriangleGroup) 0
  }
}


# Set mesh seed type states when option menu is updated
#
proc GridParameter(MESH_LAYER_TYPE,editProc) {} {
  global GridParameter

  if { [string match "SS*" $GridParameter(MESH_LAYER_TYPE)] } {
    GridParameter::setGroupState $GridParameter(seededTriangleGroup) 1
  } else {
    GridParameter::setGroupState $GridParameter(seededTriangleGroup) 0
  }
}


# Set current mesh name value
#
proc GridParameter(MESH_NAME,fillProc) {} {
  global GridParameter Model

  set GridParameter(MESH_NAME) [lindex $Model(meshNames) $GridParameter(MESH_INDEX)]
}


# Update ObjectTable gr-id from grIds when panel is opend or
# mesh has been changed in the panel
# gr --> grIds
#
proc GridParameter::updateObjectParameterId {} {
  global GridParameter Info MeshDefine Model ObjectTable

  if { [info exists MeshDefine(edited,GridParameter)] &&
       $MeshDefine(edited,GridParameter)
     } {
    return
  }

  if { ![info exists GridParameter(MESH_NAME)] } {
    set GridParameter(MESH_NAME) $Model(currentMeshName)
    set GridParameter(MESH_INDEX) $Model(currentMeshIndex)
  }

  if { $GridParameter(MESH_INDEX) == $Info(NO_INDEX) ||
       $GridParameter(MESH_INDEX) == ""
     } {
    return
  }

  foreach id $ObjectTable(ids) {

    if { ![info exists ObjectTable($id,gr)] } {
      continue
    }

    # Default values
    set ObjectTable($id,gr) $Info(NO_INDEX)
    set ObjectTable($id,excldMsh) 0

    # If grid parameter ids given, try to find one for the current
    # mesh
    set idx [lsearch $ObjectTable($id,grMshIndcs) $GridParameter(MESH_INDEX)]
    if { $idx != -1 } {
      set ObjectTable($id,gr) [lindex $ObjectTable($id,grIds) $idx]
    }

    # If meshing exclude indices given, try to find one for the current
    # mesh
    set idx [lsearch $ObjectTable($id,excldMshIndcs) $GridParameter(MESH_INDEX)]
    if { $idx != -1 } {
      set ObjectTable($id,excldMsh) 1
    } else {
      set ObjectTable($id,excldMsh) 0
    }
  }

  if { [info exists GridParameter(objectLB)] &&
       [winfo exists $GridParameter(objectLB)]
     } {
    StdPanelInit::constructObjectLBList GridParameter
    StdPanelInit::fillObjectListBox GridParameter
    #GridParameter::objectSelected $GridParameter(objectIndex)
    $GridParameter(objectLB) selection set $GridParameter(objectIndex)
  }
}


# Update ObjectTable gr-id to grIds (and grMshIndcs)
# Used after ok/apply and when mesh has been changed
# in the panel
# grIds --> gr
#
proc GridParameter::updateObjectParameterIds {} {
  global GridParameter Info MeshDefine ObjectTable

  if { [info exist GridParameter(MESH_INDEX)] &&
       $GridParameter(MESH_INDEX) != ""       &&
       $GridParameter(MESH_INDEX) != $Info(NO_INDEX)
     } {
    set msh_idx $GridParameter(MESH_INDEX)
  } else {
    set msh_idx [llength MeshDefine(mehNames)]
  }

  foreach id $ObjectTable(ids) {

    if { ![info exists ObjectTable($id,gr)] } {
      continue
    }

    #--Grid parameter ids
    set ids $ObjectTable($id,grIds)
    set indices $ObjectTable($id,grMshIndcs)
    set gid $ObjectTable($id,gr)

    if { $gid == $Info(NO_INDEX) } {
      set result [List::removeFromIndexedIdList $indices $ids $msh_idx $gid]
    } else {
      set result [List::addToIndexedIdList $indices $ids $msh_idx $gid]
    }

    set ObjectTable($id,grMshIndcs) [lindex $result 0]
    set ObjectTable($id,grIds) [lindex $result 1]

    #--Exclude from meshing indices
    set ids $ObjectTable($id,excldMshIndcs)

    if { $ObjectTable($id,excldMsh) } {
      set result [List::addToIdList $ids $msh_idx]

    } else {
      set result [List::removeFromIdList $ids $msh_idx]

    }

    set ObjectTable($id,excldMshIndcs) $result

  } ;# Foreach object

  Panel::panelDataChanged 0 GridParameter
}


# Check that QuadGrid definitions are matching after
# attaching or updating a quadgrid-type parameter
#
proc  GridParameter::checkParameterAttachment {} {
  global GridParameter Info ObjectTable

  #--Current target body
  set tid $GridParameter(targetId)

  #--If no gridable body
  if { $tid == $Info(NO_INDEX) ||
       !$ObjectTable($tid,accptStrMsh) ||
       $ObjectTable($tid,gr) == $Info(NO_INDEX)
     } {
    return
  }

  #--If no quadgrid-type parameter
  if { $GridParameter(MESH_ELEMENT_TYPE) != "Quad" } {
    return
  }

  #--Pick current body element's all sub ids (including inner
  #  boundaries!) and equivalent quadgrid n-values
  #
  set sub_ids $ObjectTable($tid,sbIds)

  set grid_ns ""
  lappend grid_ns $GridParameter(MESH_QUADGRID_N1)
  lappend grid_ns $GridParameter(MESH_QUADGRID_N2)
  lappend grid_ns $GridParameter(MESH_QUADGRID_N1)
  lappend grid_ns $GridParameter(MESH_QUADGRID_N2)

  #--Loop all objects and check if they have a quadgrid-type
  #  mesh parameter for the (current) mesh
  #
  foreach id $ObjectTable(ids) {

    #-Current body
    if { $id == $tid } {
      continue
    }

    #-No quadgrid possible for the body
    if { ![info exists ObjectTable($id,gr)] ||
         !$ObjectTable($id,accptStrMsh) ||
         $ObjectTable($id,gr) ==$Info(NO_INDEX)
       } {
      continue
    }

    set pid $ObjectTable($id,gr)  

    set etype [DataField::getFieldValue GridParameter $pid MESH_ELEMENT_TYPE]

    #-No quadgrid defined for the body
    if { $etype != "Quad" } {
      continue
    }

    #-Ok, check if body under inspection have any of same
    # elements in the current target body
    #
    foreach sid $sub_ids grid_n $grid_ns {

      set idx [lsearch $ObjectTable($id,sbIds) $sid]

      if { $idx == -1 } {
        continue
      }

      #-Overwrite the QuadGrid-N value for the
      # corresponding element

      # Get first proper variable name 
      if { $idx == 0 || $idx == 2 } {
        set fn MESH_QUADGRID_N1
      } else {
        set fn MESH_QUADGRID_N2
      }

      # The update parameter
      DataField::setFieldValue GridParameter $pid $fn $grid_n

    } ;# Foreach subelement in the current target body

  } ;# Foreach object

}


######################
### MATERIAL panel ###
######################

proc Material(checkPanel) {} {
  Material(COMPRESSIBILITY_MODEL,editProc)
}

proc Material(COMPRESSIBILITY_MODEL,fillProc) {} {

  Material(COMPRESSIBILITY_MODEL,editProc)
}


proc Material(COMPRESSIBILITY_MODEL,editProc) {} {
  global Material

  if { [Util::masksAreMatching (人) [Panel::getTargetMask Material] ] } {
    set has_heat_eq 1
  } else {
    set has_heat_eq 0
  }

  set value $Material(COMPRESSIBILITY_MODEL)

  # Note: Possible HeatCapacity field in the flow panel is
  # controlled by incompressibility when the current body
  # does not have any heat-equation
  #
  if { $value == "" || $value == "Incompressible" } {

    set Material(REFERENCE_PRESSURE,act) 0
    set Material(SPECIFIC_HEAT_RATIO,act) 0

    if { !$has_heat_eq } {
      set Material(HEAT_CAPACITY,act) 0
    }
  
  } else {
    set Material(REFERENCE_PRESSURE,act) 1
    set Material(SPECIFIC_HEAT_RATIO,act) 1

    if { !$has_heat_eq } {
      set Material(HEAT_CAPACITY,act) 1
    }
  }

  if { ![Util::panelExists Material] } {
    return
  }

  if { !$has_heat_eq } {
    Widget::configureField Material HEAT_CAPACITY
  }

  if { ![Widget::arrayWidgetExists Material COMPRESSIBILITY_MODEL] } {
    return
  }

  Widget::configureField Material REFERENCE_PRESSURE
  Widget::configureField Material SPECIFIC_HEAT_RATIO

}



####################
### SOLVER panel ###
####################

proc Solver(PROCEDURE,fillProc) {} {
  global Solver

  set vals [split $Solver(PROCEDURE) ";"]

  set Solver(LIBRARY_FILE) [lindex $vals 0]
  set Solver(FUNCTION_NAME) [lindex $vals 1]
}


proc Solver(PROCEDURE,checkProc) {} {
  global Solver

  set lb $Solver(LIBRARY_FILE)
  set fn $Solver(FUNCTION_NAME)

  set Solver(PROCEDURE) [join [list $lb $fn] ";"]

  return ""
}


proc Solver(BUBBLES,editProc) {} {
  global Solver

  if { $Solver(BUBBLES) } {
    set Solver(STABILIZE) "False"

  } else {
    set Solver(STABILIZE) "True"
  }
}


proc Solver(STABILIZE,editProc) {} {
  global Solver

  if { $Solver(STABILIZE) } {
    set Solver(BUBBLES) "False"

  } else {
    set Solver(BUBBLES) "True"
  }
}


# When Solver panel is opened
#
proc Solver(LINEAR_SYSTEM_SOLVER,fillProc) {} {
  global Solver

  Solver(LINEAR_SYSTEM_SOLVER,editProc)
}


# When Solver type OptionMenu is pressed
#
proc Solver(LINEAR_SYSTEM_SOLVER,editProc) {} {
  global Equation Solver

  if { $Solver(LINEAR_SYSTEM_SOLVER) == "Direct" } {
    #Solver::setFieldStates "True" $Solver(direct_group) 
    #Solver::setFieldStates "False" $Solver(iter_group)
    #Solver::setFieldStates "False" $Solver(multigrid_group) 
    set method_fld LINEAR_SYSTEM_DIRECT_METHOD

  } elseif { $Solver(LINEAR_SYSTEM_SOLVER) == "Iterative" } {
    #Solver::setFieldStates "True" $Solver(iter_group) 
    #Solver::setFieldStates "False" $Solver(multigrid_group) 
    #Solver::setFieldStates "False" $Solver(direct_group) 
    Solver::precond_proc
    set method_fld LINEAR_SYSTEM_ITERATIVE_METHOD

  } elseif { $Solver(LINEAR_SYSTEM_SOLVER) == "Multigrid" } {
    #Solver::setFieldStates "False" $Solver(iter_group) 
    #Solver::setFieldStates "True" $Solver(multigrid_group) 
    #Solver::setFieldStates "False" $Solver(direct_group) 
    #Solver::setFieldStates "False" LINEAR_SYSTEM_ITERATIVE_METHOD 
    Solver::precond_proc
    set method_fld LINEAR_SYSTEM_MULTIGRID_METHOD
  }

  set Solver(LINEAR_SYSTEM_METHOD) $Solver($method_fld)
  
  # Rebuild method menu if panel open
  #
  set fld LINEAR_SYSTEM_METHOD

  if { ![info exists Solver(allWidgets,$fld) ] ||
       ![winfo exists $Solver(allWidgets,$fld) ]
     } {
    return
  }

  set wdg $Solver(allWidgets,$fld)

  # Pick current MenuSelect proc
  set ms_proc [ bind $Solver(allWidgets,$fld,om) <<MenuSelect>>]

  destroy $wdg

  set mvals [DataField::getFieldProperty Solver $method_fld Limits]
  set mvals [lindex $mvals 1]

  # Build new method option-menu
  set om [Panel::createOptionMenuWidget $wdg Solver($fld) $mvals]
  $wdg configure -indicatoron 0 -width 12 
  Widget::configureS $wdg active

  bind $om <<MenuSelect>> $ms_proc

  # Save and pack rebuilt widget
  set Solver(allWidgets,$fld) $wdg
  pack $wdg
}


# When Solver method OptionMenu is pressed
#
proc Solver(LINEAR_SYSTEM_METHOD,editProc) {} {
  global Solver

  set fld LINEAR_SYSTEM_METHOD

  if { [Util::nce $Solver(LINEAR_SYSTEM_SOLVER) "Direct"] } {
    set Solver(LINEAR_SYSTEM_DIRECT_METHOD) $Solver($fld)

  } elseif { [Util::nce $Solver(LINEAR_SYSTEM_SOLVER) "Iterative"] } {
    set Solver(LINEAR_SYSTEM_ITERATIVE_METHOD) $Solver($fld)

  } elseif { [Util::nce $Solver(LINEAR_SYSTEM_SOLVER) "Multigrid"] } {
    set Solver(LINEAR_SYSTEM_MULTIGRID_METHOD) $Solver($fld)
  }
}


# When ILU-menu is used
#
proc Solver(LINEAR_SYSTEM_PRECONDITIONING,editProc) {} {
  global Solver

  Solver::precond_proc
}


# Make ILUT-threshold-entry active if preconditiong mode is "ILUT"
#
proc Solver::precond_proc {} {
  global Solver

  if { [Util::nce $Solver(LINEAR_SYSTEM_PRECONDITIONING) "ILUT"] } {
    Solver::setFieldStates "True" LINEAR_SYSTEM_ILUT_TOLERANCE
  } else {
    Solver::setFieldStates "False" LINEAR_SYSTEM_ILUT_TOLERANCE
  }
}


# end ecif_tk_panelCheck.tcl
# ********************


