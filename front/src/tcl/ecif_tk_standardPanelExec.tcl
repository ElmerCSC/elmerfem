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
#Module:    ecif_tk_panelExec.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Procedures for running "standard" panels.
#
#************************************************************************


#################################
# Selection in the OBJECT-box ###
#################################

proc BodyForce::objectSelected { selected_row {select_parameter 1} } {
  StdPanelExec::objectSelected BodyForce $selected_row $select_parameter
}

proc BodyParameter::objectSelected { selected_row {select_parameter 1} } {
  StdPanelExec::objectSelected BodyParameter $selected_row $select_parameter
}

proc Boundary::objectSelected { selected_row {select_parameter 0} } {
  upvar #0 Boundary theArray
  StdPanelExec::objectSelected Boundary $selected_row $select_parameter
}

proc BoundaryCondition::objectSelected { selected_row {select_parameter 1} } {
  global BoundaryCondition
  set BoundaryCondition(RADIATION_TARGET_BODY) ""
  StdPanelExec::objectSelected BoundaryCondition $selected_row $select_parameter
}

proc BoundaryParameter::objectSelected { selected_row {select_parameter 1} } {
  StdPanelExec::objectSelected BoundaryParameter $selected_row $select_parameter
}

proc Equation::objectSelected { selected_row {select_parameter 1} } {
  global Equation
  StdPanelExec::objectSelected Equation $selected_row $select_parameter
}

proc GridParameter::objectSelected { selected_row {select_parameter 1} } {
  global GridParameter ObjectTable
  StdPanelExec::objectSelected GridParameter $selected_row $select_parameter
  GridParameter::updateObjectBoundariesWdg
  GridParameter::checkPanel
  GridParameter::updateMeshLayerMenu
}

proc GridH::objectSelected { selected_row {select_parameter 1} } {
  StdPanelExec::objectSelected GridH $selected_row $select_parameter
}

proc InitialCondition::objectSelected { selected_row {select_parameter 1} } {
  StdPanelExec::objectSelected InitialCondition $selected_row $select_parameter
}

proc Material::objectSelected { selected_row {select_parameter 1} } {
  StdPanelExec::objectSelected Material $selected_row $select_parameter
  Material(checkPanel)
}


#-----------------------
# Common OBJECT SELECTED
#-----------------------

#-Common for all panels which have the elements list-box and have
# a parameters list-box
#
proc StdPanelExec::objectSelected { globArray selected_row {select_parameter 1} { do_selections 1} } {
  global Common Info Model ObjectTable
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  set Common(currentArray) $globArray

  if {0 == [llength [$theArray(objectLB) get 0 end]]} {
    return
  }

  # Reset data
  # ==========
  set theArray(targetMask) ""
  set theArray(objectId) ""
  set theArray(targetId) ""
  set theArray(objectName) ""

  # Set current object info
  # =======================

  #--If selected row was given in the argument
  if { $selected_row != "" &&
       $selected_row != $Info(NO_INDEX)
     } {
    set theArray(objectIndex) [lindex $selected_row 0]
    $theArray(objectLB) selection clear 0 end
    $theArray(objectLB) selection set $theArray(objectIndex) $theArray(objectIndex)
    $theArray(objectLB) see $theArray(objectIndex)

  #--Otherwise pick first selected
  } else {
    set theArray(objectIndex) [lindex [$theArray(objectLB) curselection] 0]
  }

  if { $theArray(objectIndex) == "" } {
    return
  }

  # Set selected rows
  # =================

  #--If selected row given as the argument
  if { $selected_row != "" &&
       $selected_row != $Info(NO_INDEX)
     } {
    set sindices $selected_row

    if { $do_selections } {
      foreach s_index $sindices {
        $theArray(objectLB) selection set $s_index $s_index
      }
    }

  #--Otherwise current selections in the objects list-box and
  # corresponding data-line and element-number key-part
  #
  } else {
    set sindices  [$theArray(objectLB) curselection]
  }


  #--Store and display current body name
  StdPanelExec::setCurrentBodyIdsAndNames $globArray $theArray(objectIndex)

  #--Store id and row-index for the currently selected object
  set theArray(objectId) [lindex $theArray(objectIds) $theArray(objectIndex)]

  # If panel has problem-type concept, set selected target mask
  if { $theArray(hasMaskedParameter) } {

    set trg_mask [StdPanelExec::getTargetMask $globArray $theArray(objectId)]
    set theArray(targetMask) $trg_mask
  }

  #--Set current problem (FEM) variables
  set theArray(currentVariables) [Panel::getCurrentVariables $theArray(targetMask)]

  # Check that selections match the current mask
  # ============================================
  if { [llength $sindices] > 1 } {
    StdPanelExec::checkTargetSelections $globArray \
                                        $theArray(objectLB) \
                                        $theArray(objectIds) \
                                        $theArray(targetMask)

  }


  # No boundaries
  # =============

  #--Reset values area and any param-listbox selection
  if { !$theArray(hasBoundaries) } {

    set theArray(targetId) $theArray(objectId)
    set theArray(targetIndex) $theArray(objectIndex)
  
    if { $theArray(hasParameters) && $select_parameter } {

      #-If we have problem-type dependent parameter data, we have to reset old values
      if { $theArray(hasMaskedParameter) } {

        DataField::clearLBSelection $theArray(parameterLB) $theArray(parameterIndex)
        set theArray(parameterIndex) ""
        set theArray(parameterId) ""
        set theArray(parameterName) ""

        #-Update also current parameter-mask
        if {$theArray(parameterIndex) == ""} {
          set theArray(parameterMask) ""
        } else {
          set theArray(parameterMask) $o_mask
        }

        #-Mark all parameters which match current object's mask
        #ListBox::markParameterLBByMask $globArray $o_mask 
        ListBox::markParameterLBByMask $globArray  
      }

      #--Make parameter selected (if any attached)
      #  If we work only at object-level, find and select the possible parameter definition
      #  for the object (in the case there is not yet a selected parameter!)
      set pv [string tolower $theArray(parameterType)]
      set trg_id $theArray(targetId)

      if { $ObjectTable($trg_id,$pv) == $Info(NO_INDEX) } {

        Widget::configureS $theArray(panelDetachButton) disabled
        
        if { !$theArray(hasMaskedFields) } {
          StdPanelExec::parameterSelected $globArray 0


        } else {
          #-Values area states
          StdPanelExec::clearValuesArea $globArray
          StdPanelInit::initFieldValues $globArray
          PanelCheck::execPanelFillProcs $globArray
          StdPanelExec::setValuesAreaActivity $globArray $theArray(targetMask)
          set theArray(updateAllowed) 0
        }

      } else  {
        Widget::configureS $theArray(panelDetachButton) normal
        StdPanelExec::findAndActivateParameter $globArray $trg_id
        set theArray(updateAllowed) 1
      }

    } ;# End if select parameter


    StdPanelExec::updateEquationMenu $globArray
    return

  } ;# End if no boundaries



  # Has boundaries
  # ==============

  StdPanelExec::updateBoundaryLBData $globArray

  set boundaries [List::getListRowsByIndex $theArray(boundaryLBList) \
                                           $theArray(boundaryIndices)]

  ListBox::fill $theArray(boundaryLB) $boundaries

  if { [llength $selected_row] < 2 } {
    set bndr_row 0

  } else {
    set bndr_id  [lindex $selected_row 1]
    set bndr_row [lsearch $theArray(boundaryIds) $bndr_id]

    if { $bndr_row == -1 } {
      set bndr_row 0
    }
  }

  StdPanelExec::boundarySelected $globArray $bndr_row $select_parameter
}


###########################################
### Double selection in the OBJECTS-box ###
###########################################

proc BodyForce::objectDblSelected {} {
  StdPanelExec::objectDblSelected BodyForce
}

proc BodyParameter::objectDblSelected {} {
  StdPanelExec::objectDblSelected BodyParameter
}

proc Boundary::objectDblSelected {} {
  StdPanelExec::objectDblSelected Boundary
}

proc BoundaryCondition::objectDblSelected {} {
  StdPanelExec::objectDblSelected BoundaryCondition
}

proc BoundaryParameter::objectDblSelected {} {
  StdPanelExec::objectDblSelected BoundaryParameter
}

proc Equation::objectDblSelected {} {
  StdPanelExec::objectDblSelected Equation
}

proc GridParameter::objectDblSelected {} {
  StdPanelExec::objectDblSelected GridParameter
}

proc GridH::objectDblSelected {} {
  StdPanelExec::objectDblSelected GridH
}

proc InitialCondition::objectDblSelected {} {
  StdPanelExec::objectDblSelected InitialCondition
}


#--A double-selected object is marked (toggle) and
#  then possible drawn via C++
# Original version
proc InitialCondition::objectDblSelected2 {} {
  global InitialCondition Common Info

  set source $InitialCondition(objectLB)

  if { 0 == [llength [$source get 0 end]] } {
    return
  }

  #-Current selection in the objects list-box and
  # corresponding data-line and body-number key-part
  set s_index [$source curselection]
  set hasMark [ListBox::toggleSelectionMark $source $s_index]
  set data [lindex $InitialCondition(objectList) $s_index]
  set bd_id [lindex [split $data $Info(fieldSeparator)] 1]
  set lr_id -1

  #--C++ is called
  set data "$bd_id $lr_id"
  Util::cpp_exec bodySelected $data
}


proc Material::objectDblSelected {} {

  StdPanelExec::objectDblSelected Material
}


#---------------------------
# Common OBJECT DBL-SELECTED
#---------------------------

#--A double-selected object
proc StdPanelExec::objectDblSelected {globArray} {
  global Common Model ObjectTable
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  set Common(currentArray) $globArray

  set bd_id $theArray(body1Id)
  set lr_id -1

  #--C++ is called
  set data "$bd_id $lr_id"
  Util::cpp_exec bodySelected $data
}


# Pick all boundary ids which belong to the object "objectId" and whose type
# belong to "types"
#
proc StdPanelExec::findCurrentBoundaryIds {globArray types objectId} {
  global ObjectTable
  upvar #0 $globArray theArray

  set result_ids ""

  set parentIds $objectId

  # Collect all ids hierachically
  # ============================
  while { $parentIds != "" } {

    set newParentIds ""

    # Loop all parents
    # ================
    foreach pid $parentIds {
       
      # Loop all parent's sub elements
      # ==============================
      foreach sid $ObjectTable($pid,sbIds) {
        
        if { -1 == [lsearch $theArray(allBoundaryIds) $sid] } {
          continue
        }

        # Add if not yet selected
        if { -1 == [lsearch $result_ids $sid] } {
          lappend result_ids $sid
          lappend newParentIds $sid
        }
      } ;# all sub elements

    } ;# all parents

    # Sub elements are new parents
    # ===========================
    set parentIds $newParentIds

  } ;#while parents

  set result_ids [Object::sortIdsByKeyVars \
                          $result_ids \
                          { {tp Object::codeBoundaryTypeForSorting} {tg ""} } \
                          { 0 1 } ]

  return $result_ids
}


# Create boundary listbox data when objectId is selected
#
proc StdPanelExec::updateBoundaryLBData {globArray} {
  upvar #0 $globArray theArray

  set theArray(boundaryIds) [StdPanelExec::findCurrentBoundaryIds \
                               $globArray \
                               $theArray(subObjectTypes) \
                               $theArray(objectId)]
  StdPanelInit::constructBoundaryLBList $globArray
  set theArray(targetIds) $theArray(boundaryIds)

  set theArray(boundaryIndices) ""

  set index 0
  foreach id $theArray(boundaryIds) {
    lappend theArray(boundaryIndices) $index
    incr index
  }
}


proc StdPanelExec::checkTargetSelections {globArray trg_wdg trg_ids mask} {
  global ObjectTable
  upvar #0 $globArray theArray

  if { $mask == "" } {
    return
  }

  set sindices [$trg_wdg curselection]

  if { 1 >= [llength $sindices] } {
    return
  }

  foreach idx $sindices {

    set id [lindex $trg_ids $idx]

    #set m $ObjectTable($id,msk)
    set m [StdPanelExec::getTargetMask $globArray $id]

    if { ![Util::masksAreEqual $mask $m] } {
      $trg_wdg selection clear $idx $idx
    }
  }
}



#####################################
### Selection in the BOUNDARY-box ###
#####################################

proc BoundaryParameter::boundarySelected { selected_row {select_parameter 0} {do_selections 1} } {

  StdPanelExec::boundarySelected BoundaryParameter $selected_row $select_parameter $do_selections
}

proc Boundary::boundarySelected { selected_row {select_parameter 0} {do_selections 1} } {

  StdPanelExec::boundarySelected Boundary $selected_row $select_parameter $do_selections

  Boundary::writeDiscretationData
}

proc BoundaryCondition::boundarySelected { selected_row {select_parameter 1} {do_selections 1} } {

  StdPanelExec::boundarySelected BoundaryCondition $selected_row $select_parameter $do_selections
}


# NOTE: GridParamter does not actually have a boundary box, but for 
# Quadri meshes there is a list of four boundaries!!!
#
proc GridParameter::boundarySelected { selected_row {select_parameter 1} {do_selections 1} } {

  GridParameter::checkPanel
}

proc GridH::boundarySelected { selected_row {select_parameter 1} {do_selections 1} } {

  StdPanelExec::boundarySelected GridH $selected_row $select_parameter $do_selections
  GridH::checkPanel
}


#-------------------------
# Common BOUNDARY SELECTED
#-------------------------

# Handle selections in the boundary listbox
# NOTE: boundaryIndex must be either set in the proc "StdPanelExec::boundarySelectedPre" or
# given explicitely in the argument "selected_row"
#
proc StdPanelExec::boundarySelected {globArray selected_row {select_parameter 1} {do_selections 1} } {
  global Common Info Model ObjectTable
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  set Common(currentArray) $globArray

  if { 0 == [llength [$theArray(boundaryLB) get 0 end] ] } {
    return
  }

  # Reset data
  # ==========

  set theArray(targetMask) ""
  set theArray(boundaryId) ""
  set theArray(boundaryName) ""
  
  # Set current boundary info
  # =========================

  #--If selected row given as the argument
  if { $selected_row != "" && $selected_row != $Info(NO_INDEX) } {
    set theArray(boundaryIndex) [lindex $selected_row 0]

  #--Otherwise pick first selected
  } else {
    set theArray(boundaryIndex) [lindex [$theArray(boundaryLB) curselection] 0]
  }

  if { $theArray(boundaryIndex) == "" } {
    return
  }

  set theArray(boundaryId) [lindex $theArray(boundaryIds) $theArray(boundaryIndex)]

  set theArray(targetId) $theArray(boundaryId)
  set theArray(targetIndex) $theArray(boundaryIndex)

  # If panel has problem-type concept, set selected target mask
  if { $theArray(hasMaskedTarget) } {

    set trg_mask [StdPanelExec::getTargetMask $globArray $theArray(boundaryId)]
    set theArray(targetMask) $trg_mask

    ListBox::markParameterLBByMask $globArray  
  }

  # Set selected rows
  # =================

  #--If selected row given as the argument
  if { $selected_row != "" &&
       $selected_row != $Info(NO_INDEX)
     } {
    set sindices $selected_row

    if { $do_selections } {
      foreach s_index $sindices {
        $theArray(boundaryLB) selection set $s_index $s_index
      }
    }

  #--Otherwise current selection in the elements list-box and
  # corresponding data-line and element-number key-part
  #
  } else {
    set sindices  [$theArray(boundaryLB) curselection]
  }

  set index1 [lindex $sindices 0]

  if { $index1 != "" && $index1 != $Info(NO_INDEX) } {
    $theArray(boundaryLB) see $index1
  }

  # Multiple boundaries were selected
  # =================================

  # Check that selections match the current mask
  # ============================================
  if { [llength $sindices] > 1 } {
    StdPanelExec::checkTargetSelections $globArray \
                                        $theArray(boundaryLB) \
                                        $theArray(boundaryIds) \
                                        $theArray(targetMask)

    return
  }


  # Only one boundary was selected
  # ==============================
  set trg_id $theArray(boundaryId)
  set trg_mask $theArray(targetMask)

  set theArray(boundaryName) $ObjectTable($trg_id,nm)

  #--Reset values area and any param-listbox selection
  if { $select_parameter  && $theArray(hasParameters) } {

    if { $theArray(hasMaskedParameter) } {

      DataField::clearLBSelection $theArray(parameterLB) $theArray(parameterIndex)
      set theArray(parameterIndex) ""
      set theArray(parameterId) ""
      set theArray(parameterName) ""

      #-Update also current parameter-mask
      if {$theArray(parameterIndex) == ""} {
        set theArray(parameterMask) ""
      } else {
        set theArray(parameterMask) $trg_mask
      }

      #-Mark all parameters which match current object's mask
      ListBox::markParameterLBByMask $globArray  
    }
  }

  #--Reset boundary name entry status
  # NOTE: Normally this is not relevant, because
  # boundary name is not editable 
  #
  if { $globArray == "Boundary" } {
    set theArray(boundaryName,err) 0
    set theArray(boundaryName,mod) 0
    set theArray(boundaryName,prev) $theArray(boundaryName)
    Widget::setEntryStatus $globArray boundaryName $theArray(boundaryNameWidget)
  }

  StdPanelExec::updateEquationMenu $globArray

  set theArray(updateAllowed) 1

  #--Make parameter selected (if any attached)
  if { $select_parameter && $theArray(hasParameters) } {

    set pv [string tolower $theArray(parameterType)]

    if { $ObjectTable($trg_id,$pv) == $Info(NO_INDEX) } {

      Widget::configureS $theArray(panelDetachButton) disabled
      #-Values area states
      StdPanelExec::clearValuesArea $globArray
      StdPanelInit::initFieldValues $globArray
      PanelCheck::execPanelFillProcs $globArray
      StdPanelExec::setValuesAreaActivity $globArray $trg_mask
      set theArray(updateAllowed) 0

    } else  {
      Widget::configureS $theArray(panelDetachButton) normal
      StdPanelExec::findAndActivateParameter $globArray $trg_id
    }

  #--Otherwise mark data unmodified
  } else {
    Panel::panelDataModified 0 $globArray
  }

}



############################################
### Double selection in the BOUNDARY-box ###
############################################

#--A double-selected element is marked (toggle) and
#  then drawn via C++

proc Boundary::boundaryDblSelected { {accept_body_change 0} } {

  StdPanelExec::boundaryDblSelected Boundary $accept_body_change
}

proc BoundaryParameter::boundaryDblSelected { {accept_body_change 0} } {

  StdPanelExec::boundaryDblSelected BoundaryParameter $accept_body_change
}

proc BoundaryCondition::boundaryDblSelected { {accept_body_change 0} } {

  StdPanelExec::boundaryDblSelected BoundaryCondition $accept_body_change
}

proc GridH::boundaryDblSelected { {accept_body_change 0} } {

  StdPanelExec::boundaryDblSelected GridH $accept_body_change
}


#-----------------------------
# Common BOUNDARY DBL-SELECTED
#-----------------------------

#--A double-selected element is marked (toggle) and
#  then drawn via C++
#
proc StdPanelExec::boundaryDblSelected { globArray {accept_body_change 0} } {
  global Common Info Model ObjectTable
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  set Common(currentArray) $globArray

  #--All current listbox selctions
  set s_indices [ $theArray(boundaryLB) curselection]

  set selected_row ""

  #--Find out what was the last (the new) selection in the listbox

  # Check if something was added to selections, ie listbox row
  # was selected, but the boundary is not currently marked as selected
  #
  set counter 0
  foreach id $theArray(boundaryIds) {

    set bndr_selected $ObjectTable($id,slctd)
    set src_selected [lsearch $s_indices $counter]

    if { !$bndr_selected && $src_selected != -1
       } {
      set selected_row $counter
      break
    }

    incr counter
  }

  #--NOTE: Boundary selection is done via cpp !!!
  set bndr $theArray(boundaryId)
  set bd1 $theArray(body1Id)
  set bd2 $theArray(body2Id)
  set lr1 -1
  set lr2 -1
  set update_gui 1


  set data "$bndr $bd1 $lr1 $bd2 $lr2 $accept_body_change $update_gui"

  Util::cpp_exec boundarySelected $data
}



####################################################
### Control Double selection in the BOUNDARY-box ###
####################################################

proc Boundary::boundaryCntrlDblSelected { {accept_body_change 0} } {
  StdPanelExec::boundaryCntrlDblSelected Boundary $accept_body_change
}

proc BoundaryParameter::boundaryCntrlDblSelected { {accept_body_change 0} } {
  StdPanelExec::boundaryCntrlDblSelected BoundaryParameter $accept_body_change
}

proc BoundaryCondition::boundaryCntrlDblSelected { {accept_body_change 0} } {
  StdPanelExec::boundaryCntrlDblSelected BoundaryCondition $accept_body_change
}

proc GridH::boundaryCntrlDblSelected { {accept_body_change 0} } {
  StdPanelExec::boundaryCntrlDblSelected GridH $accept_body_change
}


#---------------------------------
# Common OBJECT CNTRL-DBL-SELECTED
#---------------------------------

#--A control+double-selected element is marked (toggle) and
#  then drawn via C++
proc StdPanelExec::boundaryCntrlDblSelected { globArray {accept_body_change} } {

  # Set flag for extended selection
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_EXTEND 1

  # Call then normal double selection proc
  StdPanelExec::boundaryDblSelected $globArray $accept_body_change

  # Reset flags
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_EXTEND 0
}



######################################
### Selection in the PARAMETER-box ###
######################################

proc BodyForce::parameterSelected { selected_row } {
  StdPanelExec::parameterSelected BodyForce $selected_row
}

proc BodyParameter::parameterSelected { selected_row } {
  StdPanelExec::parameterSelected BodyParameter $selected_row
}

proc BoundaryCondition::parameterSelected { selected_row } {
  global BoundaryCondition Info Model

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  StdPanelExec::parameterSelected BoundaryCondition $selected_row

  set pid $BoundaryCondition(parameterId)

  if { $pid == $Info(NO_INDEX) } {
    return
  }

  set oid $BoundaryCondition($pid,oid)
  set mask [StdPanelExec::getTargetMask BoundaryCondition $oid]

  if { $mask != $BoundaryCondition(targetMask) } {
    #Panel::initField BoundaryCondition RADIATION_TARGET_BODY_INDEX 1
    #BoundaryCondition(RADIATION_TARGET_BODY_INDEX,fillProc)
  }

  # Check if widget must be disabled to prevent changing
  # the target body!
  #
  #BoundaryCondition(RADIATION_TARGET_BODY_NAME,checkWidgetState)
  BoundaryCondition(RADIATION_TARGET_BODY,fillProc)
}

proc BoundaryParameter::parameterSelected { selected_row } {
  StdPanelExec::parameterSelected BoundaryParameter $selected_row
}

proc Calculator::parameterSelected { selected_row } {
  StdPanelExec::parameterSelected Calculator $selected_row
}

proc Equation::parameterSelected { selected_row } {
  global Equation Model

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  Equation::initMasks

  #--Run general proc first
  StdPanelExec::parameterSelected Equation $selected_row

  set pid $Equation(parameterId)

  if {$pid == "" } {
    return
  }

  set Equation(targetMask) $Equation($pid,mask)

  #Equation::updateEquationList $Equation(equationIndex)
  #Equation::setVariableWidgetsActivity $Equation(equationIndex)
  #Equation::updateVariableList
}

proc GridParameter::parameterSelected { selected_row } {

  StdPanelExec::parameterSelected GridParameter $selected_row

  GridParameter::checkPanel
}

proc GridH::parameterSelected { selected_row } {

  StdPanelExec::parameterSelected GridH $selected_row
  GridH::checkPanel
}

proc InitialCondition::parameterSelected { selected_row } {

  StdPanelExec::parameterSelected InitialCondition $selected_row
}

proc Material::parameterSelected { selected_row } {

  StdPanelExec::parameterSelected Material $selected_row
  Material(checkPanel)
}


#--------------------------
# Common PARAMETER SELECTED
#--------------------------

# Common for all panels
# Data from (one) selected row is written to data fields.
# select_target arguments tells if the target (body or boundary) for the
# parameter should be shown selected
#
proc StdPanelExec::parameterSelected {globArray {selected_row ""} { select_target 0} } {
  global Common Info Model ObjectTable
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  set Common(currentArray) $globArray

  if { 0 == [llength [$theArray(parameterLB) get 0 end]] } {
    PanelCheck::execPanelFillProcs $globArray
    return
  }

  #--If selected row given as the argument
  if { $selected_row != "" &&
       $selected_row != $Info(NO_INDEX)
     } {
    set theArray(parameterIndex) $selected_row
    $theArray(parameterLB) selection set $selected_row $selected_row

  #-Otherwise row-index for the currently selected parameter
  } else {
    set theArray(parameterIndex) [$theArray(parameterLB) curselection]
  }

  #-If nothing was selected (or selection was cleared elsewhere!)
  if {$theArray(parameterIndex) == "" } {
    PanelCheck::execPanelFillProcs $globArray
    return
  }

  set theArray(updateAllowed) 1

  set pid [lindex $theArray(ids) $theArray(parameterIndex)]

  if { $select_target } {
    set find_first 0
    StdPanelExec::findAndActivateTarget $globArray $pid $find_first
    #return
  }

  #-Check if parameter mask matches target mask  
  set matching_masks 1

  if { $theArray(hasMaskedParameter) } {

    set tmask $theArray(targetMask)
    set pmask $theArray($pid,mask)

    if { $pmask != $tmask } {
      set matching_masks 0
    }

    set theArray(parameterMask) $tmask
    set theArray(targetMask) $tmask
  }

  #-Check if some other restriction defined against the current
  # target
  set accept_target 1

  if { [info exists theArray(hasParameterTargetCheckProc)] &&
       $theArray(hasParameterTargetCheckProc)
     } {
    set pn $globArray
    append pn ::parameterTargetCheckProc
    set accept_target [eval $pn $pid]
  }
  #-Set current parameter id and name
  set theArray(parameterId) $pid
  set theArray(parameterName) [StdPanelExec::getParameterName $globArray $pid]

  #-Store current parameter name
  StdPanelExec::storeCurrentParameterName $globArray

  #-Mark panel unmodified
  Panel::panelDataModified 0 $globArray
  set theArray(parameterName,prev) $theArray(parameterName)
  set theArray(parameterName,err) 0
  set theArray(parameterName,mod) 0

  if { [info exists theArray(parameterNameWidget)] } {
    Widget::setEntryStatus $globArray parameterName $theArray(parameterNameWidget)
  }

  #-Reset values area
  StdPanelExec::clearValuesArea $globArray

  #--Store old values
  Panel::backupFields $globArray

  Panel::resetFields $globArray

  #--Update data-field variables (in the globArray)
  StdPanelInit::initFieldValues $globArray

  # Mask matching
  # =============
  if { $accept_target && $matching_masks } {
    DataField::formDataFields $globArray $theArray($pid,data) 
    set attach_state normal    

  } else {
    DataField::formDataFields $globArray $theArray($pid,data) $tmask
    set attach_state disabled    
  }

  if { [info exists theArray((panelAttachButton)] } {
    Widget::configureS $theArray(panelAttachButton) $attach_state
  }

  #-Values area states
  PanelCheck::execPanelFillProcs $globArray
  StdPanelExec::setValuesAreaActivity $globArray $theArray(targetMask)

  # Id matching
  # ===========
  #-Check is parameter-id matches that of the current target
  set matching_pids 0
  set pv [string tolower $theArray(parameterType)]

  if { [info exists theArray(targetId)] &&
       $theArray(targetId) != $Info(NO_INDEX) &&
       $theArray(targetId) != "" &&
       $ObjectTable($theArray(targetId),$pv) == $pid
     } {
    set matching_pids 1
  }

  #-Set update button state
  if { $matching_pids } {
    #Widget::configureS $theArray(panelUpdateButton) normal

  } else {
    Widget::configureS $theArray(panelUpdateButton) disabled
  }
}



#########################################
### DblSelection in the PARAMETER-box ###
#########################################

proc BodyForce::parameterDblSelected {} {
  StdPanelExec::parameterSelected BodyForce  "" 1
}

proc BodyParameter::parameterDblSelected {} {
  StdPanelExec::parameterSelected BodyParameter  "" 1
}

proc BoundaryCondition::parameterDblSelected {} {
  global BoundaryCondition

  $BoundaryCondition(boundaryLB) selection clear 0 end
  StdPanelExec::parameterSelected BoundaryCondition  "" 1
  BoundaryCondition::boundaryDblSelected 

  #StdPanelExec::selectBoundaries BoundaryCondition $BoundaryCondition(parameterId)
}

proc BoundaryParameter::parameterDblSelected {} {
  StdPanelExec::parameterSelected BoundaryParameter  "" 1
}


proc Calculator::parameterDblSelected { selected_row } {
  StdPanelExec::parameterSelected Calculator "" 1
}

proc Equation::parameterDblSelected {} {
  global Equation Model

  if { $Model(GEOMETRY_DIMENSION) == "" } {
    return
  }

  Equation::initMasks

  #--Run general proc first
  StdPanelExec::parameterSelected Equation  "" 1

  set pid $Equation(parameterId)

  if {$pid == "" } {
    return
  }

  set Equation(targetMask) $Equation($pid,mask)
}

proc GridH::parameterDblSelected {} {
  global GridH

  $GridH(boundaryLB) selection clear 0 end
  StdPanelExec::parameterSelected GridH  "" 1
  GridH::boundaryDblSelected 

}

proc GridParameter::parameterDblSelected {} {
  StdPanelExec::parameterSelected GridParameter  "" 1
}

proc InitialCondition::parameterDblSelected {} {
  StdPanelExec::parameterSelected InitialCondition "" 1
}

proc Material::parameterDblSelected {} {
  StdPanelExec::parameterSelected Material  "" 1
}



##############################################
### CntrlDblSelection in the PARAMETER-box ###
##############################################

proc BoundaryCondition::parameterCntrlDblSelected {} {
  global BoundaryCondition

  $BoundaryCondition(boundaryLB) selection clear 0 end

  StdPanelExec::parameterSelected BoundaryCondition  "" 1

  set accept_body_change 1
  BoundaryCondition::boundaryCntrlDblSelected $accept_body_change

}



#####################
### Restore procs ###
#####################

proc BoundaryCondition::restoreProc {} {

  Util::restoreFlagValue BoundaryCondition DRAW_TARGET ""
  Util::restoreFlagValue BoundaryCondition DRAW_SOURCE DRAW_SOURCE_CAD
  Util::restoreFlagValue BoundaryCondition DRAW_SOURCE DRAW_SOURCE_MESH
}



##################
### Misc procs ###
##################

# Set model status variables after boundary condition
# update
#
proc BoundaryCondition::updateModelStatus { check_factors } {
  global BoundaryCondition Info Model ObjectTable

  # Reset flags.
  # ===========
  set Model(hasDiffuseGrayRadiation) 0

  # Do no rest if checking is turned off!
  if { $Model(Viewfactors,needsUpdate,data) != -1 } { 
    set Model(Viewfactors,needsUpdate,data) 0
  }

  # Do no rest if checking is turned off!
  if { $Model(GebhardtFactors,needsUpdate,data) != -1} {
    set Model(GebhardtFactors,needsUpdate,data) 0
  }

  # Check status
  # ============

  set arr BoundaryCondition
  set fld_em EMISSIVITY
  set fld_rd RADIATION

  set has_diffuse_gray 0              ;# Any diffuse gray --> viewfactors needed
  set gebhardtfactors_needs_update 0  ;# Flag for relevant change
  set viewfactors_needs_update 0      ;# Flag for relevant change

  foreach id $BoundaryCondition(allBoundaryIds) {

    set new_bc $Info(NO_INDEX)
    set old_bc $Info(NO_INDEX)

    if { [info exist ObjectTable($id,bc)] } {
      set new_bc $ObjectTable($id,bc)
    }

    if { [info exist ObjectTable($id,bc,old)] } {
      set old_bc $ObjectTable($id,bc,old)
    }

    #-If no bc at all
    if { $new_bc == $Info(NO_INDEX) &&
         $old_bc == $Info(NO_INDEX)
       } {
      continue
    }

    set bndr_has_diffuse_gray 0

    #-Pick new values: Emissivity, DiffuseGray
    #
    set new_val_em ""
    set new_val_rd ""

    if { $new_bc != $Info(NO_INDEX) } {
      set new_val_em [DataField::getFieldValue $arr $new_bc $fld_em]
      set new_val_rd [DataField::getFieldValue $arr $new_bc $fld_rd]
    }

    #-Pick old values: Emissivity, DiffuseGray
    # NOTE: We read old field value using also old bc-id value
    # because it may have changed also!
    #
    set old_val_em ""
    set old_val_rd ""
    set get_old_value 1

    if { $old_bc != $Info(NO_INDEX) &&
         [info exists BoundaryCondition($old_bc,data,old)]
       } {
      set old_val_em [DataField::getFieldValue $arr $old_bc $fld_em "" $get_old_value]
      set old_val_rd [DataField::getFieldValue $arr $old_bc $fld_rd "" $get_old_value]
    }

    # Model has DiffuseGray
    # =====================
    if { $new_val_rd != "" && [string equal -nocase $new_val_rd "Diffuse gray"] } {
      set bndr_has_diffuse_gray 1
      set has_diffuse_gray 1
    }

    # If do not check factors (ex. a new model)
    # =========================================
    if { !$check_factors } {
      continue
    }

    # DiffuseGray has changed
    # ==========================
    if { $bndr_has_diffuse_gray && $new_val_rd != $old_val_rd } {
      set viewfactors_needs_update 1
    }
   
    # DiffGr Emissivity has changed
    # =============================
    if { $bndr_has_diffuse_gray && $new_val_em != $old_val_em
       } {
      set gebhardtfactors_needs_update 1
    }

  }

  # ===================
  # Set new flag values
  # ===================

  # Radiation related flags
  # =======================
  if {$has_diffuse_gray} {
    set Model(hasDiffuseGrayRadiation) 1
  }
 
  # Update flag if checking is not turned off
  if { $viewfactors_needs_update &&
       $Model(Viewfactors,needsUpdate,data) != -1 
     } {
    set Model(Viewfactors,needsUpdate,data) 1
  }

  # Update flag if checking is not turned off
  if { $gebhardtfactors_needs_update  &&
       $Model(GebhardtFactors,needsUpdate,data) != -1 
     } {
    set Model(GebhardtFactors,needsUpdate,data) 1
  }

  if { $has_diffuse_gray } {
    set state normal
  } else {
    set state disabled
  }

  MenuBuild::configureMenuButtonOption Run Viewfactors state $state
  MenuBuild::configureMenuButtonOption Run GebhardtFactors state $state
}

 
################
### OK procs ###
################


proc BoundaryCondition::okProcPre { close_window} {

  set check_factors 1
  BoundaryCondition::updateModelStatus $check_factors

  if { $close_window } {
    BoundaryCondition::restoreProc
  }
}


proc Calculator::okProcPre {} {
  global Calculator

  set def_name [DataField::getFieldProperty Calculator VARIABLE InitialValue]
  
  set get_old 1
  set vars_added 0

  foreach pid $Calculator(ids) {

    set when_to_exec [DataField::getFieldValue Calculator $pid EXEC_SOLVER]

    set var_name [DataField::getFieldValue Calculator $pid VARIABLE]
    set var_name_old [DataField::getFieldValue Calculator $pid VARIABLE "" $get_old]

    set dofs [DataField::getFieldValue Calculator $pid VARIABLE_DOFS]
    set dofs_old [DataField::getFieldValue Calculator $pid VARIABLE_DOFS "" $get_old]

    # Remove previous names (if different from default name (like 'Dummy')
    #
    if { ![Util::nce $def_name $var_name_old] } {
      DataField::unsetVariableComponents $var_name_old $dofs_old 0
    }

    # Add current names (if different from default name (like 'Dummy')
    #
    if { ![Util::nce $def_name $var_name] && ![Util::nce "Never" $when_to_exec] } {
      DataField::setVariableComponents "" $var_name $dofs
      set vars_added 1
    }

  }

  if { $vars_added } {
    Panel::initFields Variables
  }
}

proc Calculator::okProcPost {} {

  SolverOrder::updateSolversAndCalculators
}


# Equation panel ok proc
#
proc Equation::okProcPre {} {
  global Common Equation EquationVariable Info ObjectTable Model

  set is_changed 0

  #--Check if parameter definitions must be updated
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }
   
    if { ![info exists ObjectTable($id,eq)] }  {
      continue
    }

    #-Old and new body masks
    set old_mask [Object::getMask $id]
    set new_mask ""

    set eid $ObjectTable($id,eq)

    if { $eid != $Info(NO_INDEX) } {
      set new_mask $Equation($eid,mask)
    }

    # If a non-empty equation was changed
    if { $old_mask != "" && $old_mask != $new_mask } {
      set is_changed 1
    }
  }

  #--Check if user wants to continue after changes
  set result ""
  if { $is_changed } {
    set msg [list "NOTE: An existing Body equation(s) was changed.\n\n" \
                  "Check parameter definitions in the Model menu! \n\n" \
                  $Info(anywayOk) ]

    if { ![Message::verifyAction $Info(noviceUser) $msg ] } {
      return 0
    }

  }

  #--Check if all bodies have an equation definition
  foreach id $ObjectTable(ids) {
    
    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    # If equation definition for the body is missing
    if {$ObjectTable($id,eq) == $Info(NO_INDEX)} {

      set Model(status) 1

      set msg [list "NOTE: A body must have an equation definition\n"  \
                    "before you can enter any other data for the body!\n\n" \
                    $Info(anywayOk) ]

      # If user wants to correct something
      if { ![Message::verifyAction $Info(noviceUser) $msg] } {
        return 0

      # Inform only on the first missing eq def.
      } else {
        break
      }
    }
  }

  # Set status ok
  set Model(status) 0

  #--Update problems and systems related data
  #
  Equation::updateEquationIndices
  Equation::constructProblemMask

  StdPanelInit::createObjectTableMasks
	
	#--Check Equation data
	# NOTE: Do not save checked equation parameters via this call!
	#
	set modified1 0
	set save_modified 0 
	Panel::checkParameters Equation modified1 $save_modified

  #--Check all mask related arries
	# NOTE: Here we should save all (other) modified paramters!
	#
  set modified2 0
	set save_modified 1 
  Panel::checkMaskedParameters modified2 $save_modified

  #--Mark panel sizes unknown (possible new fields!)
  foreach arr $Common(allPanelArries) {
    if { $arr == $Common(panelArray,EQ) } {
        continue
      }
    upvar #0 $arr theArr
    set theArr(panelSizeKnown) 0
  }

  #--Update fields in other panels if variables
  #  definition were changed
  if { $EquationVariable(dataWasChanged) } {
    
    #-Update variable fields
    Panel::initFields Variables

    #-Update Variables array
    DataField::setVariables
    
    #-Update arries for possible indexed fields
    foreach arr $Common(standardPanelArries) {
      if { $arr == $Common(panelArray,EQ) } {
        continue
      }
      if { $arr == $Common(panelArray,SL) } {
        continue
      }
      DataField::setAllFields $arr
    }
   
  }

  return 1
}


proc Equation::okProcPost {} {
  global Equation EquationVariable

  set modified 0
  Solver::checkSolverData $Equation(activeIndices) modified

  set EquationVariable(dataWasChanged) 0
}


proc GridH::okProcPre {} {
  global MeshDefine

  set MeshDefine(edited,GridH) 1
}


proc GridParameter::okProcPre {} {
  global MeshDefine ObjectTable

  set MeshDefine(edited,GridParameter) 1
}


####################
### CANCEL procs ###
####################


proc BoundaryCondition::cancelProc {} {

  BoundaryCondition::restoreProc
}


proc Equation::cancelProc {} {
  global Equation EquationVariable
  global Variables InitialCondition BoundaryCondition

  if { !$EquationVariable(dataWasChanged) } {
    return
  }

  Util::unsetArrayIdVariables EquationVariable

  # This should cancel also a possible ok in the 
  # EquationVariable-panel!
  # NOTE: We restore the data stored in special Equation variable 
  # to EquationVariable array!!
  # We do not actually have multiple EqVar parameters, but just for the
  # case...
  foreach id $Equation(equationVariables,ids) {
    set EquationVariable($id,data) $Equation(equationVariables,$id,data)
  }
  
  set Variables(initialFields) $Equation(variables,initialFields)
  set Variables(allFields) $Equation(variables,allFields)

  set InitialCondition(initialFields) $Equation(initialCondition,initialFields)
  set InitialCondition(allFields) $Equation(initialCondition,allFields)

  set BoundaryCondition(initialFields) $Equation(boundaryCondition,initialFields)
  set BoundaryCondition(allFields) $Equation(boundaryCondition,allFields)

  set BoundaryCondition(initialDirichletVnames) $Equation(boundaryCondition,initialDirichletVnames)
  set BoundaryCondition(dirichletVnames) $Equation(boundaryCondition,dirichletVnames)

  Panel::resetFields Variables

  Panel::initFields Variables

  # NOTE: Id 1 is always and the only active!!!!
  set data $EquationVariable(1,data)
  DataField::formDataFields EquationVariable $data

  # Save EquationVariable panel to restore the data in cpp-side
  # and mark it really unchanged!
  # NOTE EquationPanel save will set wasChanged to true, so the order
  # is important here!
  set inform_front 1
  set update_equation_panel 0
  EquationVariable::panelSave $inform_front $update_equation_panel
  set EquationVariable(dataWasChanged) 0

  set modified 0
  Solver::checkSolverData $Equation(activeIndices) modified

  Panel::panelDataChanged 0 Equation
}

 
#################
### ADD procs ###
#################

# Equation name is used as the paramter name, so
# it must be updted before the existence of the
# parameter (name) is checked!
#
proc Calculator::addProcPre {} {
  global Calculator

  #set Calculator(parameterName) [string trim $Calculator(EQUATION)]
}


# Equation pre processing before add
proc Equation::addProcPre {} {
  global Equation

  # Form equation target-mask
  Equation::updateVariableTypeTargetMasks 
  Equation::updateTargetMask

  set has_equations ""

  # If this is not a problen name variable
  foreach fld $Equation(allFields) {

    if { "" == [DataField::getFieldProperty Equation $fld EquationIndex "" 0] } {
      continue
    }

    if { $Equation($fld) == "True" } {
      set has_equations "ok"
      break;
    }
  }

  return $has_equations
}


# Equation post processing
#
# Nothing currently!
proc Equation::addProcPost {} {
  global Equation
 
}


####################
### UPDATE procs ###
####################


proc Calculator::updateProcPre {} {
  global Calculator

  #set Calculator(parameterName) $Calculator(EQUATION)
}


# Equation pre processing before update
proc Equation::updateProcPre {} {
  global Equation

  set has_equations ""

  foreach fld $Equation(allFields) {

    # If this is not a problen name variable
    if { "" == [DataField::getFieldProperty Equation $fld EquationIndex "" 0] } {
      continue
    }

    if { $Equation($fld) == "True" } {
      set has_equations "ok"
      break;
    }
  }

  return $has_equations
}


# Procedure updates the list which stores problem-types by bodies when
# an equation is updated in the equation-panel
#
proc Equation::updateProcPost {row_index} {
  global Equation ObjectTable

  set pid $Equation(parameterId)
  
  # Form equation target-mask
  #
  Equation::updateVariableTypeTargetMasks 
  Equation::updateTargetMask

  set Equation($pid,mask) $Equation(targetMask)

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    if { [info exists ObjectTable($id,eq)] &&
         $ObjectTable($id,eq) == $pid
       } {
      set ObjectTable($id,msk) $Equation(targetMask)
    }
  }
}


# Add radiation target body info
# NOTE: This is used when we can have multiple target bodies in the
# result list
# This needs also RADIATION_BOUNDARIES variables which storee the 
# corresponditn boundary tags
#
proc  BoundaryCondition::updateRadiationTargetInfo {target_ids {upd_parameter 1} } {
  global BoundaryCondition Info ObjectTable

  if { ![info exists BoundaryCondition(RADIATION)] ||
       $BoundaryCondition(RADIATION) == "" ||
       ![string equal -nocase $BoundaryCondition(RADIATION) "Diffuse gray"]
     } {
    return
  }

  set pid $BoundaryCondition(parameterId)

  foreach sid $target_ids {

    set mode 1
    set id $sid

    # Detach!
    if { $sid < 0 } {
      set mode -1
      set id [expr -1 * $id]
    }

    set bndr_tag $ObjectTable($id,tg)
    set bd_tag $BoundaryCondition(RADIATION_TARGET_BODY)

    set list1 $BoundaryCondition(RADIATION_BOUNDARIES)
    set list2 $BoundaryCondition(RADIATION_TARGET_BODIES)


    set idx [lsearch $list1 $bndr_tag]

    #-Attach
    if { $mode == 1 } {
      set list1 [lreplace $list1 $idx $idx $bndr_tag]
      set list2 [lreplace $list2 $idx $idx $bd_tag]

    #-Detach
    } else {
      set list1 [lreplace $list1 $idx $idx]
      set list2 [lreplace $list2 $idx $idx]
    }

    # Sort list by boundary tags
    set sort_list ""
    foreach l1 $list1 l2 $list2 {
      lappend sort_list [list $l1 $l2]
    }

    set sort_list [lsort -index 0 -integer $sort_list]

    set list1 ""
    set list2 ""

    foreach row $sort_list {

      lappend list1 [lindex $row 0]
      lappend list2 [lindex $row 1]
    }

    set BoundaryCondition(RADIATION_BOUNDARIES) $list1
    set BoundaryCondition(RADIATION_TARGET_BODIES) $list2
 
    if { $upd_parameter } {

      set an BoundaryCondition
      set fn1 RADIATION_BOUNDARIES
      set fn2 RADIATION_TARGET_BODIES

      set len [llength $list1]
      set desc "($len)"

      DataField::setFieldValue $an $pid $fn1 $list1 $desc
      DataField::setFieldValue $an $pid $fn2 $list2 $desc
    }
  }
}


####################
### DELETE procs ###
####################

# Procedure should update all the model data which is somehow
# connected to the problem-type deleted
#
proc Equation::deleteProc {sel_index} {
  global Equation

  StdPanelExec::deleteParameterById Equation $Equation(parameterId)
}


####################
### ATTACH procs ###
####################

# Procedure updates the list which stores problem-types by bodies when
# an equation is applied to a body in the equation-panel
#
proc Equation::attachProc {obj_sel_ids} {
  global Equation ObjectTable

  set pid $Equation(parameterId)

  set Equation(targetMask) $Equation($pid,mask)

  foreach oid $obj_sel_ids {
    set ObjectTable($oid,msk) $Equation(targetMask)
  }
} 


####################
### DETACH procs ###
#################### 

proc Equation::detachProc {elm_indices box_indices} {

}


###############################
### PANEL button procedures ###
###############################

#-Get globArray which corresponds the active window
#
proc StdPanelExec::getPanelArray { } {
  global Common
  global InitialCondition BoundaryCondition 
  global Equation BodyForce Material
  global GridParameter GridH

  if { $Common(currentStandardPanel) != "" } {
    #return $Common(currentStandardPanel)
  }

  #--Find toplevel-window for the panel

  #-Widget in focus
  set w [focus -displayof . ]

  #-Widget's toplevel (ie. panel)
  set w [winfo toplevel $w]

  #-Id for the panel window
  set id [winfo atom $w]

  #--Find array-name corresponding the panel
  if {$id == $InitialCondition(winId) } {
    return InitialCondition

  } elseif {$id == $BoundaryCondition(winId) } {
    return BoundaryCondition

  } elseif {$id == $Equation(winId) } {
    return Equation

  } elseif {$id == $BodyForce(winId) } {
    return BodyForce

  } elseif {$id == $Material(winId) } {
    return Material

  } elseif {$id == $GridParameter(winId) } {
    return GridParameter

  } elseif {$id == $GridH(winId) } {
    return GridH
  }

  #--If nothing was found
  return ""
}


############# 
# OK-button #
#############

#-For data panels
#-Data will be saved to C++-variables!!!
#-Result data consists of element data and parameter data.
# This data is tranferred as a combined list into
# cpp-environment, where it is split and handled separately
# Also the ids of the possibly deleted parameters are sent to cpp.
#
proc StdPanelExec::panelOk { {globArray ""} {inform_front 1} {close_panel 1} } {
  global Common Info Model

  if { $globArray == "" } {
    set globArray [StdPanelExec::getPanelArray ]
  }

  upvar #0 $globArray theArray

  if { $globArray == $Common(panelArray,EQ) } {
    $theArray(entryInfoLabel) configure -text "Checking model data..." -fg red
    Util::doWait
  }

  #--Turn auto add on, if none of the selected objecs has a paramter
  #  attached
  set all_must_be 1
  set has_unattached [StdPanelExec::selectedTargetsAreUnattached $globArray $all_must_be]

  if { $theArray(dataModified) && $has_unattached } {
    set theArray(autoAdd) 1
  } else {
    set theArray(autoAdd) 0
  }

  if { ![Panel::verifySave $globArray] } {
    return
  }

  # Close also all "hanging" entry panels
  if { $close_panel } {
    Panel::closeEntryPanels $globArray
    set theArray(procedureEntryPanelIds) ""
    set theArray(tableEntryPanelIds) ""
  }

  Panel::panelDataModified 0 $globArray

  #---If data was not updated 
  if {$theArray(dataChanged) == 0} {

    if { $close_panel } {
      destroy $theArray(winName)
    }

    if { $globArray == $Common(panelArray,BC) && $close_panel } {
      BoundaryCondition::restoreProc
    }

    return
  }

  StdPanelExec::compressParameterIds $globArray

  #---Panel specific presave processing
  switch $globArray {

    BoundaryCondition {
      #-reset DiffuseGray flags etc.
      BoundaryCondition::okProcPre $close_panel
    }

    Calculator {
      Calculator::okProcPre
    }

    Equation {
      #-check if equations are all well defined
      if { ![Equation::okProcPre] } {
        return 0
      }
    }

    GridH {

      GridH::okProcPre
    }

    GridParameter {

      GridParameter::okProcPre
    }

    default {
    }
  }

  Panel::panelDataChanged 0 $globArray

  #---Call save proc  
  if { ![StdPanelExec::panelSave $globArray $inform_front] } {
    return
  }

  #--Panel specific postsave processing
  switch $globArray {

    Calculator {
      Calculator::okProcPost
    }

    Equation {
      Equation::okProcPost
    }

  }

  if { $globArray == "Equation" } {
    $theArray(entryInfoLabel) configure -text "" -fg gray
  }

  if { $close_panel } {
    Panel::cancel $theArray(winName)
    #Panel::unsetFields $globArray
    Panel::unsetAllFieldData $globArray
    }
}


################
# APPLY-Button #
################

proc StdPanelExec::panelApply {globArray} {

  set inform_front 1
  set close_panel 0

  StdPanelExec::panelOk $globArray $inform_front $close_panel
}


#################
# CANCEL-Button #
#################

#-A "reset" for data panels
#
proc StdPanelExec::panelCancel {w {globArray ""} } {

  if { $globArray == "" } {
    set globArray [StdPanelExec::getPanelArray]
  }

  upvar #0 $globArray theArray 

  if { ![Panel::verifyCancel $globArray] } {
    return
  }

  # Close also all "hanging" entry panels
  Panel::closeEntryPanels $globArray
  set theArray(procedureEntryPanelIds) ""
  set theArray(tableEntryPanelIds) ""

  StdPanelExec::resetDataLists $globArray
  StdPanelExec::resetUpdateLists $globArray
  set theArray(nextNewParameterId) $theArray(nextNewParameterId,old)

  Panel::cancel $w

  #---Panel specific processing
  switch $globArray {
    Equation {
      #Panel::unsetFields $globArray
      Panel::unsetAllFieldData $globArray
      Equation::cancelProc
      return
    }
    BoundaryCondition {
      BoundaryCondition::cancelProc
    }
  }

  #Panel::unsetFields $globArray
  Panel::unsetAllFieldData $globArray
}


################################
# NEXT-button and PREV-buttons #
################################


# Procedure moves to next or previous Panel equation
# direction: +-1
#
proc StdPanelExec::panelEqMove {globArray direction {update 1} } {
  upvar #0 $globArray theArray

  set old_idx [lsearch $theArray(currentEquationIndices) $theArray(equationIndex)]

	set new_idx [expr $old_idx + $direction]

	# If acceptable new value, set new index and update panel
	if { $new_idx >= 0 && $new_idx <= [llength $theArray(currentEquationIndices)] } {
		set theArray(equationIndex) [lindex $theArray(currentEquationIndices) $new_idx]

    # NOTE: This should normally be called!!!.
    # (Solver panel is an exception, see Solver::panelEqMove)
    #
    if { $update } {
      StdPanelCreate::constructPanelValuesArea $globArray $theArray(equationIndex)
    }
	}
	
  StdPanelExec::setEqMoveButtonsState $globArray

  # Tell if equation was actually changed
  #
  if { $old_idx == $new_idx } {
    return 0
  } else {
    return 1
  }
}


proc StdPanelExec::setEqMoveButtonsState { globArray} {
  upvar #0 $globArray theArray
  
  set last [expr [llength $theArray(currentEquationIndices)] - 1]
  set idx [lsearch $theArray(currentEquationIndices) $theArray(equationIndex)]

	# Set button states
	if { $idx == 0 } {
		Widget::configureButton $theArray(prevEqButton) disabled
	} else {
		Widget::configureButton $theArray(prevEqButton) normal
	}

	if { $idx == $last } {
		Widget::configureButton $theArray(nextEqButton) disabled
	} else {
		Widget::configureButton $theArray(nextEqButton) normal
	}
}


# Procedure moves to next or previous Panel page
# direction: +-1
#
proc StdPanelExec::panelPageMove {globArray direction} {
  upvar #0 $globArray theArray

	set new_val [expr $theArray(panelPageNbr) + $direction]

	# If acceptable new value, update panel
	if { $new_val >= 1 && $new_val <= $theArray(maxNofPanelPages) } {
		set theArray(panelPageNbr) $new_val

		# We force the construction of the panel although there is no
		# equation change!!!
	  set force_construct 1
		StdPanelCreate::constructPanelValuesArea $globArray $theArray(equationIndex) $force_construct
	}
	
	# Set button states
	if { $theArray(panelPageNbr) == 1 } {
		Widget::configureButton $theArray(prevPageButton) disabled
	} else {
		Widget::configureButton $theArray(prevPageButton) normal
	}

	if { $theArray(panelPageNbr) == $theArray(maxNofPanelPages) } {
		Widget::configureButton $theArray(nextPageButton) disabled
	} else {
		Widget::configureButton $theArray(nextPageButton) normal
	}

  # Tell if page was actually moved
  #
  if { $new_val == $theArray(panelPageNbr) } {
    return 1
  } else {
    return 0
  }
}



###############
# KEEP-button #
###############

# Procedure takes entry data into "keep", so that the user
# can select an other parameter, copy some data and come to the
# original situation with Back-button
#
proc StdPanelExec::keepParameter {globArray} {
  global Info ObjectTable
  upvar #0 $globArray theArray

	#--Form the parameter line from panel's active fields
	set theArray(keep,objectIndex) $theArray(objectIndex)
	set theArray(keep,boundaryIndex) $theArray(boundaryIndex)
	set theArray(keep,parameterIndex) $theArray(parameterIndex)

  set theArray(keep,dataDirty) $theArray(dataDirty)
  set theArray(keep,dataModified) $theArray(dataModified)
  set theArray(keep,dataChanged) $theArray(dataChanged)

	set theArray(keep,name) $theArray(parameterName)
  set theArray(keep,data) [DataField::formParameterLine $globArray]

  Widget::configureS $theArray(panelBackButton) normal

  #---Panel specific processing
  switch $globArray {
  }
}


###############
# BACK-button #
###############

# Come back to the situation what existed when the Keep-button was
# pressed
#
proc StdPanelExec::backParameter {globArray} {
  global Info ObjectTable
  upvar #0 $globArray theArray

	if { $theArray(keep,objectIndex) >= 0 } {
		$theArray(objectLB) selection clear 0 end
		StdPanelExec::objectSelected $globArray $theArray(keep,objectIndex)
	}

	if { $theArray(keep,boundaryIndex) >= 0 } {
		$theArray(boundaryLB) selection clear 0 end
		StdPanelExec::boundarySelected $globArray $theArray(keep,boundaryIndex)
	}

	if { $theArray(keep,parameterIndex) >= 0 } {
		$theArray(parameterLB) selection clear 0 end
		StdPanelExec::parameterSelected $globArray $theArray(keep,parameterIndex)
	}

  #-Reset values area and field values
  StdPanelExec::clearValuesArea $globArray
  Panel::resetFields $globArray

	#--Fill panel fields with the data which was on the panel
	#  when Keep-button was pressed
	set theArray(parameterName) $theArray(keep,name)
  DataField::formDataFields $globArray $theArray(keep,data) 

  #-Values area states
  PanelCheck::execPanelFillProcs $globArray
  StdPanelExec::setValuesAreaActivity $globArray $theArray(targetMask)

 #--Data status
  set theArray(dataDirty) $theArray(keep,dataDirty)
  set theArray(dataModified) $theArray(keep,dataModified)
  set theArray(dataChanged) $theArray(keep,dataChanged)

	Panel::panelDataModified $theArray(dataModified) $globArray

  # After update data is clean, no need to update!
  #Widget::configureS $theArray(panelUpdateButton) disabled

  Widget::configureS $theArray(panelBackButton) disabled
	
  #---Panel specific processing
  switch $globArray {
  }
}


##############
# ADD-Button #
##############

# Procedure construct a parameter entry from the data field values
# into parameter list-box (and related data strutures).
#
proc StdPanelExec::addParameter {globArray} {
  global Info
  upvar #0  $globArray theArray 
  
  #--If the panel has problem-type concept, an object must
  #  be selected before adding a new parameter
  if { $theArray(hasMaskedParameter) } {

    if {$theArray(targetIndex) == "" } {

      set Info(messageIcon) warning
      Message::dispOkMessage "Select first an object!" \
                    "$Info(FRONT_NAME) message" \
                    $theArray(winName)
      return 0
    } 
  }

  #--Panel specific pre processing
  switch $globArray {

    Calculator {
      Calculator::addProcPre
    }
    
    Equation {
      if { "ok" != [Equation::addProcPre] } {
        set Info(messageIcon) warning
        Message::dispOkMessage "No equations selected!" 
        return 0
      }
    }

    default {
    }
  }

  #---Check name
  set theArray(parameterName) [string trim $theArray(parameterName)]

  if { ![Panel::checkParameterName $globArray $Info(NO_INDEX) $theArray(parameterName) 1] } {
    return 0
  }

  #--Check that all data values are valid.
  if { 0 == [PanelCheck::checkFields $globArray] } {
    return 0
  }

  #--Run other field check procs
  set theArray(mode) "ADD"
  PanelCheck::execPanelFillProcs $globArray
  set theArray(mode) ""

  #--Check that new parameter line is ok (not empty etc.)
  if { ![StdPanelExec::checkNewParameterLine $globArray ADD] } {
    return 0
  }

  #--Store old values
  Panel::backupFields $globArray

  #--Form the parameter line from panel's active fields
  set new_data [DataField::formParameterLine $globArray]

  # Set field modified states (to Ok state in fact)
  StdPanelExec::setValuesAreaStatus $globArray

  #--Reset parameter name entry status
  set theArray(parameterName,prev) $theArray(parameterName)
  set theArray(parameterName,err) 0
  set theArray(parameterName,mod) 0
  
  if { [info exists theArray(parameterNameWidget)] } {
    Widget::setEntryStatus $globArray parameterName $theArray(parameterNameWidget)
  }

  #--New parameter id
  set theArray(parameterId) $theArray(nextNewParameterId)
  incr theArray(nextNewParameterId)

  #--Add new parameter
  set pid $theArray(parameterId)

  lappend theArray(ids) $pid
  set theArray($pid,data) $new_data
  set theArray($pid,name) $theArray(parameterName)

  if { [info exists theArray(targetId)] } {
    set theArray($pid,oid) $theArray(targetId)
  }

  #--Update list-box stuff
  set new_line [StdPanelInit::formParameterLBRow $globArray $pid]
  lappend theArray(parameterLBList) $new_line
  $theArray(parameterLB) insert end $new_line
  $theArray(parameterLB) selection clear 0 end

  #--Update parameter's equation-type list
  # NOTE: In equation panel we copy this parameter mask to the target body,
  # it is not used to control Equation panel state!
  #
  if { $theArray(hasMaskedTarget) } {
    set theArray($pid,mask) $theArray(targetMask)
  } else {
    set theArray($pid,mask) ""
  }

  # NOTE: The following makes the added row selected, but it 
  #"really" doesn't work:  selection is displayed, but it is not "curselection"!!!
  $theArray(parameterLB) selection set end end
  set theArray(parameterIndex) [expr [$theArray(parameterLB) size] - 1]

  #--Mark panel updated
  set theArray(dataModified) 0

  #--Mark panel updateable
  set theArray(updateAllowed) 1

  #--Mark panel touched and changed (we added!)
  set theArray(dataChanged) 1
  set theArray(dataDirty) 1

  #--Attach added parameter
  # Note This proc will also update data* flags!
  #
  set trg_ids [StdPanelExec::attachParameter $globArray 0]

  #--Panel specific post processing
  switch $globArray {

    Equation {
      Equation::addProcPost
    }

    GridParameter {
      GridParameter::checkParameterAttachment
    }

    default {
    }
  }

  Widget::configureS $theArray(panelDeleteButton) normal
  Widget::configureS $theArray(panelUpdateButton) normal
  Widget::configure1 $theArray(panelAddButton) -fg $Info(en_nmod_color)

  return 1
}


#################
# UPDATE-button #
#################

# Procedure updates a parameter list-box entry with the new values entered
# in the data fields. 
#
proc StdPanelExec::updateParameter {globArray} {
  global Common Info ObjectTable
  upvar #0 $globArray theArray 

  #--Panel specific pre processing
  switch $globArray {

    BoundaryCondition {
      #BoundaryCondition::updateRadiationTargetInfo $theArray(boundaryId) 0
    }

    Equation {
      if { "ok" != [Equation::updateProcPre] } {

        set Info(messageIcon) warning
        Message::dispOkMessage "No equations selected for update!" 
        return 0
      }
    }

    default {
    }
  }

  #--Parameters box must be selected (but only one row at a time!!!)
  set sindex [$theArray(parameterLB) curselection]

  if {$sindex == "" } {

    set Info(messageIcon) warning
    Message::dispOkMessage "Nothing was selected for update!" \
                "$Info(FRONT_NAME) message" \
                $theArray(winName)
    return 0
  }
  
  #--Check that data values are valid.
  if { 0 == [PanelCheck::checkFields $globArray] } {
    StdPanelExec::setValuesAreaStatus $globArray
    return 0
  }

  #-Run other field check procs
  PanelCheck::execPanelFillProcs $globArray


  #--Parameter id
  set pid $theArray(parameterId)

  #---Check name
  set theArray(parameterName) [string trim $theArray(parameterName)]

  if { ![Panel::checkParameterName $globArray $pid $theArray(parameterName) 0] } {
    return 0
  }

  #--Store current param's name
  StdPanelExec::storeCurrentParameterName $globArray
  
  set theArray(parameterName,prev) $theArray(parameterName)
  set theArray(parameterName,err) 0
  set theArray(parameterName,mod) 0

  if { [info exists theArray(parameterNameWidget)] } {
    Widget::setEntryStatus $globArray parameterName $theArray(parameterNameWidget)
  }

  #--Check that new parameter line is ok (not empty etc.)
  if { ![StdPanelExec::checkNewParameterLine $globArray UPDATE] } {
    return 0
  }

  #--Store old values
  Panel::backupFields $globArray

  #--Form the parameter data line
  set ndata [DataField::formParameterLine $globArray]

  # Set field modified states (to Ok state in fact)
  StdPanelExec::setValuesAreaStatus $globArray

  #--Update data
  set theArray($pid,data) $ndata

  #--Form updated lb-row

  set show_mark 0

  # Note: Parameter matching (the *-marker) is concluded
  # from the Attach-button state here!
  #
  #-Add marker
  if { [info exists theArray(panelAttachButton)] &&
       "normal" == [$theArray(panelAttachButton) cget -state]
     } {
    set show_mark 1
  }

  set nrow [StdPanelInit::formParameterLBRow $globArray $pid $show_mark]

  set theArray(parameterLBList) [lreplace $theArray(parameterLBList) $sindex $sindex $nrow]
  ListBox::updateRow $theArray(parameterLB) $sindex $nrow

  #--Update target list

  # Name for parameter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  # Find target indices for the parameter
  set counter 0
  set box_indices ""

  foreach id $theArray(allTargetIds) {

    if { $ObjectTable($id,$pv) == $pid } {
      lappend box_indices $counter
    }
    incr counter
  }

  #--Mark data updated
  set theArray(dataDirty) 0
  set theArray(dataModified) 0

  #--Mark data changed
  set theArray(dataChanged) 1

  #--"Attach" updated parameter (name was possibly changed!)
  # Note This proc will also update data* flags!
  #
  set trg_ids [StdPanelExec::attachParamId $globArray "" $box_indices $pid]

  #--Panel specific post processing
  switch $globArray {

    BoundaryCondition {
      #BoundaryCondition::updateRadiationTargetInfo $theArray(boundaryId)
    }
 
    Equation {
      Equation::updateProcPost $sindex
    }

    GridParameter {
      GridParameter::checkParameterAttachment
    }

  }

  # After update data is clean, no need to update!
  Widget::configureS $theArray(panelUpdateButton) disabled

  return 1
}
 

#################
# DELETE-button #
#################

#-Procedure deletes selected parameter row.
# All possible applications to "elements" are detached.
#
proc StdPanelExec::deleteParameter {globArray } {
  global Common Info
  upvar #0 $globArray theArray 

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  #--Selected parameter row
  set cur_sel_index [$theArray(parameterLB) curselection]

  if {$cur_sel_index == "" } {

    set Info(messageIcon) warning
    Message::dispOkMessage "Nothing was selected to be deleted!" \
                "$Info(FRONT_NAME) message" \
                $theArray(winName)
    return
  }


  #--Pick current parameter's mask, so that we can make a
  #  similar parameter selected after delete
  #
  set pid $theArray(parameterId)
  set current_mask $theArray($pid,mask)

  #--Find new selection
  # NOTE: Do this before delete, so we still have
  # pid in the theArray(ids) lists!
  set npid [StdPanelExec::getNextParameterInstance $globArray $current_mask $pid]

  #--Delete selected parameter row from parameter lists
  # NOTE: Equation panel is handled separately
  # owing to the gloabl effect of the equation definition!

  if { $globArray == $Common(panelArray,EQ) } {
    Equation::deleteProc $cur_sel_index 

  } else {
    StdPanelExec::deleteParameterById $globArray $pid
  }
  
  if { $theArray(parameterId) == [expr $theArray(nextNewParameterId) - 1] } {
    incr theArray(nextNewParameterId) -1
  }

  if { $npid != $Info(NO_INDEX) } {
    set sindex [lsearch $theArray(ids) $npid]
  } else {
    set sindex 0
  }

  set theArray(parameterIndex) $sindex

  set theArray(parameterIsDeleted) 1

  # If no parameters left
  if { 0 == [llength $theArray(ids)] } {

    Widget::configureS $theArray(panelDeleteButton) disabled
    Widget::configureS $theArray(panelUpdateButton) disabled
    set theArray(parameterIndex) ""

  } else {
    StdPanelExec::parameterSelected $globArray $theArray(parameterIndex)
  }

  #--Mark data changed
  set theArray(dataChanged) 1
}


#################
# ATTACH-button #
#################

# Procedure applies the selected parameter to the selected bodies/elments 
# Return value: target object ids
#
proc StdPanelExec::attachParameter {globArray {overwrite 1} } {
  global Info ObjectTable
  upvar #0  $globArray theArray 

  if { !$theArray(hasTarget) } {
    return ""
  }

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]
  
  #--If no parameter is active, do nothing
  if {$theArray(parameterIndex) == ""} {

    set Info(messageIcon) warning
    Message::dispOkMessage "Select first a $theArray(defaultName) to attach!" \
                  "$Info(FRONT_NAME) message" \
                  $theArray(winName)
    return ""
  }

  set object_index [$theArray(objectLB) curselection]

  set trg_ids ""
  set box_indices ""

  #--Store the selection indices
  if {$theArray(hasBoundaries)} {
    set box_indices [$theArray(boundaryLB) curselection]

  } else {
    set box_indices [$theArray(objectLB) curselection]
  }

  #--If nothing was selected!
  if { ($object_index == "") && ($box_indices == "") } {

    set Info(messageIcon) warning
    Message::dispOkMessage "Select first an object or a element!" \
                  "$Info(FRONT_NAME) message" \
                  $theArray(winName)
    return ""
  }

  #---Either objectsLB or elementsLB was selected (but not both!!)
  if { $box_indices == "" } {

    #-Object was selected, but no elements --> means all elements which
    # don't have a parameter defined!
    if {$theArray(hasBoundaries)} {

      set counter 0
      foreach eid $theArray(targetIds) {

        if { $Info(NO_INDEX) == $ObjectTable($eid,$pv) } {
          lappend box_indices $counter
          lappend trg_ids $eid
        }
        incr counter
      }
    }

  # Some targets were selected explicitely
  } else {
    
    foreach idx $box_indices {
      lappend trg_ids [lindex $theArray(targetIds) $idx]
    }
  }

  #--Update parameter-id value for each selected bodies/elements.
  set pid $theArray(parameterId)

  # Check that parameter has a parent object id
  #
  if { $theArray($pid,oid) == $Info(NO_INDEX) } {
    set theArray($pid,oid) [lindex $trg_ids 0]
    set modified 0
    set new_data [DataField::checkParameterData $globArray $pid modified $theArray(targetMask)]
    set new_mask $theArray(targetMask)
    set theArray($pid,data) $new_data
    set theArray($pid,mask) $new_mask
  }

  set att_trg_ids [StdPanelExec::attachParamId $globArray $trg_ids $box_indices $pid $overwrite]

  # Mark data changed
  set theArray(dataChanged) 1

  # Attach was done, panel is "untouched"
  set theArray(dataDirty) 0

  Widget::configureS $theArray(panelDetachButton) normal

  #--Panel specific processing
  switch $globArray {
    
    BoundaryCondition {
      #BoundaryCondition::updateRadiationTargetInfo $att_trg_ids
    }

    Equation {
      Equation::attachProc $att_trg_ids 
    }

    GridParameter {
      GridParameter::checkParameterAttachment
    }

    default {
    }
  }

  return $att_trg_ids
}


#################
# DETACH-button #
#################

# Procedure detaches the parameters from selected rows
# in the objects/elements list-box data.
#
proc StdPanelExec::detachParameter {globArray} {
  global Info ObjectTable
  upvar #0 $globArray theArray

  if {$theArray(hasBoundaries)} {
    set box_indices [$theArray(boundaryLB) curselection]
  } else {
    set box_indices [$theArray(objectLB) curselection]
  }

  if {$box_indices == "" } {
    return
  }

  # Find object (ObjectTable) ids
  foreach idx $box_indices {
    lappend obj_ids [lindex $theArray(targetIds) $idx]
  }

  set trg_ids [StdPanelExec::attachParamId $globArray $obj_ids $box_indices $Info(NO_INDEX)]

  set theArray(dataChanged) 1

  Widget::configureS $theArray(panelDetachButton) disabled

  #---Panel specific processing
  switch $globArray {

    BoundaryCondition {
      #BoundaryCondition::updateRadiationTargetInfo $trg_ids
    }

    Equation {
      Equation::detachProc $obj_ids $box_indices
    }
  }
}


#############
# SAVE proc #
#############

#-Standard panel save proc
#
#-Data will be saved to C++-variables!!!
#-Result data consist of element update data and parameter
# update data. This data is tranferred as a combined list into
# cpp-environment, where is is split and handled separately
# Also the ids of the possibly deleted parameters are sent to cpp.
#
proc StdPanelExec::panelSave { globArray {inform_front 1} } {
  global Model
  upvar #0 $globArray theArray

  #---Save data

  #-Create new backup vars
  if { [info exists theArray(ids,old)] } {
    Util::unsetArrayIdVariables $globArray $theArray(ids,old) *,old

    set theArray(ids,old) $theArray(ids)

    foreach id $theArray(ids) {
      StdPanelExec::createParameterBackupData $globArray $id
    }
  }

  #-Save model data in compressed mode
  set modified 0
  Panel::compressParameters $globArray modified

  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Util::cpp_exec $theArray(panelOkProc)

  #-In workspace we use uncompressed (including empty fields etc.) data
  Panel::uncompressParameters $globArray 

  #---Set/reset flags etc.
  StdPanelExec::resetUpdateLists $globArray
 
  return 1
}


proc StdPanelExec::createParameterBackupData { globArray id} {
  upvar #0 $globArray theArray
  
  foreach vn {data name oid mask} {
    if { [info exists theArray($id,$vn)] } {
      set theArray($id,$vn,old) $theArray($id,$vn)
    }
  }
}


#-Delete one parameter by parameter id
# All possible applications to "elements" are detached.
#
proc StdPanelExec::deleteParameterById {globArray pid} {
  global Info ObjectTable
  upvar #0  $globArray theArray 

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  #--Indices for the rows to be updated
  set obj_ids ""
  set box_indices ""

  foreach id $theArray(allTargetIds) {

    if { $pid == $ObjectTable($id,$pv) } {
      lappend obj_ids $id

      set idx [lsearch $theArray(targetIds) $id]

      if { $idx != -1 } {
        lappend box_indices $idx
      }
    }
  }

  #--Delete parameter-id value from relevant object/element lists.
  StdPanelExec::attachParamId $globArray $obj_ids $box_indices $Info(NO_INDEX)

  #--Delete parameter data
  StdPanelExec::deleteActiveParameterData $globArray $pid

  set index [lsearch $theArray(ids) $pid]

  set theArray(ids) [lreplace $theArray(ids) $index $index]

  #--Delete rows from list-box related lists
  if { [info exists theArray(parameterLBList)] } {
    set theArray(parameterLBList) [lreplace $theArray(parameterLBList) $index $index]
  }

  if { [info exists theArray(parameterLB)] &&
       [winfo exists $theArray(parameterLB)] 
     } {
    $theArray(parameterLB) delete $index
  }

  set theArray(parameterIndex) ""
  set theArray(dataChanged) 1
}


# Delete active variables for the paramter "pid2 in the paramter
# array
# NOTE: This does NOT delete ",old" values!!!
#
proc StdPanelExec::deleteActiveParameterData {globArray pid} {
  upvar #0 $globArray theArray

  foreach vn {data name oid mask} {
    catch { unset theArray($pid,$vn) }
  }
}


# Procedure check if selected targets (bodies or elements) have a parameter
# attached or not
#
proc StdPanelExec::selectedTargetsAreUnattached {globArray {all_must_be 1} } {
  global Common Info ObjectTable
  upvar #0  $globArray theArray 
  
  if { !$theArray(hasTarget) } {
    return 0
  }

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  set object_index [$theArray(objectLB) curselection]
  set obj_ids ""
  set box_indices ""

  #--Store the selection indices
  if { $theArray(hasBoundaries) } {
    set box_indices [$theArray(boundaryLB) curselection]

  } else {
    set box_indices [$theArray(objectLB) curselection]
  }
  

  #--If nothing was selected!
  if { ($object_index == "") && ($box_indices == "") } {
    return 0
  }

  set has_unattached 0
  set has_all_unattached 1

  # We have elements
  # ================
  if {$theArray(hasBoundaries)} {

    #-Object was selected, but no boundaries --> select all boundaries which
    # don't have a parameter defined!
    if { $box_indices == "" } {
      
      set indices ""
      set index 0
      foreach id $theArray(boundaryIds) {
        lappend indices
        incr index
      }
 
    } else {
      set indices $box_indices
    }

    foreach index $indices {

      set oid [lindex $theArray(boundaryIds) $index]

      if { ![info exists ObjectTable($oid,$pv)] ||
           $ObjectTable($oid,$pv) == $Info(NO_INDEX)
         } {
        set has_unattached 1

      } else {
        set has_unattached_all 0
      }
    }

  # We have only bodies
  # ===================
  } else {
    set oid $theArray(objectId)

    if { ![info exists ObjectTable($oid,$pv)] ||
         $ObjectTable($oid,$pv) == $Info(NO_INDEX)
       } {
      set has_unattached 1
    } else {
      set has_unattached_all 0
    }
  }

  # Nothing unattached
  if {!$has_unattached} {
    return 0
  }

  # Not all unattached
  if { $all_must_be && !$has_all_unattached } {
    return 0
  }

  # Ok, "suitably" unattched
  return 1
}

   

########################
# Values area handling #
########################

# Procedure checks if the object which is selcted in the objects
# list-box has 'real' (ie. not just itself) subelements
# 
proc StdPanelExec::selectedObjectHasElements {globArray selected_row} {
  global ObjectTable
  upvar #0 $globArray theArray

  #-If we can decide already on panel-level
  if {!$theArray(hasBoundaries)} {
    return 0
  }

  set oid [lindex $theArray(objectIds) $selected_row]

  if { $ObjectTable($oid,sbIds) == "" } {
    return 0

  } else {
    return 1
  }
}


# Pocedure finds the problem-type mask for the selected object
#
proc StdPanelExec::getSelectedObjectMask {globArray selected_row} {
  global Info ObjectTable 
  upvar #0 $globArray theArray

  set id [lindex $theArray(objectIds) $selected_row]

  #return $ObjectTable($id,msk)
  return [StdPanelExec::getTargetMask $globArray $id]

}


# Pocedure finds the problem-type mask for given target
# of the array defined in var "globArray"
# NOTE: This should be used when "running" the panels
# instead of reading directly ObjecTable(..,msk) variable
#
# When initializing the mask in the StdPanelInit-panel it is
# ok to use directly ObjectTable-mask, because we are there
# still building the object masks!
#
proc StdPanelExec::getTargetMask {globArray target_id} {
  global Equation Info ObjectTable 
  upvar #0 $globArray theArray

  # 
  if { $globArray == "GridH" } {
    # Vertices and edges are different!
    #if { $ObjectTable($target_id,tp) == "V" } {
    #  set mask "¤bcD"
    #} else {
    #  set mask "¤bcN"
    #}
    # Only geometry mask
    set mask $Equation(dimgType)

  } elseif { $globArray == "GridParameter" } {
    # This makes the parameter to match only one body!
    set mask1 [Object::getTypeMask $target_id]
    set mask2 [Object::getMask $target_id]
    #set mask [Util::combineMasks $mask1 $mask2]
    # Only geometry mask
    set mask $Equation(dimgType)

  # NOT in use!
  } elseif { $globArray == "BB_BoundaryCondition" } {
    set mask1 [Object::getTypeMask $ObjectTable($target_id,prId)]
    set mask2 [Object::getMask $target_id]
    set mask [Util::combineMasks $mask1 $mask2]

  } else {
    set mask [Object::getMask $target_id]
  }

  return $mask
}


# Pocedure finds the problem-type mask for given object
#
proc StdPanelExec::getObjectMask {id} {

  set mask1 [Object::getTypeMask $id]
  set mask2 [Object::getMask $id]

  return [Util::combineMasks $mask1 $mask2]
}


# Reset the value area's widgets into their initial values
#
proc StdPanelExec::clearValuesArea {globArray} {
  global Info Common
  upvar #0 $globArray theArray

  # TESTING if we really need this proc!!!
  # Apparently not!!!
  return

  foreach fld $theArray(allFields) {

    # Only displayed fields are interesting
    if { !$theArray($fld,dsp) } {
      continue
    }
    
    set ival [DataField::getInitialValue $globArray $fld]
    set theArray($fld) $ival

    set wtype [DataField::getFieldProperty $globArray $fld WidgetType]

    if { $wtype == "Text" } {
      set wdg $theArray(allWidgets,$fld)
      $wdg delete 0.0 end
      $wdg insert 0.0 $theArray($fld)
    }

    if { [Panel::isGroupField $globArray $fld]} {

      foreach gf $Common($theArray(parameterType),$fld,group) {

        set ival [DataField::getInitialValue $globArray $gf]
        set theArray($gf) $ival
      }
    }

  }
}


# Select panel variables by problem-type or set panel-specific
# default variables if there is nothing special for the current
# problem-type.
#
proc StdPanelExec::setValuesAreaActivity {globArray target_mask {forced_state ""} {msg "" } } {
  global Common Info Model $globArray
  upvar #0 $globArray theArray

  if { [info exists theArray(equationIndex)] } {
    set eq_index $theArray(equationIndex)

  } else {
    set eq_index ""
  }

  if { [info exists theArray(objectId)] } {
    set eq_ids [Object::getEquationIds $theArray(objectId)]
  } else {
    set eq_ids ""
  }

  set problem_mask [Panel::getProblemMask $globArray]

  set variables ""

  # What to do if no panel mask defined
  #
  if { $theArray(hasMaskedFields) && !$theArray(acceptEmptyMask) } {
    set accept_empty_mask 0
  } else {
    set accept_empty_mask 1
  }

  if { $target_mask == "" &&
       !$accept_empty_mask
     } {
    set empty_target_mask_deactivates 1
  } else {
    set empty_target_mask_deactivates 0
  }

  # Collect all active fields in the panel
  set panelFields ""
  set activityOnlyByMask ""

  foreach fld $theArray(allFields) {

    # Only problem fields are interesting here
    if { ![info exists theArray($fld,prb)] || !$theArray($fld,prb) } {
      continue
    }

    # Group fields
    # ============
    if { [Panel::isGroupField $globArray $fld] } {


      foreach gf [Panel::getGroupFields $globArray $fld] {

        set theArray($gf,act) $theArray($fld,act)

        # For group fields we have to check memebers separately, cause
        # parent possibly has not mask at all!
        if { ![DataField::fieldProblemMaskMatches $globArray $gf $problem_mask] } {
          set theArray($gf,act) 0
          continue
        }

        lappend panelFields $gf
      }

    # Single field
    # ============
    } else {
      lappend panelFields $fld
    }

  } ;# Collect panel fields

  
  # 1. Loop all case fields and check activity by mask
  # --------------------------------------------------
  #
  foreach fld $panelFields {

    # Only displayed fields
    if { ![info exists theArray($fld,dsp)] || !$theArray($fld,dsp) } {
      continue
    }

    # If some group field is non-packed!
    if { ![info exists theArray(allWidgets,$fld)] ||
         ![winfo exists $theArray(allWidgets,$fld)]
       } {
      continue
    }

    # Check activity by mask
    #
    if { $theArray(hasMaskedFields) &&
         [DataField::getFieldProperty $globArray $fld ActivityOnlyByMask $eq_index]
       } {

      if { ($target_mask == "" && $empty_target_mask_deactivates) ||
            ![DataField::fieldTargetMaskMatches $globArray $fld $target_mask] ||
            $forced_state == "disabled" 
          } {
        set theArray($fld,act) 0

      } else {
          set theArray($fld,act) 1
      }
    }
  } ;# for all fields checking activity by mask


  # 2. Loop all case fields and check activity by parent
  # ----------------------------------------------------
  foreach fld $panelFields {

    # Only displayed fields
    if { ![info exists theArray($fld,dsp)] || !$theArray($fld,dsp) } {
      continue
    }

    # If some group field is non-packed!
    if { ![info exists theArray(allWidgets,$fld)] ||
         ![winfo exists $theArray(allWidgets,$fld)]
       } {
      continue
    }

    # Check if some parents control field's activity
    #
    set is_abp [DataField::isActivatedByParents $globArray $fld $eq_index $eq_ids $problem_mask]

    if { $is_abp == -1 } {
      set theArray($fld,act) 0

    } elseif { $is_abp == 1 } {
      set theArray($fld,act) 1
    }
  } ;# for all fields checking activity by parent


  # Finally, set field activity
  # ---------------------------
  #
  foreach fld $panelFields {

    # Only displayed fields
    if { ![info exists theArray($fld,dsp)] || !$theArray($fld,dsp) } {
      continue
    }

    # If some group field is non-packed!
    if { ![info exists theArray(allWidgets,$fld)] ||
         ![winfo exists $theArray(allWidgets,$fld)]
       } {
      continue
    }

    Widget::setFieldStatus $globArray $fld
    Widget::configureField $globArray $fld

  } ;# for all fields set activity
 
}



# Set field status in the values area (modified, error etc)
#
proc StdPanelExec::setValuesAreaStatus {globArray {forced_status ""} } {
  global Common Info Model $globArray
  upvar #0 $globArray theArray

  set variables ""

  # Collect all active fields in the panel
  set panelFields ""

  foreach fld $theArray(allFields) {

    # Only displayed fields
    if { !$theArray($fld,dsp) } {
      continue
    }

    # Group fields
    if { [Panel::isGroupField $globArray $fld] } {

      foreach gf [Panel::getGroupFields $globArray $fld] {

        lappend panelFields $gf
      }

    # Single field
    } else {
      lappend panelFields $fld
    }
    
  }

  # Loop all active fields
  foreach fld $panelFields {

   # If some group field is non-packed!
   if { ![info exists theArray(allWidgets,$fld)] ||
         ![winfo exists $theArray(allWidgets,$fld)]
       } {
      continue
    }

    set wdg $theArray(allWidgets,$fld)

    if { [catch {set wstate [$wdg cget -state] }] } {
      continue
    }

    if { $wstate == "disabled" } {
      continue
    }

    set wtype [DataField::getFieldProperty $globArray $fld WidgetType]

    if { $forced_status != "" } {
      set theArray($fld,mod) $forced_status
    }

    #--Browsable directory widget
    if { $wtype == "BrowsableDirectory" } {
      Widget::setEntryStatus $globArray $fld $wdg

    #--Browsable file widget
    } elseif { $wtype == "BrowsableFile" } {
      Widget::setEntryStatus $globArray $fld $wdg

    #--Entry widget
    } elseif { $wtype == "Entry" } {
      Widget::setEntryStatus $globArray $fld $wdg
 
    #--Button widget
    } elseif { $wtype == "Button" } {
      #Widget::setButtonStatus $globArray $fld $wdg

    #--Checkbox widget
    } elseif { $wtype == "CheckBox" } {
      Widget::setCheckBoxStatus $globArray $fld $wdg

    #--Option menu widget
    } elseif { $wtype == "OptionMenu" } {
      Widget::setOptionMenuStatus $globArray $fld $wdg
 
    #--Radio-button widget
    } elseif { $wtype == "RadioButton" } {
      #Widget::setRadioButtonStatus $globArray $fld $wdg

    #--Text widgets
    } elseif { $wtype == "Text" } {

    #--Default (set state if supported)
    } else {
      # Do nothing
    }

    #--Procedure check-box
    if { [info exists theArray(allWidgets,$fld,procedure)] } {
      set wdg $theArray(allWidgets,$fld,procedure)
      Widget::setCheckBoxStatus $globArray $fld $wdg
    }

    #--File check-box
    if { [info exists theArray(allWidgets,$fld,file)] } {
      set wdg $theArray(allWidgets,$fld,file)
      Widget::setCheckBoxStatus $globArray $fld $wdg
    }
  }
}


# Procedure updates currentBody name(s) for the selected object
#
proc StdPanelExec::setCurrentBodyIdsAndNames {globArray selected_row} {
  global Info ObjectTable
  upvar #0 $globArray theArray

  set id [lindex $theArray(objectIds) $selected_row]

  # Body
  if { $ObjectTable($id,tp) == "B" } {
    set theArray(body1Id) $id
    set theArray(body2Id) $Info(NO_INDEX)

  # BodyLayer
  } elseif { $ObjectTable($id,tp) == "BL" } {
    set theArray(body1Id) $Info(NO_INDEX)
    set theArray(body2Id) $Info(NO_INDEX)

  # BodyPair
  } elseif { $ObjectTable($id,tp) == "BP" } {
    set theArray(body1Id)  $ObjectTable($id,pr1Id)
    set theArray(body2Id)  $ObjectTable($id,pr2Id)
  }

  set theArray(objectName) $ObjectTable($id,nm)
}


# Pocedure updates currentBoundary name for the selected element
#
proc StdPanelExec::setCurrentBoundaryName {globArray selected_row} {
  upvar #0 $globArray theArray

  set theArray(boundaryName) "Unknown!"
}


# Procedure updates parameter name in the entry field
# into currentParameterName variable
#
proc StdPanelExec::updateCurrentParameterName {globArray widget} {
  upvar #0 $globArray theArray

  # NOTE: Do not trim spaces here!!!
  set theArray(parameterName) [$widget get]
}


# Procedure updates parameter name in the entry field
#
proc StdPanelExec::storeCurrentParameterName {globArray} {
  upvar #0 $globArray theArray

  if { $theArray(parameterIndex) != "" } {

    # NOTE: No blank parameter name accepted!
    # This is mainly to guarantee that no blank bc-names are
    # delivered to Solver and that finally would create problems in ElmerPost
    #
    set theArray(parameterName) [string trim $theArray(parameterName) ]

    set pid $theArray(parameterId)
    set new_n [string trim $theArray(parameterName)]
    set old_n $theArray($pid,name)

    set theArray($pid,name) $new_n
 
    if {$new_n != $old_n} {
      set theArray(dataChanged) 1
    }
  }
}


proc StdPanelExec::markParamRowUpdated {globArray row} {
  upvar #0  $globArray theArray 

  set old_row [lindex $theArray(parameterUpdateList) $row]
  set new_row [lreplace $old_row 0 0 1]

  set theArray(parameterUpdateList) \
        [lreplace $theArray(parameterUpdateList) $row $row $new_row]
}


proc StdPanelExec::setNormalTangentialLabels {globArray field_name nt_flag} {
  global $globArray 
  upvar #0 $globArray theArray

  foreach i {1 2 3} c_s {N T T} {
    set FV $field_name
    append FV _$i
    set index [expr $i - 1]

    set label_wdg $theArray(allWidgets,label,$FV)

    set label_string [DataField::getFieldProperty $globArray $FV Label]

    if { $nt_flag == "True" } {
      set coord_string "-$c_s"
    } else {
      set coord_string [Util::getCoordinateLabel $i]
      set coord_string "-$coord_string"
    }

    set unit_string  [DataField::getFieldProperty $globArray $FV UnitString]

    set f_lbl "$label_string$coord_string $unit_string"

    Widget::configure1 $label_wdg -text $f_lbl
  }
}



#######################################
#### Data and ListBoxLists updates ####
#######################################

# Check if there is an applied parameter for a selected object/element
# If found: select the corresponding parameter row and make it active
# If not found: select the LAST possible parameter row and make it active
#
proc StdPanelExec::findAndActivateParameter {globArray target_id} {
  global Info ObjectTable
  upvar #0 $globArray theArray

  if { [info exists theArray(parameterLB)] &&
       [winfo exists $theArray(parameterLB)]
     } {
    $theArray(parameterLB) selection clear 0 end
  }

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  set pid $ObjectTable($target_id,$pv)

  # If no parameter set
  if { $pid == $Info(NO_INDEX) } {

    if { $theArray(hasMaskedParameter) } {

      set mask $theArray(targetMask)
      set npid [StdPanelExec::getNextParameterInstance $globArray $mask]

      if { $npid != $Info(NO_INDEX) } {
        set pid $npid
      }
    }
  }

  # Nothing was really found!
  if {$pid == $Info(NO_INDEX)} {
    return
  }

  set prow [lsearch $theArray(ids) $pid]

  execNsProc  $globArray parameterSelected $prow
  #StdPanelExec::parameterSelected $globArray $prow

  #Widget::configureS $theArray(panelUpdateButton) normal

  return


  set pname $theArray($pid,name)
  set pdata $theArray($pid,data)

  #--Update param-info-parameters
  set theArray(parameterIndex) $prow
  set theArray(parameterId) $pid
  set theArray(parameterName) $pname

  $theArray(parameterLB) selection clear 0 end

  #---Pick data from the selected parameter
  if {$prow != -1 } {
    DataField::setDataFields $globArray $pdata
    DataField::setAndShowLBSelection $theArray(parameterLB) $prow
  }
}


# Find and select the target object for the parameter pid
#
proc StdPanelExec::findAndActivateTarget {globArray pid {find_first 1} } {
  global Info ObjectTable
  upvar #0 $globArray theArray

  if { $pid == "" || $pid == $Info(NO_INDEX) } {
    return
  }

  #--Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  #--Pick current target id
  set current_trg_id $theArray(targetId)

  #--If we should select first available or nothing
  #  is currently selected
  #
  if { $find_first ||
       $current_trg_id == "" ||
       $current_trg_id == $Info(NO_INDEX)
     } {
    set current_trg_id 1
  }

  #--Find next available target for pid
  set trg_id ""

  set count [llength $ObjectTable(ids)]
  set index [lsearch $ObjectTable(ids) $current_trg_id]
  set counter 0

  while { $counter < $count } {

    incr index

    # If round end
    if { $index >= $count } {
      set index 0
    }

    set id [lindex $ObjectTable(ids) $index]

    if { [info exists ObjectTable($id,$pv)] &&
         $ObjectTable($id,$pv) == $pid
       } {
      set trg_id $id
      break
    }

    incr counter

  } ;#end while


  # If nothing found
  if { $trg_id == "" } {
    return
  }

  # Find parent (first) for object
  set parent_id [StdPanelExec::findObjectParent $trg_id]

  if { $parent_id == $Info(NO_INDEX) } {
    return
  }

  set data [lsearch $theArray(objectIds) $parent_id]

  if { $data == -1 } {
    return
  }
 
  lappend data $trg_id

  # Simulate object selection!
  set select_parameter 1
  execNsProc $globArray objectSelected [list $data $select_parameter]
}


# Not quite working?
# MVe 23.09.00
proc StdPanelExec::findObjectParent { obj_id } {
  global Info Model ObjectTable

  if { $obj_id == $Info(NO_INDEX) } {
    return $Info(NO_INDEX)
  }

  # Boundary
  if { [info exists ObjectTable($obj_id,prId)] &&
       $ObjectTable($obj_id,prId) != "" &&
       $ObjectTable($obj_id,prId) != $Info(NO_INDEX)
     } {
    return $ObjectTable($obj_id,prId)

  } elseif { [info exists ObjectTable($obj_id,pr1Id)] } {

    set tp $ObjectTable($obj_id,tp)
    set id1 $ObjectTable($obj_id,pr1Id)
    set id2 $ObjectTable($obj_id,pr2Id)
    return [Object::findIdByParentIds $tp $id1 $id2]

  # Edge (in 3D)
  } elseif { "E" == $ObjectTable($obj_id,tp) } {
    set face_id [StdPanelExec::findSubObjectParent $obj_id]
    return [StdPanelExec::findObjectParent $face_id]

  # Vertex
  } elseif { "V" == $ObjectTable($obj_id,tp) } {

    # 3D vertex
    if { $Model(GEOMETRY_DIMENSION) == "3D" } {

      set edge_id [StdPanelExec::findSubObjectParent $obj_id]
      set face_id [StdPanelExec::findSubObjectParent $edge_id]
      return [StdPanelExec::findObjectParent $face_id]
    
    # 2D vertex
    } else {
      set edge_id [StdPanelExec::findSubObjectParent $obj_id]
      return [StdPanelExec::findObjectParent $edge_id]
    }

  } else {
    return $Info(NO_INDEX)
  }
}


# Find sub object parent for a subobject like vertex or edge (in 3D)
#
proc StdPanelExec::findSubObjectParent { obj_id } {
  global Info ObjectTable

  if { $obj_id == $Info(NO_INDEX) } {
    return $Info(NO_INDEX)
  }

  foreach id $ObjectTable(ids) {

    # No subids
    if { ![info exists ObjectTable($id,sbIds)] ||
         $ObjectTable($id,sbIds) == ""
       } {
      continue
    }

    # Check if sub-ids in this object match the obj_id    
    foreach sid $ObjectTable($id,sbIds) {
      
      # This object is the sub-object-parent
      if { $sid == $obj_id } {
        return $id
      }
    }
  }

  return $Info(NO_INDEX)
}


proc StdPanelExec::getNextParameterInstance {globArray mask {current_pid ""} } {
  global Info
  upvar #0 $globArray theArray

  # If no current reference given, select last possible
  # =============================
  if { $current_pid == "" } {
    
    set npid $Info(NO_INDEX)

    foreach id $theArray(ids) {

      if { $theArray($id,mask) == $mask } {
        set npid $id
        ;#break ;# Uncomment this if you want the first!
      }
    }

    return $npid
  }
  
  # If current reference given, find "nearest"
  # ==========================
  set available_ids ""
  set index 0
  set counter 0

  foreach id $theArray(ids) {
    
    if { $theArray($id,mask) == $mask } {
      
      # Index of current (so we try to keep the current index!)
      if { $id == $current_pid} {
        set index $counter

      # NOTE: We do not add current to this list!
      } else {
        lappend available_ids $id
      }

      incr counter
    }
  }

  if { $available_ids == "" } {
    return $Info(NO_INDEX)
  }

  # If current was the last, correct index!
  set max_index [llength $available_ids]
  incr max_index -1

  if { $index > $max_index } {
    set index $max_index
  }
  set npid [lindex $available_ids $index]

  return $npid
}


# Make boundaries with pid attached as selected (via cpp!)
#
proc StdPanelExec::selectBoundaries {globArray pid} {
  global Info ObjectTable Model
  upvar #0 $globArray theArray

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  set Common(currentArray) $globArray

  set first 1
  set count 0
  set id_data ""

  # Find all objects where parameter-id (pid) is attached
  #
  foreach id $ObjectTable(ids) {

    if { ![info exists ObjectTable($id,$pv)] ||
         $ObjectTable($id,$pv) != $pid
       } {
      set ObjectTable($id,slctd) 0
      continue
    }

    if { $ObjectTable($id,slctd) } {
      continue
    }

    set oid $theArray(objectId)

    if { $ObjectTable($oid,tp) == "BP" } {
      set bd1_id $ObjectTable($oid,pr1Id)
      set bd2_id $ObjectTable($oid,pr2Id)

    } else {
      set bd1_id $oid
      set bd2_id $Info(NO_INDEX)
    }

    append id_data " $bd1_id $bd2_id $id"

    incr count
  }

  # Collect data
  set data $count
  append data $id_data

  # Flags: accept body change, do not update Gui
  append data " 1 1"

  # First reset all selections
  Util::cpp_exec resetAllBoundarySelections

  # Set extend mode
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_EXTEND 1

  if { $count > 0 } {
    Util::cpp_exec boundariesSelected $data
  } else {
    ListBox::markSelectedBoundaries
  }

  # Reset flags
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_EXTEND 0
}


proc StdPanelExec::updateEquationMenu {globArray} {
  global Common Equation Info
  upvar #0 $globArray theArray
  
  if { !$theArray(panelsByEquation) } {
    return
  }

  if { $Equation(activeIndices) == "" } {
    return
  }

  set cindex $theArray(equationIndex)

  if { $cindex == "" || $cindex == $Info(NO_INDEX) } {
    return
  }

  # Select equations in the menu
  #
  if { $globArray == $Common(panelArray,EQ) } {
    set eindices $Equation(ownIndices)

  } else {
    set eindices [Equation::getMatchingEquationIndices $theArray(targetMask)]
  }

  # If no equations defined, no parameter adding allowed!
  #
  if { [info exists theArray(panelAddButton)] &&
       [winfo exists $theArray(panelAddButton)]
     } {

    if { $eindices == "" } {
      $theArray(panelAddButton) configure -state disabled
    } else {
      $theArray(panelAddButton) configure -state normal
    }
  }

  # If we have some problem indices (body has equation) and
  # if current problem index is not any more active, make
  # the first of the active indices the current
  #
  if { $eindices != "" && -1 == [lsearch $eindices $cindex] } {
    set eindex [lindex $eindices 0]
    set theArray(equationIndex) $eindex
    set theArray(equationLabel) [Equation::getFieldPropertyByIndex $eindex EquationLabel]

    set theArray(currentEquationIndices) $eindices

    StdPanelCreate::constructPanelValuesArea $globArray $theArray(equationIndex)
  }

  # Put all active equations for the body into the menu
  StdPanelCreate::createEquationMenu $globArray $eindices

  StdPanelExec::setEqMoveButtonsState $globArray
}


proc StdPanelExec::equationMenuSelectionProc { globArray {eindex ""} } {
  global Common Equation
  upvar #0 $globArray theArray
  
  StdPanelCreate::constructPanelValuesArea $globArray $eindex

  # Make selected equation visible in the equation panel's
  # equation listbox
  #
  if { $globArray == $Common(panelArray,EQ) } {
    set ix [lsearch $Equation(ownIndices) $Equation(equationIndex)]
    $Equation(allWidgets,EQUATION_LIST) see $ix
  }

  # Set Next/Prev equation buttons state
  #
  StdPanelExec::setEqMoveButtonsState $globArray
}


#-Procedures sets data lists into their previous state after a cancel etc.
#
proc StdPanelExec::resetDataLists {globArray} {
  global ObjectTable
  upvar #0 $globArray theArray

  set theArray(objectIndex) $theArray(objectIndex,old)
  set theArray(boundaryIndex) $theArray(boundaryIndex,old)
  set theArray(parameterIndex) $theArray(parameterIndex,old)

  #-Parameters related stuff
  if { $theArray(hasParameters) } {

    set pv [string tolower $theArray(parameterType)]

    #-Here we restore applied parameter ids
    foreach id $theArray(allTargetIds) {
      set ObjectTable($id,$pv) $ObjectTable($id,$pv,old)
    }
    
    # Delete all active parameter data
    # ================================
    # NOTE: old-values are not deleted with these, because
    # variable names are given exactly, not like data* !!!
    #
    Util::unsetArrayIdVariables $globArray $theArray(ids) data
    Util::unsetArrayIdVariables $globArray $theArray(ids) name
    Util::unsetArrayIdVariables $globArray $theArray(ids) oid
    Util::unsetArrayIdVariables $globArray $theArray(ids) mask

    # Restore old values
    # ==================
    set theArray(ids) $theArray(ids,old)

    foreach pid $theArray(ids) {
 
      foreach vn { data name oid mask} {
 
        if { [info exist theArray($pid,$vn,old)] } {
          set theArray($pid,$vn) $theArray($pid,$vn,old)
        }
      }
    }

  }

}


#-Procedures resets update-indicator lists into initial state (a zero-list)
#
proc StdPanelExec::resetUpdateLists {globArray} {
  upvar #0 $globArray theArray

  set theArray(dataChanged) 0
  set theArray(parameterIsDeleted) 0
}


#-Procedure updates all necessary lists when a parameter is attached
# to objects/elements
# When pid (param-id) = Info(NO_INDEX), it is deatched!
# Return attached target ids: id > 0 attched, id < 0 detached
#
proc StdPanelExec::attachParamId {globArray target_ids box_indices pid {overwrite 1} } {
  global Info ObjectTable
  upvar #0 $globArray theArray

  if { $target_ids == "" } {
    return ""
  }

  set attach_target_ids ""

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  # NOTE: We use 'target'-variable-names here for convenience
  # This is mainly because Tcl-list-commands (lreplace etc)
  # don't update "in situ", but result a new list, which
  # must be stored by a "set" command

  set lb_wdg ""
  set target_list ""

  # Select proper list for update
  if {!$theArray(hasBoundaries)} {

    if { [info exists theArray(objectLB)] &&
         [winfo exists $theArray(objectLB)] 
       } {
      set lb_wdg $theArray(objectLB)
      append target_list "(objectLBList)"
     }
    
  } else {

    if { [info exists theArray(boundaryLB)] &&
         [winfo exists $theArray(boundaryLB)] 
       } {
      set lb_wdg $theArray(boundaryLB)
      append target_list "(boundaryLBList)"
     }
  }

  set was_attached 0

  #---Loop all target_ids
  set counter 0
  foreach id $target_ids {

    # If we do not overwrite an existing def.
    # NOTE: This check is needed for automatic attachment
    # when adding a parameter
    if { !$overwrite &&
         $ObjectTable($id,$pv) != "" &&
         $ObjectTable($id,$pv) != $Info(NO_INDEX)
       } {
      continue
    }

    #--Attach parameter

    set was_attached 1

    # Attached
    if { $pid != $Info(NO_INDEX) } {
      lappend attach_target_ids $id

    # Detached
    } else {
      lappend attach_target_ids [expr -1 * $id]
    }

    set ObjectTable($id,$pv) $pid
  }
  
  # If some target list is active
  # =============================
  if { $target_list != "" } {

    set targetList $globArray
    append targetList $target_list
  
    upvar #0 $targetList theTargetList

    #---Loop all box_indices
    foreach idx $box_indices {
  
      if { [info exists theTargetList] } {

        set oid [lindex $theArray(targetIds) $idx]

        set lb_row [StdPanelInit::formTargetLBRow $globArray $oid]

        set theTargetList [lreplace $theTargetList $idx $idx $lb_row]
      }

      # If no listbox visible
      if { $lb_wdg == "" } {
        continue
      }
    
      #--Update listbox
      set lb_selections [$lb_wdg curselection]

    
      # Check if listbox-row was selected 
      if { -1 != [lsearch $lb_selections $idx] } {
        set select true
      } else {
        set select false
      }

      # Update listbox-row
      $lb_wdg delete $idx $idx
      $lb_wdg insert $idx $lb_row
    
      if {$select} {
        $lb_wdg selection set $idx $idx
      }
    }
  } ;# if target list active


  if {$was_attached} {

    #--Mark data updated
    set theArray(dataDirty) 0

    #--Mark data changed
    set theArray(dataChanged) 1
  }

  return $attach_target_ids
}
    

############
### Misc ###
############


# Check that the new parameter line being added/updated is ok
#
proc StdPanelExec::checkNewParameterLine {globArray mode} {
  global Info
  upvar #0 $globArray theArray

  
  # Pick only active fields
  #
  set tmp $Info(keepInactiveFields)
  set Info(keepInactiveFields) 0
  set data [DataField::formActiveParameterLine $globArray  "allFields" $theArray(targetMask)]
  set Info(keepInactiveFields) $tmp

  # Data ok (not empty)
  # ==================
  if { $data != "" } {
    return 1
  }

  # Active data would be empty
  # ==========================
  if { $mode == "ADD" } {
    set msg [list "This would add an empty parameter set!\n"]

  } elseif { $mode == "UPDATE" } {
    set msg [list "This update would create an empty parameter set.\n"  \
                  "Delete the parameter, if you want to destroy it!" ] 
  } else {
    set msg [list "This would create an empty parameter set!\n"]
  }

  set Info(messageIcon) error
  Message::dispOkMessage $msg "$Info(FRONT_NAME) message" $theArray(winName)

  return 0
}


proc StdPanelExec::detachInconsistentAttachments { globArray pid modified_flag} {
  global Info ObjectTable
  upvar #0 $globArray theArray

  upvar $modified_flag modified

  set modified 0

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  set pmask $theArray($pid,mask)

  foreach id $theArray(allTargetIds) {

    if { ![info exists ObjectTable($id,$pv)] || 
         $ObjectTable($id,$pv) != $pid
       } {
      continue
    }

    #set omask $ObjectTable($id,msk)
    set omask [StdPanelExec::getTargetMask $globArray $id]

    if { ![Util::masksAreEqual $pmask $omask] } {

      set modified 1  
      set ObjectTable($id,$pv) $Info(NO_INDEX)

      set msg "***WARNING $ObjectTable($id,nm): "
      append msg "$theArray(defaultName) $pid detached!"
      Message::showMessage $msg 
    }
  }

}
        

# Make parameter ids consequtive ie. 1 2 3 ...
# Relevant when parameters have been deleted!
#
proc StdPanelExec::compressParameterIds {globArray} {
  global Info Model ObjectTable
  upvar #0 $globArray theArray

  if {!$theArray(hasParameters)} {
    return
  }

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]


  # Form updated parameter data list where parameter
  # ids are replaced with the "compressed" ids
  #
  set new_ids ""
  set new_param_list ""

  # Update also current parameter id!
  if { [info exists theArray(parameterId)] } {
    set old_current_id $theArray(parameterId)
  } else {
    set old_current_id $Info(NO_INDEX)
  }

  set compressed 0
  set new_current_id $Info(NO_INDEX)

  set counter 1
  foreach old_id $theArray(ids) {
  
    # New id
    set new_id $counter

    if { $new_id != $old_id } {

      # Mark modified!
      set compressed 1

      set theArray($new_id,data) $theArray($old_id,data)
      set theArray($new_id,name) $theArray($old_id,name)
      set theArray($new_id,mask) $theArray($old_id,mask)

      if { $theArray(hasTarget) } {
        set theArray($new_id,oid) $theArray($old_id,oid)
      }

      #--Check if old name was an "automatic" name, if so, it will
      #  be also "renumbered"!!!

      set old_def_name [Panel::defaultParameterName $globArray $old_id]
      set new_def_name [Panel::defaultParameterName $globArray $new_id]

      if { [string equal $theArray($old_id,name) $old_def_name] } {
        set theArray($new_id,name) $new_def_name
      }

      #--Check if this was current
      if { $old_current_id == $old_id } {
        set new_current_id $new_id
      }

      #--Remove old-id vars
      # NOTE: This is save, because old_id is NOT needed anymore
      StdPanelExec::deleteActiveParameterData $globArray $old_id

      #--Update target array id-pointers
      foreach id $theArray(allTargetIds) {

        if { $ObjectTable($id,$pv) == $old_id } {
          set ObjectTable($id,$pv) $new_id
        }
      }
    }

    #--Add new id to ids-list
    lappend new_ids $counter
    lappend new_param_list [StdPanelInit::formParameterLBRow $globArray $counter]

    incr counter
  }

  if { !$compressed } {
    return
  }

  set theArray(ids) $new_ids
  set theArray(parameterLBList) $new_param_list
  set theArray(nextNewParameterId) $counter

  if { $new_current_id != $Info(NO_INDEX) } {
    set theArray(parameterId) $new_current_id
  }

  #--Update also list-box lists
  if { [info exists theArray(objectLB)]   &&
       [winfo exists $theArray(objectLB)] 
  } {

    StdPanelInit::constructLBLists $globArray
    StdPanelInit::fillLBLists $globArray

    if { $theArray(hasBoundaries) } {
      set tmp $theArray(boundaryIndex)
      StdPanelExec::objectSelected $globArray $theArray(objectIndex)
      set theArray(boundaryIndex) $tmp
      StdPanelExec::boundarySelected $globArray $theArray(boundaryIndex)
    } else {
      StdPanelExec::objectSelected $globArray $theArray(objectIndex)
    }
  }
}   


proc StdPanelExec::removeGridDefinitions { globArray new_index_array } {
  global ObjectTable
  upvar #0  $globArray theArray 
  upvar $new_index_array theIndices 

  set pv [string tolower $theArray(parameterType)]

  #-Parameter ids varname 
  set g_var $pv
  append g_var "Ids"

  #-Mesh indices varname
  set m_var $pv
  append m_var "MshIndcs"

  foreach id $ObjectTable(ids) {

    if { ![info exist ObjectTable($id,$g_var)] ||
         $ObjectTable($id,$g_var) == ""
       } {
      continue
    }

    set new_gids ""
    set new_mids ""

    set gids $ObjectTable($id,$g_var)
    set mids $ObjectTable($id,$m_var)

    set index 0
    foreach gid $gids  mid $mids {

      set new_index $theIndices($mid)

      # Remove
      if { $new_index == $Info(NO_INDEX) } {
        continue
      }

      # Copy
      lappend new_gids $gid
      lappend new_mids $new_index
    }

    set ObjectTable($id,$g_var) $new_gids
    set ObjectTable($id,$m_var) $new_mids
  }
}


# Proc marks a parameter row updated
#
proc StdPanelExec::parameterRowWasUpdated {globArray row} {
  upvar #0  $globArray theArray 

  # We store old parameter row for restoring (after cancel etc..)
  # in the parameter uppdate-list. Note, this is done only once
  # to store the situation when the panel was opened!
  if {$theArray(dataChanged) == 1} {
    return
  }

  set old_pr_row [$theArray(parameterLB) get $row $row]
  
  set old_upd_row [lindex $theArray(parameterUpdateList) $row]
  set new_upd_row [list 1 $old_pr_row]

  set theArray(parameterUpdateList) [lreplace $theArray(parameterUpdateList) \
                                  $row $row $new_upd_row]
}




###########################
### EQUATION panel procs ###
###########################

proc Equation::checkVariablesInParameters {eindex removed_var_names} {
  global Equation

  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  set act_names_var [DataField::getFieldProperty Equation $ename EquationVars]
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]

  foreach pid $Equation(ids) {

    set act_names [Equation::getActiveVariableNames $eindex $pid]

    set still_act_names ""

    foreach name $act_names {

      if { "" == [Equation::getMatchingVariableName $name $removed_var_names] } {
        lappend still_act_names $name
      }
    }

    set still_act_names [join $still_act_names $sep]

    DataField::setFieldValue Equation $pid $act_names_var $still_act_names
  }
}


# Create the Equation problem mask (based on masks of all bodies)
#
proc Equation::constructProblemMask {} {
  global Common Equation Model ObjectTable

  set mask ""

  # Generic model level masks
  if { [info exists Equation(dimgType)] && $Equation(dimgType) != "" } {
    append mask "($Equation(dimgType))"
  }

  if { [info exists Equation(dimsType)] && $Equation(dimsType) != "" } {
    append mask "($Equation(dimsType))"
  }

  #--Loop all bodies
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    # Get mask for the body
    set body_mask [Object::getMask $id]

    if { $body_mask == "" } {
      continue
    }

    # If part of the body mask is not yet in the total mask
    # add it
    foreach sub_mask $body_mask {

      if { ![Util::patternMatchesString "($sub_mask)" $mask] } {
        append mask "($sub_mask)"
      }
    }
  }
 
  set Equation(problemMask) [string trim $mask]
}


# Return list of system variable names and their activity indicators like:
# { {"Flow Solution" 1} } for navier-Stokes system or
# { {"Nitrogen" 1} {"Ogygen" 0} ...} for Advection-Diffusion system etc
#
# Argument equation_index: equation to be inspected
#          equation_ids: equation parameter ids to be searched 
#
proc Equation::findSystemNames {eindex equation_ids} {
  global Equation

  if { ![info exists Equation(ids)] || $Equation(ids) == "" } {
    return ""
  }

  #-Field name for the equation itself
  #
  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  # If UserDefined equation has been removed!
  if { $ename == "" } {
    return ""
  }

  #-Get all potential variable (subsystem) names for the equation
  #

  set all_names [Equation::getAllVariableNames $eindex]

  #-Init systems list
  set system_list ""

  if { $all_names != "" } {
    foreach all_name $all_names {
      lappend system_list [list $all_name 0]
    }
  } else {
    lappend system_list [list "" 0]
  }

  #-Loop all equations ids
  #
  foreach pid $equation_ids {

    #-Equation parameter is not applied
    if { ![Object::parameterIsApplied $pid "eq" "B"] } {
      continue
    }
    
    #-Equation field is not active
    if { "True" != [DataField::getFieldValue Equation $pid $ename] } {
      continue
    }

    #-Equation is True and no sub-vars
    # NOTE: We add the equation name (base name) below!
    #
    if { $all_names == "" } {
      set system_list [list [list "" 1]]
    }

    #-Update all active variable slots in the result
    #
    set act_names [Equation::getActiveVariableNames $eindex $pid]

    foreach act_name $act_names {
      
      set idx [List::nocaseSearch $all_names $act_name]
      set row [lindex $system_list $idx]
      set act [lindex $row 1]

      if { $act == 0 } {
        set act 1
      }

      set row [list $act_name $act]

      set system_list [lreplace $system_list $idx $idx $row]
    }
  }

  #-Add base name to sub-var names

  # Base name like "Flow Solution" for NS or "" for AdvDiff
  set base [DataField::getFieldProperty Equation $ename Variable]

  set result ""
  foreach system $system_list {

    set name [lindex $system 0]
    set true_val [lindex $system 1]

    set name [string trim "$base $name"]

    lappend result [list $name $true_val]
  }

  return $result
}


# Find active variable names for the equation parameter
#
proc Equation::getActiveVariableNames {eindex {pid ""} {only_for_true_eq 1} } {
  global Equation EquationVariable

  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  #--If equation itself must be true!
  #
  #-If equation parameter id given, read result from it
  if { $pid != "" } {
    set eq_true_val [DataField::getFieldValue Equation $pid $ename]
  
  #-Otherwise use the current field value
  } else {
    set eq_true_val $Equation($ename)
  }

  if { $only_for_true_eq &&
       $eq_true_val != "True"
     } {
    return ""
  }

  #--Next read possible sub-vars
  
  #-Variable which stores active sub-var names
  set act_var_name [DataField::getFieldProperty Equation $ename EquationVars]

  #-If equation parameter id given, read sub-vars from it
  if { $pid != "" } {
    set act_names [DataField::getFieldValue Equation $pid $act_var_name]
  
  #-Otherwise use the current field value
  } elseif { [info exists Equation($act_var_name)] } {
    set act_names $Equation($act_var_name)

  } else {
    set act_names ""
  }

  #-Make the result list
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]

  return [split $act_names $sep]
}


proc Equation::getAllVariableNames {eindex} {
  global Equation EquationVariable

  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  #--Base name like "Flow Solution"
  set base_name [DataField::getFieldProperty Equation $ename Variable]

  #--All possible sub-var names
  set pid $EquationVariable(parameterId)

  set all_names [DataField::getFieldValue EquationVariable $pid $ename]

  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep "" 0]

  return [split $all_names $sep]
}


# Get mask for the equation by the index
#
proc Equation::getEquationMask { index } {
  global Equation

  # This can be set already here!
  set extMask "$Equation(dimgType) $Equation(dimsType) $Equation(bcType)"
  set eqMask [Equation::getFieldPropertyByIndex $index EquationMask]

  return "$eqMask $extMask"
}


proc Equation::getFieldPropertyByIndex { equation_index field_name
                                         {system_index ""}
                                         {show_error_msg 1}
                                         {found_flag_var ""}
                                       } {
  global Equation
  
  set found 0

  foreach vname $Equation(allFields) {

    set idx [DataField::getFieldProperty Equation $vname EquationIndex "" 0 found]

    if { !$found || $idx != $equation_index } {
      continue
    }

    return [DataField::getFieldProperty Equation $vname $field_name \
                                        $system_index \
                                        $show_error_msg \
                                        $found_flag_var ]
  }
}


proc Equation::getMatchingEquationIndices { mask } {
  global Equation
    
  set pindices ""

  foreach fld $Equation(allFields) {

    set eindex [DataField::getFieldProperty Equation $fld EquationIndex "" 0]

    if { $eindex == "" } {
      continue
    }

    set emask [DataField::getFieldProperty Equation $fld EquationMask]

    set match_found 0

    foreach emsk $emask {

      if { 1 == [Util::patternMatchesString $emsk $mask] } {
        set match_found 1
        break
      }
    }

    if { $match_found } {
        lappend pindices $eindex
    }

  }

  return $pindices
}


proc Equation::getMatchingVariableName {test_name test_name_list} {
  
  if { $test_name_list == "" } {
    return ""
  }
  
  set test_name [string tolower $test_name]

  foreach name $test_name_list {
    
    set name_ [string tolower $name]

    if { $test_name == $name_ } {
      return $name
    }
  }

  return ""
}


# Check is an equation panel field refers to an Equation
# name like HEAT_EQUATION
#
proc Equation::isEquation {fld} {
  global Equation

  if { "" != [DataField::getFieldProperty Equation $fld EquationIndex "" 0] } {
    return 1
  } else {
    return 0
  }
}


proc Equation::setUserDefinedEquationMask {eq_name} {
  global Equation

  set msk [DataField::getFieldProperty Equation $eq_name EquationMask]

  set msk [string triml $msk "("]
  set msk [string trimr $msk ")"]

  if { $Equation($eq_name) == "True" } {

    set Equation(userDefinedEquationType) [Util::combineMasks $msk $Equation(userDefinedEquationType)]
  } else {
    set Equation(userDefinedEquationType) [Util::removeMask $Equation(userDefinedEquationType) $msk]
  }

  Equation::updateTargetMask
}


proc Equation::setUserDefinedFieldsActivity {equation} {
  global Equation

  set eq_mask [DataField::getFieldProperty Equation $equation EquationMask]
  
  if { $Equation($equation) == "True" } {
    set activity 1
  } else {
    set activity 0
  }
      
  foreach fld $Equation(allFields) {

    if { 1 != [DataField::getFieldProperty Equation $fld IsUserDefinedField "" 0] ||
         ![Equation::isEquation $fld]
       } {
      continue
    }

    if { $Equation($fld,dsp) &&
         [DataField::fieldTargetMaskMatches Equation $fld $eq_mask]
       } {
      set Equation($fld,act) $activity
    }

  }
}


proc Equation::setVariableWidgetsActivity {equation_index} {
  global Equation

  set ename [Equation::getFieldPropertyByIndex $equation_index EquationField]

  set is_multi_var [DataField::getFieldProperty Equation $ename IsMultiVar "" 0]

  if { !$is_multi_var } {
    set Equation(VARIABLE_BUTTON,act) 0
    set Equation(VARIABLE_LIST,act) 0

  } else {
    set Equation(VARIABLE_BUTTON,act) 1
  
    if { $Equation($ename) != "True" } {
    set Equation(VARIABLE_LIST,act) 0

    } else {
      set Equation(VARIABLE_LIST,act) 1
    }
  }

  Widget::configureField Equation VARIABLE_BUTTON

  # NOTE: Listbox does not have a -state option!
  Widget::configureField Equation VARIABLE_LIST
}


proc Equation::setVariableWidgetsVisibility { {eindex ""}} {
  global Equation

  if { $eindex == "" } {
    if { [info exists Equation(equationIndex)] } {
      set eindex $Equation(equationIndex)
    } else {
      return
    }
  }

  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  set is_multi_var [DataField::getFieldProperty Equation $ename IsMultiVar "" 0]
  
  set all_names [Equation::getAllVariableNames $eindex]

  #-Hide variable-list widget if not multivar or
  # there is only one variable!!!
  #
  if { !$is_multi_var ||
       [llength $all_names] < 2
     } {
    DataField::setFieldProperty Equation VARIABLE_LIST Display 0
    set Equation(VARIABLE_LIST,scr) 0

  #-Show variable-list widget
  } else {
    DataField::setFieldProperty Equation VARIABLE_LIST Display 1
    set Equation(VARIABLE_LIST,scr) 1
  }
}


# Find active equations from body equation ids
# NOTE: Only "main" equation system indices are stored here!
# AdvectionDiffusion type subsystems (like ADVECTION_DIFFUSION_EQUATION_vars==Oxygen;Nitroge)
# are not indexed!
#
proc Equation::updateEquationIndices {} {
  global Info Solver Equation ObjectTable

  set Equation(activeIndices) ""

  set body_eq_ids ""

  # Loop all bodies
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    set eq_id [Object::getParameterId $id "eq"]

    if { $eq_id == $Info(NO_INDEX) } {
      continue
    }

    if { -1 == [lsearch $body_eq_ids $eq_id] } {
      lappend body_eq_ids $eq_id
    }
  }

  foreach eq_id $body_eq_ids {

    foreach fld $Equation(allFields) {

      set eindex [DataField::getFieldProperty Equation $fld EquationIndex "" 0]

      if { $eindex == "" } {
        continue
      }
      
      if { "True" != [DataField::getFieldValue Equation $eq_id $fld] } {
        continue
      }

      if { -1 == [lsearch $Equation(activeIndices) $eindex] } {
        lappend Equation(activeIndices) $eindex
      }
    }
  }

  set Equation(activeIndices) [lsort $Equation(activeIndices)]
}


# Update Equation panel equation list
#
proc Equation::updateEquationList {} {
  global Info Equation

  set Equation(equationList) ""

  foreach eidx $Equation(ownIndices) fld $Equation(ownNames) {

    set eq_row ""

    if { "True" == $Equation($fld) } {
      append eq_row $Info(selectionBoxTrueMarker)
    } else {
      append eq_row $Info(selectionBoxFalseMarker)
    }

    append eq_row " "
    append eq_row [Equation::getFieldPropertyByIndex $eidx Label]

    lappend Equation(equationList) $eq_row
  }

  set current [$Equation(allWidgets,EQUATION_LIST) curselection]

  if { [info exists Equation(allWidgets,EQUATION_LIST)] &&
       [winfo exists $Equation(allWidgets,EQUATION_LIST)]
     } {
    ListBox::fill $Equation(allWidgets,EQUATION_LIST) $Equation(equationList)
  }
}


# Update Equation array active-variables field (like ADVECTION_DIFFUSION_EQUATION_vars)
# when EquationVariables panel is saved
#
proc Equation::updateEquationVariables { old_all_names added_all_names removed_all_names } {
  global Equation equationVariable

  set eindex $Equation(equationIndex)
  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]

  #-Add new vars to the current case
  #
  if { $added_all_names != "" } {

    set act_names [Equation::getActiveVariableNames $eindex]

    if { $act_names == "" &&
         1 == [llength $old_all_names]
       } {
      set act_names $old_all_names
    }
    
    foreach new_name $added_all_names {
      lappend act_names $new_name
    }

    set act_names_var [DataField::getFieldProperty Equation $ename EquationVars]
    set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]
    set Equation($act_names_var) [join $act_names $sep]
  }

  #-Check variable list's visibility for the equation
  Equation::setVariableWidgetsVisibility $eindex

  #-We reconstruct panel although there is no equation change!
  set force_construct 1
  StdPanelCreate::constructPanelValuesArea Equation $eindex $force_construct

  #-Check all equation parameters for removed variables
  Equation::checkVariablesInParameters $eindex $removed_all_names
}



# Procedure updates equation-panel problem-type mask
#
proc Equation::updateTargetMask {} {
  global Equation Info

  set sep $Info(maskSepMark)

  #--Flow mask
  set m ""
  append m $Equation(flow_type)
  append m $sep
  append m $Equation(flow_laminarType)
  append m $sep
  append m $Equation(flow_turbulentType)
  set Equation(flowType) $m

  #--Heat mask
  set m ""
  append m $Equation(heat_type)
  append m $sep
  append m $Equation(heat_transferType)
  append m $sep
  append m $Equation(heat_phaseType)
  set Equation(heatType) $m

  #--Stress mask
  set m ""
  append m  $Equation(stress_type)
  append m $sep
  append m  $Equation(stress_mechanicalType)
  append m $sep
  append m  $Equation(stress_thermalType)
  set Equation(stressType) $m

  #--Advection-Diffusion mask
  set m ""
  append m  $Equation(diffusion_type)
  append m $sep
  append m  $Equation(diffusion_transferType)
  append m $sep
  set Equation(diffusionType) $m

  #--Total mask (excluding DIM masks!!!)
  set m $Equation(flowType)
  append m $sep
  append m $Equation(heatType)
  append m $sep
  append m $Equation(stressType)
  append m $sep
  append m $Equation(diffusionType)
  append m $sep
  append m $Equation(userDefinedEquationType)
  append m $sep
  append m $Equation(dimsType)

  #--Add user defined equations variable masks
  foreach fld $Equation(allFields) {
    if { [info exists Equation($fld,variableType)] } {
      append m $sep
      append m $Equation($fld,variableType)
    }
  }

  set m [string trim $m]

  # Replace separator with space, so that mask
  # is easy to split into a list of submask
  regsub -all ($sep)+ $m " " m

  set Equation(targetMask) [string trim $m]

  #--Possible values area state update
  if { $Equation(updateValuesArea) } {
    StdPanelExec::setValuesAreaActivity Equation $Equation(targetMask)
  }

} 


# Update Equation panel variable list for the
# selected equation
#
proc Equation::updateVariableList { {eindex ""} } {
  global Info Equation EquationVariable

  if { $eindex == "" } {
    set eindex $Equation(equationIndex)
  }

  set ename [Equation::getFieldPropertyByIndex $eindex EquationField]
  set eval $Equation($ename)

  #-Sub-var variable (field name)
  set act_names_var [DataField::getFieldProperty Equation $ename EquationVars]

  #-Find currently active variable names
  # NOTE: We are here really interested only on the
  # variable names, equation true-value is not relevant(?!!!)
  #
  set pid ""
  set only_for_true_eq 0
  set act_names [Equation::getActiveVariableNames $eindex $pid $only_for_true_eq]

  #-All possible variable names
  set all_names [Equation::getAllVariableNames $eindex]

  # If we have only one possible sub-var
  # ====================================
  #
  if { [llength $all_names] == 1 } {
    
    #-Turn (the only variable) off
    if { $eval != "True" } {
      set Equation($act_names_var) ""

    #-Turn on
    } else {
      set Equation($act_names_var) [lindex $all_names 0]
    }

    return
  }

  # If an other than current equation was updated
  # =============================================
  if { $eindex != $Equation(equationIndex) } {
    return
  }

  if { !$Equation(VARIABLE_LIST,dsp) } {
    return
  }

  # Current equation was changed, update also the listbox list
  # ==========================================================

  set variableList ""
  set variableLBList ""


  # Loop all potential variable names
  #
  foreach var_name $all_names {

    set var_row ""

    # Mark variable to be "True" (selected) and append to the
    # *_vars list
    #
    if { "" != [Equation::getMatchingVariableName $var_name $act_names] } {

      append var_row $Info(selectionBoxTrueMarker)

      lappend variableList $var_name
      
    } else {
      append var_row $Info(selectionBoxFalseMarker)
    }

    append var_row $var_name
    lappend variableLBList $var_row
  }

  #-Variables list in parameter
  set sep [DataField::getFieldProperty EquationVariable $ename FieldDataSep]
  set Equation($act_names_var) [join $variableList $sep]

  #-Variables listbox list
  ListBox::fill $Equation(allWidgets,VARIABLE_LIST) $variableLBList
}


# Update Equation variable type mask list for the whole problem
# Read all equation parameters
# NOTE: This is not actually needed, because totalProblemMask
# is constructed from body masks, and body masks are constructed
# from applied equation parameter target masks!!!
#
proc Equation::updateVariableTypeProblemMasks {} {
  global Info Equation EquationVariable

  foreach ename $Equation(allFields) {

    set eindex [DataField::getFieldProperty Equation $ename EquationIndex "" 0]
    
    if { $eindex == "" } {
      continue
    }

    set Equation($ename,variableType) ""

    #-Only actual multi-variable equations need this mask type!
    # (We can skip here Navier-Stokes, heat Equation etc.)
    #
    if { ![DataField::getFieldProperty Equation $ename IsMultiVar] } {
      continue
    }

    set emask [DataField::getFieldProperty Equation $ename EquationMask]
    set emask [string trimleft [string trimright $emask ")"] "("]
    
    set all_act_names ""

    foreach pid $Equation(ids) {
      
      #-Equation parameter is not applied
      if { ![Object::parameterIsApplied "eq" $pid "B"] } {
        continue
      }

      #-Currently active variable names ("True")
      set act_names [Equation::getActiveVariableNames $eindex $pid]
     
      foreach act_name $act_names {

        if { -1 == [List::nocaseSearch $all_act_names $act_name] } {
          lappend all_act_names $act_name
        }
      }

    } ;# For all parameters

    # Add all active variable masks for the equation
    foreach act_name $all_act_names {

      # Append sub-var mask to the mask list
      set msk $emask
      append msk [DataField::formVariableMask $act_name]
      append Equation($ename,variableType) " $msk"
    }

  } ;# For all equations

}


# Update Equation variable type mask list for the current target
# Read all current equation fields
#
proc Equation::updateVariableTypeTargetMasks {} {
  global Info Equation EquationVariable

  foreach ename $Equation(allFields) {

    set eindex [DataField::getFieldProperty Equation $ename EquationIndex "" 0]
    
    if { $eindex == "" } {
      continue
    }

    set Equation($ename,variableType) ""

    #-Only actual multi-variable equations need this mask type!
    # (We can skip here Navier-Stokes, heat Equation etc.)
    #
    if { ![DataField::getFieldProperty Equation $ename IsMultiVar] } {
      continue
    }

    set emask [DataField::getFieldProperty Equation $ename EquationMask]
    set emask [string trimleft [string trimright $emask ")"] "("]
    
    #-Currently active variable names ("True")
    set act_names [Equation::getActiveVariableNames $eindex]

    # Add all active variable masks for the equation
    foreach act_name $act_names {

      # Append sub-var mask to the mask list
      set msk $emask
      append msk [DataField::formVariableMask $act_name]
      append Equation($ename,variableType) " $msk"
    }

  } ;# For all equations
}


# Check if variable name shoud be used in the equation name
# "Advection Diffusion Equation Oxygen"
#
proc Equation::useVariableNameInEquationName { ename } {
  global Equation EquationVariable

  if { ![Equation::isEquation $ename] } {
    return 0
  }

  set eindex [DataField::getFieldProperty Equation $ename EquationIndex]

  set names [Equation::getAllVariableNames $eindex]

  # If only one name defined, do not use the variable name in equation
  # 
  if { [llength $names] > 1 } {
    return 1

  } else {
    return 0
  }
}
 

# end ecif_tk_panelExec.tcl
# ********************


