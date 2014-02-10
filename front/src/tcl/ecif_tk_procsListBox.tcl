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
#Module:    ecif_tk_procsListBox.tcl
#Language:  Tcl
#Date:      05.10.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Listbox handling procedures
#
#************************************************************************

###########################
### List-Box procedures ###
###########################


#-Fill list-box widget with dataList-elements.
proc ListBox::fill {lbox dataList} {

  set start_pos [lindex [$lbox yview] 0]

  eval {$lbox delete 0 end}
  eval {$lbox insert end} $dataList

  $lbox yview moveto $start_pos
}


# Set the problem type indicator at the beginning of a parameter row
#
proc ListBox::markParameterLBByMask { globArray } {
  global Common Info ObjectTable
    upvar #0 $globArray theArray

  
  set index 0
  set new_list ""

  set tmask $theArray(targetMask)

  foreach id $theArray(ids) {

    # Get paramter mask via parameter's object
    #set mask [StdPanelExec::getTargetMask $globArray $theArray($id,oid)]
    set mask $theArray($id,mask)

    set accept_target 1

    if { [info exists theArray(hasParameterTargetCheckProc)] &&
         $theArray(hasParameterTargetCheckProc)
       } {
      set pn $globArray
      append pn "::parameterTargetCheckProc"
      set accept_target [eval $pn $id]
    }

    if { $accept_target && $mask == $tmask } {
      set mark $Info(paramMarked)

    } else {
      set mark $Info(paramUnmarked)
    }

    set row [lindex $theArray(parameterLBList) $index]
    set row [string replace $row 0 0 $mark]
    lappend new_list $row

    incr index
  }

  set theArray(parameterLBList) $new_list

  set selections [$theArray(parameterLB) curselection]

  ListBox::fill $theArray(parameterLB) $new_list

  foreach idx $selections {
    $theArray(parameterLB) selection set $idx
  }
}


proc ListBox::markSelectedBoundaries { {source_list "" } {current_array ""} } {
  global Info ObjectTable Common ModelFlags

  if { $current_array == "" } {
    set current_array $Common(currentArray)
  }

  # If no active panel specified
  if {$current_array == ""} {
    return 
  }

  upvar #0 $current_array theArray

  if { ![info exists theArray(boundaryLB)] ||
       ![winfo exists $theArray(boundaryLB)]
     } {
    return
  }

  if { $source_list == "" } {
    set source_list [ $theArray(boundaryLB) get 0 end]
  }

  set lbr 0
  set mark $Info(selectedMark)

  set current_selections [$theArray(boundaryLB) curselection]

  foreach old_row $source_list {

    if { $mark == [string index $old_row 0] } {
      set has_mark 1
    } else {
      set has_mark 0
    }

    #-The row without the selection mark
    set clean_row [string trimleft $old_row $mark]

    #-Check if row should be selected
    set bndr_id [lindex $theArray(boundaryIds) $lbr]
    set selected $ObjectTable($bndr_id,slctd)

    #-Insert selection mark in the beginning
    if { $selected } {
      set new_row $mark
      append new_row $clean_row

    #-No selection mark
    } else {
      set new_row $clean_row
    }

    if { $selected && !$has_mark ||
         !$selected && $has_mark
       } {
      ListBox::updateRow $theArray(boundaryLB) $lbr $new_row $selected
    }

    incr lbr
  }

  # Restore current selections
  foreach lbr $current_selections {
    $theArray(boundaryLB) selection set $lbr
  }
}


# Make a body selected in the objects list-box (from cpp)
# NOTE mode-variable not currently in use!!!
#
proc ListBox::selectBody { bd1_id lr1_id bd2_id lr2_id {mode ""} } {
  global Common Info ObjectTable

  set glob_array $Common(currentArray)

  # If currently no active panel
  if { $glob_array == "" } {
    return
  }

  upvar #0 $glob_array theArray

  if { ![info exists theArray(objectLB)] ||
       ![winfo exists $theArray(objectLB)]
     } {
    return
  }

  # Bodypair
  if { $bd1_id != $Info(NO_INDEX) &&
       $bd2_id != $Info(NO_INDEX) &&
       $bd1_id != $bd2_id
     } {
    
    set obj_id [Object::findIdByParentIds "BP" $bd1_id $bd2_id]
    set lr_id $lr1_id

  # First body
  } elseif { $bd1_id != $Info(NO_INDEX) } {
    set obj_id $bd1_id
    set lr_id $lr1_id

  # Second body
  } else {
    set obj_id $bd2_id
    set lr_id $lr2_id
  }
  
  # Try the body first
  set obj_index [lsearch $theArray(objectIds) $obj_id]

  # Ok, perhaps objects are layers (like GridParameter)
  if { $obj_index == -1 } {

    set obj_index [lsearch $theArray(objectIds) $lr_id]
    
    # Nothing found
    if { $obj_index == -1 } {
      return
    }
  }

  $theArray(objectLB) selection clear 0 end
  $theArray(objectLB) selection set $obj_index
  $theArray(objectLB) see $obj_index

  # Exec array object selected proc
  execNsProc $glob_array objectSelected $obj_index
}


# Selection in the elements list box (from cpp)
#
# NOTE: Work only with arries which have the following array variables defined:
#
# boundaryLb, boundaryIds, boundaryIndex
#
#
proc ListBox::selectBoundary {bndr_id bd1_id lr1_id bd2_id lr2_id {extend 0} } {
  global Common ModelFlags ObjectTable

  set glob_array $Common(currentArray)
  set trg_id $bndr_id

  # If currently no active panel
  if { $glob_array == "" } {
    return 
  }

  upvar #0 $glob_array theArray

  if { ![info exists theArray(boundaryLB)] ||
       ![winfo exists $theArray(boundaryLB)]
     } {
    return
  }


  if { $glob_array == $Common(panelArray,BC) } {
    if { [info exists ObjectTable($bndr_id,grpId)] } {
      set trg_id $ObjectTable($bndr_id,grpId)
    }
  }

  # Select first correct body, if not
  # yet current
  if { $theArray(body1Id) != $bd1_id ||
       $theArray(body2Id) != $bd2_id
     } {
    ListBox::selectBody $bd1_id $lr1_id $bd2_id $lr2_id
  }


  # Clear old selections, if not extended selection mode
  if { !($extend && $theArray(canMultiSelectBoundaries)) } {
     $theArray(boundaryLB) selection clear 0 end
  }
 
  #-Boundary listbox row to select
  set bndr_lb_row [lsearch $theArray(boundaryIds) $trg_id]
  $theArray(boundaryLB) selection set $bndr_lb_row

  # If not currently marked selected
  # ================================
  if { !$ObjectTable($trg_id,slctd) } {

    set $ObjectTable($trg_id,slctd) 1

  # Currently marked selected
  # =========================
  } else {
    set $ObjectTable($trg_id,slctd) 0
    #$theArray(boundaryLB) selection clear $bndr_lb_row
  }

  # Is focus has changed  
  if { $theArray(boundaryIndex) != $bndr_lb_row } {

    # Exec proper array boundary selected proc
    execNsProc $glob_array boundarySelected $bndr_lb_row
  }
}



#-Adds or removes (toggles) a selection mark (like #)
# at the end of of a listbox entry.
# Return new status (0 = doesn't have a mark, 1 = has a mark)
#
proc ListBox::toggleSelectionMark {source row {select 1} } {
  global Info
 
  set mark $Info(selectedMark)

  set data [join [$source get $row $row]]
  set data [string trimright $data]

  set mark_index [string first $mark $data]
  if {$mark_index == -1} {
    append data $mark
      set hasMark 1
  } else {
    set data [string trimright $data $mark]
      set hasMark 0
  }

  $source delete $row
  $source insert $row $data
  
  if { $select } {
    $source selection set $row
  }

  return $hasMark
}


# Updates one listbox row
proc ListBox::updateRow {lbox index row {set_selected 1} } {

  set start_pos [lindex [$lbox yview] 0]

  $lbox delete $index $index
  $lbox insert $index $row

  $lbox yview moveto $start_pos

  if { $set_selected } {
    $lbox selection set $index
  }
}


# end ecif_tk_procsListBox.tcl
# ********************
