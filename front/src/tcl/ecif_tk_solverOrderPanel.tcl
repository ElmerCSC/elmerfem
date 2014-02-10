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
#Module:    ecif_tk_solverOrderPanel.tcl
#Language:  Tcl
#Date:      20.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the solving order
#
#************************************************************************
 

#------Active solvers and solving order definition  proc------


# This procedure displays the solver solving order definition panel
#
proc SolverOrder::openPanel { } {
  global Info Common Solver SolverOrder SolverSystem
  upvar #0 Solver theArray

  set w $SolverOrder(winName)
  set wgeom $SolverOrder(winGeometry)

  set id [winfo atom $w]
  set Solver(orderWinId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow SolverOrder $id $SolverOrder(winTitle) $wgeom] } {
    return
  }  
   
  set this $w
  toplevel $w
  focus $w 

  wm title $w $SolverOrder(winTitle)
  wm geometry $w $wgeom

  set SolverOrder(dataChanged) 0
  set SolverOrder(dataModified) 0

  # Init vars and store old values
  SolverOrder::initData

  # Ups, no data!
  if { $SolverOrder(ids) == "" } {
    Message::dispOkMessage "No active solvers or calculators!"
    Panel::cancel $w
    return
  }

  foreach id $SolverOrder(ids) {
    set val $SolverOrder($id,SOLVING_ORDER)
    set SolverOrder($id,SOLVING_ORDER,old) $val
  }

  set SolverOrder(sortedIds,old) $SolverOrder(sortedIds)

  #-----WIDGET CREATION
  #-Data is organized into two frames.
  #-Solver names and solving order numbers

  # All active systems will be inserted  into this frame 
  frame $w.f1 -relief groove -bd 2
  frame $w.f11  ;# labels
  frame $w.f12  ;# check boxes

  #-Ok Cancel Apply Buttons
  frame $w.fB ;#--Buttons

  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)
 
  #---Solver name and order entry fields (create and pack)

  set SolverOrder(widgets) ""

  set m $w.f1

  # Form widget(s) in sorted order
  #
  set counter 0
  set max_label_len 0

  foreach id $SolverOrder(sortedIds) {

    set pid $SolverOrder($id,pid)
    set eq_index $SolverOrder($id,equationIndex)
    set ord $SolverOrder($id,SOLVING_ORDER)

    incr counter
    set frm [frame $m.f$counter]
    set lbl_frm [frame $frm.lf]
    set wdg_frm [frame $frm.wf]
    
    #-Calculators
    if { $SolverOrder($id,type) == "CL" } {

      set name [DataField::getFieldValue Calculator $pid EQUATION]

    #-Solvers
    } else {

      set name [DataField::getFieldValue Solver $pid EQUATION]
      
      set ename [Equation::getFieldPropertyByIndex $eq_index EquationField]
        
      if { [Equation::useVariableNameInEquationName $ename] } {
        append name " "
        append name [DataField::getFieldValue Solver $pid VARIABLE]
      }
    }
    
    set len [string length $name]

    if { $len > $max_label_len } {
      set max_label_len $len
    }

    set lbl [label $lbl_frm.l$counter -text $name]

    set order_var "SolverOrder($id,SOLVING_ORDER)"
    set wdg [entry $wdg_frm.e$counter -textvariable $order_var -width 5]
    bind $wdg <KeyRelease> "Panel::panelDataChanged 1  SolverOrder $wdg {%A %K}"
    bind $wdg <KeyPress-Return> "+SolverOrder::checkPanelData"

    # We use this to order the widget rows in the panel
    lappend SolverOrder(widgets) [list $ord $id $frm $lbl_frm $wdg_frm $lbl $wdg]

    #pack $lbl $wdg -side left -padx $fpx1
    #pack $frm -side top -pady $fpy2
  }

  SolverOrder::packWidgets $max_label_len

  pack $w.f1 -side top  -padx $fpx3 -pady $fpy3 -expand 1 -fill both

  #-Buttons (create and pack)

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "SolverOrder::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "SolverOrder::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command SolverOrder::panelApply \
                                 -state $ap]

  button $w.fB.default -text Default -command SolverOrder::panelDefault -state normal

  focus $ok_btn
  set SolverOrder(applyButton)  $ap_btn
  set SolverOrder(cancelButton) $cn_btn


  pack $ok_btn $cn_btn $ap_btn $w.fB.default -side left -padx $fpx1 

  pack $w.fB -side top  -padx $fpx3 -pady $fpy3
  
}
# End proc solverOrderPanel


proc SolverOrder::initData { {set_to_defaults 0} } {
  global Calculator Info Solver SolverOrder

  set SolverOrder(ids) ""
  set id 0

  if { [info exists Solver(ids)] } {
    set solver_ids $Solver(ids)
  } else {
    set solver_ids ""
  }

  # Read solver data
  # ================
  foreach pid $solver_ids {

    set act [DataField::getFieldValue Solver $pid ACTIVE]

    if { $act != "True" } {
      continue
    }
    
    incr id
    
    lappend SolverOrder(ids) $id
      
    set eq_nm [DataField::getFieldValue Solver $pid EQUATION]
    set eq_var [DataField::fieldNameSifToGui $eq_nm]
    set eq_idx [DataField::getFieldProperty Equation $eq_var EquationIndex]

    if { $set_to_defaults } {
      set ord [DataField::getInitialValue Solver SOLVING_ORDER $eq_idx]

    } else {
      set odr [DataField::getFieldValue Solver $pid SOLVING_ORDER]
    }

    set SolverOrder($id,pid) $pid
    set SolverOrder($id,type) "SL"
    set SolverOrder($id,SOLVING_ORDER) $odr
    set SolverOrder($id,equationIndex) $eq_idx
  }

  # Read calculator data
  # ====================
  set large_value 1000000

  if { [info exists Calculator(ids)] } {
    set calculator_ids $Calculator(ids)
  } else {
    set calculator_ids ""
  }

  foreach pid $calculator_ids {

    set act [DataField::getFieldValue Calculator $pid ACTIVE]

    if { $act != "True" } {
      continue
    }
    
    incr id
    lappend SolverOrder(ids) $id
  
    if { $set_to_defaults } {
      set ord $large_value
      incr large_value

    } else {
      set odr [DataField::getFieldValue Calculator $pid SOLVING_ORDER]
    }

    set SolverOrder($id,pid) $pid
    set SolverOrder($id,type) "CL"
    set SolverOrder($id,SOLVING_ORDER) $odr
    set SolverOrder($id,equationIndex) $Info(NO_INDEX)
  }

  SolverOrder::updateSolvingOrder
}


# Form a consequtive index series (1,2,3,..) from
# the solving-order numbers
#
proc SolverOrder::updateSolvingOrder {} {
  global SolverOrder

  set SolverOrder(sortedIds) ""
  set large_value 100000

  # This is a list of (id,solving-order-nbr)-pairs
  #
  set sort_list ""
  
  foreach id $SolverOrder(ids) {

    set ord [string trim $SolverOrder($id,SOLVING_ORDER)]

    if { $ord == "" } {
      set ord $large_value
      incr large_value
    }

    lappend sort_list [list $id $ord]
  }

  # This is the same, but sorted by solving-order number.
  # Now the row-index gives the new order number and first
  # number in the (sorted) pairs is the id
  #
  set sorted_list [lsort -index 1 -real $sort_list]

  # Finally form a consequtive sequence (1 2 3 ...) from
  # the order numers which can be non-unique!
  #
  set counter 1
  foreach pair $sorted_list {

    set id [lindex $pair 0]

    lappend SolverOrder(sortedIds) $id

    # Data was modified!
    if { $SolverOrder($id,SOLVING_ORDER) != $counter } {
      Panel::panelDataChanged 1 SolverOrder
    }

    set SolverOrder($id,SOLVING_ORDER) $counter

    incr counter
  }
}


proc SolverOrder::updateSolversAndCalculators {} {

  SolverOrder::initData
  SolverOrder::panelSave
}


proc SolverOrder::packWidgets { {label_width ""} } {
  global SolverOrder Solver Info

  # Set current order number to the 1. column
  set new_list ""

  foreach item $SolverOrder(widgets) {
    set id [lindex $item 1]
    set item [lreplace $item 0 0 $SolverOrder($id,SOLVING_ORDER)]
    lappend new_list $item
  }

  set SolverOrder(widgets) $new_list

  # Sort list according to the order number 
  set SolverOrder(widgets) [lsort -index 0 -real $SolverOrder(widgets)]

  # Pack widget in new order
  foreach item $SolverOrder(widgets) {

    set frm [lindex $item 2]
    set lbl_frm [lindex $item 3]
    set ent_frm [lindex $item 4]
    set lbl [lindex $item 5]
    set ent [lindex $item 6]

    if { [winfo exists $frm] } {
      pack forget $frm
    }
    
    if { $label_width != "" } {
      $lbl configure -width $label_width
    }

    pack $frm -side top -pady $Info(framePadY2) -expand 1 -fill both -anchor w
    pack $lbl_frm $ent_frm -side left -padx $Info(framePadX1) -anchor w -fill x -expand 0
    pack $lbl -side left -padx $Info(framePadX1) -anchor w -expand 0
    pack $ent -side left -padx $Info(framePadX1) -anchor w
  }
}


proc SolverOrder::panelSave {} {
  global Info Solver SolverOrder

  Panel::panelDataChanged 0 SolverOrder 
  Panel::panelDataModified 0 SolverOrder

  set Model(Front,needsUpdate) 1

  #--Store old values
  foreach id $SolverOrder(ids) {
    set SolverOrder($id,SOLVING_ORDER,old) $SolverOrder($id,SOLVING_ORDER)
  }

  set update_solvers 0
  set update_calculators 0

  #--Update solver and calulator paramters
  foreach id $SolverOrder(ids) {
    
    set ord $SolverOrder($id,SOLVING_ORDER)
    set pid $SolverOrder($id,pid)

    if { $SolverOrder($id,type) == "SL" } {
      DataField::setFieldValue Solver $pid SOLVING_ORDER $ord
      set update_solvers 1

    } else {
      DataField::setFieldValue Calculator $pid SOLVING_ORDER $ord
      set update_calculators 1
    }
  }

  #--Save solver data
  if { $update_solvers } {
    Solver::panelSave
  }

  #--Save calculator panel
  if { $update_calculators } {
    StdPanelExec::panelSave "Calculator"
  }
}


proc SolverOrder::panelOk {w} {
  global Solver Model SolverOrder

  #---No changes
  if { !$SolverOrder(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![SolverOrder::checkPanelData] } {
    return
  }

  #---Save data
  SolverOrder::panelSave

  Panel::cancel $w
}


proc SolverOrder::panelApply {} {
  global Info Solver SolverOrder

  #---No changes
  if { !$SolverOrder(dataChanged) } {
    return
  }

  if { ![SolverOrder::checkPanelData] } {
    return
  }
   
  SolverOrder::panelSave
}


proc SolverOrder::panelDefault {} {
  global Info Solver SolverOrder SolverSystem

  SolverOrder::initData

  Panel::panelDataChanged 1  SolverOrder
}


proc SolverOrder::panelCancel {w} {
  global Info Solver SolverOrder

  if { ![Panel::verifyCancel SolverOrder] } {
    return
  }

  #---Reset into old values
  foreach id $SolverOrder(ids) {
    set SolverOrder($id,SOLVING_ORDER) $SolverOrder($id,SOLVING_ORDER,old)
  }

  set SolverOrder(sortedIds) $SolverOrder(sortedIds,old)

  Panel::cancel $w
}


# Return 1 = ok, -1 = error
#
proc SolverOrder::checkPanelData {} {
  global SolverOrder
  #-Check solver order fields

  set min_value -1000000000

  # Replace empty values with minimum values
  #
  foreach id $SolverOrder(ids) {

    set number [string trim $SolverOrder($id,SOLVING_ORDER)]

    if { $number == "" } {
      set SolverOrder($id,SOLVING_ORDER) $min_value
      incr min_value
    }
  }

  # Form sorted ids
  SolverOrder::updateSolvingOrder

  SolverOrder::packWidgets

  return 1
}


# end ecif_tk_solverOrderPanel.tcl
# ********************
