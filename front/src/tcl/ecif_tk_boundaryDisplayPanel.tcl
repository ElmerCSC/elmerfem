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
#Module:    ecif_tk_boundaryDisplayPanel.tcl
#Language:  Tcl
#Date:      14.04.00
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for selecting boundaries to be displayed
#
#************************************************************************


#--Select boundaries to be displayed
#
proc BoundaryDisplay::openPanel {} {
  global BoundaryDisplay Info ObjectTable Model
  upvar #0 BoundaryDisplay theArray

  set w $BoundaryDisplay(winName)
  set wgeom $BoundaryDisplay(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set BoundaryDisplay(winId) $id

  set Info(thisWindow) $w
  set this $w

  if { 1 == [Util::checkPanelWindow BoundaryDisplay $id "Select boundaries" $wgeom] } {
    raise $BoundaryDisplay(winName)
    focus -force $BoundaryDisplay(winName)
    return
  }  

  set BoundaryDisplay(dataChanged) 0
  set BoundaryDisplay(dataModified) 0

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $BoundaryDisplay(winTitle)
  wm geometry $w $wgeom

  # Pick boundary object ids
  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    set BoundaryDisplay(ids) [Object::getIds "F"]
  } else {
    set BoundaryDisplay(ids) [Object::getIds "E"]
  }

  # Current display status
  foreach id $BoundaryDisplay(ids) {
    set BoundaryDisplay($id) [Object::getDisplayed $id]
    set BoundaryDisplay($id,old) [Object::getDisplayed $id]
    set BoundaryDisplay($id,slctd) [Object::getSelected $id]
  }

  #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 [expr 3 * $Info(framePadY3)]

  #---Outer frames
  set f1 [frame $w.f1]     ;#outer frame
  set f11 [frame $w.f1.f1] ;#-check box + boundary name
  set f12 [frame $w.f1.f2] ;#-All None button area
  set f2 [frame $w.f2] ;    #Apply+Ok+cancel buttons frame


  #-Boundaries listbox
  set wdg [ listbox $f11.lb -relief sunken \
             -selectmode browse -exportselection 0 \
             -height 20 -width 40 -font $Info(tableFont) \
             -xscrollcommand [list $f11.sx set] \
             -yscrollcommand [list $f11.sy set] ]

  scrollbar $f11.sx -orient horizontal -command [list $wdg xview]
  scrollbar $f11.sy -orient vertical -command [list $wdg yview]

  set BoundaryDisplay(boundaryLB) $wdg
  
  BoundaryDisplay::update

  pack $f11.sx -side bottom -fill x -expand 0 
  pack $f11.lb -side left -fill both -expand 1
  pack $f11.sy -side left -fill y -expand 0

  set BoundaryDisplay(boundaryLB) $wdg

  bind $wdg <ButtonRelease-1> "BoundaryDisplay::setDisplayMode {ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <ButtonRelease-2> "BoundaryDisplay::setDisplayMode {ButtonRelease-2 %x %y %X %Y}"
  bind $wdg <ButtonRelease-3> "BoundaryDisplay::setDisplayMode {ButtonRelease-3 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-1> "BoundaryDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-2> "BoundaryDisplay::setDisplayMode {Control-ButtonRelease-2 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-3> "BoundaryDisplay::setDisplayMode {Control-ButtonRelease-3 %x %y %X %Y}"
  bind $wdg <KeyPress-Return> "BoundaryDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <KeyPress-space> "BoundaryDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"

  bind $wdg <Double-1> "BoundaryDisplay::select"
  bind $wdg <Control-Double-1> "BoundaryDisplay::select"
  bind $wdg <Control-KeyPress-Return> "BoundaryDisplay::select"
  bind $wdg <Control-KeyPress-space> "BoundaryDisplay::select"

  #-All, None buttons
  set all_btn [button $f12.all -text All -command "BoundaryDisplay::all" \
                            -width 5]
  set none_btn [button $f12.none -text None -command "BoundaryDisplay::none" \
                            -width 5]

  set select_btn [button $f12.select -text Select -command "BoundaryDisplay::select" \
                            -width 5]

  pack $all_btn $none_btn -side top -expand 0 -pady $fpy3
  pack $select_btn -side top -expand 0 -pady [expr 2 * $fpy3]

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  #-Buttons
  set ok_btn [button $f2.ok -text OK -command "BoundaryDisplay::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "BoundaryDisplay::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f2.apply -text Apply -command "BoundaryDisplay::panelApply" \
                               -state $ap]
  
  focus $ok_btn
  set BoundaryDisplay(applyButton) $ap_btn
  set BoundaryDisplay(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 0 -padx $fpx2
  #pack $ok_btn $cn_btn -side left -expand 1 -padx $fpx1
  
  pack $f1 -side top -anchor w -expand 1 -fill both -padx $fpx2 -pady $fpy2
  pack $f11 -side left -expand 1 -fill both -padx $fpx2
  pack $f12 -side left -expand 0 -padx $fpx2
  pack $f2 -side top -expand 0 -padx $fpx2 -pady $fpy2
}


proc BoundaryDisplay::formBoundaryLBRow {id} {
  global BoundaryDisplay Info ObjectTable

  if { $BoundaryDisplay($id) } {
    set row $Info(selectionBoxTrueMarker)

  } else {
    set row $Info(selectionBoxFalseMarker)
  }

  append row "("
  append row $ObjectTable($id,tg)
  append row ")"

  append row $ObjectTable($id,nm)

  if { $BoundaryDisplay($id,slctd) } {

    append row " \[#\] "
  }

  return $row
}


proc BoundaryDisplay::all {} {
  global BoundaryDisplay ObjectTable

  set changed 0

  set index 0

  foreach id $BoundaryDisplay(ids) {
    
    if { $BoundaryDisplay($id) != 1 } {
      set changed 1

      set BoundaryDisplay($id) 1

      set ObjectTable($id,dspl) $BoundaryDisplay($id)

      set row [BoundaryDisplay::formBoundaryLBRow $id]
      ListBox::updateRow $BoundaryDisplay(boundaryLB) $index $row
    }

    incr index
  }

  $BoundaryDisplay(boundaryLB) selection clear 0 end

  if {$changed} {
    Panel::panelDataChanged 1 BoundaryDisplay
    BoundaryDisplay::updateDisplay
  }
}


proc BoundaryDisplay::none {} {
  global BoundaryDisplay ObjectTable

  set changed 0

  set index 0

  foreach id $BoundaryDisplay(ids) {
    
    if { $BoundaryDisplay($id) != 0 } {

      set changed 1

      set BoundaryDisplay($id) 0

      set ObjectTable($id,dspl) $BoundaryDisplay($id)

      set row [BoundaryDisplay::formBoundaryLBRow $id]

      ListBox::updateRow $BoundaryDisplay(boundaryLB) $index $row

    }

    incr index
  }

  $BoundaryDisplay(boundaryLB) selection clear 0 end

  if {$changed} {
    Panel::panelDataChanged 1 BoundaryDisplay
    BoundaryDisplay::updateDisplay
  }
}


proc BoundaryDisplay::select {} {
  global BoundaryDisplay ObjectTable

  set sindex [$BoundaryDisplay(boundaryLB) curselection]

  if { $sindex == "" || $sindex < 0 } {
    return
  }

  set id [lindex $BoundaryDisplay(ids) $sindex]

  set prntId $ObjectTable($id,prId)

  set bd1 $ObjectTable($prntId,pr1Id)
  set bd2 $ObjectTable($prntId,pr2Id)
  set lr1 -1
  set lr2 -1
  set bndr $id
  set accept_body_change 1
  set update_gui 1

  #--NOTE: Apply boundary selection via cpp !!!
  #
  set data "$bndr $bd1 $lr1 $bd2 $lr2 $accept_body_change $update_gui"
  Util::cpp_exec boundarySelected $data

  $BoundaryDisplay(boundaryLB) selection set $sindex
}


#
proc BoundaryDisplay::setDisplayMode {einfo} {
  global BoundaryDisplay ObjectTable

  set key [lindex $einfo 0]

  # Selection coordinates
  #-relative box coordinates in pixels
  set x [lindex $einfo 1]
  set y [lindex $einfo 2]

  #-absolute screen coordinates in pixels
  set X [lindex $einfo 3]
  set Y [lindex $einfo 4]

  #-Button-1
  if { $key == "ButtonRelease-1" } {
    
    #-'Near' the X-mark
    if { $x != "" && $x < 25 } {
      set key "ButtonRelease-3"

    #-Otherwise nothing for Button-1
    } else {
      return
    }
  }

  #-Make clicked row active, manipulate clicked
  #
  if { $key == "ButtonRelease-2" ||
       $key == "ButtonRelease-3" 
     } {
    set sindex [$BoundaryDisplay(boundaryLB) index @$x,$y]
    set windex [$BoundaryDisplay(boundaryLB) index @$x,$y]
    set do_select 1

  #-Keep currently selected row, manipulate currently selected
  } elseif { $key == "Control-ButtonRelease-1"
     } {
    set sindex [$BoundaryDisplay(boundaryLB) curselection]
    set windex [$BoundaryDisplay(boundaryLB) curselection]
    set do_select 0

  #-Keep currently selected row, manipulate clicked
  } elseif { $key == "Control-ButtonRelease-2" ||
             $key == "Control-ButtonRelease-3" 
     } {
    set sindex [$BoundaryDisplay(boundaryLB) curselection]
    set windex [$BoundaryDisplay(boundaryLB) index @$x,$y]
    set do_select 0
  }

  if { $windex == "" || $windex < 0 } {
    return
  }

  if { $do_select } {
    $BoundaryDisplay(boundaryLB) selection clear 0 end
  }

  set id [lindex $BoundaryDisplay(ids) $windex]

  set dmode $BoundaryDisplay($id)

  # Toggle selection
  if { $dmode == 1 } {
    set BoundaryDisplay($id) 0
  } else {
    set BoundaryDisplay($id) 1
  }

  set ObjectTable($id,dspl) $BoundaryDisplay($id)

  set row [BoundaryDisplay::formBoundaryLBRow $id]

  ListBox::updateRow $BoundaryDisplay(boundaryLB) $windex $row 0

  # NOTE: When listbox has focus, then Return/Space seem to
  # advance the selected row, this fixes this!
  #
  if { $sindex != "" && $sindex >= 0 } {
    $BoundaryDisplay(boundaryLB) selection set $sindex
    $BoundaryDisplay(boundaryLB) activate $sindex
  }

  Panel::panelDataChanged 1 BoundaryDisplay

  BoundaryDisplay::updateDisplay
}


proc BoundaryDisplay::updateDisplay {} {
  global BoundaryDisplay ObjectTable

  # Refresh boundary display
  Util::updateGui
  Util::cpp_exec boundaryDisplayPanelOk
}


proc BoundaryDisplay::panelOk {w} {
  global BoundaryDisplay

  BoundaryDisplay::panelApply
  
  Panel::cancel $w
}


proc BoundaryDisplay::panelApply {} {
  global BoundaryDisplay ObjectTable

  #---No changes
  if { !$BoundaryDisplay(dataChanged) } {
    return
  }

  if { -1 == [BoundaryDisplay::panelCheck] } {
    return
  }
  
  foreach id $BoundaryDisplay(ids) {
    set ObjectTable($id,dspl) $BoundaryDisplay($id)
    set BoundaryDisplay($id,old) $BoundaryDisplay($id)
  }

  Panel::panelDataChanged 0 BoundaryDisplay 
  Panel::panelDataModified 0 BoundaryDisplay 

  BoundaryDisplay::updateDisplay
}


proc BoundaryDisplay::panelCancel {w} {
  global BoundaryDisplay ObjectTable

  if {$BoundaryDisplay(dataChanged)} {

    foreach id $BoundaryDisplay(ids) {
      set BoundaryDisplay($id) $BoundaryDisplay($id,old)
    }

    # Refresh display
    BoundaryDisplay::panelApply
  }

  Panel::cancel $w
}


proc BoundaryDisplay::panelCheck {} {
  global BoundaryDisplay

  #-Ok
  return 1
}


proc BoundaryDisplay::update {} {
  global BoundaryDisplay

  set boundary_list ""

  foreach id $BoundaryDisplay(ids) {
    set row [BoundaryDisplay::formBoundaryLBRow $id]
    lappend boundary_list $row
  }

  ListBox::fill $BoundaryDisplay(boundaryLB) $boundary_list
}


# NOTE: This is called from cpp via
# Interface::setBoundarySelectionMode
# when a geometry boundary is selected
#
proc BoundaryDisplay::setSelectionMode {id mode do_update} {
  global BoundaryDisplay

  set idx [lsearch $BoundaryDisplay(ids) $id]

  if { -1 == $idx } return

  set BoundaryDisplay($id,slctd) $mode

  if { !$do_update } return

  BoundaryDisplay::update
  
  # Check that lb-row is visible (relevant when
  # selected from geometry!)
  #
  set y_pos [$BoundaryDisplay(boundaryLB) yview]
  set sz [expr 1.0 * [$BoundaryDisplay(boundaryLB) size] - 1]

  set cur_pos [expr $idx / $sz]
  set y_start [lindex $y_pos 0]
  set y_end [lindex $y_pos 1]

  if { $cur_pos < $y_start || $cur_pos > $y_end } {
    $BoundaryDisplay(boundaryLB) see $idx
  }
}


#ecif_boundaryDisplayPanel.tcl

