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
#Module:    ecif_tk_vertexDisplayPanel.tcl
#Language:  Tcl
#Date:      21.05.03
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for selecting vertices to be displayed
#
#************************************************************************


#--Select vertices to be displayed
#
proc VertexDisplay::openPanel {} {
  global VertexDisplay Info ObjectTable Model
  upvar #0 VertexDisplay theArray

  set w $VertexDisplay(winName)
  set wgeom $VertexDisplay(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set VertexDisplay(winId) $id

  set Info(thisWindow) $w
  set this $w

  if { 1 == [Util::checkPanelWindow VertexDisplay $id "Select vertices" $wgeom] } {
    raise $VertexDisplay(winName)
    focus -force $VertexDisplay(winName)
    return
  }  

  set VertexDisplay(dataChanged) 0
  set VertexDisplay(dataModified) 0

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $VertexDisplay(winTitle)
  wm geometry $w $wgeom

  # Pick vertex object ids
  set VertexDisplay(ids) [Object::getIds "V"]

  # Current display status
  foreach id $VertexDisplay(ids) {
    set VertexDisplay($id) [Object::getDisplayed $id]
    set VertexDisplay($id,old) [Object::getDisplayed $id]
    set VertexDisplay($id,slctd) [Object::getSelected $id]
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
  set f11 [frame $w.f1.f1] ;#-check box + vertex name
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

  set VertexDisplay(vertexLB) $wdg
  
  VertexDisplay::update

  pack $f11.sx -side bottom -fill x -expand 0 
  pack $f11.lb -side left -fill both -expand 1
  pack $f11.sy -side left -fill y -expand 0

  set VertexDisplay(vertexLB) $wdg

  bind $wdg <ButtonRelease-1> "VertexDisplay::setDisplayMode {ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <ButtonRelease-2> "VertexDisplay::setDisplayMode {ButtonRelease-2 %x %y %X %Y}"
  bind $wdg <ButtonRelease-3> "VertexDisplay::setDisplayMode {ButtonRelease-3 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-1> "VertexDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-2> "VertexDisplay::setDisplayMode {Control-ButtonRelease-2 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-3> "VertexDisplay::setDisplayMode {Control-ButtonRelease-3 %x %y %X %Y}"
  bind $wdg <KeyPress-Return> "VertexDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <KeyPress-space> "VertexDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"

  bind $wdg <Double-1> "VertexDisplay::select"
  bind $wdg <Control-Double-1> "VertexDisplay::select"
  bind $wdg <Control-KeyPress-Return> "VertexDisplay::select"
  bind $wdg <Control-KeyPress-space> "VertexDisplay::select"

  #-All, None buttons
  set all_btn [button $f12.all -text All -command "VertexDisplay::all" \
                            -width 5]
  set none_btn [button $f12.none -text None -command "VertexDisplay::none" \
                            -width 5]

  set select_btn [button $f12.select -text Select -command "VertexDisplay::select" \
                            -width 5]

  pack $all_btn $none_btn -side top -expand 0 -pady $fpy3
  pack $select_btn -side top -expand 0 -pady [expr 2 * $fpy3]

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  #-Buttons
  set ok_btn [button $f2.ok -text OK -command "VertexDisplay::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "VertexDisplay::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f2.apply -text Apply -command "VertexDisplay::panelApply" \
                               -state $ap]
  
  focus $ok_btn
  set VertexDisplay(applyButton) $ap_btn
  set VertexDisplay(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 0 -padx $fpx2
  #pack $ok_btn $cn_btn -side left -expand 1 -padx $fpx1
  
  pack $f1 -side top -anchor w -expand 1 -fill both -padx $fpx2 -pady $fpy2
  pack $f11 -side left -expand 1 -fill both -padx $fpx2
  pack $f12 -side left -expand 0 -padx $fpx2
  pack $f2 -side top -expand 0 -padx $fpx2 -pady $fpy2
}


proc VertexDisplay::formVertexLBRow {id} {
  global VertexDisplay Info ObjectTable

  if { $VertexDisplay($id) } {
    set row $Info(selectionBoxTrueMarker)

  } else {
    set row $Info(selectionBoxFalseMarker)
  }

  append row "("
  append row $ObjectTable($id,tg)
  append row ")"

  append row $ObjectTable($id,nm)

  if { $VertexDisplay($id,slctd) } {

    append row " \[#\] "
  }

  return $row
}


proc VertexDisplay::all {} {
  global VertexDisplay ObjectTable

  set changed 0

  set index 0

  foreach id $VertexDisplay(ids) {
    
    if { $VertexDisplay($id) != 1 } {
      set changed 1

      set VertexDisplay($id) 1

      set ObjectTable($id,dspl) $VertexDisplay($id)

      set row [VertexDisplay::formVertexLBRow $id]
      ListBox::updateRow $VertexDisplay(vertexLB) $index $row
    }

    incr index
  }

  $VertexDisplay(vertexLB) selection clear 0 end

  if {$changed} {
    Panel::panelDataChanged 1 VertexDisplay
    VertexDisplay::updateDisplay
  }
}


proc VertexDisplay::none {} {
  global VertexDisplay ObjectTable

  set changed 0

  set index 0

  foreach id $VertexDisplay(ids) {
    
    if { $VertexDisplay($id) != 0 } {

      set changed 1

      set VertexDisplay($id) 0

      set ObjectTable($id,dspl) $VertexDisplay($id)

      set row [VertexDisplay::formVertexLBRow $id]

      ListBox::updateRow $VertexDisplay(vertexLB) $index $row

    }

    incr index
  }

  $VertexDisplay(vertexLB) selection clear 0 end

  if {$changed} {
    Panel::panelDataChanged 1 VertexDisplay
    VertexDisplay::updateDisplay
  }
}


proc VertexDisplay::select {} {
  global Info ObjectTable VertexDisplay

  set sindex [$VertexDisplay(vertexLB) curselection]

  if { $sindex == "" || $sindex < 0 } {
    return
  }

  set id [lindex $VertexDisplay(ids) $sindex]
  set bd1 -1
  set lr1 -1
  set bd2 -1
  set lr2 -1

  set bndr_id [lindex $ObjectTable($id,prIds) 0]
  if { $bndr_id != "" } {
    set bndr_pr_id $ObjectTable($bndr_id,prId)

    if { $bndr_pr_id > $Info(NO_INDEX) } {
      set bd1 $ObjectTable($bndr_pr_id,pr1Id)
      set bd2 $ObjectTable($bndr_pr_id,pr2Id)
    }
  }

  set accept_body_change 1
  set update_gui 1

  #--NOTE: Apply vertex selection via cpp !!!
  #
  set data "$id $bd1 $lr1 $bd2 $lr2 $id $accept_body_change $update_gui"
  Util::cpp_exec boundarySelected $data

  $VertexDisplay(vertexLB) selection set $sindex
}


#
proc VertexDisplay::setDisplayMode {einfo} {
  global VertexDisplay ObjectTable

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
    set sindex [$VertexDisplay(vertexLB) index @$x,$y]
    set windex [$VertexDisplay(vertexLB) index @$x,$y]
    set do_select 1

  #-Keep currently selected row, manipulate currently selected
  } elseif { $key == "Control-ButtonRelease-1"
     } {
    set sindex [$VertexDisplay(vertexLB) curselection]
    set windex [$VertexDisplay(vertexLB) curselection]
    set do_select 0

  #-Keep currently selected row, manipulate clicked
  } elseif { $key == "Control-ButtonRelease-2" ||
             $key == "Control-ButtonRelease-3" 
     } {
    set sindex [$VertexDisplay(vertexLB) curselection]
    set windex [$VertexDisplay(vertexLB) index @$x,$y]
    set do_select 0
  }

  if { $windex == "" || $windex < 0 } {
    return
  }

  if { $do_select } {
    $VertexDisplay(vertexLB) selection clear 0 end
  }

  set id [lindex $VertexDisplay(ids) $windex]

  set dmode $VertexDisplay($id)

  # Toggle selection
  if { $dmode == 1 } {
    set VertexDisplay($id) 0
  } else {
    set VertexDisplay($id) 1
  }

  set ObjectTable($id,dspl) $VertexDisplay($id)

  set row [VertexDisplay::formVertexLBRow $id]

  ListBox::updateRow $VertexDisplay(vertexLB) $windex $row 0

  # NOTE: When listbox has focus, then Return/Space seem to
  # advance the selected row, this fixes this!
  #
  if { $sindex != "" && $sindex >= 0 } {
    $VertexDisplay(vertexLB) selection set $sindex
    $VertexDisplay(vertexLB) activate $sindex
  }

  Panel::panelDataChanged 1 VertexDisplay

  VertexDisplay::updateDisplay
}


proc VertexDisplay::updateDisplay {} {
  global VertexDisplay ObjectTable

  # Refresh vertex display
  Util::updateGui
  Util::cpp_exec vertexDisplayPanelOk
}


proc VertexDisplay::panelOk {w} {
  global VertexDisplay

  VertexDisplay::panelApply
  
  Panel::cancel $w
}


proc VertexDisplay::panelApply {} {
  global VertexDisplay ObjectTable

  #---No changes
  if { !$VertexDisplay(dataChanged) } {
    return
  }

  if { -1 == [VertexDisplay::panelCheck] } {
    return
  }
  
  foreach id $VertexDisplay(ids) {
    set ObjectTable($id,dspl) $VertexDisplay($id)
    set VertexDisplay($id,old) $VertexDisplay($id)
  }

  Panel::panelDataChanged 0 VertexDisplay 
  Panel::panelDataModified 0 VertexDisplay 

  VertexDisplay::updateDisplay
}


proc VertexDisplay::panelCancel {w} {
  global VertexDisplay ObjectTable

  if {$VertexDisplay(dataChanged)} {

    foreach id $VertexDisplay(ids) {
      set VertexDisplay($id) $VertexDisplay($id,old)
    }

    # Refresh display
    VertexDisplay::panelApply
  }

  Panel::cancel $w
}


proc VertexDisplay::panelCheck {} {
  global VertexDisplay

  #-Ok
  return 1
}


proc VertexDisplay::update {} {
  global VertexDisplay

  set start_pos [lindex [$VertexDisplay(vertexLB) yview] 0]

  set vertex_list ""

  foreach id $VertexDisplay(ids) {
    set row [VertexDisplay::formVertexLBRow $id]
    lappend vertex_list $row
  }

  ListBox::fill $VertexDisplay(vertexLB) $vertex_list

  $VertexDisplay(vertexLB) yview moveto $start_pos
}


# NOTE: This is called from cpp via
# Interface::setVertexSelectionMode
# when a geometry vertex is selected
#
proc VertexDisplay::setSelectionMode {id mode do_update} {
  global VertexDisplay

  set idx [lsearch $VertexDisplay(ids) $id]

  if { -1 == $idx } return

  set VertexDisplay($id,slctd) $mode

  if { !$do_update } return

  VertexDisplay::update

  # Check that lb-row is visible (relevant when
  # selected from geometry!)
  #
  set y_pos [$VertexDisplay(vertexLB) yview]

  set sz [expr 1.0 * [$VertexDisplay(vertexLB) size] - 1]
  set cur_pos [expr $idx / $sz]

  set y_start [lindex $y_pos 0]
  set y_end [lindex $y_pos 1]

  if { $cur_pos < $y_start || $cur_pos > $y_end } {
    $VertexDisplay(vertexLB) see $idx
  }
}


#ecif_VertexDisplayPanel.tcl

