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
#Module:    ecif_tk_bodyDisplayPanel.tcl
#Language:  Tcl
#Date:      13.01.03
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for selecting bodies to be displayed
#
#************************************************************************


#--Select bodies to be displayed
#
proc BodyDisplay::openPanel {} {
  global BodyDisplay Info ObjectTable Model
  upvar #0 BodyDisplay theArray

  set w $BodyDisplay(winName)
  set wgeom $BodyDisplay(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set BodyDisplay(winId) $id

  set Info(thisWindow) $w
  set this $w

  if { 1 == [Util::checkPanelWindow BodyDisplay $id "Select bodies" $wgeom] } {
    raise $BodyDisplay(winName)
    focus -force $BodyDisplay(winName)
    return
  }  

  set BodyDisplay(dataChanged) 0
  set BodyDisplay(dataModified) 0

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $BodyDisplay(winTitle)
  wm geometry $w $wgeom

  # Pick body object ids
  set BodyDisplay(ids) [Object::getIds "B"]

  # Current display status
  foreach id $BodyDisplay(ids) {
    set BodyDisplay($id) [Object::getDisplayed $id]
    set BodyDisplay($id,old) [Object::getDisplayed $id]
    #set BodyDisplay($id,slctd) [Object::getSelected $id]
    set BodyDisplay($id,slctd) 0
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
  set f11 [frame $w.f1.f1] ;#-check box + body name
  set f12 [frame $w.f1.f2] ;#-All None button area
  set f2 [frame $w.f2] ;    #Apply+Ok+cancel buttons frame


  #-Bodies listbox
  set wdg [ listbox $f11.lb -relief sunken \
             -selectmode browse -exportselection 0 \
             -height 20 -width 40 -font $Info(tableFont) \
             -xscrollcommand [list $f11.sx set] \
             -yscrollcommand [list $f11.sy set] ]

  scrollbar $f11.sx -orient horizontal -command [list $wdg xview]
  scrollbar $f11.sy -orient vertical -command [list $wdg yview]

  set BodyDisplay(bodyLB) $wdg 
  
  BodyDisplay::update

  BodyDisplay::setListColors $wdg

  pack $f11.sx -side bottom -fill x -expand 0 
  pack $f11.lb -side left -fill both -expand 1
  pack $f11.sy -side left -fill y -expand 0

  set BodyDisplay(bodyLB) $wdg

  bind $wdg <ButtonRelease-1> "BodyDisplay::setDisplayMode {ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <ButtonRelease-2> "BodyDisplay::setDisplayMode {ButtonRelease-2 %x %y %X %Y}"
  bind $wdg <ButtonRelease-3> "BodyDisplay::setDisplayMode {ButtonRelease-3 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-1> "BodyDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-2> "BodyDisplay::setDisplayMode {Control-ButtonRelease-2 %x %y %X %Y}"
  bind $wdg <Control-ButtonRelease-3> "BodyDisplay::setDisplayMode {Control-ButtonRelease-3 %x %y %X %Y}"
  bind $wdg <KeyPress-Return> "BodyDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"
  bind $wdg <KeyPress-space> "BodyDisplay::setDisplayMode {Control-ButtonRelease-1 %x %y %X %Y}"

  bind $wdg <Double-1> "BodyDisplay::select"
  bind $wdg <Control-Double-1> "BodyDisplay::select"
  bind $wdg <Control-KeyPress-Return> "BodyDisplay::select"
  bind $wdg <Control-KeyPress-space> "BodyDisplay::select"

  #-All, None buttons
  set all_btn [button $f12.all -text All -command "BodyDisplay::all" \
                            -width 3]
  set none_btn [button $f12.none -text None -command "BodyDisplay::none" \
                            -width 3]

  set select_btn [button $f12.select -text Select -command "BodyDisplay::select" \
                            -width 5]

  pack $all_btn $none_btn -side top -expand 0 -pady $fpy3
  pack $select_btn -side top -expand 0 -pady [expr 2 * $fpy3]

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  #-Buttons
  set ok_btn [button $f2.ok -text OK -command "BodyDisplay::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "BodyDisplay::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f2.apply -text Apply -command "BodyDisplay::panelApply" \
                               -state $ap]
  
  focus $ok_btn
  set BodyDisplay(applyButton) $ap_btn
  set BodyDisplay(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 0 -padx $fpx2
  #pack $ok_btn $cn_btn -side left -expand 1 -padx $fpx1
  
  pack $f1 -side top -anchor w -expand 1 -fill both -padx $fpx2 -pady $fpy2
  pack $f11 -side left -expand 1 -fill both -padx $fpx2
  pack $f12 -side left -expand 0 -padx $fpx2
  pack $f2 -side top -expand 0 -padx $fpx2 -pady $fpy2
}


proc BodyDisplay::formBodyLBRow {id} {
  global BodyDisplay Info ObjectTable

  if { $BodyDisplay($id) } {
    set row $Info(selectionBoxTrueMarker)

  } else {
    set row $Info(selectionBoxFalseMarker)
  }

  append row "("
  append row $ObjectTable($id,tg)
  append row ")"
  append row $ObjectTable($id,nm)

  if { $BodyDisplay($id,slctd) } {
    append row " \[#\] "
  }

  return $row
}


proc BodyDisplay::all {} {
  global BodyDisplay ObjectTable

  set changed 0

  set index 0

  foreach id $BodyDisplay(ids) {
    
    if { $BodyDisplay($id) != 1 } {
      set changed 1

      set BodyDisplay($id) 1

      set ObjectTable($id,dspl) $BodyDisplay($id)

      set row [BodyDisplay::formBodyLBRow $id]
      ListBox::updateRow $BodyDisplay(bodyLB) $index $row
    }

    incr index
  }

  BodyDisplay::setListColors $BodyDisplay(bodyLB)

  $BodyDisplay(bodyLB) selection clear 0 end

  if {$changed} {
    Panel::panelDataChanged 1 BodyDisplay
    BodyDisplay::updateDisplay
  }
}


proc BodyDisplay::none {} {
  global BodyDisplay ObjectTable

  set changed 0

  set index 0

  foreach id $BodyDisplay(ids) {
    
    if { $BodyDisplay($id) != 0 } {

      set changed 1

      set BodyDisplay($id) 0

      set ObjectTable($id,dspl) $BodyDisplay($id)

      set row [BodyDisplay::formBodyLBRow $id]

      ListBox::updateRow $BodyDisplay(bodyLB) $index $row

    }

    incr index
  }

  BodyDisplay::setListColors $BodyDisplay(bodyLB)

  $BodyDisplay(bodyLB) selection clear 0 end

  if {$changed} {
    Panel::panelDataChanged 1 BodyDisplay
    BodyDisplay::updateDisplay
  }
}


proc BodyDisplay::select {} {
  global BodyDisplay ObjectTable

  set sindex [$BodyDisplay(bodyLB) curselection]

  if { $sindex == "" || $sindex < 0 } {
    return
  }

  set bd_id [lindex $BodyDisplay(ids) $sindex]
  set lr_id -1

  #--NOTE: Apply body selection via cpp !!!
  #
  set data "$bd_id $lr_id"
  Util::cpp_exec bodySelected $data

  $BodyDisplay(bodyLB) selection set $sindex
}


#
proc BodyDisplay::setDisplayMode {einfo} {
  global BodyDisplay ObjectTable

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
    set sindex [$BodyDisplay(bodyLB) index @$x,$y]
    set windex [$BodyDisplay(bodyLB) index @$x,$y]
    set do_select 1

  #-Keep currently selected row, manipulate currently selected
  } elseif { $key == "Control-ButtonRelease-1"
     } {
    set sindex [$BodyDisplay(bodyLB) curselection]
    set windex [$BodyDisplay(bodyLB) curselection]
    set do_select 0

  #-Keep currently selected row, manipulate clicked
  } elseif { $key == "Control-ButtonRelease-2" ||
             $key == "Control-ButtonRelease-3" 
     } {
    set sindex [$BodyDisplay(bodyLB) curselection]
    set windex [$BodyDisplay(bodyLB) index @$x,$y]
    set do_select 0
  }

  if { $windex == "" || $windex < 0 } {
    return
  }

  if { $do_select } {
    $BodyDisplay(bodyLB) selection clear 0 end
  }

  set id [lindex $BodyDisplay(ids) $windex]

  set dmode $BodyDisplay($id)

  # Toggle selection
  if { $dmode == 1 } {
    set BodyDisplay($id) 0
  } else {
    set BodyDisplay($id) 1
  }

  set ObjectTable($id,dspl) $BodyDisplay($id)

  set row [BodyDisplay::formBodyLBRow $id]

  ListBox::updateRow $BodyDisplay(bodyLB) $windex $row 0

  if { $do_select } {
    BodyDisplay::setListColors $BodyDisplay(bodyLB) $sindex
  }

  # NOTE: When listbox has focus, then Return/Space seem to
  # advance the selected row, this fixes this!
  #
  if { $sindex != "" && $sindex >= 0 } {
    $BodyDisplay(bodyLB) selection set $sindex
    $BodyDisplay(bodyLB) activate $sindex
  }

  Panel::panelDataChanged 1 BodyDisplay

  BodyDisplay::updateDisplay
}


proc BodyDisplay::updateDisplay {} {
  global BodyDisplay ObjectTable

  # Refresh body display
  Util::updateGui
  Util::cpp_exec bodyDisplayPanelOk
}


proc BodyDisplay::panelOk {w} {
  global BodyDisplay

  BodyDisplay::panelApply
  
  Panel::cancel $w
}


proc BodyDisplay::panelApply {} {
  global BodyDisplay ObjectTable

  #---No changes
  if { !$BodyDisplay(dataChanged) } {
    return
  }

  if { -1 == [BodyDisplay::panelCheck] } {
    return
  }
  
  foreach id $BodyDisplay(ids) {
    set ObjectTable($id,dspl) $BodyDisplay($id)
    set BodyDisplay($id,old) $BodyDisplay($id)
  }

  Panel::panelDataChanged 0 BodyDisplay 
  Panel::panelDataModified 0 BodyDisplay 

  BodyDisplay::updateDisplay
}


proc BodyDisplay::panelCancel {w} {
  global BodyDisplay ObjectTable

  if {$BodyDisplay(dataChanged)} {

    foreach id $BodyDisplay(ids) {
      set BodyDisplay($id) $BodyDisplay($id,old)
    }

    # Refresh display
    BodyDisplay::panelApply
  }

  Panel::cancel $w
}


proc BodyDisplay::panelCheck {} {
  global BodyDisplay

  #-Ok
  return 1
}

# Set listbox row colors by body colors
#
proc BodyDisplay::setListColors { wdg  {idx ""} } {
  global BodyDisplay ObjectTable

  set idx 0
  foreach id $BodyDisplay(ids) {
    $wdg itemconfigure $idx -selectbackground $ObjectTable($id,clr)
    incr idx
  }
  return

  if { $idx != "" } {
    set id [lindex $BodyDisplay(ids) $idx]
    $wdg itemconfigure $idx -bg $ObjectTable($id,clr)
    $wdg itemconfigure $idx -fg $ObjectTable($id,clr)
    $wdg itemconfigure $idx -selectforeground $ObjectTable($id,clr)
  }
  return

  if { $idx != "" } {
    set id [lindex $BodyDisplay(ids) $idx]
    $wdg itemconfigure $idx -bg $ObjectTable($id,clr)
    $wdg itemconfigure $idx -fg white

  } else {
    set idx 0
    foreach id $BodyDisplay(ids) {
      $wdg itemconfigure $idx -bg $ObjectTable($id,clr)
      $wdg itemconfigure $idx -fg white
      incr idx
    }
  }

}


proc BodyDisplay::update {} {
  global BodyDisplay

  set body_list ""

  foreach id $BodyDisplay(ids) {
    set row [BodyDisplay::formBodyLBRow $id]
    lappend body_list $row
  }

  ListBox::fill $BodyDisplay(bodyLB) $body_list
}


proc BodyDisplay::setSelectionMode {id } {
  global BodyDisplay
  
  # Not working correctly, so not in use 
  # currently! (MVe 09.05.03)
  return

  set idx [lsearch $BodyDisplay(ids) $id]

  if { -1 == $idx } return
  
  # NOTE: The selection stat is simply toggled, because
  # model does not keep any selection mode info for bodies
  # We are just reacting to mouse-clicks here and toggle the
  # selection marker in this panel
  #
  # Toggle selection
  if { $BodyDisplay($id,slctd) } {
    set BodyDisplay($id,slctd) 0
  } else {
    set BodyDisplay($id,slctd) 1
  }

  BodyDisplay::update

  $BodyDisplay(bodyLB) yview $idx
}


#ecif_bodyDisplayPanel.tcl

