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
#Module:    ecif_tk_labelDisplayPanel.tcl
#Language:  Tcl
#Date:      10.02.00
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for selecting object labels to be displayed
#
#************************************************************************


#--Select labels to be displayed
#
proc LabelDisplay::openPanel {} {
  global Info LabelDisplay Model ModelFlags
  upvar #0 LabelDisplay theArray

  set w $LabelDisplay(winName)
  set wgeom $LabelDisplay(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set LabelDisplay(winId) $id

  set Info(thisWindow) $w
  set this $w

  if { 1 == [Util::checkPanelWindow LabelDisplay $id "Labels" $wgeom] } {
    raise $LabelDisplay(winName)
    focus -force $LabelDisplay(winName)
    return
  }  

  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    set LabelDisplay(targets) "FACE EDGE VERTEX"
    set LabelDisplay(texts) "Boundary Edge Vertex"

  } else {
    set LabelDisplay(targets) "EDGE VERTEX"
    set LabelDisplay(texts) "Boundary Vertex"
  }

  #---Store old values
  foreach trg $LabelDisplay(targets) {
    set ModelFlags(LABEL_DISPLAY_$trg,old) $ModelFlags(LABEL_DISPLAY_$trg)
  }

  set LabelDisplay(dataChanged) 0
  set LabelDisplay(dataModified) 0

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $LabelDisplay(winTitle)
  wm geometry $w $wgeom

  #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)

  #---Outer frames
  set f1 [frame $w.f1]     ;#outer frame
  set f11 [frame $w.f1.f1] ;#-check box + body name and color frame
  set f12 [frame $w.f1.f2] ;#-All None button area
  set f2 [frame $w.f2] ;    #Apply+Ok+cancel buttons frame
  
  #---Label checkbuttons
  foreach trg $LabelDisplay(targets) txt $LabelDisplay(texts) {
    
    set var_name "ModelFlags(LABEL_DISPLAY_$trg)"

    set f [frame $f11.f$trg]
    set cb [checkbutton $f.cb$id -indicatoron 1 -variable $var_name -text $txt]
    $cb configure -command "LabelDisplay::setDisplayMode"

    set state normal

    # For mesh only geometry, we do not have
    # boundary labels!
    if { !$ModelFlags(GEOMETRY_TYPE_CAD) } {

      if { $trg == "FACE" || $trg == "EDGE" } {
        set state disabled
      }
    }

    $cb configure -state $state

    bind $cb <ButtonRelease-1> "Panel::panelDataChanged 1 LabelDisplay $cb {%A %K}"

    pack $cb -side left -anchor w
    pack $f -side top -anchor w -expand 1
  }

  #-All, None buttons
  set all_btn [button $f12.all -text All -command "LabelDisplay::all" \
                            -width 3]
  set none_btn [button $f12.none -text None -command "LabelDisplay::none" \
                            -width 3]

  pack $all_btn $none_btn -side top -expand 1 -pady $fpy2

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  #-Buttons
  set ok_btn [button $f2.ok -text OK -command "LabelDisplay::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "LabelDisplay::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f2.apply -text Apply -command "LabelDisplay::panelApply" \
                               -state $ap]
  
  focus $ok_btn
  set LabelDisplay(applyButton) $ap_btn
  set LabelDisplay(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1
  
  pack $f1 -side top -anchor w -expand 1 -padx $fpx2 -pady $fpy2
  pack $f1.f1 -side left -expand 1 -padx $fpx2
  pack $f1.f2 -side left -expand 1 -padx $fpx2
  pack $f2 -side top -expand 1 -padx $fpx2 -pady $fpy2
}


proc LabelDisplay::all {} {
  global LabelDisplay ModelFlags

  set changed 0

  foreach trg $LabelDisplay(targets) {
    
    if { $ModelFlags(LABEL_DISPLAY_$trg) != 1 } {
      set changed 1
    }

    set ModelFlags(LABEL_DISPLAY_$trg) 1
  }

  if {$changed} {
    Panel::panelDataChanged 1 LabelDisplay
    LabelDisplay::setDisplayMode
  }
}


proc LabelDisplay::none {} {
  global LabelDisplay ModelFlags

  set changed 0

  foreach trg $LabelDisplay(targets) {
    
    if { $ModelFlags(LABEL_DISPLAY_$trg) != 0 } {
      set changed 1
    }

    set ModelFlags(LABEL_DISPLAY_$trg) 0
  }

  if {$changed} {
    Panel::panelDataChanged 1 LabelDisplay
    LabelDisplay::setDisplayMode
  }
}


proc LabelDisplay::setDisplayMode {} {
  global LabelDisplay

  # Refresh label display
  Util::updateGui

  foreach trg $LabelDisplay(targets) {  
    Util::setFlagValue LABEL_DISPLAY LABEL_DISPLAY_$trg
  }
}


proc LabelDisplay::panelOk {w} {
  global LabelDisplay

  LabelDisplay::panelApply

  Panel::cancel $w
}


proc LabelDisplay::panelApply {} {
  global LabelDisplay ModelFlags

  #---No changes
  if { !$LabelDisplay(dataChanged) } {
    return
  }

  if { -1 == [ LabelDisplay::panelCheck] } {
    return
  }

  Panel::panelDataChanged 0 LabelDisplay 
  Panel::panelDataModified 0 LabelDisplay 

  #---Store old values
  foreach trg $LabelDisplay(targets) {
    set ModelFlags(LABEL_DISPLAY_$trg,old) $ModelFlags(LABEL_DISPLAY_$trg)
  }

  LabelDisplay::setDisplayMode
}


proc  LabelDisplay::panelCancel {w} {
  global  LabelDisplay ModelFlags

  if {$LabelDisplay(dataChanged)} {

    # Restore previous values
    foreach trg $LabelDisplay(targets) {
      set ModelFlags(LABEL_DISPLAY_$trg) $ModelFlags(LABEL_DISPLAY_$trg,old)
    }

    # Refresh display
    LabelDisplay::panelApply
  }

  Panel::cancel $w
}


proc LabelDisplay::panelCheck {} {
  global LabelDisplay

  #-Ok
  return 1
}


#ecif_labelDisplayPanel.tcl

