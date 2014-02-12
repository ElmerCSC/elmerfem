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
#Module:    ecif_tk_meshSelectPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for selecting a mesh for the model
#
#************************************************************************



proc MeshSelect::openPanel {} {
  # Global variables
  global Info MeshSelect Model ModelProperty ModelInfo

  set w .selectMesh
  set wtitle "Select mesh"
  set wgeom +420+120

  toplevel $w
  focus $w
  set this $w

  wm title $w $wtitle
  wm geometry $w $wgeom 

  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  set MeshSelect(meshName) ""
  set MeshSelect(meshH) ""
  set MeshSelect(nofNodes) ""
  set MeshSelect(nofElements) ""

  # Frames
  frame $w.f1  ;# Mesh names frame
  frame $w.f2  ;# Mesh info frame
  frame $w.fB  ;# Buttons

  # Create and pack mesh names frame
  # ================================
  set m $w.f1

  set bg  $Info(nonActiveBg)
  
  set lb_frame [frame $m.lbf]

  set lb [listbox $lb_frame.lb  -height 8 -width 35 \
          -xscrollcommand [list $lb_frame.sb_x set]  \
          -yscrollcommand [list $lb_frame.sb_y set]  \
          -bg $bg ]
  set sb_x [scrollbar $lb_frame.sb_x -orient horizontal \
                                     -command [list $lb_frame.lb xview] ]
  set sb_y [scrollbar $lb_frame.sb_y -orient vertical \
                                     -command [list $lb_frame.lb yview] ]

  set MeshSelect(meshNamesLB) $lb

  ListBox::fill $MeshSelect(meshNamesLB) $Model(meshNames)

  bind $MeshSelect(meshNamesLB) <ButtonRelease-1> "+MeshSelect::displayMeshInfo"

  pack $lb_frame -expand 1 -fill both
  pack $sb_x -side bottom -fill x -expand 0
  pack $lb -side left -fill both -expand 1 
  pack $sb_y -side left -fill y -expand 0 

  
  # Create and pack mesh info frame
  # ===============================
  
  set lwid  15
  set ewid1 32
  set ewid2 15

  # Mesh name
  # ---------
  set m [frame $w.f2.name]

  #-Label
  frame $m.lf
  label $m.lf.l -text "Mesh name:" -width $lwid
  pack $m.lf.l -side top -anchor w

  #-Entry
  frame $m.ef
  entry $m.ef.e  -textvariable MeshSelect(meshName) -width $ewid1 \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -side top -anchor w -padx $fpx1 -pady $fpy3 -expand 1

  pack $m.lf $m.ef -side left -anchor w -padx $fpx1 -pady $fpy1 -expand 1
  pack $m -side top -anchor w


  # Mesh-h
  # ------
  set m [frame $w.f2.mesh]

  #-Label
  frame $m.lf
  label $m.lf.l -text "Mesh parameter:" -width $lwid
  pack $m.lf.l -side top -anchor w

  #-Entry
  frame $m.ef
  entry $m.ef.e  -textvariable MeshSelect(meshH) -width $ewid2 \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -side top -anchor w -padx $fpx1 -pady $fpy3 -expand 1

  pack $m.lf $m.ef -side left -anchor w -padx $fpx1 -pady $fpy1 -expand 1
  pack $m -side top -anchor w


  # Nof nodes
  # ---------
  set m [frame $w.f2.nodes]

  #-Label
  frame $m.lf
  label $m.lf.l -text "Nof nodes:" -width $lwid
  pack $m.lf.l -side top -anchor w

  #-Entry
  frame $m.ef
  entry $m.ef.e  -textvariable MeshSelect(nofNodes) -width $ewid2 \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -side top -anchor w -padx $fpx1 -pady $fpy3 -expand 1

  pack $m.lf $m.ef -side left -anchor w -padx $fpx1 -pady $fpy1 -expand 1
  pack $m -side top -anchor w


  # Nof elements
  # ------------
  set m [frame $w.f2.elems]

  #-Label
  frame $m.lf
  label $m.lf.l -text "Nof elements:" -width $lwid
  pack $m.lf.l -side top -anchor w

  #-Entry
  frame $m.ef
  entry $m.ef.e  -textvariable MeshSelect(nofElements) -width $ewid2 \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -side top -anchor w -padx $fpx1 -pady $fpy1 -expand 1

  pack $m.lf $m.ef -side left -anchor w -padx $fpx1 -pady $fpy3 -expand 1
  pack $m -side top -anchor w


  # Create and pack buttons
  # =======================
  set m $w.fB

  #-Ok button

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $m.ok -text OK -command "MeshSelect::panelOk $this"]
  set cn_btn [button $m.cancel -text Cancel -command "MeshSelect::panelCancel $this" \
                               -state $ca]
  set ap_btn [button $m.apply -text Apply -command "MeshSelect::panelApply $this" \
                              -state $ap]

  set MeshSelect(applyButton) $ap_btn
  set MeshSelect(cancelButton) $cn_btn

  focus $ok_btn
  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1

  pack $w.f1 $w.f2 $w.fB -side top -expand 1 -padx $fpx1 -pady $fpy3

  # Set current mesh selected
  set index [lsearch $Model(meshNames) $ModelProperty(CURRENT_MESH_NAME)]

  if { $index != -1 } {
    $MeshSelect(meshNamesLB) selection set $index
  } else {
    $MeshSelect(meshNamesLB) selection set 0
  }

  MeshSelect::displayMeshInfo
}


proc MeshSelect::panelOk {w {do_cancel 1}} {
  global MeshSelect Model ModelProperty UserSetting

  # Load new mesh if we have some and if name changed
  # or no mesh yet!
  if { [llength $Model(meshNames)] > 0 } {
  
    if { $ModelProperty(CURRENT_MESH_NAME) != $MeshSelect(meshName) ||
         !$Model(Mesh,exists)
       } {
  
      # Model file needs update!
      set Model(Front,needsUpdate) 1

      set ModelProperty(CURRENT_MESH_NAME) $MeshSelect(meshName)

      ModelProperty::setCurrentMeshDirectory
      Util::cpp_exec modelPropertyPanelOk

      MenuExec::unloadMesh

      if {$UserSetting(AUTO_LOAD_MESH)} {
        MenuExec::loadMesh
      }
    }
  }

  # If we should close the panel
  if { $do_cancel } {
    Panel::cancel $w
  }
}


proc MeshSelect::panelApply {w} {

  MeshSelect::panelOk $w 0
}


proc MeshSelect::panelCancel {w} {

  Panel::cancel $w
}


proc MeshSelect::displayMeshInfo {} {
  global MeshSelect Info Model ModelProperty MeshSelect

  set index [$MeshSelect(meshNamesLB) curselection]

  if { $index == "" } {
    return
  }

  set name [lindex $Model(meshNames) $index]

  # Name changed fromt what is current in the model!
  if { $name != $ModelProperty(CURRENT_MESH_NAME) } {
    Panel::panelDataChanged 1 MeshSelect
  } else {
    Panel::panelDataChanged 0 MeshSelect
  }

  set MeshSelect(meshName) $name

  # Read mesh-h from mesh info
  # ==========================
  set MeshSelect(meshH) ""

  set fname1 [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                        $Info(meshDirectoryName) \
                        $name \
                        $Info(meshInfoFileName)]

  if { ![catch {set ch1 [open $fname1]}] } {

    while { ![eof $ch1] } {

      set tmp [split [gets $ch1] "="]

      # Check if keyword found
      if { "h-value" == [lindex $tmp 0] } {
        set MeshSelect(meshH) [string trim [lindex $tmp 1]]
        break
      }
    }
    
    catch { close $ch1 }
  }


  # Read nof nodes etc. from mesh header
  # ====================================

  set MeshSelect(nofNodes) ""
  set MeshSelect(nofElements) ""

  set fname2 [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                        $Info(meshDirectoryName) \
                        $name \
                        "mesh.header"]

  if { ![catch {set ch2 [open $fname2]}] } {
    set info2 [split [gets $ch2]]

    catch { close $ch2 }

    set MeshSelect(nofNodes) [lindex $info2 0]
    set MeshSelect(nofElements) [lindex $info2 1]
  }

}


#End ecif_meshSelectPanel.tcl

