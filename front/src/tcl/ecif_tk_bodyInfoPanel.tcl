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
#Module:    ecif_tk_bodyInfoPanel.tcl
#Language:  Tcl
#Date:      05.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for general body information.
#
#************************************************************************


#------Bodies info panel proc------
#
# This procedure displays the bodies info panel
#
proc BodyInfo::openPanel {} {
   global BodyInfo Info ObjectTable

  set w $BodyInfo(winName)
  set wgeom $BodyInfo(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set BodyInfo(infoWinId) $id

  #--Store last created window info in Info
  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow BodyInfo $id $BodyInfo(winTitle) $wgeom] } {
    return
  }  

  toplevel $w 
  focus $w 

  #--Window properties
  wm title $w $BodyInfo(winTitle)
  wm geometry $w $wgeom

  set this $w

  set BodyInfo(ids) ""

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    lappend BodyInfo(ids) $id
  }

  #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  set bg  $Info(nonActiveBgLight)

  #---Outer frames
  set f1 [frame $w.f1] ;#List box and colr/name entries common frame
  set f11 [frame $w.f1.f1] ;#List box frame
  set f12 [frame $w.f1.f2] ;#Color/name entries frame
  set f2 [frame $w.f2] ;#Other entries common frame
  set f3 [frame $w.f3] ;#Apply+Ok+cancel buttons frame
  #
  set lb [listbox $f11.lb -selectmode browse -height 5 -width 20 \
      -yscrollcommand [list $f11.sb_y set]]
  set sb_y [scrollbar $f11.sb_y -orient vertical  -command [list $f11.lb yview] ]

  set names ""
  foreach id $BodyInfo(ids) {
    lappend names $ObjectTable($id,nm)
  }

  ListBox::fill $lb $names

  pack $lb $sb_y -side left -fill y

  set BodyInfo(listBox) $lb
  bind $lb <ButtonRelease-1> BodyInfo::setPanelContents 
  
  # Body color entry
  set ce [entry $f12.cb -width 5 -state disabled ]
  set BodyInfo(bodyColorEntry) $ce

  # Name entry
  set ne [entry $f12.ne -width 20 -textvariable BodyInfo(name) -state disabled -bg $bg]
  pack $ce $ne -side left -padx $fpx2 -pady $fpy2

  # Set first body selected
  $lb selection set 0 0
  BodyInfo::setPanelContents

  # Other fields
  # Frame
  set fX [frame $f2.fX]
  set fY [frame $f2.fY]
  set fZ [frame $f2.fZ]
  set fME [frame $f2.fME]

  # Labels
  set X_lbl  [label $fX.l -width 25 -text "Body X dimensions  \[m\]: " ]
  set Y_lbl  [label $fY.l -width 25 -text "Body Y dimensions  \[m\]: " ]
  set Z_lbl  [label $fZ.l -width 25 -text "Body Z dimensions  \[m\]: " ]
  set ME_lbl [label $fME.l -width 25 -text "Nof mesh elements: " ]

  # Entries
  set minX_e [entry $fX.minX -width 16 -textvariable BodyInfo(minX) -state disabled -bg $bg]
  set minY_e [entry $fY.minY -width 16 -textvariable BodyInfo(minY) -state disabled -bg $bg]
  set minZ_e [entry $fZ.minZ -width 16 -textvariable BodyInfo(minZ) -state disabled -bg $bg]

  set maxX_e [entry $fX.maxX -width 16 -textvariable BodyInfo(maxX) -state disabled -bg $bg]
  set maxY_e [entry $fY.maxY -width 16 -textvariable BodyInfo(maxY) -state disabled -bg $bg]
  set maxZ_e [entry $fZ.maxZ -width 16 -textvariable BodyInfo(maxZ) -state disabled -bg $bg]

  set sizeX_e [entry $fX.sizeX -width 16 -textvariable BodyInfo(sizeX) -state disabled -bg $bg]
  set sizeY_e [entry $fY.sizeY -width 16 -textvariable BodyInfo(sizeY) -state disabled -bg $bg]
  set sizeZ_e [entry $fZ.sizeZ -width 16 -textvariable BodyInfo(sizeZ) -state disabled -bg $bg]

  set nofME_e [entry $fME.nofME -width 16 -textvariable BodyInfo(nofMeshElements) -state disabled -bg $bg]

  pack $X_lbl $minX_e $maxX_e $sizeX_e -side left -padx $fpx2 
  pack $Y_lbl $minY_e $maxY_e $sizeY_e -side left -padx $fpx2 
  pack $Z_lbl $minZ_e $maxZ_e $sizeZ_e -side left -padx $fpx2 

  pack $ME_lbl $nofME_e -side left -padx $fpx2 -pady $fpy2

  pack $fX $fY $fZ $fME -side top -pady $fpy2 -anchor w

  set ok_btn [button $f3.ok -text OK -command "BodyInfo::panelOk $this"]

  focus $ok_btn

  pack $ok_btn -side left -expand 1 -padx $fpx1
  
  # Pack frames
  #
  pack $f1 -side top  -anchor w -expand 1 -padx $fpx2 -pady $fpy2
  pack $f2 $f3 -side top -expand 1 -padx $fpx2 -pady $fpy2
  pack $f11 $f12 -side left -padx $fpx3
}


# Body info panel bind proc: reacts to list box selection
proc BodyInfo::setPanelContents {} {
  global BodyInfo ObjectTable

  set index [$BodyInfo(listBox) curselection]

  set id [lindex $BodyInfo(ids) $index ]

  # Name
  set BodyInfo(name) $ObjectTable($id,nm)

  # Color
  set BodyInfo(color) $ObjectTable($id,clr)

  $BodyInfo(bodyColorEntry) configure -bg $BodyInfo(color) 

  # Size
  set BodyInfo(minX)  $ObjectTable($id,mnX)
  set BodyInfo(minY)  $ObjectTable($id,mnY)
  set BodyInfo(minZ)  $ObjectTable($id,mnZ)

  set BodyInfo(maxX)  $ObjectTable($id,mxX)
  set BodyInfo(maxY)  $ObjectTable($id,mxY)
  set BodyInfo(maxZ)  $ObjectTable($id,mxZ)

  set BodyInfo(sizeX) [expr $BodyInfo(maxX) - $BodyInfo(minX)]
  set BodyInfo(sizeY) [expr $BodyInfo(maxY) - $BodyInfo(minY)]
  set BodyInfo(sizeZ) [expr $BodyInfo(maxZ) - $BodyInfo(minZ)]

  # Nof mesh elements
  set BodyInfo(nofMeshElements) $ObjectTable($id,nofMshElm)
}
 

# Body names panel bind proc: reacts to list box selection
proc BodyInfo::setBodyColorEntry {index color} {
  global BodyInfo

  #-Color entry background color
  $BodyInfo(bodyColorEntry) configure -bg $color 
}


proc BodyInfo::panelOk {w} {
  Panel::cancel $w
}

proc BodyInfo::panelCancel {w} {
  Panel::cancel $w
}

# end ecif_tk_bodyInfo.tcl
# ********************
