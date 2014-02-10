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
#Module:    ecif_tk_inputFilenfoPanel.tcl
#Language:  Tcl
#Date:      23.03.03
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for matc,solver kwrds, definition etc. information.
#
#************************************************************************
 
#------Input file info proc------
#
proc InputFileInfo::openPanel { } {
  # Global variables
  global Info InputFileInfo

  set w $InputFileInfo(winName)
  set wgeom $InputFileInfo(winGeometry)

  set id [winfo atom $w]
  set InputFileInfo(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow InputFileInfo $id $InputFileInfo(winTitle) $wgeom] } {
    return
  }  
 
  toplevel $w
  focus $w
  set this $w

  wm title $w $InputFileInfo(winTitle)
  wm geometry $w $wgeom 

  set InputFileInfo(files) ""
  set InputFileInfo(currentFileInfo) ""

  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  # Info textbox
  # ------------
  frame $w.fT
  set m $w.fT
  label $m.l -anchor nw -text "Input files read"

  set wdg [listbox $m.lb -relief groove -height 15 -width 60 \
                        -bg white -selectborderwidth 0 \
                        -font $Info(entryFont) \
                        -xscrollcommand [list $m.sx set] \
                        -yscrollcommand [list $m.sy set] ]

  bind $wdg <ButtonRelease-1> "InputFileInfo::selectInfo"
  bind $wdg <Double-1> "InputFileInfo::browseFile"

  set InputFileInfo(infoLB) $wdg


  scrollbar $m.sx -orient horizontal \
    -command [list $m.text xview]
  scrollbar $m.sy -orient vertical  \
    -command [list $m.text yview]

  InputFileInfo::updateInfo

  pack $w.fT -expand 1 -fill both
  pack $m.l -side top -anchor w
  pack $m.sx -side bottom -anchor n -fill x -expand 0
  pack $m.lb -side left -anchor w -fill both -expand 1
  pack $m.sy -side left -anchor n -fill y  -expand 0


  # Buttons
  # -------
  frame $w.fB

  #-Ok button
  set ok_btn [button $w.fB.ok -text OK -command "InputFileInfo::panelOk $this"]
  set browse_btn [button $w.fB.browse -text Browse -command "InputFileInfo::browseFile" -state disabled]
  set edit_btn [button $w.fB.edit -text Edit -command "InputFileInfo::editFile" -state disabled]
  set load_btn [button $w.fB.load -text Load -command "InputFileInfo::loadFile" -state disabled]

  set InputFileInfo(browseButton) $browse_btn
  set InputFileInfo(editButton) $edit_btn
  set InputFileInfo(loadButton) $load_btn

  focus $ok_btn

  pack $w.fB -expand 0
  pack $ok_btn $browse_btn $edit_btn $load_btn -side left -expand 1 -padx $fpx1
}


proc InputFileInfo::updateInfo {} {
  global Info InputFileInfo

  $InputFileInfo(infoLB) delete 0 end
  set InputFileInfo(files) ""
  set idx -1

  # Definition (edf) files
  # ----------------------
  #
  # Header line
  incr idx 1
  $InputFileInfo(infoLB) insert end "          *** Definition files ***"
  $InputFileInfo(infoLB) itemconfigure $idx -bg $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectforeground black

  foreach row $Info(definitionFiles) {
    set fn [lindex $row 0]
    set mt [lindex $row 1]
    set mts [Util::timeClock2Str $mt]
    if { 0 && $mts != "" } {
      $InputFileInfo(infoLB) insert end "$fn ($mts)"
    } else {
      $InputFileInfo(infoLB) insert end "$fn"
    }
    incr idx
    lappend InputFileInfo(files) [list $idx $fn "DEF" $mt]

  }

  # Separator line (lf)
  incr idx 1
  $InputFileInfo(infoLB) insert end ""
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground white

  # SolverKeyword files
  # -------------------
  #
  # Header line
  incr idx 1
  $InputFileInfo(infoLB) insert end "          *** Solver keyword files ***"
  $InputFileInfo(infoLB) itemconfigure $idx -bg $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectforeground black

  foreach row $Info(solverKeywordFiles) {
    set fn [lindex $row 0]
    set mt [lindex $row 1]
    set mts [Util::timeClock2Str $mt]
    if { 0 && $mts != "" } {
      $InputFileInfo(infoLB) insert end "$fn  ($mts)"
    } else {
      $InputFileInfo(infoLB) insert end "$fn"
    }
    incr idx
    lappend InputFileInfo(files) [list $idx $fn "KWD" $mt]
  }

  # Separator line (lf)
  incr idx 1
  $InputFileInfo(infoLB) insert end ""
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground white

  # Matc files
  # ----------
  #
  incr idx 1
  $InputFileInfo(infoLB) insert end "          *** Matc files ***"
  $InputFileInfo(infoLB) itemconfigure $idx -bg $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectforeground black

  foreach row $Info(matcFiles) {
    set fn [lindex $row 0]
    set mt [lindex $row 1]
    set mts [Util::timeClock2Str $mt]
    if { 0 && $mts != "" } {
      $InputFileInfo(infoLB) insert end "$fn  ($mts)"
    } else {
      $InputFileInfo(infoLB) insert end "$fn"
    }
    incr idx
    lappend InputFileInfo(files) [list $idx $fn "MTC" $mt]
  }

  # Separator line (lf)
  incr idx 1
  $InputFileInfo(infoLB) insert end ""
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground white

  # Color files
  # -----------
  #
  incr idx 1
  $InputFileInfo(infoLB) insert end "          *** Color files ***"
  $InputFileInfo(infoLB) itemconfigure $idx -bg $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectbackground $Info(nonActiveBg)
  $InputFileInfo(infoLB) itemconfigure $idx -selectforeground black

  foreach row $Info(colorFiles) {
    set fn [lindex $row 0]
    set mt [lindex $row 1]
    set mts [Util::timeClock2Str $mt]
    if { 0 && $mts != "" } {
      $InputFileInfo(infoLB) insert end "$fn  ($mts)"
    } else {
      $InputFileInfo(infoLB) insert end "$fn"
    }
    incr idx
    lappend InputFileInfo(files) [list $idx $fn "CLR" $mt]
  }

  $InputFileInfo(infoLB) see 0
}


proc InputFileInfo::selectInfo {} {
  global InputFileInfo

  $InputFileInfo(browseButton) configure -state disabled
  $InputFileInfo(editButton) configure -state disabled
  $InputFileInfo(loadButton) configure -state disabled

  #-Selected info row
  set index [$InputFileInfo(infoLB) curselection]
  
  if { $index < 0 || $index == "" } return

  set is_file 0

  foreach row $InputFileInfo(files) {
    
    set idx [lindex $row 0]

    if { $idx == $index } {
      set is_file 1
      set InputFileInfo(currentFileInfo) [lrange $row 1 end]
      break
    }
  }

  if { $is_file == 0 } return
  
  $InputFileInfo(browseButton) configure -state normal
  $InputFileInfo(editButton) configure -state normal
  $InputFileInfo(loadButton) configure -state normal
}


proc InputFileInfo::browseFile {} {
  global InputFileInfo
  
  set fname [lindex $InputFileInfo(currentFileInfo) 0]

  if { $fname == "" } return

  FileBrowser::browse .winInputFileBrowse $fname
}


proc InputFileInfo::editFile {} {
  global InputFileInfo

  set fname [lindex $InputFileInfo(currentFileInfo) 0]

  if { $fname == "" } return

  FileBrowser::browse .winInputFileEdit $fname 1 1
}


proc InputFileInfo::loadFile {} {
  global InputFileInfo

  set fname [lindex $InputFileInfo(currentFileInfo) 0]
  set ftype [lindex $InputFileInfo(currentFileInfo) 1]

  switch $ftype {
    "DEF" {UserDefined::readDefinitionFile $fname}
    "KWD" {UserDefined::readSolverKeywordsFile $fname}
    "MTC" {Util::cpp_exec readMatcFile $fname}
    "CLR" {Util::cpp_exec readColorFile $fname}
  }
}


proc InputFileInfo::panelOk {w} {
  global InputFileInfo

  Panel::cancel $w
}


# end ecif_tk_InputFileInfoPanel.tcl
# ********************
