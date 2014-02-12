#/*****************************************************************************
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

#*******************************************************************************
#*
#*     Model file read utility routines
#*
#*******************************************************************************
#*
#*                     Author:       Juha Ruokolainen
#*
#*                    Address: CSC - IT Center for Science Ltd.
#*                                Keilaranta 14, P.O. BOX 405
#*                                  02101 Espoo, Finland
#*                                  Tel. +358 0 457 2723
#*                                Telefax: +358 0 457 2302
#*                              EMail: Juha.Ruokolainen@csc.fi
#*
#*                       Date: 26 Sep 1995
#*
#*                Modified by:
#*
#*       Date of modification:
#*
#*******************************************************************************

set OUTPS      0
set PSFileName ""

set DVector "Displacement"
set DScale  "1"
set DCycles "1"
set DMode   "1"
set DFrame  "20"
set DLoop   "Animate"
set DCmd    ""

proc mode_animate {} {
  global DCycles DVector DMode DScale DFrame DLoop DCmd

  if { $DLoop == "Animate" } {
    set DLoop "Stop"
    math function f(d,t) import nodes {  _f=d(0:2,time(t-1)) };
    math n=nodes;

    do t 1 [@ $DCycles*$DFrame] {
      if { $DLoop == "Animate" } {
        math nodes=n;
        return
      }
      math t=$t/$DFrame*2*pi
      math nodes=n+sin(t)*f($DVector,$DMode)*$DScale
      timestep $DMode
      update; display;
     if { $DCmd != "" } { catch [eval $DCmd]; }
    }
  }
  math nodes=n;
  set DLoop "Animate"
}

proc ModeDisplay {} {
    set w .modedisplay

    if { [winfo exists $w] } {
      wm iconify $w
      wm deiconify $w
      return
    } else {
      toplevel $w
      place_window $w
      wm title $w "Mode Display"
    }

    label $w.title -text "Mode Display"
    pack $w.title

    label $w.sp1 -text ""
    pack $w.sp1 -side top

    set DVector "Displacement"
    frame $w.arrow
    label $w.arrow.label -text "Displacement Variable: "
    button $w.arrow.but -textvariable DVector \
               -command { set DVector [make_vector_list]; }

    pack $w.arrow -side top
    pack $w.arrow.label -side left
    pack $w.arrow.but -side left -fill x


    frame $w.file -relief ridge
    label $w.file.mlab -text "Select Mode: "
    entry $w.file.mode -width 10 -textvariable DMode
    label $w.file.dlab -text "Disp scale: "
    entry $w.file.disp -width 10 -textvariable DScale
    label $w.file.flab -text "Frames/Cycle: "
    entry $w.file.frms -width 10 -textvariable DFrame
    label $w.file.clab -text "Cycles: "
    entry $w.file.cycl -width 10 -textvariable DCycles
    pack $w.file.mlab $w.file.mode $w.file.dlab $w.file.disp
    pack $w.file.flab $w.file.frms $w.file.clab $w.file.cycl
# -side left

    label $w.file.sp4 -text ""
    pack $w.file.sp4 -side top
    pack $w.file -side top -expand 1 -fill both

    frame $w.command
    label $w.command.lab -text "Do after frame:"
    entry $w.command.cmd -textvariable DCmd -width 30
    pack $w.command.lab -side left
    pack $w.command.cmd -side left -expand 1
    pack $w.command

    frame $w.anim -relief ridge
    button $w.anim.but -textvariable DLoop -command { mode_animate; } -relief ridge
    pack $w.anim.but -side left -fill both -expand 1
    pack $w.anim -side top -expand 1 -fill both

    frame $w.but
    button  $w.but.exit -text "Close" -command "destroy $w"
    pack $w.but $w.but.exit -side right
}
