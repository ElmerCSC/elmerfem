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
#* Timestep control
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
#*                       Date: 08 Jan 1996
#*
#*                Modified by:
#*
#*       Date of modification:
#*
#*******************************************************************************
#
# User level command for setting current timestep
#

set time_cmd ""
set time_loop "Loop"

proc timestep { args } {
     global time_tstep IsosurfaceRecompute time_current

     set Usage "Usage: timestep time\n\n"

     check_args timestep $Usage 0 0 opt opt_val 0 1 arg_val $args

     if { ![info exists arg_val] } { return $Usage }

     set sid [array startsearch arg_val]

     while { [array anymore arg_val $sid] != 0 } {

         set n [array nextelement arg_val $sid]
         set val $arg_val($n)
         set n [string range $n 1 [string length $n]]
     }

     array donesearch arg_val $sid

     set time_tstep $arg_val(0)

     set IsosurfaceRecompute 1
     c_TimeStep [@ $arg_val(0)-1]

     math time_current=[@ $time_tstep-1]
     math {
        if ( exists("times") ) { time_str=sprintf( "set time_current {Simulation timestep: %g, Simulation time: %g sec}",times(1,time_current) times(2,time_current) ); tcl(time_str); } } 

     UpdateObject
     display
}

#
# Loop over timesteps
#
# 09 Jan 1996
#
proc time_set_loop { {start 1} {end 0} {inc 1} {cmd ""} { count 1 } } {

   global time_loop time_edit time_step NumberOfTimesteps time_str Tt

   if { $end < 1 } { set end $NumberOfTimesteps }

   if { $time_loop == "Loop" } {
        set time_loop "Stop"
        do i 1 $count  {
            for { set t $start } { $t <= $end } { set t [@ $t+$inc] } {
               if { $time_loop  == "Loop" } return

#              math time_current=[@ $t-1]
#              math { time_str=sprintf( "set Tt Time:%gs", times(2,time_current) ); tcl(time_str); } 
#              teksti 0.4 -0.9 0 $Tt 20

               timestep $t;
               update; UpdateDisplay;
               if { $cmd != "" } { catch [eval $cmd]; }
            }
        }
   }

   set time_loop "Loop"
}

#
# Create a timestep editor toplevel window for a user to use.
#
# 09 Jan 1996
#
proc time_Edit { } {
   global time_edit time_tstep time_loop time_min time_max time_inc time_cmd NumberOfTimesteps
   global time_count time_current

   set time_edit .time_edit

   if { [winfo exists $time_edit] } {
       destroy  $time_edit.title   $time_edit.looplab $time_edit.loop    $time_edit.command $time_edit.stime
       destroy  $time_edit.buttons $time_edit.loopbut $time_edit.current $time_edit.dummy2  $time_edit.dummy3
   } else {
       toplevel $time_edit
       place_window $time_edit
   }


#  wm geometry $time_edit 300x250
   wm title $time_edit "Time Step Control"

   label $time_edit.title -text "Timestep Control"
   pack $time_edit.title -side top -expand 1 -fill both

   frame $time_edit.loop -bd 10 -relief ridge -bg lightblue

   frame $time_edit.looplab
   label $time_edit.looplab.label -text "Looping Controls"
   pack $time_edit.looplab.label
   set time_loop "Loop"

   set time_min 1
   label $time_edit.loop.minlab -text "Min: "
   entry $time_edit.loop.min -width 4 -textvariable time_min

   set time_max $NumberOfTimesteps
   label $time_edit.loop.maxlab -text "Max: "
   entry $time_edit.loop.max -width 4 -textvariable time_max

   set time_inc 1
   label $time_edit.loop.inclab -text "Inc: "
   entry $time_edit.loop.inc -width 4 -textvariable time_inc

   set time_count 1
   label $time_edit.loop.countlab -text "Loop Count: "
   entry $time_edit.loop.count -width 4 -textvariable time_count

   pack $time_edit.loop.minlab $time_edit.loop.min -side left -expand 1
   pack $time_edit.loop.maxlab $time_edit.loop.max -side left -expand 1
   pack $time_edit.loop.inclab $time_edit.loop.inc -side left -expand 1
   pack $time_edit.loop.countlab $time_edit.loop.count -side left -expand 1

   frame $time_edit.command
   label $time_edit.command.lab -text "Do after frame:"
   entry $time_edit.command.cmd -textvariable time_cmd -width 50
   pack $time_edit.command.lab -side left
   pack $time_edit.command.cmd -side left -expand 1

   frame $time_edit.stime -bd 10 -relief ridge -bg lightblue

#
# slider for time
#
   label $time_edit.stime.lab -text "Set timestep: "
   entry $time_edit.stime.ent -width 4 -textvariable time_tstep
   bind $time_edit.stime.ent <Return> { timestep $time_tstep }

#    slider $time_edit.stime.scl -orient horizontal -variable time_tstep -command "timestep" \
#      -from 1 -to $NumberOfTimesteps -resol 1 -tick [@ floor($NumberOfTimesteps/8)]
 
    pack $time_edit.stime.lab -side left
    pack $time_edit.stime.ent -side left
#    pack $time_edit.stime.scl -side left -fill x -expand 1
 
    frame $time_edit.buttons
 
    button $time_edit.buttons.close -text "Close" -command "destroy $time_edit"
    button $time_edit.buttons.help  -text "Help"  -command "help edit.html#Timestep"
 
    frame $time_edit.loopbut
    button $time_edit.loopbut.loop -textvariable time_loop -bd 6 -command \
     { time_set_loop $time_min $time_max $time_inc $time_cmd $time_count }
 
    pack  $time_edit.loopbut.loop -side left -expand 1 -fill x
 
    pack  $time_edit.buttons.help  -side right
    pack  $time_edit.buttons.close -side right
 
    pack $time_edit.stime -expand 1 -fill both -side top
 
    label $time_edit.current -textvariable time_current -relief sunken -bd 6
    pack $time_edit.current -side top -expand 1 -fill both
 
    label $time_edit.dummy2 -text " "
    pack $time_edit.dummy2 -side top
 
    pack $time_edit.looplab -expand 1 -fill both -side top
    pack $time_edit.loop -expand 1 -fill both -side top
    pack $time_edit.command -side top
 
    label $time_edit.dummy3 -text " "
    pack $time_edit.dummy3 -side top
 
    pack $time_edit.loopbut -side top -expand 1 -fill x
    pack  $time_edit.buttons -side top -fill x -expand 1
}

