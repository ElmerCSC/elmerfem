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
#* Initialize some global TCL variables, and include other TCL/TK commands
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

set tcl_interactive 1

if { [info exists env(ELMER_POST_HOME)] } {

   set ELMER_POST_HOME $env(ELMER_POST_HOME)

} else {

   set ELMER_POST_HOME "/mnt/mds/csc/jpr/SRC/ELMER/PostProcessor/"

}

frame .top -relief raised -bd 1
pack  .top -side top -fill both

frame .bot -relief raised -bd 1
pack  .bot -side top -fill both

message .top.msg -width 3i -text "Send Interrupt" -font -Adobe-Times-Medium-R-Normal-*-140-*
pack .top.msg -expand 1 -fill both -padx 5m -pady 5m

button .bot.button -bg gray -text "Send Interrupt" -command "SendInterrupt"
pack .bot.button -expand 1 -fill both -padx 5m -pady 5m
