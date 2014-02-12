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
#
#*******************************************************************************
#*
#*  Argument & option checking utility.
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

#
# Check command arguments given in string args
# Parameters:
#
# Command: Command name
#
# Usage: Usage string
#
# OptCount: Number of options possible for the command
#
# IgnoreOptions: if set ignore unknown options (?)
#
# Opt: array of Options Opt(i,name) gives name of option i,
#       Opt(i,args) gives number of arguments the option i expects
#
# Val(option name): returns the option argument string if any or string "given"
#
# Minp,Maxp: minimum & maximum number of parameters for the command
#
# Arg(): returns the argument strings if any
#
# args: command string given by user
#
#
proc check_args { Command Usage OptCount IgnoreOptions Opt Val Minp Maxp Arg args } {
    upvar $Opt Options 
    upvar $Val OptValues
    upvar $Arg ArgValues

    set i 1
    set params 0

    set args [split [string trim $args "{}"] " "]
    set match 0

    while { $i <= [llength $args] } {

       set str [string tolower [lindex $args [@ $i-1]]]

       if { [string index $str 0] == "-" && ![regexp {[0-9]} [string index $str 1]] } {
           set match 0

           do j 1 $OptCount {
              set opt $Options($j,name)

              if { [string match $str* $opt] != 0 } {
                  incr match
                  set copt $j
                  if { $match > 1 } break
              }
           }

           if { !$IgnoreOptions && $match == 0 } { return -code error "$Command: Unknown option \[$str\].\n\n$Usage" }
       } else {
           set ArgValues($params) $str
           incr params

           if { $params > $Maxp } { return -code error "$Command: too many parameters.\n\n$Usage" }
       }

       if { !$IgnoreOptions && $match  > 1 } { return -code error "$Command: option \[$str\] not unique.\n\n$Usage" }

       if { $match == 1 } {
           if { $Options($copt,args) > 0 } {
               if { !$IgnoreOptions && $i >= [llength $args] } {
                   return -code error "$Command: no value given for option \[$Options($copt,name)\].\n\n$Usage"
               }
               set OptValues($Options($copt,name)) [lindex $args $i];
               incr i
           } else {
               set OptValues($Options($copt,name)) given
           }
       }

       incr i
    }

    if { $params < $Minp } { return -code error "$Command: too few parameters.\n\n$Usage" }
}
