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
#Module:    ecif_tk_systemInfoPanel.tcl
#Language:  Tcl
#Date:      04.06.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for system and graphics environment information.
#
#************************************************************************

##################
### SystemInfo ###
##################

# This procedure displays the model statistic-info screen
#
proc SystemInfo::openPanel { } {
  # Global variables
  global Info 

  set w .infowin
  set SystemInfo(winName) $w
  set SystemInfo(winTitle) "System environment info"
  set wgeom +420+120

  set id [winfo atom $w]
  set SystemInfo(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow SystemInfo $id $SystemInfo(winTitle) $wgeom] } {
    return
  }  

  toplevel $w
  focus $w 

  wm title $w $SystemInfo(winTitle)
  wm geometry $w $wgeom 


  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  set bg  $Info(nonActiveBgLight)

  #---Data fields 

  set Info(currentDirectory) [pwd]

  set fieldInfos {
    { "User name: "               USER }
    { "Computer name: "           HOST }
    { "Machine: "                 machine }
    { "Operating system: "        os }
    { "OS version: "              osVersion }
    { "Tcl version: "             tclVersion }
    { "Tcl level: "               tclLevel }  
    { "Tk version: "              tkVersion }
    { "Tk level: "                tkLevel }
    { "Tcl library path: "        tclLibraryPath }
    { "Tk library path: "         tkLibraryPath }  
    { "Elmer Front version: "     FRONT_VERSION_NBR }
    { "Elmer Front script path: " frontScriptPath }
    { "Elmer Home: "              ELMER_HOME }
    { "Elmer Front Home: "        ELMER_FRONT_HOME }
    { "Working directory: "       workingDirectory }
  }

  if { $Info(ELMER_FRONT_DEBUG_TCL) } {
    lappend fieldInfos { "Current directory: " currentDirectory }
  }

  # Create lables and entries
  set index 1
  foreach fi $fieldInfos {

    set m [frame $w.f$index]

    set lbl [lindex $fi 0]
    set var [lindex $fi 1]
    
    # Create label and entry
    label $m.l  -width 18 -text $lbl
    entry $m.e  -width 55 -textvariable Info($var) -justify left -state disabled -bg $bg

    pack $m -side top -expand 0 -fill x -padx $fpx1 -pady $fpy1
    pack $m.l $m.e -side left -anchor w -in $m 

    incr index
  }

  # Buttons
  frame $w.fB 

  #-Ok button
  set ok_btn [button $w.fB.ok -text OK -command "Panel::cancel $SystemInfo(winName)"]
  focus $ok_btn

  pack $ok_btn
  pack $w.fB -side top -expand 0 -fill x -padx $fpx3 -pady $fpy3
}


# Dummy proc
proc SystemInfo::panelOk {w} {
}



####################
### RendererInfo ###
####################

# This procedure displays the model statistic-info screen
#
proc RendererInfo::openPanel { } {
  # Global variables
  global Info RendererInfo

  set w .infowin
  set RendererInfo(winName) $w
  set RendererInfo(winTitle) "Graphics environment info"
  set wgeom +420+120

  set id [winfo atom $w]
  set RendererInfo(infoWinId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow RendererInfo $id $RendererInfo(winTitle) $wgeom] } {
    return
  }  
  
  toplevel $w
  focus $w 

  wm title $w $RendererInfo(winTitle)
  wm geometry $w $wgeom 

  #-We set each row of data into a separate frame.
  frame $w.f1   ;# Line width granularity
  frame $w.f2   ;# Line width range
  frame $w.fB   ;# Buttons

  set wid1 12; set wid2 18; set wid3 24; set wid4 30
  #-A header text for each data field.
  label $w.f1.l  -width $wid2 -anchor w -text "Line width range: "
  label $w.f2.l  -width $wid2 -anchor w -text "Line width granularity: "

  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  set bg  $Info(nonActiveBgLight)

  #---Data fields (as labels)
  set exp 0
  set start 1
  set variable_count 2

  set variables {
    LINE_WIDTH_RANGE
    LINE_WIDTH_GRANULARITY
  } 

  set space        {15 7}
  set justify_left {1  1}

  # Create lables and entries
  for {set index 0} {$index < $variable_count} {incr index} {
    set m $w.f[expr $start + $index]
    set s  [lindex $space $index]
    set j_left [lindex $justify_left $index]
    set vars [lindex $variables $index]
    
    pack $m -side top -expand $exp -fill x -padx $fpx1 -pady $fpy1
    pack $m.l -side left -anchor w -in $m 

    set anch c
    # Should we left justify data
    if {$j_left} { set just left}

    # Create entries
    set count 0
    foreach v $vars {
      incr count
      set val [set RendererInfo($v)]
      entry $m.d$count -width $s -textvariable RendererInfo($v) -justify $just -state disabled -bg $bg 
      pack $m.d$count -side left -in $m -padx $fpx1
    }
  }

  #-Ok button
  set ok_btn [button $w.fB.ok -text OK -command "Panel::cancel $RendererInfo(winName)"]
  focus $ok_btn

  pack $ok_btn
  pack $w.fB -side top -expand $exp -fill x -padx $fpx3 -pady $fpy3
}

# Dummy proc
proc RendererInfo::panelOk {w} {
}

########################
### DEBUGGING TCL/TK ###
########################


# ================
# Sourcelist Panel
# ================

# This procedure displays all source files, you can also reload
# these files!
#
proc sourceListPanel { } {
  # Global variables
  global Info 

  set w .sourcewin
  set SourceList(winName) $w
  set SourceList(winTitle) "Source file list"
  set wgeom +420+120

  set id [winfo atom $w]
  set SourceList(infoWinId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow SourceList $id $SourceList(winTitle) $wgeom] } {
    return
  }  
 
  toplevel $w
  focus $w 

  wm title $w $SourceList(winTitle)
  wm geometry $w $wgeom 

  frame $w.fLB
  frame $w.fB

  set m $w.fLB

  listbox $m.sourceList -relief sunken -selectmode extended -exportselection 0 \
    -yscrollcommand [list $m.sy set] \
    -height 30 -width 55 -font $Info(entryFont)

  scrollbar $m.sy -orient vertical -command [list $m.sourceList yview]

  pack $m.sourceList -side left -anchor n -fill x -expand 1
  pack $m.sy -side left -anchor n -fill y 

  foreach pair $Info(sourceScripts) {
    set nam [lindex $pair 0]
    set fil [lindex $pair 1]
    set tim [lindex $pair 2]
 
    # Mark modified
    if { $tim != [file mtime $fil] } {
      append nam "       ***"
    }

    $m.sourceList insert end $nam
  }
 
  #-We set each row of data into a separate frame.
  #-Ok button
  set ok_btn [button $w.fB.ok -text OK -command "Panel::cancel $SourceList(winName)"]

  #-Load buttons
  set ld_sel_btn [button $w.fB.load_sel -text "Load selected" -command "sourceListLoadSelected $m.sourceList"]
  set ld_mod_btn [button $w.fB.load_mod -text "Load modified" -command "sourceListLoadModified $m.sourceList"]
  set ld_all_btn [button $w.fB.load_all -text "Load all"      -command "sourceListLoadAll $m.sourceList"]

  focus $ok_btn

  pack $ok_btn $ld_sel_btn $ld_mod_btn $ld_all_btn -side left -expand 1 -fill x

  pack $w.fLB $w.fB -side top -expand 1 -fill x -padx 10 -pady 10
}


# Dummy proc
proc SourceList::panelOk {w} {
}


# Can be used without opening the widget!
#
proc updateSource {} {
  sourceListLoadModified
}


proc sourceListLoadSelected {lb_wdg} {
  global Info

  set indices [$lb_wdg curselection]

  if { $indices == "" } {
    return
  }

  foreach index $indices {
    set src_def [lindex $Info(sourceScripts) $index]
    set src_file [lindex $src_def 1]

    set mtime [file mtime $src_file]
    set src_def [lreplace $src_def 2 2 $mtime]

    set Info(sourceSripts) [lreplace $Info(sourceScripts) $index $index $src_def]

    source $src_file
    Message::showMessage "Source loaded: $src_file"
  }

  sourceListRefill $lb_wdg
}


proc sourceListLoadModified { {lb_wdg ""} } {
  global Info

  set new_list ""

  set index 0
  foreach src_def $Info(sourceScripts) {

    set src_file [lindex $src_def 1]
    set src_time [lindex $src_def 2]

    set mtime [file mtime $src_file]

    # Check if modified time has changed
    if { $mtime != $src_time } {
      set src_def [lreplace $src_def 2 2 $mtime]
      set src_file [lindex $src_def 1]
      source $src_file
      Message::showMessage "Source loaded: $src_file"
    }

    lappend new_list $src_def
    incr index
  }

  set Info(sourceScripts) $new_list

  sourceListRefill $lb_wdg
}


proc sourceListLoadAll { {lb_wdg ""} } {
  global Info

  set new_list ""

  set index 0
  foreach src_def $Info(sourceScripts) {

    set src_file [lindex $src_def 1]

    set mtime [file mtime $src_file]
    set src_def [lreplace $src_def 2 2 $mtime]

    source $src_file
    Message::showMessage "Source loaded: $src_file"

    lappend new_list $src_def
    incr index
  }

  set Info(sourceScripts) $new_list

  sourceListRefill $lb_wdg

}


proc sourceListRefill {lb_wdg} {
  global Info

  if { $lb_wdg == "" } {
    return
  }

  $lb_wdg delete 0 end

  foreach src_def $Info(sourceScripts) {

    set src_name [lindex $src_def 0]

    $lb_wdg insert end $src_name
  }

}


# =================
# Tcl command panel
# =================

# This procedure displays an entry and text output widget for
# evaluatin tcl commands
#
proc tclCommandPanel {} {
  global Info TclCommand

  set id 0

  while { [winfo exists .commandwin$id] } {
    incr id
  }

  set w .commandwin$id
  set winName $w
  set winTitle "Tcl command entry"
  set wgeom +420+120
 
  toplevel $w
  focus $w 

  wm title $w $winTitle
  wm geometry $w $wgeom 

  if { ![info exists TclCommand(cmdListIndex,$id)] } {
    set TclCommand(cmdListIndex,$id) -1
  }
   
  frame $w.fT
  frame $w.fE
  frame $w.fB

  # Command output text widget
  # ==========================
  set m $w.fT
  set wdg [text $m.cmdOutput  -relief sunken  \
                              -height 20 -width 55 -setgrid true -wrap word \
                              -font $Info(entryFont) \
                              -xscrollcommand [list $m.xscroll set] \
                              -yscrollcommand [list $m.yscroll set] ]

  scrollbar $m.xscroll -orient horizontal -command [list $wdg xview]
  scrollbar $m.yscroll -orient vertical -command [list $wdg yview]

  set TclCommand(cmdOutputWidget,$id) $wdg 

  pack $m.xscroll -side bottom -fill x
  pack $m.yscroll -side right -fill y
  pack $wdg -side left -fill both -expand true

  # Command entry widget
  # ====================
  set m $w.fE
  set wdg [entry $m.cmdEntry -width 45 -textvariable TclCommand(cmd,$id) -font $Info(entryFont)]

  set TclCommand(cmdWidget,$id) $wdg 

  focus $wdg

  bind $wdg <Control-x> "TclCommand::panelOk $id $winName"

  bind $wdg <Control-c> "tclCommandOutputClear $id"

  bind $wdg <KeyPress-Return> "tclCommandEval $id 1"
  bind $wdg <Control-Return> "tclCommandEval $id 0"

  bind $wdg <Control-t> "tclCommandCmdClear $id"
  bind $wdg <Control-d> "tclCommandCmdDelete $id"
  bind $wdg <Control-r> "tclCommandCmdReset $id 1"
  bind $wdg <Control-R> "tclCommandCmdReset $id 2"
  bind $wdg <Up> "tclCommandCmdUp $id"
  bind $wdg <Down> "tclCommandCmdDown $id"

  bind $wdg <Double-1> "tclCommandCmdUp $id; $TclCommand(cmdWidget,$id) selection clear"
  bind $wdg <Control-1> "tclCommandEval $id 1"
  bind $wdg <Double-3> "tclCommandCmdDown $id; $TclCommand(cmdWidget,$id) selection clear"

  pack $wdg -side left -fill both -expand true

  # Buttons
  # =======

  #-Ok button
  set ok_btn [button $w.fB.ok -text OK -command "TclCommand::panelOk $id $winName"]

  #-Clear button
  set clear_btn [button $w.fB.clear -text "Clear out" -command "tclCommandOutputClear $id"]
  
  #-Reset button  
  set reset_btn [button $w.fB.reset -text Reset -command "tclCommandCmdReset $id"]

  #focus $ok_btn

  pack $ok_btn $reset_btn $clear_btn -side left -anchor c  

  pack $w.fT $w.fE $w.fB -side top -expand 1 -fill x -padx 10 -pady 10
}
 

# Dummy proc
proc TclCommand::panelOk {id w} {
  global TclCommand

  #unset TclCommand(cmd,$id)

  destroy $w
}

proc tclCommandOutputClear {id} {
  global TclCommand

  $TclCommand(cmdOutputWidget,$id) delete 0.0 end
}


proc tclCommandCmdReset {id {level 1}} {
  global TclCommand
  
  # Clear command list
  set TclCommand(cmdList,$id) ""
  set TclCommand(cmdListIndex,$id) -1

  if { $level < 2 } {
    return
  }

  # Clear command entry
  tclCommandCmdClear $id

  # Clear output
  tclCommandOutputClear $id


}


proc tclCommandCmdClear {id} {
  global TclCommand

  set TclCommand(cmd,$id) ""
}


proc tclCommandCmdDelete {id} {
  global TclCommand

  set tmp $TclCommand(cmdList,$id)
  set idx $TclCommand(cmdListIndex,$id)

  set TclCommand(cmdList,$id) [lreplace $tmp $idx $idx]

  tclCommandCmdUp $id
}


proc tclCommandCmdUp {id} {
  global TclCommand

  if { $TclCommand(cmdListIndex,$id) > 0 } {
    incr TclCommand(cmdListIndex,$id) -1
  }

  set idx $TclCommand(cmdListIndex,$id)
  set tmp [lindex $TclCommand(cmdList,$id) $idx]

  set TclCommand(cmd,$id) $tmp
}


proc tclCommandCmdDown {id} {
  global TclCommand

  if { $TclCommand(cmdListIndex,$id) < \
       [expr [llength $TclCommand(cmdList,$id)] - 1]
     } {
    incr TclCommand(cmdListIndex,$id)
  }

  set idx $TclCommand(cmdListIndex,$id)
  set tmp [lindex $TclCommand(cmdList,$id) $idx]

  set TclCommand(cmd,$id) $tmp
}


proc tclCommandEval {id {add_to_list 1} } {
  global TclCommand

  $TclCommand(cmdWidget,$id) selection clear 

  set TclCommand(cmd,$id) [string trim $TclCommand(cmd,$id)]

  if { $TclCommand(cmd,$id) == "" } {
    return
  }

  if { [catch { set eval_res [eval "$TclCommand(cmd,$id)" ] } msg] } {
    set res $msg
  } else {
    set res $eval_res
  }

  append show_res "*** $TclCommand(cmd,$id) "

  if { $res == "" } {
    append show_res "no result"
  } else {
    append show_res "$res"
  }

  append show_res "\n"

  $TclCommand(cmdOutputWidget,$id) insert end $show_res
  $TclCommand(cmdOutputWidget,$id) see end

  if { !$add_to_list } {
    return
  }

  lappend TclCommand(cmdList,$id) $TclCommand(cmd,$id)

  set TclCommand(cmdListIndex,$id) [llength $TclCommand(cmdList,$id)]
  incr TclCommand(cmdListIndex,$id) -1
}


# Print array variable info
#
proc print_array_info { arr pattern print_values } {
  global $arr
  upvar #0 $arr theArr

  if { ![array exists $arr] } {
    return "array does not exist!"
  }

  if { $pattern == "" } {
    set values [array get $arr]
  } else {
    set values [array get $arr $pattern]
  }

  set array_rows ""

  set len [llength $values]
  set counter 0

  while { $counter < $len } {
    set pair ""
    lappend pair [lindex $values $counter]; incr counter
    lappend pair [lindex $values $counter]; incr counter
    lappend array_rows $pair
  }

  set array_rows [lsort $array_rows]

  set len [llength $array_rows]
  append result "  ($len elements)\n"

  foreach row $array_rows {

    # Name
    set res [lindex $row 0]

    # Value
    if {$print_values} {
      append res " = "
      set val [lindex $row 1]
      if { $val != "" } {
        append res $val
      }
    }

    append res "\n"

    append result $res
  }

  return $result
}



# Command window Tcl utiliy commands
# 
# ===================================================================
# pv: print variable                      ( var [component] )       #
# an: print names in a array              ( arrayname [varpattern]) # 
# pa: print array names and their values  ( arrayname [varpattern]) #
# ===================================================================


# Print global variable "part1"
# If "part2" given, we have an array!
#
proc pv { part1 { part2 ""} } {
  upvar #0 $part1 thePart1

  append result "\n"
  if { $part2 == "" } {
    append result $thePart1
  } else {
    append result $thePart1($part2)
  }

  return $result
}


# Print array names
#
proc an { arr {pattern ""} } {
  print_array_info $arr $pattern 0
}


# Print array names and values
#
proc pa { arr {pattern ""} } {
  print_array_info $arr $pattern 1
}


# Output linefeed
#
proc lf {} {
  return "\n"
}

# Get array field property
#
proc gfp {globArray fld prop {sys_index ""} } {
  return [DataField::getFieldProperty $globArray $fld $prop $sys_index]
}


# Show current stack till level "depth"
#
proc showStack  { {depth 12} } {
  global Info
  
  if { !$Info(ELMER_FRONT_DEBUG_TCL_PRINT) } {
    return
  }

  set trc "Level 1"

  foreach lvl {-1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12} {

    if { ![catch {set trc [info level $lvl] }] } {
      MSG "$lvl: $trc"
    }

    if {$depth == -$lvl} {
      break
    }
  }

  MSG "$trc"
}


# Set write-trace for an array field
#
proc setTrace {globArray fld} {
  upvar #0 $globArray theArray

  trace variable theArray($fld) w _trace
}


# Remove op-type trace for an array field
#
proc removeTrace {globArray fld {op "w"} } {
  upvar #0 $globArray theArray

  trace vdelete theArray($fld) $op _trace
}


# Tracing proc for an array field tracing
# n1 array name
# n2 field name
# op traced operation (w(rite), r(ead))
#
proc _trace {n1 n2 op} {
  upvar #0 $n1 theArr
  
  MSG "----"

  set val "?"
  catch {set val $theArr($n2) }
  MSG "$n1\($n2\) \[$op\]: \[val=$val\]"

  foreach lvl {-1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12} {

    if { ![catch {set trc [info level $lvl] }] } {
      MSG "$lvl: $trc"
    }
  }

  MSG "----"
}

# end ecif_tk_rendererInfoPanel.tcl
# ********************
