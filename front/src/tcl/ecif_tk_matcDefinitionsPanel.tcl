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
#Module:    ecif_tk_matcDefinitionsPanel.tcl
#Language:  Tcl
#Date:      07.02.03
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining MATC functions and varaibles
#
#************************************************************************


#--Set MATC definitions
#
proc MatcDefinitions::openPanel {} {
  global MatcDefinitions Info ObjectTable Model
  upvar #0 MatcDefinitions theArray

  set w $MatcDefinitions(winName)
  set wgeom $MatcDefinitions(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set MatcDefinitions(winId) $id

  set Info(thisWindow) $w
  set this $w

  if { 1 == [Util::checkPanelWindow MatcDefinitions $id $MatcDefinitions(winTitle) $wgeom] } {
    raise $MatcDefinitions(winName)
    focus -force $MatcDefinitions(winName)
    return
  }  

  set MatcDefinitions(dataChanged) 0
  set MatcDefinitions(dataModified) 0

  set MatcDefinitions(appendMode) 1
  set MatcDefinitions(syncLog) 1
  set MatcDefinitions(wrapEntry) 0

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $MatcDefinitions(winTitle)
  wm geometry $w $wgeom


  #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #---Outer frames
  set f1 [frame $w.f1]     ;#Outer frame
  set f11 [frame $w.f1.f1] ;#-Listbox area
  set f12 [frame $w.f1.f2] ;#-Listbox control buttons area
  set f2  [frame $w.f2]    ;#Entry area
  set f21 [frame $w.f2.f1] ;#-Entry text area
  set f22 [frame $w.f2.f2] ;#-Entry control buttons area
  set f3  [frame $w.f3]    ;#Entry history area
  set f31 [frame $w.f3.f1] ;#-History text area
  set f32 [frame $w.f3.f2] ;#-History control buttons area
  set f4  [frame $w.f4]    ;#Result area
  set f41 [frame $w.f4.f1] ;#-Result text area
  set f42 [frame $w.f4.f2] ;#-Result control buttons area
  set f5  [frame $w.f5]    ;#Apply+Ok+cancel buttons frame

  set box_wid 50
  set btn_wid 7

  #-Definitions listbox
  #--------------------
  set m $f11
  label $m.l  -anchor nw -text "Current matc functions and variables defined by the user"
  set wdg [ listbox $m.lb -relief sunken \
             -selectmode extended -exportselection 0 \
             -height 8 -width $box_wid -font $Info(tableFont) \
             -xscrollcommand [list $m.sx set] \
             -yscrollcommand [list $m.sy set] ]

  set MatcDefinitions(defsLB) $wdg

  scrollbar $m.sx -orient horizontal -command [list $wdg xview]
  scrollbar $m.sy -orient vertical -command [list $wdg yview]

  MatcDefinitions::updateDefinitions

  pack $m.l -side top -anchor w
  pack $m.sx -side bottom -anchor n -fill x -expand 0
  pack $m.lb -side left -anchor w -fill both -expand 1
  pack $m.sy -side left -anchor n -fill y  -expand 0

  bind $wdg <Button-1> "MatcDefinitions::selectDefinition %x %y"
  bind $wdg <Double-1> "MatcDefinitions::editDefinition"
  bind $wdg <ButtonRelease-1> "MatcDefinitions::checkDefSelections"

  #-Listbox control buttons and checkboxes
  #
  set m $f12

  set append_box  [checkbutton $m.replace -text "Appnd" \
                   -variable MatcDefinitions(appendMode)   \
                   -offvalue 0 -onvalue 1 -state disabled   \
                   -indicatoron 1 -selectcolor $Info(select_color) ]
 
  set load_btn   [button $m.load -text "Load..." -command "MatcDefinitions::loadDefinitions" -width $btn_wid -state normal]
  set save_btn   [button $m.save -text "Save..." -command "MatcDefinitions::saveDefinitions" -width $btn_wid -state disabled]
  set edit_btn   [button $m.edit -text Edit -command "MatcDefinitions::editDefinition" -width $btn_wid -state disabled]
  set delete_btn [button $m.delete -text Delete -command "MatcDefinitions::deleteDefinition" -width $btn_wid -state disabled]

  set MatcDefinitions(saveButton) $save_btn
  set MatcDefinitions(appendBox) $append_box
  set MatcDefinitions(editButton) $edit_btn
  set MatcDefinitions(deleteButton) $delete_btn

  pack $load_btn $save_btn $append_box $edit_btn $delete_btn -side top -expand 1 -anchor w -padx $fpx1 -pady $fpy1


  #-Entry textbox
  #--------------
  set m $f21
  label $m.l -anchor nw -text "Enter a matc expression or definition"
  set wdg [text $m.text -relief groove -height 9 -width $box_wid \
                        -bg white -state normal \
                        -font $Info(entryFont) \
                        -setgrid true -wrap none \
                        -xscrollcommand [list $m.sx set] \
                        -yscrollcommand [list $m.sy set] ]

  bind $wdg <Control-Return> "MatcDefinitions::evalEntry"
  bind $wdg <Control-e> "MatcDefinitions::evalEntry"

  scrollbar $m.sx -orient horizontal \
    -command [list $m.text xview]
  scrollbar $m.sy -orient vertical  \
    -command [list $m.text yview]

  set MatcDefinitions(entryTB) $wdg

  pack $m.l -side top -anchor w
  pack $m.sx -side bottom -anchor n -fill x -expand 0
  pack $m.text -side left -anchor w -fill both -expand 1
  pack $m.sy -side left -anchor n -fill y  -expand 0

  #-Entry control buttons
  #
  set m $f22
  set eval_btn  [button $m.eval -text Eval -command "MatcDefinitions::evalEntry" -width $btn_wid -state normal]
  set clear_btn [button $m.clear -text Clear -command "MatcDefinitions::clearEntry" -width $btn_wid -state normal]

  set sync  [checkbutton $m.sync -text "Sync" -justify left \
               -variable MatcDefinitions(syncLog)    \
               -offvalue 0 -onvalue 1 -state normal   \
               -indicatoron 1 -selectcolor $Info(select_color) ]

  set wrap  [checkbutton $m.wrap -text "Wrap" -justify left \
               -variable MatcDefinitions(wrapEntry)    \
               -offvalue 0 -onvalue 1 -state normal   \
               -command MatcDefinitions::setEntryWrap  \
               -indicatoron 1 -selectcolor $Info(select_color) ]

  set clearLogs_btn  [button $m.clearLogs -text "Clear log" -command "MatcDefinitions::clearLogs" -width $btn_wid -state normal]
  set delLogRow_btn  [button $m.delLogRow -text "Del row" -command "MatcDefinitions::deleteLogRow" -width $btn_wid -state normal]
  
  set sep [frame $m.sep]

  pack $eval_btn $clear_btn -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy1
  pack $wrap -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy1
  pack $sep -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy3
  pack $clearLogs_btn -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy1
  pack $delLogRow_btn -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy1
  pack $sync -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy1


  #-Result log listbox
  #-------------------
  set m $f31
  label $m.l  -anchor nw -text "Result log"
  set wdg [ listbox $m.lb -relief sunken \
             -selectmode browse -exportselection 0 \
             -height 4 -width $box_wid -font $Info(tableFont) \
             -xscrollcommand [list $m.sx set] \
             -yscrollcommand [list $m.sy set] ]

  bind $wdg <Button-1> "MatcDefinitions::selectResultLog %x %y"

  set MatcDefinitions(resultLogLB) $wdg

  scrollbar $m.sx -orient horizontal -command [list $wdg xview]
  scrollbar $m.sy -orient vertical -command [list $wdg yview]

  pack $m.l -side top -anchor w
  pack $m.sx -side bottom -anchor n -fill x -expand 0
  pack $m.lb -side left -anchor w -fill both -expand 1
  pack $m.sy -side left -anchor n -fill y  -expand 0

  #-Result log control buttons
  #
  set m $f32
  set show_btn  [button $m.show -text Show -command "MatcDefinitions::showResult" -width $btn_wid -state normal]
  set funcs_btn  [button $m.help -text "Functions" -command "MatcDefinitions::functionList" -width $btn_wid  -state normal]
  #set clear_btn [button $m.clear -text Clear -command "MatcDefinitions::clearResult" -width $btn_wid -state normal]

  pack $show_btn $funcs_btn -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy2


  #-Entry log listbox
  #------------------
  set m $f41
  label $m.l  -anchor nw -text "Entry log"

  set wdg [ listbox $m.lb -relief sunken \
             -selectmode browse -exportselection 0 \
             -height 4 -width $box_wid -font $Info(tableFont) \
             -xscrollcommand [list $m.sx set] \
             -yscrollcommand [list $m.sy set] ]

  bind $wdg <Button-1> "MatcDefinitions::selectEntryLog %x %y"

  set MatcDefinitions(entryLogLB) $wdg

  scrollbar $m.sx -orient horizontal -command [list $wdg xview]
  scrollbar $m.sy -orient vertical -command [list $wdg yview]

  pack $m.l -side top -anchor w
  pack $m.sx -side bottom -anchor n -fill x -expand 0
  pack $m.lb -side left -anchor w -fill both -expand 1
  pack $m.sy -side left -anchor n -fill y  -expand 0

  #-Entry log control buttons
  #
  set m $f42
  set insert_btn  [button $m.insert -text Insert -command "MatcDefinitions::copyFromEntryLog 0" -width $btn_wid -state normal]
  set replace_btn [button $m.replace -text Replace -command "MatcDefinitions::copyFromEntryLog 1" -width $btn_wid -state normal]
  
  pack $insert_btn $replace_btn -side top -anchor w -expand 1 -padx $fpx1 -pady $fpy2
  

  #-Ok etc. buttons
  #----------------
  set m $f5

  if { $Model(hasMatcDefinitions) } {
    set upd normal
  } else {
    set upd disabled
  }

  set ca $Info(defaultCancelState)

  set ok_btn [button $m.ok -text OK -command "MatcDefinitions::panelOk $this"]
  set cn_btn [button $m.cancel -text Cancel -command "MatcDefinitions::panelCancel $this" \
                                -state $ca]
  set upd_btn [button $m.apply -text "Update geometry" -command "MatcDefinitions::updateCadGeometry" \
                               -state $upd]
  
  focus $ok_btn
  set MatcDefinitions(updateGeometryButton) $upd_btn
  #set MatcDefinitions(cancelButton) $cn_btn

  pack $ok_btn $upd_btn -side left -expand 1 -padx $fpx1
  #pack $ok_btn $cn_btn $upd_btn -side left -expand 1 -padx $fpx1
  #pack $ok_btn $cn_btn -side left -expand 1 -padx $fpx1
  
  pack $f1  -side top  -expand 1 -padx $fpx2 -pady $fpy2 -fill both
  pack $f11 -side left -expand 1 -padx $fpx1 -fill both -anchor w
  pack $f12 -side left -expand 0 -padx $fpx1 -anchor w
  pack $f2  -side top  -expand 1 -padx $fpx2 -pady $fpy2 -fill both
  pack $f21 -side left -expand 1 -padx $fpx1 -fill both -anchor w
  pack $f22 -side left -expand 0 -padx $fpx1 -anchor w
  pack $f3  -side top  -expand 1 -padx $fpx2 -pady $fpy2 -fill both
  pack $f31 -side left -expand 1 -padx $fpx1 -fill both -anchor w
  pack $f32 -side left -expand 0 -padx $fpx1 -anchor w
  pack $f4  -side top  -expand 1 -padx $fpx2 -pady $fpy2 -fill both
  pack $f41 -side left -expand 1 -padx $fpx1 -fill both -anchor w
  pack $f42 -side left -expand 0 -padx $fpx1 -anchor w
  pack $f5  -side top  -expand 0 -padx $fpx2 -pady $fpy2
}


# ======================#
# Definitions box procs #
# ======================#


# Load-button proc
#
proc MatcDefinitions::loadDefinitions {} {
  global MatcDefinitions

  set types {
    { {Matc} {.txt} }
    { {All} {*} }
  }
  set title "Load matc definitions"

  #--Display dialog and get path
  set path [ MenuExec::openFile $title $types ]

  if {$path == ""} {
    return
  }

  Util::cpp_exec readMatcFile $path

  MatcDefinitions::updateDefinitions

  MatcDefinitions::checkDefSelections
}


# Save-button proc
#
proc MatcDefinitions::saveDefinitions {} {
  global MatcDefinitions

  set types {
    { {Matc} {.txt} }
    { {All} {*} }
  }

  if { $MatcDefinitions(appendMode) } {
    set title "Save the selected matc definitions (Append mode)"
    set mode "a"
  } else {
    set title "Save the selected matc definitions (Replace mode)"
    set mode "w"
  }

  #--Display dialog and get path
  set path [MenuExec::saveFile $title $types $MatcDefinitions(defaultSaveFile)]

  if {$path == ""} {
    return
  }
  
  set MatcDefinitions(defaultSaveFile) $path

  #if { [catch { set fid [open $path $mode] } msg] } {
  #  return ""
  #}

  # Get listbox size
  #  
  set size [$MatcDefinitions(defsLB) size]
  set i 0

  # Store each selected listbox row to the file
  #
  set defs ""
  while { $i < $size } {

    if { [$MatcDefinitions(defsLB) selection includes $i] } {
      
      # Pick selected row
      set def [string trim [join [$MatcDefinitions(defsLB) get $i $i]]]
      
      lappend defs $def

      #puts $fid $def
    }
    incr i
  }

  #close $fid
  Util::cpp_exec updateMatcFile [list $path $mode $defs]
}


# Update definitions listbox
#
proc MatcDefinitions::updateDefinitions {} {
  global MatcDefinitions

  set defs [MatcDefinitions::getDefinitions]

  ListBox::fill $MatcDefinitions(defsLB) $defs
}


# Get list of current user defined matc-variables and functions
#
proc MatcDefinitions::getDefinitions {{add_dsign 0}} {
  global MatcDefinitions

  set defs ""

  # Get user defined variables
  #
  set vars [Util::getMatcVariables]
  foreach nm $vars {

    set val [string trim [Util::doMatc $nm]]

    set clean_val ""
    
    # Remove extra whitespace
    foreach v $val {
      append clean_val [string trim $v]
      append clean_val " "
    }

    if { $add_dsign } {
      lappend defs "\$$nm = $clean_val"
    } else {
      lappend defs "$nm = $clean_val"
    }

  }

  # Get user defined variables
  #
  set funcs [Util::getMatcFunctions]
  foreach nm $funcs {

    set fdef [string trim [Util::doMatc "help(\"$nm\")"]]

    if { $add_dsign } {
      lappend defs "\$function $fdef"

    } else {
      lappend defs "function $fdef"
    }
  }
  
  return $defs
}


# Event proc for definitons listbox selection
# Nothing useful currently!
#
proc MatcDefinitions::selectDefinition {lb_x lb_y} {
  global MatcDefinitions

  #-Selected ("clicked") definition row
  set index [$MatcDefinitions(defsLB) index @$lb_x,$lb_y]
  
  if { $index < 0 || $index == "" } return

  set row [string trim [join [$MatcDefinitions(defsLB) get $index $index]]]

  set row_info [Util::getMatcDefInfo $row]
#MSG "info=$row_info"

}


# Check selection status in the definitons listbox
#
proc MatcDefinitions::checkDefSelections {} {
  global MatcDefinitions

  set indices [$MatcDefinitions(defsLB) curselection]
  set scount [llength $indices]

  if { $scount == 0 } {
    $MatcDefinitions(saveButton) configure -state disabled
    $MatcDefinitions(appendBox) configure -state disabled
    $MatcDefinitions(deleteButton) configure -state disabled
  } else {
    $MatcDefinitions(saveButton) configure -state normal
    $MatcDefinitions(appendBox) configure -state normal
    $MatcDefinitions(deleteButton) configure -state normal
  }

  if { $scount == 0 || $scount > 1 } {
    $MatcDefinitions(editButton) configure -state disabled
  } else {
    $MatcDefinitions(editButton) configure -state normal
  }
}


# Edit-button command
#
proc MatcDefinitions::editDefinition {} {
  global MatcDefinitions

  set index [$MatcDefinitions(defsLB) curselection]
  
  # If multiple selection, cannot edit!
  if { 1 < [llength $index] } {
    return
  }

  # Ilegal index
  if { $index < 0 } {
    return
  }

  set row [join [$MatcDefinitions(defsLB) get $index $index]]
  
  $MatcDefinitions(entryTB) delete 0.0 end
  $MatcDefinitions(entryTB) insert 0.0 $row
}


# Delete-button command
#
proc MatcDefinitions::deleteDefinition {} {
  global MatcDefinitions

  set indices [$MatcDefinitions(defsLB) curselection]

  set count [llength $indices]

  if { $count == 0 || ($count == 1 && $indices < 0) } {
    return
  }

  # Loop from the last selected row upwards
  set idx $count
  while { $idx > 0 } {

    incr idx -1
    set index [lindex $indices $idx]

    set row [join [$MatcDefinitions(defsLB) get $index $index]]

    # Delete Matc definition
    #    
    set def_info [Util::getMatcDefInfo $row]
    set name [lindex $def_info 0]
    set type [lindex $def_info 1]
    
    if { $name != "" } {
      if { $type == "V" } {
        Util::doMatc "delete(\"$name\")"
      } elseif { $type == "F" } {
        Util::doMatc "funcdel(\"$name\")"
      }
    }

    # Delete listbox entry
    #
    $MatcDefinitions(defsLB) delete $index $index

    # Delete Matc object
    #
  }

  MatcDefinitions::checkDefSelections
}


# Function-button command (show list of all matc-functions)
#
proc MatcDefinitions::functionList {} {
  global MatcDefinitions
  
  set result [string trim [Util::doMatc "help"]]
  
  set help_call MatcDefinitions::functionHelp
  
  # Display function list text a separate window
  #  
  TextBrowser::browse .matcFuncs "Currently defined matc functions" $result 0 1 50 15 $help_call "Show help" 
}


# Help-button command in functions textbox 
#
proc MatcDefinitions::functionHelp {name} {
  global MatcDefinitions
  
  if { $name == "" } {
    set result "Select first a function name!"
  } else {
    set arg "help(\"$name\")"
    set result [string trim [Util::doMatc $arg 0]]
  }

  # Display in text a separate window
  #  
  TextBrowser::browse .matcHelp "Matc function help" $result 0 0 50 15
}



# ================#
# Entry box procs #
# ================#

# Eval-button command
#
proc MatcDefinitions::evalEntry {} {
  global MatcDefinitions
  
  # Get argument
  #  
  set ew $MatcDefinitions(entryTB)
  set arg [string trim [$ew get 0.0 end]]
  
  if { $arg == "" } {
    return
  }

  #--Evaluate
  #
  set result [string trim [Util::doMatc $arg]]

  #--Form a nice result
  #  
  set result [split $result]

  set result2 ""
  foreach res $result {
    set res [string trim $res]

    if { $res == "" } continue

    append result2 "$res "
  }

  #--Add expression/definition to the result log
  #
  # Was an expression
  if { $result2 != "" } {
    $MatcDefinitions(resultLogLB) insert 0 $result2
  
  # Was a definition (variable or function)
  } else {
    $MatcDefinitions(resultLogLB) insert 0 "DEF: $arg"
    MatcDefinitions::updateDefinitions
  }

  $MatcDefinitions(resultLogLB) selection clear 0 end
  $MatcDefinitions(resultLogLB) selection set 0 0
  $MatcDefinitions(resultLogLB) see 0

  #--Add expression/definition to the entry log
  #
  $MatcDefinitions(entryLogLB) insert 0 $arg

  $MatcDefinitions(entryLogLB) selection clear 0 end
  $MatcDefinitions(entryLogLB) selection set 0 0
  $MatcDefinitions(entryLogLB) see 0

}


# Clear entry textbox
#
proc MatcDefinitions::clearEntry {} {
  global MatcDefinitions
  
  set ew $MatcDefinitions(entryTB)

  $ew delete 0.0 end
}


# Set entry textbox's wrap-mode
#
proc MatcDefinitions::setEntryWrap {} {
  global MatcDefinitions

  set ew $MatcDefinitions(entryTB)

  if { $MatcDefinitions(wrapEntry) } {
    $ew configure -wrap word
  } else {
    $ew configure -wrap none
  }
}


# Clear result and entry log lisboxes
#
proc MatcDefinitions::clearLogs {} {
  global MatcDefinitions
  
  $MatcDefinitions(entryLogLB) delete 0 end
  $MatcDefinitions(resultLogLB) delete 0 end
}


# Delete a row in result and entry log listboxes 
# NOTE: Delete always corresponding rows in both listboxes
#
proc MatcDefinitions::deleteLogRow {} {
  global MatcDefinitions
  
  set index -1

  if { $index < 0 } {
    set index [$MatcDefinitions(entryLogLB) curselection]
  }
  
  if { $index < 0 } {
    set index [$MatcDefinitions(resultLogLB) curselection]
  }

  if { $index < 0 } return

  $MatcDefinitions(entryLogLB) delete $index $index
  $MatcDefinitions(entryLogLB) selection set $index $index

  $MatcDefinitions(resultLogLB) delete $index $index
  $MatcDefinitions(resultLogLB) selection set $index $index
}



# =====================#
# Result log box procs #
# =====================#

# Select-event proc for result log listbox
#
proc MatcDefinitions::selectResultLog {lb_x lb_y} {
  global MatcDefinitions

  #-Selected ("clicked") row
  set index [$MatcDefinitions(resultLogLB) index @$lb_x,$lb_y]
  
  if { $index < 0 || $index == "" } return

  set wdg $MatcDefinitions(entryLogLB)
  $wdg selection clear 0 end
  
  if { !$MatcDefinitions(syncLog) } return
  
  # If sync on, select corresponding row in entry log
  #
  $wdg selection set $index $index
  $wdg see $index
}


# Show-button command
#
# Show a result log row in a separate textbox
#
proc MatcDefinitions::showResult {} {
  global MatcDefinitions
  
  set index [$MatcDefinitions(resultLogLB) curselection]
  
  if { $index < 0 } {
    return
  }

  set row [join [$MatcDefinitions(resultLogLB) get $index $index]]

  # Display in text a separate window
  #  
  TextBrowser::browse .matcSHow "Matc result log display" $row 0 0 50 15
}


# ====================#
# Entry log box procs #
# ====================#

# Select-event proc for entry log listbox
#
proc MatcDefinitions::selectEntryLog {lb_x lb_y} {
  global MatcDefinitions

  #-Selected ("clicked") definition row
  set index [$MatcDefinitions(entryLogLB) index @$lb_x,$lb_y]
  
  if { $index < 0 || $index == "" } return

  set wdg $MatcDefinitions(resultLogLB)
  $wdg selection clear 0 end
  
  if { !$MatcDefinitions(syncLog) } return
  
  # If sync on, select corresponding row in result log
  #
  $wdg selection set $index $index
  $wdg see $index
}


# Copy-button command
#
# Copy an entry log row into entry textbox
#
proc MatcDefinitions::copyFromEntryLog {do_clear} {
  global MatcDefinitions

  set index [$MatcDefinitions(entryLogLB) curselection]

  if {$index < 0 } {
    return
  }

  set row [join [$MatcDefinitions(entryLogLB) get $index $index]]
  
  if { $do_clear } {
    $MatcDefinitions(entryTB) delete 0.0 end
  }

  $MatcDefinitions(entryTB) insert insert $row
}



# ===================#
# Panel button procs #
# ===================#

proc MatcDefinitions::panelOk {w} {
  global MatcDefinitions

  Panel::cancel $w
}


proc MatcDefinitions::updateCadGeometry {} {
  global MatcDefinitions ObjectTable
  
  MenuExec::verifyUpdateCadGeometry
}


proc MatcDefinitions::panelCheck {} {
  global MatcDefinitions

  #-Ok
  return 1
}

#ecif_MatcDefinitionsPanel.tcl

