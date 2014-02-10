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
#Module:    ecif_tk_processTablePanel.tcl
#Language:  Tcl
#Date:      28.01.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for monitoring ELMER processes
#
#************************************************************************
 
#------ProcessTable definitions  proc------

###################
### BUILD PANEL ###
###################

#This procedure displays the process-table panel
#
proc ProcessTable::openPanel { } {
  global Info Common ProcessTable ProcessIds ProcessListedIds

  set w .processTableWin
  set ProcessTable(winName) $w
  set ProcessTable(winTitle) "Process table"
  set wgeom +420+120

  set id [winfo atom $w]
  set ProcessTable(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow ProcessTable $id $ProcessTable(winTitle) $wgeom] } {
    return
  }  

  toplevel $w
  focus $w 

  set this $w

  wm title $w $ProcessTable(winTitle)
  wm geometry $w $wgeom

  set globArray "ProcessTable"

  set ProcessListedIds ""
  set ProcessTable(stepInfoFrozen) 0

  #-----WIDGET CREATION
  #-Data is organized into two principal frames:

  #----Process table
  frame $w.fP -bd 2 -relief groove

  #----Ok, Cancel etc. Buttons
  frame $w.fB 
  
  set m $w.fP

  #--Headers
  frame $m.lf
  frame $m.lf.lf1
  frame $m.lf.lf2

  # Table header
  #label $m.lf.lf1.l  -anchor nw -text "Processes:"

  # Column headers
  # ==============

  append hdr [format %-10s Process]
  append hdr [format %-4s  Nbr]
  append hdr [format %-6s  PID]
  append hdr [format %-8s  Prior.]
  append hdr [format %-10s State]
  append hdr [format %-14s Start]
  append hdr [format %-16s End]
  append hdr [format %-9s  Duration]
  append hdr [format %-1s  W]
  append hdr [format %-1s  " "]
  append hdr [format %-1s  E]

  label $m.lf.lf2.l  -text $hdr -anchor w -font $Info(processTableFont) 

  # Process Listbox
  # ===============

  frame $m.lbf
  frame $m.lbf.f1
  frame $m.lbf.f2

  listbox $m.lb -width 80 -height 4 \
                -relief sunken -selectmode browse \
                -xscrollcommand [list $m.lb_scrollbar_x set] \
                -yscrollcommand [list $m.lb_scrollbar_y set] \
                -exportselection 0 -font $Info(processTableFont)
  scrollbar $m.lb_scrollbar_x -orient horizont -command [list $m.lb xview]
  scrollbar $m.lb_scrollbar_y -orient vertical -command [list $m.lb yview]
  set ProcessTable(processLB) $m.lb
  bind $m.lb <ButtonRelease-1> "ProcessTable::info_"


  # State buttons
  # ==============

  frame $m.btnf
  frame $m.sep1
  frame $m.sep2
  frame $m.sep3

  button $m.kill     -text "Kill" -width 6 -command ProcessTable::kill  
  button $m.suspend  -text "Suspend" -width 6 -command ProcessTable::suspend  
  button $m.resume   -text "Resume" -width 6 -command ProcessTable::resume  
  button $m.delete   -text "Delete" -width 6 -command ProcessTable::delete
  button $m.browse   -text "Browse" -width 6 -command ProcessTable::browse
  set wdg [menubutton $m.priority -text "Priority" -relief raised -indicatoron 1]
  bind $wdg <ButtonPress-1> "ProcessTable::priority $wdg"

  # General process info area
  # =========================
  frame $m.genInfof
  frame $m.genInfof.f1
  frame $m.genInfof.f2

  text $m.genInfo -width 80 -height 8 -wrap none \
            -relief sunken  \
            -xscrollcommand [list $m.genInfo_sb_x set] \
            -yscrollcommand [list $m.genInfo_sb_y set] \
            -font $Info(processTableFont)
  scrollbar $m.genInfo_sb_x -orient horizont -command [list $m.genInfo xview]
  scrollbar $m.genInfo_sb_y -orient vertical -command [list $m.genInfo yview]
  set ProcessTable(generalInfoArea) $m.genInfo

  # Process step info area
  # ======================
  frame $m.stepInfof
  frame $m.stepInfof.f1
  frame $m.stepInfof.f2

	text $m.stepInfo -width 80 -height 4 -wrap none \
            -relief sunken  \
            -xscrollcommand [list $m.stepInfo_sb_x set] \
            -yscrollcommand [list $m.stepInfo_sb_y set] \
            -font $Info(processTableFont)
  scrollbar $m.stepInfo_sb_x -orient horizont -command [list $m.stepInfo xview]
  scrollbar $m.stepInfo_sb_y -orient vertical -command [list $m.stepInfo yview]
  set ProcessTable(stepInfoArea) $m.stepInfo
  
  #
  #-----WIDGET PACKING
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #-Packing principal frames
  pack $w.fP -side top -anchor nw -expand 1 -fill y -padx $fpx1 -pady $fpy1
  pack $w.fB -side top -padx $fpx3 -pady $fpy3 -expand 0


  #---Processes
  set m $w.fP
  pack $m.lf $m.lbf $m.genInfof $m.stepInfof $m.btnf -side top  -anchor w -fill both -expand 0

  #--Headers
  pack $m.lf.lf1 $m.lf.lf2  -side top  -anchor w
  pack $m.lf.lf2.l -side left -anchor w  

  #--List box
  pack $m.lbf.f1 $m.lbf.f2 -side top -fill x
  pack $m.lb $m.lb_scrollbar_y -in $m.lbf.f1 -side left -fill y 
  #pack $m.lb_scrollbar_x -in $m.lbf.f2 -side bottom -fill x

  #--Command Buttons
  pack $m.suspend $m.resume -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.sep1 -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.kill -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.sep2 -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.delete -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.sep3 -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.browse -side left -in $m.btnf -anchor w -padx $fpx1 -pady $fpy1
  pack $m.priority -side left -in $m.btnf -anchor w -padx $fpx3 -pady $fpy1

  # Priority setting only available in Win32!
  if { $Info(platform) != "windows" } {
    Widget::configureS $m.priority disabled
  }
 
  #--General info area
  pack $m.genInfof.f1 $m.genInfof.f2 -side top -fill x
  pack $m.genInfo $m.genInfo_sb_y -in $m.genInfof.f1 -side left -fill y
  pack $m.genInfo_sb_x -in $m.genInfof.f2 -side bottom -fill x

  #--Step info area
  pack $m.stepInfof.f1 $m.stepInfof.f2 -side top -fill x
  pack $m.stepInfo $m.stepInfo_sb_y -in $m.stepInfof.f1 -side left -fill y
  pack $m.stepInfo_sb_x -in $m.stepInfof.f2 -side bottom -fill x

  #pack $m -side top  -anchor nw -fill x -expand 1 -padx $fpx1 -pady $fpy1

  #---Freeze and Ok buttons
  #set ok_btn [button $w.fB.ok -text OK -command "ProcessTable::panelOk $this"]
  #pack $ok_btn -side left -padx $fpx1 

  set freeze_btn [button $w.fP.freezeStepInfo -width 6 -text Freeze -command ProcessTable::freezeStepInfo]
  pack $freeze_btn -side left -in $m.btnf -anchor w -padx $fpx2 -pady $fpy1
  set ProcessTable(freezeStepInfoButton) $freeze_btn

  set ok_btn [button $w.fP.ok -text OK -command "ProcessTable::panelOk $this"]
  pack $ok_btn -side left -in $m.btnf -anchor w -padx $fpx2 -pady $fpy1

  focus $ok_btn   


  #-Init process list listbox
  set ProcessTable(entryAdded) 0
  ProcessTable::updateProcessList 0
  ProcessTable::info_
  ProcessTable::refresh
}


######################
### PANEL COMMANDS ###
######################

# Suspend (pause) a running process
#
proc ProcessTable::suspend {} {
  global ProcessTable ProcessListedIds

  set index [$ProcessTable(processLB) curselection]

  if { $index == "" } {
    return
  }

  set nbr [lindex $ProcessListedIds $index]

  if {$ProcessTable($nbr,state) != "running"} {
    return
  }

  if { [Process::suspend $nbr] } {

    set ProcessTable($nbr,updated) 0
    set ProcessTable($nbr,state) "suspended"
    set ProcessTable($nbr,stopTime) [clock seconds]

    ProcessTable::updateProcessList $index
  }

}


# Resume (continue) a suspended process
#
proc ProcessTable::resume {} {
  global Info ProcessTable ProcessIds ProcessListedIds

  set index [$ProcessTable(processLB) curselection]

  if { $index == "" } {
    return
  }

  set nbr [lindex $ProcessListedIds $index]

  if { $ProcessTable($nbr,state) != "suspended" } {
    return
  }

  if { [Process::resume $nbr] } {

    set ProcessTable($nbr,updated) 0
    set ProcessTable($nbr,state) "running"
    set ProcessTable($nbr,lastStartTime) [clock seconds]

    ProcessTable::updateProcessList $index
  }

}


# Stop a running/suspended process
# NOTE: Number is given ex. when killing from
# browser window
#
proc ProcessTable::kill { {nbr ""} } {
  global Info ProcessTable ProcessListedIds

  set index_given false

  # If no process number given, use selected row
  if { $nbr == "" } {
    set index [$ProcessTable(processLB) curselection]
 
    if { $index == "" } {
      return
    }
    
    set index_given true
    set nbr [lindex $ProcessListedIds $index]
  }

  if { $ProcessTable($nbr,terminated) } {
    return
  }

  append p_info "Process: $ProcessTable($nbr,name) "
  append p_info "Nbr: $nbr "
  append p_info "PID: $ProcessTable($nbr,pid) "
  set msg [list "NOTE: $p_info will be terminated!\n\n"  \
           $Info(anywayOk) ]

  if { ![ Message::verifyAction $Info(advancedUser) $msg] } {
    return
  }

  # Exec kill comand
  # ================
  if { [Process::kill $nbr] } {

    set ProcessTable($nbr,updated) 0
    set ProcessTable($nbr,state) "killed"
    set ProcessTable($nbr,terminated) 1
    set ProcessTable($nbr,stopTime) [clock seconds]

    if {$index_given} {
      ProcessTable::updateProcessList $index
    }
  }

}


# Delete a process from the list
#
proc ProcessTable::delete {} {
  global Info ProcessTable ProcessIds ProcessListedIds

  set index [$ProcessTable(processLB) curselection]

  if { $index == "" } {
    return
  }

  set nbr [lindex $ProcessListedIds $index]

  set id_index [lsearch $ProcessIds $nbr]

  set ProcessListedIds [lreplace $ProcessListedIds $index $index]

  append p_info "Process: $ProcessTable($nbr,name) "
  append p_info "Nbr: $nbr "
  append p_info "PID: $ProcessTable($nbr,pid) "


  # NOTE: An already stopped process is not confirmed!

  # Exec Kill and delete
  # ====================
  if { !$ProcessTable($nbr,terminated) } {
    set msg [list "NOTE: $p_info will be terminated and removed from the list!\n\n"  \
                  $Info(anywayOk) ]

    if { ![ Message::verifyAction $Info(advancedUser) $msg ] } {
      return
    }

    # Stop process
    if { ![Process::kill $nbr] } {
      return
    }
  }

  # Delete process data
  Util::doWait 50

  ProcessTable::deleteEntry $nbr $id_index

  set sindex [llength $ProcessListedIds]
  incr sindex -1

  ProcessTable::updateProcessList $sindex
}


proc ProcessTable::browse {} {
  global ProcessTable ProcessListedIds

  set index [$ProcessTable(processLB) curselection]

  if { $index == "" } {
    return
  }

  set nbr [lindex $ProcessListedIds $index]

  if { $ProcessTable($nbr,logfile) == "" } {
    return
  }

  MenuExec::browseProcessLog $nbr ""
}


proc ProcessTable::priority {menu_btn} {
  global ProcessTable ProcessListedIds Info

  set index [$ProcessTable(processLB) curselection]

  if { $index == "" } {
    return
  }

  set nbr [lindex $ProcessListedIds $index]

  if { $ProcessTable($nbr,terminated) } {
    return
  }

  if { [winfo exist $menu_btn.menu] } {
    destroy $menu_btn.menu
  }

  # Build selection menu for the process
  set wdg [menu $menu_btn.menu -title "Select priority" -tearoff 0]

  set cur_priority $ProcessTable($nbr,priority)

  $wdg add radiobutton -label "Favor GUI" \
    -variable cur_priority \
    -value $Info(ElmerRun,lowPriority) \
    -command "Process::setPriorityLevel $nbr $Info(ElmerRun,lowPriority)"

  $wdg add radiobutton -label "Favor process"  \
    -variable cur_priority \
    -value $Info(ElmerRun,normalPriority) \
    -command "Process::setPriorityLevel $nbr $Info(ElmerRun,normalPriority)"

  $menu_btn configure -menu $wdg


  ProcessTable::updateProcessList
  #$menu_btn invoke
}


### Freeze/Unfreeze button proc ###

proc ProcessTable::freezeStepInfo {} {
  global ProcessTable

  if { $ProcessTable(stepInfoFrozen) == 1 } {
    set ProcessTable(stepInfoFrozen) 0
    $ProcessTable(freezeStepInfoButton) configure -text "Freeze" -fg black
  } else {
    set ProcessTable(stepInfoFrozen) 1
    $ProcessTable(freezeStepInfoButton) configure -text "Unfreeze" -fg red
  }
}


### OK proc ###

proc ProcessTable::panelOk {w} {
  Panel::cancel $w
}



##############################
### PROCESS TABLE HANDLING ###
##############################

proc ProcessTable::addEntry { pid runnumber priority
                              show_in_list check_all_done
                              process_tag  logfile
                              starttime
                              {process_info ""}
                            } {

  global ProcessTable ProcessIds Info ModelProperty

  set ProcessTable($runnumber,pid) $pid

  set ProcessTable($runnumber,priority) $priority

  set ProcessIds [linsert $ProcessIds 0 $runnumber]

  # NOTE: Channel-value is set from cpp-side after process has been
  # started!
  #
  if { ![info exists ProcessTable($runnumber,channel)] } {
    set ProcessTable($runnumber,channel) ""
  }

  if { ![info exists ProcessTable($runnumber,browserDelay)] } {
    set ProcessTable($runnumber,browserDelay) $Info(browser,mediumDelay)
  }

  if { $logfile != "" } {
    set tmp $process_info
    set process_info "Logfile:     $logfile\n"
    append process_info $tmp
  }

  set ProcessTable($runnumber,exists) 1
  set ProcessTable($runnumber,removeable) 0
  set ProcessTable($runnumber,show) $show_in_list
  set ProcessTable($runnumber,checkAllDone) $check_all_done
  set ProcessTable($runnumber,model) $ModelProperty(MODEL_NAME)
  set ProcessTable($runnumber,name) $process_tag
  set ProcessTable($runnumber,logfile) $logfile
  set ProcessTable($runnumber,startTime) $starttime
  set ProcessTable($runnumber,lastStartTime) $starttime
  set ProcessTable($runnumber,stopTime) ""
  set ProcessTable($runnumber,stopTimeDisplayed) ""
  set ProcessTable($runnumber,duration) ""
  set ProcessTable($runnumber,state) "running"
  set ProcessTable($runnumber,selected) 1
  set ProcessTable($runnumber,statusDone) 0
  set ProcessTable($runnumber,terminated) 0
  set ProcessTable($runnumber,updated) 0
  set ProcessTable($runnumber,processInfo) $process_info
  set ProcessTable($runnumber,activeSystem) -1
  set ProcessTable($runnumber,activeSubSystem) -1

  # NOTE: These are needed  currently only for solver step-info output
  set ProcessTable($runnumber,activeIds) ""
  set ProcessTable($runnumber,currentTime) ""
  set ProcessTable($runnumber,currentStep1) ""
  set ProcessTable($runnumber,currentStep2) ""
  set ProcessTable($runnumber,currentStep3) ""
  set ProcessTable($runnumber,printStep1) -1
  set ProcessTable($runnumber,printStep2) -1
  set ProcessTable($runnumber,printStep3) -1

  Util::unsetArrayVariables ProcessTable $runnumber,stepInfo,*
  set ProcessTable($runnumber,stepInfo,nofRows) 0

  Util::unsetArrayVariables ProcessTable $runnumber,solver*
  set ProcessTable($runnumber,solver,procNames) ""
  set ProcessTable($runnumber,solver,procIndices) ""
  set ProcessTable($runnumber,solver,procSubIndices) ""
  set ProcessTable($runnumber,solver,nofWarnings) 0
  set ProcessTable($runnumber,solver,nofErrors) 0
  set ProcessTable($runnumber,solver,outputFiles) ""
  set ProcessTable($runnumber,solver,postFiles) ""
  set ProcessTable($runnumber,solver,resultIndex) $Info(NO_INDEX)
  set ProcessTable($runnumber,solver,resultSubIndex) $Info(NO_INDEX)

  # Inform about new entry!
  set ProcessTable(entryAdded) 1
}


proc ProcessTable::deleteEntry {runnumber ids_index } {
  global ProcessTable ProcessIds

  set ProcessIds [lreplace $ProcessIds $ids_index $ids_index]

  if { $ProcessTable($runnumber,logfile) != "" } {
    set msg ""
    catch { [file delete -force $ProcessTable($runnumber,logfile)] msg }

    if { $msg != "" } {
      Message::showMessage $msg
    } else {
      Message::showMessage "Logfile for the process $runnumber ($ProcessTable($runnumber,logfile) deleted!"
    } 
  }
  
  # NOTE: We do not delete all fields automatically
  # Existence related fields are kept unless clients
  # have flagged them removable
  #
  catch {
    unset ProcessTable($runnumber,pid)
    unset ProcessTable($runnumber,exists)
    unset ProcessTable($runnumber,priority)
    unset ProcessTable($runnumber,show)
    unset ProcessTable($runnumber,checkAllDone)
    unset ProcessTable($runnumber,model)
    unset ProcessTable($runnumber,name)
    unset ProcessTable($runnumber,logfile)
    unset ProcessTable($runnumber,browserDelay)
    unset ProcessTable($runnumber,startTime)
    unset ProcessTable($runnumber,lastStartTime)
    unset ProcessTable($runnumber,stopTime) 
    unset ProcessTable($runnumber,stopTimeDisplayed)
    unset ProcessTable($runnumber,statusDone)
    unset ProcessTable($runnumber,selected)
    unset ProcessTable($runnumber,duration)
    unset ProcessTable($runnumber,updated)
    unset ProcessTable($runnumber,processInfo)
    unset ProcessTable($runnumber,activeSystem)
    unset ProcessTable($runnumber,activeSubSystem)

    unset ProcessTable($runnumber,activeIds)
    unset ProcessTable($runnumber,currentTime)
    unset ProcessTable($runnumber,currentStep1)
    unset ProcessTable($runnumber,currentStep2)
    unset ProcessTable($runnumber,currentStep3)
    unset ProcessTable($runnumber,printStep1)
    unset ProcessTable($runnumber,printStep2)
    unset ProcessTable($runnumber,printStep3)

    Util::unsetArrayVariables ProcessTable $runnumber,stepInfo,*
    
    # Solver specific
    Util::unsetArrayVariables ProcessTable $runnumber,solver*
  }

  # By default we leave these
  if { $ProcessTable($runnumber,removeable) } {
    catch {
      unset ProcessTable($runnumber,state)
      unset ProcessTable($runnumber,terminated)
      unset ProcessTable($runnumber,removeable)
    }
  }
}


proc ProcessTable::updateProcessList { {selection_index ""} } {
  global Equation ProcessTable ProcessIds ProcessListedIds SolverSystem

  $ProcessTable(processLB) delete 0 end

  set ProcessListedIds ""

  # Fill listbox
  # ============
  foreach nbr $ProcessIds {

    # If process in not be shown in  the list
    if { !$ProcessTable($nbr,show) } {
      continue
    }
    
    lappend ProcessListedIds $nbr


    # Format start and duration
    # =========================
    #set startT [clock format $ProcessTable($nbr,startTime) -format "%d.%m %H:%M:%S"]
    set startT [clock format $ProcessTable($nbr,startTime) -format "%d.%m %H:%M"]
    set runT   [Util::timeSec2Str $ProcessTable($nbr,duration)]
 
    # Format displayed stop time
    # ==========================
    # Non finished process
    if { $ProcessTable($nbr,state) != "all done" &&
         $ProcessTable($nbr,state) != "failed" && 
         $ProcessTable($nbr,state) != "killed"
       } {
      set stopT  $ProcessTable($nbr,stopTimeDisplayed)

    # Finished prosess
    } else {
      set stopT  [clock format $ProcessTable($nbr,stopTime)  -format "%d.%m %H:%M"]
    }

    switch $ProcessTable($nbr,priority) {
      LOW_PRIORITY          { set priority LOW }
      BELOW_NORMAL_PRIORITY { set priority LOW }
      NORMAL_PRIORITY       { set priority NORM }
      ABOVE_NORMAL_PRIORITY { set priority NORM }
      HIGH_PRIORITY         { set priority HIGH }
      default               { set priority "" }
    }

    # Form one data row
    # =================
    set r ""
    append r [format "%-10s" $ProcessTable($nbr,name)]
    append r [format "%-4s"  $nbr]
    append r [format "%-6s"  $ProcessTable($nbr,pid)]
    append r [format "%-8s"  $priority]
    append r [format "%-10s" [string toupper $ProcessTable($nbr,state)]]

    append r [format "%-14s" $startT]
    append r [format "%-12s" $stopT]
    append r [format "%-13s"  $runT]

    if { $ProcessTable($nbr,solver,nofWarnings) > 0 } {
      append r [format "%-1s"  "*"]
    } else {
      append r " "
    }
    
    append r " "

    if { $ProcessTable($nbr,solver,nofErrors) > 0 } {
      append r [format "%-1s"  "*"]
    } else {
      append r " "
    }

    # Mark updated
    set ProcessTable($nbr,updated) 1

    $ProcessTable(processLB) insert end $r

  } ;# Fill listbox

  if { $selection_index != "" } {
    set sidx $selection_index
  } else {
    set sidx 0
  }

  $ProcessTable(processLB) selection set $sidx $sidx
}

# This proc update ProcessTable's generic and step info area widgets
#
proc ProcessTable::info_ {} {
  global ProcessTable ProcessListedIds

  # Check that table is visible!
  #
  if { ![info exist ProcessTable(processLB)] ||
       ![winfo exist $ProcessTable(processLB)]
     } {
    return
  }

  set msg ""
  catch {
    set index [$ProcessTable(processLB) curselection]

    set iw $ProcessTable(generalInfoArea)

    #--Pick the current top of the window x,y-position
    set xpos [lindex [$iw xview] 0]
    set ypos [lindex [$iw yview] 0]

    #--Update info window contents (first delete all, then insert
    #  updated text!)

    $iw delete 0.0 end

    if { $index == "" } {
      return
    }

    set nbr [lindex $ProcessListedIds $index]
  
    $iw insert 0.0 $ProcessTable($nbr,processInfo)

    #--Note: The previous insert command resets the window to the beginning
    #  of the text, so we have to restore the position...
    $iw xview moveto $xpos
    $iw yview moveto $ypos

    #--Update step info contents
    $ProcessTable(stepInfoArea) delete 0.0 end

    for {set row 1 } { $row <= $ProcessTable($nbr,stepInfo,nofRows) } { incr row } {
      
      if { ![info exists ProcessTable($nbr,stepInfo,$row)] } continue

      $ProcessTable(stepInfoArea) insert 0.0 "$ProcessTable($nbr,stepInfo,$row)\n"
    }

  } msg ;# Catch message

  if { $msg != "" } {
    ;#MSG $msg
  }
}


# ========================================================
# ========================================================
# 
# Update process status info (only solver steps currently)
# 
# ========================================================
# ========================================================
#
# NOTE: Calls ProcessTable::info_ (=show info) if process table visible
#
proc ProcessTable::updateStatusInfo {visible} {
  global ProcessTable ProcessIds
  
  set current_nbr -1
  set showInfo 0

  foreach nbr $ProcessIds {
    
    # Killed process
    if { [string equal -nocase $ProcessTable($nbr,state) "killed"] } {
      set ProcessTable($nbr,statusDone) 1
    }

    # If process status done
    if { $ProcessTable($nbr,statusDone) } {
      continue
    }

    # If process not be shown
    if { !$ProcessTable($nbr,show) } {
      continue
    }

    #--No step info
    if { ![info exist ProcessTable($nbr,currentStep1)] } {
      continue
    }
      
    # Check if process info needs update
    # ==================================

    if { !( [string equal -nocase $ProcessTable($nbr,state) "running"] ||
            [string equal -nocase $ProcessTable($nbr,state) "all done"]
          )
       } {
      continue
    }
      
    set current_nbr $nbr

    # A process found, may show info
    #
    if { $visible || !$ProcessTable(stepInfoFrozen) } {
      set showInfo 1
    }

    # Find correct print-step
    # =======================
    
    set unprinted_found 0
    set outstep1 -1
    set outstep2 -1

    # Some info from a previous step may not be printed yet, so
    # we must find the oldest non-printed Step1 and Step2 and these
    # are then used as keys for all systems when step to be printed is selected
    #
    # NOTE: We loop systems 'backwards' in order to pick the last
    # system which have non-printed system steps!
    #
    set ids $ProcessTable($nbr,activeIds)
    set len [llength $ids]

    for {set i 0} {$i < $len} {incr i} {

      set id [lindex $ids end-$i]

      if { ![info exists ProcessTable($nbr,$id,outputSteps)]  ||
           0 == [llength $ProcessTable($nbr,$id,outputSteps)]
         } {
        continue
      }
      
      # Latest printed steps
      set pstep1 $ProcessTable($nbr,$id,printStep1)
      set pstep2 $ProcessTable($nbr,$id,printStep2)
      set pstep3 $ProcessTable($nbr,$id,printStep3)

      # Output steps
      foreach steps $ProcessTable($nbr,$id,outputSteps) {
      
        set ostep1 [lindex $steps 0]
        set ostep2 [lindex $steps 1]
        set ostep3 [lindex $steps 2]
        
        if { ($ostep1 == $pstep1 && $ostep2 == $pstep2 && $ostep3 > $pstep3) ||
             ($ostep1 == $pstep1 && $ostep2 > $pstep2) ||
             ($ostep1 > $pstep1)
           } {
          set unprinted_found 1
          set outstep1  $ostep1
          set outstep2  $ostep2
          break
        }
      }

      if { $unprinted_found } {
        break
      }

      if { $unprinted_found } {
        break
      }
    }

    if { !$unprinted_found } {

      # There is just now no need to modify output (ex to remove remaining
      # activity (*) markers
      if { ![string equal -nocase $ProcessTable($nbr,state) "all done"] } {
        return
      }
    }

    # Form output line
    # ================
    if { $unprinted_found } {

      set si ""

      if { $ProcessTable($nbr,currentTime) != "" } {
        set is_transient 1
      } else {
        set is_transient 0
      }

    
      # Transient
      # ---------
      if { $is_transient } {
        #append si [format "%4d" $ProcessTable($nbr,currentStep1)]
        append si [format "%4d" $outstep1]
        append si " "

      # Steady state
      # ------------
      } else {
        append si "ssi "
      }

      #--Find latest Step3 for each system using Outstep1 and Outstep2 as keys
      #
      foreach id $ProcessTable($nbr,activeIds) {

        if { ![info exists ProcessTable($nbr,$id,outputSteps)]  ||
             0 == [llength $ProcessTable($nbr,$id,outputSteps)]
           } {
          continue
        }

        if { $id == $ProcessTable($nbr,activeSystem) } {
          set is_active 1
        } else {
          set is_active 0
        }
      
        set ProcessTable($nbr,$id,printStep1) $outstep1
        set ProcessTable($nbr,$id,printStep2) $outstep2

        set outstep3 -1
        set index -1
        set item -1
        
        # Find latest outstep3
        #
        foreach steps $ProcessTable($nbr,$id,outputSteps) {
        
          incr item

          set step1 [lindex $steps 0]
          set step2 [lindex $steps 1]
          set step3 [lindex $steps 2]
        
          if { ($step1 == $outstep1) &&
               ($step2 == $outstep2) &&
               ($step3 > $outstep3)
             } {
            set outstep3 $step3
            set index $item
          }
        }

        set ProcessTable($nbr,$id,printStep3) $outstep3
   
        set ps1 $ProcessTable($nbr,$id,printStep1)
        set ps2 $ProcessTable($nbr,$id,printStep2)
        set ps3 $ProcessTable($nbr,$id,printStep3)

        # NOTE: We have to know when systems are still printing something
        # for the current step in order NOT to copy old timestep/ssi info
        # for the new row. Steps are interpreted somewhat differently for
        # transient/steady case, so we need different conditions here.
        # Ok, this is based more on experimental work than first class
        # logic, but it seems to work! However, this mess should be cleared out!!!
        #
        set is_in_current 0

        if { $is_transient } {
          if { $ps1 == $ProcessTable($nbr,currentStep1) &&
               ($ps2 > 1 || $ps3 > 0)
             } {
             set is_in_current 1
          }
        } else {
          if { $ps1 == $ProcessTable($nbr,currentStep1) &&
               $ps2 >= 1 &&
               $ps3 >= 1
             } {
            set is_in_current 1
          }
        }

        # If nothing found for the system, but system already in the current step, 
        # show previous iteration info!
        #
        if { $is_in_current &&
             ( $index < 0 || "" == [lindex $ProcessTable($nbr,$id,outputData) $index] )
           } {
          set index $ProcessTable($nbr,$id,lastOutputIndex)
        }
        
        # Pick system output info
        #
        if { $index != -1 } {
          set sdata [lindex $ProcessTable($nbr,$id,outputData) $index]        
          set ProcessTable($nbr,$id,lastOutputIndex) $index
        } else {
          set sdata ""
        }
        
        # Only an active system name (with asterix)  is displayed without
        # any data!
        #
        if { !$is_active && $sdata == "" } {
          continue
        }

        # Mark active (=currently solving) system with asterix
        #
        set sn $ProcessTable($nbr,$id,systemName)

        if { $is_active } {
          append si "$sn*"
        } else {
          append si "$sn "
        }
      
        # Add system's data
        #
        append si $sdata

      } ;# Each active system
       
    } ;# Unprinted data found

    
    # Store info-row count and info
    #
    set old_count $ProcessTable($nbr,stepInfo,nofRows)
    
    # Add new output data
    #
    if { $unprinted_found } {
      set ProcessTable($nbr,stepInfo,nofRows) $outstep1
      set ProcessTable($nbr,stepInfo,$outstep1) $si
    }

    set new_count $ProcessTable($nbr,stepInfo,nofRows)

    if { $old_count == $new_count &&
         [string equal -nocase $ProcessTable($nbr,state) "all done"] } {
      set ProcessTable($nbr,statusDone) 1
    }

    
    # Cosmetic hack: remove possible remaining activity markers from the old row
    #
    if { $old_count < $new_count &&
         [info exists ProcessTable($nbr,stepInfo,$old_count)]
       } {
      set row $ProcessTable($nbr,stepInfo,$old_count)
      set row [string map  {"*" " "} $row]
      set ProcessTable($nbr,stepInfo,$old_count) $row
    }

    # Cosmetic hack: if all done, remove possible remaining activity markers from the last row
    #
    if { $ProcessTable($nbr,statusDone) &&
         [info exists ProcessTable($nbr,stepInfo,$new_count)]
       } {
      set row $ProcessTable($nbr,stepInfo,$new_count)
      set row [string map  {"*" " "} $row]
      set ProcessTable($nbr,stepInfo,$new_count) $row
    }


  } ;# Each process nbr


  if { $showInfo } {
    ProcessTable::info_

  } else {
    ;#Nothing
  }

}


proc ProcessTable::refresh {} {
  global Info ProcessTable 

  set visible 0

  if { [info exist ProcessTable(winName)] &&
       [winfo exist $ProcessTable(winName)]
     } {
    
    set visible 1
      
    if { $ProcessTable(entryAdded) } {
      set index 0
      set ProcessTable(entryAdded) 0
    } else {    
      set index [$ProcessTable(processLB) curselection]
    }
    
    ProcessTable::updateProcessList $index
  }

  ProcessTable::updateStatusInfo $visible
  
  after $Info(processTableRefreshInterval) ProcessTable::refresh
}


# Remove all nonactive processes from process table
# if they are not shown in the process list panel
#
proc ProcessTable::trim {} {
  global ProcessTable ProcessIds

  set current_nbrs $ProcessIds
   
  set index 0
  foreach nbr $current_nbrs {
    
    if { !$ProcessTable($nbr,show) &&
         $ProcessTable($nbr,terminated)
       } {
      ProcessTable::deleteEntry $nbr $index
    }
    incr index
  }
}


#####################
### PROCESS INFOS ###
#####################

# Form Solver process general info field (output files etc.)
#
proc ProcessTable::formSolverInfo {} {
  global Info Model ModelProperty Status

  set md $ModelProperty(MODEL_DIRECTORY,absolute)
  set mn $ModelProperty(MODEL_NAME)
  set pn $ModelProperty(PROBLEM_NAME)

  set pinfo ""
  append pinfo "Model:        $md $mn $pn \n"
  append pinfo "Steps:        [string trim $Status(timesteps)] \n"
  append pinfo "Input file:   [MenuExec::getSolverInputFile] \n"
  append pinfo "Meshes:       [string trim $Model(activeMeshNames)] \n"
  append pinfo "Output files: [MenuExec::getMeshRelatedFiles $Model(activeMeshNames) OUTPUT_FILE 0] \n"
  append pinfo "Post files:   [MenuExec::getMeshRelatedFiles $Model(activeMeshNames) POST_FILE 0] \n"
  return $pinfo
}


# Form Mesh process info field
#
proc ProcessTable::formMeshInfo {} {
  global ModelProperty Status
  set pinfo ""
  append pinfo "Model name:  $ModelProperty(MODEL_NAME) \n"
  #append pinfo "Mesh-H:      $Status(meshParameter) \n"
  return $pinfo
}



########################
### PROCESS HANDLING ###
########################

proc Process::start {number} {
  global Info ProcessTable

  #-WIN32
  if { $Info(platform) == "windows" } {
    Util::cpp_exec processStart $number
  #-Unix
  } else {
    set pid $ProcessTable($number,pid)
    exec kill -CONT $pid
  }

  return 1
}


proc Process::kill {number} {
  global Info ProcessTable

  set processTable($number,state) "killing"

  #-WIN32
  if { $Info(platform) == "windows" } {
    Util::cpp_exec processStop $number
  #-Unix
  } else {
    set pid $ProcessTable($number,pid)
    exec kill -KILL $pid
  }
  
  #set pid $Info(ElmerRun,pipe)
  #fconfigure $pid -blocking 0
  #catch {close $pid} msg
  #catch {puts $pid "Interrupted\n"; close $pid} msg
  #puts $Info(ElmerRun,pipe) 
  #Util::cpp_exec processStop

  return 1
}


proc Process::suspend {number} {
  global Info ProcessTable

  #-WIN32-NT
  if { $Info(platform) == "windows" } {
    Util::cpp_exec processSuspend $number
  #-Unix
  } else {
    set pid $ProcessTable($number,pid)
    exec kill -STOP $pid
  }

  return 1
}


proc Process::resume {number} {
  global Info ProcessTable

  #-WIN32
  if { $Info(platform)== "windows" } {
    Util::cpp_exec processResume $number
  #-Unix
  } else {
    set pid $ProcessTable($number,pid)
    exec kill -CONT $pid
  }
  return 1
}


proc Process::setPriorityLevel {number priority} {
  global ProcessTable Info

  # Stored changed priority level
  if { $priority == $ProcessTable($number,priority) } {
    return
  } else {
    set ProcessTable($number,priority) $priority
  }

  #-WIN32
  if { $Info(platform) == "windows" } {
    Util::cpp_exec processSetPriorityLevel "$number  $priority"

  #-Unix (not implemented!)
  } else {
  }

  ProcessTable::updateProcessList

  return 1
}


# Proc checks if a external Elmer process with
# sequence number 'nbr' exists
#
proc Process::exists_old {number} {
  global Info ProcessTable

  #-WIN32
  # NOTE: In Win32 we update the global
  # variable Process(number,exists) via cpp!
  if { $Info(platform) == "windows" } {
    Util::cpp_exec processExists $number

  #-Unix
  } else {
   set pid $ProcessTable($number,pid)

   set ProcessTable($number,exists) 0
   # Check process state from system's process "file"
   if { [file exists /proc/$pid] } {

     if { [file size /proc/$pid] > 0 } {
       set ProcessTable($number,exists) 1
     }
   }

  }

  return $ProcessTable($number,exists)
}


# Proc checks if a external Elmer process with
# sequence number 'nbr' exists
#
# NOTE Modified by Juha Ruokolainen to work also
# in Linux, 28.09.00 MVe
#
proc Process::exists {number} {
  global Info ProcessTable

  #-WIN32
  # NOTE: In Win32 we update the global
  # variable Process(number,exists) via cpp!
  if { $Info(platform) == "windows" } {
    Util::cpp_exec processExists $number

  #-Unix
  } else {
    set pid $ProcessTable($number,pid)

    set ProcessTable($number,exists) 0

#--Tested in various Unixes by JPR, 1/2001

    #-Check process state from system's process "file"
    if { [file exists /proc/$pid] } {

      if { [file exists /proc/$pid/stat] } {

        set fp [open "/proc/$pid/stat" RDONLY]
        set str [string trim [gets $fp]]
        close $fp
        scan $str "%d %s %s" ppid cmd stat

        if { $stat != "Z" } {
          set ProcessTable($number,exists) 1
        }

      } else {
        if { ![catch { set stat [exec ps -p $pid | grep $pid] } msg] } {
          if { [lsearch $stat "<defunct>"] < 0 } {
            set ProcessTable($number,exists) 1
          }
        }
#       if { [file size /proc/$pid] > 0 } {
#         set ProcessTable($number,exists) 1
#       }
      }

    #-Check process state by calling ps
    } else {
      if { ![catch { set stat [exec ps -p $pid | grep $pid] } msg] } {
        if { [lsearch $stat "<defunct>"] < 0 } {
          set ProcessTable($number,exists) 1
        }
      }
    }
  } ;#-end Unix test

  return $ProcessTable($number,exists)
}


# Mark a process removeable, ie. all field can
# removed when it it deleted from the table
#
proc Process::markRemoveable {nbr} {
  global ProcessTable

  if { [info exist ProcessTable($nbr,state)] } {
    set ProcessTable($nbr,removeable) 1
  }
}


proc Process::monitor {nbr} {
  global Info ProcessTable

  #---Process active
  if { [Process::exists $nbr] } {

    # Update duration
    # ===============
    if { $ProcessTable($nbr,state) == "running" } {
      
      set currentT [clock seconds]

      set ProcessTable($nbr,duration) \
            [expr $ProcessTable($nbr,duration) + \
                  $currentT - \
                  $ProcessTable($nbr,lastStartTime)
            ]
      set ProcessTable($nbr,lastStartTime) $currentT
    }

    after $Info(processTableMonitorInterval) "Process::monitor $nbr"

  #---Process inactive,  but still in process-table
  } elseif {[info exists ProcessTable($nbr,state)] } {

    set ProcessTable($nbr,updated) 0

    # Process ended normally or failed
    # ================================
    if { $ProcessTable($nbr,state) != "killed" } {

      set ProcessTable($nbr,terminated) 1

      # Write stop info
      set ProcessTable($nbr,stopTime) [clock seconds]

      if { ![info exist ProcessTable($nbr,solver,currentStep1)] ||
           $ProcessTable($nbr,solver,currentStep1) == ""
         } { 
        set ProcessTable($nbr,stopTimeDisplayed) [clock seconds]
      }

      #--If process should send "ALL DONE" for succes
      if { $ProcessTable($nbr,checkAllDone) } {

        if { [Process::succesful $nbr] } {
          set ProcessTable($nbr,state) "all done"

        } else {
          set ProcessTable($nbr,state) "failed"
          set msg "***NOTE Process $nbr ($ProcessTable($nbr,name)) failed at "
          append msg [clock format $ProcessTable($nbr,stopTime)  -format "%d.%m %H:%M"]
          Message::showMessage $msg $Info(remMsgColor)
        }
      
      #--Otherwise just blank state name
      } else {
        set ProcessTable($nbr,state) ""
      }
 
    }
    
  #---No more valid process
  } else {
    return
  }
}


proc Process::succesful {nbr} {
  global ProcessTable

  if { [info exists ProcessTable($nbr,logfile) ] &&
       [file exists $ProcessTable($nbr,logfile)]
     } {

    set ch_id [open $ProcessTable($nbr,logfile)]

    set line2 ""

    #-Read till the end
    while {![eof $ch_id]} {

      set count [gets $ch_id line]

      if { $count > 0 } {
        set line2 $line
      }

      if { -1 < [string first "ALL DONE" $line2] } {
        catch { close $ch_id} 
        return 1
      }
    }

    catch { close $ch_id} 
    
    return 0        
  }
 
  # No logfile found
  return 1
}


# Proc parses each inputline from solver output
#
proc Process::solverBrowseParser {nbr datalines } {
  global Equation Info ProcessTable SolverSystem

# NOTE: Catch is used in the case some the user defined solvers etc.
# are not handled correctly!!!
#
catch {
  set lines [split $datalines "\n"]

  foreach line $lines {
    set time_step 0
    set coupled_step 0
    set steady_step 0
    set err 0
    set warn 0
    set eq_step 0  ;# Any equation step
    set ns_step 0
    set ke_step 0
    set dc_step 0
    set st_step 0
    set ad_step 0
    set ud_step 0 ;# User defined proc
    set result_norm 0
    set relative_change 0
    set lin_iter_done 0
    set lin_iter_cycle 0
    set lin_iter_limit 0.0

    # Try to pick data
    # ================

    #-Errors
    if { -1 != [string first "ERROR:" $line] } { 
      
      set err 1
      set err_msg [string trim $line]

    #-Warnings
    } elseif { -1 != [string first "WARNING:" $line] } { 
      
      set warn 1
      set warn_msg [string trim $line]


    #-Transient time step message
    } elseif { -1 != [string first "Time:" $line] } {
      
      set time_step 1
      set current_step [lindex $line end-1]
      set current_time [lindex $line end]

    #-Transient coupled system iteration step message
    } elseif { -1 != [string first "Coupled system iteration:" $line] } {
      
      set coupled_step 1
      set current_step [lindex $line end]

    #-Steady state time message
    } elseif { -1 != [string first "Steady state iteration:" $line] } { 

      set steady_step 1
      set current_step [lindex $line end]

    #-Some  equation iteration (procedure name in the beginning!)
    #
    } elseif { [string match -nocase "* ITERATION *" $line] } { 
      set eq_step 1

      set proc_name [lindex [split $line ":"] 0]
      set proc_idx [lsearch $ProcessTable($nbr,solver,procNames) $proc_name]

      if { $proc_idx != -1 } {
        set idx [lindex $ProcessTable($nbr,solver,procIndices) $proc_idx]
        set sidx [lindex $ProcessTable($nbr,solver,procSubIndices) $proc_idx]
      
      # Procedure name for the equation is not yet stored in ProcessTable's list
      #
      # NOTE: For predefiend Elmer equations (DC, NS etc) we can use known iteration
      # message to make the link between equation and procedure
      #
      # For user defined procs we have to use the Procedure-field in Solver parameters
      # (this field is not available for predefiend equation which do not use
      #  dynamic libraries (dlls) for the implementation!
      
      #-Navier-Stokes iteration
      } elseif { [string equal -nocase "FlowSolve" $proc_name] } { 
        set ns_step 1
        set idx [DataField::getFieldProperty Equation NAVIER-STOKES EquationIndex]
        set sidx 0

      #-KE-Turbulence Kinetic Energy iteration
      } elseif { [string equal -nocase "KinEnergySolve" $proc_name] } { 
        set ke_step 1
        set idx [DataField::getFieldProperty Equation KE_TURBULENCE EquationIndex]
        set sidx 0

      #-KE-Turbulence Kinetic Energy Dissipation iteration
      } elseif { [string equal -nocase "KinDissipationSolve" $proc_name] } { 
        set ke_step 1
        set idx [DataField::getFieldProperty Equation KE_TURBULENCE EquationIndex]
        set sidx 1

      #-Temperature iteration
      } elseif { [string equal -nocase "HeatSolve" $proc_name] } { 
        set dc_step 1
        set idx [DataField::getFieldProperty Equation HEAT_EQUATION EquationIndex]
        set sidx 0

      #-StressAnalysis iteration
      } elseif { [string equal -nocase "StressSolve" $proc_name] } { 
        set st_step 1
        set idx [DataField::getFieldProperty Equation STRESS_ANALYSIS EquationIndex]
        set sidx 0

      #-AdvectionDiffusion iteration
      } elseif { [string equal -nocase "AdvectionDiffusion" $proc_name] } { 
        set ad_step 1
        set idx [DataField::getFieldProperty Equation ADVECTION_DIFFUSION_EQUATION EquationIndex]
        set sidx 0

      #-Some User Defined proc, use Solver data to get equation index
      } else {
        set ud_step
        set idx [Solver::getEquationIndexByProcedureName $proc_name]
        set sidx 0
      }
      
    #-Result norm
    } elseif { -1 != [string first "Result Norm" $line] } { 
      set result_norm 1
      set result_norm_val [Util::format "%.2e" [lindex $line end] ]

    #-Relative change
    } elseif { -1 != [string first "Relative Change" $line] } { 
      set relative_change 1
      set relative_change_val [Util::format "%.2e" [lindex $line end] ]

    #-Linear iteration cycle
    } elseif { [string is integer -strict [lindex $line 0] ] } { 
      set lin_iter_done 1
      set lin_iter_cycle [lindex $line 0]
      set lin_iter_limit [Util::format "%.2e" [lindex $line 1] ]
    }


    # Handle picked data
    # ==================

    # System Id for status info display
    #
    set id1 $ProcessTable($nbr,solver,resultIndex)
    set id2 $ProcessTable($nbr,solver,resultSubIndex)
    set id "$id1.$id2"

    #--Error
		if {$err} {
      Process::addErrorOrWarning $nbr $err_msg
      incr ProcessTable($nbr,solver,nofErrors) 

    #--Warning
    } elseif {$warn} {
      Process::addErrorOrWarning $nbr $warn_msg    
      incr ProcessTable($nbr,solver,nofWarnings) 

    #--Time step
    } elseif {$time_step} {
      
      # Remove possible "," separators
      set current_step [string trim $current_step ","] 
      set current_step [string trim $current_step " "] 

      set current_time [string trim $current_time ","] 
      set current_time [string trim $current_time " "] 

      # Store values
      set ProcessTable($nbr,currentStep1) $current_step
      set ProcessTable($nbr,currentTime) $current_time

      set ProcessTable($nbr,stopTimeDisplayed) "Step: "
      append ProcessTable($nbr,stopTimeDisplayed) [format "%6d" $current_step]

    #--Coupled system iteration
    } elseif {$coupled_step} {

      # Remove possible "," separators
      set current_step [string trim $current_step ","] 
      set current_step [string trim $current_step " "] 

      # Store value
      set ProcessTable($nbr,currentStep2) $current_step

    #--Steady state iteration
    } elseif {$steady_step} {

      # Remove possible "," separators
      set current_step [string trim $current_step ","] 
      set current_step [string trim $current_step " "] 

      # Store value
      # NOTE: "timestep" is always the same as steady state iteration!
      set ProcessTable($nbr,currentStep1) $current_step
      set ProcessTable($nbr,currentStep2) $current_step

    #---Result norm
    } elseif {$result_norm} {
        
      # If some system already created results
      #
      if { $ProcessTable($nbr,solver,resultIndex) != $Info(NO_INDEX) } {
        set ProcessTable($nbr,solver,$id,resultNorm) $result_norm_val
      }

    #---Relative change
    } elseif {$relative_change} {
        
      # If some system already created results
      #
      if { $ProcessTable($nbr,solver,resultIndex) != $Info(NO_INDEX) } {
        set ProcessTable($nbr,solver,$id,relativeChange) $relative_change_val
      }

    #---Iteration step info
    #
    } elseif { $eq_step } {

      #-Store first still unstored proc-info
      #
      if { $ns_step || $ke_step || $dc_step || $st_step || $ad_step || $ud_step } {
        lappend $ProcessTable($nbr,solver,procNames) $proc_name
        lappend $ProcessTable($nbr,solver,procIndices) $idx
        lappend $ProcessTable($nbr,solver,procSubIndices) $sidx
      }

      # Store current linearization iteration
      set ProcessTable($nbr,currentStep3) [lindex $line end]
      
      # Store current solver ids info
      # NOTE: activeId for status info display is constructed from these!!!###!!!
      #
      set ProcessTable($nbr,solver,resultIndex) $idx
      set ProcessTable($nbr,solver,resultSubIndex) $sidx
      
      # New id
      #
      set id1 $ProcessTable($nbr,solver,resultIndex)
      set id2 $ProcessTable($nbr,solver,resultSubIndex)
      set id "$id1.$id2"
			
			# If system not yet in the process table, add it to active ids
			#
			if { ![info exists ProcessTable($nbr,$id,outputSteps)] } {

				lappend ProcessTable($nbr,activeIds) $id

        set sys_names [Equation::getFieldPropertyByIndex $idx ProcessTableName]
        set sn [lindex $sys_names $sidx]
        set ProcessTable($nbr,$id,systemName) $sn
 
				set ProcessTable($nbr,$id,outputSteps) ""
				set ProcessTable($nbr,$id,outputData) ""
        set ProcessTable($nbr,$id,printStep1) 0
        set ProcessTable($nbr,$id,printStep2) 0
        set ProcessTable($nbr,$id,printStep3) 0

        set ProcessTable($nbr,$id,lastOutputIndex) -1
      }

      # If new linearization iteration cycle started, update coupled/steady state step
      # and mark new cycle by setting step3 to zero
      #
      if { $ProcessTable($nbr,currentStep3) == 1 } {
        set step1 $ProcessTable($nbr,currentStep1)
        set step2 $ProcessTable($nbr,currentStep2)

        lappend ProcessTable($nbr,$id,outputSteps) [list $step1 $step2 0] 
        lappend ProcessTable($nbr,$id,outputData) "" 
        set ProcessTable($nbr,activeSystem) $id
        set ProcessTable($nbr,solver,$id,linIterCycle) 0
        set ProcessTable($nbr,solver,$id,linIterLimit) 0.0

      }

    #---Linear iteration cycle info
    #
    } elseif {$lin_iter_done} {
      set idx $ProcessTable($nbr,solver,resultIndex)
      set sidx $ProcessTable($nbr,solver,resultSubIndex)
      set ProcessTable($nbr,solver,$id,linIterCycle) $lin_iter_cycle
      set ProcessTable($nbr,solver,$id,linIterLimit) $lin_iter_limit

    } else {
      ;# Nothing
    }

    # Update solved step info
    # =======================
    #
    if {$relative_change} {
      set idx $ProcessTable($nbr,solver,resultIndex)
      set sidx $ProcessTable($nbr,solver,resultSubIndex)

      # Steps
      set step1 $ProcessTable($nbr,currentStep1)
      set step2 $ProcessTable($nbr,currentStep2)
      set step3 $ProcessTable($nbr,currentStep3)
      
      # Step info: "result norm (relative change, nof iterations)"
      set si ""
      append si [format "%3d" $step2]
      append si ":"
      append si [format "%2d" $step3]
        
      set data1 $ProcessTable($nbr,solver,$id,resultNorm)
      set data2 $ProcessTable($nbr,solver,$id,relativeChange)

      set data3 $ProcessTable($nbr,solver,$id,linIterCycle)
      
      append si " "
      append si [Util::format "% 1.3e" $data1]
      append si "("
      append si [Util::format "%1.2e" $data2]
      append si [Util::format "%4d" $data3]
      append si ")"

      append si "  "

      # Update step numbers list
      lappend ProcessTable($nbr,$id,outputSteps) [list $step1 $step2 $step3]

      # Update step infos list
      lappend ProcessTable($nbr,$id,outputData) $si
      
      # We do not need too old step info, 3 latest steps should be enough
      # to understand to the state of iterations
      #
      set len [llength $ProcessTable($nbr,$id,outputSteps)]

      if {$len > 3} {
        set ProcessTable($nbr,$id,outputSteps) [lrange $ProcessTable($nbr,$id,outputSteps) end-2 end]
        set ProcessTable($nbr,$id,outputData)  [lrange $ProcessTable($nbr,$id,outputData) end-2 end]
      }
    }

  }

} ;# Catch whole solverBrowseParser proc

  return $datalines
}


# Add an errror or a warning line to the process info list
#
proc Process::addErrorOrWarning {nbr msg} {
  global ProcessTable

  set count [expr $ProcessTable($nbr,solver,nofErrors) + $ProcessTable($nbr,solver,nofWarnings)]

  # Max nof printed warnings is 20
  if { $count <= 20 } {
        
    # One more printed
    if { $count < 20 } {
      set txt  "$msg\n"
    # No more printed message (this would the 21.)
    } else { 
      set txt "Over 20 error or warnings! ...\n"
    }

    set ProcessTable($nbr,processInfo) "$ProcessTable($nbr,processInfo)$msg\n"

  }
}


# end ecif_tk_processTablePanel.tcl
# ********************
