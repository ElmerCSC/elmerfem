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
#Module:    ecif_tk_postFileSelectPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for selecting current ElmerPost file to view
#
#************************************************************************


proc PostFileSelect::openPanel {} {
  # Global variables
  global Info Model ModelProperty PostFileSelect

  set w $PostFileSelect(winName)
  set wgeom $PostFileSelect(winGeometry)

  set id [winfo atom $w]
  set PostFileSelect(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow PostFileSelect $id $PostFileSelect(winTitle) $wgeom] } {
    return
  }  

  set PostFileSelect(dataChanged) 0 
  set PostFileSelect(dataModified) 0

  toplevel $w
  focus $w
  set this $w

  wm title $w $PostFileSelect(winTitle)
  wm geometry $w $wgeom 

  PostFileSelect::updateFileList

  if { $PostFileSelect(postFile) == "" } {
    set PostFileSelect(postFile) [MenuExec::getSolverResultFile 0 POST_FILE]
  }

  if { $PostFileSelect(addDirectory) == "" } {
    set PostFileSelect(addDirectory) $ModelProperty(RESULTS_DIRECTORY,absolute)
  }

  # Backup vars
  foreach var {fileList postFile addDirectory} {
    set PostFileSelect($var,old) $PostFileSelect($var)
  }

  set PostFileSelect(fileIndex) ""
  set PostFileSelect(fileName) ""
  set PostFileSelect(fileEquations) ""
  set PostFileSelect(fileTime) ""
  set PostFileSelect(fileStatus) ""
  set PostFileSelect(fileComment) ""

  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  # Frames
  frame $w.f1    ;# Postfiles namelist and cmd button frame
  frame $w.f11   ;# Namelist frame
  frame $w.f12   ;# Cmd button frame
  frame $w.f2    ;# Postfile info frame
  frame $w.fB    ;# Buttons


  # Namelist and cmd buttons frames
  # ===============================

  pack $w.f11 $w.f12 -in $w.f1 -side left -padx $fpx3 -pady $fpy3

  # Namelist frame
  # --------------
  set m $w.f11

  set bg  $Info(nonActiveBg)
  
  set lb_frame [frame $m.lbf]

  set lb [listbox $lb_frame.lb  -height 10 -width 40 \
          -xscrollcommand [list $lb_frame.sb_x set]  \
          -yscrollcommand [list $lb_frame.sb_y set]  \
          -bg $bg ]
  set sb_x [scrollbar $lb_frame.sb_x -orient horizontal \
                                     -command [list $lb_frame.lb xview] ]
  set sb_y [scrollbar $lb_frame.sb_y -orient vertical \
                                     -command [list $lb_frame.lb yview] ]

  set PostFileSelect(fileNamesLB) $lb
  bind $PostFileSelect(fileNamesLB) <ButtonRelease-1> "+PostFileSelect::displayFileInfo"


  pack $lb_frame -side left -anchor e -expand 1 -fill both
  pack $sb_x -side bottom -fill x -expand 0
  pack $lb -side left -anchor e -fill both -expand 1 
  pack $sb_y -side right -fill y -expand 0 

  # Command buttons frame
  # ---------------------
  frame $w.f121
  frame $w.f122

  set m $w.f121
  set wdg [button $m.add  -text "Add..." -width 7 -command PostFileSelect::add]
  set wdg [button $m.drop -text "Drop"  -width 7 -command PostFileSelect::drop]
  pack $m.add $m.drop -side top -padx $fpx3 -pady $fpy3

  set m $w.f122
  set wdg [button $m.post -text "Open in ElmerPost" -command PostFileSelect::post]
  Widget::configureS $wdg disabled

  pack $m.post

  set PostFileSelect(postButton) $wdg

  pack $w.f121 -in $w.f12 -side top -pady $fpy3
  pack $w.f122 -in $w.f12 -side top -anchor s -pady $fpy3
 
  
  # File info frame
  # ===============
  
  set lwid 12
  set ewid 40

  # File name
  # ---------
  set m [frame $w.f2.name]
  frame $m.lf
  frame $m.ef

  #-Label
  label $m.lf.l -text "Post file:" -width $lwid
  pack $m.lf.l

  #-Entry
  entry $m.ef.e  -textvariable PostFileSelect(fileName) -width $ewid \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -padx $fpx1 -pady $fpy3

  pack $m.lf $m.ef -side left -anchor w -padx $fpx1 -pady $fpy1
  pack $m -side top -anchor w

  # File equations
  # --------------
  set m [frame $w.f2.equations]
  frame $m.lf
  frame $m.ef

  #-Label
  label $m.lf.l -text "Equations:" -width $lwid
  pack $m.lf.l

  #-Entry
  entry $m.ef.e  -textvariable PostFileSelect(fileEquations) -width $ewid \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -side top -padx $fpx1 -pady $fpy3

  pack $m.lf $m.ef -side left -padx $fpx1 -pady $fpy1
  pack $m -side top -anchor w


  # File time
  # ---------
  set m [frame $w.f2.time]
  frame $m.lf
  frame $m.ef

  #-Label
  label $m.lf.l -text "Modified time:" -width $lwid
  pack $m.lf.l

  #-Entry
  entry $m.ef.e  -textvariable PostFileSelect(fileTime) -width $ewid \
                 -font $Info(entryFont) -state disabled -bg $bg
  pack $m.ef.e -side top -padx $fpx1 -pady $fpy3

  pack $m.lf $m.ef -side left -padx $fpx1 -pady $fpy1
  pack $m -side top -anchor w

  # File run status
  # ---------------
  set m [frame $w.f2.status]
  frame $m.lf
  frame $m.ef

  #-Label
  label $m.lf.l -text "Last run status:" -width $lwid
  pack $m.lf.l

  #-Entry
  entry $m.ef.e  -textvariable PostFileSelect(fileStatus) -width $ewid \
                 -font $Info(entryFont) -state disabled -bg $bg

  pack $m.ef.e -side top -padx $fpx1 -pady $fpy3

  pack $m.lf $m.ef -side left -padx $fpx1 -pady $fpy1
  pack $m -side top -anchor w

  # Comment
  # -------
  set m [frame $w.f2.comment]
  frame $m.lf
  frame $m.ef

  #-Label
  label $m.lf.l -text "Comment:" -width $lwid
  pack $m.lf.l -anchor w -fill x

  #-Entry
  set wdg [entry $m.ef.e  -textvariable PostFileSelect(fileComment) -width $ewid \
                          -font $Info(entryFont)]
  bind  $wdg <KeyRelease> "+Panel::panelDataChanged 1 PostFileSelect $wdg {%A %K}"

  pack $m.ef.e -side top -padx $fpx1 -pady $fpy3

  pack $m.lf -side left -anchor w -padx $fpx1 -pady $fpy1
  pack $m.ef -side left -anchor w -padx $fpx1 -pady $fpy1
  pack $m -side top -anchor w


  # Create and pack buttons
  # =======================
  set m $w.fB

  #-Ok button

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $m.ok -text OK -command "PostFileSelect::panelOk $this"]
  set cn_btn [button $m.cancel -text Cancel -command "PostFileSelect::panelCancel $this" \
                               -state $ca]
  set ap_btn [button $m.apply -text Apply -command "PostFileSelect::panelApply $this" \
                              -state $ap]

  focus $ok_btn
  set PostFileSelect(applyButton)  $ap_btn
  set PostFileSelect(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1


  # Pack parent frames
  # ==================
  pack $w.f1 -side top -anchor w -fill x -expand 1 -padx $fpx3 -pady $fpy3
  pack $w.f2 -side top -anchor w -fill x -expand 1 -padx $fpx3 -pady $fpy3
  pack $w.fB -side top -padx $fpx3 -pady $fpy3


  # Init panel
  # ==========

  PostFileSelect::updateFileNamesLB
}



proc PostFileSelect::add {} {
  global PostFileSelect

  set types {
    { {ElmerPost} {.ep .Ep .EP} }
    { {All} {*} }
  }
  set title "Open ElmerPost File"

  #--Display dialog and get path
  set path [ MenuExec::openFile $title $types $PostFileSelect(addDirectory) ]

  if { $path == "" } {
    return
  }

  set PostFileSelect(addDirectory) [file dirname $path]

  set pfname $path

  # Check if file already in the list
  foreach finfo $PostFileSelect(fileList) {

    set fname [lindex $finfo 0]
    if { $fname == $pfname } {
      return
    }
  } 

  set pf_info [PostFileSelect::getDefaultInfoEntry $pfname]

  lappend PostFileSelect(fileList) $pf_info

  PostFileSelect::updateFileNamesLB

  Panel::panelDataChanged 1 PostFileSelect
}


proc PostFileSelect::drop {} {
  global PostFileSelect

  set index $PostFileSelect(fileIndex)

  set PostFileSelect(fileList) \
    [lreplace $PostFileSelect(fileList) $index $index]

  PostFileSelect::updateFileNamesLB

  Panel::panelDataChanged 1 PostFileSelect
}


proc PostFileSelect::post {} {
  global PostFileSelect Info

  set server_per_file 0

  set index $PostFileSelect(fileIndex)

  set finfo [lindex $PostFileSelect(fileList) $index]

  set fserver [PostFileSelect::getInfoItem $finfo SERVER]
  set fname [PostFileSelect::getInfoItem $finfo FNAME]
  set fcommand [PostFileSelect::getInfoItem $finfo COMMAND]
  
  #-NOTE: server per file!
  if { $server_per_file } {
    incr PostFileSelect(serverId)
    set fserver "ElmerPostS-$PostFileSelect(serverId)"

  #-NOTE: common server
  } else {
    set fserver $Info(Results,tclServer)
  }
    
  #--Check if post server exists
  set rc [MenuExec::prepare_ElmerPost_exec $fname $fserver fcommand]

  # Use existing server
  # ==================
  if { $rc == 1 } {
    MenuExec::sendTclCommands $fserver $fcommand
    Message::showMessage "$fserver updated: $fcommand" 
    return 
  }

  # Start new server
  # ================

  set finfo [PostFileSelect::setInfoItem $finfo SERVER $fserver]
  set finfo [PostFileSelect::setInfoItem $finfo COMMAND $fcommand]

  set PostFileSelect(fileList) \
    [lreplace $PostFileSelect(fileList) $index $index $finfo]

  MenuExec::prepare_elmer_exec Results [list $fname]
  return
}


proc PostFileSelect::panelOk {w {do_cancel 1}} {
  global PostFileSelect Model 

  if { $PostFileSelect(fileIndex) != "" } {
    
    #--Update comment
    set index $PostFileSelect(fileIndex)
    set finfo [lindex $PostFileSelect(fileList) $index]
    set finfo [PostFileSelect::setInfoItem $finfo COMMENT $PostFileSelect(fileComment)]
    set PostFileSelect(fileList) [lreplace $PostFileSelect(fileList) \
                                     $index $index $finfo]
  }

  #--Backup vars, must be done after comment update!
  #  Meaningul only for Apply
  #
  if { !$do_cancel } {
    foreach var {fileList postFile addDirectory} {
      set PostFileSelect($var,old) $PostFileSelect($var)
    }
  }

  Panel::panelDataChanged 0 PostFileSelect 

  # If we should close the panel
  if { $do_cancel } {
    Panel::cancel $w
  }
}


proc PostFileSelect::panelApply {w} {

  PostFileSelect::panelOk $w 0
}


proc PostFileSelect::panelCancel {w} {
  global PostFileSelect

  # Restore vars
  foreach var {fileList postFile addDirectory} {
    set PostFileSelect($var) $PostFileSelect($var,old)
  }

  Panel::cancel $w
}


proc PostFileSelect::displayFileInfo {} {
  global PostFileSelect Model ProcessTable

  set PostFileSelect(fileIndex) [$PostFileSelect(fileNamesLB) curselection]

  set index $PostFileSelect(fileIndex)

  if { $index == "" } {
    Widget::configureS $PostFileSelect(postButton) disabled
    return
  }

  set finfo [lindex $PostFileSelect(fileList) $index]

  set fname [PostFileSelect::getInfoItem $finfo FNAME]
  set feqns [PostFileSelect::getInfoItem $finfo EQUATIONS]
  set frnbr [PostFileSelect::getInfoItem $finfo RUNNUMBER]
  set fcomm [PostFileSelect::getInfoItem $finfo COMMENT]

  if { [info exists ProcessTable($frnbr,state)] } {
    set fstatus [string toupper $ProcessTable($frnbr,state)]
    append fstatus "  (process $frnbr)"
  } else {
    set fstatus "Unknown"
  }

  if { [file exists $fname] } {
    set ftime [Util::timeClock2Str [file mtime $fname]]
    set PostFileSelect(postFile) $fname
    Widget::configureS $PostFileSelect(postButton) normal
  } else {
    set ftime "File does not exist"
    set PostFileSelect(postFile) ""
    Widget::configureS $PostFileSelect(postButton) disabled
  }


  set fn [file tail $fname]
  set mn [lindex [file split [file dirname $fname] ] end]
  set dn [file join $mn $fn]

  set PostFileSelect(fileName) $dn
  set PostFileSelect(fileEquations) $feqns
  set PostFileSelect(fileTime) $ftime
  set PostFileSelect(fileStatus) $fstatus
  set PostFileSelect(fileComment) $fcomm
}


proc PostFileSelect::updateFileNamesLB {} {
  global PostFileSelect

  # Form file name list
  set PostFileSelect(fileNames) ""

  foreach finfo $PostFileSelect(fileList) {
    lappend PostFileSelect(fileNames) [PostFileSelect::getInfoItem $finfo FNAME]
  }

  $PostFileSelect(fileNamesLB) delete 0 end

  ListBox::fill $PostFileSelect(fileNamesLB) $PostFileSelect(fileNames)

  # Set current postfile selected
  set index [lsearch $PostFileSelect(fileNames) $PostFileSelect(postFile)]

  if { $index != -1 } {
    $PostFileSelect(fileNamesLB) selection set $index
  } else {
    $PostFileSelect(fileNamesLB) selection set 0
  }

  PostFileSelect::displayFileInfo
}



# Update result file list for ElmerPost
# If a solver runnumber given, add current postFile
# to the list
# Otherwise just update the list
# Removes all unaccessable (deleted) postfiles from
# the list
#
proc PostFileSelect::updateFileList { {runnumber ""} } {
  global PostFileSelect   Model ProcessTable

  set pf_names ""

  # Post files for the current solver run
  if { $runnumber != "" &&
       [info exists ProcessTable($runnumber,solver,postFiles)]
     } {
    set pf_names $ProcessTable($runnumber,solver,postFiles)
  }

  # Loop files and update their info in the file display list
  foreach pf_name $pf_names { 

    set equations [join [MenuExec::getSolverResultFileEquations $pf_name]]

    # Default entry
    set pf_info [PostFileSelect::getDefaultInfoEntry $pf_name $runnumber $equations]

    # Check if name already stored
    # ============================
    set index 0
    foreach finfo $PostFileSelect(fileList) {

      set fname [PostFileSelect::getInfoItem $finfo FNAME]

      # Remove existing entry from the current position and
      # insert in the beginning
      if { $pf_name != "" && $pf_name == $fname } {

        # Remove from the current position
        set PostFileSelect(fileList) \
            [lreplace $PostFileSelect(fileList) $index $index]

        # Copy current data
        set pf_info $finfo
               
        set pf_info [PostFileSelect::setInfoItem $pf_info RUNNUMBER $runnumber]
        set pf_info [PostFileSelect::setInfoItem $pf_info EQUATIONS $equations]
               
        break
      }

      incr index
    }

    # Insert in the beginning
    set PostFileSelect(fileList) [linsert $PostFileSelect(fileList) \
                                              0 $pf_info]

  }

  # Finally, check all files in the list
  # ====================================
  set index 0
  foreach finfo $PostFileSelect(fileList) {

    set fname [PostFileSelect::getInfoItem $finfo FNAME]
    set frnbr [PostFileSelect::getInfoItem $finfo RUNNUMBER]

    # If the file does not exist (any more) and a solver
    # is not running for it, remove file from the list  
    if { ![file exists $fname] &&
         ( ![info exists ProcessTable($frnbr,exists)] ||
           $ProcessTable($frnbr,state) != "running"
         )
      } {
      set PostFileSelect(fileList) \
          [lreplace $PostFileSelect(fileList) $index $index]

     incr index -1
    }

    incr index
  }
}



#################################
### INFO ENTRY HANDLING PROCS ###
#################################


proc PostFileSelect::getDefaultInfoEntry {
                        fname 
                        {runnumber ""}
                        {equations ""}
                        {server ""}
                        {command ""}
                        {time ""}
                        {run_status ""}
                        {comment ""}
                        } {

  return [list $fname $runnumber $equations $server $command  \
               $time $run_status $comment]
}


proc PostFileSelect::getInfoItem { info item_name } {
  global Info

  set idx [PostFileSelect::getInfoItemIndex $item_name]

  if { $idx == $Info(NO_INDEX) } {
    return ""
  } else {
    return [lindex $info $idx]
  }
}


proc PostFileSelect::setInfoItem { info item_name value} {
  global Info

  set idx [PostFileSelect::getInfoItemIndex $item_name]

  if { $idx == $Info(NO_INDEX) } {
    return $info
  } else {
    return [lreplace $info $idx $idx $value]
  }
}


proc PostFileSelect::getInfoItemIndex {item_name} {
  global Info

  set idx $Info(NO_INDEX)

  switch [string trim $item_name] {

  FNAME       { set idx 0 }
  RUNNUMBER   { set idx 1 }
  EQUATIONS   { set idx 2 }
  SERVER      { set idx 3 }
  COMMAND     { set idx 4 }
  TIME        { set idx 5 }
  RUNSTATUS   { set idx 6 }
  COMMENT     { set idx 7 }
  }

  return $idx
}


#End ecif_postFileSelectPanel.tcl

