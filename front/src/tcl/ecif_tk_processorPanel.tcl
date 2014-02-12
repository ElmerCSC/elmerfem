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
#Module:    ecif_tk_processorDefPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining parallel processing options
#
#************************************************************************
 
#------Parallel processing definitions  proc------


proc Processor::openPanel { } {
  ## This procedure displays the model level definition panel
  ## Global variables
  global Info Common Processor Solver SolverSystem

  set w $Processor(winName)
  set wgeom $Processor(winGeometry)

  set id [winfo atom $w]
  set Processor(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow Processor $id $Processor(winTitle) $wgeom] } {
    return
  }  

  set Processor(dataModified) 0
  set Processor(dataChanged) 0
  set Processor(datav) 0

  toplevel $w
  set this $w
  focus $w 

  wm title $w $Processor(winTitle)
  wm geometry $w $wgeom 

  Panel::initFields Processor
  Panel::backupFields Processor
  
  #-----WIDGET CREATION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)
  #-Data is organized into two frames.
  #-Nof processors related data
  frame $w.f1 -bd 2 -relief groove 
  #-Buttons
  frame $w.fB ;#--Buttons
  
  set m $w.f1
  frame $m.lf
  label $m.lf.l  -anchor nw -text "Settings for parallel processing:"
  pack $m.lf.l -side top

  #-Nof processors
  frame $m.ef
  label $m.ef.l -anchor w -text "Number of processors:"
  set wdg [entry $m.ef.e  -textvariable Processor(NOF_PROCESSORS) -width 10 \
                          -font $Info(entryFont)]

  #-If no system-equation defined, no data saving!
  if {$SolverSystem(nofActiveIndices) == 0} {
    $wdg configure -state disabled
  }

  set Processor(NOF_PROCESSORS,old) $Processor(NOF_PROCESSORS)
  set Processor(NOF_PROCESSORS,prev) $Processor(NOF_PROCESSORS)
  set Processor(NOF_PROCESSORS,err) 0
  set Processor(NOF_PROCESSORS,mod) 0

  bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 Processor $wdg {%A %K}"
  Widget::setEntryBindings non_standard Processor NOF_PROCESSORS $wdg

  pack $m.ef.l  $m.ef.e -side left -padx $fpx1
  pack $m.lf $m.ef -side top -padx $fpx1 -pady $fpy1 -expand 1

  #-Buttons

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "Processor::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "Processor::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command Processor::panelApply \
                                 -state $ap]

  focus $ok_btn
  set Processor(applyButton)  $ap_btn
  set Processor(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx1 
  
  pack $w.f1 $w.fB -side top  -expand 1 -padx $fpx2 -pady $fpy3
  
}


proc Processor::panelSave {} {
  global Processor Model

  #--Store old values
  Panel::backupFields Processor

  Panel::panelDataChanged 0 Processor 
  Panel::panelDataModified 0 Processor

  set Model(Front,needsUpdate) 1
  set Model(Solver,needsUpdate) 1

  #---Construct transfer data
  foreach var $Processor(allFields) {
    lappend data $Processor($var)
  }

  #--Write data into model
  Util::cpp_exec processorPanelOk $data
}


proc Processor::panelOk {w} {
  global Processor

  #---No changes
  if { !$Processor(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![Processor::checkPanelData] } {
    return
  }
  
  #---Save data
  Processor::panelSave

  Panel::cancel $w
}


proc Processor::panelApply {} {
  global Processor

  #---No changes
  if { !$Processor(dataChanged) } {
    return
  }

  if { ![Processor::checkPanelData] } {
    return
  }

  Processor::panelSave
}


proc Processor::panelCancel {w} {
  global Processor

  if { ![Panel::verifyCancel Processor] } {
    return
  }

  #---Reset into old values
  Panel::restoreFields Processor

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
# 
proc Processor::checkPanelData {} {
  global Info Processor Model

  #-Check nof processors data
  if { ![Processor::checkNofProcessors] } {
    return 0
  }

  # Data ok
  return 1
}


# Return 1 = ok, 0 = error
#
proc Processor::checkNofProcessors {} {
  global Info Processor Model

  set n $Processor(NOF_PROCESSORS)
  if { $n < 0 || $n > $Processor(max_nof_processors) } {
    set msg "Incorrect number of processors ( should be > 0 \n"
    append msg " and no larger than $Processor(max_nof_processors) "

    set Info(messageIcon) error
    Message::dispOkMessage $msg "$Info(FRONT_NAME) error message" $Processor(winName)
    return 0
  }

  return 1
}


# end ecif_tk_processorDefPanel.tcl
# ********************
