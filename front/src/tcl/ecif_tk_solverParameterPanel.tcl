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
#Module:    ecif_tk_solverDefPanel.tcl
#Language:  Tcl
#Date:      20.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the solver parameters
#
#************************************************************************
 

# This procedure displays the solver definition panel
#
proc Solver::openPanel { } {
  global Common Equation Info Model Solver SolverSystem
  upvar #0 Solver theArray

  set w $Solver(winName)
  set wgeom $Solver(winGeometry)

  set id [winfo atom $w]
  set Solver(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow Solver $id $Solver(winTitle) $wgeom] } {
    return
  }  

  #---DATA INITIALIZATION

  set Solver(dataChanged) 0
  set Solver(dataModified) 0

  if { ![info exists Solver(ids)] } {
    set Solver(ids) ""
  }
     
  #--Store old values
  foreach id $Solver(ids) {
    set Solver($id,data,old) $Solver($id,data)
    set Solver($id,name,old) $Solver($id,name)
  }

  MenuExec::updateRadiationStatus

  #-Select first system if no earlier selections
  if { ( ![info exists Solver(equationIndex)] ||
         $Solver(equationIndex) == -1         || 
         $Solver(equationIndex) == ""
       ) ||
       ( ![info exists Solver(subsystemIndex)] ||
         $Solver(subsystemIndex) == -1         || 
         $Solver(subsystemIndex) == ""
       )
     } {

    # Defaults
    #
    set Solver(equationIndex) 0
    set Solver(subsystemIndex) 0
    
    # If some indices, pick first of them
    #
    foreach index $SolverSystem(indices) {
    
      set Solver(equationIndex) $index
      set Solver(subsystemIndex) 0
      set Solver(systemId) [lindex $SolverSystem($index,ids) 0]
   
      break
    }
  }

  #-If we have no meshes, mark panel inactive
  #
  if { $Model(meshNames) == "" || $Solver(equationIndex) == -1 } {
    set Solver(panelIsActive) 0

  } else {
    set Solver(panelIsActive) 1

    set Solver(VARIABLE) [lindex $SolverSystem($Solver(equationIndex),subsystemNames) \
                           $Solver(subsystemIndex)]
  }

  #-NOTE: This is needed for StdPanelCreate::createValuesArea !!!
  #
  if { $Equation(activeIndices) != "" } {
    set Solver(equationIndices) $Equation(activeIndices)
  } else {
    set Solver(equationIndices) $Equation(indices)
  }

  set Solver(currentEquationIndices) $Solver(equationIndices)

  #-Check that current system is still active!
  #
  if { -1 == [lsearch $Solver(equationIndices) $Solver(equationIndex)] } {
    set Solver(equationIndex) [lindex $Solver(equationIndices) 0]
    set Solver(subsystemIndex) 0
  }

  set Solver(lastEquationIndex) ""
  set Solver(lastSubsystemIndex) ""
  set Solver(lastSystemId) ""
  
  #-Init values
  #
  set Solver(equationLabel) [Equation::getFieldPropertyByIndex $Solver(equationIndex) EquationLabel]
  set Solver(equationLabel,prev) $Solver(equationLabel)
  set Solver(valuesHdr) "$Solver(equationLabel) $Solver(valuesHdrBase)"

  Panel::initFields Solver "" 0 $Solver(equationIndex)

  if { !$Solver(panelSizeKnown) } {
    StdPanelCreate::calcPanelSize Solver
    set Solver(panelSizeKnown) 1
  }

  # Set nbr of the current equation
  set Solver(panelEqNbr) [expr 1 + [lsearch $Solver(equationIndices) $Solver(equationIndex)]]

  #---WIDGET CREATION
  toplevel $w
  focus $w 

  wm title $w $Solver(winTitle)
  wm geometry $w $wgeom
 
  Solver::createAndPackWidgets $w

  # Construct values area
  # =====================
  StdPanelCreate::constructPanelValuesArea Solver
  Solver::setSteadyCoupledLabel

  if { $Equation(activeIndices) == "" } {
    Solver::setFieldStates False $Solver(allFields)

  } else {
    Panel::backupFields Solver
  }

  Solver::updateFields
  Solver::updateVariablesMenu

  #--Set state of Prev/Next equation buttons
  #
  StdPanelExec::setEqMoveButtonsState Solver
  
  #--Set set of Prev/Next page buttons
  #
	set Solver(maxNofPanelPages) $Solver($Solver(equationIndex),nofPanelPages)
 
 	if { $Solver(panelPageNbr) > 1 } {
		$Solver(prevPageButton) configure -state normal
	}

	if { $Solver(panelPageNbr) < $Solver(maxNofPanelPages) } {
		$Solver(nextPageButton) configure -state normal
	}

  # Set field label bindings for right-button help
  #Widget::setLabelBindings Solver
}


proc Solver::createAndPackWidgets { w } {
  global Info Model Solver

  set this $w

  set fpx1 $Info(framePadX1)
  set fpx2 $Info(framePadX2)
  set fpx3 $Info(framePadX3)
  set fpy1 $Info(framePadY1)
  set fpy2 $Info(framePadY2)
  set fpy3 $Info(framePadY3)

  #-Data is organized into four frames.

  # OptionMenu for Equations
  # ========================
  frame $w.f1

  #-Label
  label $w.f1.l -text "Equation:" -width 12 -justify left
  pack $w.f1.l -side left -anchor w -fill x

  #-Option menu
  set eq_labels ""
  foreach idx $Solver(equationIndices) {
    lappend eq_labels [Equation::getFieldPropertyByIndex $idx EquationLabel]
  }

  set wdg $w.f1.om
  set om [Panel::createOptionMenuWidget $w.f1.om Solver(equationLabel) $eq_labels]
  $wdg configure -indicatoron 0 -width 30 
  $wdg configure -disabledforeground black

  set cmd "PanelCheck::checkOptionMenu $om $wdg Solver equationLabel "
  # Append check-proc argument to the check-call
  append cmd "\"Solver::equationMenuSelectionProc\"  {%s}"
  bind $om <<MenuSelect>> $cmd

  set Solver(solverEquationWidget) $wdg
  pack $wdg -side left -anchor w -fill x

  #-Next/Previous eqiation buttons
  set m $w.f1

  set prev_btn [button $m.prev -text "<<" -command "Solver::panelEqMove -1 " \
                                     -state disabled -disabledforeground $Info(nonActiveFg) ]

  set next_btn [button $m.next -text ">>" -command "Solver::panelEqMove 1" \
                                      -state disabled -disabledforeground $Info(nonActiveFg) ]

  set Solver(prevEqButton) $prev_btn
  set Solver(nextEqButton) $next_btn

  pack $prev_btn -side left -anchor e -fill x -padx 10
  pack $next_btn -side left -anchor w -fill x 


  # Variable and mesh option menus + update button
  # ==============================================
  frame $w.f2       ;#-bd 2 -relief groove
  frame $w.f2.f1    ;#-Variable + Mesh
  frame $w.f2.f1.fV ;#-Variable label + Option Menu
  frame $w.f2.f1.fM ;#-Mesh label + Option Menu
  frame $w.f2.f2    ;#-Update button

  #-Variable label text
  label $w.f2.f1.fV.l -text "Variable:" -width 12 -justify left
  pack $w.f2.f1.fV.l -side left -anchor w -fill x

  #-Variable option menu
  set wdg $w.f2.f1.fV.om
  set om [Panel::createOptionMenuWidget $w.f2.f1.fV.om Solver(VARIABLE) ""]
  $wdg configure -indicatoron 0 -width 30 
  $wdg configure -state disabled
  $wdg configure -disabledforeground black
  bind $om <<MenuSelect>> "Solver::setSubsystemIndex"
  set Solver(systemVariableWidget) $wdg
  pack $wdg -side left -anchor w -fill x

  #-Mesh label text
  label $w.f2.f1.fM.l -text "Mesh:" -width 12 -justify left
  pack $w.f2.f1.fM.l -side left -anchor w -fill x

  #-Mesh option menu
  set wdg $w.f2.f1.fM.om
  set om [Panel::createOptionMenuWidget $w.f2.f1.fM.om Solver(MESH) $Model(meshNames)]
  $wdg configure -indicatoron 0 -width 30 
  #$wdg configure -state disabled
  $wdg configure -disabledforeground black

  if { $Solver(MESH) != "" && [llength $Model(meshNames)] < 2 } {
    #$wdg configure -relief groove
    #$wdg configure -state disabled
  }

  bind $om <<MenuSelect>> "PanelCheck::checkOptionMenu $om $wdg Solver MESH \"\" {%s}"
  set Solver(solverMeshWidget) $wdg
  pack $wdg -side left -anchor w -fill x
  
  #-Update button
  set wdg [button $w.f2.f2.update -text "Update" -command Solver::update]
  pack $wdg -side top -anchor w -padx $fpx3
  $wdg configure -state disabled 
  set Solver(panelUpdateButton) $wdg


  # Solver values area
  # ==================
  frame $w.f3    ;#-bd 2 -relief groove 

  #-Create Values frame
  set wdg [frame $w.f3.fValues -relief groove -bd 2]
  set Solver(valuesFrame) $wdg
  set Solver(valuesParentFrame) $w.f3

  # Fixed sized for values area if multipanel!
  if { $Solver(panelsByEquation) } {
    $wdg configure -width $Solver(valuesFrame,width)
    $wdg configure -height $Solver(valuesFrame,height)
    pack propagate $wdg 0
  }

  # Buttons
  # =======
  frame $w.fS ;#--Status info
  frame $w.fB ;#--Buttons


  #-----WIDGET PACKING

  # Equation option menu frame
  # ==========================
  set m $w.f1
  pack $m -side top -anchor w -fill x -expand 1 -padx $fpx2 -pady $fpy2 

  # Variable and Mesh OptionMenus + Update button
  # =============================
  pack $w.f2 -side top -anchor nw -fill x -expand 1 -padx $fpx2 -pady $fpy2 
  pack $w.f2.f1 -side left -fill x -anchor w
  pack $w.f2.f2 -side right -padx $fpx3
  pack $w.f2.f1.fV $w.f2.f1.fM -side top -fill x -anchor w -pady $fpy2

  # Values area
  # ===========
  pack $w.f3 -side top -anchor w -fill both -expand 1 -padx $fpx1 -pady $fpy1
  pack $w.f3.fValues -fill both -expand 1 -anchor w -padx $fpx1 -pady $fpy1

  # Buttons (create and pack)
  # =========================

	# Info in the "status" bar
	#
  label $w.fS.l -text " " -width 20 -height 3 -justify left
  set Solver(entryInfoLabel) $w.fS.l


  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "Solver::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "Solver::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command Solver::panelApply  \
                                 -state $ap]

  set prev_btn [button $w.fB.prev -text "<<" -command "StdPanelExec::panelPageMove Solver -1" \
                                  -state disabled -disabledforeground $Info(nonActiveFg) ]

  set next_btn [button $w.fB.next -text ">>" -command "StdPanelExec::panelPageMove Solver 1" \
                                  -state disabled -disabledforeground $Info(nonActiveFg) ]
 
  focus $ok_btn
  set Solver(applyButton)  $ap_btn
  set Solver(cancelButton) $cn_btn
	set Solver(prevPageButton) $prev_btn
	set Solver(nextPageButton) $next_btn

  #pack $w.fS.l -side left -anchor w -expand 1 -fill x
  pack $w.fS.l -side left -anchor w 
  
  pack $next_btn $prev_btn -side right -anchor e -fill x -padx $fpx1 
  pack $ok_btn $cn_btn $ap_btn -side left -anchor w -padx $fpx1 

  pack $w.fS $w.fB -side left -anchor w -expand 1
  #pack $w.fB -side top  -padx $fpx3 -pady $fpy3
}


proc Solver::setSteadyCoupledLabel {} {
  global Solver

  if { ![info exists Solver(allWidgets,STEADY_STATE_CONVERGENCE_TOLERANCE)] ||
       ![winfo exists $Solver(allWidgets,STEADY_STATE_CONVERGENCE_TOLERANCE)]
     } {
    return
  }

  set wdg $Solver(allWidgets,label,STEADY_STATE_CONVERGENCE_TOLERANCE)

  if { [Timestep::isSteadyState] } {
    $wdg configure -text "Steady state tolerance:"

  } else {
    $wdg configure -text "Coupled eq. tolerance:"
  }

}


# Default value for a solver mesh, if
# no mesh defined yet
#
proc Solver::getDefaultMesh {} {
  global Model

  # Pick first of the defined meshes (if any)
  return [lindex $Model(meshNames) 0]
}


proc Solver::acceptEquationChange {} {
  global Equation Solver

  if { $Solver(dataModified) } {
    if { ![Panel::verifyParameter Solver] } {
      return 0
    } else {
      Panel::panelDataModified 0 Solver
    }
  }

  return 1
}


proc Solver::equationMenuSelectionProc {} {
  global Equation Solver
  
  # Nothing to do!
  #
  if { $Solver(equationLabel) == $Solver(equationLabel,prev) } {
    return
  }

  # Check what to do with modified data
  #
  if { $Solver(dataModified) } {
    if { ![Panel::verifyParameter Solver] } {
      set Solver(equationLabel) $Solver(equationLabel,prev)
      return
    } else {
      Panel::panelDataModified 0 Solver
    }
  }

  set Solver(targetMask) [Equation::getEquationMask $Solver(equationIndex)]

  Solver::updateFields

	# Set button states
  StdPanelExec::setEqMoveButtonsState Solver
}


# Read current panel data into a Solver parameter using
# "id" as the key
#
proc Solver::updateParameter {id} {
  global Solver

  if { [info exist Solver(equationIndex)] } {
    set sys_index $Solver(equationIndex)
  } else {
    set sys_index 0
  }

  set Solver($id,data) [DataField::formParameterLine Solver $sys_index] 
  Panel::backupFields Solver

  StdPanelExec::setValuesAreaStatus Solver ""
}


# Procedure moves to next or previous Panel equation
# direction: +-1
#
# NOTE: We do not update the Solver panel via StdPanelExec::panelEqMove, because
# StdPanelCreate::createValuesArea (called in StdPanelExec::panelEqMove) cannot be called
# in the 'standard' way for Solver panel (owing to 'Variables' field etc.).
# It is called in Solver::updateFields!!!
#
proc Solver::panelEqMove {direction} {
  global Solver
  
  # Check what to do with modified data
  #
  if { $Solver(dataModified) } {
    if { ![Panel::verifyParameter Solver] } {
      return
    } else {
      Panel::panelDataModified 0 Solver
    }
  }

  # NOTE: Solver panel is not yet update here!
  #
  set update 0
  StdPanelExec::panelEqMove Solver $direction $update

  set Solver(targetMask) [Equation::getEquationMask $Solver(equationIndex)]

  # This will update values area!
  #
  Solver::updateFields $Solver(equationIndex) 
}


# Update field values
#
# Activated e.g when equation menu is used
#
# If equation index is not given as an argument, check current equation
# by the equation label (updated typically by the equation menu command!)
#
# NOTE: Equation index  is needed when using Prev/Next button, because these
# update first the equation index (+/- 1) and panel data must be written after that!!!
#
proc Solver::updateFields { { eq_index "" } } {
  global Equation Solver SolverSystem Info

  # If equation index is not given, check which equation is currently
  # on the panel by using the equation label
  #
  if { $eq_index == "" } {
    foreach fld $Equation(allFields) {
    
      if { ![Equation::isEquation $fld] } {
        continue
      }

      if { $Solver(equationLabel) == [DataField::getFieldProperty Equation $fld EquationLabel "" 0] } {
        set Solver(equationIndex) [DataField::getFieldProperty Equation $fld EquationIndex]
        break;
      }
    }
  
  # If equation index given, set new label
  #
  } else {
    set Solver(equationLabel) [Equation::getFieldPropertyByIndex $eq_index EquationLabel]
  }


  set SN_PREV  $Solver(lastEquationIndex)
  set SSN_PREV $Solver(lastSubsystemIndex)
  set sid_PREV $Solver(lastSystemId)

  set SN $Solver(equationIndex)

  # If main system changes, reset subsystem values
  #
  if { $SN != $SN_PREV } {
    set Solver(subsystemIndex) 0
    set Solver(VARIABLE) [lindex $SolverSystem($SN,subsystemNames) 0 ]
    set Solver(valuesHdr) "$Solver(equationLabel) $Solver(valuesHdrBase)"
  }

  # Update system id
  #
  set SSN $Solver(subsystemIndex)
  set Solver(systemId) [lindex $SolverSystem($SN,ids) $SSN]

  set sid $Solver(systemId)


  # Rebuild values area when system is changed
  # ==========================================
  if { $sid_PREV != $sid } {

    Panel::resetFields Solver

    Panel::initFields Solver "" 0 $Solver(equationIndex)

    #-Set field values
    DataField::formDataFields Solver $Solver($sid,data)
    
    #-If different equation type
    #
    if { $SN != $SN_PREV } {
      StdPanelCreate::constructPanelValuesArea Solver $Solver(equationIndex)

    } else {
      PanelCheck::execPanelFillProcs Solver
      StdPanelExec::setValuesAreaActivity Solver $Solver(targetMask)
    }

    #-If different equation type
    if { $SN != $SN_PREV } {
      #StdPanelCreate::packValuesArea $Solver(valuesFrame) Solver
    }

    Panel::backupFields Solver

  # Or check current data
  # ====================
  } else {
    if { ![Solver::checkData] } {
      return
    }
  }

  # Set widget states
  # =================
  foreach fld $Solver(allFields) {
    Widget::configureField Solver $fld
  }

  if { $sid_PREV == $sid } {
    return
  }

  # Update last system indicators
  #
  set Solver(lastEquationIndex) $SN
  set Solver(lastSubsystemIndex) $SSN
  set Solver(lastSystemId) $sid

  set Solver(equationLabel,prev) $Solver(equationLabel)

  # If no need to update variables' or mesh name
  if { $SN == $SN_PREV } {
    return
  }

  # Update variables menu
  #
  Solver::updateVariablesMenu

  # Check mesh name
  #
  if { $Solver(MESH) == "" } {
    set Solver(MESH) [Solver::getDefaultMesh]
  }
}


proc Solver::updateVariablesMenu {} {
  global Solver SolverSystem

  set SN  $Solver(equationIndex)
  set SSN $Solver(subsystemIndex)

  # Update variables option menu
  # ============================
  set names ""
  set indices ""
  set counter 0
  foreach act $SolverSystem($SN,actives) \
          nm  $SolverSystem($SN,subsystemNames) {

    if {$act} {
      lappend names $nm
      lappend indices $counter
    }

    incr counter
  }

  if { $indices == "" } {
    set indices 0
  }

  set index [lsearch $indices $SSN]

  if { $index == -1 } {
    set index 0
  }

  set Solver(VARIABLE) [lindex $names $index]
  set Solver(subsystemIndex) [lindex $indices $index]

  set wdg $Solver(systemVariableWidget)
  destroy $wdg

  set om [Panel::createOptionMenuWidget $wdg Solver(VARIABLE) $names]
  $wdg configure -indicatoron 0 -width 20
  $wdg configure -disabledforeground black

  if { [llength $names] < 2 } {
    $wdg configure -relief groove
    $wdg configure -state disabled
  }

  bind $om <<MenuSelect>> "Solver::setSubsystemIndex"
  set Solver(systemVariableWidget) $wdg

  pack $wdg
}


# Solver specific values area clean proc
#
proc Solver::destroyValuesArea {} {
  global Solver
  
  # NOTE: We use this method for Solver, because
  # the standard method causes noticeable flickering
  # (Most of the field are always kept, and flickering is very
  # clear in this case)
  #
  foreach fld $Solver(allFields) {

    if { 1 || !$Solver($fld,dsp) } {

      # We destroy each fields frame
      #
      if { [info exist Solver($fld,frame)] &&
           [winfo exist $Solver($fld,frame)] } {

        destroy $Solver($fld,frame)
      }
    }
  }
}


# This is called when the variables OptionMenu is used
#
proc Solver::setSubsystemIndex {} {
  global Solver SolverSystem

  # Let option menu widget to update
  Util::doWait 50

  set SN $Solver(equationIndex)
  set names $SolverSystem($SN,subsystemNames)

  set index [ List::nocaseSearch $names $Solver(VARIABLE) ]

  if { $index == -1 } {
    return
  }

  set Solver(subsystemIndex) $index
 
  if { $Solver(subsystemIndex) != $Solver(lastSubsystemIndex) } {
    Solver::updateFields    
  }
}


# Procedure sets entry states and values according to the state
# of the control-value variable
# Value "True" is considered active value
# 
proc Solver::setFieldStates {control_value fields } {
  global Solver 

  if { $control_value == "True" } {
    set active 1
  } else {
    set active 0
  }
  
  # If Solver panel is not open we do not know
  # the active equation here!
  #
  if { ![info exists Solver(winName)] ||
       ![winfo exist $Solver(winName)]
     } {
    
    foreach fld $fields {
      set Solver($fld,act) $active
    }

    return
  }

  set SN $Solver(equationIndex)

  foreach fld $fields {

    # If state is activated set into initial value
    if {$active} {
      Panel::initField Solver $fld 0 $SN
      set Solver($fld,act) 1

    } else {
      set Solver($fld,act) 0
    }

  }

  foreach fld $fields {
    Widget::configureField Solver $fld
  }
}


# Update button proc
#
proc Solver::update {} {
  global Solver

  if { !$Solver(dataModified) } {
    return
  }

  if { ![Solver::checkPanelData] } {
    return
  }

  Solver::updateParameter $Solver(systemId)

  Panel::panelDataModified 0 Solver
  Panel::panelDataChanged 1 Solver
}


proc Solver::panelSave { {inform_front 1} } {
  global Info Solver Model

  #--Store old values
  foreach id $Solver(ids) {
    set Solver($id,data,old) $Solver($id,data)
    set Solver($id,name,old) $Solver($id,name)
  }

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  set Model(Solver,needsUpdate) 1

  Panel::panelDataChanged 0 Solver 
  Panel::panelDataModified 0 Solver

  #-Store model (file) data in compressed mode
  set sys_index 0
  set modified 0
  foreach id $Solver(ids) {
    set tmp [DataField::compressParameterOutputData Solver $id modified $sys_index]
    set Solver($id,data,uncompressed) $Solver($id,data)
    set Solver($id,data) $tmp
  }

  Solver::updateActiveMeshInfo
  Solver::updateActiveSolvingOrderList

  # Update (multigrid related part of the) object masks!
  StdPanelInit::createObjectTableMasks

  Util::cpp_exec solverParameterPanelOk
}


proc Solver::panelOk {w} {
  global Solver

  #---Check that all changes are applied
  if { ![Panel::verifySave Solver] } {
    return
  }

  #---No changes
  if { !$Solver(dataChanged) } {
    Panel::cancel $w; return
  }

  if {$Solver(equationIndex) == ""} {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![Solver::checkPanelData] } {
    return
  }

  # Update current solver
  Solver::updateParameter $Solver(systemId)

  #---Save data
  Solver::panelSave
  Panel::cancel $w
}


proc Solver::panelApply {} {
  global Info Solver

  #---Check that all changes are applied
  if { ![Panel::verifySave Solver] } {
    return
  }

  #---No changes
  if { !$Solver(dataChanged) } {
    return
  }

  if {$Solver(equationIndex) == ""} {
    return
  }

  #---Error in data
  if {  ![Solver::checkPanelData] } {
    return
  }

  # Update current solver
  Solver::updateParameter $Solver(systemId)
   
  #---Save data
  Solver::panelSave
}


proc Solver::panelCancel {w} {
  global Info Solver

  if { ![Panel::verifyCancel Solver] } {
    return
  }

  #---Reset into old values
  foreach id $Solver(ids) {
    set Solver($id,data) $Solver($id,data,old)
    set Solver($id,name) $Solver($id,name,old)
  }

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc Solver::checkPanelData {} {
  global Solver

  #-Check solver number fields
  if { ![Solver::checkData] } {
    StdPanelExec::setValuesAreaStatus Solver
    return 0
  }

  return 1
}


# Return 1 = ok, 0 = error
#
proc Solver::checkData {} {
  global Info Solver Model

  set totmsg ""

  if { [Util::nce $Solver(LINEAR_SYSTEM_SOLVER) "Multigrid"] } {
  
    if { $Solver(MG_SMOOTHER) == "" } {
      set Solver(MG_SMOOTHER) $Solver(LINEAR_SYSTEM_ITERATIVE_METHOD)
    }

    if { $Solver(MG_MAX_ITERATIONS) == "" } {
      set Solver(MG_MAX_ITERATIONS) $Solver(LINEAR_SYSTEM_MAX_ITERATIONS)
    }

    if { $Solver(MG_CONVERGENCE_TOLERANCE) == "" } {
      set Solver(MG_CONVERGENCE_TOLERANCE) $Solver(LINEAR_SYSTEM_CONVERGENCE_TOLERANCE)
    }
  }

  #-Check numeric variables
  foreach fld $Solver(allFields) {

    # Check only current target fields
    # (All fields for Solver in fact!)
    if { !$Solver($fld,act) } {
      continue
    }
    
    if { ![info exist Solver($fld)] } {
      continue
    }

    set msg [PanelCheck::checkFieldValue Solver $fld]

    #-Error messages are collected
    if {$msg != ""} {
 
      set Solver($fld,err) 1

      set label [DataField::getFieldProperty Solver $fld Label]

      lappend totmsg "\n"
      lappend totmsg [string toupper $fld]
      lappend totmsg "\n$msg\n"

    #-Field Ok
    } else {
      set Solver($fld,err) 0
    }
  }

  #-Ok (no msg-variable was created!)
  if { $totmsg == "" } {
    return 1

  #-Error message is displayed
  } else {
    set totmsg [linsert $totmsg 0 $Info(fieldMsg)]

    set Info(messageIcon) error
    Message::dispOkMessage $totmsg  "$Info(FRONT_NAME) message!" $Solver(winName)
    return 0
  }

  return 1
}


#====================#
# Solver helper proc #
#====================#

# Get equation sif name
#
proc Solver::getEquationSifName { system_index } {
  
  return [Equation::getFieldPropertyByIndex $system_index SifName]
}


# Get solver variable's name
#
proc Solver::getEquationVariable { system_index } {

  set names [Equation::getAllVariableNames $system_index]

  # If only one name defined, use that name, otherwise return
  # empty result to indicate that name cannot be resolved
  if { 1 == [llength $names] } {
    return $names

  } else {
    return ""
  }
}


# Get solver variable's degrees of freedom
# Values: prefined list of integer values
# - one value  --> fixed
# - three values --> 1D, 2D and 3D values
#
proc Solver::getEquationVariableDofs { system_index } {
  global Model

  set values_list [Equation::getFieldPropertyByIndex $system_index VariableDofs "" 0]
  set values_list [split $values_list ";"]

  set tot_dof 0

  foreach values $values_list {

    #-No dimension dependence
    #
    if { 1 == [llength $values] } {
      set dof $values

    #-Depends on simulation dimension (like NS and SA)
    #
    } elseif { $Model(SIMULATION_DIMENSION) == "1D" } {
      set dof [lindex $values 0]

    } elseif { $Model(SIMULATION_DIMENSION) == "2D" } {
      set dof [lindex $values 1]

    } else {
      set dof [lindex $values end]
    }

    incr tot_dof $dof
  }

  if { $tot_dof > 0 } {
    return $tot_dof

  } else {
    return ""
  }

}


proc Solver::getMeshIndex {solver_id} {
  global Model Solver

  set mn [DataField::getFieldValue Solver $solver_id MESH]

  return [lsearch $Model(meshNames) $mn]
}


# Calculate nof active meshes (ie. nof different mesh names
# in use in active solvers
# Update globals "Model(activeMeshIndices), Model(activeMeshNames), Model(nofActiveMeshes)"
#
proc Solver::updateActiveMeshInfo {} {
  global Info Model Solver

  set Model(nofActiveMeshes) 0
  set Model(activeMeshIndices) ""
  set Model(activeMeshNames) ""

  # If no meshes available!
  if { $Model(meshNames) == "" } {

    # Update model
    Util::cpp_exec readActiveMeshIndices
    Interface::setStatusMeshes

    return
  }

  set mesh_indices ""
  set mesh_names ""

  set using_current_mesh 0
  set implicitely_using_current_mesh 0

  foreach sid $Solver(ids) {

    set is_active [string tolower [DataField::getFieldValue Solver $sid ACTIVE]]

    if  { $is_active != "true" } {
      continue
    }

    set mname [DataField::getFieldValue Solver $sid MESH]
    set mindex [lsearch $Model(meshNames) $mname]

    if { $mname != "" && $mindex != - 1 && -1 == [lsearch $mesh_names $mname] } {

      lappend mesh_indices $mindex
      lappend mesh_names $mname

      if { $mindex == $Model(currentMeshIndex) } {
        set using_current_mesh 1
      }

    } else {
      set implicitely_using_current_mesh 1
    }
  }

  set Model(nofActiveMeshes) [llength $mesh_names]
  set Model(activeMeshIndices) $mesh_indices
  set Model(activeMeshNames) $mesh_names

  if { $Model(currentMeshIndex) != $Info(NO_INDEX) &&
       $implicitely_using_current_mesh &&
       !$using_current_mesh
     } {
    incr Model(nofActiveMeshes)
    lappend Model(activeMeshIndices) $Model(currentMeshIndex)
    lappend Model(activeMeshNames) [lindex $Model(meshNames) $Model(currentMeshIndex)]
  }

  # Update model
  Util::cpp_exec readActiveMeshIndices
  Interface::setStatusMeshes
}


proc Solver::setLinearizationType {} {
}

proc Solver::setStabilizationType {} {
}


#============================#
# SolverSystem related procs #
#============================#

# Update all solver related data based on all active systems
# Argument "data_modified" is a reference variable to indicate if
# data was modified
#
proc Solver::checkSolverData { {active_system_indices ""} {data_modified ""}} {
  global Equation Solver SolverSystem

  if { $data_modified != "" } {
    upvar $data_modified modified
  }

  set modified 0

  # If no solvers (equations) yet
  #
  if { ![info exists SolverSystem(indices)] ||
       $SolverSystem(indices) == ""
     } {
    return
  }

  if { $active_system_indices == "" } {
    set active_system_indices $Equation(activeIndices)
  }
         
  #-Build first system-info table based on current
  # equations
  set system_infos ""

  # Loop all (main) solver systems
  # ------------------------------
  foreach SN $SolverSystem(indices) {

    #-If no active equation for the system
    if { -1 == [lsearch $active_system_indices $SN] } {
      continue
    }

    #-The equation name being processed
    set sys_name $SolverSystem($SN,outName)

    # NOTE: Equation need not be True valued (active)
    # Only those solvers we are not related to any existing
    # variables (like Oxygen, if it was deleted from the list
    # of AdvDiff variables in EquationVariable) will be removed
    # Other non-active solvers are kept, for the case that they
    # would be activated later by activating the corresponding equation
    #
    set var_infos [Equation::findSystemNames $SN $Equation(ids)]

    #-Subsystem (variable) infos
    #
    foreach row $var_infos {

      set var_name [lindex $row 0]
      set active_flaq [lindex $row 1]
      
      # Creatae new system-info row
      lappend system_infos [list $sys_name $var_name $active_flaq] 

    }

  }

  #-Update solvers and compress ids
  Solver::updateSolvers $system_infos modified

  #-Form SolverSystem data based on updated solvers
  Solver::formSolverSystemData
  
  #-Update solving order
  SolverOrder::updateSolversAndCalculators  

  #-Save data
  if {$modified} {
    Solver::panelSave

  } else {
    Solver::updateActiveSolvingOrderList
  }
}


proc Solver::findSolverIdsForEquation {equation_id} {
  global Info Equation Solver SolverSystem

  set sids ""

  foreach sn $SolverSystem(activeIndices) {

    set sub_names [Equation::findSystemNames $sn $equation_id]

    set index 0

    foreach sub_name $sub_names {
      set active [lindex $sub_name 1]

      if { $active } {
        lappend sids [lindex $SolverSystem($sn,ids) $index]
      }

      incr index
    }
  }

  return [List::getUniqueIdList $sids]
}


# Forms (updates) SolverSystem-data by looping all Solver parameters
#
proc Solver::formSolverSystemData {} {
  global Equation Solver SolverSystem Info

  Solver::initSolverSystemData

  foreach id $Solver(ids) {
    
    set active [DataField::getFieldValue Solver $id ACTIVE]
    set name [DataField::getFieldValue Solver $id EQUATION]
    set subname [DataField::getFieldValue Solver $id VARIABLE]

    # Find system by the name and update data
    foreach index $SolverSystem(indices) {

      if { $SolverSystem($index,outName) == $name } {

        lappend SolverSystem($index,ids) $id
        lappend SolverSystem($index,actives) $active
        lappend SolverSystem($index,subsystemNames) $subname

        if {$active} {

          lappend Solver(activeIds) $id

          if { -1 == [lsearch $SolverSystem(activeIndices) $index] } {
            lappend SolverSystem(activeIndices) $index
            incr SolverSystem(nofActiveIndices)
          }
        }

        break
      }
    }
  }
}


# Initialize SolverSystem table with predefined data
#
proc Solver::initSolverSystemData {} {
  global Equation Solver SolverSystem

  set SolverSystem(indices) ""
  set SolverSystem(activeIndices) ""
  set SolverSystem(nofActiveIndices) 0

  set Solver(activeIds) ""

  foreach fld $Equation(allFields) {

    if { ![Equation::isEquation $fld] } {
      continue
    }
       
    set index [DataField::getFieldProperty Equation $fld EquationIndex]

    lappend SolverSystem(indices) $index

    set SolverSystem($index,name)  $fld
    set SolverSystem($index,label) [DataField::getFieldProperty Equation $fld EquationLabel]
    set SolverSystem($index,outName)  [DataField::getFieldProperty Equation $fld SifName]
    set SolverSystem($index,subsystemNames) ""
    set SolverSystem($index,ids) ""     ;#NOTE: These are solver parameter ids : Solver(1,...) etc.
    set SolverSystem($index,actives) "" ;#NOTE: Which of the ids are active
  }
}



# Update solver parameters (Solver(id,..) based on
# system-infos which are created from equations
# Argument "data_modified" is a reference variable to indicate
# if data was modified
#
proc Solver::updateSolvers {system_infos data_modified} {
  global Equation Info Solver SolverSystem
  upvar $data_modified modified

  set modified 0

  # Create new solvers and remove those solver which
  # do not have any valid VARIABLE defined (ex. when an 
  # adv-diff variable was deleted)
  set inforce_ids ""

  # Loop all solver infos
  # =====================
  foreach info $system_infos {

    set info_name   [lindex $info 0]
    set info_var    [lindex $info 1]
    set info_act    [lindex $info 2]

    if { $info_act == "" } {
      set info_act 0
    }

    set ename [DataField::fieldNameSifToGui $info_name]
    set eindex [DataField::getFieldProperty Equation $ename EquationIndex]
    set emask [Equation::getEquationMask $eindex]

    set found 0

    # Check all existing solvers
    # ==========================
    foreach id $Solver(ids) {

      set name [DataField::getFieldValue Solver $id EQUATION]
      set var [DataField::getFieldValue Solver $id VARIABLE]
      set act [DataField::getFieldValue Solver $id ACTIVE]

      # If this is one of the solvers listed in "system_info"
      # ====================================================
      if { $info_name == $name && $info_var == $var } {

        # Check parameter
        set pmodified 0

        set Solver($id,data) [DataField::checkParameterData Solver $id pmodified  $emask $eindex]

        if { $pmodified } {
          set modified 1
        }

        #--This will be active solver
        if { $info_act } {

          DataField::setFieldValue Solver $id ACTIVE True
          
          #-If this was inactive, mark modified
          if { $act == "False" } {
            set modified 1
          }

        #--This will be inactive
        } else {

          DataField::setFieldValue Solver $id ACTIVE False

          #-If this was active, mark modified
          if { $act == "True" } {
            set modified 1
          }
        }
    
        #--Add to the actives' list
        set found 1
        lappend inforce_ids $id
        break
      }
    }

    # Existing solver was found (for the info)
    # -------------------------
    if { $found } {
      continue
    }

    # Create a new (default) solver!
    # ------------------------------
      
    set modified 1

    set new_id $Solver(nextNewParameterId)

    incr Solver(nextNewParameterId)

    lappend inforce_ids $new_id

    set data [Panel::createDefaultParameterLine Solver $eindex]

    set Solver($new_id,data) $data
    set Solver($new_id,name) "Solver$new_id"
    set Solver($new_id,bd1) $Info(NO_INDEX)
    set Solver($new_id,bd2) $Info(NO_INDEX)

    DataField::setFieldValue Solver $new_id VARIABLE $info_var
    
    # NOTE, we add also passive solver, if ex, Oxygen for AddDiff
    # was set initially inactive, its solver is needed as an inactive
    # instance, otherwise list sizes will not be correct in SolverSystem!
    #
    if { $info_act } {
      DataField::setFieldValue Solver $new_id ACTIVE True
    } else {
      DataField::setFieldValue Solver $new_id ACTIVE False
    }
    
    # Set correct field activity etc. for the new solver
    #
    set pmodified 0
    set Solver($new_id,data) [DataField::checkParameterData Solver $new_id pmodified  $emask $eindex]

  } ;#end all infos

  # Remove non-active solvers
  # =========================
  foreach id $Solver(ids) {

    if { -1 == [lsearch $inforce_ids $id] } {

      Util::unsetArrayIdVariables Solver $id
      set modified 1
    }
  }

  set Solver(ids) $inforce_ids

  Panel::compressParameterIds Solver
}


# Form list of equation indices sorted by solving order
#
proc Solver::updateActiveSolvingOrderList {} {
  global Equation Solver SolverSystem

  # This is a list (system-index, previous-solving-order-nbr)-pairs
  #
  set sort_list ""

  foreach id $Solver(activeIds) {
    set order [DataField::getFieldValue Solver $id SOLVING_ORDER]
    lappend sort_list [list $id $order]
  }

  # This is the same, but sorted by solving-order number.
  # Now the row-index gives the new order number and first
  # number in the (sorted) pairs is the solver id
  #
  set sorted_list [lsort -index 1 -real $sort_list]

  set SolverSystem(activeSolvingOrderIndices) ""

  foreach pair $sorted_list {

    set id [lindex $pair 0]

    set eq_var [DataField::fieldNameSifToGui [DataField::getFieldValue Solver $id EQUATION]]

    set idx [DataField::getFieldProperty Equation $eq_var EquationIndex]

    lappend SolverSystem(activeSolvingOrderIndices) $idx
  }
}


# Update varaible dofs (after coordinate changes etc)
#
proc Solver::updateVariableDofs { {is_modified ""} } {
  global Solver

  if { $is_modified != "" } {
    upvar $is_modified modified
  }

  set modified 0

  foreach id $Solver(ids) {

    set sname [DataField::getFieldValue Solver $id EQUATION]
    set ename  [DataField::fieldNameSifToGui $sname]
    set eindex [DataField::getFieldProperty Equation $ename EquationIndex]

    set old_dofs [DataField::getFieldValue Solver $id VARIABLE_DOFS]
    set new_dofs [Solver::getEquationVariableDofs $eindex]

    if { $old_dofs != $new_dofs } {
      DataField::setFieldValue Solver $id VARIABLE_DOFS $new_dofs
      set modified 1
    }
  }
}


proc Solver::getEquationIndexByProcedureName { proc_name } {
  global Info Solver

  foreach id $Solver(ids) {

    set proc_info [DataField::getFieldValue Solver $id PROCEDURE]

    set proc_info [lindex [split $proc_info $Info(dataListSeparator)] 0]

    if { $proc_info == $proc_name } {
      
      set eq_name [DataField::getFieldValue Solver $id EQUATION]
      
      set eq_name [DataField::fieldNameSifToGui $eq_name]

      set idx [DataField::getFieldProperty Equation $eq_name EquationIndex]

      if { $idx != "" } {
        return $idx
      }
    }
  }

  return $Info(NO_INDEX)
}


# end ecif_tk_solverDefPanel.tcl
# ********************
