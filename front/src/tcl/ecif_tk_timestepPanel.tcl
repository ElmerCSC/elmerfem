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
#Module:    ecif_tk_timestepDefPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the timestepping properties
#
#************************************************************************
 
#------Timestep definitions  proc------


#This procedure displays the model level definition panel
#
proc Timestep::openPanel { } {
  global Info Common Timestep  Solver Model

  set w .timestepWin
  set Timestep(winName) $w
  set Timestep(winTitle) "Timestep settings"
  set wgeom +420+120

  set id [winfo atom $w]
  set Timestep(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow Timestep $id $Timestep(winTitle) $wgeom] } {
    return
  }  

  set Timestep(dataChanged) 0
  set Timestep(dataModified) 0
  set Timestep(dataDirty) 0

  toplevel $w
  focus $w 
  set this $w

  wm title $w $Timestep(winTitle)
  wm geometry $w $wgeom

  Panel::resetFields Timestep

  set id $Timestep(parameterId)
  DataField::formDataFields Timestep $Timestep($id,data)

  set Timestep(parameterIndex) [lsearch $Timestep(ids) $id]
  set Timestep(parameterName) $Timestep($id,name)
  set Timestep(parameterName,prev) $Timestep($id,name)

  set Timestep(parameterName,err) 0
  set Timestep(parameterName,mod) 0

  # Form parameter name list for timestep-set option menu
  set Timestep(parameterNameList) ""
  foreach id $Timestep(ids) {
    lappend Timestep(parameterNameList) $Timestep($id,name)
  }

  # Map old timestepping method values
  Timestep::checkTimesteppingMethodValue

  # Backup field values
  Panel::backupFields Timestep

  # Store parameters values
  set Timestep(ids,old) $Timestep(ids)
  foreach id $Timestep(ids) {
    set Timestep($id,data,old) $Timestep($id,data)
    set Timestep($id,name,old) $Timestep($id,name)
  }

  # Purely transient vars
  set Timestep(transientGroup,vars) { 
    TIMESTEPPING_METHOD
    TIMESTEP_SIZES
    TIMESTEP_INTERVALS
    OUTPUT_INTERVALS
    timestepEntry
  }

  # Purely steady state vars
  set Timestep(steadyStateGroup,vars) { 
    #STEADY_STATE_OUTPUT_INTERVAL
    #STEADY_STATE_MAX_ITERATIONS
  }


  # ===============
  # WIDGET CREATION
  # ===============

  set bg $Info(nonActiveBg)

  #-Data is organized into four principal frames:

  #----Timesteps set selection
  frame $w.f1 

  #----Time integration related settings
  frame $w.f2 
  #--Time dependency type (transient/stready settings)
  frame $w.f2.f1 ;#-bd 2 -relief groove
  #--Time stepping: integration method and stready state settings
  frame $w.f2.f2

  #----Timesteps entry and listbox
  frame $w.f3 
  #--Timesteps: Intervals and times text window frame
  frame $w.f3.f1 
  #--Timesteps: Total info fields
  frame $w.f3.f2

  #----Buttons
  frame $w.fB 


  # Timestep set selection (and add/delete etc. buttons)
  # ======================
  set m $w.f1

  frame $m.lf 
  frame $m.omf 
  frame $m.ef 
  frame $m.bf 
  
  label $m.lf.l -text "Active set: " -anchor e

  set Timestep(timestepSetMenuFrame) $m.omf

  # Create timestep set selection OptionMenu widget
  # NOTE: This also set the bind-proc!
  #
  set wdg [Timestep::createTimestepSetOptionMenu $m.omf]

  label $m.ef.l -text "Name: "

  set wdg [entry $m.ef.e -textvariable Timestep(parameterName) -width 20 -font $Info(entryFont)]
  bind $wdg <KeyRelease> "+Panel::panelDataModified 1 Timestep $wdg {%A %K}"
  bind $wdg <KeyRelease> "+Widget::entryKeyRelease Timestep parameterName $wdg {%A %K}"
  bind $wdg <KeyPress-Escape> "+Widget::entryKeyPress-Escape Timestep parameterName $wdg {%A %K}"
  set Timestep(parameterNameWidget) $wdg

  # Add timestep button
  # ===================
  set wdg [button $m.bf.add    -text "Add" -width 5 -command Timestep::addTimestepSet]
  bind $wdg <ButtonPress> "+Panel::panelDataModified 0 Timestep $wdg"
  bind $wdg <ButtonPress> "+Panel::panelDataChanged 1 Timestep $wdg"
  set Timestep(panelAddButton) $wdg

  # Update timestep button
  # ======================
  set wdg [button $m.bf.update -text "Update" -width 5 -command Timestep::updateTimestepSet]
  bind $wdg <ButtonPress> "+Panel::panelDataModified 0 Timestep $wdg"
  bind $wdg <ButtonPress> "+Panel::panelDataChanged 1 Timestep $wdg"
  set Timestep(panelUpdateButton) $wdg
  $wdg configure -state disabled

  # Delete timestep button
  # ======================
  set wdg [button $m.bf.delete -text "Delete" -width 5 -command Timestep::deleteTimestepSet]
  if { [llength $Timestep(ids)] == 0 } {
    $wdg configure -state disabled
  }
  bind $wdg <ButtonPress> "+Panel::panelDataModified 0 Timestep $wdg"
  bind $wdg <ButtonPress> "+Panel::panelDataChanged 1 Timestep $wdg"
  set Timestep(panelDeleteButton) $wdg
  

  # Time integration related settings
  # =================================
    
  # Time dependency (transiet/steady)
  # ---------------
  set m [frame $w.f2.f1.tt -bd 2 -relief groove]

  frame $m.lf
  label $m.lf.l  -anchor nw -text "Time dependency: "

  # Add help event for the label using field name TIME_DEPENDENCY
  set wdg $m.lf.l
  set fld_hn TIME_DEPENDENCY
  bind $wdg <ButtonPress-3> "Widget::genericButton-3 Timestep $fld_hn $wdg"

  frame $m.rf
  radiobutton $m.rf.rb1  -text "Steady state"  -value "Steady State"
  radiobutton $m.rf.rb2  -text "Transient"     -value "Transient"

  #-set button variables etc.
  foreach i {1 2 } {
    $m.rf.rb$i configure -anchor w \
                        -variable Timestep(SIMULATION_TYPE) -indicatoron 1\
                        -command "Timestep::setTimestepStatus"  \
                        -selectcolor $Info(select_color)

    set wdg $m.rf.rb$i
    lappend Timestep(allWidgets,SIMULATION_TYPE) $wdg

    bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 Timestep $wdg {%A %K}"
  }

  # Steady state parameters
  # -----------------------
  set m [frame $w.f2.f2.ss -bd 2 -relief groove]

  frame $m.lf
  label $m.lf.l  -anchor nw -text "Steady state settings:"

  set Timestep(steady_coupled,group,label) $m.lf.l
  
  set fld STEADY_STATE_MAX_ITERATIONS
  frame $m.ef1
  label $m.ef1.l -text "Maximum number of iterations:" -width 25 -anchor w
  set Timestep(allWidgets,label,$fld) $m.ef1.l

  entry $m.ef1.e -textvariable Timestep($fld) -width 10 \
                 -font $Info(entryFont)
  set wdg $m.ef1.e

  # NOTE: This field is currently in both groups, MVe 09.03.99
  #
  lappend Timestep(allWidgets,$fld) $wdg
  lappend Timestep(allWidgets,steadyStateGroup) $wdg
  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg

  bind $wdg <KeyRelease> "Panel::panelDataModified 1 Timestep $wdg {%A %K}"
  Widget::setEntryBindings "non_standard" Timestep $fld $wdg

  set fld STEADY_STATE_OUTPUT_INTERVAL
  frame $m.ef2
  label $m.ef2.l -text "Output interval:" -width 25 -anchor w
  set Timestep(allWidgets,label,$fld) $m.ef2.l

  entry $m.ef2.e -textvariable Timestep($fld) -width 10 \
                 -font $Info(entryFont)
  set wdg $m.ef2.e
  set Timestep(steady_coupled,output_interval,entry) $wdg

  lappend Timestep(allWidgets,$fld) $wdg
  lappend Timestep(allWidgets,steadyStateGroup) $wdg

  bind $wdg <KeyRelease> "Panel::panelDataModified 1 Timestep $wdg {%A %K}"
  Widget::setEntryBindings "non_standard" Timestep $fld $wdg


  # Time integration method
  # -----------------------

  set m [frame $w.f2.f2.im -bd 2 -relief groove]

  frame $m.f1 
  frame $m.f2 
  frame $m.f1.lf
  frame $m.f1.omf
  
  set fld TIMESTEPPING_METHOD
  set Timestep(methodParentFrame) $w.f2.f2.im
  set Timestep(methodFrame) $w.f2.f2.im.f2
  set Timestep(allWidgets,$fld) ""

  set FV $fld
  set f_var Timestep(TIMESTEPPING_METHOD)
  set lbl [DataField::getFieldProperty Timestep $FV Label]
  set values [lindex [DataField::getFieldProperty Timestep $FV Limits] 1]
  set ww [DataField::getFieldProperty Timestep $FV WidgetWidth]

  label $m.f1.lf.l -anchor nw -text $lbl
  set Timestep(allWidgets,label,$fld) $m.f1.lf.l

  set wdg $m.f1.omf.om
  set om [Panel::createOptionMenuWidget $wdg  $f_var $values]

  $m.f1.omf.om configure -activebackground $Info(optionMenuBg) \
                         -indicatoron 0  -width $ww -state normal

  bind $om <<MenuSelect>> "PanelCheck::checkOptionMenu $om $wdg Timestep $FV Timestep::setTimesteppingMethod {%s}"

  #bind $om <<MenuSelect>> "Timestep::setTimesteppingMethod"

  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg
  
  #--Set actual method widgets
  Timestep::setTimesteppingMethod
  
  # Timesteps entry and listbox
  # ===========================
  set m $w.f3.f1

  frame $m.lf
  frame $m.lf.lf1
  frame $m.lf.lf2
  label $m.lf.lf1.l  -anchor nw -text "Timestep entry:"

  # Add help event for the label using field name TIMESTEP_ENTRY
  set wdg $m.lf.lf1.l
  set fld_hn TIMESTEP_ENTRY
  bind $wdg <ButtonPress-3> "Widget::genericButton-3 Timestep $fld_hn $wdg"

  label $m.lf.lf2.l1  -text "Step\nsize \[s\]" -width 10 -anchor e
  label $m.lf.lf2.l2  -text "Steps in\ninterval" -width 12 -anchor c
  label $m.lf.lf2.l3  -text "Output\ninterval" -width 8 -anchor c
  label $m.lf.lf2.l4  -text "Cumulatives" -width 40 -anchor c

  # Timestep entry
  # --------------
  frame $m.ef

  # NOTE This entry is with table-font!
  set wdg [entry $m.ef.e -textvariable Timestep(timestepEntry) -width 55 \
                         -font $Info(tableFont)]

  set Timestep(timestepEntryWidget) $wdg

  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg

  bind $wdg <KeyRelease> "Panel::panelDataDirty 1 Timestep $wdg {%A %K}"
  bind $wdg <KeyPress-Return> "Timestep::addTimestep"
  bind $wdg <KeyPress-Return> "+Panel::panelDataModified 1 Timestep"
  
  # Timesteps listbox
  # -----------------
  frame $m.lbf

  listbox $m.lbf.lb -width 55 -height 8 \
                    -relief sunken -selectmode extended \
                    -yscrollcommand [list $m.lbf.scrollbar_y set] \
                    -exportselection 0 -font $Info(tableFont)
  scrollbar $m.lbf.scrollbar_y -orient vertical -command [list $m.lbf.lb yview]
  set Timestep(timestepListLB) $m.lbf.lb

  bind $Timestep(timestepListLB) <Double-1> "+Timestep::copyTimestepListRow"

  # Add, Delete etc. buttons
  # ================
  frame $m.lbf.btnf

  button $m.lbf.btnf.add    -text "Add"     -width 8 -command Timestep::addTimestep  
  button $m.lbf.btnf.insert -text "Insert"  -width 8 -command Timestep::insertTimestep
  button $m.lbf.btnf.update -text "Update"  -width 8 -command Timestep::updateTimestep
  button $m.lbf.btnf.delete -text "Delete"  -width 8 -command Timestep::deleteTimestep
  
  #lappend Timestep(allWidgets,transientGroup) $m.lbf.lb
  #lappend Timestep(allWidgets,transientGroup) $m.lbf.scrollbar_y

  set wdg $m.lbf.btnf.add
  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg

  set wdg $m.lbf.btnf.insert
  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg

  set wdg $m.lbf.btnf.update
  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg

  set wdg $m.lbf.btnf.delete
  lappend Timestep(allWidgets,transientGroup) $wdg
  lappend Timestep(allWidgets,transientGroupFixed) $wdg

 
  # Totals (info fields)
  # ====================
  set m $w.f3.f2
  frame $m.lf
  label $m.lf.l  -anchor nw -text "Totals:"
  frame $m.ef
  entry $m.ef.e1 -textvariable Timestep(end_time) -width 12 \
                 -state disabled -font $Info(entryFont) -justify right -bg $bg
  entry $m.ef.e2 -textvariable  Timestep(nof_timesteps) -width 12 \
                 -state disabled -font $Info(entryFont) -justify right  -bg $bg
  entry $m.ef.e3 -textvariable Timestep(nof_output_steps) -width 12 \
                 -state disabled -font $Info(entryFont) -justify right  -bg $bg

  label $m.ef.l4 -text "End time:"

  entry $m.ef.e4 -textvariable Timestep(end_time2) -width 15 \
                 -state disabled -font $Info(entryFont) -justify right -bg $bg

 
  Timestep::setTimestepStatus

  # ==============
  # WIDGET PACKING
  # ==============

  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #-Packing principal frames
  pack $w.f1 -side top -padx $fpx1 -pady $fpy3
  pack $w.f2 -side top -anchor nw -fill x -expand 1 -padx $fpx1 -pady $fpy3
  pack $w.f3 -side top -anchor nw -fill x -expand 1 -padx $fpx1 -pady $fpy3
  pack $w.fB -side top -padx $fpx3 -pady $fpy3

  #-Frame internals

  #---Timestep sets
  set m $w.f1

  pack $m.lf $m.omf $m.ef $m.bf -side left  -anchor w -padx $fpx1 -pady $fpy3
  pack $m.lf.l $m.omf.m $m.ef.l $m.ef.e -side left  -anchor w -pady $fpy3
  pack $m.bf.add $m.bf.update $m.bf.delete -side left  -anchor w -padx $fpx1 -pady $fpy3

  #---Time integration
  pack $w.f2.f1 $w.f2.f2 -side left -padx $fpx1 -pady $fpy1

  #-Time type
  set m $w.f2.f1.tt

  pack $m -side top  -anchor center -fill both -padx $fpx3 -pady $fpy1
  pack $m.lf -side top  
  pack $m.lf.l -side top -anchor nw  
  pack $m.rf -side top  
  pack $m.rf.rb1 -side top -anchor nw  
  pack $m.rf.rb2 -side top -anchor nw  

  #-Steady state iterations
  set m $w.f2.f2.ss

  pack $m.lf -side top  
  pack $m.lf.l -side left -anchor nw  

  pack $m.ef1 $m.ef2 -side top   

  pack $m.ef1.l -side left -padx $fpx1 -pady $fpy1  
  pack $m.ef1.e -side left -padx $fpx1 -pady $fpy1  

  pack $m.ef2.l -side left -padx $fpx1 -pady $fpy1  
  pack $m.ef2.e -side left -padx $fpx1 -pady $fpy1  

  pack $m -side top  -anchor nw -fill x -expand 1 -padx $fpx2 -pady $fpy2

  #-Time integration method
  set m $w.f2.f2.im

  pack $m.f1 $m.f2 -side top
  pack $m.f1.lf $m.f1.omf -side left
  pack $m.f1.lf.l  $m.f1.omf.om -side left -anchor nw  
  pack $m -side top  -anchor nw -fill x -expand 1 -padx $fpx2 -pady $fpy2


  #---Timesteps
  #-Times and intervals
  set m $w.f3.f1
  pack $m.lf $m.ef $m.lbf -side top  -anchor w

  pack $m.lf.lf1 $m.lf.lf2  -side top  -anchor w
  pack $m.lf.lf1.l -side left -anchor w  
  pack $m.lf.lf2.l1 $m.lf.lf2.l2 $m.lf.lf2.l3 $m.lf.lf2.l4 \
        -side left -anchor w    

  pack $m.ef.e       
  set Timestep(timestepEntryWidget) $m.ef.e

  pack  $m.lbf.lb $m.lbf.scrollbar_y $m.lbf.btnf -side left -fill y     
  pack  $m.lbf.btnf.add \
        $m.lbf.btnf.insert \
        $m.lbf.btnf.update \
        $m.lbf.btnf.delete \
        -side top -anchor w -padx $fpx1 -pady $fpy1
  pack $m -side top  -anchor nw -fill x -expand 1 -padx $fpx1 -pady $fpy1

  #-Timesteps totals (info)
  set m $w.f3.f2

  pack $m.lf -side top  
  pack $m.lf.l -side left -anchor nw  
  pack $m.ef -side top   
  pack $m.ef.e1 $m.ef.e2 $m.ef.e3 -side left -padx $fpx2
  pack $m.ef.l4 -side left -anchor e -padx $fpx1
  pack $m.ef.e4 -side left -anchor e 
  pack $m -side top  -anchor nw -fill x -expand 1 -padx $fpx2 -pady $fpy3

  #---Apply/Ok buttons (create and pack)

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "Timestep::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "Timestep::panelCancel $this" -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command Timestep::panelApply -state $ap]

  focus $ok_btn
  set Timestep(applyButton)  $ap_btn
  set Timestep(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx1 
  
  #-Init timestep list and construct list-boxlist
  Timestep::initTimestepList
  Timestep::updateTimestepListData


  # Pick current packing size!
  set wdg $Timestep(methodParentFrame)
  Util::updateGui
  set wid [winfo reqwidth $wdg]
  set hig [winfo reqheight $wdg]

  $wdg configure -height $hig
  $wdg configure -width  $wid
  pack propagate $wdg 0

  # At least one timestep set is needed
  if { [llength $Timestep(ids)] == 1 } {
    $Timestep(panelDeleteButton) configure -state disabled
  }

  # Set field label bindings for right-button help
  Widget::setLabelBindings Timestep
}


# Map old timestepping method values
#
proc Timestep::checkTimesteppingMethodValue {} {
  global Timestep

  set id $Timestep(parameterId)

  set method [DataField::getFieldValue Timestep $id TIMESTEPPING_METHOD]

  if { $method == "" } {
    return
  }

  set Timestep(dataModified) 1

  switch $method {

    "Implicit Euler" { 
      DataField::setFieldValue Timestep $id TIMESTEPPING_METHOD "Newmark"
      DataField::setFieldValue Timestep $id NEWMARK_BETA 1.0
    }

    "Explicit Euler" { 
      DataField::setFieldValue Timestep $id TIMESTEPPING_METHOD "Newmark"
      DataField::setFieldValue Timestep $id NEWMARK_BETA 0.0
    }

    "Crank Nicolson" { 
      DataField::setFieldValue Timestep $id TIMESTEPPING_METHOD "Newmark"
      DataField::setFieldValue Timestep $id NEWMARK_BETA 0.5
    }

    default {
      set Timestep(dataModified) 0
    }
  }
}


proc Timestep::setTimesteppingMethod {} {
  global Timestep

  foreach wdg $Timestep(allWidgets,TIMESTEPPING_METHOD) {
    destroy $wdg
  }

  set Timestep(allWidgets,TIMESTEPPING_METHOD) ""
  set Timestep(allWidgets,transientGroup) $Timestep(allWidgets,transientGroupFixed)

  switch $Timestep(TIMESTEPPING_METHOD) {
    "Newmark" { Timestep::setTimesteppingMethod_Newmark }
    "BDF"     { Timestep::setTimesteppingMethod_BDF }
  }

  set Timestep(TIMESTEPPING_METHOD,prev) $Timestep(TIMESTEPPING_METHOD)
}


proc Timestep::setTimesteppingMethod_Newmark {} {
  global Info Timestep

  Timestep(setNewmarkBeta)

  set Timestep(NEWMARK_BETA,act) 1
  set Timestep(BDF_ORDER,act) 0

  set m $Timestep(methodFrame)

  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [frame $m.f1]
  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [frame $m.f2]

  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [frame $m.f2.lf]
  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [frame $m.f2.ef]

  set sc $Info(select_color)

  #---Selection radiobuttons
  foreach i { 1 2 3 } \
          t { "Implicit Euler" "Explicit Euler" "Crank Nicolson" } \
          v { 1.0 0.0 0.5 } {
    set wdg [radiobutton $m.f1.rb$i -text $t -value $v \
                                    -variable Timestep(NewmarkBeta) \
                                    -selectcolor $sc]
    pack $wdg -side left

    lappend Timestep(allWidgets,TIMESTEPPING_METHOD) $wdg
    lappend Timestep(allWidgets,transientGroup) $wdg

    bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 Timestep $wdg {%A %K}"
    bind $wdg <ButtonRelease-1> "+Timestep(getNewmarkBeta) $wdg"
  }

  #---Newmark beta entry
  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [label $m.f2.lf.l -text "Newmark beta:"]
  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [entry $m.f2.ef.e -textvariable Timestep(NEWMARK_BETA)]

  set wdg $m.f2.ef.e
  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) $wdg
  lappend Timestep(allWidgets,transientGroup) $wdg

  set Timestep(allWidgets,NEWMARK_BETA) $wdg
  Timestep(setNewmarkBetaStatus)

  bind $wdg <KeyRelease> "Panel::panelDataModified 1 Timestep $wdg {%A %K}"
  bind $wdg <KeyRelease> "+Timestep(setNewmarkBeta)"
  Widget::setEntryBindings "non_standard" Timestep NEWMARK_BETA $wdg

  #---Pack widgets
  pack $m.f1 $m.f2 -side top
  pack $m.f2.lf $m.f2.ef -side left
  pack $m.f2.lf.l $m.f2.ef.e -side left
}


# This proc is called when NewmarkBeta entry widget
# is modified
#
proc Timestep(setNewmarkBeta) {} {
  global Timestep

  # Set reference variable value for the radiobuttons
  if { $Timestep(NEWMARK_BETA) == 0 } {
    set Timestep(NewmarkBeta) "0.0"

  } elseif { $Timestep(NEWMARK_BETA) == 0.5 } {
    set Timestep(NewmarkBeta) "0.5"

  } elseif { $Timestep(NEWMARK_BETA) == 1.0 } {
    set Timestep(NewmarkBeta) "1.0"

  } else {
    set Timestep(NewmarkBeta) $Timestep(NEWMARK_BETA)
  }
}


# This set status for the NewmarkBeta entry field
#
proc Timestep(setNewmarkBetaStatus) {} {
  global Timestep

  # Mark entry modified status
  if { $Timestep(NEWMARK_BETA) != $Timestep(NEWMARK_BETA,old) } {
    set Timestep(NEWMARK_BETA,mod) 1
  } else {
    set Timestep(NEWMARK_BETA,mod) 0
  }

  Widget::setEntryStatus Timestep NEWMARK_BETA $Timestep(allWidgets,NEWMARK_BETA)
}



# This proc is called when one of the NewmarkBeta
# radiobuttons is presseed
#
proc Timestep(getNewmarkBeta) {wdg} {
  global Timestep

  # Read radiobutton value
  set Timestep(NEWMARK_BETA) [$wdg cget -value]

  Timestep(setNewmarkBetaStatus)
}


proc Timestep::setTimesteppingMethod_BDF {} {
  global Info Timestep

  set Timestep(BDF_ORDER,act) 1
  set Timestep(NEWMARK_BETA,act) 0

  set m $Timestep(methodFrame)

  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) [frame $m.f1]

  set sc $Info(select_color)

  # Label
  set wdg [label $m.f1.l -text "BDF Order: " -width 10]
  lappend Timestep(allWidgets,TIMESTEPPING_METHOD) $wdg
  pack $wdg -side left

  #---Selection radiobuttons
  foreach i { 1 2 3 4 5} {
    set wdg [radiobutton $m.f1.rb$i -text $i -value $i \
                                    -variable Timestep(BDF_ORDER) \
                                    -selectcolor $sc]
    pack $wdg -side left

    lappend Timestep(allWidgets,TIMESTEPPING_METHOD) $wdg
    lappend Timestep(allWidgets,transientGroup) $wdg

    bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 Timestep $wdg {%A %K}"
  }

  #---Pack widgets
  pack $m.f1 -side top
}


# Utility
#
proc Timestep::isTransient {} {
  global Timestep

  if { $Timestep(SIMULATION_TYPE) == "Transient" } {
    return 1
  } else {
    return 0
  }
}


proc Timestep::setTimestepStatus {} {
  global Timestep Solver Info

  set bg $Info(nonActiveBg)

  foreach w $Timestep(allWidgets,steadyStateGroup) {
    if { "Entry" == [winfo class $w] } {
      Widget::configureEntry $w disabled
    } else {
      $w configure -state disabled -fg gray
    }
  }

  foreach w $Timestep(allWidgets,transientGroup) {
    if { "Entry" == [winfo class $w] } {
      Widget::configureEntry $w disabled
    } else {
      $w configure -state disabled -fg gray
    }
  }

  #--Steady/transient settings

  #-Steady state is on
  if { $Timestep(SIMULATION_TYPE) == "Steady State"} {

    foreach w $Timestep(allWidgets,steadyStateGroup) {
      if { "Entry" == [winfo class $w] } {
        Widget::configureEntry $w normal
      } else {
        $w configure -state normal -fg black
      }
    }

  #-Transient is on
  } else {

    foreach w $Timestep(allWidgets,transientGroup) {
      if { "Entry" == [winfo class $w] } {
        Widget::configureEntry $w normal
      } else {
        $w configure -state normal -fg black
      }
    }

    if { [Widget::arrayWidgetExists Timestep NEWMARK_BETA] } {
      Timestep(setNewmarkBetaStatus)
    }

    Timestep::updateTimestepListAndTotals
  }

  set wdg $Timestep(steady_coupled,group,label)

  if { [Timestep::isSteadyState] } {
    $wdg configure -text "Steady state settings:"
  } else {
    $wdg configure -text "Coupled equation settings:"
  }

} 


# Get simulation type
#
proc Timestep::isSteadyState {} {
  global Timestep

  if { [info exists Timestep(SIMULATION_TYPE)] &&
       $Timestep(SIMULATION_TYPE) == "Transient"
     } {
    return 0

  } else {
    return 1
  }
}


# Set SteadyState/Coupled Equation group labels etc.
# widget: target widget
# option: option to configure (like: -text)
# values: list of possible values, index:  0 <--> transient 1 <--> steady state
#
proc Timestep::steadyCoupledConfigure {widget option values} {
  global Timestep

  if { $Timestep(SIMULATION_TYPE) == "Transient" } {
    $widget configure $option [lindex $values 0]
  } else {
    $widget configure $option [lindex $values 1]
  }
}


# ==================
# Timestep SET procs
# ==================

# Add button proc for timestep set
#
proc Timestep::addTimestepSet {} {
  global Info Timestep

  if { ![Timestep::checkTransientData] } {
    return 0
  }

  # New id and name
  set pid $Timestep(nextNewParameterId)
  set name $Timestep(parameterName)

  # Check new name
  if { ![Panel::checkParameterName Timestep $pid $name 1] } {
    return 0
  }

  # Note: Name may have been set to a default!
  set name $Timestep(parameterName)

  # Update id and index stuff
  set Timestep(parameterId) $pid
  lappend Timestep(ids) $pid
  set Timestep(parameterIndex) [expr [llength $Timestep(ids)] - 1 ]

  # Create new parameter
  set data_size [llength $Timestep(timestepList)]
  set Timestep(TIMESTEP_SIZES,dataSize) $data_size
  set Timestep(TIMESTEP_INTERVALS,dataSize) $data_size
  set Timestep(OUTPUT_INTERVALS,dataSize) $data_size

  DataField::formNonStandardParameter Timestep $pid $name

  incr Timestep(nextNewParameterId)

  # Update selection menu
  lappend Timestep(parameterNameList) $name
  set wdg [Timestep::createTimestepSetOptionMenu $Timestep(timestepSetMenuFrame)]
  pack $wdg

  # Make current "selected"
  Panel::panelDataModified 0 Timestep
  Timestep::setTimestepSet

  $Timestep(panelDeleteButton) configure -state normal

  return 1
}   


# Delete button proc for timestep set
proc Timestep::deleteTimestepSet {} {
  global Info Timestep

  set pid $Timestep(parameterId)
  set index $Timestep(parameterIndex)

  # Delete parameter data
  # ================================
  # NOTE: old-values are not deleted with these, because
  # variable names are given exactly, not like data* !!!
  #
  Util::unsetArrayIdVariables Timestep $pid data
  Util::unsetArrayIdVariables Timestep $pid name
  Util::unsetArrayIdVariables Timestep $pid bd1
  Util::unsetArrayIdVariables Timestep $pid bd2

  set Timestep(ids) [lreplace $Timestep(ids) $index $index]
  set Timestep(parameterNameList) [lreplace $Timestep(parameterNameList) $index $index]

  if { $Timestep(parameterIndex) > 0 } {
    incr Timestep(parameterIndex) -1
  }

  # Update selection menu
  set name [lindex $Timestep(parameterNameList) $Timestep(parameterIndex)]
  set Timestep(parameterName) $name
  set wdg [Timestep::createTimestepSetOptionMenu $Timestep(timestepSetMenuFrame)]
  pack $wdg

  # Make current "selected"
  Panel::panelDataModified 0 Timestep
  Timestep::setTimestepSet

  if { [llength $Timestep(ids)] == 1 } {
    $Timestep(panelDeleteButton) configure -state disabled
  }

  return 1

}   


# Update button proc for timestep set
#
proc Timestep::updateTimestepSet {} {
  global Info Timestep

  if { ![Timestep::checkTransientData] } {
    return 0
  }

  set pid $Timestep(parameterId)
  set index $Timestep(parameterIndex)

  # Check possibly updated name
  set name $Timestep(parameterName)

  if { ![Panel::checkParameterName Timestep $pid $name 0] } {
    return 0
  }

  # Note: Name may have changed to default!
  set name $Timestep(parameterName)

  # Update parameter
  set data_size [llength $Timestep(timestepList)]
  set Timestep(TIMESTEP_SIZES,dataSize) $data_size
  set Timestep(TIMESTEP_INTERVALS,dataSize) $data_size
  set Timestep(OUTPUT_INTERVALS,dataSize) $data_size

  DataField::formNonStandardParameter Timestep $pid $name

  # Update selection menu
  set tmp $Timestep(parameterNameList)
  set Timestep(parameterNameList) [lreplace $tmp $index $index $name]
  set wdg [Timestep::createTimestepSetOptionMenu $Timestep(timestepSetMenuFrame)]
  pack $wdg

  # Make current "selected"
  Panel::panelDataModified 0 Timestep
  Timestep::setTimestepSet

  return 1
}   

# Timestep set selected from the option menu
#
proc Timestep::selectTimestepSet {} {
  global Timestep

  # If same reselected do nothing!!
  #
  if { $Timestep(parameterName) == $Timestep(parameterName,prev) } {
    return
  }

  # If current was changed but not updated, verify from user
  #
  if { $Timestep(dataModified) && ![Panel::verifyParameter Timestep] } {
    set Timestep(parameterName) $Timestep(parameterName,prev)
    return
  }

  Timestep::setTimestepSet
}


# Write selected timestep set into panel
#
proc Timestep::setTimestepSet {} {
  global Timestep

  if { $Timestep(ids) == "" } {
    return
  }

  set index [lsearch $Timestep(parameterNameList) $Timestep(parameterName)]

  set Timestep(parameterName,prev) $Timestep(parameterName)

  Panel::panelDataChanged 1 Timestep
  Panel::panelDataModified 0 Timestep

  set Timestep(parameterId) [lindex $Timestep(ids) $index]
  set Timestep(parameterIndex) $index

  foreach id $Timestep(ids) {
    DataField::setFieldValue Timestep $id ACTIVE False
  }

  # Make the selective the active one!
  set pid $Timestep(parameterId)

  DataField::setFieldValue Timestep $pid ACTIVE True

  #set Timestep(parameterName) $Timestep($pid,name)

  set Timestep(parameterName,err) 0
  set Timestep(parameterName,mod) 0
  Widget::setEntryStatus Timestep parameterName $Timestep(parameterNameWidget)

  Panel::resetFields Timestep

  DataField::formDataFields Timestep $Timestep($pid,data)

  # Update timestep listbox and data
  Timestep::setTimesteppingMethod
  Timestep::initTimestepList
  Timestep::updateTimestepListData

  Timestep::setTimestepStatus
}   


# Create new timestep set option menu, when timestep sets have been
# Added/Deleted/Renamed (updated)
#
proc Timestep::createTimestepSetOptionMenu {om_frame} {
  global Info Timestep

  set fld_var parameterName
  set menu_var Timestep(parameterName)

  set menu_values $Timestep(parameterNameList)

  set wdg $om_frame.m

  destroy $wdg

  set om [Panel::createOptionMenuWidget $wdg $menu_var $menu_values]

  $wdg configure -indicatoron 0 -width 10 -anchor w \
                 -activebackground $Info(optionMenuBg)

  bind $om <<MenuSelect>> "PanelCheck::checkOptionMenu $om $wdg Timestep $fld_var Timestep::selectTimestepSet {%s}"

  return $wdg
}


# ====================
# Timestep ENTRY procs
# ====================

proc Timestep::initTimestepList {} {
  global Timestep

  set Timestep(timestepList) ""
  set Timestep(timestepEntry) ""

  # Check that data exists and is consistent
  if { $Timestep(TIMESTEP_SIZES) == "" } {
    set Timestep(TIMESTEP_INTERVALS) ""
    set Timestep(OUTPUT_INTERVALS) ""
    return
  }

  foreach size  $Timestep(TIMESTEP_SIZES)      \
          steps $Timestep(TIMESTEP_INTERVALS)  \
          outs  $Timestep(OUTPUT_INTERVALS) {

    # data and placeholders for cumulatives
    set data_row "$size $steps $outs 0.0 0 0"
    # add to the listbox list
    lappend Timestep(timestepList) $data_row
  }
}


# Udpates timestepList, totals and timestepLB 
#
proc Timestep::updateTimestepListData {} {
  global Timestep

  Timestep::updateTimestepListAndTotals

  ListBox::fill $Timestep(timestepListLB) $Timestep(timestepList)

  $Timestep(timestepListLB) see end
}


proc Timestep::updateTimestepListAndTotals {} {
  global Timestep

  set Timestep(end_time) ""
  set Timestep(end_time2) ""
  set Timestep(nof_timesteps)  ""
  set Timestep(nof_output_steps) ""

  if { ![Timestep::isTransient]  ||
       ![info exists Timestep(timestepList)]
     } {
    return
  }

  set data_list ""
  set cumul_steps 0
  set cumul_time 0.0
  set cumul_outs 0

  foreach data_row $Timestep(timestepList) {
    set size [lindex $data_row 0]
    set steps [lindex $data_row 1]
    set outs [lindex $data_row 2]

    #
    set cumul_time [expr $cumul_time + ($steps * $size)]
    incr cumul_steps $steps

    if { $outs != 0 } {
      incr cumul_outs [expr $steps / $outs]
    }

    set data_row ""

    append data_row [Util::formatNumber $size 11 1]
    append data_row [format "%6d"   $steps]
    append data_row [format "%6d"   $outs]

    append data_row [Util::formatNumber $cumul_time 15 1]
    append data_row [format "%9d"  $cumul_steps]
    append data_row [format "%8d"  $cumul_outs]

    lappend data_list $data_row
  }

  set Timestep(timestepList) $data_list

  set Timestep(end_time) $cumul_time

  set msg ""
  catch { [set Timestep(end_time2) [Util::timeSec2HourStrExact $Timestep(end_time)] ] msg}

  set Timestep(nof_timesteps)  $cumul_steps
  set Timestep(nof_output_steps) $cumul_outs
}
 

# Add button proc for a timestep entry, append
# new row at the end timestepList
#
proc Timestep::addTimestep {} {
  global Timestep

  Panel::panelDataModified 1 Timestep
  Panel::panelDataDirty 0 Timestep

  if { ![Timestep::checkTimestepEntry] } {
    return
  }

  $Timestep(timestepEntryWidget) icursor 0
  $Timestep(timestepEntryWidget) selection range 0 end

  set new_row $Timestep(timestepEntry)
  append $new_row "0.0 0 0"
  lappend Timestep(timestepList) $new_row

  Timestep::updateTimestepListData
}


# Delete button proc for a timestep entry, delete 
# current selection from the timestepList
#
proc Timestep::deleteTimestep {} {
  global Timestep Info

  Panel::panelDataModified 1 Timestep

  set indices [$Timestep(timestepListLB) curselection]

  if {$indices == ""} {

    set Info(messageIcon) error
    Message::dispOkMessage {"Select first the row to be deleted!"} \
                   "$Info(FRONT_NAME) message" \
                   $Timestep(winName)
    return
  }

  set reversed_indices ""

  foreach index $indices {
    set reversed_indices [linsert $reversed_indices 0 $index]
  }

  foreach index $reversed_indices {
    set Timestep(timestepList) [lreplace $Timestep(timestepList) \
                                      $index $index ]
  }

  Timestep::updateTimestepListData
}


# Insert button proc for a timrstep entry, insert a new row into
# timestepList above the current selection
#
proc Timestep::insertTimestep {} {
  global Timestep Info

  Panel::panelDataModified 1 Timestep
  Panel::panelDataDirty 0 Timestep

  set index [$Timestep(timestepListLB) curselection]

  if {$index == ""} {

    set Info(messageIcon) error
    Message::dispOkMessage {"Select first the row above which data is inserted!"} \
                   "$Info(FRONT_NAME) message" \
                   $Timestep(winName)
    return
  }

  if { [llength $index ] > 1} {

    set Info(messageIcon) error
    Message::dispOkMessage {"Select only one row!"} \
                   "$Info(FRONT_NAME) message" \
                   $Timestep(winName)
    return
  }

  if { ![Timestep::checkTimestepEntry] } {
    return
  }

  set new_row $Timestep(timestepEntry)
  append $new_row "0.0 0 0"

  set Timestep(timestepList) [linsert $Timestep(timestepList) \
                                      $index $new_row]

  Timestep::updateTimestepListData
}


# Updates button proc for a timestep entry, update the current
# selection in the timestepList
#
proc Timestep::updateTimestep {} {
  global Timestep Info

  Panel::panelDataModified 1 Timestep
  Panel::panelDataDirty 0 Timestep

  set index [$Timestep(timestepListLB) curselection]

  if {$index == ""} {

    set Info(messageIcon) error
    Message::dispOkMessage {"Select first the row to be updated!"} \
                  "$Info(FRONT_NAME) message" \
                  $Timestep(winName)
    return
  }

  if { [llength $index ] > 1} {

    set Info(messageIcon) error
    Message::dispOkMessage {"Select only one row!"} \
                   "$Info(FRONT_NAME) message" \
                   $Timestep(winName)
    return
  }

  if { ![Timestep::checkTimestepEntry] } {
    return
  }
  
  set new_row $Timestep(timestepEntry)
  append $new_row "0.0 0 0"
  set Timestep(timestepList) [lreplace $Timestep(timestepList) \
                                      $index $index $new_row]

  Timestep::updateTimestepListData
}


# Updates system according to the Steady/Transient state
#
proc Timestep::applyTimestepType {} {
   global Timestep

  #-Stationary problem, Initial condition menu --> inactive
  if {$Timestep(SIMULATION_TYPE) == "Steady State"} {
    #MenuBuild::configurePanelState Model initial disabled
 
  #-Transient problem, Initial condition menu --> active
  } elseif {$Timestep(SIMULATION_TYPE) == "Transient"} {
    MenuBuild::configurePanelState Model initial normal
  }
}


# Currently only time-type is set on model level
# Here we have to only set initial-condition menu-state
#
proc Timestep::updateProblemMask {changedType} {
  global Timestep

  switch $changedType {
    timestepType  Timestep::applyTimestepType
  }
} 


proc Timestep::panelSave { {inform_front 1} } {
  global Info Timestep Model

  if { [llength $Timestep(ids)] == 0 } {
    Panel::initFields Timestep "" 1
    set Timestep(nextNewParameterId) 1
    Timestep::addTimestepSet
  }

  set panel $Timestep(parameterType)

  #--Store old field values
  Panel::backupFields Timestep

  Panel::compressParameterIds Timestep

  #---Store parameters values
  set Timestep(ids,old) $Timestep(ids)
  foreach id $Timestep(ids) {
    set Timestep($id,data,old) $Timestep($id,data)
    set Timestep($id,name,old) $Timestep($id,name)
  }

  #--Set timstepping data size
  if { $Timestep(SIMULATION_TYPE) == "Transient" } {
    
    # NOTE: This must be done, if panel save proc
    # is called outside!
    if { ![info exist Timestep(timestepList)] } {
      Timestep::initTimestepList
    }
 
    set data_size [llength $Timestep(timestepList)]
    set Timestep(TIMESTEP_SIZES,dataSize) $data_size
    set Timestep(TIMESTEP_INTERVALS,dataSize) $data_size
    set Timestep(OUTPUT_INTERVALS,dataSize) $data_size
  }

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  set Model(Solver,needsUpdate) 1

  Panel::panelDataModified 0 Timestep
  Panel::panelDataChanged 0 Timestep 
  Panel::panelDataDirty 0 Timestep
  
  Util::cpp_exec timestepPanelOk

  #---Update Solver parameter panel
  Solver::setSteadyCoupledLabel
}


### OK proc ###

proc Timestep::panelOk {w} {
  global Timestep Info

  #---Check that all changes are applied
  if { ![Panel::verifySave Timestep] } {
    return
  }

  set Timestep(dataModified) 0
  set Timestep(dataDirty) 0

  #---No changes
  if { !$Timestep(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error
  if { ![Timestep::checkPanelData] } {
    return
  }

  #---Check that data makes sense

  #-Transient case
  if { $Timestep(SIMULATION_TYPE) == "Transient" } {

    if { $Timestep(TIMESTEP_SIZES) == "" } {
      set msg [list "No timestep entries were defined for a transient problem!\n" \
                    "Model cannot be solved!\n\n" \
                    $Info(anywayOk) ]

      set Info(messageIcon) warning
      set result [Message::dispCancelOkMessage \
                      $msg \
                      "$Info(FRONT_NAME) error message" \
                      $Timestep(winName) ]

      if { $result == "cancel" } {
        return
      }
    }
  }
  
  #-Steady case  
  if { $Timestep(SIMULATION_TYPE) == "Steady State" } {
    if { $Timestep(STEADY_STATE_OUTPUT_INTERVAL) == "" } {
      set msg [list "NOTE: No output interval was defined for a steady state problem!\n\n" \
                    $Info(anywayOk) ]

      set Info(messageIcon) warning
      set result [Message::dispCancelOkMessage \
                      $msg \
                      "$Info(FRONT_NAME) error message" \
                      $Timestep(winName) ]
      if { $result == "cancel" } {
        return
      }
    }
  }

  #---Save data
  Timestep::panelSave
  Panel::cancel $w
}


proc Timestep::panelApply {} {
  global Info Timestep

  #---Check that all changes are applied
  if { ![Panel::verifySave Timestep] } {
    return
  }

  #---No changes
  if { !$Timestep(dataChanged) } {
    return
  }

  #---Check that all changes are applied
  if { ![Panel::verifySave Timestep] } {
    return
  }

  #---Error in data
  if { ![Timestep::checkPanelData]} {
    return
  }

  Timestep::updateProblemMask timestepType

  Timestep::panelSave
}


proc Timestep::panelCancel {w} {
  global Info Timestep 

  if { ![Panel::verifyCancel Timestep] } {
    return
  }

  #---Reset curent fields into old values
  Panel::restoreFields Timestep

  #---Reset parameters into old values
  set Timestep(ids) $Timestep(ids,old)
  foreach id $Timestep(ids) {
    set Timestep($id,data) $Timestep($id,data,old)
    set Timestep($id,name) $Timestep($id,name,old)
  }

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc Timestep::checkPanelData {} {

  #-Check time stepping data
  set code [Timestep::checkTransientData]
  if {$code == 0} {
    return 0
  }

  #-Check steady state data
  set code [Timestep::checkSteadyStateData]
  if {$code == 0} {
    return 0
  }
 
  # Data ok
  return 1
}


# Check Transient data
# Return 1 = ok, 0 = error
#
proc Timestep::checkTransientData {} {
  global Timestep

  if { $Timestep(SIMULATION_TYPE) == "Steady State" } {
    return 1
  }

  #--Update variables from timestepList
  set Timestep(TIMESTEP_SIZES) ""
  set Timestep(TIMESTEP_INTERVALS) ""
  set Timestep(OUTPUT_INTERVALS) ""

  foreach data_row $Timestep(timestepList) {

    lappend Timestep(TIMESTEP_SIZES) [lindex $data_row 0]
    lappend Timestep(TIMESTEP_INTERVALS) [lindex $data_row 1]
    lappend Timestep(OUTPUT_INTERVALS) [lindex $data_row 2]

  }

  return 1
}


# Check Steady State data
# Return 1 = ok, 0 = error
#
proc Timestep::checkSteadyStateData {} {
  global Timestep Info

  if { $Timestep(SIMULATION_TYPE) == "Transient" } {
    return 1
  }
  
  #-Output steps in the interval
  #
  set outs [ string trim $Timestep(STEADY_STATE_OUTPUT_INTERVAL) ]

  set msg ""
  if { $outs != "" } {
    set outs [DataField::checkNumber $outs "%- 15d" msg]
  }

  if {$msg != "" } {
    set Info(messageIcon) error
    Message::dispOkMessage "Steady state output interval:\n$msg" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  if {$outs != "" && $outs < 0} {
    set Info(messageIcon) error
    Message::dispOkMessage "Steady state output interval cannot be < 0" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  set Timestep(STEADY_STATE_OUTPUT_INTERVAL) $outs

  #-Max steady iterations
  #
  set iters [ string trim $Timestep(STEADY_STATE_MAX_ITERATIONS) ]
  set msg ""
  if { $iters != "" } {
    set iters [DataField::checkNumber $iters "%- 15d" msg]
  }

  if {$msg != "" } {
    set Info(messageIcon) error
    Message::dispOkMessage "Steady state max iterations:\n$msg" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  if { $iters == "" || $iters <= 0 } {
    set Info(messageIcon) error
    Message::dispOkMessage "Steady state max iterations must be > 0" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  set Timestep(STEADY_STATE_MAX_ITERATIONS) $iters

  return 1
}


# Return 1 = ok, 0 = error
#
proc Timestep::checkTimestepEntry {} {
  global Info Timestep Model

  # Pick timestep entry and check it
  set values [string trimleft $Timestep(timestepEntry)]

  #-Three values are needed!
  if {3 != [llength $values] } {
    set Info(messageIcon) error
    Message::dispOkMessage "Enter three values!" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
      return 0
  }

  set size  [lindex $values 0]
  set steps [lindex $values 1]
  set outs  [lindex $values 2]

  #-Timestep size
   set value [DataField::checkNumber $size "%- g" msg]

  if {$msg != "" } {
    set Info(messageIcon) error
    Message::dispOkMessage $msg \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }
  if {$value <= 0} {
    set Info(messageIcon) error
    Message::dispOkMessage "Timestep size must be > 0" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  #-Timesteps in the interval
  set value [DataField::checkNumber $steps "%- 15d" msg]
  if {$msg != "" } {
    set Info(messageIcon) error
    Message::dispOkMessage $msg \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }
  if {$value <= 0} {
    set Info(messageIcon) error
    Message::dispOkMessage "Nof timesteps in the interval must be > 0" \
                "$Info(FRONT_NAME) error message" \
                $Timestep(winName)
    return 0
  }

  if {$value >= 100000} {
    set Info(messageIcon) error
    Message::dispOkMessage "Nof timesteps in the interval must be < 100000" \
                "$Info(FRONT_NAME) error message" \
                $Timestep(winName)
    return 0
  }

  #-Output steps in the interval
  set value [DataField::checkNumber $outs "%- 15d" msg]
  if {$msg != "" } {
    set Info(messageIcon) error
    Message::dispOkMessage $msg \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  if {$value <= 0} {
    set Info(messageIcon) error
    Message::dispOkMessage "Output step interval must be > 0" \
                  "$Info(FRONT_NAME) error message" \
                  $Timestep(winName)
    return 0
  }

  if {$value >= 100000} {
    set Info(messageIcon) error
    Message::dispOkMessage "Output step interval must be < 100000" \
                "$Info(FRONT_NAME) error message" \
                $Timestep(winName)
    return 0
  }

  return 1
}


proc Timestep::copyTimestepListRow {} {
  global Timestep Info

  #Timestep(timestepEntryWidget)
  #Timestep(timestepListLB)

  set index [$Timestep(timestepListLB) curselection]

  if {$index == ""} {
    return
  }

  set list_row [join [$Timestep(timestepListLB) get $index $index] ]

  set old_row [$Timestep(timestepEntryWidget) get]
  set old_row [string trim $old_row]

  if { $old_row != "" } {
    set msg "NOTE: Current entry data will be replaced by the selected row!"
    if { ![Message::verifyAction $Info(noviceUser) $msg ok warning] } {
      return
    }
  }

  set new_row ""
  foreach i {0 1 2} {
    append new_row [lindex $list_row $i]
    append new_row " "
  }

  $Timestep(timestepEntryWidget) delete 0 end
  $Timestep(timestepEntryWidget) insert 0 $new_row

}


proc Timestep::defineTimesteps {} {
}

proc Timestep::defineTimestepSize {} {
}


# end ecif_tk_timestepDefPanel.tcl
# ********************
