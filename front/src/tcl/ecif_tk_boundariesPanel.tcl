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
#Module:    ecif_tk_boundariesPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining new (mesh) boundaries and boundary names
#
#************************************************************************
 
#------Panel for editing boundaries (split, combine)
#

###########################
### Screen construction ###
### and initialization  ###
###########################

# Screen is divided into three main parts:
# 1: Objects listbox (fObjects), Elements listbox (fElements)
# 2: Split, Combine etc. buttons (fButtons1)
# 3: Ok and cancel etc. buttons (fButtons2)


proc Boundary::openPanel {} {
  global Boundary Info Common Model ModelFlags ObjectTable

  set w $Boundary(winName)
  set wgeom $Boundary(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set Boundary(winId) $id

  #--Store last created window info in Info
  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow Boundary $id $Boundary(winTitle) $wgeom] } {
    return
  }  

 if { $ModelFlags(GEOMETRY_TYPE_CAD) ||
      ([llength $Model(meshNames)] > 1)
    } {
    set Boundary(acceptEdit) 0
    
  } else {
    set Boundary(acceptEdit) 1
  }

  set Boundary(dataChanged) 0
  set Boundary(dataModified) 0
  
  set Common(currentArray) Boundary

  set this $w
  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $Boundary(winTitle)
  wm geometry $w $wgeom

  set Boundary(objectName) ""

  # Currently we can edit only Mesh geometry 
  #Util::setFlagValue GEOMETRY_TYPE GEOMETRY_TYPE_CAD 0
  #Util::setFlagValue GEOMETRY_TYPE GEOMETRY_TYPE_MESH 1

  # We are editing boundary objects
  if {$Model(GEOMETRY_DIMENSION) == "2D"} {
    Util::setFlagValue DRAW_TARGET DRAW_TARGET_EDGES 1 

  } else {
    Util::setFlagValue DRAW_TARGET DRAW_TARGET_SURFACES 1
  }

  Util::setFlagValue SELECT_METHOD SELECT_METHOD_SINGLE 1
  Util::setFlagValue SELECT_METHOD SELECT_METHOD_BY_NEIGHBOR 1
  Util::setFlagValue SELECT_MODE SELECT_MODE_EXTEND 1
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_TOGGLE 1

  Boundary::setEditBoundaries 1

  set Boundary(currentSplitCombineIndex) 0
  set Boundary(maxSplitCombineIndex) 0

  #----WIDGET DEFINITION
  #
  #---Outer frames
  frame $w.f1 
  frame $w.f2 
  frame $w.bf1 
  frame $w.bf2 
  frame $w.bf3 
  frame $w.bf4 

  #-Object list-boxes frame 
  frame $w.f1.fObjectBox
  frame $w.f1.fObjectBox.fObjects 
  frame $w.f1.fBoundaryBox
  frame $w.f1.fBoundaryBox.fBoundaries

  #-Body and bounadry names frame
  frame $w.f2.fNames 

  #-Buttons-1 frame (Split Combine Update)
  frame $w.bf1.fButtons 

  #-Buttons-2 frame (Select method: Select by normal, Select by plane, Tolerances)
  frame $w.bf2.fButtons 

  #-Buttons-3 frame (Select, Undo, Redo)
  frame $w.bf3.fButtons 

  #-Buttons-4 frame (Cancel, Ok)
  frame $w.bf4.fButtons 

  #---Frame internals
  #-Object list-box
  StdPanelCreate::createObjectBoxArea $w.f1.fObjectBox Boundary
  
  #-Element list-box
  StdPanelCreate::createBoundaryBoxArea $w.f1.fBoundaryBox Boundary
  
  #-Body name field
  Boundary::createNamesArea $w.f2.fNames
  
  #-Buttons-1
  Boundary::createButtons1Area $w.bf1.fButtons
  #-Buttons-2
  Boundary::createButtons2Area $w.bf2.fButtons
  #-Buttons-3
  Boundary::createButtons3Area $w.bf3.fButtons
  #-Buttons-4
  Boundary::createButtons4Area $w.bf4.fButtons $this
  #
  #-end WIDGET DEFINITION

  #---WIDGET BINDINGS and DATA INITIALIZATION
  Boundary::packPanel $w
  StdPanelInit::createWidgetBindings Boundary
  StdPanelInit::initPanelData Boundary
  #Boundary::packPanel $w

  Boundary::backupData
}
#-end openPanel proc


proc Boundary::packPanel {w} {
  global Boundary Info

  #---All frames
  set fpx0 0
  set fpy0 0
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)
  set expnd 1
  set fill both
  set anch w

  #NOTE: We give highest priority for button2 frame (=f4) to keep Ok/Cancel
  #buttons visible
  #---All frames
  set frames [list $w.bf4 $w.bf3 $w.bf2 $w.bf1 $w.f2 $w.f1 ]
  foreach fr $frames {
    pack $fr -side bottom -expand $expnd -fill $fill -padx $fpx2 -pady $fpy2
  }

  #---Individual frames
  #-Object list-box frame
  pack $w.f1.fObjectBox -side left -expand $expnd -fill $fill -anchor $anch \
                        -pady $fpy2

  #-Element list-box frame
  pack $w.f1.fBoundaryBox -side left -expand $expnd -fill $fill -anchor $anch \
                          -pady $fpy2

  #-Body name frame
  pack $w.f2.fNames -side left -expand $expnd -fill $fill -anchor $anch \
                    -pady $fpy2

  #-Buttons-1 frame
  pack $w.bf1.fButtons -expand $expnd -fill $fill -pady $fpy1

  #-Buttons-2 frame
  pack $w.bf2.fButtons -expand $expnd  -pady $fpy1
  
  #-Buttons-3 frame
  pack $w.bf3.fButtons -expand $expnd -pady $fpy1

  #-Buttons-4 frame
  pack $w.bf4.fButtons -expand $expnd -pady $fpy1

  #---Frame internals
  #-Object list-box
  StdPanelCreate::packObjectBoxArea $w.f1.fObjectBox Boundary

  #-Element list-box
  StdPanelCreate::packBoundaryBoxArea $w.f1.fBoundaryBox Boundary

  #-Body name 
  Boundary::packNamesArea $w.f2.fNames

  #-Buttons-1
  Boundary::packButtons1Area $w.bf1.fButtons

  #-Buttons-2
  Boundary::packButtons2Area $w.bf2.fButtons

  #-Buttons-3
  Boundary::packButtons3Area $w.bf3.fButtons

  #-Buttons-4
  Boundary::packButtons4Area $w.bf4.fButtons

  #-end SCREEN PACKING
}


##########################
#### Helper functions ####
##########################

#--Body name area
#--Creating
proc Boundary::createNamesArea {frame} {
  global Boundary Common Info

  # Body name 
  label $frame.body_name_label -text "Body name:" -width 10
  entry $frame.body_name -textvariable Boundary(objectName) -font $Info(entryFont) -state disabled -bg $Info(nonActiveBg)

  # Boundary name 
  label $frame.bndr_name_label -text "Boundary name: " -width 15
  entry $frame.bndr_name_entry -textvariable Boundary(boundaryName) -font $Info(entryFont)

  set wdg $frame.bndr_name_entry
  set  Boundary(boundaryNameWidget) $wdg

  # Set bindings
  set Boundary(boundaryName,err) 0
  set Boundary(boundaryName,mod) 0
  #bind $wdg <KeyRelease> "Boundary::applyBoundaryNameEntryContents"
  #bind $wdg <KeyRelease> "+Panel::panelDataChanged 1 Boundary $wdg {%A %K}"
  bind $wdg <KeyRelease> "+Panel::panelDataModified 1 Boundary $wdg {%A %K}"
  bind $wdg <KeyRelease> "+Widget::entryKeyRelease Boundary boundaryName $wdg {%A %K}"
  bind $wdg <KeyPress-Escape> "+Widget::entryKeyPress-Escape Boundary boundaryName $wdg {%A %K}"

}

#--Packing
proc Boundary::packNamesArea {frame} {
  global Boundary Common Info

  #-Body name entry
  pack $frame.body_name_label  -side left -pady 1  -padx 3 -anchor w
  pack $frame.body_name  -side left -pady 1  -padx 3 -anchor w 

  #-Boundary name entry
  pack $frame.bndr_name_label  -side left -pady 0  -padx 3 -anchor w 
  pack $frame.bndr_name_entry  -side left -pady 0  -padx 3 -anchor w 
}


#-----Buttons area

#--Creating
# Slit/Combine buttons
#
proc Boundary::createButtons1Area {frame} {
  global Boundary Common Info ModelFlags

  if { $Boundary(acceptEdit) } {
    set state normal
  } else {
    set state disabled
  }

  # Split and Combine buttons
  #
  set wdg [button $frame.split -text Split -command "Boundary::split_" \
                               -state $state]
  bind $wdg <ButtonPress-1> "Panel::panelDataChanged 1 Boundary $wdg {%A %K}"


  set wdg [button $frame.combine -text Combine -command "Boundary::combine" \
                                 -state $state]
  bind $wdg <ButtonPress-1> "Panel::panelDataChanged 1 Boundary $wdg {%A %K}"


  # Split/Combine Undo/Redo buttons
  #
  set wdg [button $frame.sac_undo -text "Undo" \
                                  -command "Boundary::splitCombineUndo" \
                                  -state disabled]
  
  set Boundary(splitCombineUndoButton) $wdg

  set wdg [button $frame.sac_redo -text "Redo" \
                                  -command "Boundary::splitCombineRedo" \
                                  -state disabled]

  set Boundary(splitCombineRedoButton) $wdg


  # Update boundary data button
  #
  set wdg [button $frame.update -text "Update" \
                                -command "Boundary::update" \
                                -state disabled]

  set Boundary(panelUpdateButton) $wdg

  bind $wdg <ButtonPress> "+Panel::panelDataChanged 1 Boundary $wdg {%A %K}"
  bind $wdg <ButtonPress> "+Panel::panelDataModified 0 Boundary $wdg {%A %K}"



  #-Store split/combine buttons to panel specific array
  set Boundary(buttons1) [list $frame.split $frame.combine \
                               $frame.sac_undo $frame.sac_redo]

}


#--Packing
# Split/Combine buttons
#
proc Boundary::packButtons1Area {frame} {
  global Boundary Common Info

  set expnd 0
  set fill both

  # Split and Combine buttons
  set px 5; set py 2
  pack $frame.split    -side left -pady $py -padx $px -anchor w
  pack $frame.combine  -side left -pady $py -padx $px -anchor w

  frame $frame.separator
  pack $frame.separator -side left -pady $py -padx $px -anchor w

  pack $frame.sac_undo -side left -pady $py -padx $px -anchor w
  pack $frame.sac_redo -side left -pady $py -padx $px -anchor w

  pack $frame.update -side right -pady $py -padx $px -anchor e
}


#--Creating
# Select method buttons, eometry discretization buttons
#
proc Boundary::createButtons2Area {frame} {
  global Boundary Common Info ModelFlags
  
  set b_wid1  9
  set b_wid2  9
  set b_wid3  12
  set e_wid   5
  set e_wid2  12
  set l_wid  19
  set l_wid2 14
  set s_wid  13
  set s_len 150

  frame $frame.buttons

  #---Selection method and mode buttons
  #
  frame $frame.buttonsS
  frame $frame.buttonsS.methods -relief groove -bd 2
  frame $frame.buttonsS.modes -relief groove -bd 2

  if { $Boundary(acceptEdit) } {
    set state normal
    set bg white
  } else {
    set state disabled
    set bg $Info(nonActiveBg)
  }

  set dsc_state disabled
  set dsc_bg $Info(nonActiveBg)

  if { $ModelFlags(GEOMETRY_TYPE_CAD) == 0 } {
    set fn_state disabled
    set fn_bg $Info(nonActiveBg)
  } else {
    set fn_state normal
    set fn_bg white
  }

  #-Method buttons
  set m $frame.buttonsS.methods

  label $m.label -text "Selection method:"

  radiobutton $m.by_neighbor -width $b_wid1  -text "By neighbor" \
    -variable ModelFlags(SELECT_METHOD_BY_NEIGHBOR) -value 1      \
    -command "Boundary::setSelectMethodByNeighbor" \
    -state $state -anchor w

  radiobutton $m.by_normal -width $b_wid1  -text "By normal" \
    -variable ModelFlags(SELECT_METHOD_BY_NORMAL) -value 1      \
    -command "Boundary::setSelectMethodByNormal" \
    -state $state -anchor w

  radiobutton $m.by_plane -width $b_wid1 -text "By plane" \
    -variable ModelFlags(SELECT_METHOD_BY_PLANE) -value 1      \
    -command "Boundary::setSelectMethodByPlane" \
    -state $state -anchor w

  radiobutton $m.all -width $b_wid1  -text "All" \
    -variable ModelFlags(SELECT_METHOD_ALL) -value 1      \
    -command "Boundary::setSelectMethodAll" \
    -state $state -anchor w


  #---Selection mode buttons
  #
  set m $frame.buttonsS.modes

  label $m.label -text "Selection mode:"

  radiobutton $m.extend -width $b_wid2 -text "Extend" -indicatoron 1 \
    -variable ModelFlags(SELECT_MODE_EXTEND) -value 1      \
    -command "Boundary::setSelectModeExtend" \
    -state $state -anchor w

  radiobutton $m.reduce -width $b_wid2 -text "Reduce" -indicatoron 1 \
    -variable ModelFlags(SELECT_MODE_REDUCE) -value 1      \
    -command "Boundary::setSelectModeReduce" \
    -state $state -anchor w

  radiobutton $m.toggle -width $b_wid2 -text "Toggle" -indicatoron 1 \
    -variable ModelFlags(SELECT_MODE_TOGGLE) -value 1      \
    -command "Boundary::setSelectModeToggle" \
    -state $state -anchor w

  #-Discretization buttons and entries
  #
  frame $frame.buttonsD -relief groove -bd 2
  frame $frame.buttonsD.types 
  frame $frame.buttonsD.entries

  label $frame.buttonsD.label -text "Discretization settings:"

  set m $frame.buttonsD.entries
  
  frame $m.fT
  frame $m.fV
  frame $m.fF

  label $m.fT.label -width 16 -text "Delta types (H,N,-):"
  entry $m.fT.eT -width 13  -textvariable Boundary(DELTA_TYPES) \
    -font $Info(entryFont) -state $dsc_state -bg $dsc_bg

  label $m.fV.label -width 16 -text "Delta values:"
  entry $m.fV.eV -width 13 -textvariable Boundary(DELTA_VALUES) \
    -font $Info(entryFont) -state $dsc_state -bg $dsc_bg

  label $m.fF.label -width 16 -text "Fixed N in meshing"
  entry $m.fF.eF -width 13 -textvariable Boundary(FIXED_N) \
    -font $Info(entryFont) -state $fn_state -bg $fn_bg

  set wdg $m.fT.eT
  set Boundary(allWidgets,DELTA_TYPES) $wdg
  bind $wdg <KeyRelease> "Panel::panelDataModified 1 Boundary $wdg {%A %K}"
  #bind $wdg <KeyRelease> "Panel::panelDataChanged 1 Boundary"
  Widget::setEntryBindings "standard" Boundary DELTA_TYPES $m.fT.eT
  set Boundary(DELTA_TYPES,err) 0
  set Boundary(DELTA_TYPES,mod) 0

  set wdg $m.fV.eV
  set Boundary(allWidgets,DELTA_VALUES) $wdg
  bind $wdg <KeyRelease> "Panel::panelDataChanged 1 Boundary $wdg {%A %K}"
  bind $wdg <KeyRelease> "Panel::panelDataModified 1 Boundary $wdg {%A %K}"
  Widget::setEntryBindings "standard" Boundary DELTA_VALUES $wdg
  set Boundary(DELTA_VALUES,err) 0
  set Boundary(DELTA_VALUES,mod) 0

  set wdg $m.fF.eF
  set Boundary(allWidgets,FIXED_N) $wdg
  bind $wdg <KeyRelease> "Panel::panelDataChanged 1 Boundary $wdg {%A %K}"
  bind $wdg <KeyRelease> "Panel::panelDataModified 1 Boundary $wdg {%A %K}"
  Widget::setEntryBindings "standard" Boundary FIXED_N $wdg
  set Boundary(FIXED_N,err) 0
  set Boundary(FIXED_N,mod) 0

  #---Tolerance scales
  #
  frame $frame.tolerances
  frame $frame.tolerances.scales
  frame $frame.tolerances.scales.normal
  frame $frame.tolerances.scales.distance
  frame $frame.tolerances.nextActive

  set m $frame.tolerances.scales.normal
  set resol 0.25 
  label $m.label -width $l_wid -text "Normal tol. (Deg 0-$Info(maxNormalTolerance)):"

  entry $m.entry -width $e_wid  -textvariable Info(normalTolerance) \
    -font $Info(entryFont) -state $state -bg $bg

  scale $m.scale  -orient horizontal -width $s_wid      \
                  -variable Info(normalTolerance)       \
                  -from 0 -to $Info(maxNormalTolerance) \
                  -digits 4 -showvalue 0                \
                  -resolution $resol                    \
                  -font $Info(entryFont) -state $state


  set m $frame.tolerances.scales.distance
  set resol 0.01
  label $m.label -width $l_wid -text "Distance tol. (% 0-$Info(maxDistanceTolerance)):"

  entry $m.entry -width $e_wid  -textvariable Info(distanceTolerance) \
    -font $Info(entryFont) -state $state  -bg $bg

  scale $m.scale  -orient horizontal -width $s_wid        \
                  -variable Info(distanceTolerance)       \
                  -from 0 -to $Info(maxDistanceTolerance) \
                  -digits 3 -showvalue 0                  \
                  -resolution $resol                      \
                  -font $Info(entryFont) -state $state

  set m $frame.tolerances.nextActive

  label $m.label  -width $l_wid -text "Next active:"

  entry $m.entry  -width $e_wid2 \
                  -textvariable Info(NEXT_ACTIVE_SELECTION_TOLERANCE) \
                  -font $Info(entryFont) -state $state -bg $bg
}


#--Packing
proc Boundary::packButtons2Area {frame} {
  global Boundary

  set expnd 1
  set fill none

  pack $frame.buttons $frame.tolerances -side top -pady 2

  pack $frame.buttonsS  -side left -anchor w
  pack $frame.buttonsD -side left -padx 6m -anchor e
  
  #---Method and mode buttons
  #
  pack $frame.buttonsS.methods -side left -anchor w
  pack $frame.buttonsS.modes -side left -padx 2m -anchor e

  set m $frame.buttonsS.methods
  pack $m.label \
       $m.by_neighbor $m.by_normal $m.by_plane $m.all \
       -side top -anchor w -padx 1m

  set m $frame.buttonsS.modes
  pack $m.label \
       $m.extend  $m.reduce $m.toggle \
       -side top -anchor w -padx 1m
  
  #---Discretization buttons and entries
  #
  pack $frame.buttonsD.label $frame.buttonsD.entries \
       -side top
  
  #-Label + entry 
  set m $frame.buttonsD.entries

  pack $m.fT $m.fV $m.fF -side top -anchor w -pady 1m

  pack $m.fT.label $m.fT.eT -side left -anchor w -padx 1m
  pack $m.fV.label $m.fV.eV -side left -anchor w -padx 1m
  pack $m.fF.label $m.fF.eF -side left -anchor w -padx 1m

  #---Tolerance scales and entries
  #
  pack $frame.tolerances.scales $frame.tolerances.nextActive -side left

  set m $frame.tolerances.scales
  pack $m.normal $m.distance -side top -pady 1m
 
  set m $frame.tolerances.scales.normal
  pack $m.label $m.entry $m.scale -side left -padx 1m

  set m $frame.tolerances.scales.distance
  pack $m.label $m.entry $m.scale -side left -padx 1m

  set m $frame.tolerances.nextActive
  pack $m.label $m.entry -side top -anchor n -pady 1m
}


#--Creating
# Select Do/Undo buttons
#
proc Boundary::createButtons3Area {frame} {
  global Boundary ModelFlags Info

  if { $Boundary(acceptEdit) } {
    set state normal
  } else {
    set state disabled
  }

  set wdg [button $frame.select_do -text "Select" \
            -command "Boundary::selectMeshBoundaryElements do" \
            -state $state]


  set wdg [button $frame.select_undo -text "Undo" \
            -command "Boundary::selectMeshBoundaryElements undo" \
            -state disabled]

  set Boundary(selectUndoButton) $wdg

  set wdg [button $frame.select_redo -text "Redo" \
            -command "Boundary::selectMeshBoundaryElements redo" \
            -state disabled]

  set Boundary(selectRedoButton) $wdg
}


#--Packing
# Select Do/Undo buttons
#
proc Boundary::packButtons3Area {frame} {
  global Boundary Common Info

  set expnd 1
  set fill none
  pack $frame.select_do $frame.select_undo $frame.select_redo -side left -padx 4m 
}


#--Creating
proc Boundary::createButtons4Area {frame this} {
  global Boundary Common Info Model

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $frame.ok -text OK -command "Boundary::panelOk $this"]
  set cn_btn [button $frame.cancel -text Cancel -command "Boundary::panelCancel $this" \
                                   -state $ca]
  set ap_btn [button $frame.apply -text Apply -command "Boundary::panelApply" \
                                  -state $ap]

  focus $ok_btn

  set Boundary(applyButton)  $ap_btn
  set Boundary(cancelButton) $cn_btn
}


#--Packing
proc Boundary::packButtons4Area {frame} {
  global Common Info

  set expnd 1
  set fill none
  pack $frame.ok $frame.cancel $frame.apply -side left -padx 4m -pady 2m
}



################################
# BOUNDADRY SLITTING/COMBINING #
################################

# Combine boundaries
#
proc Boundary::combine {} {
  global Info Boundary ObjectTable

  # Save names etc.
  Boundary::saveNames

  # Combine boundaries
  set data ""
  append data $Boundary(body1Id)
  append data " "
  append data $Boundary(body2Id)

	set Boundary(boundaryIndex) ""

  Util::cpp_exec combineBoundaries $data
	
  # Update panel with the newly set data
  Boundary::panelUpdate

  incr Boundary(currentSplitCombineIndex)

  Boundary::splitCombineUpdate
}


# Split boundaries
#
proc Boundary::split_ {} {
  global Info Boundary ObjectTable

  # Save names.
  Boundary::saveNames

  # Split boundary 
  set data ""
  append data $Boundary(body1Id)
  append data " "
  append data $Boundary(body2Id)

  Util::cpp_exec splitBoundary $data

  # Update panel with the newly set data
  Boundary::panelUpdate

  incr Boundary(currentSplitCombineIndex)

  Boundary::splitCombineUpdate
}


# Update split/combine data after slpit/combine
#
proc Boundary::splitCombineUpdate {} {
  global Boundary ModelFlags

  # Boundary geometry was edited
  if { $Boundary(currentSplitCombineIndex) > 0 } {
    set ModelFlags(GEOMETRY_EDITED_BOUNDARIES) 1
  } 

  if { $Boundary(maxSplitCombineIndex) < 
       $Boundary(currentSplitCombineIndex)
     } {
    set Boundary(maxSplitCombineIndex) $Boundary(currentSplitCombineIndex)
  }

  if { $Boundary(currentSplitCombineIndex) == 0 } {
    $Boundary(splitCombineUndoButton) configure -state disabled
  } else {
    $Boundary(splitCombineUndoButton) configure -state normal
  }

  if { $Boundary(currentSplitCombineIndex) == $Boundary(maxSplitCombineIndex) } {
    $Boundary(splitCombineRedoButton) configure -state disabled
  } else { 
    $Boundary(splitCombineRedoButton) configure -state normal
  }
}


# Redo boundaries split/combine
#
proc Boundary::splitCombineRedo {} {
  global Boundary

  # Save names.
  Boundary::saveNames

  if { $Boundary(currentSplitCombineIndex) == 
       $Boundary(maxSplitCombineIndex)
     } {
    Message::showMessage "Nothing to redo!"
    return
  }

  Util::cpp_exec splitCombineBoundariesRedo

  incr Boundary(currentSplitCombineIndex)

  Boundary::splitCombineUpdate

  # Update panel with the newly set data
  Boundary::panelUpdate
}


# Undo boundaries split/combine
#
proc Boundary::splitCombineUndo {} {
  global Boundary

  # Save names.
  Boundary::saveNames

  if { $Boundary(currentSplitCombineIndex) == 0 } {
    Message::showMessage "Nothing to undo!"
    return
  }
 
  Util::cpp_exec splitCombineBoundariesUndo

  # Restore values
  incr Boundary(currentSplitCombineIndex) -1

  Boundary::splitCombineUpdate

  # Update panel with the newly set data
  Boundary::panelUpdate
}


# Update panel with the newly set data
# after split/combine or undo/redo
#
proc Boundary::panelUpdate {} {
  global Boundary

  StdPanelInit::initPanelData Boundary
  $Boundary(boundaryLB) see end

  if { $Boundary(currentSplitCombineIndex) > 0 } {
    Panel::panelDataChanged 1 Boundary
  } 

}


#########################
# BOUNDARY DISCRETATION #
#########################

proc Boundary::writeDiscretationData {} {
  global Boundary Model ModelFlags ObjectTable

  if { $ModelFlags(GEOMETRY_TYPE_CAD) == 0 } {
    return
  }

  Widget::configureEntry $Boundary(allWidgets,FIXED_N) Normal
  
  set id $Boundary(boundaryId)

  set cmp_cnt $ObjectTable($id,nofCmp)
  
  set tps $ObjectTable($id,dscTp)
  set dvals $ObjectTable($id,dscU)
  set useFN $ObjectTable($id,useFN)

  if { $tps != "" } {
    set st normal
  } else {
    set st disabled
  }

  #--Format delta values
  #
  set vals ""
  set idx 0
  foreach val $dvals {
    set tp [string index $tps $idx]
    incr idx

    if { $tp == "N" } {
      lappend vals [expr round($val)]
    } else {
      lappend vals $val
    }
  }

  #--Insert spaces between type characters
  set types ""
  for {set idx 0} {$idx < [string length $tps]} {incr idx} {
    append types [string index $tps $idx]
    append types " "
  }
  set types [string trim $types]

  # Set  field data
  #
  set Boundary(DELTA_TYPES) $types
  set Boundary(DELTA_VALUES) $vals
  set Boundary(FIXED_N) $useFN

  Widget::configureEntry $Boundary(allWidgets,DELTA_TYPES) $st
  Widget::configureEntry $Boundary(allWidgets,DELTA_VALUES) $st
}


proc Boundary::readDiscretationData {} {
  global Boundary Info Model ModelFlags ObjectTable

  if { $ModelFlags(GEOMETRY_TYPE_CAD) == 0 } {
    return
  }

  set id $Boundary(boundaryId)

  if { $ObjectTable($id,dscTp) != "" } {
    set is_discretized 1
  } else {
    set is_discretized 0
  }

  set types $Boundary(DELTA_TYPES)
  
  #--Check delta types
  # Remove spaces from string
  #
  set tps ""
  foreach chr $types {
    if { [string equal -nocase $chr "H"] ||
         [string equal -nocase $chr "N"] ||
         [string equal -nocase $chr "-"] 
       } {
      append tps $chr
    } elseif { $chr != " " } {
      set Info(messageIcon) error
      Message::dispOkMessage "Illegal value ($chr).  Delta types should be H,N or -" \
                  "$Info(FRONT_NAME) error message" \
                  $Boundary(winName)
      set Boundary(DELTA_TYPES,err) 1
      return 0
    }
  }
  
  #--Check item counts
  #
  set cmp_cnt $ObjectTable($id,nofCmp)

  set cnt [string length $tps]
  if { $is_discretized && $cnt != $cmp_cnt } {
    set Info(messageIcon) error
    Message::dispOkMessage "The nof delta types ($cnt) does not match the nof geometry components ($cmp_cnt)" \
             "$Info(FRONT_NAME) error message" \
              $Boundary(winName)
    set Boundary(DELTA_TYPES,err) 1
    return 0
  }

  set cnt [llength $Boundary(DELTA_VALUES)]
  if { $is_discretized && $cnt != $cmp_cnt } {
    set Info(messageIcon) error
    Message::dispOkMessage "The nof delta values ($cnt) does not match the nof geometry components ($cmp_cnt)" \
             "$Info(FRONT_NAME) error message" \
              $Boundary(winName)
    set Boundary(DELTA_VALUES,err) 1
    return 0
  }

  set cnt [llength $Boundary(FIXED_N)]
  if { $cnt != $cmp_cnt } {
    set Info(messageIcon) error
    Message::dispOkMessage "The nof fixed N meshing flags ($cnt) does not match the nof geometry components ($cmp_cnt)" \
             "$Info(FRONT_NAME) error message" \
              $Boundary(winName)
    set Boundary(FIXED_N,err) 1
    return 0
  }
  
  #--Check delta values
  #
  set vals ""
  set idx 0
  foreach val $Boundary(DELTA_VALUES) {
    
    set tp [string index $tps $idx]
    incr idx
    
    if { $chr != "-" && $val <= 0 } {
      set Info(messageIcon) error
      Message::dispOkMessage "Illegal value ($val).  Delta values should be positive for N and H" \
               "$Info(FRONT_NAME) error message" \
                $Boundary(winName)
      set Boundary(DELTA_VALUES,err) 1
      return 0
    }

    if { $chr == "N" } {
      lappend vals [expr round($val)]
    } else {
      lappend vals $val
    }
  }
  set Boundary(DELTA_VALUES) $vals

  #--Check fixed N flags
  #
  set vals ""
  set idx 0
  foreach val $Boundary(FIXED_N) {

    if { ($val != 0) && ($val != 1) } {
      set Info(messageIcon) error
      Message::dispOkMessage "Illegal value ($val).  Fixed N meshing flags should be 0 or 1" \
               "$Info(FRONT_NAME) error message" \
                $Boundary(winName)
      set Boundary(FIXED_N,err) 1
      return 0
    }

    lappend vals $val
  }
  set Boundary(FIXED_N) $vals

  # Store data in ObjectTable
  #
  if {  $is_discretized } {
    set ObjectTable($id,dscTp) $tps
    set ObjectTable($id,dscU) $Boundary(DELTA_VALUES)
  }

  set ObjectTable($id,useFN) $Boundary(FIXED_N)

  # Data Ok
  #
  return 1
}



##################
# SAVE etc procs #
##################


# Update button proc
#
proc Boundary::update {} {
  global Boundary

  if { ![Boundary::applyBoundaryNameEntryContents] } {
    return 0
  }

  return 1
}


# Save boundary names
#
proc Boundary::saveNames {} {
  global Boundary Model ObjectTable

  # Update name
  foreach id $Boundary(allTargetIds) {
  
    set name $ObjectTable($id,nm)

    set old_len [string length $name]
    set name [string trim $name]
    set new_len [string length $name]

    if { $name == "" } {
      set name "Boundary$id"
    }

    if { $old_len != $new_len || 
         $new_len == 0
       } {
      set ObjectTable($id,nm) $name
    }
  }

  Util::cpp_exec boundaryPanelNamesOk
}


proc Boundary::panelSave { {inform_front 1} } {
  global Boundary Model ModelFlags

  Boundary::saveNames

  # Boundary geometry was edited
  if { 1 || $Boundary(currentSplitCombineIndex) > 0 } {
    Util::cpp_exec boundariesPanelOk
  } 

  # Reset indices
  set Boundary(currentSplitCombineIndex) 0
  set Boundary(maxSplitCombineIndex) 0

  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }
  
  Panel::panelDataChanged 0 Boundary 
  Panel::panelDataModified 0 Boundary 
}


proc Boundary::panelOk {w} {
  global Boundary

  if { ![Panel::verifySave Boundary] } {
    return
  }

  #---No changes
  if { !$Boundary(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![Boundary::panelCheck] } {
    return
  }

  Boundary::panelSave

  Boundary::setEditBoundaries 0
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_TOGGLE 1
  Util::cpp_exec resetAllBoundarySelections
  Util::cpp_exec stopEditMeshBoundaries "ok"
  Util::cpp_exec storeBoundaryNames

  Panel::cancel $w
}


proc Boundary::panelApply {} {
  global Boundary

  #---No changes
  if { !$Boundary(dataModified) && !$Boundary(dataChanged) } {
    return
  }

  #---Error in data
  if { ![Boundary::panelCheck] } {
    return
  }

  Boundary::update

  Boundary::panelSave

  Boundary::backupData

  Util::cpp_exec resetAllBoundarySelections
  Util::cpp_exec stopEditMeshBoundaries "ok"
  Util::cpp_exec storeBoundaryNames
}


proc Boundary::panelCancel {w} {
  global Boundary 

  if { ![Panel::verifyCancel Boundary] } {
    return
  }

  Panel::cancel $w

  Panel::panelDataChanged 0 Boundary 
  Panel::panelDataModified 0 Boundary 

  Boundary::setEditBoundaries 0
  Util::setFlagValue SELECT_OBJECTS SELECT_OBJECTS_TOGGLE 1

  Boundary::restoreData

  Util::cpp_exec resetAllBoundarySelections
  Util::cpp_exec stopEditMeshBoundaries "cancel"
  Util::cpp_exec restoreBoundaryNames
}


proc Boundary::panelCheck {} {
  global Boundary Info ObjectTable
   
  # Check boundary names
  #
  foreach id1 $Boundary(allTargetIds) {

    foreach id2 $Boundary(allTargetIds) {

      if { $id2 <= $id1 } {
        continue
      }

      set name1 $ObjectTable($id1,nm)
      set name2 $ObjectTable($id2,nm)

      set tag1 $ObjectTable($id1,tg)
      set tag2 $ObjectTable($id2,tg)

      if { [string equal -nocase $name1 $name2] } {

        set msg [list "NOTE: Boundaries $tag1 and $tag2 have similar names:\n" \
                      "$name1 \n" \
                      "$name2 \n\n" \
                      $Info(anywayOk)]

        if { ![Message::verifyAction $Info(powerUser) $msg cancel question] } {
          return 0
        }
      }
    } ;# for id2
  } ;# for id1

  if { ![Boundary::readDiscretationData] } {
    return 0
  }

  # Ok
  return 1
}


###############
# OTHER procs #
###############


proc Boundary::backupData {} {
  global Boundary ModelFlags ObjectTable

  # The following data is relevant only for ad geometries
  #
  if { $ModelFlags(GEOMETRY_TYPE_CAD) == 0 } {
    return
  }

  foreach id $Boundary(boundaryIds) {
    set ObjectTable($id,dscTp,old) $ObjectTable($id,dscTp)
    set ObjectTable($id,dscU,old) $ObjectTable($id,dscU)
    set ObjectTable($id,useFN,old) $ObjectTable($id,useFN)
  }
}


proc Boundary::restoreData {} {
  global Boundary ModelFlags ObjectTable
  
  # The following data is relevant only for ad geometries
  #
  if { $ModelFlags(GEOMETRY_TYPE_CAD) == 0 } {
    return
  }

  foreach id $Boundary(boundaryIds) {
    set ObjectTable($id,dscTp) $ObjectTable($id,dscTp,old)
    set ObjectTable($id,dscU) $ObjectTable($id,dscU,old)
    set ObjectTable($id,useFN) $ObjectTable($id,useFN,old)
  }
}


# Boundary name entry field bind proc
# Trims trailing spaces before storing data
#
proc Boundary::applyDiscretationSettings {} {
  global Boundary Info ObjectTable

  set fields { DELTA_TYPES DELTA_VALUES FIXED_N}
  foreach fld $fields {
    set Boundary($fld,err) 0
    set Boundary($fld,mod) 0
    Widget::setEntryStatus Boundary $fld $Boundary(allWidgets,$fld)
  }

  return 1
}



# Boundary name entry field bind proc
# Trims trailing spaces before storing data
#
proc Boundary::applyBoundaryNameEntryContents {} {
  global Boundary Info ObjectTable

  set bid $Boundary(boundaryId)

  set name $Boundary(boundaryName) 

  # No empty field accepted!
  if { $name == "" } {
    set name "Boundary$bid"
    $Boundary(bounadryNameWidget) insert 0 $name
  }

  set Boundary(boundaryName) $name

  #-NOTE: Do not trim in the entry field, only when storing!
  set name [string trim $name]

  set ObjectTable($bid,nm) $name
  set Boundary(boundaryName,prev) $name
  set Boundary(boundaryName,old) $name

  
  #-Update boundary name in the listbox
  StdPanelInit::constructBoundaryLBList Boundary

  set sindex [$Boundary(boundaryLB) curselection]
  set nrow [StdPanelInit::formTargetLBRow Boundary $bid]
  ListBox::updateRow $Boundary(boundaryLB) $sindex $nrow

  set index $Boundary(currentSplitCombineIndex)

  set ObjectTable($bid,Nm,old,$index) $ObjectTable($bid,nm)

  set Boundary(boundaryName,err) 0
  set Boundary(boundaryName,mod) 0

  set Boundary(updateAllowed) 0
  Panel::panelDataModified 0 Boundary

  Widget::setEntryStatus Boundary boundaryName $Boundary(boundaryNameWidget)

  return 1
}


proc Boundary::setEditBoundaries {state} {

  Util::cpp_exec rendererSetEditBoundaries $state
}


proc Boundary::selectMeshBoundaryElements {type } {
  global Boundary Info ModelFlags

   # Update first tolerances in the model
   Util::cpp_exec readSelectionTolerances

   # Then select elements
   Util::cpp_exec selectMeshBoundaryElements $type

   set Info(lastAppliedSelectionMethod) $ModelFlags(SELECT_METHOD)
   set Info(lastAppliedSelectionMode) $ModelFlags(SELECT_MODE)

   $Boundary(selectRedoButton) configure -state normal
   $Boundary(selectUndoButton) configure -state normal
}


proc Boundary::setSelectMethodByNeighbor {} {
  Boundary::setNextActiveSelectionToleranceInfo SELECT_METHOD_BY_NEIGHBOR
  Util::setFlagValue SELECT_METHOD SELECT_METHOD_BY_NEIGHBOR 1
}


proc Boundary::setSelectMethodByNormal {} {
  Boundary::setNextActiveSelectionToleranceInfo SELECT_METHOD_BY_NORMAL
  Util::setFlagValue SELECT_METHOD SELECT_METHOD_BY_NORMAL 1
}


proc Boundary::setSelectMethodByPlane {} {
  Boundary::setNextActiveSelectionToleranceInfo SELECT_METHOD_BY_PLANE
  Util::setFlagValue SELECT_METHOD SELECT_METHOD_BY_PLANE 1
}


proc Boundary::setSelectMethodAll {} {
  Boundary::setNextActiveSelectionToleranceInfo SELECT_METHOD_ALL
  Util::setFlagValue SELECT_METHOD SELECT_METHOD_ALL 1
}


proc Boundary::setSelectModeExtend {} {
  Util::setFlagValue SELECT_MODE SELECT_MODE_EXTEND
}


proc Boundary::setSelectModeReduce {} {
  Util::setFlagValue SELECT_MODE SELECT_MODE_REDUCE
}


proc Boundary::setSelectModeToggle {} {
  Util::setFlagValue SELECT_MODE SELECT_MODE_TOGGLE
}


proc Boundary::setNextActiveSelectionToleranceInfo {selection_method} {
  global Info ModelFlags

  if { $Info(lastAppliedSelectionMethod) == $selection_method } {
      set Info(NEXT_ACTIVE_SELECTION_TOLERANCE) \
          $Info(nextActiveSelectionTolerance)
  } else {
      set Info(NEXT_ACTIVE_SELECTION_TOLERANCE) ""
  }
}


# end ecif_tk_boundariesPanel.tcl
# ********************

