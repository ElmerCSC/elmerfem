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
#Module:    ecif_tk_modelInfoPanel.tcl
#Language:  Tcl
#Date:      23.10.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for general model information.
#
#************************************************************************
 
#------Model info proc------
#
# This procedure displays the model statistic-info screen
# NOTE: This is an info panel, but Cad input file and Mesh input file
# names are editable here. These names are however just for info, changing
# them does not affect anything in the model. This feature is provided
# mainly for the possibilty to keep this info sensible when model files are
# moved from a machine/platform to an other.
#
proc ModelInfo::openPanel { } {
  # Global variables
  global Info Model ModelProperty ModelInfo

  set w $ModelInfo(winName)
  set wgeom $ModelInfo(winGeometry)

  set id [winfo atom $w]
  set ModelInfo(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow ModelInfo $id $ModelInfo(winTitle) $wgeom] } {
    return
  }  
 
  # Create local value for this panel (for backup etc.)
  set ModelInfo(cadPath)  $Model(cadPath)
  set ModelInfo(meshPath) $Model(meshPath)
  set ModelInfo(created)  $Model(created)

  set ModelInfo(dataChanged) 0
  set ModelInfo(dataModified) 0

  toplevel $w
  focus $w
  set this $w

  wm title $w $ModelInfo(winTitle)
  wm geometry $w $wgeom 

  #-Frame padding parameters
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  set Model(nofMeshElements_) [expr $Model(nofMeshElements) - $Model(nofMeshSplittedElements)]

  #-Entry descriptions:
  # Lbl, Lbl width, Array, Variable, Editable flag, Width, Left justify flag

  set x_lbl [Util::getCoordinateLabel 1]
  set y_lbl [Util::getCoordinateLabel 2]
  set z_lbl [Util::getCoordinateLabel 3]

  set fieldInfos ""
  lappend fieldInfos "\"Cad input file:\"             12  Model cadPath             1 55 1"
  lappend fieldInfos "\"Mesh input file:\"            12  Model meshPath            1 55 1"
  lappend fieldInfos "\"Model created:\"              12  Model created             1 55 1"
  lappend fieldInfos "\"Model file:\"                 12  Model EMF_PATH            0 55 1"
  lappend fieldInfos "\"Model directory:\"            12  ModelProperty MODEL_DIRECTORY,absolute   0 55 1"
  lappend fieldInfos "\"Include path:\"               12  ModelProperty INCLUDE_PATH,absolute      0 55 1"
  lappend fieldInfos "\"Results dir:\"                12  ModelProperty RESULTS_DIRECTORY,absolute 0 55 1"
  lappend fieldInfos "\"Log dir:\"                    12  ModelProperty LOG_DIRECTORY,absolute     0 55 1"
  lappend fieldInfos "\"Nof bodies:\"                 24  Model nofBodies           0 20 0"
  lappend fieldInfos "\"Nof boundaries:\"             24  Model nofElements         0 20 0"
  lappend fieldInfos "\"Nof outer boundaries:\"       24  Model nofOuterBoundaries  0 20 0"
  lappend fieldInfos "\"Nof inner boundaries:\"       24  Model nofInnerBoundaries  0 20 0"
  lappend fieldInfos "\"Nof vertices:\"               24  Model nofVertices         0 20 0"
  lappend fieldInfos "\"Model  $x_lbl dimensions \[m\]:\"  24  Model {minX maxX}         0 20 0"
  lappend fieldInfos "\"Model  $y_lbl dimensions \[m\]:\"  24  Model {minY maxY}         0 20 0"
  lappend fieldInfos "\"Model  $z_lbl dimensions \[m\]:\"  24  Model {minZ maxZ}         0 20 0"
  lappend fieldInfos "\"Minimum edge size:\"          24  Model minEdgeSize         0 20 0"
  lappend fieldInfos "\"Current mesh:\"               12  ModelProperty CURRENT_MESH_DIRECTORY,absolute   0 55 1"
  lappend fieldInfos "\"Nof mesh nodes:\"             24  Model nofMeshNodes        0 20 0"
  lappend fieldInfos "\"Nof mesh elements:\"          24  Model nofMeshElements_    0 20 0"
  lappend fieldInfos "\"Nof mesh boundary elements:\" 24  Model nofMeshBndrElements 0 20 0"

  # Create labels and variable entries
  # ==================================

  set index 0
  foreach fi $fieldInfos {

    incr index

    #--Frame 
    set m [frame $w.f$index]

    #--Label
    set txt [lindex $fi 0]
    set wid [lindex $fi 1]
    label $m.l -width $wid -text $txt -anchor w

    # Pack frame and label
    pack $m -side top -expand 0 -fill x -padx $fpx1 -pady $fpy2
    pack $m.l -side left -anchor w -in $m 

    #--Entries
    set arr   [lindex $fi 2]
    set vars  [lindex $fi 3]
    set edit  [lindex $fi 4]
    set wid   [lindex $fi 5]
    set ljust [lindex $fi 6]

    if {$edit} {
      set state normal
      set bgcolor white
      set relief $Info(activeRelief)
    } else {
      set state disabled
      set bgcolor $Info(nonActiveBgLight)
      set relief $Info(nonActiveRelief)
    }

    if {$ljust} {
      set just left
    } else {
      set just center
    }


    set count 0
    foreach var $vars {
      incr count

      set txtvar $arr

      append txtvar "($var)"    

      set wdg [entry $m.d$count -width $wid -relief $relief -textvariable $txtvar \
                                -justify $just -state $state  -bg $bgcolor ]
      
      set ModelInfo(allWidgets,$var) $wdg

      if { $state == "normal" } {
        bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 ModelInfo"
      }
                        
      pack $m.d$count -side left -in $m -padx $fpx1
    }
  }


  # Create buttons
  # ===============

  frame $w.fB

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  #-Ok button
  set ok_btn [button $w.fB.ok -text OK -command "ModelInfo::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "ModelInfo::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command "ModelInfo::panelApply" \
                                 -state $ap]

  focus $ok_btn
  set ModelInfo(applyButton)  $ap_btn
  set ModelInfo(cancelButton) $cn_btn

  pack $w.fB
  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1
}


# This does not actually save anything, just updates
# modified states
# Changes are stored directly to Model-array vars
#
proc ModelInfo::panelSave { {inform_front 1} } {
  global ModelInfo Model

  # Backup values
  set ModelInfo(cadPath)  $Model(cadPath)
  set ModelInfo(meshPath) $Model(meshPath)

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 ModelInfo 
  Panel::panelDataModified 0 ModelInfo 

  Util::updateMainWindowTitle
}


proc ModelInfo::panelOk {w} {
  global ModelInfo

  #---No changes
  if { !$ModelInfo(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![ModelInfo::panelCheck] } {
    return
  }

  ModelInfo::panelSave
  Panel::cancel $w
}

proc ModelInfo::panelApply {} {
  global ModelInfo

  #---No changes
  if { !$ModelInfo(dataChanged) } {
    return
  }

  #---Error in data
  if { ![ModelInfo::panelCheck] } {
    return
  }

  ModelInfo::panelSave
}

proc ModelInfo::panelCancel {w} {
  global ModelInfo Model

  if { ![Panel::verifyCancel ModelInfo] } {
    return
  }

  # Restore values
  set Model(cadPath)  $ModelInfo(cadPath)
  set Model(meshPath) $ModelInfo(meshPath)
  set Model(created)  $ModelInfo(created)

  Panel::cancel $w
}


proc ModelInfo::panelCheck {} {
  return 1
}

# end ecif_tk_modelInfoPanel.tcl
# ********************
