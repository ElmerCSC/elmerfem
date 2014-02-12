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
#Module:    ecif_tk_coordinateDefPanel.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining the coordiante system
#
#************************************************************************


proc Coordinate::openPanel {} {
  ## This procedure displays the coordiante definition panel
  ## Global variables
  global Info Common Model Coordinate oldValues
  global platform
  upvar #0 Coordinate theArray

  set w $Coordinate(winName)
  set wgeom $Coordinate(winGeometry)

  set id [winfo atom $w]
  set Coordinate(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow Coordinate $id $Coordinate(winTitle) $wgeom] } {
    return
  }  

  set Coordinate(dataChanged) 0
  set Coordinate(dataModified) 0

  set this $w
  toplevel $w
  focus $w

  wm title $w $Coordinate(winTitle)
  wm geometry $w $wgeom 

  Panel::resetFields Coordinate

  Panel::initFields Coordinate

  set id $Coordinate(parameterId)
  DataField::formDataFields Coordinate $Coordinate($id,data)

  Panel::backupFields Coordinate
  
  set Info(coordiantePanelApplied) 0

  #-----WIDGET CREATION
  #-We organise data into four frames. These frames are then
  #-packed below each other and devided internally if needed.
  frame $w.f1 ;#--Coordinate geometry 
  frame $w.f1.cd  ;#--Coordinate related stuff
  frame $w.f1.cd.1 -bd 2 -relief groove ;#--Coordinate systems
  frame $w.f1.cd.2 -bd 2 -relief groove ;#--Coordinate mapping
  frame $w.fB ;#--Buttons
  
  # Coordinate system radiobuttons
  # ==============================
  #
  set fld COORDINATE_SYSTEM
  set m $w.f1.cd.1
  frame $m.lf
  label $m.lf.l  -anchor nw -text "Coordinate system: "

  # Add help event for the label
  set wdg $m.lf.l
  bind $wdg <ButtonPress-3> "Widget::genericButton-3 Coordinate $fld $wdg"

  frame $m.rf
  radiobutton $m.rf.r1  -text "Cartesian 2D"        -value "Cartesian 2D" 
  radiobutton $m.rf.r2  -text "Axi Symmetric"       -value "Axi Symmetric" 
  radiobutton $m.rf.r3  -text "Cylindric Symmetric" -value "Cylindric Symmetric" 
  radiobutton $m.rf.r4  -text "Polar 2D"            -value "Polar 2D" 
  radiobutton $m.rf.r5  -text "Cartesian 3D"        -value "Cartesian 3D" 
  radiobutton $m.rf.r6  -text "Cylindrical"           -value "Cylindrical"
  radiobutton $m.rf.r7  -text "Polar 3D"            -value "Polar 3D" 
  
  set Coordinate(allWidgets,$fld) ""

  #-Set button variables etc.
  set coordinate_system_ids {1 2 3 4 5 6 7}
  set is2D  {1 1 1 1 0 0 0}

  foreach i $coordinate_system_ids is_2d $is2D {

    set wdg $m.rf.r$i

   # Add help event for the radiobutton, using text as field name!
    set wfld [DataField::fieldNameSifToGui [$wdg cget -text]]
    bind $wdg <ButtonPress-3> "Widget::genericButton-3 Coordinate $wfld $wdg"

    $wdg configure -anchor w -state normal \
                   -variable Coordinate($fld) \
                   -indicatoron 1 \
                   -command "Coordinate::setCoordinateStatus" \
                   -selectcolor $Info(select_color)

    # Disable wrong dimensions
    if { $Model(GEOMETRY_DIMENSION) == "2D" && !$is_2d ||
         $Model(GEOMETRY_DIMENSION) == "3D" && $is_2d  
       } {
        $wdg configure -state disabled
    }
    
    lappend Coordinate(allWidgets,$fld) $wdg
    bind $wdg <ButtonRelease-1> "Panel::panelDataChanged 1 Coordinate $wdg {%A %K}"
  }

  #-Labels according to different coordinate types
  set label_index [Coordinate::getLabelIndex]
  set labels [lindex $Coordinate(coordinateLabels) $label_index]

  # Coordinate mapping entries
  # ==========================
  #
  set fld COORDINATE_MAPPING
  set m $w.f1.cd.2
  frame $m.lf
  label $m.lf.l  -anchor nw -text "Coordinate mapping:"

  # Add help event for the label
  set wdg $m.lf.l
  bind $wdg <ButtonPress-3> "Widget::genericButton-3 Coordinate $fld $wdg"

  frame $m.ef

  label $m.ef.l1  -text [lindex $labels 0] -width 5
  entry $m.ef.e1  -textvariable Coordinate(COORDINATE_MAPPING_1) -width 2 -font $Info(entryFont)

  label $m.ef.l2  -text [lindex $labels 1] -width 5
  entry $m.ef.e2  -textvariable Coordinate(COORDINATE_MAPPING_2) -width 2 -font $Info(entryFont)

  #-third entry needs special handling according to the dimension!
  label $m.ef.l3  -text [lindex $labels 2] -width 5

  set state normal
  set bg_color white

  if { $Model(SIMULATION_DIMENSION) != "3D" } {
    set state disabled
    set bg_color gray
  }
  entry $m.ef.e3  -textvariable Coordinate(COORDINATE_MAPPING_3) -width 2 \
                  -bg $bg_color -state $state -font $Info(entryFont)


  # Store widgets and set bindings
  set Coordinate(allWidgets,label,$fld) ""
  set Coordinate(allWidgets,$fld) ""

  foreach i { 1 2 3 } {
    set lbl $m.ef.l$i
    set wdg $m.ef.e$i
    
    set var $fld
    append var "_$i"

    lappend Coordinate(allWidgets,label,$fld) $lbl
    lappend Coordinate(allWidgets,$fld) $wdg

    set Coordinate(allWidgets,$var) $wdg

    bind $wdg <KeyRelease> "Panel::panelDataChanged 1 Coordinate $wdg {%A %K}"
    Widget::setEntryBindings non_standard Coordinate $var $wdg
  }

 
  # Set values
  #
  foreach i {0 1 2} n {1 2 3} {

    set var $fld
    append var "_$n"
    
    set Coordinate($var) [lindex $Coordinate($fld) $i]

    # NOTE: This is needed because mapping vars are not "normal"
    set Coordinate($var,err) 0
    set Coordinate($var,mod) 0

    set Coordinate($var,prev) $Coordinate($var)
    set Coordinate($var,err,prev) 0
    set Coordinate($var,mod,prev) 0
  }
  

  #---WIDGET PACKING
  set fpx $Info(framePadX1)
  set fpy $Info(framePadY1)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #-----Geometry & time widgets packing
  #pack $w.f1 -side top  -anchor nw -fill x -padx $fpx -pady $fpy
  pack $w.f1 -side top  -anchor nw -fill x -padx $fpx3 -pady $fpy3
  pack $w.f1.cd -side left  -anchor nw -fill x -padx $fpx -pady $fpy

  #---Coordinates

  # Coordinate systems
  set m $w.f1.cd.1
  pack $m -side top  -anchor nw -fill x -padx $fpx -pady $fpy
  pack $m.lf -side top  
  pack $m.lf.l -side top -anchor nw  
  pack $m.rf -side top  
  foreach i $coordinate_system_ids {
    pack $m.rf.r$i -side top -anchor nw
  }  

  # Mapping
  set m $w.f1.cd.2
  pack $m -side top  -anchor nw -fill x -padx $fpx -pady $fpy
  pack $m.lf -side top  
  pack $m.lf.l -side top -anchor nw  
  pack $m.ef -side top  
  pack $m.ef.l1 $m.ef.e1 -side left -anchor nw  
  pack $m.ef.l2 $m.ef.e2 -side left -anchor nw  
  pack $m.ef.l3 $m.ef.e3 -side left -anchor nw  

  #-----Buttons packing widgets packing
  pack $w.fB -side top  -padx $fpx -pady $fpy
  #-----Apply, Ok and cancel buttons creating and packing

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "Coordinate::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "Coordinate::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply  -text Apply -command Coordinate::panelApply \
                                  -state $ap]

  focus $ok_btn
  set Coordinate(applyButton)  $ap_btn
  set Coordinate(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -padx $fpx 

  #-----Initialization
  #-Nothing so far
 
  
}


proc Coordinate::panelSave { {inform_front 1} } {
  global Info Coordinate Model

  set panel $Coordinate(parameterType)

  #--Store old values
  Panel::backupFields Coordinate

  set sdim1 $Model(SIMULATION_DIMENSION)
  Coordinate::applyCoordinateDimension $Coordinate(COORDINATE_SYSTEM)
  set sdim2 $Model(SIMULATION_DIMENSION)

  # If simulation dimension was changed
  #
  if { $sdim1 != "" && $sdim1 != $sdim2 } {

    StdPanelInit::updateEquationDataAndMasks 
    Equation::setOtherMasks
    StdPanelInit::createObjectTableMasks
    Equation::constructProblemMask

    # Solver variable dofs!
    set modified 0
    Solver::updateVariableDofs modified
    if { $modified } {
      Solver::panelSave
    }

    # Calculator variable dofs!
    set modified 0
    Calculator::updateVariableDofs modified
    if { $modified } {
      Calculator::okProcPre
      StdPanelExec::panelSave Calculator
    }

    Panel::checkMaskedParameters modified
  }

  #--Set vector data size
  set Coordinate(COORDINATE_MAPPING,dataSize) 3

  #--Form parameter data
  set Coordinate(ids) 1
  DataField::formNonStandardParameter Coordinate 1 "Coordinates1"

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  set Model(Solver,needsUpdate) 1

  Panel::panelDataChanged 0 Coordinate 
  Panel::panelDataModified 0 Coordinate 
  Util::cpp_exec coordinatePanelOk
}


proc Coordinate::panelOk {w} {
  global Coordinate
    
  #---No changes
  if { !$Coordinate(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![Coordinate::checkPanelData] } {
    return
  }

  #---Ok 
  Coordinate::panelSave
  Panel::cancel $w
}


proc Coordinate::panelApply {} {
  global Info Coordinate

  #---No changes
  if { !$Coordinate(dataChanged) } {
    return
  }

  if { ![Coordinate::checkPanelData] } {
    return
  }
   
  Coordinate::panelSave

  # This is needed to make mapping fields in
  # correct status in the panel
  Coordinate::setCoordinateStatus
}


proc Coordinate::getLabelIndex {} {
  global Coordinate

  switch $Coordinate(COORDINATE_SYSTEM) {
    "Cartesian 2D"        -
    "Cartesian 3D"        {return 0}
    "Axi Symmetric"       -
    "Cylindric Symmetric" -
    "Cylindrical"         {return 1}
    "Polar 2D"            -
    "Polar 3D"            {return 2}
  }

  #-Default if not found!
  return 0
}


proc Coordinate::updateLabels {} {
  global Coordinate

  set index [Coordinate::getLabelIndex]
  set Coordinate(coordinateLabelIndex) $index
  set labels [lindex $Coordinate(coordinateLabels) $index]

  foreach wid $Coordinate(allWidgets,label,COORDINATE_MAPPING) lbl $labels {
    set text $lbl
    append text ":"
    $wid configure -text $text
  }
}

proc Coordinate::panelCancel {w} {
  global Coordinate

  if { ![Panel::verifyCancel Coordinate] } {
    return
  }

  #---Reset into old values
  Panel::restoreFields Coordinate

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc Coordinate::checkPanelData {} {
  global Info Coordinate Model

  # Coordinate mapping
  # ==================
  # we must have a list of {1 2 3} (in any order)!
  set mappings ""
  set is_ok 1

  # Note: We check from end to start to display
  # the possible same-value error for the smaller
  # number, because this is sensible for 2D, where
  # the 3rd coordinate should always be 3 when it is
  # not editable
  #
  foreach i {3 2 1} {
    set var COORDINATE_MAPPING_$i
    set val [string trim $Coordinate($var)]
    set Coordinate($var) $val

    # NOTE: This is needed because mapping vars are not "normal"
    set Coordinate($var,err) 0
    set Coordinate($var,mod) 0
    set is_ok 1

    # Illegal value or same values
    #
    if { ($val < 1 || $val > 3) ||
         -1 != [lsearch $mappings $val]    
      } {
      set Coordinate($var,err) 1
      set is_ok 0
    }

    Widget::setFieldStatus Coordinate $var
  
    if {$is_ok} {
      lappend mappings $val
    }
  }

  #--If error
  if { 3 != [llength $mappings] } {

    set Info(messageIcon) error
    Message::dispOkMessage {"Incorrect coordinate mapping!"} \
                  "$Info(FRONT_NAME) error message" \
                  $Coordinate(winName) 

    return 0
  }

  set Coordinate(COORDINATE_MAPPING) "" 

  #--If mapping data ok
  foreach i {1 2 3} {
 
    set var COORDINATE_MAPPING_$i

    # NOTE: This is needed because mapping vars are not "normal"
    set Coordinate($var,prev) $Coordinate($var)
    set Coordinate($var,err,prev) $Coordinate($var,err)
    set Coordinate($var,mod,prev) $Coordinate($var,mod)

    lappend Coordinate(COORDINATE_MAPPING) $Coordinate($var)
  }

  # Ok
  return 1
}


# Set coordinate panel widgets state according
# to the model geometry dimension
#
proc Coordinate::setCoordinateStatus {} {
  global Coordinate

  Coordinate::updateLabels

  switch $Coordinate(COORDINATE_SYSTEM) {
    "Cartesian 2D"        -
    "Polar 2D"            -
    "Axi Symmetric"       { set activate_3 0 }
    "Cylindric Symmetric" -
    "Cartesian 3D"        -
    "Cylindrical"         -
    "Polar 3D"            { set activate_3 1 }
  }

  set wdg [lindex $Coordinate(allWidgets,COORDINATE_MAPPING) 2]

  if { $activate_3 } {
    $wdg configure -state normal -bg white

  } else {
    $wdg configure -state disabled -bg gray

    # Note: We have to force a legal mapping if
    # dimension is changed, in practice this means
    # 2D and CylindricSymmetric changes
    set Coordinate(COORDINATE_MAPPING_3) 3

    if { $Coordinate(COORDINATE_MAPPING_1) <= \
         $Coordinate(COORDINATE_MAPPING_2)
       } {
      set Coordinate(COORDINATE_MAPPING_1) 1
      set Coordinate(COORDINATE_MAPPING_2) 2

    } else { 
      set Coordinate(COORDINATE_MAPPING_1) 2
      set Coordinate(COORDINATE_MAPPING_2) 1
    }

  }
}


# Updates system according to the simulation dimension value
#
proc Coordinate::applyCoordinateDimension {coordinate_system} {
  global Coordinate Model

  switch $coordinate_system {
    "Cartesian 2D"        -
    "Polar 2D"            -
    "Axi Symmetric"       { set Model(SIMULATION_DIMENSION) "2D" }
    "Cylindric Symmetric" -
    "Cartesian 3D"        -
    "Cylindrical"         -
    "Polar 3D"            { set Model(SIMULATION_DIMENSION) "3D" }
  }

  set Coordinate(coordinateLabelIndex) [Coordinate::getLabelIndex]
}


# Set default coordinate system
#
proc Coordinate::setDefaultCoordinateSystem {} {
  global Coordinate Model

  if { $Model(GEOMETRY_DIMENSION) == "2D" } {
    set Coordinate(COORDINATE_SYSTEM) "Cartesian 2D"

  } else {
    set Coordinate(COORDINATE_SYSTEM) "Cartesian 3D"
  }
    
  # Store new value in the parameter
  set id $Coordinate(parameterId)
  DataField::setFieldValue Coordinate $id COORDINATE_SYSTEM $Coordinate(COORDINATE_SYSTEM)

  execNsProc Coordinate panelSave 0
}



# end ecif_tk_coordinateDefPanel.tcl
# ********************
