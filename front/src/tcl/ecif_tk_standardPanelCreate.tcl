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
#Module:    ecif_tk_panelCreate.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A procedure for creating a "standard" panel.
#           These are panel with the common structure where we 
#           bodies, [elements], paremeters and data values for
#           the parameters. Examples are initialConditons, 
#           boundaryCondtions, materials etc.
#
#************************************************************************
 
#------Standard panel construction proc------
#

###########################
### Screen construction ###
### and initialization  ###
###########################

# Screen is divided into four main parts:
# 1: Objects listbox (fObjects),
#    Boundaries listbox (fBoundaries)
#    Parameter listbox area (fParam)
# 2: Apply,Remove; New,Back,Add,Updqate,Delete params handling buttons (fButtons1)
# 3: Problem Buttons are (not all panels have this)
# 4: Data values (fValues)
# 5: Ok and cancel etc. buttons (fButtons2)


proc StdPanelCreate::openPanel {globArray} {
  global $globArray Info Common Equation SolverSystem
  upvar #0 $globArray theArray

  set panel $theArray(parameterType)

  set wgeom $theArray(winGeometry)

  set w $theArray(winName)
  set id [winfo atom $w]
  set theArray(winId) $id
  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow $globArray $id $theArray(winTitle) $wgeom] } {
    return
  }  

  set theArray(entryInfo) ""

  # Variables and Problem indices
  # =============================
  set theArray(targetMask) ""
  Panel::initFields $globArray

  #-Equation panel uses its own problem variables
  if { $globArray == $Common(panelArray,EQ) } {
    set theArray(equationIndices) $Equation(ownIndices)

  #-All other panels use equation attachment based variables!
  } else {
    set theArray(equationIndices) $Equation(activeIndices)
  }

  if { $theArray(equationIndices) != "" } {
    set theArray(currentEquationIndices) $theArray(equationIndices)
  } else {
    set theArray(currentEquationIndices) $Equation(indices)
  }

  # If no problems, we display them all, because
  # the panel is inactive anyway!
  #
  if { $theArray(equationIndices) == "" } {
    set theArray(equationIndices) $Equation(indices)
  }
    
  # Current problem index
  # =====================
  #-Check that last current problem index is still active
  if { [info exist theArray(lastEquationIndex)]  &&
       -1 != [lsearch $theArray(equationIndices) $theArray(lastEquationIndex)]
     } {
    set theArray(equationIndex) $theArray(lastEquationIndex)
    
  #-Otherwise take the first active
  } else {
    set theArray(equationIndex) [lindex $theArray(equationIndices) 0]
  }

  set theArray(lastEquationIndex) $Info(NO_INDEX)


  # Panel variables and labels
  # ==========================
  if { $theArray(panelsByEquation) } {
    set problem_label [Equation::getFieldPropertyByIndex $theArray(equationIndex) EquationLabel]
    append problem_label " "
    Panel::selectDisplayFields $globArray
  
  } else {
    set problem_label ""
  }

  set theArray(objectName) ""
  set theArray(parameterName) ""
  set theArray(currentEntryField) ""
  set theArray(valuesHdr) $problem_label$theArray(valuesHdrBase)

  #--Find panel size parameters (an array of nofFrames)
  # (per problem if needed)
  # =========================
  if { !$theArray(panelSizeKnown) } {
    StdPanelCreate::calcPanelSize $globArray
    set theArray(panelSizeKnown) 1
  }

  #--Window properties
  set this $w
  toplevel $w
  focus $w 
  wm title $w $theArray(winTitle)
  wm geometry $w $wgeom
  
  #---MAIN for the panel
  StdPanelCreate::initPanelData $globArray
  StdPanelCreate::createPanelWidgets $globArray $w
  Widget::setPanelBindings $globArray $w
  StdPanelInit::createWidgetBindings $globArray
  StdPanelInit::initPanelData $globArray
  StdPanelCreate::packPanel $globArray $w
  StdPanelCreate::constructPanelValuesArea $globArray  $theArray(equationIndex)

  #--Prev/Next equation buttons
  if { $theArray(panelsByEquation) } {
    StdPanelExec::setEqMoveButtonsState $globArray
  }

  #--Prev/Next page buttons
  if { $theArray(panelsByEquation) } {
		set theArray(maxNofPanelPages) $theArray($theArray(equationIndex),nofPanelPages)
  } else {
		set theArray(maxNofPanelPages) $theArray(nofPanelPages)
	}

	if { $theArray(panelPageNbr) > 1 } {
		$theArray(prevPageButton) configure -state normal
	}

	if { $theArray(panelPageNbr) < $theArray(maxNofPanelPages) } {
		$theArray(nextPageButton) configure -state normal
	}

}
#-end openPanel proc


# Calc panel max size (nof valueareas and value area size)
#
proc StdPanelCreate::calcPanelSize {globArray} {
  upvar #0 $globArray theArray

  StdPanelCreate::setNofValuesAreaFrames $globArray

  if { $theArray(panelsByEquation) } {
    StdPanelCreate::calcValuesAreaSize $globArray
  }
}


# Call this for any array-specific processing before
# creating panel widgets
#
proc StdPanelCreate::initPanelData {globArray} {
  global Equation
  upvar #0 $globArray theArray

  # Currently nothing!
  # =================
  return

  # A trial to handle HeatCapacity field, which is in Heat and in ComprFlow
  if { $globArray == "Material" } {

    set eidx [DataField::getFieldProperty Equation HEAT_EQUATION EquationIndex]

    if { 1 == [lsearch $theArray(equationIndices) $eidx] } {
      DataField::setFieldProperty $globArray HEAT_CAPACITY Mask "(¤HT)"
      DataField::setFieldProperty $globArray HEAT_CAPACITY ActivityOnlyByMask 1

    } else {
      DataField::setFieldProperty $globArray HEAT_CAPACITY Mask "(¤FL)"
      DataField::setFieldProperty $globArray HEAT_CAPACITY ActivityOnlyByMask 0
    }
  }
}

proc StdPanelCreate::createPanelWidgets {globArray top_level} {
  global Info Common $globArray
  upvar #0 $globArray theArray

  set w $top_level

  #----WIDGET DEFINITION
  
  #---Outer frames
  frame $w.f1 ;# List boxes
  frame $w.f2 ;# Update etc. buttons
  frame $w.f3 ;# Equation selection area
  frame $w.f4 ;# Values area
  frame $w.f5 ;# Ok etc. buttons

  #-Object list-boxes frame 
  frame $w.f1.fObjectBox
  frame $w.f1.fObjectBox.fObjects 
  frame $w.f1.fBoundaryBox
  frame $w.f1.fBoundaryBox.fBoundaries
  frame $w.f1.fParamBox
  frame $w.f1.fParamBox.fParams

  #-Buttons-1 frame
  frame $w.f2.fButtons1 

  #-Problem menu frame
  frame $w.f3.fEquationMenu

  #-Values frame
  set wdg [frame $w.f4.fValues -relief groove -bd 2]
  set theArray(valuesFrame) $wdg
  set theArray(valuesParentFrame) $w.f4

  # Fixed sized for values area if multipanel!
  if { $theArray(panelsByEquation) } {
    $wdg configure -width $theArray(valuesFrame,width)
    $wdg configure -height $theArray(valuesFrame,height)
    pack propagate $wdg 0
  }

  #-Buttons-2 frame
  frame $w.f5.fButtons2 

  #---Frame internals
  #-Object listbox
  StdPanelCreate::createObjectBoxArea $w.f1.fObjectBox $globArray

  #-Boundary listbox
  StdPanelCreate::createBoundaryBoxArea $w.f1.fBoundaryBox $globArray

  #-Parameters
  StdPanelCreate::createParamBoxArea $w.f1.fParamBox $globArray

  #-Buttons-1
  StdPanelCreate::createButtons1Area $w.f2.fButtons1 $globArray

  #-Equation menu
  StdPanelCreate::createEquationMenuArea $w.f3.fEquationMenu $globArray

  #-Values 
  StdPanelCreate::createValuesArea $w.f4.fValues $globArray

  #-Buttons-2
  StdPanelCreate::createButtons2Area $w.f5.fButtons2 $globArray $w
  
  #-end WIDGET DEFINITION
}


proc StdPanelCreate::packPanel { globArray w } {
  global Info
  upvar #0 $globArray theArray

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
  if { $theArray(panelsByEquation) } {
    set frames [list $w.f5 $w.f4 $w.f3 $w.f2 $w.f1 ]
  } else {
    set frames [list $w.f5 $w.f4 $w.f2 $w.f1 ]
  }
  foreach fr $frames {
    pack $fr -side bottom -expand $expnd -fill $fill -padx $fpx1 -pady $fpy1
  }

  #---Individual frames
  #-Object list-box frame
  pack $w.f1.fObjectBox -side left -padx 2m -expand $expnd -fill $fill -anchor $anch \
                        -padx $fpx0 -pady $fpy0

  #-Boundary list-box frame
  pack $w.f1.fBoundaryBox -side left -padx 2m -expand $expnd -fill $fill -anchor $anch \
                          -padx $fpx0 -pady $fpy0

  #-Params list-box frame
  pack $w.f1.fParamBox -side left -padx 2m -expand $expnd -fill $fill -anchor $anch \
                       -padx $fpx0 -pady $fpy0

  #-Buttons-1 frame (Attach etc, Body etc. names)
  pack $w.f2.fButtons1 -expand $expnd -fill $fill \
                       -padx $fpx0 -pady $fpy0

  #-Problem buttons frame
  pack $w.f3.fEquationMenu -padx $fpx0 -pady $fpy3

  #-Values frames
  pack $w.f4.fValues -expand $expnd -fill $fill -anchor $anch \
                     -padx $fpx1 -pady $fpy1

  #-Buttons-2 frame
  pack $w.f5.fButtons2 -expand $expnd  -padx $fpx0 -pady $fpy0
  

  #---Frame internals
  #-Object list-box
  StdPanelCreate::packObjectBoxArea $w.f1.fObjectBox $globArray

  #-Element list-box
  StdPanelCreate::packBoundaryBoxArea $w.f1.fBoundaryBox $globArray

  #-Params list-box
  StdPanelCreate::packParamBoxArea $w.f1.fParamBox $globArray

  #-Buttons-1
  StdPanelCreate::packButtons1Area $w.f2.fButtons1 $globArray

  #-Equation menu
  StdPanelCreate::packEquationMenuArea $w.f3.fEquationMenu $globArray

  #-Values
  StdPanelCreate::packValuesArea $w.f4.fValues $globArray

  #-Buttons-2
  StdPanelCreate::packButtons2Area $w.f5.fButtons2 $globArray

  #-end SCREEN PACKING
}


##########################
#### Helper functions ####
##########################


# Objects listbox area
# =======================

#--Creating
proc StdPanelCreate::createObjectBoxArea {frame globArray {wh ""} {ww ""}} {
  upvar #0 $globArray theArray
  global Common Info

  if { $wh == "" } {
    set wh $Info(standardLbHeight)
  }

  if { $ww == "" } {
    if {$theArray(hasBoundaries)} {
      set ww $Info(standardObjectLbWidth3)
    } else {
      set ww $Info(standardObjectLbWidth2)
    }
  }

  #---Objects area internals (a listbox with v-scrollbar)
  set m $frame.fObjects

  label $m.l -text $theArray(objectHdr)

  if { [info exists theArray(canMultiSelectObjects)] &&
       $theArray(canMultiSelectObjects)
     } {
    set select_mode extended
    #set select_mode multiple
  } else {
    set select_mode browse
  }


  listbox $m.objects -relief sunken -selectmode $select_mode -exportselection 0 \
    -yscrollcommand [list $m.sy set] \
    -height $wh -width $ww
    
  scrollbar $m.sy -orient vertical -command [list $m.objects yview]
  
  set theArray(objectLB) $m.objects
}

#--Packing
proc StdPanelCreate::packObjectBoxArea {frame globArray} {
  upvar #0 $globArray theArray
  global Common Info

  set expnd 1
  set fill both
  pack $frame.fObjects -side top -expand $expnd -fill $fill

  #--Objects frame area internals
  set mo $frame.fObjects
  pack $mo.l -side top -anchor nw -expand $expnd -fill $fill  
  pack $mo.objects -side left -anchor n -fill $fill -expand $expnd
  pack $mo.sy -side left -anchor n -fill y 
}


# Boundary listbox area
# =====================

#--Creating
proc StdPanelCreate::createBoundaryBoxArea {frame globArray {wh ""} {ww ""}} {
  upvar #0 $globArray theArray
  global Common Info

  #-First we check if boundary listbox belongs to the panel!
  if {!$theArray(hasBoundaries)} {
    return
  }

  if { $wh == "" } {
    set wh $Info(standardLbHeight)
  }

  if { $ww == "" } {
    set ww $Info(standardBoundaryLbWidth)
  }

  #---Boundary area internals (a listbox with v-scrollbar)
  if { [info exist theArray(canMultiSelectBoundaries)] &&
       $theArray(canMultiSelectBoundaries)
     } {
    set select_mode extended
    #set select_mode multiple
  } else {
    set select_mode browse
  }

  set m $frame.fBoundaries

  label $m.l -text $theArray(boundaryHdr)

  listbox $m.boundaries -relief sunken -selectmode $select_mode -exportselection 0 \
    -yscrollcommand [list $m.sy set]  \
    -height $wh -width $ww

  scrollbar $m.sy -orient vertical -command [list $m.boundaries yview]

  set theArray(boundaryLB) $m.boundaries
}


#--Packing
proc StdPanelCreate::packBoundaryBoxArea {frame globArray} {
  upvar #0 $globArray theArray
  global Common Info

  #-First we check if boundary listbox belongs to the panel!
  if {!$theArray(hasBoundaries)} {
    return
  }

  set expnd 1
  set fill both
  pack $frame.fBoundaries -side top -expand $expnd -fill $fill

  #--Boundary frame area internals
  set m $frame.fBoundaries
  pack $m.l -side top -anchor nw -expand $expnd  -fill $fill
  pack $m.boundaries -side left -anchor n -fill $fill -expand $expnd
  pack $m.sy -side left -anchor n -fill y 
}


# Parameters listbox area
# =======================
 
#--Creating
proc StdPanelCreate::createParamBoxArea {frame globArray {wh ""} {ww ""} } {
  upvar #0 $globArray theArray
  global Common Info

  if { $wh == "" } {
    set wh $Info(standardLbHeight)
  }

  if { $ww == "" } {
    set ww $Info(standardParameterLbWidth)
  }

  set f $frame.fParams

  label $f.l -text $theArray(paramsHdr)

  listbox $f.param -relief sunken -selectmode browse \
    -yscrollcommand [list $f.sy set] -exportselection 0 \
    -height $wh -width $ww 

  scrollbar $f.sy -orient vertical -command [list $f.param yview]

  set theArray(parameterLB) $f.param
}

#--Packing
proc StdPanelCreate::packParamBoxArea {frame globArray} {
  upvar #0 $globArray theArray
  global Common Info

  set expnd 1
  set fill both
  pack $frame.fParams -side top -expand $expnd -fill $fill
  set mc $frame.fParams
  pack $mc.l -side top -anchor nw -expand $expnd -fill $fill
  pack $mc.param -side left -fill x -expand $expnd
  pack $mc.sy -side right -fill y 
}


# Button1 area (Attach etc. Body, Boundary and Parameter names)
# ==========================

#--Creating
proc StdPanelCreate::createButtons1Area {frame globArray} {
  global Common Info ModelFlags
  upvar #0 $globArray theArray

  frame $frame.sf1
  frame $frame.sf2

  # Attach and Detach buttons
  button $frame.sf1.attach -text Attach -command "StdPanelExec::attachParameter $globArray"
  button $frame.sf1.detach -text Detach -command "StdPanelExec::detachParameter $globArray"

  # Empty separator frame
  frame $frame.sf1.separator 

  # Keep, Back, Add, Update and Delete buttons
  button $frame.sf1.keep -text Keep -command "StdPanelExec::keepParameter $globArray"
  button $frame.sf1.back -text Back -state disabled -command "StdPanelExec::backParameter $globArray"
  button $frame.sf1.add -text Add -command "StdPanelExec::addParameter $globArray"
  button $frame.sf1.update -text Update -command "StdPanelExec::updateParameter $globArray"
  button $frame.sf1.delete -text Delete -command "StdPanelExec::deleteParameter $globArray"

  set wdg $frame.sf1.keep
  set theArray(panelNewButton) $wdg

  set wdg $frame.sf1.back
  set theArray(panelBackButton) $wdg

  set wdg $frame.sf1.attach;  bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  set theArray(panelAttachButton) $wdg

  set wdg $frame.sf1.detach;  bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  set theArray(panelDetachButton) $wdg

  set wdg $frame.sf1.add;     bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  set theArray(panelAddButton) $wdg

  set wdg $frame.sf1.update;  bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  $wdg configure -state disabled
  set theArray(panelUpdateButton) $wdg

  set wdg $frame.sf1.delete;  bind $wdg <ButtonPress-1> "+Panel::panelDataChanged 1 $globArray $wdg {%A %K}"

  if { ![info exists theArray(ids)] ||
       0 == [llength $theArray(ids)]
     } {
    $wdg configure -state disabled
  }

  set theArray(panelDeleteButton) $wdg

  # These updates make data "unmodified"
  set wdg $frame.sf1.add;     bind $wdg <ButtonRelease-1> "+Panel::panelDataModified 0 $globArray $wdg {%A %K}"
  set wdg $frame.sf1.update;  bind $wdg <ButtonRelease-1> "+Panel::panelDataModified 0 $globArray $wdg {%A %K}"

  #-Store buttons to panel specific array
  set m $frame.sf1
  set theArray(buttons1) [list $m.attach $m.detach $m.add $m.update $m.delete]

  # Body name
  # =========

  # Body name variable
  set bn_var $globArray
  append bn_var (objectName)

  # Object (body) name is normally non-editable
  label $frame.sf2.oname_label -text "Body: "
  entry $frame.sf2.oname -textvariable $bn_var -font $Info(entryFont) \
                         -state disabled -width 15 -background $Info(nonActiveBg)
  frame $frame.sf2.separator1
  set theArray(bodyNameWidget) $frame.sf2.oname

  # Boundary (element) name is normally non-editable
  # =======================
  # Boundary name variable
  set el_var $globArray
  append el_var (boundaryName)

  if { $theArray(hasBoundaries) } {
    label $frame.sf2.bname_label -text "Boundary: "
    entry $frame.sf2.bname -textvariable $el_var -font $Info(entryFont) \
                           -state disabled  -width 15 -background $Info(nonActiveBg)
    frame $frame.sf2.separator2
    set theArray(boundaryNameWidget) $frame.sf2.bname
  }

  # Parameter name is editable
  # ==========================
  label $frame.sf2.pname_label -text "Name: "

  if { $theArray(hasBoundaries) } {
    set wid 18
  } else {
    set wid 25
  }

  # Parameter name variable
  set pn_var $globArray
  append pn_var (parameterName)

  entry $frame.sf2.pname -textvariable $pn_var -font $Info(entryFont)  -width $wid

  set wdg $frame.sf2.pname
  set theArray(parameterNameWidget) $wdg
  set theArray(parameterName,err) 0
  set theArray(parameterName,mod) 0

  # Set bindings
  set pn "parameterName"
  bind $wdg <KeyRelease> "+StdPanelExec::updateCurrentParameterName $globArray $wdg"
  bind $wdg <KeyRelease> "+Panel::panelDataModified 1 $globArray $wdg"
  bind $wdg <KeyRelease> "+Widget::entryKeyRelease $globArray $pn $wdg {%A %K}"
  bind $wdg <FocusIn> "+Panel::setProcAndTableButtonStates $globArray disabled disabled"
  bind $wdg <KeyPress-Escape> "+Widget::entryKeyPress-Escape $globArray $pn $wdg {%A %K}"
}


#--Packing
proc StdPanelCreate::packButtons1Area {frame globArray} {
  global Common Info
  upvar #0 $globArray theArray

  set expnd 0
  set fill both
  set px 5; set py 4

  # Attach and Detach buttons
  pack $frame.sf1.attach  -side left -pady $py -padx $px -anchor w
  pack $frame.sf1.detach -side left -pady $py -padx $px -anchor w

  # Empty separator frame
  pack $frame.sf1.separator -side left -pady $py -padx $px -anchor w

  # New, Back, Add, Upadate and Delete buttons
  pack $frame.sf1.delete -side right -pady $py -padx $px -anchor e
  pack $frame.sf1.update -side right -pady $py -padx $px -anchor e
  pack $frame.sf1.add  -side right -pady $py -padx $px -anchor e
  pack $frame.sf1.back  -side right -pady $py -padx $px -anchor e
  pack $frame.sf1.keep  -side right -pady $py -padx $px -anchor e

  #-Object (body) name
  pack $frame.sf2.oname_label  -side left -pady 1  -padx 3 -anchor w 
  pack $frame.sf2.oname  -side left -pady 1  -padx 3 -anchor w 
  pack $frame.sf2.separator1  -side left -pady 1  -padx 1 -anchor w 

  #-Boundary (element) name
  if {$theArray(hasBoundaries)} {
    pack $frame.sf2.bname_label  -side left -pady 1  -padx 3 -anchor w 
    pack $frame.sf2.bname  -side left -pady 1  -padx 3 -anchor w 
    pack $frame.sf2.separator2  -side left -pady 1  -padx 1 -anchor w 
  }

  # parameter name
  pack $frame.sf2.pname -side right -pady 1 -padx 2 -anchor w
  pack $frame.sf2.pname_label -side right -pady 1 -padx 2 -anchor w 
  
  pack $frame.sf1 $frame.sf2 -side top -fill $fill -expand $expnd -padx 3
  #pack $frame.sf1 $frame.sf2 -side top
}


# Equation selection area
# =======================

#--Creating
proc StdPanelCreate::createEquationMenuArea {frame globArray} {
  global Info
  upvar #0 $globArray theArray

  if { !$theArray(panelsByEquation) } {
    return
  }

  set m [frame $frame.frb]

  set theArray(allWidgets,equationMenuFrame) $m
  
  #--Label
  label $m.l -text "Parameters for equation: " -width 20 -justify left
  pack $m.l -side left -anchor w -fill x

  #--Create OptionMenu for Equations
  StdPanelCreate::createEquationMenu $globArray $theArray(equationIndices)

  #--Next/Previous equation buttons
  set prev_btn [button $m.prev -text "<<" -command "StdPanelExec::panelEqMove $globArray -1 " \
                                     -state disabled -disabledforeground $Info(nonActiveFg) ]

  set next_btn [button $m.next -text ">>" -command "StdPanelExec::panelEqMove $globArray 1" \
                                      -state disabled -disabledforeground $Info(nonActiveFg) ]

  set theArray(prevEqButton) $prev_btn
  set theArray(nextEqButton) $next_btn

  pack $next_btn -side right -anchor w -fill x 
  pack $prev_btn -side right -anchor w -fill x -padx 10
}
 


proc StdPanelCreate::createEquationMenu { globArray indices} {
  global Equation
  upvar #0 $globArray theArray

  set m $theArray(allWidgets,equationMenuFrame)

  set eqlabels ""

  foreach idx $indices {
    lappend eqlabels [Equation::getFieldPropertyByIndex $idx EquationLabel]
  }

  #-No equations for the body
  if { $indices == "" } {
    set theArray(equationLabel) ""
  
  #-Pick last visited (if applicable) or first active
  } else {

    if { -1 != [lsearch $indices $theArray(lastEquationIndex)] } { 
      set index $theArray(lastEquationIndex)
    } else {
      set index  [lindex $theArray(equationIndices) 0]
    }
    
    set theArray(equationLabel) [Equation::getFieldPropertyByIndex $index EquationLabel]
  }
  
  set el_var $globArray
  append el_var "(equationLabel)"

  set wdg $m.om

  if { [winfo exists $wdg] } {
    destroy $wdg
  }

  set om [Panel::createOptionMenuWidget $m.om $el_var $eqlabels]
  $wdg configure -indicatoron 0 -width 30 
  #$wdg configure -state disabled
  $wdg configure -disabledforeground black
  set cmd "PanelCheck::checkOptionMenu $om $wdg $globArray equationLabel "
  append cmd "\"StdPanelExec::equationMenuSelectionProc $globArray \" {%s}"
  #bind $om <<MenuSelect>> "PanelCheck::checkOptionMenu $om $wdg $globArray equationLabel \"\" {%s}"
  bind $om <<MenuSelect>> $cmd
  set theArray(equationWidget) $wdg

  pack $wdg -side left -anchor w -fill x
}


#--Packing
proc StdPanelCreate::packEquationMenuArea {frame globArray} {
  global Common Info Equation
  upvar #0 $globArray theArray

  if { !$theArray(panelsByEquation) } {
    return
  }

  set expnd 0
  set fill both

  pack $frame.frb -side top -fill $fill -expand $expnd
  set m $frame.frb

  pack $m
}


# Values area
# ===========

#--Creating
#
proc StdPanelCreate::createValuesArea {frame globArray} {
  global Common Equation Info 
  upvar #0 $globArray theArray

  if { $theArray(panelsByEquation) &&
       [info exists theArray(equationIndex)]
     } {
    set eq_index $theArray(equationIndex)

  } else {
    set eq_index ""
  }

  # NOTE: This is a harmful init !!!???
  Panel::initValuesAreaFields $globArray "" 0 $eq_index

  #-Select display fields
  Panel::selectDisplayFields $globArray

  set m $frame
	set mp $theArray(panelPageNbr)

  #--Create all values sub-frames

  #-Parameter file field frame (INCLUDE field, can be left free of entries!)
  if { ![info exists m.valuesSubFrame0] &&
       ![winfo exists $m.valuesSubFrame0]
     } {
    frame $m.valuesSubFrame0
  }

  #--Data field sub-frames
  if { $eq_index != "" } {
    set nof_sub_frms $theArray($eq_index,$mp,nofSubPanels)
  } else {
    set nof_sub_frms $theArray($mp,nofSubPanels)
  }

  for {set i 1} {$i <= $nof_sub_frms} {incr i} {

    if { ![info exist m.valuesSubFrame$i] &&
         ![winfo exist $m.valuesSubFrame$i]
       } {
      frame $m.valuesSubFrame$i -relief groove -bd 2
    }
  }
  
  #--Common values frame header
  if { ![info exists m.l] &&
       ![winfo exist  $m.l]
     } {
    label $m.l -text $theArray(valuesHdr)
  } else {
    $m.l configure -text $theArray(valuesHdr)
  }

  #--Create fields into proper sub-frame
  set unit_sp "  "
  set fill_sp "          "

  set panel_mask [Panel::getProblemMask $globArray]

  foreach fld $theArray(allFields) {

    #-Only displayed fields!!!
    if { !$theArray($fld,dsp) } {
      continue
    }

    #-Sub-frame for the field
    #
    if { $theArray(panelsByEquation) } {
      set sp_index [DataField::getSubPanel $globArray $fld $eq_index]

    } else {
      set sp_index [DataField::getSubPanel $globArray $fld]
    }
		
  	#-Select only fields which belong to the current Panel Page
    #
		# NOTE: Sub panel-0 fields (Include, Clear, Edit etc. ) are displayed in all
		# main panels!
		#
		if { $sp_index > 0 } {
			set mp_indices [DataField::getPanelPage $globArray $fld $eq_index]

			if { -1 == [lsearch $mp_indices $mp] } {
				continue
			}
		}

    set sub_frame valuesSubFrame$sp_index
    set ff $m.$sub_frame

    #-Create each field into its own frame in the proper subframe !!
    # NOTE: Create  only if the subframe widget exists!!
    #
    if { [info exists ff] && 
         [winfo exists $ff] &&
         ![info exists ff.f$fld] && 
         ![winfo exists $ff.f$fld]
       } {
      set fld_frame [frame $ff.f$fld]
      StdPanelCreate::createField $globArray "" $fld $fld_frame $panel_mask $eq_index
    }

    if { [Panel::isSimpleDataField $globArray $fld] } {
      Widget::setFieldStatus $globArray $fld
    }
  }
}


#--Packing
proc StdPanelCreate::packValuesArea {frame globArray} {
  global Common Info
  upvar #0 $globArray theArray

  set f $frame
	set mp $theArray(panelPageNbr)

  #--Parameter file etc. field frame
  pack $f.valuesSubFrame0 -side top -anchor nw -expand 0 -fill both  

  #--Pack common data field label
  pack $f.l -side top -anchor n -pady 5

  #--Data field sub-frames
  if { $theArray(panelsByEquation) } {
    set nof_sub_frms $theArray($theArray(equationIndex),$mp,nofSubPanels)

  } else {
    set nof_sub_frms $theArray($mp,nofSubPanels)
  }

  for {set i 1} {$i <= $nof_sub_frms} {incr i} {
    pack $f.valuesSubFrame$i -side left -anchor nw -expand 0 -fill both
  }

  # We are ready!
  # =============
  set theArray(valuesAreaPacked) 1

  # Comment this for testing!
  return
  
  #========================================#
  # NOTE: The following is for testing!!!! #
  # To se how fields are packed when       # 
  # created turn packing off in the        #
  # createField proc and use this for pack!#
  #========================================#

  #--Pack each field into proper sub_frame
  foreach fld $theArray(allFields) {

    if { !$theArray($fld,prb) } {
      continue
    }
#print "packing loop for panel = $globArray and field = $fld "
    set index [DataField::getSubPanel $globArray $fld $theArray(equationIndex)]
    set sub_frame valuesSubFrame$index
    set ff $f.$sub_frame
    #StdPanelCreate::packField $globArray $fld $ff.f$fld 1m 3m 
    StdPanelCreate::packField $globArray $fld $ff.f$fld 1 3 
  }
}


proc StdPanelCreate::packField {globArray fld frame fieldSep procedureSep} {
  global Common
  upvar #0 $globArray theArray

  set f $frame

  set isGroup [Panel::isGroupField $globArray $fld]

  if {$isGroup} {
    set panel $theArray(parameterType)
    set GroupFields $Common($panel,$fld,group)
  } else {
    set GroupFields {""}
  }

  pack $f -side top -pady $fieldSep -fill x -anchor n

  #-The following packing order is essential, if we want to
  # pack procedure checkbuttons to the right of the entry fields !!!***!!!
  pack $f.fc -side right -pady $fieldSep -fill x -anchor e
  pack $f.fe -side right -pady $fieldSep -fill x -anchor w
  pack $f.fl -side left -pady $fieldSep -fill x -anchor w

  set ornt [DataField::getFieldProperty $globArray $fld WidgetPacking]
  if {$ornt == "Horizontal"} {
    set wornt left 
  } elseif {$ornt == "Vertical" } {
    set wornt top
  }

  #-Pack individual group elements (x,y,z)
  foreach GF $GroupFields {
    if {$isGroup} {
      set FV $GF
    # A single entry (we don't have a group!)
    } else {
      set FV $fld
    }

    #MSG "packing field:  $fld   GF: $GF   FV: $FV"
    #MSG "Hi, look I'm being built!"

    pack $f.fl.l -side $wornt
    pack $f.fe.e -side $wornt
  }

  #-Pack procedure check-button if defined
  if {1 == [DataField::getFieldProperty $globArray $fld Procedurable]} {
    pack $f.fc.c -side $wornt -padx $procedureSep 
  }
}


#--Destroying Values Area fields
#
proc StdPanelCreate::destroyValuesArea {frame globArray} {
  global $globArray Common Info
  upvar #0 $globArray theArray

   set m $frame
	 set mp $theArray(panelPageNbr)

  #--Common label
  #destroy $m.l

  #--Destroy all values sub-frames

  #-We do not destroy this area, because no changes!
  #destroy $m.valuesSubFrame0

  #--Data field sub-frames
  if { $theArray(panelsByEquation) } {
    set nof_sub_frms $theArray($theArray(equationIndex),$mp,nofSubPanels)
  } else {
    set nof_sub_frms $theArray($mp,nofSubPanels)
  }

  if { $globArray == $Common(panelArray,EQ) } {
    set fstart 2
  } else {
    set fstart 1
  }

  for {set i $fstart} {$i <= 10} {incr i} {
    
    if { [winfo exists $m.valuesSubFrame$i] } {
      destroy $m.valuesSubFrame$i
    }
  }
}


# Button2 area (Ok etc.)
# ======================

#--Creating
proc StdPanelCreate::createButtons2Area {frame globArray this} {
  global Common Info Model 
  upvar #0 $globArray theArray

  frame $frame.f1
  frame $frame.f2
  
	# Entry info in the "status" bar
	#
  label $frame.f1.l -text " " -width 20 -height 3 -justify left
  set theArray(entryInfoLabel) $frame.f1.l

	# Buttons
	#
  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $frame.f2.ok -text OK -command "StdPanelExec::panelOk $globArray"]
  set cn_btn [button $frame.f2.cancel -text Cancel -command "StdPanelExec::panelCancel $this $globArray" \
                                      -state $ca ]
  set ap_btn [button $frame.f2.apply -text Apply -command "StdPanelExec::panelApply $globArray" \
                                     -state $ap ]

  set prev_btn [button $frame.f2.prev -text "<<" -command "StdPanelExec::panelPageMove $globArray -1 " \
                                     -state disabled  -disabledforeground $Info(nonActiveFg) ]

  set next_btn [button $frame.f2.next -text ">>" -command "StdPanelExec::panelPageMove $globArray 1" \
                                      -state disabled -disabledforeground $Info(nonActiveFg)]

  focus $ok_btn
  set theArray(applyButton)  $ap_btn
  set theArray(cancelButton) $cn_btn
  set theArray(prevPageButton) $prev_btn
  set theArray(nextPageButton) $next_btn

}


#--Packing
proc StdPanelCreate::packButtons2Area {frame globArray} {
  global Common Info
  upvar #0 $globArray theArray

  set expnd 1
  set fill none

  pack $frame.f1 $frame.f2 -side left -anchor w -expand 1 -fill x

  pack $frame.f1.l -side left -anchor w -expand 1 -fill x
  pack $frame.f2.ok $frame.f2.cancel $frame.f2.apply -side left -padx 4m 
  pack $frame.f2.next $frame.f2.prev -side right -padx 2m 


}


# NOTE globArray(equationIndex) is updated when this proc is called!
#
proc StdPanelCreate::constructPanelValuesArea {globArray {problem_index ""} {construct_always 0} } {
  global Common Equation
  upvar #0 $globArray theArray

  # If no problem index given, search problem using current equation label!
  #
  if { $problem_index == "" } {

    foreach fld $Equation(allFields) {

      if { ![Equation::isEquation $fld] } {
        continue
      }

      set eq_lbl [DataField::getFieldProperty Equation $fld EquationLabel]

      if { $eq_lbl == $theArray(equationLabel) } {
        set problem_index [DataField::getFieldProperty Equation $fld EquationIndex]
        break
      }
    }
  }

  if { $problem_index == "" || 
       ($problem_index == $theArray(lastEquationIndex) &&
        !$construct_always
       )
     } {
    return
  }

  #-Set header
  if { $theArray(panelsByEquation) } {
    set eq_label [Equation::getFieldPropertyByIndex $problem_index EquationLabel]
    set theArray(valuesHdr) "$eq_label $theArray(valuesHdrBase)"

  } else {
    set eq_label ""
    set theArray(valuesHdr) "$theArray(valuesHdrBase)"
  }

  set theArray(equationLabel) $eq_label
  set theArray(equationIndex) $problem_index
  set theArray(lastEquationIndex) $theArray(equationIndex)

  # Array specific processing
  #
  if { $globArray == "Equation" } {
    Equation::setVariableWidgetsVisibility $problem_index
  }

  set theArray(valuesAreaPacked) 0

  #-Destroy current frame contents
  #
  StdPanelCreate::destroyValuesArea $theArray(valuesFrame) $globArray

  #-Create new contents
  #
  StdPanelCreate::createValuesArea $theArray(valuesFrame) $globArray
  PanelCheck::execPanelFillProcs $globArray
  StdPanelExec::setValuesAreaActivity $globArray $theArray(targetMask)
  StdPanelCreate::packValuesArea $theArray(valuesFrame) $globArray

  # Set label binding for field help
  #
  Widget::setLabelBindings $globArray
}



#======================================#
# Values area and field creation procs #
#======================================#


#--Find maximum size for the values area
#
proc StdPanelCreate::calcValuesAreaSize {globArray} {
  global $globArray Common Info
  upvar #0 $globArray theArray

  #-We use a test frame!
  if { [winfo exists .testf] } {
    destroy .testf
  }

  set f [frame .testf]

  set max_w 0
  set max_h 0

  set has_pages 0

  if { [info exists theArray(panelPageNbr)] } {
    set has_pages 1
  }

  #-Store current values before testing
  if { $has_pages } {
    set current_page_nbr $theArray(panelPageNbr)
  }
  set current_eindex $theArray(equationIndex)

  #-Select indices to test
  if { $theArray(panelsByEquation) } {
    set eindices $theArray(equationIndices)

  } else {
    set eindices 0
  }

  #---Loop all active equations and calculate size
  #   for the values area
  #
  foreach i $eindices {

    if { $theArray(panelsByEquation) } {
      set theArray(equationIndex) $i
    }

    if { $has_pages } {

      set theArray(panelPageNbr) 1

      if { $theArray(panelsByEquation) } {
		    set maxNofPanelPages $theArray($theArray(equationIndex),nofPanelPages)
      } else {
		    set maxNofPanelPages $theArray(nofPanelPages)
      }
    }

    #--Loop all pages (if exists)
    #
    while { 1 } {

      #-Create new contents
      Panel::initFields $globArray "" 0 $i
      StdPanelCreate::createValuesArea $f $globArray
      StdPanelCreate::packValuesArea $f $globArray

      # Pick current packing size and update max-sizes!
      Util::updateGui

      set cur_w [winfo reqwidth $f]
      set cur_h [winfo reqheight $f]

      if { $cur_w > $max_w } {
        set max_w $cur_w
      }

      if { $cur_h > $max_h } {
        set max_h $cur_h
      }
 
      #-Destroy current frame contents
      StdPanelCreate::destroyValuesArea $f $globArray
      
      # Move to next page if exists
      #
      if { !$has_pages } {
        break;
      } else {
        incr theArray(panelPageNbr)

        if { $theArray(panelPageNbr) > $maxNofPanelPages } {
          break;
        }
      }

    } ;# For each page

  } ;# For each equation

  # NOTE: Set small marginals to the estimated sizes
  # to avoid "compressed" fields especially at the bottom
  #
  set theArray(valuesFrame,width) [expr 1.01 * $max_w]
  set theArray(valuesFrame,height) [expr 1.025 * $max_h]

  #-Restore current state
  if { $has_pages } {
    set theArray(panelPageNbr) $current_page_nbr
  }

  if { $theArray(panelsByEquation) } {
    set theArray(equationIndex) $current_eindex
  }
}


# Find nof frames in values area for each equation
# Updates array like Material("sys-index","main-panel-nbr",nofSubPanels)
#
proc StdPanelCreate::setNofValuesAreaFrames {globArray} {
  global Equation
  upvar #0 $globArray theArray

  if { $theArray(panelsByEquation) } {
    set eindices $theArray(equationIndices)
  } else {
    set eindices 0
  }

  foreach eindex $eindices {
		
		for {set i 0} { $i <= 10 } { incr i} {
			set max_frame($i) 0
		}

    if { $theArray(panelsByEquation) } {
      set panel_mask [Equation::getEquationMask $eindex]
    } else {
      set panel_mask ""
    }
		
		set panel_page_nbrs ""

    foreach fld $theArray(allFields) {
		
			# Count nof main panels
			# =====================
			
			# NOTE: A field can be in many panel pages (for convenience)
			if { $theArray(panelsByEquation) } {
        set panel_nbrs [DataField::getPanelPage $globArray $fld $eindex]
        set frame_nbr [DataField::getSubPanel $globArray $fld $eindex]
      } else {
        set panel_nbrs [DataField::getPanelPage $globArray $fld]
        set frame_nbr [DataField::getSubPanel $globArray $fld]
      }
			
			# Add to the list of panel page nbrs
			foreach nbr $panel_nbrs {
				if { -1 == [lsearch $panel_page_nbrs $nbr] } {
					lappend panel_page_nbrs $nbr
				}
			}
			
			# Count nof sub panels
			# ====================

      # Field must belong to the problem
      if { !$theArray($fld,prb) } {
        continue
      }

      if { ![DataField::fieldProblemMaskMatches $globArray $fld $panel_mask] } {
        continue
      }

			foreach nbr $panel_nbrs {
				if { $frame_nbr > $max_frame($nbr) } {
					set max_frame($nbr) $frame_nbr
				}
			}
    }

    if { $theArray(panelsByEquation) } {
      set theArray($eindex,nofPanelPages) [llength $panel_page_nbrs]
    } else {
      set theArray(nofPanelPages) [llength $panel_page_nbrs]
    }
		
		foreach nbr $panel_page_nbrs {
			if { $theArray(panelsByEquation) } {
				set theArray($eindex,$nbr,nofSubPanels) $max_frame($nbr)
			} else {
				set theArray($nbr,nofSubPanels) $max_frame($nbr)
			}
		}
   
  }
}


#===========================#
# Generic field create call #
#===========================#

# Create and pack one field for globArray into frame-argument
# If field is part of a group, then group-argument is non-empty
# Field itself can again be a group entry or a single entry
#
proc StdPanelCreate::createField {globArray group fld fld_frame panelMask eq_index} {
  global Common Info Model $globArray 
  global unit_sp fill_sp selected_$fld 
  upvar #0 $globArray theArray 

  set is_group_entry [Panel::isGroupField $globArray $fld]
  set is_member_entry [Panel::isGroupMemberField $globArray $fld]

  set fieldIndent 3
  set fieldSep 1
 
  set filePadX 3
  set procPadX 3

  set f $fld_frame

  # Screen pad amounts
  set fieldPadX [DataField::getFieldProperty $globArray $fld ScreenPadX $eq_index]
  set fieldPadY [DataField::getFieldProperty $globArray $fld ScreenPadY $eq_index]

  # If not given, apply default, but not
  # if we have a member entry (they are packed
  # denser)

  if { $fieldPadX == "" } {
    if {!$is_member_entry} {
      set fieldPadX $fieldIndent
    } else {
      set fieldPadX 0
    }
  }

  if { $fieldPadY == "" } {
    if {!$is_member_entry} {
      set fieldPadY $fieldSep
    } else {
      set fieldPadY 0
    }
  }

  #--We also collect all field widgets (per group) into
  #  panel and field specific array-variables to be referenced later.

  #--If group entry: recursively using there frames:
  #  On top the frame for the group master field
  #  Below it two frames, Indent-frame and members frame area

  if {$is_group_entry} {

    set groupFields [Panel::getGroupFields $globArray $fld]
    set separateMemberArea [DataField::getFieldProperty $globArray $fld ScreenMemberArea $eq_index]
    set memberIndent [DataField::getFieldProperty $globArray $fld ScreenIndentMembers $eq_index]
    set memberPadX [DataField::getFieldProperty $globArray $fld ScreenPadXMembers $eq_index]
    set memberPadY [DataField::getFieldProperty $globArray $fld ScreenPadYMembers $eq_index]

    if { $memberPadX == "" } {
      set memberPadX 0
    }

    if { $memberPadY == "" } {
      set memberPadY 0
    }

    pack $f -side top -padx $fieldPadX -pady $fieldPadY -fill x -anchor w
    
    set masterFrame [frame $f.master_area]
    set memberFrame [frame $f.member_area]

    if {$separateMemberArea} {
      pack $masterFrame $memberFrame -side top  -fill x -anchor w
    } else {
      #NOTE: We have this strange packing order, because the
      #      natural -left order does not pack these horizontally!!!
      pack $memberFrame $masterFrame -side right -anchor w
    }

    frame $memberFrame.indent
    frame $memberFrame.members

    # Group members indent (dummy label with text as blanks)
    if { $memberIndent > 0 } {
      set member_indent [string range "                    " 1 $memberIndent]
      label $memberFrame.indent.dummy -text $member_indent 
      pack $memberFrame.indent.dummy -side left
    }

    pack $memberFrame.indent $memberFrame.members -side left -padx $memberPadX -pady $memberPadY \
                                                  -fill x -anchor w

    set fm $memberFrame.members

    # Create all fields in the group
    # ==============================
    foreach gv $groupFields {

      if { ![DataField::fieldProblemMaskMatches $globArray $gv $panelMask] } {
        continue
      }

      frame $fm.f$gv

      StdPanelCreate::createField $globArray $fld $gv $fm.f$gv $panelMask $eq_index
    }

    # Back to master's frame before group master continues itself
    set f $masterFrame
  }
  
  #--Three separate frames for label, entries and possible proc/file/use check-box
  frame $f.fl
  frame $f.fe
  frame $f.fc 
  
  #--Store entry field's frame-widget name into a panel specific variable, so
  # it can be referenced later like: globInit(TEMPERATURE,wframe) 
  set theArray($fld,frame) $f
  set theArray($fld,wframe) $f.fe.e$fld

  set FV $fld
  set fieldExpand 0
  
  #---Pick widget and field types
  set wtype   [DataField::getFieldProperty $globArray $FV WidgetType $eq_index]

  #---Field's label 
  set label_string [DataField::getFieldProperty $globArray $FV Label $eq_index]
  set label_width [DataField::getFieldProperty $globArray $FV LabelWidth $eq_index]
  set field_indent [DataField::getFieldProperty $globArray $FV ScreenIndent $eq_index]

  # Create indent * spaces
  set label_indent [string range "                    " 1 $field_indent]

  set label_by_coord [DataField::getFieldProperty $globArray $FV LabelByCoordinate $eq_index]
  set coord_string ""

  # If field is labelled by the coordinate index
  #
  if { $label_by_coord } {
    set coord_index [DataField::getFieldProperty $globArray $FV CoordinateIndex $eq_index]
    set coord_string [Util::getCoordinateLabel $coord_index]

    if { $coord_string != "" } {
      set coord_string "-$coord_string"
    }
  }

  set unit_string  [DataField::getFieldProperty $globArray $FV UnitString $eq_index]
  set f_lbl "$label_string$coord_string $unit_string"
  set f_lbl_indented $label_indent$f_lbl

  label $f.fl.l -text $f_lbl_indented -anchor w -width $label_width

  lappend label_widgets $f.fl.l
  array set $globArray [list allWidgets,label,$FV $f.fl.l]

  #---Field's state
  if { !$theArray($fld,act) } {
    set wdg_state disabled
  } else {
    set wdg_state normal
  }
  set wdg_info [list $wdg_state]

  if {$group == ""} {
    set checkProcName PanelCheck::single
    set checkProcArgument $fld

  } else {
    set checkProcName PanelCheck::group
    set checkProcArgument "$group $fld"
  }
  set cp_info [list $checkProcName $checkProcArgument]

  #---Create widget
  set fld_wdg ""

  # Directory entry with browser button
  # ===================================
  if { $wtype == "BrowsableDirectory" } {
    set fld_wdg [StdPanelCreate::createBrowsableDirectoryField $globArray $FV $f.fe $f.fc $cp_info $wdg_info]

  # File entry with browser button and Abs checkbox
  # ===============================================
  } elseif { $wtype == "BrowsableFile" } {
    set fld_wdg [StdPanelCreate::createBrowsableFileField $globArray $FV $f.fe $f.fc $cp_info $wdg_info]

  # Button
  # ======
  } elseif { $wtype == "Button" } {
    set fld_wdg [StdPanelCreate::createButtonField $globArray $FV $f.fe $cp_info $wdg_info]

  # CheckButton ("Yes/No" button)
  # ===========
  } elseif { $wtype == "CheckButton" } {
    set fld_wdg [StdPanelCreate::createCheckButtonField $globArray $FV $f.fe $cp_info $wdg_info]

  # CheckBox 
  # ========
  } elseif { $wtype == "CheckBox" } {
    set fld_wdg [StdPanelCreate::createCheckBoxField $globArray $FV $f.fe $cp_info $wdg_info]

  # Entry
  # =======
  } elseif { $wtype == "Entry" } {
    set fld_wdg [StdPanelCreate::createEntryField $globArray $FV $f.fe $cp_info $wdg_info]

  # File entry with browser button and Use checkbox
  # ===============================================
  } elseif { $wtype == "IncludeFile" } {
    set fld_wdg [StdPanelCreate::createIncludeFileField $globArray $FV $f.fe $f.fc $cp_info $wdg_info]

  # Listbox 
  # ========
  } elseif { $wtype == "ListBox" } {
    set fld_wdg [StdPanelCreate::createListBoxField $globArray $FV $f.fe $cp_info $wdg_info]

  # OptionMenu
  # ===========
  } elseif { $wtype == "OptionMenu" } {
    set fld_wdg [StdPanelCreate::createOptionMenuField $globArray $FV $f.fe $cp_info $wdg_info]
    set fieldExpand 1

  # RadioButton
  # ===========
  } elseif { $wtype == "RadioButton" } {
    set rb_info [list $f_lbl $group]
    set fld_wdg [StdPanelCreate::createRadioButtonField $globArray $FV $f.fe $cp_info $wdg_info $rb_info]
    set fieldIndent $field_indent

  # RadioButtonOptionMenu
  # =====================
  } elseif { $wtype == "RadioButtonOptionMenu" } {
    set rb_info [list $f_lbl $group]
    set fld_wdg [StdPanelCreate::createRadioButtonOptionMenuField $globArray $FV $f.fe $cp_info $wdg_info $rb_info]

  # SelectionBox
  # ============
  } elseif { $wtype == "SelectionBox" } {
    set fld_wdg [StdPanelCreate::createSelectionBoxField $globArray $FV $f.fe $cp_info $wdg_info]

  # Text
  # =======
  } elseif { $wtype == "Text" } {
    set fld_wdg [StdPanelCreate::createTextField $globArray $FV $f.fe $cp_info $wdg_info]

  # Dummy (for a screen grouper or for a pure label)
  # =====
  # NOTE We need a dummy entry as a pair for the pure label (clumsy, yes)!!!
  #
  } elseif { $wtype == "Label" ||
             $wtype == "Separator" || 
             $wtype == "GroupDataAndScreen" ||
             $wtype == "GroupScreen" 
            } {
    set fld_wdg [StdPanelCreate::createDummyField $globArray $FV $f.fe $cp_info $wdg_info]

  } ;#End widget type


  #--WIDGET PACKING

  #-Packing the widget frame
  set orientation [DataField::getFieldProperty $globArray $FV WidgetPacking $eq_index]
  set anch [DataField::getFieldProperty $globArray $FV WidgetAnchor  $eq_index 0]
  
  #--Field's frame
  #
  if { $wtype == "Label" } {

    if { $fieldPadY < 0 } {
      set fieldPadY 0
    } elseif { $fieldPadY == 0 } {
      set fieldPadY 8
    }

    pack $f -side top -padx $fieldPadX -pady $fieldPadY -fill x -anchor c -expand 0

  } elseif { $wtype == "Separator" } {

    if { $fieldPadY < 0 } {
      set fieldPadY 0
    } elseif { $fieldPadY == 0 } {
      set fieldPadY 2
    }

    pack $f -side top -padx $fieldPadX -pady $fieldPadY -fill x -anchor c -expand 0
    
    # NOTE: This is needed to prevent 'Separator' text to be displayed!
    #
    $f.fl.l configure -text " "

  } elseif {$orientation == "Horizontal"} {
    pack $f -side left -padx $fieldPadX -pady $fieldPadY -fill x -anchor w -expand 0

  } else {
    pack $f -side top -padx $fieldPadX -pady $fieldPadY -fill x -anchor w -expand 0
  }

  #--Field's label (other than radiobuttons)
  #
  if {$wtype != "RadioButton" } {
    
    # Center label field
    if { $wtype == "Label" } {
      pack $f.fl -side right -fill x -expand 1 -anchor c
      pack $f.fl.l -side right -fill x -expand 1 -anchor c 

    # This will pack the field itself at the right
    } elseif { $anch == "" || $anch == "e" } {
      pack $f.fl -side left -fill x -expand 1 -anchor w
      pack $f.fl.l -side left -fill x -expand 1 -anchor w 

    # Leave room at the left for the field itself
    } else {
      pack $f.fl -side left -anchor w
      pack $f.fl.l -side left -anchor w
    }
  }

  #--Field itself
  #
  if { $fld_wdg != "" } {

    set theArray(allWidgets,$FV) $fld_wdg

    if {$wtype == "ListBox" } {
      pack $f.fe -side top -fill x -expand 0 -anchor w
      pack $fld_wdg -side top -fill x -expand 0 -anchor w

    } else {
      # Field to the right
      if { $anch == "" || $anch == "e" } {
        pack $f.fe -side right -fill x -expand 0 -anchor w

      # Field to the left
      } else {
        pack $f.fe -side left -fill x -expand 0 -anchor w
      }
      pack $fld_wdg -side left -fill x -expand 0 -anchor w
    }
  }
}


#=====================#
# Browsable directory #
#=====================#

proc StdPanelCreate::createBrowsableDirectoryField {globArray FV entry_frame cb_frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set padX 3

  #--File entry field
  set wdg [StdPanelCreate::createEntryField $globArray $FV $entry_frame $cp_info $wdg_info]

  #--Directory browser button
  if {1 == [DataField::getFieldProperty $globArray $FV UseAppendMode "" 0]} {
    set browse_wdg  [button $cb_frame.bb -text "..." -width 1 -height 0 \
          -command "Util::fillDirectoryEntry $wdg \
          -array $globArray -panel $globArray -field $FV -append 1" ]

  } else {
    set browse_wdg  [button $cb_frame.bb -text "..." -width 1 -height 0 \
        -command "Util::fillDirectoryEntry $wdg \
        -array $globArray -panel $globArray -field $FV -append 0" ]
  }

  set theArray(allWidgets,$FV,Browse) $browse_wdg

  pack $browse_wdg -in $cb_frame -side right -padx $padX 

  pack $cb_frame -side right -anchor w

  return $wdg
}


#================#
# Browsable file #
#================#

proc StdPanelCreate::createBrowsableFileField {globArray FV entry_frame cb_frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set padX 3

  set has_abs 0
  set has_use 0

  #--File entry field
  set wdg [StdPanelCreate::createEntryField $globArray $FV $entry_frame $cp_info $wdg_info]

  #--File browser button
  if {1 == [DataField::getFieldProperty $globArray $FV HasAbsCheckBox "" 0]} {
    set browse_wdg  [button $cb_frame.bb -text "..." -width 1 -height 0 \
          -command "Util::fillFileEntry $wdg \
                    -array $globArray -absFld $FV,Abs -keepExtFld $FV,Abs -panel $globArray -field $FV "]
  } else {
    set browse_wdg  [button $cb_frame.bb -text "..." -width 1 -height 0 \
          -command "Util::fillFileEntry $wdg \
                    -array $globArray -panel $globArray -field $FV "]
  }

  set theArray(allWidgets,$FV,Browse) $browse_wdg

  #--Abs check box field
  if {1 == [DataField::getFieldProperty $globArray $FV HasAbsCheckBox "" 0]} {

    set has_abs 1
    set f_var $globArray
    append f_var "($FV,Abs)"

    set abs_wdg [checkbutton $cb_frame.ab -text "Abs" -variable $f_var \
             -command "PanelCheck::singleAbsButton $globArray $FV"]

    set theArray(allWidgets,$FV,Abs) $abs_wdg
  }

  #--Use check box
  if {1 == [DataField::getFieldProperty $globArray $FV HasUseCheckBox "" 0]} {

    set has_use 1
    set f_var $globArray
    append f_var "($FV,Use)"

    set use_wdg [checkbutton $cb_frame.ub -text "Use" -variable $f_var \
                 -command "PanelCheck::singleUseButton $globArray $FV"]

    set theArray(allWidgets,$FV,Use) $use_wdg

    if { $theArray(hasParameters) } {
      bind $use_wdg <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $use_wdg {%A %K}"
    } else {
      bind $use_wdg <ButtonRelease-1> "Panel::panelDataChanged 1 $globArray $use_wdg {%A %K}"
    }
  }

  if { $has_abs && $has_use } {
    pack $use_wdg $abs_wdg $browse_wdg -in $cb_frame -side right -padx $padX 
  } elseif { $has_abs } {
    pack $abs_wdg $browse_wdg -in $cb_frame -side right -padx $padX 
  } elseif { $has_use } {
    pack $use_wdg $browse_wdg -in $cb_frame -side right -padx $padX 
  } else {
    pack $browse_wdg -in $cb_frame -side right -padx $padX 
  }

  pack $cb_frame -side right -anchor w

  return $wdg
}


#========#
# Button #
#========#

proc StdPanelCreate::createButtonField {globArray FV frame cp_info wdg_info} {
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  set checkProcArgument [lindex $cp_info 1]

  append checkProcName Button

  set b_state [lindex $wdg_info 0]

  set wl [DataField::getFieldProperty $globArray $FV WidgetLabel]
  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]
  set wh [DataField::getFieldProperty $globArray $FV WidgetHeight]
  set wcp [DataField::getFieldProperty $globArray $FV WidgetCommandProc]

  set wdg [ button $frame.e -state $b_state \
              -text $wl -width $ww -height $wh \
              -command "$checkProcName $globArray $checkProcArgument" ]

  if { $wcp != "" } {
    $wdg configure -command "$wcp $globArray"
  }

  #set ipath $Info(imagePath)
  #set btn_img [ image create photo buttonImage -file [file join $ipath scaledown.gif] ]
  #$wdg configure -image $btn_img
  #$wdg configure -width 8 -height 8 

  return $wdg
}


#==========#
# CheckBox #
#==========#

proc StdPanelCreate::createCheckBoxField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName CheckBox
  set checkProcArgument [lindex $cp_info 1]

  set b_state [lindex $wdg_info 0]

  set wl [DataField::getFieldProperty $globArray $FV WidgetLabel]
  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]
  set wcp [DataField::getFieldProperty $globArray $FV WidgetCommandProc]

  #-Find value for button when it is selected.
  # (this should be defined as a Set-type with two values (on/off))
  set value [lindex [DataField::getFieldProperty $globArray $FV Limits] 1]

  if { $value == "" } {
    set on_val "True"
    set off_val "False"
  } else {
    set on_val [lindex $value 0]
    set off_val [lindex $value 1]
  }

  set b_var $globArray
  append b_var "($FV)"

  set wdg [ checkbutton $frame.e -anchor w \
              -variable $b_var -offvalue $off_val -onvalue $on_val \
              -command "$checkProcName $globArray $checkProcArgument" \
              -indicatoron 1  -width $ww -state $b_state -selectcolor $Info(select_color) ]
  
  # Add label to the widget
  if { $wl != "" } {
    $wdg configure -text $wl
  }

  if { $wcp != "" } {
    $wdg configure -command "$wcp $globArray"
  }

  if { $theArray(hasParameters) } {
    bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"
  } else {
    bind $wdg <ButtonRelease-1> "Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  }

  return $wdg
}


#=============#
# CheckButton #
#=============#

proc StdPanelCreate::createCheckButtonField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName CheckButton
  set checkProcArgument [lindex $cp_info 1]

  set b_state [lindex $wdg_info 0]

  set b_var $globArray
  append b_var "($FV)"

  set wdg [ checkbutton $frame.e -indicatoron 0  -textvariable $b_var \
              -state $b_state -selectcolor $Info(select_color) \
              -variable $b_var -offvalue "N" -onvalue "Y" \
              -command "$checkProcName $globArray $checkProcArgument" \
              -width 0 -height 0 ]

  $wdg configure -width 1 -height 1

  if { $theArray(hasParameters) } {
    bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"
  } else {
    bind $wdg <ButtonRelease-1> "Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  }

  return $wdg
}


#=======#
# Dummy #
#=======#

proc StdPanelCreate::createDummyField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]
  set wh [DataField::getFieldProperty $globArray $FV WidgetHeight]

  set wdg [ frame $frame.e -width $ww -height $wh ]

  return $wdg
}


#=======#
# Entry #
#=======#

proc StdPanelCreate::createEntryField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName Entry
  
  set e_state [lindex $wdg_info 0]

  if { ![Panel::isSimpleDataField $globArray $FV] } {
    set e_relief sunken
  } else {
    set e_relief groove
  }

  if {$e_state != "Normal"} {
    set bg_color $Info(nonActiveBg)
  } else {
    set bg_color white
  }

  set e_var $globArray
  append e_var "($FV)"

  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]

  set wdg [ entry $frame.e -relief $e_relief -textvariable $e_var -width $ww \
              -fg black -bg $bg_color -state $e_state -font $Info(entryFont) ]
   
  if { $theArray(hasParameters) } {
    bind $wdg <KeyRelease> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"
  } else {
    bind $wdg <KeyRelease> "Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  }

  Widget::setEntryBindings "standard" $globArray $FV $wdg
  
  return $wdg
}


#==============#
# Include file #
#==============#

proc StdPanelCreate::createIncludeFileField {globArray FV entry_frame cb_frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName IncludeFile
  set checkProcArgument [lindex $cp_info 1]

  set padX 3

  set has_abs 0
  set has_use 0

  #--File entry field
  set wdg [StdPanelCreate::createEntryField $globArray $FV $entry_frame $cp_info $wdg_info]
  bind $wdg <Double-1> "$checkProcName $globArray $checkProcArgument Double-1"

  #--File browser button
  set browse_wdg  [button $cb_frame.bb -text "..." -width 1 -height 0 \
                    -command "Util::fillFileEntry $wdg -panel $globArray -dropInc 1"]

  set theArray(allWidgets,$FV,Browse) $browse_wdg

  #--Abs check box field
  if {1 == [DataField::getFieldProperty $globArray $FV HasAbsCheckBox "" 0]} {

    set has_abs 1
    set f_var $globArray
    append f_var "($FV,Abs)"

    set abs_wdg [checkbutton $cb_frame.ab -text "Abs" -variable $f_var \
             -command "PanelCheck::singleAbsButton $globArray $FV"]

    set theArray(allWidgets,$FV,Abs) $abs_wdg
  }

  #--Use check box
  if {1 == [DataField::getFieldProperty $globArray $FV HasUseCheckBox "" 0]} {

    set has_use 1
    set f_var $globArray
    append f_var "($FV,Use)"

    set use_wdg [checkbutton $cb_frame.ub -text "Use" -variable $f_var \
                 -command "PanelCheck::singleUseButton $globArray $FV"]

    set theArray(allWidgets,$FV,Use) $use_wdg

    if { $theArray(hasParameters) } {
      bind $use_wdg <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $use_wdg {%A %K}"
    } else {
      bind $use_wdg <ButtonRelease-1> "Panel::panelDataChanged 1 $globArray $use_wdg {%A %K}"
    }
  }

  if { $has_abs && $has_use } {
    pack $use_wdg $abs_wdg $browse_wdg -in $cb_frame -side right -padx $padX 
  } elseif { $has_abs } {
    pack $abs_wdg $browse_wdg -in $cb_frame -side right -padx $padX 
  } elseif { $has_use } {
    pack $use_wdg $browse_wdg -in $cb_frame -side right -padx $padX 
  } else {
    pack $browse_wdg -in $cb_frame -side right -padx $padX 
  }

   pack $cb_frame -side right -anchor e

  return $wdg
}


#=========#
# ListBox #
#=========#

proc StdPanelCreate::createListBoxField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName ListBox
  set checkProcArgument [lindex $cp_info 1]

  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]
  set wh [DataField::getFieldProperty $globArray $FV WidgetHeight]

  set wdg [listbox $frame.e -relief sunken \
             -selectmode browse -exportselection 0 \
             -height $wh -width $ww \
             -yscrollcommand [list $frame.sy$FV set] ]


  scrollbar $frame.sy$FV -orient vertical -command [list $wdg yview]
  pack $frame.sy$FV -side right -fill y  

  bind $wdg <Double-1> "$checkProcName $globArray $checkProcArgument Double-1"
  bind $wdg <ButtonRelease-1> "$checkProcName $globArray $checkProcArgument ButtonRelease-1"

  return $wdg
}


#============#
# OptionMenu #
#============#

proc StdPanelCreate::createOptionMenuField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName OptionMenu
  set checkProcArgument [lindex $cp_info 1]

  set m_state [lindex $wdg_info 0]

  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]

  #-Find value for the menu
  # (this should be defined as a Set-type)
  set values [lindex [DataField::getFieldProperty $globArray $FV Limits] 1]

  # If values is a list of lists
  if { "\{" == [string index [string trim $values] 0] } {
    set values [lindex $values 0]
  }

  # If values is a variable reference
  if { "\$" == [string index [string trim $values] 0] } {
    set values [ DataField::evalInitialValue $values]
  }
  #set theArray($FV) [DataField::getInitialValue $globArray $FV]

  set m_var $globArray
  append m_var "($FV)"

  set wdg $frame.e_om
  set om [Panel::createOptionMenuWidget $wdg $m_var $values]
  set theArray(allWidgets,$FV,om) $om

  set check_proc_call "\"$checkProcName $globArray $checkProcArgument\""

  bind $om <<MenuSelect>> "PanelCheck::checkOptionMenu $om $wdg $globArray $FV $check_proc_call {%s}"

  $wdg configure -activebackground $Info(optionMenuBg) \
                 -indicatoron 0  -width $ww -state $m_state

  # These work incorrectly (launch always when menu opened)
  #bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"
  # These doe not work (does not launch when something is selected)
  #bind $om <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"

  return $wdg
}


#=============#
# RadioButton #
#=============#

proc StdPanelCreate::createRadioButtonField {globArray FV frame cp_info wdg_info rb_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName RadioButton
  set checkProcArgument [lindex $cp_info 1]

  set b_state [lindex $wdg_info 0]

  set b_lbl [lindex $rb_info 0]
  set b_group [lindex $rb_info 1]
  set checkProcArgument [lindex $cp_info 1]
 
  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]

  #-Use radio-button group parent variable as the button-group variable
  set var $globArray
  append var "($b_group)"

  #-If needed, set initial value for the variable (from grpup parent initial value)
  if { ![info exists $var ] } {
    set $var [DataField::getInitialValue $globArray $b_group]
  }

  #-Value for button when it is selected.
  # (this should be defined as a Constant-type) 
  set value [lindex [DataField::getFieldProperty $globArray $FV Limits] 1]

  set wdg [radiobutton $frame.e -anchor w \
            -indicatoron 1 -width $ww -state $b_state -text $b_lbl \
            -selectcolor $Info(select_color) -variable $var -value $value \
            -command "$checkProcName $globArray $checkProcArgument"]

  if { $theArray(hasParameters) } {
    bind $wdg <ButtonRelease-1> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"
  } else {
    bind $wdg <ButtonRelease-1> "Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  }

  return $wdg
}


#======================#
# RadiButtonOptionMenu #
#======================#

proc StdPanelCreate::createRadioButtonOptionMenuField {globArray FV frame cp_info wdg_info rb_info} {
  global Info
  upvar #0 $globArray theArray

  # RadioButton and OptionMenu fields
  set rb_fld [DataField::getFieldProperty $globArray $FV RadioButtonField "" 0]
  set om_fld [DataField::getFieldProperty $globArray $FV OptionMenuField "" 0]

  if { $rb_fld == "" } {
    set rb_fld [lindex $rb_info 1]
  }

  if { $om_fld == "" } {
    set om_fld $FV
  }

  # RadioButton and OptionMenu values
  set limits [DataField::getFieldProperty $globArray $FV Limits]

  set rb_value [lindex [lindex $limits 1] 0] ;# RadioButton value
  set om_list [lindex [lindex $limits 1] 1]  ;# OptionMenu list

  #-Radio button
  set rb_limits [list "Set" $rb_value]
  set rb_info [list "" $rb_fld]
  set rb_cp_info [list [lindex $cp_info 0] $rb_fld]

  set lims [DataField::getFieldProperty $globArray $rb_fld Limits] ;# Backup Limits
  DataField::setFieldProperty $globArray $rb_fld Limits $rb_limits ;#-Needed for constructing the wdg!
  set rb_wdg [StdPanelCreate::createRadioButtonField $globArray $rb_fld $frame $rb_cp_info $wdg_info $rb_info]
  $rb_wdg configure -state normal
  DataField::setFieldProperty $globArray $rb_fld Limits $lims ;# Restore Limits!

  #-Option menu
  set om_limits [list "Set" [list $om_list]]
  set om_cp_info [list [lindex $cp_info 0] $om_fld]

  set lims [DataField::getFieldProperty $globArray $om_fld Limits] ;# Backup Limits
  DataField::setFieldProperty $globArray $om_fld Limits $om_limits ;#-Needed for constructing the wdg!
  set om_wdg [StdPanelCreate::createOptionMenuField $globArray $om_fld $frame $om_cp_info $wdg_info]
  DataField::setFieldProperty $globArray $om_fld Limits $lims ;# Restore Limits!
  
  pack $frame -side right -fill x -expand 0 -anchor w
  pack $rb_wdg $om_wdg -side left -expand 0 -anchor w
  
  #set theArray(allWidgets,$rb_fld) $rb_wdg
  #set theArray(allWidgets,$om_fld) $om_wdg

  set theArray(allWidgets,$FV) [list $rb_wdg $om_wdg]
  
  return ""
}


#==============#
# SelectionBox # (a listbox with selection indicators)
#==============#

proc StdPanelCreate::createSelectionBoxField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set checkProcName [lindex $cp_info 0]
  append checkProcName SelectionBox
  set checkProcArgument [lindex $cp_info 1]

  set wl [DataField::getFieldProperty $globArray $FV WidgetLabel]
  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]
  set wh [DataField::getFieldProperty $globArray $FV WidgetHeight]

  set wdg [ listbox $frame.e -relief sunken \
              -selectmode browse -exportselection 0 \
              -height $wh -width $ww -font $Info(tableFont) \
              -yscrollcommand [list $frame.sy set]]

  scrollbar $frame.sy -orient vertical -command [list $wdg yview]

  set theArray(allWidgets,$FV) $wdg

  # WidgetLabel (if defined) is put above the box!
  if { $wl != "" } {
    label $frame.l -text $wl
    pack $frame.l -side top
  }

  pack $frame.sy -side right -fill y  

  #-Find values to be dislayed
  set values ""
  set use_fill_proc 0

  # (this should be defined as a Set-type)
  set values [lindex [DataField::getFieldProperty $globArray $FV Limits] 1]

  # If values is a list of lists
  if { "\{" == [string index [string trim $values] 0] } {
    set values [lindex $values 0]

  # If values is a variable reference
  } elseif { "\$" == [string index [string trim $values] 0] } {
    set values [ DataField::evalInitialValue $values]

  # If values is by a procedure
  } elseif { "\@" == [string index [string trim $values] 0] } {
    set use_fill_proc 1
    #set values [ DataField::evalInitialValue $values]
  }

  #-If list values given
  if { $values != "" } {
    #-Call a proc to fill the values
    if { $use_fill_proc } {
      DataField::evalInitialValue $values
    #-Use data list
    } else {
      ListBox::fill $wdg $values
    }
  }

  set check_proc_call "\"$checkProcName $globArray $checkProcArgument\""

  bind $wdg <ButtonRelease-1> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call ButtonRelease-1 %x %y"
  bind $wdg <ButtonRelease-2> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call ButtonRelease-3 %x %y"
  bind $wdg <ButtonRelease-3> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call ButtonRelease-3 %x %y"
  bind $wdg <Control-ButtonRelease-1> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call Control-ButtonRelease-1 %x %y"
  bind $wdg <Control-ButtonRelease-2> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call Control-ButtonRelease-3 %x %y"
  bind $wdg <Control-ButtonRelease-3> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call Control-ButtonRelease-3 %x %y"
  bind $wdg <Double-1> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call Control-ButtonRelease-3 %x %y"

  bind $wdg <KeyPress-Return> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call ButtonRelease-1"
  bind $wdg <KeyPress-space> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call ButtonRelease-3"
  bind $wdg <Control-KeyPress-Return> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call Control-ButtonRelease-1"
  bind $wdg <Control-KeyPress-space> "PanelCheck::checkSelectionBox $wdg $globArray $FV $check_proc_call Control-ButtonRelease-3"
  
  return $wdg
}


#======#
# Text #
#======#

proc StdPanelCreate::createTextField {globArray FV frame cp_info wdg_info} {
  global Info
  upvar #0 $globArray theArray

  set t_state [lindex $wdg_info 0]

  if {$t_state != "Normal"} {
    set bg_color $Info(nonActiveBg)
  } else {
    set bg_color white
  }

  set ww [DataField::getFieldProperty $globArray $FV WidgetWidth]
  set wh [DataField::getFieldProperty $globArray $FV WidgetHeight]

  set wdg [text $frame.e -relief groove -height $wh -width $ww \
                            -bg $bg_color -state $t_state \
                            -font $Info(entryFont)]

  $wdg insert 0.0 $theArray($FV)

  if { $theArray(hasParameters) } {
    bind $wdg <KeyRelease> "Panel::panelDataModified 1 $globArray $wdg {%A %K}"
  } else {
    bind $wdg <KeyRelease> "Panel::panelDataChanged 1 $globArray $wdg {%A %K}"
  }

  return $wdg
}



# end ecif_tk_panelCreate.tcl
# ********************
