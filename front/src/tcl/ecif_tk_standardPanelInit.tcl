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
#Module:    ecif_tk_panelInit.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Initialization procedures for "standard" panels.
#
#************************************************************************


################################
#### Definition of bindings ####
################################

proc StdPanelInit::createWidgetBindings {globArray} {
  upvar #0 $globArray p
  global Common Info

  set selected_row $Info(NO_INDEX)

  if { [info exist p(objectLB)] && [winfo exists $p(objectLB)] } {
    bind $p(objectLB) <ButtonRelease-1> "execNsProc $globArray objectSelected $selected_row "
    bind $p(objectLB) <Double-1> "execNsProc $globArray objectDblSelected"
  }

  set p(realButtonPress) 1
  set p(realButtonRelease) 1

  if { [info exist p(boundaryLB)] && [winfo exists $p(boundaryLB)] } {
    bind $p(boundaryLB) <ButtonRelease-1> "execNsProc $globArray boundarySelected $selected_row "
    bind $p(boundaryLB) <Double-1> "execNsProc $globArray boundaryDblSelected"
    bind $p(boundaryLB) <Control-Double-1> "execNsProc $globArray boundaryCntrlDblSelected"
  }

  if { [info exist p(parameterLB)] && [winfo exists $p(parameterLB)] } {
    bind $p(parameterLB) <ButtonRelease-1> "execNsProc $globArray parameterSelected $selected_row "
    bind $p(parameterLB) <Double-1> "execNsProc $globArray parameterDblSelected"
    bind $p(parameterLB) <Control-Double-1> "execNsProc $globArray parameterCntrlDblSelected"
  }
}



##############################################
#### Data and ListBoxLists initialization ####
##############################################

proc StdPanelInit::disablePanel { globArray } {
  upvar #0 $globArray theArray

  set theArray(addAllowed) 0
  set theArray(updateAllowed) 0
  
  set btns {
    panelAddButton panelUpdateButton panelDeleteButton
    cancelButton applyButton
  }

  foreach btn $btns {

    if { [info exists theArray($btn)]  &&
         [winfo exists $theArray($btn)]
       } {
      Widget::configureS $theArray($btn) disabled
    }
  }
}


#-Initialize panel data.
#
proc StdPanelInit::initPanelData {globArray {reset 1}} {
  global Model ObjectTable
  upvar #0 $globArray theArray

  set theArray(dataDirty) 0
  set theArray(dataChanged) 0
  set theArray(dataModified) 0
  set theArray(autoAdd) 0
  set theArray(updateAllowed) 0

  if { $ObjectTable(ids) == "" } {
    StdPanelInit::disablePanel $globArray
    return
  }

  set Model(DB,willNeedUpdate) 0

  #---Construct all necessary lists from cpp-delivered strings.
  StdPanelInit::createPanelData $globArray

  #---Panel specific pre initialization
  switch $globArray {

    BoundaryCondition {
      #Util::cpp_exec resetBoundarySelections
    }

    GridH {
      GridH::preInitData
    }

    GridParameter {
      GridParameter::preInitData
    }
  }

  #---Fill list boxes etc.
  StdPanelInit::fillPanelData $globArray

  #---Select possible previously selected object, element and parameter
  #
  # NOTE: We have to store previous selections before
  # we apply any selection, because these
  # the global variables are changed during the selections
  #
  set obj_row $theArray(objectIndex)
  set bndr_row $theArray(boundaryIndex)
  set param_row $theArray(parameterIndex)
  set theArray(objectIndex,old) $obj_row
  set theArray(boundaryIndex,old) $bndr_row
  set theArray(parameterIndex,old) $param_row 

  if { [info exists theArray(objectLB)] && $obj_row != "" } {

    $theArray(objectLB) selection clear 0 end
    set rc [execNsProc $globArray objectSelected $obj_row]

    # If array's ns-proc not defined, exec generic proc
    if { !$rc } {
      StdPanelExec::objectSelected $globArray $obj_row
    }
  }

  if { [info exists theArray(boundaryLB)] && $bndr_row != "" } {
    $theArray(boundaryLB) selection clear 0 end
    StdPanelExec::boundarySelected $globArray $bndr_row
  }

  if { [info exists theArray(parameterLB)] && $param_row != "" } {
    $theArray(parameterLB) selection clear 0 end
    StdPanelExec::parameterSelected $globArray $param_row
  }


  #---Get possible pending selections
  #cpp_setSelectionsToGui

  #---Panel specific post initialization
  switch $globArray {

    BoundaryCondition {
      BoundaryCondition::postInitData
    }

    Equation {
      Equation::postInitData
    }

    GridParameter {
      GridParameter::postInitData
    }
  }
}


#-Create panel data
#
proc StdPanelInit::createPanelData {globArray} {
  global Info
  upvar #0 $globArray theArray

  set OS $Info(objectSeparator)
  set FS $Info(fieldSeparator)

  set theArray(mode) ""

  if { ![info exists theArray(objectIndex)] } {
    set theArray(objectIndex) 0
  }

  set theArray(targetIndex) $theArray(objectIndex)
  
  if { $theArray(hasBoundaries) } {
    
    if { ![info exists theArray(boundaryIndex)] } {
      set theArray(boundaryIndex) 0
      set theArray(targetIndex) 0
    }

    set theArray(targetIndex) $theArray(boundaryIndex)

  } else {
    set theArray(boundaryIndex) ""
  }

  if { ![info exists theArray(parameterIndex)] } {
    set theArray(parameterIndex) ""
  }

  StdPanelInit::createObjectAndBoundaryLists $globArray
  
  if { $theArray(hasParameters) } {
    StdPanelInit::createParameterLists $globArray
  }
}


#-Create object and boundary data lists data
#
proc StdPanelInit::createObjectAndBoundaryLists {globArray} {
  global Model ObjectTable
  upvar #0 $globArray theArray

  set tp $theArray(parameterType)

  set theArray(objectTypes) ""
  set theArray(subObjectTypes) ""
  set must_have_subs 0

  #-Equation etc.: Body level stuff
  if { $tp == "EQ" || $tp == "BF" || $tp == "IC" || $tp == "MT"} {
    set theArray(objectTypes) "B"

  #-Body level stuff (body parameters)
  } elseif { $tp == "BodyP"} {
    set theArray(objectTypes) "B"

  #-Boundaries
  # NOTE: Do not accept Intra layer boundaries (ILB)!
  #
  } elseif { $tp == "BA" || $tp == "BndrP" } {
    set theArray(objectTypes) "B BP"

    if { $Model(GEOMETRY_DIMENSION) == "3D" } {
      set theArray(subObjectTypes) { {F {-ILB}} }
    } else {
      set theArray(subObjectTypes) { {E {-ILB}} }
    }
    set must_have_subs 1

  #-Boundary condition
  # NOTE: Do not accept Intra layer boundaries (ILB)!
  #
  } elseif { $tp == "BC" } {
    set theArray(objectTypes) "B BP"
    set theArray(subObjectTypes) { {FG {-ILB}} {EG {-ILB}} {VG {-ILB}} }
    set must_have_subs 1

  #-Mesh structure (GridParameter)
  } elseif { $tp == "GR" } {
    set theArray(objectTypes) { {BL {CSD OPN}} }

  #-Mesh density (GridH)
  } elseif { $tp == "GH" } {
    set theArray(objectTypes) { {B {CSD OPN}} {BP} }
    set theArray(subObjectTypes) "F E V"
    set must_have_subs 1
  }

  # Object lists (always created)
  # ============
  StdPanelInit::findAllObjectIds $globArray $theArray(objectTypes) $must_have_subs
  set theArray(allObjectIds) [Object::sortIdsByParentTags $theArray(allObjectIds)]
  set theArray(objectIds) $theArray(allObjectIds)

  StdPanelInit::constructObjectLBList $globArray

  # Objects are targets
  # -------------------
  if { !$theArray(hasBoundaries) } {
    set theArray(allTargetIds) $theArray(allObjectIds)
    set theArray(targetIds) $theArray(allObjectIds)
  }

  if { $theArray(hasBoundaries) } {
    StdPanelInit::findAllBoundaryIds $globArray $theArray(subObjectTypes)
    set theArray(allTargetIds) $theArray(allBoundaryIds)
    set theArray(targetIds) ""
    set theArray(boundaryIds) ""
    StdPanelInit::constructBoundaryLBList $globArray
  }

  #-Here we store old parameter attach values
  if { $theArray(hasParameters) } {

    set pv [string tolower $tp]
 
    foreach id $ObjectTable(ids) {
      if { [info exists ObjectTable($id,$pv)] } {
        set ObjectTable($id,$pv,old) $ObjectTable($id,$pv)
      }
    }
  }

}


# Create all parameter data list
#
proc StdPanelInit::createParameterLists {globArray} {
  global Common Info ObjectTable
  upvar #0 $globArray theArray

  #-Check if panel has paramter ids!
  if { ![info exists theArray(ids)] } {
    return
  }

  # Delete old backup data
  # ======================
  Util::unsetArrayIdVariables $globArray $theArray(ids) *,old

  # Create new backup data
  # =====================
  set theArray(ids,old) $theArray(ids)

  # All other variables but mask
  foreach pid $theArray(ids) {
    foreach vn {data name oid} {

      if { [info exists theArray($pid,$vn)] } {
        set theArray($pid,$vn,old) $theArray($pid,$vn)
      }
    }
  }

  # Masks need special handling!
  foreach id $theArray(ids) {

    #-Parameter get masks by objects, but not Equation!
    #
    if { $globArray != $Common(panelArray,EQ) } {

      if { $theArray(hasMaskedParameter) } {

        set oid $theArray($id,oid)

        if { $oid != $Info(NO_INDEX) } {
          set theArray($id,mask) [StdPanelExec::getTargetMask $globArray $oid]
        } else {
          set theArray($id,mask) ""
        }

      } else {
        set theArray($id,mask) ""
      }

    }
  
    set theArray($id,mask,old) $theArray($id,mask)
  }

  # Listbox list
  StdPanelInit::constructParameterLBList $globArray
}


proc StdPanelInit::constructLBLists {globArray} {
  upvar #0 $globArray theArray

  if { $theArray(hasBodies) } {
    StdPanelInit::constructObjectLBList $globArray
  }

  if { $theArray(hasBoundaries) } {
    StdPanelInit::constructBoundaryLBList $globArray
  }

  if { $theArray(hasParameters) } {
    StdPanelInit::constructParameterLBList $globArray
  }
}


proc StdPanelInit::constructObjectLBList {globArray} {
  upvar #0 $globArray theArray

  #---Objects (body or body-pair)
  set theArray(objectLBList) ""

  if {$theArray(hasBoundaries)} {
    set as_target 0
  } else {
    set as_target 1
  }

  foreach id $theArray(objectIds) {
    set row [StdPanelInit::formObjectLBRow $globArray $id $as_target]
    lappend theArray(objectLBList) $row 
  }
}


proc StdPanelInit::constructBoundaryLBList {globArray} {
  upvar #0 $globArray theArray

  #---Objects (body or body-pair)
  set theArray(boundaryLBList) ""

  set as_target 1

  foreach id $theArray(boundaryIds) {
    set row [StdPanelInit::formObjectLBRow $globArray $id $as_target]
    lappend theArray(boundaryLBList) $row 
  }
}



proc StdPanelInit::constructBoundaryLBList2 {globArray} {
  upvar #0 $globArray theArray

  #---Objects (body or body-pair)
  set theArray(boundaryLBList) ""

  set as_target 1

  foreach id $theArray(boundaryIds) {
    set tp [Object::getType $id]
    if { $tp != "F" && $tp != "FG" } continue
    set row [StdPanelInit::formObjectLBRow $globArray $id $as_target]
    lappend theArray(boundaryLBList) $row 
  }
  
  foreach id $theArray(boundaryIds) {
    set tp [Object::getType $id]
    if { $tp != "E" && $tp != "EG" } continue
    set row [StdPanelInit::formObjectLBRow $globArray $id $as_target]
    lappend theArray(boundaryLBList) $row 
  }

  foreach id $theArray(boundaryIds) {
    set tp [Object::getType $id]
    if { $tp != "V" && $tp != "VG" } continue
    set row [StdPanelInit::formObjectLBRow $globArray $id $as_target]
    lappend theArray(boundaryLBList) $row 
  }
}


# Create parameter listbox list
#
proc StdPanelInit::constructParameterLBList {globArray} {
  upvar #0 $globArray theArray

  #-Ids and data for the parameters themselves
  set theArray(parameterLBList) ""

  foreach id $theArray(ids) {
    lappend theArray(parameterLBList) [StdPanelInit::formParameterLBRow $globArray $id] 
  }
}


#-Set panel data values 
#
proc StdPanelInit::fillPanelData {globArray} {
  global Info
  upvar #0 $globArray theArray

  set tp $theArray(parameterType)
 
  set theArray(objectId) ""
  set theArray(objectName) ""
  set theArray(targetMask) ""

  set theArray(boundaryIds) ""
  set theArray(boundaryIndices) ""

  set theArray(parameterId) ""
  set theArray(parameterName) ""
  set theArray(parameterMask) ""

  set theArray(body1Id) ""
  set theArray(body2Id) ""

  set theArray(tableEntryPanelIds) ""
  set theArray(procedureEntryPanelIds) ""

  if { $theArray(hasBodies) } {
    StdPanelInit::fillObjectListBox $globArray
  }

  if { $theArray(hasBoundaries) } {
    StdPanelInit::fillBoundaryListBox $globArray
  }

  if { !$theArray(hasParameters) } {
    return
  }

  StdPanelInit::fillParameterListBox $globArray

  #---Set initial values into arrays's fields
  StdPanelInit::initFieldValues $globArray

  #---Other initial values
  set theArray(nextNewParameterId,old) $theArray(nextNewParameterId)
}


proc StdPanelInit::fillLBLists {globArray} {
  upvar #0 $globArray theArray

  if { $theArray(hasBodies) } {
    ListBox::fill $theArray(objectLB) $theArray(objectLBList)
  }

  if { $theArray(hasBoundaries) } {
    ListBox::fill $theArray(boundaryLB) $theArray(boundaryLBList)
  }

  if { $theArray(hasParameters) } {
    ListBox::fill $theArray(parameterLB) $theArray(parameterLBList)
  }
}


proc StdPanelInit::fillObjectListBox {globArray} {
  upvar #0 $globArray theArray

  ListBox::fill $theArray(objectLB) $theArray(objectLBList)
}


proc StdPanelInit::fillBoundaryListBox {globArray} {
  upvar #0 $globArray theArray

  ListBox::fill $theArray(boundaryLB) $theArray(boundaryLBList)
}


proc StdPanelInit::fillParameterListBox {globArray} {
  upvar #0 $globArray theArray

  #-NOTE: Parameters are set unmarked when panel is opened
  set new_list ""

  foreach row $theArray(parameterLBList) {
    set row [string replace $row 0 0 " "]
    lappend new_list $row
  }

  set theArray(parameterLBList) $new_list

  ListBox::fill $theArray(parameterLB) $theArray(parameterLBList)

}


# Pick all objects ids whose type belong to "types"
#
# NOTE: Updates array variable theArray(allObjectIds)
#
proc StdPanelInit::findAllObjectIds {globArray types {must_have_subs 0} } {
  global ObjectTable
  upvar #0 $globArray theArray

  set theArray(allObjectIds) ""

  foreach type $types {

    set otp [lindex $type 0] ;# Object type
    set ocl [lindex $type 1] ;# Object class (ie. subtype like VIR(tual))
    
    set ids [Object::getIds $otp $ocl]

    foreach id $ids {

      if { $must_have_subs &&
           $ObjectTable($id,sbIds) == ""
         } {
        continue
      }

      lappend theArray(allObjectIds) $id
    }
  }

}


# Pick all boundary ids whose type belong to "types"
#
proc StdPanelInit::findAllBoundaryIds { globArray types } {
  global ObjectTable
  upvar #0 $globArray theArray

  set theArray(allBoundaryIds) ""

  foreach type $types {

    set btp [lindex $type 0] ;# Boundary type
    set bcl [lindex $type 1] ;# Boundary class (ie. subtype like VIR(tual))

    set ids [Object::getIds $btp $bcl]

    foreach id $ids {
      lappend theArray(allBoundaryIds) $id
    }
  }
}


# Cover proc to select activate target list format
#
proc StdPanelInit::formTargetLBRow {globArray obj_id} {
  upvar #0 $globArray theArray

  set as_target 1

  StdPanelInit::formObjectLBRow $globArray $obj_id $as_target
}


proc StdPanelInit::formObjectLBRow {globArray obj_id as_target} {
  global Info Model ObjectTable
  upvar #0 $globArray theArray

  set otp $ObjectTable($obj_id,tp)

  set row "("

  # Vertex/Edge as subelement
  # =========================
  if { $otp == "V" || ( $otp == "E" && $Model(GEOMETRY_DIMENSION) == "3D" )
     } {
    append row $otp
    append row $ObjectTable($obj_id,tg)

  # Vertex group as subelement
  # ==========================
  } elseif { $otp == "VG" } {
    if { 1 < [llength $ObjectTable($obj_id,mbrIds)] } {
      append row "VG"

    } else {
      append row "V"
    }
    append row $ObjectTable($obj_id,tg)

  # Edge group as subelement
  # ========================
  } elseif { $otp == "EG" && $Model(GEOMETRY_DIMENSION) == "3D" } {
    if { 1 < [llength $ObjectTable($obj_id,mbrIds)] } {
      append row "EG"

    } else {
      append row "E"
    }
    append row $ObjectTable($obj_id,tg)

  # Bodypair
  # ========
  } elseif { $otp == "BP" } {

    set id1 $ObjectTable($obj_id,pr1Id)
    set id2 $ObjectTable($obj_id,pr2Id)

    append row $ObjectTable($id1,tg)
    append row ","
    append row $ObjectTable($id2,tg)

  # Others
  # ======
  } else {
    append row $ObjectTable($obj_id,tg)
  }

  append row ") "
  append row $ObjectTable($obj_id,nm)

  # If this is not a target object table, we can
  # stop here!
  #
  if {!$as_target} {
    return $row
  }

  if { !$theArray(hasParameters) } {
    return $row
  }

  # Name for paramter variable like bf, eq etc
  set pv [string tolower $theArray(parameterType)]

  if { ![info exists ObjectTable($obj_id,$pv)] ||
        $ObjectTable($obj_id,$pv) == $Info(NO_INDEX)
     } {
    return $row
  }

  set pid $ObjectTable($obj_id,$pv)

  set pname $theArray($pid,name)

  append row " - "
  append row "\[$pid\] "
  append row $pname

  return $row
}


# Procedure forms one parameter listbox line like: "(1) Inflow"
#
proc StdPanelInit::formParameterLBRow { globArray pid { show_mark 1} } {
  global Info
  upvar #0 $globArray theArray

  #  For equation-type dependent parameters, we also mark, that parameter
  #  fits current object's equation type
  if { $theArray(hasMaskedParameter) && $show_mark } {

    set pline $Info(paramMarked)($pid)

  } else {
    set pline $Info(paramUnmarked)($pid)
  }

  append pline " "
  append pline [StdPanelExec::getParameterName $globArray $pid]
}


proc StdPanelExec::getParameterName {globArray pid} {
  upvar #0 $globArray theArray

  if { [info exists theArray($pid,name)] } {
    return $theArray($pid,name)

  } else {
    return ""
  }
}


# Create equation parameter masks from equation data.
#
# These masks are used as target mask (for bodies when
# an eqaution is attached)
#
# Update data according to check proc result and store
#  new data if it is modified.
#
proc StdPanelInit::updateEquationDataAndMasks {} {
  global Equation ObjectTable Info

  # Equation masks
  # ==============
  Equation::initMasks
  
  #--Form "normal" masks

  # Turn off values area update (avoid loop!)
  set flag $Equation(updateValuesArea)
  set Equation(updateValuesArea) 0

  set modified 0

  foreach id $Equation(ids) {

    set Equation(targetMask) ""

    Panel::resetFields Equation

    StdPanelInit::initFieldValues Equation

    DataField::formDataFields Equation $Equation($id,data)

    # Set variableType masks by equations
    # -----------------------------------
    Equation::updateVariableTypeTargetMasks

    # Predefined equations
    # --------------------
    Equation(NAVIER-STOKES,fillProc)
    Equation(HEAT_EQUATION,fillProc)
    Equation(STRESS_ANALYSIS,fillProc)
    Equation(ADVECTION_DIFFUSION_EQUATION,fillProc)

    # Possible user defined equations
    # -------------------------------
    foreach fld $Equation(allFields) {

      if { 1 != [DataField::getFieldProperty Equation $fld IsUserDefinedField "" 0] ||
           ![Equation::isEquation $fld]
         } {
        continue
      }

      Equation::setUserDefinedEquationMask $fld
    }


    # Update parameter data
    # ---------------------
    set tmp $Equation($id,data)

    set Equation($id,data) [DataField::formParameterLine Equation]
    set Equation($id,mask) $Equation(targetMask)

    if { $tmp != $Equation($id,data) } {
      set modified 1
    }

    Panel::unsetAllFieldData Equation
  }

  # Save data if modified
  # =====================
  if {$modified} {
    StdPanelExec::panelSave Equation 1
  }

  # Restore original flag value
  set Equation(updateValuesArea) $flag
} 


# Procedure creates a sparse query list: Arguments are largest id,
# list of existing ids and their corresponding values.
# Size of the list is max_id + 1 <--> ids should be [0...max_id]
# (or some subset of these)
#
proc StdPanelInit::constructSparseQueryList {max_id ids_list values_list} {

  set result ""
  set counter 0

  while {$counter <= $max_id } {
    lappend result ""
    incr counter
  }

  foreach id $ids_list value $values_list {
    set result [lreplace $result $id $id $value] 
  }

  return $result
}


proc StdPanelInit::setObjectTableSubIds {objType subType {keep_current 0} } {
  global ObjectTable

  #---Loop object table to find searched objects of the "objType"
  foreach oid $ObjectTable(ids) {

    if { $ObjectTable($oid,tp) != $objType } {
      continue
    }
    
    # If current sub-ids should be accpeted if defined
    #
    if { $keep_current && $ObjectTable($oid,sbIds) != "" } {
      continue
    }

    set subIds ""
    set subTgs ""

    #--Loop all subobjects and find they rows using their id as the key
    foreach stg $ObjectTable($oid,sbTgs) {
     
      #-Loop again object table to find the subobject row
      foreach sid $ObjectTable(ids) {

        if { $ObjectTable($sid,tp) != $subType ||
             $ObjectTable($sid,tg) != $stg
           } {
          continue
        }

        lappend subIds $sid
        lappend subTgs $stg
      }
    }

    set ObjectTable($oid,sbIds) $subIds
    set ObjectTable($oid,sbTgs) $subTgs
  }
}


# Construct ObjectTable masks for bodies and their subelemnts
#
proc StdPanelInit::createObjectTableMasks {} {
  global Model

  StdPanelInit::createObjectTableBodyMasks
  StdPanelInit::createObjectTableSubElementMasks
}


# Construc ObjectTable body masks from equation masks
#
proc StdPanelInit::createObjectTableBodyMasks {} {
  global Equation Info Model ModelFlags ObjectTable

  # Base mask from geometry
  set gm "$Equation(dimgType) $Equation(dimsType)"

  # Bodies
  # ======
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    set eqid $ObjectTable($id,eq)
    set em ""

    #-If equation defined for the body group
    if { $eqid != $Info(NO_INDEX) } {
      set em $Equation($eqid,mask)
    }

    set ObjectTable($id,msk) [Util::combineMasks $gm $em]
  }

  # BodyPairs from Bodies
  # =====================
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "BP" } {
      continue
    }

    set o1 $ObjectTable($id,pr1Id)
    set o2 $ObjectTable($id,pr2Id)

    set m1 $ObjectTable($o1,msk)
    set m2 $ObjectTable($o2,msk)
    set m [Util::combineMasks $m1 $m2]

    set ObjectTable($id,msk) $m
  }

}


# Construc ObjectTable subelement masks by processing objects and
# their possible subobjects hierachically
#
proc StdPanelInit::createObjectTableSubElementMasks {} {
  global Equation Info Model ObjectTable

  # Reset subelement masks
  # ======================
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) == "B"  ||
         $ObjectTable($id,tp) == "BP"
       } {
      continue
    }

    set ObjectTable($id,msk) ""
  }

  # Body and BodyPair subelements (Faces or Edges)
  # =============================
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B"  &&
         $ObjectTable($id,tp) != "BP"  
       } {
      continue
    }
    
    set m1 $ObjectTable($id,msk)
    foreach sid $ObjectTable($id,sbIds) {

      if { [info exists ObjectTable($sid,msk)] } {
        set m2 $ObjectTable($sid,msk)
        set m [Util::combineMasks $m1 $m2]
        set ObjectTable($sid,msk) $m
      }
    }

  }

  # Edges from Faces
  # ================
  if { $Model(GEOMETRY_DIMENSION) == "3D" } {

    foreach id $ObjectTable(ids) {

      if { $ObjectTable($id,tp) != "F" &&
           $ObjectTable($id,tp) != "FG" 
         } {
        continue
      }

      set m1 $ObjectTable($id,msk)

      foreach sid $ObjectTable($id,sbIds) {

        if { [info exists ObjectTable($sid,msk)] } {
          set m2 $ObjectTable($sid,msk)
          set m [Util::combineMasks $m1 $m2]
          set ObjectTable($sid,msk) $m
        }
      }

    }
  }

  # Vertices from Edges
  # ===================
  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "E" &&
         $ObjectTable($id,tp) != "EG"
       } {
      continue
    }
    
    set m1 $ObjectTable($id,msk)

    foreach sid $ObjectTable($id,sbIds) {

      if { [info exists ObjectTable($sid,msk)] } {
        set m2 $ObjectTable($sid,msk)
        set m [Util::combineMasks $m1 $m2]
        set ObjectTable($sid,msk) $m
      }
    }

  }

  # Add BC masks for Boundaries
  # ===========================
  foreach id $ObjectTable(ids) {

    set otp $ObjectTable($id,tp)

    if { $otp != "FG" && $otp != "EG" && $otp != "VG" } {
      continue
    }

    # Face group
    if { $otp == "FG" } {
      append ObjectTable($id,msk) " 刃cD 刃cN3"

    # Edge group
    } elseif { $otp == "EG" } {
      append ObjectTable($id,msk) " 刃cD 刃cN2"

    # Vertex group
    } else {
       append ObjectTable($id,msk) " 刃cD"
    }
    
  }
    
}


#====================================#
#  BoundaryCondition panel specific  #
#====================================#

proc BoundaryCondition::postInitData {} {
  global BoundaryCondition Model ModelFlags

  set Model(Emissivity,hasChanged) 0

  # We set display mode (temporarily) to edges/surfaces
  
  Util::setFlagValue SELECT_METHOD SELECT_METHOD_SINGLE 1

  set BoundaryCondition(DRAW_TARGET,prev) $ModelFlags(DRAW_TARGET)
  set BoundaryCondition(DRAW_SOURCE_CAD,prev) $ModelFlags(DRAW_SOURCE_CAD)
  set BoundaryCondition(DRAW_SOURCE_MESH,prev) $ModelFlags(DRAW_SOURCE_MESH)

  if { $Model(GEOMETRY_DIMENSION) == "2D" } {
    set BoundaryCondition(DRAW_TARGET,cur) DRAW_TARGET_EDGES
    Util::setFlagValue DRAW_TARGET DRAW_TARGET_EDGES 1

  } else {
    set BoundaryCondition(DRAW_TARGET,cur) DRAW_TARGET_SURFACES
    Util::setFlagValue DRAW_TARGET DRAW_TARGET_SURFACES 1
  }

  if { $ModelFlags(GEOMETRY_TYPE_CAD) } {
    set BoundaryCondition(DRAW_SOURCE_CAD,cur) 1
    set BoundaryCondition(DRAW_SOURCE_MESH,cur) 0
    Util::setFlagValue DRAW_SOURCE DRAW_SOURCE_CAD 1
    Util::setFlagValue DRAW_SOURCE DRAW_SOURCE_MESH 0
 }

}


#===========================#
#  Equation panel specific  #
#===========================#
#NOTE: Equation panel needs some special handling, fex.
#      parameters can be added although there is no object selected.


# Set equation to be shown in Equation-panel equation option menu
#
proc Equation::setOwnEquationIndices {} {
  global Equation

  set Equation(ownIndices) ""
  set Equation(ownNames) ""

  foreach fld $Equation(allFields) {

    # Is this an equation name?
    set eidx [DataField::getFieldProperty Equation $fld EquationIndex "" 0]

    if { $eidx == "" } {
      continue
    }

    # Should it be displayed in the Equation panel?
    if { 1 != [DataField::getFieldProperty Equation $fld IncludeInEquationPanel "" 0] } {
      continue
    }

    lappend Equation(ownNames) $fld
    lappend Equation(ownIndices) $eidx
  }

  # Check that active indices are legal
  if { [info exists Equation(activeIndices)] } {
    
    set tmp ""

    foreach idx $Equation(activeIndices) {

      if { "" != [Equation::getFieldPropertyByIndex $idx EquationField] } {
        lappend tmp $idx
      }
    }

    set Equation(activeIndices) $tmp
  }

}


proc Equation::postInitData {} {
  global Equation EquationVariable
  global Variables InitialCondition BoundaryCondition

  set Equation(currentEquationIndex) ""
  set Equation(updateValuesArea) 1

  #--We need this for handling Equation(id,mask) array
  set Equation(parameterLBList,old) $Equation(parameterLBList) 

  #--We need these to cancel an ok in the eq. variables panel!!!
  #  ===========================================================
  set Equation(equationVariables,ids) $EquationVariable(ids)

  foreach id $EquationVariable(ids) {
    set Equation(equationVariables,$id,data) $EquationVariable($id,data)
  }

  set Equation(variables,initialFields) $Variables(initialFields)
  set Equation(variables,allFields) $Variables(allFields)

  set Equation(initialCondition,initialFields) $InitialCondition(initialFields)
  set Equation(initialCondition,allFields) $InitialCondition(allFields)

  set Equation(boundaryCondition,initialFields) $BoundaryCondition(initialFields)
  set Equation(boundaryCondition,allFields) $BoundaryCondition(allFields)

  set Equation(boundaryCondition,initialDirichletVnames) $BoundaryCondition(initialDirichletVnames)
  set Equation(boundaryCondition,dirichletVnames) $BoundaryCondition(dirichletVnames)


  #--Init var values and system masks and body mask
  Equation::initMasks
}


proc Equation::initMasks {} {
  global Equation Model

  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    set Equation(dimgType) 千imG3
  } else {
    set Equation(dimgType) 千imG2
  }

  if { $Model(SIMULATION_DIMENSION) == "3D" } {
    set Equation(dimsType) 千imS3
  } else {
    set Equation(dimsType) 千imS2
  }

  set Equation(flow_type) ""
  set Equation(flow_laminarType) ""
  set Equation(flow_turbulentType) ""
  set Equation(flowType) ""

  set Equation(heat_type) ""
  set Equation(heat_transferType) ""
  set Equation(heat_phaseType) ""
  set Equation(heatType) ""

  set Equation(stress_type) ""
  set Equation(stress_mechanicalType) ""
  set Equation(stress_thermalType) ""
  set Equation(stressType) ""

  set Equation(diffusion_type) ""
  set Equation(diffusion_transferType) ""
  set Equation(diffusionType) ""

  set Equation(bcType) ""

  set Equation(userDefinedEquationType) ""

  Equation::initVariableTypeMasks
}


proc Equation::initVariableTypeMasks {} {
  global Equation

  foreach fld $Equation(allFields) {

    set eindex [DataField::getFieldProperty Equation $fld EquationIndex "" 0]
    
    if { $eindex == "" } {
      continue
    }

    set Equation($fld,variableType) ""
  }
}

proc Equation::setOtherMasks {} {
  global Equation Common Info Model

  # Dimension masks
  # ===============
  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    set dimgMask " (千imG3)"
  } else {
    set dimgMask " (千imG2)"
  }

  if { $Model(SIMULATION_DIMENSION) == "3D" } {
    set dimsMask " (千imS3)"
  } else {
    set dimsMask " (千imS2)"
  }

  # Boundary condition mask
  # =======================
  set bcMask "(刃cD) (刃cN2) (刃cN3)"
 
  set Equation(dimgType) $dimgMask
  set Equation(dimsType) $dimsMask
  set Equation(bcType) $bcMask
  
}


#========================#
#  GridH panel specific  #
#========================#
proc GridH::preInitData {} {
  global GridH Model
 
  GridH::updateObjectParameterId 
}


#================================#
#  GridParameter panel specific  #
#================================#
proc GridParameter::preInitData {} {

  GridParameter::updateObjectParameterId 
}


proc GridParameter::postInitData {} {

  GridParameter::checkPanel 
}


#-Set initial values for the fields in the argument
# panel array
#
proc StdPanelInit::initFieldValues {globArray} {
  global Common
  upvar #0 $globArray theArray

  set panel $theArray(parameterType)

  foreach fld $theArray(allFields) {

    # List boxes, buttons, labels etc. field are skipped here!
    if { ![Panel::isSimpleDataField $globArray  $fld] } {
      set theArray($fld,err) 0
      set theArray($fld,mod) 0
      #set theArray($fld,act) 0
      continue
    }

    #-Special sub widgets
    set sub_flds {use abs proc file}
    foreach sfld $sub_flds {
      if { [info exists theArray($fld,$sfld)] } {
        set theArray($fld,$sfld,act) 1
      }
    }

    set theArray($fld,variables) ""
    set theArray($fld,dataSize) ""

    #---Group field
    if { [Panel::isGroupField $globArray $fld] } {

      #-Radiobutton group data is stored in the parent, so
      # it needs an intial value. 
      if { [Panel::isRadioButtonField $globArray $fld] } {
        StdPanelInit::initFieldValue $globArray $fld
      }

      set glist $Common($panel,$fld,group)

      foreach g_fld $glist {

        if { ![Panel::isSimpleDataField $globArray  $g_fld] } {
          continue
        }

        StdPanelInit::initFieldValue $globArray $g_fld
      }

    #---Single field
    } else {
      StdPanelInit::initFieldValue $globArray $fld
    }

  }
}


# Set value for a single field
#
proc StdPanelInit::initFieldValue {globArray fld} {
  upvar #0 $globArray theArray

  # Bacic init
  Panel::initField $globArray $fld


  # If this is a group main widget, check its state
  if { [Panel::isGroupField $globArray $fld]  &&
       [info exists theArray(allWidgets,$fld)] &&
       [winfo exists $theArray(allWidgets,$fld)] 
     } {

    # If this widget type has a "state"
    set woptions [$theArray(allWidgets,$fld) configure]
    if { -1 != [lsearch $woptions "state"] &&
         "disabled" == [$theArray(allWidgets,$fld) cget -state]
       } {
        set theArray($fld) ""

    } else {
      set theArray($fld) [DataField::getInitialValue $globArray $fld]
    }

  # Otherwise just use the initial value
  } else {
    set theArray($fld) [DataField::getInitialValue $globArray $fld]
  }

  set procedurable [DataField::getFieldProperty $globArray $fld Procedurable]
  if {!$procedurable} {
    set theArray($fld,proc) 0
  }

  set tableable [DataField::getFieldProperty $globArray $fld Tableable]
  if {!$tableable} {
    set theArray($fld,table) 0
  }

}

# end ecif_tk_panelInit.tcl
# ********************
