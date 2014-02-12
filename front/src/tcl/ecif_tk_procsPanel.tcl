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
#Module:    ecif_tk_procsPanel.tcl
#Language:  Tcl
#Date:      29.06.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Generic proc for panels
#
#************************************************************************


# Get panel (array's) menu name (used mainly for
# user friendly messages!)
#
proc Panel::panelNameGuiToMenu {globArray} {
  global Common
  upvar #0 $globArray theArray

  if { [info exist theArray(menuName)] } {
    return $theArray(menuName)

  } else {
    return $globArray
  }
}


# =====================
# Fields handling stuff
# =====================

# Backup fields
# =============


# Store old values and reset states for argument field
#
proc Panel::backupField {globArray fld} {
  upvar #0 $globArray theArray

  set theArray($fld,old) $theArray($fld)
  set theArray($fld,prev) $theArray($fld)

  set theArray($fld,err) 0
  set theArray($fld,mod) 0
}

# Store old values and reset states for all case fields
#
proc Panel::backupFields {globArray {fields ""} } {
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {

    # Backup only current target fields
    if { ![info exists theArray($fld,trg)] || !$theArray($fld,trg) } {
      continue
    }
    
    # Backup only output fields, other are just
    # reset to normal status
    #
    if { ![Panel::isOutputField $globArray $fld] } {
      set theArray($fld,err) 0
      set theArray($fld,mod) 0
      continue
    }

    Panel::backupField $globArray $fld

    if { [info exists theArray($fld,file)] } {
      Panel::backupField $globArray $fld,file
    }

    if { [info exists theArray($fld,proc)] } {
      Panel::backupField $globArray $fld,proc
    }

    if { [info exists theArray($fld,table)] } {
      Panel::backupField $globArray $fld,table
    }
      
  }
}


# Reset fields 
# ============

# Reset one field
#
proc Panel::resetField {globArray fld} {
  global $globArray
  upvar #0 $globArray theArray

  set theArray($fld,err) 0
  set theArray($fld,mod) 0

  set theArray($fld,proc) 0
  set theArray($fld,table) 0

  # Remove all possible data value type fields
  catch { unset theArray($fld,dataSize) }
  catch { unset theArray($fld,variables) }
  catch { unset theArray($fld,nofEntries) }
  catch { unset theArray($fld,valueData) }
  catch { unset theArray($fld,procData) }
  catch { unset theArray($fld,procData,dataSize) }
  catch { unset theArray($fld,procData,variables) }
  catch { unset theArray($fld,procData,nofEntries) }
  catch { unset theArray($fld,tableData) }
  catch { unset theArray($fld,tableData,dataSize) }
  catch { unset theArray($fld,tableData,variables) }
  catch { unset theArray($fld,tableData,nofEntries) }

  # Value-data is the "default" data type
  #
  set theArray($fld) ""
  set theArray($fld,valueData) ""
  set theArray($fld,firstInstance) 1
}


# Reset list of fields
#
proc Panel::resetFields {globArray {fields ""} } {
	global Common
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {

    Panel::resetField $globArray $fld

    # Init also group members!
    if { [Panel::isGroupField $globArray $fld] } {

      foreach gfld $Common($theArray(parameterType),$fld,group) {
        Panel::resetField $globArray $gfld
      }
    }
  }
}


# Init fields
# ===========


# Init single field
#
proc Panel::initField {globArray fld {force_initial 0} {system_index ""} } {
  global Common Model
  upvar #0 $globArray theArray

  set always_initialize [DataField::getFieldProperty $globArray $fld AlwaysInitialize $system_index]

  if { $force_initial || $always_initialize } {
    Panel::resetField $globArray $fld
  }

  set prb_mask [Panel::getProblemMask $globArray]
  set trg_mask [Panel::getTargetMask $globArray]

  set def_dsp [DataField::getFieldProperty $globArray $fld Display $system_index]

  #-No data entry if no model defined (except UserSetting)!
  if { $Model(GEOMETRY_DIMENSION) == "" &&
       $globArray != $Common(panelArray,US)
     } {
    set def_act 0

  #-Otherwise field property value decides
  } else {
    set def_act [DataField::getFieldProperty $globArray $fld InitiallyActive $system_index]
  }

  # NOTE: 'Positive' attitude here. Panels should
  # turn off these flags when needed!
  #
  set theArray($fld,prb) 1 ;# Belongs to the problem
  set theArray($fld,trg) 1 ;# Belongs to the current case
  set theArray($fld,out) 1 ;# Is output field

  # Field is displayed (packed) by default, if not
  # explicitely defined something else
  if { !$def_dsp } {
    set theArray($fld,dsp) 0
    set theArray($fld,scr) 0
  } else {
    set theArray($fld,dsp) 1
    set theArray($fld,scr) 1
  }

  # NOTE: Field is active (def_act=1), if not
  # explicitely defined something else
  #
  if { !$def_act } {
    set theArray($fld,act) 0
  } else {
    set theArray($fld,act) 1
  }

  set theArray($fld,err) 0 ;# Error status
  set theArray($fld,mod) 0 ;# Modified status

  if {1 == [DataField::getFieldProperty $globArray $fld HasAbsCheckBox "" 0]} {
    set theArray($fld,Abs,err) 0
    set theArray($fld,Abs,mod) 0
  }

  if {1 == [DataField::getFieldProperty $globArray $fld HasUseCheckBox "" 0]} {
    set theArray($fld,Use,err) 0
    set theArray($fld,Use,mod) 0
  }

  #--If field does not belong to the problem
  if { $prb_mask != "" &&
       ![DataField::fieldProblemMaskMatches $globArray $fld $prb_mask]
     } {
    set theArray($fld,prb) 0
    set theArray($fld,dsp) 0
    set theArray($fld,trg) 0
    set theArray($fld,scr) 0
    set theArray($fld,act) 0
  }

  #--If field does not belong to current case
  if { $theArray($fld,prb) &&
       $trg_mask != "" &&
       ![DataField::fieldTargetMaskMatches $globArray $fld $trg_mask]
     } {
    set theArray($fld,trg) 0
    set theArray($fld,act) 0
  }

  #--Field is not a screen field
  if { $theArray($fld,prb) &&
       ![Panel::isScreenField $globArray $fld]
     } {
    set theArray($fld,scr) 0
    set theArray($fld,dsp) 0
  }

  # If any initializing allowed
  #
  if { $force_initial != -1 } {

    # Pick proper default value
    if { [info exist theArray($fld,default)] &&
         $theArray($fld,default) != ""
       } {
      set def_value $theArray($fld,default)

    } else {
      set def_value [DataField::getInitialValue $globArray $fld $system_index]
    }
    
    set drop_initial [DataField::getFieldProperty $globArray $fld DropInitialValued]

    # Initialize field value
    #
    # NOTE: Do not initialize a non-forced field if would be dropped
    # as initial valued. Otherwise the model file would always seem
    # to changed when loaded!!!
    #
    if { $force_initial ||
         $always_initialize ||
         ( !$drop_initial && 
           (![info exists theArray($fld)] || $theArray($fld) == "")
         )
       } {
      set theArray($fld) $def_value
      set theArray($fld,valueData) $def_value
    }
  }

  # If field not defined, stop!
  #  
  if { ![info exists theArray($fld)] } {
    return
  }

  set theArray($fld,old) $theArray($fld)
  set theArray($fld,prev) $theArray($fld)
}


# Init fields
#
proc Panel::initFields {globArray {fields ""} {force_initial 0} {system_index ""} } {
  global Common
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {

    Panel::initField $globArray $fld $force_initial $system_index

    # Init also group members!
    if { [Panel::isGroupField $globArray $fld] } {

      foreach gfld $Common($theArray(parameterType),$fld,group) {
        Panel::initField $globArray $gfld $force_initial $system_index
      }
    }

  }
}


# Init values area fields
#
proc Panel::initValuesAreaFields {globArray {fields ""} {force_initial 0} {system_index ""} } {
  global Common
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  # We skip frame area 0 here, for Equation also
  # the equalion list frmame (1)
  #
  if { $globArray == $Common(panelArray,EQ) } {
    set min_frame 2
  } else {
    set min_frame 1
  }

  foreach fld $fields {

    set fld_frame [DataField::getFieldProperty $globArray $fld SubPanel $system_index]

    if { $fld_frame < $min_frame } {
      continue
    }

    Panel::initField $globArray $fld $force_initial $system_index

    # Init also group members!
    if { [Panel::isGroupField $globArray $fld] } {

      foreach gfld $Common($theArray(parameterType),$fld,group) {

        Panel::initField $globArray $gfld $force_initial $system_index
      }
    }

  }
}


# Init new (added) fields in the panel
#
proc Panel::initNewPanelFields {globArray  {force_initial 0} {system_index ""} } {
  global $globArray
  upvar #0 $globArray theArray

  # Collect panel fields
  set fields ""
  foreach fld $theArray(allFields) {

    # Init only new fields
    if { ![info exist theArray($fld,trg)] } {
      lappend fields $fld
    }
  }

  Panel::initFields $globArray $fields $force_initial $system_index
}


# Init fields in the currently visible panel (pane)
#
proc Panel::initCurrentPanelFields {globArray  {force_initial 0} {system_index ""} } {
  global $globArray
  upvar #0 $globArray theArray

  # Collect panel fields
  set fields ""
  foreach fld $theArray(allFields) {

    # Init all problem fields
    if { $theArray($fld,trg) } {
      lappend fields $fld
    }
  }

  Panel::initFields $globArray $fields $force_initial $system_index
}


# Init all panel fields (for all panes)
#
proc Panel::initAllPanelFields {globArray  {force_initial 0} {system_index ""} } {
  global $globArray
  upvar #0 $globArray theArray

  # Collect panel fields
  set fields ""
  foreach fld $theArray(allFields) {

    # Init all problem fields
    if { $theArray($fld,prb) } {
      lappend fields $fld
    }
  }

  Panel::initFields $globArray $fields $force_initial $system_index
}


# Restore fields
# ==============

# Restore current target fields
#
proc Panel::restoreFields {globArray {fields ""} } {
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {

    if { !$theArray($fld,trg) } {
      continue
    }

    set theArray($fld) $theArray($fld,old)
  }
}


# Unset array fields
# =====================
#
proc Panel::unsetFields {globArray {fields ""} } {
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {

    Panel::unsetField $globArray "$fld"
    Panel::unsetField $globArray "$fld,old"
    Panel::unsetField $globArray "$fld,prev"

    Panel::unsetField $globArray "$fld,file"
    Panel::unsetField $globArray "$fld,file,old"
    Panel::unsetField $globArray "$fld,file,prev"

    Panel::unsetField $globArray "$fld,procedure"
    Panel::unsetField $globArray "$fld,procedure,old"
    Panel::unsetField $globArray "$fld,procedure,prev"

    Panel::unsetField $globArray "$fld,prb"
    Panel::unsetField $globArray "$fld,trg"
    Panel::unsetField $globArray "$fld,dsp"
    Panel::unsetField $globArray "$fld,scr"
    Panel::unsetField $globArray "$fld,act"
    Panel::unsetField $globArray "$fld,out"

    Panel::unsetField $globArray "$fld,err"
    Panel::unsetField $globArray "$fld,mod"

    Panel::unsetField $globArray "$fld,dataSize" 
    Panel::unsetField $globArray "$fld,valueData" 
    Panel::unsetField $globArray "$fld,procData"
    Panel::unsetField $globArray "$fld,tableData"
    Panel::unsetField $globArray "$fld,variables" 
    Panel::unsetField $globArray "$fld,wframe"
  }
}


proc Panel::unsetAllFieldData {globArray {fields ""} } {
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {
    Panel::unsetField $globArray $fld
    Panel::unsetField $globArray "$fld,*"
  }
}


# Unset one array field
# (remove all variables related to array variable "fld")
#
proc Panel::unsetField {globArray fld } {
  upvar #0 $globArray theArray

  Util::unsetArrayVariables $globArray $fld
}


# Set field states (mod etc.)
#
proc Panel::setFieldStates {globArray {fields ""}} {
  upvar #0 $globArray theArray

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {
    Widget::configureField $globArray $fld
  }
}



# Field select procs
####################


# Set proper fields to be displayed for the
# current problem (system) on the panel. 
# Problem mask is used as the criterium.
#
# NOTE: Display of the field is controlled by the globArray(fld,dsp)
# variable which is updated here!!!
#
proc Panel::selectDisplayFields {globArray} {
  global Equation Info
  upvar #0 $globArray theArray

  if { $theArray(panelsByEquation) &&
       $theArray(hasMaskedFields)
     } {

    if { ![info exists theArray(equationIndex)] ||
         $theArray(equationIndex) == "" ||
         $theArray(equationIndex) < 0
       } {
      return
    }

    set theArray(displayMask) [Equation::getEquationMask $theArray(equationIndex)]

  } else {
    set theArray(displayMask) ""
  }

  #--Select fields which match display mask
  set pnames ""

  foreach fld $theArray(allFields) {

    set theArray($fld,dsp) 0

    # Only problem screen fields can be
    # display fields!
    if { !$theArray($fld,prb) ||
         !$theArray($fld,scr) 
       } {
      continue
    }

    # Skip index variables like DIFFUSION (only their derivatives
    # like <Oxygen>DIFFUSION are included!)
    #
    if { "" != [DataField::getFieldProperty $globArray $fld IndexVariable] } {
      continue
    }
    
    # Display field
    if { [DataField::fieldDisplayMaskMatches $globArray $fld $theArray(displayMask)] } {
      set theArray($fld,dsp) 1
    }

  }
}


# Select fields for the array starting from the given
# source-fields list  and the given problem-mask
#
proc Panel::selectFields {globArray source_fields {problem_mask ""} } {
  global Info
  upvar #0 $globArray theArray

  set fnames ""

  foreach fld $source_fields {

    #--Should we pick the field
    if { $problem_mask != "" &&
         ![DataField::fieldProblemMaskMatches $globArray $fld $problem_mask]
       } {
      continue
    }

    #--Is it an active (normal, non-disabled) field
    if { $problem_mask != "" &&
         ![DataField::fieldTargetMaskMatches $globArray $fld $problem_mask]
       } {
      set theArray($fld,act) 0
    } else {
      set theArray($fld,act) 1
    }

    #--Add field to the result list
    lappend fnames $fld
#MSG "Array=$globArray fld=$fld"
  }

  return $fnames
}


# Select variables for the case using the given mask as the criterium
# NOTE: Inactive fields are also selected, if Info(keepInactiveFields)
# flag is on
#
proc Panel::selectFields2 {globArray fld_names problem_mask} {
  global Info
  upvar #0 $globArray theArray
 
  set vnames ""

  foreach fld $fld_names {

    if { $problemMask != "" &&
         ![DataField::fieldProblemMaskMatches $globArray $fld $problem_mask]
       } {
      set active  0
    } else {
      set active 1
    }

    set theArray($fld,act) $active

    if { [Panel::isGroupField $globArray $fld] } {
      Panel::setGroupMembersValue $globArray $fld "active" $active
    }

    # Keep also inactive fields if the flag is on
    #
    if { $theArray($fld,act) || $Info(keepInactiveFields) } {
      lappend vnames $fld
    }
  }

  return $vnames
}


proc Panel::getProblemMask {globArray} {
  global Common Equation Info

  #-Equation panel accepts all defined fields, otherwise
  # we could not add new equations!
  if { $globArray == $Common(panelArray,EQ) } {
    #return $Equation(totalMask)
    return ""

  #-In other panels we are more restrictive

  # GridParameter specific
  } elseif { $globArray == $Common(panelArray,GR) } {
    return $Equation(dimgType)

  # GridH (density) specific
  } elseif { $globArray == $Common(panelArray,GH) } {
    return $Equation(dimgType)

  # Solver
  } elseif { $globArray == $Common(panelArray,SL) } {
      return $Equation(problemMask)
    
  # Other arries
  } else {

    # These arries follow equation mask
    if { 1 == [Util::getArrayValue $globArray hasMaskedTarget] } {
      #return $Equation(problemMask)
      return ""
    
    # These arries could have some array specific value defined, but
    # normally the default value blank is returned here!
    } else {
      return [Util::getArrayValue $globArray problemMask]
    }

  }
}


proc Panel::getTargetMask {globArray} {
  global Common Equation Info
  upvar #0 $globArray theArray

  #-In equation panel we accept all defined fields, otherwise
  # we could not add new equations!
  if { $globArray == $Common(panelArray,EQ) } {
    #return $Equation(problemMask)
    return $Equation(targetMask)

  # GridParameter
  } elseif { $globArray == $Common(panelArray,GR) } {
    return $Equation(dimgType)

  # GridH
  } elseif { $globArray == $Common(panelArray,GH) } {
    return $Equation(dimgType)
    
  #-In other panels we are more restrictive
  } else {	
    return [Util::getArrayValue $globArray targetMask]
  }
}


proc Panel::getSystemIndex {globArray pid} {
  global Common SolverSystem
  upvar #0 $globArray theArray

  if { $globArray != $Common(panelArray,SL) } {
    return ""
  }

  # Find system index by equation name
  #
  set eq_name [DataField::getFieldValue Solver $pid EQUATION]

  foreach idx $SolverSystem(indices) {
    
    if { [string equal -nocase $eq_name $SolverSystem($idx,outName)] } {
      return $idx
    }
  }

  # Nothing found, Error for Solver!
  if { $globArray == $Common(panelArray,SL) } {
    return -1
  
  # Other panels,ok
  } else {
    return ""
  }
}


# =========================
# Parameters handling stuff
# =========================

# Check all input parameters in all arries given in the argument
# Remove non-input fields from each parameter
#
proc Panel::checkInputParameters { array_list vnames_name} {

  # Loop argument arries
  # ====================
  foreach globArray $array_list {
    upvar #0 $globArray theArray

    if { ![info exists theArray(ids)] ||
         $theArray(ids) == ""
       } {
      continue
    }

    # Loop all parameter entries
    # ==========================
    foreach id $theArray(ids) {
      
      if { ![info exists theArray($id,data)] } continue

      set new_data [DataField::compressParameterInputData $globArray $id $vnames_name]
      set theArray($id,data) $new_data
    }
  }
}


# Check all masked parameters
#
proc Panel::checkMaskedParameters {modified_flag {save_modified_data 1} } {
  global Common
  upvar $modified_flag modified
    
  Panel::checkParameters $Common(maskedParameters) modified $save_modified_data
}
 

# Check all parameters in arries given in the array_list argument
# Argument "data_modified" is a reference variable to indicate if data was modifeid
#
# NOTE: Normally modified data will be saved, but if this proc is called ex.
# during Equation::ok processing for checking Equation-parameters themselves, it
# may  not make sense to save here, because Equation will save its own stuff
# anyway
#
proc Panel::checkParameters { array_list data_modified {save_modified_data 1} { force_initial 0} } {
  global Common Info Equation ObjectTable
  upvar $data_modified modified

  set modified 0
  set detached 0

  # Loop argument arries
  # ====================
  foreach globArray $array_list {
    upvar #0 $globArray theArray

    if { ![info exists theArray(ids)] ||
         $theArray(ids) == ""
       } {
      continue
    }
#MSG "array=$globArray"
    # Loop all parameter entries
    # ==========================
    set new_param_list ""
    set new_problem_list ""

    set amodified 0

    foreach pid $theArray(ids) {

      set new_mask ""

      #-Pick mask for the parameters from parent objects
      if { $theArray(hasMaskedParameter) } {
        set oid $theArray($pid,oid)      
        set new_mask [StdPanelExec::getTargetMask $globArray $oid]

      #-Or just use current mask (if exists)
      } elseif { [info exists theArray($pid,mask)] } {
        set new_mask $theArray($pid,mask)
      }

      #-Form checked parameter data line
      set pmodified 0

      set system_index [Panel::getSystemIndex $globArray $pid]

      #set new_data $theArray($pid,data)

      # If something wrong with data, remove it!
      #
	    if { $system_index == -1 } {
		    set new_data ""
        Message::showMessage "WARNING! Incorrect parameter definition ($globArray $pid) removed!" $Info(wrnMsgColor)
      
      # Data Ok
      # NOTE: DataField::checkParameterData takes time!!!
      #
	    } else {
        #set new_data $theArray($pid,data)
		    set new_data [DataField::checkParameterData $globArray $pid pmodified $new_mask $system_index $force_initial]
	    }

      if { $pmodified } {
        set amodified 1
      }

      # If result was an empty parameter, delete it
      if { $new_data == "" } {
        StdPanelExec::deleteParameterById $globArray $pid
        set dmodified 1

      # Otherwise update (possibly) modified parameter
      } else {

  			set dmodified 0
        set theArray($pid,data) $new_data

        # NOTE: Do not set Equation mask here! It is NOT
        # concluded from body masks!!!
        #
        if { $theArray(hasMaskedParameter) } {
          set theArray($pid,mask) $new_mask
          StdPanelExec::detachInconsistentAttachments $globArray $pid dmodified
        }
      }
      
      if { $dmodified } {
        set amodified 1
        set detached 1
      }
    } ;# For each pid

    # Update modified array
    # =====================
    if { $amodified && $save_modified_data } {
      set modified 1
      StdPanelExec::compressParameterIds $globArray

      StdPanelExec::panelSave $globArray
    }
  } ; # For each array

  if {$detached} {
    set msg [list "NOTE: Equation definition(s) were changed and inconsistent parameter\n" \
                  "attachments were detached. Check warnings in the message area and\n" \
                  "and correct parameter attachments (Model menu)!"]

    set Info(messageIcon) warning
    Message::dispOkMessage $msg
  }

  # Make current parameter selected if panel was visible
  if { [Util::panelExists $globArray] } {
    if { [info exists theArray(parameterId)] } {
      StdPanelExec::parameterSelected $globArray
    }
  }

}


# Make parameter ids consequtive ie. 1 2 3 ...
# Relevant when parameters have been deleted!

# NOTE: This is different from corresponding
# proc in StdPanelExec:: !!!
# These parameters do not have targets like bodies or 
# boundaries
# Used for: Solver, Timestep etc.
#
proc Panel::compressParameterIds {globArray} {
  global Info
  upvar #0 $globArray theArray

  # Update also current parameter id!
  if { [info exists theArray(parameterId)] } {
    set old_current_id $theArray(parameterId)
  } else {
    set old_current_id $Info(NO_INDEX)
  }

  set new_current_id $Info(NO_INDEX)

  # NOTE: Sort as integers, This is important, because we delete
  # the old-id data if id-number was "compressed", and for this
  # it is important that we go through the old-ids in increading
  # order!!!
  #
  set theArray(ids) [lsort -integer $theArray(ids)]
  set new_ids ""
 
  # Form updated parameter data list where parameter
  # ids are replaced with the "compressed" ids
  #
  set counter 1
  foreach old_id $theArray(ids) {
  
    # New id
    set new_id $counter
    
    if { $new_id != $old_id } {

      set theArray($new_id,data) $theArray($old_id,data)
      set theArray($new_id,name) $theArray($old_id,name)
      
      #--Check if this was current
      if { $old_current_id == $old_id } {
        set new_current_id $new_id
      }

      #--Check if old name was an "automatic" name, if so, it will
      #  be also "renumbered"!!!
      #
      set old_def_name [Panel::defaultParameterName $globArray $old_id]
      set new_def_name [Panel::defaultParameterName $globArray $new_id]

      if { [string equal $theArray($old_id,name) $old_def_name] } {
        set theArray($new_id,name) $new_def_name
      }

      # Remove old-id vars

      Util::unsetArrayIdVariables $globArray $old_id
    }

    lappend new_ids $new_id

    incr counter
  }

  set theArray(ids) $new_ids
  set theArray(nextNewParameterId) $counter
  
  if { $new_current_id != $Info(NO_INDEX) } {
    set theArray(parameterId) $new_current_id
  }
}   


# Remove non-output fields from all parameters in the array before saving to cpp
# NOTE: Do not loose (...,uncompressed) data before using uncompress proc!!!
#
proc Panel::compressParameters {globArray modified_flag} {
  upvar #0 $globArray theArray
  upvar $modified_flag modified

  #-In cpp we use compressed (excluding empty fields etc.) data
  foreach id $theArray(ids) {
    set amodified 0
  
    set tmp [DataField::compressParameterOutputData $globArray $id amodified]
#print "comp=$tmp"
    set theArray($id,data,uncompressed) $theArray($id,data)
    set theArray($id,data) $tmp

    if {$amodified} {
      set modified 1
    }
  }
}


# Restore all fields for all parameters in the array after saving to cpp
# NOTE: Do not loose (...,uncompressed) data before using this proc!!!
#
proc Panel::uncompressParameters {globArray} {
  upvar #0 $globArray theArray

  #-In workspace we use uncompressed (including empty fields etc.) data
  foreach id $theArray(ids) {
    set theArray($id,data) $theArray($id,data,uncompressed)
    unset theArray($id,data,uncompressed)
  }
}


proc Panel::createDefaultParameterLine {globArray {system_index ""} } {
  global $globArray Info
  upvar #0 $globArray theArray

  #-Data separator
  set sep $Info(dispListSeparatorOut)

  Panel::initFields $globArray "" 1 $system_index

  set result ""  

  # Create field values for the current target
  #
  foreach fld $theArray(allFields) {

    if { !$theArray($fld,trg) } {
      continue
    }

    set tmp [DataField::formInitialParameterLineValue $globArray $fld $system_index]

    if { $tmp == "" } {
      continue
    }

    append result $tmp$sep
  }

  set result [string trimright $result $sep]

  return $result
}





# Set array variable value for field's group members
#
proc Panel::setGroupMembersValue {globArray fld variable_name value} {
  global Common
  upvar #0 $globArray theArray

  if { ![info exists Common($theArray(parameterType),$fld,group)] } {
    return
  }

  foreach gv $Common($theArray(parameterType),$fld,group) {

    set theArray($gv,$variable_name) $value
  }
}


# Return needed-level:
# 0 <--> missing value: accepted
# 1 <--> missing value: warning
# 2 <--> missing value: error
#
proc Panel::getFieldNeededLevel {globArray fld} {
  upvar #0 $globArray theArray

  set ndl [DataField::getFieldProperty $globArray $fld NeededLevel]

  # No NeededLevel defined
  # Missing value accepted
  if { $ndl == "" || $ndl == 0 } {
    return 0
  }

  # If field's NeededLevel value is <= 2 and a parameter file given
  # Missing value accepted
  if { $ndl <= 2 && [Panel::isActiveField $globArray INCLUDE] } {
    return 0

  # Missing value ==> warning!
  } elseif { $ndl == 1 } {
    return 1

  # Missing value ==> error!
  } else {
    return 2
  }
}



# Construct a list pair of the currently matching
# (by object mask) problem variables
# 1.pair  list of internal names like "COORDINATE_1"
# 2.pair  list of external names like "Coordinate 1"
#
proc Panel::getCurrentVariables {object_mask {accept_array 0} } {
  global Common Model Timestep Variables

  set var_list1 ""
  set var_list2 ""

  set dim $Model(SIMULATION_DIMENSION)
  
  if { $Timestep(SIMULATION_TYPE) == "Transient" } {
    #set object_mask [Util::combineMasks $object_mask $Common(timeMask)]
    set object_mask [Util::combineMasks $object_mask ¤T]
  }

  foreach fld $Variables(allFields) {

    # Pick only variables for the current case
		#
    if { !$Variables($fld,trg) } {
      continue
    }

    if { $object_mask != "" &&
         ![DataField::fieldTargetMaskMatches Variables $fld $object_mask] } {
      continue
    }

		# Pick only variables who have proper dimension
		# NOTE: This is used to drop non-scalar variables from table entries!
		#
		# Check if field is an array
		if { "Array" == [DataField::getFieldProperty Variables $fld FieldDataType] } {
			set is_array 1
		} else {
			set is_array 0
		}

		if { $is_array && !$accept_array } {
			continue
		}

    set name [DataField::getSifName Variables $fld]

    # Note: Not in use!
    set coord_index [DataField::getFieldProperty Variables $fld CoordinateIndex]

    if { $coord_index > 0 } {

      if { $coord_index == 3 && $dim == "2D" } {
        ;#continue
      } else {
        ;#append name " $coord_index"
      }
    }
    #

    lappend var_list1 $fld
    lappend var_list2 $name
  }

  return [list $var_list1 $var_list2]
}



# Create an option menu button
#
proc Panel::createOptionMenuWidget {widget fld values } {

  append om "tk_optionMenu $widget $fld "

  # Make a 'flat' list of values: "A" "B" "C" ...
  #
  set values [join $values "\" \""]
  set values "\"$values\""
  append om $values

  # Create widget
  return [eval $om]
}


# CLEAR ENTRY button
#
proc Panel::clearEntryButtonCommandProc {globArray} {
  global Info
  upvar #0 $globArray theArray

  set fld $theArray(currentEntryField)

  if { $fld == "" } {
    return
  }

  set label [string toupper [DataField::getFieldProperty $globArray $fld Label]]

  set msg "Are you sure to delete the contents of the $label entry!"

  if { $theArray($fld) != "" && ![Message::verifyAction $Info(noviceUser) $msg question] } {
    return
  }

  # Current is proc
  if { $theArray($fld,proc) } {
    set theArray($fld,proc) 0
    set theArray($fld,procData) ""
 
  # Current is table
  } elseif { $theArray($fld,table) } {
    set theArray($fld,table) 0
    set theArray($fld,tableData) ""
  
  # Current is scalar (value)
  } else {
    set theArray($fld,valueData) ""
  }

  # Next is always the scalar (value)!
  set theArray($fld) $theArray($fld,valueData)

  # Check need of an initial value
  set force_initial 0
  Panel::initField $globArray $fld $force_initial

}


# CLEAR PANEL button
# Clears all fields in current panel
# E.g. all Heat Equation initial condition fields
#
proc Panel::clearPanelButtonCommandProc {globArray} {
  global Info
  upvar #0 $globArray theArray

  set msg "Are you sure to delete the contents of all entry fields in the panel!"

  if { ![Message::verifyAction $Info(noviceUser) $msg ok question] } {
    return
  }

  Panel::initCurrentPanelFields $globArray 1
}


# CLEAR PANELS button
# Clears all fields in all panes
# E.g. Navier-Stokes, Heat Equation etc. initial condition fields
#
proc Panel::clearPanelsButtonCommandProc {globArray} {
  global Info
  upvar #0 $globArray theArray

  set msg "Are you sure to delete the contents of all entries in all $globArray panels!"

  if { ![Message::verifyAction $Info(advancedUser) $msg ok question] } {
    return
  }

  Panel::initAllPanelFields $globArray 1

  if { $theArray(hasParameters) } {
    set theArray(parameterName) ""
  }
}


# EDIT button
#
proc Panel::editButtonCommandProc {globArray} {
  upvar #0 $globArray theArray

  set fld $theArray(currentEntryField)

  if { $fld == "" } {
    return
  }

  set wdg $theArray(allWidgets,$fld)

  Widget::entryDouble-1 $globArray $fld $wdg
}


# ================================
# Table and Procedure button stuff
# ================================

proc Panel::procedureCheckBoxCommandProc {globArray} {
  global Info
  upvar #0 $globArray theArray

  if { $theArray(currentEntryField) != "" } {
    set fld $theArray(currentEntryField)
  } else {
    return
  }

  # Store normal data, if it was active
  if { $theArray($fld,proc) && 
       (![info exists theArray($fld,table)] || !$theArray($fld,table))
     } {
    set theArray($fld,valueData) $theArray($fld)
  }

  # Proc checked on
  if { $theArray($fld,proc) } {
    set theArray($fld,table) 0

     if { [info exist theArray($fld,procData)] } {
      set theArray($fld) $theArray($fld,procData)
      set theArray($fld,dataSize) $theArray($fld,procData,dataSize)
      set theArray($fld,variables) $theArray($fld,procData,variables)

    } else {
      set theArray($fld) ""
      catch { unset theArray($fld,dataSize) }
      catch { unset theArray($fld,variables) }
      Widget::entryDouble-1 $globArray $fld $theArray(allWidgets,$fld)
    }


  # Proc checked off
  } else {

    if { [info exist theArray($fld,valueData)] } {
      set theArray($fld) $theArray($fld,valueData)

    } else {
      set theArray($fld) ""
    }

    catch { unset theArray($fld,dataSize) }
    catch { unset theArray($fld,variables) }
  }
}


proc Panel::tableCheckBoxCommandProc {globArray} {
  global Info
  upvar #0 $globArray theArray

  if { $theArray(currentEntryField) != "" } {
    set fld $theArray(currentEntryField)
  } else {
    return
  }

  # Store normal data, if it was active
  if { $theArray($fld,table) && 
       (![info exists theArray($fld,proc)] || !$theArray($fld,proc))
     } {
    set theArray($fld,valueData) $theArray($fld)
  }

  # Table checked on
  if { $theArray($fld,table) } {
    set theArray($fld,proc) 0

    if { [info exist theArray($fld,tableData)] } {
      set theArray($fld) $theArray($fld,tableData)
      set theArray($fld,dataSize) $theArray($fld,tableData,dataSize)
      set theArray($fld,variables) $theArray($fld,tableData,variables)

    } else {
      set theArray($fld) ""

      catch { unset theArray($fld,dataSize) }
      catch { unset theArray($fld,variables) }

      Widget::entryDouble-1 $globArray $fld $theArray(allWidgets,$fld)
    }

  # Table checked off
  } else {

    if { [info exist theArray($fld,valueData)] } {
      set theArray($fld) $theArray($fld,valueData)

    } else {
      set theArray($fld) ""
    }

    catch { unset theArray($fld,dataSize) }
    catch { unset theArray($fld,variables) }
  }
}


# Used eg. when parameter name field is entered
#
proc Panel::setProcAndTableButtonStates {globArray proc_state table_state } {
  upvar #0 $globArray theArray

  if { $theArray(hasProcButton) } {
    Widget::configureS $theArray(allWidgets,PROCEDURE_BUTTON) $proc_state
  }

  if { $theArray(hasTableButton) } {
    Widget::configureS $theArray(allWidgets,TABLE_BUTTON) $table_state
  }
}

proc Panel::editArray {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  #-Create qualified name for the result
  # NOTE This is needed because dispDataEntry proc
  # uses button commands which are in their own scope!
  namespace eval MyNs {
    variable result ""
  }

  Widget::dispDataEntry ::MyNs::result $theArray($fld) "Entry data window" "Enter values"
  tkwait variable ::MyNs::result

  set entry_result $::MyNs::result
 
  #-If cancel button was pressed
  if { $entry_result == "!cancel!" } {
    return
  }

  #-Set new value
  $wdg delete 0 end
  $wdg insert 0 $entry_result
}



# Proc creates a data entry panel for an entry field
#
proc Panel::editTable {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  set theArray(isActive) 0
  set theArray(hasFocus) 0

  #-Create a data entry widget
  set de_id [TableEntryPanel::create $globArray $fld Panel::tableEntryPanelCallback]
  
  # Store entry panel id
  lappend theArray(tableEntryPanelIds) $de_id

}


# Proc is the entry field update callback for a data entry panel
#
proc Panel::tableEntryPanelCallback {de_id globArray fld} {
  global Info
  upvar #0 $globArray theArray

  set entry_result [set TableEntryPanel::result$de_id]
  set entry_list [split $entry_result $Info(dataListSeparator)]
   
  if { 1 == [set TableEntryPanel::updated$de_id] } {
    Panel::panelDataModified 1 $globArray ""
  }

  #--Update field data
  #
  set sz1 [set TableEntryPanel::dataSize1$de_id]
  set sz2 [set TableEntryPanel::dataSize2$de_id]
  set ez  [set TableEntryPanel::nofEntries$de_id]

  set vn  [set TableEntryPanel::variableName$de_id]

  if { $sz1 == "N" } {
    set sz1 [llength [lindex $entry_list 0]]
  }

  if { $sz2 == "N" } {
    set sz2 [llength $entry_list]
  }

  if { $vn == "none" } {
    set vn ""
  }  

  set theArray($fld) $entry_result
  set theArray($fld,tableData) $entry_result

  set theArray($fld,dataSize) [list $sz1 $sz2]
  set theArray($fld,tableData,dataSize) [list $sz1 $sz2]

  set theArray($fld,nofEntries) $ez
  set theArray($fld,tableData,nofEntries) $ez

  set theArray($fld,variables) $vn
  set theArray($fld,tableData,variables) $vn

  #--Set new values to entry
  set wdg $theArray(allWidgets,$fld)
  if {[winfo exists $wdg]} {
    $wdg delete 0 end
    $wdg insert 0 $entry_result
  }

  #--Clear the data entry instance
  TableEntryPanel::delete $de_id

  #--Remove entry panel id from the list of active panels
  set tmp ""
  foreach id $theArray(tableEntryPanelIds) {
    if { $id != $de_id } {
      lappend tmp $id
    }
  }
  set theArray(tableEntryPanelIds) $tmp

}


# Proc creates a procedure entry panel for an entry field
#
proc Panel::editProcedure {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  set theArray(isActive) 0
  set theArray(hasFocus) 0

  #-Create a procedure entry widget and get the result
  set pe_id [ProcedureEntryPanel::create $globArray $fld Panel::procedureEntryPanelCallback]

  # Store procedure panel id
  lappend theArray(procedureEntryPanelIds) $pe_id

}


# Proc is the panel entry field update callback for a procedure entry panel
#
proc Panel::procedureEntryPanelCallback { pe_id globArray fld} {
  global Info
  upvar #0 $globArray theArray

  set entry_result [set ProcedureEntryPanel::result$pe_id]

  if { 1 == [set ProcedureEntryPanel::updated$pe_id] } {
    Panel::panelDataModified 1 $globArray ""
  }

  set theArray($fld) $entry_result
  set theArray($fld,procData) $entry_result

  set theArray($fld,variables) [set ProcedureEntryPanel::variableName$pe_id]

  if { [string equal -nocase "none" $theArray($fld,variables)] } {
    set theArray($fld,variables) ""
  }

  set theArray($fld,procData,variables) [set ProcedureEntryPanel::variableName$pe_id]

  set theArray($fld,dataSize) [list 1 1]
  set theArray($fld,procData,dataSize) [list 1 1]

  set theArray($fld,nofEntries) 1
  set theArray($fld,procData,nofEntries) 1

  #-Set new value to entry
  set wdg $theArray(allWidgets,$fld)
  if {[winfo exists $wdg]} {
    $wdg delete 0 end
    $wdg insert 0 $entry_result
  }

  #-Clear the procedure entry panel
  ProcedureEntryPanel::delete $pe_id

  #--Remove entry panel id from the list of active panels
  set tmp ""
  foreach id $theArray(procedureEntryPanelIds) {
    if { $id != $pe_id } {
      lappend tmp $id
    }
  }
  set theArray(procedureEntryPanelIds) $tmp

}

# ============================
# PANEL MODIFIED/CHANGED PROCS
# ============================

# Dirty when: panel data entry was touched but not "used"
# Updaters should reset this flag!!!
# NOTE: This is used for entries like timestep-entry, where you
# have to add data with "sub-add" buttons, which is different from
# parameter add-button!
#
proc Panel::panelDataDirty { is_dirty globArray {widget ""} {event_info ""} } {

  upvar #0 $globArray theArray

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }

  if { $widget != "" &&
       "disabled" == [$widget cget -state]
     } {
    return
  }

  set theArray(dataDirty) $is_dirty

}



# Panel data was modified but not "applied"
# Updaters should reset this flag!!!
#
proc Panel::panelDataModified { is_modified globArray {widget ""} {event_info ""} } {
  global Common Info Model
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" &&
       $globArray != $Common(panelArray,US)
     } {
    return
  }

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }

  if { $widget != "" &&
       "disabled" == [$widget cget -state]
     } {
    return
  }

  set theArray(dataModified) $is_modified

  # Set state for the panel Update-button, if such
  # exists and the array has something to update!

  # NOTE: If argument is the button itself, do not
  # change it state (to disabled), because it would block
  # -command option scripts for the button!!!

  set has_update_button 0
  set has_add_button 0

  if { [info exists theArray(panelUpdateButton)] &&
       [winfo exists $theArray(panelUpdateButton)] &&
       $widget != $theArray(panelUpdateButton)
     } {
    set has_update_button 1
  }

  if { [info exists theArray(panelAddButton)] &&
       [winfo exists $theArray(panelAddButton)] &&
       $widget != $theArray(panelAddButton)
     } {
    set has_add_button 1
  }

  # Is update or add currently allowed in the panel
  set accept_update 1
  set accept_add 1

  if { [info exists theArray(updateAllowed)] &&
       !$theArray(updateAllowed)
     } {
    set accept_update 0
  }

  if { [info exists theArray(addAllowed)] &&
       !$theArray(addAllowed)
     } {
    set accept_add 0
  }

  set has_data 1
  if { [info exist theArray(ids)] && $theArray(ids) == ""} {
    set has_data 1
  }

  # Set Update-button state
  #
  if { $has_update_button && $has_data } {
    
    # If Update-button should be highlighted
    #
    if { $is_modified && $accept_update } {
      Widget::configure1 $theArray(panelUpdateButton) -state normal
      #Widget::configure1 $theArray(panelUpdateButton) -relief raised
      Widget::configure1 $theArray(panelUpdateButton) -fg $Info(en_mod_color)

    } else {
      Widget::configure1 $theArray(panelUpdateButton) -state disabled
      #Widget::configure1 $theArray(panelUpdateButton) -relief raised
      Widget::configure1 $theArray(panelUpdateButton) -fg $Info(en_nmod_color)
    }
  }

  # Set Add-button state
  #
  if { $has_add_button && $has_data } {

    # If Add-button should be highlighted
    #
    if { $is_modified && !$accept_update && $accept_add } {
      Widget::configure1 $theArray(panelAddButton) -fg $Info(en_mod_color)

    } else {
      Widget::configure1 $theArray(panelAddButton) -fg $Info(en_nmod_color)
    }
  }

  if { $is_modified && !$has_update_button } {
    set is_changed 1
  } else {
    set is_changed $theArray(dataChanged)
  }

  Panel::panelDataChanged $is_changed $globArray $widget $event_info
}


# Modified data was "applied"
# Savers should reset this flag!!!
#
proc Panel::panelDataChanged {is_changed globArray {widget ""} {event_info ""} } {
  global Common Info Model
  upvar #0 $globArray theArray

  if { $Model(GEOMETRY_DIMENSION) == "" &&
       $globArray != $Common(panelArray,US)
     } {
    return
  }

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }

  if { $widget != "" &&
       "disabled" == [$widget cget -state]
     } {
    return
  }

  # Set changed/modified status
  set theArray(dataChanged) $is_changed

  if { !$theArray(dataChanged) } {
    #set theArray(dataModified) 0
  }

  # Mark window updated
  #
  if { [info exists theArray(winName)] &&
       [info exists theArray(winTitle)]  &&
       [winfo exists $theArray(winName)] 
     } {

    #--Panel(s) modfied and update(s) done ("deep modify")
    if { $theArray(dataChanged) && $theArray(dataModified) } {
      wm title $theArray(winName) "$theArray(winTitle)(*)"

    #--Panel(s) modfied and update(s) done ("deep modify")
    } elseif { $theArray(dataChanged) } {
      wm title $theArray(winName) "$theArray(winTitle)*"

    #--Visible panel (pane) was modified, but no update yet ("ligh modify")
    } elseif { $theArray(dataModified) } {
      wm title $theArray(winName) "$theArray(winTitle)(*)"

    #--No unsaved changes
    } else {
      wm title $theArray(winName) "$theArray(winTitle)"
    }
  }

  # Activate/Deactivate cancel-button
  # NOTE: currently defaultCancelState is normal!, MVe 04.10.00
  #
  if { [info exists theArray(cancelButton)] &&
       [winfo exists $theArray(cancelButton)]
     } {

    if {$is_changed} {
      Widget::configureS $theArray(cancelButton) normal

    } else {
      Widget::configureS $theArray(cancelButton) $Info(defaultCancelState)
    }
  }

  # Activate/Deactivate apply-button
  if { [info exists theArray(applyButton)] &&
       [winfo exists $theArray(applyButton)]
     } {

    if { $theArray(dataModified) || $theArray(dataChanged) } {
      Widget::configureS $theArray(applyButton) normal
      #Widget::configure1 $theArray(applyButton) -fg $Info(en_mod_color)

    } else {
      Widget::configureS $theArray(applyButton) $Info(defaultApplyState)
      #Widget::configure1 $theArray(applyButton) -fg $Info(en_nmod_color)
    }
  }
}


# Check is any field is modified
#
proc Panel::panelHasModifiedStatus {globArray} {
  upvar #0 $globArray theArray

  foreach fld $theArray(allFields) {

    # Only case variables are interesting
    if { !$theArray($fld,trg) } {
      continue
    }
    
    if { $theArray($fld,mod) } {
      return 1
    }
  }

  return 0
}


# Check is any field is in error state
#
proc Panel::panelHasErrorStatus {globArray} {
  upvar #0 $globArray theArray

  foreach fld $theArray(allFields) {

    # Only case variables are interesting
    if { !$theArray($fld,trg) } {
      continue
    }
    
    if { $theArray($fld,err) } {
      return 1
    }
  }

  return 0
}


# Check if panel data values were changed
#
proc Panel::panelValuesChanged {globArray} {
  upvar #0 $globArray theArray

  foreach fld $theArray(allFields) {

    # Only case variables are interesting
    if { !$theArray($fld,trg) } {
      continue
    }

    if { [Panel::isNumericField $globArray $fld] } {

      if { [expr $theArray($fld)] != [expr $theArray($fld,old)] } {
        return 1
      } 
    } else {
      if { $theArray($fld)] != $theArray($fld,old) } {
        return 1
      }
    }
  }

  return 0
}


# ===========
# OTHER STUFF
# ===========

# OK in a window.
# This is a general version which is currently
# just the same as cancel!
#
proc Panel::ok { arguments } {
  global Info

  set argc [llength $arguments]
  if { $argc == 0} {
    destroy $Info(thisWindow)
  } else {
    destroy $arguments
  }
}


# Check that new parameter name does not already exists or
# what to do when it is empty
#
proc Panel::checkParameterName {globArray pid name add_mode} {
  global Info
  upvar #0  $globArray theArray 

  # Default ("autmatic") name
  set default_name [Panel::defaultParameterName $globArray]

  # Empty name
  # ==========
  if { $name == "" } {

    set msg [list "NOTE: Parameter name is missing!\n" \
                  "Press Ok to use the default name: \n\n$default_name" ]

    if { ![Message::verifyAction $Info(noviceUser) $msg ] } {
      return 0

    } else {
      set theArray(parameterName) $default_name
      return 1
    }
  }

  set is_like_default [Panel::isLikeDefaultParameterName $globArray $name]
  set name_exists 0

  # Check against existing names
  # ============================
  foreach id $theArray(ids) {
    
    if { $id != $pid && $name == $theArray($id,name) } {
      set name_exists 1
      break
    }
    
  }

  # What to do, if name already exist
  # =================================
  if {$name_exists} {

    # Current name is like adefault name and we are adding
    # parameter, just generate the autoamtic name!!!
    #
    if { $add_mode && $is_like_default } {
      set theArray(parameterName) $default_name
      return 1


    #-ERROR: Name already exists, is is NOT accepted
    #
    } else {
      
      if { [info exists theArray(parameterEntryName)] } {
        set pname $theArray(parameterEntryName)
      } else {
        set pname "A parameter"
      }

      set msg [list "ERROR: $pname already exists with the name: \n\n$name"]
      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return 0
    }
  }

  # All ok
  # ======
  return 1
}


# Form default parameter name
#
proc Panel::defaultParameterName {globArray {pid ""} } {
  upvar #0  $globArray theArray 

  if { $pid == "" } {
    set pid $theArray(nextNewParameterId)
  }

  return $theArray(defaultName)$pid
}


# Check if parameter name is like a autoamtically
# generated name (like Contraint123)
#
proc Panel::isLikeDefaultParameterName {globArray name} {
  upvar #0  $globArray theArray 

  set name [string trim $name]

  # Form test-pattern like Constraint1234...
  #
  set pat $theArray(defaultName)
  append pat "\[0-9\]*"

  # Check if there is nothing else but numbers
  # after the default parameter name
  #
  regsub $pat $name "" rest

  if { $rest == "" } {
    return 1
  } else {
    return 0
  }

}



# Check if it is ok to select a new parameter, when
# data in the panel was modified, but not updated
#
proc Panel::verifyParameter {globArray} {
  global Info
  upvar #0 $globArray theArray 

  set msg [list "Panel data has been modified but is not updated!\n\n"  \
                $Info(anywayOk) ]

  set Info(messageIcon) warning

  if { "cancel" == [Message::dispOkCancelMessage $msg] } {
    return 0
  }

  return 1
}


#-Verifies that Save (OK or Apply) is all right if the
# panel data has been "modified" but not applied
# Returns:
# 0 <--> do not continue
# 1 <--> it is ok to continue
#
proc Panel::verifySave {globArray} {
  global Info
  upvar #0 $globArray theArray 

  #--Verify first that there are no non-applied
  #  data entry sub-panels for the panel
  #
  if { ![Panel::verifyEntryPanels $globArray] } {
    return 0
  }

  #--Panel name to be used in messages
  if { [info exists theArray(winTitle)] } {
    set pname $theArray(winTitle)
  } else {
    set pname $globArray
  }

  # NOTE: UPDATE has higher priority than ADD!
  # ==========================================

  set updateDone 0
  set addDone    0

  #--If panel data was modified, but not updated --> add or update!!!
  if { [info exists theArray(dataModified)] && $theArray(dataModified) } {

    # Check update
    if { !$updateDone                               &&
         [info exists theArray(panelUpdateButton)]  &&
         [winfo exists $theArray(panelUpdateButton)]
       } {
      set wdg $theArray(panelUpdateButton)
      set cmd [$wdg cget -command]
      set sta [$wdg cget -state]

      set rc ""

      # Activate update command if button active!
      if { $sta == "normal" && $cmd != "" } {
        set rc [eval $cmd]
        set updateDone 1
      }


      # Id update button returned something and
      # it was equal to zero, update not succesful!
      #
      if { $rc != "" && !$rc } {
        return 0
      }
    }

    # Check next auto add, will be don only if
    # update was not done!
    if { !$updateDone                             &&
         [info exists theArray(autoAdd)]          &&
         $theArray(autoAdd)                       &&
         [info exists theArray(panelAddButton)]   &&
         [winfo exists $theArray(panelAddButton)]
       } {

      set wdg $theArray(panelAddButton)
      set cmd [$wdg cget -command]
      set sta [$wdg cget -state]

      set rc ""

      # Activate add command if button active!
      if { $sta == "normal" && $cmd != "" } {
        set rc [eval $cmd]
        set addDone 1
      }

      # Id add button returned something and
      # it was equal to zero, update not succesful!
      #
      if { $rc != "" && !$rc } {
        return 0
      }
    }

  }

  #--If non applied data (like non-attached added parameters) inform
  #  (novice) user
  #
  if { $Info(userLevel) == $Info(noviceUser) } {
    
    if { $theArray(dataModified) ||
         ([info exists theArray(dataDirty)] && $theArray(dataDirty))
       } {

      set msg [list "NOTE: $pname\n\n" \
                    "Panel data was changed, but the new data was not used!\n\n" \
                    $Info(anywayOk) ]

      set Info(messageIcon) question

      if { "cancel" == [Message::dispCancelOkMessage $msg] } {
        return 0
      } 
    }
  }
 
  return 1
}


#-Verifies that cancel is all right if the
# panel data has changed ("modified" or "changed")
# Returns:
# 0 <--> do not cancel
# 1 <--> it is ok to cancel
#
proc Panel::verifyCancel {globArray} {
  upvar #0 $globArray theArray 
  global Info

  # NEW VERSION
  # Do not verify, just accepts cancel!!!
  # =====================================
  Panel::panelDataChanged 0 $globArray
  return 1


  # OLD VERSION
  # Verify chnages!!!
  # =================
  #--Verify first possibly only modified data
  if { ![verifySave $globArray] } {
    return 0
  }

  if { [info exists theArray(winTitle)] } {
    set pname $theArray(winTitle)
  } else {
    set pname $globArray
  }

  #--Then verify applied, but not saved data
  if { [info exists theArray(dataChanged)] } {
    if { $theArray(dataChanged) } {
      set msg [list "NOTE: $pname\nPanel data was changed, but not saved!\n\n" \
                     $Info(anywayOk) ]

      set Info(messageIcon) warning
      if { "cancel" == [ Message::dispOkCancelMessage $msg ] } {
        return 0
      }
    }
  }

  return 1
}


#-Verifies that cancel/ok/apply is all right if there
# are entry panels opnen which have "modified" but not applied
# Returns:
# 0 <--> do not continue
# 1 <--> it is ok to continue
#
proc Panel::verifyEntryPanels {globArray} {
  upvar #0 $globArray theArray 
  global Info

  # Check table entry panels
  if { [info exists theArray(tableEntryPanelIds)] } {

    foreach id $theArray(tableEntryPanelIds) {

      if { [info exist TableEntryPanel::updated$id] &&
           [set TableEntryPanel::updated$id]
         } {

        set fld [set TableEntryPanel::fieldName$id]
        set msg [list "NOTE: Data entry for $fld was changed, but not applied!\n\n" \
                       $Info(anywayOk) ]

        set Info(messageIcon) warning

        if { "cancel" == [ Message::dispOkCancelMessage $msg ] } {
          return 0
        }
      }
    }
  }

  # Check procedure entry panels
  if { [info exists theArray(procedureEntryPanelIds)] } {

    foreach id $theArray(procedureEntryPanelIds) {

      if { [set ProcedureEntryPanel::updated$id] } {

        set fld [set ProcedureEntryPanel::fieldName$id]
        set msg [list "NOTE: Procedure entry for $fld was changed, but not applied!\n\n" \
                       $Info(anywayOk) ]

        set Info(messageIcon) warning

        if { "cancel" == [ Message::dispOkCancelMessage $msg ] } {
          return 0
        }
      }
    }
  }

  return 1
}



#-Closes all open entry panels for the main panel
#
proc Panel::closeEntryPanels {globArray} {
  upvar #0 $globArray theArray 
  global Info

  if { [info exists theArray(tableEntryPanelIds)] } {
    foreach id $theArray(tableEntryPanelIds) {
      TableEntryPanel::cancel $id
    }
  }

  if { [info exists theArray(procedureEntryPanelIds)] } {
    foreach id $theArray(procedureEntryPanelIds) {
      ProcedureEntryPanel::cancel $id
    }
  }
}


# CANCEL in a window
# Destroys the window
#
proc Panel::cancel { {arguments ""} } {
  global Info gMW

  if { $arguments == ""} {
    destroy $Info(thisWindow)
  } else {
    global $arguments
    destroy $arguments
  }

  set tmp ""

  foreach row $Info(windowList) {

    set w [winfo atomname [lindex $row 1]]

    if {[winfo exists $w]} {
      lappend tmp $row
    }
  }

  set Info(windowList) $tmp

  MenuBuild::createWindowListMenu $gMW(Window,windowList)
}


# ============#
# FIELD TYPES #
# ============#

# Checks if the array field is currently active
#
proc Panel::isActiveField { globArray fld } {
  upvar #0 $globArray theArray

  if { ![info exists theArray($fld,act)] ||
       !$theArray($fld,act)
     } {
    return 0

  } else {
    return 1
  }
}


# Check is the array field is a data entry which must be checked
#
proc Panel::isCheckableDataField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray
  
  set wclass [DataField::getFieldProperty $globArray $fld WidgetClass ]
  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  set dsp [DataField::getFieldProperty $globArray $fld Display]

  #-When IS NOT checkable entry
  if { !$dsp ||
       $wclass == "Parent" ||
       ( $wtype != "BrowsableFile" &&
         $wtype != "Entry" &&
         $wtype != "IncludeFile" &&
         $wtype != "Text")
     } {
    return 0
  }

  return 1
}


# Check if the array field is a data (fileable data)
#
proc Panel::isSimpleDataField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]

  #-When NOT data entry
  if { $wtype  == "GroupScreen"  ||
       $wtype  == "Label"        ||
       $wtype  == "Separator"    ||
       $wtype  == "Button"       || 
       $wtype  == "ListBox"      ||
       $wtype  == "SelectionBox"
     } {
    return 0
  }
  
  return 1
}


# Check is the array field is a grouping entry
#
proc Panel::isGroupField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  #-When IS group entry
  if { $wtype == "GroupData" ||
       $wtype == "GroupDataAndScreen" ||
       $wtype == "GroupScreen"
     } {
    return 1
  }
  
  return 0
}


# Check is the array field is a grouping entry
#
proc Panel::isGroupMemberField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set wclass [DataField::getFieldProperty $globArray $fld WidgetType]

  #-When IS group member
  if { $wclass == "MemberOutputField" ||
       $wclass == "MemeberField"
     } {
    return 1
  }

  return 0
}



# Check if the array field should be read from model file
#
proc Panel::isInputField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set fclass [DataField::getFieldProperty $globArray $fld FieldClass ]

  #-When NOT data entry
  if { $fclass == "OutputOnly" } {
    return 0
  }

  return 1
}


# Check if the array field should be output to cpp side 
#
proc Panel::isOutputField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set wclass [DataField::getFieldProperty $globArray $fld WidgetClass]

  if { ![Panel::isSimpleDataField $globArray $fld] } {
    return 0
  }

  if { $wclass == "PanelField"  ||
       $wclass == "WorkField" 
     } {
    return 0
  }

  # Output status can ve controlled also by
  # the save flag
  if { [info exist theArray($fld,save)] &&
       !$theArray($fld,save)
     } {
    return 0
  }
  
  return 1  
}



# Check if the array field is a screen entry (affecting panels)
#
proc Panel::isScreenField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set wclass [DataField::getFieldProperty $globArray $fld WidgetClass]
  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  set dsp [DataField::getFieldProperty $globArray $fld Display]

  #-When NOT screen entry
  if { !$dsp ||
       $wclass == "ParentField" ||
       $wclass == "WorkField" ||
       $wtype  == "GroupData"
     } {
    return 0
  }

  return 1
}


# Check if the array field is a sif-field (SolverInputFile field)
#
proc Panel::isSifField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  # Must be at least an output field
  if { ![Panel::isOutputField $globArray $fld] } {
    return 0
  }

  if { 1 == [DataField::getFieldProperty $globArray $fld ExcludeFromSif] } {
    return 0
  }

  if { $globArray != $Common(panelArray,EQ) } {
    return 1
  }

  # Equation panel need special handling, because mask-system is not
  # controlling field's activity, so we use this to control equation
  # fields output to sif!!!
  
  # NOTE: We check field's target (panel activity) mask against problem mask!!!
  #
  set tmask $Equation(problemMask)
  if { ![DataField::fieldTargetMaskMatches Equation $fld $tmask] } {
    return 0
  }

  return 1
}


# Check if the array field is a radio-button type
#
proc Panel::isRadioButtonField {globArray fld} {
  global Common $globArray
  upvar #0 $globArray theArray

  set wclass [DataField::getFieldProperty $globArray $fld WidgetClass ]
  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]

  #-When IS radio-button type entry
  if {  $wclass  == "RadioButtonsContainer" ||
        $wtype  == "RadioButton"
     } {
    return 1
  }

  return 0
}


# ===================#
# FIELD ACCESS PROCS #
# ===================#

# Pick group fields for a group widget
#
proc Panel::getGroupFields {globArray fld} {
  global Common $globArray 
  upvar #0 $globArray theArray 

  set panel $theArray(parameterType)
  return $Common($panel,$fld,group)
}


# Pick group labels for a group widget
# Returns "" if label not defined
#
proc Panel::getGroupLabels {globArray fld} {
  global Common $globArray Coordinate
  upvar #0 $globArray theArray 

  set coord_index [DataField::getFieldProperty $globArray $fld CoordinateIndex]
  set result ""

  if {$coord_index > 0} {
    catch { set result [lindex $Common(CoordinateSymbols) $Coordinate(coordinateLabelIndex)] }
  } else {
    catch { set result [DataField::getFieldProperty $globArray $fld Label] }
  }

  return $result
}


# Pick proper label for a field.
# Returns "" if label not defined
#
proc Panel::getCommonLabel {globArray fld is_group} {
  global Common $globArray Coordinate
  upvar #0 $globArray theArray 

  set label [DataField::getFieldProperty $globArray $fld Label]
  set coord_index [DataField::getFieldProperty $globArray $fld CoordinateIndex]

  # If variable is coordinate dependent (like velocity)
  if {$coord_index > 0 && $is_group == 0} {
    catch { set symbols [lindex $Common(CoordinateSymbols) $Coordinate(coordinateLabelIndex)]  
            set symbol [lindex $symbols $coord_index]
            set label [append label $symbol]
          }
  }

  return $label
}



# ============
# OLD VERSIONS
# ============

# Panel data was modified but not "applied"
# Updaters should reset this flag!!!
#
proc Panel::panelDataModified_old { is_modified globArray {widget ""} {event_info ""} } {
  upvar #0 $globArray theArray

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }

  if { $widget != "" &&
       "disabled" == [$widget cget -state]
     } {
    return
  }

  set theArray(dataModified) $is_modified

  # Set state for the panel Update-button, if such
  # exists and the array has something to update!

  # NOTE: If argument is the button itself, do not
  # change it state (to disabled), because it would block
  # -command option scripts for the button!!!
  set has_button 0
  if { [info exists theArray(panelUpdateButton)] &&
       [winfo exists $theArray(panelUpdateButton)] &&
       $widget != $theArray(panelUpdateButton)
     } {
    set has_button 1
  }

  # Is update currently allowed in the panel
  set accept_update 1
  if { [info exists theArray(updateAllowed)] &&
       !$theArray(updateAllowed)
     } {
    set accept_update 0
  }

  set has_data 1
  if { [info exist theArray(ids)] && $theArray(ids) == ""} {
    set has_data 0
  }

  if { $has_button && $has_data } {

    if { $is_modified && $accept_update } {
      Widget::configure1 $theArray(panelUpdateButton) -state normal
      Widget::configure1 $theArray(panelUpdateButton) -relief raised

    } else {
      Widget::configure1 $theArray(panelUpdateButton) -state disabled
      Widget::configure1 $theArray(panelUpdateButton) -relief raised
    }
  }

  if { $is_modified && $accept_update } {
    set is_changed 1
  } else {
    set is_changed $theArray(dataChanged)
  }

  Panel::panelDataChanged $is_changed $globArray $widget $event_info
}


# ===========
# Modified data was "applied"
# Savers should reset this flag!!!
#
proc Panel::panelDataChanged_old {is_changed globArray {widget ""} {event_info ""} } {
  global Info
  upvar #0 $globArray theArray

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }

  if { $widget != "" &&
       "disabled" == [$widget cget -state]
     } {
    return
  }

  # Set changed/modified status
  set theArray(dataChanged) $is_changed

  if { !$theArray(dataChanged) } {
    set theArray(dataModified) 0
  }

  # Mark window updated
  #
  if { [info exists theArray(winName)] &&
       [info exists theArray(winTitle)]  &&
       [winfo exists $theArray(winName)] 
     } {
    if {$is_changed} {
      wm title $theArray(winName) "$theArray(winTitle)*"
    } else {
      wm title $theArray(winName) "$theArray(winTitle)"
    }
  }

  # Activate/Deactivate cancel-button
  if { [info exists theArray(cancelButton)] &&
       [winfo exists $theArray(cancelButton)]
     } {
    if {$is_changed} {
      Widget::configureS $theArray(cancelButton) normal
    } else {
      Widget::configureS $theArray(cancelButton) $Info(defaultCancelState)
    }
  }

  # Activate/Deactivate apply-button
  if { [info exists theArray(applyButton)] &&
       [winfo exists $theArray(applyButton)]
     } {
    if {$is_changed} {
      Widget::configureS $theArray(applyButton) normal
    } else {
      Widget::configureS $theArray(applyButton) $Info(defaultApplyState)
    }
  }
}


# end ecif_tk_procsPanel.tcl
# ********************
