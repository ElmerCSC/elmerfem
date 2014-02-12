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
#Module:    ecif_tk_procsDataField.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Data field handling procedures
#
#************************************************************************

# Inform abouyt fields emoved from parameters
#
proc DataField::informAboutRemovedFields { {msg_lst ""} } {
  global Model

  if { $Model(removedFieldsList) == "" } {
    return
  }
  
  set msg ""

  foreach row $msg_lst {
    lappend msg $row
  }

  foreach row $Model(removedFieldsList) {
    lappend msg $row
  }

  Message::dispOkMessage $msg
}



# Define all active fields for the argument arries
# by expanding indexed fields!
#
# Indices are defined in the "indexing" array,like
# EquationVariable(...)
#
# Start: theArray(initialFields)
# Result: theArray(allFields)
#
proc DataField::setAllFields {arries} {
  global BoundaryCondition Common Equation EquationVariable Info

  set update_variables 0
  set bc_extras ""  ;# These field names will be removed from BC-array
  set ic_extras ""  ;# These field names will be removed from IC-array

  # Loop all arries (which could have indexed variables)
  # ===============
  foreach arr $arries {
    upvar #0 $arr theArray
    
    #-Variables-array does not need any expanding!
    #
    if { $arr == "Variables" } {
      continue
    }

    #-If no initial fields available!
    #
    if { ![info exists theArray(initialFields)] } {
      return
    }
    
    set total_mask [Panel::getProblemMask $arr]
    set has_expanded 0
    set new_fields $theArray(initialFields)

    # Loop all fields in the array
    # ============================
    foreach fld $theArray(initialFields) {

      #-Check if some indexing is in use
      set index_type [ DataField::getFieldProperty $arr $fld IndexType ]

      if { $index_type == "" } {
        continue
      }

      set index_var [ DataField::getFieldProperty $arr $fld IndexVariable ]

      #-Data values separator in the index var
      set index_sep [ DataField::getFieldProperty EquationVariable $index_var FieldDataSep ]

      #-Get indices from the indexing var in the array

      set pid $EquationVariable(parameterId)
      set the_indices [DataField::getFieldValue EquationVariable $pid $index_var]

      if { $the_indices == "" } {
        continue
      }

      #---Expand field
      set has_expanded 1

      set indices [split $the_indices $index_sep]
      
      set index_also_mask [DataField::getFieldProperty $arr $fld IndexAlsoMask]
      set label [DataField::getFieldProperty $arr $fld Label]
      
      # We use Sif Name as the base varaible name,if it given, otherwise
      # we use field's label!
      #
      set found 0
      set sif_label [DataField::getFieldProperty $arr $fld SifName "" 0 found]

      if { !$found } {
        set sif_label $label
      }

      if { 1 == [llength $indices] } {
        set number_mask 0
      } else {
        set number_mask 1
      }

      #-Loop all indexing names
      #
      foreach index $indices {
        
        # Indexing variable's gui name
        #
        set index_ [DataField::variableNameSifToGui $index]
                  
        if { $index_type == "Pre" } {
          set new_var "<$index_>"
          append new_var $fld
          set new_label "$index $label"
          set sif_name "$index $sif_label"

        } else {
          set new_var $fld
          append new_var "<$index_>"
          set new_label "$label Of $index"
          set sif_name "$sif_label Of $index"
        }
        
        set new_widget_class [DataField::getFieldProperty $arr $fld ChildWidgetClass]

        #--Create own mask for the indexed variable
        
        # First pick the "template" field mask
        set new_mask [DataField::getFieldProperty $arr $fld Mask]

        set mask1 [lindex $new_mask 0]
        set mask2 [lindex $new_mask 1]

        # If an indexed variable should also add indexed mask components
        # like (¤AD#Oxygen .. or ¤UD010#SoundP ... for user defined fields)
        #
        if { $index_also_mask } {

          set msk2 ""
 
          if { $mask2 != "" } {
            set len [string length $mask2]
            set pos 0
            set add_var_msk 0  
            set check_next_ch 0

            while { $pos < $len } {
              
              set ch [string index $mask2 $pos]
              
              incr pos

              # If not a part teminator, add right away
              if { $ch != ")" } {
                append msk2 $ch
              }
                
              # NOTE: Variable name is added only for proper equation masks, which start
              # with a capital letter, not dimension or geometry related masks
              # like (¤dimS2) or (¤bcD) whic always should start with a small letter!
              #
              if { $check_next_ch } {
                if { $ch == [string toupper $ch] } {
                  set add_var_msk 1
                }
                set check_next_ch 0
              }

              # If a mask starter, mark that next character should be checked
              if { $ch == "¤" } {
                set check_next_ch 1
              }

              # Part terminator, check if variable mask should
              # be inserted
              if { $ch == ")" } {

                # Ok, insert variable mask !!!
                if { $number_mask && $add_var_msk } {

                  append msk2 [DataField::formVariableMask $index]
                }

                set add_var_msk 0
                append msk2 $ch
              }

            } ;# while mask2
          } ;# If mask2 != ""

          set mask2 $msk2

          set new_mask $mask1

          if { $mask2 != "" } {
            lappend new_mask $mask2
          }
        } ; # Index also mask

        DataField::setFieldProperty $arr $new_var TemplateVariable $fld

        # Set some field properties for the new variable
        #
        # NOTE: These properties must be set, otherwise they are inherited
        # (possibly incorrectly!) from the parent field owing to the
        # "TemplateVariable" property being set below!!!
        #
        DataField::setFieldProperty $arr $new_var Mask $new_mask
        DataField::setFieldProperty $arr $new_var WidgetClass $new_widget_class
        DataField::setFieldProperty $arr $new_var Label $new_label
        DataField::setFieldProperty $arr $new_var SifName $sif_name
        DataField::setFieldProperty $arr $new_var IndexType ""      ;# Important to nullify!
        DataField::setFieldProperty $arr $new_var IndexVariable ""  ;# Important to nullify!

        lappend new_fields $new_var

        # If new field's label is empty, it will be an equation-variable type field in
        # InitialCondition or BoundaryCondition arries and we must remove any instances
        # of these variables which may have come via SetVariableComponents-procedure
        # The same concerns the variable BoundaryCondition(initialDirichletNames)
        if { $label == "" } {
          if { $arr == $Common(panelArray,BC) } {
            lappend bc_extras [string toupper $index_ ]
          } elseif { $arr == "InitialCondition" } {
            lappend ic_extras [string toupper $index_ ]
          }
        }
        
      } ;#Each index
 
    } ;#Each fld 
 
    # Store expanded field names
    set theArray(allFields) $new_fields

  } ;#Each array

  
  # Remove possible duplicate fields
  # NOTE: Actually the user may have unintentionally
  # defiend fields in other equations wich are similar to
  # ex. AdvectionDiffusion variables, but we cannot differentiate between
  # that, sorry!!!###!!!
  #
  if { $bc_extras != "" } {
    global BoundaryCondition
    set BoundaryCondition(allFields) [
      List::removeFromList $BoundaryCondition(allFields) $bc_extras]
  }

  if { $ic_extras != "" } {
    global InitialCondition
    set InitialCondition(allFields) [
      List::removeFromList $InitialCondition(allFields) $bc_extras]
  }

}


# Form mask for a  variable name (= field index)
#
# NOTE: We have to use a special stop mark (like "^") to
# prevent partial matches like AD#species to AD#species2
# With AD#species^ and AD#species2^ there is no partial match!
#
proc DataField::formVariableMask {variable_name} {
  global Info

  #-Replace spaces with something else
  set vmn [string map {" " "_"} [string tolower $variable_name]]
  
  #-Form mask (like #species£ )
  set mask $Info(variableMaskStart)
  append mask $vmn
  append mask $Info(variableMaskStop)

  return $mask
}


# Create Equation Variable based fields in the following arries:
# Variables and possibly in InitialCondition, BoundaryCondition
#
proc DataField::setVariables {} {
  global Equation

  # Add names from predefined "indexed" equation (like AdvDiff and from
  # user defined "indexed" equations
  #
  foreach ename $Equation(allFields) {

    #-Check if this an equation name field
    #
    if { ![Equation::isEquation $ename] } {
      continue
    }

    if { ![DataField::getFieldProperty Equation $ename IsMultiVar] } {
      continue
    }

    set eindex [DataField::getFieldProperty Equation $ename EquationIndex]

    #-Get all variable names
    set evars [Equation::getAllVariableNames $eindex]

    #-Create new Variables-array fields from equation variables
    foreach var $evars {
      DataField::setVariableComponents $ename $var
    }
  }

  Panel::initFields Variables
}


# Add variable components based on variable name and degrees of freedom
# For example: SoundP and dof=2 will create variable components
# "SoundP 1" and "SoundP 2"
#
# Components are added to: Variables, InitialCondition and BoundaryCondition
# NOTE: Only Real equation vars are added to Boundary and InitialCondition panels
#       Calculator and Exported equation variables are not!
#
# NOTE: mask: equation level mask. It is the display mask!
# The variable mask (activity mask) is build from it and variable name
#
#
proc DataField::setVariableComponents { eq_name varname { dofs ""} { mask "" } } {
  global BoundaryCondition  InitialCondition UserDefined Variables

  set dofs_list [split $dofs ";"]
  set base_list ""
  set emask $mask

  #-Set variable properites
  # =======================
  #
  # NOTE: These properties are picked from equation field and where they are specially
  # named with the prefix "Var_", so that they are not mixed with
  # equation field's own properties!
  #
  set fieldProperties ""

  foreach prop_name $UserDefined(equationVarFieldKeywords) {
    
    if { ![UserDefined::convertPropertyNameToGuiFormat $prop_name gui_prop] } {
      continue
    }

    lappend fieldProperties $gui_prop
  }

  #-Set property values
  # ===================
  set propertyValues ""

  if { $eq_name != "" } {

    set ebase [DataField::getFieldProperty Equation $eq_name VariableComponents "" 0]
    set edofs [DataField::getFieldProperty Equation $eq_name VariableDofs "" 0]
    set emask [DataField::getFieldProperty Equation $eq_name EquationMask]

    set dofs_list [split $edofs ";"]
    set base_list [split $ebase ";"]

    # NOTE: These prperties are specially named in the equation field, so
    # we do not mix them with the properties of the equation field itself!
    #
    foreach prop $fieldProperties {
      
      set values [DataField::getFieldProperty Equation $eq_name "Var_$prop" "" 0]
      set values [split $values ";"]
      lappend propertyValues $values
    }
  }

  # NOTE: Currently index is always 0, but this could be used
  # to handle multi-name variables like (Pressure, Velocity) for
  # the 'Flow Solution'.
  #
  set index 0

  foreach dofs $dofs_list {

    # If more than one dof-value --> simulation dimension
    # dependency (like Velocity etc.)
    #
    if { 1 != [llength $dofs] } {
      set dim_dep 1

      if { 2 == [llength $dofs] } {
        set dofs [linsert $dofs 0 1]
      }

    } else {
      set dim_dep 0
    }
    
    set dof_max [lindex $dofs end]
    set idx 1
    set base [lindex $base_list $index]

    while { $idx <= $dof_max }  {

      #--We number components only if dof > 1
      #
      if { $dof_max > 1 } {
        set cn "$base $varname $idx"
      } else {
        set cn "$base $varname"
      }
      
      # For the case of an empty base!  
      set cn [string trim $cn]
      set cn_base [string trim "$base $varname"]

      #--Variable Gui (field) name
      set fn [DataField::fieldNameSifToGui $cn]

      #--Variable Sif name
      DataField::setFieldProperty "" $fn SifName $cn

      #--Set equation variable properties
      if { $eq_name != "" } {

        # Predefined (fixed) property values
        # ==================================
        DataField::setFieldProperty "" $fn WidgetType Entry
        DataField::setFieldProperty "" $fn FieldValueType Real

        # User entered properties
        # =======================
        #-Set actively set (non-empty) property values for the field
        #
        # NOTE: We use (of course) normal property names here!
        #
        foreach prop $fieldProperties values $propertyValues {
          
          # NOTE: Normally we have only one value in the list
          #
          set value [lindex $values $index]

          if { $value != "" }  {
            DataField::setFieldProperty "" $fn $prop $value
          }
        }

        # System set properties
        # =====================
        #
        # NOTE:
        # dsp_msk is equation level
        # act_msk is variable level
        #
        set dsp_mask $emask

        set act_mask [string trimright $emask ")"]
        append act_mask [DataField::formVariableMask $varname]
        append act_mask ")"
        
        # NOTE: We replace the possible user set (0/1) indicator
        # with the actual coordinate (dof) index
        #
        set is_coord [DataField::getFieldProperty "" $fn LabelByCoordinate]
        
        if { $is_coord } {
          DataField::setFieldProperty "" $fn CoordinateIndex $idx
        }
          
        # If dimension dependent field
        if { $dim_dep } {

          if { $idx <= [lindex $dofs 0] } {
            set act_mask $act_mask\&((¤dimS1)|(¤dimS2)|(¤dimS3))


          } elseif { $idx <= [lindex $dofs 1] } {
            set act_mask $act_mask\&((¤dimS2)|(¤dimS3))

          } else {
            set act_mask $act_mask\&(¤dimS3)
          }
        }
    
        # Variable field's properties
        # ===============================
        
        # Label
        # -----
        # If no user defined label, we create the default label from
        # the variable name

        set lbl [DataField::getFieldProperty "" $fn Label]
        
        #-Coordinate will be added to the base field label
        #
        if { $is_coord } {
        
          if { $lbl == "" } {
            set lbl $cn_base
          }

        #-No coordinated label, but we may have to use the Dof to
        # label the variables
        #
        } else {
    
          # If no label given, we use the variable Sifname (which
          # is already indexed) for the label
          #
          if { $lbl == "" } {
            set lbl $cn

          } elseif { $dof_max > 1 } {
            set lbl [string trimright $lbl " $idx"]
            set lbl "$lbl $idx"
          }
        }

        DataField::setFieldProperty "" $fn Label $lbl

        # Mask
        # ----
        # NOTE: Mask is combined with the possibly existing mask (like
        # for Displacement of Stress Analysis and Displacement of the
        # user defined non-elastic equation)
        # NOTE: Mask is set directly for BoundaryCondition and InitialCondition panels

        DataField::addFieldMaskItem BoundaryCondition $fn [list "|$dsp_mask" "|$act_mask"]
        DataField::addFieldMaskItem InitialCondition $fn [list "|$dsp_mask" "|$act_mask"]
        DataField::addFieldMaskItem Variables $fn [list "|$dsp_mask" "|$act_mask"]

        # Add Dirichlet-bc mask
        set m ""
        lappend m {}
        lappend m &(¤bcD)
        DataField::addFieldMaskItem BoundaryCondition $fn $m

      } ;# if equation
    

      #--If field not yet installed, add it
      if { -1 == [lsearch $Variables(initialFields) $fn] } {

        lappend Variables(initialFields) $fn

        if { $eq_name != "" } {
          lappend InitialCondition(initialFields) $fn
          lappend BoundaryCondition(initialFields) $fn
          lappend BoundaryCondition(initialDirichletVnames) $fn
        }
      }

      incr idx

    } ;# for 1...max_dof

  } ;# For all dofs/base pairs

  set Variables(allFields) [lsort $Variables(initialFields)]
  
  # If this was for a real equation
  #
  if { $eq_name != "" } {
    #set InitialCondition(allFields) $InitialCondition(initialFields)
    #set BoundaryCondition(allFields) $BoundaryCondition(initialFields)
    #set BoundaryCondition(dirichletVnames) $BoundaryCondition(initialDirichletVnames)
  }
}


# Remove  variable components based on variable name and degrees of freedom
# For example: SoundP and dof=2 will remove variable components
# "SoundP 1" and "SoundP 2"
#
# Components are removed from: Variables, InitialCondition and BoundaryCondition
# NOTE: Only "equation" vars are in Boundary and InitialCondition panels
#       Calculator variables are not!
#
proc DataField::unsetVariableComponents { varname dofs {is_equation_var 1} } {
  global Variables InitialCondition BoundaryCondition

  # If more than one dof-value --> simulation dimension
  # dependency (like Velocity etc.)
  #
  if { 1 != [llength $dofs] } {
    set dim_dep 1
  } else {
    set dim_dep 0
  }

  set dof_max [lindex $dofs end]
  set idx 1

  while { $idx <= $dof_max }  {

    #--We number components only if dof > 1
    if { $dof_max > 1 } {
      set cn "$varname $idx"
    } else {
      set cn "$varname"
    }

    #--Variable field name
    set fn [DataField::fieldNameSifToGui $cn]

    #-Variables array
    set index [lsearch $Variables(initialFields) $fn]
    if { $index != -1 } {
      set Variables(initialFields) [lreplace $Variables(initialFields) $index $index]
    }

    if { $is_equation_var } {

      #-InitialCondition array
      set index [lsearch $InitialCondition(initialFields) $fn]
      if { $index != -1 } {
        set InitialCondition(initialFields) [lreplace $InitialCondition(initialFields) $index $index]
      }

      #-BoundaryCondition array
      set index [lsearch $BoundaryCondition(initialFields) $fn]
      if { $index != -1 } {
        set BoundaryCondition(initialFields) [lreplace $BoundaryCondition(initialFields) $index $index]
      }

      set index [lsearch $BoundaryCondition(initialDirichletVnames) $fn]
      if { $index != -1 } {
        set BoundaryCondition(initialDirichletVnames) [lreplace $BoundaryCondition(initialDirichletVnames) $index $index]
      }
    }

    incr idx
  }

  set Variables(allFields) [lsort $Variables(initialFields)]
  
  if { $is_equation_var } {
    #set InitialCondition(allFields) $InitialCondition(initialFields)
    #set BoundaryCondition(allFields) $BoundaryCondition(initialFields)
    #set BoundaryCondition(dirichletVnames) $BoundaryCondition(initialDirichletVnames)
  }
}



###########################
### GENERAL FIELD PROCS ###
###########################


proc DataField::textWidget2Data {globArray fld text_wdg {trim_empty_lines 0} } {
  global Info
  upvar #0 $globArray theArray

 
  if { ![winfo exists $text_wdg] } {
    return
  }

  #---Split text entry into lines
  set data [$text_wdg get 0.0 end]
  #set values [split [string trimright $data \n] \n]
  set values [split $data \n]

  #---Make a list like "Name1;Name2;Name3"
  set sep $Info(dataListSeparator) 
  set result ""

  foreach val $values {
    set val [string trimleft $val]
    set val [string trimright $val]

    # Possibly skip empty lines
    if { $trim_empty_lines && $val == "" } {
      continue
    }

    append result "$val$sep"
  }
  
  set result [string trim $result $sep]

#print "data=$data values=$values result=$result"

  set theArray($fld) $result
}


proc DataField::data2TextWidget { globArray fld text_wdg {trim_empty_lines 1}} {
  global Info
  upvar #0 $globArray theArray

  if { ![info exists theArray($fld)] } {
    return
  }

  if { ![winfo exists $text_wdg] } {
    return
  }

  set sep $Info(dataListSeparator) 

  set values [split $theArray($fld) $sep]

  set result ""
  foreach val $values {
    append val \n
    append result $val
  }
  
  set result [string trimright $result \n]
  $text_wdg delete 0.0 end
  $text_wdg insert 0.0 $result
}


# Copy field propeties from panel2(fld2) to panel1(fld1)
#
proc DataField::copyFieldProperties {array1 fld1 array2 fld2} {
  global Common

  foreach pname $Common(allFieldProperties) {

    set pval [ DataField::getFieldProperty $array2 $fld2 $pname ]

    DataField::setFieldProperty $array1 $fld1 $pname $pval
  }
}


############################
### FIELD PROPERTY PROCS ###
############################


# Change field gui name to sif name
# Like "VELOCITY_1" --> Velocity 1"
#
proc DataField::fieldNameGuiToSif { gui_name } {

  set sif_name ""
  set next_to_upper 1

  set len [string length $gui_name]
  set i 0

  while { $i < $len } {

    set c [string index $gui_name $i]
    incr i

    if { $next_to_upper } {
      set c [string toupper $c]
    } else {
      set c [string tolower $c]
  }

    # Like "Heat Equation"
    if { $c == "_" } {
      set c " "
      set next_to_upper 1

    # Like "Navier-Stokes"
    } elseif { $c == "-" } {
      set next_to_upper 1

    # Like "My_Name"
    } elseif { $c == "~" } {
      set c "_"
      set next_to_upper 1

    } else {
      set next_to_upper 0
    }

    append sif_name $c
  }

  return $sif_name

}


# Change field sif name to gui name
# Like "Velocity 1" --> "VELOCITY_1"
#
proc DataField::fieldNameSifToGui { sif_name } {

  set gui_name ""

  set len [string length $sif_name]
  set i 0

  while { $i < $len } {

    set c [string index $sif_name $i]
    incr i

    set c [string toupper $c]

    if { $c == " " } { 
      set c "_"

    } elseif { $c == "_" } {
      set c  "~"
    }

    append gui_name $c
  }

  return $gui_name
}


# Change variable (index) gui name to sif name
# Like "oxygen" --> "OxyGen" (it it was given by the user
# in that format)
#
proc DataField::variableNameGuiToSif { gui_name } {
  global EquationVariable

  set sif_name [DataField::fieldNameGuiToSif $gui_name]

  # Get actual variable name as it is entered by the user!!
  #
  return [Equation::getMatchingVariableName $sif_name $EquationVariable(allVariableNames)]
}


# Change variable (index) sif name to gui name
# Like "Some Quantity" --> "some_quantity"
#
proc DataField::variableNameSifToGui { sif_name } {

  set gui_name [DataField::fieldNameSifToGui $sif_name]
  
  # Gui variable (index) names are in lower case!!!
  return [string tolower $gui_name]
}



#--------------------------#
# Get field property procs #
#--------------------------#

# Pick the property value for the field in the panel
# from Common array
#
proc DataField::getFieldProperty { globArray field property
                                   {system_index ""}
                                   {show_error_msg 1}
                                   {found_flag_var ""}
                 } {
  global Common

  if { $found_flag_var != "" } {
    upvar $found_flag_var found_flag
  }

  set found_flag 1

  # Use a generic definition
  if { $globArray == "" || $globArray == "#"} {
    set panel "#"

  # Use an array specific definition
  } else {
    upvar #0 $globArray theArray
    set panel $theArray(parameterType)
  }

  if { $system_index == "" } {

    if { [info exists theArray(panelsByEquation)] &&
         $theArray(panelsByEquation) &&
         [info exists theArray(equationIndex)] &&
         -1 == [lsearch $Common(propertiesNotByEquation) $property]
       } {
      set system_index $theArray(equationIndex)
    }

  } elseif { $system_index == - 1 } {
    set system_index ""
  }

  # If direct call succeds 
  # ======================
  if { [info exists Common($panel,$field,$property) ] } {
    set values $Common($panel,$field,$property)
    set value [DataField::getValueBySystemIndex $globArray $field $property $values $system_index]
    return [string trim $value]
  }

  # Try next if a template variable definition exists
  # =================================================
  if { [info exists Common($panel,$field,TemplateVariable) ] } {
    set template $Common($panel,$field,TemplateVariable)

    if { 2 == [llength $template] } {
      set arr [lindex $template 0]
    } else {
      set arr $globArray
    }

    set fld [lindex $template end]

    set fnd_flag 0
    set value [DataField::getFieldProperty $arr $fld $property $system_index 0 fnd_flag]
    set value [string trim $value]

    # Ok, something found via template variable
    # =========================================
    if { $fnd_flag } {
      return $value
    }
  }

  # Try next generic property for the field
  # =======================================
  if { [info exists Common(#,$field,$property) ] } {
    set values $Common(#,$field,$property)
    set value [DataField::getValueBySystemIndex $globArray $field $property $values $system_index]
    return [string trim $value]
  }

  # Try next generic property for the panel
  # =======================================
  if { [info exists Common($panel,default,$property) ] } {
    set values $Common($panel,default,$property)
    set value [DataField::getValueBySystemIndex $globArray $field $property $values $system_index]
    return [string trim $value]
  }

  # Use generic property
  # ====================
  if { [info exists Common(#,default,$property) ] } {
    set value [DataField::getDefaultProperty $globArray $field $property $system_index]
    return [string trim $value]
  }

  # ERROR, property not defined!
  # ============================
  if { $show_error_msg } {
    set msg "Unknown field property!\n\n"
    append msg "Panel: $panel\nField: $field\nProperty: $property\n\n"
    set Info(messageIcon) error
    Message::dispOkMessage $msg
    showStack 12
  }

  # No result in fact, but return blank value!
  # =========================================
  set found_flag 0
  return ""
}


proc DataField::getDefaultProperty { globArray field property system_index} {
  global Common

  if { $property != "WidgetWidth"  &&
       $property != "Mask"         &&
       $property != "FieldFormat"  &&
       $property != "Procedurable" &&
       $property != "Tableable"    &&
       $property != "Variabled"
     } {
    return $Common(#,default,$property)
  }

  set is_entry 0
  set is_real 0
  set is_integer 0
  set is_logical 0
  set is_string 0
  set is_file 0

  set fv_type [DataField::getFieldProperty $globArray $field FieldValueType $system_index]

  if {"Real" == $fv_type } {
    set is_real 1
  } elseif { "Integer" == $fv_type } {
    set is_integer 1
  } elseif { "Logical" == $fv_type } {
    set is_logical 1
  } elseif { "String" == $fv_type } {
    set is_string 1
  } elseif { "File" == $fv_type } {
    set is_file 1
  }


  set wd_type [DataField::getFieldProperty $globArray $field WidgetType $system_index]

  if { "Entry" == $wd_type } {
    set is_entry 1
  }


  # WidgetWidth
  # ===========
  if { $property == "WidgetWidth" } {
    if {$is_entry} { 
      return $Common(#,default,WidgetWidth)
    } else {
      return 1
    }

  # FieldFormat
  # ===========
  } elseif { $property == "FieldFormat" } {
    #-Real Entry
    if { $is_entry && $is_real } { 
      return $Common(#,default,FieldFormat)
    #-Integer Entry
    } elseif { $is_entry && $is_integer } { 
      return "%d"
    #-Other
    } else {
      return "%s"
    }

  # Prcedurable
  # ===========
  } elseif { $property == "Procedurable" } {
    if { !$is_entry || !$is_real } { 
      return 0

    } else {
      return $Common(#,default,Procedurable)
    }

  # Tableable
  # =========
  } elseif { $property == "Tableable" } {
    if { !$is_entry || !$is_real } { 
      return 0
    } else {
      return $Common(#,default,Tableable)
    }

  # Variabled
  # =========
  } elseif { $property == "Variabled" } {
    if { !$is_entry || !($is_real || $is_integer) } { 
      return 0
    } else {
      return $Common(#,default,Variabled)
    }
  }
}


proc DataField::getValueBySystemIndex {globArray field property values system_index} {
  global Info

  # Get "raw" value
  if { $system_index == -2 } {
    return $values
  }

  set index_found 0

  # If system index given
  # =====================
  if { $system_index != "" } {

    foreach ival $values {

      set is_index 0
      set sindex [string trim [lindex $ival 0]]

      if { "#" == [string index $sindex 0] } {
        set sindex [string range $sindex 1 end]
        set is_index 1
      }
    
      # Indexed value
      if { $is_index } {
    
        if { $sindex == $system_index } {
          set value [join [lrange $ival 1 end]]
#MSG "index=$system_index ival: values=$values  value=$value"
          set index_found 1
          break
        }
      }
    }
  }

  # If no system index given or found
  # =================================
  if { !$index_found } {
    
    # Check if first non-blank is "#" <--> no non-indexed value given
    set tmp [string trimleft $values "\{"]
    set tmp [string trim $tmp]

    # So, no non-indexed value given!
    if { "#" == [string index  $tmp 0] } {
      return [DataField::getDefaultProperty $globArray $field $property ""]
    
    # Ok, pick the non-indexed value (should be first in the list!)
    } else {
      set tmp [lindex [split $values "#"] 0]
      set value [string trimright [string trimright $tmp] "\{"]
    }
  }

  # Remove string separators fron the value
  # NOTE: We cannot store quoted strings in model file parameters!
  #set value [string trim $value "\""]

  return [string trim $value]
}


proc DataField::getFieldSize { globArray field } {
  global Model

  # Default values
  set size1 1
  set size2 1
  
  # Check if field is an array
  if { "Array" == [DataField::getFieldProperty $globArray $field FieldDataType] } {
    set is_array 1
  } else {
    set is_array 0
  }

  # Symbol DIM in the definition is replaced with this value
  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    set dim 3
  } else {
    set dim 2
  }

  set size [DataField::getFieldProperty $globArray $field "FieldDataSize"]

  # Check if size is a symbolic expression containing parenthesis
  set lpc 0
  set rpc 0
  set pos -1

  for {set i 0} {$i < [string length $size]} {incr i} {

    if { [string index $size $i] == "(" } {
      incr lpc
    }
  
    if { [string index $size $i] == ")" } {
      incr rpc
    }

    if { $lpc > 0 && $lpc == $rpc } {
      set pos  $i
      break
    }
  }

  # If parenthesis found, split at the "pos"
  if { $pos > 0 } {
    set sz1 [string range $size 0 $pos]
    incr pos
    set sz2 [string range $size $pos end]
  
  # Otherwise split at the first blank
  } else {
    set sz1 [lindex $size 0]
    set sz2 [lindex $size 1]
  }
  
  set sz1 [string toupper $sz1]
  set sz2 [string toupper $sz2]

  set dc1 [regsub -all "DIM" $sz1 "$dim" sz1]
  set dc2 [regsub -all "DIM" $sz2 "$dim" sz2]
  
  # NOTE: If size given by DIM and field is not marked as an array, we add
  # 1 into row/column sizes, because a scalar is then also  accepted in additon to a
  # possible 'tensor' type

  if { $sz1 != "" } {
    if { $dc1 > 0 && !$is_array } {
      set size1 [list 1 [expr $sz1]]
    } elseif { $sz1 == "N" } {
      set size1 $sz1
    } else {
      set size1 [list [expr $sz1]]
    }
  }

  if { $sz2 != "" } {
    if { $dc2 > 0 && !$is_array } {
      set size2 [list 1 [expr $sz2]]
    } elseif { $sz2 == "N" } {
      set size2 $sz2]
    } else {
      set size2 [list [expr $sz2]]
    }
  }

  return [list $size1 $size2]
}


# Reformats field data as list, which is
# suitable fex. for TableEntry
#
proc DataField::getFieldDataAsList {globArray fld} {
  global Info $globArray
  upvar #0 $globArray theArray

  set avn $globArray
  append avn "($fld)"

  if { ![info exists $avn] } {
    return ""
  }
  
  if { -1 != [string first $Info(dataListSeparator) $theArray($fld)] } {
    set values [split $theArray($fld) $Info(dataListSeparator)]
  } else {
    set values [split $theArray($fld) " "]
  }

  return $values
}


# Get initial value for the array var
# (possibly for the system defined by "system_index")
#
proc DataField::getInitialValue {globArray var {system_index ""} } {
  global $globArray
  upvar #0 $globArray theArray

  set ival [DataField::getFieldProperty $globArray $var InitialValue $system_index]

  #-If no value defined
  if { $ival == "" } {
    return ""
  }
  
  set ival [string trim [string trim $ival] "\""]

  #-Final value
  return [DataField::evalInitialValue $ival $system_index]
}


proc DataField::evalInitialValue {initial_value {system_index ""} } {

  set initial_value [Util::stringRemoveQuotes $initial_value]

  # If default value is a variable (possibly array) name
  if {"$" == [string index $initial_value 0] } {
    # Pick possible matrix name and make it global
    set def_name [split $initial_value "("]
    set def_name [string range [lindex $def_name 0] 1 end]

    if { $def_name != "" } {
      global $def_name
    }

    set ival [subst $initial_value]
    return $ival

  # If default value is a proc name
  #
  } elseif {"@" == [string index $initial_value 0] } {

    set def_proc [string range $initial_value 1 end]

    #-If some arguments, first must be "system_index"!!
    if { "" != [info args $def_proc] } {
      set ival [eval $def_proc $system_index]

    #-No argument
    } else {
      set ival [eval $def_proc]
    }

    return $ival

  # Normal data
  } else {
    return $initial_value
  }
}


proc DataField::getPanelPage {globArray fld {problem_index 0} } {
  upvar #0 $globArray theArray

  return [DataField::getFieldProperty $globArray $fld PanelPage $problem_index]
}


proc DataField::getSubPanel {globArray fld {problem_index 0} } {
  upvar #0 $globArray theArray

  return [DataField::getFieldProperty $globArray $fld SubPanel $problem_index]
}


proc DataField::getSifName {globArray fld {problem_index 0} } {
  upvar #0 $globArray theArray

  set sif_name [DataField::getFieldProperty $globArray $fld SifName $problem_index]

  if { $sif_name == "" } {
    set sif_name [DataField::fieldNameGuiToSif $fld]
  }

  return $sif_name
}


# Get mask for the field
# For group fields,mask is picked for the member field
# based on group field's value!
# mindex = 0 --> state mask, mindex = 1 --> displayed mask
#
proc DataField::getFieldMask {globArray fld mindex} {
  global Common
  upvar #0 $globArray theArray 
  
  if { $mindex < 0 || $mindex > 1 } {
    return ""
  }
      
  set f_mask ""

  set pt $theArray(parameterType)

  #-Pick gield's mask by mask index
  set f_mask [lindex [DataField::getFieldProperty $globArray $fld Mask] $mindex]
#MSG "fld=$fld m=$f_mask"
  # Single field
  # ============
  if { ![Panel::isGroupField $globArray $fld] } {
    return $f_mask
  }

  # Group field
  # ===========
  # NOTE we suppouse that group field values are binary (two valued) sets!
  set groupFlds $Common($pt,$fld,group)

  foreach gf $groupFlds {
   
    # Pick group member value which is valid when member is turned on
    set value [lindex [DataField::getFieldProperty $globArray $gf Limits] 1]
    
    if { [info exist theArray($fld)] && $value == $theArray($fld) } {
      set f_mask [lindex [DataField::getFieldProperty $globArray $gf Mask] $mindex]
      break
    }
  }

  return $f_mask
}



#--------------------------#
# Set field property procs #
#--------------------------#

proc DataField::setFieldProperty_old {globArray field property value} {
  global Common

  if { $globArray == "" || $globArray == "#"} {
    set panel "#"

  } else {
    upvar #0 $globArray theArray
    set panel $theArray(parameterType)
  }

  set Common($panel,$field,$property) $value
}


# Set field property value for an array var
# (possibly for the system defined by "system_index")
#
proc DataField::setFieldProperty {globArray fld property value {system_index ""} } {
  global Common
  upvar #0 $globArray theArray

  if { $globArray == "" || $globArray == "#"} {
    set panel "#"

  } else {
    upvar #0 $globArray theArray
    set panel $theArray(parameterType)
  }

  if { $system_index == "" } {
    set Common($panel,$fld,$property) $value
    return
  }

  set value [string trim $value]
  set new_value [list "#$system_index" $value]
  set new_values ""
  set old_found 0
  
  # Check if new value is already stored
  #
  set values [DataField::getFieldProperty $globArray $fld $property -2]

  foreach vals $values {
    
    set idx [lindex $vals 0]

    if { $idx != "#$system_index" } {
      lappend new_values $vals

    } else {
      set old_found 1  
      lappend new_values $new_value
    }
  }

  # Add if not stored
  #
  if { !$old_found } {
    lappend new_values $new_value
  }

  set Common($panel,$fld,$property) $new_values
}


proc DataField::setFieldProperties1 {globArray field properties values} {

  foreach property $properties value $values {
    DataField::setFieldProperty $globArray $field $property $value
  }
}


proc DataField::setFieldProperties2 {globArray field property_value_pairs} {

  foreach property_value_pair $property_value_pairs {
    set property  [lindex $property_value_pair 0]
    set value     [lindex $property_value_pair 1]
    DataField::setFieldProperty $globArray $field $property $value
  }
}


# Set default values for some field properties if needed
#
proc DataField::setDefaultFieldProperties {array fields} {

  foreach fld $fields {

    set wt [DataField::getFieldProperty $array $fld WidgetType]
    
    #-If widget type not yet set
    if { $wt == "" } {
      set fvt [DataField::getFieldProperty $array $fld FieldValueType]

      if { $fvt == "Logical" } {
        DataField::setFieldProperty $array $fld WidgetType "CheckBox"
      } else {
        DataField::setFieldProperty $array $fld WidgetType "Entry"
      }
    }
  }
}


# Set initial value for the array var
# (possibly for the system defined by "system_index")
#
proc DataField::setInitialValue {globArray var value {system_index ""} } {
  global $globArray

  DataField::setFieldProperty $globArray $var "InitialValue" $value $system_index
}



#------------------#
# Field mask procs #
#------------------#

proc DataField::addFieldMaskItem {globArray vnames item} {

  foreach fn $vnames {

    #--Pick field Mask
    set pm [DataField::getFieldProperty $globArray $fn "Mask"]

    set pm0  [lindex $pm 0]
    set pm1  [lindex $pm 1]

    set itm0 [lindex $item 0]
    set itm1 [lindex $item 1]

    if { $itm0 != "" } {
      if { $pm0 == "" } {
        set itm0 [string trim $itm0 "|"]
        set itm0 [string trim $itm0 "&"]
      }

      append pm0 $itm0
    }

    if { $itm1 != "" } {
      if { $pm1 == "" } {
        set itm1 [string trim $itm1 "|"]
        set itm1 [string trim $itm1 "&"]
      }

      append pm1 $itm1
    }

    #--Collect result from updated components
    set pm $pm0

    if { $pm1 != "" } {
      lappend pm $pm1
    }

    #--Update field Mask
    DataField::setFieldProperty $globArray $fn "Mask" $pm
  }
}


# Check if a field should be included on panel data fields
#
proc DataField::fieldProblemMaskMatches {globArray fld mask} {
  upvar #0 $globArray theArray

  # Empty mask matches any field!
  if { $mask == "" } {
    return 1
  }

  set fp_mask [DataField::getFieldMask $globArray $fld 0]

  return [Util::masksAreMatching $fp_mask $mask]
}


# Check if a field should be displayed (packed) on the panel
#
proc DataField::fieldDisplayMaskMatches {globArray fld mask} {
  upvar #0 $globArray theArray

  # Empty mask matches any field!
  if { $mask == "" } {
    return 1
  }

  set fd_mask [DataField::getFieldMask $globArray $fld 0]

  return [Util::masksAreMatching $fd_mask $mask]
}


# Check if a (displayed) field should be active (not disabled) on the panel
#
proc DataField::fieldTargetMaskMatches {globArray fld mask} {
  upvar #0 $globArray theArray

  # Empty mask matches any field!
  if { $mask == "" } {
    return 1
  }

  set ft_mask [DataField::getFieldMask $globArray $fld 1]
    
  return [Util::masksAreMatching $ft_mask $mask]
}

  
##################################
### Data validiating routines ###
##################################

#-All check-procedures return empty string if valued checked is OK 
# and a string like "value for a number" if ERROR
#
# formattedvalue: 'global' variable for the result
# actualsize: 'global' variable for the actual size
# value: entered value
# frmt: format to be used for the result
# valuetype: 
# "nbr" (numeric), "file" (file name), "dir" (directory name), "str" (string)
# "proc" (procedure name), "matc" (matc script), "lib" (library name), "func" (function name), "var" (variable name)
# needed_level: 0=empty accepted, 1=empty warned, 2=empty is error!
#-Returns possible error message --> "" is OK. 
#
proc DataField::checkValue { 
      formattedvalue
      actualsize
      value
      frmt
      {valuetype "nbr"}
      {needed_level 0}
      {range ""}
      {valueset ""}
      {var_size 0} 
      {data_size 1} 
      {accept_scalar 1}
      {min_max_sizes ""}
     } {

  global Info
  upvar $formattedvalue result_value
  upvar $actualsize result_size

  set result_value ""
  set result_size ""

  # Trim and remove surrounding quotes
  #
  set value [Util::stringTRQ $value]

  # Value is a Matc-command, apply it first
  # NOTE: precision set to 8 digits
  #
  if { $valuetype == "matc" } {
    set to_do "format($Info(defaultMatcPrecision),\"rowcol\");$value"
    set value [Util::doMatc $to_do]
    set valuetype "nbr"
  }

  #-case: empty field
  if { $value == "" } {
    
    # Blank accepted
    if {$needed_level == 0 } {
      return ""

    # Blank warned (no warning currently!)
    } elseif { $needed_level == 1 } {
      return ""

    # Blank is ERROR!
    } else {
      return "missing value!"
    }
  }


  #-Number or a procedure/file name
  #-NOTE: message used as an UPVAR

  # Number
  # ======
  if { $valuetype == "nbr" } {
    set check_result [DataField::checkNumbers $value $frmt message $range $valueset \
                                 $var_size $data_size $accept_scalar $min_max_sizes]
    set result_value [lindex $check_result 0]
    set result_size  [lrange $check_result 1 end] ;# data-size + nof-entries

  # File name
  # =========
  } elseif { $valuetype == "file" } {
    set result_value [DataField::checkFileName $value $frmt message]

  # Directory name
  # ==============
  } elseif { $valuetype == "dir" } {
    set result_value [DataField::checkDirectoryName $value $frmt message]

  # String
  # ======
  } elseif { $valuetype == "str" } {
    set result_value [DataField::checkString $value $frmt message]

  # Procedure name
  # ==============
  } elseif { $valuetype == "proc" } {
    set result_value [DataField::checkProcedureName $value $frmt message]

  # Library name
  # ============
  } elseif { $valuetype == "lib" } {
    set result_value [DataField::checkLibraryName $value $frmt message]

  # Function name
  # =============
  } elseif { $valuetype == "func" } {
    set result_value [DataField::checkFunctionName $value $frmt message]

  # Variable name
  # =============
  } elseif { $valuetype == "var" } {
    set result_value [DataField::checkVariableName $value $frmt message]

  } else {
    Message::showMessage "UNKNOWN value type ($valuetype) for value ($value) in DataField::checkValue"
    set message ""
    set result_value ""
  }

  return $message
}


# Checks vector type data field
# formattedvalue: result as formatted (by data type)
# maxElements: max number of elements in the vector
# other arguments as in proc checkValue.
#-Returns possible error message
#
proc DataField::checkVector {
          formattedvalue
          value
           frmt
          {datatype 0}
          {maxElements 3}
          {needed_level 0} } {

  upvar $formattedvalue result

  set result ""

  #set value [string trim $value " \t\n"]
  set value [string trim $value " "]

  #-empty field
  if { $value == "" } {

    # Blank accepted
    if {$needed_level == 0 } {
      return ""

    # Blank warned (no warning given currently!)
    } elseif { $needed_level == 1 } {
      return ""

    # Blank is ERROR!
    } else {
      return "missing value!"
    }
  }

  #-procedure name
  if {$datatype == 1} {
    set result [DataField::checkProcedureName $value $frmt message]
    return $message
  }

  #-vector
#print "vector: $value"
  set count 0 ; set result ""
  set list1 [split $value ","]
  foreach el1 $list1 {
    set list2 [split $el1 " "]
      foreach el2 $list2 {
      if {$el2 == ""} {
        continue
      }
      incr count
      if {$count > $maxElements} {
        return "Too many elements in vector!"
      }
      set res_val [DataField::checkNumber $el2 frmt message]
      if {$message == ""} {
        set result "$result$res_val, "
      } else {
        return $messsage
      }
    }
  }

  #-Format result and return ok.
  for {set i 0} { $i < [expr $maxElements - $count]} {incr i} {
    set result "$result$defaultvalue, "
  }

  #-Trim away extra spaces
  set result [string trimright $result " ,"]
  return $message
}


#-Checks numbers (integers, floats etc)
#-Data type can be: scalar or array
#-Possible error msg is put into 'message' variable
#-Returns a list: { {formatted result} {actual size} }
#
proc DataField::checkNumbers { values
                    frmt
                    message
                    {range ""}
                    {valueset ""}
                    {var_size 0}
                    {data_size 1}
                    {accept_scalar 1}
                    {min_max_sizes ""}
                  } {
  global Info
  upvar $message result_msg

  set result_msg ""

  #--Accepted column and row sizes
  #  -----------------------------
  #
  set accepted_d1 [lindex [lindex $data_size 0] end] ;# Nof columns <--> entrySize
  set accepted_d2 [lindex [lindex $data_size 1] end] ;# Nof rows    <--> nofEntries

  if { $accepted_d1 == "" } {
    set accepted_d1 1
  }

  if { $accepted_d2 == "" } {
    set accepted_d2 1
  }

  # If a table with variable values (entry values in one row, any number of rows)
  if { $var_size != 0 } {
    set accepted_d1 [expr $var_size + ($accepted_d1 * $accepted_d2)]
    set accepted_d2 "N"
  }

  #--Actual column and row sizes
  #  ---------------------------
  # Split values by possible list-separator (pick off the possible
  # last trailing first!)
  #
  set values [string trim $values $Info(arrayDataSeparator)]
  set values [split $values $Info(arrayDataSeparator)]

  #-Actual dim2: nof of separated entries
  # -----------
  set actual_d2 [llength $values]

  #-Evaluate first possible field expression
  set slot_value [lindex $values 0]
 
  if { [DataField::evaluateFieldValue $slot_value eval_value] } {
    set values [lreplace $values 0 0 $eval_value ]
  }

  #-Actual dim1: size of the first entry
  # -----------
  set actual_d1 [llength [lindex $values 0]]

  #-Possible minimum/maximum sizes
  #
  set min_sz1 [lindex $min_max_sizes 0]
  set max_sz1 [lindex $min_max_sizes 1]

  set min_sz2 [lindex $min_max_sizes 2]
  set max_sz2 [lindex $min_max_sizes 3]

  
  #--Start checking
  #  --------------
  set size_ok 0

  #--If scalar is ok
  if {$accept_scalar} {
    if { ( $actual_d1 == 1 && $actual_d2 == 1 )
         ||
         ( $actual_d1 == [expr $var_size + 1] && 
           ($actual_d2 == 1 || $accepted_d2 == "N")
         )
       } {
      set size_ok 1
    }
  }

  #-Check if a "flat" data matches expected dimensions with the stride
  #
  if { !$size_ok && $actual_d2 == 1 } {

    if { $min_sz1 != "" && $actual_d1 < $min_sz1 } {
      set result_msg "Nof values ($actual_d1) should be >= $min_sz1"
      return ""
    }

    if { $max_sz1 != "" && $actual_d1 > $max_sz1 } {
      set result_msg "Nof values ($actual_d1) should be <= $max_sz1"
      return ""
    }

    set data_size $actual_d1
    set multiple_info [list $data_size $accepted_d1 $accepted_d2 actual_d1 actual_d2]

    set result_msg [DataField::checkSizeParameter "Size for data" "" "" $multiple_info]

    if { $result_msg != "" } {
      return ""
    }

    set size_ok 1
  }


  #-Check that all entries are of equal size
  #
  if {!$size_ok} {

    if { $min_sz2 != "" && $actual_d2 < $min_sz2 } {
      set result_msg "Nof values per row ($actual_d2) should be >= $min_sz2"
      return ""
    }

    if { $max_sz2 != "" && $actual_d2 > $max_sz2 } {
      set result_msg "Nof values per row ($actual_d2) should be <= $max_sz2"
      return ""
    }

    for {set i 0 } {$i < $actual_d2} {incr i } {

      #-Evaluate first possible expression
      set slot_value [lindex $values $i]

      if { [DataField::evaluateFieldValue $slot_value eval_value] } {
        set values [lreplace $values 0 0 $eval_value ]
      }
      
      #-Actual value
      set actual_d1 [llength [lindex $values $i]]

      set index [expr 1 + $i]

      set result_msg [DataField::checkSizeParameter "Size for entry($index)" $accepted_d1 $actual_d1]

      if { $result_msg != "" } {
        return ""
      }
    }
  }

  #-Check that nof entries is ok (rows, size2)
  #
  if {!$size_ok} {

    set result_msg [DataField::checkSizeParameter "Nof entries " $accepted_d2 $actual_d2]

    if { $result_msg != "" } {
      return ""
    }
  }

  #--Format data
  #
  set values [join $values] ;# Flat data
  set data_length [llength $values] ;# Total data size
  
  #-Row separator flag
  if { $actual_d2 == 1 } {
    set add_separator 0
  } else {
    set add_separator 1
  }

  set start_pos 0
  set end_pos [expr $actual_d1 - 1]

  set counter 0

  set result ""

  #-Loop data and format (row by row)
  #
  while { $end_pos < $data_length } {

    #-Append separator after the first "row"
    if { $add_separator && $counter > 0 } {
      append result "$Info(arrayDataSeparator) "
    }

    #-Take one "row"
    set values1 [lrange $values $start_pos $end_pos]

    set result1 ""

    # Check each number
    # =================
    foreach value $values1 {
      set res_value [DataField::checkNumber $value $frmt result_msg $range $valueset]

      if {$result_msg != "" } {
        return ""
      }

      set res_value [string trim $res_value]
      append result1 $res_value
      append result1 " "
    }

    #-Append formatted row to the final result
    append result $result1

    incr start_pos $actual_d1
    incr end_pos $actual_d1
    incr counter
  }

  #--Ok 
  set dsize ""

  if { $var_size > 0 } {
    set dsize $data_size
    set nof_entries $actual_d2
  } else {
    lappend dsize $actual_d1
    lappend dsize $actual_d2
    set nof_entries 1
  }

  return [list $result $dsize $nof_entries]
}


# Evaluates a field expression if not simple number
# Returns: 1 <--> evaluated, 0 <--> not evaluated
#
proc DataField::evaluateFieldValue { value result_value } {
  upvar $result_value result

  set result $value

  if { [string is double $value] } {
    return 0
  }

  if { ![catch {set value [expr $value]} ] } {

    if { $value < 10000 } {

      set decs [expr round(log10($value))]
      set pwr [expr 6 - $decs]

      set value [expr pow(10,-$pwr) * round(pow(10,$pwr)*$value)]
    }

    set result $value
   
    return 1
  }
   
  return 0
}


#-Checks numbers (integers, floats etc)
#-Possible error msg is put into 'message' variable
#-Returns a formatted result.
proc DataField::checkNumber {value frmt message {range ""} {valueset ""} } {
  upvar $message msg

  set msg ""
  
  # Unformatted input value
  #
  set input_value $value

  if { ![catch {set test_value [expr $value]} ] } {

    set value $test_value

    # Use scientific format for large real numbers
    #
    if { ([string match -nocase "*g" $frmt] != 0) &&
         ([expr abs($value)] >= 10000)
     } {
      set frmt "% -12.4e"
    }
  }

  #-Make formatted result value and check also a 
  # possible incorrect number format
  #
  if { [catch {set result [format $frmt $value] } dummy] } {

    # Not correct!
    set msg $dummy
    return $value
  }


  
  # If input value was already formatted in e-format do no reformat!
  #
  # NOTE: E-format is reformatted, this way user can reformat an e-formatted!!!
  #
  if { -1 != [string first "e" $input_value] } {
    set result $input_value
  }

  #-Within range? []
  if {$range != ""} {

    set range [string trim $range]

    set min_open 0
    set max_open 0

    # Is interval is open at left
    if { "(" == [string index $range 0] } {
      set min_open 1
    }

    # Is interval is open at right
    if { ")" == [string index $range end] } {
      set max_open 1
    }

    # Remove possible interval markes
    set range [string trim $range "\(\)\[\]"]

    set min_val [ lindex $range 0 ]
    set max_val [ lindex $range 1 ]

    if { $min_val != "-Inf"} {
      if { $min_open && $result <= $min_val ||
           !$min_open && $result < $min_val
         } {
        set msg "Number should be larger than $min_val"
        return $result
      }
    } 

    if { $max_val != "Inf"} {
      if { $max_open && $result >= $max_val ||
           !$max_open && $result > $max_val
         } {
        set msg "Number should be smaller than $max_val"
        return $result
      }
    } 
  }

  #-Among correct valueset?
  if { $valueset != "" } {
    set ok 0

    foreach value $valueset {
      if { $result == $value } {
        set ok 1
        break
      }
    }

    if { $ok != 1 } {
      set msg "Number ($result) not among correct values ($valueset)!"
      return $result
    } 
  }

  #-Ok 
  return $result
}


# Format numbers (scalar and arries)
#
# Returns formatted result
#
proc DataField::formatNumbers {values frmt} {
  global Info

  # Data size
  set values [split $values $Info(arrayDataSeparator)]
  set d1 [llength [lindex $values 0]]
  set d2 [llength $values]

  set counter 0

  set result ""

  while { $counter < $d2 } {
  
    #-Append separator after the first "row"
    if { $counter > 0 } {
      append result "$Info(arrayDataSeparator)"
    }

    #-Take one "row"
    set values1 [lindex $values $counter]
    set result1 ""
    
    #-Format number in the "row"
    foreach value $values1 {
      #set res_value [format $frmt $value]
      set msg ""
      set res_value [DataField::checkNumber $value $frmt msg]
      set res_value [string trim $res_value]
      append result1 $res_value
      append result1 " "
    }

    #-Append formatted row to the final result
    append result [string trim $result1]

    incr counter
  }

  #-Ok 
  return $result
}


# Check that given size is ok
# return non-null error message if NOT ok
# If "multiple_info" is given, we accept a "flat" entry
# if size (stride)  and multiple-count match with the data.
# "multiple_info" should be a list of five items:
# 3 value arguments (numbers): total length, accepted-d1, accepted-d2,
# 2 reference arguments (names) for results: actuald1, actuald2
#
proc DataField::checkSizeParameter { size_name expected_size given_size { multiple_info ""} } {
  global Model

  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    set dim 3
  } else {
    set dim 2
  }

  # If expected size is based on:
  #-Dimension
  if { $expected_size == "DIM" } {
    set needed_size $dim
  #-Predefined value
  } else {
    set needed_size $expected_size
  }

  # If some multiple count accepted
  #
  if { $multiple_info != "" } {

    set data_length [lindex $multiple_info 0]

    set accepted_multiple1 [lindex $multiple_info 1]
    set accepted_multiple2 [lindex $multiple_info 2]

    upvar [lindex $multiple_info 3] actual_multiple1
    upvar [lindex $multiple_info 4] actual_multiple2

    # If accepted multiple is based on dim:
    if { $accepted_multiple1 == "DIM" } {
      set accepted_multiple1 $dim
    }
    if { $accepted_multiple2 == "DIM" } {
      set accepted_multiple2 $dim
    }
  }
  
  # Free index, ok
  #
  if { $expected_size == "N" } {
    return ""
  }

  # Check if we have an acceptable multiple count in data
  if { $multiple_info != "" } {

    set residual 0

    set multiple_msg "does not match expected size: ($accepted_multiple1 , $accepted_multiple2)"

    if { $accepted_multiple1 != "N" && $accepted_multiple1 != 0 } {
      set residual [expr $data_length % $accepted_multiple1]
      set actual_multiple1 $accepted_multiple1
      set actual_multiple2 [expr $data_length / $accepted_multiple1]

    } elseif { $accepted_multiple2 != "N" && $accepted_multiple2 != 0 } {
      set residual [expr $data_length % $accepted_multiple2]
      set actual_multiple1 [expr $data_length / $accepted_multiple2]
      set actual_multiple2 $accepted_multiple2
    }

    if { $residual != 0 } {
      return "Residual=$residual\n$multiple_msg" 
    }


    if { ($accepted_multiple1 != "N" && $accepted_multiple1 != $actual_multiple1) ||
         ($accepted_multiple2 != "N" && $accepted_multiple2 != $actual_multiple2)
       } { 
      return "Data size ($actual_multiple1 , $actual_multiple2):\n$multiple_msg"
    }

    return ""
  }

  #-If size is not ok
  if { $given_size != $needed_size } {
    set msg "$size_name: Expected ($needed_size) , Given ($given_size)"
    return $msg
  }

  return ""
}


#-Checks MATC script
#-NOTE: Field start with "$"
#-Possible error msg is put into 'message' variable
#-Format is not relevant here!
#
proc DataField::checkMatcScript {value frmt message} {
  global Info
  upvar $message msg

  set result ""
  set msg ""

  set value [string trim $value]
  
  Util::doMath [string range $value 1 end]

  return $result
}


#-Checks procedure-type field
#-NOTE: Currently only demanded that name starts with a letter!
#-Possible error msg is put into 'message' variable
#-Format is not relevatn here!
#
proc DataField::checkProcedureName {value frmt message} {
  global Info
  upvar $message msg

  set result ""
  set msg ""

  set value [string trim $value]
  #set value [string trim $value $Info(groupDataSeparator)]

  set values [DataField::splitProcedureEntryValue $value]

  set nof_values [llength $values]

  if { $nof_values < 2} {
    set msg "NOTE: At least library and function names must be given!\n"
    append msg "General format is:\nLIBRARY-FILE-name; FUNCTION-name; VARIABLE-name; SOURCE-name"
    return $result
  }

  if { $nof_values < 3 } {
    lappend values [list ""]
  }

  if { $nof_values < 4 } {
    lappend values [list ""]
  }

  if { $nof_values > 4} {
    set msg "too many strings in \[$values\] for a procedure name.\n"
    append msg "Should be:\nLIBRAY-FILE-name; FUNCTION-name; VARIABLE-name; SOURCE-name"
    return $result
  }

  set val [lindex $values 1]
  if { [string equal -nocase -length 5 "\$MATC" $val] } {
    set is_matc_proc 1
  } else {
    set is_matc_proc 0
  }

  set valuetypes { "Library" "Function" "Variable" "File" }
  
  if { $is_matc_proc } {  
    set needs {0 1 0 0}
  } else {
    set needs {1 1 0 0}
  }

  set checked_value ""

  foreach val $values type $valuetypes need $needs {
    
    set val [join $val]
    set val [string trim $val "\""]
    set val [string trimleft $val]

    # Library and function names must be given
    #
    if { $val == "" && $need } {
      set msg "ERROR: $type name missing!"
      return $result
    }

    # Check a non-empty value
    #
    if { $val != "" } {
      set msg ""

      set check_proc "check"
      append check_proc $type
      append check_proc "Name"

      set res_val [$check_proc $val "" msg]

      if { $msg != "" } {
        return $result
      }

    # Add empty value
    #
    } else {
      set res_val ""
    } 

    append checked_value $res_val
    append checked_value $Info(groupDataSeparator)
  }

  set checked_value [string trimright $checked_value $Info(groupDataSeparator)]


  #-Set into Library;File;... format
  set result [DataField::trimStringGroupData $checked_value]

  return $result
}


#-Checks filename-type field
#-NOTE: Currently only changes \ --> /
#-Possible error msg is put into 'message' variable
#-Returns a formatted result.
#
proc DataField::checkFileName {value frmt message} {
  global Info
  upvar $message msg

  set msg ""

  set val [Util::stringTrim $value]

  set val [string map { \\ / } $val]

  # Check that no quotes inside the name!  
  set val [string trim $value "\""]

  if { -1 != [string first "\"" $val] } { 
    set msg "Illegal file name (quote inside name): $val"
  }
  
  #-Format not in real use!
  set result $val

  return $result
}


#-Checks directoryname-type field
#-NOTE: Currently only changes \ --> /
#-Possible error msg is put into 'message' variable
#-Returns a formatted result.
#
proc DataField::checkDirectoryName {value frmt message} {
  global Info
  upvar $message msg

  set msg ""

  set val [Util:trimString $value]

  set val [string map { \\ / } $val]

  # Check that no quotes inside the name!  
  set val [string trim $value "\""]

  if { -1 != [string first "\"" $val] } { 
    set msg "Illegal directory name (quote inside name): $val"
  }
  
  #-Format not in real use!
  set result $val

  return $result
}


# Quotes inside a string are not accepted!
#
proc DataField::checkString {value frmt message} {
  global Info
  upvar $message msg

  set value [Util::stringTrim $value]
    
  set msg ""

  # Check that no quotes inside the string!  
  set val [string trim $value "\""]

  if { -1 != [string first "\"" $val] } { 
    set msg "Illegal string value (quote inside string): $val"
  }

  return $value
}


#-Checks function name value
#-NOTE: Currently only demanded that name starts with a letter!
#-Possible error msg is put into 'message' variable
#-Returns a formatted result.
#
proc DataField::checkFunctionName {value frmt message} {
  global Info
  upvar $message msg

  set msg ""
  
  set val [Util::stringTrim $value]

  if { $val == "" } {
    return $val
  }

  #-A Matc function definition, it must be have been earlier!
  if { [string equal -nocase -length 5 "\$MATC" $val] } {
    return $val
  }

  #-Must start with a letter
  if { ![regexp {[a-zA-Z$]} [string index $val 0]] } {
    set msg "Invalid function name: Must start with a letter!"
    return $val
  }

  if { [regexp {[^a-zA-Z0-9_$]} $val] } {
    set msg "Invalid function name: Contains illegal characters or spaces!"
    return $val
  }

  return $val
}


#-Checks library name value
#-NOTE: No checking, value just trimmed!
# MVe 18.11.99
#-Possible error msg is put into 'message' variable
#-Returns a formatted result.
#
proc DataField::checkLibraryName {value frmt message} {
  global Info
  upvar $message msg

  set msg ""

  #-Trim input value
  set val [Util::stringTrim $value]

  return $val
}


#-Checks variable name (in procedure definition) value
#-Possible error msg is put into 'message' variable
#-Returns a formatted result.
#
proc DataField::checkVariableName {value frmt message} {
  global Info
  upvar $message msg

  set msg ""

  set val [Util::stringTrim $value]

  #-Must start with a letter
  if { ![regexp {^[a-zA-Z]+.*} $val] } {
    set msg "$val is not a proper variable name"
  }

  if { "none" == [string tolower $val] } {
    set val ""
  }

  #-Format not in real use!
  set result $val

  return $result
}


proc DataField::splitProcedureEntryValue {value} {
  global Info

  #-Data given in "Library" "File" format
  if { "\"" == [string index $value 0] } {
    set values [split $value]

  #-Data given in Library;File format
  } else  {
    set values [split $value $Info(groupDataSeparator)]
  }

  return $values
}


# -Trims leading and trailing blanks from "list of strings"
#  List separator is the character in the variable Info(groupDataSeparator)
#-Returns a formatted result.
#
proc DataField::trimStringGroupData {value} {
  global Info

  set values [split $value $Info(groupDataSeparator)]
  set result ""
  set list_count [llength $values]

  set counter 0
  foreach val $values {
    incr counter
    set val [string trimright [string trimleft $val]]
    append result $val
    if {$counter < $list_count} {
      append result $Info(groupDataSeparator)
    }
  }

  return $result
}


# Checks floating point-type field
# NOTE: 10E-2 etc. currently not accepted!
# Uses regular expressions
proc DataField::checkFloat {value} {

  #-(Optional +,-)(123 or 123.123)(Must end with number)
  set code [regexp {^[+-]?([0-9]+|[0-9]*[.][0-9]+)$} $value]
  if {0 == $code} {
    return "not a number"
  }

  return ""
}


# Checks integer-type field
# Uses regular expressions
proc DataField::checkInteger {value} {

  #-(Optional +,-)(1213)(Must end with number)
  set code [regexp {^[+-]?[0-9]+$} $value]
  if {0 == $code} {
    return "not a number"
  }

  return ""
}



################################
#### Fields <--> Data lines ####
################################


# Get field value (like TEMPERATURE)
# Reference argument "active_var" tells if field is active
# Description is returned in the reference argument "desc_var" if that
# is given
# Returns field's value part
#
proc DataField::getFieldValue {globArray id var {active_var "" } {old_value 0 } {desc_var ""} } {
  global Info
  upvar #0 $globArray theArray

  if { $active_var != "" } {
    upvar $active_var active
  }
  set active 0

  
  if { $desc_var != "" } {
    upvar $desc_var desc
  }
  set desc ""
  
  #-Normal field
  set pattern1 $var
  append pattern1 "=="

  #-File name
  set pattern2 $var
  append pattern2 "=:"

  #-Proc name
  set pattern3 $var
  append pattern3 "=."

  set field_data ""

  if { $old_value } {

    if { ![info exists theArray($id,data,old)] } {
      return ""
    } else {
      set fields [split $theArray($id,data,old) $Info(dispListSeparatorIn)]
    }

  } else {

    if { ![info exists theArray($id,data)] } {
      return ""
    } else {
      set fields [split $theArray($id,data) $Info(dispListSeparatorIn)]
    }
  }

  foreach fld $fields {

    if { $Info(inactiveMark) == [string index $fld 0] } { 
      set active 0
    } else {
      set active 1
    }

    set fld [string trimleft $fld $Info(inactiveMark)]

    if { [regexp $pattern1 $fld] } {
      regsub $pattern1 $fld "" field_data
      break

    } elseif {  [regexp $pattern2 $fld] } {
      regsub $pattern2 $fld "" field_data
      break

    } elseif {  [regexp $pattern3 $fld] } {
      regsub $pattern3 $fld "" field_data
      break
    }
  }

  set result ""
  set desc [DataField::extractVariableAndSizeDescription result $field_data]

  return [string trim $result $Info(dataListSeparator)] 
}


# Set new value for a field, replace also with a new description
# if given
# NOTE: Works only for active fields if no inactive marker set in front
# of the argument "fld_name" (like -TEMPERATURE)
# Does not return any value
#
proc DataField::setFieldValue {globArray id field_name new_value { new_desc ""} } {
  global Info
  upvar #0 $globArray theArray

  #--This will contain everything before data field
  #  itself
  set pattern $field_name
  append pattern "=="

  #--Split data into fields
  set fields [split $theArray($id,data) $Info(dispListSeparatorIn)]

  set index 0
  set found 0

  #--Find old field value
  foreach fld $fields {

    if { [regexp $pattern $fld] } {
      regsub $pattern $fld "" old_field_data
      set found 1
      break
    }

    incr index 
  }

  #--An existing field
  #
  if {$found} {

    set new_field $pattern

    #-Pick old description (we are not interested in the old data here,
    # but procedure needs this argument!
    #
    set old_value ""
    set old_desc [DataField::extractVariableAndSizeDescription old_value $old_field_data]

    #-Set possible description
    if { $new_desc != "" } {
      append new_field $new_desc

    } else {    
      append new_field $old_desc
    }
    
    append new_field $new_value

    #--Replace with the updated field
    set fields [lreplace $fields $index $index $new_field]

  #--New field
  } else {
    set idx [lsearch $theArray(allFields) $field_name]

    if { $idx == -1 } {
      return
    }

    # Very simple currently!!!
    # MVe 23.09.00.
    set new_field "$field_name==$new_value"
    set fields [linsert $fields $idx $new_field]

  }
#MSG "fields=$fields"
  #--Update parameter
  set theArray($id,data) [join $fields $Info(dispListSeparatorIn)]
}


#-Procedure formats a parameter id-number for display
#
proc DataField::formatParamId {globArray pid} {
  upvar #0 $globArray theArray
  global Info
#  set param_format "%0$theArray(parameterIdLen)d"
  set param_format "%d"
  return [format $param_format $pid]
}


# Procedure forms a parameter for a "non-standard" panel like
# Coordinate, Constant, Solver ... from current field values
#
proc DataField::formNonStandardParameter { globArray pid parameter_name } {
  global Info
  upvar #0 $globArray theArray

  set DS $Info(dispListSeparatorOut)
  set FS $Info(fieldSeparator)

  #-Name, Body1Id, Body2Id
  set theArray($pid,name) $parameter_name
  set theArray($pid,bd1) $Info(NO_INDEX)
  set theArray($pid,bd2) $Info(NO_INDEX)


  #-Data 
  set data [DataField::formParameterLine $globArray]

  set theArray($pid,data) $data
}


# Procedure forms a parameter for a "non-standard" panel like
# Coordinate, Constant, Solver ... from initial values
#
proc DataField::formDefaultNonStandardParameter { globArray pid } {
  global Info
  upvar #0 $globArray theArray

  #-Name, Body1Id, Body2Id and data
  set theArray($pid,name) $globArray$pid
  set theArray($pid,bd1)  $Info(NO_INDEX)
  set theArray($pid,bd2)  $Info(NO_INDEX)
  set theArray($pid,data) [Panel::createDefaultParameterLine $globArray]
}


# Check parameter data against given problem mask
#
# This is used ex. when equation definitions have been changed
# and all parameters must be checked against the new problem mask!
#
# Removes obsolate fields and initializes added fields
# Argument "data_modified" is a reference variable to indicate if data was modifeid
#
# NOTE: This takes time!!!
#
proc DataField::checkParameterData {globArray pid data_modified {mask ""} {system_index ""} {force_initial 0} } {
  global Info 
  upvar #0 $globArray theArray
  upvar $data_modified modified

  set modified 0

  set old_data $theArray($pid,data)

  set restore_tmask 0

  if { [info exist theArray(targetMask)] } {
    set tmask $theArray(targetMask)
    set restore_tmask 1
  }

  # Set array variable
  set theArray(targetMask) $mask

  Panel::unsetAllFieldData $globArray

  Panel::resetFields $globArray

  Panel::initFields $globArray $theArray(allFields) $force_initial $system_index

  DataField::formDataFields $globArray $old_data $mask

  PanelCheck::execPanelFillProcs $globArray

  set new_data [DataField::formParameterLine $globArray $system_index]


  # WARNING: This simple test makes Model as updated(*)
  # nearly always, should be better!!!
  #
  if { $new_data != $old_data } {
    set modified 1
  }

  # NOTE: This was commented because when Apply is done
  # on the Equation panel, then ex. VARIABLE_LIST,act is
  # deleted and check procs cannot work any more!!!###!!!
  # MVe 20.11.02
  #Panel::unsetAllFieldData $globArray

  if { $restore_tmask } {
    set theArray(targetMask) $tmask
  }

  return $new_data
}


# Check if this field has activity parents whose value should
# make this field active! (Like "Diffuse gray" RADIATION for EMISSIVITY etc.)
#
# Return:
#  1 <--> activated
#  0 <--> do not know
# -1 <--> deactivated
#
# NOTE: If a field has multiple activity parents, this proc should should NOT normlly
# be called as "parent_flds" given, because then the order of parents could decide the
# activity decision. So, normally this argument should be left blank and let the proc use
# all parents at the same time!!!
#
proc DataField::isActivatedByParents { globArray fld eq_index param_ids problem_mask { parent_flds ""} } {
  global Info
  upvar #0 $globArray theArray 

#MSG "globArray=$globArray field=$fld eq-index=$eq_index param-ids=$param_ids"

  set activity_parents [DataField::getFieldProperty $globArray $fld ActivityParents $eq_index]

  #--No activity info defined
  if { $activity_parents == "" } {
    return 0
  }

  if { $parent_flds == "" } {
    set parent_infos $activity_parents

  } else {
    set parent_infos $parent_flds
  }

  set activity_values [DataField::getFieldProperty $globArray $fld ActivityValues $eq_index]

  set some_active 0
  set some_nonactive 0
  
  # Check all parents
  # -----------------
  #
  foreach parent_info $parent_infos values $activity_values {

    set rc [DataField::isActivatedByParent_1 $globArray $fld $eq_index $param_ids $problem_mask $parent_info $values]
    
    #-At least one is  non-active
    #
    if { $rc == -1 } {
      set some_nonactive 1
    }
    
    #-At least one is active
    #
    if { $rc == 1 } {
      set some_active 1
    }
  }
  
  #-One active parent 'wins'
  #
  if { $some_active } {
    return 1
  
  #-Next one non-active parent 'wins'
  #
  } elseif {$some_nonactive } {
    return -1

  #-Cannot decide!
  #
  } else {
    return 0
  }
}


# Check if this field has an activity parent whose value should
# make this field active! (Like "Diffuse gray" RADIATION for EMISSIVITY etc.)
#
# NOTE: The parent fld and activity values given!!!
#
# Return:
#  1 <--> activated
#  0 <--> do not know
# -1 <--> deactivated
#
proc DataField::isActivatedByParent_1 { globArray fld eq_index param_ids problem_mask parent_info activity_values} {
  global Info Model

  # Pick parent array and field name
  #
  if { 1 < [llength $parent_info] } {
    set parent_arr [lindex $parent_info 0]
    set parent_fld [lindex $parent_info 1]
    set own_panel 0
  } else {
    set parent_arr $globArray
    set parent_fld [lindex $parent_info 0]
    set own_panel 1
  }

  upvar #0 $parent_arr theArr 

  set parent_value ""
  set pvt ""
  
  # Special 'panels'
  # ----------------
  if { [string equal -nocase $parent_arr "Model"] } {
    
    if { $parent_fld == "GEOMETRY_DIMENSION" ||
         $parent_fld == "SIMULATION_DIMENSION"
       } {
      set parent_value [string index $Model($parent_fld) 0]
      set pvt "INTEGER"
    }

  # Field's own panel
  # -----------------
  } elseif { $own_panel } {
    
    # If parent field is not active!
    if { [info exists theArr($parent_fld,act)] && 
         !$theArr($parent_fld,act)
       } {
      return -1
    } else {
      set parent_value $theArr($parent_fld)
    }

  # Other panel (= Equation)
  # ------------------------
  } else {

    foreach pid $param_ids {

      set parent_value [DataField::getFieldValue $parent_arr $pid $parent_fld]
    }
  }

  set values $activity_values

  #--If no criteria values defined, just the activity of the parent
  #  decides first!
  #
  if { $values == "" } {
    if { ![DataField::fieldProblemMaskMatches $parent_arr $parent_fld $problem_mask] } {
      return -1
    }
  }

  #--If no value for the parent (yet) available!
  #
  if { $own_panel && ![info exists theArr($parent_fld)] } {
    return 0
  }
  
  if { $pvt == "" } {
    set pvt [DataField::getFieldProperty $parent_arr $parent_fld FieldValueType]
  }
  
  # Set default activity value
  #
  if { $activity_values == "" } {

    if { [string equal -nocase $pvt "LOGICAL"] } {
      set values True
    } elseif { [string equal -nocase $pvt "INTEGER"] } {
      set values ~0
    } elseif { [string equal -nocase $pvt "REAL"] } {
      set values ~0.0
    } else {
      set values ~
    }
  }

#MSG "fld=$fld parent-arr=$parent_arr parent-fld=$parent_fld values=$values"

  #--Loop each criteria value and check if field could be active
  #
  foreach value $values {
          
    #-If criterium value is a NEGATION (eg. for COMPRESSIBILITY_MODEL)!
    #
    if { $Info(negationMark) == [string index $value 0] } {

      set value [string range $value 1 end]
      
      # Numeric value
      #
      if { [string equal -nocase $pvt "LOGICAL"] || 
           [string equal -nocase $pvt INTEGER"]  || 
           [string equal -nocase $pvt "REAL"]
         } {

        if { $parent_value != $value } {
          return 1
        }

      # String
      #
      } else {
        if { ![string equal $parent_value $value] } {
          return 1
        }
      }

    #-Normal, POSITIVE criterium value
    #
    } else {
    
      # Numeric value
      #
      if { [string equal -nocase $pvt "LOGICAL"] || 
           [string equal -nocase $pvt INTEGER"]  || 
           [string equal -nocase $pvt "REAL"]
         } {
        if { $parent_value == $value } {
          return 1
        }

      # String
      #
      } else {
        if { [string equal $parent_value $value] } {
          return 1
        }
      }
    }
  }

  return -1
}



# Check if this field has activity slaves whose activity should be set based
# on current field value
# NOTE: This proc is launched typically from field check procs
#
#
proc DataField::checkActivitySlaves { globArray fld problem_mask} {
  upvar #0 $globArray theArray

  if { [info exists theArray(equationIndex)] } {
    set eq_index $theArray(equationIndex)
  } else {
    set eq_index ""
  }
  
  if { [info exists theArray(objectId)] } {
    set eq_ids [Object::getEquationIds $theArray(objectId)]
  } else {
    set eq_ids ""
  }

  set slave_infos [DataField::getFieldProperty $globArray $fld ActivitySlaves $eq_index]

  #--No activity info defined
  if { $slave_infos == "" } {
    return
  }


  foreach slave_info $slave_infos {
    
    if { 1 < [llength $slave_info] } {
      set slave_arr [lindex $slave_info 0]  
      set slave_fld [lindex $slave_info 1]
    } else {
      set slave_arr $globArray  
      set slave_fld [lindex $slave_info 0]
    }

    upvar #0 $slave_arr theArr

    set val [DataField::isActivatedByParents $slave_arr $slave_fld $eq_index $eq_ids $problem_mask]

    if { $val == 1 } {
      set theArr($slave_fld,act) 1
    } elseif { $val == -1 } {
      set theArr($slave_fld,act) 0
    }

    Widget::configureField $slave_arr $slave_fld
    
    # Recursive call for the slave!
    #
    DataField::checkActivitySlaves $slave_arr $slave_fld $problem_mask
  }
}


proc DataField::addActivityParent {globArray fld parent_info {eq_index ""} } {

  set parent_infos [DataField::getFieldProperty $globArray $fld ActivityParents $eq_index]

  # Pick new parent-info array and field
  #
  if { 1 < [llength $parent_info] } {
    set parent_info_arr [lindex $parent_info 0]
    set parent_info_fld [lindex $parent_info 1]
  } else {
    set parent_info_arr $globArray
    set parent_info_fld [lindex $parent_info 0]
  }

  # Check that parent-info does not already exist
  #
  foreach pi $parent_infos {

    if { 1 < [llength $pi] } {
      set pi_arr [lindex $pi 0]
      set pi_fld [lindex $pi 1]
    } else {
      set pi_arr $globArray
      set pi_fld [lindex $pi 0]
    }

    if { $pi_arr == $parent_info_arr &&
         $pi_fld == $parent_info_fld
       } {
      return
    }
  }
  
  # Append new info to the old list and store new list
  #
  lappend parent_infos $parent_info
  DataField::setFieldProperty $globArray $fld ActivityParents $parent_infos $eq_index
}


proc DataField::addActivitySlave {globArray parent slave_info {eq_index ""} } {

  set slave_infos [DataField::getFieldProperty $globArray $parent ActivitySlaves $eq_index]

  # Pick new slave-info array and field
  #
  if { 1 < [llength $slave_info] } {
    set slave_info_arr [lindex $slave_info 0]
    set slave_info_fld [lindex $slave_info 1]
  } else {
    set slave_info_arr $globArray
    set slave_info_fld [lindex $slave_info 0]
  }

  # Check that slave-info does not already exist
  #
  foreach si $slave_infos {

    if { 1 < [llength $si] } {
      set si_arr [lindex $si 0]
      set si_fld [lindex $si 1]
    } else {
      set si_arr $globArray
      set si_fld [lindex $si 0]
    }

    if { $si_arr == $slave_info_arr &&
         $si_fld == $slave_info_fld
       } {
      return
    }
  }
  
  # Append new info to the old list and store new list
  #
  lappend slave_infos $slave_info
  DataField::setFieldProperty $globArray $parent ActivitySlaves $slave_infos $eq_index
}


# Remove non-input fields from input parameter line: FieldClass OutputOnly
#
proc DataField::compressParameterInputData {globArray pid vnames_name} {
  global Info 
  upvar #0 $globArray theArray

  if { ![info exists theArray($vnames_name)] } {
    return
  }

  #---Split data row into a list using separators like '|'
  set sep $Info(dispListSeparatorIn)
  set val_list [split $theArray($pid,data) $sep]


  set new_data "" 

  # Loop all data fields
  # ====================
  foreach val $val_list {

    # Remove possible inactive indicator
    set clean_val [string trimleft $val $Info(inactiveMark)]

    #-Pick variable name and value part
    set tmp [split $clean_val "="]

    #-Field name
    set fld [lindex $tmp 0]

    #-Variable type
    if { ![Panel::isInputField $globArray $fld] } {
      continue
    } else {
      append new_data $val$sep
    }
  }

  #-Final data
  #
  if { $new_data != "" } {
    set new_data [string trimright $new_data $sep]
  }

if { $globArray == "Solver" } {
#MSG "new-data=$new_data"
}

   return $new_data
}


# Remove unneeded fields from output parameter line: empty, False and None values
# This is before saving the parameter
# Argument "data_modified" is a reference variable to indicate if data was modifeid
#
proc DataField::compressParameterOutputData {globArray pid data_modified {system_index ""} } {
  global Info 
  upvar #0 $globArray theArray
  upvar $data_modified modified

  set modifed 0

  #---Split data row into a list using separators like '|'
  set sep $Info(dispListSeparatorIn)
  set val_list [split $theArray($pid,data) $sep]


  set new_data "" 

  #set count [expr [llength $val_list] - 1]
  set count [llength $val_list]
  set index  0
  set index2 1


  # Loop all data fields
  # ====================
  while { $index < $count } {

    set val [lindex $val_list $index]

    # Remove possible inactive indicator
    if { [string index $val 0] == $Info(inactiveMark) } {
      set clean_val [string trimleft $val $Info(inactiveMark)]
    } else {
      set clean_val $val
    }

    #-Pick variable name and value part
    set tmp [split $clean_val "="]

    #-Field name
    set val_var [lindex $tmp 0]
    
    #-If field is not to be output
    if { ![Panel::isOutputField $globArray $val_var] } {
      set modified 1
      incr index
      incr index2
      continue
    }

    #-Field value
    # This for normal data (== value flag)
    if { "" == [lindex $tmp 1] } {
      set val_val [lindex $tmp 2]

    # This is for file and proc (=: or =. value flag)
    } else {
      set val_val [string range [lindex $tmp 1] 1 end]
    }

    set val_val [string trim [join $val_val] $Info(dataListSeparator)]

    set is_single_entry 1

    #-Peak next variable name, if it is the same as the
    # current, then an empty field could be a "blocking empty"
    # (see below)
    #
    if { $index2 < $count } {

      set val2 [lindex $val_list $index2]

      # Remove possible inactive indicator
      if { [string index $val2 0] == $Info(inactiveMark) } {
        set clean_val2 [string trimleft $val2 $Info(inactiveMark)]
      } else {
        set clean_val2 $val2
      }

      #-Pick variable name and value part
      set tmp2 [split $clean_val2 "="]

      #-Field name
      set val_var2 [lindex $tmp2 0]

      if { $val_var2 == $val_var } {
        set is_single_entry 0
      }
    }

    #-Increment variable indices
    incr index
    incr index2

    # Should we drop the field because it is in initial value
    # (Like Incompressible for COMPRESSIBILITY)
    #
    set drop_initial_valued 0

    set ival [DataField::getFieldProperty $globArray $val_var InitialValue $system_index]

    if { [DataField::getFieldProperty $globArray $val_var DropInitialValued] } {

      if { $val_val == $ival } {
        set drop_initial_valued 1
      }
    }

    #-Find if data-var can be dropped
    #
    #-NOTE: We do not drop "empty" entry unless it is the only
    # entry for the variable in the data line, because we must
    # block "inactive" proc and table entries to become
    # active by accident (as next non-empty instances when a
    # parameter is reactivated!)
    # However, we never drop the ACTIVE-field!
    # NOTE: We do not drop an empty single value if its initial
    # value is non-empty, because it is now blocking the initial
    # value to be set on automatically!
    #
    # Drop value
    if { $val_var != "ACTIVE" &&
         $is_single_entry &&
         ($drop_initial_valued ||
          ($val_val == "" && $ival == "")
         )
       } {
      set modified 1
      continue
    }

    #-Add the variable to the result
    #
    append new_data $val$sep
  }

  #-Final data
  if { $new_data == "" } {
    return ""

  } else {
    set return [string trimright $new_data $sep]
  }
}


#-Procedure reads data from a parameter line and
# fills data area fields with the data
#
proc DataField::setDataFields {globArray data_line} {
  upvar #0 $globArray theArray

  Panel::resetFields $globArray

  #--Update data-field variables
  StdPanelInit::initFieldValues $globArray

  DataField::formDataFields $globArray $data_line
}


proc DataField::clearLBSelection {list_box {selection_row "" } } {

  if {$selection_row != ""} {
    $list_box selection clear $selection_row $selection_row
  } else {
    $list_box selection clear 0 end
  }
}


proc DataField::setAndShowLBSelection {list_box selection_row} {

  if {$selection_row != ""} {
    $list_box selection set $selection_row 
    $list_box see $selection_row 
  }
}


#########################################
#########################################
### A. Data fields --> Parameter line ### 
#########################################
#########################################

# Procedure forms a parameter dataline from globArray's
# actice variables, variables-list variable name given as
# argument
# If problem_mask != "" it is used to drop variables which
# don't match the mask.
# NOTE: Only active fields are picked.
# NOTE: No parameter data update actually done!!!
#       This proc is used be only to check
#       if the parameter line would contain any active data
#
proc DataField::formActiveParameterLine { globArray vnames_var {problem_mask ""} } {
  global Info Common
  upvar #0 $globArray theArray

  set sep $Info(dispListSeparatorOut)

  set parameter_line ""
  foreach var $theArray($vnames_var) {

    #---Skip non-data entries and non-active entries
    if { ![Panel::isOutputField $globArray $var] ||
         !$theArray($var,act)
       } {
      continue
    }

    #---If problem mask given, check that variable matches the mask
    if { $problem_mask != "" } {

      # If mask does not match
      if { ![DataField::fieldTargetMaskMatches $globArray $var $problem_mask] } {
          continue
      }
    }

    #--If field in initial value would be inactive in output
    if { [DataField::getFieldProperty $globArray $var DropInitialValued] } {

      set initial_value [DataField::getInitialValue $globArray $var]

      if { $initial_value == $theArray($var) } {
        continue
      }
    }

    set fld_data [DataField::formFieldParameterLineData $globArray $var]

    if { $fld_data == "" } {
      continue
    }

    # Like: TEMP=100|
    append fld_data $sep


    append parameter_line $fld_data
  }

  set parameter_line [string trimright $parameter_line $sep]

  return $parameter_line
}


# Procedure forms a parameter dataline from globArray's
# field variables
# Inclusion of the inactive fields is controlled by the
# global flag variable Info(keepInactiveFields) (0/1)
#   
proc DataField::formParameterLine {globArray {system_index ""} } {
  global Info Common
  upvar #0 $globArray theArray

  set sep $Info(dispListSeparatorOut)

  set parameter_line ""

  foreach fld $theArray(allFields) {

    #---Skip non-data entries 
    if { ![Panel::isOutputField $globArray $fld] } {
      continue
    }

    #---If inactives are not kept, skip field
    #
    if { !$theArray($fld,act) &&
         !$Info(keepInactiveFields)
       } {
      continue
    }

    set fld_data [DataField::formFieldParameterLineData $globArray $fld $system_index]

    if { $fld_data == "" } {
      continue
    }

    # Like: TEMP=100|
    append fld_data $sep

    append parameter_line $fld_data
  }

  set parameter_line [string trimright $parameter_line $sep]

  return $parameter_line
}


proc DataField::formFieldParameterLineData {globArray fld {system_index ""} } {
  global Info

  upvar #0 $globArray theArray 

  #-Data separators
  set sep $Info(dispListSeparatorOut)
  set extra1 $Info(dispListSeparatorExtra1)
  set extra2 $Info(dispListSeparatorExtra2)

  #-Data item types:
  # Normal value like: TEMP==123.12                         (==)
  # Procedure name like: TEMP=."Temp_modules" "Temp1_proc"  (=.)
  # File name like: TEMP=:Filename.txt                      (=:)

  set flags ""
  set datas ""
  
  #---Ddata type flags
  
  #-Current type

  set is_value 1
  set is_proc 0
  set is_table 0

  if { [info exists theArray($fld,proc)] && $theArray($fld,proc) } {
    set is_proc 1
  }

  if { [info exists theArray($fld,table)] && $theArray($fld,table) } {
    set is_table 1
  }

  #-Existing types
  set iv [DataField::getInitialValue $globArray $fld $system_index]
  set value_exists 0
  set proc_exists 0
  set table_exists 0

  if { [info exists theArray($fld,valueData)] && $theArray($fld,valueData) != "" } {
    set value_exists 1
  }

  if { [info exists theArray($fld,procData)] && $theArray($fld,procData) != "" } {
    set proc_exists 1
  }

  if { [info exists theArray($fld,tableData)] && $theArray($fld,tableData) != "" } {
    set table_exists 1
  }

  #---Handle active data type
  #-Procedure
  if { $is_proc } {
    lappend flags $Info(procFlag)
    lappend datas procData
 
  #-Table
  } elseif { $is_table } {
    lappend flags $Info(valueFlag)
    lappend datas tableData

  #-Value data (normal scalar value)
  # NOTE: An empty scalar (value) field is removed only if its initial value is empty and
  # no other value types are present, otherwise it is kept as overwritting an initial value 
  # of for blocking the other types to become automatically active as
  # the next non-empty instance!
  } else {
    if { $theArray($fld) != "" || $iv != "" || $proc_exists || $table_exists } {

      set theArray($fld,valueData) $theArray($fld)
      set is_value 1
      set value_exists 1

      lappend flags $Info(valueFlag)
      lappend datas valueData
    }
  }

  #---If field value is NOT for procedure, but proc data exists, save it also
  if { !$is_proc  && $proc_exists } {
    lappend flags $Info(procFlag)
    lappend datas procData
  }

  #---If field value is NOT for table, but table data exists, save it also
  if { !$is_table && $table_exists } {
    lappend flags $Info(tableFlag)
    lappend datas tableData
  }

  #---If normal value exists, but is not current, store it also
  if { !$is_value && $value_exists } {
    lappend flags $Info(valueFlag)
    lappend datas valueData
  }

  # Skip completely empty fields
  if { $datas == "" } {
    return ""
  }

  set values_list ""

  # Form all field ids and values into separate lists
  # =================================================
  foreach flag $flags data $datas {

    # Field ids
    # =========
    if { $theArray($fld,act) } {
      set fld_id $fld$extra1$flag$extra2

    } else {
      set fld_id $Info(inactiveMark)$fld$extra1$flag$extra2
    }

    lappend fld_ids $fld_id

    # Form var name (for getting the actual data)
    if { $data != "" } {
      set vn $fld,$data
    } else {
      set vn $fld
    }  
  
    set is_prc 0
    set is_tbl 0

    if { $data == "procData" } {
      set is_prc 1
    } elseif { $data == "tableData" } {
      set is_tbl 1
    }

    # Field values
    # ============
    #-Get possible description
    set fld_desc [DataField::formVariableAndSizeDescription $globArray $vn]

    #-Group-type procedure name is picked from the first entry field
    if { [Panel::isGroupField $globArray $fld] } {

      set fld_val [DataField::formGroupParameterLineValue $globArray $fld $is_prc]

    } else {
      set fld_val [DataField::formSingleParameterLineValue $globArray $fld $vn 0]
    }

    append values_list $fld_val

    # Empty field
    # ------------------
    if { $fld_val == "" } {
      lappend fld_values ""
 
    # Data containing field (active or storable inactive field)
    # ---------------------
    } else {
      lappend fld_values $fld_desc$fld_val
    }
  }
  
  # Not any data values for the field!
  #
  if { $values_list == "" } {
    return ""
  }

  # Form final field data
  # ======================
  set fld_data ""

  foreach fld_id $fld_ids val $fld_values {
    append fld_data $fld_id$val$sep
  }

  set fld_data [string trimright $fld_data $sep]

  return $fld_data
}


# Form dataline (parameter) representation of a single field
# NOTE: var_item is the actual array item like TEMPERATURE,table
#
proc DataField::formSingleParameterLineValue {globArray var var_item isFile} {
  global Info
  upvar #0 $globArray theArray

  if { ![info exists theArray($var_item)] } {
    return ""
  }

  # Normal field
  # ============
  if { !$isFile } {
    return $theArray($var_item)
  }

  # File name field
  # ===============

  set sep $Info(dataListSeparator) 
  set result ""

  set values [split $theArray($var_item) $sep]

  foreach val $values {

    set val [string trim $val "\""]
    append result $val
    append result $sep
  }

  set result [string trimright $result $sep]

  return $result
}


proc DataField::formGroupParameterLineValue {globArray var isProcedure} {
  upvar #0 $globArray theArray
  global Info Common
  set panel $theArray(parameterType)

  set glist $Common($panel,$var,group)

  #-Procedure name. It is read from the first group entry-field
  if {$isProcedure} {
    set result $theArray([lindex $glist 0])
    return $result
  }

  #-Normal values. If all fields are empty, "result" is set to an empty list
  set hasData 0

  #-Radio-button values is stored into the group parent itself
  if { [Panel::isRadioButtonField $globArray $var] } {
    set fld $theArray($var)
    if {$fld != ""} {
      set hasData 1
      lappend result $fld
    }

  #-Otherwise read all group member widgets and make a data list
  } else {
    foreach v $glist {
      if { ![Panel::isOutputField $globArray $v] } {
        continue
      }
      set fld [string trim $theArray($v)]
      if {$fld != ""} {
        set hasData 1
      }
      lappend result $fld
    }
  }

  if {$hasData} {
    set result [join $result $Info(groupDataSeparator)]
    #--Put separator also at the end of the data
    append result $Info(groupDataSeparator)
  } else {
    set result ""
  }

  return $result
}


# Form field's default (=initial value) contribution to a paramter
#
proc DataField::formInitialParameterLineValue {globArray fld {system_index ""} } {
  global $globArray Info
  upvar #0 $globArray theArray

  #-Data separators
  set sep $Info(dispListSeparatorOut)
  set extra1 $Info(dispListSeparatorExtra1)
  set extra2 $Info(dispListSeparatorExtra2)

  #-Initial value
  set val [DataField::getInitialValue $globArray $fld $system_index]
  set act [DataField::getFieldProperty $globArray $fld InitiallyActive $system_index]

  # Set possible inactive marker for the field
  if {!$act} {
    set act_mark $Info(inactiveMark)
  } else {
    set act_mark ""
  }

  #-Type
  if { "File" == [DataField::getFieldProperty $globArray $fld FieldValueType] } {
    set flag $Info(fileFlag)
  } else {
    set flag $Info(valueFlag)
  }

  #-Result
  return $act_mark$fld$extra1$flag$extra2$val
}


# Calculates the size of a data field
# NOTE: Currently just calculates the list length
#       This need changes when we have table!!!
#
#proc calcFieldSize {globArray fld} {
#  global Info Common
#  upvar #0 $globArray theArray
#  return [llength $theArray($fld)]
#}

# Ok this should work for arries and tables etc
#
proc DataField::calcFieldSize {field_data} {
  global Info

  set values [split $field_data $Info(arrayDataSeparator)]

  #-Column size from first row
  set dim1 [llength [lindex $values 0] ]

  #-Row size 
  # NOTE: Second index in Solver Input size-field!!!
  set dim2 [llength $values]

  if { $dim2 > 1 } {
    set result "$dim1 $dim2"
  } else {
    set result $dim1
  }

  return $result
}


#########################################
#########################################
### B. Parameter line --> Data fields ###
#########################################
#########################################


# Fill array data fields (like InitialCondition(TEMPERATURE) from
# parameter dataline
# NOTE: Var-name indexing by system-id is used mainly for solvers!
#
proc DataField::formDataFields {globArray data_line {problem_mask "" } } {
  upvar #0 $globArray theArray
  global Common Info Model

  set panel $theArray(parameterType)

  #-1. Split data row into a list using separators like '|'
  set disp_sep $Info(dispListSeparatorIn)
  set data [split $data_line $disp_sep]

  #-2. Pick the size indicator (like (2 2) ) from the data
  # and collect cleaned data into values-list.
  set var_val_list ""
  set actives_list ""

  foreach val $data {
    
    if { $Info(inactiveMark) == [string index $val 0] } {
      lappend actives_list 0
    } else {
      lappend actives_list 1
    }

    # Remove possible inactive mark from the beginning
    set clean_val [string trimleft $val $Info(inactiveMark)]

    # Pick the possible file size indicator
    set field_size [DataField::extractFieldSize cleaned_result $clean_val]

    # Append to the "pure" list of var-name+value
    lappend var_val_list $cleaned_result
  }


  #-3. Pick off object identifier, possible variable description and 
  # put the data-value into the field variable.
  # Check also if value is for a procedure name.
  set var_index 0

  # Reset fields
  #
  #Panel::resetFields $globArray $theArray(allFields)

  set len_before [llength $var_val_list]
  set len_after -1

  while { $var_val_list != ""  && $len_after < $len_before } {

    set len_before [llength $var_val_list]

    #--Loop array fields
    foreach var $theArray(allFields) {
    
      if { ![Panel::isOutputField $globArray $var]} {
        incr var_index
        continue
      }
    
      # Drop inactive field
      #
      if { $problem_mask != "" &&
           !$Info(keepInactiveFields) &&
           ![DataField::fieldTargetMaskMatches $globArray $var $problem_mask]
         } {
        incr var_index
        continue
      }

      if { [catch {

          if { [Panel::isGroupField $globArray $var] } {
            set pos_index [DataField::fillGroupFromDispList $var_val_list $actives_list $globArray $var $var_index]
          } else {
            set pos_index [DataField::fillFieldFromDispList $var_val_list $actives_list $globArray $var $var_index]
          }
    
        } err_msg ] && $Info(showDataFormatErrorMsg) } {
    
        set msg "Error in forming data for ($globArray) $var because: \n\n"
        append msg $err_msg
        print "$msg"

        set Info(showDataFormatErrorMsg) 1
      }

      incr var_index
    
      #--Remove found field from search list

      #-If field not found
      if { $pos_index == -1 } {
        continue
      }
      
      #-Remove matched items from the lists
      set var_val_list [lreplace $var_val_list $pos_index $pos_index]
      set actives_list [lreplace $actives_list $pos_index $pos_index]

      set theArray($var,firstInstance) 0
      set theArray($var,prev) $theArray($var)
      set theArray($var,old) $theArray($var)
    
      #-If field needs special formatting
      if { [info exists Common($panel,$var,data2WidgetProc)] } {
        Util::execArrayProc $globArray $var,data2WidgetProc
      }
    
    } ;# Foreach var

    set len_after [llength $var_val_list]

  } ;# While var_val_list


  # If we should inform about removed fields!
  # =========================================
  if { $Info(informAboutRemovedFields) && $var_val_list != "" } {

    foreach pair $var_val_list act $actives_list {
      
      if { !$act } {
        continue
      }

      set var [lindex [split $pair "="] 0]
      
      if { ![Panel::isOutputField $globArray $var] } {
        continue
      }

      # Some specially named fields do not tell anything
      # to the user!
      #
      set var2 [string toupper $var]
      if { ![string equal $var $var2] } {
        continue
      }

      set arr [Panel::panelNameGuiToMenu $globArray]
      set var [DataField::fieldNameGuiToSif $var]
      
      lappend Model(removedFieldsList) "$arr:    $var\n"
    }
  }

}


#-Find an element from a list based on TEMP== type key.
# Keys are suppoused to be in the beginning of each list item.
# Finally trim away the key-part from identified element.
# Example:
# TEMPERATURE== means a normal value     (see Info(valueFlag)
# TEMPERATURE=: means a filename         (see Info(fileFlag)
# TEMPERATURE=. means a procedure name   (see Info(procFlag)
#
proc DataField::fillFieldFromDispList { var_val_list actives_list globArray var var_index } {
  global Common Info
  upvar #0 $globArray theArray

  set panel $theArray(parameterType)

  #-Extra fill characters around separators (normally spaces )
  set extra1 $Info(dispListSeparatorExtra1)
  set extra2 $Info(dispListSeparatorExtra2)

  set procedurable [DataField::getFieldProperty $globArray $var Procedurable]
  set tableable    [DataField::getFieldProperty $globArray $var Tableable]

  set variable_id $var
  # Add "=" at the end to prevent partial matching like:
  # HEAT_FLUX_BC matching too early with HEAT_FLUX
  append variable_id "="

  set index [lsearch -regexp $var_val_list ^$extra1$variable_id]

  #---Field not found
  if {$index == -1} {
    return -1
  }

  #---Insert values into the global field variable
  #
  set data [lindex $var_val_list $index]
  set active  [lindex $actives_list $index]

  # Handle Include field separately because it is a panel level field 
  # (no mask!) and has of the Use-checkbox etc.
  #
  if { $var == "INCLUDE" } {
    set theArray(INCLUDE,Use,act) 1

    if { $active } {
      set theArray(INCLUDE,act) 1
      set theArray(INCLUDE,Use) 1
    } else {
      set theArray(INCLUDE,act) 0
      set theArray(INCLUDE,Use) 0
    }
  }

  # If activity is set for the field
  #
  if { [info exists theArray($var,act) ] } {
    set active $theArray($var,act)
  }

  #-Variable-id (like TEMP) is trimmed away.
  regsub $extra1$variable_id $data {} rest

  #-Trim fill-character(s) extra1 before value-type indicator
  set rest [string trimleft $rest $extra1]

  #-Add "=" character back at front of the data, so that the following flags
  # can match
  set rest "=$rest"
  set flag ""
  set is_file 0
  set is_proc 0
  set is_matc 0
  set is_table 0

  #-Is this a numeric field
  set ftype [DataField::getFieldProperty $globArray $var FieldValueType]
  set is_numeric 0

  if { $ftype == "Real" || $ftype == "Integer" } {
    set is_numeric 1
  }

  #-Is value a file, procedure or matc-script(==$)?
  if { $Info(fileFlag) == [string range $rest 0 1] } {
    set is_file 1
    set flag $Info(fileFlag)

  } elseif { $Info(procFlag) == [string range $rest 0 1] } {
    set is_proc 1
    set flag $Info(procFlag)

  } elseif { $Info(matcFlag) == [string range $rest 0 2] } {
    set is_matc 1
    set flag ""
  }

  set rest [string trimleft  $rest $Info(valueFlag)$flag$extra2]
  set rest [string trimright $rest $extra2]

  set vn $var

  #-Pick possible variable and size description
  set description [DataField::extractVariableAndSizeDescription cleaned $rest]

  # Check if field is a table (or array)
  if { $tableable && !$is_proc && !$is_matc &&
       ( $description != "" ||
         ($rest != "" && $is_numeric && ![string is double -strict $rest])
       ) } {
    set is_table 1
  }

  #-Now we can set data type
  if { $is_proc } {
    set dtype procData

  } elseif { $is_table } {
    set dtype tableData

  } elseif { $is_matc } {
    set dtype valueData

  } else {
    set dtype valueData
  }

  set check 0
  setVariableAndSizeValues $globArray $vn $dtype $description $check

  set rest $cleaned

  # Format numeric data 
  # ===================
  if { $rest != "" && $is_numeric && !($is_file || $is_proc || $is_matc) } {

    set frmt [DataField::getFieldProperty $globArray $var FieldFormat]

    if {$frmt != ""} {
      if { [catch  {set rest [DataField::formatNumbers $rest $frmt] } err_msg ] } {
        set msg "Invalid data: $rest \nfor ($globArray) in field: $var\n\n"
        append msg $err_msg

        set Info(messageIcon) error
        Message::dispOkMessage $msg
      }
    }
  }

  set rest [string trim $rest]

  set theArray($vn,$dtype) $rest

  #-First occurance of a field defines the active type and field value
  #
  if { $theArray($vn,firstInstance) } {
    set theArray($vn) $rest
    set theArray($vn,act) $active
    set theArray($vn,proc) $is_proc
    set theArray($vn,table) $is_table

    if { [info exist theArray($vn,$dtype,dataSize)] } {
      set theArray($vn,dataSize) $theArray($vn,$dtype,dataSize)
    }

    if { [info exist theArray($vn,$dtype,variables)] } {
      set theArray($vn,variables) $theArray($vn,$dtype,variables)
    }

    if { [info exist theArray($vn,$dtype,nofEntries)] } {
      set theArray($vn,nofEntries) $theArray($vn,$dtype,nofEntries)
    }
  }

  return $index
}


proc DataField::fillGroupFromDispList { var_val_list actives_list globArray var var_index } {
  global Common Info $globArray
  upvar #0 $globArray theArray

  set panel $theArray(parameterType)

  set extra1 $Info(dispListSeparatorExtra1)
  set extra2 $Info(dispListSeparatorExtra2)
  set variable_id $var
  set index [lsearch -regexp $var_val_list ^$extra1$variable_id]

  #---Field not found
  if {$index == -1} {
    return -1
  }

  #---Insert values into the global field variable
  set data [lindex $var_val_list $index]
  set active  [lindex $actives_list $index]

  # Handle Include field separately because it is a panel level field 
  # (no mask!) and has of the Use-checkbox etc.
  #
  if { $var == "INCLUDE" } {
    set theArray(INCLUDE,Use,act) 1

    if { $active } {
      set theArray(INCLUDE,act) 1
      set theArray(INCLUDE,Use) 1
    } else {
      set theArray(INCLUDE,act) 0
      set theArray(INCLUDE,Use) 0
    }
  }

  # If activity is set for the field
  #
  if { [info exists theArray($var,act) ] } {
    set active $theArray($var,act)
  }

  #-Variable-id (like TEMP=) is trimmed away.
  regsub $extra1$variable_id $data {} rest

  #-Trim extra1 before value-type indicator
  set rest [string trimleft $rest $extra1]

  set flag ""
  set is_file 0
  set is_proc 0
  set is_matc 0
  set is_table 0

  #-Is value a file or procedure name?
  if { $Info(fileFlag) == [string range $rest 0 1] } {
    set is_file 1
    set flag $Info(fileFlag)

  } elseif { $Info(procFlag) == [string range $rest 0 1] } {
    set is_proc 1
    set flag $Info(procFlag)

  } elseif { $Info(matcFlag) == [string range $rest 0 2] } {
    set is_matc 1
    set flag ""
  }

  set rest [string trimleft  $rest $Info(valueFlag)$flag$extra2]
  set rest [string trimright $rest $extra2]


  #---If value is a file or procedure name, put it to the first field of the group
  if { $is_file || $is_proc } {

    set vn [lindex $Common($panel,$var,group) 0]

    set rest [string trim $rest]

    if { $theArray($vn,firstInstance) } {
      set theArray($vn) $rest
      set theArray($vn,proc) $is_proc
    }

    #-Store possible procedure/table values separately
    if { $is_proc } {
      set theArray($vn,procData) $rest

      if { [info exist theArray($vn,dataSize)] } {
        set theArray($vn,procData,dataSize) $theArray($vn,dataSize)
      }

      if { [info exist theArray($vn,variables)] } {
        set theArray($vn,procData,variables) $theArray($vn,variables)
      }
    }

    return $index
  }

  #---Otherwise split data-line into group's entry fields

  #-Take last separator character away from the end of the data
  set len [string last $Info(groupDataSeparator) $rest]
  incr len -1
  set rest [string range $rest 0 $len]


  #-If we have a radio-button group, we update only group's parent's value
  if { [Panel::isRadioButtonField $globArray $var] } {

    set vn $var

    set theArray($vn) $rest
    set theArray($vn,act) $active

    Util::execArrayProc $globArray "$var,groupEditProc" 
    return $index
  }   

  #-Otherwise loop data-list and group-variable list in parallel
  set dlist [split $rest $Info(groupDataSeparator)]
  set glist $Common($panel,$var,group)

  #-Append empty-lists at the end of dlist-variable
  # so that all data-fields are set a value or blanked
  #
  set app_cnt [expr [llength $glist] - [llength $dlist]]
  for {set i 0} {$i < $app_cnt} {incr i} {
    lappend dlist ""
  }

  foreach fv $dlist vn $glist {

    if {$fv != ""} {
      set frmt [DataField::getFieldProperty $globArray $vn FieldFormat]
      if {$frmt != ""} {
        set fv [format $frmt $fv]
      }
    }

    #-First occurance of a field defines the active type and field value
    #
    if { $theArray($vn,firstInstance) } {
      set theArray($vn) $fv
      set theArray($vn,act) $active
      set theArray($vn,proc) $is_proc
      set theArray($vn,table) $is_table
    }

    #-Set data type
    if { $is_proc } {
      set dtype procData
    } elseif { $is_table } {
      set dtype tableData
    } else {
      set dtype valueData
    }

    set theArray($vn,$dtype) $fv

    # Store other possible data components
    if { $is_proc || $is_table } {

      if { [info exist theArray($vn,dataSize)] } {
        set theArray($vn,$dtype,dataSize) $theArray($vn,dataSize)
      }

      if { [info exist theArray($vn,variables)] } {
        set theArray($vn,$dtype,variables) $theArray($vn,variables)
      }
    }

  }

  return $index
}


# Pick variable description from parameter field
#
proc DataField::pickVariableAndSizeDescription {field_data} {
  global Info

  set pos1 [string first $Info(leftFieldSizeSep) $field_data ]
  set pos2 [string first $Info(rightFieldSizeSep) $field_data ]

  set dstart [incr pos1]
  set dend [incr pos2 -1]


  set var_description [string range $field_data $dstart $dend]

  return $var_description
}


# Extract variable description like (Velocity 1) from parameter line
# Return description part and put cleaned (= size removed) data
# into result
#
proc DataField::extractVariableAndSizeDescription {result field_data} {
  global Info
  upvar $result cleaned_result

  set field_data [string trim $field_data " \t"]

  set cleaned_result $field_data
   
  if { $field_data == "" } {
    set cleaned_result ""
    return ""
  }

  # Size descriptor must be in the beginning!
  #
  set pos1 [string first $Info(leftFieldSizeSep) $field_data ]

  if { $pos1 != 0 } {
    return ""
  }

  set pos2 [string first $Info(rightFieldSizeSep) $field_data ]

  set end1 $pos1; incr end1 -1
  set start2 $pos2; incr start2 

  incr pos1; incr pos2 -1

  set var_desc [string trim [string range $field_data $pos1 $pos2]]

  set tmp [string range $field_data 0 $end1]
  append  tmp [string range $field_data $start2 end]
  set cleaned_result $tmp

  return $var_desc
}


# Form display format description for the field
# like (Velocity 1)
#
proc DataField::formVariableAndSizeDescription {globArray fld} {
  global Info

  upvar #0 $globArray theArray

  if { ![info exists theArray($fld,dataSize)] } {
    return ""
  }

  if { ![info exists theArray($fld,variables)] } {
    set theArray($fld,variables) ""
  }

  set description ""
 
  set sz1 [lindex $theArray($fld,dataSize) 0]
  set sz2 [lindex $theArray($fld,dataSize) 1]

  if { ($sz1 == "" || $sz1 == 1) && ($sz2 == "" || $sz2 == 1) } {
    set add_size 0
  } else {
    set add_size 1
  }


  if { $theArray($fld,variables) != "" && $add_size } {
    set description_sep $Info(dataListSeparator)
  } else {
    set description_sep ""
  }

  # If variable or size different from 1 1, form description
  #
  if { $theArray($fld,variables) != ""  || $add_size } {
    append description "("
    append description $theArray($fld,variables)
    append description $description_sep
    if {$add_size} {
      append description $theArray($fld,dataSize)
    }
    append description ")"
  }

  return $description
}


# Set field variable and data size info from description
#
proc DataField::setVariableAndSizeValues {globArray var data_type description check } {
  global Info $globArray
  upvar #0 $globArray theArray

  if { $data_type != "" } {
    set vn "$var,$data_type"
  } else {
    set vn $var
  }

#  if { ![info exists theArray($vn)] ||
#       $theArray($vn)  == ""
#     } {
#    return
#  }

  #---Default values
  if { ![info exists theArray($vn)] || 
       $theArray($vn) == "" ||
       ![info exist theArray($vn,dataSize)]
     } {
    set theArray($vn,dataSize) [list 1 1]
    set theArray($vn,entrySize) 1
    set theArray($vn,variables) ""
    set theArray($vn,variablesSize) 0
  }

  #---If no description given, use default values
  if { $description == "" } {
    return ""
  }

  #---Parse description
  set desc_list [split $description $Info(dataListSeparator)]

  set variables ""
  set dsize ""

  set nof_items [llength $desc_list]
  set nof_variables 0

  set counter 0
  # Pick all variables and the possible size info
  foreach item $desc_list {

    incr counter

    set item [string trim $item]

    # Alphabetic <--> VARIABLE NAME info
    if { [catch {expr [join "$item 1" *]} msg] } {

      if {$check} {
        set msg ""
        if {![DataField::isValidVariable $globArray $item]} {
          append msg "Incorrect variable name ($item) for field\n"
          append msg [string toupper $var]
          return $msg
        }
      }

      if { $nof_variables > 0 } {
        append variables "$Info(dataListSeparator) "
      }

      append variables $item

      incr nof_variables          

    # Numeric <--> SIZE info
    } else {

      if {$check} {
        set msg ""
        if { $counter != $nof_items } {
          append msg "Incorrect data size description for field\n"
          append msg [string toupper $var]
          append msg "\n\nEntry size must the last data item in description!"
          return $msg
        }
      }

      # Pick size description
      set dsize $item
    }
  }

  set theArray($vn,variables) $variables

  set len [llength $dsize]

  # If we dont have variables
  if { $len > 2 } {
    append msg "Incorrect data size description for field\n"
    append msg [string toupper $var]
    append msg "\n\nEntry size is at most 2 numbers!"
    return $msg
  }

  # Data size at least (1 1)
  lappend dsize 1
  lappend dsize 1
  set dsize [lrange $dsize 0 1]

  set theArray($vn,dataSize) $dsize

  set theArray($vn,variablesSize) [llength $theArray($vn,variables)]

  set vz $theArray($vn,variablesSize)

  set s1 [lindex $dsize 0]
  set s2 [lindex $dsize 1]
  set theArray($vn,entrySize) [expr $vz + ($s1 * $s2)]

  return ""
}


proc DataField::isValidVariable {globArray variable} {
  upvar #0 $globArray theArray

  set external_names [lindex $theArray(currentVariables) 1]

  set index [lsearch -exact $external_names $variable]

  if { $index == -1 } {
    return 0
  }

  return 1
}


# Extract size description from parameter line
# Return field size and put cleaned (= size removed) data
# into result
#
proc DataField::extractFieldSize_old {result field_data} {
  global Info
  upvar $result cleaned_result

  set pos1 [string first $Info(leftFieldSizeSep) $field_data ]
  set pos2 [string first $Info(rightFieldSizeSep) $field_data ]

  set end1 $pos1; incr end1 -1
  set start2 $pos2; incr start2 

  incr pos1; incr pos2 -1

  set field_size [string range $field_data $pos1 $pos2]
  set tmp [string range $field_data 0 $end1]
  append  tmp [string range $field_data $start2 end]
  set cleaned_result $tmp

  return $field_size
}


# Extract size description from parameter line
# Return field size and put cleaned (= size removed) data
# into result
# NOTE: This version does not remove a variable description
# by accident MVe 11.02.99 !!!
#
proc DataField::extractFieldSize {result field_data} {
  global Info
  upvar $result cleaned_result

  set cleaned_result $field_data

  set pos1 [string first $Info(leftFieldSizeSep) $field_data ]
  set pos2 [string first $Info(rightFieldSizeSep) $field_data ]

  if { $pos1 != 0 } {
    return ""
  }

  set end1 $pos1; incr end1 -1
  set start2 $pos2; incr start2 

  incr pos1; incr pos2 -1

  set field_size [string trim [string range $field_data $pos1 $pos2] ]

  # Check that field-size is purely numerical, so it cannot
  # be a variable description
  set msg ""
  foreach item $field_size {

    if { [catch {expr $item} msg] } {
      return ""
    }
  }

  # Ok, we really have a size description
  #
  set tmp [string range $field_data 0 $end1]
  append  tmp [string range $field_data $start2 end]

  set cleaned_result $tmp

  return $field_size
}


# end ecif_tk_procsDataField.tcl
# ********************
