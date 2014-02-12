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
#Module:    ecif_tk_procsUserDefined.tcl
#Language:  Tcl
#Date:      15.01.01
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  User defined equation handling
#
#************************************************************************

global UserDefined


#================================================================#
#                     Definition tables                          #
#================================================================#

# Accepted generic keywords
# =========================
set UserDefined(genericKeywords) {
    "Include" "Equation" "Panel" "End File"
}

# Accepted equation section keywords
# ==================================
set UserDefined(equationKeywords) {
    "Variable" "Variable Components" "Variable DOFs"
    "Exported Variable" "Exported Variable DOFs"
    "Solver Procedure"
    "Process Table Name" "Status Line Name"
    "Panel"
    "End Equation"
}

# Accepted panel section keywords
# ===============================
set UserDefined(panelKeywords) {
    "Include"
    "Field"
    "End Field"
    "End Panel"
}

# Accepted field property keywords
# ================================
set UserDefined(panelFieldKeywords) {
    "Type" "Size" "Display"
    "Widget" "Widget Width" "Widget Height"
    "Sif Name"
    "Label" "Label Width" "Unit Label"
    "Initial Value" "Always Initialize" 
    "Drop Initial Valued" "Always Output"
    "Activity Parent" "Activity Parents"
    "Activity Value" "Activity Values"
    "Can Be Procedure" "Can Be Table" "Can Has Variable"
    "PanelPage" "Panel Page" 
    "PanelFrame" "Panel Frame" 
    "SubPanel" "Sub Panel" 
    "Screen PadX" "Screen PadY"
    "Index Type"
    "Values"
    "Add Absolute Box"
    "Field Target"
}

# Accepted equation variable property keywords
# ============================================
# NOTE: These properties can be given to the variable which
# is solved in the equation. They are given in the main
# equation block
#
set UserDefined(equationVarFieldKeywords) {
    "Display"
    "Widget" "Widget Width" "Widget Height"
    "Label" "Label Width" "Label By Coordinate" "Unit Label"
    "Initial Value" "Always Initialize"
    "Activity Parent" "Activity Parents"
    "Activity Value" "Activity Values"
    "Can Be Procedure" "Can Be Table" "Can Has Variable" 
    "PanelPage" "Panel Page" 
    "PanelFrame" "Panel Frame" 
    "SubPanel" "Sub Panel" 
    "Screen PadX" "Screen PadY"
    "Index Type"
}
  
# Add field keywords to keywords accepted in
# Equation and in Panel blocks

#-Panel block
foreach fld_keyword $UserDefined(panelFieldKeywords) {
  lappend UserDefined(panelKeywords) $fld_keyword
}

#-Equation block
foreach fld_keyword $UserDefined(equationVarFieldKeywords) {
  lappend UserDefined(equationKeywords) $fld_keyword
}


# Accepted panel names (in compressed format!) and
# their gui equivalent formats
#
# 1. group: Menu based names
# 2. group: Sif based names
# 3. group: Alternative (all imaginable) versions
#
# 1. Panel file name (in compressed format)
# 2. Gui array name
#
set UserDefined(panelNameTable) {

    {Equations Equation}
    {BodyForces BodyForce}
    {BodyParameters BodyParameter}
    {BoundaryConditions BoundaryCondition}
    {BoundaryParameters BoundaryParameter}
    {InitialConditions InitialCondition}
    {Materials Material}
    {ModelParameters ModelParameter}
    {PhysicalConstants Constant}
    {SimulationParameters SimulationParameter}
    {SolverControl SolverControl}
    {SolverSettings Solver}

    {Constant Constant}
    {Equation Equation}
    {BodyForce BodyForce}
    {BoundaryCondition BoundaryCondition}
    {InitialCondition InitialCondition}
    {Material Material}
    {Simulation SimulationParameter}
    {Solver Solver}

    {EquationParameter Equation}
    {EquationParameters Equation}
    {EquationSettings Equation}
    {BodyForceParameter BodyForce}
    {BodyForceParameters BodyForce}
    {BodyForceSettings BodyForce}
    {BodyParameter BodyParameter}
    {BodyParameterSettings BodyParameter}
    {BoundaryConditionParameter BoundaryCondition}
    {BoundaryConditionParameters BoundaryCondition}
    {BoundaryConditionSettings BoundaryCondition}
    {BoundaryParameter BoundaryParameter}
    {BoundaryParameterSettings BoundaryParameter}
    {InitialConditionParameter InitialCondition}
    {InitialConditionParameters InitialCondition}
    {InitialConditionSettings InitialCondition}
    {MaterialParameter Material}
    {MaterialParameters Material}
    {MaterialParameterSettings Material}
    {Model ModelParameter}
    {ModelParameter ModelParameter}
    {ModelParameterSettings ModelParameter}
    {PhysicalConstant Constant}
    {PhysicalConstantSettings Constant}
    {SimulationParameter SimulationParameter}
    {SimulationParameterSettings SimulationParameter}
    {Solvers Solver}
    {SolverParameter Solver}
    {SolverParameters Solver}
    {SolverParameterSettings Solver}
    {SolverControlParameter SolverControl}
    {SolverControlParameters SolverControl}
    {SolverControlSettings SolverControl}

}


# Accepted property names in compressed format and
# their gui equivalent formats
#
# NOTE: Accepted keywords must be first defined by tre two variables above:
#
# UserDefined(panelFieldKeywords)
# UserDefined(equationVarFieldKeywords)
#
# List list just defines their gui-equivalents!!!
#
# 1. Definition file name (in compressed format)
# 2. Gui property name
#
set UserDefined(propertyKeywordTable) {
    {Variable             Variable}
    {VariableDofs         VariableDofs}
    {VariableComponents   VariableComponents}
    {ExportedVariable     ExportedVariable}
    {ExportedVariableDofs ExportedVariableDofs}
    {VariableDofs         VariableDofs}
    {SolverProcedure      Procedure}
    {StatusLineName       StatusLineName}
    {ProcessTableName     ProcessTableName}
    {Widget               WidgetType}
    {WidgetWidth          WidgetWidth}
    {WidgetHeight         WidgetHeight}
    {Type                 FieldValueType}
    {Size                 FieldDataSize}
    {SifName              SifName}
    {Label                Label}
    {LabelWidth           LabelWidth}
    {LabelByCoordinate    LabelByCoordinate}
    {UnitLabel            UnitString}
    {PanelPage            PanelPage}
    {PanelFrame           SubPanel}
    {SubPanel             SubPanel}
    {ScreenPadX           ScreenPadX}
    {ScreenPadY           ScreenPadY}
    {CanBeProcedure       Procedurable}
    {CanBeTable           Tableable}
    {CanHasVariable       Variabled}
    {IndexType            IndexType}
    {InitialValue         InitialValue}
    {AlwaysInitialize     AlwaysInitialize}
    {DropInitialValued    DropInitialValued}
    {AlwaysOutput         AlwaysOutput}
    {ActivityParent       ActivityParents}
    {ActivityParents      ActivityParents}
    {ActivityValue        ActivityValues}
    {ActivityValues       ActivityValues}
    {Values               Limits}
    {Display              Display}
    {AddAbsoluteBox       HasAbsCheckBox}
    {FieldTarget          FieldTarget}
}

#================================================================#
#================================================================#


proc UserDefined::openDefinitionFile {filename {dir ""} } {
  global Info

  # Try first to open with pure filename
  #
  set path1 $filename
  set path2 $filename

  set cmsg ""
  if { [ catch {set in [open $path1 "r"]} cmsg ] } {
    if { $dir != "" } {
      set path2 [file join $dir $filename]
    }
  }      

  # Then add possible dir and try again
  #
  if { [ catch {set in [open $path2 "r"]} cmsg ] } {

    Message::showMessage ""
    Message::showMessage "WARNING*** $cmsg" $Info(wrnMsgColor)
    Message::showMessage ""
    
    #set msg "Cannot open user's equation definition file: $filename"
    #set Info(messageIcon) error
    #Message::dispOkCancelMessage $msg

    return -1
  }

  return $in
}



#========================#
#   readDefinitionFile   #
#========================#
#
# Main call for user defined equation
#
proc UserDefined::readDefinitionFile { filename } {
  global Equation Info Model
  
  # =========
  # Open file
  # =========
  
  set in [UserDefined::openDefinitionFile $filename]

  if { $in < 0 } {
    return 0
  }

  set lastEqId 0    ;# Last equation id to be used
  set newLastEqId 0 ;# Updated last equation id when read equations has been added
  set eqDefIds ""   ;# List of equation ids
  
  # eqDefs          ;# Equation definition array, stores generic data for equation. Keys: eq-id
  # panelDefs       ;# Panel definition array, stores field data fore panels. Keys: eq-id, panel name, field name
  # panelNames      ;# Panel names array, needed for looping panelDefs. Keys: eq-id, panel name, field name
  # fieldNames      ;# Field names array, needed for looping panelDefs, Keys: eq-id, panel name, field name

  Message::showMessage "Reading definition file:  $filename" $Info(remMsgColor)

  if { ![UserDefined::readDefinitionData \
            $filename $in $lastEqId  \
            newLastEqId eqDefIds \
            eqDefs panelDefs     \
            panelNames fieldNames]
     } {

    Message::showMessage ""
    Message::showMessage "ERROR IN DEFINITION FILE. NO DEFINITIONS APPLIED!" $Info(errMsgColor)
    Message::showMessage ""
    set Model(newDefinitions) ""
    close $in
    return 0
  }

  close $in

  if { ![UserDefined::addDefinitionData $eqDefIds eqDefs panelDefs panelNames fieldNames] } {
    Message::showMessage ""
    Message::showMessage "ERROR IN DEFINITION FILE. NO DEFINITIONS APPLIED!" $Info(errMsgColor)
    Message::showMessage ""
    set Model(newDefinitions) ""
    return 0
  }

  Message::showMessage "***All definitions read from file $filename"
  Message::showMessage ""

  # If new definitions were succesfully read
  #
  if { $Model(newDefinitions) != "" } {

    foreach line $Model(newDefinitions) {
      lappend Model(allDefinitions) $line
    }

    set Model(hasUserDefinitions) 1
  }

  # Store filename in info-list
  Util::addToFileList Info(definitionFiles) $filename

  return 1
}


#========================#
#   readDefinitionData   #
#========================#
#
# Read definition file contents
#
# NOTE: proc UserDefined::addToDefinitionLog adds definition
# lines to the global arraies:
# Model(allDefinitions)
# Model(newDefinitions)
#
proc UserDefined::readDefinitionData { filename channel_id lastEqId
                                       newLastEqId eqDefinitionIds
                                       eqDefinitions panelDefinitions
                                       panelNameLists fieldNameLists
                                       {currentEquation "" }
                                       {currentPanel "" }
                                       {currentField "" }
                                       {currentPanelVar "" }
                                       {currentFieldVar "" }
                                      } {
  global UserDefined
  global Equation Info Model

  upvar $newLastEqId newLastEq
  upvar $eqDefinitionIds eqDefIds
  upvar $eqDefinitions eqDefs
  upvar $panelDefinitions panelDefs
  upvar $panelNameLists panelNames
  upvar $fieldNameLists fieldNames

  set eq $lastEqId
  set newLastId $lastEqId

  set defIndentLevel 0


  # This sort order should quarantee that "longer" keywords are checked
  # first ("Variable DOFs" before "Variable" etc.)
  #
  set generic_keywords $UserDefined(genericKeywords)
  set generic_keywords [lsort -ascii -decreasing $generic_keywords]

  set equation_keywords $UserDefined(equationKeywords)
  set equation_keywords [lsort -ascii -decreasing $equation_keywords]

  set panel_keywords $UserDefined(panelKeywords)
  set panel_keywords [lsort -ascii -decreasing $panel_keywords]

  set filename [Util::makeTclFilename [Util::stringRemoveQuotes $filename]]
  
  set UserDefined(currentFileName) $filename

  # =========
  # Read file
  # =========
  
  set current_equation $currentEquation
  set current_panel $currentPanel
  set current_field $currentField
  set current_panel_var $currentPanelVar
  set current_field_var $currentFieldVar

  set err_file "ERROR! When reading definition file   \"$filename\"\n\n"
  set wrn_file "Warning! When reading definition file   \"$filename\"\n\n"

  while {1} {

    set line [Util::readInputFileLine $channel_id]

    set line [Util::stringTrim $line]

    set err_line "\nin the input line:   $line\n\n"

    # No more data (eof)
    # ==================
    if { $line == "" } {
      UserDefined::addToDefinitionLog "End File !********* ($filename)" 0
      return 1
    }
    
    # Read next (accepted) keyword
    # ============================
    set keyword ""
    set value ""

    #-Generic keyword
    # ---------------
    if { $current_equation == "" && $current_panel == "" } {
      set err_block "for general definitions\n"  
      set rc [UserDefined::readAcceptedKeyword $generic_keywords $line keyword value]

    #-Equation keyword
    # ----------------
    } elseif { $current_panel == "" } {
      set err_block "for Equation \"$current_equation\"\n"
      set rc [UserDefined::readAcceptedKeyword $equation_keywords $line keyword value]

    #-Panel keyword
    # -------------
    } elseif { $current_panel != "" } {
      set err_block "for Panel \"$current_panel\"\n"
      set rc [UserDefined::readAcceptedKeyword $panel_keywords $line keyword value]
    }

    #-Unknown keyword!
    # ===============
    if {!$rc} {

      if { $keyword != "" } {
        set err "Unknown value   $value  for $keyword\n"
      } else {
        set err "Unknown keyword\n"
      }
      
      set Info(messageIcon) error
      Message::dispOkMessage [list $err_file $err $err_line $err_block]
      return 0
    }

    # Check section changes
    # =====================

    #-New equation started
    # --------------------
    if { 0 == [string compare -nocase $keyword "Equation"] } {
      
      set current_equation [Util::stringRemoveQuotes $value]
      set current_panel ""
      set current_field ""
      
      if { $current_equation == "" } {
        set err "\nEquation name missing after Equation keyword!\n"
  
        set Info(messageIcon) error
        Message::dispOkMessage [list $err_file $err]
        return 0
      }
      
      incr eq
      set newLastEq $eq

      lappend eqDefIds $eq      
      set eqDefs($eq,NAME) $current_equation
      set eqDefs($eq,ExportedVariable) ""
      set eqDefs($eq,ExportedVariableDofs) ""
      set panelNames($eq) ""
      
      set defIndentLevel 0
      UserDefined::addToDefinitionLog " "
      UserDefined::addToDefinitionLog $line $defIndentLevel
      incr defIndentLevel 1
      continue
    
    #-Equation ended
    # --------------
    } elseif { 0 == [string compare -nocase $keyword "End Equation"] } {
      
      set current_equation ""
      set current_panel ""
      set current_field ""

      set defIndentLevel 0
      UserDefined::addToDefinitionLog $line $defIndentLevel
      continue

    #-New panel started
    # -----------------
    } elseif { 0 == [string compare -nocase $keyword "Panel"] } {

      set current_panel [Util::stringRemoveQuotes $value]

      if { $current_panel == "" } {
        set err "\nPanel name missing after Panel keyword!\n"
  
        set Info(messageIcon) error
        Message::dispOkMessage [list $err_file $err]
        return 0
      }

      set name_tbl $UserDefined(panelNameTable)

      if { ![Util::convertPanelNameToGuiFormat $name_tbl $current_panel current_panel_var] } {
        set err "Unknown panel name:  $current_panel\n"
  
        set Info(messageIcon) error
        Message::dispOkMessage [list $err_file $err $err_line]
        return 0
      }

      set current_field ""

      # If this an "all-equation" panel (ie. a panel outside 
      # any equation definition)
      if { $current_equation == "" } {
  
        incr eq
        set newLastEq $eq

        lappend eqDefIds $eq
        set eqDefs($eq,NAME) ""
        set panelNames($eq) ""
        set defIndentLevel 0
      }

      set fieldNames($eq,$current_panel_var) ""

      #-Add to the panel name list for the current equation
      # if not yet there!
      if { -1 == [lsearch $panelNames($eq) $current_panel_var] } {
        lappend panelNames($eq) $current_panel_var
      }
      
      UserDefined::addToDefinitionLog " "
      UserDefined::addToDefinitionLog $line $defIndentLevel
      incr defIndentLevel 1
      continue

    #-Panel ended
    # -----------
    } elseif { 0 == [string compare -nocase $keyword "End Panel"] } {
      
      set current_panel ""
      set current_field ""

      incr defIndentLevel -1
      UserDefined::addToDefinitionLog $line $defIndentLevel
      continue
    
    #-Panel level Include file (make a recursive call!)
    # ------------------------
    } elseif { $current_panel != "" &&
               0 == [string compare -nocase $keyword "Include"]
             } {

      set fn $value
      set base_dir [file dirname $filename]
      set inc_in [UserDefined::openDefinitionFile $fn $base_dir]
      
      if { $inc_in > 0 } {
        set ceqn $current_equation
        set cpnl $current_panel
        set cfld $current_field
        set cpvar $current_panel_var
        set cfvar $current_field_var

        set rc [UserDefined::readDefinitionData $fn $inc_in $eq newLastEq eqDefIds eqDefs panelDefs panelNames fieldNames $ceqn $cpnl $cfld $cpvar $cfvar]
        set eq $newLastEq
      }

      continue

    #-New field started
    # -----------------
    } elseif { $current_panel != "" && 
               0 == [string compare -nocase $keyword "Field"]
             } {
      set value  [Util::stringRemoveQuotes $value]

      if { $value == "" } {
        set err "\nField name missing after Field keyword!\n"

        set Info(messageIcon) error
        Message::dispOkMessage [list $err_file $err]
        return 0
      }

      set current_field $value
      set current_field_var [DataField::fieldNameSifToGui $current_field]

      #-Add to the field name list for the current panel
      # if not yet there!
      if { -1 == [lsearch $fieldNames($eq,$current_panel_var) $current_field_var] } {

        lappend fieldNames($eq,$current_panel_var) $current_field_var

        #UserDefined::setDefaultFieldProperties panelDefs $eq $current_panel_var $current_field_var ""
      }

      UserDefined::addToDefinitionLog $line $defIndentLevel
      continue

    #-Field ended (optional)
    # ----------
    } elseif { 0 == [string compare -nocase $keyword "End Field"] } {

      set current_field ""
      continue

    #-File ended (optional)
    # ----------
    } elseif { 0 == [string compare -nocase $keyword "End File"] } {

      #MSG "User definition file  \"$filename\"  read succesfully!"

      set current_equation ""
      set current_panel ""
      set current_field ""
      continue
          
    #-Property line
    # -------------
    } else {
      incr defIndentLevel 1
      UserDefined::addToDefinitionLog $line $defIndentLevel
      incr defIndentLevel -1
    }


    # Handle keyword
    # ==============
    #-Keyword variable name in the transfer array
    if { ![UserDefined::convertPropertyNameToGuiFormat $keyword gui_keyword] } {
      
      set err "\nIllegal keyword: $keyword\n"

      set Info(messageIcon) error
      Message::dispOkMessage [list $err_file $err $err_line $err_block]
      return 0
    }

   
    #-Panel level keyword
    # -------------------
    #
    # An entry like: panelDefs(1,InitialCondition,SOUNDP) <-- {LabelWidth 15}
    #
    if { $current_panel != "" } {

      if { $current_field == "" } {
        set err "\nField name missing, cannot store any value for property:  $keyword\n"
  
        set Info(messageIcon) error
        Message::dispOkMessage [list $wrn_file $err $err_line $err_block]
        return 0
      
      } else {
        lappend panelDefs($eq,$current_panel_var,$current_field_var) [list $gui_keyword $value]
      }

    #-Equation level keyword
    # ----------------------
    #
    # An entry like eqDefs(1,VariableDofs) <-- 2 3;1
    #
    } elseif { $current_equation != "" } {
      
      # There may be multiple entries for these keywords and they
      # are stored in lists!
      # NOTE: Nofs Dofs is kept the same as for variables although there were
      # more or less dofs-definitions!
      #
      if { $gui_keyword == "ExportedVariable" ||
           $gui_keyword == "ExportedVariableDofs"
         } {
        
        # Add possibly missing default Dofs value (1)
        #
        if { $gui_keyword == "ExportedVariable"  &&
             ( [llength $eqDefs($eq,ExportedVariableDofs)] <
               [llength $eqDefs($eq,ExportedVariable)]
             )
           } {
            lappend eqDefs($eq,ExportedVariableDofs) 1
        
        # Remove possibly repetitious Dofs values
        #
        } elseif { $gui_keyword == "ExportedVariableDofs"  &&
                   ( [llength $eqDefs($eq,ExportedVariableDofs)] >=
                     [llength $eqDefs($eq,ExportedVariable)]
                   )
                 } {
          set eqDefs($eq,ExportedVariableDofs) \
              [lreplace $eqDefs($eq,ExportedVariableDofs) end end]
        }
        
        lappend eqDefs($eq,$gui_keyword) $value
      
      # Other 'normal' equation level keywords
      #
      } else {
        set eqDefs($eq,$gui_keyword) $value
      }

    #-Generic keyword (Inlcude etc.)
    # ---------------
    } else {

      # Include file (make a recursive call!)
      # ------------
      if { 0 == [string compare -nocase $keyword "Include"] } {
        set fn $value
        set rc [UserDefined::readDefinitionData $fn $eq newLastEq eqDefIds eqDefs panelDefs panelNames fieldNames]
        set eq $newLastEq
      }
    }

  } ;#while(1)


  return 1
}


#=======================#
#   addDefinitionData   #
#=======================#
#
# Add definition data read from the file
#
# Returns: 1 <--> ok, 0 <--> error
#
proc UserDefined::addDefinitionData { eqDefIds
                                      eqDefinitions
                                      panelDefinitions
                                      panelNameLists
                                      fieldNameLists } {
  global Equation Info
    
  upvar $eqDefinitions eqDefs
  upvar $panelDefinitions panelDefs
  upvar $panelNameLists panelNames
  upvar $fieldNameLists fieldNames

  # ===================
  # Loop  equations ids
  # ===================
  foreach eq $eqDefIds {

    if { ![info exists eqDefs($eq,NAME)] } {
      continue
    }

    set en [DataField::fieldNameSifToGui $eqDefs($eq,NAME)]
   
    # Equation level data
    # ===================
    if { $en != "" } {

      #-Add a new equation
      # ------------------
      if { -1 == [lsearch $Equation(initialFields) $en] } {
        Message::showMessage "  Adding equation: $en"  
        
        if { ![UserDefined::addEquation $eq eqDefs] } {
          return 0
        }
    
      #-Update an existing equation
      # ---------------------------
      } else {
        Message::showMessage "  Updating equation: $en"  
        if { ![UserDefined::updateEquation $eq eqDefs] } {
          return 0
        }
      }
    }
    
    # Panel level data
    # ================

    #-An equation specific field
    if { $en != "" } {
      set emsk $eqDefs($eq,MASK)
      set eidx $eqDefs($eq,INDEX)
      set ename $eqDefs($eq,NAME)
    
    #-An "all-equations" field
    } else {
      set emsk ""
      set eidx ""
      set ename ""
    }

    #-Loop all panels
    foreach panel $panelNames($eq) {

      if { $eidx == "" } {
        Message::showMessage "  Panel: $panel"
      } else {
        Message::showMessage "      Panel: $panel"
      }

      foreach field $fieldNames($eq,$panel) {

        #Message::showMessage "Adding field: $panel $field emask=$emsk"

        if { ![UserDefined::addPanelField $eq $panel $field panelDefs $emsk $eidx $ename] } {
          return 0
        }
      }
    }
  }

  return 1
}


#=================#
#   addEquation   #
#=================#
#
# Add new equation
#
# Returns: 1 --> ok, 0 --> error
#
proc UserDefined::addEquation { eq defArray } {
  global Common Equation EquationVariable Solver Info
  upvar $defArray def

  set is_new_eq 1

  #-Check equation data
  # ===================
  if { ![UserDefined::checkEquationData $is_new_eq $eq def] } {
    return 0
  }

  #-Equation name
  # =============
  set sif_name [Util::stringRemoveQuotes $def($eq,NAME)]
  set ename [DataField::fieldNameSifToGui $sif_name]

  if { [info exists def($eq,LABEL)] } {
    set lbl [Util::stringRemoveQuotes $def($eq,LABEL)]
  } else {
    set lbl $sif_name
  }

  #-Find equation index
  # ===================
  set eq_index -1
  set uidx 0

  foreach fld $Equation(initialFields) {

    #-Find last used equation index and user defined eq. index
    #
    set idx [DataField::getFieldProperty Equation $fld EquationIndex "" 0]
    
    #-If this is an equation system field (ie. it has the EquationIndex
    # property value)
    #
    if { $idx != "" } {
      if { $idx > $eq_index } {
        set eq_index $idx
      }
    }
    
    #-User defined equation counter
    #  
    if { 1 == [DataField::getFieldProperty Equation $fld IsUserDefinedField "" 0] } {
      incr uidx
    }
  }

  #-Indices, ids and the mask
  # =========================
  incr eq_index
  incr uidx

  # Build a mask like ¤UD001
  #
  set ustr [format "%03d" $uidx]
  set eq_mask "(¤UD$ustr)"

  set def($eq,INDEX) $eq_index
  set def($eq,MASK) $eq_mask

  #-Add Exported variables
  UserDefined::addExportedEquationVariables $eq def


  #-Set Equation array properties
  # =============================

  #-'Equation' field
  #
  DataField::setFieldProperty Equation $ename Label $lbl
  DataField::setFieldProperty Equation $ename Display 0
  DataField::setFieldProperty Equation $ename FieldValueType Logical
  DataField::setFieldProperty Equation $ename Limits {Set {True ""}}
  DataField::setFieldProperty Equation $ename InitialValue False
  DataField::setFieldProperty Equation $ename UnitString ""
  DataField::setFieldProperty Equation $ename PanelPage 1
  DataField::setFieldProperty Equation $ename SubPanel 1
  DataField::setFieldProperty Equation $ename Mask ""

  DataField::setFieldProperty Equation $ename EquationIndex $eq_index
  DataField::setFieldProperty Equation $ename EquationLabel $lbl
  DataField::setFieldProperty Equation $ename EquationMask $eq_mask
  DataField::setFieldProperty Equation $ename EquationField $ename
  DataField::setFieldProperty Equation $ename EquationVars "$ename\_vars"
  DataField::setFieldProperty Equation $ename IsMultiVar 1
  DataField::setFieldProperty Equation $ename SifName $sif_name

  DataField::setFieldProperty Equation $ename IsUserDefinedField 1
  DataField::setFieldProperty Equation $ename IncludeInEquationPanel 1
  DataField::setFieldProperty Equation $ename SystemInfoIndices 0

  #-'Equation'_vars field
  #
  DataField::setFieldProperty Equation "$ename\_vars" Display 0
  DataField::setFieldProperty Equation "$ename\_vars" FieldValueType String
  DataField::setFieldProperty Equation "$ename\_vars" ExcludeFromSif 1

  DataField::addFieldMaskItem Equation VARIABLE_LIST [list "|$eq_mask" "|$eq_mask"]


  #-Set data in other arries
  # ========================

  #-EquationVariable
  #
  lappend EquationVariable(initialFields) $ename
  DataField::setFieldProperty EquationVariable $ename Display 0
  DataField::setFieldProperty EquationVariable $ename FieldValueType String
  DataField::setFieldProperty EquationVariable $ename FieldDataSep ";"
  DataField::setFieldProperty EquationVariable $ename InitialValue ""

  #-Solver
  #
  set proc [join $def($eq,Procedure) "\;"]
  DataField::setInitialValue Solver PROCEDURE $proc $eq_index
  DataField::setInitialValue Solver SOLVING_ORDER $eq_index $eq_index

  #-Common lists
  #
  lappend Common(allFieldNames) $ename
  lappend Common(allFieldNames) "$ename\_vars"

  lappend Equation(indices) $eq_index
  lappend Equation(initialFields) $ename
  lappend Equation(initialFields) "$ename\_vars"


  #-Set some default values for further processing
  # ==============================================
 
  #-If equation variable defined, default VariableDofs = 1
  #
  if { [info exists def($eq,Variable)] &&
       ![info exist def($eq,VariableDofs)]
     } {
    set def($eq,VariableDofs) 1
  }

  #-Possible equation variable components
  #
  if { ![info exist def($eq,VariableComponents)] } {
    set def($eq,VariableComponents) ""
  }

  #-Main window equation status line name
  #
  if { ![info exist def($eq,StatusLineName)] } {
    set def($eq,StatusLineName) "UD$uidx"
  }

  #-Process table status line name
  #
  if { ![info exist def($eq,ProcessTableName)] } {
    set def($eq,ProcessTableName) $def($eq,StatusLineName)
  }

  #-Try to add equation properties
  #
  if { ![UserDefined::setEquationFieldProperties $eq def $ename] } {
    return 0
  }

  return 1
}



#====================#
#   updateEquation   #
#====================#
#
# Update an existing equation
#
# Returns: 1 --> ok, 0 --> error
#
proc UserDefined::updateEquation { eq defArray } {
  upvar $defArray def

  set is_new_eq 0

  #-Check equation data
  # ===================
  if { ![UserDefined::checkEquationData $is_new_eq $eq def] } {
    return 0
  }

  set ename [DataField::fieldNameSifToGui [Util::stringRemoveQuotes $def($eq,NAME)]]

  set eindex [DataField::getFieldProperty Equation $ename EquationIndex "" 0]

  if { $eindex == "" } {
    return 0
  }

  set emask [DataField::getFieldProperty Equation $ename EquationMask]
  set def($eq,MASK) $emask
  set def($eq,INDEX) $eindex
  
  #-Add Exported variables
  UserDefined::addExportedEquationVariables $eq def

  #-Solver data
  #
  if { [info exists def($eq,SolverProcedure)] } {
    set proc [join $def($eq,SolverProcedure) "\;"]
    DataField::setInitialValue Solver SolverProcedure $proc $eq_index
  }

  #-Try to add equation properties
  #
  if { ![UserDefined::setEquationFieldProperties $eq def $ename] } {
    return 0
  }

  return 1
}


#======================#
#   checkEquationData  #
#======================#
#
# Check new equation data
#
# Returns: 1 --> ok, 0 --> error
#
proc UserDefined::checkEquationData { is_new_eq eq defArray } {
  global Info UserDefined
  upvar $defArray def

  set m0 "ERROR! When reading definition file   \"$UserDefined(currentFileName)\"\n\n"

  if { $is_new_eq && ![info exist def($eq,Variable)] } {

    set m1 "Equation variable name missing!\n\n"
    set m2 "equation:      $def($eq,NAME)\n\n"

    set Info(messageIcon) error
    ##Message::dispOkMessage [list $m0 $m1 $m2]
    ##return 0
  }

  # Check that no more than three integer values in Dofs
  #
  if { [info exist def($eq,VariableDofs)] } {
    set dofs [split $def($eq,VariableDofs)]
    set err 0

    if { 3 < [llength $dofs] } {
      set err 1
    }

    foreach dof $dofs {
      if { ![string is integer $dof] } {
        set err 1
        break
      }
    }

    if { $err } {
      set m1 "Illegal Variable DOFS value:   $def($eq,VariableDofs)\n\n"
      set m2 "equation:      $def($eq,NAME)\n\n"

      set Info(messageIcon) error
      Message::dispOkMessage [list $m0 $m1 $m2]
      return 0
    }

    # If two values (2D and 3D), add default 1D value 1
    #
    if { 2 == [llength $def($eq,VariableDofs)] } {
      set def($eq,VariableDofs) [linsert $def($eq,VariableDofs) 0 1]
    }
  }

  # Check that no more than three integer values in ExportedVariableDofs
  #
  if { [info exist def($eq,ExportedVariableDofs)] } {

    set checked_dofs ""

    foreach exported_dofs $def($eq,ExportedVariableDofs) {

      set dofs [split $exported_dofs]
      set err 0

      if { 3 < [llength $dofs] } {
        set err 1
      }

      foreach dof $dofs {
        if { ![string is integer $dof] } {
          set err 1
          break
        }
      }

      if { $err } {
        set m1 "Illegal Exported Variable DOFS value:   $exported_dofs\n\n"
        set m2 "equation:      $def($eq,NAME)\n\n"

        set Info(messageIcon) error
        Message::dispOkMessage [list $m0 $m1 $m2]
        return 0
      }

      # If two values (2D and 3D), add default 1D value 1
      #
      if { 2 == [llength $dofs] } {
        set dofs [linsert $dofs 0 1]
      }

      lappend checked_dof $dofs

    }

    set def($eq,ExportedVariableDofs) $checked_dofs
  }


  if { $is_new_eq && ![info exist def($eq,Procedure)] } {

    set m1 "Solver procedure missing!\n\n"
    set m2 "equation:      $def($eq,NAME)\n\n"

    set Info(messageIcon) error
    Message::dispOkMessage [list $m0 $m1 $m2]
    return 0
  }

  if { [info exist def($eq,Procedure)] &&
       2 != [llength $def($eq,Procedure)]
     } {

    set m1 "Illegal Solver procedure definition:   $def($eq,Procedure)\n\n"
    set m2 "equation:      $def($eq,NAME)\n\n"

    set Info(messageIcon) error
    Message::dispOkMessage [list $m0 $m1 $m2]
    return 0
  }

  return 1
}



#==================================#
#   addExportedEquationVariables   #
#==================================#
#
proc UserDefined::addExportedEquationVariables {eq defArray} {
  global Info UserDefined
  upvar $defArray def

  if { ![info exist def($eq,ExportedVariable)] } {
    return
  }
  
  # Add possibly missing (in practice the last) Dofs values
  #
  while { [llength $def($eq,ExportedVariableDofs)] < 
          [llength $def($eq,ExportedVariable)]
        } {
    lappend def($eq,ExportedVariableDofs) 1
  }

  foreach var $def($eq,ExportedVariable) dofs $def($eq,ExportedVariableDofs) {
    DataField::setVariableComponents "" $var $dofs $def($eq,MASK)
  }
}


#================================#
#   setEquationFieldProperties   #
#================================#
#
# NOTE: These properties are set to the equation field-variable
# like NAVIER-STOKES. They are used mainly for InitialCondition and
# BoundaryCondition panels where equation variables (like Velocity)
# are displayed. These properties are not by system (because they
# are already attached to a system!), but many of them can be
# by component (like Velocity and Pressure for Navier-Stokes)
#
# Returns: 1 <--> ok, 0 <--> error
#
proc UserDefined::setEquationFieldProperties { eq defArray ename} {
  global Info UserDefined
  upvar $defArray def
  
  set m0 "ERROR! When reading definition file   \"$UserDefined(currentFileName)\"\n\n"

  #-Simple (non-list ) properties for equation field
  # ================================================
  #
  set equation_properties {
    Variable StatusLineName ProcessTableName
  }

  foreach prop $equation_properties {
    
    if { [info exists def($eq,$prop)] } {

      set val $def($eq,$prop)
      
      if { ![UserDefined::checkFieldPropertyValue $prop $val checked_val] } {
        set m1 "Illegal field property value:  $val \n\n"
        set m2 "property:      $prop \n"
        set m3 "equation:      $def($eq,NAME)\n\n"

        set Info(messageIcon) error
        Message::dispOkMessage [list $m0 $m1 $m2 $m3]
        return 0
      }

      DataField::setFieldProperty Equation $ename $prop $checked_val
    }
  }
  
  #-Field properties which can be component-wise
  # ============================================

  set field_properties ""

  #-These are for the equation field, no prefix!
  #
  lappend field_properties [list "" VariableDofs]
  lappend field_properties [list "" VariableComponents]
  
  # NOTE: !!!
  #-These properties are stored using a property prefix "Var_" , so that they are not mixed
  # with equation fields own properties!
  #
  foreach prop_name $UserDefined(equationVarFieldKeywords) {
    
    if { ![UserDefined::convertPropertyNameToGuiFormat $prop_name gui_prop] } {
      continue
    }

    lappend field_properties [list "Var_" $gui_prop]
  }
  
  # Store property values
  # =====================
  
  foreach prop_info $field_properties {
    
    #-Pick prefix and property name
    set prefix [lindex $prop_info 0]
    set prop [lindex $prop_info 1]

    set value ""

    if { [info exists def($eq,$prop)] } {

      set values [split $def($eq,$prop) ";"]

      foreach val $values {

        set val [Util::stringTRQ $val]

        if { ![UserDefined::checkFieldPropertyValue $prop $val checked_val] } {
          set m1 "Illegal field property value:  $val \n\n"
          set m2 "property:      $prop\n"
          set m3 "equation:      $def($eq,NAME)\n\n"

          set Info(messageIcon) error
          Message::dispOkMessage [list $m0 $m1 $m2 $m3]
          return 0
        }

        append value $checked_val
        append value ";"
      }

      set value [string trim $value ";"]

      #-Store value
      DataField::setFieldProperty Equation $ename $prefix$prop $value
    }
  }

  # EquationVariable related stuff
  # ===============================
  #
  set var [DataField::getFieldProperty Equation $ename Variable "" 0]
  set is_multi_var [DataField::getFieldProperty Equation $ename IsMultiVar]

  if { "" != [DataField::getFieldProperty Equation $ename VariableComponents "" 0] } {
    set has_components 1
  } else {
    set has_components 0
  }

  #-For a multi-var equation without components, the variable name is in use
  # in EquationVariable
  #
  if { $is_multi_var && $var != "" && !$has_components } {

    DataField::setFieldProperty EquationVariable $ename InitialValue $var

    # NOTE: User defined variables are always added to the list of
    # equation's variables, indepenntly what is previously defined.
    # This is different compared to the predefined multi-var equation
    # like AdvectionDiffusion
    #
    DataField::setFieldProperty EquationVariable $ename AlwaysAddInitialValue 1

    #-Create possible variable component names (like "SoundP 1" SoundP 2")
    DataField::setVariableComponents $ename $var

  } else {
    DataField::setFieldProperty EquationVariable $ename InitialValue ""
  }

  #-For a multi-var equation, if no variable components, we do not use any Variable
  # in Equation (compare how AdvDiff is handled!)
  #
  if { $is_multi_var &&
       "" == [DataField::getFieldProperty Equation $ename VariableComponents "" 0]
     } {
    DataField::setFieldProperty Equation $ename Variable ""
  }

  return 1
}



#===================#
#   addPanelField   #
#===================#
#
# Add panel field
#
# Returns: 1 <--> ok, 0 <--> error
#
proc UserDefined::addPanelField { eq panel fld defArray eq_mask eq_index eq_name} {
  global Common Equation Info UserDefined
  upvar $defArray def
  upvar #0 $panel theArray

  set m0 "ERROR! When reading definition file   \"$UserDefined(currentFileName)\"\n\n"

  if { -1 == [lsearch $theArray(initialFields) $fld] } {
    set is_new_fld 1
  } else {
    set is_new_fld 0
  }

  UserDefined::setDefaultFieldProperties $defArray $eq $panel $fld $eq_index

  # Set property values
  # ===================

  #-Some properties defined
  if { [info exist def($eq,$panel,$fld)] } {
    set prop_pairs $def($eq,$panel,$fld)

  #-Only field name given!
  } else {
    set prop_pairs ""
  }

  #-Check and set property values
  #
  foreach prop_pair $prop_pairs {

    set prop [lindex $prop_pair 0]

    set val  [lindex $prop_pair 1]
    set val [Util::stringRemoveQuotes $val]

    #-Check property value
    #
    if { ![UserDefined::checkFieldPropertyValue $prop $val checked_val] } {
      set m1 "Illegal field property value:  $val \n\n"
      set m2 "property:     $prop \n"
      set m3 "field:          [DataField::fieldNameGuiToSif $fld]\n"
      set m4 "panel:        $panel\n"

      if { $eq_name != "" } {
        set m5 "equation:   $eq_name\n\n"
      } else {
        set m5 "\n"
      }
 
      set Info(messageIcon) error
      Message::dispOkMessage [list $m0 $m1 $m2 $m3 $m4 $m5]
      return 0
    }

    # NOTE: Own version of setFieldProperty is used here!
    #    
    UserDefined::setFieldProperty $panel $fld $prop $checked_val $eq_index
    
    # NOTE: For a new field, we set also the non-indexed property value!
    #
    if { $is_new_fld } {
      UserDefined::setFieldProperty $panel $fld $prop $checked_val ""
    }
  }


  # Check some necessary property values
  # ====================================

  #-Widget type
  #
  if { "undef" == [DataField::getFieldProperty $panel $fld WidgetType] } {

    set vt [DataField::getFieldProperty $panel $fld FieldValueType]

    if { $vt == "Logical" } {
      DataField::setFieldProperty $panel $fld WidgetType CheckBox

    } else {
      DataField::setFieldProperty $panel $fld WidgetType Entry
    }
  }
      
  #-Field label is needed
  #
  if { "" == [DataField::getFieldProperty $panel $fld Label] } {
    set lbl [DataField::fieldNameGuiToSif $fld]
    DataField::setFieldProperty $panel $fld Label $lbl
  }

  #-In Equation panel subframe-1 is reserved for the selection box!
  #
  if { $panel == "Equation" } {
    set pf [DataField::getFieldProperty Equation $fld SubPanel]

    if { $pf == 1 } {
      DataField::setFieldProperty Equation $fld SubPanel 2
    }
  }

  # Field mask
  # ==========
  #

  # If field has an activity parent for the equation, we do not add equation mask to the 
  # target mask!!!
  #
  set fld_act_parents [DataField::getFieldProperty $panel $fld ActivityParents $eq_index]

  #-Field belongs to specific equations (not an All equation field)!
  # -----------------------------------
  #
  if { $eq_index != "" } {
        
    set fd_mask [DataField::getFieldMask $panel $fld 0]
    set ft_mask [DataField::getFieldMask $panel $fld 1]
    
    #-Display mask
    #
    if { $fd_mask == "" } {
      set fd_mask $eq_mask

    } elseif  { ![Util::masksAreMatching $eq_mask $fd_mask] } {
      set fd_mask "($fd_mask)"
      append fd_mask "|$eq_mask"
    }
    
    #-Activity mask
    #
    if { $ft_mask == "" && $fld_act_parents == "" } {
      set ft_mask $eq_mask

    } elseif { ![Util::masksAreMatching $eq_mask $ft_mask] } {
      set ft_mask "($ft_mask)"

      if { $fld_act_parents == "" } {
        append ft_mask "|$eq_mask"
      }
    }

    DataField::setFieldProperty $panel $fld Mask [list $fd_mask $ft_mask]

  #-Field belongs to All equations!
  # -------------------------------
  #
  } else {
    set fd_mask ""
    set ft_mask ""
    DataField::setFieldProperty $panel $fld Mask [list $fd_mask $ft_mask]    
  }


  # Some postprocessing for the field
  # =================================

  set wt [DataField::getFieldProperty $panel $fld WidgetType $eq_index]
  set ww [DataField::getFieldProperty $panel $fld WidgetWidth $eq_index]
  set fvt [DataField::getFieldProperty $panel $fld FieldValueType $eq_index]

  set lims [DataField::getFieldProperty $panel $fld Limits $eq_index]
  set ival [DataField::getFieldProperty $panel $fld InitialValue $eq_index]

  #-Evaluate a numeric initial value
  #
  if { $ival != "" &&
       ($fvt == "Real" || $fvt == "Integer" || $fvt == "Logical")
     } {
    
    set msg ""

    if { [DataField::evaluateFieldValue $ival res_ival] } {
      DataField::setFieldProperty $panel $fld InitialValue $res_ival $eq_index
    }
  }

  #-Set option menu properties
  if { $wt == "OptionMenu" } {
    
    # Values list
    if { $lims != "Unlimited" } {
      
      # Add "Set" identifier in front of the limits set, however remove
      # first the existing Set token possibly added by the system
      #
      set lims_set "Set [string trimleft [string trimleft $lims "Set "] ]"
      DataField::setFieldProperty $panel $fld Limits $lims_set $eq_index

      # Initial value (set first list value, if no iv defined)
      if { "" == [DataField::getFieldProperty $panel $fld InitialValue $eq_index] } {
        DataField::setFieldProperty $panel $fld InitialValue [lindex [join $lims] 0] $eq_index
      }

      # Widget width
      if { $ww <= 1 } {

        foreach val [join $lims] {
          set len [string length $val]

          if { $len > $ww } {
            set ww $len
          }
        }
      }
    }
  }


  #-Set browsable file properties
  if { $wt == "BrowsableFile" } {

    # Widget width
    if { $ww <= 1 } {
      set ww 15
    }

    if { $fvt == "Real" } {
      set fvt "File"
    }
  }
  
  
  # Set some (checked) property values
  #
  DataField::setFieldProperty $panel $fld WidgetWidth $ww $eq_index
  DataField::setFieldProperty $panel $fld FieldValueType $fvt $eq_index

  if { $is_new_fld } {
    DataField::setFieldProperty $panel $fld WidgetWidth $ww
    DataField::setFieldProperty $panel $fld FieldValueType $fvt
  }

  # Add field to field lists (if not yet there)
  # ========================
  if { $is_new_fld } {
    lappend theArray(initialFields) $fld
  }

  if { -1 == [lsearch $Common(allFieldNames) $fld] } {
    lappend Common(allFieldNames) $fld
  }

  return 1

};# End addPanelField




#==================#
# Helper functions #
#==================#

# Form "tight" list from list like:
# "Velocity ; Pressure " --> "Velocity;Pressure"
#
proc UserDefined::getVariableComponents {name_list} {
  
  set name_list [split $name_list ";"]

  set names ""

  foreach name $name_list {

    set name [Util::stringTRQ $name]
    
    append names $name
    append names ";"
  }

  return [string trim $names ";"]
}


# Set new value for a field property and check if equation index
# should be used
#
proc UserDefined::setFieldProperty {panel fld property value eq_index} {
  global Common
  upvar #0 $panel thePanel

  # Check if equation index should be used
  #
  if { $eq_index != "" && $eq_index >= 0 } {
   
    # If panel is not by equation or the property is skipped from
    # equationed handling (like WidgetType or Mask)
    #
    if { !$thePanel(panelsByEquation) ||
         -1 != [lsearch $Common(propertiesNotByEquation) $property]
       } {
      set eq_index ""
    }
  }
  
  # These properties need special handling when setting the data
  # (values are lists etc.)
  #
  if { [string equal -nocase $property "ActivityParents"] } {

    set parent_info $value
    
    # Add activity parent info
    DataField::addActivityParent $panel $fld $parent_info $eq_index

    # Add also activity slave info if parent is in the same panel (ie. parent is a field name)
    if { 1 == [llength $parent_info] } {
      DataField::addActivitySlave $panel $parent_info $fld $eq_index
    }

  # Otherwise set propeprty values normally
  #
  } else {
    DataField::setFieldProperty $panel $fld $property $value $eq_index
  }

}


# Set some default property values
#
proc UserDefined::setDefaultFieldProperties {defArray eq panel fld eq_index} {
  global Common
  upvar $defArray def
  upvar #0 $panel theArray

  #-Check is this a new field
  #
  if { -1 == [lsearch $theArray(initialFields) $fld] } {
    set is_new_fld 1
  } else {
    set is_new_fld 0
  }

  # New field: mark some properties undefined, so that we can
  # later set correct default values for them
  #
  if { $is_new_fld } {
    DataField::setFieldProperty $panel $fld WidgetType "undef"
  }

  #-We do not have procedure/table button in these panels (Solver,Constant)!
  #
  if { $panel == $Common(panelArray,SL) || $panel == $Common(panelArray,CN) } {

    lappend def($eq,$panel,$fld) [list Procedure 0]
    lappend def($eq,$panel,$fld) [list Table 0]
    
    # New Solver field: the default is non-packed and inactive!
    #
    if { $panel == $Common(panelArray,SL) } {

      # Generic non-indexed properties (hidden, non-active!)
      if { $is_new_fld && $eq_index != ""} {
        DataField::setFieldProperty Solver $fld Display 0
        DataField::setFieldProperty Solver $fld InitiallyActive 0
      }

      # If equation-id given, make active for the equation!
      if { $eq_index != "" } {
        DataField::setFieldProperty Solver $fld Display 1 $eq_index
        DataField::setFieldProperty Solver $fld InitiallyActive 1 $eq_index

      } else {
        DataField::setFieldProperty Solver $fld Display 1
        DataField::setFieldProperty Solver $fld InitiallyActive 1
      }
    }
  }
}


# Convert user entered property name to Gui format
#
proc UserDefined::convertPropertyNameToGuiFormat {property gui_property} {
  global UserDefined
  upvar $gui_property gui_prop

  foreach pair $UserDefined(propertyKeywordTable) {

    set in_file [lindex $pair 0]
    set in_gui  [lindex $pair 1]
    
    # Remove quotes and all spaces
    set property [Util::stringRemoveQuotes $property]
    set property [string map {" " ""} $property]

    if { [string equal -nocase $in_file $property] } {
      set gui_prop $in_gui
      return 1
    }
  }

  # Error!
  return 0
}


# Check that user has entered an accepted field property value
#
# NOTE: Property values are here already converted from those entered by
# the user to the intenal format!
#
proc UserDefined::checkFieldPropertyValue {property value checked_value} {
  upvar $checked_value checked_val

  # Accepted property values
  # NOTE: Empty accepted value means that a matched value is ok
  # as it is!
  #
  set pairs {
    {WidgetType         { {{CheckBox Entry BrowsableFile Label Separator} ""} {{SelectionList} OptionMenu} }  }
    {FieldValueType     { {{Real Integer Logical String File Procedure} ""} }  }
    {PanelPage          { {{1 2 3 4 5 6 7 8 9 10} ""} }  }
    {SubPanel           { {{1 2 3} ""} }  }
    {IndexType          { {{Prefix Pre} Pre} {{Suffix Suf Post} Post} }  }
    {LabelByCoordinate  { {{1 True} 1} {{0 False} 0} }  }
    {Procedurable       { {{1 True} 1} {{0 False} 0} }  }
    {Tableable          { {{1 True} 1} {{0 False} 0} }  }
    {Variabled          { {{1 True} 1} {{0 False} 0} }  }
    {Coordinate         { {{1 True} 1} {{0 False} 0} }  }
    {Display            { {{1 True} 1} {{0 False} 0} }  }
    {HasAbsCheckBox     { {{1 True} 1} {{0 False} 0} }  }
    {AlwaysOutput       { {{1 True} 1} {{0 False} 0} }  }
    {DropInitialValued  { {{1 True} 1} {{0 False} 0} }  }
    {FieldTarget        { {{Solver} ""} }  }
  }

  # Check if property is listed in checkable property value list
  foreach pair $pairs {

    set prop [lindex $pair 0]

    if { $property != $prop } {
      continue
    }
    
    # Check property value
    # ====================

    #-These property value cannot have internal blanks!
    set val [Util::stringRemoveQuotes $value]
    set val [Util::makeArrayName $val]

    #-Accepted values and their qui equivalent
    set check_pairs [lindex $pair 1]

    #-Loop all accepted values list
    foreach check_pair $check_pairs {

      set accepted_vals [lindex $check_pair 0]

      set idx [lsearch $accepted_vals $val]

      #-If value is in the list
      if { $idx != -1 } {
  
        set gui_val [lindex $check_pair 1]
        
        # If gui-value is empty, we use "cleaned" value itself
        if { $gui_val == "" } {
          set checked_val $val

        #-Otherewise we use given gui-value
        } else {
          set checked_val $gui_val
        }

        return 1
      }
    }
    
    #-No accepted value was found for the property
    return 0
  }

  # These properties need some special checking/processing
  # ======================================================
  #

  #-Activity parent can be a single field or a list of arr+field
  #
  if { [string equal -nocase "ActivityParents" $property] } {
    set value [Util::stringTrim $value]

    # If given as a list ( {Equation "Field NAme"}
    #
    if { "\{"== [string index $value 0] } {
      set value [join $value]
      set arr_value [lindex $value 0]
      set fld_value [lindex $value 1]

    # Single field name
    #
    } else {
      set arr_value ""
      set fld_value $value
    }
    
    set value ""

    set arr_value [Util::stringTrim $arr_value]
    set fld_value [Util::stringTrim $fld_value]

    # Store as a list
    #
    if { $arr_value != "" } {
      lappend value $arr_value
      lappend value [DataField::fieldNameSifToGui $fld_value]

    # Store as a single field name
    #
    } else {
      set value [DataField::fieldNameSifToGui $fld_value]
    }
  }


  # This property was not in the checkable list, value
  # is always ok!
  #
  set checked_val $value

  return 1
}


# Read one predefined keyword , it must be in the "accepted_keywords" list
#
# Updates argument variables "kword" and "kvalue"
# Returns: 1 --> Ok, 0 --> Not an accepted keyword
#
proc UserDefined::readAcceptedKeyword { accepted_keywords line kword kvalue } {
  upvar $kword kwrd
  upvar $kvalue kval
  
  foreach kw $accepted_keywords {

    set val ""

    if { [UserDefined::readKeyword $line $kw val] } {
      set kwrd $kw
      set kval [Util::stringTrim $val]
      set kval [Util::stringRemoveQuotes $kval 0]

#MSG "kwd=$kw line=$line kval=$kval"

      return 1
    }
  }

#MSG "Illegal keyword $kwrd"
  
  return 0
}


# Try to read one keyword entry named in argument "keyword"
# 
# Updates argument variable "value"
# Returns: 1 -> Ok, 0 --> No keyword found
#
proc UserDefined::readKeyword {line keyword value} {
  upvar $value kval

  set kval ""

  # Try first something separted by equal (=) sign
  #
  if { -1 != [string first "=" $line] } {
    set tmp [split $line "="]
    set kwd [string trim [lindex $tmp 0]]
    set res_value [string trim [lindex $tmp 1] ]

    if { [string equal -nocase $keyword $kwd] } {
      set kval $res_value
      return 1
    }
  }

  # Try next the whole keyword in quotes
  #
  if { "\"" == [string index $line 0] } {
    set tmp [string range $line 1 end]
    set idx [string first "\"" $tmp]
    set idx1 [expr $idx - 1]
    set idx2 [expr $idx + 1]
    set kwd [string trim [string range $tmp 0 $idx1] ]
    set res_value [string trim [string range $tmp $idx2 end] ]

    if { [string equal -nocase $keyword $kwd] } {
      set kval $res_value
      return 1
    }
  }

  # Try next whole keyword, but without any quotes
  #
  if { 1 == [regsub -nocase $keyword $line "" res_value] } {

    # NOTE: We do not accept partial match here!
    # (Like "Variablee" for "Variable"
    #
    if { $res_value != "" &&
         " " != [string index $res_value 0]
       } {
      return 0
    }

    set kval $res_value
    return 1
  }

  # Keyword is single word, so no chance!
  #
  if { 1 == [llength $keyword] } {
    return 0
  }

  # Next try each piece of keyword (in case of extra blanks in the line!)
  foreach kwp $keyword {
    
    if { ![regsub -nocase $keyword $res_value "" res_value] } {
      return 0
    }
  }

  set kval $res_value

#MSG "ok: $keyword   $kval"

  return 1
}


# Store definition information to global arries for
# later reference
#
proc UserDefined::addToDefinitionLog { line {indent_level 0} } {
  global Model

  set is [expr 2 * $indent_level]
  set indent [string repeat " " $is]

  lappend Model(newDefinitions) $indent$line
}


############################
# SOLVER.KEYWODSDS file stuff
############################

# Load all possible (updated) SOLVER.KEYWORD files
#
proc UserDefined::loadSolverKeywordsFiles {} {
	global Info Model

	#--Possible ELMER_HOME/lib version
	set kwd_file "$Info(ELMER_HOME)/lib/SOLVER.KEYWORDS"

	if { [file exists $kwd_file] } {
		set mtm [file mtime $kwd_file]
		if { $kwd_file != $Info(file_SKWD_lib) || $mtm != $Info(mtime_SKWD_lib) } {
			set Info(file_SKWD_lib) $kwd_file
			set Info(mtime_SKWD_lib) $mtm
			UserDefined::readSolverKeywordsFile $kwd_file "lib"
		}
	}

	#--Possible ELMER_FRONT_PREFIX/lib version
	set kwd_file "$Info(ELMER_FRONT_INSTALL_LIB)/SOLVER.KEYWORDS"

	if { [file exists $kwd_file] } {
		set mtm [file mtime $kwd_file]
		if { $kwd_file != $Info(file_SKWD_lib) || $mtm != $Info(mtime_SKWD_lib) } {
			set Info(file_SKWD_lib) $kwd_file
			set Info(mtime_SKWD_lib) $mtm
			UserDefined::readSolverKeywordsFile $kwd_file "lib"
		}
	}

	#--Possible ELMER_USER_HOME version
	set kwd_file "$Info(ELMER_USER_HOME)/SOLVER.KEYWORDS"

	if { [file exists $kwd_file] } {
		set mtm [file mtime $kwd_file]
		if { $kwd_file != $Info(file_SKWD_user) || $mtm != $Info(mtime_SKWD_user) } {
			set Info(file_SKWD_user) $kwd_file
			set Info(mtime_SKWD_user) $mtm
			UserDefined::readSolverKeywordsFile $kwd_file "user"
		}
	}

	#--Possible model's directory version
	if { $Model(EMF_FILE) != "" } {
		set kwd_file "$Model(EMF_OPEN_DIR)/SOLVER.KEYWORDS"
		if { [file exists $kwd_file] } {
			set mtm [file mtime $kwd_file]
			if { $kwd_file != $Info(file_SKWD_model) || $mtm != $Info(mtime_SKWD_model) } {
				set Info(file_SKWD_model) $kwd_file
				set Info(mtime_SKWD_model) $mtm
				UserDefined::readSolverKeywordsFile $kwd_file "model"
			}
		}
	}
}


# Read one SOLVER.KEYWORDS file
#
proc UserDefined::readSolverKeywordsFile { filename {type ""} } {
  global Info SKWD
	  
	set sep ":"

	Util::unsetArray SKWD

  # =========
  # Open file
  # =========

  set cmsg ""
  if { [ catch {set in [open $filename "r"]} cmsg ] } {

    Message::showMessage ""
    Message::showMessage "WARNING*** $cmsg" $Info(wrnMsgColor)
    Message::showMessage ""
    
    return 0
  }

  #Message::showMessage "Reading SOLVER.KEYWORDS file:  $filename" $Info(remMsgColor)
  
  # Read type definitions
  # =====================
  while {1} {

    set line [Util::readInputFileLine $in]

    set line [Util::stringTrim $line]

    set err_line "\nin the input line:   $line\n\n"

    # No more data (eof)
    # ------------------
    if { $line == "" } {
      break
    }
    
    # Read next keyword
    # -----------------
		set section ""
    set keyword ""
    set type ""
		set sec ""
		
		set items [split $line $sep]

		set section [string tolower [lindex $items 0]]

		if { 3 == [llength $items] } {
			set type [string tolower [lindex $items 1] ]
		}

		# Illegal type!
		#
		if { $type != "real" &&
		     $type != "integer" &&
		     $type != "logical" &&
		     $type != "string" &&
		     $type != "file" &&
		     $type != "procedure"
			 } {
			continue
		}

		set keyword [string tolower [Util::stringTRQ [lindex $items end]]]

		switch $section {

			"body" {set sec "BO"}
			"bodyforce" {set sec "BF"}
			"boundary" {set sec "BA"}
			"bc" {set sec "BC"}
			"constants" {set sec "CN"}
			"equation" {set sec "EQ"}
			"ic" {set sec "IC"}
			"material" {set sec "MT"}
 			"simulation" {set sec "SI"}
			"solver" {set sec "SL"}
		}

#MSG "$section $sec $type $keyword"
		
		# Unknown section!
		#
		if { $sec == "" } {
			continue
		}

		if { $type == "" } {
			set SKWD($sec,$keyword) 0
		} else {
			set SKWD($sec,$keyword) 1
		}

  } ;#while(1)

	close $in 

  Message::showMessage "***SOLVER.KEYWORDS file loaded:  $filename" $Info(remMsgColor)

  # Store filename in info-list
  Util::addToFileList Info(solverKeywordFiles) $filename

	return 1
}


# end ecif_tk_procsUserDefined.tcl
# *************************
