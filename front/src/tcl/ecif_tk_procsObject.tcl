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
#Module:    ecif_tk_procsObject.tcl
#Language:  Tcl
#Date:      20.12.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Object handling procs
#
#************************************************************************

# ========= #
# GET procs #
# ========= #


proc Object::getIds { {object_type ""} {object_classes ""} } {
  global ObjectTable

  set result ""

  set otype [string tolower $object_type]

  if { $otype == "" || $otype == "a" || $otype == "all"} {
    return $ObjecTable(ids)

  } elseif { $otype == "b" || $otype == "body" } {
    set type "B"

  } elseif { $otype == "bl" || $otype == "bodylayer" } {
    set type "BL"

  } elseif { $otype == "bp" || $otype == "bodypair" } {
    set type "BP"

  } elseif { $otype == "bg" || $otype == "boundarygroup" } {
    set type "BG"

  } elseif { $otype == "fg" || $otype == "facegroup" } {
    set type "FG"

  } elseif { $otype == "eg" || $otype == "edgegroup" } {
    set type "EG"

  } elseif { $otype == "vg" || $otype == "vertexgroup" } {
    set type "VG"

  } elseif { $otype == "f" || $otype == "face" } {
    set type "F"

  } elseif { $otype == "e" || $otype == "edge" } {
    set type "E"

  } elseif { $otype == "v" || $otype == "vertex" } {
    set type "V"

  } else { 
    MSG "WARNING: Illegal object type: $object_type"
    return ""
  }

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != $type } {
      continue
    }
    
    set found 1

    if { $object_classes != "" } {

      set found 0
      foreach cls $object_classes {

        # Negation condition
        # Like: -ILB meaning not an intra layer boundray
        #
        if { "-" == [string index $cls 0] } {
          if { $ObjectTable($id,cl) != $cls } {
            set found 1
            break
          }
        
        # Normal condition
        #
        } else {
          if { $ObjectTable($id,cl) == $cls } {
            set found 1
            break
          }
        }
      }
    }

    if { $found == 0 } {
      continue
    }

    lappend result $id
  }

  return $result
}


proc Object::getEquationIds {oid} {
  global ObjectTable

  if { ![info exists ObjectTable($oid,tp)] } {
    return ""
  }

  set eids ""
  set tp $ObjectTable($oid,tp)

  if { $tp == "B" } {

    if { [info exists ObjectTable($oid,eq)] } {
      set eids $ObjectTable($oid,eq)
    }

  } elseif { $tp == "BP" } {

    # Corrected: MVe 02.06.2004
    #set oid1 $ObjectTable($oid,bd1Id)
    #set oid2 $ObjectTable($oid,bd2Id)
    set oid1 $ObjectTable($oid,pr1Id)
    set oid2 $ObjectTable($oid,pr2Id)
    
    if { [info exists ObjectTable($oid1,eq)] } {
      lappend eids $ObjectTable($oid1,eq)
    }

    if { [info exists ObjectTable($oid2,eq)] } {
      lappend eids $ObjectTable($oid2,eq)
    }

  }

  return $eids
}


proc Object::getDisplayed {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,dspl)] } {
    return $ObjectTable($id,dspl)

  } else {
    return false
  }
}


proc Object::getClass {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,cl)] } {
    return $ObjectTable($id,cl)

  } else {
    return ""
  }
}


proc Object::getColor {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,clr)] } {
    return $ObjectTable($id,clr)

  } else {
    return ""
  }
}



proc Object::getMask {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,msk)] } {
    return $ObjectTable($id,msk)

  } else {
    return ""
  }
}


proc Object::getTypeMask {id} {
  global Info ObjectTable

  if { ![info exists ObjectTable($id,tp)] } {
    return ""

  } else {
    set m ""
    append m $Info(maskStart)
    append m $id
    append m $ObjectTable($id,tp)

    return $m
  }
}


proc Object::getName {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,nm)] } {
    return $ObjectTable($id,nm)

  } else {
    return ""
  }

}


proc Object::getParent {id} {
  global Info ObjectTable

  if { [info exists ObjectTable($id,prId)] } {
    return $ObjectTable($id,prId)

  } else {
    return $Info(NO_INDEX)
  }

}


proc Object::getParents {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,pr1Id)] &&
       [info exists ObjectTable($id,pr2Id)]
     } {
    return "$ObjectTable($id,pr1Id) $ObjectTable($id,pr2Id)"

  } else {
    return ""
  }
}


proc Object::getParameterId {id pv} {
  global Info ObjectTable

  if { [info exists ObjectTable($id,$pv)] } {
    return $ObjectTable($id,$pv) 

  } else {
    return $Info(NO_INDEX)
  }
}


proc Object::getSelected {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,slctd)] } {
    return $ObjectTable($id,slctd)

  } else {
    return false
  }
}


proc Object::getSubIds {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,sbIds)] } {
    return $ObjectTable($id,sbIds)

  } else {
    return ""
  }
}


proc Object::getSubTags {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,sbTgs)] } {
    return $ObjectTable($id,sbTgs)

  } else {
    return ""
  }
}


proc Object::getTag {id} {
  global Info ObjectTable

  if { [info exists ObjectTable($id,tg)] } {
    return $ObjectTable($id,tg)

  } else {
    return $Info(NO_INDEX)
  }
}


proc Object::getType {id} {
  global ObjectTable

  if { [info exists ObjectTable($id,tp)] } {
    return $ObjectTable($id,tp)

  } else {
    return ""
  }
}


# ========= #
# SET procs #
# ========= #

proc Object::setClass {id value} {
  global ObjectTable

  set ObjectTable($id,cl) $value
}


proc Object::setColor {id value} {
  global ObjectTable

  set ObjectTable($id,clr) $value
}


proc Object::setDisplayed {id value} {
  global ObjectTable

  set ObjectTable($id,slctd) $value
}


proc Object::setMask {id value} {
  global ObjectTable

  set ObjectTable($id,msk) $value
}


proc Object::setName {id value} {
  global ObjectTable

  set ObjectTable($id,nm) $value
}


proc Object::setParent {id value} {
  global ObjectTable

  set ObjectTable($id,prId) $value
}


proc Object::setParents {id value1 value2} {
  global ObjectTable

  set ObjectTable($id,pr1Id) $value1
  set ObjectTable($id,pr2Id) $value2
}


proc Object::setSelected {id value} {
  global ObjectTable

  set ObjectTable($id,slctd) $value
}


proc Object::setSubIds {id value} {
  global ObjectTable

  set ObjectTable($id,sbIds) $value
}


proc Object::setSubTags {id value} {
  global ObjectTable

  set ObjectTable($id,sbTgs) $value
}


proc Object::setType {id value} {
  global ObjectTable

  set ObjectTable($id,tp) $value
}


# ===========
# OTHER procs
# ===========


proc Object::findIdByTag { otype tags } {
  global Info ObjectTable

  foreach otp $otypes {
    foreach id $ObjectTable(ids) {

      if { $ObjectTable($id,tp) == $otp &&
           $ObjectTable($id,tg) == $tag
         } {
        return $id
      }
    }
  }

  # Not found!
  return $Info(NO_INDEX)
}


proc Object::findIdByParentIds { otype id1 id2 } {
  global Info ObjectTable

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) == $otype  &&
         $ObjectTable($id,pr1Id) == $id1 &&
         $ObjectTable($id,pr2Id) == $id2
       } {
      return $id
    }
  }

  # Not found!
  return $Info(NO_INDEX)
}

# Check if a boundary is between layers in a same body
#
# NOTE: These boundaries are not shown in BC-panel, tehy are
# relevant only for meshing!
#
proc Object::isIntraLayerBoundary {id } {
  global ObjectTable Model

  set tp $ObjectTable($id,tp)

  if { $Model(GEOMETRY_DIMENSION) == "3D" } {
    if { $tp != "F" && $tp != "FG" } {
      return 0
    }
  } elseif { $Model(GEOMETRY_DIMENSION) == "2D" } {
    if { $tp != "E" && $tp != "EG" } {
      return 0
    }
  }

  set prId $ObjectTable($id,prId)
  set blIds $ObjectTable($prId,blIds)

  if { 2 != [llength $blIds] } {
    return 0
  }

  foreach blId $blIds {
    if { -1 == [lsearch $ObjectTable($blId,sbIds) $id] } {
      return 0
    }
  }

  return 1
}


proc Object::parameterIsApplied {pid pv object_types } {
  global ObjectTable

  foreach id $ObjectTable(ids) {

    foreach tp $object_types {

      if { $ObjectTable($id,tp) == $tp &&
           $pid == [Object::getParameterId $id $pv]
         } {
        return 1
      }
    }
  }

  return 0
}


proc Object::setBodyPairNames { } {
  global Info ObjectTable

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "BP" } continue
    
    set nm ""
    append nm $ObjectTable($ObjectTable($id,pr1Id),nm)
    append nm "-"
    append nm $ObjectTable($ObjectTable($id,pr2Id),nm)

    set ObjectTable($id,nm) $nm
  }
}


# Mark the boundary selection status in the
# ObjectTable(id,slctd) flags
#
proc Object::setBoundarySelectionMode {bndr_id {mode ""}} {
  global Common ObjectTable

  set glob_array $Common(currentArray)
  set trg_id $bndr_id

  # If currently no active panel
  if {$glob_array == ""} {
    return ""
  }

  if { $glob_array == $Common(panelArray,BC) } {
    if { [info exists ObjectTable($bndr_id,grpId)] } {
      set trg_id $ObjectTable($bndr_id,grpId)
    }
  }

  global $glob_array
  upvar #0 $glob_array theArray

  #-Toggle old mode value
  if { $mode == "" } {

    set value $ObjectTable($trg_id,slctd) 

    if { $value == 1 } {
      set new_value 0

    } else {
      set new_value 1
    }

  #-Mode set explicitely
  } else {
    set new_value $mode
  }

  set ObjectTable($trg_id,slctd) $new_value
}


proc Object::sortIdsByParentTags { ids } {
  global Info ObjectTable

  set ids_table ""

  foreach id $ids {

    set row $id

    if { [info exist ObjectTable($id,pr1Id)] } {
      lappend row $ObjectTable($ObjectTable($id,pr1Id),tg)
    } else {
      lappend row $Info(NO_INDEX)
    }

    if { [info exist ObjectTable($id,pr2Id)] } {
      lappend row $ObjectTable($ObjectTable($id,pr2Id),tg)
    } else {
      lappend row $Info(NO_INDEX)
    }

    lappend ids_table $row

  }

  # Parent2 is the secondary key
  set ids_table [lsort -integer -index 2 $ids_table]

  # Parent1 is the primary key
  set ids_table [lsort -integer -index 1 $ids_table]

  set result ""

  foreach item $ids_table {

    lappend result [lindex $item 0]
  }

  return $result
}


# NOTE: key_entries is list of {var-name var-code-proc} entries
# where var-code-proc may be empty
#
proc Object::sortIdsByKeyVars { ids key_entries numeric_key_flags } {
  global Info ObjectTable

  if { $key_entries == "" } {
    return $ids
  }


  #-Form list of rows where we have the id, key1, key2, etc
  set ids_table ""
  foreach id $ids {

    set row $id

    foreach key_entry $key_entries {

      set key_var [string tolower [lindex $key_entry 0]]

      set key_code_proc [lindex $key_entry 1]

      if { [info exist ObjectTable($id,$key_var)] } {

        set value $ObjectTable($id,$key_var)

        if { $key_code_proc != "" } {
          lappend row [eval $key_code_proc $ObjectTable($id,$key_var)]
        } else {
          lappend row $ObjectTable($id,$key_var)
        }

      } else {
        lappend row $Info(NO_INDEX)
      }
     
    }

    lappend ids_table $row
  }

  #-Sort by each key, starting from the lowest level
  # until to the primary (first) key-var
  set index [llength $key_entries]

  while { $index > 0 } {

    set is_numeric [lindex $numeric_key_flags [expr $index - 1]]

    if { $is_numeric } {
      set ids_table [lsort -integer -index $index $ids_table]

    } else {
      set ids_table [lsort -index $index $ids_table]
    }

    incr index -1
  }

  set result ""

  #--Pick finally the ids
  foreach item $ids_table {
    lappend result [lindex $item 0]
  }

  return $result
}


# Recode boundary object type for sorting
#
proc Object::codeBoundaryTypeForSorting {tp} {

  # Face
  if { $tp == "F" || $tp == "FG" } {
    set tp "a"
  # Edge
  } elseif { $tp == "E" || $tp == "EG" } {
    set tp "b"
  # Vertex
  } else {
    set tp "c"
  }

  return $tp
}



