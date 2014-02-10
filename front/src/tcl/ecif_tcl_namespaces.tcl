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
#Module:    ecif_tcl_namespaces.tcl
#Language:  Tcl
#Date:      28.06.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Script defines all namespaces for Front Tcl/Tk
#
#************************************************************************

# Executes an 'namespace' proc like: InitialCondition::objectSelected
# Returns 1 if proc found, otherwise returns 0
#
proc execNsProc {nsName procName {arguments ""}} {

  set ns_proc_name $nsName
  append ns_proc_name "::$procName"

  #-If proc exists
  if {"" != [info commands $ns_proc_name] } {

    if { $arguments == "" } {
      eval $ns_proc_name
    } else {
      eval $ns_proc_name $arguments
    }
    return 1

  #-Otherwise
  } else {
MSG "Proc $ns_proc_name not found!"
    return 0
  }
}


# =================================
# Namespace variable handling procs
# =================================

proc getNsVar {NS VN} {
  set n $NS::$VN
  return [set $n]
}


proc getNsIdVar {NS ID VN} {
  set n $NS::$ID,$VN
  return [set $n]
}


proc setNsIdVar {NS ID VN val} {
  set n $NS::$ID,$VN
  set rc [eval { set $n $val }]
}

proc setNsVar {NS VN val} {
  set n $NS::$VN
  set rc [eval { set $n $val }]
}


proc updateNsListRow { NS VN row_index new_row} {
  
  set old_list [getNsVar $NS $VN]
  set new_list [lreplace $old_list $row_index $row_index $new_row]
  setNsVar $NS $VN $new_list
}
    

proc setNsIdVar {NS ID VN val} {
  set n $NS::$ID,$VN
  set rc [eval { set $n $val }]
}


# ==========
# Namespaces
# ==========

namespace eval About {
  namespace export openPanel
  namespace export panelOk
}

namespace eval BodyDisplay {
  namespace export all
  namespace export none
  namespace export select
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export setDisplayMode
}

namespace eval BodyForce {
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export parameterDblSelected
  namespace export parameterSelected
}

namespace eval BodyParameter {
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export parameterDblSelected
  namespace export parameterSelected
}

namespace eval BodyInfo {
  namespace export openPanel
  namespace export panelCancel
  namespace export panelOk
  namespace export setBodyColorEntry
  namespace export setPanelContents
}

namespace eval BodyProperty {
  namespace export applyBodyNameEntryContents
  namespace export getColorName
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export panelSave
  namespace export selectNextBody
  namespace export setBodyColor
  namespace export setBodyName
  namespace export setColorContents
  namespace export setNameAndColorContents
  namespace export setNameContents
}

namespace eval Boundary {
  namespace export applyBoundaryNameEntryContents
  namespace export boundaryCntrlDblSelected
  namespace export boundaryDblSelected
  namespace export boundarySelected
  namespace export combine
  namespace export createButtons1Area
  namespace export createButtons2Area
  namespace export createButtons3Area
  namespace export createButtons4Area
  namespace export createNamesArea
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export openPanel
  namespace export packButtons1Area
  namespace export packButtons2Area
  namespace export packButtons3Area
  namespace export packButtons4Area
  namespace export packNamesArea
  namespace export packPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export panelSave
  namespace export panelUpdate
  namespace export saveNames
  namespace export selectMeshBoundaryElements
  namespace export setEditBoundaries
  namespace export setNextActiveSelectionToleranceInfo
  namespace export setSelectMethodAll
  namespace export setSelectMethodByNeighbor
  namespace export setSelectMethodByNormal
  namespace export setSelectMethodByPlane
  namespace export setSelectModeExtend
  namespace export setSelectModeReduce
  namespace export setSelectModeToggle
  namespace export splitCombineRedo
  namespace export splitCombineUndo
  namespace export splitCombineUpdate
  namespace export split_
}

namespace eval BoundaryCondition {
  namespace export boundaryCntrlDblSelected
  namespace export boundaryDblSelected
  namespace export boundarySelected
  namespace export cancelProc
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export okProc
  namespace export parameterDblSelected
  namespace export parameterSelected
  namespace export parameterTargetCheckProc
  namespace export postInitData
  namespace export restoreProc
  namespace export updateModelStatus
  namespace export updateRadiationTargetInfo
}

namespace eval BoundaryDisplay {
  namespace export all
  namespace export none
  namespace export select
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export setBoundaryDisplayMode
  namespace export setDisplayMode
}

namespace eval BoundaryParameter {
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export parameterDblSelected
  namespace export parameterSelected
}

namespace eval CheckBoxList {
  variable id 0
  namespace export all
  namespace export apply
  namespace export cancel
  namespace export close
  namespace export create
  namespace export fillListBox
  namespace export none
  namespace export ok
  namespace export save
}

namespace eval Calculator {
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
}

namespace eval Constant {
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
}

namespace eval Coordinate {
  namespace export applyCoordinateDimension
  namespace export checkPanelData
  namespace export getLabelIndex
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
  namespace export setCoordinateStatus
  namespace export setSymmetryPlaneStatus
  namespace export updateLabels
}

namespace eval DataField {
  namespace export addFieldMaskItem
  namespace export calcFieldSize
  namespace export checkFileName
  namespace export checkFloat
  namespace export checkFunctionName
  namespace export checkInteger
  namespace export checkModuleName
  namespace export checkNumber
  namespace export checkNumbers
  namespace export checkParameterData
  namespace export checkProcedureName
  namespace export checkSizeParameter
  namespace export checkString
  namespace export checkValue
  namespace export checkVariableName
  namespace export checkVector
  namespace export clearLBSelection
  namespace export compressParameterInputData
  namespace export compressParameterOutputData
  namespace export copyFieldProperties
  namespace export data2TextWidget
  namespace export evalInitialValue
  namespace export extractFieldSize
  namespace export extractFieldSizeNew
  namespace export extractVariableAndSizeDescription
  namespace export fieldDisplayMaskMatches
  namespace export fieldProblemMaskMatches
  namespace export fieldTargetMaskMatches
  namespace export fillFieldFromDispList
  namespace export fillGroupFromDispList
  namespace export formDataFields
  namespace export formDefaultNonStandardParameter
  namespace export formFieldParameterLineData
  namespace export formGroupParameterLineValue
  namespace export formInitialParameterLineValue
  namespace export formNonStandardParameter
  namespace export formParameterLine
  namespace export formSingleParameterLineValue
  namespace export formVariableAndSizeDescription
  namespace export formatNumbers
  namespace export formatParamId
  namespace export getFieldDataAsList
  namespace export getFieldMask
  namespace export getFieldProperty
  namespace export getFieldSizes
  namespace export getFieldValue
  namespace export getInitialValue
  namespace export getInitiallyActive
  namespace export getNavierStokesDOFs
  namespace export getStressAnalysisDOFs
  namespace export isValidVariable
  namespace export pickVariableAndSizeDescription
  namespace export selectArrayFields
  namespace export setAndShowLBSelection
  namespace export setDataFields
  namespace export setFieldProperties1
  namespace export setFieldProperties2
  namespace export setFieldProperty
  namespace export setFieldValue
  namespace export setVariableAndSizeValues
  namespace export splitProcedureEntryValue
  namespace export textWidget2Data
  namespace export trimStringGroupData
  namespace export updatePanelValuesArea
}

namespace eval Datafile {
  namespace export checkRestartFile
  namespace export checkRestartPosition
  namespace export getDefaultMeshInputFile
  namespace export getDefaultOutputFile
  namespace export getDefaultPostFile
  namespace export getDefaultSolverInputFile
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export panelSave
}

namespace eval Equation {
  namespace export addEquationIndex
  namespace export addProcPost
  namespace export addProcPre
  namespace export attachProc
  namespace export cancelProc
  namespace export constructProblemMask
  namespace export deleteProc
  namespace export detachProc
  namespace export getMatchingEquationIndices
  namespace export initMasks
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export okProc
  namespace export parameterDblSelected
  namespace export parameterSelected
  namespace export postInitData
  namespace export setEquationAndProblemMasks
  namespace export updateEquationIndices
  namespace export updateTargetMask
  namespace export updateProcPost
  namespace export updateProcPre
}

namespace eval EquationVariable {
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export panelSave
}

namespace eval FileBrowser {
  variable id 0
  namespace export browse
  namespace export closed
  namespace export monitor
  namespace export readHandler
  namespace export refresh
  namespace export save
  namespace export toggleUpdateMode
}

namespace eval GridH {
  namespace export boundaryCntrlDblSelected
  namespace export boundaryDblSelected
  namespace export boundarySelected
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export okProc
  namespace export parameterDblSelected
  namespace export parameterSelected
  namespace export preInitData
  namespace export updateObjectParameterId
  namespace export updateObjectParameterIds
}

namespace eval GridParameter {
  namespace export checkPanel
  namespace export checkParameterAttachment
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export okProc
  namespace export parameterDblSelected
  namespace export parameterSelected
  namespace export postInitData
  namespace export preInitData
  namespace export setGroupState
  namespace export updateObjectBoundariesWdg
  namespace export updateObjectParameterId
  namespace export updateObjectParameterIds
}

namespace eval InitialCondition {
  namespace export objectDblSelected
  namespace export objectDblSelected2
  namespace export objectSelected
  namespace export parameterDblSelected
  namespace export parameterSelected
}

namespace eval InputFileInfo {
  namespace export openPanel
  namespace export panelOk
}

namespace eval Interface {
  namespace export activateFileMenus
  namespace export applyAllParameterData
  namespace export applyBodyData
  namespace export applyBoundaryData
  namespace export applyModelData
  namespace export applyModelFlags
  namespace export applyModelGeometryDimension
  namespace export applyModelStatus
  namespace export applyObjectTableData
  namespace export applyStatsData
  namespace export checkmeshInfoTs
  namespace export configureButtonOption
  namespace export configureButtons
  namespace export configureMenuButton
  namespace export configureMenuButtons
  namespace export createBodiesMenus
  namespace export debug
  namespace export displayBodySelectPanel
  namespace export displayBoundarySelectPanel
  namespace export displayErrorMsg
  namespace export displayLabelSelectPanel
  namespace export getColorName
  namespace export getNextActiveSelectionTolerance
  namespace export getParallelInfo
  namespace export initEquationedPanels
  namespace export markSelectedBoundaries
  namespace export printFlags
  namespace export rendererReset
  namespace export resetObjectTable
  namespace export saveModelPropertyData
  namespace export selectBody
  namespace export selectBoundaries
  namespace export selectBoundary
  namespace export selectBoundary2
  namespace export setBoundarySelectionMode
  namespace export setCurrentMeshF
  namespace export setCurrentMeshH
  namespace export setCurrentTimestamp
  namespace export setInitialMeshH
  namespace export setInitialState
  namespace export setMeshEdited
  namespace export setMeshExists
  namespace export setModelFilePath
  namespace export setModelHasElmerMesh
  namespace export setModelHasMeshParameter
  namespace export setNeedsUpdate
  namespace export setParameterFieldValueState
  namespace export setRotatePriority
  namespace export setStatusBodyForces
  namespace export setStatusEquations
  namespace export setStatusInitialConditions
  namespace export setStatusInnerBoundaryConditions
  namespace export setStatusMaterials
  namespace export setStatusMeshes
  namespace export setStatusOuterBoundaryConditions
  namespace export setStatusTimestamps
  namespace export setWasUpdated
  namespace export showMessage
}

namespace eval LabelDisplay {
  namespace export all
  namespace export none
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export setDisplayMode
}

namespace eval List {
  namespace export addToIdList
  namespace export addToIndexedIdList
  namespace export getListRowsByIndex
  namespace export getListRowsByRange
  namespace export getUniqueIdList
  namespace export makeIndicesFromRange
  namespace export removeFromIdList
  namespace export removeFromIndexedIdList
}

namespace eval ListBox {
  namespace export fill
  namespace export markParameterLBByMask
  namespace export markSelectedBoundaries
  namespace export selectBody
  namespace export selectBoundary
  namespace export setBoundarySelectionMode
  namespace export toggleSelectionMark
  namespace export updateRow
}


namespace eval MatcDefinitions {
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export selectMatcDef
}
namespace eval Material {
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export parameterDblSelected
  namespace export parameterSelected
  namespace export setValuesAreaStates
}

namespace eval MenuBuild {
  namespace export DisplayMenu
  namespace export EditMenu
  namespace export FileMenu
  namespace export HelpMenu
  namespace export MeshMenu
  namespace export ModelMenu
  namespace export ProblemMenu
  namespace export RunMenu
  namespace export SolverMenu
  namespace export WindowMenu
  namespace export addMenuBarItem
  namespace export addMenuCascade
  namespace export addMenuCommand
  namespace export applyUserLevel
  namespace export configureButtonOption
  namespace export configureButtons
  namespace export configureFileLoadMeshMenus
  namespace export configureMenuButtonOption
  namespace export configureMenuButtons
  namespace export configurePanelState
  namespace export createBodyListMenu
  namespace export createMenuBarItem
  namespace export createWindowListMenu
  namespace export deleteMenuBarItem
  namespace export initMenuInfo
  namespace export insertMenuBarItem
}

namespace eval MenuExec {
  namespace export applySettings
  namespace export applySettingsToCaseDirectories
  namespace export askModelOk
  namespace export browseLogFile
  namespace export browseModelFile
  namespace export browseProcessLog
  namespace export browseSolverInputFile
  namespace export case_info
  namespace export checkCadData
  namespace export checkDBDirectory
  namespace export checkDBName
  namespace export checkDirectory
  namespace export checkEquationMeshes
  namespace export checkGebhardtFactors
  namespace export checkMeshData
  namespace export checkModelFile
  namespace export checkOpenPanels
  namespace export checkPostFile
  namespace export checkProcessEnvironment
  namespace export checkSolverData
  namespace export checkSolverInputFile
  namespace export checkUserLevel
  namespace export checkViewfactors
  namespace export cifExit
  namespace export closeAllWindows 
  namespace export closePanel
  namespace export closeUnmodifiedWindows 
  namespace export do_elmer_exec
  namespace export editSolverInputFile
  namespace export existsTclServer
  namespace export formElmerPostCommand
  namespace export formElmerPostCommandFile
  namespace export getCurrentDirectory
  namespace export getElmerPostHeader
  namespace export getMeshRelatedFile
  namespace export getMeshRelatedFiles
  namespace export getSolverInputFile
  namespace export getSolverMeshHeaderFile
  namespace export getSolverMeshName
  namespace export getSolverMeshDirectories
  namespace export getSolverMeshDirectory
  namespace export getSolverResultDirectory
  namespace export getSolverResultFile
  namespace export getSolverResultFileEquations
  namespace export getSolverResultFiles
  namespace export helpContents
  namespace export helpIndex
  namespace export initStatusFields
  namespace export loadMesh
  namespace export loadNewModel
  namespace export monitorMesh
  namespace export openCadFile
  namespace export openFile
  namespace export openMeshDialog
  namespace export openMeshFile
  namespace export openModelFile
  namespace export prepare_ElmerPost_exec
  namespace export prepare_elmer_exec
  namespace export presetModelData
  namespace export removeCadGeometry
  namespace export removeInactiveParameterData
  namespace export rendererDisplay
  namespace export rendererMove
  namespace export rendererReset
  namespace export rendererRotate
  namespace export rendererScale
  namespace export rendererSetRotatePriorities
  namespace export resetArrayData
  namespace export resetModelFlags
  namespace export resetObjectTable
  namespace export reset_elmer_exec
  namespace export saveCheckedModelFile
  namespace export saveElmerMeshFile
  namespace export saveElmerPostMeshFile
  namespace export saveFile
  namespace export saveModelFile
  namespace export saveModelFileAs
  namespace export saveSolverInputFile
  namespace export saveThetisMeshFile
  namespace export sendTclCommands
  namespace export send_WIN32
  namespace export setModelFile
  namespace export setRotatePriority
  namespace export setWorkingDirectory
  namespace export showCaseDirectoriesMessage
  namespace export show_results
  namespace export stopFrontTask
  namespace export storeModelInDB
  namespace export storeModelInDBAs
  namespace export storeModelInDB_old
  namespace export unloadMesh
  namespace export updateModelData
  namespace export updateRadiationStatus 
}

namespace eval MeshDefine {
  namespace export applyMeshDelete
  namespace export applyMeshNew
  namespace export applyMeshRename
  namespace export callPanel
  namespace export checkF
  namespace export checkH
  namespace export checkMeshName
  namespace export checkMeshUpdate
  namespace export checkPanelData
  namespace export meshDelete
  namespace export meshDisplay
  namespace export meshGenerate
  namespace export meshNameModified
  namespace export meshNew
  namespace export meshRename
  namespace export meshSelected
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
  namespace export setAllButtonsState
  namespace export setCommandButtonsState
  namespace export setCurrentF
  namespace export setCurrentH
  namespace export setDefaultF
  namespace export setDefaultH
  namespace export setEditButtonsState
  namespace export updateMeshNamesLB
  namespace export updateNofNodesAndElements
}

namespace eval MeshSelect {
  namespace export displayMeshInfo
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
}

namespace eval Message {
  namespace export clearMessageArea
  namespace export dispCancelOkMessage
  namespace export dispCancelYesNoMessage
  namespace export dispMessage
  namespace export dispOkCancelMessage
  namespace export dispOkMessage
  namespace export dispYesNo
  namespace export dispYesNoCancelMessage
  namespace export showMessage
  namespace export verifyAction
}

namespace eval ModelInfo {
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export panelSave
}

namespace eval ModelParameter {
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
}

namespace eval ModelProperty {
  namespace export applyDirectoryValueType
  namespace export applyModelDirectoryValue
  namespace export checkCaseDirectoryValues
  namespace export checkEntry
  namespace export checkModelDirectory
  namespace export checkModelDirectoryEntry
  namespace export checkModelNames
  namespace export checkNameForFile
  namespace export findModelMeshNames
  namespace export getMeshDirectory
  namespace export initCaseDirectoryVariables
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export panelSave
  namespace export resetFieldStates
  namespace export setCurrentMeshDirectory
  namespace export setDefaultDirectoryValue
  namespace export setModelCaseName
  namespace export updateCaseDirectoryVariables
  namespace export updateGridVariable
  namespace export updateMeshNamesAndIndices
  namespace export valueTypeMenuKey2Var
  namespace export valueTypeMenuVar2Key
  namespace export verifyMeshNameChange
}

namespace eval Object {
  namespace export codeBoundaryTypeForSorting
  namespace export findIdByParentIds
  namespace export findIdByTag
  namespace export getColor
  namespace export getDisplayed
  namespace export getIds
  namespace export getMask
  namespace export getName
  namespace export getParameterId
  namespace export getParent
  namespace export getParents
  namespace export getSelected
  namespace export getSubIds
  namespace export getSubTags
  namespace export getTag
  namespace export getType
  namespace export getTypeMask
  namespace export parameterIsApplied
  namespace export setColor
  namespace export setDisplayed
  namespace export setMask
  namespace export setName
  namespace export setParent
  namespace export setParents
  namespace export setSelected
  namespace export setSubIds
  namespace export setSubTags
  namespace export setType
  namespace export sortIdsByKeyVars
  namespace export sortIdsByParentTags
}

namespace eval Panel {
  namespace export backupFieldList
  namespace export backupFields
  namespace export cancel
  namespace export checkInputParameters
  namespace export checkMaskedParameters
  namespace export checkParameterName
  namespace export checkParameters
  namespace export clearEntryButtonCommandProc
  namespace export clearPanelButtonCommandProc
  namespace export clearPanelsButtonCommandProc
  namespace export closeEntryPanels
  namespace export compressParameterIds
  namespace export compressParameters
  namespace export createDefaultParameterLine
  namespace export createOptionMenuWidget
  namespace export defaultParameterName
  namespace export editArray
  namespace export editButtonCommandProc
  namespace export editProcedure
  namespace export editTable
  namespace export getCommonLabel
  namespace export getCurrentVariables
  namespace export getGroupFields
  namespace export getGroupLabels
  namespace export getProblemMask
  namespace export initAllPanelFields
  namespace export initCurrentPanelFields
  namespace export initField
  namespace export initFields
  namespace export isActiveField
  namespace export isCheckableDataField
  namespace export isDataField
  namespace export isGroupField
  namespace export isGroupMemberField
  namespace export isInputField
  namespace export isOutputField
  namespace export isRadioButtonField
  namespace export isScreenField
  namespace export isLikeDefaultParameterName
  namespace export ok
  namespace export panelDataChanged
  namespace export panelDataDirty
  namespace export panelDataModified
  namespace export panelValuesChanged
  namespace export procedureCheckBoxCommandProc
  namespace export procedureEntryPanelCallback
  namespace export resetField
  namespace export resetFields
  namespace export restoreFields
  namespace export selectFields
  namespace export selectProblemFields
  namespace export setProcAndTableButtonStates
  namespace export tableCheckBoxCommandProc
  namespace export tableEntryPanelCallback
  namespace export uncompressParameters
  namespace export unsetIdVariables
  namespace export unsetField
  namespace export unsetFields
  namespace export verifyCancel
  namespace export verifyEntryPanels
  namespace export verifyParameter
  namespace export verifySave
}


namespace eval PanelCheck {
  namespace export checkField
  namespace export checkFieldValue
  namespace export checkFields
  namespace export checkOptionMenu
  namespace export execPanelGroupEditProcs
  namespace export groupButton
  namespace export groupCheckBox
  namespace export groupCheckButton
  namespace export groupProcedureButton
  namespace export groupRadioButton
  namespace export setNormalTangentialLabels
  namespace export singleButton
  namespace export singleCheckBox
  namespace export singleCheckButton
  namespace export singleFileButton
  namespace export singleOptionMenu
  namespace export singleProcedureButton
  namespace export singleRadioButton
  namespace export singleText
}

namespace eval PostFileSelect {
  namespace export add
  namespace export displayFileInfo
  namespace export drop
  namespace export getDefaultInfoEntry
  namespace export getInfoItem
  namespace export getInfoItemIndex
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export post
  namespace export setInfoItem
  namespace export updateFileList
  namespace export updateFileNamesLB
}

namespace eval ProcedureEntryPanel {
  variable id 0
  variable insertMarker "Y = 1.0"
  namespace export apply
  namespace export cancel
  namespace export checkData
  namespace export checkVariable
  namespace export compileSource
  namespace export create
  namespace export delete
  namespace export editSource
  namespace export lineParser
  namespace export ok
  namespace export save
  namespace export updated
}

namespace eval Process {
  namespace export addWarning
  namespace export exists
  namespace export kill
  namespace export markRemoveable
  namespace export monitor
  namespace export resume
  namespace export setPriorityLevel
  namespace export solverBrowseParser
  namespace export start
  namespace export succesful
  namespace export suspend
}

namespace eval ProcessTable {
  namespace export addEntry
  namespace export browse
  namespace export delete
  namespace export deleteEntry
  namespace export formMeshInfo
  namespace export formSolverInfo
  namespace export info_
  namespace export kill
  namespace export openPanel
  namespace export panelOk
  namespace export priority
  namespace export refresh
  namespace export resume
  namespace export suspend
  namespace export trim
  namespace export updateProcessList
}

namespace eval Processor {
  namespace export checkNofProcessors
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
}

namespace eval RendererInfo {
  namespace export openPanel
  namespace export panelOk
}

namespace eval Screen {
  namespace export setParameters
}


namespace eval SimulationParameter {
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
}

namespace eval Solver {
  namespace export applyMiniListBox
  namespace export checkData
  namespace export checkData2
  namespace export checkPanelData
  namespace export checkSolverData
  namespace export checkSolverData_old
  namespace export checkSolverData_old2
  namespace export createAndPackParamRow
  namespace export createAndPackSystemButton
  namespace export destroyCurrentListBox
  namespace export direct_proc
  namespace export doLboxWin
  namespace export findSolverIdsForEquation
  namespace export findSubsystemNames
  namespace export findSubsystemNames_old
  namespace export findSystemNames
  namespace export formSolverSystemData
  namespace export getDefaultMesh
  namespace export getMeshIndex
  namespace export initPanelRowWidgets
  namespace export initSolverSystemData
  namespace export iter_proc
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
  namespace export setComboBoxWidgetState
  namespace export setGroupStateAndValues
  namespace export setLinearization2GroupState
  namespace export setLinearizationType
  namespace export setStabilizationType
  namespace export setSubsystemIndex
  namespace export setUseNewtonGroupState
  namespace export setWidgetState
  namespace export update
  namespace export updateActiveMeshInfo
  namespace export updateFields
  namespace export updateParameter
  namespace export updateSolvers
  namespace export updateSolvingOrder
}

namespace eval SolverControl {
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
}

namespace eval SolverOrder {
  namespace export checkPanelData
  namespace export openPanel
  namespace export packWidgets
  namespace export panelApply
  namespace export panelCancel
  namespace export panelDefault
  namespace export panelOk
  namespace export panelSave
}


namespace eval StdPanelCreate {
  namespace export calcPanelSize
  namespace export calcValuesAreaSize
  namespace export constructPanelValuesArea
  namespace export createBoundaryBoxArea
  namespace export createButtons1Area
  namespace export createButtons2Area
  namespace export createField
  namespace export createObjectBoxArea
  namespace export createPanelWidgets
  namespace export createParamBoxArea
  namespace export createEquationButtonsArea
  namespace export createValuesArea
  namespace export destroyValuesArea
  namespace export openPanel
  namespace export packBoundaryBoxArea
  namespace export packButtons1Area
  namespace export packButtons2Area
  namespace export packField
  namespace export packObjectBoxArea
  namespace export packPanel
  namespace export packParamBoxArea
  namespace export packEquationButtonsArea
  namespace export packValuesArea
  namespace export setNofValuesAreaFrames
}

namespace eval StdPanelExec {
  namespace export addParameter
  namespace export attachParamId
  namespace export attachParameter
  namespace export boundaryCntrlDblSelected
  namespace export boundaryDblSelected
  namespace export boundarySelected
  namespace export boundarySelectedPre
  namespace export checkTargetSelections
  namespace export clearValuesArea
  namespace export compressParameterIds
  namespace export createParameterBackupData
  namespace export deleteActiveParameterData
  namespace export deleteParameter
  namespace export deleteParameterById
  namespace export detachInconsistentAttachments
  namespace export detachParameter
  namespace export execStandardPanel
  namespace export findAndActivateParameter
  namespace export findAndActivateTarget
  namespace export findCurrentBoundaryIds
  namespace export findObjectParent
  namespace export findSubObjectParent
  namespace export getNextParameterInstance
  namespace export getObjectMask
  namespace export getPanelArray
  namespace export getSelectedObjectMask
  namespace export getTargetMask
  namespace export markParamRowUpdated
  namespace export objectDblSelected
  namespace export objectSelected
  namespace export objectSelectedPre
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
  namespace export parameterRowWasUpdated
  namespace export parameterSelected
  namespace export parameterSelectedPre
  namespace export removeGridDefinitions
  namespace export resetDataLists
  namespace export resetUpdateLists
  namespace export selectBoundaries
  namespace export selectedObjectHasElements
  namespace export selectedTargetsAreUnattached
  namespace export setActiveEquationButtons
  namespace export setCurrentBodyIdsAndNames
  namespace export setCurrentBoundaryName
  namespace export setValuesAreaStates
  namespace export storeCurrentParameterName
  namespace export updateBoundaryLBData
  namespace export updateCurrentParameterName
  namespace export updateParameter
}

namespace eval StdPanelInit {
  namespace export constructLBLists
  namespace export constructBoundaryLBList
  namespace export constructObjectLBList
  namespace export constructParameterLBList
  namespace export constructSparseQueryList
  namespace export createEquationAndObjectMasks
  namespace export createObjectAndBoundaryLists
  namespace export createObjectTableBodyMasks
  namespace export createObjectTableMasks
  namespace export createObjectTableSubElementMasks
  namespace export createPanelData
  namespace export createParameterLists
  namespace export createWidgetBindings
  namespace export fillLBLists
  namespace export fillBoundaryListBox
  namespace export fillObjectListBox
  namespace export fillPanelData
  namespace export fillParameterListBox
  namespace export findAllBoundaryIds
  namespace export findAllObjectIds
  namespace export formObjectLBRow
  namespace export formParameterLBRow
  namespace export formTargetLBRow
  namespace export initFieldValue
  namespace export initFieldValues
  namespace export initPanelData
  namespace export setObjectTableSubIds
  namespace export updateEquationDataAndMasks
}

namespace eval SystemInfo {
  namespace export openPanel
  namespace export panelOk
}

namespace eval TableEntry {
  variable id 0
  variable sortMarker "XYZ"
  namespace export addListEntry
  namespace export copyListEntry
  namespace export create
  namespace export delete
  namespace export deleteListEntry
  namespace export getDataList
  namespace export insertListEntry
  namespace export toggleLineNumbers
  namespace export updateDataList
  namespace export updateListEntry
}

namespace eval TableEntryPanel {
  variable id 0
  variable sortMarker "XYZ"
  namespace export apply
  namespace export cancel
  namespace export checkDataList
  namespace export checkDataSize
  namespace export checkEntry
  namespace export checkVariable
  namespace export create
  namespace export delete
  namespace export ok
  namespace export save
}

namespace eval TextBrowser {
  variable id 0
  namespace export browse
  namespace export closed
  namespace export save
}

namespace eval Timestep {
  namespace export addTimestep
  namespace export addTimestepSet
  namespace export applyTimestepType
  namespace export checkPanelData
  namespace export checkSteadyStateData
  namespace export checkTimestepEntry
  namespace export checkTimesteppingMethodValue
  namespace export checkTransientData
  namespace export copyTimestepListRow
  namespace export createTimestepSetOptionMenu
  namespace export defineTimestepSize
  namespace export defineTimesteps
  namespace export deleteTimestep
  namespace export deleteTimestepSet
  namespace export initTimestepList
  namespace export insertTimestep
  namespace export isTransient
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
  namespace export selectTimestepSet
  namespace export setTimestepSet
  namespace export setTimestepStatus
  namespace export setTimesteppingMethod
  namespace export setTimesteppingMethod_BDF
  namespace export setTimesteppingMethod_Newmark
  namespace export steadyCoupledConfigure
  namespace export updateProblemMask
  namespace export updateTimestep
  namespace export updateTimestepListAndTotals
  namespace export updateTimestepListData
  namespace export updateTimestepSet
}

namespace eval UserDefined {
  namespace export readDefinitionFile
  namespace export addEquation
  namespace export addField
}

namespace eval UserSetting {
  namespace export applyUserLevel
  namespace export checkPanelData
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelOk
  namespace export panelSave
  namespace export saveForSession
  namespace export saveInModel
  namespace export saveInSettingsFile
}

namespace eval Util {
  namespace export checkAndCreateDirectory2B
  namespace export checkAndCreateDirectory3B
  namespace export checkPanelWindow
  namespace export combineMasks
  namespace export compareTimestamps
  namespace export cpp_exec
  namespace export createDirectory
  namespace export doWait
  namespace export execArrayProc
  namespace export execProc
  namespace export fillDirectoryEntry
  namespace export findNextSubPattern
  namespace export formatNumber
  namespace export getArrayValue
  namespace export getCoordinateLabel
  namespace export getSubdirectoryList
  namespace export isSimpleFilename
  namespace export isWritableFile
  namespace export makePathAbsolute
  namespace export masksAreEqual
  namespace export masksAreMatching
  namespace export normalizeVector
  namespace export patternMatchesString
  namespace export patternMatchesString_org
  namespace export restoreFlagValue
  namespace export setCurrentTimestamp
  namespace export setFlagValue
  namespace export setMainWindowTitle
  namespace export setTearOffTitle
  namespace export timeClock2Str
  namespace export timeSec2HourStr
  namespace export timeSec2Str
  namespace export timeStr2Clock
  namespace export unsetArrayVariables
  namespace export unsetArrayIdVariables
  namespace export updateGui
  namespace export updateFlagNameValues
  namespace export updateMainWindowTitle
}

namespace eval VertexDisplay {
  namespace export all
  namespace export none
  namespace export select
  namespace export openPanel
  namespace export panelApply
  namespace export panelCancel
  namespace export panelCheck
  namespace export panelOk
  namespace export setVertexDisplayMode
  namespace export setDisplayMode
}

namespace eval Widget {
  namespace export bindEditingEvent
  namespace export checkBox_generic
  namespace export configure1
  namespace export configureCheckField
  namespace export configureEntryField
  namespace export configureS
  namespace export configureSBR
  namespace export dispDataEntry 
  namespace export dispDirectoryEntry 
  namespace export entryDouble-1
  namespace export entryEnter
  namespace export entryEvent
  namespace export entryFocusIn
  namespace export entryFocusOut
  namespace export entryKeyPress
  namespace export entryKeyPress-Down
  namespace export entryKeyPress-Escape
  namespace export entryKeyPress-Return
  namespace export entryKeyPress-Up
  namespace export entryKeyRelease
  namespace export entryLeave
  namespace export isEditingEvent
  namespace export panelFocusIn
  namespace export panelFocusIn
  namespace export panelFocusOut
  namespace export panelFocusOut
  namespace export radioButton_generic
  namespace export setCheckBoxBindings
  namespace export setCheckBoxStatus
  namespace export setEntryBindings
  namespace export setEntryStatus
  namespace export setOptionMenuStatus
  namespace export setPanelBindings
  namespace export setRadioButtonBindings
  namespace export setRadioButtonStatus
  namespace export setStatus
}



# For debugging
# =============

namespace eval SourceList {
  namespace export panelOk
}

namespace eval TclCommand {
  namespace export panelOk
}




# end ecif_tcl_namespaces.tcl
# ********************
