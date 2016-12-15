!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *                               EXPERIMENTAL
! *
! *  Module for exporting parallel results in Xdmf/HDF5 file format
! *
! ******************************************************************************
! *
! *  Authors: Mikko Lyly
! *
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 11 Feb 2011
! *
! *****************************************************************************/
!
!------------------------------------------------------------------------------
SUBROUTINE XdmfWriter(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE HDF5
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: MaxFields = 1000
  INTEGER :: NofScalarFields, NofVectorFields
  CHARACTER(LEN=MAX_NAME_LEN) :: ScalarFieldNames(MaxFields)
  CHARACTER(LEN=MAX_NAME_LEN) :: VectorFieldNames(MaxFields)
  CHARACTER(LEN=MAX_NAME_LEN) :: BaseFileName
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: PEs, MyPE, i, ierr
  INTEGER, ALLOCATABLE :: NofNodes(:), NofElements(:), NofStorage(:), itmp(:)
  INTEGER(HID_T) :: file_id, plist_id, fill_type_id
  LOGICAL :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: Str
  REAL(KIND=dp) :: RealTime, StartTime, TotalTime
  INTEGER:: Order203(3), Order510(10), Order820(20)
  REAL(KIND=dp), ALLOCATABLE :: TimeValues(:), dtmp(:)
  TYPE(Variable_t), POINTER :: Var
  INTEGER :: Counter = 1
  SAVE TimeValues, Counter
!------------------------------------------------------------------------------
  StartTime = RealTime()
  Mesh => GetMesh()
  PEs = ParEnv % PEs
  MyPE = ParEnv % MyPE + 1

  ! Save simulation time history:
  !-------------------------------
  IF(Counter == 1) THEN
     ALLOCATE(TimeValues(1))
     TimeValues(1) = GetTime()
  ELSE
     ALLOCATE(dtmp(SIZE(TimeValues)))
     dtmp = TimeValues
     DEALLOCATE(TimeValues)
     ALLOCATE(TimeValues(SIZE(dtmp)+1))
     TimeValues(1:SIZE(dtmp)) = dtmp
     TimeValues(SIZE(dtmp)+1) = GetTime()
     DEALLOCATE(dtmp)
  END IF

  ! Determine the base file name and field variables:
  !---------------------------------------------------
  BaseFileName = ListGetString(Solver % Values, 'base file name', Found)
  IF(.NOT.Found) BaseFileName = 'results'
  CALL INFO('XdmfWriter', 'Base file name: '//TRIM(BaseFileName))

  CALL FindScalarFields(NofScalarFields, ScalarFieldNames)
  DO i = 1, NofScalarFields
     CALL INFO('XdmfWriter', 'Scalar field: '//TRIM(ScalarFieldNames(i)))
  END DO
  
  CALL FindVectorFields(NofVectorFields, VectorFieldNames)
  DO i = 1, NofVectorFields
     CALL INFO('XdmfWriter', 'Vector field: '//TRIM(VectorFieldNames(i)))
  END DO

  ! Set up node permutation vectors for quadratic elements:
  !---------------------------------------------------------
  Order203(:) = (/ 1,3,2 /)
  Order510(:) = (/ 1,2,4,3,5,9,8,7,6,10/)
  Order820(:) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16 /)

  ! Determine Nof nodes, Nof elements and mixed xdmf storage size for all PEs:
  !----------------------------------------------------------------------------
  ALLOCATE(itmp(PEs))

  ALLOCATE(NofNodes(PEs), NofElements(PEs), NofStorage(PEs))

  NofNodes = 0; itmp = 0;  itmp(MyPE) = Mesh % NumberOfNodes
  CALL MPI_ALLREDUCE(itmp, NofNodes, PEs, MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)

  NofElements = 0; itmp = 0; itmp(MyPE) = GetNofActive()
  CALL MPI_ALLREDUCE(itmp, NofElements, PEs, MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)

  NofStorage = 0; itmp = 0; itmp(MyPE) = GetElementStorageSize()
  CALL MPI_ALLREDUCE(itmp, NofStorage, PEs, MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)

  DEALLOCATE(itmp)  

  ! Initialize and create/open the hdf5 file collectively:
  !--------------------------------------------------------
  CALL h5open_f(ierr)
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
  CALL h5pset_fapl_mpio_f(plist_id, ELMER_COMM_WORLD, MPI_INFO_NULL, ierr)

  IF(Counter == 1) THEN
     CALL h5fcreate_f(TRIM(BaseFileName)//'.h5', H5F_ACC_TRUNC_F, file_id, ierr, access_prp = plist_id)
  ELSE
     CALL h5fopen_f(TRIM(BaseFileName)//'.h5', H5F_ACC_RDWR_F, file_id, ierr, access_prp = plist_id)
  END IF

  CALL h5pclose_f(plist_id, ierr)

  ! Determine output precision for reals:
  !---------------------------------------
  fill_type_id = H5T_NATIVE_DOUBLE
  IF(ListGetLogical(Solver % Values, 'single precision', Found)) fill_type_id = H5T_NATIVE_REAL
  IF(.NOT.Found) fill_type_id = H5T_NATIVE_DOUBLE

  ! Write nodes, elements and part numbers (first time only):
  !-----------------------------------------------------------
  IF(Counter == 1) THEN
     CALL WriteNodes(file_id, PEs, MyPE, NofNodes, fill_type_id)
     CALL WriteElements(file_id, PEs, MyPE, NofStorage)
     CALL WriteParts(file_id, PEs, MyPE, NofNodes, fill_type_id)
  END IF

  ! Export field variables:
  !-------------------------
  DO i = 1, NofScalarFields
     CALL WriteScalars(file_id, PEs, MyPE, NofNodes, ScalarFieldNames(i), fill_type_id)
  END DO

  DO i = 1, NofVectorFields
     CALL WriteVectors(file_id, PEs, MyPE, NofNodes, VectorFieldNames(i), fill_type_id)
  END DO

  ! Rewrite the xdmf-file:
  !------------------------
  IF(MyPE == 1) CALL WriteXdmfFile(PEs, NofNodes, NofElements, &
       NofStorage, NofScalarFields, ScalarFieldNames, NofVectorFields, &
       VectorFieldNames, BaseFileName, fill_type_id)

  ! Finalize:
  !-----------
  CALL h5fclose_f(file_id, ierr)
  DEALLOCATE(NofElements, NofNodes, NofStorage)
  TotalTime = RealTime() - StartTime
  WRITE(Str, *) TotalTime
  CALL INFO('XdmfWriter', 'Total write time (REAL): '//TRIM(ADJUSTL(Str)))
  Counter = Counter + 1

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE FindScalarFields(NofScalarFields, ScalarFieldNames)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: NofScalarFields
    CHARACTER(LEN=MAX_NAME_LEN), INTENT(OUT) :: ScalarFieldNames(:)

    INTEGER :: id
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: LHS, RHS, Tmp1
    TYPE(Variable_t), POINTER :: Variable

    NofScalarFields = 0

    DO id = 1, SIZE(ScalarFieldNames)
       WRITE(Tmp1, *) id

       WRITE(LHS, '(A)') 'scalar field '//TRIM(ADJUSTL(Tmp1))

       RHS = ListGetString(Solver % Values, TRIM(LHS), Found)

       IF(.NOT.Found) CYCLE

       Variable => VariableGet(Solver % Mesh % Variables, TRIM(RHS))

       IF(.NOT.ASSOCIATED(Variable)) THEN
          CALL INFO('XdmfWriter', 'Bad scalar field: '//TRIM(RHS))
          CYCLE
       END IF

       NofScalarFields = NofScalarFields + 1
       ScalarFieldNames(NofScalarFields) = TRIM(RHS)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE FindScalarFields
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE FindVectorFields(NofVectorFields, VectorFieldNames)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: NofVectorFields
    CHARACTER(LEN=MAX_NAME_LEN), INTENT(OUT) :: VectorFieldNames(:)

    INTEGER :: id
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: LHS, RHS, Tmp1
    TYPE(Variable_t), POINTER :: Variable

    NofVectorFields = 0

    DO id = 1, SIZE(VectorFieldNames)
       WRITE(Tmp1, *) id

       WRITE(LHS, '(A)') 'vector field '//TRIM(ADJUSTL(Tmp1))

       RHS = ListGetString(Solver % Values, TRIM(LHS), Found)

       IF(.NOT.Found) CYCLE

       Variable => VariableGet(Solver % Mesh % Variables, TRIM(RHS))

       IF(.NOT.ASSOCIATED(Variable)) THEN
          ! Try componentwise:
          !--------------------
          WRITE(Tmp1, '(A)') TRIM(RHS)//' 1'
          Variable => VariableGet(Solver % Mesh % Variables, TRIM(Tmp1))
          IF(.NOT.ASSOCIATED(Variable)) THEN
             CALL INFO('XdmfWriter', 'Bad vector field: '//TRIM(RHS))
             CYCLE
          END IF
       END IF

       NofVectorFields = NofVectorFields + 1
       VectorFieldNames(NofVectorFields) = TRIM(RHS)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE FindVectorFields
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  INTEGER FUNCTION GetXdmfCode(ElementCode) RESULT(XdmfCode)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: ElementCode, XdmfCode
    
    SELECT CASE(ElementCode)
    CASE(202)
       XdmfCode = 2 ! polyline
    CASE(303)
       XdmfCode = 4 ! linear triangle
    CASE(404)
       XdmfCode = 5 ! linear quadrilateral
    CASE(504)
       XdmfCode = 6 ! linear tetrahedron
    CASE(510)
       XdmfCode = 38 ! quadratic tetrahedron
    CASE(808)
       XdmfCode = 9 ! linear hexahedron
    CASE(820)
       XdmfCode = 48 ! quadratic hexahedron
    CASE DEFAULT
       XdmfCode = -1 ! not supported, yet
    END SELECT

!XDMF_NOTOPOLOGY     0x0
!XDMF_POLYVERTEX     0x1
!XDMF_POLYLINE       0x2
!XDMF_POLYGON        0x3
!XDMF_TRI            0x4
!XDMF_QUAD           0x5
!XDMF_TET            0x6
!XDMF_PYRAMID        0x7
!XDMF_WEDGE          0x8
!XDMF_HEX            0x9
!XDMF_EDGE_3         0x0022
!XDMF_TRI_6          0x0024
!XDMF_QUAD_8         0x0025
!XDMF_TET_10         0x0026
!XDMF_PYRAMID_13     0x0027
!XDMF_WEDGE_15       0x0028
!XDMF_WEDGE_18       0x0029
!XDMF_HEX_20         0x0030
!XDMF_HEX_24         0x0031
!XDMF_HEX_27         0x0032
!XDMF_MIXED          0x0070
!XDMF_2DSMESH        0x0100
!XDMF_2DRECTMESH     0x0101
!XDMF_2DCORECTMESH   0x0102
!XDMF_3DSMESH        0x1100
!XDMF_3DRECTMESH     0x1101
!XDMF_3DCORECTMESH   0x1102

!------------------------------------------------------------------------------
  END FUNCTION GetXdmfCode
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
  INTEGER FUNCTION GetElementStorageSize() RESULT(StorageSize)
!------------------------------------------------------------------------------
    INTEGER :: i, StorageSize
    TYPE(Element_t), POINTER :: Element
    INTEGER :: XdmfCode

    StorageSize = 0
    DO i = 1, GetNofActive()
       Element => GetActiveElement(i)
       IF(.NOT.ASSOCIATED(Element)) CYCLE
       XdmfCode = GetXdmfCode(Element % Type % ElementCode)
       IF(XdmfCode < 0) CYCLE ! unknown: skip this element
       StorageSize = StorageSize + 1
       IF(XdmfCode == 2) StorageSize = StorageSize + 1 ! polyline
       StorageSize = StorageSize + GetElementNofNodes()
    END DO

!------------------------------------------------------------------------------
  END FUNCTION GetElementStorageSize
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteNodes(file_id, PEs, MyPE, NofNodes, fill_type_id)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: file_id, PEs, MyPE, NofNodes(:)
    INTEGER(HID_T) :: fill_type_id

    INTEGER :: i, ierr
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER(HID_T) :: dset_id(PEs), filespace, memspace, plist_id
    REAL(KIND=dp), ALLOCATABLE :: data(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Str

    ! Create the datasets collectively:
    !-----------------------------------
    DO i = 1, PEs
       dims(1) = 3
       dims(2) = NofNodes(i)

       WRITE(Str, *) i
       WRITE(Str, '(A)') 'nodes_'//TRIM(ADJUSTL(Str))

       CALL h5screate_simple_f(2, dims, filespace, ierr)
       CALL h5dcreate_f(file_id, TRIM(ADJUSTL(Str)), fill_type_id, filespace, dset_id(i), ierr)
       CALL h5sclose_f(filespace, ierr)
    END DO

    ! Write the data independently:
    !-------------------------------
    ALLOCATE(data(3, NofNodes(MyPE)))

    dims(1) = SIZE(data, 1)
    dims(2) = SIZE(data, 2)

    DO i = 1, dims(2)
       data(1, i) = Mesh % Nodes % x(i)
       data(2, i) = Mesh % Nodes % y(i)
       data(3, i) = Mesh % Nodes % z(i)
    END DO
    
    CALL h5screate_simple_f(2, dims, memspace, ierr)
    CALL h5dget_space_f(dset_id(MyPE), filespace, ierr)
    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
    
    CALL h5dwrite_f(dset_id(MyPE), H5T_NATIVE_DOUBLE, data, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Finalize:
    !-----------
    DEALLOCATE(data)
    CALL h5pclose_f(plist_id, ierr)
    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)
   
    DO i = 1, PEs
       CALL h5dclose_f(dset_id(i), ierr)
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE WriteNodes
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteElements(file_id, PEs, MyPE, NofStorage)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: file_id, PEs, MyPE, NofStorage(:)

    INTEGER :: i, j, k, ierr, XdmfCode
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER(HID_T) :: dset_id(PEs), filespace, memspace, plist_id
    INTEGER, ALLOCATABLE :: data(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Str
    TYPE(Element_t), POINTER :: Element

    ! Create the datasets collectively:
    !-----------------------------------
    DO i = 1, PEs
       dims(1) = 1
       dims(2) = NofStorage(i)

       WRITE(Str, *) i
       WRITE(Str, '(A)') 'elements_'//TRIM(ADJUSTL(Str))

       CALL h5screate_simple_f(2, dims, filespace, ierr)
       CALL h5dcreate_f(file_id, TRIM(ADJUSTL(Str)), H5T_NATIVE_INTEGER, filespace, dset_id(i), ierr)
       CALL h5sclose_f(filespace, ierr)
    END DO

    ! Write the data independently:
    !-------------------------------
    ALLOCATE(data(1, NofStorage(MyPE)))

    dims(1) = SIZE(data, 1)
    dims(2) = SIZE(data, 2)

    j = 0
    DO i = 1, GetNofActive()
       Element => GetActiveElement(i)
       IF(.NOT.ASSOCIATED(Element)) CYCLE
       XdmfCode = GetXdmfCode(Element % Type % ElementCode)
       IF(XdmfCode < 0) CYCLE ! unknown: skip this element

       j = j + 1
       data(1, j) = XdmfCode

       IF(XdmfCode == 2) THEN
          j = j + 1 ! polyline: nof nodes
          data(1, j) = GetElementNofNodes()
       END IF

       DO k = 1, GetElementNofNodes()
          j = j + 1

          ! Permuted C-style numbering
          SELECT CASE(Element % Type % ElementCode)
          CASE(203)
             data(1, j) = Element % NodeIndexes(Order203(k)) - 1
          CASE(510)
             data(1, j) = Element % NodeIndexes(Order510(k)) - 1
          CASE(820)
             data(1, j) = Element % NodeIndexes(Order820(k)) - 1
          CASE DEFAULT
             data(1, j) = Element % NodeIndexes(k) - 1
          END SELECT

       END DO
    END DO

    IF(j /= NofStorage(MyPE)) CALL Fatal('XdmfWriter', 'Bad element numbering')
    
    CALL h5screate_simple_f(2, dims, memspace, ierr)
    CALL h5dget_space_f(dset_id(MyPE), filespace, ierr)
    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    CALL h5dwrite_f(dset_id(MyPE), H5T_NATIVE_INTEGER, data, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Finalize:
    !-----------
    DEALLOCATE(data)
    CALL h5pclose_f(plist_id, ierr)
    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)
   
    DO i = 1, PEs
       CALL h5dclose_f(dset_id(i), ierr)
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE WriteElements
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteParts(file_id, PEs, MyPE, NofNodes, fill_type_id)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: file_id, PEs, MyPE, NofNodes(:)
    INTEGER(HID_T) :: fill_type_id

    INTEGER :: i, ierr
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER(HID_T) :: dset_id(PEs), filespace, memspace, plist_id
    REAL(KIND=dp), ALLOCATABLE :: data(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Str

    ! Create the datasets collectively:
    !-----------------------------------
    DO i = 1, PEs
       dims(1) = 1
       dims(2) = NofNodes(i)

       WRITE(Str, *) i
       WRITE(Str, '(A)') 'part_number_'//TRIM(ADJUSTL(Str))

       CALL h5screate_simple_f(2, dims, filespace, ierr)
       CALL h5dcreate_f(file_id, TRIM(ADJUSTL(Str)), fill_type_id, filespace, dset_id(i), ierr)
       CALL h5sclose_f(filespace, ierr)
    END DO

    ! Write the data independently:
    !-------------------------------
    ALLOCATE(data(1, NofNodes(MyPE)))

    dims(1) = SIZE(data, 1)
    dims(2) = SIZE(data, 2)

    data = MyPE
    
    CALL h5screate_simple_f(2, dims, memspace, ierr)
    CALL h5dget_space_f(dset_id(MyPE), filespace, ierr)
    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    CALL h5dwrite_f(dset_id(MyPE), H5T_NATIVE_DOUBLE, data, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Finalize:
    !-----------
    DEALLOCATE(data)
    CALL h5pclose_f(plist_id, ierr)
    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)
   
    DO i = 1, PEs
       CALL h5dclose_f(dset_id(i), ierr)
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE WriteParts
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteScalars(file_id, PEs, MyPE, NofNodes, ScalarFieldName, fill_type_id)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: file_id, PEs, MyPE, NofNodes(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: ScalarFieldName
    INTEGER(HID_T) :: fill_type_id

    INTEGER :: i, j, ierr
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER(HID_T) :: dset_id(PEs), filespace, memspace, plist_id
    REAL(KIND=dp), ALLOCATABLE :: data(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Str, Tmp1, Tmp2
    TYPE(Variable_t), POINTER :: Var

    ! Create the datasets collectively:
    !-----------------------------------
    DO i = 1, PEs
       dims(1) = 1
       dims(2) = NofNodes(i)
       
       WRITE(Tmp1, *) i
       WRITE(Tmp2, *) Counter
       
       WRITE(Str, '(A)') TRIM(ScalarFieldName)//'_'//TRIM(ADJUSTL(Tmp2))//'_'//TRIM(ADJUSTL(Tmp1))
       
       CALL h5screate_simple_f(2, dims, filespace, ierr)
       CALL h5dcreate_f(file_id, TRIM(ADJUSTL(Str)), fill_type_id, filespace, dset_id(i), ierr)
       CALL h5sclose_f(filespace, ierr)
    END DO
    
    ! Write the data independently:
    !-------------------------------
    ALLOCATE(data(1, NofNodes(MyPE)))

    dims(1) = SIZE(data, 1)
    dims(2) = SIZE(data, 2)
    
    Var => VariableGet(Solver % Mesh % Variables, TRIM(ScalarFieldName))
    IF(.NOT.ASSOCIATED(Var)) CALL INFO('XdmfWriter', 'Scalar not found')
    
    DO i = 1, dims(2)
       j = Var % Perm(i)
       data(1, i) = Var % Values(j)
    END DO
    
    CALL h5screate_simple_f(2, dims, memspace, ierr)
    CALL h5dget_space_f(dset_id(MyPE), filespace, ierr)
    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    CALL h5dwrite_f(dset_id(MyPE), H5T_NATIVE_DOUBLE, data, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    
    ! Finalize:
    !-----------
    DEALLOCATE(data)
    CALL h5pclose_f(plist_id, ierr)
    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)
    
    DO i = 1, PEs
       CALL h5dclose_f(dset_id(i), ierr)
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE WriteScalars
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteVectors(file_id, PEs, MyPE, NofNodes, VectorFieldName, fill_type_id)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: file_id, PEs, MyPE, NofNodes(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VectorFieldName
    INTEGER(HID_T) :: fill_type_id

    INTEGER :: i, j, k, ierr, dofs
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER(HID_T) :: dset_id(PEs), filespace, memspace, plist_id
    REAL(KIND=dp), ALLOCATABLE :: data(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Str, Tmp1, Tmp2
    TYPE(Variable_t), POINTER :: Var

    ! Create the datasets collectively:
    !-----------------------------------
    DO i = 1, PEs
       dims(1) = 3
       dims(2) = NofNodes(i)
       
       WRITE(Tmp1, *) i
       WRITE(Tmp2, *) Counter
       
       WRITE(Str, '(A)') TRIM(VectorFieldName)//'_'//TRIM(ADJUSTL(Tmp2))//'_'//TRIM(ADJUSTL(Tmp1))
       
       CALL h5screate_simple_f(2, dims, filespace, ierr)
       CALL h5dcreate_f(file_id, TRIM(ADJUSTL(Str)), fill_type_id, filespace, dset_id(i), ierr)
       CALL h5sclose_f(filespace, ierr)
    END DO
    
    ! Write the data independently:
    !-------------------------------
    ALLOCATE(data(3, NofNodes(MyPE)))

    dims(1) = SIZE(data, 1)
    dims(2) = SIZE(data, 2)

    data = 0.0d0
    
    Var => VariableGet(Solver % Mesh % Variables, TRIM(VectorFieldName))

    IF(ASSOCIATED(Var)) THEN
       ! Storage type: multiple DOFs per node:
       !---------------------------------------
       dofs = Var % DOFs
       DO i = 1, dims(2)
          j = Var % Perm(i)
          DO k = 1, MIN(3, dofs)
             data(k, i) = Var % Values(dofs*(j - 1) + k)
          END DO
       END DO
    ELSE
       ! Storage type: multiple scalars per node:
       !-----------------------------------------
       DO k = 1, 3
          WRITE(Tmp1, *) k
          WRITE(Tmp2, '(A)') TRIM(VectorFieldName)//' '//TRIM(ADJUSTL(Tmp1))
          Var => VariableGet(Solver % Mesh % Variables, TRIM(Tmp2))
          IF(.NOT.ASSOCIATED(Var)) CYCLE
          dofs = Var % DOFs
          IF(dofs /= 1) CYCLE
          DO i = 1, dims(2)
             j = Var % Perm(i)
             data(k, i) = Var % Values(j)
          END DO
       END DO
    END IF

    CALL h5screate_simple_f(2, dims, memspace, ierr)
    CALL h5dget_space_f(dset_id(MyPE), filespace, ierr)
    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    CALL h5dwrite_f(dset_id(MyPE), H5T_NATIVE_DOUBLE, data, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Finalize:
    !-----------
    DEALLOCATE(data)
    CALL h5pclose_f(plist_id, ierr)
    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)
    
    DO i = 1, PEs
       CALL h5dclose_f(dset_id(i), ierr)
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE WriteVectors
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  LOGICAL FUNCTION Dmp(fid, indent, str) RESULT(ok)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: fid, indent
    CHARACTER(LEN=*) :: str
    LOGICAL :: ok

    INTEGER :: i
    CHARACTER(LEN=MAX_NAME_LEN) :: line

    WRITE(line, '(A)') REPEAT(' ', indent)//TRIM(ADJUSTL(str))//CHAR(10)
    
    WRITE(fid) TRIM(line)

    ok = .TRUE. ! TODO: Return false if WRITE fails
!------------------------------------------------------------------------------
  END FUNCTION Dmp
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
  SUBROUTINE WriteXdmfFile(PEs, NofNodes, NofElements, &
       NofStorage, NofScalarFields, ScalarFieldNames, &
       NofVectorFields, VectorFieldNames, BaseFileName, &
       fill_type_id)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: PEs, NofScalarFields, NofVectorFields
    INTEGER :: NofElements(:), NofNodes(:), NofStorage(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: ScalarFieldNames(:), VectorFieldNames(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: BaseFileName
    INTEGER(HID_T) :: fill_type_id

    CHARACTER(LEN=MAX_NAME_LEN) :: Tmp1, Tmp2, Tmp3, Tmp4, Tmp5, Tmp6
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    CHARACTER(LEN=MAX_NAME_LEN) :: H5FileName
    INTEGER :: i, j, k
    LOGICAL :: ok

    ! Initialize:
    !-------------
    FileName = TRIM(ADJUSTL(BaseFileName))//'.xmf'
    H5FileName = TRIM(ADJUSTL(BaseFileName))//'.h5'

    OPEN(UNIT=10, FILE=TRIM(FileName), FORM='unformatted', ACCESS='stream', STATUS='unknown')

    ok = Dmp(10, 0, '<?xml version="1.0" ?>')
    ok = Dmp(10, 0, '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    ok = Dmp(10, 0, '')
    ok = Dmp(10, 0, '<Xdmf>')
    ok = Dmp(10, 2, '<Domain>')

    ok = Dmp(10, 4, '<Grid Name="mesh" GridType="Collection" CollectionType="Temporal">')

    DO j = 1, Counter
       WRITE(Tmp1, *) j;             WRITE(Tmp1, '(A)') TRIM(ADJUSTL(Tmp1))
       WRITE(Tmp2, *) TimeValues(j); WRITE(Tmp2, '(A)') TRIM(ADJUSTL(Tmp2))

       ok = Dmp(10, 6, '<Grid Name="mesh_'//TRIM(Tmp1)//'" GridType="Collection" CollectionType="Spatial">')
       ok = Dmp(10, 8, '<Time Value="'//TRIM(Tmp2)//'" />')

       DO i = 1, PEs
          WRITE(Tmp1, *) i;              WRITE(Tmp1, '(A)') TRIM(ADJUSTL(Tmp1))
          WRITE(Tmp2, *) NofElements(i); WRITE(Tmp2, '(A)') TRIM(ADJUSTL(Tmp2))
          WRITE(Tmp3, *) NofStorage(i);  WRITE(Tmp3, '(A)') TRIM(ADJUSTL(Tmp3))
          WRITE(Tmp4, *) NofNodes(i);    WRITE(Tmp4, '(A)') TRIM(ADJUSTL(Tmp4))
          WRITE(Tmp5, *) j;              WRITE(Tmp5, '(A)') TRIM(ADJUSTL(Tmp5))
          WRITE(Tmp6, *) TimeValues(j);  WRITE(Tmp6, '(A)') TRIM(ADJUSTL(Tmp6))
          
          ! Init part:
          !------------
          ok = Dmp(10, 8, '<Grid Name="mesh_'//TRIM(Tmp5)//'_'//TRIM(Tmp1)//'">')
          ok = Dmp(10, 10, '<Time Value="'//TRIM(Tmp6)//'" />')
          
          ! Write elements:
          !-----------------
          ok = Dmp(10, 10, '<Topology Type="Mixed" NumberOfElements="'//TRIM(Tmp2)//'">')
          ok = Dmp(10, 12, '<DataItem Format="HDF" DataType="Int" Dimensions="'//TRIM(Tmp3)//'">')
          ok = Dmp(10, 14, TRIM(H5FileName)//':/elements_'//TRIM(Tmp1))
          ok = Dmp(10, 12, '</DataItem>')
          ok = Dmp(10, 10, '</Topology>')
          
          ! Write nodes:
          !--------------
          ok = Dmp(10, 10, '<Geometry Type="XYZ">')
          IF(fill_type_id == H5T_NATIVE_DOUBLE) THEN
             ok = Dmp(10, 12, '<DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="'//TRIM(Tmp4)//' 3">')
          ELSE
             ok = Dmp(10, 12, '<DataItem Format="HDF" DataType="Float" Precision="4" Dimensions="'//TRIM(Tmp4)//' 3">')
          END IF
          ok = Dmp(10, 14, TRIM(H5FileName)//':/nodes_'//TRIM(Tmp1))
          ok = Dmp(10, 12, '</DataItem>')
          ok = Dmp(10, 10, '</Geometry>')
          
          ! Write part number:
          !--------------------
          ok = Dmp(10, 10, ' <Attribute Name="part_number" AttributeType="Scalar" Center="Node">')
          IF(fill_type_id == H5T_NATIVE_DOUBLE) THEN
             ok = Dmp(10, 12, ' <DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="'//TRIM(Tmp4)//' 1">')
          ELSE
             ok = Dmp(10, 12, ' <DataItem Format="HDF" DataType="Float" Precision="4" Dimensions="'//TRIM(Tmp4)//' 1">')
          END IF
          ok = Dmp(10, 14, TRIM(H5FileName)//':/part_number_'//TRIM(Tmp1))
          ok = Dmp(10, 12, '</DataItem>')
          ok = Dmp(10, 10, '</Attribute>')
          
          ! Write scalar fields:
          !----------------------
          DO k = 1, NofScalarFields
             ok = Dmp(10, 10, ' <Attribute Name="'//TRIM(ScalarFieldNames(k))//'" AttributeType="Scalar" Center="Node">')
             IF(fill_type_id == H5T_NATIVE_DOUBLE) THEN
                ok = Dmp(10, 12, ' <DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="'//TRIM(Tmp4)//' 1">')
             ELSE
                ok = Dmp(10, 12, ' <DataItem Format="HDF" DataType="Float" Precision="4" Dimensions="'//TRIM(Tmp4)//' 1">')
             END IF
             ok = Dmp(10, 14, TRIM(H5FileName)//':/'//TRIM(ScalarFieldNames(k))//'_'//TRIM(Tmp5)//'_'//TRIM(Tmp1))
             ok = Dmp(10, 12, '</DataItem>')
             ok = Dmp(10, 10, '</Attribute>')
          END DO
          
          ! Write vector fields:
          !----------------------
          DO k = 1, NofVectorFields
             ok = Dmp(10, 10, ' <Attribute Name="'//TRIM(VectorFieldNames(k))//'" AttributeType="Vector" Center="Node">')
             IF(fill_type_id == H5T_NATIVE_DOUBLE) THEN
                ok = Dmp(10, 12, ' <DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="'//TRIM(Tmp4)//' 3">')
             ELSE
                ok = Dmp(10, 12, ' <DataItem Format="HDF" DataType="Float" Precision="4" Dimensions="'//TRIM(Tmp4)//' 3">')
             END IF
             ok = Dmp(10, 14, TRIM(H5FileName)//':/'//TRIM(VectorFieldNames(k))//'_'//TRIM(Tmp5)//'_'//TRIM(Tmp1))
             ok = Dmp(10, 12, '</DataItem>')
             ok = Dmp(10, 10, '</Attribute>')
          END DO
          
          ! Finalize part:
          !----------------
          ok = Dmp(10, 8, '</Grid>') ! part
       END DO

       ok = Dmp(10, 6, '</Grid>') ! spatial collection

    END DO

    ok = Dmp(10, 4, '</Grid>') ! temporal collection
    ok = Dmp(10, 2, '</Domain>')
    ok = Dmp(10, 0, '</Xdmf>')

    CLOSE(10)

!------------------------------------------------------------------------------
  END SUBROUTINE WriteXdmfFile
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE XdmfWriter
!------------------------------------------------------------------------------
