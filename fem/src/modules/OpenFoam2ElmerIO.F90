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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 26.2.2018
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for reading data from OpenFOAM cell centers and interpolating that to
!> continuous Elmer field.
!------------------------------------------------------------------------------
SUBROUTINE OpenFoam2ElmerFit( Model,Solver,dt,TransientSimulation )
  
  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  
! local variables
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var
  TYPE(Mesh_t), POINTER :: Mesh, OFMesh
  CHARACTER(LEN=MAX_NAME_LEN) :: FileName, DirName, BaseDir, OFfile
  INTEGER :: i, NoDir, IOStatus, PassiveCoord
  LOGICAL :: Found, Visited = .FALSE., GotData
  REAL(KIND=dp) :: MinF, MaxF, MeanF, Coeff, Norm
  REAL(KIND=dp), POINTER :: RhsVector(:), WeightVector(:), OfField(:)
  TYPE(Matrix_t), POINTER :: Amat
  
  SAVE OFMesh, NoDir, Visited
  

  CALL Info('OpenFoam2ElmerFit','-----------------------------------------', Level=4 )
  CALL Info('OpenFoam2ElmerFit','Projecting OpenFOAM data to Elmer mesh',Level=4) 
    
  ! The variable containing the field contributions
  !------------------------------------------------------------------------
  Params => GetSolverParams()
  Mesh => GetMesh()

  Var => Solver % Variable
  IF(.NOT. ASSOCIATED( Var ) ) THEN
    CALL Fatal('OpenFoam2ElmerFit','Variable must exist for the solver!')
  END IF

  Amat => Solver % Matrix
  RhsVector => Amat % Rhs
  ALLOCATE( WeightVector( SIZE( RhsVector ) ) )  
  
  ! If we visit this the second time, then destroy the structures that were saved last time.
  IF( Visited ) THEN
    CALL ReleaseMesh( OFMesh )
    DEALLOCATE( OFMesh ) 
  END IF
  OFMesh => AllocateMesh()
  
  
  ! If the Elmer mesh has different dimension we may make a simple
  ! dimensional reduction for the OpenFOAM mesh.
  !-------------------------------------------------------------------------
  PassiveCoord = ListGetInteger( Solver % Values,'Passive OpenFOAM Coordinate',Found )
  IF( .NOT. Found .AND. Mesh % MeshDim < 3 ) THEN
    CALL Warn('OpenFOAM2Elmer','Dimension of Elmer mesh is reduced, and OpenFOAM not?!')
  END IF  
  
  BaseDir = GetString( Params,'OpenFOAM Directory',Found)
  IF( Found ) THEN
    CALL Info('OpenFoam2ElmerFit','Using given > OpenFOAM Directory < : '//TRIM(BaseDir),Level=6)
  ELSE
    CALL Fatal('OpenFoam2ElmerFit','> OpenFOAM Directory < must exist for the solver!')
  END IF
   
  ! If the blocks do not exist find them 
  ! When they are stored as keywords the user may give them also manually
  IF( .NOT. ListCheckPresent(Params,'OpenFOAM Mesh 1') ) THEN
    CALL OpenFOAMBlocks()
  END IF
  NoDir = NINT( ParallelReduction(1.0_dp * NoDir ) )
  CALL Info('OpenFOAM2ElmerFit','Number of active OpenFOAM blocks: '//TRIM(I2S(NoDir)),Level=5)

  
  CALL DefaultInitialize()  
  WeightVector = 0.0_dp
   
  DO i = 1, NoDir    
    IF( NoDir > 1 ) THEN
      CALL Info('OpenFOAM2ElmerFit','****************************',Level=8)
      CALL Info('OpenFOAM2ElmerFit','Treating OpenFOAM block: '//TRIM(I2S(i)),Level=5)
    END IF
    DirName = ListGetString(Params,'OpenFOAM Mesh '//TRIM(I2S(i)),Found)
    IF(.NOT. Found ) CALL Fatal('OpenFoam2ElmerFit','Could not find keyword: '//TRIM(DirName))
    
    FileName = TRIM(DirName)//'C'    
    CALL Info('OpenFoam2ElmerFit','Reading OpenFOAM cell center from file: '//TRIM(FileName),Level=5)

    CALL CreateFOAMMesh(FileName,OFMesh)
    ALLOCATE( OFField(OFMesh % NumberOfNodes)  ) 
    
    CALL ReadFOAMField(Filename, GotData)
    
    IF( GotData ) THEN
      MinF = MINVAL( OFField )
      MaxF = MAXVAL( OFField )
      MeanF = SUM( OFField ) / SIZE( OFField )    
      
      WRITE( Message,'(A,ES12.5)') 'Minimum field value: ',MinF
      CALL Info('OpenFoam2ElmerFit',Message,Level=6)
      WRITE( Message,'(A,ES12.5)') 'Maximum field value: ',MaxF
      CALL Info('OpenFoam2ElmerFit',Message,Level=6)
      WRITE( Message,'(A,ES12.5)') 'Average field value: ',MeanF
      CALL Info('OpenFoam2ElmerFit',Message,Level=6)      

      CALL DataAssembly()
    END IF
      
    ! We can only have one OpenFOAM mesh at a time, hence release the structures if we have a second mesh.
    IF( i < NoDir ) CALL ReleaseMesh( OFMesh )
    DEALLOCATE( OFField ) 
  END DO

  Coeff = ListGetCReal( Params,'Fit Coefficient',Found )
  IF( .NOT. Found ) Coeff = 1.0_dp

  !PRINT *,'rhs sum:',SUM( RhsVector )
  !PRINT *,'weight sum:',SUM( WeightVector )
  !PRINT *,'diag sum:',SUM( Amat % Values( Amat % Diag ) )

  RhsVector = Coeff * RhsVector 
  Amat % Values( Amat % Diag ) = Coeff * WeightVector 
    
  CALL DiffusionAssembly()

  CALL DefaultFinishBulkAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()
  
  Norm = DefaultSolve( )
  
  Visited = .TRUE.
    
  CALL Info('OpenFoam2ElmerFit','All done', Level=4 )
  CALL Info('OpenFoam2ElmerFit','-----------------------------------------', Level=4 )  
  
  
CONTAINS 
  
  SUBROUTINE OpenFOAMBlocks( )
    
    CHARACTER(LEN=MAX_NAME_LEN) :: DirCommand
    LOGICAL :: FileExists
    INTEGER, PARAMETER :: InFileUnit = 28
    
    
    NoDir = 0

#ifdef __INTEL_COMPILER
    ! Fortran standard states that inquiry for a file returns true if the queried entity is a file
    INQUIRE( Directory = TRIM(BaseDir), Exist = FileExists )
#else
    INQUIRE( File = TRIM(BaseDir), Exist = FileExists )
#endif
    IF(.NOT. FileExists ) THEN
      CALL Fatal('OpenFoam2ElmerFit','OpenFOAM directory does not exist: '//TRIM(BaseDir))
    END IF

    DirName = TRIM(BaseDir)//'/0/'
#ifdef __INTEL_COMPILER
    INQUIRE( Directory = DirName, Exist = FileExists )
#else
    INQUIRE( File = DirName, Exist = FileExists )
#endif
    IF(.NOT. FileExists ) THEN
      CALL Fatal('OpenFoam2ElmerFit','OpenFOAM mesh does not exist: '//TRIM(DirName))
    END IF
    
    FileName = TRIM(DirName)//'C'
    CALL Info('OpenFoam2ElmerFit','Inquire file: '//TRIM(FileName),Level=12)
    INQUIRE( File = FileName, Exist = FileExists )
   
    IF( FileExists ) THEN
      CALL Info('OpenFoam2ElmerFit','Using OpenFOAM centers in: '//TRIM(FileName),Level=10)
      CALL ListAddString( Params, 'OpenFOAM Mesh 1', DirName, .FALSE.)
      NoDir = 1
      RETURN
    END IF

    DirCommand = 'ls -d '//TRIM(DirName)//'*/ > OpenFOAMBlocks.txt' 
    CALL Info('OpenFoam2ElmerFit','Performing command: '//TRIM(DirCommand),Level=12)
    CALL SystemCommand( DirCommand )

    OPEN(InFileUnit,File='OpenFOAMBlocks.txt',IOStat=IOstatus)
    IF(IOStatus /= 0 ) THEN
      CALL Fatal('OpenFoam2ElmerFit','Could not open file: OpenFOAMBlocks.txt')
    END IF
     
    DO
      READ(InFileUnit,'(A)',IOStat = IOStatus) DirName
      IF( IOStatus /= 0 ) EXIT
      FileName = TRIM(DirName)//'C'
      CALL Info('OpenFoam2ElmerFit','Inquire file: '//TRIM(FileName),Level=12)
      INQUIRE( File = FileName, Exist = FileExists )
      IF( FileExists ) THEN
        NoDir = NoDir + 1
        CALL Info('OpenFoam2ElmerFit','Using OpenFOAM centers in: '//TRIM(FileName),Level=10)
        CALL ListAddString( Params, 'OpenFOAM Mesh '//TRIM(I2S(NoDir)), DirName, .FALSE.)
      ELSE
        CALL Info('OpenFoam2ElmerFit','No OpenFOAM center in: '//TRIM(DirName),Level=12)
      END IF
    END DO
    CLOSE(InFileUnit)
    
    IF( NoDir == 0 ) THEN
      CALL Fatal('OpenFoam2ElmerFit','No OpenFOAM mesh blocks found!')
    ELSE
      CALL Info('OpenFoam2ElmerFit','Number of OpenFOAM blocks: '//TRIM(I2S(NoDir)),Level=10)
    END IF
    
  END SUBROUTINE OpenFOAMBlocks
    


  
  !------------------------------------------------------------------------
  !> Open file in OpenFOAM format and read the cell centers from there.
  !-------------------------------------------------------------------------
  SUBROUTINE CreateFOAMMesh( Filename, Mesh ) 
    
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t), POINTER :: Mesh
   
    INTEGER :: line,i,j,k,n
    REAL(KIND=dp) :: x,y,z
    INTEGER :: NumberOfNodes, IOStatus
    INTEGER, PARAMETER :: InFileUnit = 28
    CHARACTER(LEN=:), ALLOCATABLE :: ReadStr
    
    ALLOCATE( Mesh % Nodes )
    Mesh % NumberOfBulkElements = 0
    Mesh % NumberOfBoundaryElements = 0
    
    
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::ReadStr)                   
    
    OPEN(InFileUnit,FILE = Filename, STATUS='old', IOSTAT=IOstatus)
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('OpenFoam2ElmerFit','Could not open file for reading: '//TRIM(FileName))
    END IF
    
    CALL Info('OpenFoam2ElmerFit','Reading data points from file: '//TRIM(FileName),Level=6)
    
    j = 0
    DO Line = 1, 100
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Warn('OpenFoam2ElmerFit','End of file after '//TRIM(I2S(Line))//' lines')
        EXIT
      END IF

      j =  INDEX( ReadStr,'internalField',.TRUE.) 
      IF( j > 0 ) EXIT
    END DO

    IF( j == 0 ) THEN
      CALL Warn('OpenFoam2ElmerFit','Could not find > internalField < in header!')
    ELSE
      CALL Info('OpenFoam2ElmerFit','internalField found at line: '//TRIM(I2S(Line)),Level=7)    
    END IF

    j = INDEX( ReadStr,'nonuniform',.TRUE.)
    IF( j == 0 ) THEN
      CALL Warn('OpenFoam2ElmerFit','This routine only knows how to read nonuniform lists!')
      RETURN
    END IF

    
    READ(InFileUnit,*,IOSTAT=IOStatus) NumberOfNodes    
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('OpenFoam2ElmerFit','Could not read number of nodes!')
    END IF
    CALL Info('OpenFoam2ElmerFit','Number of OpenFOAM nodes: '&
        //TRIM(I2S(NumberOfNodes)))

    i = ListGetInteger(Params,'Number of cells',Found)
    IF( i > 0 .AND. i < NumberOfNodes ) THEN
      NumberOfNodes = i
      CALL Info('OpenFoam2ElmerFit','Limiting number of OpenFOAM nodes: '&
          //TRIM(I2S(NumberOfNodes)))
    END IF

    
    n = NumberOfNodes    
    ALLOCATE( Mesh % Nodes % x(n), &          
        Mesh % Nodes % y(n), &
        Mesh % Nodes % z(n) )

    Mesh % Nodes % x(1:n) = 0.0_dp
    Mesh % Nodes % y(1:n) = 0.0_dp
    Mesh % Nodes % z(1:n) = 0.0_dp
    
    Mesh % NumberOfNodes = NumberOfNodes
         
    ! This is just empty left paranthesis
    READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
   
    DO i=1,n
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('OpenFoam2ElmerFit','Could not read coordinate line: '//TRIM(I2S(i)))
      END IF
      
      j =  INDEX( ReadStr,'(',.TRUE.) 
      IF( j == 0 ) THEN
        CALL Fatal('OpenFoam2ElmerFit',&
            'Expecting a paranthesis at the start of OpenFOAM line: '//TRIM(I2S(i)))
      END IF
      k =  INDEX( ReadStr,')',.TRUE.) 
      IF( k == 0 ) THEN
        CALL Fatal('OpenFoam2ElmerFit',&
            'Expecting a paranthesis at the end of OpenFOAM line: '//TRIM(I2S(i)))
      END IF
      
      READ( ReadStr(j+1:k-1),*,IOSTAT=IOStatus ) x,y,z
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('OpenFoam2ElmerFit','Could not read coordinate values: '//TRIM(I2S(i)))
      END IF
      Mesh % Nodes % x(i) = x
      Mesh % Nodes % y(i) = y
      Mesh % Nodes % z(i) = z
    END DO
    CLOSE( InFileUnit ) 

    CALL Info('OpenFoam2ElmerFit','Created temporal OpenFOAM mesh just for nodes',Level=8)

    !PRINT *,'range x:',MINVAL( Mesh % Nodes % x ), MAXVAL( Mesh % Nodes % x ) 
    !PRINT *,'range y:',MINVAL( Mesh % Nodes % y ), MAXVAL( Mesh % Nodes % y ) 
    !PRINT *,'range z:',MINVAL( Mesh % Nodes % z ), MAXVAL( Mesh % Nodes % z ) 
   
    IF( PassiveCoord == 1 ) THEN
      Mesh % Nodes % x = 0.0_dp
    ELSE IF( PassiveCoord == 2 ) THEN
      Mesh % Nodes % y = 0.0_dp
    ELSE IF( PassiveCoord == 3 ) THEN
      Mesh % Nodes % z = 0.0_dp
    END IF
    
  END SUBROUTINE CreateFOAMMesh



  !------------------------------------------------------------------------
  !> Open file in OpenFOAM format result file and read the data from there.
  !-------------------------------------------------------------------------
  SUBROUTINE ReadFOAMField( Filename, GotData )
    
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    LOGICAL :: GotData

    CHARACTER(LEN=MAX_NAME_LEN) :: TFileName, TSuffix   
    INTEGER :: line,i,j,k,n,nstep
    REAL(KIND=dp) :: val
    INTEGER :: NumberOfNodes, IOStatus
    INTEGER, PARAMETER :: InFileUnit = 28
    CHARACTER(LEN=:), ALLOCATABLE :: ReadStr
    LOGICAL :: IsScalar

    GotData = .FALSE.
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::ReadStr)
    
    TSuffix = ListGetString( Params,'OpenFOAM field',Found ) 
    IF( .NOT. Found ) THEN
      CALL Fatal('OpenFoam2ElmerFit','Give OpenFOAM field name!')
    END IF
    
    nstep = ListGetInteger( Params,'OpenFOAM Timestep',Found ) 
    IF(Found ) CALL Info('OpenFoam2ElmerFit','Replacing 0 with timestep '//TRIM(I2S(nstep)))

    j = INDEX( FileName,'/0/')
    k = len_trim( FileName ) 
    
    WRITE(TFileName,'(A)') FileName(1:j)//TRIM(I2S(nstep))//FileName(j+2:k-1)//TRIM(TSuffix)
    
    OPEN(InFileUnit,FILE = TFilename, STATUS='old', IOSTAT=IOstatus)
    IF( IOStatus /= 0 ) THEN
      CALL Warn('OpenFoam2ElmerFit','Could not open file for reading: '//TRIM(TFileName))
      RETURN
    END IF
    
    CALL Info('OpenFoam2ElmerFit','Reading data field from file: '//TRIM(TFileName),Level=6)

    j = 0
    DO Line = 1, 100
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Warn('OpenFoam2ElmerFit','End of file after '//TRIM(I2S(Line))//' lines')
        GOTO 10
      END IF

      j =  INDEX( ReadStr,'internalField',.TRUE.) 
      IF( j > 0 ) THEN
        k = INDEX( ReadStr,'<vector>',.TRUE.)
        IF( k > 0 ) THEN
          CALL Fatal('OpenFoam2ElmerFit','Currently implemented only for scalar fields!')         
        END IF
        EXIT
      END IF
    END DO
    
    IF( j == 0 ) THEN
      CALL Warn('OpenFoam2ElmerFit','Could not find > internalField < in header!')
      GOTO 10
    ELSE
      CALL Info('OpenFoam2ElmerFit','internalField found at line: '//TRIM(I2S(Line)),Level=7)    
    END IF

    j = INDEX( ReadStr,'nonuniform',.TRUE.)
    IF( j == 0 ) THEN
      CALL Warn('OpenFoam2ElmerFit','This routine only knows how to read nonuniform lists!')
      GOTO 10
    END IF
    
      
    READ(InFileUnit,*,IOSTAT=IOStatus) n
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('OpenFoam2ElmerFit','Could not read number of nodes!')
    END IF
    CALL Info('OpenFoam2ElmerFit','Number of OpenFOAM nodes: '&
        //TRIM(I2S(n)),Level=10)
    
    ! This is just empty left paranthesis
    READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
   
    DO i=1,OFMesh % NumberOfNodes
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('OpenFoam2ElmerFit','Could not read coordinate line: '//TRIM(I2S(i)))
      END IF

      READ( ReadStr,*,IOSTAT=IOStatus ) val
      OFField(i) = val
             
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('OpenFoam2ElmerFit','Could not read field values: '//TRIM(I2S(i)))
      END IF
    END DO

    
    CALL Info('OpenFoam2ElmerFit','Read data from OpenFOAM mesh region',Level=7)
    GotData = .TRUE.
    
    ! PRINT *,'range f:',MINVAL( OFField ), MAXVAL( OFField )

10  CLOSE( InFileUnit )
    
  END SUBROUTINE ReadFOAMField

 

  SUBROUTINE DataAssembly()
    
    REAL(KIND=dp) :: GlobalCoords(3), LocalCoords(3), val, weight, u, v, w, DetJ
    INTEGER :: i,j,t, n, ElementIndex
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: ElementNodes
    REAL(KIND=dp), POINTER :: Basis(:)
    LOGICAL :: Stat
    
    
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n) )

    ElementIndex = 0 
    DO t = 1, OFMesh % NumberOfNodes
      GlobalCoords(1) = OFMesh % Nodes % x(t)
      GlobalCoords(2) = OFMesh % Nodes % y(t)
      GlobalCoords(3) = OFMesh % Nodes % z(t)
      
      val = OFField(t)

      CALL LocateParticleInMeshOctree( ElementIndex, GlobalCoords, LocalCoords )
      
      IF( ElementIndex == 0 ) CYCLE

      Element => Mesh % Elements( ElementIndex )
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      CALL GetElementNodes(ElementNodes,Element)
      
      u = LocalCoords(1)
      v = LocalCoords(2)
      w = LocalCoords(3)
      
      stat = ElementInfo( Element, ElementNodes, U, V, W, DetJ, Basis )
      
      DO i = 1,n
        j = Var % Perm( NodeIndexes(i) )
        
        IF( j == 0 ) CYCLE
        
        weight = Basis(i)
        
        RhsVector( j ) = RhsVector( j ) + weight * val
        WeightVector( j ) = WeightVector( j ) + weight  
        
      END DO
    END DO

    DEALLOCATE( Basis )
    
  END SUBROUTINE DataAssembly


  !------------------------------------------------------------------------
  ! Assemble the matrix equation 
  !-------------------------------------------------------------------------
  SUBROUTINE DiffusionAssembly()
    
    INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
    INTEGER :: i,j,p,q,k,t,n,istat,active,dim
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    CHARACTER(LEN=MAX_NAME_LEN) :: DiffusivityName
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: Coeff, detJ, val, DiffMatrix(3,3)
    REAL(KIND=dp), POINTER :: Hwrk(:,:,:) => Null()
    REAL(KIND=dp), POINTER :: DataDiffusivity(:,:,:)
    TYPE(Matrix_t), POINTER :: StiffMatrix 
    LOGICAL :: stat, GlobalDiffuse, LocalDiffuse, Visited = .FALSE.
    TYPE(ValueList_t), POINTER :: Material
    
    
    SAVE Visited, Nodes, STIFF, FORCE, Basis, dBasisdx, DataDiffusivity

    ! Assembly the diffusion part used for regularization
    !----------------------------------------------------------
    Coeff = GetCReal( Solver % Values,'Diffusion Coefficient',GlobalDiffuse)
    DiffusivityName = GetString( Solver % Values,'Diffusivity Name',Found )
    IF(.NOT. Found ) DiffusivityName = 'Heat Conductivity'
    LocalDiffuse = .FALSE.

    active = GetNOFActive()
    StiffMatrix => Solver % Matrix

    dim = CoordinateSystemDimension()

    IF(.NOT. Visited) THEN
      Visited = .TRUE.
      N = Solver % Mesh % MaxElementNodes 
      ALLOCATE( Basis(n), dBasisdx(n, 3), FORCE(N), STIFF(N,N), &
          DataDiffusivity( 3,3,N ), STAT=istat )
      IF( istat /= 0) CALL Fatal('OpenFoam2ElmerFit','Allocation error in DiffusionAssembly!')
    END IF
    

    DO t=1,active

      Element => GetActiveElement(t)
      n = GetElementNOFNodes(Element)
      Indexes => Element % NodeIndexes
      
      CALL GetElementNodes( Nodes, Element )
      STIFF = 0.0d0
      FORCE = 0.0d0
      
      IF( .NOT. GlobalDiffuse ) THEN
        Material => GetMaterial()
        CALL ListGetRealArray( Material,DiffusivityName,Hwrk,n,Indexes,LocalDiffuse)
        IF( LocalDiffuse ) THEN
          DataDiffusivity = 0.0d0
          IF ( SIZE(Hwrk,1) == 1 ) THEN
            DO i=1,3
              DataDiffusivity( i,i,1:n ) = Hwrk( 1,1,1:n )
            END DO
          ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
            DO i=1,MIN(3,SIZE(Hwrk,1))
              DataDiffusivity(i,i,1:n) = Hwrk(i,1,1:n)
            END DO
          ELSE
            DO i=1,MIN(3,SIZE(Hwrk,1))
              DO j=1,MIN(3,SIZE(Hwrk,2))
                DataDiffusivity( i,j,1:n ) = Hwrk(i,j,1:n)
              END DO
            END DO
          END IF
        END IF
      END IF

      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )
      DO k=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis, dBasisdx )
        
        ! Compute the local conductivity tensor
        ! -------------------------------
        IF( LocalDiffuse ) THEN
          DO p=1,dim
            DO q=1,dim
              DiffMatrix(p,q) = SUM( DataDiffusivity(p,q,1:n) * Basis(1:n) )
            END DO
          END DO
        END IF

        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        IF( GlobalDiffuse ) THEN
          DO i=1,n
            DO j=1,n
              STIFF(i,j) = STIFF(i,j) + IP % s(k) * DetJ * &
                  Coeff * SUM( dBasisdx(i,1:dim) * dBasisdx(j,1:dim) ) 
            END DO
          END DO
        ELSE IF( LocalDiffuse ) THEN                        
          DO i=1,n
            DO j=1,n
              STIFF(i,j) = STIFF(i,j) + IP % s(k) * DetJ * &
                  SUM(MATMUL(DiffMatrix(1:dim,1:dim), dBasisdx(j,1:dim)) * dBasisdx(i,1:dim)) 
            END DO
          END DO
        END IF

      END DO
      
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO

  END SUBROUTINE DiffusionAssembly

!------------------------------------------------------------------------------
  

END SUBROUTINE OpenFoam2ElmerFit
