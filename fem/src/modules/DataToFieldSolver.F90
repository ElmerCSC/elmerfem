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
! *  Authors: Peter Råback, Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.06.2011
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Solves by Galerkin method the problem where a continuous field is 
!> fitted to data. The data may be created in advance or it may be 
!> given as a property of discrete particles. The data has a contribution
!> on the r.h.s. of the equation only. Regularization (i.e. diffusion) may be
!> added to reduce noise from the fitting of the data. Also data may be used
!> only selectively as defined by some mask variable. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE DataToFieldSolver( Model,Solver,dt,TransientSimulation )
  
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
  REAL(KIND=dp), POINTER :: WeightVector(:),ForceVector(:),MaskVector(:),FieldVector(:)
  INTEGER, POINTER :: WeightPerm(:),ForcePerm(:),MaskPerm(:),FieldPerm(:)
  REAL(KIND=dp) :: Norm, MinMaskVal, MaxMaskVal, GlobalWeight
  LOGICAL :: GivenNormalize, NodalNormalize, Found, Found2, Mask, MaskDiffusion, &
      ConstantWeightSum, UseLogScale, RevertLogScale, ContinueWithBC, GlobalNormalize
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, FileName


  CALL Info('DataToFieldSolver','-----------------------------------------', Level=4 )
  CALL Info('DataToFieldSolver','Resolving field from given data',Level=4) 

  
  Var => Solver % Variable
  FieldVector => Var % Values
  FieldPerm => Var % Perm 
  CALL Info('DataToFieldSolver','Fitting to variable: '//TRIM(Var % Name),Level=6)

  ! The variable containing the field contributions
  !------------------------------------------------------------------------
  Params => GetSolverParams()
  VarName = GetString( Params,'Data Field Name',Found)
  IF(.NOT. Found ) VarName = GetString( Params,'Target Variable',Found)
  IF(.NOT. Found ) THEN
    CALL Fatal('DataToFieldSolver','> Target Variable < must exist for the solver!')
  END IF

  Var => VariableGet(Solver % Mesh % Variables, VarName )
  IF( ASSOCIATED( Var ) ) THEN    
    ForceVector => Var % Values
    ForcePerm => Var % Perm
  ELSE
    CALL Fatal('DataToFieldSolver','Variable not present:'//TRIM(VarName))      
  END IF

  ! If normalization is requested then need the vector of weights as well
  !------------------------------------------------------------------------
  GlobalWeight = GetCReal( Params,'Weight Coefficient',GlobalNormalize )

  GivenNormalize = GetLogical( Params,'Normalize by Given Weight',Found)
  IF(.NOT. GivenNormalize ) THEN
    GivenNormalize = GetLogical( Params,'Normalize Data by Weight',Found)
  END IF
  IF( GivenNormalize ) THEN
    VarName = GetString( Params,'Weight Field Name',Found)
    IF(.NOT. Found) VarName = GetString( Params,'Weight Variable',Found)
    IF(.NOT. Found ) THEN
      CALL Fatal('DataToFieldSolver','> Weight Variable < must exist for the solver!')
    END IF    
    Var => VariableGet(Solver % Mesh % Variables, VarName )
    IF( ASSOCIATED( Var ) ) THEN    
      WeightVector => Var % Values
      WeightPerm => Var % Perm
    ELSE
      CALL Fatal('DataToFieldSolver','Variable not present: '//TRIM(VarName))      
    END IF
    CALL Info('DataToFieldSolver','Normalizing source using: '//TRIM(VarName),Level=6) 

    ConstantWeightSum = GetLogical( Params,'Set Constant Weight Sum',Found)
  END IF
  
  NodalNormalize = GetLogical( Params,'Normalize by Nodal Weight',Found)    
  IF( NodalNormalize ) THEN
    CALL Info('DataToFieldSolver','Normalizing source using nodal weight',Level=6) 
  END IF

  IF( GivenNormalize .AND. NodalNormalize ) THEN
    CALL Fatal('DataToFieldSolver','Normalization cannot be both given and nodal!')
  END IF


  ! The variable containing the field contributions
  !------------------------------------------------------------------------
  MaskDiffusion = .FALSE.
  VarName = ListGetString( Params,'Mask Field Name', Mask )
  IF(.NOT. Mask) VarName = ListGetString( Params,'Mask Variable',Mask )
  IF( Mask ) THEN
    Var => VariableGet(Solver % Mesh % Variables, VarName )
    IF( ASSOCIATED( Var ) ) THEN    
      MaskVector => Var % Values
      MaskPerm => Var % Perm
    ELSE
      CALL Fatal('DataToFieldSolver','Variable not present: '//TRIM(VarName))      
    END IF
    CALL Info('DataToFieldSolver','Masking source using: '//TRIM(VarName),Level=6) 

    MaxMaskVal = ListGetCReal( Params,'Max Mask Value',Found ) 
    IF(.NOT. Found) MaxMaskVal = HUGE( MaxMaskVal )
    MinMaskVal = ListGetCReal( Params,'Min Mask Value',Found2 )
    IF(.NOT. Found2 ) THEN
      IF(.NOT. Found ) THEN
        MinMaskVal = 0.0_dp
      ELSE
        MinMaskVal = -HUGE(MinMaskVal) 
      END IF
    END IF      
    MaskDiffusion = GetLogical( Params,'Mask Diffusion',Found )
    IF( MaskDiffusion ) THEN
      CALL Info('DataToFieldSolver','Masking diffusion terms',Level=6) 
    END IF
  END IF


  ! If the data is retrieved from ascii table then map that into mesh
  ! Only the r.h.s. is assemblied here. 
  !------------------------------------------------------------------------
  Filename = GetString( Params,'Point Data Filename',Found)
  IF( Found ) THEN
    CALL AsciiPointsToMesh()
  END IF

  UseLogScale = GetLogical( Params,'Logarithmic Fitting',Found)
  IF( UseLogScale ) THEN
    CALL Info('DataToFieldSolver','Using logarithmic scale when fitting data!',Level=6)
  END IF

  ! Make the continue BCs using BC assembly / or the modification of bulk assembly
  ContinueWithBC = GetLogical( Params,'Continue BC With Boundary Assembly',Found )
  IF( ContinueWithBC ) THEN
    CALL Info('DataToFieldSolver','Using boundary assembly to continue the solution at BCs',Level=6)
  END IF
  

  ! Create the matrix equation with r.h.s. data and regularization
  !------------------------------------------------------------------------
  CALL DefaultInitialize()  

  CALL BulkAssembly()

  CALL DefaultFinishBulkAssembly()

  IF( ContinueWithBC ) THEN
    CALL BoundaryAssembly()
  END IF

  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()

  ! Solver the matrix 
  !------------------------------------------------------------------------

  ! If we compute in log scale and revert back then compute the change only 
  ! after reverting back to real scale.
  RevertLogScale = GetLogical( Params,'Revert Logarithmic Fitting',Found)
  IF( RevertLogScale ) THEN
    CALL Info('DataToFieldSolver','Reverting from logarithmic scale!',Level=6)
    CALL ListAddLogical( Solver % Values,'Skip Compute Nonlinear Change',.TRUE.)
    Norm = DefaultSolve()
    FieldVector = EXP( FieldVector ) 
    CALL ComputeChange( Solver, .FALSE. )
  ELSE
    Norm = DefaultSolve( )
  END IF

  CALL Info('DataToFieldSolver','All done', Level=4 )
  CALL Info('DataToFieldSolver','-----------------------------------------', Level=4 )
  
  
CONTAINS 
  

  !------------------------------------------------------------------------
  !> Go through list of given points and add their contribution to FE mesh.
  !-------------------------------------------------------------------------
  SUBROUTINE AsciiPointsToMesh()
    
    INTEGER :: i,j,n,dim,No,ElementIndex=0
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: SqrtElementMetric, Weight,u,v,w,LocalCoords(3),val,InputData(4)
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    LOGICAL :: AllocationsDone = .FALSE., Stat
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: GlobalCoords(3)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER :: Success, DataColumn, IOUnit


    SAVE :: AllocationsDone, ElementIndex, Basis, dBasisdx, ElementNodes
    
    Mesh => Solver % Mesh

    n = Mesh % MaxElementNodes
    ALLOCATE(ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
        Basis(n), dBasisdx(n,3) )

    
    GlobalCoords(3) = 0.0_dp
    dim = CoordinateSystemDimension()
    
    DataColumn = ListGetInteger(Params,'Point Data Column',Found )
    IF( .NOT. Found ) DataColumn = dim + 1

    OPEN(NEWUNIT=IOUnit,FILE = Filename, IOSTAT=Success)
    IF( Success /= 0 ) THEN
      CALL Fatal('DataToFieldSolver','Could not open file for reading: '//TRIM(FileName))
    ELSE
      CALL Info('DataToFieldSolver','Reading data points from file: '//TRIM(FileName),Level=6)
    END IF

    No = 0 
    DO WHILE(.TRUE.) 
      READ(IOUnit,*,IOSTAT=Success) InputData(1:DataColumn)
      IF( Success /= 0 ) THEN
        CALL Info('DataToFieldSolver','End of file after '//TRIM(I2S(No))//' values')
        IF( No == 0 ) CALL Fatal('DataToFieldSolver','Could not read any points from file!')
        EXIT
      END IF

      No = No + 1

      GlobalCoords(1:dim) = InputData(1:dim)
      val = InputData(DataColumn)

      CALL LocateParticleInMeshOctree( ElementIndex, GlobalCoords, LocalCoords )

      IF( ElementIndex == 0 ) CYCLE

      CurrentElement => Mesh % Elements( ElementIndex )
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      CALL GetElementNodes(ElementNodes,CurrentElement)

      u = LocalCoords(1)
      v = LocalCoords(2)
      w = LocalCoords(3)

      stat = ElementInfo( CurrentElement, ElementNodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx )

      DO i = 1,n
        j = FieldPerm( NodeIndexes(i) )

        IF( j == 0 ) CYCLE

        ! As the weight should be proportional to the particle amount rather than
        ! element volume the weight is not multiplied with local element size!
        ! Note that the weight could be also ~1/r^2 from the nodes etc.
        !-------------------------------------------------------------------------
        weight = Basis(i)
        
        ForceVector( j ) = ForceVector( j ) + weight * val
        IF( GivenNormalize ) THEN
          WeightVector( j ) = WeightVector( j ) + weight  
        END IF

      END DO      
    END DO

    CLOSE(IOUnit) 

    CALL Info('DataToFieldSolver','Done reading data points',Level=12)

  END SUBROUTINE AsciiPointsToMesh
   


  !------------------------------------------------------------------------
  ! Assemble the matrix equation 
  !-------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
    
    INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
    INTEGER :: i,j,p,q,j2,j3,k,t,n,istat,active,BoundaryNodes,dim,MaskActive
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryName, DiffusivityName
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), NodalWeight(:)
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: Coeff, detJ, WeightCorr, val, Wrhs, Wmat, DiffMatrix(3,3)
    REAL(KIND=dp), POINTER :: MatValues(:), MatRhs(:)
    REAL(KIND=dp), POINTER :: Hwrk(:,:,:) => Null()
    REAL(KIND=dp), POINTER :: DataDiffusivity(:,:,:)
    INTEGER, POINTER :: MatDiag(:)
    TYPE(Matrix_t), POINTER :: StiffMatrix 
    LOGICAL :: stat, GlobalDiffuse, LocalDiffuse, Visited = .FALSE.
    TYPE(ValueList_t), POINTER :: Material
    
    
    SAVE Visited, Nodes, STIFF, FORCE, Basis, dBasisdx, NodalWeight, &
        BoundaryPerm, BoundaryNodes, DataDiffusivity

    ! Assembly the diffusion part used for regularization
    !----------------------------------------------------------
    Coeff = GetCReal( Solver % Values,'Diffusion Coefficient',GlobalDiffuse)
    LocalDiffuse = .FALSE.

    DiffusivityName = GetString( Solver % Values,'Diffusivity Name',Found )
    IF(.NOT. Found ) DiffusivityName = 'Data Diffusivity'
    LocalDiffuse = .FALSE.

    active = GetNOFActive()
    StiffMatrix => Solver % Matrix
    MatValues => StiffMatrix % Values
    MatRhs => StiffMatrix % Rhs
    MatDiag => StiffMatrix % Diag
    IF(.NOT. ASSOCIATED( StiffMatrix ) ) THEN
      CALL Fatal('DataToFieldSolver','StiffMatrix not associated!')
    END IF

    dim = CoordinateSystemDimension()

    IF(.NOT. Visited) THEN
      Visited = .TRUE.
      N = Solver % Mesh % MaxElementNodes 
      ALLOCATE( Basis(n), dBasisdx(n, 3), FORCE(N), STIFF(N,N), &
          DataDiffusivity( 3,3,N ), STAT=istat )
      IF( istat /= 0) CALL Fatal('DataToFieldSolver','Allocation error 1 in BulkAssembly!')
      
      n = StiffMatrix % NumberOfRows
      ALLOCATE( NodalWeight( n ), STAT=istat )
      IF( istat /= 0) CALL Fatal('DataToFieldSolver','Allocation error 2 in BulkAssembly!')

      IF( .NOT. ContinueWithBC ) THEN
        N = Solver % Mesh % NumberOfNodes
        ALLOCATE( BoundaryPerm(n) )
        BoundaryPerm = 0
        BoundaryNodes = 0
        BoundaryName = ComponentName( Solver % Variable )
        BoundaryName = TRIM(BoundaryName)//' continue'
        CALL MakePermUsingMask( CurrentModel,Solver,Solver % Mesh,BoundaryName, &
            .FALSE., BoundaryPerm, BoundaryNodes )
        IF( BoundaryNodes > 0 ) THEN
          WRITE( Message,'(A,I0)') 'Nodes with BC > '// TRIM( BoundaryName ) //' <  true: ',BoundaryNodes
          CALL Info('DataToFieldSolver',Message,Level=6)
        END IF
      END IF

    END IF
    
    NodalWeight = 0.0_dp

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
        
        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        DO i=1,n

          IF( .NOT. ContinueWithBC .AND. BoundaryNodes > 0 ) THEN
            IF( BoundaryPerm( Indexes(i) ) > 0 ) CYCLE
          END IF
          
          ! Compute the rowsum that is used in the normalization
          !-----------------------------------------------------
          val = IP % s(k) * DetJ * Basis(i)        
          j = Indexes(i)
          NodalWeight( j ) = NodalWeight( j ) + val


          ! Compute the local conductivity tensor
          ! -------------------------------
          IF( LocalDiffuse ) THEN
            DO p=1,dim
              DO q=1,dim
                DiffMatrix(p,q) = SUM( DataDiffusivity(p,q,1:n) * Basis(1:n) )
              END DO
            END DO            
          END IF

          ! This condition should remove the diffusion for proper data, and 
          ! use it only for outlier data.
          !----------------------------------------------------------------
          IF( MaskDiffusion ) THEN
            j2 = Indexes(i)
            IF( ASSOCIATED( MaskPerm ) ) j2 = MaskPerm(j2)
            val = MaskVector( j2 )
            IF( .NOT. (val < MinMaskVal .OR. val > MaxMaskVal ) ) CYCLE
          END IF
          
          ! This condition removes the natural boundary condition that would 
          ! try to fix the normal gradient of the field to zero.
          ! Does not seem to work though...
          !--------------------------------------------------------------------

          IF( GlobalDiffuse ) THEN
            DO j=1,n
              STIFF(i,j) = STIFF(i,j) + IP % s(k) * DetJ * &
                  Coeff * SUM( dBasisdx(i,1:dim) * dBasisdx(j,1:dim) ) 
            END DO
          ELSE IF( LocalDiffuse ) THEN            
            DO j=1,n
              STIFF(i,j) = STIFF(i,j) + IP % s(k) * DetJ * &
                  SUM(MATMUL(DiffMatrix(1:dim,1:dim), dBasisdx(j,1:dim)) * dBasisdx(i,1:dim)) 
            END DO
          END IF

        END DO
      END DO
      
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO

    MaskActive = 0
    IF( GivenNormalize .AND. ConstantWeightSum ) THEN
      !-----------------------------------------------------------------------
      ! Set the weight to the diagonal i.e. make the mass matrix contribution  
      ! The data is normalized so that if it would be constant it would yield the 
      ! same equation as the normal one and the weights would also be constant.
      ! This way diffusion will not depend on the amount of data, whether
      ! that is desirable, or not, I don't know. 
      !-----------------------------------------------------------------------
      WeightCorr = SUM( NodalWeight) / SUM( WeightVector )
    ELSE
      WeightCorr = 1.0_dp
    END IF

    DO i=1,Solver % Mesh % NumberOfNodes

      ! If some values are not to be trusted don't include them in the allocation process
      !----------------------------------------------------------------------------------
      IF( Mask ) THEN
        j = i
        IF( ASSOCIATED( MaskPerm ) ) j = MaskPerm(i)
        IF( j > 0 ) THEN
          val = MaskVector(j)
          IF( val < MinMaskVal .OR. val > MaxMaskVal ) THEN
            MaskActive = MaskActive + 1
            CYCLE
          END IF
        END IF
      END IF
      
      ! field to be solved for
      j = FieldPerm(i)
      IF( j == 0 ) CYCLE

      ! force vector on the r.h.s.
      j2 = i
      IF( ASSOCIATED( ForcePerm ) ) j2 = ForcePerm(j2)     
      IF( j2 > 0 ) THEN
        val = ForceVector(j2)
        IF( UseLogScale ) val = LOG( val )
      ELSE
        val = 0.0_dp
      END IF

      IF( GivenNormalize ) THEN
        j3 = i
        IF( ASSOCIATED( WeightPerm ) ) j3 = WeightPerm(j3)        
        IF( j3 > 0 ) THEN
          Wmat = WeightCorr * WeightVector(j3)
        ELSE
          Wmat = 0.0_dp
        END IF
        Wrhs = WeightCorr * val
        
      ELSE IF( NodalNormalize ) THEN
        Wmat = NodalWeight(i)        
        Wrhs = val

      ELSE
        Wmat = NodalWeight(i)
        Wrhs = NodalWeight(i) * val
      END IF

      IF( GlobalNormalize ) THEN
        Wmat = GlobalWeight * Wmat
        Wrhs = GlobalWeight * Wrhs
      END IF

      k = MatDiag(j) 
      MatValues( k ) = MatValues( k ) + Wmat
      MatRhs(j) = MatRhs(j) + Wrhs
    END DO


    IF( Mask ) THEN
      WRITE( Message,'(A,I0,A,I0,A)') 'Mask is active for ',MaskActive,&
          ' nodes (out of ',SIZE(MaskVector),')'
      CALL Info('DataToFieldSolver',Message,Level=4)
    END IF

  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
  SUBROUTINE BoundaryAssembly()
!------------------------------------------------------------------------------
    INTEGER :: t,n,np
    TYPE(Element_t), POINTER :: Element,ParentElement
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: CondName

    CondName = ComponentName( Solver % Variable )
    CondName = TRIM( CondName ) //' continue'

    DO t=1,GetNOFBoundaryElements()
      
      Element => GetBoundaryElement(t)
      IF(.NOT. ActiveBoundaryElement(Element)) CYCLE
      
      BC => GetBC( Element ) 
      IF(.NOT. GetLogical( BC, CondName ,Found)) CYCLE

      n  = GetElementNOFNodes(Element)
      
      ParentElement => Element % BoundaryInfo % Left
      IF( .NOT. ASSOCIATED( ParentElement ) ) THEN
        CALL Fatal('DataToFieldSolver','Could not find parent element!')
      ELSE IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
        CALL Fatal('DataToFieldSolver','This does not make sense for internal BC!')
      END IF

      np  = GetElementNOFNodes(ParentElement)

      CALL BoundaryLocalMatrix( Element, ParentElement, n, np )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BoundaryLocalMatrix( Element, ParentElement, n, np )
!------------------------------------------------------------------------------
    INTEGER :: n, np
    TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(np,np), FORCE(np)
    REAL(KIND=dp), POINTER :: A(:,:),M(:,:)    
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), s, DetJ,u,v,w
    REAL(KIND=dp) :: ParentBasis(np),ParentdBasisdx(np,3),DiffMatrix(3,3)
    REAL(KIND=dp) :: Nrm(3),Coeff
    REAL(KIND=dp) :: DataDiffusivity(3,3,np)
    REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
    LOGICAL :: Stat,Found
    INTEGER :: i,j,p,q,t,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes, ParentNodes
    LOGICAL :: LocalDiffuse, GlobalDiffuse
    TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()    
    
    Coeff = GetCReal( Solver % Values,'Diffusion Coefficient',GlobalDiffuse)
    LocalDiffuse = .FALSE.

    CALL GetElementNodes( Nodes,Element )
    CALL GetElementNodes( ParentNodes,ParentElement )

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    
    IF( .NOT. GlobalDiffuse ) THEN
      Material => GetMaterial( ParentElement )
      CALL ListGetRealArray( Material,'Data Diffusivity',Hwrk,np,&
          ParentElement % NodeIndexes,LocalDiffuse)
      IF( LocalDiffuse ) THEN
        DataDiffusivity = 0.0d0
        IF ( SIZE(Hwrk,1) == 1 ) THEN
          DO i=1,3
            DataDiffusivity( i,i,1:np ) = Hwrk( 1,1,1:np )
          END DO
        ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Hwrk,1))
            DataDiffusivity(i,i,1:np) = Hwrk(i,1,1:np)
          END DO
        ELSE
          DO i=1,MIN(3,SIZE(Hwrk,1))
            DO j=1,MIN(3,SIZE(Hwrk,2))
              DataDiffusivity( i,j,1:np ) = Hwrk(i,j,1:np)
            END DO
          END DO
        END IF
      END IF
    END IF

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      u =  IP % U(t)
      v =  IP % V(t)
      w =  IP % W(t)

      stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )
      
      s = IP % s(t) * DetJ
      Nrm = NormalVector(Element,Nodes,u,v, .TRUE.)

      CALL GetParentUVW(Element, n, ParentElement, np, U, V, W, Basis)

      stat = ElementInfo(ParentElement, ParentNodes, U, V, W, detJ, &
          ParentBasis, ParentdBasisdx)


      ! Compute the local conductivity tensor
      ! -------------------------------
      IF( LocalDiffuse ) THEN
        DO p=1,dim
          DO q=1,dim
            DiffMatrix(p,q) = SUM( DataDiffusivity(p,q,1:np) * ParentBasis(1:np) )
          END DO
        END DO
      END IF
      

      DO p=1,np
        DO q=1,np
          IF( GlobalDiffuse ) THEN
            STIFF(p,q) = STIFF(p,q) - s * ParentBasis(p)* &
                Coeff * SUM( ParentdBasisDx(q,1:dim) * Nrm(1:dim) )
          ELSE IF( LocalDiffuse ) THEN
            STIFF(p,q) = STIFF(p,q) - s * ParentBasis(p) * &
                SUM( MATMUL(DiffMatrix, dBasisdx(q,:)) * Nrm(1:dim) )
          END IF
        END DO
      END DO
    END DO

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=ParentElement )

!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryLocalMatrix
!------------------------------------------------------------------------------

END SUBROUTINE DataToFieldSolver



!-------------------------------------------------------------------
!> Default initialization for the primary solver.
!-------------------------------------------------------------------
SUBROUTINE DataToFieldSolver_init( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE Lists

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  
! local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: GivenNormalize, NodalNormalize, Found, HaveFile
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  

  Params => GetSolverParams()

  HaveFile = ListCheckPresent( Params,'Point Data Filename')

  GivenNormalize = GetLogical( Params,'Normalize by Given Weight',Found)
  IF(.NOT. GivenNormalize) GivenNormalize = GetLogical( Params,'Normalize Data by Weight',Found)

  ! If the file is given then the two following fields will be created internally
  ! and one can allocate them here, if not given otherwise. 
  !-----------------------------------------------------------------------------
  VarName = GetString( Params,'Target Variable',Found)
  IF(.NOT. Found .AND. HaveFile ) THEN
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable',Params),'Point Data')
    CALL ListAddString( Params,'Target Variable','Point Data')
    CALL Info('DataToFieldSolver_init','Creating > Point Data < as exported variable')
  END IF

  IF( GivenNormalize .AND. HaveFile ) THEN
    VarName = GetString( Params,'Weight Variable',Found)
    IF(.NOT. Found ) THEN
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable',Params),'Point Weight')
      CALL ListAddString( Params,'Weight Variable','Point Weight')
      CALL Info('DataToFieldSolver_init','Creating > Point Weight < as exported variable')
    END IF
  END IF

  ! If the field is not given, create it
  ! A different default name if point data set is used.
  !-----------------------------------------------------
  VarName = GetString( Params,'Variable',Found)
  IF(.NOT. Found ) THEN
    IF( HaveFile ) THEN
      CALL ListAddString( Params,'Variable','Point Fit')     
    ELSE
      CALL ListAddString( Params,'Variable','Fit')
    END IF
  END IF


  CALL ListAddInteger( Params, 'Time derivative order', 0 )
  
  ! Add some cheap linear system defaults: bicgstab + none (diagonal when scaled)
  !------------------------------------------------------------------------
  CALL ListAddNewString(Params,'Linear System Solver','Iterative')
  CALL ListAddNewString(Params,'Linear System Iterative Method','bicgstab')
  CALL ListAddNewString(Params,'Linear System Preconditioning','none')
  CALL ListAddNewInteger(Params,'Linear System Max Iterations',1000)
  CALL ListAddNewInteger(Params,'Linear System Residual Output',20)
  CALL ListAddNewConstReal(Params,'Linear System Convergence Tolerance',1.0e-10_dp)
  
END SUBROUTINE DataToFieldSolver_init
