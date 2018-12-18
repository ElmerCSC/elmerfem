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
! *  Original Date: 5.9.2013
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
!> Subroutine that helps to give current sources for closed coils. Normally there
!> is the difficulty that the potential in such a coil should be discontinuous. 
!> In this solver two different potential fields are solved (say a left and a right 
!> field). The two fields have different boundary conditions that are automatically 
!> set on the bulk nodes so that a current is induced. The solution settles typically
!> on the other side such that the union of the solutions is always reasonable.
!> The main assumption is that the coil axis is aligned with the z-axis. 
!
!> The module includes a special user function that returns the potential 
!> at the good side. 
!
!> \ingroup Solvers
!------------------------------------------------------------------------------


SUBROUTINE CoilSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params 
  INTEGER :: dim
  LOGICAL :: Found, CalcCurr

  dim = CoordinateSystemDimension()
  Params => GetSolverParams()

  IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable','-nooutput CoilTmp')
  END IF

  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
      'CoilPot')

  IF( GetLogical( Params,'Coil Conductivity Fix', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
        'CoilFix')
  END IF
    
  IF( GetLogical( Params,'Coil Closed', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
        'CoilPotB')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
        'PotSelect')
   
    IF( GetLogical( Params,'Save Coil Set', Found ) ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
          'CoilSetB')
    END IF
  END IF

  IF( GetLogical( Params,'Save Coil Set', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
        'CoilSet')
  END IF
  
  IF( GetLogical( Params,'Save Coil Index', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
        'CoilIndex')
  END IF

  CalcCurr = GetLogical( Params,'Calculate Coil Current',Found )
  IF( .NOT. Found ) CalcCurr = .TRUE.
  IF( CalcCurr ) THEN
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable',Params),&
        'CoilCurrent[CoilCurrent:'//TRIM(I2S(dim))//']')
  END IF
    
! Loads are needed to compute the induced currents in a numerically optimal way
  CALL ListAddLogical( Params,'Calculate Loads',.TRUE.)

END SUBROUTINE CoilSolver_init



!------------------------------------------------------------------------------
SUBROUTINE CoilSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm, s
  INTEGER :: i,j,k,n, nb, nd, t, active, iter, Part, sgn, nsize, CoilParts, &
       MaxNonlinIter, dim,dimi,ierr, NoCoils, MaxNoCoils
  INTEGER, POINTER :: Perm(:), Set(:),TargetBodies(:)
  INTEGER, ALLOCATABLE, TARGET :: SetA(:), SetB(:)
  TYPE(Matrix_t), POINTER :: StiffMatrix
  REAL(KIND=dp), POINTER :: ForceVector(:)
  TYPE(Variable_t), POINTER :: PotVar, FixVar, SolVar, FluxVar, LoadVar, DistVar
  TYPE(Variable_t), POINTER :: PotVarA,PotVarB,PotSelect,CoilIndexVar,CoilSetVar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params, CoilList 
  REAL(KIND=dp) :: CoilCrossSection,InitialCurrent, Coeff, val, x0
  REAL(KIND=dp), ALLOCATABLE :: DesiredCoilCurrent(:), DesiredCurrentDensity(:)
  LOGICAL :: Found, CoilClosed, CoilAnisotropic, UseDistance, FixConductivity, &
      NormalizeCurrent, FitCoil, SelectNodes, CalcCurr, NarrowInterface
  LOGICAL, ALLOCATABLE :: GotCurr(:), GotDens(:)
  REAL(KIND=dp) :: CoilCenter(3), CoilNormal(3), CoilTangent1(3), CoilTangent2(3), &
      MinCurr(3),MaxCurr(3),TmpCurr(3)
  INTEGER, ALLOCATABLE :: CoilIndex(:)


 !------------------------------------------------------------------------------

  CALL Info('CoilSolver','--------------------------------------')
  CALL Info('CoilSolver','Solving current distribution in a coil')
  CALL Info('CoilSolver','--------------------------------------')

  Params => GetSolverParams()

  nsize = SIZE( Solver % Variable % Values ) 
  StiffMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % Rhs
  Mesh => Solver % Mesh
  SolVar => Solver % Variable
  Perm => SolVar % Perm

  dim = CoordinateSystemDimension()

  ! The following definiotions are same for all coils, even IF there would be many 
  !--------------------------------------------------------------------------------
  FixConductivity = GetLogical( Params,'Coil Conductivity Fix', Found ) 

  CoilAnisotropic = GetLogical( Params,'Coil Anisotropic', Found )

  UseDistance = GetLogical( Params,'Use Wall Distance',Found)
  IF( UseDistance ) THEN
     DistVar => VariableGet( Mesh % Variables,'Wall Distance' )
     IF( .NOT. ASSOCIATED( DistVar ) ) THEN
        CALL Fatal('CoilSolver','> Wall Distance < not associated!')
     END IF
  END IF

  CoilClosed = GetLogical( Params,'Coil Closed', Found )
  IF( CoilClosed ) THEN
    CoilParts = 2
  ELSE
    IF( .NOT. ListGetLogicalAnyBC( Model,'Coil Start') ) THEN
      CALL Info('CoilSolver','Assuming coil that is not closed',Level=3)
      CALL Fatal('CoilSolver','> Coil Start < must be defined on some BC')
    END IF
    IF( .NOT. ListGetLogicalAnyBC( Model,'Coil End') ) THEN
      CALL Info('CoilSolver','Assuming coil that is not closed',Level=3)
      CALL Fatal('CoilSolver','> Coil End < must be defined on some BC')
    END IF
    CoilParts = 1
  END IF

  NarrowInterface = GetLogical( Params,'Narrow Interface',Found )
  
  CalcCurr = GetLogical( Params,'Calculate Coil Current',Found )
  IF( .NOT. Found ) CalcCurr = .TRUE.

  NormalizeCurrent = ListGetLogical( Params,'Normalize Coil Current',Found ) 

  IF( FixConductivity ) THEN
    IF( CoilAnisotropic ) THEN
      MaxNonlinIter = 3
    ELSE
      MaxNonlinIter = 2
    END IF
  ELSE
    MaxNonlinIter = 1
  END IF

  PotVarA => VariableGet( Mesh % Variables,'CoilPot' )
  ALLOCATE( SetA(nsize) )
  SetA = 0

  IF( .NOT. ASSOCIATED( PotVarA ) ) THEN
    CALL Fatal('CoilSolver','CoilPot not associated!')
  END IF
  IF( CoilParts == 2 ) THEN
    ALLOCATE( SetB(nsize) )
    SetB = 0 

    PotVarB => VariableGet( Mesh % Variables,'CoilPotB' )
    IF( .NOT. ASSOCIATED( PotVarB ) ) THEN
      CALL Fatal('CoilSolver','CoilPotB not associated!')
    END IF
    PotSelect => VariableGet( Mesh % Variables,'PotSelect' )
    IF( .NOT. ASSOCIATED( PotSelect ) ) THEN
      CALL Fatal('CoilSolver','PotSelect not associated!')
    END IF
  END IF

  IF( FixConductivity ) THEN
    FixVar => VariableGet( Mesh % Variables,'CoilFix' )
    IF( .NOT. ASSOCIATED( FixVar ) ) THEN
      CALL Fatal('CoilSolver','CoilFix not associated!')
    END IF
    FixVar % Values = 1.0_dp
  END IF
    
  ! Get the loads
  LoadVar => VariableGet( Mesh % Variables,&
      TRIM(SolVar % Name)//' Loads' )
  IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
    CALL Fatal('CoilSolver','> '//TRIM(SolVar % Name)//' < Loads not associated!')
  END IF

  MaxNoCoils = MAX( 1, Model % NumberOfComponents ) 
  ALLOCATE( DesiredCoilCurrent( MaxNoCoils ), DesiredCurrentDensity(MaxNoCoils), &
      GotCurr( MaxNoCoils ), GotDens( MaxNoCoils) )
  DesiredCoilCurrent = 0.0_dp
  GotCurr = .FALSE.
  GotDens = .FALSE.


  ! These are different for different coils, would there be many
  !-----------------------------------------------------------------------
  NoCoils = 0
  SelectNodes = .FALSE.
  DO i=1,Model % NumberOfComponents + 1
    IF( i <= Model % NumberOfComponents ) THEN
      CoilList => Model % Components(i) % Values

      IF(.NOT. ListCheckPresent( CoilList,'Coil Type' ) ) CYCLE      
      TargetBodies => ListGetIntegerArray( CoilList,'Master Bodies',Found )
      IF( .NOT. Found ) TargetBodies => ListGetIntegerArray( CoilList,'Body',Found )
      IF( .NOT. Found ) CALL Fatal('CoilSolver','Coil fitting requires > Master Bodies <') 

      CALL Info('CoilSolver','Treating coil in Component: '//TRIM(I2S(i)),Level=7)

      IF(.NOT. ALLOCATED( CoilIndex ) ) THEN
        ALLOCATE( CoilIndex( Mesh % NumberOfNodes ) )
        CoilIndex = 0
      END IF
      NoCoils = NoCoils + 1
      SelectNodes = .TRUE.
      CALL DefineCoilCenter( CoilCenter, CoilList, TargetBodies )
      CALL DefineCoilParameters( CoilNormal, CoilTangent1, CoilTangent2, &
          CoilList, TargetBodies )
    ELSE
      IF( NoCoils > 0 ) EXIT
      NoCoils = 1
      CoilList => Params
      CALL DefineCoilCenter( CoilCenter, Params )
      CALL DefineCoilParameters( CoilNormal, CoilTangent1, CoilTangent2, Params )
    END IF

    ! Choose nodes where the Dirichlet values are set. 
    IF( CoilClosed ) THEN
      IF( NarrowInterface ) THEN
        Set => SetA
        CALL ChooseFixedBulkNodesNarrow(Set,1,SelectNodes)
        Set => SetB
        CALL ChooseFixedBulkNodesNarrow(Set,2,SelectNodes)
      ELSE
        Set => SetA
        CALL ChooseFixedBulkNodes(Set,1,SelectNodes)
        Set => SetB
        CALL ChooseFixedBulkNodes(Set,2,SelectNodes)
      END IF
    ELSE
      Set => SetA
      CALL ChooseFixedEndNodes(Set)
    END IF

    DesiredCoilCurrent(NoCoils) = ListGetCReal( CoilList,'Desired Coil Current',Found )
    IF(.NOT. Found ) DesiredCoilCurrent(NoCoils) = 1.0_dp
    GotCurr(NoCoils) = Found

    DesiredCurrentDensity(NoCoils) = ListGetCReal( CoilList,'Desired Current Density',Found)
    IF(.NOT. Found ) DesiredCurrentDensity(NoCoils) = 1.0_dp
    GotDens(NoCoils) = Found


    ! If we know the coil cross section we can relate the wishes in total 
    ! current through cross section and current density. 
    CoilCrossSection = ListGetCReal( CoilList,'Coil Cross Section',Found ) 
    IF( Found ) THEN
      IF( GotCurr(NoCoils) ) THEN
        DesiredCurrentDensity(NoCoils) = DesiredCoilCurrent(NoCoils) / CoilCrossSection
        GotDens(NoCoils) = .TRUE.
      ELSE IF( GotDens(NoCoils) ) THEN
        DesiredCoilCurrent(NoCoils) = DesiredCurrentDensity(NoCoils) * CoilCrossSection 
        GotCurr(NoCoils) = .TRUE.
      END IF
    END IF
    
  END DO

  CALL Info('CoilSolver','Coil system consists of '//TRIM(I2S(NoCoils))//' coils',Level=7)


  ! Count the fixing nodes just for information 
  Set => SetA
  CALL CountFixingNodes(Set,1)
  IF( CoilClosed ) THEN
    Set => SetB
    CALL CountFixingNodes(Set,2)
  END IF
  
       
  DO iter=1,MaxNonlinIter
      
    IF( iter > 1 ) THEN
      CALL Info('CoilSolver','Fixing the conductivity field')
      
      CALL DefaultInitialize()
      
      Active = GetNOFActive()
      DO t=1,Active
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()           
        CALL LocalFixMatrix(  Element, n, nd )
      END DO
      
      ! Solve the fixing field
      !--------------------------
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.) 
      Norm = DefaultSolve()
      FixVar % Values = SolVar % Values
    END IF

    
    CALL Info('CoilSolver','Computing the dummy potential field')

    ! For closed coils the solution is computed in two parts
    DO Part = 1,CoilParts
      
      ! For closed coils the solution is computed in two parts
      IF( Part == 1 ) THEN        
        Set => SetA
        PotVar => PotVarA
      ELSE
        Set => SetB
        PotVar => PotVarB
      END IF
      
      CALL DefaultInitialize()

      Active = GetNOFActive()          
      DO t=1,Active
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        CALL LocalPotMatrix(  Element, n, nd )
      END DO
      
      ! This routine is needed to save the BulkValues
      ! Bulk values are needed to compute the currents
      CALL DefaultFinishBulkAssembly()
      CALL DefaultFinishAssembly()

      ! Set the potential values on the nodes
      ! Values 0/1 are used as the potential is later scaled. 
      DO i = 1,nsize
        sgn = Set(i)
        IF( sgn == 0 ) CYCLE
        
        IF( sgn > 0 ) THEN
          val = 1.0_dp
        ELSE
          val = 0.0_dp
        END IF
        
        !s = StiffMatrix % Values(StiffMatrix % Diag(i))
        !ForceVector(i) = val * s

        CALL UpdateDirichletDof( Solver % Matrix, i, val )
                
        !CALL ZeroRow( StiffMatrix,i )
        !CALL SetMatrixElement( StiffMatrix,i,i,1.0d0*s )     
      END DO
      
      ! If we use narrow strategy we need to cut the connections in the bulk values
      ! between the two different Dirichlet conditions. Otherwise the load computation
      ! will produce crap.
      IF( NarrowInterface ) THEN
        CALL CutInterfaceConnections( StiffMatrix, Set )
      END IF

      ! Only Default dirichlet conditions activate the BCs above!
      CALL DefaultDirichletBCs()

      
      ! Solve the potential field
      !--------------------------
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.FALSE.) 

      ! Here solve actually the potential
      !--------------------------------------
      Norm = DefaultSolve()
      
      PotVar % Values = SolVar % Values

      ! Get the nodal loads i.e. nodal currents in this case
      LoadVar => VariableGet( Mesh % Variables,&
          TRIM(SolVar % Name)//' Loads' )
      IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
        CALL Fatal('CoilSolver','> '//TRIM(SolVar % Name)//' < Loads not associated!')
      END IF

      CALL ScalePotential()
        
    END DO
  END DO


  ! Compute the current
  !---------------------------------
  IF( CalcCurr ) THEN
    CALL ListAddLogical( Params,'Calculate Loads',.FALSE.)

    MinCurr = 0.0_dp
    MaxCurr = 0.0_dp

    DO Part = 1,CoilParts

      ! The solution is computed separately for each component
      !--------------------------------------------------------

      DO dimi = 1,dim
        CALL Info('CoilSolver','Computing current component: '//TRIM(I2S(dimi)),Level=6)
        
        CALL DefaultInitialize()
        Active = GetNOFActive()          
        DO t=1,Active
          Element => GetActiveElement(t)
          n  = GetElementNOFNodes()
          nd = GetElementNOFDOFs()           
          CALL LocalFluxMatrix(  Element, n, nd, dimi )
        END DO
        
        ! Solve the flux in direction dimi
        !--------------------------------------
        CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.) 
        Norm = DefaultSolve()

        ! If the coil is computed in two parts then 
        ! for the other part pick only the values which are on the better half
        FluxVar => VariableGet( Mesh % Variables,'CoilCurrent '//TRIM(I2S(dimi)) )
        IF( .NOT. ASSOCIATED( FluxVar ) ) THEN
          CALL Fatal('CoilSolver','CoilCurrent not associated!')
        END IF

        MinCurr(dimi) = MINVAL( SolVar % Values ) 
        MaxCurr(dimi) = MAXVAL( SolVar % Values )

        FluxVar % Values = SolVar % Values
      END DO

      IF( ParEnv % PEs > 1 ) THEN
        TmpCurr = MinCurr
        CALL MPI_ALLREDUCE(MinCurr,TmpCurr,3,MPI_DOUBLE_PRECISION,MPI_MIN,ELMER_COMM_WORLD,ierr)
        MinCurr = TmpCurr
        TmpCurr = MaxCurr
        CALL MPI_ALLREDUCE(MaxCurr,TmpCurr,3,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
        MaxCurr = TmpCurr
      END IF

      WRITE( Message,'(A,3ES12.4)') 'Minimum current components: ',MinCurr
      CALL Info('CoilSolver',Message,Level=7)
      WRITE( Message,'(A,3ES12.4)') 'Maximum current components: ',MaxCurr
      CALL Info('CoilSolver',Message,Level=7)
    END DO

    IF( NormalizeCurrent ) THEN
      CALL NormalizeCurrentDensity() 
    END IF
    CALL ListAddLogical( Params,'Calculate Loads',.TRUE.)
  END IF
    

  ! Some optional postprocessing mainly for debugging purposes
  CoilIndexVar => VariableGet( Mesh % Variables,'CoilIndex' )
  IF( ASSOCIATED( CoilIndexVar ) ) THEN
    IF( .NOT. ALLOCATED( CoilIndex ) ) THEN
      CALL Warn('CoilSolver','CoilIndex requested for saving but it does not exist!')
    ELSE
      DO i=1,Mesh % NumberOfNodes
        j = CoilIndexVar % Perm(i)
        IF( j > 0 ) CoilIndexVar % Values(j) = 1.0_dp * CoilIndex(i)
      END DO
    END IF
  END IF

  ! Some optional postprocessing mainly for debugging purposes
  CoilSetVar => VariableGet( Mesh % Variables,'CoilSet' )
  IF( ASSOCIATED( CoilSetVar ) ) THEN
    DO i=1,Mesh % NumberOfNodes
      j = CoilSetVar % Perm(i)
      IF( j > 0 ) CoilSetVar % Values(j) = 1.0_dp * SetA(j) 
    END DO
  END IF

  IF( CoilParts == 2 ) THEN
    CoilSetVar => VariableGet( Mesh % Variables,'CoilSetB' )
    IF( ASSOCIATED( CoilSetVar ) ) THEN
      DO i=1,Mesh % NumberOfNodes
        j = CoilSetVar % Perm(i)
        IF( j > 0 ) CoilSetVar % Values(j) = 1.0_dp * SetB(j) 
      END DO
    END IF
  END IF

  ! Finally, always use the primary variable for testing convergence in
  ! coupled system level etc.
  Solver % Variable % Values = PotVarA % Values
  
  CALL Info('CoilSolver','All done',Level=7)
  CALL Info('CoilSolver','--------------------------------------')
   
 

CONTAINS 


!------------------------------------------------------------------------------
  SUBROUTINE CutInterfaceConnections( A, Set )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A             
    INTEGER, POINTER :: Set(:)
!------------------------------------------------------------------------------
    INTEGER, POINTER  CONTIG :: Cols(:),Rows(:)
    REAL(KIND=dp), POINTER  CONTIG :: Values(:)

    INTEGER :: i,j,k,n,cnt,s1,s2,dia

!------------------------------------------------------------------------------     
     n = A % NumberOfRows
     Rows   => A % Rows
     Cols   => A % Cols

     IF(.NOT. ASSOCIATED( A % BulkValues ) ) THEN
       CALL Fatal('CutInterfaceConnections','Dont have bulk values!')
     END IF
     Values => A % BulkValues

     cnt = 0
     
     DO i=1,n
       s1 = Set(i)
       IF( s1 == 0 ) CYCLE
       dia = A % Diag(i)
       DO j=Rows(i),Rows(i+1)-1
         k = Cols(j)
         s2 = Set(k)
         IF( s1 * s2 < 0 ) THEN
           ! The diagonal needs also to be compensated for the cutted connections.
           ! For Laplace operator the row sum is known to be zero. 
           A % BulkValues(dia) = A % BulkValues(dia) + A % BulkValues(j)
           A % BulkValues(j) = 0.0_dp
           cnt = cnt + 1
         END IF
       END DO
     END DO
     
     CALL Info('CutInterfaceConnections','Number of connections cut: '//TRIM(I2S(cnt)),Level=7)
     
!------------------------------------------------------------------------------
   END SUBROUTINE CutInterfaceConnections
!------------------------------------------------------------------------------


  

  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil.
  !----------------------------------------------------------------------------  
  SUBROUTINE DefineCoilCenter(CoilCenter, Params, TargetBodies)
    REAL(KIND=dp) :: CoilCenter(3)
    TYPE(ValueList_t), POINTER :: Params 
    INTEGER, POINTER, OPTIONAL :: TargetBodies(:)

    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: DetJ,r(3),s
    INTEGER :: e,t,i,j,n,Active
    LOGICAL :: stat,Found,CoilCenterSet
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Volume,Center(3),SerTmp(4),ParTmp(4),ierr
    REAL(KIND=dp), POINTER :: HelperArray(:,:)

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )


    ! If there will be many coils then mark no nodes associated to them
    IF( PRESENT( TargetBodies ) ) THEN
      Active = GetNOFActive()
      DO e=1,Active
        Element => GetActiveElement(e)
        n  = GetElementNOFNodes()        
        IF( ALL( TargetBodies /= Element % BodyId ) ) CYCLE
        CoilIndex( Element % NodeIndexes ) = NoCoils
      END DO
    END IF


    HelperArray => ListGetConstRealArray( Params, 'Coil Center', CoilCenterSet)
    IF( CoilCenterSet ) THEN
      Center(1:3) = HelperArray(1:3,1)
    ELSE
      Center(1) = ListGetCReal( Params,'Coil x0',CoilCenterSet)
      Center(2) = ListGetCReal( Params,'Coil y0',Found)
      CoilCenterSet = CoilCenterSet .OR. Found
      Center(3) = ListGetCReal( Params,'Coil z0',Found)
      CoilCenterSet = CoilCenterSet .OR. Found
    END IF
      
    IF( CoilCenterSet ) THEN
      CoilCenter = Center
      CALL Info('CoilSolver','Coil center defined by user',Level=20)
      RETURN
    END IF

    ! If coil center not given by user then compute the center of the coil
    !---------------------------------------------------------------------
    Volume = 0.0_dp
    Active = GetNOFActive()
    DO e=1,Active
      Element => GetActiveElement(e)
      n  = GetElementNOFNodes()
      
      IF( PRESENT( TargetBodies ) ) THEN
        IF( ALL( TargetBodies /= Element % BodyId ) ) CYCLE
      END IF

      CALL GetElementNodes( Nodes, Element )
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints(Element)
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )
        
        r(1) = SUM(Nodes % x(1:n)  *Basis(1:n))
        r(2) = SUM(Nodes % y(1:n) * Basis(1:n))
        r(3) = SUM(Nodes % z(1:n) * Basis(1:n))
        
        s = IP % s(t) * detJ
        
        Volume = Volume + s
        Center = Center + s * r 
      END DO
    END DO

    
    IF( ParEnv % PEs > 1 ) THEN
      SerTmp(1:3) = Center
      SerTmp(4) = Volume
      CALL MPI_ALLREDUCE(SerTmp,ParTmp,4,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
      Center = ParTmp(1:3)
      Volume = ParTmp(4)
    END IF

    IF( Volume < EPSILON( Volume ) ) CALL Fatal('DefineCoilCenter','Coil has no volume!')

    CoilCenter = Center / Volume
    
    WRITE( Message,'(A,ES12.4)') 'Coil volume:',Volume
    CALL Info('CoilSolver',Message,Level=7)

    WRITE( Message,'(A,3ES12.4)') 'Coil center:',CoilCenter
    CALL Info('CoilSolver',Message,Level=7)
    
  END SUBROUTINE DefineCoilCenter



  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil.
  !----------------------------------------------------------------------------  
  SUBROUTINE DefineCoilParameters(CoilNormal, CoilTangent1, CoilTangent2, &
      Params, TargetBodies )
    REAL(KIND=dp) :: CoilNormal(3), CoilTangent1(3), CoilTangent2(3)
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, POINTER, OPTIONAL :: TargetBodies(:)

    REAL(KIND=dp), POINTER :: HelperArray(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: DetJ,r(3),s,CoilTangentTmp(3)
    INTEGER :: e,t,i,j,n,Active
    LOGICAL :: stat,Found
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Imoment(9), EigVec(3,3), EigVal(3), ParTmp(9), ierr
    REAL(KIND=dp) :: EigWrk(20)
    INTEGER :: EigInfo, Three

    FitCoil = GetLogical( Params,'Fit Coil',Found )
    IF(.NOT. Found ) FitCoil = .TRUE. 

    HelperArray => ListGetConstRealArray( Params, 'Coil normal', Found)
    IF(.NOT. FitCoil .or. Found) THEN
      IF(Found) THEN
        CoilNormal(1:3) = HelperArray(1:3,1)
      ELSE
        CoilNormal = 0.0_dp
        CoilNormal(3) = 1.0_dp
      END IF
      CALL TangentDirections(CoilNormal, CoilTangent1, CoilTangent2)
      RETURN
    END IF

    CALL Info('DefineCoilParametes','Fitting the coil by maximizing inertia',Level=7)

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )

    Imoment = 0.0_dp

    Active = GetNOFActive()

    DO e=1,Active
      Element => GetActiveElement(e)
      
      IF( PRESENT( TargetBodies ) ) THEN
        IF( ALL( TargetBodies /= Element % BodyId ) ) CYCLE
      END IF

      n  = GetElementNOFNodes()
      
      CALL GetElementNodes( Nodes, Element )
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints(Element)
      DO t=1,IP % n
          ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )
        
        r(1) = SUM(Nodes % x(1:n)  *Basis(1:n))
        r(2) = SUM(Nodes % y(1:n) * Basis(1:n))
        r(3) = SUM(Nodes % z(1:n) * Basis(1:n))
        
        s = IP % s(t) * detJ
        r = r - CoilCenter
        DO i=1,3
          Imoment(3*(i-1)+i) = Imoment(3*(i-1)+i) + s * SUM( r**2 )
          DO j=1,3
            Imoment(3*(i-1)+j) = Imoment(3*(i-1)+j) - s * r(i) * r(j)
          END DO
        END DO
      END DO
    END DO

    IF( ParEnv % PEs > 1 ) THEN
      CALL MPI_ALLREDUCE(Imoment,ParTmp,9,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
      Imoment = ParTmp
    END IF

    DO i=1,3
      DO j=1,3
        EigVec(i,j) = Imoment(3*(i-1)+j)
      END DO
      EigVec(i,i) = EigVec(i,i) - 1.0_dp
    END DO

    EigInfo = 0
    Three = 3
    
    CALL DSYEV( 'V','U', Three, EigVec, Three, EigVal, EigWrk, SIZE(EigWrk), EigInfo )

    IF (EigInfo /= 0) THEN 
      CALL Fatal( 'CoilSolver', 'DSYEV cannot generate eigen basis')
    END IF

    WRITE( Message,'(A,3ES12.4)') 'Coil inertia eigenvalues:',EigVal
    CALL Info('CoilSolver',Message,Level=10)

    CoilNormal = EigVec(:,3)
    CoilTangent1 = EigVec(:,1)
    CoilTangent2 = EigVec(:,2)

    IF( -MINVAL( CoilTangent1 ) > MAXVAL( CoilTangent1 ) ) THEN
      CoilTangent1 = -CoilTangent1 
    END IF
    IF( -MINVAL( CoilTangent2 ) > MAXVAL( CoilTangent2 ) ) THEN
      CoilTangent2 = -CoilTangent2 
    END IF

    WRITE( Message,'(A,3ES12.4)') 'Coil axis normal:',CoilNormal
    CALL Info('CoilSolver',Message,Level=10)
    WRITE( Message,'(A,3ES12.4)') 'Coil tangent1:',CoilTangent1
    CALL Info('CoilSolver',Message,Level=10)
    WRITE( Message,'(A,3ES12.4)') 'Coil tangent2:',CoilTangent2
    CALL Info('CoilSolver',Message,Level=10)
    
  END SUBROUTINE DefineCoilParameters



  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil.
  !----------------------------------------------------------------------------  
  SUBROUTINE ChooseFixedBulkNodes( Set, SetNo, SelectNodes )
    
    INTEGER :: SetNo
    INTEGER, POINTER :: Set(:)
    LOGICAL :: SelectNodes

    LOGICAL :: Mirror 
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: x,y,z,x0,y0,dy
    REAL(KIND=dp) :: MinCoord(3),MaxCoord(3),r(3),rp(3),ParTmp(3),ierr
    INTEGER :: i,j,k,ioffset
    LOGICAL :: Found

    CALL Info('CoilSolver','Choosing fixing nodes for set: '//TRIM(I2S(SetNo)))

    Mirror = ( SetNo == 2 )

    Mesh => Solver % Mesh

    ! The maximum coordinate difference for an acceptable node
    ! The larger the value the more there will be nodes in the set.
    ! There should be enough nodes so that the BC is good one,
    ! but not too many either. 

    MinCoord = HUGE( MinCoord )
    MaxCoord = -HUGE( MaxCoord )
    ioffset = 10 * NoCoils 


    DO i=1,Mesh % NumberOfNodes
      IF( Perm(i) == 0 ) CYCLE

      IF( SelectNodes ) THEN
        IF( CoilIndex(i) /= NoCoils ) CYCLE
      END IF

      r(1) = Mesh % Nodes % x(i)
      r(2) = Mesh % Nodes % y(i)
      r(3) = Mesh % Nodes % z(i)

      ! Move to coil origin
      r = r - CoilCenter

      IF( mirror ) r = -r
      
      ! Coordinate projected to coil coordinates
      rp(1) = SUM( CoilTangent1 * r ) 
      rp(2) = SUM( CoilTangent2 * r ) 
      rp(3) = SUM( CoilNormal * r ) 
      
      DO j=1,3
        MinCoord(j) = MIN( MinCoord(j), rp(j) )
        MaxCoord(j) = MAX( MaxCoord(j), rp(j) ) 
      END DO
    END DO

    IF( ParEnv % PEs > 1 ) THEN
      CALL MPI_ALLREDUCE(MinCoord,ParTmp,3,MPI_DOUBLE_PRECISION,MPI_MIN,ELMER_COMM_WORLD,ierr)
      MinCoord = ParTmp
      CALL MPI_ALLREDUCE(MaxCoord,ParTmp,3,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
      MaxCoord = ParTmp
    END IF

    dy = ListGetCReal( Params,'Coil Bandwidth',Found)
    IF(.NOT. Found ) THEN
      dy = 0.2 * ( MaxCoord(2) - MinCoord(2) )
    END IF


    
    DO i=1,Mesh % NumberOfNodes

      IF( SelectNodes ) THEN
        IF( CoilIndex(i) /= NoCoils ) CYCLE
      END IF

      j = Perm(i)
      IF( j == 0 ) CYCLE
      
      r(1) = Mesh % Nodes % x(i)
      r(2) = Mesh % Nodes % y(i)
      r(3) = Mesh % Nodes % z(i)

      r = r - CoilCenter
      IF( mirror ) r = -r

      ! Coordinate projected to coil coordinates
      rp(1) = SUM( CoilTangent1 * r ) 
      rp(2) = SUM( CoilTangent2 * r ) 
      rp(3) = SUM( CoilNormal * r ) 

      IF( SetNo == 1 ) THEN
        ! This is used to determine "left" and "right" side of the coil
        PotSelect % Values( PotSelect % Perm(i) ) = rp(1)
      END IF

      IF( ABS( rp(2) ) > dy ) CYCLE

      ! Values with abs 1 indicate the narrow band that is omitted when computing the currents
      ! Values with abs 2 indicate the wide band
      IF( rp(1) > 0 ) THEN
        IF( rp(2) > dy / 2 ) THEN
          Set(j) = 2 + ioffset
        ELSE IF( rp(2) > 0.0 ) THEN
          Set(j) = 1 + ioffset
        ELSE IF( rp(2) > -dy / 2 ) THEN
          Set(j) = -1 - ioffset
        ELSE 
          Set(j) = -2 - ioffset
        END IF
      END IF

    END DO


  END SUBROUTINE ChooseFixedBulkNodes




  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil. Narrow version.
  !----------------------------------------------------------------------------  
  SUBROUTINE ChooseFixedBulkNodesNarrow( Set, SetNo, SelectNodes )
    
    INTEGER :: SetNo
    INTEGER, POINTER :: Set(:)
    LOGICAL :: SelectNodes

    LOGICAL :: Mirror 
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: x,y,z,x0,y0
    REAL(KIND=dp) :: r(3),rp(3),MinCut, MaxCut, CutDist(27)
    INTEGER :: t,i,j,k,n,ioffset
    LOGICAL :: Found, Hit
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: Indexes(:)


    CALL Info('CoilSolver','Choosing fixing nodes for set: '//TRIM(I2S(SetNo)))

    Mirror = ( SetNo == 2 )

    Mesh => Solver % Mesh

    ioffset = 10 * NoCoils 


    DO t=1,Mesh % NumberOfBulkElements

      Element => Mesh % Elements(t)
      Indexes => Element % NodeIndexes
      n = Element % Type % NumberOfNodes

      ! Study the elements belonging to the coil under study
      IF( SelectNodes ) THEN
        IF( ANY( CoilIndex(Indexes) /= NoCoils ) ) CYCLE
      END IF


      Hit = .TRUE.

      DO k = 1, n
        i = Indexes(k)

        j = Perm(i)
        IF( j == 0 ) THEN
          Hit = .FALSE.
          CYCLE
        END IF
          
        r(1) = Mesh % Nodes % x(i)
        r(2) = Mesh % Nodes % y(i)
        r(3) = Mesh % Nodes % z(i)

        r = r - CoilCenter
        IF( mirror ) r = -r

        ! Coordinate projected to coil coordinates
        rp(1) = SUM( CoilTangent1 * r ) 
        rp(2) = SUM( CoilTangent2 * r ) 
        rp(3) = SUM( CoilNormal * r ) 

        IF( SetNo == 1 ) THEN
          ! This is used to determine "left" and "right" side of the coil
          PotSelect % Values( PotSelect % Perm(i) ) = rp(1)
        END IF

        ! This element can not be an the interface as it is on the wrong side
        IF( rp(1) < 0 ) THEN
          Hit = .FALSE.
          CYCLE
        END IF

        CutDist(k) = rp(2)
      END DO
      
      IF(.NOT. Hit) CYCLE

      MaxCut = MAXVAL( CutDist(1:n) )
      MinCut = MINVAL( CutDist(1:n) )

      IF( MaxCut >= -EPSILON( MaxCut) .AND. MinCut <= EPSILON( MinCut)  ) THEN
        DO k=1,n
          i = Indexes(k)
          j = Perm(i)
          IF( CutDist(k) > 0.0_dp ) THEN
            Set(j) = 2 + ioffset
          ELSE
            Set(j) = -2 - ioffset
          END IF
        END DO
      END IF
    END DO

  END SUBROUTINE ChooseFixedBulkNodesNarrow


  

  ! Choose end nodes as assingled by "Coil Start" and "Coil End" flags.
  ! The result of this imitate the previous routine in order to be able
  ! to use the same way to computed the resulting currents.
  !--------------------------------------------------------------------
  SUBROUTINE ChooseFixedEndNodes( Set )
    
    INTEGER, POINTER :: Set(:)

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:)
    INTEGER :: t,i,j,k,nminus,nplus,ioffset
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found

    Set = 0
    ioffset = 10 * NoCoils

    Mesh => Solver % Mesh

    
    DO t=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       
       Element => Mesh % Elements(t)
       Model % CurrentElement => Element
       Indexes => Element % NodeIndexes

       IF( ANY( Perm( Indexes ) == 0 ) ) CYCLE

       BC => GetBC(Element)
       IF( ListGetLogical( BC,'Coil Start',Found ) ) THEN
          Set( Perm( Indexes ) ) = 2 + ioffset
       ELSE IF( ListGetLogical( BC,'Coil End',Found ) ) THEN
          Set( Perm( Indexes ) ) = -2 - ioffset
       END IF
    END DO
          
  END SUBROUTINE ChooseFixedEndNodes


  ! Count the nodes that will be used to set the fixing nodes 
  !----------------------------------------------------------
  SUBROUTINE CountFixingNodes( Set, SetNo )

    INTEGER, POINTER :: Set(:)
    INTEGER :: SetNo   
    INTEGER :: i,j,nplus,nminus
  

    ! Just count the nodes set 
    IF( ParEnv % PEs == 1 ) THEN
      nplus = COUNT( Set > 0  ) 
      nminus = COUNT( Set < 0  ) 
    ELSE
      nplus = 0
      nminus = 0
      DO i=1,nsize
        IF( Solver % Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1) &
            /= ParEnv % MyPe ) CYCLE
        IF( Set(i) > 0 ) THEN
          nplus = nplus + 1
        ELSE IF( Set(i) < 0 ) THEN
          nminus = nminus + 1
        END IF
      END DO
      nplus = NINT( ParallelReduction( 1.0_dp * nplus ) ) 
      nminus = NINT( ParallelReduction( 1.0_dp * nminus ) ) 
    END IF

    CALL Info('CoilSolver','Set'//TRIM(I2S(SetNo))//' : '&
        //TRIM(I2S(nplus))//' +nodes and ' &
        //TRIM(I2S(nminus))//' -nodes')

    IF( nplus == 0 .OR. nminus == 0 ) THEN
      CALL Warn('CoilSolver','Cannot set Dirichlet conditions with this set')
    END IF

  END SUBROUTINE CountFixingNodes




  ! Assembly the potential equation related to the dummy potential 
  !------------------------------------------------------------------------------
  SUBROUTINE LocalPotMatrix( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,GradAtIp(3), &
        FixAtIp, AbsGradAtIp, CondAtIp(3), DistGradAtIp(3),AbsDistGradAtIp,AbsCondAtIp
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), NodalPot(n), &
        NodalFix(n), NodalDist(n),DotProd, ElCond(n)
    LOGICAL :: Stat,Found, GotElCond
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: Material
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    Material => GetMaterial( Element ) 
    ElCond(1:n) = GetReal( Material, 'Electric Conductivity', GotElCond ) 
   
    IF( CoilParts == 1 ) THEN
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    ELSE      
      IF( MINVAL( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes)) ) > 0.0_dp ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( iter > 1 ) THEN
      NodalFix(1:n) = FixVar % Values( Perm( Element % NodeIndexes(1:n) ) )
    END IF

    IF( UseDistance ) THEN
      NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes(1:n) ) )
    END IF

    ! Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )

      CondAtIp = 1.0_dp
      IF( iter > 1 ) THEN
        FixAtIp = SUM( Basis(1:n) * NodalFix(1:n) )
        IF( CoilAnisotropic ) THEN
          DO i=1,3
            GradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalPot(1:n) )
          END DO
          AbsGradAtIp = SQRT( SUM( GradAtIp ** 2 ) )
          CondAtIp = ABS( GradAtIp ) / AbsGradAtIp
        END IF
      ELSE
        FixAtIp = 1.0_dp
      END IF

                
      ! If distance from the tangential walls is used then
      ! remove from conductivity the normal component.
      ! The normal vector is defined as the gradient of distance.
      ! As it can in center of coil be less than zero because it is 
      ! poorly defined, normalize the vector only if its larger. 
      ! At center the model could therefore be more isotropic.
      !-------------------------------------------------------------
      IF( UseDistance ) THEN
        DO i=1,3
          DistGradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalDist(1:n) )
        END DO
        AbsDistGradAtIp = SQRT( SUM( DistGradAtIp ** 2 ) )

        AbsCondAtIp = SQRT(SUM( CondAtIP**2 ))
        DotProd = SUM( CondAtIp * DistGradAtIp ) / ( AbsDistGradAtIp * SQRT(3.0) )

        IF( AbsDistGradAtIp > 1.0_dp ) THEN
          DistGradAtIp = DistGradAtIp / AbsDistGradAtIp
        END IF

        CondAtIp = CondAtIp - DotProd * DistGradAtIp
      END IF

      IF( GotElCond ) THEN
        CondAtIP = SUM( Basis(1:n) * ElCond(1:n) ) * FixAtIp * CondAtIp
      ELSE
        CondAtIP = FixAtIp * CondAtIp        
      END IF
        
      
      Weight = IP % s(t) * DetJ

      ! diffusion term (Cond*grad(u),grad(v)):
      ! -----------------------------------
      DO p=1,nd
        DO q=1,nd
          DO i=1,dim
            STIFF(p,q) = STIFF(p,q) + Weight * &
                CondAtIp(i) * dBasisdx(p,i) * dBasisdx(q,i) 
          END DO
        END DO
      END DO
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalPotMatrix
!------------------------------------------------------------------------------


! Assembly the equation for the fixing coefficient
!------------------------------------------------------------------------------
  SUBROUTINE LocalFixMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: diff_coeff(n), AbsGradAtIp, Weight, Dreg, AbsCondAtIp
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,GradAtIp(3),DistGradAtIp(3)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd),NodalPot(nd),NodalDist(nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    IF( CoilParts == 1 ) THEN
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    ELSE
      IF( MINVAL( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes)) ) > 0.0_dp ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( UseDistance ) THEN
       NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes ) )
    END IF
    
    diff_coeff(1:n) = GetReal(Params,'CFix Diffusion',Found)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      ! By construction the r.h.s. should be one
      !------------------------------------------
      LoadAtIP = 1.0_dp

      DO i=1,3
         GradAtIP(i) = SUM( dBasisdx( 1:nd, i) * NodalPot(1:nd) )
      END DO

      IF( UseDistance ) THEN
         DO i=1,3
            DistGradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalDist(1:n) )
         END DO
         GradAtIp = GradAtIp - SUM( GradAtIp * DistGradAtIp ) * DistGradAtIp 
      END IF

      AbsGradAtIp = SQRT( SUM( GradAtIp**2) )      

      Dreg = SUM( Basis(1:nd) * Diff_Coeff(1:nd) )
      
      Weight = IP % s(t) * DetJ

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
             Dreg * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * AbsGradAtIp * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalFixMatrix
!------------------------------------------------------------------------------


! Assembly the equation for the flux component
!------------------------------------------------------------------------------
  SUBROUTINE LocalFluxMatrix( Element, n, nd, dimi  )
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dimi
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),Weight,DetJ,LoadAtIP,FixAtIp, &
         GradAtIp(3),AbsGradAtIp,CondAtIp(3),DistGradAtIp(3),AbsDistGradAtIp, &
         AbsCondAtIp,DotProd
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd),NodalPot(nd),&
         NodalFix(nd),NodalDist(nd),Elcond(nd)
    LOGICAL :: Stat, GotElCond
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: Material

    SAVE Nodes
!------------------------------------------------------------------------------


    Material => GetMaterial( Element ) 
    ElCond(1:n) = GetReal( Material, 'Electric Conductivity', GotElCond ) 
    
    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    IF( CoilParts == 1 ) THEN
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    ELSE
      IF( MINVAL( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes)) ) > 0.0_dp ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( FixConductivity ) THEN
      NodalFix(1:n) = FixVar % Values( Perm( Element % NodeIndexes ) )
    END IF

    IF( UseDistance ) THEN
      NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes ) )
    END IF

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )


      DO i=1,3
        GradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalPot(1:n) )
      END DO
      
      CondAtIp = 1.0_dp
      IF( FixConductivity ) THEN
        FixAtIp = SUM( Basis(1:n) * NodalFix(1:n) )
        IF( CoilAnisotropic ) THEN
          AbsGradAtIp = SQRT( SUM( GradAtIp ** 2 ) )
          CondAtIp = ABS( GradAtIp ) / AbsGradAtIp
        END IF
      ELSE
        FixAtIp = 1.0_dp
      END IF

      IF( UseDistance ) THEN
        DO i=1,3
          DistGradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalDist(1:n) )
        END DO
        AbsDistGradAtIp = SQRT( SUM( DistGradAtIp ** 2 ) )

        AbsCondAtIp = SQRT(SUM( CondAtIP**2 ))
        DotProd = SUM( CondAtIp * DistGradAtIp ) / ( AbsDistGradAtIp * SQRT(3.0) )

        IF( AbsDistGradAtIp > 1.0_dp ) THEN
          DistGradAtIp = DistGradAtIp / AbsDistGradAtIp
        END IF

        CondAtIp = CondAtIp - DotProd * DistGradAtIp
      END IF

      
      IF( GotElCond ) THEN
        CondAtIP = SUM( Basis(1:n) * ElCond(1:n) ) * FixAtIp * CondAtIp
      ELSE
        CondAtIP = FixAtIp * CondAtIp        
      END IF
  
      
      LoadAtIp = CondAtIp(dimi) * GradAtIp(dimi) 

      Weight = IP % s(t) * DetJ

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

!    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalFluxMatrix
!------------------------------------------------------------------------------


! Scale the potential such that the total current is the desired one.
! The scaling uses the feature where the nodal current is computed from
! r=Ax-b where the components of r are the nodal currents corresponding to
! the enforces Dirichlet conditions.
!----------------------------------------------------------------------------- 
  SUBROUTINE ScalePotential() 
    
    REAL(KIND=dp) :: InitialCurrent,possum, negsum, sumerr
    INTEGER :: i,j,k,Coil,nsize,posi,negi


    nsize = SIZE( LoadVar % Perm ) 

    DO Coil = 1, NoCoils 

      ! Evaluate the positive and negative currents
      ! As the total current vanishes these should be roughly the same 
      ! but with different signs. 
      possum = 0.0_dp
      negsum = 0.0_dp
      posi = 0
      negi = 0

      DO i=1,nsize
        j = LoadVar % Perm(i) 
        IF( j == 0 ) CYCLE

        IF( ParEnv % PEs > 1 ) THEN
          IF( Solver % Matrix % ParallelInfo % Neighbourlist(j) % Neighbours(1) &
              /= ParEnv % MyPe ) CYCLE
        END IF
        
        IF( NoCoils > 1 ) THEN
          IF( CoilIndex(i) /= Coil ) CYCLE
        END IF

        ! Note that the narrow part of the gap is omitted and here only 
        ! currents related to the outher parts of the gap (with |sgn|=2) 
        ! are accounted for. 
        sgn = Set(j)
        IF( MODULO(sgn,10) == 2  ) THEN
          possum = possum + LoadVar % Values(j)
          posi = posi + 1
        ELSE IF( MODULO(sgn,10)-10 == -2 ) THEN
          negsum = negsum + LoadVar % Values(j)
          negi = negi + 1
        END IF
      END DO
    
      IF( ParEnv % PEs > 1 ) THEN
        possum = ParallelReduction( possum ) 
        negsum = ParallelReduction( negsum ) 
      END IF

      WRITE (Message,'(A,I6,ES12.4)') 'Positive coil currents:',posi,possum
      CALL Info('CoilSolver',Message,Level=12)

      WRITE (Message,'(A,I6,ES12.4)') 'Negative coil currents:',negi,negsum
      CALL Info('CoilSolver',Message,Level=12)
    
      IF( ABS( possum ) < EPSILON( possum ) ) THEN
        CALL Fatal('CoilSolver','No positive current sources on coil end!')
      END IF
      IF( ABS( negsum ) < EPSILON( negsum ) ) THEN
        CALL Fatal('CoilSolver','No negative current sources on coil end!')
      END IF            
      
      sumerr = 2.0 * ABS( ABS( possum ) - ABS( negsum ) )  / (ABS( possum ) + ABS( negsum ) ) 
      WRITE (Message,'(A,ES12.4)') 'Discrepancy of start and end coil currents: ',sumerr 
      CALL Info('CoilSolver',Message,Level=7)

      IF( sumerr > 0.5 ) THEN
        CALL Warn('CoilSolver','Positive and negative sums differ too much!')
      END IF

      InitialCurrent = ( possum - negsum ) / 2.0
      
      WRITE( Message,'(A,ES12.4)') 'Initial coil current for coil '&
          //TRIM(I2S(Coil))//':',InitialCurrent
      CALL Info('CoilSolver',Message,Level=5)
    

      ! Scale the potential such that the current is as desired
      ! The current scales linearly with the potential.
      Coeff = DesiredCoilCurrent(Coil) / InitialCurrent

      WRITE( Message,'(A,ES12.4)') 'Coil potential multiplier:',Coeff
      CALL Info('CoilSolver',Message,Level=5)

      IF( NoCoils == 1 ) THEN
        PotVar % Values = Coeff * PotVar % Values     
      ELSE
        DO i=1,nsize
          IF( CoilIndex(i) == Coil ) THEN
            j = PotVar % Perm(i) 
            IF( j == 0 ) CYCLE            
            PotVar % Values( j ) = Coeff * PotVar % Values( j )
          END IF
        END DO
      END IF
    END DO

  END SUBROUTINE ScalePotential


  
! Normalize the current density to a given length
! When using this feature the potential can no longer be used
! to obtain the same effective current.   
!------------------------------------------------------------------------------
  SUBROUTINE NormalizeCurrentDensity()
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dimi
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Weight,DetJ,LocalCurr(3),AbsCurr,ScaleCurr,&
        TotCurr,TotVol,TargetDensity
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),NodalCurr(:,:)
    LOGICAL :: Stat
    INTEGER :: i,elem,t,Coil
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes
!------------------------------------------------------------------------------

    CALL Info('CoilSolver','Normalizing current density to a constant value')

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n), NodalCurr(3,n) )
    
    FluxVar => VariableGet( Mesh % Variables,'CoilCurrent')
    IF( .NOT. ASSOCIATED( FluxVar ) ) THEN
      CALL Fatal('CoilSolver','CoilCurrent not associated!')
    END IF
    

    DO Coil = 1, NoCoils 

      ! Current is already scaled, but we don't know the desired current density
      ! Let's assume that the average current density is a good candidate
      !--------------------------------------------------------------------------
      IF( .NOT. GotDens(Coil) .AND. GotCurr(Coil) ) THEN
        
        TotVol = 0.0_dp
        TotCurr = 0.0_dp
        NodalCurr = 0.0_dp
        LocalCurr = 0.0_dp
        
        Active = GetNOFActive()          
        DO elem=1,Active
          Element => GetActiveElement(elem)
          n  = GetElementNOFNodes()
          nd = GetElementNOFDOFs()           
          
          IF( NoCoils > 1 ) THEN
            IF( ANY( CoilIndex( Element % NodeIndexes ) /= Coil ) ) CYCLE
          END IF

          CALL GetElementNodes( Nodes )
          
          DO dimi=1,dim
            NodalCurr(dimi,1:n) = FluxVar % Values( dim*(Perm( Element % NodeIndexes )-1) + dimi )
          END DO
          
          ! Numerical integration:
          !----------------------
          IP = GaussPoints( Element )
          DO t=1,IP % n
            stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )
            
            DO dimi=1,dim
              LocalCurr(dimi) = SUM( Basis(1:n) * NodalCurr(dimi,1:n) )
            END DO
            AbsCurr = SQRT( SUM( LocalCurr(1:dim)**2 ) )
            
            Weight = IP % s(t) * DetJ
            
            TotVol = TotVol + Weight
            TotCurr = TotCurr + Weight * AbsCurr
          END DO
        END DO

        IF( ParEnv % PEs > 1 ) THEN
          TotCurr = ParallelReduction( TotCurr ) 
          TotVol = ParallelReduction( TotVol ) 
        END IF
        
        TargetDensity = TotCurr / TotVol
        
        WRITE( Message,'(A,ES12.4)') 'Average current density:',TargetDensity
        CALL Info('CoilSolver',Message)
      ELSE
        TargetDensity = DesiredCurrentDensity(Coil)
        WRITE( Message,'(A,ES12.4)') 'Desired current density:',TargetDensity
      END IF


      
      ! Now perform the normalization to the target value
      DO i = 1, Solver % Mesh % NumberOfNodes
        IF( NoCoils > 1 ) THEN
          IF( CoilIndex(i) /= Coil ) CYCLE
        END IF
        j = FluxVar % Perm(i)
        IF( j == 0 ) CYCLE
      
        LocalCurr(1:dim) = FluxVar % Values( dim*(j-1)+1: dim*(j-1)+dim )
        AbsCurr = SQRT( SUM( LocalCurr(1:dim) ** 2 ) )
        IF( AbsCurr > TINY( AbsCurr ) ) THEN
          ScaleCurr = TargetDensity / AbsCurr 
          FluxVar % Values( dim*(j-1)+1: dim*(j-1)+dim ) = &
              ScaleCurr * LocalCurr(1:dim)
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE NormalizeCurrentDensity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CoilSolver
!------------------------------------------------------------------------------



! Returns the coil potential in the case of closed coil if it should be rather 
! used as a potential value than a current. The potential consists of two parts
! and is not continuous but it's gradient should be continuous. 
!------------------------------------------------------------------------------

FUNCTION CoilPotential( Model, n, t ) RESULT(f)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f

  LOGICAL :: Visited = .FALSE., Found
  TYPE( Variable_t), POINTER :: PotA, PotB, PotP, PotSelect => NULL()
  TYPE(Element_t), POINTER :: Element, PrevElement => NULL()
  REAL(KIND=dp) :: xmin

  SAVE Visited, PotA, PotB, PotP, PotSelect, PrevElement

  IF( .NOT. Visited ) THEN
    PotA => VariableGet( Model % Mesh % Variables,'CoilPot' )
    PotP => PotA

    PotB => VariableGet( Model % Mesh % Variables,'CoilPotB' )
    PotSelect => VariableGet( Model % Mesh % Variables,'PotSelect' )

    Visited = .TRUE.
  END IF

  IF( ASSOCIATED( PotSelect ) ) THEN
    Element => Model % CurrentElement
    IF( .NOT. ASSOCIATED( Element, PrevElement ) ) THEN
      ! One could use as well max or mean, for example
      ! Consistancy is most important    
      xmin = MINVAL( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes ) ) ) 
      IF( xmin > 0.0 ) THEN
        PotP => PotB
      ELSE
        PotP => PotA
      END IF
    END IF
  END IF

  f = PotP % Values( PotP % Perm(n) )

END FUNCTION CoilPotential


! Returns the potential such that it is normalized to the desired current density.
!---------------------------------------------------------------------------------
FUNCTION CoilPotentialNormalized( Model, n, t ) RESULT(f)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f

  LOGICAL :: Visited = .FALSE., Found, Stat
  TYPE( Variable_t), POINTER :: PotA, PotB, PotP, PotSelect => NULL()
  TYPE( Nodes_t) :: Nodes
  INTEGER :: i,j,m
  TYPE(Element_t), POINTER :: Element, PrevElement => NULL()
  REAL(KIND=dp) :: xmin, DesiredCurrentDensity, u, v, w, NormCoeff, detJ, gradPot(3)
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), NodalPot(:)
  TYPE(GaussIntegrationPoints_t) :: IP

  SAVE Visited, PotA, PotB, PotP, PotSelect, Nodes, &
      PrevElement, xmin, DesiredCurrentDensity, Basis, dBasisdx, NodalPot, &
      GradPot, NormCoeff

  IF( .NOT. Visited ) THEN

    PotB => VariableGet( Model % Mesh % Variables,'CoilPotB' )
    PotA => VariableGet( Model % Mesh % Variables,'CoilPot' )
    PotP => PotA
    PotSelect => VariableGet( Model % Mesh % Variables,'PotSelect' )
   
    DO i=1,Model % NumberOfSolvers
      DesiredCurrentDensity = ListGetCReal( Model % Solvers(i) % Values,&
          'Desired Current Density',Found )
      IF( Found ) EXIT
    END DO
    IF(.NOT. Found ) DesiredCurrentDensity = 1.0_dp

    n = Model % Solver % Mesh % MaxElementNodes
    ALLOCATE( Basis(n), dBasisdx(n,3), NodalPot(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    Visited = .TRUE.
  END IF

  Element => Model % CurrentElement

  IF( .NOT. ASSOCIATED( Element, PrevElement ) ) THEN
    m = GetElementNOFNodes()
        
    Nodes % x(1:m) = Model % Mesh % Nodes % x(Element % NodeIndexes)
    Nodes % y(1:m) = Model % Mesh % Nodes % y(Element % NodeIndexes)
    Nodes % z(1:m) = Model % Mesh % Nodes % z(Element % NodeIndexes)

    IP = GaussPoints( Element )
    u = SUM( IP % U ) / IP % n
    v = SUM( IP % V ) / IP % n
    w = SUM( IP % W ) / IP % n

    stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )

    IF( ASSOCIATED( PotSelect ) ) THEN
      ! One could use as well max or mean, for example
      ! Consistancy is most important
      xmin = MINVAL( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes ) ) ) 
      IF( xmin > 0.0_dp ) THEN
        PotP => PotB
      ELSE
        PotP => PotA
      END IF
    END IF

    NodalPot(1:m) = PotP % Values( PotP % Perm(Element % NodeIndexes) )
    
    DO i=1,3
      GradPot(i) = SUM( dBasisdx(1:m,i) * NodalPot(1:m) )
    END DO
    
    NormCoeff = DesiredCurrentDensity / SQRT( SUM(GradPot**2) )     
    PrevElement => Element
  END IF

  f = NormCoeff * PotP % Values( PotP % Perm(n) )

END FUNCTION CoilPotentialNormalized


