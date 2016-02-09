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
  LOGICAL :: Found  

  dim = CoordinateSystemDimension()
  Params => GetSolverParams()

  IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable','-nooutput CoilTmp')
  END IF
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
      'CoilPot')
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
      'CoilFix')

  IF( GetLogical( Params,'Coil Closed', Found ) ) THEN
     CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
          'CoilPotB')
  END IF

  CALL ListAddString( Params,&
       NextFreeKeyword('Exported Variable',Params),&
       'CoilCurrent[CoilCurrent:'//TRIM(I2S(dim))//']')

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
       MaxNonlinIter, dim,dimi,ierr
  INTEGER, POINTER :: Perm(:), Set(:)
  INTEGER, ALLOCATABLE, TARGET :: SetA(:), SetB(:)
  TYPE(Matrix_t), POINTER :: StiffMatrix
  REAL(KIND=dp), POINTER :: ForceVector(:)
  TYPE(Variable_t), POINTER :: PotVar, FixVar, SolVar, FluxVar, LoadVar, DistVar
  TYPE(Variable_t), POINTER :: PotVarA,PotVarB
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params 
  REAL(KIND=dp) :: possum, negsum, DesiredCoilCurrent, DesiredCurrentDensity, &
      CoilCrossSection,InitialCurrent, Coeff, val, x0
  LOGICAL :: Found, CoilClosed, CoilAnisotropic, UseDistance, FixConductivity, &
      NormalizeCurrent, GotCurr, GotDens, FitCoil
  REAL(KIND=dp) :: CoilCenter(3), CoilNormal(3), CoilTangent1(3), CoilTangent2(3), &
      MinCurr(3),MaxCurr(3),TmpCurr(3)


 !------------------------------------------------------------------------------

  CALL Info('CoilSolver','--------------------------------------')
  CALL Info('CoilSolver','Solving current distribution in a coil')
  CALL Info('CoilSolver','--------------------------------------')

  Params => GetSolverParams()

  nsize = SIZE( Solver % Variable % Values ) 
  ALLOCATE( SetA(nsize), SetB(nsize) )

  StiffMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % Rhs
  Mesh => Solver % Mesh
  SolVar => Solver % Variable
  Perm => SolVar % Perm

  dim = CoordinateSystemDimension()

  FixConductivity = GetLogical( Params,'Coil Conductivity Fix', Found ) 

  CoilAnisotropic = GetLogical( Params,'Coil Anisotropic', Found )


  CALl DefineCoilCenter( CoilCenter )

  FitCoil = GetLogical( Params,'Fit Coil',Found ) 
  IF(.NOT. Found ) FitCoil = .TRUE. 

  IF( FitCoil ) THEN
    CALL DefineCoilParameters( CoilNormal, CoilTangent1, CoilTangent2)
  ELSE
    CoilNormal   = 0.0_dp; CoilNormal(3)   = 1.0_dp
    CoilTangent1 = 0.0_dp; CoilTangent1(1) = 1.0_dp
    CoilTangent2 = 0.0_dp; CoilTangent2(2) = 1.0_dp
  END IF


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

  DesiredCoilCurrent = ListGetCReal( Params,'Desired Coil Current',GotCurr )
  IF(.NOT. GotCurr ) DesiredCoilCurrent = 1.0_dp

  DesiredCurrentDensity = ListGetCReal( Params,'Desired Current Density',GotDens)
  IF(.NOT. GotDens ) DesiredCurrentDensity = 1.0_dp

  ! If we know the coil cross section we can relate the wishes in total 
  ! current through cross section and current density. 
  CoilCrossSection = ListGetCReal( Params,'Coil Cross Section',Found ) 
  IF( Found ) THEN
    IF( GotCurr ) THEN
      DesiredCurrentDensity = DesiredCoilCurrent / CoilCrossSection
      GotDens = .TRUE.
    ELSE IF( GotDens ) THEN
      DesiredCoilCurrent = DesiredCurrentDensity * CoilCrossSection 
      GotCurr = .TRUE.
    END IF
  END IF

  NormalizeCurrent = ListGetLogical( Params,'Normalize Coil Current',Found ) 

  MaxNonlinIter = GetInteger( Params,'Nonlinear System Max Iterations',Found )
  IF(.NOT. Found ) THEN
    IF( FixConductivity ) THEN
      IF( CoilAnisotropic ) THEN
        MaxNonlinIter = 3
      ELSE
        MaxNonlinIter = 2
      END IF
    ELSE
      MaxNonlinIter = 1
    END IF
  END IF

  IF( MaxNonlinIter > 1 .AND. .NOT. FixConductivity ) THEN
    CALL Info('CoilSolver','No conductivity fixing, setting nonlinear iterations to one!',Level=6)
  END IF

  PotVarA => VariableGet( Mesh % Variables,'CoilPot' )
  IF( .NOT. ASSOCIATED( PotVarA ) ) THEN
    CALL Fatal('CoilSolver','CoilPot not associated!')
  END IF
  IF( CoilParts == 2 ) THEN
    PotVarB => VariableGet( Mesh % Variables,'CoilPotB' )
    IF( .NOT. ASSOCIATED( PotVarB ) ) THEN
      CALL Fatal('CoilSolver','CoilPotB not associated!')
    END IF
  END IF

  FixVar => VariableGet( Mesh % Variables,'CoilFix' )
  IF( .NOT. ASSOCIATED( FixVar ) ) THEN
    CALL Fatal('CoilSolver','CoilFix not associated!')
  END IF

  ! Get the loads
  LoadVar => VariableGet( Mesh % Variables,&
      TRIM(SolVar % Name)//' Loads' )
  IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
    CALL Fatal('CoilSolver','> '//TRIM(SolVar % Name)//' < Loads not associated!')
  END IF
  

  ! Choose nodes where the Dirichlet values are set. 
  IF( CoilClosed ) THEN
    Set => SetA
    CALL ChooseFixedBulkNodes(Set,1)
    Set => SetB
    CALL ChooseFixedBulkNodes(Set,2)
    x0 = ListGetCReal( Params,'Coil x0',Found )
  ELSE
    Set => SetA
    CALL ChooseFixedEndNodes(Set)
  END IF

       
  DO iter=1,MaxNonlinIter
      
    IF( iter == 1 ) THEN
      FixVar % Values = 1.0_dp
    ELSE
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
     
      !PRINT *,'Fix Range:',MINVAL( FixVar % Values), MAXVAL( FixVar % Values )
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
        
        s = StiffMatrix % Values(StiffMatrix % Diag(i))
        ForceVector(i) = val * s
        CALL ZeroRow( StiffMatrix,i )
        CALL SetMatrixElement( StiffMatrix,i,i,1.0d0*s )     
      END DO
      
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
      
      ! Evaluate the positive and negative currents
      ! As the total current vanishes these should be roughly the same 
      ! but with different signs. 
      possum = 0.0_dp
      negsum = 0.0_dp
      DO i=1,nsize
        sgn = Set(i)
        
        IF( ParEnv % PEs > 1 ) THEN
          IF( Solver % Matrix % ParallelInfo % Neighbourlist(i) % Neighbours(1) &
              /= ParEnv % MyPe ) CYCLE
        END IF

        ! Note that the narrow part of the gap is omitted and here only 
        ! currents related to the outher parts of the gap (with |sgn|=2) 
        ! are accounted for. 
        IF( sgn == 2 ) THEN
          possum = possum + LoadVar % Values(i)
        ELSE IF( sgn == -2 ) THEN
          negsum = negsum + LoadVar % Values(i)
        END IF
      END DO

      IF( ParEnv % PEs > 1 ) THEN
        possum = ParallelReduction( possum ) 
        negsum = ParallelReduction( negsum ) 
      END IF
   
      InitialCurrent = ( possum - negsum ) / 2.0

      WRITE( Message,'(A,ES12.4)') 'Initial coil current:',InitialCurrent
      CALL Info('CoilSolver',Message,Level=10)

      ! Scale the potential such that the current is as desired
      ! The current scales linearly with the potential.
      Coeff = DesiredCoilCurrent / InitialCurrent

      WRITE( Message,'(A,ES12.4)') 'Coil current multiplier:',Coeff
      CALL Info('CoilSolver',Message,Level=10)

      PotVar % Values = Coeff * PotVar % Values     
    END DO
  END DO

  ! Compute the current
  !---------------------------------
  MinCurr = 0.0_dp
  MaxCurr = 0.0_dp

  DO Part = 1,CoilParts
    
    ! The solution is computed separately for each component
    !--------------------------------------------------------
    DO dimi = 1,dim
      CALL Info('CoilSolver','Computing the current components')
      
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
      CALL MPI_ALLREDUCE(MinCurr,TmpCurr,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      MinCurr = TmpCurr
      TmpCurr = MaxCurr
      CALL MPI_ALLREDUCE(MaxCurr,TmpCurr,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
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


  DEALLOCATE( SetA, SetB ) 

  CALL Info('CoilSolver','All done',Level=7)
   
 

CONTAINS 

  

  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil.
  !----------------------------------------------------------------------------  
  SUBROUTINE DefineCoilCenter(CoilCenter)
    REAL(KIND=dp) :: CoilCenter(3)

    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: DetJ,r(3),s
    INTEGER :: e,t,i,j,n,Active
    LOGICAL :: stat,Found,CoilCenterSet
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Volume,Center(3),SerTmp(4),ParTmp(4),ierr
    

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )

    Volume = 0.0_dp
    Center = 0.0_dp

    Center(1) = ListGetCReal( Params,'Coil x0',CoilCenterSet)
    Center(2) = ListGetCReal( Params,'Coil y0',Found)
    CoilCenterSet = CoilCenterSet .OR. Found
    Center(3) = ListGetCReal( Params,'Coil z0',Found)
    CoilCenterSet = CoilCenterSet .OR. Found

    IF( CoilCenterSet ) THEN
      CoilCenter = Center
      CALL Info('CoilSolver','Coil center defined by user',Level=20)
      RETURN
    END IF

    ! If coil center not given by user then compute the center of the coil
    !---------------------------------------------------------------------
    Active = GetNOFActive()
    DO e=1,Active
      Element => GetActiveElement(e)
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
        
        Volume = Volume + s
        Center = Center + s * r 
      END DO
    END DO

    IF( ParEnv % PEs > 1 ) THEN
      SerTmp(1:3) = Center
      SerTmp(4) = Volume
      CALL MPI_ALLREDUCE(SerTmp,ParTmp,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      Center = ParTmp(1:3)
      Volume = ParTmp(4)
    END IF

    CoilCenter = Center / Volume
    
    WRITE( Message,'(A,ES12.4)') 'Coil volume:',Volume
    CALL Info('CoilSolver',Message,Level=7)

    WRITE( Message,'(A,3ES12.4)') 'Coil center:',CoilCenter
    CALL Info('CoilSolver',Message,Level=7)

    ! Add this also to the Simulation section as the UDF needs it
    CALL ListAddConstReal( Model % Simulation,'Coil x0',CoilCenter(1) )
    CALL ListAddConstReal( Model % Simulation,'Coil y0',CoilCenter(2) )
    CALL ListAddConstReal( Model % Simulation,'Coil z0',CoilCenter(3) )
    
  END SUBROUTINE DefineCoilCenter



  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil.
  !----------------------------------------------------------------------------  
  SUBROUTINE DefineCoilParameters(CoilNormal, CoilTangent1, CoilTangent2 )
    REAL(KIND=dp) :: CoilNormal(3), CoilTangent1(3), CoilTangent2(3)

    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: DetJ,r(3),s
    INTEGER :: e,t,i,j,n,Active
    LOGICAL :: stat,Found
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Imoment(9), EigVec(3,3), EigVal(3), ParTmp(9), ierr
    REAL(KIND=dp) :: EigWrk(20)
    INTEGER :: EigInfo, Three
    

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )

    Imoment = 0.0_dp

    Active = GetNOFActive()

    DO e=1,Active
      Element => GetActiveElement(e)
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
      CALL MPI_ALLREDUCE(Imoment,ParTmp,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
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

    ! In order to reproduce the stanard cartesian directions for the simple cases
    ! the sign of the tangent vectors is checked.
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
  SUBROUTINE ChooseFixedBulkNodes( Set, SetNo )
    
    INTEGER :: SetNo
    INTEGER, POINTER :: Set(:)

    LOGICAL :: Mirror 
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: x,y,z,x0,y0,rad2deg,fii,dfii,dy
    REAL(KIND=dp) :: MinCoord(3),MaxCoord(3),r(3),rp(3),ParTmp(3),ierr
    INTEGER :: i,j,k,nminus,nplus
    LOGICAL :: Found

    CALL Info('CoilSolver','Choosing fixing nodes for set: '//TRIM(I2S(SetNo)))

    Mirror = ( SetNo == 2 )

    rad2deg = 180.0_dp / PI

    Mesh => Solver % Mesh

    ! The angle of acceptable nodes 
    ! 90 degs effectively chooses the right half
    dfii = ListGetCReal( Params,'Coil dfii',Found)
    IF(.NOT. Found ) dfii = 90.0

    ! The maximum coordinate difference for an acceptable node
    ! The larger the value the more there will be nodes in the set.
    ! There should be enough nodes so that the BC is good one,
    ! but not too many either. 

    MinCoord = HUGE( MinCoord )
    MaxCoord = -HUGE( MaxCoord )

    DO i=1,Mesh % NumberOfNodes
      IF( Perm(i) == 0 ) CYCLE

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
      CALL MPI_ALLREDUCE(MinCoord,ParTmp,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      MinCoord = ParTmp
      CALL MPI_ALLREDUCE(MaxCoord,ParTmp,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      MaxCoord = ParTmp
    END IF

    dy = ListGetCReal( Params,'Coil Bandwidth',Found)
    IF(.NOT. Found ) THEN
      dy = 0.2 * ( MaxCoord(2) - MinCoord(2) )
    END IF

    Set = 0
    
    DO i=1,Mesh % NumberOfNodes
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

      fii = rad2deg * ATAN2( rp(1), rp(2) )
      IF( fii > 180.0 ) fii = fii - 360.0
      
      IF( ABS( fii ) > dfii ) CYCLE
      
      IF( ABS( rp(1) ) > dy ) CYCLE

      ! Values with abs 1 indicate the narrow band that is omitted when computing the currents
      ! Values with abs 2 indicate the wide band
      IF( rp(1) > 0 ) THEN
        Set(j) = 1
        IF( rp(1) > dy / 2 ) Set(j) = 2
      ELSE
        Set(j) = -1
        IF( y < -dy / 2 ) Set(j) = -2 
      END IF      
    END DO
    
    ! Just count the nodes set 
    IF( ParEnv % PEs == 1 ) THEN
      nplus = COUNT( Set > 0 ) 
      nminus = COUNT( Set < 0 ) 
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
      CALL Fatal('CoilSolver','Cannot set Dirichlet conditions with this set')
    END IF


  END SUBROUTINE ChooseFixedBulkNodes


  ! Choose end nodes as assingled by "Coil Start" and "Coil End" flags.
  ! The result of this imitate the previous routine in order to be able
  ! to use the same way to computed the resulting currents.
  !--------------------------------------------------------------------
  SUBROUTINE ChooseFixedEndNodes( Set )
    
    INTEGER, POINTER :: Set(:)

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:)
    INTEGER :: t,i,j,k,nminus,nplus
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found

    Set = 0

    Mesh => Solver % Mesh
    
    DO t=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       
       Element => Mesh % Elements(t)
       Model % CurrentElement => Element
       Indexes => Element % NodeIndexes

       IF( ANY( Perm( Indexes ) == 0 ) ) CYCLE

       BC => GetBC(Element)
       IF( ListGetLogical( BC,'Coil Start',Found ) ) THEN
          Set( Perm( Indexes ) ) = 2
       ELSE IF( ListGetLogical( BC,'Coil End',Found ) ) THEN
          Set( Perm( Indexes ) ) = -2 
       END IF
    END DO
          
    ! Just count the nodes set 
    nplus = COUNT( Set > 0 ) 
    nminus = COUNT( Set < 0 ) 

    IF( ParEnv % PEs > 1 ) THEN
      nplus = NINT( ParallelReduction( 1.0_dp * nplus ) ) 
      nminus = NINT( ParallelReduction( 1.0_dp * nminus ) ) 
    END IF

    CALL Info('CoilSolver','Found '&
        //TRIM(I2S(nplus))//' start nodes and ' &
        //TRIM(I2S(nminus))//' end nodes')

    IF( nplus == 0 .OR. nminus == 0 ) THEN
      CALL Fatal('CoilSolver','Cannot set Dirichlet conditions with this set')
    END IF

  END SUBROUTINE ChooseFixedEndNodes



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
         NodalFix(n), NodalDist(n),DotProd
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
      IF( MINVAL( Nodes % x(1:n)) - x0 > 0.0 ) THEN
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

      IF( iter > 1 ) THEN
         FixAtIp = SUM( Basis(1:n) * NodalFix(1:n) )
         IF( CoilAnisotropic ) THEN
            DO i=1,3
               GradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalPot(1:n) )
            END DO
            AbsGradAtIp = SQRT( SUM( GradAtIp ** 2 ) )
            CondAtIp = FixAtIp * ABS( GradAtIp ) / AbsGradAtIp
         ELSE
            CondAtIp = FixAtIp
         END IF
      ELSE
         CondAtIp = 1.0_dp
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

         IF( ANY( CondAtIp < 0 ) ) THEN
            PRINT *,'DistGrad',DistGradAtIp,AbsDistGradAtIp
            PRINT *,'Cond',CondAtIp,DotProd,AbsCondAtIp
         END IF

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
      IF( MINVAL( Nodes % x(1:n)) - x0 > 0.0 ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( ALL( ABS( NodalPot(1:n) ) < 1.0e-8 ) ) THEN
      PRINT *,'NodalPot',NodalPot
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

!    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
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
         NodalFix(nd),NodalDist(nd)
    LOGICAL :: Stat
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
      IF( MINVAL( Nodes % x(1:n)) - x0 > 0.0 ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( iter > 1 ) THEN
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
         GradAtIP(i) = SUM( dBasisdx( 1:nd, i) * NodalPot(1:nd) )
      END DO
      AbsGradAtIp = SQRT( SUM( GradAtIp**2) )            

      IF( iter == 1 ) THEN
        CondAtIp = 1.0_dp
      ELSE
        FixAtIP = SUM( Basis( 1:nd) * NodalFix(1:nd) )
        IF( CoilAnisotropic ) THEN
          CondAtIp = FixAtIp * ABS( GradAtIp ) / AbsGradAtIp
        ELSE
          CondAtIp = FixAtIp 
        END IF
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


! Normalize the current density to a given length
!------------------------------------------------------------------------------
  SUBROUTINE NormalizeCurrentDensity( )
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dimi
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight,DetJ,LocalCurr(3),AbsCurr,ScaleCurr,TotCurr,TotVol
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),NodalCurr(:,:)
    LOGICAL :: Stat
    INTEGER :: i,elem,t
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
    
    ! Current is already scaled, but we don't know the desired current density
    ! Let's assume that the average current density is a good candidate
    !--------------------------------------------------------------------------
    IF( .NOT. GotDens .AND. GotCurr ) THEN

      TotVol = 0.0_dp
      TotCurr = 0.0_dp
      NodalCurr = 0.0_dp
      LocalCurr = 0.0_dp

      Active = GetNOFActive()          
      DO elem=1,Active
        Element => GetActiveElement(elem)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()           
        
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

      DesiredCurrentDensity = TotCurr / TotVol
      
      WRITE( Message,'(A,ES12.4)') 'Average current density:',DesiredCurrentDensity
      CALL Info('CoilSolver',Message)
    END IF
    

    DO i=1, SIZE(FluxVar % Values ) / dim
      LocalCurr(1:dim) = FluxVar % Values( dim*(i-1)+1: dim*(i-1)+dim )
      AbsCurr = SQRT( SUM( LocalCurr(1:dim) ** 2 ) )
      IF( AbsCurr > TINY( AbsCurr ) ) THEN
        ScaleCurr = DesiredCurrentDensity / AbsCurr 
        FluxVar % Values( dim*(i-1)+1: dim*(i-1)+dim ) = &
            ScaleCurr * LocalCurr(1:dim)
      END IF
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

FUNCTION UnitedPotential( Model, n, t ) RESULT(f)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f

  LOGICAL :: Visited = .FALSE., Found
  TYPE( Variable_t), POINTER :: PotA, PotB
  TYPE( Nodes_t), POINTER :: Nodes
  INTEGER :: i,j,m
  TYPE(Element_t), POINTER :: Element, PrevElement => NULL()
  REAL(KIND=dp) :: xmin, r(3), rp(3), CoilCenter(3), CoilTangent(3)
  REAL(KIND=dp), POINTER :: WrkPntr(:,:) => NULL()
  LOGICAL :: GotCoilTangent

  SAVE Visited, PotA, PotB, Nodes, CoilCenter, CoilTangent, GotCoilTangent, &
      PrevElement, xmin

  IF( .NOT. Visited ) THEN
    PotB => VariableGet( Model % Mesh % Variables,'CoilPotB' )
    PotA => VariableGet( Model % Mesh % Variables,'CoilPot' )
    Nodes => Model % Mesh % Nodes

    WrkPntr => ListGetConstRealArray( Model % Simulation,'Coil Center', Found ) 
    IF( Found ) THEN
      CoilCenter(1:3) = WrkPntr(1:3,1)
    ELSE
      CoilCenter(1) = ListGetCReal( Model % Simulation,'Coil x0',Found )
      CoilCenter(2) = ListGetCReal( Model % Simulation,'Coil y0',Found )
      CoilCenter(3) = ListGetCReal( Model % Simulation,'Coil z0',Found )
    END IF
     
    WrkPntr => ListGetConstRealArray( Model % Simulation,'Coil Tangent', GotCoilTangent )
    IF( GotCoilTangent ) THEN
      CoilTangent(1:3) = WrkPntr(1:3,1)
    ELSE
      CoilTangent = 0.0_dp
      CoilTangent(1) = 1.0_dp
    END IF
   
    Visited = .TRUE.
  END IF

  Element => Model % CurrentElement

  IF( .NOT. ASSOCIATED( Element, PrevElement ) ) THEN
    m = GetElementNOFNodes()
    
    ! One could use as well max or mean, for example
    ! Consistancy is most important
    
    xmin = HUGE( xmin ) 
    DO i=1,m
      j = Element % NodeIndexes(i)
      r(1) = Model % Mesh % Nodes % x(j) 
      IF( GotCoilTangent ) THEN
        r(2) = Model % Mesh % Nodes % y(j) 
        r(3) = Model % Mesh % Nodes % z(j) 
        r = r - CoilCenter
        rp(1) = SUM( CoilTangent * r ) 
      ELSE
        rp(1) = r(1) - CoilCenter(1)
      END IF
      xmin = MIN( xmin, rp(1) )
    END DO
  END IF

  IF( xmin > 0.0 ) THEN
    f = PotB % Values( PotB % Perm(n) )
  ELSE
    f = PotA % Values( PotA % Perm(n) )
  END IF

END FUNCTION UnitedPotential


FUNCTION UnitedPotentialNormalized( Model, n, t ) RESULT(f)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f

  LOGICAL :: Visited = .FALSE., Found, Stat
  TYPE( Variable_t), POINTER :: PotA, PotB
  TYPE( Nodes_t) :: Nodes
  INTEGER :: i,j,m
  TYPE(Element_t), POINTER :: Element, PrevElement => NULL()
  REAL(KIND=dp) :: xmin, r(3), rp(3), CoilCenter(3), CoilTangent(3), &
      DesiredCurrentDensity, u, v, w, NormCoeff, detJ, gradPot(3)
  REAL(KIND=dp), POINTER :: WrkPntr(:,:) => NULL()
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), NodalPot(:)
  LOGICAL :: GotCoilTangent, UsePotB
  TYPE(GaussIntegrationPoints_t) :: IP

  SAVE Visited, PotA, PotB, Nodes, CoilCenter, CoilTangent, GotCoilTangent, &
      PrevElement, xmin, DesiredCurrentDensity, Basis, dBasisdx, NodalPot, &
      UsePotB

  IF( .NOT. Visited ) THEN
    PotB => VariableGet( Model % Mesh % Variables,'CoilPotB' )
    PotA => VariableGet( Model % Mesh % Variables,'CoilPot' )

    WrkPntr => ListGetConstRealArray( Model % Simulation,'Coil Center', Found ) 
    IF( Found ) THEN
      CoilCenter(1:3) = WrkPntr(1:3,1)
    ELSE
      CoilCenter(1) = ListGetCReal( Model % Simulation,'Coil x0',Found )
      CoilCenter(2) = ListGetCReal( Model % Simulation,'Coil y0',Found )
      CoilCenter(3) = ListGetCReal( Model % Simulation,'Coil z0',Found )
    END IF
     
    WrkPntr => ListGetConstRealArray( Model % Simulation,'Coil Tangent', GotCoilTangent )
    IF( GotCoilTangent ) THEN
      CoilTangent(1:3) = WrkPntr(1:3,1)
    ELSE
      CoilTangent = 0.0_dp
      CoilTangent(1) = 1.0_dp
    END IF
   
    DesiredCurrentDensity = ListGetCReal( Model % Simulation,'Desired Current Density',Found)
    IF(.NOT. Found ) DesiredCurrentDensity = 1.0_dp

    n = Model % Solver % Mesh % MaxElementNodes
    ALLOCATE( Basis(n), dBasisdx(n,3), NodalPot(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    Visited = .TRUE.
  END IF

  Element => Model % CurrentElement

  IF( .NOT. ASSOCIATED( Element, PrevElement ) ) THEN
    m = GetElementNOFNodes()
    
    ! One could use as well max or mean, for example
    ! Consistancy is most important
    
    Nodes % x(1:m) = Model % Mesh % Nodes % x(Element % NodeIndexes)
    Nodes % y(1:m) = Model % Mesh % Nodes % y(Element % NodeIndexes)
    Nodes % z(1:m) = Model % Mesh % Nodes % z(Element % NodeIndexes)

    xmin = HUGE( xmin ) 
    DO i=1,m
      r(1) = Nodes % x(i) 
      IF( GotCoilTangent ) THEN
        r(2) = Nodes % y(i) 
        r(3) = Nodes % z(i) 
        r = r - CoilCenter
        rp(1) = SUM( CoilTangent * r ) 
      ELSE
        rp(1) = r(1) - CoilCenter(1)
      END IF
      xmin = MIN( xmin, rp(1) )    
    END DO

    IP = GaussPoints( Element )
    u = SUM( IP % U ) / IP % n
    v = SUM( IP % V ) / IP % n
    w = SUM( IP % W ) / IP % n

    stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )
    UsePotB = ( xmin > 0.0 )
    IF( UsePotB ) THEN
      NodalPot(1:m) = PotB % Values( PotB % Perm(Element % NodeIndexes) )
    ELSE
      NodalPot(1:m) = PotA % Values( PotA % Perm(Element % NodeIndexes) )
    END IF
    
    DO i=1,3
      GradPot(i) = SUM( dBasisdx(1:m,i) * NodalPot(1:m) )
    END DO
    
    NormCoeff = DesiredCurrentDensity / SQRT( SUM(GradPot**2) )     
    PrevElement => Element
  END IF

  IF( UsePotB ) THEN
    f = NormCoeff * PotB % Values( PotB % Perm(n) )
  ELSE
    f = NormCoeff * PotA % Values( PotA % Perm(n) )
  END IF

END FUNCTION UnitedPotentialNormalized


FUNCTION PotentialNormalized( Model, n, t ) RESULT(f)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f

  LOGICAL :: Visited = .FALSE., Found, Stat
  TYPE( Variable_t), POINTER :: PotA
  TYPE( Nodes_t) :: Nodes
  INTEGER :: i,j,m
  TYPE(Element_t), POINTER :: Element, PrevElement => NULL()
  REAL(KIND=dp) :: DesiredCurrentDensity, u, v, w, NormCoeff, detJ, gradPot(3)
  REAL(KIND=dp), POINTER :: WrkPntr(:,:) => NULL()
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), NodalPot(:)
  TYPE(GaussIntegrationPoints_t) :: IP

  SAVE Visited, PotA, Nodes, PrevElement, DesiredCurrentDensity, Basis, dBasisdx, NodalPot

  IF( .NOT. Visited ) THEN
    PotA => VariableGet( Model % Mesh % Variables,'CoilPot' )

    DesiredCurrentDensity = ListGetCReal( Model % Simulation,'Desired Current Density',Found)
    IF(.NOT. Found ) DesiredCurrentDensity = 1.0_dp

    n = Model % Solver % Mesh % MaxElementNodes
    ALLOCATE( Basis(n), dBasisdx(n,3), NodalPot(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    Visited = .TRUE.
  END IF

  Element => Model % CurrentElement

  IF( .NOT. ASSOCIATED( Element, PrevElement ) ) THEN
    m = GetElementNOFNodes()
    
    ! One could use as well max or mean, for example
    ! Consistancy is most important
    
    Nodes % x(1:m) = Model % Mesh % Nodes % x(Element % NodeIndexes)
    Nodes % y(1:m) = Model % Mesh % Nodes % y(Element % NodeIndexes)
    Nodes % z(1:m) = Model % Mesh % Nodes % z(Element % NodeIndexes)

    IP = GaussPoints( Element )
    u = SUM( IP % U ) / IP % n
    v = SUM( IP % V ) / IP % n
    w = SUM( IP % W ) / IP % n

    stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )
    NodalPot(1:m) = PotA % Values( PotA % Perm(Element % NodeIndexes) )
    
    DO i=1,3
      GradPot(i) = SUM( dBasisdx(1:m,i) * NodalPot(1:m) )
    END DO
    
    NormCoeff = DesiredCurrentDensity / SQRT( SUM(GradPot**2) )     
    PrevElement => Element
  END IF

  f = NormCoeff * PotA % Values( PotA % Perm(n) )

END FUNCTION PotentialNormalized
