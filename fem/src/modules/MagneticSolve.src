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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  Solver for the MHD Maxwell equations (or the induction equation).
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE MagneticSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

    USE Maxwell
    USE MaxwellAxiS
    USE MaxwellGeneral
    USE DefUtils
    USE Differentials
!------------------------------------------------------------------------------

    IMPLICIT NONE

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: dt

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,l,n,nd,t,iter,LocalNodes,k1,k2,istat

     TYPE(ValueList_t),POINTER :: Material, Equation, BF, BC, SolverParams
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,Gravity(3), &
                      Tdiff,Normal(3),s,r,NewtonTol,NonlinearTol

     TYPE(Variable_t), POINTER :: ElectricSol, ForceSol, ElecFieldSol
     INTEGER :: NSDOFs,NewtonIter,NonlinearIter, dim

     REAL(KIND=dp), POINTER :: MagneticField(:),ElectricCurrent(:), &
      Work(:,:), M1(:),M2(:),M3(:),E1(:),E2(:),E3(:), &
      ForceVector(:), divB(:),ExB(:), LrF(:), LrF1(:), LrF2(:), LrF3(:), &
      ElecField(:), Field1(:), Field2(:), Field3(:)

     LOGICAL :: Stabilize,NewtonLinearization = .FALSE.,GotForceBC,GotIt

     INTEGER :: body_id,bf_id,eq_id
     INTEGER, POINTER :: NodeIndexes(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: VelocityVarName, VelocityVarOne
     CHARACTER(LEN=MAX_NAME_LEN) :: VelocityVarTwo, VelocityVarDre
     CHARACTER(LEN=MAX_NAME_LEN) :: ForceVarName

     LOGICAL :: AllocationsDone = .FALSE., FreeSurfaceFlag, UserDefinedVelo
     LOGICAL :: CalculateMagneticForce = .FALSE.

     REAL(KIND=dp),ALLOCATABLE:: MASS(:,:),STIFF(:,:), LoadVector(:,:),FORCE(:), &
      Conductivity(:),Mx(:),My(:),Mz(:),U(:),V(:),W(:),Alpha(:),Beta(:), &
         Permeability(:),ExBx(:),ExBy(:),ExBz(:),B1(:),B2(:),B3(:),MU(:),MV(:),MW(:)

     SAVE Mx,My,Mz,U,V,W,MASS,STIFF,LoadVector, &
       FORCE,ElementNodes,Alpha,Beta, &
         Conductivity, AllocationsDone,LocalNodes, &
           Permeability, divB, ExBx,ExBy,ExBz,B1,B2,B3,MU,MV,MW, &
            VelocityVarOne, VelocityVarTwo, VelocityVarDre


!------------------------------------------------------------------------------
     REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime


!------------------------------------------------------------------------------
!    Get variables needed for solving the system
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     dim = CoordinateSystemDimension()

     ElectricSol   => VariableGet( Model % Variables, 'Electric Current' )
     IF ( ASSOCIATED(ElectricSol) ) ElectricCurrent => ElectricSol % Values

     ForceVarName = GetString( Model % Solver % Values, &
          'Force Variable Name', GotIt )
     IF ( .NOT. GotIt )  ForceVarName = 'Lorentz Force'

     ForceSol => VariableGet( Model % Variables, ForceVarName )
     IF ( ASSOCIATED( ForceSol ) ) LrF => ForceSol % Values

     ElecFieldSol => VariableGet( Model % Variables, 'E Field' )
     IF ( ASSOCIATED( ElecFieldSol ) ) ElecField => ElecFieldSol % Values

     LocalNodes = COUNT( Solver % Variable % Perm>0 )

     UNorm = Solver % Variable % Norm
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------

     IF ( .NOT.AllocationsDone ) THEN
       N = Solver % Mesh % MaxElementDOFs

       ALLOCATE( U(N), V(N), W(N),    &
                 MU(N), MV(N), MW(N), &
                 MX(N), MY(N), MZ(N), &
                 ExBx(N),ExBy(N),ExBz(N), &
                 ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 Conductivity( N ),  &
                 Permeability( N ),  &
                 FORCE( 3*N ), &
                 MASS(  3*N,3*N ),  &
                 STIFF( 3*N,3*N ),  &
                 divB(Model % NumberOfNodes), &
                 B1(LocalNodes), &
                 B2(LocalNodes), &
                 B3(LocalNodes), &
                 LoadVector( 3,N ), Alpha( N ), Beta( N ),STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'MagneticSolve', 'Memory allocation error.' )
       END IF

       VelocityVarName = GetString( Solver % Values, &
            'Velocity Solution Name', GotIt )
       IF ( .NOT. GotIt )  VelocityVarName = 'Velocity'

       VelocityVarOne = TRIM(VelocityVarName) // ' 1'
       VelocityVarTwo = TRIM(VelocityVarName) // ' 2'
       VelocityVarDre = TRIM(VelocityVarName) // ' 3'

       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------

     SolverParams => GetSolverParams()

     Stabilize = GetLogical( SolverParams,'Stabilize',GotIt )
     IF ( .NOT.GotIt ) Stabilize = .TRUE.

     NonlinearTol = GetConstReal( SolverParams, &
        'Nonlinear System Convergence Tolerance' )

     NewtonTol = GetConstReal( SolverParams, &
        'Nonlinear System Newton After Tolerance' )

     NewtonIter = GetInteger( SolverParams,  &
        'Nonlinear System Newton After Iterations' )

     NonlinearIter = GetInteger( SolverParams, &
        'Nonlinear System Max Iterations' )

!------------------------------------------------------------------------------
!    Check if free surfaces present
!------------------------------------------------------------------------------
     FreeSurfaceFlag = .FALSE.
     DO i=1,Model % NumberOfBCs
       FreeSurfaceFlag = FreeSurfaceFlag.OR. &
          GetLogical( Model % BCs(i) % Values, 'Free Surface', GotIt )
       IF ( FreeSurfaceFlag ) EXIT
     END DO
!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'MagneticSolve', ' ', Level=4 )
       CALL Info( 'MagneticSolve', ' ', Level=4 )
       CALL Info( 'MagneticSolve', &
             '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'Magnetic induction iteration: ', iter
       CALL Info( 'MagneticSolve', Message, Level=4 )
       CALL Info( 'MagneticSolve', &
             '-------------------------------------', Level=4 )
       CALL Info( 'MagneticSolve', ' ', Level=4 )
!------------------------------------------------------------------------------
       CALL DefaultInitialize()

       DO t=1,Solver % NumberOFActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Model % NumberOfBulkElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
                       
           CALL Info( 'MagneticSolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         Element => GetActiveElement(t)
!
!------------------------------------------------------------------------------
         n  = GetElementNOFNodes()
!         nd = GetElementNOFDOFs()
         NodeIndexes => Element % NodeIndexes

         Equation => GetEquation()

         UserDefinedVelo = .FALSE.
         UserDefinedVelo = GetLogical(Equation,'User Defined Velocity', gotIt)
         
         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

         Material => GetMaterial()

         Permeability(1:n) = GetReal(Material, 'Magnetic Permeability'   )
         Conductivity(1:n) = GetReal(Material, 'Electrical Conductivity',GotIt )
         IF( GotIt ) THEN
           CALL Warn('MagneticSolve','Use electric conductivity instead of electrical')
         ELSE
           Conductivity(1:n) = GetReal(Material, 'Electric Conductivity' )
         END IF

         Mx = GetReal( Material, 'Applied Magnetic Field 1', Gotit )
         My = GetReal( Material, 'Applied Magnetic Field 2', Gotit )
         Mz = GetReal( Material, 'Applied Magnetic Field 3', Gotit )

         ExBx = 0.0d0
         ExBy = 0.0d0
         ExBz = 0.0d0
! If you want to use time-domain solution for high-frequendy part,
! leave external field out. Better to use frequency-domain solver!
!#if 1
         CALL GetScalarLocalSolution( ExBx, 'Magnetic Flux Density 1' )
         CALL GetScalarLocalSolution( ExBy, 'Magnetic Flux Density 2' )
         CALL GetScalarLocalSolution( ExBz, 'Magnetic Flux Density 3' )
!#endif

         Mx(1:n) = Mx(1:n) + ExBx(1:n)
         My(1:n) = My(1:n) + ExBy(1:n)
         Mz(1:n) = Mz(1:n) + ExBz(1:n)

         U  = 0.0d0
         V  = 0.0d0
         W  = 0.0d0
         MU = 0.0d0
         MV = 0.0d0
         MW = 0.0d0
! For high-f part (in time-domain), leave velocity contribution out.
!#if 1
         IF ( UserDefinedVelo ) THEN     
           ! check for given constant velocity
           U(1:n) = GetReal( Material, 'MHD Velocity 1', gotIt )
           V(1:n) = GetReal( Material, 'MHD Velocity 2', gotIt )
           W(1:n) = GetReal( Material, 'MHD Velocity 3', gotIt )
         ELSE
           CALL GetScalarLocalSolution( U, VelocityVarOne )
           CALL GetScalarLocalSolution( V, VelocityVarTwo )
           CALL GetScalarLocalSolution( W, VelocityVarDre )

           CALL GetScalarLocalSolution( MU, 'Mesh Velocity 1' )
           CALL GetScalarLocalSolution( MV, 'Mesh Velocity 2' )
           CALL GetScalarLocalSolution( MW, 'Mesh Velocity 3' )
         END IF
!#endif
         U = U - MU
         V = V - MV
         W = W - MW
!------------------------------------------------------------------------------
!        Set body forces
!------------------------------------------------------------------------------
         BF => getBodyForce()
 
         LoadVector = 0.0D0
         IF ( ASSOCIATED(BF) ) THEN
           LoadVector(1,1:n) = LoadVector(1,1:n) + GetReal( &
                   BF, 'Magnetic Bodyforce 1', GotIt )
           LoadVector(2,1:n) = LoadVector(2,1:n) + GetReal( &
                   BF, 'Magnetic Bodyforce 2', GotIt )
           LoadVector(3,1:n) = LoadVector(3,1:n) + GetReal( &
                   BF, 'Magnetic Bodyforce 3', GotIt )
         END IF
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
            CALL MaxwellCompose( &
                MASS,STIFF,FORCE, &
                    LoadVector,Conductivity*Permeability,Mx,My,Mz,U,V,W, &
                        Element,n,ElementNodes )
          ELSE IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
            CALL MaxwellAxiSCompose( &
                MASS,STIFF,FORCE, &
                    LoadVector,Conductivity*Permeability,Mx,My,Mz,U,V,W, &
                       Element,n,ElementNodes )
         ELSE 
            CALL MaxwellGeneralCompose( &
                MASS,STIFF,FORCE, &
                    LoadVector,Conductivity*Permeability,Mx,My,Mz,U,V,W, &
                        Element,n,ElementNodes )
         END IF

!------------------------------------------------------------------------------
!        If time dependent simulation, add mass matrix to global 
!        matrix and global RHS vector
!------------------------------------------------------------------------------
         IF ( TransientSimulation ) CALL Default1stOrderTime(MASS,STIFF,FORCE)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO

      CALL DefaultFinishBulkAssembly()
      CALL Info( 'MagneticSolve', 'Assembly done', Level=4 )

      at = CPUTime() - at
      st = CPUTime()

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      DO t = 1,Solver % Mesh % NumberOFBoundaryElements

        Element => GetBoundaryElement(t)
        IF ( GetElementFamily()==1 .OR. .NOT. ActiveBoundaryElement() ) CYCLE

        BC => GetBC()
        IF (.NOT. ASSOCIATED(BC) )  CYCLE

!        nd = GetElementNOFDOFs()
        n  = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

!------------------------------------------------------------------------------
!           (at the moment the following is done...)
!           BC: \tau \cdot n = \alpha n +  @\beta/@t + F
!------------------------------------------------------------------------------
        LoadVector = 0.0D0
        Alpha      = 0.0D0
        Beta       = 0.0D0

        GotForceBC = GetLogical( BC, 'Magnetic Force BC',gotIt )
        IF ( GotForceBC ) THEN
!------------------------------------------------------------------------------
!         normal force BC: \tau\cdot n = \alpha n
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!         tangential force BC:
!         \tau\cdot n = @\beta/@t (tangential derivative of something)
!------------------------------------------------------------------------------
!         force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------
          LoadVector(1,1:n) =  GetReal( BC, 'Magnetic Force 1',GotIt )
          LoadVector(2,1:n) =  GetReal( BC, 'Magnetic Force 2',GotIt )
          LoadVector(3,1:n) =  GetReal( BC, 'Magnetic Force 3',GotIt )

!------------------------------------------------------------------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
            CALL MaxwellBoundary( STIFF,FORCE, &
                LoadVector,Alpha,Beta,Element,n,ElementNodes )
          ELSE
            CALL MaxwellGeneralBoundary( STIFF,FORCE, &
                LoadVector,Alpha,Beta,Element,n,ElementNodes )
          END IF

!------------------------------------------------------------------------------
!         Update global matrices from local matrices
!------------------------------------------------------------------------------
          IF ( TransientSimulation ) THEN
            MASS = 0.0d0
            CALL Default1stOrderTime(MASS,STIFF,FORCE)
          END IF

          CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
        END IF
      END DO
!------------------------------------------------------------------------------

      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm
      UNorm = DefaultSolve()

      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE( Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info('MagneticSolve', Message, LEVEL=4 )
      WRITE( Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info('MagneticSolve', Message, LEVEL=4 )

!------------------------------------------------------------------------------
      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2 * ABS(PrevUNorm-UNorm)/(UNorm + PrevUNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm     : ',UNorm
      CALL Info( 'MagneticSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'MagneticSolve', Message, Level=4 )

      IF ( RelativeChange < NewtonTol .OR. &
             iter > NewtonIter ) NewtonLinearization = .TRUE.

     IF ( RelativeChange < NonLinearTol ) EXIT
    END DO

    M1 => Solver % Variable %  Values(1::3)
    M2 => Solver % Variable %  Values(2::3)
    M3 => Solver % Variable %  Values(3::3)

    IF ( ASSOCIATED( ElectricSol ) ) THEN
      E1 => ElectricCurrent(1::3)
      E2 => ElectricCurrent(2::3)
      E3 => ElectricCurrent(3::3)

      ! Compute the electric current from the magnetic field
      !------------------------------------------------------
      IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
        CALL AxiSCurl( M1,M2,M3,E1,E2,E3,Solver % Variable % Perm)
      ELSE
        CALL Curl( M1,M2,M3,E1,E2,E3,Solver % Variable % Perm )
      END IF
      CALL InvalidateVariable( Model % Meshes, Solver % Mesh,'Electric Current' )
    END IF

!------------------------------------------------------------------------------
!   Compute the Lorentz force if requested
!------------------------------------------------------------------------------

    IF ( ASSOCIATED( ForceSol ) ) THEN
       LrF1 => LrF(1::3)
       LrF2 => LrF(2::3)
       LrF3 => LrF(3::3)

       CALL LorentzForceNodal( LrF1,LrF2,LrF3,M1,M2,M3,Solver % Variable % Perm )
    END IF

!------------------------------------------------------------------------------
!   Compute the electric field if requested
!------------------------------------------------------------------------------

    IF ( ASSOCIATED( ElecFieldSol ) ) THEN

       IF ( ASSOCIATED( ElectricSol ) ) THEN

          Field1 => ElecField(1::3)
          Field2 => ElecField(2::3)
          Field3 => ElecField(3::3)

          CALL ComputeNodalField( Field1, Field2, Field3, &
               M1, M2, M3, E1, E2, E3, Solver % Variable % Perm )

       ELSE
          CALL Warn( 'MagneticSolve', &
               'Need also current as exported variable to compute electric field' )
       END IF

    END IF

!------------------------------------------------------------------------------

  CONTAINS

    SUBROUTINE ComputeNodalField( Field1, Field2, Field3, B1,B2,B3, Ji1, Ji2, Ji3, Reorder )

      IMPLICIT NONE

      REAL(KIND=dp) :: B1(:), B2(:), B3(:)
      REAL(KIND=dp) :: Ji1(:), Ji2(:), Ji3(:)
      REAL(KIND=dp) :: Field1(:), Field2(:), Field3(:)
      INTEGER :: Reorder(:)


      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
      TYPE(ValueList_t), POINTER :: Material

      REAL(KIND=dp), POINTER :: EConductivity(:), MPermeability(:)
      REAL(KIND=dp), POINTER :: ExBx(:), ExBy(:), ExBz(:)
      REAL(KIND=dp), POINTER :: ExMx(:), ExMy(:), ExMz(:)
      REAL(KIND=dp), POINTER :: TotBx(:), TotBy(:), TotBz(:)

      REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
      REAL(KIND=dp) :: EF(3)

      INTEGER, POINTER :: NodeIndexes(:),Visited(:)
      INTEGER :: p,q,i,t,n, M

      LOGICAL :: gotit
 
!------------------------------------------------------------------------------

     ALLOCATE( Visited(Model % NumberOfNodes) )

     M = Solver % Mesh % MaxElementNodes

     ALLOCATE( Nodes % x( M ), &
               Nodes % y( M ), &
               Nodes % z( M ) )

     ALLOCATE( EConductivity( M ), MPermeability( M ),  &
               ExBx( M ), ExBy( M ), ExBz( M ), &
               ExMx( M ), ExMy( M ), ExMz( M ), &
               TotBx( M ), TotBy( M ), TotBz( M ) )


     Visited = 0

     Field1 = 0.0d0
     Field2 = 0.0d0
     Field3 = 0.0d0

     DO t = 1, Solver % NumberOfActiveElements

        Element => Solver % Mesh % Elements( Solver % ActiveElements(t) )
        NodeIndexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes

        Nodes % x(1:n) = Solver % Mesh % Nodes % x( NodeIndexes )
        Nodes % y(1:n) = Solver % Mesh % Nodes % y( NodeIndexes )
        Nodes % z(1:n) = Solver % Mesh % Nodes % z( NodeIndexes )

        IF ( MINVAL( Reorder( NodeIndexes ) ) > 0 ) THEN

           Material => GetMaterial( Element )
           EConductivity(1:n) = GetReal( Material, 'Electrical Conductivity', &
                gotit, Element )
           IF( GotIt ) THEN 
             CALL Warn('MagenticSolve','Use electric conductivity instead of electrical')
           ELSE
             EConductivity(1:n) = GetReal( Material, 'Electric Conductivity', &
                gotit, Element )
           END IF

           MPermeability(1:n) = GetReal( Material, 'Magnetic Permeability', &
                gotit, Element )

!------------------------------------------------------------------------------

           ExMx = 0.0d0
           ExMy = 0.0d0
           ExMz = 0.0d0

           ExMx(1:n) = GetReal( Material, 'Applied Magnetic Field 1', Gotit )
           ExMy(1:n) = GetReal( Material, 'Applied Magnetic Field 2', Gotit )
           ExMz(1:n) = GetReal( Material, 'Applied Magnetic Field 3', Gotit )

           ExBx = 0.0d0
           ExBy = 0.0d0
           ExBz = 0.0d0

           CALL GetScalarLocalSolution( ExBx, 'Magnetic Flux Density 1' )
           CALL GetScalarLocalSolution( ExBy, 'Magnetic Flux Density 2' )
           CALL GetScalarLocalSolution( ExBz, 'Magnetic Flux Density 3' )

           ExMx(1:n) = ExMx(1:n) + ExBx(1:n)
           ExMy(1:n) = ExMy(1:n) + ExBy(1:n)
           ExMz(1:n) = ExMz(1:n) + ExBz(1:n)

           TotBx(1:n) = ExMx(1:n) + B1( Reorder( NodeIndexes ) )
           TotBy(1:n) = ExMy(1:n) + B2( Reorder( NodeIndexes ) )
           TotBz(1:n) = ExMz(1:n) + B3( Reorder( NodeIndexes ) )

!------------------------------------------------------------------------------

           U  = 0.0d0
           V  = 0.0d0
           W  = 0.0d0
           MU = 0.0d0
           MV = 0.0d0
           MW = 0.0d0

           IF ( UserDefinedVelo ) THEN     
              ! check for given constant velocity
              U(1:n) = GetReal( Material, 'MHD Velocity 1', gotIt )
              V(1:n) = GetReal( Material, 'MHD Velocity 2', gotIt )
              W(1:n) = GetReal( Material, 'MHD Velocity 3', gotIt )
           ELSE
              CALL GetScalarLocalSolution( U, VelocityVarOne )
              CALL GetScalarLocalSolution( V, VelocityVarTwo )
              CALL GetScalarLocalSolution( W, VelocityVarDre )

              CALL GetScalarLocalSolution( MU, 'Mesh Velocity 1' )
              CALL GetScalarLocalSolution( MV, 'Mesh Velocity 2' )
              CALL GetScalarLocalSolution( MW, 'Mesh Velocity 3' )
           END IF

           U = U - MU
           V = V - MV
           W = W - MW

!------------------------------------------------------------------------------

           DO p=1,n

              q = Reorder( NodeIndexes(p) )

              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 EF(1) = V(p)* TotBz(p) - W(p)*TotBy(p)
                 EF(2) = W(p)* TotBx(p) - U(p)*TotBz(p)
                 EF(3) = U(p)* TotBy(p) - V(p)*TotBx(p)

              ELSE IF ( CurrentCoordinateSystem()  == CylindricSymmetric ) THEN

                 CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,&
                      Nodes % x(p), Nodes % y(p), Nodes % z(p) )
 
                 EF(1) = W(p)*TotBy(p) - V(p)*TotBz(p)
                 EF(2) = U(p)*TotBz(p) - W(p)*TotBx(p)
                 ! You might want to use SI units for the azimuthal component,
                 ! if you compute electric field at nodal points and symmetry axis,
                 ! otherwise you divide by zero.
!#ifdef SI_UNITS
!                 EF(3) = V(p)*TotBx(p) - U(p)*TotBy(p)
!#else
                 IF ( SqrtMetric > 1.0d-10 ) THEN
                    EF(3) = ( V(p)*TotBx(p) - U(p)*TotBy(p) ) / SqrtMetric
                 ELSE
                    EF(3) = 0.d0
                 END IF
!#endif

              ELSE
                 CALL Warn( 'MagneticSolve', & 
                      'Unsupported coordinate system in computing electric field' )

                 EF = 0.0d0
              END IF

              Field1(q) = Field1(q) + &
                   Ji1(q) / ( EConductivity(p) * MPermeability(p) ) - EF(1)
              Field2(q) = Field2(q) + &
                   Ji2(q) / ( EConductivity(p) * MPermeability(p) ) - EF(2)
              Field3(q) = Field3(q) + &
                   Ji3(q) / ( EConductivity(p) * MPermeability(p) ) - EF(3)

              Visited(q) = Visited(q) + 1
           
           END DO
        END IF
     END DO

     DO i = 1, Model % NumberOfNodes
         IF ( Visited(i) > 1 ) THEN
            Field1(i) = Field1(i) / Visited(i)
            Field2(i) = Field2(i) / Visited(i)
            Field3(i) = Field3(i) / Visited(i)
         END IF
      END DO

      DEALLOCATE( Visited )
      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( EConductivity, MPermeability )
      DEALLOCATE( ExBx, ExBy, ExBz )
      DEALLOCATE( ExMx, ExMy, ExMz )
      DEALLOCATE( TotBx, TotBy, TotBz )

      PRINT*, 'And finishing...'

    END SUBROUTINE ComputeNodalField

!------------------------------------------------------------------------------

   SUBROUTINE LorentzForceNodal( LrF1,LrF2,LrF3,B1,B2,B3,Reorder )

     IMPLICIT NONE
     REAL(KIND=dp) :: B1(:),B2(:),B3(:)
     REAL(KIND=dp) :: LrF1(:),LrF2(:),LrF3(:), Lorentz(3)
     REAL(KIND=dp), POINTER :: Density(:)
     INTEGER :: Reorder(:)

     TYPE(Element_t), POINTER :: Element
     TYPE(Nodes_t) :: Nodes 

     LOGICAL :: Stat, Averaged

     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag

     INTEGER, POINTER :: NodeIndexes(:),Visited(:)
     INTEGER :: p,q,i,t,n, mat_id

     REAL(KIND=dp) :: u,v,w

!------------------------------------------------------------------------------

     ALLOCATE( Visited(Model % NumberOfNodes) )

     n = Model % Mesh % MaxElementNodes
     ALLOCATE(Nodes % x(n),Nodes % y(n),Nodes % z(n))
     ALLOCATE( Density(n) )

     Visited = 0

     LrF1 = 0.0d0
     LrF2 = 0.0d0
     LrF3 = 0.0d0

     DO t=1,Model % NumberOfBulkElements

        Element => Model % Elements(t)
        Bf_id = GetBodyForceId( Element, Stat )
        IF ( .NOT. Stat )  CYCLE

        Averaged = .FALSE.
        IF ( .NOT. GetLogical( Model % BodyForces(Bf_id) % Values, &
             'Lorentz Force', Stat ) )  THEN
           IF ( GetLogical( Model % BodyForces(bf_id) % Values, &
                'Averaged Lorentz Force', Stat ) ) THEN
              Averaged = .TRUE.
           ELSE
              CYCLE
           END IF
        END IF

        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        Nodes % x(1:n) = Model % Nodes % x( NodeIndexes )
        Nodes % y(1:n) = Model % Nodes % y( NodeIndexes )
        Nodes % z(1:n) = Model % Nodes % z( NodeIndexes )


        IF ( Averaged ) THEN

           mat_id = GetMaterialId( Element )
           Density(1:n) = GetReal( Model % Materials(mat_id) % Values , 'Density', &
                Stat, Element )

           CompressibilityFlag = ListGetString( Model % Materials(mat_id) % Values, &
                'Compressibility Model', Stat )
           IF ( Stat .AND. CompressibilityFlag /= 'incompressible' )  &
                CALL Fatal( 'MagneticSolve', &
                'Averaged Lorentz Force implemented only for incompressible flow' )

        END IF

        IF ( MINVAL(Reorder(NodeIndexes)) > 0 ) THEN

           DO p=1,n

              q = Reorder(NodeIndexes(p))
              u = Element % TYPE % NodeU(p)
              v = Element % TYPE % NodeV(p)

              IF ( Element % TYPE % DIMENSION == 3 ) THEN
                 w = Element % TYPE % NodeW(p)
              ELSE
                 w = 0.0D0
              END IF

! Call LorentzForce from here (and zero B_e, B_ac for high-f part)
! For r < 1.0d-10, avoid the 1/r term in ComputeLorentz (if CylindricSymmetric)
! by using SI_UNITS
              Lorentz = LorentzForce( Element,Nodes,u,v,w,n )
              IF ( Averaged ) THEN
                 LrF1(q) = LrF1(q) + Lorentz(1) / Density(p)
                 LrF2(q) = LrF2(q) + Lorentz(2) / Density(p)
                 LrF3(q) = LrF3(q) + Lorentz(3) / Density(p)
              ELSE
                 LrF1(q) = LrF1(q) + Lorentz(1)
                 LrF2(q) = LrF2(q) + Lorentz(2)
                 LrF3(q) = LrF3(q) + Lorentz(3)
              END IF

              Visited(q) = Visited(q) + 1
           
           END DO
        END IF
      END DO

      DO i=1,Model % NumberOfNodes
         IF ( Visited(i) > 1 ) THEN
            LrF1(i) = LrF1(i) / Visited(i)
            LrF2(i) = LrF2(i) / Visited(i)
            LrF3(i) = LrF3(i) / Visited(i)
         END IF
      END DO

      DEALLOCATE( Visited )
      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( Density )

    END SUBROUTINE LorentzForceNodal

!------------------------------------------------------------------------------
  END SUBROUTINE MagneticSolver
!------------------------------------------------------------------------------
