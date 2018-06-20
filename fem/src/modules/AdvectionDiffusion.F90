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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Ville Savolainen, Antti Pursula
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{

!------------------------------------------------------------------------------
!> Initialization of the main solver: AdvectionDiffusionSolver
!------------------------------------------------------------------------------
   SUBROUTINE AdvectionDiffusionSolver_init( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
     TYPE(Model_t), TARGET :: Model    !< All model information (mesh, materials, BCs, etc...)
     REAL(KIND=dp) :: Timestep         !< Timestep size for time dependent simulations
     LOGICAL :: TransientSimulation    !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Params

     Params => GetSolverParams()
     CALL ListAddInteger( Params,'Time Derivative Order', 1 )

   END SUBROUTINE AdvectionDiffusionSolver_init



!------------------------------------------------------------------------------
!>  Advection-diffusion equation solver for scalar fields. 
!------------------------------------------------------------------------------
   SUBROUTINE AdvectionDiffusionSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------

     USE SolverUtils
     USE Differentials
     USE DefUtils

! Will need this for Density later
     USE MaterialModels
! Need these for mass conservation check
     USE Integration

     IMPLICIT NONE
!------------------------------------------------------------------------------
 
     TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
     TYPE(Model_t), TARGET :: Model    !< All model information (mesh, materials, BCs, etc...)
     REAL(KIND=dp) :: Timestep         !< Timestep size for time dependent simulations
     LOGICAL :: TransientSimulation    !< Steady state or transient simulation
 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,m,n,pn,t,tmax,iter,istat,bf_id,CoordinateSystem, outbody
 
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     TYPE(Nodes_t)   :: ElementNodes, ParentNodes
     TYPE(Element_t),POINTER :: CurrentElement, Parent
     TYPE(Mesh_t), POINTER :: Mesh

     REAL(KIND=dp) :: Norm,RelativeChange
     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: Stabilize = .FALSE., Bubbles, GotIt, GotIt2, AbsoluteMass = .FALSE.
     LOGICAL :: AllocationsDone = .FALSE., ScaledToSolubility = .FALSE.
     LOGICAL :: ErrorWritten

     TYPE(ValueList_t), POINTER :: Material, BC, Eq

     INTEGER, POINTER :: SpeciesPerm(:), MeshPerm(:)
     REAL(KIND=dp), POINTER :: Species(:),ForceVector(:), MeshVelocity(:), Hwrk(:,:,:)
 
     REAL(KIND=dp), ALLOCATABLE :: LocalMassMatrix(:,:), SoretDiffusivity(:), &
       NEConst(:), LocalStiffMatrix(:,:),Load(:),Diffusivity(:,:,:), &
                   C0(:),C1(:),CT(:),C2(:,:,:),LocalForce(:), TimeForce(:)
     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, HeatSolName, PhiSolName, ConvectName
! Use C1, C2 as in users guide:
! Relative mass units: C1=Density, C2=Density*Diff, C0 = 0, Ct=C1
! Absolute mass units: C1=1, C2=Diff, C0 = div v, Ct=1
! Used to be for C_I, C_V
! Add C, when previous solution needed for linearization
     TYPE(Variable_t), POINTER :: TempSol,PhiSol,FlowSol,MeshSol

     INTEGER, POINTER :: TempPerm(:),PotentialPerm(:),FlowPerm(:)
     INTEGER :: NSDOFs,NonlinearIter,body_id,lbody,rbody,eq_id,MDOFs
     REAL(KIND=dp) :: Relax, dt

! For the moment assume a linear system
! LocalMassMatrix also added, since t-dep. possible
     REAL(KIND=dp) :: SpecificHeatRatio, ReferencePressure, Ratio, MaxSol, IonCharge
! Ratio => Ratio of solubilities, used in flux jump boundary condition

     REAL(KIND=dp), POINTER :: Temperature(:),EPotential(:),FlowSolution(:)
! We need also Density, but it must be calculated for Compressible
! We read in also gamma & c_p in *.sif
     REAL(KIND=dp), ALLOCATABLE :: U(:), V(:), W(:), LocalTemperature(:), &
         LocalPotential(:), Density(:), HeatCapacity(:), Pressure(:),GasConstant(:), &
         SpeciesTransferCoeff(:), SExt(:), MU(:), MV(:), MW(:)

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,totat,st,totst,at0
#else
     REAL(KIND=dp) :: at,totat,st,totst,CPUTime,at0,RealTime
#endif
     REAL(KIND=dp) :: Nrm(3), Nrm2(3)
     REAL(KIND=dp), ALLOCATABLE :: BackupForceVector(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, StabilizeFlag
     INTEGER :: CompressibilityModel
     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName
     CHARACTER(LEN=MAX_NAME_LEN) :: ConcentrationUnits
! We will read velocity, temperature and electrical potential in
! For the moment assume steady-state velocity, temperature and potential
! Then ElmerSolver reads them in automatically
! If time-dependent, have to somewhere:
! CALL LoadRestartFile( RestartFile,k,CurrentModel )

! Variables needed for the mass conservation check
     REAL(KIND=dp) :: Mass, PreviousMass
     INTEGER, DIMENSION( Model %  NumberOfBulkElements) :: ElementList

     SAVE LocalMassMatrix,LocalStiffMatrix,Load,C0,C1,CT,C2,Diffusivity, &
         TimeForce, LocalForce, ElementNodes,AllocationsDone, &
         U, V, W, LocalTemperature,LocalPotential,Density, HeatCapacity, Pressure, &
         GasConstant, SpeciesTransferCoeff, Hwrk, MU, MV, MW, Mass, &
         ParentNodes, SExt, SoretDiffusivity, NEConst

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     CALL Info( 'AdvectionDiffusion', ' ', Level=4 )
     CALL Info( 'AdvectionDiffusion', &
         '-------------------------------------', Level=4 )
     CALL Info( 'AdvectionDiffusion', 'Solving for species transport', Level=4 )
     CALL Info( 'AdvectionDiffusion', &
         '-------------------------------------', Level=4 )

     CoordinateSystem = CurrentCoordinateSystem()
     Mesh => Solver % Mesh

     Species     => Solver % Variable % Values
     SpeciesPerm => Solver % Variable % Perm

     IF ( ALL( SpeciesPerm == 0) ) RETURN

! These are available in Model after LoadRestartFile

     TempSol => VariableGet( Mesh % Variables, 'Temperature' )
     IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm    => TempSol % Perm
       Temperature => TempSol % Values
     END IF

     PhiSol => VariableGet( Mesh % Variables, 'Potential' )
     IF ( ASSOCIATED( PhiSol ) ) THEN
       PotentialPerm => PhiSol % Perm
       EPotential    => PhiSol % Values
     END IF

     FlowSol => VariableGet( Mesh % Variables, 'Flow Solution' )
     IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm     => FlowSol % Perm
       NSDOFs       =  FlowSol % DOFs
       FlowSolution => FlowSol % Values
     END IF

     MeshSol => VariableGet( Mesh % Variables, 'Mesh Velocity' )
     NULLIFY( MeshVelocity )
     IF ( ASSOCIATED(MeshSol ) ) THEN
       MDOFs    =  MeshSol % DOFs
       MeshPerm => MeshSol % Perm
       MeshVelocity => MeshSol % Values
     END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

     Norm = Solver % Variable % Norm
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Mesh % MaxElementNodes

! Solve species weakly coupled, allocate accordingly
       ALLOCATE( U( N ),  V( N ), W( N ), &
                 MU( N ),MV( N ),MW( N ), &
                 LocalTemperature( N ),   &
                 LocalPotential( N ),     &
                 Pressure( N ),           &
                 Density( N ),            &
                 ElementNodes % x( N ),   &
                 ElementNodes % y( N ),   &
                 ElementNodes % z( N ),   &
                 ParentNodes % x( N ),   &
                 ParentNodes % y( N ),   &
                 ParentNodes % z( N ),   &
                 LocalForce( 2*N ),         &
                 TimeForce( 2*N ),         &
                 LocalMassMatrix( 2*N,2*N ),  &
                 LocalStiffMatrix( 2*N,2*N ), &
                 Load( N ), Diffusivity( 3,3,N ), &
                 SoretDiffusivity( N ),   &
                 NEConst( N ),            &
                 HeatCapacity( N ),       &
                 GasConstant( N ),        &
                 SpeciesTransferCoeff( N ), &
                 SExt( N ),              &
                 C0( N ), C1( N ), CT( N ), C2( 3,3,N ),STAT=istat )
 
       NULLIFY( HWrk)

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'AdvectionDiffusion', 'Memory allocation error.' )
       END IF

! Add parallel solution check here

       Mass = 0.0d0

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     Bubbles = .FALSE.
     Stabilize = .FALSE.
     
     StabilizeFlag = GetString( Solver % Values, 'Stabilization Method', GotIt )
     IF( GotIt ) THEN
       SELECT CASE(StabilizeFlag)
       CASE('stabilized')
         Stabilize = .TRUE.
       CASE('bubbles')
         Bubbles = .TRUE.
       CASE DEFAULT
         CALL Fatal('AdvectionDiffusion','Unknown stabilization method: '//TRIM(StabilizeFlag))
       END SELECT
     ELSE
       Stabilize = ListGetLogical( Solver % Values,'Stabilize',GotIt )
       Bubbles = ListGetLogical( Solver % Values,'Bubbles',GotIt2 )

       IF( .NOT. ( GotIt .OR. GotIt2 ) ) THEN
         CALL Info('AdvectionDiffusion','Defaulting stabilization to bubbles',Level=10)
         Bubbles = .TRUE.
       ELSE IF( Stabilize .AND. Bubbles ) THEN
         CALL Fatal('AdvectionDiffusion','Choose either Bubbles or Stabilize!')       
       END IF
     END IF
       
     NonlinearIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Max Iterations',GotIt )
     IF ( .NOT.GotIt ) NonlinearIter = 1

     Relax = GetCReal( Solver % Values, &
         'Nonlinear System Relaxation Factor',GotIt )
     IF ( .NOT.GotIt ) Relax = 1

     EquationName = ListGetString( Solver % Values, 'Equation' )

!------------------------------------------------------------------------------

     dt = Timestep

!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     DO iter=1,NonlinearIter
       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'AdvectionDiffusion', &
             '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'Iteration step: ', iter
       CALL Info( 'AdvectionDiffusion', Message, Level=4 )
       CALL Info( 'AdvectionDiffusion', &
             '-------------------------------------', Level=4 )

!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       body_id = -1
       NULLIFY(Material)
!------------------------------------------------------------------------------
!      Do the assembly for bulk elements
!------------------------------------------------------------------------------

       tmax = Solver % NumberOfActiveElements
       CALL Info( 'AdvectionDiffusion','Bulk Assembly')
 
       DO t = 1, tmax

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', &
               INT(100.0 - 100.0 * (tmax-t) / (1.0*tmax)), ' % done'           
           CALL Info( 'AdvectionDiffusion', Message, Level=5 )             
           at0 = RealTime()
         END IF
         
         CurrentElement => GetActiveElement(t)
!
!------------------------------------------------------------------------------
         IF ( CurrentElement % BodyId /= body_id ) THEN
!------------------------------------------------------------------------------
           body_id = CurrentElement % Bodyid    
           eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation')
           Eq => Model % Equations(eq_id) % Values

           ConvectionFlag = ListGetString( Eq, &
               TRIM(ComponentName(Solver % Variable)) //' Convection', GotIt )
           IF ( .NOT. GotIt ) &
             ConvectionFlag = ListGetString( Eq,'Convection', GotIt )
           
           ScaledToSolubility = .FALSE.
           ConcentrationUnits = ListGetString( Eq, 'Concentration Units', GotIt )
           IF ( .NOT.GotIt ) AbsoluteMass = .FALSE.
           IF (ConcentrationUnits == 'absolute mass') THEN
              AbsoluteMass = .TRUE.
           ELSE IF (ConcentrationUnits == 'mass to max solubility' ) THEN
              AbsoluteMass = .TRUE.
              ScaledToSolubility = .TRUE.
           ELSE
              AbsoluteMass = .FALSE.
           END IF

           k = ListGetInteger( Model % Bodies( body_id) % Values, 'Material' )
           Material => Model % Materials(k) % Values

           HeatSolName = ListGetString( Material, &
               'Temperature Field Variable', GotIt )
           IF ( Gotit ) THEN 
             TempSol => VariableGet( Mesh % Variables, &
                 TRIM( HeatSolName ) )
             IF ( ASSOCIATED( TempSol ) ) THEN
               TempPerm     => TempSol % Perm
               Temperature  => TempSol % Values
             ELSE
               WRITE( Message, * ) 'No temperature variable ' &
                   // TRIM( HeatSolName ) // ' available'
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           END IF

           PhiSolName = ListGetString( Material, &
               TRIM(ComponentName(Solver % Variable)) // &
               ' EPotential Field Variable', GotIt )
           IF ( .NOT. GotIt ) &
             PhiSolName = ListGetString( Material, &
                 'EPotential Field Variable', GotIt )
           IF ( Gotit ) THEN
             PhiSol => VariableGet( Mesh % Variables, &
                 TRIM( PhiSolName ) )
             IF ( ASSOCIATED( PhiSol ) ) THEN
               PotentialPerm => PhiSol % Perm
               EPotential    => PhiSol % Values
             ELSE
               WRITE( Message, * ) 'No electrical potential variable ' &
                   // TRIM( PhiSolName ) // ' available'
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           END IF

           ConvectName = ListGetString( Eq, &
               TRIM(ComponentName(Solver % Variable)) // &
               ' Convection Field Variable', GotIt )
           IF ( .NOT. GotIt ) &
             ConvectName = ListGetString( Material, &
                 'Convection Field Variable', GotIt )
           IF ( GotIt ) THEN
             FlowSol => VariableGet( Mesh % Variables, &
                 TRIM( ConvectName ) )
             IF ( ASSOCIATED( FlowSol ) ) THEN
               FlowPerm     => FlowSol % Perm
               NSDOFs       =  FlowSol % DOFs
               FlowSolution => FlowSol % Values
             ELSE
               WRITE( Message, * ) 'No convection  variable ' // &
                   TRIM( ConvectName ) // ' available'
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           END IF

!------------------------------------------------------------------------------
           CompressibilityFlag = ListGetString( Material, &
               'Compressibility Model', GotIt)
           IF ( .NOT.GotIt ) CompressibilityModel = Incompressible

           SELECT CASE( CompressibilityFlag )

             CASE( 'incompressible' )
             CompressibilityModel = Incompressible

             CASE( 'user defined 1' )
             CompressibilityModel = UserDefined1

             CASE( 'user defined 2' )
             CompressibilityModel = UserDefined2

             CASE( 'perfect gas equation 1' )
             CompressibilityModel = PerfectGas1

             CASE( 'perfect gas equation 2' )
             CompressibilityModel = PerfectGas2

             CASE( 'perfect gas equation 3' )
             CompressibilityModel = PerfectGas3

           CASE DEFAULT
             CompressibilityModel = Incompressible
           END SELECT
!------------------------------------------------------------------------------
         END IF
!------------------------------------------------------------------------------

         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes
 
!------------------------------------------------------------------------------
!        Get element nodal coordinates
!------------------------------------------------------------------------------
         ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
         IF ( ASSOCIATED( TempSol ) ) THEN
           LocalTemperature(1:n) = Temperature( TempPerm(NodeIndexes) )
         ELSE
           LocalTemperature(1:n) = ListGetReal(  Material, &
                'Reference Temperature',n,NodeIndexes,GotIt )
         ENDIF
         IF ( ASSOCIATED( PhiSol ) ) THEN
           LocalPotential(1:n) = EPotential( PotentialPerm(NodeIndexes) )
         ELSE
           LocalPotential(1:n) = ListGetReal(  Material, &
                'Reference EPotential',n,NodeIndexes,GotIt )
         ENDIF
!------------------------------------------------------------------------------
!        Get system parameters for element nodal points
!------------------------------------------------------------------------------

         CALL ListGetRealArray( Material,  &
             TRIM(ComponentName(Solver % Variable)) // &
                     ' Diffusivity', Hwrk, n, NodeIndexes )

         Diffusivity = 0.0d0
         IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,3
             Diffusivity( i,i,1:n ) = Hwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Hwrk,1))
             Diffusivity(i,i,1:n) = Hwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Hwrk,1))
             DO j=1,MIN(3,SIZE(Hwrk,2))
               Diffusivity( i,j,1:n ) = Hwrk(i,j,1:n)
             END DO
           END DO
         END IF

! This is also different for each species (in the same carrier gas)

         SoretDiffusivity(1:n) = ListGetReal( Material, &
              TRIM(ComponentName(Solver % Variable)) // &
              ' Soret Diffusivity', n, NodeIndexes, GotIt )
         IF ( .NOT. GotIt )  SoretDiffusivity = 0.0d0

         IonCharge = ListGetConstReal( Material, &
              TRIM(ComponentName(Solver % Variable)) // &
              ' Ion Charge', GotIt )
         IF ( .NOT. GotIt ) THEN
           NEConst   = 0._dp
           IonCharge = 0._dp
         ELSE
           NEConst(1:n) = IonCharge / LocalTemperature(1:n) * &
               ListGetConstReal( Model % Constants,'Unit Charge') / &
               ListGetConstReal( Model % Constants,'Boltzmann Constant')
         END IF

         HeatCapacity(1:n) = ListGetReal( Material,'Heat Capacity', &
             n,NodeIndexes,GotIt )

!------------------------------------------------------------------------------
!      Previous solution for element nodal points
!------------------------------------------------------------------------------
!       C(1:n) = Species( SpeciesPerm(NodeIndexes) )
!------------------------------------------------------------------------------

! We need also Density, need to go through all that jazz, if incompressible
         IF ( CompressibilityModel >= PerfectGas1 .AND. &
             CompressibilityModel <= PerfectGas3 ) THEN
!------------------------------------------------------------------------------
! Read Specific Heat Ratio
!------------------------------------------------------------------------------
           SpecificHeatRatio = ListGetConstReal( Material, &
               'Specific Heat Ratio', GotIt )
           IF ( .NOT.GotIt ) SpecificHeatRatio = 1.4d0
!------------------------------------------------------------------------------
! For an ideal gas, \gamma, c_p and R are really a constant
! GasConstant is an array only since HeatCapacity formally is
!------------------------------------------------------------------------------
           GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) * &
               HeatCapacity(1:n) / SpecificHeatRatio
         ELSE
           Density(1:n) = ListGetReal( Material,'Density',n,NodeIndexes )
         END IF
! Read p_0
!------------------------------------------------------------------------------
         ReferencePressure = 0.0_dp
         IF ( CompressibilityModel /= Incompressible ) THEN
           ReferencePressure = ListGetConstReal( Material, &
               'Reference Pressure', GotIt)
         END IF

! etc.
         Load = 0.d0
         Pressure = 0.d0

!------------------------------------------------------------------------------
!        Check for convection model
!------------------------------------------------------------------------------         
         
         IF ( ConvectionFlag == 'constant' ) THEN
           U = ListGetReal( Eq, TRIM(ComponentName(Solver % Variable)) // &
               ' Convection Velocity 1',n,NodeIndexes,GotIt )
           IF ( .NOT. GotIt ) &
               U = ListGetReal( Material,'Convection Velocity 1',n,NodeIndexes, GotIt)
           V = ListGetReal( Eq, TRIM(ComponentName(Solver % Variable)) // &
               ' Convection Velocity 2',n,NodeIndexes,GotIt )
           IF ( .NOT. GotIt ) &
             V = ListGetReal( Material,'Convection Velocity 2',n,NodeIndexes, GotIt)
           W = ListGetReal( Eq, TRIM(ComponentName(Solver % Variable)) // &
                ' Convection Velocity 3',n,NodeIndexes,GotIt )
           IF ( .NOT. GotIt ) &
             W = ListGetReal( Material,'Convection Velocity 3',n,NodeIndexes, GotIt)
         ELSE IF ( ConvectionFlag == 'computed' ) THEN
           
           IF( .NOT. ASSOCIATED( FlowSol ) ) THEN
             CALL Fatal('AdvectionDiffusion','Give > Convection Field Variable <')
           END IF
           
           DO i=1,n
             k = FlowPerm(NodeIndexes(i))
             IF ( k > 0 ) THEN

!------------------------------------------------------------------------------
               SELECT CASE( CompressibilityModel )
!------------------------------------------------------------------------------
                 CASE( PerfectGas1,PerfectGas2,PerfectGas3 )
                 Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
                 Density(i)  = Pressure(i) / &
                       ( GasConstant(i) * LocalTemperature(i) )
!------------------------------------------------------------------------------
               END SELECT
!------------------------------------------------------------------------------

               SELECT CASE( NSDOFs )
                 CASE(3)
                 U(i) = FlowSolution( NSDOFs*k-2 )
                 V(i) = FlowSolution( NSDOFs*k-1 )
                 W(i) = 0.0D0

               CASE(4)
                 U(i) = FlowSolution( NSDOFs*k-3 )
                 V(i) = FlowSolution( NSDOFs*k-2 )
                 W(i) = FlowSolution( NSDOFs*k-1 )
               END SELECT
           
             ELSE
               U(i) = 0.0d0
               V(i) = 0.0d0
               W(i) = 0.0d0
             END IF
           END DO

         ELSE  ! no convection
           U = 0.0d0
           V = 0.0d0
           W = 0.0d0
         END IF

! Set C1 to density or 1.0 if there's convection or Nernst-Planck migration
         IF ( ConvectionFlag == 'constant' .OR. ConvectionFlag == 'computed' &
              .OR. IonCharge /= 0.0d0 ) THEN
            IF (.NOT. AbsoluteMass) THEN
               C1 = Density
            ELSE
               C1 = 1.d0
            END IF
         ELSE
           C1 = 0.0d0
         END IF

         IF (.NOT. AbsoluteMass) THEN
            CT = Density
            DO i=1,3
               DO j=1,3
                  C2(i,j,:) = Density * Diffusivity(i,j,:)
               END DO
            END DO
         ELSE
            CT = 1.0d0
            DO i=1,3
               DO j=1,3
                  C2(i,j,:) = Diffusivity(i,j,:)
               END DO
            END DO
         END IF

         MU  = 0.0d0
         MV  = 0.0d0
         MW  = 0.0d0
         IF ( ASSOCIATED( MeshVelocity ) ) THEN
            DO i=1,n
              IF ( MeshPerm( NodeIndexes(i) ) > 0 ) THEN
                 MU(i) = MeshVelocity( MDOFs*(MeshPerm(NodeIndexes(i))-1)+1 )
                 MV(i) = MeshVelocity( MDOFs*(MeshPerm(NodeIndexes(i))-1)+2 )
                 IF ( MDOFs > 2 ) THEN
                    MW(i) = MeshVelocity( MDOFs*(MeshPerm(NodeIndexes(i))-1)+3 )
                 END IF
              END IF
            END DO
         END IF

         C0 = Density
         IF ( ScaledToSolubility ) THEN
            MaxSol = ListGetConstReal( Material, &
              TRIM(ComponentName(Solver % Variable)) // &
                 ' Maximum Solubility', GotIt )
            IF ( .NOT. GotIT ) THEN
               WRITE( Message, * ) 'Maximum solubility not defined in body : ', &
                    CurrentElement % BodyId
               CALL Fatal( 'AdvectionDiffusion', Message )
            END IF
            C0 = C0 / MaxSol
         END IF
             

!------------------------------------------------------------------------------
!      Get element local matrix, and rhs vector
!------------------------------------------------------------------------------
! Add body forces here, if any
         bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
           Values, 'Body Force', GotIt, 1, Model % NumberOfBodyForces )

         IF ( GotIt ) THEN
!------------------------------------------------------------------------------
!          Given species source
!------------------------------------------------------------------------------
! Not multiplied by density => Absolute mass source of c_n [kg/m3s]
!
! Take into account scaling of units if source given in physical units

           IF ( ScaledToSolubility .AND. ListGetLogical( Model % &
               BodyForces(bf_id) % Values, 'Physical Units', GotIt ) ) THEN

             Ratio = ListGetConstReal( Material, &
              TRIM(ComponentName(Solver % Variable)) // &
                 ' Maximum Solubility', GotIt )
             IF ( .NOT. GotIT ) THEN
               WRITE( Message, * ) 'Maximum solubility not defined in body : ', &
                   CurrentElement % BodyId
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           ELSE
             Ratio = 1.0d0
           END IF
           Load(1:n) = Load(1:n) + &
               ListGetReal( Model % BodyForces(bf_id) % Values,  &
               TRIM(ComponentName(Solver % Variable)) // &
                 ' Diffusion Source',n,NodeIndexes,gotIt ) / Ratio
         END IF

!------------------------------------------------------------------------------
         IF ( CoordinateSystem == Cartesian ) THEN
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveCompose( &
               LocalMassMatrix, LocalStiffMatrix, LocalForce, Load, &
               CT, C0, C1, C2, LocalTemperature, LocalPotential, U(1:n), V(1:n), W(1:n), MU, MV, MW, &
               SoretDiffusivity, NEConst, (CompressibilityModel /= Incompressible), &
               AbsoluteMass,Stabilize, Bubbles, CurrentElement, n, ElementNodes )           
         ELSE
           CALL DiffuseConvectiveGenCompose( &
               LocalMassMatrix, LocalStiffMatrix, LocalForce, Load, &
               CT, C0, C1, C2, LocalTemperature, LocalPotential, U, V, W, MU, MV, MW, &
               SoretDiffusivity, NEConst, (CompressibilityModel /= Incompressible), &
               AbsoluteMass,Stabilize,CurrentElement, n, ElementNodes )
         END IF

!------------------------------------------------------------------------------
!        If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
         TimeForce  = 0.0_dp
         IF ( TransientSimulation ) THEN
!------------------------------------------------------------------------------
!          NOTE: This will replace LocalStiffMatrix and LocalForce with the
!                combined information...
!------------------------------------------------------------------------------
           CALL Default1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
             LocalForce )
         END IF
!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
         IF ( Bubbles .AND. ( ConvectionFlag == 'computed' .OR. &
              ConvectionFlag == 'constant' ) ) THEN
           CALL Condensate( N, LocalStiffMatrix,  LocalForce, TimeForce )
         END IF

         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

!------------------------------------------------------------------------------
       END DO     !  Bulk elements
!------------------------------------------------------------------------------

       CALL DefaultFinishBulkAssembly()


!------------------------------------------------------------------------------
!     Mixed bulk - boundary element assembly
!
!     This is needed for the flux condition g_1 / g_2 = beta over a boundary
!------------------------------------------------------------------------------

       IF ( ScaledToSolubility ) THEN
         CALL Info( 'AdvectionDiffusion', 'Mixed bulk-boundary assembly')


         DO t=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           CurrentElement => Mesh % Elements(t)
           GotIt = .FALSE.
           DO i=1,Model % NumberOfBCs
             GotIt = CurrentElement % BoundaryInfo % Constraint == &
                 Model % BCs(i) % Tag 
             IF( GotIt ) EXIT
           END DO
           IF( .NOT. GotIt ) CYCLE

           BC => Model % BCs(i) % Values
           IF ( .NOT. ListGetLogical( BC, TRIM(ComponentName(Solver % Variable)) &
               // ' Solubility Change Boundary', GotIt ) )  CYCLE

!------------------------------------------------------------------------------
!             Set the current element pointer in the model structure to
!             reflect the element being processed
!------------------------------------------------------------------------------
           Model % CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
           n = CurrentElement % TYPE % NumberOfNodes
           NodeIndexes => CurrentElement % NodeIndexes
           
           ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
           ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

!------------------------------------------------------------------------------
!             Get normal target body. The parent element from the opposite
!             direction is used for normal derivative calculation
!------------------------------------------------------------------------------
           body_id = ListGetInteger( BC,'Normal Target Body', GotIt )
           
           lbody = 0
           IF ( ASSOCIATED( CurrentElement % BoundaryInfo % Left ) ) &
               lbody = CurrentElement % BoundaryInfo % Left % BodyId
           rbody = 0
           IF ( ASSOCIATED( CurrentElement % BoundaryInfo % Right ) ) &
               rbody = CurrentElement % BoundaryInfo % Right % BodyId
!------------------------------------------------------------------------------
!             If normal target body not defined check which direction is used
!             by the NormalVector function
!------------------------------------------------------------------------------

           IF ( .NOT. GotIt ) THEN
!------------------------------------------------------------------------------
!               Force normal to point to lbody
!------------------------------------------------------------------------------
             body_id = lbody
             outbody = CurrentElement % BoundaryInfo % OutBody
             CurrentElement % BoundaryInfo % OutBody = body_id
             
             Nrm = NormalVector( CurrentElement, ElementNodes, &
                 CurrentElement % TYPE % NodeU(1), CurrentElement % TYPE % NodeV(1), &
                 .TRUE. )
!------------------------------------------------------------------------------
!               Check which direction NormalVector chooses
!------------------------------------------------------------------------------
             CurrentElement % BoundaryInfo % OutBody = outbody
             Nrm2 = NormalVector( CurrentElement, ElementNodes, &
                 CurrentElement % TYPE % NodeU(1), CurrentElement % TYPE % NodeV(1), &
                 .TRUE. )
!------------------------------------------------------------------------------
!               Change body_id if Nrm and Nrm2 point to different directions
!------------------------------------------------------------------------------
             IF ( SUM( Nrm(1:3) * Nrm2(1:3) ) < 0 ) body_id = rbody
           END IF
!------------------------------------------------------------------------------

           k = ListGetInteger( Model % Bodies( body_id ) % Values, &
               'Material', minv=1, maxv=Model % NumberOfMaterials )
           Material => Model % Materials(k) % Values                
           
           Ratio = ListGetConstReal( Material, &
               TRIM(ComponentName(Solver % Variable)) // &
               ' Maximum Solubility', GotIt )
           
           IF ( .NOT. GotIT ) THEN
             WRITE( Message, * ) 'No maximum solubility defined for material : ', k
             CALL Fatal( 'AdvectionDiffusion', Message )
           END IF
           
           IF ( lbody == body_id ) THEN
             k = ListGetInteger( Model % Bodies( rbody ) % Values, 'Material', &
                 minv=1,maxv=Model % NumberOfMaterials )
             Parent => CurrentElement % BoundaryInfo % Right
           ELSE
             k = ListGetInteger( Model % Bodies(lbody) % Values, 'Material', &
                 minv=1,maxv=Model % NumberOfMaterials )
             Parent => CurrentElement % BoundaryInfo % Left
           END IF
           Material => Model % Materials(k) % Values                
           Ratio = ListGetConstReal( Material, &
               TRIM(ComponentName(Solver % Variable)) // &
               ' Maximum Solubility', GotIt ) / Ratio
           
           IF ( .NOT. GotIT ) THEN
             WRITE( Message, * ) 'No maximum solubility defined for material : ', k
             CALL Fatal( 'AdvectionDiffusion', Message )
           END IF
           Ratio = Ratio - 1.0d0
           
!------------------------------------------------------------------------------
!            Get the diffusivity tensor
!------------------------------------------------------------------------------
           CALL ListGetRealArray( Material,  &
               TRIM(ComponentName(Solver % Variable)) // &
               ' Diffusivity', Hwrk, n, NodeIndexes )
           
           Diffusivity = 0.0d0
           IF ( SIZE(Hwrk,1) == 1 ) THEN
             DO m=1,3
               Diffusivity( m,m,1:n ) = Hwrk( 1,1,1:n )
             END DO
           ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
             DO m=1,MIN(3,SIZE(Hwrk,1))
               Diffusivity(m,m,1:n) = Hwrk(m,1,1:n)
             END DO
           ELSE
             DO m=1,MIN(3,SIZE(Hwrk,1))
               DO j=1,MIN(3,SIZE(Hwrk,2))
                 Diffusivity( m,j,1:n ) = Hwrk(m,j,1:n)
               END DO
             END DO
           END IF
           
           pn = Parent % TYPE % NumberOfNodes
           
           ParentNodes % x(1:pn) = Mesh % Nodes % x(Parent % NodeIndexes)
           ParentNodes % y(1:pn) = Mesh % Nodes % y(Parent % NodeIndexes)
           ParentNodes % z(1:pn) = Mesh % Nodes % z(Parent % NodeIndexes)
!------------------------------------------------------------------------------
!             Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
           IF ( CoordinateSystem == Cartesian ) THEN
             CALL DiffuseConvectiveBBoundary( LocalStiffMatrix, Parent, pn, &
                 ParentNodes, Ratio, CurrentElement, n, ElementNodes )
           ELSE
             CALL DiffuseConvectiveGenBBoundary(LocalStiffMatrix, Parent, &
                 pn, ParentNodes, Ratio, CurrentElement,n , ElementNodes ) 
           END IF
!------------------------------------------------------------------------------
!             Update global matrices from local matrices
!------------------------------------------------------------------------------
           IF ( TransientSimulation ) THEN
             LocalMassMatrix = 0.0d0
             LocalForce = 0.0d0
             CALL Add1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
                 LocalForce,dt,pn,1,SpeciesPerm(Parent % NodeIndexes),Solver )
           END IF
           
           CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
               ForceVector, LocalForce, pn, 1, SpeciesPerm(Parent % NodeIndexes) )
!------------------------------------------------------------------------------
         END DO   ! Boundary - bulk element assembly
!------------------------------------------------------------------------------
       END IF  ! If ScaledToSolubility

!------------------------------------------------------------------------------
!     Boundary element assembly
!------------------------------------------------------------------------------
       CALL Info( 'AdvectionDiffusion','Boundary Assembly')
       
       ErrorWritten = .FALSE.
       DO t=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         CurrentElement => Mesh % Elements(t)

         DO i=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN

!------------------------------------------------------------------------------
!             Set the current element pointer in the model structure to
!             reflect the element being processed
!------------------------------------------------------------------------------
              Model % CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
              n = CurrentElement % TYPE % NumberOfNodes
              NodeIndexes => CurrentElement % NodeIndexes

              IF ( ANY( SpeciesPerm( NodeIndexes ) <= 0 ) ) CYCLE

              ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
              ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

              SpeciesTransferCoeff = 0.0D0
              SExt = 0.0D0
              Load = 0.0D0

              BC => Model % BCs(i) % Values

              SpeciesTransferCoeff(1:n) = ListGetReal( BC, &
                  'Mass Transfer Coefficient', n, NodeIndexes, GotIt )
              IF ( GotIt ) THEN
                SExt(1:n) = ListGetReal( BC, &
                    'External Concentration', n, NodeIndexes, GotIt )
                
                IF ( .NOT. AbsoluteMass .OR. ScaledToSolubility ) THEN
                  IF ( .NOT. ErrorWritten ) THEN
                    CALL Error( 'AdvectionDiffusion', '--------------------' )
                    CALL Error( 'AdvectionDiffusion', &
                        'Mass transfer coefficient possible to use only with absolute mass concentrations' )
                    CALL Error( 'AdvectionDiffusion', &
                        'Ignoring mass transfer BC' )
                    CALL Error( 'AdvectionDiffusion', '--------------------' )
                    ErrorWritten = .TRUE.
                  END IF
                  SExt = 0.0d0
                  SpeciesTransferCoeff = 0.0d0
                END IF
              ELSE
                SExt(1:n) = 0.0d0
              END IF
              
!------------------------------------------------------------------------------
!           BC: -D@c/@n = \alpha(C - Cext)
!------------------------------------------------------------------------------
              DO j=1,n
                Load(j) = Load(j) + SpeciesTransferCoeff(j) * SExt(j)
              END DO
               
!------------------------------------------------------------------------------
!             BC: j_n=-\rho*\alpha*@c/@n = g
!------------------------------------------------------------------------------
               
              IF ( ScaledToSolubility .AND. &
                  ListGetLogical( BC, 'Physical Units', GotIt ) ) THEN
                
                Ratio = ListGetConstReal( Material, &
                    TRIM(ComponentName(Solver % Variable)) // &
                    ' Maximum Solubility', GotIt )
                IF ( .NOT. GotIT ) THEN
                  WRITE( Message, * ) 'No maximum solubility defined in body : ', &
                      CurrentElement % BodyId
                  CALL Fatal( 'AdvectionDiffusion', Message )
                END IF
              ELSE
                Ratio = 1.0d0
              END IF
              Load(1:n) = Load(1:n) + &
                  ListGetReal( BC, &
                  TRIM(ComponentName(Solver % Variable)) // &
                  ' Flux', n,NodeIndexes,gotIt ) / Ratio
!------------------------------------------------------------------------------
!             Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
              IF ( CoordinateSystem == Cartesian ) THEN
                CALL DiffuseConvectiveBoundary( LocalStiffMatrix,LocalForce, &
                    Load,SpeciesTransferCoeff,CurrentElement,n,ElementNodes )
              ELSE
                CALL DiffuseConvectiveGenBoundary(LocalStiffMatrix,LocalForce,&
                    Load,SpeciesTransferCoeff,CurrentElement,n,ElementNodes ) 
              END IF
!------------------------------------------------------------------------------
!             Update global matrices from local matrices
!------------------------------------------------------------------------------
              IF ( TransientSimulation ) THEN
                LocalMassMatrix = 0.0d0
                CALL Add1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
                    LocalForce,dt,n,1,SpeciesPerm(NodeIndexes),Solver )
              END IF
              
              CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
                  ForceVector, LocalForce, n, 1, SpeciesPerm(NodeIndexes) )
!------------------------------------------------------------------------------
          END IF ! of currentelement bc == bcs(i)
        END DO ! of i=1,model bcs
      END DO   ! Boundary element assembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    FinishAssemebly must be called after all other assembly steps, but before
!    Dirichlet boundary settings. Actually no need to call it except for
!    transient simulations.
!------------------------------------------------------------------------------
      CALL DefaultFinishBoundaryAssembly()
      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
!    Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
      CALL Info( 'AdvectionDiffusion', 'Assembly done')

      at = CPUTime() - at
      st = CPUTime()

!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------

      Norm = DefaultSolve()

      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info( 'AdvectionDiffusion', Message, Level=5 )
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info( 'AdvectionDiffusion', Message, Level=5 )
!------------------------------------------------------------------------------
      RelativeChange = Solver % Variable % NonlinChange 

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'AdvectionDiffusion', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'AdvectionDiffusion', Message, Level=4 )

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT

!------------------------------------------------------------------------------
    END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Finally, check if integration of species density over volume requested.
! This is really old - don't know really what is done
!------------------------------------------------------------------------------
    IF( ListCheckPresent( Model % Simulation, 'Species Density') ) THEN        
      PreviousMass = Mass
      Mass = VolumeIntegrate( Model, Solver % ActiveElements, 'Species Density' )      
      PRINT *,'Species Mass: ',Mass
      IF ( TransientSimulation )  THEN
        PRINT *,'Mass Gain: ',Mass - PreviousMass
      END IF
    END IF

!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
!> Diffuse-convective local matrix computing (cartesian coordinates)
!>  Returns element local matrices and RSH vector for diffusion-convection
!>  equation. 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveCompose( MassMatrix,StiffMatrix,ForceVector,  &
       LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,Temperature,EPotential, &
       Ux,Uy,Uz,MUx,MUy,MUz,SoretD,NEC,Compressible,AbsoluteMass, &
       Stabilize,UseBubbles,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT:
!
!  REAL(KIND=dp) :: NodalCT,NodalC0,NodalC1
!     INPUT: Coefficient of the time derivative term, 0 degree term, and
!            the convection term respectively
!
!  REAL(KIND=dp) :: NodalC2(:,:,:)
!     INPUT: Nodal values of the diffusion term coefficient tensor
!
!
!  REAL(KIND=dp) :: Temperature
!     INPUT: Temperature from previous iteration, needed if we model
!            phase change
!
!  REAL(KIND=dp) :: EPotential
!     INPUT: Electrical potential from previous iteration
!
!  REAL(KIND=dp) :: SoretD
!     INPUT: Soret Diffusivity D_t : j_t = D_t grad(T)
!
!  REAL(KIND=dp) :: NEC
!     INPUT: Nernst-Einstein mobility const N_p : j_p = D N_p C grad(phi)
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!******************************************************************************

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Ux,Uy,Uz,MUx,MUy,MUz,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: Temperature(:), SoretD(:)
     REAL(KIND=dp) :: EPotential(:), NEC(:)
     REAL(KIND=dp) :: NodalC0(:),NodalC1(:),NodalCT(:),NodalC2(:,:,:),dT

     LOGICAL :: Stabilize,UseBubbles,Compressible,AbsoluteMass

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: Velo(3),dVelodx(3,3),Force

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,Tau,x,y,z

     REAL(KIND=dp) :: SorD, GradTemp(3), SoretForce
     REAL(KIND=dp) :: GradPhi(3), NEConst

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s,u,v,w,DivVelo

     REAL(KIND=dp) :: C0,C00,C1,CT,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)
     REAL(KIND=dp) :: NodalCThermal(n), CThermal
     REAL(KIND=dp) :: NodalCNP(n), CNP

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat,Convection,ConvectOrNPAndStabilize,Bubbles,ThermalDiffusion,NPDiffusion

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     c = dim + 1

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0
     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. Stabilize .AND. UseBubbles ) THEN
        NBasis = 2*n
        Bubbles = .TRUE.
     END IF
     
     ThermalDiffusion = .FALSE.
     IF ( ANY( ABS( SoretD(1:n) ) > AEPS ) ) THEN
        ThermalDiffusion = .TRUE. 
        NodalCThermal = NodalC0(1:n)
     END IF
     NPDiffusion = .FALSE.
     IF ( ANY( ABS( NEC(1:n) ) > AEPS ) ) THEN
        NPDiffusion = .TRUE.
        NodalCNP = NodalC0(1:n)
     END IF
     NodalC0 = 0.0d0   ! this is the way it works, maybe could do better some time

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!    If there is no convection term we don t need stabilization.
!------------------------------------------------------------------------------
     ConvectOrNPAndStabilize = .FALSE.
     IF ( Stabilize .AND. ( Convection .OR. NPDiffusion ) ) THEN
       ConvectOrNPAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
     END IF

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
             Basis,dBasisdx,ddBasisddx,ConvectOrNPAndStabilize,Bubbles )

       s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       C0 = SUM( NodalC0(1:n) * Basis(1:n) )
       C1 = SUM( NodalC1(1:n) * Basis(1:n) )
       CT = SUM( NodalCT(1:n) * Basis(1:n) )


!------------------------------------------------------------------------------
!      Coefficient of the diffusion term & it s derivatives at the
!      integration point
!------------------------------------------------------------------------------
       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO
!------------------------------------------------------------------------------
!      If there's no convection term and no Nernst-Planck migration,
!      we don't need the velocities, and
!      also no need for stabilization
!------------------------------------------------------------------------------
       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
          Convection = .TRUE.
!------------------------------------------------------------------------------
!         Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
          Velo = 0.0D0
          Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
          Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
          IF ( dim > 2 ) Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )

          IF ( Compressible .AND. AbsoluteMass ) THEN
            dVelodx = 0.0D0
            DO i=1,3
              dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
              dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
              IF ( dim > 2 ) dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
            END DO

            DivVelo = 0.0D0
            DO i=1,dim
              DivVelo = DivVelo + dVelodx(i,i)
            END DO
            C0 = DivVelo
          END IF

!------------------------------------------------------------------------------
!         For Nernst-Planck, D_ij * NEC * dphi/dx_j is like -u_i
!------------------------------------------------------------------------------
          IF ( NPDiffusion ) THEN
             NEConst = SUM ( NEC(1:n)*Basis(1:n) )
             GradPhi = 0.0d0
             DO i = 1, dim
                GradPhi(i) = SUM( dBasisdx(1:n,i) * EPotential(1:n) )
             END DO
             DO i = 1, dim
                DO j = 1, dim
                   Velo(i) = Velo(i) - C2(i,j) * NEConst * GradPhi(j)
                END DO
             END DO
          END IF

          IF ( Stabilize ) THEN
!------------------------------------------------------------------------------
!           Stabilization parameter Tau
!------------------------------------------------------------------------------
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

!#if 1
            Pe  = MIN( 1.0D0, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )

            Tau = 0.0D0
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF
!#else
!            C00 = C0
!            IF ( DT /= 0.0d0 ) C00 = C0 + CT / DT
!
!            Pe1 = 0.0d0
!            IF ( C00 /= 0.0d0 ) THEN
!              Pe1 = 2 * ABS(C2(1,1)) / ( mK * C00 * hK**2 )
!              Pe1 = C00 * hK**2 * MAX( 1.0d0, Pe1 )
!            ELSE
!              Pe1 = 2 * ABS(C2(1,1)) / mK
!            END IF
!
!            Pe2 = 0.0d0
!            IF ( C2(1,1) /= 0.0d0 ) THEN
!              Pe2 = ( mK * C1 * VNorm * hK ) / ABS(C2(1,1))
!              Pe2 = 2 * ABS(C2(1,1)) * MAX( 1.0d0, Pe2 ) / mK
!            ELSE
!              Pe2 = 2 * hK * C1 * VNorm
!            END IF
!
!            Tau = hk**2 / ( Pe1 + Pe2 )
!#endif
!------------------------------------------------------------------------------

            DO i=1,dim
              DO j=1,dim
                DO k=1,dim
                  dC2dx(i,j,k) = SUM( NodalC2(i,j,1:n)*dBasisdx(1:n,k) )
                END DO
              END DO
            END DO

!------------------------------------------------------------------------------
!           Compute residual & stablization vectors
!------------------------------------------------------------------------------
            DO p=1,N
              SU(p) = C0 * Basis(p)
              DO i = 1,dim
                SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SU(p) = SU(p) - C2(i,j) * ddBasisddx(p,i,j)
                  SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                END DO
              END DO

              SW(p) = C0 * Basis(p)
              DO i = 1,dim
                SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SW(p) = SW(p) - C2(i,j) * ddBasisddx(p,i,j)
                  SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                END DO
              END DO
            END DO
          END IF
        END IF

!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
        DO p=1,NBasis
        DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = CT * Basis(q) * Basis(p)
          A = C0 * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          DO i=1,dim
            DO j=1,dim
              A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
            END DO
          END DO

          IF ( Convection ) THEN
!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
            DO i=1,dim
              A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO
!------------------------------------------------------------------------------
!           Next we add the stabilization...
!------------------------------------------------------------------------------
            IF ( Stabilize ) THEN
              A = A + Tau * SU(q) * SW(p)
              M = M + Tau * CT * Basis(q) * SW(p)
            END IF
          END IF

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
          MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
        END DO
        END DO

!------------------------------------------------------------------------------
!       The righthand side...
!------------------------------------------------------------------------------
!       Force at the integration point
!------------------------------------------------------------------------------
        Force = SUM( LoadVector(1:n)*Basis(1:n) )

!------------------------------------------------------------------------------
        DO p=1,NBasis
          Load = Basis(p)
          IF ( ConvectOrNPAndStabilize ) Load = Load + Tau * SW(p)
          ForceVector(p) = ForceVector(p) + s * Force * Load
        END DO

!------------------------------------------------------------------------------
!     Add Soret diffusivity if necessary
!     -div( rho D_t grad(T)) 
!------------------------------------------------------------------------------

        IF ( ThermalDiffusion ) THEN

           CThermal = SUM( NodalCThermal(1:n) * Basis(1:n) )

           GradTemp = 0.0d0
           DO i = 1, dim
              GradTemp(i) = SUM( dBasisdx(1:n,i) * Temperature(1:n) )
           END DO
           SorD = SUM( Basis(1:n) * SoretD(1:n) )

           DO p=1,NBasis
              IF ( ConvectOrNPAndStabilize ) THEN
                 Load = Tau * SW(p)
              ELSE
                 Load = 1.0d0
              END IF

              SoretForce = CThermal * SorD * SUM( GradTemp(1:dim) * dBasisdx(p,1:dim) )
              ForceVector(p) = ForceVector(p) - s * SoretForce * Load
           END DO
            
        END IF

     END DO


!------------------------------------------------------------------------------
   END SUBROUTINE DiffuseConvectiveCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Returns element local matrices for a discontinuous flux boundary conditions
!>  of diffusion equation.
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveBBoundary( BoundaryMatrix, Parent, &
               pn, ParentNodes, Ratio, Element, n, Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!      OUTPUT: coefficient matrix if equations
!
!  TYPE(Element_t) :: Parent
!       INPUT: Structure describing the boundary elements parent
!
!  INTEGER :: pn
!       INPUT: Number of parent element nodes
!
!  TYPE(Nodes_t) :: ParentNodes
!       INPUT: Parent element node coordinates
!
!  REAL(KIND=dp) :: Ratio
!       INPUT: The ratio of maximal solubilities - 1 (defining the 
!             measure of discontinuity)
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number  of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, ParentNodes
     TYPE(Element_t), POINTER :: Element, Parent
     REAL(KIND=dp) :: BoundaryMatrix(:,:), Ratio
     INTEGER :: n, pn

     REAL(KIND=dp) :: ddBasisddx(n,3,3), ParentdBasisdx(pn,3)
     REAL(KIND=dp) :: Basis(n), ParentBasis(pn)
     REAL(KIND=dp) :: dBasisdx(n,3), SqrtElementMetric

     REAL(KIND=dp) :: Diff(3,3)
     REAL(KIND=dp) :: u, v, w, s, x(n), y(n), z(n), Normal(3), FluxVector(3)
     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)

     INTEGER :: ParentNodeIndexes(n)
     INTEGER :: i, t, q, p, N_Integ, j

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
           Basis,dBasisdx,ddBasisddx,.FALSE. )

       s = SqrtElementMetric * S_Integ(t)

       Normal = Normalvector( Element, Nodes, u, v, .TRUE. )
!------------------------------------------------------------------------------
!      Need parent element basis functions for calculating normal derivatives
!------------------------------------------------------------------------------
       DO i = 1,n
         DO j = 1,pn
           IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
             x(i) = Parent % TYPE % NodeU(j)
             y(i) = Parent % TYPE % NodeV(j)
             z(i) = Parent % TYPE % NodeW(j)
             ParentNodeIndexes(i) = j
             EXIT
           END IF
         END DO
       END DO

       u = SUM( Basis(1:n) * x(1:n) )
       v = SUM( Basis(1:n) * y(1:n) )
       w = SUM( Basis(1:n) * z(1:n) )

       stat = ElementInfo( Parent, ParentNodes,u, v, w, SqrtElementMetric, &
           ParentBasis, ParentdBasisdx, ddBasisddx, .FALSE. )

       FluxVector = 0.0d0

       DO i = 1, 3
         DO j = 1, 3
           Diff(i,j) = SUM( Diffusivity(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO q = 1, pn
         DO j = 1, 3
           FluxVector(j) = SUM( Diff(j,1:3) * ParentdBasisdx(q,1:3) )
         END DO
         DO i = 1, n
           p = ParentNodeIndexes(i)
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + Ratio * &
               s * Basis(i) * SUM( FluxVector(1:3) * Normal(1:3) )
         END DO
       END DO

     END DO

   END SUBROUTINE DiffuseConvectiveBBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Returns element local matrices and RSH vector for boundary conditions
!>  of diffusion convection equation.
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveBoundary( BoundaryMatrix,BoundaryVector, &
               LoadVector,NodalAlpha,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: coefficient matrix if equations
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT: coefficient of the force term
!
!  REAL(KIND=dp) :: NodalAlpha
!     INPUT: coefficient for temperature dependent term
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!   INTEGER :: n
!       INPUT: Number  of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LoadVector(:),NodalAlpha(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: u,v,w,s
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx,ddBasisddx,.FALSE. )

       s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis )
       Alpha = SUM( NodalAlpha(1:n)*Basis )

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
              s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
   END SUBROUTINE DiffuseConvectiveBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Diffuse-convective local matrix computing for general euclidian coordinates.
!>  Returns element local matrices and RSH vector for diffusion-convection
!>  equation.
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveGenCompose( MassMatrix,StiffMatrix,ForceVector, &
       LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,Temperature,EPotential, &
         Ux,Uy,Uz,MUx,MUy, MUz,SoretD,NEC, Compressible,AbsoluteMass, &
             Stabilize,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT:
!
!  REAL(KIND=dp) :: NodalCT,NodalC0,NodalC1
!     INPUT: Coefficient of the time derivative term, 0 degree term, and the
!             convection term respectively
!
!  REAL(KIND=dp) :: NodalC2(:,:,:)
!     INPUT: Nodal values of the diffusion term coefficient tensor
!

!  REAL(KIND=dp) :: Temperature
!     INPUT: Temperature from previous iteration, needed if we model
!            phase change
!
!  REAL(KIND=dp) :: EPotential
!     INPUT: Electrical potential from previous iteration
!
!  REAL(KIND=dp) :: SoretD
!     INPUT: Soret Diffusivity D_t : j_t = D_t grad(T)
!
!  REAL(KIND=dp) :: NEC
!     INPUT: Nernst-Einstein mobility const N_p : j_p = D N_p C grad(phi)
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!          used only if coefficient of the convection term (C1) is nonzero
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:) :: ForceVector,Ux,Uy,Uz,MUx,MUy,MUz,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: NodalC0(:),NodalC1(:),NodalCT(:),NodalC2(:,:,:)
     REAL(KIND=dp) :: Temperature(:), SoretD(:), dT
     REAL(KIND=dp) :: EPotential(:), NEC(:), dPhi

     LOGICAL :: Stabilize,Compressible,AbsoluteMass

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: Velo(3),Force

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,C00,Tau,Delta,x,y,z

     REAL(KIND=dp) :: SorD, GradTemp(3), SoretForce
     REAL(KIND=dp) :: GradPhi(3), NEConst

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis

     REAL(KIND=dp) :: s,u,v,w,DivVelo,dVelodx(3,3)

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     REAL(KIND=dp) :: C0,CT,C1,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)
     REAL(KIND=dp) :: NodalCThermal(n), CThermal
     REAL(KIND=dp) :: NodalCNP(n), CNP

     LOGICAL :: stat,CylindricSymmetry,Convection,ConvectOrNPAndStabilize,Bubbles
     LOGICAL :: ThermalDiffusion,NPDiffusion

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

     CylindricSymmetry = (CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                  CurrentCoordinateSystem() == AxisSymmetric)

     IF ( CylindricSymmetry ) THEN
       dim = 3
     ELSE
       dim = CoordinateSystemDimension()
     END IF
     n = element % TYPE % NumberOfNodes

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. Stabilize ) THEN
        NBasis = 2*n
        Bubbles = .TRUE.
     END IF
     
     ThermalDiffusion = .FALSE.
     IF ( ANY( ABS( SoretD(1:n) ) > AEPS ) ) THEN
        ThermalDiffusion = .TRUE. 
        NodalCThermal = NodalC0(1:n)
     END IF
     NPDiffusion = .FALSE.
     IF ( ANY( ABS( NEC(1:n) ) > AEPS ) ) THEN
        NPDiffusion = .TRUE.
        NodalCNP = NodalC0(1:n)
     END IF
     NodalC0 = 0.0d0   ! this is the way it works, maybe could do better some time

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n
 
!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!    If there is no convection term we don't need stabilization.
!------------------------------------------------------------------------------
     ConvectOrNPAndStabilize = .FALSE.
     IF ( Stabilize .AND. (ANY(NodalC1 /= 0.0D0) .OR. NPDiffusion) ) THEN
       ConvectOrNPAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
     END IF

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
             Basis,dBasisdx,ddBasisddx,ConvectOrNPAndStabilize,Bubbles )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( nodes % x(1:n)*Basis(1:n) )
         y = SUM( nodes % y(1:n)*Basis(1:n) )
         z = SUM( nodes % z(1:n)*Basis(1:n) )
       END IF

       CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )

       s = SqrtMetric * SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms at the
!      integration point
!------------------------------------------------------------------------------
       C0 = SUM( NodalC0(1:n)*Basis(1:n) )
       CT = SUM( NodalCT(1:n)*Basis(1:n) )
       C1 = SUM( NodalC1(1:n)*Basis(1:n) )
!------------------------------------------------------------------------------
!     Compute effective heatcapacity, if modelling phase change,
!     at the integration point.
!     NOTE: This is for heat equation only, not generally for diff.conv. equ.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!      Coefficient of the diffusion term & its derivatives at the
!      integration point
!------------------------------------------------------------------------------
       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SQRT(Metric(i,i)) * SQRT(Metric(j,j)) * &
                SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO
 
!------------------------------------------------------------------------------
!      If there's no convection term we don't need the velocities, and
!      also no need for stabilization
!------------------------------------------------------------------------------
       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
         Convection = .TRUE.
!------------------------------------------------------------------------------
!        Velocity and pressure (deviation) from previous iteration
!        at the integration point
!------------------------------------------------------------------------------
         Velo = 0.0D0
         Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
         Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
         IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
           Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )
         END IF

         IF ( Compressible .AND. AbsoluteMass ) THEN

           dVelodx = 0.0D0
           DO i=1,3
             dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
             dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
             IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) &
               dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
           END DO
  
           DivVelo = 0.0D0
           DO i=1,dim
             DivVelo = DivVelo + dVelodx(i,i)
           END DO
           IF ( CurrentCoordinateSystem() >= Cylindric .AND. &
                CurrentCoordinateSystem() <= AxisSymmetric ) THEN
! Cylindrical coordinates
             DivVelo = DivVelo + Velo(1)/x
           ELSE
! General coordinate system
             DO i=1,dim
               DO j=i,dim
                 DivVelo = DivVelo + Velo(j)*Symb(i,j,i)
               END DO
             END DO
           END IF
           C0 = DivVelo
         END IF

!------------------------------------------------------------------------------
!         For Nernst-Planck, D_ij * NEC * dphi/dx_j is like -u_i
!------------------------------------------------------------------------------
          IF ( NPDiffusion ) THEN
             NEConst = SUM ( NEC(1:n)*Basis(1:n) )
             GradPhi = 0.0d0
             DO i = 1, dim
                GradPhi(i) = SUM( dBasisdx(1:n,i) * EPotential(1:n) )
             END DO
             DO i = 1, dim
                DO j = 1, dim
                   Velo(i) = Velo(i) - C2(i,j) * NEConst * GradPhi(j)
                END DO
             END DO
          END IF

!------------------------------------------------------------------------------
!          Stabilization parameters...
!------------------------------------------------------------------------------
         IF ( Stabilize ) THEN
!          VNorm = SQRT( SUM(Velo(1:dim)**2) )
 
           Vnorm = 0.0D0
           DO i=1,dim
              Vnorm = Vnorm + Velo(i)*Velo(i) / Metric(i,i)
           END DO
           Vnorm = SQRT( Vnorm )
 
!#if 1
           Pe = MIN(1.0D0,mK*hK*C1*VNorm/(2*ABS(C2(1,1))))

           Tau = 0.0D0
           IF ( VNorm /= 0.0D0 ) THEN
             Tau = hK * Pe / (2 * C1 * VNorm)
           END IF
!#else
!            C00 = C0
!            IF ( dT > 0 ) C00 = C0 + CT
!
!            Pe1 = 0.0d0
!            IF ( C00 > 0 ) THEN
!              Pe1 = 2 * ABS(C2(1,1)) / ( mK * C00 * hK**2 )
!              Pe1 = C00 * hK**2 * MAX( 1.0d0, Pe1 )
!            ELSE
!              Pe1 = 2 * ABS(C2(1,1)) / mK
!            END IF
!
!            Pe2 = 0.0d0
!            IF ( C2(1,1) /= 0.0d0 ) THEN
!              Pe2 = ( mK * C1 * VNorm * hK ) / ABS(C2(1,1))
!              Pe2 = 2*ABS(C2(1,1)) * MAX( 1.0d0, Pe2 ) / mK
!            ELSE
!              Pe2 = 2 * hK * C1 * VNorm
!            END IF
!
!            Tau = hk**2 / ( Pe1 + Pe2 )
!#endif

!------------------------------------------------------------------------------
           DO i=1,dim
             DO j=1,dim
               DO k=1,3
                 dC2dx(i,j,k) = SQRT(Metric(i,i))*SQRT(Metric(j,j))* &
                      SUM(NodalC2(i,j,1:n)*dBasisdx(1:n,k))
               END DO
             END DO
           END DO
!------------------------------------------------------------------------------
!          Compute residual & stabilization weight vectors
!------------------------------------------------------------------------------
           DO p=1,n
             SU(p) = C0 * Basis(p)
             DO i = 1,dim
               SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
               IF ( Element % TYPE % BasisFunctionDegree <= 1 ) CYCLE

               DO j=1,dim
                 SU(p) = SU(p) - C2(i,j) * ddBasisddx(p,i,j)
                 SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                 DO k=1,dim
                   SU(p) = SU(p) + C2(i,j) * Symb(i,j,k) * dBasisdx(p,k)
                   SU(p) = SU(p) - C2(i,k) * Symb(k,j,j) * dBasisdx(p,i)
                   SU(p) = SU(p) - C2(k,j) * Symb(k,j,i) * dBasisdx(p,i)
                 END DO
               END DO
             END DO

             SW(p) = C0 * Basis(p)

             DO i = 1,dim
               SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
               IF ( Element % TYPE % BasisFunctionDegree <= 1 ) CYCLE

               DO j=1,dim
                 SW(p) = SW(p) - C2(i,j) * ddBasisddx(p,i,j)
                 SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                 DO k=1,dim
                   SW(p) = SW(p) + C2(i,j) * Symb(i,j,k) * dBasisdx(p,k)
                   SW(p) = SW(p) - C2(i,k) * Symb(k,j,j) * dBasisdx(p,i)
                   SW(p) = SW(p) - C2(k,j) * Symb(k,j,i) * dBasisdx(p,i)
                 END DO
               END DO
             END DO
           END DO
         END IF
       END IF
!------------------------------------------------------------------------------
!      Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,NBasis
       DO q=1,NBasis
!------------------------------------------------------------------------------
!        The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
         M = CT * Basis(q) * Basis(p)
         A = C0 * Basis(q) * Basis(p)
         DO i=1,dim
           DO j=1,dim
             A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
           END DO
         END DO

         IF ( Convection ) THEN
           DO i=1,dim
             A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
           END DO

!------------------------------------------------------------------------------
!        Next we add the stabilization...
!------------------------------------------------------------------------------
           IF ( Stabilize ) THEN
             A = A + Tau * SU(q) * SW(p)
             M = M + Tau * CT * Basis(q) * SW(p)
           END IF
         END IF

         StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
         MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
       END DO
       END DO

!------------------------------------------------------------------------------
!      Force at the integration point
!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis(1:n) )

!------------------------------------------------------------------------------
!      The righthand side...
!------------------------------------------------------------------------------
       DO p=1,NBasis
         Load = Basis(p)

         IF ( ConvectOrNPAndStabilize ) THEN
           Load = Load + Tau * SW(p)
         END IF

         ForceVector(p) = ForceVector(p) + s * Load * Force
       END DO

!------------------------------------------------------------------------------
!     Add Soret diffusivity if necessary
!     -div( rho D_t grad(T)) 
!------------------------------------------------------------------------------

        IF ( ThermalDiffusion ) THEN

           CThermal = SUM( NodalCThermal(1:n) * Basis(1:n) )

           GradTemp = 0.0d0
           IF ( CurrentCoordinateSystem() >= Cylindric .AND. &
               CurrentCoordinateSystem() <= AxisSymmetric ) THEN

             DO i = 1, dim
               GradTemp(i) = SUM( dBasisdx(1:n,i) * Temperature(1:n) )
             END DO
           ELSE
             CALL Error( 'AdvectionDiffusion', &
                 'Thermal diffusion not implemented for this coordinate system. Ignoring it.' )
           END IF

           SorD = SUM( Basis(1:n) * SoretD(1:n) )

           DO p=1,NBasis
              IF ( ConvectOrNPAndStabilize ) THEN
                 Load = Tau * SW(p)
              ELSE
                 Load = 1.0d0
              END IF

              SoretForce = CThermal * SorD * SUM( GradTemp(1:dim) * dBasisdx(p,1:dim) )
              ForceVector(p) = ForceVector(p) - s * SoretForce * Load
           END DO
            
        END IF

     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE DiffuseConvectiveGenCompose
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Return element local matrices for a discontinuous flux boundary conditions
!>  of diffusion equation in general coordinate system: 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveGenBBoundary( BoundaryMatrix, Parent, &
               pn, ParentNodes, Ratio, Element, n, Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!      OUTPUT: coefficient matrix if equations
!
!  TYPE(Element_t) :: Parent
!       INPUT: Structure describing the boundary elements parent
!
!  INTEGER :: pn
!       INPUT: Number of parent element nodes
!
!  TYPE(Nodes_t) :: ParentNodes
!       INPUT: Parent element node coordinates
!
!  REAL(KIND=dp) :: Ratio
!       INPUT: The ratio of maximal solubilities - 1 (defining the 
!             measure of discontinuity)
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number  of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!******************************************************************************

     TYPE(Nodes_t)   :: Nodes, ParentNodes
     TYPE(Element_t), POINTER :: Element, Parent
     REAL(KIND=dp) :: BoundaryMatrix(:,:), Ratio
     INTEGER :: n, pn

     REAL(KIND=dp) :: ddBasisddx(n,3,3), ParentdBasisdx(pn,3)
     REAL(KIND=dp) :: Basis(n), ParentBasis(pn)
     REAL(KIND=dp) :: dBasisdx(n,3), SqrtElementMetric

     REAL(KIND=dp) :: Diff(3,3)
     REAL(KIND=dp) :: xpos, ypos, zpos
     REAL(KIND=dp) :: u, v, w, s, x(n), y(n), z(n), Normal(3), FluxVector(3)
     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)

     INTEGER :: ParentNodeIndexes(n)
     INTEGER :: i, t, q, p, N_Integ, j

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
           Basis,dBasisdx,ddBasisddx,.FALSE. )

       s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         xpos = SUM( Nodes % x(1:n)*Basis )
         ypos = SUM( Nodes % y(1:n)*Basis )
         zpos = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( xpos,ypos,zpos )
       END IF

       Normal = Normalvector( Element, Nodes, u, v, .TRUE. )
!------------------------------------------------------------------------------
!      Need parent element basis functions for calculating normal derivatives
!------------------------------------------------------------------------------
       DO i = 1,n
         DO j = 1,pn
           IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
             x(i) = Parent % TYPE % NodeU(j)
             y(i) = Parent % TYPE % NodeV(j)
             z(i) = Parent % TYPE % NodeW(j)
             ParentNodeIndexes(i) = j
             EXIT
           END IF
         END DO
       END DO

       u = SUM( Basis(1:n) * x(1:n) )
       v = SUM( Basis(1:n) * y(1:n) )
       w = SUM( Basis(1:n) * z(1:n) )

       stat = ElementInfo( Parent, ParentNodes,u, v, w, SqrtElementMetric, &
           ParentBasis, ParentdBasisdx, ddBasisddx, .FALSE. )

       FluxVector = 0.0d0
       DO i = 1, 3
         DO j = 1, 3
           Diff(i,j) = SUM( Diffusivity(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO q = 1, pn
         DO j = 1, 3
           FluxVector(j) = SUM( Diff(j,1:3) * ParentdBasisdx(q,1:3) )
         END DO
         DO i = 1, n
           p = ParentNodeIndexes(i)
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + Ratio * &
               s * Basis(i) * SUM( FluxVector(1:3) * Normal(1:3) )
         END DO
       END DO

     END DO
     
   END SUBROUTINE DiffuseConvectiveGenBBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for boundary conditions
!>  of diffusion convection equation in general euclidian coordinates.
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveGenBoundary( BoundaryMatrix,BoundaryVector, &
              LoadVector,NodalAlpha,Element,n,Nodes)
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: coefficient matrix if equations
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT: coefficient of the force term
!
!  REAL(KIND=dp) :: NodalAlpha
!     INPUT: coefficient for temperature dependent term
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:)
     REAL(KIND=dp) :: LoadVector(:),NodalAlpha(:)
     TYPE(Nodes_t)    :: Nodes
     TYPE(Element_t),POINTER  :: Element

     INTEGER :: n
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     INTEGER :: i,t,q,p,N_Integ

     LOGICAL :: stat,CylindricSymmetry

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
 
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n
 
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx,ddBasisddx,.FALSE. )

       s =  S_Integ(t) * SqrtElementMetric
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
       END IF
!------------------------------------------------------------------------------
!      Basis function values at the integration point
!------------------------------------------------------------------------------
       Alpha = SUM( NodalAlpha(1:n)*Basis )
       Force = SUM( LoadVector(1:n)*Basis )

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
               s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
  END SUBROUTINE DiffuseConvectiveGenBoundary
!------------------------------------------------------------------------------

END SUBROUTINE AdvectionDiffusionSolver
!------------------------------------------------------------------------------

!> \}
