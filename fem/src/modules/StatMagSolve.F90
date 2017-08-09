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
! *  Authors: Juha Ruokolainen, Ville Savolainen, Jussi Heikonen,
! *           Peter R�back, Antti Pursula
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 May 2000
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  Solver for static or harmonic magnetic field with vector potential
!>  formulation. Note that this solver has a long history and all choices
!>  do not reflect the current best practices.
!
!>  Instead of this solver use MagnetoDynamics2D for 2D and cylindrically
!>  symmetric problems, and MagnetoDynamics for 3D problems.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE StatMagSolver( Model, Solver, dt, Transient )
!------------------------------------------------------------------------------

    USE DefUtils
    USE Differentials
    
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver
    REAL (KIND=DP) :: dt 
    LOGICAL :: Transient
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    TYPE(Matrix_t),POINTER :: StiffMatrix
    REAL (KIND=DP), POINTER :: ForceVector(:),MVP(:),MFD(:), PhaseAngle(:)
    REAL (KIND=DP), POINTER :: A(:), Br(:), Bz(:), Bp(:), Brim(:), Bzim(:), Babs(:), &
        Joule(:), absJoule(:), Ax(:), Ay(:), Az(:)
    
    TYPE(ValueList_t),POINTER :: Material
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t),POINTER :: CurrentElement

    TYPE(Variable_t), POINTER :: MagneticSol, BVar, TempVar

    INTEGER, POINTER :: NodeIndexes(:),MagneticPerm(:)

    REAL (KIND=DP), ALLOCATABLE :: Reluctivity(:), CurrentDensity(:), &
        LocalMassMatrix(:,:),LocalStiffMatrix(:,:),LocalForce(:), &
        Ap(:),Permeability(:), Conductivity(:), Ae(:),VecLoadVector(:,:)

#ifdef USE_ISO_C_BINDINGS
    REAL (KIND=DP) :: UNorm,RelativeChange, &
        at,at0,AngularFrequency, jc, jre, jim, &
        TotalHeating, DesiredHeating, TotalVolume, PermeabilityOfVacuum
#else
    REAL (KIND=DP) :: UNorm,RelativeChange, &
        at,at0,RealTime,CPUTime, AngularFrequency, jc, jre, jim, &
        TotalHeating, DesiredHeating, TotalVolume, PermeabilityOfVacuum
#endif

    INTEGER :: body_id, eq_id, bf_id, LocalNodes, NonlinearIter
    INTEGER :: t,n,k,istat,i,iter,q,j,dofs,dim

    LOGICAL :: GotIt, HarmonicSimulation, CalculateJouleHeating, &
        CalculateMagneticFlux, CalculateMagneticFluxAbs, AllocationsDone = .FALSE.

    SAVE LocalMassMatrix,LocalStiffMatrix,CurrentDensity,LocalForce, &
        ElementNodes,Reluctivity,AllocationsDone,Ap,MFD, &
        Permeability, Conductivity, Ae, VecLoadVector, Joule, &
        absJoule, PhaseAngle, Br, Bz, Bp, Brim, Bzim

 !------------------------------------------------------------------------------
!    Get variables needed for solving the system
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

    CALL Info( 'StatMagSolve', ' ', Level=4 )
    CALL Info( 'StatMagSolve','------------------------------------------------', Level=4 )
    CALL Info( 'StatMagSolve', 'Starting simulation', Level=4 )

    dim = Solver % Mesh % MeshDim

    PermeabilityOfVacuum = GetConstReal( Model % Constants,'Permeability Of Vacuum',GotIt)
    IF(.NOT. GotIt) PermeabilityOfVacuum = PI * 4.0d-7 

    MagneticSol => Solver % Variable 
    MagneticPerm  => MagneticSol % Perm
    MVP => MagneticSol % Values

    LocalNodes = COUNT( MagneticPerm > 0 )
    IF ( LocalNodes <= 0 ) RETURN

    StiffMatrix => Solver % Matrix
    ForceVector => StiffMatrix % RHS
    dofs = Solver % Variable % DOFs

    HarmonicSimulation = ListGetLogical( Solver % Values, &      
              'Harmonic Simulation',gotIt )
    IF (.NOT.gotIt) HarmonicSimulation = ( dofs == 2 )

    IF(HarmonicSimulation) THEN
      CALL Info('StatMagSolver','Assuming harmonic simulation')
    ELSE
      CALL Info('StatMagSolver','Assuming steady-state simulation')
    END IF

    IF( Dim == 3 ) THEN
      IF ( dofs /= 3 ) THEN
        CALL Fatal('StatMagSolver','In 3D there must be three components for the vector potential!')
        CalculateJouleHeating = .FALSE.      
      END IF
      CALL Warn('StatMagSolver','This solver does not really fulfill the Coulomb gauge, use Whitney solver in 3D!')
    ELSE
      IF(HarmonicSimulation) THEN
        IF( dofs /= 2 ) THEN
          CALL Fatal('StatMagSolver','For harmonic simulation there must be two components for the vector potential!')
        END IF
        Solver % Matrix % COMPLEX = .TRUE.
        AngularFrequency = GetAngularFrequency()
      ELSE 
        IF( dofs /= 1) THEN
          CALL Fatal('StatMagSolver','For 2d steady-state cases there should be just one component for vector potential!')
        END IF
      END IF
    END IF

    CalculateMagneticFlux = &
        ListGetLogical( Solver % Values, 'Calculate Magnetic Flux', GotIt )
    CalculateMagneticFluxAbs = &
        ListGetLogical( Solver % Values, 'Calculate Magnetic Flux Abs', GotIt )
    CalculateJouleHeating = &
        ListGetLogical( Solver % Values, 'Calculate Joule Heating', GotIt )
    IF( CalculateJouleHeating .AND. .NOT. HarmonicSimulation ) THEN
      CALL Warn('StatMagSolver','Computation of Joule heating is only relevan for harmonic systems!')
      CalculateJouleHeating = .FALSE.
    END IF


!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------

    IF (.NOT. AllocationsDone) THEN

      N = Model % MaxElementNodes

      ALLOCATE(Reluctivity(N), &
          Permeability(N), &
          Conductivity(N), &
          ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          CurrentDensity(N), &
          PhaseAngle(N), &
          LocalMassMatrix(dofs*N,dofs*N), &
          LocalStiffMatrix(dofs*N,dofs*N), &
          LocalForce(dofs*N), &
          Ap(n), &
          Ae(n), &
          VecLoadVector(3,n), & 
          STAT=istat)

!------------------------------------------------------------------------------
!    Add magnetic flux density to variables
!------------------------------------------------------------------------------        
       
      IF(CalculateMagneticFlux) THEN
        TempVar => VariableGet( Solver % Mesh % Variables,'Magnetic Flux Density')
        MFD => TempVar % Values
        Br => MFD(1::3) ! Bx
        Bz => MFD(2::3) ! By
        Bp => MFD(3::3) ! Bz

        IF(CalculateMagneticFluxAbs) THEN
          TempVar => VariableGet( Solver % Mesh % Variables,'Magnetic Flux Density_abs')
          Babs => TempVar % Values
        END IF

        IF( HarmonicSimulation ) THEN
          ALLOCATE( Brim( Model%NumberOfNodes ), Bzim( Model%NumberOfNodes ), STAT=istat )
          IF ( istat /= 0 ) CALL Fatal( 'StatMagSolve', 'Memory allocation error.' )
        END IF
      END IF

      IF ( CalculateJouleHeating )  THEN
        TempVar => VariableGet( Solver % Mesh % Variables,'Joule Heating')
        Joule => TempVar % Values

        TempVar => VariableGet( Solver % Mesh % Variables,'Joule Field')
        absJoule => TempVar % Values
      END IF

      AllocationsDone = .TRUE.
    END IF

!--------------------------------------------------------------

    NonlinearIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Max Iterations', GotIt )
    IF(.NOT. GotIt) NonlinearIter = 1

!---------------------------------------------------------------

    DO iter=1,NonlinearIter
  
       WRITE( Message, '(A,I0)' ) 'Magnetic Field Iteration: ', iter
       CALL Info( 'StatMagSolve', Message, Level=4)
      
       at  = CPUTime()
       at0 = RealTime()

       CALL DefaultInitialize()
       
       DO t=1, Solver % NumberOfActiveElements 
         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
               (Model % NumberOfBulkElements-t) / &
               (1.0*Model % NumberOfBulkElements)), ' % done'
           
           CALL Info( 'StatMagSolve', Message, Level=5 )
           at0 = RealTime()
         END IF
     
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where the equations
!        should be calculated
!------------------------------------------------------------------------------

         CurrentElement => GetActiveElement(t)
         NodeIndexes => CurrentElement % NodeIndexes
         n  = GetElementNOFNodes()

         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))
         
         body_id = CurrentElement % BodyId
         k = ListGetInteger( Model % Bodies(body_id) % Values, 'Material', &
                  minv=1, maxv=Model % NumberOFMaterials )
         Material => Model % Materials(k) % Values
         
         Permeability(1:n) = ListGetReal(Material, &
             'Relative Permeability',n,NodeIndexes,GotIt)
         IF( GotIt ) THEN
	   Permeability(1:n) = PermeabilityOfVacuum * Permeability(1:n)
         ELSE
           Permeability(1:n) = ListGetReal(Material, &
             'Magnetic Permeability',n,NodeIndexes,GotIt)
           IF(.NOT. GotIt) Permeability = PermeabilityOfVacuum
         END IF

         Reluctivity(1:n) = 1.0 / Permeability(1:n)
         
         IF(HarmonicSimulation) THEN
           Conductivity(1:n) = ListGetReal(Material, &
               'Electrical Conductivity',n,NodeIndexes,GotIt)
           IF( GotIt ) THEN
             CALL Warn('StatMagSolve','Use > Electric Conductivity < instead of electrial')
           ELSE
             Conductivity(1:n) = ListGetReal(Material, &
                 'Electric Conductivity',n,NodeIndexes)
           END IF
         END IF
         
!------------------------------------------------------------------------------
!        Set body forces (applied current densities)
!------------------------------------------------------------------------------
  
         bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
             'Body Force',gotIt, minv=1, maxv=Model % NumberOFBodyForces )


         IF( dim < 3) THEN
           IF ( bf_id > 0  ) THEN
             CurrentDensity(1:n) = ListGetReal( &
                 Model % BodyForces(bf_id) % Values,'Current Density',n,NodeIndexes,GotIt )
           ELSE 
             CurrentDensity(1:n) = 0.0d0
           END IF
           
           IF(HarmonicSimulation) THEN
             IF(bf_id > 0) THEN
               PhaseAngle(1:n) = ListGetReal( &
                   Model % BodyForces(bf_id) % Values,'Current Phase Angle',n,NodeIndexes,GotIt )
             ELSE
               PhaseAngle(1:n) = 0.0d0
             END IF
           END IF
         END IF

         IF ( dim == 3) THEN
           VecLoadVector=0.0_dp
           IF ( bf_id > 0  ) THEN
             VecLoadVector(1,1:n) = VecLoadVector(1,1:n) + ListGetReal( &
                 Model % BodyForces(bf_id) % Values, &
                 'Current Density 1',n,NodeIndexes,gotIt )
             
             VecLoadVector(2,1:n) = VecLoadVector(2,1:n) + ListGetReal( &
                 Model % BodyForces(bf_id) % Values, &
                 'Current Density 2',n,NodeIndexes,gotIt )
             
             VecLoadVector(3,1:n) = VecLoadVector(3,1:n) + ListGetReal( &
                 Model % BodyForces(bf_id) % Values, &
                 'Current Density 3',n,NodeIndexes,gotIt )
           END IF
         END IF

!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
       IF ( dim < 3 ) THEN
         IF(HarmonicSimulation) THEN
           CALL HarmMagAxisCompose( &
               LocalStiffMatrix,LocalForce,CurrentDensity,PhaseAngle,Reluctivity, &
               Conductivity,AngularFrequency,CurrentElement,n,ElementNodes )
         ELSE
           CALL StatMagAxisCompose( &
               LocalMassMatrix,LocalStiffMatrix,LocalForce, &
               CurrentDensity,Reluctivity,Ap,CurrentElement,n,ElementNodes )
         END IF
       ELSE
         ! Note: this formulation is far from general
          CALL StatMagCartesianCompose( &
               LocalStiffMatrix, LocalForce, VecLoadVector, Reluctivity, &
               CurrentElement, n, ElementNodes )
       END IF

       CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
       
     END DO
      
     CALL Info( 'StatMagSolve', 'Assembly done', Level=4 )

     CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
     CALL DefaultDirichletBCs()
     
     CALL Info( 'StatMagSolve', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
     UNorm = DefaultSolve()

!---------------------------------------------------------------------
! Compute the magnetic flux density from the potential: B = curl A
! Here A = A_phi e_phi and B = B_rho e_rho + B_z e_z due to axisymmetry
!---------------------------------------------------------------------

     IF(CalculateMagneticFlux) THEN
       
       MFD = 0.0d0
       IF(HarmonicSimulation) THEN
         A => MVP(2::2)     
         CALL AxiSCurl(Bp,Bp,A,Brim,Bzim,Bp,MagneticPerm)
         
         A => MVP(1::2)
         CALL AxiSCurl(Bp,Bp,A,Br,Bz,Bp,MagneticPerm)
         
         DO i=1, Model%NumberofNodes
           j = MagneticPerm(i)
           IF(j > 0) THEN
             Br(j) = SQRT( Br(j)**2 + Brim(j)**2 )
             Bz(j) = SQRT( Bz(j)**2 + Bzim(j)**2 )
             Brim(j) = 0.0d0
             Bzim(j) = 0.0d0
           END IF
         END DO
       ELSE
         IF ( dim < 3 ) THEN
           A => MVP
           CALL AxiSCurl(Bp,Bp,A,Br,Bz,Bp,MagneticPerm)
         ELSE
           Ax => MVP(1::3)
           Ay => MVP(2::3)
           Az => MVP(3::3)
           
           Br = 0.0d0
           Bp = 0.0d0
           Bz = 0.0d0
           CALL Curl( Ax, Ay, Az, Br, Bz, Bp, MagneticPerm )
         END IF
       END IF
       
       CALL InvalidateVariable( Model % Meshes, Solver % Mesh, &
         'Magnetic Flux Density')
              
       IF( CalculateMagneticFluxAbs ) THEN
         DO i=1,SIZE(Babs)
           Babs(i) = SQRT( Br(i)**2 + Bz(i)**2 + Bp(i)**2 )
         END DO

         CALL InvalidateVariable( Model % Meshes, Solver % Mesh, &
             'Magnetic Flux Density_abs')
       END IF
     END IF
     

     IF( Solver % Variable % NonlinConverged == 1 ) THEN
       WRITE( Message,'(A,I0,A)' ) 'Convergence after ',iter,' iterations'
       CALL Info( 'StatMagSolve', Message, Level = 4 )
       EXIT
     END IF

   END DO


!---------------------------------------------------------------------
! In case of harmonic simulation it is possible to compute the Joule losses
!---------------------------------------------------------------------

   IF(HarmonicSimulation .AND. CalculateJouleHeating) THEN
    
     jc = 0.5_dp * AngularFrequency**2
     
     DO i=1, Model % NumberofNodes
       j = MagneticPerm(i) 
       IF(j > 0) THEN
         jre = MVP(2*j-1)
         jim = MVP(2*j)
         absJoule(j) = jc * (jre*jre+jim*jim)
       END IF
     END DO
     
     TotalHeating = 0.0d0
     TotalVolume = 0.0d0

    
     DO t=1, Solver % NumberOfActiveElements
       
       CurrentElement => GetActiveElement(t)    
       NodeIndexes => CurrentElement % NodeIndexes
       n  = GetElementNOFNodes()
       
       ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)
       
       body_id = CurrentElement % BodyId
       k = ListGetInteger( Model % Bodies(body_id) % Values, 'Material', &
               minv=1, maxv=Model % NumberOFMaterials )
       Material => Model % Materials(k) % Values
       
       Conductivity(1:n) = ListGetReal(Material, &
           'Electrical Conductivity',n,NodeIndexes,GotIt)
       IF( GotIt ) THEN
         CALL Warn('StatMagSolve','Use electric conductivity instead of electrical')
       ELSE
         Conductivity(1:n) = ListGetReal(Material, &
             'Electric Conductivity',n,NodeIndexes)
       END IF
       
       Ae(1:n) = absJoule(MagneticPerm(NodeIndexes(1:n)))
       CALL JouleIntegrate(Ae,Conductivity,TotalHeating,TotalVolume,&
           CurrentElement,n,ElementNodes )
       
       DO i=1,n
         j = MagneticPerm(NodeIndexes(i))
         IF(j > 0) THEN
           jc = Conductivity(i) * absJoule(j)
           IF(jc > Joule(j)) Joule(j) = jc
         END IF
       END DO
       
     END DO
   
     DesiredHeating = ListGetConstReal( Solver % Values, 'Desired Heating Power',gotIt)
     IF(.NOT. GotIt) DesiredHeating = ListGetConstReal( Solver % Values, 'Power Control',gotIt)
     IF(gotIt .AND. TotalHeating > 0.0d0) THEN
       absJoule = (DesiredHeating/TotalHeating) * absJoule
       Joule = (DesiredHeating/TotalHeating) * Joule
     END IF
     
     WRITE(Message,'(A,ES15.4)') 'Joule Heating (W): ',TotalHeating
     CALL Info('StatMagSolve',Message,Level=4)
     CALL ListAddConstReal( Model % Simulation, 'res: Joule heating',TotalHeating)
   END IF

   CALL Info( 'StatMagSolve', 'All done, exiting',Level=4)
   CALL Info( 'StatMagSolve','------------------------------------------------',Level=4)


CONTAINS

!------------------------------------------------------------------------------
!> Subroutine for computing local matrices for static magnetic field
!> in cylindrical coordinates with axisymmetry.
!  Author:       Jussi Heikonen
!------------------------------------------------------------------------------
  SUBROUTINE StatMagAxisCompose( &
      MassMatrix,StiffMatrix,ForceVector,LoadVector,NodalReluctivity, &
      Ap,Element,n,Nodes)

!------------------------------------------------------------------------------
!
!  REAL (KIND=DP) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL (KIND=DP) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL (KIND=DP) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL (KIND=DP) :: LoadVector(:)
!     INPUT:
!
!  REAL (KIND=DP) :: NodalReluctivity(:)
!     INPUT: Nodal values of relucitivity ( 1 / permeability )
!
!  REAL (KIND=DP) :: Ap(:)
!     INPUT: Vector potential from the previous iteration for computing
!            the "reluctivity" in a nonlinear material

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
      USE Types
      USE Integration
      USE ElementDescription

      IMPLICIT NONE
     
      REAL (KIND=DP),TARGET :: MassMatrix(:,:),StiffMatrix(:,:),&
          ForceVector(:)
      REAL (KIND=DP) :: NodalReluctivity(:), Reluctivity
      REAL (KIND=DP) :: LoadVector(:),Ap(:)

      INTEGER :: n

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

      REAL (KIND=DP) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
      REAL (KIND=DP) :: SqrtElementMetric
      REAL (KIND=DP) :: Force,r,Br,Bz,Babs,mat
      REAL (KIND=DP), POINTER :: A(:,:),M(:,:),Load(:)
      INTEGER :: DIM,t,i,j,p,q
      REAL (KIND=DP) :: s,u,v,w
      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
      INTEGER :: N_Integ
      REAL (KIND=DP), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,&
          S_Integ
      LOGICAL :: stat
!------------------------------------------------------------------------------
	  
     DIM = 2

     ForceVector = 0.0D0
     MassMatrix  = 0.0D0
     StiffMatrix = 0.0D0

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
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx,ddBasisddx,.FALSE. )

       r = SUM( Basis(1:n) * Nodes%x(1:n))
       s = SqrtElementMetric * S_Integ(t)
       
!------------------------------------------------------------------------------
!     Values at integration point
!------------------------------------------------------------------------------

       Force = SUM( LoadVector(1:n) * Basis(1:n) )
       Reluctivity = SUM( NodalReluctivity(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
       DO p=1,N
         DO q=1,N
!------------------------------------------------------------------------------
!      The equation for the vector potential
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!      Mass matrix:
!------------------------------------------------------------------------------
!           MassMatrix(p,q) = MassMatrix(p,q) + Basis(p)*Basis(q)*s*r

!------------------------------------------------------------------------------
!      Stiffness matrix:
!------------------------------

           mat = r*dBasisdx(p,1)*dBasisdx(q,1) + &
               r*dBasisdx(p,2)*dBasisdx(q,2) + &
               Basis(p)*dBasisdx(q,1) + &
               Basis(q)*dBasisdx(p,1) + &
               Basis(p)*Basis(q)/r 
           mat = mat * Reluctivity * s

           StiffMatrix(p,q) = StiffMatrix(p,q) + mat

         END DO
       END DO

!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------
       DO p=1,N
         ForceVector(p) = ForceVector(p) + Force*Basis(p)*r*s
       END DO

     END DO

   END SUBROUTINE StatMagAxisCompose


!------------------------------------------------------------------------------
!> Subroutine for computing local matrices for harmonic magnetic field
!> in cylindrical coordinates with axisymmetry.
! Author: Peter R�back
!------------------------------------------------------------------------------
   SUBROUTINE HarmMagAxisCompose( &
       StiffMatrix,ForceVector,CurrentDensity,NodalAngle,NodalReluctivity, NodalConductivity,&
       Wang,Element,n,Nodes)

!------------------------------------------------------------------------------
     USE Types
     USE Integration
     USE ElementDescription
     
     IMPLICIT NONE
     
     REAL (KIND=DP),TARGET :: StiffMatrix(:,:), ForceVector(:)
     REAL (KIND=DP) :: NodalReluctivity(:), NodalAngle(:), Reluctivity, &
         NodalConductivity(:), Conductivity, Angle
     REAL (KIND=DP) :: CurrentDensity(:),Wang
     
     INTEGER :: n
     
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL (KIND=DP) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
     REAL (KIND=DP) :: SqrtElementMetric     
     REAL (KIND=DP) :: Force,r,a11,a21,a12,a22
     REAL (KIND=DP), POINTER :: A(:,:),M(:,:),Load(:)
     
     INTEGER :: DIM,t,i,j,p,q     
     REAL (KIND=DP) :: s,u,v,w
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ
     REAL (KIND=DP), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,&
         S_Integ
     
     LOGICAL :: stat

!------------------------------------------------------------------------------

     DIM = 2

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0

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
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx,ddBasisddx,.FALSE. )

       r = SUM( basis*nodes%x(1:n))
       s = SqrtElementMetric * S_Integ(t)
       
!------------------------------------------------------------------------------
!     Force at integration point
!------------------------------------------------------------------------------
       Force = 0.0D0
       Force = SUM( CurrentDensity(1:n)*Basis )

       Reluctivity = SUM( NodalReluctivity(1:n)*Basis(1:n) )
       Conductivity = SUM( NodalConductivity(1:n)*Basis(1:n) )
       Angle = (PI/180.0d0) * SUM( NodalAngle(1:n)*Basis(1:n) )

!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
       DO p=1,N
         DO q=1,N
!------------------------------------------------------------------------------
!      The equation for the vector potential
!      Stiffness matrix:
!------------------------------

           a11 = dBasisdx(p,1)*dBasisdx(q,1)*r + &
               dBasisdx(p,2)*dBasisdx(q,2)*r + &
               Basis(p)*dBasisdx(q,1) + &
               Basis(q)*dBasisdx(p,1) + &
               Basis(p)*Basis(q)/r 
           a11 = Reluctivity * s * a11
           a22 = a11

           a21 = -Conductivity * wang * s * r * Basis(q) * Basis(p) 
           a12 = -a21

           StiffMatrix(2*p-1,2*q-1) = StiffMatrix(2*p-1,2*q-1) + a11
           StiffMatrix(2*p,2*q)     = StiffMatrix(2*p,2*q) + a22

           StiffMatrix(2*p,2*q-1) = StiffMatrix(2*p,2*q-1) + a21
           StiffMatrix(2*p-1,2*q) = StiffMatrix(2*p-1,2*q) + a12

         END DO
       END DO

!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------
       DO p=1,N

         ForceVector(2*p-1) = ForceVector(2*p-1) + Force * COS(Angle) * Basis(p) * r * s
         ForceVector(2*p)   = ForceVector(2*p) + Force * SIN(Angle) * Basis(p) * r * s         

       END DO

     END DO

   END SUBROUTINE HarmMagAxisCompose



!------------------------------------------------------------------------------
   SUBROUTINE JouleIntegrate( &
       NodalField,NodalConductivity,TotalHeating,TotalVolume,Element,n,Nodes)

!------------------------------------------------------------------------------

     USE Types
     USE Integration
     USE ElementDescription
     
     IMPLICIT NONE
     
     REAL (KIND=DP) :: NodalConductivity(:), NodalField(:)
     REAL (KIND=DP) :: TotalHeating, TotalVolume, Conductivity, Field
     INTEGER :: n
     
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL (KIND=DP) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
     REAL (KIND=DP) :: SqrtElementMetric     
     INTEGER :: DIM,t,i,j,p,q     
     REAL (KIND=DP) :: r,s,u,v,w
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ
     REAL (KIND=DP), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     
     LOGICAL :: stat

!------------------------------------------------------------------------------

     DIM = 2

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
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx,ddBasisddx,.FALSE. )

       r = SUM( basis*nodes%x(1:n))
       s = SqrtElementMetric * S_Integ(t)
       
!------------------------------------------------------------------------------
!     Force at integration point
!------------------------------------------------------------------------------

       Field = SUM( NodalField(1:n)*Basis(1:n) )
       Conductivity = SUM( NodalConductivity(1:n)*Basis(1:n) )

       DO p=1,N
         TotalVolume = TotalVolume + 2.0d0 * PI * r * s * Basis(p) 
         TotalHeating = TotalHeating + 2.0d0 * PI * r * s * Basis(p) * Field * Conductivity
       END DO

     END DO

   END SUBROUTINE JouleIntegrate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for computing local matrices for static magnetic field
!> in 3D cartesian coordinates. Makes some heavy assumptions.
! Author:       Antti Pursula
!------------------------------------------------------------------------------
    SUBROUTINE StatMagCartesianCompose( &
        StiffMatrix,ForceVector,LoadVector,NodalReluctivity, &
        Element,n,Nodes)

!------------------------------------------------------------------------------
!
!
!  REAL (KIND=DP) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL (KIND=DP) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL (KIND=DP) :: LoadVector(:)
!     INPUT: The source current density
!
!  REAL (KIND=DP) :: NodalReluctivity(:)
!     INPUT: Nodal values of relucitivity ( 1 / permeability )
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
      USE Types
      USE Integration
      USE ElementDescription

      IMPLICIT NONE
     
      REAL (KIND=DP),TARGET :: StiffMatrix(:,:), ForceVector(:)
      REAL (KIND=DP) :: NodalReluctivity(:)
      REAL (KIND=DP) :: LoadVector(:,:)

      INTEGER :: n

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
      REAL (KIND=DP) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
      REAL (KIND=DP) :: SqrtElementMetric
      REAL (KIND=DP) :: Force(3), mat, Reluctivity
      REAL (KIND=DP) :: s,u,v,w
      REAL (KIND=DP), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
      INTEGER :: N_Integ
      INTEGER :: t,p,q
      LOGICAL :: stat

!------------------------------------------------------------------------------

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0

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
     DO t = 1, N_Integ

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
!     Values at integration point
!------------------------------------------------------------------------------

       DO i = 1, 3
          Force(i) = SUM( LoadVector(i,1:n) * Basis(1:n) )
       END DO

       Reluctivity = SUM( NodalReluctivity(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
       DO p=1,N
         DO q=1,N
!------------------------------------------------------------------------------
!      The equation for the vector potential
!------------------------------------------------------------------------------
!      Stiffness matrix:
!------------------------------

           mat = SUM( dBasisdx(p,1:3)*dBasisdx(q,1:3) )
           mat = mat * Reluctivity * s

           StiffMatrix(3*p-2,3*q-2) = StiffMatrix(3*p-2,3*q-2) + mat
           StiffMatrix(3*p-1,3*q-1) = StiffMatrix(3*p-1,3*q-1) + mat
           StiffMatrix(3*p,3*q) = StiffMatrix(3*p,3*q) + mat

         END DO
       END DO

!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------
       DO p=1,N

         ForceVector(3*p-2) = ForceVector(3*p-2) + Force(1)*Basis(p)*s
         ForceVector(3*p-1) = ForceVector(3*p-1) + Force(2)*Basis(p)*s
         ForceVector(3*p) = ForceVector(3*p) + Force(3)*Basis(p)*s

      END DO

     END DO

   END SUBROUTINE StatMagCartesianCompose
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE StatMagSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization for the primary solver: StatMagSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE StatMagSolver_Init( Model, Solver, dt, Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found, Calculate, CalculateAbs
    INTEGER :: Dim,GivenDim
!------------------------------------------------------------------------------
    Params => GetSolverParams()
    Dim = CoordinateSystemDimension()

    IF( GetLogical( Params,'Harmonic Simulation',Found ) ) THEN
      CALL ListAddInteger( Params,'Variable Dofs',2 )
    END IF

    Calculate = ListGetLogical(Params,'Calculate Magnetic Flux',Found)
    CalculateAbs = ListGetLogical( Params, 'Calculate Magnetic Flux Abs',Found)
  
    IF( CalculateAbs .AND. .NOT. Calculate ) THEN
      CALL Warn('StatMagSolver_init','Cannot compute Abs without computing field')
    END IF

    IF( Calculate ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs 3 Magnetic Flux Density' )
      IF( CalculateAbs ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
            'Magnetic Flux Density_abs' )	
      END IF
    END IF
    Calculate = ListGetLogical( Params, 'Calculate Joule Heating', Found )
    IF( Calculate ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          'Joule Heating' )
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          'Joule Field' )
    END IF

  END SUBROUTINE StatMagSolver_Init
