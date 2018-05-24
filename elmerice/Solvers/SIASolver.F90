!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 23 Jully 2009 
! * 
! *****************************************************************************
!> Solver to inquire the velocity and isotropic pressure from the SIA solution
!> Exported Variables SIAFlow, dof=dim
!> Needs to first compute the Depth and the FreeSurfGrad using FlowDepth Solver 
! *****************************************************************************
SUBROUTINE SIASolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the SIA Flow solution !
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
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
  TYPE(Element_t),POINTER :: CurrentElement, Element
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, Grad1Sol, Grad2Sol, &
                               DepthSol, VeloSol

  LOGICAL :: AllocationsDone = .FALSE., Found, LimitVelocity,UnFoundFatal=.TRUE.
  
  INTEGER :: i, j, n, m, t, istat, DIM, p, Indexes(128), COMP, constrainedvelocities
  INTEGER, POINTER :: Permutation(:), VeloPerm(:), &
       DepthPerm(:), GradSurface1Perm(:), GradSurface2Perm(:), &
       NodeIndexes(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Depth(:), GradSurface1(:), &
                            GradSurface2(:), Velocity(:), PrevVelo(:,:)
  REAL(KIND=dp) :: Norm, cn, dd 

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalGravity(:), NodalViscosity(:), NodalDensity(:), &
           NodalDepth(:), NodalSurfGrad1(:), NodalSurfGrad2(:), &
           NodalU(:), NodalV(:)
  REAL(KIND=dp) :: VelocityCutOff, VeloAbs, VelocityDirection(2) 

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
       

  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName
  SAVE NodalGravity, NodalViscosity, NodalDensity, &
           NodalDepth, NodalSurfGrad1, NodalSurfGrad2, &
           NodalU, NodalV
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values


  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     WRITE(SolverName, '(A)') 'SIASolver'
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalGravity, &
                       NodalViscosity, NodalDensity, NodalDepth, &
                       NodalSurfGrad1, NodalSurfGrad2, NodalU, NodalV )

     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), &
          NodalGravity(N), NodalDensity(N), NodalViscosity(N), &
          NodalDepth(N), NodalSurfGrad1(N), NodalSurfGrad2(N), &
          NodalU(N), NodalV(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF


     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()

  VeloSol => VariableGet( Solver % Mesh % Variables, 'SIAFlow',UnFoundFatal=UnFoundFatal)
  Velocity => VeloSol % Values
  VeloPerm => VeloSol % Perm
  PrevVelo => veloSol % PrevValues
  
  DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
  Depth => DepthSol % Values
  DepthPerm => DepthSol % Perm
  
  Grad1Sol => VariableGet( Solver % Mesh % Variables, 'FreeSurfGrad1',UnFoundFatal=UnFoundFatal)
  GradSurface1 => Grad1Sol % Values
  GradSurface1Perm => Grad1Sol % Perm

  IF (dim > 2) THEN
     Grad2Sol => VariableGet( Solver % Mesh % Variables, 'FreeSurfGrad2',UnFoundFatal=UnFoundFatal)
     GradSurface2 => Grad2Sol % Values
     GradSurface2Perm => Grad2Sol % Perm
  END IF

  VelocityCutOff = ListGetConstReal( Solver % Values, 'Velocity Cutoff', Found)
  IF (.NOT.Found) THEN
     LimitVelocity = .FALSE.
     CALL INFO(SolverName,'No Velocity Cutoff has been set',Level=1)
  ELSE
     LimitVelocity = .TRUE.
     WRITE(Message, '(A,E10.4)') 'Velocity Cutoff set to ', VelocityCutOff
     CALL INFO(SolverName,Message, Level=1)
  END IF


     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS


  ! Loop over the velocity components and pressure 
  ! If DIM = 2 u, w, p
  ! If DIM = 3 u, v, w, p
  !-----------------------------------------------
  DO  COMP = 1, DIM+1

! No non-linear iteration, no time dependency  
  VariableValues = 0.0d0
  Norm = Solver % Variable % Norm



  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes

     ! Read the gravity in the Body Force Section 
     BodyForce => GetBodyForce()
     NodalGravity = 0.0_dp
     IF ( ASSOCIATED( BodyForce ) ) THEN
           IF (DIM==2) THEN 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
           ELSE 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
           END IF
     END IF
     
     ! Read the Viscosity eta, density, and exponent m in Material Section
     ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
     Material => GetMaterial()

     NodalDensity = 0.0D0
     NodalDensity(1:n) = ListGetReal( &
         Material, 'Density', n, NodeIndexes, Found )

     NodalViscosity = 0.0D0
     NodalViscosity(1:n) = ListGetReal( &
         Material, 'Viscosity', n, NodeIndexes, Found )
     
     cn = ListGetConstReal( Material, 'Viscosity Exponent',Found)

     ! Get the Nodal value of Depth, FreeSurfGrad1 and FreeSurfGrad2
     NodalDepth(1:n) = Depth(DepthPerm(NodeIndexes(1:n)))
     NodalSurfGrad1(1:n) = GradSurface1(GradSurface1Perm(NodeIndexes(1:n)))
     NodalSurfGrad2 = 0.0D0
     IF (DIM==3) NodalSurfGrad2(1:n) = GradSurface2(GradSurface2Perm(NodeIndexes(1:n)))

     IF (COMP==1) THEN     ! u
        CALL LocalMatrixUV (  STIFF, FORCE, Element, n, NodalGravity, &
        NodalDensity, NodalViscosity, NodalDepth, NodalSurfGrad1, &
        NodalSurfGrad2, cn, COMP)

     ELSE IF (COMP==DIM) THEN  ! w
        NodalU(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+1)
        NodalV = 0.0D0
        IF (DIM==3) NodalV(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+2)
        CALL LocalMatrixW (  STIFF, FORCE, Element, n, NodalU, NodalV ) 

     ELSE IF (COMP==DIM+1) THEN ! p
        CALL LocalMatrixP (  STIFF, FORCE, Element, n )

     ELSE               ! v if dim=3
        CALL LocalMatrixUV (  STIFF, FORCE, Element, n, NodalGravity, &
        NodalDensity, NodalViscosity, NodalDepth, NodalSurfGrad1, &
        NodalSurfGrad2, cn, COMP)

     END IF

     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  ! Neumann conditions only for w and p
  IF (COMP .GE. DIM) THEN
  DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     Element => GetBoundaryElement(t)
     IF ( GetElementFamily() == 1 ) CYCLE
     NodeIndexes => Element % NodeIndexes
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()
     STIFF = 0.0D00
     FORCE = 0.0D00

     IF (COMP==DIM) THEN
     ! only for the surface nodes
        dd = SUM(ABS(Depth(Depthperm(NodeIndexes(1:n)))))
        IF (dd < 1.0e-6) THEN
           NodalU(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+1)
           NodalV = 0.0D0
           IF (DIM==3) NodalV(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+2)
           CALL LocalMatrixBCW (  STIFF, FORCE, Element, n, NodalU, NodalV ) 
        END IF
     ELSE IF (COMP==DIM+1) THEN
            CALL LocalMatrixBCP(  STIFF, FORCE, Element, n, NodalDensity, &
                    NodalGravity )
     END IF
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  END IF

  CALL DefaultFinishAssembly()

  ! Dirichlet 
     CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          ComponentName('SIAFlow',COMP), 1,1, Permutation )

  CALL EnforceDirichletConditions( Solver, Solver % Matrix , Solver % Matrix % RHS)

  !Solve the system
  Norm = DefaultSolve()

  ! Save the solution on to the right variable
  constrainedvelocities = 0
  DO i = 1, Model % Mesh % NumberOfNodes
     IF (VeloPerm(i)>0) THEN
        Velocity ((DIM+1)*(VeloPerm(i)-1) + COMP) = VariableValues(Permutation(i)) 
        IF (LimitVelocity .AND. (COMP==(DIM-1))) THEN
           VeloAbs = (Velocity ((DIM+1)*(VeloPerm(i)-1) + COMP))**2
           IF (COMP==2) VeloAbs = VeloAbs + (Velocity ((DIM+1)*(VeloPerm(i)-1) + COMP - 1))**2
           VeloAbs = SQRT(VeloAbs)
           IF (VeloAbs > VelocityCutOff) THEN                    
              constrainedvelocities = constrainedvelocities + 1
              DO j= 1,DIM-1
                 Velocity ((DIM+1)*(VeloPerm(i)-1) + j)  = VelocityCutOff * Velocity ((DIM+1)*(VeloPerm(i)-1) + j) / VeloAbs
              END DO
           END IF
              
        END IF
     END IF
  END DO

END DO ! Loop p

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUV(  STIFF, FORCE, Element, n, gravity, &
           Density, Viscosity, LocalDepth, SurfGrad1, SurfGrad2, cm, cp)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
                     Viscosity(:), LocalDepth(:), SurfGrad1(:), SurfGrad2(:)
    INTEGER :: n, cp
    REAL(KIND=dp) :: cm
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    REAL(KIND=dp) :: g, rho, eta, d, dhdx, dhdy, dU2dz2, nn, detadz
    REAL(KIND=dp) :: grad1, grad2                            
    LOGICAL :: Stat
    INTEGER :: t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    nn = 1.0/cm

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO
       ! Point value
       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )
       eta = SUM( Viscosity(1:n) * Basis(1:n) )
       detadz = SUM( Viscosity(1:n) * dBasisdx(1:n,dim) )
       grad1 = SUM( SurfGrad1(1:n) * Basis(1:n) )
       grad2 = SUM( SurfGrad2(1:n) * Basis(1:n) )
       d = SUM( LocalDepth(1:n) * Basis(1:n) )


       IF (cp==1) dU2dz2 = nn * (rho * g / eta)**nn * grad1 * (1.0 + d * detadz / eta)
       IF (cp==2) dU2dz2 = nn * (rho * g / eta)**nn * grad2 * (1.0 + d * detadz / eta)

       ! Non linear case
       IF (nn > 1.0) THEN
               dU2dz2 = dU2dz2 * (d * SQRT(grad1**2.0 + grad2**2.0))**(nn-1.0)
       END IF
       
       FORCE(1:n) = FORCE(1:n) - dU2dz2 * IP % s(t) * detJ * Basis(1:n) 

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUV
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixW(  STIFF, FORCE, Element, n, VeloU, VeloV)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), VeloU(:), VeloV(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ, &
                     dU2dxz, dV2dyz
    LOGICAL :: Stat
    INTEGER :: t, p,q , DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    DIM = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .TRUE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO

       dU2dxz = SUM(VeloU(1:n)*ddBasisddx(1:n,1,DIM))
       dV2dyz = 0.0d0
       IF (DIM==3) dV2dyz = SUM(VeloV(1:n)*ddBasisddx(1:n,2,3))
       

       FORCE(1:n) = FORCE(1:n) + (dU2dxz + dV2dyz) * IP % s(t) * detJ * Basis(1:n) 

    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixW

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixP(  STIFF, FORCE, Element, n)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixP
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCW(  STIFF, FORCE, Element, n, VeloU, VeloV )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), veloU(:), veloV(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ, Normal(3), grad, dUdx, dVdy  
    LOGICAL :: Stat
    INTEGER :: t, DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    DIM = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       dUdx = SUM( VeloU(1:n) * dBasisdx(1:n,1) )
       dVdy = 0.0e0
       IF (DIM==3) dVdy = SUM( VeloV(1:n) * dBasisdx(1:n,2) )

       grad = - (dUdx + dVdy) 

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCW
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCP(  STIFF, FORCE, Element, n, Density, & 
                      Gravity)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), density(:), Gravity(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), rho, g, grad
    LOGICAL :: Stat
    INTEGER :: t, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )

       grad = - rho * g 

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCP
!------------------------------------------------------------------------------
END SUBROUTINE SIASolver
!------------------------------------------------------------------------------
!> Allows to have the SIA Velocity as primary variables (and not exported one)
!> Allow then to have access to the Previous values to construct a stable 
!> time discretization sheme. 
SUBROUTINE SIAVariable( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
END SUBROUTINE SIAVariable
!------------------------------------------------------------------------------


