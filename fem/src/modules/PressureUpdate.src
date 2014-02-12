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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Mika  Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 Dec 2003
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  A pressure computation for the Krylov accelerated consistent splitting scheme.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE PressureSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Types
  USE CRSMatrix
  USE ParallelUtils  
  USE Lists
  USE SparIterSolve

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found, Parallel, OutflowBC, BlockPreconditioning
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC

  REAL(KIND=dp) :: Norm, Visc
  INTEGER :: n, nb, nd, t, istat, active, j, i, dim, k
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), r(:), TmpVec(:), &
      CVelo(:,:), PVelo(:,:), PPres(:,:)

  TYPE(Matrix_t), POINTER :: M
  REAL(KIND=dp), POINTER :: Mx(:),Mb(:),Mr(:)
  INTEGER, ALLOCATABLE :: Indexes(:)

  TYPE(Variable_t), POINTER :: VeloVar, DivVar, PresVar



  SAVE STIFF, LOAD, FORCE, r, AllocationsDone, Indexes, TmpVec, CVelo, PVelo, PPres
  !------------------------------------------------------------------------------
 
  ! Traction computations on boundary
  TYPE(Element_t),POINTER :: Parent
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), ALLOCATABLE :: TotalForce(:,:), TotalArea(:), mu(:) 
  REAL(KIND=dp) :: Traction(3), Area
  INTEGER, ALLOCATABLE :: ParentIndexes(:)
  TYPE(Nodes_t) :: ElementNodes, ParentNodes  
  INTEGER, POINTER :: FlowPerm(:)
  LOGICAL :: CalculateTraction
  INTEGER :: np, ndp, nlen
  CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryName

  SAVE TotalForce, TotalArea, ParentIndexes, mu
  !-------------------------------------------------------------------------------


!------------------------------------------------------------------------------

  Parallel = ParEnv % PEs > 1
  !PRINT *, 'This run is parallel...', Parallel
  BlockPreconditioning = GetLogical( GetSolverParams(), 'Block Preconditioning', Found )
  IF ( .NOT. Found )  BlockPreconditioning = .TRUE.

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  dim = CoordinateSystemDimension()
  Mesh => GetMesh()

  VeloVar => VariableGet( Mesh % Variables, "VelocityTot" )
  DivVar => VariableGet( Mesh % Variables, "Divergence" ) 
  if (BlockPreconditioning) then
     PresVar => VariableGet( Mesh % Variables, "Flow" )
     IF ( .NOT. ASSOCIATED(PresVar) ) &
          CALL Fatal( 'PressureSolver', 'The coupled flow variable Flow was not found' )     
  end if

  IF ( .NOT. ASSOCIATED(VeloVar) ) &
      CALL Fatal( 'PressureUpdate', 'Velocity Variable VelocityTot was not found' )
  IF ( .NOT. ASSOCIATED(DivVar) ) &
      CALL Fatal( 'PressureUpdate', 'Divergence Variable was not found' )


  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     t = Solver % Matrix % NumberOfRows
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), Indexes(N), r(t), &
         CVelo(dim,n), PVelo(dim,n), PPres(n,2), &
         TotalForce(Model % NumberOfBCs,3), &
         TotalArea(Model % NumberOfBCs), & 
         ParentIndexes(n), &    
         mu(n), &
         STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PressureUpdate', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

   !System assembly:
   !----------------
   Active = GetNOFActive()
   CALL DefaultInitialize()
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementDOFs( Indexes )
      nb = GetElementNOFBDOFs()


      ! Get previous elementwise velocity iterate:
      !-------------------------------------------
      if (BlockPreconditioning) then
         DO i=1,dim
            PVelo(i,1:nd) = 0.0d0 ! This should work here since previous velocity is divergence-free  
            CVelo(i,1:nd) = VeloVar % Values( &
                 VeloVar % DOFs * (VeloVar % Perm(Indexes(1:nd))-1)+i)            
         end DO
         PPres(1:nd,1) = DivVar % Values( DivVar % Perm(Indexes(1:nd)) )
         PPres(1:nd,2) = PresVar % Values( &
              PresVar % DOFs*(PresVar % &
              Perm(Indexes(1:nd))-1)+dim+1)     
      else

         DO i=1,dim
            CVelo(i,1:nd) = VeloVar % Values( &
                 VeloVar % DOFs * (VeloVar % Perm(Indexes(1:nd))-1)+i)

            PVelo(i,1:nd) = VeloVar % PrevValues( &
                 VeloVar % DOFs * (VeloVar % Perm(Indexes(1:nd))-1)+i,1)
         END DO

         PPres(1:nd,1) = DivVar % Values( DivVar % Perm(Indexes(1:nd)) )
         PPres(1:nd,2) = Solver % Variable % PrevValues( &
              Solver % Variable % Perm( Indexes(1:nd) ),1)

      end if

      LOAD = 0.0d0
      !BodyForce => GetBodyForce()
      !IF ( ASSOCIATED(BodyForce) ) &
      !   Load(1:n) = GetReal( BodyForce, 'Source', Found )

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd, nd+nb )
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO

   CALL DefaultFinishAssembly()
   

   Active = GetNOFBoundaryElements()
   !PRINT *, Active
  
   DO t=1, Active           !Solver % Mesh % NumberOfBoundaryElements   
     Element => GetBoundaryElement(t)
     IF ( .NOT. ActiveBoundaryElement() ) CYCLE

     n  = GetElementNOFNodes()
     nd = GetElementDOFs( Indexes )

     IF ( GetElementFamily() == 1 ) CYCLE

     BC => GetBC()
     IF ( ASSOCIATED( BC ) ) THEN
       OutflowBC= ListGetLogical( BC, 'Outflow Boundary', Found ) 
       IF (OutFlowBC) THEN

         DO i=1,dim
           CVelo(i,1:nd) = VeloVar % Values( &
               VeloVar % DOFs*(VeloVar % Perm(Indexes(1:nd))-1)+i)

           PVelo(i,1:nd) = VeloVar % PrevValues( &
               VeloVar % DOFs*(VeloVar % Perm(Indexes(1:nd))-1)+i,1)
         END DO


         DO i=1,nd

           j = Solver % Variable % Perm( Indexes(i) )
           k = DivVar % Perm( Indexes(i) )


           CALL ZeroRow( Solver % Matrix, j)
           CALL SetMatrixElement( Solver % Matrix, j, j, 1.0d0 )

           if (BlockPreconditioning) then
              IF ( i > n ) THEN
                 Solver % Matrix % RHS(j) = 0.0d0
              ELSE
                 Solver % Matrix % RHS(j) = PresVar % Values( PresVar % DOFs * (PresVar % &
                      Perm(Indexes(i))-1)+dim+1) - DivVar % Values(k)
              END IF              
           else
              IF ( i > n ) THEN
                 Solver % Matrix % RHS(j) = 0.0d0
              ELSE
                 Solver % Matrix % RHS(j) = Solver % Variable % PrevValues(j,1) - DivVar % Values(k)
              END IF
              !Solver % Matrix % RHS(j) = 1.0d0
           end if
         END DO

       END IF
     END IF
   END DO


   !CALL DefaultDirichletBCs()


   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()


   IF ( .NOT. BlockPreconditioning) THEN
      !-------------------------------------------------------------------------
      ! Compute tractions if desired
      !-------------------------------------------------------------------------
      FlowPerm => Solver % Variable % Perm
      TotalForce = 0.0d0
      TotalArea = 0.0d0
      Active = GetNOFBoundaryElements()

      DO t=1, Active

         Element => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement() ) CYCLE

         !------------------------------------------------------------------------------
         ! Skip degenerate element types
         !------------------------------------------------------------------------------
         IF ( Element % TYPE % ElementCode == 101 ) CYCLE
         IF ( (dim > 2) .AND. Element % TYPE % ElementCode == 202 ) CYCLE    
         !------------------------------------------------------------------------------

         n  = GetElementNOFNodes()
         nd = GetElementDOFs( Indexes )

         DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == &
                 Model % BCs(i) % Tag ) THEN

               IF ( .NOT. ListGetLogical( Model % BCs(i) % Values, &
                    'Calculate Fluidic Force', Found ) ) CYCLE

               Parent => Element % BoundaryInfo % Left
               Found = ASSOCIATED( Parent )
               IF (Found) Found = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

               IF ( .NOT. Found) THEN
                  Parent => ELement % BoundaryInfo % Right
                  Found = ASSOCIATED( Parent )
                  IF (Found) Found = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

                  IF ( .NOT. Found )  CALL Fatal( 'PressureUpdate', &
                       'No parent element can be found for given boundary element' )
               END IF

               np = Parent % TYPE % NumberOfNodes
               ndp = GetElementDOFs( ParentIndexes, Parent )

               k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
                    'Material' )
               Material => Model % Materials(k) % Values

               mu(1:np) = ListGetReal( Material, 'Viscosity', &
                    np, Parent % NodeIndexes(1:np) )   


               DO j=1,dim
                  CVelo(j,1:ndp) = VeloVar % Values( &
                       VeloVar % DOFs*(VeloVar % Perm(ParentIndexes(1:ndp))-1)+j)
               END DO
               PPres(1:ndp,1) = Solver % Variable % Values( Solver % Variable % Perm( ParentIndexes(1:ndp) ) )


               Traction = 0.0d0
               CALL SurfaceForceIntegration(Element, Parent, Traction, Area, &
                    CVelo, PPres, mu, np, ndp)

               TotalForce(i,1:3) = TotalForce(i,1:3) + Traction(1:3)
               TotalArea(i) = TotalArea(i) + Area

            END IF
         END DO
      END DO


      DO k=1, Model % NumberOfBCs
         IF ( .NOT. ListGetLogical(Model % BCs(k) % Values,'Calculate Fluidic Force', Found ) ) CYCLE
         IF( Model % NumberOfBCs < 10 ) THEN
            WRITE( BoundaryName, '("bc ",I1)') k
         ELSE
            WRITE( BoundaryName, '("bc ",I2)') k
         END IF

         nlen = LEN_TRIM(BoundaryName)

         CALL Info('ForceCompute','Forces on Boundary '//BoundaryName(1:nlen),Level=4 )
         WRITE( Message, '("Fluidic Force (X,Y,Z):", 3ES17.6E2)') TotalForce(k,1:3)
         CALL Info( 'ForceCompute', Message, Level=4 ) 
         WRITE( Message, '("Contact Area:   ", ES17.6E2)') TotalArea(k)
         CALL Info( 'ForceCompute', Message, Level=4 )


         CALL ListAddConstReal( Model % Simulation, &
              'res: contact force 1 '//BoundaryName(1:nlen), TotalForce(k,1) )
         CALL ListAddConstReal( Model % Simulation, &
              'res: contact force 2 '//BoundaryName(1:nlen), TotalForce(k,2) )
         IF ( DIM > 2 )  CALL ListAddConstReal( Model % Simulation, &
              'res: contact force 3 '//BoundaryName(1:nlen), TotalForce(k,3) )

         CALL ListAddConstReal( Model % Simulation, & 
              'res: contact force area '//BoundaryName(1:nlen), TotalArea(k) )

      END DO

   END IF


CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd, ntot )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP, s, Ddivu
    LOGICAL :: Stat
    INTEGER :: i, t, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      s = IP % s(t) * DetJ

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      Ddivu = 0.0d0
      DO i=1,dim
        Ddivu = Ddivu + ( SUM( Cvelo(i,1:nd) * dBasisdx(1:nd,i) ) - &
            SUM( Pvelo(i,1:nd) * dBasisdx(1:nd,i) ) ) / dt 
      END DO

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      DO p=1,n
        DO q=1,n
          stiff(p,q) = stiff(p,q) + s * SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )
        END DO
        FORCE(p) = FORCE(p) - s * Ddivu * Basis(p)
      END DO

      !STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + s * &
      !       MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      !FORCE(1:nd) = FORCE(1:nd) +  s * LoadAtIP * Basis(1:nd)

    END DO
    
    Force(1:n) = Force(1:n) - MATMUL( Stiff(1:n, 1:n), PPres(1:n,1) ) + &
        MATMUL( Stiff(1:n, 1:n), Ppres(1:n,2) )
    

    ! Eliminate the bubble and edge degrees of freedom if any
    
    DO i = n+1,ntot
      FORCE(i)   = 0.0d0
      STIFF(i,:) = 0.0d0
      STIFF(:,i) = 0.0d0
      STIFF(i,i) = 1.0d0
    END DO
    

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!----------------------------------------------------------------------------------
  SUBROUTINE SurfaceForceIntegration(Element, Parent, Traction, Area, &
     Velo, Pres, Viscosity, np, ndp)
!----------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element, Parent  
  REAL(kind=dp) :: Traction(3), Area, Velo(:,:), &
      Viscosity(:), Pres(:,:)
  INTEGER :: np, ndp
!----------------------------------------------------------------------------------
  LOGICAL :: stat
  INTEGER :: N_Integ, t, i, dim
  REAL(kind=dp) :: u, v, w, detJ, s, Basis(ndp), dBasisdx(ndp,3), Normal(3), &
      Visc, Lambda, ReGrad(3,3), ImGrad(3,3), ReDiv, ImDiv, ReD(3,3), ImD(3,3), &
      ReP, ImP, tmpmat(3,ndp), tmp = 0.0d0, r1, r2, r3, r
  REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

  TYPE(Nodes_t) :: Nodes, ParentNodes
  SAVE Nodes, ParentNodes

  SAVE tmp
!----------------------------------------------------------------------------------
  dBasisdx = 0.0d0
  Traction = 0.0d0
  Area = 0.0d0
  dim = CoordinateSystemDimension()
  CALL GetElementNodes( Nodes, Element )
!----------------------------------------------------------------------------------
! The normal vector to the boundary element
!----------------------------------------------------------------------------------
  SELECT CASE( Element % TYPE % NumberOfNodes )
  CASE(2)
    Normal = Normalvector(Element, Nodes, 0.0d0, 0.0d0, .TRUE.)    
  CASE( 3 )
    Normal = Normalvector(Element, Nodes, 0.3d0, 0.3d0, .TRUE.)          
  CASE( 4 )
    Normal = Normalvector(Element, Nodes, 0.0d0, 0.0d0, .TRUE.) 
  END SELECT

  CALL GetElementNodes( ParentNodes, Parent )
!----------------------------------------------------------------------------------
! Compute field derivatives at the integration points of the parent element
!----------------------------------------------------------------------------------
  IntegStuff = Gausspoints(Parent)
  U_Integ => IntegStuff % u
  V_Integ => IntegStuff % v
  W_Integ => IntegStuff % w
  S_Integ => IntegStuff % s
  N_Integ =  IntegStuff % n 

!------------------------------------------------------------------------------
  DO t=1,N_Integ
!------------------------------------------------------------------------------
    u = U_Integ(t)
    v = V_Integ(t)
    w = W_Integ(t)
!------------------------------------------------------------------------------
!   Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
    stat = ElementInfo( Parent, ParentNodes, u, v, w, &
        detJ, Basis, dBasisdx )

    s = detJ * S_Integ(t)

    Visc = SUM( Viscosity(1:np) * Basis(1:np) ) 
    ReP = SUM(Pres(1:ndp,1) * Basis(1:ndp))

    tmpmat(1,1:ndp) = Velo(1,1:ndp)
    tmpmat(2,1:ndp) = Velo(2,1:ndp)
    IF ( dim > 2 ) THEN
       tmpmat(3,1:ndp) = Velo(3,1:ndp)
    ELSE
       tmpmat(3,1:ndp) = 0.0d0
    END IF

    ReGrad = MATMUL( tmpmat, dBasisdx )

    ReD = Visc * ( ReGrad + TRANSPOSE(ReGrad) )

    Traction(1:3) = Traction(1:3) - ReP * Normal(1:3) + MATMUL( ReD, Normal)

  END DO
  
  Area = ElementArea(Solver % Mesh, Element, Element % TYPE % NumberOfNodes)

  Traction = -Area/N_Integ * Traction


!----------------------------------------------------------------------------------
END SUBROUTINE SurfaceForceIntegration
!------------------------------------------------------------------------------------



!------------------------------------------------------------------------------
END SUBROUTINE PressureSolver
!------------------------------------------------------------------------------
