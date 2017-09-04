!/*****************************************************************************
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
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.11.2005
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
!>  Computes the curvature from the level set function. Additional diffusion may 
!>  be added in order to limit the singular curvature peaks.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetCurvature( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     USE SolverUtils
     USE Integration

     IMPLICIT NONE
!------------------------------------------------------------------------------ 
     TYPE(Model_t), TARGET :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: Timestep
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     TYPE(Nodes_t) :: ElementNodes, ParentNodes
     TYPE(Element_t),POINTER :: Element, Parent
     TYPE(Variable_t), POINTER :: SurfSol

     INTEGER :: i,j,k,n,pn,t,istat,bf_id,CoordinateSystem
     REAL(KIND=dp) :: Norm, Coeff, Diff, Alpha, Val, Delta
     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER, POINTER :: CurvPerm(:), SurfPerm(:)
     REAL(KIND=dp), POINTER :: Curvature(:),ForceVector(:), Curv(:),Surface(:) 
     REAL(KIND=dp), ALLOCATABLE :: LocalStiffMatrix(:,:),LocalForce(:),Surf(:)
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,totat,st,totst
#else
     REAL(KIND=dp) :: at,totat,st,totst,CPUTime
#endif
     CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName
     LOGICAL :: GotIt, Stat, AllocationsDone = .FALSE.

     SAVE LocalStiffMatrix,LocalForce, ElementNodes,ParentNodes,Surf,AllocationsDone

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     CoordinateSystem = CurrentCoordinateSystem()

     Curvature => Solver % Variable % Values
     CurvPerm => Solver % Variable % Perm
     IF ( SIZE( Curvature ) == 0 ) RETURN

     LevelSetVariableName = ListGetString(Solver % Values,'LevelSet Variable',GotIt) 
     IF(GotIt) THEN
       SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(LevelSetVariableName) )
     ELSE  
       SurfSol => VariableGet( Solver % Mesh % Variables, 'Surface' )
     END IF
     Surface => Surfsol % Values
     SurfPerm => SurfSol % Perm

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     Norm = Solver % Variable % Norm

     Diff = ListGetConstReal(Solver % Values,'Curvature Diffusion',GotIt)

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Solver % Mesh % MaxElementNodes

       ALLOCATE( ElementNodes % x( N ),   &
           ElementNodes % y( N ),   &
           ElementNodes % z( N ),   &
           ParentNodes % x( N ),   &
           ParentNodes % y( N ),   &
           ParentNodes % z( N ),   &
           LocalForce( N ),       &
           LocalStiffMatrix( N, N ), &
           Surf( N ), &
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'CurvatureSolve', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     at = CPUTime()
     
     CALL Info( 'LevelSetCurvature','-------------------------------------', Level=4 )
     CALL Info( 'LevelSetCurvature','Solving Level set curvature', Level=4 )
     CALL Info( 'LevelSetCurvature','-------------------------------------', Level=4 )

     CALL DefaultInitialize()

!------------------------------------------------------------------------------
!      Do the assembly for bulk elements
!------------------------------------------------------------------------------

     DO t=1,Solver % NumberOfActiveElements
       
       Element => Solver % Mesh % Elements(Solver % ActiveElements(t))
       n = Element % TYPE % NumberOfNodes
       NodeIndexes => Element % NodeIndexes
       Model % CurrentElement => Element
 
       ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

       Surf = Surface( SurfPerm(NodeIndexes) )

       CALL LocalMatrix( LocalStiffMatrix, LocalForce, &
           Surf, Element, n, ElementNodes )
         
       CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

     END DO     
     CALL DefaultFinishBulkAssembly()

!------------------------------------------------------------------------------
!      Do the assembly for boundary elements
!------------------------------------------------------------------------------

     DO t=Solver % Mesh % NumberOfBulkElements + 1, &
         Solver % Mesh % NumberOfBulkElements + &
         Solver % Mesh % NumberOfBoundaryElements
      
       Element => Solver % Mesh % Elements(t)
!------------------------------------------------------------------------------
       DO i=1,Model % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN

           n = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           Model % CurrentElement => Element
          
           IF ( ANY( CurvPerm(NodeIndexes) <= 0 ) ) CYCLE         
           
           IF ( .NOT. ListGetLogical(Model % BCs(i) % Values, &
               'Levelset Curvature BC',gotIt) ) CYCLE

           ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
           
           Parent => Element % BoundaryInfo % Left
          
           stat = ASSOCIATED( Parent )
           IF ( stat ) stat = stat .AND. ALL(CurvPerm(Parent % NodeIndexes) > 0)
          
           IF ( .NOT. stat ) THEN
             Parent => Element % BoundaryInfo % Right
             
            stat = ASSOCIATED( Parent )
            IF ( stat ) stat = ALL(CurvPerm(Parent % NodeIndexes(1:pn)) > 0)
            
            IF ( .NOT. stat )  THEN
              CALL Warn( 'LevelSetCurvature', &
                  'No curvature solution available for specified boundary' )
              CYCLE
            END IF
          END IF
          
          pn = Parent % TYPE % NumberOfNodes           
          ParentNodes % x(1:pn) = Solver % Mesh % Nodes % x(Parent % NodeIndexes)
          ParentNodes % y(1:pn) = Solver % Mesh % Nodes % y(Parent % NodeIndexes)
          ParentNodes % z(1:pn) = Solver % Mesh % Nodes % z(Parent % NodeIndexes)

          Surf(1:pn) = Surface(SurfPerm(Parent % NodeIndexes))
          
!------------------------------------------------------------------------------
!             Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
          
          CALL LocalBoundary( LocalStiffMatrix, LocalForce,  &
              Surf, Element, Parent, n, pn, ElementNodes, ParentNodes )
!------------------------------------------------------------------------------
!             Update global matrices from local matrices
!------------------------------------------------------------------------------
          CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
!------------------------------------------------------------------------------
           END IF   
        END DO    
      END DO   
!------------------------------------------------------------------------------

     CALL DefaultFinishAssembly()
     
     at = CPUTime() - at
     WRITE(Message,'(a,F8.2)') 'Assembly done in time (s):',at
     CALL Info( 'LevelSetCurvature',Message, Level=4 )

!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
     st = CPUTime()
     CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, &
         Curvature, Norm, 1, Solver )
     
     st = CPUTIme()-st
     WRITE(Message,'(a,F8.2)') 'Solution done in time (s):',st
     CALL Info( 'LevelSetCurvature',Message, Level=4 )

!------------------------------------------------------------------------------

     Coeff = ListGetConstReal(Solver % Values,'Curvature Coefficient',GotIt) 
     IF(GotIt) THEN
       Curvature = Coeff * Curvature
     END IF

     Alpha = ListGetConstReal(Model % Simulation,'Levelset Bandwidth',GotIt) 
     IF(.NOT. GotIt) Alpha = ListGetConstReal(Solver % Values,'Levelset Bandwidth',GotIt)      
     IF(GotIt) THEN
       DO i=1,SIZE(CurvPerm) 
         j = CurvPerm(i)
         k = SurfPerm(i)
         IF(j == 0 .OR. k == 0) CYCLE
         Val = Surface(k)

         IF( Val < -Alpha) THEN
           Delta = 0.0d0
         ELSE IF(Val > Alpha) THEN
           Delta = 0.0d0
         ELSE
           Delta = (1.0d0 + COS( (Val/Alpha) * PI ) ) / (2.0d0 * Alpha)
         END IF
         
         Curvature(j) = Delta * Curvature(j)
       END DO
     END IF


CONTAINS

   SUBROUTINE LocalMatrix( StiffMatrix,ForceVector, Surf, &
       Element, n, Nodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Surf
     REAL(KIND=dp), DIMENSION(:,:) :: StiffMatrix
     INTEGER :: n
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ
     REAL(KIND=dp) :: Grad(3),GradAbs,A,B,xpos
     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ
     REAL(KIND=dp) :: s,u,v,w
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     Grad = 0.0d0
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
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
       
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx)

       s = detJ * S_Integ(t)

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( Coordinates == AxisSymmetric .OR. Coordinates == CylindricSymmetric ) THEN
         xpos = SUM( Nodes % x(1:n) * Basis(1:n) )
         s = xpos * s
       END IF

!------------------------------------------------------------------------------

       DO i=1,dim
         Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
       END DO
       GradAbs = SQRT( SUM( Grad * Grad ) )
       IF ( GradAbs > 10*AEPS ) THEN
         Grad = Grad / GradAbs
       END IF

!------------------------------------------------------------------------------
!      Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,n

         DO q=1,n
           A = Basis(q) * Basis(p)
           A = A + Diff * SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )

           StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
         END DO

         B = - SUM( dBasisdx(p, 1:dim) * Grad(1:dim))  
         ForceVector(p) = ForceVector(p) + s * B        

       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix


!------------------------------------------------------------------------------
   SUBROUTINE LocalBoundary( BoundaryMatrix, BoundaryVector, &
        Surf, Element, Parent, n, pn, ElementNodes, ParentNodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), Surf(:)
     TYPE(Nodes_t)   :: ElementNodes, ParentNodes
     TYPE(Element_t), POINTER :: Element, Parent
     INTEGER :: n, pn

     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), Normal(3), Grad(3), GradAbs, detJ, A
     REAL(KIND=dp) :: ParentBasis(pn), ParentdBasisdx(pn,3)
     REAL(KIND=dp) :: u,v,w,s,x(n),y(n),z(n),xpos, Force
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,j,k,t,p,q,N_Integ,dim

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0d0
     BoundaryMatrix = 0.0d0
     dim = CoordinateSystemDimension()
     Grad = 0.0d0

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
!     Basis function values & derivates at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,ElementNodes,u,v,w,detJ,Basis,dBasisdx )
       s = S_Integ(t) * detJ

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        xpos = SUM( ElementNodes % x(1:n)*Basis(1:n) )
        s = xpos * s
      END IF

      Normal = Normalvector( Element, ElementNodes, u, v, .TRUE. )
      
      !------------------------------------------------------------------------------
      ! Need parent element basis etc., for computing normal derivatives on boundary.
      !------------------------------------------------------------------------------
      
      DO i = 1,n
        DO j = 1,pn
          IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            x(i) = Parent % TYPE % NodeU(j)
            y(i) = Parent % TYPE % NodeV(j)
            z(i) = Parent % TYPE % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      
      u = SUM( Basis(1:n) * x(1:n) )
      v = SUM( Basis(1:n) * y(1:n) )
      w = SUM( Basis(1:n) * z(1:n) )
      
      stat = ElementInfo( Parent, ParentNodes, u, v, w, detJ, ParentBasis, ParentdBasisdx )
          
      DO j = 1,DIM
        Grad(j) = SUM( ParentdBasisdx(1:pn,j) * Surf(1:pn) )
      END DO
      GradAbs = SQRT( SUM(Grad(1:dim) * Grad(1:dim)) )

      IF ( GradAbs > 10*AEPS ) THEN         
        Grad = Grad / GradAbs
      END IF

!------------------------------------------------------------------------------
      
      DO p=1,n
        DO q=1,n
          A = Diff * SUM( Normal(1:dim) * dBasisdx(q,1:dim)) * Basis(p)
          BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + s * A
        END DO
      END DO

      Force = SUM( Normal(1:dim) * Grad(1:dim) )
      DO q=1,N
        BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
      END DO

    END DO

  END SUBROUTINE LocalBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE LevelSetCurvature
!------------------------------------------------------------------------------
