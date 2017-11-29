!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
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


!------------------------------------------------------------------------------
!> Initialization for the primary solver: SmitcSolver
!------------------------------------------------------------------------------
  SUBROUTINE SmitcSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Solver_t) :: Solver
    TYPE(Model_t)  :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 3 )
      CALL ListAddString( SolverParams, 'Variable', 'Deflection' )
    END IF
    


    
    CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )
!------------------------------------------------------------------------------
  END SUBROUTINE SmitcSolver_Init
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
!>  Solve the Reissner-Mindlin equations i.e. displacement equations for
!> elastic plates. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
 SUBROUTINE SmitcSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Solver_t):: Solver
     TYPE(Model_t) :: Model
 
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: k,n,nd,t,bf_id,mat_id,istat,DOFs
 
     TYPE(Element_t),POINTER :: Element
     TYPE(Nodes_t) :: ElementNodes
     REAL(KIND=dp) :: Norm
     TYPE(ValueList_t), POINTER :: Material, BodyForce
 
     LOGICAL :: AllocationsDone = .FALSE., HoleCorrection, &
         got_mat_id, got_bf_id, NeglectSprings, EigenOrHarmonic, Found

     INTEGER, POINTER, SAVE ::  Indexes(:)
     INTEGER :: MaxIter, iter
     TYPE(ValueList_t), POINTER :: SolverParams
     
     REAL(KIND=dp), ALLOCATABLE :: &
                     STIFF(:,:), Load(:), Load2(:), FORCE(:), &
                     Poisson(:), Thickness(:), Young(:), Tension(:), &
                     MASS(:,:), DAMP(:,:), Density(:), &
                     DampingCoef(:), HoleFraction(:), HoleSize(:), SpringCoef(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: HoleType
     LOGICAL :: GotIt, GotHoleType
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,st
#else
     REAL(KIND=dp) :: at,st,CPUTime
#endif
     SAVE STIFF, MASS, Load, Load2, FORCE, ElementNodes, &
          Poisson, Density, Young, Thickness, Tension, AllocationsDone, &
          DAMP, DampingCoef, HoleFraction, HoleSize, SpringCoef

!
!    Allocate some permanent storage, this is done first time only
!    -------------------------------------------------------------
     DOFs = Solver % Variable % DOFs

     IF ( .NOT. AllocationsDone ) THEN
       N = Solver % Mesh % MaxElementDOFs
       ALLOCATE( Indexes( N ),      &
                 FORCE( DOFs*N ),      &
                 STIFF( DOFs*N, DOFs*N ), &
                 MASS( DOFs*N, DOFs*N ),  &
                 DAMP( DOFs*N, DOFs*N ),  &
                 Load( N ), Load2(N), Poisson( N ), Young( N ),   &
                 Density ( N ), Thickness( N ), DampingCoef( N ), &
                 Tension( N ), HoleFraction( N ), HoleSize( N ),  &
                 SpringCoef( N ), STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal('SmitcSolver','Memory allocation error')
       END IF
       AllocationsDone = .TRUE.
     END IF

!
!    Do some additional initialization, and go for it
!    ------------------------------------------------
     CALL Info( 'SmitcSolver', '--------------------------------------------------',Level=4 )
     CALL Info( 'SmitcSolver', 'Solving the Reissner-Mindlin equations for plates',Level=4 )     
     CALL Info( 'SmitcSolver', '--------------------------------------------------',Level=4 )

     SolverParams => GetSolverParams()
     
     
     EigenOrHarmonic = EigenOrHarmonicAnalysis() &
         .OR. ListGetLogical( SolverParams,'Harmonic Mode',Found ) 

     CALL DefaultStart()     

     MaxIter = GetInteger( SolverParams, &
         'Nonlinear System Max Iterations',GotIt )
     IF ( .NOT. GotIt ) MaxIter = 1

     DO iter=1,MaxIter
    
       at = CPUTime()
       CALL DefaultInitialize()

       !
       ! These keywords enable that the use of a second parameter set for the
       ! same elements where the material properties are given in an additional
       ! body. May be used to model microphone and its backplate, for example:
       ! --------------------------------------------------------------------
       mat_id = ListGetInteger( SolverParams, 'Material Index',got_mat_id, &
           minv=1, maxv=Model % NumberOFMaterials )

       Material => NULL()
       IF(got_mat_id) THEN
         Material => Model % Materials(mat_id) % Values
       END IF

       bf_id = ListGetInteger( SolverParams, 'Body Force Index', &
           got_bf_id,  minv=1, maxv=Model % NumberOFBodyForces )

       BodyForce => NULL()
       IF(got_bf_id) THEN
         BodyForce => Model % Materials(mat_id) % Values
       END IF
       HoleCorrection = ListGetLogical( SolverParams, &
           'Hole Correction',gotIt )

       !
       !    Do the assembly:
       !    ----------------
       DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t)

         n = GetElementDOFs(Indexes)
         n = GetElementNOFNodes()
         CALL GetElementNodes(ElementNodes)

         IF(.NOT. got_bf_id) BodyForce => GetBodyForce()

         Load = 0.0d0
         Load2 = 0.0d0

         ! There may be three forces, which should be introduced
         ! in the following order
         ! ------------------------------------------------------
         IF(ASSOCIATED(BodyForce)) THEN
           Load(1:n) = GetReal( BodyForce, 'Pressure',  GotIt)
           IF(GotIt) THEN
             Load2(1:n) = GetReal( BodyForce, 'Pressure B', GotIt )
             IF(GotIt) THEN
               Load(1:n) = Load(1:n) + Load2(1:n)
               Load2(1:n) = GetReal( BodyForce, 'Pressure C', GotIt)
               IF(GotIt) Load(1:n) = Load(1:n) + Load2(1:n)
             END IF
           END IF
         END IF

         IF (.NOT. got_mat_id) Material => GetMaterial()

         Poisson(1:n) = GetReal( Material, 'Poisson ratio' )
         Density(1:n) = GetReal( Material, 'Density' )
         Thickness(1:n) = GetReal( Material,'Thickness' )
         Young(1:n) = GetReal( Material,'Youngs modulus' )
         Tension(1:n) = GetReal( Material, 'Tension', GotIt )

         ! In some cases it is preferable that the damping and spring
         ! coefficients are related to the body force.
         ! ----------------------------------------------------------
         IF ( ASSOCIATED(BodyForce) ) THEN
           DampingCoef(1:n) = GetReal( BodyForce, 'Damping', GotIt )
         ELSE
           GotIt = .FALSE.
         END IF
         IF (.NOT. GotIt) DampingCoef(1:n) = &
             GetReal( Material, 'Damping', GotIt )
         IF (.NOT. GotIt ) DampingCoef(1:n) = 0.0d0

         IF ( ASSOCIATED(BodyForce)) THEN
           SpringCoef(1:n) = GetReal( BodyForce, 'Spring', GotIt )
         ELSE
           GotIt = .FALSE.
         END IF
         IF (.NOT. GotIt) SpringCoef(1:n) = &
             GetReal( Material, 'Spring', GotIt )
         IF (.NOT. GotIt ) SpringCoef(1:n) = 0.0d0

         IF(HoleCorrection) THEN
           HoleType = GetString(Material,'Hole Type',GotHoleType)
           IF(GotHoleType) THEN
             HoleSize(1:n) = GetReal( Material, 'Hole Size' )
             HoleFraction(1:n) = GetReal( Material, 'Hole Fraction' )
           END IF
         END IF

         ! Get element local matrix, and rhs vector
         !-----------------------------------------
         CALL LocalMatrix(  STIFF, DAMP, MASS, FORCE, Load, &
             Element,n, DOFs, ElementNodes, DampingCoef, SpringCoef )

         IF( TransientSimulation ) &
             CALL Default2ndOrderTime( MASS,DAMP,STIFF,FORCE )

         ! Update global matrix and rhs vector from local matrix & vector
         !---------------------------------------------------------------
         CALL DefaultUpdateEquations( STIFF, FORCE )

         IF ( EigenOrHarmonic ) THEN
           CALL DefaultUpdateMass( MASS )
           CALL DefaultUpdateDamp( DAMP )
         END IF
       END DO
       CALL DefaultFinishBulkAssembly()

       ! No Flux BCs
       CALL DefaultFinishBoundaryAssembly()
       CALL DefaultFinishAssembly()
       
       !------------------------------------------------------------------------------

       ! Dirichlet boundary conditions
       !------------------------------
       CALL DefaultDirichletBCs()

       at = CPUTime() - at

       WRITE (Message,*) 'Assembly (s): ',at
       CALL Info('SmitcSolver',Message,Level=4)

       !------------------------------------------------------------------------------

       ! Solve the system and we are done
       !---------------------------------
       st = CPUTime()
       Norm =  DefaultSolve()

       st = CPUTime() - st
       WRITE (Message,*) 'Solve (s): ',st
       CALL Info('SmitcSolver',Message,Level=4)

       IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
     END DO
     
     CALL DefaultFinish()
       
!------------------------------------------------------------------------------
 
   CONTAINS

!------------------------------------------------------------------------------
     SUBROUTINE LocalMatrix( STIFF, DAMP, MASS, &
          Force, Load, Element, n, DOFs, Nodes, DampingCoef, SpringCoef )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), DAMP(:,:), &
            MASS(:,:), Force(:), Load(:), DampingCoef(:), SpringCoef(:)
       INTEGER :: n, DOFs
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), &
                           Curvature(3,100), ShearStrain(2,100), &
                           Ematrix(3,3), Gmatrix(2,2), Tmatrix(2,2)
       REAL(KIND=dp) :: detJ,U,V,W,S,Kappa,rho,h,qeff
       REAL(KIND=dp) :: Pressure, DampCoef, WinklerCoef
       LOGICAL :: Stat
       INTEGER :: i,j,p,q,t
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------

       Force = 0.0d0
       STIFF = 0.0d0
       DAMP = 0.0d0
       MASS = 0.0d0
       Curvature = 0.0d0
       ShearStrain = 0.0d0
!
!      Numerical integration:
!      ----------------------

       IntegStuff = GaussPoints(Element,3)

       DO t = 1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)

!
!        Basis function values & derivatives at the integration point
!        -------------------------------------------------------------
         stat = ElementInfo(Element,Nodes,U,V,W,detJ,Basis,dBasisdx)

         S = S * detJ

         Pressure = SUM( Load(1:n)*Basis(1:n) )
         h = SUM( Thickness(1:n)*Basis(1:n) )
         DampCoef = SUM( DampingCoef(1:n) * Basis(1:n) )
         WinklerCoef = SUM( SpringCoef(1:N) * Basis(1:N) )

         IF(HoleCorrection .AND. GotHoleType) THEN
           CALL PerforatedElasticity(Ematrix,Gmatrix, &
               Poisson,Young,Thickness,HoleFraction,HoleSize,Basis,n)
           qeff = SUM(HoleFraction(1:n)*Basis(1:n))
           rho = (1.0d0-qeff) * SUM( Density(1:n)*Basis(1:n) )
           Tmatrix = 0.0d0
           Tmatrix(1,1) = SQRT(1.0d0-qeff**2) * SUM( Tension(1:n)*Basis(1:n) ) * h
           Tmatrix(2,2) = Tmatrix(1,1)
         ELSE
           CALL IsotropicElasticity(Ematrix,Gmatrix, &
               Poisson,Young,Thickness,Basis,n)
           rho = SUM( Density(1:n)*Basis(1:n) )
           Tmatrix = 0.0d0
           Tmatrix(1,1) = SUM( Tension(1:n)*Basis(1:n) ) * h
           Tmatrix(2,2) = SUM( Tension(1:n)*Basis(1:n) ) * h
         END IF


!
!        The degrees-of-freedom are  (u_x, u_y, u_z, r_x, r_y, r_z)
!        where u_i is the displacement [m] and r_i the rotation [1]
!        ----------------------------------------------------------

!        Bending stiffness:
!        ------------------
         DO p=1,n
            Curvature(1,3*p-1) = dBasisdx(p,1)
            Curvature(2,3*p  ) = dBasisdx(p,2)
            Curvature(3,3*p-1) = dBasisdx(p,2)
            Curvature(3,3*p  ) = dBasisdx(p,1)
         END DO

         CALL AddInnerProducts(STIFF,Ematrix,Curvature,3,3*n,s)

!        In-plane stiffness:
!        -------------------

!
!        Shear stiffness:
!        ----------------
         CALL CovariantInterpolation(ShearStrain, &
              Basis, Nodes % x(1:n),Nodes % y(1:n),U,V,n)

         CALL ShearCorrectionFactor(Kappa, h, &
              Nodes % x(1:n), Nodes % y(1:n), n)

         DO p=1,n
            ShearStrain(1,3*p-2) = dBasisdx(p,1)
            ShearStrain(2,3*p-2) = dBasisdx(p,2)
         END DO 

         CALL AddInnerProducts(STIFF, &
              Gmatrix,ShearStrain,2,3*n,Kappa*s)
!
!        Tensile stiffness:
!        ------------------
         ShearStrain = 0.0d0
         DO p=1,n
            ShearStrain(1,3*p-2) = dBasisdx(p,1)
            ShearStrain(2,3*p-2) = dBasisdx(p,2)
         END DO 

         CALL AddInnerProducts(STIFF,Tmatrix,ShearStrain,2,3*n,s)

!        Spring Coeffficient:
!        -------------------
         DO p = 1,n
            i = DOFs*(p-1)+1
            DO q = 1,n
               j = DOFs*(q-1)+1
               STIFF(i,j) = STIFF(i,j) &
                    + WinklerCoef * Basis(p) * Basis(q) * s
            END DO
         END DO

!
!        Load vector:
!        ------------
         DO p=1,n
            i = DOFs*(p-1)+1
            Force(i) = Force(i) + Pressure * Basis(p) * s
         END DO
!
!        Mass matrix:
!        ------------
         DO p = 1,n
            i = DOFs*(p-1)+1
            DO q = 1,n
               j = DOFs*(q-1)+1
               MASS(i,j) = MASS(i,j) + rho * h * Basis(p) * Basis(q) * s
            END DO
         END DO

!
!        Damping matrix:
!        ---------------
         DO p = 1,n
            i = DOFs*(p-1)+1
            DO q = 1,n
               j = DOFs*(q-1)+1
               DAMP(i,j) = DAMP(i,j) + DampCoef * Basis(p) * Basis(q) * s
            END DO
         END DO

!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix


!==============================================================================


     SUBROUTINE IsotropicElasticity(Ematrix, &
          Gmatrix,Poisson,Young,Thickness,Basis,n)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Ematrix(:,:), Gmatrix(:,:), Basis(:)
     REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:)
     REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw
     INTEGER :: n
!------------------------------------------------------------------------------
       Euvw = SUM( Young(1:n)*Basis(1:n) )
       Puvw = SUM( Poisson(1:n)*Basis(1:n) )
       Tuvw = SUM( Thickness(1:n)*Basis(1:n) )
       Guvw = Euvw/(2.0d0*(1.0d0 + Puvw))

       Ematrix = 0.0d0
       Ematrix(1,1) = 1.0d0
       Ematrix(1,2) = Puvw
       Ematrix(2,1) = Puvw
       Ematrix(2,2) = 1.0d0
       Ematrix(3,3) = (1.0d0-Puvw)/2.0d0

       Ematrix = Ematrix* Euvw * (Tuvw**3) / (12.0d0*(1.0d0-Puvw**2))

       Gmatrix = 0.0d0
       Gmatrix(1,1) = Guvw*Tuvw
       Gmatrix(2,2) = Guvw*Tuvw
!------------------------------------------------------------------------------
     END SUBROUTINE IsotropicElasticity


!------------------------------------------------------------------------------
!> The elastic model for perforated plates is taken directly from
!> M. Pedersen, W. Olthuis, P. BergWald:
!> 'On the mechanical behavior of thin perforated plates and their application
!>  in silicon condenser microphones', Sensors and Actuators A 54 (1996) 499-504.
! The model in verified in the special assignment of Jani Paavilainen 
!------------------------------------------------------------------------------

     SUBROUTINE PerforatedElasticity(Ematrix, &
         Gmatrix,Poisson,Young,Thickness,HoleFraction, &
         HoleSize, Basis,n)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Ematrix(:,:), Gmatrix(:,:), Basis(:), HoleSize(:)
     REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:), HoleFraction(:)
     REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw, q, a, b, k, sq

     INTEGER :: n
!------------------------------------------------------------------------------
       Euvw = SUM( Young(1:n)*Basis(1:n) )
       Puvw = SUM( Poisson(1:n)*Basis(1:n) )
       Tuvw = SUM( Thickness(1:n)*Basis(1:n) )

       q = SUM( HoleFraction(1:n)*Basis(1:n) )
       a = SUM( HoleSize(1:n)*Basis(1:n))
       sq = SQRT(q)
       b = 2*a/sq

       IF(Tuvw > b-2*a) THEN
         k = (Tuvw-0.63_dp*(b-2*a)) * (b-2*a)**3 / 3.0d0
       ELSE
         k = ((b-2*a)-0.63_dp*Tuvw) * Tuvw**3 / 3.0d0
       END IF

       Ematrix = 0.0d0
       Ematrix(1,1) = (1.0d0-sq)/(1.0d0-Puvw**2) + 0.5d0*sq * (1.0d0-sq)**2
       Ematrix(2,2) = Ematrix(1,1)
       Ematrix(1,2) = Puvw * (1.0d0-sq)/(1.0-Puvw**2)
       Ematrix(2,1) = Ematrix(1,2)
       Ematrix(3,3) = 0.5d0*(1.0d0-sq)/(1.0d0+Puvw) + &
           1.5d0*k*sq*(1-sq)/(b*(1.0d0+Puvw)*Tuvw**3)

       Guvw = Euvw * Ematrix(3,3) ! * 2.0d0? 

       Ematrix = Ematrix * Euvw * (Tuvw**3) / 12.0d0

       Gmatrix = 0.0d0
       Gmatrix(1,1) = Guvw * Tuvw
       Gmatrix(2,2) = Gmatrix(1,1)
!------------------------------------------------------------------------------
     END SUBROUTINE PerforatedElasticity

!==============================================================================


     SUBROUTINE ShearCorrectionFactor(Kappa,Thickness,x,y,n)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Kappa,Thickness,x(:),y(:)
       INTEGER :: n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: x21,x32,x43,x13,x14,y21,y32,y43,y13,y14, &
            l21,l32,l43,l13,l14,alpha,h
!------------------------------------------------------------------------------
       Kappa = 1.0d0
       SELECT CASE(n)
          CASE(3)
             alpha = 0.20d0
             x21 = x(2)-x(1)
             x32 = x(3)-x(2)
             x13 = x(1)-x(1)
             y21 = y(2)-y(1)
             y32 = y(3)-y(2)
             y13 = y(1)-y(1)
             l21 = SQRT(x21**2 + y21**2)
             l32 = SQRT(x32**2 + y32**2)
             l13 = SQRT(x13**2 + y13**2)
             h = MAX(l21,l32,l13)
             Kappa = (Thickness**2)/(Thickness**2 + alpha*(h**2))
          CASE(4)
             alpha = 0.10d0
             x21 = x(2)-x(1)
             x32 = x(3)-x(2)
             x43 = x(4)-x(3)
             x14 = x(1)-x(4)
             y21 = y(2)-y(1)
             y32 = y(3)-y(2)
             y43 = y(4)-y(3)
             y14 = y(1)-y(4)
             l21 = SQRT(x21**2 + y21**2)
             l32 = SQRT(x32**2 + y32**2)
             l43 = SQRT(x43**2 + y43**2)
             l14 = SQRT(x14**2 + y14**2)
             h = MAX(l21,l32,l43,l14)
             Kappa = (Thickness**2)/(Thickness**2 + alpha*(h**2))
          CASE DEFAULT
            CALL WARN('SmitcSolver','Illegal number of nodes for Smitc elements')
          END SELECT
!------------------------------------------------------------------------------
     END SUBROUTINE ShearCorrectionFactor


!==============================================================================


     SUBROUTINE AddInnerProducts(A,B,C,m,n,s)
!------------------------------------------------------------------------------
!      Performs the operation
!
!         A = A + C' * B * C * s
!
!      with
!
!         Size( A ) = n x n
!         Size( B ) = m x m
!         Size( C ) = m x n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: A(:,:),B(:,:),C(:,:),s
       INTEGER :: m,n
!------------------------------------------------------------------------------
       INTEGER :: i,j,k,l
!------------------------------------------------------------------------------
       DO i=1,n
          DO j=1,n
             DO k=1,m
                DO l=1,m
                   A(i,j) = A(i,j) + C(k,i)*B(k,l)*C(l,j) * s
                END DO
             END DO
          END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE AddInnerProducts


!==============================================================================


     SUBROUTINE CovariantInterpolation(ShearStrain,Basis,X,Y,U,V,n)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: ShearStrain(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: detJ,Jmat(2,2),invJ(2,2),ShearRef(2,100)
       REAL(KIND=dp) :: Tau(2),Sdofs(100)
       INTEGER :: j

       SELECT CASE(n)

!      The SMITC3 element
!      ==================
       CASE(3)
          CALL Jacobi3(Jmat,invJ,detJ,x,y)
          ShearRef = 0.0d0
          ShearStrain = 0.0d0

!         Compute the shear-dofs for edge 12:
!         ===================================
          Tau(1) = 1.0d0
          Tau(2) = 0.0d0
          Tau = (/ 1.0d0, 0.0d0/)

          Sdofs = 0.0d0
          Sdofs(2) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(3) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
          Sdofs(5) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(6) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
           
          DO j = 1,9
             ShearRef(1,j) = ShearRef(1,j) + (1.0d0-V)*Sdofs(j)
             ShearRef(2,j) = ShearRef(2,j) + (U)*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 23:
!         ===================================
          Tau(1) = -1.0d0/SQRT(2.0d0)
          Tau(2) =  1.0d0/SQRT(2.0d0)

          Sdofs = 0.0d0
          Sdofs(5) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/SQRT(2.0d0)
          Sdofs(6) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/SQRT(2.0d0)
          Sdofs(8) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/SQRT(2.0d0)
          Sdofs(9) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/SQRT(2.0d0)

          DO j = 1,9
             ShearRef(1,j) = ShearRef(1,j) + (-V)*Sdofs(j)
             ShearRef(2,j) = ShearRef(2,j) + (U)*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 31:
!         ===================================
          Tau(1) =  0.0d0
          Tau(2) = -1.0d0

          Sdofs = 0.0d0
          Sdofs(2) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(3) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
          Sdofs(8) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(9) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0

          DO j = 1,9
             ShearRef(1,j) = ShearRef(1,j) + ( -V )*Sdofs(j)
             ShearRef(2,j) = ShearRef(2,j) + (-1.0d0+U)*Sdofs(j)
          END DO

!         Compute the final reduced shear strain
!         ======================================
          ShearStrain(1:2,1:9) = MATMUL(invJ,ShearRef(1:2,1:9))


!      The SMITC4 element
!      ==================
       CASE(4)
          ShearRef = 0.0d0
          ShearStrain = 0.0d0

!         Compute the shear-dofs for edge 12:
!         ===================================
          Tau(1) = 1.0d0
          Tau(2) = 0.0d0

          CALL Jacobi4(Jmat,invJ,detJ,0.0d0,-1.0d0,x,y)
          
          Sdofs = 0.0d0
          Sdofs(2) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(3) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(5) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(6) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,12
             ShearRef(1,j) = ShearRef(1,j) + (1-V)/4.0d0*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 23:
!         ===================================
          Tau(1) = 0.0d0
          Tau(2) = 1.0d0

          CALL Jacobi4(Jmat,invJ,detJ,1.0d0,0.0d0,x,y)

          Sdofs = 0.0d0
          Sdofs(5) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(6) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(8) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(9) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,12
             ShearRef(2,j) = ShearRef(2,j) + (1+U)/4.0d0*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 34:
!         ===================================
          Tau(1) = -1.0d0
          Tau(2) =  0.0d0

          CALL Jacobi4(Jmat,invJ,detJ,0.0d0,1.0d0,x,y)

          Sdofs = 0.0d0
          Sdofs(8)  = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(9)  = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(11) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(12) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,12
             ShearRef(1,j) = ShearRef(1,j) + (-1-V)/4.0d0*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 41:
!         ===================================
          Tau(1) =  0.0d0
          Tau(2) = -1.0d0

          CALL Jacobi4(Jmat,invJ,detJ,-1.0d0,0.0d0,x,y)

          Sdofs = 0.0d0
          Sdofs(2)  = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(3)  = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(11) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(12) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,12
             ShearRef(2,j) = ShearRef(2,j) + (-1+U)/4.0d0*Sdofs(j)
          END DO

!         Compute the final reduced shear strain
!         ======================================
          CALL Jacobi4(Jmat,invJ,detJ,U,V,x,y)
          ShearStrain(1:2,1:12) = MATMUL(invJ,ShearRef(1:2,1:12))

       CASE DEFAULT
         CALL WARN('SmitcSolver','Illegal number of nodes for Smitc elements.')

       END SELECT
!------------------------------------------------------------------------------
     END SUBROUTINE CovariantInterpolation


!==============================================================================


     SUBROUTINE Jacobi3(Jmat,invJ,detJ,x,y)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Jmat(:,:),invJ(:,:),detJ,x(:),y(:)
!------------------------------------------------------------------------------
       Jmat(1,1) = x(2)-x(1)
       Jmat(2,1) = x(3)-x(1)
       Jmat(1,2) = y(2)-y(1)
       Jmat(2,2) = y(3)-y(1)

       detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

       invJ(1,1) =  Jmat(2,2)/detJ
       invJ(2,2) =  Jmat(1,1)/detJ
       invJ(1,2) = -Jmat(1,2)/detJ
       invJ(2,1) = -Jmat(2,1)/detJ
!------------------------------------------------------------------------------
     END SUBROUTINE Jacobi3


!==============================================================================


     SUBROUTINE Jacobi4(Jmat,invJ,detJ,xi,eta,x,y)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Jmat(:,:),invJ(:,:),detJ,xi,eta,x(:),y(:)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: dNdxi(4), dNdeta(4)
       INTEGER :: i

       dNdxi(1) = -(1-eta)/4.0d0
       dNdxi(2) =  (1-eta)/4.0d0
       dNdxi(3) =  (1+eta)/4.0d0
       dNdxi(4) = -(1+eta)/4.0d0
       dNdeta(1) = -(1-xi)/4.0d0
       dNdeta(2) = -(1+xi)/4.0d0
       dNdeta(3) =  (1+xi)/4.0d0
       dNdeta(4) =  (1-xi)/4.0d0
       
       Jmat = 0.0d0
       DO i=1,4
          Jmat(1,1) = Jmat(1,1) + dNdxi(i)*x(i)
          Jmat(1,2) = Jmat(1,2) + dNdxi(i)*y(i)
          Jmat(2,1) = Jmat(2,1) + dNdeta(i)*x(i)
          Jmat(2,2) = Jmat(2,2) + dNdeta(i)*y(i)
       END DO

       detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

       invJ(1,1) = Jmat(2,2)/detJ
       invJ(2,2) = Jmat(1,1)/detJ
       invJ(1,2) = -Jmat(1,2)/detJ
       invJ(2,1) = -Jmat(2,1)/detJ
!------------------------------------------------------------------------------
     END SUBROUTINE Jacobi4

!==============================================================================


   END SUBROUTINE SmitcSolver
!------------------------------------------------------------------------------
