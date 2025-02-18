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
! * Module for the solution of reduced dimensional Navier-Stokes equation.
! * Here we assume that the average velocity of the film defines the full velocity
! * profile and hence the depth direction may be eliminated from the flow solution.
! * The intended use is film/channel flow but also cylindrical pipe flow is implemented.
! * This fills the gap between full Navier-Stokes solver and reduced dimensional
! * Reynolds solver. 
! *
! * The module is compatible with p-bubbles and/or p2/p1 elements, e.g.
! * 303b1    - triangle with one bubble
! * 303e1b1  - p2/p1 triangle!
! * 404b4    - quad with four bubbles
! * 404e1b1  - q2/q1 quad
! *
! * This module has been defived from a historical Navier-Stokes solver and
! * updated much later for problems involving challel flows.
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 27.10.2022
! *
!/*****************************************************************************/


!------------------------------------------------------------------------------
SUBROUTINE FilmFlowSolver_init0( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  LOGICAL :: Found, Serendipity
  TYPE(ValueList_t), POINTER :: Params 
  Params => GetSolverParams()

  Serendipity = GetLogical( GetSimulation(), 'Serendipity P Elements', Found)
  IF(.NOT.Found) Serendipity = .TRUE.

  CALL ListAddNewInteger(Params, 'Time derivative Order', 1)
  
  IF(Serendipity) THEN
    CALL ListAddNewString(Params,'Element','p:1 -line b:1 -tri b:1 -quad b:3')
  ELSE
    CALL ListAddNewString(Params,'Element','p:1 -line b:1 -tri b:1 -quad b:4')
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE FilmFlowSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE FilmFlowSolver_init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params 
  LOGICAL :: Found
  INTEGER :: mdim
  CHARACTER(*), PARAMETER :: Caller = 'FilmFlowSolver_init'
!------------------------------------------------------------------------------ 
  Params => GetSolverParams() 
  
  mdim = ListGetInteger(Params,'Model Dimension',Found)
  IF(.NOT. Found ) THEN
    CALL Fatal(Caller,'Give "Model Dimension" i.e. the dimension of N-S equation!')    
  END IF
   
  IF ( mdim == 2 ) THEN
    CALL ListAddNewString(Params, 'Variable', &
        'Flow[FilmVelocity:2 FilmPressure:1]')
  ELSE IF( mdim == 1 ) THEN
    CALL ListAddNewString(Params, 'Variable', &
        'Flow[FilmSpeed:1 FilmPressure:1]')
  ELSE
    CALL Fatal(Caller,'This module does not make sense in dim: '//I2S(mdim))    
  END IF

  ! Study only velocity components in linear system
  CALL ListAddNewInteger(Params, 'Nonlinear System Norm DOFs', mdim )

  ! This should be true to incompressible flows where pressure level is not uniquely determined
  CALL ListAddNewLogical(Params, 'Relative Pressure Relaxation', .TRUE. )

  ! It makes sense to eliminate the bubbles to save memory and time
  CALL ListAddNewLogical(Params, 'Bubbles in Global System', .FALSE.)

  ! Use global mass matrix in time integration
  CALL ListAddNewLogical(Params, 'Global Mass Matrix', .TRUE.)

!------------------------------------------------------------------------------ 
END SUBROUTINE FilmFlowSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE FilmFlowSolver( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Found, Convect, CSymmetry
  TYPE(Element_t),POINTER :: Element
  INTEGER :: i,n, nb, nd, t, istat, dim, mdim, BDOFs=1,Active,iter,maxiter,CoupledIter
  REAL(KIND=dp) :: Norm = 0, mingap
  TYPE(ValueList_t), POINTER :: Params, BodyForce, Material, BC
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), &
      FORCE(:), rho(:), gap(:), gap0(:), mu(:), ac(:), Pres(:), Velocity(:,:), MASS(:,:),&
      PrevPressure(:), FsiRhs(:,:), PrevGap(:)
  LOGICAL :: GradP, LateralStrain, GotAc, SurfAc, UsePrevGap
  TYPE(Variable_t), POINTER :: pVar, thisVar
  INTEGER :: GapDirection
  REAL(KIND=dp) :: GapFactor
  CHARACTER(*), PARAMETER :: Caller = 'FilmFlowSolver'
  LOGICAL :: Debug
  
  SAVE STIFF, MASS, LOAD, FORCE, rho, ac, gap, gap0, mu, Pres, Velocity, &
      PrevPressure, AllocationsDone, pVar, GotAc, SurfAC, FsiRhs, PrevGap
!------------------------------------------------------------------------------

  CALL Info(Caller,'Computing reduced dimensional Navier-Stokes equations!')

  Mesh => GetMesh()
  Element => GetActiveElement(1)

  dim = CoordinateSystemDimension()
    
  Params => GetSolverParams()
  thisVar => Solver % Variable
  
  mdim = ListGetInteger( Params,'Model Dimension',UnFoundFatal=.TRUE.)
  Convect = GetLogical( Params, 'Convect', Found )
  GradP = GetLogical( Params, 'GradP Discretization', Found ) 
  LateralStrain = GetLogical( Params,'Lateral Strain',Found )
  mingap = ListGetCReal( Params,'Min Gap Height',Found )
  IF(.NOT. Found) mingap = TINY(mingap)
  GotAC = ListCheckPresentAnyMaterial( Model,'Artificial Compressibility')

  UsePrevGap = ListGetLogical( Params,'Use Gap Average',Found )
  
  CoupledIter = GetCoupledIter()
  
  GapDirection = 0
  GapFactor = ListGetCReal( Params,'Gap Addition Factor',Found )
  IF( Found ) THEN  
    GapDirection = mdim+1
    IF( ABS( GapFactor ) > 1.0_dp ) THEN
      CALL Warn(Caller,'"Gap Addition Factor" greater to unity does not make sense!')
    END IF
  END IF
  
  CSymmetry = ListGetLogical( Params,'Axi Symmetric',Found )
  IF(.NOT. Found ) THEN 
    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) 
  END IF

  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    CALL Info(Caller,'Dimension of Navier-Stokes equation: '//I2S(mdim))
    CALL Info(Caller,'Dimension of coordinate system: '//I2S(dim))

    n = (mdim+1)*(Mesh % MaxElementDOFs+4*BDOFs)  ! just big enough for elemental arrays
    ALLOCATE( FORCE(n), LOAD(mdim+2,n), STIFF(n,n), MASS(n,n), &
        rho(n), ac(n), gap(n), gap0(n), mu(n), Pres(n), Velocity(mdim+1,n), STAT=istat )
    Velocity = 0.0_dp
    IF ( istat /= 0 ) THEN
      CALL Fatal( Caller, 'Memory allocation error.' )
    END IF
    IF( GradP ) THEN
      CALL Info(Caller,'"Gradp Discretization" is set True',Level=10)
    END IF     

    pVar => VariableGet( Mesh % Variables,'FilmPressure')
    IF ( .NOT. ASSOCIATED(pVar) ) THEN
      CALL Fatal( Caller, 'Could not find required field "FilmPressure"!')
    END IF
    n = SIZE(pVar % Values) 
    IF( GotAC ) THEN
      CALL Info(Caller,'Using artificial compressibility for FSI emulation!') 
      ALLOCATE(PrevPressure(n),FsiRhs(2,n))
      PrevPressure = 0.0_dp
      FsiRhs = 0.0_dp
      SurfAc = ListGetLogical( Params,'Surface Compressibility',Found )
    END IF
    IF(UsePrevGap) THEN
      ALLOCATE(PrevGap(n))
      PrevGap = 0.0_dp
    END IF
    AllocationsDone = .TRUE.
  END IF

  IF(GotAc) THEN  
    PrevPressure = pVar % Values
    FsiRhs = 0.0_dp
  END IF
   
  maxiter = ListGetInteger( Params,'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  DO iter=1,maxiter
    
    !Initialize the system and do the assembly:
    !----------------
    CALL DefaultInitialize()

    Newton = GetNewtonActive()
    
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      ! Volume forces:
      !---------------
      BodyForce => GetBodyForce()
      LOAD = 0.0d0
      IF ( ASSOCIATED(BodyForce) ) THEN
        Load(1,1:n) = GetReal( BodyForce, 'Flow Bodyforce 1', Found )
        IF(mdim==2) Load(2,1:n) = GetReal( BodyForce, 'Flow Bodyforce 2', Found )
        Load(mdim+1,1:n) = GetReal( BodyForce, 'Normal Velocity', Found )
        Load(mdim+2,1:n) = GetReal( BodyForce, 'Fsi Velocity', Found )
      END IF

      ! Material parameters:
      !---------------------
      Material => GetMaterial()
      rho(1:n) = GetReal( Material, 'Density' )
      mu(1:n)  = GetReal( Material, 'Viscosity' )
      gap(1:n) = GetReal( Material, 'Gap Height' )

      WHERE(gap(1:n) < mingap )
        gap(1:n) = mingap
      END WHERE

      IF(UsePrevGap) THEN
        IF(CoupledIter == 1) THEN
          PrevGap(pVar % Perm(Element % NodeIndexes)) = gap(1:n)
        END IF
        gap0(1:n) = PrevGap(pVar % Perm(Element % NodeIndexes))
      END IF
        
      
      IF( GotAC ) THEN
        ac(1:n) = GetReal( Material,'Artificial Compressibility',Found )
        Pres(1:n) = PrevPressure(pVar % Perm(Element % NodeIndexes))
      END IF
            
      ! Get previous elementwise velocity iterate:
      ! Note: pressure is the dim+1 component here!
      !-------------------------------------------
      IF ( Convect ) CALL GetVectorLocalSolution( Velocity )

      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(  MASS, STIFF, FORCE, LOAD, rho, gap, gap0, mu, &
          ac, Velocity, Pres, Element, n, nd, nd+nb, &
          dim, mdim )
      
      IF ( nb>0 ) THEN
        CALL LCondensate( nd, nb, mdim, STIFF, FORCE )
      END IF

      IF ( Transient ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF
      
      ! Update global matrix and rhs vector from local matrix & vector:
      !----------------------------------------------------------------
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
    CALL DefaultFinishBulkAssembly()


    IF( GotAC ) THEN
      BLOCK
        REAL(KIND=dp) :: sorig, sfsi, coeff
        sorig = SUM(FsiRhs(1,:))
        sfsi = SUM(FsiRhs(2,:))
        coeff = 1.0_dp
        IF(sfsi /= 0.0) coeff = sorig / sfsi            
        IF(sfsi < sorig) coeff = 1.0_dp

        ! Just report incoming and outgoing total fluxes
        PRINT *,'RHSComp:',sorig,sfsi,coeff
      END BLOCK
    END IF
      
    DO t=1, Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE

      n = GetElementNOFNodes()      
      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE

      rho(1:n) = GetParentMatProp( 'Density', Element, Found )
      mu(1:n)  = GetParentMatProp( 'Viscosity', Element, Found )
      gap(1:n) = GetParentMatProp( 'Gap Height', Element, Found )

      WHERE(gap(1:n) < mingap )
        gap(1:n) = mingap
      END WHERE

      IF(UsePrevGap) THEN
        gap0(1:n) = PrevGap(pVar % Perm(Element % NodeIndexes))
      END IF
      
      DO i=1,mdim
        Load(i,1:n) = GetReal( BC, 'Pressure '//I2S(i), Found ) 
      END DO
      Load(mdim+1,1:n) = GetReal( BC, 'Mass Flux', Found )
      
      CALL LocalBoundaryMatrix(  MASS, STIFF, FORCE, Load, rho, gap, mu, &
          Element, n, dim, mdim )

      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
    
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    
    Norm = DefaultSolve()
    
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT    
  END DO


  CALL DefaultFinish()

  BLOCK
    REAL(KIND=dp), POINTER :: Comp(:)

    IF( InfoActive(10) ) THEN
      n = SIZE(pVar % Values)
      DO i=1, Solver % Variable % dofs
        Comp => Solver % Variable % Values(i::Solver % Variable % Dofs)      
        CALL VectorValuesRange(Comp,n,'Velocity '//I2S(i))       
      END DO
      CALL VectorValuesRange(pVar % Values,n,'Pressure')       
      IF(GotAc) THEN
        CALL VectorValuesRange(PrevPressure,n,'Pressure0')       
        CALL VectorValuesRange(pVar % Values - PrevPressure,n,'PressureDiff '//I2S(CoupledIter))       
      END IF
    END IF      
  END BLOCK
  
  CALL Info(Caller,'All done',Level=12)

  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(  MASS, STIFF, FORCE, LOAD, Nodalrho, NodalGap, &
      NodalGap0, Nodalmu, NodalAC, NodalVelo, NodalPres, Element, n, nd, &
      ntot, dim, mdim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), NodalAC(:), Nodalrho(:), &
        NodalGap(:), NodalGap0(:), NodalPres(:), NodalVelo(:,:)
    INTEGER :: dim, mdim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3)
    REAL(KIND=dp) :: DetJ,LoadAtIP(mdim+2),Velo(mdim), VeloGrad(mdim,mdim), gapGrad(mdim) 
    REAL(KIND=dp), POINTER :: A(:,:),F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q, geomc
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, pres, gap, gap0, gap2, ac, s, s0, s1, MinPres
    LOGICAL :: Visited = .FALSE.
    
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes, Visited, MinPres
!------------------------------------------------------------------------------

    
    CALL GetElementNodes( Nodes )

    IF( GapDirection > 0 ) THEN
      SELECT CASE( GapDirection )
      CASE(1)
        Nodes % x(1:n) = Nodes % x(1:n) + GapFactor * NodalGap(1:n)
      CASE(2)
        Nodes % y(1:n) = Nodes % y(1:n) + GapFactor * NodalGap(1:n)
      CASE(3)
        Nodes % z(1:n) = Nodes % z(1:n) + GapFactor * NodalGap(1:n)
      END SELECT
      ! Does this have an effect?
      !PRINT *,'GapFactor:',GapFactor * NodalGap(1:n)
    END IF

    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0
    gapGrad = 0.0_dp
        
    ! To my understanding we want to include the gap height to weight
    IF( Csymmetry ) THEN
      geomc = 2
    ELSE
      geomc = 1
    END IF
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element, PReferenceElement = .TRUE. )

    !IP = GaussPoints( Element )

    IF(.NOT. Visited) THEN
      CALL Info(Caller,'Number of integration points: '//I2S(IP % n))
      MinPres = ListGetConstReal( Params,'Min FilmPressure',Found )
      IF(.NOT. Found) MinPres = -HUGE(MinPres)
      Visited = .TRUE.
    END IF
    
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
         IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       s1 = s 
       s0 = s 
       
       
       ! Material parameters at the integration point:
       !----------------------------------------------      
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       gap = SUM( Basis(1:n) * NodalGap(1:n) ) 
       gap0 = SUM( Basis(1:n) * NodalGap0(1:n) ) 
       
       DO i=1,mdim
         gapGrad(i) = SUM( NodalGap(1:nd) * dBasisdx(1:nd,i) )
       END DO
       
       ! Previous velocity at the integration point:
       !--------------------------------------------
       Velo = MATMUL( NodalVelo(1:mdim,1:nd), Basis(1:nd) )
       VeloGrad = MATMUL( NodalVelo(1:mdim,1:nd), dBasisdx(1:nd,1:mdim) )
       
       IF( GotAC ) THEN
         Pres = SUM( NodalPres(1:n) * Basis(1:n) )
         Pres = MAX(MinPres,Pres)
         ac = SUM( NodalAC(1:n) * Basis(1:n) ) / dt
         !IF(.NOT. SurfAc) ac = ac * gap
       END IF
       
       ! The source term at the integration point:
       !------------------------------------------
       DO i=1,mdim+2
         LoadAtIP(i) = SUM( Basis(1:n) * LOAD(i,1:n) )
       END DO

       IF ( Convect .AND. Newton ) THEN
         LoadAtIp(1:mdim) = LoadAtIp(1:mdim) + rho * MATMUL(VeloGrad(1:mdim,1:mdim),Velo(1:mdim))
       END IF
       LoadAtIp(mdim+1) = geomc * LoadAtIp(mdim+1) 
       ! Fsi velocity
       LoadAtIp(mdim+2) = geomc * LoadAtIp(mdim+2) 
         
       ! Finally, the elemental matrix & vector:
       !----------------------------------------       
       DO p=1,ntot
         DO q=1,ntot
           i = (mdim+1) * (p-1) + 1
           j = (mdim+1) * (q-1) + 1
           A => STIFF(i:i+mdim,j:j+mdim)
           M => MASS(i:i+mdim,j:j+mdim)

           DO i=1,mdim
             IF( Transient ) THEN
               M(i,i) = M(i,i) + s * rho * Basis(q) * Basis(p)
             END IF

             DO j = 1,mdim
               IF( LateralStrain ) THEN 
                 A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                 A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
               END IF
                 
               IF ( Convect ) THEN
                 A(i,i) = A(i,i) + s * rho * Velo(j) * dBasisdx(q,j) * Basis(p)
                 IF ( Newton ) THEN
                    A(i,j) = A(i,j) + s * rho * VeloGrad(i,j) * Basis(q) * Basis(p)
                 END IF
               END IF
             END DO
             
             ! This is the Poisseille flow resistance
             IF(UsePrevGap) THEN
               ! This takes the analytical average when going from 1/d_0^2 to 1/d^2. 
               gap2 = gap*gap0
             ELSE
               gap2 = gap**2
             END IF

             IF( CSymmetry ) THEN
               A(i,i) = A(i,i) + s * ( 8 * mu / gap2 ) * Basis(q) * Basis(p)  
             ELSE
               A(i,i) = A(i,i) + s * ( 12 * mu / gap2 ) * Basis(q) * Basis(p)  
             END IF
               
             ! Note that here the gap height must be included in the continuity equation
             IF( GradP ) THEN
               A(i,mdim+1) = A(i,mdim+1) + s * dBasisdx(q,i) * Basis(p)
               A(mdim+1,i) = A(mdim+1,i) - s * gap * rho * Basis(q) * dBasisdx(p,i)               
             ELSE
               A(i,mdim+1) = A(i,mdim+1) - s * Basis(q) * dBasisdx(p,i)
               A(mdim+1,i) = A(mdim+1,i) + s * gap * rho * dBasisdx(q,i) * Basis(p) & 
                   + geomc * s * rho * Basis(q) * gapGrad(i) * Basis(p)
             END IF
           END DO
             
           ! This is the implicit term in artificial compressibility for FSI coupling
           ! applied to thin film flow. 
           ! Div(u) + (c/dt)*p^(m) = (c/dt)*p^(m-1)
           ! See Raback et al., CFD Eccomas 2001.
           ! "FLUID-STRUCTURE INTERACTION BOUNDARY CONDITIONS BY ARTIFICIAL COMPRESSIBILITY".
           IF(GotAC) A(mdim+1,mdim+1) = A(mdim+1,mdim+1) + ac * s * rho * Basis(q) * Basis(p)              
         END DO
         
         i = (mdim+1) * (p-1) + 1
         F => FORCE(i:i+mdim)
         
         ! This is the explit term in artificial compressibility for FSI coupling
         IF( GotAC ) F(mdim+1) = F(mdim+1) + ac * s * rho * Basis(p) * Pres         

         ! Body force for velocity components and pressure
         F(1:mdim+1) = F(1:mdim+1) + s * rho * Basis(p) * LoadAtIp(1:mdim+1)

         ! Additional body force from FSI velocity
         F(mdim+1) = F(mdim+1) - s * rho * Basis(p) * LoadAtIp(mdim+2) 
       END DO

       ! These are just recorded in order to study the total forced
       ! and induced (by FSI coupling) fluxes. 
       IF(GotAC) THEN
         FsiRhs(1,ThisVar % Perm(Element % NodeIndexes)) = &
             FsiRhs(1,ThisVar % Perm(Element % NodeIndexes))  + &
             s * rho * LoadAtIp(mdim+1) * Basis(1:n)             
         
         FsiRhs(2,ThisVar % Perm(Element % NodeIndexes)) = &
             FsiRhs(2,ThisVar % Perm(Element % NodeIndexes))  + &
             s * rho * LoadAtIp(mdim+2) * Basis(1:n)             
       END IF
     END DO
     
   ! for p2/p1 elements set Dirichlet constraint for unused dofs,
   ! EliminateDirichlet will get rid of these:
   !-------------------------------------------------------------
    DO p = n+1,ntot
      i = (mdim+1) * p
      FORCE(i)   = 0.0d0
      MASS(:,i)  = 0.0d0
      MASS(i,:)  = 0.0d0
      STIFF(i,:) = 0.0d0
      STIFF(:,i) = 0.0d0
      STIFF(i,i) = 1.0d0
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalBulkMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalBoundaryMatrix(  MASS, STIFF, FORCE, NodalLoad, Nodalrho, NodalGap, &
      Nodalmu, Element, n, dim, mdim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), NodalLoad(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), NodalGap(:)
    INTEGER :: dim, mdim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Load(mdim+1)
    REAL(KIND=dp), POINTER :: F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, p, q, geomc
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, gap, ac, s
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0
    
    ! To my understanding we want to include the gap height to weight
    IF( Csymmetry ) THEN
      geomc = 2
    ELSE
      geomc = 1
    END IF
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
         IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ
              
       ! Material parameters at the integration point:
       !----------------------------------------------      
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       gap = SUM( Basis(1:n) * NodalGap(1:n) ) 

       DO i=1,mdim+1
         Load(i) = SUM( Basis(1:n) * NodalLoad(i,1:n) ) 
       END DO
         
       ! Finally, the elemental matrix & vector:
       !----------------------------------------       
       DO p=1,n
         i = (mdim+1) * (p-1) + 1
         F => FORCE(i:i+mdim)
         F = F + s * Basis(p) * Load
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalBoundaryMatrix
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
    SUBROUTINE LCondensate( N, nb, dim, K, F )
!------------------------------------------------------------------------------
      USE LinearAlgebra
      INTEGER :: N, nb, dim
      REAL(KIND=dp) :: K(:,:),F(:), Kbb(Nb*dim,Nb*dim), &
       Kbl(nb*dim,n*(dim+1)),Klb(n*(dim+1),nb*dim),Fb(nb*dim)

      INTEGER :: m, i, j, l, p, Cdofs((dim+1)*n), Bdofs(dim*nb)

      m = 0
      DO p = 1,n
        DO i = 1,dim+1
          m = m + 1
          Cdofs(m) = (dim+1)*(p-1) + i
        END DO
      END DO
      
      m = 0
      DO p = 1,nb
        DO i = 1,dim
          m = m + 1
          Bdofs(m) = (dim+1)*(p-1) + i + n*(dim+1)
        END DO
      END DO

      Kbb = K(Bdofs,Bdofs)
      Kbl = K(Bdofs,Cdofs)
      Klb = K(Cdofs,Bdofs)
      Fb  = F(Bdofs)

      CALL InvertMatrix( Kbb,Nb*dim )

      F(1:(dim+1)*n) = F(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
      K(1:(dim+1)*n,1:(dim+1)*n) = &
           K(1:(dim+1)*n,1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )
!------------------------------------------------------------------------------
    END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
END SUBROUTINE FilmFlowSolver
!------------------------------------------------------------------------------
