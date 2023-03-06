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
  
  IF(Serendipity) THEN
    CALL ListAddNewString(Params,'Element', &
        'p:1 -line b:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4')
  ELSE
    CALL ListAddNewString(Params,'Element', &
        'p:1 -line b:1 -tri b:1 -tetra b:1 -quad b:4 -brick b:8 -prism b:4 -pyramid b:4')
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
  INTEGER :: i,j,k,n, nb, nd, t, istat, dim, mdim, BDOFs=1,Active,iter,maxiter
  REAL(KIND=dp) :: Norm = 0, PrevNorm, RelC, mingap
  TYPE(ValueList_t), POINTER :: Params, BodyForce, Material, BC
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), &
      FORCE(:), rho(:), gap(:), mu(:), ac(:), NormalVelo(:), Pres(:), Velocity(:,:), MASS(:,:),&
      PseudoPressure(:)  
  LOGICAL :: GradP, LateralStrain, GotAc, SurfAc
  TYPE(Variable_t), POINTER :: pVar
  CHARACTER(*), PARAMETER :: Caller = 'FilmFlowSolver'
 
  SAVE STIFF, MASS, LOAD, FORCE, rho, ac, gap, mu, NormalVelo, Pres, Velocity, &
      PseudoPressure, AllocationsDone, pVar, GotAc, SurfAC
!------------------------------------------------------------------------------

  CALL Info(Caller,'Computing reduced dimensional Navier-Stokes equations!')

  Mesh => GetMesh()
  Element => GetActiveElement(1)

  dim = CoordinateSystemDimension()
    
  Params => GetSolverParams()

  mdim = GetInteger( Params,'Model Dimension')
  Convect = GetLogical( Params, 'Convect', Found )
  GradP = GetLogical( Params, 'GradP Discretization', Found ) 
  LateralStrain = GetLogical( Params,'Lateral Strain',Found )
  mingap = ListGetCReal( Params,'Min Gap Height',Found )
  IF(.NOT. Found) mingap = TINY(mingap)
  GotAC = ListCheckPresentAnyMaterial( Model,'Artificial Compressibility')

  CSymmetry = ListGetLogical( Params,'Axi Symmetric',Found )
  IF(.NOT. Found ) THEN 
    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) 
  END IF
    
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    CALL Info(Caller,'Dimension of Navier-Stokes equation: '//I2S(mdim))
    CALL Info(Caller,'Dimension of coordinate system: '//I2S(dim))

    n = (mdim+1)*(Mesh % MaxElementDOFs+BDOFs)  ! just big enough for elemental arrays
    ALLOCATE( FORCE(n), LOAD(n,4), STIFF(n,n), MASS(n,n), &
        rho(n), ac(n), gap(n), mu(n), NormalVelo(n), Pres(n), Velocity(mdim+1,n), STAT=istat )
    Velocity = 0.0_dp
    IF ( istat /= 0 ) THEN
      CALL Fatal( Caller, 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
    IF( GradP ) THEN
      CALL Info(Caller,'"Gradp Discretization" is set True',Level=10)
    END IF     
    IF( GotAC ) THEN
      pVar => VariableGet( Mesh % Variables,'FilmPressure')
      IF ( .NOT. ASSOCIATED(pVar) ) THEN
        CALL Fatal( Caller, 'Could not find required field "FilmPressure"!')
      END IF
      ALLOCATE(PseudoPressure(SIZE(pVar % Values)))
      PseudoPressure = 0.0_dp

      SurfAc = ListGetLogical( Params,'Surface Compressibility',Found )
    END IF
  END IF

  IF(GotAc) THEN  
    PseudoPressure = pVar % Values
    WRITE(Message,'(A,T25,E15.4)') 'PseudoPressure mean: ',&
        SUM(PseudoPressure)/SIZE(PseudoPressure)
    CALL Info(Caller,Message,Level=5)
  END IF
   
  maxiter = ListGetInteger( Params, &
      'Nonlinear System Max Iterations',Found,minv=1)
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

      !IF(t==1) PRINT *,'Element:',Element % TYPE % ElementCode, n, nd, nb
      
      ! Volume forces:
      !---------------
      BodyForce => GetBodyForce()
      LOAD = 0.0d0
      IF ( ASSOCIATED(BodyForce) ) THEN
        Load(1:n,1) = GetReal( BodyForce, 'Flow Bodyforce 1', Found )
        Load(1:n,2) = GetReal( BodyForce, 'Flow Bodyforce 2', Found )
        Load(1:n,3) = GetReal( BodyForce, 'Flow Bodyforce 3', Found )
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

      IF( GotAC ) THEN
        ac(1:n) = GetReal( Material,'Artificial Compressibility',Found )
        Pres(1:n) = PseudoPressure(pVar % Perm(Element % NodeIndexes))
      END IF

      NormalVelo(1:n) = GetReal( Material,'Normal Velocity',Found )
            
      ! Get previous elementwise velocity iterate:
      ! Note: pressure is the dim+1 component here!
      !-------------------------------------------
      IF ( Convect ) CALL GetVectorLocalSolution( Velocity )
      
      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(  MASS, STIFF, FORCE, LOAD, rho, gap, mu, &
          ac, NormalVelo, Velocity, Pres, Element, n, nd, nd+nb, &
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

  CALL Info(Caller,'All done',Level=12)

  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(  MASS, STIFF, FORCE, LOAD, Nodalrho, NodalGap, &
      Nodalmu, NodalAC, NodalNormalVelo, NodalVelo, NodalPres, Element, n, nd, &
      ntot, dim, mdim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), NodalNormalVelo(:), NodalAC(:), Nodalrho(:), &
        NodalGap(:), NodalPres(:), NodalVelo(:,:)
    INTEGER :: dim, mdim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3), &
        DetJ,LoadAtIP(dim+1),Velo(dim), VeloGrad(dim,dim), gapGrad(dim),meangap
    REAL(KIND=dp), POINTER :: A(:,:), F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q, geomc
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, normalvelo, pres, gap, ac, s, s0, s1, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0

    gapGrad = 0.0_dp
    meangap = SUM( NodalGap(1:n) ) / n
       
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
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

       !rho = rho * gap
       
       DO i=1,dim
         gapGrad(i) = SUM( NodalGap(1:nd) * dBasisdx(1:nd,i) )
       END DO
       
       normalVelo = SUM( Basis(1:n) * NodalNormalVelo(1:n) )
       
       ! Previous velocity at the integration point:
       !--------------------------------------------
       Velo = MATMUL( NodalVelo(1:mdim,1:nd), Basis(1:nd) )
       VeloGrad = MATMUL( NodalVelo(1:mdim,1:nd), dBasisdx(1:nd,1:mdim) )

       IF( GotAC ) THEN
         Pres = SUM( NodalPres(1:n) * Basis(1:n) )
         ac = SUM( NodalAC(1:n) * Basis(1:n) ) / dt
         IF(SurfAc) ac = ac / gap
       END IF
         
       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP(1:mdim+1) = MATMUL( Basis(1:n), LOAD(1:n,1:mdim+1) )
       IF ( Convect .AND. Newton ) THEN
         LoadAtIp(1:mdim) = LoadAtIp(1:mdim) + rho * MATMUL(VeloGrad,Velo)
       END IF

       ! To my understanding we want to include the gap height to weight
       IF( Csymmetry ) THEN
!         s = s * (meangap**2)
         geomc = 2
       ELSE
!         s = s * meangap
         geomc = 1
       END IF
       
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
             IF( CSymmetry ) THEN
               A(i,i) = A(i,i) + s * ( 8 * mu / gap**2 ) * Basis(q) * Basis(p)  
             ELSE
               A(i,i) = A(i,i) + s * ( 12 * mu / gap**2 ) * Basis(q) * Basis(p)  
             END IF
               
             ! Note that here the gap height must be included in the continuity equation
             IF( GradP ) THEN
               A(i,mdim+1) = A(i,mdim+1) + s * dBasisdx(q,i) * Basis(p)
               A(mdim+1,i) = A(mdim+1,i) - s * rho * Basis(q) * dBasisdx(p,i)               
             ELSE
               A(i,mdim+1) = A(i,mdim+1) - s * Basis(q) * dBasisdx(p,i)
               A(mdim+1,i) = A(mdim+1,i) + s * gap * rho * dBasisdx(q,i) * Basis(p) & 
                   + geomc * s * rho * Basis(q) * gapGrad(i) * Basis(p)
!                   + 1.0*(geomc*s/gap) * rho * Basis(q) * gapGrad(i) * Basis(p)

               !A(i,dim+1) = A(i,dim+1) - s * Basis(q) * dBasisdx(p,i)
               !A(dim+1,i) = A(dim+1,i) + s * rho * dBasisdx(q,i) * Basis(p)
               

             END IF

             ! This is the flow into the domain from side
             A(mdim+1,i) = A(mdim+1,i) + geomc * s * NormalVelo * rho * Basis(q) * Basis(p)
               
             ! This is the artificial compressibility for FSI coupling
             IF(GotAC) A(mdim+1,mdim+1) = ac * s * rho * Basis(q) * Basis(p)              
           END DO
         END DO
         i = (mdim+1) * (p-1) + 1
         F => FORCE(i:i+mdim)

         F = F + s * LoadAtIP * Basis(p) 
         IF( GotAC ) F = F + ac * s * rho * Basis(p) * Pres
       END DO
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
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Load(dim+1)
    REAL(KIND=dp), POINTER :: A(:,:),F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q, geomc
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, normalvelo, pres, gap, ac, s, c
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0

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
         
       ! To my understanding we want to include the gap height to weight
       IF( Csymmetry ) THEN
         geomc = 2
       ELSE
         geomc = 1
       END IF
       
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
