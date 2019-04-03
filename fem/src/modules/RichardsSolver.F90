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
! *****************************************************************************
!
!******************************************************************************
! *
! *  Richards equation for porous flow
! *
! *****************************************************************************
! *
! *  Authors: Peter Råback, Serge-�tienne Parent
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 9 Apr 2004
! *
! ****************************************************************************


!--------------------------------------------------------------------
!> This module contains porosity models that are computed at the 
!> integration point given the value of matric suction, psy.
!> More materials models could be added.
!--------------------------------------------------------------------
MODULE PorousMaterials

  USE Types
  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE

  INTEGER, PARAMETER :: POROSITY_DEFAULT=0, &
      POROSITY_VAN_GENUCHTEN = 1, &
      POROSITY_BROOKS_COREY = 2 

CONTAINS


  FUNCTION HydraulicConductivity(Element,Material,Basis,elempsy) RESULT ( kw )
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Basis(:),elempsy(:)
    REAL(KIND=dp) :: psy, kw
    REAL(KIND=dp) :: kwsat
    REAL(KIND=dp) :: a, m, n ! van Genuchten
    REAL(KIND=dp) :: aev, lambda ! Brooks and Corey
    REAL(KIND=dp), ALLOCATABLE :: nodalkw(:)
    TYPE(ValueList_t), POINTER :: PrevMaterial => NULL()
    TYPE(Element_t), POINTER :: PrevElement => NULL()
    CHARACTER(LEN=MAX_NAME_LEN) :: PorosityModel
    INTEGER :: nnodes, PorosityModelIndex
    LOGICAL :: Found, SameParameters 
    
    SAVE PrevMaterial,PrevElement,nnodes,PorosityModelIndex,kwsat,a,n,m,aev,lambda,nodalkw 
    
    ! Parameters for the certain model are assumed to be 
    ! constant within a material. This is merely to save a little
    ! bit of time. 
    !--------------------------------------------------------
    SameParameters = .FALSE.
    IF( ASSOCIATED(Material, PrevMaterial) ) THEN
      IF( PorosityModelIndex == POROSITY_VAN_GENUCHTEN .OR. &
          PorosityModelIndex == POROSITY_BROOKS_COREY ) THEN
        SameParameters = .TRUE.
      ELSE IF( ASSOCIATED(Element,PrevElement)) THEN
        SameParameters = .TRUE.
      END IF
    ELSE
      PrevMaterial => Material 
      PorosityModel = GetString( Material,'Porosity Model',Found)
      IF( PorosityModel == 'van genuchten') THEN
        PorosityModelIndex = POROSITY_VAN_GENUCHTEN
      ELSE IF( PorosityModel == 'brooks and corey') THEN 
        PorosityModelIndex = POROSITY_BROOKS_COREY 
      ELSE
        PorosityModelIndex = POROSITY_DEFAULT                
      END IF      
      IF( .NOT. ALLOCATED( nodalkw ) ) THEN
        ALLOCATE( nodalkw(CurrentModel % Mesh % MaxElementNodes))
        nodalkw = 0.0_dp
      END IF
    END IF


    IF( PorosityModelIndex == POROSITY_VAN_GENUCHTEN ) THEN
      IF(.NOT. SameParameters) THEN
        kwsat = GetConstReal( Material,'Saturated Hydraulic Conductivity')
        a = GetConstReal( Material,'van Genuchten Alpha')
        n = GetConstReal( Material,'van Genuchten N')
        m = GetConstReal( Material,'van Genuchten M')
      END IF
      nnodes = Element % TYPE % NumberOfNodes
      psy = SUM( Basis(1:nnodes) * elemPsy(1:nnodes) )
      IF( psy <= 0.0_dp ) THEN
        kw = kwsat
      ELSE
        kw=kwsat*(1-(a*psy)**(n*m)*(1+(a*psy)**n)**(-m))**2*(1+(a*psy)**n)**(-m/2)
      END IF

    ELSE IF( PorosityModelIndex == POROSITY_BROOKS_COREY ) THEN
      IF(.NOT. SameParameters) THEN
        kwsat = GetConstReal( Material,'Saturated Hydraulic Conductivity')
        aev = GetConstReal( Material,'Brooks and Corey Air entry value')
        lambda = GetConstReal( Material,'Brooks and Corey Lambda')
      END IF
      nnodes = Element % TYPE % NumberOfNodes
      psy = SUM( Basis(1:nnodes) * elemPsy(1:nnodes) )
      IF( psy <= aev ) THEN
        kw = kwsat
      ELSE
        kw=kwsat*(psy/aev)**(-2-3*lambda)
      END IF

    ELSE
      nnodes = Element % TYPE % NumberOfNodes
      nodalkw = ListGetReal( Material,'Hydraulic Conductivity',&
          nnodes,Element % NodeIndexes)
      kw = SUM( Basis(1:nnodes) * nodalkw(1:nnodes) )
    END IF
       
  END FUNCTION HydraulicConductivity



  FUNCTION WaterContent( Element, Material, Basis, elemPsy ) RESULT ( teta )
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Basis(:),ElemPsy(:)
    REAL(KIND=dp) :: psy, teta
    REAL(KIND=dp) :: tetaR, tetaS
    REAL(KIND=dp) :: a, m, n ! van Genuchten
    REAL(KIND=dp) :: aev, lambda ! Brooks and Corey
    REAL(KIND=dp), ALLOCATABLE :: Nodalteta(:)
    TYPE(ValueList_t), POINTER :: PrevMaterial => NULL()
    TYPE(Element_t), POINTER :: PrevElement => NULL()
    CHARACTER(LEN=MAX_NAME_LEN) :: PorosityModel
    INTEGER :: PorosityModelIndex, nnodes
    LOGICAL :: SameParameters, Found

    SAVE PrevMaterial,PrevElement,nnodes,PorosityModelIndex,nodalteta,tetaR,tetaS,a,n,m,aev,lambda 
    
    ! Parameters for the certain model are assumed to be 
    ! constant within a material
    !--------------------------------------------------------
    SameParameters = .FALSE.
    IF( ASSOCIATED(Material, PrevMaterial) ) THEN
      IF( PorosityModelIndex == POROSITY_VAN_GENUCHTEN .OR. &
          PorosityModelIndex == POROSITY_BROOKS_COREY ) THEN
        SameParameters = .TRUE.
      ELSE IF( ASSOCIATED(Element,PrevElement)) THEN
        SameParameters = .TRUE.
      END IF
    ELSE
      PrevMaterial => Material 
      PorosityModel = GetString( Material,'Porosity Model',Found)
      IF( PorosityModel == 'van genuchten') THEN
        PorosityModelIndex = POROSITY_VAN_GENUCHTEN
      ELSE IF( PorosityModel == 'brooks and corey') THEN 
        PorosityModelIndex = POROSITY_BROOKS_COREY 
      ELSE
        PorosityModelIndex = POROSITY_DEFAULT                
        IF( .NOT. ALLOCATED( nodalteta ) ) THEN
          ALLOCATE( nodalteta(CurrentModel % Mesh % MaxElementNodes))
          nodalteta = 0.0_dp
        END IF
      END IF
    END IF


    IF( PorosityModelIndex == POROSITY_VAN_GENUCHTEN ) THEN
      IF(.NOT. SameParameters) THEN
        tetaR = GetConstReal( Material,'Residual Water Content')
        tetaS = GetConstReal( Material,'Saturated Water Content')
        a = GetConstReal( Material,'van Genuchten Alpha')
        n = GetConstReal( Material,'van Genuchten N')
        m = GetConstReal( Material,'van Genuchten M')
      END IF

      nnodes = Element % TYPE % NumberOfNodes
      psy = SUM( Basis(1:nnodes) * elemPsy(1:nnodes) )

      IF( psy <= 0.0 ) THEN
        teta = tetaS
      ELSE
        teta=tetaR+(tetaS-tetaR)/((1+(a*psy)**n)**m)
      END IF

    ELSE IF( PorosityModelIndex == POROSITY_BROOKS_COREY ) THEN
      IF(.NOT. SameParameters) THEN
        tetaR = GetConstReal( Material,'Residual Water Content')
        tetaS = GetConstReal( Material,'Saturated Water Content')
        aev = GetConstReal( Material,'Brooks and Corey Air entry value')
        lambda = GetConstReal( Material,'Brooks and Corey Lambda')
      END IF

      nnodes = Element % TYPE % NumberOfNodes
      psy = SUM( Basis(1:nnodes) * elemPsy(1:nnodes) )

      IF( psy <= aev ) THEN
        teta = tetaS
      ELSE
        teta=tetaR+(tetaS-tetaR)*(psy/aev)**(-lambda)
      END IF

    ELSE  
      nnodes = Element % TYPE % NumberOfNodes
      nodalteta = ListGetReal( Material,'Water Content',&
          nnodes,Element % NodeIndexes)
      teta = SUM( Basis(1:nnodes) * nodalteta(1:nnodes) )
    END IF
    
  END FUNCTION WaterContent

END MODULE PorousMaterials



!------------------------------------------------------------------------------
!> Initialization for the primary solver: RichardsSolver
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE RichardsSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

    USE DefUtils
    USE Types

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
    
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: UseDG, Found, Calculate 

    Params => GetSolverParams()
    UseDG = GetLogical( Params,'Discontinuous Galerkin',Found)
    
    IF( UseDG ) THEN
      CALL ListAddString( Params,'Exported Variable 1',&
          'Nodal PressureHead')
    END IF
    
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
            'Matric Suction' )

    ! The sign of total head is different that the default convention of diffusion equations
    IF( .NOT. ListCheckPresent( Params, 'Limiter Load Sign Negative') ) THEN
      CALL ListAddLogical( Params,'Limiter Load Sign Negative',.TRUE.) 	
    END IF
    
  END SUBROUTINE RichardsSolver_Init


!------------------------------------------------------------------------------
!> Solves the Richards equation for porous flow.
!> A variably saturated formulation is used. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE RichardsSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils
    USE PorousMaterials
    USE Types
   
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC, BodyForce, Material, SolverParams    
    TYPE(Element_t), POINTER :: Element
    TYPE( Element_t ), POINTER :: Faces(:)
    TYPE(Nodes_t) :: ElementNodes
    LOGICAL :: AllocationsDone = .FALSE., Found, Stat, Bubbles, UseDG
    INTEGER :: NumberOfFAces
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), TimeForce(:), &
        ElemHead(:),ElemPrevHead(:),ElemMatric(:),ElemPrevMatric(:),&
        ElemSource(:),ElemFlux(:)
    REAL(KIND=dp), POINTER :: Hcoord(:), TotalHead(:), MatricSuction(:)
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: TotalHeadPerm(:)

    LOGICAL :: SubroutineVisited = .FALSE.,InitSolution, ResetRelax
    INTEGER :: i,j,n, nd, t, istat, iter, dim, Active, NonlinearIterMax, &
        ActiveCoordinate
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: at,st,totst,totat
#else
    REAL(KIND=dp) :: at,st,totst,totat,CPUTime
#endif
    REAL(KIND=dp) :: Norm, Relax
    TYPE(Variable_t), POINTER :: Var
    
    
    SAVE MASS, STIFF, FORCE, TimeForce, &
        ElemHead, ElemPrevHead, ElemMatric,ElemPrevMatric,ElemSource, ElemFlux, &
        AllocationsDone,SubroutineVisited
    
     !----------------------------------------------------------------------------
    WRITE(Message,'(A,A)') 'RichardsSolver for variable '// &
        TRIM( Solver % Variable % Name )
    CALL INFO('RichardsSolver',Message,Level=3)
    
    Mesh => GetMesh()
    DIM = CoordinateSystemDimension()
    
    ! Initialize & allocate some permanent storage, this is done first time only:
    !----------------------------------------------------------------------------
    IF ( .NOT. AllocationsDone ) THEN
      N = 2 * MAX(Mesh % MaxElementDOFs, Mesh % MaxElementNodes )
      ALLOCATE( FORCE(N), MASS(n,n), STIFF(N,N), TimeForce(N), &
          ElemHead(N),ElemPrevHead(N),ElemMatric(N),ElemPrevMatric(N), &
          ElemSource(N), ElemFlux(N), &
          STAT = istat )
      
      IF ( istat /= 0 ) CALL FATAL('RichardsSolver','Memory allocation error.' )
      AllocationsDone = .TRUE.
    END IF
    
    !------------------------------------------------------------------------------
    !    Read physical and numerical constants and initialize 
    !------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    
    Bubbles = GetLogical( SolverParams, 'Bubbles', Found )
    
    ActiveCoordinate = GetInteger( SolverParams,'Active Coordinate',Found)
    IF(.NOT. Found) ActiveCoordinate = DIM
    IF( ActiveCoordinate == 1 ) THEN
      Hcoord => Solver % Mesh % Nodes % x
    ELSE IF( ActiveCoordinate == 2 ) THEN
      Hcoord => Solver % Mesh % Nodes % y
    ELSE
      Hcoord => Solver % Mesh % Nodes % z
    END IF

    !------------------------------------------------------------------------------
    ! Field variables 
    !------------------------------------------------------------------------------
    TotalHead => Solver % Variable % Values
    TotalHeadPerm => Solver % Variable % Perm

    Var => VariableGet( Solver % Mesh % Variables,'Matric Suction')
    MatricSuction => Var % Values

    
    UseDG = GetLogical( SolverParams,'Discontinuous Galerkin',Found)
    IF( UseDG ) THEN
      CALL Info('RichardsSolver','Using DG for discretization') 
      IF ( DIM == 2 ) THEN
        Faces => Mesh % Edges
        NumberOfFaces = Mesh % NumberOfEdges
      ELSE
        Faces => Mesh % Faces
        NumberOfFaces = Mesh % NumberOfFaces
      END IF
    END IF
    
    NonlinearIterMax = GetInteger( SolverParams, &
        'Nonlinear System Max Iterations',Found )
    IF ( .NOT.Found ) THEN
      CALL WARN('RichardsSolver','No >Nonlinear System Max Iterations< found. Setting 1')
      NonlinearIterMax = 1
    END IF
        
    !------------------------------------------------------------------------------
    !       non-linear system iteration loop
    !------------------------------------------------------------------------------
    
    totst = 0; totat = 0;
    DO iter=1,NonlinearIterMax
      
      at  = CPUTime()
      
      CALL Info( 'RichardsSolver', ' ', Level=6 )
      CALL Info( 'RichardsSolver', '-------------------------------------',Level=4 )
      WRITE( Message,'(A,I4,A,I4)') &
          'Nonlinear iteration no.', iter,' of',NonlinearIterMax
      CALL Info( 'RichardsSolver', Message, Level=4 )
      CALL Info( 'RichardsSolver', '-------------------------------------',Level=4 )
      CALL Info( 'RichardsSolver', ' ', Level=6 )
      CALL Info( 'RichardsSolver', 'Starting Assembly...', Level=6 )
      
      InitSolution = .FALSE.
      IF( iter == 1 .AND. .NOT. SubroutineVisited ) THEN
        InitSolution = GetLogical( SolverParams,'Saturated Initial Guess',Found )
      END IF
      
      IF( InitSolution ) THEN
        Relax = GetConstReal(SolverParams,&
            'Nonlinear System Relaxation Factor',Found)
        IF(Found) THEN
          CALL ListAddConstReal(SolverParams,&
              'Nonlinear System Relaxation Factor',1.0_dp)
          ResetRelax = .TRUE.
        ELSE
          ResetRelax = .FALSE.  
        END IF
      END IF
      
      CALL DefaultInitialize()
      
      !------------------------------------------------------------------------------
      !  Bulk Assembly
      !------------------------------------------------------------------------------
      Active = GetNOFActive()
      DO t = 1, Active  
        !------------------------------------------------------------------------------
        ! assign pointers and get number of nodes in element
        ! The material parameters are defined in the library but the pointer to 
        ! the correct material must be set here. 
        !------------------------------------------------------------------------------  
        Element => GetActiveElement( t )
        n = GetElementNOfNodes( Element )           
        nd = GetElementNOFDOFs()
        
        Material => GetMaterial()
        
        !------------------------------------------------------------------------------
        ! the body force (r.h.s) = source
        !------------------------------------------------------------------------------         
        BodyForce => GetBodyForce( Element )           
        ElemSource(1:n) = GetReal( BodyForce,'Richards Source', Found )    
        
        CALL GetScalarLocalSolution( ElemHead )
        IF( Transient ) THEN
          CALL GetScalarLocalSolution( ElemPrevHead, tStep = -1 )
        END IF
        ElemMatric(1:n) = Hcoord(Element % NodeIndexes ) - ElemHead(1:n)
        ElemPrevMatric(1:n) = Hcoord(Element % NodeIndexes) - ElemPrevHead(1:n)

        CALL LocalBulkMatrix( MASS, STIFF, FORCE, &
            ElemMatric, ElemPrevMatric, ElemSource, Element, n, nd ) 
        
        TimeForce  = 0.0_dp
        IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
        
        IF (  Bubbles ) THEN
          CALL Condensate( N, STIFF, FORCE, TimeForce )
        END IF
        
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO

      CALL DefaultFinishBulkAssembly( )

      
      !------------------------------------------------------------------------------
      !  Boundary Assembly
      !------------------------------------------------------------------------------
      DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE
       
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        
        BC => GetBC()
        IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
        
        ElemFlux(1:n) = GetReal(BC,'Richards Flux',Found)
        IF(.NOT. Found) CYCLE
        
        CALL LocalBoundaryMatrix( MASS, STIFF, FORCE, &
            ElemFlux, Element, n )
        
        IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE )

        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO
      
      CALL DefaultFinishAssembly()
      CALL Info( 'RichardsSolver', 'Assembly done', Level=6 )
      
      CALL DefaultDirichletBCs()
      CALL Info( 'RichardsSolver', 'Dirichlet conditions done', Level=6 )
      
      !------------------------------------------------------------------------------
      !     Solve the system and check for convergence
      !------------------------------------------------------------------------------
      at = CPUTime() - at
      st = CPUTime()
      
      ! Solve the system:
      !------------------
      Norm = DefaultSolve()
      
      IF(InitSolution .AND. ResetRelax ) THEN
        CALL ListAddConstReal(SolverParams,&
            'Nonlinear System Relaxation Factor',Relax)
      END IF
      
      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info( 'RichardsSolver', Message, Level=4 )
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info( 'RichardsSolver', Message, Level=4 )
      
      IF ( Solver % Variable % NonlinConverged == 1 ) THEN
        WRITE(Message,'(A,I6,A,I6,A)') &
            'Nonlinear iteration converged after ', iter, &
            ' out of max ',NonlinearIterMax,' iterations'
        CALL INFO('RichardsSolver',Message, Level=4)
        EXIT
      END IF
      
      ! Update matric suction which may be used in the nonlinear iteration as 
      !----------------------------------------------------------------------
      DO i = 1, Solver % Mesh % NumberOfNodes
        j = TotalHeadPerm( i )
        IF( j == 0 ) CYCLE
        MatricSuction(j) = Hcoord(i) - TotalHead(j)
      END DO

      SubroutineVisited = .TRUE.
      
    END DO ! End of nonlinear iteration loop
    !----------------------------------------------
    
    IF( Solver % Variable % NonlinConverged == 0 ) THEN
      CALL WARN('RichardsSolver','Maximum nonlinear iterations reached, but system not converged')       
    END IF




  CONTAINS


!------------------------------------------------------------------------------------
!> The water content derivative is computed using the real differential. This way the 
!> differential of water content with time will be consistent.
!------------------------------------------------------------------------------------
    FUNCTION WaterContentDerivative( Element, Material, Basis, elemMatric, elemPrevMatric ) RESULT ( dtetadmatric )
      
      TYPE(Element_t), POINTER :: Element
      TYPE(ValueList_t), POINTER :: Material
      REAL(KIND=dp) :: Basis(:),ElemMatric(:),ElemPrevMatric(:)
      REAL(KIND=dp) :: matric, prevmatric, dtetadmatric, &
          teta, prevteta, eps = 1.0d-6
      
      teta = WaterContent( Element, Material, Basis, elemmatric ) 

      n = Element % TYPE % NumberOfNodes 
      matric = SUM( Basis(1:n) * elemMatric(1:n) ) 
      prevmatric = SUM( Basis(1:n) * elemPrevMatric(1:n) )

      IF( ABS( matric - prevmatric ) > eps ) THEN
        prevteta = WaterContent( Element, Material, Basis, elemPrevMatric ) 
        dtetadmatric = (teta - prevteta) / (matric - prevmatric)
      ELSE
        prevteta = WaterContent( Element, Material, Basis, elemMatric - eps ) 
        dtetadmatric = (teta - prevteta) / eps 
      END IF
      
    END FUNCTION WaterContentDerivative
    
    

!------------------------------------------------------------------------------      
     SUBROUTINE LocalBulkMatrix(MASS, STIFF, FORCE, &
         ElemMatric, ElemPrevMatric, ElemSource, Element, n, nd)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), &
           ElemMatric(:), ElemPrevMatric(:), ElemSource(:)
       INTEGER :: n, nd
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3)
       REAL(KIND=dp) :: x,y,z,detJ,SqrtMetric
       REAL(KIND=dp) :: U,V,W,S,ipMatric,ipPrevMatric,ipCond,&
           ipSource,ipContent
       LOGICAL :: Stat
       INTEGER :: i,j,p,q,t,dim,Nbasis, CoordSys
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       SAVE Nodes
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       CoordSys = CurrentCoordinateSystem()
!------------------------------------------------------------------------------
       FORCE = 0.0d0
       STIFF = 0.0d0
       MASS  = 0.0d0
       CALL GetElementNodes( Nodes, Element )

!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       NBasis = n
       IF ( Bubbles ) THEN
         NBasis = 2 * n
         IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints2 )
       ELSE
         NBasis = nd
         IntegStuff = GaussPoints( Element )
       END IF
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
             Basis, dBasisdx, Bubbles = Bubbles )
         
         S = S * detJ

         IF ( CoordSys /= Cartesian ) THEN
           X = SUM( Nodes % X(1:n) * Basis(1:n) )
           S = S * X
         END IF

         ipCond = HydraulicConductivity(Element,Material,Basis,elemMatric)

         IF( Transient .AND. .NOT. InitSolution ) THEN
           ipContent = WaterContentDerivative(Element,Material,Basis,elemMatric,elemPrevMatric)
         ELSE
           ipContent = 0.0_dp
         END IF
 
         ipSource = SUM( ElemSource(1:n) *  Basis(1:n) )
         
!------------------------------------------------------------------------------
!        The Richards equation
!------------------------------------------------------------------------------
         DO p=1,Nbasis
           DO q=1,Nbasis
             MASS(p,q)  = MASS(p,q) + s * ipContent * Basis(q) * Basis(p)
             STIFF(p,q) = STIFF(p,q) - s * ipCond * &
                 SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )
           END DO           
           FORCE(p) = FORCE(p) + s * ipSource * Basis(p)
         END DO
         
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalBulkMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
     SUBROUTINE LocalBoundaryMatrix(MASS, STIFF, FORCE, &
         ElemFlux, Element, n )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), ElemFlux(:)
       INTEGER :: n, nd
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: x,y,z,detJ,SqrtMetric
       REAL(KIND=dp) :: U,V,W,S,ipFlux
       LOGICAL :: Stat
       INTEGER :: i,j,p,q,t,dim,CoordSys
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       SAVE Nodes
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       CoordSys = CurrentCoordinateSystem()
!------------------------------------------------------------------------------
       FORCE = 0.0d0
       STIFF = 0.0d0
       MASS  = 0.0d0
       CALL GetElementNodes( Nodes, Element )
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------

       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
             Basis, dBasisdx )
         
         S = S * detJ
         
         IF ( CoordSys /= Cartesian ) THEN
           X = SUM( Nodes % X(1:n) * Basis(1:n) )
           s = s * X
         END IF
         
         ipFlux = SUM( Basis(1:n) * ElemFlux(1:n) )
         
         DO p=1,n
           FORCE(p) = FORCE(p) + s * ipFlux * Basis(p)
         END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalBoundaryMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   END SUBROUTINE RichardsSolver
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
!> Initialization for the primary solver: RichardsPostprocess
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE RichardsPostprocess_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
    INTEGER :: dim
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 1 )
      CALL ListAddString( SolverParams, 'Variable','-nooutput flux_temp' )
      IF(dim == 2) THEN
        CALL ListAddString(SolverParams, 'Exported Variable 1','Richards Flux[Richards Flux:2]')
      ELSE IF(dim == 3) THEN
        CALL ListAddString(SolverParams, 'Exported Variable 1','Richards Flux[Richards Flux:3]')
      END IF
    END IF
    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )

    ! Add linear system defaults: cg+diagonal
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
      CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
      CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
      CALL ListAddString(SolverParams,'Linear System Preconditioning','diagonal')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
      CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
      CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
      CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-10_dp)

!------------------------------------------------------------------------------
  END SUBROUTINE RichardsPostprocess_Init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Solve the flux resulting from the solution of the Richards equation using
!>  the Galerkin method. The solver is somewhat complicated because the components
!>  of the flux (f_x,f_y,f_z) are computed one-by-one.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE RichardsPostprocess( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE CoordinateSystems
  USE DefUtils
  USE PorousMaterials

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i,j,dim,DOFs,ActiveCoordinate
  LOGICAL :: Found, ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry
  REAL(KIND=dp) :: Unorm, Totnorm, FluxMultiplier
  REAL(KIND=dp), POINTER CONTIG :: ForceVector(:,:), SaveRHS(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: FluxSol
  
 
  CALL Info( 'RichardsPostprocess', '-------------------------------------',Level=4 )
  CALL Info( 'RichardsPostprocess','Computing the flux',Level=4 )
  CALL Info( 'RichardsPostprocess', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  
  SolverParams => GetSolverParams()

  ActiveCoordinate = GetInteger( SolverParams,'Active Coordinate',Found)
  IF(.NOT. Found) ActiveCoordinate = dim

!-------------------------------------------------------------------------------
! Use an auxiliary variable to store all the dimensions
!-------------------------------------------------------------------------------
  VarName = GetString(SolverParams,'Flux Result Variable',Found )
  IF(.NOT. Found) VarName = 'Richards Flux'
  FluxSol => VariableGet( Solver % Mesh % Variables,VarName)
  IF(ASSOCIATED(FluxSol)) THEN
    Dofs = FluxSol % DOFs
    IF(Dofs /= DIM) THEN
      CALL Fatal('RichardsPostprocess','The flux should have DOFs equal to DIM')
    END IF
  ELSE
    CALL Fatal('RichardsPostprocess','Flux Result Variable is missing: '//TRIM(VarName))      
  END IF
  
  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  
  VarName = GetString(SolverParams,'Target Variable',Found )
  IF(.NOT. Found) VarName = 'Total Head'

  at0 = RealTime()
  
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', Found )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
    Solver % Matrix % rhs = 0.0_dp
  ELSE
    CALL DefaultInitialize()
  END IF

  ! We need DIM r.h.s. vectors, allocated DIM-1 additional ones  
  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),Dofs-1))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS
  
  CALL BulkAssembly()
  CALL DefaultFinishBulkAssembly()

  CALL DefaultFinishAssembly()
  
  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'RichardsPostprocess', Message, Level=5 )

!------------------------------------------------------------------------------     

  TotNorm = 0.0_dp
  DO i=1,Dofs
    IF(i==1) THEN
      Solver % Matrix % RHS => SaveRHS
    ELSE
      Solver % Matrix % RHS => ForceVector(:,i-1)
    END IF
    UNorm = DefaultSolve()
    TotNorm = TotNorm + Unorm ** 2
    DO j=1,Solver % Matrix % NumberOfRows
      FluxSol % Values(DOFs*(j-1)+i) = Solver % Variable % Values(j)
    END DO
  END DO

! This may be used to multiply the resulting flux to a more convenient value range
  FluxMultiplier = GetConstReal( SolverParams,'Flux Multiplier',Found) 
  IF( Found ) THEN
    FluxSol % Values = FluxMultiplier * FluxSol % Values
  END IF


  DEALLOCATE( ForceVector )  
  Solver % Matrix % RHS => SaveRHS
  TotNorm = SQRT(TotNorm)
  Solver % Variable % Norm = Totnorm

!------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'RichardsPostprocess', Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( 'RichardsPostprocess', Message, Level=4 )
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,grad(3),ipCond,ipMatric,ipZ,coeff,detJ,z
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: ElemMatric(:),ElemHead(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    
    SAVE Nodes
    
    n = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(dim,n) )
    ALLOCATE( ElemMatric(n), ElemHead(n), Basis(n), dBasisdx(n,3) )

    DO elem = 1,Solver % NumberOFActiveElements
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      NodeIndexes => Element % NodeIndexes
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      Material => GetMaterial()
      
      CALL GetScalarLocalSolution( elemMatric,'matric suction')
      CALL GetScalarLocalSolution( elemHead, VarName )
      

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )
        
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        ipCond = HydraulicConductivity(Element,Material,Basis,elemMatric)
        Grad(1:dim) = MATMUL( elemHead(1:nd), dBasisdx(1:nd,1:dim) )
        
        DO i=1,dim
          FORCE(i,1:nd) = FORCE(i,1:nd) - &
              Weight * ipCond * Grad(i) * Basis(1:nd)
        END DO
      END DO

!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------
      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        Solver % Matrix % RHS => SaveRHS
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      ELSE
        CALL DefaultUpdateForce( FORCE(1,1:nd) )       
      END IF

      DO i=2,Dofs
        Solver % Matrix % RHS => ForceVector(:,i-1)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO
    END DO
    
    DEALLOCATE( elemMatric, STIFF, FORCE, Basis, dBasisdx )
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE RichardsPostprocess
!------------------------------------------------------------------------------


