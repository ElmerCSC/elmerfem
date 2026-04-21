!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Module for solving 1dof Ordinary Differential Equation (ODE) solver on moving
! *  coordinates. The mesh is assumed to be extruded and we assume constant draw
! *  velocity. Hence we can march the ODE at material point in time. 
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 4.12.2019
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{
 

!------------------------------------------------------------------------------
!> Initialization for the primary solver
!------------------------------------------------------------------------------
SUBROUTINE MarchingODESolver_init( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  
  TYPE(Model_t) :: Model    
  REAL(KIND=dp) :: dt       
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'MarchingODESolver_init'
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams()

  ! By construction time derivative order must be one!
  CALL ListAddNewInteger( Params, 'Time derivative order', 1 )

  ! Global mass matrix allows us to use library routines for time integration
  CALL ListAddNewLogical( Params, 'Use Global Mass Matrix', .TRUE. )

  ! We initialize our own matrix structures
  CALL ListAddNewLogical( Params, 'No Matrix', .TRUE. )
    
!------------------------------------------------------------------------------
END SUBROUTINE MarchingODESolver_init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This solver combines a simple solver for ODE in time with marching along
!> structured mesh.
!------------------------------------------------------------------------------
SUBROUTINE MarchingODESolver( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MeshUtils, ONLY : DetectExtrudedStructure, MarkBCNodes
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
  LOGICAL :: Transient              !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'MarchingODESolver'
  LOGICAL :: Found
  REAL(KIND=dp) :: Norm, Change, dz, dtime, velo, NonLinTol, Beta, &
      Hparam, dth, time, mincons, sumcons
  INTEGER :: t,i,j,n,m,iter,MaxIter,TimeOrder,BotNodes,layer,dtn,dti,NoActive
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Element_t), POINTER :: Element
  INTEGER, POINTER :: BotPointer(:), UpPointer(:)
  INTEGER, POINTER :: BotPerm(:),InvPerm(:),PrevInvPerm(:),MaskPerm(:),SingleIndex(:),Node2DG(:)
  INTEGER, ALLOCATABLE :: ParentElem(:),DGIndexes(:)
  INTEGER :: NumberOfLayers, NoBCNodes, dofs, subt, maxsubt
  TYPE(Variable_t), POINTER :: ExtVar, Var3D, AddVar
  TYPE(ValueList_t), POINTER :: Material
  LOGICAL :: MaskExist, ParabolicModel, RequireBC, DoTransient, AnyDG, VectorSource
  REAL(KIND=dp), POINTER :: Coord(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TimeMethod, VarName, str
  LOGICAL, ALLOCATABLE :: BCNode(:)
  REAL(KIND=dp), POINTER :: xvec(:),xivec(:),dxvec(:),x0vec(:),&
      fvec(:),rvec(:),cvec(:),f0vec(:),r0vec(:),c0vec(:)
  LOGICAL :: HaveF, HaveC, HaveR, UseInternalVals, SetMin, SetSum  
  LOGICAL, SAVE :: Initialized = .FALSE.
!------------------------------------------------------------------------------

  SAVE :: BotPointer, UpPointer, BotPerm, InvPerm, PrevInvPerm, ParentElem, &
      SingleIndex, MaskPerm, MaskExist, NumberOfLayers, ExtVar, BotNodes, &
      TimeMethod, RequireBC, UseInternalVals, AddVar, AnyDG, HParam, SetMin, SetSum, mincons, &
      xivec, dxvec, xvec, fvec, rvec, cvec, &
      x0vec, f0vec, r0vec, c0vec, Coord
  
  CALL Info(Caller,'-----------------------------------------------------',Level=6)
  CALL Info(Caller,'Solving ODE on moving coordinates in structured mesh',Level=4)

  Params => GetSolverParams()
  Mesh => Solver % Mesh
  PSolver => Solver
  
  ! In principle we could use the same solver to advect many fields on
  ! same extruded mesh. Hence these are not saved. 
  Var3D => Solver % Variable
  IF(.NOT. ASSOCIATED(Var3D)) THEN
    CALL Fatal(Caller,'Solver variable not associated!')
  END IF
    
  VarName = TRIM( Var3D % Name )  
  CALL Info(Caller,'Working with variable: '//TRIM(VarName),Level=7)

  dofs = Var3D % dofs
  CALL Info(Caller,'Number of components in variable: '//I2S(dofs),Level=7)
  IF(dofs<1) CALL Fatal(Caller,'Invalid number of components!')
    
  IF( .NOT. Initialized ) THEN
    ! Just check what type of elements we have 
    NoActive = GetNOFActive()
    DO t=1,NoActive
      Element => GetActiveElement(t)
      i = Element % TYPE % ElementCode
      IF(t == 1) THEN
        j = i
      ELSE IF(i /= j) THEN
        EXIT
      END IF
    END DO
    IF( i /= j ) THEN
      CALL Info(Caller,'There are at least elements of type: '//I2S(i)//' and '//I2S(j))
    ELSE
      CALL Info(Caller,'All elements are of type: '//I2S(i),Level=12)
    END IF

  
    CALL Info(Caller,'Initializing structured mesh and ODE structures',Level=6)
    CALL Info(Caller,'Solving for variable: '//TRIM(VarName),Level=6)

    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = ExtVar, &
        BotNodePointer = BotPointer, UpNodePointer = UpPointer, &
        NumberOfLayers = NumberOfLayers, MaskVar = Var3D )
        
    CALL Info(Caller,'Number of element layers: '//I2S(NumberOfLayers),Level=7)

    MaskExist = ASSOCIATED( ExtVar % Perm ) 
    IF( MaskExist ) THEN
      CALL Info(Caller,'We have a mask of size:'//I2S(MAXVAL(ExtVar % Perm)),Level=7)
    ELSE
      CALL Info(Caller,'No mask associated to solver',Level=20)
    END IF

    IF( MaskExist ) MaskPerm => ExtVar % Perm
    Coord => ExtVar % Values

    IF(.NOT. ASSOCIATED(Coord) ) THEN
      CALL Fatal(Caller,'Coord is not associated!')
    END IF
    
    WRITE(Message,'(A,2ES12.3)') 'Extruded Coordinate Range:',MINVAL(Coord),MAXVAL(Coord)
    CALL Info(Caller,Message,Level=8)
       
    dz = MAXVAL(Coord)-MINVAL(Coord)
    IF( dz < EPSILON(dz) ) THEN
      CALL Fatal(Caller,'We cannot march when dz is ~zero!')
    END IF
    
    AnyDG = ListGetLogicalAnySolver( Model,'Discontinuous Galerkin')

    UseInternalVals = ListGetLogical( Params,'Use Internal Values', Found ) 
    str = ListGetString( Params,'Additional Internal Variable', Found ) 
    IF(Found) THEN
      AddVar => VariableGet( Mesh % Variables, str )
      CALL Info(Caller,'Using additional internal variable: '//TRIM(AddVar % Name))      
    END IF

    mincons = ListGetCReal( Params,'Minimum Cons',SetMin )
    SetSum = ListGetLogical( Params,'Enforce Unity Sum',Found )    
    
    
    ! It is not trivial to know to which element a node belongs to.
    ! This structure is needed when we want to know the DG field value of a given node.
    ! Only if we also have "Discontinuous Bodies" within this active domain will this be
    ! uniquely defined. 
    !-----------------------------------------------------------------------------------
    IF( MaskExist .OR. AnyDG ) THEN
      CALL Info(Caller,'Creating inverse node parent look-up table',Level=7)
      ALLOCATE( ParentElem(Mesh % NumberOfNodes) )
      ParentElem = 0

      NoActive = GetNOFActive()
      DO t=1,NoActive
        Element => GetActiveElement(t)
        n = Element % TYPE % NumberOfNodes
        IF( MaskExist ) THEN
          IF( ANY( MaskPerm(Element % NodeIndexes) == 0 ) ) CYCLE
        END IF
        DO i=1,n
          j = Element % NodeIndexes(i)
          IF( ParentElem(j) == 0 ) ParentElem(j) = Element % ElementIndex
        END DO
      END DO
    END IF
    
    Hparam = ( MAXVAL( Coord ) - MINVAL( Coord ) ) / NumberOfLayers
    WRITE(Message,'(A,ES12.3)') 'Constant mesh parameter: ',Hparam
    CALL Info(Caller,Message,Level=7)
    
    ! We may choose only to apply the ODE to BC nodes
    RequireBC = ListGetLogical( Params,'Apply BCs Only',Found )
    IF( RequireBC ) THEN
      CALL MarkBCNodes( Mesh,BCNode,NoBCNodes)
      IF(NoBCNodes == 0 ) RequireBC = .FALSE.
      CALL Info(Caller,'Number of BC nodes: '//I2S(NoBCNodes),Level=7)
    END IF
    
    ! Create the permutation using the bottom layer
    BotNodes = 0
    ALLOCATE( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    DO i=1,Mesh % NumberOfNodes
      ! The variable to be marched does not exist at the 1st layer
      IF( Var3D % Perm(i) == 0 ) CYCLE
      
      IF( RequireBC ) THEN
        IF( .NOT. BcNode(i) ) CYCLE
      END IF
      
      j = i
      IF( MaskExist ) THEN
        j = MaskPerm(i)
        IF( j == 0 ) CYCLE
      END IF
        
      ! This is not at the bottom
      IF(BotPointer(j) /= i) CYCLE

      ! Ok, check that also the next node would be at BC
      ! The 1st layer may be an exception, check for 2nd too

      ! The variable to be marched does not exist 
      IF( Var3D % Perm(UpPointer(j)) == 0 ) CYCLE

      ! We want BC node but don't have one
      IF( RequireBC ) THEN
        IF( .NOT. BcNode(UpPointer(j)) ) CYCLE
      END IF      
      
      BotNodes = BotNodes + 1
      BotPerm(i) = BotNodes
    END DO
    n = BotNodes
    ALLOCATE( InvPerm(n), PrevInvPerm(n), SingleIndex(1) )

    CALL Info(Caller,'Number of bottom nodes: '//I2S(n),Level=7)
    
    ! Allocate some vectors to study convergence 
    m = dofs * n
    ALLOCATE( xvec(m), fvec(m), rvec(m), cvec(m), f0vec(m), r0vec(m), &
        c0vec(m), xivec(m), dxvec(m), x0vec(m) )
       
    Initialized = .TRUE.
    CALL Info(Caller,'Initialization done',Level=10)
  END IF
  
  ! The variable on the layer
  n = BotNodes
  m = dofs * n
  
  ! We just use one parameter to define the timestepping.
  ! This defines how the coefficients are to be evaluated. 
  Beta = ListGetCReal( Params,'Newmark Beta',Found )
  IF(.NOT. Found ) THEN
    ! Default timestepping is impicit euler
    Beta = 1.0_dp
    TimeMethod = ListGetString( Params, 'Timestepping Method',Found )
    IF( Found ) THEN      
      IF( TimeMethod == 'implicit euler' ) THEN
        Beta = 1.0_dp
      ELSE IF( TimeMethod == 'explicit euler' ) THEN
        Beta = 0.0_dp
      ELSE IF( TimeMethod == 'crank-nicolsen' ) THEN
        Beta = 0.5_dp
      ELSE IF( TimeMethod == 'newmark' ) THEN
        Beta = ListGetCReal( Params,'Newmark Beta',UnfoundFatal=.TRUE. )
      END IF
    END IF
  END IF
  
  ParabolicModel = ListGetLogical( Params,'Parabolic Model',Found )
  IF( ParabolicModel ) THEN
    CALL Info(Caller,'Using parabolic growth model',Level=7)
  END IF
  
  velo = GetCReal( Model % Simulation,'Draw Velocity',Found )
  IF(.NOT. Found ) velo = GetCReal( Params,'Draw Velocity',Found )
  IF(.NOT. Found ) THEN
    CALL Fatal(Caller,'>Draw Velocity< is needed for marching solver!')
  END IF
    
  MaxIter = GetInteger( Params,'Nonlinear System Max Iterations',Found )
  IF(.NOT. Found) MaxIter = 1
  CALL Info(Caller,'Number of nonlinear iterations set to: '//I2S(MaxIter),Level=7)
  
  NonLinTol = GetCReal( Params,'Nonlinear System Convergence Tolerance',Found )

  Material => FirstExtrudedMaterial()
  
  xvec = 0.0_dp; fvec = 0.0_dp; rvec = 0.0_dp; cvec = 0.0_dp
  x0vec = 0.0_dp; f0vec = 0.0_dp; r0vec = 0.0_dp; c0vec = 0.0_dp
  xivec = 0.0_dp; dxvec = 0.0_dp

  HaveF = ListCheckPresent( Material,TRIM(VarName)//': Source')
  HaveR = ListCheckPresent( Material,TRIM(VarName)//': Reaction Coefficient')
  HaveC = ListCheckPresent( Material,TRIM(VarName)//': Time Derivative Coefficient')

  VectorSource = .FALSE.
  IF(dofs > 1 ) THEN
    IF(HaveF .OR. HaveR .OR. HaveC ) THEN
      VectorSource = .TRUE.
    ELSE            
      HaveF = ListCheckSuffix( Material,': Source')
      HaveR = ListCheckSuffix( Material,': Reaction Coefficient')
      HaveC = ListCheckSuffix( Material,': Time Derivative Coefficient')
    END IF
  END IF
    
  IF( HaveR ) THEN
    CALL Fatal(Caller,'Code some more to account for reaction term!')
  ELSE IF( .NOT. HaveC ) THEN
    CALL Info(Caller,'By default "Time Derivative Coefficient" one will be used',Level=7)
  END IF

  dtn = 0  
  IF( Transient ) THEN
    ! timestep defined by velocity & mesh parameter
    dth = Hparam / velo 
    dtn = NINT(dt / dth ) 
    IF( dtn >= NumberOfLayers ) THEN
      CALL Info(Caller,'Timestep so large than we can use steady algo!')
      dtn = 0
    ELSE    
      IF( ABS( dt/dth - dtn ) > 0.01 ) THEN
        PRINT *,'Mesh parameter:',Hparam
        PRINT *,'Draw Velocity:',velo
        PRINT *,'Suggested timesteps:',dt,dth
        CALL Fatal(Caller,'Timesteps are not matching')        
      ELSE
        CALL Info(Caller,'Number of marching steps for each timestep: '//I2S(dtn),Level=5)
      END IF
    END IF
  END IF
  DoTransient = ( dtn > 0 ) 
    
  !------------------------------------------------------------------------
  ! This is a counter for optional case where externally given timestep is
  ! a multitude of internally preferred timestep. 
  dti = 1
1 CONTINUE



  ! First layer (= 0) is determined by the initial conditions (=boundary conditions)
  !-----------------------------------------------------------------------------------
  DO i=1,Mesh % NumberOfNodes
    j = BotPerm(i)
    IF( j > 0 ) THEN
      InvPerm(j) = i
    END IF
  END DO
  PrevInvPerm = InvPerm    

  
  IF( ANY(InvPerm == 0 ) ) THEN
    CALL Fatal(Caller,'Number of nodes has InvPerm undefined: '//I2S(COUNT(InvPerm==0)))
  END IF
  
  CALL GetLayerValues( xvec ) 
  
  ! 0-values at values at the previous layer
  ! The 1st layer cannot really change since it is the BC. 
  CALL GetCoefficients(Set0=.TRUE., xlayer = xvec)
  x0vec = xvec
      
  maxsubt = ListGetInteger( Params,'Timestep Divisions',Found)
  IF(.NOT. Found) maxsubt = 1
  IF( maxsubt > 1 .AND. .NOT. UseInternalVals ) THEN
    CALL Fatal(caller,'We can only have substeps when we use internal values!')
  END IF
  
  
  DO layer=1,NumberOfLayers

    CALL Info(Caller,'Solving for layer: '//I2S(layer),Level=12)
            
    ! Find the next level of nodes, and remember the previous one. 
    PrevInvPerm = InvPerm    
    IF( MaskExist ) THEN
      InvPerm = UpPointer(MaskPerm(PrevInvPerm))
    ELSE
      InvPerm = UpPointer(PrevInvPerm)
    END IF
    IF( InvPerm(1) == PrevInvPerm(1) ) THEN
      CALL Fatal(Caller,'InvPerm is the same on different layers: '//I2S(InvPerm(1)))
    END IF
          
    ! This sets the timestep assuming that all nodes are extruded equally.
    ! Hence this only applied to cartesian drawing. 
    IF( MaskExist ) THEN
      dz = Coord(MaskPerm(InvPerm(1))) - Coord(MaskPerm(PrevInvPerm(1)))
    ELSE
      dz = Coord(InvPerm(1)) - Coord(PrevInvPerm(1))
    END IF
    dtime = dz / velo

    WRITE(Message,'(A,2ES12.3)') 'Layer thickness and timestep: ',dz,dtime
    CALL Info(Caller,Message,Level=12)

    IF(dtime < EPSILON(dtime) ) THEN
      CALL Fatal(Caller,'Cannot march if timestep is ~zero')
    END IF

      ! xi is the value of x at the previous iterate of this layer
#if 1    
    CALL GetLayerValues( xivec ) 
#else
    xivec = xvec
#endif
    !PRINT *,'xivec:',xivec(1)
   
    DO subt=1,maxsubt
      IF( maxsubt > 1 ) THEN
        CALL Info(Caller,'Solving for sub timestep: '//I2S(subt),Level=15)
      END IF


      !IF(subt > 1 ) THEN
      !  xivec = xvec
      !END IF
        
      
      ! We may have iteration if the ODE is nonlinear.
      ! However, more typical could be to iterate over coupled systems. 
      DO iter = 1, MaxIter 
        IF( MaxIter > 1 ) THEN
          CALL Info(Caller,'Nonlinear iteration: '//I2S(iter),Level=20)
          xivec = xvec
        END IF

        CALL GetCoefficients(Beta, xlayer = xvec )

        ! We don't really need any linear solver for this as there is no coupling among dofs.      

        !PRINT *,'x0vec:',x0vec(1),dtime,fvec(1),Beta        
        IF( HaveC ) THEN
          xvec = x0vec + ( dtime / maxsubt ) * fvec / cvec
        ELSE          
          xvec = x0vec + ( dtime / maxsubt ) * fvec
        END IF
        
        IF(SetMin) THEN
          xvec = MAX(xvec,mincons)
        END IF
        IF(SetSum) THEN
          DO i=1,n
            sumcons = SUM(xvec(dofs*(i-1)+1:dofs*i))
            ! The scaling get a little bit more complex when we want to maintain the minimum cuts.
            IF(SetMin) THEN
              xvec(dofs*(i-1)+1:dofs*i) = (1-dofs*mincons)/(sumcons-dofs*mincons)*xvec(dofs*(i-1)+1:dofs*i)
            ELSE
              xvec(dofs*(i-1)+1:dofs*i) = xvec(dofs*(i-1)+1:dofs*i) / sumcons
            END IF
          END DO
        END IF
          
        ! For the 1st iteration the error corresponds to the error with respect to previous solution.
        ! For 2nd etc. iteration the error is of the nonlinear iteration. 
        Norm = SQRT(SUM(xvec*xvec))
        dxvec = xvec-xivec
        Change = SQRT(SUM(dxvec*dxvec)) / Norm

        ! This must be in the loop since we may have dependence on some field value that has changed!
        IF(.NOT. UseInternalVals ) THEN
          CALL SetLayerValues( xvec )
        END IF

        IF( Change < NonLinTol ) EXIT
      END DO
      
      IF(UseInternalVals .AND. subt == maxsubt ) THEN
        CALL SetLayerValues( xvec )
      END IF
      
    
      IF( InfoActive(20) ) THEN
        PRINT *,'Layer:',layer,dtime,Norm,Change      
      END IF
      
      ! 0-values at values at the previous layer
      IF( layer < NumberOfLayers .OR. subt < maxsubt ) THEN
        IF( DoTransient ) THEN        
          ! The only way to have transient simulation is that the timestep is such
          ! that we take exactly one extruded layer. Then the initial value is the
          ! next starting value of the next layer. We have to do some back-and-forth
          ! stuff to have true previous timestep starting values for the coefficients. 
          IF(.NOT. UseInternalVals ) THEN
            CALL SetLayerValues( xivec )
          END IF
          
          ! Set the coefficients at the start at the timestep, hence Set0 = .TRUE. 
          CALL GetCoefficients(Set0=.TRUE.,xlayer=xivec)
          x0vec = xivec
          
          IF(.NOT. UseInternalVals ) THEN
            CALL SetLayerValues( xvec )
          END IF
        ELSE
          ! For steady state the initial value is the final value of this layer. 
          CALL GetCoefficients(Set0=.TRUE.,xlayer=xvec)
          x0vec = xvec
        END IF
      END IF

    END DO
      
  END DO

  IF( dti < dtn ) THEN
    dti = dti + 1
    CALL Info(Caller,'Taking marching step: '//I2S(dti))
    GOTO 1
  END IF
  
  CALL Info(Caller,'All done',Level=5)
  CALL Info(Caller,'-----------------------------------------------------',Level=6)

  
CONTAINS

  ! Finds pointer to the extruded material.
  ! Note that we assume that this is the only material being extruded.
  !--------------------------------------------------------------------------------
  FUNCTION FirstExtrudedMaterial() RESULT ( Material ) 
    TYPE(Valuelist_t), POINTER :: Material
    TYPE(Element_t), POINTER :: Element
    
    Material => NULL()
    Element => GetActiveElement(1) 
    IF( ASSOCIATED( Element ) ) THEN
      Material => GetMaterial(Element)
    END IF
    IF(.NOT. ASSOCIATED( Material ) ) THEN
      CALL Fatal(Caller,'Could not set material for extruded domain!')
    END IF
    
  END FUNCTION FirstExtrudedMaterial
  

  ! Back-substitute values on one single layer back to the distributed field on the
  ! finite element mesh
  !-----------------------------------------------------------------------------------
  SUBROUTINE SetLayerValues( xlayer )
    REAL(KIND=dp) :: xlayer(:)
    INTEGER :: j
    
    IF(dofs == 1 ) THEN
      IF( ParabolicModel ) THEN
        Var3D % Values(Var3D % Perm(InvPerm)) = SQRT( 2 * xlayer(1:m) ) 
      ELSE      
        Var3D % Values(Var3D % Perm(InvPerm)) = xlayer(1:m)
      END IF
    ELSE
      DO j=1,dofs
        IF( ParabolicModel ) THEN
          Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j) = SQRT( 2 * xlayer(j:m:dofs))
        ELSE
          Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j) = xlayer(j:m:dofs)
        END IF
      END DO
    END IF
        
  END SUBROUTINE SetLayerValues

  ! Reverse of the previous
  !-----------------------------------------------------------------------------------
  SUBROUTINE GetLayerValues( xlayer )
    REAL(KIND=dp) :: xlayer(:)
    INTEGER :: j
    
    IF( dofs == 1 ) THEN
      IF( ParabolicModel ) THEN
        xlayer(1:n) = 0.5_dp * Var3D % Values(Var3D % Perm(InvPerm))**2 
      ELSE
        xlayer(1:n) = Var3D % Values(Var3D % Perm(InvPerm))       
      END IF
    ELSE
      DO j=1,dofs
        IF( ParabolicModel ) THEN
          xlayer(j:m:dofs) = 0.5_dp * Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j)**2 
        ELSE
          xlayer(j:m:dofs) = Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j)       
        END IF
      END DO     
    END IF
  END SUBROUTINE GetLayerValues

  
  ! Creates the local matrix equation for the ODY before time integration.
  ! We use namespace in order to allow the same solver to be used for
  ! several fields. The equation is of type:
  ! c*du/dt + r*u = f. 
  !-----------------------------------------------------------------------
  SUBROUTINE GetCoefficients( q, Set0, xlayer )
    REAL(KIND=DP), OPTIONAL :: q
    LOGICAL, OPTIONAL :: Set0
    REAL(KIND=dp), OPTIONAL :: xlayer(:)
    INTEGER :: i,j,k,l
    TYPE(Element_t), POINTER :: Element
    LOGICAL, SAVE :: Visited = .FALSE.
    TYPE(ValueHandle_t), POINTER, SAVE :: Source_h(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: xloc(:)
    REAL(KIND=dp), POINTER :: farray(:,:) => NULL()
    INTEGER, SAVE :: cdofs 
    INTEGER :: interp
    REAL(KIND=dp) :: addv(2), qadd
    
    
    IF( UseInternalVals ) THEN
      IF(.NOT. PRESENT(xlayer) ) THEN
        CALL Fatal(Caller,'"xlayer" must be provided as a parameter!')
      END IF

      
      IF(.NOT. Visited) THEN
        cdofs = dofs
        IF(ASSOCIATED(AddVar)) THEN
          cdofs = cdofs + 1
        END IF
        ALLOCATE(xloc(cdofs))
        CALL Info(Caller,'Get ODE coefficients using '//I2S(cdofs)//' internal variables!')
        IF( dofs == 1 .OR. VectorSource ) THEN
          ALLOCATE(Source_h(1))
          CALL ListInitElementKeyword( Source_h(1),'Material',TRIM(VarName)//': Source', &
              DummyCount = cdofs )
        ELSE
          ALLOCATE(Source_h(dofs))
          DO j=1,dofs
            CALL ListInitElementKeyword( Source_h(j),'Material',TRIM(VarName)//' '//I2S(j)//': Source', &
                DummyCount = cdofs )
          END DO
        END IF
        Visited = .TRUE.
      END IF
      
      DO i=1,n
        j = InvPerm(i)        
        IF(j==0) CALL Fatal('GetCoefficients','We should have positive index!')
        Element => Mesh % Elements( ParentElem(j) )        
        Model % CurrentElement => Element

        ! field values of the additional variable (most likely temperature)
        IF(ASSOCIATED(AddVar)) THEN
          addv = 0.0_dp
          DO interp=1,2
            IF( AddVar % TYPE == Variable_on_nodes_on_elements ) THEN
              DO l=1,Element % TYPE % NumberOfNodes
                IF(Element % NodeIndexes(l) == j) k = Element % DGIndexes(l)
              END DO
            ELSE
              k = j
            END IF
            
            addv(interp) = AddVar % Values(AddVar % Perm(k))
            IF(maxsubt == 1 ) EXIT

            IF(interp == 1) THEN
              j = PrevInvPerm(i)        
            END IF
          END DO

          IF(maxsubt == 1 ) THEN
            qadd = 1.0_dp
          ELSE
            qadd = 1.0_dp * subt / maxsubt
          END IF
          
          xloc(1) = qadd * addv(1) + (1-qadd) * addv(2)
          k = 1
        ELSE
          k = 0
        END IF

        ! the field values of this solver
        DO j=1,dofs
          xloc(k+j) = xlayer(dofs*(i-1)+j)
        END DO
               
        IF( VectorSource ) THEN
          j = 0          
          qadd = ListGetElementReal( Source_h(1), Element = Element, Found = Found, DummyVals = xloc, Rdim=j, Rtensor=farray )
          DO j=1,dofs
            l = dofs*(i-1)+j
            fvec(l:l) = farray(j,1)
          END DO
        ELSE
          DO j=1,dofs
            l = dofs*(i-1)+j
            fvec(l:l) = ListGetElementReal( Source_h(j), Element = Element, Found = Found, DummyVals = xloc )
          END DO
        END IF
      END DO
      
    ELSE IF(MaskExist .OR. AnyDG ) THEN
      k = 1
      DO i=1,n
        j = InvPerm(i)
        SingleIndex(1) = j
        
        IF(j==0) CALL Fatal('GetCoefficients','We should have positive index!')
        Element => Mesh % Elements( ParentElem(j) )
        
        Model % CurrentElement => Element

        !PRINT *,'ParentIndex:',ParentElem(j), j
        !PRINT *,'Nodes:',Element % NodeIndexes
        !PRINT *,'DGs:',Element % DGIndexes

        IF(dofs == 1 ) THEN
          IF( HaveF ) THEN
            fvec(i:i) = ListGetReal( Material,&
                TRIM(VarName)//': Source',k,SingleIndex(1:1))
          END IF
          IF( HaveR ) THEN
            rvec(i:i) = ListGetReal( Material,&
                TRIM(VarName)//': Reaction Coefficient',k,SingleIndex(1:1) )
          END IF
          IF( HaveC ) THEN
            cvec(i:i) = ListGetReal( Material,&
                TRIM(VarName)//': Time Derivative Coefficient',k,SingleIndex(1:1) )
          END IF
        ELSE
          DO j=1,dofs
            l = dofs*(i-1)+j
            IF( HaveF ) THEN
              fvec(l:l) = ListGetReal( Material,&
                  TRIM(VarName)//' '//I2S(j)//': Source',k,SingleIndex(1:1), Found )
            END IF
            IF( HaveR ) THEN
              rvec(l:l) = ListGetReal( Material,&
                  TRIM(VarName)//' '//I2S(j)//': Reaction Coefficient',k,SingleIndex(1:1), Found )
            END IF
            IF( HaveC ) THEN
              cvec(l:l) = ListGetReal( Material,&
                  TRIM(VarName)//' '//I2S(j)//': Time Derivative Coefficient',k,SingleIndex(1:1), Found)
              IF(.NOT. Found) cvec(l:l) = 1.0_dp
            END IF
          END DO
        END IF

        IF(ALL(SingleIndex(1) /= Element % NodeIndexes ) ) STOP
      END DO
    ELSE
      IF( dofs == 1 ) THEN
        IF( HaveF ) THEN
          fvec(1:n) = ListGetReal( Material,&
              TRIM(VarName)//': Source',n,InvPerm )
        END IF
        IF( HaveR ) THEN
          rvec(1:n) = ListGetReal( Material,&
              TRIM(VarName)//': Reaction Coefficient',n,InvPerm )
        END IF
        IF( HaveC ) THEN
          cvec(1:n) = ListGetReal( Material,&
              TRIM(VarName)//': Time Derivative Coefficient',n,invPerm )
        END IF
      ELSE
        DO j=1,dofs
          IF( HaveF ) THEN
            fvec(j:m:dofs) = ListGetReal( Material,&
                TRIM(VarName)//' '//I2S(j)//': Source',n,InvPerm, Found )
          END IF
          IF( HaveR ) THEN
            rvec(j:m:dofs) = ListGetReal( Material,&
                TRIM(VarName)//' '//I2S(j)//': Reaction Coefficient',n,InvPerm, Found )
          END IF
          IF( HaveC ) THEN
            cvec(j:m:dofs) = ListGetReal( Material,&
                TRIM(VarName)//' '//I2S(j)//': Time Derivative Coefficient',n,InvPerm, Found)
            IF(.NOT. Found) cvec(j:m:dofs) = 1.0_dp
          END IF
        END DO
      END IF
        
    END IF

    !PRINT *,'fvec',fvec(1:n)
    !PRINT *,'rvec',rvec(1:n)
    !PRINT *,'cvec',cvec(1:n)
    
    ! When using different integration we may need to access the
    ! value of coefficients at previous mesh layer.
    IF( PRESENT( Set0 ) ) THEN
      IF( Set0 ) THEN
        IF( HaveF ) f0vec = fvec
        IF( HaveR ) r0vec = rvec
        IF( HaveC ) c0vec = cvec
        RETURN
      END IF
    END IF
    
    IF(.NOT. PRESENT( q ) ) RETURN
    IF( ABS(q-1.0_dp) < EPSILON(q) ) RETURN

    IF( HaveF ) fvec = q * fvec + (1-q) * f0vec
    IF( HaveR ) rvec = q * rvec + (1-q) * r0vec
    IF( HaveC ) cvec = q * cvec + (1-q) * c0vec   
    
  END SUBROUTINE GetCoefficients

!------------------------------------------------------------------------------
END SUBROUTINE MarchingODESolver
!------------------------------------------------------------------------------
!> \}
