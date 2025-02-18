!------------------------------------------------------------------------------
!> Calculate lumped fields for vector helmholtz type solver.
!------------------------------------------------------------------------------
 SUBROUTINE VectorHelmholtzLumping(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   USE DefUtils

   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Solver_t) :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------
   TYPE(Variable_t), POINTER :: EVar, PotVar
   TYPE(Element_t), POINTER :: Element
   INTEGER :: i, j, t, k, vdofs, soln, Active, iMode, jMode=0
   TYPE(Solver_t), POINTER :: pSolver
   REAL(KIND=dp) :: mu0inv, eps0, Omega
   CHARACTER(LEN=MAX_NAME_LEN) :: Pname
   LOGICAL :: Found, stat, InitHandles
   TYPE(Mesh_t), POINTER :: Mesh
   COMPLEX(KIND=dp), ALLOCATABLE :: IntPoynt(:,:), IntCurr(:,:), IntVolt(:,:), IntPot(:,:), IntEl(:,:)   
   REAL(KIND=dp), ALLOCATABLE :: IntCenter(:,:), IntWeight(:), IntNorm(:)   
   INTEGER, ALLOCATABLE :: CenterNode(:)
   INTEGER, POINTER :: PotPerm(:)
   REAL(KIND=dp), POINTER :: PotVals(:)
   INTEGER :: EdgeBasisDegree, NoModes
   LOGICAL :: PiolaVersion, NodalMode, EdgeMode
   TYPE(ValueList_t), POINTER :: SolverParams, BC 
   LOGICAL :: Visited = .FALSE.
   CHARACTER(*), PARAMETER :: Caller = 'VectorHelmholtzLumping'

   SAVE IntCenter, IntWeight, IntCurr, IntPot, IntVolt, IntPoynt, IntEl, IntNorm, &
       CenterNode, jMode, NoModes 

   
!-------------------------------------------------------------------------------------------

   CALL Info(Caller,'',Level=6 )
   CALL Info(Caller,'----------------------------------------------------------',Level=6 )
   CALL Info(Caller,'Computing derived fields for electromagnetic wave equation!',Level=4 )
   
   SolverParams => GetSolverParams()

   soln = ListGetInteger( SolverParams,'Primary Solver Index', Found) 
   IF( soln == 0 ) THEN
     CALL Fatal(Caller,'We should know > Primary Solver Index <')
   END IF

   ! Pointer to primary solver
   pSolver => Model % Solvers(soln)

   Mesh => GetMesh()
   
   Omega = GetAngularFrequency(pSolver % Values,Found=Found)
   IF(.NOT. Found) CALL Fatal(Caller,'We need angular frequency!')
   
   Found = .FALSE.
   IF( ASSOCIATED( Model % Constants ) ) THEN
     mu0inv = 1.0_dp / GetConstReal( Model % Constants,'Permeability of Vacuum', Found )
   END IF
   IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
   
   Found = .FALSE.
   IF( ASSOCIATED( Model % Constants ) ) THEN
     eps0 = GetConstReal ( Model % Constants,'Permittivity of Vacuum', Found )
   END IF
   IF(.NOT. Found ) eps0 = 8.854187817d-12   

   IF(.NOT. Visited ) THEN
     NoModes = 0
     DO i=1,Model % NumberOfBCs
       j = ListGetInteger( Model % BCs(i) % Values,'Constraint Mode', Found )
       NoModes = MAX(NoModes, j)
     END DO
     ALLOCATE( IntPoynt(NoModes,NoModes), IntEl(NoModes,Nomodes), IntCurr(NoModes,NoModes), &
         IntPot(NoModes,Nomodes), IntVolt(NoModes,NoModes), IntCenter(NoModes,3), &
         IntWeight(NoModes), CenterNode(NoModes), IntNorm(NoModes) )
     IntPoynt = 0.0_dp
     IntEl = 0.0_dp
     IntCurr = 0.0_dp
     IntVolt = 0.0_dp
     IntPot = 0.0_dp
   END IF
   jMode = jMode + 1

   IF(jMode > NoModes ) THEN
     CALL Fatal(Caller,'The lumping was already called "NoModes" times!')
   END IF

   potVar => VariableGet( Mesh % Variables,'Potential', ThisOnly = .TRUE.)
   
   ! One should be able to toggle between using nodal or edge basis.
   ! The results are not quite the same but should hopefully be close.
   NodalMode = ListGetLogical(SolverParams,'Nodal Target Field',Found )   
   EdgeMode = .NOT. NodalMode

   IF( NodalMode ) THEN  
     CALL Info(Caller,'Assuming electric field living on nodes')

     evar => VariableGet( Mesh % Variables, 'Electric field e', ThisOnly = .TRUE.)
     IF(.NOT. ASSOCIATED( evar ) ) THEN
       evar => VariableGet( Mesh % Variables, 'Electric field', ThisOnly = .TRUE.)
     END IF
     IF(.NOT. ASSOCIATED(evar) ) THEN
       CALL Fatal(Caller,'Could not find nodal electric field!')
     END IF
     IF( evar % dofs /= 6 ) THEN
       CALL Fatal(Caller,'Nodal mode assumes exactly 6 dofs, not '//I2S(evar % dofs))
     END IF
   ELSE
     CALL Info(Caller,'Assuming electric field to live in Hcurl')
     Pname = ListGetString( SolverParams,'Target Variable', Found )
     IF(Found ) THEN
       Evar => VariableGet( pSolver % Mesh % Variables, Pname ) 
     ELSE
       Evar => pSolver % Variable
       Pname = getVarName(pSolver % Variable)
     END IF
     IF(.NOT. ASSOCIATED(evar) ) THEN
       CALL Fatal(Caller,'Could not find electric field living on edges!')
     END IF
     
     CALL Info(Caller,'Name of target variable: '//TRIM(pName),Level=10)

     ! Inherit the solution basis from the primary solver
     vDOFs = Evar % DOFs
     IF( vDofs /= 2 ) CALL Fatal(Caller,'Primary field should have two components!')

     CALL EdgeElementStyle(pSolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
     IF (PiolaVersion) CALL Info(Caller,'Using Piola transformed finite elements',Level=8)
   END IF
  
   
   Active = GetNOFBoundaryElements()

   InitHandles = .TRUE.
   IntCenter = 0.0_dp
   IntWeight = 0.0_dp
   IntNorm = 0.0_dp
   
   DO t=1,Active
     Element => GetBoundaryElement(t)
     BC => GetBC()
     IF (.NOT. ASSOCIATED(BC) ) CYCLE

     SELECT CASE(GetElementFamily())
     CASE(1)
       CYCLE
     CASE(2)
       k = GetBoundaryEdgeIndex(Element,1); Element => Mesh % Edges(k)
     CASE(3,4)
       k = GetBoundaryFaceIndex(Element)  ; Element => Mesh % Faces(k)
     END SELECT
     
     iMode = ListGetInteger( BC,'Constraint Mode',Found )
     IF( iMode == 0 ) CYCLE       

     CALL LocalIntegBC(BC,Element,InitHandles )
   END DO

   IF( ParEnv % PEs > 1 ) THEN
     DO i=1,NoModes
       IntPoynt(jMode,i) = ParallelReduction(IntPoynt(jMode,i))
       IntEl(jMode,i) = ParallelReduction(IntEl(jMode,i))
       IntCurr(jMode,i) = ParallelReduction(IntCurr(jMode,i))
       IntWeight(i) = ParallelReduction(IntWeight(i))
       IntNorm(i) = ParallelReduction(IntNorm(i))
     END DO
   END IF

   PRINT *,'Elfield before norm through port:',IntEl(jMode,:)
   PRINT *,'Elfield norm:',IntNorm(:)

   ! Why sqrt(2) ? 
   IntEl(jMode,:) = IntEl(jMode,:) / IntNorm(:)
   IntEl(jMode,jMode) = IntEl(jMode,jMode) - 1.0_dp

   PRINT *,'Elfield itse through port:',IntEl(jMode,:)
   PRINT *,'Elfield area port:',IntWeight(jMode)
   PRINT *,'Elfield Energy through port:',IntPoynt(jMode,:)
   PRINT *,'Elfield Current through port:',IntCurr(jMode,:)
   IF( ASSOCIATED(PotVar) ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       DO i=1,NoModes
         IntPot(jMode,i) = ParallelReduction(IntPot(jMode,i))
       END DO
     END IF
     IntPot(jMode,:) = IntPot(jMode,:) / IntWeight(:)
     PRINT *,'Average pot on port:',IntPot(jMode,:)
   END IF
   
   CALL CenterPortLoc()
   
   IF( ASSOCIATED( potVar ) ) THEN
     CALL Info(Caller,'Computing voltage using "Potential" as path indicator')
     potVals => PotVar % Values
     PotPerm => PotVar % Perm        
     CALL EdgeVoltageIntegral()
   ELSE
     CALL Info(Caller,'Cannot compute voltage as no "Potential" is present')
   END IF
     
   IF( jMode == NoModes .AND. ParEnv % MyPe == 0 ) THEN
     CALL Info(Caller,'Writing results on final visit!')

     OPEN (10, FILE="Poynt_re.dat")
     DO i=1,NoModes
       WRITE(10,*) REAL(IntPoynt(i,:))
     END DO
     CLOSE(10) 
     OPEN (10, FILE="Poynt_im.dat")
     DO i=1,NoModes
       WRITE(10,*) AIMAG(IntPoynt(i,:))
     END DO
     CLOSE(10)
     OPEN (10, FILE="El_re.dat")
     DO i=1,NoModes
       WRITE(10,*) REAL(IntEl(i,:))
     END DO
     CLOSE(10) 
     OPEN (10, FILE="El_im.dat")
     DO i=1,NoModes
       WRITE(10,*) AIMAG(IntEl(i,:))
     END DO
     CLOSE(10)
     OPEN (10, FILE="Curr_re.dat")
     DO i=1,NoModes
       WRITE(10,*) REAL(IntCurr(i,:))
     END DO
     CLOSE(10) 
     OPEN (10, FILE="Curr_im.dat")
     DO i=1,NoModes
       WRITE(10,*) AIMAG(IntCurr(i,:))
     END DO
     CLOSE(10) 
     IF( ASSOCIATED( PotVar ) ) THEN
       OPEN (10, FILE="dPot_re.dat")
       DO i=1,NoModes
         WRITE(10,*) REAL(IntPot(i,:))
       END DO
       CLOSE(10) 
       OPEN (10, FILE="dPot_im.dat")
       DO i=1,NoModes
         WRITE(10,*) AIMAG(IntPot(i,:))
       END DO
       CLOSE(10)
       OPEN (10, FILE="Volt_re.dat")
       DO i=1,NoModes
         WRITE(10,*) REAL(IntVolt(i,:))
       END DO
       CLOSE(10) 
       OPEN (10, FILE="Volt_im.dat")
       DO i=1,NoModes
         WRITE(10,*) AIMAG(IntVolt(i,:))
       END DO
       CLOSE(10)
     END IF
   END IF
        
   Visited = .TRUE.
   CALL Info(Caller,'All done for now!',Level=20)   
   
   
CONTAINS

  
!-----------------------------------------------------------------------------
  SUBROUTINE LocalIntegBC( BC, Element, InitHandles )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: B, Zs, L(3), muinv, TemGrad(3), BetaPar, jn, eps, &
        e_ip(3), e_ip_norm, e_ip_tan(3), f_ip_tan(3), imu, phi
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),WBasis(:,:),RotWBasis(:,:), e_local(:,:), phi_local(:,:)
    REAL(KIND=dp) :: weight, DetJ, Normal(3), cond, u, v, w, x, y, z
    LOGICAL :: Stat, Found
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, m, np, p, q, ndofs, n, nd   
    TYPE(Nodes_t), SAVE :: Nodes, ParentNodes
    LOGICAL :: AllocationsDone = .FALSE.
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueHandle_t), SAVE :: MagLoad_h, ElRobin_h, MuCoeff_h, Absorb_h, TemRe_h, TemIm_h
    TYPE(ValueHandle_t), SAVE :: TransferCoeff_h, ElCurrent_h, RelNu_h, CondCoeff_h, CurrDens_h, EpsCoeff_h
     
    SAVE AllocationsDone, WBasis, RotWBasis, Basis, dBasisdx, e_local, phi_local

    ndofs = evar % dofs
    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDOFs
      ALLOCATE( WBasis(m,3), RotWBasis(m,3), Basis(m), dBasisdx(m,3), e_local(ndofs,m), phi_local(2,m) )      
      AllocationsDone = .TRUE.
    END IF
 
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( ElRobin_h,'Boundary Condition','Electric Robin Coefficient',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MagLoad_h,'Boundary Condition','Magnetic Boundary Load', InitIm=.TRUE.,InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( Absorb_h,'Boundary Condition','Absorbing BC')
      CALL ListInitElementKeyword( TemRe_h,'Boundary Condition','TEM Potential')
      CALL ListInitElementKeyword( TemIm_h,'Boundary Condition','TEM Potential Im')
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( TransferCoeff_h,'Boundary Condition','Electric Transfer Coefficient',InitIm=.TRUE.)
      CALL ListInitElementKeyword( ElCurrent_h,'Boundary Condition','Electric Current Density',InitIm=.TRUE.)
      CALL ListInitElementKeyword( CurrDens_h,'Body Force','Current Density', InitIm=.TRUE.,InitVec3D=.TRUE.)      
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
      InitHandles = .FALSE.
    END IF

    imu = CMPLX(0.0_dp, 1.0_dp)

    n = Element % TYPE % NumberOfNodes
    CALL GetElementNodes( Nodes, Element )    
    Parent => GetBulkElementAtBoundary(Element)
    IF(.NOT. ASSOCIATED( Parent ) ) THEN
      CALL Fatal(Caller,'Model lumping requires parent element!')
    END IF

    e_local = 0.0_dp

    IF( ASSOCIATED( PotVar ) ) THEN        
      phi_local = 0.0_dp
      CALL GetVectorLocalSolution( phi_local, uelement = Element, uvariable = potvar )
    END IF
    
    IF( NodalMode ) THEN
      CALL GetVectorLocalSolution( e_local, uelement = Element, uvariable = evar, Found=Found)
    ELSE
      CALL GetElementNodes( ParentNodes, Parent )    
      CALL GetVectorLocalSolution( e_local, UElement = Parent, Uvariable=eVar, &
          uSolver=pSolver, Found=Found)
      np = n * MAXVAL(Solver % Def_Dofs(GetElementFamily(Parent),:,1))
      np = 0
      nd = GetElementNOFDOFs(Parent,uSolver=pSolver)
    END IF
    IF(.NOT. Found) THEN
      CALL Fatal(Caller,'Could not find field data on boundary!?')
    END IF
            
    Normal = NormalVector(Element, Nodes, Check=.TRUE.)
    
    ! Numerical integration:
    !-----------------------
    IF( NodalMode ) THEN
      IP = GaussPoints(Element)
    ELSE
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)
    END IF

    DO t=1,IP % n  
      
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )              
      weight = IP % s(t) * detJ

      B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )
            
      ! Get material properties from parent element.
      !----------------------------------------------
      muinv = ListGetElementComplex( MuCoeff_h, Basis, Parent, Found, GaussPoint = t )      
      IF( Found ) THEN
        muinv = muinv * mu0inv
      ELSE
        muinv = mu0inv
      END IF

      eps = ListGetElementComplex( EpsCoeff_h, Basis, Parent, Found, GaussPoint = t )      
      IF( Found ) THEN
        eps = eps * eps0
      ELSE
        eps = eps0
      END IF

      Cond = ListGetElementReal( CondCoeff_h, Basis, Parent, Found, GaussPoint = t )

      Zs = imu * Omega / (B * muinv)
      
      IntWeight(iMode) = IntWeight(iMode) + weight
      x = SUM(Basis(1:n) * Nodes % x(1:n)) 
      y = SUM(Basis(1:n) * Nodes % y(1:n)) 
      z = SUM(Basis(1:n) * Nodes % z(1:n)) 
      
      IntCenter(iMode,1) = IntCenter(iMode,1) + weight * x
      IntCenter(iMode,2) = IntCenter(iMode,2) + weight * y
      IntCenter(iMode,3) = IntCenter(iMode,3) + weight * z

      L = ListGetElementComplex3D( MagLoad_h, Basis, Element, Found, GaussPoint = t )

      TemGrad = CMPLX( ListGetElementRealGrad( TemRe_h,dBasisdx,Element,Found), &
          ListGetElementRealGrad( TemIm_h,dBasisdx,Element,Found) )
      L = L + TemGrad
      
      B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )

      L = L / (2*SQRT(2.0_dp)*B)

      IF( NodalMode ) THEN
        DO i=1,3
          e_ip(i) = CMPLX( SUM( Basis(1:n) * e_local(i,1:n) ), SUM( Basis(1:n) * e_local(i+3,1:n) ) )
        END DO
      ELSE        
        ! In order to get the normal component of the electric field we must operate on the
        ! parent element. The surface element only has tangential components. 
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp
        CALL FindParentUVW( Element, n, Parent, Parent % TYPE % NumberOfNodes, U, V, W, Basis )
        IF (GetElementFamily(Element) == 2) THEN
          stat = EdgeElementInfo(Parent, ParentNodes, u, v, w, detF = detJ, &
              Basis = Basis, EdgeBasis = Wbasis, RotBasis = RotWBasis, dBasisdx = dBasisdx, &
              BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
        ELSE    
          stat = ElementInfo( Parent, ParentNodes, u, v, w, detJ, Basis, dBasisdx, &
              EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver )
        END IF
        e_ip(1:3) = CMPLX(MATMUL(e_local(1,np+1:nd),WBasis(1:nd-np,1:3)), MATMUL(e_local(2,np+1:nd),WBasis(1:nd-np,1:3)))       
      END IF
            
      e_ip_norm = SUM(e_ip*Normal)
      e_ip_tan = e_ip - e_ip_norm * Normal

      IntNorm(iMode) = IntNorm(iMode) + weight * ABS( SUM( L * CONJG(L) ) )
      IntEl(jMode,iMode) = IntEl(jMode,iMode) + weight * SUM(e_ip_tan * CONJG(L) ) 
      
      IntPoynt(jMode,iMode) = IntPoynt(jMode,iMode) + weight * &
          0.5_dp * SUM(e_ip_tan * CONJG(e_ip_tan) ) / Zs
      IntCurr(jMode,iMode) = IntCurr(jMode,iMode) + weight * &
          ( imu * Omega * Eps + cond ) * e_ip_norm 

      ! If potential is given compute the integral over the potential on the electrode to get
      ! the avarage potential at the surface.
      IF( ASSOCIATED( PotVar ) ) THEN        
        CALL GetVectorLocalSolution( phi_local, uelement = Element, uvariable = potvar )
        phi = CMPLX( SUM( Basis(1:n) * phi_local(1,1:n) ), SUM( Basis(1:n) * phi_local(2,1:n) ) )
        IntPot(jMode,iMode) = IntPot(jMode,iMode) + weight * phi
      END IF

    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalIntegBC
!------------------------------------------------------------------------------


  ! Find a node closest to the port center.
  !-----------------------------------------
  SUBROUTINE CenterPortLoc()
    
    INTEGER :: i,iMode,mini
    REAL(KIND=dp) :: mindist2, dist2, Coord0(3), Coord1(3)    

    PRINT *,'IntWeight:',IntWeight
    
    ! Calculate the center nodes for each mode.
    IF( ANY( IntWeight < EPSILON(dist2) ) ) THEN
      PRINT *,'IntWeight:',IntWeight
      CALL Fatal(Caller,'Some weight is zero!')
    END IF

    DO i=1,3
      IntCenter(:,i) = IntCenter(:,i) / IntWeight(:)
    END DO

        
    ! Find the node closest to the center for each port
    ! As each port is planar the node with minimum distance
    ! should hopefully lie on the port as well.  
    CenterNode = 0
    DO iMode=1,NoModes
      Coord0 = IntCenter(iMode,:)
      PRINT *,'Center:',Coord0
            
      mindist2 = HUGE(mindist2)
      DO i=1,Mesh % NumberOfNodes
        Coord1(1) = Mesh % Nodes % x(i)
        Coord1(2) = Mesh % Nodes % y(i)
        Coord1(3) = Mesh % Nodes % z(i)        
        dist2 = SUM((Coord0-Coord1)**2)
        IF(dist2 < mindist2 ) THEN
          mindist2 = dist2
          mini = i
        END IF
      END DO

      PRINT *,'Minimum distance:',iMode,mini,SQRT(mindist2)
      
      CenterNode(iMode) = mini      
    END DO
    
  END SUBROUTINE CenterPortLoc
    

  ! Perform line integral from the center of port to ground. The ground is defined as a node having
  ! the smallest potential that can be reached following the edges. The integral is taking over this
  ! route. 
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE EdgeVoltageIntegral()

    TYPE(Matrix_t), POINTER :: NodeGraph
    COMPLEX :: gradv(3), Circ
    REAL(KIND=dp) :: pot, minpot, EdgeVector(3), s
    INTEGER :: iMode, nsteps, sgn, i, j, k, kmin, imin, i1, i2, j1, j2, l, lp, n0
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Edge
    
    ! Create a graph for node-to-edge connectivity
    !----------------------------------------------
    NodeGraph => AllocateMatrix()
    NodeGraph % FORMAT = MATRIX_LIST         
    DO i = Mesh % NumberOfEdges, 1, -1
      Edge => Mesh % Edges(i)
      DO j=1, Edge % TYPE % NumberOfNodes 
        CALL List_AddToMatrixElement( NodeGraph % ListMatrix,Edge % NodeIndexes(j),i,1.0_dp )
      END DO
    END DO
    CALL List_ToCRSMatrix(NodeGraph)
    PRINT *,'Nonzeros per row NodeGraph:',1.0_dp * SIZE(NodeGraph % Values) / NodeGraph % NumberOfRows
    n0 = Mesh % NumberOfNodes
    
    DO iMode=1,NoModes
      i = CenterNode(iMode)
      minpot = PotVals(PotPerm(i))
      nsteps = 0
      Circ = 0.0_dp
      
      PRINT *,'Starting pot:',iMode, minpot, IntPot(jMode,iMode) 
      
      DO WHILE(.TRUE.)
        kmin = 0
        ! Among the edges related to node "i" find the one that has the steepest
        ! potential descent.
        DO j = NodeGraph % Rows(i),NodeGraph % Rows(i+1)-1
          k = NodeGraph % Cols(j)
          Edge => Mesh % Edges(k)
          NodeIndexes => Edge % NodeIndexes
          DO l=1,2
            IF(NodeIndexes(l) == i) CYCLE
            lp = PotPerm(NodeIndexes(l))
            IF(lp == 0) CYCLE
            pot = PotVals(lp)
            IF( pot < minpot ) THEN
              kmin = k
              imin  = NodeIndexes(l)
              minpot = pot              
            END IF
          END DO
        END DO

        ! When no smaller potential is found we are done.
        IF( kmin == 0 ) EXIT

        nsteps = nsteps + 1
        
        ! The edge to integrate over
        k = kmin

        ! Edge that goes to the minimum value
        Edge => Mesh % Edges(k)

        i1 = Edge % NodeIndexes(1)
        i2 = Edge % NodeIndexes(2)
                
        EdgeVector(1) = Mesh % Nodes % x(i2) - Mesh % Nodes % x(i1)
        EdgeVector(2) = Mesh % Nodes % y(i2) - Mesh % Nodes % y(i1)
        EdgeVector(3) = Mesh % Nodes % z(i2) - Mesh % Nodes % z(i1)

        ! Integration length and direction
        s = SQRT(SUM(EdgeVector**2))

        ! If we do the path integral in the wrong direction compared to definiotion of edge switch the sign
        sgn = 1
        IF(i /= i1 ) sgn = -sgn 

        IF( NodalMode ) THEN
          j1 = eVar % Perm(i1)
          j2 = eVar % Perm(i2)          
          DO k=1,3
            gradv(k) = CMPLX(eVar % Values(6*(j1-1)+k) + eVar % Values(6*(j2-1)+k),&
                eVar % Values(6*(j1-1)+3+k) + eVar % Values(6*(j2-1)+3+k) ) / 2
          END DO
          Circ = Circ + sgn * SUM(gradv*EdgeVector)
        ELSE        
          ! Check the sign if the direction based on global edge direction rules
          IF( ParEnv % PEs > 1 ) THEN                            
            i1 = Mesh % ParallelInfo % GlobalDOFs(i1)             
            i2 = Mesh % ParallelInfo % GlobalDOFs(i2)             
          END IF
          IF( i1 < i2) sgn = -sgn      
                                      
          j = eVar % Perm(n0 + k)
          Circ = Circ + s * sgn * CMPLX( eVar % Values(2*j-1),eVar % Values(2*j) )
        END IF

        ! Continue from the end point
        i = imin
      END DO

      PRINT *,'Path integral:',iMode, minpot, nsteps, Circ
      IntVolt(jMode,iMode) = Circ
    END DO

    CALL FreeMatrix(NodeGraph)
    
  END SUBROUTINE EdgeVoltageIntegral
    
!------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzLumping
!------------------------------------------------------------------------

