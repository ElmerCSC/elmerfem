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
! *  Original Date: 30. April 2010
! * 
! *****************************************************************************
SUBROUTINE AdjointSSA_GradientSolver_init0(Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: Name

  Name = ListGetString( Solver % Values, 'Equation',UnFoundFatal=.TRUE.)
  CALL ListAddNewString( Solver % Values,'Variable',&
          '-nooutput '//TRIM(Name)//'_var')
  CALL ListAddLogical(Solver % Values, 'Optimize Bandwidth',.FALSE.)

END SUBROUTINE AdjointSSA_GradientSolver_init0
! *****************************************************************************
SUBROUTINE AdjointSSA_GradientSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!   Compute the Gradient of user defined cost functions with respect to various
!   input parameters of the SSA model:
!
!    
!     OUTPUT is : (OPTIONAL) DJDBeta (gradient wr to slip coeff)
!                 (OPTIONAL) DJDZs (gradient wr to Zs)      
!                 (OPTIONAL) DJDZb (gradient wr to Zb) 
!                 (OPTIONAL) DJDRho (gradient wr to Mean Density)
!                 (OPTIONAL) DJDEta (gradient wr to Mean Viscosity)
!
!     BE careful: - change of variable (for example slip coeff = 10^alpha) has to
!     be taken care by the user in the .sif
!                 - by default the gradient is reset to 0 here; set "Reset DJD... = Logical False" 
!     if part of the gradient has already been computed before 
!
!
!     INPUT PARAMETERS are (in addition to the SSA required INPUTS):
!
!      In solver section:
!             Flow Solution Name = String (default SSAVelocity)
!             Adjoint Solution Name = String (default Adjoint)
!
!             Compute DJDBeta = Logical (default False)
!             Reset DJDBeta = Logical (default True)
!
!             Compute DJDZs = Logical (default False)
!             Reset DJDZs = Logical (default True)
!
!             Compute DJDZb = Logical (default False)  
!             Reset DJDZb = Logical (default True)
!
!             Compute DJDRho = Logical (default False)
!             Reset DJDRho = Logical (default True)
!
!             Compute DJDEta = Logical (default False)
!             Reset DJDEta = Logical (default True)
!                
!
!
!******************************************************************************
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
  TYPE(Solver_t), POINTER :: DSolver
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, BC
  TYPE(Variable_t), POINTER ::  ZsSol, ZbSol, &
                               VeloSolN,VeloSolD,&
                               DJDBetaSol,DJDZsSol,DJDZbSol,DJDRhoSol,DJDEtaSol

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt, CalvingFront 
  LOGICAL :: Newton

  INTEGER :: i, n, m, t, istat, DIM, p, STDOFs
  INTEGER :: NonlinearIter, NewtonIter, iter, other_body_id
          
  INTEGER, POINTER :: ZsPerm(:), ZbPerm(:), &
       VeloNPerm(:),VeloDPerm(:),&
       DJDBetaPerm(:),DJDZsPerm(:),DJDZbPerm(:), DJDRhoPerm(:),DJDEtaPerm(:),&
       NodeIndexes(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: Zs(:), Zb(:)
  REAL(KIND=dp), POINTER :: DJDBeta(:),DJDZs(:),DJDZb(:),DJDRho(:),DJDEta(:),VelocityN(:),VelocityD(:)
                            
  REAL(KIND=dp) :: UNorm, cn, dd, NonlinearTol, NewtonTol, MinSRInv, rhow, sealevel, &
                   PrevUNorm, relativeChange,minv

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalGravity(:), NodalViscosity(:), NodalDensity(:), &
           NodalZs(:), NodalZb(:), NodalGM(:),NodalBed(:),  &
           NodalU(:), NodalV(:), NodalBeta(:),LocalLinVelo(:),&
           Nodalbetab(:),Nodalzsb(:),Nodalzbb(:),NodalRhob(:),NodalEtab(:),&
           NodalEtaDer(:),NodalBetaDer(:)

  INTEGER :: iFriction
  REAL(KIND=dp) :: fm,U0
  CHARACTER(LEN=MAX_NAME_LEN) :: Friction
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='AdjointSSA_GradientSolver'
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: at, at0
#else
    REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif 
  LOGICAL :: SEP ! Sub-element parametrization for Grounding line
  INTEGER :: GLnIP ! number of Integ. Points for GL Sub-element parametrization
  TYPE(Variable_t), POINTER :: GMSol,BedrockSol

  LOGICAL , SAVE :: ComputeDJDBeta,ComputeDJDZs,ComputeDJDZb,ComputeDJDRho,ComputeDJDEta,Reset
  LOGICAL :: HaveBetaDer,HaveEtaDer
  CHARACTER(LEN=MAX_NAME_LEN), SAVE :: NeumannSolName,AdjointSolName,SName
  INTEGER,SAVE :: SolverInd

  SAVE rhow,sealevel
  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName, ElementNodes
  SAVE NodalGravity, NodalViscosity, NodalDensity, &
           NodalZs, NodalZb, NodalGM,NodalBed,  &
           NodalU, NodalV, NodeIndexes, NodalBeta,LocalLinVelo, &
           Nodalbetab,Nodalzsb,NodalZbb,NodalRhob,NodalEtab,&
           NodalEtaDer,NodalBetaDer
  SAVE STDOFs

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()




!!!!!!!!!!! get Solver Variables
  SolverParams => GetSolverParams()
  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

    NeumannSolName =  GetString( SolverParams,'Flow Solution Name', Found)
    IF(.NOT.Found) THEN        
            CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Solver<')
            CALL WARN(SolverName,'Taking default value >SSAVelocity<')
            WRITE(NeumannSolName,'(A)') 'SSAVelocity'
    END IF
  !! SSA Solution
     VeloSolN => VariableGet( Solver % Mesh % Variables, TRIM(NeumannSolName), UnFoundFatal=.TRUE.   )
     STDOFs = VeloSolN % DOFs

     ! Get the Direct solver
     DO i=1,Model % NumberOfSolvers
        if (TRIM(NeumannSolName) == TRIM(Model % Solvers(i) % Variable % Name)) exit
     End do
     if (i.eq.(Model % NumberOfSolvers+1)) CALL FATAL(SolverName,'Could not find Flow Solver Equation Name')
     SolverInd=i

     AdjointSolName =  GetString( SolverParams,'Adjoint Solution Name', Found)
     IF(.NOT.Found) THEN        
         CALL WARN(SolverName,'Keyword >Adjoint Solution Name< not found in section >Solver<')
         CALL WARN(SolverName,'Taking default value >Adjoint<')
         WRITE(AdjointSolName,'(A)') 'Adjoint'
     END IF
     ComputeDJDBeta =  GetLogical( SolverParams,'Compute DJDBeta', Found)
     IF(.NOT.Found) ComputeDJDBeta=.False.
     ComputeDJDZs =  GetLogical( SolverParams,'Compute DJDZs', Found)
     IF(.NOT.Found) ComputeDJDZs=.False.
     ComputeDJDZb =  GetLogical( SolverParams,'Compute DJDZb', Found)
     IF(.NOT.Found) ComputeDJDZb=.False.
     ComputeDJDRho =  GetLogical( SolverParams,'Compute DJDRho', Found)
     IF(.NOT.Found) ComputeDJDRho=.False.
     ComputeDJDEta =  GetLogical( SolverParams,'Compute DJDEta', Found)
     IF(.NOT.Found) ComputeDJDEta=.False.


     ! Get some constants
     rhow = ListGetConstReal( Model % Constants, 'Water Density', UnFoundFatal=.TRUE.)

     sealevel = ListGetConstReal( Model % Constants, 'Sea Level', UnFoundFatal=.TRUE. )

     ! Allocate
     N = Model % MaxElementNodes
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalGravity, &
                       NodalViscosity, NodalDensity,  &
                       NodalZb, NodalZs, NodalGM,NodalBed, NodalU, NodalV, &
                       NodalBeta,LocalLinVelo, Nodalbetab,Nodalzsb,NodalRhob,&
                       NodalEtab,NodalZbb,&
                       NodalBetaDer,NodalEtaDer,&
                       ElementNodes % x, &
                       ElementNodes % y, ElementNodes % z )

     ALLOCATE( FORCE(STDOFs*N), LOAD(N), STIFF(STDOFs*N,STDOFs*N), &
          NodalGravity(N), NodalDensity(N), NodalViscosity(N), &
          NodalZb(N), NodalZs(N) , NodalGM(N),NodalBed(N),&
          NodalU(N), NodalV(N), NodalBeta(N), LocalLinVelo(N),&
          Nodalbetab(N),NodalZsb(N), NodalRhob(N),NodalEtab(N),&
          NodalZbb(N),&
          NodalBetaDer(N),NodalEtaDer(N),&
          ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N), &
           STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

! Get Info from Direct Solver Params
  DSolver => Model % Solvers(SolverInd)
!  Sub - element GL parameterisation
  SEP=GetLogical( DSolver % Values, 'Sub-Element GL parameterization',GotIt)
  IF (.NOT.GotIt) SEP=.False.
  IF (SEP) THEN
     GLnIP=ListGetInteger( DSolver % Values, &
           'GL integration points number',UnFoundFatal=.TRUE. )
     WRITE(Message,'(a,i0,a)') 'Sub-Element GL parameterization using ',GLnIP,' IPs'
     CALL INFO(SolverName,TRIM(Message),level=4)
     GMSol => VariableGet( Solver % Mesh % Variables, 'GroundedMask',UnFoundFatal=.TRUE. )
     BedrockSol => VariableGet( Solver % Mesh % Variables, 'bedrock',UnFoundFatal=.TRUE. )
  END IF

  ZbSol => VariableGet( Solver % Mesh % Variables, 'Zb', UnFoundFatal=.TRUE. )
  Zb => ZbSol % Values
  ZbPerm => ZbSol % Perm

  ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs', UnFoundFatal=.TRUE.  )
  Zs => ZsSol % Values
  ZsPerm => ZsSol % Perm

  !! SSA Solution
  VeloSolN => VariableGet( Solver % Mesh % Variables, NeumannSolName, UnFoundFatal=.TRUE.  )
  VelocityN => VeloSolN % Values
  VeloNPerm => VeloSolN % Perm

  !! Adjoint Solution
  VeloSolD => VariableGet( Solver % Mesh % Variables, AdjointSolName, UnFoundFatal=.TRUE. )
  VelocityD => VeloSolD % Values
  VeloDPerm => VeloSolD % Perm

  IF (ComputeDJDBeta) Then
     SName =  GetString( SolverParams,'DJDBeta Name', Found)
     IF(.NOT.Found) THEN        
            CALL WARN(SolverName,'Keyword >DJDBeta Name< not found in section >Solver<')
            CALL WARN(SolverName,'Taking default value >DJDBeta<')
            WRITE(SName,'(A)') 'DJDBeta'
            CALL ListAddString(  SolverParams, 'DJDBeta Name', TRIM(SName))
      END IF
      DJDBetaSol => VariableGet( Solver % Mesh % Variables, TRIM(SName) ,UnFoundFatal=.TRUE. )
      DJDBeta => DJDBetaSol % Values
      DJDBetaPerm => DJDBetaSol % Perm

      Reset =  GetLogical( SolverParams,'Reset DJDBeta', Found)
      if (Reset.OR.(.NOT.Found)) DJDBeta = 0.0
  End if

  IF (ComputeDJDZs) Then
       DJDZsSol => VariableGet( Solver % Mesh % Variables, 'DJDZs',UnFoundFatal=.TRUE. )
       DJDZs => DJDZsSol % Values
       DJDZsPerm => DJDZsSol % Perm
          
       Reset =  GetLogical( SolverParams,'Reset DJDZs', Found)
       if (Reset.OR.(.NOT.Found)) DJDZs = 0.0
  End if
  IF (ComputeDJDZb) Then
       DJDZbSol => VariableGet( Solver % Mesh % Variables, 'DJDZb', UnFoundFatal=.TRUE. )
       DJDZb => DJDZbSol % Values
       DJDZbPerm => DJDZbSol % Perm
          
       Reset =  GetLogical( SolverParams,'Reset DJDZb', Found)
       if (Reset.OR.(.NOT.Found)) DJDZb = 0.0
  End if

  IF (ComputeDJDRho) Then
        DJDRhoSol => VariableGet( Solver % Mesh % Variables, 'DJDRho', UnFoundFatal=.TRUE. )
        DJDRho => DJDRhoSol % Values
        DJDRhoPerm => DJDRhoSol % Perm

        Reset =  GetLogical( SolverParams,'Reset DJDRho', Found)
        if (Reset.OR.(.NOT.Found)) DJDRho = 0.0
  End if
  IF (ComputeDJDEta) Then
        DJDEtaSol => VariableGet( Solver % Mesh % Variables, 'DJDEta' ,UnFoundFatal=.TRUE.)
        DJDEta => DJDEtaSol % Values
        DJDEtaPerm => DJDEtaSol % Perm

        Reset =  GetLogical( SolverParams,'Reset DJDEta', Found)
        if (Reset.OR.(.NOT.Found)) DJDEta = 0.0
  END IF

  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes(Element)

     NodeIndexes => Element % NodeIndexes

 ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (STDOFs == 1) THEN !1D SSA
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (STDOFs == 2) THEN !2D SSA
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                STDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

     ! Read the gravity in the Body Force Section 
     BodyForce => GetBodyForce()
     NodalGravity = 0.0_dp
     IF ( ASSOCIATED( BodyForce ) ) THEN
           IF (STDOFs==1) THEN 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 2', n, NodeIndexes, UnFoundFatal=.TRUE.)
           ELSE 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 3', n, NodeIndexes, UnFoundFatal=.TRUE.)
           END IF
     END IF

     ! Read the Viscosity eta, density, and exponent m in MMaterial Section
     ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
     Material => GetMaterial(Element)
     cn = ListGetConstReal( Material, 'Viscosity Exponent',UnFoundFatal=.TRUE.)
     MinSRInv = ListGetConstReal( Material, 'Critical Shear Rate',UnFoundFatal=.TRUE.)

     NodalDensity=0.0_dp
     NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n,NodeIndexes,UnFoundFatal=.TRUE.)

     NodalViscosity=0.0_dp
     NodalViscosity(1:n) = ListGetReal( Material, 'SSA Mean Viscosity',n, NodeIndexes,UnFoundFatal=.TRUE.)
     NodalEtaDer(1:n) = ListGetReal( Material, 'SSA Mean Viscosity Derivative',n, NodeIndexes,Found=HaveEtaDer)

     Friction = ListGetString(Material, 'SSA Friction Law', UnFoundFatal=.TRUE.)

    SELECT CASE(Friction)
       CASE('linear')
         iFriction = 1
         fm = 1.0_dp
       CASE('weertman')
        iFriction = 2
       CASE('regularized coulomb')
        iFriction = 3
       CASE DEFAULT
         CALL FATAL(SolverName,'Friction should be linear or Weertman or regularized coulomb')
   END SELECT


   NodalBeta(1:n) = ListGetReal( Material, 'SSA Friction Parameter',n, NodeIndexes,UnFoundFatal=.TRUE.)
   NodalBetaDer(1:n) = ListGetReal( Material, 'SSA Friction Parameter Derivative',n, NodeIndexes,Found=HaveBetaDer)
   IF (iFriction > 1) THEN
        fm = ListGetConstReal( Material, 'SSA Friction Exponent', UnFoundFatal=.TRUE. )

        LocalLinVelo = 0.0_dp
        LocalLinVelo(1:n) = ListGetReal(Material, 'SSA Friction Linear Velocity', n, NodeIndexes,UnFoundFatal=.TRUE.)
   END IF

   IF (iFriction == 3) THEN
        U0 = ListGetConstReal( Material, 'SSA Friction Threshold Velocity', Found, UnFoundFatal=.TRUE.)
   END IF


   IF (SEP) THEN
     NodalGM(1:n)=GMSol%Values(GMSol%Perm(NodeIndexes(1:n)))
     NodalBed(1:n)=BedrockSol%Values(BedrockSol%Perm(NodeIndexes(1:n)))
   ENDIF

   ! Get the Nodal value of Zb and Zs
   NodalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
   NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))

   ! Previous Velocity 
   NodalU(1:n) = VelocityN(STDOFs*(VeloNPerm(NodeIndexes(1:n))-1)+1)
   NodalV = 0.0
   IF (STDOFs.EQ.2) NodalV(1:n) = VelocityN(STDOFs*(VeloNPerm(NodeIndexes(1:n))-1)+2)
      
   CALL LocalMatrixUVSSA (  STIFF, FORCE, Element, n, ElementNodes, NodalGravity, &
        NodalDensity, NodalViscosity, NodalEtaDer,HaveEtaDer,NodalZb, NodalZs, &
        NodalU, NodalV, NodalBeta,NodalBetaDer,HaveBetaDer,iFriction,fm,LocalLinVelo, cn, &
        NodalGM,NodalBed,SEP,GLnIP,sealevel,rhow,&
        MinSRInv , STDOFs,&
        Nodalbetab,nodalzsb,nodalzbb,nodalrhob,nodaletab)

   IF (ComputeDJDBeta) &
        DJDBeta(DJDBetaPerm(NodeIndexes(1:n)))=DJDBeta(DJDBetaPerm(NodeIndexes(1:n)))+Nodalbetab(1:n)
   IF (ComputeDJDZs) &
        DJDZs(DJDZsPerm(NodeIndexes(1:n)))=DJDZs(DJDZsPerm(NodeIndexes(1:n)))+Nodalzsb(1:n)
   IF (ComputeDJDZb) &
        DJDZb(DJDZbPerm(NodeIndexes(1:n)))=DJDZb(DJDZbPerm(NodeIndexes(1:n)))+Nodalzbb(1:n)
   IF (ComputeDJDRho) &
        DJDRho(DJDRhoPerm(NodeIndexes(1:n)))=DJDRho(DJDRhoPerm(NodeIndexes(1:n)))+Nodalrhob(1:n)
   IF (ComputeDJDEta) &
        DJDEta(DJDEtaPerm(NodeIndexes(1:n)))=DJDEta(DJDEtaPerm(NodeIndexes(1:n)))+Nodaletab(1:n)

  END DO
  
!  
! Neumann condition
!
  DO t=1,GetNOFBoundaryElements()
     BoundaryElement => GetBoundaryElement(t)


     IF ( .NOT. ActiveBoundaryElement(BoundaryElement,Solver) ) CYCLE
     IF ( GetElementFamily() == 1 ) CYCLE

     NodeIndexes => BoundaryElement % NodeIndexes
     !IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE


     n = GetElementNOFNodes(BoundaryElement)
     FORCE = 0.0e0
     STIFF = 0.0e0

 ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (STDOFs == 1) THEN
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (STDOFs == 2) THEN
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA with SSA var DOFs=',&
                STDOFs, '. Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

     BC => GetBC()
     IF (.NOT.ASSOCIATED( BC ) ) CYCLE

! Find the nodes for which 'Calving Front' = True             
     CalvingFront=.False. 
     CalvingFront = ListGetLogical( BC, 'Calving Front', GotIt )
     IF (CalvingFront) THEN
        NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))
        NodalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
     
       ! Need to access Parent Element to get Material properties
        other_body_id = BoundaryElement % BoundaryInfo % outbody
        IF (other_body_id < 1) THEN ! only one body in calculation
          ParentElement => BoundaryElement % BoundaryInfo % Right
          IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
        ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
          ParentElement => BoundaryElement %  BoundaryInfo % Right
          IF (ParentElement % BodyId == other_body_id) ParentElement =>  BoundaryElement % BoundaryInfo % Left
        END IF

        ! Read Density in the Material Section
        Material => GetMaterial(ParentElement)

        NodalDensity=0.0_dp
        NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,UnFoundFatal=.TRUE.)

        ! Read the gravity in the Body Force Section 
        BodyForce => GetBodyForce(ParentElement)
        NodalGravity = 0.0_dp
        IF ( ASSOCIATED( BodyForce ) ) THEN
           IF (STDOFs==1) THEN 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 2', n, NodeIndexes,UnFoundFatal=.TRUE.)
           ELSE 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 3', n, NodeIndexes,UnFoundFatal=.TRUE.)
           END IF
        END IF

        CALL LocalMatrixBCSSA(  STIFF, FORCE, BoundaryElement, n, ElementNodes,&
               NodalDensity, NodalGravity, NodalZb, NodalZs, rhow, sealevel , &
               nodalzsb,nodalzbb,nodalrhob)

        IF (ComputeDJDZs) &
             DJDZs(DJDZsPerm(NodeIndexes(1:n)))=DJDZs(DJDZsPerm(NodeIndexes(1:n)))+Nodalzsb(1:n)
        IF (ComputeDJDZb) &
             DJDZb(DJDZbPerm(NodeIndexes(1:n)))=DJDZb(DJDZbPerm(NodeIndexes(1:n)))+Nodalzbb(1:n)
        IF (ComputeDJDRho) &
             DJDRho(DJDRhoPerm(NodeIndexes(1:n)))=DJDRho(DJDRhoPerm(NodeIndexes(1:n)))+Nodalrhob(1:n)
     END IF
  END DO

  RETURN
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUVSSA(  STIFF, FORCE, Element, n, Nodes, gravity, &
           Density, Viscosity,NodalEtaDer,HaveEtaDer, LocalZb, LocalZs, LocalU, &
           LocalV, LocalBeta,NodalBetaDer,HaveBetaDer,iFriction,fm,LocalLinVelo, cm,&
           NodalGM,NodalBed,SEP,GLnIP,sealevel,rhow,&
           MinSRInv, STDOFs , &
            nodalbetab,nodalzsb,nodalzbb,nodalrhob,nodaletab )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
                     Viscosity(:), LocalZb(:), LocalZs(:), &
                     LocalU(:), LocalV(:) , LocalBeta(:),LocalLinVelo(:), &
                     nodalbetab(:),nodalzbb(:),nodalzsb(:),nodalrhob(:),nodaletab(:),&
                     NodalEtaDer(:),NodalBetaDer(:)
    LOGICAL :: HaveEtaDer,HaveBetaDer
    INTEGER :: n, cp , STDOFs
    INTEGER :: iFriction
    REAL(KIND=dp) :: cm,fm
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Newton
    REAL(KIND=dp) :: NodalGM(:),NodalBed(:)
    REAL(KIND=dp) :: sealevel,rhow
    LOGICAL :: SEP
    INTEGER :: GLnIP
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    REAL(KIND=dp) :: g, rho, eta, h, dhdx, dhdy , muder, bedrock,Hf
    REAL(KIND=dp) :: gradS(2),gradSb(2),Slip,  A(2,2),  Exx, Eyy, Exy, Ezz, Ee, MinSRInv                            
    REAL(KIND=dp) :: beta,Slip2,Velo(2),LinVelo,ub
    REAL(KIND=dp) :: betab,hb,rhob,etab
    REAL(KIND=dp) :: Id2
    LOGICAL :: Stat, NewtonLin
    INTEGER :: i, j, t, p, q 
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: PartlyGroundedElement

    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------


    nodalbetab = 0.0_dp
    nodalzsb = 0.0_dp
    nodalzbb = 0.0_dp
    nodalrhob = 0.0_dp
    nodaletab = 0.0_dp

    IF (SEP) THEN
     PartlyGroundedElement=(ANY(NodalGM(1:n).GE.0._dp).AND.ANY(NodalGM(1:n).LT.0._dp))
     IF (PartlyGroundedElement) THEN
        IP = GaussPoints( Element , np=GLnIP )
     ELSE
        IP = GaussPoints( Element )
     ENDIF
    ELSE
     IP = GaussPoints( Element )
    ENDIF

    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

! Needed Integration Point value

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )
       eta = SUM( Viscosity(1:n) * Basis(1:n) )
       gradS = 0._dp
       gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
       if (STDOFs == 2) gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
       h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n) )
       beta = SUM( LocalBeta(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
! In the non-linear case, effective viscosity       
       Id2=1.0_dp
       IF (cm.NE.1.0_dp) THEN
           Exx = SUM(LocalU(1:n)*dBasisdx(1:n,1))
           Eyy = 0.0
           Exy = 0.0
           IF (STDOFs.EQ.2) THEN
              Eyy = SUM(LocalV(1:n)*dBasisdx(1:n,2))
              Ezz = -Exx - Eyy
              Exy = SUM(LocalU(1:n)*dBasisdx(1:n,2))
              Exy = 0.5*(Exy + SUM(LocalV(1:n)*dBasisdx(1:n,1)))
              Ee = 0.5*(Exx**2.0 + Eyy**2.0 + Ezz**2.0) + Exy**2.0
              !Ee = SQRT(Ee)
           ELSE
              !Ee = ABS(Exx)
              Ee = Exx * Exx
           END IF
           muder = eta * 0.5 * (2**cm) * ((cm-1.0)/2.0) *  Ee**((cm-1.0)/2.0 - 1.0)
           IF (sqrt(Ee) < MinSRInv) then
                Ee = MinSRInv*MinSRInv
                muder = 0.0_dp
           Endif
           Id2 =  0.5 * (2**cm) * Ee**((cm-1.0)/2.0)
       END IF 

       betab = 0.0
       hb = 0.0
       rhob = 0.0
       etab = 0.0
       gradsb = 0.0

       A = 0.0_dp
       DO p=1,n
         DO q=1,n
         A(1,1) = 2.0*dBasisdx(q,1)*dBasisdx(p,1)  
           IF (STDOFs.EQ.2) THEN
           A(1,1) = A(1,1) + 0.5*dBasisdx(q,2)*dBasisdx(p,2)
           A(1,2) = dBasisdx(q,2)*dBasisdx(p,1) + &
                             0.5*dBasisdx(q,1)*dBasisdx(p,2)
           A(2,1) = dBasisdx(q,1)*dBasisdx(p,2) + &
                             0.5*dBasisdx(q,2)*dBasisdx(p,1)
           A(2,2) = 2.0*dBasisdx(q,2)*dBasisdx(p,2) +&
                             0.5*dBasisdx(q,1)*dBasisdx(p,1)  
           END IF

           DO i=1,STDOFs
             betab = betab +  Basis(q) * Basis(p) * IP % S(t) * detJ *&
                     (- VelocityN(STDOFs*(VeloNPerm(NodeIndexes(q))-1)+i) * &
                        VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i))

             DO j=1,STDOFs
                etab = etab + 2.0 * h * Id2 * A(i,j) * IP % S(t) * detJ *&
                       (- VelocityN(STDOFs*(VeloNPerm(NodeIndexes(q))-1)+j) * &
                       VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i))

                hb = hb + 2.0 * eta * Id2 * A(i,j)  * IP % S(t) * detJ *&
                     (- VelocityN(STDOFs*(VeloNPerm(NodeIndexes(q))-1)+j) * &
                        VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i))
             END DO  !j
           END DO !i

         END DO !q

         DO i=1,STDOFs
             rhob = rhob - h*g*gradS(i) * IP % s(t) * detJ * Basis(p) *&
                        VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i)
             hb = hb - rho*g*gradS(i) * IP % s(t) * detJ * Basis(p) *&
                        VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i)
             gradsb(i) = gradsb(i) - rho * g * h * IP % s(t) * detJ * Basis(p) *&
                        VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i)
         END DO !i
       END DO !p

       IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1 !linear
       IF (iFriction > 1) THEN
           LinVelo = SUM( LocalLinVelo(1:n) * Basis(1:n) )
           Velo = 0.0_dp
           Velo(1) = SUM(LocalU(1:n) * Basis(1:n))
           IF (STDOFs == 2) Velo(2) = SUM(LocalV(1:n) * Basis(1:n))
           ub = SQRT(Velo(1)*Velo(1)+Velo(2)*Velo(2))
           IF (ub < LinVelo) then 
              ub = LinVelo
           ENDIF
           IF (iFriction == 2) THEN !Weertman
             betab = betab * ub**(fm-1.0_dp)           
           ELSE IF (iFriction == 3) THEN !regularized coulomb
             betab = betab * ub**(fm-1.0_dp) / (ub + U0)**fm           
           END IF
       END IF

       IF (SEP) THEN
        IF (ALL(NodalGM(1:n).LT.0._dp)) THEN
           betab=0._dp
        ELSE IF (PartlyGroundedElement) THEN
           bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
           Hf= rhow * (sealevel-bedrock) / rho
           if (h.lt.Hf) betab=0._dp
        END IF
       END IF

       IF (HaveBetaDer) THEN
         nodalbetab(1:n)=nodalbetab(1:n)+betab*Basis(1:n)*NodalBetaDer(1:n)
       ELSE
         nodalbetab(1:n)=nodalbetab(1:n)+betab*Basis(1:n)
       ENDIF

       nodalzbb(1:n)=nodalzbb(1:n)-hb*Basis(1:n)
       nodalzsb(1:n)=nodalzsb(1:n)+hb*Basis(1:n)
       Do i=1,STDOFS
        nodalzsb(1:n)=nodalzsb(1:n)+gradsb(i)*dBasisdx(1:n,i)
       End do
       nodalrhob(1:n)=nodalrhob(1:n)+rhob*Basis(1:n)

       IF (HaveEtaDer) THEN
         nodaletab(1:n)=nodaletab(1:n)+etab*Basis(1:n)*NodalEtaDer(1:n)
       ELSE
         nodaletab(1:n)=nodaletab(1:n)+etab*Basis(1:n)
       ENDIF

    END DO !IP

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSA
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSSA(  STIFF, FORCE, Element, n, ENodes, Density, & 
                      Gravity, LocalZb, LocalZs, rhow, &
                      sealevel,nodalzsb,nodalzbb,nodalrhob)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) ::  ENodes
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:),  density(:), Gravity(:), LocalZb(:),&
                         LocalZs(:),rhow, sealevel,&
                         nodalzsb(:),nodalzbb(:),nodalrhob(:)
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), rhoi, g, alpha, h, h_im,norm
    REAL(KIND=dp) :: alphab,zsb,zbb,rhob
    LOGICAL :: Stat
    INTEGER :: t, i
    TYPE(GaussIntegrationPoints_t) :: IP

!------------------------------------------------------------------------------
    nodalzsb=0.0
    nodalzbb=0.0
    nodalrhob=0.0

! The front force is a concentrated nodal force in 1D-SSA and
! a force distributed along a line in 2D-SSA    

! 1D-SSA Case : concentrated force at each nodes
    IF (STDOFs==1) THEN  !1D SSA but should be 2D problem (does elmer work in 1D?)
      DO i = 1, n
         g = ABS( Gravity(i) )
         rhoi = Density(i)
         h = LocalZs(i)-LocalZb(i) 
         h_im=max(0._dp,sealevel-LocalZb(i))
         alpha=0.5 * g * (rhoi * h**2.0 - rhow * h_im**2.0)

         alphab=VelocityD(STDOFs*(VeloDPerm(NodeIndexes(i))-1)+1)
         nodalzsb(i)=+alphab*2.0*h*rhoi*g*0.5
         nodalzbb(i)=-alphab*2.0*h*rhoi*g*0.5
         if ((sealevel-LocalZb(i)).GT.0._dp) then
           nodalzbb(i)=nodalzbb(i)+alphab*2.0*h_im*rhow*g*0.5
         endif
         nodalrhob(i)=alphab * 0.5 * g *  h**2.0
      END DO

! 2D-SSA Case : force distributed along the line       
! This will work in DIM=3D only if working with Extruded Mesh and Preserve
! Baseline as been set to True to keep the 1D-BC 
    ELSE IF (STDOFs==2) THEN

          IP = GaussPoints( Element )
          DO t=1,IP % n
             stat = ElementInfo( Element, ENodes, IP % U(t), IP % V(t), &
                 IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
 
             g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
             rhoi = SUM( Density(1:n) * Basis(1:n) )
             h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n))
             h_im = max(0.0_dp , SUM( (sealevel-LocalZb(1:n)) * Basis(1:n)) )
             alpha=0.5 * g * (rhoi * h**2.0 - rhow * h_im**2.0)

! Normal in the (x,y) plane
             Normal = NormalVector( Element, ENodes, IP % U(t), IP % V(t), .TRUE.)
             norm=SQRT(normal(1)**2.0+normal(2)**2.0)
             Normal(1) = Normal(1)/norm
             Normal(2) = Normal(2)/norm

             rhob=0._dp
             zbb=0.0_dp
             zsb=0.0_dp
             DO p=1,n
                DO i=1,STDOFs
                   alphab=VelocityD(STDOFs*(VeloDPerm(NodeIndexes(p))-1)+i)*&
                          Normal(i) * IP % s(t) * detJ * Basis(p)
                   rhob=rhob+alphab * 0.5 * g *  h**2.0
                   zsb=zsb+alphab*2.0*h*rhoi*g*0.5
                   zbb=zbb-alphab*2.0*h*rhoi*g*0.5
                   if (SUM( (sealevel-LocalZb(1:n)) * Basis(1:n)).GT.0.0_dp) Then
                      zbb=zbb+alphab*2.0*h_im*rhow*g*0.5
                   endif
                END DO !p
             END DO !q

             nodalrhob(1:n)=nodalrhob(1:n)+rhob*Basis(1:n)
             nodalzsb(1:n)=nodalzsb(1:n)+zsb*Basis(1:n)
             nodalzbb(1:n)=nodalzbb(1:n)+zbb*Basis(1:n)
          END DO !IP

    ELSE   

      CALL FATAL('SSASolver-SSABasalSolver','Do not work for STDOFs <> 1 or 2')

    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCSSA

END SUBROUTINE AdjointSSA_GradientSolver