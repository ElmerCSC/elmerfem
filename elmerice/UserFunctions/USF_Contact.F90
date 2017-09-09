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
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *  Date Modifications:2007/10/25. Gael Durand
! * 
! *****************************************************************************
!>  tests if resting ice becomes floating ice (GroundedMask from 0 or 1 to -1)
!>  
!>  Return a friction coefficient
!>  -from sliding weertman for resting ice (GroundedMask = 0 or 1)
!>  -of 0 for floating ice (GroundedMask = -1)
!>  2014 : Introduce 3 different ways of defining the grounding line (Mask = 0)
!>  Last Grounded ; First Floating ; Discontinuous

FUNCTION SlidCoef_Contact ( Model, nodenumber, y) RESULT(Bdrag)

  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  TYPE(variable_t), POINTER :: TimeVar, NormalVar, VarSurfResidual, GroundedMaskVar, HydroVar, DistanceVar, FrictionVar
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: Element, CurElement, BoundaryElement
  TYPE(Nodes_t), SAVE :: Nodes

  REAL(KIND=dp), POINTER :: NormalValues(:), ResidValues(:), GroundedMask(:), Hydro(:), &
       Distance(:), FrictionValues(:)
  REAL(KIND=dp) :: Bdrag, t, told, thresh
  REAL(KIND=dp), ALLOCATABLE :: Normal(:), Fwater(:), Fbase(:)

  INTEGER, POINTER :: NormalPerm(:), ResidPerm(:), GroundedMaskPerm(:), HydroPerm(:), &
       DistancePerm(:), FrictionPerm(:)
  INTEGER :: nodenumber, ii, DIM, GL_retreat, n, tt, Nn, jj, MSum, ZSum

  LOGICAL :: FirstTime = .TRUE., GotIt, Yeschange, GLmoves, Friction,UnFoundFatal=.TRUE.

  REAL (KIND=dp) ::  y, relChange, relChangeOld, Sliding_Budd, Sliding_Weertman, Friction_Coulomb

  REAL(KIND=dp) :: comp, cond, TestContact
  CHARACTER(LEN=MAX_NAME_LEN) :: USF_Name='SlidCoef_Contact', Sl_Law, GLtype, FrictionVarName

  SAVE FirstTime, yeschange, told, GLmoves, thresh, GLtype, TestContact
  SAVE DIM, USF_Name, Normal, Fwater, Fbase, relChangeOld, Sl_Law
  SAVE FrictionVar, FrictionValues, FrictionPerm, BC

!----------------------------------------------------------------------------

! Real time import
  Timevar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)

! GroundedMask import
  GroundedMaskVar => VariableGet( Model % Mesh % Variables, 'GroundedMask',UnFoundFatal=UnFoundFatal)
  GroundedMask => GroundedMaskVar % Values
  GroundedMaskPerm => GroundedMaskVar % Perm
  
  relchange = Model % Solver % Variable % NonLinChange

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First time step for the First time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (FirstTime) THEN
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     n = Model % MaxElementNodes
     told = t
     
! means have the possibility to change
     yesChange = .TRUE.
     ALLOCATE( Normal(DIM), Fwater(DIM), Fbase(DIM) )
     
     relChangeOld = relChange
    
     ! choice of the Sliding Law
     BoundaryElement => Model % CurrentElement
     BC => GetBC(BoundaryElement)
     
     Sl_Law = GetString( BC, 'Sliding Law', GotIt )
     IF (.NOT.Gotit) THEN
        CALL FATAL(USF_Name,'No "Sliding law" Name given')
     END IF

     GLtype = GetString( BC, 'Grounding Line Definition', GotIt )
     IF (.NOT.Gotit) THEN
        GLtype = 'last grounded'
        CALL Info(USF_Name, 'Grounded Line Defined as the last Grounded point', Level=3)
     ELSE
        WRITE(Message, '(A,A)') 'Grounding Line Defined as ', GLtype
        CALL Info(USF_Name, Message, Level=3)
     END IF
     
     ! Possibility to fix the grounding line, default is a moving Grounding Line
     GLmoves = GetLogical( BC, 'Grounding line moves', GotIt )
     IF (.NOT.GotIt) THEN
        GLmoves = .TRUE.
     END IF
     IF (GLmoves) THEN
        CALL Info(USF_Name, 'GL may move by default', Level=3)
        CALL Info(USF_Name, 'If you want to fix the Grounding Line, put the keyword "Grounding line moves" to False', Level=3)
     ELSE
        CALL Info(USF_Name, 'GL will be fixed', Level=3)
     END IF
     
     TestContact = GetConstReal( BC, 'Test Contact Tolerance', GotIt )
     IF (.NOT.Gotit) THEN
        TestContact = 1.0e-3     
        CALL Info(USF_Name, 'Contact will be tested for a tolerance of 1.0e-3', Level=3)
     ELSE
        WRITE(Message, '(A,e14.8)') 'Contact tested for a tolerance of ', TestContact
        CALL Info(USF_Name, Message, Level=3)
     END IF
     
     ! Possibility to avoid detachement from nodes that are too far inland from the Grounding line
     ! Uses the DistanceSolver
     ! Default is non possible detachment
     thresh = GetConstReal( BC, 'non detachment inland distance', GotIt )
     IF (.NOT.GotIt) THEN
        thresh = -10000.0_dp
        CALL INFO( USF_Name, 'far inland nodes have the possibility to detach by default', Level=3)
        CALL INFO( USF_Name, 'to avoid detachment (when bedrock is well below sea level),', Level=3)
        CALL INFO( USF_Name, 'use the keyword "non detachment inland distance" to the distance you wish', Level=3)
        CALL INFO( USF_Name, 'This works with the DistanceSolver', Level=3)
     ELSE
        CALL INFO( USF_Name, 'far inland nodes will not detach', level=3)
     END IF
  ENDIF
  
  IF(Sl_Law(1:10) == 'prescribed') THEN
     FrictionVarName = GetString( BC, 'Friction Variable Name', GotIt )
     IF(.NOT. GotIt) CALL Fatal(USF_Name, 'Prescribed friction requested but no &
          "Friction Variable Name" found!')
     FrictionVar => VariableGet(Model % Mesh % Variables, FrictionVarName)
     IF(ASSOCIATED(FrictionVar)) THEN
        FrictionValues => FrictionVar % Values
        FrictionPerm => FrictionVar % Perm
     ELSE
        WRITE(Message, '(A,A)') 'Unable to find variable: ',FrictionVarName
        CALL Fatal(USF_Name, Message)
     END IF
  END IF


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First time step for a New time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF ( t > told ) THEN
     told = t
     yesChange = .TRUE.
     relChangeOld = relChange
  END IF
  
  ! to use the non detachment possibility when a grounded node is too far from the grounding line
  ! and positioned on a well below sea level bedrock
  IF (thresh.GT.0.0) THEN
     DistanceVar => VariableGet( Model % Mesh % Variables, 'Distance',UnFoundFatal=UnFoundFatal)
     Distance => DistanceVar % Values
     DistancePerm => DistanceVar % Perm
  END IF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look at the convergence of the FlowSolver.
  ! If relative change < TestContact, test if traction occurs. To apply one time
  !
  ! Only to release contact between bed and ice as hydrostatic pressure is higher than normal stress
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Normal = 0.0_dp
  Fwater = 0.0_dp
  Fbase = 0.0_dp

  IF ( (relChange.NE.relChangeOld) .AND. (relchange.GT.0.0_dp) .AND. & 
       &            (relchange.LT.TestContact) .AND. (yesChange) .AND. GLmoves ) THEN
     ! Change the basal condition just once per timestep
     yesChange = .FALSE.

     CALL Info(USF_name,'FLOW SOLVER HAS SLIGHTLY CONVERGED: look for new basal conditions', Level=3)

     VarSurfResidual => VariableGet( Model % Mesh % Variables, 'Flow Solution Loads',UnFoundFatal=UnFoundFatal)
     ResidPerm => VarSurfResidual  % Perm
     ResidValues => VarSurfResidual % Values

     NormalVar => VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
     
     !Force exerted by the water, computed for each good boundary nodes (whatever on the bed or floating)
     !From GetHydrostaticLoads
     
     HydroVar => VariableGet( Model % Mesh % Variables, 'Fw',UnFoundFatal=UnFoundFatal)
     Hydro => HydroVar % Values
     HydroPerm => HydroVar % Perm
     
     ! Retreat of the Grounding line if Hydro loads higher than residual values
     GL_retreat = 0

     CurElement => Model % CurrentElement
     DO tt = 1, Model % NumberOfBoundaryElements

        Element => GetBoundaryElement(tt)
        IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
        n = GetElementNOFNodes(Element)

        CALL GetElementNodes(Nodes, Element)

        IF (ANY(GroundedMaskPerm(Element % NodeIndexes(1:n))==0)) CYCLE
        DO ii = 1,n

           Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
           ! the grounded mask is not defined here
           IF (Nn==0) CYCLE
           IF (GroundedMask(Nn) < -0.5_dp) CYCLE
           
           jj = Element % NodeIndexes(ii)
           
           ! comparison between water load and reaction
           
           Normal = NormalValues(DIM*(NormalPerm(jj)-1)+1 : DIM*NormalPerm(jj))
           Fwater = Hydro(DIM*(HydroPerm(jj)-1)+1 : DIM*HydroPerm(jj))
           Fbase = ResidValues((DIM+1)*(ResidPerm(jj)-1)+1 : (DIM+1)*ResidPerm(jj)-1)

           ! comparison between water pressure and bed action
           comp = ABS( SUM( Fwater * Normal ) ) - ABS( SUM( Fbase * Normal ) )
           
           
           IF (comp .GE. 0.0_dp) THEN
              IF (thresh.LE.0.0_dp) THEN
                 GroundedMask(Nn) = -1.0_dp
                 GL_retreat = GL_retreat + 1
                 PRINT *, 'Retreat of the Grounding Line : '
                 PRINT *, Nodes % x(ii), Nodes % y(ii), Nodes % z(ii)
              ELSE
                 IF ( Distance(DistancePerm(Element % NodeIndexes(ii))).LE.thresh ) THEN
                    GroundedMask(Nn) = -1.0_dp
                    GL_retreat = GL_retreat + 1
                    PRINT *, 'Retreat of the Grounding Line : '
                    PRINT *, Nodes % x(ii), Nodes % y(ii), Nodes % z(ii)
                 END IF
              END IF
           END IF
        END DO
        
     END DO
     Model % CurrentElement => CurElement
     
     ! with the previous step
     ! Some 0 (Grounding line) may have been replaced by -1 (floating nodes)
     ! here replacement of some 1 by 0's
     
     IF (GL_retreat.GT.0) THEN
        CurElement => Model % CurrentElement
        DO tt = 1, Model % NumberOfBoundaryElements
           
           Element => GetBoundaryElement(tt)
           IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
           n = GetElementNOFNodes(Element)
           
           CALL GetElementNodes(Nodes, Element)
           MSum = 0
           ZSum = 0
           
           IF (ANY(GroundedMaskPerm(Element % NodeIndexes(1:n))==0)) CYCLE
           DO ii = 1,n
              
              Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
              ! the grounded mask is not defined here
              IF (Nn==0) CYCLE
              MSum = MSum + INT(GroundedMask(Nn))
              IF (GroundedMask(Nn)==0.0_dp) ZSum = ZSum + 1
              
           END DO
           
           IF (MSum+ZSum .LT. n) THEN
              DO ii=1,n
                 Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
                 IF (Nn==0) CYCLE
                 
                 IF (GroundedMask(Nn)==1.0_dp) THEN
                    GroundedMask(Nn)=0.0_dp
                 END IF

              END DO
           END IF
           
        END DO
        Model % CurrentElement => CurElement
     END IF
  END IF
  
  relChangeOld = relChange  
  
  
  IF (GroundedMaskPerm(nodenumber) > 0) THEN
  ! for the bottom surface, where the GroundedMask is defined
     cond = GroundedMask(GroundedMaskPerm(nodenumber))
     
  ! Definition of the Grounding line in term of friction 
  ! If GLtype = Last Grounded -> Bdrag = 0 if Nodal Mask < 0
  ! If GLtype = First Floating -> Bdrag = 0 if Nodal Mask <=0
  ! If GLtype = Discontinuous -> Bdrag = 0 if at one node of the element Mask<=0

     Friction = .FALSE.
     SELECT CASE(GLtype)
     CASE('last grounded')
        IF (cond > -0.5) Friction = .TRUE.  
     CASE('first floating')
        IF (cond > 0.5) Friction = .TRUE. 
     CASE('discontinuous')
        BoundaryElement => Model % CurrentElement
        IF (ALL(GroundedMask(GroundedMaskPerm(BoundaryElement % NodeIndexes))>-0.5)) Friction = .TRUE. 
     CASE DEFAULT
        WRITE(Message, '(A,A)') 'GL type not recognised ', GLtype 
        CALL FATAL( USF_Name, Message)
     END SELECT

     IF (Friction) THEN
        ! grounded node
        SELECT CASE(Sl_law)
        CASE ('weertman')
           Bdrag = Sliding_weertman(Model, nodenumber, y)
        CASE ('budd')
           Bdrag = Sliding_Budd(Model, nodenumber, y)
        CASE ('coulomb')
           Bdrag = Friction_Coulomb(Model, nodenumber, y)
        CASE('prescribed beta2')
           Bdrag = FrictionValues(FrictionPerm(nodenumber))**2.0_dp
        CASE('prescribed power')
           Bdrag = 10.0_dp**FrictionValues(FrictionPerm(nodenumber))
        CASE DEFAULT
           WRITE(Message, '(A,A)') 'Sliding law not recognised ',Sl_law
           CALL FATAL( USF_Name, Message)
        END SELECT
     ELSE
        ! floating node
        Bdrag = 0.0_dp
     END IF
  ELSE
     ! for other surfaces, typically for lateral surfaces within buttressing experiments
     Bdrag = Sliding_weertman(Model, nodenumber, y)
  END IF
END FUNCTION SlidCoef_Contact



