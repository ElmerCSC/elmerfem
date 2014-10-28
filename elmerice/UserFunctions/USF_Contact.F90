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
FUNCTION SlidCoef_Contact ( Model, nodenumber, y) RESULT(Bdrag)

  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  TYPE(variable_t), POINTER :: TimeVar, NormalVar, VarSurfResidual, GroundedMaskVar, HydroVar, FlowVariable, DistanceVar
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: Element, CurElement, BoundaryElement
  TYPE(Nodes_t), SAVE :: Nodes

  REAL(KIND=dp), POINTER :: NormalValues(:), ResidValues(:), GroundedMask(:), Hydro(:), FlowValues(:), Distance(:)
  REAL(KIND=dp) :: Bdrag, t, told, test(3), thresh
  REAL(KIND=dp), ALLOCATABLE :: Normal(:), Fwater(:), Fbase(:), NormalTest(:)

  INTEGER, POINTER :: NormalPerm(:), ResidPerm(:), GroundedMaskPerm(:), HydroPerm(:), FlowPerm(:), DistancePerm(:)
  INTEGER :: nodenumber, ii, DIM, GL_retreat, n, tt, Nn, jj, nnn, MSum, ZSum

  LOGICAL :: FirstTime = .TRUE., FirstTimeTime = .TRUE., GotIt, Yeschange, GLmoves

  REAL (KIND=dp) ::  y, relChange, relChangeOld, Sliding_Budd, Sliding_Weertman, Friction_Coulomb, C, m

  REAL(KIND=dp) :: NonLinIter, comp, cond
  CHARACTER*20 :: USF_Name='SlidCoef_Contact', Sl_Law

  SAVE FirstTime, yeschange, told, GLmoves, thresh
  SAVE DIM, USF_Name, Normal, Fwater, Fbase, relChangeOld, Sl_Law

!----------------------------------------------------------------------------

! Real time import
  Timevar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)

! GroundedMask import
  GroundedMaskVar => VariableGet( Model % Mesh % Variables, 'GroundedMask')
  IF ( ASSOCIATED( GroundedMaskVar ) ) THEN
    GroundedMask => GroundedMaskVar % Values
    GroundedMaskPerm => GroundedMaskVar % Perm
  ELSE
    CALL FATAL( USF_Name, 'need to get GroundedMask')
  END IF

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

    ! Possiblity to fix the grounding line, default is a moving Grounding Line
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
  IF (thresh.GT.0) THEN
    DistanceVar => VariableGet( Model % Mesh % Variables, 'Distance')
    IF ( ASSOCIATED( DistanceVar ) ) THEN
      Distance => DistanceVar % Values
      DistancePerm => DistanceVar % Perm
    ELSE
      CALL FATAL( USF_Name, 'need to get DistanceSolver for the use of "non detachment inland distance"' )
    END IF
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look at the convergence of the FlowSolver.
  ! If relative change < 1.e-3, test if traction occurs. To apply one time
  !
  ! Only to release contact between bed and ice as hydrostatic pressure is higher than normal stress
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Normal = 0.0_dp
  Fwater = 0.0_dp
  Fbase = 0.0_dp

  IF ( (relChange.NE.relChangeOld) .AND. (relchange.GT.0.0_dp) .AND. (relchange.LT.1.0e-3_dp) .AND. (yesChange) .AND. GLmoves ) THEN
    ! Change the basal condition just once per timestep
    yesChange = .FALSE.

    CALL Info(USF_name,'FLOW SOLVER HAS SLIGHTLY CONVERGED: look for new basal conditions', Level=3)

    VarSurfResidual => VariableGet( Model % Mesh % Variables, 'Flow Solution Loads' )
    IF ( ASSOCIATED( VarSurfResidual ) ) THEN
      ResidPerm => VarSurfResidual  % Perm
      ResidValues => VarSurfResidual % Values
    ELSE
      WRITE(Message, '(A)') '> Flow Solution Loads< not associated'
      CALL FATAL( USF_Name, Message)
    END IF

    NormalVar => VariableGet(Model % Variables,'Normal Vector')
    IF ( ASSOCIATED( NormalVar ) ) THEN
      NormalPerm => NormalVar % Perm
      NormalValues => NormalVar % Values
    ELSE
      WRITE(Message, '(A)') '>Normal Vector< not associated'
      CALL FATAL(USF_Name, Message)
    END IF
    
    !Force exerted by the water, computed for each good boundary nodes (whatever on the bed or floating)
    !From GetHydrostaticLoads

    HydroVar => VariableGet( Model % Mesh % Variables, 'Fw')
    IF ( ASSOCIATED( HydroVar ) ) THEN
      Hydro => HydroVar % Values
      HydroPerm => HydroVar % Perm
    ELSE
      WRITE(Message, '(A)') '>Fw< not associated'
      CALL FATAL( USF_Name, Message)
    END IF

    ! Retreat of the Grounding line if Hydro loads higher than residual values
    GL_retreat = 0

    CurElement => Model % CurrentElement
    DO tt = 1, Model % NumberOfBoundaryElements

      Element => GetBoundaryElement(tt)
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      n = GetElementNOFNodes(Element)

      CALL GetElementNodes(Nodes, Element)

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
        n = GetElementNOFNodes(Element)

        CALL GetElementNodes(Nodes, Element)
        MSum = 0
        ZSum = 0

        DO ii = 1,n

          Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
          ! the grounded mask is not defined here
          IF (Nn==0) CYCLE
          MSum = MSum + GroundedMask(Nn)
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
  ! END relChange.NE.relChangeOld .AND. relchange.GT.0.0_dp .AND. relchange.LT.1.0e-3_dp .AND. yesChange 

  relChangeOld = relChange  

  !!!!!!!!!!!!!!!!!!!!!!!!
  ! Weertman sliding law
  !!!!!!!!!!!!!!!!!!!!!!!!

  ! Bdrag = 0 if GroundedMask=-1
  ! Bdrag = Sliding_law if GroundedMask=0 or 1

  IF (GroundedMaskPerm(nodenumber) > 0) THEN
  ! for the bottom surface, where the GroundedMask is defined

    cond = GroundedMask(GroundedMaskPerm(nodenumber))

    IF (cond >= -0.5_dp) THEN
       ! grounded node
       SELECT CASE(Sl_law)
       CASE ('weertman')
          Bdrag = Sliding_weertman(Model, nodenumber, y)
       CASE ('budd')
          Bdrag = Sliding_Budd(Model, nodenumber, y)
       CASE ('coulomb')
          Bdrag = Friction_Coulomb(Model, nodenumber, y)
       CASE DEFAULT
          WRITE(Message, '(A,A)') 'Sliding law not recognised ',Sl_law
          CALL FATAL( USF_Name, Message)
       END SELECT
!      IF (Sl_Law == 'weertman') Bdrag = Sliding_weertman(Model, nodenumber, y)
!      IF (Sl_Law == 'budd')     Bdrag = Sliding_Budd(Model, nodenumber, y)
!      IF (Sl_Law == 'coulomb')  Bdrag = Friction_Coulomb(Model, nodenumber, y)
    ELSE
    ! floating node
      Bdrag = 0.0_dp
    END IF

  ELSE
  ! for other surfaces, typically for lateral surfaces within buttressing experiments
    Bdrag = Sliding_weertman(Model, nodenumber, y)

  END IF

END FUNCTION SlidCoef_Contact



