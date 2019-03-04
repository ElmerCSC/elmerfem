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
! *  Module for simulating electro-kinetic flows with heat transfer
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Thomas Zwinger
! *  Email:   Juha.Ruokolainen@csc.fi, Thomas.Zwinger@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 13 Apr 2005
! *
! *****************************************************************************/
!
!******************************************************************************
! *        file for routines used for electro-kinetic flow simulation including
! *        heat transfer.
! *        file contains following functions:
! *
! *        helmholtz_smoluchowski1/2/3: Elmer Interface for component-whise
! *                                     computation of HS-velocity to function
! *                                     helmholtz_smoluchowski_comp
! *
! *        helmholtz_smoluchowski_comp:  calculates components of the velocity
! *                                     given by the electroosmotic slip
! *                                     condition at the boundary;
! *
! *        helmholtz_smoluchowski: calculates absolute tangential value of the
! *                                Helmholtz-Smoluchowski velocity for
! *                                electro-osmotic velocity slip condition
! *
! *        getJouleHeat: computes inductive heat source as a function of 
! *                       electric field. Needs conductivity as material 
! *                       parameter input
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
!> Computes Helmholtz Smoluchowski velocity in x-direction.
!> \ingroup UDF
!----------------------------------------------------------------------------------
FUNCTION helmholtz_smoluchowski1( Model, NodeNumber, dummyargument) RESULT(hs_velocity1)
  USE DefUtils
  IMPLICIT NONE
! external variables
! ------------------
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) :: dummyargument, hs_velocity1

  INTERFACE 
     FUNCTION helmholtz_smoluchowski_comp( Model, NodeNumber, direction) RESULT(hs_velocity_comp)
       USE DefUtils
       IMPLICIT NONE
       ! external variables
       ! ------------------
       TYPE(Model_t) :: Model
       INTEGER :: NodeNumber, direction
       REAL(KIND=dp) :: hs_velocity_comp
     END FUNCTION helmholtz_smoluchowski_comp
  END INTERFACE

  ! ---------------------
  ! HS-velocity component
  ! ---------------------
  hs_velocity1 = helmholtz_smoluchowski_comp( Model, NodeNumber, 1)
END FUNCTION helmholtz_smoluchowski1


!----------------------------------------------------------------------------------
!> Computes Helmholtz Smoluchowski velocity in y-direction.
!> \ingroup UDF
!----------------------------------------------------------------------------------
FUNCTION helmholtz_smoluchowski2( Model, NodeNumber, dummyargument) RESULT(hs_velocity2)
  USE DefUtils
  IMPLICIT NONE
! external variables
! ------------------
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) :: dummyargument, hs_velocity2

  INTERFACE 
     FUNCTION helmholtz_smoluchowski_comp( Model, NodeNumber, direction) RESULT(hs_velocity_comp)
       USE DefUtils
       IMPLICIT NONE
       ! external variables
       ! ------------------
       TYPE(Model_t) :: Model
       INTEGER :: NodeNumber, direction
       REAL(KIND=dp) :: hs_velocity_comp
     END FUNCTION helmholtz_smoluchowski_comp
  END INTERFACE

  ! ---------------------
  ! HS-velocity component
  ! ---------------------
  hs_velocity2 = helmholtz_smoluchowski_comp( Model, NodeNumber, 2)
END FUNCTION helmholtz_smoluchowski2


!----------------------------------------------------------------------------------
!> Computes Helmholtz Smoluchowski velocity in z-direction.
!> \ingroup UDF
!----------------------------------------------------------------------------------
FUNCTION helmholtz_smoluchowski3( Model, NodeNumber, dummyargument) RESULT(hs_velocity3)

  USE DefUtils
  IMPLICIT NONE
! external variables
! ------------------
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) :: dummyargument, hs_velocity3

  INTERFACE 
     FUNCTION helmholtz_smoluchowski_comp( Model, NodeNumber, direction) RESULT(hs_velocity_comp)
       USE DefUtils
       IMPLICIT NONE
       ! external variables
       ! ------------------
       TYPE(Model_t) :: Model
       INTEGER :: NodeNumber, direction
       REAL(KIND=dp) :: hs_velocity_comp
     END FUNCTION helmholtz_smoluchowski_comp
  END INTERFACE

  ! ---------------------
  ! HS-velocity component
  ! ---------------------
  hs_velocity3 = helmholtz_smoluchowski_comp( Model, NodeNumber, 3)
END FUNCTION helmholtz_smoluchowski3



!----------------------------------------------------------------------------------
!> Computes Helmholtz Smoluchowski velocity in the tangential direction.
!> \ingroup UDF
!----------------------------------------------------------------------------------
FUNCTION helmholtz_smoluchowski( Model, NodeNumber, dummyargument) RESULT(hs_velocity_tang)
  USE DefUtils
  IMPLICIT NONE
! external variables
! ------------------
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) :: dummyargument, hs_velocity_tang

  INTERFACE 
     FUNCTION helmholtz_smoluchowski_comp( Model, NodeNumber, direction) RESULT(hs_velocity_comp)
       USE DefUtils
       IMPLICIT NONE
       ! external variables
       ! ------------------
       TYPE(Model_t) :: Model
       INTEGER :: NodeNumber, direction
       REAL(KIND=dp) :: hs_velocity_comp
     END FUNCTION helmholtz_smoluchowski_comp
  END INTERFACE

  ! ---------------------
  ! HS-velocity component
  ! ---------------------
  hs_velocity_tang = helmholtz_smoluchowski_comp( Model, NodeNumber, 0)
END FUNCTION helmholtz_smoluchowski



!----------------------------------------------------------------------------------
!> Computes component of Helmholtz Smoluchowski velocity for given direction. Not usually called alone.
!----------------------------------------------------------------------------------
FUNCTION helmholtz_smoluchowski_comp( Model, NodeNumber, direction) RESULT(hs_velocity_comp)
  USE DefUtils
  IMPLICIT NONE
! external variables
! ------------------
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber, direction
  REAL(KIND=dp) :: hs_velocity_comp
! internal variables
! ------------------
  INTEGER :: i, N, DIM, istat, body_id, material_id, eq_id, other_body_id,&
       bc_id, Nboundary, BoundaryElementNode, coordinateDirection, &
       dummyIntegerArray(1)
  INTEGER, TARGET :: NodeNumberTarget(1)
  INTEGER, DIMENSION(:), POINTER :: NodeNumberArray
  REAL(KIND=dp) :: dummyArray(1), viscosity, dielectricConstant, vacuumPerm, zetapotential,&
       electricField(3), tang_electricField(3), norm_electricField,&
       U, V, W, Normal(3), Tangent(3), Tangent2(3), hs_velocity(3), eoMobility
  CHARACTER(LEN=MAX_NAME_LEN) :: ElectricFieldMethod, ValueName
  LOGICAL :: GotIt, ElectricFieldExists, FirstTime=.TRUE., SkipThisTime=.FALSE.,& 
       CalcMobility
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement  
  TYPE(Nodes_t) :: Nodes
  TYPE(Variable_t), POINTER :: EFieldComp
! remember this
! -------------
  SAVE Nodes, DIM, FirstTime, vacuumPerm, SkipThisTime

! provide arrays to let the ListGetReal/Integer function 
! work with just a nodal value to be read in
! -------------------------------------------------------
  NodeNumberTarget(1) = NodeNumber
  NodeNumberArray => NodeNumberTarget
!---------------------End of Header section ------------------------------------
  !-----------------------------------------------------------------
  !    return zero in simulation setup
  !-----------------------------------------------------------------
  hs_velocity_comp = 0.0d00
  eoMobility = 0.0d00
  SkipThisTime = ListGetLogical(Model % Simulation, 'Initializaton Phase', GotIt)  
  IF (SkipThisTime) RETURN

  !-----------------------------------------------------------------
  ! things to be done for the first time only
  !-----------------------------------------------------------------
  IF (FirstTime) THEN
     DIM = CoordinateSystemDimension()

     IF( DIM  == 1 ) THEN
       CALL Fatal('electrokinetics (helmholtz_smoluchowski_comp)',&
           'The electrokinetic BCs do not make sense in 1D')
     END IF

     IF( direction == 0 .AND. Dim > 2) THEN
       CALL Fatal('electrokinetics (helmholtz_smoluchowski_comp)',&
           'The tangential velocity may only be defined in 2D')
     END IF

     N = Model % MaxElementNodes 
     ALLOCATE(Nodes % x(N), Nodes % y(N), Nodes % z(N),&
          STAT = istat)
     IF (istat /= 0) THEN
        CALL FATAL('electrokinetics (helmholtz_smoluchowski_comp)','Allocations failed')
     END IF
     vacuumPerm = ListGetConstReal( Model % Constants, 'Permittivity of Vacuum', GotIt )
     IF (.NOT. GotIt) THEN
        CALL WARN('electrokinetics (helmholtz_smoluchowski_comp)',&
             'No value for >Permittivity of Vacuum< found in section Constants')
        CALL WARN('electrokinetics (helmholtz_smoluchowski_comp)',&
             '         Using default SI value 8.8542E-12')
        vacuumPerm = 8.8542D-12
     END IF
     FirstTime = .FALSE.
  END IF

  !-------------------------------------------------------------------------
  ! get some information upon active boundary element and its parent element
  !-------------------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL('electrokinetics (helmholtz_smoluchowski_comp)','No boundary element found')
  END IF
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  ! just to be on the save side, check again
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     CALL FATAL('electrokinetics (helmholtz_smoluchowski_comp)','No parent element found for boudnary element')
  END IF

  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF (.NOT. ASSOCIATED(ParentMaterial)) THEN
     CALL FATAL('electrokinetics (helmholtz_smoluchowski_comp)','No material values could be found')
  END IF
  eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation', &
       minv=1,maxv=Model % NumberOfEquations )

  bc_id = GetBCId(BoundaryElement)
  BC => GetBC(BoundaryElement)
  !-------------------------------------------
  ! Get normal of the boundary element at node
  !-------------------------------------------
  Nboundary = BoundaryElement % TYPE % NumberOfNodes
  DO BoundaryElementNode=1,Nboundary
     IF ( NodeNumber == BoundaryElement % NodeIndexes(BoundaryElementNode) ) EXIT
  END DO
  U = BoundaryElement % TYPE % NodeU(BoundaryElementNode)
  V = BoundaryElement % TYPE % NodeV(BoundaryElementNode)
  Nodes % x(1:Nboundary) = Model % Nodes % x(BoundaryElement % NodeIndexes)
  Nodes % y(1:Nboundary) = Model % Nodes % y(BoundaryElement % NodeIndexes)
  Nodes % z(1:Nboundary) = Model % Nodes % z(BoundaryElement % NodeIndexes)
  Normal = NormalVector( BoundaryElement, Nodes, U, V,.TRUE. )

  !-----------------------------------------------------------------
  ! get material parameters from boundary
  !-----------------------------------------------------------------  
  dummyArray = ListGetReal( BC, 'Zeta Potential', 1, NodeNumberArray, CalcMobility )

  IF( CalcMobility ) THEN
    zetapotential = dummyArray(1)

    ! get material parameters from bulk 
    !-----------------------------------------------------------------
    Model % CurrentElement => ParentElement ! we need this in case of the material parameter being a function
    dummyArray  = ListGetReal( ParentMaterial, 'Viscosity', 1, NodeNumberArray, GotIt)
    IF (.NOT. GotIt) THEN
      WRITE(Message,'(a,i0)' )'No viscosity found in material section no. ', material_id
      CALL FATAL( 'electrokinetics (helmholtz_smoluchowski)',Message )
    ELSE
      viscosity = dummyArray(1)
    END IF
    dummyArray= ListGetReal( ParentMaterial, 'Relative Permittivity', 1, NodeNumberArray, GotIt )
    IF (.NOT. GotIt) THEN
      WRITE(Message,'(a,i0)' )&
          'No keyword >Relative Permittivity< found in material section no. ', material_id
      CALL Info('electrokinetics (helmholtz_smoluchowski)',Message,Level=4)
      CALL Info('electrokinetics (helmholtz_smoluchowski)','setting to 1', Level=4)
      dielectricConstant = 1.0d00
    ELSE
      dielectricConstant = dummyArray(1)
    END IF

    ! calculate mobility
    !-----------------------------------------------------------------
    eoMobility = (zetapotential * dielectricConstant * vacuumPerm)/ viscosity
  ELSE
    dummyArray = ListGetReal( BC, 'EO Mobility', 1, NodeNumberArray, GotIt )
    IF (.NOT.GotIt) THEN
      WRITE(Message,'(a,i0)' )'Neither zeta potential nor EO mobility found in boundary condition no. ', bc_id
      CALL Info('electrokinetics (helmholtz_smoluchowski)',Message,Level=4)
      CALL Info('electrokinetics (helmholtz_smoluchowski)','setting EO mobility to 0', Level=4)
      eoMobility = 0.0d00
    ELSE
      eoMobility = dummyArray(1)
    END IF
  END IF

  !-----------------------------------------------------------------
  ! get electric field
  !-----------------------------------------------------------------
  ElectricFieldMethod = ListGetString( Model % Equations(eq_id) % &
         Values, 'Electric Field', GotIt )
  IF ( .NOT.GotIt ) THEN 
    WRITE(Message,'(a,i0)' )'No external electric field defined for Equation no. ', eq_id
    CALL Info('electrokinetics (helmholtz_smoluchowski)',Message, level=4)
    ElectricFieldExists = .FALSE.
    hs_velocity = 0.0D00
  ELSE
    ElectricFieldExists = .TRUE.
    IF ( ElectricFieldMethod == 'constant') THEN ! read field components from Material section
      ElectricFieldExists = .FALSE.
      ElectricField = 0.0_dp
      
      DO i=1,dim
        WRITE(ValueName,'(a,i2)' ) 'Electric Field',i
        dummyArray =  ListGetReal( ParentMaterial,ValueName, 1, NodeNumberArray, GotIt )
        IF( GotIt ) THEN
          electricField(i) = dummyArray(1) 
          ElectricFieldExists = .TRUE.
        END IF
      END DO
      
      IF (.NOT. ElectricFieldExists) THEN
        WRITE(Message,'(a,i0)' )'No component for >Electric Field {1,2,3}< found in Material ',&
            material_id, ' although defined as constant'
        CALL WARN('electrokinetics (helmholtz_smoluchowski)',Message)
      END IF
    ELSE IF ( ElectricFieldMethod == 'computed') THEN ! get Electric Field from Electrostatic Solver
      electricField = 0.0_dp      
      DO i=1,dim
        WRITE(ValueName,'(a,i2)' ) 'Electric Field',i
        EFieldComp => VariableGet( Model % Variables, ValueName ) 
        IF (ASSOCIATED(EFieldComp)) THEN
          electricField(i) = EFieldComp%Values(EFieldComp%Perm(NodeNumber))
          ElectricFieldExists = .TRUE.
        END IF
      END DO
    ELSE         
      WRITE(Message,'(a,a,a,i0)' ) 'Unknown entry, ', ElectricFieldMethod,&
          ',for keyword >Electric Field< for Equation no. ', eq_id
      CALL WARN('electrokinetics (helmholtz_smoluchowski)',Message)
      ElectricFieldExists = .FALSE.
      hs_velocity = 0.0D00        
    END IF
  END IF
  Model % CurrentElement => BoundaryElement  ! restore correct pointer

  !---------------------------------------------
  ! compute the Helmholtz-Smoluchowski velocity
  !---------------------------------------------
  IF (ElectricFieldExists) THEN     
    norm_electricField = SUM(electricField(1:3)*Normal(1:3))
    tang_electricField = electricField - norm_electricField * Normal
    hs_velocity = tang_electricField * eoMobility
    
    IF( direction == 0 ) THEN
      CALL TangentDirections( Normal,Tangent,Tangent2)
      hs_velocity_comp = SUM(hs_velocity * Tangent )
    ELSE
      hs_velocity_comp = hs_velocity(direction)
    END IF
  END IF
  
END FUNCTION helmholtz_smoluchowski_comp



!------------------------------------------------------------------------------
!> Joule heat source as a function of electric field
!> This subroutine is basically obsolete. A more accurate version is built
!> inside the Differentials for the field Potential.
!> \deprecated Is this used any more?
!> \ingroup UDF
!-------------------------------------------------------------------------------
FUNCTION getJouleHeat( Model, NodeNumber, realDummy ) RESULT(jouleHeat)
  USE Types
  USE Lists
  USE CoordinateSystems
!
  IMPLICIT NONE
! external variables
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) :: realDummy, jouleHeat
! internal variables and parameters
  INTEGER, TARGET :: NodeNumberTarget(1)
  INTEGER, DIMENSION(:), POINTER :: NodeNumberArray
  REAL(KIND=dp) :: dummyArray(1), elConductivity, density, electricField(3), electricField_square
  INTEGER :: i, body_id, material_id, eq_id, DIM
  CHARACTER(LEN=MAX_NAME_LEN) :: ElectricFieldMethod,ValueName
  LOGICAL :: GotIt, ElectricFieldExists
  TYPE(Element_t), POINTER :: Parent
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Variable_t), POINTER :: EFieldComp


  CALL Warn('GetJouleHeat','Deprecated function!')

  DIM = CoordinateSystemDimension()
  jouleHeat = 0.0_dp

! provide arrays to let the ListGetReal/Integer function 
! work with just a nodal value to be read in
! -------------------------------------------------------
  NodeNumberTarget(1) = NodeNumber
  NodeNumberArray => NodeNumberTarget

  !-----------------------------------------------------------------
  ! get element info
  !-----------------------------------------------------------------
  body_id = Model % CurrentElement % BodyId
  eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation')
  material_id = ListGetInteger( Model % Bodies(body_id) % Values, 'Material')
  Material => Model % Materials(material_id) % Values
  IF (.NOT. ASSOCIATED(Material)) THEN
     WRITE(Message, '(a,i3)' )'No Material found for body-id ', body_id
     CALL WARN('electrokinetics (getJouleHeat)',Message)
     RETURN
  END IF
  !-----------------------------------------------------------------
  ! get material parameters
  !-----------------------------------------------------------------
  dummyArray = ListGetReal( Material,&
      'Electric Conductivity', 1, NodeNumberArray, GotIt )
  IF (.NOT. GotIt) THEN
    WRITE(Message, '(a,i0)' )&
        'No >Electric Conductivity< found in Material section ',material_id
    CALL WARN('electrokinetics (getJouleHeat)',Message)
    RETURN
  ELSE
    elConductivity = dummyArray(1)
  END IF

  dummyArray = ListGetReal( Material,'Density', 1, NodeNumberArray, GotIt )
  IF (.NOT. GotIt) THEN
    WRITE(Message, '(a,i0)' )'No >Density< found in Material section ',&
        material_id 
    CALL WARN('electrokinetics (getJouleHeat)',Message)
    CALL WARN('electrokinetics (getJouleHeat)','setting reference density to 1')
    density = 1.0D00
  ELSE
    density = dummyArray(1) 
  END IF

  !-----------------------------------------------------------------
  ! get electric field
  !-----------------------------------------------------------------
  ElectricFieldMethod = ListGetString( Model % Equations(eq_id) % &
      Values, 'Electric Field', GotIt )
  IF ( .NOT.GotIt ) THEN 
    WRITE(Message,'(a,i3)' )'No external electric field defined for Equation no.', eq_id
    CALL Info('electrokinetics (helmholtz_smoluchowski)',Message, level=4)
    ElectricFieldExists = .FALSE.
  ELSE
    ElectricFieldExists = .TRUE.
    IF ( ElectricFieldMethod == 'constant') THEN ! read field components from Material section
      ElectricFieldExists = .FALSE.
      ElectricField = 0.0_dp
      
      DO i=1,dim
        WRITE(ValueName,'(a,i2)' ) 'Electric Field',i
        dummyArray =  ListGetReal( Material,ValueName, 1, NodeNumberArray, GotIt )
        IF( GotIt ) THEN
          electricField(i) = dummyArray(1) 
          ElectricFieldExists = .TRUE.
        END IF
      END DO
      
      IF (.NOT. ElectricFieldExists) THEN
        WRITE(Message,'(a,i3)' )'No component for >Electric Field {1,2,3}< found in Material',&
            material_id, ' although defined as constant'
        CALL WARN('electrokinetics (getJouleHeat)',Message)
      END IF
    ELSE IF ( ElectricFieldMethod == 'computed') THEN ! get Electric Field from Electrostatic Solver
      electricField = 0.0_dp      
      DO i=1,dim
        WRITE(ValueName,'(a,i2)' ) 'Electric Field',i
        EFieldComp => VariableGet( Model % Variables, ValueName ) 
        IF (ASSOCIATED(EFieldComp)) THEN
          electricField(i) = EFieldComp%Values(EFieldComp%Perm(NodeNumber))
          ElectricFieldExists = .TRUE.
        END IF
      END DO
    ELSE         
      WRITE(Message,'(a,a,a,i3)' ) 'Unknown entry, ', ElectricFieldMethod,&
          ',for keyword >Electric Field< for Equation no.', eq_id
      CALL WARN('electrokinetics (getJouleHeat)',Message)
      ElectricFieldExists = .FALSE.
    END IF
  END IF


  IF (ElectricFieldExists) THEN
    electricField_square = SUM(electricField(1:3)*electricField(1:3))
    jouleHeat = elConductivity * electricField_square / density
  END IF

END FUNCTION getJouleHeat
