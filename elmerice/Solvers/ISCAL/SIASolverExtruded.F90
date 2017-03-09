!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
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
!/******************************************************************************
! *
! *  SIASolver: Solver to inquire the velocity and pressure from 
! *                   the SIA solution
! *  Exported Variables SIAFlow, dof=dim
! *
! *
! ******************************************************************************
! *
! *  Authors: Josefin Ahlkrona
! *  Email:   ahlkrona@math.uni-kiel.de
! *  Web:     http://www.csc.fi/elmer
! *
! *  Original Date: 25 June 2013 
! * 
! *****************************************************************************/



! *****************************************************************************
SUBROUTINE SIASolverExtruded( Model,Solver,dt,TransientSimulation )
  !DEC$ATTRIBUTES DLLEXPORT :: SIASolver
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  Solve the SIA Flow solution !
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
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, Grad1Sol, Grad2Sol, &
       VeloSol, FreeSurfSol, GradGrad1Sol, GradGrad2Sol, BGrad1Sol, BGrad2Sol, &
       TempGrad1Sol, TempGrad2Sol, TempSol, NormalSol, ViscositySol, dAdTSol, &
       SSAVelocity, SSAgradSol, SSAygradSol, HybridVariable,FSVeloVariable

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt, SetArrheniusFactor, ConstTemp, AllocationsBetaDone = .FALSE.

  INTEGER :: i, n, m, t, istat, DIM, p, COMP, layercounter, bc, bottomindex
  INTEGER, POINTER :: Permutation(:), VeloPerm(:), &
       GradSurface1Perm(:), GradSurface2Perm(:), &
       GradBed1Perm(:), GradBed2Perm(:), &
       FreeSurfPerm(:), GradGradSurface1Perm(:), GradGradSurface2Perm(:), &
       TempGrad1Perm(:), TempGrad2Perm(:), TempPerm(:), NormalPerm(:), &
       ViscosityPerm(:), dAdTPerm(:), SSAVelocityPerm(:), SSAgradPerm(:), &
       SSAygradPerm(:), HybridPerm(:),FSVeloPerm(:)

  INTEGER, SAVE :: Timestep

  TYPE(Variable_t), POINTER :: TimeVar

  REAL(KIND=dp), POINTER :: VariableValues(:), GradSurface1(:), &
       GradSurface2(:), Velocity(:), PrevVelo(:,:), FreeSurf(:), &
       GradGradSurface1(:),GradGradSurface2(:),  GradBed1(:), GradBed2(:), &
       Temperature(:), TempGrad1(:), TempGrad2(:), Normal(:), ViscosityVar(:), &
       dAdT(:), SSAVelocityValues(:), SSAgrad(:), SSAygrad(:), HybridValues(:), &
       FSVeloValues(:)

  REAL(KIND=dp) :: Norm, SurfGrad1,SurfGrad2, Surf,Gravity(3), &
       Integral,SurfGrad1Grad1,SurfGrad2Grad1,BedGrad1, BedGrad2, &
       SurfGrad1Grad2,SurfGrad2Grad2, Viscosity, n_powerlaw, Slip1, Slip2, &
       Position, PositionUnder

  REAL(KIND=dp), ALLOCATABLE :: SIAxVelocity(:), SIAyVelocity(:), SIApressure(:), &
       A3hminusz3(:),dvxdx(:),dvydy(:),threeA3hminusz2(:),intdvxdx(:), &
       dvssadx(:), dvssady(:), intdvxdxssa(:)

  REAL(KIND=dp), ALLOCATABLE :: ArrheniusFactor(:), dArrheniusFactordT(:), dAdxhminusz3(:), &
       dAdyhminusz3(:)

  REAL(KIND=dp), ALLOCATABLE :: Beta1(:), Beta2(:), Beta1Node(:), &
       Beta2Node(:)

  REAL(KIND=dp) :: uB1, uB2, uB3, up !basal velocities
    REAL(KIND=dp) :: BasalVeloX, BasalVeloY !basal velocities
      REAL(KIND=dp) :: Tangent(3),Tangent2(3)
  REAL(KIND=dp) :: VELO_LIMIT, veloratio !maximum absolute velocity
  LOGICAL :: FOUND_VELO_LIMIT


  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, TemperatureName, ViscosityFlag, &
       BedrockName, SurfaceName

  SAVE  AllocationsDone, DIM, SolverName

  REAL(KIND=dp) :: g, rho, nGlen, nGleninv, siatime

  REAL(KIND=dp) :: Temp, R, Tlimit, A1, A2, Q1, Q2, vx, vy, vz

  TYPE(Element_t),POINTER :: Element

  INTEGER, ALLOCATABLE :: LayerOfNode(:)

  LOGICAL :: BedrockData, HybridSolution


  !-----------Variables needed for integration --------------------------------
  TYPE(Solver_t), POINTER :: PSolver
  INTEGER :: j,k,l,Dofs,dof,nsize,TopNodes,BotNodes
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
  LOGICAL :: Initialized = .FALSE.,GotVar
  REAL(KIND=dp) :: dx,Level,q
  REAL(KIND=dp), POINTER :: Coord(:)
  TYPE(Variable_t), POINTER :: Var, OldVar

  !REAL(KIND=dp) :: CPUTime,RealTime

  !---------------------------------------------------------------------------
  LOGICAL :: ComputeSurface = .FALSE.
  TYPE(Variable_t), POINTER :: ThickSol 
  INTEGER, POINTER :: ThickPerm(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: ThickName
  REAL(KIND=dp), POINTER :: Thick(:)

  !-----Needed for transient stuff etc.

  SAVE Material, ArrheniusFactor,  dArrheniusFactordT, A3hminusz3, &
       BotPointer, TopPointer, UpPointer, DownPointer, SIApressure, &
       Coord, nsize,  threeA3hminusz2, SurfGrad1,SurfGrad2, &
       SurfGrad1Grad1,SurfGrad2Grad1,BedGrad1, BedGrad2, Surf, dvxdx, &
       dvydy, intdvxdx,dAdxhminusz3,dAdyhminusz3, SetArrheniusFactor, &
       LayerOfNode, Beta1Node, Beta2Node, Beta1,Beta2, Slip1, Slip2, &
       dvssadx,dvssady, intdvxdxssa


	SolverName = "SIA_Solver"

  siatime=CPUTime()

  !------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  N = Solver % Mesh % MaxElementDOFs !needed right now only for beta

WRITE(*,*) 'sia no bulkd elems', Solver % Mesh % NumberOfBulkElements
WRITE(*,*) 'sia no bnd elems', Solver % Mesh % NumberOfBoundaryElements
WRITE(*,*)  'sia active elements',Solver % NumberOfActiveElements

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     IF (AllocationsDone) DEALLOCATE(SIAxVelocity, SIAyVelocity,&
          SIApressure,A3hminusz3,dvxdx,dvydy,threeA3hminusz2,intdvxdx,&
          ArrheniusFactor, dArrheniusFactordT,dAdxhminusz3,dAdyhminusz3, &
          LayerOfNode, dvssadx, dvssady,intdvxdxssa)

     ALLOCATE(SIAxVelocity(Model % Mesh % NumberOfNodes), &
          SIAyvelocity(Model % Mesh % NumberOfNodes), &
          SIApressure(Model % Mesh % NumberOfNodes), &
          A3hminusz3(Model % Mesh % NumberOfNodes), &
          dvxdx(Model % Mesh % NumberOfNodes), &
          dvydy(Model % Mesh % NumberOfNodes), &
          dvssadx(Model % Mesh % NumberOfNodes), &
          dvssady(Model % Mesh % NumberOfNodes), &
          threeA3hminusz2(Model % Mesh % NumberOfNodes), &
          intdvxdx(Model % Mesh % NumberOfNodes),&
          intdvxdxssa(Model % Mesh % NumberOfNodes),&
          ArrheniusFactor(Model % Mesh % NumberOfNodes), &
          dArrheniusFactordT(Model % Mesh % NumberOfNodes),&
          dAdxhminusz3(Model % Mesh % NumberOfNodes), &
          LayerOfNode(Model % Mesh % NumberOfNodes), &
          dAdyhminusz3(Model % Mesh % NumberOfNodes), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  !    Get Timestep
  TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
  Timestep = NINT(Timevar % Values(1))


  !------------------------------------------------------------------------------
  !   Initialize the pointers to top and bottom nodes, needed for the integration 
  !------------------------------------------------------------------------------
  IF( .NOT. Initialized ) THEN

     ! Choose active direction coordinate and set corresponding unit vector
     !---------------------------------------------------------------------
     PSolver => Solver
     CALL DetectExtrudedStructure( Solver % Mesh, PSolver, Var, &
          TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
          UpNodePointer = UpPointer, DownNodePointer = DownPointer )

     Coord => Var % Values
     nsize = SIZE( Coord )
     Initialized = .TRUE.

  END IF

  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  VeloSol => VariableGet( Solver % Mesh % Variables, 'SIAFlow' )
  IF (ASSOCIATED(veloSol)) THEN
     Velocity => VeloSol % Values
     VeloPerm => VeloSol % Perm
     PrevVelo => VeloSol % PrevValues
  ELSE
     CALL FATAL(SolverName,'Could not find variable >SIAFlow<')
  END IF


  !Get name of surface
  SurfaceName=GetString( Solver % Values, 'Surface Name', gotIt )


  !get free surface
  FreeSurfSol => VariableGet( Solver % Mesh % Variables, TRIM(SurfaceName))

  IF (ASSOCIATED(FreeSurfSol)) THEN
     FreeSurf => FreeSurfSol % Values
     FreeSurfPerm => FreeSurfSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find variable >'//TRIM(SurfaceNAme)//'<')
  END IF

  !get free surface gradient

  Grad1Sol => VariableGet( Solver % Mesh % Variables, TRIM(SurfaceName)//' grad 1')

  IF  (ASSOCIATED(Grad1Sol)) THEN
     GradSurface1 => Grad1Sol % Values
     GradSurface1Perm => Grad1Sol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find a variable >'//TRIM(SurfaceNAme)//' grad 1<')
  END IF
  
  HybridSolution=GetLogical( Solver % Values, 'Hybrid Solution', gotIt )
  IF ( HybridSolution ) THEN
  
  SSAgradSol => VariableGet( Solver % Mesh % Variables, TRIM('SSAVelocity 1')//' grad 1')

  IF  (ASSOCIATED(SSAgradSol)) THEN
     SSAgrad => SSAgradSol % Values
     SSAgradPerm => SSAgradSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find a variable >'//TRIM('SSAVelocity 1')//' grad 1<')
  END IF

  END IF

  IF (dim > 2) THEN
    IF ( HybridSolution ) THEN
      SSAygradSol => VariableGet( Solver % Mesh % Variables, TRIM('SSAVelocity 2')//' grad 2')

     IF  (ASSOCIATED(SSAygradSol)) THEN
      SSAygrad => SSAygradSol % Values
      SSAygradPerm => SSAygradSol % Perm
     ELSE
      CALL FATAL(SolverName,'Could not find a variable >'//TRIM('SSAVelocity 2')//' grad 2<')
     END IF
    END IF

     Grad2Sol => VariableGet( Solver % Mesh % Variables, TRIM(SurfaceName)//' grad 2')

     IF  (ASSOCIATED(Grad2Sol)) THEN
        GradSurface2 => Grad2Sol % Values
        GradSurface2Perm => Grad2Sol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find a variable >'//TRIM(SurfaceNAme)//' grad 2<')
     END IF

  END IF

  !get second derivatives of surface

  GradGrad1Sol => VariableGet( Solver % Mesh % Variables, TRIM(SurfaceName)//' grad 1 grad')
  IF  (ASSOCIATED(GradGrad1Sol)) THEN
     GradGradSurface1 => GradGrad1Sol % Values
     GradGradSurface1Perm => GradGrad1Sol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find a variable >'//TRIM(SurfaceName)//' grad 1 grad<')
  END IF

  IF (dim > 2) THEN
     GradGrad2Sol => VariableGet( Solver % Mesh % Variables, TRIM(SurfaceName)//' grad 2 grad')
     IF  (ASSOCIATED(GradGrad2Sol)) THEN
        GradGradSurface2 => GradGrad2Sol % Values
        GradGradSurface2Perm => GradGrad2Sol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find a variable >'//TRIM(SurfaceName)//' grad 2 grad<')
     END IF
  END IF
  
  
      FSVeloVariable => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
     IF ( ASSOCIATED( FSVeloVariable ) ) THEN
        FSVeloPerm    => FSVeloVariable % Perm
        FSVeloValues  => FSVeloVariable % Values
     ELSE
        CALL Fatal( 'Approximation updater','Cannot find variable <Flow Solution>' )
     END IF

!get bedrock and bedrock gradients
  BedrockData=GetLogical( Solver % Values, 'Bedrock Data', gotIt )

  IF ( BedrockData ) THEN
     BedrockName=GetString( Solver % Values, 'Bedrock Name', gotIt )
     BGrad1Sol => VariableGet( Solver % Mesh % Variables, 'FreeBedGrad1')
     IF (.NOT. ASSOCIATED(BGrad1Sol)) THEN
        BGrad1Sol => VariableGet( Solver % Mesh % Variables, TRIM(BedrockName)//' grad 1')
     END IF
     IF (ASSOCIATED(BGrad1Sol)) THEN
        GradBed1 => BGrad1Sol % Values
        GradBed1Perm => BGrad1Sol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find variable >'//TRIM(BedrockName)//' grad 1<')
     END IF

     IF (dim > 2) THEN
        BGrad2Sol => VariableGet( Solver % Mesh % Variables, 'FreeBedGrad2')  
        IF (.NOT. ASSOCIATED(BGrad2Sol)) THEN
           BGrad2Sol => VariableGet( Solver % Mesh % Variables, TRIM(BedrockName)//' grad 2')
        END IF
        IF (ASSOCIATED(BGrad2Sol)) THEN
           GradBed2 => BGrad2Sol % Values
           GradBed2Perm => BGrad2Sol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >'//TRIM(BedrockName)//' grad 2<')
        END IF
     END IF
  ELSE
      WRITE( Message, * ) 'Bedrock data zero'
     CALL Info('SIASolver',Message,Level=4)
  END IF



  ComputeSurface=GetLogical( Solver % Values, 'Compute Surface From Thickness', gotIt )

  IF (ComputeSurface) THEN
    ! WRITE(*,*) 'compute surface'
     ThickName=GetString( Solver % Values, 'Thickness Name', gotIt )
     ThickSol => VariableGet( Solver % Mesh % Variables, TRIM(ThickName))
     IF (ASSOCIATED(ThickSol)) THEN
        Thick => ThickSol % Values
        ThickPerm => ThickSol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find thickness variable')
     END IF

     DO i = 1,Model % Mesh % NumberOfNodes
        IF (DIM==2) THEN
           FreeSurf(FreeSurfPerm(i))=Thick(ThickPerm(i))+ &
                Model % Nodes % y(BotPointer(i))
        ELSEIF (DIM==3) THEN
           FreeSurf(FreeSurfPerm(i))=Thick(ThickPerm(i))+ &
                Model % Nodes % z(BotPointer(i))
        END IF
     END DO

  END IF

  VELO_LIMIT=GetConstReal( Solver % Values, 'Velocity Cutoff', FOUND_VELO_LIMIT)

  VariableValues = 0.0d0
  Norm = Solver % Variable % Norm

  IF (DIM==2) THEN 
     g = abs(GetConstReal(Model % BodyForces(1) % Values, "Flow BodyForce 2", Found))
  ELSE
     g = abs(GetConstReal(Model % BodyForces(1) % Values, "Flow BodyForce 3", Found))
  END IF

  !------------------------------------------------------------------------------
  !   Viscosity reading
  !------------------------------------------------------------------------------

  !Material => GetMaterial() 
  Material =>  Model % Materials(1) % Values

  ViscosityFlag = ListGetString( Model % Materials(1) % Values ,'Viscosity Model', GotIt) 

  SELECT CASE( ViscosityFlag )

  CASE ('glen') !! need updating for non-isothermal
     CALL INFO('SIASolver','viscosity model: glen', Level=3)
     ConstTemp=.TRUE.

     CALL INFO('SIASolver','Setting Arrhenius Factor', Level = 5)
     ArrheniusFactor = GetConstReal(Material,'Arrhenius Factor', GotIt)
     dArrheniusFactordT = 0.0_dp

     nGlen = GetConstReal(Model % Materials(1) % Values, "Glen Exponent", Found)

  CASE('power law')
     CALL INFO('SIASolver','viscosity model: power law', Level=5)
     ConstTemp = GetLogical(Material, 'Isothermal', GotIt)

     n_powerlaw = GetConstReal(Material, 'Viscosity Exponent', Found)
     nGlen=AINT(100.0*1.0/n_powerlaw)/100.0 !Got weird errors before
     
    ! WRITE(*,*) 'tjoho'

     IF (ConstTemp) THEN
        Viscosity = GetConstReal(Material, 'Viscosity', Found)
        ArrheniusFactor = 1.0/(2.0*Viscosity**nGlen)
        dArrheniusFactordT = 0.0_dp
            ! WRITE(*,*) 'const temp'
     ELSE !Variable temperature
             !WRITE(*,*) 'variable temp'
        DO i = 1, Model % Mesh % NumberOfNodes
           Viscosity = ListGetRealAtNode(Material,'Viscosity', i, GotIt)
           ArrheniusFactor(i) = 1.0/(2.0*Viscosity**nGlen)
           dArrheniusFactordT(i) =  ListGetRealAtNode(Material, &
                'dArrheniusFactordT', i, GotIt)           
        END DO

 	TempGrad1Sol => VariableGet(  Solver % Mesh % Variables, 'Temp grad 1')
	IF (ASSOCIATED(TempGrad1Sol)) THEN
           TempGrad1 => TempGrad1Sol % Values
           TempGrad1Perm => TempGrad1Sol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Temp grad 1<')
        END IF

        IF (DIM==3) THEN
           TempGrad2Sol => VariableGet(  Solver % Mesh % Variables, 'Temp grad 2')
           IF (ASSOCIATED(TempGrad2Sol)) THEN
              TempGrad2 => TempGrad2Sol % Values
              TempGrad2Perm => TempGrad2Sol % Perm
           ELSE
              CALL FATAL(SolverName,'Could not find variable >Temp grad 2<')
           END IF
        END IF
     END IF
  END SELECT

  rho = GetConstReal(Model % Materials(1) % Values, "Density", Found)

  !  nGleninv = GetConstReal(Model % Materials(1) % Values, "Viscosity Exponent", Found)
  !  nGlen=1.0_dp/nGleninv

  !--------------------------------------------------------------
  ! Compute SIA-solution
  !--------------------------------------------------------------

!WRITE(*,*) 'Starting integration'
  !------------------------------------------------------------------------------
  !   Integral of 2A*(rho*g)^n*(h-z)^n
  !------------------------------------------------------------------------------
  !Integration from bottom to z-level, trapezoidal method
  DO i=1,nsize
     IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
        l = i !start at the bottom of the line
        A3hminusz3(l) = 0.0_dp !Integral is zero at the bottom
        DO k=1,nsize

           !              Surf = FreeSurf(FreeSurfPerm(TopPointer(l)))
           IF (dim == 2) THEN
              Surf= Model % Nodes % y(TopPointer(l))
              Position =  Model % Nodes % y(l)
              PositionUnder=Model % Nodes % y(DownPointer(l))
           ELSEIF (dim ==3) THEN
              Surf= Model % Nodes % z(TopPointer(l))
              Position = Model % Nodes % z(l)
              PositionUnder= Model % Nodes % z(DownPointer(l))
           END IF

           IF( k > 1 ) THEN !above bottom

              dx = (Coord(l) - Coord(DownPointer(l)))

              A3hminusz3(l)=A3hminusz3(DownPointer(l))+0.5_dp*dx*( &
                   2*ArrheniusFactor(l)*(rho*g)**nGlen*& 
                   (Surf - Position)**nGlen + &
                   2*ArrheniusFactor(DownPointer(l))*(rho*g)**nGlen*&
                   (Surf - PositionUnder)**nGlen)
           END IF

           IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
           l = UpPointer(l)   !step up    

        END DO
     END IF
  END DO

!WRITE(*,*) 'Starting integration 2'
  !!------------------------------------------------------------------------------
  !   Integral of 2*dAdx*(rho*g)^3*(h-z)^3
  !------------------------------------------------------------------------------


  IF (.NOT. ConstTemp) THEN

     !Integration from bottom to z-level, trapezoidal method
     DO i=1,nsize
        IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
           l = i !start at the bottom of the line
           dAdxhminusz3(l) = 0.0_dp !Integral is zero at the bottom
           DO k=1,nsize

              IF (dim == 2) THEN
                 Surf= Model % Nodes % y(TopPointer(l))
                 Position =  Model % Nodes % y(l)
                 PositionUnder=Model % Nodes % y(DownPointer(l))
              ELSEIF (dim ==3) THEN
                 Surf= Model % Nodes % z(TopPointer(l))
                 Position = Model % Nodes % z(l)
                 PositionUnder= Model % Nodes % z(DownPointer(l))
              END IF

              IF( k > 1 ) THEN !above bottom
                 dx = (Coord(l) - Coord(DownPointer(l)))
                 dAdxhminusz3(l)=dAdxhminusz3(DownPointer(l))+0.5_dp*dx*( &
                      2*dArrheniusFactordT(l)*TempGrad1(TempGrad1Perm(l))*(rho*g)**nGlen*& 
                      (Surf - Position)**nGlen + &
                      2*dArrheniusFactordT(DownPointer(l))*TempGrad1(TempGrad1Perm(DownPointer(l)))*(rho*g)**nGlen*&
                      (Surf - PositionUnder)**nGlen)
              END IF

              IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
              l = UpPointer(l)   !step up    

           END DO
        END IF
     END DO
  ELSE
     dAdxhminusz3 = 0.0_dp 
  END IF

!WRITE(*,*) 'Starting integration 3'
  !------------------------------------------------------------------------------
  !   integral of 3*2A*(rho*g)^3*(h-z)^2
  !------------------------------------------------------------------------------

  !Integration from bottom to z-level, trapezoidal method
  DO i=1,nsize
     IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
        l = i !start at the bottom of the line
        threeA3hminusz2(l) = 0.0_dp !Integral is zero at the bottom
        DO k=1,nsize

           !Surf = FreeSurf(FreeSurfPerm(TopPointer(l)))
           IF (dim == 2) THEN
              Surf= Model % Nodes % y(TopPointer(l))
              Position =  Model % Nodes % y(l)
              PositionUnder=Model % Nodes % y(DownPointer(l))
           ELSEIF (dim ==3) THEN
              Surf= Model % Nodes % z(TopPointer(l))
              Position = Model % Nodes % z(l)
              PositionUnder= Model % Nodes % z(DownPointer(l))
           END IF

           IF( k > 1 ) THEN !above bottom
              dx = (Coord(l) - Coord(DownPointer(l)))
              threeA3hminusz2(l)=threeA3hminusz2(DownPointer(l))+3.0_dp*0.5_dp*dx*( &
                   2*ArrheniusFactor(l)*(rho*g)**nGlen*& 
                   (Surf - Position)**(nGlen-1) + &
                   2*ArrheniusFactor(DownPointer(l))*(rho*g)**nGlen*&
                   (Surf - PositionUnder)**(nGlen-1))
           END IF

           IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
           l = UpPointer(l)   !step up    

        END DO
     END IF
  END DO

!WRITE(*,*) 'pressure'
  !------------------------------------------------------------------------------
  !   Pressure
  !------------------------------------------------------------------------------

  SIApressure = 0.0_dp 
  DO i=1,nsize

     IF (dim == 2) THEN
        Surf= Model % Nodes % y(TopPointer(i))
        Position =  Model % Nodes % y(i)
     ELSEIF (dim ==3) THEN
        Surf= Model % Nodes % z(TopPointer(i))
        Position = Model % Nodes % z(i)
     END IF

     SIApressure(i)=rho*g*(Surf - Position)
  END DO

!WRITE(*,*) 'integration 4'

  IF (dim == 3) THEN!--------------------------------------------------------------------------------------

     !------------------------------------------------------------------------------
     !   Integral of 2*dAdy*(rho*g)^3*(h-z)^3
     !------------------------------------------------------------------------------


     IF (.NOT. ConstTemp) THEN
        !Integration from bottom to z-level, trapezoidal method
        DO i=1,nsize
           IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
              l = i !start at the bottom of the line
              dAdyhminusz3(l) = 0.0_dp !Integral is zero at the bottom
              DO k=1,nsize

                 ! Surf = FreeSurf(FreeSurfPerm(TopPointer(l)))
                 Surf =  Model % Nodes % z(TopPointer(l))

                 IF( k > 1 ) THEN !above bottom
                    dx = (Coord(l) - Coord(DownPointer(l)))
                    dAdyhminusz3(l)=dAdyhminusz3(DownPointer(l))+0.5_dp*dx*( &
                         2*dArrheniusFactordT(l)*TempGrad2(TempGrad2Perm(l))*(rho*g)**nGlen*& 
                         (Surf - Model % Nodes % z(l))**nGlen + &
                         2*dArrheniusFactordT(DownPointer(l))*TempGrad2(TempGrad2Perm(DownPointer(l)))*(rho*g)**nGlen*&
                         (Surf - Model % Nodes % z(DownPointer(l)))**nGlen)
                 END IF

                 IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
                 l = UpPointer(l)   !step up    

              END DO
           END IF
        END DO
     ELSE !constant temperature
        dAdyhminusz3 = 0.0_dp 
     END IF
  ELSE !2D case
     dAdyhminusz3 = 0.0_dp 
  END IF

  !------------------------------------------------------------------------------
  !   dvxdx, dvydy
  !------------------------------------------------------------------------------
!WRITE(*,*) 'dvxdx'

  dvxdx=0.0
  dvydy=0.0
  BedGrad1=0.0_dp 
  BedGrad2=0.0_dp

  SurfGrad2=0.0
  SurfGrad2Grad2=0.0
  SurfGrad2Grad1=0.0
  SurfGrad1Grad2=0.0

  DO i = 1, Model % Mesh % NumberOfNodes

     IF (dim == 2) THEN
        Surf= Model % Nodes % y(TopPointer(i))
        Position =  Model % Nodes % y(i)
     ELSEIF (dim ==3) THEN
        Surf= Model % Nodes % z(TopPointer(i))
        Position = Model % Nodes % z(i)
     END IF

     SurfGrad1 = GradSurface1(GradSurface1Perm(TopPointer(i)))
     SurfGrad1Grad1=GradGradSurface1(DIM*(GradSurface1Perm(TopPointer(i))-1)+1)
     
     IF ( HybridSolution ) THEN
      dvssadx(i)= SSAgrad(SSAgradPerm(BotPointer(i)))
      dvssady(i)= 0.0!SSAgrad(SSAgradPerm(BotPointer(i)))
     ELSE 
      dvssadx(i)= 0.0
      dvssady(i)= 0.0
     END IF

     IF (dim ==3) THEN
        SurfGrad2= GradSurface2(GradSurface2Perm(TopPointer(i)))
        SurfGrad2Grad2=GradGradSurface2(DIM*(GradSurface2Perm(TopPointer(i))-1)+2)
        SurfGrad2Grad1=GradGradSurface2(DIM*(GradSurface2Perm(TopPointer(i))-1)+1)
        SurfGrad1Grad2=GradGradSurface1(DIM*(GradSurface1Perm(TopPointer(i))-1)+2)
     END IF

     IF (BedrockData) THEN
        BedGrad1 = GradBed1(GradBed1Perm(BotPointer(i)))
        IF (dim ==3) THEN
           BedGrad2 = GradBed2(GradBed2Perm(BotPointer(i)))
        END IF
     END IF

     dvxdx(i)=  -SurfGrad1*(2.0_dp*SurfGrad1*SurfGrad1Grad1 + 2.0_dp*SurfGrad2*SurfGrad2Grad1)**((nGlen-1.0)/2.0)*A3hminusz3(i) & 
          -SurfGrad1Grad1*SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*A3hminusz3(i)  &
          -SurfGrad1* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*SurfGrad1*threeA3hminusz2(i) &
          + SurfGrad1* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)&
          *BedGrad1*2.0*ArrheniusFactor(i)*(rho*g)**nGlen*(Surf - Model % Nodes % z(i))**nGlen  & !leibniz term
          -SurfGrad1* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*dAdxhminusz3(i)

     IF (dim ==3) THEN
        dvydy(i)= -SurfGrad2*(2.0_dp*SurfGrad1*SurfGrad1Grad2 + 2.0_dp*SurfGrad2*SurfGrad2Grad2)**((nGlen-1.0)/2.0)*A3hminusz3(i) & 
             -SurfGrad2Grad2*SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*A3hminusz3(i)  &
             -SurfGrad2* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*SurfGrad2*threeA3hminusz2(i) &
             + SurfGrad2* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)&
             * BedGrad2*2.0*ArrheniusFactor(i)*(rho*g)**nGlen*(Surf - Model % Nodes % z(i))**nGlen  & !leibniz term
             -SurfGrad2* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*dAdyhminusz3(i)
     END IF

  END DO

!WRITE(*,*) 'dvzdz'
  !------------------------------------------------------------------------------
  !   z-Velocity
  !------------------------------------------------------------------------------
  intdvxdx = 0.0_dp
  intdvxdxssa= 0.0_dp
  !Integration from bottom to z-level, trapezoidal method
  DO i=1,nsize
     IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
        l = i !start at the bottom of the line
        intdvxdx(l) = 0.0_dp !Integral is zero at the bottom
        DO k=1,nsize

           IF( k > 1 ) THEN !above bottom
              dx = (Coord(l) - Coord(DownPointer(l)))
              intdvxdx(l)=intdvxdx(DownPointer(l))+0.5_dp*dx*( &
                   dvxdx(l) +  dvxdx(DownPointer(l)) + dvydy(l) +  dvydy(DownPointer(l)))
              intdvxdxssa(l)=intdvxdxssa(DownPointer(l))+0.5_dp*dx*( &
                   dvssadx(BotPointer(l))+ dvssadx(BotPointer(l)) + &
                   dvssady(BotPointer(l)) +  dvssady(BotPointer(l)))
           END IF

           IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
           l = UpPointer(l)   !step up

        END DO

     END IF
  END DO



  !First a try run to get what bc has the slip coefficient, if it has it
  bottomindex=0
  DO bc=1,Model % NumberOfBCs
     BasalVeloX=ListGetRealAtNode( Model % BCs(bc) % Values, 'Basal Velocity X', 1, GotIt)
     IF (.NOT. GotIt) THEN
       Slip1 = ListGetRealAtNode( Model % BCs(bc) % Values, 'Slip Coefficient 2', 1, GotIt)
     END IF
     IF (GotIt) THEN
        bottomindex = bc
     END IF
  END DO


  IF (bottomindex .NE. 0) THEN !we have a slip-coefficient set  , wo we wanna do a bc

     ! We need the normal vector to do the bc computation
     NormalSol => VariableGet( Solver % Mesh % Variables, 'Normal Vector' )
     IF (ASSOCIATED(NormalSol)) THEN
        Normal => NormalSol % Values
        NormalPerm => NormalSol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find variable >NormalSol<')
     END IF

  END IF !bottomindex .NE. 0

  !------------------------------------------------------------------------------
  ! Add sliding and save the solution on the right variable
  !------------------------------------------------------------------------------
  open(unit=119, file='dhdx.txt',STATUS='REPLACE')
  open(unit=120, file='viskositetsia.txt',STATUS='REPLACE')


  IF (dim == 2) THEN

     DO i = 1, Model % Mesh % NumberOfNodes
        IF (VeloPerm(i)>0) THEN
           SurfGrad1 = GradSurface1(GradSurface1Perm(TopPointer(i)))
           SurfGrad2= 0.0_dp   

           vx= -SurfGrad1* ABS(SurfGrad1)**(nGlen-1.0)*A3hminusz3(i)
           vy= -intdvxdx(i)

           Velocity ((DIM+1)*(VeloPerm(i)-1) + 1) = vx
           Velocity ((DIM+1)*(VeloPerm(i)-1) + 2) = vy
           Velocity ((DIM+1)*(VeloPerm(i)-1) + 3) = SIApressure(i)

           WRITE(119,*) SurfGrad1

           Surf= Model % Nodes % y(TopPointer(i))
           Position =  Model % Nodes % y(i)


           WRITE(120,*) 1.0/(2.0*ArrheniusFactor(i)*(rho*g*(Surf - Position)*ABS(SurfGrad1))**2.0)
        END IF
     END DO



  ELSE IF (dim == 3) THEN

     DO i = 1, Model % Mesh % NumberOfNodes
        IF (VeloPerm(i)>0) THEN


           !------------------------------------------------------------------------------

           SurfGrad1 = GradSurface1(GradSurface1Perm(TopPointer(i)))
           SurfGrad2= GradSurface2(GradSurface2Perm(TopPointer(i)))
           Surf = Model % Nodes % z(TopPointer(i))

           !------------------------------------------------------------------------------
           ! Sliding
           !------------------------------------------------------------------------------

           !Zero if not at bottom AND sliding
           uB1=0
           uB2=0
           uB3=0
           IF (bottomindex .NE. 0) THEN  !sliding
              
              j = BotPointer(i) !sliding is at bottom but added everywhere in the column
              BasalVeloX=ListGetRealAtNode( Model % BCs(bottomindex) % Values, 'Basal Velocity X', j, GotIt)
              BasalVeloY=ListGetRealAtNode( Model % BCs(bottomindex) % Values, 'Basal Velocity Y', j, GotIt)
             
             
             up = ListGetRealAtNode( Model % BCs(bottomindex) % Values, 'External Pressure', j, GotIt)
             up=0
              IF (GotIt) THEN
                 !Slipz = ListGetRealAtNode( Model % BCs(bottomindex) % Values, &
                 !  'Slip Coefficient 1', j, GotIt)
                   
               !IF (GotIt) THEN
               !  WRITE(*,*) 'need t implement seaspring'
               !END IF
              
                 CALL TangentDirections((/ Normal(DIM*(NormalPerm(j)-1)+1), &
                   Normal(DIM*(NormalPerm(j)-1)+2),Normal(DIM*(NormalPerm(j)-1)+3)/), Tangent, Tangent2 ) 
                uB1=BasalVeloX
                uB2=BasalVeloY
                uB3=0.0
	
              	!uB1=BasalVeloX*Tangent2(1)+BasalVeloY*Tangent(1)
              	!uB2=BasalVeloX*Tangent2(2)+BasalVeloY*Tangent(2)
              	!uB3=BasalVeloX*Tangent2(3)+BasalVeloY*Tangent(3)
           	
              ELSE


              Slip1 = ListGetRealAtNode( Model % BCs(bottomindex) % Values, &
                   'Slip Coefficient 2', j, GotIt)

              Slip2 = ListGetRealAtNode( Model % BCs(bottomindex) % Values, &
                   'Slip Coefficient 3', j, GotIt)
                   

              CALL LinearSliding( Slip1, Slip2, Normal(DIM*(NormalPerm(j)-1)+1), &
                   Normal(DIM*(NormalPerm(j)-1)+2),Normal(DIM*(NormalPerm(j)-1)+3), &
                   - rho*g*SurfGrad1*(Surf - Model % Nodes % z(j)), - rho*g*SurfGrad2*(Surf - Model % Nodes % z(j)), &
   rho*g*(Surf - Model % Nodes % z(j)), &
                   uB1,uB2,uB3) 
                   END IF !if not slip coeff
           END IF !if at bottom


!uB1= 0.0
!uB2 = 0.0
!uB3 = 0.0
           vx=  uB1-SurfGrad1* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*A3hminusz3(i)
           vy = uB2 -SurfGrad2* SQRT(SurfGrad1**2.0 + SurfGrad2**2.0)**(nGlen-1.0)*A3hminusz3(i)
           vz = uB3 -intdvxdx(i) -intdvxdxssa(i)
        
      HybridVariable => VariableGet( Solver % Mesh % Variables, 'ApproximationHybrid' )
     IF ( ASSOCIATED( HybridVariable ) ) THEN
        HybridPerm    => HybridVariable % Perm
        HybridValues  => HybridVariable % Values
     ELSE
        CALL Fatal( 'FlowSolveSIAFS','Cannot find variable <ApproximationHybrid>' )
     END IF
     
     HybridValues(HybridPerm(i))=&
     SQRT((uB1**2.0+uB2**2.0+uB3**2.0+(intdvxdxssa(i))**2.0)&
     /(vx**2.0+vy**2.0+vz**2.0))
           
	   IF (FOUND_VELO_LIMIT) THEN !limit the velocity magnitude
              IF ((vx**2.0+vy**2.0+vz**2.0) > VELO_LIMIT**2.0 ) THEN
                        veloratio = VELO_LIMIT**2.0/(vx**2.0+vy**2.0+vz**2.0)
                         vx=vx*veloratio
                         vy=vy*veloratio
                         vz=vz*veloratio
              END IF
           END IF

           Velocity ((DIM+1)*(VeloPerm(i)-1) + 1) = vx
           Velocity ((DIM+1)*(VeloPerm(i)-1) + 2) = vy
           Velocity ((DIM+1)*(VeloPerm(i)-1) + 3) = vz
           Velocity ((DIM+1)*(VeloPerm(i)-1) + 4) =  SIApressure(i) +up
           
 
     
      Velocity ((DIM+1)*(VeloPerm(i)-1) + 4) = FSVeloValues ((DIM+1)*(FSVeloPerm(i)-1) + 4)

           Position =  Model % Nodes % z(i)

           WRITE(120,*) 1.0/(2.0*ArrheniusFactor(i)*(rho*g*(Surf - Position)*SQRT(SurfGrad1**2.0+SurfGrad2**2.0))**2.0)


        END IF
     END DO

  END IF
  
   

  siatime=CPUTime()-siatime

  open (unit=135, file="TimingSIA.dat",POSITION='APPEND')

  WRITE(135,*)  'At ', RealTime(), ' ******** Timestep = ', Timestep, ' ********'

  WRITE(135,*)  'SIA Solver time: ', siatime
  WRITE(135,*) '***************************************************************'
  WRITE(135,*) '                                                               '

  close(135)


  close(119)


  !------------------------------------------------------------------------------
END SUBROUTINE SIASolverExtruded
!------------------------------------------------------------------------------

