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
! *  Authors: F. GILLET-CHAULET
! *  Email:   fabien.gillet-cahulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: Nov. 2017
! * 
! *****************************************************************************
!!! Apply the flotation criterion to update Zb and Zs from the ice thickness H
!!!    Execute this solver on the bottom boundary where the thickness equation is solved
!!!    IF the mesh is vertically extruded export Zs on the top surface.
!
!   OUTPUT Variables:
!     Zb
!     Zs
!     DZbDt (optional, if variable found and transient simulation)
!     DZsDt (optional, if variable found and transient simulation)
!   
!   INPUT Variable:
!     H
!     bedrock (optional)
!
! PARAMETERS:
!   Constants: 
!     zsea
!     rhow
!   Material:
!      SSA Mean Density
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  Flotation_init( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !--------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: ZbName,ZsName,HName
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Flotation'
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: GotIt
  
  SolverParams => Solver % Values 

  ZbName = GetString(SolverParams, 'Bottom Surface Name', GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO(SolverName, 'Bottom Surface Name not found - using default Zb', level=3)
    CALL ListAddString(SolverParams,'Bottom Surface Name','Zb')
  END IF

  ZsName = GetString(SolverParams, 'Top Surface Name', GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO(SolverName, 'Top Surface Name not found - using default Zs', level=3)
    CALL ListAddString(SolverParams,'Top Surface Name','Zs')
  END IF

  HName = GetString(SolverParams, 'Thickness Variable Name', GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO(SolverName, 'Thickness Variable Name not found - using default H', level=3)
    CALL ListAddString(SolverParams,'Thickness Variable Name','H')
  END IF

END SUBROUTINE  Flotation_init
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Flotation( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Solver_t),POINTER :: PSolver
  TYPE(Variable_t),POINTER :: Var
  TYPE(Variable_t),POINTER :: ZbVar,ZsVar
  TYPE(Variable_t),POINTER :: HVar,BedVar
  TYPE(Variable_t),POINTER :: GLMask,HafVar
  TYPE(Variable_t),POINTER :: sftgif,sftgrf,sftflf
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t),POINTER :: BodyForce,Material, Params
  TYPE(Nodes_t),SAVE :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IP

  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: Density
  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: GL
  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: MinH,NodalH
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:)
  REAL(KIND=dp) :: zsea,rhow,rhoi
  REAL(KIND=dp) :: H,zb,zs,bedrock,Hf
  REAL(KIND=dp), PARAMETER :: EPS=EPSILON(1.0)
  REAL(KIND=dp) :: FillValue
  REAL(KIND=dp) :: flarea,grarea,area,IParea,detJ


  INTEGER,DIMENSION(:),POINTER,SAVE :: BotPointer,TopPointer,UpPointer
  INTEGER, POINTER,SAVE :: NodeIndexes(:)
  INTEGER :: GroundedNode
  INTEGER :: ActiveDirection
  INTEGER :: t,i,n,kk,ll
  INTEGER :: topnode
  INTEGER :: Active
  INTEGER :: Eindex
  INTEGER :: GlnIP

  LOGICAL,SAVE :: Initialized = .FALSE.
  LOGICAL,SAVE :: ExtrudedMesh=.False.
  LOGICAL :: Found,GotIt
  LOGICAL :: BoundarySolver
  LOGICAL :: ComputeIceMasks,LimitedSolution,IceFree
  LOGICAL :: stat
  LOGICAL :: SEP

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Flotation'
  CHARACTER(LEN=MAX_NAME_LEN) :: ZbName,ZsName,HName

!------------------------------------------------------------------------------

  Mesh => Model % Mesh

  Params => Solver % Values

  BoundarySolver = ( Solver % ActiveElements(1) > Model % Mesh % NumberOfBulkElements )

!!! get required variables Zb,Zs,H
  ZbName = ListGetString(Params, 'Bottom Surface Name', UnFoundFatal=.TRUE.)
  zbVar => VariableGet( Model % Mesh % Variables, ZbName,UnFoundFatal=.TRUE.)
  
  ZsName = ListGetString(Params, 'Top Surface Name', UnFoundFatal=.TRUE.)
  zsVar => VariableGet( Model % Mesh % Variables, ZsName,UnFoundFatal=.TRUE.)

  HName = ListGetString(Params, 'Thickness Variable Name', UnFoundFatal=.TRUE.)
  HVar => VariableGet( Model % Mesh % Variables, HName, UnFoundFatal=.TRUE.)

!!
!! get optional variables GLMAsk,bedrock
  GLMAsk => VariableGet( Model % Mesh % Variables, 'GroundedMask')
  IF (.NOT.ASSOCIATED(GLMAsk)) THEN
    Message='GroundedMask not found'
    CALL INFO(SolverName,Message,level=5)
  ELSE
    IF ((ParEnv % PEs>1).AND.(.NOT.ASSOCIATED(Solver%Matrix))) &
       CALL FATAL(SolverName,'Solver%Matrix should be associated to update GLMask')
  END IF
  BedVar => VariableGet( Model % Mesh % Variables, 'bedrock')
  IF (.NOT.ASSOCIATED(BedVar)) THEN
     Message='bedrock not found'
     CALL INFO(SolverName,Message,level=5)
  END IF
  ! height above flotation
  HafVar => VariableGet( Model % Mesh % Variables, 'Haf')
  IF (.NOT.ASSOCIATED(HafVar)) THEN
     Message='<Haf> not found; do not compute height above flotation'
     CALL INFO(SolverName,Message,level=5)
  END IF

  !! compute ice area farctions
  ComputeIceMasks = ListGetLogical(Params,"compute ice area fractions",Gotit)
  LimitedSolution = .FALSE.
  IF (ComputeIceMasks) THEN
     FillValue = ListGetConstReal(Params,"Ice free mask values",Gotit)
     IF (.NOT.Gotit) FillValue=0._dp
          
    sftgif => VariableGet( Model % Mesh % Variables, 'sftgif',UnFoundFatal=.TRUE.)
    sftgif  % Values = 0._dp
    IF (sftgif % TYPE /= Variable_on_elements) &
           CALL FATAL(SolverName,"sftgif type should be on_elements")
    sftgrf => VariableGet( Model % Mesh % Variables, 'sftgrf',UnFoundFatal=.TRUE.)
    sftgrf  % Values = 0._dp
    IF (sftgrf % TYPE /= Variable_on_elements) &
           CALL FATAL(SolverName,"sftgrf type should be on_elements")
    sftflf => VariableGet( Model % Mesh % Variables, 'sftflf',UnFoundFatal=.TRUE.)
    sftflf  % Values = 0._dp
    IF (sftflf % TYPE /= Variable_on_elements) &
           CALL FATAL(SolverName,"sftflf type should be on_elements")

    ! change number of IPs for partially grounded elements
    GLnIP=ListGetInteger( Params,'GL integration points number',SEP)

    ! check if we have a limited solution
    ! internal Elmer limiters...
    LimitedSolution=ListCheckPresentAnyBodyForce(CurrentModel, TRIM(HName)//" Lower Limit")
    IF (.NOT.LimitedSolution) &
      LimitedSolution=ListCheckPresentAnyMaterial(CurrentModel, "Min "//TRIM(HName))
  ENDIF

!!! Do some initialisation/allocation
  IF ((.NOT.Initialized) .OR. Mesh % Changed) THEN

    ActiveDirection = ListGetInteger(Params,'Active Coordinate',ExtrudedMesh)
    IF (ExtrudedMesh) THEN
      ! Choose active direction coordinate and set corresponding unit vector
      !---------------------------------------------------------------------

      PSolver => Solver
      CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, BotNodePointer = BotPointer , &
                                TopNodePointer = TopPointer, UpNodePointer = UpPointer)
    END IF

    IF (Initialized) deallocate(Density,GL,MinH,NodalH,Basis)
    
    N=Model % MaxElementNodes
    allocate(Density(N),GL(N),MinH(N),NodalH(N),Basis(N))

    Initialized = .TRUE.
  END IF
!!

 zsea = ListGetCReal( Model % Constants, 'Sea Level', UnFoundFatal=.TRUE. )
 rhow = ListGetCReal( Model % Constants, 'water density', UnFoundFatal=.TRUE. )

 IF (ASSOCIATED(GLMask)) GLMask%Values = -1.0

   IF (BoundarySolver) THEN
     Active = GetNOFBoundaryElements()
   ELSE
     Active = Solver % Mesh % NumberOfBulkElements
   ENDIF

   IF (ASSOCIATED(HafVar)) HafVar%Values = 0._dp

   Do t=1,Active

    IF (BoundarySolver) THEN
      Element => GetBoundaryElement(t,Solver)
    ELSE
      Element => Solver % Mesh % Elements(t)
      CurrentModel % CurrentElement => Element
    ENDIF

    Eindex = Element%ElementIndex
    n = GetElementNOFNodes(Element)
    NodeIndexes => Element % NodeIndexes

    Material => GetMaterial(Element)
    BodyForce => GetBodyForce(Element)

    NodalH(1:n) = HVar%Values(HVar%Perm(NodeIndexes(1:n)))

    IF (ComputeIceMasks) THEN
      kk=sftgif % Perm(Eindex)
      IF (kk>0) sftgif % Values ( kk ) = 1._dp
    END IF

    ! If H is limited check if all element is at lower limit...
    IceFree=.FALSE.
    IF (LimitedSolution) THEN
      MinH = 0._dp
      ! check for limited solution;
      MinH = ListGetConstReal(BodyForce,TRIM(HName)//' Lower Limit', Gotit)
      IF (.NOT.Gotit) &
       MinH = ListGetConstReal(Material,'Min '//TRIM(HName), Gotit)
      IF (.NOT.GotIt) CALL FATAL(SolverName,TRIM(HName)//" not found...but was supposed to be limited")
      IF (ALL((NodalH(1:n)-MinH(1:n)) < EPS)) IceFree=.TRUE.
    END IF

    Density(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,UnFoundFatal=.TRUE.)

    GroundedNode=0
    GL=-1
    Do i=1,n 

       H=NodalH(i)
       rhoi=Density(i)
       zb=zsea-H*rhoi/rhow
      ! if bedrock defined check flotation criterion
       IF(ASSOCIATED(BedVar)) THEN
         bedrock=BedVar%Values(BedVar%Perm(NodeIndexes(i)))
         IF (zb < bedrock) THEN
           zb=bedrock
           GL(i)=1
           GroundedNode=GroundedNode+1
         END IF
         IF (ASSOCIATED(HafVar)) THEN
           Hf=MAX(zsea-bedrock,0._dp)*rhow/rhoi
           HafVar%Values(HafVar%Perm(NodeIndexes(i)))=H-Hf
         END IF
       END IF

       zs=zb+H
       ZbVar%Values(ZbVar%Perm(NodeIndexes(i)))=zb
       ZsVar%Values(ZsVar%Perm(NodeIndexes(i)))=zs

       ! Export Zs on top surface if required
       IF (ExtrudedMesh) THEN
          topnode=TopPointer(NodeIndexes(i))
          ZsVar%Values(ZsVar%Perm(topnode))=zs
       END IF

    End do
    IF (ASSOCIATED(GLMask)) THEN
       IF ((GroundedNode > 0).AND.(GroundedNode < n)) THEN
           WHERE(GL > 0._dp) GL=0._dp
       END IF
       GLMask%Values(GLMask%Perm(NodeIndexes(1:n)))= GL(1:n) * ABS(GLMask%Values(GLMask%Perm(NodeIndexes(1:n))))
    END IF


    IF (ComputeIceMasks) THEN
      IF (GroundedNode == 0) THEN
          !floating element
          IF (sftgrf % Perm(Eindex) > 0 ) &
           sftgrf % Values ( sftgrf % Perm(Eindex)) = 0._dp
          IF (sftflf % Perm(Eindex) > 0 ) &
           sftflf % Values ( sftflf % Perm(Eindex)) = 1._dp
      ELSEIF (GroundedNode < n) THEN
         !partly grounded
         CALL GetElementNodes( ElementNodes, Element )
         IF (SEP .AND. GLnIP > 0 ) THEN
           IP = GaussPoints( Element ,np=GLnIP )
         ELSE
           IP = GaussPoints( Element )
         ENDIF
         area=0._dp
         flarea=0._dp
         grarea=0._dp
         DO ll=1,IP % n
           stat = ElementInfo( Element, ElementNodes, IP % U(ll), IP % V(ll), &
             IP % W(ll),  detJ, Basis )
           rhoi=SUM(Density(1:n)*Basis(1:n))
           H=SUM(NodalH(1:n)*Basis(1:n))
           bedrock=SUM(BedVar%Values(BedVar%Perm(NodeIndexes(1:n)))*Basis(1:n))
           Hf=(zsea-bedrock)*rhow/rhoi
           IParea=detJ*IP % s(ll)
           IF (H < Hf) THEN
             flarea=flarea+IParea
           ELSE
             grarea=grarea+IParea
           ENDIF
           area=area+IParea
         END DO

         IF (sftgrf % Perm(Eindex) > 0 ) &
           sftgrf % Values ( sftgrf % Perm(Eindex)) = grarea/area
         IF (sftflf % Perm(Eindex) > 0 ) &
           sftflf % Values ( sftflf % Perm(Eindex)) = flarea/area
      ELSE
         !grounded
         IF (sftgrf % Perm(Eindex) > 0 ) &     
           sftgrf % Values ( sftgrf % Perm(Eindex)) = 1._dp
         IF (sftflf % Perm(Eindex) > 0 ) &
           sftflf % Values ( sftflf % Perm(Eindex)) = 0._dp
      END IF
      IF (IceFree) THEN
        IF (sftgif % Perm(Eindex) > 0) &
          sftgif % Values ( sftgif % Perm(Eindex)) = FillValue
        IF (sftgrf % Perm(Eindex) > 0) &
          sftgrf % Values ( sftgrf % Perm(Eindex)) = FillValue
        IF (sftflf % Perm(Eindex) > 0) &
          sftflf % Values ( sftflf % Perm(Eindex)) = FillValue
      END IF
    ENDIF
 End Do

 IF (ASSOCIATED(GLMask).AND.( ParEnv % PEs>1 )) CALL ParallelSumVector( Solver % Matrix, GLMask%Values ,1 )

!------------------------------------------------------------------------------
END SUBROUTINE Flotation
!------------------------------------------------------------------------------
