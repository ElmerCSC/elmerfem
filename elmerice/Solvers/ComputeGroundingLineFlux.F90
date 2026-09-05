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
!/*****************************************************************************/
! *  A solver to compute the flux accross the grounding line.
! *    to execute on a bottom boundary of a 3D vertically extruded mesh
! *    where the contact problem for the grounding line is solved.
! * 
! *
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: F. Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: Jan. 2026
! *
! *****************************************************************************
!------------------------------------------------------------------------------
SUBROUTINE ComputeGroundingLineFlux_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  LOGICAL :: Found
  INTEGER :: i

!! default solver variable; the total grounding line flux
  IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable',&
          '-nooutput -global tendligroundf')
  END IF

!! default extrusion direction.  
  IF( .NOT. ListCheckPresent( Solver % Values,'Active Coordinate') ) THEN
      CALL ListAddInteger( Solver % Values,'Active Coordinate',3)
  END IF

!! add elem variable ligroundf  
  i = 1
  DO WHILE(.TRUE.)
    IF(ListCheckPresent(Solver % Values , "Exported Variable "//i2s(i))) THEN
      i=i+1
    ELSE
      EXIT
    END IF
  END DO

  CALL ListAddString(Solver % Values , "Exported Variable "//i2s(i), &
          "-elem ligroundf" )


END SUBROUTINE ComputeGroundingLineFlux_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE ComputeGroundingLineFlux( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE MeshUtils
  USE MeshTransform, ONLY : DetectExtrudedStructure
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t),POINTER :: Psolver
  TYPE(Element_t),POINTER :: BottomElement
  TYPE(Element_t),POINTER :: TopFace,BottomFace,lFace
  TYPE(Element_t),POINTER :: Parent,TopParent
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: EFluxVar,GMask,FlowVar,AreaVar
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  TYPE(Nodes_t),SAVE :: FaceNodes,ElementNodes
  REAL(KIND=dp), ALLOCATABLE,SAVE :: NodalGM(:)
  REAL(KIND=dp), ALLOCATABLE,SAVE :: Basis(:)
  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: cell_area
  REAL(KIND=dp),DIMENSION(3) :: Normal,Flow
  REAL(KIND=dp) :: FillValue
  REAL(KIND=dp) :: U,V,W,detJ
  REAL(KIND=dp) :: gFlux,GlobalFlux,ParVal
  REAL(KIND=dp) :: Norm,PrevNorm,Change
  INTEGER, POINTER,SAVE :: BotPointer(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: DIM,NoDofs
  INTEGER :: n,ngl
  INTEGER :: m
  INTEGER :: EIndex
  INTEGER :: tt,ff,kk,jj
  INTEGER :: nlayers
  INTEGER :: ProjCoord
  INTEGER :: Active
  INTEGER,SAVE :: NumberOfLayers
  LOGICAL,SAVE :: Initialized=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: Caller="ComputeGLFlux"
  LOGICAL :: GotIt
  LOGICAL :: stat
  LOGICAL :: DoProj


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialisation : get paramters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SolverParams => GetSolverParams()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Restricted to 3D vertically extruded
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Mesh => GetMesh()
  IF (Mesh % MeshDim .NE. 3) &
     CALL Fatal(Caller,"Can't handle but 3D mesh, sorry.")

  DIM = CoordinateSystemDimension()
  IF (DIM .NE. 3) &
     CALL Fatal(Caller,"Can't handle but 3D , sorry.")

  IF (.NOT.Initialized) THEN 

     !! Allocation     
     M = CurrentModel % MaxElementNodes
     ALLOCATE(NodalGM(M),Basis(M))

     Active = GetNOFActive(Solver)
     ALLOCATE(cell_area(Active))

     !! Get mesh faces
     IF (.NOT.ASSOCIATED(Mesh % Faces)) &
        CALL FindMeshFaces3D(Mesh)

     !! Get Extruded structure
     Psolver => Solver
     CALL DetectExtrudedStructure( Mesh, PSolver, &
          BotNodePointer = BotPointer, NumberOfLayers = NumberOfLayers)


     !! Compute the element area;
     AreaVar => VariableGet(Mesh%Variables,'bc_area')
     IF (ASSOCIATED(AreaVar)) THEN
       IF (AreaVar % TYPE /= Variable_on_elements) &
         CALL FATAL(Caller,"bc_area type should be on_elements")
     END IF

     !!  we are interested by the projected area; so ProjCoord should be 3..
     ProjCoord = ListGetInteger(Solver % Values,"Active Coordinate",DoProj)
     cell_area=0._dp
     DO tt=1,Solver % NumberOfActiveElements
        BottomElement => GetActiveElement(tt)
        n = GetElementNOFNodes(BottomElement)
        NodeIndexes => BottomElement % NodeIndexes

        CALL GetElementNodes( ElementNodes, BottomElement )
        IF (DoProj) THEN
          SELECT CASE (ProjCoord)
            CASE (1)
              ElementNodes % x(:)=0._dp
            CASE (2)
              ElementNodes % y(:)=0._dp
            CASE (3)
              ElementNodes % z(:)=0._dp
            CASE DEFAULT
              CALL FATAL(Caller, &
                            "Wrong active coordinate"//I2S(ProjCoord))

           END SELECT
         END IF

         IntegStuff = GaussPoints( BottomElement )
         DO kk=1,IntegStuff % n
            U = IntegStuff % u(kk)
            V = IntegStuff % v(kk)
            W = IntegStuff % w(kk)
            stat = ElementInfo(BottomElement,ElementNodes,U,V,W,detJ,Basis)
            cell_area(tt)=cell_area(tt)+detJ*IntegStuff % s(kk)
         END DO
         IF (ASSOCIATED(AreaVar)) THEN
           EIndex=BottomElement % ElementIndex
           IF (ASSOCIATED(AreaVar % Perm)) EIndex=AreaVar % Perm (EIndex)
           IF (EIndex > 0) AreaVar % Values ( EIndex ) = cell_area(tt)
           !PRINT *,AreaVar % Values ( EIndex )
         END IF

     END DO


     Initialized=.TRUE.
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! get required variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT:
  !  - Velocity : assume 'Flow Solution' for now
  FlowVar => VariableGet(Mesh%Variables,'Flow Solution',UnfoundFatal=.TRUE.)
  NoDofs = FlowVar % DOFs

  !  - GroundedMask : assume 'groundedmask' for now
  !    +1=grounded; 0=grounding line (last grounded); -1: Floating
  GMask => VariableGet(Mesh%Variables,'groundedmask',UnfoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(GMask % Perm)) &
     CALL FATAL(Caller,"GMask % Perm not associated")

  ! OUTPUT:
  !  - grounding line flux: 'ligroundf' should be on_element
  EFluxVar => VariableGet(Mesh%Variables,'ligroundf',UnfoundFatal=.TRUE.)
  IF (EFluxVar % TYPE /= Variable_on_elements) &
     CALL FATAL(Caller,"ligroundf type should be on_elements")

  ! optional FillValue for elements not containing the grouning line
  FillValue=ListGetCReal( SolverParams,'FillValue',GotIt)
  IF (GotIt) THEN
    EFluxVar % Values = FillValue
  ELSE
    EFluxVar % Values = 0._dp
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Start the integration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  GlobalFlux=0._dp
  ! go through the active elements, i.e. the bottom boundary elements
  DO tt=1,Solver % NumberOfActiveElements
    BottomElement => GetActiveElement(tt)

    ! restricted to linear elements for now.....
    IF (BottomElement % TYPE % BasisFunctionDegree>1) &
      CALL Fatal(Caller,"Can't handle but linear elements, sorry.")

    n = GetElementNOFNodes(BottomElement)
    NodeIndexes => BottomElement % NodeIndexes

    ! get the nodal groundedmask
    IF(ANY(GMask % Perm ( NodeIndexes(1:n)) == 0)) &
       CALL Fatal(Caller,"GMask % Perm is 0; it should not!!")
    NodalGM(1:n) = GMask % Values ( GMask % Perm ( NodeIndexes(1:n) ) )

    ! the "true" grounding line is within an element
    ! with at list two nodes with mask=0
    ! and one floating nodes (mask=-1)
    ngl=COUNT(NodalGM(1:n) == 0)
    IF (ngl < 2) CYCLE
    IF (.NOT.ANY(NodalGM(1:n).LT.0._dp)) CYCLE


    ! Get the parent element; i.e. the 3D bulk element
    ! parent of the actual bottom face
    IF (.NOT.ASSOCIATED(BottomElement % BoundaryInfo)) &
      CALL Fatal(Caller,"BottomElement % BoundaryInfo not associated !!")

    Parent => BottomElement % BoundaryInfo % Left
    IF (.NOT.ASSOCIATED(Parent)) &
       Parent => BottomElement % BoundaryInfo % Right
    IF (.NOT.ASSOCIATED(Parent)) &
       CALL Fatal(Caller,"Element has no parent !")

    ! Vertical integration 
    BottomFace => BottomElement
    nlayers = 0
    gflux=0._dp
    DO WHILE (ASSOCIATED(Parent))
      nlayers = nlayers + 1

      ! Go through the vertical faces
      ! restricted to vertcially extruded mesh
      !  => face 1 and 2 are the horizontal faces
      !  => others are vertical faces
      DO ff=3,Parent % TYPE % NumberOfFaces
        lFace => Mesh % Faces(Parent % FaceIndexes(ff))
        IF (.NOT.ASSOCIATED(lFace)) &
           CALL Fatal(Caller,"Face not associated")
        n = GetElementNOFNodes(lFace)
        NodeIndexes => lFace % NodeIndexes

        ! A face is a vertical extrusion of an Edge
        ! check the GroundedMask values corresponding to the edge
        IF(ANY(GMask % Perm ( BotPointer(NodeIndexes(1:n))) == 0)) &
           CALL Fatal(Caller,"Bottom GMask % Perm is 0; it should not!!")
        ! Edge is GL if all GM=0
        NodalGM(1:n) = GMask % Values ( GMask % Perm ( BotPointer(NodeIndexes(1:n))))
        IF ( ANY( abs(NodalGM(1:n)) .GT. AEPS ) ) CYCLE

        ! Integrate the flux accross the face
        CALL GetElementNodes(FaceNodes, lFace)

        IntegStuff = GaussPoints( lFace )
        DO kk=1,IntegStuff % n
          U = IntegStuff % u(kk)
          V = IntegStuff % v(kk)
          W = IntegStuff % w(kk)
          stat = ElementInfo(lFace,FaceNodes,U,V,W,detJ,Basis)
          ! normal will point out of the parent
          ! -NormalVector will point to the inside, i.e. towards the floating node
          Normal = 0._dp
          Normal = -NormalVector(lFace,FaceNodes,U,V,.FALSE.,Parent)
          Flow=0._dp
          DO jj=1,DIM
             Flow(jj) = SUM( FlowVar % Values ( NoDofs*(FlowVar % Perm ( NodeIndexes(1:n) ) - 1) + jj ) * Basis(1:n))
          END DO
          gflux=gflux + detJ * IntegStuff % s(kk) * SUM(Normal*Flow)
        END DO
        
      END DO

      !! The Top Face
      !! assume the bottom face is always the 1rst?
      TopFace => Mesh%Faces(Parent % FaceIndexes(2))

      !! sanity check that it's not the actual bottom face...
      IF (BottomFace % ElementIndex == TopFace % ElementIndex) &
         CALL Fatal(Caller,"Bottom and top Faces are the same")

      IF (.NOT.ASSOCIATED(TopFace)) &
         CALL Fatal(Caller,"Face not associated!!")
      IF (.NOT.ASSOCIATED(TopFace % BoundaryInfo)) &
         CALL Fatal(Caller,"Face % BoundaryInfo not associated!!")

      !! get the next parent just above
      TopParent => TopFace % BoundaryInfo % Left
      IF (ASSOCIATED(TopParent)) THEN
        IF (TopParent % ElementIndex == Parent % ElementIndex) &
           TopParent => TopFace % BoundaryInfo % Right
      ELSE 
        TopParent => TopFace % BoundaryInfo % Right
      END IF
      IF (ASSOCIATED(TopParent)) THEN
         Parent => TopParent
         BottomFace => TopFace
      ELSE
         Parent => NULL()
      END IF   

    END DO

    !! sanity check; we should always find the same number of extruded layers..
    IF (nlayers.NE.NumberOfLayers) &
       CALL Fatal(Caller,"Wrong number of layers")

    ! assign the value to the 2D element variable
    EIndex=BottomElement % ElementIndex
    IF (ASSOCIATED(EFluxVar % Perm)) EIndex=EFluxVar % Perm (EIndex)
    IF (EIndex > 0) EFluxVar % Values ( EIndex ) = gflux / cell_area(tt)
    GlobalFlux=GlobalFlux+gflux

  END DO

  !------------------------------------------------------------------------------
  ! Parallel reduction if required
  !-----------------------------------------------------------------------------
  ParVal=ParallelReduction(GlobalFlux,0)
  !------------------------------------------------------------------------------
  ! For consistency checks one may print out a value imitating ComputeChange
  !------------------------------------------------------------------------------
  PrevNorm = Solver % Variable % Norm
  Norm = ABS( ParVal )

  Change = ABS( Norm-PrevNorm )
  IF(Norm + PrevNorm > 0.0) THEN
    Change = Change * 2.0_dp/ (Norm+PrevNorm)
  END IF

  WRITE( Message, '(a,g15.8,g15.8,a)') &
          'SS (ITER=1) (NRM,RELC): (',ParVal, Change,&
          ' ) :: '//TRIM( "GlFlux" )
  CALL Info( 'ComputeChange', Message, Level=3 )

  !------------------------------------------------------------------------------
  ! The solver varaible is a global variable with the total flux
  !------------------------------------------------------------------------------
  Solver % Variable % Values = ParVal
  Solver % Variable % Norm = Norm

END SUBROUTINE ComputeGroundingLineFlux
