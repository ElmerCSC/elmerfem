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
!
!/******************************************************************************
! *
! *  Authors: fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: May 2022
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Module containing utility routines to compute fluxes e.g. at the GL
!--------------------------------------------------------------------------------
      MODULE ComputeFluxUtils
      USE DefUtils

      IMPLICIT NONE
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Compute grounding line flux for 2D applications
!!    the mean flow velocity should be the main solver variable (e.g. ssavelocity) 
!!    required variable groundemask,H              
!!    output variable ligroundf: the element mean GL Flux
!!       
!!    GL flux is the sum of all fluxes from the GL_Edges of a partially
!!     grounded element, i.e. which contain the True GL.
!!    a GL_Edge has 2 Groundedmask==0 (so limited to linear elements for now)
!!    a partially grounded element has at least one GL_Edge and one floating node
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ComputeGLFlux_2D(Solver,FillValue)
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp), OPTIONAL :: FillValue

      Type(Mesh_t), POINTER :: Mesh
      Type(Variable_t), POINTER :: GMask,FlowVar,HVar,EFluxVar
      Type(Element_t), POINTER ::  Element,Edge
      TYPE(Nodes_t),SAVE :: EdgeNodes,ElementNodes
      INTEGER :: tt,ii,jj,kk
      INTEGER :: M,n,nEdges
      INTEGER :: ngl
      INTEGER :: DIM
      INTEGER, POINTER :: NodeIndexes(:)
      LOGICAL, SAVE :: FirstTime=.TRUE.
      CHARACTER(LEN=MAX_NAME_LEN) :: Caller="ComputeGLFlux_2D"

      REAL(KIND=dp), ALLOCATABLE,SAVE :: NodalGM(:)
      REAL(KIND=dp), ALLOCATABLE,SAVE :: Basis(:)
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: U,V,W,detJ
      REAL(KIND=dp),DIMENSION(3) :: Normal,Flow
      REAL(KIND=dp) :: Flux,area,H
      INTEGER :: EIndex
      LOGICAL :: stat

      IF (FirstTime) THEN
         M = CurrentModel % MaxElementNodes
         ALLOCATE(NodalGM(M), Basis(M))
         FirstTime=.FALSE.
      END IF

      !! get required variables
      EFluxVar => VariableGet(Solver%Mesh%Variables,'ligroundf',UnfoundFatal=.TRUE.)
      IF (EFluxVar % TYPE /= Variable_on_elements) &
           CALL FATAL(Caller,"ligroundf type should be on_elements")

      !! min horizontal velocity; SSA velocity in general....
      FlowVar => Solver % Variable
      DIM = FlowVar % DOFs
      IF (DIM /= 2) &
         CALL Fatal(Caller,"Can't handle but 2D flow variable, sorry")
      !! grounded mask
      GMask =>  VariableGet(Solver%Mesh%Variables,'GroundedMask',UnfoundFatal=.TRUE.)

      !! thickness
      HVar =>  VariableGet(Solver%Mesh%Variables,'H',UnfoundFatal=.TRUE.)

      Mesh => Solver % Mesh

      ! Edges will be required....
      CALL FindMeshEdges(Mesh,FindFaces=.FALSE.)

      IF (PRESENT(FillValue)) THEN
        EFluxVar % Values = FillValue
      ELSE
        EFluxVar % Values = 0._dp
      END IF

      DO tt=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(tt)
 
        IF (Element % TYPE % BasisFunctionDegree>1) &
           CALL Fatal(Caller,"Can't handle but linear elements, sorry.")

        n = GetElementNOFNodes(Element)
        NodeIndexes => Element % NodeIndexes

        NodalGM(1:n) = GMask % Values ( GMask % Perm ( NodeIndexes(1:n) ) )

        ! we have an edge GL if at least 2 nodes are GL
        ! and we have a least one floating node
        ngl=COUNT(NodalGM == 0)
        IF (ngl < 2) CYCLE
        IF (.NOT.ANY(NodalGM(1:n).LT.0._dp)) CYCLE

        !! compute total flux from edges
        nEdges = Element % Type % NumberOfEdges
        Flux = 0._dp
        DO ii=1,nEdges
          Edge => Mesh%Edges(Element % EdgeIndexes(ii))
          n = GetElementNOFNodes(Edge)
          NodeIndexes => Edge % NodeIndexes
          NodalGM(1:n) = GMask % Values ( GMask % Perm ( NodeIndexes(1:n) ) )
          ! Edge is GL if all GM=0
          IF ( ANY( abs(NodalGM(1:n)) .GT. AEPS ) ) CYCLE

          CALL GetElementNodes(EdgeNodes, Edge)

          IntegStuff = GaussPoints( Edge )
          DO kk=1,IntegStuff % n
             U = IntegStuff % u(kk)
             V = IntegStuff % v(kk)
             W = IntegStuff % w(kk)
             stat = ElementInfo(Edge,EdgeNodes,U,V,W,detJ,Basis)

             ! normal will poingt ou of the parent
             Normal = -NormalVector(Edge,EdgeNodes,U,V,.FALSE.,Element)

             DO jj=1,DIM
                Flow(jj) = SUM( FlowVar % Values ( DIM*(FlowVar % Perm ( NodeIndexes(1:n) ) - 1) + jj ) * Basis(1:n))
             END DO
             H = SUM (HVar % Values ( HVar % Perm ( NodeIndexes(1:n) ) ) * Basis(1:n))

             Flux = Flux + detJ * IntegStuff % s(kk) * h * SUM(Normal(1:DIM)*Flow(1:DIM))
          END DO
        END DO

        CALL GetElementNodes(ElementNodes, Element)
        ! compute element area
        IntegStuff = GaussPoints( Element )
        area=0._dp
        DO kk=1,IntegStuff % n
          U = IntegStuff % u(kk)
          V = IntegStuff % v(kk)
          W = IntegStuff % w(kk)
          stat = ElementInfo(Element,ElementNodes,U,V,W,detJ,Basis)
          area=area+IntegStuff % s(kk)*detJ
        END DO

        ! mean flux
        EIndex=Element % ElementIndex
        IF (ASSOCIATED(EFluxVar % Perm)) EIndex=EFluxVar % Perm (EIndex)
        IF (EIndex > 0) EFluxVar % Values ( EIndex ) = Flux / area

      END DO
              
      END SUBROUTINE ComputeGLFlux_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Compute Calving Front flux for 2D applications
!!    the mean flow velocity should be the main solver variable (e.g. ssavelocity)
!!    required variable H
!!    output variable : calving_front_flux the parent element mean Calving Front flux
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ComputeCalvingFrontFlux_2D(Solver,FillValue)
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp), OPTIONAL :: FillValue

      Type(Variable_t), POINTER :: FlowVar,HVar,CFluxVar
      TYPE(ValueList_t), POINTER :: BC
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      TYPE(Element_t), POINTER :: Element,Parent
      TYPE(Nodes_t),SAVE :: ElementNodes,ParentNodes

      REAL(KIND=dp), ALLOCATABLE,SAVE :: Basis(:)
      REAL(KIND=dp) :: U,V,W,detJ
      REAL(KIND=dp) :: Normal(3),Flow(3)
      REAL(KIND=dp) :: CalvingFlux,H
      REAL(KIND=dp) :: area

      INTEGER,POINTER :: NodeIndexes(:)
      INTEGER :: t,i,j,k
      INTEGER :: n
      INTEGER :: M
      INTEGER :: EIndex,NofActive
      INTEGER :: DIM

      LOGICAL :: CalvingFront
      LOGICAL :: Found
      LOGICAL :: stat
      LOGICAL,SAVE :: FirstTime=.TRUE.
      LOGICAL,ALLOCATABLE :: VisitedParent(:)

      CHARACTER(LEN=MAX_NAME_LEN) :: Caller="ComputeCalvingFrontFlux_2D"

      IF (FirstTime) THEN
         M = CurrentModel % MaxElementNodes
         ALLOCATE(Basis(M))
         FirstTime=.FALSE.
      END IF


      !! get required variables
      CFluxVar => VariableGet(Solver%Mesh%Variables,'calving_front_flux',UnfoundFatal=.TRUE.)
      IF (CFluxVar % TYPE /= Variable_on_elements) &
       CALL FATAL(Caller,"calving_front_flux type should be on_elements")

      IF (PRESENT(FillValue)) THEN
        CFluxVar % Values = FillValue
      ELSE
        CFluxVar % Values = 0._dp
      END IF

      NofActive = Solver % Mesh % NumberOfBulkElements
      ALLOCATE(VisitedParent(NofActive))
      VisitedParent=.FALSE.

      !! min horizontal velocity; SSA velocity in general....
      FlowVar => Solver % Variable
      DIM = FlowVar % DOFs
      IF (DIM /= 2) &
         CALL Fatal(Caller,"Can't handle but 2D flow variable, sorry")

      !! thickness
      HVar =>  VariableGet(Solver%Mesh%Variables,'H',UnfoundFatal=.TRUE.)


      DO t = 1,GetNOFBoundaryElements()
         Element => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement(Element,Solver) ) CYCLE
         IF ( GetElementFamily(Element) == 1 ) CYCLE

         BC => GetBC()
         IF ( .NOT. ASSOCIATED(BC) ) CYCLE
         CalvingFront=.FALSE.
         CalvingFront=ListGetLogical(BC,'Calving Front', Found)
         IF (.NOT.CalvingFront) CYCLE

         n = GetElementNOFNodes(Element)
         NodeIndexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes )

         IntegStuff = GaussPoints( Element )
         CalvingFlux = 0._dp
         DO i=1,IntegStuff % n
           U = IntegStuff % u(i)
           V = IntegStuff % v(i)
           W = IntegStuff % w(i)
           stat = ElementInfo(Element,ElementNodes,U,V,W,detJ,Basis)

           Normal=0._dp
           Normal = NormalVector( Element,ElementNodes,u,v,.TRUE. )

           Flow=0._dp
           DO j=1,DIM
              Flow(j) = SUM( FlowVar % Values(DIM*(FlowVar % Perm(NodeIndexes(1:n))-1)+j) * Basis(1:n) )
           END DO

           H = SUM (HVar % Values ( HVar % Perm ( NodeIndexes(1:n) ) ) * Basis(1:n))

           CalvingFlux=CalvingFlux+&
                       H*SUM(Normal * Flow)*detJ*IntegStuff % s(i)
         END DO

         ! attribute the flux to the active parent
         Parent => Element % BoundaryInfo % Right
         IF (ASSOCIATED(Parent)) THEN
           ! we can have a passive/active BC=> we attribute the flux to
           ! the active parent
           IF (CheckPassiveElement( Parent )) &
             Parent => Element % BoundaryInfo % Left
         ELSE
           Parent => Element % BoundaryInfo % Left
         END IF
         IF (.NOT.ASSOCIATED(Parent)) CYCLE
         IF (ParEnv % myPe .NE. Parent % partIndex) CYCLE

         CALL GetElementNodes(ParentNodes,Parent)
         ! compute element area
         IntegStuff = GaussPoints( Parent )
         area=0._dp
         DO k=1,IntegStuff % n
            U = IntegStuff % u(k)
            V = IntegStuff % v(k)
            W = IntegStuff % w(k)
            stat = ElementInfo(Parent,ParentNodes,U,V,W,detJ,Basis)
            area=area+IntegStuff % s(k)*detJ
          END DO

        ! mean flux
        EIndex=Parent % ElementIndex
        IF (ASSOCIATED(CFluxVar % Perm)) EIndex=CFluxVar % Perm (EIndex)
        IF (EIndex > 0) THEN
           IF (VisitedParent(Parent % ElementIndex)) THEN
             CFluxVar % Values ( EIndex ) = CFluxVar % Values ( EIndex ) + CalvingFlux / area
           ELSE 
             CFluxVar % Values ( EIndex ) = CalvingFlux / area
             IF (Parent % ElementIndex.GT.NofActive) &
                CALL FATAL(Caller,"Pb with VisitedParent allocated size")
             VisitedParent(Parent % ElementIndex)=.TRUE.
           END IF
        END IF

      END DO

      DEALLOCATE(VisitedParent)

      End

      END MODULE
