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
! *  Authors: F. Gillet-Chaulet (IGE-France)
! *  Web:     http://elmerice.elmerfem.org
! *  Original Date: 04/2019
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conservative FE projection of element variable to a nodal variable
!   Required Solver parameters:
!      Elemental Variable Name = String
!      Nodal Variable Name = String
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE CELL_TO_NODE( Model,Solver,dt,TransientSimulation )
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
      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="CELL_TO_NODE"
      TYPE(ValueList_t), POINTER :: SolverParams
      TYPE(Variable_t),POINTER :: EVar,NVar
      TYPE(Element_t), POINTER :: Element
      TYPE(Nodes_t),SAVE :: ElementNodes
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp),allocatable,SAVE :: weight(:),NodalVar(:)
      REAL(KIND=dp),allocatable,SAVE :: Basis(:), dBasisdx(:,:)
      REAL(KIND=dp) :: U,V,W,SqrtElementMetric
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER, POINTER :: Perm(:)
      INTEGER :: i,t
      INTEGER :: ne,n
      INTEGER :: EIndex
      CHARACTER (len=100) :: EVarName,NVarName
      LOGICAL :: Parallel
      LOGICAL, SAVE :: FirstTime=.TRUE.
      LOGICAL :: stat

      IF (.NOT.ASSOCIATED(Solver % Variable)) &
        CALL FATAL(SolverName,'Solver should have a variable associated')
      Perm => Solver % Variable % Perm

      IF (FirstTime) THEN
        ALLOCATE( weight( Model % NumberOfNodes ),&
                  NodalVar ( Model % NumberOfNodes ),&
                  Basis(Model % MaxElementNodes),&
                  dBasisdx(Model % MaxElementNodes,3))
        FirstTime=.FALSE.
      END IF

! get parameters
      SolverParams => GetSolverParams()

      EVarName = ListGetString(SolverParams,'Elemental Variable Name',UnFoundFatal=.TRUE.)
      NVarName = ListGetString(SolverParams,'Nodal Variable Name',UnFoundFatal=.TRUE.)

! check if this is a paralell run
      Parallel=(ParEnv % PEs > 1)

! get variables
      EVar => VariableGet( Model % Mesh % Variables,TRIM(EVarName),UnFoundFatal=.TRUE.)
      IF(EVar % TYPE /= Variable_on_elements) &
        CALL FATAL(SolverName,'Wrong variable type; use -elem ')

      NVar => VariableGet( Model % Mesh % Variables,TRIM(NVarName),UnFoundFatal=.TRUE.)
      IF(NVar % TYPE /= Variable_on_nodes) &
        CALL FATAL(SolverName,'Wrong variable type; should be nodal')

      NodalVar=0._dp
      weight=0._dp

      ne=GetNOFActive()
      DO t = 1,ne
         Element => GetActiveElement(t)
         EIndex = Element % ElementIndex

         NodeIndexes => Element % NodeIndexes
         IF ( ANY(Perm(NodeIndexes) == 0) ) CYCLE

         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )

         IntegStuff = GaussPoints( Element )

         DO i=1,IntegStuff % n
            U = IntegStuff % u(i)
            V = IntegStuff % v(i)
            W = IntegStuff % w(i)

            stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                        Basis,dBasisdx )

            NodalVar(Perm(NodeIndexes(1:n)))=&
               NodalVar(Perm(NodeIndexes(1:n)))+ &
               EVAR % Values(EVar % Perm(EIndex))*SqrtElementMetric*IntegStuff % s(i) * Basis(1:n)

            weight(Perm(NodeIndexes(1:n)))=weight(Perm(NodeIndexes(1:n)))+&
                      SqrtElementMetric*IntegStuff % s(i) * Basis(1:n)
         END DO
       END DO

       IF (Parallel) THEN
         CALL ParallelSumVector(Solver % Matrix, NodalVar)
         CALL ParallelSumVector(Solver % Matrix, weight)
       END IF

       DO i = 1, Model % NumberOfNodes
        IF ( ABS( weight(Perm(i)) ) > 0.0D0 ) THEN
          NVar % Values (NVar % Perm(i)) = NodalVar(Perm(i)) / weight(Perm(i))
        END IF
      END DO


      END SUBROUTINE CELL_TO_NODE
