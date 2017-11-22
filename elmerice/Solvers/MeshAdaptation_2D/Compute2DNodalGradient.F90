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
! * 
! *****************************************************************************
SUBROUTINE Compute2DNodalGradient( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Nodes_t),SAVE   :: ElementNodes
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: Solution, VSol
  TYPE(GaussIntegrationPoints_t) :: IP

  LOGICAL,SAVE :: AllocationsDone = .FALSE.
  LOGICAL,SAVE :: FEConsistent=.TRUE.
  LOGICAL :: UnfoundFatal=.TRUE.
  LOGICAL :: Found
  LOGICAL :: stat

  LOGICAL,ALLOCATABLE,SAVE :: ActiveNode(:)
  
  INTEGER :: STDOFs

  INTEGER :: i,j,k,l,t
  INTEGER :: n,m
  INTEGER :: istat
  INTEGER :: proc,ierr

  INTEGER, POINTER :: Perm(:), VPerm(:)
  INTEGER, POINTER,SAVE :: NodeIndexes(:)

  REAL(KIND=dp), POINTER :: Values(:), VValues(:)
  REAL(KIND=dp) :: U,V,W
  REAL(KIND=dp) :: dsdx(2)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: grad(:),weight(:)
  REAL(KIND=dp),allocatable,save :: Basis(:),dBasisdx(:,:),ddBasisddx(:,:,:)
  REAL(KIND=dp) :: detJ

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Compute2DNodalGradient'
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: VarName

  TYPE lbuff_t
        INTEGER, ALLOCATABLE :: buff(:)
        REAL(KIND=dp), ALLOCATABLE :: values(:)
  END TYPE lbuff_t
  INTEGER, POINTER :: nlist(:)
  TYPE(lbuff_t), ALLOCATABLE :: n_index(:)
  REAL(KIND=dp), ALLOCATABLE :: nbuff(:)
  INTEGER, ALLOCATABLE :: n_count(:), gbuff(:), n_comp(:)
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

!------------------------------------------------------------------------------
  Solution => Solver % Variable
  Perm   => Solution % Perm
  Values => Solution % Values
  STDOFs =  Solution % DOFs 

     !------------------------------------------------------------------------------
     !    Get variables needed for solution
     !------------------------------------------------------------------------------
      VarName =  ListGetString( Solver % Values,'Variable Name', UnFoundFatal=UnFoundFatal)
      VSol => VariableGet( Solver % Mesh % Variables, trim(VarName), UnFoundFatal=UnFoundFatal )
      VValues => VSol % Values
      VPerm => VSol % Perm

      FEConsistent = GetLogical(Solver % Values,'FE consistent average',Found)
      IF (.NOT.Found) FEConsistent=.True.

     !--------------------------------------------------------------
     !Allocate some permanent storage, this is done first time only:
     !--------------------------------------------------------------
      IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
        N = Model % MaxElementNodes
        M = Model % Mesh % NumberOfNodes
        IF (AllocationsDone) DEALLOCATE(ElementNodes % x, &
                       ElementNodes % y, ElementNodes % z, &
                       Basis,dBasisdx,ddBasisddx,&
                       grad, &
                       ActiveNode, &
                       weight)

        ALLOCATE(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N), &
             Basis(N),dBasisdx(N,3),ddBasisddx(N,3,3),&
              grad(STDOFs*M), &
              ActiveNode(M), &
              weight(M),&
           STAT=istat )
        IF ( istat /= 0 ) THEN
          CALL Fatal( SolverName, 'Memory allocation error.' )
        END IF
        AllocationsDone = .TRUE.
        CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
      END IF

      ActiveNode=.False.
      grad=0._dp
      weight=0.0

   !!!!!! Compute gradient elementwise
      Do t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes

        ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (STDOFs == 1) THEN !1D 
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (STDOFs == 2) THEN !2D 
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute slope with DOFs=',&
                  STDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message)
        END IF

        Do i=1,n
           k=NodeIndexes(i)
           ActiveNode(k)=.True.

           IF (FEConsistent) THEN
              IP = GaussPoints(Element)
              DO j=1,IP % n
                stat = ElementInfo(Element, ElementNodes, IP % U(j), &
                                           IP % v(j), IP % W(j), detJ, &
                                           Basis,dBasisdx, ddBasisddx,.FALSE. )

                dsdx(1)=SUM( VValues(VPerm(NodeIndexes(1:n))) * dBasisdx(1:n,1) )
                IF (STDOFs == 2) &
                   dsdx(2)=SUM( VValues(VPerm(NodeIndexes(1:n))) * dBasisdx(1:n,2) )

                weight(k)=weight(k)+IP%s(j)*detJ*Basis(i)

                grad(STDOFs*(k-1)+1) = grad(STDOFs*(k-1)+1) + &
                                       dsdx(1)*IP%s(j)*detJ*Basis(i)
                IF (STDOFs == 2) grad(STDOFs*(k-1)+2) =grad(STDOFs*(k-1)+2) + &
                                                dsdx(2)*IP%s(j)*detJ*Basis(i)
              END DO !IPs
           ELSE

              U=Element % TYPE % NodeU(i)
              V=Element % TYPE % NodeV(i)
              W=0.0

              stat = ElementInfo( Element, ElementNodes, U, V, &
                     W,  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

              weight(k)=weight(k)+1.0_dp


             grad(STDOFs*(k-1)+1) = grad(STDOFs*(k-1)+1) + &
                     SUM( VValues(VPerm(NodeIndexes(1:n))) * dBasisdx(1:n,1) )

             IF (STDOFs == 2) grad(STDOFs*(k-1)+2) =grad(STDOFs*(k-1)+2) + &
                            SUM( VValues(VPerm(NodeIndexes(1:n))) * dBasisdx(1:n,2) )
           END IF
        End do !elements nodes

      End do !elements

      IF (ParEnv % PEs>1 ) THEN  !! if parallel need to sum values at interfaces
                                     !! here is a copy of what is done in
                                     !SolverUtils.src to average boundary normals
         ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
         n_count = 0

         DO i=1,Solver%Mesh % NumberOfNodes
            IF (.NOT.ActiveNode(i)) CYCLE
            IF (.NOT.Solver % Mesh % ParallelInfo % INTERFACE(i) ) CYCLE
  
            nlist => Solver%Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
               k = nlist(j)+1
               IF ( k-1 == ParEnv % myPE ) CYCLE
               n_count(k) = n_count(k)+1
            END DO
         END DO
         DO i=1,ParEnv % PEs
            IF ( n_count(i)>0 ) &
               ALLOCATE( n_index(i) % buff(n_count(i)), &
                         n_index(i) % values((STDOFs+1)*n_count(i)) )
         END DO

         n_count = 0
         DO i=1,Model % NumberOfNodes
            IF (.NOT.ActiveNode(i)) CYCLE
            IF (.NOT.Solver % Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

            nlist =>Solver% Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
               k = nlist(j)+1
               IF ( k-1 == ParEnv % myPE ) CYCLE
               n_count(k) = n_count(k)+1
               n_index(k) % buff(n_count(k)) = Solver%Mesh % Parallelinfo % &
                  GlobalDOFs(i)
               n_index(k) % values((STDOFs+1)*(n_count(k)-1)+1)=grad(STDOFs*(i-1)+1)
               if (STDOFS.EQ.2) &
                  n_index(k) % Values((STDOFs+1)*(n_count(k)-1)+2)=grad(STDOFs*(i-1)+2)
               n_index(k) % values((STDOFs+1)*(n_count(k)-1)+STDOFs+1)=weight(i)
            END DO
         END DO

         DO i=1,ParEnv % PEs
            IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
              CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                  900, MPI_COMM_WORLD, ierr )
               IF ( n_count(i)>0 ) THEN
                 CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, MPI_COMM_WORLD, ierr )
                 CALL MPI_BSEND( n_index(i) % values, (STDOFs+1)*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, MPI_COMM_WORLD, ierr )
               END IF
            END IF
        END DO
        DO i=1,ParEnv % PEs
           IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % values)

           IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
              CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                    900, MPI_COMM_WORLD, status, ierr )
              IF ( n>0 ) THEN
                 proc = status(MPI_SOURCE)
                 ALLOCATE( gbuff(n), nbuff((STDOFs+1)*n) )
                 CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                     901, MPI_COMM_WORLD, status, ierr )

                 CALL MPI_RECV( nbuff, (STDOFs+1)*n, MPI_DOUBLE_PRECISION, proc, &
                     902, MPI_COMM_WORLD, status, ierr )

                 DO j=1,n
                    k = SearchNodeL(Solver% Mesh % ParallelInfo, gbuff(j), Solver%Mesh % NumberOfNodes )

                    IF ( k>0 ) THEN
                      !n_comp(k) = n_comp(k)+1
                      ! IF ( l>0 ) THEN
                      grad(STDOFs*(k-1)+1)=grad(STDOFs*(k-1)+1)+nbuff((STDOFs+1)*(j-1)+1)
                      IF (STDOFs.EQ.2) &
                         grad(STDOFs*(k-1)+2)=grad(STDOFs*(k-1)+2)+nbuff((STDOFs+1)*(j-1)+2)
                      weight(k)=weight(k)+nbuff((STDOFs+1)*(j-1)+STDOFs+1)
                      ! END IF
                   END IF
                 END DO
                 DEALLOCATE(gbuff, nbuff)
               END IF
            END IF
         END DO
         DEALLOCATE( n_index, n_count )
           !DEALLOCATE(n_comp)

       END IF !end do parallel reduction

       Do t=1,Solver % Mesh % NumberOfNodes
          IF (.NOT.ActiveNode(t)) CYCLE 
          Values(STDOFs*(Perm(t)-1)+1) = grad(STDOFs*(t-1)+1) / weight(t)
          IF (STDOFs == 2) Values(STDOFs*(Perm(t)-1)+2)=grad(STDOFs*(t-1)+2) / weight(t)
       END DO

!!! As there is no equation solved this will allow to interpolate the solution
!from meshes to meshes if VariableGet is used
          CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
                    Solver % Variable % Name )

!------------------------------------------------------------------------------
END SUBROUTINE Compute2DNodalGradient
!------------------------------------------------------------------------------

