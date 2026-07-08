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
! *****************************************************************************
! * COUPLING WITH PDAF (Parallel Data Assimilation Framework)
! *  https://pdaf.awi.de/trac/wiki/WikiStart
! *  Code adapted from PDAF templates and examples
! *****************************************************************************
! * Routines used to initialised the local domain analyses in local filters
! * LSEIK/LETKF/LESTKF/LNETF/LKNETF and the 3DEnVar and hybrid 3DVAR
! * 
! * init_n_domains_pdaf: 
! *  - Initialise the  number of local analysis domains on the local PE
! * init_dim_l_pdaf:
! *  - Initialise the dimension of the local model  state on the current analysis
! *   domain.
! ***************************************************************************** 
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)
  USE DefUtils
  USE mod_statevector_pdaf, ONLY : n_fields,sfields,DomainToNode

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains

  TYPE(Mesh_t),POINTER :: Mesh
  INTEGER :: node,j
  LOGICAL :: IsDomain

! ************************************
! *** Initialize number of domains ***
! ************************************
! * A domain is a node where et least one state variable is Active
! * Initialise the DomainToNode permutation table
  Mesh => GetMesh()

  n_domains_p = 0
  DO node=1,Mesh%NumberOfNodes

    IsDomain=.FALSE.
    DO j=1,n_fields
      IsDomain=IsDomain.OR.(sfields(j)%NodeToState(node) > 0)
    END DO  
    IF (.NOT.IsDomain) CYCLE

    n_domains_p = n_domains_p + 1
    DomainToNode(n_domains_p) = node
  END DO

END SUBROUTINE init_n_domains_pdaf

SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE DefUtils
  USE PDAF, &                  ! Routine to provide local indices to PDAF
       ONLY: PDAFlocal_set_indices
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l
  USE mod_statevector_pdaf, ONLY : n_fields,sfields,DomainToNode,mfields,n_mfields,&
                                   ReducedMask
  USE mod_io

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  CHARACTER(*), PARAMETER :: Caller="init_dim_l_pdaf"
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
  INTEGER :: node,j,cnt
  INTEGER :: mask
  TYPE(Mesh_t),POINTER :: Mesh


  Mesh => GetMesh()
! ****************************************
! *** Initialize local state dimension ***
! ****************************************
 !* get mesh node corresponding to the current domain
  node=DomainToNode(domain_p)

 !* Count number of State variables active at the given node
 !* If the state variable has a max, exclude the variable from the analysis if Max(max) < 0.
  dim_l = 0
  DO j=1,n_fields
   IF (sfields(j)%NodeToState(node) > 0) THEN
      mask=sfields(j)%mask
      IF (mask > 0) THEN
         IF (mask > n_mfields) &
            CALL PDAF_FATAL(Caller,'maks > n_mask for variable '//TRIM(sfields(j)%name))
         IF (.NOT.(ReducedMask(mfields(mask)%NodeToState(node) + mfields(mask)%off) < 0._dp)) &
            dim_l=dim_l + 1
      ELSE
         dim_l=dim_l + 1
      END IF        
   END IF
  END DO

! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************
  !* Global coordinates of local analysis domain
  coords_l(1)=Mesh % Nodes % x(node)
  coords_l(2)=Mesh % Nodes % y(node)
  coords_l(3)=Mesh % Nodes % z(node)


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************
  ! Allocate array
  ALLOCATE(id_lstate_in_pstate(dim_l))
  cnt = 0
   DO j=1,n_fields
     IF (sfields(j)%NodeToState(node) > 0) THEN
       mask=sfields(j)%mask
       IF (mask > 0) THEN
         IF (.NOT.(ReducedMask(mfields(mask)%NodeToState(node) + mfields(mask)%off) < 0._dp)) THEN
           cnt=cnt+1
           id_lstate_in_pstate(cnt)=sfields(j)%off + sfields(j)%NodeToState(node)
         END IF
       ELSE
          cnt=cnt+1
          id_lstate_in_pstate(cnt)=sfields(j)%off + sfields(j)%NodeToState(node)        
       ENDIF
     END IF
   END DO

  ! Provide the index vector to PDAF
  CALL PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  DEALLOCATE(id_lstate_in_pstate)

END SUBROUTINE init_dim_l_pdaf

