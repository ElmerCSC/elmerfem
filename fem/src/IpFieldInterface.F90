!------------------------------------------------------------------------------
!> The explicit interface to Ip2DgFieldInElement, which is an external
!> procedure in InterpolateMeshToMesh.F90.
!>
!> That routine has to stay external: it needs Integration and
!> ElementDescription, while one of its callers is Lists.F90, which sits below
!> both in the module hierarchy and hence cannot USE a module that hosts it.
!> What it does not need is a private copy of the interface in each caller.
!> Five such copies existed, and the sixth copy of this kind -- the one for
!> HierarchicPToLagrange -- is what silently drifted from its definition and
!> killed high-order p-element output. One copy, here, keeps them honest.
!------------------------------------------------------------------------------
MODULE IpFieldInterface

  USE Types

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE Ip2DgFieldInElement( Mesh, Element, nip, fip, ndg, fdg )
      USE Types
      IMPLICIT NONE
      TYPE(Mesh_t) :: Mesh
      TYPE(Element_t), TARGET :: Element
      INTEGER :: nip, ndg
      REAL(KIND=dp) :: fip(:), fdg(:)
    END SUBROUTINE Ip2DgFieldInElement
  END INTERFACE

!------------------------------------------------------------------------------
END MODULE IpFieldInterface
!------------------------------------------------------------------------------
