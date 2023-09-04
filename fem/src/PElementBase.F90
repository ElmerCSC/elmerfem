!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Mikko Byckling, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20 Aug 2004
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!>  Module defining p element basis functions. All p basis (and related) 
!>  functions are defined here, as well as few helper routines for determining 
!>  if an element is p element. For mappings related to p elements see 
!>  module PElementMaps.
!-----------------------------------------------------------------------------

MODULE PElementBase
  USE PElementMaps
  USE Messages
  USE Types, ONLY : dp, Element_t, Mesh_t
  IMPLICIT NONE

  CONTAINS

    ! 1D ELEMENTS

!------------------------------------------------------------------------------
    FUNCTION LineNodalPBasis(node, u) RESULT(value)
!------------------------------------------------------------------------------
!
!  DESCRIPTION:
!     Nodal basis of line element at point (u)
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of nodal function to calculate, node = {1,2}
!
!    REAL(KIND=dp) :: u
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of lines nodal function i at point u, i.e.
!       value = N_i(u)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u
      REAL (KIND=dp) :: value

      SELECT CASE(node)
      CASE (1)
         value = (1d0-u)/2d0
      CASE (2)
         value = (1d0+u)/2d0
      CASE DEFAULT
         value = 0.0_dp
         CALL Fatal('PElementBase::LineNodalPBasis', 'Unknown node for line')
      END SELECT
    END FUNCTION LineNodalPBasis


    ! As previous except obtain all values at once.
    SUBROUTINE LineNodalPBasisAll(u, phi) 

      IMPLICIT NONE

      REAL (KIND=dp), INTENT(IN) :: u
      REAL (KIND=dp) :: phi(:)
      REAL(Kind=dp), PARAMETER :: c = 1.0_dp/2.0_dp      
      INTEGER, PARAMETER :: usgn(2) = [-1,1]
      
      phi(1:2) = c*(1+usgn*u)
      
    END SUBROUTINE LineNodalPBasisAll


    
!------------------------------------------------------------------------------
!>     Derivative of line elements nodal basis at point (u).
!------------------------------------------------------------------------------
    FUNCTION dLineNodalPBasis(node, u) RESULT(grad)
!------------------------------------------------------------------------------
!
!    INTEGER :: node
!      INPUT: number of nodal function to calculate, node = {1,2}
!
!    REAL(KIND=dp) :: u
!      INPUT: point at which to evaluate function derivative
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of derivative of lines nodal function i at point u, i.e.
!       value = dN_i(u)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u
      REAL (KIND=dp) :: grad

      SELECT CASE(node)
      CASE (1)
         grad = -1d0/2
      CASE (2)
         grad = 1d0/2
      CASE DEFAULT
         grad = 0.0_dp
         CALL Fatal('PElementBase::dLineNodalPBasis', 'Unknown node for line')
      END SELECT
    END FUNCTION dLineNodalPBasis


    ! As previous except obtain all values at once.
    SUBROUTINE dLineNodalPBasisAll(u, gradphi) 

      IMPLICIT NONE

      REAL (KIND=dp), INTENT(IN) :: u
      REAL (KIND=dp) :: gradphi(:,:)
      REAL(Kind=dp), PARAMETER :: c = 1.0_dp/2.0_dp      
      INTEGER, PARAMETER :: usgn(2) = [-1,1]
      
      gradphi(1:2,1) = c*(usgn)
      
    END SUBROUTINE dLineNodalPBasisAll

    
!------------------------------------------------------------------------------
!>     2nd derivative of line elements nodal basis at point (u).
!------------------------------------------------------------------------------
    FUNCTION ddLineNodalPBasis(node, u) RESULT(grad)
!------------------------------------------------------------------------------
!
!    INTEGER :: node
!      INPUT: number of nodal function to calculate, node = {1,2}
!
!    REAL(KIND=dp) :: u
!      INPUT: point at which to evaluate function derivative
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of derivative of lines nodal function i at point u, i.e.
!       value = dN_i(u)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u
      REAL (KIND=dp) :: grad

      grad = 0
    END FUNCTION ddLineNodalPBasis


!------------------------------------------------------------------------------
!>     Bubble function i of line element.
!------------------------------------------------------------------------------
    FUNCTION LineBubblePBasis(i, u, invertEdge) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i
!      INPUT: index of bubble function to calculate, i = {2,...}
!
!    REAL(KIND=dp) :: u
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert this edge or not. Used in calculation of edge 
!      boundary values for 2d element. If direction of bubble function is 
!      inverted parameter of phi function is varied from [1,-1] in stead of 
!      the usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of lines bubble function i at point u, i.e.
!       value = N_i^(0)(u)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: u
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp) :: phiPar, value
      LOGICAL :: invert
      
      ! Check if line basis has been inverted (not by default)
      invert = .FALSE.
      IF (PRESENT( invertEdge )) invert = invertEdge

      phiPar = u
      IF (invert) phiPar = -phiPar

      value = Phi(i,phipar)
    END FUNCTION LineBubblePBasis



!------------------------------------------------------------------------------
!>     Derivative of bubble function i of line element.
!------------------------------------------------------------------------------
    FUNCTION dLineBubblePBasis(i, u, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i
!      INPUT: index of bubble function derivative to calculate, , i = {1,...}
!
!    REAL(KIND=dp) :: u
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert this edge or not. Used in calculation of edge 
!      boundary values for 2d element. If direction of bubble function is 
!      inverted parameter of phi function is varied from [1,-1] in stead of 
!      the usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of derivative of lines bubble function i at point u, i.e.
!       value = dN_i^(0)(u)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: u
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp) :: phiPar, grad 
      LOGICAL :: invert
      
      ! Check if line basis has been inverted (not by default)
      invert = .FALSE.
      IF (PRESENT( invertEdge )) invert = invertEdge
      
      phiPar = u
      IF (invert) phiPar = -phiPar

      grad = dPhi(i,phiPar)
    END FUNCTION dLineBubblePBasis



!------------------------------------------------------------------------------
!>     2nd derivative of bubble function i of line element.
!------------------------------------------------------------------------------
    FUNCTION ddLineBubblePBasis(i, u, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i
!      INPUT: index of bubble function derivative to calculate, , i = {1,...}
!
!    REAL(KIND=dp) :: u
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert this edge or not. Used in calculation of edge 
!      boundary values for 2d element. If direction of bubble function is 
!      inverted parameter of phi function is varied from [1,-1] in stead of 
!      the usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of 2nd derivative of lines bubble function i at point u
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: u
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp) :: phiPar, grad 
      LOGICAL :: invert
      
      ! Check if line basis has been inverted (not by default)
      invert = .FALSE.
      IF (PRESENT( invertEdge )) invert = invertEdge
      
      phiPar = u
      IF (invert) phiPar = -phiPar
      
      grad = ddPhi(i,phiPar)
    END FUNCTION ddLineBubblePBasis


    ! 2D ELEMENTS
!------------------------------------------------------------------------------
!>     Quadrilateral nodal basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION QuadNodalPBasis(node, u, v) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of quadrilaterals nodal function to calculate
!        node = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of quadrilaterals nodal function at point (u,v), i.e.
!       value = N_i(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: value
      
      value = 0
      ! By local edge, calculate value of nodal function
      SELECT CASE(node)
      CASE (1)
         value = (1-u)*(1-v)/4
      CASE (2)
         value = (1+u)*(1-v)/4
      CASE (3)
         value = (1+u)*(1+v)/4
      CASE (4)
         value = (1-u)*(1+v)/4
      CASE DEFAULT
         CALL Fatal('PElementBase::QuadNodalPBasis', 'Unknown node for quadrilateral')
      END SELECT
    END FUNCTION QuadNodalPBasis


    ! As previous except obtain all values at once.
    SUBROUTINE QuadNodalPBasisAll(u, v, phi) 

      IMPLICIT NONE

      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: phi(:)
      REAL(Kind=dp), PARAMETER :: c = 1.0_dp/4.0_dp      
      INTEGER, PARAMETER :: usgn(4) = [-1,1,1,-1]
      INTEGER, PARAMETER :: vsgn(4) = [-1,-1,1,1]
      
      phi(1:4) = c*(1+usgn*u)*(1+vsgn*v)
      
    END SUBROUTINE QuadNodalPBasisAll
       

!------------------------------------------------------------------------------
!>     Gradient of quadrilateral nodal basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION dQuadNodalPBasis(node, u, v) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of derivative of quadrilateral s nodal function to 
!        calculate, node = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function derivative
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals nodal function at point (u,v),
!       i.e. value = dN_i(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp), DIMENSION(2) :: grad
      
      grad = 0
      ! By local edge, calculate value of nodal function
      SELECT CASE(node)
      CASE (1)
         grad(1) = -(1-v)/4
         grad(2) = -(1-u)/4
      CASE (2)
         grad(1) =  (1-v)/4
         grad(2) = -(1+u)/4
      CASE (3)
         grad(1) =  (1+v)/4
         grad(2) =  (1+u)/4
      CASE (4)
         grad(1) = -(1+v)/4
         grad(2) =  (1-u)/4
      CASE DEFAULT
         CALL Fatal('PElementBase::dQuadNodalPBasis', 'Unknown node for quadrilateral')
      END SELECT
    END FUNCTION dQuadNodalPBasis
!------------------------------------------------------------------------------


!>     2nd derivatives of quadrilateral nodal basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION ddQuadNodalPBasis(node, u, v) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of derivative of quadrilateral s nodal function to 
!        calculate, node = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function derivative
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of quadrilaterals nodal function at point (u,v),
!       i.e. value = dN_i(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp), DIMENSION(2,2) :: grad
      
      ! By local edge, calculate value of nodal function
      grad = 0
      SELECT CASE(node)
      CASE (1,3)
         grad(1,2) =  1
         grad(2,1) =  1
      CASE (2,4)
         grad(1,2) = -1
         grad(2,1) = -1
      CASE DEFAULT
         CALL Fatal('PElementBase::ddQuadNodalPBasis', 'Unknown node for quadrilateral')
      END SELECT
      grad = grad/4
    END FUNCTION ddQuadNodalPBasis


    ! As previous except obtain all values at once 
    SUBROUTINE dQuadNodalPBasisAll(u, v, gradphi) 
      IMPLICIT NONE

      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: gradphi(:,:)
      
      INTEGER, PARAMETER :: usgn(4) = [-1,1,1,-1]
      INTEGER, PARAMETER :: vsgn(4) = [-1,-1,1,1]
      REAL(Kind=dp), PARAMETER :: c = 1.0_dp/4.0_dp      
     
      gradphi(1:4,1) = c*(usgn)*(1+vsgn*v)
      gradphi(1:4,2) = c*(1+usgn*u)*(vsgn)
      
    END SUBROUTINE dQuadNodalPBasisAll



!  --- start serendipity quad ---
    
!------------------------------------------------------------------------------
!>     Quadrilateral edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION SD_QuadEdgePBasis(edge, i, u, v, invertEdge) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of quadrilaterals edge function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of quadrilaterals edge function i at point (u,v), i.e.
!       value = N_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN) :: edge, i
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: value
      LOGICAL :: invert 
      
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      value = 0
      ! By local edge, calculate value of edge function
      SELECT CASE(edge)
      CASE (1)
         IF (.NOT. invert) THEN
            value = 1d0/2*(1-v)*Phi(i,u)
         ELSE
            value = 1d0/2*(1-v)*Phi(i,-u)
         END IF
      CASE (2)
         IF (.NOT. invert) THEN
            value = 1d0/2*(1+u)*Phi(i,v)
         ELSE
            value = 1d0/2*(1+u)*Phi(i,-v)
         END IF
      CASE (3)
         IF (.NOT. invert) THEN
            value = 1d0/2*(1+v)*Phi(i,u)
         ELSE
            value = 1d0/2*(1+v)*Phi(i,-u)
         END IF
      CASE (4)
         IF (.NOT. invert) THEN
            value = 1d0/2*(1-u)*Phi(i,v)
         ELSE 
            value = 1d0/2*(1-u)*Phi(i,-v)
         END IF
      CASE DEFAULT
         CALL Fatal('PElementBase::QuadEdgePBasis', 'Unknown edge for quadrilateral')
      END SELECT
    END FUNCTION SD_QuadEdgePBasis

!------------------------------------------------------------------------------
!>     2nd derivatives of quadrilateral edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION SD_ddQuadEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of quadrilaterals edge function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of quadrilaterals edge function i at point (u,v), i.e.
!       grad = dN_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN) :: edge, i
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp), DIMENSION(2,2) :: grad
      LOGICAL :: invert

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      grad = 0
      ! By local edge, calculate value of edge function
      SELECT CASE(edge)
      CASE (1)
         IF (.NOT. invert) THEN
            grad(1,1) =  (1-v)*ddPhi(i,u)
            grad(1,2) = -dPhi(i,u)
            grad(2,2) = 0
         ELSE 
            grad(1,1) =  (1-v)*ddPhi(i,-u)
            grad(1,2) =  dPhi(i,-u)
            grad(2,2) = 0
         END IF
      CASE (2)
         IF (.NOT. invert) THEN
            grad(1,1) = 0
            grad(1,2) = dPhi(i,v)
            grad(2,2) = (1+u)*ddPhi(i,v)
         ELSE 
            grad(1,1) = 0
            grad(1,2) =-dPhi(i,-v)
            grad(2,2) = (1+u)*ddPhi(i,-v)
         END IF
      CASE (3)
         IF (.NOT. invert) THEN
            grad(1,1) = (1+v)*ddPhi(i,u)
            grad(1,2) = dPhi(i,u)
            grad(2,2) = 0
         ELSE
            grad(1,1) = (1+v)*ddPhi(i,-u)
            grad(1,2) =-dPhi(i,-u)
            grad(2,2) = 0
         END IF
      CASE (4)
         IF (.NOT. invert) THEN
            grad(1,1) = 0
            grad(1,2) =-dPhi(i,v)
            grad(2,2) = (1-u)*ddPhi(i,v)
         ELSE
            grad(1,1) = 0
            grad(1,2) = dPhi(i,-v)
            grad(2,2) = (1-u)*ddPhi(i,-v)
         END IF
      CASE DEFAULT
         CALL Fatal('PElementBase::ddQuadEdgePBasis', 'Unknown edge for quadrilateral')
      END SELECT
      grad = grad/2
      grad(2,1) = grad(1,2)
    END FUNCTION SD_ddQuadEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of quadrilateral edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION SD_dQuadEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of quadrilaterals edge function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals edge function i at point (u,v), i.e.
!       grad = dN_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN) :: edge, i
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp), DIMENSION(2) :: grad
      LOGICAL :: invert

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      grad = 0
      ! By local edge, calculate value of edge function
      SELECT CASE(edge)
      CASE (1)
         IF (.NOT. invert) THEN
            grad(1) = 1d0/2*(1-v)*dPhi(i,u)
            grad(2) = -1d0/2*Phi(i,u)
         ELSE 
            grad(1) = -1d0/2*(1-v)*dPhi(i,-u)
            grad(2) = -1d0/2*Phi(i,-u)
         END IF
      CASE (2)
         IF (.NOT. invert) THEN
            grad(1) = 1d0/2*Phi(i,v)
            grad(2) = 1d0/2*(u+1)*dPhi(i,v)
         ELSE 
            grad(1) = 1d0/2*Phi(i,-v)
            grad(2) = -1d0/2*(u+1)*dPhi(i,-v)
         END IF
      CASE (3)
         IF (.NOT. invert) THEN
            grad(1) = 1d0/2*(1+v)*dPhi(i,u)
            grad(2) = 1d0/2*Phi(i,u)
         ELSE
            grad(1) = -1d0/2*(1+v)*dPhi(i,-u)
            grad(2) = 1d0/2*Phi(i,-u)
         END IF
      CASE (4)
         IF (.NOT. invert) THEN
            grad(1) = -1d0/2*Phi(i,v)
            grad(2) = 1d0/2*(1-u)*dPhi(i,v)
         ELSE
            grad(1) = -1d0/2*Phi(i,-v)
            grad(2) = -1d0/2*(1-u)*dPhi(i,-v)
         END IF
      CASE DEFAULT
         CALL Fatal('PElementBase::dQuadEdgePBasis', 'Unknown edge for quadrilateral')
      END SELECT
    END FUNCTION SD_dQuadEdgePBasis


!------------------------------------------------------------------------------
!>     Quadrilateral bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION SD_QuadBubblePBasis(i,j,u,v,localNumbers) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of edge function, (i,j) = {(2,2),(2,3),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of quarilateral to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of quadrilateral faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of quadrilaterals bubble function (i,j) at point (u,v), 
!       i.e. value = N_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(4)
      REAL (KIND=dp) :: La, Lb, Lc, value

      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
         value = Phi(i,u)*Phi(j,v)
         RETURN
      END IF
      
      ! Numbering present, so use it
      La = QuadL(localNumbers(1),u,v)
      Lb = QuadL(localNumbers(2),u,v)
      Lc = QuadL(localNumbers(4),u,v)

      ! Calculate value of function from general form
      value = Phi(i,Lb-La)*Phi(j,Lc-La)
    END FUNCTION SD_QuadBubblePBasis


!------------------------------------------------------------------------------
!>     Gradient of quadrilateral bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION SD_dQuadBubblePBasis(i,j,u,v, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of edge function, (i,j) = {(2,2),(2,3),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of quarilateral to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of quadrilateral faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals bubble function (i,j) at point (u,v), 
!       i.e. grad = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(4)
      REAL(Kind=dp) :: La, Lb, Lc
      REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, grad
      
      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
         grad(1) = dPhi(i,u)*Phi(j,v)
         grad(2) = Phi(i,u)*dPhi(j,v)
         RETURN
      END IF

      ! Numbering present, so use it
      La = QuadL(localNumbers(1),u,v)
      Lb = QuadL(localNumbers(2),u,v)
      Lc = QuadL(localNumbers(4),u,v)
      dLa = dQuadL(localNumbers(1),u,v)
      dLb = dQuadL(localNumbers(2),u,v)
      dLc = dQuadL(localNumbers(4),u,v)

      grad = dPhi(i,Lb-La)*(dLb-dLa)*Phi(j,Lc-La) + &
           Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc-dLa)
    END FUNCTION SD_dQuadBubblePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of quadrilateral bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION SD_ddQuadBubblePBasis(i,j,u,v, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of edge function, (i,j) = {(2,2),(2,3),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of quarilateral to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of quadrilateral faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of quadrilaterals bubble function (i,j) at point (u,v), 
!       i.e. grad = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(4)
      REAL(Kind=dp) :: La, Lb, Lc
      REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc
      REAL(Kind=dp), DIMENSION(2,2) :: grad
      
      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
         grad(1,1) = ddPhi(i,u)*Phi(j,v)
         grad(1,2) = dPhi(i,u)*dPhi(j,v)
         grad(2,1) = dPhi(i,u)*dPhi(j,v)
         grad(2,2) = Phi(i,u)*ddPhi(j,v)
         RETURN
      END IF

      ! Numbering present, so use it
      La = QuadL(localNumbers(1),u,v)
      Lb = QuadL(localNumbers(2),u,v)
      Lc = QuadL(localNumbers(4),u,v)

      dLa = dQuadL(localNumbers(1),u,v)
      dLb = dQuadL(localNumbers(2),u,v)
      dLc = dQuadL(localNumbers(4),u,v)

      grad(1,1) = ddPhi(i,Lb-La)*(dLb(1)-dLa(1))**2*Phi(j,Lc-La)
      grad(1,1) = grad(1,1) + dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,1) = grad(1,1) + dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLC(1)-dLa(1))
      grad(1,1) = grad(1,1) + Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLC(1)-dLa(1))**2

      grad(1,2) = ddPhi(i,Lb-La)*(dLb(1)-dLa(1))*(dLb(2)-dLa(2))*Phi(j,Lc-La)
      grad(1,2) = grad(1,2) + dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(1,2) = grad(1,2) + dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,2) = grad(1,2) + Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(1)-dLa(1))*(dLc(2)-dLa(2))

      grad(2,1) = grad(1,2)

      grad(2,2) = ddPhi(i,Lb-La)*(dLb(2)-dLa(2))**2*Phi(j,Lc-La)
      grad(2,2) = grad(2,2) + dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,2) = grad(2,2) + dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,2) = grad(2,2) + Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(2)-dLa(2))**2
    END FUNCTION SD_ddQuadBubblePBasis


!  --- end serendipity quad ---



    
!------------------------------------------------------------------------------
!>     Quadrilateral edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION QuadEdgePBasis(edge, i, u, v, invertEdge) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of quadrilaterals edge function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of quadrilaterals edge function i at point (u,v), i.e.
!       value = N_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN) :: edge, i
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: value, PhiPar, La, Lb, Na, Nb
      INTEGER :: nodes(2)
      LOGICAL :: invert 
      
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      value = 0
      ! Parameter validity check
      IF (edge < 1 .OR. edge > 4) THEN
         CALL Fatal('PElementBase::QuadEdgePBasis','Unknown edge for quad.')
      END IF

      ! Get nodes of edge
      nodes(1:2) = getQuadEdgeMap(edge)

      ! Bilinear nodal functions
      Na = QuadNodalPBasis(nodes(1),u,v)
      Nb = QuadNodalPBasis(nodes(2),u,v)

      ! Affine functions for edge direction
      La = QuadL(nodes(1),u,v)
      Lb = QuadL(nodes(2),u,v)

      ! For inverted edges swap direction
      phiPar = Lb-La
      IF (invert) phiPar = -PhiPar

      ! Get value of edge function
      value = Na*Nb*varPhi(i,phiPar)
    END FUNCTION QuadEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of quadrilateral edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION dQuadEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of quadrilaterals edge function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals edge function i at point (u,v), i.e.
!       grad = dN_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v
      LOGICAL, OPTIONAL :: invertEdge

      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp) :: Na,Nb,La,Lb, vPhi, PhiPar, dVPhi(2)
      REAL (KIND=dp), DIMENSION(2) :: dNa, dNb, dLa, dLb, grad,dPhiPar
      INTEGER :: nodes(2), swap

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 4) THEN
         CALL Fatal('PElementBase::dQuadPyraEdgePBasis','Unknown edge for quad.')
      END IF

      ! Get nodes of edge
      nodes(1:2) = getQuadEdgeMap(edge)

      ! Trilinear nodal functions and their derivatives
      Na  = QuadNodalPBasis(nodes(1),u,v)
      Nb  = QuadNodalPBasis(nodes(2),u,v)

      dNa = dQuadNodalPBasis(nodes(1),u,v)
      dNb = dQuadNodalPBasis(nodes(2),u,v)

      ! Affine functions and their derivatives for edge direction
      La  = QuadL(nodes(1),u,v)
      Lb  = QuadL(nodes(2),u,v)

      dLa = dQuadL(nodes(1),u,v)
      dLb = dQuadL(nodes(2),u,v)

      PhiPar = Lb-La
      dPhiPar = dLb-dLa

      IF(Invert) THEN
        PhiPar = -PhiPar
        dPhiPar = -dPhiPar
      END IF

      ! Get value of edge function
      vPhi  = varPhi(i,Phipar)
      dvPhi = dvarPhi(i,Phipar)*dPhiPar
      grad = dNa*Nb*vPhi + Na*dNb*vPhi + Na*Nb*dvPhi
    END FUNCTION dQuadEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of quadrilateral edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION ddQuadEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of quadrilaterals edge function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals edge function i at point (u,v), i.e.
!       grad = dN_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v
      LOGICAL, OPTIONAL :: invertEdge

      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp), DIMENSION(2) :: dNa, dNb, dLa, dLb, dPhiPar
      REAL (KIND=dp) :: Na, Nb, La, Lb, vPhi, PhiPar, dVPhi(2)
      REAL (KIND=dp) :: ddNa(2,2), ddNb(2,2), ddVPhi(2,2),grad(2,2)
      INTEGER :: nodes(2), swap,  p, q

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 4) THEN
         CALL Fatal('PElementBase::dQuadEdgePBasis','Unknown edge for quad.')
      END IF

      ! Get nodes of edge
      nodes(1:2) = getQuadEdgeMap(edge)

      ! Trilinear nodal functions and their derivatives
      Na  = QuadNodalPBasis(nodes(1),u,v)
      Nb  = QuadNodalPBasis(nodes(2),u,v)

      dNa = dQuadNodalPBasis(nodes(1),u,v)
      dNb = dQuadNodalPBasis(nodes(2),u,v)

      ddNa = ddQuadNodalPBasis(nodes(1),u,v)
      ddNb = ddQuadNodalPBasis(nodes(2),u,v)

      ! For inverted edges swap direction
      IF (invert) THEN
        swap=nodes(1); nodes(1)=nodes(2); nodes(2)=swap
      END IF

      ! Affine functions and their derivatives for edge direction
      La  = QuadL(nodes(1),u,v)
      Lb  = QuadL(nodes(2),u,v)

      dLa = dQuadL(nodes(1),u,v)
      dLb = dQuadL(nodes(2),u,v)

      PhiPar = Lb-La
      dPhiPar = dLb-dLa

      ! Get value of edge function
      vPhi  = VarPhi(i,Phipar)
      dvPhi = dVarPhi(i,Phipar)*dPhiPar
      DO p=1,2
        DO q=p,2
          ddVPhi(p,q) = ddVarPhi(i,PhiPar)*dPhiPar(p)*dPhiPar(q)
        END DO
      END DO

!     grad = dNa*Nb*vPhi + Na*dNb*vPhi + Na*Nb*dVPhi
      grad(1,1) = ddNa(1,1)*Nb*Vphi + dNa(1)*dNb(1)*vPhi + dNa(1)*Nb*dVPhi(1) + &
                  dNa(1)*dNb(1)*Vphi + Na*ddNb(1,1)*vPhi + Na*dNb(1)*dVPhi(1) + &
                  dNa(1)*Nb*dVphi(1) + Na*dNb(1)*dvPhi(1) + Na*Nb*ddVPhi(1,1)

      grad(1,2) = ddNa(1,2)*Nb*Vphi + dNa(2)*dNb(1)*vPhi + dNa(2)*Nb*dVPhi(1) + &
                  dNa(1)*dNb(2)*Vphi + Na*ddNb(1,2)*vPhi + Na*dNb(2)*dVPhi(1) + &
                  dNa(1)*Nb*dVphi(2) + Na*dNb(1)*dvPhi(2) + Na*Nb*ddVPhi(1,2)

      grad(2,2) = ddNa(2,2)*Nb*Vphi + dNa(2)*dNb(2)*vPhi + dNa(2)*Nb*dVPhi(2) + &
                  dNa(2)*dNb(2)*Vphi + Na*ddNb(2,2)*vPhi + Na*dNb(2)*dVPhi(2) + &
                  dNa(2)*Nb*dVphi(2) + Na*dNb(2)*dvPhi(2) + Na*Nb*ddVPhi(2,2)

      grad(2,1) = grad(1,2)
    END FUNCTION ddQuadEdgePBasis


!------------------------------------------------------------------------------
!>     Quadrilateral bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION QuadBubblePBasis(i,j,u,v,localNumbers) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of edge function, (i,j) = {(2,2),(2,3),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of quarilateral to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of quadrilateral faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of quadrilaterals bubble function (i,j) at point (u,v), 
!       i.e. value = N_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(4)
      REAL (KIND=dp) :: La, Lb, Lc, value, Pa, Pb

      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
        value  = Phi(i+2,u)*Phi(j+2,v)
        RETURN
      END IF
      
      ! Numbering present, so use it
      La = QuadL(localNumbers(1),u,v)
      Lb = QuadL(localNumbers(2),u,v)
      Lc = QuadL(localNumbers(4),u,v)

      Pa = QuadNodalPBasis(LocalNumbers(1),u,v)
      Pb = QuadNodalPBasis(LocalNumbers(3),u,v)

      ! Calculate value of function from general form
      value = Pa*Pb*LegendreP(i,Lb-La)*LegendreP(j,Lc-La)
    END FUNCTION QuadBubblePBasis


!------------------------------------------------------------------------------
!>     Gradient of quadrilateral bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION dQuadBubblePBasis(i,j,u,v, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of edge function, (i,j) = {(2,2),(2,3),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of quarilateral to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of quadrilateral faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals bubble function (i,j) at point (u,v), 
!       i.e. grad = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(4)
      REAL(Kind=dp) :: La, Lb, Lc, Pa, Pb, Legi, Legj
      REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, grad, dPa, dPb, dLegi, dLegj

      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
         grad(1) = dPhi(i+2,u)*Phi(j+2,v)
         grad(2) = Phi(i+2,u)*dPhi(j+2,v)
         RETURN
      END IF

      ! Numbering present, so use it
      La = QuadL(localNumbers(1),u,v)
      Lb = QuadL(localNumbers(2),u,v)
      Lc = QuadL(localNumbers(4),u,v)

      dLa = dQuadL(localNumbers(1),u,v)
      dLb = dQuadL(localNumbers(2),u,v)
      dLc = dQuadL(localNumbers(4),u,v)

      Pa = QuadNodalPBasis(localNumbers(1),u,v)
      Pb = QuadNodalPBasis(localNumbers(3),u,v)

      dPa = dQuadNodalPBasis(localNumbers(1),u,v)
      dPb = dQuadNodalPBasis(localNumbers(3),u,v)

      Legi = LegendreP(i,Lb-La)
      Legj = LegendreP(j,Lc-La)

      dLegi = dLegendreP(i,Lb-La)*(dLb-dLa)
      dLegj = dLegendreP(j,Lc-La)*(dLc-dLa)

      grad = dPa*Pb*Legi*Legj + Pa*dPb*Legi*Legj + &
             Pa*Pb*dLegi*Legj + Pa*Pb*Legi*dLegj
    END FUNCTION dQuadBubblePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of quadrilateral bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION ddQuadBubblePBasis(i,j,u,v, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of edge function, (i,j) = {(2,2),(2,3),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of quarilateral to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of quadrilateral faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of quadrilaterals bubble function (i,j) at point (u,v), 
!       i.e. grad = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(4)
      INTEGER :: p,q, local(4)
      REAL(Kind=dp) :: La, Lb, Lc, Legi, Legj, Pa, Pb
      REAL(Kind=dp), DIMENSION(2,2) :: ddLegi, ddLegj, grad, ddPa, ddPb
      REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, dLegi, dLegj, dPa, dPb
      
      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
         grad(1,1) = ddPhi(i+2,u)*Phi(j+2,v)
         grad(1,2) = dPhi(i+2,u)*dPhi(j+2,v)
         grad(2,1) = dPhi(i+2,u)*dPhi(j+2,v)
         grad(2,2) = Phi(i+2,u)*ddPhi(j+2,v)
         RETURN
      END IF

      ! Numbering present, so use it
      La = QuadL(localNumbers(1),u,v)
      Lb = QuadL(localNumbers(2),u,v)
      Lc = QuadL(localNumbers(4),u,v)

      dLa = dQuadL(localNumbers(1),u,v)
      dLb = dQuadL(localNumbers(2),u,v)
      dLc = dQuadL(localNumbers(4),u,v)

      Pa = QuadNodalPBasis(localNumbers(1),u,v)
      Pb = QuadNodalPBasis(localNumbers(3),u,v)

      dPa = dQuadNodalPBasis(localNumbers(1),u,v)
      dPb = dQuadNodalPBasis(localNumbers(3),u,v)

      ddPa = ddQuadNodalPBasis(localNumbers(1),u,v)
      ddPb = ddQuadNodalPBasis(localNumbers(3),u,v)

      Legi = LegendreP(i,Lb-La)
      Legj = LegendreP(j,Lc-La)

      dLegi = dLegendreP(i,Lb-La)*(dLb-dLa)
      dLegj = dLegendreP(j,Lc-La)*(dLc-dLa)

      DO p=1,2
        DO q=p,2
          ddLegi(p,q) = ddLegendreP(i,Lb-La)*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
          ddLegj(p,q) = ddLegendreP(j,Lc-La)*(dLc(p)-dLa(p))*(dLc(q)-dLa(q))
        END DO
      END DO

!     grad = dPa*Pb*Legi*Legj + Pa*dPb*Legi*Legj + &
!            Pa*Pb*dLegi*Legj + Pa*Pb*Legi*dLegj

      grad = 0
      grad(1,1) = grad(1,1) + ddPa(1,1)*Pb*Legi*Legj + dPa(1)*dPb(1)*Legi*Legj + &
                  dPa(1)*Pb*dLegi(1)*Legj + dPa(1)*Pb*Legi*dLegj(1)

      grad(1,1) = grad(1,1) + dPa(1)*dPb(1)*Legi*Legj + Pa*ddPb(1,1)*Legi*Legj + &
                  Pa*dPb(1)*dLegi(1)*Legj + Pa*dPb(1)*Legi*dLegj(1)

      grad(1,1) = grad(1,1) + dPa(1)*Pb*dLegi(1)*Legj + Pa*dPb(1)*dLegi(1)*Legj + &
                  Pa*Pb*ddLegi(1,1)*Legj + Pa*Pb*dLegi(1)*dLegj(1)

      grad(1,1) = grad(1,1) + dPa(1)*Pb*Legi*dLegj(1) + Pa*dPb(1)*Legi*dLegj(1) + &
                  Pa*Pb*dLegi(1)*dLegj(1) + Pa*Pb*Legi*ddLegj(1,1)

      grad(1,2) = grad(1,2) + ddPa(1,2)*Pb*Legi*Legj + dPa(2)*dPb(1)*Legi*Legj + &
                  dPa(2)*Pb*dLegi(1)*Legj + dPa(2)*Pb*Legi*dLegj(1)

      grad(1,2) = grad(1,2) + dPa(1)*dPb(2)*Legi*Legj + Pa*ddPb(1,2)*Legi*Legj + &
                  Pa*dPb(2)*dLegi(1)*Legj + Pa*dPb(2)*Legi*dLegj(1)

      grad(1,2) = grad(1,2) + dPa(1)*Pb*dLegi(2)*Legj + Pa*dPb(1)*dLegi(2)*Legj + &
                  Pa*Pb*ddLegi(1,2)*Legj + Pa*Pb*dLegi(2)*dLegj(1)

      grad(1,2) = grad(1,2) + dPa(1)*Pb*Legi*dLegj(2) + Pa*dPb(1)*Legi*dLegj(2) + &
                  Pa*Pb*dLegi(1)*dLegj(2) + Pa*Pb*Legi*ddLegj(1,2)

      grad(2,2) = grad(2,2) + ddPa(2,2)*Pb*Legi*Legj + dPa(2)*dPb(2)*Legi*Legj + &
                  dPa(2)*Pb*dLegi(2)*Legj + dPa(2)*Pb*Legi*dLegj(2)

      grad(2,2) = grad(2,2) + dPa(2)*dPb(2)*Legi*Legj + Pa*ddPb(2,2)*Legi*Legj + &
                  Pa*dPb(2)*dLegi(2)*Legj + Pa*dPb(2)*Legi*dLegj(2)

      grad(2,2) = grad(2,2) + dPa(2)*Pb*dLegi(2)*Legj + Pa*dPb(2)*dLegi(2)*Legj + &
                  Pa*Pb*ddLegi(2,2)*Legj + Pa*Pb*dLegi(2)*dLegj(2)

      grad(2,2) = grad(2,2) + dPa(2)*Pb*Legi*dLegj(2) + Pa*dPb(2)*Legi*dLegj(2) + &
                  Pa*Pb*dLegi(2)*dLegj(2) + Pa*Pb*Legi*ddLegj(2,2)

      grad(2,1) = grad(1,2)

    END FUNCTION ddQuadBubblePBasis


!------------------------------------------------------------------------------
!>     Defines linear functions for quadrilateral nodes. These are used in 
!>     calculation of changing parameters for bubbles of quadrilateral if 
!>     directional function values are requested. 
!------------------------------------------------------------------------------
    PURE FUNCTION QuadL(which, u, v) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: which
!      INPUT: Node in which calculate function value, which = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of quadrilaterals linear nodal function at point (u,v), 
!       i.e. value = N_i^l(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: which
      REAL(Kind=dp), INTENT(IN) :: u, v
      REAL(Kind=dp) :: value

      value = 0.0_dp
      SELECT CASE (which)
      CASE (1)
         value = (2-u-v)/2
      CASE (2)
         value = (2+u-v)/2
      CASE (3)
         value = (2+u+v)/2
      CASE (4)
         value = (2-u+v)/2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::QuadL', 'Unknown helper function L for quad')
#endif
      END SELECT
    END FUNCTION QuadL


!------------------------------------------------------------------------------
!>     Defines gradients of linear functions for quadrilateral nodes. For use 
!>     see QuadL. 
!------------------------------------------------------------------------------
    PURE FUNCTION dQuadL(which, u, v) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: which
!      INPUT: Node in which calculate function value, which = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of quadrilaterals linear nodal function at point (u,v), 
!       i.e. value = dN_i^l(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: which
      REAL(Kind=dp), INTENT(IN) :: u, v
      REAL(Kind=dp) :: grad(2)
      
      SELECT CASE (which)
      CASE (1)
         grad(1:2) = [ -1d0/2, -1d0/2 ]
      CASE (2)
         grad(1:2) = [  1d0/2, -1d0/2 ]
      CASE (3)
         grad(1:2) = [  1d0/2,  1d0/2 ]
      CASE (4)
         grad(1:2) = [ -1d0/2,  1d0/2 ]
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::dQuadL', 'Unknown helper function dL for quad')
#endif
      END SELECT
    END FUNCTION dQuadL


!------------------------------------------------------------------------------
!>     Triangle nodal basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION TriangleNodalPBasis(node, u, v) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of triangles nodal function to calculate
!        node = {1,2,3}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of triangles nodal function at point (u,v), i.e.
!       value = N_i(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: value

      value = 0
      SELECT CASE(node)
      CASE (1)
         value = (1-u-v/SQRT(3d0))/2
      CASE (2)
         value = (1+u-v/SQRT(3d0))/2
      CASE (3)
         value = v/SQRT(3d0)
      CASE DEFAULT
         CALL Fatal('PElementBase::TriangleNodalPBasis', 'Unknown node for triangle')
      END SELECT 
    END FUNCTION TriangleNodalPBasis



    SUBROUTINE TriangleNodalPBasisAll(u, v, phi) 
      IMPLICIT NONE
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: phi(:)
      REAL(KIND=dp), PARAMETER :: half=1.0_dp/2, c3=1.0_dp/SQRT(3.0_dp)
      
      phi(1) = half*(1-u-c3*v)
      phi(2) = half*(1+u-c3*v)
      phi(3) = c3*v
    END SUBROUTINE TriangleNodalPBasisAll


    SUBROUTINE TriangleNodalLBasisAll(u, v, phi) 
      IMPLICIT NONE
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: phi(:)
      
      phi(1) = 1.0_dp-u-v
      phi(2) = u
      phi(3) = v
    END SUBROUTINE TriangleNodalLBasisAll
    
!------------------------------------------------------------------------------
!>     Gradient of triangle nodal basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION dTriangleNodalPBasis(node, u, v) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of triangles nodal function to calculate
!        node = {1,2,3}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of triangles nodal function at point (u,v), i.e.
!       grad = dN_i(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v
      ! Return value
      REAL (KIND=dp), DIMENSION(2) :: grad

      grad = 0
      SELECT CASE(node)
      CASE (1)
         grad(1) = -1d0/2
         grad(2) = -SQRT(3d0)/6
      CASE (2)
         grad(1) = 1d0/2
         grad(2) = -SQRT(3d0)/6
      CASE (3)
         grad(1) = 0
         grad(2) = SQRT(3d0)/3
      CASE DEFAULT
         CALL Fatal('PElementBase::dTriangleNodalPBasis', 'Unknown node for triangle')
      END SELECT
    END FUNCTION dTriangleNodalPBasis
!------------------------------------------------------------------------------

!>     2nd derivatives of triangle nodal basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION ddTriangleNodalPBasis(node, u, v) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of triangles nodal function to calculate
!        node = {1,2,3}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of triangles nodal function at point (u,v), i.e.
!       grad = dN_i(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v
      ! Return value
      REAL (KIND=dp), DIMENSION(2,2) :: grad

      grad = 0
    END FUNCTION ddTriangleNodalPBasis


    SUBROUTINE dTriangleNodalPBasisAll(u, v, gradphi)
      IMPLICIT NONE
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: gradphi(:,:)
      REAL(KIND=dp), PARAMETER :: half=1.0_dp/2, c6=SQRT(3.0_dp)/6.0_dp

      gradphi(1,1) = -half
      gradphi(1,2) = -c6
      gradphi(2,1) = half
      gradphi(2,2) = -c6
      gradphi(3,1) = 0
      gradphi(3,2) = 2*c6
    END SUBROUTINE dTriangleNodalPBasisAll

    SUBROUTINE dTriangleNodalLBasisAll(u, v, gradphi)
      IMPLICIT NONE
      REAL (KIND=dp), INTENT(IN) :: u,v
      REAL (KIND=dp) :: gradphi(:,:)

      gradphi(1,1) = -1.0_dp
      gradphi(1,2) = -1.0_dp
      gradphi(2,1) = 1.0_dp
      gradphi(2,2) = 0.0_dp
      gradphi(3,1) = 0.0_dp
      gradphi(3,2) = 1.0_dp
    END SUBROUTINE dTriangleNodalLBasisAll

    
!------------------------------------------------------------------------------
!>     Triangle edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION TriangleEdgePBasis(edge, i, u, v, invertEdge) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of triangles edge function to calculate
!        edge = {1,2,3}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of triangles edge function i at point (u,v), i.e.
!       value = N_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      LOGICAL, INTENT(IN), OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v

      REAL (KIND=dp) :: L1, L2, L3, value
      LOGICAL :: invert

      ! Check if edge needs to be inverted. The default is not inverted edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      value = 0
      SELECT CASE(edge)
      CASE (1)
         L1 = TriangleNodalPBasis(1,u,v)
         L2 = TriangleNodalPBasis(2,u,v)

         ! Invert edge for parity if needed
         IF (.NOT. invert) THEN
            value = L1*L2*varPhi(i,L2-L1)
         ELSE
            value = L1*L2*varPhi(i,L1-L2)
         END IF
      CASE (2)
         L2 = TriangleNodalPBasis(2,u,v)
         L3 = TriangleNodalPBasis(3,u,v)

         ! Invert edge for parity if needed
         IF (.NOT. invert) THEN
            value = L2*L3*varPhi(i,L3-L2)
         ELSE
            value = L2*L3*varPhi(i,L2-L3)
         END IF
      CASE (3)
         L1 = TriangleNodalPBasis(1,u,v)
         L3 = TriangleNodalPBasis(3,u,v)

         ! Invert edge for parity if needed
         IF (.NOT. invert) THEN
            value = L1*L3*varPhi(i,L1-L3)
         ELSE
            value = L1*L3*varPhi(i,L3-L1)
         END IF
      CASE DEFAULT
         CALL Fatal('PElementBase::TriangleEdgePBasis', 'Unknown edge for triangle')
      END SELECT 
    END FUNCTION TriangleEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of triangle edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION dTriangleEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of triangles edge function to calculate
!        edge = {1,2,3}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of triangles edge function i at point (u,v), i.e.
!       grad = dN_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u,v
      LOGICAL, OPTIONAL :: invertEdge
      ! Return value
      REAL (KIND=dp), DIMENSION(2) :: grad
      ! Variables
      REAL (KIND=dp) :: L1, L2, L3, L3_L2, L1_L3, vPhi
      LOGICAL :: invert

      
      ! Check if edge needs to be inverted. The default is not inverted edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      grad = 0
      SELECT CASE(edge)
      CASE (1)
         L1 = TriangleNodalPBasis(1,u,v)
         L2 = TriangleNodalPBasis(2,u,v)
         
         ! Invert edge 
         IF (.NOT. invert) THEN
            vPhi = varPhi(i,u)
            grad(1) = -1d0/2*L2*vPhi+1d0/2*L1*vPhi+L1*L2*dVarPhi(i,u) 
            grad(2) = -SQRT(3d0)/6*L2*vPhi-SQRT(3d0)/6*L1*vPhi
         ELSE
            vPhi = varPhi(i,-u)
            grad(1) = -1d0/2*L2*vPhi+1d0/2*L1*vPhi-L1*L2*dVarPhi(i,-u) 
            grad(2) = -SQRT(3d0)/6*L2*vPhi-SQRT(3d0)/6*L1*vPhi   
         END IF
      CASE (2)
         L2 = TriangleNodalPBasis(2,u,v)
         L3 = TriangleNodalPBasis(3,u,v)

         ! Invert edge
         IF (.NOT. invert) THEN 
            L3_L2 = L3-L2
            vPhi = varPhi(i,L3_L2)
            grad(1) = 1d0/2*L3*vPhi-1d0/2*L2*L3*dVarPhi(i,L3_L2)
            grad(2) = -SQRT(3d0)/6*L3*vPhi+SQRT(3d0)/3*L2*vPhi+SQRT(3d0)/2*L2*L3*dVarPhi(i,L3_L2)
         ELSE
            L3_L2 = L2-L3
            vPhi = varPhi(i,L3_L2)
            grad(1) = 1d0/2*L3*vPhi+1d0/2*L2*L3*dVarPhi(i,L3_L2)
            grad(2) = -SQRT(3d0)/6*L3*vPhi+SQRT(3d0)/3*L2*vPhi-SQRT(3d0)/2*L2*L3*dVarPhi(i,L3_L2)
         END IF
      CASE (3)
         L1 = TriangleNodalPBasis(1,u,v)
         L3 = TriangleNodalPBasis(3,u,v)

         ! Invert edge
         IF (.NOT. invert) THEN
            L1_L3 = L1-L3
            vPhi = varPhi(i,L1_L3)
            grad(1) = -1d0/2*L3*vPhi-1d0/2*L1*L3*dVarPhi(i,L1_L3)
            grad(2) = -SQRT(3d0)/6*L3*vPhi+SQRT(3d0)/3*L1*vPhi-SQRT(3d0)/2*L1*L3*dVarPhi(i,L1_L3)
         ELSE
            L1_L3 = L3-L1
            vPhi = varPhi(i,L1_L3)
            grad(1) = -1d0/2*L3*vPhi+1d0/2*L1*L3*dVarPhi(i,L1_L3)
            grad(2) = -SQRT(3d0)/6*L3*vPhi+SQRT(3d0)/3*L1*vPhi+SQRT(3d0)/2*L1*L3*dVarPhi(i,L1_L3)
         END IF
      CASE DEFAULT
         CALL Fatal('PElementBase::dTriangleEdgePBasis', 'Unknown edge for triangle')
      END SELECT 
    END FUNCTION dTriangleEdgePBasis
      

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>     2nd derivatives of triangle edge basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION ddTriangleEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of triangles edge function to calculate
!        edge = {1,2,3}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of triangles edge function i at point (u,v), i.e.
!       grad = dN_i^{edge}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u,v
      LOGICAL, OPTIONAL :: invertEdge
      ! Return value
      REAL (KIND=dp), DIMENSION(2,2) :: grad
      ! Variables
      REAL (KIND=dp) :: LA,LB,dLAu,dLAv,dLBu,dLBv,swap,dL1u,dL2u,dL3u,dL1v,dL2v,dL3v,s
      REAL (KIND=dp) :: dd,dv,vp,varArg
      LOGICAL :: invert
      REAL(KIND=dp), PARAMETER :: half=1.0_dp/2, c3=1.0_dp/SQRT(3.0_dp)

      
      ! Check if edge needs to be inverted. The default is not inverted edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      grad = 0

      dL1u = -half 
      dL2u =  half 
      dL3u =  0

      dL1v = -half*c3
      dL2v = -half*c3
      dL3v =  c3

      SELECT CASE(edge)
      CASE (1)
         dLAu = dL1u
         dLBu = dL2u
         dLAv = dL1v
         dLBv = dL2v
         LA = TriangleNodalPBasis(1,u,v)
         LB = TriangleNodalPBasis(2,u,v)
      CASE(2)
         dLAu = dL2u
         dLBu = dL3u
         dLAv = dL2v
         dLBv = dL3v
         LA = TriangleNodalPBasis(2,u,v)
         LB = TriangleNodalPBasis(3,u,v)
      CASE(3)
         dLAu = dL3u
         dLBu = dL1u
         dLAv = dL3v
         dLBv = dL1v
         LA = TriangleNodalPBasis(3,u,v)
         LB = TriangleNodalPBasis(1,u,v)
      CASE DEFAULT
         CALL Fatal('PElementBase::dTriangleEdgePBasis', 'Unknown edge for triangle')
      END SELECT 

      s = 1; varArg = LB-LA
      IF(Invert) THEN
        s = -1; varArg = -varArg
      END IF

      dd = ddVarPhi(i,varArg)
      vp = VarPhi(i,varArg)
      dv = s*dVarPhi(i,varArg)

      grad(1,1) = LA*LB*(dLBu-dLAu)**2*dd + &
        2*(dLAu*LB+LA*dLBu)*(dLBu-dLAu)*dv + (dLAu*dLBu+dLAu*dLBu)*vp

      grad(1,2) = LA*LB*(dLBu-dLAu)*(dLBv-dLAv)*dd + &
                  (dLAv*LB+LA*dLBv)*(dLBu-dLAu)*dv + &
                  (dLAu*LB+LA*dLBu)*(dLBv-dLAv)*dv + &
                     (dLAu*dLBv+dLAv*dLBu)*vp

      grad(2,2) = LA*LB*(dLBv-dLAv)**2*dd + &
        2*(dLAv*LB+LA*dLBv)*(dLBv-dLAv)*dv + (dLAv*dLBv+dLAv*dLBv)*vp

      grad(2,1) = grad(1,2)
    END FUNCTION ddTriangleEdgePBasis
      

!------------------------------------------------------------------------------
!>     Triangle bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION TriangleBubblePBasis(j,n,u,v,localNumbers) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of bubble function, (i,j) = {(0,0),(0,1),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of triangle to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of triangular faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of triangles bubble function (i,j) at point (u,v), 
!       i.e. value = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: j, n
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(3)
      REAL (KIND=dp) :: La, Lb, Lc, value
      INTEGER :: local(3)
      
      ! If local numbering present, use it
      IF (PRESENT(localNumbers)) THEN
         local(1:3) = localNumbers(1:3)
      ELSE
      ! Local numbering not present. Use default numbering
         local(1:3) = [ 1,2,3 ]
      END IF

      La = TriangleNodalPBasis(local(1),u,v)
      Lb = TriangleNodalPBasis(local(2),u,v)
      Lc = TriangleNodalPBasis(local(3),u,v)
 
      value = La*Lb*Lc*((Lb-La)**j)*((2*Lc-1)**n)
    END FUNCTION TriangleBubblePBasis


    FUNCTION TriangleEBubblePBasis(i,j,u,v,localNumbers) RESULT(value)
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: i,j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(3)
      REAL (KIND=dp) :: La, Lb, Lc, value
      INTEGER :: local(3)
      
      ! If local numbering present, use it
      IF (PRESENT(localNumbers)) THEN
         local(1:3) = localNumbers(1:3)
      ELSE
      ! Local numbering not present. Use default numbering
         local(1:3) = [ 1,2,3 ]
      END IF

      La = TriangleNodalPBasis(local(1),u,v)
      Lb = TriangleNodalPBasis(local(2),u,v)
      Lc = TriangleNodalPBasis(local(3),u,v)
 
      value = La*Lb*Lc*LegendreP(i,Lb-La)*LegendreP(j,2*Lc-1)
    END FUNCTION TriangleEBubblePBasis


!------------------------------------------------------------------------------
!>     Gradient of triangles bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION dTriangleBubblePBasis(j,n,u,v,localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of bubble function, (i,j) = {(0,0),(0,1),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of triangle to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of triangular faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2)
!       gradient of triangles bubble function (i,j) at point (u,v), 
!       i.e. grad = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: j, n
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(3)
      ! Variables
      REAL (KIND=dp) :: La, Lb, Lc, Lc_1, Lb_Laj, Lc_1n
      REAL (KIND=dp), DIMENSION(2) :: dLa, dLb, dLc, grad
      INTEGER :: local(3)

      ! If local numbering present, use it
      IF (PRESENT(localNumbers)) THEN
         local(1:3) = localNumbers(1:3)
      ELSE
         ! Local numbering not present. Use default numbering
         local(1:3) = [ 1,2,3 ]
      END IF

      La = TriangleNodalPBasis(local(1),u,v)
      Lb = TriangleNodalPBasis(local(2),u,v)
      Lc = TriangleNodalPBasis(local(3),u,v)
      dLa = dTriangleNodalPBasis(local(1),u,v)
      dLb = dTriangleNodalPBasis(local(2),u,v)
      dLc = dTriangleNodalPBasis(local(3),u,v)

      Lb_Laj = toExp(Lb-La,j)
      Lc_1n = toExp(2*Lc-1,n)

      ! Calculate value of function from general form
      grad = dLa*Lb*Lc*Lb_Laj*Lc_1n + La*dLb*Lc*Lb_Laj*Lc_1n + &
           La*Lb*dLc*Lb_Laj*Lc_1n + La*Lb*Lc*j*toExp(Lb-La,j-1)*(dLb-dLa)*Lc_1n + &
           La*Lb*Lc*Lb_Laj*n*toExp(2*Lc-1,n-1)*(2d0*dLc)
    END FUNCTION dTriangleBubblePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of triangles bubble basis at point (u,v).
!------------------------------------------------------------------------------
    FUNCTION ddTriangleBubblePBasis(j,n,u,v,localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j
!      INPUT: index of bubble function, (i,j) = {(0,0),(0,1),...}
!
!    REAL(KIND=dp) :: u,v
!      INPUT: point at which to evaluate function
! 
!    INTEGER, OPTIONAL :: localNumbers(4) 
!      INPUT: local numbering of triangle to define direction of bubble
!        function. Used with 3d element boundary integrals to give correct 
!        directions to bubble functions of triangular faces. 
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(2,2)
!       gradient of triangles bubble function (i,j) at point (u,v), 
!       i.e. grad = dN_{m(i,j)}^{0}(u,v)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=dp) :: grad(2,2)

      ! Parameters
      INTEGER, INTENT(IN) :: j, n
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(3)
      ! Variables
      REAL (KIND=dp), DIMENSION(2) :: dLa, dLb, dLc, dLb_Laj, dLc_1n
      REAL (KIND=dp) :: La, Lb, Lc, Lc_1, Lb_Laj, Lc_1n, ddLb_Laj, ddLc_1n
      INTEGER :: local(3)

      ! If local numbering present, use it
      IF (PRESENT(localNumbers)) THEN
         local(1:3) = localNumbers(1:3)
      ELSE
         ! Local numbering not present. Use default numbering
         local(1:3) = [ 1,2,3 ]
      END IF

      La = TriangleNodalPBasis(local(1),u,v)
      Lb = TriangleNodalPBasis(local(2),u,v)
      Lc = TriangleNodalPBasis(local(3),u,v)

      dLa = dTriangleNodalPBasis(local(1),u,v)
      dLb = dTriangleNodalPBasis(local(2),u,v)
      dLc = dTriangleNodalPBasis(local(3),u,v)

      Lb_Laj = toExp(Lb-La,j)
      Lc_1n = toExp(2*Lc-1,n)

      dLb_Laj = j*toExp(Lb-La,j-1) * (dLb-dLa)
      dLc_1n  = n*toExp(2*Lc-1,n-1) * 2*dLc

      ddLb_Laj = j*(j-1)*toExp(Lb-La,j-2)
      ddLc_1n = n*(n-1)*toExp(2*Lc-1,n-2)

      ! Calculate value of function from general form
!     grad = dLa*Lb*Lc*Lb_Laj*Lc_1n + La*dLb*Lc*Lb_Laj*Lc_1n + &
!          La*Lb*dLc*Lb_Laj*Lc_1n + La*Lb*Lc*j*toExp(Lb-La,j-1)*(dLb-dLa)*Lc_1n + &
!          La*Lb*Lc*Lb_Laj*n*toExp(2*Lc-1,n-1)*(2d0*dLc)
!     value = La*Lb*Lc*((Lb-La)**j)*((2*Lc-1)**n)

      grad(1,1) = &
         dLa(1)*(dLb(1)*(Lc*Lb_Laj*Lc_1n)+dLc(1)*(Lb*Lb_Laj*Lc_1n)+dLb_Laj(1)*(Lb*Lc*Lc_1n)+dLc_1n(1)*(Lb*Lc*Lb_Laj)) + &
         dLb(1)*(dLa(1)*(Lc*Lb_Laj*Lc_1n)+dLc(1)*(La*Lb_Laj*Lc_1n)+dLb_Laj(1)*(La*Lc*Lc_1n)+dLc_1n(1)*(La*Lc*Lb_Laj)) + &
         dLc(1)*(dLa(1)*(Lb*Lb_Laj*Lc_1n)+dLb(1)*(La*Lb_Laj*Lc_1n)+dLb_Laj(1)*(La*Lb*Lc_1n)+dLc_1n(1)*(La*Lb*Lb_Laj)) + &
         dLb_Laj(1)*(dLa(1)*(Lb*Lc*Lc_1n)+dLb(1)*(La*Lc*Lc_1n)+dLc(1)*(La*Lb*Lc_1n)+dLc_1n(1)*(La*Lb*Lc)) + &
         dLc_1n(1)*(dLa(1)*(Lb*Lc*Lb_Laj)+dLb(1)*(La*Lc*Lb_Laj)+dLc(1)*(La*Lb*Lb_Laj)+dLb_Laj(1)*(La*Lb*Lc)) + &
         ddLb_Laj*(La*Lb*Lc*Lc_1n)*(dLb(1)-dLa(1))**2 + ddLc_1n*(La*Lb*Lc*Lb_Laj)*4*dLc(1)**2

      grad(1,2) = &
        dLa(1)*(dLb(2)*(Lc*Lb_Laj*Lc_1n)+dLc(2)*(Lb*Lb_Laj*Lc_1n)+dLb_Laj(2)*(Lb*Lc*Lc_1n)+dLc_1n(2)*(Lb*Lc*Lb_Laj)) + &
        dLb(1)*(dLa(2)*(Lc*Lb_Laj*Lc_1n)+dLc(2)*(La*Lb_Laj*Lc_1n)+dLb_Laj(2)*(La*Lc*Lc_1n)+dLc_1n(2)*(La*Lc*Lb_Laj)) + &
        dLc(1)*(dLa(2)*(Lb*Lb_Laj*Lc_1n)+dLb(2)*(La*Lb_Laj*Lc_1n)+dLb_Laj(2)*(La*Lb*Lc_1n)+dLc_1n(2)*(La*Lb*Lb_Laj)) + &
        dLb_Laj(1)*(dLa(2)*(Lb*Lc*Lc_1n)+dLb(2)*(La*Lc*Lc_1n)+dLc(2)*(La*Lb*Lc_1n)+dLc_1n(2)*(La*Lb*Lc)) + &
        dLc_1n(1)*(dLa(2)*(Lb*Lc*Lb_Laj)+dLb(2)*(La*Lc*Lb_Laj)+dLc(2)*(La*Lb*Lb_Laj)+dLb_Laj(2)*(La*Lb*Lc)) + &
        ddLb_Laj*(La*Lb*Lc*Lc_1n)*(dLb(1)-dLa(1))*(dLb(2)-dLa(2)) + ddLc_1n*(La*Lb*Lc*Lb_Laj)*4*dLc(1)*dLc(2)

      grad(2,2) = &
        dLa(2)*(dLb(2)*(Lc*Lb_Laj*Lc_1n)+dLc(2)*(Lb*Lb_Laj*Lc_1n)+dLb_Laj(2)*(Lb*Lc*Lc_1n)+dLc_1n(2)*(Lb*Lc*Lb_Laj)) + &
        dLb(2)*(dLa(2)*(Lc*Lb_Laj*Lc_1n)+dLc(2)*(La*Lb_Laj*Lc_1n)+dLb_Laj(2)*(La*Lc*Lc_1n)+dLc_1n(2)*(La*Lc*Lb_Laj)) + &
        dLc(2)*(dLa(2)*(Lb*Lb_Laj*Lc_1n)+dLb(2)*(La*Lb_Laj*Lc_1n)+dLb_Laj(2)*(La*Lb*Lc_1n)+dLc_1n(2)*(La*Lb*Lb_Laj)) + &
        dLb_Laj(2)*(dLa(2)*(Lb*Lc*Lc_1n)+dLb(2)*(La*Lc*Lc_1n)+dLc(2)*(La*Lb*Lc_1n)+dLc_1n(2)*(La*Lb*Lc)) + &
        dLc_1n(2)*(dLa(2)*(Lb*Lc*Lb_Laj)+dLb(2)*(La*Lc*Lb_Laj)+dLc(2)*(La*Lb*Lb_Laj)+dLb_Laj(2)*(La*Lb*Lc)) + &
        ddLb_Laj*(La*Lb*Lc*Lc_1n)*(dLb(2)-dLa(2))**2 + ddLc_1n*(La*Lb*Lc*Lb_Laj)*4*dLc(2)**2

      grad(2,1) = grad(1,2)
    END FUNCTION ddTriangleBubblePBasis


    FUNCTION dTriangleEBubblePBasis(i,j,u,v,localNumbers) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i, j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(3)
      ! Variables
      REAL (KIND=dp) :: La, Lb, Lc, Lc_1, Legi, Legj
      REAL (KIND=dp), DIMENSION(2) :: dLa, dLb, dLc, grad
      INTEGER :: local(3)

      ! If local numbering present, use it
      IF (PRESENT(localNumbers)) THEN
         local(1:3) = localNumbers(1:3)
      ELSE
         ! Local numbering not present. Use default numbering
         local(1:3) = [ 1,2,3 ]
      END IF

      La = TriangleNodalPBasis(local(1),u,v)
      Lb = TriangleNodalPBasis(local(2),u,v)
      Lc = TriangleNodalPBasis(local(3),u,v)
      dLa = dTriangleNodalPBasis(local(1),u,v)
      dLb = dTriangleNodalPBasis(local(2),u,v)
      dLc = dTriangleNodalPBasis(local(3),u,v)

      Legi = LegendreP(i,Lb-La)
      Legj = LegendreP(j,2*Lc-1)

      ! Calculate value of function from general form
      grad = dLa*Lb*Lc*Legi*Legj + La*dLb*Lc*Legi*Legj + &
           La*Lb*dLc*Legi*Legj + La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb-dLa)*Legj + &
           La*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*(2*dLc)
    END FUNCTION dTriangleEBubblePBasis


    FUNCTION ddTriangleEBubblePBasis(i,j,u,v,localNumbers) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i, j
      REAL (KIND=dp), INTENT(IN) :: u,v
      INTEGER, OPTIONAL :: localNumbers(3)
      ! Variables
      REAL (KIND=dp), DIMENSION(2) :: dLa, dLb, dLc
      REAL (KIND=dp) :: La, Lb, Lc, Lc_1, Legi, Legj, grad(2,2)
      INTEGER :: local(3)

      REAL(KIND=dp) :: r,s,t
      REAL(KIND=dp) :: Lab,dLab(2),ddLab(2,2)
      REAL(KIND=dp) :: Labc,dLabc(2),ddLabc(2,2)
      REAL(KIND=dp) :: G1,dG1(2), ddG1(2,2)
      REAL(KIND=dp) :: G2,dG2(2), ddG2(2,2)
      REAL(KIND=dp) :: G12,dG12(2), ddG12(2,2)
      REAL(KIND=dp) :: LabcG12, dLabcG12(2), ddLabcG12(2,2)
      INTEGER :: p,q

      ! If local numbering present, use it
      IF (PRESENT(localNumbers)) THEN
         local(1:3) = localNumbers(1:3)
      ELSE
         ! Local numbering not present. Use default numbering
         local(1:3) = [ 1,2,3 ]
      END IF

      La = TriangleNodalPBasis(local(1),u,v)
      Lb = TriangleNodalPBasis(local(2),u,v)
      Lc = TriangleNodalPBasis(local(3),u,v)

      dLa = dTriangleNodalPBasis(local(1),u,v)
      dLb = dTriangleNodalPBasis(local(2),u,v)
      dLc = dTriangleNodalPBasis(local(3),u,v)

      Lab = La*Lb
      dLab = dLa*Lb + La*dLb
      DO p=1,2
        DO q=p,2
          ddLab(p,q) = dLa(p)*dLb(q) + dLa(q)*dLb(p)
        END DO        
      END DO        
 
      Labc = Lab*Lc
      dLabc = dLab*Lc + Lab*dLc
      DO p=1,2
        DO q=p,2
          ddLabc(p,q) = ddLab(p,q)*Lc + dLab(p)*dLc(q) + dLab(q)*dLc(p)
        END DO        
      END DO        

      G1 = LegendreP(i,Lb-La)
      G2 = LegendreP(j,2*Lc-1)

      dG1 = dLegendreP(i,Lb-La)*(dLb-dLa)
      dG2 = dLegendreP(j,2*Lc-1)*2*dLc

      r = ddLegendreP(i,Lb-La)
      s = ddLegendreP(j,2*Lc-1)
      DO p=1,2
        DO q=p,2
          ddG1(p,q) = r * (dLb(p)-dLa(p))*(dLb(q)-dLa(q))
          ddG2(p,q) = s * 4*dLc(p)*dLc(q)
        END DO
      END DO

      G12 = G1*G2
      dG12 = dG1*G2 + G1*dG2
      DO p=1,2
        DO q=p,2
          ddG12(p,q) = ddG1(p,q)*G2 + dG1(p)*dG2(q) + dG1(q)*dG2(p) + G1*ddG2(p,q)
        END DO
      END DO

      ! dd(La*Lb*Lc*G1*G2)
      DO p=1,2
        DO q=p,2
          grad(p,q)=ddLabc(p,q)*G12+dLabc(p)*dG12(q)+dLabc(q)*dG12(p)+Labc*ddG12(p,q)
        END DO
      END DO

      grad(2,1)=grad(1,2)
    END FUNCTION ddTriangleEBubblePBasis


    ! 3D ELEMENTS

    FUNCTION BrickNodalPBasis(node, u, v, w) RESULT(value)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(Kind=dp), INTENT(IN) :: u,v,w

      ! Variables
      REAL(Kind=dp) :: value

      SELECT CASE (node)
      CASE (1)
         value = (1-u)*(1-v)*(1-w)
      CASE (2)
         value = (1+u)*(1-v)*(1-w)
      CASE (3)
         value = (1+u)*(1+v)*(1-w)
      CASE (4)
         value = (1-u)*(1+v)*(1-w)
      CASE (5)
         value = (1-u)*(1-v)*(1+w)
      CASE (6)
         value = (1+u)*(1-v)*(1+w)
      CASE (7)
         value = (1+u)*(1+v)*(1+w)
      CASE (8)
         value = (1-u)*(1+v)*(1+w)
      CASE DEFAULT
         CALL Fatal('PElementBase::BrickNodalPBasis','Unknown node for brick')
      END SELECT
      value = value/8
    END FUNCTION BrickNodalPBasis

    ! As previous except obtain all nodal lvalues at once. 
    SUBROUTINE BrickNodalPBasisAll(u, v, w, phi) 
      IMPLICIT NONE
      
      ! Parameters
      REAL(Kind=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: phi(:)
      
      REAL(Kind=dp), PARAMETER :: c = 1.0_dp/8.0_dp      
      INTEGER, PARAMETER :: usgn(8) = [-1,1,1,-1,-1,1,1,-1]
      INTEGER, PARAMETER :: vsgn(8) = [-1,-1,1,1,-1,-1,1,1]
      INTEGER, PARAMETER :: wsgn(8) = [-1,-1,-1,-1,1,1,1,1]
            
      phi(1:8) = c*(1+usgn*u)*(1+vsgn*v)*(1+wsgn*w)
      
    END SUBROUTINE BrickNodalPBasisAll

    
    FUNCTION dBrickNodalPBasis(node, u, v, w) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(Kind=dp), INTENT(IN) :: u,v,w

      ! Variables
      REAL(Kind=dp) :: grad(3)

      SELECT CASE (node)
      CASE (1)
         grad(1) = -(1-v)*(1-w) 
         grad(2) = -(1-u)*(1-w)
         grad(3) = -(1-u)*(1-v)
      CASE (2)
         grad(1) =  (1-v)*(1-w) 
         grad(2) = -(1+u)*(1-w)
         grad(3) = -(1+u)*(1-v)
      CASE (3)
         grad(1) =  (1+v)*(1-w) 
         grad(2) =  (1+u)*(1-w)
         grad(3) = -(1+u)*(1+v)
      CASE (4)
         grad(1) = -(1+v)*(1-w) 
         grad(2) =  (1-u)*(1-w)
         grad(3) = -(1-u)*(1+v)
      CASE (5)
         grad(1) = -(1-v)*(1+w) 
         grad(2) = -(1-u)*(1+w)
         grad(3) =  (1-u)*(1-v)
      CASE (6)
         grad(1) =  (1-v)*(1+w) 
         grad(2) = -(1+u)*(1+w)
         grad(3) =  (1+u)*(1-v)
      CASE (7)
         grad(1) =  (1+v)*(1+w) 
         grad(2) =  (1+u)*(1+w)
         grad(3) =  (1+u)*(1+v)
      CASE (8)
         grad(1) = -(1+v)*(1+w) 
         grad(2) =  (1-u)*(1+w)
         grad(3) =  (1-u)*(1+v)
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickNodalPBasis','Unknown node for brick')
      END SELECT
      grad = grad/8
    END FUNCTION dBrickNodalPBasis


    FUNCTION ddBrickNodalPBasis(node, u, v, w) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(Kind=dp), INTENT(IN) :: u,v,w

      ! Variables
      REAL(Kind=dp) :: grad(3,3)

      grad = 0
      SELECT CASE(node)
      CASE (1)
        grad(1,2) =  (1-w)
        grad(1,3) =  (1-v)
        grad(2,3) =  (1-u)
      CASE (2)
        grad(1,2) = -(1-w) 
        grad(1,3) = -(1-v)
        grad(2,3) =  (1+u)
      CASE (3)
        grad(1,2) =  (1-w) 
        grad(1,3) = -(1+v)
        grad(2,3) = -(1+u)
      CASE (4)
        grad(1,2) = -(1-w) 
        grad(1,3) =  (1+v)
        grad(2,3) = -(1-u)
      CASE (5)
        grad(1,2) =  (1+w) 
        grad(1,3) = -(1-v)
        grad(2,3) = -(1-u)
      CASE (6)
        grad(1,2) = -(1+w) 
        grad(1,3) =  (1-v)
        grad(2,3) = -(1+u)
      CASE (7)
        grad(1,2) =  (1+w) 
        grad(1,3) =  (1+v)
        grad(2,3) =  (1+u)
      CASE (8)
        grad(1,2) = -(1+w) 
        grad(1,3) = -(1+v)
        grad(2,3) =  (1-u)
      CASE DEFAULT
        CALL Fatal('PElementBase::dBrickNodalPBasis','Unknown node for brick')
      END SELECT
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
      grad = grad / 8
    END FUNCTION ddBrickNodalPBasis



    ! As previous except obtain all nodal values at once. 
    SUBROUTINE dBrickNodalPBasisAll(u, v, w, gradphi) 
      IMPLICIT NONE

      ! Parameters
      REAL(Kind=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: gradphi(:,:)

      REAL(Kind=dp), PARAMETER :: c = 1.0_dp/8.0_dp            
      INTEGER, PARAMETER :: usgn(8) = [-1,1,1,-1,-1,1,1,-1]
      INTEGER, PARAMETER :: vsgn(8) = [-1,-1,1,1,-1,-1,1,1]
      INTEGER, PARAMETER :: wsgn(8) = [-1,-1,-1,-1,1,1,1,1]
            
      gradphi(1:8,1) = c*(usgn)*(1+vsgn*v)*(1+wsgn*w)
      gradphi(1:8,2) = c*(1+usgn*u)*(vsgn)*(1+wsgn*w)
      gradphi(1:8,3) = c*(1+usgn*u)*(1+vsgn*v)*(wsgn)

    END SUBROUTINE dBrickNodalPBasisAll


! --- start serendipity brick

!------------------------------------------------------------------------------
!>     Brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION SD_BrickEdgePBasis(edge, i , u, v, w, invertEdge) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of bricks edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks edge function i at point (u,v,w), i.e.
!       value = N_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      LOGICAL, OPTIONAL :: invertEdge
      
      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp) :: phipar, value

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      ! For inverted edges parameter inside phi function changes mark
      SELECT CASE(edge)
      ! Xi
      CASE (1,3,5,7)
         phiPar = u
      ! Eta
      CASE (2,4,6,8)
         phiPar = v
      ! Zeta
      CASE (9,10,11,12)
         phiPar = w
      END SELECT

      IF (invert) THEN
         phiPar = -phiPar
      END IF

      value = 0
      SELECT CASE(edge)
      CASE (1)
         value = 1d0/4*Phi(i,phiPar)*(1-v)*(1-w)
      CASE (2)
         value = 1d0/4*Phi(i,phiPar)*(1+u)*(1-w)
      CASE (3)
         value = 1d0/4*Phi(i,phiPar)*(1+v)*(1-w)
      CASE (4)
         value = 1d0/4*Phi(i,phiPar)*(1-u)*(1-w)
      CASE (5)
         value = 1d0/4*Phi(i,phiPar)*(1-v)*(1+w)
      CASE (6)
         value = 1d0/4*Phi(i,phiPar)*(1+u)*(1+w)
      CASE (7)
         value = 1d0/4*Phi(i,phiPar)*(1+v)*(1+w)
      CASE (8)
         value = 1d0/4*Phi(i,phiPar)*(1-u)*(1+w)
      CASE (9)
         value = 1d0/4*Phi(i,phiPar)*(1-u)*(1-v)
      CASE (10)
         value = 1d0/4*Phi(i,phiPar)*(1+u)*(1-v)
      CASE (11)
         value = 1d0/4*Phi(i,phiPar)*(1+u)*(1+v)
      CASE (12)
         value = 1d0/4*Phi(i,phiPar)*(1-u)*(1+v)
      CASE DEFAULT
         CALL Fatal('PElementBase::BrickEdgePBasis','Unknown edge for brick')
      END SELECT 
    END FUNCTION SD_BrickEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION SD_dBrickEdgePBasis(edge, i , u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of bricks edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp) :: phiU, phiV, phiW, phiPar 
      REAL (KIND=dp), DIMENSION(3) :: grad
      
      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      ! For inverted edges parameter inside phi function changes mark
      SELECT CASE(edge)
      ! Xi
      CASE (1,3,5,7)
         phiPar = u
      ! Eta
      CASE (2,4,6,8)
         phiPar = v
      ! Zeta
      CASE (9,10,11,12)
         phiPar = w
      END SELECT

      IF (invert) THEN
         phiPar = -phiPar
      END IF

      grad = 0
      SELECT CASE(edge)
      CASE (1)
         phiU = Phi(i,phiPar)
         grad(1) = 1d0/4*dPhi(i,phiPar)*(1-v)*(1-w)
         grad(2) = -1d0/4*phiU*(1-w)
         grad(3) = -1d0/4*phiU*(1-v)
      CASE (2)
         phiV = Phi(i,phiPar)
         grad(1) = 1d0/4*phiV*(1-w)
         grad(2) = 1d0/4*dPhi(i,phiPar)*(1+u)*(1-w)
         grad(3) = -1d0/4*phiV*(1+u)
      CASE (3)
         phiU = Phi(i,phiPar)
         grad(1) = 1d0/4*dPhi(i,phiPar)*(1+v)*(1-w)
         grad(2) = 1d0/4*phiU*(1-w)
         grad(3) = -1d0/4*phiU*(1+v)
      CASE (4)
         phiV = Phi(i,phiPar)
         grad(1) = -1d0/4*phiV*(1-w)
         grad(2) = 1d0/4*dPhi(i,phiPar)*(1-u)*(1-w)
         grad(3) = -1d0/4*phiV*(1-u)
      CASE (5)
         phiU = Phi(i,phiPar)
         grad(1) = 1d0/4*dPhi(i,phiPar)*(1-v)*(1+w)
         grad(2) = -1d0/4*phiU*(1+w)
         grad(3) = 1d0/4*phiU*(1-v)
      CASE (6)
         phiV = Phi(i,phiPar)
         grad(1) = 1d0/4*phiV*(1+w)
         grad(2) = 1d0/4*dPhi(i,phiPar)*(1+u)*(1+w)
         grad(3) = 1d0/4*phiV*(1+u)
      CASE (7)
         phiU = Phi(i,phiPar)
         grad(1) = 1d0/4*dPhi(i,phiPar)*(1+v)*(1+w)
         grad(2) = 1d0/4*phiU*(1+w)
         grad(3) = 1d0/4*phiU*(1+v)
      CASE (8)
         phiV = Phi(i,phiPar)
         grad(1) = -1d0/4*phiV*(1+w)
         grad(2) = 1d0/4*dPhi(i,phiPar)*(1-u)*(1+w)
         grad(3) = 1d0/4*phiV*(1-u)
      CASE (9)
         phiW = Phi(i,phiPar)
         grad(1) = -1d0/4*phiW*(1-v)
         grad(2) = -1d0/4*phiW*(1-u)
         grad(3) = 1d0/4*dPhi(i,phiPar)*(1-u)*(1-v)
      CASE (10)
         phiW = Phi(i,phiPar)
         grad(1) = 1d0/4*phiW*(1-v)
         grad(2) = -1d0/4*phiW*(1+u)
         grad(3) = 1d0/4*dPhi(i,phiPar)*(1+u)*(1-v)
      CASE (11)
         phiW = Phi(i,phiPar)
         grad(1) = 1d0/4*phiW*(1+v)
         grad(2) = 1d0/4*phiW*(1+u)
         grad(3) = 1d0/4*dPhi(i,phiPar)*(1+u)*(1+v)
      CASE (12)
         phiW = Phi(i,phiPar)
         grad(1) = -1d0/4*phiW*(1+v)
         grad(2) = 1d0/4*phiW*(1-u)
         grad(3) = 1d0/4*dPhi(i,phiPar)*(1-u)*(1+v)
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickEdgePBasis','Unknown edge for brick')
      END SELECT 

      ! Finally add derivative of dPhi s inner function to gradient
      ! if edge was inverted (=multiply by -1)
      IF (invert) THEN
         SELECT CASE(edge)
         ! Xi
         CASE (1,3,5,7)
            grad(1) = -grad(1)
         ! Eta
         CASE (2,4,6,8)
            grad(2) = -grad(2)
         ! Zeta
         CASE (9,10,11,12)
            grad(3) = -grad(3)
         END SELECT
      END IF
    END FUNCTION SD_dBrickEdgePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION SD_ddBrickEdgePBasis(edge, i , u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of bricks edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3,3)
!       gradient of bricks edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp) :: phiU, phiV, phiW, phiPar 
      REAL (KIND=dp), DIMENSION(3,3) :: grad
      
      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      ! For inverted edges parameter inside phi function changes mark
      SELECT CASE(edge)
      ! Xi
      CASE (1,3,5,7)
         phiPar = u
      ! Eta
      CASE (2,4,6,8)
         phiPar = v
      ! Zeta
      CASE (9,10,11,12)
         phiPar = w
      END SELECT
      IF (invert) phiPar = -phiPar

      grad = 0
      SELECT CASE(edge)
      CASE (1)
         grad(1,1) = ddPhi(i,phiPar)*(1-v)*(1-w)
         grad(1,2) = -dPhi(i,phiPar)*(1-w)
         grad(1,3) = -dPhi(i,phiPar)*(1-v)
         grad(2,3) = Phi(i,phiPar)
      CASE (2)
         grad(1,2) = dPhi(i,phiPar)*(1-w)
         grad(1,3) = -Phi(i,phiPar)
         grad(2,2) = ddPhi(i,phiPar)*(1+u)*(1-w)
         grad(2,3) = -dPhi(i,phiPar)*(1+u)
      CASE (3)
         grad(1,1) =  ddPhi(i,phiPar)*(1+v)*(1-w)
         grad(1,2) =  dPhi(i,phiPar)*(1-w)
         grad(1,3) = -dPhi(i,phiPar)*(1+v)
         grad(2,3) = -Phi(i,phiPar)
      CASE (4)
         grad(1,2) = -dPhi(i,phiPar)*(1-w)
         grad(1,3) = Phi(i,phiPar)
         grad(2,2) = ddPhi(i,phiPar)*(1-u)*(1-w)
         grad(2,3) = -dPhi(i,phiPar)*(1-u)
      CASE (5)
         grad(1,1) = ddPhi(i,phiPar)*(1-v)*(1+w)
         grad(1,2) = -dPhi(i,phiPar)*(1+w)
         grad(1,3) =  dPhi(i,phiPar)*(1-v)
         grad(2,3) = -Phi(i,phiPar)
      CASE (6)
         grad(1,2) = dPhi(i,phiPar)*(1+w)
         grad(1,3) = Phi(i,phiPar)
         grad(2,2) = ddPhi(i,phiPar)*(1+u)*(1+w)
         grad(2,3) = dPhi(i,phiPar)*(1+u)
      CASE (7)
         grad(1,1) = ddPhi(i,phiPar)*(1+v)*(1+w)
         grad(1,2) = dPhi(i,phiPar)*(1+w)
         grad(1,3) = dPhi(i,phiPar)*(1+v)
         grad(2,3) = Phi(i,phiPar)
      CASE (8)
         grad(1,2) = -dPhi(i,phiPar)*(1+w)
         grad(1,3) = -Phi(i,phiPar)
         grad(2,2) = ddPhi(i,phiPar)*(1-u)*(1+w)
         grad(2,3) = dPhi(i,phiPar)*(1-u)
      CASE (9)
         grad(1,2) = Phi(i,phiPar)
         grad(1,3) = -dPhi(i,phiPar)*(1-v)
         grad(2,3) = -dPhi(i,phiPar)*(1-u)
         grad(3,3) = ddPhi(i,phiPar)*(1-u)*(1-v)
      CASE (10)
         grad(1,2) = -Phi(i,phiPar)
         grad(1,3) = dPhi(i,phiPar)*(1-v)
         grad(2,3) = -dPhi(i,phiPar)*(1+u)
         grad(3,3) = ddPhi(i,phiPar)*(1+u)*(1-v)
      CASE (11)
         grad(1,2) = Phi(i,phiPar)
         grad(1,3) = dPhi(i,phiPar)*(1+v)
         grad(2,3) = dPhi(i,phiPar)*(1+u)
         grad(3,3) = ddPhi(i,phiPar)*(1+u)*(1+v)
      CASE (12)
         grad(1,2) = -Phi(i,phiPar)
         grad(1,3) = -dPhi(i,phiPar)*(1+v)
         grad(2,3) = dPhi(i,phiPar)*(1-u)
         grad(3,3) = ddPhi(i,phiPar)*(1-u)*(1+v)
      CASE DEFAULT
         CALL Fatal('PElementBase::ddBrickEdgePBasis','Unknown edge for brick')
      END SELECT 

      ! Finally add derivative of dPhi s inner function to gradient
      ! if edge was inverted (=multiply by -1)
      IF (invert) THEN
         SELECT CASE(edge)
         ! Xi
         CASE (1,3,5,7)
            grad(1,2) = -grad(1,2)
            grad(1,3) = -grad(1,3)
         ! Eta
         CASE (2,4,6,8)
            grad(1,2) = -grad(1,2)
            grad(2,3) = -grad(2,3)
         ! Zeta
         CASE (9,10,11,12)
            grad(1,3) = -grad(1,3)
            grad(2,3) = -grad(2,3)
         END SELECT
      END IF
      grad = grad / 4
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION SD_ddBrickEdgePBasis

!------------------------------------------------------------------------------    
!>     Brick face basis at point (u,v,w)
!------------------------------------------------------------------------------    
    FUNCTION SD_BrickFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(value)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of bricks face function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square face to define direction of face
!        function. Default numbering is that defined for face in PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks face function m(i,j) at point (u,v,w), i.e.
!       value = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      INTEGER, DIMENSION(4) :: local
      REAL (KIND=dp) :: La, Lb, Lc, Lh, value

      ! If local numbering not present use default numbering
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getBrickFaceMap(face)
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set parameters for face value calculation
      La = BrickL(local(1),u,v,w)
      Lb = BrickL(local(2),u,v,w)
      Lc = BrickL(local(4),u,v,w)

      SELECT CASE(face)
      CASE (1)
         Lh = (1-w)
      CASE (2)
         Lh = (1+w)
      CASE (3)
         Lh = (1-v)
      CASE (4)
         Lh = (1+u)
      CASE (5)
         Lh = (1+v)
      CASE (6)
         Lh = (1-u)
      CASE DEFAULT
         CALL Fatal('PElementBase::BrickFacePBasis','Unknown face for brick')
      END SELECT

      ! Calculate value of function from general form
      value = 1d0/2*Lh*Phi(i,Lb-La)*Phi(j,Lc-La)
    END FUNCTION SD_BrickFacePBasis


!------------------------------------------------------------------------------
!>     Gradient of brick face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_dBrickFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of bricks face function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square face to define direction of face
!        function. Default numbering is that defined for face in PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers

      ! Variables
      INTEGER, DIMENSION(4) :: local
      REAL (KIND=dp) :: La, Lb, Lc, Lh, phiI, phiJ
      REAL (KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLh, grad
      
      ! If local numbering not present use default numbering
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getBrickFaceMap(face)
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set parameters for face calculation
      La = BrickL(local(1),u,v,w)
      Lb = BrickL(local(2),u,v,w)
      Lc = BrickL(local(4),u,v,w)
      dLa = dBrickL(local(1),u,v,w)
      dLb = dBrickL(local(2),u,v,w)
      dLc = dBrickL(local(4),u,v,w)

      SELECT CASE(face)
      CASE(1)
         Lh = (1-w)
         dLh = [ 0d0,0d0,-1d0 ]
      CASE(2)
         Lh = (1+w)
         dLh = [ 0d0,0d0,1d0 ]
      CASE(3)
         Lh = (1-v)
         dLh = [ 0d0,-1d0,0d0 ]
      CASE(4)
         Lh = (1+u)
         dLh = [ 1d0, 0d0, 0d0]
      CASE(5)
         Lh = (1+v)
         dLh = [ 0d0,1d0,0d0 ]
      CASE(6)
         Lh = (1-u)
         dLh = [ -1d0,0d0,0d0 ]
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickFacePBasis','Unknown face for brick')
      END SELECT

      ! Calculate value of gradient from general form
      grad = 0
      phiI = Phi(i,Lb-La)
      phiJ = Phi(j,Lc-La)

      grad = 1d0/2*(dLh*phiI*phiJ+ Lh*dPhi(i,Lb-La)*(dLb-dLa)*phiJ + &
           Lh*phiI*dPhi(j,Lc-La)*(dLc-dLa))
    END FUNCTION SD_dBrickFacePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of brick face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_ddBrickFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of bricks face function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square face to define direction of face
!        function. Default numbering is that defined for face in PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers

      ! Variables
      INTEGER, DIMENSION(4) :: local
      REAL (KIND=dp) :: La, Lb, Lc, Lh, phiI, phiJ, grad(3,3)
      REAL (KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLh

      ! If local numbering not present use default numbering
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getBrickFaceMap(face)
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set parameters for face calculation
      La = BrickL(local(1),u,v,w)
      Lb = BrickL(local(2),u,v,w)
      Lc = BrickL(local(4),u,v,w)

      dLa = dBrickL(local(1),u,v,w)
      dLb = dBrickL(local(2),u,v,w)
      dLc = dBrickL(local(4),u,v,w)

      SELECT CASE(face)
      CASE(1)
         Lh = (1-w)
         dLh = [ 0d0,0d0,-1d0 ]
      CASE(2)
         Lh = (1+w)
         dLh = [ 0d0,0d0,1d0 ]
      CASE(3)
         Lh = (1-v)
         dLh = [ 0d0,-1d0,0d0 ]
      CASE(4)
         Lh = (1+u)
         dLh = [ 1d0, 0d0, 0d0]
      CASE(5)
         Lh = (1+v)
         dLh = [ 0d0,1d0,0d0 ]
      CASE(6)
         Lh = (1-u)
         dLh = [ -1d0,0d0,0d0 ]
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickFacePBasis','Unknown face for brick')
      END SELECT

      ! Calculate value of gradient from general form
      grad = 0
      phiI = Phi(i,Lb-La)
      phiJ = Phi(j,Lc-La)

!     grad = (dLh*phiI*phiJ + Lh*dPhi(i,Lb-La)*(dLb-dLa)*phiJ + &
!                  Lh*phiI*dPhi(j,Lc-La)*(dLc-dLa)) / 2

      grad(1,1) = grad(1,1) + dLh(1)*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*Phi(j,Lc-La)
      grad(1,1) = grad(1,1) + dLh(1)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,1) = grad(1,1) + dLh(1)*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*Phi(j,Lc-La)
      grad(1,1) = grad(1,1) + Lh*ddPhi(i,Lb-La)*(dLb(1)-dLa(1))**2*Phi(j,Lc-La)
      grad(1,1) = grad(1,1) + Lh*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,1) = grad(1,1) + dLh(1)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,1) = grad(1,1) + Lh*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,1) = grad(1,1) + Lh*Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(1)-dLa(1))**2

      grad(1,2) = grad(1,2) + dLh(1)*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*Phi(j,Lc-La)
      grad(1,2) = grad(1,2) + dLh(1)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(1,2) = grad(1,2) + dLh(2)*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*Phi(j,Lc-La)
      grad(1,2) = grad(1,2) + Lh*ddPhi(i,Lb-La)*(dLb(1)-dLa(1))*(dLb(2)-dLa(2))*Phi(j,Lc-La)
      grad(1,2) = grad(1,2) + Lh*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(1,2) = grad(1,2) + dLh(2)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,2) = grad(1,2) + Lh*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,2) = grad(1,2) + Lh*Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(1)-dLa(1))*(dLc(2)-dLa(2))

      grad(1,3) = grad(1,3) + dLh(1)*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*Phi(j,Lc-La)
      grad(1,3) = grad(1,3) + dLh(1)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(3)-dLa(3))
      grad(1,3) = grad(1,3) + dLh(3)*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*Phi(j,Lc-La)
      grad(1,3) = grad(1,3) + Lh*ddPhi(i,Lb-La)*(dLb(1)-dLa(1))*(dLb(3)-dLa(3))*Phi(j,Lc-La)
      grad(1,3) = grad(1,3) + Lh*dPhi(i,Lb-La)*(dLb(1)-dLa(1))*dPhi(j,Lc-La)*(dLC(3)-dLa(3))
      grad(1,3) = grad(1,3) + dLh(3)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,3) = grad(1,3) + Lh*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*dPhi(j,Lc-La)*(dLc(1)-dLa(1))
      grad(1,3) = grad(1,3) + Lh*Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(1)-dLa(1))*(dLc(3)-dLa(3))

      grad(2,2) = grad(2,2) + dLh(2)*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*Phi(j,Lc-La)
      grad(2,2) = grad(2,2) + dLh(2)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,2) = grad(2,2) + dLh(2)*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*Phi(j,Lc-La)
      grad(2,2) = grad(2,2) + Lh*ddPhi(i,Lb-La)*(dLb(2)-dLa(2))**2*Phi(j,Lc-La)
      grad(2,2) = grad(2,2) + Lh*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLC(2)-dLa(2))
      grad(2,2) = grad(2,2) + dLh(2)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,2) = grad(2,2) + Lh*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,2) = grad(2,2) + Lh*Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(2)-dLa(2))**2

      grad(2,3) = grad(2,3) + dLh(2)*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*Phi(j,Lc-La)
      grad(2,3) = grad(2,3) + dLh(2)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(3)-dLa(3))
      grad(2,3) = grad(2,3) + dLh(3)*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*Phi(j,Lc-La)
      grad(2,3) = grad(2,3) + Lh*ddPhi(i,Lb-La)*(dLb(2)-dLa(2))*(dLb(3)-dLa(3))*Phi(j,Lc-La)
      grad(2,3) = grad(2,3) + Lh*dPhi(i,Lb-La)*(dLb(2)-dLa(2))*dPhi(j,Lc-La)*(dLC(3)-dLa(3))
      grad(2,3) = grad(2,3) + dLh(3)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,3) = grad(2,3) + Lh*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*dPhi(j,Lc-La)*(dLc(2)-dLa(2))
      grad(2,3) = grad(2,3) + Lh*Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(2)-dLa(2))*(dLc(3)-dLa(3))

      grad(3,3) = grad(3,3) + dLh(3)*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*Phi(j,Lc-La)
      grad(3,3) = grad(3,3) + dLh(3)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(3)-dLa(3))
      grad(3,3) = grad(3,3) + dLh(3)*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*Phi(j,Lc-La)
      grad(3,3) = grad(3,3) + Lh*ddPhi(i,Lb-La)*(dLb(3)-dLa(3))**2*Phi(j,Lc-La)
      grad(3,3) = grad(3,3) + Lh*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*dPhi(j,Lc-La)*(dLC(3)-dLa(3))
      grad(3,3) = grad(3,3) + dLh(3)*Phi(i,Lb-La)*dPhi(j,Lc-La)*(dLc(3)-dLa(3))
      grad(3,3) = grad(3,3) + Lh*dPhi(i,Lb-La)*(dLb(3)-dLa(3))*dPhi(j,Lc-La)*(dLc(3)-dLa(3))
      grad(3,3) = grad(3,3) + Lh*Phi(i,Lb-La)*ddPhi(j,Lc-La)*(dLc(3)-dLa(3))**2

      grad = grad / 2
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION SD_ddBrickFacePBasis


!------------------------------------------------------------------------------
!>    Brick bubble basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION SD_BrickBubblePBasis(i, j, k, u, v, w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(2,2,2),(2,2,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      ! Result
      REAL (KIND=dp) :: value

      value = Phi(i,u)*Phi(j,v)*Phi(k,w)
    END FUNCTION SD_BrickBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of brick bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_dBrickBubblePBasis(i, j, k, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(2,2,2),(2,2,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      ! Variables
      REAL (KIND=dp) :: phiU, phiV, phiW
      REAL (KIND=dp), DIMENSION(3) :: grad

      grad = 0
      phiU = Phi(i,u) 
      phiV = Phi(j,v)
      phiW = Phi(k,w)
      grad(1) = dPhi(i,u)*phiV*phiW
      grad(2) = phiU*dPhi(j,v)*phiW
      grad(3) = phiU*phiV*dPhi(k,w)
    END FUNCTION SD_dBrickBubblePBasis

!------------------------------------------------------------------------------
!>    2nd derivatives of brick bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_ddBrickBubblePBasis(i, j, k, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(2,2,2),(2,2,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      ! Variables
      REAL (KIND=dp) :: phiU, phiV, phiW
      REAL (KIND=dp), DIMENSION(3,3) :: grad

      grad = 0
      phiU = Phi(i,u) 
      phiV = Phi(j,v)
      phiW = Phi(k,w)

      grad(1,1) = ddPhi(i,u)*PhiV*phiW
      grad(1,2) = dPhi(i,u)*dPhi(j,v)*phiW
      grad(1,3) = dPhi(i,u)*PhiV*dPhi(k,w)

      grad(2,2) = PhiU*ddPhi(j,v)*phiW
      grad(2,3) = PhiU*dPhi(j,v)*dPhi(k,w)

      grad(3,3) = PhiU*PhiV*ddPhi(k,w)

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION SD_ddBrickBubblePBasis

! --- end serendipity brick

          
!------------------------------------------------------------------------------
!>     Brick edge basis at point (u,v,w). Compatible with pyramidal edge basis.
!------------------------------------------------------------------------------
    FUNCTION BrickEdgePBasis(edge, i , u, v, w, invertEdge) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of bricks edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks edge function i at point (u,v,w), i.e.
!       value = N_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u, v, w
      LOGICAL, OPTIONAL :: invertEdge
      
      ! Variables
      LOGICAL :: invert
      INTEGER :: local(2)
      REAL(KIND=dp) :: Pa, Pb, La, Lb, phiPar, value

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 12) THEN
         CALL Fatal('PElementBase::BrickPyraEdgePBasis','Unknown edge for brick')
      END IF

      ! Get nodes of edge
      local(1:2) = getBrickEdgeMap(edge)

      ! Bilinear nodal functions
      Pa = BrickNodalPBasis(local(1),u,v,w)
      Pb = BrickNodalPBasis(local(2),u,v,w)

      ! Affine functions for edge direction
      La = BrickL(local(1),u,v,w)
      Lb = BrickL(local(2),u,v,w)

      PhiPar = Lb-La
      IF ( invert ) PhiPar = -PhiPar

      ! Get value of edge function
      value = Pa*Pb*varPhi(i,PhiPar)
    END FUNCTION BrickEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION dBrickEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of bricks edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      LOGICAL, OPTIONAL :: invertEdge

      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp) :: Pa,Pb,La,Lb, vPhi, PhiPar, dVPhi(3)
      REAL (KIND=dp), DIMENSION(3) :: dPa, dPb, dLa, dLb, grad,dPhiPar
      INTEGER :: local(2), swap

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 12) THEN
         CALL Fatal('PElementBase::dBrickPyraEdgePBasis','Unknown edge for brick')
      END IF

      ! Get local of edge
      local(1:2) = getBrickEdgeMap(edge)

      ! Trilinear nodal functions and their derivatives
      Pa  = BrickNodalPBasis(local(1),u,v,w)
      Pb  = BrickNodalPBasis(local(2),u,v,w)

      dPa = dBrickNodalPBasis(local(1),u,v,w)
      dPb = dBrickNodalPBasis(local(2),u,v,w)

      ! Affine functions and their derivatives for edge direction
      La  = BrickL(local(1),u,v,w)
      Lb  = BrickL(local(2),u,v,w)

      dLa = dBrickL(local(1),u,v,w)
      dLb = dBrickL(local(2),u,v,w)

      PhiPar = Lb-La
      dPhiPar = dLb-dLa

      IF ( invert ) THEN
        PhiPar = -PhiPar
        dPhiPar = -dPhiPar
      END IF

      ! Get value of edge function
      VPhi  = VarPhi(i,Phipar)
      dVPhi = dVarPhi(i,Phipar)*dPhiPar
      grad = dPa*Pb*vPhi + Pa*dPb*vPhi + Pa*Pb*dVPhi
    END FUNCTION dBrickEdgePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION ddBrickEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of bricks edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      LOGICAL, OPTIONAL :: invertEdge

      REAL(KIND=dp) :: grad(3,3)

      ! Variables
      LOGICAL :: invert
      REAL (KIND=dp) :: s,Pa,Pb,La,Lb, vPhi, PhiPar, ddPa(3,3), ddPb(3,3)
      REAL (KIND=dp) :: dPa(3), dPb(3), dLa(3), dLb(3), dPhiPar(3), dvPhi(3), ddvPhi(3,3)
      INTEGER :: local(2),p,q

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 12) THEN
         CALL Fatal('PElementBase::dBrickPyraEdgePBasis','Unknown edge for brick')
      END IF

      ! Get local of edge
      local = getBrickEdgeMap(edge)

      ! Trilinear nodal functions and their derivatives
      Pa   = BrickNodalPBasis(local(1),u,v,w)
      Pb   = BrickNodalPBasis(local(2),u,v,w)

      dPa  = dBrickNodalPBasis(local(1),u,v,w)
      dPb  = dBrickNodalPBasis(local(2),u,v,w)

      ddPa = ddBrickNodalPBasis(local(1),u,v,w)
      ddPb = ddBrickNodalPBasis(local(2),u,v,w)

      ! Affine functions and their derivatives for edge direction
      La  = BrickL(local(1),u,v,w)
      Lb  = BrickL(local(2),u,v,w)

      dLa = dBrickL(local(1),u,v,w)
      dLb = dBrickL(local(2),u,v,w)

      PhiPar = Lb-La
      dPhiPar = dLb-dLa

      IF ( invert ) THEN
        PhiPar = -PhiPar
        dPhiPar = -dPhiPar
      END IF

      ! Get value of edge function
      vPhi  = varPhi(i,PhiPar)
      dvPhi = dVarPhi(i,PhiPar)*dPhiPar
      s = ddVarPhi(i,PhiPar)
      ddVphi = 0
      DO p=1,3
        DO q=p,3
          ddVPhi(p,q) = s*dPhiPar(p)*dPhiPar(q)
        END DO
      END DO

!     grad(1) = dPa(1)*Pb*vPhi + Pa*dPb(1)*vPhi + Pa*Pb*dvPhi(1)
#if 1
      grad = 0
      grad(1,1) = grad(1,1) + ddPa(1,1)*Pb*vPhi + dPa(1)*dPb(1)*vPhi + dPa(1)*Pb*dvPhi(1)
      grad(1,1) = grad(1,1) + dPa(1)*dPb(1)*vPhi + Pa*ddPb(1,1)*vPhi + Pa*dPb(1)*dvPhi(1)
      grad(1,1) = grad(1,1) + dPa(1)*Pb*dvPhi(1) + Pa*dPb(1)*dvPhi(1) + Pa*Pb*ddvPhi(1,1)

      grad(1,2) = grad(1,2) + ddPa(1,2)*Pb*vPhi + dPa(2)*dPb(1)*vPhi + dPa(2)*Pb*dvPhi(1)
      grad(1,2) = grad(1,2) + dPa(1)*dPb(2)*vPhi + Pa*ddPb(1,2)*vPhi + Pa*dPb(2)*dvPhi(1)
      grad(1,2) = grad(1,2) + dPa(1)*Pb*dvPhi(2) + Pa*dPb(1)*dvPhi(2) + Pa*Pb*ddvPhi(1,2)

      grad(1,3) = grad(1,3) + ddPa(1,3)*Pb*vPhi + dPa(3)*dPb(1)*vPhi + dPa(3)*Pb*dvPhi(1)
      grad(1,3) = grad(1,3) + dPa(1)*dPb(3)*vPhi + Pa*ddPb(1,3)*vPhi + Pa*dPb(3)*dvPhi(1)
      grad(1,3) = grad(1,3) + dPa(1)*Pb*dvPhi(3) + Pa*dPb(1)*dvPhi(3) + Pa*Pb*ddvPhi(1,3)

      grad(2,2) = grad(2,2) + ddPa(2,2)*Pb*vPhi + dPa(2)*dPb(2)*vPhi + dPa(2)*Pb*dvPhi(2)
      grad(2,2) = grad(2,2) + dPa(2)*dPb(2)*vPhi + Pa*ddPb(2,2)*vPhi + Pa*dPb(2)*dvPhi(2)
      grad(2,2) = grad(2,2) + dPa(2)*Pb*dvPhi(2) + Pa*dPb(2)*dvPhi(2) + Pa*Pb*ddvPhi(2,2)

      grad(2,3) = grad(2,3) + ddPa(2,3)*Pb*vPhi + dPa(3)*dPb(2)*vPhi + dPa(3)*Pb*dvPhi(2)
      grad(2,3) = grad(2,3) + dPa(2)*dPb(3)*vPhi + Pa*ddPb(2,3)*vPhi + Pa*dPb(3)*dvPhi(2)
      grad(2,3) = grad(2,3) + dPa(2)*Pb*dvPhi(3) + Pa*dPb(2)*dvPhi(3) + Pa*Pb*ddvPhi(2,3)

      grad(3,3) = grad(3,3) + ddPa(3,3)*Pb*vPhi + dPa(3)*dPb(3)*vPhi + dPa(3)*Pb*dvPhi(3)
      grad(3,3) = grad(3,3) + dPa(3)*dPb(3)*vPhi + Pa*ddPb(3,3)*vPhi + Pa*dPb(3)*dvPhi(3)
      grad(3,3) = grad(3,3) + dPa(3)*Pb*dvPhi(3) + Pa*dPb(3)*dvPhi(3) + Pa*Pb*ddvPhi(3,3)
#else
      BLOCK
        REAL(KIND=dp) :: f(3), df(3,3), ddf(3,3,3)

        f=[Pa, Pb, vPhi]
        df(1,:) = dPa
        df(2,:) = dPb
        df(3,:) = dVphi
        ddf(1,:,:)=ddPa; ddf(2,:,:)=ddPb; ddf(3,:,:)=ddVphi 

        grad = Product2ndDerivatives(3,f,df,ddf,3,0)
      END BLOCK
#endif

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddBrickEdgePBasis


!------------------------------------------------------------------------------    
!>     Brick face basis at point (u,v,w)
!------------------------------------------------------------------------------    
    FUNCTION BrickFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(value)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of bricks face function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square face to define direction of face
!        function. Default numbering is that defined for face in PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks face function m(i,j) at point (u,v,w), i.e.
!       value = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      INTEGER, DIMENSION(4) :: local
      REAL (KIND=dp) :: La, Lb, Lc, Lh, value, Pa,Pb

      ! If local numbering not present use default numbering
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getBrickFaceMap(face)
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set parameters for face value calculation
      La = BrickL(local(1),u,v,w)
      Lb = BrickL(local(2),u,v,w)
      Lc = BrickL(local(4),u,v,w)

      ! Calculate value of function from general form
      Pa = BrickNodalPBasis(local(1),u,v,w)
      Pb = BrickNodalPBasis(local(3),u,v,w)
      value = Pa*Pb*LegendreP(i,Lb-La)*LegendreP(j,Lc-La)
    END FUNCTION BrickFacePBasis


!------------------------------------------------------------------------------
!>     Gradient of brick face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dBrickFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of bricks face function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square face to define direction of face
!        function. Default numbering is that defined for face in PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers

      ! Variables
      INTEGER, DIMENSION(4) :: local
      REAL (KIND=dp) :: La, Lb, Lc, Lh, phiI, phiJ, Pa, Pb
      REAL (KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLh, grad, &
                dPa, dPb, dPhiI, dPhiJ
      
      ! If local numbering not present use default numbering
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getBrickFaceMap(face)
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set parameters for face calculation
      La = BrickL(local(1),u,v,w)
      Lb = BrickL(local(2),u,v,w)
      Lc = BrickL(local(4),u,v,w)

      dLa = dBrickL(local(1),u,v,w)
      dLb = dBrickL(local(2),u,v,w)
      dLc = dBrickL(local(4),u,v,w)

      ! Calculate value of gradient from general form
      Pa = BrickNodalPBasis(local(1),u,v,w)
      Pb = BrickNodalPBasis(local(3),u,v,w)

      dPa = dBrickNodalPBasis(local(1),u,v,w)
      dPb = dBrickNodalPBasis(local(3),u,v,w)

      PhiI = LegendreP(i,Lb-La)
      PhiJ = LegendreP(j,Lc-La)

      dPhiI = dLegendreP(i,Lb-La)*(dLb-dLa)
      dPhiJ = dLegendreP(j,Lc-La)*(dLc-dLa)

      grad = dPa*Pb*PhiI*PhiJ + Pa*dPb*PhiI*PhiJ + &
             Pa*Pb*dPhiI*PhiJ + Pa*Pb*PhiI*dPhiJ
    END FUNCTION dBrickFacePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of brick face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddBrickFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of bricks face function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square face to define direction of face
!        function. Default numbering is that defined for face in PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of bricks face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE 

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers

      REAL(KIND=dp) :: grad(3,3)

      ! Variables
      INTEGER :: local(4), p, q
      REAL (KIND=dp) :: La, Lb, Lc, Lh, phiI, phiJ, Pa, Pb, s, t, ddLegi, ddLegj, &
          dPa(3), dPb(3), ddPa(3,3), ddPb(3,3), dPhiI(3), dPhiJ(3), ddPhiI(3,3), ddPhiJ(3,3)
      REAL (KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLh, ds, dt

      ! If local numbering not present use default numbering
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getBrickFaceMap(face)
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set parameters for face calculation
      La  = BrickL(local(1),u,v,w)
      Lb  = BrickL(local(2),u,v,w)
      Lc  = BrickL(local(4),u,v,w)

      dLa = dBrickL(local(1),u,v,w)
      dLb = dBrickL(local(2),u,v,w)
      dLc = dBrickL(local(4),u,v,w)

      Pa = BrickNodalPBasis(local(1),u,v,w)
      Pb = BrickNodalPBasis(local(3),u,v,w)

      dPa = dBrickNodalPBasis(local(1),u,v,w)
      dPb = dBrickNodalPBasis(local(3),u,v,w)

      ddPa = ddBrickNodalPBasis(local(1),u,v,w)
      ddPb = ddBrickNodalPBasis(local(3),u,v,w)

      s = Lb-La
      ds = dLb-dLa

      t = Lc-La
      dt = dLc-dLa

      ! Calculate value of gradient from general form
      phiI = LegendreP(i,s)
      phiJ = LegendreP(j,t)

      dPhiI = dLegendreP(i,s)*ds
      dPhiJ = dLegendreP(j,t)*dt
      
      ddLegi = ddLegendreP(i,s)
      ddLegj = ddLegendreP(j,t)
      DO p=1,3
        DO q=p,3
          ddPhiI(p,q) = ddLegi*ds(p)*ds(q)
          ddPhiJ(p,q) = ddLegj*dt(p)*dt(q)
        END DO
      END DO

      BLOCK
        REAL(KIND=dp) :: f(4), df(4,3), ddf(4,3,3)

        f = [Pa,Pb,PhiI,PhiJ]
        df(1,:) = dPa
        df(2,:) = dPb
        df(3,:) = dPhiI
        df(4,:) = dPhiJ
        ddf(1,:,:)=ddPa; ddf(2,:,:)=ddPb; ddf(3,:,:)=ddPhiI; ddf(4,:,:)=ddPhiJ

        grad = Product2ndDerivatives(4,f,df,ddf,3,0)
      END BLOCK
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddBrickFacePBasis


!------------------------------------------------------------------------------
!>    Brick bubble basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION BrickBubblePBasis(i, j, k, u, v, w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(2,2,2),(2,2,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      ! Result
      REAL (KIND=dp) :: value

      value=Phi(i+2,u)*Phi(j+2,v)*Phi(k+2,w)
    END FUNCTION BrickBubblePBasis


!------------------------------------------------------------------------------
!>    Brick bubble basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION dBrickBubblePBasis(i, j, k, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(2,2,2),(2,2,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u, v, w
      ! Result
      REAL (KIND=dp) :: grad(3)

      grad(1) = dPhi(i+2,u)*Phi(j+2,v)*Phi(k+2,w)
      grad(2) = Phi(i+2,u)*dPhi(j+2,v)*Phi(k+2,w)
      grad(3) = Phi(i+2,u)*Phi(j+2,v)*dPhi(k+2,w)
    END FUNCTION dBrickBubblePBasis
!------------------------------------------------------------------------------


!>    Brick bubble basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION ddBrickBubblePBasis(i, j, k, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k)
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of bricks bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL(KIND=dp), INTENT(IN) :: u, v, w
      ! Result
      REAL(KIND=dp) :: grad(3,3)

      grad=0
      grad(1,1) = ddPhi(i+2,u)*Phi(j+2,v)*Phi(k+2,w)
      grad(1,2) = dPhi(i+2,u)*dPhi(j+2,v)*Phi(k+2,w)
      grad(1,3) = dPhi(i+2,u)*Phi(j+2,v)*dPhi(k+2,w)
      grad(2,2) = Phi(i+2,u)*ddPhi(j+2,v)*Phi(k+2,w)
      grad(2,3) = Phi(i+2,u)*dPhi(j+2,v)*dPhi(k+2,w)
      grad(3,3) = Phi(i+2,u)*Phi(j+2,v)*ddPhi(k+2,w)

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddBrickBubblePBasis


    PURE FUNCTION BrickL(which, u, v, w) RESULT(value)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v,w 
      ! Variables
      REAL(KIND=dp) :: value
      
      value = 0
      SELECT CASE(which)
      CASE (1)
         value = (3-u-v-w)/2d0
      CASE (2)
         value = (3+u-v-w)/2d0
      CASE (3)
         value = (3+u+v-w)/2d0
      CASE (4)
         value = (3-u+v-w)/2d0
      CASE (5)
         value = (3-u-v+w)/2d0
      CASE (6)
         value = (3+u-v+w)/2d0
      CASE (7)
         value = (3+u+v+w)/2d0
      CASE (8)   
         value = (3-u+v+w)/2d0
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::BrickL','Unknown function L for brick')
#endif
      END SELECT
    END FUNCTION BrickL

    PURE FUNCTION dBrickL(which, u, v, w) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v,w 
      ! Variables
      REAL(KIND=dp) :: grad(3)
      
      grad = 0
      SELECT CASE(which)
      CASE (1)
         grad(1) = -1d0/2 
         grad(2) = -1d0/2
         grad(3) = -1d0/2
      CASE (2)
         grad(1) =  1d0/2 
         grad(2) = -1d0/2
         grad(3) = -1d0/2
      CASE (3)
         grad(1) =  1d0/2 
         grad(2) =  1d0/2
         grad(3) = -1d0/2
      CASE (4)
         grad(1) = -1d0/2 
         grad(2) =  1d0/2
         grad(3) = -1d0/2
      CASE (5)
         grad(1) = -1d0/2 
         grad(2) = -1d0/2
         grad(3) =  1d0/2
      CASE (6)
         grad(1) =  1d0/2 
         grad(2) = -1d0/2
         grad(3) =  1d0/2
      CASE (7)
         grad(1) =  1d0/2 
         grad(2) =  1d0/2
         grad(3) =  1d0/2
      CASE (8)   
         grad(1) = -1d0/2 
         grad(2) =  1d0/2
         grad(3) =  1d0/2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickL','Unknown function dL for brick')
#endif
      END SELECT
    END FUNCTION dBrickL

    
!------------------------------------------------------------------------------
!>     Tetrahedron nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION TetraNodalPBasis(node, u, v, w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of tetrahedral nodal function to calculate
!        node = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of tetras nodal function at point (u,v,w), i.e.
!       value = N_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Result
      REAL (KIND=dp) :: value

      value = 0
      SELECT CASE(node)
         CASE(1)
            value = 1d0/2*(1-u-v/SQRT(3d0)-w/SQRT(6d0))
         CASE(2)
            value = 1d0/2*(1+u-v/SQRT(3d0)-w/SQRT(6d0))
         CASE(3)
            value = SQRT(3d0)/3*(v-w/SQRT(8d0))
         CASE(4)
            value = SQRT(3d0/8d0)*w
         CASE DEFAULT 
            CALL Fatal('PElementBase::TetraNodalPBasis','Unknown node for tetrahedron')
         END SELECT
    END FUNCTION TetraNodalPBasis

    
    SUBROUTINE TetraNodalPBasisAll(u, v, w, phi) 
      IMPLICIT NONE
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      REAL (KIND=dp) :: phi(:)
      REAL(KIND=dp), PARAMETER :: half = 1.0_dp/2.0_dp, &
          c3 = 1.0_dp/SQRT(3.0_dp), c6 = 1.0_dp/SQRT(6.0_dp), c8 = SQRT(3.0_dp/8.0_dp)

      phi(1) = half*(1-u-c3*v-c6*w)
      phi(2) = half*(1+u-c3*v-c6*w)
      phi(3) = c3*v - half*c6*w 
      phi(4) = c8*w
    END SUBROUTINE TetraNodalPBasisAll

    SUBROUTINE TetraNodalLBasisAll(u, v, w, phi) 
      IMPLICIT NONE
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      REAL (KIND=dp) :: phi(:)
      phi(1) = 1.0_dp-u-v-w
      phi(2) = u
      phi(3) = v
      phi(4) = w
    END SUBROUTINE TetraNodalLBasisAll

    
    
!------------------------------------------------------------------------------
!>     Gradient of tetrahedrons nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dTetraNodalPBasis(node, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of tetrahedral nodal function to calculate
!        node = {1,2,3,4}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras nodal function at point (u,v,w), i.e.
!       grad = dN_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: node
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Result
      REAL (KIND=dp), DIMENSION(3) :: grad

      grad = 0
      SELECT CASE(node)
      CASE(1)
         grad(1)=-1d0/2
         grad(2)=-SQRT(3d0)/6
         grad(3)=-SQRT(6d0)/12
      CASE(2)
         grad(1)=1d0/2
         grad(2)=-SQRT(3d0)/6
         grad(3)=-SQRT(6d0)/12
      CASE(3)
         grad(1)=0
         grad(2)=SQRT(3d0)/3
         grad(3)=-SQRT(6d0)/12
      CASE(4)
         grad(1)=0
         grad(2)=0
         grad(3)=SQRT(6d0)/4
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraNodalPBasis','Unknown node for tetrahedron')
      END SELECT
    END FUNCTION dTetraNodalPBasis


    SUBROUTINE dTetraNodalPBasisAll(u, v, w, gradphi )
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: gradphi(:,:)
      REAL(KIND=dp), PARAMETER :: half = 1.0_dp/2.0_dp, &
          c3 = 1.0_dp/SQRT(3.0_dp), c6 = 1.0_dp/SQRT(6.0_dp), c8 = SQRT(3.0_dp/8.0_dp)
      
      gradphi(1,1) = -half
      gradphi(1,2) = -c3*half
      gradphi(1,3) = -c6*half
      gradphi(2,1) = half
      gradphi(2,2) = -c3*half
      gradphi(2,3) = -c6*half
      gradphi(3,1) = 0
      gradphi(3,2) = c3
      gradphi(3,3) = -c6*half
      gradphi(4,1) = 0
      gradphi(4,2) = 0
      gradphi(4,3) = c8
    END SUBROUTINE dTetraNodalPBasisAll

    SUBROUTINE dTetraNodalLBasisAll(u, v, w, gradphi )
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: gradphi(:,:)
      
      gradphi(1,1:3) = -1.0_dp
      gradphi(2:4,1:3) = 0.0_dp
      gradphi(2,1) = 1.0_dp
      gradphi(3,2) = 1.0_dp
      gradphi(4,3) = 1.0_dp
    END SUBROUTINE dTetraNodalLBasisAll

    
!------------------------------------------------------------------------------    
!>     Tetra edge basis at point (u,v,w)
!------------------------------------------------------------------------------    
    FUNCTION TetraEdgePBasis(edge, i, u, v, w, tetratype) RESULT(value)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of tetras edge function to calculate
!        edge = {1,2,..,12}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: tetratype
!      INPUT: Type of tetrahedron. Defines type of tetrahedron and thus
!        direction for some edge basis functions. Default is 1, tetratype={1,2}
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of tetras edge function i at point (u,v,w), i.e.
!       value = N_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      INTEGER, INTENT(IN), OPTIONAL :: tetratype 
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      
      ! Variables
      INTEGER :: t
      REAL (KIND=dp) :: L1, L2, L3, L4, value
      
      ! Use default type (1) if type not present
      t = 1 
      IF (PRESENT(tetratype)) t = tetratype

      value = 0
      SELECT CASE(edge)
      CASE(1)
         L1=TetraNodalPBasis(1,u,v,w)
         L2=TetraNodalPBasis(2,u,v,w) 
         value = L1*L2*varPhi(i,L2-L1)
      CASE(2)
         L2=TetraNodalPBasis(2,u,v,w)
         L3=TetraNodalPBasis(3,u,v,w) 
         
         ! Choose correct edge function by type
         SELECT CASE(t)
         CASE (1)
            value = L2*L3*varPhi(i,L3-L2)
         CASE (2) 
            value = L2*L3*varPhi(i,L2-L3)
         CASE DEFAULT
            CALL Fatal('PElementBase::TetraEdgePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE(3)
         L1=TetraNodalPBasis(1,u,v,w)
         L3=TetraNodalPBasis(3,u,v,w) 
         value = L1*L3*varPhi(i,L3-L1)
      CASE(4)
         L1=TetraNodalPBasis(1,u,v,w)
         L4=TetraNodalPBasis(4,u,v,w) 
         value = L1*L4*varPhi(i,L4-L1)
      CASE(5)
         L2=TetraNodalPBasis(2,u,v,w)
         L4=TetraNodalPBasis(4,u,v,w) 
         value = L2*L4*varPhi(i,L4-L2)
      CASE(6)
         L3=TetraNodalPBasis(3,u,v,w)
         L4=TetraNodalPBasis(4,u,v,w) 
         value = L3*L4*varPhi(i,L4-L3)
      CASE DEFAULT 
         CALL Fatal('PElementBase::TetraEdgePBasis','Unknown edge for tetrahedron')
      END SELECT
    END FUNCTION TetraEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of tetrahedrons edge basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dTetraEdgePBasis(edge, i, u, v, w, tetratype) RESULT(grad)
!******************************************************************************
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of tetras edge function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: tetratype
!      INPUT: Type of tetrahedron. Defines type of tetrahedron and thus
!        direction for some edge basis functions. Default is 1, tetratype={1,2}
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      INTEGER, INTENT(IN), OPTIONAL :: tetratype 
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      
      ! Variables
      INTEGER :: t
      REAL (KIND=dp) :: Lb, La, vPhi
      REAL (KIND=dp), DIMENSION(3) :: grad, dLb_La, dLb, dLa
      
      ! Use default type (1) if type not present
      t = 1 
      IF (PRESENT(tetratype)) t = tetratype

      grad = 0
      ! Set parameters for gradient 
      SELECT CASE(edge)
      CASE(1)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(2,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(2,u,v,w) 
      CASE(2)
         ! Choose correct edge function by type
         SELECT CASE(t)
         ! Type 1 tetrahedron, edge 4:(2->3)
         CASE (1)
            La = TetraNodalPBasis(2,u,v,w)
            Lb = TetraNodalPBasis(3,u,v,w)
            dLa = dTetraNodalPBasis(2,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w) 
         ! Type 2 tetrahedron, edge 4:(3->2)
         CASE (2) 
            La = TetraNodalPBasis(3,u,v,w)
            Lb = TetraNodalPBasis(2,u,v,w)
            dLa = dTetraNodalPBasis(3,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w) 
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraEdgePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE(3)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(3,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(3,u,v,w) 
      CASE(4)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
      CASE(5)
         La = TetraNodalPBasis(2,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(2,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
      CASE(6)
         La = TetraNodalPBasis(3,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(3,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraEdgePBasis','Unknown edge for tetrahedron')
      END SELECT

      ! Calculate gradient from given parameters
      ! General form for tetra edge gradients is 
      ! 
      ! Grad(Le) = 
      ! Grad(La)*Lb*varPhi(i,Lb-La) + La*Grad(Lb)*varPhi(i,Lb-La) + 
      ! La*Lb*dVarPhi(i,Lb-La)*Grad(Lb-La) 
      dLb_La = dLb-dLa

      vPhi = varPhi(i, Lb-La)
      grad = dLa*Lb*vPhi + La*dLb*vPhi + La*Lb*dVarPhi(i,Lb-La)*dLb_La
    END FUNCTION dTetraEdgePBasis
      

!------------------------------------------------------------------------------
!>     2nd derivatives of tetrahedrons edge basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddTetraEdgePBasis(edge, i, u, v, w, tetratype) RESULT(grad)
!******************************************************************************
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of tetras edge function to calculate
!        edge = {1,2,..,6}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: tetratype
!      INPUT: Type of tetrahedron. Defines type of tetrahedron and thus
!        direction for some edge basis functions. Default is 1, tetratype={1,2}
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      INTEGER, INTENT(IN), OPTIONAL :: tetratype 
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      
      ! Variables
      INTEGER :: t
      REAL (KIND=dp) :: Lb, La, vPhi, grad(3,3)
      REAL (KIND=dp), DIMENSION(3) :: dLb_La, dLb, dLa
      
      ! Use default type (1) if type not present
      t = 1 
      IF (PRESENT(tetratype)) t = tetratype

      ! Set parameters for gradient 
      SELECT CASE(edge)
      CASE(1)
         La  = TetraNodalPBasis(1,u,v,w)
         Lb  = TetraNodalPBasis(2,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(2,u,v,w) 
      CASE(2)
         ! Choose correct edge function by type
         SELECT CASE(t)
         ! Type 1 tetrahedron, edge 4:(2->3)
         CASE (1)
            La  = TetraNodalPBasis(2,u,v,w)
            Lb  = TetraNodalPBasis(3,u,v,w)
            dLa = dTetraNodalPBasis(2,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w) 
         ! Type 2 tetrahedron, edge 4:(3->2)
         CASE (2) 
            La  = TetraNodalPBasis(3,u,v,w)
            Lb  = TetraNodalPBasis(2,u,v,w)
            dLa = dTetraNodalPBasis(3,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w) 
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraEdgePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE(3)
         La  = TetraNodalPBasis(1,u,v,w)
         Lb  = TetraNodalPBasis(3,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(3,u,v,w) 
      CASE(4)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
      CASE(5)
         La = TetraNodalPBasis(2,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(2,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
      CASE(6)
         La = TetraNodalPBasis(3,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(3,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraEdgePBasis','Unknown edge for tetrahedron')
      END SELECT

      ! Calculate gradient from given parameters
      ! General form for tetra edge gradients is 
      ! 
      ! Grad(Le) = 
      ! Grad(La)*Lb*varPhi(i,Lb-La) + La*Grad(Lb)*varPhi(i,Lb-La) + 
      ! La*Lb*dVarPhi(i,Lb-La)*Grad(Lb-La) 

      dLb_La = dLb-dLa
      vPhi = varPhi(i, Lb-La)
!     grad = dLa*Lb*vPhi + La*dLb*vPhi + La*Lb*dVarPhi(i,Lb-La)*dLb_La

      grad = 0
      grad(1,1) = grad(1,1) + dLa(1)*(dLb(1)*vPhi + Lb*dVarPhi(i,Lb-La)*dLb_La(1))
      grad(1,1) = grad(1,1) + dLb(1)*(dLa(1)*vPhi + La*dVarPhi(i,Lb-La)*dLb_La(1))
      grad(1,1) = grad(1,1) + dLa(1)*Lb*dVarPhi(i,Lb-La)*dLb_La(1)
      grad(1,1) = grad(1,1) + La*dLb(1)*dVarPhi(i,Lb-La)*dLb_La(1)
      grad(1,1) = grad(1,1) + La*Lb*ddVarPhi(i,Lb-La)*dLb_La(1)**2

      grad(1,2) = grad(1,2) + dLa(1)*(dLb(2)*vPhi + Lb*dVarPhi(i,Lb-La)*dLb_La(2))
      grad(1,2) = grad(1,2) + dLb(1)*(dLa(2)*vPhi + La*dVarPhi(i,Lb-La)*dLb_La(2))
      grad(1,2) = grad(1,2) + dLa(2)*Lb*dVarPhi(i,Lb-La)*dLb_La(1)
      grad(1,2) = grad(1,2) + La*dLb(2)*dVarPhi(i,Lb-La)*dLb_La(1)
      grad(1,2) = grad(1,2) + La*Lb*ddVarPhi(i,Lb-La)*dLb_La(2)*dLb_La(1)

      grad(1,3) = grad(1,3) + dLa(1)*(dLb(3)*vPhi + Lb*dVarPhi(i,Lb-La)*dLb_La(3))
      grad(1,3) = grad(1,3) + dLb(1)*(dLa(3)*vPhi + La*dVarPhi(i,Lb-La)*dLb_La(3))
      grad(1,3) = grad(1,3) + dLa(3)*Lb*dVarPhi(i,Lb-La)*dLb_La(1)
      grad(1,3) = grad(1,3) + La*dLb(3)*dVarPhi(i,Lb-La)*dLb_La(1)
      grad(1,3) = grad(1,3) + La*Lb*ddVarPhi(i,Lb-La)*dLb_La(3)*dLb_La(1)

      grad(2,2) = grad(2,2) + dLa(2)*(dLb(2)*vPhi + Lb*dVarPhi(i,Lb-La)*dLb_La(2))
      grad(2,2) = grad(2,2) + dLb(2)*(dLa(2)*vPhi + La*dVarPhi(i,Lb-La)*dLb_La(2))
      grad(2,2) = grad(2,2) + dLa(2)*Lb*dVarPhi(i,Lb-La)*dLb_La(2)
      grad(2,2) = grad(2,2) + La*dLb(2)*dVarPhi(i,Lb-La)*dLb_La(2)
      grad(2,2) = grad(2,2) + La*Lb*ddVarPhi(i,Lb-La)*dLb_La(2)**2

      grad(2,3) = grad(2,3) + dLa(2)*(dLb(3)*vPhi + Lb*dVarPhi(i,Lb-La)*dLb_La(3))
      grad(2,3) = grad(2,3) + dLb(2)*(dLa(3)*vPhi + La*dVarPhi(i,Lb-La)*dLb_La(3))
      grad(2,3) = grad(2,3) + dLa(3)*Lb*dVarPhi(i,Lb-La)*dLb_La(2)
      grad(2,3) = grad(2,3) + La*dLb(3)*dVarPhi(i,Lb-La)*dLb_La(2)
      grad(2,3) = grad(2,3) + La*Lb*ddVarPhi(i,Lb-La)*dLb_La(3)*dLb_La(2)

      grad(3,3) = grad(3,3) + dLa(3)*(dLb(3)*vPhi + Lb*dVarPhi(i,Lb-La)*dLb_La(3))
      grad(3,3) = grad(3,3) + dLb(3)*(dLa(3)*vPhi + La*dVarPhi(i,Lb-La)*dLb_La(3))
      grad(3,3) = grad(3,3) + dLa(3)*Lb*dVarPhi(i,Lb-La)*dLb_La(3)
      grad(3,3) = grad(3,3) + La*dLb(3)*dVarPhi(i,Lb-La)*dLb_La(3)
      grad(3,3) = grad(3,3) + La*Lb*ddVarPhi(i,Lb-La)*dLb_La(3)**2

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddTetraEdgePBasis
      

!------------------------------------------------------------------------------    
!>     Tetra face basis at point (u,v,w)
!------------------------------------------------------------------------------    
    FUNCTION TetraFacePBasis(face, i, j, u, v, w, tetratype) RESULT(value)
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of tetras face function to calculate
!        edge = {1,2,3,4}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: tetratype
!      INPUT: Type of tetrahedron. Defines type of tetrahedron and thus
!        direction for some face basis functions. Default is 1, tetratype={1,2}
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of tetras face function m(i,j) at point (u,v,w), i.e.
!       value = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      INTEGER, INTENT(IN), OPTIONAL :: tetratype
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      INTEGER :: t
      REAL (KIND=dp) :: L1, L2, L3, L4, value
      
      ! Use default type (1) if type not present
      t = 1 
      IF (PRESENT(tetratype)) t = tetratype

      value = 0
      SELECT CASE(face)
      CASE (1)
         L1 = TetraNodalPBasis(1,u,v,w)
         L2 = TetraNodalPBasis(2,u,v,w)
         L3 = TetraNodalPBasis(3,u,v,w)
            
         SELECT CASE(t)   
         CASE (1)
            value = L1*L2*L3*LegendreP(i,L2-L1)*LegendreP(j,2*L3-1)
         CASE (2)
            value = L1*L2*L3*LegendreP(i,L3-L1)*LegendreP(j,2*L2-1)
         CASE DEFAULT
            CALL Fatal('PElementBase::TetraFacePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE (2)
         L1 = TetraNodalPBasis(1,u,v,w)
         L2 = TetraNodalPBasis(2,u,v,w)
         L4 = TetraNodalPBasis(4,u,v,w)
         value = L1*L2*L4*LegendreP(i,L2-L1)*LegendreP(j,2*L4-1)
      CASE (3)
         L2=TetraNodalPBasis(2,u,v,w)
         L3=TetraNodalPBasis(3,u,v,w)
         L4=TetraNodalPBasis(4,u,v,w)

         SELECT CASE(t)
         CASE (1)
            value = L2*L3*L4*LegendreP(i,L3-L2)*LegendreP(j,2*L4-1)
         CASE (2) 
            value = L2*L3*L4*LegendreP(i,L2-L3)*LegendreP(j,2*L4-1)
         CASE DEFAULT
            CALL Fatal('PElementBase::TetraFacePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE (4)
         L1=TetraNodalPBasis(1,u,v,w)
         L3=TetraNodalPBasis(3,u,v,w)
         L4=TetraNodalPBasis(4,u,v,w)
         value = L1*L3*L4*LegendreP(i,L3-L1)*LegendreP(j,2*L4-1)
      CASE DEFAULT 
         CALL Fatal('PElementBase::TetraFacePBasis','Unknown face for tetrahedron')
      END SELECT

    END FUNCTION TetraFacePBasis


!------------------------------------------------------------------------------
!>     Gradient of tetra face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dTetraFacePBasis(face, i, j, u, v, w, tetratype) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of tetras face function to calculate
!        edge = {1,2,..,4}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: tetratype
!      INPUT: Type of tetrahedron. Defines type of tetrahedron and thus
!        direction for some face basis functions. Default is 1, tetratype={1,2}
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      INTEGER, INTENT(IN), OPTIONAL :: tetratype
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      INTEGER :: t
      REAL (KIND=dp) :: La, Lb, Lc, Legi, Legj
      REAL (KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLb_La, dLc_1, grad
      
      ! Use default type (1) if type not present
      t = 1 
      IF (PRESENT(tetratype)) t = tetratype

      SELECT CASE(face)
      CASE (1)
         SELECT CASE (t)
         CASE (1)
            La = TetraNodalPBasis(1,u,v,w)
            Lb = TetraNodalPBasis(2,u,v,w)
            Lc = TetraNodalPBasis(3,u,v,w)

            dLa = dTetraNodalPbasis(1,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w)
            dLc = dTetraNodalPBasis(3,u,v,w)
         CASE (2)
            La = TetraNodalPBasis(1,u,v,w)
            Lb = TetraNodalPBasis(3,u,v,w)
            Lc = TetraNodalPBasis(2,u,v,w)

            dLa = dTetraNodalPbasis(1,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w)
            dLc = dTetraNodalPBasis(2,u,v,w)
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraFacePBasis','Unknown type for tetrahedron')       
         END SELECT
      CASE (2)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(2,u,v,w)
         Lc = TetraNodalPBasis(4,u,v,w)

         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(2,u,v,w)
         dLc = dTetraNodalPBasis(4,u,v,w)
      CASE (3)
         SELECT CASE(t)
         ! Type 1 tetrahedron: Face 4:(2,3,4)
         CASE (1)
            La = TetraNodalPBasis(2,u,v,w)
            Lb = TetraNodalPBasis(3,u,v,w)
            Lc = TetraNodalPBasis(4,u,v,w)

            dLa = dTetraNodalPBasis(2,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w)
            dLc = dTetraNodalPBasis(4,u,v,w)
         ! Type 2 tetrahedron: Face 4:(3,2,4)
         CASE (2)
            La = TetraNodalPBasis(3,u,v,w)
            Lb = TetraNodalPBasis(2,u,v,w)
            Lc = TetraNodalPBasis(4,u,v,w)

            dLa = dTetraNodalPBasis(3,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w)
            dLc = dTetraNodalPBasis(4,u,v,w)
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraFacePBasis','Unknown type for tetrahedron')
         END SELECT
         ! Derivative of second inner function is equal for both types
      CASE (4)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(3,u,v,w)
         Lc = TetraNodalPBasis(4,u,v,w)

         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(3,u,v,w)
         dLc = dTetraNodalPBasis(4,u,v,w)
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraFacePBasis','Unknown face for tetrahedron')
      END SELECT

      Legi = LegendreP(i, Lb-La)
      Legj = LegendreP(j, 2*Lc-1)

      ! Calculate gradient from given parameters 
      grad = dLa*Lb*Lc*Legi*Legj + La*dLb*Lc*Legi*Legj + La*Lb*dLc*Legi*Legj + &
         La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb-dLa)*Legj + La*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc
    
    END FUNCTION dTetraFacePBasis

!------------------------------------------------------------------------------
!>     2nd derivatives of tetra face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddTetraFacePBasis(face, i, j, u, v, w, tetratype) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of tetras face function to calculate
!        edge = {1,2,..,4}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: tetratype
!      INPUT: Type of tetrahedron. Defines type of tetrahedron and thus
!        direction for some face basis functions. Default is 1, tetratype={1,2}
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: face, i, j
      INTEGER, INTENT(IN), OPTIONAL :: tetratype
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      INTEGER :: t
      REAL (KIND=dp) :: La, Lb, Lc, Legi, Legj, grad(3,3)
      REAL (KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLb_La, dLc_1
      
      ! Use default type (1) if type not present
      t = 1 
      IF (PRESENT(tetratype)) t = tetratype

      SELECT CASE(face)
      CASE (1)
         SELECT CASE (t)
         CASE (1)
            La = TetraNodalPBasis(1,u,v,w)
            Lb = TetraNodalPBasis(2,u,v,w)
            Lc = TetraNodalPBasis(3,u,v,w)
            dLa = dTetraNodalPbasis(1,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w)
            dLc = dTetraNodalPBasis(3,u,v,w)
         CASE (2)
            La = TetraNodalPBasis(1,u,v,w)
            Lb = TetraNodalPBasis(3,u,v,w)
            Lc = TetraNodalPBasis(2,u,v,w)
            dLa = dTetraNodalPbasis(1,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w)
            dLc = dTetraNodalPBasis(2,u,v,w)
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraFacePBasis','Unknown type for tetrahedron')       
         END SELECT
      CASE (2)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(2,u,v,w)
         Lc = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(2,u,v,w)
         dLc = dTetraNodalPBasis(4,u,v,w)
      CASE (3)
         SELECT CASE(t)
         ! Type 1 tetrahedron: Face 4:(2,3,4)
         CASE (1)
            La = TetraNodalPBasis(2,u,v,w)
            Lb = TetraNodalPBasis(3,u,v,w)
            Lc = TetraNodalPBasis(4,u,v,w)
            dLa = dTetraNodalPBasis(2,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w)
            dLc = dTetraNodalPBasis(4,u,v,w)
         ! Type 2 tetrahedron: Face 4:(3,2,4)
         CASE (2)
            La = TetraNodalPBasis(3,u,v,w)
            Lb = TetraNodalPBasis(2,u,v,w)
            Lc = TetraNodalPBasis(4,u,v,w)
            dLa = dTetraNodalPBasis(3,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w)
            dLc = dTetraNodalPBasis(4,u,v,w)
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraFacePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE (4)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(3,u,v,w)
         Lc = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(3,u,v,w)
         dLc = dTetraNodalPBasis(4,u,v,w)
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraFacePBasis','Unknown face for tetrahedron')
      END SELECT

      
      Legi = LegendreP(i, Lb-La)
      Legj = LegendreP(j, 2*Lc-1)

      ! Calculate gradient from given parameters 
!     dLb_La = dLb-dLa
!     dLc_1 = 2*dLc

!     grad = dLa*Lb*Lc*Legi*Legj + La*dLb*Lc*Legi*Legj + La*Lb*dLc*Legi*Legj + &
!          La*Lb*Lc*dLegendreP(i,Lb-La)*dLb_La*Legj + La*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*dLc_1

      grad = 0
      grad(1,1) = grad(1,1) + dLa(1)*(dLb(1)*Lc*Legi*Legj + Lb*dLc(1)*Legi*Legj + &
            Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                              Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) )
      grad(1,1) = grad(1,1) + dLb(1)*(dLa(1)*Lc*Legi*Legj + La*dLc(1)*Legi*Legj + &
            La*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                              La*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) )
      grad(1,1) = grad(1,1) + dLc(1)*(dLa(1)*Lb*Legi*Legj + La*dLb(1)*Legi*Legj + &
            La*Lb*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                              La*Lb*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) )
      grad(1,1) = grad(1,1) + dLa(1)*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
            La*dLb(1)*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
               La*Lb*dLc(1)*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                  La*Lb*Lc*ddLegendreP(i,Lb-La)*(dLb(1)-dLa(1))**2*Legj + &
                       La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*dLegendreP(j,2*Lc-1)*2*dLc(1)
      grad(1,1) = grad(1,1) + dLa(1)*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
            La*dLb(1)*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                La*Lb*dLc(1)*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                    La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                         La*Lb*Lc*Legi*ddLegendreP(j,2*Lc-1)*4*dLc(1)**2

      grad(2,2) = grad(2,2) + dLa(2)*(dLb(2)*Lc*Legi*Legj + Lb*dLc(2)*Legi*Legj + &
            Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) )
      grad(2,2) = grad(2,2) + dLb(2)*(dLa(2)*Lc*Legi*Legj + La*dLc(2)*Legi*Legj + &
            La*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                La*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) )
      grad(2,2) = grad(2,2) + dLc(2)*(dLa(2)*Lb*Legi*Legj + La*dLb(2)*Legi*Legj + &
            La*Lb*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                La*Lb*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) )
      grad(2,2) = grad(2,2) + dLa(2)*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
            La*dLb(2)*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                La*Lb*dLc(2)*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                    La*Lb*Lc*ddLegendreP(i,Lb-La)*(dLb(2)-dLa(2))**2*Legj + &
                         La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*dLegendreP(j,2*Lc-1)*2*dLc(2)
      grad(2,2) = grad(2,2) + dLa(2)*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
            La*dLb(2)*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
                 La*Lb*dLc(2)*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
                     La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
                          La*Lb*Lc*Legi*ddLegendreP(j,2*Lc-1)*4*dLc(2)**2
    
      grad(3,3) = grad(3,3) + dLa(3)*(dLb(3)*Lc*Legi*Legj + Lb*dLc(3)*Legi*Legj + &
            Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
               Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(3,3) = grad(3,3) + dLb(3)*(dLa(3)*Lc*Legi*Legj + La*dLc(3)*Legi*Legj + &
            La*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
               La*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(3,3) = grad(3,3) + dLc(3)*(dLa(3)*Lb*Legi*Legj + La*dLb(3)*Legi*Legj + &
            La*Lb*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
               La*Lb*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(3,3) = grad(3,3) + dLa(3)*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
            La*dLb(3)*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
                La*Lb*dLc(3)*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
                     La*Lb*Lc*ddLegendreP(i,Lb-La)*(dLb(3)-dLa(3))**2*Legj + &
                           La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*dLegendreP(j,2*Lc-1)*2*dLc(3)
      grad(3,3) = grad(3,3) + dLa(3)*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) + &
            La*dLb(3)*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) + &
                 La*Lb*dLc(3)*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) + &
                     La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*dLegendreP(j,2*Lc-1)*2*dLc(3) + &
                          La*Lb*Lc*Legi*ddLegendreP(j,2*Lc-1)*4*dLc(3)**2


      grad(1,2) = grad(1,2) + dLa(1)*(dLb(2)*Lc*Legi*Legj + Lb*dLc(2)*Legi*Legj + &
            Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
               Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) )
      grad(1,2) = grad(1,2) + dLb(1)*(dLa(2)*Lc*Legi*Legj + La*dLc(2)*Legi*Legj + &
            La*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
               La*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) )
      grad(1,2) = grad(1,2) + dLc(1)*(dLa(2)*Lb*Legi*Legj + La*dLb(2)*Legi*Legj + &
            La*Lb*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
               La*Lb*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) )
      grad(1,2) = grad(1,2) + dLa(2)*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
           La*dLb(2)*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                La*Lb*dLc(2)*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                    La*Lb*Lc*ddLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*(dLb(2)-dLa(2))*Legj + &
                          La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*dLegendreP(j,2*Lc-1)*2*dLc(2)
      grad(1,2) = grad(1,2) + dLa(2)*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
            La*dLb(2)*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                  La*Lb*dLc(2)*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                       La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                              La*Lb*Lc*Legi*ddLegendreP(j,2*Lc-1)*4*dLc(1)*dLc(2)
    
      grad(1,3) = grad(1,3) + dLa(1)*(dLb(3)*Lc*Legi*Legj + Lb*dLc(3)*Legi*Legj + &
            Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
                Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(1,3) = grad(1,3) + dLb(1)*(dLa(3)*Lc*Legi*Legj + La*dLc(3)*Legi*Legj + &
            La*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
                La*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(1,3) = grad(1,3) + dLc(1)*(dLa(3)*Lb*Legi*Legj + La*dLb(3)*Legi*Legj + &
            La*Lb*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
                La*Lb*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(1,3) = grad(1,3) + dLa(3)*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
            La*dLb(3)*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                 La*Lb*dLc(3)*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*Legj + &
                     La*Lb*Lc*ddLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*(dLb(3)-dLa(3))*Legj + &
                          La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(1)-dLa(1))*dLegendreP(j,2*Lc-1)*2*dLc(3)
      grad(1,3) = grad(1,3) + dLa(3)*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
            La*dLb(3)*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                 La*Lb*dLc(3)*Legi*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                     La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*dLegendreP(j,2*Lc-1)*2*dLc(1) + &
                           La*Lb*Lc*Legi*ddLegendreP(j,2*Lc-1)*4*dLc(1)*dLc(3)

      grad(2,3) = grad(2,3) + dLa(2)*(dLb(3)*Lc*Legi*Legj + Lb*dLc(3)*Legi*Legj + &
            Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
               Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(2,3) = grad(2,3) + dLb(2)*(dLa(3)*Lc*Legi*Legj + La*dLc(3)*Legi*Legj + &
            La*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
               La*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(2,3) = grad(2,3) + dLc(2)*(dLa(3)*Lb*Legi*Legj + La*dLb(3)*Legi*Legj + &
            La*Lb*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*Legj + &
               La*Lb*Legi*dLegendreP(j,2*Lc-1)*2*dLc(3) )
      grad(2,3) = grad(2,3) + dLa(3)*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
            La*dLb(3)*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                 La*Lb*dLc(3)*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*Legj + &
                     La*Lb*Lc*ddLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*(dLb(3)-dLa(3))*Legj + &
                          La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(2)-dLa(2))*dLegendreP(j,2*Lc-1)*2*dLc(3)
      grad(2,3) = grad(2,3) + dLa(3)*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
            La*dLb(3)*Lc*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
                  La*Lb*dLc(3)*Legi*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
                       La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb(3)-dLa(3))*dLegendreP(j,2*Lc-1)*2*dLc(2) + &
                            La*Lb*Lc*Legi*ddLegendreP(j,2*Lc-1)*4*dLc(2)*dLc(3)
    
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddTetraFacePBasis


!------------------------------------------------------------------------------
!>    Tetra bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION TetraBubblePBasis(i, j, k, u, v, w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,0),(0,0,1),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of tetras bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL (KIND=dp) :: L1, L2, L3, L4, value

      L1 = TetraNodalPBasis(1,u,v,w)
      L2 = TetraNodalPBasis(2,u,v,w)
      L3 = TetraNodalPBasis(3,u,v,w)
      L4 = TetraNodalPBasis(4,u,v,w)
      value = L1*L2*L3*L4*LegendreP(i,L2-L1)*LegendreP(j,2*L3-1)*LegendreP(k,2*L4-1)
    END FUNCTION TetraBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of tetra bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dTetraBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,0),(0,0,1),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL (KIND=dp) :: L1, L2, L3, L4, a, b, c
      REAL (KIND=dp), DIMENSION(3) :: grad, grad1, dL1, dL2, dL3, dL4, da, db, dc

      grad = 0
      L1 = TetraNodalPBasis(1,u,v,w)
      L2 = TetraNodalPBasis(2,u,v,w)
      L3 = TetraNodalPBasis(3,u,v,w)
      L4 = TetraNodalPBasis(4,u,v,w)

      dL1 = dTetraNodalPBasis(1,u,v,w)
      dL2 = dTetraNodalPBasis(2,u,v,w)
      dL3 = dTetraNodalPBasis(3,u,v,w)
      dL4 = dTetraNodalPBasis(4,u,v,w)

      a = LegendreP(i,L2-L1)
      b = LegendreP(j,2*L3-1)
      c = LegendreP(k,2*L4-1)

      da = dLegendreP(i,L2-L1)*(dL2-dL1)
      db = dLegendreP(j,2*L3-1)*2*dL3
      dc = dLegendreP(k,2*L4-1)*2*dL4

      ! Gradients of tetrahedral bubble basis functions 
      grad = (dL1*L2*L3*L4 + L1*dL2*L3*L4 + L1*L2*dL3*L4 + L1*L2*L3*dL4)*a*b*c + &
                   L1*L2*L3*L4*(da*b*c + a*db*c + a*b*dc)
    END FUNCTION dTetraBubblePBasis


!------------------------------------------------------------------------------
!>    2nd derivatives of tetra bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddTetraBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,0),(0,0,1),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of tetras bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i, j, k
      REAL (KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL (KIND=dp) :: L1, L2, L3, L4, L2_L1, L3_1, L4_1, a, b, c 
      REAL (KIND=dp), DIMENSION(3,3) :: grad

      REAL(KIND=dp) :: dL1(3),dL2(3),dL3(3),dL4(3), r,s,t
      REAL(KIND=dp) :: L12,dL12(3),ddL12(3,3)
      REAL(KIND=dp) :: L34,dL34(3),ddL34(3,3)
      REAL(KIND=dp) :: L1234,dL1234(3),ddL1234(3,3)
      REAL(KIND=dp) :: G1,dG1(3), ddG1(3,3)
      REAL(KIND=dp) :: G2,dG2(3), ddG2(3,3)
      REAL(KIND=dp) :: G3,dG3(3), ddG3(3,3)
      REAL(KIND=dp) :: G12,dG12(3), ddG12(3,3)
      REAL(KIND=dp) :: G123, dG123(3), ddG123(3,3)
      REAL(KIND=dp) :: L1234G123, dL1234G123(3), ddL1234G123(3,3)
      INTEGER :: p,q

!     value = L1*L2*L3*L4*LegendreP(i,L2-L1)*LegendreP(j,2*L3-1)*LegendreP(k,2*L4-1)

      L1 = TetraNodalPBasis(1,u,v,w)
      L2 = TetraNodalPBasis(2,u,v,w)
      L3 = TetraNodalPBasis(3,u,v,w)
      L4 = TetraNodalPBasis(4,u,v,w)

      dL1 = dTetraNodalPBasis(1,u,v,w)
      dL2 = dTetraNodalPBasis(2,u,v,w)
      dL3 = dTetraNodalPBasis(3,u,v,w)
      dL4 = dTetraNodalPBasis(4,u,v,w)

      L12 = L1*L2
      dL12 = dL1*L2 + L1*dL2
      DO p=1,3
        DO q=p,3
          ddL12(p,q) = dL1(p)*dL2(q) + dL1(q)*dL2(p)
        END DO        
      END DO        
 
      L34 = L3*L4
      dL34 = dL3*L4 + L3*dL4
      DO p=1,3
        DO q=p,3
          ddL34(p,q) = dL3(p)*dL4(q) + dL3(q)*dL4(p)
        END DO        
      END DO        

      L1234 = L12*L34
      dL1234 = dL12*L34 + L12*dL34
      DO p=1,3
        DO q=p,3
          ddL1234(p,q) = ddL12(p,q)*L34 + dL12(p)*dL34(q) + dL12(q)*dL34(p) + L12*ddL34(p,q)
        END DO        
      END DO        

      G1 = LegendreP(i,L2-L1)
      G2 = LegendreP(j,2*L3-1)
      G3 = LegendreP(k,2*L4-1)

      dG1 = dLegendreP(i,L2-L1)*(dL2-dL1)
      dG2 = dLegendreP(j,2*L3-1)*2*dL3
      dG3 = dLegendreP(k,2*L4-1)*2*dL4

      r = ddLegendreP(i,L2-L1)
      s = ddLegendreP(j,2*L3-1)
      t = ddLegendreP(k,2*L4-1)
      DO p=1,3
        DO q=p,3
          ddG1(p,q) = r * (dL2(p)-dL1(p))*(dL2(q)-dL1(q))
          ddG2(p,q) = s * 4*dL3(p)*dL3(q)
          ddG3(p,q) = t * 4*dL4(p)*dL4(q)
        END DO
      END DO

      G12 = G1*G2
      dG12 = dG1*G2 + G1*dG2
      DO p=1,3
        DO q=p,3
          ddG12(p,q) = ddG1(p,q)*G2 + dG1(p)*dG2(q) +dG1(q)*dG2(p) + G1*ddG2(p,q)
        END DO
      END DO

      G123 = G12 * G3
      dG123 = dG12*G3 + G12*dG3
      DO p=1,3
        DO q=p,3
         ddG123(p,q) = ddG12(p,q)*G3 + dG12(p)*dG3(q) + dG12(q)*dG3(p) + G12*ddG3(p,q)
        END DO
      END DO

      ! dd(L1*L2*L3*L4*G1*G2*G3)
      DO p=1,3
        DO q=p,3
          grad(p,q)=ddL1234(p,q)*G123+dL1234(p)*dG123(q)+dL1234(q)*dG123(p)+L1234*ddG123(p,q)
        END DO
      END DO

      grad(2,1)=grad(1,2)
      grad(3,1)=grad(1,3)
      grad(3,2)=grad(2,3)
    END FUNCTION ddTetraBubblePBasis


!------------------------------------------------------------------------------
!>     Wedge nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION WedgeNodalPBasis(node, u, v, w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of wedge nodal function to calculate
!        node = {1,2,3,4,5,6}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges nodal function at point (u,v,w), i.e.
!       value = N_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Result
      REAL(KIND=dp) :: value 

      value = 0
      SELECT CASE (node)
      CASE (1)
         value = WedgeL(1,u,v)*(1-w)
      CASE (2)
         value = WedgeL(2,u,v)*(1-w)
      CASE (3)
         value = WedgeL(3,u,v)*(1-w)
      CASE (4)
         value = WedgeL(1,u,v)*(1+w)
      CASE (5)   
         value = WedgeL(2,u,v)*(1+w)
      CASE (6)
         value = WedgeL(3,u,v)*(1+w)
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeNodalPBasis','Unknown node for wedge')
      END SELECT
      value = value/2
    END FUNCTION WedgeNodalPBasis


    SUBROUTINE WedgeNodalPBasisAll(u, v, w, phi) 
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: phi(:)
      REAL(KIND=dp) :: tri(3),line(2)
      REAL(KIND=dp), PARAMETER :: half = 1.0_dp/2.0_dp, c3 = 1.0_dp/SQRT(3.0_dp)
                 
      tri(1) = half*(1d0-u-c3*v)
      tri(2) = half*(1d0+u-c3*v)
      tri(3) = c3*v

      line(1) = half*(1-w)
      line(2) = half*(1+w)
      
      phi(1:3) = line(1)*tri
      phi(4:6) = line(2)*tri
    END SUBROUTINE WedgeNodalPBasisAll

    SUBROUTINE WedgeNodalLBasisAll(u, v, w, phi) 
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: phi(:)
      REAL(KIND=dp) :: tri(3),line(2)
      REAL(KIND=dp), PARAMETER :: half = 1.0_dp/2.0_dp
                 
      tri(1) = 1.0_dp-u-v
      tri(2) = u
      tri(3) = v

      line(1) = half*(1-w)
      line(2) = half*(1+w)
      
      phi(1:3) = line(1)*tri
      phi(4:6) = line(2)*tri
    END SUBROUTINE WedgeNodalLBasisAll

    
    
!------------------------------------------------------------------------------
!>     Gradient of wedges nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dWedgeNodalPBasis(node, u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of wedge nodal function to calculate
!        node = {1,2,..,6}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges nodal function at point (u,v,w), i.e.
!       grad = dN_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: signW, dL(3), L, grad(3) 

      grad = 0
      SELECT CASE(node)
      CASE (1,2,3)
         signW = -1
      CASE (4,5,6)
         signW = 1
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeNodalPBasis','Unknown node for wedge')
      END SELECT

      ! Calculate gradient from the general form
      dL(1:3) = dWedgeL(node,u,v)
      L = WedgeL(node,u,v)
      grad(1) = 1d0/2*dL(1)*(1+signW*w)  
      grad(2) = 1d0/2*dL(2)*(1+signW*w)
      grad(3) = signW*1d0/2*L
    END FUNCTION dWedgeNodalPBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of wedges nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddWedgeNodalPBasis(node, u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of wedge nodal function to calculate
!        node = {1,2,..,6}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3,3)
!       gradient of wedges nodal function at point (u,v,w), i.e.
!       grad = dN_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: signW, dL(3), L, grad(3,3)

      grad = 0
      SELECT CASE(node)
      CASE (1,2,3)
         signW = -1
      CASE (4,5,6)
         signW = 1
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeNodalPBasis','Unknown node for wedge')
      END SELECT

      ! Calculate gradient from the general form
      dL = dWedgeL(node,u,v)
!     grad(1) = 1d0/2*dL(1)*(1+signW*W)
!     grad(2) = 1d0/2*dL(2)*(1+signW*W)
!     grad(3) = signW*1d0/2*L

      grad = 0
      grad(1,3) = dL(1)*signW/2
      grad(2,3) = dL(2)*signW/2
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddWedgeNodalPBasis


    SUBROUTINE dWedgeNodalPBasisAll(u, v, w, gradphi) 
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: gradphi(:,:)
      REAL(KIND=dp) :: tri(3),line(2),gradtri(3,2),gradline(2)

      REAL(KIND=dp), PARAMETER :: half=1.0_dp/2, c3=1.0_dp/SQRT(3.0_dp)

      tri(1) = half*(1d0-u-c3*v)
      tri(2) = half*(1d0+u-c3*v)
      tri(3) = c3*v

      line(1) = half*(1-w)
      line(2) = half*(1+w)
      
      gradtri(1,1) = -half
      gradtri(1,2) = -half*c3
      gradtri(2,1) = half
      gradtri(2,2) = -half*c3
      gradtri(3,1) = 0
      gradtri(3,2) = c3

      gradline(1) = -half
      gradline(2) = half
          
      gradphi(1:3,1:2) = gradtri * line(1)
      gradphi(4:6,1:2) = gradtri * line(2)
      
      gradphi(1:3,3) = tri * gradline(1)
      gradphi(4:6,3) = tri * gradline(2)
      
    END SUBROUTINE dWedgeNodalPBasisAll

    SUBROUTINE dWedgeNodalLBasisAll(u, v, w, gradphi) 
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      REAL(KIND=dp), INTENT(OUT) :: gradphi(:,:)
      REAL(KIND=dp) :: tri(3),line(2),gradtri(3,2),gradline(2)
      REAL(KIND=dp), PARAMETER :: half=1.0_dp/2.0_dp

      tri(1) = 1.0_dp-u-v
      tri(2) = u
      tri(3) = v

      line(1) = half*(1-w)
      line(2) = half*(1+w)

      gradtri(1,1) = -1.0_dp
      gradtri(1,2) = -1.0_dp
      gradtri(2,1) = 1.0_dp
      gradtri(2,2) = 0.0_dp
      gradtri(3,1) = 0.0_dp
      gradtri(3,2) = 1.0_dp

      gradline(1) = -half
      gradline(2) = half
          
      gradphi(1:3,1:2) = gradtri * line(1)
      gradphi(4:6,1:2) = gradtri * line(2)
      
      gradphi(1:3,3) = tri * gradline(1)
      gradphi(4:6,3) = tri * gradline(2)
      
    END SUBROUTINE dWedgeNodalLBasisAll


! --- start serendipity wedge

!------------------------------------------------------------------------------      
!>     Wedge edge basis at point (u,v,w)
!------------------------------------------------------------------------------      
    FUNCTION SD_WedgeEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(value)
!------------------------------------------------------------------------------      
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of wedges edge function to calculate
!        edge = {1,2,..,9}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges edge function i at point (u,v,w), i.e.
!       value = N_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp) :: parW, La, Lb, tmp, value
      LOGICAL :: invert
      
      ! Edge is not inverted by default 
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      value = 0
      ! Set sign of w for edges 1,..,6
      SELECT CASE(edge)
      CASE (1,2,3)
         parW = -w
      CASE (4,5,6)
         parW = w
      END SELECT

      SELECT CASE(edge)
      CASE (1,4)
         La = WedgeL(1,u,v)
         Lb = WedgeL(2,u,v)
      CASE (2,5)
         La = WedgeL(2,u,v)
         Lb = WedgeL(3,u,v)
      CASE (3,6)
         La = WedgeL(3,u,v)
         Lb = WedgeL(1,u,v)
      CASE (7,8,9)
         ! Invert edge if needed
         IF (invert) THEN
            parW = -w
         ELSE 
            parW = w
         END IF

         value = WedgeL(edge-6,u,v)*Phi(i,parW)
         RETURN
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeEdgePBasis','Unknown edge for wedge')
      END SELECT

      ! Swap parameters for inverted edges
      IF (invert) THEN
         tmp = La
         La = Lb
         Lb = tmp
      END IF

      ! Calculate value from general form for edges 1-6
      value = 1d0/2*La*Lb*varPhi(i,Lb-La)*(1+parW)
    END FUNCTION SD_WedgeEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of wedge edge basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_dWedgeEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of wedges edge function to calculate
!        edge = {1,2,..,9}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp) :: parW, La, Lb, phiI, tmp, &
           grad(3), dLa(3), dLb(3), dW(3), dtmp(3) 
      LOGICAL :: invert 
      
      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      grad = 0
      ! Set sign of w and derivative of w
      dW = 0
      SELECT CASE(edge)
      CASE (1,2,3)
         parW = -w
         dW(3) = -1
      CASE (4,5,6)
         parW = w
         dW(3) = 1
      END SELECT

      SELECT CASE(edge)
      CASE (1,4)
         La = WedgeL(1,u,v)
         Lb = WedgeL(2,u,v)
         dLa(1:3) = dWedgeL(1,u,v)
         dLb(1:3) = dWedgeL(2,u,v)
      CASE (2,5)
         La = WedgeL(2,u,v)
         Lb = WedgeL(3,u,v)
         dLa(1:3) = dWedgeL(2,u,v)
         dLb(1:3) = dWedgeL(3,u,v)
      CASE (3,6)
         La = WedgeL(3,u,v)
         Lb = WedgeL(1,u,v)
         dLa(1:3) = dWedgeL(3,u,v)
         dLb(1:3) = dWedgeL(1,u,v)
      CASE (7,8,9)
         ! Invert edge if needed
         IF (invert) THEN
            parW = -w
            dW(3) = -1
         ELSE 
            parW = w
            dW(3) = 1
         END IF

         phiI = Phi(i,parW)
         dLa(1:3) = dWedgeL(edge-6,u,v)
         
         ! Calculate value of edge function and return
         grad(1) = dLa(1)*phiI
         grad(2) = dLa(2)*phiI
         grad(3) = WedgeL(edge-6,u,v)*dPhi(i,parW)*dW(3)
         RETURN
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeEdgePBasis','Unknown edge for wedge')
      END SELECT

      ! Swap parameters for inverted edges
      IF (invert) THEN
         tmp = La
         La = Lb
         Lb = tmp
         dtmp(1:3) = dLa
         dLa(1:3) = dLb
         dLb(1:3) = dtmp
      END IF      
      
      phiI = varPhi(i,Lb-La)

      ! Calculate value of function from general form
      grad = (dLa*Lb*phiI*(1+parW) + La*dLb*phiI*(1+parW) + &
           La*Lb*dVarPhi(i,Lb-La)*(dLb-dLa)*(1+parW) + La*Lb*phiI*dW)/2
    END FUNCTION SD_dWedgeEdgePBasis

!------------------------------------------------------------------------------
!>     2nd derivatives of wedge edge basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_ddWedgeEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of wedges edge function to calculate
!        edge = {1,2,..,9}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp) :: parW, La, Lb, phiI, tmp, &
           grad(3,3), dLa(3), dLb(3), dW(3), dtmp(3) 
      LOGICAL :: invert 
      
      INTEGER :: p,q
      REAL(KIND=dp) :: f(4), df(4,3), ddf(4,3,3), ddPhiI

      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      grad = 0
      ! Set sign of w and derivative of w
      dW = 0
      SELECT CASE(edge)
      CASE (1,2,3)
         parW = -w
         dW(3) = -1
      CASE (4,5,6)
         parW = w
         dW(3) = 1
      END SELECT

      SELECT CASE(edge)
      CASE (1,4)
         La = WedgeL(1,u,v)
         Lb = WedgeL(2,u,v)
         dLa(1:3) = dWedgeL(1,u,v)
         dLb(1:3) = dWedgeL(2,u,v)
      CASE (2,5)
         La = WedgeL(2,u,v)
         Lb = WedgeL(3,u,v)
         dLa(1:3) = dWedgeL(2,u,v)
         dLb(1:3) = dWedgeL(3,u,v)
      CASE (3,6)
         La = WedgeL(3,u,v)
         Lb = WedgeL(1,u,v)
         dLa(1:3) = dWedgeL(3,u,v)
         dLb(1:3) = dWedgeL(1,u,v)
      CASE (7,8,9)
         ! Invert edge if needed
         IF (invert) THEN
            parW = -w
            dW(3) = -1
         ELSE 
            parW = w
            dW(3) = 1
         END IF

         phiI = Phi(i,parW)
         dLa(1:3) = dWedgeL(edge-6,u,v)
         
         ! Calculate value of edge function and return
!        grad(1) = dLa(1)*phiI
!        grad(2) = dLa(2)*phiI
!        grad(3) = WedgeL(edge-6,u,v)*dPhi(i,parW)*dW(3)

         grad=0
         grad(1,3) = dLa(1)*dPhi(i,parW)*dW(3)
         grad(3,1) = grad(1,3)
         grad(2,3) = dLa(2)*dPhi(i,parW)*dW(3)
         grad(3,2) = grad(2,3)
         grad(3,3) = WedgeL(edge-6,u,v)*ddPhi(i,parW)

         RETURN
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeEdgePBasis','Unknown edge for wedge')
      END SELECT

      ! Swap parameters for inverted edges
      IF (invert) THEN
         tmp = La
         La = Lb
         Lb = tmp
         dtmp(1:3) = dLa
         dLa(1:3) = dLb
         dLb(1:3) = dtmp
      END IF      
      
      phiI = varPhi(i,Lb-La)

      ! Calculate value of function from general form
!     grad = (dLa*Lb*phiI*(1+parW) + La*dLb*phiI*(1+parW) + &
!              La*Lb*dVarPhi(i,Lb-La)*(dLb-dLa)*(1+parW) + La*Lb*phiI*dW)/2


      f(1)=La; f(2)=Lb; f(3)=PhiI; f(4)=1+parW
      df(1,:) = dLa
      df(2,:) = dLb
      df(3,:) = dVarPhi(i,Lb-La)*(dLb-dLa)
      df(4,:) = dW

      ddPhiI = ddVarPhi(i,Lb-La)
      ddf = 0
      DO p=1,3
        DO q=p,3
          ddf(3,p,q) = ddPhiI*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
        END DO
      END DO

      grad = Product2ndDerivatives(4,f,df,ddf,3,0)/2

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)

    END FUNCTION SD_ddWedgeEdgePBasis


!------------------------------------------------------------------------------
!>     Wedge face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_WedgeFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of wedges face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges face function m(i,j) at point (u,v,w), i.e.
!       value = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      REAL(KIND=dp) :: La, Lb, Lc, Lha, Lhc, parW, value
      INTEGER :: local(4)

      ! If local numbering not present, use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local = 0
         local(1:4) = getWedgeFaceMap(face)
      ! Numbering was present. Use it
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF

      ! Set sign of w for faces 1 and 2
      SELECT CASE (face)
      CASE (1)
         parW = -w
      CASE (2)
         parW = w
      END SELECT
      
      ! Get value of face function
      value = 0
      SELECT CASE(face)
         CASE (1,2)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            Lc = WedgeL(local(3),u,v)
            value = 1d0/2*(1+parW)*LegendreP(i,Lb-La)*LegendreP(j,2*Lc-1)*La*Lb*Lc
         CASE (3,4,5)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            Lha = WedgeH(local(1),w)
            Lhc = WedgeH(local(4),w)
            value = La*Lb*varPhi(i,Lb-La)*Phi(j,Lhc-Lha)
         CASE DEFAULT
            CALL Fatal('PElementBase::WedgeFacePBasis','Unknown face for wedge')
      END SELECT
    END FUNCTION SD_WedgeFacePBasis


!------------------------------------------------------------------------------
!>     Gradient of wedge face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_dWedgeFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of wedges face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      REAL(KIND=dp) :: La, Lb, Lc, Lha, Lhc, Legi, Legj, parW 
      REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLha, dLhc, dW, grad
      INTEGER :: local(4)

      ! If local numbering not present, use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getWedgeFaceMap(face)
      ! Numbering was present, use it
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF
      
      ! Set sign of w and its derivative for faces 1 and 2
      dW = 0
      SELECT CASE (face)
      CASE (1)
         parW = -w
         dW(3) = -1
      CASE (2)
         parW = w
         dW(3) = 1
      END SELECT

      ! Get value of face function
      grad = 0
      SELECT CASE(face)
         CASE (1,2)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            Lc = WedgeL(local(3),u,v)
            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            dLc = dWedgeL(local(3),u,v)
            
            ! Precalculate values of legenre functions
            Legi = LegendreP(i,Lb-La)
            Legj = LegendreP(j,2d0*Lc-1)

            ! Get value of gradient
            grad = 1d0/2*dLa*Lb*Lc*Legi*Legj*(1+parW)+&
                 1d0/2*La*dLb*Lc*Legi*Legj*(1+parW) +&
                 1d0/2*La*Lb*dLc*Legi*Legj*(1+parW) +&
                 1d0/2*La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb-dLa)*Legj*(1+parW) +&
                 1d0/2*La*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*(2*dLc)*(1+parW) +&
                 1d0/2*La*Lb*Lc*Legi*Legj*dW
         CASE (3,4,5)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            
            Lha = WedgeH(local(1),w)
            Lhc = WedgeH(local(4),w)
            dLha = dWedgeH(local(1),w)
            dLhc = dWedgeH(local(4),w)

            Legi = varPhi(i,Lb-La)
            Legj = Phi(j,Lhc-Lha)

            grad = dLa*Lb*Legi*Legj + La*dLb*Legi*Legj + &
                 La*Lb*dVarPhi(i,Lb-La)*(dLb-dLa)*Legj + &
                 La*Lb*Legi*dPhi(j,Lhc-Lha)*(dLhc-dLha)
         CASE DEFAULT
            CALL Fatal('PElementBase::dWedgeFacePBasis','Unknown face for wedge')
      END SELECT
    END FUNCTION SD_dWedgeFacePBasis


!------------------------------------------------------------------------------
!>    Wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
!>     2nd derivatives of wedge face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_ddWedgeFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of wedges face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLha, dLhc, dW
      REAL(KIND=dp) :: La, Lb, Lc, Lha, Lhc, Legi, Legj, parW, grad(3,3)
      INTEGER :: p,q,local(4)
      REAL(KIND=dp) :: f(6), df(6,3), ddf(6,3,3), ddPi, ddPj

      ! If local numbering not present, use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getWedgeFaceMap(face)
      ! Numbering was present, use it
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF
      
      ! Set sign of w and its derivative for faces 1 and 2
      dW = 0
      SELECT CASE (face)
      CASE (1)
         parW = -w
         dW(3) = -1
      CASE (2)
         parW = w
         dW(3) = 1
      END SELECT

      ! Get value of face function
      grad = 0
      SELECT CASE(face)
         CASE (1,2)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            Lc = WedgeL(local(3),u,v)

            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            dLc = dWedgeL(local(3),u,v)
            
            ! Precalculate values of legenre functions
            Legi = LegendreP(i,Lb-La)
            Legj = LegendreP(j,2d0*Lc-1)

            ! Get value of gradient
!           grad = 1d0/2*dLa*Lb*Lc*Legi*Legj*(1+parW)+&
!                1d0/2*La*dLb*Lc*Legi*Legj*(1+parW) +&
!                1d0/2*La*Lb*dLc*Legi*Legj*(1+parW) +&
!                1d0/2*La*Lb*Lc*dLegendreP(i,Lb-La)*(dLb-dLa)*Legj*(1+parW) +&
!                1d0/2*La*Lb*Lc*Legi*dLegendreP(j,2*Lc-1)*(2*dLc)*(1+parW) +&
!                1d0/2*La*Lb*Lc*Legi*Legj*dW


            f(1) = La; f(2)=Lb; f(3)=Lc; f(4)=Legi; f(5)=Legj; f(6)=1+parW
            df(1,:) = dLa
            df(2,:) = dLb
            df(3,:) = dLc
            df(4,:) = dLegendreP(i,Lb-La)*(dLb-dLa)
            df(5,:) = dLegendreP(j,2*Lc-1)*2*dLc
            df(6,:) = dW

            ddf = 0
            ddPi = ddLegendreP(i,Lb-La)
            ddPj = ddLegendreP(j,2*Lc-1)
            DO p=1,3
              DO q=p,3
                ddf(4,p,q) = ddPi*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
                ddf(5,p,q) = ddPj*4*dLc(p)*dLc(q)
              END DO
            END DO

            grad = Product2ndDerivatives(6,f,df,ddf,3,0)/2
            grad(2,1) = grad(1,2)
            grad(3,1) = grad(1,3)
            grad(3,2) = grad(2,3)
         CASE (3,4,5)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            
            Lha = WedgeH(local(1),w)
            Lhc = WedgeH(local(4),w)
            dLha = dWedgeH(local(1),w)
            dLhc = dWedgeH(local(4),w)

            Legi = varPhi(i,Lb-La)
            Legj = Phi(j,Lhc-Lha)

!           grad = dLa*Lb*Legi*Legj + La*dLb*Legi*Legj + &
!                La*Lb*dVarPhi(i,Lb-La)*(dLb-dLa)*Legj + &
!                La*Lb*Legi*dPhi(j,Lhc-Lha)*(dLhc-dLha)

           f(1) = La; f(2)=Lb; f(3)=Legi; f(4)=Legj
           df(1,:) = dLa
           df(2,:) = dLb
           df(3,:) = dVarPhi(i,Lb-La)*(dLb-dLa)
           df(4,:) = dPhi(j,Lhc-Lha)*(dLhc-dLha)
           ddPi = ddVarPhi(i,Lb-La)
           ddPj = ddPhi(j,Lhc-Lha)
           ddf = 0
           DO p=1,3
             DO q=p,3
               ddf(3,p,q) = ddPi*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
               ddf(4,p,q) = ddPj*(dLhc(p)-dLha(p))*(dLhc(q)-dLha(q))
             END DO
           END DO

           grad = Product2ndDerivatives(4,f,df,ddf,3,0)
           grad(2,1) = grad(1,2)
           grad(3,1) = grad(1,3)
           grad(3,2) = grad(2,3)

         CASE DEFAULT
            CALL Fatal('PElementBase::dWedgeFacePBasis','Unknown face for wedge')
      END SELECT
    END FUNCTION SD_ddWedgeFacePBasis


!------------------------------------------------------------------------------
!>    Wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_WedgeBubblePBasis(i,j,k,u,v,w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,2),(0,0,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i,j,k
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: L1,L2,L3,value
      
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      ! Get value of bubble function
      value = L1*L2*L3*LegendreP(i,L2-L1)*LegendreP(j,2d0*L3-1)*Phi(k,w)
    END FUNCTION SD_WedgeBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_dWedgeBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,2),(0,0,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i,j,k
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: L1,L2,L3,Legi,Legj,phiW
      REAL(KIND=dp), DIMENSION(3) :: dL1, dL2, dL3, dW, grad
      
      ! Initialize derivative of w
      dW = [ 0,0,1 ]
      ! Values of function L
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      dL1 = dWedgeL(1,u,v)
      dL2 = dWedgeL(2,u,v)
      dL3 = dWedgeL(3,u,v)

      Legi = LegendreP(i,L2-L1)
      Legj = LegendreP(j,2d0*L3-1)
      phiW = Phi(k,w)

      grad = dL1*L2*L3*Legi*Legj*phiW + L1*dL2*L3*Legi*Legj*phiW +&
           L1*L2*dL3*Legi*Legj*phiW + L1*L2*L3*dLegendreP(i,L2-L1)*(dL2-dL1)*Legj*phiW +&
           L1*L2*L3*Legi*dLegendreP(j,2d0*L3-1)*(2d0*dL3)*phiW +&
           L1*L2*L3*Legi*Legj*dPhi(k,w)*dW
    END FUNCTION SD_dWedgeBubblePBasis

!------------------------------------------------------------------------------
!>    2nd derivatives of wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION SD_ddWedgeBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,2),(0,0,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i,j,k
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: dL1, dL2, dL3, dW
      REAL(KIND=dp) :: L1,L2,L3,Legi,Legj,phiW,grad(3,3)

      INTEGER :: p,q
      REAL(KIND=dp) :: f(6), df(6,3), ddf(6,3,3), ddPi, ddPj, ddPk
      
      ! Initialize derivative of w
      dW = [ 0,0,1 ]
      ! Values of function L
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      dL1 = dWedgeL(1,u,v)
      dL2 = dWedgeL(2,u,v)
      dL3 = dWedgeL(3,u,v)

      Legi = LegendreP(i,L2-L1)
      Legj = LegendreP(j,2d0*L3-1)
      phiW = Phi(k,w)

!     grad = dL1*L2*L3*Legi*Legj*phiW + L1*dL2*L3*Legi*Legj*phiW +&
!          L1*L2*dL3*Legi*Legj*phiW + L1*L2*L3*dLegendreP(i,L2-L1)*(dL2-dL1)*Legj*phiW +&
!          L1*L2*L3*Legi*dLegendreP(j,2d0*L3-1)*(2d0*dL3)*phiW +&
!          L1*L2*L3*Legi*Legj*dPhi(k,w)*dW

      f(1)=L1; f(2)=L2; f(3)=L3; f(4)=Legi; f(5)=Legj; f(6)=phiW
      df(1,:) = dL1
      df(2,:) = dL2
      df(3,:) = dL3
      df(4,:) = dLegendreP(i,L2-L1)*(dL2-dL1)
      df(5,:) = dLegendreP(j,2*L3-1)*2*dL3
      df(6,:) = dPhi(k,w)*dW

      ddPi = ddLegendreP(i,L2-L1)
      ddPj = ddLegendreP(j,2*L3-1)
      ddPk = ddPhi(k,w)
      ddf = 0
      DO p=1,3
        DO q=p,3
          ddf(4,p,q) = ddPi*(dL2(p)-dL1(p))*(dL2(q)-dL1(q))
          ddf(5,p,q) = ddPj*4*dL3(p)*dL3(q)
          ddf(6,p,q) = ddPk*dW(p)*dW(q)
        END DO
      END DO

      grad = Product2ndDerivatives(6,f,df,ddf,3,0)
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)

    END FUNCTION SD_ddWedgeBubblePBasis

! --- end serendipity wedge

    
!------------------------------------------------------------------------------      
!>     Wedge edge basis at point (u,v,w)
!------------------------------------------------------------------------------      
    FUNCTION WedgeEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(value)
!------------------------------------------------------------------------------      
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of wedges edge function to calculate
!        edge = {1,2,..,9}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges edge function i at point (u,v,w), i.e.
!       value = N_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp) :: La, Lb, tmp, value, Pa, Pb, PhiPar
      LOGICAL :: invert
      INTEGER :: local(2)
      
      ! Edge is not inverted by default 
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      local = getWedgeEdgeMap(edge)
      Pa = WedgeNodalPBasis(local(1),u,v,w)
      Pb = WedgeNodalPBasis(local(2),u,v,w)

      SELECT CASE(edge)
      CASE (1,4)
         La = WedgeL(1,u,v)
         Lb = WedgeL(2,u,v)
         PhiPar = Lb-La
      CASE (2,5)
         La = WedgeL(2,u,v)
         Lb = WedgeL(3,u,v)
         PhiPar = Lb-La
      CASE (3,6)
         La = WedgeL(3,u,v)
         Lb = WedgeL(1,u,v)
         PhiPar = Lb-La
      CASE (7,8,9)
         PhiPar = w
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeEdgePBasis','Unknown edge for wedge')
      END SELECT

      ! Swap parameters for inverted edges
      IF (invert) PhiPar = -PhiPar

      ! Calculate value from general form for edges
      value = Pa*Pb*varPhi(i,PhiPar)
    END FUNCTION WedgeEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of wedge edge basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dWedgeEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of wedges edge function to calculate
!        edge = {1,2,..,9}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp) :: parW, La, Lb, phiI, tmp, PhiPar, dPhiPar(3), dPhiI(3), &
           grad(3), dLa(3), dLb(3), dW(3), dtmp(3), Pa, Pb, dPa(3), dPb(3) 
      LOGICAL :: invert 
      INTEGER :: local(2)
      
      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      local = getWedgeEdgeMap(edge)

      Pa  = WedgeNodalPBasis(local(1),u,v,w)
      Pb  = WedgeNodalPBasis(local(2),u,v,w)

      dPa = dWedgeNodalPBasis(local(1),u,v,w)
      dPb = dWedgeNodalPBasis(local(2),u,v,w)

      SELECT CASE(edge)
      CASE (1,4)
         La  = WedgeL(1,u,v)
         Lb  = WedgeL(2,u,v)
         dLa = dWedgeL(1,u,v)
         dLb = dWedgeL(2,u,v)
         PhiPar = Lb-La; dPhiPar = dLb-dLa
      CASE (2,5)
         La  = WedgeL(2,u,v)
         Lb  = WedgeL(3,u,v)
         dLa = dWedgeL(2,u,v)
         dLb = dWedgeL(3,u,v)
         PhiPar = Lb-La; dPhiPar = dLb-dLa
      CASE (3,6)
         La  = WedgeL(3,u,v)
         Lb  = WedgeL(1,u,v)
         dLa = dWedgeL(3,u,v)
         dLb = dWedgeL(1,u,v)
         PhiPar = Lb-La; dPhiPar = dLb-dLa
      CASE (7,8,9)
         PhiPar = w; dPhiPar = [0,0,1];
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeEdgePBasis','Unknown edge for wedge')
      END SELECT

      ! Swap parameters for inverted edges
      IF (invert) THEN
         PhiPar = -PhiPar; dPhiPar=-dPhiPar
      END IF      
      
      phiI  = varPhi(i,PhiPar)
      dphiI = dvarPhi(i,PhiPar)*dPhiPar

      ! Calculate value of function from general form
      grad = dPa*Pb*PhiI + Pa*dPb*PhiI + Pa*Pb*dPhiI
    END FUNCTION dWedgeEdgePBasis


!------------------------------------------------------------------------------
!>     2nd derivatives of wedge edge basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddWedgeEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of wedges edge function to calculate
!        edge = {1,2,..,9}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp) :: parW, Pa, Pb, dPa(3), dPb(3), ddPa(3,3), ddPb(3,3), &
           La, Lb, phiI, tmp, grad(3,3), dLa(3), dLb(3), dW(3), dtmp(3), &
           PhiPar, dPhiPar(3), dPhiI(3), ddPhiI(3,3), s
      LOGICAL :: invert 
      
      INTEGER :: p,q, local(2)

      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      local = getWedgeEdgeMap(edge)

      Pa  = WedgeNodalPBasis(local(1),u,v,w)
      Pb  = WedgeNodalPBasis(local(2),u,v,w)

      dPa = dWedgeNodalPBasis(local(1),u,v,w)
      dPb = dWedgeNodalPBasis(local(2),u,v,w)

      ddPa = ddWedgeNodalPBasis(local(1),u,v,w)
      ddPb = ddWedgeNodalPBasis(local(2),u,v,w)

      SELECT CASE(edge)
      CASE (1,4)
         La  = WedgeL(1,u,v)
         Lb  = WedgeL(2,u,v)
         dLa = dWedgeL(1,u,v)
         dLb = dWedgeL(2,u,v)
         PhiPar = Lb-La; dPhiPar = dLb-dLa
      CASE (2,5)
         La  = WedgeL(2,u,v)
         Lb  = WedgeL(3,u,v)
         dLa = dWedgeL(2,u,v)
         dLb = dWedgeL(3,u,v)
         PhiPar = Lb-La; dPhiPar = dLb-dLa
      CASE (3,6)
         La  = WedgeL(3,u,v)
         Lb  = WedgeL(1,u,v)
         dLa = dWedgeL(3,u,v)
         dLb = dWedgeL(1,u,v)
         PhiPar = Lb-La; dPhiPar = dLb-dLa
      CASE (7,8,9)
         PhiPar = w; dPhiPar = [0,0,1]
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeEdgePBasis','Unknown edge for wedge')
      END SELECT

      ! Swap parameters for inverted edges
      IF (invert) THEN
        PhiPar = -PhiPar; dPhiPar = -dPhiPar
      END IF      
      
      phiI  = varPhi(i,PhiPar)
      dphiI = dvarPhi(i,PhiPar)*dPhiPar

      s = ddVarPhi(i,PhiPar)
      DO p=1,3
        DO q=p,3
          ddphiI(p,q) = s*dPhiPar(p)*dPhiPar(q)
        END DO
      END DO

      ! Calculate value of function from general form
!     grad = dPa*Pb*PhiI + Pa*dPb*PhiI + Pa*Pb*dPhiI

      BLOCK
        REAL(KIND=dp) :: f(3), df(3,3), ddf(3,3,3)

        f = [Pa,Pb,PhiI]
        df(1,:) = dPa
        df(2,:) = dPb
        df(3,:) = dPhiI
        ddf(1,:,:) = ddPa
        ddf(2,:,:) = ddPb
        ddf(3,:,:) = ddPhiI

        grad = Product2ndDerivatives(3,f,df,ddf,3,0)
      END BLOCK

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)

    END FUNCTION ddWedgeEdgePBasis


!------------------------------------------------------------------------------
!>     Wedge face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION WedgeFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of wedges face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges face function m(i,j) at point (u,v,w), i.e.
!       value = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      REAL(KIND=dp) :: La, Lb, Lc, Ld, Lha, Lhc, parW, value, Pa, Pb
      INTEGER :: local(4)

      ! If local numbering not present, use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
        local(1:4) = getWedgeFaceMap(face)
      ! Numbering was present. Use it
      ELSE
        local(1:4) = localNumbers(1:4)
      END IF


      ! Set sign of w for faces 1 and 2
      SELECT CASE (face)
      CASE (1)
         parW = -w
      CASE (2)
         parW =  w
      END SELECT

      ! Get value of face function
      value = 0
      SELECT CASE(face)
      CASE (1,2)
        La = WedgeL(local(1),u,v)
        Lb = WedgeL(local(2),u,v)
        Lc = WedgeL(local(3),u,v)
        Ld = (1+parW)/2
        value = LegendreP(i,Lb-La)*LegendreP(j,2*Lc-1)*La*Lb*Lc*Ld
      CASE (3,4,5)
        Pa = WedgeNodalPBasis(local(1),u,v,w)
        Pb = WedgeNodalPBasis(local(3),u,v,w)

        La  = WedgeL(local(1),u,v)
        Lb  = WedgeL(local(2),u,v)

        Lha = WedgeH(local(1),w)
        Lhc = WedgeH(local(4),w)

        value = Pa*Pb*LegendreP(i,Lb-La)*LegendreP(j,Lhc-Lha)
      CASE DEFAULT
        CALL Fatal('PElementBase::WedgeFacePBasis','Unknown face for wedge')
      END SELECT
    END FUNCTION WedgeFacePBasis


!------------------------------------------------------------------------------
!>     Gradient of wedge face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dWedgeFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of wedges face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      REAL(KIND=dp) :: La, Lb, Lc, Lha, Lhc, Ld, Legi, Legj, parW ,& 
                 Pa,Pb,dPa(3),dPb(3)
      REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLd,  &
                dLha, dLhc, dW, dLegi, dLegj, grad
      INTEGER :: local(4)

      ! If local numbering not present, use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getWedgeFaceMap(face)
      ! Numbering was present, use it
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF
      
      ! Set sign of w and its derivative for faces 1 and 2
      dW = 0
      SELECT CASE (face)
      CASE (1)
         parW = -w
         dW(3) = -1
      CASE (2)
         parW = w
         dW(3) = 1
      END SELECT

      ! Get value of face function
      grad = 0
      SELECT CASE(face)
         CASE (1,2)
            La  = WedgeL(local(1),u,v)
            Lb  = WedgeL(local(2),u,v)
            Lc  = WedgeL(local(3),u,v)
            Ld  = (1+ParW)/2

            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            dLc = dWedgeL(local(3),u,v)
            dLd = dW/2
            
            ! Precalculate values of legenre functions
            Legi = LegendreP(i,Lb-La)
            Legj = LegendreP(j,2*Lc-1)

            dLegi = dLegendreP(i,Lb-La)*(dLb-dLa)
            dLegj = dLegendreP(j,2*Lc-1)*2*dLc

            ! Get value of gradient
            grad = dLa*Lb*Lc*Ld*Legi*Legj + La*dLb*Lc*Ld*Legi*Legj + &
                   La*Lb*dLc*Ld*Legi*Legj + La*Lb*Lc*dLd*Legi*Legj + &
                   La*Lb*Lc*Ld*dLegI*Legj + La*lb*Lc*Ld*Legi*dLegj
         CASE (3,4,5)
            La  = WedgeL(local(1),u,v)
            Lb  = WedgeL(local(2),u,v)

            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            
            Lha  = WedgeH(local(1),w)
            Lhc  = WedgeH(local(4),w)

            dLha = dWedgeH(local(1),w)
            dLhc = dWedgeH(local(4),w)

            Pa  = WedgeNodalPBasis(local(1),u,v,w)
            Pb  = WedgeNodalPBasis(local(3),u,v,w)

            dPa = dWedgeNodalPBasis(local(1),u,v,w)
            dPb = dWedgeNodalPBasis(local(3),u,v,w)

            Legi = LegendreP(i,Lb-La)
            Legj = LegendreP(j,Lhc-Lha)

            dLegi = dLegendreP(i,Lb-La)*(dLb-dLa)
            dLegj = dLegendreP(j,Lhc-Lha)*(dLhc-dLha)

            grad = dPa*Pb*Legi*Legj + Pa*dPb*Legi*Legj + &
                   Pa*Pb*dLegi*Legj + Pa*Pb*Legi*dLegj
         CASE DEFAULT
            CALL Fatal('PElementBase::dWedgeFacePBasis','Unknown face for wedge')
      END SELECT
    END FUNCTION dWedgeFacePBasis


!------------------------------------------------------------------------------
!>    Wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
!>     2nd derivatives of wedge face basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddWedgeFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of wedges face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, DIMENSION(4), OPTIONAL :: localNumbers
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dLha, dLhc, dW
      REAL(KIND=dp) :: La, Lb, Lc, Lha, Lhc, Legi, Legj, parW, grad(3,3)
      REAL(KIND=dp) :: ddPi, ddPj, dLegi(3), dLegj(3), ddLegi(3,3), ddLegj(3,3)
      REAL(KIND=dp) :: Pa, Pb, dPa(3), dPb(3), ddPa(3,3), ddPb(3,3), s,t,ds(3),dt(3)
      INTEGER :: p,q,local(4)

      ! If local numbering not present, use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getWedgeFaceMap(face)
      ! Numbering was present, use it
      ELSE
         local(1:4) = localNumbers(1:4)
      END IF
      
      ! Set sign of w and its derivative for faces 1 and 2
      dW = 0
      SELECT CASE (face)
      CASE (1)
         parW = -w
         dW(3) = -1
      CASE (2)
         parW = w
         dW(3) = 1
      END SELECT

      ! Get value of face function
      grad = 0
      SELECT CASE(face)
         CASE (1,2)
            La = WedgeL(local(1),u,v)
            Lb = WedgeL(local(2),u,v)
            Lc = WedgeL(local(3),u,v)

            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            dLc = dWedgeL(local(3),u,v)
            
            ! Precalculate values of legenre functions
            Legi = LegendreP(i,Lb-La)
            Legj = LegendreP(j,2*Lc-1)

            ! Get value of gradient
           BLOCK
              REAL(KIND=dp) :: f(6), df(6,3), ddf(6,3,3)

              f = [La, Lb, Lc, Legi, Legj, 1+parW]
              df(1,:) = dLa
              df(2,:) = dLb
              df(3,:) = dLc
              df(4,:) = dLegendreP(i,Lb-La)*(dLb-dLa)
              df(5,:) = dLegendreP(j,2*Lc-1)*2*dLc
              df(6,:) = dW

              ddPi = ddLegendreP(i,Lb-La)
              ddPj = ddLegendreP(j,2*Lc-1)
              ddf = 0
              DO p=1,3
                DO q=p,3
                  ddf(4,p,q) = ddPi*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
                  ddf(5,p,q) = ddPj*4*dLc(p)*dLc(q)
                END DO
              END DO
              grad = Product2ndDerivatives(6,f,df,ddf,3,0)/2
            END BLOCK
            grad(2,1) = grad(1,2)
            grad(3,1) = grad(1,3)
            grad(3,2) = grad(2,3)
         CASE (3,4,5)
            La  = WedgeL(local(1),u,v)
            Lb  = WedgeL(local(2),u,v)

            dLa = dWedgeL(local(1),u,v)
            dLb = dWedgeL(local(2),u,v)
            
            Lha = WedgeH(local(1),w)
            Lhc = WedgeH(local(4),w)

            dLha = dWedgeH(local(1),w)
            dLhc = dWedgeH(local(4),w)

            Pa  = WedgeNodalPBasis(local(1),u,v,w)
            Pb  = WedgeNodalPBasis(local(3),u,v,w)

            dPa = dWedgeNodalPBasis(local(1),u,v,w)
            dPb = dWedgeNodalPBasis(local(3),u,v,w)

            ddPa = ddWedgeNodalPBasis(local(1),u,v,w)
            ddPb = ddWedgeNodalPBasis(local(3),u,v,w)

            s  = Lb-La
            ds = dLb-dLa

            t  = Lhc-Lha
            dt = dLhc-dLha

            Legi = LegendreP(i,s)
            Legj = LegendreP(j,t)

            dLegi = dLegendreP(i,s)*ds
            dLegj = dLegendreP(j,t)*dt

            ddPi = ddLegendreP(i,s)
            ddPj = ddLegendreP(j,t)
            DO p=1,3
              DO q=p,3
                ddLegi(p,q) = ddPi*ds(p)*ds(q)
                ddLegj(p,q) = ddPj*dt(p)*dt(q)
              END DO
            END DO

            BLOCK
              REAL(KIND=dp) :: f(4), df(4,3), ddf(4,3,3)

              f = [Pa, Pb, Legi, Legj]
              df(1,:) = dPa
              df(2,:) = dPb
              df(3,:) = dLegi
              df(4,:) = dLegj
              ddf(1,:,:) = ddPa
              ddf(2,:,:) = ddPb
              ddf(3,:,:) = ddLegi
              ddf(4,:,:) = ddLegj

              grad = Product2ndDerivatives(4,f,df,ddf,3,0)
            END BLOCK

            grad(2,1) = grad(1,2)
            grad(3,1) = grad(1,3)
            grad(3,2) = grad(2,3)

         CASE DEFAULT
            CALL Fatal('PElementBase::dWedgeFacePBasis','Unknown face for wedge')
      END SELECT
    END FUNCTION ddWedgeFacePBasis


!------------------------------------------------------------------------------
!>    Wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION WedgeBubblePBasis(i,j,k,u,v,w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,2),(0,0,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of wedges bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i,j,k
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: L1,L2,L3,L4,value,s,t
      
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      s = L2-L1
      t = 2*L3-1

      ! Get value of bubble function
      value = L1*L2*L3*LegendreP(i,s)*LegendreP(j,t)*Phi(k+2,w)
    END FUNCTION WedgeBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dWedgeBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,2),(0,0,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i,j,k
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: dL1, dL2, dL3, grad
      REAL(KIND=dp) :: L1,L2,L3,Legi,Legj,Legk,dLegi(3),dLegj(3),dLegk(3), &
                       s,t,ds(3),dt(3), L4, dL4(3), value
      
      ! Values of function L
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      dL1 = dWedgeL(1,u,v)
      dL2 = dWedgeL(2,u,v)
      dL3 = dWedgeL(3,u,v)

      s = L2-L1
      ds = dL2-dL1

      t = 2*L3-1
      dt = 2*dL3

      Legi = LegendreP(i,s)
      Legj = LegendreP(j,t)
      Legk = Phi(k+2,w)

      dLegi = dLegendreP(i,s)*ds
      dLegj = dLegendreP(j,t)*dt
      dLegk = dPhi(k+2,w)*[0,0,1]

      grad = dL1*L2*L3*Legi*Legj*Legk + L1*dL2*L3*Legi*Legj*Legk + &
             L1*L2*dL3*Legi*Legj*Legk + L1*L2*L3*dLegi*Legj*Legk + &
             L1*L2*L3*Legi*dLegj*Legk + L1*L2*L3*Legi*Legj*dLegk
    END FUNCTION dWedgeBubblePBasis


!------------------------------------------------------------------------------
!>    2nd derivatives of wedge bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddWedgeBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,2),(0,0,3),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of wedges bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: i,j,k
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: dL1, dL2, dL3, dW
      REAL(KIND=dp) :: L1,L2,L3,Legi,Legj,Legk,phiW,grad(3,3)

      INTEGER :: p,q
      REAL(KIND=dp) :: ddPi, ddPj, ddPk
      REAL(KIND=dp) :: s,t,ds(3),dt(3), dLegi(3),dLegj(3),dLegk(3)
      
      ! Initialize derivative of w
      dW = [0,0,1]

      ! Values of function L
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      dL1 = dWedgeL(1,u,v)
      dL2 = dWedgeL(2,u,v)
      dL3 = dWedgeL(3,u,v)

      s = L2-L1
      ds = dL2-dL1

      t = 2*L3-1
      dt = 2*dL3

      Legi = LegendreP(i,s)
      Legj = LegendreP(j,t)
      Legk = Phi(k+2,w)

      dLegi = dLegendreP(i,s)*ds
      dLegj = dLegendreP(j,t)*dt
      dLegk = dPhi(k+2,w)*dw

      BLOCK
        REAL(KIND=dp) :: f(6), df(6,3), ddf(6,3,3)

        f=[L1,L2,L3,Legi,Legj,Legk]
        df(1,:) = dL1
        df(2,:) = dL2
        df(3,:) = dL3
        df(4,:) = dLegi
        df(5,:) = dLegj
        df(6,:) = dLegk

        ddPi = ddLegendreP(i,s)
        ddPj = ddLegendreP(j,t)
        ddPk = ddPhi(k+2,w)
        ddf = 0
        DO p=1,3
          DO q=p,3
            ddf(4,p,q) = ddPi*ds(p)*ds(q)
            ddf(5,p,q) = ddPj*dt(p)*dt(q)
            ddf(6,p,q) = ddPk*dw(p)*dw(q)
          END DO
        END DO
        grad = Product2ndDerivatives(6,f,df,ddf,3,0)
      END BLOCK

      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)

    END FUNCTION ddWedgeBubblePBasis


    PURE FUNCTION WedgeL(which, u, v) RESULT(value)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v
      ! Result
      REAL(KIND=dp) :: value 
      
      value = 0
      SELECT CASE(which)
      CASE (1,4)
         value = 1d0/2*(1d0-u-v/SQRT(3d0))
      CASE (2,5)
         value = 1d0/2*(1d0+u-v/SQRT(3d0))
      CASE (3,6)
         value = SQRT(3d0)/3*v
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeL','Unknown linear function L for wedge')
#endif
      END SELECT
    END FUNCTION WedgeL

    PURE FUNCTION WedgeH(which, w) RESULT(value)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: w
      ! Result
      REAL(KIND=dp) :: value 
      
      value = 0
      SELECT CASE(which)
      CASE (1,2,3)
         value = -w/2
      CASE (4,5,6)
         value = w/2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeH','Unknown linear function H for wedge')
#endif
      END SELECT
    END FUNCTION WedgeH

    PURE FUNCTION dWedgeL(which, u, v) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v
      ! Result
      REAL(KIND=dp) :: grad(3) 

      grad = 0 
      SELECT CASE(which)
      CASE (1,4)
         grad(1) = -1d0/2
         grad(2) = -SQRT(3d0)/6
      CASE (2,5)
         grad(1) = 1d0/2
         grad(2) = -SQRT(3d0)/6
      CASE (3,6)
         grad(2) = SQRT(3d0)/3
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeL','Unknown linear function dL for wedge')
#endif
      END SELECT
    END FUNCTION dWedgeL

    PURE FUNCTION dWedgeH(which, w) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: w
      ! Result
      REAL(KIND=dp) :: grad(3) 

      grad = 0 
      SELECT CASE(which)
      CASE (1,2,3)
         grad(3) = -1d0/2
      CASE (4,5,6)
         grad(3) = 1d0/2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeH','Unknown linear function dH for wedge')
#endif
      END SELECT
    END FUNCTION dWedgeH


!------------------------------------------------------------------------------
!>     Pyramid nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION PyramidNodalPBasis(node, u, v, w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of wedge nodal function to calculate
!        node = {1,2,3,4,5}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of pyramids nodal function at point (u,v,w), i.e.
!       value = N_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: value, s, sq2 = SQRT(2.0_dp)

      s = w/sq2
      SELECT CASE(node)
      CASE (1)
         value = (1-u-v-s+u*v/(1-s))/4
      CASE (2)
         value = (1+u-v-s-u*v/(1-s))/4
      CASE (3)
         value = (1+u+v-s+u*v/(1-s))/4
      CASE (4)
         value = (1-u+v-s-u*v/(1-s))/4
      CASE (5)
         value = s
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidNodalPBasis','Unknown node for pyramid')
      END SELECT
    END FUNCTION PyramidNodalPBasis


!------------------------------------------------------------------------------
!>     Gradient of pyramids nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dPyramidNodalPBasis(node, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of pyramid nodal function to calculate
!        node = {1,2,3,4,5}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids nodal function at point (u,v,w), i.e.
!       grad = dN_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: grad
      REAL(KIND=dp) :: s, sq2=SQRT(2.0_dp)


      grad = 0
      s = w/sq2
      SELECT CASE(node)
      CASE (1)
         grad(1) = -1 + v/(1-s)
         grad(2) = -1 + u/(1-s)
         grad(3) = -1 + u*v/(1-s)**2
      CASE (2)
         grad(1) =  1 - v/(1-s)
         grad(2) = -1 - u/(1-s)
         grad(3) = -1 - u*v/(1-s)**2
      CASE (3)
         grad(1) =  1 + v/(1-s)
         grad(2) =  1 + u/(1-s)
         grad(3) = -1 + u*v/(1-s)**2
      CASE (4)
         grad(1) = -1 - v/(1-s)
         grad(2) =  1 - u/(1-s)
         grad(3) = -1 - u*v/(1-s)**2
      CASE (5)
         grad(3) = 4
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidNodalPBasis','Unknown node for pyramid')
      END SELECT
      grad = grad / 4
      grad(3) = grad(3)/sq2
    END FUNCTION dPyramidNodalPBasis


!------------------------------------------------------------------------------
!>     Gradient of pyramids nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddPyramidNodalPBasis(node, u, v, w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: node
!      INPUT: number of pyramid nodal function to calculate
!        node = {1,2,3,4,5}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids nodal function at point (u,v,w), i.e.
!       grad = dN_i(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: s, sq2=SQRT(2.0_dp), grad(3,3)

      grad = 0
      s = w/sq2
      SELECT CASE(node)
      CASE (1)
         grad(1,1) = 0
         grad(1,2) = 1/(1-s)
         grad(1,3) = v/(1-s)**2

         grad(2,1) = 1/(1-s)
         grad(2,2) = 0
         grad(2,3) = u/(1-s)**2

         grad(3,1) = v/(1-s)**2
         grad(3,2) = u/(1-s)**2
         grad(3,3) = 2*u*v/(1-s)**3
      CASE (2)
         grad(1,1) =  0
         grad(1,2) =  -1/(1-s)
         grad(1,3) =  -v/(1-s)**2

         grad(2,1) = -1/(1-s)
         grad(2,2) = 0
         grad(2,3) = -u/(1-s)**2

         grad(3,1) = -v/(1-s)**2
         grad(3,2) = -u/(1-s)**2
         grad(3,3) = -2*u*v/(1-s)**3
      CASE (3)
         grad(1,1) = 0
         grad(1,2) = 1/(1-s)
         grad(1,3) = v/(1-s)**2

         grad(2,1) = 1/(1-s)
         grad(2,2) = 0
         grad(2,3) = u/(1-s)**2

         grad(3,1) = v/(1-s)**2
         grad(3,2) = u/(1-s)**2
         grad(3,3) = 2*u*v/(1-s)**3
      CASE (4)
         grad(1,1) = 0
         grad(1,2) = -1/(1-s)
         grad(1,3) = -v/(1-s)**2

         grad(2,1) = -1/(1-s)
         grad(2,2) =  0
         grad(2,3) = -u/(1-s)**2

         grad(3,1) = -v/(1-s)**2
         grad(3,2) = -u/(1-s)**2
         grad(3,3) = -2*u*v/(1-s)**3
      CASE (5)
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidNodalPBasis','Unknown node for pyramid')
      END SELECT
      grad = grad / 4
      grad(3,:) = grad(3,:)/sq2
      grad(:,3) = grad(:,3)/sq2
    END FUNCTION ddPyramidNodalPBasis



!-----------------------------------------------------------------------------
!>     Pyramid edge basis at point (u,v,w)
!-----------------------------------------------------------------------------
    FUNCTION PyramidEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(value)
!-----------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of pyramids edge function to calculate
!        edge = {1,2,..,8}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of pyramids edge function i at point (u,v,w), i.e.
!       value = N_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, INTENT(IN), OPTIONAL :: invertEdge
      ! Variables
      INTEGER :: local(2)
      LOGICAL :: invert
      REAL(KIND=dp) :: La, Lb, Pa, Pb, phiPar, value, s, sq2=SQRT(2.0_dp)
            
      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      local = getPyramidEdgeMap(edge)

      Pa = PyramidNodalPBasis(local(1),u,v,w)
      Pb = PyramidNodalPBasis(local(2),u,v,w)

      SELECT CASE(Edge)
      CASE(1,2,3,4)
        La = PyramidL(local(1),u,v)
        Lb = PyramidL(local(2),u,v)
      CASE(5,6,7,8)
        La = PyramidTL(local(1),u,v,w)
        Lb = PyramidTL(local(2),u,v,w)
      END SELECT

      s = w/sq2
      ! Invert edge if needed
      PhiPar = Lb-La
      IF (invert) phiPar = -phiPar

      ! Calculate value of edge function
      value = Pa*Pb*varPhi(i,phiPar)
    END FUNCTION PyramidEdgePBasis


!-------------------------------------------------------------------------------
!>     Gradient of pyramid edge basis at point (u,v,w)
!-------------------------------------------------------------------------------
    FUNCTION dPyramidEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!-------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of pyramids edge function to calculate
!        edge = {1,2,..,8}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, INTENT(IN), OPTIONAL :: invertEdge
      ! Variables
      REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dPa, dPb, dPhiPar, grad
      REAL(KIND=dp) :: La, Lb, Pa, Pb, phiPar, vPhiI, sq2=SQRT(2.0_dp)
      LOGICAL :: invert
      INTEGER :: local(2)

      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      local = getPyramidEdgeMap(edge)

      Pa  = PyramidNodalPBasis(local(1), u,v,w)
      Pb  = PyramidNodalPBasis(local(2), u,v,w)
      dPa = dPyramidNodalPBasis(local(1), u,v,w)
      dPb = dPyramidNodalPBasis(local(2), u,v,w)

      SELECT CASE(edge)
      CASE(1,2,3,4)
        La = PyramidL(local(1),u,v)
        Lb = PyramidL(local(2),u,v)

        dLa = dPyramidL(local(1),u,v)
        dLb = dPyramidL(local(2),u,v)
      CASE(5,6,7,8)
        La = PyramidTL(local(1),u,v,w)
        Lb = PyramidTL(local(2),u,v,w)

        dLa = dPyramidTL(local(1),u,v,w)
        dLb = dPyramidTL(local(2),u,v,w)
      END SELECT

      ! Invert edge if needed
      PhiPar = Lb-La
      dPhiPar = dLb-dLa
      IF (invert) THEN
         phiPar  = -phiPar
         dPhiPar = -dPhiPar
      END IF

      vPhiI = varPhi(i,phiPar)
      grad = dPa*Pb*vPhiI + Pa*dPb*vPhiI + Pa*Pb*dVarPhi(i,phiPar)*dPhiPar
    END FUNCTION dPyramidEdgePBasis

!-------------------------------------------------------------------------------
!>     Gradient of pyramid edge basis at point (u,v,w)
!-------------------------------------------------------------------------------
    FUNCTION ddPyramidEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
!-------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: edge
!      INPUT: number of pyramids edge function to calculate
!        edge = {1,2,..,8}
!
!    INTEGER :: i
!      INPUT: index of edge function, i = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    LOGICAL, OPTIONAL :: invertEdge
!      INPUT: whether to invert edge or not. If this flag is set to true
!        edge changing parameter of edge function is varied from [1,-1] in
!        stead of usual [-1,1].
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids edge function i at point (u,v,w), i.e.
!       grad = dN_i^{edge}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: edge, i
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      LOGICAL, INTENT(IN), OPTIONAL :: invertEdge
      ! Variables
      INTEGER :: p,q, local(2)
      REAL(KIND=dp) :: La, Lb, Pa, Pb, phiPar, vPhiI
      REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dPa, dPb, dPhiPar, dvPhiI
      REAL(KIND=dp), DIMENSION(3,3) :: ddPa, ddPb, ddvPhiI, grad
      LOGICAL :: invert

      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      local = getPyramidEdgeMap(edge)

      Pa   = PyramidNodalPBasis(local(1),u,v,w)
      Pb   = PyramidNodalPBasis(local(2),u,v,w)

      dPa  = dPyramidNodalPBasis(local(1),u,v,w)
      dPb  = dPyramidNodalPBasis(local(2),u,v,w)

      ddPa = ddPyramidNodalPBasis(local(1),u,v,w)
      ddPb = ddPyramidNodalPBasis(local(2),u,v,w)

      SELECT CASE(edge)
      CASE(1,2,3,4)
        La  = PyramidL(local(1),u,v)
        Lb  = PyramidL(local(2),u,v)

        dLa = dPyramidL(local(1),u,v)
        dLb = dPyramidL(local(2),u,v)
      CASE(5,6,7,8)
        La  = PyramidTL(local(1),u,v,w)
        Lb  = PyramidTL(local(2),u,v,w)

        dLa = dPyramidTL(local(1),u,v,w)
        dLb = dPyramidTL(local(2),u,v,w)
      END SELECT
 
      ! Invert edge if needed
      PhiPar = Lb-La
      dPhiPar = dLb-dLa
      IF (invert) THEN
         phiPar  = -phiPar
         dPhiPar = -dPhiPar
      END IF

      vPhiI = varPhi(i,phiPar)
      dvPhiI = dvarPhi(i,phiPar)*dPhiPar
      DO p=1,3
        DO q=p,3
          ddvPhiI(p,q) = ddvarPhi(i,phiPar)*dPhiPar(p)*dPhiPar(q)
        END DO
      END DO

!     grad = dPa*Pb*vPhiI + Pa*dPb*vPhiI + Pa*Pb*dVarPhi(i,phiPar)*dPhiPar
#if 1
      grad = 0
      grad(1,1) = grad(1,1) + ddPa(1,1)*Pb*vPhiI+dPa(1)*dPb(1)*vPhiI+dPa(1)*Pb*dvPhiI(1)
      grad(1,1) = grad(1,1) + dPa(1)*dPb(1)*vPhiI+Pa*ddPb(1,1)*vPhiI+Pa*dPb(1)*dvPhiI(1)
      grad(1,1) = grad(1,1) + dPa(1)*Pb*dvPhiI(1)+Pa*dPb(1)*dvPhiI(1)+Pa*Pb*ddvPhiI(1,1)

      grad(1,2) = grad(1,2) + ddPa(1,2)*Pb*vPhiI+dPa(2)*dPb(1)*vPhiI+dPa(2)*Pb*dvPhiI(1)
      grad(1,2) = grad(1,2) + dPa(1)*dPb(2)*vPhiI+Pa*ddPb(1,2)*vPhiI+Pa*dPb(2)*dvPhiI(1)
      grad(1,2) = grad(1,2) + dPa(1)*Pb*dvPhiI(2)+Pa*dPb(1)*dvPhiI(2)+Pa*Pb*ddvPhiI(1,2)
      grad(2,1) = grad(1,2)

      grad(1,3) = grad(1,3) + ddPa(1,3)*Pb*vPhiI+dPa(3)*dPb(1)*vPhiI+dPa(3)*Pb*dvPhiI(1)
      grad(1,3) = grad(1,3) + dPa(1)*dPb(3)*vPhiI+Pa*ddPb(1,3)*vPhiI+Pa*dPb(3)*dvPhiI(1)
      grad(1,3) = grad(1,3) + dPa(1)*Pb*dvPhiI(3)+Pa*dPb(1)*dvPhiI(3)+Pa*Pb*ddvPhiI(1,3)
      grad(3,1) = grad(1,3)

      grad(2,2) = grad(2,2) + ddPa(2,2)*Pb*vPhiI+dPa(2)*dPb(2)*vPhiI+dPa(2)*Pb*dvPhiI(2)
      grad(2,2) = grad(2,2) + dPa(2)*dPb(2)*vPhiI+Pa*ddPb(2,2)*vPhiI+Pa*dPb(2)*dvPhiI(2)
      grad(2,2) = grad(2,2) + dPa(2)*Pb*dvPhiI(2)+Pa*dPb(2)*dvPhiI(2)+Pa*Pb*ddvPhiI(2,2)

      grad(2,3) = grad(2,3) + ddPa(2,3)*Pb*vPhiI+dPa(3)*dPb(2)*vPhiI+dPa(3)*Pb*dvPhiI(2)
      grad(2,3) = grad(2,3) + dPa(2)*dPb(3)*vPhiI+Pa*ddPb(2,3)*vPhiI+Pa*dPb(3)*dvPhiI(2)
      grad(2,3) = grad(2,3) + dPa(2)*Pb*dvPhiI(3)+Pa*dPb(2)*dvPhiI(3)+Pa*Pb*ddvPhiI(2,3)
      grad(3,2) = grad(2,3)

      grad(3,3) = grad(3,3) + ddPa(3,3)*Pb*vPhiI+dPa(3)*dPb(3)*vPhiI+dPa(3)*Pb*dvPhiI(3)
      grad(3,3) = grad(3,3) + dPa(3)*dPb(3)*vPhiI+Pa*ddPb(3,3)*vPhiI+Pa*dPb(3)*dvPhiI(3)
      grad(3,3) = grad(3,3) + dPa(3)*Pb*dvPhiI(3)+Pa*dPb(3)*dvPhiI(3)+Pa*Pb*ddvPhiI(3,3)
#else
      BLOCK
        REAL(KIND=dp) :: f(3), df(3,3), ddf(3,3,3)

        f(1)=Pa; f(2)=Pb; f(3)=vPhiI
        df(1,:)=dPa; df(2,:)=dPb; df(3,:)=dvPhiI
        ddf(1,:,:)=ddPa; ddf(2,:,:)=ddPb; ddf(3,:,:)=ddvPhiI

        grad = Product2ndDerivatives(3,f,df,ddf,3,0)
      END BLOCK
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
#endif

    END FUNCTION ddPyramidEdgePBasis

	
!-----------------------------------------------------------------------------	
!>     Pyramid face basis at point (u,v,w)
!-----------------------------------------------------------------------------	
    FUNCTION PyramidFacePBasis(face, i, j, u, v, w, localNumbers ) RESULT(value)
!-----------------------------------------------------------------------------	
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of pyramids face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of pyramids face function m(i,j) at point (u,v,w), i.e.
!       value = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)
      ! Variables
      REAL(KIND=dp) :: Pa, Pb, Pc, La, Lb, Lc, value, s, varg,swap
      INTEGER :: local(4)
      
      ! If local numbering not present use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getPyramidFaceMap(face)
      ELSE
         local(1:4) = localNumbers
      END IF

      SELECT CASE(face)
      CASE (1)
         Pa = PyramidNodalPBasis(local(1),u,v,w)
         Pb = PyramidNodalPBasis(local(3),u,v,w)

         La = PyramidL(local(1),u,v)
         Lb = PyramidL(local(2),u,v)
         Lc = PyramidL(local(4),u,v)
         value = Pa*Pb*LegendreP(i,Lb-La)*LegendreP(j,Lc-La)

      CASE (2,3,4,5)
         Pa = PyramidNodalPBasis(local(1),u,v,w)
         Pb = PyramidNodalPBasis(local(2),u,v,w)
         Pc = PyramidNodalPBasis(local(3),u,v,w)

         La = PyramidTL(local(1),u,v,w)
         Lb = PyramidTL(local(2),u,v,w)
         Lc = PyramidTL(local(3),u,v,w)

         value = Pa*Pb*Pc*LegendreP(i,Lb-La)*LegendreP(j,2*Lc-1)
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidFacePBasis','Unknown face for pyramid')
      END SELECT
    END FUNCTION PyramidFacePBasis

!-----------------------------------------------------------------------------	
!>     Gradient of pyramid face basis at point (u,v,w)
!-----------------------------------------------------------------------------		
    FUNCTION dPyramidFacePBasis(face, i, j, u, v, w, localNumbers )  RESULT(grad)
!-----------------------------------------------------------------------------	
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of pyramids face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)
      ! Variables
      REAL(KIND=dp) :: Pa, Pb, Pc, La, Lb, Lc, legI, legJ, s, swap
      REAL(KIND=dp), DIMENSION(3) :: dPa, dPb, dPc, dLa, dLb, dLc, dLegI,dLegJ,grad
      INTEGER :: local(4)

      ! If local numbering not present use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getPyramidFaceMap(face)
      ELSE
         local(1:4) = localNumbers
      END IF

      ! Calculate value of function by face
      s = w/SQRT(2.0_dp)

      SELECT CASE (face)
      ! Square face
      CASE(1)
         Pa  = PyramidNodalPBasis(local(1),u,v,w)
         Pb  = PyramidNodalPBasis(local(3),u,v,w)
         dPa = dPyramidNodalPBasis(local(1),u,v,w)
         dPb = dPyramidNodalPBasis(local(3),u,v,w)

         La  = PyramidL(local(1),u,v)
         Lb  = PyramidL(local(2),u,v)
         Lc  = PyramidL(local(4),u,v)

         dLa = dPyramidL(local(1),u,v)
         dLb = dPyramidL(local(2),u,v)
         dLc = dPyramidL(local(4),u,v)

         legI = LegendreP(i,Lb-La)
         legJ = LegendreP(j,Lc-La)
         grad = dPa*Pb*legI*legJ + Pa*dPb*legI*legJ +&
              Pa*Pb*dLegendreP(i,Lb-La)*(dLb-dLa)*legJ + &
              Pa*Pb*legI*dLegendreP(j,Lc-La)*(dLc-dLa)
      CASE (2,3,4,5)
      ! Triangle face
         Pa  = PyramidNodalPBasis(local(1),u,v,w)
         Pb  = PyramidNodalPBasis(local(2),u,v,w)
         Pc  = PyramidNodalPBasis(local(3),u,v,w)

         dPa = dPyramidNodalPBasis(local(1),u,v,w)
         dPb = dPyramidNodalPBasis(local(2),u,v,w)
         dPc = dPyramidNodalPBasis(local(3),u,v,w)

         La = PyramidTL(local(1),u,v,w)
         Lb = PyramidTL(local(2),u,v,w)
         Lc = PyramidTL(local(3),u,v,w)

         dLa = dPyramidTL(local(1),u,v,w)
         dLb = dPyramidTL(local(2),u,v,w)
         dLc = dPyramidTL(local(3),u,v,w)

         legI  = LegendreP(i,Lb-La)
         legJ  = LegendreP(j,2*Lc-1)
         dLegI = dLegendreP(i,Lb-La)*(dLb-dLa)
         dLegJ = dLegendreP(j,2*Lc-1)*2*dLc

         grad = dPa*Pb*Pc*LegI*LegJ + Pa*dPb*Pc*LegI*LegJ + &
            Pa*Pb*dPc*LegI*LegJ + Pa*Pb*Pc*dLegI*LegJ + Pa*Pb*Pc*LegI*dLegJ
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidFacePBasis','Unknown face for pyramid')
      END SELECT
    END FUNCTION dPyramidFacePBasis

	
!-----------------------------------------------------------------------------	
!>     Gradient of pyramid face basis at point (u,v,w)
!-----------------------------------------------------------------------------		
    FUNCTION ddPyramidFacePBasis(face, i, j, u, v, w, localNumbers )  RESULT(grad)
!-----------------------------------------------------------------------------	
!
!  ARGUMENTS:
!    INTEGER :: face
!      INPUT: number of pyramids face function to calculate
!        edge = {1,2,..,5}
!
!    INTEGER :: i,j
!      INPUT: index of face function, i,j = {2,3,...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!    INTEGER, OPTIONAL :: localNumber(4)
!      INPUT: local numbering of square or triangle face to define direction 
!        of face basis function. Default numbering is that defined for face in 
!        PElementMaps
! 
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids face function m(i,j) at point (u,v,w), i.e.
!       grad = N_{m(i,j)}^{face}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: face, i, j
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)
      ! Variables
      INTEGER :: p,q
      REAL(KIND=dp) :: Pa, Pb, Pc, La, Lb, Lc, legI, legJ, s, swap, grad(3,3)
      REAL(KIND=dp), DIMENSION(3) :: dPa, dPb, dPc, dLa, dLb, dLc, dLegI,dLegJ
      REAL(KIND=dp), DIMENSION(3,3) :: ddPa, ddPb, ddPc, ddLegI,ddLegJ
      INTEGER :: local(4)

      ! If local numbering not present use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getPyramidFaceMap(face)
      ELSE
         local(1:4) = localNumbers
      END IF

      ! Calculate value of function by face
      s = w/SQRT(2.0_dp)

      SELECT CASE (face)
      ! Square face
      CASE(1)
         Pa  = PyramidNodalPBasis(local(1),u,v,w)
         Pb  = PyramidNodalPBasis(local(3),u,v,w)
         dPa = dPyramidNodalPBasis(local(1),u,v,w)
         dPb = dPyramidNodalPBasis(local(3),u,v,w)

         ddPa = ddPyramidNodalPBasis(local(1),u,v,w)
         ddPb = ddPyramidNodalPBasis(local(3),u,v,w)

         La  = PyramidL(local(1),u,v)
         Lb  = PyramidL(local(2),u,v)
         Lc  = PyramidL(local(4),u,v)

         dLa = dPyramidL(local(1),u,v)
         dLb = dPyramidL(local(2),u,v)
         dLc = dPyramidL(local(4),u,v)

         legI = LegendreP(i,Lb-La)
         legJ = LegendreP(j,Lc-La)

         dLegI = dLegendreP(i,Lb-La)*(dLb-dLa)
         dLegJ = dLegendreP(j,Lc-La)*(dLc-dLa)

         ddLegI = ddLegendreP(i,Lb-La)
         ddLegJ = ddLegendreP(j,Lc-La)

         DO p=1,3
           DO q=p,3
             ddLegI(p,q) = ddLegI(p,q)*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
             ddLegJ(p,q) = ddLegJ(p,q)*(dLc(p)-dLa(p))*(dLc(q)-dLa(q))
           END DO
         END DO

!        grad = dPa*Pb*legI*legJ + Pa*dPb*legI*legJ + Pa*Pb*dLegI*legJ + Pa*Pb*legI*dLegJ

#if 1
         grad = 0
         grad(1,1) = grad(1,1) + ddPa(1,1)*Pb*LegI*LegJ + dPa(1)*dPb(1)*LegI*LegJ + &
                                 dPa(1)*Pb*dLegI(1)*LegJ + dPa(1)*Pb*LegI*dLegJ(1)
         grad(1,1) = grad(1,1) + dPa(1)*dPb(1)*LegI*LegJ + Pa*ddPb(1,1)*LegI*LegJ + &
                                 Pa*dPb(1)*dLegI(1)*LegJ + Pa*dPb(1)*LegI*dLegJ(1)
         grad(1,1) = grad(1,1) + dPa(1)*Pb*dLegI(1)*LegJ + Pa*dPb(1)*dLegI(1)*LegJ + &
                                 Pa*Pb*ddLegI(1,1)*LegJ + Pa*Pb*dLegI(1)*dLegJ(1)
         grad(1,1) = grad(1,1) + dPa(1)*Pb*LegI*dLegJ(1) + Pa*dPb(1)*LegI*dLegJ(1) + &
                                 Pa*Pb*dLegI(1)*dLegJ(1) + Pa*Pb*LegI*ddLegJ(1,1)

         grad(1,2) = grad(1,2) + ddPa(1,2)*Pb*LegI*LegJ + dPa(2)*dPb(1)*LegI*LegJ + &
                                 dPa(2)*Pb*dLegI(1)*LegJ + dPa(2)*Pb*LegI*dLegJ(1)
         grad(1,2) = grad(1,2) + dPa(1)*dPb(2)*LegI*LegJ + Pa*ddPb(1,2)*LegI*LegJ + &
                                 Pa*dPb(2)*dLegI(1)*LegJ + Pa*dPb(2)*LegI*dLegJ(1)
         grad(1,2) = grad(1,2) + dPa(1)*Pb*dLegI(2)*LegJ + Pa*dPb(1)*dLegI(2)*LegJ + &
                                 Pa*Pb*ddLegI(1,2)*LegJ + Pa*Pb*dLegI(2)*dLegJ(1)
         grad(1,2) = grad(1,2) + dPa(1)*Pb*LegI*dLegJ(2) + Pa*dPb(1)*LegI*dLegJ(2) + &
                                 Pa*Pb*dLegI(1)*dLegJ(2) + Pa*Pb*LegI*ddLegJ(1,2)
         grad(2,1) = grad(1,2)

         grad(2,2) = grad(2,2) + ddPa(2,2)*Pb*LegI*LegJ + dPa(2)*dPb(2)*LegI*LegJ + &
                                 dPa(2)*Pb*dLegI(2)*LegJ + dPa(2)*Pb*LegI*dLegJ(2)
         grad(2,2) = grad(2,2) + dPa(2)*dPb(2)*LegI*LegJ + Pa*ddPb(2,2)*LegI*LegJ + &
                                 Pa*dPb(2)*dLegI(2)*LegJ + Pa*dPb(2)*LegI*dLegJ(2)
         grad(2,2) = grad(2,2) + dPa(2)*Pb*dLegI(2)*LegJ + Pa*dPb(2)*dLegI(2)*LegJ + &
                                 Pa*Pb*ddLegI(2,2)*LegJ + Pa*Pb*dLegI(2)*dLegJ(2)
         grad(2,2) = grad(2,2) + dPa(2)*Pb*LegI*dLegJ(2) + Pa*dPb(2)*LegI*dLegJ(2) + &
                                 Pa*Pb*dLegI(2)*dLegJ(2) + Pa*Pb*LegI*ddLegJ(2,2)


         grad(1,3) = grad(1,3) + ddPa(1,3)*Pb*LegI*LegJ + dPa(3)*dPb(1)*LegI*LegJ + &
                                 dPa(3)*Pb*dLegI(1)*LegJ + dPa(3)*Pb*LegI*dLegJ(1)
         grad(1,3) = grad(1,3) + dPa(1)*dPb(3)*LegI*LegJ + Pa*ddPb(1,3)*LegI*LegJ + &
                                 Pa*dPb(3)*dLegI(1)*LegJ + Pa*dPb(3)*LegI*dLegJ(1)
         grad(1,3) = grad(1,3) + dPa(1)*Pb*dLegI(3)*LegJ + Pa*dPb(1)*dLegI(3)*LegJ + &
                                 Pa*Pb*ddLegI(1,3)*LegJ + Pa*Pb*dLegI(3)*dLegJ(1)
         grad(1,3) = grad(1,3) + dPa(1)*Pb*LegI*dLegJ(3) + Pa*dPb(1)*LegI*dLegJ(3) + &
                                 Pa*Pb*dLegI(1)*dLegJ(3) + Pa*Pb*LegI*ddLegJ(1,3)
         grad(3,1) = grad(1,3)

         grad(2,3) = grad(2,3) + ddPa(2,3)*Pb*LegI*LegJ + dPa(3)*dPb(2)*LegI*LegJ + &
                                 dPa(3)*Pb*dLegI(2)*LegJ + dPa(3)*Pb*LegI*dLegJ(2)
         grad(2,3) = grad(2,3) + dPa(2)*dPb(3)*LegI*LegJ + Pa*ddPb(2,3)*LegI*LegJ + &
                                 Pa*dPb(3)*dLegI(2)*LegJ + Pa*dPb(3)*LegI*dLegJ(2)
         grad(2,3) = grad(2,3) + dPa(2)*Pb*dLegI(3)*LegJ + Pa*dPb(2)*dLegI(3)*LegJ + &
                                 Pa*Pb*ddLegI(2,3)*LegJ + Pa*Pb*dLegI(3)*dLegJ(2)
         grad(2,3) = grad(2,3) + dPa(2)*Pb*LegI*dLegJ(3) + Pa*dPb(2)*LegI*dLegJ(3) + &
                                 Pa*Pb*dLegI(2)*dLegJ(3) + Pa*Pb*LegI*ddLegJ(2,3)
         grad(3,2) = grad(2,3)

         grad(3,3) = grad(3,3) + ddPa(3,3)*Pb*LegI*LegJ + dPa(3)*dPb(3)*LegI*LegJ + &
                                 dPa(3)*Pb*dLegI(3)*LegJ + dPa(3)*Pb*LegI*dLegJ(3)
         grad(3,3) = grad(3,3) + dPa(3)*dPb(3)*LegI*LegJ + Pa*ddPb(3,3)*LegI*LegJ + &
                                 Pa*dPb(3)*dLegI(3)*LegJ + Pa*dPb(3)*LegI*dLegJ(3)
         grad(3,3) = grad(3,3) + dPa(3)*Pb*dLegI(3)*LegJ + Pa*dPb(3)*dLegI(3)*LegJ + &
                                 Pa*Pb*ddLegI(3,3)*LegJ + Pa*Pb*dLegI(3)*dLegJ(3)
         grad(3,3) = grad(3,3) + dPa(3)*Pb*LegI*dLegJ(3) + Pa*dPb(3)*LegI*dLegJ(3) + &
                                 Pa*Pb*dLegI(3)*dLegJ(3) + Pa*Pb*LegI*ddLegJ(3,3)
#else
         BLOCK
           REAL(KIND=dp) :: f(4), df(4,3), ddf(4,3,3)

           f(1)=Pa; f(2)=Pb; f(3)=LegI; f(4)=LegJ
           df(1,:)=dPa; df(2,:)=dPb; df(3,:)=dLegI; df(4,:)=dLegJ

           ddf(1,:,:)=ddPa; ddf(2,:,:)=ddPb;
           ddf(3,:,:)=ddLegI; ddf(4,:,:)=ddLegJ

           grad = Product2ndDerivatives(4,f,df,ddf,3,0)
         END BLOCK
         grad(2,1) = grad(1,2)
         grad(3,1) = grad(1,3)
         grad(3,2) = grad(2,3)
#endif

      CASE (2,3,4,5)
      ! Triangle face
         Pa  = PyramidNodalPBasis(local(1),u,v,w)
         Pb  = PyramidNodalPBasis(local(2),u,v,w)
         Pc  = PyramidNodalPBasis(local(3),u,v,w)

         dPa = dPyramidNodalPBasis(local(1),u,v,w)
         dPb = dPyramidNodalPBasis(local(2),u,v,w)
         dPc = dPyramidNodalPBasis(local(3),u,v,w)

         ddPa = ddPyramidNodalPBasis(local(1),u,v,w)
         ddPb = ddPyramidNodalPBasis(local(2),u,v,w)
         ddPc = ddPyramidNodalPBasis(local(3),u,v,w)

         La  = PyramidTL(local(1),u,v,w)
         Lb  = PyramidTL(local(2),u,v,w)
         Lc  = PyramidTL(local(3),u,v,w)

         dLa = dPyramidTL(local(1),u,v,w)
         dLb = dPyramidTL(local(2),u,v,w)
         dLc = dPyramidTL(local(3),u,v,w)

         legI   = LegendreP(i,Lb-La)
         legJ   = LegendreP(j,2*Lc-1)
         dLegI  = dLegendreP(i,Lb-La)*(dLb-dLa)
         dLegJ  = dLegendreP(j,2*Lc-1)*2*dLc
         ddLegI(:,:) = ddLegendreP(i,Lb-La)
         ddLegJ(:,:) = ddLegendreP(j,2*Lc-1)
         DO p=1,3
           DO q=p,3
             ddLegI(p,q) = ddLegI(p,q)*(dLb(p)-dLa(p))*(dLb(q)-dLa(q))
             ddLegJ(p,q) = ddLegJ(p,q)*4*dLc(p)*dLc(q)
           END DO
         END DO

!        grad = dPa*Pb*Pc*LegI*LegJ + Pa*dPb*Pc*LegI*LegJ + &
!          Pa*Pb*dPc*LegI*LegJ + Pa*Pb*Pc*dLegI*LegJ + Pa*Pb*Pc*LegI*dLegJ

         BLOCK
           REAL(KIND=dp) :: f(5), df(5,3), ddf(5,3,3)

           f(1)=Pa; f(2)=Pb; f(3)=Pc; f(4)=LegI; f(5)=LegJ
           df(1,:)=dPa; df(2,:)=dPb; df(3,:)=dPc; df(4,:)=dLegI; df(5,:)=dLegJ

           ddf(1,:,:)=ddPa; ddf(2,:,:)=ddPb; ddf(3,:,:)=ddPc;
           ddf(4,:,:)=ddLegI; ddf(5,:,:)=ddLegJ

           grad = Product2ndDerivatives(5,f,df,ddf,3,0)
         END BLOCK
         grad(2,1) = grad(1,2)
         grad(3,1) = grad(1,3)
         grad(3,2) = grad(2,3)

      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidFacePBasis','Unknown face for pyramid')
      END SELECT
    END FUNCTION ddPyramidFacePBasis

	
	
!------------------------------------------------------------------------------
!>    Pyramid bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION PyramidBubblePBasis(i,j,k,u,v,w) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,0),(0,0,1),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: value
!       value of pyramids bubble function (i,j,k) at point (u,v,w), 
!       i.e. value = N_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: i, j, k 
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: value, s

      ! Calculate value of function
      s = w / SQRT(2._dp)
      value = PyramidNodalPBasis(1,u,v,w)*PyramidNodalPBasis(3,u,v,w)* &
       PyramidNodalPBasis(5,u,v,w)*LegendreP(i,u)*LegendreP(j,v)*LegendreP(k,2*s-1)

!     value = PyramidNodalPBasis(1,u,v,w)*PyramidNodalPBasis(3,u,v,w)* &
!          PyramidNodalPBasis(5,u,v,w)*LegendreP(i,u/(1-w/SQRT(2d0)))* &
!          LegendreP(j,v/(1-w/SQRT(2d0)))*LegendreP(k,w/SQRT(2d0))
    END FUNCTION PyramidBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of pyramid bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION dPyramidBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,0),(0,0,1),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      ! Parameters
      INTEGER, INTENT(IN) :: i, j, k 
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: P1, P3, P5, legI, legJ, legK, s
      REAL(KIND=dp), DIMENSION(3) :: dP1, dP3, dP5, dLegI, dLegJ, &
           dLegK, grad

      s = w/SQRT(2._dp)
      P1  = PyramidNodalPBasis(1,u,v,w)
      P3  = PyramidNodalPBasis(3,u,v,w)
      P5  = PyramidNodalPBasis(5,u,v,w)

      dP1 = dPyramidNodalPBasis(1,u,v,w)
      dP3 = dPyramidNodalPBasis(3,u,v,w)
      dP5 = dPyramidNodalPBasis(5,u,v,w)

      legI = LegendreP(i,u)
      legJ = LegendreP(j,v)
      legK = LegendreP(k,2*s-1)

      dLegI=0; dLegJ=0; dLegK=0
      dLegI(1) = dLegendreP(i,u)
      dLegJ(2) = dLegendreP(j,v)
      dLegK(3) = dLegendreP(k,2*s-1)*2/SQRT(2.0_dp)

      ! Calculate value of gradient
      grad = dP1*P3*P5*legI*legJ*legK + P1*dP3*P5*legI*legJ*legK + &
             P1*P3*dP5*legI*legJ*legK + P1*P3*P5*dLegI*legJ*legK + &
             P1*P3*P5*legI*dLegJ*legK + P1*P3*P5*legI*legJ*dLegK
    END FUNCTION dPyramidBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of pyramid bubble basis at point (u,v,w)
!------------------------------------------------------------------------------
    FUNCTION ddPyramidBubblePBasis(i,j,k,u,v,w) RESULT(grad)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER :: i,j,k
!      INPUT: index of bubble function, (i,j,k) = {(0,0,0),(0,0,1),...}
!
!    REAL(KIND=dp) :: u,v,w
!      INPUT: point at which to evaluate function
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: grad(3)
!       gradient of pyramids bubble function (i,j,k) at point (u,v,w), 
!       i.e. grad = dN_{m(i,j,k)}^{0}(u,v,w)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      ! Parameters
      INTEGER, INTENT(IN) :: i, j, k 
      REAL(KIND=dp), INTENT(IN) :: u,v,w
      ! Variables
      REAL(KIND=dp) :: P1, P3, P5, legI, legJ, legK, s
      REAL(KIND=dp), DIMENSION(3) :: dP1, dP3, dP5, dLegI, dLegJ, dLegK
      REAL(KIND=dp), DIMENSION(3,3) :: ddP1, ddP3, ddP5, ddLegI, ddLegJ, ddLegK, grad

      s = w/SQRT(2._dp)
      P1  = PyramidNodalPBasis(1,u,v,w)
      P3  = PyramidNodalPBasis(3,u,v,w)
      P5  = PyramidNodalPBasis(5,u,v,w)

      dP1 = dPyramidNodalPBasis(1,u,v,w)
      dP3 = dPyramidNodalPBasis(3,u,v,w)
      dP5 = dPyramidNodalPBasis(5,u,v,w)

      ddP1 = ddPyramidNodalPBasis(1,u,v,w)
      ddP3 = ddPyramidNodalPBasis(3,u,v,w)
      ddP5 = ddPyramidNodalPBasis(5,u,v,w)

      legI = LegendreP(i,u)
      legJ = LegendreP(j,v)
      legK = LegendreP(k,2*s-1)

      dLegI=0; dLegJ=0; dLegK=0
      dLegI(1) = dLegendreP(i,u)
      dLegJ(2) = dLegendreP(j,v)
      dLegK(3) = dLegendreP(k,2*s-1)*2/SQRT(2.0_dp)

      ddLegI=0; ddLegJ=0; ddLegK=0
      ddLegI(1,1) = ddLegendreP(i,u)
      ddLegJ(2,2) = ddLegendreP(j,v)
      ddLegK(3,3) = ddLegendreP(k,2*s-1)*2

      BLOCK
        REAL(KIND=dp) :: f(6), df(6,3), ddf(6,3,3)

        f(1)=P1; f(2)=P3; f(3)=P5; f(4)=LegI; f(5)=LegJ; f(6)=LegK
        df(1,:)=dP1; df(2,:)=dP3; df(3,:)=dP5;
        df(4,:)=dLegI; df(5,:)=dLegJ; df(6,:)=dLegK
        ddf(1,:,:)=ddP1; ddf(2,:,:)=ddP3; ddf(3,:,:)=ddP5;
        ddf(4,:,:)=ddLegI; ddf(5,:,:)=ddLegJ; ddf(6,:,:)=ddLegK

        grad = Product2ndDerivatives(6,f,df,ddf,3,0)
      END BLOCK
      grad(2,1) = grad(1,2)
      grad(3,1) = grad(1,3)
      grad(3,2) = grad(2,3)
    END FUNCTION ddPyramidBubblePBasis

    ! Define affine coordinates for pyramid square face
    PURE FUNCTION PyramidL(which, u, v) RESULT(value)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v
      ! Variables
      REAL(KIND=dp) :: value
      
      SELECT CASE (which)
      CASE (1)
         value = ((1-u)+(1-v))/2 
      CASE (2)
         value = ((1+u)+(1-v))/2
      CASE (3)
         value = ((1+u)+(1+v))/2
      CASE (4)
         value = ((1-u)+(1+v))/2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidL','Unknown affine coordinate for square face')
#endif
      END SELECT
    END FUNCTION PyramidL

    PURE FUNCTION dPyramidL(which, u, v) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v
      ! Variables
      REAL(KIND=dp) :: grad(3)

      SELECT CASE (which)
      CASE (1)
         grad = [-1d0/2,-1d0/2,0d0 ]
      CASE (2)
         grad = [ 1d0/2,-1d0/2,0d0 ]
      CASE (3)
         grad = [ 1d0/2, 1d0/2,0d0 ]
      CASE (4)
         grad = [-1d0/2, 1d0/2,0d0 ]
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidL','Unknown affine coordinate for square face')
#endif
      END SELECT
    END FUNCTION dPyramidL


    PURE FUNCTION PyramidTL(which, u, v, w) RESULT(value)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v,w 
      ! Variables
      REAL(KIND=dp) :: value,s
      
      value = 0
      s = w/SQRT(2.0_dp)
      SELECT CASE(which)
      CASE (1)
         value = (2-u-v-s)/2
      CASE (2)
         value = (2+u-v-s)/2
      CASE (3)
         value = (2+u+v-s)/2
      CASE (4)
         value = (2-u+v-s)/2
      CASE (5)
         value = s
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidTL','Unknown function L for brick')
#endif
      END SELECT
    END FUNCTION PyramidTL

    PURE FUNCTION dPyramidTL(which, u, v, w) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v,w 
      ! Variables
      REAL(KIND=dp) :: grad(3),s
      
      s = w/SQRT(2.0_dp)
      grad = 0
      SELECT CASE(which)
      CASE (1)
         grad(1) = -1
         grad(2) = -1
         grad(3) = -1
      CASE (2)
         grad(1) =  1
         grad(2) = -1
         grad(3) = -1
      CASE (3)
         grad(1) =  1
         grad(2) =  1
         grad(3) = -1
      CASE (4)
         grad(1) = -1
         grad(2) =  1
         grad(3) = -1
      CASE (5)
         grad(3) =  2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidTL','Unknown function L for brick')
#endif
      END SELECT
      grad = grad/2
      grad(3) = grad(3)/SQRT(2._dp)
    END FUNCTION dPyramidTL


!------------------------------------------------------------------------------
!>    Phi function value at point x. Phi is defined as (Szabo & Babuska: Finite
!>    Element Analysis, p.38). 
!------------------------------------------------------------------------------
    PURE FUNCTION Phi(i,x) RESULT(value)
!------------------------------------------------------------------------------
!
!  DESCRIPTION:
!    
!    Phi(i,x)=SQRT(1/(2*(2*i-1)))(P(i,x)-P(i-2,x)), i=2,3,...
!
!    where P(i,x) are legendre polynomials,
!
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: i
!      INPUT: parameter of phi
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate phi
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of phi function i at point x i.e
!       value = Phi(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: x
      ! Return value
      REAL (KIND=dp) :: value

#ifdef DEBUG_PBASIS
      IF (i < 2) THEN
         CALL Fatal('PElementBase::Phi','Phi(i,x) not defined for i<2')
      END IF
#endif

      ! Check if value is available by direct calculation
      IF (i <= 20) THEN
         ! Get value from the varphi function
         value = varPhi(i,x)*(1-x**2)/4
      ELSE 
         value = SQRT(1d0/(2*(2*i-1)))*(LegendreP(i,x)-LegendreP(i-2,x))
      END IF
    END FUNCTION Phi
   
 
!------------------------------------------------------------------------------
!>    Derivative of phi function value at point x.   
!>    Phi,(i,x)=SQRT(1/(2*(2*i-1)))(P,(i,x)-P,(i-2,x)), i=2,3,... 
!>    where P,(i,x) are derivatives of legendre polynomials.
!------------------------------------------------------------------------------
    PURE FUNCTION dPhi(i,x) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: i
!      INPUT: parameter of phi
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate phi
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of derivated phi function i at point x i.e
!       value = Phi,(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: x
      ! Return value
      REAL (KIND=dp) :: value

#ifdef DEBUG_PBASIS
      IF (i < 2) THEN
         CALL Fatal('PElementBase::dPhi','dPhi(i,x) not defined for i<2')
      END IF
#endif

      ! 20 first derivatives of phi functions are precalculated
      ! They are all generated with Maple 
      SELECT CASE (i)
      CASE (2)
         value = sqrt(0.6D1) * x / 0.2D1      
      CASE (3)
         value = sqrt(0.10D2) * dble(3 * x ** 2 - 1) / 0.4D1
      CASE (4)
         value = sqrt(0.14D2) * dble(x) * dble(5 * x ** 2 - 3) / 0.4D1
      CASE (5)
         value = 0.3D1 / 0.16D2 * sqrt(0.2D1) * dble(35 * x ** 4 - 30 * x ** 2 + 3)
      CASE (6)
         value = sqrt(0.22D2) * dble(x) * dble(63 * x ** 4 - 70 * x ** 2 +& 
              15) / 0.16D2
      CASE (7)
         value = sqrt(0.26D2) * dble(231 * x ** 6 - 315 * x ** 4 + 105 * &
              x ** 2 - 5) / 0.32D2
      CASE (8)    
         value = sqrt(0.30D2) * dble(x) * dble(429 * x ** 6 - 693 * x ** 4 +& 
              315 * x ** 2 - 35) / 0.32D2
      CASE (9)
         value = sqrt(0.34D2) * dble(6435 * x ** 8 - 12012 * x ** 6 + 6930 * &
              x ** 4 - 1260 * x ** 2 + 35) / 0.256D3
      CASE (10)
         value = sqrt(0.38D2) * dble(x) * dble(12155 * x ** 8 - 25740 * x ** 6 &
              + 18018 * x ** 4 - 4620 * x ** 2 + 315) / 0.256D3
      CASE (11)   
         value = sqrt(0.42D2) * dble(46189 * x ** 10 - 109395 * x ** 8 + 90090 &
              * x ** 6 - 30030 * x ** 4 + 3465 * x ** 2 - 63) / 0.512D3
      CASE (12)
         value = sqrt(0.46D2) * dble(x) * dble(88179 * x ** 10 - 230945 * x ** 8 + & 
              218790 * x ** 6 - 90090 * x ** 4 + 15015 * x ** 2 - 693) / 0.512D3
      CASE (13)
         value = 0.5D1 / 0.2048D4 * sqrt(0.2D1) * dble(676039 * x ** 12 - 1939938 * &
              x ** 10 + 2078505 * x ** 8 - 1021020 * x ** 6 + 225225 * x ** 4 - 18018 * &
              x ** 2 + 231)
      CASE (14)
         value = 0.3D1 / 0.2048D4 * sqrt(0.6D1) * dble(x) * dble(1300075 * x ** 12 &
              - 4056234D0 * x ** 10 + 4849845D0 * x ** 8 - 2771340 * x ** 6 + 765765 * x ** 4 &
              - 90090D0 * x ** 2 + 3003)
      CASE (15)
         value = sqrt(0.58D2) * dble(5014575D0 * x ** 14 - 16900975D0 * x ** 12 &
              + 22309287D0 * x ** 10 - 14549535D0 * x ** 8 + 4849845D0 * x ** 6 - 765765 * &
              x ** 4 + 45045 * x ** 2 - 429) / 0.4096D4
      CASE (16)
         value = sqrt(0.62D2) * dble(x) * dble(9694845D0 * x ** 14 - 35102025D0 &
              * x ** 12 + 50702925D0 * x ** 10 - 37182145D0 * x ** 8 + 14549535D0 * x ** 6 - & 
              2909907D0 * x ** 4 + 255255 * x ** 2 - 6435) / 0.4096D4
      CASE (17)
         value = sqrt(0.66D2) * dble(300540195D0 * x ** 16 - 1163381400D0 * x ** 14 &
              + 1825305300D0 * x ** 12 - 1487285800D0 * x ** 10 + 669278610D0 * x** 8 - & 
              162954792D0 * x ** 6 + 19399380D0 * x ** 4 - 875160 * x ** 2 + 6435) / 0.65536D5
      CASE (18)
         value = sqrt(0.70D2) * dble(x) * dble(583401555D0 * x ** 16 - 2404321560D0 * &
              x ** 14 + 4071834900D0 * x ** 12 - 3650610600D0 * x ** 10 + 1859107250D0 * & 
              x ** 8 - 535422888D0 * x ** 6 + 81477396D0 * x ** 4 - 5542680 * x ** 2 + & 
              109395) / 0.65536D5
      CASE (19)      
         value = sqrt(0.74D2) * dble(2268783825D0 * x ** 18 - 9917826435D0 * x ** 16 + &
              18032411700D0 * x ** 14 - 17644617900D0 * x ** 12 + 10039179150D0 * x ** 10 -& 
              3346393050D0 * x ** 8 + 624660036D0 * x ** 6 - 58198140D0 * x ** 4 + 2078505D0 * &
              x ** 2 - 12155) / 0.131072D6
      CASE (20)    
         value = sqrt(0.78D2) * dble(x) * dble(4418157975D0 * x ** 18 - 20419054425D0 & 
              * x ** 16 + 39671305740D0 * x ** 14 - 42075627300D0 * x ** 12 + &
              26466926850D0 * x ** 10 - 10039179150D0 * x ** 8 + 2230928700D0 * x ** 6 - &
              267711444D0 * x ** 4 + 14549535D0 * x ** 2 - 230945) / 0.131072D6
         ! If no precalculated value available generate value of function
#ifdef DEBUG_PBASIS
      CASE DEFAULT 
         value = SQRT(1d0/(2*(2*i-1)))*(dLegendreP(i,x)-dLegendreP(i-2,x))
#endif
      END SELECT
    END FUNCTION dPhi

!------------------------------------------------------------------------------
!>    2nd derivative of phi function value at point x.   
!>    Phi,(i,x)=SQRT(1/(2*(2*i-1)))(P,(i,x)-P,(i-2,x)), i=2,3,... 
!>    where P,(i,x) are derivatives of legendre polynomials.
!------------------------------------------------------------------------------
    PURE FUNCTION ddPhi(i,x) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: i
!      INPUT: parameter of phi
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate phi
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of derivated phi function i at point x i.e
!       value = Phi,(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: x
      ! Return value
      REAL (KIND=dp) :: value

#ifdef DEBUG_PBASIS
      IF (i < 2) THEN
         CALL Fatal('PElementBase::dPhi','dPhi(i,x) not defined for i<2')
      END IF
#endif

      ! 20 first derivatives of phi functions are precalculated
      ! They are all generated with Maple 
      SELECT CASE (i)
      CASE (2)
         value = sqrt(0.6D1) / 0.2D1      
      CASE (3)
         value = sqrt(0.10D2) * 3 * 2*x / 0.4D1
      CASE (4)
         value = sqrt(0.14D2) * (5 * 3*x ** 2 - 3) / 0.4D1
      CASE (5)
         value = 0.3D1 / 0.16D2 * sqrt(0.2D1) * (35 * 4*x ** 3 - 30 * 2*x)
      CASE (6)
         value = sqrt(0.22D2) * (63 * 5*x ** 4 - 70 * 3*x ** 2 + 15) / 0.16D2
      CASE (7)
         value = sqrt(0.26D2) * (6*231*x**5 - 4*315*x**3 + 2*105*x) / 0.32D2
      CASE (8)    
         value = sqrt(0.30D2) * (7*429*x**6 - 5*693*x**4 + 3*315*x**2 - 35) / 0.32D2
      CASE (9)
         value = sqrt(0.34D2) * (8*6435*x**7 - 6*12012*x**5  + 4*6930*x**3 - 2*1260*x) / 0.256D3
      CASE (10)
         value = sqrt(0.38D2) * (9*12155*x**8 - 7*25740*x**6 &
              + 5*18018*x**4 - 3*4620*x**2 + 315) / 0.256D3
      CASE (11)   
         value = sqrt(0.42D2) * dble(10*46189*x**9 - 8*109395*x**7 + 6*90090 &
              *x**5 - 4*30030*x**3 + 2*3465*x) / 0.512D3
      CASE (12)
         value = sqrt(0.46D2) * (11*88179*x**10 - 9*230945*x**8 + &
              7*218790*x**6 - 5*90090*x**4 + 3*15015*x**2 - 693) / 0.512D3
      CASE (13)
         value = 0.5D1 / 0.2048D4 * sqrt(0.2D1) * (12*676039*x**11 - 10*1939938 * &
              x**9 + 8*2078505*x**7 - 6*1021020*x**5 + 4*225225*x**3 - 2*18018*x)
      CASE (14)
         value = 0.3D1 / 0.2048D4 * sqrt(0.6D1) * (13*1300075 * x ** 12 &
              - 11*4056234D0 * x ** 10 + 9*4849845D0 * x ** 8 - 7*2771340 * x ** 6 + 5*765765 * x ** 4 &
              - 3*90090D0 * x ** 2 + 3003)
      CASE (15)
         value = sqrt(0.58D2) * (14*5014575D0*x**13 - 12*16900975D0*x**11 &
              + 10*22309287D0*x**9 - 8*14549535D0*x**7 + 6*4849845D0*x**5 - 4*765765 * &
              x**3 + 2*45045*x) / 0.4096D4
      CASE (16)
         value = sqrt(0.62D2) * (15*9694845D0 * x ** 14 - 13*35102025D0 &
              * x ** 12 + 11*50702925D0 * x ** 10 - 9*37182145D0 * x ** 8 + 7*14549535D0 * x ** 6 - & 
              5*2909907D0 * x ** 4 + 3*255255 * x ** 2 - 6435) / 0.4096D4
      CASE (17)
         value = sqrt(0.66D2) * (16*300540195D0*x**15 - 14*1163381400D0*x**13 &
              + 12*1825305300D0*x**11 - 10*1487285800D0*x**9 + 8*669278610D0*x**7 - & 
              6*162954792D0*x**5 + 4*19399380D0*x**3 - 2*875160*x) / 0.65536D5
      CASE (18)
         value = sqrt(0.70D2) * (17*583401555D0 * x ** 16 - 15*2404321560D0 * &
              x ** 14 + 13*4071834900D0 * x ** 12 - 11*3650610600D0 * x ** 10 + 9*1859107250D0 * & 
              x ** 8 - 7*535422888D0 * x ** 6 + 5*81477396D0 * x ** 4 - 3*5542680 * x ** 2 + & 
              109395) / 0.65536D5
      CASE (19)      
         value = sqrt(0.74D2) * (18*2268783825D0*x**17 - 16*9917826435D0*x**15 + &
              14*18032411700D0*x**13 - 12*17644617900D0*x** 11 + 10*10039179150D0*x**9 -& 
              8*3346393050D0*x**7 + 6*624660036D0*x**5 - 4*58198140D0*x**3 + 2*2078505D0 * &
              x) / 0.131072D6
      CASE (20)    
         value = sqrt(0.78D2) * (19*4418157975D0 * x ** 18 - 17*20419054425D0 & 
              * x ** 16 + 15*39671305740D0 * x ** 14 - 13*42075627300D0 * x ** 12 + &
              11*26466926850D0 * x ** 10 - 9*10039179150D0 * x ** 8 + 7*2230928700D0 * x ** 6 - &
              5*267711444D0 * x ** 4 + 3*14549535D0 * x ** 2 - 230945) / 0.131072D6
         ! If no precalculated value available generate value of function
      CASE DEFAULT 
#ifdef DEBUG_PBASIS
         PRINT*,'Legendre phi: ', i
         STOP 'no ddph > 20'
!        value = SQRT(1d0/(2*(2*i-1)))*(dLegendreP(i,x)-dLegendreP(i-2,x))
#endif
      END SELECT
    END FUNCTION ddPhi



    
!------------------------------------------------------------------------------
!>    varPhi function value at point x. Phi is defined as (Szabo & Babuska: Finite
!>    Element Analysis, p. 103).     
!>    Phi(i,x)=1/4*(1-x^2)*varPhi(i,x)
!------------------------------------------------------------------------------
    PURE FUNCTION varPhi(i,x) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: i
!      INPUT: parameter of varPhi
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate varPhi
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of varPhi function i at point x i.e
!       value = varPhi(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: x
      REAL (KIND=dp), PARAMETER :: dx = 1E-11
      REAL (KIND=dp) :: value
      
      ! 20 first varphi functions are precalculated
      ! These are all generated with Maple
      SELECT CASE(i)
#ifdef DEBUG_PBASIS
      CASE (:1)
         CALL Fatal('PElementBase::varPhi','varPhi not defined for i<2')
#endif
      CASE (2)
         value = -SQRT(0.6D1)
      CASE (3)
         value = -x * SQRT(0.10D2)
      CASE (4)
         value = -DBLE(5 * x ** 2 - 1) * SQRT(0.14D2) / 0.4D1
      CASE (5)
         value = -0.3D1 / 0.4D1 * DBLE(7 * x ** 2 - 3) * DBLE(x) * SQRT(0.2D1)
      CASE (6)
         value = -DBLE(21 * x ** 4 - 14 * x ** 2 + 1) * SQRT(0.22D2) / 0.8D1
      CASE (7)
         value = -DBLE(33 * x ** 4 - 30 * x ** 2 + 5) * DBLE(x) * SQRT(0.26D2) / 0.8D1
      CASE (8)
         value = -DBLE(429 * x ** 6 - 495 * x ** 4 + 135 * x ** 2 - 5) & 
              * SQRT(0.30D2) / 0.64D2
      CASE (9)
         value = -DBLE(715 * x ** 6 - 1001 * x ** 4 + 385 * x ** 2 - 35) & 
              * DBLE(x) * SQRT(0.34D2) / 0.64D2
      CASE (10)
         value = -DBLE(2431 * x ** 8 - 4004 * x ** 6 + 2002 * x ** 4 - 308 & 
              * x ** 2 + 7) * SQRT(0.38D2) / 0.128D3
      CASE (11)
         value = -DBLE(4199 * x ** 8 - 7956 * x ** 6 + 4914 * x ** 4 - 1092 &
              * x ** 2 + 63) * DBLE(x) * SQRT(0.42D2) / 0.128D3
      CASE (12)
         value = -DBLE(29393 * x ** 10 - 62985 * x ** 8 + 46410 * x ** 6 - &
              13650 * x ** 4 + 1365 * x ** 2 - 21) * SQRT(0.46D2) / 0.512D3
      CASE (13)      
         value = -0.5D1 / 0.512D3 * DBLE(52003D0 * x ** 10 - 124355D0 * x ** 8 &
              + 106590D0 * x ** 6 - 39270D0 * x ** 4 + 5775 * x ** 2 - 231) * & 
              DBLE(x) * SQRT(0.2D1)
      CASE (14)
         value = -0.3D1 / 0.1024D4 * DBLE(185725D0 * x ** 12 - 490314D0 * x ** &
              10 + 479655D0 * x ** 8 - 213180D0 * x ** 6 + 42075 * x ** 4 - 2970 & 
              * x** 2 + 33) * SQRT(0.6D1)
      CASE (15)
         value = -DBLE(334305D0 * x ** 12 - 965770D0 * x ** 10 + 1062347D0 * x ** &
              8 - 554268D0 * x ** 6 + 138567D0 * x ** 4 - 14586 * x ** 2 + 429) * &
              DBLE(x) * SQRT(0.58D2) / 0.1024D4
      CASE (16)
         value = -DBLE(9694845D0 * x ** 14 - 30421755D0 * x ** 12 + 37182145D0 * &
              x ** 10 - 22309287D0 * x ** 8 + 6789783D0 * x ** 6 - 969969D0 * x ** 4 + &
              51051D0 * x ** 2 - 429) * SQRT(0.62D2) / 0.16384D5
      CASE (17)
         value = -DBLE(17678835D0 * x ** 14 - 59879925D0 * x ** 12 + 80528175D0 * & 
              x ** 10 - 54679625D0 * x ** 8 + 19684665D0 * x ** 6 - 3594591D0 * x ** 4 &
              + 285285D0 * x ** 2 - 6435) * dble(x) * SQRT(0.66D2) / 0.16384D5
      CASE (18)
         value = -DBLE(64822395D0 * x ** 16 - 235717800D0 * x ** 14 + 345972900D0 * & 
              x ** 12 - 262462200D0 * x ** 10 + 109359250D0 * x ** 8 - 24496472 * &
              x ** 6 + 2662660D0 * x ** 4 - 108680 * x ** 2 + 715) * SQRT(0.70D2) / 0.32768D5
      CASE (19)
         value = -DBLE(119409675D0 * x ** 16 - 463991880D0 * x ** 14 + 738168900D0 & 
              * x ** 12 - 619109400D0 * x ** 10 + 293543250D0 * x ** 8 - 78278200D0 * &
              x ** 6 + 10958948D0 * x ** 4 - 680680D0 * x ** 2 + 12155) * DBLE(x) * & 
              SQRT(0.74D2) / 0.32768D5
      CASE (20)
         value = -DBLE(883631595D0 * x ** 18 - 3653936055D0 * x ** 16 + 6263890380D0 &
              * x ** 14 - 5757717420D0 * x ** 12 + 3064591530D0 * x ** 10 - 951080130D0 & 
              * x ** 8 + 164384220D0 * x ** 6 - 14090076D0 * x ** 4 + 459459 * &
              x ** 2 - 2431) * SQRT(0.78D2) / 0.131072D6
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         IF (x==1 .OR. x==-1) THEN
            ! TEMP SOLUTION!
            ! Try to interpolate value of function
            value = ((4*Phi(i,(x-dx))/(1-(x-dx)**2))+(4*Phi(i,(x+dx))/(1-(x+dx)**2)))/2
         ELSE 
            value = 4*Phi(i,x)/(1-x**2)
         END IF
#endif
      END SELECT
      
    END FUNCTION varPhi


!------------------------------------------------------------------------------
!>    Derivative of varPhi function at point x.
!------------------------------------------------------------------------------
    PURE FUNCTION dVarPhi(i,x) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: i
!      INPUT: parameter of varPhi
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate varPhi
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of derivative of varPhi function i at point x i.e
!       value = dVarPhi(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: x
      REAL (KIND=dp), PARAMETER :: dx = 0.001 !1E-10
      REAL (KIND=dp) :: value, vp, vm

      ! 
      SELECT CASE(i)
#ifdef DEBUG_PBASIS
      CASE (:1)
         CALL Fatal('PElementBase::dVarPhi','dVarPhi not defined for i<2')
#endif
      CASE (2)
         value = 0      
      CASE (3)
         value = -SQRT(0.10D2)
      CASE (4)
         value = -0.5D1 / 0.2D1 * x * SQRT(0.14D2)
      CASE (5)
         value = -0.63D2 / 0.4D1 * x ** 2 * SQRT(0.2D1) + 0.9D1 / 0.4D1 * SQRT(0.2D1)
      CASE (6)
         value = -0.7D1 / 0.2D1 * DBLE(x) * DBLE(3 * x ** 2 - 1) * SQRT(0.22D2)
      CASE (7)
         value = -0.165D3 / 0.8D1 * x ** 4 * SQRT(0.26D2) + 0.45D2 / 0.4D1 * & 
              x ** 2 * SQRT(0.26D2) - 0.5D1 / 0.8D1 * SQRT(0.26D2)
      CASE (8)
         value = -0.9D1 / 0.32D2 * DBLE(x) * DBLE(143 * x ** 4 - 110 * x **2 + 15) & 
              * SQRT(0.30D2)
      CASE (9)
         value = -0.5005D4 / 0.64D2 * x ** 6 * SQRT(0.34D2) + 0.5005D4 / 0.64D2 * &
              x ** 4 * SQRT(0.34D2) - 0.1155D4 / 0.64D2 * x ** 2 * SQRT(0.34D2) + 0.35D2 & 
              / 0.64D2 * SQRT(0.34D2)
      CASE (10)
         value = -0.11D2 / 0.16D2 * DBLE(x) * DBLE(221 * x ** 6 - 273 * x ** 4 & 
              + 91 * x ** 2 - 7) * SQRT(0.38D2)
      CASE (11)
         value = -0.37791D5 / 0.128D3 * x ** 8 * SQRT(0.42D2) + 0.13923D5 / &
              0.32D2 * x ** 6 * SQRT(0.42D2) - 0.12285D5 / 0.64D2 * x ** 4 * SQRT(0.42D2) &
              + 0.819D3 / 0.32D2 * x ** 2 * SQRT(0.42D2) - 0.63D2 / 0.128D3 * SQRT(0.42D2)
      CASE (12)
         value = -0.65D2 / 0.256D3 * DBLE(x) * DBLE(2261 * x ** 8 - 3876 * &
              x ** 6 + 2142 * x ** 4 - 420 * x ** 2 + 21) * SQRT(0.46D2)
      CASE (13)
         value = -0.2860165D7 / 0.512D3 * x ** 10 * SQRT(0.2D1) + 0.5595975D7 &
              / 0.512D3 * x ** 8 * SQRT(0.2D1) - 0.1865325D7 / 0.256D3 * x **6 * &
              SQRT(0.2D1) + 0.490875D6 / 0.256D3 * x ** 4 * SQRT(0.2D1) - 0.86625D5 &
              / 0.512D3 * x ** 2 * SQRT(0.2D1) + 0.1155D4 / 0.512D3 * SQRT(0.2D1)
      CASE (14)
         value = -0.45D2 / 0.256D3 * DBLE(x) * DBLE(37145 * x ** 10 - 81719* &
              x ** 8 + 63954 * x ** 6 - 21318 * x ** 4 + 2805 * x ** 2 - 99) * SQRT(0.6D1)
      CASE (15)
         value = -0.4345965D7 / 0.1024D4 * x ** 12 * SQRT(0.58D2) + 0.5311735D7 & 
              / 0.512D3 * x ** 10 * SQRT(0.58D2) - 0.9561123D7 / 0.1024D4 * x ** 8 * &
              SQRT(0.58D2) + 0.969969D6 / 0.256D3 * x ** 6 * SQRT(0.58D2) - 0.692835D6 &
              / 0.1024D4 * x ** 4 * SQRT(0.58D2) + 0.21879D5 / 0.512D3 * x ** 2 * & 
              SQRT(0.58D2) - 0.429D3 / 0.1024D4 * SQRT(0.58D2)
      CASE (16)
         value = -0.119D3 / 0.8192D4 * DBLE(x) * DBLE(570285D0 * x ** 12 - 1533870D0 * &
              x ** 10 + 1562275D0 * x ** 8 - 749892D0 * x ** 6 + 171171 * x ** 4 - 16302 * &
              x ** 2 + 429) * SQRT(0.62D2)
      CASE (17)
         value = -0.265182525D9 / 0.16384D5 * x ** 14 * SQRT(0.66D2) + 0.778439025D9 & 
              / 0.16384D5 * x ** 12 * SQRT(0.66D2) - 0.885809925D9 / 0.16384D5 * & 
              x ** 10 * SQRT(0.66D2) + 0.492116625D9 / 0.16384D5 * x ** 8 * SQRT(0.66D2) - &
              0.137792655D9 / 0.16384D5 * x ** 6 * SQRT(0.66D2) + 0.17972955D8 / 0.16384D5 * &
              x ** 4 * SQRT(0.66D2) - 0.855855D6 / 0.16384D5 * x ** 2 * SQRT(0.66D2) + &
              0.6435D4 / 0.16384D5 * SQRT(0.66D2)
      CASE (18)
         value = -0.19D2 / 0.2048D4 * DBLE(x) * DBLE(3411705D0 * x ** 14 - 10855425D0 * &
              x ** 12 + 13656825D0 * x ** 10 - 8633625D0 * x ** 8 + 2877875D0 * x ** 6 - &
              483483D0 * x ** 4 + 35035 * x ** 2 - 715) * SQRT(0.70D2)
      CASE (19)
         value = -0.2029964475D10 / 0.32768D5 * x ** 16 * SQRT(0.74D2) + 0.869984775D9 &
              / 0.4096D4 * x ** 14 * SQRT(0.74D2) - 0.2399048925D10 &
              / 0.8192D4 * x ** 12 * SQRT(0.74D2) + 0.851275425D9 / 0.4096D4 * x** 10 * &
              SQRT(0.74D2) - 0.1320944625D10 / 0.16384D5 * x ** 8 * SQRT(0.74D2) + &
              0.68493425D8 / 0.4096D4 * x ** 6 * SQRT(0.74D2) - 0.13698685D8 / 0.8192D4 * &
              x ** 4 * SQRT(0.74D2) + 0.255255D6 / 0.4096D4 * x ** 2 * SQRT(0.74D2) - &
              0.12155D5 / 0.32768D5 * SQRT(0.74D2)
      CASE (20)
         value = -0.63D2 / 0.65536D5 * DBLE(x) * DBLE(126233085D0 * x ** 16 - &
              463991880D0 * x ** 14 + 695987820D0 * x ** 12 - 548354040D0 * x ** 10 + &
              243221550D0 * x ** 8 - 60386040D0 * x ** 6 + 7827820D0 * x ** 4 - 447304D0 * &
              x ** 2 + 7293) * SQRT(0.78D2)
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         IF (x==1 .OR. x==-1) THEN
            ! TEMP SOLUTION
            ! Try to interpolate value of function
            vp = 4*((1-(x+dx)**2)*dPhi(i,(x+dx))+2*(x+dx)*Phi(i,(x+dx)))/(1-(x+dx)**2)**2
            vm = 4*((1-(x-dx)**2)*dPhi(i,(x-dx))+2*(x-dx)*Phi(i,(x-dx)))/(1-(x-dx)**2)**2
            value = (vp+vm)/2
            ! WRITE (*,*) i,x,vp,vm,value
         ELSE
            value = 4*((1-x**2)*dPhi(i,x)+2*x*Phi(i,x))/(1-x**2)**2
         END IF
#endif
      END SELECT
    END FUNCTION dVarPhi

!------------------------------------------------------------------------------
!>    2nd derivatives of varPhi function at point x.
!------------------------------------------------------------------------------
    PURE FUNCTION ddVarPhi(i,x) RESULT(value)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: i
!      INPUT: parameter of varPhi
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate varPhi
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of derivative of varPhi function i at point x i.e
!       value = dVarPhi(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      REAL (KIND=dp), INTENT(IN) :: x
      REAL (KIND=dp) :: value, vp, vm
      REAL (KIND=dp), PARAMETER :: dx = 0.001 !1E-10

      ! 
      SELECT CASE(i)
#ifdef DEBUG_PBASIS
      CASE (:1)
         CALL Fatal('PElementBase::dVarPhi','dVarPhi not defined for i<2')
#endif
      CASE (2)
         value = 0      
      CASE (3)
         value = 0
      CASE (4)
         value = -0.5D1 / 0.2D1 * SQRT(0.14D2)
      CASE (5)
         value = -0.63D2 / 0.4D1 * 2*x * SQRT(0.2D1)
      CASE (6)
         value = -0.7D1 / 0.2D1 * (3 * 3*x ** 2 - 1) * SQRT(0.22D2)
      CASE (7)
         value = -0.165D3 / 0.8D1 * 4*x**3 * SQRT(0.26D2) + 0.45D2 / 0.4D1 * & 
              2*x * SQRT(0.26D2)
      CASE (8)
         value = -0.9D1 / 0.32D2 * (143 * 5*x**4 - 3*110*x**2 + 15) * SQRT(0.30D2)
      CASE (9)
         value = -0.5005D4 / 0.64D2 * 6*x**5 * SQRT(0.34D2) + 0.5005D4 / 0.64D2 * &
              4*x**3 * SQRT(0.34D2) - 0.1155D4 / 0.64D2 * 2*x*SQRT(0.34D2)
      CASE (10)
         value = -0.11D2 / 0.16D2 * (221 * 7*x* 6 - 273*5*x**4 & 
              + 91 * 3*x**2 - 7) * SQRT(0.38D2)
      CASE (11)
         value = -0.37791D5 / 0.128D3 * 8*x**7 * SQRT(0.42D2) + 0.13923D5 / &
              0.32D2 * 6*x**5 * SQRT(0.42D2) - 0.12285D5 / 0.64D2 * 4*x**3 * SQRT(0.42D2) &
              + 0.819D3 / 0.32D2 * 2*x * SQRT(0.42D2)
      CASE (12)
         value = -0.65D2 / 0.256D3 * (2261 * 9*x**8 - 3876 * &
              7*x**6 + 2142*5*x**4 - 420 * 3*x**2 + 21) * SQRT(0.46D2)
      CASE (13)
         value = -0.2860165D7 / 0.512D3 * 10*x**9 * SQRT(0.2D1) + 0.5595975D7 &
              / 0.512D3 * 8*x**7 * SQRT(0.2D1) - 0.1865325D7 / 0.256D3 * 6*x**5 * &
              SQRT(0.2D1) + 0.490875D6 / 0.256D3 * 4*x**3 * SQRT(0.2D1) - 0.86625D5 &
              / 0.512D3 * 2*x * SQRT(0.2D1)
      CASE (14)
         value = -0.45D2 / 0.256D3 * (37145 * 11*x**10 - 81719* &
              9*x ** 8 + 63954 * 7*x ** 6 - 21318 * 5*x ** 4 + 2805 * 3*x ** 2 - 99) * SQRT(0.6D1)
      CASE (15)
         value = -0.4345965D7 / 0.1024D4 * 12*x ** 11 * SQRT(0.58D2) + 0.5311735D7 & 
              / 0.512D3 * 10*x ** 9 * SQRT(0.58D2) - 0.9561123D7 / 0.1024D4 * 8*x ** 7 * &
              SQRT(0.58D2) + 0.969969D6 / 0.256D3 * 6*x ** 5 * SQRT(0.58D2) - 0.692835D6 &
              / 0.1024D4 * 4*x ** 3 * SQRT(0.58D2) + 0.21879D5 / 0.512D3 * 2*x * & 
              SQRT(0.58D2)
      CASE (16)
         value = -0.119D3 / 0.8192D4 * (570285D0 * 13*x ** 12 - 1533870D0 * &
              11*x ** 10 + 1562275D0 * 9*x ** 8 - 749892D0 * 7*x ** 6 + 171171 * 5*x ** 4 - 16302 * &
              3*x ** 2 + 429) * SQRT(0.62D2)
      CASE (17)
         value = -0.265182525D9 / 0.16384D5 * 14*x ** 13 * SQRT(0.66D2) + 0.778439025D9 & 
              / 0.16384D5 * 12*x ** 11 * SQRT(0.66D2) - 0.885809925D9 / 0.16384D5 * & 
              10*x ** 9 * SQRT(0.66D2) + 0.492116625D9 / 0.16384D5 * 8*x ** 7 * SQRT(0.66D2) - &
              0.137792655D9 / 0.16384D5 * 6*x ** 5 * SQRT(0.66D2) + 0.17972955D8 / 0.16384D5 * &
              4*x ** 3 * SQRT(0.66D2) - 0.855855D6 / 0.16384D5 * 2*x * SQRT(0.66D2)
      CASE (18)
         value = -0.19D2 / 0.2048D4 * (3411705D0 * 15*x ** 14 - 10855425D0 * &
              13*x ** 12 + 13656825D0 * 11*x ** 10 - 8633625D0 * 9*x ** 8 + 2877875D0 * 7*x ** 6 - &
              483483D0 * 5*x ** 4 + 35035 * 3*x ** 2 - 715) * SQRT(0.70D2)
      CASE (19)
         value = -0.2029964475D10 / 0.32768D5 * 16*x ** 15 * SQRT(0.74D2) + 0.869984775D9 &
              / 0.4096D4 * 14*x ** 13 * SQRT(0.74D2) - 0.2399048925D10 &
              / 0.8192D4 * 12*x ** 11 * SQRT(0.74D2) + 0.851275425D9 / 0.4096D4 * 10*x** 9 * &
              SQRT(0.74D2) - 0.1320944625D10 / 0.16384D5 * 8*x ** 7 * SQRT(0.74D2) + &
              0.68493425D8 / 0.4096D4 * 6*x ** 5 * SQRT(0.74D2) - 0.13698685D8 / 0.8192D4 * &
              4*x ** 3 * SQRT(0.74D2) + 0.255255D6 / 0.4096D4 * 2*x * SQRT(0.74D2)
      CASE (20)
         value = -0.63D2 / 0.65536D5 * (126233085D0 * 17*x ** 16 - &
              463991880D0 * 15*x ** 14 + 695987820D0 * 13*x ** 12 - 548354040D0 * 11*x ** 10 + &
              243221550D0 * 9*x ** 8 - 60386040D0 * 7*x ** 6 + 7827820D0 * 5*x ** 4 - 447304D0 * &
              3*x ** 2 + 7293) * SQRT(0.78D2)
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         stop ' ddvarphi > 20 ?'
         IF (x==1 .OR. x==-1) THEN
            ! TEMP SOLUTION
            ! Try to interpolate value of function
            vp = 4*((1-(x+dx)**2)*dPhi(i,(x+dx))+2*(x+dx)*Phi(i,(x+dx)))/(1-(x+dx)**2)**2
            vm = 4*((1-(x-dx)**2)*dPhi(i,(x-dx))+2*(x-dx)*Phi(i,(x-dx)))/(1-(x-dx)**2)**2
            value = (vp+vm)/2
            ! WRITE (*,*) i,x,vp,vm,value
         ELSE
            value = 4*((1-x**2)*dPhi(i,x)+2*x*Phi(i,x))/(1-x**2)**2
            value = 4*((1-x**2)*dPhi(i,x)+2*x*Phi(i,x))/(1-x**2)**2
         END IF
#endif
      END SELECT
    END FUNCTION ddVarPhi

    
!------------------------------------------------------------------------------
!>    Function LegendreP returns value of l,th Legendre polynomial
!>    for point x.
!
!>    Value of legendre polynomial is precalculated for l=<20 and calculated from 
!>    recursion for l>20,
!
!>    P(i+1,x)=1/(1+i)*((2*i+1)*x*P(i,x)+i*P(i-1,x)), 
!>    where P(0,x)=1, P(1,x)=x.
!------------------------------------------------------------------------------
    PURE RECURSIVE FUNCTION LegendreP(l,x) RESULT(value) 
!------------------------------------------------------------------------------
! 
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: l
!      INPUT: parameter of Legendre polynomial
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate Legendre polynomial
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of legendre polynomial l at point x i.e
!       value = P(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: l
      REAL (KIND=dp), INTENT(IN) :: x
      ! Return value
      REAL (KIND=dp) :: value
      
      ! Internal variables
      REAL (KIND=dp) :: P_l_1, P_l, PT
      INTEGER :: k 

      ! First 20 Legendre polynomials are precalculated. These are 
      ! all generated with Maple
      SELECT CASE (l)
#ifdef DEBUG_PBASIS
      CASE (:-1)
         CALL Fatal('PElementBase::LegendreP','LegendreP not defined for l < 0')
#endif
      CASE (0)
         value = 1      
      CASE (1)
         value = x
      CASE (2)
         value = -0.1D1 / 0.2D1 + 0.3D1 / 0.2D1 * x ** 2
      CASE (3)
         value = 0.5D1 / 0.2D1 * x ** 3 - 0.3D1 / 0.2D1 * x
      CASE (4)
         value = 0.3D1 / 0.8D1 + 0.35D2 / 0.8D1 * x ** 4 - 0.15D2 / 0.4D1 * x ** 2
      CASE (5)
         value = 0.63D2 / 0.8D1 * x ** 5 - 0.35D2 / 0.4D1 * x ** 3 + 0.15D2 / 0.8D1 * x
      CASE (6)
         value = -0.5D1 / 0.16D2 + 0.231D3 / 0.16D2 * x ** 6 - 0.315D3 / 0.16D2 * & 
              x ** 4 + 0.105D3 / 0.16D2 * x ** 2
      CASE (7)
         value = 0.429D3 / 0.16D2 * x ** 7 - 0.693D3 / 0.16D2 * x ** 5 + 0.315D3 / &
              0.16D2 * x ** 3 - 0.35D2 / 0.16D2 * x
      CASE (8)
         value = 0.35D2 / 0.128D3 + 0.6435D4 / 0.128D3 * x ** 8 - 0.3003D4 &
              / 0.32D2 * x ** 6 + 0.3465D4 / 0.64D2 * x ** 4 - 0.315D3 / 0.32D2 * x ** 2
      CASE (9)
         value = 0.12155D5 / 0.128D3 * x ** 9 - 0.6435D4 / 0.32D2 * x ** 7 + & 
              0.9009D4 / 0.64D2 * x ** 5 - 0.1155D4 / 0.32D2 * x ** 3 + 0.315D3 / 0.128D3 * x
      CASE (10)
         value = -0.63D2 / 0.256D3 + 0.46189D5 / 0.256D3 * x ** 10 - 0.109395D6 / &
              0.256D3 * x ** 8 + 0.45045D5 / 0.128D3 * x ** 6 - 0.15015D5/ 0.128D3 * &
              x ** 4 + 0.3465D4 / 0.256D3 * x ** 2
      CASE (11)
         value = 0.88179D5 / 0.256D3 * x ** 11 - 0.230945D6 / 0.256D3 * x ** 9 + &
              0.109395D6 / 0.128D3 * x ** 7 - 0.45045D5 / 0.128D3 * x ** 5 + 0.15015D5 / &
              0.256D3 * x ** 3 - 0.693D3 / 0.256D3 * x
      CASE (12)
         value = 0.231D3 / 0.1024D4 + 0.676039D6 / 0.1024D4 * x ** 12 - 0.969969D6 / &
              0.512D3 * x ** 10 + 0.2078505D7 / 0.1024D4 * x ** 8 - 0.255255D6 / 0.256D3 * &
              x ** 6 + 0.225225D6 / 0.1024D4 * x ** 4 - 0.9009D4 / 0.512D3 * x ** 2
      CASE (13)
         value = 0.1300075D7 / 0.1024D4 * x ** 13 - 0.2028117D7 / 0.512D3 * x ** 11 + &
              0.4849845D7 / 0.1024D4 * x ** 9 - 0.692835D6 / 0.256D3 * x ** 7 + 0.765765D6 / &
              0.1024D4 * x ** 5 - 0.45045D5 / 0.512D3 * x ** 3 + 0.3003D4 / 0.1024D4 * x
      CASE (14)
         value = -0.429D3 / 0.2048D4 + 0.5014575D7 / 0.2048D4 * x ** 14 - 0.16900975D8 / &
              0.2048D4 * x ** 12 + 0.22309287D8 / 0.2048D4 * x ** 10 - 0.14549535D8 / 0.2048D4 * &
              x ** 8 + 0.4849845D7 / 0.2048D4 * x ** 6 - 0.765765D6 / 0.2048D4 * x ** 4 + &
              0.45045D5 / 0.2048D4 * x ** 2
      CASE (15)     
         value = 0.9694845D7 / 0.2048D4 * x ** 15 - 0.35102025D8 / 0.2048D4 * x ** 13 +&
              0.50702925D8 / 0.2048D4 * x ** 11 - 0.37182145D8 / 0.2048D4 * x ** 9 + &
              0.14549535D8 / 0.2048D4 * x ** 7 - 0.2909907D7 / 0.2048D4 * x ** 5 + & 
              0.255255D6 / 0.2048D4 * x ** 3 - 0.6435D4 / 0.2048D4 * x
      CASE(16)
         value = 0.6435D4 / 0.32768D5 + 0.300540195D9 / 0.32768D5 * x ** 16  - &
              0.145422675D9 / 0.4096D4 * x ** 14 + 0.456326325D9 / 0.8192D4 * x ** 12 -&
              0.185910725D9 / 0.4096D4 * x ** 10 + 0.334639305D9 / 0.16384D5 * x ** 8 - &
              0.20369349D8 / 0.4096D4 * x ** 6 + 0.4849845D7 / 0.8192D4 * x ** 4 - &
              0.109395D6 / 0.4096D4 * x ** 2
      CASE (17)
         value = 0.583401555D9 / 0.32768D5 * x ** 17 - 0.300540195D9 / 0.4096D4 *&
              x ** 15 + 0.1017958725D10 / 0.8192D4 * x ** 13 - 0.456326325D9 / 0.4096D4 *&
              x ** 11 + 0.929553625D9 / 0.16384D5 * x ** 9 - 0.66927861D8 / 0.4096D4 * &
              x ** 7 + 0.20369349D8 / 0.8192D4 * x ** 5 - 0.692835D6 / 0.4096D4 * x ** 3 +&
              0.109395D6 / 0.32768D5 * x
      CASE (18)
         value = -0.12155D5 / 0.65536D5 + 0.2268783825D10 / 0.65536D5 * x ** 18 - &
              0.9917826435D10 / 0.65536D5 * x ** 16 + 0.4508102925D10 / 0.16384D5 * &
              x ** 14 - 0.4411154475D10 / 0.16384D5 * x ** 12 + 0.5019589575D10 / &
              0.32768D5 * x ** 10 - 0.1673196525D10 / 0.32768D5 * x ** 8 + 0.156165009D9 &
              / 0.16384D5 * x ** 6 - 0.14549535D8 / 0.16384D5 * x ** 4 + 0.2078505D7 / &
              0.65536D5 * x ** 2
      CASE (19)
         value = 0.4418157975D10 / 0.65536D5 * x ** 19 - 0.20419054425D11 / &
              0.65536D5 * x ** 17 + 0.9917826435D10 / 0.16384D5 * x ** 15 - &
              0.10518906825D11 / 0.16384D5 * x ** 13 + 0.13233463425D11 / 0.32768D5 * &
              x ** 11 - 0.5019589575D10 / 0.32768D5 * x ** 9 + 0.557732175D9 / 0.16384D5 * &
              x ** 7 - 0.66927861D8 / 0.16384D5 * x ** 5 + 0.14549535D8 / 0.65536D5 * &
              x ** 3 - 0.230945D6 / 0.65536D5 * x
      CASE (20)
         value = 0.46189D5 / 0.262144D6 + 0.34461632205D11 / 0.262144D6 * x ** 20 - &
              0.83945001525D11 / 0.131072D6 * x ** 18 + 0.347123925225D12 / 0.262144D6 *&
              x ** 16 - 0.49589132175D11 / 0.32768D5 * x ** 14 + 0.136745788725D12 / &
              0.131072D6 * x ** 12 - 0.29113619535D11 / 0.65536D5 * x ** 10 + &
              0.15058768725D11 / 0.131072D6 * x ** 8 - 0.557732175D9 / 0.32768D5 * &
              x ** 6 + 0.334639305D9 / 0.262144D6 * x ** 4 - 0.4849845D7 / 0.131072D6 * x ** 2
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         ! Generate n:th Legendre polynomial
         
         ! Initialize first two legendre functions 
         ! P(19,x)=... 
         P_l_1=LegendreP(19,x)
         ! P(20,x)=...
         P_l=LegendreP(20,x)
       
         ! Generate (k+1):th legendre polynomial
         DO k=20,(l-1)
            PT = (1d0/(k+1))*((2*k+1)*x*P_l - k*P_l_1)
            ! Advance to next legendre polynomial
            P_l_1=P_l
            P_l=PT
         END DO
     
         value = P_l
#endif
      END SELECT
      ! Value contains LegendreP(l,x)
    END FUNCTION LegendreP


!------------------------------------------------------------------------------
!>    Function dLegendreP returns value of derivative of l:th Legendre polynomial
!>    at point x.
!
!>    Value of legendre polynomial is precalculated for l=<20 and calculated from 
!>    recursion for l>20,
!
!>    P,(l+1,x)=x*P,(l,x)+(l+1)*P(l,x), 
!>    where P,(0,x)=0, P(1,x)=1.
!------------------------------------------------------------------------------
    PURE RECURSIVE FUNCTION dLegendreP(l,x) RESULT(value)
!------------------------------------------------------------------------------
! 
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: l
!      INPUT: parameter of Legendre polynomial
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate Legendre polynomial
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of derivative of legendre polynomial l at point x i.e
!       value = P,(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: l
      REAL (KIND=dp), INTENT(IN) :: x
      ! Return value
      REAL (KIND=dp) :: value
      
      ! Internal variables
      REAL (KIND=dp) :: P_l, dP_l, dPT
      INTEGER :: k 

      SELECT CASE(l)
#ifdef DEBUG_PBASIS
      CASE (:-1)
         CALL Fatal('PElementBase::dLegendreP','dLegendreP not defined for l < 0')
#endif
      CASE (0)
         value = 0
      CASE (1)
         value = 1
      CASE (2)
         value = 3 * x
      CASE (3)
         value = 0.15D2 / 0.2D1 * x ** 2 - 0.3D1 / 0.2D1
      CASE (4)
         value = 0.35D2 / 0.2D1 * x ** 3 - 0.15D2 / 0.2D1 * x
      CASE (5)
         value = 0.315D3 / 0.8D1 * x ** 4 - 0.105D3 / 0.4D1 * x ** 2 + 0.15D2 / 0.8D1
      CASE (6)
         value = 0.693D3 / 0.8D1 * x ** 5 - 0.315D3 / 0.4D1 * x ** 3 + 0.105D3 &
              / 0.8D1 * x
      CASE (7)     
         value = 0.3003D4 / 0.16D2 * x ** 6 - 0.3465D4 / 0.16D2 * x ** 4 + &
              0.945D3 / 0.16D2 * x ** 2 - 0.35D2 / 0.16D2
      CASE (8)     
         value = 0.6435D4 / 0.16D2 * x ** 7 - 0.9009D4 / 0.16D2 * x ** 5 + &
              0.3465D4 / 0.16D2 * x ** 3 - 0.315D3 / 0.16D2 * x
      CASE (9)   
         value = 0.109395D6 / 0.128D3 * x ** 8 - 0.45045D5 / 0.32D2 * x ** 6 +&
              0.45045D5 / 0.64D2 * x ** 4 - 0.3465D4 / 0.32D2 * x ** 2 + 0.315D3 / 0.128D3
      CASE (10)
         value = 0.230945D6 / 0.128D3 * x ** 9 - 0.109395D6 / 0.32D2 * x ** 7 + &
              0.135135D6 / 0.64D2 * x ** 5 - 0.15015D5 / 0.32D2 * x ** 3 + &
              0.3465D4 / 0.128D3 * x
      CASE (11)
         value = 0.969969D6 / 0.256D3 * x ** 10 - 0.2078505D7 / 0.256D3 * x ** 8 +&
              0.765765D6 / 0.128D3 * x ** 6 - 0.225225D6 / 0.128D3 * x ** 4 + 0.45045D5 / &
              0.256D3 * x ** 2 - 0.693D3 / 0.256D3
      CASE (12) 
         value = 0.2028117D7 / 0.256D3 * x ** 11 - 0.4849845D7 / 0.256D3 * x ** 9 + &
              0.2078505D7 / 0.128D3 * x ** 7 - 0.765765D6 / 0.128D3 * x ** 5 + &
              0.225225D6 / 0.256D3 * x ** 3 - 0.9009D4 / 0.256D3 * x
      CASE (13) 
         value = 0.16900975D8 / 0.1024D4 * x ** 12 - 0.22309287D8 / 0.512D3 * &
              x ** 10 + 0.43648605D8 / 0.1024D4 * x ** 8 - 0.4849845D7 / 0.256D3 * &
              x ** 6 + 0.3828825D7 / 0.1024D4 * x ** 4 - 0.135135D6 / 0.512D3 * x ** 2 +&
              0.3003D4 / 0.1024D4
      CASE (14)   
         value = 0.35102025D8 / 0.1024D4 * x ** 13 - 0.50702925D8 / 0.512D3 * &
              x ** 11 + 0.111546435D9 / 0.1024D4 * x ** 9 - 0.14549535D8 / 0.256D3 * &
              x ** 7 + 0.14549535D8 / 0.1024D4 * x ** 5 - 0.765765D6 / 0.512D3 * &
              x ** 3 + 0.45045D5 / 0.1024D4 * x
      CASE (15)     
         value = 0.145422675D9 / 0.2048D4 * x ** 14 - 0.456326325D9 / 0.2048D4 * &
              x ** 12 + 0.557732175D9 / 0.2048D4 * x ** 10 - 0.334639305D9 / 0.2048D4 * &
              x ** 8 + 0.101846745D9 / 0.2048D4 * x ** 6 - 0.14549535D8 / 0.2048D4 * &
              x ** 4 + 0.765765D6 / 0.2048D4 * x ** 2 - 0.6435D4 / 0.2048D4
      CASE (16)
         value = 0.300540195D9 / 0.2048D4 * x ** 15 - 0.1017958725D10 / 0.2048D4 * &
              x ** 13 + 0.1368978975D10 / 0.2048D4 * x ** 11 - 0.929553625D9 / 0.2048D4 * &
              x ** 9 + 0.334639305D9 / 0.2048D4 * x ** 7 - 0.61108047D8 / 0.2048D4 * &
              x ** 5 + 0.4849845D7 / 0.2048D4 * x ** 3 - 0.109395D6 / 0.2048D4 * x
      CASE (17)      
         value = 0.9917826435D10 / 0.32768D5 * x ** 16 - 0.4508102925D10 / 0.4096D4 * &
              x ** 14 + 0.13233463425D11 / 0.8192D4 * x ** 12 - &
              0.5019589575D10 / 0.4096D4 * x ** 10 + 0.8365982625D10 / 0.16384D5 * x ** 8 -&
              0.468495027D9 / 0.4096D4 * x ** 6 + 0.101846745D9 / 0.8192D4 * x ** 4 - &
              0.2078505D7 / 0.4096D4 * x ** 2 + 0.109395D6 / 0.32768D5
      CASE (18)
         value = 0.20419054425D11 / 0.32768D5 * x ** 17 - 0.9917826435D10 / &
              0.4096D4 * x ** 15 + 0.31556720475D11 / 0.8192D4 * x ** 13 - &
              0.13233463425D11 / 0.4096D4 * x ** 11 + 0.25097947875D11 / 0.16384D5 * &
              x ** 9 - 0.1673196525D10 / 0.4096D4 * x ** 7 + 0.468495027D9 / 0.8192D4 * &
              x ** 5 - 0.14549535D8 / 0.4096D4 * x ** 3 + 0.2078505D7 / 0.32768D5 * x
      CASE (19)
         value = 0.83945001525D11 / 0.65536D5 * x ** 18 - 0.347123925225D12 / &
              0.65536D5 * x ** 16 + 0.148767396525D12 / 0.16384D5 * x ** 14 - &
              0.136745788725D12 / 0.16384D5 * x ** 12 + 0.145568097675D12 / 0.32768D5 * &
              x ** 10 - 0.45176306175D11 / 0.32768D5 * x ** 8 + 0.3904125225D10 / &
              0.16384D5 * x ** 6 - 0.334639305D9 / 0.16384D5 * x ** 4 + 0.43648605D8 / &
              0.65536D5 * x ** 2 - 0.230945D6 / 0.65536D5
      CASE (20)   
         value = 0.172308161025D12 / 0.65536D5 * x ** 19 - 0.755505013725D12 / &
              0.65536D5 * x ** 17 + 0.347123925225D12 / 0.16384D5 * x ** 15 - & 
              0.347123925225D12 / 0.16384D5 * x ** 13 + 0.410237366175D12 / 0.32768D5 * &
              x ** 11 - 0.145568097675D12 / 0.32768D5 * x ** 9 + 0.15058768725D11 / &
              0.16384D5 * x ** 7 - 0.1673196525D10 / 0.16384D5 * x ** 5 + 0.334639305D9 /&
              0.65536D5 * x ** 3 - 0.4849845D7 / 0.65536D5 * x
#ifdef DEBUG_PBASIS
      CASE DEFAULT
         ! Generate derivative of n:th Legendre polynomial

         ! Initialize derivative of legendre polynomial for l=20
         ! P,(20,x)=... 
         dP_l = dLegendreP(20,x)

         ! Generate derivative of (k+1):th legendre polynomial
         DO k=20,(l-1)
            ! P(k,x)=...
            P_l=LegendreP(k,x)
            dPT = x*dP_l+(k+1)*P_l
            ! Advance to next legendre polynomial
            dP_l=dPT
         END DO
     
         value = dP_l
#endif
      END SELECT
      ! Value now contains P,(l,x)
    END FUNCTION dLegendreP

!------------------------------------------------------------------------------
!>    Function ddLegendreP returns value of 2nd derivative of l:th Legendre polynomial
!>    at point x.
!
!>    Value of legendre polynomial is precalculated for l=<20 and calculated from 
!>    recursion for l>20,
!
!>    P,(l+1,x)=x*P,(l,x)+(l+1)*P(l,x), 
!>    where P,(0,x)=0, P(1,x)=1.
!------------------------------------------------------------------------------
    PURE RECURSIVE FUNCTION ddLegendreP(l,x) RESULT(value)
!------------------------------------------------------------------------------
! 
!  ARGUMENTS:
!    INTEGER, INTENT(IN) :: l
!      INPUT: parameter of Legendre polynomial
!
!    REAL(Kind=dp), INTENT(IN) :: x
!      INPUT: point at which to evaluate Legendre polynomial
!
!  FUNCTION VALUE:
!    REAL(Kind=dp) :: value
!       value of derivative of legendre polynomial l at point x i.e
!       value = P,(i,x)
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: l
      REAL (KIND=dp), INTENT(IN) :: x
      ! Return value
      REAL (KIND=dp) :: value
      
      ! Internal variables
      REAL (KIND=dp) :: P_l, dP_l, dPT
      INTEGER :: k 

      SELECT CASE(l)
#ifdef DEBUG_PBASIS
      CASE (:-1)
         CALL Fatal('PElementBase::dLegendreP','dLegendreP not defined for l < 0')
#endif
      CASE (0)
         value = 0
      CASE (1)
         value = 0
      CASE (2)
         value = 3
      CASE (3)
         value = 0.15D2 / 0.2D1 * 2*x
      CASE (4)
         value = 0.35D2 / 0.2D1 * 3*x**2 - 0.15D2 / 0.2D1
      CASE (5)
         value = 0.315D3 / 0.8D1 * 4*x ** 3 - 0.105D3 / 0.4D1 * 2*x
      CASE (6)
         value = 0.693D3 / 0.8D1 * 5*x ** 4 - 0.315D3 / 0.4D1 * 3*x ** 2 + 0.105D3/0.8D1
      CASE (7)     
         value = 0.3003D4 / 0.16D2 * 6*x ** 5 - 0.3465D4 / 0.16D2 * 4*x ** 3 + &
              0.945D3 / 0.16D2 * 2*x
      CASE (8)     
         value = 0.6435D4 / 0.16D2 * 7*x ** 6 - 0.9009D4 / 0.16D2 * 5*x ** 4 + &
              0.3465D4 / 0.16D2 * 3*x ** 2 - 0.315D3 / 0.16D2
      CASE (9)   
         value = 0.109395D6 / 0.128D3 * 8*x ** 7 - 0.45045D5 / 0.32D2 * 6*x ** 5 +&
              0.45045D5 / 0.64D2 * 4*x ** 3 - 0.3465D4 / 0.32D2 * 2*x
      CASE (10)
         value = 0.230945D6 / 0.128D3 * 9*x ** 8 - 0.109395D6 / 0.32D2 * 7*x ** 6 + &
              0.135135D6 / 0.64D2 * 5*x ** 4 - 0.15015D5 / 0.32D2 * 3*x ** 2 + &
              0.3465D4 / 0.128D3
      CASE (11)
         value = 0.969969D6 / 0.256D3 * 10*x ** 9 - 0.2078505D7 / 0.256D3 * 8*x ** 7 +&
              0.765765D6 / 0.128D3 * 6*x ** 5 - 0.225225D6 / 0.128D3 * 4*x ** 3 + 0.45045D5 / &
              0.256D3 * 2*x
      CASE (12) 
         value = 0.2028117D7 / 0.256D3 * 11*x ** 10 - 0.4849845D7 / 0.256D3 * 9*x ** 8 + &
              0.2078505D7 / 0.128D3 * 7*x ** 6 - 0.765765D6 / 0.128D3 * 5*x ** 4 + &
              0.225225D6 / 0.256D3 * 3*x ** 2 - 0.9009D4 / 0.256D3
      CASE (13) 
         value = 0.16900975D8 / 0.1024D4 * 12*x ** 11 - 0.22309287D8 / 0.512D3 * &
              10*x ** 9 + 0.43648605D8 / 0.1024D4 * 8*x ** 7 - 0.4849845D7 / 0.256D3 * &
              6*x ** 5 + 0.3828825D7 / 0.1024D4 * 4*x ** 3 - 0.135135D6 / 0.512D3 * 2*x
      CASE (14)   
         value = 0.35102025D8 / 0.1024D4 * 13*x ** 12 - 0.50702925D8 / 0.512D3 * &
              11*x ** 10 + 0.111546435D9 / 0.1024D4 * 9*x ** 8 - 0.14549535D8 / 0.256D3 * &
              7*x ** 6 + 0.14549535D8 / 0.1024D4 * 5*x ** 4 - 0.765765D6 / 0.512D3 * &
              3*x ** 2 + 0.45045D5 / 0.1024D4
      CASE (15)     
         value = 0.145422675D9 / 0.2048D4 * 14*x ** 13 - 0.456326325D9 / 0.2048D4 * &
              12*x ** 11 + 0.557732175D9 / 0.2048D4 * 10*x ** 9 - 0.334639305D9 / 0.2048D4 * &
              8*x ** 7 + 0.101846745D9 / 0.2048D4 * 6*x ** 5 - 0.14549535D8 / 0.2048D4 * &
              4*x ** 3 + 0.765765D6 / 0.2048D4 * 2*x
      CASE (16)
         value = 0.300540195D9 / 0.2048D4 * 15*x ** 14 - 0.1017958725D10 / 0.2048D4 * &
              13*x ** 12 + 0.1368978975D10 / 0.2048D4 * 11*x ** 10 - 0.929553625D9 / 0.2048D4 * &
              9*x ** 8 + 0.334639305D9 / 0.2048D4 * 7*x ** 6 - 0.61108047D8 / 0.2048D4 * &
              5*x ** 4 + 0.4849845D7 / 0.2048D4 * 3*x ** 2 - 0.109395D6 / 0.2048D4
      CASE (17)      
         value = 0.9917826435D10 / 0.32768D5 * 16*x ** 15 - 0.4508102925D10 / 0.4096D4 * &
              14*x ** 13 + 0.13233463425D11 / 0.8192D4 * 12*x ** 11 - &
              0.5019589575D10 / 0.4096D4 * 10*x ** 9 + 0.8365982625D10 / 0.16384D5 * 8*x ** 7 -&
              0.468495027D9 / 0.4096D4 * 6*x ** 5 + 0.101846745D9 / 0.8192D4 * 4*x ** 3 - &
              0.2078505D7 / 0.4096D4 * 2*x
      CASE (18)
         value = 0.20419054425D11 / 0.32768D5 * 17*x ** 16 - 0.9917826435D10 / &
              0.4096D4 * 15*x ** 14 + 0.31556720475D11 / 0.8192D4 * 13*x ** 12 - &
              0.13233463425D11 / 0.4096D4 * 11*x ** 10 + 0.25097947875D11 / 0.16384D5 * &
              9*x ** 8 - 0.1673196525D10 / 0.4096D4 * 7*x ** 6 + 0.468495027D9 / 0.8192D4 * &
              5*x ** 4 - 0.14549535D8 / 0.4096D4 * 3*x ** 2 + 0.2078505D7 / 0.32768D5
      CASE (19)
         value = 0.83945001525D11 / 0.65536D5 * 18*x ** 17 - 0.347123925225D12 / &
              0.65536D5 * 16*x ** 15 + 0.148767396525D12 / 0.16384D5 * 14*x ** 13 - &
              0.136745788725D12 / 0.16384D5 * 12*x ** 11 + 0.145568097675D12 / 0.32768D5 * &
              10*x ** 9 - 0.45176306175D11 / 0.32768D5 * 8*x ** 7 + 0.3904125225D10 / &
              0.16384D5 * 6*x ** 5 - 0.334639305D9 / 0.16384D5 * 4*x ** 3 + 0.43648605D8 / &
              0.65536D5 * 2*x
      CASE (20)   
         value = 0.172308161025D12 / 0.65536D5 * 19*x ** 18 - 0.755505013725D12 / &
              0.65536D5 * 17*x ** 16 + 0.347123925225D12 / 0.16384D5 * 15*x ** 14 - & 
              0.347123925225D12 / 0.16384D5 * 13*x ** 12 + 0.410237366175D12 / 0.32768D5 * &
              11*x ** 10 - 0.145568097675D12 / 0.32768D5 * 9*x ** 8 + 0.15058768725D11 / &
              0.16384D5 * 7*x ** 6 - 0.1673196525D10 / 0.16384D5 * 5*x ** 4 + 0.334639305D9 /&
              0.65536D5 * 3*x ** 2 - 0.4849845D7 / 0.65536D5
      CASE DEFAULT
#ifdef DEBUG_PBASIS
         ! Generate derivative of n:th Legendre polynomial

         STOP 'No 2nd derivative for Legendre > 20'

         ! Initialize derivative of legendre polynomial for l=20
         ! P,(20,x)=... 
         dP_l = dLegendreP(20,x)

         ! Generate derivative of (k+1):th legendre polynomial
         DO k=20,(l-1)
            ! P(k,x)=...
            P_l=LegendreP(k,x)
            dPT = x*dP_l+(k+1)*P_l
            ! Advance to next legendre polynomial
            dP_l=dPT
         END DO
     
         value = dP_l
#endif
      END SELECT
      ! Value now contains P,(l,x)
    END FUNCTION ddLegendreP

    ! Function value = x^n
    PURE FUNCTION toExp(x,n) RESULT(value)
      IMPLICIT NONE
      
      REAL(KIND=dp), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: n
      REAL(KIND=dp) :: value

      ! Anything to 0 is 1
      IF (n == 0) THEN
         value = 1
      ! 0 to anything not 0 is 0
      ELSE IF (x == 0) THEN
         value = 0
      ELSE 
         value = x**n
      END IF
    END FUNCTION toExp


    RECURSIVE FUNCTION Product2ndDerivatives(n,g,dg,ddg,dim,level) RESULT(grad)
      IMPLICIT NONE

      INTEGER :: n,dim,level
      REAL(KIND=dp) :: g(:), dg(:,:), ddg(:,:,:), grad(dim,dim)

      INTEGER :: i,j,p,q
      REAL(KIND=dp) :: h,dh(dim),ddh(dim,dim),l,dl(dim),ddl(dim,dim)

!     df1*f2 + f1*df2
!     ddf1(p,q)*f2 + df1(i)*df2(j) + df1(j)*df2(i) + f1*ddf2(p,q)

      j = 0
      DO i=1,n-1,2
        h  = g(i)
        l  = g(i+1)
        dh = dg(i,:)
        dl = dg(i+1,:)
        ddh = ddg(i,:,:)
        ddl = ddg(i+1,:,:)
        j = j + 1
        g(j) = h*l
        dg(j,:) = dh*l + h*dl

        ddg(j,1,1) = ddh(1,1)*l + 2*dh(1)*dl(1) + h*ddl(1,1)
        ddg(j,1,2) = ddh(1,2)*l + dh(1)*dl(2) + dh(2)*dl(1) + h*ddl(1,2)
        ddg(j,2,2) = ddh(2,2)*l + 2*dh(2)*dl(2) + h*ddl(2,2)
        IF(dim>2) THEN
          ddg(j,1,3) = ddh(1,3)*l + dh(1)*dl(3) + dh(3)*dl(1) + h*ddl(1,3)
          ddg(j,2,3) = ddh(2,3)*l + dh(2)*dl(3) + dh(3)*dl(2) + h*ddl(2,3)
          ddg(j,3,3) = ddh(3,3)*l + 2*dh(3)*dl(3) + h*ddl(3,3)
        END IF
      END DO
      IF(mod(n,2) /= 0) THEN
        h  = g(j)
        l  = g(n)
        dh = dg(j,:)
        dl = dg(n,:)
        ddh = ddg(j,:,:)
        ddl = ddg(n,:,:)
        g(j) = h*l
        dg(j,:) = dh*l + h*dl

        ddg(j,1,1) = ddh(1,1)*l + 2*dh(1)*dl(1) + h*ddl(1,1)
        ddg(j,1,2) = ddh(1,2)*l + dh(1)*dl(2) + dh(2)*dl(1) + h*ddl(1,2)
        ddg(j,2,2) = ddh(2,2)*l + 2*dh(2)*dl(2) + h*ddl(2,2)
        IF(dim>2) THEN
          ddg(j,1,3) = ddh(1,3)*l + dh(1)*dl(3) + dh(3)*dl(1) + h*ddl(1,3)
          ddg(j,2,3) = ddh(2,3)*l + dh(2)*dl(3) + dh(3)*dl(2) + h*ddl(2,3)
          ddg(j,3,3) = ddh(3,3)*l + 2*dh(3)*dl(3) + h*ddl(3,3)
        END IF
      END IF

      IF(n/2>1) ddh=Product2ndDerivatives(n/2,g,dg,ddg,dim,level+1)
      IF(Level==0) THEN
        grad = 0
        DO p=1,dim
          DO q=p,dim
            grad(p,q) = ddg(1,p,q)
          END DO
        END DO
      END IF

    END FUNCTION Product2ndDerivatives


END MODULE PElementBase

!> \}
