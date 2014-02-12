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

      value = phi(i,phiPar)
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
         value = 1d0/4*(1-u)*(1-v)
      CASE (2)
         value = 1d0/4*(1+u)*(1-v)
      CASE (3)
         value = 1d0/4*(1+u)*(1+v)
      CASE (4)
         value = 1d0/4*(1-u)*(1+v)
      CASE DEFAULT
         CALL Fatal('PElementBase::QuadNodalPBasis', 'Unknown node for quadrilateral')
      END SELECT
    END FUNCTION QuadNodalPBasis


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
         grad(1) = -1d0/4*(1-v)
         grad(2) = -1d0/4*(1-u)
      CASE (2)
         grad(1) = 1d0/4*(1-v)
         grad(2) = -1d0/4*(1+u)
      CASE (3)
         grad(1) = 1d0/4*(1+v)
         grad(2) = 1d0/4*(1+u)
      CASE (4)
         grad(1) = -1d0/4*(1+v)
         grad(2) = 1d0/4*(1-u)
      CASE DEFAULT
         CALL Fatal('PElementBase::dQuadNodalPBasis', 'Unknown node for quadrilateral')
      END SELECT
    END FUNCTION dQuadNodalPBasis


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
    END FUNCTION QuadEdgePBasis

!------------------------------------------------------------------------------
!>     Quadrilateral edge basis at point (u,v), which is compatible with
!>     basis for pyramidal 3d element square face
!------------------------------------------------------------------------------
    FUNCTION QuadPyraEdgePBasis(edge, i, u, v, invertEdge) RESULT(value)
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
      
      ! Parameters 
      INTEGER, INTENT(IN) :: edge, i
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v
      
      ! Variables
      INTEGER :: nodes(2),tmp
      REAL (KIND=dp) :: Na, Nb, La, Lb, value
      LOGICAL :: invert 
      
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge
      
      ! Check edge parameter validity
      IF (edge < 1 .OR. edge > 4) THEN
         CALL Fatal('PElementBase::QuadPyraEdgePBasis', 'Unknown edge for quadrilateral')
      END IF

      nodes(1:2) = getQuadEdgeMap(edge)
      ! Get bilinear nodal function values
      Na = QuadNodalPBasis(nodes(1),u,v)
      Nb = QuadNodalPBasis(nodes(2),u,v)
      
      ! Invert edge direction if needed
      IF (invert) THEN
         tmp = nodes(1)
         nodes(1) = nodes(2)
         nodes(2) = tmp
      END IF

      ! Get affine function values for edge direction
      La = QuadL(nodes(1),u,v)
      Lb = QuadL(nodes(2),u,v)

      value = Na*Nb*varPhi(i,Lb-La)
    END FUNCTION QuadPyraEdgePBasis

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
    END FUNCTION dQuadEdgePBasis

!------------------------------------------------------------------------------
!>     Gradient of quadrilateral edge basis at point (u,v) which is
!>     compatible with pyramidal 3d element square face edges.
!------------------------------------------------------------------------------
    FUNCTION dQuadPyraEdgePBasis(edge, i, u, v, invertEdge) RESULT(grad)
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
      LOGICAL, OPTIONAL :: invertEdge
      REAL (KIND=dp), INTENT(IN) :: u,v

      ! Variables
      INTEGER :: nodes(2), tmp
      REAL (KIND=dp) :: Na, Nb, La, Lb, vPhi
      REAL (KIND=dp), DIMENSION(2) :: dNa, dNb, dLa, dLb, grad
      LOGICAL :: invert

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      IF (edge < 1 .OR. edge > 4) THEN
         CALL Fatal('PElementBase::dQuadEdgePBasis', 'Unknown edge for quadrilateral')
      END IF
      
      nodes(1:2) = getQuadEdgeMap(edge)
      ! Get bilinear nodal function values and their gradients
      Na = QuadNodalPBasis(nodes(1),u,v)
      Nb = QuadNodalPBasis(nodes(2),u,v)
      dNa = dQuadNodalPBasis(nodes(1),u,v)
      dNb = dQuadNodalPBasis(nodes(2),u,v)

      ! Invert edge direction if needed
      IF (invert) THEN
         tmp = nodes(1)
         nodes(1) = nodes(2)
         nodes(2) = tmp
      END IF

      ! Get affine function values and their gradients for edge direction
      La = QuadL(nodes(1),u,v)
      Lb = QuadL(nodes(2),u,v)
      dLa = dQuadL(nodes(1),u,v)
      dLb = dQuadL(nodes(2),u,v)

      vPhi = varPhi(i,Lb-La)

      ! Calculate value of gradient from general form
      grad = 0
      grad = dNa*Nb*vPhi + Na*dNb*vPhi + Na*Nb*dVarPhi(i,Lb-La)*(dLb-dLa)
    END FUNCTION dQuadPyraEdgePBasis


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
      REAL(Kind=dp) :: La, Lb, Lc
      REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, grad
      
      ! Calculate value of function without direction and return
      ! if local numbering not present
      IF (.NOT. PRESENT(localNumbers)) THEN
         grad = 0
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
    END FUNCTION dQuadBubblePBasis


!------------------------------------------------------------------------------
!>     Defines linear functions for quadrilateral nodes. These are used in 
!>     calculation of changing parameters for bubbles of quadrilateral if 
!>     directional function values are requested. 
!------------------------------------------------------------------------------
    FUNCTION QuadL(which, u, v) RESULT(value)
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

      SELECT CASE (which)
      CASE (1)
         value = (2d0-u-v)/2d0
      CASE (2)
         value = (2d0+u-v)/2d0
      CASE (3)
         value = (2d0+u+v)/2d0
      CASE (4)
         value = (2d0-u+v)/2d0
      CASE DEFAULT
         value = 0.0_dp
         CALL Fatal('PElementBase::QuadL', 'Unknown helper function L for quad')
      END SELECT
    END FUNCTION QuadL


!------------------------------------------------------------------------------
!>     Defines gradients of linear functions for quadrilateral nodes. For use 
!>     see QuadL. 
!------------------------------------------------------------------------------
    FUNCTION dQuadL(which, u, v) RESULT(grad)
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
         grad(1:2) = [ 1d0/2, -1d0/2 ]
      CASE (3)
         grad(1:2) = [ 1d0/2, 1d0/2 ]
      CASE (4)
         grad(1:2) = [ -1d0/2, 1d0/2 ]
      CASE DEFAULT
         CALL Fatal('PElementBase::dQuadL', 'Unknown helper function dL for quad')
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
         value = 1d0/2*(1-u-v/SQRT(3d0))
      CASE (2)
         value = 1d0/2*(1+u-v/SQRT(3d0))
      CASE (3)
         value = v/SQRT(3d0)
      CASE DEFAULT
         CALL Fatal('PElementBase::TriangleNodalPBasis', 'Unknown node for triangle')
      END SELECT 
    END FUNCTION TriangleNodalPBasis


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

    ! 3D ELEMENTS

    FUNCTION BrickNodalPBasis(node, u, v, w) RESULT(value)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(Kind=dp), INTENT(IN) :: u,v,w

      ! Variables
      REAL(Kind=dp) :: value

      value = 0
      
      SELECT CASE (node)
      CASE (1)
         value = 1d0/8*(1-u)*(1-v)*(1-w)
      CASE (2)
         value = 1d0/8*(1+u)*(1-v)*(1-w)
      CASE (3)
         value = 1d0/8*(1+u)*(1+v)*(1-w)
      CASE (4)
         value = 1d0/8*(1-u)*(1+v)*(1-w)
      CASE (5)
         value = 1d0/8*(1-u)*(1-v)*(1+w)
      CASE (6)
         value = 1d0/8*(1+u)*(1-v)*(1+w)
      CASE (7)
         value = 1d0/8*(1+u)*(1+v)*(1+w)
      CASE (8)
         value = 1d0/8*(1-u)*(1+v)*(1+w)
      CASE DEFAULT
         CALL Fatal('PElementBase::BrickNodalPBasis','Unknown node for brick')
      END SELECT
    END FUNCTION BrickNodalPBasis

    FUNCTION dBrickNodalPBasis(node, u, v, w) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: node
      REAL(Kind=dp), INTENT(IN) :: u,v,w

      ! Variables
      REAL(Kind=dp) :: grad(3)

      grad = 0

      SELECT CASE (node)
      CASE (1)
         grad(1) = -1d0/8*(1-v)*(1-w) 
         grad(2) = -1d0/8*(1-u)*(1-w)
         grad(3) = -1d0/8*(1-u)*(1-v)
      CASE (2)
         grad(1) = 1d0/8*(1-v)*(1-w) 
         grad(2) = -1d0/8*(1+u)*(1-w)
         grad(3) = -1d0/8*(1+u)*(1-v)
      CASE (3)
         grad(1) = 1d0/8*(1+v)*(1-w) 
         grad(2) = 1d0/8*(1+u)*(1-w)
         grad(3) = -1d0/8*(1+u)*(1+v)
      CASE (4)
         grad(1) = -1d0/8*(1+v)*(1-w) 
         grad(2) = 1d0/8*(1-u)*(1-w)
         grad(3) = -1d0/8*(1-u)*(1+v)
      CASE (5)
         grad(1) = -1d0/8*(1-v)*(1+w) 
         grad(2) = -1d0/8*(1-u)*(1+w)
         grad(3) = 1d0/8*(1-u)*(1-v)
      CASE (6)
         grad(1) = 1d0/8*(1-v)*(1+w) 
         grad(2) = -1d0/8*(1+u)*(1+w)
         grad(3) = 1d0/8*(1+u)*(1-v)
      CASE (7)
         grad(1) = 1d0/8*(1+v)*(1+w) 
         grad(2) = 1d0/8*(1+u)*(1+w)
         grad(3) = 1d0/8*(1+u)*(1+v)
      CASE (8)
         grad(1) = -1d0/8*(1+v)*(1+w) 
         grad(2) = 1d0/8*(1-u)*(1+w)
         grad(3) = 1d0/8*(1-u)*(1+v)
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickNodalPBasis','Unknown node for brick')
      END SELECT
    END FUNCTION dBrickNodalPBasis


!------------------------------------------------------------------------------
!>     Brick edge basis at point (u,v,w).
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
    END FUNCTION BrickEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION dBrickEdgePBasis(edge, i , u, v, w, invertEdge) RESULT(grad)
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
    END FUNCTION dBrickEdgePBasis


!------------------------------------------------------------------------------
!>     Brick edge basis at point (u,v,w). Compatible with pyramidal edge basis.
!------------------------------------------------------------------------------
    FUNCTION BrickPyraEdgePBasis(edge, i , u, v, w, invertEdge) RESULT(value)
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
      REAL (KIND=dp) :: Na, Nb, La, Lb, phiPar, value
      INTEGER :: nodes(2)

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 12) THEN
         CALL Fatal('PElementBase::BrickPyraEdgePBasis','Unknown edge for brick')
      END IF

      ! Get nodes of edge
      nodes(1:2) = getBrickEdgeMap(edge)
      ! Bilinear nodal functions
      Na = BrickNodalPBasis(nodes(1),u,v,w)
      Nb = BrickNodalPBasis(nodes(2),u,v,w)
      ! Affine functions for edge direction
      La = BrickL(nodes(1),u,v,w)
      Lb = BrickL(nodes(2),u,v,w)

      ! For inverted edges swap direction
      IF (invert) THEN
         phiPar = La-Lb
      ELSE
         phiPar = Lb-La
      END IF

      ! Get value of edge function
      value = Na*Nb*varPhi(i,phiPar)
    END FUNCTION BrickPyraEdgePBasis


!------------------------------------------------------------------------------
!>     Gradient of brick edge basis at point (u,v,w).
!------------------------------------------------------------------------------
    FUNCTION dBrickPyraEdgePBasis(edge, i, u, v, w, invertEdge) RESULT(grad)
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
      REAL (KIND=dp) :: Na,Nb,La,Lb, vPhi
      REAL (KIND=dp), DIMENSION(3) :: dNa, dNb, dLa, dLb, grad
      INTEGER :: nodes(2), tmp

      ! By default do not invert edges
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge      

      ! Parameter validity check
      IF (edge < 1 .OR. edge > 12) THEN
         CALL Fatal('PElementBase::dBrickPyraEdgePBasis','Unknown edge for brick')
      END IF

      ! Get nodes of edge
      nodes(1:2) = getBrickEdgeMap(edge)
      ! Bilinear nodal functions and their derivatives
      Na = BrickNodalPBasis(nodes(1),u,v,w)
      Nb = BrickNodalPBasis(nodes(2),u,v,w)
      dNa = dBrickNodalPBasis(nodes(1),u,v,w)
      dNb = dBrickNodalPBasis(nodes(2),u,v,w)

      ! For inverted edges swap direction
      IF (invert) THEN
         tmp=nodes(1)
         nodes(1) = nodes(2)
         nodes(2) = tmp
      END IF

      ! Affine functions and their derivatives for edge direction
      La = BrickL(nodes(1),u,v,w)
      Lb = BrickL(nodes(2),u,v,w)
      dLa = dBrickL(nodes(1),u,v,w)
      dLb = dBrickL(nodes(2),u,v,w)

      vPhi = varPhi(i,Lb-La)

      ! Get value of edge function
      grad = dNa*Nb*vPhi + Na*dNb*vPhi + Na*Nb*dVarPhi(i,Lb-La)*(dLb-dLa)
    END FUNCTION dBrickPyraEdgePBasis

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
    END FUNCTION dBrickFacePBasis


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

      value = Phi(i,u)*Phi(j,v)*Phi(k,w)
    END FUNCTION BrickBubblePBasis


!------------------------------------------------------------------------------
!>    Gradient of brick bubble basis at point (u,v,w)
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
    END FUNCTION dBrickBubblePBasis

    FUNCTION BrickL(which, u, v, w) RESULT(value)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v,w 
      ! Variables
      REAL(KIND=dp) :: value
      
      value = 0
      SELECT CASE(which)
      CASE (1)
         value = (3d0-u-v-w)/2d0
      CASE (2)
         value = (3d0+u-v-w)/2d0
      CASE (3)
         value = (3d0+u+v-w)/2d0
      CASE (4)
         value = (3d0-u+v-w)/2d0
      CASE (5)
         value = (3d0-u-v+w)/2d0
      CASE (6)
         value = (3d0+u-v+w)/2d0
      CASE (7)
         value = (3d0+u+v+w)/2d0
      CASE (8)   
         value = (3d0-u+v+w)/2d0
      CASE DEFAULT
         CALL Fatal('PElementBase::BrickL','Unknown function L for brick')
      END SELECT
    END FUNCTION BrickL

    FUNCTION dBrickL(which, u, v, w) RESULT(grad)
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
         grad(1) = 1d0/2 
         grad(2) = -1d0/2
         grad(3) = -1d0/2
      CASE (3)
         grad(1) = 1d0/2 
         grad(2) = 1d0/2
         grad(3) = -1d0/2
      CASE (4)
         grad(1) = -1d0/2 
         grad(2) = 1d0/2
         grad(3) = -1d0/2
      CASE (5)
         grad(1) = -1d0/2 
         grad(2) = -1d0/2
         grad(3) = 1d0/2
      CASE (6)
         grad(1) = 1d0/2 
         grad(2) = -1d0/2
         grad(3) = 1d0/2
      CASE (7)
         grad(1) = 1d0/2 
         grad(2) = 1d0/2
         grad(3) = 1d0/2
      CASE (8)   
         grad(1) = -1d0/2 
         grad(2) = 1d0/2
         grad(3) = 1d0/2
      CASE DEFAULT
         CALL Fatal('PElementBase::dBrickL','Unknown function dL for brick')
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
         dLb_La(1) = 1
         dLb_La(2) = 0
         dLb_La(3) = 0
      CASE(2)
         ! Choose correct edge function by type
         SELECT CASE(t)
         ! Type 1 tetrahedron, edge 4:(2->3)
         CASE (1)
            La = TetraNodalPBasis(2,u,v,w)
            Lb = TetraNodalPBasis(3,u,v,w)
            dLa = dTetraNodalPBasis(2,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w) 
            dLb_La(1) = -1d0/2
            dLb_La(2) = SQRT(3d0)/2
            dLb_La(3) = 0 
         ! Type 2 tetrahedron, edge 4:(3->2)
         CASE (2) 
            La = TetraNodalPBasis(3,u,v,w)
            Lb = TetraNodalPBasis(2,u,v,w)
            dLa = dTetraNodalPBasis(3,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w) 
            dLb_La(1) = 1d0/2
            dLb_La(2) = -SQRT(3d0)/2
            dLb_La(3) = 0
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraEdgePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE(3)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(3,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(3,u,v,w) 
         dLb_La(1) = 1d0/2
         dLb_La(2) = SQRT(3d0)/2 
         dLb_La(3) = 0
      CASE(4)
         La = TetraNodalPBasis(1,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
         dLb_La(1) = 1d0/2
         dLb_La(2) = SQRT(3d0)/6
         dLb_La(3) = SQRT(6d0)/3
      CASE(5)
         La = TetraNodalPBasis(2,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(2,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
         dLb_La(1) = -1d0/2
         dLb_La(2) = SQRT(3d0)/6
         dLb_La(3) = SQRT(6d0)/3
      CASE(6)
         La = TetraNodalPBasis(3,u,v,w)
         Lb = TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(3,u,v,w)
         dLb = dTetraNodalPBasis(4,u,v,w) 
         dLb_La(1) = 0
         dLb_La(2) = -SQRT(3d0)/3
         dLb_La(3) = SQRT(6d0)/3
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraEdgePBasis','Unknown edge for tetrahedron')
      END SELECT

      ! Calculate gradient from given parameters
      ! General form for tetra edge gradients is 
      ! 
      ! Grad(Le) = 
      ! Grad(La)*Lb*varPhi(i,Lb-La) + La*Grad(Lb)*varPhi(i,Lb-La) + 
      ! La*Lb*dVarPhi(i,Lb-La)*Grad(Lb-La) 

      vPhi = varPhi(i, Lb-La)
      grad = dLa*Lb*vPhi+La*dLb*vPhi+La*Lb*dVarPhi(i,Lb-La)*dLb_La
    END FUNCTION dTetraEdgePBasis
      

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
         L1=TetraNodalPBasis(1,u,v,w)
         L2=TetraNodalPBasis(2,u,v,w)
         L3=TetraNodalPBasis(3,u,v,w)
            
         SELECT CASE(t)   
         CASE (1)
            value = L1*L2*L3*LegendreP(i,L2-L1)*LegendreP(j,2*L3-1)
         CASE (2)
            value = L1*L2*L3*LegendreP(i,L3-L1)*LegendreP(j,2*L2-1)
         CASE DEFAULT
            CALL Fatal('PElementBase::TetraFacePBasis','Unknown type for tetrahedron')
         END SELECT
      CASE (2)
         L1=TetraNodalPBasis(1,u,v,w)
         L2=TetraNodalPBasis(2,u,v,w)
         L4=TetraNodalPBasis(4,u,v,w)
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

      grad = 0
      SELECT CASE(face)
      CASE (1)
         SELECT CASE (t)
         CASE (1)
            La=TetraNodalPBasis(1,u,v,w)
            Lb=TetraNodalPBasis(2,u,v,w)
            Lc=TetraNodalPBasis(3,u,v,w)
            dLa = dTetraNodalPbasis(1,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w)
            dLc = dTetraNodalPBasis(3,u,v,w)
            ! Set up derivative of first inner function
            dLb_La(1) =  1d0 
            dLb_La(2) =  0d0
            dLb_La(3) =  0d0
            ! Derivative of second inner function
            dLc_1(1) =  0d0
            dLc_1(2) =  2d0*SQRT(3d0)/3
            dLc_1(3) =  -SQRT(6d0)/6
         CASE (2)
            La=TetraNodalPBasis(1,u,v,w)
            Lb=TetraNodalPBasis(3,u,v,w)
            Lc=TetraNodalPBasis(2,u,v,w)
            dLa = dTetraNodalPbasis(1,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w)
            dLc = dTetraNodalPBasis(2,u,v,w)
            ! Set up derivative of first inner function
            dLb_La(1) = 1d0/2 
            dLb_La(2) = SQRT(3d0)/2 
            dLb_La(3) = 0 
            ! Derivative of second inner function
            dLc_1(1) = 1d0 
            dLc_1(2) = -SQRT(3d0)/3 
            dLc_1(3) = -SQRT(6d0)/6 
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraFacePBasis','Unknown type for tetrahedron')       
         END SELECT
      CASE (2)
         La=TetraNodalPBasis(1,u,v,w)
         Lb=TetraNodalPBasis(2,u,v,w)
         Lc=TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(2,u,v,w)
         dLc = dTetraNodalPBasis(4,u,v,w)
         dLb_La(1) =  1d0 
         dLb_La(2) =  0d0
         dLb_La(3) =  0d0
         dLc_1(1) =  0d0
         dLc_1(2) =  0d0
         dLc_1(3) =  SQRT(6d0)/2
      CASE (3)
         SELECT CASE(t)
         ! Type 1 tetrahedron: Face 4:(2,3,4)
         CASE (1)
            La=TetraNodalPBasis(2,u,v,w)
            Lb=TetraNodalPBasis(3,u,v,w)
            Lc=TetraNodalPBasis(4,u,v,w)
            dLa = dTetraNodalPBasis(2,u,v,w)
            dLb = dTetraNodalPBasis(3,u,v,w)
            dLc = dTetraNodalPBasis(4,u,v,w)
            dLb_La(1) =  -1d0/2 
            dLb_La(2) =  SQRT(3d0)/2
            dLb_La(3) =  0d0
         ! Type 2 tetrahedron: Face 4:(3,2,4)
         CASE (2)
            La=TetraNodalPBasis(3,u,v,w)
            Lb=TetraNodalPBasis(2,u,v,w)
            Lc=TetraNodalPBasis(4,u,v,w)
            dLa = dTetraNodalPBasis(3,u,v,w)
            dLb = dTetraNodalPBasis(2,u,v,w)
            dLc = dTetraNodalPBasis(4,u,v,w)
            dLb_La(1) =  1d0/2 
            dLb_La(2) =  -SQRT(3d0)/2
            dLb_La(3) =  0d0
         CASE DEFAULT
            CALL Fatal('PElementBase::dTetraFacePBasis','Unknown type for tetrahedron')
         END SELECT
         ! Derivative of second inner function is equal for both types
         dLc_1(1) =  0d0
         dLc_1(2) =  0d0
         dLc_1(3) =  SQRT(6d0)/2
      CASE (4)
         La=TetraNodalPBasis(1,u,v,w)
         Lb=TetraNodalPBasis(3,u,v,w)
         Lc=TetraNodalPBasis(4,u,v,w)
         dLa = dTetraNodalPBasis(1,u,v,w)
         dLb = dTetraNodalPBasis(3,u,v,w)
         dLc = dTetraNodalPBasis(4,u,v,w)
         dLb_La(1) =  1d0/2 
         dLb_La(2) =  SQRT(3d0)/2
         dLb_La(3) =  0d0
         dLc_1(1) =  0d0
         dLc_1(2) =  0d0
         dLc_1(3) =  SQRT(6d0)/2
      CASE DEFAULT 
         CALL Fatal('PElementBase::dTetraFacePBasis','Unknown face for tetrahedron')
      END SELECT
      
      Legi = LegendreP(i, Lb-La)
      Legj = LegendreP(j, 2*Lc-1)

      ! Calculate gradient from given parameters 
      grad = dLa*Lb*Lc*Legi*Legj + La*dLb*Lc*Legi*Legj + La*Lb*dLc*Legi*Legj + &
           La*Lb*Lc*dLegendreP(i,Lb-La)*dLb_La*Legj + La*Lb*Lc*Legi*dLegendreP(j,2*Lc-1) * dLc_1
    
    END FUNCTION dTetraFacePBasis


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

      value = 0
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
      REAL (KIND=dp) :: L1, L2, L3, L4, L2_L1, L3_1, L4_1, a, b, c 
      REAL (KIND=dp), DIMENSION(3) :: grad

      grad = 0
      L1 = TetraNodalPBasis(1,u,v,w)
      L2 = TetraNodalPBasis(2,u,v,w)
      L3 = TetraNodalPBasis(3,u,v,w)
      L4 = TetraNodalPBasis(4,u,v,w)
      L2_L1 = L2 - L1
      L3_1 = 2*L3 - 1
      L4_1 = 2*L4 - 1
      a = LegendreP(i,L2_L1)
      b = LegendreP(j,L3_1)
      c = LegendreP(k,L4_1)

      ! Gradients of tetrahedral bubble basis functions 
      grad(1) = -1d0/2*L2*L3*L4*a*b*c + 1d0/2*L1*L3*L4*a*b*c + L1*L2*L3*L4*dLegendreP(i,L2_L1)*b*c
      grad(2) = -SQRT(3d0)/6*L2*L3*L4*a*b*c - SQRT(3d0)/6*L1*L3*L4*a*b*c + SQRT(3d0)/3*L1*L2*L4*a*b*c &
           + 2*SQRT(3d0)/3*L1*L2*L3*L4*a*dLegendreP(j,L3_1)*c
      grad(3) = -SQRT(6d0)/12*L2*L3*L4*a*b*c - SQRT(6d0)/12*L1*L3*L4*a*b*c - SQRT(6d0)/12*L1*L2*L4*a*b*c &
           + SQRT(6d0)/4*L1*L2*L3*a*b*c - SQRT(6d0)/6*L1*L2*L3*L4*a*dLegendreP(j,L3_1)*c &
           + SQRT(6d0)/2*L1*L2*L3*L4*a*b*dLegendreP(k,L4_1)
    END FUNCTION dTetraBubblePBasis


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
         value = 1d0/2*WedgeL(1,u,v)*(1d0-w)
      CASE (2)
         value = 1d0/2*WedgeL(2,u,v)*(1d0-w)
      CASE (3)
         value = 1d0/2*WedgeL(3,u,v)*(1-w)
      CASE (4)
         value = 1d0/2*WedgeL(1,u,v)*(1d0+w)
      CASE (5)   
         value = 1d0/2*WedgeL(2,u,v)*(1d0+w)
      CASE (6)
         value = 1d0/2*WedgeL(3,u,v)*(1+w)
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeNodalPBasis','Unknown node for wedge')
      END SELECT
    END FUNCTION WedgeNodalPBasis


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
      grad = 1d0/2*dLa*Lb*phiI*(1d0+parW) + 1d0/2*La*dLb*phiI*(1d0+parW) + &
           1d0/2*La*Lb*dVarPhi(i,Lb-La)*(dLb-dLa)*(1+parW) + 1d0/2*La*Lb*phiI*dW
    END FUNCTION dWedgeEdgePBasis


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
    END FUNCTION dWedgeFacePBasis


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
      REAL(KIND=dp) :: L1,L2,L3,value
      
      L1 = WedgeL(1,u,v)
      L2 = WedgeL(2,u,v)
      L3 = WedgeL(3,u,v)

      ! Get value of bubble function
      value = L1*L2*L3*LegendreP(i,L2-L1)*LegendreP(j,2d0*L3-1)*Phi(k,w)
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
    END FUNCTION dWedgeBubblePBasis

    FUNCTION WedgeL(which, u, v) RESULT(value)
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
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeL','Unknown linear function L for wedge')
      END SELECT
    END FUNCTION WedgeL

    FUNCTION WedgeH(which, w) RESULT(value)
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
      CASE DEFAULT
         CALL Fatal('PElementBase::WedgeH','Unknown linear function H for wedge')
      END SELECT
    END FUNCTION WedgeH

    FUNCTION dWedgeL(which, u, v) RESULT(grad)
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
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeL','Unknown linear function dL for wedge')
      END SELECT
    END FUNCTION dWedgeL

    FUNCTION dWedgeH(which, w) RESULT(grad)
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
      CASE DEFAULT
         CALL Fatal('PElementBase::dWedgeH','Unknown linear function dH for wedge')
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
      REAL(KIND=dp) :: Ta, Tb, value

      SELECT CASE(node)
      CASE (1)
         Ta = PyramidT(0,u,w)
         Tb = PyramidT(0,v,w)
         value = Ta*Tb*(1-w/SQRT(2d0))
      CASE (2)
         Ta = PyramidT(1,u,w)
         Tb = PyramidT(0,v,w)
         value = Ta*Tb*(1-w/SQRT(2d0))
      CASE (3)
         Ta = PyramidT(1,u,w)
         Tb = PyramidT(1,v,w)
         value = Ta*Tb*(1-w/SQRT(2d0))
      CASE (4)
         Ta = PyramidT(0,u,w)
         Tb = PyramidT(1,v,w)
         value = Ta*Tb*(1-w/SQRT(2d0))
      CASE (5)
         value = w/SQRT(2d0)
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
      REAL(KIND=dp) :: Ta, Tb, Tc, Temp(2)
      REAL(KIND=dp), DIMENSION(3) :: dTa, dTb, dTc, grad

      grad = 0
      dTa = 0
      dTb = 0
      SELECT CASE(node)
      CASE (1)
         Ta = PyramidT(0,u,w)
         Tb = PyramidT(0,v,w)
         temp = dPyramidT(0,u,w)
         dTa(1) = temp(1)
         dTa(3) = temp(2)
         dTb(2:3) = dPyramidT(0,v,w)
      CASE (2)
         Ta = PyramidT(1,u,w)
         Tb = PyramidT(0,v,w)
         temp = dPyramidT(1,u,w)
         dTa(1) = temp(1)
         dTa(3) = temp(2)
         dTb(2:3) = dPyramidT(0,v,w)
      CASE (3)
         Ta = PyramidT(1,u,w)
         Tb = PyramidT(1,v,w)
         temp = dPyramidT(1,u,w)
         dTa(1) = temp(1)
         dTa(3) = temp(2)
         dTb(2:3) = dPyramidT(1,v,w)
      CASE (4)
         Ta = PyramidT(0,u,w)
         Tb = PyramidT(1,v,w)
         temp = dPyramidT(0,u,w)
         dTa(1) = temp(1)
         dTa(3) = temp(2)
         dTb(2:3) = dPyramidT(1,v,w)
      CASE (5)
         ! Calculate value of gradient and return
         grad(3) = 1d0/SQRT(2d0)
         RETURN
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidNodalPBasis','Unknown node for pyramid')
      END SELECT

      ! Set parameter Tc
      Tc = (1-w/SQRT(2d0))
      dTc = 0
      dTc(3) = -1d0/SQRT(2d0) 

      ! Calculate value of gradient from general form
      grad = dTa*Tb*Tc + Ta*dTb*Tc + Ta*Tb*dTc
    END FUNCTION dPyramidNodalPBasis


	
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
      REAL(KIND=dp) :: Pa, Pb, phiPar, value
      LOGICAL :: invert
            
      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      ! Set parameters
      value = 0
      SELECT CASE(edge)
      CASE (1)
         phiPar = u
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(2,u,v,w)
      CASE (2)
         phiPar = v
         Pa = PyramidNodalPBasis(2,u,v,w)
         Pb = PyramidNodalPBasis(3,u,v,w)
      CASE (3)
         phiPar = u
         Pa = PyramidNodalPBasis(4,u,v,w)
         Pb = PyramidNodalPBasis(3,u,v,w)
      CASE (4)
         phiPar = v
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(4,u,v,w)
      CASE (5)
         phiPar = u/2+v/2+w/SQRT(2d0)
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
      CASE (6)
         phiPar = -u/2+v/2+w/SQRT(2d0)
         Pa = PyramidNodalPBasis(2,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
      CASE (7)
         phiPar = -u/2-v/2+w/SQRT(2d0)
         Pa = PyramidNodalPBasis(3,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
      CASE (8)
         phiPar = u/2-v/2+w/SQRT(2d0)
         Pa = PyramidNodalPBasis(4,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidEdgePBasis','Unknown edge for pyramid')
      END SELECT

      ! Invert edge if needed
      IF (invert) THEN
         phiPar = -phiPar
      END IF

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
      REAL(KIND=dp) :: Pa, Pb, phiPar, vPhiI
      REAL(KIND=dp), DIMENSION(3) :: dPa, dPb, dPhiPar, grad
      LOGICAL :: invert
            
      ! Edge is not inverted by default
      invert = .FALSE.
      IF (PRESENT(invertEdge)) invert = invertEdge

      ! Set parameters
      grad = 0
      dPhiPar = 0
      SELECT CASE(edge)
      CASE (1)
         phiPar = u
         dPhiPar(1) = 1 
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(2,u,v,w)
         dPa = dPyramidNodalPBasis(1,u,v,w)
         dPb = dPyramidNodalPBasis(2,u,v,w)
      CASE (2)
         phiPar = v
         dPhiPar(2) = 1
         Pa = PyramidNodalPBasis(2,u,v,w)
         Pb = PyramidNodalPBasis(3,u,v,w)
         dPa = dPyramidNodalPBasis(2,u,v,w)
         dPb = dPyramidNodalPBasis(3,u,v,w)
      CASE (3)
         phiPar = u
         dPhiPar(1) = 1
         Pa = PyramidNodalPBasis(4,u,v,w)
         Pb = PyramidNodalPBasis(3,u,v,w)
         dPa = dPyramidNodalPBasis(4,u,v,w)
         dPb = dPyramidNodalPBasis(3,u,v,w)
      CASE (4)
         phiPar = v
         dPhiPar(2) = 1
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(4,u,v,w)
         dPa = dPyramidNodalPBasis(1,u,v,w)
         dPb = dPyramidNodalPBasis(4,u,v,w)
      CASE (5)
         phiPar = u/2+v/2+w/SQRT(2d0)
         dPhiPar(1) = 1d0/2
         dPhiPar(2) = 1d0/2
         dPhiPar(3) = SQRT(2d0)/2
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
         dPa = dPyramidNodalPBasis(1,u,v,w)
         dPb = dPyramidNodalPBasis(5,u,v,w)
      CASE (6)
         phiPar = -u/2+v/2+w/SQRT(2d0)
         dPhiPar(1) = -1d0/2
         dPhiPar(2) = 1d0/2
         dPhiPar(3) = SQRT(2d0)/2
         Pa = PyramidNodalPBasis(2,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
         dPa = dPyramidNodalPBasis(2,u,v,w)
         dPb = dPyramidNodalPBasis(5,u,v,w)
      CASE (7)
         phiPar = -u/2-v/2+w/SQRT(2d0)
         dPhiPar(1) = -1d0/2
         dPhiPar(2) = -1d0/2
         dPhiPar(3) = SQRT(2d0)/2
         Pa = PyramidNodalPBasis(3,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
         dPa = dPyramidNodalPBasis(3,u,v,w)
         dPb = dPyramidNodalPBasis(5,u,v,w)
      CASE (8)
         phiPar = u/2-v/2+w/SQRT(2d0)
         dPhiPar(1) = 1d0/2
         dPhiPar(2) = -1d0/2
         dPhiPar(3) = SQRT(2d0)/2
         Pa = PyramidNodalPBasis(4,u,v,w)
         Pb = PyramidNodalPBasis(5,u,v,w)
         dPa = dPyramidNodalPBasis(4,u,v,w)
         dPb = dPyramidNodalPBasis(5,u,v,w)
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidEdgePBasis','Unknown edge for pyramid')
      END SELECT

      ! Invert edge if needed
      IF (invert) THEN
         phiPar = -phiPar
         dPhiPar = -dPhiPar
      END IF

      vPhiI = varPhi(i,phiPar)

      grad = dPa*Pb*vPhiI + Pa*dPb*vPhiI + Pa*Pb*dVarPhi(i,phiPar)*dPhiPar
    END FUNCTION dPyramidEdgePBasis

	
!-----------------------------------------------------------------------------	
!>     Pyramid face basis at point (u,v,w)
!-----------------------------------------------------------------------------	
    FUNCTION PyramidFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(value)
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
      REAL(KIND=dp) :: Pa, Pb, Pc, La, Lb, Lc, value
      INTEGER :: local(4)
      
      ! If local numbering not present use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getPyramidFaceMap(face)
      ELSE
         local(1:4) = localNumbers
      END IF
      
      SELECT CASE(face)
      CASE (1)
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(3,u,v,w)
         La = PyramidL(local(1),u,v)
         Lb = PyramidL(local(2),u,v)
         Lc = PyramidL(local(4),u,v)
         value = Pa*Pb*varPhi(i,Lb-La)*varPhi(j,Lc-La)
      CASE (2,3,4,5)
         Pa = PyramidNodalPBasis(local(1),u,v,w)
         Pb = PyramidNodalPBasis(local(2),u,v,w)
         Pc = PyramidNodalPBasis(local(3),u,v,w)
         value = Pa*Pb*Pc*LegendreP(i,Pb-Pa)*LegendreP(j,2*Pc-1)
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidFacePBasis','Unknown face for pyramid')
      END SELECT
    END FUNCTION PyramidFacePBasis

!-----------------------------------------------------------------------------	
!>     Gradient of pyramid face basis at point (u,v,w)
!-----------------------------------------------------------------------------		
    FUNCTION dPyramidFacePBasis(face, i, j, u, v, w, localNumbers) RESULT(grad)
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
      REAL(KIND=dp) :: Pa, Pb, Pc, La, Lb, Lc, legI, legJ
      REAL(KIND=dp), DIMENSION(3) :: dPa, dPb, dPc, dLa, dLb, dLc, grad
      INTEGER :: local(4)

      ! If local numbering not present use default numbers
      IF (.NOT. PRESENT(localNumbers)) THEN
         local(1:4) = getPyramidFaceMap(face)
      ELSE
         local(1:4) = localNumbers
      END IF

      ! Calculate value of function by face
      SELECT CASE (face)
      ! Square face
      CASE(1)
         Pa = PyramidNodalPBasis(1,u,v,w)
         Pb = PyramidNodalPBasis(3,u,v,w)
         dPa = dPyramidNodalPBasis(1,u,v,w)
         dPb = dPyramidNodalPBasis(3,u,v,w)
         La = PyramidL(local(1),u,v)
         Lb = PyramidL(local(2),u,v)
         Lc = PyramidL(local(4),u,v)
         dLa = dPyramidL(local(1),u,v)
         dLb = dPyramidL(local(2),u,v)
         dLc = dPyramidL(local(4),u,v)
         legI = varPhi(i,Lb-La)
         legJ = varPhi(j,Lc-La)
         grad = dPa*Pb*legI*legJ + Pa*dPb*legI*legJ +&
              Pa*Pb*dVarPhi(i,Lb-La)*(dLb-dLa)*legJ + &
              Pa*Pb*legI*dVarPhi(j,Lc-La)*(dLc-dLa)
      ! Triangle face
      CASE (2,3,4,5)
         Pa = PyramidNodalPBasis(local(1),u,v,w)
         Pb = PyramidNodalPBasis(local(2),u,v,w)
         Pc = PyramidNodalPBasis(local(3),u,v,w)
         dPa = dPyramidNodalPBasis(local(1),u,v,w)
         dPb = dPyramidNodalPBasis(local(2),u,v,w)
         dPc = dPyramidNodalPBasis(local(3),u,v,w)
         legI = LegendreP(i,Pb-Pa)
         legJ = LegendreP(j,2*Pc-1)
         grad = dPa*Pb*Pc*legI*legJ + Pa*dPb*Pc*legI*legJ + &
              Pa*Pb*dPc*legI*legJ + &
              Pa*Pb*Pc*dLegendreP(i,Pb-Pa)*(dPb-dPa)*legJ + &
              Pa*Pb*Pc*legI*dLegendreP(j,2*Pc-1)*(2*dPc)
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidFacePBasis','Unknown face for pyramid')
      END SELECT
    END FUNCTION dPyramidFacePBasis

	
	
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
      REAL(KIND=dp) :: value

      ! Calculate value of function
      value = PyramidNodalPBasis(1,u,v,w)*PyramidNodalPBasis(3,u,v,w)* &
           PyramidNodalPBasis(5,u,v,w)*LegendreP(i,u/(1-w/SQRT(2d0)))* &
           LegendreP(j,v/(1-w/SQRT(2d0)))*LegendreP(k,w/SQRT(2d0))
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
      REAL(KIND=dp) :: P1, P3, P5, legI, legJ, legK
      REAL(KIND=dp), DIMENSION(3) :: dP1, dP3, dP5, dLegIPar, dLegJPar, &
           dLegKPar, grad

      P1 = PyramidNodalPBasis(1,u,v,w)
      P3 = PyramidNodalPBasis(3,u,v,w)
      P5 = PyramidNodalPBasis(5,u,v,w)
      dP1 = dPyramidNodalPBasis(1,u,v,w)
      dP3 = dPyramidNodalPBasis(3,u,v,w)
      dP5 = dPyramidNodalPBasis(5,u,v,w)
      legI = LegendreP(i,u/(1-w/SQRT(2d0)))
      dLegIPar = [ 1d0/(1-w/SQRT(2d0)), 0d0, u*SQRT(2d0)/(2*(1-w/SQRT(2d0))**2) ]
      legJ = LegendreP(j,v/(1-w/SQRT(2d0)))
      dLegJPar = [ 0d0, 1d0/(1-w/SQRT(2d0)), v*SQRT(2d0)/(2*(1-w/SQRT(2d0))**2) ]
      legK = LegendreP(k,w/SQRT(2d0))
      dLegKPar = [ 0d0, 0d0, 1d0/SQRT(2d0) ]

      ! Calculate value of gradient
      grad = 0
      grad = dP1*P3*P5*legI*legJ*legK + P1*dP3*P5*legI*legJ*legK + &
           P1*P3*dP5*legI*legJ*legK + &
           P1*P3*P5*dLegendreP(i,u/(1-w/SQRT(2d0)))*dLegIPar*legJ*legK + &
           P1*P3*P5*legI*dLegendreP(j,v/(1-w/SQRT(2d0)))*dLegJPar*legK + &
           P1*P3*P5*legI*legJ*dLegendreP(k,w/SQRT(2d0))*dLegKPar
    END FUNCTION dPyramidBubblePBasis

    FUNCTION PyramidT(which, c, t) RESULT(value)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: c,t
      ! Variables
      REAL(KIND=dp) :: value

      SELECT CASE(which)
      CASE (0)
         value = ((1-t/SQRT(2d0))-c)/(2*(1-t/SQRT(2d0)))
      CASE (1)
         value = ((1-t/SQRT(2d0))+c)/(2*(1-t/SQRT(2d0)))
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidT','Unknown function T for pyramid')
      END SELECT
    END FUNCTION PyramidT

    FUNCTION dPyramidT(which, c, t) RESULT(grad)
      IMPLICIT NONE

      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: c,t
      ! Variables
      REAL(KIND=dp), DIMENSION(2) :: grad

      SELECT CASE(which)
      CASE (0)
         grad(1) = -1d0/(2-t*SQRT(2d0))
         grad(2) = -SQRT(2d0)/(2*(2-t*SQRT(2d0)))+ &
              (1-t*SQRT(2d0)/2-c)*SQRT(2d0)/(2-t*SQRT(2d0))**2
      CASE (1)
         grad(1) = 1d0/(2-t*SQRT(2d0))
         grad(2) = -SQRT(2d0)/(2*(2-t*SQRT(2d0)))+ &
              (1-t*SQRT(2d0)/2+c)*SQRT(2d0)/(2-t*SQRT(2d0))**2
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidT','Unknown function dT for pyramid')
      END SELECT
    END FUNCTION dPyramidT

    ! Define affine coordinates for pyramid square face

    FUNCTION PyramidL(which, u, v) RESULT(value)
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
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidL','Unknown affine coordinate for square face')
      END SELECT
    END FUNCTION PyramidL

    FUNCTION dPyramidL(which, u, v) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      REAL(KIND=dp), INTENT(IN) :: u,v
      ! Variables
      REAL(KIND=dp) :: grad(3)

      SELECT CASE (which)
      CASE (1)
         grad = [ -1d0/2,-1d0/2,0d0 ]
      CASE (2)
         grad = [ 1d0/2,-1d0/2,0d0 ]
      CASE (3)
         grad = [ 1d0/2,1d0/2,0d0 ]
      CASE (4)
         grad = [ -1d0/2,1d0/2,0d0 ]
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidL','Unknown affine coordinate for square face')
      END SELECT
    END FUNCTION dPyramidL


!------------------------------------------------------------------------------
!>    Phi function value at point x. Phi is defined as (Szabo & Babuska: Finite
!>    Element Analysis, p.38). 
!------------------------------------------------------------------------------
    FUNCTION Phi(i,x) RESULT(value)
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

      IF (i < 2) THEN
         CALL Fatal('PElementBase::Phi','Phi(i,x) not defined for i<2')
      END IF

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
    FUNCTION dPhi(i,x) RESULT(value)
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

      IF (i < 2) THEN
         CALL Fatal('PElementBase::dPhi','dPhi(i,x) not defined for i<2')
      END IF

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
      CASE DEFAULT 
         value = SQRT(1d0/(2*(2*i-1)))*(dLegendreP(i,x)-dLegendreP(i-2,x))
      END SELECT
    END FUNCTION dPhi

    
!------------------------------------------------------------------------------
!>    varPhi function value at point x. Phi is defined as (Szabo & Babuska: Finite
!>    Element Analysis, p. 103).     
!>    Phi(i,x)=1/4*(1-x^2)*varPhi(i,x)
!------------------------------------------------------------------------------
    FUNCTION varPhi(i,x) RESULT(value)
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
      CASE (:1)
         CALL Fatal('PElementBase::varPhi','varPhi not defined for i<2')
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
      CASE DEFAULT
         IF (x==1 .OR. x==-1) THEN
            ! TEMP SOLUTION!
            ! Try to interpolate value of function
            value = ((4*Phi(i,(x-dx))/(1-(x-dx)**2))+(4*Phi(i,(x+dx))/(1-(x+dx)**2)))/2
         ELSE 
            value = 4*Phi(i,x)/(1-x**2)
         END IF
      END SELECT
      
    END FUNCTION varPhi


!------------------------------------------------------------------------------
!>    Derivative of varPhi function at point x.
!------------------------------------------------------------------------------
    FUNCTION dVarPhi(i,x) RESULT(value)
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
      CASE (:1)
         CALL Fatal('PElementBase::dVarPhi','dVarPhi not defined for i<2')
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
      END SELECT
    END FUNCTION dVarPhi

    
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
    RECURSIVE FUNCTION LegendreP(l,x) RESULT(value) 
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
      CASE (:-1)
         CALL Fatal('PElementBase::LegendreP','LegendreP not defined for l < 0')
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
    RECURSIVE FUNCTION dLegendreP(l,x) RESULT(value)
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
      CASE (:-1)
         CALL Fatal('PElementBase::dLegendreP','dLegendreP not defined for l < 0')
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
      END SELECT
      ! Value now contains P,(l,x)
    END FUNCTION dLegendreP

    ! Function value = x^n

    FUNCTION toExp(x,n) RESULT(value)
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

END MODULE PElementBase

!> \}