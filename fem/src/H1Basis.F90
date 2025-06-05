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
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 31 May 2017
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!>  Module defining vectorized p element basis functions and mappings.
!-----------------------------------------------------------------------------

#include "../config.h"

MODULE H1Basis
  USE Messages
  USE Types, ONLY : dp, VECTOR_BLOCK_LENGTH
  
  ! Module contains vectorized version of FE basis
  ! functions for selected elements

  INTEGER, PARAMETER :: H1Basis_MaxPElementEdgeNodes=2, &
          H1Basis_MaxPElementEdges = 12, &
          H1Basis_MaxPElementFaceNodes = 4, &
          H1Basis_MaxPElementFaces = 6
CONTAINS

  SUBROUTINE H1Basis_GetEdgeDirection(ecode, nedges, globalind, direction)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ecode, nedges
    INTEGER, DIMENSION(:), POINTER CONTIG, INTENT(IN) :: globalind
    INTEGER, DIMENSION(H1Basis_MaxPElementEdgeNodes, &
                       H1Basis_MaxPElementEdges), INTENT(INOUT) :: direction
    INTEGER :: i, tmp

    CALL H1Basis_GetEdgeMap(ecode, direction)

!DIR$ LOOP COUNT MAX=12
    DO i=1,nedges
      ! Swap local indices if needed to set the global direction (increasing)
      IF (globalind(direction(2,i))<globalind(direction(1,i))) THEN
        tmp = direction(1,i)
        direction(1,i)=direction(2,i)
        direction(2,i)=tmp
      END IF
    END DO
  END SUBROUTINE H1Basis_GetEdgeDirection
  
  SUBROUTINE H1Basis_GetTetraEdgeDirection(ttype, direction)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ttype
    INTEGER, DIMENSION(H1Basis_MaxPElementEdgeNodes, &
                       H1Basis_MaxPElementEdges), TARGET, INTENT(INOUT) :: direction
    ! Local variables
    INTEGER, DIMENSION(:), POINTER CONTIG :: flatdir

    ! Remap array bounds
    flatdir(1:H1Basis_MaxPElementEdgeNodes*H1Basis_MaxPElementEdges) => direction
    SELECT CASE (ttype)
    CASE (1)
      flatdir(1:12) = [ 1,2, &
                        2,3, &
                        1,3, &
                        1,4, &
                        2,4, &
                        3,4 ]
    CASE (2)
      flatdir(1:12) = [ 1,2, &
                        3,2, &
                        1,3, &
                        1,4, &
                        2,4, &
                        3,4 ]      
    CASE DEFAULT
      CALL Fatal('H1Basis_GetTetraEdgeDirection','Unknown tetra type')
    END SELECT
  END SUBROUTINE H1Basis_GetTetraEdgeDirection
  
  SUBROUTINE H1Basis_GetFaceDirection(ecode, nfaces, globalind, direction)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ecode, nfaces
    INTEGER, DIMENSION(:), POINTER CONTIG, INTENT(IN) :: globalind
    INTEGER, DIMENSION(H1Basis_MaxPElementFaceNodes, &
                       H1Basis_MaxPElementFaces), INTENT(INOUT) :: direction
    INTEGER :: face, tmp, tmp1, tmp2, i, minI

    CALL H1Basis_GetFaceMap(ecode, direction)

!DIR$ LOOP COUNT MAX=6
    DO face=1,nfaces
      IF (direction(H1Basis_MaxPElementFaceNodes, face) == 0) THEN
        ! Triangle face (last local index 0)
        
        ! Global face direction is consistent when the local indices are chosen such 
        ! that the global index is sorted

        ! Sort globalind(direction(:,face))
        ! (1,2) -> (1,2) or (2,1)
        IF (globalind(direction(1,face))>globalind(direction(2,face))) THEN
          tmp = direction(1,face)
          direction(1,face)=direction(2,face)
          direction(2,face)=tmp
        END IF
        ! (1,3) -> (1,3) or (3,1)
        IF (globalind(direction(1,face))>globalind(direction(3,face))) THEN
          tmp = direction(1,face)
          direction(1,face)=direction(3,face)
          direction(3,face)=tmp
        END IF
        ! (2,3) -> (2,3) or (3,2)
        IF (globalind(direction(2,face))>globalind(direction(3,face))) THEN
          tmp = direction(2,face)
          direction(2,face)=direction(3,face)
          direction(3,face)=tmp
        END IF
      ELSE
        ! Square face
        ! Find index of minimum global index (A)
        minI = 1
        DO i=2,4
          IF (globalind(direction(i,face))<globalind(direction(minI,face))) minI=i
        END DO
        
        ! In place rotate face to left or right cyclically to move minI as first element
        ! sqface(1:4)=cshift(direction(1:4,face),minI-1)
        SELECT CASE(minI-1)
        CASE(0)
          ! Nothing to do!
        CASE(1)
          tmp = direction(1,face)
          ! arr(i) -> arr(i+1)
          DO i=1,3
            direction(i,face)=direction(i+1,face)
          END DO
          direction(4,face)=tmp
        CASE(2)
          ! (1,3)->(3,1)
          tmp1 = direction(1,face)
          direction(1,face)=direction(3,face)
          direction(3,face)=tmp1
          ! (2,4)->(4,2)
          tmp2 = direction(2,face)
          direction(2,face)=direction(4,face)
          direction(4,face)=tmp2
        CASE(3)
          tmp = direction(4,face)
          ! arr(i+1) -> arr(i)
          DO i=3,1,-1
            direction(i+1,face)=direction(i,face)
          END DO
          direction(1,face)=tmp
        END SELECT

        ! Find second smallest index B next to global minimum (either 2 or 4)
        IF (globalind(direction(2,face))>globalind(direction(4,face))) THEN
          tmp=direction(2,face)
          direction(2,face)=direction(4,face)
          direction(4,face)=tmp
        END IF

        ! Let C be the remaining index next to the global minimum A and 
        ! D the index opposite of A -> [A,B,D,C] forms a globally consistent ordering
      END IF
    END DO
  END SUBROUTINE H1Basis_GetFaceDirection

  SUBROUTINE H1Basis_GetTetraFaceDirection(ttype, direction)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ttype
    INTEGER, DIMENSION(H1Basis_MaxPElementFaceNodes, &
                       H1Basis_MaxPElementFaces), TARGET, INTENT(INOUT) :: direction
    ! Local variables
    INTEGER, DIMENSION(:), POINTER CONTIG :: flatdir

    ! Remap array bounds
    flatdir(1:H1Basis_MaxPElementFaceNodes*H1Basis_MaxPElementFaces) => direction
    SELECT CASE (ttype)
    CASE (1)
      flatdir(1:16) = [ 1,2,3,0, & 
                        1,2,4,0, &
                        2,3,4,0, &
                        1,3,4,0 ]
    CASE (2)
      flatdir(1:16) = [ 1,3,2,0, &
                        1,2,4,0, &
                        3,2,4,0, &
                        1,3,4,0 ]
    CASE DEFAULT
      CALL Fatal('H1Basis_GetTetraFaceDirection','Unknown tetra type')
    END SELECT
  END SUBROUTINE H1Basis_GetTetraFaceDirection

  SUBROUTINE H1Basis_GetEdgeMap(ecode, map)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ecode
    INTEGER, DIMENSION(H1Basis_MaxPElementEdgeNodes, &
                       H1Basis_MaxPElementEdges), TARGET, INTENT(OUT) :: map
    ! Local variables
    INTEGER, DIMENSION(:), POINTER CONTIG :: flatmap

    ! Remap array bounds
    flatmap(1:H1Basis_MaxPElementEdgeNodes*H1Basis_MaxPElementEdges) => map

    SELECT CASE (ecode)
    CASE(202)
      ! Line edge mappings
      flatmap(1:2) = [ 1,2 ]
    CASE(303)
      ! Triangle edge mappings
      flatmap(1:6) = [ 1,2, &
                       2,3, &
                       3,1 ]
    CASE(404)
      ! Quad edge mappings
      flatmap(1:8) = [ 1,2,&
                       2,3,&
                       4,3,&
                       1,4 ]
    CASE(504)
      ! Tetra edge mapping (for completeness, not needed for enforcing parity)
      CALL H1Basis_GetTetraEdgeDirection(1, map)

    CASE(605)
      flatmap(1:16) = [ 1,2,&
                        2,3,&
                        4,3,&
                        1,4,&
                        1,5,&
                        2,5,&
                        3,5,&
                        4,5 ]

    CASE(706)
      flatmap(1:18) = [ 1,2,&
                        2,3,&
                        3,1,&
                        4,5,&
                        5,6,&
                        6,4,&
                        1,4,&
                        2,5,&
                        3,6 ]
    CASE(808)
      flatmap(1:24) = [ 1,2,&
                        2,3,&
                        4,3,&
                        1,4,&
                        5,6,&
                        6,7,&
                        8,7,&
                        5,8,&
                        1,5,&
                        2,6,&
                        3,7,&
                        4,8 ]

    CASE DEFAULT
      CALL Fatal('H1Basis_GetEdgeMap','Not fully implemented yet!')
    END SELECT
  END SUBROUTINE H1Basis_GetEdgeMap

  SUBROUTINE H1Basis_GetFaceMap(ecode, map)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ecode
    INTEGER, DIMENSION(H1Basis_MaxPElementFaceNodes, &
                       H1Basis_MaxPElementFaces), TARGET, INTENT(OUT) :: map
    ! Local variables
    INTEGER, DIMENSION(:), POINTER CONTIG :: flatmap

    ! Remap array bounds
    flatmap(1:H1Basis_MaxPElementFaceNodes*H1Basis_MaxPElementFaces) => map

    SELECT CASE (ecode)
    CASE(303)
      ! 2D boundary element in a 3D mesh
      flatmap(1:4) = [ 1,2,3,0 ]
    CASE(404)
      ! 2D boundary element in a 3D mesh
      flatmap(1:4) = [ 1,2,3,4 ]
    CASE(504)
      ! Tetra face mapping (for completeness, not needed for enforcing parity)
      CALL H1Basis_GetTetraFaceDirection(1, map)

    CASE(605)
      ! Pyramid face mappings
      flatmap(1:20) = [ 1,2,3,4, &
                        1,2,5,0, &
                        2,3,5,0, &
                        3,4,5,0, &
                        4,1,5,0 ]

    CASE(706)
      flatmap(1:20) = [ 1,2,3,0, &
                        4,5,6,0, &
                        1,2,5,4, &
                        2,3,6,5, &
                        3,1,4,6 ]
    CASE(808)
      flatmap(1:24) = [ 1,2,3,4, &
                        5,6,7,8, &
                        1,2,6,5, &
                        2,3,7,6, &
                        4,3,7,8, &
                        1,4,8,5 ]
    CASE DEFAULT
      CALL Fatal('H1Basis_GetFaceMap','Not fully implemented yet!')
    END SELECT
  END SUBROUTINE H1Basis_GetFaceMap
  
  SUBROUTINE H1Basis_LineNodal(nvec, u, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    REAL(Kind=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: j
!DIR$ ASSUME_ALIGNED u:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,nbasis+1) = c*(1-u(j))
      ! Node 2
      fval(j,nbasis+2) = c*(1+u(j))
    END DO
    nbasis = nbasis + 2
  END SUBROUTINE H1Basis_LineNodal

  SUBROUTINE H1Basis_dLineNodal(nvec, u, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    REAL(Kind=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: j
!DIR$ ASSUME_ALIGNED u:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) = -c
      grad(j,nbasis+2,1) =  c
    END DO
    nbasis = nbasis + 2
  END SUBROUTINE H1Basis_dLineNodal

  SUBROUTINE H1Basis_LineBubbleP(nvec, u, pmax, nbasismax, fval, nbasis, invertEdge)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    LOGICAL, OPTIONAL :: invertEdge

    ! Local variables
    LOGICAL :: invert
    INTEGER :: p, j
!DIR$ ASSUME_ALIGNED u:64, fval:64

    ! Check if line basis has been inverted (not by default)
    invert = .FALSE.
    IF (PRESENT( invertEdge )) invert = invertEdge

    IF (invert) THEN
      DO p=1,pmax-1
        !_ELMER_OMP_SIMD 
        DO j=1,nvec
           fval(j,nbasis+p) = H1Basis_Phi(p+1,-u(j))
        END DO
      END DO
    ELSE
      DO p=1,pmax-1
        !_ELMER_OMP_SIMD 
        DO j=1,nvec
          fval(j,nbasis+p) = H1Basis_Phi(p+1,u(j))
        END DO
      END DO
    END IF

    ! nbasis = nbasis+(pmax-2)+1
    nbasis = nbasis+pmax-1
  END SUBROUTINE H1Basis_LineBubbleP

  SUBROUTINE H1Basis_dLineBubbleP(nvec, u, pmax, nbasismax, grad, nbasis, invertEdge)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    LOGICAL, OPTIONAL :: invertEdge

    ! Local variables
    LOGICAL :: invert
    INTEGER :: p, j
!DIR$ ASSUME_ALIGNED u:64, grad:64

    ! Check if line basis has been inverted (not by default)
    invert = .FALSE.
    IF (PRESENT( invertEdge )) invert = invertEdge

    IF (invert) THEN
      DO p=1,pmax-1
        ! First coordinate (xi)
        !_ELMER_OMP_SIMD
        DO j=1,nvec
          grad(j,nbasis+p,1) = H1Basis_dPhi(p+1,-u(j))
        END DO
      END DO
    ELSE
      DO p=1,pmax-1
        ! First coordinate (xi)
        !_ELMER_OMP_SIMD
        DO j=1,nvec
          grad(j,nbasis+p,1) = H1Basis_dPhi(p+1,u(j))
        END DO
      END DO
    END IF

    ! nbasis = nbasis+(pmax-2)+1
    nbasis = nbasis+pmax-1
  END SUBROUTINE H1Basis_dLineBubbleP

  SUBROUTINE H1Basis_TriangleNodalP(nvec, u, v, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    
    INTEGER :: j
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0)
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      fval(j,nbasis+1) = c*(1-u(j)-d*v(j))
      fval(j,nbasis+2) = c*(1+u(j)-d*v(j))
      fval(j,nbasis+3) = d*v(j)
    END DO
    nbasis = nbasis + 3
  END SUBROUTINE H1Basis_TriangleNodalP

  FUNCTION H1Basis_TriangleL(node, u, v) RESULT(fval)
    IMPLICIT NONE

    ! Parameters 
    INTEGER, INTENT(IN) :: node
    REAL(KIND=dp), INTENT(IN) :: u,v
    ! Result
    REAL(KIND=dp) :: fval
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) &
    !_ELMER_OMP _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1)
      fval = c*(1-u-d*v)
    CASE(2)
      fval = c*(1+u-d*v)
    CASE(3)
      fval = d*v
    END SELECT
  END FUNCTION H1Basis_TriangleL

  SUBROUTINE H1Basis_dTriangleNodalP(nvec, u, v, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER :: j
    REAL(KIND=dp), PARAMETER :: half=1._dp/2, a=1/SQRT(3._dp), b=-a/2
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) = -half
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,1) =  half
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,1) = 0
    END DO
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = b
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,2) = b
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,2) = a
    END DO
    nbasis = nbasis + 3
  END SUBROUTINE H1Basis_dTriangleNodalP
  
  FUNCTION H1Basis_dTriangleL(node) RESULT(grad)
    IMPLICIT NONE

    ! Parameters 
    INTEGER, INTENT(IN) :: node
    ! REAL(KIND=dp), INTENT(IN) :: u,v
    ! Result
    REAL(KIND=dp) :: grad(2)
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1)
      grad = c*[REAL(-1,dp), REAL(-d,dp)]
    CASE(2)
      grad = c*[REAL(1,dp), REAL(-d,dp)]
    CASE(3)
      grad = [REAL(0,dp), REAL(d,dp)]
    END SELECT
  END FUNCTION H1Basis_dTriangleL
  
  SUBROUTINE H1Basis_TriangleEdgeP(nvec, u, v, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    ! For each edge
    DO i=1,3
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb)
        DO k=1,nvec
          La = H1Basis_TriangleL(edgedir(1,i),u(k),v(k))
          Lb = H1Basis_TriangleL(edgedir(2,i),u(k),v(k))
          fval(k, nbasis+j) = La*Lb*H1Basis_varPhi(j+1,Lb-La)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO

  END SUBROUTINE H1Basis_TriangleEdgeP

  SUBROUTINE H1Basis_dTriangleEdgeP(nvec, u, v, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, vPhi, dVPhi, dLa(2), dLb(2)
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    ! For each edge
    DO i=1,3
      dLa = H1Basis_dTriangleL(edgedir(1,i))
      dLb = H1Basis_dTriangleL(edgedir(2,i))

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, vPhi, dVPhi)
        DO k=1,nvec
          La = H1Basis_TriangleL(edgedir(1,i),u(k),v(k))
          Lb = H1Basis_TriangleL(edgedir(2,i),u(k),v(k))
          vPhi = H1Basis_varPhi(j+1,Lb-La)
          dVPhi = H1Basis_dVarPhi(j+1,Lb-La)

          grad(k, nbasis+j, 1) = dLa(1)*Lb*vPhi+ &
                  La*dLb(1)*vPhi + La*Lb*dVPhi*(dLb(1)-dLa(1))
          grad(k, nbasis+j, 2) = dLa(2)*Lb*vPhi+ &
                  La*dLb(2)*vPhi + La*Lb*dVPhi*(dLb(2)-dLa(2))
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO

  END SUBROUTINE H1Basis_dTriangleEdgeP

  SUBROUTINE H1Basis_TriangleBubbleP(nvec, u, v, pmax, nbasismax, fval, nbasis, localnumbers)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, INTENT(IN), DIMENSION(3), OPTIONAL :: localnumbers
      
    ! Variables
    INTEGER :: i,j,k
    REAL (KIND=dp) :: La, Lb, Lc
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    IF (.NOT. PRESENT(localnumbers)) THEN
      ! Triangle bubble basis
      DO i = 0,pmax-3
        DO j = 1,pmax-i-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc)
          DO k=1,nvec
            La = H1Basis_TriangleL(1,u(k),v(k))
            Lb = H1Basis_TriangleL(2,u(k),v(k))
            Lc = H1Basis_TriangleL(3,u(k),v(k))
          
            fval(k,nbasis+j) = La*Lb*Lc*(H1Basis_PowInt((Lb-La),i))*&
                    H1Basis_PowInt((2*Lc-1),j-1)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-3) + 1
        nbasis = nbasis + MAX(pmax-i-2,0)
      END DO
    ELSE
      ! Triangular face basis
      DO i=0,pmax-3
        DO j=1,pmax-i-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc)
          DO k=1,nvec
            La=H1Basis_TriangleL(localnumbers(1),u(k),v(k))
            Lb=H1Basis_TriangleL(localnumbers(2),u(k),v(k))
            Lc=H1Basis_TriangleL(localnumbers(3),u(k),v(k))
            
            fval(k,nbasis+j) = La*Lb*Lc*H1Basis_LegendreP(i,Lb-La)* &
                                          H1Basis_LegendreP(j-1,2*Lc-1)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-3 + 1)
        nbasis = nbasis + MAX(pmax-i-2,0)
      END DO
    END IF
  END SUBROUTINE H1Basis_TriangleBubbleP

  SUBROUTINE H1Basis_dTriangleBubbleP(nvec, u, v, pmax, nbasismax, grad, nbasis, localnumbers)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, INTENT(IN), DIMENSION(3), OPTIONAL :: localnumbers

    ! Variables
    REAL (KIND=dp) :: La, Lb, Lc, Lc_1, Lb_Lai, Lc_1n, a, b, Lb_La,&
            dLa(2), dLb(2), dLc(2)
    INTEGER :: i,j,k
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = -SQRT(3D0)/6, e = SQRT(3D0)/3
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    IF (.NOT. PRESENT(localnumbers)) THEN
      ! Triangle bubble basis
      DO i = 0,pmax-3
        DO j = 1,pmax-i-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Lb_Lai, Lc_1n)
          DO k=1,nvec
            La = H1Basis_TriangleL(1,u(k),v(k))
            Lb = H1Basis_TriangleL(2,u(k),v(k))
            Lc = H1Basis_TriangleL(3,u(k),v(k))
            
            Lb_Lai = H1Basis_PowInt((Lb-La), i)
            Lc_1n = H1Basis_PowInt((2*Lc-1), j-1)
            
            ! Calculate value of function from general form
            grad(k,nbasis+j,1) = -c*Lb*Lc*Lb_Lai*Lc_1n + La*c*Lc*Lb_Lai*Lc_1n + &
                    La*Lb*Lc*i*(H1Basis_PowInt((Lb-La),i-1))*Lc_1n
            grad(k,nbasis+j,2) = d*Lb*Lc*Lb_Lai*Lc_1n + La*d*Lc*Lb_Lai*Lc_1n + &
                    La*Lb*e*Lb_Lai*Lc_1n + &
                    La*Lb*Lc*Lb_Lai*(j-1)*(H1Basis_PowInt((2*Lc-1),j-2))*2*e
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-3) + 1
        nbasis = nbasis + MAX(pmax-i-2,0)
      END DO
    ELSE
      ! Triangular face basis
      dLa = H1Basis_dTriangleL(localnumbers(1))
      dLb = H1Basis_dTriangleL(localnumbers(2))
      dLc = H1Basis_dTriangleL(localnumbers(3))

      DO i=0,pmax-3
        DO j=1,pmax-i-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Lb_La, Lc_1, a, b)
          DO k=1,nvec
            La=H1Basis_TriangleL(localnumbers(1),u(k),v(k))
            Lb=H1Basis_TriangleL(localnumbers(2),u(k),v(k))
            Lc=H1Basis_TriangleL(localnumbers(3),u(k),v(k))
            Lb_La=Lb-La
            Lc_1=2*Lc-1
            a=H1Basis_LegendreP(i,Lb_La)
            b=H1Basis_LegendreP(j-1,Lc_1)

            grad(k,nbasis+j,1) = dLa(1)*Lb*Lc*a*b + La*dLb(1)*Lc*a*b + &
                    La*Lb*dLc(1)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(i,Lb_La)*(dLb(1)-dLa(1))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(j-1,Lc_1)*2*dLc(1)
            grad(k,nbasis+j,2) = dLa(2)*Lb*Lc*a*b + La*dLb(2)*Lc*a*b + &
                    La*Lb*dLc(2)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(i,Lb_La)*(dLb(2)-dLa(2))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(j-1,Lc_1)*2*dLc(2)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-3 + 1)
        nbasis = nbasis + MAX(pmax-i-2,0)
      END DO
    END IF

  END SUBROUTINE H1Basis_dTriangleBubbleP

  FUNCTION H1Basis_PowInt(x,j) RESULT(powi)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: j
    
    INTEGER :: i
    REAL(KIND=dp) :: powi
    !_ELMER_OMP_DECLARE_SIMD _ELMER_LINEAR_REF(x) UNIFORM(j) NOTINBRANCH
    
    powi=REAL(1,dp)
    DO i=1,j
      powi=powi*x
    END DO

  END FUNCTION H1Basis_PowInt

  SUBROUTINE H1Basis_QuadNodal(nvec, u, v, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    REAL(Kind=dp), PARAMETER :: c = 1D0/4D0
    INTEGER :: j
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,1) = c*(1-u(j))*(1-v(j))
      ! Node 2
      fval(j,2) = c*(1+u(j))*(1-v(j))
      ! Node 3
      fval(j,3) = c*(1+u(j))*(1+v(j))
      ! Node 4
      fval(j,4) = c*(1-u(j))*(1+v(j))
    END DO
    nbasis = nbasis + 4
  END SUBROUTINE H1Basis_QuadNodal

  SUBROUTINE H1Basis_dQuadNodal(nvec, u, v, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    REAL(Kind=dp), PARAMETER :: c = 1D0/4D0
    INTEGER :: j
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) = -c*(1-v(j))
      grad(j,nbasis+2,1) =  c*(1-v(j))
      grad(j,nbasis+3,1) =  c*(1+v(j))
      grad(j,nbasis+4,1) = -c*(1+v(j))
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = -c*(1-u(j))
      grad(j,nbasis+2,2) = -c*(1+u(j))
      grad(j,nbasis+3,2) =  c*(1+u(j))
      grad(j,nbasis+4,2) =  c*(1-u(j))
    END DO
    
    nbasis = nbasis + 4
  END SUBROUTINE H1Basis_dQuadNodal

! --- start serendipity quad


  SUBROUTINE H1Basis_SD_QuadEdgeP(nvec, u, v, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir
    REAL(KIND=dp), PARAMETER :: c = 1/2D0

    REAL(KIND=dp) :: La, Lb
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    ! For each edge
    DO i=1,4
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb)
        DO k=1,nvec
          La = H1Basis_QuadL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_QuadL(edgedir(2,i), u(k), v(k))

          fval(k, nbasis+j) = c*(La+Lb-1)*H1Basis_Phi(j+1, Lb-La)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_SD_QuadEdgeP

  SUBROUTINE H1Basis_SD_dQuadEdgeP(nvec, u, v, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Phi, dPhi, dLa(2), dLb(2)
    REAL(KIND=dp), PARAMETER :: c = 1/2D0
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64
    
    ! For each edge
    DO i=1,4
      dLa = H1Basis_dQuadL(edgedir(1,i))
      dLb = H1Basis_dQuadL(edgedir(2,i))

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Phi, dPhi)
        DO k=1,nvec
          La = H1Basis_QuadL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_QuadL(edgedir(2,i), u(k), v(k))
          
          Phi = H1Basis_Phi(j+1, Lb-La)
          dPhi = H1Basis_dPhi(j+1,Lb-La)

          grad(k, nbasis+j, 1) = c*((dLa(1)+dLb(1))*Phi + &
                  (La+Lb-1)*dPhi*(dLb(1)-dLa(1)))
          grad(k, nbasis+j, 2) = c*((dLa(2)+dLb(2))*Phi + &
                  (La+Lb-1)*dPhi*(dLb(2)-dLa(2)))
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_SD_dQuadEdgeP

  SUBROUTINE H1Basis_SD_QuadBubbleP(nvec, u, v, pmax, nbasismax, fval, nbasis, localNumbers)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)

    INTEGER :: i,j,k
    REAL(KIND=dp) :: La, Lb, Lc
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    ! Calculate value of function without direction and return
    ! if local numbering not present
    IF (.NOT. PRESENT(localNumbers)) THEN
      DO i=2,(pmax-2)
        DO j=1,(pmax-i)-1
          !_ELMER_OMP_SIMD
          DO k=1,nvec
            fval(k,nbasis+j) = H1Basis_Phi(i,u(k))*H1Basis_Phi(j+1,v(k))
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-1,0)
      END DO
    ELSE
      DO i=2,(pmax-2)
        DO j=1,(pmax-i)-1
          !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lc)
          DO k=1,nvec
            ! Directed quad bubbles
            La = H1Basis_QuadL(localNumbers(1),u(k),v(k))
            Lb = H1Basis_QuadL(localNumbers(2),u(k),v(k))
            Lc = H1Basis_QuadL(localNumbers(4),u(k),v(k))

            ! Calculate value of function from the general form
            fval(k,nbasis+j) = H1Basis_Phi(i,Lb-La)*H1Basis_Phi(j+1,Lc-La)
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-1,0)
      END DO
    END IF
  END SUBROUTINE H1Basis_SD_QuadBubbleP

  SUBROUTINE H1Basis_SD_dQuadBubbleP(nvec, u, v, pmax, nbasismax, grad, nbasis, localNumbers)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)

    INTEGER :: i,j,k
    REAL(KIND=dp) :: La, Lb, Lc
    REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, dLbdLa, dLcdLa
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    IF (.NOT. PRESENT(localNumbers)) THEN
      DO i=2,(pmax-2)
        DO j=1,(pmax-i)-1
          !_ELMER_OMP_SIMD
          DO k=1,nvec
            ! First coordinate (Xi)
            grad(k,nbasis+j,1) = H1Basis_dPhi(i,u(k))*H1Basis_Phi(j+1,v(k))
            ! Second coordinate (Eta)
            grad(k,nbasis+j,2) = H1Basis_Phi(i,u(k))*H1Basis_dPhi(j+1,v(k))
          END DO
        END DO
        nbasis = nbasis+MAX((pmax-i)-1,0)
      END DO
    ELSE
      ! Numbering present, so use it
      dLa = H1Basis_dQuadL(localNumbers(1))
      dLb = H1Basis_dQuadL(localNumbers(2))
      dLc = H1Basis_dQuadL(localNumbers(4))

      dLBdLa = dLb-dLa
      dLcdLa = dLc-dLa

      DO i=2,(pmax-2)
        DO j=1,(pmax-i)-1
          !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lc)
          DO k=1,nvec
            La = H1Basis_QuadL(localNumbers(1),u(k),v(k))
            Lb = H1Basis_QuadL(localNumbers(2),u(k),v(k))
            Lc = H1Basis_QuadL(localNumbers(4),u(k),v(k))

            grad(k,nbasis+j,1) = H1Basis_dPhi(i,Lb-La)*(dLbdLa(1))*H1Basis_Phi(j+1,Lc-La) + &
                    H1Basis_Phi(i,Lb-La)*H1Basis_dPhi(j+1,Lc-La)*(dLcdLa(1))
            grad(k,nbasis+j,2) = H1Basis_dPhi(i,Lb-La)*(dLbdLa(2))*H1Basis_Phi(j+1,Lc-La) + &
                    H1Basis_Phi(i,Lb-La)*H1Basis_dPhi(j+1,Lc-La)*(dLcdLa(2))
          END DO
        END DO
        nbasis = nbasis+MAX((pmax-i)-1,0)
      END DO
    END IF
  END SUBROUTINE H1Basis_SD_dQuadBubbleP


! --- end serendipity quad


  SUBROUTINE H1Basis_QuadEdgeP(nvec, u, v, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir
    REAL(KIND=dp), PARAMETER :: c = 1/2D0*2

    REAL(KIND=dp) :: La, Lb, N(VECTOR_BLOCK_LENGTH,nbasismax), Na, Nb
    INTEGER :: i,j,k, nnb, node1, node2
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    ! For each edge
    DO i=1,4
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          La = H1Basis_QuadL(node1, u(k), v(k))
          Lb = H1Basis_QuadL(node2, u(k), v(k))

          Na = fval(k,node1)
          Nb = fval(k,node2)

          fval(k, nbasis+j) = c*Na*Nb*H1Basis_varPhi(j+1, Lb-La)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_QuadEdgeP

  SUBROUTINE H1Basis_dQuadEdgeP(nvec, u, v, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Phi, dPhi, dLa(2), dLb(2), N(VECTOR_BLOCK_LENGTH,nbasismax), &
           dN(VECTOR_BLOCK_LENGTH,nbasismax,3), Na, Nb, dNa(2), dNb(2)
    REAL(KIND=dp), PARAMETER :: c = 1/2D0*2
    INTEGER :: i,j,k,nnb, node1, node2
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64
    
    nnb = 0; CALL H1Basis_QuadNodal(nvec, u, v, nbasismax, n, nnb )

    ! For each edge
    DO i=1,4
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      dLa = H1Basis_dQuadL(node1)
      dLb = H1Basis_dQuadL(node2)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Phi, dPhi, Na, Nb, dNa, dNb)
        DO k=1,nvec
          La = H1Basis_QuadL(node1, u(k), v(k))
          Lb = H1Basis_QuadL(node2, u(k), v(k))
          
          Na  = N(k,node1)
          Nb  = N(k,node2)
          dNa = grad(k,node1,1:2)
          dNb = grad(k,node2,1:2)

          Phi = H1Basis_varPhi(j+1, Lb-La)
          dPhi = H1Basis_dvarPhi(j+1,Lb-La)

          grad(k, nbasis+j, 1:2) = c*(dNa*Nb*Phi + Na*dNb*Phi + &
                    Na*Nb*dPhi*(dLb-dLa))
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_dQuadEdgeP

  SUBROUTINE H1Basis_QuadBubbleP(nvec, u, v, pmax, nbasismax, fval, nbasis, localNumbers)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)

    INTEGER :: i,j,k
    REAL(KIND=dp) :: La, Lb, Lc, Pa, Pb
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64


    ! Calculate value of function without direction and return
    ! if local numbering not present
    IF (.NOT. PRESENT(localNumbers)) THEN
      DO i=2,pmax
        DO j=1,pmax-1
          !_ELMER_OMP_SIMD
          DO k=1,nvec
            fval(k,nbasis+j) = H1Basis_Phi(i,u(k))*H1Basis_Phi(j+1,v(k))
          END DO
        END DO
        nbasis = nbasis + pmax - 1
      END DO
    ELSE
      DO i=0,pmax-2
        DO j=1,pmax-1
          !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lc,Pa,Pb)
          DO k=1,nvec
            ! Directed quad bubbles
            La = H1Basis_QuadL(localNumbers(1),u(k),v(k))
            Lb = H1Basis_QuadL(localNumbers(2),u(k),v(k))
            Lc = H1Basis_QuadL(localNumbers(4),u(k),v(k))

            Pa = fval(k,localNumbers(1))
            Pb = fval(k,localNumbers(3))

            ! Calculate value of function from the general form
            fval(k,nbasis+j) = Pa*Pb*H1Basis_LegendreP(i,Lb-La)*H1Basis_LegendreP(j-1,Lc-La)
          END DO
        END DO
        nbasis = nbasis + pmax - 1
      END DO
    END IF
  END SUBROUTINE H1Basis_QuadBubbleP

  SUBROUTINE H1Basis_dQuadBubbleP(nvec, u, v, pmax, nbasismax, grad, nbasis, localNumbers)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, INTENT(IN), OPTIONAL :: localNumbers(4)

    INTEGER :: i,j,k, nnb
    REAL(KIND=dp) :: La, Lb, Lc, Legi, Legj, Pa, Pb
    REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, dLbdLa, dLcdLa, dLegi,dLegj, dPa,dPb
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: N
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    IF (.NOT. PRESENT(localNumbers)) THEN
      DO i=2,pmax
        DO j=1,pmax-1
          !_ELMER_OMP_SIMD
          DO k=1,nvec
            ! First coordinate (Xi)
            grad(k,nbasis+j,1) = H1Basis_dPhi(i,u(k))*H1Basis_Phi(j+1,v(k))
            ! Second coordinate (Eta)
            grad(k,nbasis+j,2) = H1Basis_Phi(i,u(k))*H1Basis_dPhi(j+1,v(k))
          END DO
        END DO
        nbasis = nbasis + pmax - 1
      END DO
    ELSE
      nnb = 0; CALL H1Basis_QuadNodal(nvec, u, v, nbasismax, n, nnb )

      ! Numbering present, so use it
      dLa = H1Basis_dQuadL(localNumbers(1))
      dLb = H1Basis_dQuadL(localNumbers(2))
      dLc = H1Basis_dQuadL(localNumbers(4))

      DO i=0,pmax-2
        DO j=1,pmax-1
          !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lc,Legi,Legj,dLegi,dLegj,Pa,Pb,dPa,dPb)
          DO k=1,nvec
            La = H1Basis_QuadL(localNumbers(1),u(k),v(k))
            Lb = H1Basis_QuadL(localNumbers(2),u(k),v(k))
            Lc = H1Basis_QuadL(localNumbers(4),u(k),v(k))

            Pa = N(k,localNumbers(1))
            Pb = N(k,localNumbers(3))

            dPa = grad(k,localNumbers(1),1:2)
            dPb = grad(k,localNumbers(3),1:2)

            Legi = H1Basis_LegendreP(i,Lb-La)
            Legj = H1Basis_LegendreP(j-1,Lc-La)

            dLegi = H1Basis_dLegendreP(i,Lb-La)*(dLb-dLa)
            dLegj = H1Basis_dLegendreP(j-1,Lc-La)*(dLc-dLa)

            grad(k,nbasis+j,1:2) = dPa*Pb*Legi*Legj + Pa*dPb*Legi*Legj + &
                                   Pa*Pb*dLegi*Legj + Pa*Pb*Legi*dLegj
          END DO
        END DO
        nbasis = nbasis + pmax - 1
      END DO
    END IF
  END SUBROUTINE H1Basis_dQuadBubbleP

  FUNCTION H1Basis_QuadL(node, u, v) RESULT(fval)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: node
    REAL(KIND=dp), INTENT(IN) :: u, v
    REAL(KIND=dp) :: fval
    REAL(KIND=dp), PARAMETER :: c = 1/2D0
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) &
    !_ELMER_OMP _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) NOTINBRANCH
    
    SELECT CASE (node)
    CASE (1)
      fval = c*(2-u-v)
    CASE (2)
      fval = c*(2+u-v)
    CASE (3)
      fval = c*(2+u+v)
    CASE (4)
      fval = c*(2-u+v)
    END SELECT
  END FUNCTION H1Basis_QuadL

  FUNCTION H1Basis_dQuadL(node) RESULT(grad)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: node
    ! REAL(KIND=dp), INTENT(IN) :: u, v
    REAL(KIND=dp) :: grad(2)
    REAL(KIND=dp), PARAMETER :: c = 1/2D0
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE (1)
      grad(1:2) = c*[-1,-1 ]
    CASE (2)
      grad(1:2) = c*[ 1,-1 ]
    CASE (3)
      grad(1:2) = c*[ 1, 1 ]
    CASE (4)
      grad(1:2) = c*[-1, 1 ]
    END SELECT
  END FUNCTION H1Basis_dQuadL

  SUBROUTINE H1Basis_TetraNodalP(nvec, u, v, w, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER :: j

    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0), &
            e = 1D0/SQRT(6D0), f = 1D0/SQRT(8D0), g = SQRT(3D0)*f
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,nbasis+1) = c*(1-u(j)-d*v(j)-e*w(j))
      ! Node 2
      fval(j,nbasis+2) = c*(1+u(j)-d*v(j)-e*w(j))
      ! Node 3
      fval(j,nbasis+3) = d*(v(j)-f*w(j))
      ! Node 4
      fval(j,nbasis+4) = g*w(j)
    END DO
    nbasis = nbasis + 4
  END SUBROUTINE H1Basis_TetraNodalP

  FUNCTION H1Basis_TetraL(node, u, v, w) RESULT(fval)
    IMPLICIT NONE

    ! Parameters 
    INTEGER, INTENT(IN) :: node
    REAL(KIND=dp), INTENT(IN) :: u,v,w
    ! Result
    REAL(KIND=dp) :: fval
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0), &
            e = 1D0/SQRT(6D0), f = 1D0/SQRT(8D0), g = SQRT(3D0)*f
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) &
    !_ELMER_OMP _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) _ELMER_LINEAR_REF(w) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1)
      fval = c*(1-u-d*v-e*w)
    CASE(2)
      fval = c*(1+u-d*v-e*w)
    CASE(3)
      fval = d*(v-f*w)
    CASE(4)
      fval = g*w      
    END SELECT
  END FUNCTION H1Basis_TetraL
  
  SUBROUTINE H1Basis_dTetraNodalP(nvec, u, v, w, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    ! Variables
    INTEGER, INTENT(IN) :: nbasismax
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER :: j
    REAL(KIND=dp), PARAMETER :: half=1d0/2, a = 1/SQRT(3._dp), b=-a/2, c = -SQRT(6._dp)/12, &
                  d = -3*c
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) =-half
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,1) = half
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,1) = 0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+4,1) = 0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = b
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,2) = b
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,2) = a
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+4,2) = 0
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,3) = c
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,3) = c
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,3) = c
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+4,3) = d
    END DO
    nbasis = nbasis + 4
  END SUBROUTINE H1Basis_dTetraNodalP
  
  FUNCTION H1Basis_dTetraL(node) RESULT(grad)
    IMPLICIT NONE

    ! Parameters 
    INTEGER, INTENT(IN) :: node
    ! Result
    REAL(KIND=dp) :: grad(3)
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0), &
            e = 1D0/SQRT(6D0), f = 1D0/SQRT(8D0), g = SQRT(3._dp)*f
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1)
      grad = c*[-1D0, -d, -e]
    CASE(2)
      grad = c*[1D0, -d, -e]
    CASE(3)
      grad = d*[0D0, 1D0, -f]
    CASE(4)
      grad = [0D0, 0D0, g]
    END SELECT
  END FUNCTION H1Basis_dTetraL

  SUBROUTINE H1Basis_TetraEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each edge
    DO i=1,6
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb)
        DO k=1,nvec
          La = H1Basis_TetraL(edgedir(1,i), u(k), v(k), w(k))
          Lb = H1Basis_TetraL(edgedir(2,i), u(k), v(k), w(k))
          fval(k, nbasis+j) = La*Lb*H1Basis_varPhi(j+1,Lb-La)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_TetraEdgeP
  
  SUBROUTINE H1Basis_dTetraEdgeP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, vPhi, dVPhi, dLa(3), dLb(3)
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each edge
    DO i=1,6
      dLa = H1Basis_dTetraL(edgedir(1,i))
      dLb = H1Basis_dTetraL(edgedir(2,i))

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, vPhi, dVPhi)
        DO k=1,nvec
          La = H1Basis_TetraL(edgedir(1,i),u(k),v(k),w(k))
          Lb = H1Basis_TetraL(edgedir(2,i),u(k),v(k),w(k))
          vPhi = H1Basis_varPhi(j+1,Lb-La)
          dVPhi = H1Basis_dVarPhi(j+1,Lb-La)

          grad(k, nbasis+j, 1) = dLa(1)*Lb*vPhi + &
                  La*dLb(1)*vPhi + La*Lb*dVPhi*(dLb(1)-dLa(1))
          grad(k, nbasis+j, 2) = dLa(2)*Lb*vPhi+ &
                  La*dLb(2)*vPhi + La*Lb*dVPhi*(dLb(2)-dLa(2))
          grad(k, nbasis+j, 3) = dLa(3)*Lb*vPhi+ &
                  La*dLb(3)*vPhi + La*Lb*dVPhi*(dLb(3)-dLa(3))
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_dTetraEdgeP

  SUBROUTINE H1Basis_TetraFaceP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each face
    DO i=1,4
      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-j-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc)
          DO l=1,nvec
            La=H1Basis_TetraL(facedir(1,i),u(l),v(l),w(l))
            Lb=H1Basis_TetraL(facedir(2,i),u(l),v(l),w(l))
            Lc=H1Basis_TetraL(facedir(3,i),u(l),v(l),w(l))
    
            fval(l,nbasis+k) = La*Lb*Lc*H1Basis_LegendreP(j,Lb-La)* &
                                        H1Basis_LegendreP(k-1,2*Lc-1)
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-3 + 1)
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_TetraFaceP
  
  SUBROUTINE H1Basis_dTetraFaceP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, dLa(3), dLb(3), dLc(3), Lb_La, Lc_1, a, b
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each face
    DO i=1,4
      dLa = H1Basis_dTetraL(facedir(1,i))
      dLb = H1Basis_dTetraL(facedir(2,i))
      dLc = H1Basis_dTetraL(facedir(3,i))

      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-j-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Lb_La, Lc_1, a, b)
          DO l=1,nvec
            La=H1Basis_TetraL(facedir(1,i),u(l),v(l),w(l))
            Lb=H1Basis_TetraL(facedir(2,i),u(l),v(l),w(l))
            Lc=H1Basis_TetraL(facedir(3,i),u(l),v(l),w(l))
            Lb_La=Lb-La
            Lc_1=2*Lc-1
            a=H1Basis_LegendreP(j,Lb_La)
            b=H1Basis_LegendreP(k-1,Lc_1)

            grad(l,nbasis+k,1) = dLa(1)*Lb*Lc*a*b + La*dLb(1)*Lc*a*b + &
                    La*Lb*dLc(1)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(j,Lb_La)*(dLb(1)-dLa(1))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(k-1,Lc_1)*2*dLc(1)
            grad(l,nbasis+k,2) = dLa(2)*Lb*Lc*a*b + La*dLb(2)*Lc*a*b + &
                    La*Lb*dLc(2)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(j,Lb_La)*(dLb(2)-dLa(2))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(k-1,Lc_1)*2*dLc(2)
            grad(l,nbasis+k,3) = dLa(3)*Lb*Lc*a*b + La*dLb(3)*Lc*a*b + &
                    La*Lb*dLc(3)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(j,Lb_La)*(dLb(3)-dLa(3))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(k-1,Lc_1)*2*dLc(3)
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-3 + 1)
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dTetraFaceP
  
  SUBROUTINE H1Basis_TetraBubbleP(nvec, u, v, w, pmax, nbasismax, fval, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    
    ! Variables
    INTEGER :: i, j, k, l
    ! Variables
    REAL(KIND=dp) :: L1, L2, L3, L4, L2_L1, L3_1, L4_1
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=0,pmax-4
      DO j=0,pmax-i-4
        DO k=1,pmax-i-j-3
          !_ELMER_OMP_SIMD PRIVATE(L1,L2,L3,L4,L2_L1,L3_1,L4_1)
          DO l=1,nvec
            L1 = H1Basis_TetraL(1,u(l),v(l),w(l))
            L2 = H1Basis_TetraL(2,u(l),v(l),w(l))
            L3 = H1Basis_TetraL(3,u(l),v(l),w(l))
            L4 = H1Basis_TetraL(4,u(l),v(l),w(l))
            L2_L1 = L2-L1
            L3_1 = 2*L3-1
            L4_1 = 2*L4-1

            fval(l,nbasis+k) = L1*L2*L3*L4*&
                    H1Basis_LegendreP(i,L2_L1)*&
                    H1Basis_LegendreP(j,L3_1)*&
                    H1Basis_LegendreP(k-1,L4_1)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-j-4) + 1
        nbasis = nbasis + MAX(pmax-i-j-3,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_TetraBubbleP

  SUBROUTINE H1Basis_dTetraBubbleP(nvec, u, v, w, pmax, nbasismax, grad, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis

    ! Variables
    INTEGER :: i, j, k, l
    REAL(KIND=dp) :: L1, L2, L3, L4, L2_L1, L3_1, L4_1, a, b, c 
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    DO i=0,pmax-4
      DO j=0,pmax-i-4
        DO k=1,pmax-i-j-3
          !_ELMER_OMP_SIMD PRIVATE(L1,L2,L3,L4,L2_L1,L3_1,L4_1,a,b,c)
          DO l=1,nvec
            L1 = H1Basis_TetraL(1,u(l),v(l),w(l))
            L2 = H1Basis_TetraL(2,u(l),v(l),w(l))
            L3 = H1Basis_TetraL(3,u(l),v(l),w(l))
            L4 = H1Basis_TetraL(4,u(l),v(l),w(l))
            L2_L1 = L2 - L1
            L3_1 = 2*L3 - 1
            L4_1 = 2*L4 - 1
            a = H1Basis_LegendreP(i,L2_L1)
            b = H1Basis_LegendreP(j,L3_1)
            c = H1Basis_LegendreP(k-1,L4_1)

            ! Gradients of tetrahedral bubble basis functions 
            grad(l,nbasis+k,1) = -1d0/2*L2*L3*L4*a*b*c + &
                    1d0/2*L1*L3*L4*a*b*c + &
                    L1*L2*L3*L4*H1Basis_dLegendreP(i,L2_L1)*b*c
            grad(l,nbasis+k,2) = -SQRT(3d0)/6*L2*L3*L4*a*b*c - &
                    SQRT(3d0)/6*L1*L3*L4*a*b*c + &
                    SQRT(3d0)/3*L1*L2*L4*a*b*c &
                    + 2*SQRT(3d0)/3*L1*L2*L3*L4*a*H1Basis_dLegendreP(j,L3_1)*c
            grad(l,nbasis+k,3) = -SQRT(6d0)/12*L2*L3*L4*a*b*c - &
                    SQRT(6d0)/12*L1*L3*L4*a*b*c - &
                    SQRT(6d0)/12*L1*L2*L4*a*b*c &
                    + SQRT(6d0)/4*L1*L2*L3*a*b*c - &
                    SQRT(6d0)/6*L1*L2*L3*L4*a*H1Basis_dLegendreP(j,L3_1)*c &
                    + SQRT(6d0)/2*L1*L2*L3*L4*a*b*H1Basis_dLegendreP(k-1,L4_1)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-j-4) + 1
        nbasis = nbasis + MAX(pmax-i-j-3,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dTetraBubbleP


! --- start serendipity wedge

  
  SUBROUTINE H1Basis_SD_WedgeEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Na, Nb
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each triangle face edge
    DO i=1,6
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_WedgeL(edgedir(2,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          fval(k, nbasis+j) = c*La*Lb*H1Basis_varPhi(j+1,Lb-La)*(1+Na+Nb)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
      
    ! For each square face edge
    DO i=7,9
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Na, Nb)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          fval(k, nbasis+j) = La*H1Basis_Phi(j+1, Nb-Na)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_SD_WedgeEdgeP
  
  SUBROUTINE H1Basis_SD_dWedgeEdgeP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Na, Nb, Phi, dPhi, vPhi, dVPhi, &
          dLa(3), dLb(3), dNa(3), dNb(3), NaNb
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each triangle face edge
    DO i=1,6
      dLa = H1Basis_dWedgeL(edgedir(1,i))
      dLb = H1Basis_dWedgeL(edgedir(2,i))
      dNa = H1Basis_dWedgeH(edgedir(1,i))
      dNb = H1Basis_dWedgeH(edgedir(2,i))
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, vPhi, dVPhi, NaNb)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_WedgeL(edgedir(2,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          vPhi=H1Basis_varPhi(j+1,Lb-La)
          dVPhi=H1Basis_dVarPhi(j+1,Lb-La)
          NaNb=1+Na+Nb
          grad(k, nbasis+j,1) = c*dLa(1)*Lb*vPhi*NaNb + &
                c*La*dLb(1)*vPhi*NaNb + c*La*Lb*dVPhi*(dLb(1)-dLa(1))*NaNb
          grad(k, nbasis+j,2) = c*dLa(2)*Lb*vPhi*NaNb + &
                c*La*dLb(2)*vPhi*NaNb + c*La*Lb*dVPhi*(dLb(2)-dLa(2))*NaNb
          grad(k, nbasis+j,3) = c*La*Lb*vPhi*(dNb(3)+dNa(3))
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
      
    ! For each square face edge
    DO i=7,9
      dLa = H1Basis_dWedgeL(edgedir(1,i))
      dNa = H1Basis_dWedgeH(edgedir(1,i))
      dNb = H1Basis_dWedgeH(edgedir(2,i))
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Na, Nb, Phi, dPhi)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          Phi = H1Basis_Phi(j+1, Nb-Na)
          dPhi = H1Basis_dPhi(j+1,Nb-Na)
          
          grad(k, nbasis+j,1) = dLa(1)*Phi+La*dPhi*(dNb(1)-dNa(1))
          grad(k, nbasis+j,2) = dLa(2)*Phi+La*dPhi*(dNb(2)-dNa(2))
          grad(k, nbasis+j,3) = dLa(3)*Phi+La*dPhi*(dNb(3)-dNa(3)) 
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_SD_dWedgeEdgeP
  
  SUBROUTINE H1Basis_SD_WedgeFaceP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Na, Nc
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l
    LOGICAL :: nonpermuted
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! Triangle faces
    DO i=1,2      
      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-j-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Na)
          DO l=1,nvec
            La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
            Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
            Lc = H1Basis_WedgeL(facedir(3,i), u(l), v(l))
            Na = H1Basis_WedgeH(facedir(1,i), w(l))

            fval(l,nbasis+k) = c*(1+2*Na)*La*Lb*Lc* &
                                            H1Basis_LegendreP(j, Lb-La)* &
                                            H1Basis_LegendreP(k-1, 2*Lc-1)
          END DO
        END DO
        ! nbasis = nbasis + (p-j-3) + 1
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
    ! Square faces
    DO i=3,5
      ! First and second node must form an edge in upper or lower triangle
      nonpermuted = (facedir(1,i) >= 1 .AND. facedir(1,i) <= 3 .AND.&
                     facedir(2,i) >= 1 .AND. facedir(2,i) <= 3) .OR. &
                    (facedir(1,i) >= 4 .AND. facedir(1,i) <= 6 .AND.&
                     facedir(2,i) >= 4 .AND. facedir(2,i) <= 6)
      
      DO j=2,pmax(i)-2
        DO k=1,pmax(i)-j-1
          IF (nonpermuted) THEN
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc)
            DO l=1,nvec
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(4,i), w(l))
              fval(l,nbasis+k) = La*Lb*H1Basis_varPhi(j, Lb-La)* &
                                         H1Basis_Phi(k+1, Nc-Na)
            END DO
          ELSE
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc)
            DO l=1,nvec
              ! Numbering permuted, switch indices j,k and nodes 2 and 4
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(4,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(2,i), w(l))
              fval(l,nbasis+k) = La*Lb*H1Basis_varPhi(k+1, Lb-La)* &
                                         H1Basis_Phi(j, Nc-Na)
            END DO
          END IF
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_WedgeFaceP
  
  SUBROUTINE H1Basis_SD_dWedgeFaceP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Na, Nc, LegPLbLa, LegP2Lc1, &
            cNa, vPhi, Phi, dLa(3), dLb(3), dLc(3), dNa(3), dNc(3)
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l
    LOGICAL :: nonpermuted
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! Triangle faces
    DO i=1,2
      dLa = H1Basis_dWedgeL(facedir(1,i))
      dLb = H1Basis_dWedgeL(facedir(2,i))
      dLc = H1Basis_dWedgeL(facedir(3,i))
      dNa = H1Basis_dWedgeH(facedir(1,i))
      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-j-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Na, LegPLbLa, LegP2Lc1, cNa)
          DO l=1,nvec
            La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
            Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
            Lc = H1Basis_WedgeL(facedir(3,i), u(l), v(l))
            Na = H1Basis_WedgeH(facedir(1,i), w(l))

            LegPLbLa = H1Basis_LegendreP(j, Lb-La)
            LegP2Lc1 = H1Basis_LegendreP(k-1, 2*Lc-1)
            cNa = c*(1+2*Na)

            grad(l,nbasis+k,1) = cNa*dLa(1)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*dLb(1)*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*dLc(1)*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(1)-dLa(1))*LegP2Lc1 + &
                    cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k-1, 2*Lc-1)*(2*dLc(1))
            grad(l,nbasis+k,2) = cNa*dLa(2)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*dLb(2)*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*dLc(2)*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(2)-dLa(2))*LegP2Lc1 + &
                    cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k-1, 2*Lc-1)*(2*dLc(2))
            grad(l,nbasis+k,3) = 2*c*dNa(3)*La*Lb*Lc*LegPLbLa*LegP2Lc1
          END DO
        END DO
        ! nbasis = nbasis + (p-j-3) + 1
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
    ! Square faces
    DO i=3,5
      ! First and second node must form an edge in upper or lower triangle
      nonpermuted = (facedir(1,i) >= 1 .AND. facedir(1,i) <= 3 .AND.&
                     facedir(2,i) >= 1 .AND. facedir(2,i) <= 3) .OR. &
                    (facedir(1,i) >= 4 .AND. facedir(1,i) <= 6 .AND.&
                     facedir(2,i) >= 4 .AND. facedir(2,i) <= 6)
      
      IF (nonpermuted) THEN
        dLa = H1Basis_dWedgeL(facedir(1,i))
        dLb = H1Basis_dWedgeL(facedir(2,i))
        dNa = H1Basis_dWedgeH(facedir(1,i))
        dNc = H1Basis_dWedgeH(facedir(4,i))
      ELSE
        dLa = H1Basis_dWedgeL(facedir(1,i))
        dLb = H1Basis_dWedgeL(facedir(4,i))
        dNa = H1Basis_dWedgeH(facedir(1,i))
        dNc = H1Basis_dWedgeH(facedir(2,i))
      END IF

      DO j=2,pmax(i)-2
        DO k=1,pmax(i)-j-1
          IF (nonpermuted) THEN
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc,vPhi,Phi)
            DO l=1,nvec
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(4,i), w(l))
              
              vPhi = H1Basis_varPhi(j, Lb-La)
              Phi = H1Basis_Phi(k+1, Nc-Na)

              grad(l,nbasis+k,1) = dLa(1)*Lb*vPhi*Phi + La*dLb(1)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(j, Lb-La)*(dLb(1)-dLa(1))*Phi
              grad(l,nbasis+k,2) = dLa(2)*Lb*vPhi*Phi + La*dLb(2)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(j, Lb-La)*(dLb(2)-dLa(2))*Phi
              grad(l,nbasis+k,3) = La*Lb*vPhi*H1Basis_dPhi(k+1, Nc-Na)*(dNc(3)-dNa(3))
            END DO
          ELSE
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc,vPhi,Phi)
            DO l=1,nvec
              ! Numbering permuted, switch indices j,k and nodes 2 and 4
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(4,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(2,i), w(l))
              
              vPhi = H1Basis_varPhi(k+1, Lb-La)
              Phi = H1Basis_Phi(j, Nc-Na)

              grad(l,nbasis+k,1) = dLa(1)*Lb*vPhi*Phi + La*dLb(1)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(k+1, Lb-La)*(dLb(1)-dLa(1))*Phi
              grad(l,nbasis+k,2) = dLa(2)*Lb*vPhi*Phi + La*dLb(2)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(k+1, Lb-La)*(dLb(2)-dLa(2))*Phi
              grad(l,nbasis+k,3) = La*Lb*vPhi*H1Basis_dPhi(j, Nc-Na)*(dNc(3)-dNa(3))
            END DO
          END IF
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_dWedgeFaceP
  
  SUBROUTINE H1Basis_SD_WedgeBubbleP(nvec, u, v, w, pmax, nbasismax, fval, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    
    INTEGER :: i, j, k, l
    REAL(KIND=dp) :: L1, L2, L3, L2_L1, L3_1
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=0,pmax-5
      DO j=0,pmax-5-i
        DO k=1,pmax-4-i-j
          !_ELMER_OMP_SIMD PRIVATE(L1, L2, L3, L2_L1, L3_1)
          DO l=1,nvec
            L1 = H1Basis_WedgeL(1,u(l),v(l))
            L2 = H1Basis_WedgeL(2,u(l),v(l))
            L3 = H1Basis_WedgeL(3,u(l),v(l))
            L2_L1 = L2-L1
            L3_1 = 2*L3-1

            ! Get value of bubble function
            fval(l,nbasis+k) = L1*L2*L3*H1Basis_LegendreP(i,L2_L1)*&
                    H1Basis_LegendreP(j,L3_1)*H1Basis_Phi(k+1,w(l))
          END DO
        END DO
        ! nbasis = nbasis + (pmax-3-i-j)-2+1
        nbasis = nbasis + MAX(pmax-4-i-j,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_WedgeBubbleP

  SUBROUTINE H1Basis_SD_dWedgeBubbleP(nvec, u, v, w, pmax, nbasismax, grad, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    
    ! Parameters
    INTEGER :: i,j,k,l
    ! Variables
    REAL(KIND=dp) :: L1,L2,L3,Legi,Legj,phiW,L2_L1,L3_1
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64      

    DO i=0,pmax-5
      DO j=0,pmax-5-i
        DO k=1,pmax-4-i-j
          !_ELMER_OMP_SIMD PRIVATE(L1, L2, L3, Legi, Legj, phiW, L2_L1, L3_1)
          DO l=1,nvec

            ! Values of function L
            L1 = H1Basis_WedgeL(1,u(l),v(l))
            L2 = H1Basis_WedgeL(2,u(l),v(l))
            L3 = H1Basis_WedgeL(3,u(l),v(l))
            
            L2_L1 = L2-L1
            L3_1 = 2d0*L3-1

            Legi = H1Basis_LegendreP(i,L2_L1)
            Legj = H1Basis_LegendreP(j,L3_1)
            phiW = H1Basis_Phi(k+1,w(l))
            
            grad(l,nbasis+k,1) = ((-1d0/2)*L2*L3*Legi*Legj + &
                    L1*(1d0/2)*L3*Legi*Legj +&
                    L1*L2*L3*H1Basis_dLegendreP(i,L2_L1)*Legj)*phiW
            grad(l,nbasis+k,2) = ((-SQRT(3d0)/6)*L2*L3*Legi*Legj + &
                    L1*(-SQRT(3d0)/6)*L3*Legi*Legj +&
                    L1*L2*(SQRT(3D0)/3)*Legi*Legj + &
                    L1*L2*L3*Legi*H1Basis_dLegendreP(j,L3_1)*(2d0*SQRT(3D0)/3))*phiW 
            grad(l,nbasis+k,3) = L1*L2*L3*Legi*Legj*H1Basis_dPhi(k+1,w(l))
          END DO
        END DO
        ! nbasis = nbasis + (pmax-3-i-j)-2+1
        nbasis = nbasis + MAX(pmax-4-i-j,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_dWedgeBubbleP


! --- end serendipity wedge


  SUBROUTINE H1Basis_WedgeNodalP(nvec, u, v, w, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    ! Result
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER :: j

    REAL(KIND=dp), PARAMETER :: c = 1D0/4D0, d = 1D0/SQRT(3D0), &
            e = SQRT(3d0)/6D0
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,nbasis+1) = c*(1-u(j)-d*v(j))*(1-w(j))
      ! Node 2
      fval(j,nbasis+2) = c*(1+u(j)-d*v(j))*(1-w(j))
      ! Node 3
      fval(j,nbasis+3) = e*v(j)*(1-w(j))
      ! Node 4
      fval(j,nbasis+4) = c*(1-u(j)-d*v(j))*(1+w(j))
      ! Node 5
      fval(j,nbasis+5) = c*(1+u(j)-d*v(j))*(1+w(j))
      ! Node 6
      fval(j,nbasis+6) = e*v(j)*(1+w(j))
    END DO
    nbasis = nbasis + 6
  END SUBROUTINE H1Basis_WedgeNodalP

  SUBROUTINE H1Basis_dWedgeNodalP(nvec, u, v, w, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    
    ! Variables
    INTEGER :: j
    REAL(KIND=dp), PARAMETER :: c = 1D0/4D0, d = 1D0/SQRT(3D0), &
            e = SQRT(3d0)/6D0, f = SQRT(3d0)/12D0
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) = -c*(1-w(j))
      grad(j,nbasis+2,1) =  c*(1-w(j))  
      grad(j,nbasis+3,1) =  0
      grad(j,nbasis+4,1) = -c*(1+w(j))
      grad(j,nbasis+5,1) =  c*(1+w(j))
      grad(j,nbasis+6,1) =  0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2)=  -f*(1-w(j))
      grad(j,nbasis+2,2) = -f*(1-w(j))
      grad(j,nbasis+3,2) =  e*(1-w(j))
      grad(j,nbasis+4,2) = -f*(1+w(j))
      grad(j,nbasis+5,2) = -f*(1+w(j))
      grad(j,nbasis+6,2) =  e*(1+w(j))
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,3) = -c*(1-u(j)-d*v(j))
      grad(j,nbasis+2,3) = -c*(1+u(j)-d*v(j))
      grad(j,nbasis+3,3) = -e*v(j)
      grad(j,nbasis+4,3) =  c*(1-u(j)-d*v(j))
      grad(j,nbasis+5,3) =  c*(1+u(j)-d*v(j))
      grad(j,nbasis+6,3) =  e*v(j)
    END DO
        nbasis = nbasis + 6
  END SUBROUTINE H1Basis_dWedgeNodalP

  FUNCTION H1Basis_WedgeL(node, u, v) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: node
    REAL(KIND=dp), INTENT(IN) :: u,v
    ! Result
    REAL(KIND=dp) :: fval
    REAL(KIND=dp), PARAMETER :: a=1/SQRT(3.0_dp)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) &
    !_ELMER_OMP _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) NOTINBRANCH
    
    SELECT CASE(node)
    CASE (1,4)
      fval = (1-u-a*v)/2
    CASE (2,5)
      fval = (1+u-a*v)/2
    CASE (3,6)
      fval = a*v
    END SELECT
  END FUNCTION H1Basis_WedgeL

  FUNCTION H1Basis_dWedgeL(node) RESULT(grad)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: node
    ! Result
    REAL(KIND=dp) :: grad(3)
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0, d = 1D0/SQRT(3D0)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1,4)
      grad = c*[REAL(-1,dp), REAL(-d,dp), 0D0]
    CASE(2,5)
      grad = c*[REAL(1,dp), REAL(-d,dp), 0D0]
    CASE(3,6)
      grad =   [REAL(0,dp), REAL(d,dp), 0D0]
    END SELECT    
  END FUNCTION H1Basis_dWedgeL

  FUNCTION H1Basis_WedgeH(node, w) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: node
    REAL(KIND=dp), INTENT(IN) :: w
    ! Result
    REAL(KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) _ELMER_LINEAR_REF(w) NOTINBRANCH
    
    SELECT CASE(node)
    CASE (1,2,3)
      fval = -w/2
    CASE (4,5,6)
      fval =  w/2
    END SELECT
  END FUNCTION H1Basis_WedgeH

  FUNCTION H1Basis_dWedgeH(node) RESULT(grad)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: node
    ! Result
    REAL(KIND=dp) :: grad(3)
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE (1,2,3)
      grad = [0D0, 0D0, -c]
    CASE (4,5,6)
      grad = [0D0, 0D0,  c]
    END SELECT
  END FUNCTION H1Basis_dWedgeH
  
  SUBROUTINE H1Basis_WedgeEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Na, Nb
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l, node1, node2
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each triangle face edge
    DO i=1,6
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          La = H1Basis_WedgeL(node1, u(k), v(k))
          Lb = H1Basis_WedgeL(node2, u(k), v(k))
          Na = fval(k,node1)
          Nb = fval(k,node2)
          fval(k, nbasis+j) = Na*Nb*H1Basis_varPhi(j+1,Lb-La)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
      
    ! For each square face edge
    DO i=7,9
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          La = H1Basis_WedgeH(node1, w(k))
          Lb = H1Basis_WedgeH(node2, w(k))
          Na = fval(k,node1)
          Nb = fval(k,node2)
          fval(k, nbasis+j) = Na*Nb*H1Basis_varPhi(j+1, Lb-La)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_WedgeEdgeP
  
  SUBROUTINE H1Basis_dWedgeEdgeP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Na, Nb, Phi, dPhi, vPhi, dVPhi, &
          dLa(3), dLb(3), dNa(3), dNb(3), N(VECTOR_BLOCK_LENGTH,nbasismax)
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k, node1, node2, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    nnb = 0; CALL H1Basis_WedgeNodalP(nvec, u, v, w, nbasismax, n, nnb )

    ! For each triangle face edge
    DO i=1,6
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      dLa = H1Basis_dWedgeL(node1)
      dLb = H1Basis_dWedgeL(node2)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, dNa, dNb, vPhi, dVPhi)
        DO k=1,nvec
          La = H1Basis_WedgeL(node1, u(k), v(k))
          Lb = H1Basis_WedgeL(node2, u(k), v(k))

          Na = N(k,node1)
          Nb = N(k,node2)
          dNa = grad(k,node1,:)
          dNb = grad(k,node2,:)

          vPhi = H1Basis_varPhi(j+1,Lb-La)
          dVPhi = H1Basis_dVarPhi(j+1,Lb-La)
          ! fval(k, nbasis+j-1) = c*La*Lb*H1Basis_varPhi(j,Lb-La)*(1+Na+Nb)

          grad(k, nbasis+j,:) = dNa*Nb*vPhi + Na*dNb*vPhi + Na*Nb*dvPhi*(dLb-dLa)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
      
    ! For each square face edge
    DO i=7,9
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      dLa = H1Basis_dWedgeH(node1)
      dLb = H1Basis_dWedgeH(node2)
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, dNa, dNb, vPhi, dvPhi)
        DO k=1,nvec
          La = H1Basis_WedgeH(node1, w(k))
          Lb = H1Basis_WedgeH(node2, w(k))

          vPhi = H1Basis_varPhi(j+1, Lb-La)
          dvPhi = H1Basis_dvarPhi(j+1, Lb-La)

          Na = N(k,node1)
          Nb = N(k,node2)
          dNa = grad(k,node1,:)
          dNb = grad(k,node2,:)

          grad(k, nbasis+j,:) = dNa*Nb*vPhi + Na*dNb*vPhi + Na*Nb*dvPhi*(dLb-dLa)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_dWedgeEdgeP
  
  SUBROUTINE H1Basis_WedgeFaceP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Na, Nb, Lha, Lhc
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l, node1, node2, node3, node4
    LOGICAL :: nonpermuted
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! Triangle faces
    DO i=1,2      
      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-j-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Na)
          DO l=1,nvec
            La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
            Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
            Lc = H1Basis_WedgeL(facedir(3,i), u(l), v(l))
            Na = H1Basis_WedgeH(facedir(1,i), w(l))

            fval(l,nbasis+k) = c*(1+2*Na)*La*Lb*Lc* &
                                            H1Basis_LegendreP(j, Lb-La)* &
                                            H1Basis_LegendreP(k-1, 2*Lc-1)
          END DO
        END DO
        ! nbasis = nbasis + (p-j-3) + 1
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
    ! Square faces
    DO i=3,5
      ! First and second node must form an edge in upper or lower triangle
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      node4 = facedir(4,i)

      nonpermuted = (node1 >= 1 .AND. node1 <= 3 .AND.&
                     node2 >= 1 .AND. node2 <= 3) .OR. &
                    (node1 >= 4 .AND. node1 <= 6 .AND.&
                     node2 >= 4 .AND. node2 <= 6)

      
      DO j=0,pmax(i)-2
        DO k=1,pmax(i)-1
          IF (nonpermuted) THEN
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nb,Lha,Lhc)
            DO l=1,nvec
              La = H1Basis_WedgeL(node1, u(l), v(l))
              Lb = H1Basis_WedgeL(node2, u(l), v(l))

              Lha = H1Basis_WedgeH(node1, w(l))
              Lhc = H1Basis_WedgeH(node4, w(l))

              Na = fval(l,node1)
              Nb = fval(l,node3)

              fval(l,nbasis+k) = Na*Nb*H1Basis_LegendreP(j, Lb-La)*H1Basis_LegendreP(k-1,Lhc-Lha)
            END DO
          ELSE
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nb,Lha,Lhc)
            DO l=1,nvec
              ! Numbering permuted, switch indices j,k and nodes 2 and 4
              La = H1Basis_WedgeL(node1, u(l), v(l))
              Lb = H1Basis_WedgeL(node4, u(l), v(l))

              Lha = H1Basis_WedgeH(node1, w(l))
              Lhc = H1Basis_WedgeH(node2, w(l))

              Na = fval(l,node1)
              Nb = fval(l,node3)

              fval(l,nbasis+k) = Na*Nb*H1Basis_LegendreP(k-1, Lb-La)*H1Basis_LegendreP(j,Lhc-Lha)
            END DO
          END IF
        END DO
        nbasis = nbasis + MAX(pmax(i)-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_WedgeFaceP
  
  SUBROUTINE H1Basis_dWedgeFaceP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Lha, Lhc, Na, Nb, Nc, LegPLbLa, LegP2Lc1, &
            cNa, vPhi, Phi, dLa(3), dLb(3), dLc(3), dLhc(3), dLha(3), dNc(3), &
                 Legj, Legk, dLegj(3), dLegk(3), dNa(3), dNb(3), N(VECTOR_BLOCK_LENGTH, nbasismax) 
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l, node1, node2, node3, node4, nnb
    LOGICAL :: nonpermuted
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    nnb = 0; CALL H1Basis_WedgeNodalP(nvec, u, v, w, nbasismax, n, nnb)

    ! Triangle faces
    DO i=1,2
      dLa = H1Basis_dWedgeL(facedir(1,i))
      dLb = H1Basis_dWedgeL(facedir(2,i))
      dLc = H1Basis_dWedgeL(facedir(3,i))
      dNa = H1Basis_dWedgeH(facedir(1,i))
      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-j-2
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Na, LegPLbLa, LegP2Lc1, cNa)
          DO l=1,nvec
            La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
            Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
            Lc = H1Basis_WedgeL(facedir(3,i), u(l), v(l))
            Na = H1Basis_WedgeH(facedir(1,i), w(l))

            LegPLbLa = H1Basis_LegendreP(j, Lb-La)
            LegP2Lc1 = H1Basis_LegendreP(k-1, 2*Lc-1)
            cNa = c*(1+2*Na)

            grad(l,nbasis+k,1) = cNa*dLa(1)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*dLb(1)*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*dLc(1)*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(1)-dLa(1))*LegP2Lc1 + &
                    cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k-1, 2*Lc-1)*(2*dLc(1))
            ! + 2*c*dNa(1)*La*Lb*Lc*LegPLbLa*LegP2Lc1 ! Term always zero
            grad(l,nbasis+k,2) = cNa*dLa(2)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*dLb(2)*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*dLc(2)*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(2)-dLa(2))*LegP2Lc1 + &
                    cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k-1, 2*Lc-1)*(2*dLc(2))
            ! + 2*c*dNa(2)*La*Lb*Lc*LegPLbLa*LegP2Lc1 ! Term always zero
            grad(l,nbasis+k,3) = 2*c*dNa(3)*La*Lb*Lc*LegPLbLa*LegP2Lc1
            ! Rest of the terms always zero
            ! cNa*dLa(3)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
            ! cNa*La*dLb(3)*Lc*LegPLbLa*LegP2Lc1 + &
            ! cNa*La*Lb*dLc(3)*LegPLbLa*LegP2Lc1 + &
            ! cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(3)-dLa(3))*LegP2Lc1 + &
            ! cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k, 2*Lc-1)*(2*dLc(3))
          END DO
        END DO
        ! nbasis = nbasis + (p-j-3) + 1
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO

    ! Square faces
    DO i=3,5

      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      node4 = facedir(4,i)

      ! First and second node must form an edge in upper or lower triangle
      nonpermuted = (node1 >= 1 .AND. node1 <= 3 .AND.&
                     node2 >= 1 .AND. node2 <= 3) .OR. &
                    (node1 >= 4 .AND. node1 <= 6 .AND.&
                     node2 >= 4 .AND. node2 <= 6)
      
      IF (nonpermuted) THEN
        dLa  = H1Basis_dWedgeL(node1)
        dLb  = H1Basis_dWedgeL(node2)
        dLha = H1Basis_dWedgeH(node1)
        dLhc = H1Basis_dWedgeH(node4)
      ELSE
        dLa  = H1Basis_dWedgeL(node1)
        dLb  = H1Basis_dWedgeL(node4)
        dLha = H1Basis_dWedgeH(node1)
        dLhc = H1Basis_dWedgeH(node2)
      END IF

      DO j=0,pmax(i)-2
        DO k=1,pmax(i)-1
          IF (nonpermuted) THEN
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lha,Lhc,Na,Nb,dNa,dNb,Legj,Legk,dLegj,dLegk)
            DO l=1,nvec
              La = H1Basis_WedgeL(node1, u(l), v(l))
              Lb = H1Basis_WedgeL(node2, u(l), v(l))

              Lha = H1Basis_WedgeH(node1, w(l))
              Lhc = H1Basis_WedgeH(node4, w(l))

              Na = N(l,node1)
              Nb = N(l,node3)

              dNa = grad(l,node1,:)
              dNb = grad(l,node3,:)

              Legj = H1Basis_LegendreP(j,Lb-La)
              Legk = H1Basis_LegendreP(k-1,Lhc-Lha)

              dLegj = H1Basis_dLegendreP(j,Lb-La)*(dLb-dLa)
              dLegk = H1Basis_dLegendreP(k-1,Lhc-Lha)*(dLhc-dLha)

              grad(l,nbasis+k,:) = dNa*Nb*Legj*Legk + Na*dNb*Legj*Legk + &
                     Na*Nb*dLegj*Legk + Na*Nb*Legj*dLegk
            END DO
          ELSE
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lha,Lhc,Na,Nb,dNa,dNb,Legj,Legk,dLegj,dLegk)
            DO l=1,nvec
              La = H1Basis_WedgeL(node1, u(l), v(l))
              Lb = H1Basis_WedgeL(node4, u(l), v(l))

              Lha = H1Basis_WedgeH(node1, w(l))
              Lhc = H1Basis_WedgeH(node2, w(l))

              Na = N(l,node1)
              Nb = N(l,node3)

              dNa = grad(l,node1,:)
              dNb = grad(l,node3,:)

              Legj = H1Basis_LegendreP(k-1,Lb-La)
              Legk = H1Basis_LegendreP(j,Lhc-Lha)

              dLegj = H1Basis_dLegendreP(k-1,Lb-La)*(dLb-dLa)
              dLegk = H1Basis_dLegendreP(j,Lhc-Lha)*(dLhc-dLha)

              grad(l,nbasis+k,:) = dNa*Nb*Legj*Legk + Na*dNb*Legj*Legk + &
                     Na*Nb*dLegj*Legk + Na*Nb*Legj*dLegk
            END DO
          END IF
        END DO
        nbasis = nbasis + MAX(pmax(i)-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dWedgeFaceP
  

  SUBROUTINE H1Basis_WedgeBubbleP(nvec, u, v, w, pmax, nbasismax, fval, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    
    INTEGER :: i, j, k, l
    REAL(KIND=dp) :: L1, L2, L3, s,t
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=0,pmax-3
      DO j=0,pmax-3-i
        DO k=1,pmax-1
          !_ELMER_OMP_SIMD PRIVATE(L1, L2, L3, s,t)
          DO l=1,nvec
            L1 = H1Basis_WedgeL(1,u(l),v(l))
            L2 = H1Basis_WedgeL(2,u(l),v(l))
            L3 = H1Basis_WedgeL(3,u(l),v(l))

            s = L2-L1
            t = 2*L3-1

            ! Get value of bubble function
            fval(l,nbasis+k) = L1*L2*L3*H1Basis_LegendreP(i,s)*&
                    H1Basis_LegendreP(j,t)*H1Basis_Phi(k+1,w(l))
          END DO
        END DO
        ! nbasis = nbasis + (pmax-3-i-j)-2+1
        nbasis = nbasis + MAX(pmax-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_WedgeBubbleP

  SUBROUTINE H1Basis_dWedgeBubbleP(nvec, u, v, w, pmax, nbasismax, grad, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    
    ! Parameters
    INTEGER :: i,j,k,l
    ! Variables
    REAL(KIND=dp) :: L1,L2,L3,Legi,Legj, Legk, dLegi(3), dLegj(3), dLegk(3), &
                  s,t,ds(3),dt(3), dL1(3),dL2(3),dL3(3)
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64      

    dL1 = H1Basis_dWedgeL(1)
    dL2 = H1Basis_dWedgeL(2)
    dL3 = H1Basis_dWedgeL(3)

    DO i=0,pmax-3
      DO j=0,pmax-3-i
        DO k=1,pmax-1
          !_ELMER_OMP_SIMD PRIVATE(L1, L2, L3, Legi, Legj, Legk, dLegi, dLegj, dLegk, s,t,ds,dt)
          DO l=1,nvec

            ! Values of function L
            L1 = H1Basis_WedgeL(1,u(l),v(l))
            L2 = H1Basis_WedgeL(2,u(l),v(l))
            L3 = H1Basis_WedgeL(3,u(l),v(l))
            
            s = L2-L1
            t = 2*L3-1
            ds = dL2-dL1
            dt = 2*dL3

            Legi = H1Basis_LegendreP(i,s)
            Legj = H1Basis_LegendreP(j,t)
            Legk = H1Basis_Phi(k+1,w(l))

            dLegi = H1Basis_dLegendreP(i,s)*ds
            dLegj = H1Basis_dLegendreP(j,t)*dt
            dLegk = H1Basis_dPhi(k+1,w(l))*[0,0,1]

            grad(l,nbasis+k,:) = dL1*L2*L3*Legi*Legj*Legk + L1*dL2*L3*Legi*Legj*Legk + &
                    L1*L2*dL3*Legi*Legj*Legk + L1*L2*L3*dLegi*Legj*Legk + &
                    L1*L2*L3*Legi*dLegj*Legk + L1*L2*L3*Legi*Legj*dLegk
          END DO
        END DO
        ! nbasis = nbasis + (pmax-3-i-j)-2+1
        nbasis = nbasis + MAX(pmax-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dWedgeBubbleP

!------------------------------------------------------------------------------
!>     Pyramid nodal basis at points (u,v,w)
!------------------------------------------------------------------------------
  SUBROUTINE H1Basis_PyramidNodalP(nvec, u, v, w, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    ! Result
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER :: j
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: s, sq2 = SQRT(2.0_dp)
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

      !_ELMER_OMP_SIMD
      DO j=1,nvec
        s = w(j) / sq2
        fval(j,nbasis+1) = (1-u(j)-v(j)-s+u(j)*v(j)/(1-s))/4
        fval(j,nbasis+2) = (1+u(j)-v(j)-s-u(j)*v(j)/(1-s))/4
        fval(j,nbasis+3) = (1+u(j)+v(j)-s+u(j)*v(j)/(1-s))/4
        fval(j,nbasis+4) = (1-u(j)+v(j)-s-u(j)*v(j)/(1-s))/4
        fval(j,nbasis+5) =  s
      END DO
      nbasis = nbasis + 5
    END SUBROUTINE H1Basis_PyramidNodalP


!------------------------------------------------------------------------------
!>     Gradient of pyramids nodal basis at point (u,v,w)
!------------------------------------------------------------------------------
  SUBROUTINE H1Basis_dPyramidNodalP(nvec, u, v, w, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    ! Result
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER :: j
    REAL(KIND=dp) :: s, sq2=SQRT(2.0_dp), sq42=4*SQRT(2._dp)
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64

      !_ELMER_OMP_SIMD
      DO j=1,nvec
        s = v(j)/(1-w(j)/sq2)
        grad(j,nbasis+1,1) = (-1 + s)/4
        grad(j,nbasis+2,1) = ( 1 - s)/4
        grad(j,nbasis+3,1) = ( 1 + s)/4
        grad(j,nbasis+4,1) = (-1 - s)/4
        grad(j,nbasis+5,1) = 0
      END DO

      !_ELMER_OMP_SIMD
      DO j=1,nvec
        s = u(j)/(1-w(j)/sq2)
        grad(j,nbasis+1,2) = (-1 + s)/4
        grad(j,nbasis+2,2) = (-1 - s)/4
        grad(j,nbasis+3,2) = ( 1 + s)/4
        grad(j,nbasis+4,2) = ( 1 - s)/4
        grad(j,nbasis+5,2) = 0
      END DO

      !_ELMER_OMP_SIMD
      DO j=1,nvec
        s = u(j)*v(j)/(1-w(j)/sq2)**2
        grad(j,nbasis+1,3) = (-1 + s)/sq42
        grad(j,nbasis+2,3) = (-1 - s)/sq42
        grad(j,nbasis+3,3) = (-1 + s)/sq42
        grad(j,nbasis+4,3) = (-1 - s)/sq42
        grad(j,nbasis+5,3) = 1/sq2
      END DO
      nbasis = nbasis + 5
  END SUBROUTINE H1Basis_dPyramidnodalP

    ! Define affine coordinates for pyramid square face
    FUNCTION H1Basis_PyramidL(which, u, v) RESULT(value)
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
    END FUNCTION H1Basis_PyramidL

    FUNCTION H1Basis_dPyramidL(which) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      ! Variables
      REAL(KIND=dp) :: grad(3)

      SELECT CASE (which)
      CASE (1)
         grad = [-1d0/2,-1d0/2, 0d0 ]
      CASE (2)
         grad = [ 1d0/2,-1d0/2, 0d0 ]
      CASE (3)
         grad = [ 1d0/2, 1d0/2, 0d0 ]
      CASE (4)
         grad = [-1d0/2, 1d0/2, 0d0 ]
      CASE DEFAULT
         CALL Fatal('PElementBase::dPyramidL','Unknown affine coordinate for square face')
      END SELECT
    END FUNCTION H1Basis_dPyramidL


    FUNCTION H1Basis_PyramidTL(which, u, v, w) RESULT(value)
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
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidTL','Unknown function L for brick')
      END SELECT
    END FUNCTION H1Basis_PyramidTL

    FUNCTION H1Basis_dPyramidTL(which) RESULT(grad)
      IMPLICIT NONE
      
      ! Parameters
      INTEGER, INTENT(IN) :: which
      ! Variables
      REAL(KIND=dp) :: grad(3),s
      
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
      CASE DEFAULT
         CALL Fatal('PElementBase::PyramidTL','Unknown function L for brick')
      END SELECT
      grad = grad/2
      grad(3) = grad(3)/SQRT(2._dp)
    END FUNCTION H1Basis_dPyramidTL



  SUBROUTINE H1Basis_PyramidEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: n
    REAL(KIND=dp) :: La, Lb, Na, Nb
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l, node1, node2, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=1,4 ! square face
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          Na = fval(k,node1)
          Nb = fval(k,node2)
          La = H1Basis_PyramidL(node1,u(k),v(k))
          Lb = H1Basis_PyramidL(node2,u(k),v(k))

          fval(k,nbasis+j) = Na*Nb*H1Basis_varPhi(j+1,Lb-La)
        END DO
      END DO
      nbasis = nbasis + MAX(pmax(i)-1,0)
    END DO

    DO i=5,8 ! triangular faces
      Node1 = edgedir(1,i)
      Node2 = edgedir(2,i)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          Na = fval(k,node1)
          Nb = fval(k,node2)
          La = H1Basis_PyramidTL(node1,u(k),v(k),w(k))
          Lb = H1Basis_PyramidTL(node2,u(k),v(k),w(k))
          fval(k,nbasis+j) = Na*Nb*H1Basis_varPhi(j+1,Lb-La)
        END DO
      END DO
      nbasis = nbasis + MAX(pmax(i)-1,0)
    END DO
  END SUBROUTINE H1Basis_PyramidEdgeP


  SUBROUTINE H1Basis_dPyramidEdgeP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: n
    REAL(KIND=dp) :: La, Lb, dLa(3), dLb(3), Na, Nb, dNa(3), dNb(3), Phi, dPhi(3)
    INTEGER :: i,j,k,l, node1, node2, nnb
!!!!DIR$ ASSUME_ALIGNED u:64, v:64, w:64

    REAL(KIND=dp) :: s, sq2=SQRT(2.0_dp)

    nnb=0; CALL H1Basis_PyramidNodalP(nvec, u, v, w, nbasismax, n, nnb)

    DO i=1,4 ! square face
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)

      dLa = H1Basis_dPyramidL(node1)
      dLb = H1Basis_dPyramidL(node2)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, dNa, dNb, Phi, dPhi)
        DO k=1,nvec
          Na = N(k,node1)
          Nb = N(k,node2)

          dNa = grad(k,node1,:)
          dNb = grad(k,node2,:)

          La = H1Basis_PyramidL(node1,u(k),v(k))
          Lb = H1Basis_PyramidL(node2,u(k),v(k))

          Phi = H1Basis_varPhi(j+1,Lb-La)
          dPhi = H1Basis_dvarPhi(j+1,Lb-La)*(dLb-dLa)

          grad(k,nbasis+j,:) = dNa*Nb*Phi + Na*dNb*Phi + Na*Nb*dPhi
        END DO
      END DO
      nbasis = nbasis + MAX(pmax(i)-1,0)
    END DO

    DO i=5,8 ! triangular faces
      Node1 = edgedir(1,i)
      Node2 = edgedir(2,i)

      dLa = H1Basis_dPyramidTL(node1)
      dLb = H1Basis_dPyramidTL(node2)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, dNa, dNb, Phi, dPhi)
        DO k=1,nvec
          Na = N(k,node1)
          Nb = N(k,node2)
          dNa = grad(k,node1,:)
          dNb = grad(k,node2,:)

          La = H1Basis_PyramidTL(node1,u(k),v(k),w(k))
          Lb = H1Basis_PyramidTL(node2,u(k),v(k),w(k))

          Phi = H1Basis_varPhi(j+1,Lb-La)
          dPhi = H1Basis_dvarPhi(j+1,Lb-La)*(dLb-dLa)

          grad(k,nbasis+j,:) = dNa*Nb*Phi + Na*dNb*Phi + Na*Nb*dPhi
        END DO
      END DO
      nbasis = nbasis + MAX(pmax(i)-1,0)
    END DO
  END SUBROUTINE H1Basis_dPyramidEdgeP



  SUBROUTINE H1Basis_PyramidFaceP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: n
    REAL(KIND=dp) :: La, Lb, Lc, Pa, Pb, Pc
    INTEGER :: i,j,k,l, node1, node2, node3, node4, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=1,1 ! square face
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      node4 = facedir(4,i)

      DO j=0,pmax(i)-2
        DO k=1,pmax(i)-1
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Pa, Pb)
          DO l=1,nvec
            Pa = fval(l,node1)
            Pb = fval(l,node3)

            La = H1Basis_PyramidL(node1,u(l),v(l))
            Lb = H1Basis_PyramidL(node2,u(l),v(l))
            Lc = H1Basis_PyramidL(node4,u(l),v(l))
            fval(l,nbasis+k) = Pa*Pb*H1Basis_LegendreP(j,Lb-La)*H1Basis_LegendreP(k-1,Lc-La)
          END DO
        END DO
        nbasis = nbasis + pmax(i)-1
      END DO
    END DO

    DO i=2,5 ! tri face
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-2-j
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Pa, Pb, Pc)
          DO l=1,nvec
            Pa = fval(l,node1)
            Pb = fval(l,node2)
            Pc = fval(l,node3)

            La = H1Basis_PyramidTL(node1,u(l),v(l),w(l))
            Lb = H1Basis_PyramidTL(node2,u(l),v(l),w(l))
            Lc = H1Basis_PyramidTL(node3,u(l),v(l),w(l))
            fval(l,nbasis+k) = Pa*Pb*Pc*H1Basis_LegendreP(j,Lb-La)*H1Basis_LegendreP(k-1,2*Lc-1)
          END DO
        END DO
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_PyramidFaceP


  SUBROUTINE H1Basis_dPyramidFaceP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: n
    REAL(KIND=dp) :: La, Lb, Lc, Pa, Pb, Pc, Legi, Legj
    REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dPa, dPb, dPc, dLegi, dLegj
    INTEGER :: i,j,k,l, node1, node2, node3, node4, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64

    nnb = 0; CALL H1Basis_PyramidNodalP(nvec,u,v,w,nbasismax,n,nnb)

    DO i=1,1 ! square face
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      node4 = facedir(4,i)

      dLa = H1Basis_dPyramidL(node1)
      dLb = H1Basis_dPyramidL(node2)
      dLc = H1Basis_dPyramidL(node4)

      DO j=0,pmax(i)-2
        DO k=1,pmax(i)-1
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Pa, Pb, dPa, dPb, Legi, Legj, dLegi, dLegj)
          DO l=1,nvec
            Pa = n(l,node1)
            Pb = n(l,node3)

            dPa = grad(l,node1,:)
            dPb = grad(l,node3,:)

            La = H1Basis_PyramidL(node1,u(l),v(l))
            Lb = H1Basis_PyramidL(node2,u(l),v(l))
            Lc = H1Basis_PyramidL(node4,u(l),v(l))

            Legi = H1Basis_LegendreP(j,Lb-La)
            Legj = H1Basis_LegendreP(k-1,Lc-La)

            dLegi = H1Basis_dLegendreP(j,Lb-La)*(dLb-dLa)
            dLegj = H1Basis_dLegendreP(k-1,Lc-La)*(dLc-dLa)
 
            grad(l,nbasis+k,:) = dPa*Pb*Legi*Legj + Pa*dPb*Legi*Legj + &
                       Pa*Pb*dLegi*Legj + Pa*Pb*Legi*dLegj
          END DO
        END DO
        nbasis = nbasis + pmax(i)-1
      END DO
    END DO

    DO i=2,5 ! tri face
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
 
      dLa = H1Basis_dPyramidTL(node1)
      dLb = H1Basis_dPyramidTL(node2)
      dLc = H1Basis_dPyramidTL(node3)

      DO j=0,pmax(i)-3
        DO k=1,pmax(i)-2-j
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Pa, Pb, Pc, dPa, dPb, dPc, Legi, Legj, dLegi, dLegj)
          DO l=1,nvec
            Pa = n(l,node1)
            Pb = n(l,node2)
            Pc = n(l,node3)

            dPa = grad(l,node1,:)
            dPb = grad(l,node2,:)
            dPc = grad(l,node3,:)

            La = H1Basis_PyramidTL(node1,u(l),v(l),w(l))
            Lb = H1Basis_PyramidTL(node2,u(l),v(l),w(l))
            Lc = H1Basis_PyramidTL(node3,u(l),v(l),w(l))

            Legi = H1Basis_LegendreP(j,Lb-La)
            Legj = H1Basis_LegendreP(k-1,2*Lc-1)

            dLegi = H1Basis_dLegendreP(j,Lb-La)*(dLb-dLa)
            dLegj = H1Basis_dLegendreP(k-1,2*Lc-1)*2*dLc
 
            grad(l,nbasis+k,:) = dPa*Pb*Pc*Legi*Legj + Pa*dPb*Pc*Legi*Legj + &
               Pa*Pb*dPc*Legi*Legj + Pa*Pb*Pc*dLegi*Legj + Pa*Pb*Pc*Legi*dLegj
          END DO
        END DO
        nbasis = nbasis + MAX(pmax(i)-j-2,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dPyramidFaceP

  SUBROUTINE H1Basis_PyramidBubbleP(nvec, u, v, w, pmax, nbasismax, fval, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis

    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: n
    REAL(KIND=dp) :: La, Lb, Lc, Pa, Pb, Pc, Legi, Legj, Legk, s
    REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dPa, dPb, dPc, dLegi, dLegj
    INTEGER :: i,j,k,l, node1, node2, node3, node4, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! Calculate value of function
    DO i=0,pmax-3
      DO j=0,pmax-i-3
        DO k=1,pmax-i-j-2
          nbasis = nbasis + 1
          !_ELMER_OMP_SIMD PRIVATE(Pa, Pb, Pc, Legi, Legj, Legk,s)
          DO l=1,nvec
             s = w(l) / SQRT(2._dp)
             Pa = fval(l,1)
             Pb = fval(l,3)
             Pc = fval(l,5)
             Legi = H1Basis_LegendreP(i,u(l))
             Legj = H1Basis_LegendreP(j,v(l))
             Legk = H1Basis_LegendreP(k-1,2*s-1)
             fval(l,nbasis) = Pa*Pb*Pc*Legi*Legj*Legk
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE H1Basis_PyramidBubbleP

  SUBROUTINE H1Basis_dPyramidBubbleP(nvec, u, v, w, pmax, nbasismax, grad, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis

    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax) :: n
    REAL(KIND=dp) :: La, Lb, Lc, Pa, Pb, Pc, Legi, Legj, Legk, s
    REAL(KIND=dp), DIMENSION(3) :: dLa, dLb, dLc, dPa, dPb, dPc, dLegi, dLegj, dLegk
    INTEGER :: i,j,k,l, node1, node2, node3, node4, nnb
    REAL(KIND=dp), PARAMETER :: a = 1/SQRT(2._dp)
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64


    nnb = 0; CALL H1Basis_PyramidNodalP(nvec,u,v,w,nbasismax,n,nnb)
    ! Calculate value of function
    DO i=0,pmax-3
      DO j=0,pmax-i-3
        DO k=1,pmax-i-j-2
          nbasis = nbasis + 1
          !_ELMER_OMP_SIMD PRIVATE(Pa, Pb, Pc, Legi, Legj, Legk, dPa, dPb, dPc, dLegi, dLegj, dLegk, s)
          DO l=1,nvec
             s = a * w(l)
             Pa = N(l,1)
             Pb = N(l,3)
             Pc = N(l,5)

             dPa = grad(l,1,:)
             dPb = grad(l,3,:)
             dPc = grad(l,5,:)

             Legi = H1Basis_LegendreP(i,u(l))
             Legj = H1Basis_LegendreP(j,v(l))
             Legk = H1Basis_LegendreP(k-1,2*s-1)

             dLegi = 0; dLegj=0; dLegk=0
             dLegi(1) = H1Basis_dLegendreP(i,u(l))
             dLegj(2) = H1Basis_dLegendreP(j,v(l))
             dLegk(3) = H1Basis_dLegendreP(k-1,2*s-1)*2*a

             grad(l,nbasis,:) = dPa*Pb*Pc*Legi*Legj*Legk + Pa*dPb*Pc*Legi*Legj*Legk + &
                                Pa*Pb*dPc*Legi*Legj*Legk + Pa*Pb*Pc*dLegi*Legj*Legk + &
                                Pa*Pb*Pc*Legi*dLegj*Legk + Pa*Pb*Pc*Legi*Legj*dLegk
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE H1Basis_dPyramidBubbleP


  SUBROUTINE H1Basis_BrickNodal(nvec, u, v, w, nbasismax, fval, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(KIND=dp) :: fval(VECTOR_BLOCK_LENGTH,nbasismax)
    INTEGER, INTENT(INOUT) :: nbasis
    
    REAL(Kind=dp), PARAMETER :: c = 1D0/8D0
    INTEGER :: j
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,nbasis+1) = c*(1-u(j))*(1-v(j))*(1-w(j))
      ! Node 2
      fval(j,nbasis+2) = c*(1+u(j))*(1-v(j))*(1-w(j))
      ! Node 3
      fval(j,nbasis+3) = c*(1+u(j))*(1+v(j))*(1-w(j))
      ! Node 4
      fval(j,nbasis+4) = c*(1-u(j))*(1+v(j))*(1-w(j))
      ! Node 5
      fval(j,nbasis+5) = c*(1-u(j))*(1-v(j))*(1+w(j))
      ! Node 6
      fval(j,nbasis+6) = c*(1+u(j))*(1-v(j))*(1+w(j))
      ! Node 7
      fval(j,nbasis+7) = c*(1+u(j))*(1+v(j))*(1+w(j))
      ! Node 8
      fval(j,nbasis+8) = c*(1-u(j))*(1+v(j))*(1+w(j))
    END DO
    nbasis = nbasis + 8
  END SUBROUTINE H1Basis_BrickNodal

  SUBROUTINE H1Basis_dBrickNodal(nvec, u, v, w, nbasismax, grad, nbasis)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp) :: grad(VECTOR_BLOCK_LENGTH,nbasismax,3)
    INTEGER, INTENT(INOUT) :: nbasis
    REAL(Kind=dp), PARAMETER :: c = 1D0/8D0
    INTEGER :: j
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) = -c*(1-v(j))*(1-w(j))
      grad(j,nbasis+2,1) =  c*(1-v(j))*(1-w(j))
      grad(j,nbasis+3,1) =  c*(1+v(j))*(1-w(j))
      grad(j,nbasis+4,1) = -c*(1+v(j))*(1-w(j))
      grad(j,nbasis+5,1) = -c*(1-v(j))*(1+w(j))
      grad(j,nbasis+6,1) =  c*(1-v(j))*(1+w(j))
      grad(j,nbasis+7,1) =  c*(1+v(j))*(1+w(j))
      grad(j,nbasis+8,1) = -c*(1+v(j))*(1+w(j))
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = -c*(1-u(j))*(1-w(j))
      grad(j,nbasis+2,2) = -c*(1+u(j))*(1-w(j))
      grad(j,nbasis+3,2) =  c*(1+u(j))*(1-w(j))
      grad(j,nbasis+4,2) =  c*(1-u(j))*(1-w(j))
      grad(j,nbasis+5,2) = -c*(1-u(j))*(1+w(j))
      grad(j,nbasis+6,2) = -c*(1+u(j))*(1+w(j))
      grad(j,nbasis+7,2) =  c*(1+u(j))*(1+w(j))
      grad(j,nbasis+8,2) =  c*(1-u(j))*(1+w(j))
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,3) = -c*(1-u(j))*(1-v(j))
      grad(j,nbasis+2,3) = -c*(1+u(j))*(1-v(j))
      grad(j,nbasis+3,3) = -c*(1+u(j))*(1+v(j))
      grad(j,nbasis+4,3) = -c*(1-u(j))*(1+v(j))
      grad(j,nbasis+5,3) =  c*(1-u(j))*(1-v(j))
      grad(j,nbasis+6,3) =  c*(1+u(j))*(1-v(j))
      grad(j,nbasis+7,3) =  c*(1+u(j))*(1+v(j))
      grad(j,nbasis+8,3) =  c*(1-u(j))*(1+v(j))
    END DO
        nbasis = nbasis + 8
  END SUBROUTINE H1Basis_dBrickNodal

  FUNCTION H1Basis_BrickL(node, u, v, w) RESULT(fval)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: node
    REAL(KIND=dp), INTENT(IN) :: u, v, w
    REAL(KIND=dp) :: fval
    REAL(KIND=dp), PARAMETER :: c = 1/2D0
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) _ELMER_LINEAR_REF(u) &
    !_ELMER_OMP _ELMER_LINEAR_REF(v) _ELMER_LINEAR_REF(w) NOTINBRANCH
    
    SELECT CASE (node)
    CASE (1)
      fval = c*(3-u-v-w)
    CASE (2)
      fval = c*(3+u-v-w)
    CASE (3)
      fval = c*(3+u+v-w)
    CASE (4)
      fval = c*(3-u+v-w)
    CASE (5)
      fval = c*(3-u-v+w)
    CASE (6)
      fval = c*(3+u-v+w)
    CASE (7)
      fval = c*(3+u+v+w)
    CASE (8)
      fval = c*(3-u+v+w)
    END SELECT
  END FUNCTION H1Basis_BrickL

  FUNCTION H1Basis_dBrickL(node) RESULT(grad)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: node
    ! REAL(KIND=dp), INTENT(IN) :: u, v
    REAL(KIND=dp) :: grad(3)
    REAL(KIND=dp), PARAMETER :: c = 1/2D0
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE (1)
      grad(1:3) = c*[-1,-1,-1 ]
    CASE (2)
      grad(1:3) = c*[ 1,-1,-1 ]
    CASE (3)
      grad(1:3) = c*[ 1, 1,-1 ]
    CASE (4)
      grad(1:3) = c*[-1, 1,-1 ]
    CASE (5)
      grad(1:3) = c*[-1,-1, 1 ]
    CASE (6)
      grad(1:3) = c*[ 1,-1, 1 ]
    CASE (7)
      grad(1:3) = c*[ 1, 1, 1 ]
    CASE (8)
      grad(1:3) = c*[-1, 1, 1 ]
    END SELECT
  END FUNCTION H1Basis_dBrickL

  SUBROUTINE H1Basis_BrickEdgeL(edge, u, v, w, La, Lb)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: edge
    REAL(KIND=dp), INTENT(IN) :: u, v, w
    REAL(KIND=dp), INTENT(OUT) :: La, Lb
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(edge) _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) &
    !_ELMER_OMP _ELMER_LINEAR_REF(w) _ELMER_LINEAR_REF(La) _ELMER_LINEAR_REF(Lb) NOTINBRANCH
    
    SELECT CASE(edge)
    CASE (1)
      La=1-v
      Lb=1-w
    CASE (2)
      La=1+u
      Lb=1-w
    CASE (3)
      La=1+v
      Lb=1-w
    CASE (4)
      La=1-u
      Lb=1-w
    CASE (5)
      La=1-v
      Lb=1+w
    CASE (6)
      La=1+u
      Lb=1+w
    CASE (7)
      La=1+v
      Lb=1+w
    CASE (8)
      La=1-u
      Lb=1+w
    CASE (9)
      La=1-u
      Lb=1-v
    CASE (10)
      La=1+u
      Lb=1-v
    CASE (11)
      La=1+u
      Lb=1+v
    CASE (12)
      La=1-u
      Lb=1+v
    END SELECT
  END SUBROUTINE H1Basis_BrickEdgeL

  SUBROUTINE H1Basis_dBrickEdgeL(edge, dLa, dLb)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: edge
    REAL(KIND=dp), INTENT(OUT) :: dLa(3), dLb(3)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(edge) NOTINBRANCH
    
    SELECT CASE(edge)
    CASE (1)
      dLa=[ 0,-1, 0 ] ! 1-v
      dLb=[ 0, 0,-1 ] ! 1-w
    CASE (2)
      dLa=[ 1, 0, 0 ] ! 1+u
      dLb=[ 0, 0,-1 ] ! 1-w
    CASE (3)
      dLa=[ 0, 1, 0 ] ! 1+v
      dLb=[ 0, 0,-1 ] ! 1-w
    CASE (4)
      dLa=[-1, 0, 0 ] ! 1-u
      dLb=[ 0, 0,-1 ] ! 1-w
    CASE (5)
      dLa=[ 0,-1, 0 ] ! 1-v
      dLb=[ 0, 0, 1 ] ! 1+w
    CASE (6)
      dLa=[ 1, 0, 0 ] ! 1+u
      dLb=[ 0, 0, 1 ] ! 1+w
    CASE (7)
      dLa=[ 0, 1, 0 ] ! 1+v
      dLb=[ 0, 0, 1 ] ! 1+w
    CASE (8)
      dLa=[-1, 0, 0 ] ! 1-u
      dLb=[ 0, 0, 1 ] ! 1+w
    CASE (9)
      dLa=[-1, 0, 0 ] ! 1-u
      dLb=[ 0,-1, 0 ] ! 1-v
    CASE (10)
      dLa=[ 1, 0, 0 ] ! 1+u
      dLb=[ 0,-1, 0 ] ! 1-v
    CASE (11)
      dLa=[ 1, 0, 0 ] ! 1+u
      dLb=[ 0, 1, 0 ] ! 1+v
    CASE (12)
      dLa=[-1, 0, 0 ] ! 1-u
      dLb=[ 0, 1, 0 ] ! 1+v
    END SELECT
  END SUBROUTINE H1Basis_dBrickEdgeL

! --- start serendipity brick

  
  SUBROUTINE H1Basis_SD_BrickEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp), PARAMETER :: c = 1/4D0
    REAL(KIND=dp) :: La, Lb, Aa, Ba
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each edge
    DO i=1,12
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Aa, Ba)
        DO k=1,nvec
          La=H1Basis_BrickL(edgedir(1,i), u(k), v(k), w(k))
          Lb=H1Basis_BrickL(edgedir(2,i), u(k), v(k), w(k))
          CALL H1Basis_BrickEdgeL(i, u(k), v(k), w(k), Aa, Ba)
          fval(k, nbasis+j) = c*H1Basis_Phi(j+1, Lb-La)*Aa*Ba
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_SD_BrickEdgeP
  
  SUBROUTINE H1Basis_SD_dBrickEdgeP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp), PARAMETER :: c = 1/4D0
    REAL(KIND=dp) :: La, Lb, Aa, Ba, Phi, dPhi, dLa(3), &
            dLb(3), dAa(3), dBa(3)
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each edge
    DO i=1,12
      dLa = H1Basis_dBrickL(edgedir(1,i))
      dLb = H1Basis_dBrickL(edgedir(2,i))
      CALL H1Basis_dBrickEdgeL(i, dAa, dBa)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Aa, Ba, Phi, dPhi)
        DO k=1,nvec
          La=H1Basis_BrickL(edgedir(1,i), u(k), v(k), w(k))
          Lb=H1Basis_BrickL(edgedir(2,i), u(k), v(k), w(k))
          CALL H1Basis_BrickEdgeL(i, u(k), v(k), w(k), Aa, Ba)

          Phi = H1Basis_Phi(j+1, Lb-La)
          dPhi = H1Basis_dPhi(j+1, Lb-La)

          grad(k,nbasis+j,1) = c*dPhi*(dLb(1)-dLa(1))*Aa*Ba +&
                  c*Phi*dAa(1)*Ba + c*Phi*Aa*dBa(1)
          grad(k,nbasis+j,2) = c*dPhi*(dLb(2)-dLa(2))*Aa*Ba +&
                  c*Phi*dAa(2)*Ba + c*Phi*Aa*dBa(2)
          grad(k,nbasis+j,3) = c*dPhi*(dLb(3)-dLa(3))*Aa*Ba +&
                  c*Phi*dAa(3)*Ba + c*Phi*Aa*dBa(3)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_SD_dBrickEdgeP
  
  SUBROUTINE H1Basis_SD_BrickFaceP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Ld
    REAL(KIND=dp), PARAMETER :: a = 1D0/4
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each face
    DO i=1,6
      DO j=2,pmax(i)
        DO k=1,pmax(i)-j-1
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Ld)
          DO l=1,nvec
            La=H1Basis_BrickL(facedir(1,i), u(l), v(l), w(l))
            Lb=H1Basis_BrickL(facedir(2,i), u(l), v(l), w(l))
            Lc=H1Basis_BrickL(facedir(3,i), u(l), v(l), w(l))
            Ld=H1Basis_BrickL(facedir(4,i), u(l), v(l), w(l))
            
            fval(l, nbasis+k) = (a*(La+Lb+Lc+Ld)-1)*H1Basis_Phi(j, Lb-La) &
                               *H1Basis_Phi(k+1, Ld-La)
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_BrickFaceP
  
  SUBROUTINE H1Basis_SD_dBrickFaceP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Ld, dLa(3), dLb(3), dLc(3), dLd(3), &
          PhiLbLa, PhiLdLa, LaLbLcLd
    REAL(KIND=dp), PARAMETER :: a=1D0/4D0
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each face
    DO i=1,6
      dLa=H1Basis_dBrickL(facedir(1,i))
      dLb=H1Basis_dBrickL(facedir(2,i))
      dLc=H1Basis_dBrickL(facedir(3,i))
      dLd=H1Basis_dBrickL(facedir(4,i))
      DO j=2,pmax(i)
        DO k=1,pmax(i)-j-1
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Ld, PhiLbLa, PhiLdLa, LaLbLcLd)
          DO l=1,nvec
            La=H1Basis_BrickL(facedir(1,i), u(l), v(l), w(l))
            Lb=H1Basis_BrickL(facedir(2,i), u(l), v(l), w(l))
            Lc=H1Basis_BrickL(facedir(3,i), u(l), v(l), w(l))
            Ld=H1Basis_BrickL(facedir(4,i), u(l), v(l), w(l))
            
            PhiLbLa=H1Basis_Phi(j, Lb-La)
            PhiLdLa=H1Basis_Phi(k+1, Ld-La)
            LaLbLcLd=a*(La+Lb+Lc+Ld)-1

            grad(l, nbasis+k, 1)=a*(dLa(1)+dLb(1)+dLc(1)+dLd(1))*PhiLbLa*PhiLdLa + &
                  LaLbLcLd*H1Basis_dPhi(j,Lb-La)*(dLb(1)-dLa(1))*PhiLdLa + &
                  LaLbLcLd*PhiLbLa*H1Basis_dPhi(k+1, Ld-La)*(dLd(1)-dLa(1))
            grad(l, nbasis+k, 2)=a*(dLa(2)+dLb(2)+dLc(2)+dLd(2))*PhiLbLa*PhiLdLa + &
                  LaLbLcLd*H1Basis_dPhi(j,Lb-La)*(dLb(2)-dLa(2))*PhiLdLa + &
                  LaLbLcLd*PhiLbLa*H1Basis_dPhi(k+1, Ld-La)*(dLd(2)-dLa(2))
            grad(l, nbasis+k, 3)=a*(dLa(3)+dLb(3)+dLc(3)+dLd(3))*PhiLbLa*PhiLdLa + &
                  LaLbLcLd*H1Basis_dPhi(j,Lb-La)*(dLb(3)-dLa(3))*PhiLdLa + &
                  LaLbLcLd*PhiLbLa*H1Basis_dPhi(k+1, Ld-La)*(dLd(3)-dLa(3))
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_dBrickFaceP
  
  SUBROUTINE H1Basis_SD_BrickBubbleP(nvec, u, v, w, pmax, nbasismax, fval, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis

    INTEGER :: i, j, k, l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each bubble calculate value of basis function and its derivative
    ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,pmax
    DO i=2,pmax-4
      DO j=2,pmax-i-2
        DO k=1,pmax-i-j-1
          !_ELMER_OMP_SIMD
          DO l=1,nvec
            fval(l,nbasis+k) = H1Basis_Phi(i,u(l))*H1Basis_Phi(j,v(l))*H1Basis_Phi(k+1,w(l))
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_BrickBubbleP

  SUBROUTINE H1Basis_SD_dBrickBubbleP(nvec, u, v, w, pmax, nbasismax, grad, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis

    ! Variables
    INTEGER :: i, j, k, l
    REAL(KIND=dp) :: phiU, phiV, phiW
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each bubble calculate value of basis function and its derivative
    ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,pmax
    DO i=2,pmax-4
      DO j=2,pmax-i-2
        DO k=1,pmax-i-j-1
          !_ELMER_OMP_SIMD PRIVATE(phiU, phiV, phiW)
          DO l=1,nvec
            phiU = H1Basis_Phi(i,u(l))
            phiV = H1Basis_Phi(j,v(l))
            phiW = H1Basis_Phi(k+1,w(l))
            grad(l,nbasis+k,1) = H1Basis_dPhi(i,u(l))*phiV*phiW
            grad(l,nbasis+k,2) = phiU*H1Basis_dPhi(j,v(l))*phiW
            grad(l,nbasis+k,3) = phiU*phiV*H1Basis_dPhi(k+1,w(l))
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_SD_dBrickBubbleP


! --- end serendipity brick
  
  SUBROUTINE H1Basis_BrickEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir

    REAL(KIND=dp) :: La, Lb, Na, Nb, Lp(nvec,8)
    INTEGER :: i,j,k,l, node1, node2, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=1,8
      !_ELMER_OMP_SIMD
      DO k=1,nvec
        Lp(k,i) = H1Basis_BrickL(i, u(k), v(k), w(k))
      END DO
    END DO

    ! For each edge
    DO i=1,12
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)
      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          Na = fval(k, node1)
          Nb = fval(k, node2)
          La = Lp(k, node1)
          Lb = Lp(k, node2)
          fval(k, nbasis+j) = Na*Nb*H1Basis_varPhi(j+1, Lb-La)
        END DO
      END DO
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_BrickEdgeP
  
  SUBROUTINE H1Basis_dBrickEdgeP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, edgedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: edgedir
    REAL(KIND=dp) :: La, Lb, Na, Nb, Phi, dPhi, dLa(3), Lp(nvec,8), &
            dLb(3), dNa(3), dNb(3), N(VECTOR_BLOCK_LENGTH, nbasismax)
    INTEGER :: i,j,k, node1, node2, nnb
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    nnb = 0; CALL H1Basis_BrickNodal(nvec, u, v, w, nbasismax, n, nnb )

    DO i=1,8
      !_ELMER_OMP_SIMD
      DO k=1,nvec
        Lp(k,i) = H1Basis_BrickL(i, u(k), v(k), w(k))
      END DO
    END DO

    ! For each edge
    DO i=1,12
      node1 = edgedir(1,i)
      node2 = edgedir(2,i)

      dLa = H1Basis_dBrickL(node1)
      dLb = H1Basis_dBrickL(node2)

      DO j=1,pmax(i)-1
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, Phi, dPhi)
        DO k=1,nvec
          La = Lp(k,node1)
          Lb = Lp(k,node2)

          Na = N(k,node1)
          Nb = N(k,node2)
          dNa = grad(k,node1,:)
          dNb = grad(k,node2,:)

          Phi = H1Basis_varPhi(j+1, Lb-La)
          dPhi = H1Basis_dvarPhi(j+1, Lb-La)

          grad(k,nbasis+j,1) = dPhi*(dLb(1)-dLa(1))*Na*Nb + &
                  Phi*dNa(1)*Nb + Phi*Na*dNb(1)
          grad(k,nbasis+j,2) = dPhi*(dLb(2)-dLa(2))*Na*Nb + &
                  Phi*dNa(2)*Nb + Phi*Na*dNb(2)
          grad(k,nbasis+j,3) = dPhi*(dLb(3)-dLa(3))*Na*Nb + &
                  Phi*dNa(3)*Nb + Phi*Na*dNb(3)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
  END SUBROUTINE H1Basis_dBrickEdgeP
  
  SUBROUTINE H1Basis_BrickFaceP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Ld, Na, Nb, Lp(nvec,8)
    INTEGER :: i,j,k,l, node1, node2, node3, node4
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=1,8
      !_ELMER_OMP_SIMD
      DO k=1,nvec
        Lp(k,i) = H1Basis_BrickL(i, u(k), v(k), w(k))
      END DO
    END DO

    ! For each face
    DO i=1,6
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      node4 = facedir(4,i)
      DO j=0,pmax(i)-2
        DO k=1,pmax(i)-1
          !_ELMER_OMP_SIMD PRIVATE(Na,Nb,La, Lb, Ld)
          DO l=1,nvec
            Na = fval(l,node1)
            Nb = fval(l,node3)
            La = Lp(l,node1)
            Lb = Lp(l,node2)
            Ld = Lp(l,node4)
            fval(l, nbasis+k) = Na*Nb*H1Basis_LegendreP(j, Lb-La)*H1Basis_LegendreP(k-1, Ld-La)
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + pmax(i)-1
      END DO
    END DO
  END SUBROUTINE H1Basis_BrickFaceP
  
  SUBROUTINE H1Basis_dBrickFaceP(nvec, u, v, w, pmax, nbasismax, grad, nbasis, facedir)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, DIMENSION(:) CONTIG, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis
    INTEGER, DIMENSION(:,:) CONTIG, INTENT(IN) :: facedir

    REAL(KIND=dp) :: La, Lb, Lc, Ld, dLa(3), dLb(3), dLc(3), dLd(3), Lp(nvec,8), &
          Na,Nb,dNa(3),dNb(3), N(VECTOR_BLOCK_LENGTH, nbasismax), PhiU, PhiV, dPhiU, dPhiV
    INTEGER :: i,j,k,l, nnb, node1, node2, node3, node4
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    nnb = 0; CALL H1Basis_BrickNodal(nvec, u, v, w, nbasismax, n, nnb )
    DO i=1,8
      !_ELMER_OMP_SIMD
      DO k=1,nvec
        Lp(k,i) = H1Basis_BrickL(i, u(k), v(k), w(k))
      END DO
    END DO

    ! For each face
    DO i=1,6
      node1 = facedir(1,i)
      node2 = facedir(2,i)
      node3 = facedir(3,i)
      node4 = facedir(4,i)
      dLa = H1Basis_dBrickL(node1)
      dLb = H1Basis_dBrickL(node2)
      dLc = H1Basis_dBrickL(node3)
      dLd = H1Basis_dBrickL(node4)

      DO j=0,pmax(i)-2
        DO k=1,pmax(i)-1
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Ld, PhiU, PhiV, Na, Nb, dNa, dNb)
          DO l=1,nvec

            La = Lp(l,node1)
            Lb = Lp(l,node2)
            Ld = Lp(l,node4)
            
            PhiU = H1Basis_LegendreP(j, Lb-La)
            PhiV = H1Basis_LegendreP(k-1, Ld-La)

            dPhiU = H1Basis_dLegendreP(j, Lb-La)
            dPhiV = H1Basis_dLegendreP(k-1, Ld-La)

            Na = N(l,node1)
            Nb = N(l,node3)

            dNa = grad(l,node1,:)
            dNb = grad(l,node3,:)

            grad(l, nbasis+k, :) = dNa*Nb*PhiU*PhiV + Na*dNb*PhiU*PhiV + &
                Na*Nb*dPhiU*(dLb-dLa)*PhiV + Na*Nb*PhiU*dPhiV*(dLd-dLa)
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + pmax(i)-1
      END DO
    END DO
  END SUBROUTINE H1Basis_dBrickFaceP
  
  SUBROUTINE H1Basis_BrickBubbleP(nvec, u, v, w, pmax, nbasismax, fval, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER, INTENT(INOUT) :: nbasis

    INTEGER :: i, j, k, l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each bubble calculate value of basis function and its derivative
    ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,pmax
    DO i=2,pmax
      DO j=2,pmax
        DO k=1,pmax-1
          !_ELMER_OMP_SIMD
          DO l=1,nvec
            fval(l,nbasis+k) = H1Basis_Phi(i,u(l))*H1Basis_Phi(j,v(l))*H1Basis_Phi(k+1,w(l))
          END DO
        END DO
        nbasis = nbasis + pmax-1
      END DO
    END DO
  END SUBROUTINE H1Basis_BrickBubbleP

  SUBROUTINE H1Basis_dBrickBubbleP(nvec, u, v, w, pmax, nbasismax, grad, nbasis)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nvec
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u, v, w
    INTEGER, INTENT(IN) :: pmax
    INTEGER, INTENT(IN) :: nbasismax
    REAL(KIND=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER, INTENT(INOUT) :: nbasis

    ! Variables
    INTEGER :: i, j, k, l
    REAL(KIND=dp) :: phiU, phiV, phiW
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! For each bubble calculate value of basis function and its derivative
    ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,pmax
    DO i=2,pmax
      DO j=2,pmax
        DO k=1,pmax-1
          !_ELMER_OMP_SIMD PRIVATE(phiU, phiV, phiW)
          DO l=1,nvec
            phiU = H1Basis_Phi(i,u(l))
            phiV = H1Basis_Phi(j,v(l))
            phiW = H1Basis_Phi(k+1,w(l))

            grad(l,nbasis+k,1) = H1Basis_dPhi(i,u(l))*phiV*phiW
            grad(l,nbasis+k,2) = phiU*H1Basis_dPhi(j,v(l))*phiW
            grad(l,nbasis+k,3) = phiU*phiV*H1Basis_dPhi(k+1,w(l))
          END DO
        END DO
        nbasis = nbasis+pmax-1
      END DO
    END DO
  END SUBROUTINE H1Basis_dBrickBubbleP

  PURE FUNCTION H1Basis_Phi(k, x) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: x
    ! Return value
    REAL (KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(k) _ELMER_LINEAR_REF(x) NOTINBRANCH
    
    ! Phi function values (autogenerated to Horner form)
    SELECT CASE(k)
    CASE(2)
      fval = -SQRT(0.6D1) / 0.4D1 + SQRT(0.6D1) * x ** 2 / 0.4D1
    CASE(3)
      fval = (-SQRT(0.10D2) / 0.4D1 + SQRT(0.10D2) * x ** 2 / 0.4D1) * x
    CASE(4)
      fval = SQRT(0.14D2) / 0.16D2 + (-0.3D1 / 0.8D1 * SQRT(0.14D2) &
              + 0.5D1 / 0.16D2 * SQRT(0.14D2) * x ** 2) * x ** 2
    CASE(5)
      fval = (0.9D1 / 0.16D2 * SQRT(0.2D1) + (-0.15D2 / 0.8D1 * SQRT(0.2D1)&
              + 0.21D2 / 0.16D2 * SQRT(0.2D1) * x ** 2) * x ** 2) * x
    CASE(6)
      fval = -SQRT(0.22D2) / 0.32D2 + (0.15D2 / 0.32D2 * SQRT(0.22D2) + &
              (-0.35D2 / 0.32D2 * SQRT(0.22D2) + 0.21D2 / 0.32D2 * &
              SQRT(0.22D2) * x ** 2) * x ** 2) * x ** 2
    CASE(7)
      fval = (-0.5D1 / 0.32D2 * SQRT(0.26D2) + (0.35D2 / 0.32D2 * SQRT(0.26D2) &
              + (-0.63D2 / 0.32D2 * SQRT(0.26D2) + 0.33D2 / 0.32D2 * &
              SQRT(0.26D2) * x ** 2) * x ** 2) * x ** 2) * x
    CASE(8)
      fval = 0.5D1 / 0.256D3 * SQRT(0.30D2) + (-0.35D2 / 0.64D2 * SQRT(0.30D2) &
              + (0.315D3 / 0.128D3 * SQRT(0.30D2) + (-0.231D3 / 0.64D2 * SQRT(0.30D2) &
              + 0.429D3 / 0.256D3 * SQRT(0.30D2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(9)
      fval = (0.35D2 / 0.256D3 * SQRT(0.34D2) + (-0.105D3 / 0.64D2 * SQRT(0.34D2) +&
              (0.693D3 / 0.128D3 * SQRT(0.34D2) + (-0.429D3 / 0.64D2 * SQRT(0.34D2) +&
              0.715D3 / 0.256D3 * SQRT(0.34D2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
    CASE(10)
      fval = -0.7D1 / 0.512D3 * SQRT(0.38D2) + (0.315D3 / 0.512D3 * SQRT(0.38D2) + &
              (-0.1155D4 / 0.256D3 * SQRT(0.38D2) + (0.3003D4 / 0.256D3 * SQRT(0.38D2) +&
              (-0.6435D4 / 0.512D3 * SQRT(0.38D2) + 0.2431D4 / 0.512D3 * SQRT(0.38D2) &
              * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(11)
      fval = (-0.63D2 / 0.512D3 * SQRT(0.42D2) + (0.1155D4 / 0.512D3 * SQRT(0.42D2) + &
              (-0.3003D4 / 0.256D3 * SQRT(0.42D2) + (0.6435D4 / 0.256D3 * SQRT(0.42D2) +&
              (-0.12155D5 / 0.512D3 * SQRT(0.42D2) + 0.4199D4 / 0.512D3 * &
              SQRT(0.42D2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
    CASE(12)
      fval = 0.21D2 / 0.2048D4 * SQRT(0.46D2) + (-0.693D3 / 0.1024D4 * SQRT(0.46D2) + &
              (0.15015D5 / 0.2048D4 * SQRT(0.46D2) + (-0.15015D5 / 0.512D3 * SQRT(0.46D2) +&
              (0.109395D6 / 0.2048D4 * SQRT(0.46D2) + (-0.46189D5 / 0.1024D4 * SQRT(0.46D2) +&
              0.29393D5 / 0.2048D4 * SQRT(0.46D2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) &
              * x ** 2) * x ** 2
    CASE(13)
      fval = (0.1155D4 / 0.2048D4 * SQRT(0.2D1) + (-0.15015D5 / 0.1024D4 * SQRT(0.2D1) +&
              (0.225225D6 / 0.2048D4 * SQRT(0.2D1) + (-0.182325D6 / 0.512D3 * SQRT(0.2D1) +&
              (0.1154725D7 / 0.2048D4 * SQRT(0.2D1) + (-0.440895D6 / 0.1024D4 * SQRT(0.2D1) +&
              0.260015D6 / 0.2048D4 * SQRT(0.2D1) * x ** 2) * x ** 2) * x ** 2) * x ** 2) &
              * x ** 2) * x ** 2) * x
    CASE(14)
      fval = -0.99D2 / 0.4096D4 * SQRT(0.6D1) + (0.9009D4 / 0.4096D4 * SQRT(0.6D1) + &
              (-0.135135D6 / 0.4096D4 * SQRT(0.6D1) + (0.765765D6 / 0.4096D4 * SQRT(0.6D1) + &
              (-0.2078505D7 / 0.4096D4 * SQRT(0.6D1) + (0.2909907D7 / 0.4096D4 * SQRT(0.6D1) +&
              (-0.2028117D7 / 0.4096D4 * SQRT(0.6D1) + 0.557175D6 / 0.4096D4 * SQRT(0.6D1) &
              * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_Phi

  PURE FUNCTION H1Basis_dPhi(k, x) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: x
    ! Return value
    REAL (KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(k) _ELMER_LINEAR_REF(x) NOTINBRANCH
    
    ! Phi function values (autogenerated to Horner form)
    SELECT CASE(k)
    CASE(2)
      fval = x * SQRT(0.6D1) / 0.2D1
    CASE(3)
      fval = DBLE(3 * x ** 2 - 1) * SQRT(0.10D2) / 0.4D1
    CASE(4)
      fval = DBLE(5 * x ** 2 - 3) * DBLE(x) * SQRT(0.14D2) / 0.4D1
    CASE(5)
      fval = 0.3D1 / 0.16D2 * DBLE(3 + (35 * x ** 2 - 30) * x ** 2) * SQRT(0.2D1)
    CASE(6)
      fval = DBLE(15 + (63 * x ** 2 - 70) * x ** 2) * DBLE(x) * SQRT(0.22D2) / 0.16D2
    CASE(7)
      fval = DBLE(-5 + (105 + (231 * x ** 2 - 315) * x ** 2) * x ** 2) * &
              SQRT(0.26D2) / 0.32D2
    CASE(8)
      fval = DBLE(-35 + (315 + (429 * x ** 2 - 693) * x ** 2) * x ** 2) * &
              DBLE(x) * SQRT(0.30D2) / 0.32D2
    CASE(9)
      fval = DBLE(35 + (-1260 + (6930 + (6435 * x ** 2 - 12012) * x ** 2) &
              * x ** 2) * x ** 2) * SQRT(0.34D2) / 0.256D3
    CASE(10)
      fval = DBLE(315 + (-4620 + (18018 + (12155 * x ** 2 - 25740) * x ** 2) * &
              x ** 2) * x ** 2) * DBLE(x) * SQRT(0.38D2) / 0.256D3
    CASE(11)
      fval = DBLE(-63 + (3465 + (-30030 + (90090 + (46189 * x ** 2 - 109395) * &
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * SQRT(0.42D2) / 0.512D3
    CASE(12)
      fval = DBLE(-693 + (15015 + (-90090 + (218790 + (88179 * x ** 2 - 230945) *&
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * DBLE(x) * SQRT(0.46D2) / 0.512D3
    CASE(13)
      fval = 0.5D1 / 0.2048D4 * DBLE(231 + (-18018 + (225225 + (-1021020 + &
              (2078505 + (676039 * x ** 2 - 1939938) * x ** 2) * x ** 2) * &
              x ** 2) * x ** 2) * x ** 2) * SQRT(0.2D1)
    CASE(14)
      fval = 0.3D1 / 0.2048D4 * DBLE(3003 + (-90090 + (765765 + (-2771340 + &
              (4849845 + (1300075 * x ** 2 - 4056234) * x ** 2) * x ** 2) * &
              x ** 2) * x ** 2) * x ** 2) * DBLE(x) * SQRT(0.6D1)
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_dPhi

  PURE FUNCTION H1Basis_VarPhi(k, x) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: x
    ! Return value
    REAL (KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(k) _ELMER_LINEAR_REF(x) NOTINBRANCH

    SELECT CASE(k)
    CASE(2)
      fval = -SQRT(0.6D1)
    CASE(3)
      fval = -x * SQRT(0.10D2)
    CASE(4)
      fval = -DBLE(5 * x ** 2 - 1) * SQRT(0.14D2) / 0.4D1
    CASE(5)
      fval = -0.3D1 / 0.4D1 * DBLE(7 * x ** 2 - 3) * DBLE(x) * SQRT(0.2D1)
    CASE(6)
      fval = -DBLE(1 + (21 * x ** 2 - 14) * x ** 2) * SQRT(0.22D2) / 0.8D1
    CASE(7)
      fval = -DBLE(5 + (33 * x ** 2 - 30) * x ** 2) * DBLE(x) * &
              SQRT(0.26D2) / 0.8D1
    CASE(8)
      fval = -DBLE(-5 + (135 + (429 * x ** 2 - 495) * x ** 2) * x ** 2) *&
              SQRT(0.30D2) / 0.64D2
    CASE(9)
      fval = -DBLE(-35 + (385 + (715 * x ** 2 - 1001) * x ** 2) * x ** 2) * &
              DBLE(x) * SQRT(0.34D2) / 0.64D2
    CASE(10)
      fval = -DBLE(7 + (-308 + (2002 + (2431 * x ** 2 - 4004) * x ** 2) * &
              x ** 2) * x ** 2) * SQRT(0.38D2) / 0.128D3
    CASE(11)
      fval = -DBLE(63 + (-1092 + (4914 + (4199 * x ** 2 - 7956) * x ** 2) *&
              x ** 2) * x ** 2) * DBLE(x) * SQRT(0.42D2) / 0.128D3
    CASE(12)
      fval = -DBLE(-21 + (1365 + (-13650 + (46410 + (29393 * x ** 2 - 62985) *&
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * SQRT(0.46D2) / 0.512D3
    CASE(13)
      fval = -0.5D1 / 0.512D3 * DBLE(-231 + (5775 + (-39270 + (106590 + &
              (52003 * x ** 2 - 124355) * x ** 2) * x ** 2) * x ** 2) * &
              x ** 2) * DBLE(x) * SQRT(0.2D1)
    CASE(14)
      fval = -0.3D1 / 0.1024D4 * DBLE(33 + (-2970 + (42075 + (-213180 + &
              (479655 + (185725 * x ** 2 - 490314) * x ** 2) * x ** 2) *&
              x ** 2) * x ** 2) * x ** 2) * SQRT(0.6D1)
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_VarPhi

  PURE FUNCTION H1Basis_dVarPhi(k, x) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: x
    ! Return value
    REAL (KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(k) _ELMER_LINEAR_REF(x) NOTINBRANCH
    
    SELECT CASE(k)
    CASE(2)
      fval = 0
    CASE(3)
      fval = -SQRT(0.10D2)
    CASE(4)
      fval = -0.5D1 / 0.2D1 * SQRT(0.14D2) * x
    CASE(5)
      fval = -0.9D1 / 0.4D1 * SQRT(0.2D1) * DBLE(7 * x ** 2 - 1)
    CASE(6)
      fval = -0.7D1 / 0.2D1 * SQRT(0.22D2) * DBLE(x) * DBLE(3 * x ** 2 - 1)
    CASE(7)
      fval = -0.5D1 / 0.8D1 * SQRT(0.26D2) * DBLE(1 + (33 * x ** 2 - 18) * &
              x ** 2)
    CASE(8)
      fval = -0.9D1 / 0.32D2 * SQRT(0.30D2) * DBLE(x) * DBLE(15 + (143 * x ** 2 -&
              110) * x ** 2)
    CASE(9)
      fval = -0.35D2 / 0.64D2 * SQRT(0.34D2) * DBLE(-1 + (33 + (143 * x ** 2 - 143) &
              * x ** 2) * x ** 2)
    CASE(10)
      fval = -0.11D2 / 0.16D2 * SQRT(0.38D2) * DBLE(x) * DBLE(-7 + (91 + (221 * &
              x ** 2 - 273) * x ** 2) * x ** 2)
    CASE(11)
      fval = -0.9D1 / 0.128D3 * SQRT(0.42D2) * DBLE(7 + (-364 + (2730 + (4199 * &
              x ** 2 - 6188) * x ** 2) * x ** 2) * x ** 2)
    CASE(12)
      fval = -0.65D2 / 0.256D3 * SQRT(0.46D2) * DBLE(x) * DBLE(21 + (-420 + (2142 +&
              (2261 * x ** 2 - 3876) * x ** 2) * x ** 2) * x ** 2)
    CASE(13)
      fval = -0.385D3 / 0.512D3 * SQRT(0.2D1) * DBLE(-3 + (225 + (-2550 + (9690 + &
              (7429 * x ** 2 - 14535) * x ** 2) * x ** 2) * x ** 2) * x ** 2)
    CASE(14)
      fval = -0.45D2 / 0.256D3 * DBLE(x) * SQRT(0.6D1) * DBLE(-99 + (2805 + (-21318 + &
              (63954 + (37145 * x ** 2 - 81719) * x ** 2) * x ** 2) * x ** 2) * x ** 2)
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_dVarPhi

  PURE FUNCTION H1Basis_LegendreP(k, x) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: x
    ! Return value
    REAL (KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(k) _ELMER_LINEAR_REF(x) NOTINBRANCH
    
    SELECT CASE(k)
    CASE(0)
      fval = 1
    CASE(1)

      fval = x
    CASE(2)

      fval = -0.1D1 / 0.2D1 + 0.3D1 / 0.2D1 * x ** 2
    CASE(3)

      fval = (-0.3D1 / 0.2D1 + 0.5D1 / 0.2D1 * x ** 2) * x
    CASE(4)

      fval = 0.3D1 / 0.8D1 + (-0.15D2 / 0.4D1 + 0.35D2 / 0.8D1 * x ** 2) * &
              x ** 2
    CASE(5)

      fval = (0.15D2 / 0.8D1 + (-0.35D2 / 0.4D1 + 0.63D2 / 0.8D1 * x ** 2) * &
              x ** 2) * x
    CASE(6)

      fval = -0.5D1 / 0.16D2 + (0.105D3 / 0.16D2 + (-0.315D3 / 0.16D2 + 0.231D3 / &
              0.16D2 * x ** 2) * x ** 2) * x ** 2
    CASE(7)
      fval = (-0.35D2 / 0.16D2 + (0.315D3 / 0.16D2 + (-0.693D3 / 0.16D2 + 0.429D3 / &
              0.16D2 * x ** 2) * x ** 2) * x ** 2) * x
    CASE(8)                              
      fval = 0.35D2 / 0.128D3 + (-0.315D3 / 0.32D2 + (0.3465D4 / 0.64D2 + &
              (-0.3003D4 / 0.32D2 + 0.6435D4 / 0.128D3 * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(9)
      fval = (0.315D3 / 0.128D3 + (-0.1155D4 / 0.32D2 + (0.9009D4 / 0.64D2 + &
              (-0.6435D4 / 0.32D2 + 0.12155D5 / 0.128D3 * x ** 2) * x ** 2) * &
              x ** 2) * x ** 2) * x
    CASE(10)
      fval = -0.63D2 / 0.256D3 + (0.3465D4 / 0.256D3 + (-0.15015D5 / 0.128D3 + &
              (0.45045D5 / 0.128D3 + (-0.109395D6 / 0.256D3 + 0.46189D5 / 0.256D3 *&
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(11)
      fval = (-0.693D3 / 0.256D3 + (0.15015D5 / 0.256D3 + (-0.45045D5 / 0.128D3 +&
              (0.109395D6 / 0.128D3 + (-0.230945D6 / 0.256D3 + 0.88179D5 / 0.256D3 * &
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
    CASE(12)
      fval = 0.231D3 / 0.1024D4 + (-0.9009D4 / 0.512D3 + (0.225225D6 / 0.1024D4 + &
              (-0.255255D6 / 0.256D3 + (0.2078505D7 / 0.1024D4 + (-0.969969D6 / 0.512D3 +&
              0.676039D6 / 0.1024D4 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(13)
      fval = (0.3003D4 / 0.1024D4 + (-0.45045D5 / 0.512D3 + (0.765765D6 / 0.1024D4 + &
              (-0.692835D6 / 0.256D3 + (0.4849845D7 / 0.1024D4 + (-0.2028117D7 / 0.512D3 + &
              0.1300075D7 / 0.1024D4 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
              x ** 2) * x
    CASE(14)
      fval = -0.429D3 / 0.2048D4 + (0.45045D5 / 0.2048D4 + (-0.765765D6 / 0.2048D4 +&
              (0.4849845D7 / 0.2048D4 + (-0.14549535D8 / 0.2048D4 + (0.22309287D8 / &
              0.2048D4 + (-0.16900975D8 / 0.2048D4 + 0.5014575D7 / 0.2048D4 * x ** 2) *&
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_LegendreP

  PURE FUNCTION H1Basis_dLegendreP(k, x) RESULT(fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: x
    ! Return value
    REAL (KIND=dp) :: fval
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(k) _ELMER_LINEAR_REF(x) NOTINBRANCH
    
    SELECT CASE(k)
    CASE(0)
      fval = 0
    CASE(1)
      fval = 1
    CASE(2)
      fval = 3 * x
    CASE(3)
      fval = 0.15D2 / 0.2D1 * x ** 2 - 0.3D1 / 0.2D1
    CASE(4)
      fval = 0.5D1 / 0.2D1 * DBLE(7 * x ** 2 - 3) * DBLE(x)
    CASE(5)
      fval = 0.15D2 / 0.8D1 + (-0.105D3 / 0.4D1 + 0.315D3 / 0.8D1 * x ** 2) *&
              x ** 2
    CASE(6)
      fval = 0.21D2 / 0.8D1 * DBLE(5 + (33 * x ** 2 - 30) * x ** 2) * DBLE(x)
    CASE(7)
      fval = -0.35D2 / 0.16D2 + (0.945D3 / 0.16D2 + (-0.3465D4 / 0.16D2 + 0.3003D4 / &
              0.16D2 * x ** 2) * x ** 2) * x ** 2
    CASE(8)
      fval = 0.9D1 / 0.16D2 * DBLE(-35 + (385 + (715 * x ** 2 - 1001) * x ** 2) * x ** 2) *&
              DBLE(x)
    CASE(9)
      fval = 0.315D3 / 0.128D3 + (-0.3465D4 / 0.32D2 + (0.45045D5 / 0.64D2 + (-0.45045D5 / &
              0.32D2 + 0.109395D6 / 0.128D3 * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(10)
      fval = 0.55D2 / 0.128D3 * DBLE(63 + (-1092 + (4914 + (4199 * x ** 2 - 7956) * x ** 2) * &
              x ** 2) * x ** 2) * DBLE(x)
    CASE(11)
      fval = -0.693D3 / 0.256D3 + (0.45045D5 / 0.256D3 + (-0.225225D6 / 0.128D3 + &
              (0.765765D6 / 0.128D3 + (-0.2078505D7 / 0.256D3 + 0.969969D6 / 0.256D3 *&
              x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(12)
      fval = 0.39D2 / 0.256D3 * DBLE(-231 + (5775 + (-39270 + (106590 + (52003 * x ** 2 - &
              124355) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * DBLE(x)
    CASE(13)
      fval = 0.3003D4 / 0.1024D4 + (-0.135135D6 / 0.512D3 + (0.3828825D7 / 0.1024D4 + &
              (-0.4849845D7 / 0.256D3 + (0.43648605D8 / 0.1024D4 + (-0.22309287D8 / 0.512D3 + &
              0.16900975D8 / 0.1024D4 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE(14)
      fval = 0.105D3 / 0.1024D4 * DBLE(429 + (-14586 + (138567 + (-554268 + (1062347 + (334305 *&
              x ** 2 - 965770) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * DBLE(x)
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_dLegendreP

  ! To be deprecated
  
  ! WARNING: this is not a barycentric triangle
  SUBROUTINE H1Basis_TriangleNodal(nvec, u, v, nbasismax, fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER :: j

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,1) = 1-u(j)-v(j)
      ! Node 2
      fval(j,2) = u(j)
      ! Node 3
      fval(j,3) = v(j)
    END DO
  END SUBROUTINE H1Basis_TriangleNodal

  ! WARNING: this is not a barycentric triangle
  SUBROUTINE H1Basis_dTriangleNodal(nvec, u, v, nbasismax, grad)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER :: j

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,1) = -1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,1) =  1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,1) =  0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,2) = -1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,2) =  0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,2) =  1
    END DO
  END SUBROUTINE H1Basis_dTriangleNodal
  
  ! WARNING: this is not a barycentric tetra
  SUBROUTINE H1Basis_TetraNodal(nvec, u, v, w, nbasismax, fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    INTEGER :: j

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,1) = 1-u(j)-v(j)-w(j)
      ! Node 2
      fval(j,2) = u(j)
      ! Node 3
      fval(j,3) = v(j)
      ! Node 4
      fval(j,4) = w(j)
    END DO
  END SUBROUTINE H1Basis_TetraNodal

  ! WARNING: this is not a barycentric tetra
  SUBROUTINE H1Basis_dTetraNodal(nvec, u, v, w, nbasismax, grad)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    INTEGER :: j

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,1) = -1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,1) =  1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,1) =  0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,4,1) =  0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,2) = -1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,2) =  0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,2) =  1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,4,2) =  0
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,3) = -1
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,3) =  0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,3) =  0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,4,3) =  1
    END DO
  END SUBROUTINE H1Basis_dTetraNodal

  ! WARNING: this is not a barycentric wedge
  SUBROUTINE H1Basis_WedgeNodal(nvec, u, v, w, nbasismax, fval)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax), INTENT(INOUT) :: fval
    REAL(Kind=dp), PARAMETER :: c = 1D0/2
    INTEGER :: j

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,1) = c*(1-u(j)-v(j))*(1-w(j))
      ! Node 2
      fval(j,2) = c*u(j)*(1-w(j))
      ! Node 3
      fval(j,3) = c*v(j)*(1-w(j))
      ! Node 4
      fval(j,4) = c*(1-u(j)-v(j))*(1+w(j))
      ! Node 5
      fval(j,5) = c*u(j)*(1+w(j))
      ! Node 6
      fval(j,6) = c*v(j)*(1+w(j))
    END DO
  END SUBROUTINE H1Basis_WedgeNodal

  ! WARNING: this is not a barycentric wedge
  SUBROUTINE H1Basis_dWedgeNodal(nvec, u, v, w, nbasismax, grad)
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: nvec
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH), INTENT(IN) :: u,v,w
    INTEGER, INTENT(IN) :: nbasismax
    ! Variables
    REAL(Kind=dp), DIMENSION(VECTOR_BLOCK_LENGTH,nbasismax,3), INTENT(INOUT) :: grad
    REAL(Kind=dp), PARAMETER :: c = 1D0/2
    INTEGER :: j

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,1) = -c*(1-w(j))
      grad(j,2,1) =  c*(1-w(j))
      grad(j,3,1) =  0
      grad(j,4,1) = -c*(1+w(j))
      grad(j,5,1) =  c*(1+w(j))
      grad(j,6,1) =  0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,2) = -c*(1-w(j))
      grad(j,2,2) =  0
      grad(j,3,2) =  c*(1-w(j))
      grad(j,4,2) = -c*(1+w(j))
      grad(j,5,2) =  0
      grad(j,6,2) =  c*(1+w(j))
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,3) = -c*(1-u(j)-v(j))
      grad(j,2,3) = -c*u(j)
      grad(j,3,3) = -c*v(j)
      grad(j,4,3) =  c*(1-u(j)-v(j))
      grad(j,5,3) =  c*u(j)
      grad(j,6,3) =  c*v(j)
    END DO
  END SUBROUTINE H1Basis_dWedgeNodal

END MODULE H1Basis
