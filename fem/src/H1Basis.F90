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
! *  Authors: Mikko Byckling
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
  USE Types, ONLY : dp, VECTOR_BLOCK_LENGTH
  USE Messages
  
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
      grad(j,nbasis+2,1) = c
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
      DO p=2,pmax
        !_ELMER_OMP_SIMD 
        DO j=1,nvec
           fval(j,nbasis+p-1) = H1Basis_Phi(p,-u(j))
        END DO
      END DO
    ELSE
      DO p=2,pmax
        !_ELMER_OMP_SIMD 
        DO j=1,nvec
          fval(j,nbasis+p-1) = H1Basis_Phi(p,u(j))
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
      DO p=2,pmax
        ! First coordinate (xi)
        !_ELMER_OMP_SIMD
        DO j=1,nvec
          grad(j,nbasis+p-1,1) = H1Basis_dPhi(p,-u(j))
        END DO
      END DO
    ELSE
      DO p=2,pmax
        ! First coordinate (xi)
        !_ELMER_OMP_SIMD
        DO j=1,nvec
          grad(j,nbasis+p-1,1) = H1Basis_dPhi(p,u(j))
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
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1) = -1d0/2
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,1) = 1d0/2
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,1) = REAL(0,dp)
    END DO
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = -SQRT(3d0)/6
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,2) = -SQRT(3d0)/6
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,2) = SQRT(3d0)/3
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
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb)
        DO k=1,nvec
          La = H1Basis_TriangleL(edgedir(1,i),u(k),v(k))
          Lb = H1Basis_TriangleL(edgedir(2,i),u(k),v(k))

          fval(k, nbasis+j-1) = La*Lb*H1Basis_varPhi(j,Lb-La)
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

      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, vPhi, dVPhi)
        DO k=1,nvec
          La = H1Basis_TriangleL(edgedir(1,i),u(k),v(k))
          Lb = H1Basis_TriangleL(edgedir(2,i),u(k),v(k))
          vPhi = H1Basis_varPhi(j,Lb-La)
          dVPhi = H1Basis_dVarPhi(j,Lb-La)

          grad(k, nbasis+j-1, 1) = dLa(1)*Lb*vPhi+ &
                  La*dLb(1)*vPhi +&
                  La*Lb*dVPhi*(dLb(1)-dLa(1))
          grad(k, nbasis+j-1, 2) = dLa(2)*Lb*vPhi+ &
                  La*dLb(2)*vPhi +&
                  La*Lb*dVPhi*(dLb(2)-dLa(2))
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
        DO j = 0,pmax-i-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc)
          DO k=1,nvec
            La = H1Basis_TriangleL(1,u(k),v(k))
            Lb = H1Basis_TriangleL(2,u(k),v(k))
            Lc = H1Basis_TriangleL(3,u(k),v(k))
          
            fval(k,nbasis+j+1) = La*Lb*Lc*(H1Basis_PowInt((Lb-La),i))*&
                    H1Basis_PowInt((2D0*Lc-1),j)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-3) + 1
        nbasis = nbasis + MAX(pmax-i-2,0)
      END DO
    ELSE
      ! Triangular face basis
      DO i=0,pmax-3
        DO j=0,pmax-i-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc)
          DO k=1,nvec
            La=H1Basis_TriangleL(localnumbers(1),u(k),v(k))
            Lb=H1Basis_TriangleL(localnumbers(2),u(k),v(k))
            Lc=H1Basis_TriangleL(localnumbers(3),u(k),v(k))
            
            fval(k,nbasis+j+1) = La*Lb*Lc*H1Basis_LegendreP(i,Lb-La)* &
                                          H1Basis_LegendreP(j,2*Lc-1)
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
        DO j = 0,pmax-i-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Lb_Lai, Lc_1n)
          DO k=1,nvec
            La = H1Basis_TriangleL(1,u(k),v(k))
            Lb = H1Basis_TriangleL(2,u(k),v(k))
            Lc = H1Basis_TriangleL(3,u(k),v(k))
            
            Lb_Lai = H1Basis_PowInt((Lb-La), i)
            Lc_1n = H1Basis_PowInt((2D0*Lc-1), j)
            
            ! Calculate value of function from general form
            grad(k,nbasis+j+1,1) = -c*Lb*Lc*Lb_Lai*Lc_1n + La*c*Lc*Lb_Lai*Lc_1n + &
                    La*Lb*Lc*i*(H1Basis_PowInt((Lb-La),i-1))*Lc_1n
            grad(k,nbasis+j+1,2) = d*Lb*Lc*Lb_Lai*Lc_1n + La*d*Lc*Lb_Lai*Lc_1n + &
                    La*Lb*e*Lb_Lai*Lc_1n + &
                    La*Lb*Lc*Lb_Lai*j*(H1Basis_PowInt((2D0*Lc-1),j-1))*2D0*e
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
        DO j=0,pmax-i-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Lb_La, Lc_1, a, b)
          DO k=1,nvec
            La=H1Basis_TriangleL(localnumbers(1),u(k),v(k))
            Lb=H1Basis_TriangleL(localnumbers(2),u(k),v(k))
            Lc=H1Basis_TriangleL(localnumbers(3),u(k),v(k))
            Lb_La=Lb-La
            Lc_1=2*Lc-1
            a=H1Basis_LegendreP(i,Lb_La)
            b=H1Basis_LegendreP(j,Lc_1)

            grad(k,nbasis+j+1,1) = dLa(1)*Lb*Lc*a*b + La*dLb(1)*Lc*a*b + &
                    La*Lb*dLc(1)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(i,Lb_La)*(dLb(1)-dLa(1))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(j,Lc_1)*2*dLc(1)
            grad(k,nbasis+j+1,2) = dLa(2)*Lb*Lc*a*b + La*dLb(2)*Lc*a*b + &
                    La*Lb*dLc(2)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(i,Lb_La)*(dLb(2)-dLa(2))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(j,Lc_1)*2*dLc(2)
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
      grad(j,nbasis+2,1) = c*(1-v(j))
      grad(j,nbasis+3,1) = c*(1+v(j))
      grad(j,nbasis+4,1) = -c*(1+v(j))
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = -c*(1-u(j))
      grad(j,nbasis+2,2) = -c*(1+u(j))
      grad(j,nbasis+3,2) = c*(1+u(j))
      grad(j,nbasis+4,2) = c*(1-u(j))
    END DO
    
    nbasis = nbasis + 4
  END SUBROUTINE H1Basis_dQuadNodal

  SUBROUTINE H1Basis_QuadEdgeP(nvec, u, v, pmax, nbasismax, fval, nbasis, edgedir)
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
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb)
        DO k=1,nvec
          La = H1Basis_QuadL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_QuadL(edgedir(2,i), u(k), v(k))

          fval(k, nbasis+j-1) = c*(La+Lb-1)*H1Basis_Phi(j, Lb-La)
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

    REAL(KIND=dp) :: La, Lb, Phi, dPhi, dLa(2), dLb(2)
    REAL(KIND=dp), PARAMETER :: c = 1/2D0
    INTEGER :: i,j,k
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64
    
    ! For each edge
    DO i=1,4
      dLa = H1Basis_dQuadL(edgedir(1,i))
      dLb = H1Basis_dQuadL(edgedir(2,i))

      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Phi, dPhi)
        DO k=1,nvec
          La = H1Basis_QuadL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_QuadL(edgedir(2,i), u(k), v(k))
          
          Phi = H1Basis_Phi(j, Lb-La)
          dPhi = H1Basis_dPhi(j,Lb-La)

          grad(k, nbasis+j-1, 1) = c*((dLa(1)+dLb(1))*Phi + &
                  (La+Lb-1)*dPhi*(dLb(1)-dLa(1)))
          grad(k, nbasis+j-1, 2) = c*((dLa(2)+dLb(2))*Phi + &
                  (La+Lb-1)*dPhi*(dLb(2)-dLa(2)))
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
    REAL(KIND=dp) :: La, Lb, Lc
!DIR$ ASSUME_ALIGNED u:64, v:64, fval:64

    ! Calculate value of function without direction and return
    ! if local numbering not present
    IF (.NOT. PRESENT(localNumbers)) THEN
      DO i=2,(pmax-2)
        DO j=2,(pmax-i)
          !_ELMER_OMP_SIMD
          DO k=1,nvec
            fval(k,nbasis+j-1) = H1Basis_Phi(i,u(k))*H1Basis_Phi(j,v(k))
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-1,0)
      END DO
    ELSE
      DO i=2,(pmax-2)
        DO j=2,(pmax-i)
          !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lc)
          DO k=1,nvec
            ! Directed quad bubbles
            La = H1Basis_QuadL(localNumbers(1),u(k),v(k))
            Lb = H1Basis_QuadL(localNumbers(2),u(k),v(k))
            Lc = H1Basis_QuadL(localNumbers(4),u(k),v(k))

            ! Calculate value of function from the general form
            fval(k,nbasis+j-1) = H1Basis_Phi(i,Lb-La)*H1Basis_Phi(j,Lc-La)
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-1,0)
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

    INTEGER :: i,j,k
    REAL(KIND=dp) :: La, Lb, Lc
    REAL(Kind=dp), DIMENSION(2) :: dLa, dLb, dLc, dLbdLa, dLcdLa
!DIR$ ASSUME_ALIGNED u:64, v:64, grad:64

    IF (.NOT. PRESENT(localNumbers)) THEN
      DO i=2,(pmax-2)
        DO j=2,(pmax-i)
          !_ELMER_OMP_SIMD
          DO k=1,nvec
            ! First coordinate (Xi)
            grad(k,nbasis+j-1,1) = H1Basis_dPhi(i,u(k))*H1Basis_Phi(j,v(k))
            ! Second coordinate (Eta)
            grad(k,nbasis+j-1,2) = H1Basis_Phi(i,u(k))*H1Basis_dPhi(j,v(k))
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
        DO j=2,(pmax-i)
          !_ELMER_OMP_SIMD PRIVATE(La,Lb,Lc)
          DO k=1,nvec
            La = H1Basis_QuadL(localNumbers(1),u(k),v(k))
            Lb = H1Basis_QuadL(localNumbers(2),u(k),v(k))
            Lc = H1Basis_QuadL(localNumbers(4),u(k),v(k))

            grad(k,nbasis+j-1,1) = H1Basis_dPhi(i,Lb-La)*(dLbdLa(1))*H1Basis_Phi(j,Lc-La) + &
                    H1Basis_Phi(i,Lb-La)*H1Basis_dPhi(j,Lc-La)*(dLcdLa(1))
            grad(k,nbasis+j-1,2) = H1Basis_dPhi(i,Lb-La)*(dLbdLa(2))*H1Basis_Phi(j,Lc-La) + &
                    H1Basis_Phi(i,Lb-La)*H1Basis_dPhi(j,Lc-La)*(dLcdLa(2))
          END DO
        END DO
        nbasis = nbasis+MAX((pmax-i)-1,0)
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
      fval = c*(2d0-u-v)
    CASE (2)
      fval = c*(2d0+u-v)
    CASE (3)
      fval = c*(2d0+u+v)
    CASE (4)
      fval = c*(2d0-u+v)
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
            e = 1D0/SQRT(6D0), f = 1D0/SQRT(8D0)
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    !_ELMER_OMP_SIMD
    DO j=1,nvec
      ! Node 1
      fval(j,nbasis+1) = c*(1-u(j)-d*v(j)-e*w(j))
      ! Node 2
      fval(j,nbasis+2) = c*(1+u(j)-d*v(j)-e*w(j))
      ! Node 3
      fval(j,nbasis+3) = SQRT(3d0)/3*(v(j)-f*w(j))
      ! Node 4
      fval(j,nbasis+4) = SQRT(3d0/8d0)*w(j)
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
            e = 1D0/SQRT(6D0), f = 1D0/SQRT(8D0)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) &
    !_ELMER_OMP _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) _ELMER_LINEAR_REF(w) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1)
      fval = c*(1-u-d*v-e*w)
    CASE(2)
      fval = c*(1+u-d*v-e*w)
    CASE(3)
      fval = SQRT(3d0)/3*(v-f*w)
    CASE(4)
      fval = SQRT(3d0/8d0)*w      
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
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,1)=-1d0/2
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,1)=1d0/2
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,1)=REAL(0,dp)
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+4,1)=REAL(0,dp)
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2)=-SQRT(3d0)/6
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,2)=-SQRT(3d0)/6
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,2)=SQRT(3d0)/3
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+4,2)=REAL(0,dp)
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,3)=-SQRT(6d0)/12
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+2,3)=-SQRT(6d0)/12
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+3,3)=-SQRT(6d0)/12
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+4,3)=SQRT(6d0)/4
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
            e = 1D0/SQRT(6D0), f = 1D0/SQRT(8D0)
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) NOTINBRANCH
    
    SELECT CASE(node)
    CASE(1)
      grad = c*[-1D0, -d, -e]
    CASE(2)
      grad = c*[1D0, -d, -e]
    CASE(3)
      grad = d*[0D0, 1D0, -f]
    CASE(4)
      grad = [0D0, 0D0, SQRT(3d0/8d0)]
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
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb)
        DO k=1,nvec
          La = H1Basis_TetraL(edgedir(1,i), u(k), v(k), w(k))
          Lb = H1Basis_TetraL(edgedir(2,i), u(k), v(k), w(k))

          fval(k, nbasis+j-1) = La*Lb*H1Basis_varPhi(j,Lb-La)
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

      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, vPhi, dVPhi)
        DO k=1,nvec
          La = H1Basis_TetraL(edgedir(1,i),u(k),v(k),w(k))
          Lb = H1Basis_TetraL(edgedir(2,i),u(k),v(k),w(k))
          vPhi = H1Basis_varPhi(j,Lb-La)
          dVPhi = H1Basis_dVarPhi(j,Lb-La)

          grad(k, nbasis+j-1, 1) = dLa(1)*Lb*vPhi+ &
                  La*dLb(1)*vPhi +&
                  La*Lb*dVPhi*(dLb(1)-dLa(1))
          grad(k, nbasis+j-1, 2) = dLa(2)*Lb*vPhi+ &
                  La*dLb(2)*vPhi +&
                  La*Lb*dVPhi*(dLb(2)-dLa(2))
          grad(k, nbasis+j-1, 3) = dLa(3)*Lb*vPhi+ &
                  La*dLb(3)*vPhi +&
                  La*Lb*dVPhi*(dLb(3)-dLa(3))
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
        DO k=0,pmax(i)-j-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc)
          DO l=1,nvec
            La=H1Basis_TetraL(facedir(1,i),u(l),v(l),w(l))
            Lb=H1Basis_TetraL(facedir(2,i),u(l),v(l),w(l))
            Lc=H1Basis_TetraL(facedir(3,i),u(l),v(l),w(l))
    
            fval(l,nbasis+k+1) = La*Lb*Lc*H1Basis_LegendreP(j,Lb-La)* &
                                        H1Basis_LegendreP(k,2*Lc-1)
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
        DO k=0,pmax(i)-j-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Lb_La, Lc_1, a, b)
          DO l=1,nvec
            La=H1Basis_TetraL(facedir(1,i),u(l),v(l),w(l))
            Lb=H1Basis_TetraL(facedir(2,i),u(l),v(l),w(l))
            Lc=H1Basis_TetraL(facedir(3,i),u(l),v(l),w(l))
            Lb_La=Lb-La
            Lc_1=2*Lc-1
            a=H1Basis_LegendreP(j,Lb_La)
            b=H1Basis_LegendreP(k,Lc_1)

            grad(l,nbasis+k+1,1) = dLa(1)*Lb*Lc*a*b + La*dLb(1)*Lc*a*b + &
                    La*Lb*dLc(1)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(j,Lb_La)*(dLb(1)-dLa(1))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(k,Lc_1)*2*dLc(1)
            grad(l,nbasis+k+1,2) = dLa(2)*Lb*Lc*a*b + La*dLb(2)*Lc*a*b + &
                    La*Lb*dLc(2)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(j,Lb_La)*(dLb(2)-dLa(2))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(k,Lc_1)*2*dLc(2)
            grad(l,nbasis+k+1,3) = dLa(3)*Lb*Lc*a*b + La*dLb(3)*Lc*a*b + &
                    La*Lb*dLc(3)*a*b + &
                    La*Lb*Lc*H1Basis_dLegendreP(j,Lb_La)*(dLb(3)-dLa(3))*b + &
                    La*Lb*Lc*a*H1Basis_dLegendreP(k,Lc_1)*2*dLc(3)
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
        DO k=0,pmax-i-j-4
          !_ELMER_OMP_SIMD PRIVATE(L1,L2,L3,L4,L2_L1,L3_1,L4_1)
          DO l=1,nvec
            L1 = H1Basis_TetraL(1,u(l),v(l),w(l))
            L2 = H1Basis_TetraL(2,u(l),v(l),w(l))
            L3 = H1Basis_TetraL(3,u(l),v(l),w(l))
            L4 = H1Basis_TetraL(4,u(l),v(l),w(l))
            L2_L1 = L2-L1
            L3_1 = 2*L3-1
            L4_1 = 2*L4-1

            fval(l,nbasis+k+1) = L1*L2*L3*L4*&
                    H1Basis_LegendreP(i,L2_L1)*&
                    H1Basis_LegendreP(j,L3_1)*&
                    H1Basis_LegendreP(k,L4_1)
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
        DO k=0,pmax-i-j-4
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
            c = H1Basis_LegendreP(k,L4_1)

            ! Gradients of tetrahedral bubble basis functions 
            grad(l,nbasis+k+1,1) = -1d0/2*L2*L3*L4*a*b*c + &
                    1d0/2*L1*L3*L4*a*b*c + &
                    L1*L2*L3*L4*H1Basis_dLegendreP(i,L2_L1)*b*c
            grad(l,nbasis+k+1,2) = -SQRT(3d0)/6*L2*L3*L4*a*b*c - &
                    SQRT(3d0)/6*L1*L3*L4*a*b*c + &
                    SQRT(3d0)/3*L1*L2*L4*a*b*c &
                    + 2*SQRT(3d0)/3*L1*L2*L3*L4*a*H1Basis_dLegendreP(j,L3_1)*c
            grad(l,nbasis+k+1,3) = -SQRT(6d0)/12*L2*L3*L4*a*b*c - &
                    SQRT(6d0)/12*L1*L3*L4*a*b*c - &
                    SQRT(6d0)/12*L1*L2*L4*a*b*c &
                    + SQRT(6d0)/4*L1*L2*L3*a*b*c - &
                    SQRT(6d0)/6*L1*L2*L3*L4*a*H1Basis_dLegendreP(j,L3_1)*c &
                    + SQRT(6d0)/2*L1*L2*L3*L4*a*b*H1Basis_dLegendreP(k,L4_1)
          END DO
        END DO
        ! nbasis = nbasis + (pmax-i-j-4) + 1
        nbasis = nbasis + MAX(pmax-i-j-3,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dTetraBubbleP

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
      fval(j,nbasis+1) = c*(1d0-u(j)-d*v(j))*(1d0-w(j))
      ! Node 2
      fval(j,nbasis+2) = c*(1d0+u(j)-d*v(j))*(1d0-w(j))
      ! Node 3
      fval(j,nbasis+3) = e*v(j)*(1-w(j))
      ! Node 4
      fval(j,nbasis+4) = c*(1d0-u(j)-d*v(j))*(1d0+w(j))
      ! Node 5
      fval(j,nbasis+5) = c*(1d0+u(j)-d*v(j))*(1d0+w(j))
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
      grad(j,nbasis+2,1) = c*(1-w(j))  
      grad(j,nbasis+3,1) = REAL(0, dp)
      grad(j,nbasis+4,1) = -c*(1+w(j))
      grad(j,nbasis+5,1) = c*(1+w(j))
      grad(j,nbasis+6,1) = REAL(0, dp)
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2)= -f*(1-w(j))
      grad(j,nbasis+2,2) = -f*(1-w(j))
      grad(j,nbasis+3,2) = e*(1-w(j))
      grad(j,nbasis+4,2) = -f*(1+w(j))
      grad(j,nbasis+5,2) = -f*(1+w(j))
      grad(j,nbasis+6,2) = e*(1+w(j))
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,3) = -c*(1d0-u(j)-d*v(j))
      grad(j,nbasis+2,3) = -c*(1d0+u(j)-d*v(j))
      grad(j,nbasis+3,3) = -e*v(j)
      grad(j,nbasis+4,3) = c*(1d0-u(j)-d*v(j))
      grad(j,nbasis+5,3) = c*(1d0+u(j)-d*v(j))
      grad(j,nbasis+6,3) = e*v(j)
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
    !_ELMER_OMP_DECLARE_SIMD UNIFORM(node) &
    !_ELMER_OMP _ELMER_LINEAR_REF(u) _ELMER_LINEAR_REF(v) NOTINBRANCH
    
    SELECT CASE(node)
    CASE (1,4)
      fval = 1d0/2*(1d0-u-v/SQRT(3d0))
    CASE (2,5)
      fval = 1d0/2*(1d0+u-v/SQRT(3d0))
    CASE (3,6)
      fval = SQRT(3d0)/3*v
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
      grad = [REAL(0,dp), REAL(d,dp), 0D0]
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
      fval = -w/2d0
    CASE (4,5,6)
      fval = w/2d0
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
      grad = [0D0, 0D0, c]
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
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each triangle face edge
    DO i=1,6
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_WedgeL(edgedir(2,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          fval(k, nbasis+j-1) = c*La*Lb*H1Basis_varPhi(j,Lb-La)*(1+Na+Nb)
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
      nbasis = nbasis + pmax(i) - 1
    END DO
      
    ! For each square face edge
    DO i=7,9
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Na, Nb)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          fval(k, nbasis+j-1) = La*H1Basis_Phi(j, Nb-Na)
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
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Na, Nb, vPhi, dVPhi, NaNb)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Lb = H1Basis_WedgeL(edgedir(2,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          vPhi=H1Basis_varPhi(j,Lb-La)
          dVPhi=H1Basis_dVarPhi(j,Lb-La)
          NaNb=1+Na+Nb
          ! fval(k, nbasis+j-1) = c*La*Lb*H1Basis_varPhi(j,Lb-La)*(1+Na+Nb)
          grad(k, nbasis+j-1,1) = c*dLa(1)*Lb*vPhi*NaNb + &
                c*La*dLb(1)*vPhi*NaNb + c*La*Lb*dVPhi*(dLb(1)-dLa(1))*NaNb
                ! + c*La*Lb*vPhi*(dNb(1)+dNa(1)) ! Last term zero
          grad(k, nbasis+j-1,2) = c*dLa(2)*Lb*vPhi*NaNb + &
                c*La*dLb(2)*vPhi*NaNb + c*La*Lb*dVPhi*(dLb(2)-dLa(2))*NaNb
                ! c*La*Lb*vPhi*(dNb(2)+dNa(2)) ! Last term zero
          grad(k, nbasis+j-1,3) = c*La*Lb*vPhi*(dNb(3)+dNa(3))
                ! c*dLa(3)*Lb*vPhi*NaNb + ! First term zero
                ! c*La*dLb(3)*vPhi*NaNb + ! Second term zero
                ! c*La*Lb*dVPhi*(dLb(3)-dLa(3))*NaNb ! Third term zero
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
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Na, Nb, Phi, dPhi)
        DO k=1,nvec
          La = H1Basis_WedgeL(edgedir(1,i), u(k), v(k))
          Na = H1Basis_WedgeH(edgedir(1,i), w(k))
          Nb = H1Basis_WedgeH(edgedir(2,i), w(k))
          Phi = H1Basis_Phi(j, Nb-Na)
          dPhi = H1Basis_dPhi(j,Nb-Na)
          
          grad(k, nbasis+j-1,1) = dLa(1)*Phi+La*dPhi*(dNb(1)-dNa(1))
          grad(k, nbasis+j-1,2) = dLa(2)*Phi+La*dPhi*(dNb(2)-dNa(2))
          grad(k, nbasis+j-1,3) = dLa(3)*Phi+La*dPhi*(dNb(3)-dNa(3)) 
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

    REAL(KIND=dp) :: La, Lb, Lc, Na, Nc
    REAL(KIND=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: i,j,k,l
    LOGICAL :: nonpermuted
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! Triangle faces
    DO i=1,2      
      DO j=0,pmax(i)-3
        DO k=0,pmax(i)-j-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Na)
          DO l=1,nvec
            La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
            Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
            Lc = H1Basis_WedgeL(facedir(3,i), u(l), v(l))
            Na = H1Basis_WedgeH(facedir(1,i), w(l))

            fval(l,nbasis+k+1) = c*(1+2*Na)*La*Lb*Lc* &
                                            H1Basis_LegendreP(j, Lb-La)* &
                                            H1Basis_LegendreP(k, 2*Lc-1)
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
        DO k=2,pmax(i)-j          
          IF (nonpermuted) THEN
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc)
            DO l=1,nvec
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(4,i), w(l))
              fval(l,nbasis+k-1) = La*Lb*H1Basis_varPhi(j, Lb-La)* &
                                         H1Basis_Phi(k, Nc-Na)
            END DO
          ELSE
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc)
            DO l=1,nvec
              ! Numbering permuted, switch indices j,k and nodes 2 and 4
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(4,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(2,i), w(l))
              fval(l,nbasis+k-1) = La*Lb*H1Basis_varPhi(k, Lb-La)* &
                                         H1Basis_Phi(j, Nc-Na)
            END DO
          END IF
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
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
        DO k=0,pmax(i)-j-3
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Na, LegPLbLa, LegP2Lc1, cNa)
          DO l=1,nvec
            La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
            Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
            Lc = H1Basis_WedgeL(facedir(3,i), u(l), v(l))
            Na = H1Basis_WedgeH(facedir(1,i), w(l))

            LegPLbLa = H1Basis_LegendreP(j, Lb-La)
            LegP2Lc1 = H1Basis_LegendreP(k, 2*Lc-1)
            cNa = c*(1+2*Na)

            grad(l,nbasis+k+1,1) = cNa*dLa(1)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*dLb(1)*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*dLc(1)*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(1)-dLa(1))*LegP2Lc1 + &
                    cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k, 2*Lc-1)*(2*dLc(1))
            ! + 2*c*dNa(1)*La*Lb*Lc*LegPLbLa*LegP2Lc1 ! Term always zero
            grad(l,nbasis+k+1,2) = cNa*dLa(2)*Lb*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*dLb(2)*Lc*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*dLc(2)*LegPLbLa*LegP2Lc1 + &
                    cNa*La*Lb*Lc*H1Basis_dLegendreP(j, Lb-La)*(dLb(2)-dLa(2))*LegP2Lc1 + &
                    cNa*La*Lb*Lc*LegPLbLa*H1Basis_dLegendreP(k, 2*Lc-1)*(2*dLc(2))
            ! + 2*c*dNa(2)*La*Lb*Lc*LegPLbLa*LegP2Lc1 ! Term always zero
            grad(l,nbasis+k+1,3) = 2*c*dNa(3)*La*Lb*Lc*LegPLbLa*LegP2Lc1
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
        DO k=2,pmax(i)-j          
          IF (nonpermuted) THEN
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc,vPhi,Phi)
            DO l=1,nvec
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(2,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(4,i), w(l))
              
              vPhi = H1Basis_varPhi(j, Lb-La)
              Phi = H1Basis_Phi(k, Nc-Na)

              ! fval(l,nbasis+k-1) = La*Lb*H1Basis_varPhi(j, Lb-La)* &
              !                            H1Basis_Phi(k, Nc-Na)
              grad(l,nbasis+k-1,1) = dLa(1)*Lb*vPhi*Phi+ &
                      La*dLb(1)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(j, Lb-La)*(dLb(1)-dLa(1))*Phi
              ! La*Lb*vPhi*H1Basis_dPhi(k, Nc-Na)*(dNc(1)-dNa(1)) ! Term always zero
              grad(l,nbasis+k-1,2) = dLa(2)*Lb*vPhi*Phi+ &
                      La*dLb(2)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(j, Lb-La)*(dLb(2)-dLa(2))*Phi
              ! La*Lb*vPhi*H1Basis_dPhi(k, Nc-Na)*(dNc(2)-dNa(2)) ! Term always zero
              grad(l,nbasis+k-1,3) = La*Lb*vPhi*H1Basis_dPhi(k, Nc-Na)*(dNc(3)-dNa(3))
              ! Rest of the terms always zero
              ! dLa(3)*Lb*vPhi*Phi+ &
              ! La*dLb(3)*vPhi*Phi + &
              ! La*Lb*H1Basis_dVarPhi(j, Lb-La)*(dLb(3)-dLa(3))*Phi
            END DO
          ELSE
            !_ELMER_OMP_SIMD PRIVATE(La,Lb,Na,Nc)
            DO l=1,nvec
              ! Numbering permuted, switch indices j,k and nodes 2 and 4
              La = H1Basis_WedgeL(facedir(1,i), u(l), v(l))
              Lb = H1Basis_WedgeL(facedir(4,i), u(l), v(l))
              Na = H1Basis_WedgeH(facedir(1,i), w(l))
              Nc = H1Basis_WedgeH(facedir(2,i), w(l))
              
              vPhi = H1Basis_varPhi(k, Lb-La)
              Phi = H1Basis_Phi(j, Nc-Na)

              ! fval(l,nbasis+k-1) = La*Lb*H1Basis_varPhi(k, Lb-La)* &
              !                            H1Basis_Phi(j, Nc-Na)
              grad(l,nbasis+k-1,1) = dLa(1)*Lb*vPhi*Phi+ &
                      La*dLb(1)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(k, Lb-La)*(dLb(1)-dLa(1))*Phi
              ! La*Lb*vPhi*H1Basis_dPhi(j, Nc-Na)*(dNc(1)-dNa(1)) ! Term always zero
              grad(l,nbasis+k-1,2) = dLa(2)*Lb*vPhi*Phi+ &
                      La*dLb(2)*vPhi*Phi + &
                      La*Lb*H1Basis_dVarPhi(k, Lb-La)*(dLb(2)-dLa(2))*Phi
              ! La*Lb*vPhi*H1Basis_dPhi(j, Nc-Na)*(dNc(2)-dNa(2)) ! Term always zero
              grad(l,nbasis+k-1,3) = La*Lb*vPhi*H1Basis_dPhi(j, Nc-Na)*(dNc(3)-dNa(3))
              ! Rest of the terms always zero
              ! dLa(3)*Lb*vPhi*Phi+ &
              ! La*dLb(3)*vPhi*Phi + &
              ! La*Lb*H1Basis_dVarPhi(k, Lb-La)*(dLb(3)-dLa(3))*Phi
            END DO
          END IF
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
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
    REAL(KIND=dp) :: L1, L2, L3, L2_L1, L3_1
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    DO i=0,pmax-5
      DO j=0,pmax-5-i
        DO k=2,pmax-3-i-j
          !_ELMER_OMP_SIMD PRIVATE(L1, L2, L3, L2_L1, L3_1)
          DO l=1,nvec
            L1 = H1Basis_WedgeL(1,u(l),v(l))
            L2 = H1Basis_WedgeL(2,u(l),v(l))
            L3 = H1Basis_WedgeL(3,u(l),v(l))
            L2_L1 = L2-L1
            L3_1 = 2d0*L3-1

            ! Get value of bubble function
            fval(l,nbasis+k-1) = L1*L2*L3*H1Basis_LegendreP(i,L2_L1)*&
                    H1Basis_LegendreP(j,L3_1)*H1Basis_Phi(k,w(l))
          END DO
        END DO
        ! nbasis = nbasis + (pmax-3-i-j)-2+1
        nbasis = nbasis + MAX(pmax-4-i-j,0)
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
    REAL(KIND=dp) :: L1,L2,L3,Legi,Legj,phiW,L2_L1,L3_1
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, grad:64      

    DO i=0,pmax-5
      DO j=0,pmax-5-i
        DO k=2,pmax-3-i-j
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
            phiW = H1Basis_Phi(k,w(l))
            
            grad(l,nbasis+k-1,1) = ((-1d0/2)*L2*L3*Legi*Legj + &
                    L1*(1d0/2)*L3*Legi*Legj +&
                    L1*L2*L3*H1Basis_dLegendreP(i,L2_L1)*Legj)*phiW
            grad(l,nbasis+k-1,2) = ((-SQRT(3d0)/6)*L2*L3*Legi*Legj + &
                    L1*(-SQRT(3d0)/6)*L3*Legi*Legj +&
                    L1*L2*(SQRT(3D0)/3)*Legi*Legj + &
                    L1*L2*L3*Legi*H1Basis_dLegendreP(j,L3_1)*(2d0*SQRT(3D0)/3))*phiW 
            grad(l,nbasis+k-1,3) = L1*L2*L3*Legi*Legj*H1Basis_dPhi(k,w(l))
          END DO
        END DO
        ! nbasis = nbasis + (pmax-3-i-j)-2+1
        nbasis = nbasis + MAX(pmax-4-i-j,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dWedgeBubbleP

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
      grad(j,nbasis+2,1) = c*(1-v(j))*(1-w(j))
      grad(j,nbasis+3,1) = c*(1+v(j))*(1-w(j))
      grad(j,nbasis+4,1) = -c*(1+v(j))*(1-w(j))
      grad(j,nbasis+5,1) = -c*(1-v(j))*(1+w(j))
      grad(j,nbasis+6,1) = c*(1-v(j))*(1+w(j))
      grad(j,nbasis+7,1) = c*(1+v(j))*(1+w(j))
      grad(j,nbasis+8,1) = -c*(1+v(j))*(1+w(j))
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,2) = -c*(1-u(j))*(1-w(j))
      grad(j,nbasis+2,2) = -c*(1+u(j))*(1-w(j))
      grad(j,nbasis+3,2) = c*(1+u(j))*(1-w(j))
      grad(j,nbasis+4,2) = c*(1-u(j))*(1-w(j))
      grad(j,nbasis+5,2) = -c*(1-u(j))*(1+w(j))
      grad(j,nbasis+6,2) = -c*(1+u(j))*(1+w(j))
      grad(j,nbasis+7,2) = c*(1+u(j))*(1+w(j))
      grad(j,nbasis+8,2) = c*(1-u(j))*(1+w(j))
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,nbasis+1,3) = -c*(1-u(j))*(1-v(j))
      grad(j,nbasis+2,3) = -c*(1+u(j))*(1-v(j))
      grad(j,nbasis+3,3) = -c*(1+u(j))*(1+v(j))
      grad(j,nbasis+4,3) = -c*(1-u(j))*(1+v(j))
      grad(j,nbasis+5,3) = c*(1-u(j))*(1-v(j))
      grad(j,nbasis+6,3) = c*(1+u(j))*(1-v(j))
      grad(j,nbasis+7,3) = c*(1+u(j))*(1+v(j))
      grad(j,nbasis+8,3) = c*(1-u(j))*(1+v(j))
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
      fval = c*(3D0-u-v-w)
    CASE (2)
      fval = c*(3D0+u-v-w)
    CASE (3)
      fval = c*(3D0+u+v-w)
    CASE (4)
      fval = c*(3D0-u+v-w)
    CASE (5)
      fval = c*(3D0-u-v+w)
    CASE (6)
      fval = c*(3D0+u-v+w)
    CASE (7)
      fval = c*(3D0+u+v+w)
    CASE (8)
      fval = c*(3D0-u+v+w)
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
  
  SUBROUTINE H1Basis_BrickEdgeP(nvec, u, v, w, pmax, nbasismax, fval, nbasis, edgedir)
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
      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Aa, Ba)
        DO k=1,nvec
          La=H1Basis_BrickL(edgedir(1,i), u(k), v(k), w(k))
          Lb=H1Basis_BrickL(edgedir(2,i), u(k), v(k), w(k))
          CALL H1Basis_BrickEdgeL(i, u(k), v(k), w(k), Aa, Ba)
          
          fval(k, nbasis+j-1) = c*H1Basis_Phi(j, Lb-La)*Aa*Ba
        END DO
      END DO
      ! nbasis = nbasis + (pmax(i)-2) + 1
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

      DO j=2,pmax(i)
        !_ELMER_OMP_SIMD PRIVATE(La, Lb, Aa, Ba, Phi, dPhi)
        DO k=1,nvec
          La=H1Basis_BrickL(edgedir(1,i), u(k), v(k), w(k))
          Lb=H1Basis_BrickL(edgedir(2,i), u(k), v(k), w(k))
          CALL H1Basis_BrickEdgeL(i, u(k), v(k), w(k), Aa, Ba)

          Phi = H1Basis_Phi(j, Lb-La)
          dPhi = H1Basis_dPhi(j, Lb-La)

          grad(k,nbasis+j-1,1) = c*dPhi*(dLb(1)-dLa(1))*Aa*Ba +&
                  c*Phi*dAa(1)*Ba +&
                  c*Phi*Aa*dBa(1)
          grad(k,nbasis+j-1,2) = c*dPhi*(dLb(2)-dLa(2))*Aa*Ba +&
                  c*Phi*dAa(2)*Ba +&
                  c*Phi*Aa*dBa(2)
          grad(k,nbasis+j-1,3) = c*dPhi*(dLb(3)-dLa(3))*Aa*Ba +&
                  c*Phi*dAa(3)*Ba +&
                  c*Phi*Aa*dBa(3)
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

    REAL(KIND=dp) :: La, Lb, Lc, Ld
    REAL(KIND=dp), PARAMETER :: a=1D0/4D0
    INTEGER :: i,j,k,l
!DIR$ ASSUME_ALIGNED u:64, v:64, w:64, fval:64

    ! For each face
    DO i=1,6
      DO j=2,pmax(i)
        DO k=2,pmax(i)-j
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Ld)
          DO l=1,nvec
            La=H1Basis_BrickL(facedir(1,i), u(l), v(l), w(l))
            Lb=H1Basis_BrickL(facedir(2,i), u(l), v(l), w(l))
            Lc=H1Basis_BrickL(facedir(3,i), u(l), v(l), w(l))
            Ld=H1Basis_BrickL(facedir(4,i), u(l), v(l), w(l))
            
            fval(l, nbasis+k-1) = (a*(La+Lb+Lc+Ld)-1)*H1Basis_Phi(j, Lb-La) &
                                                     *H1Basis_Phi(k, Ld-La)
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
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
        DO k=2,pmax(i)-j
          !_ELMER_OMP_SIMD PRIVATE(La, Lb, Lc, Ld, PhiLbLa, PhiLdLa, LaLbLcLd)
          DO l=1,nvec
            La=H1Basis_BrickL(facedir(1,i), u(l), v(l), w(l))
            Lb=H1Basis_BrickL(facedir(2,i), u(l), v(l), w(l))
            Lc=H1Basis_BrickL(facedir(3,i), u(l), v(l), w(l))
            Ld=H1Basis_BrickL(facedir(4,i), u(l), v(l), w(l))
            
            PhiLbLa=H1Basis_Phi(j, Lb-La)
            PhiLdLa=H1Basis_Phi(k, Ld-La)
            LaLbLcLd=a*(La+Lb+Lc+Ld)-1

            grad(l, nbasis+k-1, 1)=a*(dLa(1)+dLb(1)+dLc(1)+dLd(1))*PhiLbLa*PhiLdLa + &
                  LaLbLcLd*H1Basis_dPhi(j,Lb-La)*(dLb(1)-dLa(1))*PhiLdLa + &
                  LaLbLcLd*PhiLbLa*H1Basis_dPhi(k, Ld-La)*(dLd(1)-dLa(1))
            grad(l, nbasis+k-1, 2)=a*(dLa(2)+dLb(2)+dLc(2)+dLd(2))*PhiLbLa*PhiLdLa + &
                  LaLbLcLd*H1Basis_dPhi(j,Lb-La)*(dLb(2)-dLa(2))*PhiLdLa + &
                  LaLbLcLd*PhiLbLa*H1Basis_dPhi(k, Ld-La)*(dLd(2)-dLa(2))
            grad(l, nbasis+k-1, 3)=a*(dLa(3)+dLb(3)+dLc(3)+dLd(3))*PhiLbLa*PhiLdLa + &
                  LaLbLcLd*H1Basis_dPhi(j,Lb-La)*(dLb(3)-dLa(3))*PhiLdLa + &
                  LaLbLcLd*PhiLbLa*H1Basis_dPhi(k, Ld-La)*(dLd(3)-dLa(3))
          END DO
        END DO
        ! nbasis = nbasis + (pmax(i)-j-2) + 1
        nbasis = nbasis + MAX(pmax(i)-j-1,0)
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
    DO i=2,pmax-4
      DO j=2,pmax-i-2
        DO k=2,pmax-i-j
          !_ELMER_OMP_SIMD
          DO l=1,nvec
            fval(l,nbasis+k-1) = H1Basis_Phi(i,u(l))*H1Basis_Phi(j,v(l))*H1Basis_Phi(k,w(l))
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-j-1,0)
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
    DO i=2,pmax-4
      DO j=2,pmax-i-2
        DO k=2,pmax-i-j
          !_ELMER_OMP_SIMD PRIVATE(phiU, phiV, phiW)
          DO l=1,nvec
            phiU = H1Basis_Phi(i,u(l))
            phiV = H1Basis_Phi(j,v(l))
            phiW = H1Basis_Phi(k,w(l))
            grad(l,nbasis+k-1,1) = H1Basis_dPhi(i,u(l))*phiV*phiW
            grad(l,nbasis+k-1,2) = phiU*H1Basis_dPhi(j,v(l))*phiW
            grad(l,nbasis+k-1,3) = phiU*phiV*H1Basis_dPhi(k,w(l))
          END DO
        END DO
        nbasis = nbasis+MAX(pmax-i-j-1,0)
      END DO
    END DO
  END SUBROUTINE H1Basis_dBrickBubbleP

  FUNCTION H1Basis_Phi(k, x) RESULT(fval)
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
!!!      CASE(15)
!!!        fval = (-0.429D3 / 0.4096D4 * sqrt(0.58D2) + (0.15015D5 / 0.4096D4 * sqrt(0.58D2) +&
!!!                (-0.153153D6 / 0.4096D4 * sqrt(0.58D2) + (0.692835D6 / 0.4096D4 * sqrt(0.58D2) +&
!!!                (-0.1616615D7 / 0.4096D4 * sqrt(0.58D2) + (0.2028117D7 / 0.4096D4 * sqrt(0.58D2) +&
!!!                (-0.1300075D7 / 0.4096D4 * sqrt(0.58D2) + 0.334305D6 / 0.4096D4 * sqrt(0.58D2) &
!!!                * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
!!!      CASE(16)
!!!        fval = 0.429D3 / 0.65536D5 * sqrt(0.62D2) + (-0.6435D4 / 0.8192D4 * sqrt(0.62D2) +&
!!!                (0.255255D6 / 0.16384D5 * sqrt(0.62D2) + (-0.969969D6 / 0.8192D4 * sqrt(0.62D2) +&
!!!                (0.14549535D8 / 0.32768D5 * sqrt(0.62D2) + (-0.7436429D7 / 0.8192D4 * &
!!!                sqrt(0.62D2) + (0.16900975D8 / 0.16384D5 * sqrt(0.62D2) + (-0.5014575D7 / 0.8192D4 &
!!!                * sqrt(0.62D2) + 0.9694845D7 / 0.65536D5 * sqrt(0.62D2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
!!!      CASE(17)
!!!        fval = (0.6435D4 / 0.65536D5 * sqrt(0.66D2) + (-0.36465D5 / 0.8192D4 * sqrt(0.66D2) +&
!!!                (0.969969D6 / 0.16384D5 * sqrt(0.66D2) + (-0.2909907D7 / 0.8192D4 * sqrt(0.66D2) +&
!!!                (0.37182145D8 / 0.32768D5 * sqrt(0.66D2) + (-0.16900975D8 / 0.8192D4 * &
!!!                sqrt(0.66D2)  + (0.35102025D8 / 0.16384D5 * sqrt(0.66D2) + &
!!!                (-0.9694845D7 / 0.8192D4 * sqrt(0.66D2) + 0.17678835D8 / 0.65536D5 * sqrt(0.66D2) &
!!!                * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
!!!      CASE(18)
!!!        fval = -0.715D3 / 0.131072D6 * sqrt(0.70D2) + (0.109395D6 / 0.131072D6 * sqrt(0.70D2) +&
!!!                (-0.692835D6 / 0.32768D5 * sqrt(0.70D2) + (0.6789783D7 / 0.32768D5 * &
!!!                sqrt(0.70D2) + (-0.66927861D8 / 0.65536D5 * sqrt(0.70D2) + (0.185910725D9 / &
!!!                0.65536D5 * sqrt(0.70D2) + (-0.152108775D9 / 0.32768D5 * sqrt(0.70D2) + &
!!!                (0.145422675D9 / 0.32768D5 * sqrt(0.70D2) + (-0.300540195D9 / 0.131072D6 * &
!!!                sqrt(0.70D2) + 0.64822395D8 / 0.131072D6 * sqrt(0.70D2) * x ** 2) * x ** 2) &
!!!                * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
!!!      CASE(19)
!!!        fval = (-0.12155D5 / 0.131072D6 * sqrt(0.74D2) + (0.692835D6 / 0.131072D6 * &
!!!                sqrt(0.74D2) + (-0.2909907D7 / 0.32768D5 * sqrt(0.74D2) + (0.22309287D8 /&
!!!                0.32768D5 * sqrt(0.74D2) + (-0.185910725D9 / 0.65536D5 * sqrt(0.74D2) + &
!!!                (0.456326325D9 / 0.65536D5 * sqrt(0.74D2) + (-0.339319575D9 / 0.32768D5 * &
!!!                sqrt(0.74D2) + (0.300540195D9 / 0.32768D5 * sqrt(0.74D2) + (-0.583401555D9 / &
!!!                0.131072D6 * sqrt(0.74D2) + 0.119409675D9 / 0.131072D6 * sqrt(0.74D2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x
!!!      CASE(20)
!!!        fval = 0.2431D4 / 0.524288D6 * sqrt(0.78D2) + (-0.230945D6 / 0.262144D6 * sqrt(0.78D2) + &
!!!                (0.14549535D8 / 0.524288D6 * sqrt(0.78D2) + (-0.22309287D8 / 0.65536D5 * &
!!!                sqrt(0.78D2) + (0.557732175D9 / 0.262144D6 * sqrt(0.78D2) + (-0.1003917915D10 / &
!!!                0.131072D6 * sqrt(0.78D2) + (0.4411154475D10 / 0.262144D6 * sqrt(0.78D2) + &
!!!                (-0.1502700975D10 / 0.65536D5 * sqrt(0.78D2) + (0.9917826435D10 / 0.524288D6 *&
!!!                sqrt(0.78D2) + (-0.2268783825D10 / 0.262144D6 * sqrt(0.78D2) + 0.883631595D9 / &
!!!                0.524288D6 * sqrt(0.78D2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) &
!!!                * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_Phi

  FUNCTION H1Basis_dPhi(k, x) RESULT(fval)
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
!!!      CASE(15)
!!!        fval = dble(-429 + (45045 + (-765765 + (4849845 + (-14549535 + &
!!!                (22309287 + (5014575 * x ** 2 - 16900975) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * sqrt(0.58D2) / 0.4096D4
!!!      CASE(16)
!!!        fval = dble(-6435 + (255255 + (-2909907 + (14549535 + (-37182145 + &
!!!                (50702925 + (9694845 * x ** 2 - 35102025) * x ** 2) * x ** 2) &
!!!                * x ** 2) * x ** 2) * x ** 2) * x ** 2) * dble(x) * sqrt(0.62D2) / 0.4096D4
!!!      CASE(17)
!!!        fval = dble(6435 + (-875160 + (19399380 + (-162954792 + (669278610 + &
!!!                (-1487285800 + (1825305300 + (300540195 * x ** 2 - 1163381400) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
!!!                sqrt(0.66D2) / 0.65536D5
!!!      CASE(18)
!!!        fval = dble(109395 + (-5542680 + (81477396 + (-535422888 + (1859107250 + &
!!!                (-3650610600 + (4071834900 + (583401555 * x ** 2 - 2404321560) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
!!!                dble(x) * sqrt(0.70D2) / 0.65536D5
!!!      CASE(19)
!!!        fval = dble(-12155 + (2078505 + (-58198140 + (624660036 + (-3346393050 + &
!!!                (10039179150 + (-17644617900 + (18032411700 + (2268783825 * x ** 2 -&
!!!                9917826435) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * sqrt(0.74D2) / 0.131072D6
!!!      CASE(20)
!!!        fval = dble(-230945 + (14549535 + (-267711444 + (2230928700 + (-10039179150 +&
!!!                (26466926850 + (-42075627300 + (39671305740 + (4418157975 * x ** 2 -&
!!!                20419054425) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * dble(x) * sqrt(0.78D2) / 0.131072D6
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_dPhi

  FUNCTION H1Basis_VarPhi(k, x) RESULT(fval)
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
!!!      CASE(15)
!!!        fval = -dble(429 + (-14586 + (138567 + (-554268 + (1062347 + (334305 *&
!!!                x ** 2 - 965770) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * dble(x) * sqrt(0.58D2) / 0.1024D4
!!!      CASE(16)
!!!        fval = -dble(-429 + (51051 + (-969969 + (6789783 + (-22309287 + &
!!!                (37182145 + (9694845 * x ** 2 - 30421755) * x ** 2) * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * sqrt(0.62D2) / 0.16384D5
!!!      CASE(17)
!!!        fval = -dble(-6435 + (285285 + (-3594591 + (19684665 + (-54679625 + &
!!!                (80528175 + (17678835 * x ** 2 - 59879925) * x ** 2) * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * dble(x) * sqrt(0.66D2) / 0.16384D5
!!!      CASE(18)
!!!        fval = -dble(715 + (-108680 + (2662660 + (-24496472 + (109359250 + &
!!!                (-262462200 + (345972900 + (64822395 * x ** 2 - 235717800) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * sqrt(0.70D2) / 0.32768D5
!!!      CASE(19)
!!!        fval = -dble(12155 + (-680680 + (10958948 + (-78278200 + (293543250 + &
!!!                (-619109400 + (738168900 + (119409675 * x ** 2 - 463991880) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * dble(x) * sqrt(0.74D2) / 0.32768D5
!!!      CASE(20)
!!!        fval = -dble(-2431 + (459459 + (-14090076 + (164384220 + (-951080130 + &
!!!                (3064591530 + (-5757717420 + (6263890380 + (883631595 * x ** 2 - &
!!!                3653936055) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * sqrt(0.78D2) / 0.131072D6
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_VarPhi

  FUNCTION H1Basis_dVarPhi(k, x) RESULT(fval)
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
!!!      CASE(15)
!!!        fval = -0.13D2 / 0.1024D4 * sqrt(0.58D2) * dble(33 + (-3366 + (53295 + &
!!!                (-298452 + (735471 + (334305 * x ** 2 - 817190) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2)
!!!      CASE(16)
!!!        fval = -0.119D3 / 0.8192D4 * sqrt(0.62D2) * dble(x) * dble(429 + (-16302 + &
!!!                (171171 + (-749892 + (1562275 + (570285 * x ** 2 - 1533870) * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2)
!!!      CASE(17)
!!!        fval = -0.45D2 / 0.16384D5 * sqrt(0.66D2) * dble(-143 + (19019 + (-399399 + &
!!!                (3062059 + (-10935925 + (19684665 + (5892945 * x ** 2 - 17298645) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2)
!!!      CASE(18)
!!!        fval = -0.19D2 / 0.2048D4 * sqrt(0.70D2) * dble(x) * dble(-715 + (35035 + &
!!!                (-483483 + (2877875 + (-8633625 + (13656825 + (3411705 * x ** 2 - &
!!!                10855425) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2)
!!!      CASE(19)
!!!        fval = -0.85D2 / 0.32768D5 * sqrt(0.74D2) * dble(143 + (-24024 + (644644 + &
!!!                (-6446440 + (31081050 + (-80120040 + (112896420 + (23881935 * x ** 2 - &
!!!                81880920) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
!!!                x ** 2)
!!!      CASE(20)
!!!        fval = -0.63D2 / 0.65536D5 * sqrt(0.78D2) * dble(x) * dble(7293 + (-447304 + &
!!!                (7827820 + (-60386040 + (243221550 + (-548354040 + (695987820 + &
!!!                (126233085 * x ** 2 - 463991880) * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2)
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_dVarPhi

  FUNCTION H1Basis_LegendreP(k, x) RESULT(fval)
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
!!!      CASE(15)
!!!        fval = (-0.6435D4 / 0.2048D4 + (0.255255D6 / 0.2048D4 + (-0.2909907D7 / 0.2048D4 + &
!!!                (0.14549535D8 / 0.2048D4 + (-0.37182145D8 / 0.2048D4 + (0.50702925D8 / &
!!!                0.2048D4 + (-0.35102025D8 / 0.2048D4 + 0.9694845D7 / 0.2048D4 * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
!!!      CASE(16)
!!!        fval = 0.6435D4 / 0.32768D5 + (-0.109395D6 / 0.4096D4 + (0.4849845D7 / 0.8192D4 + &
!!!                (-0.20369349D8 / 0.4096D4 + (0.334639305D9 / 0.16384D5 + (-0.185910725D9 / &
!!!                0.4096D4 + (0.456326325D9 / 0.8192D4 + (-0.145422675D9 / 0.4096D4 + &
!!!                0.300540195D9 / 0.32768D5 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) &
!!!                * x ** 2) * x ** 2) * x ** 2
!!!      CASE(17)
!!!        fval = (0.109395D6 / 0.32768D5 + (-0.692835D6 / 0.4096D4 + (0.20369349D8 / 0.8192D4 + &
!!!                (-0.66927861D8 / 0.4096D4 + (0.929553625D9 / 0.16384D5 + (-0.456326325D9 / &
!!!                0.4096D4 + (0.1017958725D10 / 0.8192D4 + (-0.300540195D9 / 0.4096D4 + &
!!!                0.583401555D9 / 0.32768D5 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
!!!      CASE(18)
!!!        fval = -0.12155D5 / 0.65536D5 + (0.2078505D7 / 0.65536D5 + (-0.14549535D8 / &
!!!                0.16384D5 + (0.156165009D9 / 0.16384D5 + (-0.1673196525D10 / 0.32768D5 + &
!!!                (0.5019589575D10 / 0.32768D5 + (-0.4411154475D10 / 0.16384D5 + &
!!!                (0.4508102925D10 / 0.16384D5 + (-0.9917826435D10 / 0.65536D5 + &
!!!                0.2268783825D10 / 0.65536D5 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
!!!      CASE(19) 
!!!        fval = (-0.230945D6 / 0.65536D5 + (0.14549535D8 / 0.65536D5 + (-0.66927861D8 / &
!!!                0.16384D5 + (0.557732175D9 / 0.16384D5 + (-0.5019589575D10 / 0.32768D5 + &
!!!                (0.13233463425D11 / 0.32768D5 + (-0.10518906825D11 / 0.16384D5 + &
!!!                (0.9917826435D10 / 0.16384D5 + (-0.20419054425D11 / 0.65536D5 + &
!!!                0.4418157975D10 / 0.65536D5 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x
!!!      CASE(20)
!!!        fval = 0.46189D5 / 0.262144D6 + (-0.4849845D7 / 0.131072D6 + (0.334639305D9 / &
!!!                0.262144D6 + (-0.557732175D9 / 0.32768D5 + (0.15058768725D11 / 0.131072D6 +&
!!!                (-0.29113619535D11 / 0.65536D5 + (0.136745788725D12 / 0.131072D6 + &
!!!                (-0.49589132175D11 / 0.32768D5 + (0.347123925225D12 / 0.262144D6 + &
!!!                (-0.83945001525D11 / 0.131072D6 + 0.34461632205D11 / 0.262144D6 * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2
    CASE DEFAULT
      ! TODO: Handle error somehow
    END SELECT
  END FUNCTION H1Basis_LegendreP

  FUNCTION H1Basis_dLegendreP(k, x) RESULT(fval)
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
!!!      CASE(15)
!!!        fval = -0.6435D4 / 0.2048D4 + (0.765765D6 / 0.2048D4 + (-0.14549535D8 / 0.2048D4 +&
!!!                (0.101846745D9 / 0.2048D4 + (-0.334639305D9 / 0.2048D4 + (0.557732175D9 / &
!!!                0.2048D4 + (-0.456326325D9 / 0.2048D4 + 0.145422675D9 / 0.2048D4 * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
!!!      CASE(16)
!!!        fval = 0.17D2 / 0.2048D4 * dble(-6435 + (285285 + (-3594591 + (19684665 + (-54679625 +&
!!!                (80528175 + (17678835 * x ** 2 - 59879925) * x ** 2) * x ** 2) * x ** 2) * &
!!!                x ** 2) * x ** 2) * x ** 2) * dble(x)
!!!      CASE(17)
!!!        fval = 0.109395D6 / 0.32768D5 + (-0.2078505D7 / 0.4096D4 + (0.101846745D9 / 0.8192D4 +&
!!!                (-0.468495027D9 / 0.4096D4 + (0.8365982625D10 / 0.16384D5 + (-0.5019589575D10 /&
!!!                0.4096D4 + (0.13233463425D11 / 0.8192D4 + (-0.4508102925D10 / 0.4096D4 + &
!!!                0.9917826435D10 / 0.32768D5 * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2
!!!      CASE(18)
!!!        fval = 0.171D3 / 0.32768D5 * dble(12155 + (-680680 + (10958948 + (-78278200 + &
!!!                (293543250 + (-619109400 + (738168900 + (119409675 * x ** 2 - 463991880) *&
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * dble(x)
!!!      CASE(19)
!!!        fval = -0.230945D6 / 0.65536D5 + (0.43648605D8 / 0.65536D5 + (-0.334639305D9 / 0.16384D5 +&
!!!                (0.3904125225D10 / 0.16384D5 + (-0.45176306175D11 / 0.32768D5 + &
!!!                (0.145568097675D12 / 0.32768D5 + (-0.136745788725D12 / 0.16384D5 + &
!!!                (0.148767396525D12 / 0.16384D5 + (-0.347123925225D12 / 0.65536D5 + &
!!!                0.83945001525D11 / 0.65536D5 * x ** 2) * x ** 2) * x ** 2) * x ** 2) *&
!!!                x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2
!!!    CASE(20)
!!!      fval = 0.105D3 / 0.65536D5 * dble(-46189 + (3187041 + (-63740820 + (573667380 + &
!!!              (-2772725670 + (7814045070 + (-13223768580 + (13223768580 + (1641030105 *&
!!!              x ** 2 - 7195285845) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * x ** 2) * &
!!!              x ** 2) * x ** 2) * x ** 2) * dble(x)
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
      grad(j,1,1) = -1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,1) = 1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,1) = 0D0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,2) = -1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,2) = 0D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,2) = 1D0
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
      grad(j,1,1) = -1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,1) = 1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,1) = 0D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,4,1) = 0D0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,2) = -1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,2) = 0D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,2) = 1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,4,2) = 0D0
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,3) = -1D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,2,3) = 0D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,3,3) = 0D0
    END DO
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,4,3) = 1D0
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
    REAL(Kind=dp), PARAMETER :: c = 1D0/2D0
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
    REAL(Kind=dp), PARAMETER :: c = 1D0/2D0
    INTEGER :: j

    ! First coordinate (xi)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,1) = -c*(1-w(j))
      grad(j,2,1) = c*(1-w(j))
      grad(j,3,1) = 0D0
      grad(j,4,1) = -c*(1+w(j))
      grad(j,5,1) = c*(1+w(j))
      grad(j,6,1) = 0D0
    END DO
    
    ! Second coordinate (eta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,2) = -c*(1-w(j))
      grad(j,2,2) = 0D0
      grad(j,3,2) = c*(1-w(j))
      grad(j,4,2) = -c*(1+w(j))
      grad(j,5,2) = 0D0
      grad(j,6,2) = c*(1+w(j))
    END DO
    
    ! Third coordinate (zeta)
    !_ELMER_OMP_SIMD
    DO j=1,nvec
      grad(j,1,3) = -c*(1-u(j)-v(j))
      grad(j,2,3) = -c*u(j)
      grad(j,3,3) = -c*v(j)
      grad(j,4,3) = c*(1-u(j)-v(j))
      grad(j,5,3) = c*u(j)
      grad(j,6,3) = c*v(j)
    END DO
  END SUBROUTINE H1Basis_dWedgeNodal

END MODULE H1Basis
