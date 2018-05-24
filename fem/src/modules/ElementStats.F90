!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Module Authors: Peter Råback
! *  Email:          Peter.Raback@csc.fi
! *  Web:            http://www.csc.fi/elmer
! *  Address:        CSC - IT Center for Science Ltd.
! *
! *  Original Date: 2.8.2016
! *
! *****************************************************************************/
!-----------------------------------------------------------------------------
!> A solver that reports the staticstics of the finite element mesh.
!------------------------------------------------------------------------------
SUBROUTINE ElementStats( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: i, j, t, t0, n, m, is_bc
  INTEGER :: FirstElem, LastElem

  REAL(KIND=dp) :: WarnSize, MinF(3), MaxF(3), AveF(3), F, dF
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: ElementCount(1000)
  CHARACTER(MAX_NAME_LEN) :: Str, OperName
  LOGICAL :: Found, DoCohorts, DoingCohorts
  INTEGER :: NoCohorts, OperNo
  INTEGER, ALLOCATABLE :: CohortCount(:,:)
!------------------------------------------------------------------------------
  
  CALL Info('ElementStats','------------------------------',Level=6)
  CALL Info('ElementStats','Studying elements of the mesh')

  Params => GetSolverParams()
  Mesh => GetMesh()

  ! Do cohorst only if requested
  DoCohorts = GetLogical( Params,'Create Histogram',Found ) 
  IF( DoCohorts ) THEN
    NoCohorts = GetInteger( Params,'Histogram Intervals',Found ) 
    IF( .NOT. Found ) NoCohorts = 10 
    ALLOCATE( CohortCount( NoCohorts, 3 ) )
    CohortCount = 0
  END IF

  ! Go through bulk and boundary elements separately
  DO is_bc = 0, 1

    IF( is_bc == 0 ) THEN
      t0 = 0
      FirstElem = 1
      LastElem = Mesh % NumberOfBulkElements 
      Str = 'bulk'
    ELSE
      t0 = Mesh % NumberOfBulkElements
      FirstElem = 1
      LastElem = Mesh % NumberOfBoundaryElements 
      Str = 'boundary'
    END IF
    
    IF( FirstElem > LastElem ) CYCLE
    
    CALL Info('ElementStats','Going through '//TRIM(Str)//' elements')
    
    ! Initialize separately for bulk and boundary
    WarnSize = EPSILON( WarnSize ) 
    MinF = HUGE( MinF )
    MaxF = -HUGE( MaxF ) 
    AveF = 0.0_dp
    DoingCohorts = .FALSE.
    ElementCount = 0

    
100 DO t=FirstElem, LastElem
      
      Element => Mesh % Elements(t0+t)
      Model % CurrentElement => Element
      
      n = GetElementNOFNodes(Element)
      m = GetElementCode(Element)

      IF(.NOT. DoingCohorts ) THEN
        ElementCount(m) = ElementCount(m) + 1
      END IF
              
      DO OperNo = 1,3 

        SELECT CASE(OperNo)
        CASE( 1 ) 
          CALL ElementSizeStudy( Element, n, f )
          
        CASE(2)
          CALL ElementSkewStudy( Element, n, f )
          
        CASE(3)
          CALL ElementRatioStudy( Element, n, f )
          
        END SELECT
        
        IF( .NOT. DoingCohorts ) THEN        
          MinF( OperNo ) = MIN( f, MinF( OperNo ) ) 
          MaxF( OperNo ) = MAX( f, MaxF( OperNo ) )
          AveF( OperNo ) = f + AveF( OperNo )         
        ELSE
          i = CEILING( NoCohorts * (f-MinF(OperNo))/(MaxF(OperNo)-MinF(OperNo)) )
          i = MAX( 1, MIN( i, NoCohorts ) )
          CohortCount(i,OperNo) = CohortCount(i,OperNo) + 1
        END IF
      END DO
    END DO

    ! If we study cohorst we need a second sweep over the elements
    IF( DoCohorts .AND. .NOT. DoingCohorts ) THEN
      DoingCohorts = .TRUE.
      GOTO 100
    END IF

    IF( LastElem > FirstElem ) THEN
      AveF = AveF / ( LastElem - FirstElem + 1 ) 
    END IF
      
    
    DO OperNo = 1,3     
      SELECT CASE(OperNo)
      CASE( 1 ) 
        OperName = 'Element size'
        
      CASE(2)
        OperName = 'Element skew (degs)'
        
      CASE(3)
        OperName = 'Element ratio'
      END SELECT

      PRINT *,'Statistics for mesh operator: '//TRIM(OperName)
      PRINT *,'Minimum values: ',MinF(OperNo)
      PRINT *,'Maximum values: ',MaxF(OperNo)
      PRINT *,'Average values: ',AveF(OperNo)
      
      IF( OperNo == 1 ) THEN
        PRINT *,'Ratio of element sizes: ',MaxF(OperNo) / MinF(OperNo)
      END IF

      IF( DoCohorts ) THEN        
      PRINT *,'Histogram:'
        dF = ( MaxF(OperNo) - MinF(OperNo) ) / NoCohorts
        DO i=1,NoCohorts
          IF( CohortCount(i,OperNo) > 0 ) THEN
            PRINT *,i,': ',CohortCount(i,OperNo),' (',MinF(OperNo)+(i-1)*dF,' to ',MinF(OperNo)+i*dF,')'
          END IF
        END DO
      END IF
    END DO
    
  END DO

  CALL Info('ElementStats','------------------------------',Level=6)


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE ElementSizeStudy( Element, n, ElemSize )
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: ElemSize
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight
    REAL(KIND=dp) :: Basis(n),DetJ
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------


    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )

    ElemSize = 0.0_dp

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ
      ElemSize = ElemSize + Weight
    END DO

    IF( ElemSize < EPSILON( ElemSize ) ) THEN
      PRINT *,'Element too small: ',Element % ElementIndex, ElemSize
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE ElementSizeStudy
!------------------------------------------------------------------------------



  !------------------------------------------------------------------------------
  SUBROUTINE ElementRatioStudy( Element, n, q )
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: q
!------------------------------------------------------------------------------
    INTEGER :: ElemFamily
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER :: i1,i2,k
    REAL(KIND=dp) :: ri(3), s2, s2min, s2max
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    
    CALL GetElementNodes( Nodes )
    ElemFamily = GetElementFamily( Element )    
    EdgeMap => GetEdgeMap( ElemFamily )
    
    s2max = 0.0_dp
    s2min = HUGE( s2min ) 
    
    DO k=1,Element % Type % NumberOfEdges
      i1 = EdgeMap(k,1)
      i2 = EdgeMap(k,2)
      
      ri(1) = Nodes % x(i2) - Nodes % x(i1)
      ri(2) = Nodes % y(i2) - Nodes % y(i1)
      ri(3) = Nodes % z(i2) - Nodes % z(i1)
      
      s2 = SQRT( SUM( ri * ri ) )

      s2max = MAX( s2max, s2 )
      s2min = MIN( s2min, s2 )
    END DO

    ! Compute the elementanl max ratio
    q = SQRT( s2max / s2min )
    
!------------------------------------------------------------------------------
  END SUBROUTINE ElementRatioStudy
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ElementSkewStudy( Element, n, alpha )
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: alpha
!------------------------------------------------------------------------------
    INTEGER :: ElemFamily
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER :: NoEdgePairs,i1,i2,j1,j2,k
    INTEGER :: EdgePairs(24)
    REAL(KIND=dp) :: Beta, MaxBeta, ri(3), rj(3), ri_len, rj_len
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    
    CALL GetElementNodes( Nodes )
    ElemFamily = GetElementFamily( Element )    
    EdgeMap => GetEdgeMap( ElemFamily )
    

    Alpha = 0.0_dp

    
    SELECT CASE ( ElemFamily )
      
    CASE(4) ! quads
      NoEdgePairs = 2
      EdgePairs(1:2*NoEdgePairs) = (/1,3,2,4/)
      ! I.e. the edge pairs to study are 1 & 3, and 2 & 4
      
    CASE(6) ! pyramids
      NoEdgePairs = 2
      EdgePairs(1:2*NoEdgePairs) = (/1,3,2,4/)
      
    CASE(7) ! wedges/prisms
      NoEdgePairs = 6
      EdgePairs(1:2*NoEdgePairs) = (/1,4,7,8,2,5,8,9,3,6,9,7/)
      
    CASE(8) ! hexas
      EdgePairs = 12
      EdgePairs(1:2*NoEdgePairs) = &
          (/1,3,2,4,5,7,6,8,1,5,9,10,2,6,10,11,3,7,11,12,4,8,12,9/)      

    CASE DEFAULT
      RETURN
            
    END SELECT
        
    MaxBeta = 0.0_dp
    
    DO k=1,NoEdgePairs
      i1 = EdgeMap(EdgePairs(2*k-1),1)
      i2 = EdgeMap(EdgePairs(2*k-1),2)
      j1 = EdgeMap(EdgePairs(2*k),1)
      j2 = EdgeMap(EdgePairs(2*k),2)
      
      ri(1) = Nodes % x(i2) - Nodes % x(i1)
      ri(2) = Nodes % y(i2) - Nodes % y(i1)
      ri(3) = Nodes % z(i2) - Nodes % z(i1)
      
      rj(1) = Nodes % x(j2) - Nodes % x(j1)
      rj(2) = Nodes % y(j2) - Nodes % y(j1)
      rj(3) = Nodes % z(j2) - Nodes % z(j1)
      
      ri_len = SQRT( SUM( ri * ri ) )
      rj_len = SQRT( SUM( rj * rj ) )
      
      IF( ri_len * rj_len > 0 ) THEN
        beta = ACOS( SUM( ri * rj ) / ( ri_len * rj_len ) )
        beta = MIN( ABS( beta ), ABS( PI - beta ) )
        MaxBeta = MAX( beta, MaxBeta )
      END IF
    END DO

    ! Transform angle to degrees
    Alpha = 180 * MaxBeta / PI 
    
!------------------------------------------------------------------------------
  END SUBROUTINE ElementSkewStudy
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE ElementStats
!------------------------------------------------------------------------------
