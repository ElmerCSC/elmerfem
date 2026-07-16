!/*****************************************************************************
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
! * A simple solver to regulate velocity such that the surface velocity never
! * points out of the surface. 
! *
! *****************************************************************************/

SUBROUTINE FixTangentVelo_init( Model,Solver,dt,Transient )
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL ListAddNewLogical( Solver % Values,'No Matrix',.TRUE.)
  CALL ListAddNewString( Solver % Values,'Variable',&
      '-nooutput -global FixTangent_var')
   
END SUBROUTINE FixTangentVelo_Init
  


SUBROUTINE FixTangentVelo( Model,Solver,dt,Transient )
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!---------------------------------------------------------------
  INTEGER :: n,i,j,k,dim
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: VeloVar, FixVar
  CHARACTER(LEN=MAX_NAME_LEN) :: str
  INTEGER, POINTER :: FixPerm(:)
  REAL(KIND=dp), POINTER :: FixVals(:)
  LOGICAL :: Visited=.FALSE., Found, Top, DownstreamOnly 
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC, Params
  INTEGER :: i1,i2,t,ActiveCoord
  REAL(KIND=dp) :: e1(3), e2(3), Amat(3,3), velem(4,4), v(3), c(3), Normal(3), Nrm
  
  SAVE Visited, Nodes, FixVar, VeloVar, FixVals, FixPerm, ActiveCoord, DownstreamOnly

  
  CALL Info('FixTangentVelo','Regulating surface velocity not to go out domain')

  IF(.NOT. Visited) THEN
    Mesh => Solver % Mesh
    Params => Solver % Values

    dim = CoordinateSystemDimension()
    IF(dim /= 3 ) CALL Fatal('FixTangentVelo','Implemented only for 3D systems!')

    ActiveCoord = ListGetInteger(Params,'Active Coordinate',Found)
    IF(.NOT. Found) ActiveCoord = dim

    DownstreamOnly = ListGetLogical(Params,'Fix Downstream Only',Found )
    
    str = ListGetString(Params,'Velocity Variable Name',Found)
    IF(.NOT. Found) str = 'flow solution'
    VeloVar => VariableGet( Mesh % Variables, str )
    IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
      CALL Fatal('FixTangentVelo','Velocity field variable does not exist: '//TRIM(str))           
    END IF

    ALLOCATE(FixPerm(Mesh % NumberOfNodes) )
    FixPerm = 0
    DO t=Mesh % NumberOfBulkElements+1,&
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      BC => GetBC(Element)   
      IF(.NOT. ASSOCIATED(BC)) CYCLE
      IF( ListGetLogical(BC,'Tangent Surface',Found) ) THEN
        FixPerm(Element % NodeIndexes) = 1
      END IF
    END DO

    n = 0
    DO i=1,Mesh % NumberOfNodes
      IF(FixPerm(i) > 0 ) THEN
        n = n+1
        FixPerm(i) = n
      END IF
    END DO

    n = ParallelReduction(n)
    CALL Info('FixTangentVelo','Total number of nodes to fix: '//I2S(n))
    IF(n==0) CALL Fatal('FixTangentVelo','No nodes to fix!')
    
    CALL VariableAddVector( Mesh % Variables, Mesh, Solver,'FixVelo', Perm = FixPerm )
    FixVar => VariableGet( Mesh % Variables,'FixVelo')    
    FixVals => FixVar % Values
    Visited = .TRUE.        
  END IF


  FixVals = 0.0_dp

  DO t=Mesh % NumberOfBulkElements+1,&
      Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    Element => Mesh % Elements(t)
    BC => GetBC(Element)   
    IF(.NOT. ASSOCIATED(BC)) CYCLE
    Top = ListGetLogical(BC,'Tangent Surface',Found)
    IF(.NOT. Top) CYCLE

    CALL GetVectorLocalSolution(velem,UElement=Element,UVariable=VeloVar)
    CALL GetElementNodes(Nodes, Element )
    Normal = NormalVector(Element, Nodes, Check=.TRUE. )
     
    n = Element % TYPE % NumberOfNodes
    IF(.NOT. (n==3 .OR. n==4)) THEN
      CALL Fatal('FixTangentVelo','Only implemented for 3 and 4 nodes!')
    END IF
    
    ! Go through each element node. So lazy coding here fetching the same stuff many times. 
    DO i=1,n
      v = velem(1:3,i)
      i1 = MODULO(i,n)+1
      i2 = MODULO(i-2,n)+1
                 
      ! Edges point to two other nodes
      e1(1) = Nodes % x(i1) - Nodes % x(i)
      e1(2) = Nodes % y(i1) - Nodes % y(i)
      e1(3) = Nodes % z(i1) - Nodes % z(i)
      
      e2(1) = Nodes % x(i2) - Nodes % x(i)
      e2(2) = Nodes % y(i2) - Nodes % y(i)
      e2(3) = Nodes % z(i2) - Nodes % z(i)
      
      ! The flow needs to be towards either node that we can have an issue.
      ! It is enough that either nodes gives postive dot product. 
      IF( DownstreamOnly ) THEN
        IF(SUM(v*e1) < 0 .AND. SUM(v*e2) < 0 ) CYCLE
      END IF
        
      Amat = 0.0_dp
      Amat(:,1) = e1
      Amat(:,2) = e2
      Amat(ActiveCoord,3) = 1.0_dp

      CALL SolveLinSys3x3( Amat, c, v ) 
      
      j = FixPerm(Element % NodeIndexes(i))
      IF( Normal(ActiveCoord) > 0.0_dp ) THEN
        FixVals(j) = MAX(FixVals(j),c(3))
      ELSE
        FixVals(j) = MIN(FixVals(j),c(3))
      END IF
    END DO
  END DO
  
  IF(InfoActive(20) ) THEN
    CALL VectorValuesRange(FixVals, SIZE(FixVals),'FixVals')
  END IF
  
  
  IF( ListGetLogical(Params,'Fix Tangent Velocity',Found ) ) THEN
    CALL ApplyTangentFix(.FALSE.)
  END IF

  Nrm = ComputeNorm(Solver,SIZE(FixVals),FixVals)
  Solver % Variable % Values = Nrm
  Solver % Variable % Norm = Nrm
  
  
CONTAINS

  SUBROUTINE ApplyTangentFix(UndoFix)
    LOGICAL :: UndoFix
    INTEGER :: sgn,i,jv,jf,vdofs

    CALL Info('FixTangentVelo','Applying the fix to current velocity field')
    
    sgn = -1
    IF(UndoFix) sgn=1

    vdofs = VeloVar % dofs        
    DO i=1, Mesh % NumberOfNodes
      jf = FixVar % Perm(i)
      jv = VeloVar % Perm(i)
      IF(jf==0 .OR. jv==0) CYCLE
      VeloVar % Values(vdofs*(jv-1)+ActiveCoord) = VeloVar % Values(vdofs*(jv-1)+ActiveCoord) +&
          sgn * FixVar % Values(jf)
    END DO
    
  END SUBROUTINE ApplyTangentFix
 
  
END SUBROUTINE FixTangentVelo
