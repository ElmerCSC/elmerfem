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
! *  Authors: Peter Råback, Juha Ruokolainen, Eelis Takala
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 5.9.2013
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
!> Subroutine that helps to give current sources for closed coils. Normally there
!> is the difficulty that the potential in such a coil should be discontinuous. 
!> In this solver two different potential fields are solved (say a left and a right 
!> field). The two fields have different boundary conditions that are automatically 
!> set on the bulk nodes so that a current is induced. The solution settles typically
!> on the other side such that the union of the solutions is always reasonable.
!> The main assumption is that the coil axis is aligned with the z-axis. 
!
!> The module includes a special user function that returns the potential 
!> at the good side. 
!
!> \ingroup Solvers
!------------------------------------------------------------------------------

SUBROUTINE WVectorFix_init( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: Params 
  INTEGER :: dim
  LOGICAL :: Found, CalcCurr, CalculateElemental, CalculateNodal

  dim = CoordinateSystemDimension()
  Params => GetSolverParams()

  CALL ListAddNewString( Params,'Variable','-nooutput WFixVarTmp')

END SUBROUTINE WVectorFix_init

!------------------------------------------------------------------------------
SUBROUTINE WVectorFix( Model,Solver,dt,TransientSimulation )
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
  REAL(KIND=dp) :: Norm!, s
  INTEGER :: active, dim, n, nd, nsize, t, j, i, varindex!,k, nb, iter, Part, sgn, CoilParts, &
!       MaxNonlinIter,dimi,ierr, NoCoils, MaxNoCoils
  INTEGER :: Dnode_index
  INTEGER, POINTER :: Perm(:)!, Set(:),TargetBodies(:)
!  INTEGER, ALLOCATABLE, TARGET :: SetA(:), SetB(:)
  TYPE(Matrix_t), POINTER :: StiffMatrix
  REAL(KIND=dp), POINTER :: ForceVector(:)
  TYPE(Variable_t), POINTER :: SolVar, GradV!, PotVar, FixVar, FluxVar, LoadVar, DistVar
!  TYPE(Variable_t), POINTER :: PotVarA,PotVarB,PotSelect,CoilIndexVar,CoilSetVar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params!, CoilList 
  TYPE(Valuelist_t), POINTER :: Material
  CHARACTER(LEN=MAX_NAME_LEN) :: FluxVarname(2)
!  REAL(KIND=dp) :: CoilCrossSection,InitialCurrent, Coeff, val, x0
!  REAL(KIND=dp), ALLOCATABLE :: DesiredCoilCurrent(:), DesiredCurrentDensity(:)
!  LOGICAL :: Found, CoilClosed, CoilAnisotropic, UseDistance, FixConductivity, &
!      NormalizeCurrent, FitCoil, SelectNodes, CalcCurr, NarrowInterface
!  LOGICAL, ALLOCATABLE :: GotCurr(:), GotDens(:)
!  REAL(KIND=dp) :: CoilCenter(3), CoilNormal(3), CoilTangent1(3), CoilTangent2(3), &
!      MinCurr(3),MaxCurr(3),TmpCurr(3)
!  INTEGER, ALLOCATABLE :: CoilIndex(:)
!  CHARACTER(LEN=MAX_NAME_LEN) :: CondName

 !------------------------------------------------------------------------------

  CALL Info('WVectorFix','----------------------------------------')
  CALL Info('WVectorFix','Fixing the coil gradient potential field')
  CALL Info('WVectorFix','----------------------------------------')

  Params => GetSolverParams()

  nsize = SIZE( Solver % Variable % Values ) 
  StiffMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % Rhs
  Mesh => Solver % Mesh
  SolVar => Solver % Variable
  Perm => SolVar % Perm

  dim = CoordinateSystemDimension()

  FluxVarname(1) = 'W Vector E'
  FluxVarname(2) = 'J Vector E'

  DO varindex = 1, 2
    GradV => VariableGet( Mesh % Variables,FluxVarname(varindex))  
    IF( .NOT. ( ASSOCIATED( GradV ) ) ) THEN
      IF (varindex .EQ. 2) CYCLE
      CALL Fatal('WVectorFix','Fixing is done currently only for '//TRIM(FluxVarname(varindex))//&
        '! I cannot find it.')
    END IF    

    CALL INFO('WVectorFix','Applying elemental div()=0 to: '//TRIM(FluxVarname(varindex)),level=3)
    
    CALL DefaultInitialize()
    Active = GetNOFActive() 
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()           
      CALL LocalCorrMatrix(  Element, n, nd )
    END DO
    IF( Parenv%myPE == parenv % pes ) Dnode_index = Element % NodeIndexes(1)

    ! Set just one Dirichlet node.
    ! In serial it can be as well the 1st one. 
    !CALL UpdateDirichletDof(Solver % Matrix, Dnode_index, 0._dp)        
      
    ! Only Default dirichlet conditions activate the BCs above!
    !CALL DefaultDirichletBCs()
    
    Norm = DefaultSolve()
    
    CALL Info('WVectorFix','Fixing '//TRIM(FluxVarname(varindex))//' to be divergence free!')
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()           
      CALL LocalCorrCurrent(  Element, n, nd )
    END DO
  END DO

CONTAINS 

  ! Assembly the potential equation related to Jfix field with divergence of
  ! elemental coil current as source.
  !-------------------------------------------------------------------------
  SUBROUTINE LocalCorrMatrix( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,L(3),Lfix(3,nd),Tcoef(nd)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), C
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: Material
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    Material => GetMaterial( Element )
    IF (.NOT. ASSOCIATED(Material)) CALL Fatal('WVectorFix','Material not found.')

    Tcoef(1:n) = GetReal(Material,'Electric Conductivity',Found)

    ! Current density at nodes
    DO i=1,n
      j = GradV % Perm( Element % DGIndexes(i) )      
      Lfix(1:3,i) = GradV % Values( 3*j-2: 3*j )
    END DO
    
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )                
      Weight = IP % s(t) * detJ

      ! Current density at integration point
      L = MATMUL(Lfix(1:3,1:n),Basis(1:n))
      C = MAXVAL((/1._dp, SUM( Tcoef(1:n) * Basis(1:n) )/))
      
      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + Weight * C * SUM(dBasisdx(q,:)*dBasisdx(p,:))
        END DO
        FORCE(p) = FORCE(p) + SUM(L * dBasisdx(p,:)) * Weight
      END DO
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalCorrMatrix
!------------------------------------------------------------------------------

  ! Correct the current elementwise to be divergence free.
  ! This is done by solving elementwise the fixing currents.
  !------------------------------------------------------------------------------
  SUBROUTINE LocalCorrCurrent( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),NodalPot(nd),DetJ,Weight,Tcoef(nd)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd,3),x(nd),C
    INTEGER :: pivot(nd)
    LOGICAL :: Stat,Found,Erroneous
    INTEGER :: i,j,k,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    Material => GetMaterial( Element )
    IF (.NOT. ASSOCIATED(Material)) CALL Fatal('WVectorFix','Material not found.')

    Tcoef(1:n) = GetReal(Material,'Electric Conductivity',Found)

    NodalPot(1:n) = Solver % Variable % Values( Perm( Element % NodeIndexes ) )             
    
    ! Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )                
      Weight = IP % s(t) * detJ

      C = MAXVAL((/1._dp, SUM( Tcoef(1:n) * Basis(1:n) )/))
      

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(p) * Basis(q) 
        END DO
        ! Source for each coordinate direction
        DO k=1,3
          FORCE(p,k) = FORCE(p,k) + Weight * Basis(p) * SUM(dBasisdx(1:n,k)*NodalPot(1:n))
        END DO
      END DO
    END DO


    CALL LUdecomp(STIFF,n,pivot,Erroneous)
    IF (Erroneous) CALL Fatal('WVectorFix', 'LU-decomposition fails')
    ! Fixing for each coordinate direction
    DO k=1,3
      x = FORCE(1:n,k)
      CALL LUSolve(n,STIFF,x,pivot)

      DO i=1,n
        j = GradV % Perm( Element % DGIndexes(i) )              
        GradV % Values( 3*(j-1)+k) = GradV % Values( 3*(j-1)+k) - x(i)
      END DO
    END DO    
    
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalCorrCurrent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE WVectorFix
!------------------------------------------------------------------------------
