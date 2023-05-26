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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Mika  Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 Dec 2003
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  A consistent splitting scheme for the incompressible Navier-Stokes equations:
!>  The computation of viscous compressibility term
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE CompressibilitySolver( Model,Solver,dt,TransientSimulation )
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
  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Found, OutflowBC, &
      NewTimeStep, ConstantBulkMatrix, ConstantBulkMatrixInUse

  TYPE(Element_t),POINTER :: Element

  INTEGER :: i, j, n, nb, nd, t, istat, dim, active, CurrentDoneTime = 0
  REAL(KIND=dp) :: Norm = 0, PrevNorm, RelC, PressureRelax, Visc

  TYPE(Variable_t), POINTER :: PressureVariable
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), FORCE(:), &
            cvelo(:,:),pvelo(:,:), ppres(:,:), MASS(:,:)

  INTEGER, ALLOCATABLE :: Indexes(:)

  TYPE(Variable_t), POINTER :: VeloVar
  TYPE(ValueList_t), POINTER :: BC, SolverParams

  SAVE MASS, STIFF, LOAD, FORCE, AllocationsDone, CVelo,PVelo, ppres, indexes, &
      CurrentDoneTime
!------------------------------------------------------------------------------

  dim = CoordinateSystemDimension()
  Mesh => GetMesh()
  SolverParams => GetSolverParams()
  ! Mesh => Solver % Mesh

  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', Found )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
       ASSOCIATED(Solver % Matrix % BulkValues)

 
  VeloVar => VariableGet( Mesh % Variables, "VelocityTot" )

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( MASS(n,n), FORCE(N), STIFF(N,N), CVelo(dim,N), &
             indexes(n), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'CompressibilitySolver', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

  ! Initialize the system and do the assembly:
  !-------------------------------------------
  Active = GetNOFActive()

  CALL DefaultInitialize(Solver, ConstantBulkMatrixInUse)

  DO t=1,Active      !Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes()
     nb = GetElementNOFBDOFs()
     nd = GetElementDOFs( Indexes )

     ! Get previous elementwise velocity iterate:
     !-------------------------------------------
     DO i=1,dim
        CVelo(i,1:nd) = VeloVar % Values( &
             VeloVar % DOFs*(VeloVar % Perm(Indexes(1:nd))-1)+i)
     END DO

     Visc = ListGetConstReal( GetMaterial(), 'Viscosity' )

     ! Get element local matrix and rhs vector:
     !-----------------------------------------
     CALL LocalMatrix(  MASS, STIFF, FORCE,  CVelo, &
            Element, n, nd, nd+nb, dim )

     ! Update global matrix and rhs vector from local matrix & vector:
     !----------------------------------------------------------------
     IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        CALL DefaultUpdateEquations( STIFF, FORCE )
     ELSE
        CALL DefaultUpdateForce( FORCE ) 
     END IF
        
  END DO

  IF ( ConstantBulkMatrix ) THEN
    CALL DefaultFinishBulkAssembly(BulkUpdate = .NOT.ConstantBulkMatrixInUse, &
        RHSUpdate = .FALSE.)
  ELSE
    CALL DefaultFinishBulkAssembly()
  END IF

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()
  
!------------------------------------------------------------------------------

  ! Solve the system:
  !------------------
  Norm = DefaultSolve()


!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  MASS, STIFF, FORCE, NodalVelo, Element, n, nd, ntot, dim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: NodalVelo(:,:)
    INTEGER :: dim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3), &
                  DetJ, divU, du(3), DPres(3), AverDiv, Ddivu
    REAL(KIND=dp), POINTER :: A(:,:), F(:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS  = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    AverDiv = 0.0d0
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       ! Divergence at the integration point:
       !--------------------------------------------
       divu = 0.0d0

       DO i=1,dim
          divu = divu + SUM( Cvelo(i,1:nd) * dBasisdx(1:nd,i) )
       END DO

       ! Finally, the elemental matrix & vector:
       !----------------------------------------
       DO p=1,n
          IF ( .NOT. ConstantBulkMatrixInUse ) THEN
             DO q=1,n
                stiff(p,q) = stiff(p,q) + s * Basis(q) * Basis(p)
             END DO
          END IF
          !FORCE(2*(p-1)+1) = FORCE(2*(p-1)+1) + s * 2.0d0 * visc * divu * Basis(p)
          !FORCE(2*(p-1)+2) = FORCE(2*(p-1)+2) + s * SUM( Du(1:dim) * dBasisdx(p,1:dim) )
          FORCE(p) = FORCE(p) + s * 2.0d0 * visc * divu * Basis(p)
       END DO

    END DO


   ! Eliminate the bubble and edge degrees of freedom if any
    
   DO i = n+1,ntot
      FORCE(i)   = 0.0d0
      STIFF(i,:) = 0.0d0
      STIFF(:,i) = 0.0d0
      STIFF(i,i) = 1.0d0      
   END DO



!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
END SUBROUTINE CompressibilitySolver
!------------------------------------------------------------------------------
