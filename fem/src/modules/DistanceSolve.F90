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
!/*****************************************************************************/
!  
!    Authors: Juha Ruokolainen
!    Email:   Juha.Ruokolainen@csc.fi
!    Web:     http://www.csc.fi/elmer
!    Address: CSC - IT Center for Science Ltd.
!             Keilaranta 14
!             02101 Espoo, Finland
!  
!    Original Date: 09 Nov 2007
!
!/*****************************************************************************/


!------------------------------------------------------------------------------
!> Solves the equation: (grad(d),grad(d))=1. The solution of this 
!> gives the closest distance a to boundary where distance is forced to zero.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE DistanceSolver( Model,Solver,dt,TransientSimulation )
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
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element

  TYPE(ValueList_t), POINTER :: SolverParams, BC

  REAL(KIND=dp) :: Pnorm,Norm,RelaxDT,x0,y0,z0,x1,y1,z1
  INTEGER :: n, nb, nd, t, i,j,k,istat, active, MaxIter
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), buf(:)

  SAVE STIFF, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  CALL Info('DistanceSolver','Using PDE based distance solver')


  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), STIFF(N,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF

     n = Mesh % NumberOfNodes
     IF ( ALL( Solver % Variable % Values == 0.0_dp ) ) &
       Solver % Variable % Values(Solver % Variable % Perm(1:n)) = &
           SQRT( Mesh % Nodes % x(1:n)**2 + &
                 Mesh % Nodes % y(1:n)**2 + &
                 Mesh % Nodes % z(1:n)**2 )

     AllocationsDone = .TRUE.
  END IF

  SolverParams => GetSolverParams()
  RelaxDT = GetCReal( SolverParams, 'Distance Pseudo dt', Found )

  MaxIter = GetInteger( SolverParams, 'Nonlinear System Max Iterations', Found )
  IF ( .NOT. Found ) MaxIter = 100

  DO i=1,Model % NumberOFBCs
    BC => Model % BCs(i) % Values
    IF ( GetLogical( BC, 'Noslip Wall BC', Found ) ) &
      CALL ListAddConstReal( BC, ComponentName(Solver % Variable), 0.0_dp )
  END DO

  DO i=1,MaxIter
    !System assembly:
    !----------------
    Active = GetNOFActive()
    CALL DefaultInitialize()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, Element, n, nd+nb )
      CALL CondensateP( nd, nb, STIFF, FORCE )
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
    CALL DefaultFinishBulkAssembly()

    ! No Flux BCs
    CALL DefaultFinishAssembly()

    CALL DefaultDirichletBCs()

    ! And finally, solve:
     !--------------------
    Norm = DefaultSolve()
    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

  CALL Info('DistanceSolver','All done')


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ, &
            SOL(nd),Grad(3),SU(nd),stab,uabs,hscale
    LOGICAL :: Stat
    INTEGER :: i,j,t,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim  = CoordinateSystemDimension()
    CALL GetScalarLocalSolution( SOL )

    hscale = GetCReal( GetSolverParams(), 'H scale', stat )

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )

      DO i=1,dim
        Grad(i) = SUM( dBasisdx(1:nd,i) * SOL(1:nd) )
      END DO

      Uabs = SQRT( SUM( Grad(1:dim)**2 ) )
      Grad(1:dim) = Grad(1:dim) / MAX(Uabs,1.d-8)
      stab = hscale*Element % hk/2

      SU = 0.0_dp
      DO i=1,nd
        DO j=1,dim
          SU(i) = SU(i) + Grad(j) * dBasisdx(i,j)
        END DO
      END DO

      DO i=1,nd
      DO j=1,nd
         STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ * &
           SUM( Grad(1:dim) * dBasisdx(j,1:dim) ) * Basis(i)

         STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ * stab * SU(i) * SU(j)

         IF ( RelaxDT > 0 ) STIFF(i,j) = STIFF(i,j) + &
          1._dp/RelaxDT*IP % s(t) * detJ * Basis(j) * Basis(i)
      END DO
      END DO
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (Basis(1:nd)+stab*SU(1:nd))

      IF ( RelaxDT > 0 ) FORCE(1:nd) = FORCE(1:nd) + &
         1._dp/RelaxDT*IP % s(t) * DetJ * Basis(1:nd)*SUM(Basis(1:nd)*SOL(1:nd))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE DistanceSolver
!------------------------------------------------------------------------------

!--------------------------------------------------------------------
!> Initialization for the primary solver
!--------------------------------------------------------------------
SUBROUTINE DistanceSolver_init( Model,Solver,dt,TransientSimulation )
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found

  SolverParams => GetSolverParams()

  IF( .NOT. ListCheckPresent(SolverParams,'Nonlinear System Convergence Tolerance') ) THEN
    CALL ListAddConstReal(SolverParams,'Nonlinear System Convergence Tolerance',1.0d-8) 
  END IF	

END SUBROUTINE DistanceSolver_init	



!------------------------------------------------------------------------------
!> Geometric variant of the DistanceSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE DistanceSolver1( Model,Solver,dt,TransientSimulation )
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
  LOGICAL ::  Found, Found1, Found2
  TYPE(Element_t),POINTER :: Element

  TYPE(ValueList_t), POINTER :: SolverParams, BC, BodyForce

  INTEGER :: n, m, nb, nnb, nd, t, i,j,k,istat, active, MaxIter, maxnode,comm
  REAL(KIND=dp) :: Pnorm,Norm,RelaxDT,TOL,x0,y0,z0,dist
  TYPE(Mesh_t), POINTER :: Mesh

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), buf(:) ,xp(:),yp(:),zp(:),bdist(:),bd(:),&
      condition(:), work(:)
  REAL(KIND=dp), POINTER :: distance(:)
  INTEGER, POINTER :: gPerm(:), ibuf(:), aperm(:),bperm(:),cperm(:)
  LOGICAL :: DummyDistance 

  SAVE STIFF, FORCE
!------------------------------------------------------------------------------

  CALL Info('DistanceSolver1','-------------------------------------------',Level=5)
  CALL Info('DistanceSolver1','Using geometric distance solver',Level=5)
  CALL Info('DistanceSolver1','-------------------------------------------',Level=5)
  CALL ResetTimer('DistanceSolver1')

  Mesh => GetMesh()

  CALL Info('DistanceSolver1','Working with mesh: '&
      //TRIM(Mesh % Name),Level=10)
  CALL Info('DistanceSolver1','Solving for variable: '&
      //TRIM(Solver % variable % Name ),Level=10)

  n = Mesh % NumberOfNodes
  ALLOCATE( aperm(n), bperm(n), bdist(n) )
  aperm = 0; bperm = 0; bdist=0

  m = Mesh % MaxElementNodes
  ALLOCATE( condition(m), Work(m) )
  condition = -1.0_dp; Work = 0.0_dp

  SolverParams => GetSolverParams()

  DummyDistance = GetLogical( SolverParams,'Dummy Distance Computation',Found ) 


  nb = 0
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF ( .NOT.ASSOCIATED(Element) ) CYCLE

     BodyForce => GetBodyForce()
     IF(ASSOCIATED(BodyForce)) THEN
       nd = GetElementNOFNodes(Element)
       Work(1:nd) = ListGetReal(BodyForce, TRIM(Solver % Variable % Name), &
           nd, Element % NodeIndexes, Found )
       IF ( Found) THEN
          condition(1:nd) = GetReal(BodyForce, TRIM(Solver % Variable % Name) // " Condition", Found )        
          DO i=1,nd             
             j = Element % NodeIndexes(i)

             IF (Found .AND. condition(i) < 0.0) CYCLE

             IF ( bperm(j) == 0 ) THEN
                nb = nb + 1
                aperm(nb) = j
                bperm(j)  = nb
                bdist(j) = work(i)
             END IF
          END DO
       END IF
     END IF
  END DO


  DO t=1,Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)

    IF ( .NOT. ActiveBoundaryElement(Element) ) CYCLE
    BC => GetBC(Element)
    nd = GetElementNOFNodes(Element)

    work(1:nd) = GetReal(BC, Solver % Variable % Name, Found )
    IF ( Found .OR. GetLogical( BC, 'Noslip Wall BC', Found1 ) ) THEN
      DO i=1,nd
        j = Element % NodeIndexes(i)
        IF ( bperm(j) == 0 ) THEN
          nb = nb + 1
          aperm(nb) = j
          bperm(j)  = nb
          bdist(j) = work(i)
        END IF
      END DO
    END IF
  END DO

  
  IF( ParEnv % PEs == 1 ) THEN
    IF( nb == 0 ) THEN
      CALL Warn('DistanceSolver1','No known distances given for the distance solver!')
      RETURN
    ELSE 
      WRITE( Message,'(A,I0)') 'Number of fixed nodes on bulk: ',nb
      CALL Info('DistanceSolver1',Message,Level=8)
    END IF
  END IF


  IF ( ParEnv % PEs > 1 ) THEN
    comm = Solver % Matrix % Comm
    gPerm => Mesh % ParallelInfo % GlobalDOFs

    maxnode = NINT(ParallelReduction(MAXVAL(gPerm(aperm(1:nb)))*1._dp,2))
    ALLOCATE( ibuf(maxnode),cperm(maxnode) );

    cperm=0
    cperm(gperm(aperm(1:nb))) = 1
    CALL MPI_ALLREDUCE(cperm, ibuf, maxnode, MPI_INTEGER, MPI_SUM, comm, i)

    nnb = 0; cperm=0
    DO i=1,maxnode
     IF ( ibuf(i)/=0 ) THEN
       nnb=nnb+1
       cperm(i)=nnb
       ibuf(nnb)=ibuf(i)
     END IF
    END DO

    ALLOCATE( xp(nnb), yp(nnb), zp(nnb), bd(nnb), buf(nnb) )

    buf=0._dp
    buf(cperm(gperm(aperm(1:nb)))) = Mesh % Nodes % x(aperm(1:nb))
    CALL MPI_ALLREDUCE(buf, xp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    xp(1:nnb) = xp(1:nnb)/ibuf(1:nnb)

    buf=0._dp
    buf(cperm(gperm(aperm(1:nb)))) = Mesh % Nodes % y(aperm(1:nb))
    CALL MPI_ALLREDUCE(buf, yp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    yp(1:nnb) = yp(1:nnb)/ibuf(1:nnb)

    buf=0._dp
    buf(cperm(gperm(aperm(1:nb)))) = Mesh % Nodes % z(aperm(1:nb))
    CALL MPI_ALLREDUCE(buf, zp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    zp(1:nnb) = zp(1:nnb)/ibuf(1:nnb)

    buf=0._dp
    buf(cperm(gperm(aperm(1:nb)))) = bdist(aperm(1:nb))
    CALL MPI_ALLREDUCE(buf, bd, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    bd(1:nnb) = bd(1:nnb)/ibuf(1:nnb)

    nb = nnb
    DEALLOCATE(ibuf,buf,cperm)
  ELSE
    ALLOCATE( xp(nb),yp(nb),zp(nb),bd(nb) )

    DO i=1,nb
      xp(i) = Solver % Mesh % Nodes % x(aperm(i))
      yp(i) = Solver % Mesh % Nodes % y(aperm(i))
      zp(i) = Solver % Mesh % Nodes % z(aperm(i))
      bd(i) = bdist(aperm(i))
    END DO
  END IF

  distance => Solver % Variable % Values

  IF( DummyDistance ) THEN
    CALL distcomp0()
  ELSE
    CALL distcomp()
  END IF

  IF( ParEnv % PEs == 1 ) THEN
    WRITE( Message,'(A,ES12.3)') 'Maximum distance from given nodes: ',MAXVAL(distance)
    CALL Info('DistanceSolver1',Message,Level=8)
  END IF


  DEALLOCATE(xp,yp,zp,aperm,bperm,condition)
  Solver % Variable % Norm = SQRT(SUM(distance**2))

  CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, Solver % Variable % Name )	

  CALL CheckTimer('DistanceSolver1',Delete=.TRUE.)
  CALL Info('DistanceSolver1','All done')

CONTAINS

  SUBROUTINE distcomp0

    INTEGER :: i,j,k,cnt
    REAL(KIND=dp) :: xl,yl,zl

    distance = HUGE(1._dp)
    DO i=1,n
      k = Solver % Variable % Perm(i)
      IF ( k <= 0 ) CYCLE
 
      IF ( bperm(i) /= 0 ) THEN
          distance(k) = 0._dp
      ELSE
        xl = Mesh % Nodes % x(i)
        yl = Mesh % Nodes % y(i)
        zl = Mesh % Nodes % z(i)
        DO j=1,nb
          dist = (xp(j)-xl)**2 + (yp(j)-yl)**2 + (zp(j)-zl)**2
          distance(k) = MIN(distance(k), dist)
        END DO
      END IF
    END DO
  END SUBROUTINE distcomp0

  SUBROUTINE distcomp

    REAL(KIND=dp), ALLOCATABLE :: xxp(:),yyp(:),zzp(:),d(:),dd(:), bbd(:)
    REAL(KIND=dp) :: xl,yl,zl,dl,c(3),ldist
    INTEGER :: i,j,k,l,nnb,cnt,minl
    INTEGER, ALLOCATABLE :: dperm(:), near(:)

    nnb = n+nb
    ALLOCATE( xxp(nnb), yyp(nnb), zzp(nnb), d(nnb), dperm(nnb), &
              near(nnb), dd(nb), bbd(nnb) )

    CALL RANDOM_NUMBER(c)

    xxp(1:n) = Mesh % Nodes % x
    yyp(1:n) = Mesh % Nodes % y
    zzp(1:n) = Mesh % Nodes % z
    bbd(1:n) = bdist(1:n)

    xxp(n+1:nnb) = xp
    yyp(n+1:nnb) = yp
    zzp(n+1:nnb) = zp
    bbd(n+1:nnb) = bd
 
    xl = MAXVAL(xxp)-MINVAL(xxp)
    yl = MAXVAL(yyp)-MINVAL(yyp)
    zl = MAXVAL(zzp)-MINVAL(zzp)
    xxp = xxp + xl*(2*c(1)-1)
    yyp = yyp + yl*(2*c(2)-1)
    zzp = zzp + zl*(2*c(3)-1)

    d = SQRT(xxp**2 + yyp**2 + zzp**2)
    dperm = (/ (i,i=1,nnb) /)

    CALL SortR(nnb,dperm,d)

    j = 0
    DO i=1,nnb
      near(i)=j
      IF ( dperm(i)>n ) THEN
        j=j+1
        xp(j)=xxp(dperm(i))
        yp(j)=yyp(dperm(i))
        zp(j)=zzp(dperm(i))
        dd(j)=d(i)
        bd(j)=bbd(dperm(i))
      END IF
    END DO

    DO i=1,nnb
      j=dperm(i)
      IF ( j>n ) CYCLE

      k = Solver % Variable % Perm(j)
      IF ( k <= 0 ) CYCLE

      IF ( bperm(j)/=0 ) THEN
        distance(k) = bdist(j)
      ELSE
        dl = d(i)
        xl = xxp(j)
        yl = yyp(j)
        zl = zzp(j)
        dist = HUGE(1._dp)
        DO l=near(i)+1,nb
          IF( (dd(l)-dl)**2>dist ) EXIT
          ldist = (xp(l)-xl)**2+(yp(l)-yl)**2+(zp(l)-zl)**2
          IF ( dist>ldist) THEN
            minl = l; dist=ldist
          END IF
        END DO

        DO l=near(i),1,-1
          IF( (dd(l)-dl)**2>dist ) EXIT
          ldist = (xp(l)-xl)**2+(yp(l)-yl)**2+(zp(l)-zl)**2
          IF ( dist>ldist) THEN
            minl = l; dist=ldist
          END IF
        END DO

        distance(k) = SQRT(dist)+bd(minl)

      END  IF
    END DO
    DEALLOCATE( xxp, yyp, zzp, d, dperm, near, dd, bbd )
  END SUBROUTINE distcomp
!------------------------------------------------------------------------------
END SUBROUTINE DistanceSolver1
!------------------------------------------------------------------------------


