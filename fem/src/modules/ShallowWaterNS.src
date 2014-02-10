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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 22.12.2011
! *
! *****************************************************************************/

!-----------------------------------------------------------------------------
!>  Solve the shallow water n-s equations
!-----------------------------------------------------------------------------
SUBROUTINE ShallowWaterNSSolver( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: iter,  n, nb, nd, t, istat, active

  LOGICAL :: Found, NewtonLinearization=.FALSE., MassFluxDiscretization

  TYPE(ValueList_t), POINTER :: Params

  INTEGER :: NewtonIter,  NonlinIter
  REAL(KIND=dp):: NewtonTol, NonlinTol, RelativeChange
!------------------------------------------------------------------------------
  Params => GetSolverParams(Solver)

  NewtonIter = GetInteger(Params,'Nonlinear System Newton After Iterations',Found)
  NewtonTol  = GetConstReal(Params,'Nonlinear System Newton After Tolerance',Found)

  NonlinIter = GetInteger(Params,'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter=1

  NonlinTol = GetConstReal(Params,'Nonlinear System Convergence Tolerance',Found)

  MassFluxDiscretization = GetLogical(Params,'Mass Flux Discretization', Found )

  DO iter=1,NonlinIter
    CALL DefaultInitialize()
    CALL BulkAssembly()
    CALL DefaultFinishBulkAssembly()

    CALL BoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    RelativeChange = Solver % Variable % NonlinChange

    WRITE( Message, * ) 'Result Norm   : ',Norm
    CALL Info( 'ShallowWaterNSSolver', Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ',RelativeChange
    CALL Info( 'ShallowWaterNSSolver', Message, Level=4 )

    IF ( RelativeChange < NewtonTol .OR. iter >= NewtonIter ) &
             NewtonLinearization = .TRUE.
    IF ( RelativeChange < NonlinTol ) EXIT
  END DO

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
    INTEGER :: t,n,nd

!$omp parallel do private(Element,n,nd)
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL LocalMatrix( Element, n, nd )
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(3*nd,3*nd), STIFF(3*nd,3*nd), FORCE(3*nd), LOAD(3,nd)
    REAL(KIND=dp), POINTER :: A(:,:),M(:,:)

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), DetJ, SOL(3,nd),Depth(nd)
    LOGICAL :: Stat,isScalar,Found
    INTEGER :: i,j,k,p,q,t,l,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: s, f, g, diff, Speed, ICz, d, h, U(3), diff_n(n), &
     gradU(2,2), gradF(2), gradH(2), gradd(2), Uwind(2), wcoeff, Uwind_n(2,n)

    TYPE(ValueList_t), POINTER :: Material, BF

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    dim=CoordinateSystemDimension()
    L=dim+1

    CALL GetElementNodes( Nodes,Element )

    Material => GetMaterial(Element)

    BF => GetBodyForce(Element)
    IF (ASSOCIATED(BF)) THEN
    END IF

    CALL GetVectorLocalSolution(SOL,UElement=Element)
    CALL GetScalarLocalSolution(Depth,'Ground',UElement=Element)
    Depth = -Depth

    g = GetCReal(Material,'Graviational Acceleration', Found)

    ICz = GetCReal(Material,'Inverse Chezy Coefficient',Found)

    wcoeff = GetCReal(Material,'Wind-Stress Coefficient',Found)
    Uwind_n(1,1:n) = GetReal(Material,'Wind 1',Found,Element)
    Uwind_n(2,1:n) = GetReal(Material,'Wind 2',Found,Element) 
    f = GetCReal(Material,'Coriolis Coefficient',Found)

    diff_n(1:n) = GetReal(Material,'Viscosity Coefficient',Found,Element)

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    Uwind(1)=0.0_dp
    Uwind(2)=0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      s = IP % s(t) * DetJ

      DO i=1,L
        U(i) = SUM( SOL(i,1:nd) * Basis(1:nd) )
      END DO

      gradU(1:dim,1:dim) = MATMUL(SOL(1:dim,1:nd),dBasisdx(1:nd,1:dim))

      DO i=1,dim
        gradF(i) = SUM(SOL(L,1:nd)*dBasisdx(1:nd,i))
      END DO

      Speed = ICz**2*SQRT(SUM(U(1:dim)**2))

      h = SUM(Depth(1:n)*Basis(1:n))
      DO i=1,dim
        gradh(i)=SUM(Depth(1:n)*dBasisdx(1:n,i))
      END DO

      DO i=1,dim
        Uwind(i) = SUM(Basis(1:n)*Uwind_n(i,1:n))
      END DO

      diff = SUM(Basis(1:n)*Diff_n(1:n))

      DO p=1,nd
        i=L*(p-1)
        DO q=1,nd
          j=L*(q-1)
          M => MASS(i+1:i+L,j+1:j+L)
          A => STIFF(i+1:i+L,j+1:j+L)

          DO j=1,L
            M(j,j) = M(j,j) + s * Basis(q)*Basis(p)
          END DO

          DO j=1,dim
            A(j,j) = A(j,j) + s * diff*SUM(dBasisdx(q,:)*dBasisdx(p,:))
            DO k=1,dim
              A(j,j) = A(j,j) + s * U(k)*dBasisdx(q,k)*Basis(p)
              IF (NewtonLinearization) THEN
                A(j,k) = A(j,k) + s * Basis(q)*gradU(j,k)*Basis(p)
              END IF
            END DO
            A(j,j) = A(j,j) + s * g * Speed * Basis(q) * Basis(p)
            IF (MassFluxDiscretization) THEN
              A(j,L) = A(j,L) + s * g * dBasisdx(q,j) * Basis(p)
            ELSE
              A(j,L) = A(j,L) - s * g * Basis(q) * dBasisdx(p,j)
            END IF

            A(j,j) = A(j,j) + s * wcoeff*Basis(q)/(h+U(L))*Basis(p)

            IF (MassFluxDiscretization) THEN
              A(L,j) = A(L,j) - s * Basis(q)*(h+U(L))*dBasisdx(p,j)
              IF(NewtonLinearization) &
                A(L,L) = A(L,L) - s * U(j)*Basis(q)*dBasisdx(p,j)
            ELSE
              A(L,j) = A(L,j) + s * (dBasisdx(q,j)*(h+U(L))+Basis(q)*(gradh(j)+gradf(j)))*Basis(p)
              IF(NewtonLinearization) &
                A(L,L) = A(L,L) + s * (gradU(j,j)*Basis(q)+U(j)*dBasisdx(q,j))*Basis(p)
            END IF
          END DO
          IF (dim==2) THEN
            A(1,2) = A(1,2) + s * f * Basis(q) * Basis(p)
            A(2,1) = A(2,1) + s * f * Basis(q) * Basis(p)
          END IF
        END DO
        DO j=1,dim
          IF(NewtonLinearization) THEN
            IF (MassFluxDiscretization) THEN
              FORCE(i+L) = FORCE(i+L) - s * U(j)*U(L)*dBasisdx(p,j)
            ELSE
              FORCE(i+L) = FORCE(i+L) + s * (U(L)*gradU(j,j)+U(j)*gradf(j))*Basis(p)
            END IF

            DO k=1,dim
              FORCE(i+j) = FORCE(i+j) + s * U(k)*gradU(j,k) * Basis(p)
            END DO
          END IF
        END DO
        DO j=1,dim
          FORCE(i+j) = FORCE(i+j) + s * wcoeff*Uwind(j)/(h+U(L))*Basis(p)
        END DO
      END DO
    END DO

    j=n+1
    IF (dim==1) j=MAX(j,nd-2)

    DO i=j,nd
     q=L*(i-1)+L
     STIFF(q,:)=0
     STIFF(:,q)=0

     MASS(q,:)=0
     MASS(:,q)=0

     STIFF(q,q)=1
     FORCE(q)=0
    END DO

    IF(Transient) THEN
      CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element,USolver=Solver)
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element,USolver=Solver)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BoundaryAssembly()
!------------------------------------------------------------------------------
    INTEGER :: t,n,nd

!$omp parallel do private(Element,n,nd)
    DO t=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(t)
      IF(.NOT.ActiveBoundaryElement(Element)) CYCLE

      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL BoundaryLocalMatrix( Element, n, nd )
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BoundaryLocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(3*nd,3*nd), STIFF(3*nd,3*nd), FORCE(3*nd)
    REAL(KIND=dp), POINTER :: A(:,:),M(:,:)

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), s, DetJ
    REAL(KIND=dp) :: SOL(3,nd), Udot, U(3), h, Depth(nd), Nrm(3)
    LOGICAL :: Stat,Found,Rout
    INTEGER :: i,j,p,q,t,pind(2),dim,L
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t), SAVE :: Nodes, PNodes
!$omp threadprivate(Nodes,PNodes)
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    L = dim+1

    IF(.NOT.MassFluxDiscretization) RETURN

    BC => GetBC(Element)
    IF(.NOT.ASSOCIATED(BC)) RETURN

    CALL GetElementNodes( Nodes,Element )
    rout = .FALSE.
    IF(dim==1) THEN
     Parent => Element % BoundaryInfo % Left
     CALL GetElementNodes( PNodes,Parent )
     pind =[1,2]
     IF (Element % NodeIndexes(1)/=Parent % NodeIndexes(1)) pind=[2,1]
     rout = PNodes % x(pind(1)) > PNodes % x(pind(2))
    END IF

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    CALL GetVectorLocalSolution(SOL,UElement=Element)
    CALL GetScalarLocalSolution(Depth,'Ground',UElement=Element)
    Depth = -Depth

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      s = IP % s(t) * DetJ

      DO i=1,dim+1
        U(i)=SUM(SOL(i,1:nd)*Basis(1:nd))
      END DO
      IF (dim>1) THEN
        Nrm = NormalVector(Element,Nodes,ip % U(t), ip % V(t), .TRUE.)
        Udot = SUM(Nrm(1:dim)*U(1:dim))
      ELSE IF(rout) THEN
        Udot =  U(1)
      ELSE
        Udot = -U(1)
      END IF

      h = SUM(Basis(1:n)*Depth(1:n))

      DO p=1,nd
        i=L*(p-1)
        DO q=1,nd
          j=L*(q-1)
          A => STIFF(i+1:i+L,j+1:j+L)
          A(L,L) = A(L,L) + s*Udot*Basis(q)*Basis(p)
        END DO
        FORCE(i+L) = FORCE(i+L) - s*Udot*h*Basis(p)
      END DO
    END DO

    IF(Transient) THEN
      CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE ShallowWaterNSSolver
!------------------------------------------------------------------------------
