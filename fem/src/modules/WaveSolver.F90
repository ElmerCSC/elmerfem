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
! *  Solvers: WaveSolver 
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: ~2013
! *  Modified by: Peter Råback
! *  Modification date: 16 Feb 2015
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Solves the transient wave equation using nodal basis functions and Galerkin 
!> or bubble stabilezed formulation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WaveSolver( Model,Solver,dt,TransientSimulation )
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
  REAL(KIND=dp) :: Norm, Pave
  INTEGER :: n, nb, nd, t, active
  LOGICAL :: Found
!------------------------------------------------------------------------------
  
  CALL DefaultStart()
  
   !System assembly:
   !----------------
   CALL DefaultInitialize()
   Active = GetNOFActive()
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)
      CALL LocalMatrix(  Element, n, nd+nb )
   END DO

   CALL DefaultFinishBulkAssembly()

   Active = GetNOFBoundaryElements()
   DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element)
        nb = GetElementNOFBDOFs(Element)
        CALL LocalMatrixBC(  Element, n, nd+nb )
      END IF
   END DO

   CALL DefaultFinishBoundaryAssembly()

   CALL DefaultFinishAssembly()
   CALL DefaultDirichletBCs()


   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()

   IF( GetLogical( Solver % Values,'Set Average To Zero',Found ) ) THEN
     Pave = SUM( Solver % Variable % Values) / &
         SIZE( Solver % Variable % Values ) 
     Solver % Variable % Values = Solver % Variable % Values - 0.5 * Pave
   END IF

   CALL DefaultFinish()
   

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
!    REAL(KIND=dp) :: diff_coeff(n), conv_coeff(n),react_coeff(n), &
!                     time_coeff(n), D,C,R, rho,Velo(3,n),a(3), Weight
    REAL(KIND=dp) :: speed(n), density(n), damping(n), reaction(n), &
        vel, rho, att, react, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n), DAMP(nd,nd)
    LOGICAL :: Stat,Found,ParametersSet = .FALSE.
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BodyForce, Material

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes, ParametersSet, LoadAtIp, vel, rho, att, react
!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    DAMP = 0._dp
    LOAD = 0._dp

    IF(.NOT. ParametersSet ) THEN
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) &
          Load(1:n) = GetReal( BodyForce, 'sound source', Found )

      Material => GetMaterial()
      speed(1:n) = GetReal(Material,'sound speed',Found)
      density(1:n) = GetReal(Material,'density',Found)    
      damping(1:n) = GetReal(Material,'sound damping',Found)
      reaction(1:n) = GetReal(Material,'sound reaction',Found)
    END IF

!    Velo = 0._dp
!   DO i=1,dim
!      Velo(i,1:n)=GetReal(GetMaterial(),'a '//TRIM(I2S(i)),Found)
!    END DO

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

!      rho = SUM(Basis(1:n)*time_coeff(1:n))
!      a = MATMUL(Velo(:,1:n),Basis(1:n))
!      D = SUM(Basis(1:n)*dens_coeff(1:n))
!      C = SUM(Basis(1:n)*conv_coeff(1:n))
!      R = SUM(Basis(1:n)*react_coeff(1:n))

      IF( .NOT. ParametersSet ) THEN
        LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
        rho = SUM(Basis(1:n)*density(1:n))
        vel = SUM(Basis(1:n)*speed(1:n))
        att = SUM(Basis(1:n)*damping(1:n))
        react = SUM(Basis(1:n)*reaction(1:n))
        ParametersSet = .TRUE.
      END IF

      Weight = IP % s(t) * DetJ

      ! (D*grad(u),grad(v))  (diffusion term)
      ! -----------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          (1.0_dp/rho) * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      DO p=1,nd
        DO q=1,nd
          !  (c^2 grad(u), grad(v))
          ! (advection term (C*grad(u),v)) 
          ! -----------------------------------
          !STIFF (p,q) = STIFF(p,q) + Weight * &
          !   C * SUM(a(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! reaction term (R*u,v)
          ! -----------------------------------
          STIFF(p,q) = STIFF(p,q) + Weight * react * Basis(q) * Basis(p)
          

          ! time d^2u/dt^2,v
          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight  * &
              (1.0_dp/ (rho * vel**2) ) * Basis(q) * Basis(p)

          ! damping term - check the scaling
          !---------------------------------
          DAMP(p,q) = DAMP(p,q) + Weight * &
              att * Basis(q) * Basis(p)                        
        END DO
      END DO

      ! dens*c^2*(Q,v)
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF( TransientSimulation) THEN
       CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
    END IF

    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Speed(n), Density(n), rho, vel, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: FORCE(nd), LOAD(n), DAMP(nd,nd), MASS(nd,nd),STIFF(nd,nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    DAMP = 0._dp
    MASS = 0._dp
    LOAD = 0._dp

    IF( .NOT. GetLogical( BC,'Plane Wave BC',Found ) ) RETURN

    Density(1:n) = GetParentMatProp( 'Density', Element, Found )
    Speed(1:n) = GetParentMatProp( 'Sound Speed',Element, Found )

    IF( .NOT. Found ) RETURN

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Plane wave velo at IP 
      ! ---------------------------
      rho = SUM( Basis(1:n) * Density(1:n) )
      vel = SUM( Basis(1:n) * Speed(1:n) )

      DO p=1,nd
        DO q=1,nd
          DAMP(p,q) = DAMP(p,q) + Weight * Basis(q) * Basis(p) / (rho * vel)
        END DO
      END DO

!      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
    END DO

    IF( TransientSimulation) THEN
       CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
    END IF
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = &
         K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE WaveSolver
!------------------------------------------------------------------------------
