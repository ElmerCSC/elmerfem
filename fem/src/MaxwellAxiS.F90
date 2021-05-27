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
! *  Authors: Juha Ruokolainen, Jussi Heikonen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/



!-----------------------------------------------------------------------------
!>  Module computing MHD Maxwell equations (or the induction equation) local
!>  matrices (cartesian coordinates).
!------------------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{

MODULE MaxwellAxiS

  USE Types
  USE Integration
  USE ElementDescription

  IMPLICIT NONE

  CONTAINS


!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for the MHD Maxwell equation
!------------------------------------------------------------------------------
   SUBROUTINE MaxwellAxiSCompose  (                                     &
       MassMatrix,StiffMatrix,ForceVector,LoadVector,NodalConductivity, &
                   Mx,My,Mz,Ux,Uy,Uz,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT:
!
!  REAL(KIND=dp) :: NodalConductivity(:)
!     INPUT: Nodal values of electric conductivity (times the magnetic
!            permeability)
!
!  REAL(KIND=dp) :: Mx(:),My(:),Mz(:)
!     INPUT: Nodal values of applied magnetic field components
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp),TARGET :: MassMatrix(:,:),StiffMatrix(:,:),ForceVector(:)
     REAL(KIND=dp), DIMENSION(:) :: Ux,Uy,Uz,Mx,My,Mz
     REAL(KIND=dp) :: NodalConductivity(:),LoadVector(:,:)

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
     REAL(KIND=dp) :: SqrtElementMetric

     REAL(KIND=dp) :: Velo(3),dVelodx(3,3),Force(3),Metric(3,3),Symb(3,3,3)
     REAL(KIND=dp) :: MField(3),dMFielddx(3,3)

     REAL(KIND=dp), POINTER :: A(:,:),M(:,:),Load(:)
     REAL(KIND=dp) :: SU(n,4,4),SW(n,4,4)

     REAL(KIND=dp) :: Lambda=1.0,Re,Tau,Delta
     REAL(KIND=dp) :: VNorm,hK,mK,Conductivity,dConductivitydx(3)

     INTEGER :: i,j,k,c,p,q,t,DIM

     REAL(KIND=dp) :: s,u,v,w,r
  
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat
!------------------------------------------------------------------------------

     DIM = 3

     ForceVector = 0.0D0
     MassMatrix  = 0.0D0
     StiffMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n
 
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx )

      r = SUM( basis*nodes%x(1:n) )
      s = r * SqrtElementMetric * S_Integ(t)

!------------------------------------------------------------------------------
!     Applied magnetic field
!------------------------------------------------------------------------------
      MField    = 0.0D0
      MField(1) = SUM( Mx(1:n)*Basis )
      MField(2) = SUM( My(1:n)*Basis )
      MField(3) = SUM( Mz(1:n)*Basis )

      dMFielddx = 0.0D0
      DO i=1,3
        dMFielddx(1,i) = SUM( Mx(1:n)*dBasisdx(1:n,i) )
        dMFielddx(2,i) = SUM( My(1:n)*dBasisdx(1:n,i) )
        dMFielddx(3,i) = SUM( Mz(1:n)*dBasisdx(1:n,i) )
      END DO

!------------------------------------------------------------------------------
!     Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
      Velo = 0.0d0
      Velo(1) = SUM( Ux(1:n)*Basis )
      Velo(2) = SUM( Uy(1:n)*Basis )
      Velo(3) = SUM( Uz(1:n)*Basis )

      dVelodx = 0.0d0
      DO i=1,3
        dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
        dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
        dVelodx(3,i) = r * SUM( Uz(1:n)*dBasisdx(1:n,i) )
      END DO   
      dVelodx(3,1) = Velo(3) + dVelodx(3,1)

      Velo(3) = r * Velo(3)
   
!------------------------------------------------------------------------------
!     Force at integration point
!------------------------------------------------------------------------------
      Force = 0.0D0
      DO i=1,DIM
        Force(i) = SUM( LoadVector(i,1:n)*Basis )
      END DO
!------------------------------------------------------------------------------
!     Effective conductivity
!------------------------------------------------------------------------------
      Conductivity = SUM( NodalConductivity(1:n)*Basis )
!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
     DO p=1,N
     DO q=1,N
!------------------------------------------------------------------------------
!      The MHD Maxwell equations
!------------------------------------------------------------------------------
       i = DIM*(p-1)
       j = DIM*(q-1)
       M => MassMatrix ( i+1:i+DIM,j+1:j+DIM )
       A => StiffMatrix( i+1:i+DIM,j+1:j+DIM )
!------------------------------------------------------------------------------
!      Mass matrix:
!------------------------------------------------------------------------------
       DO i=1,DIM
         M(i,i) = M(i,i) + s * Basis(q) * Basis(p)
       END DO

!------------------------------------------------------------------------------
!      Stiffness matrix:
!------------------------------
!      Diffusive terms
!------------------------------------------------------------------------------

       DO i=1,DIM
         A(i,i) = A(i,i) + ( s * ( dBasisdx(q,1) * dBasisdx(p,1) &
          + dBasisdx(q,2) * dBasisdx(p,2) ) ) / Conductivity
         IF (i /= 2) THEN
           A(i,i) = A(i,i) + ( s * Basis(q) * Basis(p)/ (r**2) ) / Conductivity
         END IF
       END DO
!------------------------------------------------------------------------------
!    The curl(u x B) terms when nabla . u = nabla . B = 0
!------------------------------------------------------------------------------
#ifdef FULL_INDUCTION     
      DO i=1,DIM
        A(i,1) = A(i,1) - s * Basis(q) * dVelodx(i,1) * Basis(p)
        A(i,2) = A(i,2) - s * Basis(q) * dVelodx(i,2) * Basis(p)
        A(i,i) = A(i,i) + s * ( Velo(1) * dBasisdx(q,1) &
            + Velo(2) * dBasisdx(q,2) ) * Basis(p)
      END DO

      A(3,3) = A(3,3) - s * Basis(q) * Velo(1) * Basis(p) / r
      A(3,1) = A(3,1) + s * Basis(q) * Velo(3) * Basis(p) / r
#endif
     END DO
     END DO

!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------
     DO p=1,n
       Load => ForceVector( DIM*(p-1)+1 : DIM*(p-1)+DIM )

       DO i=1,DIM
         Load(i) = Load(i) + s * Force(i) * Basis(p)
       END DO

       ! nabla ^ 2 B_0 / sigma mu  goes here 

!------------------------------------------------------------------------------
!    The curl(u x B) terms
!------------------------------------------------------------------------------

       DO i=1,DIM
         Load(i) = Load(i) + s * ( MField(1) * dVelodx(i,1) &
             + MField(2) * dVelodx(i,2) &
             - Velo(1) * dMFielddx(i,1) &
             - Velo(2) * dMFielddx(i,2) ) * Basis(p)
       END DO

       Load(3) = Load(3) + s * ( MField(3) * Velo(1) - & 
           MField(1) * Velo(3) ) * Basis(p) / r 

     END DO
   END DO

 END SUBROUTINE MaxwellAxiSCompose
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for the MHD Maxwell equation
!------------------------------------------------------------------------------
 SUBROUTINE MaxwellAxiSBoundary( BoundaryMatrix,BoundaryVector,LoadVector, &
                     NodalAlpha,NodalBeta,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:,:)
!     INPUT: Nodal values force in coordinate directions
!
!  REAL(KIND=dp) :: NodalAlpha(:,:)
!     INPUT: Nodal values of force in normal direction
!
!  REAL(KIND=dp) :: NodalBeta(:,:)
!     INPUT: Nodal values of something which will be taken derivative in
!            tangential direction and added to force...
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of boundary element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

   IMPLICIT NONE

   REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:),LoadVector(:,:), &
                             NodalAlpha(:),NodalBeta(:)

   INTEGER :: n

   TYPE(Element_t),POINTER  :: Element
   TYPE(Nodes_t)    :: Nodes

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
   REAL(KIND=dp) :: SqrtElementMetric

   REAL(KIND=dp) :: u,v,w,s
   REAL(KIND=dp) :: Force(3),Alpha
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   INTEGER :: i,t,q,p,c,DIM,N_Integ

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------
   DIM = 3
   BoundaryVector = 0.0D0
   BoundaryMatrix = 0.0D0
!
!------------------------------------------------------------------------------
!  Integration stuff
!------------------------------------------------------------------------------
   IntegStuff = GaussPoints( element )
   U_Integ => IntegStuff % u
   V_Integ => IntegStuff % v
   W_Integ => IntegStuff % w
   S_Integ => IntegStuff % s
   N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!  Now we start integrating
!------------------------------------------------------------------------------
   DO t=1,N_Integ

     u = U_Integ(t)
     v = V_Integ(t)
     w = W_Integ(t)
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
     stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx )

     s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!    Add to load: tangential derivative of something
!------------------------------------------------------------------------------
     DO i=1,DIM
       Force(i) = SUM( NodalBeta(1:n)*dBasisdx(1:n,i) )
     END DO
!------------------------------------------------------------------------------
!    Add to load: given force in normal direction
!------------------------------------------------------------------------------
     Alpha = SUM( NodalAlpha(1:n)*Basis )
     Force = Force + Alpha * NormalVector( element,nodes,u,v,.TRUE. )
!------------------------------------------------------------------------------
!    Add to load: given force in coordinate directions
!------------------------------------------------------------------------------
     DO i=1,DIM
       Force(i) = Force(i) + SUM( LoadVector(i,1:n)*Basis )
     END DO

!    DO p=1,N
!      DO q=1,N
!        DO i=1,dim
!          BoundaryMatrix((p-1)*dim+i,(q-1)*dim+i) =  &
!           BoundaryMatrix((p-1)*dim+i,(q-1)*dim+i) + &
!              s * Gamma(i) * Basis(q) * Basis(p)
!        END DO
!      END DO
!    END DO

     DO q=1,N
       DO i=1,DIM
         BoundaryVector((q-1)*DIM+i) = &
             BoundaryVector((q-1)*DIM+i) + s * Basis(q) * Force(i)
       END DO
     END DO

   END DO

 END SUBROUTINE MaxwellAxiSBoundary
!------------------------------------------------------------------------------

END MODULE MaxwellAxiS

!> \}