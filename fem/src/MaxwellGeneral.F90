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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/

!-------------------------------------------------------------------------------
!>  Module computing MHD Maxwell equations (or the induction equation) local
!>  matrices.
!-------------------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{



MODULE MaxwellGeneral

  USE Integration
  USE ElementDescription

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for the MHD Maxwell equations
!------------------------------------------------------------------------------
   SUBROUTINE MaxwellGeneralCompose  (                                  &
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
!     permeability)
!
!  REAL(KIND=dp) :: Mx(:),My(:),Mz(:)
!     INPUT: Nodal values of applied magnetic field components
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!
!  TYPE(Element_t) :: Element
!     INPUT: Structure describing the element (dimension, nof nodes,
!             interpolation degree, etc, ... )
!
!  INTEGER :: n
!     INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!     INPUT: Element node coordinates
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

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     REAL(KIND=dp) :: Velo(3),dVelodx(3,3),Force(3)
     REAL(KIND=dp) :: MField(3),dMFielddx(3,3)

     REAL(KIND=dp), POINTER :: A(:,:),M(:,:),Load(:)
     REAL(KIND=dp) :: SU(n,4,4),SW(n,4,4),x,y,z

     REAL(KIND=dp) :: Lambda=1.0,Re,Tau,Delta
     REAL(KIND=dp) :: VNorm,hK,mK,Conductivity,dConductivitydx(3)

     INTEGER :: i,j,k,l,h,c,p,q,t,dim

     REAL(KIND=dp) :: s,u,v,w
  
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat
!------------------------------------------------------------------------------
     dim = 3 

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
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        X = SUM( nodes % x(1:n)*Basis )
        Y = SUM( nodes % y(1:n)*Basis )
        Z = SUM( nodes % z(1:n)*Basis )
      END IF

      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )

      s = SqrtMetric * SqrtElementMetric * S_Integ(t)
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
      Velo = 0.0D0
      Velo(1) = SUM( Ux(1:n)*Basis )
      Velo(2) = SUM( Uy(1:n)*Basis )
      Velo(3) = SUM( Uz(1:n)*Basis )

      dVelodx = 0.0D0
      DO i=1,3
        dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
        dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
        dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
      END DO
!------------------------------------------------------------------------------
!     Force at integration point
!------------------------------------------------------------------------------
      Force = 0.0D0
      DO i=1,dim
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
!      The MHD Maxwell equations (the induction equation)
!------------------------------------------------------------------------------
       i = dim*(p-1)
       j = dim*(q-1)
       M => MassMatrix ( i+1:i+dim,j+1:j+dim )
       A => StiffMatrix( i+1:i+dim,j+1:j+dim )
!------------------------------------------------------------------------------
!      Mass matrix:
!------------------------------------------------------------------------------
       DO i=1,dim
         M(i,i) = M(i,i) + s * Basis(q) * Basis(p)
       END DO
!
!------------------------------------------------------------------------------
!      Stiffness matrix:
!------------------------------
!      Diffusive terms
!------------------------------------------------------------------------------
       s = s / Conductivity
       DO i=1,dim
         DO j = 1,dim
           DO k = 1,dim
             A(i,i) = A(i,i) + s * Metric(j,k) * dBasisdx(q,k) * dBasisdx(p,j)
!            A(i,j) = A(i,j) - s * Metric(i,k) * dBasisdx(q,j) * dBasisdx(p,k) ! xxx
             DO l = 1,dim
               A(l,i) = A(l,i) - s * Metric(j,k) * dBasisdx(q,k) * Symb(i,j,l) * Basis(p)
               A(i,l) = A(i,l) + s * Metric(j,k) * Symb(l,k,i) * Basis(q) * dBasisdx(p,j)

               DO h = 1,dim
                 A(h,l) = A(h,l) - s * Metric(j,k) * Basis(q) * Basis(p) * Symb(l,k,i) * Symb(i,j,h)
               END DO

!              A(l,j) = A(l,j) + s * Metric(i,k) * dBasisdx(q,j) * Symb(i,k,l) * Basis(p)
!              A(i,l) = A(i,l) - s * Metric(i,k) * Symb(l,j,j) * Basis(q) * dBasisdx(p,k)
!              DO h = 1,dim
!                A(h,l) = A(h,l) + s * Metric(i,k) * Basis(q) * Basis(p) * Symb(l,j,j) * Symb(i,k,h)
!              END DO
             END DO
           END DO
         END DO
       END DO
       s = s * Conductivity

!

!------------------------------------------------------------------------------
!    The curl(u x B) terms
!------------------------------------------------------------------------------
       DO i=1,dim
         DO j=1,dim
!          A(i,j) = A(i,j) - s * Velo(i) * dBasisdx(q,j) * Basis(p)
!          A(i,i) = A(i,i) + s * Basis(q) * dVelodx(j,j) * Basis(p)
           A(i,j) = A(i,j) - s * Basis(q) * dVelodx(i,j) * Basis(p)
           A(i,i) = A(i,i) + s * Velo(j) * dBasisdx(q,j) * Basis(p)
           DO k=1,dim
!            A(i,k) = A(i,k) - s * Symb(k,j,j) * Velo(i) * Basis(q) * Basis(p)
!            A(i,i) = A(i,i) + s * Symb(k,j,j) * Velo(k) * Basis(q) * Basis(p)
             A(i,j) = A(i,j) - s * Symb(k,j,i) * Velo(k) * Basis(q) * Basis(p)
             A(i,k) = A(i,k) + s * Symb(k,j,i) * Velo(j) * Basis(q) * Basis(p)
           END DO
         END DO
       END DO
     END DO
     END DO

!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------

     DO p=1,N
       Load => ForceVector( dim*(p-1)+1 : dim*(p-1)+dim )

       DO i=1,dim
         Load(i) = Load(i) + s * Force(i) * Basis(p)
#if 0
         s = s / Conductivity
         DO j = 1,dim
           DO k = 1,dim
             Load(i) = Load(i) - s * Metric(j,k) * dMFielddx(i,k) * dBasisdx(p,j)
!            Load(i) = Load(i) + s * Metric(i,k) * dMFielddx(j,j) * dBasisdx(p,k) ! xxx
             DO l = 1,dim
               Load(l) = Load(l) + s * Metric(j,k) * dMFielddx(i,k) * Symb(i,j,l) * Basis(p)
               Load(i) = Load(i) - s * Metric(j,k) * Symb(l,k,i) * MField(l) * dBasisdx(p,j)
               DO h = 1,dim
                 Load(h) = Load(h) + s * Metric(j,k) * MField(l) * Basis(p) * Symb(l,k,i) * Symb(i,j,h)
               END DO

!              Load(l) = Load(l) - s * Metric(i,k) * dMFielddx(j,j) * Symb(i,k,l) * Basis(p)
!              Load(i) = Load(i) + s * Metric(i,k) * Symb(l,j,j) * MField(l) * dBasisdx(p,k)
!              DO h = 1,dim
!                Load(h) = Load(h) - s * Metric(i,k) * MField(l) * Basis(p) * Symb(l,j,j) * Symb(i,k,h)
!              END DO
             END DO
           END DO
         END DO
         s = s * Conductivity
#endif
       END DO
!------------------------------------------------------------------------------
!    The curl(u x B) terms
!------------------------------------------------------------------------------
       DO i=1,dim
         DO j=1,dim
! u\nabla.B
!          Load(i) = Load(i) + s * Velo(i) * dMFielddx(j,j) * Basis(p)
! B\nabla.u
!          Load(i) = Load(i) - s * MField(i) * dVelodx(j,j) * Basis(p)
! B.\nabla u
           Load(i) = Load(i) + s * MField(j) * dVelodx(i,j) * Basis(p)
! u.\nabla B
           Load(i) = Load(i) - s * Velo(j) * dMFielddx(i,j) * Basis(p)
           DO k=1,dim
!            Load(i) = Load(i) + s * Symb(k,j,j) * Velo(i) * MField(k) * Basis(p)
!            Load(i) = Load(i) - s * Symb(k,j,j) * Velo(k) * MField(i) * Basis(p)
             Load(i) = Load(i) + s * Symb(k,j,i) * Velo(k) * MField(j) * Basis(p)
             Load(i) = Load(i) - s * Symb(k,j,i) * Velo(j) * MField(k) * Basis(p)
           END DO
         END DO
       END DO
     END DO
   END DO 

 END SUBROUTINE MaxwellGeneralCompose
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for the MHD
!>  Maxwell equation
!------------------------------------------------------------------------------
 SUBROUTINE MaxwellGeneralBoundary( BoundaryMatrix,BoundaryVector,LoadVector, &
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

   REAL(KIND=dp) :: u,v,w,x,y,z,s
   REAL(KIND=dp) :: Force(3),Alpha
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

   INTEGER :: i,t,q,p,c,dim,N_Integ

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------
   dim = 3
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
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         X = SUM( nodes % x(1:n)*Basis )
         Y = SUM( nodes % y(1:n)*Basis )
         Z = SUM( nodes % z(1:n)*Basis )
       END IF

       CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )

       s = SqrtMetric * SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!    Add to load: tangential derivative of something
!------------------------------------------------------------------------------
     DO i=1,dim
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
     DO i=1,dim
       Force(i) = Force(i) + SUM( LoadVector(i,1:n)*Basis )
     END DO

!    DO p=1,N
!      DO q=1,N
!        DO i=1,dim
!          BoundaryMatrix((p-1)*dim+i,(q-1)*dim+i) =  &
!             BoundaryMatrix((p-1)*dim+i,(q-1)*dim+i) + &
!               s * Gamma(i) * Basis(q) * Basis(p)
!        END DO
!      END DO
!    END DO

     DO q=1,N
       DO i=1,dim
         BoundaryVector((q-1)*dim+i) = &
             BoundaryVector((q-1)*dim+i) + s * Basis(q) * Force(i)
       END DO
     END DO

   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE MaxwellGeneralBoundary
!------------------------------------------------------------------------------

END MODULE MaxwellGeneral

!> \}