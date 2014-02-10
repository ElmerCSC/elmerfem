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

!> \ingroup ElmerLib
!> \{

!> Module computing local matrices for stress computation in general coordinates.


MODULE StressGeneral

!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE Integration
  USE ElementDescription

  IMPLICIT NONE

!------------------------------------------------------------------------------
  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE StressGeneralCompose( MassMatrix, StiffMatrix,ForceVector,      &
      LoadVector,NodalYoung, NodalPoisson,NodalDensity,PlaneStress,Isotropic, &
           NodalHeatExpansion, NodalTemperature,Element,n,Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: MassMatrix(:,:), StiffMatrix(:,:),NodalHeatExpansion(:,:,:)
     REAL(KIND=dp) :: NodalTemperature(:),NodalDensity(:),LoadVector(:,:),NodalYoung(:,:,:)
     REAL(KIND=dp), DIMENSION(:) :: ForceVector,NodalPoisson

     LOGICAL :: PlaneStress, Isotropic(2)

     TYPE(Element_t) :: Element
     TYPE(Nodes_t) :: Nodes

     INTEGER :: n
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: Force(3),NodalLame1(n),NodalLame2(n),Lame1,Lame2

     REAL(KIND=dp) :: Load(3),Temperature,X,Y,Z,Density
     REAL(KIND=dp), DIMENSION(3,3) :: A,MT,HeatExpansion

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     INTEGER :: i,j,k,l,m,p,q,t,dim

     REAL(KIND=dp) :: s,u,v,w
  
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat,CylindricSymmetry
!------------------------------------------------------------------------------

     CylindricSymmetry = ( CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                           CurrentCoordinateSystem() == AxisSymmetric )

     IF ( CylindricSymmetry ) THEN
       dim = 3
     ELSE
       dim = CoordinateSystemDimension()
     END IF

     IF ( PlaneStress ) THEN
       NodalLame1(1:n) = NodalYoung(1,1,1:n) * NodalPoisson(1:n) /  &
            ( (1.0d0 - NodalPoisson(1:n)**2) )
     ELSE
       NodalLame1(1:n) = NodalYoung(1,1,1:n) * NodalPoisson(1:n) /  &
          (  (1.0d0 + NodalPoisson(1:n)) * (1.0d0 - 2.0d0*NodalPoisson(1:n)) )
     END IF

     NodalLame2(1:n) = NodalYoung(1,1,1:n)  / ( 2* (1.0d0 + NodalPoisson(1:n)) )


     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
!    
!    Integration stuff
!    
     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
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
!
!      Coordinatesystem dependent info
!
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        X = SUM( nodes % x(1:n)*Basis )
        Y = SUM( nodes % y(1:n)*Basis )
        Z = SUM( nodes % z(1:n)*Basis )
      END IF

      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
!  
      s = SqrtMetric * SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!  
!     Force at integration point
!   
      Force = 0.0D0
      DO i=1,dim
        Force(i) = SUM( LoadVector(i,1:n)*Basis )
      END DO
#if 0
!
!     NOTE: convert unit base vector to contravariant base...
!
      DO i=1,dim
        Force(i) = Force(i) / SQRT(Metric(i,i))
      END DO
#endif
!
!     Lame parameters at the integration point
!
      Lame1 = SUM( NodalLame1(1:n)*Basis )
      Lame2 = SUM( NodalLame2(1:n)*Basis )

      Density = SUM( NodalDensity(1:n)*Basis )
!
!     Temperature at the integration point
!
      Temperature = SUM( NodalTemperature(1:n)*Basis )
!
!     Heat expansion tensor values at the integration point
!
      DO i=1,3
        DO j=1,3
          HeatExpansion(i,j) = SUM( NodalHeatExpansion(i,j,:)*Basis ) * &
               SQRT(Metric(i,i)) * SQRT(Metric(j,j))
        END DO
      END DO
!
!    Loop over basis functions (of both unknowns and weights)
!
     DO p=1,N
     DO q=1,N
       MT = 0.0

       DO i=1,dim
         MT(i,i) = Density * Basis(p) * Basis(q)
       END DO

       A = 0.0
!
!      Stiffness Matrix:
! -----------------------------------
!
       DO i=1,dim
         DO j = 1,dim
           DO k = 1,dim
             A(i,k) = A(i,k) + Lame1 * Metric(i,j) * dBasisdx(q,k) * dBasisdx(p,j)

             A(i,i) = A(i,i) + Lame2 * Metric(j,k) * dBasisdx(q,k) * dBasisdx(p,j)

             A(i,j) = A(i,j) + Lame2 * Metric(i,k) * dBasisdx(q,k) * dBasisdx(p,j)

             IF ( CurrentCoordinateSystem() /= Cartesian .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
                DO l=1,dim
                   A(i,l) = A(i,l)  +  Lame1 * Metric(i,j) * Symb(l,k,k) * Basis(q) * dBasisdx(p,j)

                   A(l,k) = A(l,k)  -  Lame1 * Metric(i,j) * dBasisdx(q,k) * Symb(i,j,l) * Basis(p)

                   A(i,l) = A(i,l)  +  Lame2 * Metric(j,k) * Symb(l,k,i) * Basis(q) * dBasisdx(p,j)

                  A(l,i) = A(l,i)  -  Lame2 * Metric(j,k) * dBasisdx(q,k) * Symb(i,j,l) * Basis(p)

                   A(i,l) = A(i,l)  +  Lame2 * Metric(i,k) * Symb(l,k,j) * Basis(q) * dBasisdx(p,j)

                   A(l,j) = A(l,j)  -  Lame2 * Metric(i,k) * dBasisdx(q,k) * Symb(i,j,l) * Basis(p)
                   DO m=1,dim
                     A(l,m) = A(l,m) - Lame1 * Metric(i,j) * Symb(m,k,k) * Basis(q) * Symb(i,j,l) * Basis(p)

                     A(l,m) = A(l,m) - Lame2 * Metric(j,k) * Symb(m,k,i) * Basis(q) * Symb(i,j,l) * Basis(p)

                     A(l,m) = A(l,m) - Lame2 * Metric(i,k) * Symb(m,k,j) * Basis(q) * Symb(i,j,l) * Basis(p)
                   END DO
                END DO
             END IF

           END DO
         END DO
       END DO

!
! Add nodal matrix to element matrix
!
       IF ( CurrentCoordinateSystem() == AxisSymmetric ) THEN
         A(1,1) = A(1,1) + (Lame1 + 2*Lame2) * Metric(3,3) * Basis(q) * Basis(p)
         A(1,1) = A(1,1) + Lame1 * Symb(1,3,3) * Basis(q) * dBasisdx(p,1)
         A(2,1) = A(2,1) + Lame1 * Symb(1,3,3) * Basis(q) * dBasisdx(p,2)
         A(1,1) = A(1,1) - Lame1 * Metric(3,3) * dBasisdx(q,1) * &
                       Symb(3,3,1) * Basis(p)
         A(1,2) = A(1,2) - Lame1 * Metric(3,3) * dBasisdx(q,2) * &
                       Symb(3,3,1) * Basis(p)

         DO i=1,2
           DO j=1,2
             StiffMatrix( 2*(p-1)+i,2*(q-1)+j ) = &
                 StiffMatrix( 2*(p-1)+i,2*(q-1)+j ) + s * A(i,j)

             MassMatrix( 2*(p-1)+i,2*(q-1)+j ) = &
                 MassMatrix( 2*(p-1)+i,2*(q-1)+j ) + s * MT(i,j)
           END DO
         END DO
       ELSE
         DO i=1,dim
           DO j=1,dim
             StiffMatrix( dim*(p-1)+i,dim*(q-1)+j ) = &
                 StiffMatrix( dim*(p-1)+i,dim*(q-1)+j ) + s*A(i,j)

             MassMatrix( dim*(p-1)+i,dim*(q-1)+j ) = &
                 MassMatrix( dim*(p-1)+i,dim*(q-1)+j ) + s*MT(i,j)
           END DO
         END DO
       END IF

     END DO
     END DO
!
! The righthand side...
!
     DO p=1,N
       Load = 0.0D0
  
       DO i=1,dim
          Load(i) = Load(i) + Force(i) * Basis(p)
       END DO

       DO i=1,dim
         DO j=1,dim
           Load(i) = Load(i) + 2 * Lame2 * HeatExpansion(i,j) * &
                    Temperature * dBasisdx(p,j)
           DO k=1,3
             Load(i) = Load(i) + Metric(i,j) * Lame1 * HeatExpansion(k,k) * &
                    Temperature * dBasisdx(p,j)
           END DO
         END DO
       END DO

       IF ( CurrentCoordinateSystem() == AxisSymmetric ) THEN
         DO i=1,2
           ForceVector(2*(p-1)+i) = ForceVector(2*(p-1)+i) + s*Load(i)
         END DO
       ELSE
         DO i=1,dim
           ForceVector(dim*(p-1)+i) = ForceVector(dim*(p-1)+i) + s*Load(i)
         END DO
       END IF
     END DO

   END DO 
!------------------------------------------------------------------------------
 END SUBROUTINE StressGeneralCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE StressGeneralBoundary( BoundaryMatrix,BoundaryVector,LoadVector, &
                         NodalAlpha,NodalBeta,Element,n,Nodes )
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:)
   REAL(KIND=dp) :: LoadVector(:,:),NodalAlpha(:),NodalBeta(:)
   TYPE(Element_t),POINTER  :: Element
   TYPE(Nodes_t)    :: Nodes

   INTEGER :: n
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n)
   REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

   REAL(KIND=dp) :: u,v,w,s
   REAL(KIND=dp) :: Force(3),Alpha(3),Beta,Normal(3)
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   INTEGER :: i,j,t,q,p,dim,N_Integ

   REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),x,y,z

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

   n   = Element % Type % NumberOfNodes
   dim = Element % Type % Dimension + 1

   BoundaryVector = 0.0D0
   BoundaryMatrix = 0.0D0
!
!  Integration stuff
!
   IntegStuff = GaussPoints( element )
   U_Integ => IntegStuff % u
   V_Integ => IntegStuff % v
   W_Integ => IntegStuff % w
   S_Integ => IntegStuff % s
   N_Integ =  IntegStuff % n
!
!  Now we start integrating
!
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
!
!    Coordinatesystem dependent info
!
     IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
       x = SUM( nodes % x(1:n)*Basis )
       y = SUM( nodes % y(1:n)*Basis )
       z = SUM( nodes % z(1:n)*Basis )
     END IF

     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
!
     s = SqrtMetric * SqrtElementMetric * S_Integ(t)

     Force = 0.0D0
     DO i=1,dim
       Force(i) = SUM( LoadVector(i,1:n)*Basis )
     END DO
!
!------------------------------------------------------------------------------
!    Add to load: given force in normal direction
!------------------------------------------------------------------------------
!
     Normal = NormalVector( Element,Nodes,u,v,.TRUE. )

     Beta  = SUM( NodalBeta(1:n)*Basis )
     DO i=1,dim
       DO j=1,dim
         Force(i) = Force(i) + Beta*Metric(i,j)*Normal(j)
       END DO
     END DO
!------------------------------------------------------------------------------
!

     Alpha(1:3) = SUM( NodalAlpha(1:n) * Basis(1:n) ) * Normal
     DO p=1,N
       DO q=1,N
         DO i=1,dim
           BoundaryMatrix((p-1)*dim+i,(q-1)*dim+i) =  &
               BoundaryMatrix((p-1)*dim+i,(q-1)*dim+i) + &
                  s * Alpha(i) * Basis(q) * Basis(p)
         END DO
       END DO
     END DO

     DO q=1,N
       DO i=1,dim
         BoundaryVector((q-1)*dim+i) = BoundaryVector((q-1)*dim+i) + &
                      s * Basis(q) * Force(i)
       END DO
     END DO

   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE StressGeneralBoundary
!------------------------------------------------------------------------------

END MODULE StressGeneral
