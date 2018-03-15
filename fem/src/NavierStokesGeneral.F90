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

!------------------------------------------------------------------------------
!>  Module computing Navier-Stokes local matrices in general coordinate system 
!> (i.e. not cartesian, axisymmetric or cylindrically symmetric.
!------------------------------------------------------------------------------

MODULE NavierStokesGeneral

  USE CoordinateSystems
  USE Integration
  USE Differentials
  USE Materialmodels

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for Navier-Stokes-Equations
!>  in general Euclidian coordinate system.
!------------------------------------------------------------------------------
   SUBROUTINE NavierStokesGeneralCompose  (                               &
            MassMatrix,StiffMatrix,ForceVector,LoadVector,NodalViscosity, &
        NodalDensity,Ux,Uy,Uz,MUx,MUy,MUz,Stabilize,NewtonLinearization,Element,n,Nodes )
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
!  REAL(KIND=dp) :: NodalViscosity(:)
!     INPUT: Nodal values for viscosity (i.e. if turbulence model, or
!             power-law viscosity is used, the values vary in space)
!
!  REAL(KIND=dp) :: NodalDensity(:)
!     INPUT: nodal density values
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ?
!
!  LOGICAL :: NewtonLinearization
!      INPUT: Picard or Newton  linearization of the convection term ?
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix,LoadVector
     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Ux,Uy,Uz,MUx,MUy,MUz
     REAL(KIND=dp) :: NodalViscosity(:),NodalDensity(:)
     LOGICAL :: Stabilize,NewtonLinearization

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Velo(3),UVelo(3),dVelodx(3,3),Force(4)

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     REAL(KIND=dp), DIMENSION(4,4) :: A,Mass
     REAL(KIND=dp), DIMENSION(4) :: Load,SU(n,4,4),SW(n,4,4),LrF(3)

     REAL(KIND=dp) :: Lambda=1.0,Re,Tau,Delta,x0,y0
     REAL(KIND=dp) :: VNorm,hK,mK,Viscosity,dViscositydx(3),Density

     INTEGER :: i,j,k,l,m,c,p,q,t,dim,N_Integ

     REAL(KIND=dp) :: s,u,v,w,x,y,z
  
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
 
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat,CylindricSymmetry

     INTEGER :: IMap(3) = (/ 1,2,4 /)
!------------------------------------------------------------------------------

     CylindricSymmetry = (CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                   CurrentCoordinateSystem() == AxisSymmetric)

     IF ( CylindricSymmetry ) THEN
       dim = 3
     ELSE
       dim = CoordinateSystemDimension()
     END IF
     c = dim + 1

     N = element % Type % NumberOfNodes

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ  = IntegStuff % n
!
!------------------------------------------------------------------------------
!    Stabilization parameter mK
!------------------------------------------------------------------------------
     IF ( Stabilize ) THEN
       hK = element % hK
       mK = element % StabilizationMK
     END IF
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
!
    DO t=1,N_Integ
!
!     Integration stuff
!
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx,ddBasisddx,Stabilize )
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
!      density at the integration point
!------------------------------------------------------------------------------
       Density = SUM( NodalDensity(1:n)*Basis )

!------------------------------------------------------------------------------
!     Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
      Velo    = 0.0D0
      Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis )
      Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis )
      IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
        Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis )
      END IF

      IF ( NewtonLinearization ) THEN
        UVelo    = 0.0D0
        UVelo(1) = SUM( Ux(1:n) * Basis )
        UVelo(2) = SUM( Uy(1:n) * Basis )
        IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
          UVelo(3) = SUM( Uz(1:n) * Basis )
        END IF

        dVelodx = 0.0D0
        DO i=1,3
          dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,1) )
          dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,2) )

          IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
            dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,3) )
          END IF
        END DO
      END IF
!  
!------------------------------------------------------------------------------
!     Force at the integration point
!------------------------------------------------------------------------------
      LrF  = LorentzForce( Element,Nodes,u,v,w,n )

      Force = 0.0D0
      Force(1) = SUM( LoadVector(1,1:n)*Basis )
      Force(2) = SUM( LoadVector(2,1:n)*Basis )
      Force(3) = SUM( LoadVector(3,1:n)*Basis )
      IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
        Force(4) = SUM( LoadVector(4,1:n)*Basis )
        Force(1:3) = Force(1:3) + LrF / Density
      ELSE
        Force(1:2) = Force(1:2) + LrF(1:2) / Density
      END IF

#if 0
!
!     NOTE: convert unit base vector, to contravariant base !
!     TODO: this actually valid for orthogonal coordinates...
!
      DO i=1,3
        Force(i) = Force(i) * SQRT(Metric(i,i))
      END DO
#endif
!------------------------------------------------------------------------------
!     Effective viscosity & derivatives at integration point
!------------------------------------------------------------------------------
      Viscosity = SUM( NodalViscosity(1:n)*Basis )
      Viscosity = EffectiveViscosity( Viscosity, Density, Ux, Uy, Uz, &
            Element, Nodes, n, n, u, v, w, LocalIP=t )
!   
!------------------------------------------------------------------------------
!      Stabilization parameters Tau & Delta
!------------------------------------------------------------------------------
     IF ( Stabilize ) THEN
       DO i=1,3
         dViscositydx(i) = SUM( NodalViscosity(1:n)*dBasisdx(1:n,i) )
       END DO

       VNorm = MAX( SQRT( SUM(Velo(1:dim)**2) ), 1.0D-12 )
!
!      VNorm = 0
!      DO i=1,dim
!        VNorm = VNorm + Velo(i)*Velo(i)/Metric(i,i)
!      END DO
!      VNorm = SQRT(VNorm)

       Re = MIN( 1.0D0, Density * mK * hK * VNorm / (4 * Viscosity) )

       Tau = 0.0D0
       IF ( VNorm /= 0.0D0 ) THEN
         Tau = hK * Re / (2 * Density * VNorm)
       END IF

       Delta = Density * Lambda * VNorm * Re * hK
!
!------------------------------------------------------------------------------
!      SU will contain residual of ns-equations, SW will contain the
!      weight function terms
!------------------------------------------------------------------------------
       SU = 0.0D0
       SW = 0.0D0
       DO p=1,N
          DO i=1,dim
            DO j=1,dim
              SU(p,i,c) = SU(p,i,c) + Metric(i,j) * dBasisdx(p,j)

              SU(p,i,i) = SU(p,i,i) + Density * dBasisdx(p,j) * Velo(j)
              IF ( NewtonLinearization ) THEN
                SU(p,i,j) = SU(p,i,j) + Density * dVelodx(i,j) * Basis(p)
              END IF

              DO k =1,dim
                SU(p,i,i) = SU(p,i,i) - Viscosity * Metric(j,k) * ddBasisddx(p,j,k)

                SU(p,i,i) = SU(p,i,i) - dViscositydx(j) * Metric(j,k) * dBasisdx(p,k)

                SU(p,i,j) = SU(p,i,j) - Viscosity * Metric(i,k) * ddBasisddx(p,j,k)

                SU(p,i,j) = SU(p,i,j) - dViscositydx(j) * Metric(i,k) * dBasisdx(p,k)

                IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                  SU(p,i,k) = SU(p,i,k) + Density * Symb(k,j,i) * Basis(p) * Velo(j)
                  IF ( NewtonLinearization ) THEN
                      SU(p,i,j) = SU(p,i,j) + Density * Symb(k,j,i) * Velo(k) * Basis(p)
                  END IF

                  DO l=1,dim
                    SU(p,i,i) = SU(p,i,i) + Viscosity * Metric(j,k) * Symb(j,k,l) * dBasisdx(p,l)

                    SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j,k) * Symb(l,j,i) * dBasisdx(p,k)

                    SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j,k) * Symb(l,k,i) * dBasisdx(p,j)

                    SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j,k) * dSymb(l,j,i,k) * Basis(p)

                    SU(p,i,j) = SU(p,i,j) + Viscosity * Metric(i,k) * Symb(j,k,l) * dBasisdx(p,l)

                    SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i,k) * Symb(l,j,j) * dBasisdx(p,k)

                    SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i,k) * Symb(l,k,j) * dBasisdx(p,j)

                    SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i,k) * dSymb(l,j,j,k) * Basis(p)

                    DO m=1,dim
                      SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j,k) * Symb(m,k,i) * Symb(l,j,m) * Basis(p)

                      SU(p,i,l) = SU(p,i,l) + Viscosity * Metric(j,k) * Symb(j,k,m) * Symb(l,m,i) * Basis(p)

                      SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i,k) * Symb(m,k,j) * Symb(l,j,m) * Basis(p)

                      SU(p,i,l) = SU(p,i,l) + Viscosity * Metric(i,k) * Symb(j,k,m) * Symb(l,m,j) * Basis(p)
                    END DO
                  END DO
                END IF

              END DO
            END DO

            SW(p,i,c) = SW(p,i,c) + Density * dBasisdx(p,i)

            DO j=1,dim

              SW(p,i,i) = SW(p,i,i) + Density * dBasisdx(p,j) * Velo(j)
              DO k =1,dim
                SW(p,i,i) = SW(p,i,i) + Viscosity * Metric(j,k) * ddBasisddx(p,j,k)

                SW(p,i,i) = SW(p,i,i) + dViscositydx(j) * Metric(j,k) * dBasisdx(p,j)

                SW(p,i,j) = SW(p,i,j) + Viscosity * Metric(j,k) * ddBasisddx(p,i,k)

                SW(p,i,j) = SW(p,i,j) + dViscositydx(j) * Metric(j,k) * dBasisdx(p,i)

                IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                  SW(p,i,k) = SW(p,i,k) - Density * Symb(i,j,k) * Basis(p) * Velo(j)
                  DO l=1,dim
                    SW(p,i,i) = SW(p,i,i) - Viscosity * Metric(j,k) * Symb(j,k,l) * dBasisdx(p,l)

                    SW(p,i,j) = SW(p,i,j) - Viscosity * Metric(j,k) * Symb(i,k,l) * dBasisdx(p,l)

                    SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j,k) * Symb(i,j,l) * dBasisdx(p,k)

                    SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j,k) * Symb(i,j,l) * dBasisdx(p,k)

                    SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j,k) * Symb(i,k,l) * dBasisdx(p,j)

                    SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j,k) * Symb(j,k,l) * dBasisdx(p,i)

                    SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j,k) * dSymb(i,j,l,k) * Basis(p)

                    SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j,k) * dSymb(i,j,l,k) * Basis(p)

                    DO m=1,dim
                      SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j,k) * Symb(i,k,m) * Symb(m,j,l) * Basis(p)
 
                      SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j,k) * Symb(j,k,m) * Symb(m,i,l) * Basis(p)

                      SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j,k) * Symb(i,k,m) * Symb(m,j,l) * Basis(p)

                      SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j,k) * Symb(j,k,m) * Symb(m,i,l) * Basis(p)
                    END DO
                  END DO
                END IF

              END DO
            END DO

            IF ( CurrentCoordinateSystem() == AxisSymmetric ) THEN
              SU(p,i,3) = 0.0D0
              SW(p,i,3) = 0.0D0
            END IF

         END DO
       END DO
     END IF

!
!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
!
     DO p=1,N
     DO q=1,N
!
!------------------------------------------------------------------------------
! First plain Navier-Stokes
!------------------------------------------------------------------------------
!
!      Mass Matrix:
! -----------------------------------
!------------------------------------------------------------------------------
!
       Mass = 0.0D0
       DO i=1,dim
         Mass(i,i) = Density * Basis(q) * Basis(p)
       END DO

       A = 0.0
!
!------------------------------------------------------------------------------
!      Stiffness Matrix:
! -----------------------------------
!      Diffusive terms
!------------------------------------------------------------------------------
!
       DO i=1,dim
         DO j = 1,dim
           DO k = 1,dim
             A(i,i) = A(i,i) + Viscosity * Metric(j,k) * dBasisdx(q,k) * dBasisdx(p,j)

             A(i,j) = A(i,j) + Viscosity * Metric(i,k) * dBasisdx(q,k) * dBasisdx(p,j)

             IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
               DO l=1,dim
                 A(i,l) = A(i,l) + Viscosity * Metric(j,k) * Symb(l,k,i) * Basis(q) * dBasisdx(p,j)

                 A(i,l) = A(i,l) + Viscosity * Metric(i,k) * Symb(l,k,j) * Basis(q) * dBasisdx(p,j)

                 A(l,i) = A(l,i) - Viscosity * Metric(j,k) * dBasisdx(q,k) * Symb(i,j,l) * Basis(p)

                 A(l,j) = A(l,j) - Viscosity * Metric(i,k) * dBasisdx(q,k) * Symb(i,j,l) * Basis(p)
                 DO m=1,dim
                   A(l,m) = A(l,m) - Viscosity * Metric(j,k) * Symb(m,k,i) * Basis(q) * Symb(i,j,l) * Basis(p)

                   A(l,m) = A(l,m) - Viscosity * Metric(i,k) * Symb(m,k,j) * Basis(q) * Symb(i,j,l) * Basis(p)
                 END DO
               END DO
             END IF

           END DO
         END DO
       END DO
!
!------------------------------------------------------------------------------
!      Convection terms, Picard linearization
!------------------------------------------------------------------------------
!
       DO i=1,dim
         DO j = 1,dim
           A(i,i) = A(i,i) + Density * dBasisdx(q,j) * Velo(j) * Basis(p)
           IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             DO k=1,dim 
                A(i,k) = A(i,k) + Density * Symb(k,j,i) * Basis(q) * Velo(j) * Basis(p)
             END DO
           END IF
         END DO
       END DO
!
!------------------------------------------------------------------------------
!      Newton linearization
!------------------------------------------------------------------------------
!
       IF ( NewtonLinearization ) THEN
         DO i=1,dim
           DO j=1,dim
             A(i,j) = A(i,j) + Density * dVelodx(i,j) * Basis(q) * Basis(p)
             IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
               DO k=1,dim
                 A(i,j) = A(i,j) + Density * Symb(k,j,i) * Velo(k) * Basis(q) * Basis(p)
               END DO
             END IF
           END DO
         END DO
       END IF
!
!------------------------------------------------------------------------------
!      Pressure terms
!------------------------------------------------------------------------------
       DO i=1,dim
         DO j=1,dim
           A(i,c) = A(i,c) - Metric(i,j) * Basis(q) * dBasisdx(p,j)
           IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             DO k=1,dim
                A(k,c) = A(k,c) + Metric(i,j) * Basis(q) * Symb(i,j,k) * Basis(p)
             END DO
           END IF
         END DO
       END DO
!
!------------------------------------------------------------------------------
!      Continuity equation
!------------------------------------------------------------------------------
!
       DO i=1,dim
         A(c,i) = A(c,i) + Density * dBasisdx(q,i) * Basis(p)
         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
           DO j=1,dim
             A(c,j) = A(c,j) + Density * Symb(j,i,i) * Basis(q) * Basis(p)
           END DO
         END IF
       END DO

!
!------------------------------------------------------------------------------
!      Stabilization...
!------------------------------------------------------------------------------
!
       IF ( Stabilize ) THEN
         DO i=1,dim
           DO j=1,c
             Mass(j,i) = Mass(j,i) + Tau * Density * Basis(q) * SW(p,i,j)
             DO k=1,c
               A(j,k) = A(j,k) + Tau * SU(q,i,k) * SW(p,i,j)
             END DO
           END DO

           DO j=1,dim
             DO k=1,dim
               A(j,i) = A(j,i) + Delta * dBasisdx(q,i) * Metric(j,k) * dBasisdx(p,k)
               IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                 DO l=1,dim
                   A(l,i) = A(l,i) - Delta * dBasisdx(q,i) * Metric(j,k) * Symb(j,k,l) * Basis(p)
                   A(j,l) = A(j,l) + Delta * Symb(l,i,i) * Basis(q) * Metric(j,k) * dBasisdx(p,k)
                   DO m=1,dim
                     A(m,l) = A(m,l) - Delta * Symb(l,i,i) * Basis(q) * Metric(j,k) * Symb(j,k,m) * Basis(p)
                   END DO
                 END DO
               END IF
             END DO
           END DO
         END DO
       END IF
!
!------------------------------------------------------------------------------
! Add nodal matrix to element matrix
!------------------------------------------------------------------------------
!
       IF ( CurrentCoordinateSystem() == AxisSymmetric ) THEN
         DO i=1,3
           DO j=1,3
             StiffMatrix( 3*(p-1)+i, 3*(q-1)+j ) = &
                StiffMatrix( 3*(p-1)+i, 3*(q-1)+j ) + s * A(IMap(i),IMap(j))

             MassMatrix(  3*(p-1)+i, 3*(q-1)+j ) = &
                MassMatrix( 3*(p-1)+i, 3*(q-1)+j )  + s * Mass(IMap(i),IMap(j))
            END DO
         END DO
       ELSE
         DO i=1,c
           DO j=1,c
             StiffMatrix( c*(p-1)+i, c*(q-1)+j ) = &
                 StiffMatrix( c*(p-1)+i, c*(q-1)+j ) + s * A(i,j)

             MassMatrix(  c*(p-1)+i, c*(q-1)+j ) = &
                 MassMatrix(  c*(p-1)+i, c*(q-1)+j ) + s * Mass(i,j)
            END DO
         END DO
       END IF
 
     END DO
     END DO

!
!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------
!
     IF ( NewtonLinearization ) THEN
       DO i=1,dim
         DO j=1,dim
           Force(i) = Force(i) + dVelodx(i,j) * Uvelo(j)
           IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             DO k=1,dim
               Force(i) = Force(i) + Symb(k,j,i) * Uvelo(k) * Uvelo(j)
             END DO
           END IF
         END DO
       END DO
     END IF 

     DO p=1,N
       Load = 0.0D0
       DO i=1,c
          Load(i) = Load(i) + Density * Force(i) * Basis(p)
       END DO

       IF ( Stabilize ) THEN
         DO i=1,c
           DO j=1,c
             Load(j) = Load(j) + Tau * Density * Force(i) * SW(p,i,j)
           END DO
         END DO
       END IF

       IF ( CurrentCoordinateSystem() == AxisSymmetric ) THEN
         DO i=1,3
           ForceVector(3*(p-1)+i) = ForceVector(3*(p-1)+i) + s*Load(IMap(i))
         END DO
       ELSE
         DO i=1,c
           ForceVector(c*(p-1)+i) = ForceVector(c*(p-1)+i) + s*Load(i)
         END DO
       END IF
     END DO

   END DO 
!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesGeneralCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for Navier-Stokes-equations
!>  boundary conditions. (No velocity dependent velocity BC:s ("Newton BCs")
!>  at the moment, so BoundaryMatrix will contain only zeros at exit...)
!------------------------------------------------------------------------------
 SUBROUTINE NavierStokesGeneralBoundary( BoundaryMatrix,BoundaryVector, &
   LoadVector,NodalAlpha,NodalBeta,NodalExtPressure,NodalSlipCoeff,Element,n,Nodes )
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
       NodalAlpha(:),NodalBeta(:),NodalSlipCoeff(:,:), NodalExtPressure(:)

   TYPE(Element_t),POINTER  :: Element
   TYPE(Nodes_t)    :: Nodes

   INTEGER :: n

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------

   REAL(KIND=dp) :: Basis(n)
   REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

   REAL(KIND=dp) :: Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
   REAL(KIND=dp) :: u,v,w,s,x,y,z,SlipCoeff
   REAL(KIND=dp) :: Force(3),Alpha,Normal(3)
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   INTEGER :: i,j,t,q,p,c,dim,N_Integ

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

   IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
     dim = 3
   ELSE
     dim = CoordinateSystemDimension()
   ENDIF
   c = dim + 1

   BoundaryVector = 0.0D0
   BoundaryMatrix = 0.0D0
!
!------------------------------------------------------------------------------
!  Integration stuff
!------------------------------------------------------------------------------
!
   IntegStuff = GaussPoints( element )
   U_Integ => IntegStuff % u
   V_Integ => IntegStuff % v
   W_Integ => IntegStuff % w
   S_Integ => IntegStuff % s
   N_Integ =  IntegStuff % n
!
!------------------------------------------------------------------------------
!  Now we start integrating
!------------------------------------------------------------------------------
!
   DO t=1,N_Integ

     u = U_Integ(t)
     v = V_Integ(t)
     w = W_Integ(t)
!
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx )
!------------------------------------------------------------------------------
!    Coordinatesystem dependent info
!------------------------------------------------------------------------------
!
     IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
       x = SUM( nodes % x(1:n)*Basis )
       y = SUM( nodes % y(1:n)*Basis )
       z = SUM( nodes % z(1:n)*Basis )
     END IF

     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
!
     s = SqrtMetric * SqrtElementMetric * S_Integ(t)
!
!------------------------------------------------------------------------------
!    Add to load: tangential derivative of something
!------------------------------------------------------------------------------
     DO i=1,3
       Force(i) = SUM( NodalBeta(1:n)*dBasisdx(1:n,i) )
     END DO
!
!------------------------------------------------------------------------------
!    Add to load: given force in normal direction
!------------------------------------------------------------------------------
!
     Normal = NormalVector( Element,Nodes,u,v,.TRUE. )

     Alpha  = SUM( NodalExtPressure(1:n)*Basis )
     DO i=1,dim
       DO j=1,dim
         Force(i) = Force(i) + Alpha*Metric(i,j)*Normal(j)
       END DO
     END DO

     Alpha  = SUM( NodalAlpha(1:n)*Basis )
!
!------------------------------------------------------------------------------
!    Add to load: given force in coordinate directions
!------------------------------------------------------------------------------
!
     DO i=1,3
       Force(i) = Force(i) + SUM( LoadVector(i,1:n) )
     END DO

     IF ( ANY( NodalSlipCoeff(:,:) /= 0.0d0 ) ) THEN
       DO p=1,n
         DO q=1,n
           DO i=1,dim
             SlipCoeff = SUM( NodalSlipCoeff(i,1:n) * Basis(1:n) )
             BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) = &
               BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) + s * SlipCoeff * Basis(q) * Basis(p)
           END DO
         END DO
       END DO
     END IF

     DO q=1,N
       DO i=1,dim
         BoundaryVector( (q-1)*c+i ) = BoundaryVector( (q-1)*c+i ) + &
                     s * Basis(q) * Force(i)
       END DO
     END DO
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesGeneralBoundary
!------------------------------------------------------------------------------


END MODULE NavierStokesGeneral

!> \}
