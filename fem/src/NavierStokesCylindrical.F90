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
!>  Module computing Navier-stokes local matrices in cylindrical coordinates.
!------------------------------------------------------------------------------
MODULE NavierStokesCylindrical

  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE ElementDescription, ONLY: GetEdgeMap

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for Navier-Stokes-Equations
!>  in axisymmetric, cylindric symmetric or cylindrical coordinates.
!------------------------------------------------------------------------------
   SUBROUTINE NavierStokesCylindricalCompose  (                             &
       MassMatrix,StiffMatrix,ForceVector,LoadVector,NodalViscosity,   &
       NodalDensity,Ux,Uy,Uz,MUx,MUy,MUz,NodalPressure,NodalTemperature, &
       Convect,StabilizeFlag,Compressible,PseudoCompressible, NodalCompressibility, &
       NodalGasConstant,Porous, NodalDrag, PotentialForce, PotentialField, &
       PotentialCoefficient, MagneticForce, divDiscretization, gradPDiscretization, &
       NewtonLinearization,Element,n,Nodes )
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
!            power-law viscosity is used, the values vary in space)
!
!  REAL(KIND=dp) :: NodalDensity(:)
!     INPUT: nodal density values
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!
!  REAL(KIND=dp) :: NodalPressure(:)
!     INPUT: Nodal values of total pressure from previous iteration
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ?
!
!  LOGICAL :: Compressible
!     INPUT: Should compressible terms be added =
!
!  LOGICAL :: PseudoCompressible
!     INPUT: Should artificial compressibility be added ?
!
!  REAL(KIND=dp) :: NodalCompressibility(:)
!     INPUT: Artificial compressibility for the nodes
!
!  LOGICAL :: MagneticForce
!      INPUT: Should Lorentz force for magnetohydrodynamics be included
!
!  LOGICAL :: NewtonLinearization
!      INPUT: Picard or Newton  linearization of the convection term ?
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number  of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Ux,Uy,Uz,MUx,MUy,MUz, &
                     NodalPressure(:), NodalTemperature(:)
     REAL(KIND=dp) :: NodalViscosity(:),NodalDensity(:),LoadVector(:,:), &
         NodalCompressibility(:), NodalDrag(:,:), PotentialField(:), &
         PotentialCoefficient(:), NodalGasConstant(:)
     LOGICAL :: Convect,NewtonLinearization,Compressible,&
         PseudoCompressible, Porous, PotentialForce, &
         MagneticForce, divDiscretization, gradPDiscretization

      CHARACTER(LEN=*) :: StabilizeFlag

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Velo(3),UVelo(3),dVelodx(3,3),Force(4)

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, dNodalBasisdx(n,n,3)
     REAL(KIND=dp) :: Basis(2*n), PBasis(n), PdBasisdx(n,3),BaseP

     REAL(KIND=dp) :: GMetric(3,3),Metric(3),SqrtMetric,&
                            Symb(3,3,3),dSymb(3,3,3,3), Drag(3)

     REAL(KIND=dp), DIMENSION(4,4) :: A,Mass
     REAL(KIND=dp) :: Load(4),SU(n,4,4),SW(N,4,4),LrF(3)

     REAL(KIND=dp) :: Lambda=1.0,Re,Tau,Delta,x0,y0
     REAL(KIND=dp) :: VNorm,hK,mK,Viscosity,dViscositydx(3),Density
     REAL(KIND=dp) :: Pressure,Temperature,dTemperaturedx(3), &
         dDensitydx(3), Compress, dPressuredx(3), GasConstant

     INTEGER :: i,j,k,l,m,c,p,q,t,DIM,N_Integ,NBasis

     REAL(KIND=dp) :: s,u,v,w,x,y,z
  
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
 
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER, POINTER :: EdgeMap(:,:)

     TYPE(ElementType_t), POINTER :: LinearType, SaveType
     INTEGER :: LinearBasis, LinearCode(3:8) = (/ 303,404,504,605,706,808 /)
  
     LOGICAL :: stat,CylindricSymmetry,Stabilize,Bubbles,PBubbles,P2P1

     INTEGER :: IMap(3) = (/ 1,2,4 /)
!------------------------------------------------------------------------------

     CylindricSymmetry = ( CurrentCoordinateSystem() == CylindricSymmetric )

! Assembly assuming cylindric symmetry and then pick the entries for axi symmetric case
! Hence dimension is assumed to be zero
!--------------------------------------------------------------------------------------
     dim = 3
     c = dim + 1
     N = element % TYPE % NumberOfNodes

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     NBasis    = n
     Bubbles   = .FALSE.
     PBubbles  = .FALSE.
     P2P1      = .FALSE.
     Stabilize = StabilizeFlag == 'stabilized'
     IF ( .NOT.Stabilize .OR. Compressible ) THEN
       PBubbles = isActivePelement(Element) .AND. Element % BDOFs > 0
       IF ( PBubbles ) THEN
          NBasis = n + Element % BDOFs
       ELSE IF ( StabilizeFlag == 'bubbles' .OR. &
              Element % TYPE % BasisFunctionDegree<=1 ) THEN
          NBasis    = 2 * n
          Bubbles   = .TRUE.
       ELSE
          P2P1 = .TRUE.
       END IF
       Stabilize = .FALSE.
     END IF

     LinearBasis = 0
     IF ( P2P1 ) THEN
        SaveType => Element % TYPE
        j = GetElementFamily()
        LinearType => GetElementType(LinearCode(j))
        LinearBasis = LinearType % NumberOfNodes
     END IF

     IF ( Bubbles ) THEN
       IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints2 )
     ELSE
       IntegStuff = GaussPoints( Element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ  = IntegStuff % n

!------------------------------------------------------------------------------
!    Stabilization parameter mK
!------------------------------------------------------------------------------
     IF ( Stabilize ) THEN
       hK = element % hK
       mK = element % StabilizationMK

       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % TYPE % NodeU(p)
         v = Element % TYPE % NodeV(p)
         w = Element % TYPE % NodeW(p)
         stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
     END IF
!
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
    DO t=1,N_Integ

!------------------------------------------------------------------------------
!     Integration stuff
!------------------------------------------------------------------------------
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      IF ( P2P1 ) THEN
         Element % TYPE => LinearType
         stat = ElementInfo( Element, Nodes, u, v, w, detJ,pBasis,pdBasisdx )
         Element % TYPE => SaveType
      END IF

      stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
             Basis,dBasisdx, Bubbles=Bubbles )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      x = SUM( Nodes % x(1:n) * Basis(1:n) )
      y = SUM( Nodes % y(1:n) * Basis(1:n) )
      z = SUM( Nodes % z(1:n) * Basis(1:n) )

      CALL CoordinateSystemInfo( GMetric,SqrtMetric,Symb,dSymb,x,y,z )
      DO i=1,3
        Metric(i) = GMetric(i,i)
      END DO
      s = SqrtMetric * detJ * S_Integ(t)

!------------------------------------------------------------------------------
!     Density at the integration point
!------------------------------------------------------------------------------
      Density = SUM( NodalDensity(1:n)*Basis(1:n) )
      IF ( Compressible ) THEN
        IF ( P2P1 ) THEN
          k = LinearBasis
          Temperature = SUM( NodalTemperature(1:k)*pBasis(1:k) )
          Pressure = SUM( NodalPressure(1:k)*pBasis(1:k) )
          DO i=1,dim
            dDensitydx(i)     = SUM( NodalDensity(1:k)*pdBasisdx(1:k,i) )
            dPressuredx(i)    = SUM( NodalPressure(1:k)*pdBasisdx(1:k,i) )
            dTemperaturedx(i) = SUM( NodalTemperature(1:k)*pdBasisdx(1:k,i) )
          END DO
        ELSE
          Temperature = SUM( NodalTemperature(1:n)*Basis(1:n) )
          Pressure = SUM( NodalPressure(1:n)*Basis(1:n) )
          DO i=1,dim
            dDensitydx(i)     = SUM( NodalDensity(1:n)*dBasisdx(1:n,i) )
            dPressuredx(i)    = SUM( NodalPressure(1:n)*dBasisdx(1:n,i) )
            dTemperaturedx(i) = SUM( NodalTemperature(1:n)*dBasisdx(1:n,i) )
          END DO
        END IF

        GasConstant = SUM( Basis(1:n) * NodalGasConstant(1:n) )
        Density = Pressure / (Temperature*GasConstant)
      END IF

      IF(PseudoCompressible) THEN
        Pressure = SUM( NodalPressure(1:n) * Basis(1:n) )        
        Compress = Density * SUM(NodalCompressibility(1:n)*Basis(1:n))      
      END IF

!------------------------------------------------------------------------------
!     Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
      Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
      Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
      IF ( CylindricSymmetry ) THEN
        Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )
      ELSE
        Velo(3) = 0.0_dp
      END IF

      IF ( NewtonLinearization ) THEN
        dVelodx = 0.0D0
        DO i=1,3
          dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
          dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
          IF ( CylindricSymmetry ) THEN
            dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
          END IF
        END DO
      END IF
!  
!------------------------------------------------------------------------------
!     Force at the integration point
!------------------------------------------------------------------------------
      Lrf = 0.0d0
      IF ( MagneticForce ) THEN
         Lrf = LorentzForce( Element,Nodes,u,v,w,n )
      END IF

      Force = 0.0D0
! Use this, if the (time-averaged) Lorentz force is not divided by density
! and its phi-component is given in SI units
#ifdef LORENTZ_AVE
      Force(1) = SUM( LoadVector(1,1:n)*Basis(1:n) )/Density
      Force(2) = SUM( LoadVector(2,1:n)*Basis(1:n) )/Density
      IF( CylindricSymmetry ) THEN
        Force(3) = SUM( LoadVector(3,1:n)*Basis(1:n) )/(Density*x)
        Force(4) = SUM( LoadVector(4,1:n)*Basis(1:n) )
      ELSE 
        Force(3) = 0.0_dp
        Force(4) = SUM( LoadVector(3,1:n)*Basis(1:n) )
      END IF
#else
      Force(1) = SUM( LoadVector(1,1:n)*Basis(1:n) )
      Force(2) = SUM( LoadVector(2,1:n)*Basis(1:n) )
      IF( CylindricSymmetry ) THEN
        Force(3) = SUM( LoadVector(3,1:n)*Basis(1:n) )
        Force(4) = SUM( LoadVector(4,1:n)*Basis(1:n) )
      ELSE
        Force(3) = 0.0_dp
        Force(4) = SUM( LoadVector(3,1:n)*Basis(1:n) )
      END IF
#endif

      IF ( CylindricSymmetry ) THEN
        Force(1:3) = Force(1:3) + LrF(1:3) / Density
      ELSE
        Force(1:2) = Force(1:2) + LrF(1:2) / Density
      END IF
      
      IF ( Compressible ) Force(4) = Force(4) / Temperature

!------------------------------------------------------------------------------
!     Additional forces due to gradient forces (electrokinetic flow) and
!     viscous drag in porous media.
!------------------------------------------------------------------------------

      IF(PotentialForce) THEN        
        Force(1) = Force(1) - SUM( PotentialCoefficient(1:n) * Basis(1:n) ) * &
            SUM(  PotentialField(1:n) * dBasisdx(1:n,1) )
        Force(2) = Force(2) - SUM( PotentialCoefficient(1:n) * Basis(1:n) ) * &
            SUM(  PotentialField(1:n) * dBasisdx(1:n,2) ) / x
      END IF
      
      IF(Porous) THEN
        DO i=1,DIM
          Drag(i) = SUM( NodalDrag(i,1:n) * Basis(1:n) )
        END DO
      END IF

!------------------------------------------------------------------------------
!     Effective viscosity & derivatives at integration point
!------------------------------------------------------------------------------
      Viscosity = SUM( NodalViscosity(1:n) * Basis(1:n) )
      Viscosity = EffectiveViscosity( Viscosity, Density, Ux, Uy, Uz, &
                    Element, Nodes, n, n, u, v, w, LocalIP=t )

!------------------------------------------------------------------------------
!      Stabilization parameters Tau & Delta
!------------------------------------------------------------------------------
     IF ( Stabilize ) THEN
!------------------------------------------------------------------------------
       DO i=1,3
         dViscositydx(i) = SUM( NodalViscosity(1:n)*dBasisdx(1:n,i) )
       END DO
!------------------------------------------------------------------------------
!      VNorm = SQRT( SUM(Velo(1:dim)**2) )
 
       IF ( Convect ) THEN
         Vnorm = 0.0D0
         DO i=1,DIM
            Vnorm = Vnorm + Velo(i) * Velo(i) / Metric(i)
         END DO
         Vnorm = MAX( SQRT( Vnorm ), 1.0d-12 )
 
         Re = MIN( 1.0D0, Density * mK * hK * VNorm / (4 * Viscosity) )

         Tau = 0.0D0
         IF ( VNorm /= 0.0D0 ) THEN
           Tau = hK * Re / (2 * Density * VNorm)
         END IF

         Delta = Density * Lambda * Re * hK * VNorm
       ELSE
         Delta = 0._dp
         Tau = mK * hK**2 / ( 8 * Viscosity )
       END IF

!------------------------------------------------------------------------------
!      SU will contain residual of ns-equations, SW will contain the
!      weight function terms
!------------------------------------------------------------------------------
       SU = 0.0D0
       SW = 0.0D0
       DO p=1,N
         DO i=1,dim
           SU(p,i,c) = SU(p,i,c) + Metric(i) * dBasisdx(p,i)

           IF(Porous) THEN
             SU(p,i,i) = SU(p,i,i) + Viscosity * Drag(i) * Basis(p)
           END IF

           IF ( Convect ) THEN
             DO j=1,DIM
               SU(p,i,i) = SU(p,i,i) + Density * dBasisdx(p,j) * Velo(j)
               IF ( NewtonLinearization ) THEN
                 SU(p,i,j) = SU(p,i,j) + Density * dVelodx(i,j) * Basis(p)
               END IF
               DO k =1,dim,2
                 SU(p,i,k) = SU(p,i,k) + Density * Symb(k,j,i) * Basis(p) * Velo(j)
                 IF ( NewtonLinearization ) THEN
                   SU(p,i,j) = SU(p,i,j) + Density * Symb(k,j,i) * Velo(k) * Basis(p)
                 END IF
               END DO
             END DO
           END IF

           !
           ! diffusion
           ! ----------
           DO j=1,dim
             SU(p,i,i) = SU(p,i,i) - dViscositydx(j) * Metric(j) * dBasisdx(p,j)
             SU(p,i,j) = SU(p,i,j) - dViscositydx(j) * Metric(i) * dBasisdx(p,i)
             SU(p,i,i) = SU(p,i,i) - Viscosity * Metric(j) * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
             SU(p,i,j) = SU(p,i,j) - Viscosity * Metric(i) * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,i))
           END DO

           DO j=1,dim,2
             DO l=1,dim,2
               SU(p,i,i) = SU(p,i,i) + Viscosity * Metric(j) * Symb(j,j,l) * dBasisdx(p,l)
               SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j) * Symb(l,j,i) * dBasisdx(p,j)
               SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j) * Symb(l,j,i) * dBasisdx(p,j)
               SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j) * dSymb(l,j,i,j) * Basis(p)
               SU(p,i,j) = SU(p,i,j) + Viscosity * Metric(i) * Symb(j,i,l) * dBasisdx(p,l)
               SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i) * Symb(l,j,j) * dBasisdx(p,i)
               SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i) * Symb(l,i,j) * dBasisdx(p,j)
               SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i) * dSymb(l,j,j,i) * Basis(p)
               DO m=1,DIM,2
                 SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(j) * Symb(j,m,i) * Symb(l,j,m) * Basis(p)
                 SU(p,i,l) = SU(p,i,l) + Viscosity * Metric(j) * Symb(j,j,m) * Symb(l,m,i) * Basis(p)
                 SU(p,i,l) = SU(p,i,l) - Viscosity * Metric(i) * Symb(m,i,j) * Symb(l,j,m) * Basis(p)
                 SU(p,i,l) = SU(p,i,l) + Viscosity * Metric(i) * Symb(j,i,m) * Symb(l,m,j) * Basis(p)
               END DO
             END DO
! then -mu,_j (g^{jk} U^i_,k + g^{ik} U^j_,k)
             DO l=1,DIM,2
               SU(p,i,l) = SU(p,i,l) - dViscositydx(j) * &
                 ( Metric(j) * Basis(p) * Symb(j,l,i) + Metric(i) * Basis(p) * Symb(i,l,j) )
             END DO
           END DO

!------------------------------------------------------------------------------
           IF ( Convect ) THEN
             SW(p,i,c) = SW(p,i,c) + Density * dBasisdx(p,i)
             DO j=1,DIM
               SW(p,i,i) = SW(p,i,i) + Density * dBasisdx(p,j) * Velo(j)
               DO k =1,DIM,2
                 SW(p,i,k) = SW(p,i,k) - Density * Symb(i,j,k) * Basis(p) * Velo(j)
               END DO
             END DO
           ELSE
             SW(p,i,c) = SW(p,i,c) + dBasisdx(p,i)
           END IF

           !
           !  Diffusion
           !
           DO j=1,dim
             SW(p,i,i) = SW(p,i,i) + dViscositydx(j) * Metric(j) * dBasisdx(p,j)
             SW(p,i,j) = SW(p,i,j) + dViscositydx(j) * Metric(j) * dBasisdx(p,i)
             SW(p,i,i) = SW(p,i,i) - Viscosity * Metric(j) * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
             SW(p,i,j) = SW(p,i,j) - Viscosity * Metric(j) * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,i))
           END DO

           DO j=1,dim,2
              DO l=1,dim,2
                SW(p,i,i) = SW(p,i,i) - Viscosity * Metric(j) * Symb(j,j,l) * dBasisdx(p,l)
                SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j) * Symb(i,j,l) * dBasisdx(p,j)
                SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j) * Symb(i,j,l) * dBasisdx(p,j)
                SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j) * dSymb(i,j,l,j) * Basis(p)
                SW(p,i,j) = SW(p,i,j) - Viscosity * Metric(j) * Symb(i,j,l) * dBasisdx(p,l)
                SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j) * Symb(i,j,l) * dBasisdx(p,j)
                SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j) * Symb(j,j,l) * dBasisdx(p,i)
                SW(p,i,l) = SW(p,i,l) - Viscosity * Metric(j) * dSymb(i,j,l,j) * Basis(p)

                DO m=1,dim,2
                  SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j) * Symb(i,j,m) * Symb(m,j,l) * Basis(p)
                  SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j) * Symb(j,j,m) * Symb(m,i,l) * Basis(p)
                  SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j) * Symb(j,j,m) * Symb(m,i,l) * Basis(p)
                  SW(p,i,l) = SW(p,i,l) + Viscosity * Metric(j) * Symb(i,j,m) * Symb(m,j,l) * Basis(p)
                END DO
              END DO
! then -mu,_j g^{jk} (w_i,_k +  w_k,_i)
              DO l=1,dim,2
                SW(p,i,l) = SW(p,i,l) + dViscositydx(j) * &
                  ( Metric(j) * Basis(p) * Symb(i,j,l) + Metric(i) * Basis(p) * Symb(j,i,l) )
              END DO
           END DO

           IF ( .NOT. CylindricSymmetry ) THEN
             SU(p,i,3) = 0.0_dp
             SW(p,i,3) = 0.0_dp
           END IF
         END DO
       END DO
     END IF

!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
     DO p=1,NBasis
     DO q=1,NBasis
!
!------------------------------------------------------------------------------
! First plain Navier-Stokes
!------------------------------------------------------------------------------

       BaseP = Basis(p)
       IF ( P2P1 .AND. p<=LinearBasis ) THEN
          BaseP = PBasis(p)
       ELSE
          BaseP = Basis(p)
       END IF

!------------------------------------------------------------------------------
!      Mass Matrix:
! -----------------------------------
!------------------------------------------------------------------------------
       Mass = 0.0D0
       DO i=1,DIM
         Mass(i,i) = Density * Basis(q) * Basis(p)
       END DO

! Continuity equation
       IF ( Compressible ) THEN
         Mass(c,c) = ( Density / Pressure ) * Basis(q) * BaseP
       END IF
!------------------------------------------------------------------------------

       A = 0.0
!
!------------------------------------------------------------------------------
!      Stiffness Matrix:
!------------------------------
! Possible Porous media effects
!------------------------------------------------------------------------------
       IF(Porous) THEN
         DO i=1,DIM
           A(i,i) = A(i,i) + Viscosity * Drag(i) * Basis(q) * Basis(p)
         END DO
       END IF

! -----------------------------------
!      Diffusive terms
!------------------------------------------------------------------------------
       DO i=1,dim
         DO j = 1,dim
           A(i,i) = A(i,i) + Viscosity * Metric(j) * dBasisdx(q,j) * dBasisdx(p,j)
           A(i,j) = A(i,j) + Viscosity * Metric(i) * dBasisdx(q,i) * dBasisdx(p,j)

           IF ( Compressible ) THEN
!------------------------------------------------------------------------------
!  For compressible flows add (2/3) \mu \nabla \cdot u
!  Partial derivative terms only here
!------------------------------------------------------------------------------
             A(i,j) = A(i,j) - (2.0d0/3.0d0) * Viscosity * Metric(i) * &
                           dBasisdx(q,j) * dBasisdx(p,i)
           END IF
         END DO
       END DO

!------------------------------------------------------------------------------
       IF ( Compressible ) THEN
!------------------------------------------------------------------------------
!  For compressible flows add (2/3) \mu \nabla \cdot u
!  Terms involving Christoffel symbols
!------------------------------------------------------------------------------
         DO i=1,DIM
           A(i,1) = A(i,1) - ( 2.0d0 /3.0d0 ) * Viscosity * Metric(i) * &
                    Basis(q) * Symb(3,1,3) * dBasisdx(p,i)
         END DO

         DO j=1,DIM
           A(1,j) = A(1,j) + ( 2.0d0 / 3.0d0 ) * Viscosity * dBasisdx(q,j) * &
                    Metric(3) * Symb(3,3,1) * Basis(p)
         END DO

         A(1,1) = A(1,1) + ( 2.0d0 / 3.0d0 ) * Viscosity * Symb(3,1,3) * &
                  Basis(q) * Metric(3) * Symb(3,3,1) * Basis(p)
       END IF
!------------------------------------------------------------------------------

       IF ( .NOT. CylindricSymmetry ) THEN
         A(1,3) = A(1,3) - 2 * Viscosity * Metric(3) * Symb(3,3,1) * Basis(p) * dBasisdx(q,3)
         A(3,1) = A(3,1) + 2 * Viscosity * Metric(3) * Symb(1,3,3) * Basis(q) * dBasisdx(p,3)
         A(3,1) = A(3,1) - 2 * Viscosity * Metric(3) * Symb(1,3,3) * Basis(p) * dBasisdx(q,3)
       END IF
       A(1,1) = A(1,1) + 2 * Viscosity * Metric(3) * Basis(q) * Basis(p)
       A(3,3) = A(3,3) - 2 * Viscosity * Symb(3,1,3) * dBasisdx(q,1) * Basis(p)

!------------------------------------------------------------------------------
!      Convection terms, Picard linearization
!------------------------------------------------------------------------------

       IF ( Convect ) THEN
          DO i = 1,dim
            DO j = 1,dim
              A(i,i) = A(i,i) + Density * dBasisdx(q,j) * Velo(j) * Basis(p)
            END DO
          END DO
          A(1,3) = A(1,3) + Density * Symb(3,3,1) * Basis(q) * Velo(3) * Basis(p)
          A(3,1) = A(3,1) + Density * Symb(1,3,3) * Basis(q) * Velo(3) * Basis(p)
          A(3,3) = A(3,3) + Density * Symb(3,1,3) * Basis(q) * Velo(1) * Basis(p)
 
!------------------------------------------------------------------------------
!      Convection terms, Newton linearization
!------------------------------------------------------------------------------
!
          IF ( NewtonLinearization ) THEN
            DO i=1,dim
              DO j=1,dim
                A(i,j) = A(i,j) + Density * dVelodx(i,j) * Basis(q) * Basis(p)
              END DO
            END DO
            A(1,3) = A(1,3) + Density * Symb(3,3,1) * Basis(q) * Velo(3) * Basis(p)
            A(3,1) = A(3,1) + Density * Symb(3,1,3) * Basis(q) * Velo(3) * Basis(p)
            A(3,3) = A(3,3) + Density * Symb(1,3,3) * Basis(q) * Velo(1) * Basis(p)
          END IF
       END IF
!
!------------------------------------------------------------------------------
!      Pressure terms
!------------------------------------------------------------------------------
!
       IF ( gradpDiscretization ) THEN
         DO i=1,dim
           A(i,c) = A(i,c) + Metric(i) * dBasisdx(q,i) * Basis(p)
         END DO
       ELSE
         DO i=1,dim
           A(i,c) = A(i,c) - Metric(i) * Basis(q) * dBasisdx(p,i)
         END DO
         A(1,c) = A(1,c) + Metric(3) * Basis(q) * Symb(3,3,1) * Basis(p)  
       END IF

!
!------------------------------------------------------------------------------
!      Continuity equation
!------------------------------------------------------------------------------
       DO i=1,dim
!------------------------------------------------------------------------------
         IF ( Compressible ) THEN
           A(c,i) = A(c,i) + (Density / Pressure) * Basis(q) * &
                      dPressuredx(i) * BaseP / 2

           A(c,c) = A(c,c) + (Density / Pressure) * Velo(i)  * &
                      dBasisdx(q,i) * BaseP / 2

           A(c,i) = A(c,i) - (Density/Temperature) * &
                 Basis(q) * dTemperaturedx(i) * BaseP

           A(c,i) = A(c,i) + Density * dBasisdx(q,i) * BaseP
         ELSE
           IF ( Convect ) THEN
             IF ( gradpDiscretization ) THEN
               A(c,i) = A(c,i) - Density * Basis(q) * dBasisdx(p,i)
             ELSE
               A(c,i) = A(c,i) + Density * dBasisdx(q,i) * BaseP
             END IF
           ELSE
             IF ( gradpDiscretization ) THEN
               A(c,i) = A(c,i) - Basis(q) * dBasisdx(p,i)
             ELSE
               A(c,i) = A(c,i) + dBasisdx(q,i) * BaseP
             END IF
           END IF
         END IF
!------------------------------------------------------------------------------
       END DO

       IF ( .NOT. gradpDiscretization ) THEN
         IF ( Compressible .OR. Convect ) THEN
           A(c,1) = A(c,1) + Density * Symb(1,3,3) * Basis(q) * BaseP
         ELSE
           A(c,1) = A(c,1) + Symb(1,3,3) * Basis(q) * BaseP
         END IF
       END IF
  
!------------------------------------------------------------------------------
!      Artificial Compressibility, affects only the continuity equation
!------------------------------------------------------------------------------  

       IF(PseudoCompressible) THEN
         A(c,c) = A(c,c) + Compress * Basis(q) * Basis(p)
       END IF


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
                A(j,i) = A(j,i) + Delta * dBasisdx(q,i) * Metric(j) * dBasisdx(p,j)
                DO l=1,dim,2
                   A(l,i) = A(l,i) - Delta * dBasisdx(q,i) * Metric(j) * Symb(j,j,l) * Basis(p)

                   A(j,l) = A(j,l) + Delta * Symb(l,i,i) * Basis(q) * Metric(j) * dBasisdx(p,j)
                   DO m=1,dim,2
                      A(m,l) = A(m,l) - Delta * Symb(l,i,i) * Basis(q) * Metric(j) * Symb(j,j,m) * Basis(p)
                   END DO
                END DO
             END DO
          END DO
       END IF

!
!------------------------------------------------------------------------------
! Add nodal matrix to element matrix
!------------------------------------------------------------------------------
!
       IF ( CylindricSymmetry ) THEN
         DO i=1,c
           DO j=1,c
             StiffMatrix( c*(p-1)+i,c*(q-1)+j ) = &
                 StiffMatrix( c*(p-1)+i,c*(q-1)+j ) + s*A(i,j)

             MassMatrix(  c*(p-1)+i,c*(q-1)+j ) = &
                 MassMatrix(  c*(p-1)+i,c*(q-1)+j ) + s*Mass(i,j)
           END DO
         END DO
       ELSE
         DO i=1,3
           DO j=1,3
             StiffMatrix( 3*(p-1)+i,3*(q-1)+j ) = &
                 StiffMatrix( 3*(p-1)+i,3*(q-1)+j ) + s*A(IMap(i),IMap(j))
             
             MassMatrix(  3*(p-1)+i,3*(q-1)+j ) =  &
                MassMatrix(  3*(p-1)+i,3*(q-1)+j ) + s*Mass(IMap(i),IMap(j))
           END DO
         END DO
       END IF
 
     END DO
     END DO

!
!------------------------------------------------------------------------------
! The righthand side...
!------------------------------------------------------------------------------
!
     IF (  Convect .AND. NewtonLinearization ) THEN

       Uvelo(1) = SUM( Basis(1:n) * Ux(1:n) )
       Uvelo(2) = SUM( Basis(1:n) * Uy(1:n) )
       IF ( CylindricSymmetry ) THEN
         Uvelo(3) = SUM( Basis(1:n) * Uz(1:n) )
       ELSE
         Uvelo(3) = 0.0_dp
       END IF

       DO i=1,dim
         DO j=1,dim
           Force(i) = Force(i) + dVelodx(i,j) * Uvelo(j)
         END DO
       END DO

       IF ( CylindricSymmetry ) THEN
         Force(1) = Force(1) + Symb(3,3,1) * UVelo(3) * UVelo(3)
         Force(3) = Force(3) + Symb(3,1,3) * UVelo(1) * UVelo(3)
         Force(3) = Force(3) + Symb(1,3,3) * UVelo(3) * UVelo(1)
       END IF
     END IF 

     DO p=1,NBasis
       Load = 0.0d0
       DO i=1,c
         Load(i) = Load(i) + Density * Force(i) * Basis(p)
       END DO

       IF(PseudoCompressible) THEN
         Load(c) = Load(c) + Pressure * Basis(p) * Compress
       END IF

       IF ( Stabilize ) THEN
          DO i=1,DIM
            DO j=1,c
              Load(j) = Load(j) + Tau * Density * Force(i) * SW(p,i,j)
            END DO
          END DO
       END IF

       IF ( CylindricSymmetry ) THEN
         DO i=1,c
           ForceVector( c*(p-1)+i ) = ForceVector( c*(p-1)+i ) + s*Load(i)
         END DO
       ELSE
         DO i=1,3
           ForceVector(3*(p-1)+i) = ForceVector(3*(p-1)+i) + s*Load(IMap(i))
         END DO
       END IF
     END DO

   END DO 

   IF ( CylindricSymmetry ) THEN
     k = 4
   ELSE
     k = 3
   END IF

   IF ( P2P1 ) THEN
     j = GetElementFamily()
     EdgeMap => GetEdgeMap(j)

     DO i=j+1,j+SIZE(EdgeMap,1)
       p = EdgeMap(i-j,1)
       q = EdgeMap(i-j,2)
       StiffMatrix( k*i, : ) = 0.0d0
       MassMatrix(  k*i, : ) = 0.0d0
       ForceVector( k*i ) = 0.0d0
       StiffMatrix( k*i, k*i ) =  1.0d0
       StiffMatrix( k*i, k*p ) = -1.0d0/2.0d0
       StiffMatrix( k*i, k*q ) = -1.0d0/2.0d0
     END DO
   END IF

   IF ( PBubbles ) THEN
      DO i=n+1,nBasis
        StiffMatrix( k*i, : ) = 0.0d0
        MassMatrix(  k*i, : ) = 0.0d0
        StiffMatrix( :, k*i ) = 0.0d0
        MassMatrix(  :, k*i ) = 0.0d0
        ForceVector( k*i ) = 0.0d0
        StiffMatrix( k*i, k*i ) =  1.0d0
     END DO
   END IF

!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesCylindricalCompose
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for Navier-Stokes-equations
!>  boundary conditions. (No velocity dependent velocity BC:s ("Newton BCs")
!>  at the moment, so BoundaryMatrix will contain only zeros at exit...).
!------------------------------------------------------------------------------
 SUBROUTINE NavierStokesCylindricalBoundary( BoundaryMatrix,BoundaryVector, &
     LoadVector,NodalAlpha,NodalBeta,NodalExtPressure,NodalSlipCoeff,NormalTangential, &
         Element,n,Nodes )
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

   USE ElementUtils

   IMPLICIT NONE

   REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:),LoadVector(:,:), &
                 NodalAlpha(:),NodalBeta(:),NodalSlipCoeff(:,:), NodalExtPressure(:)

   TYPE(Element_t),POINTER :: Element
   TYPE(Nodes_t)    :: Nodes

   INTEGER :: n

   LOGICAL :: NormalTangential

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n)
   REAL(KIND=dp) :: dBasisdx(n,3),detJ

   REAL(KIND=dp) :: Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
   REAL(KIND=dp) :: u,v,w,s,x,y,z
   REAL(KIND=dp) :: TangentForce(3), Force(4),Alpha,Normal(3),Tangent(3), &
                     Tangent2(3), Vect(3), SlipCoeff
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   INTEGER :: i,j,k,l,t,q,p,c,DIM,N_Integ

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

   IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
     dim = 3
   ELSE
     dim = CoordinateSystemDimension()
   END IF
   c = dim + 1

   BoundaryVector = 0.0_dp
   BoundaryMatrix = 0.0_dp
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
     stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                Basis,dBasisdx )
!
!------------------------------------------------------------------------------
!    Coordinatesystem dependent info
!------------------------------------------------------------------------------
!
     IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
       x = SUM( Nodes % x(1:n) * Basis(1:n) )
       y = SUM( Nodes % y(1:n) * Basis(1:n) )
       z = SUM( Nodes % z(1:n) * Basis(1:n) )
     END IF

     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
     s = SqrtMetric * detJ * S_Integ(t)
!
!------------------------------------------------------------------------------
!    Add to load: tangential derivative of something
!------------------------------------------------------------------------------
!
     TangentForce = 0._dp
     DO i=1,dim
       TangentForce(i) = SUM( NodalBeta(1:n)*dBasisdx(1:n,i) )
     END DO

!------------------------------------------------------------------------------
!    Add to load: given force in coordinate directions
!------------------------------------------------------------------------------
     Force = 0.0d0
     DO i=1,c
        Force(i) = Force(i) + SUM( LoadVector(i,1:n)*Basis(1:n) )
     END DO

!------------------------------------------------------------------------------
!    Add to load: given force in normal direction
!------------------------------------------------------------------------------
     Normal = NormalVector( Element, Nodes, u,v,.TRUE. )
     Alpha = SUM( NodalExtPressure(1:n)*Basis(1:n) )
     IF ( NormalTangential ) THEN
       Force(1) = Force(1) + Alpha
     ELSE
       DO i=1,dim
          Force(i) = Force(i) + Alpha * Normal(i)
       END DO
     END IF
!------------------------------------------------------------------------------

     Alpha = SUM( NodalAlpha(1:n)*Basis(1:n) )
!

     SELECT CASE( Element % TYPE % DIMENSION )
        CASE(1)
           Tangent(1) =  Normal(2)
           Tangent(2) = -Normal(1)
           Tangent(3) =  0.0_dp
           Tangent2   =  0.0_dp
        CASE(2)
           CALL TangentDirections( Normal, Tangent, Tangent2 ) 
     END SELECT

     IF ( ANY( NodalSlipCoeff(:,:) /= 0.0d0 ) ) THEN
        DO p=1,n
          DO q=1,n
            DO i=1,DIM

              SlipCoeff = SUM( NodalSlipCoeff(i,1:n) * Basis(1:n) )
              IF ( NormalTangential ) THEN
                 SELECT CASE(i)
                    CASE(1)
                      Vect = Normal
                    CASE(2)
                      Vect = Tangent
                    CASE(3)
                      Vect = Tangent2
                 END SELECT

                 DO j=1,dim
                    DO k=1,dim
                       BoundaryMatrix( (p-1)*c+j,(q-1)*c+k ) = &
                          BoundaryMatrix( (p-1)*c+j,(q-1)*c+k ) + &
                           s * SlipCoeff * Basis(q) * Basis(p) * Vect(j) * Vect(k)
                    END DO
                 END DO
              ELSE
                  BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) = &
                     BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) + &
                         s * SlipCoeff * Basis(q) * Basis(p)
              END IF
            END DO
          END DO
        END DO
     END IF

     DO q=1,N
       DO i=1,dim
          k = (q-1)*c + i
          IF ( NormalTangential ) THEN
             SELECT CASE(i)
                CASE(1)
                  Vect = Normal
                CASE(2)
                  Vect = Tangent
                CASE(3)
                  Vect = Tangent2
             END SELECT

             DO j=1,dim
                l = (q-1)*c + j
                BoundaryVector(l) = BoundaryVector(l) +  &
                             s * Basis(q) * Force(i) * Vect(j)
             END DO
          ELSE
             BoundaryVector(k) = BoundaryVector(k) + s * Basis(q) * Force(i)
          END IF
          BoundaryVector(k) = BoundaryVector(k) - s * Alpha * dBasisdx(q,i)
          BoundaryVector(k) = BoundaryVector(k) + s * TangentForce(i) * Basis(q)
       END DO
       k = (q-1)*c + 1
       BoundaryVector(k) = BoundaryVector(k) - s * Alpha * Basis(q) * Symb(3,1,3)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesCylindricalBoundary
!------------------------------------------------------------------------------

END MODULE NavierStokesCylindrical

!> \}
