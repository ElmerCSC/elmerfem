!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
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
! ******************************************************************************/

!----------------------------------------------------------------------
!> Diffuse-convective local matrix computing in cartesian coordinates.
!> Used nowadays only by HeatSolver.
!----------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{
MODULE DiffuseConvective

  USE DefUtils
  USE Differentials
  USE MaterialModels

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for diffusion-convection
!>  equation: 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveCompose( MassMatrix,StiffMatrix,ForceVector,  &
      LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,PhaseChange,NodalTemperature, &
      Enthalpy,Ux,Uy,Uz,MUx,MUy,MUz,Nodalmu,Nodalrho,NodalPressure, &
      NodaldPressureDt, NodalPressureCoeff, Compressible, Stabilize, &
      UseBubbles, Element,n,Nodes )
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
!  REAL(KIND=dp) :: NodalCT,NodalC0,NodalC1
!     INPUT: Coefficient of the time derivative term, 0 degree term, and
!            the convection term respectively
!
!  REAL(KIND=dp) :: NodalC2(:,:,:)
!     INPUT: Nodal values of the diffusion term coefficient tensor
!
!  LOGICAL :: PhaseChange
!     INPUT: Do we model phase change here...
!
!  REAL(KIND=dp) :: NodalTemperature
!     INPUT: NodalTemperature from previous iteration
!
!  REAL(KIND=dp) :: Enthalpy
!     INPUT: Enthalpy from previous iteration, needed if we model
!            phase change
!
!  REAL(KIND=dp) :: UX(:),UY(:),UZ(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values of the viscosity
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
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

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,UX,UY,UZ,MUX,MUY,MUZ,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: NodalTemperature(:),Enthalpy(:),Nodalmu(:), &
       NodaldPressureDt(:),NodalPressure(:),NodalPressureCoeff(:),Nodalrho(:)
     REAL(KIND=dp) :: NodalC0(:),NodalC1(:),NodalCT(:),NodalC2(:,:,:)

     LOGICAL :: UseBubbles,PhaseChange,Compressible,Stabilize, VectH

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     CHARACTER(LEN=MAX_NAME_LEN) :: StabilizeFlag
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: ddBasisddx(n,3,3),dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: Velo(3),Grad(3,3),Force

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK,hScale
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,Tau,x,y,z

     REAL(KIND=dp) :: Tau_M, Tau_C, Gmat(3,3),Gvec(3), dt=0._dp, LC1, NodalVelo(4,n), &
       RM(n),LC(3,n), gradP(n), VRM(3), GradNodal(n,4,3), PVelo(3), NodalPVelo(4,n), &
       Work(3,n), Grav(3), expc, reft, temperature

     REAL(KIND=dp), POINTER :: gWrk(:,:)

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis,Order

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s,u,v,w,dEnth,dTemp,mu,DivVelo,Pressure,rho,&
                      Pcoeff, minl, maxl, SSCond

     REAL(KIND=dp) :: C0,C00,C1,CT,CL,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: Vms, Found, Transient, stat,Convection,ConvectAndStabilize,Bubbles, &
          FrictionHeat, PBubbles
     TYPE(ValueList_t), POINTER :: BodyForce, Material
     LOGICAL :: GotCondModel
     
!------------------------------------------------------------------------------

     VectH = GetLogical( GetSolverParams(), 'VectH', Found )
     IF(.NOT. Found)  VectH = .FALSE.

     StabilizeFlag = GetString( GetSolverParams(),'Stabilization Method',Found )
     Vms = StabilizeFlag == 'vms'

     Material => GetMaterial()
     GotCondModel = ListCheckPresent( Material,'Heat Conductivity Model')


     dim = CoordinateSystemDimension()
     c = dim + 1

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     CL = 0

     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. (Vms .OR. Stabilize) .AND. UseBubbles ) THEN
       PBubbles  = isActivePElement(Element) .AND. Element % BDOFs > 0
       IF ( PBubbles ) THEN
          NBasis = n + Element % BDOFs
       ELSE
        NBasis = 2*n
        Bubbles = .TRUE.
       END IF
     END IF

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!    If there is no convection term we don t need stabilization.
!------------------------------------------------------------------------------
     hScale = GetCReal( GetSolverParams(), 'H scale', Found )
     IF(.NOT.Found) hScale = 1._dp

     hK = element % hK * hScale
     mK = element % StabilizationMK

     ConvectAndStabilize = .FALSE.
     IF  ( Vms .OR. VectH ) THEN
       NodalVelo(1,1:n) = Ux(1:n)
       NodalVelo(2,1:n) = Uy(1:n)
       NodalVelo(3,1:n) = Uz(1:n)
       NodalVelo(dim+1,1:n) = NodalPressure(1:n)

       dNodalBasisdx = 0._dp
       GradNodal = 0._dp
       DO p=1,n
          u = Element % TYPE % NodeU(p)
          v = Element % TYPE % NodeV(p)
          w = Element % TYPE % NodeW(p)
          stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
 
          dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
          GradNodal(p,1:dim,1:dim) = MATMUL( NodalVelo(1:dim,1:n), dBasisdx(1:n,1:dim) )
          GradNodal(p,dim+1,1:dim) = MATMUL( NodalVelo(dim+1,1:n), dBasisdx(1:n,1:dim) )
       END DO

       NodalPvelo = 0._dp

       ! transient flag only needed for vms
       Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'
       SSCond = ListGetConstReal(GetSolverParams(), "Steady State Condition", Found)
       IF(Found .AND. SSCond > 0.0_dp) Transient = .FALSE.

       IF ( Transient ) THEN
         dt = CurrentModel % Solver % dt
         Order = MIN(CurrentModel % Solver % DoneTime,CurrentModel % Solver % Order)
 
         CALL GetVectorLocalSolution( NodalPVelo, 'Flow Solution', tStep=-1 )
 
         IF ( Order<2 ) THEN
           NodalPVelo(1:dim,1:n)=(NodalVelo(1:dim,1:n)-NodalPVelo(1:dim,1:n))/dt
         ELSE
           CALL GetVectorLocalSolution( Work, 'Flow Solution', tStep=-2 )
           NodalPVelo(1:dim,1:n)=(1.5_dp*NodalVelo(1:dim,1:n) - &
                  2._dp*NodalPVelo(1:dim,1:n)+0.5_dp*Work(1:dim,1:n) )/dt
         END IF
       END IF

       expc = GetCReal(Material,'Heat Expansion Coefficient',Found)
       reft = GetCReal(Material,'Reference Temperature',Found)
       CALL GetConstRealArray( GetConstants(), gWrk, 'Grav',Found )
       IF ( Found ) THEN
         grav = gWrk(1:3,1)*gWrk(4,1)
       ELSE
         grav    =  0.00_dp
         grav(2) = -9.81_dp
       END IF

       LC1 = 2._dp/mK
       LC(1,1:n) = Element % TYPE % NodeU(1:n)
       LC(2,1:n) = Element % TYPE % NodeV(1:n)
       LC(3,1:n) = Element % TYPE % NodeW(1:n)

       DO i=1,Element % TYPE % DIMENSION
         minl=MINVAL(LC(i,1:n))
         maxl=MAXVAL(LC(i,1:n))
         LC(i,1:n) = 2*(LC(i,1:n)-minl)/(maxl-minl)-1
       END DO
     ELSE IF ( Stabilize .AND. Convection ) THEN
       ConvectAndStabilize = .TRUE.
       hK = element % hK * hScale
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

     BodyForce => GetBodyForce()
     FrictionHeat = .FALSE.
     IF (ASSOCIATED(BodyForce))  &
       FrictionHeat = GetLogical( BodyForce, 'Friction Heat', Found )
!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
             Basis,dBasisdx, Bubbles=Bubbles )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       C0 = SUM( NodalC0(1:n) * Basis(1:n) )
       C1 = SUM( NodalC1(1:n) * Basis(1:n) )
       CT = SUM( NodalCT(1:n) * Basis(1:n) )
!------------------------------------------------------------------------------
!       Compute effective heatcapacity, if modelling phase change,
!       at the integration point.
!       NOTE: This is for heat equation only, not generally for diff.conv. equ.
!------------------------------------------------------------------------------
        IF ( PhaseChange ) THEN
          dEnth = 0.0D0
          dTemp = 0.0D0
          DO i=1,3
            dEnth = dEnth + SUM( Enthalpy(1:n) * dBasisdx(1:n,i) )**2
            dTemp = dTemp + SUM( NodalTemperature(1:n) * dBasisdx(1:n,i) )**2
          END DO

          IF( dTemp > TINY( dTemp ) ) THEN
            CL = SQRT( dEnth/dTemp )
          ELSE
            CALL Info('DiffuseConvectiveCompose',&
                'Temperature difference almost zero, cannot account for phase change!',Level=7)
            CL = 0.0
          END IF
          CT = CT + CL
        END IF
!------------------------------------------------------------------------------
!      Coefficient of the diffusion term & it s derivatives at the
!      integration point
!------------------------------------------------------------------------------
       rho = SUM( Nodalrho(1:n) * Basis(1:n) ) 

       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       IF( GotCondModel ) THEN       
         DO i=1,dim
           C2(i,i) = EffectiveConductivity( C2(i,i), rho, Element, &
               NodalTemperature, UX,UY,UZ, Nodes, n, n, u, v, w )
         END DO
       END IF
         
!------------------------------------------------------------------------------
!      If there's no convection term we don't need the velocities, and
!      also no need for stabilization
!------------------------------------------------------------------------------
       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
          Convection = .TRUE.
          IF ( PhaseChange ) C1 = C1 + CL
!------------------------------------------------------------------------------
!         Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
          Velo = 0.0D0
          Velo(1) = SUM( (UX(1:n)-MUX(1:n))*Basis(1:n) )
          Velo(2) = SUM( (UY(1:n)-MUY(1:n))*Basis(1:n) )
          IF ( dim > 2 ) Velo(3) = SUM( (UZ(1:n)-MUZ(1:n))*Basis(1:n) )

          IF ( Compressible ) THEN
            Grad = 0.0D0
            DO i=1,3
              Grad(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
              Grad(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
              IF ( dim > 2 ) Grad(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
            END DO

            Pressure = SUM( NodalPressure(1:n)*Basis(1:n) )
            DivVelo = 0.0D0
            DO i=1,dim
              DivVelo = DivVelo + Grad(i,i)
            END DO
          END IF


          IF ( Vms .OR. VectH ) THEN
            mu = GetCReal( Material, 'Viscosity', Found )
            mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
                   Element, Nodes, n, n, u,v,w,LocalIP=t )

            Grad = 0.0D0
            DO i=1,3
              Grad(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
              Grad(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
              Grad(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
            END DO
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

            Temperature = SUM(Basis(1:n)*NodalTemperature(1:n))

            DO i=1,dim
              GradP(i) = SUM( NodalPressure(1:n)*dBasisdx(1:n,i) )
            END DO

            Gmat = 0._dp
            Gvec = 0._dp
            DO i=1,dim
              DO j=1,dim
                Gvec(i) = Gvec(i) + SUM(LC(j,1:n)*dBasisdx(1:n,i))
                DO k=1,dim
                  Gmat(i,j) = Gmat(i,j) + SUM(LC(k,1:n)*dBasisdx(1:n,i)) * &
                                          SUM(LC(k,1:n)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO

            IF ( Transient ) THEN
              Tau_M = 1._dp / SQRT( SUM(C1*Velo*MATMUL(Gmat,C1*Velo)) + &
                    C2(1,1)**2 * LC1**2*SUM(Gmat*Gmat)/dim + C1**2*4/dt**2 )
            ELSE
              Tau_M = 1._dp / SQRT( SUM(C1*Velo*MATMUL(Gmat,C1*Velo)) + &
                    C2(1,1)**2 * LC1**2 * SUM(Gmat*Gmat)/dim )
            END IF

            IF(Vms) THEN
              RM = 0._dp
              DO p=1,n
                RM(p) = C0 * Basis(p)
                DO i=1,dim
                  RM(p) = RM(p) + C1 * Velo(i) * dBasisdx(p,i)
                  DO j=1,dim
                    RM(p) = RM(p) - C2(i,j)*SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                  END DO
                END DO
              END DO

              VRM = 0._dp
              DO i=1,dim
                VRM(i) = SUM(NodalPVelo(i,1:n)*Basis(1:n))
                DO j=1,dim
                  VRM(i) = VRM(i) + Velo(j) * Grad(i,j)
                  VRM(i) = VRM(i) - (mu/rho)*SUM(GradNodal(1:n,i,j)*dBasisdx(1:n,j))
                END DO
                VRM(i) = VRM(i) + GradP(i)
                VRM(i) = VRM(i) + Grav(i)*ExpC*( Temperature - RefT )
              END DO
            END IF
          END IF

          IF ( .NOT.Vms.AND.Stabilize ) THEN
!------------------------------------------------------------------------------
!           Stabilization parameter Tau
!------------------------------------------------------------------------------
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

#if 1
            Pe  = MIN( 1.0D0, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )
            Tau = 0.0D0
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF

            IF (VectH) Tau = 2*Tau_M
#else
            C00 = C0
            IF ( DT /= 0.0d0 ) C00 = C0+CT/DT

            Pe1 = 0.0d0
            IF ( C00 /= 0.0d0 ) THEN
              Pe1 = 2*ABS(C2(1,1)) / (mK*C00*hK**2)
              Pe1 = C00 * hK**2 * MAX( 1.0d0, Pe1 )
            ELSE
              Pe1 = 2 * ABS(C2(1,1)) / mK
            END IF

            Pe2 = 0.0d0
            IF ( C2(1,1) /= 0.0d0 ) THEN
              Pe2 = ( mK * C1 * VNorm * hK ) / ABS(C2(1,1))
              Pe2 = 2 * ABS(C2(1,1)) * MAX( 1.0d0, Pe2 ) / mK
              Pe  = MIN( 1.0D0, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )
            ELSE
              Pe2 = 2 * hK * C1 * VNorm
            END IF

            Tau = hk**2 / ( Pe1 + Pe2 )
#endif
!------------------------------------------------------------------------------

            DO i=1,dim
              DO j=1,dim
                DO k=1,dim
                  dC2dx(i,j,k) = SUM( NodalC2(i,j,1:n)*dBasisdx(1:n,k) )
                END DO
              END DO
            END DO

!------------------------------------------------------------------------------
!           Compute residual & stablization vectors
!------------------------------------------------------------------------------
            DO p=1,N
              SU(p) = C0 * Basis(p)
              DO i = 1,dim
                SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SU(p) = SU(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO

              SW(p) = C0 * Basis(p)
              DO i = 1,dim
                SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SW(p) = SW(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO
          END IF
        END IF

!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
        DO p=1,NBasis
        DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = CT * Basis(q) * Basis(p)
          A = C0 * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          DO i=1,dim
            DO j=1,dim
              A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
            END DO
          END DO

          IF ( Convection ) THEN
!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
            DO i=1,dim
              A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO
!------------------------------------------------------------------------------
!           Next we add the stabilization...
!------------------------------------------------------------------------------
            IF ( Vms ) THEN
              DO i=1,dim
                A = A - C1 * Tau_M * VRM(i) * dBasisdx(q,i) * Basis(p)

                A = A + C1 * Velo(i) * Tau * RM(q) * dBasisdx(p,i)
                M = M + C1 * Velo(i) * Tau * CT * Basis(q) * dBasisdx(p,i)

                A = A - C1 * Tau_M * VRM(i) * Tau * RM(q) * dBasisdx(p,i)
                M = M - C1 * Tau_M * VRM(i) * Tau * CT * Basis(q) * dBasisdx(p,i)
              END DO
            ELSE IF ( Stabilize ) THEN
              A = A + Tau * SU(q) * SW(p)
              M = M + Tau * CT * Basis(q) * SW(p)
            END IF
          END IF

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
          MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
        END DO
        END DO

!------------------------------------------------------------------------------
!       The righthand side...
!------------------------------------------------------------------------------
!       Force at the integration point
!------------------------------------------------------------------------------
        Force = SUM( LoadVector(1:n)*Basis(1:n) ) + &
          JouleHeat( Element, Nodes, u, v, w, n )

        IF ( Convection ) THEN
!         IF ( Compressible ) Force = Force - Pressure * DivVelo

          Pcoeff = SUM(NodalPressureCoeff(1:n)*Basis(1:n))
          IF ( Pcoeff /= 0.0_dp ) THEN
            Force = Force + Pcoeff * SUM(NodalDPressureDt(1:n)*Basis(1:n))
            DO i=1,dim
              Force = Force + Pcoeff*Velo(i)*SUM(NodalPressure(1:n)*dBasisdx(1:n,i))
            END DO
          END IF

          IF ( FrictionHeat ) THEN
            mu = SUM( Nodalmu(1:n) * Basis(1:n) )
            mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
                   Element, Nodes, n, n, u,v,w )
            IF ( mu > 0.0d0 ) THEN
              IF ( .NOT.Compressible ) THEN
                Grad = 0.0D0
                DO i=1,3
                  Grad(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
                  Grad(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
                  IF ( dim > 2 ) Grad(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
                END DO
              END IF
              Force = Force + 0.5d0 * mu*SecondInvariant(Velo,Grad)
            END IF
          END IF
        END IF
!------------------------------------------------------------------------------
        DO p=1,NBasis
          Load = Force * Basis(p)
          IF ( Vms ) THEN
            DO i=1,dim
              Load = Load + C1 * Velo(i) * Tau  * Force * dBasisdx(p,i)
              Load = Load - C1 * Tau_M * VRM(i) * Tau * Force * dBasisdx(p,i)
            END DO
          ELSE IF ( ConvectAndStabilize ) THEN
            Load = Load + Tau * Force * SW(p)
          END IF
          ForceVector(p) = ForceVector(p) + s * Load
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE DiffuseConvectiveCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for boundary conditions
!>  of diffusion convection equation: 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveBoundary( BoundaryMatrix,BoundaryVector, &
               LoadVector,NodalAlpha,OpenBC,NodalCond,NodalExt,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: coefficient matrix if equations
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT: coefficient of the force term
!
!  REAL(KIND=dp) :: NodalAlpha
!     INPUT: coefficient for temperature dependent term
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!   INTEGER :: n
!       INPUT: Number  of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LoadVector(:),NodalAlpha(:),NodalCond(:),NodalExt(:)
     LOGICAL :: OpenBC
     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t), POINTER :: Element

     INTEGER :: n

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ

     REAL(KIND=dp) :: U,V,W,S
     REAL(KIND=dp) :: Force,Alpha,Normal(3),Coord(3),x,y,z
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       U = U_Integ(t)
       V = V_Integ(t)
       W = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                  Basis,dBasisdx )

       S = detJ * S_Integ(t)
!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis )
       Alpha = SUM( NodalAlpha(1:n)*Basis )

       IF( OpenBC ) THEN
         x = SUM( Nodes % x(1:n)*Basis(1:n) )
         y = SUM( Nodes % y(1:n)*Basis(1:n) )
         z = SUM( Nodes % z(1:n)*Basis(1:n) )         

         Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
         Coord(1) = x
         Coord(2) = y
         Coord(3) = z
         Alpha = SUM( Basis(1:n) * NodalCond(1:n) ) *  &
             SUM( Coord * Normal ) / SUM( Coord * Coord ) 
         Force = Alpha * SUM( Basis(1:n) * NodalExt(1:n) )
       END IF

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
              s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
   END SUBROUTINE DiffuseConvectiveBoundary
!------------------------------------------------------------------------------

END MODULE DiffuseConvective

!> \}
