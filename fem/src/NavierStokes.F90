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
!>  Module computing Navier-Stokes local matrices (cartesian coordinates)
!------------------------------------------------------------------------------

MODULE NavierStokes

  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE ElementDescription!, ONLY: GetEdgeMap

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for Navier-Stokes-Equations
!> in cartesian coordinates.
!------------------------------------------------------------------------------
   SUBROUTINE NavierStokesCompose  (                                            &
       MassMatrix, StiffMatrix, ForceVector, LoadVector, Nodalmu,        &
       Nodalrho, Ux, Uy, Uz, MUx, MUy, MUz, NodalPressure, NodalTemperature,&
       Convect, StabilizeFlag, Cmodel, &
       PseudoCompressible, NodalCompressibility, NodalGasC, Porous, NodalDrag,  &
       PotentialForce, PotentialField, PotentialCoefficient, MagneticForce,     &
       Rotating, Omega, &
       DivDiscretization,  gradPDiscretization, NewtonLinearization,            &
       Element, n, Nodes )
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
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values for viscosity (i.e. if turbulence model or
!             power-law viscosity is used, the values vary in space)
!
!  REAL(KIND=dp) :: Nodalrho(:)
!     INPUT: Nodal values of density
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!
!  REAL(KIND=dp) :: NodalPressure(:)
!     INPUT: Nodal values of total pressure from previous iteration
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilization be used ?
!
!  LOGICAL :: Compressible1, Compressible2
!     INPUT: Should compressible flow terms be added ?
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
!  LOGICAL :: Rotating
!      INPUT: Is the coordinate system rotating
!  
!  REAL(KIND=dp) :: Omega(:)
!      INPUT: If previous is True, components of angular velocity
!
!  LOGICAL :: NewtonLinearization
!      INPUT: Picard or Newton linearization of the convetion term ?
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
     REAL(KIND=dp), DIMENSION(:) :: Ux,Uy,Uz,MUx,MUy,MUz,Omega
     LOGICAL :: PseudoCompressible, Porous, &
         NewtonLinearization, Convect, DivDiscretization, gradPDiscretization, &
         PotentialForce, MagneticForce, Rotating
     CHARACTER(LEN=*) :: StabilizeFlag
     REAL(KIND=dp) :: Nodalmu(:),Nodalrho(:), &
       NodalPressure(:), LoadVector(:,:), NodalTemperature(:), &
       NodalCompressibility(:), NodalDrag(:,:), PotentialField(:), &
       PotentialCoefficient(:), NodalGasC(:)

     INTEGER :: n, Cmodel

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

     REAL(KIND=dp) :: LES(6,6,n)
     TYPE(Variable_t), POINTER :: LESVar
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp),target :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp),target :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: Velo(3),UVelo(3),Grad(3,3),Force(4),Metric(3,3),Symb(3,3,3), Drag(3), Coord(3)

     REAL(KIND=dp), POINTER :: A(:,:),M(:,:),Load(:), Jac(:,:)
     REAL(KIND=dp) :: SU(n,4,4),SW(n,4,4),LrF(3), LSx(3,3), VC(6,6), B(6,3), G(3,6)

     REAL(KIND=dp) :: Lambda=1.0,Re,Tau,Delta,Re1,Re2,GasC, hScale

     REAL(KIND=dp) :: VNorm,hK,mK,mu,dmudx(3),Temperature, &
               drhodp,drhodp_n(n)

     REAL(KIND=dp) :: drhodx(3), rho, Pressure, dTemperaturedx(3), &
                      dPressuredx(3),dPrevPressuredx(3), Compress, masscoeff

     REAL(KIND=dp), TARGET :: JacM(8*n,8*n), SOL(8*n)

     REAL(KIND=dp) :: Tau_M, Tau_C, Gmat(3,3),Gvec(3), dt=0._dp, C1, NodalVelo(4,n), &
       RM(n,4,4),LC(3,n), gradP(n), PRM(4), GradNodal(n,4,3), PVelo(3), NodalPVelo(4,n), &
       Work(3,n), maxl, minl, Strain(3,3), muder0, muder, ViscConstantCondition

     INTEGER :: i,j,k,l,c,p,q,t,dim

     REAL(KIND=dp) :: s,u,v,w,volume

     INTEGER, POINTER :: EdgeMap(:,:)

     TYPE(ElementType_t), POINTER :: LinearType, SaveType
     INTEGER :: LinearBasis, LinearCode(3:8) = (/ 303,404,504,605,706,808 /)
  
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis, deg(100), Order
     REAL(KIND=dp), POINTER :: NodalC(:,:,:)
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, Transient, stat, Bubbles, PBubbles, Stabilize, &
                VMS, P2P1, Isotropic, drhodp_found, Compressible, ViscNewtonLin, &
                ViscNonnewtonian, LaplaceDiscretization,OutOfPlaneFlow

     ! local transposed values 
     REAL(KIND=dp) :: dBasisdxp(n*2)
     REAL(KIND=dp) :: muderq(n*2)
     REAL(KIND=dp),target :: StiffMatrixTrabsp(n*2*4, n*2*4)
     REAL(KIND=dp),target :: MassMatrixTrabsp(n*2*4, n*2*4)
     REAL(KIND=dp) :: JacMTrabsp(n*2*4, n*2*4)
     REAL(KIND=dp),target :: BasePVec(2*n) ! todo, change this and the assorted behavior to a pointer
     REAL(KIND=dp),POINTER :: BasePPtr(:)
     REAL(KIND=dp),POINTER :: dBasisdxPtrP(:),dBasisdxPtrQ(:)
     REAL(KIND=dp) :: DivLapConst, ComprConvConst, ViscNewtonLinConst
     REAL(KIND=dp),POINTER :: gradPDiscPtrP(:), gradPDiscPtrQ(:)
     REAL(KIND=dp) :: gradPDiscConst, CmodelConst
     REAL(KIND=dp) ::SWxSU(n), StrainS(2*n, 4)
     INTEGER :: P2P1stop, ii, jj, kk, ll
     ! end of transposed varible 

     
!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     hScale = GetCReal( Params, 'H scale', Found )
     IF ( .NOT. Found )  hScale = 1._dp

     OutOfPlaneFlow = GetLogical( Params, 'Out of Plane flow', Found)
     IF ( .NOT. Found ) OutOfPlaneFlow = .FALSE.

!#ifdef LES
!     LESVar => VariableGet( CurrentModel % Variables, 'LES' )
!     IF ( ASSOCIATED( LESVar ) ) THEN
!        k = 0
!        DO i=1,dim
!          DO j=i,dim
!             k = k + 1
!             LES(i,j,1:n) = LESVar % Values(LESVar % DOFs*(LESVar % Perm( &
!                         Element % NodeIndexes)-1)+k)
!          END DO
!        END DO
!        DO i=1,dim
!        DO j=1,i-1
!           LES(i,j,1:n) = LES(j,i,1:n)
!        END DO
!        END DO
!     END IF
!#endif

     c = dim + 1

     ForceVector = 0._dp
     MassMatrix  = 0._dp
     StiffMatrix = 0._dp

     StiffMatrixTrabsp = 0._dp
     MassMatrixTrabsp = 0._dp
     JacMTrabsp = 0._dp


     IF ( NewtonLinearization ) JacM = 0._dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis    = n
     Bubbles   = .FALSE.
     PBubbles  = .FALSE.
     P2P1 = .FALSE.
     VMS =  StabilizeFlag == 'vms'

     Material => GetMaterial()
     Isotropic=GetString( Material, 'Viscosity Model',ViscNonnewtonian )/='anisotropic'
     IF ( .NOT. Isotropic ) THEN
       CALL ListGetRealArray( Material, 'Viscosity', &
            NodalC,n,Element % NodeIndexes,Found )
     END IF

     IF( ViscNonnewtonian ) THEN
       ViscConstantCondition = GetCReal( Params, 'Newtonian Viscosity Condition',Found)
       IF( Found .AND. ViscConstantCondition > 0.0_dp ) ViscNonnewtonian = .FALSE.
     END IF

     LaplaceDiscretization = GetLogical( Params,'Laplace Discretization', Found)

     Transient = GetString(GetSimulation(),'Simulation Type',Found)=='transient'
     IF ( VMS .AND. Transient ) THEN
       dt = CurrentModel % Solver % dt
       Order = MIN(CurrentModel % Solver % DoneTime, CurrentModel % Solver % Order)
     END IF

     Compressible = Cmodel /= Incompressible
     drhodp_found = .FALSE.
     drhodp_n(1:n) = GetReal( Material, 'drho/dp', drhodp_found )

     Stabilize = StabilizeFlag == 'stabilized'
     IF ( .NOT.(Vms .OR. Stabilize) .OR. Compressible ) THEN
       PBubbles  = isActivePElement(Element) .AND. Element % BDOFs > 0
       IF ( PBubbles ) THEN
          NBasis = n + Element % BDOFs
       ELSE IF ( StabilizeFlag == 'bubbles' .OR. &
              Element % TYPE % BasisFunctionDegree<=1 ) THEN
          NBasis    = 2 * n
          Bubbles   = .TRUE.
       ELSE
          P2P1 = StabilizeFlag == 'p2/p1' .OR. StabilizeFlag == 'p2p1'
       END IF
       Stabilize = .FALSE.
     END IF

     LinearBasis = 0
     IF ( P2P1 ) THEN
        SaveType => Element % TYPE
        j = GetElementFamily()
        LinearType => GetElementType(LinearCode(j))
        LinearBasis = LinearType % NumberOFNodes
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
     N_Integ =  IntegStuff % n
!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!------------------------------------------------------------------------------
    
     hK = element % hK*hscale
     mK = element % StabilizationMK

     IF ( Stabilize ) THEN
       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % TYPE % NodeU(p)
         v = Element % TYPE % NodeV(p)
         w = Element % TYPE % NodeW(p)
         stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
     ELSE IF  ( Vms ) THEN
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
      IF ( Transient ) THEN
        CALL GetVectorLocalSolution( NodalPVelo, tStep=-1 )

        IF ( Order<2 ) THEN
          NodalPVelo(1:dim,1:n)=(NodalVelo(1:dim,1:n)-NodalPVelo(1:dim,1:n))/dt
        ELSE
          CALL GetVectorLocalSolution( Work, tStep=-2 )
          NodalPVelo(1:dim,1:n)=(1.5_dp*NodalVelo(1:dim,1:n)- &
                 2._dp*NodalPVelo(1:dim,1:n)+0.5_dp*Work(1:dim,1:n) )/dt
        END IF
      END IF

      C1 = 2._dp/mK
      LC(1,1:n) = Element % TYPE % NodeU(1:n)
      LC(2,1:n) = Element % TYPE % NodeV(1:n)
      LC(3,1:n) = Element % TYPE % NodeW(1:n)

      DO i=1,Element % TYPE % DIMENSION
        minl=MINVAL(LC(i,1:n))
        maxl=MAXVAL(LC(i,1:n))
        LC(i,1:n) = 2*(LC(i,1:n)-minl)/(maxl-minl)-1
      END DO
     END IF

     ViscNewtonLin = .FALSE.

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
      IF ( P2P1 ) THEN
         Element % TYPE => LinearType
         stat = ElementInfo( Element, Nodes, u, v, w, detJ,pBasis,pdBasisdx )
         Element % TYPE => SaveType
      END IF

      stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
              Basis, dBasisdx, Bubbles=Bubbles )

      s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!     rho at the integration point
!------------------------------------------------------------------------------
      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
      dPressuredx = 0._dp
      dTemperaturedx = 0._dp
      SELECT CASE(Cmodel)
      CASE(PerfectGas1)
        IF ( P2P1 ) THEN
          k = LinearBasis
          Pressure    = SUM( NodalPressure(1:k) * pBasis(1:k) )
          Temperature = SUM( NodalTemperature(1:k) * pBasis(1:k) )
          DO i=1,dim
            dPressuredx(i)  = SUM( NodalPressure(1:k)  * pdBasisdx(1:k,i) )
            dTemperaturedx(i)  = SUM( NodalTemperature(1:k)  * pdBasisdx(1:k,i) )
          END DO
        ELSE
          Pressure    = SUM( NodalPressure(1:n) * Basis(1:n) )
          Temperature = SUM( NodalTemperature(1:n) * Basis(1:n) )
          DO i=1,dim
            dPressuredx(i)  = SUM( NodalPressure(1:n)  * dBasisdx(1:n,i) )
            dTemperaturedx(i)  = SUM( NodalTemperature(1:n)  * dBasisdx(1:n,i) )
          END DO
        END IF

        GasC = SUM( Basis(1:n) * NodalGasC(1:n) )
        rho = Pressure / ( GasC*Temperature )
      CASE (UserDefined1,Thermal)
        DO i=1,dim
          drhodx(i)  = SUM( Nodalrho(1:n) * dBasisdx(1:n,i) )
        END DO
      CASE(UserDefined2)
        DO i=1,dim
          dPressuredx(i) = SUM( NodalPressure(1:n) * dBasisdx(1:n,i) )
        END DO
        drhodp = SUM( drhodp_n(1:n) * Basis(1:n) )
      END SELECT

      IF ( PseudoCompressible ) THEN
        Pressure = SUM( NodalPressure(1:n) * Basis(1:n) )        
        Compress = rho * SUM(NodalCompressibility(1:n)*Basis(1:n))      
      END IF

!------------------------------------------------------------------------------
!     Velocity from previous iteration (relative to mesh velocity)
!     at the integration point
!------------------------------------------------------------------------------
      Velo = 0.0_dp
      Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
      Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
      IF ( DIM > 2 .OR. OutOfPlaneFlow ) Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )

      Grad = 0.0_dp
      DO i=1,3
        Grad(1,i) = SUM( Ux(1:n) * dBasisdx(1:n,i) )
        Grad(2,i) = SUM( Uy(1:n) * dBasisdx(1:n,i) )
        IF ( DIM > 2 .OR. OutOfPlaneFlow ) Grad(3,i) = SUM( Uz(1:n) * dBasisdx(1:n,i) )
      END DO
!------------------------------------------------------------------------------
!     Force at integration point
!------------------------------------------------------------------------------
      Force = 0.0_dp
      DO i=1,c
        Force(i) = SUM( LoadVector(i,1:n) * Basis(1:n) )
      END DO

      ! The relative change in temperature is the source term 
      ! for continuity equation.
      !------------------------------------------------------
      IF ( Compressible .AND. Cmodel==PerfectGas1 ) THEN
        Force(c) = Force(c) / Temperature
      END IF


      IF ( MagneticForce ) THEN
         LrF = LorentzForce( Element,Nodes,u,v,w,n )
         Force(1:DIM) = Force(1:DIM) + Lrf(1:DIM) / rho
      END IF

      IF( Rotating ) THEN
        Coord(1) = SUM( Basis(1:n) * Nodes % x(1:n) )
        Coord(2) = SUM( Basis(1:n) * Nodes % y(1:n) )
        Coord(3) = SUM( Basis(1:n) * Nodes % z(1:n) )
        
        ! langranges formula is used to simplify the triple product
        ! omega x ( omega x coord ) = omega(omega.coord) - coord(omega.omega)

        ! This will be multiplied by density later on
        Force(1:dim) = Force(1:dim) - Omega(1:dim) * SUM( Omega * Coord )
        Force(1:dim) = Force(1:dim) + Coord(1:dim) * SUM( Omega * Omega )
      END IF


!------------------------------------------------------------------------------
!     Additional forces due to gradient forces (electrokinetic flow) and
!     viscous drag in porous media.
!------------------------------------------------------------------------------

      IF(PotentialForce) THEN
        DO i=1,DIM
          Force(i) = Force(i) - SUM( PotentialCoefficient(1:n) * Basis(1:n) ) * &
              SUM(  PotentialField(1:n) * dBasisdx(1:n,i) )
        END DO
      END IF

      IF(Porous) THEN
        DO i=1,DIM
          Drag(i) = SUM( NodalDrag(i,1:n) * Basis(1:n) )
        END DO
      END IF

      IF ( Convect .AND. NewtonLinearization ) THEN
        Uvelo = 0.0_dp
        Uvelo(1) = SUM( Basis(1:n) * Ux(1:n) )
        Uvelo(2) = SUM( Basis(1:n) * Uy(1:n) )
        IF ( DIM > 2 .OR. OutOfPlaneFlow) Uvelo(3) = SUM( Basis(1:n) * Uz(1:n) )

        DO i=1,dim
          DO j=1,dim
            Force(i) = Force(i) + Grad(i,j) * Uvelo(j) 
          END DO
        END DO
      END IF

!------------------------------------------------------------------------------
!     Effective viscosity & derivatives at integration point
!------------------------------------------------------------------------------
      IF ( isotropic ) THEN
        mu = SUM( Nodalmu(1:n) * Basis(1:n) )

        IF( ViscNonnewtonian ) THEN
          mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
              Element, Nodes, n, n, u, v, w,  muder0, LocalIP=t )

          ViscNewtonLin = NewtonLinearization .AND. muder0/=0._dp
          IF ( ViscNewtonLin )  Strain = (Grad+TRANSPOSE(Grad))/2
        END IF

      ELSE
        DO i=1,6
        DO j=1,6
          VC(i,j) = SUM(NodalC(i,j,1:n)*Basis(1:n))
        END DO
        END DO
      END IF


      IF ( Stabilize ) THEN
        DO i=1,3
          dmudx(i) = SUM( Nodalmu(1:n)*dBasisdx(1:n,i) )
        END DO
!------------------------------------------------------------------------------
!       Stabilization parameters Tau & Delta
!------------------------------------------------------------------------------
        IF ( Convect ) THEN
           VNorm = MAX( SQRT( SUM(Velo(1:DIM)**2) ), 1.0d-12 )
           Re = MIN( 1.0d0, rho * mK * hK * VNorm / (4 * mu) )
 
           Tau = hK * Re / (2 * rho * VNorm)
           Delta = rho * Lambda * Re * hK * VNorm
        ELSE
           Delta = 0._dp
           Tau   = mK * hK**2  / ( 8 * mu )
        END IF
!------------------------------------------------------------------------------
!       SU will contain residual of ns-equations (except for the time derivative
!       and force terms). SW will contain the weight function values.
!------------------------------------------------------------------------------
        SU(1:n,:,:) = 0.0D0
        SW(1:n,:,:) = 0.0D0

        DO p=1,N
          DO i=1,DIM
            SU(p,i,c) = SU(p,i,c) + dBasisdx(p,i)
            IF(Porous) THEN
              SU(p,i,i) = SU(p,i,i) + mu * Drag(i) * Basis(p)
            END IF

            IF ( Convect ) THEN
              DO j=1,DIM
                SU(p,i,i) = SU(p,i,i) + rho * dBasisdx(p,j) * Velo(j)
              END DO
            END IF

            DO j=1,DIM
              SU(p,i,i) = SU(p,i,i) - dmudx(j) * dBasisdx(p,j)
              SU(p,i,j) = SU(p,i,j) - dmudx(j) * dBasisdx(p,i)
              SU(p,i,i) = SU(p,i,i) - mu * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
              SU(p,i,j) = SU(p,i,j) - mu * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
            END DO

            IF ( Convect .AND. NewtonLinearization ) THEN
              DO j=1,DIM
                SU(p,i,j) = SU(p,i,j) + rho * Grad(i,j) * Basis(p)
              END DO
            END IF
!
!------------------------------------------------------------------------------

            IF ( Convect ) THEN
              SW(p,c,i) = SW(p,c,i) + rho * dBasisdx(p,i)
              DO j=1,dim
                SW(p,i,i) = SW(p,i,i) + rho * dBasisdx(p,j) * Velo(j)
              END DO
            ELSE
              SW(p,c,i) = SW(p,c,i) + dBasisdx(p,i)
            END IF

            DO j=1,dim
              SW(p,i,i) = SW(p,i,i) - dmudx(j) * dBasisdx(p,j)
              SW(p,j,i) = SW(p,j,i) - dmudx(j) * dBasisdx(p,i)
              SW(p,i,i) = SW(p,i,i) - mu * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
              SW(p,j,i) = SW(p,j,i) - mu * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
            END DO
          END DO
        END DO
      ELSE IF ( Vms ) THEN

        mu = mu / rho
        rho = 1._dp

        DO i=1,dim
          PVelo(i) = rho*SUM(NodalPVelo(i,1:n)*Basis(1:n))
          gradP(i) = SUM(NodalVelo(dim+1,1:n)*dBasisdx(1:n,i))
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
          Tau_M = 1._dp / SQRT( SUM(Velo*MATMUL(Gmat,Velo)) + &
            C1**2*(mu/rho)**2*SUM(Gmat*Gmat)/dim + 4/dt**2 )
        ELSE
          Tau_M = 1._dp / SQRT( SUM(Velo*MATMUL(Gmat,Velo)) + &
             C1**2 * (mu/rho)**2 * SUM(Gmat*Gmat)/dim )
        END IF
        Tau_C = 1._dp / (Tau_M*SUM(Gvec*Gvec))

!
!       VNorm = MAX( SQRT( SUM(Velo(1:dim)**2) ), 1.0d-12 )
!       Re = MIN( 1.0d0, rho * mK * hK * VNorm / (4 * mu) )
!
!       Tau = hK * Re / (2 * rho * VNorm)
!       Delta = rho * Lambda * Re * hK * VNorm
!PRINT*, hk,mk,2/c1,2._dp/sqrt(Gmat(1,1)),2._dp/sqrt(Gmat(2,2))
!PRINT*,tau,tau_m,dim
!PRINT*,'---'
!

        RM = 0._dp
        DO p=1,n
          DO i=1,dim
            DO j=1,dim
              RM(p,i,i) = RM(p,i,i) + rho*Velo(j)*dBasisdx(p,j)
              RM(p,i,i) = RM(p,i,i) - mu*SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
            END DO
            RM(p,i,dim+1) = RM(p,i,dim+1) + dBasisdx(p,i)
            RM(p,dim+1,i) = RM(p,dim+1,i) + rho*dBasisdx(p,i)
          END DO
        END DO

        PRM = 0._dp
        DO i=1,dim
          DO j=1,dim
            PRM(i) = PRM(i) + rho*Velo(j)*Grad(i,j)
            PRM(i) = PRM(i) - mu*SUM(GradNodal(1:n,i,j)*dBasisdx(1:n,j))
          END DO
          PRM(i) = PRM(i) + GradP(i)
          PRM(dim+1) = PRM(dim+1) + rho*Grad(i,i)
        END DO
      END IF


!#ifdef LES
!#define LSDelta (SQRT(2.0d0)*Element % hK)
!#define LSGamma 6

!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!    Vectroized first 
!------------------------------------------------------------------------------
!
    IF (ViscNewtonLin) THEN
      DO j=1,dim
        !$omp simd
        DO p=1,NBasis
          StrainS(p, j) = SUM(Strain(j,:)*dBasisdx(p,:))
        END DO  
      END DO  
      DO i=1,dim
        DO j=1,dim
          !DIR$ Unroll(2)
          DO q=1,NBasis
          !$omp simd
            DO p=1,NBasis
              !Jac(j,i)=Jac(j,i)
              JacMTrabsp(((j-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = JacMTrabsp(((j-1)*(NBasis)) + (p),((i-1)*(NBasis)) & 
              + (q)) + s * 2 * muder0 * 4 * StrainS(q, i) * StrainS(p, j)
            END DO
          END DO
        END DO
      END DO
    END If

    IF ( P2P1 ) THEN
      DO p = 1, NBasis
        IF (p<=LinearBasis ) THEN
          BasePVec(p) = PBasis(p)
        ELSE
          BasePVec(p) =  Basis(p)
        END IF
      END DO
    ELSE
        BasePVec =  Basis
    END IF


    !c = dim + 1
    IF ( Isotropic ) THEN
      DO i=1,dim
        DO j = 1,dim
          IF ( divDiscretization ) THEN
            dBasisdxPtrQ => dBasisdx(:,j)
            dBasisdxPtrP => dBasisdx(:,i)
          ELSE IF (.NOT.LaplaceDiscretization) THEN
            dBasisdxPtrQ => dBasisdx(:,i)
            dBasisdxPtrP => dBasisdx(:,j)
          END IF
          !DIR$ Unroll(2)
          DO q=1,NBasis           
            !$omp simd
            DO p=1,NBasis
              !A(i,i) = A(i,i)
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((i-1)*(NBasis)) + (q)) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
            END DO ! p nbasis simd

            IF ( divDiscretization .or. .NOT.LaplaceDiscretization) THEN 
              !$omp simd
              DO p=1,NBasis
                !A(i,j) = A(i,j)
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                  ((j-1)*(NBasis)) + (q)) + s * mu * dBasisdxPtrQ(q) * dBasisdxPtrP(p)
              END DO ! p nbasis simd
            END IF

            IF ( Compressible ) THEN
              !$omp simd
              DO p=1,NBasis
                !A(i,j) = A(i,j) 
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                  ((j-1)*(NBasis)) + (q)) - s * ( 2._dp / 3._dp ) * mu * dBasisdx(q,j) * dBasisdx(p,i)
              END DO ! p nbasis simd
            END IF ! compressible 
          END DO ! q basis
        END DO ! j dim
      END DO ! i dim
    END IF ! isotropic

    IF ( Convect ) THEN
      DO i=1,dim
        DO j = 1,dim
          !DIR$ Unroll(2)
          DO q=1,NBasis
            !$omp simd
            DO p=1,NBasis
              !A(i,i) = A(i,i)
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((i-1)*(NBasis)) + (q)) + s * rho * dBasisdx(q,j) * Velo(j) * Basis(p)
            END DO ! p nbasis simd

!------------------------------------------------------------------------------
!      Convection, Newton linearization
!------------------------------------------------------------------------------
            IF (NewtonLinearization ) THEN
              !$omp simd
              DO p=1,NBasis
                !A(i,j) = A(i,j) 
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                  ((j-1)*(NBasis)) + (q)) + s * rho * Grad(i,j) * Basis(q) * Basis(p)
              END DO ! p nbasis simd
            END IF ! NewtonLinearization
          END DO ! q basis
        END DO ! j dim
      END DO ! i dim
    END IF ! convect


    DO i=1,dim
      IF (P2P1) THEN ! this could be done with pointers 
        P2P1stop = LinearBasis
        IF ( gradPDiscretization  ) THEN
          gradPDiscPtrQ => PdBasisdx(:,i)
          gradPDiscPtrP => Basis(:)
          gradPDiscConst = s
        ELSE
          gradPDiscPtrQ => PBasis(:)
          gradPDiscPtrP => dBasisdx(:,i)
          gradPDiscConst = -s
        END IF
      ELSE
        P2P1stop = NBasis

        IF ( gradPDiscretization  ) THEN
          gradPDiscPtrQ => dBasisdx(:,i)
          gradPDiscPtrP => Basis(:)
          gradPDiscConst = s
        ELSE
          gradPDiscPtrQ => Basis(:)
          gradPDiscPtrP => dBasisdx(:,i)
          gradPDiscConst = -s
        END IF
      END IF ! P2P1

       !DIR$ Unroll(2)
       DO q=1,P2P1stop
         !$omp simd
         DO p=1,NBasis
           !A(i,c) = A(i,c)
           StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((c-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
            ((c-1)*(NBasis)) + (q)) + gradPDiscConst * gradPDiscPtrQ(q) * gradPDiscPtrP(p)
         END DO ! p nbasis simd
       END DO ! q NBasis
    END DO ! i dim

!
    IF ( Compressible .OR. Convect ) THEN
      ComprConvConst = rho
    ELSE
      ComprConvConst = 1
    END IF

    DO i=1,dim
      IF ( gradPDiscretization  ) THEN
        gradPDiscPtrQ => Basis(:)
        gradPDiscPtrP => dBasisdx(:,i)
        gradPDiscConst = -s * rho
      ELSE
        gradPDiscPtrQ => dBasisdx(:,i)
        gradPDiscPtrP => BasePVec(:)
        gradPDiscConst = s * ComprConvConst
        END IF
      !DIR$ Unroll(2)
      DO q=1,NBasis
        !$omp simd
        DO p=1,NBasis
          !A(c,i) = A(c,i) &
          StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
            ((i-1)*(NBasis)) + (q)) + gradPDiscConst * gradPDiscPtrQ(q) * gradPDiscPtrP(p)
        END DO ! p nbasis simd
      END DO ! q nbasis
    END DO ! i dim
    
    IF ( .NOT. gradPDiscretization .AND. (Cmodel==PerfectGas1 .OR. Cmodel==UserDefined1 .OR. Cmodel==Thermal &
      .OR. Cmodel==UserDefined2) ) THEN
      DO i=1,dim
        !DIR$ Unroll(2)
        DO q=1,NBasis
          SELECT CASE(Cmodel)
            CASE(PerfectGas1)
              !$omp simd
              DO p=1,NBasis
                !A(c,i) = A(c,i) &
                StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
                  ((i-1)*(NBasis)) + (q)) + s * ( rho / Pressure ) * Basis(q) * dPressuredx(i) * BasePVec(p) / 2
                !A(c,c) = A(c,c) &
                StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((c-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
                  ((c-1)*(NBasis)) + (q)) + s * ( rho / Pressure ) * Velo(i) * dBasisdx(q,i) * BasePVec(p) / 2
                !A(c,i) = A(c,i) &
                StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
                  ((i-1)*(NBasis)) + (q)) - s * ( rho / Temperature ) * Basis(q) * dTemperaturedx(i) *  BasePVec(p)
              END DO ! p nbasis simd
            CASE(UserDefined1, Thermal)
              !$omp simd
              DO p=1,NBasis
                !A(c,i) = A(c,i) &
                StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
                  ((i-1)*(NBasis)) + (q)) + s * drhodx(i) * Basis(q) * BasePVec(p)
              END DO ! p nbasis simd
            CASE(UserDefined2)
              !$omp simd
              DO p=1,NBasis
                !A(c,c) = A(c,c) &
                StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((c-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
                  ((c-1)*(NBasis)) + (q)) + s * drhodp*dBasisdx(q,i)*Velo(i)*BasePVec(p)/2
                !A(c,i) = A(c,i) &
                StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
                  ((i-1)*(NBasis)) + (q)) + s * drhodp*dPressuredx(i)*Basis(q)*BasePVec(p)/2
              END DO ! p nbasis simd
          END SELECT ! Cmodel
        END DO ! q nbasis
      END DO ! i dim
    END IF !gradPDiscretization

!!------------------------------
!! Possible Porous media effects
!!------------------------------------------------------------------------------
    IF(Porous) THEN
      DO i=1,dim
        !DIR$ Unroll(2)
        DO q=1,NBasis
            !$omp simd
            DO p=1,NBasis
              !A(i,i) = A(i,i)
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((i-1)*(NBasis)) + (q)) + s * mu * Drag(i) * Basis(q) * Basis(p)
            END DO ! p nbasis simd
        END DO ! q nbasis
      END DO ! i dim
    END IF
!
!
!
!!------------------------------------------------------------------------------
!!      Mass matrix:
!!------------------------------------------------------------------------------
!
    DO i=1,dim
      !DIR$ Unroll(2)
      DO q=1,NBasis
      !$omp simd
        DO p=1,NBasis
          !M(i,i) = M(i,i)
          MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
            ((i-1)*(NBasis)) + (q)) + s*rho*Basis(q)*Basis(p)
        END DO
      END DO
    END DO ! i dim

    IF ( Cmodel==PerfectGas1 .OR. Cmodel==UserDefined2) THEN 
      IF ( Cmodel==PerfectGas1 ) THEN 
        CmodelConst = ( rho / Pressure )
      ELSE IF ( Cmodel==UserDefined2) THEN
        CmodelConst = drhodp
      END IF ! Cmodel==PerfectGas1

      DO q=1,NBasis
        !$omp simd
        DO p=1,NBasis
          !M(c,c) = M(c,c) &
          MassMatrixTrabsp(((c-1)*(NBasis)) + (p),((c-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((c-1)*(NBasis)) + (p), & 
            ((c-1)*(NBasis)) + (q)) + s * CmodelConst * Basis(q) * BasePVec(p)
        END DO ! p nbasis simd
      END DO

    END IF
    
    IF (PseudoCompressible) THEN 
      DO q=1,NBasis
          !$omp simd
          DO p=1,NBasis
            !A(c,c) = A(c,c) &
            StiffMatrixTrabsp(((c-1)*(NBasis)) + (p),((c-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((c-1)*(NBasis)) + (p), &
              ((c-1)*(NBasis)) + (q)) + s * Compress * Basis(q) * BasePVec(p)
          END DO ! p nbasis simd
      END DO ! q nbasis
    END IF ! PseudoCompressible

    IF( Rotating ) THEN
      DO q=1,NBasis
        !$omp simd
        DO p=1,NBasis
          masscoeff = 2 * s * rho * Basis(q) * Basis(p)
          !A(1,2) = A(1,2)&
          StiffMatrixTrabsp(((1-1)*(NBasis)) + (p),((2-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((1-1)*(NBasis)) + (p), & 
          ((2-1)*(NBasis)) + (q)) - masscoeff * Omega(3)
          !A(2,1) = A(2,1)&
          StiffMatrixTrabsp(((2-1)*(NBasis)) + (p),((1-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((2-1)*(NBasis)) + (p), & 
          ((1-1)*(NBasis)) + (q)) + masscoeff * Omega(3)
        END DO ! p nbasis simd
        IF( dim == 3) THEN
          !$omp simd
          DO p=1,NBasis
            masscoeff = 2 * s * rho * Basis(q) * Basis(p)
            !A(1,3) = A(1,3) &
            StiffMatrixTrabsp(((1-1)*(NBasis)) + (p),((3-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((1-1)*(NBasis)) + (p), & 
              ((3-1)*(NBasis)) + (q)) + masscoeff * Omega(2)          
            !A(2,3) = A(2,3) &
            StiffMatrixTrabsp(((2-1)*(NBasis)) + (p),((3-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((2-1)*(NBasis)) + (p), & 
              ((3-1)*(NBasis)) + (q)) - masscoeff * Omega(1)
            !A(3,2) = A(3,2) &
            StiffMatrixTrabsp(((3-1)*(NBasis)) + (p),((2-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((3-1)*(NBasis)) + (p), & 
              ((2-1)*(NBasis)) + (q)) + masscoeff * Omega(1)
            !A(3,1) = A(3,1) &
            StiffMatrixTrabsp(((1-1)*(NBasis)) + (p),((3-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((1-1)*(NBasis)) + (p), & 
              ((3-1)*(NBasis)) + (q)) - masscoeff * Omega(2)
          END DO ! p nbasis simd
        END IF
      END DO ! q nbasis
    END IF ! Rotating

!!------------------------------------------------------------------------------
!!    Loop over basis functions (of both unknowns and weights)
!!    Then the stuff i havent vectorized yet
!!------------------------------------------------------------------------------
    
!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
    IF (.NOT.Isotropic) THEN
      !
      DO p=1,NBasis
        G = 0.0d0
        G(1,1) = dBasisdx(p,1)
        G(2,2) = dBasisdx(p,2)
        G(3,3) = dBasisdx(p,3)
        G(1,4) = dBasisdx(p,2)
        G(2,4) = dBasisdx(p,1)
        G(2,5) = dBasisdx(p,3)
        G(3,5) = dBasisdx(p,2)
        G(1,6) = dBasisdx(p,3)
        G(3,6) = dBasisdx(p,1)
        G = MATMUL( G, VC )
        DO q=1,NBasis
          A => StiffMatrixTrabsp( p:NBasis*c: NBasis, q:NBasis*c: NBasis )

          B = 0.0d0
          B(1,1) = dBasisdx(q,1)
          B(2,2) = dBasisdx(q,2)
          B(3,3) = dBasisdx(q,3)
          B(4,1) = dBasisdx(q,2)
          B(4,2) = dBasisdx(q,1)
          B(5,2) = dBasisdx(q,3)
          B(5,3) = dBasisdx(q,2)
          B(6,1) = dBasisdx(q,3)
          B(6,3) = dBasisdx(q,1)
          A(1:3,1:3) = A(1:3,1:3) + s * MATMUL( G, B )
        END DO
      END DO
    END IF ! .NOT.Isotropic
    
            
    IF ( Stabilize ) THEN 
      
      DO q=1,NBasis 
        DO j=1,c
          DO i=1,c
            SWxSU = 0.0
            DO k=1,dim 
              !$omp simd 
              DO p=1,NBasis
                SWxSU(p) = SWxSU(p) + (SW(p, i, k) * SU(q, k, j))
              END DO
            END DO 
            !$omp simd 
            DO p=1,NBasis
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) + s*Tau* SWxSU(p)
            END DO ! p nbasis simd
          END DO ! i c
          DO i=1,dim
            !$omp simd
            DO p=1,NBasis
              !M => MassMatrixTrabsp ( p:NBasis*c: NBasis, q:NBasis*c: NBasis )
              !M(j,i) = M(j,i) &
              MassMatrixTrabsp(((j-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((j-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) + s * Tau * rho * Basis(q) * SW(p,j,i)
            END DO ! p nbasis simd
          END DO ! q nbasis 
        END DO ! j c
        DO i=1,dim
          DO j=1,dim
            !$omp simd
            DO p=1,NBasis
              !A => StiffMatrixTrabsp( p:NBasis*c: NBasis, q:NBasis*c: NBasis )
              !A(j,i) = A(j,i) &
              StiffMatrixTrabsp(((j-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((j-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) + s * Delta * dBasisdx(q,i) * dBasisdx(p,j)
            END DO ! p nbasis simd
        END DO ! j dims
      END DO ! i dims
    END DO ! q nbasis 

    ELSE IF ( Vms ) THEN
      DO i=1,dim
        DO p=1,NBasis
          !omp simd
          DO q=1,NBasis
            !M => MassMatrixTrabsp ( p:NBasis*c: NBasis, q:NBasis*c: NBasis )
            !M(dim+1,i) = M(dim+1,i) &
            MassMatrixTrabsp(((dim+1-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((dim+1-1)*(NBasis)) + (p), & 
            ((i-1)*(NBasis)) + (q)) + s * rho*Tau_M*rho*Basis(q)*dBasisdx(p,i)
          END DO ! p nbasis simd
        END DO ! Q nbasis
        DO k=1,dim+1
          DO q=1,NBasis
            !omp simd
            DO p=1,NBasis
              !A => StiffMatrixTrabsp( p:NBasis*c: NBasis, q:NBasis*c: NBasis )
              !A(dim+1,k) = A(dim+1,k) &
              StiffMatrixTrabsp(((dim+1-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((dim+1-1)*(NBasis)) + (p), & 
              ((k-1)*(NBasis)) + (q)) + s * rho*Tau_M*RM(q,i,k)*dBasisdx(p,i)
            END DO ! p nbasis simd
          END DO ! Q nbasis
        END DO ! k dim+1
        DO j=1,dim
          DO q=1,NBasis
            !omp simd
            DO p=1,NBasis
              !M => MassMatrixTrabsp ( p:NBasis*c: NBasis, q:NBasis*c: NBasis )
              !A => StiffMatrixTrabsp( p:NBasis*c: NBasis, q:NBasis*c: NBasis )

              !A(i,j) = A(i,j) &
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) + s * rho * Tau_M * PRM(i) * Basis(q) * dBasisdx(p,j)/2
              !A(i,i) = A(i,i) &
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) - s * rho * Tau_M * PRM(j) * dBasisdx(q,j) * Basis(p)/2
              !A(i,j) = A(i,j) &
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) - s * rho * Tau_M * rho*Force(i) * Basis(q) * dBasisdx(p,j)
              !A(i,i) = A(i,i) &
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) + s * rho * Tau_M * rho*Force(j) * dBasisdx(q,j) * Basis(p)
              !A(i,j) = A(i,j) &
              StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) + s * rho*Tau_C*RM(q,dim+1,j) * dBasisdx(p,i)

              !M(i,j) = M(i,j) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) - s * rho * Tau_M * rho*Basis(q) * Grad(i,j) * Basis(p)
              !M(i,i) = M(i,i) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) + s * rho * Tau_M * rho*Basis(q) * Velo(j) * dBasisdx(p,j)
              !M(i,j) = M(i,j) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) - s * rho * Tau_M**2 * PRM(i) * rho*Basis(q) * dBasisdx(p,j)
              !M(i,i) = M(i,i) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) - s * rho * Tau_M**2 * rho*Basis(q) * PRM(j) * dBasisdx(p,j)
              !M(i,i) = M(i,i) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) + s * rho * Tau_M**2 * rho*Basis(q) * rho*Force(j) * dBasisdx(p,j)
              !M(i,j) = M(i,j) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) + s * rho * Tau_M**2 * rho*Basis(q) * rho*Force(i) * dBasisdx(p,j)
              !M(i,i) = M(i,i) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((i-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((i-1)*(NBasis)) + (q)) - s * rho * Tau_M**2 * rho*Basis(q) * PVelo(j) * dBasisdx(p,j)
              !M(i,j) = M(i,j) &
              MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q)) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
              ((j-1)*(NBasis)) + (q)) - s * rho * Tau_M**2 * PVelo(i) * rho*Basis(q) * dBasisdx(p,j)

            END DO ! p nbasis simd
          END DO ! Q nbasis
          DO k=1,dim+1
            DO q=1,NBasis
              !omp simd
              DO p=1,NBasis
                !A => StiffMatrixTrabsp( p:NBasis*c: NBasis, q:NBasis*c: NBasis )
                !A(i,k) = A(i,k) &
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((k-1)*(NBasis)) + (q)) - s * rho * Tau_M * RM(q,j,k) * Grad(i,j) * Basis(p)/2
                !A(i,k) = A(i,k) &
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((k-1)*(NBasis)) + (q)) + s * rho * Tau_M * RM(q,i,k) * Velo(j) * dBasisdx(p,j)/2
                !A(i,k) = A(i,k) &
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((k-1)*(NBasis)) + (q)) - s * rho * Tau_M**2 * RM(q,i,k) * PRM(j) * dBasisdx(p,j)/2
                !A(i,k) = A(i,k) &
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((k-1)*(NBasis)) + (q)) - s * rho * Tau_M**2 * PRM(i) * RM(q,j,k) * dBasisdx(p,j)/2
                !A(i,k) = A(i,k) &
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((k-1)*(NBasis)) + (q)) + s * rho * Tau_M**2 * RM(q,i,k) * rho*Force(j) * dBasisdx(p,j)
                !A(i,k) = A(i,k) &
                StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((k-1)*(NBasis)) + (q)) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p), & 
                ((k-1)*(NBasis)) + (q)) + s * rho * Tau_M**2 * rho*Force(i) * RM(q,j,k) * dBasisdx(p,j)
              END DO ! p nbasis simd
            END DO ! Q nbasis
          END DO ! k dim+1
        END DO ! j dim
      END DO ! i dim
    END IF ! VMS
!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------
!     LSx=0
!     IF ( ASSOCIATED( LESVar ) ) THEN
!       DO i=1,dim
!       DO j=1,dim
!          LSx(i,j) = LSx(i,j) + SUM( Basis(1:n) * LES(i,j,1:n) )
!       END DO
!       END DO
!     ELSE
!       DO i=1,dim
!       DO j=1,dim
!          DO k=1,dim
!             LSx(i,j) = LSx(i,j) + Grad(i,k) * Grad(j,k)
!          END DO
!       END DO
!       END DO
!     END IF
!#endif

     DO p=1,NBasis
       Load => ForceVector( c*(p-1)+1 : c*(p-1)+c )

       DO i=1,c
         Load(i) = Load(i) + s * rho * Force(i) * Basis(p)
       END DO

!#ifdef LES
!       DO i=1,dim
!         DO j=1,dim
!            Load(i) = Load(i)+s*rho*LSDelta**2/(2*LSGamma)*LSx(i,j)*dBasisdx(p,j)
!         END DO
!       END DO
!#endif

       IF ( PseudoCompressible ) THEN
          Load(c) = Load(c) + s * Pressure * Basis(p) * Compress
       END IF

       IF ( Stabilize ) THEN
         DO i=1,DIM
           DO j=1,c
             Load(j) = Load(j) + s * Tau * rho * Force(i) * SW(p,j,i)
           END DO
         END DO
       ELSE IF ( Vms ) THEN
         DO i=1,dim
           DO j=1,dim
             Load(i) = Load(i) - s * rho*Tau_M**2 * rho*Force(i) * rho*Force(j) * dBasisdx(p,j)
           END DO
           Load(dim+1) = Load(dim+1) + s * rho*Tau_M * rho*Force(i) * dBasisdx(p,i)
         END DO
       END IF
     END DO
   END DO 

! untranspose the matrixes once the integration is done 
    DO i=1,c
      DO j = 1,c
        DO q=1,NBasis
          DO p=1,NBasis
            StiffMatrix( (c*(p-1))+i, (c*(q-1))+j ) = StiffMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q))
            MassMatrix( (c*(p-1))+i, (c*(q-1))+j ) = MassMatrixTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q))
            JacM( (c*(p-1))+i, (c*(q-1))+j ) = JacMTrabsp(((i-1)*(NBasis)) + (p),((j-1)*(NBasis)) + (q))
          END DO
        END DO
      END DO
    END DO

   IF ( ViscNewtonLin ) THEN
     SOL=0._dp
     SOL(1:c*n:c) = Ux(1:n)
     SOL(2:c*n:c) = Uy(1:n)
     IF ( dim>2 .OR. OutOfPlaneFlow ) SOL(3:c*n:c) = Uz(1:n)
     p = c*nBasis
     StiffMatrix(1:p,1:p) = StiffMatrix(1:p,1:p)+JacM(1:p,1:p)
     FORCEvector(1:p)=Forcevector(1:p)+MATMUL(JacM(1:p,1:p),SOL(1:p))
   END IF

   ! Original P2P1 implementation: Now commented out and replaced by an alternate
   ! implementation since this appears to produce numerical oscillations. To be
   ! removed?
   IF (.FALSE.) THEN
   !IF (  P2P1 ) THEN
     j = GetElementFamily()
     EdgeMap => GetEdgeMap(j)

     DO i=j+1,j+SIZE(EdgeMap,1)
       p = EdgeMap(i-j,1)
       q = EdgeMap(i-j,2)
       StiffMatrix( c*i, : ) = 0._dp
       MassMatrix(  c*i, : ) = 0._dp
       ForceVector( c*i ) = 0._dp
       StiffMatrix( c*i, c*i ) =  1._dp
       StiffMatrix( c*i, c*p ) = -0.5_dp
       StiffMatrix( c*i, c*q ) = -0.5_dp
     END DO
   END IF
   
   ! Produce zero values for pressure at the nodes which are not needed
   ! in the lowest-order pressure interpolation. We shall return consistent
   ! values after the nonlinear iteration has terminated:
   IF (P2P1) THEN
      DO i=LinearBasis+1,nBasis
        StiffMatrix( c*i, : ) = 0._dp
        MassMatrix(  c*i, : ) = 0._dp
        StiffMatrix( :, c*i ) = 0._dp
        MassMatrix(  :, c*i ) = 0._dp
        ForceVector( c*i ) = 0._dp
        StiffMatrix( c*i, c*i ) = 1._dp
     END DO
   END IF

   IF ( PBubbles ) THEN
      DO i=n+1,nBasis
        StiffMatrix( c*i, : ) = 0._dp
        MassMatrix(  c*i, : ) = 0._dp
        StiffMatrix( :, c*i ) = 0._dp
        MassMatrix(  :, c*i ) = 0._dp
        ForceVector( c*i ) = 0._dp
        StiffMatrix( c*i, c*i ) = 1._dp
     END DO
   END IF
!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for Navier-Stokes-equations
!>  boundary conditions in cartesian coordinates.
!------------------------------------------------------------------------------
 SUBROUTINE NavierStokesBoundary( BoundaryMatrix,BoundaryVector,LoadVector,   &
    NodalAlpha, NodalBeta, NodalExtPressure, NodalSlipCoeff, NormalTangential, Element, n, Nodes )
             
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
     NodalAlpha(:),NodalBeta(:), NodalSlipCoeff(:,:), NodalExtPressure(:)

   INTEGER :: n,pn

   TYPE(Element_t),POINTER  :: Element, Parent
   TYPE(Nodes_t)    :: Nodes, ParentNodes

   LOGICAL :: NormalTangential

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
   REAL(KIND=dp) :: detJ,FlowStress(3,3),SlipCoeff
#if 0
   REAL(KIND=dp) :: PBasis(pn),PdBasisdx(pn,3)
#endif

   REAL(KIND=dp) :: u,v,w,ParentU,ParentV,ParentW,s,x(n),y(n),z(n)
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)
   REAL(KIND=dp) :: TangentForce(3),Force(3),Normal(3),Tangent(3),Tangent2(3), &
               Vect(3), Alpha, mu,Grad(3,3),Velo(3)

   REAL(KIND=dp) :: xx, yy, ydot, ydotdot, MassFlux

   INTEGER :: i,j,k,l,k1,k2,t,q,p,c,dim,N_Integ

   LOGICAL :: stat, Found

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

   dim = CoordinateSystemDimension()
   c = dim + 1
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
     stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                 Basis, dBasisdx )

     s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!    Add to load: tangential derivative of something
!------------------------------------------------------------------------------
     DO i=1,dim
       TangentForce(i) = SUM(NodalBeta(1:n)*dBasisdx(1:n,i))
     END DO

!------------------------------------------------------------------------------
!    Add to load: given force in coordinate directions
!------------------------------------------------------------------------------
     Force = 0.0d0
     DO i=1,dim
       Force(i) = Force(i) + SUM( LoadVector(i,1:n)*Basis )
     END DO

!------------------------------------------------------------------------------
!    Add to load: given force in normal direction
!------------------------------------------------------------------------------
     Normal = NormalVector( Element, Nodes, u,v,.TRUE. )

     Alpha = SUM( NodalExtPressure(1:n) * Basis )
     IF ( NormalTangential ) THEN
       Force(1) = Force(1) + Alpha
     ELSE
        DO i=1,dim
           Force(i) = Force(i) + Alpha * Normal(i)
        END DO
     END IF

!------------------------------------------------------------------------------

     Alpha = SUM( NodalAlpha(1:n) * Basis )
     MassFlux = SUM( Loadvector(4,1:n) * Basis(1:n) )

!------------------------------------------------------------------------------

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
           DO i=1,dim
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

     DO q=1,n
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
               BoundaryVector(l) = BoundaryVector(l) + &
                 s * Basis(q) * Force(i) * Vect(j)
            END DO
         ELSE
            BoundaryVector(k) = BoundaryVector(k) + s*Basis(q)*Force(i)
         END IF
         BoundaryVector(k) = BoundaryVector(k) - s * Alpha * dBasisdx(q,i)
         BoundaryVector(k) = BoundaryVector(k) + s * TangentForce(i)*Basis(q)
       END DO
       BoundaryVector(q*c) = BoundaryVector(q*c) + s * MassFlux * Basis(q)
     END DO

   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE VMSWalls( BoundaryMatrix,BoundaryVector,LayerThickness,&
    SurfaceRoughness,Nodalmu,Nodalrho,Ux,Uy,Uz,Element,n,Nodes )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Return element local matrices and RSH vector for Navier-Stokes-equations
!  boundary conditions.
!
!  ARGUMENTS:
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LayerThickness
!     INPUT: Boundary layer thickness
!
!  REAL(KIND=dp) :: SurfaceRoughness
!     INPUT: Measure of surface roughness (f.ex. 9)
!
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values of viscosity
!
!  REAL(KIND=dp) :: Nodalrho(:,:)
!     INPUT: Nodal values of density
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity from previous iteration
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
!******************************************************************************
!------------------------------------------------------------------------------

   IMPLICIT NONE

   REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
     Nodalmu(:),Nodalrho(:),Ux(:),Uy(:),Uz(:)

   REAL(KIND=dp) :: LayerThickness(:),SurfaceRoughness(:)

   INTEGER :: n

   TYPE(Element_t),POINTER  :: Element
   TYPE(Nodes_t)    :: Nodes

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
   REAL(KIND=dp) :: detJ

   REAL(KIND=dp) :: u,v,w,s,r,rho,mu,Dist,Roughness,SlipCoeff,Vect(3),vabs
   REAL(KIND=dp) :: x,y,z,h
   REAL(KIND=dp) :: Velo(3),Normal(3),Tangent(3), Tangent2(3)
   REAL(KIND=dp) :: TangentialVelocity(2), FrictionVelocity,DFX,DKERR
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   INTEGER, POINTER :: Ind(:)
   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: i,j,k,k1,k2,t,q,p,c,dim,N_Integ

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

   INTERFACE
     SUBROUTINE SOLVE_UFRIC(DENSIT,VISCOS,DIST,ROUGH,UT,UFRIC,DFX)
         DOUBLE PRECISION :: Densit, Viscos,Dist, Rough, Ut, Ufric, DFX
     END SUBROUTINE SOLVE_UFRIC
   END INTERFACE
!------------------------------------------------------------------------------

   IF ( GetString(GetSolverParams(), 'Stabilization Method' )/='vms' ) RETURN

   dim = CoordinateSystemDimension()
   c = dim + 1

   Mesh => GetMesh()
   Ind => Element % BoundaryInfo % Left % NodeIndexes
   k = 0
   x = 0; y=0; z=0
   DO i=1,SIZE(ind)
     DO j=1,n
       IF ( ind(i)==Element % NodeIndexes(j) ) EXIT
     END DO
     IF ( j>n ) THEN
       x = x+Mesh % Nodes % x(Ind(i))
       y = y+Mesh % Nodes % y(Ind(i))
       z = z+Mesh % Nodes % z(Ind(i))
       k = k + 1
     END IF
   END DO

   IF ( k>0 ) THEN
     x = x/k - SUM(Nodes % x)/n
     y = y/k - SUM(Nodes % y)/n
     z = z/k - SUM(Nodes % z)/n
   END IF
   h = SQRT(x**2 + y**2 + z**2)

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
     stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx )

     s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!    rho and viscosity at integration point
!------------------------------------------------------------------------------
     rho   = SUM( Nodalrho(1:n)*Basis(1:n) )
     mu = SUM( Nodalmu(1:n)*Basis(1:n) )
!------------------------------------------------------------------------------
!    Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
     Velo = 0.0d0
     Velo(1) = SUM( Ux(1:n)*Basis(1:n) )
     Velo(2) = SUM( Uy(1:n)*Basis(1:n) )
     IF ( dim > 2 ) Velo(3) = SUM( Uz(1:n)*Basis(1:n) )
!------------------------------------------------------------------------------
!    Normal & tangent directions
!------------------------------------------------------------------------------
     Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
 
     IF ( dim <= 2 ) THEN
       Tangent(1) =  Normal(2)
       Tangent(2) = -Normal(1)
       Tangent(3) =  0.0d0
     ELSE
       CALL TangentDirections( Normal, Tangent, Tangent2 ) 
     END IF
     TangentialVelocity(1) = SUM( Velo(1:dim) * Tangent(1:dim) )

     IF ( dim==3 ) THEN
       TangentialVelocity(2) = SUM( Velo(1:dim) * Tangent2(1:dim) )
     END IF

     Dist = SUM( LayerThickness(1:n) * Basis(1:n) )
     IF ( Dist == 0 ) Dist = h
     Roughness = SUM( SurfaceRoughness(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!    Solve friction velocity and its derivative with respect to
!    the tangential velocity:
!------------------------------------------------------------------------------
     FrictionVelocity=0
     DKERR=0
     Vabs = MAX( SQRT(SUM(Velo**2)), 1.d-8 )
     CALL Solve_UFric( rho,mu,Dist,Roughness, &
          Vabs,FrictionVelocity,DFX )
     DKERR = rho*FrictionVelocity**2

!------------------------------------------------------------------------------
    DO p=1,n
    DO q=1,n
      DO i=1,dim
        k1 = c*(p-1)+i
        k2 = c*(q-1)+i
        BoundaryMatrix(k1,k2) = BoundaryMatrix(k1,k2) - &
              s*2*mu*SUM(dBasisdx(q,:)*Normal)*Basis(p)

        BoundaryMatrix(k1,k2) = BoundaryMatrix(k1,k2) - &
              s*2*mu*SUM(dBasisdx(p,:)*Normal)*Basis(q)
      END DO
    END DO
    END DO


     DO p=1,n
       DO q=1,n
         DO i=1,dim
           k1 = c*(p-1)+i
           k2 = c*(q-1)+i
           BoundaryMatrix(k1,k2) = BoundaryMatrix(k1,k2) + &
                  s * DKERR * Basis(q) * Basis(p)
         END DO
       END DO
     END DO
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE VmsWalls
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Return the local matrices and RHS contrbution for the wall law.
!------------------------------------------------------------------------------
 SUBROUTINE NavierStokesWallLaw( BoundaryMatrix,BoundaryVector,LayerThickness,&
    SurfaceRoughness,Nodalmu,Nodalrho,Ux,Uy,Uz,Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LayerThickness
!     INPUT: Boundary layer thickness
!
!  REAL(KIND=dp) :: SurfaceRoughness
!     INPUT: Measure of surface roughness (f.ex. 9)
!
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values of viscosity
!
!  REAL(KIND=dp) :: Nodalrho(:,:)
!     INPUT: Nodal values of density
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity from previous iteration
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

   REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
     Nodalmu(:),Nodalrho(:),Ux(:),Uy(:),Uz(:)

   REAL(KIND=dp) :: LayerThickness(:),SurfaceRoughness(:)

   INTEGER :: n

   TYPE(Element_t),POINTER  :: Element
   TYPE(Nodes_t)    :: Nodes

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
   REAL(KIND=dp) :: detJ

   REAL(KIND=dp) :: u,v,w,s,r,rho,mu,Dist,Roughness,SlipCoeff,Vect(3),vabs
   REAL(KIND=dp) :: Velo(3),Normal(3),Tangent(3), Tangent2(3)
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)
   REAL(KIND=dp) :: TangentialVelocity(2), FrictionVelocity(2),DFX(2),DKERR(2)

   INTEGER :: i,j,k,k1,k2,t,q,p,c,dim,N_Integ

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

   INTERFACE
     SUBROUTINE SOLVE_UFRIC(DENSIT,VISCOS,DIST,ROUGH,UT,UFRIC,DFX)
         DOUBLE PRECISION :: Densit, Viscos,Dist, Rough, Ut, Ufric, DFX
     END SUBROUTINE SOLVE_UFRIC
   END INTERFACE
!------------------------------------------------------------------------------

   dim = CoordinateSystemDimension()
   c = dim + 1

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
     stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
           Basis,dBasisdx )

     s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!    rho and viscosity at integration point
!------------------------------------------------------------------------------
     rho   = SUM( Nodalrho(1:n)*Basis(1:n) )
     mu = SUM( Nodalmu(1:n)*Basis(1:n) )
!------------------------------------------------------------------------------
!    Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
     Velo = 0.0d0
     Velo(1) = SUM( Ux(1:n)*Basis(1:n) )
     Velo(2) = SUM( Uy(1:n)*Basis(1:n) )
     IF ( dim > 2 ) Velo(3) = SUM( Uz(1:n)*Basis(1:n) )
     vabs = MAX(1.d-9,SQRT(SUM(Velo**2)))
!------------------------------------------------------------------------------
!    Normal & tangent directions
!------------------------------------------------------------------------------
     Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
 
     IF ( dim <= 2 ) THEN
       Tangent(1) =  Normal(2)
       Tangent(2) = -Normal(1)
       Tangent(3) =  0.0d0
     ELSE
       CALL TangentDirections( Normal, Tangent, Tangent2 ) 
     END IF


     TangentialVelocity(1) = SUM( Velo(1:dim) * Tangent(1:dim) )
     IF ( TangentialVelocity(1) < 0 ) THEN
       Tangent = -Tangent
       TangentialVelocity(1) = -TangentialVelocity(1)
     END IF
     IF ( dim==3 ) THEN
       TangentialVelocity(2) = SUM( Velo(1:dim) * Tangent2(1:dim) )
       IF ( TangentialVelocity(2) < 0 ) THEN
         Tangent2 = -Tangent2
         TangentialVelocity(2) = -TangentialVelocity(2)
       END IF
     END IF

     Dist = SUM( LayerThickness(1:n) * Basis(1:n) )
     Roughness = SUM( SurfaceRoughness(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!    Solve friction velocity and its derivative with respect to
!    the tangential velocity:
!------------------------------------------------------------------------------
     FrictionVelocity=0
     DKERR=0
     IF ( TangentialVelocity(1) > 1.d-9 ) THEN
       CALL Solve_UFric( rho,mu,Dist,Roughness, &
           TangentialVelocity(1),FrictionVelocity(1),DFX(1) )
       DKERR(1) = 2.0d0 * rho * FrictionVelocity(1) / DFX(1)
     END IF

     IF ( dim==3 ) THEN
       IF ( TangentialVelocity(2) > 1.d-9 ) THEN
         CALL Solve_UFric( rho,mu,Dist,Roughness, &
             TangentialVelocity(2),FrictionVelocity(2),DFX(2) )
         DKERR(2) = 2.0d0 * rho * FrictionVelocity(2) / DFX(2)
       END IF
     END IF

!------------------------------------------------------------------------------
     DO p=1,n
       DO q=1,n
         DO i=1,dim
           DO j=1,dim
             k1 = (p-1)*c + i
             k2 = (q-1)*c + j
             BoundaryMatrix(k1,k2) = BoundaryMatrix(k1,k2) + &
               s * DKERR(1) * Tangent(i) * Tangent(j) * Basis(q) * Basis(p)
             IF ( dim==3 ) THEN
               BoundaryMatrix(k1,k2) = BoundaryMatrix(k1,k2) + &
                 s * DKERR(2) * Tangent2(i) * Tangent2(j) * Basis(q) * Basis(p)
             END IF
           END DO
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
     DO q=1,n
       DO i=1,dim
         k1 = (q-1)*c + i
         BoundaryVector(k1) = BoundaryVector(k1) + &
           s * ( DKERR(1) * TangentialVelocity(1) - &
             rho * FrictionVelocity(1)**2 ) * Tangent(i) * Basis(q)
         IF ( dim==3 ) THEN
           BoundaryVector(k1) = BoundaryVector(k1) + &
             s * ( DKERR(2) * TangentialVelocity(2) - &
               rho * FrictionVelocity(2)**2 ) * Tangent2(i) * Basis(q)
         END IF
       END DO
     END DO
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesWallLaw
!------------------------------------------------------------------------------

END MODULE NavierStokes

!> \}
