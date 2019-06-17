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
!>  Module computing local matrices for stress computation (cartesian
!>  coordinates, axisymmetric)
!------------------------------------------------------------------------------


MODULE StressLocal

!------------------------------------------------------------------------------
  USE DefUtils


  IMPLICIT NONE

!------------------------------------------------------------------------------
  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE StressCompose( MASS, DAMP, STIFF, FORCE, FORCE_im, LOAD, LOAD_im, ElasticModulus, &
     NodalPoisson, NodalDensity, PlaneStress, Isotropic,           &
     NodalPreStress, NodalPreStrain, NodalStressLoad, NodalStrainLoad,           &
     NodalHeatExpansion, NodalTemperature, Element, n, ntot, Nodes, RelIntegOrder, StabilityAnalysis, &
     GeometricStiffness, NodalDisplacement, RotateC, TransformMatrix, NodalMeshVelo, &
     NodalDamping, RayleighDamping, RayleighAlpha, RayleighBeta  )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:), MASS(:,:), DAMP(:,:), FORCE(:), LOAD(:,:)
     REAL(KIND=dp) :: FORCE_im(:), LOAD_im(:,:)
     REAL(KIND=dp) :: NodalTemperature(:),ElasticModulus(:,:,:)
     REAL(KIND=dp) :: NodalPreStress(:,:), NodalPreStrain(:,:)
     REAL(KIND=dp) :: NodalStressLoad(:,:), NodalStrainLoad(:,:)
     REAL(KIND=dp) :: NodalDisplacement(:,:), NodalHeatExpansion(:,:,:)
     REAL(KIND=dp) :: TransformMatrix(:,:), NodalMeshVelo(:,:)
     REAL(KIND=dp) :: RayleighAlpha(:), RayleighBeta(:)
     REAL(KIND=dp), DIMENSION(:) :: NodalPoisson, NodalDensity, NodalDamping

     LOGICAL :: PlaneStress, Isotropic(2), StabilityAnalysis, GeometricStiffness
     LOGICAL :: RotateC, RayleighDamping


     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: RelIntegOrder

     INTEGER :: n, ntot
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(ntot)
     REAL(KIND=dp) :: dBasisdx(ntot,3),detJ

     REAL(KIND=dp) :: LoadAtIp(3), LoadatIp_im(3), Poisson, Young, Ident(3,3)

     REAL(KIND=dp) :: M(3,3),D(3,3),HeatExpansion(3,3), A(4,4)
     REAL(KIND=dp) :: Temperature,Density, C(6,6), Damping,MeshVelo(3)
     REAL(KIND=dp) :: StressTensor(3,3), StrainTensor(3,3), InnerProd, NodalViscosity(n)
     REAL(KIND=dp) :: StressLoad(6), StrainLoad(6), PreStress(6), PreStrain(6)

     INTEGER :: i,j,k,l,p,q,t,dim,NBasis,ind(3)

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), xPhi

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry, NeedMass, NeedHeat, NeedStress, NeedHarmonic, &
         NeedPreStress, ActiveGeometricStiffness, GPA


     TYPE(ValueList_t), POINTER :: BF
   
     REAL(KIND=dp) :: GPA_Coeff(n)

     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: ndim
     LOGICAL :: Found, Incompressible, MaxwellMaterial
     REAL(KIND=dp) :: Pres, Pres0
     REAL(KIND=dp) :: PSOL(4,32), SOL(4,32), ShearModulus, Viscosity, PrevStress(3,3)
!------------------------------------------------------------------------------

     TYPE(Variable_t), POINTER, SAVE :: ve_stress

     REAL(KIND=dp), ALLOCATABLE, SAVE :: StressStore(:,:,:,:)


     dim = CoordinateSystemDimension()

     Incompressible = GetLogical( GetSolverParams(), 'Incompressible', Found )
     IF (Incompressible) THEN
       ndim = dim+1
     ELSE
       ndim = dim
     END IF

     Ident = 0._dp
     DO i=1,dim
       Ident(i,i) = 1._dp
     END DO

     CSymmetry = .FALSE.
     CSymmetry = CSymmetry .OR. CurrentCoordinateSystem() == AxisSymmetric
     CSymmetry = CSymmetry .OR. CurrentCoordinateSystem() == CylindricSymmetric

     FORCE_im = 0.0d0
     FORCE = 0.0d0
     STIFF = 0.0d0
     MASS  = 0.0d0
     DAMP  = 0.0d0

     NeedMass = ANY( NodalDensity(1:n) /= 0.0d0 )
     NeedMass = NeedMass .OR. ANY( NodalDamping(1:n) /= 0.0d0 ) .OR. RayleighDamping

     NeedHeat = ANY( NodalTemperature(1:n) /= 0.0d0 )
     NeedHarmonic = ANY( LOAD_im(:,1:n) /= 0.0d0 ) 
     NeedPreStress = ANY( NodalPreStrain(1:6,1:n) /= 0.0d0 ) 
     NeedPreStress = NeedPreStress .OR. ANY( NodalPreStress(1:6,1:n) /= 0.0d0 ) 


     BF => GetBodyForce()
     GPA = .FALSE.
     IF(ASSOCIATED(BF)) THEN
        GPA = GetLogical(BF, 'Gravitational Prestress Advection', Found )
       IF ( GPA ) THEN
         GPA_Coeff(1:n) = GetReal( BF, 'GPA Coeff', Found )
       END IF
     END IF


     !      ! Integration stuff:
     ! ------------------  
     NBasis = ntot
     IntegStuff = GaussPoints( element, RelOrder = RelIntegOrder )

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     Mesh => GetMesh()
     MaxwellMaterial = GetLogical( GetMaterial(), 'Maxwell material', Found )
     IF( MaxwellMaterial ) THEN
       ve_stress => variableget( Mesh % Variables, 've_stress' )
       IF(.NOT.ASSOCIATED(ve_stress)) THEN
         CALL Fatal( 'StressCompose', '"Maxwell material" set, but no storage space for stresses present?' )
       END IF

       i = Element % ElementIndex
       j = ve_stress % Perm( i+1 ) - ve_stress % Perm ( i )
       IF( IntegStuff % n /= j ) THEN
         PRINT *,'Inconsistent number of gauss points:',i, IntegStuff % n, j
       END IF

       NodalViscosity(1:n) = GetReal( GetMaterial(), 'Viscosity', Found )

       SOL = 0; PSOL = 0
       CALL GetVectorLocalSolution( SOL )
       CALL GetVectorLocalSolution( PSOL, tStep=-2 )
       PSOL = SOL - PSOL
       PrevStress = 0._dp
     END IF

     ! Now we start integrating:
     ! -------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------

       IF ( NeedMass ) THEN
         Density = SUM( NodalDensity(1:n)*Basis(1:n) )
         Damping = SUM( NodalDamping(1:n)*Basis(1:n) )
         DO i=1,dim
           MeshVelo(i) = SUM( NodalMeshVelo(i,1:n)*Basis(1:n) )
         END DO
       END IF


       IF ( NeedHeat ) THEN
         ! Temperature at the integration point:
         !-------------------------------------- 
         Temperature = SUM( NodalTemperature(1:ntot)*Basis(1:ntot) )
 
         ! Heat expansion tensor values at the integration point:
         !-------------------------------------------------------
         HeatExpansion = 0.0d0
         DO i=1,3
           IF ( Isotropic(2) ) THEN
              HeatExpansion(i,i) = SUM( NodalHeatExpansion(1,1,1:n)*Basis(1:n) )
           ELSE
              DO j=1,3
                HeatExpansion(i,j) = SUM( NodalHeatExpansion(i,j,1:n)*Basis(1:n) )
              END DO
           END IF
         END DO
       END IF

       IF ( Isotropic(1) ) Poisson = SUM( Basis(1:n) * NodalPoisson(1:n) )

       C = 0
       IF ( .NOT. Isotropic(1) ) THEN 
          DO i=1,SIZE(ElasticModulus,1)
            DO j=1,SIZE(ElasticModulus,2)
               C(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
            END DO
          END DO
       ELSE
          Young = SUM( Basis(1:n) * ElasticModulus(1,1,1:n) )
       END IF

       SELECT CASE(dim)
       CASE(2)
         IF ( CSymmetry ) THEN
           IF ( Isotropic(1) ) THEN
              C(1,1) = 1.0d0 - Poisson
              C(1,2) = Poisson
              C(1,3) = Poisson
              C(2,1) = Poisson
              C(2,2) = 1.0d0 - Poisson
              C(2,3) = Poisson
              C(3,1) = Poisson
              C(3,2) = Poisson
              C(3,3) = 1.0d0 - Poisson
              C(4,4) = 0.5d0 - Poisson

              C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
           END IF
           Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
           s = s * Radius
         ELSE
           IF ( Isotropic(1) ) THEN
              IF ( PlaneStress ) THEN
                 C(1,1) = 1.0d0
                 C(1,2) = Poisson
                 C(2,1) = Poisson
                 C(2,2) = 1.0d0
                 C(3,3) = 0.5d0*(1-Poisson)
 
                 C = C * Young / ( 1 - Poisson**2 )
              ELSE
                 C(1,1) = 1.0d0 - Poisson
                 C(1,2) = Poisson
                 C(2,1) = Poisson
                 C(2,2) = 1.0d0 - Poisson
                 C(3,3) = 0.5d0 - Poisson
                 C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
              END IF
           ELSE
              IF ( PlaneStress ) THEN
                C(1,1) = C(1,1) - C(1,3)*C(3,1) / C(3,3)
                C(1,2) = C(1,2) - C(1,3)*C(2,3) / C(3,3)
                C(2,1) = C(2,1) - C(1,3)*C(2,3) / C(3,3)
                C(2,2) = C(2,2) - C(2,3)*C(3,2) / C(3,3)
              ELSE
                IF ( NeedHeat ) THEN
                  HeatExpansion(1,1) = HeatExpansion(1,1) + HeatExpansion(3,3) * &
                     ( C(2,2)*C(1,3)-C(1,2)*C(2,3) ) / ( C(1,1)*C(2,2) - C(1,2)*C(2,1) )
  
                  HeatExpansion(2,2) = HeatExpansion(2,2) + HeatExpansion(3,3) * &
                     ( C(1,1)*C(2,3)-C(1,2)*C(1,3) ) / ( C(1,1)*C(2,2) - C(1,2)*C(2,1) )
                END IF
              END IF
              C(3,3) = C(4,4)
              C(1,3) = 0.0d0
              C(3,1) = 0.0d0
              C(2,3) = 0.0d0
              C(3,2) = 0.0d0
              C(4:6,:) = 0.0d0
              C(:,4:6) = 0.0d0
           END IF
         END IF

       CASE(3)
         IF ( Isotropic(1) ) THEN
            C = 0
            C(1,1) = 1.0d0 - Poisson
            C(1,2) = Poisson
            C(1,3) = Poisson
            C(2,1) = Poisson
            C(2,2) = 1.0d0 - Poisson
            C(2,3) = Poisson
            C(3,1) = Poisson
            C(3,2) = Poisson
            C(3,3) = 1.0d0 - Poisson
            C(4,4) = 0.5d0 - Poisson
            C(5,5) = 0.5d0 - Poisson
            C(6,6) = 0.5d0 - Poisson

            C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
!--------------------------------------------------------------------------------
!   Rotate elasticity tensor if required

          ELSE
            IF ( RotateC ) THEN
              CALL RotateElasticityMatrix( C, TransformMatrix, 3 )
            END IF
          END IF

       END SELECT

       ActiveGeometricStiffness = StabilityAnalysis.OR.GeometricStiffness
       IF ( ActiveGeometricStiffness ) THEN
         CALL LocalStress( StressTensor,StrainTensor,NodalPoisson,ElasticModulus, &
             NodalHeatExpansion, NodalTemperature, Isotropic,CSymmetry,PlaneStress,   &
             NodalDisplacement,Basis,dBasisdx,Nodes,dim,n,ntot )
       END IF

       IF( NeedPreStress ) THEN
         DO i=1,6
           PreStrain(i) = SUM( NodalPreStrain(i,1:n)*Basis(1:n) )
           PreStress(i) = SUM( NodalPreStress(i,1:n)*Basis(1:n) )
         END DO
         PreStress = PreStress + MATMUL( C, PreStrain  )
         
         DO i=1,6
           StrainLoad(i) = SUM( NodalStrainLoad(i,1:n)*Basis(1:n) )
           StressLoad(i) = SUM( NodalStressLoad(i,1:n)*Basis(1:n) )
         END DO
         StressLoad = StressLoad + MATMUL( C, StrainLoad )
         
         IF( .NOT. ( StabilityAnalysis .OR. GeometricStiffness ) ) THEN 
           StressTensor = 0.0d0
           StrainTensor = 0.0d0          
         END IF
         
         SELECT CASE(dim)
         CASE(2)
           IF ( Csymmetry ) THEN
             StressTensor(1,1) = StressTensor(1,1) + PreStress(1)
             StressTensor(2,2) = StressTensor(2,2) + PreStress(2)
             StressTensor(3,3) = StressTensor(3,3) + PreStress(3)
             StressTensor(1,2) = StressTensor(1,2) + PreStress(4)
             StressTensor(2,1) = StressTensor(2,1) + PreStress(4)
           ELSE
             StressTensor(1,1) = StressTensor(1,1) + PreStress(1)
             StressTensor(2,2) = StressTensor(2,2) + PreStress(2)
             StressTensor(1,2) = StressTensor(1,2) + PreStress(3)
             StressTensor(2,1) = StressTensor(2,1) + PreStress(3)
           END IF
         CASE(3)
           StressTensor(1,1) = StressTensor(1,1) + PreStress(1)
           StressTensor(2,2) = StressTensor(2,2) + PreStress(2)
           StressTensor(3,3) = StressTensor(3,3) + PreStress(3)
           StressTensor(1,2) = StressTensor(1,2) + PreStress(4)
           StressTensor(2,1) = StressTensor(2,1) + PreStress(4)
           StressTensor(2,3) = StressTensor(2,3) + PreStress(5)
           StressTensor(3,2) = StressTensor(3,2) + PreStress(5)
           StressTensor(1,3) = StressTensor(1,3) + PreStress(6)
           StressTensor(3,1) = StressTensor(3,1) + PreStress(6)
         END SELECT
       END IF

       IF(MaxwellMaterial) THEN
         Viscosity = SUM( NodalViscosity(1:n) * Basis(1:n) )
         xPhi = ViscoElasticLoad( ve_stress, t, StressLoad )
         NeedPreStress = .TRUE.
       ELSE
         xPhi = 1
       END IF

       !
       ! Loop over basis functions (of both unknowns and weights):
       ! ---------------------------------------------------------
       A = 0.0d0
       M = 0.0d0
       D = 0.0d0
       B = 0.0d0

       DO p=1,NBasis

         G = 0.0d0
         SELECT CASE(dim)
         CASE(2)
           IF ( CSymmetry ) THEN
             G(1,1) = dBasisdx(p,1)
             G(1,3) = Basis(p) / Radius
             G(1,4) = dBasisdx(p,2)
             G(2,2) = dBasisdx(p,2)
             G(2,4) = dBasisdx(p,1)
           ELSE
             G(1,1) = dBasisdx(p,1)
             G(1,3) = dBasisdx(p,2)
             G(2,2) = dBasisdx(p,2)
             G(2,3) = dBasisdx(p,1)
           END IF

         CASE(3)
           G(1,1) = dBasisdx(p,1)
           G(2,2) = dBasisdx(p,2)
           G(3,3) = dBasisdx(p,3)
           G(1,4) = dBasisdx(p,2)
           G(2,4) = dBasisdx(p,1)
           G(2,5) = dBasisdx(p,3)
           G(3,5) = dBasisdx(p,2)
           G(1,6) = dBasisdx(p,3)
           G(3,6) = dBasisdx(p,1)
         END SELECT

         LoadAtIp    = 0.0d0
         LoadAtIp_im = 0.0d0
         IF( NeedPreStress ) THEN
           DO i=1,dim
             DO j=1,6
               LoadAtIp(i) = LoadAtIp(i) + StressLoad(j) * G(i,j)
             END DO
           END DO
         END IF


         IF (.NOT. Incompressible ) G = MATMUL( G, C )

         DO q=1,NBasis
           IF ( NeedMass ) THEN
             DO i=1,dim
               M(i,i) = Density * Basis(p) * Basis(q)
               D(i,i) = Damping * Basis(p) * Basis(q)
               DO j=1,dim
                 D(i,i) = D(i,i) + Density * MeshVelo(j) * dBasisdx(q,j) * Basis(p)
               END DO
             END DO
           END IF
 
           SELECT CASE(dim)
           CASE(2)
              IF ( CSymmetry ) THEN
                 B(1,1) = dBasisdx(q,1)
                 B(2,2) = dBasisdx(q,2)
                 B(3,1) = Basis(q) / Radius
                 B(4,1) = dBasisdx(q,2)
                 B(4,2) = dBasisdx(q,1)
              ELSE
                 B(1,1) = dBasisdx(q,1)
                 B(2,2) = dBasisdx(q,2)
                 B(3,1) = dBasisdx(q,2)
                 B(3,2) = dBasisdx(q,1)
              END IF
 
           CASE(3)
              B(1,1) = dBasisdx(q,1)
              B(2,2) = dBasisdx(q,2)
              B(3,3) = dBasisdx(q,3)
              B(4,1) = dBasisdx(q,2)
              B(4,2) = dBasisdx(q,1)
              B(5,2) = dBasisdx(q,3)
              B(5,3) = dBasisdx(q,2)
              B(6,1) = dBasisdx(q,3)
              B(6,3) = dBasisdx(q,1)
           END SELECT
 
           A = 0._dp
           IF ( .NOT. Incompressible ) THEN
              A(1:3,1:3) = MATMUL( G, B ) * xPhi
           ELSE
              DO i=1,dim 
                DO j=1,dim 
                  A(i,i) = A(i,i) + Young/3 * dBasisdx(q,j) * dBasisdx(p,j)
                  A(i,j) = A(i,j) + Young/3 * dBasisdx(q,i) * dBasisdx(p,j)
                END DO
                A(i,:) = A(i,:) * xPhi

                A(i,ndim) = A(i,ndim) - Basis(q) * dBasisdx(p,i)
                A(ndim,i) = A(ndim,i) - dBasisdx(q,i) * Basis(p)
             END DO
           END IF
 
           IF( GPA ) THEN
             DO i=1,dim
               A(i,dim) = A(i,dim) + SUM(GPA_Coeff(1:n)*Basis(1:n))*dBasisdx(q,i)*Basis(p)
             END DO
           END IF

           !
           ! Add nodal matrix to element matrix:
           ! -----------------------------------
           DO i=1,ndim
             DO j=1,ndim
               STIFF( ndim*(p-1)+i,ndim*(q-1)+j ) =  &
                    STIFF( ndim*(p-1)+i,ndim*(q-1)+j ) + s*A(i,j)
             END DO
           END DO

           IF ( NeedMass .AND. (.NOT.StabilityAnalysis) ) THEN
              DO i=1,dim
                DO j=1,dim
                  MASS( ndim*(p-1)+i,ndim*(q-1)+j ) =  &
                       MASS( ndim*(p-1)+i,ndim*(q-1)+j ) + s*M(i,j)

                  DAMP( ndim*(p-1)+i,ndim*(q-1)+j ) =  &
                       DAMP( ndim*(p-1)+i,ndim*(q-1)+j ) + s*D(i,j)
                END DO
              END DO
           END IF
      
           IF ( ActiveGeometricStiffness ) THEN
             DO k = 1,dim
               InnerProd = 0.0d0
               DO i = 1,dim
                 DO j = 1,dim
                    InnerProd = InnerProd + &
                      dBasisdx(p,i) * dBasisdx(q,j) * StressTensor(i,j)
                 END DO
               END DO

               IF ( StabilityAnalysis ) THEN
                 MASS( ndim*(p-1)+k,ndim*(q-1)+k ) &
                     = MASS( ndim*(p-1)+k,ndim*(q-1)+k ) - s * InnerProd
               ELSE
                 STIFF( ndim*(p-1)+k,ndim*(q-1)+k ) &
                    = STIFF( ndim*(p-1)+k,ndim*(q-1)+k ) + s * InnerProd 
               END IF
             END DO
           END IF

         END DO

         !
         ! The (rest of the) righthand side:
         ! ---------------------------------
         DO i=1,dim
           LoadAtIp(i) = LoadAtIp(i) + &
               SUM( LOAD(i,1:n)*Basis(1:n) ) * Basis(p) + &
               SUM( LOAD(4,1:n)*Basis(1:n) ) * dBasisdx(p,i)
           IF( NeedHarmonic ) THEN
             LoadAtIp_im(i) = LoadAtIp_im(i) + &
                 SUM( LOAD_im(i,1:n)*Basis(1:n) ) * Basis(p) + &
                 SUM( LOAD_im(4,1:n)*Basis(1:n) ) * dBasisdx(p,i)
           END IF
         END DO

         IF ( NeedHeat ) THEN
           DO i=1,dim
             IF ( CSymmetry ) THEN
               DO j=1,3
                 LoadAtIp(i) = LoadAtIp(i) +  &
                   G(i,j) * HeatExpansion(j,j) * Temperature
               END DO
             ELSE
               DO j=1,dim
                 LoadAtIp(i) = LoadAtIp(i) + &
                   G(i,j) * HeatExpansion(j,j) * Temperature
               END DO
             END IF
           END DO
         END IF

         DO i=1,dim
           FORCE(ndim*(p-1)+i) = FORCE(ndim*(p-1)+i) + s*LoadAtIp(i)
           IF( NeedHarmonic ) THEN
             FORCE_im(ndim*(p-1)+i) = FORCE_im(ndim*(p-1)+i) + s*LoadAtIp_im(i)
           END IF
         END DO
      END DO
    END DO

    IF ( Incompressible ) THEN
      DO i=n+1,ntot
        j = ndim*i
        FORCE(j)   = 0._dp
        STIFF(j,:) = 0._dp
        STIFF(:,j) = 0._dp
        STIFF(j,j) = 1._dp
      END DO
    END IF

    DAMP  = ( DAMP  + TRANSPOSE(DAMP) )  / 2.0_dp
    MASS  = ( MASS  + TRANSPOSE(MASS) )  / 2.0_dp
    STIFF = ( STIFF + TRANSPOSE(STIFF) ) / 2.0_dp

    IF( RayleighDamping ) THEN
        DAMP = RayleighAlpha(1) * MASS + RayleighBeta(1) * STIFF
    END IF

CONTAINS

!------------------------------------------------------------------------------
 FUNCTION ViscoElasticLoad(ve_stress, ip, StressLoad) RESULT(xPhi)
!------------------------------------------------------------------------------
    TYPE(Variable_t) :: ve_stress
    INTEGER :: ip
    REAL(KIND=dp) :: StressLoad(6), Xphi
!------------------------------------------------------------------------------
    INTEGER :: i
    REAL(KIND=dp) :: StressTensor(3,3), PrevStress(3,3), Pres, Pres0, &
           ShearModulus

    StressTensor  = 0._dp
    CALL LocalStress( StressTensor,StrainTensor,NodalPoisson,ElasticModulus, &
         NodalHeatExpansion, NodalTemperature, Isotropic,CSymmetry,PlaneStress,   &
         PSOL,Basis,dBasisdx,Nodes,dim,n,ntot, .FALSE. )

    IF(Incompressible) THEN
      ShearModulus = Young / 3
      Pres  = SUM( Basis(1:n) * SOL(ndim,1:n) )
      Pres0 = SUM( Basis(1:n) * (SOL(ndim,1:n) - PSOL(ndim,1:n)) )
    ELSE
      Pres = 0._dp; Pres0 = 0._dp
      ShearModulus = Young / (2*(1+Poisson))
    END IF

    xPhi = 1._dp / ( 1 + ShearModulus / Viscosity * GetTimeStepSize() )

    i = dim**2*(ve_stress % perm(Element % ElementIndex) + ip - 1)
    PrevStress(1:dim,1:dim) = RESHAPE( ve_stress % values(i+1:i+dim**2), [dim,dim] )
    PrevStress = xPhi * (StressTensor + PrevStress + Pres0 * Ident) - Pres * Ident
    ve_stress % values(i+1:i+dim**2) = RESHAPE( PrevStress(1:dim,1:dim), [dim**2] )

    StressTensor  = 0._dp
    CALL LocalStress( StressTensor,StrainTensor,NodalPoisson,ElasticModulus, &
        NodalHeatExpansion, NodalTemperature, Isotropic,CSymmetry,PlaneStress,   &
        SOL,Basis,dBasisdx,Nodes,dim,n,ntot, .FALSE. )

    StressTensor = xPhi * (StressTensor - PrevStress - Pres * Ident)
            
    CALL Tensor26Vector( StressTensor, StressLoad, dim, CSymmetry )
!------------------------------------------------------------------------------
 END FUNCTION ViscoElasticLoad
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 END SUBROUTINE StressCompose
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE StressForceCompose( FORCE, FORCE_im, LOAD, LOAD_im, ElasticModulus, NodalPoisson,     &
     PlaneStress, Isotropic,NodalStressLoad, NodalStrainLoad, NodalHeatExpansion,&
     NodalTemperature, Element, n, ntot, Nodes,  RelIntegOrder, RotateC, TransformMatrix )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: FORCE(:), LOAD(:,:)
     REAL(KIND=dp) :: FORCE_im(:), LOAD_im(:,:)
     REAL(KIND=dp) :: NodalTemperature(:),ElasticModulus(:,:,:)
     REAL(KIND=dp) :: NodalStressLoad(:,:), NodalStrainLoad(:,:)
     REAL(KIND=dp) :: NodalHeatExpansion(:,:,:)
     REAL(KIND=dp) :: TransformMatrix(:,:), NodalPoisson(:)

     LOGICAL :: PlaneStress, Isotropic(2)
     LOGICAL :: RotateC

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: RelIntegOrder

     INTEGER :: n, ntot
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(ntot)
     REAL(KIND=dp) :: dBasisdx(ntot,3),detJ

     REAL(KIND=dp) :: LoadAtIp(3), LoadAtIp_im(3), Poisson, Young

     REAL(KIND=dp), DIMENSION(3,3) :: HeatExpansion
     REAL(KIND=dp) :: Temperature, C(6,6)
     REAL(KIND=dp) :: StressLoad(6), StrainLoad(6)

     INTEGER :: i,j,k,l,p,q,t,dim,NBasis,ind(3)

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6)

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry, NeedHeat

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     CSymmetry = .FALSE.
     CSymmetry = CSymmetry .OR. CurrentCoordinateSystem() == AxisSymmetric
     CSymmetry = CSymmetry .OR. CurrentCoordinateSystem() == CylindricSymmetric

     FORCE = 0.0D0
     FORCE_im = 0.0D0

     NeedHeat = ANY( NodalTemperature(1:n) /= 0.0d0 )
     !    
     ! Integration stuff:
     ! ------------------  
     NBasis = ntot
     IntegStuff = GaussPoints( element, RelOrder = RelIntegOrder )

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     !
     ! Now we start integrating:
     ! -------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx )
       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
       IF ( NeedHeat ) THEN
         ! Temperature at the integration point:
         !-------------------------------------- 
         Temperature = SUM( NodalTemperature(1:ntot)*Basis(1:ntot) )
 
         ! Heat expansion tensor values at the integration point:
         !-------------------------------------------------------
         HeatExpansion = 0.0d0
         DO i=1,3
           IF ( Isotropic(2) ) THEN
              HeatExpansion(i,i) = SUM( NodalHeatExpansion(1,1,1:n)*Basis(1:n) )
           ELSE
              DO j=1,3
                HeatExpansion(i,j) = SUM( NodalHeatExpansion(i,j,1:n)*Basis(1:n) )
              END DO
           END IF
         END DO
       END IF

       IF ( Isotropic(1) ) Poisson = SUM( Basis(1:n) * NodalPoisson(1:n) )

       C = 0
       IF ( .NOT. Isotropic(1) ) THEN 
          DO i=1,SIZE(ElasticModulus,1)
            DO j=1,SIZE(ElasticModulus,2)
               C(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
            END DO
          END DO
       ELSE
          Young = SUM( Basis(1:n) * ElasticModulus(1,1,1:n) )
       END IF

       SELECT CASE(dim)
       CASE(2)
         IF ( CSymmetry ) THEN
           IF ( Isotropic(1) ) THEN
              C(1,1) = 1.0d0 - Poisson
              C(1,2) = Poisson
              C(1,3) = Poisson
              C(2,1) = Poisson
              C(2,2) = 1.0d0 - Poisson
              C(2,3) = Poisson
              C(3,1) = Poisson
              C(3,2) = Poisson
              C(3,3) = 1.0d0 - Poisson
              C(4,4) = 0.5d0 - Poisson

              C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
           END IF
           Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
           s = s * Radius
         ELSE
           IF ( Isotropic(1) ) THEN
              IF ( PlaneStress ) THEN
                 C(1,1) = 1.0d0
                 C(1,2) = Poisson
                 C(2,1) = Poisson
                 C(2,2) = 1.0d0
                 C(3,3) = 0.5d0*(1-Poisson)
 
                 C = C * Young / ( 1 - Poisson**2 )
              ELSE
                 C(1,1) = 1.0d0 - Poisson
                 C(1,2) = Poisson
                 C(2,1) = Poisson
                 C(2,2) = 1.0d0 - Poisson
                 C(3,3) = 0.5d0 - Poisson
                 C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
              END IF
           ELSE
              IF ( PlaneStress ) THEN
                C(1,1) = C(1,1) - C(1,3)*C(3,1) / C(3,3)
                C(1,2) = C(1,2) - C(1,3)*C(2,3) / C(3,3)
                C(2,1) = C(2,1) - C(1,3)*C(2,3) / C(3,3)
                C(2,2) = C(2,2) - C(2,3)*C(3,2) / C(3,3)
              ELSE
                IF ( NeedHeat ) THEN
                  HeatExpansion(1,1) = HeatExpansion(1,1) + HeatExpansion(3,3) * &
                     ( C(2,2)*C(1,3)-C(1,2)*C(2,3) ) / ( C(1,1)*C(2,2) - C(1,2)*C(2,1) )
  
                  HeatExpansion(2,2) = HeatExpansion(2,2) + HeatExpansion(3,3) * &
                     ( C(1,1)*C(2,3)-C(1,2)*C(1,3) ) / ( C(1,1)*C(2,2) - C(1,2)*C(2,1) )
                END IF
              END IF
              C(3,3) = C(4,4)
              C(1,3) = 0.0d0
              C(3,1) = 0.0d0
              C(2,3) = 0.0d0
              C(3,2) = 0.0d0
              C(4:6,:) = 0.0d0
              C(:,4:6) = 0.0d0
           END IF
         END IF

       CASE(3)
         IF ( Isotropic(1) ) THEN
            C = 0
            C(1,1) = 1.0d0 - Poisson
            C(1,2) = Poisson
            C(1,3) = Poisson
            C(2,1) = Poisson
            C(2,2) = 1.0d0 - Poisson
            C(2,3) = Poisson
            C(3,1) = Poisson
            C(3,2) = Poisson
            C(3,3) = 1.0d0 - Poisson
            C(4,4) = 0.5d0 - Poisson
            C(5,5) = 0.5d0 - Poisson
            C(6,6) = 0.5d0 - Poisson

            C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
!--------------------------------------------------------------------------------
!   Rotate elasticity tensor if required

          ELSE
            IF ( RotateC ) THEN
                CALL RotateElasticityMatrix( C, TransformMatrix, 3 )
            END IF
          END IF

       END SELECT

       DO i=1,6
          StrainLoad(i) = SUM( NodalStrainLoad(i,1:n)*Basis(1:n) )
          StressLoad(i) = SUM( NodalStressLoad(i,1:n)*Basis(1:n) )
       END DO
       StressLoad = StressLoad + MATMUL( C, StrainLoad )

       !
       ! Loop over basis functions (of both unknowns and weights):
       ! ---------------------------------------------------------
       B = 0.0d0
       DO p=1,NBasis
         G = 0.0d0
         SELECT CASE(dim)
         CASE(2)
           IF ( CSymmetry ) THEN
             G(1,1) = dBasisdx(p,1)
             G(1,3) = Basis(p) / Radius
             G(1,4) = dBasisdx(p,2)
             G(2,2) = dBasisdx(p,2)
             G(2,4) = dBasisdx(p,1)
           ELSE
             G(1,1) = dBasisdx(p,1)
             G(1,3) = dBasisdx(p,2)
             G(2,2) = dBasisdx(p,2)
             G(2,3) = dBasisdx(p,1)
           END IF

         CASE(3)
           G(1,1) = dBasisdx(p,1)
           G(2,2) = dBasisdx(p,2)
           G(3,3) = dBasisdx(p,3)
           G(1,4) = dBasisdx(p,2)
           G(2,4) = dBasisdx(p,1)
           G(2,5) = dBasisdx(p,3)
           G(3,5) = dBasisdx(p,2)
           G(1,6) = dBasisdx(p,3)
           G(3,6) = dBasisdx(p,1)
         END SELECT

         LoadatIp = 0.0d0
         LoadatIp_im = 0.0d0
         DO i=1,dim
           DO j=1,6
             LoadAtIp(i) = LoadAtIp(i) + StressLoad(j) * G(i,j)
           END DO
         END DO

         G = MATMUL( G, C )

         !
         ! The (rest of the) righthand side:
         ! ---------------------------------
         DO i=1,dim
           LoadAtIp(i) = LoadAtIp(i) + &
               SUM( LOAD(i,1:n)*Basis(1:n) ) * Basis(p) + &
               SUM( LOAD(4,1:n)*Basis(1:n) ) * dBasisdx(p,i)
           LoadAtIp_im(i) = LoadAtIp_im(i) + &
               SUM( LOAD_im(i,1:n)*Basis(1:n) ) * Basis(p) + &
               SUM( LOAD_im(4,1:n)*Basis(1:n) ) * dBasisdx(p,i)
         END DO

         IF ( NeedHeat ) THEN
           DO i=1,dim
             IF ( CSymmetry ) THEN
               DO j=1,3
                 LoadAtIp(i) = LoadAtIp(i) +  &
                   G(i,j) * HeatExpansion(j,j) * Temperature
               END DO
             ELSE
               DO j=1,dim
                 LoadAtIp(i) = LoadAtIp(i) + &
                   G(i,j) * HeatExpansion(j,j) * Temperature
               END DO
             END IF
           END DO
         END IF

         DO i=1,dim
           FORCE(dim*(p-1)+i) = FORCE(dim*(p-1)+i) + s*LoadAtIp(i)
           FORCE_im(dim*(p-1)+i) = FORCE_im(dim*(p-1)+i) + s*LoadAtIp_im(i)
         END DO
      END DO
    END DO
!------------------------------------------------------------------------------
 END SUBROUTINE StressForceCompose
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE StressBoundary( STIFF,DAMP,FORCE,FORCE_im,LOAD,LOAD_im,NodalSpring, &
      NormalSpring, NodalDamp, NodalBeta,NodalBeta_im,NodalStress,NormalTangential, &
         Element,n,ntot,Nodes )
   USE ElementUtils
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: NodalSpring(:,:,:),NodalDamp(:),NodalBeta(:),LOAD(:,:)
   REAL(KIND=dp) :: LOAD_im(:,:),FORCE_im(:),NodalBeta_im(:)
   TYPE(Element_t),POINTER  :: Element
   TYPE(Nodes_t)    :: Nodes
   REAL(KIND=dp) :: STIFF(:,:),DAMP(:,:),FORCE(:), NodalStress(:,:)

   INTEGER :: n,ntot
   LOGICAL :: NormalTangential, NormalSpring
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(ntot)
   REAL(KIND=dp) :: dBasisdx(ntot,3),detJ

   REAL(KIND=dp) :: u,v,w,s
   REAL(KIND=dp) :: LoadAtIp(3),LoadAtIp_im(3), SpringCoeff(3,3),DampCoeff(3),Beta,Normal(3),&
                    Tangent(3), Tangent2(3), Vect(3), Vect2(3), Stress(3,3), Tf(3,3)
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

   INTEGER :: i,j,k,l,q,p,t,ii,jj,kk,dim,N_Integ, ndim

   LOGICAL :: stat, Csymm, Incompressible

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

   dim = CoordinateSystemDimension()
   Csymm = CurrentCoordinateSystem() == AxisSymmetric .OR. &
           CurrentCoordinateSystem() == CylindricSymmetric

   Incompressible = GetLogical( GetSolverParams(), 'Incompressible', stat )
   IF (Incompressible) THEN
     ndim = dim+1
   ELSE
     ndim = dim
   END IF

   STIFF = 0.0d0
   DAMP  = 0.0d0
   FORCE = 0.0D0
   FORCE_im = 0.0D0
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

     ! Basis function values & derivatives at the integration point:
     !--------------------------------------------------------------
     stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
        Basis, dBasisdx )

     s = detJ * S_Integ(t)
     IF ( Csymm ) s = s * SUM( Nodes % x(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------

     LoadAtIp = 0.0_dp
     LoadAtIp_im = 0.0_dp
     DO i=1,dim
       LoadAtIp(i) = SUM( LOAD(i,1:n)*Basis(1:n) )
       LoadAtIp_im(i) = SUM( LOAD_im(i,1:n)*Basis(1:n) )
     END DO

     Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
     IF ( NormalTangential ) THEN
       LoadAtIp(1) = LoadAtIp(1) + SUM( NodalBeta(1:n)*Basis(1:n) )
       LoadAtIp_im(1) = LoadAtIp_im(1) + SUM( NodalBeta_im(1:n)*Basis(1:n) )
     ELSE
       LoadAtIp = LoadAtIp + SUM( NodalBeta(1:n)*Basis(1:n) ) * Normal
       LoadAtIp_im = LoadAtIp_im + SUM( NodalBeta_im(1:n)*Basis(1:n) ) * Normal
     END IF

     Stress = 0.0_dp
     SELECT CASE(dim)
     CASE(2)
       Stress(1,1) = SUM( NodalStress(1,1:n)*Basis(1:n) )
       Stress(2,2) = SUM( NodalStress(2,1:n)*Basis(1:n) )
       Stress(1,2) = SUM( NodalStress(3,1:n)*Basis(1:n) )
       Stress(2,1) = SUM( NodalStress(3,1:n)*Basis(1:n) )
     CASE(3)
       Stress(1,1) = SUM( NodalStress(1,1:n)*Basis(1:n) )
       Stress(2,2) = SUM( NodalStress(2,1:n)*Basis(1:n) )
       Stress(3,3) = SUM( NodalStress(3,1:n)*Basis(1:n) )
       Stress(1,2) = SUM( NodalStress(4,1:n)*Basis(1:n) )
       Stress(2,1) = SUM( NodalStress(4,1:n)*Basis(1:n) )
       Stress(3,2) = SUM( NodalStress(5,1:n)*Basis(1:n) )
       Stress(2,3) = SUM( NodalStress(5,1:n)*Basis(1:n) )
       Stress(1,3) = SUM( NodalStress(6,1:n)*Basis(1:n) )
       Stress(3,1) = SUM( NodalStress(6,1:n)*Basis(1:n) )
     END SELECT
     LoadAtIp = LoadatIp + MATMUL( Stress, Normal )

     IF ( NormalTangential ) THEN
        Tf=0._dp
        SELECT CASE( Element % TYPE % DIMENSION )
        CASE(1)
           Tangent(1) =  Normal(2)
           Tangent(2) = -Normal(1)
           Tangent(3) =  0.0_dp
           Tangent2   =  0.0_dp
        CASE(2)
           CALL TangentDirections( Normal, Tangent, Tangent2 ) 
        END SELECT
     END IF

     DO i=1,dim
       DO j=1,dim
         SpringCoeff(i,j) = SUM(Basis(1:n)*NodalSpring(1:n,i,j))
       END DO
     END DO

     IF ( NormalTangential ) THEN
       DampCoeff(1) = SUM( NodalDamp(1:n)*Basis(1:n) )
     ELSE
       DampCoeff(1:3) = SUM( NodalDamp(1:n)*Basis(1:n))*Normal
       IF ( NormalSpring ) THEN
         SpringCoeff(1,1) = SpringCoeff(1,1)*Normal(1)
         SpringCoeff(2,2) = SpringCoeff(1,1)*Normal(2)
         SpringCoeff(3,3) = SpringCoeff(1,1)*Normal(3)
       END IF
     END IF

     DO p=1,Ntot
       DO q=1,Ntot
         DO i=1,dim
           IF ( NormalTangential ) THEN
             SELECT CASE(i)
                CASE(1)
                  Vect = Normal
                CASE(2)
                  Vect = Tangent
                CASE(3)
                  Vect = Tangent2
             END SELECT

             DO ii = 1,dim
               DO jj = 1,dim
                  k = (p-1)*ndim + ii
                  l = (q-1)*ndim + jj
                  DAMP(k,l)  = DAMP(k,l) + s * DampCoeff(i) * &
                     Vect(ii) * Vect(jj) * Basis(q) * Basis(p)

                  DO j=1,dim
                    SELECT CASE(j)
                       CASE(1)
                         Vect2 = Normal
                       CASE(2)
                         Vect2 = Tangent
                       CASE(3)
                         Vect2 = Tangent2
                    END SELECT
                    STIFF(k,l) = STIFF(k,l) + s * SpringCoeff(i,j) * &
                       Vect(ii) * Vect2(jj) * Basis(q) * Basis(p)
                  END DO
               END DO
             END DO
           ELSE
              k = (p-1)*ndim + i
              l = (q-1)*ndim + i
              DAMP(k,l)  = DAMP(k,l)  + s * DampCoeff(i) * Basis(q) * Basis(p)

              DO j=1,dim
                l = (q-1)*ndim + j
                STIFF(k,l) = STIFF(k,l) + s * SpringCoeff(i,j) * Basis(q) * Basis(p)
              END DO
           END IF
         END DO
       END DO
     END DO

     DO q=1,Ntot
       DO i=1,dim
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
               k = (q-1)*ndim + j
               FORCE(k) = FORCE(k) + &
                   s * Basis(q) * LoadAtIp(i) * Vect(j)
               FORCE_im(k) = FORCE_im(k) + &
                   s * Basis(q) * LoadAtIp_im(i) * Vect(j)
            END DO
         ELSE
            k = (q-1)*ndim + i
            FORCE(k) = FORCE(k) + s * Basis(q) * LoadAtIp(i)
            FORCE_im(k) = FORCE_im(k) + s * Basis(q) * LoadAtIp_im(i)
         END IF
       END DO
     END DO
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE StressBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE LocalStress( Stress, Strain, PoissonRatio, ElasticModulus, &
      Heatexpansion, NodalTemp, Isotropic, CSymmetry, PlaneStress,     &
      NodalDisp, Basis, dBasisdx, Nodes, dim, n, nBasis, ApplyPressure )
!------------------------------------------------------------------------------
     LOGICAL :: Isotropic(2), CSymmetry, PlaneStress  
     LOGICAL, OPTIONAL :: ApplyPressure
     INTEGER :: n,nd,dim
     INTEGER, OPTIONAL :: nBasis
     TYPE(Nodes_t) :: Nodes
     REAL(KIND=dp) :: Stress(:,:), Strain(:,:), ElasticModulus(:,:,:), &
                      HeatExpansion(:,:,:), NodalTemp(:), Temperature
     REAL(KIND=dp) :: Basis(:), dBasisdx(:,:), PoissonRatio(:), NodalDisp(:,:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,p,q,IND(9),ic
     LOGICAL :: Found, Incompressible
     REAL(KIND=dp) :: C(6,6), Young, LGrad(3,3), Poisson, S(6), &
             Pressure, Radius, HEXP(3,3)
!------------------------------------------------------------------------------

     Incompressible = GetLogical( GetSolverParams(), 'Incompressible', Found )

     Stress = 0.0d0
     Strain = 0.0d0

     nd = n
     IF ( PRESENT( nBasis ) ) nd = nBasis

     ic = dim
     IF ( CSymmetry ) ic=ic+1
!
!    Material parameters:
!    --------------------
     IF ( Isotropic(1) ) Poisson = SUM( Basis(1:n) * PoissonRatio(1:n) )

     C = 0
     IF ( Isotropic(1) ) THEN 
       Young = SUM( Basis(1:n) * ElasticModulus(1,1,1:n) )
     ELSE
       DO i=1,SIZE(ElasticModulus,1)
         DO j=1,SIZE(ElasticModulus,2)
            C(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
         END DO
       END DO
     END IF

     HEXP = 0.0_dp
     IF ( Isotropic(2) ) THEN
        DO i=1,ic
          HEXP(i,i) = SUM( Basis(1:n) * HeatExpansion(1,1,1:n) )
        END DO
     ELSE
        DO i=1,ic
          DO j=1,ic
            HEXP(i,j) = SUM( Basis(1:n) * HeatExpansion(i,j,1:n) )
          END DO
        END DO
     END IF

     Temperature = SUM( Basis(1:n) * NodalTemp(1:n) )

     SELECT CASE(dim)
     CASE(2)
       IF ( CSymmetry ) THEN
         IF ( Isotropic(1) ) THEN
            C(1,1) = 1.0d0 - Poisson
            C(1,2) = Poisson
            C(1,3) = Poisson
            C(2,1) = Poisson
            C(2,2) = 1.0d0 - Poisson
            C(2,3) = Poisson
            C(3,1) = Poisson
            C(3,2) = Poisson
            C(3,3) = 1.0d0 - Poisson
            C(4,4) = 0.5d0 - Poisson

            C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
         END IF
       ELSE
         IF ( Isotropic(1) ) THEN
            IF ( PlaneStress ) THEN
               C(1,1) = 1.0d0
               C(1,2) = Poisson
               C(2,1) = Poisson
               C(2,2) = 1.0d0
               C(3,3) = 0.5d0*(1-Poisson)

               C = C * Young / ( 1 - Poisson**2 )
             ELSE
               C(1,1) = 1.0d0 - Poisson
               C(1,2) = Poisson
               C(2,1) = Poisson
               C(2,2) = 1.0d0 - Poisson
               C(3,3) = 0.5d0 - Poisson

!              To compute Stress_zz afterwards....!
               C(4,1) = Poisson
               C(4,2) = Poisson

               C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
            END IF
         ELSE
            IF ( PlaneStress ) THEN
               C(1,1) = C(1,1) - C(1,3) * C(3,1) / C(3,3)
               C(1,2) = C(1,2) - C(1,3) * C(2,3) / C(3,3)
               C(2,1) = C(2,1) - C(2,3) * C(1,3) / C(3,3)
               C(2,2) = C(2,2) - C(2,3) * C(3,2) / C(3,3)
            ELSE
!              To compute Stress_zz afterwards....!
               C(4,1) = C(3,1)
               C(4,2) = C(3,2)
               C(4,3) = C(3,4)
            END IF
            C(3,3) = C(4,4)
            C(1,3) = 0; C(3,1) = 0
            C(2,3) = 0; C(3,2) = 0
         END IF
       END IF

     CASE(3)
       IF ( Isotropic(1) ) THEN
          C = 0
          C(1,1) = 1.0d0 - Poisson
          C(1,2) = Poisson
          C(1,3) = Poisson
          C(2,1) = Poisson
          C(2,2) = 1.0d0 - Poisson
          C(2,3) = Poisson
          C(3,1) = Poisson
          C(3,2) = Poisson
          C(3,3) = 1.0d0 - Poisson
          C(4,4) = 0.5d0 - Poisson
          C(5,5) = 0.5d0 - Poisson
          C(6,6) = 0.5d0 - Poisson

          C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
       END IF
     END SELECT
!
!    Compute strain: 
!    ---------------
     LGrad = 0._dp
     LGrad(1:dim,1:dim) = MATMUL( NodalDisp(1:dim,1:nd), dBasisdx(1:nd,1:dim) )
     Strain = ( LGrad + TRANSPOSE(LGrad) ) / 2

     IF ( CSymmetry ) THEN
       Strain(1,3) = 0.0d0
       Strain(2,3) = 0.0d0
       Strain(3,1) = 0.0d0
       Strain(3,2) = 0.0d0
       Strain(3,3) = 0.0d0

       Radius = SUM( Nodes % x(1:n) * Basis(1:n) )

       IF ( Radius /= 0.0d0 ) THEN
         Strain(3,3) = SUM(NodalDisp(1,1:nd)*Basis(1:nd))/Radius
       END IF
     END IF

     DO i=1,ic
       Strain(i,i) = Strain(i,i) - HEXP(i,i)*Temperature
     END DO

     !
     ! Compute stresses: 
     ! -----------------
     IF (Incompressible) THEN
       Stress = 2 * Young * Strain / 3
       IF (ApplyPressure) THEN
         Pressure = SUM(NodalDisp(dim+1,1:n)*Basis(1:n))
         DO j=1,dim
           Stress(j,j) = Stress(j,j) - Pressure
         END DO
       END IF
     ELSE
       CALL Strain2Stress( Stress, Strain, C, dim, CSymmetry )
     END IF

     IF ( dim==2 .AND. .NOT. CSymmetry .AND. .NOT. PlaneStress ) THEN
        S(1) = Strain(1,1)
        S(2) = Strain(2,2)
        S(3) = Strain(1,2)
        Stress(3,3) = Stress(3,3) + SUM( C(4,1:3) * S(1:3) )
     END IF
   END SUBROUTINE LocalStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Strain2Stress( Stress, Strain, C, dim, CSymmetry )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Stress(:,:), Strain(:,:), C(:,:)
     INTEGER :: dim
     LOGICAL :: CSymmetry
!------------------------------------------------------------------------------
     INTEGER :: i,j,n,p,q
     INTEGER :: i1(6), i2(6)
     REAL(KIND=dp) :: S(9), csum
!------------------------------------------------------------------------------
     S = 0.0d0
     SELECT CASE(dim)
     CASE(2)
        IF ( CSymmetry ) THEN
          n = 4
          S(1) = Strain(1,1)
          S(2) = Strain(2,2)
          S(3) = Strain(3,3)
          S(4) = Strain(1,2)*2
          i1(1:n) = [ 1,2,3,1 ]
          i2(1:n) = [ 1,2,3,2 ]
        ELSE
          n = 3
          S(1) = Strain(1,1)
          S(2) = Strain(2,2)
          S(3) = Strain(1,2)*2
          i1(1:n) = [ 1,2,1 ]
          i2(1:n) = [ 1,2,2 ]
        END IF
     CASE(3)
        n = 6
        S(1) = Strain(1,1)
        S(2) = Strain(2,2)
        S(3) = Strain(3,3)
        S(4) = Strain(1,2)*2
        S(5) = Strain(2,3)*2
        S(6) = Strain(1,3)*2
        i1(1:n) = [ 1,2,3,1,2,1 ]
        i2(1:n) = [ 1,2,3,2,3,3 ]
     END SELECT

     DO i=1,n
       p = i1(i)
       q = i2(i)
       csum = 0.0d0
       DO j=1,n
          csum = csum + C(i,j) * S(j)
       END DO
       Stress(p,q) = csum
       Stress(q,p) = csum
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Strain2Stress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Tensor26vector( X, V, dim, CSymmetry )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: X(:,:), V(:)
     INTEGER :: dim
     LOGICAL :: CSymmetry
!------------------------------------------------------------------------------
     INTEGER :: i,j,n,p,q
     INTEGER :: i1(6), i2(6)
     REAL(KIND=dp) :: S(9), csum
!------------------------------------------------------------------------------
     S = 0.0d0
     SELECT CASE(dim)
     CASE(2)
        IF ( CSymmetry ) THEN
          n = 4
          i1(1:n) = [ 1,2,3,1 ]
          i2(1:n) = [ 1,2,3,2 ]
        ELSE
          n = 3
          i1(1:n) = [ 1,2,1 ]
          i2(1:n) = [ 1,2,2 ]
        END IF
     CASE(3)
        n = 6
        i1(1:n) = [ 1,2,3,1,2,1 ]
        i2(1:n) = [ 1,2,3,2,3,3 ]
     END SELECT


     V = 0
     DO i=1,n
       p = i1(i)
       q = i2(i)
       V(i) = (X(p,q)+X(q,p))/2
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Tensor26Vector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE InputTensor( Tensor, IsScalar, Name, Material, n, NodeIndexes, Found )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
      LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      INTEGER :: i,j

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0
      IsScalar = .TRUE.

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IF( PRESENT( Found ) ) Found = Stat  
      IF ( .NOT. stat ) RETURN

      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1
      IF ( IsScalar ) THEN
        DO i=1,SIZE(Tensor,1)
          Tensor(i,i,1:n) = Hwrk(1,1,1:n)
        END DO
      ELSE
        IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,MIN(6,SIZE(HWrk,2) )
              Tensor( i,i,1:n ) = Hwrk( 1,1,1:n )
           END DO
        ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(6,SIZE(Hwrk,1))
              Tensor( i,i,1:n ) = Hwrk( i,1,1:n )
           END DO
        ELSE
          DO i=1,MIN(6,SIZE(Hwrk,1))
             DO j=1,MIN(6,SIZE(Hwrk,2))
                Tensor( i,j,1:n ) = Hwrk( i,j,1:n )
             END DO
          END DO
        END IF
      END IF

!------------------------------------------------------------------------------
   END SUBROUTINE InputTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateStressVector(C,T)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: T(:,:), C(:), CT(3,3)
    INTEGER :: i,j,p,q,r,s
    INTEGER :: I1(6) = [ 1,2,3,1,2,1 ], I2(6) = [ 1,2,3,2,3,3 ]

    !
    ! Convert stress vector to stress tensor:
    ! ----------------------------------------
    CT = 0.0d0
    DO i=1,6
      p = I1(i)
      q = I2(i)
      CT(p,q) = C(i)
      CT(q,p) = C(i)
    END DO

    !
    ! Rotate the tensor:
    ! ------------------
    CALL Rotate2IndexTensor( CT, T, 3 )

    !
    ! Convert back to vector form:
    ! ----------------------------
    DO i=1,6
      p = I1(i)
      q = I2(i)
      C(i) = CT(p,q)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateStressVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateStrainVector(C,T)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: T(:,:), C(:), CT(3,3)
    INTEGER :: i,j,p,q,r,s
    INTEGER :: I1(6) = [ 1,2,3,1,2,1 ], I2(6) = [ 1,2,3,2,3,3 ]

    !
    ! Convert strain vector to strain tensor:
    ! ---------------------------------------
    CT = 0.0d0
    C(4:6) = C(4:6)/2
    DO i=1,6
      p = I1(i)
      q = I2(i)
      CT(p,q) = C(i)
      CT(q,p) = C(i)
    END DO

    !
    ! Rotate the tensor:
    ! ------------------
    CALL Rotate2IndexTensor( CT, T, 3 )

    !
    ! Convert back to vector form:
    ! ----------------------------
    DO i=1,6
      p = I1(i)
      q = I2(i)
      C(i) = CT(p,q)
    END DO
    C(4:6) = 2*C(4:6)
!------------------------------------------------------------------------------
  END SUBROUTINE RotateStrainVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateElasticityMatrix(C,T,dim)
!------------------------------------------------------------------------------
    INTEGER :: dim
    REAL(KIND=dp) :: T(:,:), C(:,:)
!------------------------------------------------------------------------------
    SELECT CASE(dim)
    CASE(2)
      CALL RotateElasticityMatrix2D(C,T)
    CASE(3)
      CALL RotateElasticityMatrix3D(C,T)
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE RotateElasticityMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateElasticityMatrix2D(C,T)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: T(:,:), C(:,:), CT(2,2,2,2)
    INTEGER :: i,j,p,q,r,s
    INTEGER :: I1(3) = [ 1,2,1 ], I2(3) = [ 1,2,2 ]

    !
    ! Convert C-matrix to 4 index elasticity tensor:
    ! ----------------------------------------------
    CT = 0.0d0
    DO i=1,2
      p = I1(i)
      q = I2(i)
      DO j=1,2
        r = I1(j)
        s = I2(j)
        CT(p,q,r,s) = C(i,j)
        CT(p,q,s,r) = C(i,j)
        CT(q,p,r,s) = C(i,j)
        CT(q,p,s,r) = C(i,j)
      END DO
    END DO

    !
    ! Rotate the tensor:
    ! ------------------
    CALL Rotate4IndexTensor( CT, T, 2 )

    !
    ! Convert back to matrix form:
    ! ----------------------------
    DO i=1,2
      p = I1(i)
      q = I2(i)
      DO j=1,2
        r = I1(j)
        s = I2(j)
        C(i,j) = CT(p,q,r,s)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateElasticityMatrix2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateElasticityMatrix3D(C,T)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: T(:,:), C(:,:), CT(3,3,3,3)
    INTEGER :: i,j,p,q,r,s
    INTEGER :: I1(6) = [ 1,2,3,1,2,1 ], I2(6) = [ 1,2,3,2,3,3 ]

    !
    ! Convert C-matrix to 4 index elasticity tensor:
    ! ----------------------------------------------
    CT = 0.0d0
    DO i=1,6
      p = I1(i)
      q = I2(i)
      DO j=1,6
        r = I1(j)
        s = I2(j)
        CT(p,q,r,s) = C(i,j)
        CT(p,q,s,r) = C(i,j)
        CT(q,p,r,s) = C(i,j)
        CT(q,p,s,r) = C(i,j)
      END DO
    END DO

    !
    ! Rotate the tensor:
    ! ------------------
    CALL Rotate4IndexTensor( CT, T, 3 )

    !
    ! Convert back to matrix form:
    ! ----------------------------
    DO i=1,6
      p = I1(i)
      q = I2(i)
      DO j=1,6
        r = I1(j)
        s = I2(j)
        C(i,j) = CT(p,q,r,s)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateElasticityMatrix3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Rotate2IndexTensor( C, T, dim )
!------------------------------------------------------------------------------
     INTEGER :: dim
     REAL(KIND=dp) :: C(:,:),T(:,:)
!------------------------------------------------------------------------------
     INTEGER :: i,j
     REAL(KIND=dp) :: C1(dim,dim)
!------------------------------------------------------------------------------
     C1 = 0
     DO i=1,dim
       DO j=1,dim
         C1(:,i) = C1(:,i) + T(i,j)*C(:,j)
       END DO
     END DO

     C = 0
     DO i=1,dim
       DO j=1,dim
         C(i,:) = C(i,:) + T(i,j)*C1(j,:)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Rotate2IndexTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Rotate4IndexTensor( C, T, dim )
!------------------------------------------------------------------------------
     INTEGER :: dim
     REAL(KIND=dp) :: C(:,:,:,:),T(:,:)
!------------------------------------------------------------------------------
     INTEGER :: i,j
     REAL(KIND=dp) :: C1(dim,dim,dim,dim)
!------------------------------------------------------------------------------
     C1 = 0
     DO i=1,dim
       DO j=1,dim
         C1(:,:,:,i) = C1(:,:,:,i) + T(i,j)*C(:,:,:,j)
       END DO
     END DO

     C = 0
     DO i=1,dim
       DO j=1,dim
         C(:,:,i,:) = C(:,:,i,:) + T(i,j)*C1(:,:,j,:)
       END DO
     END DO

     C1 = 0
     DO i=1,dim
       DO j=1,dim
         C1(:,i,:,:) = C1(:,i,:,:) + T(i,j)*C(:,j,:,:)
       END DO
     END DO

     C = 0
     DO i=1,dim
       DO j=1,dim
         C(i,:,:,:) = C(i,:,:,:) + T(i,j)*C1(j,:,:,:)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Rotate4IndexTensor
!------------------------------------------------------------------------------

END MODULE StressLocal

!> \}
