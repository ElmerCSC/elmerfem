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
  USE Materialmodels

  IMPLICIT NONE

  INTEGER, PARAMETER :: VOIGT_I1(6) = [1,2,3,1,2,1], VOIGT_I2(6) = [1,2,3,2,3,3]

!> Tensor index pair (i,j) at 3*(i-1)+j -> the slot it is stored in, for the output
!> layout SymTensorOutputComponents describes. A slot beyond that layout's length is
!> one the layout does not carry.
  INTEGER, PARAMETER :: SYMTENSOR_IND(9) = [ 1,4,6,4,2,5,6,5,3 ]

!> Every field an elasticity solver may recover from the displacement. Both
!> solvers declare exactly this set, so ElasticityStoreEigenmode can walk it
!> without being told which of them are in play: one a solver did not declare is
!> simply not on the mesh, and is skipped.
  CHARACTER(LEN=16), PARAMETER :: STRESS_OUTPUT_FIELDS(7) = &
      [ CHARACTER(LEN=16) :: 'Stress', 'vonMises', 'Principal Stress', 'Strain', &
        'Principal Strain', 'Principal Angle', 'Tresca' ]

!------------------------------------------------------------------------------
!> Persistent state of a nodal projection. Stress and strain fields are recovered
!> from their integration point values by an L2 projection, which needs a solver
!> of its own: a scalar mass matrix, a permutation, and a hidden variable to run
!> the component solves through. All of it survives between calls, so it is kept
!> here rather than in a pile of SAVEd locals.
!------------------------------------------------------------------------------
  TYPE :: NodalProjector_t
    TYPE(Solver_t), POINTER :: PSolver => NULL()
    !> The primary solver this projection serves. Held because an assembly callback
    !> is necessarily file scope -- see NodalProjectorAssemble -- so it cannot reach
    !> the solver by host association and has nowhere else to get it. PSolver is no
    !> substitute: its Variable is the projection's hidden scalar, not the field.
    TYPE(Solver_t), POINTER :: Solver => NULL()
    INTEGER, POINTER :: Perm(:) => NULL()
    CHARACTER(LEN=MAX_NAME_LEN) :: EqName = ' '
    LOGICAL :: UseMask = .FALSE.
    LOGICAL :: Initialized = .FALSE.
    ! Solver state held aside while the projection runs, see NodalProjectorBegin.
    LOGICAL :: LimiterOn = .FALSE., ContactOn = .FALSE.
    LOGICAL :: EigenOn = .FALSE., HarmonicOn = .FALSE., ResidualOn = .FALSE.
    REAL(KIND=dp) :: Relax = 1.0_dp
    LOGICAL :: RelaxFound = .FALSE.
    LOGICAL :: Factorize = .FALSE., FoundFactorize = .FALSE.
    LOGICAL :: FreeFactorize = .FALSE., FoundFreeFactorize = .FALSE.
    LOGICAL :: SkipChange = .FALSE., FoundSkipChange = .FALSE.
  END TYPE NodalProjector_t

!------------------------------------------------------------------------------
!> What a caller of NodalProjectorAssemble supplies: the tensor, or the pair of
!> them, to be fitted at one integration point. This is the only part of a nodal
!> projection that knows any physics -- everything around it is the same whatever
!> is being projected.
!>
!> An implementation of this MUST be a file scope procedure, not an internal one.
!> gfortran carries host association through a trampoline on the stack, which marks
!> the shared object as needing an executable stack and makes the solver module
!> fail to load outright ("cannot enable executable stack as shared object
!> requires"). Check with: readelf -lW <module>.so | grep GNU_STACK -- RW is fine,
!> RWE is this bug. That is why Proj is passed rather than reached: it is the
!> callback's only route to the solver and hence to the material and the solution.
!>
!> Per-element work must not be repeated at every point; do it when t == 1.
!>
!> T2 is untouched unless the caller asked for two fields. Both arrive zeroed, so
!> a component the model does not set stays zero.
!------------------------------------------------------------------------------
  ABSTRACT INTERFACE
    SUBROUTINE ProjectedTensors_i( Proj, Element, Nodes, n, nd, t, Basis, dBasisdx, T1, T2 )
      IMPORT :: dp, Element_t, Nodes_t, NodalProjector_t
      TYPE(NodalProjector_t) :: Proj
      TYPE(Element_t), POINTER :: Element
      TYPE(Nodes_t) :: Nodes
      INTEGER :: n, nd, t
      REAL(KIND=dp) :: Basis(:), dBasisdx(:,:)
      REAL(KIND=dp) :: T1(3,3), T2(3,3)
    END SUBROUTINE ProjectedTensors_i
  END INTERFACE

!------------------------------------------------------------------------------
  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE StressCompose( MASS, DAMP, STIFF, FORCE, FORCE_im, LOAD, LOAD_im, ElasticModulus, &
     NodalPoisson, NodalDensity, PlaneStress, Isotropic,           &
     NodalPreStress, NodalPreStrain, NodalStressLoad, NodalStrainLoad,           &
     NodalHeatExpansion, NodalTemperature, Element, n, ntot, Nodes, RelIntegOrder, StabilityAnalysis, &
     GeometricStiffness, NodalDisplacement, RotateC, TransformMatrix, NodalMeshVelo, &
     NodalDamping, RayleighDamping, RayleighAlpha, RayleighBeta, EvaluateAtIP, EvaluateLoadAtIp, NeedMass)
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
     LOGICAL  :: EvaluateAtIP(3),EvaluateLoadAtIp,NeedMass
 

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t),POINTER :: Element
     INTEGER :: RelIntegOrder

     INTEGER :: n, ntot
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(ntot)
     REAL(KIND=dp) :: dBasisdx(ntot,3),detJ

     REAL(KIND=dp) :: LoadAtIp(3), LoadatIp_im(3), Poisson, Young, Ident(3,3)

     REAL(KIND=dp) :: M(3,3),D(3,3),HeatExpansion(3,3), A(4,4)
     REAL(KIND=dp) :: Temperature,Density, C(6,6), Damping,MeshVelo(3)
     REAL(KIND=dp) :: StressTensor(3,3), StrainTensor(3,3), ElasticStress(3,3), &
                      InnerProd, NodalViscosity(n)
     REAL(KIND=dp) :: StressLoad(6), StrainLoad(6), PreStress(6), PreStrain(6)
     ! The viscoelastic lag stress load, kept apart from StressLoad rather than
     ! written over it. The two are separate contributions to the same force term,
     ! and while they shared one array a Maxwell material silently discarded
     ! whatever "Stress Load" and "Strain Load" had put there: this routine set
     ! StressLoad from those keywords and ViscoElasticLoad then overwrote it,
     ! Tensor26Vector zeroing its output first. Nothing said a keyword was dropped,
     ! and NeedPreStress was forced true immediately after, so the mechanism carried
     ! on with the lag stress alone.
     REAL(KIND=dp) :: VeLoad(6)

     INTEGER :: i,j,k,l,p,q,t,dim,NBasis

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), xPhi, Ux(ntot), Uy(ntot), Uz(ntot)
     REAL(KIND=dp) :: GPA_ip, Load4_ip

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry, NeedHeat, NeedHarmonic, &
         NeedPreStress, ActiveGeometricStiffness, GPA
     TYPE(ValueHandle_t), SAVE :: BetaIP_h, EIP_h, nuIP_h, Load_h(4), Load_h_im(4)

     TYPE(ValueList_t), POINTER :: BF
   
     REAL(KIND=dp) :: GPA_Coeff(n)

     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: ndim, nve
     LOGICAL :: Found, Incompressible,  MaxwellMaterial, FirstTime = .TRUE., ReuseC
     REAL(KIND=dp) :: dt
     REAL(KIND=dp) :: PSOL(4,ntot), SOL(4,ntot), Viscosity, muder0
     CHARACTER :: DimensionString
!------------------------------------------------------------------------------

     TYPE(Variable_t), POINTER, SAVE :: ve_stress

     SAVE FirstTime, dim

     IF (FirstTime) THEN
       dim = CoordinateSystemDimension()
       IF(EvaluateAtIP(1)) &
            CALL ListInitElementKeyword( EIP_h,'Material','Youngs Modulus')
       IF(EvaluateAtIP(2)) &
            CALL ListInitElementKeyword( BetaIP_h,'Material','Heat Expansion Coefficient')
       IF(EvaluateAtIP(3)) &
            CALL ListInitElementKeyword( nuIP_h,'Material','Poisson Ratio')
       IF(EvaluateLoadAtIP) THEN
         DO I=1,DIM
           WRITE(DimensionString,'(I1)') I
           CALL ListInitElementKeyword( Load_h(I),'Body Force','Stress BodyForce '//TRIM(DimensionString))          
           CALL ListInitElementKeyword( Load_h_im(I),'Body Force','Stress BodyForce '//TRIM(DimensionString)//' im')
         END DO
         CALL ListInitElementKeyword( Load_h(4),'Body Force','Stress Pressure')
         CALL ListInitElementKeyword( Load_h_im(4),'Body Force','Stress Pressure im')
       END IF
       FirstTime = .FALSE.
     END IF
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

     IF (NeedMass) &
          NeedMass = ANY( NodalDensity(1:n) /= 0.0d0 )       
     NeedMass = NeedMass .OR. ANY( NodalDamping(1:n) /= 0.0d0 ) .OR. RayleighDamping

     NeedHeat = ANY( NodalTemperature(1:ntot) /= 0.0d0 )
     IF (.NOT. EvaluateLoadAtIP) THEN
       NeedHarmonic = ANY( LOAD_im(:,1:n) /= 0.0d0 )
     ELSE
       NeedHarmonic = .FALSE.
     END IF
     NeedPreStress = ANY( NodalPreStrain(1:6,1:n) /= 0.0d0 ) .OR. &
         ANY( NodalStrainLoad(1:6,1:n) /= 0.0d0 )
     NeedPreStress = NeedPreStress .OR. ANY( NodalPreStress(1:6,1:n) /= 0.0d0 ) .OR. &
         ANY( NodalStressLoad(1:6,1:n) /= 0.0d0 )

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

       IF(.NOT.ASSOCIATED(ve_stress % PrevValues)) THEN
         ALLOCATE(ve_stress % PrevValues(SIZE(ve_stress % values),1))
         ve_stress % PrevValues = 0._dp
       END IF

       ! The lag stress is symmetric, so only its independent components are kept.
       ! StressSolver_Init sizes the variable by the same rule; if the two ever
       ! disagree the indexing below would silently run into the neighbouring
       ! point, so say so instead.
       nve = SymTensorComponents( dim, CSymmetry )
       IF( ve_stress % DOFs /= nve ) THEN
         CALL Fatal( 'StressCompose', 'Variable "ve_stress" has '//I2S(ve_stress % DOFs)// &
             ' components per point, expected '//I2S(nve) )
       END IF

       i = Element % ElementIndex
       j = ve_stress % Perm( i+1 ) - ve_stress % Perm ( i )
       IF( IntegStuff % n /= j ) THEN
         PRINT *,'Inconsistent number of gauss points:',i, IntegStuff % n, j
       END IF

       ! The lag stress converged at the end of the previous timestep is the
       ! reference for the whole of this one, so save it once per element here
       ! instead of retesting the iteration counters at every integration point.
       IF ( GetNonlinIter()==1 .AND. GetCoupledIter()==1 ) THEN
         k = nve * ve_stress % Perm( i )
         l = k + nve * j
         ve_stress % PrevValues(k+1:l,1) = ve_stress % Values(k+1:l)
       END IF

       ! The assembly loop evaluates the isotropic material parameters and builds
       ! the elasticity matrix at each integration point anyway; those can be
       ! handed to LocalStress below rather than evaluated again. Not when the
       ! parameters are given at integration points, though: LocalStress is called
       ! without that context here and would answer differently.
       ReuseC = Isotropic(1) .AND. .NOT. EvaluateAtIP(1) .AND. .NOT. EvaluateAtIP(3)

       NodalViscosity(1:n) = GetReal( GetMaterial(), 'Viscosity', Found )

       SOL = 0; PSOL = 0
       CALL GetVectorLocalSolution( SOL )
       CALL GetVectorLocalSolution( PSOL, tStep=-1 )

       dt = GetTimeStepSize()
       SELECT CASE(dim)
       CASE(1)
         Ux = (SOL(1,1:ntot) - PSOL(1,1:ntot))/dt
         Uy = 0._dp
         Uz = 0._dp
       CASE(2)
         Ux = (SOL(1,1:ntot) - PSOL(1,1:ntot))/dt
         Uy = (SOL(2,1:ntot) - PSOL(2,1:ntot))/dt
         Uz = 0._dp
       CASE(3)
         Ux = (SOL(1,1:ntot) - PSOL(1,1:ntot))/dt
         Uy = (SOL(2,1:ntot) - PSOL(2,1:ntot))/dt
         Uz = (SOL(3,1:ntot) - PSOL(3,1:ntot))/dt
       CASE DEFAULT
         CALL Fatal( 'StressCompose', 'Unkown coordinate system dimension' ) 
       END SELECT

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
       stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------

       Density = SUM( NodalDensity(1:n)*Basis(1:n) )
       IF ( NeedMass ) THEN
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
             IF (EvaluateAtIP(2)) THEN
               HeatExpansion(i,i) = ListGetElementReal( BetaIP_h, Basis, Element, Found, GaussPoint=t)
             ELSE
               HeatExpansion(i,i) = SUM( NodalHeatExpansion(1,1,1:n)*Basis(1:n) )
             END IF
           ELSE
              DO j=1,3
                HeatExpansion(i,j) = SUM( NodalHeatExpansion(i,j,1:n)*Basis(1:n) )
              END DO
           END IF
         END DO
       END IF

       IF ( Isotropic(1) ) THEN
         IF (EvaluateAtIP(3)) THEN
           Poisson = ListGetElementReal( nuIP_h, Basis, Element, Found, GaussPoint=t)
         ELSE
           Poisson = SUM( Basis(1:n) * NodalPoisson(1:n) )
         END IF
       END IF
       
       C = 0
       IF ( .NOT. Isotropic(1) ) THEN 
          DO i=1,SIZE(ElasticModulus,1)
            DO j=1,SIZE(ElasticModulus,2)
               C(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
            END DO
          END DO
       ELSE
         IF (EvaluateAtIP(1)) THEN
           Young = ListGetElementReal( EIP_h, Basis, Element, Found, GaussPoint=t)
         ELSE
           Young = SUM( Basis(1:n) * ElasticModulus(1,1,1:n) )
         END IF
       END IF

       SELECT CASE(dim)
       CASE(2)
         IF ( CSymmetry ) THEN
           IF ( Isotropic(1) ) CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
           Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
           s = s * Radius
         ELSE
           IF ( Isotropic(1) ) THEN
              CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
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
            CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
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

       StressLoad = 0.0d0
       VeLoad = 0.0d0
       IF( NeedPreStress ) THEN
         DO i=1,6
           PreStrain(i) = SUM( NodalPreStrain(i,1:n)*Basis(1:n) )
           PreStress(i) = SUM( NodalPreStress(i,1:n)*Basis(1:n) )
         END DO
         PreStress = PreStress - MATMUL( C, PreStrain  )
         
         DO i=1,6
           StrainLoad(i) = SUM( NodalStrainLoad(i,1:n)*Basis(1:n) )
           StressLoad(i) = SUM( NodalStressLoad(i,1:n)*Basis(1:n) )
         END DO
         StressLoad = MATMUL( C, StrainLoad ) - StressLoad
         
         IF( .NOT. ActiveGeometricStiffness ) THEN 
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
         Viscosity = EffectiveViscosity( Viscosity, Density, Ux, Uy, Uz, &
            Element, Nodes, n, ntot, u, v, w,  muder0, LocalIP=t )

         ! The purely elastic response at this point, which the viscoelastic
         ! update needs alongside the lag stress carried over from the previous
         ! timestep. C, Young and Poisson above refer to this same point.
         IF( ReuseC ) THEN
           CALL LocalStress( ElasticStress,StrainTensor,NodalPoisson,ElasticModulus, &
                NodalHeatExpansion, NodalTemperature, Isotropic,CSymmetry,PlaneStress, &
                SOL, Basis, dBasisdx, Nodes, dim, n, ntot, .FALSE.,                    &
                argC = C, argYoung = Young, argPoisson = Poisson )
         ELSE
           CALL LocalStress( ElasticStress,StrainTensor,NodalPoisson,ElasticModulus, &
                NodalHeatExpansion, NodalTemperature, Isotropic,CSymmetry,PlaneStress, &
                SOL, Basis, dBasisdx, Nodes, dim, n, ntot, .FALSE. )
         END IF

         xPhi = ViscoElasticLoad( ve_stress, t, ElasticStress, VeLoad )
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

       GPA_ip   = MERGE( SUM(GPA_Coeff(1:n)*Basis(1:n)), 0._dp, GPA )
       Load4_ip = SUM( LOAD(4,1:n)*Basis(1:n) )

       DO p=1,NBasis

         CALL BuildGMatrix( G, dBasisdx, Basis, p, Radius, dim, CSymmetry )

         LoadAtIp    = 0.0d0
         LoadAtIp_im = 0.0d0
         IF( NeedPreStress ) THEN
           DO i=1,dim
             DO j=1,6
               LoadAtIp(i) = LoadAtIp(i) + ( StressLoad(j) + VeLoad(j) ) * G(i,j)
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
 
           B = 0.0d0
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
               A(i,dim) = A(i,dim) + GPA_ip*dBasisdx(q,i)*Basis(p)
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
         IF (EvaluateLoadAtIP) THEN
           DO I=1,DIM
             LoadAtIp(I) = LoadAtIp(I) &
                  + ListGetElementReal( Load_h(I), Basis, Element, Found, GaussPoint=t)* Basis(p) &
                  + ListGetElementReal( Load_h(4), Basis, Element, Found, GaussPoint=t)* dBasisdx(p,i)
             LoadAtIp_im(i) = LoadAtIp_im(i) &
                  + ListGetElementReal( Load_h_im(I), Basis, Element, Found, GaussPoint=t)* Basis(p) &
                  + ListGetElementReal( Load_h_im(4), Basis, Element, Found, GaussPoint=t)* dBasisdx(p,i)
             NeedHarmonic = NeedHarmonic .OR. Found
           END DO
         ELSE
           DO i=1,dim
             LoadAtIp(i) = LoadAtIp(i) + &
                  SUM( LOAD(i,1:n)*Basis(1:n) ) * Basis(p) + &
                  Load4_ip * dBasisdx(p,i)
             IF( NeedHarmonic ) THEN
               LoadAtIp_im(i) = LoadAtIp_im(i) + &
                    SUM( LOAD_im(i,1:n)*Basis(1:n) ) * Basis(p) + &
                    SUM( LOAD_im(4,1:n)*Basis(1:n) ) * dBasisdx(p,i)
             END IF
           END DO
         END IF

         IF ( NeedHeat ) THEN
           DO i=1,dim
             DO j=1,MERGE(3, dim, CSymmetry)
               LoadAtIp(i) = LoadAtIp(i) + G(i,j) * HeatExpansion(j,j) * Temperature
             END DO
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
   FUNCTION ViscoElasticLoad(ve_stress, ip, ElasticStress, LagLoad) RESULT(xPhi)
!------------------------------------------------------------------------------
     TYPE(Variable_t) :: ve_stress
     INTEGER :: ip
     REAL(KIND=dp) :: ElasticStress(3,3), LagLoad(6), xPhi
!------------------------------------------------------------------------------
     INTEGER :: i
     REAL(KIND=dp) :: D_new(3,3), PrevD(3,3), VeVec(6), Pres, Pres0, ShearModulus

     i = nve*(ve_stress % perm(Element % ElementIndex) + ip - 1)

     IF(Incompressible) THEN
       ShearModulus = Young / 3
       Pres  = SUM( Basis(1:n) * SOL(ndim,1:n) )
       Pres0 = SUM( Basis(1:n) * PSOL(ndim,1:n) )
     ELSE
       Pres = 0._dp; Pres0 = 0._dp
       ShearModulus = Young / (2*(1+Poisson))
     END IF

     xPhi = 1._dp / ( 1 + ShearModulus / Viscosity * dt )

     ! Lag stress from previous timestep: d = C:u - sigma_VE. It is copied aside
     ! once per element before the integration loop.
     CALL Vector62Tensor( ve_stress % PrevValues(i+1:i+nve,1), PrevD, dim, CSymmetry )

     ! RHS contribution from stored lag stress (no LocalStress call needed):
     StressTensor = xPhi * (PrevD - Pres0 * Ident)
     CALL Tensor26Vector( StressTensor, LagLoad, dim, CSymmetry )

     ! Update lag stress: d_new = (1-xPhi)*C:u + xPhi*(d_prev - p0*I) + p*I
     D_new = (1._dp - xPhi)*ElasticStress + xPhi*(PrevD - Pres0*Ident) + Pres*Ident
     CALL Tensor26Vector( D_new, VeVec, dim, CSymmetry )
     ve_stress % Values(i+1:i+nve) = VeVec(1:nve)
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

     INTEGER :: i,j,k,l,p,q,t,dim,NBasis

     REAL(KIND=dp) :: s,u,v,w, Radius, G(3,6), Load4_ip, Load4im_ip

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

     NeedHeat = ANY( NodalTemperature(1:ntot) /= 0.0d0 )
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
           IF ( Isotropic(1) ) CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
           Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
           s = s * Radius
         ELSE
           IF ( Isotropic(1) ) THEN
              CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
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
            CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
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
       StressLoad = MATMUL( C, StrainLoad ) - StressLoad

       !
       ! Loop over basis functions (of both unknowns and weights):
       ! ---------------------------------------------------------
       Load4_ip   = SUM( LOAD(4,1:n)*Basis(1:n) )
       Load4im_ip = SUM( LOAD_im(4,1:n)*Basis(1:n) )

       DO p=1,NBasis
         CALL BuildGMatrix( G, dBasisdx, Basis, p, Radius, dim, CSymmetry )

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
                Load4_ip * dBasisdx(p,i)
           LoadAtIp_im(i) = LoadAtIp_im(i) + &
                SUM( LOAD_im(i,1:n)*Basis(1:n) ) * Basis(p) + &
                Load4im_ip * dBasisdx(p,i)
         END DO


         IF ( NeedHeat ) THEN
           DO i=1,dim
             DO j=1,MERGE(3, dim, CSymmetry)
               LoadAtIp(i) = LoadAtIp(i) + G(i,j) * HeatExpansion(j,j) * Temperature
             END DO
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
      NormalSpring, NodalDamp, NormalDamp, NodalBeta,NodalBeta_im,NodalStress,NormalTangential, &
         Element,n,ntot,Nodes )
   USE ElementUtils
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: NodalSpring(:,:,:),NodalDamp(:,:,:),NodalBeta(:),LOAD(:,:)
   REAL(KIND=dp) :: LOAD_im(:,:),FORCE_im(:),NodalBeta_im(:)
   TYPE(Element_t), TARGET :: Element
   TYPE(Nodes_t)    :: Nodes
   REAL(KIND=dp) :: STIFF(:,:),DAMP(:,:),FORCE(:), NodalStress(:,:)

   INTEGER :: n,ntot
   LOGICAL :: NormalTangential, NormalSpring, NormalDamp
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(ntot)
   REAL(KIND=dp) :: dBasisdx(ntot,3),detJ

   REAL(KIND=dp) :: u,v,w,s
   REAL(KIND=dp) :: LoadAtIp(3),LoadAtIp_im(3), SpringCoeff(3,3),DampCoeff(3,3),Normal(3),&
                    Tangent(3), Tangent2(3), Vect(3), Vect2(3), Stress(3,3)
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
     stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )

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
         DampCoeff(i,j) = SUM(Basis(1:n)*NodalDamp(1:n,i,j))
       END DO
     END DO

     IF ( .NOT. NormalTangential ) THEN
       IF ( NormalSpring ) THEN
         SpringCoeff(1,1) = SpringCoeff(1,1)*Normal(1)
         SpringCoeff(2,2) = SpringCoeff(1,1)*Normal(2)
         SpringCoeff(3,3) = SpringCoeff(1,1)*Normal(3)
       END IF
       IF( NormalDamp ) THEN
         DampCoeff(1,1) = DampCoeff(1,1)*Normal(1)
         DampCoeff(2,2) = DampCoeff(1,1)*Normal(2)
         DampCoeff(3,3) = DampCoeff(1,1)*Normal(3)
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
                   DAMP(k,l) = DAMP(k,l) + s * DampCoeff(i,j) * &
                       Vect(ii) * Vect2(jj) * Basis(q) * Basis(p)
                 END DO
               END DO
             END DO
           ELSE
             k = (p-1)*ndim + i

             DO j=1,dim
               l = (q-1)*ndim + j
               STIFF(k,l) = STIFF(k,l) + s * SpringCoeff(i,j) * Basis(q) * Basis(p)
               DAMP(k,l) = DAMP(k,l) + s * DampCoeff(i,j) * Basis(q) * Basis(p)
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
      NodalDisp, Basis, dBasisdx, Nodes, dim, n, nBasis, ApplyPressure,&
      argEvaluateAtIP, argEvaluateLoadAtIP, GaussPoint, argC, argYoung,&
      argPoisson)
!------------------------------------------------------------------------------
     LOGICAL :: Isotropic(2), CSymmetry, PlaneStress  
     LOGICAL, OPTIONAL :: ApplyPressure
     INTEGER :: n,nd,dim
     INTEGER, OPTIONAL :: nBasis, GaussPoint
     TYPE(Nodes_t) :: Nodes
     REAL(KIND=dp) :: Stress(:,:), Strain(:,:), ElasticModulus(:,:,:), &
                      HeatExpansion(:,:,:), NodalTemp(:), Temperature
     REAL(KIND=dp) :: Basis(:), dBasisdx(:,:), PoissonRatio(:), NodalDisp(:,:)
     LOGICAL, OPTIONAL :: argEvaluateAtIP(3), argEvaluateLoadAtIP
     REAL(KIND=dp), OPTIONAL :: argC(:,:), argYoung, argPoisson
!------------------------------------------------------------------------------
     INTEGER :: i,j,p,q,ic
     LOGICAL :: Found, Incompressible, FirstTime=.TRUE., PreBuiltC
     ! S went with the plane strain recovery, which is now PlaneStrainStressZZ.
     REAL(KIND=dp) :: C(6,6), Young, LGrad(3,3), Poisson, &
          Pressure, Radius, HEXP(3,3), EzzC(3)
     TYPE(ValueHandle_t), SAVE :: BetaIP_h, EIP_h, nuIP_h, Load_h(4), Load_h_im(4)
     TYPE(Element_t), POINTER :: Element
     CHARACTER :: DimensionString
     LOGICAL :: EvaluateAtIP(3), EvaluateLoadAtIP     
!------------------------------------------------------------------------------

     SAVE FirstTime

     IF(PRESENT(argEvaluateAtIP)) THEN
       EvaluateAtIp = argEvaluateAtIp
     ELSE
       EvaluateAtIp = .FALSE.
     END IF

     IF(PRESENT(argEvaluateLoadAtIP)) THEN
       EvaluateLoadAtIp = argEvaluateLoadAtIp
     ELSE
       EvaluateLoadAtIp = .FALSE.
     END IF
     
     Incompressible = GetLogical( GetSolverParams(), 'Incompressible', Found )

     Element => CurrentModel % CurrentElement
     IF (FirstTime) THEN
       dim = CoordinateSystemDimension()
       IF (PRESENT(argEvaluateAtIP)) THEN
         IF(EvaluateAtIP(1)) &
              CALL ListInitElementKeyword( EIP_h,'Material','Youngs Modulus')
         IF(EvaluateAtIP(2)) &
              CALL ListInitElementKeyword( BetaIP_h,'Material','Heat Expansion Coefficient')
         IF(EvaluateAtIP(3)) &
              CALL ListInitElementKeyword( nuIP_h,'Material','Poisson Ratio')
       END IF
       IF(PRESENT(argEvaluateLoadAtIP) ) THEN
         IF(EvaluateLoadAtIP) THEN
           DO I=1,DIM
             WRITE(DimensionString,'(I1)') I
             CALL ListInitElementKeyword( Load_h(I),'Body Force','Stress BodyForce '//TRIM(DimensionString))          
             CALL ListInitElementKeyword( Load_h_im(I),'Body Force','Stress BodyForce '//TRIM(DimensionString)//' im')
           END DO
           CALL ListInitElementKeyword( Load_h(4),'Body Force','Stress Pressure')
           CALL ListInitElementKeyword( Load_h_im(4),'Body Force','Stress Pressure im')
         END IF
       END IF
       FirstTime = .FALSE.
     END IF
     
     Stress = 0.0d0
     Strain = 0.0d0

     nd = n
     IF ( PRESENT( nBasis ) ) nd = nBasis

     ic = dim
     IF ( CSymmetry ) ic=ic+1
!
!    Material parameters:
!    --------------------
!    The caller may have evaluated the isotropic parameters and built the
!    elasticity matrix at this very integration point already, in which case they
!    are handed over rather than evaluated a second time. Only the isotropic
!    matrix can be passed in: in the anisotropic case the caller has reduced its
!    own copy (rows and columns 4:6 cleared), so the out-of-plane row rebuilt
!    below could no longer be recovered from it.
     PreBuiltC = .FALSE.
     IF ( PRESENT(argC) ) PreBuiltC = Isotropic(1)

     EzzC = 0.0_dp

     IF ( PreBuiltC ) THEN
       C = argC
       Young = argYoung
       Poisson = argPoisson
     ELSE
       IF ( Isotropic(1) ) THEN
         IF (EvaluateAtIP(3)) THEN
           Poisson =  ListGetElementReal(nuIP_h, Basis, Element, Found, GaussPoint=GaussPoint)
         ELSE
           Poisson = SUM( Basis(1:n) * PoissonRatio(1:n) )
         END IF
       END IF

       C = 0
       IF ( Isotropic(1) ) THEN
         IF (EvaluateAtIP(1)) THEN
           Young = ListGetElementReal( EIP_h, Basis, Element, Found, GaussPoint=GaussPoint)
         ELSE
           Young = SUM( Basis(1:n) * ElasticModulus(1,1,1:n) )
         END IF
       ELSE
         DO i=1,SIZE(ElasticModulus,1)
           DO j=1,SIZE(ElasticModulus,2)
              C(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
           END DO
         END DO
       END IF
     END IF

     HEXP = 0.0_dp
     IF ( Isotropic(2) ) THEN
       DO i=1,ic
         IF (EvaluateAtIP(2)) THEN
           HEXP(i,i)= ListGetElementReal( BetaIP_h, Basis, Element, Found, GaussPoint=GaussPoint)
         ELSE
           HEXP(i,i) = SUM( Basis(1:n) * HeatExpansion(1,1,1:n) )
         END IF
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
         IF ( Isotropic(1) .AND. .NOT. PreBuiltC ) &
             CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
       ELSE
         IF ( Isotropic(1) ) THEN
            IF ( .NOT. PreBuiltC ) &
                CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
            IF ( .NOT. PlaneStress ) THEN
              C(4,1) = C(1,2)  ! coefficient for out-of-plane Stress_zz
              C(4,2) = C(1,2)
            END IF
         ELSE
            CALL CondensePlaneElasticityMatrix( C, PlaneStress, EzzC )
         END IF
       END IF

     CASE(3)
       IF ( Isotropic(1) .AND. .NOT. PreBuiltC ) &
           CALL BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
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

!    In two dimensions one out-of-plane component is not carried by the plane
!    system but is still determined by it, and which one depends on the
!    assumption: plane strain leaves Stress_zz to be recovered, plane stress
!    leaves Strain_zz. Neither does any virtual work here -- the plane weak form
!    cannot see them -- but both are reportable, so fill them in for output.
     IF ( dim==2 .AND. .NOT. CSymmetry ) THEN
       IF ( .NOT. PlaneStress ) THEN
         Stress(3,3) = Stress(3,3) + PlaneStrainStressZZ( C, Strain )
       ELSE
         IF ( Isotropic(1) ) THEN
           Strain(3,3) = -Poisson / ( 1.0_dp - Poisson ) * ( Strain(1,1) + Strain(2,2) )
         ELSE
           Strain(3,3) = PlaneStressStrainZZ( EzzC, Strain )
         END IF
       END IF
     END IF
   END SUBROUTINE LocalStress
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
!> Reduce a full 6x6 elasticity matrix to the three-component plane packing that
!> Strain2Stress uses when dim is 2 and CSymmetry is false, namely (11,22,12) with
!> the engineering shear.
!>
!> This is the step that separates a solver which supports anisotropy in the plane
!> from one which does not, and it is nothing to do with the assembly: a plane
!> assembly written for a general dim contracts C(1:3,1:3) with (e11,e22,2e12), so
!> handed a raw 6x6 it would read C(3,3) -- the 33 modulus -- as the shear modulus,
!> and the 11-33 couplings as normal-to-shear ones. Hence the condensation, and
!> hence extracting it here rather than letting a second solver grow its own.
!>
!> WHAT COMES BACK, and it is two different things depending on the assumption,
!> because the plane system determines one out-of-plane component without carrying
!> it:
!>
!>   plane strain -> Stress_zz is left to be recovered, and its coefficients are
!>     stashed in C(4,1:3), to be contracted with (e11,e22,2e12). EzzC is untouched.
!>   plane stress -> Strain_zz is left to be recovered, and its coefficients come
!>     back in EzzC. C(4,1:3) is untouched, and still holds the shear row.
!>
!> The C(4,1:3) stash is not a design so much as an existing convention -- the
!> isotropic path a few lines up in LocalStress fills the same slots by hand -- and
!> it is kept because changing it would alter what every caller reads for no gain
!> here. A caller wanting the plane strain coefficients reads C(4,1:3); one wanting
!> the plane stress coefficients reads EzzC; neither is meaningful under the other
!> assumption.
!>
!> KNOWN LIMITATION, pre-existing and preserved deliberately: the normal-to-shear
!> couplings C(1,4) and C(2,4) of the input are NOT carried into the reduced
!> matrix -- positions (1,3) and (2,3) are zeroed rather than loaded from them. For
!> a material whose axes are the coordinate axes those couplings vanish and nothing
!> is lost, which is why no test has ever noticed. Fixing it would change
!> StressSolve's answers, so it is recorded rather than done.
!------------------------------------------------------------------------------
   SUBROUTINE CondensePlaneElasticityMatrix( C, PlaneStress, EzzC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: C(:,:)     !< 6x6 in; reduced plane packing in C(1:3,1:3) out
     LOGICAL :: PlaneStress
     REAL(KIND=dp) :: EzzC(:)    !< plane stress only: the Strain_zz coefficients
!------------------------------------------------------------------------------
     IF ( PlaneStress ) THEN
!      Coefficients recovering the out-of-plane strain from the in-plane ones, from
!      the plane stress condition Stress_zz = 0. Taken before the condensation
!      below overwrites the out-of-plane row. The third multiplies the engineering
!      shear, which is how C is indexed.
       EzzC(1) = -C(3,1) / C(3,3)
       EzzC(2) = -C(3,2) / C(3,3)
       EzzC(3) = -C(3,4) / C(3,3)

       C(1,1) = C(1,1) - C(1,3) * C(3,1) / C(3,3)
       C(1,2) = C(1,2) - C(1,3) * C(2,3) / C(3,3)
       C(2,1) = C(2,1) - C(2,3) * C(1,3) / C(3,3)
       C(2,2) = C(2,2) - C(2,3) * C(3,2) / C(3,3)
     ELSE
!      To compute Stress_zz afterwards....!
       C(4,1) = C(3,1)
       C(4,2) = C(3,2)
       C(4,3) = C(3,4)
     END IF
     C(3,3) = C(4,4)
     C(1,3) = 0; C(3,1) = 0
     C(2,3) = 0; C(3,2) = 0
!------------------------------------------------------------------------------
   END SUBROUTINE CondensePlaneElasticityMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The out-of-plane normal stress that a plane strain system determines but does
!> not carry, from the coefficients CondensePlaneElasticityMatrix stashed in
!> C(4,1:3). Anisotropic materials only; the isotropic case has a closed form in
!> terms of the Poisson ratio and does not need the stash.
!>
!> Extracted so that the two elasticity solvers cannot drift apart on it. It is one
!> expression and it would be trivial to write twice, which is exactly how a factor
!> of two on the engineering shear below came to live in this file unnoticed until
!> a test was written that could see it.
!------------------------------------------------------------------------------
   PURE FUNCTION PlaneStrainStressZZ( C, Strain ) RESULT( Szz )
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(IN) :: C(:,:), Strain(:,:)
     REAL(KIND=dp) :: Szz
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: S(3)
!------------------------------------------------------------------------------
     S(1) = Strain(1,1)
     S(2) = Strain(2,2)
     ! The engineering shear: C is indexed for it throughout, as Strain2Stress
     ! shows by doubling this same component before contracting with C.
     S(3) = 2.0_dp * Strain(1,2)
     Szz = SUM( C(4,1:3) * S(1:3) )
!------------------------------------------------------------------------------
   END FUNCTION PlaneStrainStressZZ
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The out-of-plane normal strain that a plane stress system determines but does
!> not carry, from the coefficients CondensePlaneElasticityMatrix returned in EzzC.
!> Anisotropic materials only, for the same reason as above.
!------------------------------------------------------------------------------
   PURE FUNCTION PlaneStressStrainZZ( EzzC, Strain ) RESULT( Ezz )
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(IN) :: EzzC(:), Strain(:,:)
     REAL(KIND=dp) :: Ezz
!------------------------------------------------------------------------------
     Ezz = EzzC(1) * Strain(1,1) + EzzC(2) * Strain(2,2) + &
         EzzC(3) * 2.0_dp * Strain(1,2)
!------------------------------------------------------------------------------
   END FUNCTION PlaneStressStrainZZ
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Number of components of a stress or strain field written out for the user, in
!> the order (11,22,33,12) in two dimensions and (11,22,33,12,23,13) in three --
!> the same order truncated, so IND maps into both.
!>
!> This is deliberately NOT SymTensorComponents, which answers a different
!> question. That one gives the length of the vector the weak form contracts, and
!> so drops the out-of-plane component in the plane case, where it does no virtual
!> work. Here the question is what is *reportable*, and the out-of-plane component
!> very much is: in plane strain Stress_zz is nonzero and of the same order as
!> Stress_xx, and in plane stress Strain_zz is. Only the 23 and 13 shears are
!> identically zero in two dimensions, so 6 -> 4 is the whole of the reduction.
!------------------------------------------------------------------------------
   FUNCTION SymTensorOutputComponents( dim ) RESULT( ncomp )
!------------------------------------------------------------------------------
     INTEGER :: dim, ncomp
!------------------------------------------------------------------------------
     ncomp = MERGE( 4, 6, dim <= 2 )
!------------------------------------------------------------------------------
   END FUNCTION SymTensorOutputComponents
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The tensor index pair a stored slot corresponds to, as "11", "23" and so on,
!> for progress messages. Slot 3 is always the out-of-plane component, which under
!> axial symmetry is the hoop.
!------------------------------------------------------------------------------
   FUNCTION SymTensorComponentName( slot ) RESULT( Name )
!------------------------------------------------------------------------------
     INTEGER :: slot
     CHARACTER(LEN=2) :: Name
!------------------------------------------------------------------------------
     CHARACTER(LEN=2), PARAMETER :: Names(6) = [ '11','22','33','12','23','13' ]
!------------------------------------------------------------------------------
     Name = Names(slot)
!------------------------------------------------------------------------------
   END FUNCTION SymTensorComponentName
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The "Exported Variable" definition of a stress or strain field, naming exactly
!> the components SymTensorOutputComponents keeps, in that order. Built here so
!> that the declaration and the assembly that indexes into it cannot drift apart.
!------------------------------------------------------------------------------
   FUNCTION StressFieldDefinition( Name, dim ) RESULT( Str )
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: Name
     INTEGER :: dim
     CHARACTER(LEN=MAX_NAME_LEN) :: Str
!------------------------------------------------------------------------------
     CHARACTER(LEN=2), PARAMETER :: Comp(6) = [ 'xx','yy','zz','xy','yz','xz' ]
     INTEGER :: i, ncomp
!------------------------------------------------------------------------------
     ncomp = SymTensorOutputComponents( dim )

     Str = TRIM(Name)//'['//TRIM(Name)//'_'//Comp(1)//':1'
     DO i=2,ncomp
       Str = TRIM(Str)//' '//TRIM(Name)//'_'//Comp(i)//':1'
     END DO
     Str = TRIM(Str)//']'
!------------------------------------------------------------------------------
   END FUNCTION StressFieldDefinition
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Keep the postprocessed fields of one eigenmode.
!>
!> Each mode has its own displacement and hence its own stresses, so the
!> postprocessing runs once per mode -- but the nodal fields can only ever hold
!> the last one computed. The results therefore go into each field variable's own
!> EigenVectors, mirroring the displacement's, so that a mode written out carries
!> its stresses with it.
!>
!> Call this after the postprocessing of mode Mode has run. The fields are read
!> from the variables themselves, which is where the projection has just left
!> them, so nothing has to be threaded through: whichever fields the solver
!> declared are the ones stored.
!>
!> Imaginary is mandatory rather than OPTIONAL on purpose -- testing an absent
!> OPTIONAL reads uninitialised memory instead of failing to build. A harmonic
!> caller runs each mode twice, .FALSE. for the real part and then .TRUE. for the
!> imaginary; an eigen caller passes .FALSE. once.
!------------------------------------------------------------------------------
   SUBROUTINE ElasticityStoreEigenmode( Solver, Mesh, Mode, Imaginary )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: Mode
     LOGICAL :: Imaginary
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Var, iVar
     INTEGER :: i, k, dofs, nomodes
!------------------------------------------------------------------------------
     nomodes = Solver % NOFEigenValues

     DO i=1,SIZE(STRESS_OUTPUT_FIELDS)
       Var => VariableGet( Mesh % Variables, STRESS_OUTPUT_FIELDS(i) )
       IF ( .NOT. ASSOCIATED( Var ) ) CYCLE
       dofs = Var % DOFs

       IF ( Mode == 1 .AND. .NOT. Imaginary ) THEN
         ! Sized from the field itself rather than from the displacement: the two
         ! agree whenever both are on the same permutation, and where they would
         ! not, this is the length the assignments below actually need.
         IF ( .NOT. ASSOCIATED( Var % EigenVectors ) ) THEN
           ALLOCATE( Var % EigenVectors( nomodes, SIZE( Var % Values ) ) )
           Var % EigenVectors = 0.0_dp
         END IF
         IF ( .NOT. ASSOCIATED( Var % EigenValues ) ) THEN
           ALLOCATE( Var % EigenValues( nomodes ) )
           Var % EigenValues = 0.0_dp
         END IF
         Var % EigenValues = Solver % Variable % EigenValues

         ! The components of a multi-dof field are variables in their own right,
         ! and they are what actually gets written out, so point each at the slice
         ! it names.
         IF ( dofs > 1 ) THEN
           DO k=1,dofs
             iVar => VariableGet( Mesh % Variables, ComponentName( Var % Name, k ) )
             IF ( .NOT. ASSOCIATED( iVar ) ) CALL Fatal( 'ElasticityStoreEigenmode', &
                 'No variable associated: '//TRIM( ComponentName( Var % Name, k ) ) )
             iVar % EigenValues => Var % EigenValues
             iVar % EigenVectors => Var % EigenVectors(:,k::dofs)
           END DO
         END IF
       END IF

       IF ( Imaginary ) THEN
         Var % EigenVectors(Mode,:) = Var % EigenVectors(Mode,:) + &
             CMPLX( 0.0_dp, Var % Values, KIND=dp )
       ELSE
         Var % EigenVectors(Mode,:) = Var % Values
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE ElasticityStoreEigenmode
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Expand one node's stored components, as laid out by SymTensorOutputComponents,
!> back into a full 3x3 tensor. V is that node's slice alone. A slot the layout
!> does not carry is one that is identically zero, so it reads back as zero and
!> the tensor is fully defined either way.
!------------------------------------------------------------------------------
   SUBROUTINE OutputVector2Tensor( V, ncomp, T )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: V(:), T(3,3)
     INTEGER :: ncomp
!------------------------------------------------------------------------------
     INTEGER :: i,j,p,k
!------------------------------------------------------------------------------
     T = 0.0_dp
     p = 0
     DO i=1,3
       DO j=1,3
         p = p + 1
         k = SYMTENSOR_IND(p)
         IF ( k <= ncomp ) T(i,j) = V(k)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE OutputVector2Tensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Add one integration point's contribution to the local mass matrix of a nodal
!> projection. Every component is fitted against this same Galerkin mass, which is
!> why one matrix serves them all and they differ only in the right hand side.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorMass( Mass, Basis, nd, Weight )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Mass(:,:), Basis(:), Weight
     INTEGER :: nd
!------------------------------------------------------------------------------
     INTEGER :: p,q
!------------------------------------------------------------------------------
     DO p=1,nd
       DO q=1,nd
         Mass(p,q) = Mass(p,q) + Weight * Basis(q) * Basis(p)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorMass
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Add one integration point's contribution of a symmetric tensor to the local
!> projection right hand side, packed in the stored component order.
!>
!> Slots the layout does not carry are skipped rather than special cased. In two
!> dimensions those are the 23 and 13 shears, identically zero there; the
!> out-of-plane 33 is carried, and under axial symmetry it is the hoop, which the
!> tensor holds at (3,3) like any other out-of-plane component. That is why this
!> needs no axisymmetric branch, where the callers each used to have one.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorTensor( Force, Basis, nd, ncomp, Weight, T )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Force(:), Basis(:), Weight, T(:,:)
     INTEGER :: nd, ncomp
!------------------------------------------------------------------------------
     INTEGER :: p,i,j,k
!------------------------------------------------------------------------------
     DO p=1,nd
       DO i=1,3
         DO j=i,3
           k = SYMTENSOR_IND( 3*(i-1)+j )
           IF ( k > ncomp ) CYCLE
           Force(ncomp*(p-1)+k) = Force(ncomp*(p-1)+k) + Weight * T(i,j) * Basis(p)
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> As NodalProjectorTensor, but for a source already held as a component vector in
!> the stored order rather than as a tensor -- which is how a UMAT hands its stress
!> back.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorVector( Force, Basis, nd, ncomp, Weight, V )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Force(:), Basis(:), Weight, V(:)
     INTEGER :: nd, ncomp
!------------------------------------------------------------------------------
     INTEGER :: p,i
!------------------------------------------------------------------------------
     DO p=1,nd
       DO i=1,ncomp
         Force(ncomp*(p-1)+i) = Force(ncomp*(p-1)+i) + Weight * V(i) * Basis(p)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Add an element's local projection right hand side into the global one, the
!> components staying interleaved with stride ncomp.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorGlue( ForceG, Force, Perm, Indexes, nd, ncomp )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: ForceG(:), Force(:)
     INTEGER :: Perm(:), Indexes(:), nd, ncomp
!------------------------------------------------------------------------------
     INTEGER :: p,i,l
!------------------------------------------------------------------------------
     DO p=1,nd
       l = Perm(Indexes(p))
       DO i=1,ncomp
         ForceG(ncomp*(l-1)+i) = ForceG(ncomp*(l-1)+i) + Force(ncomp*(p-1)+i)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorGlue
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Number of independent components of a symmetric stress tensor as the weak form
!> sees it, i.e. the length of the vector that Tensor26Vector fills and
!> Vector62Tensor reads back: (11,22,12) in plane, (11,22,33,12) in axisymmetric
!> and (11,22,33,12,23,13) in three dimensions. The plane case has no out-of-plane
!> entry because BuildGMatrix has no column for one -- see
!> SymTensorOutputComponents for the layout used when writing fields out.
!------------------------------------------------------------------------------
   FUNCTION SymTensorComponents( dim, CSymmetry ) RESULT( ncomp )
!------------------------------------------------------------------------------
     INTEGER :: dim, ncomp
     LOGICAL :: CSymmetry
!------------------------------------------------------------------------------
     SELECT CASE( dim )
     CASE( 1 )
       ncomp = 1
     CASE( 2 )
       ncomp = MERGE( 4, 3, CSymmetry )
     CASE DEFAULT
       ncomp = 6
     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION SymTensorComponents
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Expand a symmetric stress vector back to tensor form, the inverse of
!> Tensor26Vector. Components the vector does not carry are zeroed, so the tensor
!> is fully defined whatever the dimension.
!------------------------------------------------------------------------------
   SUBROUTINE Vector62Tensor( V, X, dim, CSymmetry )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: V(:), X(:,:)
     INTEGER :: dim
     LOGICAL :: CSymmetry
!------------------------------------------------------------------------------
     INTEGER :: i,n,p,q
     INTEGER :: i1(6), i2(6)
!------------------------------------------------------------------------------
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
     CASE DEFAULT
        n = 1
        i1(1:n) = [ 1 ]
        i2(1:n) = [ 1 ]
     END SELECT

     X = 0
     DO i=1,n
       p = i1(i)
       q = i2(i)
       X(p,q) = V(i)
       X(q,p) = V(i)
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Vector62Tensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Tensor26vector( X, V, dim, CSymmetry )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: X(:,:), V(:)
     INTEGER :: dim
     LOGICAL :: CSymmetry
!------------------------------------------------------------------------------
     INTEGER :: i,n,p,q
     INTEGER :: i1(6), i2(6)
!------------------------------------------------------------------------------
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
              Tensor( i,i,1:n ) = Hwrk( 1,i,1:n )
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
    INTEGER :: i,p,q

    !
    ! Convert stress vector to stress tensor:
    ! ----------------------------------------
    CT = 0.0d0
    DO i=1,6
      p = VOIGT_I1(i)
      q = VOIGT_I2(i)
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
      p = VOIGT_I1(i)
      q = VOIGT_I2(i)
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
    INTEGER :: i,p,q

    !
    ! Convert strain vector to strain tensor:
    ! ---------------------------------------
    CT = 0.0d0
    C(4:6) = C(4:6)/2
    DO i=1,6
      p = VOIGT_I1(i)
      q = VOIGT_I2(i)
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
      p = VOIGT_I1(i)
      q = VOIGT_I2(i)
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

    !
    ! Convert C-matrix to 4 index elasticity tensor:
    ! ----------------------------------------------
    CT = 0.0d0
    DO i=1,6
      p = VOIGT_I1(i)
      q = VOIGT_I2(i)
      DO j=1,6
        r = VOIGT_I1(j)
        s = VOIGT_I2(j)
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
      p = VOIGT_I1(i)
      q = VOIGT_I2(i)
      DO j=1,6
        r = VOIGT_I1(j)
        s = VOIGT_I2(j)
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

!------------------------------------------------------------------------------
  SUBROUTINE BuildIsotropicC( C, Young, Poisson, dim, CSymmetry, PlaneStress )
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(OUT) :: C(6,6)
    REAL(KIND=dp), INTENT(IN)  :: Young, Poisson
    INTEGER,       INTENT(IN)  :: dim
    LOGICAL,       INTENT(IN)  :: CSymmetry, PlaneStress
!------------------------------------------------------------------------------
    C = 0
    SELECT CASE(dim)
    CASE(2)
      IF ( CSymmetry ) THEN
        C(1,1) = 1.0d0 - Poisson;  C(1,2) = Poisson;       C(1,3) = Poisson
        C(2,1) = Poisson;           C(2,2) = 1.0d0 - Poisson; C(2,3) = Poisson
        C(3,1) = Poisson;           C(3,2) = Poisson;       C(3,3) = 1.0d0 - Poisson
        C(4,4) = 0.5d0 - Poisson
        C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
      ELSE IF ( PlaneStress ) THEN
        C(1,1) = 1.0d0;  C(1,2) = Poisson
        C(2,1) = Poisson; C(2,2) = 1.0d0
        C(3,3) = 0.5d0*(1-Poisson)
        C = C * Young / ( 1 - Poisson**2 )
      ELSE
        C(1,1) = 1.0d0 - Poisson;  C(1,2) = Poisson
        C(2,1) = Poisson;           C(2,2) = 1.0d0 - Poisson
        C(3,3) = 0.5d0 - Poisson
        C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
      END IF
    CASE(3)
      C(1,1) = 1.0d0 - Poisson;  C(1,2) = Poisson;         C(1,3) = Poisson
      C(2,1) = Poisson;           C(2,2) = 1.0d0 - Poisson; C(2,3) = Poisson
      C(3,1) = Poisson;           C(3,2) = Poisson;         C(3,3) = 1.0d0 - Poisson
      C(4,4) = 0.5d0 - Poisson
      C(5,5) = 0.5d0 - Poisson
      C(6,6) = 0.5d0 - Poisson
      C = C * Young / ( (1+Poisson) * (1-2*Poisson) )
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE BuildIsotropicC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BuildGMatrix( G, dBasisdx, Basis, p, Radius, dim, CSymmetry )
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(OUT) :: G(3,6)
    REAL(KIND=dp), INTENT(IN)  :: dBasisdx(:,:), Basis(:), Radius
    INTEGER,       INTENT(IN)  :: p, dim
    LOGICAL,       INTENT(IN)  :: CSymmetry
!------------------------------------------------------------------------------
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
      G(1,1) = dBasisdx(p,1); G(2,2) = dBasisdx(p,2); G(3,3) = dBasisdx(p,3)
      G(1,4) = dBasisdx(p,2); G(2,4) = dBasisdx(p,1)
      G(2,5) = dBasisdx(p,3); G(3,5) = dBasisdx(p,2)
      G(1,6) = dBasisdx(p,3); G(3,6) = dBasisdx(p,1)
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE BuildGMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Build, or rebuild after a mesh change, the auxiliary solver that a nodal
!> projection runs through: a scalar CRS matrix over the same mesh, its own right
!> hand side and communicator, and a hidden variable that the component solves are
!> written into. The keyword namespace is left active on return; the caller clears
!> it when the projection is finished.
!>
!> Rebuilt reports that the matrix was created afresh, so that force vectors the
!> caller keeps between calls have to be reallocated to the new row count.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorSetup( Proj, Solver, NameSpace, MaskKeyword, TempName, &
       GlobalBubbles, Rebuilt, VarPerm, ReuseExisting )
!------------------------------------------------------------------------------
     TYPE(NodalProjector_t) :: Proj
     TYPE(Solver_t), TARGET :: Solver
     CHARACTER(LEN=*) :: NameSpace, MaskKeyword, TempName
     !> Bubble handling, derived by the caller from the element definition. It is
     !> deliberately not settable for the projection alone: the assembly loop walks
     !> the element's DOFs, so a matrix built without the bubbles the element
     !> actually has is glued into out of bounds and segfaults in
     !> CRS_GlueLocalMatrix. Setting "Bubbles in Global System" without a namespace
     !> still works and applies to the primary solver and the projection alike.
     LOGICAL :: GlobalBubbles
     LOGICAL, INTENT(OUT) :: Rebuilt
     !> Permutation recorded on the hidden variable. Defaults to the projection's
     !> own, which is what the component solves index with; callers that register
     !> it against the field permutation instead pass theirs.
     INTEGER, POINTER, OPTIONAL :: VarPerm(:)
     !> Adopt an existing variable of this name, and its permutation, when the mesh
     !> already carries one -- as it can after a refinement has interpolated it
     !> across. Off by default, in which case a fresh one is always made.
     LOGICAL, OPTIONAL :: ReuseExisting
!------------------------------------------------------------------------------
     LOGICAL :: Found, OptimizeBW
     REAL(KIND=dp), POINTER :: TempValues(:)
     INTEGER, POINTER :: PermForVar(:)
!------------------------------------------------------------------------------
     CALL ListSetNameSpace( NameSpace )
     Proj % Solver => Solver

     Rebuilt = ( .NOT. Proj % Initialized ) .OR. Solver % MeshChanged
     IF ( .NOT. Rebuilt ) RETURN

     ! The auxiliary solver object is reassigned rather than freed. The variable
     ! added below records it as its owner and VariableGet dereferences that owner
     ! on its interpolation path, so releasing it would leave the previous mesh's
     ! variable list pointing at freed memory.
     IF ( Proj % Initialized ) THEN
       CALL FreeMatrix( Proj % PSolver % Matrix )
     ELSE
       ALLOCATE( Proj % PSolver )
     END IF
     Proj % PSolver = Solver

     ! An earlier run may have left this variable on the mesh, in which case its
     ! permutation is adopted rather than a new one built.
     Proj % PSolver % Variable => NULL()
     IF ( PRESENT(ReuseExisting) ) THEN
       IF ( ReuseExisting ) Proj % PSolver % Variable => &
           VariableGet( Proj % PSolver % Mesh % Variables, TempName, ThisOnly=.TRUE. )
     END IF

     IF ( ASSOCIATED( Proj % PSolver % Variable ) ) THEN
       Proj % Perm => Proj % PSolver % Variable % Perm
     ELSE
       ALLOCATE( Proj % Perm( SIZE(Solver % Variable % Perm) ) )
       Proj % Perm = 0
     END IF

     OptimizeBW = GetLogical( Proj % PSolver % Values, 'Optimize Bandwidth', Found )
     IF ( .NOT. Found ) OptimizeBW = .TRUE.


     ! Restrict the projection to the bodies that asked for it, when any did.
     IF ( ListGetLogicalAnyEquation( CurrentModel, MaskKeyword ) ) THEN
       Proj % UseMask = .TRUE.
       Proj % EqName = MaskKeyword
     ELSE
       Proj % UseMask = .FALSE.
       Proj % EqName = TRIM( ListGetString( Proj % PSolver % Values,'Equation') )
     END IF

     Proj % PSolver % Matrix => CreateMatrix( CurrentModel, Solver, Solver % Mesh, &
         Proj % Perm, 1, MATRIX_CRS, OptimizeBW, Proj % EqName, GlobalBubbles=GlobalBubbles )

     ALLOCATE( Proj % PSolver % Matrix % RHS(Proj % PSolver % Matrix % NumberOfRows) )
     Proj % PSolver % Matrix % Comm = Solver % Matrix % Comm

     IF ( .NOT. ASSOCIATED( Proj % PSolver % Variable ) ) THEN
       PermForVar => Proj % Perm
       IF ( PRESENT(VarPerm) ) PermForVar => VarPerm
       ALLOCATE( TempValues(Proj % PSolver % Matrix % NumberOfRows) )
       TempValues = 0.0_dp
       CALL VariableAdd( Proj % PSolver % Mesh % Variables, Proj % PSolver % Mesh, &
           Proj % PSolver, TempName, 1, TempValues, PermForVar, Output=.FALSE. )
       Proj % PSolver % Variable => VariableGet( Proj % PSolver % Mesh % Variables, TempName )
     END IF

     Proj % Initialized = .TRUE.
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorSetup
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Hand the solver over to the projection. Several things the primary solve wants
!> are meaningless for an L2 fit and actively harmful if left on -- limiters and
!> contact conditions would be applied to a stress component, an eigen or harmonic
!> setting would send the component solves down the wrong path, and a relaxation
!> factor would damp a solve that is not iterating on anything. Each is put aside
!> here and restored by NodalProjectorEnd.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorBegin( Proj, Solver )
!------------------------------------------------------------------------------
     TYPE(NodalProjector_t) :: Proj
     TYPE(Solver_t), TARGET :: Solver
!------------------------------------------------------------------------------
     LOGICAL :: Found
     TYPE(ValueList_t), POINTER :: Params
!------------------------------------------------------------------------------
     Params => Solver % Values

     Proj % LimiterOn = ListGetLogical( Params, 'Apply Limiter', Found )
     IF ( Proj % LimiterOn ) CALL ListAddLogical( Params, 'Apply Limiter', .FALSE. )

     Proj % ContactOn = ListGetLogical( Params, 'Apply Contact BCs', Found )
     IF ( Proj % ContactOn ) CALL ListAddLogical( Params, 'Apply Contact BCs', .FALSE. )

     Proj % EigenOn = ListGetLogical( Params, 'Eigen Analysis', Found )
     IF ( Proj % EigenOn ) CALL ListAddLogical( Params, 'Eigen Analysis', .FALSE. )

     Proj % HarmonicOn = ListGetLogical( Params, 'Harmonic Analysis', Found )
     IF ( Proj % HarmonicOn ) CALL ListAddLogical( Params, 'Harmonic Analysis', .FALSE. )

     Proj % ResidualOn = ListGetLogical( Params, 'Linear System Residual Mode', Found )
     IF ( Proj % ResidualOn ) &
         CALL ListAddLogical( Params, 'Linear System Residual Mode', .FALSE. )

     Proj % Relax = GetCReal( Params, 'Nonlinear System Relaxation Factor', Proj % RelaxFound )
     IF ( .NOT. Proj % RelaxFound ) Proj % Relax = 1.0_dp
     CALL ListAddConstReal( Params, 'Nonlinear System Relaxation Factor', 1.0_dp )

     ! Every component solve runs against the same matrix and differs only in its
     ! right hand side, so the factorization is formed once and kept for the whole
     ! bracket. And a component solve is not a nonlinear iteration, so it must not
     ! feed the primary solve's convergence measure.
     Proj % Factorize = ListGetLogical( Params, 'Linear System Refactorize', &
         Proj % FoundFactorize )
     Proj % FreeFactorize = ListGetLogical( Params, 'Linear System Free Factorization', &
         Proj % FoundFreeFactorize )
     Proj % SkipChange = ListGetLogical( Params, 'Skip Compute Nonlinear Change', &
         Proj % FoundSkipChange )

     CALL ListAddLogical( Params, 'Linear System Refactorize', .FALSE. )
     CALL ListAddLogical( Params, 'Linear System Free Factorization', .FALSE. )
     CALL ListAddLogical( Params, 'Skip Compute Nonlinear Change', .TRUE. )

     Proj % PSolver % NOFEigenValues = 0
     CurrentModel % Solver => Proj % PSolver
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorBegin
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Assemble a nodal projection: walk the active elements, and at every integration
!> point add the Galerkin mass and the right hand side of whatever tensor the
!> caller supplies. One or two fields, the second being how a caller that wants
!> stress and strain gets both out of a single pass.
!>
!> GetTensors must be a file scope procedure -- see ProjectedTensors_i for why an
!> internal one cannot be used, and what it costs if tried.
!>
!> CSymmetry is passed rather than asked, because the callers do not agree on it:
!> StressSolve counts CylindricSymmetric as axisymmetric and ElasticSolve does not.
!> That disagreement is theirs to keep until someone decides it.
!>
!> The right hand side handed to DefaultUpdateEquations is a zero vector: only the
!> mass matrix is wanted from it, and each component solve overwrites the system's
!> right hand side from ForceG before running.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorAssemble( Proj, ncomp, CSymmetry, GetTensors, ForceG, ForceG2 )
!------------------------------------------------------------------------------
     TYPE(NodalProjector_t) :: Proj
     INTEGER :: ncomp
     LOGICAL :: CSymmetry
     PROCEDURE(ProjectedTensors_i) :: GetTensors
     REAL(KIND=dp) :: ForceG(:)
     REAL(KIND=dp), OPTIONAL :: ForceG2(:)
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PSolver
     TYPE(Element_t), POINTER :: Element
     TYPE(ValueList_t), POINTER :: Equation
     TYPE(Nodes_t), SAVE :: Nodes
     TYPE(GaussIntegrationPoints_t) :: IP
     INTEGER, POINTER :: Indices(:)
     REAL(KIND=dp), ALLOCATABLE :: Mass(:,:), Force(:), SForce(:), Zero(:), &
         Basis(:), dBasisdx(:,:)
     REAL(KIND=dp) :: T1(3,3), T2(3,3), detJ, Weight
     INTEGER :: nmax, elem, n, nd, t
     LOGICAL :: stat, Found, Two
!------------------------------------------------------------------------------
     PSolver => Proj % PSolver
     Two = PRESENT( ForceG2 )

     nmax = PSolver % Mesh % MaxElementDOFs
     ALLOCATE( Indices(nmax), Mass(nmax,nmax), Force(ncomp*nmax), &
         SForce(ncomp*nmax), Zero(nmax), Basis(nmax), dBasisdx(nmax,3) )
     Zero = 0.0_dp

     ForceG = 0.0_dp
     IF ( Two ) ForceG2 = 0.0_dp

     CALL InitializeToZero( PSolver % Matrix, PSolver % Matrix % RHS )
     IF( ALLOCATED(PSolver % Matrix % ConstrainedDOF) ) &
         PSolver % Matrix % ConstrainedDOF = .FALSE.
     IF( ALLOCATED(PSolver % Matrix % Dvalues) ) PSolver % Matrix % Dvalues = 0.0_dp

     DO elem = 1, PSolver % NumberOfActiveElements
       Element => GetActiveElement( elem, PSolver )

       ! Only the bodies that asked for this field, when any did.
       IF ( Proj % UseMask ) THEN
         Equation => GetEquation()
         IF ( .NOT. GetLogical( Equation, Proj % EqName, Found ) ) CYCLE
       END IF

       n = GetElementNOFNodes()
       nd = GetElementDOFs( Indices )
       CALL GetElementNodes( Nodes )

       Mass = 0.0_dp
       Force = 0.0_dp
       SForce = 0.0_dp

       IP = GaussPoints( Element )
       DO t=1,IP % n
         stat = ElementInfo( Element, Nodes, IP % u(t), IP % v(t), IP % w(t), &
             detJ, Basis, dBasisdx )

         Weight = IP % s(t) * detJ
         ! An L2 fit is a volume integral, so in axisymmetric coordinates it is
         ! weighted by the radius like any other. Without this the fit is made in
         ! the Cartesian metric and every varying component comes out wrong.
         IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )

         T1 = 0.0_dp
         T2 = 0.0_dp
         CALL GetTensors( Proj, Element, Nodes, n, nd, t, Basis, dBasisdx, T1, T2 )

         CALL NodalProjectorMass( Mass, Basis, nd, Weight )
         CALL NodalProjectorTensor( Force, Basis, nd, ncomp, Weight, T1 )
         IF ( Two ) CALL NodalProjectorTensor( SForce, Basis, nd, ncomp, Weight, T2 )
       END DO

       CALL DefaultUpdateEquations( Mass, Zero(1:nd) )

       CALL NodalProjectorGlue( ForceG, Force, Proj % Perm, Indices, nd, ncomp )
       IF ( Two ) CALL NodalProjectorGlue( ForceG2, SForce, Proj % Perm, Indices, nd, ncomp )
     END DO

     DEALLOCATE( Indices, Mass, Force, SForce, Zero, Basis, dBasisdx )
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorAssemble
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Run the projection component by component and scatter each result into the
!> nodal field. The components share the matrix and differ only in the right hand
!> side, whose slots are interleaved with stride ncomp; refactorization is already
!> suppressed for the whole bracket by NodalProjectorBegin, so the factorization is
!> formed once and reused.
!>
!> Perm indexes the nodal field, which need not be the projection's own
!> permutation -- it is the field variable's, and the two differ whenever the
!> projection was built over a mask.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorSolve( Proj, FieldName, ncomp, ForceG, Nodal, Perm )
!------------------------------------------------------------------------------
     TYPE(NodalProjector_t) :: Proj
     CHARACTER(LEN=*) :: FieldName
     INTEGER :: ncomp
     REAL(KIND=dp) :: ForceG(:), Nodal(:)
     INTEGER, POINTER :: Perm(:)
!------------------------------------------------------------------------------
     INTEGER :: i, l
     REAL(KIND=dp) :: Norm
     TYPE(Solver_t), POINTER :: PSolver
!------------------------------------------------------------------------------
     PSolver => Proj % PSolver

     DO i=1,ncomp
       CALL Info( 'NodalProjectorSolve', TRIM(FieldName)//' component '// &
           SymTensorComponentName(i), Level=5 )

       PSolver % Matrix % RHS = ForceG(i::ncomp)
       ! Cold start, which is what every caller of this did before it was shared:
       ! the previous component's solution is no kind of guess for this one.
       PSolver % Variable % Values = 0.0_dp

       Norm = DefaultSolve()

       DO l=1,SIZE( Proj % Perm )
         IF ( Proj % Perm(l) <= 0 ) CYCLE
         Nodal(ncomp*(Perm(l)-1)+i) = PSolver % Variable % Values(Proj % Perm(l))
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Give the solver back what NodalProjectorBegin put aside, and drop the keyword
!> namespace the projection was reading through.
!------------------------------------------------------------------------------
   SUBROUTINE NodalProjectorEnd( Proj, Solver )
!------------------------------------------------------------------------------
     TYPE(NodalProjector_t) :: Proj
     TYPE(Solver_t), TARGET :: Solver
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Params
!------------------------------------------------------------------------------
     Params => Solver % Values

     IF ( Proj % EigenOn ) CALL ListAddLogical( Params, 'Eigen Analysis', .TRUE. )
     IF ( Proj % HarmonicOn ) CALL ListAddLogical( Params, 'Harmonic Analysis', .TRUE. )
     ! Put back what was there, or take the keyword away again if it was not:
     ! leaving one behind changes the paths that merely test for its presence.
     IF ( Proj % RelaxFound ) THEN
       CALL ListAddConstReal( Params, 'Nonlinear System Relaxation Factor', Proj % Relax )
     ELSE
       CALL ListRemove( Params, 'Nonlinear System Relaxation Factor' )
     END IF
     IF ( Proj % LimiterOn ) CALL ListAddLogical( Params, 'Apply Limiter', .TRUE. )
     IF ( Proj % ContactOn ) CALL ListAddLogical( Params, 'Apply Contact BCs', .TRUE. )
     IF ( Proj % ResidualOn ) &
         CALL ListAddLogical( Params, 'Linear System Residual Mode', .TRUE. )

     IF ( Proj % FoundFactorize ) THEN
       CALL ListAddLogical( Params, 'Linear System Refactorize', Proj % Factorize )
     ELSE
       CALL ListRemove( Params, 'Linear System Refactorize' )
     END IF
     IF ( Proj % FoundFreeFactorize ) THEN
       CALL ListAddLogical( Params, 'Linear System Free Factorization', Proj % FreeFactorize )
     ELSE
       CALL ListRemove( Params, 'Linear System Free Factorization' )
     END IF
     IF ( Proj % FoundSkipChange ) THEN
       CALL ListAddLogical( Params, 'Skip Compute Nonlinear Change', Proj % SkipChange )
     ELSE
       CALL ListRemove( Params, 'Skip Compute Nonlinear Change' )
     END IF

     CurrentModel % Solver => Solver
     CALL ListSetNameSpace('')
!------------------------------------------------------------------------------
   END SUBROUTINE NodalProjectorEnd
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Trace of a second order tensor over its leading dim by dim block. Was one copy
!> in each solver, differing only in argument names and in spelling REAL(KIND=dp)
!> as DOUBLE PRECISION.
!------------------------------------------------------------------------------
  FUNCTION TRACE( F, dim ) RESULT(t)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F(:,:), t
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER :: i
!------------------------------------------------------------------------------
    t = 0.0_dp
    DO i=1,dim
       t = t + F(i,i)
    END DO
!------------------------------------------------------------------------------
  END FUNCTION TRACE
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Double contraction of two second order tensors over their leading N by N block,
!> A:B. Was one copy internal to StressSolve's interior residual and another named
!> DDOT_PRODUCT in ElasticSolve, identical but for the case of the loop variables.
!------------------------------------------------------------------------------
  FUNCTION DDOTPROD(A,B,N) RESULT(C)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: A(:,:), B(:,:), C
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: i,j
!------------------------------------------------------------------------------
    C = 0.0_dp
    DO i = 1,N
       DO j = 1,N
          C = C + A(i,j)*B(i,j)
       END DO
    END DO
!------------------------------------------------------------------------------
  END FUNCTION DDOTPROD
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Isotropic Lame parameters at a point, from Young's modulus and Poisson ratio.
!------------------------------------------------------------------------------
   SUBROUTINE LameParameters( Young, Poisson, PlaneStress, Lame1, Lame2 )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Young, Poisson, Lame1, Lame2
     LOGICAL :: PlaneStress
!------------------------------------------------------------------------------
     IF ( PlaneStress ) THEN
       Lame1 = Young * Poisson / ( 1.0_dp - Poisson**2 )
     ELSE
       Lame1 = Young * Poisson / ( (1.0_dp + Poisson) * (1.0_dp - 2.0_dp*Poisson) )
     END IF
     Lame2 = Young / ( 2.0_dp * (1.0_dp + Poisson) )
!------------------------------------------------------------------------------
   END SUBROUTINE LameParameters
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The geometrically nonlinear stress an adaptivity estimator wants: Green-Lagrange
!> strain, St Venant-Kirchhoff stress from it, pushed forward with the deformation
!> gradient. Note what comes back is therefore the **first Piola-Kirchhoff** stress,
!> where the small strain path returns Cauchy; the two coincide only as the
!> displacement gradient vanishes. Isotropic only, which is what ElasticSolve's
!> estimators have always been.
!------------------------------------------------------------------------------
   SUBROUTINE LargeDeflectionResidualStress( Stress, Strain, Grad, Lame1, Lame2, dim )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Stress(3,3), Strain(3,3), Grad(3,3), Lame1, Lame2
     INTEGER :: dim
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: DefG(3,3), Stress2(3,3), Identity(3,3), tr
     INTEGER :: i
!------------------------------------------------------------------------------
     Identity = 0.0_dp
     DO i=1,3
       Identity(i,i) = 1.0_dp
     END DO

     DefG = Identity + Grad
     Strain = ( TRANSPOSE(Grad) + Grad + MATMUL(TRANSPOSE(Grad),Grad) ) / 2.0_dp

     tr = 0.0_dp
     DO i=1,dim
       tr = tr + Strain(i,i)
     END DO

     Stress2 = 2.0_dp*Lame2*Strain + Lame1*tr*Identity
     Stress = MATMUL( DefG, Stress2 )
!------------------------------------------------------------------------------
   END SUBROUTINE LargeDeflectionResidualStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adaptivity error estimators for elasticity, shared by both solvers.
!>
!> This is StressSolve's version, which is the more capable of the two that
!> existed: the material arrives through InputTensor as a full tensor with an
!> Isotropic flag so anisotropy is supported, the stress comes from the shared
!> LocalStress rather than being computed inline, and elements are reached through
!> GetElementNOFDOFs/GetElementNodes so p-elements and bubbles are handled.
!>
!> These do not cover every feature of either solver -- they estimate the residual
!> of the plain elasticity operator and predate much of what the two solvers grew
!> afterwards. Sharing them does not change that; it only stops there being two
!> versions of the same partial coverage.
!>
!> Element is deliberately NOT a POINTER dummy here: Adaptive's InsideResidual
!> interface declares it as a plain TYPE(Element_t), and a POINTER dummy would read
!> the target address it is handed as though it were a pointer descriptor. That
!> silently produces a garbage element and a heap overrun rather than a diagnostic.
!------------------------------------------------------------------------------

   SUBROUTINE ElasticityBoundaryResidual( Model, Edge, Mesh, Quant, Perm, Gnorm, Indicator, LargeDeflection )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
     !> Geometrically nonlinear stress, as ElasticSolve's estimators have always
     !> used; false gives the small strain path through LocalStress. Deliberately
     !> not OPTIONAL -- an absent optional cannot be tested without PRESENT, and
     !> getting that wrong reads uninitialised memory rather than failing to build.
     LOGICAL :: LargeDeflection
     REAL(KIND=dp) :: Lame1, Lame2, ResGrad(3,3)

!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     TYPE( Mesh_t )    :: Mesh
     TYPE( Element_t ) :: Edge
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry
     INTEGER :: i,j,k,n,l,t,dim,DOFs,nd,Pn,En
     LOGICAL :: stat, Found
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Normal(3), EdgeLength
     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), dEdgeBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:), ExtPressure(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Force(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: ElasticModulus(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: LocalTemp(:), LocalHexp(:,:,:)

     REAL(KIND=dp) :: Residual(3), ResidualNorm, Area
     REAL(KIND=dp) :: ForceSolved(3), Dir(3)
     REAL(KIND=dp) :: Displacement(3)
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Grad(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)
     REAL(KIND=dp) :: Identity(3,3), YoungsAverage

     LOGICAL :: PlaneStress, Isotropic(2)=.TRUE., CSymmetry = .FALSE.
     TYPE(ValueList_t), POINTER :: Material, Equation, BodyForce, BC
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     SAVE Nodes, EdgeNodes
!------------------------------------------------------------------------------

     ! Initialize:
     ! -----------
     Gnorm = 0.0d0
     Indicator = 0.0d0

     Identity = 0.0d0
     DO i=1,3
        Identity(i,i) = 1.0d0
     END DO

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric

     dim = CoordinateSystemDimension()
     DOFs = dim

!    --------------------------------------------------
     Element => Edge % BoundaryInfo % Left

     IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     En = GetElementNOFNodes( Edge )
     CALL GetElementNodes( EdgeNodes )

     nd = GetElementNOFDOFs( Element )
     Pn = GetElementNOFNodes( Element )
     CALL GetElementNodes( Nodes, UElement=Element )

     ALLOCATE( EdgeBasis(En), dEdgeBasisdx(En,3), x(En), y(En), z(En), &
        ExtPressure(En), Basis(nd), dBasisdx(nd,3), Force(3,En), &
        NodalDisplacement(3,nd), ElasticModulus(6,6,Pn),&
        NodalPoissonRatio(Pn), LocalTemp(nd), LocalHexp(3,3,Pn) )

     LocalTemp = 0
     LocalHexp = 0

     DO l = 1,En
       DO k = 1,Pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO

     ! Integrate square of residual over boundary element:
     ! ---------------------------------------------------
     Indicator     = 0.0d0
     EdgeLength    = 0.0d0
     YoungsAverage = 0.0d0
     ResidualNorm  = 0.0d0

     BC => GetBC()
     IF ( .NOT.ASSOCIATED( BC ) ) RETURN

     ! Logical parameters:
     ! -------------------
     Equation => GetEquation( Element )
     PlaneStress = GetLogical( Equation, 'Plane Stress' ,Found )

     Material => GetMaterial( Element )
     NodalPoissonRatio(1:pn) = GetReal( &
                  Material, 'Poisson Ratio',Found, Element )
     CALL InputTensor( ElasticModulus, Isotropic(1), &
                 'Youngs Modulus', Material, Pn, Element % NodeIndexes )

     ! Given traction:
     ! ---------------
     Force = 0.0d0
     Force(1,1:En) = GetReal( BC, 'Force 1', Found )
     Force(2,1:En) = GetReal( BC, 'Force 2', Found )
     Force(3,1:En) = GetReal( BC, 'Force 3', Found )

     ! Force in normal direction:
     ! ---------------------------
     ExtPressure(1:En) = GetReal( BC, 'Normal Force', Found )

     ! If dirichlet BC for displacement in any direction given,
     ! nullify force in that direction:
     ! --------------------------------------------------------
     Dir = 1.0d0
     IF ( ListCheckPresent( BC, 'Displacement' ) )   Dir = 0
     IF ( ListCheckPresent( BC, 'Displacement 1' ) ) Dir(1) = 0
     IF ( ListCheckPresent( BC, 'Displacement 2' ) ) Dir(2) = 0
     IF ( ListCheckPresent( BC, 'Displacement 3' ) ) Dir(3) = 0

     ! Elementwise nodal solution:
     ! ---------------------------
     CALL GetVectorLocalSolution( NodalDisplacement, UElement=Element )

     ! Integration:
     ! ------------
     EdgeLength    = 0.0d0
     YoungsAverage = 0.0d0
     ResidualNorm  = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
            EdgeBasis, dEdgeBasisdx )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
   
           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )

           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

        u = SUM( EdgeBasis(1:En) * x(1:En) )
        v = SUM( EdgeBasis(1:En) * y(1:En) )
        w = SUM( EdgeBasis(1:En) * z(1:En) )

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
           Basis, dBasisdx )

        ! Stress tensor on the edge:
        ! --------------------------
        IF ( LargeDeflection ) THEN
          CALL LameParameters( SUM( ElasticModulus(1,1,1:pn) * Basis(1:pn) ), &
              SUM( NodalPoissonRatio(1:pn) * Basis(1:pn) ), PlaneStress, Lame1, Lame2 )
          ResGrad = MATMUL( NodalDisplacement(:,1:nd), dBasisdx(1:nd,:) )
          CALL LargeDeflectionResidualStress( Stress1, Strain, ResGrad, Lame1, Lame2, dim )
        ELSE
        CALL LocalStress( Stress1, Strain, NodalPoissonRatio, &
           ElasticModulus, LocalHExp, LocalTemp, &
           Isotropic, CSymmetry, PlaneStress, &
           NodalDisplacement, Basis, dBasisdx, Nodes, dim, pn, nd )
        END IF

        ! Given force at the integration point:
        ! -------------------------------------
        Residual = MATMUL( Force(:,1:En), EdgeBasis(1:En) ) - &
          SUM( ExtPressure(1:En) * EdgeBasis(1:En) ) * Normal

        ForceSolved = MATMUL( Stress1, Normal )
        Residual = Residual - ForceSolved * Dir

        EdgeLength    = EdgeLength + s
        ResidualNorm  = ResidualNorm  + s * SUM(Residual(1:DIM) ** 2)
        YoungsAverage = YoungsAverage + &
                    s * SUM( ElasticModulus(1,1,1:Pn) * Basis(1:Pn) )
     END DO

     IF ( YoungsAverage > AEPS ) THEN
        YoungsAverage = YoungsAverage / EdgeLength
        Indicator = EdgeLength * ResidualNorm / YoungsAverage
     END IF

     DEALLOCATE( EdgeBasis, dEdgeBasisdx, x, y, z, ExtPressure, Basis, &
      dBasisdx, Force, NodalDisplacement, ElasticModulus, NodalPoissonRatio, &
      LocalTemp, LocalHexp )
!------------------------------------------------------------------------------
   END SUBROUTINE ElasticityBoundaryResidual

  SUBROUTINE ElasticityEdgeResidual( Model,Edge,Mesh,Quant,Perm, Indicator, LargeDeflection )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
     !> Geometrically nonlinear stress, as ElasticSolve's estimators have always
     !> used; false gives the small strain path through LocalStress. Deliberately
     !> not OPTIONAL -- an absent optional cannot be tested without PRESENT, and
     !> getting that wrong reads uninitialised memory rather than failing to build.
     LOGICAL :: LargeDeflection
     REAL(KIND=dp) :: Lame1, Lame2, ResGrad(3,3)


     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE( Mesh_t )    :: Mesh
     TYPE( Element_t ) :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,l,n,t,dim,DOFs,En,Pn, nd
     LOGICAL :: stat, Found

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Stressi(3,3,2), Jump(3), Identity(3,3)
     REAL(KIND=dp) :: Normal(3)
     REAL(KIND=dp) :: Displacement(3)
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: YoungsAverage
     REAL(KIND=dp) :: Grad(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)

     REAL(KIND=dp), ALLOCATABLE :: LocalTemp(:), LocalHexp(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: ElasticModulus(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), Basis(:), dBasisdx(:,:)

     LOGICAL :: PlaneStress, Isotropic(2)=.TRUE., CSymmetry

     TYPE(ValueList_t), POINTER :: Material, Equation

     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp) :: Residual, ResidualNorm, EdgeLength

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     SAVE Nodes, EdgeNodes
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     dim = CoordinateSystemDimension()
     DOFs = dim

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric


     Identity = 0.0d0
     Metric   = 0.0d0
     DO i = 1,3
        Metric(i,i)   = 1.0d0
        Identity(i,i) = 1.0d0
     END DO
!
!    ---------------------------------------------
     En = GetElementNOFNodes( Edge )
     CALL GetElementNodes( EdgeNodes, Edge )

     Element => Edge % BoundaryInfo % Left
     pn = GetElementNOFNodes( Element )
     nd = GetElementNOFDOFs( Element )

     Element => Edge % BoundaryInfo % Right
     nd = MAX( nd, GetElementNOFDOFs( Element ) )
     pn = MAX( pn, GetElementNOFNodes( Element ) )

     ALLOCATE( LocalTemp(nd), LocalHexp(3,3,Pn), x(En), y(En), z(En), &
      NodalDisplacement(3,nd), ElasticModulus(6,6,pn), &
      NodalPoissonRatio(pn), EdgeBasis(En), Basis(nd), dBasisdx(nd,3) )

     LocalTemp = 0
     LocalHexp = 0

!    Integrate square of jump over edge:
!    ------------------------------------
     ResidualNorm  = 0.0d0
     EdgeLength    = 0.0d0
     Indicator     = 0.0d0
     Grad          = 0.0d0
     YoungsAverage = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Stressi = 0.0d0
        DO i = 1,2
           IF ( i==1 ) THEN
              Element => Edge % BoundaryInfo % Left
           ELSE
              Element => Edge % BoundaryInfo % Right
           END IF

           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE

           pn = GetElementNOFNodes( Element )
           nd = GetElementNOFDOFs( Element )
           CALL GetElementNodes( Nodes, Element )
           DO j = 1,en
              DO k = 1,pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
               Basis, dBasisdx )

           ! Logical parameters:
           ! -------------------
           Equation => GetEquation( Element )
           PlaneStress = GetLogical( Equation,'Plane Stress',Found )

           ! Material parameters:
           ! --------------------
           Material => GetMaterial( Element )
           NodalPoissonRatio(1:pn) = GetReal( Material, 'Poisson Ratio', Found, Element )
           CALL InputTensor( ElasticModulus, Isotropic(1), &
                         'Youngs Modulus', Material, pn, Element % NodeIndexes )

           ! Elementwise nodal solution:
           ! ---------------------------
           CALL GetVectorLocalSolution( NodalDisplacement, UElement=Element )

           ! Stress tensor on the edge:
           ! --------------------------
        IF ( LargeDeflection ) THEN
          CALL LameParameters( SUM( ElasticModulus(1,1,1:pn) * Basis(1:pn) ), &
              SUM( NodalPoissonRatio(1:pn) * Basis(1:pn) ), PlaneStress, Lame1, Lame2 )
          ResGrad = MATMUL( NodalDisplacement(:,1:nd), dBasisdx(1:nd,:) )
          CALL LargeDeflectionResidualStress( Stress1, Strain, ResGrad, Lame1, Lame2, dim )
        ELSE
           CALL LocalStress( Stress1, Strain, NodalPoissonRatio, &
              ElasticModulus, LocalHExp, LocalTemp, Isotropic, CSymmetry, PlaneStress, &
              NodalDisplacement, Basis, dBasisdx, Nodes, dim, pn, nd )
        END IF

           Stressi(:,:,i) = Stress1
        END DO

        EdgeLength  = EdgeLength + s
        Jump = MATMUL( ( Stressi(:,:,1) - Stressi(:,:,2)), Normal )
        ResidualNorm = ResidualNorm + s * SUM( Jump(1:DIM) ** 2 )

        YoungsAverage = YoungsAverage + s *  &
                    SUM( ElasticModulus(1,1,1:pn) * Basis(1:pn) )
     END DO

     YoungsAverage = YoungsAverage / EdgeLength
     Indicator = EdgeLength * ResidualNorm / YoungsAverage

     DEALLOCATE( LocalTemp, LocalHexp, x, y, z, NodalDisplacement, &
       ElasticModulus, NodalPoissonRatio, EdgeBasis, Basis, dBasisdx )
!------------------------------------------------------------------------------
   END SUBROUTINE ElasticityEdgeResidual

   SUBROUTINE ElasticityInsideResidual( Model, Element,  &
                      Mesh, Quant, Perm, Fnorm, Indicator, LargeDeflection )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
     !> Geometrically nonlinear stress, as ElasticSolve's estimators have always
     !> used; false gives the small strain path through LocalStress. Deliberately
     !> not OPTIONAL -- an absent optional cannot be tested without PRESENT, and
     !> getting that wrong reads uninitialised memory rather than failing to build.
     LOGICAL :: LargeDeflection
     REAL(KIND=dp) :: Lame1, Lame2, ResGrad(3,3)

!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE( Mesh_t )    :: Mesh
     TYPE( Element_t ) :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,m,n,nd,t,dim,DOFs,I1(6),I2(6)
     INTEGER, ALLOCATABLE :: Indexes(:)

     LOGICAL :: stat, Found

     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp) :: Density
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Damping
     REAL(KIND=dp) :: Displacement(3),Identity(3,3), YoungsAverage
     REAL(KIND=dp) :: Grad(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)
     REAL(KIND=dp) :: Energy

     REAL(KIND=dp), ALLOCATABLE :: ElasticModulus(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDamping(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: LocalHexp(:,:,:), vec(:)
     REAL(KIND=dp), ALLOCATABLE :: Stressi(:,:,:), LocalTemp(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalForce(:,:), Veloc(:,:), Accel(:,:)

     LOGICAL :: PlaneStress, CSymmetry, Isotropic(2)=.TRUE., Transient

     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual(3), ResidualNorm, Area

     TYPE(ValueList_t), POINTER :: Material, BodyForce, Equation

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     SAVE Nodes
!------------------------------------------------------------------------------
     ! Initialize:
     ! -----------
     Fnorm     = 0.0d0
     Indicator = 0.0d0

     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     dim = CoordinateSystemDimension()
     DOFs = dim 

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric

     ! Element nodal points:
     ! ---------------------
     nd = GetElementNOFDOFs()
     n  = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )

     ALLOCATE( ElasticModulus(6,6,nd), NodalDensity(n), NodalPoissonRatio(n), &
         NodalDamping(n), NodalDisplacement(3,nd), LocalHExp(3,3,n), vec(nd), &
         Stressi(3,3,nd), LocalTemp(nd), Basis(nd), dBasisdx(nd,3), &
         NodalForce(4,n), Veloc(3,nd), Accel(3,nd) )

     LocalTemp = 0
     LocalHexp = 0

     ! Logical parameters:
     ! -------------------
     equation => GetEquation()
     PlaneStress = GetLogical( Equation, 'Plane Stress',Found )

     ! Material parameters:
     ! --------------------
     Material => GetMaterial()

     CALL InputTensor( ElasticModulus, Isotropic(1), &
           'Youngs Modulus', Material, n, Element % NodeIndexes )

     NodalPoissonRatio(1:n) = GetReal( Material, 'Poisson Ratio', Found )

     ! Check for time dep.
     ! -------------------
     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Transient = .TRUE.
        Var => VariableGet( Model % Variables, 'Displacement', .TRUE. )

        nd = GetElementDOFs( Indexes )

        Veloc = 0.0d0
        Accel = 0.0d0
        DO i=1,DOFs
           Veloc(i,1:nd) = Var % PrevValues(DOFs*(Var % Perm(Indexes(1:nd))-1)+i,1)
           Accel(i,1:nd) = Var % PrevValues(DOFs*(Var % Perm(Indexes(1:nd))-1)+i,2)
        END DO
        NodalDensity(1:n) = GetReal( Material, 'Density', Found )
        NodalDamping(1:n) = GetReal( Material, 'Damping', Found )
     ELSE
        Transient = .FALSE.
     END IF

     ! Elementwise nodal solution:
     ! ---------------------------
     CALL GetVectorLocalSolution( NodalDisplacement )

     ! Body Forces:
     ! ------------
     BodyForce => GetBodyForce()

     NodalForce = 0.0d0

     IF ( ASSOCIATED( BodyForce ) ) THEN
        NodalForce(1,1:n) = NodalForce(1,1:n) + GetReal( &
            BodyForce, 'Stress BodyForce 1', Found )
        NodalForce(2,1:n) = NodalForce(1,1:n) + GetReal( &
            BodyForce, 'Stress BodyForce 2', Found )
        NodalForce(3,1:n) = NodalForce(1,1:n) + GetReal( &
            BodyForce, 'Stress BodyForce 3', Found )
     END IF

     Identity = 0.0D0
     DO i = 1,dim
        Identity(i,i) = 1.0D0
     END DO
     CSymmetry = .FALSE.

     Var => VariableGet( Model % Variables, 'Stress 1' )
     IF ( ASSOCIATED( Var ) ) THEN

       ! If stress already computed:
       ! ---------------------------
       I1(1:6) = (/ 1,2,3,1,2,1 /)
       I2(1:6) = (/ 1,2,3,2,3,3 /)
       DO i=1,6
         CALL GetScalarLocalSolution(Vec(1:nd),'Stress ' // CHAR(i+ICHAR('0')))
         Stressi(I1(i),I2(i),1:nd) = Vec(1:nd)
         Stressi(I2(i),I1(i),1:nd) = Vec(1:nd)
       END DO
     ELSE
       ! Values of the stress tensor at node points:
       ! -------------------------------------------
       DO i = 1,n
         u = Element % TYPE % NodeU(i)
         v = Element % TYPE % NodeV(i)
         w = Element % TYPE % NodeW(i)

         stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )

        IF ( LargeDeflection ) THEN
          CALL LameParameters( SUM( ElasticModulus(1,1,1:n) * Basis(1:n) ), &
              SUM( NodalPoissonRatio(1:n) * Basis(1:n) ), PlaneStress, Lame1, Lame2 )
          ResGrad = MATMUL( NodalDisplacement(:,1:nd), dBasisdx(1:nd,:) )
          CALL LargeDeflectionResidualStress( Stressi(:,:,i), Strain, ResGrad, Lame1, Lame2, dim )
        ELSE
         CALL LocalStress( Stressi(:,:,i), Strain, NodalPoissonRatio, &
                   ElasticModulus, LocalHExp, LocalTemp, Isotropic, CSymmetry, PlaneStress, &
                   NodalDisplacement, Basis, dBasisdx, Nodes, dim, n, nd )
        END IF
       END DO
     END IF

     ! Integrate square of residual over element:
     ! ------------------------------------------
     ResidualNorm = 0.0d0
     Fnorm = 0.0d0
     Area = 0.0d0
     Energy = 0.0d0
     YoungsAverage = 0.0d0

     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,u,v,w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! Residual of the diff.equation:
        ! ------------------------------
        Residual = 0.0d0
        DO i = 1,3
           Residual(i) = -SUM( NodalForce(i,1:n) * Basis(1:n) )

           IF ( Transient ) THEN
              Residual(i) = Residual(i) + SUM(NodalDensity(1:n)*Basis(1:n)) * &
                            SUM( Accel(i,1:nd) * Basis(1:nd) )
              Residual(i) = Residual(i) + SUM(NodalDamping(1:n)*Basis(1:n)) * &
                            SUM( Veloc(i,1:nd) * Basis(1:nd) )
           END IF

           DO j = 1,3
             Residual(i) = Residual(i) - SUM(Stressi(i,j,1:nd)*dBasisdx(1:nd,j))
           END DO
        END DO

!       IF ( CSymmetry ) THEN
!          DO k=1,3
!             Residual(1) = Residual(1) + ...
!          END DO
!       END IF

       ! Dual norm of the load:
       ! ----------------------
        DO i = 1,dim
           Fnorm = Fnorm + s * SUM( NodalForce(i,1:n) * Basis(1:n) ) ** 2
        END DO

        YoungsAverage = YoungsAverage + s*SUM( ElasticModulus(1,1,1:n) * Basis(1:n) )

        ! Energy:
        ! -------
        IF ( LargeDeflection ) THEN
          CALL LameParameters( SUM( ElasticModulus(1,1,1:n) * Basis(1:n) ), &
              SUM( NodalPoissonRatio(1:n) * Basis(1:n) ), PlaneStress, Lame1, Lame2 )
          ResGrad = MATMUL( NodalDisplacement(:,1:nd), dBasisdx(1:nd,:) )
          CALL LargeDeflectionResidualStress( Stress1, Strain, ResGrad, Lame1, Lame2, dim )
        ELSE
        CALL LocalStress( Stress1, Strain, NodalPoissonRatio, &
           ElasticModulus, LocalHExp, LocalTemp, Isotropic, CSymmetry, PlaneStress, &
           NodalDisplacement, Basis, dBasisdx, Nodes, dim, n, nd )
        END IF

        Energy = Energy + s*DDOTPROD(Strain,Stress1,dim) / 2.0d0

        Area = Area + s
        ResidualNorm = ResidualNorm + s * SUM( Residual(1:dim) ** 2 )
     END DO

     YoungsAverage = YoungsAverage / Area
     Fnorm = Energy
     Indicator = Area * ResidualNorm / YoungsAverage
 
     DEALLOCATE( ElasticModulus, NodalDensity, NodalPoissonRatio,  &
         NodalDamping, NodalDisplacement, LocalHExp, vec, Stressi, &
         LocalTemp, Basis, dBasisdx, NodalForce, Veloc, Accel )


!------------------------------------------------------------------------------
   END SUBROUTINE ElasticityInsideResidual

END MODULE StressLocal

!> \}
