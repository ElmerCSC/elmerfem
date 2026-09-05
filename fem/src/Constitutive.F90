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
! *  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA  02110-1301  USA
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> One constitutive law at one integration point: given the deformation state,
!> the environment and the history there, return the stress -- and, where a model
!> can form one, the consistent tangent.
!>
!> The shape follows UMAT's semantics without its thirty-five positional
!> arguments, so that a model ignores what it does not need and fields can be
!> added without touching every implementation.
!>
!> THE TWO ENTRY POINTS ARE NOT INTERCHANGEABLE, and which one a driver wants
!> depends on how its assembly is written rather than on the model:
!>
!>   ConstitutiveStress_i  -- stress from the strain measure in Point. Every model
!>     has one. This is all the postprocessing needs, and it is also what an
!>     assembly that forms its tangent by Gateaux derivatives needs, applying the
!>     law once per (integration point, basis function, component).
!>   ConstitutiveTangent_i -- stress and the consistent tangent. Optional. What an
!>     assembly in B^T C B form needs, once per integration point.
!>
!> Measured before being written, because the plan for this interface assumed one
!> call per integration point and that is not what either Elmer elasticity solver
!> does: 7.6 ns a call for the stress shape and 21.4 ns with the tangent, against
!> an element that costs 20-40 microseconds to assemble. So the indirection is
!> affordable at either call frequency; what is not affordable is reshaping an
!> assembly to suit the interface, and hence the two entry points.
!>
!> Three conventions, each chosen against an alternative and each worth stating
!> because getting one wrong is silent:
!>
!>   TENSORS, NOT VOIGT. Strain and stress are 3x3. Elmer's own kernels already
!>     pass 3x3 -- Strain2Stress takes tensors and does the Voigt packing with the
!>     engineering shear internally, LocalStress hands back tensors -- so carrying
!>     Voigt here would introduce a second convention across every caller for the
!>     benefit of the one model, UMAT, whose adapter can convert. A factor of two
!>     on a shear term is the failure mode, and it is invisible in an isotropic
!>     test. Tangent is unavoidably Voigt: a fourth-order tensor in matrix form is.
!>   PROPS ARE PRE-EVALUATED BY THE DRIVER, as a plain array, not a ValueList
!>     pointer. A model that read its constants from the sif at every integration
!>     point would spend a keyword lookup per call -- string comparison, hundreds
!>     of nanoseconds at best -- which is one to two orders of magnitude more than
!>     the entire interface overhead measured above. The driver interpolates once
!>     per point and hands the numbers over. This is UMAT's own arrangement.
!>
!>     How much may Props carry? Arbitrarily much, and a matrix as readily as a
!>     scalar: the anisotropic model below takes a whole 6x6 elasticity matrix
!>     through it, flattened in Fortran's own column-major order so that the
!>     driver packs with RESHAPE and the model indexes back with the same
!>     arithmetic. That settles the question against a second, matrix-shaped
!>     channel into the interface, and it draws the line where it belongs:
!>
!>       MATERIAL data, however much of it, goes in Props;
!>       SOLUTION data goes in MaterialPoint_t.
!>
!>     The distinction is what separates the two remaining inline branches. The
!>     anisotropic elasticity matrix is material data that merely happens to vary
!>     from point to point, so Props takes it. A mixed formulation's pressure is a
!>     component of the solution, so it belongs with the deformation state in
!>     MaterialPoint_t and no amount of Props would be the right home for it.
!>   THE MODEL DECLARES ITS STRESS MEASURE rather than the driver assuming.
!>     Elmer already has one kernel read two ways: Strain2Stress is fed
!>     infinitesimal strain and read as Cauchy by StressSolve, and Green-Lagrange
!>     strain and read as second Piola-Kirchhoff by ElasticSolve. Nothing in its
!>     signature says which, and that ambiguity is what StressMeasure removes.
!------------------------------------------------------------------------------

MODULE Constitutive

  USE Types
  USE Messages
  USE LinearAlgebra

  IMPLICIT NONE

  !> Which measure the returned stress is in. Declared by the model; converted by
  !> the driver, which is the only party that knows the deformation gradient it
  !> would push forward with.
  INTEGER, PARAMETER :: STRESS_CAUCHY = 1
  INTEGER, PARAMETER :: STRESS_PK2    = 2

  !> Which strain measure the driver has put in Point % Strain.
  INTEGER, PARAMETER :: KINEMATICS_SMALL_STRAIN     = 1
  INTEGER, PARAMETER :: KINEMATICS_LARGE_DEFLECTION = 2
  INTEGER, PARAMETER :: KINEMATICS_HENCKY           = 3

!------------------------------------------------------------------------------
!> What the driver supplies about one integration point.
!------------------------------------------------------------------------------
  TYPE :: MaterialPoint_t
    REAL(KIND=dp) :: Strain(3,3) = 0.0_dp   !> the strain measure Kinematics names
    REAL(KIND=dp) :: dStrain(3,3) = 0.0_dp  !> its increment over the step
    REAL(KIND=dp) :: F(3,3) = 0.0_dp        !> deformation gradient now
    REAL(KIND=dp) :: F0(3,3) = 0.0_dp       !> and at the start of the step
    REAL(KIND=dp) :: DetF = 1.0_dp          !> its determinant, already to hand
    INTEGER :: Kinematics = KINEMATICS_SMALL_STRAIN
    REAL(KIND=dp) :: Temperature = 0.0_dp, dTemperature = 0.0_dp
    REAL(KIND=dp) :: Time = 0.0_dp, dt = 0.0_dp
    REAL(KIND=dp) :: Coord(3) = 0.0_dp, CharLength = 0.0_dp
    !> Dim is the coordinate system dimension and CDim the number of displacement
    !> components, which differ under axial symmetry. A model needs both: the
    !> trace runs over Dim while the identity it multiplies is CDim-diagonal.
    INTEGER :: Dim = 3, CDim = 3
    LOGICAL :: PlaneStress = .FALSE., AxiSymmetric = .FALSE.
    !> A pressure carried as an unknown of the system rather than derived from the
    !> deformation, which is what Elmer's "Mixed Formulation" does for a nearly
    !> incompressible material. It sits here and not in Props because it is a
    !> component of the SOLUTION -- the one piece of a constitutive response that a
    !> model cannot compute for itself and the driver must hand over. A model that
    !> also has a compressible form checks PressureSupplied and computes its own
    !> volumetric term when it is false; the constitutive relation stays in the
    !> model either way, which is the whole point of asking the driver only for
    !> what the driver alone knows.
    REAL(KIND=dp) :: Pressure = 0.0_dp
    LOGICAL :: PressureSupplied = .FALSE.
    INTEGER :: Element = 0, IP = 0
  END TYPE MaterialPoint_t

!------------------------------------------------------------------------------
!> History at one integration point. v is read and written, v0 is the converged
!> value at the start of the step and is read only.
!>
!> This is the one place two existing mechanisms collapse into one: Maxwell keeps
!> its previous stress in an Elmer -ip variable while UMAT keeps its own statev
!> array, for the same purpose. A model that wants neither leaves n at zero.
!------------------------------------------------------------------------------
  TYPE :: MaterialState_t
    REAL(KIND=dp), POINTER :: v(:) => NULL()
    REAL(KIND=dp), POINTER :: v0(:) => NULL()
    INTEGER :: n = 0
  END TYPE MaterialState_t

!------------------------------------------------------------------------------
!> What a model returns. Tangent is filled only by ConstitutiveTangent_i.
!>
!> The type carries Tangent on both paths deliberately. Splitting it to spare the
!> stress path thirty-six unused doubles looks like an optimisation and is not:
!> passing the same type through both entry points measured 7.6 ns against 21.4,
!> so the cost is filling the tangent, not carrying it -- these are passed by
!> reference.
!------------------------------------------------------------------------------
  TYPE :: MaterialResponse_t
    REAL(KIND=dp) :: Stress(3,3) = 0.0_dp
    REAL(KIND=dp) :: Tangent(6,6) = 0.0_dp
    INTEGER :: StressMeasure = STRESS_CAUCHY
    REAL(KIND=dp) :: EnergyElastic = 0.0_dp
    REAL(KIND=dp) :: EnergyPlastic = 0.0_dp
    REAL(KIND=dp) :: EnergyViscous = 0.0_dp
    REAL(KIND=dp) :: SuggestedDtFactor = 1.0_dp   !> UMAT's pnewdt
    LOGICAL :: Failed = .FALSE.
  END TYPE MaterialResponse_t

  ABSTRACT INTERFACE
    SUBROUTINE ConstitutiveStress_i( Point, Props, State, Response )
      IMPORT :: dp, MaterialPoint_t, MaterialState_t, MaterialResponse_t
      TYPE(MaterialPoint_t), INTENT(IN) :: Point
      REAL(KIND=dp), INTENT(IN) :: Props(:)
      TYPE(MaterialState_t), INTENT(INOUT) :: State
      TYPE(MaterialResponse_t), INTENT(OUT) :: Response
    END SUBROUTINE ConstitutiveStress_i

    SUBROUTINE ConstitutiveTangent_i( Point, Props, State, Response )
      IMPORT :: dp, MaterialPoint_t, MaterialState_t, MaterialResponse_t
      TYPE(MaterialPoint_t), INTENT(IN) :: Point
      REAL(KIND=dp), INTENT(IN) :: Props(:)
      TYPE(MaterialState_t), INTENT(INOUT) :: State
      TYPE(MaterialResponse_t), INTENT(OUT) :: Response
    END SUBROUTINE ConstitutiveTangent_i

    !> Many strains at ONE point, in one call. Optional, and purely a cost
    !> measure: a model that omits it loses nothing but speed, because
    !> ConstitutiveStresses below falls back to looping the single-strain entry.
    !>
    !> It exists because of what a Gateaux assembly actually asks for. Forming
    !> the tangent by directional derivatives applies the law once per
    !> (integration point, test function, component) -- 208 times per trilinear
    !> hexahedron where a B^T C B form asks 8 times -- and all 208 share one
    !> Point, one Props and one State. Only the strain differs. Measured on that
    !> assembly, routing it through the single-strain entry point added 25% to
    !> it, with the constitutive arithmetic eliminated as the cause: gutting the
    !> isotropic law to a single assignment left the overhead unchanged. Batching
    !> gives 21 of those 25 points back.
    !>
    !> POINT % STRAIN IS IGNORED HERE, and deliberately not quietly used as the
    !> first entry: the strains are the array, all n of them, and a model that
    !> reads Point % Strain in this entry point is reading a stale value from
    !> whatever the driver last did.
    !>
    !> Stresses must be written IN FULL, every one of the nine components of
    !> every slice. The single-strain entry gets its zeroing free, from
    !> MaterialResponse_t being INTENT(OUT) with default initialisers; here there
    !> is no such type, and a model that writes only the components its dimension
    !> reaches leaves the rest carrying whatever the driver's stack held. That is
    !> not hypothetical -- it is the exact defect (unzeroed out-of-plane row,
    !> masked by a multiplication by zero, surfacing only as an intermittent NaN)
    !> that this file's own history records.
    SUBROUTINE ConstitutiveStressBatch_i( Point, Props, State, n, Strains, Stresses )
      IMPORT :: dp, MaterialPoint_t, MaterialState_t
      TYPE(MaterialPoint_t), INTENT(IN) :: Point
      REAL(KIND=dp), INTENT(IN) :: Props(:)
      TYPE(MaterialState_t), INTENT(INOUT) :: State
      INTEGER, INTENT(IN) :: n
      REAL(KIND=dp), INTENT(IN) :: Strains(3,3,n)
      REAL(KIND=dp), INTENT(OUT) :: Stresses(3,3,n)
    END SUBROUTINE ConstitutiveStressBatch_i
  END INTERFACE

!------------------------------------------------------------------------------
!> What a model registers. A driver reads this to find out what it may do, rather
!> than the model's author deciding by writing another assembly routine.
!------------------------------------------------------------------------------
  TYPE :: MaterialModel_t
    PROCEDURE(ConstitutiveStress_i), POINTER, NOPASS :: Stress => NULL()
    PROCEDURE(ConstitutiveTangent_i), POINTER, NOPASS :: Tangent => NULL()
    !> The batched form of Stress, or null. A driver does not test this: it calls
    !> ConstitutiveStresses, which dispatches. See ConstitutiveStressBatch_i.
    PROCEDURE(ConstitutiveStressBatch_i), POINTER, NOPASS :: StressBatch => NULL()
    !> True when the response is linear in the strain measure, so that applying
    !> Stress to a strain INCREMENT gives the directional derivative of the
    !> stress. That is the licence for the cheap Gateaux path, and it is a
    !> property of the model: it holds for sigma = C:eps and fails for
    !> neo-Hookean, Hencky and UMAT.
    !>
    !> This flag is why ElasticSolve has three assembly routines rather than one.
    !> LocalMatrix takes the Gateaux shortcut; NeoHookeanLocalMatrix and
    !> LocalMatrixWithUMAT exist because the shortcut is not a derivative for
    !> them. They are one idea and two escapes from it, and a driver that reads
    !> this flag needs neither escape written out again.
    LOGICAL :: StrainLinear = .FALSE.
    INTEGER :: nState = 0             !> history slots wanted per integration point
    INTEGER :: nProps = 0             !> how many pre-evaluated constants expected
    CHARACTER(LEN=MAX_NAME_LEN) :: Name = ' '
  END TYPE MaterialModel_t

  !> Props layout for IsotropicLinearStress.
  INTEGER, PARAMETER :: ISOLIN_LAME1 = 1, ISOLIN_LAME2 = 2, ISOLIN_NPROPS = 2

  !> Props layout for AnisotropicLinearStress: the 6x6 elasticity matrix in
  !> Elmer's Voigt order (11,22,33,12,23,13) with engineering shear, flattened
  !> column major, which is what RESHAPE(C,[36]) produces.
  INTEGER, PARAMETER :: ANISOLIN_C = 1, ANISOLIN_NPROPS = 36

  !> Props layout for NeoHookeanStress. The same two constants the isotropic
  !> linear law takes, read as the compressible neo-Hookean parameters: LAME2 is
  !> the shear modulus and LAME1 the one governing the volumetric term, which is
  !> used only when the driver does not supply a pressure.
  INTEGER, PARAMETER :: NEOHOOKE_LAME1 = 1, NEOHOOKE_LAME2 = 2, NEOHOOKE_NPROPS = 2

CONTAINS

!------------------------------------------------------------------------------
!> The identity that the plane and axisymmetric reductions actually want, shared
!> by every model here because getting it wrong is how a plane stress case
!> silently becomes a plane strain one.
!>
!> It is not the plain identity: it is CDim-diagonal, and carries the out-of-plane
!> entry only away from plane stress. With that entry zeroed a volumetric term
!> makes no contribution to sigma_33, which is what leaves sigma_33 zero as plane
!> stress requires. Under axial symmetry the out-of-plane direction is the hoop and
!> does carry stress, so the entry is present whatever the plane stress flag says.
!------------------------------------------------------------------------------
  PURE FUNCTION ReducedIdentity( Point ) RESULT( Identity )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp) :: Identity(3,3)
!------------------------------------------------------------------------------
    INTEGER :: i
!------------------------------------------------------------------------------
    Identity = 0.0_dp
    DO i=1,Point % CDim
      Identity(i,i) = 1.0_dp
    END DO
    IF ( Point % AxiSymmetric .OR. .NOT. Point % PlaneStress ) Identity(3,3) = 1.0_dp
!------------------------------------------------------------------------------
  END FUNCTION ReducedIdentity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Evaluate a model at n strains sharing one point, through its batched entry
!> point where it has one and by looping the single-strain one where it has not.
!>
!> This is what a driver calls. It is a plain procedure and not a pointer on the
!> model, so that the fallback is written once here rather than as an "if the
!> model has a batch" at every call site -- a model gains the batch entry point
!> without its callers changing, and loses it the same way.
!------------------------------------------------------------------------------
  SUBROUTINE ConstitutiveStresses( Model, Point, Props, State, n, Strains, Stresses )
!------------------------------------------------------------------------------
    TYPE(MaterialModel_t), INTENT(IN) :: Model
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), INTENT(IN) :: Strains(3,3,n)
    REAL(KIND=dp), INTENT(OUT) :: Stresses(3,3,n)
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t) :: P
    TYPE(MaterialResponse_t) :: Response
    INTEGER :: k
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Model % StressBatch ) ) THEN
      CALL Model % StressBatch( Point, Props, State, n, Strains, Stresses )
      RETURN
    END IF

    ! The fallback. The point is copied once rather than per strain, and it is
    ! copied rather than aliased because Strain is the one field that varies and
    ! the caller's Point is INTENT(IN).
    P = Point
    DO k=1,n
      P % Strain = Strains(:,:,k)
      CALL Model % Stress( P, Props, State, Response )
      Stresses(:,:,k) = Response % Stress
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ConstitutiveStresses
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Isotropic linear elasticity: sigma = 2 mu eps + lambda tr(eps) I.
!>
!> Read as Cauchy from an infinitesimal strain and as second Piola-Kirchhoff from
!> a Green-Lagrange one -- the same arithmetic either way, which is why the
!> measure is taken from the kinematics rather than fixed here.
!>
!> The identity is not the plain one. It is CDim-diagonal, and carries the
!> out-of-plane entry only away from plane stress, which is how the plane stress
!> reduction enters: with it zeroed, the volumetric term makes no contribution to
!> sigma_33, and sigma_33 is then zero as plane stress requires. Under axial
!> symmetry the out-of-plane direction is the hoop and does carry stress, so the
!> entry is present whatever the plane stress flag says.
!------------------------------------------------------------------------------
  SUBROUTINE IsotropicLinearStress( Point, Props, State, Response )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    TYPE(MaterialResponse_t), INTENT(OUT) :: Response
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Identity(3,3), tr
    INTEGER :: i
!------------------------------------------------------------------------------
    Identity = ReducedIdentity( Point )

    tr = 0.0_dp
    DO i=1,Point % Dim
      tr = tr + Point % Strain(i,i)
    END DO

    Response % Stress = 2.0_dp * Props(ISOLIN_LAME2) * Point % Strain + &
        Props(ISOLIN_LAME1) * tr * Identity

    Response % StressMeasure = MERGE( STRESS_PK2, STRESS_CAUCHY, &
        Point % Kinematics /= KINEMATICS_SMALL_STRAIN )
!------------------------------------------------------------------------------
  END SUBROUTINE IsotropicLinearStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The batched form of the above. Identity and the two constants are hoisted out
!> of the loop; the arithmetic per strain is otherwise the same expression in the
!> same association, which is what makes the two forms agree to the bit rather
!> than merely to a tolerance.
!>
!> TwoMu is the one thing worth checking rather than assuming: the single-strain
!> form evaluates 2 * mu * eps left to right, so hoisting 2 * mu is exact -- and
!> exact for the ordinary reason that a factor of two is a power of two, not for
!> any general licence to rearrange.
!------------------------------------------------------------------------------
  SUBROUTINE IsotropicLinearStressBatch( Point, Props, State, n, Strains, Stresses )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), INTENT(IN) :: Strains(3,3,n)
    REAL(KIND=dp), INTENT(OUT) :: Stresses(3,3,n)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Identity(3,3), Lame1, TwoMu, tr
    INTEGER :: i, k, d
!------------------------------------------------------------------------------
    Identity = ReducedIdentity( Point )
    Lame1 = Props(ISOLIN_LAME1)
    TwoMu = 2.0_dp * Props(ISOLIN_LAME2)
    d = Point % Dim

    DO k=1,n
      tr = 0.0_dp
      DO i=1,d
        tr = tr + Strains(i,i,k)
      END DO
      Stresses(:,:,k) = TwoMu * Strains(:,:,k) + Lame1 * tr * Identity
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE IsotropicLinearStressBatch
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The model record for the above, for a driver to select and interrogate.
!------------------------------------------------------------------------------
  FUNCTION IsotropicLinearModel() RESULT( Model )
!------------------------------------------------------------------------------
    TYPE(MaterialModel_t) :: Model
!------------------------------------------------------------------------------
    Model % Stress => IsotropicLinearStress
    Model % StressBatch => IsotropicLinearStressBatch
    Model % StrainLinear = .TRUE.
    Model % nState = 0
    Model % nProps = ISOLIN_NPROPS
    Model % Name = 'isotropic linear'
!------------------------------------------------------------------------------
  END FUNCTION IsotropicLinearModel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Anisotropic linear elasticity: sigma = C : eps, with C given in full as the
!> pre-evaluated Props array.
!>
!> TWO PACKINGS, and which one applies is decided by Dim rather than chosen here.
!> In three dimensions C is the full 6x6 in Elmer's Voigt order. In two it is the
!> reduced three-component plane packing (11,22,12), and the caller is required to
!> have condensed it -- CondensePlaneElasticityMatrix in StressLocal is what does
!> that, moving the shear modulus C(4,4) into slot 3 and clearing the out-of-plane
!> couplings, plus a static condensation of the out-of-plane row under plane stress.
!>
!> That requirement cannot be checked from here and getting it wrong is silent, so
!> it is worth being explicit about the failure mode: handed a RAW 6x6 with Dim 2,
!> this routine reads C(3,3) -- the 33 modulus -- as the shear modulus, and the
!> 11-33 couplings as normal-to-shear ones. On the material in the anisotropic
!> tests that is 1346 in place of 385, wrong by three and a half times, and wrong in
!> the stiffness rather than only in the output.
!>
!> Axial symmetry DOES arrive here, and takes the three-dimensional packing: Dim is
!> the dimension of the state of stress, which is three there, and the components
!> are ordered (r, z, phi) with the hoop at index 3 -- so C is read in the Voigt
!> order (rr, zz, hoop, rz, z-hoop, r-hoop) it was written in.
!>
!> That was not always true and the history is worth keeping, because the failure it
!> describes is silent. ElasticSolve's assemblies used to order the axisymmetric
!> components (r, phi, z), hoop at index 2, while this file, its own postprocessor,
!> StressSolve and the UMAT interface all used (r, z, phi). An isotropic law cannot
!> see the difference -- the trace and the identity are blind to which axis is which
!> -- but for an anisotropic C it is exactly C(2,2) against C(3,3), so the case was
!> refused outright rather than answered wrongly. Aligning the conventions was the
!> whole of the fix; the permutation adapter that had been anticipated for the
!> elasticity matrix turned out not to be needed at all.
!------------------------------------------------------------------------------
  SUBROUTINE AnisotropicLinearStress( Point, Props, State, Response )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    TYPE(MaterialResponse_t), INTENT(OUT) :: Response
!------------------------------------------------------------------------------
    !> Elmer's Voigt order as index pairs, the same tables Strain2Stress builds --
    !> the full six in three dimensions, the leading three in the reduced plane
    !> packing, where slot 3 is the shear and there is no out-of-plane row.
    INTEGER, PARAMETER :: I1(6) = [ 1,2,3,1,2,1 ], I2(6) = [ 1,2,3,2,3,3 ]
    INTEGER, PARAMETER :: I1P(3) = [ 1,2,1 ], I2P(3) = [ 1,2,2 ]
    REAL(KIND=dp) :: S(6), csum
    INTEGER :: i, j, p, q, n
!------------------------------------------------------------------------------
    ! The off-diagonal entries are engineering shear, doubled here because C is
    ! indexed for it. A factor of two lost on these is invisible in any isotropic
    ! test, which is why they are written out rather than looped.
    SELECT CASE ( Point % Dim )
    CASE ( 2 )
      n = 3
      S(1) = Point % Strain(1,1)
      S(2) = Point % Strain(2,2)
      S(3) = 2.0_dp * Point % Strain(1,2)
    CASE ( 3 )
      n = 6
      S(1) = Point % Strain(1,1)
      S(2) = Point % Strain(2,2)
      S(3) = Point % Strain(3,3)
      S(4) = 2.0_dp * Point % Strain(1,2)
      S(5) = 2.0_dp * Point % Strain(2,3)
      S(6) = 2.0_dp * Point % Strain(1,3)
    CASE DEFAULT
      CALL Fatal( 'AnisotropicLinearStress', &
          'Material anisotropy implemented for two and three dimensions only' )
    END SELECT

    DO i=1,n
      IF ( n == 3 ) THEN
        p = I1P(i)
        q = I2P(i)
      ELSE
        p = I1(i)
        q = I2(i)
      END IF
      csum = 0.0_dp
      DO j=1,n
        csum = csum + Props(ANISOLIN_C - 1 + 6*(j-1) + i) * S(j)
      END DO
      Response % Stress(p,q) = csum
      Response % Stress(q,p) = csum
    END DO

    Response % StressMeasure = MERGE( STRESS_PK2, STRESS_CAUCHY, &
        Point % Kinematics /= KINEMATICS_SMALL_STRAIN )
!------------------------------------------------------------------------------
  END SUBROUTINE AnisotropicLinearStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The batched form of the above. The packing decision, the index tables and the
!> Fatal are hoisted out of the loop; the contraction itself is the same sum in
!> the same order, so the two forms agree to the bit.
!>
!> THE ZEROING IS NOT INCIDENTAL. The single-strain form gets it free, from
!> MaterialResponse_t being INTENT(OUT) with a default initialiser, and in two
!> dimensions it needs it: the loop writes four of the nine components and the
!> out-of-plane row and column are never touched. Here nothing zeroes them for
!> us. Left out, this routine would hand the assembly the driver's stack, in the
!> one entry a plane assembly multiplies by a dBasisdx that is zero -- silent
!> until the garbage happens to carry a NaN pattern, which is precisely the
!> defect valgrind found in this assembly once already.
!------------------------------------------------------------------------------
  SUBROUTINE AnisotropicLinearStressBatch( Point, Props, State, n, Strains, Stresses )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), INTENT(IN) :: Strains(3,3,n)
    REAL(KIND=dp), INTENT(OUT) :: Stresses(3,3,n)
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: I1(6) = [ 1,2,3,1,2,1 ], I2(6) = [ 1,2,3,2,3,3 ]
    INTEGER, PARAMETER :: I1P(3) = [ 1,2,1 ], I2P(3) = [ 1,2,2 ]
    REAL(KIND=dp) :: S(6), csum
    INTEGER :: i, j, k, p, q, m, P1(6), P2(6)
!------------------------------------------------------------------------------
    SELECT CASE ( Point % Dim )
    CASE ( 2 )
      m = 3
      P1(1:3) = I1P
      P2(1:3) = I2P
    CASE ( 3 )
      m = 6
      P1 = I1
      P2 = I2
    CASE DEFAULT
      CALL Fatal( 'AnisotropicLinearStressBatch', &
          'Material anisotropy implemented for two and three dimensions only' )
    END SELECT

    Stresses = 0.0_dp

    DO k=1,n
      ! Engineering shear on the off-diagonals, doubled because C is indexed for
      ! it. A factor of two lost here is invisible in any isotropic test.
      IF ( m == 3 ) THEN
        S(1) = Strains(1,1,k)
        S(2) = Strains(2,2,k)
        S(3) = 2.0_dp * Strains(1,2,k)
      ELSE
        S(1) = Strains(1,1,k)
        S(2) = Strains(2,2,k)
        S(3) = Strains(3,3,k)
        S(4) = 2.0_dp * Strains(1,2,k)
        S(5) = 2.0_dp * Strains(2,3,k)
        S(6) = 2.0_dp * Strains(1,3,k)
      END IF

      DO i=1,m
        p = P1(i)
        q = P2(i)
        csum = 0.0_dp
        DO j=1,m
          csum = csum + Props(ANISOLIN_C - 1 + 6*(j-1) + i) * S(j)
        END DO
        Stresses(p,q,k) = csum
        Stresses(q,p,k) = csum
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AnisotropicLinearStressBatch
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The model record for the above.
!------------------------------------------------------------------------------
  FUNCTION AnisotropicLinearModel() RESULT( Model )
!------------------------------------------------------------------------------
    TYPE(MaterialModel_t) :: Model
!------------------------------------------------------------------------------
    Model % Stress => AnisotropicLinearStress
    Model % StressBatch => AnisotropicLinearStressBatch
    Model % StrainLinear = .TRUE.
    Model % nState = 0
    Model % nProps = ANISOLIN_NPROPS
    Model % Name = 'anisotropic linear'
!------------------------------------------------------------------------------
  END FUNCTION AnisotropicLinearModel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compressible neo-Hookean elasticity, returning the second Piola-Kirchhoff
!> stress of the current iterate:
!>
!>   S = p C^-1 + mu ( I - C^-1 ),     C = F^T F
!>
!> This is the first model here that is NOT linear in its strain measure, so
!> StrainLinear is false and the Gateaux shortcut is not a derivative for it --
!> which is exactly why the calling solver keeps a separate assembly routine for
!> it. Nothing in this file depends on that; the flag is how a driver finds out.
!>
!> The pressure arrives one of two ways, and the split is the reason this model was
!> the one the interface could not absorb until MaterialPoint_t grew a field:
!>
!>   SUPPLIED -- Elmer's "Mixed Formulation" carries pressure as an unknown of the
!>     system, so it is solution data and only the driver can produce it. It comes
!>     in through Point % Pressure.
!>   COMPUTED -- otherwise the pressure is the material's own volumetric response
!>     to the dilatation, lambda/2 (J^2 - 1), which is a constitutive relation and
!>     so belongs here rather than in the driver that used to hold it.
!>
!> Written with (J-1)(J+1) rather than (J^2-1) deliberately: that is the form the
!> driver had, and the two are not the same in floating point.
!------------------------------------------------------------------------------
  SUBROUTINE NeoHookeanStress( Point, Props, State, Response )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    TYPE(MaterialResponse_t), INTENT(OUT) :: Response
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: InvC(3,3), Identity(3,3), Pres
!------------------------------------------------------------------------------
    IF ( Point % PressureSupplied ) THEN
      Pres = Point % Pressure
    ELSE
      Pres = Props(NEOHOOKE_LAME1)/2.0_dp * (Point % DetF - 1.0_dp) * &
          (Point % DetF + 1.0_dp)
    END IF

    Identity = ReducedIdentity( Point )

    ! The inverse of the right Cauchy-Green tensor. Inverted over Dim, the
    ! dimension of the state of stress, and not over the three rows it is stored
    ! in -- the out-of-plane row is not part of the plane system.
    InvC = MATMUL( TRANSPOSE(Point % F), Point % F )
    CALL InvertMatrix( InvC, Point % Dim )

    Response % Stress = Pres * InvC + Props(NEOHOOKE_LAME2) * (Identity - InvC)
    Response % StressMeasure = STRESS_PK2
!------------------------------------------------------------------------------
  END SUBROUTINE NeoHookeanStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The model record for the above.
!------------------------------------------------------------------------------
  FUNCTION NeoHookeanModel() RESULT( Model )
!------------------------------------------------------------------------------
    TYPE(MaterialModel_t) :: Model
!------------------------------------------------------------------------------
    Model % Stress => NeoHookeanStress
    !> Not linear in the strain measure: see the note on the type's own field.
    Model % StrainLinear = .FALSE.
    Model % nState = 0
    Model % nProps = NEOHOOKE_NPROPS
    Model % Name = 'neo-Hookean'
!------------------------------------------------------------------------------
  END FUNCTION NeoHookeanModel
!------------------------------------------------------------------------------

END MODULE Constitutive
!------------------------------------------------------------------------------
