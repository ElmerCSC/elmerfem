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
! *  Some common subroutines for the solvers of solid mechanics
! *
! *  Authors: Mika Malinen, Mikko Lyly
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Oct 8, 2020
! *
! *****************************************************************************/

MODULE SolidMechanicsUtils

  USE DefUtils
  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Integrate and assemble the local stiffness matrix corresponding to the 
!> one-dimensional Timoshenko beam equations. The local DOFs always 
!> correspond to the displacement components along the tangent direction and the
!> principal axes of the cross section. The transformation to global DOFs is done
!> within this subroutine. The stiffness matrix K corresponding to the global 
!> DOFs is thus obtained as K = R^T k R and the RHS vector F is obtained as 
!> F = R^T f.
!------------------------------------------------------------------------------
  SUBROUTINE BeamStiffnessMatrix(Element, n, nd, nb, TransientSimulation, &
      MassAssembly, HarmonicAssembly, LargeDeflection, LocalSol, RHSForce, &
      CombineWithShell, ApplyRotation, DrillingDOFs)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, nb
    LOGICAL, INTENT(IN) :: TransientSimulation
    LOGICAL, INTENT(IN) :: MassAssembly               ! To activate mass matrix integration
    LOGICAL, OPTIONAL, INTENT(IN) :: HarmonicAssembly ! To activate the global mass matrix updates
    LOGICAL, OPTIONAL, INTENT(IN) :: LargeDeflection  ! To activate nonlinear terms
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: LocalSol(:,:) ! The previous solution iterate
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: RHSForce(:)  ! Local RHS vector corresponding to external loads
    LOGICAL, OPTIONAL, INTENT(IN) :: CombineWithShell    ! Set .TRUE. if the caller is the shell solver 
    LOGICAL, OPTIONAL, INTENT(IN) :: ApplyRotation    ! Rotate DOFs in the context of shell analysis
    LOGICAL, OPTIONAL, INTENT(IN) :: DrillingDOFs     ! Assume drilling DOFs in the context of shell analysis
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes, LocalNodes
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Found, Stat
    LOGICAL :: NonlinAssembly, RotationNeeded
    LOGICAL :: ApplyOffset
    LOGICAL :: DampingBetaWarning = .FALSE.
    
    INTEGER :: DOFs
    INTEGER :: i, t, p, q
    INTEGER :: i0, p0, q0
    INTEGER :: MomentFreeAxis

    REAL(KIND=dp), POINTER :: ArrayPtr(:,:) => NULL()
    REAL(KIND=dp), POINTER :: StiffBlock(:,:), MassBlock(:,:), DampBlock(:,:)
    REAL(KIND=dp), DIMENSION(3), PARAMETER :: ZBasis = (/ 0.0d0, 0.0d0, 0.1d1 /)

    REAL(KIND=dp), TARGET :: Mass(6*nd,6*nd), Stiff(6*nd,6*nd), Damp(6*nd,6*nd)
    REAL(KIND=dp) :: Force(6*nd)
    REAL(KIND=dp) :: RBlock(3,3), R(6*nd,6*nd)
    REAL(KIND=dp) :: Basis(nd), dBasis(nd,3), DetJ, Weight
    REAL(KIND=dp) :: Youngs_Modulus(n), Shear_Modulus(n), Area(n), Density(n)
    REAL(KIND=dp) :: Form_Factor(n)
    REAL(KIND=dp) :: Torsional_Constant(n) 
    REAL(KIND=dp) :: Area_Moment_2(n), Area_Moment_3(n)
    REAL(KIND=dp) :: Offset_2, Offset_3
    REAL(KIND=dp) :: Mass_Inertia_Moment(n), Damping(n), RayleighBeta(n)
    REAL(KIND=dp) :: Load(3,n), f(3)
    REAL(KIND=dp) :: PrevSolVec(6*nd)
    REAL(KIND=dp) :: E, A, G, rho, DampCoef, FormFact
    REAL(KIND=dp) :: EA, GA, MOI, Mass_per_Length 
    REAL(KIND=dp) :: E_diag(3)

    REAL(KIND=dp) :: p1(3), p2(3), e1(3), e2(3), e3(3)
    REAL(KIND=dp) :: L, Norm

    SAVE Nodes, LocalNodes, DampingBetaWarning
!------------------------------------------------------------------------------
    IF (n > 2) CALL Fatal('BeamStiffnessMatrix', &
        'Only 2-node background meshes supported currently')

    DOFs = 6
!    dim = CoordinateSystemDimension()

    CALL GetElementNodes(Nodes)

    Mass  = 0.0_dp
    Stiff = 0.0_dp
    Damp = 0.0_dp
    Force = 0.0_dp

    IF (PRESENT(RHSForce)) RHSForce = 0.0d0
    IF (PRESENT(LargeDeflection)) THEN
      NonlinAssembly = LargeDeflection
    ELSE
      NonlinAssembly = .FALSE.
    END IF
    IF (NonlinAssembly) THEN
      IF (.NOT. PRESENT(LocalSol)) CALL Fatal('BeamStiffnessMatrix', &
          'Previous solution iterate needed')
      DO i=1,DOFs
        PrevSolVec(i:DOFs*(nd-nb):DOFs) = LocalSol(i,1:(nd-nb))
      END DO  
    END IF


    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) THEN
      !
      ! Force components refer to the basis of the global frame:
      !
      Load(1,1:n) = GetReal(BodyForce, 'Body Force 1', Found)
      Load(2,1:n) = GetReal(BodyForce, 'Body Force 2', Found)
      Load(3,1:n) = GetReal(BodyForce, 'Body Force 3', Found)
    ELSE
      Load = 0.0_dp
    END IF

    Material => GetMaterial()
    Youngs_Modulus(1:n) = GetReal(Material, 'Youngs Modulus', Found)
    IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Youngs Modulus needed')
    Shear_Modulus(1:n) = GetReal(Material, 'Shear Modulus', Found)
    IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Shear Modulus needed')
    Form_Factor(1:n) = GetReal(Material, 'Shear Correction Factor', Found)
    IF (.NOT. Found) Form_Factor(1:n) = 1.0_dp
    Area(1:n) = GetReal(Material, 'Cross Section Area', Found)
    IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Cross Section Area needed')
    Torsional_Constant(1:n) = GetReal(Material, 'Torsional Constant', Found)
    IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Torsional Constant needed')

    ! If we don't give the moment of area component-wise, it is assumed to be the same
    Area_Moment_2(1:n) = GetReal(Material, 'Second Moment of Area', Found)
    IF( Found ) THEN
      Area_Moment_3(1:n) = Area_Moment_2(1:n)
    ELSE
      Area_Moment_2(1:n) = GetReal(Material, 'Second Moment of Area 2', Found)
      IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Second Moment of Area 2 needed')
      Area_Moment_3(1:n) = GetReal(Material, 'Second Moment of Area 3', Found)
      IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Second Moment of Area 3 needed')
    END IF
      
    IF (MassAssembly) THEN
      Density(1:n) = GetReal(Material, 'Density', Found)
      IF (.NOT. Found) CALL Fatal('BeamStiffnessMatrix', 'Density needed')
      Damping(1:n) = GetReal(Material, 'Rayleigh Damping Alpha', Found)
      RayleighBeta = GetReal(Material, 'Rayleigh Damping Beta', Found)
      IF (Found .AND. .NOT.DampingBetaWarning) THEN
        CALL Warn('BeamStiffnessMatrix', 'Only mass-proportional damping, neglecting Rayleigh Damping Beta = ...')
        DampingBetaWarning = .TRUE.
      END IF
    END IF

    !
    ! Compute the tangent vector e1 to the beam axis:
    !
    p1(1) = Nodes % x(1)
    p1(2) = Nodes % y(1)
    p1(3) = Nodes % z(1)
    p2(1) = Nodes % x(2)
    p2(2) = Nodes % y(2)
    p2(3) = Nodes % z(2)
    e1 = p2 - p1
    L = SQRT(SUM(e1(:)**2))
    e1 = 1.0_dp/L * e1
    !
    ! Cross section parameters are given with respect to a local frame. 
    ! Determine its orientation:
    !
    ArrayPtr => ListGetConstRealArray(Material, 'Director', Found)
    IF (Found) THEN
      e3 = 0.0d0
      DO i=1,SIZE(ArrayPtr,1)
        e3(i) = ArrayPtr(i,1)
      END DO
      Norm = SQRT(SUM(e3(:)**2))
      e3 = 1.0_dp/Norm * e3
      IF (ABS(DOT_PRODUCT(e1,e3)) > 100.0_dp * AEPS) CALL Fatal('BeamStiffnessMatrix', &
          'Director should be orthogonal to the beam axis')
      e2 = CrossProduct(e3, e1)
    ELSE
      ArrayPtr => ListGetConstRealArray(Material, 'Principal Direction 2', Found)
      IF (Found) THEN
        e2 = 0.0d0
        DO i=1,SIZE(ArrayPtr,1)
          e2(i) = ArrayPtr(i,1)
        END DO
        Norm = SQRT(SUM(e2(:)**2))
        e2 = 1.0_dp/Norm * e2     
        e3 = CrossProduct(e1, e2)
      ELSE
        IF (ANY( ABS(Area_Moment_3(1:n) - Area_Moment_2(1:n)) > 5.0_dp * AEPS )) THEN
          CALL Fatal('BeamStiffnessMatrix', &
              'The default principal directions need the same moments of area')
        END IF
        !e2 = -ZBasis
        !e3 = CrossProduct(e1, e2)
        CALL TangentDirections( e1, e2, e3 ) 
      END IF
      IF (ABS(DOT_PRODUCT(e1,e2)) > 100.0_dp * AEPS) CALL Fatal('BeamStiffnessMatrix', &
          'Principal Direction 2 should be orthogonal to the beam axis')
    END IF

 
    !
    ! Allocate an additional variable so as to write nodes data with respect to
    ! the local frame.
    !
    IF (.NOT. ASSOCIATED(LocalNodes % x)) THEN
      ALLOCATE(LocalNodes % x(n), LocalNodes % y(n), LocalNodes % z(n) ) 
      LocalNodes % NumberOfNodes = n
      LocalNodes % y(:) = 0.0_dp
      LocalNodes % z(:) = 0.0_dp
    END IF
    LocalNodes % x(1) = 0.0d0
    LocalNodes % x(2) = L

    !-----------------------
    ! Numerical integration:
    !-----------------------
    IF (.NOT. IsActivePElement(Element) .AND. nd > n) THEN
      IP = GaussPoints(Element, 3)
    ELSE
      IP = GaussPoints(Element)
    END IF

    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo(Element, LocalNodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasis)

      ! Create a bubble if the element is the standard 2-node element:
      IF (.NOT. IsActivePElement(Element) .AND. nd > n) THEN
        Basis(n+1) = Basis(1) * Basis(2)
        dBasis(3,:) = dBasis(1,:) * Basis(2) + Basis(1) * dBasis(2,:)
      END IF

      !------------------------------------------
      ! The model data at the integration point:
      !------------------------------------------
      f(1) = SUM(Basis(1:n) * Load(1,1:n))
      f(2) = SUM(Basis(1:n) * Load(2,1:n))
      f(3) = SUM(Basis(1:n) * Load(3,1:n))      

      ! TO DO: Add option to give the applied moment load

      E = SUM(Basis(1:n) * Youngs_Modulus(1:n))
      G = SUM(Basis(1:n) * Shear_Modulus(1:n))
      FormFact = SUM(Basis(1:n) * Form_Factor(1:n))
      A = SUM(Basis(1:n) * Area(1:n))

      E_diag(1) = G * SUM(Basis(1:n) * Torsional_Constant(1:n))
      E_diag(2) = E * SUM(Basis(1:n) * Area_Moment_2(1:n))
      E_diag(3) = E * SUM(Basis(1:n) * Area_Moment_3(1:n)) 

      IF (MassAssembly) THEN
        rho = SUM(Basis(1:n) * Density(1:n))
        MOI = rho/E * (E_diag(2) + E_diag(3))
        Mass_per_Length = rho * A
        DampCoef = SUM(Damping(1:n) * Basis(1:n))
      END IF

      GA = FormFact*G*A
      EA = E*A

      ! TO DO: Add option to give shear correction factors

      Weight = IP % s(t) * DetJ

      DO p=1,nd
        p0 = (p-1)*DOFs
        DO q=1,nd
          q0 = (q-1)*DOFs
          StiffBlock => Stiff(p0+1:p0+DOFs,q0+1:q0+DOFs)
          MassBlock => Mass(p0+1:p0+DOFs,q0+1:q0+DOFs)
          DampBlock => Damp(p0+1:p0+DOFs,q0+1:q0+DOFs)
          !
          ! (Du',v'):
          !
          StiffBlock(1,1) = StiffBlock(1,1) + &
              EA * dBasis(q,1) * dBasis(p,1) * Weight
          StiffBlock(2,2) = StiffBlock(2,2) + &
              GA * dBasis(q,1) * dBasis(p,1) * Weight
          StiffBlock(3,3) = StiffBlock(3,3) + &
              GA * dBasis(q,1) * dBasis(p,1) * Weight
  
          IF (MassAssembly) THEN
            MassBlock(1,1) = MassBlock(1,1) + &
                Mass_per_Length * Basis(q) * Basis(p) * Weight
            MassBlock(2,2) = MassBlock(2,2) + &
                Mass_per_Length * Basis(q) * Basis(p) * Weight
            MassBlock(3,3) = MassBlock(3,3) + &
                Mass_per_Length * Basis(q) * Basis(p) * Weight
            
            DampBlock(1,1) = DampBlock(1,1) + &
                DampCoef * Mass_per_Length * Basis(q) * Basis(p) * Weight
            DampBlock(2,2) = DampBlock(2,2) + &
                DampCoef * Mass_per_Length * Basis(q) * Basis(p) * Weight
            DampBlock(3,3) = DampBlock(3,3) + &
                DampCoef * Mass_per_Length * Basis(q) * Basis(p) * Weight
          END IF

          IF (q > n) CYCLE
          !
          ! -(D theta x t,v'):
          !
          StiffBlock(2,6) = StiffBlock(2,6) - &
              GA * Basis(q) * dBasis(p,1) * Weight
          StiffBlock(3,5) = StiffBlock(3,5) + &
              GA * Basis(q) * dBasis(p,1) * Weight
        END DO
        
        Force(p0+1) = Force(p0+1) + Weight * DOT_PRODUCT(f,e1)* Basis(p)
        Force(p0+2) = Force(p0+2) + Weight * DOT_PRODUCT(f,e2)* Basis(p)
        Force(p0+3) = Force(p0+3) + Weight * DOT_PRODUCT(f,e3)* Basis(p)

        IF (p > n) CYCLE

        DO q=1,nd
          q0 = (q-1)*DOFs
          StiffBlock => Stiff(p0+1:p0+DOFs,q0+1:q0+DOFs)
          MassBlock => Mass(p0+1:p0+DOFs,q0+1:q0+DOFs)
          !
          ! -(D u',psi x t):
          !
          StiffBlock(5,3) = StiffBlock(5,3) + &
              GA * Basis(p) * dBasis(q,1) * Weight
          StiffBlock(6,2) = StiffBlock(6,2) - &
              GA * Basis(p) * dBasis(q,1) * Weight

          IF (q > n) CYCLE

          !
          ! (E theta',psi') + (D theta x t,psi x t):
          !
          StiffBlock(4,4) = StiffBlock(4,4) + &
              E_diag(1) * dBasis(q,1) * dBasis(p,1) * Weight
          StiffBlock(5,5) = StiffBlock(5,5) + &
              E_diag(2) * dBasis(q,1) * dBasis(p,1) * Weight + &
              GA * Basis(p) * Basis(q) * Weight
          StiffBlock(6,6) = StiffBlock(6,6) + &
              E_diag(3) * dBasis(q,1) * dBasis(p,1) * Weight + &
              GA * Basis(p) * Basis(q) * Weight

          IF (MassAssembly) THEN
            MassBlock(4,4) = MassBlock(4,4) + MOI * Basis(q) * Basis(p) * Weight
            MassBlock(5,5) = MassBlock(5,5) + rho/E * E_diag(2) * &
                Basis(q) * Basis(p) * Weight
            MassBlock(6,6) = MassBlock(6,6) + rho/E * E_diag(3) * &
                Basis(q) * Basis(p) * Weight
          END IF

        END DO
      END DO
    END DO

    CALL BeamCondensate(nd-nb, nb, DOFs, 3, Stiff, Force)
    
    IF (PRESENT(CombineWithShell)) THEN
      IF (CombineWithShell) THEN

        Offset_3 = GetConstReal(Material, 'Beam Axis Offset 3', ApplyOffset)
        Offset_2 = GetConstReal(Material, 'Beam Axis Offset 2', Found)
        ApplyOffset = Found .OR. ApplyOffset
        IF (ApplyOffset) THEN
          R = 0.0d0
          DO i=1,nd-nb
            i0 = (i-1)*DOFs
            R(i0+1,i0+1) = 1.0d0
            R(i0+1,i0+5) = Offset_3
            R(i0+1,i0+6) = -Offset_2
            R(i0+2,i0+2) = 1.0d0
            R(i0+3,i0+3) = 1.0d0
            R(i0+4,i0+4) = 1.0d0
            R(i0+5,i0+5) = 1.0d0
            R(i0+6,i0+6) = 1.0d0
          END DO
          DOFs = (nd-nb)*DOFs
          Stiff(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)), &
              MATMUL(Stiff(1:DOFs,1:DOFs),R(1:DOFs,1:DOFs)))
          Force(1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)),Force(1:DOFs))

          IF (MassAssembly) &
              Mass(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)), &
              MATMUL(Mass(1:DOFs,1:DOFs),R(1:DOFs,1:DOFs)))
          DOFs = 6
        END IF
        
        IF (PRESENT(ApplyRotation)) THEN
          RotationNeeded = ApplyRotation
        ELSE
          RotationNeeded = .TRUE.
        END IF

        IF (PRESENT(DrillingDOFs)) THEN
          IF (DrillingDOFs) RotationNeeded = .FALSE.
        END IF
        
        IF (RotationNeeded) THEN
          !
          ! Switch to rotation variables which conform with the rotated moments - M x d:
          !
          R = 0.0d0
          DO i=1,nd-nb
            i0 = (i-1)*DOFs
            R(i0+1,i0+1) = 1.0d0
            R(i0+2,i0+2) = 1.0d0
            R(i0+3,i0+3) = 1.0d0
            R(i0+4,i0+5) = 1.0d0
            R(i0+5,i0+4) = -1.0d0
            R(i0+6,i0+6) = 1.0d0
          END DO
          DOFs = (nd-nb)*DOFs
          Stiff(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)), &
              MATMUL(Stiff(1:DOFs,1:DOFs),R(1:DOFs,1:DOFs)))
          Force(1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)),Force(1:DOFs))

          IF (MassAssembly) &
              Mass(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)), &
              MATMUL(Mass(1:DOFs,1:DOFs),R(1:DOFs,1:DOFs)))
          DOFs = 6
        END IF

        !
        ! The moment around the director is not compatible with the shell model.
        ! Remove its contribution:
        !
        MomentFreeAxis = GetInteger(Material, 'Moment-free Axis', Found)
        IF (.NOT. Found) MomentFreeAxis = 3
          
        DO p=1,nd-nb
          i = (p-1)*DOFs + 3 + MomentFreeAxis
          Stiff(i,:) = 0.0d0
          Stiff(:,i) = 0.0d0
          !
          ! If the beam axis lies on the mid-surface, assembling a zero
          ! diagonal entry would be consistent. However, if the beam continues
          ! outside the mid-surface, a row of zeroes can cause troubles.
          ! As a remedy, we now give a minimal rigidity against a moment load.
          ! TO DO: Develop a more consistent strategy to handle this trouble
          !
          Stiff(i,i) = 1.0d1 * AEPS
          Force(i) = 0.0d0
          Mass(i,:) = 0.0d0
          Mass(:,i) = 0.0d0
        END DO
      END IF
    END IF

    !
    ! Build the transformation matrix in order to switch to the global DOFs
    !
    DOFs = 6
    R = 0.0d0
    RBlock(1,1:3) = e1(1:3)
    RBlock(2,1:3) = e2(1:3)
    RBlock(3,1:3) = e3(1:3)
    DO i=1,nd-nb
      i0 = (i-1)*DOFs
      R(i0+1:i0+3,i0+1:i0+3) =  RBlock(1:3,1:3)
      R(i0+4:i0+6,i0+4:i0+6) =  RBlock(1:3,1:3)
    END DO

    !-------------------------------------------------------
    ! Transform to the global DOFs:
    !-------------------------------------------------------
    DOFs = (nd-nb)*DOFs
    Stiff(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)), &
        MATMUL(Stiff(1:DOFs,1:DOFs),R(1:DOFs,1:DOFs)))
    Force(1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)),Force(1:DOFs))

    IF (PRESENT(RHSForce)) RHSForce(1:DOFs) = Force(1:DOFs)
    IF (NonlinAssembly) Force(1:DOFs) = Force(1:DOFs) - &
        MATMUL(Stiff(1:DOFs,1:DOFs), PrevSolVec(1:DOFs))

    IF (MassAssembly) THEN
      Mass(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(R(1:DOFs,1:DOFs)), &
          MATMUL(Mass(1:DOFs,1:DOFs),R(1:DOFs,1:DOFs)))
      IF (TransientSimulation) THEN
        CALL Default2ndOrderTime(Mass, Damp, Stiff, Force)
      ELSE IF (PRESENT(HarmonicAssembly)) THEN
        IF (HarmonicAssembly) CALL DefaultUpdateMass(Mass)
      END IF
    END IF

    CALL DefaultUpdateEquations(Stiff, Force)
!------------------------------------------------------------------------------
  END SUBROUTINE BeamStiffnessMatrix
!------------------------------------------------------------------------------

    
!------------------------------------------------------------------------------
  SUBROUTINE BeamCondensate(n, nb, dofs, dim, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n    ! Nodes after condensation
    INTEGER, INTENT(IN) :: nb   ! The number of bubble basis functions
    INTEGER, INTENT(IN) :: dofs ! DOFs per node
    INTEGER, INTENT(IN) :: dim  ! The first dim fields have bubbles
    REAL(KIND=dp), INTENT(INOUT) :: K(:,:)          ! The stiffness matrix
    REAL(KIND=dp), INTENT(INOUT) :: F(:)            ! The RHS vector
    REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: F1(:) ! Some other RHS vector
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Kbl(nb*dim,n*dofs), Kbb(nb*dim,nb*dim), Fb(nb*dim)
    REAL(KIND=dp) :: Klb(n*dofs,nb*dim)
    
    INTEGER :: i, m, p, Cdofs(dofs*n), Bdofs(dim*nb)
!------------------------------------------------------------------------------
    
    Cdofs(1:n*dofs) = (/ (i, i=1,n*dofs) /)

    m = 0
    DO p = 1,nb
      DO i = 1,dim
        m = m + 1
        Bdofs(m) = dofs*(n+p-1) + i
      END DO
    END DO

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Cdofs)
    Klb = K(Cdofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb*dim )

    F(1:dofs*n) = F(1:dofs*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    K(1:dofs*n,1:dofs*n) = &
        K(1:dofs*n,1:dofs*n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )

    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:dofs*n) = F1(1:dofs*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE BeamCondensate
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
!> Perform the operation
!>
!>    A = A + C' * B * C * s
!>
!> with
!>
!>    Size( A ) = n x n
!>    Size( B ) = m x m
!>    Size( C ) = m x n
!------------------------------------------------------------------------------
  SUBROUTINE StrainEnergyDensity(A, B, C, m, n, s)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: A(:,:)
    REAL(KIND=dp), INTENT(IN) :: B(:,:), C(:,:)
    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp), INTENT(IN) :: s
!------------------------------------------------------------------------------
    A(1:n,1:n) = A(1:n,1:n) + s * MATMUL(TRANSPOSE(C(1:m,1:n)),MATMUL(B(1:m,1:m),C(1:m,1:n))) 
!------------------------------------------------------------------------------
  END SUBROUTINE StrainEnergyDensity
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
  SUBROUTINE Jacobi3(Jmat, invJ, detJ, x, y)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(OUT) :: Jmat(:,:), invJ(:,:), detJ
    REAL(KIND=dp), INTENT(IN) :: x(:), y(:)
!------------------------------------------------------------------------------
    Jmat(1,1) = x(2)-x(1)
    Jmat(2,1) = x(3)-x(1)
    Jmat(1,2) = y(2)-y(1)
    Jmat(2,2) = y(3)-y(1)

    detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

    invJ(1,1) =  Jmat(2,2)/detJ
    invJ(2,2) =  Jmat(1,1)/detJ
    invJ(1,2) = -Jmat(1,2)/detJ
    invJ(2,1) = -Jmat(2,1)/detJ
!------------------------------------------------------------------------------
  END SUBROUTINE Jacobi3
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Jacobi4(Jmat, invJ, detJ, xi, eta, x, y)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(OUT) :: Jmat(:,:), invJ(:,:), detJ
    REAL(KIND=dp), INTENT(IN) :: xi, eta, x(:), y(:)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dNdxi(4), dNdeta(4)
    INTEGER :: i
!------------------------------------------------------------------------------
    dNdxi(1) = -(1-eta)/4.0d0
    dNdxi(2) =  (1-eta)/4.0d0
    dNdxi(3) =  (1+eta)/4.0d0
    dNdxi(4) = -(1+eta)/4.0d0
    dNdeta(1) = -(1-xi)/4.0d0
    dNdeta(2) = -(1+xi)/4.0d0
    dNdeta(3) =  (1+xi)/4.0d0
    dNdeta(4) =  (1-xi)/4.0d0
       
    Jmat = 0.0d0
    DO i=1,4
      Jmat(1,1) = Jmat(1,1) + dNdxi(i)*x(i)
      Jmat(1,2) = Jmat(1,2) + dNdxi(i)*y(i)
      Jmat(2,1) = Jmat(2,1) + dNdeta(i)*x(i)
      Jmat(2,2) = Jmat(2,2) + dNdeta(i)*y(i)
    END DO

    detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

    invJ(1,1) = Jmat(2,2)/detJ
    invJ(2,2) = Jmat(1,1)/detJ
    invJ(1,2) = -Jmat(1,2)/detJ
    invJ(2,1) = -Jmat(2,1)/detJ
!------------------------------------------------------------------------------
  END SUBROUTINE Jacobi4
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ShearCorrectionFactor(Kappa, Thickness, x, y, n, StabParam)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(OUT) :: Kappa
    REAL(KIND=dp), INTENT(IN) :: Thickness, x(:), y(:)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: StabParam
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: x21,x32,x43,x13,x14,y21,y32,y43,y13,y14, &
        l21,l32,l43,l13,l14,alpha,h
    REAL(KIND=dp) :: StabPar
!------------------------------------------------------------------------------
    IF (PRESENT(StabParam)) THEN
      StabPar = StabParam
    ELSE
      StabPar = 1.0d0
    END IF

    Kappa = 1.0d0
    SELECT CASE(n)
    CASE(3)
      alpha = 0.20d0 * StabPar
      x21 = x(2)-x(1)
      x32 = x(3)-x(2)
      x13 = x(1)-x(1)
      y21 = y(2)-y(1)
      y32 = y(3)-y(2)
      y13 = y(1)-y(1)
      l21 = SQRT(x21**2 + y21**2)
      l32 = SQRT(x32**2 + y32**2)
      l13 = SQRT(x13**2 + y13**2)
      h = MAX(l21,l32,l13)
      Kappa = (Thickness**2)/(Thickness**2 + alpha*(h**2))
    CASE(4)
      alpha = 0.10d0 * StabPar
      x21 = x(2)-x(1)
      x32 = x(3)-x(2)
      x43 = x(4)-x(3)
      x14 = x(1)-x(4)
      y21 = y(2)-y(1)
      y32 = y(3)-y(2)
      y43 = y(4)-y(3)
      y14 = y(1)-y(4)
      l21 = SQRT(x21**2 + y21**2)
      l32 = SQRT(x32**2 + y32**2)
      l43 = SQRT(x43**2 + y43**2)
      l14 = SQRT(x14**2 + y14**2)
      h = MAX(l21,l32,l43,l14)
      Kappa = (Thickness**2)/(Thickness**2 + alpha*(h**2))
    CASE DEFAULT
      CALL Fatal('ShearCorrectionFactor','Illegal number of nodes for Smitc elements: '//I2S(n))
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE ShearCorrectionFactor
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE IsotropicElasticity(Ematrix, Gmatrix, Poisson, Young, Thickness,&
      Basis, n)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Ematrix(:,:), Gmatrix(:,:), Basis(:)
    REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:)
    REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw
    INTEGER :: n
!------------------------------------------------------------------------------
    Euvw = SUM( Young(1:n)*Basis(1:n) )
    Puvw = SUM( Poisson(1:n)*Basis(1:n) )
    Tuvw = SUM( Thickness(1:n)*Basis(1:n) )
    Guvw = Euvw/(2.0d0*(1.0d0 + Puvw))

    Ematrix = 0.0d0
    Ematrix(1,1) = 1.0d0
    Ematrix(1,2) = Puvw
    Ematrix(2,1) = Puvw
    Ematrix(2,2) = 1.0d0
    Ematrix(3,3) = (1.0d0-Puvw)/2.0d0

    Ematrix = Ematrix* Euvw * (Tuvw**3) / (12.0d0*(1.0d0-Puvw**2))
    
    Gmatrix = 0.0d0
    Gmatrix(1,1) = Guvw*Tuvw
    Gmatrix(2,2) = Guvw*Tuvw
!------------------------------------------------------------------------------
  END SUBROUTINE IsotropicElasticity
!------------------------------------------------------------------------------

  
END MODULE SolidMechanicsUtils

