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

!------------------------------------------------------------------------------
 MODULE MGDynMaterialUtils
!------------------------------------------------------------------------------
 USE DefUtils
 IMPLICIT NONE

 INTERFACE GetPermittivity
   MODULE PROCEDURE GetPermittivityR, GetPermittivityC
 END INTERFACE

 CONTAINS
!------------------------------------------------------------------------------
  FUNCTION GetElectricConductivityTensor(Element, n, Part, &
                   CoilBody,CoilType) RESULT (Tcoef)  
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), SAVE, POINTER :: Cwrk(:,:,:) => NULL()
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, i, j
    TYPE(Valuelist_t), POINTER :: Material
    REAL(KIND=dp) :: Tcoef(3,3,n)
    CHARACTER(LEN=*):: CoilType
    CHARACTER(LEN=2) :: Part
    LOGICAL :: Found
    LOGICAL :: CoilBody
!$OMP THREADPRIVATE(Cwrk)

    Tcoef=0._dp
    Material => GetMaterial( Element )
    IF ( ASSOCIATED(Material) ) THEN
      IF (Part=='re') THEN 
        CALL ListGetRealArray( Material, &
             'Electric Conductivity', Cwrk, n, Element % NodeIndexes, Found )
      ELSE
        CALL ListGetRealArray( Material, &
             'Electric Conductivity im', Cwrk, n, Element % NodeIndexes, Found )
      END IF 
      IF (Found) THEN
        IF (SIZE(Cwrk,1) == 1 .AND. SIZE(Cwrk,2) == 1) THEN
          DO i=1,3
            Tcoef(i,i,1:n) = Cwrk(1,1,1:n)
          END DO
        ELSE
          IF ( SIZE(Cwrk,1) == 1 ) THEN
            DO i=1,MIN(3, SIZE(Cwrk,2))
              Tcoef( i,i,1:n ) = Cwrk( 1,i,1:n )
            END DO
          ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
            DO i=1,MIN(3,SIZE(Cwrk,1))
              Tcoef(i,i,1:n) = Cwrk(i,1,1:n)
            END DO
          ELSE
            DO i=1,MIN(3,SIZE(Cwrk,1))
              DO j=1,MIN(3,SIZE(Cwrk,2))
                Tcoef( i,j,1:n ) = Cwrk(i,j,1:n)
              END DO
            END DO
          END IF
        END IF
      END IF
    END IF

    IF (CoilBody) THEN 
      SELECT CASE (CoilType)
      CASE ('stranded')
        !Tcoef(1,1,1:n) = 0._dp
        !Tcoef(2,2,1:n) = 0._dp
      CASE ('foil winding')
        Tcoef(1,1,1:n) = 0._dp
      END SELECT
    END IF
 
!------------------------------------------------------------------------------
  END FUNCTION GetElectricConductivityTensor
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 
  FUNCTION GetCMPLXElectricConductivityTensor(Element, n, CoilBody, CoilType) &
                  RESULT (TCoef) 
!------------------------------------------------------------------------------ 
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: TCoef(3,3,n)
    REAL(KIND=dp) :: TCoefRe(3,3,n), TCoefIm(3,3,n)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, i, j
    LOGICAL :: CoilBody
    CHARACTER(LEN=*) :: CoilType
    
    TCoef=0._dp
    TCoefRe=0._dp
    TCoefIm=0._dp
    TCoefRe = GetElectricConductivityTensor(Element,n,'re',CoilBody,CoilType)
    TCoefIm = GetElectricConductivityTensor(Element,n,'im',CoilBody,CoilType)
    DO i=1,3
       DO j=1,3
          Tcoef( i,j,1:n ) = CMPLX( TcoefRe( i,j,1:n ), TCoefIm( i,j,1:n ), KIND=dp)
       END DO
    END DO

!------------------------------------------------------------------------------ 
  END FUNCTION GetCMPLXElectricConductivityTensor
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------
  FUNCTION GetPermeabilityTensor(Element, n, Part) &
                  RESULT (mu)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), SAVE, POINTER :: Cwrk(:,:,:) => NULL()
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, i, j
    TYPE(Valuelist_t), POINTER :: Material
    REAL(KIND=dp) :: mu(3,3,n)
    CHARACTER(LEN=2) :: Part
    LOGICAL :: Found
!$OMP THREADPRIVATE(Cwrk)

    mu=0._dp
    Material => GetMaterial( Element )
    IF ( ASSOCIATED(Material) ) THEN
      IF (Part=='re') THEN 
        CALL ListGetRealArray( Material, &
             'Relative Permeability', Cwrk, n, Element % NodeIndexes, Found )
      ELSE
        CALL ListGetRealArray( Material, &
             'Relative Permeability im', Cwrk, n, Element % NodeIndexes, Found )
      END IF 
      IF (Found) THEN
        IF (SIZE(Cwrk,1) == 1 .AND. SIZE(Cwrk,2) == 1) THEN
          DO i=1,3
            mu(i,i,1:n) = Cwrk(1,1,1:n)
          END DO
        ELSE
          IF ( SIZE(Cwrk,1) == 1 ) THEN
            DO i=1,MIN(3, SIZE(Cwrk,2))
              mu( i,i,1:n ) = Cwrk( 1,i,1:n )
            END DO
          ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
            DO i=1,MIN(3,SIZE(Cwrk,1))
              mu(i,i,1:n) = Cwrk(i,1,1:n)
            END DO
          ELSE
            DO i=1,MIN(3,SIZE(Cwrk,1))
              DO j=1,MIN(3,SIZE(Cwrk,2))
                mu( i,j,1:n ) = Cwrk(i,j,1:n)
              END DO
            END DO
          END IF
        END IF
      END IF
    END IF
!------------------------------------------------------------------------------
  END FUNCTION GetPermeabilityTensor
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------
  FUNCTION GetTensor(Element, n, tsize, varname, Part, Found) &
                  RESULT (T)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), POINTER :: Cwrk(:,:,:) => NULL()
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, i, j, slen, tsize 
    TYPE(Valuelist_t), POINTER :: Material
    REAL(KIND=dp) :: T(tsize,tsize,n)
    CHARACTER(LEN=2) :: Part
    CHARACTER(LEN=*) :: varname
    LOGICAL, OPTIONAL :: Found
!$OMP THREADPRIVATE(Cwrk)

    IF (.NOT. ASSOCIATED(Element)) CALL Fatal ('GetTensor', 'Element not associated')
    T=0._dp
    Material => GetMaterial( Element )
    IF ( ASSOCIATED(Material) ) THEN
      slen = LEN_TRIM(varname)
      IF (Part=='re') THEN 
        CALL ListGetRealArray( Material, &
          varname(1:slen), Cwrk, n, Element % NodeIndexes, Found )
      ELSE
        CALL ListGetRealArray( Material, &
          varname(1:slen)//' im', Cwrk, n, Element % NodeIndexes, Found )
      END IF 
      IF (Found) THEN
        IF (SIZE(Cwrk,1) == 1 .AND. SIZE(Cwrk,2) == 1) THEN
          DO i=1,tsize
            T(i,i,1:n) = Cwrk(1,1,1:n)
          END DO
        ELSE
          IF ( SIZE(Cwrk,1) == 1 ) THEN
            DO i=1,MIN(tsize, SIZE(Cwrk,2))
              T( i,i,1:n ) = Cwrk( 1,i,1:n )
            END DO
          ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
            DO i=1,MIN(tsize,SIZE(Cwrk,1))
              T(i,i,1:n) = Cwrk(i,1,1:n)
            END DO
          ELSE
            DO i=1,MIN(tsize,SIZE(Cwrk,1))
              DO j=1,MIN(tsize,SIZE(Cwrk,2))
                T( i,j,1:n ) = Cwrk(i,j,1:n)
              END DO
            END DO
          END IF
        END IF
      END IF
    END IF
!------------------------------------------------------------------------------
  END FUNCTION GetTensor
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------
  FUNCTION GetCMPLXTensor(Element, n, tsize, varname, Found) &
                  RESULT (T)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, i, j, slen, tsize 
    COMPLEX(KIND=dp) :: T(tsize,tsize,n)
    REAL(KIND=dp) :: TRe(tsize,tsize,n), TIm(tsize,tsize,n)
    CHARACTER(LEN=2) :: Part
    CHARACTER(LEN=*) :: varname
    LOGICAL, OPTIONAL :: Found
    LOGICAL :: FoundRe, FoundIm

    T=0._dp
    TRe=0._dp
    TIm=0._dp
    TRe = GetTensor(Element,n,tsize,varname,'re',FoundRe)
    TIm = GetTensor(Element,n,tsize,varname,'im',FoundIm)
    Found = FoundRe .OR. FoundIm
    DO i=1,tsize
       DO j=1,tsize
          T( i,j,1:n ) = CMPLX( REAL(TRe( i,j,1:n )), TIm( i,j,1:n ), KIND=dp)
       END DO
    END DO
!------------------------------------------------------------------------------
  END FUNCTION GetCMPLXTensor
!------------------------------------------------------------------------------ 


!-------------------------------------------------------------------
  FUNCTION Get2x2MatrixInverse(M) &
   RESULT (Minv)
!-------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: M(2,2), Minv(2,2)
    REAL(KIND=dp) :: det, a, b, c, d 

    a = M(1,1); b = M(1,2); c=M(2,1); d=M(2,2)
    Minv=0._dp
  
    IF ( ABS(a) <= TINY(a) .AND. ABS(b) <= TINY(b) .AND. &
         ABS(c) <= TINY(c) .AND. ABS(d) <= TINY(d)         ) RETURN
    det = a*d-b*c
    IF (ABS(det) <= TINY(det)) CALL Fatal('Get2x2MatrixInverse', 'Determinant is zero! This should not happen...') 
    
    Minv(1,1) =  1/det * d
    Minv(1,2) = -1/det * b 
    Minv(2,1) = -1/det * c 
    Minv(2,2) =  1/det * a 
    
!-------------------------------------------------------------------
  END FUNCTION Get2x2MatrixInverse
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  FUNCTION Get2x2TensorInverse(T, n) &
    RESULT (Tinv)
!-------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: T(2,2,n), Tinv(2,2,n)
    INTEGER :: i, n

    DO i = 1, n
      Tinv(:,:,i) = Get2x2MatrixInverse(T(:,:,i))
    END DO

!-------------------------------------------------------------------
  END FUNCTION Get2x2TensorInverse
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  FUNCTION Get2x2CMPLXMatrixInverse(M) &
   RESULT (Minv)
!-------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: M(2,2), Minv(2,2)
    COMPLEX(KIND=dp) :: det, a, b, c, d 
    REAL(KIND=dp) :: r

    a = M(1,1); b = M(1,2); c=M(2,1); d=M(2,2)
    Minv=0._dp
  
    IF ( ABS(a) <= TINY(r) .AND. ABS(b) <= TINY(r) .AND. &
         ABS(c) <= TINY(r) .AND. ABS(d) <= TINY(r)         ) RETURN
    det = a*d-b*c
    IF (ABS(det) <= TINY(r)) CALL Fatal('Get2x2MatrixInverse', 'Determinant is zero! This should not happen...') 
    
    Minv(1,1) =  1/det * d
    Minv(1,2) = -1/det * b 
    Minv(2,1) = -1/det * c 
    Minv(2,2) =  1/det * a 
    
!-------------------------------------------------------------------
  END FUNCTION Get2x2CMPLXMatrixInverse
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  FUNCTION Get2x2CMPLXTensorInverse(T, n) &
    RESULT (Tinv)
!-------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: T(2,2,n), Tinv(2,2,n)
    INTEGER :: i, n

    DO i = 1, n
      Tinv(:,:,i) = Get2x2CMPLXMatrixInverse(T(:,:,i))
    END DO

!-------------------------------------------------------------------
  END FUNCTION Get2x2CMPLXTensorInverse
!-------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   USE CircuitUtils
   IMPLICIT NONE
   TYPE(Mesh_t), POINTER, SAVE :: Mesh
   TYPE(Element_t), POINTER :: Element
   TYPE(Valuelist_t), POINTER :: CompParams
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
   TYPE(Variable_t), POINTER, SAVE :: RotMvar !, alphavecvar 
   REAL(KIND=dp), POINTER, SAVE :: ConstArray(:,:)
   REAL(KIND=dp) :: Origin(3), alpha_ref(3), beta_ref(3)
   REAL(KIND=dp) :: x(3), r(3), xref(3)
   REAL(KIND=dp) :: C, S, t
   LOGICAL, SAVE :: visited = .FALSE.
   TYPE(Nodes_t), SAVE :: Nodes

   LOGICAL :: GotIt

   IF(.NOT. visited) THEN
     visited = .TRUE.
     Mesh => GetMesh()
     RotMvar => VariableGet( Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF

!     alphavecvar => VariableGet( Mesh % Variables, 'Alpha Vector E')
!     IF(.NOT. ASSOCIATED(alphavecvar)) THEN
!       CALL Fatal('RotMSolver()','Alpha Vector E variable not found')
!     END IF
   END IF

   RotM = 0._dp

   CompParams => GetComponentParams(Element)
   CALL GetConstRealArray( CompParams, ConstArray, 'Rotation Matrix Origin', GotIt )
   IF (GotIt) THEN
     IF (SIZE(Origin) /= 3) CALL Fatal('GetElementRotM', 'Rotation Matrix Origin needs three components!')
     Origin = ConstArray(1:3,1)
     CALL GetConstRealArray( CompParams, ConstArray, 'Rotation Matrix Alpha Reference', GotIt )
     IF (.NOT. GotIt) CALL Fatal('GetElementRotM', 'Rotation Matrix Origin set but Rotation Matrix Alpha Reference not found!')
     IF (SIZE(alpha_ref) /= 3) CALL Fatal('GetElementRotM', 'Rotation Matrix Alpha Reference needs three components!')
     alpha_ref = ConstArray(1:3,1)

     CALL GetConstRealArray( CompParams, ConstArray, 'Rotation Matrix Beta Reference', GotIt )
     IF (.NOT. GotIt) CALL Fatal('GetElementRotM', 'Rotation Matrix Origin set but Rotation Matrix Beta Reference not found!')
     IF (SIZE(beta_ref) /= 3) CALL Fatal('GetElementRotM', 'Rotation Matrix Beta Reference needs three components!')
     beta_ref = ConstArray(1:3,1)

     CALL GetElementNodes( Nodes )
     DO j = 1, n
       x(1) = Nodes % x(j)
       x(2) = Nodes % y(j)
       x(3) = Nodes % z(j)

       r = beta_ref/SQRT(SUM(beta_ref**2.)) ! Normalize the rotation axis

       xref = x - Origin ! take reference according to the origin of rotation
       xref = xref - SUM(xref * r) * r ! project xref to the rotation plane (beta_ref is the normal)
       C = SUM(xref * alpha_ref)/SQRT(SUM(xref**2.))/SQRT(SUM(alpha_ref**2.)) ! cosine of the angle of rotation
       S = SQRT(1-C**2) ! sine of the angle of rotation
       t = 1-C

       RotM(1,1,j) = t*r(1)**2+C
       RotM(1,2,j) = t*r(1)*r(2)-S*r(3)
       RotM(1,3,j) = t*r(1)*r(3)+S*r(2)

       RotM(2,1,j) = t*r(1)*r(2)+S*r(3)
       RotM(2,2,j) = t*r(2)**2+C
       RotM(2,3,j) = t*r(2)*r(3)-S*r(1)

       RotM(3,1,j) = t*r(1)*r(3)-S*r(2)
       RotM(3,2,j) = t*r(2)*r(3)+S*r(1)
       RotM(3,3,j) = t*r(3)**2+C
     END DO
     ! This is debug stuff (remove later if necessary)
!     DO j = 1, n
!       DO k=1,RotMvar % DOFs
!         RotMvar % Values(RotMvar % DOFs*(&
!           RotMvar % Perm(Element % DGIndexes(j))-1)+k) = RotM(ind1(k),ind2(k),j) 
!       END DO
!
!       IF (ASSOCIATED(alphavecvar)) THEN
!         x=(/1,1,1/)
!         x=MATMUL(RotM(:,:,j),x)
!         DO k=1,alphavecvar % DOFs
!           alphavecvar % Values( alphavecvar % DOFs*(alphavecvar % Perm( &
!                 Element % DGIndexes(j))-1)+k) = x(k)
!         END DO
!       END IF
!     END DO
   ELSE
     DO j = 1, n
       DO k=1,RotMvar % DOFs
         RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
               RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
       END DO
     END DO
   END IF

!       IF (SUM(Origin) == 0) THEN
!         x(1) = 0
!         x(2) = 1
!         x(3) = 0
!
!         r = beta_ref/SQRT(SUM(beta_ref**2.)) ! Normalize the rotation axis
!
!         xref = x - Origin ! take reference according to the origin of rotation
!         xref = xref - xref * beta_ref ! project xref to the rotation plane (beta_ref is the normal)
!         C = SUM(xref * alpha_ref)/SQRT(SUM(xref**2.))/SQRT(SUM(alpha_ref**2.)) ! cosine of the angle of rotation
!         S = SQRT(1-C**2) ! sine of the angle of rotation
!         t = 1+C
!
!         RotM(1,1,1) = t*r(1)**2+C
!         RotM(1,2,1) = t*r(1)*r(2)-S*r(3)
!         RotM(1,3,1) = t*r(1)*r(3)+S*r(2)
!
!         RotM(2,1,1) = t*r(1)*r(2)+S*r(3)
!         RotM(2,2,1) = t*r(2)**2+C
!         RotM(2,3,1) = t*r(2)*r(3)-S*r(1)
!
!         RotM(3,1,1) = t*r(1)*r(3)-S*r(2)
!         RotM(3,2,1) = t*r(2)*r(3)+S*r(1)
!         RotM(3,3,1) = t*r(3)**2+C
!
!         x = MATMUL(RotM(:,:,1),x)
!         print *, "RotM", RotM
!         print *, "x", x
!      END IF
!   END IF

       
!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetPermittivityR(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE., Warned = .FALSE.
    REAL(KIND=dp) :: Pvacuum
    SAVE FirstTime, Warned, Pvacuum
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Pvacuum = GetConstReal( CurrentModel % Constants, &
              'Permittivity of Vacuum', Found )
      IF (.NOT. Found) Pvacuum = 8.854187817d-12
      FirstTime = .FALSE.
    END IF

    Acoef(1:n) = GetReal( Material, 'Relative Permittivity', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Pvacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permittivity', Found )
    END IF

    IF( .NOT. Found ) THEN
      IF(.NOT. Warned ) THEN
        CALL Warn('GetPermittivity','Permittivity not defined in material, defaulting to that of vacuum')
        Warned = .TRUE.
      END IF
      Acoef(1:n) = Pvacuum
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetPermittivityR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetPermittivityC(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    COMPLEX(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, Found_im, FirstTime = .TRUE., Warned = .FALSE.
    REAL(KIND=dp) :: Pvacuum
    SAVE FirstTime, Warned, Pvacuum
    REAL(KIND=dp), PARAMETER :: im  = (0._dp, 1._dp)
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Pvacuum = GetConstReal( CurrentModel % Constants, &
              'Permittivity of Vacuum', Found )
      IF (.NOT. Found) Pvacuum = 8.854187817d-12
      FirstTime = .FALSE.
    END IF

    Acoef(1:n) = GetReal( Material, 'Relative Permittivity', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Pvacuum * Acoef(1:n)
      Acoef(1:n) = Acoef(1:n) + im * Pvacuum * & 
              GetReal(Material,'Relative Permittivity  im', Found_im )
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permittivity', Found )
      Acoef(1:n) = Acoef(1:n) + im * &
                 GetReal(Material,'Permittivity  im', Found_im )
    END IF

    IF( .NOT. Found ) THEN
      IF(.NOT. Warned ) THEN
        CALL Warn('GetPermittivity','Permittivity not defined in material, defaulting to that of vacuum')
        Warned = .TRUE.
      END IF
      Acoef(1:n) = Pvacuum
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetPermittivityC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This gets keywords for a loss model that is needed by FourierLossSolver and
!> MagnetoDynamicsCalcFields. There is an old and new format. 
!------------------------------------------------------------------------------
  SUBROUTINE GetLossExponents(vList,FreqPower,FieldPower,Ncomp,OldKeywords)
!------------------------------------------------------------------------------
    TYPE(ValueList_t),POINTER :: vList
    REAL(KIND=dp) :: FreqPower(:), FieldPower(:)
    INTEGER :: Ncomp
    LOGICAL, OPTIONAL :: OldKeywords
    
    REAL(KIND=dp), POINTER :: WrkArray(:,:)
    LOGICAL :: Found
    INTEGER :: icomp
    CHARACTER(*), PARAMETER :: Caller = 'GetLossExponents'

    IF( PRESENT(OldKeywords) ) THEN
      IF( OldKeywords ) THEN      
        CALL Info('GetLossExponents','Using old keyword format',Level=20)
        FreqPower(1) = GetCReal( vList,'Harmonic Loss Linear Frequency Exponent',Found )
        IF( .NOT. Found ) FreqPower(1) = 1.0_dp
        
        FreqPower(2) = GetCReal( vList,'Harmonic Loss Quadratic Frequency Exponent',Found )
        IF( .NOT. Found ) FreqPower(2) = 2.0_dp
        
        FieldPower(1) = GetCReal( vList,'Harmonic Loss Linear Exponent',Found ) 
        IF( .NOT. Found ) FieldPower(1) = 2.0_dp
        FieldPower(1) = FieldPower(1) / 2.0_dp
        
        FieldPower(2) = GetCReal( vList,'Harmonic Loss Quadratic Exponent',Found ) 
        IF( .NOT. Found ) FieldPower(2) = 2.0_dp
        FieldPower(2) = FieldPower(2) / 2.0_dp    
        RETURN
      ELSE
        CALL Info('GetLossExponents','Using new keyword format',Level=20)    
      END IF
    END IF
      
    WrkArray => ListGetConstRealArray( vList,'Harmonic Loss Frequency Exponent',Found )    
    IF( Found ) THEN 
      IF( SIZE( WrkArray,1 ) < Ncomp ) THEN
        CALL Fatal(Caller,'> Harmonic Loss Frequency Exponent < too small')
      END IF
      FreqPower(1:Ncomp) = WrkArray(1:Ncomp,1)
    ELSE       
      DO icomp = 1, Ncomp
        FreqPower(icomp) = GetCReal( vList,'Harmonic Loss Frequency Exponent '//I2S(icomp) )
      END DO
    END IF

    WrkArray => ListGetConstRealArray( vList,'Harmonic Loss Field Exponent',Found )
    IF( Found ) THEN
      IF( SIZE( WrkArray,1 ) < Ncomp ) THEN
        CALL Fatal(Caller,'> Harmonic Loss Field Exponent < too small')
      END IF
      FieldPower(1:Ncomp) = WrkArray(1:Ncomp,1)        
    ELSE
      DO icomp = 1, Ncomp
        FieldPower(icomp) = GetCReal( vList,'Harmonic Loss Field Exponent '//I2S(icomp) )
      END DO
    END IF    
    
  END SUBROUTINE GetLossExponents

  
!------------------------------------------------------------------------------
 END MODULE MGDynMaterialUtils
!------------------------------------------------------------------------------
