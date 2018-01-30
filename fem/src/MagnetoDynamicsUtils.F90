!------------------------------------------------------------------------------
 MODULE MGDynMaterialUtils
!------------------------------------------------------------------------------
 USE DefUtils
 IMPLICIT NONE

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
    CHARACTER(LEN=MAX_NAME_LEN):: CoilType
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
         IF ( SIZE(Cwrk,1) == 1 ) THEN
            DO i=1,3
               Tcoef( i,i,1:n ) = Cwrk( 1,1,1:n )
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
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    
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
         IF ( SIZE(Cwrk,1) == 1 ) THEN
            DO i=1,3
               mu( i,i,1:n ) = Cwrk( 1,1,1:n )
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
         IF ( SIZE(Cwrk,1) == 1 ) THEN
            DO i=1,tsize
               T( i,i,1:n ) = Cwrk( 1,1,1:n )
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
 END MODULE MGDynMaterialUtils
!------------------------------------------------------------------------------
