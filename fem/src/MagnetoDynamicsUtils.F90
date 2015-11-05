!------------------------------------------------------------------------------
 MODULE MGDynMaterialUtils
!------------------------------------------------------------------------------
 USE DefUtils

 CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE GetElectricConductivityTensor(Tcoef, Element, n, Part, &
                   CoilBody,CoilType) 
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    REAL(KIND=dp), POINTER :: Cwrk(:,:,:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, i, j
    TYPE(Valuelist_t), POINTER :: Material
    REAL(KIND=dp) :: Tcoef(3,3,n)
    CHARACTER(LEN=MAX_NAME_LEN):: CoilType
    CHARACTER(LEN=2) :: Part
    LOGICAL :: Found
    LOGICAL :: CoilBody

    Tcoef=0._dp
    NULLIFY( Cwrk )
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
        Tcoef(1,1,1:n) = 0._dp
        Tcoef(2,2,1:n) = 0._dp
      CASE ('foil winding')
        Tcoef(1,1,1:n) = 0._dp
      END SELECT
    END IF
 
!------------------------------------------------------------------------------
  END SUBROUTINE GetElectricConductivityTensor
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 
  SUBROUTINE GetCMPLXElectricConductivityTensor(TCoef, Element, n, CoilBody, CoilType)
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
    CALL GetElectricConductivityTensor(TCoefRe,Element,n,'re',CoilBody,CoilType)
    CALL GetElectricConductivityTensor(TCoefIm,Element,n,'im',CoilBody,CoilType)
    DO i=1,3
       DO j=1,3
          Tcoef( i,j,1:n ) = CMPLX( REAL(TcoefRe( i,j,1:n )), TCoefIm( i,j,1:n ), KIND=dp)
       END DO
    END DO

!------------------------------------------------------------------------------ 
  END SUBROUTINE GetCMPLXElectricConductivityTensor
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------
 END MODULE MGDynMaterialUtils
!------------------------------------------------------------------------------
