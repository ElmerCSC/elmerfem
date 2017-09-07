!-------------------------------------------------------------
! Gives the x(or y) component of tangential reaction force
! multiplied by a friction coefficient. 
!
! The subroutine is needed bacause when we use the model for 
! friction force where the nodal loads in normal direction is
! multiplied by a friction coefficient we end in problems for
! normal-tangential coordinate system. There we need to eliminate
! the tangential reaction force from the loads. 
!-------------------------------------------------------------
  FUNCTION TangentForce( Model, n, t ) RESULT(f)
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: t,f
!-------------------------------------------------------------
    INTEGER :: i,j,LoadDofs
    REAL(KIND=dp) :: Normal(3), Force(3),NormalForce,Mu
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Variable_t), POINTER :: LoadVar
    INTEGER, POINTER :: LoadPerm(:)
    REAL(KIND=dp), POINTER :: LoadValues(:)
    LOGICAL :: Visited = .FALSE.

    SAVE Visited, ElementNodes, LoadPerm, LoadDofs, LoadValues


    IF( .NOT. Visited ) THEN
      LoadVar => VariableGet( Model % Variables,'Displacement Loads')
      IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
        CALL Fatal('TangentForce','No > Displacement Loads < given!')
      END IF
      LoadDofs = LoadVar % Dofs
      IF( LoadDofs /= 2 ) THEN
        CALL Fatal('TangentForce','Implemented only for 2D!')
      END IF

      LoadValues => LoadVar % Values
      LoadPerm => LoadVar % Perm
    END IF

    Mu = 0.5_dp

    Element => Model % CurrentElement       
    CALL GetElementNodes( ElementNodes, Element )
    Normal = NormalVector( Element, ElementNodes, 0.0d0, 0.0d0 )

    Force = 0.0_dp
    DO i=1,LoadDofs
      j = LoadPerm(n)
      IF( j == 0 ) CYCLE
      Force(i) = LoadValues( LoadDofs * (j-1) + i )
    END DO

    ! This sign convention seems to be consistent for outer normals
    NormalForce = -SUM( Force * Normal ) 

    f = Mu * NormalForce

  END FUNCTION TangentForce
!-------------------------------------------------------------

