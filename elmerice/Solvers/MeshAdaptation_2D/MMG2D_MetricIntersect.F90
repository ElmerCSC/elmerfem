!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Author: F. Gillet-Chaulet (IGE)
! *  Email:  fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:    http://elmerice.elmerfem.org
! *
! *  Original Date: 13-07-2017, 
! *****************************************************************************
!> Compute the intersection of 2 metrics
!>  
!------------------------------------------------------------------------------
   SUBROUTINE MMG2D_MetricIntersect_Init( Model,Solver,dt,TransientSimulation )
   USE DefUtils
   IMPLICIT NONE
   !------------------------------------------------------------------------------
   TYPE(Solver_t), TARGET :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
   !--------------------------------------------------------------------------
   CHARACTER(LEN=MAX_NAME_LEN) :: Name,TensorName,ExportName
   TYPE(ValueList_t), POINTER :: SolverParams
   TYPE(Variable_t), POINTER :: TensorVariable
   LOGICAL :: GotIt

   SolverParams => Solver % Values 

   Name = ListGetString( SolverParams, 'Equation',GotIt)
   IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
        CALL ListAddString( SolverParams,'Variable',&
           '-nooutput '//TRIM(Name)//'_var')
   ENDIF

   IF(.NOT. ListCheckPresent(SolverParams,'Optimize Bandwidth')) &
        CALL ListAddLogical(SolverParams,'Optimize Bandwidth',.FALSE.)
   
   TensorName = ListGetString( SolverParams, 'Metric Variable Name',  UnFoundFatal=.TRUE. )
   TensorVariable => VariableGet( Solver % Mesh % Variables, TensorName, ThisOnly=.TRUE. )
   IF (.NOT.ASSOCIATED(TensorVariable)) THEN
      ExportName=NextFreeKeyword('Exported Variable',SolverParams)
      CALL ListAddString( SolverParams,TRIM(ExportName),'-nooutput '//TRIM(TensorName))
      CALL ListAddInteger( SolverParams,TRIM(ExportName)//' DOFs',3)
   ENDIF

   END SUBROUTINE MMG2D_MetricIntersect_Init
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE MMG2D_MetricIntersect( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL ::  TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='MMG_MetricAniso'

    TYPE(ValueList_t),POINTER :: SolverParams

    TYPE(Variable_t), POINTER :: TensorVariable,Tensor1Variable,Tensor2Variable
    REAL(KIND=dp),DIMENSION(:),POINTER :: TensorValues,Tensor1Values,Tensor2Values
    INTEGER,DIMENSION(:),POINTER :: TensorPerm,Tensor1Perm,Tensor2Perm
    CHARACTER(LEN=MAX_NAME_LEN) :: TensorName,Tensor1Name,Tensor2Name

    REAL(KIND=dp),DIMENSION(2,2) :: Metric,M,M1,M2,N
    REAL(KIND=dp),DIMENSION(2,2) :: M1I,Nev,NevI
    REAL(KIND=dp),DIMENSION(2) :: Nei
    REAL(KIND=dp) :: lambda,mu

    INTEGER :: node,k,i,j
    INTEGER, DIMENSION(3),PARAMETER :: CI=(/1,2,1/),CJ=(/1,2,2/)

    LOGICAL :: UnFoundFatal=.TRUE.,GotIt

!------------------------------------------------------------------------------
     SolverParams => GetSolverParams(Solver)

    TensorName = ListGetString( SolverParams, 'Metric Variable Name',  UnFoundFatal=UnFoundFatal )    
    TensorVariable => VariableGet( Solver % Mesh % Variables, TensorName, UnFoundFatal=UnFoundFatal ) 
    TensorPerm => TensorVariable % Perm
    TensorValues => TensorVariable % Values
    IF (TensorVariable % DOFs /= 3) &
          CALL Fatal( SolverName, 'Bad dimension of Tensor Variable; should be 3' )

    Tensor1Name = ListGetString( SolverParams, 'Metric 1 Variable Name',  UnFoundFatal=UnFoundFatal )    
    Tensor1Variable => VariableGet( Solver % Mesh % Variables, Tensor1Name, UnFoundFatal=UnFoundFatal ) 
    Tensor1Perm => Tensor1Variable % Perm
    Tensor1Values => Tensor1Variable % Values
    IF (Tensor1Variable % DOFs /= 3) &
          CALL Fatal( SolverName, 'Bad dimension of Tensor Variable; should be 3' )

    Tensor2Name = ListGetString( SolverParams, 'Metric 2 Variable Name',  UnFoundFatal=UnFoundFatal )    
    Tensor2Variable => VariableGet( Solver % Mesh % Variables, Tensor2Name, UnFoundFatal=UnFoundFatal ) 
    Tensor2Perm => Tensor2Variable % Perm
    Tensor2Values => Tensor2Variable % Values
    IF (Tensor2Variable % DOFs /= 3) &
          CALL Fatal( SolverName, 'Bad dimension of Tensor Variable; should be 3' )


    DO node=1,Solver % Mesh % NumberOfNodes
       Do k=1,3
         M1(CI(k),CJ(k))=Tensor1Values(3*(Tensor1Perm(node)-1)+k)
         M2(CI(k),CJ(k))=Tensor2Values(3*(Tensor2Perm(node)-1)+k)
       End do
       k=3
       M1(CJ(k),CI(k))=M1(CI(k),CJ(k))
       M2(CJ(k),CI(k))=M2(CI(k),CJ(k))
       CALL Inv2D(M1,M1I)
       N=MATMUL(M1I(:,:),M2(:,:))
       CALL Eigen2D(N,Nei,Nev)
       M=0._dp
       Do k=1,2
         lambda=0._dp
         mu=0._dp
         Do i=1,2
          Do j=1,2 
            lambda=lambda+M1(i,j)*Nev(i,k)*Nev(j,k)
            mu=mu+M2(i,j)*Nev(i,k)*Nev(j,k)
          End do
         End do
         M(k,k)=max(lambda,mu)
       End do
       CALL Inv2D(Nev,NevI)
       Metric=0._dp
       Metric(:,:)=MATMUL(TRANSPOSE(NevI(:,:)),MATMUL(M(:,:),NevI(:,:)))
       Do k=1,3
          TensorValues(3*(TensorPerm(node)-1)+k)=Metric(CI(k),CJ(k))
       End Do
    END DO
      
CONTAINS
! Return eigenvalues (lambda) of 2x2 matrice a
! Return eigenvectors (columns of e) of 2x2 matrice a
      SUBROUTINE Eigen2D(a,lambda,e)
      implicit none
      REAL(KIND=dp),intent(in) :: a(2,2)
      REAL(KIND=dp),intent(out) :: lambda(2)
      REAL(KIND=dp),intent(out) :: e(2,2)

      REAL(KIND=dp) :: T,D,en
      REAL(KIND=dp),parameter :: Zero=1.0e-8  !Zero=100*AEPS
      INTEGER :: i

      T=a(1,1)+a(2,2)
      D=a(1,1)*a(2,2)-a(1,2)*a(2,1)

      lambda(1)=0.5_dp*T+sqrt(0.25_dp*T*T-D)
      lambda(2)=0.5_dp*T-sqrt(0.25_dp*T*T-D)


      e=0._dp
     ! only 1 double eigenvalue (ONLY append when interscting isotropic
     ! metric??) - TO CHECK
      IF (abs((lambda(1)-lambda(2))/T).LT.Zero) THEN
         e(1,1)=1._dp
         e(2,2)=1._dp
      ELSEIF (abs(a(2,1)/T).GT.Zero) THEN
       !first eigenvector
         e(1,1)=lambda(1)-a(2,2)
         e(2,1)=a(2,1)
        !2nd eigenvector
         e(1,2)=lambda(2)-a(2,2)
         e(2,2)=a(2,1)
      ELSE IF (abs(a(1,2)/T).GT.Zero) THEN
       !first eigenvector
         e(1,1)=a(1,2)
         e(2,1)=lambda(1)-a(1,1)
        !2nd eigenvector
         e(1,2)=a(1,2)
         e(2,2)=lambda(2)-a(1,1)
      ELSE
         lambda(1)=a(1,1)
         lambda(2)=a(2,2)
         e(1,1)=1._dp
         e(2,2)=1._dp
      END IF
      
      DO i=1,2
         en=SQRT(SUM(e(:,i)*e(:,i)))
         e(:,i)=e(:,i)/en
      END DO
      END SUBROUTINE Eigen2D
        
! Return inverse (e) of 2x2 matrice a
      SUBROUTINE Inv2D(a,e)
      implicit none
      REAL(KIND=dp),intent(in) :: a(2,2)
      REAL(KIND=dp),intent(out) :: e(2,2)

      REAL(KIND=dp) :: D
      REAL(KIND=dp),parameter :: Zero=TINY(1.0_dp)

      D=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      IF (abs(D).LE.Zero) THEN
         PRINT *,a(1,1),a(2,1),a(1,2),a(2,2),D
         CALL FATAL('Inv2D','2D Matrix has no inverse')
      END IF

      e(1,1)=a(2,2)
      e(2,2)=a(1,1)
      e(1,2)=-a(1,2)
      e(2,1)=-a(2,1)
      e(:,:)=e(:,:)/D

      END SUBROUTINE Inv2D
!------------------------------------------------------------------------------
      END SUBROUTINE MMG2D_MetricIntersect
!------------------------------------------------------------------------------
