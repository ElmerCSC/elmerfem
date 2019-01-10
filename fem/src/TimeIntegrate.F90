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
! *  Original Date: 09 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!> Time integration schemes for first and second order partial differential
!> equation with in time. 
!------------------------------------------------------------------------------
MODULE TimeIntegrate

   USE Types

   IMPLICIT NONE

CONTAINS

!>
!------------------------------------------------------------------------------

   SUBROUTINE RungeKutta( N, dt, MassMatrix, StiffMatrix, &
                   Force, PrevSolution, CurrSolution )
!------------------------------------------------------------------------------
    INTEGER :: N
    REAL(KIND=dp) :: Force(:),PrevSolution(:),dt
    REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:),CurrSolution(:)

!------------------------------------------------------------------------------
    INTEGER :: i,j,NB1,NB2

    REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
    NB1 = SIZE( StiffMatrix,1 )
    NB2 = SIZE( StiffMatrix,2 )

    DO i=1,NB1
      s = 0.0_dp
      DO j=1,n
            s = s + (1/dt) * MassMatrix(i,j) * PrevSolution(j) - &
                 StiffMatrix(i,j) * CurrSolution(j)
      END DO

      DO j=1,NB2
         StiffMatrix(i,j) = (1/dt)*MassMatrix(i,j)
      END DO
      Force(i) = Force(i) + s
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE RungeKutta
!------------------------------------------------------------------------------

!> Apply newmark beta scheme to local matrix equation.
!> This is used also for Implicit Euler with beta=1.0 and
!> Crank-Nicolson with beta=0.5.
!------------------------------------------------------------------------------

   SUBROUTINE NewmarkBeta( N, dt, MassMatrix, StiffMatrix, &
                   Force, PrevSolution, Beta )
!------------------------------------------------------------------------------
    INTEGER :: N
    REAL(KIND=dp) :: Force(:),PrevSolution(:),dt
    REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:),Beta

!------------------------------------------------------------------------------
    INTEGER :: i,j,NB1,NB2

    REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
    NB1 = SIZE( StiffMatrix,1 )
    NB2 = SIZE( StiffMatrix,2 )

    DO i=1,NB1
      s = 0.0_dp
      DO j=1,n
            s = s + (1/dt) * MassMatrix(i,j) * PrevSolution(j) - &
                      (1-Beta) * StiffMatrix(i,j) * PrevSolution(j)
      END DO

      DO j=1,NB2
         StiffMatrix(i,j) = Beta * StiffMatrix(i,j) + (1/dt)*MassMatrix(i,j)
      END DO
      Force(i) = Force(i) + s
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NewmarkBeta
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------   
! Adams-Bashforth and Adams-Moulton time-integration schemes implemented
! by Gong Cheng, 2017. 
!   
!> Apply Adams-Bashforth method to local matrix equation.
!> This method stores prev_stiff*prevSol in Element % propertydata and 
!> the corresponding forcing terms.
!> PredCorrOrder=1, Explicit Euler
!> PredCorrOrder=2, 2nd order Adams-Bashforth
!------------------------------------------------------------------------------
   SUBROUTINE AdamsBashforth( N, dt, MassMatrix, StiffMatrix, &
                   Force, PrevSolution, zeta, PredCorrOrder)
!------------------------------------------------------------------------------
    INTEGER :: N    ! Size of the unknowns
    REAL(KIND=dp) :: Force(:),PrevSolution(:),dt, zeta
    REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:)
    TYPE(Element_t), POINTER :: Element
    TYPE(elementdata_t), POINTER :: tempRes
    LOGICAL :: GotIt
    INTEGER :: PredCorrOrder

!------------------------------------------------------------------------------
     INTEGER :: i,j,NB1,NB2

     REAL(KIND=dp) :: s_curr, m_curr, residual, curr_res, preForce
!------------------------------------------------------------------------------
     NB1 = SIZE( StiffMatrix,1 )
     NB2 = SIZE( StiffMatrix,2 ) 

     Element => CurrentModel % CurrentElement
     IF (.NOT. ASSOCIATED(Element % propertydata)) THEN
       ALLOCATE( Element % propertydata )
       ALLOCATE( Element % propertydata % values(NB1*2) )  
     END IF


     DO i=1,NB1
       s_curr = 0.0_dp
       m_curr = 0.0_dp
       DO j=1,n
         s_curr = s_curr + StiffMatrix(i,j) * PrevSolution(j) 
         m_curr = m_curr + (1/dt) * MassMatrix(i,j) * PrevSolution(j) 
       END DO

       curr_res = - s_curr

       IF ( PredCorrOrder == 1 ) THEN
         residual = curr_res
         preForce = Force(i)
       ELSE
         residual = Element % propertydata % values(i)
         preForce = Element % propertydata % values(i+NB1)        
       END IF
       Element % propertydata % values(i) = curr_res
       Element % propertydata % values(i+NB1) = Force(i)

      Force(i) = Force(i) - s_curr  + m_curr + 0.5_dp * zeta * (Force(i) - preForce) + &
                      0.5_dp * zeta * (curr_res - residual)

       DO j=1,NB2
         StiffMatrix(i,j) = (1/dt) * MassMatrix(i,j)
       END DO

     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE AdamsBashforth
!------------------------------------------------------------------------------


!> Apply Adams-Moulton(Corrector) method to local matrix equation.
!>
!> A two steps method with second order accuracy in time.
!> This method is only used in the correction phase of adaptive timestepping,
!> PrevSolution(:,2) -- the true solution from previous correction step H^{n-1}
!> PrevSolution(:,1) -- the corrector \tilde{H^n}
!>
!> This method can only be used in the corrector phase of Predictor-Corrector 
!> scheme and just after Adams-Bashforth method with the same order, otherwise
!> the residual at n-1 step will be incorrect.
!> PredCorrOrder=1, Implicit Euler
!> PredCorrOrder=2, 2nd order Adams-Moulton

!------------------------------------------------------------------------------
   SUBROUTINE AdamsMoulton( N, dt, MassMatrix, StiffMatrix, &
       Force, PrevSolution, PredCorrOrder)
!------------------------------------------------------------------------------
     INTEGER :: N
     REAL(KIND=dp) :: Force(:),PrevSolution(:,:),dt
     REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:)
     INTEGER :: PredCorrOrder
!------------------------------------------------------------------------------
     INTEGER :: i,j,NB1,NB2
     TYPE(Element_t), POINTER :: Element

     REAL(KIND=dp) :: s_curr, m_curr, residual, preForce
!------------------------------------------------------------------------------
     NB1 = SIZE( StiffMatrix,1 )
     NB2 = SIZE( StiffMatrix,2 ) 

     Element => CurrentModel % CurrentElement

     IF (.NOT. ASSOCIATED(Element % propertydata)) THEN
       CALL Fatal( 'AdamsMoulton', &
           'Adams-Moulton method must be executed after Adams-Bashforth method!')
     END IF

     DO i=1,NB1
       s_curr = 0.0_dp
       m_curr = 0.0_dp
       DO j=1,n
         s_curr = s_curr + StiffMatrix(i,j) * PrevSolution(j,1) 
         m_curr = m_curr + (1/dt) * MassMatrix(i,j) * PrevSolution(j,2) 
       END DO
         
       DO j=1,NB2
         StiffMatrix(i,j) =   (1/dt) * MassMatrix(i,j)
       END DO

       residual = Element % propertydata % values(i)
       preForce = Element % propertydata % values(i+NB1)        
       IF ( PredCorrOrder == 1 ) THEN
         Force(i) =  Force(i) + m_curr - s_curr  
       ELSE
         Force(i) =  0.5_dp * (Force(i) + preForce) + m_curr + 0.5_dp * (-s_curr + residual)
       END IF

     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE AdamsMoulton
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Apply constant timestep BDF integration scheme to elementwise matrix entry.
!------------------------------------------------------------------------------
   SUBROUTINE BDFLocal( N, dt, MassMatrix, StiffMatrix, &
                   Force, PrevSolution, Order )
!------------------------------------------------------------------------------

     INTEGER :: N, Order
     REAL(KIND=dp) :: Force(:),PrevSolution(:,:),dt
     REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:)


!------------------------------------------------------------------------------
     INTEGER :: i,j,NB1,NB2

     REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
     NB1 = SIZE(StiffMatrix,1)
     NB2 = SIZE(StiffMatrix,2)

     SELECT CASE( Order)
     CASE(1)
       DO i=1,NB1
         s = 0.0_dp
         DO j=1,N
           s = s + (1.0d0/dt)*MassMatrix(i,j)*PrevSolution(j,1)
         END DO
         Force(i) = Force(i) + s

         DO j=1,NB2
           StiffMatrix(i,j) = (1.0d0/dt)*MassMatrix(i,j) + StiffMatrix(i,j)
         END DO
       END DO

     CASE(2)
       DO i=1,NB1
         s = 0.0_dp
         DO j=1,N
             s = s + (1.0_dp/dt)*MassMatrix(i,j) * ( &
                   2.0_dp*PrevSolution(j,1) - 0.5_dp*PrevSolution(j,2) )
         END DO
         Force(i) = Force(i) + s

         DO j=1,NB2
           StiffMatrix(i,j) = (1.5d0/dt)*MassMatrix(i,j) + StiffMatrix(i,j)
         END DO
       END DO

     CASE(3)
       DO i=1,NB1
         s = 0.0_dp
         DO j=1,N
           s = s + (1.0d0/dt)*MassMatrix(i,j) * ( &
             3._dp*PrevSolution(j,1) &
             - (3._dp/2._dp)*PrevSolution(j,2) &
             + (1._dp/3._dp)*PrevSolution(j,3) )
         END DO
         Force(i) = Force(i) + s

         DO j=1,NB2
           StiffMatrix(i,j) = (11._dp/(6._dp*dt))*MassMatrix(i,j) &
              + StiffMatrix(i,j)
         END DO
       END DO

     CASE(4)
       DO i=1,NB1
         s = 0.0_dp
         DO j=1,N
           s = s + (1._dp/dt)*MassMatrix(i,j) * ( &
             4._dp*PrevSolution(j,1) &
             - 3._dp*PrevSolution(j,2) &
             + (4._dp/3._dp)*PrevSolution(j,3) &
             - (1._dp/4._dp)*PrevSolution(j,4) )
         END DO
         Force(i) = Force(i) + s

         DO j=1,NB2
           StiffMatrix(i,j) = (25.0_dp/(12._dp*dt))*MassMatrix(i,j) &
             + StiffMatrix(i,j)
         END DO
       END DO

     CASE(5)
       DO i=1,NB1
         s = 0.0_dp
         DO j=1,N
           s = s + (1._dp/dt)*MassMatrix(i,j) * ( &
             5._dp*PrevSolution(j,1)    &
             - 5._dp*PrevSolution(j,2)  &
             + (10._dp/3._dp)*PrevSolution(j,3)  &
             - (5._dp/4._dp)*PrevSolution(j,4)   &
             + (1._dp/5._dp)*PrevSolution(j,5) )
         END DO
         Force(i) = Force(i) + s

         DO j=1,NB2
           StiffMatrix(i,j) = (137._dp/(60._dp*dt)) * MassMatrix(i,j) &
              + StiffMatrix(i,j)
         END DO
       END DO

     CASE DEFAULT
        WRITE( Message, * ) 'Invalid order BDF', Order
        CALL Fatal( 'BDFLocal', Message )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE BDFLocal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply variable timestep BDF integration scheme to elementwise matrix entry.
!------------------------------------------------------------------------------
   SUBROUTINE VBDFLocal( N, dts, MassMatrix, StiffMatrix, &
                   Force, PrevSolution, Order )
!------------------------------------------------------------------------------

     INTEGER :: N, Order
     REAL(KIND=dp) :: Force(:),PrevSolution(:,:),dts(:)
     REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:)
     REAL(KIND=dp) :: a(4)

!------------------------------------------------------------------------------
     INTEGER :: i,j,k,NB1,NB2
     REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
     NB1 = SIZE(StiffMatrix,1)
     NB2 = SIZE(StiffMatrix,2)

     a = 0.0_dp

     a(1) = 1.0_dp / Dts(1)
     a(2) = -1.0_dp / Dts(1)
     IF(Order >= 2) THEN
       a(1) = a(1) + 1.0_dp / (Dts(1)+Dts(2)) 
       a(2) = a(2) - (1.0_dp + Dts(1)/Dts(2)) / (Dts(1)+Dts(2)) 
       a(3) = (Dts(1)/Dts(2)) / (Dts(1)+Dts(2)) 
     END IF
     IF(Order >= 3) THEN
       a(1) = a(1) + 1.0_dp / (Dts(1)+Dts(2)+Dts(3)) 
       a(2) = a(2) - (1.0_dp + Dts(1)/Dts(2)*(1.0+(Dts(1)+Dts(2))/(Dts(2)+Dts(3)))) / (Dts(1)+Dts(2)+Dts(3)) 
       a(3) = a(3) + (Dts(1)/Dts(2)*(1.0+(Dts(1)+Dts(2))/(Dts(2)+Dts(3))) + &
               Dts(1)/Dts(3)*(Dts(1)+Dts(2))/(Dts(2)+Dts(3)) ) / (Dts(1)+Dts(2)+Dts(3)) 
       a(4) = -(Dts(1)/Dts(3))*(Dts(1)+Dts(2))/(Dts(2)+Dts(3)) / (Dts(1)+Dts(2)+Dts(3))
     END IF
     IF(Order > 3) THEN
       CALL Warn('VBDFLocal','Variable timestep BDF implemented only to order 3')
     END IF

     DO i=1,NB1
       s = 0.0_dp
       DO k=1,MIN(Order,3)
         DO j=1,N
           s = s - a(k+1)*MassMatrix(i,j) * PrevSolution(j,k)
         END DO
       END DO 
       Force(i) = Force(i) + s

       DO j=1,NB2
         StiffMatrix(i,j) = a(1)*MassMatrix(i,j) + StiffMatrix(i,j)
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE VBDFLocal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply second order Bossok time integration scheme to the elementwise matrix
!> entry.
!------------------------------------------------------------------------------
   SUBROUTINE Bossak2ndOrder( N, dt, MassMatrix, DampMatrix, StiffMatrix, &
                        Force, X,V,A,Alpha )
!------------------------------------------------------------------------------

     INTEGER :: N
     REAL(KIND=dp) :: Force(:),X(:),V(:),A(:),dt
     REAL(KIND=dp) :: Alpha,Beta, Gamma
     REAL(KIND=dp) :: MassMatrix(:,:),DampMatrix(:,:),StiffMatrix(:,:)

!------------------------------------------------------------------------------
     INTEGER :: i,j,n1,n2

     REAL(KIND=dp) :: s, aa
!------------------------------------------------------------------------------
     n1 = MIN( N, SIZE(StiffMatrix,1) )
     n2 = MIN( N, SIZE(StiffMatrix,2) )

     Gamma = 0.5d0 - Alpha
     Beta = (1.0d0 - Alpha)**2 / 4.0d0
     DO i=1,n1
       s = 0.0d0
       DO j=1,n2
         s = s + ( (1.0d0 - Alpha) / (Beta*dt**2) ) * MassMatrix(i,j) * X(j)
         s = s + ( (1.0d0 - Alpha) / (Beta*dt)) * MassMatrix(i,j) * V(j)
         s = s - ( (1.0d0 - Alpha) * (1.0d0 - 1.0d0 / (2.0d0*Beta)) + Alpha )*&
                              MassMatrix(i,j) * A(j)

         s = s + ( Gamma / (Beta*dt) ) * DampMatrix(i,j) * X(j)
         s = s + ( Gamma/Beta - 1.0d0) * DampMatrix(i,j) * V(j)
         s = s - ((1.0d0 - Gamma) + Gamma * (1.0d0 - 1.0d0 / (2.0d0*Beta))) * &
                          dt * DampMatrix(i,j) * A(j)

         StiffMatrix(i,j) = StiffMatrix(i,j) +  &
           ( (1.0d0 - Alpha) / (Beta*dt**2) ) * MassMatrix(i,j) + &
                  (Gamma / (Beta*dt)) * DampMatrix(i,j)
       END DO 
       Force(i) = Force(i) + s
     END DO 
   END SUBROUTINE Bossak2ndOrder
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply fractional step integration scheme.
!------------------------------------------------------------------------------
   SUBROUTINE FractionalStep( N, dt, MassMatrix, StiffMatrix, &
                   Force, PrevSolution, Beta, Solver )
!------------------------------------------------------------------------------
     USE Types
     USE Lists

     TYPE(Solver_t) :: Solver

     INTEGER :: N
     REAL(KIND=dp) :: Force(:),PrevSolution(:),dt, fsstep, fsTheta, fsdTheta, &
                      fsAlpha, fsBeta, MassCoeff, ForceCoeff
     REAL(KIND=dp) :: MassMatrix(:,:),StiffMatrix(:,:),Beta

!------------------------------------------------------------------------------
     INTEGER :: i,j,NB

     REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
     NB = SIZE( StiffMatrix,1 )

! Check out the coefficients related to each step of the fractional stepping

     fsstep   = ListGetConstReal( Solver % Values, 'fsstep')
     fsTheta  = ListGetConstReal( Solver % Values, 'fsTheta')
     fsdTheta = ListGetConstReal( Solver % Values, 'fsdTheta')
     fsAlpha  = ListGetConstReal( Solver % Values, 'fsAlpha')
     fsBeta   = ListGetConstReal( Solver % Values, 'fsBeta')

     SELECT CASE( INT(fsstep) )     
       CASE(1)
        MassCoeff = fsAlpha * fsTheta
        ForceCoeff = fsBeta * fsTheta 
       CASE(2)
        MassCoeff = fsBeta * fsdTheta
        ForceCoeff = fsAlpha * fsdTheta
       CASE(3)
        MassCoeff = fsAlpha * fsTheta
        ForceCoeff = fsBeta * fsTheta
     END SELECT

     DO i=1,NB
       s = 0.0d0
       DO j=1,N
         s = s + (1.0d0/dt) * MassMatrix(i,j) * PrevSolution(j) - &
              ForceCoeff * StiffMatrix(i,j) * PrevSolution(j)
       END DO
       Force(i) = Force(i) + s

       DO j=1,NB
           StiffMatrix(i,j) = MassCoeff * StiffMatrix(i,j) + (1.0d0/dt)*MassMatrix(i,j)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE FractionalStep
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply 2nd order Newmark time integration scheme.
!------------------------------------------------------------------------------
   SUBROUTINE Newmark2ndOrder( N, dt, MassMatrix, DampMatrix, StiffMatrix, &
                        Force, PrevSol0,PrevSol1, Avarage )
!------------------------------------------------------------------------------

     INTEGER :: N
     REAL(KIND=dp) :: Force(:),PrevSol0(:),PrevSol1(:),dt
     LOGICAL :: Avarage
     REAL(KIND=dp) :: MassMatrix(:,:),DampMatrix(:,:),StiffMatrix(:,:)

!------------------------------------------------------------------------------
     INTEGER :: i,j

     REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
     IF ( Avarage ) THEN 
       DO i=1,N
         s = 0.0d0
         DO j=1,N
           s = s - ((1/dt**2)*MassMatrix(i,j) - (1/(2*dt))*DampMatrix(i,j) + &
                       StiffMatrix(i,j) / 3 ) * PrevSol0(j)

           s = s + ((2/dt**2)*MassMatrix(i,j) - StiffMatrix(i,j) / 3) *  &
                               PrevSol1(j)

           StiffMatrix(i,j) = StiffMatrix(i,j) / 3 +  &
                        (1/dt**2)*MassMatrix(i,j) + (1/(2*dt))*DampMatrix(i,j)
         END DO
         Force(i) = Force(i) + s
       END DO
     ELSE 
       DO i=1,N
         s = 0.0d0
         DO j=1,N
           s = s - ((1/dt**2)*MassMatrix(i,j) - (1/(2*dt))*DampMatrix(i,j)) * &
                                     PrevSol0(j)

           s = s + (2/dt**2)*MassMatrix(i,j) * PrevSol1(j)

           StiffMatrix(i,j) = StiffMatrix(i,j) +  &
                        (1/dt**2)*MassMatrix(i,j) + (1/(2*dt))*DampMatrix(i,j)
         END DO
         Force(i) = Force(i) + s
       END DO
     END IF
   END SUBROUTINE Newmark2ndOrder
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
! These are similar subroutines as above except they operate on the full
! matrix level and assumes that A % Values includes stiffness matrix
! and A % MassValues includes the mass matrix.
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE RungeKutta_CRS( dt, Matrix, Force, PrevSolution, CurrSolution )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: Force(:)
    REAL(KIND=dp) :: PrevSolution(:),CurrSolution(:),dt
    TYPE(Matrix_t), POINTER :: Matrix
!------------------------------------------------------------------------------
    INTEGER :: i,j,n
    REAL(KIND=dp), POINTER :: Stiff(:),Mass(:),MassL(:)
    REAL(KIND=dp), POINTER :: PrevResidual(:)
    INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp) :: su,mu,uj,ui
!------------------------------------------------------------------------------
    
    n = Matrix % NumberOfRows
    Rows   => Matrix % Rows
    Cols   => Matrix % Cols
    Stiff => Matrix % Values
    Mass => Matrix % MassValues
    
    ! For true 2nd order accuracy of Crank-Nicolsen for nonlinear
    ! problems the old matrix and rhs need to be considered.
    !-------------------------------------------------------------
    IF( ASSOCIATED( Matrix % BulkResidual ) ) THEN
      ! Ax-f from previous iteration
      PrevResidual => Matrix % BulkResidual

      IF( ASSOCIATED( Matrix % MassValuesLumped ) ) THEN
        MassL => Matrix % MassValuesLumped
        !$omp parallel do private(j,ui)
        DO i=1,n
          j = Matrix % Diag(i)
          ui = PrevSolution(i)
          Force(i) = PrevResidual(i) + (1.0/dt) * MassL(i) * ui
        END DO
        !$omp end parallel do
        Stiff = 0
        Stiff(Matrix % Diag) = Stiff(Matrix % Diag) + (1.0d0/dt) * MassL(i)
      ELSE
        !$omp parallel do private(j,uj,mu)     
        DO i=1,n
          mu = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            uj = PrevSolution(Cols(j))
            mu = mu + Mass(j) * uj
          END DO
          Force(i) = PrevResidual(i) + (1.0/dt)*mu
        END DO
        !$omp end parallel do
        Stiff = (1.0d0/dt) * Mass
      END IF
    ELSE

      IF( ASSOCIATED( Matrix % MassValuesLumped ) ) THEN
        MassL => Matrix % MassValuesLumped
        !$omp parallel do private(j,uj,mu)
        DO i=1,n
          su = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            uj = CurrSolution(Cols(j))
            su = su + Stiff(j) * uj
          END DO
          mu = MassL(i)*PrevSolution(i)
          Force(i) = Force(i) - su + (1.0/dt)*mu
        END DO
        !$omp end parallel do
        Stiff = 0
        Stiff(Matrix % Diag) = Stiff(Matrix % Diag) + (1.0/dt) * MassL
      ELSE
        !$omp parallel do private(j,uj,su,mu)     
        DO i=1,n
          su = 0.0_dp
          mu = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            uj = CurrSolution(Cols(j))
            su = su + Stiff(j) * uj
            uj = PrevSolution(Cols(j))
            mu = mu + Mass(j) * uj
          END DO
          Force(i) = Force(i) - su + (1.0/dt)*mu
        END DO
        !$omp end parallel do
        Stiff = (1.0d0/dt) * Mass
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE RungeKutta_CRS
!------------------------------------------------------------------------------

!> Apply newmark beta scheme to CRS_Matrix.
!> This is used also for Implicit Euler with beta=1.0 and
!> Crank-Nicolson with beta=0.5.
!------------------------------------------------------------------------------
   SUBROUTINE NewmarkBeta_CRS( dt, Matrix, Force, PrevSolution, Beta )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: Force(:)
    REAL(KIND=dp) :: PrevSolution(:),dt
    REAL(KIND=dp) :: Beta
    TYPE(Matrix_t), POINTER :: Matrix
!------------------------------------------------------------------------------
    INTEGER :: i,j,n
    REAL(KIND=dp), POINTER :: Stiff(:),Mass(:),MassL(:)
    REAL(KIND=dp), POINTER :: PrevResidual(:)
    INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp) :: su,mu,uj,ui
!------------------------------------------------------------------------------
    
    n = Matrix % NumberOfRows
    Rows   => Matrix % Rows
    Cols   => Matrix % Cols
    Stiff => Matrix % Values
    Mass => Matrix % MassValues
    
    ! For true 2nd order accuracy of Crank-Nicolsen for nonlinear
    ! problems the old matrix and rhs need to be considered.
    !-------------------------------------------------------------
    IF( ASSOCIATED( Matrix % BulkResidual ) ) THEN
      ! Ax-f from previous iteration
      PrevResidual => Matrix % BulkResidual

      IF( ASSOCIATED( Matrix % MassValuesLumped ) ) THEN
        MassL => Matrix % MassValuesLumped
        !$omp parallel do private(j,ui)
        DO i=1,n
          j = Matrix % Diag(i)
          ui = PrevSolution(i)
          Force(i) = Beta * Force(i) - (1-Beta) * PrevResidual(i) + &
              (1.0/dt) * MassL(i) * ui
        END DO
        !$omp end parallel do
        Stiff = Beta * Stiff
        Stiff(Matrix % Diag) = Stiff(Matrix % Diag) + (1.0d0/dt) * MassL(i)
      ELSE
        !$omp parallel do private(j,uj,mu)     
        DO i=1,n
          mu = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            uj = PrevSolution(Cols(j))
            mu = mu + Mass(j) * uj
          END DO
          Force(i) = Beta * Force(i) - (1-Beta) * PrevResidual(i) + &
              (1.0/dt)*mu
        END DO
        !$omp end parallel do
        Stiff = Beta * Stiff + (1.0d0/dt) * Mass
      END IF
    ELSE

      IF( ASSOCIATED( Matrix % MassValuesLumped ) ) THEN
        MassL => Matrix % MassValuesLumped
        !$omp parallel do private(j,uj,mu)
        DO i=1,n
          su = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            uj = PrevSolution(Cols(j))
            su = su + Stiff(j) * uj
          END DO
          mu = MassL(i)*PrevSolution(i)
          Force(i) = Force(i) - (1-Beta)*su + (1.0/dt)*mu
        END DO
        !$omp end parallel do
        Stiff = Beta * Stiff
        Stiff(Matrix % Diag) = Stiff(Matrix % Diag) + (1.0/dt) * MassL
      ELSE
        !$omp parallel do private(j,uj,su,mu)     
        DO i=1,n
          su = 0.0_dp
          mu = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            uj = PrevSolution(Cols(j))
            su = su + Stiff(j) * uj
            mu = mu + Mass(j) * uj
          END DO
          Force(i) = Force(i) - (1-Beta)*su + (1.0/dt)*mu
        END DO
        !$omp end parallel do
        Stiff = Beta * Stiff + (1.0d0/dt) * Mass
      END IF
    END IF


!------------------------------------------------------------------------------
  END SUBROUTINE NewmarkBeta_CRS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply constant timestep BDF integration scheme to CRS Matrix.
!------------------------------------------------------------------------------
  SUBROUTINE BDF_CRS( dt, Matrix, Force, PrevSolution, Order )
!------------------------------------------------------------------------------

    TYPE(Matrix_t), POINTER :: Matrix
    INTEGER :: Order
    REAL(KIND=dp) :: Force(:),PrevSolution(:,:),dt
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n
    REAL(KIND=dp), POINTER :: Stiff(:),Mass(:),MassL(:)
    INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp) :: su,mu,uj,ui,a(6)
!------------------------------------------------------------------------------
    
    n = Matrix % NumberOfRows
    Rows   => Matrix % Rows
    Cols   => Matrix % Cols
    Stiff => Matrix % Values
    Mass => Matrix % MassValues
    IF( Matrix % Lumped ) THEN
      MassL => Matrix % MassValuesLumped
    END IF

    a = 0.0_dp
    SELECT CASE( Order)
    CASE(1)
      a(1) = 1.0_dp
      a(2) = -1.0_dp
    CASE(2)
      a(1) = 1.5_dp
      a(2) = -2.0_dp
      a(3) = 0.5_dp
    CASE(3)
      a(1) = 11._dp/6_dp
      a(2) = -3._dp
      a(3) = 3._dp/2._dp
      a(4) = -1._dp/3._dp
    CASE(4)
      a(1) = 25._dp/12_dp
      a(2) = -4._dp
      a(3) = 3._dp
      a(4) = -4._dp/3._dp
      a(5) = 1._dp/4._dp
    CASE(5)
      a(1) = 137._dp/60_dp
      a(2) = -5._dp
      a(3) = 5._dp
      a(4) = -10._dp/3._dp
      a(5) = 5._dp/4._dp
      a(6) = -1._dp/5._dp
    CASE DEFAULT
      CALL Fatal('BDF_CRS','Constant timestep BDF implemented only to order 5')
    END SELECT

    a = a / dt

    IF( Matrix % Lumped ) THEN
      !$omp parallel do private(j,k,ui,mu)
      DO i=1,n
        DO k=1,Order
          ui = PrevSolution(i,k)
          Force(i) = Force(i) - MassL(i) * a(k+1) * ui
        END DO
        j = Matrix % Diag(i)      
        Stiff(j) = Stiff(j) + a(1) * MassL(i)
      END DO
      !$omp end parallel do
    ELSE
      !$omp parallel do private(j,k,uj,mu)
      DO i=1,n
        mu = 0.0_dp
        DO j=Rows(i),Rows(i+1)-1
          DO k=1,Order
            uj = PrevSolution(Cols(j),k)
            mu = mu - Mass(j) * a(k+1) * uj  
          END DO
          Stiff(j) = Stiff(j) + a(1) * Mass(j)
        END DO
        Force(i) = Force(i) + mu
      END DO
      !$omp end parallel do
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE BDF_CRS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply variable timestep BDF integration scheme to CRS matrix.
!------------------------------------------------------------------------------
   SUBROUTINE VBDF_CRS( dts, Matrix, Force, PrevSolution, Order )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: Matrix
     INTEGER :: Order
     REAL(KIND=dp) :: Force(:),PrevSolution(:,:),dts(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n
     REAL(KIND=dp), POINTER :: Stiff(:),Mass(:),MassL(:)
     INTEGER, POINTER :: Cols(:),Rows(:)
     REAL(KIND=dp) :: su,mu,uj,ui,a(4)
!------------------------------------------------------------------------------
    
     n = Matrix % NumberOfRows
     Rows   => Matrix % Rows
     Cols   => Matrix % Cols
     Stiff => Matrix % Values
     Mass => Matrix % MassValues
     IF( Matrix % Lumped ) THEN
       MassL => Matrix % MassValuesLumped
     END IF
          
     a = 0.0_dp

     a(1) = 1.0_dp / Dts(1)
     a(2) = -1.0_dp / Dts(1)
     IF(Order >= 2) THEN
       a(1) = a(1) + 1.0_dp / (Dts(1)+Dts(2)) 
       a(2) = a(2) - (1.0_dp + Dts(1)/Dts(2)) / (Dts(1)+Dts(2)) 
       a(3) = (Dts(1)/Dts(2)) / (Dts(1)+Dts(2)) 
     END IF
     IF(Order >= 3) THEN
       a(1) = a(1) + 1.0_dp / (Dts(1)+Dts(2)+Dts(3)) 
       a(2) = a(2) - (1.0_dp + Dts(1)/Dts(2)*(1.0+(Dts(1)+Dts(2))/(Dts(2)+Dts(3)))) / (Dts(1)+Dts(2)+Dts(3)) 
       a(3) = a(3) + (Dts(1)/Dts(2)*(1.0+(Dts(1)+Dts(2))/(Dts(2)+Dts(3))) + &
               Dts(1)/Dts(3)*(Dts(1)+Dts(2))/(Dts(2)+Dts(3)) ) / (Dts(1)+Dts(2)+Dts(3)) 
       a(4) = -(Dts(1)/Dts(3))*(Dts(1)+Dts(2))/(Dts(2)+Dts(3)) / (Dts(1)+Dts(2)+Dts(3))
     END IF

     ! The constant timestep BDF is implemented up to 5th degree hence higher order
     ! schemes are possible even though this scheme cannot handle them. 
     IF(Order > 3) THEN
       CALL Warn('VBDF_CRS','Variable timestep BDF implemented only to order 3')
     END IF

     IF( Matrix % Lumped ) THEN
      !$omp parallel do private(j,k,ui,mu)
       DO i=1,n
         DO k=1,MIN(Order,3)
           ui = PrevSolution(i,k)
           Force(i) = Force(i) - MassL(i) * a(k+1) * ui
         END DO
         j = Matrix % Diag(i)      
         Stiff(j) = Stiff(j) + a(1) * MassL(i)
       END DO
       !$omp end parallel do
     ELSE
       !$omp parallel do private(j,k,uj,mu)
       DO i=1,n
         mu = 0.0_dp
         DO j=Rows(i),Rows(i+1)-1
           DO k=1,MIN(Order,3)
             uj = PrevSolution(Cols(j),k)
             mu = mu - Mass(j) * a(k+1) * uj  
           END DO
         END DO
         Force(i) = Force(i) + mu
       END DO
       !$omp end parallel do
       Stiff = Stiff + a(1) * Mass
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE VBDF_CRS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply fractional step integration scheme for a CRS matrix equation.
!------------------------------------------------------------------------------
   SUBROUTINE FractionalStep_CRS( dt, Matrix, Force, PrevSolution, Solver )
!------------------------------------------------------------------------------
     USE Types
     USE Lists

     TYPE(Solver_t) :: Solver    
     TYPE(Matrix_t), POINTER :: Matrix
     REAL(KIND=dp) :: Force(:),PrevSolution(:),dt
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: fsstep, fsTheta, fsdTheta, &
                      fsAlpha, fsBeta, MassCoeff, ForceCoeff
     REAL(KIND=dp),POINTER :: Mass(:),Stiff(:),MassL(:)
     INTEGER :: n, i,j
     INTEGER, POINTER :: Rows(:),Cols(:)
     REAL(KIND=dp) :: uj,ui,su,mu
!------------------------------------------------------------------------------

! Check out the coefficients related to each step of the fractional stepping

     fsstep   = ListGetConstReal( Solver % Values, 'fsstep')
     fsTheta  = ListGetConstReal( Solver % Values, 'fsTheta')
     fsdTheta = ListGetConstReal( Solver % Values, 'fsdTheta')
     fsAlpha  = ListGetConstReal( Solver % Values, 'fsAlpha')
     fsBeta   = ListGetConstReal( Solver % Values, 'fsBeta')

     SELECT CASE( INT(fsstep) )     
       CASE(1)
        MassCoeff = fsAlpha * fsTheta
        ForceCoeff = fsBeta * fsTheta 
       CASE(2)
        MassCoeff = fsBeta * fsdTheta
        ForceCoeff = fsAlpha * fsdTheta
       CASE(3)
        MassCoeff = fsAlpha * fsTheta
        ForceCoeff = fsBeta * fsTheta
     END SELECT

     n = Matrix % NumberOfRows
     Rows   => Matrix % Rows
     Cols   => Matrix % Cols
     Stiff => Matrix % Values
     Mass => Matrix % MassValues
 

     IF( ASSOCIATED( Matrix % MassValuesLumped ) ) THEN
       MassL => Matrix % MassValuesLumped
       !$omp parallel do private(j,uj,ui,su)
       DO i=1,n
         su = 0.0_dp
         DO j=Rows(i),Rows(i+1)-1
           uj = PrevSolution(Cols(j))
           su = su + Stiff(j) * uj
         END DO
         j = Matrix % Diag(i)
         ui = PrevSolution(Cols(j))
         
         Force(i) = Force(i) - ForceCoeff * su + (1._dp/dt)*MassL(i)*ui
         Stiff(j) = MassCoeff * Stiff(j) + (1.0d0/dt) * MassL(i)
       END DO
       !$omp end parallel do
     ELSE
       !$omp parallel do private(j,uj,su,mu)     
       DO i=1,n
         su = 0.0_dp
         mu = 0.0_dp
         DO j=Rows(i),Rows(i+1)-1
           uj = PrevSolution(Cols(j))
           su = su + Stiff(j) * uj
           mu = mu + Mass(j) * uj
         END DO
         Force(i) = Force(i) - ForceCoeff * su + (1.0/dt)*mu
       END DO
       !$omp end parallel do
       Stiff = MassCoeff * Stiff + (1.0d0/dt) * Mass
     END IF

!------------------------------------------------------------------------------
   END SUBROUTINE FractionalStep_CRS
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
END MODULE TimeIntegrate
!------------------------------------------------------------------------------

!> \}
