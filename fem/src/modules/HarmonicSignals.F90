!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Author:    Eelis Takala(Trafotek Oy) 
! *  Email:     eelis.takala@trafotek.fi
! *  Web:       http://www.trafotek.fi
! *  Addresses: Trafotek Oy
! *             Kaarinantie 700
! *             Turku
! *
! *  Original Date: May 2017
! *
! *****************************************************************************/
 

MODULE HarmUtils
  IMPLICIT NONE
  
  CONTAINS
    
    !------------------------------------------------------------------
    ! This function returns the harmonic sum of a list of 
    ! sinus curves: sum_k( A_k * sin(k * x))
    !
    ! Args:
    !   fundamental_f       :: fundamental frequency
    !   amplitudes(:, 2)    :: list of harmonic numbers (first column)
    !                          and their amplitudes (second column)
    !   time                :: time
    !
    ! Returns:
    !   
    ! -----------------------------------------------------------------
    FUNCTION SinSum(fundamental_f, amplitudes, t, phase) RESULT(sumA)
      USE DefUtils
      IMPLICIT NONE
      
      REAL(KIND=dp) :: fundamental_f, fundamental_omega, t
      REAL(KIND=dp) :: amplitudes(:,:)
      REAL(KIND=dp) :: sumA
      REAL(KIND=dp) :: phase
      
      INTEGER :: i
      
      sumA = 0
      fundamental_omega = 2 * PI * fundamental_f
      DO i = 1, SIZE(amplitudes(:,1))
        sumA = sumA + amplitudes(i, 2) * sqrt(2._dp) * sin(amplitudes(i,1)  * fundamental_omega * t + phase)
      END DO
      
    END FUNCTION SinSum
    
    FUNCTION CosSum(fundamental_f, amplitudes, t, phase) RESULT(sumA)
      USE DefUtils
      IMPLICIT NONE
      
      REAL(KIND=dp) :: fundamental_f, fundamental_omega, t
      REAL(KIND=dp) :: amplitudes(:,:)
      REAL(KIND=dp) :: sumA
      REAL(KIND=dp) :: phase
      
      INTEGER :: i
      
      sumA = 0
      fundamental_omega = 2 * PI * fundamental_f
      DO i = 1, SIZE(amplitudes(:,1))
        sumA = sumA + amplitudes(i, 2) * sqrt(2._dp) * cos(amplitudes(i,1)  * fundamental_omega * t + phase)
      END DO
      
    END FUNCTION CosSum
    
END MODULE HarmUtils


!!!!!!!!!!!!!!!!! ALL PHASES
FUNCTION source( model, n, time ) RESULT(current)
  USE DefUtils
  USE HarmUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  TYPE(ValueList_t), POINTER :: BF
  LOGICAL :: Found, amplitude_fade_in
  INTEGER :: n, i
  REAL(KIND=dp) :: time, current, Isum, freq, coeff, phase
  REAL(KIND=dp), POINTER :: Amplitudes(:,:)=>NULL()
  
  BF => Model % BodyForces(1) % Values
  
  CALL ListPushNameSpace(TRIM(ListGetActiveName())//':')
  
  freq = GetConstReal(BF, 'Frequency', Found)
  IF (.NOT. FOUND) CALL FATAL('source', ListGetActiveName()//': Frequency not found in Body Force 1 section.')
  
  CALL GetConstRealArray( BF, Amplitudes, 'Harmonic Content', Found)
  IF (.NOT. FOUND) CALL FATAL('source', ListGetActiveName()//': Harmonic Content not found in Body Force 1 section.')
  
  phase = GetConstReal(BF, 'Phase', Found)
  IF (.NOT. FOUND) CALL FATAL('source', ListGetActiveName()//': phase not found in Body Force 1 section.')
  
  amplitude_fade_in = GetLogical(BF, 'Amplitude Fade In', Found)
  IF (.NOT. Found) amplitude_fade_in = .TRUE.
  
!  Isum = CosSum(freq, Amplitudes, time, phase)
  Isum = SinSum(freq, Amplitudes, time, phase)
  
  DO i = 1, SIZE(amplitudes(:,1))
    WRITE(Message,*) TRIM(ListGetActiveName())//': Amplitudes = ', Amplitudes(i,1), ' ', Amplitudes(i,2)
    CALL Info('source', Message, Level=5 )
  END DO
  WRITE(Message,*) TRIM(ListGetActiveName())//': Current = ', Isum
  CALL Info('source', Message, Level=5 )
  WRITE(Message,*) TRIM(ListGetActiveName())//': Frequency = ', freq
  CALL Info('source', Message, Level=5 )
  WRITE(Message,*) TRIM(ListGetActiveName())//': Phase = ', phase/PI*180._dp
  CALL Info('source', Message, Level=5 )
 
  CALL ListPopNameSpace()
  
  coeff = 1._dp
  IF (amplitude_fade_in .AND. 2*pi*freq*time <= pi/2) coeff = 4._dp*freq*time 
  
  CALL Info('source', Message, Level=5 )
  
  current = coeff * Isum
 
END FUNCTION source
