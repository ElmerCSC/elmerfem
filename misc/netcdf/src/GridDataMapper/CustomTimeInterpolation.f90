!------------------------------------------------------------------------------
! Vili Forsell
! Created: 15.6.2011
! Last Modified: 16.6.2011
!------------------------------------------------------------------------------
! This module is the customizable file to implement optional functions for
! time interpolation between two points.
! HOW TO USE:
! - add a case calling your interpolation function (can be implemented within this file)
! - add a matching string into Solver Input File in variable "Time Interpolation Method"
! - uses linear interpolation by default
!------------------------------------------------------------------------------
MODULE CustomTimeInterpolation
  USE DefUtils, ONLY: dp, MAX_NAME_LEN
  IMPLICIT NONE

  ! If different parameters needed, adjust GridDataMapper.f90 interpolation call and add your function name here
!  INTERFACE ChooseTimeInterpolation
!    MODULE PROCEDURE ChooseTimeInterpolation
!  END INTERFACE

  PRIVATE :: SeasonalInterpolation

  CONTAINS

     !--------------- ChooseTimeInterpolation() -----------
     !--- Chooses the time interpolation method on basis of variable "Time Interpolation Method" in the SIF
     FUNCTION ChooseTimeInterpolation( time, left, right, method, give_info ) RESULT( NetCDF_value )
     !-----------------------------------------------------

       USE NetCDFInterpolate, ONLY: LinearInterpolation
       USE DefUtils, ONLY: dp,MAX_NAME_LEN
       USE Messages
       IMPLICIT NONE

       !--- Input parameters
       REAL(KIND=dp), INTENT(IN) :: time, & ! The time index as a real value
                                    left(:), & ! The closest left NetCDF point of form (time index, NetCDF value)
                                    right(:) ! The closest right NetCDF point of form (time index, NetCDF value)
       CHARACTER(*), INTENT(IN) :: method ! The string used for selecting the method
       ! It's contents come from the SIF variable "Time Interpolation Method"
       LOGICAL, INTENT(IN) :: give_info

       !--- Output parameter
       REAL(KIND=dp) :: NetCDF_value ! The NetCDF value obtained on basis of interpolating in terms of time

       !--- Checks
       IF ( size(left) .NE. size(right) ) CALL Fatal('CustomTimeInterpolation','Input vector sizes do not match')

       !--- Selects the implemented method
       SELECT CASE (method) ! Selects the chosen interpolation method
         CASE ("seasonal")
           IF ( give_info ) WRITE (Message,'(A)') 'Seasonal interpolation in effect'
           NetCDF_value = SeasonalInterpolation(time,left,right)
         CASE ("linear")
           IF ( give_info ) WRITE (Message,'(A)') 'Linear interpolation in effect'
           NetCDF_value = LinearInterpolation(time,left,right)
         CASE DEFAULT ! Defaults to linear interpolation
           IF ( give_info ) WRITE (Message,'(A)') 'Default time interpolation in effect (linear)'
           NetCDF_value = LinearInterpolation(time,left,right)
       END SELECT

       IF ( give_info ) CALL Info('ChooseTimeInterpolation',Message)

     END FUNCTION ChooseTimeInterpolation

     !-------------- SinusoidalInterpolation() -------------'
     ! An example time interpolation function, 
     ! which applies sinusoidal behaviour between the two input points
     FUNCTION SeasonalInterpolation(time,left,right) RESULT( val )
     !------------------------------------------------------

       ! --- PREFACE ---
       ! Assume NetCDF data (and all grid points) are June temperatures and we want to consider also winter.
       ! First, we fit a sinusoid between the interval (with NetCDF value constant) to get the general drop in temperature.
       ! Second, we adjust the sinusoid values with the value difference between the two grid points
       ! to take linear (assumption) temperature change between Junes into account.
       ! Hence, we have a general interpolation for temperatures between the June seasons.
       USE DefUtils, ONLY: PI
       USE Messages ! Elmer messaging interface for error messages, etc.
       IMPLICIT NONE
       REAL(KIND=dp), INTENT(IN) :: time, left(:), right(:) ! Input variables
       REAL(KIND=dp) :: val ! Output variable
       REAL(KIND=dp) :: frequency, amplitude, phase, t ! Sinusoid parameters
       INTEGER :: alloc_stat
       REAL(KIND=dp), ALLOCATABLE :: diff(:) ! For allocating memory on basis of input size

       ! Allocate memory on basis of input and check that it didn't run out
       ALLOCATE ( diff(size(left)), STAT = alloc_stat )
       IF ( alloc_stat .NE. 0 ) THEN
         CALL Fatal( 'Sinusoidal Interpolation', 'Memory ran out' ) ! This automatically terminates the program
       END IF

       diff(:) = right(:) - left(:) ! For time and value differences in the uniform NetCDF grid

       ! 1) Fit a sinusoid between the interval:
       phase = PI ! The part from pi to 2*pi is below the y axis, so pi ~ t_left and 2pi ~ t_right
       frequency = 0.5 ! Half a sine per interval
       t = (time - left(1))/diff(1) ! time is between left and right, and t is where on sine we are, between 0 and 1
       amplitude = 30 ! Assume 30 degree difference between june and december

       val = amplitude*SIN(2*PI*t*frequency + phase)

       ! 2) Adjust with the value difference at current time point and we're ready
       val = val + t*diff(2)
!       WRITE (*,*) 'Value ', val

     END FUNCTION SeasonalInterpolation

END MODULE CustomTimeInterpolation
