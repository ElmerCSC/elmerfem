!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                           Main prog
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
PROGRAM sico2elmer
  IMPLICIT NONE
      
  INTEGER :: &
       noption=999,  i, j, k, flag = 1, ndata, istat, imax,jmax,kcmax,ktmax,krmax, gotit
  REAL ::&
       deform, Dx=0.0e0, time_gr, H_R_gr, YEAR_SEC = 31556926.0
  CHARACTER :: &
       runname * 5, ergnum * 5, rname * 5, enum * 5, grid
  LOGICAL :: &
       dataread=.FALSE.
  INTEGER, ALLOCATABLE ::&
       maske_gr(:,:), n_cts_gr(:,:)
  REAL, ALLOCATABLE :: &
       xi_gr(:), eta_gr(:),&
       zs_gr(:,:), zm_gr(:,:),zb_gr(:,:), zb0_gr(:,:),&
       z_c_gr(:,:,:), z_t_gr(:,:,:), z_r_gr(:,:,:),&
       H_c_gr(:,:), H_t_gr(:,:),&
       vx_c_gr(:,:,:), vy_c_gr(:,:,:), vz_c_gr(:,:,:),&
       vx_t_gr(:,:,:), vy_t_gr(:,:,:), vz_t_gr(:,:,:),&
       temp_c_gr(:,:,:), temp_r_gr(:,:,:), temp_t_gr(:,:,:),&
       temph_c_gr(:,:,:), temph_t_gr(:,:,:),&
       age_c_gr(:,:,:), age_t_gr(:,:,:), omega_t_gr(:,:,:),am_perp_gr(:,:),&
       qx_gr(:,:), qy_gr(:,:), Q_bm_gr(:,:), Q_tld_gr(:,:)

  SAVE dataread, runname

! inquire run identifyer and allocate fields arcording to information in log-file
!--------------------------------------------------------------------------------
  WRITE( *, '(A)', ADVANCE = 'YES')' '
  WRITE( *, '(A)', ADVANCE = 'YES')'This is sico2elmer'
  WRITE( *, '(A)', ADVANCE = 'YES')'******************'
  WRITE( *, '(A)', ADVANCE = 'YES')' '
  WRITE( *, '(A)', ADVANCE = 'NO') 'Input of run identifyer (5 characters): '
  READ (*,'(A5)',ADVANCE = 'YES') runname
  WRITE( *, '(A,A)', ADVANCE = 'YES') 'chosen run identifyer: ', runname
  CALL readlog_c(runname, imax,jmax,kcmax,ktmax,krmax,deform,Dx,gotit)
  IF (gotit .EQ. 0) THEN
     WRITE(6,'(A)', ADVANCE = 'YES') 'Error occurred while opening log-file!'
     STOP
  ELSE
     WRITE( *, '(A)', ADVANCE = 'YES') 'Log read successful'
!     WRITE( *, '(I4,I4,I4,I4,F8.4,F8.4)', ADVANCE = 'YES') imax,jmax,kcmax,ktmax,deform,dx
  END IF
  ALLOCATE( &
       maske_gr(0:imax,0:jmax), n_cts_gr(0:imax,0:jmax),&
       xi_gr(0:imax), eta_gr(0:jmax),&
       z_c_gr(0:imax,0:jmax,0:kcmax),&
       z_t_gr(0:imax,0:jmax,0:ktmax),&
       z_r_gr(0:imax,0:jmax,0:krmax),&
       zs_gr(0:imax,0:jmax), zm_gr(0:imax,0:jmax),&
       zb_gr(0:imax,0:jmax), zb0_gr(0:imax,0:jmax),&
       H_c_gr(0:imax,0:jmax),H_t_gr(0:imax,0:jmax),&
       vx_c_gr(0:imax,0:jmax,0:kcmax),&
       vy_c_gr(0:imax,0:jmax,0:kcmax),&
       vz_c_gr(0:imax,0:jmax,0:kcmax),&
       vx_t_gr(0:imax,0:jmax,0:ktmax),&
       vy_t_gr(0:imax,0:jmax,0:ktmax),&
       vz_t_gr(0:imax,0:jmax,0:ktmax),&
       temp_c_gr(0:imax,0:jmax,0:kcmax),&
       temph_c_gr(0:imax,0:jmax,0:kcmax),&
       age_c_gr(0:imax,0:jmax,0:kcmax),&
       omega_t_gr(0:imax,0:jmax,0:ktmax),&
       temp_t_gr(0:imax,0:jmax,0:ktmax),&
       temph_t_gr(0:imax,0:jmax,0:ktmax),&
       age_t_gr(0:imax,0:jmax,0:ktmax),&
       temp_r_gr(0:imax,0:jmax,0:krmax),&
       qx_gr(0:imax,0:jmax), qy_gr(0:imax,0:jmax),&
       Q_bm_gr(0:imax,0:jmax), Q_tld_gr(0:imax,0:jmax),&
       am_perp_gr(0:imax,0:jmax),&
       STAT=istat )
  
  IF ( istat /= 0 ) THEN
     WRITE( *, '(A)', ADVANCE = 'YES') 'Error in allocation of memory!'
     STOP
  ELSE
     WRITE( *, '(A)', ADVANCE = 'YES') 'Allocation of arrays done'
  END IF
  
! Main Menu
!--------------------------------------------------------------------------------
  DO 
     WRITE( *, '(A)', ADVANCE = 'YES') ' '
     WRITE( *, '(A)', ADVANCE = 'YES') ' '
     WRITE( *, '(A)', ADVANCE = 'YES') 'Options:' 
     WRITE( *, '(A)', ADVANCE = 'YES') ' ' 
     WRITE( *, '(A)', ADVANCE = 'YES') ' (1) Read SICOPOLIS timeslice-file'
     WRITE( *, '(A)', ADVANCE = 'YES') ' (2) Output of timestep-grid for ELMERPOST'
     WRITE( *, '(A)', ADVANCE = 'YES') ' (3) Output of timestep-data for ELMERPOST'
     WRITE( *, '(A)', ADVANCE = 'YES') ' (4) Output of timestep for input to ELMER SOLVER' 
     WRITE( *, '(A)', ADVANCE = 'YES') ' (5) Output of timestep-data in ASCII-format' 
     WRITE( *, '(A)', ADVANCE = 'YES') ' ' 
     WRITE( *, '(A)', ADVANCE = 'YES') ' (0) Quit'
     WRITE( *, '(A)', ADVANCE = 'YES') ' '
     WRITE( *, '(A)', ADVANCE = 'NO') 'Your choice: '
     READ ( *, '(I1)', ADVANCE = 'YES') noption
     WRITE( *, '(A)', ADVANCE = 'YES') ' '
     WRITE( *, '(A)', ADVANCE = 'YES') ' '
     
     IF ((noption == 1).and.(.NOT.dataread)) THEN
        WRITE( *, '(A,A)', ADVANCE = 'YES') 'Reading timeslice for run: ', runname
        CALL ReadData(runname, ergnum, &
             imax,jmax,kcmax,ktmax,krmax,deform,& 
             maske_gr, n_cts_gr,&
             time_gr, xi_gr, eta_gr, zs_gr, zm_gr,zb_gr, zb0_gr,&
             z_c_gr, z_t_gr, z_r_gr,&
             H_c_gr, H_t_gr, H_R_gr, vx_c_gr, vy_c_gr, vz_c_gr,&
             vx_t_gr, vy_t_gr, vz_t_gr, temp_c_gr, temp_r_gr, temp_t_gr,&
             temph_c_gr, temph_t_gr, am_perp_gr, age_c_gr, age_t_gr, omega_t_gr,&
             qx_gr, qy_gr, Q_bm_gr, Q_tld_gr)
        dataread = .TRUE.
     ELSE IF ((noption == 2).and.(dataread)) THEN
        WRITE( *, '(A,A)', ADVANCE = 'YES') 'Writing timestep-grid (ElmerPost) for run: ', runname
        CALL postgrid(xi_gr, eta_gr, z_c_gr, z_t_gr, Dx, imax, jmax, kcmax, ktmax,&
             runname, ergnum, maske_gr, 1)
     ELSE IF ((noption == 3).and.(dataread)) THEN
        WRITE( *, '(A,A)', ADVANCE = 'YES') 'Writing timestep-data (ElmerPost) for run: ', runname
        CALL elmerdata(imax, jmax, kcmax, ktmax, z_c_gr, z_t_gr,&
             vx_c_gr, vy_c_gr, vz_c_gr, age_c_gr, temp_c_gr, vx_t_gr, vy_t_gr,&
             vz_t_gr, temp_t_gr,  age_t_gr, omega_t_gr, Q_bm_gr,  Q_tld_gr,&
             am_perp_gr, qx_gr, qy_gr, n_cts_gr, maske_gr, runname, ergnum, flag)
     ELSE IF ((noption == 4).and.(dataread)) THEN
        WRITE( *, '(A,A)', ADVANCE = 'YES') 'Writing timestep-grid (ElmerSolver) for run: ', runname
        call pregrid(xi_gr, eta_gr, z_c_gr, z_t_gr, imax, jmax, kcmax, ktmax, runname, ergnum,&
             maske_gr, DX, flag)
     ELSE IF ((noption == 5).and.(dataread)) THEN
        WRITE( *, '(A,A)', ADVANCE = 'YES') 'Output of timestep-data in ASCII-format for run: ',runname 
        CALL  asciidata(xi_gr,eta_gr,IMAX,JMAX,KCMAX,KTMAX,&
             z_c_gr, z_t_gr, vx_c_gr,vy_c_gr,vz_c_gr,age_c_gr,temph_c_gr,&
             vx_t_gr,vy_t_gr,vz_t_gr,temph_t_gr,age_t_gr, omega_t_gr,&
             Q_bm_gr, Q_tld_gr,am_perp_gr,qx_gr,qy_gr,n_cts_gr,&
             maske_gr,runname,ergnum,1)
     ELSE IF (noption == 0) THEN
        DO i=1,5
           WRITE( *, '(A)', ADVANCE = 'YES') ' '
        END DO
        WRITE( *, '(A)', ADVANCE = 'YES') 'Thank you for using sico2elmer!'
        WRITE( *, '(A)', ADVANCE = 'YES') 'Bye'
        STOP
     ELSE
        WRITE( *, '(A)', ADVANCE = 'NO') 'Could not perform command due to '
        IF (noption > 5) THEN
           WRITE( *, '(A)', ADVANCE = 'YES') 'invalid selection!'
        ELSE
           WRITE( *, '(A)', ADVANCE = 'YES') 'missing timeslice data! Select option (1) first'
        END IF
     END IF
     ! clear screen
!     PAUSE
     DO i=1,5
        WRITE( *, '(A)', ADVANCE = 'YES') ' '
     END DO
  END DO

CONTAINS
!==============================================================================
  SUBROUTINE ReadData(runname, ergnum, &
             imax,jmax,kcmax,ktmax,krmax,deform,& 
             maske_gr, n_cts_gr,&
             time_gr, xi_gr, eta_gr, zs_gr, zm_gr,zb_gr, zb0_gr,&
             z_c_gr, z_t_gr, z_r_gr,&
             H_c_gr, H_t_gr, H_R_gr, vx_c_gr, vy_c_gr, vz_c_gr,&
             vx_t_gr, vy_t_gr, vz_t_gr, temp_c_gr, temp_r_gr, temp_t_gr,&
             temph_c_gr, temph_t_gr, am_perp_gr, age_c_gr, age_t_gr, omega_t_gr,&
             qx_gr, qy_gr, Q_bm_gr, Q_tld_gr)
    IMPLICIT NONE
!   external variables:
!   ------------------------------------------------------------------------
    CHARACTER :: &
         runname * 5, ergnum * 2
    INTEGER ::&
         imax,jmax,kcmax,ktmax,krmax 
    REAL :: &
         time_gr, deform, H_R_gr
    INTEGER  ::&
       maske_gr(0:imax,0:jmax), n_cts_gr(0:imax,0:jmax)
    REAL  ::&
         xi_gr(0:imax), eta_gr(0:jmax),&
         zs_gr(0:imax,0:jmax), zm_gr(0:imax,0:jmax),zb_gr(0:imax,0:jmax), zb0_gr(0:imax,0:jmax),&
         z_c_gr(0:imax,0:jmax,0:kcmax), z_t_gr(0:imax,0:jmax,0:ktmax), z_r_gr(0:imax,0:jmax,0:krmax),&
         H_c_gr(0:imax,0:jmax), H_t_gr(0:imax,0:jmax),&
         vx_c_gr(0:imax,0:jmax,0:kcmax), vy_c_gr(0:imax,0:jmax,0:kcmax), vz_c_gr(0:imax,0:jmax,0:kcmax),&
         vx_t_gr(0:imax,0:jmax,0:ktmax), vy_t_gr(0:imax,0:jmax,0:ktmax), vz_t_gr(0:imax,0:jmax,0:ktmax),&
         temp_c_gr(0:imax,0:jmax,0:kcmax), temp_r_gr(0:imax,0:jmax,0:krmax), temp_t_gr(0:imax,0:jmax,0:ktmax),&
         temph_c_gr(0:imax,0:jmax,0:kcmax), temph_t_gr(0:imax,0:jmax,0:ktmax),&
         age_c_gr(0:imax,0:jmax,0:kcmax), age_t_gr(0:imax,0:jmax,0:ktmax), omega_t_gr(0:imax,0:jmax,0:ktmax),&
         am_perp_gr(0:imax,0:jmax),&
         qx_gr(0:imax,0:jmax), qy_gr(0:imax,0:jmax), Q_bm_gr(0:imax,0:jmax), Q_tld_gr(0:imax,0:jmax)
!   internal variables:
!   ------------------------------------------------------------------------
    CHARACTER :: &
         ergfile * 11
    INTEGER :: &
         ios, kmax=0, i, j, k
    REAL :: &
         ScalarDummy
     REAL, PARAMETER :: &
          YEAR_SEC = 3.1556926e07, BETA = 8.70e-04
    INTEGER, ALLOCATABLE :: &
         TwoDimDummyI(:,:)
    REAL, ALLOCATABLE :: &
         OneDimDummyX(:), OneDimDummyY(:), TwoDimDummy(:,:),&
         ThreeDimDummyR(:,:,:), ThreeDimDummyT(:,:,:), ThreeDimDummyC(:,:,:)
    LOGICAL :: &
         FirstTime = .TRUE.

    SAVE ScalarDummy, TwoDimDummy, TwoDimDummyI, ThreeDimDummyR, ThreeDimDummyT,& 
         ThreeDimDummyC, OneDimDummyX, OneDimDummyY, FirstTime

    ! allocate some stuff (first time only)
    !---------------------------------------------------------------
    IF (FirstTime) THEN
       ALLOCATE(&
            OneDimDummyX(0:imax),&
            OneDimDummyY(0:jmax),&
            TwoDimDummy(0:jmax,0:imax),&
            TwoDimDummyI(0:jmax,0:imax),&
            ThreeDimDummyR(0:krmax,0:jmax,0:imax),&
            ThreeDimDummyT(0:ktmax,0:jmax,0:imax),&
            ThreeDimDummyC(0:kcmax,0:jmax,0:imax),&
            STAT=istat )
       
       IF ( istat /= 0 ) THEN
          WRITE( *, '(A)', ADVANCE = 'YES') 'Error in allocation of memory!'
          STOP
       ELSE
          WRITE( *, '(A)', ADVANCE = 'YES') 'Allocation of dummy arrays done'
       END IF 
       FirstTime = .FALSE.
    END IF

    ! inquire timeslice and open file
    !--------------------------------
    WRITE( *, '(A)', ADVANCE = 'NO') 'Enter number of timeslice (may start with 0 if < 10): '
    READ (*,'(A)', ADVANCE = 'YES') ergnum
    
    ergfile = runname//ergnum//'.erg'
    
    WRITE( *, '(A,A)', ADVANCE = 'YES') 'Atempting to open timeslice-file ', ergfile

    OPEN(UNIT=10, IOSTAT=ios, FILE=ergfile, STATUS='old', FORM='unformatted')
    IF (ios /= 0) THEN
       WRITE(6,'(A,A,A)') 'Error occurred while opening timeslice-file ', ergfile,'!'
       STOP
    ELSE
       WRITE( *, '(A,A)', ADVANCE = 'YES') 'Reading from ', ergfile
    END IF

    ! read in stuff
    !-------------
    !time
    READ(10) ScalarDummy
    time_gr = ScalarDummy/YEAR_SEC ! sec -> a
    !x-coords
    READ(10) OneDimDummyX
    xi_gr(0:imax) = 1.0e-03 * OneDimDummyX(0:imax) ! m -> km
    !y-coords
    READ(10) OneDimDummyY
    eta_gr(0:jmax) = 1.0e-03 * OneDimDummyY(0:jmax) ! m -> km
    ! glaciation mask
    read(10) TwoDimDummyI
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertFieldI2(imax, jmax, TwoDimDummyI, maske_gr) ! 1 -> 1
    ! poly-thermal condition mask
    read(10) TwoDimDummyI
    WRITE( *, '(A)', ADVANCE = 'NO') '.' 
    CALL ConvertFieldI2(imax, jmax, TwoDimDummyI, n_cts_gr) ! 1 -> 1
    ! z-coords 
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, zs_gr, 1.0e-03) !  m -> km
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, zm_gr, 1.0e-03) !  m -> km
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, zb_gr, 1.0e-03) !  m -> km
    ! depths
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, H_c_gr, 1.0e-03) !  m -> km
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, H_t_gr, 1.0e-03) !  m -> km    
    read(10) ScalarDummy
    H_R_gr = 1.0e-03 * ScalarDummy! m --> km
    ! convert heights
!    CALL ElevationColdD(imax, jmax, kcmax, deform, xi_gr, eta_gr,  z_c_gr)
    CALL ElevationCold(imax, jmax, kcmax, deform, zm_gr, H_c_gr,  z_c_gr)
    CALL ElevationTemp(imax, jmax, ktmax, zb_gr, H_t_gr, z_t_gr)
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    ! velocity components cold region
    read(10) ThreeDimDummyC
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, kcmax, ThreeDimDummyC, vx_c_gr, YEAR_SEC) ! (m/s) --> (m/yr)
    read(10) ThreeDimDummyC
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, kcmax, ThreeDimDummyC, vy_c_gr, YEAR_SEC) ! (m/s) --> (m/yr)
    read(10) ThreeDimDummyC
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, kcmax, ThreeDimDummyC, vz_c_gr, YEAR_SEC) ! (m/s) --> (m/yr)    
    ! velocity components temperate region
    read(10) ThreeDimDummyT
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, ktmax, ThreeDimDummyT, vx_t_gr, YEAR_SEC) ! (m/s) --> (m/yr)
    read(10) ThreeDimDummyT
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, ktmax, ThreeDimDummyT, vy_t_gr, YEAR_SEC) ! (m/s) --> (m/yr)
    read(10) ThreeDimDummyT
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, ktmax, ThreeDimDummyT, vz_t_gr, YEAR_SEC) ! (m/s) --> (m/yr)
    ! temperature field cold region (Celsius and homologe temperature)
    read(10) ThreeDimDummyC
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, kcmax, ThreeDimDummyC, temp_c_gr, 1.0e0) ! in C
    CALL HomologousTempC(imax, jmax, kcmax, temp_c_gr, H_c_gr,  deform, BETA, temph_c_gr)  ! with respect to pressure melting point
    ! water content in temperate layer
    read(10) ThreeDimDummyT
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, ktmax, ThreeDimDummyT, omega_t_gr, 1.0e0) ! dimensionless
    ! temperature field bedrock 
    read(10) ThreeDimDummyR(0:krmax,0:jmax,0:imax)
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, krmax, ThreeDimDummyR, temp_r_gr,  1.0e0) ! in C
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, Q_bm_gr, YEAR_SEC) ! m3/(m2*s) --> m3/(m2*yr)
!!$    DO i = 0, imax 
!!$       WRITE( *,'(I):', ADVANCE = 'NO') i
!!$       DO j = 0, jmax
!!$          WRITE( *,'(F8.4) ', ADVANCE = 'NO') Q_bm_gr(j,i)
!!$       END DO
!!$       WRITE( *,'(A)', ADVANCE = 'YES') 'END'
!!$    END DO
    ! water drainage rate from the temperated region
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, Q_tld_gr, YEAR_SEC) ! m3/(m2*s) --> m3/(m2*yr)
    ! ice volume flux through CTS 
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, am_perp_gr, YEAR_SEC) ! (m/s) --> (m/yr)
    ! volume-flux components
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, qx_gr, 1.0e-03) ! m2/s --> 1000 m2/yr
    read(10) TwoDimDummy
    WRITE( *, '(A)', ADVANCE = 'NO') '.'
    CALL ConvertField2(imax, jmax, TwoDimDummy, qy_gr, 1.0e-03) ! m2/s --> 1000 m2/yr
    read(10) ThreeDimDummyC
    WRITE( *, '(A)', ADVANCE = 'NO') '....'
    CALL ConvertField3(imax, jmax, kcmax, ThreeDimDummyC, age_c_gr, 1.0e0/(1.0e03*YEAR_SEC)) ! s --> kyr
    read(10) ThreeDimDummyT
    WRITE( *, '(A)', ADVANCE = 'YES') '....'
    CALL ConvertField3(imax, jmax, ktmax, ThreeDimDummyT, age_t_gr, 1.0e0/(1.0e03*YEAR_SEC)) ! s --> kyr
    WRITE( *, '(A)', ADVANCE = 'YES') 'Read in completed'
    ! close file
    !-----------
    CLOSE(10, STATUS='keep')
!==============================================================================
  END SUBROUTINE ReadData
  !--------------------------------------------------------------------------
  !   Writes 2d arrays read in from file into graphical output arrays
  !--------------------------------------------------------------------------
  SUBROUTINE ConvertField2(imax, jmax, in, out, factor)
    IMPLICIT NONE
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax
    REAL ::&
        in(0:jmax,0:imax), out(0:imax,0:jmax), factor !in(0:(imax+1)*(jmax+1)-1)
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j
    
    DO i = 0, imax
       DO j = 0, jmax
!          out(i,j) = factor * in(i*(jmax+1) + j)
          out(i,j) = factor * in(j,i)
       END DO
    END DO
  END SUBROUTINE ConvertField2
!==============================================================================
  !--------------------------------------------------------------------------
  !   Writes 2d integer arrays read in from file into graphical output arrays
  !--------------------------------------------------------------------------
  SUBROUTINE ConvertFieldI2(imax, jmax, in, out)
    IMPLICIT NONE
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax
    INTEGER ::&
         in(0:jmax,0:imax), out(0:imax,0:jmax)
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j
    
    DO i = 0, imax
       DO j = 0, jmax
          out(i,j) = in(j,i)
       END DO
    END DO
  END SUBROUTINE ConvertFieldI2
!==============================================================================
  !--------------------------------------------------------------------------
  !   Writes 3d arrays read in from file into graphical output arrays
  !--------------------------------------------------------------------------
  SUBROUTINE ConvertField3(imax, jmax, kmax, in, out, factor)
    IMPLICIT NONE
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax, kmax
    REAL ::&
         in(0:kmax,0:jmax,0:imax), out(0:imax,0:jmax,0:kmax), factor
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j, k
    
    DO i = 0, imax
       DO j = 0, jmax
          DO k=0, kmax
             out(i,j,k) =  factor * in(k,j,i)
          END DO
       END DO
    END DO
  END SUBROUTINE ConvertField3
!==============================================================================
  !--------------------------------------------------------------------------
  !   Transforms Celsius temperature into pressure meltinghomologous temperature
  !--------------------------------------------------------------------------
  SUBROUTINE HomologousTempC(imax, jmax, kcmax, in, depth, deform, BETA, out)
    IMPLICIT NONE
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax, kcmax
    REAL ::&
         in(0:imax,0:jmax,0:kcmax), out(0:imax,0:jmax,0:kcmax),&
         depth(0:imax,0:jmax)
    REAL :: &
         BETA, deform
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j, k
    REAL :: &
         ea, eaz_c_quotient, zeta_c, eaz_c

    ea = exp(deform)
    DO k=0, kcmax
       zeta_c = real(k)/real(kcmax)
       eaz_c = exp(DEFORM*zeta_c)
       eaz_c_quotient =(eaz_c-1.0)/(ea-1.0)
       DO i = 0, imax
          DO j = 0, jmax
             out(i,j,k) =  in(i,j,k)&
                  - ( -1000.0*BETA*depth(i,j)*(1.0-eaz_c_quotient) )
          END DO
       END DO
    END DO
  END SUBROUTINE HomologousTempC

!==============================================================================
  !--------------------------------------------------------------------------
  !   Transforms Celsius temperature into pressure meltinghomologous temperature
  !--------------------------------------------------------------------------
  SUBROUTINE HomologousTempT(imax, jmax, ktmax, in, depth1, depth2, BETA, out)
    IMPLICIT NONE
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax, ktmax
    REAL ::&
         in(0:imax,0:jmax,0:ktmax), out(0:imax,0:jmax,0:ktmax),&
         depth1(0:imax,0:jmax), depth2(0:imax,0:jmax)
    REAL :: &
         BETA
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j, k         
    REAL :: &
         zeta_t 

    DO k=0, ktmax
       zeta_t  = real(k)/real(ktmax)
       DO i = 0, imax
          DO j = 0, jmax
             out(i,j,k) =  in(i,j,k)&
                  - ( -1000.0*BETA*(depth1(i,j)+depth2(i,j))*(1.0-zeta_t) )
          END DO
       END DO
    END DO
  END SUBROUTINE HomologousTempT
!==============================================================================
  !--------------------------------------------------------------------------
  !   Transforms normalized into real elevations in cold layer
  !--------------------------------------------------------------------------
  SUBROUTINE ElevationTemp(imax, jmax, ktmax, zb_gr, H_t_gr, z_t_gr)
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax, ktmax
    REAL ::&
         zm_gr(0:imax,0:jmax), zb_gr(0:imax,0:jmax),&
         H_t_gr(0:imax,0:jmax), z_t_gr(0:imax,0:jmax,0:ktmax)
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j, k
    REAL :: &
         zeta_t
    
    DO k=0, ktmax
       zeta_t  = real(k)/real(ktmax)
       DO i = 0, imax
          DO j = 0, jmax
              z_t_gr(i,j,k) = zb_gr(i,j) +H_t_gr(i,j)*zeta_t
          END DO
       END DO
    END DO
  END SUBROUTINE ElevationTemp
!==============================================================================
  !--------------------------------------------------------------------------
  !   Transforms normalized into real elevations in cold layer
  !--------------------------------------------------------------------------
  SUBROUTINE ElevationColdD(imax, jmax, kcmax, deform, xi_gr, eta_gr,  z_c_gr)
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax, kcmax
    REAL ::&
         xi_gr(0:imax), eta_gr(0:jmax), z_c_gr(0:imax,0:jmax,0:kcmax)
    REAL :: &
         deform
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j, k
    DO i = 0, imax
       DO j = 0, jmax
          DO k=0, kcmax
!             z_c_gr(i,j,k) =  xi_gr(i) * eta_gr(j)
             z_c_gr(i,j,k) = 0.0e0
          END DO
       END DO
    END DO
  END SUBROUTINE ElevationColdD

  SUBROUTINE ElevationCold(imax, jmax, kcmax, deform, zm_gr, H_c_gr,  z_c_gr)
    ! external variables
    !-------------------
    INTEGER ::&
         imax, jmax, kcmax
    REAL ::&
         zm_gr(0:imax,0:jmax), H_c_gr(0:imax,0:jmax), z_c_gr(0:imax,0:jmax,0:kcmax)
    REAL :: &
         deform
    ! internal variables
    !-------------------
    INTEGER ::&
         i, j, k
    REAL :: &
         ea, eaz_c_quotient, zeta_c, eaz_c
    
    ea = exp(deform)
    DO k=0, kcmax
       zeta_c = real(k)/real(kcmax)
       eaz_c = exp(DEFORM*zeta_c)
       eaz_c_quotient =(eaz_c-1.0e0)/(ea-1.0e0)
       DO i = 0, imax
          DO j = 0, jmax
             z_c_gr(i,j,k) = zm_gr(i,j) + H_c_gr(i,j)*eaz_c_quotient
          END DO
       END DO
    END DO
  END SUBROUTINE ElevationCold
!==============================================================================
  SUBROUTINE  ReadLog(runname, imax,jmax,kcmax,ktmax,krmax,deform,Dx)
    IMPLICIT NONE
!   external variables:
!   ------------------------------------------------------------------------

    INTEGER ::&
         imax,jmax,kcmax,ktmax,krmax 
    CHARACTER :: &
       runname * 5
    REAL ::&
         deform, Dx, rDummy
!   internal variables:
!    ------------------------------------------------------------------------
    CHARACTER :: &
         logdat * 9, chtrash
    INTEGER :: &
         ios
     
    logdat = runname//'.log'
    WRITE( *, '(A,A)', ADVANCE = 'YES') 'Atempting to open log-file ', logdat
    OPEN(UNIT=10, iostat=ios, file=logdat, status='old')
    IF (ios.ne.0) THEN
       WRITE(6,'(A)') 'Error occurred while opening log-file!'
       STOP
    ELSE
       WRITE( *, '(A,A)', ADVANCE = 'YES') 'Reading from ', logdat
    END IF    
!    READ(10,'(a7,i4)') chtrash
!    READ(10,'(a7,i4)') chtrash
!    READ(10,'(a7,i4)') chtrash
    READ(10,'(a7,i4)') chtrash, imax
    READ(10,'(a7,i4)') chtrash, jmax
    READ(10,'(a7,i4)') chtrash, kcmax
    READ(10,'(a7,i4)') chtrash, ktmax
    READ(10,'(a7,i4)') chtrash, krmax
    READ(10,'(a)') chtrash
    READ(10,'(a3,e9.2)')  chtrash, deform  
    READ(10,'(a,e9.2)') chtrash
    READ(10,'(a12,e9.2)') chtrash, rDummy
    READ(10,'(a12,e9.2)') chtrash, rDummy
    READ(10,'(a)') chtrash
    READ(10,'(a9,e9.2)') chtrash, Dx
    CLOSE(10)
    WRITE( *, '(A)', ADVANCE = 'YES')
    WRITE( *, '(A)', ADVANCE = 'YES') 'The following parameters have been read in:'
    WRITE( *, '(A,i4)', ADVANCE = 'YES') '   imax=', imax
    WRITE( *, '(A,i4)', ADVANCE = 'YES') '   jmax=', jmax
    WRITE( *, '(A,i4)', ADVANCE = 'YES') '  kcmax=', kcmax
    WRITE( *, '(A,i4)', ADVANCE = 'YES') '  ktmax=', ktmax
    WRITE( *, '(A,i4)', ADVANCE = 'YES') '  krmax=', krmax
    WRITE( *, '(A,f9.2)', ADVANCE = 'YES') ' deform=', deform
    WRITE( *, '(A,f9.2)', ADVANCE = 'YES') '     Dx=', Dx
    WRITE( *, '(A)', ADVANCE = 'YES') ' '
    WRITE( *, '(A)', ADVANCE = 'YES')
  END SUBROUTINE ReadLog
!==============================================================================
END PROGRAM sico2elmer
!==============================================================================
!==============================================================================
