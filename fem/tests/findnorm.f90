PROGRAM findnorm

 IMPLICIT NONE

 CHARACTER(LEN=512) :: line, cpustring, filename
 INTEGER :: Noargs, i,j,n, success
 REAL(8) :: target_nrm, target_eps = 1d-5, f, norm, CPUT, REALT, cumul_time

 NoArgs = IARGC() 
 CALL GetArg(1,filename)
 IF (NoArgs >= 2) THEN
   CALL GetArg(2,cpustring)
   READ(cpustring, *) cumul_time
 END IF

 OPEN( 1, File = filename, STATUS='old' )
 DO WHILE(.TRUE.)
   READ(1,'(a)', END=100, ERR=100 ) line
   i = INDEX( line, 'END TEST CASE' )
   IF (i>0 ) THEN
     i = INDEX( line, 'NRM=' )
     IF(i>0) READ(line(i+4:), * ) target_nrm

     i = INDEX( line, 'EPS=' )
     IF(i>0) READ(line(i+4:), * ) target_eps
   END IF

   i = INDEX(line, '(NRM,RELC)')
   IF(i>0) THEN
     DO j=i+10,LEN_TRIM(line)
       IF( line(j:j) == '(' ) THEN
         READ( line(j+1:), *, END=10,ERR=10 ) f
         IF(f/=0) norm = f
10       CONTINUE
       END IF
     END DO
   END IF

   i = INDEX(line, 'Check NRM')
   IF(i>0) THEN
     DO j=i+10,LEN_TRIM(line)
       IF( line(j:j) == '(' ) THEN
         READ( line(j+1:), *, END=20,ERR=20 ) f
         IF(f/=0) norm = f
20       CONTINUE
       END IF
     END DO
   END IF

   i = INDEX(line, '(CPU,REAL')
   IF (i>0) THEN
     READ(line(i+11:), *,END=30,ERR=30) REALT, CPUT
30   CONTINUE
   END IF
 END DO

100 CONTINUE

 Success = COMPARE( norm, target_nrm, target_eps );
 IF ( noargs < 2 ) THEN
   IF (Success<=0) THEN
     n = LEN_TRIM(line)
     j = n
     DO i=n,1,-1
       IF(ICHAR(line(i:i))==10.OR.ICHAR(line(i:i))==13) j=j-1
     END DO
     WRITE(6,'(a,g0.6,a,g0.6)' ) '[FAILED] Computed NRM=', norm,' Target NRM=', target_nrm
    ELSE
     WRITE(6,'(I1)') Success
    END IF
 ELSE
   WRITE(6, '(g10.4)') cumul_time + CPUT
 END IF


  FLUSH(6)

CONTAINS

FUNCTION Compare( norm1, norm2, eps ) RESULT(success)
  REAL(8) :: norm1, norm2, eps
  INTEGER :: success
 
  success = 0

  if ( norm1 /= -1 )then
    if ( eps < 0 ) then
       if ( norm2 < norm1 ) then
          success = 1
          return;
       else
          return;
       end if
    else  if ( 2 * abs(norm1-norm2) / (norm1+norm2) < eps ) then
       success = 1
       return
    else
       return
    end if
  else 
    success = 1
    return
  end if
END FUNCTION Compare

END PROGRAM findnorm
