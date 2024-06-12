*  Routine:    PIVOUT - Parallel version of ARPACK UTILITY ROUTINE IVOUT
*
*  Purpose:    Integer vector output routine.
*
*  Usage:      CALL PIVOUT (COMM, LOUT, N, IX, IDIGIT, IFMT)
*
*  Arguments
*     COMM   - MPI Communicator for the processor grid
*     N      - Length of array IX. (Input)
*     IX     - Integer array to be printed. (Input)
*     IFMT   - Format to be used in printing array IX. (Input)
*     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*\SCCS Information: 
* FILE: ivout.F   SID: 1.2   DATE OF SID: 3/19/97   RELEASE: 1
*
*-----------------------------------------------------------------------
*
      SUBROUTINE PIVOUT (COMM, LOUT, N, IX, IDIGIT, IFMT)
*
      include  'mpif.h'
*
*     .. MPI VARIABLES AND FUNCTIONS ..
*     .. Variable Declaration ..
      integer    COMM, MYID, IERR
*
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IX(*), N, IDIGIT, LOUT
      CHARACTER  IFMT*(*)
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, NDIGIT, K1, K2, LLL
      CHARACTER*80 LINE
*     ...
*     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
*
*     ..
*     .. Executable Statements ..
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
*     Determine processor configuration
*
      call MPI_COMM_RANK( comm, myid, ierr )
*
*     .. Only Processor 0 will write to file LOUT ..
*
      IF ( MYID .EQ. 0 ) THEN
*
      LLL = MIN ( LEN ( IFMT ), 80 )
      DO 1 I = 1, LLL
          LINE(I:I) = '-'
    1 CONTINUE
*
      DO 2 I = LLL+1, 80
          LINE(I:I) = ' '
    2 CONTINUE
*
      WRITE ( LOUT, 2000 ) IFMT, LINE(1:LLL)
 2000 FORMAT ( /1X, A  /1X, A )
*
      IF (N .LE. 0) RETURN
      NDIGIT = IDIGIT
      IF (IDIGIT .EQ. 0) NDIGIT = 4
*
*=======================================================================
*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
*=======================================================================
*
      IF (IDIGIT .LT. 0) THEN
*
      NDIGIT = -IDIGIT
      IF (NDIGIT .LE. 4) THEN
         DO 10 K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   10    CONTINUE
*
      ELSE IF (NDIGIT .LE. 6) THEN
         DO 30 K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
   30    CONTINUE
*
      ELSE IF (NDIGIT .LE. 10) THEN
         DO 50 K1 = 1, N, 5
            K2 = MIN0(N,K1+4)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
   50    CONTINUE
*
      ELSE
         DO 70 K1 = 1, N, 3
            K2 = MIN0(N,K1+2)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
   70    CONTINUE
      END IF
*
*=======================================================================
*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
*=======================================================================
*
      ELSE
*
      IF (NDIGIT .LE. 4) THEN
         DO 90 K1 = 1, N, 20
            K2 = MIN0(N,K1+19)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   90    CONTINUE
*
      ELSE IF (NDIGIT .LE. 6) THEN
         DO 110 K1 = 1, N, 15
            K2 = MIN0(N,K1+14)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
  110    CONTINUE
*
      ELSE IF (NDIGIT .LE. 10) THEN
         DO 130 K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
  130    CONTINUE
*
      ELSE
         DO 150 K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
  150    CONTINUE
      END IF
      END IF
      WRITE (LOUT,1004)
 
      ENDIF
*
 1000 FORMAT(1X,I4,' - ',I4,':',20(1X,I5))
 1001 FORMAT(1X,I4,' - ',I4,':',15(1X,I7))
 1002 FORMAT(1X,I4,' - ',I4,':',10(1X,I11))
 1003 FORMAT(1X,I4,' - ',I4,':',7(1X,I15))
 1004 FORMAT(1X,' ')
*
      RETURN
      END
