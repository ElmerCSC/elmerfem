c\BeginDoc
c
c\Name: pslarnv
c
c Message Passing Layer: MPI
c
c\Description:
c
c  Parallel Version of ARPACK utility routine slarnv
c
c  PSLARNV returns a vector of n (nloc) random real numbers from a uniform or
c  normal distribution. It is assumed that X is distributed across a 1-D array 
c  of processors ( nprocs < 1000 )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid
c
c  IDIST   (input) INTEGER
c          Specifies the distribution of the random numbers:
c          = 1:  uniform (0,1)
c          = 2:  uniform (-1,1)
c          = 3:  normal (0,1)
c
c  ISEED   (input/output) INTEGER array, dimension (4)
c          On entry, the seed of the random number generator; the array
c          elements must be between 0 and 4095, and ISEED(4) must be
c          odd.
c          On exit, the seed is updated.
c
c  N       (input) INTEGER
c          The number of random numbers to be generated.
c
c  X       (output) Real array, dimension (N)
c          The generated random numbers.
c
c\Author: Kristi Maschhoff
c
c\Details
c
c  Simple parallel version of LAPACK auxiliary routine slarnv 
c  for X distributed across a 1-D array of processors.
c  This routine calls the auxiliary routine SLARNV to generate random
c  real numbers from a uniform (0,1) distribution. Output is consistent
c  with serial version. 
c
c\SCCS Information: 
c FILE: larnv.F   SID: 1.4   DATE OF SID: 04/16/99   
c
c-----------------------------------------------------------------------
c
      subroutine pslarnv( comm, idist, iseed, n, x )
c
      include  'mpif.h'
c
c     .. MPI VARIABLES AND FUNCTIONS ..
      integer   comm
c     ..
c     .. Scalar Arguments ..
      integer			idist, n
c     ..
c     .. Array Arguments ..
      integer			iseed( 4 )
      Real			
     &                  x( * )
c     ..
c     .. External Subroutines ..
      external			slarnv
c     ..
c     .. Executable Statements ..
c
      call slarnv ( idist, iseed, n, x )
c
      return
      end
