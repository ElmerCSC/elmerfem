c\BeginDoc
c
c\Name: psnorm2
c
c Message Passing Layer: MPI
c
c\Description:
c
c\Usage:
c  call psnorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
c
c-----------------------------------------------------------------------
c
      Real function psnorm2 ( comm, n, x, inc )
c
      include   'mpif.h'
c
c     %---------------%
c     | MPI Variables |
c     %---------------%
c
      integer    comm, ierr
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer      n, inc
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Real
     &             x(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      Real
     &             max, buf, zero
      parameter    ( zero = 0.0 )
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Real       
     &             snrm2
      External     snrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      psnorm2 = snrm2( n, x, inc)
c
      buf = psnorm2
      call MPI_ALLREDUCE( buf, max, 1, MPI_REAL,
     &                    MPI_MAX, comm, ierr )
      if ( max .eq. zero ) then
         psnorm2 = zero
      else
         buf = (psnorm2/max)**2.0
         call MPI_ALLREDUCE( buf, psnorm2, 1, MPI_REAL,
     &                       MPI_SUM, comm, ierr )
         psnorm2 = max * sqrt(abs(psnorm2))
      endif
c
c     %----------------%
c     | End of psnorm2 |
c     %----------------%
c
      return
      end
