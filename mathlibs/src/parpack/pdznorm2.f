c\BeginDoc
c
c\Name: pdznorm2
c
c Message Passing Layer: MPI
c
c\Description:
c
c\Usage:
c  call pdznorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
c
c-----------------------------------------------------------------------
c
      Double precision function pdznorm2 ( comm, n, x, inc )
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
      Complex*16
     &             x(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      Double precision
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
      Double precision       
     &             dznrm2
      External     dznrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      pdznorm2 = dznrm2( n, x, inc)
c
      buf = pdznorm2
      call MPI_ALLREDUCE( buf, max, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, comm, ierr )
      if ( max .eq. zero ) then
         pdznorm2 = zero
      else
         buf = (pdznorm2/max)**2.0
         call MPI_ALLREDUCE( buf, pdznorm2, 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, comm, ierr )
         pdznorm2 = max * sqrt(abs(pdznorm2))
      endif
c
c     %-----------------%
c     | End of pdznorm2 |
c     %-----------------%
c
      return
      end
