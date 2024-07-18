c\BeginDoc
c
c\Name: pdnorm2
c
c Message Passing Layer: MPI
c
c\Description:
c
c\Usage:
c  call pdnorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
c
c-----------------------------------------------------------------------
c
      Double precision function pdnorm2 ( comm, n, x, inc )
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
      Double precision
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
     &             dnrm2
      External     dnrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      pdnorm2 = dnrm2( n, x, inc)
c
      buf = pdnorm2
      call MPI_ALLREDUCE( buf, max, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, comm, ierr )
      if ( max .eq. zero ) then
         pdnorm2 = zero
      else
         buf = (pdnorm2/max)**2.0
         call MPI_ALLREDUCE( buf, pdnorm2, 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, comm, ierr )
         pdnorm2 = max * sqrt(abs(pdnorm2))
      endif
c
c     %----------------%
c     | End of pdnorm2 |
c     %----------------%
c
      return
      end
