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
! *  Authors: Juhani Kataja
! *  Email:   juhani.kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Aug 24, 2023
! *
! *****************************************************************************/


SUBROUTINE VolumeProbe( Model, Solver, dt, TransientSimulation )
  USE DefUtils

  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
! Local variables
  
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueHandle_t), POINTER :: probe_h(:)
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp), allocatable :: integral(:)
  integer :: active, cdim, n, nd, ngp, t
  TYPE(Element_t), POINTER :: elm
  LOGICAL :: Found, stat, parallel
  CHARACTER(LEN=MAX_NAME_LEN) :: probe_names
  integer, allocatable :: biases(:)

  Mesh => GetMesh()

  probe_names = ListGetString(Solver % values, 'Body force probe names', found)

  ! Chop 'probe names' string with white-space delimiters
  if (found) then
    block
      integer :: bias, numsubs, lenstr, counter, addbias


      lenstr = len(trim(probe_names))
      store: do while (.true.)

        bias = 1
        numsubs = 0
        counter = 0 ! to force exit at maxiter

        remove_left_ws: do while (probe_names(bias:bias) == " ")
          bias = bias + 1
        end do remove_left_ws
        if (bias >= lenstr) then ! empty string full of ws
          exit store
        end if
        numsubs = 1 ! at least 1 substring

        if (allocated(biases)) then
          biases(1) = bias
        end if

        scanstr: do while (.true.)
          remove_left_nws: do while (probe_names(bias:bias) /= " " .and. bias <= lenstr)
            bias = bias+1
          end do remove_left_nws

          if (bias >= lenstr) then ! no more nws
            numsubs = numsubs + 1
            if (allocated(biases)) then
              biases(numsubs) = lenstr+1 ! last one is length of the string
            end if  
            exit scanstr
          end if

          hog_whitespace: do while (probe_names(bias:bias) == " " .and. bias <= lenstr)
            bias = bias + 1
          end do hog_whitespace

          numsubs = numsubs + 1
          if (allocated(biases)) then
            biases(numsubs) = bias
          end if

          ! sanity quit
          counter = counter + 1
          if (counter .ge. 10) then
            exit scanstr
          end if
        end do scanstr

        if (.not. allocated(biases)) then
          allocate(biases(numsubs))
          allocate(probe_h(numsubs-1))
          allocate(integral(numsubs-1))
        else
          exit store
        end if
      end do store

      do counter = 1, numsubs-1
        integral(counter)=0.0

        call ListInitElementKeyword( probe_h(counter), 'Body Force', &
            trim(probe_names(biases(counter):biases(counter+1)-1)) )
      end do

    end block
  end if

  cdim = CoordinateSystemDimension()
  Active = GetNOFActive(Solver)

  do t=1,active
    elm => GetActiveElement(t)
    IP = GaussPoints(elm)
    ngp = IP % n

    n  = GetElementNOFNodes(elm)
    nd = GetElementNOFDOFs(elm) + GetElementNOFBDOFs(elm)

    CALL GetElementNodes( Nodes, UElement=elm )

   block  ! actually collect the integral here
     REAL(KIND=dp) :: Basis(nd), detj, integrand_value
     integer :: nd_n, m, quad_n
       do quad_n = 1,ngp
         stat = ElementInfo( elm, Nodes, IP % U(quad_n), IP % V(quad_n), IP % W(quad_n), detj, Basis)
         do m = 1,size(probe_h)
           integrand_value = ListGetElementReal( probe_h(m), Basis, elm, Found, gausspoint=ngp)
           IF( Found ) THEN
             integral(m) = integral(m) + integrand_value*detj*IP%s(quad_n)
           END IF
         end do
       end do
   end block

  end do
  
 Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Solver % Mesh % SingleMesh ) 

 if (Parallel) then
   do t = 1,size(integral)
     integral(t) = parallelreduction(integral(t))
   end do
 end if

 do t = 1,size(integral)
   CALL ListAddConstReal( Model % Simulation, 'res: volume probe ' // & 
       trim(probe_names(biases(t):biases(t+1)-1)), integral(t) )
 end do

END SUBROUTINE VolumeProbe
