!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!!!
!!!  modules to interface the Scattered2DDataInterpolator Solver with the 
!!!    external c libraries nn and csa (http://code.google.com/p/nn-c http://code.google.com/p/csa-c)
        module point_type
           USE ISO_C_BINDING
            TYPE,BIND(C) :: POINT
              REAL(C_DOUBLE) :: x,y,z
            END TYPE POINT
        end module point_type

        module NearestNeighbour
        USE ISO_C_BINDING
        USE point_type
  
       
        enum , BIND(C) 
            ENUMERATOR :: SIBSON, NON_SIBSONIAN
        end enum
        common /nn_rule/nn_rule
        INTEGER(KIND(SIBSON)) :: nn_rule
        BIND(C) :: /nn_rule/

        INTERFACE
           SUBROUTINE nnpi_interpolate_points(nin,pin,wmin,nout,pout) BIND(C)
           USE ISO_C_BINDING
           USE point_type
           implicit none
           INTEGER(C_INT), VALUE :: nin,nout
           REAL(C_DOUBLE),VALUE :: wmin
           TYPE(POINT) :: pin(nin),pout(nout)
          END SUBROUTINE nnpi_interpolate_points
          
          SUBROUTINE csa_interpolate_points(nin,pin,nout,pout,nppc, k) BIND(C)
           USE ISO_C_BINDING
           USE point_type
           implicit none
           INTEGER(C_INT), VALUE :: nin,nout,nppc,k
           TYPE(POINT) :: pin(nin),pout(nout)
          END SUBROUTINE csa_interpolate_points

          SUBROUTINE lpi_interpolate_points(nin,pin,nout,pout) BIND(C)
           USE ISO_C_BINDING
           USE point_type
           implicit none
           INTEGER(C_INT), VALUE :: nin,nout
           TYPE(POINT) :: pin(nin),pout(nout)
          END SUBROUTINE lpi_interpolate_points
        END INTERFACE
        end module NearestNeighbour
