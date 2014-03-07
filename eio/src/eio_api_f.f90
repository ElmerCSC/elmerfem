!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write
! * to the Free Software Foundation, Inc., 51 Franklin Street,
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Sami Ilvonen
! *  Email:   sami.ilvonen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 26 Feb 2014
! *
! *****************************************************************************/

! Module for wrapping selected parts of eio library using iso_c_binding 
! interface


MODULE EIOFortranAPI
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE

     SUBROUTINE eio_init(info) BIND(c, name="eio_init")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_init

     SUBROUTINE eio_init_parallel(procs, me, info) &
          & BIND(c, name="eio_init_parallel")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: procs, me, info
     END SUBROUTINE eio_init_parallel

     SUBROUTINE eio_close(info) BIND(c, name="eio_close")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_close

     SUBROUTINE eio_create_model(directory, info) &
          & BIND(c, name="eio_create_model")
       USE, INTRINSIC :: ISO_C_BINDING
       CHARACTER(kind=c_char) :: directory(*)
       INTEGER(c_int) :: info
     END SUBROUTINE eio_create_model

     SUBROUTINE eio_open_model(directory, info) BIND(c, name="eio_open_model")
       USE, INTRINSIC :: ISO_C_BINDING
       CHARACTER(kind=c_char) :: directory(*)
       INTEGER(c_int) :: info
     END SUBROUTINE eio_open_model

     SUBROUTINE eio_close_model(info) BIND(c, name="eio_close_model")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_close_model

     SUBROUTINE eio_create_geometry(info) BIND(c, name="eio_create_geometry")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_create_geometry

     SUBROUTINE eio_open_geometry(info) BIND(c, name="eio_open_geometry")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_open_geometry

     SUBROUTINE eio_close_geometry(info) BIND(c, name="eio_close_geometry")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_close_geometry

     SUBROUTINE eio_set_geometry_description(bodyC, boundaryC, outerC, &
          & innerC, vertexC, loopC, maxLooplen, info) &
          & BIND(c, name="eio_set_geometry_description")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: bodyC, boundaryC, outerC, innerC, vertexC, loopC,&
            & maxLooplen, info
     END SUBROUTINE eio_set_geometry_description

     SUBROUTINE eio_get_geometry_description(bodyC, boundaryC, outerC, &
          & innerC, vertexC, loopC, maxLooplen, info) &
          & BIND(c, name="eio_get_geometry_description")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: bodyC, boundaryC, outerC, innerC, vertexC, loopC, &
            & maxLooplen, info
     END SUBROUTINE eio_get_geometry_description

     SUBROUTINE eio_set_geometry_body(tag, meshControl, loopC, &
          & loops, info) BIND(c, name="eio_set_geometry_body")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, meshControl, loopC, loops, info
     END SUBROUTINE eio_set_geometry_body

     SUBROUTINE eio_get_geometry_body(tag, meshControl, loopC, &
          & loops, info) BIND(c, name="eio_get_geometry_body")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, meshControl, loopC, loops, info
     END SUBROUTINE eio_get_geometry_body

     SUBROUTINE eio_set_geometry_body_loop(tag, field, nodes, info) &
          & BIND(c, name="eio_set_geometry_body_loop")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, field, nodes, info
     END SUBROUTINE eio_set_geometry_body_loop

    SUBROUTINE eio_open_mesh(dir, info) &
        BIND(C, name="eio_open_mesh")
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(C_CHAR), INTENT(IN) :: dir(*)
        INTEGER(c_int) :: info
    END SUBROUTINE eio_open_mesh

     SUBROUTINE eio_close_mesh(info) BIND(c, name="eio_close_mesh")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: info
     END SUBROUTINE eio_close_mesh

     SUBROUTINE eio_set_mesh_bndry_element(tag, boundary, leftElement, &
          & rightElement, TYPEDEF, nodes, info) &
          & BIND(c, name="eio_set_mesh_bndry_element")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, boundary, leftElement, rightElement, &
            & TYPEDEF, nodes(*), info
     END SUBROUTINE eio_set_mesh_bndry_element

     SUBROUTINE eio_get_mesh_bndry_element(tag, part, boundary, leftElement, &
          & rightElement, typedef, nodes, coord, info) &
          & BIND(c, name="eio_get_mesh_bndry_element")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, part, boundary, leftElement, rightElement, &
            & typedef, nodes(*), info
       REAL(c_double) :: coord(*)
     END SUBROUTINE eio_get_mesh_bndry_element

     SUBROUTINE eio_get_mesh_description(nodeCount, elementCount, &
          & boundaryElementCount, usedElementTypes, elementTypeTags, &
          & elementCountByType, info) &
          & BIND(c, name="eio_get_mesh_description")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: nodeCount, elementCount, boundaryElementCount, &
            & usedElementTypes, elementTypeTags(*), elementCountByType(*), info
     END SUBROUTINE eio_get_mesh_description

    SUBROUTINE eio_get_mesh_element_conns(tag, part, body, typedef, &
                pdofs, nodes, info) &
        BIND(c, name="eio_get_mesh_element_conns")
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT) :: tag, part, body, typedef, pdofs(*), nodes(*), info
    END SUBROUTINE eio_get_mesh_element_conns

    SUBROUTINE eio_get_mesh_nodes(tags, coord, info) &
          & BIND(c, name="eio_get_mesh_nodes")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tags(*), info
       REAL(c_double) :: coord(*)
     END SUBROUTINE eio_get_mesh_nodes

     SUBROUTINE eio_get_part_description(sharedNodeCount, info) &
          & BIND(c, name="eio_get_part_description")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: sharedNodeCount, info
     END SUBROUTINE eio_get_part_description

     SUBROUTINE eio_set_part_node(tag, typedef, coord, partcount, &
          & parts, info) BIND(c, name="eio_set_part_node")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, typedef, partcount, parts(*), info
       REAL(c_double) :: coord(*)
     END SUBROUTINE eio_set_part_node

     SUBROUTINE eio_get_part_node(tag, constraint, coord, partcount, &
          & parts, info) BIND(c, name="eio_get_part_node")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(c_int) :: tag, constraint, partcount, parts(*), info
       REAL(c_double) :: coord(*)
     END SUBROUTINE eio_get_part_node

  END INTERFACE
END MODULE EIOFortranAPI
