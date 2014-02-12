/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/*******************************************************************************
 *
 * Type & structure definitions for objects & geometry. This is really the
 * definition of the structure of ElmerPost.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 26 Sep 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

#ifdef MODULE_GEOMETRY
#   define GEO_EXT
#else
#   define GEO_EXT extern
#endif

#define GEO_TRIANGLE_BLOCK_SIZE 4096
#define GEO_VERTEX_BLOCK_SIZE   4096

#define FLOAT float

typedef struct group_s
{
    struct group_s *Next;
    int status,Open;
    char *Name;
} group_t;

/*
 * Triangle
 */
typedef struct
{
    int   v[3];    /* vertex pointers */
    FLOAT Fu[3];   /* triangle normal */
    FLOAT u[3][3]; /* vertex normals  */

    struct element_s *Element;

    int Count;
    logical_t Edge[3];  /* beginning of an edge flag */
} triangle_t;

/*
 * list of faces connected to a vertex
 */
typedef struct vertex_face_s
{
    int Face;
    struct vertex_face_s *Next;
} vertex_face_t;

/*
 *  vertex def's
 */
typedef struct vertex_s
{
    FLOAT x[3];
    vertex_face_t *Faces;
    logical_t ElementModelNode;
} vertex_t;

typedef struct
{
    FLOAT x[3],y[3],z[3];
    FLOAT u[3],v[3],w[3];
    FLOAT c[3],f[3];
} polygon_t;

typedef enum
{
    line_style_line,line_style_cylinder
} line_style_t;

typedef struct line_s
{
    float x[3];
    float y[3];
    float z[3];
    float f[3],c[3];
} line_t;

/*
 * Edges of elements are hold in an array of lists that are
 * indexed by smallest numbered vertex of a particular edge.
 */
typedef enum
{
    edge_style_all, edge_style_free
} edge_style_t;

typedef struct edge_list_s
{
    struct edge_list_s *Next;
    int Entry,Count;
    struct element_s *Element;
} edge_list_t;
    
typedef struct edge_s
{
     edge_list_t *EdgeList; 
} edge_t;

/*
 *  geometry def's
 */
typedef struct geometry_s
{
    struct Geometry_s *Next;
    char *Name;

    triangle_t *Triangles;
    int TriangleCount,MaxTriangleCount;

    vertex_t *Vertices;
    int VertexCount,MaxVertexCount;

    edge_t *Edges;

    double Scale;
    vertex_t MinMax[2];
} geometry_t;

GEO_EXT geometry_t Geometry;

typedef struct data_s
{
    int a;    
} data_t;


void geo_free_groups( group_t *groups );
int geo_add_vertex( geometry_t *geometry, vertex_t *vertex );
void geo_free_edge_tables( geometry_t *geometry );
void geo_free_vertex_face_tables( geometry_t *geometry );
int geo_add_triangle( geometry_t *geometry, triangle_t *triangle );
void geo_triangle_normal( geometry_t *geom,triangle_t *triangle );
void geo_vertex_normals( geometry_t *geometry, double Ang ) ;
