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
 *     Geometry manipulation utilities.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address:   Keilaranta 14, P.O. BOX 405
 *                                 02101 Espoo, Finland
 *                                 Tel. +358 0 457 2723
 *                               Telefax: +358 0 457 2302
 *                             EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 26 Sep 1995
 *
 * Modification history:
 *
 * 27 Sep 1995: added geo_add_edge - routine, see below, Juha Ruokolainen
 *
 ******************************************************************************/

#define MODULE_GEOMETRY

#include <elmerpost.h>

#include <tcl.h>
extern Tcl_Interp *TCLInterp;

void geo_free_groups( group_t *groups )
{
    group_t *grp;
    int gid = 0;
    char str[64];

    for( ; groups != NULL; gid++ )
    {
        grp = groups->Next;

        sprintf( str, "GroupStatus(%d)", gid );
        Tcl_UnlinkVar( TCLInterp, str ); 

        free( groups );
        groups = grp;
    }
    Tcl_SetVar( TCLInterp, "NumberOfGroups", "0", TCL_GLOBAL_ONLY );
}

int geo_group_id( element_model_t *model, char *name,int open )
{
   group_t *group,*prevgroup;
   int groupid = 0;

   static char str[128];

   prevgroup  = model->Groups;
   for( group = model->Groups; group!=NULL; group=group->Next ) 
   {
      if ( strcmp(group->Name, name) == 0 ) break;
      groupid++;
      prevgroup = group;
   }

   if ( !group )
   {
      if ( !prevgroup )
         group = model->Groups = (group_t *)malloc(sizeof(group_t));
      else
         group = prevgroup->Next = (group_t *)malloc(sizeof(group_t));

       sprintf( str, "GroupStatus(%d)", groupid );
       Tcl_LinkVar( TCLInterp,str,(char *)&group->status,TCL_LINK_INT );

       group->status = 1;
       group->Next = NULL;
       group->Name = (char *)malloc(strlen(name)+1);
       strcpy(group->Name,name);
   }

   group->Open = open;

   return groupid;
}

/*******************************************************************************
 *
 *     Name:         geo_add_vertex( geometry_t *, vertex_t * )
 *
 *     Purpose:      Add a vertex to triangulation description
 *
 *     Parameters: 
 *
 *         Input:    (geometry_t *) pointer to geometry structure to modify
 *                   (vertex_t   *) pointer to vertex to add
 *
 *         Output:   (geometry_t *) is modified
 *
 *   Return value:   TRUE if success, FALSE if malloc() fails
 *
 ******************************************************************************/
int geo_add_vertex( geometry_t *geometry, vertex_t *vertex )
{
    vertex_t *ptr;

    if ( geometry->VertexCount >= geometry->MaxVertexCount )
    {
        geometry->MaxVertexCount += GEO_VERTEX_BLOCK_SIZE;
        ptr = (vertex_t *)realloc(geometry->Vertices,geometry->MaxVertexCount*sizeof(vertex_t));

        if ( !ptr )
        {
            fprintf( stderr, "FATAL: geo_add_vertex: Can't allocate memory.\n" );
            return -1;
        }

        geometry->Vertices = ptr;
    }

    vertex->Faces = NULL;
    geometry->Vertices[geometry->VertexCount++] = *vertex;

    return (geometry->VertexCount-1);
}

/*******************************************************************************
 *
 *     Name:         geo_add_edge( geometry_t *, edge_t *,element_t * )
 *
 *     Purpose:      Add an edge to edge table
 *
 *     Parameters: 
 *          Input:    (geometry_t *) pointer to geometry structure to modify
 *                    (edge_t   *) pointer to edge to add
 *                    (element_t   *) pointer to parent element
 *
 *         Output:   (geometry_t *) is modified
 *
 *   Return value:   TRUE if success, FALSE if malloc() fails
 *
 ******************************************************************************/
int geo_add_edge( geometry_t *geometry,int v0,int v1,element_t *elm )
{
    int swap;
    edge_list_t *ptr;

    if ( !GlobalOptions.VolumeEdges && elm->ElementType->ElementCode>=500 ) return TRUE;

    if ( v0 > v1 ) { swap = v0; v0 = v1; v1 = swap; }

    if ( v0 >= geometry->VertexCount || v1 >= geometry->VertexCount )
    {
       fprintf( stderr, "geo_add_edge: vertex number [%d,%d] out of range\n", v0,v1);
       return FALSE;
    }


    for( ptr=geometry->Edges[v0].EdgeList; ptr != NULL; ptr=ptr->Next )
    {
       if ( v1 == ptr->Entry ) break;
    }

    if ( ptr != NULL )
    {
        if ( elm->ElementType->ElementCode < 300 ) 
        {
           ptr->Count = -1;
           ptr->Element = elm;
        } else if ( ptr->Count > 0 ) ptr->Count++;
    } else
    {
        ptr = (edge_list_t *)malloc( sizeof(edge_list_t) );
        if ( !ptr )
        {
            fprintf( stderr, "geo_add_edge: FATAL: can't alloc memory.\n" );
            return FALSE;
        }
        if ( elm->ElementType->ElementCode < 300 ) 
          ptr->Count = -1;
        else
          ptr->Count =  1;
        ptr->Entry = v1;
        ptr->Element = elm;
        ptr->Next = geometry->Edges[v0].EdgeList;

        geometry->Edges[v0].EdgeList = ptr;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:         geo_free_edge_tables( geometry_t * )
 *
 *     Purpose:      free geometry edge tables
 *
 *     Parameters: 
 *          Input:    (geometry_t *) pointer to geometry structure to modify
 *
 *         Output:   (geometry_t *) is modified
 *
 *   Return value:   void
 *
 ******************************************************************************/
void geo_free_edge_tables( geometry_t *geometry )
{
    edge_list_t *ptr,*ptr1;
    int i;

    if ( geometry->Edges )
    {
        for( i=0; i<geometry->VertexCount; i++ )
        {
            for( ptr=geometry->Edges[i].EdgeList; ptr != NULL; )
            {
                ptr1 = ptr->Next;
                free( ptr );
                ptr = ptr1;
            }
            geometry->Edges[i].EdgeList = NULL;
        }
        free( geometry->Edges );
        geometry->Edges = NULL;
    }
}

/*******************************************************************************
 *
 *     Name:         geo_free_vertex_face_tables( geometry_t * )
 *
 *     Purpose:      free geometry vertex face tables
 *
 *     Parameters: 
 *          Input:    (geometry_t *) pointer to geometry structure to modify
 *
 *         Output:   (geometry_t *) is modified
 *
 *   Return value:   void
 *
 ******************************************************************************/
void geo_free_vertex_face_tables( geometry_t *geometry )
{
    vertex_face_t *face,*face1;
    vertex_t *vertex;
    int i;

    for( i=0; i<geometry->VertexCount; i++ )
    {
        vertex = &geometry->Vertices[i];
        for( face=vertex->Faces; face != NULL; )
        {
            face1 = face->Next;
            free( face );
            face = face1;
        }
        vertex->Faces = NULL;
    }
}

/*******************************************************************************
 *
 *     Name:         geo_add_triangle( geometry_t *, triangle_t * )
 *
 *     Purpose:      Add a  triangle to  triangulation description.  The vertices
 *                   to which (triangle_t *) points must exist before this routine
 *                   is called.
 *
 *     Parameters: 
 *
 *         Input:    (geometry_t *) pointer to geometry structure to modify
 *                   (triangle_t *) pointer to triangle ito add
 *
 *         Output:   (geometry_t *) is modified
 *
 *   Return value:    TRUE if success, FALSE if malloc() fails
 *
 ******************************************************************************/
int geo_add_triangle( geometry_t *geometry, triangle_t *triangle )
{
    triangle_t *ptr,*tri;
    vertex_t *vertex;
    vertex_face_t *face;

    int i,j,v0=triangle->v[0],v1=triangle->v[1],v2=triangle->v[2],swap,w0,w1,w2;

    if ( v0 >= geometry->VertexCount || v1 >= geometry->VertexCount || v2 >= geometry->VertexCount)
    {
        fprintf( stderr, "geo_add_triangle: vertex number [%d,%d,%d] out of range\n", v0,v1,v2);
        return FALSE;
    }

    if ( v0 > v1 ) { swap=v0; v0=v1; v1=swap; }
    if ( v1 > v2 ) { swap=v1; v1=v2; v2=swap; }
    if ( v0 > v1 ) { swap=v0; v0=v1; v1=swap; }

    vertex = &geometry->Vertices[v0];
    for( face=vertex->Faces; face != NULL; face=face->Next )
    {
        tri = &geometry->Triangles[face->Face];

        if ( v0 == tri->v[0] ) 
           if ( v1 == tri->v[1] ) 
              if ( v2 == tri->v[2] )
              {
                if ( triangle->Element->ElementType->ElementCode < 500 )
                {
                   tri->Count = - 1;
                   tri->Element = triangle->Element;
                } else if ( tri->Count > 0 ) tri->Count++;
                return TRUE;
              }
    }

    if ( geometry->TriangleCount >= geometry->MaxTriangleCount )
    {
        geometry->MaxTriangleCount += GEO_TRIANGLE_BLOCK_SIZE;
        ptr = (triangle_t *)realloc( geometry->Triangles,geometry->MaxTriangleCount*sizeof(triangle_t) );
        if ( !ptr )
        {
            fprintf( stderr, "FATAL: geo_add_triangle: Can't allocate memory.\n" );
            return FALSE;
        }

        geometry->Triangles = ptr;
    }

    if ( GlobalOptions.SurfaceSides ) {
      if ( triangle->Edge[0] ) geo_add_edge( geometry,triangle->v[0],triangle->v[1],triangle->Element );
      if ( triangle->Edge[1] ) geo_add_edge( geometry,triangle->v[1],triangle->v[2],triangle->Element );
      if ( triangle->Edge[2] ) geo_add_edge( geometry,triangle->v[2],triangle->v[0],triangle->Element );
    }

    triangle->v[0] = v0;
    triangle->v[1] = v1;
    triangle->v[2] = v2;

    for( j=0; j<3; j++ )
    {
        vertex = &geometry->Vertices[triangle->v[j]];

        face = (vertex_face_t *)malloc( sizeof(vertex_face_t) );
        if ( !face )
        {
            fprintf( stderr, "geo_add_triangle: FATAL: malloc() failed.\n" );
            return FALSE;
        }

        face->Face = geometry->TriangleCount;
        face->Next = vertex->Faces;
        vertex->Faces = face;
    }

    if ( triangle->Element->ElementType->ElementCode < 500 )
      triangle->Count = -1;
    else
      triangle->Count =  1;
    geometry->Triangles[geometry->TriangleCount++] = *triangle;
 
    return TRUE;
}

/*******************************************************************************
 *
 *     Name:         geo_triangle_normal
 *
 *     Purpose:      generate triangle face normal by cross product
 *
 *     Parameters:
 *
 *         Input:    (geometry_t *) holding triangulation information
 *                   (triangle_t *) the triangle for which the normal
 *                                  is computed
 *                                  (this could be index to geometry,hmm)
 *
 *         Output:   (triangle_t *) is modified
 *
 *   Return value:   void
 *
 ******************************************************************************/
void geo_triangle_normal( geometry_t *geom,triangle_t *triangle )
{
    int i;

    double ax,ay,az,bx,by,bz,u,v,w,s,uu,vv,ww;

    vertex_t *Vert = geom->Vertices;

    int v0 = triangle->v[0];
    int v1 = triangle->v[1];
    int v2 = triangle->v[2];

    if ( v0 >= geom->VertexCount || v1 >= geom->VertexCount || v2 >= geom->VertexCount)
    {
        fprintf( stderr, "geo_triangle_normal: vertex number [%d,%d,%d] out of range\n", v0,v1,v2);
        return;
    }

    ax = Vert[triangle->v[1]].x[0] - Vert[triangle->v[0]].x[0];
    ay = Vert[triangle->v[1]].x[1] - Vert[triangle->v[0]].x[1];
    az = Vert[triangle->v[1]].x[2] - Vert[triangle->v[0]].x[2];

    bx = Vert[triangle->v[2]].x[0] - Vert[triangle->v[0]].x[0];
    by = Vert[triangle->v[2]].x[1] - Vert[triangle->v[0]].x[1];
    bz = Vert[triangle->v[2]].x[2] - Vert[triangle->v[0]].x[2];

    u = ay*bz - az*by;
    v = az*bx - ax*bz;
    w = ax*by - ay*bx;
    s = 1.0 / sqrt(u*u + v*v + w*w);

    u = u * s;
    v = v * s;
    w = w * s;

    for( i=0; i<3; i++ )
    {
        triangle->u[i][0] = u;
        triangle->u[i][1] = v;
        triangle->u[i][2] = w;
    }

    triangle->Fu[0] = u;
    triangle->Fu[1] = v;
    triangle->Fu[2] = w;
}


/*******************************************************************************
 *
 *     Name:         geo_vertex_normals( geometry_t *, double )
 *
 *     Purpose:      Generete triangle vertex normals from face normals,
 *                   by averaging face normals connected to same vertex
 *                   and within a given angle.
 *
 *     Parameters:
 *
 *         Input:    (geometry_t *) pointer to geometry structure to modify
 *                   (double)       threshold angle
 *
 *         Output:   (geometry_t *) is modified
 *
 *  Return value:    void
 *
 ******************************************************************************/
void geo_vertex_normals( geometry_t *geometry, double Ang ) 
{
    vertex_t *vertex,*p = geometry->Vertices;

    vertex_face_t *Faces;

    triangle_t *t = geometry->Triangles;

    double CU,CV,CW,u,v,w,s,A;

    int i,j,k;

    A  = cos( M_PI*Ang / 180.0 );

    for( i=0; i<geometry->TriangleCount; i++ ) geo_triangle_normal( geometry,&t[i] );

    for( i=0; i<geometry->TriangleCount; i++ ) /* loop through faces */
    {
        u = t[i].Fu[0];
        v = t[i].Fu[1];
        w = t[i].Fu[2];

        for( j=0; j<3; j++ )                   /* trough face vertices */
        {
            t[i].u[j][0] = u;
            t[i].u[j][1] = v;
            t[i].u[j][2] = w;

            vertex = &p[t[i].v[j]];

            /*
             * loop trough faces connected to a vertex
             */
            for( Faces=vertex->Faces; Faces != NULL; Faces=Faces->Next )
            {
                if ( Faces->Face == i ) continue;

                CU = t[Faces->Face].Fu[0];
                CV = t[Faces->Face].Fu[1];
                CW = t[Faces->Face].Fu[2];

                s = u*CU + v*CV + w*CW;

                if ( ABS(s) >= A )
                {
                    if ( s < 0.0 )
                    {
                        t[i].u[j][0] -= CU;
                        t[i].u[j][1] -= CV;
                        t[i].u[j][2] -= CW;
                    } else
                    {
                        t[i].u[j][0] += CU;
                        t[i].u[j][1] += CV;
                        t[i].u[j][2] += CW;
                    }
                }
            }
        }
    }

    /*
     *   normalize
     */
    for( i=0; i<geometry->TriangleCount; i++ )
        for( j=0; j<3; j++ )
        {
            u = t[i].u[j][0];
            v = t[i].u[j][1];
            w = t[i].u[j][2];
            s = sqrt( u*u + v*v + w*w );
 
            if ( s > 1.0E-10 )
            {
                t[i].u[j][0] /= s;
                t[i].u[j][1] /= s;
                t[i].u[j][2] /= s;
            }
        }
}
