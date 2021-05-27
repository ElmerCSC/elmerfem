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
 * Action routines for the mesh visual class.
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
 *
 * Modification history:
 *
 * 28 Sep 1995, modified vis_initialize_mesh_visual to set the VisualName field
 *             of the visual_type structure
 *
 * Juha R.
 *
 ******************************************************************************/

#include "../elmerpost.h"


/******************************************************************************
 * 
 * Parameter sructure definitios for mesh visual class
 *
 ******************************************************************************/
static char *mesh_style_names[] =
{
   "none", "line", "surf", "line_and_surf", NULL
};

typedef struct mesh_s
{
    scalar_t   *ColorData;

    mesh_style_t Style;

    edge_style_t EdgeStyle;

    material_t *Material;
    material_t *EdgeMaterial;
    colormap_t *ColorMap;

    int LineQuality;
    double LineWidth;

    line_style_t LineStyle;
 
    logical_t NodeNumbers;
} mesh_t;


/*******************************************************************************
 *
 *     Name:         vis_polygon
 *
 *     Purpose:      
 *                  
 *
 *     Parameters:
 *
 *         Input:    (polygon_t *) polygon
 *
 *         Output:   graphics
 *
 *   Return value:   void
 *
 ******************************************************************************/
void vis_polygon( polygon_t *poly)
{
    int n=3;

    gra_poly3( n,poly->x,poly->y,poly->z,poly->u,poly->v,poly->w,poly->c );
}

/*******************************************************************************
 *
 *     Name:         vis_triangle
 *
 *     Purpose:      
 *                  
 *
 *     Parameters:
 *
 *         Input:    (triangle_t *) triangle
 *                   (vertex_t   *) vertex array
 *                   (double     *) quantity to use color the edges (or NULL)
 *                   (double,double ) CScl,CAdd are constatnt to scale color
 *                                  range (0-1)
 *                   (line_style_t) line style, either line_style_line or
 *                                  line_style_cylider
 *
 *         Output:   graphics
 *
 *   Return value:   void
 *
 ******************************************************************************/
void vis_triangle
   (
      triangle_t *t,vertex_t *v,double *color,double CScl,double CAdd
   )
{
    float x[3][3],n[3][3],c[3];
    int j,k;

    for( j=0; j<3; j++ )
    {
        k = t->v[j];
        x[j][0] = v[k].x[0];
        x[j][1] = v[k].x[1];
        x[j][2] = v[k].x[2];

        n[j][0] = t->u[j][0];
        n[j][1] = t->u[j][1];
        n[j][2] = t->u[j][2];

        if ( color ) c[j] = CScl * ( color[k]-CAdd );
    }

    gra_triangle( x,n,c );
}

/*******************************************************************************
 *
 *     Name:         vis_draw_edge
 *
 *     Purpose:      draw element edge as line or cylinder
 *
 *     Parameters:
 *
 *         Input:    (vertex_t   *) vertex array
 *                   (int,int)      edge vertex indices
 *                   (double *)      color quantity
 *                   (double,double ) CScl,CAdd are constatnt to scale color
 *                                  range (0-1)
 *                   (line_style_t) line style, either line_style_line or
 *                                  line_style_cylider
 *                   (double)       line width, cylinder radius.
 *                   
 *
 *         Output:   graphics
 *
 *   Return value:   void
 *
 ******************************************************************************/
static void vis_draw_edge(vertex_t *vertex,int v0,int v1,double *color,double CScl,
                       double CAdd,line_style_t style,double width)
{
    double c0=0.0,c1=1.0;

    /*
     *  if color function given scale values to 0-1
     */
    if ( color )
    {
        c0 = CScl*( color[v0] - CAdd );
        c1 = CScl*( color[v1] - CAdd );
    }
     
    gra_line( vertex[v0].x,c0,vertex[v1].x,c1,style,width );
}

/*******************************************************************************
 *
 *     Name:        vis_mesh
 *
 *     Purpose:     draw mesh as lines or surface, with color codes or not
 *
 *     Parameters: 
 *
 *         Input:   (geometry_t *) triangles to draw
 *                  (mesh_t     *) mesh display parameters
 *                  (double) 
 *
 *         Output:  graphics
 *   
 *   Return value:  if mouse interaction is going on, and time used exeeds
 *                  given value (TooLong1,2) FALSE, otherwise TRUE
 *
 ******************************************************************************/
static int vis_mesh( geometry_t *geometry, element_model_t *model, mesh_t *Mesh,double dt )
{
    scalar_t *ColorData = Mesh->ColorData;

    vertex_t *v = geometry->Vertices;

    edge_list_t *edge;

    double width = Mesh->LineWidth*0.005,CScl=1.0,CAdd=0.0,*C=NULL;

    int i,j,quick,N=geometry->VertexCount,NT=geometry->TriangleCount;

    element_t *elements = model->Elements;

    static char str[100];

    if ( !GlobalOptions.StereoMode )
      if ( Mesh->Material->Diffuse[3]  < 1.0 )
      {
          if ( GlobalPass != 0 ) return TRUE;
      } else if ( GlobalPass == 0 )
      {
          return TRUE;
      }

    quick  = (Mesh->Style == mesh_style_line && Mesh->LineStyle == line_style_line );
    quick &= ~Mesh->NodeNumbers;
    quick |= epMouseDown && epMouseDownTakesTooLong;

    gra_set_material( Mesh->Material );

    if ( quick && !(epMouseDown && epMouseDownTakesTooLong) )
    {
        gra_line_width( Mesh->LineWidth );
    } else {
        gra_line_width( 1.0 );
    }

    if ( ColorData && ColorData->f )
    {
        C = ColorData->f;

        CAdd = ColorData->min;
        if ( ABS(ColorData->max - ColorData->min)>0.0 )
            CScl = 1.0 / (ColorData->max - ColorData->min);
        else
            CScl = 1.0;

        gra_set_colormap( Mesh->ColorMap );
    } else gra_set_colormap( NULL );

    if ( quick )
    {
        gra_set_material( Mesh->EdgeMaterial );
        gra_beg_lines();

        for( i=0; i<N; i++ )
        {
            if ( v[i].ElementModelNode )
            {
                for( edge=geometry->Edges[i].EdgeList; edge != NULL; edge = edge->Next )
                {
                    if ( edge->Element && !edge->Element->DisplayFlag ) continue;
                    if ( Mesh->EdgeStyle == edge_style_all || ABS(edge->Count) == 1 )
                    {
                         vis_draw_edge( v,i,edge->Entry,C,CScl,CAdd,line_style_line,width );
                    }
                }
            } else break;

            if ( epMouseDown && (i & 128) )
            {
                if ( RealTime() - dt > TooLong2 )
                    if ( ++epMouseDownTakesTooLong > 3 )
                    {
                        gra_end_lines();
                        return FALSE;
                    } else dt = RealTime();
            }

            if ( BreakLoop ) break;
        }

        gra_end_lines();
        gra_set_material( Mesh->Material );

        return TRUE;
    }

    if ( Mesh->NodeNumbers )
    {
       for( i=0; i<N; i++ )
       {
          if ( v[i].ElementModelNode )
          {
            glRasterPos3f( v[i].x[0],v[i].x[1],v[i].x[2]);
            sprintf( str, "N %d", i );
            PrintString( str );
          }
       }
    }

    if ( Mesh->Style == mesh_style_surf || Mesh->Style == mesh_style_line_and_surf )
    {
        triangle_t *t = geometry->Triangles;

        gra_begin( GRA_TRIANGLES );
        for( i=0; i<NT; i++,t++ )
        {
            if ( t->Element && !t->Element->DisplayFlag ) continue;

            if ( Mesh->EdgeStyle != edge_style_all && t->Count > 1 ) continue;

            vis_triangle( t,v,C,CScl,CAdd );

            if ( epMouseDown && (i & 8) )
                if ( RealTime() - dt > TooLong1 )
                {
                    epMouseDownTakesTooLong++;
                    gra_end();
                    return FALSE;
                }

            if ( BreakLoop ) break;
        }
        gra_end();
    }

    if ( Mesh->Style == mesh_style_line_and_surf ) gra_set_colormap( NULL );

    if ( Mesh->Style == mesh_style_line || Mesh->Style == mesh_style_line_and_surf )
    {

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glTranslatef(0.0,0.0,0.005);
        glMatrixMode(GL_MODELVIEW);

        gra_set_material( Mesh->EdgeMaterial );

        if ( Mesh->LineStyle == line_style_cylinder )
        {
            gra_sphere_quality( Mesh->LineQuality );

            for( i=0; i<N; i++ )
            {
                if ( v[i].ElementModelNode )
                {
                    for( edge=geometry->Edges[i].EdgeList; edge != NULL; edge = edge->Next )
                    {
                        if ( edge->Element && !edge->Element->DisplayFlag ) continue;
                        if ( Mesh->EdgeStyle == edge_style_all || ABS(edge->Count) == 1 )
                        {
                             j = edge->Entry;
                             vis_draw_edge( v,i,j,C,CScl,CAdd,line_style_cylinder,width );
                        }
                   }
                } else break;

                if ( epMouseDown && (i & 8) )
                    if ( RealTime() - dt > TooLong1 )
                    {
                        epMouseDownTakesTooLong++;
                        glMatrixMode(GL_PROJECTION);
                        glPopMatrix();
                        glMatrixMode( GL_MODELVIEW );
                        return FALSE;
                    }

                if ( BreakLoop ) break;
            }
        } else
        {
            gra_line_width( Mesh->LineWidth );
            gra_beg_lines();
            for( i=0; i<N; i++ )
            {
                if ( v[i].ElementModelNode )
                {
                    for( edge=geometry->Edges[i].EdgeList; edge != NULL; edge = edge->Next )
                    {
                        if ( edge->Element && !edge->Element->DisplayFlag ) continue;
                        if ( Mesh->EdgeStyle == edge_style_all || ABS(edge->Count) == 1 )
                        {
                             vis_draw_edge( v,i,edge->Entry,C,CScl,CAdd,line_style_line,width );
                        }
                    }
                } else break;

                if ( epMouseDown && (i & 32) )
                    if ( RealTime() - dt > TooLong1 )
                    {
                        epMouseDownTakesTooLong++;
                        gra_end_lines();
                        glMatrixMode(GL_PROJECTION);
                        glPopMatrix();
                        glMatrixMode( GL_MODELVIEW );
                        return FALSE;
                    }

                if ( BreakLoop ) break;
            }
            gra_end_lines();
        }
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode( GL_MODELVIEW );
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        vis_mesh_alloc
 *
 *     Purpose:     allocate memory for mesh_t structure
 *
 *     Parameters: 
 *
 *         Input:   none
 *
 *         Output:  none
 *   
 *   Return value:  pointer to allocated memory
 *
 ******************************************************************************/
static mesh_t *vis_mesh_alloc()
{
     mesh_t *mesh = (mesh_t *)calloc(sizeof(mesh_t),1);

     if ( !mesh )
     {
         fprintf( stderr, "vis_mesh_alloc: FATAL: can't alloc a few bytes of memory\n" );
     }

     return mesh;
}

/*******************************************************************************
 *
 *     Name:        vis_mesh_delete
 *
 *     Purpose:     free memory associated with mesh_t structure
 *
 *     Parameters: 
 *
 *         Input:   (mesh_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_mesh_delete(mesh_t *mesh)
{
    if ( mesh ) free( mesh );
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_mesh_visual
 *
 *     Purpose:     Register "Mesh" visual type
 *
 *     Parameters: 
 *
 *         Input:   none
 *
 *         Output:  none
 *   
 *   Return value:  vis_add_visual_type (malloc success probably)...
 *
 ******************************************************************************/
int vis_initialize_mesh_visual()
{
    static char *visual_name = "Mesh";
    visual_type_t VisualDef; 

    static mesh_t mesh;

    static visual_param_t MeshParams[] =
    {
        { "ColorData",     "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "Style",         "%d", 0, VIS_VISUAL_PARAM_INT,     mesh_style_line, 0.0, NULL },
        { "EdgeStyle",     "%d", 0, VIS_VISUAL_PARAM_INT,     edge_style_all,  0.0, NULL },
        { "Material",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { "EdgeMaterial",  "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultEdgeMaterial },
        { "ColorMap",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "LineWidth",     "%lf", 0, VIS_VISUAL_PARAM_FLOAT,  0, 1.0, NULL },
        { "LineQuality",   "%d", 0, VIS_VISUAL_PARAM_INT,     1, 0.0, NULL },
        { "LineStyle",     "%d", 0, VIS_VISUAL_PARAM_INT,     0, 0.0, NULL },
        { "NodeNumbers",   "%d", 0, VIS_VISUAL_PARAM_LOGICAL, 0, 0.0, NULL },
        { NULL, NULL, 0, 0, 0, 0.0, NULL }
    };

    int n = 0;

    MeshParams[n++].Offset = (char *)&mesh.ColorData     - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.Style         - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.EdgeStyle     - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.Material      - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.EdgeMaterial  - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.ColorMap      - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.LineWidth     - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.LineQuality   - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.LineStyle     - (char *)&mesh;
    MeshParams[n++].Offset = (char *)&mesh.NodeNumbers   - (char *)&mesh;

    VisualDef.VisualName    = visual_name;

    VisualDef.RealizeVisual = (int   (*)()) vis_mesh;
    VisualDef.AllocParams   = (void *(*)()) vis_mesh_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_mesh_delete;
    VisualDef.VisualParams  = MeshParams;

    return vis_add_visual_type( &VisualDef );
}
