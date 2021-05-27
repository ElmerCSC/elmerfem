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
 * Action routines for the sphere visual class.
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
 *                       Date: 28 Sep 1995
 *
 * Modification history:
 *
 * 28 Sep 1995, modified vis_initialize_sphere_visual to set the VisualName
 *              field of the visual_type structure
 *
 * 29 Sep 1995, modified vis_initialize_sphere_visual to set new field
 *              VisualParams of the visual_type structure, and removed
 *              vis_spehre_set_param routine, the program now uses the
 *              visual_params_t structure and the routine vis_set_param
 *              in file visual.c
 *
 * Juha R.
 *
 ******************************************************************************/

#include "../elmerpost.h"
/******************************************************************************
 * 
 * Parameter sructure definitios for sphere visual class
 *
 ******************************************************************************/
typedef struct shpere_s
{
    scalar_t *ColorData;
    scalar_t *RadiusData;
 
    double RadiusScale;

    scalar_t *ThresholdData;
    double Floor,Ceiling;

    material_t *Material;
    colormap_t *ColorMap;

    int Quality;
} sphere_t;

/*******************************************************************************
 *
 *     Name:        vis_sphere
 *
 *     Purpose:     draw sphere as lines or surface, with color coded or not
 *
 *     Parameters: 
 *
 *         Input:   (geometry_t *)   geometry description
 *                  (sphere_t     *) sphere display parameters
 *                  (double) 
 *
 *         Output:  graphics
 *   
 *   Return value:  if mouse interaction is going on, and time used exeeds
 *                  given value (TooLong1,2) FALSE, otherwise TRUE
 *
 ******************************************************************************/
static int vis_sphere( geometry_t *geometry, element_model_t *model, sphere_t *Sphere,double dt )
{
    scalar_t *RadiusData = Sphere->RadiusData;
    scalar_t *ColorData  = Sphere->ColorData;
    scalar_t *ThresholdData = Sphere->ThresholdData;

    vertex_t *v = geometry->Vertices;

    vertex_face_t *face;

    int i,j,quick,N=geometry->VertexCount;

    double *C=NULL, *R=NULL, *T=NULL;
    double CScl=1.0,CAdd=0.0,RScl=1.0,RAdd=0.0,Rad=1.0,Col=1.0;

    if ( !GlobalOptions.StereoMode )
      if ( Sphere->Material->Diffuse[3]  < 1.0 )
      {
          if ( GlobalPass != 0 ) return TRUE;
      } else if ( GlobalPass == 0 )
      {
          return TRUE;
      }

    if ( ColorData && ColorData->f )
    {
        C = ColorData->f;
        CAdd = ColorData->min;
        CScl = 1.0/(ColorData->max - ColorData->min);

        gra_set_colormap( Sphere->ColorMap );
    } else gra_set_colormap( NULL );

    if ( RadiusData && RadiusData->f )
    {
        R = RadiusData->f;
        RAdd = RadiusData->min;
        RScl = 1.0/(RadiusData->max - RadiusData->min);
    }
    else { RAdd = 0.0; RScl = 0.05*Sphere->RadiusScale; }

    if ( ThresholdData && ThresholdData->f ) T = ThresholdData->f;

    quick = (epMouseDown && epMouseDownTakesTooLong) || (Sphere->Quality < 0);

    if ( quick )
    {
        gra_begin( GRA_POINTS );
#if 0
        gra_polygon_mode( GRA_LINE );
        gra_sphere_quality( 1 );
#endif
    } else {
       gra_sphere_quality( Sphere->Quality );
       if ( Sphere->Quality == 0 ) gra_polygon_mode( GRA_LINE );
    }

    gra_set_material( Sphere->Material );

    for( i=0; i<N; i++ )
    {
        if ( v[i].ElementModelNode )
        {
            for( face=v[i].Faces; face!=NULL; face=face->Next )
            {
               if ( geometry->Triangles[face->Face].Element->DisplayFlag ) break;
            }
#if 1
            if ( v[i].Faces && !face ) continue;
#else
            if ( !face ) continue;
#endif

            if ( T && (T[i]<Sphere->Floor || T[i]>Sphere->Ceiling) ) continue;

            if ( R ) 
                Rad = RScl*(R[i]-RAdd)+0.05;
            else Rad = 1.0;

            Rad *= 0.05*Sphere->RadiusScale;

            if ( C ) Col = CScl*(C[i] - CAdd); else Col = CScl;

            if ( quick ) {
               gra_point( v[i].x[0],v[i].x[1],v[i].x[2],Col,Rad );
            } else {
               gra_sphere( v[i].x[0],v[i].x[1],v[i].x[2],Col,Rad );
            }

        } else break;

        if ( epMouseDown && (i & 8) )
        {
            if ( !epMouseDownTakesTooLong )
            {
                if ( RealTime() - dt > TooLong1 )
                {
                    if ( quick ) gra_end();
                    gra_polygon_mode( GRA_FILL );
                    ++epMouseDownTakesTooLong;
                    return FALSE;
                }
            }
            else if ( RealTime() - dt > TooLong2 )
                if ( ++epMouseDownTakesTooLong > 3 )
                {
                    gra_end();
                    gra_polygon_mode( GRA_FILL );
                    return FALSE;
                } else dt = RealTime();
        }

        if ( BreakLoop ) break;
    }

    if ( quick ) gra_end();
    gra_polygon_mode( GRA_FILL );

    return TRUE;
}


/*******************************************************************************
 *
 *     Name:        vis_sphere_alloc
 *
 *     Purpose:     allocate memory for sphere_t structure
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
static sphere_t *vis_sphere_alloc()
{
     sphere_t *sphere = (sphere_t *)calloc(sizeof(sphere_t),1);

     if ( !sphere )
     {
         fprintf( stderr, "vis_sphere_alloc: FATAL: can't alloc a few bytes of memory\n" );
     }

     return sphere;
}

/*******************************************************************************
 *
 *     Name:        vis_sphere_delete
 *
 *     Purpose:     free memory associated with sphere_t structure
 *
 *     Parameters: 
 *
 *         Input:   (sphere_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_sphere_delete(sphere_t *sphere)
{
    if ( sphere ) free( sphere );
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_sphere_visual
 *
 *     Purpose:     Register "Spheres" visual type
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
int vis_initialize_sphere_visual()
{
    static char *visual_name = "Spheres";
    visual_type_t VisualDef;

    static sphere_t sphere;

    int n = 0;

    static visual_param_t SphereParams[] =
    {
        { "ColorData",     "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "RadiusData",    "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "RadiusScale",   "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "ThresholdData", "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 1.0, NULL },
        { "Floor",         "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 0.0, NULL },
        { "Ceiling",       "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "Material",      "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { "ColorMap",      "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "Quality",       "%d",  0, VIS_VISUAL_PARAM_INT,     1, 0.0, NULL },
        { NULL, NULL, 0, 0, 0,  0.0, NULL }
    };

    SphereParams[n++].Offset = (char *)&sphere.ColorData     - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.RadiusData    - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.RadiusScale   - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.ThresholdData - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.Floor         - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.Ceiling       - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.Material      - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.ColorMap      - (char *)&sphere;
    SphereParams[n++].Offset = (char *)&sphere.Quality       - (char *)&sphere;

    VisualDef.VisualName    =  visual_name;
    VisualDef.RealizeVisual = (int   (*)()) vis_sphere;
    VisualDef.AllocParams   = (void *(*)()) vis_sphere_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_sphere_delete;
    VisualDef.VisualParams  = SphereParams;

    return vis_add_visual_type( &VisualDef );
}
