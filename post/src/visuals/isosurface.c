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
 * Action routines for visual class IsoSurface.
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
 * 28 Sep 1995, modified vis_initialize_isosurface_visual to set the
 *              VisualName field of the visual_type structure
 *
 * Juha R.
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include "../elements/elements.h"
/******************************************************************************
 * 
 * Parameter sructure definitios for Isosurface visual class
 *
 ******************************************************************************/


typedef struct isosurface_s
{
    mesh_style_t Style;

    scalar_t *ContourData;
    scalar_t *ColorData;

    scalar_t *NormalData[3];

    int NofLevels;
    double *Levels;

    material_t *Material;
    colormap_t *ColorMap;

    int Recompute;

    int NPolygons;
    int MaxPolygons;
    polygon_t *Polygons;

    int LineQuality;
    double LineWidth;

    line_style_t LineStyle;
} isosurface_t;

#define FEPS 1.0E-9

/*******************************************************************************
 *
 *     Name:        vis_get_isosurface
 *
 *     Purpose:     Extract a surface with given threshold 
 *
 *     Parameters: 
 *
 *         Input:   (triangle_t *)
 *                  (vertex_t   *)
 *                  (double *,double,double)  color quantity, and scales => 0,1 
 *                  (double *)             surface quantity
 *                  (double)               threshold value
 *
 *         Output:  (line_t *)             place to store the line
 *   
 *   Return value:  number of points generated, line exists if (n>=2)
 *
 ******************************************************************************/
static int vis_get_isosurface
   (
            element_model_t *model, element_t *element, vertex_t *vertices, double *C, double *F, 
     double *NX, double *NY, double *NZ, polygon_t *Polygons, int nlevels, double *levels,double CScl,double CAdd
   )
{
    static double x[ELM_MAX_ELEMENT_NODES];
    static double y[ELM_MAX_ELEMENT_NODES];
    static double z[ELM_MAX_ELEMENT_NODES];
    static double f[ELM_MAX_ELEMENT_NODES];
    static double c[ELM_MAX_ELEMENT_NODES];
    static double u[ELM_MAX_ELEMENT_NODES];
    static double v[ELM_MAX_ELEMENT_NODES];
    static double w[ELM_MAX_ELEMENT_NODES];

    int i,j,n,*T=element->Topology;
    element_type_t *elmt = element->ElementType;

    if ( elmt->IsoSurface ) 
    {
        for( i=0; i<elmt->NumberOfNodes; i++ )
        {
            x[i] = vertices[T[i]].x[0];
            y[i] = vertices[T[i]].x[1];
            z[i] = vertices[T[i]].x[2];

            f[i] = F[T[i]];
            if ( C )
                c[i] = CScl*(C[T[i]]-CAdd);
             else
                c[i] = 0.0;

            if ( NX )
            {
                u[i] = NX[T[i]];
                v[i] = NY[T[i]];
                w[i] = NZ[T[i]];
            }
        }

        n = 0;
        for( i=0; i<nlevels; i++ )
        {
            n += (*elmt->IsoSurface)( levels[i],f,c,x,y,z,u,v,w,&Polygons[n] );
        } 

        return n;
    }

    return 0;
}

/*******************************************************************************
 *
 *     Name:        vis_isosurfaces
 *
 *     Purpose:     Draw surface given data, and threshold values
 *
 *     Parameters: 
 *
 *         Input:   (geometry_t *)  
 *                  (isosurface_t *) surface display params
 *                  (double)            time used
 *
 *         Output:  graphics
 *   
 *   Return value:  if mouse interaction is going on, and time used exeeds
 *                  given value (TooLong1,2) FALSE, otherwise TRUE
 *
 ******************************************************************************/
static int vis_isosurface( geometry_t *geometry, element_model_t *model, isosurface_t *IsoSurface, double dt )
{
    double CScl,CAdd,*Levels = IsoSurface->Levels;
    double *C=NULL, *F=NULL, *NX = NULL, *NY = NULL, *NZ = NULL;

    scalar_t *ColorData   = IsoSurface->ColorData;
    scalar_t *ContourData = IsoSurface->ContourData;

    int i,j,k,l,m,n,quick;

    double width = IsoSurface->LineWidth*0.005, color[3];

    element_t *elements = model->Elements;

    if ( !ContourData || !ContourData->f ) return TRUE;

    if ( !GlobalOptions.StereoMode )
      if ( IsoSurface->Material->Diffuse[3]  < 1.0 )
      {
          if ( GlobalPass != 0 ) return TRUE;
      } else if ( GlobalPass == 0 )
      {
          return TRUE;
      }

    F = ContourData->f;

    if ( ColorData && ColorData->f )
    {
        CScl = 1.0/(ColorData->max - ColorData->min);
        C = ColorData->f;
        CAdd = ColorData->min;

        gra_set_colormap( IsoSurface->ColorMap );
    } else gra_set_colormap( NULL );

    if ( IsoSurface->NormalData[0] && IsoSurface->NormalData[0]->f )
    {
        NX = IsoSurface->NormalData[0]->f;
        NY = IsoSurface->NormalData[1]->f;
        NZ = IsoSurface->NormalData[2]->f;
    } else
    {
        gra_loff();
    }

    quick = (IsoSurface->Style == mesh_style_line);
    quick |= epMouseDown && epMouseDownTakesTooLong;

    if ( !quick && (IsoSurface->LineStyle == line_style_cylinder) )
    {
        gra_sphere_quality( IsoSurface->LineQuality );
    }

    if ( quick && !(epMouseDown && epMouseDownTakesTooLong) )
    {
        gra_line_width( IsoSurface->LineWidth );
    } else {
        gra_line_width( 1.0 );
    }

    gra_set_material( IsoSurface->Material );

    if ( quick ) gra_polygon_mode( GRA_LINE );
    gra_begin( GRA_TRIANGLES );

    if ( IsoSurface->Recompute )
    {
        IsoSurface->NPolygons = 0;
        for( i=0; i<model->NofElements; i++ )
        {
            if ( !elements[i].DisplayFlag ) continue;

            if ( IsoSurface->MaxPolygons-IsoSurface->NPolygons < 1200 )
            {
                   IsoSurface->MaxPolygons += 1200;
                   IsoSurface->Polygons = (polygon_t *)realloc(
                              IsoSurface->Polygons,IsoSurface->MaxPolygons*sizeof(polygon_t)
                           );
            }

            n = vis_get_isosurface(
                              model,&elements[i],geometry->Vertices,C,F,NX,NY,NZ,
               &IsoSurface->Polygons[IsoSurface->NPolygons],IsoSurface->NofLevels,Levels,CScl,CAdd
                       );

            for( j=0; j<n; j++ ) vis_polygon( &IsoSurface->Polygons[IsoSurface->NPolygons+j] );

            IsoSurface->NPolygons += n;

            if ( epMouseDown )
            {
                if ( quick )
                {
                    if ( (i & 32) && (RealTime() - dt > TooLong2) )
                        if ( ++epMouseDownTakesTooLong > 3 )
                        {
                            gra_end();
                            gra_lon();
                            gra_polygon_mode( GRA_FILL );
                            return FALSE;
                        } else dt = RealTime();
               } else
               {
                   if ( RealTime() - dt > TooLong1 )
                   {
                       gra_end();
                       gra_lon();
  
                       epMouseDownTakesTooLong++;
                       return FALSE;
                   }
               }
            }
            if ( BreakLoop ) break;
        }
        IsoSurface->Recompute = FALSE;
    } else
    {
       for( j=0; j<IsoSurface->NPolygons; j++ )
       {
           vis_polygon( &IsoSurface->Polygons[j] );

           if ( (j & 64) && epMouseDown )
           {
               if ( quick )
               {
                   if ( RealTime() - dt > TooLong2 )
                       if ( ++epMouseDownTakesTooLong > 3 )
                       {
                           gra_end();
                           gra_lon();
                           gra_polygon_mode( GRA_FILL );
 
                           return FALSE;
                       } else dt = RealTime();
              } else
              {
                  if ( RealTime() - dt > TooLong1 )
                  {
                      gra_end();
                      gra_lon();
 
                      epMouseDownTakesTooLong++;
                      return FALSE;
                  }
              }
           }
           if ( BreakLoop ) break;
       }
    }

    gra_end();
    gra_lon();
    if ( quick ) gra_polygon_mode( GRA_FILL );

    return TRUE;
}


/*******************************************************************************
 *
 *     Name:        vis_isosurface_alloc
 *
 *     Purpose:     allocate memory for isosurface_t structure
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
static isosurface_t *vis_isosurface_alloc()
{
     isosurface_t *isosurface = (isosurface_t *)calloc(sizeof(isosurface_t),1);

     if ( !isosurface )
     {
         fprintf( stderr, "vis_isosurface_alloc: FATAL: can't alloc a few bytes of memory.\n" );
     }

     return isosurface;
}

/*******************************************************************************
 *
 *     Name:        vis_isosurface_delete
 *
 *     Purpose:     free memory associated with isosurface_t structure
 *
 *     Parameters: 
 *
 *         Input:   (isosurface_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_isosurface_delete(isosurface_t *isosurface)
{
    if ( isosurface ) {
       if ( isosurface -> MaxPolygons > 0 && isosurface -> Polygons ) free( isosurface->Polygons ); 
       free( isosurface );
    }
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_isosurface_visual
 *
 *     Purpose:     Register "Contour Line" visual type
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
int vis_initialize_isosurface_visual()
{
    static char *visual_name = "Isosurfaces";
    visual_type_t VisualDef; 

    static isosurface_t isosurface;

    static visual_param_t IsosurfaceParams[] =
    {
        { "Style",         "%d", 0, VIS_VISUAL_PARAM_INT, 0, 0.0, NULL },
        { "ContourData",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ColorData",     "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "NormalData0",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "NormalData1",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "NormalData2",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "NofLevels",     "%d", 0, VIS_VISUAL_PARAM_INT,     0, 0.0, NULL },
        { "Levels",        "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "Material",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { "ColorMap",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "LineWidth",     "%lf",0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "LineQuality",   "%d", 0, VIS_VISUAL_PARAM_INT,     1, 0.0, NULL },
        { "LineStyle",     "%d", 0, VIS_VISUAL_PARAM_INT,line_style_line, 0.0, NULL },
        { "Recompute",     "%d", 0, VIS_VISUAL_PARAM_INT,  TRUE, 0.0, NULL },
        { NULL, NULL, 0, 0, 0, 0.0, NULL }
    };

    int n = 0;

    IsosurfaceParams[n++].Offset = (char *)&isosurface.Style          - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.ContourData    - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.ColorData      - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.NormalData[0]  - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.NormalData[1]  - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.NormalData[2]  - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.NofLevels      - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.Levels         - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.Material       - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.ColorMap       - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.LineWidth      - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.LineQuality    - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.LineStyle      - (char *)&isosurface;
    IsosurfaceParams[n++].Offset = (char *)&isosurface.Recompute      - (char *)&isosurface;

    VisualDef.VisualName    = visual_name;

    VisualDef.RealizeVisual = (int   (*)()) vis_isosurface;
    VisualDef.AllocParams   = (void *(*)()) vis_isosurface_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_isosurface_delete;
    VisualDef.VisualParams  = IsosurfaceParams;

    return vis_add_visual_type( &VisualDef );
}
