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
 * Action routines for arrow visual class.
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
 *
 * Modification history:
 *
 * 28 Sep 1995, modified vis_initialize_arrow_visual to set the VisualName field
 *              of the visual_type structure
 *
 * Juha R.
 *
 ******************************************************************************/

#include "../elmerpost.h"


/******************************************************************************
 * 
 * Parameter sructure definitios for arrow visual class
 *
 ******************************************************************************/
typedef struct arrow_s
{
    scalar_t *VectorData[4];           /* arrow direction quantity */
    scalar_t *ColorData;               /* arrow coloring quantity */
    scalar_t *LengthData;              /* arrow length quantity */
    scalar_t *ThresholdData;           /* arrow thresholding quantity */

    arrow_style_t Style;               /* arrow style, stick or arrow */

    double Floor,Ceiling;               /* threshold min & max */

    double RadiusScale,LengthScale;
    logical_t EqualLength,EqualRadius;

    material_t *Material;
    colormap_t *ColorMap;

    int LineQuality;
    line_style_t LineStyle;
} arrow_t;


/*******************************************************************************
 *
 *     Name:        vis_arrow
 *
 *     Purpose:     draw arrows given vector field
 *
 *     Parameters:
 *
 *         Input:   (geometry_t *) geometry description
 *                  (arrow_t    *) arrow display parameters
 *                  (double)       real time for interaction
 *
 *         Output:  graphics
 *
 *   Return value:  if mouse interaction is going on, and time used exeeds
 *                  given value (TooLong1,2) FALSE, otherwise TRUE
 *
 ******************************************************************************/
static int vis_arrow( geometry_t *geometry, element_model_t *model, arrow_t *Arrows,double dt )
{
    line_style_t line_style   = Arrows->LineStyle;
    arrow_style_t arrow_style = Arrows->Style;
  
    vertex_t *v = geometry->Vertices;

    int i,N=geometry->VertexCount,quick;

    float x[3];

    vertex_face_t *face;

    element_t *elements = model->Elements;

    double LScl,CScl,CAdd,vx,vy,vz,co;
    double R,RadiusScale,L,LengthScale;

    if ( !Arrows->VectorData ||
         !Arrows->VectorData[1]->f || !Arrows->VectorData[2]->f || !Arrows->VectorData[3]->f
        )
    {
/*
        fprintf( stderr, "vis_arrows: no arrow data\n" );
*/
        return TRUE;
    }

    if ( !GlobalOptions.StereoMode )
      if ( Arrows->Material->Diffuse[3]  < 1.0 )
      {
          if ( GlobalPass != 0 ) return TRUE;
      } else if ( GlobalPass == 0 )
      {
          return TRUE;
      }

    LengthScale = Arrows->LengthScale / 10.0;
    RadiusScale = Arrows->RadiusScale / 10.0;

    if ( Arrows->LengthData && Arrows->LengthData->f )
    {
        LScl = 0.05;
        if ( ABS( Arrows->LengthData->max - Arrows->LengthData->min ) > 1.0E-10 )
            LScl =  1.0 / (Arrows->LengthData->max - Arrows->LengthData->min);
    }

    if ( Arrows->ColorData && Arrows->ColorData->f )
    {
        CAdd = Arrows->ColorData->min;
        CScl = 1.0 / ( Arrows->ColorData->max - Arrows->ColorData->min );
        gra_set_colormap( Arrows->ColorMap );
    } else gra_set_colormap( NULL );

    quick = epMouseDown && epMouseDownTakesTooLong;

    if ( quick )
    {
      line_style  = line_style_line;
      arrow_style = arrow_style_stick;
      gra_line_width( 1.0 );
    } else
    {
      if ( line_style == line_style_cylinder )
        gra_sphere_quality( Arrows->LineQuality );
      else
        gra_line_width( Arrows->RadiusScale );
    }

    gra_set_material( Arrows->Material );

    for( i=0; i<N; i++,v++ )
    {
        if ( !v->ElementModelNode ) continue;

         for( face=v->Faces; face!=NULL; face=face->Next )
         {
           if ( geometry->Triangles[face->Face].Element->DisplayFlag ) break;
         }
#if 1
         if ( v->Faces && !face ) continue;
#else
         if ( !face ) continue;
#endif
        if ( Arrows->ThresholdData && Arrows->ThresholdData->f )
        {
            if ( Arrows->ThresholdData->f[i] < Arrows->Floor ||
                 Arrows->ThresholdData->f[i] > Arrows->Ceiling ) continue;
        }

        vx = Arrows->VectorData[1]->f[i];
        vy = Arrows->VectorData[2]->f[i];
        vz = Arrows->VectorData[3]->f[i];

        L  = sqrt( vx*vx + vy*vy + vz*vz );
        if ( L < 1.0E-10 ) continue;

        L = LengthScale / L;
        R = LengthScale * RadiusScale;

        if ( Arrows->ColorData && Arrows->ColorData->f )
        {
            co = CScl*(Arrows->ColorData->f[i]-CAdd);
        }

        x[0] = L*vx;
        x[1] = L*vy;
        x[2] = L*vz;

        if ( !Arrows->EqualLength && Arrows->LengthData && Arrows->LengthData->f )
        {
            x[0] *= LScl*Arrows->LengthData->f[i];
            x[1] *= LScl*Arrows->LengthData->f[i];
            x[2] *= LScl*Arrows->LengthData->f[i];

            if ( !Arrows->EqualRadius )
            {
                R *= LScl*Arrows->LengthData->f[i];
            }
        }

        gra_arrow( v->x,x,co,line_style,arrow_style,R );

        if ( epMouseDown && (i & 8) )
        {
            if ( !epMouseDownTakesTooLong )
            {
                if ( RealTime() - dt > TooLong1 )
                {
                    ++epMouseDownTakesTooLong;
                    return FALSE;
                }
            }
            else if ( RealTime() - dt > TooLong2 )
                if ( ++epMouseDownTakesTooLong > 3 )
                {
                    return FALSE;
                } else dt = RealTime();
        }

        if ( BreakLoop ) break;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        vis_arrow_alloc
 *
 *     Purpose:     allocate memory for arrow_t structure
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
static arrow_t *vis_arrow_alloc()
{
     arrow_t *arrow = (arrow_t *)calloc(sizeof(arrow_t),1);

     if ( !arrow )
     {
         fprintf( stderr, "vis_arrow_alloc: FATAL: can't alloc a few bytes of memory\n" );
     }

     return arrow;
}

/*******************************************************************************
 *
 *     Name:        vis_arrow_delete
 *
 *     Purpose:     free memory associated with arrow_t structure
 *
 *     Parameters: 
 *
 *         Input:   (arrow_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_arrow_delete(arrow_t *arrow)
{
    if ( arrow ) free( arrow );
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_arrow_visual
 *
 *     Purpose:     Register arrow visual type
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
int vis_initialize_arrow_visual()
{
    static char *visual_name = "Arrows";

    visual_type_t VisualDef; 

    static arrow_t arrow;

    static visual_param_t ArrowParams[] =
    {
        { "VectorData0",   "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "VectorData1",   "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "VectorData2",   "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "VectorData3",   "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ColorData",     "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "LengthData",    "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ThresholdData", "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "Style",         "%d",  0, VIS_VISUAL_PARAM_INT,     arrow_style_arrow, 0.0, NULL },
        { "Floor",         "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 0.0, NULL },
        { "Ceiling",       "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "RadiusScale",   "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "LengthScale",   "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "EqualLength",   "%c",  0, VIS_VISUAL_PARAM_LOGICAL, 1, 0.0, NULL },
        { "EqualRadius",   "%c",  0, VIS_VISUAL_PARAM_LOGICAL, 0, 0.0, NULL },
        { "Material",      "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { "ColorMap",      "%s",  0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "LineQuality",   "%d",  0, VIS_VISUAL_PARAM_INT,     1, 0.0, NULL },
        { "LineStyle",     "%d",  0, VIS_VISUAL_PARAM_INT,     line_style_line, 0.0, NULL },
        { NULL, NULL, 0, 0, 0,0.0, NULL }
    };

    int n = 0;

    ArrowParams[n++].Offset = (char *)&arrow.VectorData[0] - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.VectorData[1] - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.VectorData[2] - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.VectorData[3] - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.ColorData     - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.LengthData    - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.ThresholdData - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.Style         - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.Floor         - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.Ceiling       - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.RadiusScale   - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.LengthScale   - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.EqualLength   - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.EqualRadius   - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.Material      - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.ColorMap      - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.LineQuality   - (char *)&arrow;
    ArrowParams[n++].Offset = (char *)&arrow.LineStyle     - (char *)&arrow;

    VisualDef.VisualName    = visual_name;
    VisualDef.RealizeVisual = (int   (*)()) vis_arrow;
    VisualDef.AllocParams   = (void *(*)()) vis_arrow_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_arrow_delete;
    VisualDef.VisualParams  = ArrowParams;

    return vis_add_visual_type( &VisualDef );
}
