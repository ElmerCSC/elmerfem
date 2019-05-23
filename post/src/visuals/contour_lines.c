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
 * Action routines for visual classe ContourLines.
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
 * 28 Sep 1995, modified vis_initialize_contour_lines_visual to set the
 *              VisualName field of the visual_type structure
 *
 * Juha R.
 *
 ******************************************************************************/

#include "../elmerpost.h"

/******************************************************************************
 * 
 * Parameter sructure definitios for Contour Lines visual class
 *
 ******************************************************************************/

typedef struct contour_lines_s
{
    scalar_t *ContourData;
    scalar_t *ColorData;

    int NofLevels;
    double *Levels;

    material_t *Material;
    colormap_t *ColorMap;

    int LineQuality;
    double LineWidth;

    line_style_t LineStyle;
} contour_lines_t;

/*******************************************************************************
 *
 *     Name:        vis_get_isolines
 *
 *     Purpose:     Extract isolines with given threshold 
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
static int vis_get_isolines
   (
    element_model_t *model, element_t *element, vertex_t *vertices, line_t *Lines,
      double *C, double *F, int nlevels, double *levels,double CScl,double CAdd
   )
{
    static double x[ELM_MAX_ELEMENT_NODES];
    static double y[ELM_MAX_ELEMENT_NODES];
    static double z[ELM_MAX_ELEMENT_NODES];
    static double f[ELM_MAX_ELEMENT_NODES];
    static double c[ELM_MAX_ELEMENT_NODES];

    int i,j,n,*T=element->Topology;
    element_type_t *elmt = element->ElementType;

    if ( elmt->IsoLine ) 
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
        }

        n = 0;
        for( i=0; i<nlevels; i++ )
        {
            n += (*elmt->IsoLine)( levels[i],f,c,x,y,z,&Lines[n] );
        } 

        return n;
    }

    return 0;
}


/*******************************************************************************
 *
 *     Name:        vis_draw_line
 *
 *     Purpose:     Draw  a line
 *
 *     Parameters: 
 *
 *         Input:   (line_t *)  
 *                  (int )     line/solid
 *                  (double)    width of line
 *
 *         Output:  graphics
 *   
 *   Return value:   void
 *
 ******************************************************************************/
static void vis_draw_line( line_t *line, int quick, double width )
{
    float x[2][3],c0,c1;

    x[0][0] = line->x[0]; 
    x[0][1] = line->y[0]; 
    x[0][2] = line->z[0]; 

    x[1][0] = line->x[1]; 
    x[1][1] = line->y[1]; 
    x[1][2] = line->z[1]; 

    c0 = line->c[0];
    c1 = line->c[1];

    if ( quick )
    {
        gra_line( x[0],c0,x[1],c1,line_style_line,width );
    } else
    {
        gra_line( x[0],c0,x[1],c1,line_style_cylinder,width );
        gra_sphere( x[0][0],x[0][1],x[0][2],c0, width );
        gra_sphere( x[1][0],x[1][1],x[1][2],c1, width );
    }
}

/*******************************************************************************
 *
 *     Name:        vis_contour_lines
 *
 *     Purpose:     Draw contour lines given data, and threshold values
 *
 *     Parameters: 
 *
 *         Input:   (geometry_t *)  
 *                  (contour_lines_t *) contour line display params
 *                  (double)            time used
 *
 *         Output:  graphics
 *   
 *   Return value:  if mouse interaction is going on, and time used exeeds
 *                  given value (TooLong1,2) FALSE, otherwise TRUE
 *
 ******************************************************************************/
static int vis_contour_lines( geometry_t *geometry, element_model_t *model, contour_lines_t *ContourLines,double dt )
{
    double CScl,CAdd,*Levels=ContourLines->Levels;
    double *C=NULL, *F=NULL;

    scalar_t *ColorData   = ContourLines->ColorData;
    scalar_t *ContourData = ContourLines->ContourData;

    int i,j,n,quick;

    double width = ContourLines->LineWidth*0.005;

    static line_t Lines[1000];

    element_t *elements = model->Elements;

    if ( !ContourData || !ContourData->f ) return TRUE;
    
    if ( !GlobalOptions.StereoMode )
      if ( ContourLines->Material->Diffuse[3]  < 1.0 )
      {
          if ( GlobalPass != 0 ) return TRUE;
      } else if ( GlobalPass == 0 )
      {
          return TRUE;
      }

    F = ContourData->f;

    if ( ColorData && ColorData->f )
    {
        CAdd = ColorData->min;
        CScl = 1.0/(ColorData->max - ColorData->min);

        C = ColorData->f;

        gra_set_colormap( ContourLines->ColorMap );
    } else gra_set_colormap( NULL );

    quick  = ContourLines->LineStyle == line_style_line;
    quick |= epMouseDown && epMouseDownTakesTooLong;

    if ( !quick && (ContourLines->LineStyle == line_style_cylinder) )
    {
        gra_sphere_quality( ContourLines->LineQuality );
    }

    if ( quick && !(epMouseDown && epMouseDownTakesTooLong) )
    {
        gra_line_width( ContourLines->LineWidth );
    } else {
        gra_line_width( 1.0 );
    }

    gra_set_material( ContourLines->Material );

    if ( quick ) gra_beg_lines();

    for( i=0; i<model->NofElements; i++ )
    {

        if ( !elements[i].DisplayFlag ) continue;

        n = vis_get_isolines
         (
           model,&elements[i],geometry->Vertices,Lines,C,F,ContourLines->NofLevels,Levels,CScl,CAdd
         );
        for( j=0; j<n; j++ ) vis_draw_line( &Lines[j], quick, width );

        if ( epMouseDown )
        {
            if ( quick )
            {
                if ( (i & 32) && (RealTime() - dt > TooLong2) )
                    if ( ++epMouseDownTakesTooLong > 3 )
                    {
                        gra_end_lines();
                        return FALSE;
                    } else dt = RealTime();
            } else
            {
                if ( RealTime() - dt > TooLong1 )
                {
                    epMouseDownTakesTooLong++;
                    return FALSE;
                }
            }
        }

        if ( BreakLoop ) break;
    }

    if ( quick ) gra_end_lines();

    return TRUE;
}


/*******************************************************************************
 *
 *     Name:        vis_contour_lines_alloc
 *
 *     Purpose:     allocate memory for contour_lines_t structure
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
static contour_lines_t *vis_contour_lines_alloc()
{
     contour_lines_t *contour_line = (contour_lines_t *)calloc(sizeof(contour_lines_t),1);

     if ( !contour_line )
     {
         fprintf( stderr, "vis_contour_lines_alloc: FATAL: can't alloc a few bytes of memory\n" );
     }

     return contour_line;
}

/*******************************************************************************
 *
 *     Name:        vis_contour_lines_delete
 *
 *     Purpose:     free memory associated with contour_lines_t structure
 *
 *     Parameters: 
 *
 *         Input:   (contour_lines_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_contour_lines_delete(contour_lines_t *contour_line)
{
    if ( contour_line ) free( contour_line );
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_contour_lines_visual
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
int vis_initialize_contour_line_visual()
{
    static char *visual_name = "Contour Lines";
    visual_type_t VisualDef; 

    static contour_lines_t contours;

    static visual_param_t ContourLineParams[] =
    {
        { "ContourData",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ColorData",     "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "NofLevels","%d", 0, VIS_VISUAL_PARAM_INT,     0, 0.0, NULL },
        { "Levels",        "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "Material",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { "ColorMap",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "LineWidth",     "%lf", 0, VIS_VISUAL_PARAM_FLOAT,   0, 1.0, NULL },
        { "LineQuality",   "%d", 0, VIS_VISUAL_PARAM_INT,     1, 0.0, NULL },
        { "LineStyle",     "%d", 0, VIS_VISUAL_PARAM_INT,line_style_line, 0.0, NULL },
        { NULL, NULL, 0, 0, 0, 0.0, NULL }
    };

    int n = 0;

    ContourLineParams[n++].Offset = (char *)&contours.ContourData    - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.ColorData      - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.NofLevels      - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.Levels         - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.Material       - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.ColorMap       - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.LineWidth      - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.LineQuality    - (char *)&contours;
    ContourLineParams[n++].Offset = (char *)&contours.LineStyle      - (char *)&contours;

    VisualDef.VisualName    = visual_name;

    VisualDef.RealizeVisual = (int   (*)()) vis_contour_lines;
    VisualDef.AllocParams   = (void *(*)()) vis_contour_lines_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_contour_lines_delete;
    VisualDef.VisualParams  = ContourLineParams;

    return vis_add_visual_type( &VisualDef );
}
