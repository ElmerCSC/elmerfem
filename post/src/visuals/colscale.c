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
 * Action routines for the colscale visual class.
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
 *                       Date: 15 Jan 1996
 *
 *
 * Modification history:
 *
 *
 ******************************************************************************/

#include "../elmerpost.h"


/******************************************************************************
 * 
 * Parameter sructure definitios for colscale visual class
 *
 ******************************************************************************/

typedef struct colscale_s
{
    double Length,Thickness,XPosition,YPosition,FontSize;
    int Entries,Decimals,Style,FontColor;

    scalar_t   *ColorData;
    colormap_t *ColorMap;
    material_t *Material;
} colscale_t;


/*******************************************************************************
 *
 *     Name:        vis_colscale
 *
 *     Purpose:     draw colscale as lines or surface, with color codes or not
 *
 *     Parameters: 
 *
 *         Input:   (geometry_t *) triangles to draw
 *                  (colscale_t     *) colscale display parameters
 *                  (double) 
 *
 *         Output:  graphics
 *   
 *   Return value:  if mouse interaction is going on, and time used exeeds
 *                  given value (TooLong1,2) FALSE, otherwise TRUE
 *
 ******************************************************************************/
static int vis_colscale( geometry_t *geometry, element_model_t *model,
                    colscale_t *ColorScale,double dt )
{
    scalar_t *ColorData = ColorScale->ColorData;

    static char str[120],fmt[120];
    int i,j,n;
 
    double CScl,CAdd,*C,y,xl;

    GLboolean clip[6], lights_on;

    float coords[4][3];

    if ( epMouseDown && epMouseDownTakesTooLong ) return TRUE;

#if 0
    if ( !(ColorData && ColorData->f) ) return TRUE;
#endif

    for( i=0; i<6; i++ )
    {
      glGetBooleanv( GL_CLIP_PLANE0+i, &clip[i] );
      if ( clip[i] ) glDisable( GL_CLIP_PLANE0+i );
    }

    glMatrixMode( GL_MODELVIEW );
    gra_push_matrix();
    gra_load_identity(); 

    glMatrixMode(GL_PROJECTION);
    gra_push_matrix();
    gra_load_identity();

    glGetBooleanv( GL_LIGHTING, &lights_on );
    glDisable( GL_LIGHTING );

    gra_polygon_mode( GRA_FILL );

    if ( !ColorData || !ColorData->f )
    {
       C = ColorData->f;

       CAdd = ColorData->min;
       if ( ABS(ColorData->max - ColorData->min)>0.0 )
         CScl = 1.0 / (ColorData->max - ColorData->min);
       else
         CScl = 1.0;
    } else {
       CAdd = 0.0;
       CScl = 1.0;
    }

    gra_set_colormap( ColorScale->ColorMap );
    gra_set_material( ColorScale->Material );

    n = ColorScale->Entries;
    gra_begin( GRA_QUADS );

    if ( ColorScale->Style ==  0 )
    {
      for( i=0; i<ColorScale->ColorMap->NumberOfEntries-1; i++  )
      {
         y = ColorScale->Length * i / (ColorScale->ColorMap->NumberOfEntries-1.0) + ColorScale->YPosition;
         coords[0][0] = ColorScale->XPosition;
         coords[1][0] = ColorScale->XPosition + ColorScale->Thickness;
         coords[2][0] = ColorScale->XPosition + ColorScale->Thickness;
         coords[3][0] = ColorScale->XPosition;

         coords[0][1] = y;
         coords[1][1] = y;
         coords[2][1] = y + ColorScale->Length/(ColorScale->ColorMap->NumberOfEntries-1.0);
         coords[3][1] = y + ColorScale->Length/(ColorScale->ColorMap->NumberOfEntries-1.0);

         coords[0][2] = 0.0;        
         coords[1][2] = 0.0;
         coords[2][2] = 0.0;
         coords[3][2] = 0.0;
   
         gra_flat_quad( coords, i / (ColorScale->ColorMap->NumberOfEntries-1.0) );
      }

      gra_end();

      glColor3f( ((ColorScale->FontColor & 0xff0000)>>16)/255.0,
                 ((ColorScale->FontColor & 0xff00)>>8)/255.0,
                    (ColorScale->FontColor & 0xff)/255.0 );

      sprintf( fmt, "%%.%de", ColorScale->Decimals );

      for( i=0; i<n; i++  )
      {
         sprintf( str, fmt, (ColorData->max - ColorData->min)*i/(n-1.0) + ColorData->min );
         y = ColorScale->Length * i / (n-1.0) + ColorScale->YPosition;

         xl = (ColorScale->Decimals + 9) * ColorScale->FontSize / (double)GraphicsXSize;
#ifndef WIN32
         if ( CurrentXFont ) {
           xl = (XTextWidth( CurrentXFont, str, strlen(str) ) + 125.0 ) / (double)GraphicsXSize;
         }
#endif
         glRasterPos3f( ColorScale->XPosition - xl,y,0.0 );
         PrintString( str );

         xl = ColorScale->FontSize / GraphicsXSize;
         gra_beg_lines();
           glVertex3f( ColorScale->XPosition-xl,y,0.01 );
           glVertex3f( ColorScale->XPosition,y,0.01 );
         gra_end_lines();
      }

      if ( ColorScale->ColorData->name ) {
         xl =  strlen(ColorScale->ColorData-name) * ColorScale->FontSize / (double)GraphicsXSize;
#ifndef WIN32
         if ( CurrentXFont ) {
           xl = (XTextWidth( CurrentXFont, str, strlen(ColorScale->ColorData->name)))/(double)GraphicsXSize;
         }
#endif
         glRasterPos3f( ColorScale->XPosition - xl,
             ColorScale->YPosition-ColorScale->Thickness,0.0 );

         PrintString( ColorScale->ColorData->name );
      }

    } else {

      for( i=0; i<ColorScale->ColorMap->NumberOfEntries; i++  )
      {
         y = ColorScale->Length * i / (ColorScale->ColorMap->NumberOfEntries-1.0) + ColorScale->XPosition;

         coords[0][0] = y;
         coords[1][0] = y + ColorScale->Length/(ColorScale->ColorMap->NumberOfEntries-1.0);
         coords[2][0] = y + ColorScale->Length/(ColorScale->ColorMap->NumberOfEntries-1.0);
         coords[3][0] = y;

         coords[0][1] = ColorScale->YPosition;
         coords[1][1] = ColorScale->YPosition;
         coords[2][1] = ColorScale->YPosition + ColorScale->Thickness;
         coords[3][1] = ColorScale->YPosition + ColorScale->Thickness;

         coords[0][2] = 0.0;        
         coords[1][2] = 0.0;
         coords[2][2] = 0.0;
         coords[3][2] = 0.0;
   
         gra_flat_quad( coords, i / (ColorScale->ColorMap->NumberOfEntries-1.0) );
      }

      gra_end();

      glColor3f( ((ColorScale->FontColor & 0xff0000)>>16)/255.0,
                 ((ColorScale->FontColor & 0xff00)>>8)/255.0,
                    (ColorScale->FontColor & 0xff)/255.0 );

      sprintf( fmt, "%%-8.%dg",ColorScale->Decimals );

      for( i=0; i<n; i++  )
      {
         sprintf( str, fmt, (ColorData->max - ColorData->min)*i/(n-1.0) + ColorData->min );

         y = ColorScale->Length * i / (n-1.0) + ColorScale->XPosition;
         glRasterPos3f( y,ColorScale->YPosition+ColorScale->Thickness+0.05,0.0 );
         PrintString( str );

         gra_beg_lines();
           glVertex3f( y,ColorScale->YPosition+ColorScale->Thickness,0.01 );
           glVertex3f( y,ColorScale->YPosition+ColorScale->Thickness+0.025,0.01 );
         gra_end_lines();
      }

      if ( ColorScale->ColorData->name ) {
         glRasterPos3f( ColorScale->XPosition,
             ColorScale->YPosition-ColorScale->Thickness,0.0 );

         PrintString( ColorScale->ColorData->name );
      }
    }

    if ( lights_on ) glEnable( GL_LIGHTING );
    gra_pop_matrix();

    glMatrixMode( GL_MODELVIEW );
    gra_pop_matrix();

    for( i=0; i<6; i++ )
      if ( clip[i] ) glEnable( GL_CLIP_PLANE0+i );

    return TRUE;
}



/*******************************************************************************
 *
 *     Name:        vis_colscale_alloc
 *
 *     Purpose:     allocate memory for colscale_t structure
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
static colscale_t *vis_colscale_alloc()
{
     colscale_t *colscale = (colscale_t *)calloc(sizeof(colscale_t),1);

     if ( !colscale )
     {
         fprintf( stderr, "vis_colscale_alloc: FATAL: can't alloc a few bytes of memory\n" );
     }

     return colscale;
}

/*******************************************************************************
 *
 *     Name:        vis_colscale_delete
 *
 *     Purpose:     free memory associated with colscale_t structure
 *
 *     Parameters: 
 *
 *         Input:   (colscale_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_colscale_delete(colscale_t *colscale)
{
    if ( colscale ) free( colscale );
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_colscale_visual
 *
 *     Purpose:     Register "ColorScale" visual type
 *
 *     Parameters: 
 *
 *         Input:   none
 *
 *         Output:  none
 *   
 *   Return value:  vis_add_visual_type (malloc success propably)...
 *
 ******************************************************************************/
int vis_initialize_colscale_visual()
{
    static char *visual_name = "ColorScale";
    visual_type_t VisualDef; 

    static colscale_t colscale;


    static visual_param_t ColorScaleParams[] =
    {
        { "Entries",       "%d", 0, VIS_VISUAL_PARAM_INT,  6, 0.0, NULL },
        { "Decimals",      "%d", 0, VIS_VISUAL_PARAM_INT,  2, 0.0, NULL },
        { "Style",         "%d", 0, VIS_VISUAL_PARAM_INT,  0,  0.0, NULL },
        { "XPosition",     "%d", 0, VIS_VISUAL_PARAM_FLOAT,  0,  0.8, NULL },
        { "YPosition",     "%d", 0, VIS_VISUAL_PARAM_FLOAT,  0, -0.8, NULL },
        { "Length",        "%d", 0, VIS_VISUAL_PARAM_FLOAT,  0,  1.5, NULL },
        { "Thickness",     "%d", 0, VIS_VISUAL_PARAM_FLOAT,  0,  0.1, NULL },
        { "Font Color",    "%d", 0, VIS_VISUAL_PARAM_INT,  0xffffff,  0.0, NULL },
        { "Font Size",     "%d", 0, VIS_VISUAL_PARAM_FLOAT, 0,  17.0, NULL },
        { "ColorData",     "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0,  NULL },
        { "ColorMap",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "Material",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { NULL, NULL, 0, 0, 0, 0.0, NULL }
    };

    int n = 0;

    ColorScaleParams[n++].Offset = (char *)&colscale.Entries     - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.Decimals    - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.Style       - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.XPosition   - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.YPosition   - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.Length      - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.Thickness   - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.FontColor   - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.FontSize    - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.ColorData   - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.ColorMap    - (char *)&colscale;
    ColorScaleParams[n++].Offset = (char *)&colscale.Material    - (char *)&colscale;

    VisualDef.VisualName    = visual_name;

    VisualDef.RealizeVisual = (int   (*)()) vis_colscale;
    VisualDef.AllocParams   = (void *(*)()) vis_colscale_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_colscale_delete;
    VisualDef.VisualParams  = ColorScaleParams;

    return vis_add_visual_type( &VisualDef );
}
