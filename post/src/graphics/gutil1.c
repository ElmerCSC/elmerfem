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
 * Graphics utilities part 1. OpenGL is called heavily by routines in this file.
 *
 *******************************************************************************
 *
 *                     Author: Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                             EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 26 Sep 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 * $Id: gutil1.c,v 1.4 2003/02/06 09:37:50 jpr Exp $
 *
 * $Log: gutil1.c,v $
 * Revision 1.4  2003/02/06 09:37:50  jpr
 * *** empty log message ***
 *
 * Revision 1.3  2001/07/02 11:21:44  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/07/31 13:36:56  jpr
 *
 * Added id, started log.
 *
 *
 ******************************************************************************/

#include "../elmerpost.h"

static colormap_t *CurrentColormap = NULL;

void vsetcolori( int ncolor )
{
    float r,g,b;

    if ( CurrentColormap )
    {
       r = CurrentColormap->Values[ncolor].r / 255.0;
       g = CurrentColormap->Values[ncolor].g / 255.0;
       b = CurrentColormap->Values[ncolor].b / 255.0;
       glColor4f( r,g,b,1.0 );
    } 
}

void vsetcolor( double color )
{
    float r,g,b;
    int ncolor;

    if ( CurrentColormap )
    {
       ncolor = (CurrentColormap->NumberOfEntries-1) * color + 0.5;
       r = CurrentColormap->Values[ncolor].r / 255.0;
       g = CurrentColormap->Values[ncolor].g / 255.0;
       b = CurrentColormap->Values[ncolor].b / 255.0;
       glColor4f( r,g,b,1.0 );
    } 
}

/*******************************************************************************
 *
 *     Name:        gra_set_material
 *
 *     Purpose:     set current material
 *
 *     Parameters: 
 *
 *         Input:   (material_t *) pointer to material
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void gra_set_material( material_t *Material )
{
    static material_t *CurrentMaterial = NULL;

    if ( !Material ) return;

/*
 *  if ( Material == CurrentMaterial && !Material->Changed ) return;
 */

    glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, Material->Shininess );
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, Material->Specular );

#if 0
    if ( Material->Specular[0] != 0.0 ||
         Material->Specular[1] != 0.0 ||
         Material->Specular[2] != 0.0
       ) glEnable( GL_NORMALIZE ); else glDisable( GL_NORMALIZE );
#endif
 
    glColor4fv( Material->Diffuse );

    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    if ( Material->Diffuse[3] < 1 )
    {
        glDisable( GL_DEPTH_TEST );
        glEnable( GL_BLEND );
    } else
    {
        glEnable( GL_DEPTH_TEST );
    }

    Material->Changed = FALSE;
    CurrentMaterial = Material;
}

/*******************************************************************************
 *
 *     Name:        gra_set_colormap
 *
 *     Purpose:     set current colormap
 *
 *     Parameters: 
 *
 *         Input:   (colormap_t *) pointer to colormap
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void gra_set_colormap( colormap_t *ColorMap )
{
    static rgb_t v = { 255,255,255 };
    int N;
    static colormap_t ctr;
    colormap_t *ct = ColorMap;

    if ( !ct || !ct->Values || ct->NumberOfEntries <= 1 )
    {
/*        { glDisable( GL_TEXTURE_1D ); return } */
          
          ct = &ctr;
          ct->Changed = TRUE;
          ct->Values = &v;
          ct->NumberOfEntries = 1;
    }

    glEnable( GL_TEXTURE_1D );

    if ( ct == CurrentColormap && !ct->Changed ) return;

    glGetIntegerv( GL_MAX_TEXTURE_SIZE, &N );

    if ( ct->NumberOfEntries > N )
    {
        fprintf( stderr, "gra_set_colormap: Colormap size [%d] too "
                         "big for this OpenGL-implementation\n", ct->NumberOfEntries );
        return;
    }

    glTexParameterf( GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexParameterf( GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameterf( GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP );

    glTexEnvf( GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE );

    glTexImage1D( GL_TEXTURE_1D, 0, 3, ct->NumberOfEntries, 0,
             GL_RGB, GL_UNSIGNED_BYTE, ct->Values );

    glEnable( GL_TEXTURE_1D );

    ct->Changed = FALSE;
    CurrentColormap = ct;
}

colormap_t *gra_get_colormap()
{
   return CurrentColormap;
}

/*******************************************************************************
 *
 *     Name:        gra_triangle
 *
 *     Purpose:     draw one triangle, glBegin() should have been called,
 *                  glEnd() should be called...
 *
 *     Parameters: 
 *
 *         Input:   (float [3][3]) vertex coordinates
 *                  (float [3][3]) normal vectors
 *                  (float [3])    vertex colors
 *
 *         Output:  graphics
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void gra_output_triangle( float coords[3][3],float normal[3][3],int ncolor )
{
    int i,j;
    double x0,y0,z0,x1,y1,z1,n0,n1,n2;


    x0 = coords[1][0] - coords[0][0];
    y0 = coords[1][1] - coords[0][1];
    z0 = coords[1][2] - coords[0][2];

    x1 = coords[2][0] - coords[0][0];
    y1 = coords[2][1] - coords[0][1];
    z1 = coords[2][2] - coords[0][2];

    n0 = y0*z1 - y1*z0;
    n1 = z0*x1 - z1*x0;
    n2 = x0*y1 - x1*y0;

    vsetcolori( ncolor );

    if ( n0*normal[0][0] + n1*normal[0][1] + n2*normal[0][2] > 0 ) {
       glNormal3fv( normal[0] );
       glVertex3fv( coords[0] );

       glNormal3fv( normal[1] );
       glVertex3fv( coords[1] );

       glNormal3fv( normal[2] );
       glVertex3fv( coords[2] );
    } else {
       glNormal3fv( normal[2] );
       glVertex3fv( coords[2] );

       glNormal3fv( normal[1] );
       glVertex3fv( coords[1] );

       glNormal3fv( normal[0] );
       glVertex3fv( coords[0] );
    }

    glNormal3f( 0.0,0.0,1.0 );
}

static void gra_internal_triangle( float coords[3][3],float normal[3][3],float color[3], int ncolor )
{
     float ncoord[3][3], nnorma[3][3], ccolor[3];
     double x[2], y[2], z[2], u[2], v[2], w[2], c[2];

     double c0=color[0],c1=color[1],c2=color[2],t,K;
     int i,n=0,k;

     int S, S0, S1, S2;

     if ( ncolor >= CurrentColormap->NumberOfEntries-1 ) 
     {
        gra_output_triangle( coords, normal, ncolor );
        return;
     }

     K = (double)(ncolor+1) / (double)CurrentColormap->NumberOfEntries;

     S0 = color[0] > K;
     S1 = color[1] > K;
     S2 = color[2] > K;
     S = S0 + S1 + S2;

     if ( S==0 )
     {
        gra_output_triangle( coords, normal, ncolor );
        return;
     }

     if ( S==3 )
     {
        gra_internal_triangle( coords, normal, color, ncolor+1 );
        return;
     }

     n = 0;
     if (S0 ^ S1)
     {
         t = (K - c0) / (c1 - c0);
 
         x[n] = t * ( coords[1][0] - coords[0][0] ) + coords[0][0];
         y[n] = t * ( coords[1][1] - coords[0][1] ) + coords[0][1];
         z[n] = t * ( coords[1][2] - coords[0][2] ) + coords[0][2];
         u[n] = t * ( normal[1][0] - normal[0][0] ) + normal[0][0];
         v[n] = t * ( normal[1][1] - normal[0][1] ) + normal[0][1];
         w[n] = t * ( normal[1][2] - normal[0][2] ) + normal[0][2];
         n++;
     }

     if (S0 ^ S2)
     {
         t = (K - c0) / (c2 - c0);
 
         x[n] = t * ( coords[2][0] - coords[0][0] ) + coords[0][0];
         y[n] = t * ( coords[2][1] - coords[0][1] ) + coords[0][1];
         z[n] = t * ( coords[2][2] - coords[0][2] ) + coords[0][2];
         u[n] = t * ( normal[2][0] - normal[0][0] ) + normal[0][0];
         v[n] = t * ( normal[2][1] - normal[0][1] ) + normal[0][1];
         w[n] = t * ( normal[2][2] - normal[0][2] ) + normal[0][2];
         n++;
     }

     if (S1 ^ S2)
     {
         t = (K - c1) / (c2 - c1);
 
         x[n] = t * ( coords[2][0] - coords[1][0] ) + coords[1][0];
         y[n] = t * ( coords[2][1] - coords[1][1] ) + coords[1][1];
         z[n] = t * ( coords[2][2] - coords[1][2] ) + coords[1][2];
         u[n] = t * ( normal[2][0] - normal[1][0] ) + normal[1][0];
         v[n] = t * ( normal[2][1] - normal[1][1] ) + normal[1][1];
         w[n] = t * ( normal[2][2] - normal[1][2] ) + normal[1][2];
         n++;
     }

     k = 3;
     if ( (S0 ^ S1) && (S0 ^ S2) ) k = 0;
     if ( (S0 ^ S1) && (S1 ^ S2) ) k = 1;
     if ( (S0 ^ S2) && (S1 ^ S2) ) k = 2;

     for( i=0; i<2; i++ )
     {
        ncoord[i][0] = x[i];
        ncoord[i][1] = y[i];
        ncoord[i][2] = z[i];
        nnorma[i][0] = u[i];
        nnorma[i][1] = v[i];
        nnorma[i][2] = w[i];
        ccolor[i] = K;
     }

     ncoord[2][0] = coords[k][0];
     ncoord[2][1] = coords[k][1];
     ncoord[2][2] = coords[k][2];
     nnorma[2][0] = normal[k][0];
     nnorma[2][1] = normal[k][1];
     nnorma[2][2] = normal[k][2];
     ccolor[2]    = color[k];
     if ( color[k] > K )
         gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );
     else
         gra_output_triangle( ncoord, nnorma, ncolor );

     switch( k ) 
     {
        case 2:
          for( i=0; i<3; i++ ) 
          {
             ncoord[2][i] = coords[0][i];
             nnorma[2][i] = normal[0][i];
          }
          ccolor[2] = color[0];
          if ( color[k] > K )
              gra_output_triangle( ncoord, nnorma, ncolor );
          else
              gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );

          for( i=0; i<3; i++ ) 
          {
             ncoord[0][i] = coords[1][i];
             nnorma[0][i] = normal[1][i];
          }
          ccolor[0] = color[1];
          if ( color[k] > K )
              gra_output_triangle( ncoord, nnorma, ncolor );
          else
              gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );
        break;

        case 1:
          for( i=0; i<3; i++ ) 
          {
             ncoord[2][i] = coords[0][i];
             nnorma[2][i] = normal[0][i];
          }
          ccolor[2] = color[0];
          if ( color[k] > K )
              gra_output_triangle( ncoord, nnorma, ncolor );
          else
              gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );

          for( i=0; i<3; i++ ) 
          {
             ncoord[0][i] = coords[2][i];
             nnorma[0][i] = normal[2][i];
          }
          ccolor[0] = color[2];
          if ( color[k] > K )
              gra_output_triangle( ncoord, nnorma, ncolor );
          else
              gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );
        break;

        case 0:
          for( i=0; i<3; i++ ) 
          {
             ncoord[2][i] = coords[1][i];
             nnorma[2][i] = normal[1][i];
          }
          ccolor[2] = color[1];
          if ( color[k] > K )
              gra_output_triangle( ncoord, nnorma, ncolor );
          else
              gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );

          for( i=0; i<3; i++ ) 
          {
             ncoord[0][i] = coords[2][i];
             nnorma[0][i] = normal[2][i];
          }
          ccolor[0] = color[2];
          if ( color[k] > K )
              gra_output_triangle( ncoord, nnorma, ncolor );
          else
              gra_internal_triangle( ncoord, nnorma, ccolor, ncolor+1 );
        break;
     }
}

void gra_triangle( float coords[3][3],float normal[3][3],float color[3] )
{
   int i;

   if ( GlobalOptions.OutputPS ) {
      gra_internal_triangle( coords, normal, color, 0 );
   } else {
      for( i=0; i<3; i++ )
      {
         glTexCoord1f( color[i] );
         glNormal3fv( normal[i] );
         glVertex3fv( coords[i] );
      }
   }
   glNormal3f( 0.0, 0.0, 1.0 );
}

/*******************************************************************************
 *
 *     Name:        gra_poly
 *
 *     Purpose:     draw one polygon
 *
 *     Parameters: 
 *
 *         Input:   int n
 *                  (float *) x,y,z vertices
 *                  (float *) u,v,w normals
 *                  (float *) c color
 *
 *         Output:  graphics
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void gra_poly3( int n,float *x,float *y,float *z,float *u,float *v,float *w,float *c)
{
   float coords[3][3], normal[3][3];
   int i;

   for( i=0; i<n; i++ )
   {
       normal[i][0] = u[i];
       normal[i][1] = v[i];
       normal[i][2] = w[i];

       coords[i][0] = x[i];
       coords[i][1] = y[i];
       coords[i][2] = z[i];
   }
   gra_triangle( coords,normal,c );
   glNormal3f( 0.0, 0.0, 0.0 );
}


/*******************************************************************************
 *
 *     Name:        gra_flat_quad
 *
 *     Purpose:     draw one flat quad, glBegin() should have been called,
 *                  glEnd() should be called...
 *
 *     Parameters: 
 *
 *         Input:   (float [4][3]) vertex coordinates
 *                  (float [1])    face color
 *
 *         Output:  graphics
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void gra_flat_quad( float coords[4][3],double color )
{
    int j;

    glNormal3f( 0.0,0.0,1.0 );

    if ( GlobalOptions.OutputPS )
       vsetcolor( color );
    else
       glTexCoord1f( color );

    glVertex3fv( coords[0] );
    glVertex3fv( coords[1] );
    glVertex3fv( coords[2] );
    glVertex3fv( coords[3] );
}

/*******************************************************************************
 *
 *     Name:         gra_line
 *
 *     Purpose:      draw a line with line style given, if style is
 *                   line_style_line glBegin() should have been called,
 *                   glEnd() should be called...
 *
 *     Parameters:
 *
 *         Input:    (float [3])    coordinates of the first endpoint
 *                   (double)       color of the first endpoint (0-1)
 *                   (float [3])    coordinates of the second endpoint
 *                   (double)       color of the second endpoint
 *                   (line_style_t) line_style_line or line_style_cylinder
 *                   (double)       line width, cylinder radius
 *
 *         Output:   graphics
 *
 *   Return value:   void
 *
 ******************************************************************************/
void gra_line(float *x0,double c0,float *x1,double c1,line_style_t style,double width)
{
    if ( style == line_style_line )
    {
        if ( GlobalOptions.OutputPS ) 
           vsetcolor( c0 );
        else
           glTexCoord1f( c0 );  
        glVertex3fv(  x0 );

        if ( GlobalOptions.OutputPS ) 
           vsetcolor( c1 );
        else
           glTexCoord1f( c1 );
        glVertex3fv(  x1 );
    } else
    {
        if ( GlobalOptions.OutputPS ) vsetcolor( c1 );
        gra_cylinder( x0[0],x0[1],x0[2],c0,x1[0],x1[1],x1[2],c1,width );
    }
}


/*******************************************************************************
 *
 *     Name:        gra_arrow
 *
 *     Purpose:     draw one arrow
 *
 *     Parameters: 
 *
 *         Input:   (float [3]) coordinate of the arrow base
 *                  (float [3]) vector components x,y,z
 *                  (double)    color (0-1)
 *                  (line_style_t)  line_style_line or line_style_cylinder
 *                  (arrow_style_t) arrow_style_stick or arrow_style_arrow
 *                  (double)         cylinder radius, line width
 *
 *         Output:  graphics
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void gra_arrow(float From[3],float Vector[3],double color,line_style_t line_style,arrow_style_t arrow_style,double width)
{
    double x  = Vector[0], y  = Vector[1], z  = Vector[2];
    double xo = From[0],   yo = From[1],   zo = From[2];

    double xa = xo + 0.65*x;
    double ya = yo + 0.65*y;
    double za = zo + 0.65*z;

    float r,g,b;
    int ncolor;

    x += xo;
    y += yo;
    z += zo;

    if ( GlobalOptions.OutputPS ) vsetcolor( color );

    if ( line_style == line_style_cylinder )
    {
        if ( arrow_style == arrow_style_arrow )
        {
            gra_cylinder( xo,yo,zo,color,xa,ya,za,color,width );
            gra_cone( xa,ya,za,color,2.0*width,x,y,z,color,0.0 );
        } else {
            gra_cylinder( xo,yo,zo,color,x,y,z,color,width );
        }
    } else
    {
          glDisable( GL_LIGHTING );
          glTexCoord1f( color ); 

         glBegin(GL_LINE_STRIP);
            glVertex3f( xo,yo,zo );
            glVertex3f( x,y,z );
         glEnd();

         if ( arrow_style == arrow_style_arrow )
         {
              double xp,yp,zp;

              xp = 0.75 * ( x - xo ) + xo;
              yp = 0.75 * ( y - yo ) + yo;
              zp = 0.75 * ( z - zo ) + zo;
              xa = xp + ( y - yp );
              ya = yp - ( x - xp );
              za = zp;

              glBegin(GL_LINE_STRIP);
                 glVertex3f(x,y,z);
                 glVertex3f(xa,ya,za);
                 glVertex3f(xp,yp,zp);
                 xa = xp - ( y - yp );
                 ya = yp + ( x - xp );
                 za = zp;
                 glVertex3f(xa,ya,za);
                 glVertex3f(x,y,z);
              glEnd();

              glBegin(GL_LINE_STRIP);
                 xa = xp + ( z - zp );
                 ya = yp;
                 za = zp - ( x - xp );
                 glVertex3f(x,y,z);
                 glVertex3f(xa,ya,za);
                 glVertex3f(xp,yp,zp);
                 xa = xp - ( y - yp );
                 ya = yp;
                 za = zp + ( x - xp );
                 glVertex3f(xa,ya,za);
                 glVertex3f(x,y,z);
              glEnd();
        }
        glEnable( GL_LIGHTING );
    }
}
