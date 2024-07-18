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
 * Graphics utilities Part II. OpenGL is being called.
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
 * $Id: gutil2.c,v 1.3 2003/02/06 09:37:50 jpr Exp $
 *
 * $Log: gutil2.c,v $
 * Revision 1.3  2003/02/06 09:37:50  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/07/31 13:36:58  jpr
 *
 * Added id, started log.
 *
 *
 ******************************************************************************/
#include "../elmerpost.h"

static int GRA_SphereQuality=3,GRA_SphereInitDone=FALSE;

#define GRA_MAX_QUALITY 32

static double GRA_CosA[4*GRA_MAX_QUALITY];
static double GRA_SinA[4*GRA_MAX_QUALITY];
static double GRA_CosB[4*GRA_MAX_QUALITY];
static double GRA_SinB[4*GRA_MAX_QUALITY];

/*******************************************************************************
 *
 *     Name:          gra_sphere_quality
 *
 *     Purpose:       initialize internal cos & sin tables with given resolution
 *
 *     Parameters:
 *
 *         Input:     (int)  quality factor
 *
 *         Output:    none
 *
 *   Return value:    void
 *
 ******************************************************************************/
void gra_sphere_quality(int quality)
{
    int i,j,Division;
    double a;

    if ( GRA_SphereInitDone )
        if ( GRA_SphereQuality == quality ) return;

    if ( quality > GRA_MAX_QUALITY ) quality = GRA_MAX_QUALITY;
    else if ( quality < 1 ) quality = 1;

    GRA_SphereQuality = quality;

    Division = 4*quality;

    for( i=0; i<Division; i++ )
    {
        a = 2.0*M_PI*i/(double)Division;
        GRA_CosA[i] = cos(a);
        GRA_SinA[i] = sin(a);
    }
    GRA_CosA[i] = GRA_CosA[0];
    GRA_SinA[i] = GRA_SinA[0];

    for( i=-Division/2,j=0; i<=Division/2; i+=2, j++ )
    {
        a = M_PI*i/(double)Division;
        GRA_CosB[j] = cos(a);
        GRA_SinB[j] = sin(a);
    }

    GRA_SphereInitDone = TRUE;
}

/*******************************************************************************
 *
 *     Name:          gra_vector_to_angles
 *
 *     Purpose:       return angles about x and y axis to rotate x-axis
 *                    to given vector
 *
 *     Parameters:
 *
 *         Input:     (double vx,vy,vz)     first end point
 *                    (double *px,*py,*pz)  second end point
 *
 *         Output:    (double *px,*py)      computed angles
 *
 *   Return value:    void
 *
 ******************************************************************************/
void gra_vector_to_angles(double vx,double vy,double vz,double *px,double *py,double *pz,double r)
{
     double a,b,x,y,z,RadToDeg=180.0/M_PI;

     x = (*px-vx)*r;
     y = (*py-vy)*r;
     z = (vz-*pz)*r;

     b = atan2(y,z);

     r = y*sin(b) + z*cos(b);
     a = atan2(r,x);

     *px = RadToDeg*b;
     *py = RadToDeg*a;
     *pz = 0.0;
}

/*******************************************************************************
 *
 *     Name:          gra_sphere
 *
 *     Purpose:       draw sphere centered at a point, with given color and
 *                    radius
 *
 *     Parameters:
 *
 *         Input:     (double x,y,z)  sphere center point
 *                    (double f)      color
 *                    (double r)      radius
 *
 *         Output:    graphics
 *
 *   Return value:    void
 *
 ******************************************************************************/
void gra_sphere(double x,double y,double z,double f,double r)
{
    int Division = 4*GRA_SphereQuality;
    int i,j;

    if ( !GRA_SphereInitDone ) gra_sphere_quality( GRA_SphereQuality );

    glPushMatrix();

    glTranslatef(x,y,z);
    glScalef(r,r,r);

    for( i=0; i<Division/2; i++ )
    {
        glBegin(GL_QUAD_STRIP);
            x = GRA_CosB[i];
            y = 0.0;
            z = GRA_SinB[i];

            glTexCoord1d( f );
            glNormal3f(-x,-y,-z);
            glVertex3f(x,y,z);

            x = GRA_CosB[i+1];
            y = 0.0;
            z = GRA_SinB[i+1];

            glTexCoord1d( f );
            glNormal3d( -x,-y,-z );
            glVertex3d( x,y,z );

            for( j=1; j<=Division; j++ )
            {
                x = GRA_CosA[j]*GRA_CosB[i];
                y = GRA_SinA[j]*GRA_CosB[i];
                z = GRA_SinB[i];

                glTexCoord1d( f );
                glNormal3d( -x,-y,-z );
                glVertex3d( x,y,z );

                x = GRA_CosA[j]*GRA_CosB[i+1];
                y = GRA_SinA[j]*GRA_CosB[i+1];
                z = GRA_SinB[i+1];

                glTexCoord1d( f );
                glNormal3d( -x,-y,-z );
                glVertex3d( x,y,z );
            }
        glEnd();
    }

    glPopMatrix();

    glNormal3d( 0.0,0.0,1.0 );
}

/*******************************************************************************
 *
 *     Name:          gra_point
 *
 *     Purpose:       draw point  with given color and  radius
 *
 *     Parameters:
 *
 *         Input:     (double x,y,z)  sphere center point
 *                    (double f)      color
 *                    (double r)      radius
 *
 *         Output:    graphics
 *
 *   Return value:    void
 *
 ******************************************************************************/
void gra_point(double x,double y,double z,double f,double r)
{
    glTexCoord1d( f );
    glVertex3d( x,y,z ); 
}

/*******************************************************************************
 *
 *     Name:          gra_cylinder
 *
 *     Purpose:       draw cylinder between points given
 *
 *     Parameters:
 *
 *         Input:     (double x0,y0,z0)  first end point coordinates
 *                    (double f0)        first end point color
 *                    (double x1,y1,z1)) second end point coordinates
 *                    (double f1)        second end point color
 *                    (double R)         cylinder radius
 *
 *         Output:    graphics
 *
 *   Return value:    void
 *
 ******************************************************************************/
void gra_cylinder(double x0,double y0,double z0,double f0,double x1,double y1,double z1,double f1,double R)
{
    double ax,ay,az,s;

    int i,j;

    if ( !GRA_SphereInitDone ) gra_sphere_quality( GRA_SphereQuality );

    ax = x1;
    ay = y1;
    az = z1;
    x1 = ax - x0;
    y1 = ay - y0;
    z1 = az - z0;
    s = sqrt( x1*x1+y1*y1+z1*z1 );
    if ( s < 1.0e-10 ) return;

    gra_vector_to_angles( x0,y0,z0,&ax,&ay,&az,1/s );

    glPushMatrix();

    glTranslated( x0,y0,z0 );

    glRotated( ax,1.0,0.0,0.0 );
    glRotated( ay,0.0,1.0,0.0 );

    glScaled( s,R,R );

    glBegin(GL_QUAD_STRIP);

        glTexCoord1d( f0 );
        glNormal3d( 0.0,0.0,1.0 );
        glVertex3d( 0.0,0.0,1.0 );

        glTexCoord1d( f1 );
        glVertex3d( 1.0,0.0,1.0 );

        for( j=1; j<=GRA_SphereQuality*4;j++ )
        {
            y0 = GRA_SinA[j];
            z0 = GRA_CosA[j];

            glTexCoord1d( f0 );
            glNormal3d( 0.0,y0,z0 );
            glVertex3d( 0.0,y0,z0 );

            glTexCoord1d( f1 );
            glVertex3d( 1.0,y0,z0 );
        }

    glEnd();

    glPopMatrix();

    glNormal3d( 0.0,0.0,1.0 );
}

/*******************************************************************************
 *
 *     Name:          gra_cone
 *
 *     Purpose:       draw a cone between points given
 *
 *     Parameters:
 *
 *         Input:     (double x0,y0,z0)  first end point coordinates
 *                    (double f0)        first end point color
 *                    (double R0)        first end point radius
 *                    (double x1,y1,z1)) second end point coordinates
 *                    (double f1)        second end point color
 *                    (double R1)        second end point radius
 *
 *         Output:    graphics
 *
 *   Return value:    void
 *
 ******************************************************************************/
void gra_cone(double x0,double y0,double z0,double f0,double R0,double x1,double y1,double z1,double f1,double R1)
{
    double ax,ay,az,r,s;
    int i,j;

    if ( !GRA_SphereInitDone ) gra_sphere_quality(GRA_SphereQuality);

    ax = x1;
    ay = y1;
    az = z1;

    x1 = ax - x0;
    y1 = ay - y0;
    z1 = az - z0;
    r = sqrt( x1*x1+y1*y1+z1*z1 );
    if ( r < 1.0E-10 ) return;

    gra_vector_to_angles( x0,y0,z0,&ax,&ay,&az,1/r );

    glPushMatrix();

    glTranslated(x0,y0,z0);

    glRotated( ax,1.0,0.0,0.0 );
    glRotated( ay,0.0,1.0,0.0 );

    glScaled( r,1.0,1.0 );
    
    s = sqrt(1+(R0-R1)*(R0-R1));
    s = 1/s;

#define NORMAL(x,y,z) glNormal3d(s*(x),s*(y),s*(z))
    
    glBegin(GL_QUAD_STRIP);

        NORMAL( R0-R1,0.0,1.0 );
        glTexCoord1d( f0 );
        glVertex3d( 0.0,0.0,R0 );

        NORMAL( R0-R1,0.0,1.0 );
        glTexCoord1d( f1 );
        glVertex3d( 1.0,0.0,R1 );

        for( j=1; j<=GRA_SphereQuality*4;j++ )
        {
            y0 = GRA_SinA[j];
            z0 = GRA_CosA[j];

            NORMAL( R0-R1,y0,z0 );

            glTexCoord1d( f0 );
            glVertex3d( 0.0,R0*y0,R0*z0 );

            glTexCoord1d( f1 );
            glVertex3d( 1.0,R1*y0,R1*z0 );
        }

    glEnd();

#undef NORMAL

    glPopMatrix();

    glNormal3d( 0.0,0.0,1.0 );
}
