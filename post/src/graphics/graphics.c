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
 * Graphics main module & some utilities
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
 *                       Date: 27 Sep 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 * 
 * $Id: graphics.c,v 1.4 2003/02/06 09:37:50 jpr Exp $
 *
 * $Log: graphics.c,v $
 * Revision 1.4  2003/02/06 09:37:50  jpr
 * *** empty log message ***
 *
 * Revision 1.3  2001/07/06 11:55:24  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/07/31 13:36:54  jpr
 *
 * Added id, started log.
 *
 *
 ******************************************************************************/

#define MODULE_GRAPHICS

#include "../elmerpost.h"

void gra_clip_plane( int plane, double equation[] )
{
    glClipPlane( (int)GL_CLIP_PLANE0+plane, equation );
    glEnable( (int)GL_CLIP_PLANE0+plane );
}

void gra_disable_clip( int plane )
{
    glDisable( (int)GL_CLIP_PLANE0+plane );
}

void gra_init()
{
    GLfloat pos[4] = { 0.0,0.0,1.0,0.0 };

    GLfloat amb[4] = { 0.2,0.2,0.2,1.0 };
    GLfloat spc[4] = { 1.0,1.0,1.0,1.0 };
    GLfloat dif[4] = { 0.8,0.8,0.8,1.0 };

    glPixelStorei( GL_PACK_ALIGNMENT,   1 );
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    gluPerspective( 30.0,1.0,1.0,20.0 );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glColorMaterial( GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE );
    glEnable( GL_COLOR_MATERIAL );

    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE,1 );
    glLightModelf( GL_LIGHT_MODEL_LOCAL_VIEWER,1.0 );
    glEnable( GL_LIGHTING );

    glLightfv( GL_LIGHT0, GL_POSITION, pos );
    glLightfv( GL_LIGHT0, GL_DIFFUSE,  dif );
    glLightfv( GL_LIGHT0, GL_SPECULAR, spc );
    glLightfv( GL_LIGHT0, GL_AMBIENT,  amb );

    glEnable( GL_LIGHT0 );

    gluLookAt( 0.0,0.0,5.0,0.0,0.0,0.0,0.0,1.0,0.0 );

    glEnable( GL_NORMALIZE );
    glEnable( GL_DEPTH_TEST );

    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glDisable( GL_BLEND );

    glCullFace( GL_BACK );
    glDisable( GL_CULL_FACE ); 
    glDepthRange( 0.0,1.0 );
    glDepthMask( GL_TRUE );
    glDepthFunc( GL_LEQUAL );

    // Enable anti-aliasing:
#if 0
    glEnable( GL_BLEND );
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
#endif
}

void gra_draw_viewport_frame()
{
    float Color[4];

    GLboolean DepthEnabled;

    glLineWidth( 1.0 );

    glGetBooleanv( GL_DEPTH_TEST,&DepthEnabled );
    glDisable( GL_DEPTH_TEST );

    glPushMatrix();
    glLoadIdentity(); 

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity(); 

    glGetFloatv( GL_CURRENT_COLOR,Color );

    glColor4f( 0.0,1.0,0.0,1.0 );

    glBegin( GL_LINE_LOOP );
    glVertex2f(-1.0,-1.0 );
    glVertex2f( 1.0,-1.0 );
    glVertex2f( 1.0, 1.0 );
    glVertex2f(-1.0, 1.0 );
    glEnd();

    glColor4fv( Color );
 
    glPopMatrix();

    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
    if ( DepthEnabled ) glEnable( GL_DEPTH_TEST );
}

void gra_set_projection( camera_proj_t Proj,double FieldAngle,
	double lx,double hx,double ly,double hy,double Near,double Far, int frame )
{
    int ilx     = (GraphicsXSize-1)*lx + 0.5;
    int iwidth  = (GraphicsXSize-1)*(hx-lx) + 0.5;

    int ily     = (GraphicsYSize-1)*ly + 0.5;
    int iheight = (GraphicsYSize-1)*(hy-ly) + 0.5;

    double Aspect = (double)iwidth/(double)iheight;

    glViewport( ilx,ily,iwidth,iheight );

    if ( frame ) gra_draw_viewport_frame();

    glMatrixMode( GL_PROJECTION );

    glLoadIdentity();

    if ( Proj == camera_proj_ortho )
        glOrtho( -Aspect,Aspect,-1.0,1.0,Near,Far );
    else
        gluPerspective( FieldAngle,Aspect,Near,Far );

    glMatrixMode( GL_MODELVIEW );
}

void gra_bbox( double XMin,double XMax,double YMin,double YMax,double ZMin,double ZMax )
{
    GLboolean TexEnabled,DepthEnabled;
    void gra_set_colormap( colormap_t * );
    colormap_t *gra_get_colormap(), *SaveColorMap;

#if 0
    glGetBooleanv(GL_TEXTURE_1D,&TexEnabled);
    if (  TexEnabled )  glDisable(GL_TEXTURE_1D);
#else
    SaveColorMap = gra_get_colormap();
    gra_set_colormap( NULL );
#endif

#if 0
    glGetBooleanv(GL_DEPTH_TEST,&DepthEnabled);
    if ( DepthEnabled ) glDisable(GL_DEPTH_TEST);
#endif
    glLineWidth( 1.0 );

    glDisable( GL_LIGHTING );
    glBegin(GL_LINE_LOOP);
    glVertex3f( XMin, YMin, ZMin );
    glVertex3f( XMax, YMin, ZMin );
    glVertex3f( XMax, YMax, ZMin );
    glVertex3f( XMin, YMax, ZMin );
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3f( XMin, YMin, ZMax );
    glVertex3f( XMax, YMin, ZMax );
    glVertex3f( XMax, YMax, ZMax );
    glVertex3f( XMin, YMax, ZMax );
    glEnd();

    glBegin(GL_LINES);
    glVertex3f( XMin, YMin, ZMin );
    glVertex3f( XMin, YMin, ZMax );

    glVertex3f( XMax, YMin, ZMin );
    glVertex3f( XMax, YMin, ZMax );

    glVertex3f( XMax, YMax, ZMin );
    glVertex3f( XMax, YMax, ZMax );

    glVertex3f( XMin, YMax, ZMin );
    glVertex3f( XMin, YMax, ZMax );
    glEnd();

    glBegin(GL_LINES);
    glVertex3f( 0.0,0.0,0.0 );
    glVertex3f( 0.5,0.0,0.0 );

    glVertex3f( 0.6,-0.1,0.0 );
    glVertex3f( 0.8, 0.1,0.0 );

    glVertex3f( 0.6, 0.1,0.0 );
    glVertex3f( 0.8,-0.1,0.0 );
    glEnd();

    glBegin(GL_LINES);
    glVertex3f( 0.0,0.0,0.0 );
    glVertex3f( 0.0,0.5,0.0 );

    glVertex3f(-0.1, 0.7,0.0 );
    glVertex3f( 0.0, 0.6,0.0 );

    glVertex3f( 0.0, 0.6,0.0 );
    glVertex3f( 0.1, 0.7,0.0 );

    glVertex3f( 0.0, 0.6,0.0 );
    glVertex3f( 0.0, 0.5,0.0 );
    glEnd();

    glBegin(GL_LINES);
    glVertex3f( 0.0,0.0,0.0 );
    glVertex3f( 0.0,0.0,0.5 );

    glVertex3f( -0.1, 0.1, 0.6 );
    glVertex3f(  0.1, 0.1, 0.6 );

    glVertex3f(  0.1, 0.1, 0.6 );
    glVertex3f( -0.1,-0.1, 0.6 );

    glVertex3f( -0.1,-0.1, 0.6 );
    glVertex3f(  0.1,-0.1, 0.6 );
    glEnd();
    glEnable( GL_LIGHTING );

#if 0
    if ( TexEnabled )   glEnable( GL_TEXTURE_1D );
#else
    gra_set_colormap( SaveColorMap );
#endif

#if 0
    if ( DepthEnabled ) glEnable( GL_DEPTH_TEST );
#endif
}

