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
 *     Clipping routines for MATC graphics.
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
 *                       Date: 30 May 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/
/**************************************************************************
*
*   by Otto Pesonen
*
*   Started        08-SEP-88 / OP
*   Last edited    16-FEB-89 / OP
*
*   Version:       0.0
*   File   :       clip.c
*
*   Clip a polygon against any bounding box
*   An example program.
*
************************************o*************************************/


/*
 * $Id: clip.c,v 1.2 2005/05/27 12:26:19 vierinen Exp $ 
 *
 * $Log: clip.c,v $
 * Revision 1.2  2005/05/27 12:26:19  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:31  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

/**************************************************************************
   Clip a 2D-polygon against clipping box
   If all points are inside the clipping box the polygon is trivially
   accepted.
   Otherwise it is clipped in turn against single clipping edge.

   RETURN: Number of points in the polygon inside the box

   NOTICE! The new number of points may be GRATER than the in origin!
           The theoretical number of points newer execeeds twice 
           the number in origin. 
************************************o*************************************/
int clip_poly(n,x,y)
int *n;                            /* Number of points in polygon      */
    double *x,*y;                  /* Coordinate arrays of the polygon */
{
  int last,                        /* Last point code (in/out?)     */
      this;                        /* Current point code (in/out?)   */
  int i,j,k;                       /* Arrays index (orig, new & temp) */
  int nn;                          /* Internal counter of points       */
  int eg;                          /* Current edge for clipping         */
  int icode;                       /* Sum of points inside bounding box */
  double xx,yy;                    /* Current point coordinates, i:th   */
  double cx,cy;                    /* Coordinates of the clipped point  */ 
  double lx,ly;                    /* Coordinates of the previous point */
  double dx,dy;                    /* Length of the line segment        */

  nn = *n;
  icode = 0;
  last  = FALSE;

  for(i=0; i<nn ; i++)
    if (x[i]<=CL_XMAX && x[i]>=CL_XMIN &&
        y[i]<=CL_YMAX && y[i]>=CL_YMIN) icode++;
  if(icode == nn) return(nn);      /* All points inside the box! */

  for( eg=0 ; eg<4 ; eg++)         /* Loop over all the edges */
  {
    x[nn] = x[0];                  /* Handle the first point */
    y[nn] = y[0];                  /* this eases the program a lot! */
    nn++;

    for(i=j=0 ; i<nn ; i++)        /* Loop over the points! */
    {
      xx=x[i];
      yy=y[i];

      this = FALSE;
      switch(eg)                   /* Code the current point */
      { 
        case 0: if (yy <= CL_YMAX) this = TRUE; break;
        case 1: if (yy >= CL_YMIN) this = TRUE; break;
        case 2: if (xx <= CL_XMAX) this = TRUE; break;
        case 3: if (xx >= CL_XMIN) this = TRUE; break;
      }

      if(i>0)
      {
        if((this && !last)  ||  (last && !this))
        {
          dx = xx-lx;
          dy = yy-ly;
          switch(eg)               /* Cut against edge */
          {
            case 0:
              cx = lx + ((CL_YMAX-ly)*dx)/dy;
              cy = CL_YMAX;
              break;       
            case 1:
              cx = lx + ((CL_YMIN-ly)*dx)/dy;
              cy = CL_YMIN;
              break;       
            case 2:
              cy = ly + ((CL_XMAX-lx)*dy)/dx;
              cx = CL_XMAX;
              break;       
            case 3:
              cy = ly + ((CL_XMIN-lx)*dy)/dx;
              cx = CL_XMIN;
              break;       
          }
        }

        if(last)                   /* Decide to store the point(s) */
          if(this) 
            { x[j]=xx; y[j]=yy; j++; }
          else
            { x[j]=cx; y[j]=cy; j++; }
        else
          if(this) 
          {
            if(j+2 > i)            /* Too many points whitout shift */
            {
              for(k=nn; k>i ; k--) { x[k]=x[k-1]; y[k]=y[k-1]; }
              nn++;
              i++;                 /* Update pointer AND the limit  */
            } 
            x[j]=cx; y[j]=cy; j++; /* Store two points */
            x[j]=xx; y[j]=yy; j++;
          }
      }
      else
        if(this)                   /* Store the first point */
          { x[j]=xx; y[j]=yy; j++; }

      lx = xx;
      ly = yy;
      last = this;                 /* Save the last point */

    }
    *n = j;
    if(!j) return(j);              /* No need to proceed */
    nn = j;
  }
  *n = j;
  return(j);
}

#define NULL_EDGE   0
#define LEFT_EDGE   1
#define RIGHT_EDGE  2
#define BOTTOM_EDGE 4
#define TOP_EDGE    8

void clip_code(x,y,c) double x,y; int *c;
{
  *c = NULL_EDGE;

  if (x < CL_XMIN)
    *c = LEFT_EDGE;
  else if (x > CL_XMAX)
    *c = RIGHT_EDGE;
  
  if (y < CL_YMIN)
    *c |= BOTTOM_EDGE;
  else if (y > CL_YMAX)
    *c |= TOP_EDGE;
}

/*
 * Polyline clipping routine. Return value is number of points inside
 * clipping limits when first linesegment not totally inside
 * the limits is found (the last one is included in the count if it is
 * part of a visible line segment). One should call this routine repeatedly
 * until all the linesegments in the polyline are handled. Edge crossing
 * points of the last linesegment are stored in place of originals.
 */
int clip_line(n,x,y)
int *n;
double *x, *y;
{
  double xx, yy,
         px, py;

  int i,c,c1,c2;

  px = x[0]; 
  py = y[0];
  clip_code(px,py,&c1); 
  for(i = 1; i < *n; i++)
  {
    clip_code(x[i],y[i],&c2);
    if (c1 || c2)
    {
      while(c1 || c2)
      {
        if (c1 & c2)
        {
          *n = i;
          return *n;
        }

        c = c1 ? c1 : c2;

        if (c & LEFT_EDGE)
        {
          yy = py+(y[i]-py)*(CL_XMIN-px)/(x[i]-px);
          xx = CL_XMIN;
        }
        else if (c & RIGHT_EDGE)
        {
          yy = py+(y[i]-py)*(CL_XMAX-px)/(x[i]-px);
          xx = CL_XMAX;
        }
        else if (c & BOTTOM_EDGE)
        {
          xx = px+(x[i]-px)*(CL_YMIN-py)/(y[i]-py);
          yy = CL_YMIN;
        }
        else if (c & TOP_EDGE)
        {
          xx = px+(x[i]-px)*(CL_YMAX-py)/(y[i]-py);
          yy = CL_YMAX;
        }
     
        if (c == c1)
        {
          x[i-1] = px = xx; 
          y[i-1] = py = yy; 
          clip_code(xx,yy,&c1);
        }
        else
        {
          x[i] = xx; 
          y[i] = yy; 
          clip_code(xx,yy,&c2);
        }
      }
      *n = i + 1;
      return *n; 
    }
    px = x[i]; 
    py = y[i];
    c1 = c2;
  }
  return *n;
}
