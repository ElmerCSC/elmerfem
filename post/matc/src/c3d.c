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
 *     MATC surface display routines.
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
/***********************************************************************
|
|  Written      ??-???-88 / JPR
|  Last edited  16-FEB-89 / OP
|
*  File: c3d
^
|  USING THE C3D (preliminary) 
|
|  Full usage will be included later!
|
^**********************************************************************/


/*
 * $Id: c3d.c,v 1.2 2005/05/27 12:26:14 vierinen Exp $ 
 *
 * $Log: c3d.c,v $
 * Revision 1.2  2005/05/27 12:26:14  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:30  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

#define C3D_MASK 9
#define C3D_HALFMASK (1 << (C3D_MASK - 1))

/***********************************************************************
  Typedefinitions for panel tree
***********************************************************************/
struct node
{
  int x, y, z, d;
};
 
struct element
{
  struct node *node[4];
  int d, z; 
};

struct el_tree
{
  struct el_tree *left, *right;
  struct element *entry;
};

void C3D_Contour(double *, int, int);
void C3D_Levels(int); 

void C3D_Show_Tri(int *, int *, int *d);
void C3D_Show_Elem(struct element *el);
void C3D_SelCol(int);

void C3D_Show_El_Tree(struct el_tree *head);
void C3D_Add_El_Tree(struct el_tree *head, struct el_tree *add);

void C3D_Pcalc(int, int, int, int, int, int, int *, int *, int *,int  *);
void C3D_AreaFill(int, int, int *, int *);

VARIABLE *c3d_gc3d(var) VARIABLE *var;
{
  C3D_Contour(MATR(var), NROW(var), NCOL(var));
  return (VARIABLE *)NULL;
}  

VARIABLE *c3d_gc3dlevels(var) VARIABLE *var;
{
  C3D_Levels((int)*MATR(var));
  return (VARIABLE *)NULL;
}

/***********************************************************************
      Globals needed for C3D All start with the prefix c3d_
***********************************************************************/
static int c3d_clevels     = 10,
    c3d_perspective = FALSE;

/***********************************************************************
void C3D_Contour(matrix, nr, nc) double *matrix; int nr, nc;
?  Main function to be used.
|  This function displays a matrix seen from the current C3D-viewpoint
|  and with number of levels selected.
|
|  NOTICE! CALL FROM FORTRAN:
|
|    CALL C3D_Contour( a , %val( DIMY ) , %val( DIMX ) )
|
|    The a should be of type double precision (or REAL*8 in VAX/VMS)
|    and the dimensions are in reverse order
!  No error messages are generated
&  C3D...
***********************************************************************/
void C3D_Contour(double *matrix, int nr, int nc)
{
  double xmin, xmax, ymin, ymax, dmin, dmax;
  double x, y, z, d, xs, ys, zs;

  struct el_tree *tr, *trb, *element_el_tree = NULL;
  struct element *el, *elements = NULL;
  struct node *nodes = NULL;

  GMATRIX gm;

  static GMATRIX ident = 
         { 
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
         };

  int i, j, k, n;

  nodes = (struct node *)ALLOCMEM(nr * nc * sizeof(struct node));

  xmin = ymin = dmin =  1.0e20;
  xmax = ymax = dmax = -1.0e20;

  for(i = 0, n = 0; i < nr; i++) 
    for(j = 0; j < nc; j++, n++) {
      z = matrix[n]; dmin = min(dmin, z);  dmax = max(dmax, z);
    }

  for(i = 0, n = 0; i < nr; i++) 
    for(j = 0; j < nc; j++, n++)
    {
      x = 2 * (double)i / (double)nr - 1;
      y = 2 * (double)j / (double)nc - 1;
      d = (matrix[n] - dmin) / (dmax - dmin);
      z = 2 * d - 1;
      gra_mtrans(x, y, z, &xs, &ys, &zs);
      xs *= (1 << 20);
      ys *= (1 << 20);
      zs *= (1 << 20);
      nodes[n].x = xs;
      nodes[n].y = ys;
      nodes[n].z = zs;
      nodes[n].d = (c3d_clevels * d + 1) * (1 << C3D_MASK);
      xmin = min(xmin, xs);
      xmax = max(xmax, xs);
      ymin = min(ymin, ys);
      ymax = max(ymax, ys);
    }

  for(i = 0, n = 0; i < nr; i++) 
    for(j = 0; j < nc; j++, n++)
    {
      nodes[n].x = 4095 * (nodes[n].x-xmin) / (xmax-xmin); 
      nodes[n].y = 4095 * (nodes[n].y-ymin) / (ymax-ymin); 
    }

  elements = (struct element *)
             ALLOCMEM((nr - 1) * (nc - 1) * sizeof(struct element));

  trb = (struct el_tree *)
             ALLOCMEM((nr -1) * (nc - 1) * sizeof(struct el_tree));

  for(i = 0, n = 0; i < nr - 1; i++)
    for(j = 0; j < nc - 1; j++, n++)  {
      tr = &trb[n];
      el = tr -> entry = &elements[n];
      el -> node[0] = &nodes[nc*i + j];
      el -> node[1] = &nodes[nc*(i+1) + j];
      el -> node[2] = &nodes[nc*(i+1) + j + 1];
      el -> node[3] = &nodes[nc*i + j + 1];
      el -> d = 0; el -> z = 0;
      for(k = 0; k < 4; k++)  {
        el -> d += el -> node[k] -> d;
        el -> z += el -> node[k] -> z;
      }
      el -> d = (el -> d + 2) >> 2;
      tr -> left = tr -> right = NULL;
      if (element_el_tree == NULL)
        element_el_tree = tr;
      else
        C3D_Add_El_Tree(element_el_tree, tr);
    }

  GRA_GETMATRIX(gm);
  GRA_SETMATRIX(ident);
  GRA_WINDOW(0.0,4096.0,0.0,4096.0,-1.0,1.0);

  C3D_Show_El_Tree(element_el_tree);

  if (gra_state.pratio>0)
    GRA_PERSPECTIVE(gra_state.pratio);
  GRA_SETMATRIX(gm);
  GRA_FLUSH();

  FREEMEM(elements);
  FREEMEM(trb);
  FREEMEM(nodes);
}

void C3D_Levels(levels) int levels; 
/***********************************************************************
|  Set the number of colorlevels to be used in coloring
|  use 0 for singlecolor hidden surface plot
!  Number of levels must be between 0 - 13
^**********************************************************************/
{
  if (levels >= 0) 
    c3d_clevels = levels;
  else
    error("C3D_Levels: level parameter negative, not changed.\n");
}

void C3D_Persp(tf) int tf; 
/***********************************************************************
|  Use perspective or orthogonal transformation ?
!  Number of levels must be between 0 - 13
^**********************************************************************/
{
  c3d_perspective = tf;
}

void C3D_Add_El_Tree(head, add) struct el_tree *head, *add;
/***********************************************************************
|  Add a single element into the current leaf of the element tree
|  Only internal use
^**********************************************************************/
{
  if (add -> entry -> z > head -> entry -> z) {
    if (head -> left == NULL) {
      head -> left = add;
    } 
    else {
      C3D_Add_El_Tree(head -> left, add);
    }
  }
  else if (add -> entry -> z < head -> entry -> z) {
    if (head -> right == NULL) {
      head -> right = add;
    }
    else {
      C3D_Add_El_Tree(head -> right, add);
    }
  }
  else {
    add -> left = head -> left;
    head -> left = add;
  }
}   

void C3D_Show_El_Tree(head) struct el_tree *head;
/***********************************************************************
|  Display the contents of the current leaf of the element tree
|  Only internal use
^**********************************************************************/
{
  if (head == NULL) return;
  C3D_Show_El_Tree(head->left);
  C3D_Show_Elem(head->entry);
  C3D_Show_El_Tree(head->right);
}

void C3D_Free_El_Tree(head) struct el_tree *head;
/***********************************************************************
|  Free the memory allocated to the current leaf of the element tree
|  Only internal use
^**********************************************************************/
{
  if (head == NULL) return;
  C3D_Free_El_Tree(head->left);
  C3D_Free_El_Tree(head->right);
  FREEMEM(head);
}

int C3D_Convex_Test(x, y) int x[], y[];
/***********************************************************************
|  Test four-node ploygon 
|  Only internal use
^**********************************************************************/
{
   int a1, a2, at1, at2, amax, aind;

   amax = abs(y[0]*(x[2]-x[1])+y[1]*(x[0]-x[2])+y[2]*(x[1]-x[0]));
   aind = 3;

   at1 = abs(y[2]*(x[0]-x[3])+y[3]*(x[2]-x[0])+y[0]*(x[3]-x[2]));

   a1 = amax + at1;

   if (at1 > amax) { 
     amax = at1; aind = 1; 
   }

   at1 = abs(y[1]*(x[3]-x[2])+y[2]*(x[1]-x[3])+y[3]*(x[2]-x[1]));

   if (at1 > amax) { 
     amax = at1; aind = 0;
   }
 
   at2 = abs(y[3]*(x[1]-x[0])+y[0]*(x[3]-x[1])+y[1]*(x[0]-x[3]));

   if (at2 > amax) { 
     aind = 2;
   }

   a2 = at1 + at2;

   if (a1 == a2) return -1;

   return aind;
}

void C3D_Show_Elem(el) struct element *el;
/***********************************************************************
|  Display a single element by coloring and transformation
|  Only internal use
^**********************************************************************/
{
  int x[5], y[5], z[5], zz, xp, yp;
  int xi[5], yi[5], i;
  int d[5], zp, col[3];

  Point p[5];

  for(i = 0; i != 4; i++)
  {
    x[i] = el -> node[i] -> x;
    y[i] = el -> node[i] -> y;
    d[i] = el -> node[i] -> d;
  }

  if ((d[0]>>C3D_MASK) == (d[1]>>C3D_MASK) && 
      (d[0]>>C3D_MASK) == (d[2]>>C3D_MASK) && 
      (d[0]>>C3D_MASK) == (d[3]>>C3D_MASK))
  {
    C3D_SelCol(d[0]>>C3D_MASK);
    C3D_AreaFill(4, TRUE, x, y); 
    return;
  }

  switch(C3D_Convex_Test(x,y))
  {
    case 0:
      C3D_Show_Tri(x, y, d);
      xi[0]  = x[2]; xi[1]  = x[3]; xi[2]  = x[0];
      yi[0]  = y[2]; yi[1]  = y[3]; yi[2]  = y[0];
      col[0] = d[2]; col[1] = d[3]; col[2] = d[0];
      C3D_Show_Tri(xi, yi, col);
      break;

    case 1:
      C3D_Show_Tri(&x[1], &y[1], &d[1]);
      xi[0]  = x[0]; xi[1]  = x[1]; xi[2]  = x[3];
      yi[0]  = y[0]; yi[1]  = y[1]; yi[2]  = y[3];
      col[0] = d[0]; col[1] = d[1]; col[2] = d[3];
      C3D_Show_Tri(xi, yi, col);
      break;

    case 2:
      C3D_Show_Tri(x, y, d);
      xi[0]  = x[2]; xi[1]  = x[3]; xi[2]  = x[0];
      yi[0]  = y[2]; yi[1]  = y[3]; yi[2]  = y[0];
      col[0] = d[2]; col[1] = d[3]; col[2] = d[0];
      C3D_Show_Tri(xi, yi, col);
      break;

    case 3:
      C3D_Show_Tri(&x[1], &y[1], &d[1]);
      xi[0]  = x[0]; xi[1]  = x[1]; xi[2]  = x[3];
      yi[0]  = y[0]; yi[1]  = y[1]; yi[2]  = y[3];
      col[0] = d[0]; col[1] = d[1]; col[2] = d[3];
      C3D_Show_Tri(xi, yi, col);
      break;

    default:
      xp = 0; yp = 0;
      for(i = 0; i != 4; i++)
      {
        xp += x[i]; yp += y[i];
      }
      xp = (xp + 2) >> 2; yp = (yp + 2) >> 2;
      zp = el -> d;

      xi[0]  = x[0]; xi[1]  = x[1]; xi[2] = xp;
      yi[0]  = y[0]; yi[1]  = y[1]; yi[2] = yp;
      col[0] = d[0]; col[1] = d[1]; col[2] = zp;
      C3D_Show_Tri(xi, yi, col);

      xi[0]  = x[1]; xi[1]  = x[2]; 
      yi[0]  = y[1]; yi[1]  = y[2]; 
      col[0] = d[1]; col[1] = d[2];
      C3D_Show_Tri(xi, yi, col);
 
      xi[0]  = x[2]; xi[1]  = x[3];
      yi[0]  = y[2]; yi[1]  = y[3]; 
      col[0] = d[2]; col[1] = d[3];
      C3D_Show_Tri(xi, yi, col);
  
      xi[0]  = x[3]; xi[1]  = x[0];
      yi[0]  = y[3]; yi[1]  = y[0]; 
      col[0] = d[3]; col[1] = d[0];
      C3D_Show_Tri(xi, yi, col);
      break;
  }

  p[0].x = (int)(x[0]+0.5); p[0].y = (int)(y[0]+0.5); p[0].z = 0.0;
  p[1].x = (int)(x[1]+0.5); p[1].y = (int)(y[1]+0.5); p[1].z = 0.0;
  p[2].x = (int)(x[2]+0.5); p[2].y = (int)(y[2]+0.5); p[2].z = 0.0;
  p[3].x = (int)(x[3]+0.5); p[3].y = (int)(y[3]+0.5); p[3].z = 0.0;
  p[4].x = (int)(x[0]+0.5); p[4].y = (int)(y[0]+0.5); p[4].z = 0.0;
  GRA_COLOR(1);
  GRA_POLYLINE(5, p); 
}

void C3D_Show_Tri(x, y, d) int x[3], y[3], d[3];
/***********************************************************************
|  Only internal use
^**********************************************************************/
{
  int xx[128], yy[128], dd[128], px[7], py[7];
  int i, j, k, n = 0;

  if ((d[0] >> C3D_MASK) == (d[1] >> C3D_MASK) && 
      (d[0] >> C3D_MASK) == (d[2] >> C3D_MASK))
  {
    C3D_SelCol(d[0] >> C3D_MASK); 
    C3D_AreaFill(3, FALSE, x, y);
    return;
  }

  C3D_Pcalc(x[0],y[0],d[0],x[1],y[1],d[1],&i,&xx[n],&yy[n],&dd[n]); n += i;
  C3D_Pcalc(x[1],y[1],d[1],x[2],y[2],d[2],&i,&xx[n],&yy[n],&dd[n]); n += i;
  C3D_Pcalc(x[2],y[2],d[2],x[0],y[0],d[0],&i,&xx[n],&yy[n],&dd[n]); n += i;

  for(i = 0; i < 2; i++)
  {
    xx[n+i] = xx[i];
    yy[n+i] = yy[i]; 
    dd[n+i] = dd[i];
  }

  for(i = 0, k = 0; i < n - 2; k = 0, i++)
  {
    px[k] = xx[i];   py[k++] = yy[i];
    px[k] = xx[i+1]; py[k++] = yy[i+1];
 
    if (dd[i] == dd[i+1]) {
      i++; px[k] = xx[i+1]; py[k++] = yy[i+1];
    }

    for(j = n - 1; j > i; j--) {
      if (dd[i] == dd[j]) {
        if (dd[j-1] == dd[j]) {
          px[k] = xx[j-1];   py[k++] = yy[j-1];
        }
        px[k] = xx[j];   py[k++] = yy[j];
        px[k] = xx[j+1]; py[k++] = yy[j+1];
        if (dd[j] == dd[j+1]) {
          j++; px[k] = xx[j+1]; py[k++] = yy[j+1];
        }
        break;
      }
    }

    if (k > 2) {
      C3D_SelCol(dd[i]);
      C3D_AreaFill(k, FALSE, px, py);
    }  

  }  
}

void C3D_Pcalc(x0, y0, d0, x1, y1, d1, n, xx, yy, dd) 
/**/ int x0, y0, d0, x1, y1, d1, *n, xx[], yy[], dd[];
/***********************************************************************
|  Only internal use
^**********************************************************************/
{
  int x, y, deltax, deltay, deltad;
  int i, d, ds;

  *n = abs((d1 >> C3D_MASK) - (d0 >> C3D_MASK));

  xx[0] = x0; yy[0] = y0; dd[0] = d0 >> C3D_MASK;
  (*n)++;

  if (*n != 1) {

    deltad = (d1 >= d0 ? 1 : (-1));
    ds = d0 & ((1 << C3D_MASK)-1);

    if (d1 > d0) {
      ds = (1 << C3D_MASK) - ds;
    }

    d = abs(d1 - d0);

    if (x1 > x0) {
      deltax = ((x1 - x0) << (C3D_MASK << 1)) / d;
      deltax >>= C3D_MASK;
      x = x0 + ((ds * deltax + C3D_HALFMASK) >> C3D_MASK);
    }
    else {
      deltax = ((x0 - x1) << (C3D_MASK << 1)) /  d;
      deltax >>= C3D_MASK;
      x = x0 - ((ds * deltax + C3D_HALFMASK) >> C3D_MASK);
      deltax = -deltax;
    }

    if (y1 > y0) {
      deltay = ((y1 - y0) << (C3D_MASK << 1)) / d;
      deltay >>= C3D_MASK;
      y = y0 + ((ds * deltay + C3D_HALFMASK) >> C3D_MASK);
    }
    else {
      deltay = ((y0 - y1) << (C3D_MASK << 1)) / d;
      deltay >>= C3D_MASK;
      y = y0 - ((ds * deltay + C3D_HALFMASK) >> C3D_MASK);
      deltay = -deltay;
    }

    for(i = 1; i != *n; i++) {
      dd[i] = dd[i-1] + deltad;
      xx[i] = x; yy[i] = y;
      x = x + deltax; y = y + deltay;
    }
  }
}

void C3D_SelCol(col) int col;
/***********************************************************************
|  Only internal use
^**********************************************************************/
{
  GRA_COLOR(++col);
}

void C3D_AreaFill(n, border, x, y) int x[], y[], n, border;
/***********************************************************************
|  Only internal use
^**********************************************************************/
{
  int i, j;
  Point p[8];

  while(n > 0 && x[n-1] == x[0] && y[n-1] == y[0]) n--;

  for(i = 0; i < n; i++) {
    p[i].x = (int)(x[i]+0.5); p[i].y = (int)(y[i]+0.5); p[i].z = 0.0;
  }

  for(i = 0; i < n - 1; i++)
    if (p[i].x == p[i+1].x && p[i].y == p[i+1].y) {
      for(j = i + 1; j < n - 1; j++) {
        p[j].x = p[j+1].x; p[j].y = p[j+1].y;
      }
      n--;
    }

  if (n > 2) {
    GRA_AREAFILL(n, p);
    if (border)
    {
      p[n].x = p[0].x;
      p[n].y = p[0].y;
      p[n].z = 0;
      GRA_COLOR(1);
      GRA_POLYLINE(n+1, p);
    }
  }
}
