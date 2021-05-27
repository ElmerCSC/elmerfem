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
 *  Element model type definitions etc.
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
 *                       Date: 20 Sep 1995
 *
 * Modification history:
 *
 * 28 Sep 1995, modified element_defs_t and element_type_t structures to hold
 *              list of element_types instead of an array
 * Juha R
 *
 *******************************************************************************/

#if !defined(ELEMENTS_H)

#define ELEMENTS_H

#define TRUE  1
#define FALSE 0

#ifdef MODULE_ELEMENTS
#define ELM_EXT 
#else
#define ELM_EXT extern
#endif

#define AEPS 1.0E-12

#define ELM_NULL_ELEMENT          -1

#define ELM_MAX_ELEMENT_TYPES    100
#define ELM_MAX_ELEMENT_CODE     999

#define ELM_MAX_ELEMENT_NODES 27  /* TODO:  FIX THIS WHEN YOU'VE GOT MORE... */

#define MAX_GROUP_IDS 8

typedef struct element_s
{
   struct element_type_s *ElementType;
   signed char DisplayFlag;
   int *Topology;
   signed char GroupIds[MAX_GROUP_IDS];
} element_t;

typedef struct element_model_s
{
   double *NodeArray;
   element_t *Elements;

   group_t *Groups;

   int NofNodes, NofElements, NofTimesteps;
} element_model_t;

ELM_EXT element_model_t ElementModel;
ELM_EXT element_t *Elements;

typedef struct element_type_s
{
    struct element_type_s *Next;

    char *ElementName;     /* One line description of the element */
    int   ElementCode;     /* Numeric code for the element        */

    double *NodeU;         /* node u coordinates */
    double *NodeV;         /* node v coordinates */
    double *NodeW;         /* node w coordinates */

    int NumberOfNodes;     /* number of nodes */

    /*
     * function to give value of a variable (f), given coordinates (u,v,w)
     */
    double (*FunctionValue)( double *f,double u,double v,double w );

    /*
     * function to give value of first partial derivate in (u) of a variable (f),
     * given coordinates (u,v,w)
     */
    double (*PartialU)( double *f,double u,double v,double w );

    /*
     * function to give value of first partial derivate in (v) of a variable (f),
     * given coordinates (u,v,w)
     */
    double (*PartialV)( double *f,double u,double v,double w );

    /*
     * function to give value of first partial derivate in (w) of a variable (f),
     * given coordinates (u,v,w)
     */
    double (*PartialW)( double *f,double u,double v,double w );

    /*
     * function to give value of second partial derivates of a variable (f),
     * given coordinates (u,v,w)
     */
    double (*SecondPartials)( double *f,double u,double v,double w,double *Values );

    /*
     * Trianglulate the element given node coordinates. Return value is 1 for
     * success, 0 for failure.
     */
    int (*Triangulate)( geometry_t *,element_t *,element_t * );

    /*
     * Check if a point is inside element boundaries, and return element coordinates
     * of the point if it is.
     */
    int (*PointInside)
         (
                       double *nodex, double *nodey, double *nodez,
                double x, double y, double z, double *u,double *v,double *w
         );

    /*
     * Isoline extraction for an element.
     */
    int (*IsoLine)
        (
           double K, double *F, double *C, double *nx, double *ny, double *nz, line_t *line
        );

    /*
     * Isosurface extraction for element.
     */
    int (*IsoSurface)
         (
               double K, double *F, double *C, double *nx, double *ny,double *nz,
                      double *nu,double *nv,double *nw,polygon_t *poly
         );

} element_type_t;

/*
 *  Element type definitions
 */
typedef struct element_defs_s
{
    element_type_t *ElementTypes;
    int NumberOfTypes;
} element_defs_t;

ELM_EXT element_defs_t ElementDefs;

#ifdef MODULE_ELEMENTS

int ElmBrickFace[6][9] = 
{
    { 0,1,2,3, 8, 9,10,11,20 },
    { 4,5,6,7,16,17,18,19,21 },
    { 0,1,5,4, 8,13,16,12,22 },
    { 3,2,6,7,10,14,18,15,24 },
    { 0,3,7,4,11,15,19,12,25 },
    { 1,2,6,5, 9,14,17,13,23 }
};

int ElmWedgeFace[5][8] = 
{
    { 0, 1, 4, 3,  6, 13,  9, 12 },
    { 0, 2, 5, 3,  8, 14, 11, 12 },
    { 1, 2, 5, 4,  7, 14, 10, 13 },
    { 0, 1, 2, 6,  7,  8,  0,  0 },
    { 3, 4, 5, 9, 10, 11,  0,  0 } 
};

int ElmTetraFace[4][7] = 
{
    { 0, 1, 2, 4, 5, 6, 10 },
    { 0, 1, 3, 4, 8, 7, 11 },
    { 1, 2, 3, 5, 9, 8, 12 },
    { 0, 2, 3, 6, 9, 7, 13 }
};

int ElmTetraFaceCubic[4][10] = 
{
    { 0, 1, 2, 4, 5, 6, 7, 8, 9,16 },
    { 0, 1, 3, 4, 5,11,14,13,10,17 },
    { 1, 2, 3, 6, 7,12,15,14,11,18 },
    { 0, 2, 3, 9, 8,12,15,13,10,19 }
};

#else

extern int ElmBrickFace[6][9];
extern int ElmWedgeFace[5][8];
extern int ElmTetraFace[4][7];
extern int ElmTetraFaceCubic[4][10];

#endif

void lu_mtrinv( double *, int );
int elm_initialize_element_types();
#endif
