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

/******************************************************************************
 *
 *
 *
 ******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 02 Jun 1997
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "../../config.h"

#ifdef WIN32
double drand48();
#endif

#include <sys/types.h>

#ifdef MODULE_MAIN
#define EXT
#else
#define EXT extern
#endif

#ifndef MIN
#define MIN(x,y) ((x)>(y)?(y):(x))
#endif

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

#define ABS(s) ((s)>0?(s):(-(s)))

#define TRUE  1
#define FALSE 0

#define MAX_GEOM_ELEM  8000
#define MAX_GEOM_NODES 8000
#define MAX_BOXES       512

double Fvalue_4node( double *,double, double );
double ElementOfArea_4node( double *,double *,double *,double, double );

EXT double U_Integ1d[32],S_Integ1d[32],U_Integ[128],V_Integ[128],S_Integ[128];
EXT int N_Integ,N_Integ1d,N_Integ3;

EXT double ShapeFunctionMatrix[16][16],ShapeFunctionMatrix4[4][4],
           ShapeFunctionMatrix3[3][3], ShapeFunctionMatrix2[2][2];

EXT double U_Integ3[128],V_Integ3[128],S_Integ3[128];

EXT double XMin,XMax,YMin,YMax,ZMin,ZMax; 
EXT char str[512];

typedef double Matrix_t[3][3];

typedef struct
{
    double XMin,YMin,ZMin,XMax,YMax,ZMax;
} BBox_t;

typedef struct
{
    double x,y,z;
} Point_t;

typedef struct
{
    double Radius;
    Point_t CenterPoint;
} Sphere_t;

typedef struct
{
    /*
     * You can have a circle segment (area between RMin,RMax) 
     */
    double RMin,RMax;

    /*
     * define the rotation axis by direction vector and a point in the axis.
     * the point MUST be the centerpoint of the circle.
     */
    Point_t Axis;
    Point_t CenterPoint; 

    int IdentMatrix;
    Matrix_t RotationMatrix;
} Circle_t;

typedef struct
{
    double Radius,Length;

    /*
     * define the rotation axis by direction vector and a point in the axis.
     * the point MUST be the centerpoint of the cylinder.
     */
    Point_t Axis;
    Point_t CenterPoint;

    int IdentMatrix;
    Matrix_t RotationMatrix;
} Cylinder_t;

typedef struct
{
    /*
     * these define a quadratic equation r^2 = eqn[2]*z^2 + eqn[1]*z + eqn[0]
     */
    double RMin,RMax,ZMin,ZMax,REqn[3];

    /*
     * define the rotation axis by direction vector and a point in the axis
     */
    Point_t Axis;
    Point_t CenterPoint;

    int IdentMatrix;
    Matrix_t RotationMatrix;
} RotQuadric_t;

typedef struct
{
    int NumberOfSides;
    Point_t *Normals;
    Point_t *Vertices;
} Polygon_t;

typedef struct
{
    double  PolyFactors[6][16];
    double BezierFactors[6][16];
} BiCubic_t;

typedef struct
{
    double  PolyFactors[6][9];
    double BezierFactors[6][9];
} BiQuadratic_t;

typedef struct
{
    double  PolyFactors[6][4];
} BiLinear_t;

typedef struct
{
    double  PolyFactors[6][3];
} Triangle_t;

typedef struct
{
    double  PolyFactors[6][2];
} Linear_t;


#define MAX_GEOMETRY_TYPES   9

#define GEOMETRY_LINE        0
#define GEOMETRY_TRIANGLE    1
#define GEOMETRY_BILINEAR    2
#define GEOMETRY_BICUBIC     3
#define GEOMETRY_BIQUADRATIC 4
#define GEOMETRY_SPHERE      5
#define GEOMETRY_CYLINDER    6
#define GEOMETRY_CIRCLE      7
#define GEOMETRY_ROTQUADRIC  8

typedef struct GeometryList
{
    struct GeometryList *Next;

    double ViewFactor;
    struct Geometry *Entry;
} GeometryList_t;

#define GEOMETRY_FLAG_PLANE 1
#define GEOMETRY_FLAG_LEAF  2
#define GEOMETRY_FLAG_USED_INTEG 4

typedef struct Geometry
{
    int GeometryType;
    BBox_t BBox;

    union
    {
         Triangle_t *Triangle;
         BiQuadratic_t *BiQuadratic;
         BiLinear_t *BiLinear;
         Linear_t *Linear;
         BiCubic_t *BiCubic;
         Polygon_t *Polygon;
         Sphere_t *Sphere;
         Circle_t *Circle;
         Cylinder_t *Cylinder;
         RotQuadric_t *RotQuadric;
    } GeometryEntry;

    int Flags;
    int N;
    double Area,B,E,M;
	
    GeometryList_t *Link;
    struct Geometry *Left,*Right;
}  Geometry_t;

#define Triangle    GeometryEntry.Triangle
#define BiQuadratic GeometryEntry.BiQuadratic
#define BiLinear    GeometryEntry.BiLinear
#define Linear      GeometryEntry.Linear
#define BiCubic     GeometryEntry.BiCubic
#define Polygon     GeometryEntry.Polygon
#define Sphere      GeometryEntry.Sphere
#define Circle      GeometryEntry.Circle
#define Cylinder    GeometryEntry.Cylinder
#define RotQuadric  GeometryEntry.RotQuadric

typedef struct VolumeBounds
{
   struct VolumeBounds *Left,*Right;

   BBox_t BBox;
   int n;
   int *Elements;
} VolumeBounds_t;

EXT VolumeBounds_t *VolumeBounds;
EXT int NVolumeBounds,*VolumeTested;

EXT int NGeomElem,NElements;

EXT int (*RayHit[MAX_GEOMETRY_TYPES])(Geometry_t *,double,double,double,double,double,double,double);
EXT int RayHitGeometry(double,double,double,double,double,double );
void InitRayTracer(double);
void InitVolumeBounds( int,int,Geometry_t * );

void GetMatrixToRotateVectorToZAxis(double x,double y,double z,Matrix_t Matrix,int *Ident);
void RotateVector(double *x,double *y,double *z,Matrix_t Matrix);

EXT Geometry_t *Geometry,*Elements,*RTElements;

void BiCubicMonomialToBezier(double *MonomialFactors,double *BezierFactors);
void BiCubicBezierToMonomial(double *MonomialFactors,double *BezierFactors);

#define BiCubicValue(u,v,N) \
         (((((N[15]*u+N[14])*u+N[13])*u+N[12])*v+((N[11]*u+N[10])*u+N[9])*u+N[8])*v+ \
                ((N[7]*u+N[6])*u+N[5])*u+N[4])*v+((N[3]*u+N[2])*u+N[1])*u+N[0]

#define BiCubicPartialU(u,v,N) \
     (((((3*N[15]*v+3*N[11])*v+3*N[7])*v+3*N[3])*u+((2*N[14]*v+2*N[10])*v+2*N[6])*v+ \
                     2*N[2])*u+((N[13]*v+N[9])*v+N[5])*v+N[1])

#define BiCubicPartialV(u,v,N) \
    (((((3*N[15]*u+3*N[14])*u+3*N[13])*u+3*N[12])*v+((2*N[11]*u+2*N[10])*u+2*N[9])*u+ \
                    2*N[8])*v+((N[7]*u+N[6])*u+N[5])*u+N[4])

double BiCubicEofA(double U,double V,double *X,double *Y,double *Z);
double BiCubicArea(Geometry_t *);
double BiCubicLength(Geometry_t *,int);
double BiCubicIntegrateDiffToArea(Geometry_t *,double,double,double,double,double,double);
void BiCubicSubdivide(Geometry_t *,int,int);
int BiCubicIsAPlane(double *,double *,double *);
void BiCubicComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int Level,int );

void BiQuadraticMonomialToBezier(double *MonomialFactors,double *BezierFactors);
void BiQuadraticBezierToMonomial(double *MonomialFactors,double *BezierFactors);

#define BiQuadraticValue(u,v,N) (((N[8]*u+N[7])*u+N[6])*v+(N[5]*u+N[4])*u+N[3])*v+(N[2]*u+N[1])*u+N[0]
#define BiQuadraticPartialU(u,v,N) ((2*N[8]*u+N[7])*v+(2*N[5]*u+N[4]))*v+2*N[2]*u+N[1]
#define BiQuadraticPartialV(u,v,N) ((2*N[8]*u+2*N[7])*u+2*N[6])*v+(N[5]*u+N[4])*u+N[3]

double BiQuadraticEofA(double U,double V,double *X,double *Y,double *Z);
double BiQuadraticArea(Geometry_t *);
double BiQuadraticLength(Geometry_t *,int);
double BiQuadraticIntegrateDiffToArea(Geometry_t *,double,double,double,double,double,double);
void BiQuadraticSubdivide(Geometry_t *,int,int);
int BiQuadraticIsAPlane(double *,double *,double *);
void BiQuadraticComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int Level,int );

#define BiLinearValue(u,v,N) ((N[3]*u+N[2])*v+N[1]*u+N[0])
#define BiLinearPartialU(u,v,N) (N[3]*v+N[1])
#define BiLinearPartialV(u,v,N) (N[3]*u+N[2])

double BiLinearEofA(double U,double V,double *X,double *Y,double *Z);
double BiLinearArea(Geometry_t *);
double BiLinearLength(Geometry_t *,int);
double BiLinearIntegrateDiffToArea(Geometry_t *,double,double,double,double,double,double);
void BiLinearSubdivide(Geometry_t *,int,int);
void BiLinearComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int,int );

#define TriangleValue(u,v,N) (N[2]*v+N[1]*u+N[0])
#define TrianglePartialU(u,v,N) (N[1])
#define TrianglePartialV(u,v,N) (N[2])

double TriangleEofA(double U,double V,double *X,double *Y,double *Z);
double TriangleArea(Geometry_t *);
double TriangleLength(Geometry_t *,int);
double TriangleIntegrateDiffToArea(Geometry_t *,double,double,double,double,double,double);
void TriangleSubdivide(Geometry_t *,int,int);
void TriangleComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int,int );

#define LinearValue(u,N) (N[1]*u+N[0])
#define LinearPartialU(u,N) (N[1])

double LinearEofA(double,double, double *,double *, double *);
double LinearArea(Geometry_t *);
double LinearLength(Geometry_t *,int);
double LinearIntegrateDiffToArea(Geometry_t *,double,double,double,double,double,double);
void LinearSubdivide(Geometry_t *,int,int);
void LinearComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int,int );
void LinearComputeRadiatorFactors(Geometry_t *GA,double, double, double, int );
void BiLinearComputeRadiatorFactors(Geometry_t *GA,double, double, double, int );
void TriangleComputeRadiatorFactors(Geometry_t *GA,double, double, double, int );

void elm_4node_quad_shape_functions(double B[4][4]); 

static double FunctionValue( Geometry_t *Geom,double U,double V,int N )
{
	switch( Geom->GeometryType )
    {
  	   case GEOMETRY_LINE:
		   return LinearValue(U,Geom->Linear->PolyFactors[N]);

  	   case GEOMETRY_TRIANGLE:
		   return TriangleValue(U,V,Geom->Triangle->PolyFactors[N]);

  	   case GEOMETRY_BILINEAR:
		   return BiLinearValue(U,V,Geom->BiLinear->PolyFactors[N]);

  	   case GEOMETRY_BICUBIC:
		   return BiCubicValue(U,V,Geom->BiCubic->PolyFactors[N]); 

  	   case GEOMETRY_BIQUADRATIC:
		   return BiQuadraticValue(U,V,Geom->BiQuadratic->PolyFactors[N]);

   }

   return 0.0;
}

EXT double (*AreaCompute[MAX_GEOMETRY_TYPES])(Geometry_t *);
EXT void (*ViewFactorCompute[MAX_GEOMETRY_TYPES])(Geometry_t *,Geometry_t *,int, int);
EXT void (*RadiatorFactorsCompute[MAX_GEOMETRY_TYPES])(Geometry_t *,double,double,double,int);
EXT void (*Subdivide[MAX_GEOMETRY_TYPES])(Geometry_t *,int,int);
EXT double (*IntegrateDiffToArea[MAX_GEOMETRY_TYPES])(Geometry_t *,double,double,double,double,double,double);

void OutputGeometry( Geometry_t *,int );
void LinearSolveGaussSeidel( Geometry_t *,int,double *);

EXT double AreaEPS,FactorEPS,RayEPS;
EXT int hits, Nrays;

typedef struct CRSRows
{
   struct CRSRows *Head;
   struct CRSRows *Next;
   int col;
   double value;
} CRSRows_t;

typedef struct CRSMatrix {
   int N;
   CRSRows_t *Rows;
} CRSMatrix_t;

EXT CRSMatrix_t CRSFactors;
