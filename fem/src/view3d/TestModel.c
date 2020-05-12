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

#include <ViewFactors.h>

/*******************************************************/
/*                    TEST SECTION                     */
/*******************************************************/

void elm_16node_quad_shape_functions(double a[16][16]);
void elm_4node_quad_shape_functions(double a[4][4]);

static double GeomNodes[MAX_GEOM_NODES][3],GeomNorms[MAX_GEOM_NODES][3];
static int    GeomElem[MAX_GEOM_ELEM][16],Type[MAX_GEOM_ELEM];

extern double ShapeFunctionMatrix[16][16],ShapeFunctionMatrix4[4][4],ShapeFunctionMatrix3[3][3];

void MakeTestModelCubic()
{
    double a,r,PI=2*acos(0.0),XMin,XMax,YMin,YMax,ZMin,ZMax;

    int i,j,k,n,CN=6*3,RN=4*3;

    elm_16node_quad_shape_functions( ShapeFunctionMatrix );
    elm_4node_quad_shape_functions(  ShapeFunctionMatrix4 );

    k = 0;
    for( i=0; i<=CN; i++ )
    for( j=0; j<=RN; j++,k++ )
    {
        a = 2.0*PI*i/CN;
        r = j/(double)RN;

        GeomNodes[k][0] = 0.25*cos(a);
        GeomNodes[k][1] = 0.25*sin(a);
        GeomNodes[k][2] = r;
    }

    for( i=0; i<k; i++ )
    {
        GeomNorms[i][0] = -GeomNodes[i][0];
        GeomNorms[i][1] = -GeomNodes[i][1];
        GeomNorms[i][2] =  0.0;

        r = sqrt( GeomNorms[i][0]*GeomNorms[i][0] + GeomNorms[i][1]*GeomNorms[i][1] );

        GeomNorms[i][0] /= r;
        GeomNorms[i][1] /= r;
        GeomNorms[i][2] /= r;
    }

    for( j=0; j<=CN; j++ )
    for( i=0; i<=RN; i++,k++ )
    {
        a = 2.0*PI*j/CN;
        r = 0.25*i/RN;

        GeomNodes[k][0] = r*cos(a);
        GeomNodes[k][1] = r*sin(a);
        GeomNodes[k][2] = 0.0;

        GeomNorms[k][0] = 0.0;
        GeomNorms[k][1] = 0.0;
        GeomNorms[k][2] = 1.0;
    }

    for( j=0; j<=CN; j++ )
    for( i=0; i<=RN; i++,k++ )
    {
        a = 2.0*PI*j/CN;
        r = 0.25*i/RN;

        GeomNodes[k][0] = r * cos(a);
        GeomNodes[k][1] = r * sin(a);
        GeomNodes[k][2] = 1.0;

        GeomNorms[k][0] =  0.0;
        GeomNorms[k][1] =  0.0;
        GeomNorms[k][2] = -1.0;
    }

    for( j=0; j<=CN; j++ )
    for( i=0; i<=RN; i++,k++ )
    {
        a = 2.0*PI*j/CN;
        r = 0.25*(i+10)/(RN+10);

        GeomNodes[k][0] = r * cos(a);
        GeomNodes[k][1] = r * sin(a);
        GeomNodes[k][2] = 1.0/2.0;

        GeomNorms[k][0] = 0.0;
        GeomNorms[k][1] = 0.0;
        GeomNorms[k][2] = 1.0;
    }

    for( j=0; j<=CN; j++ )
    for( i=0; i<=RN; i++,k++ )
    {
        a = 2.0*PI*j/CN;
        r = 0.25*(i+10)/(RN+10);

        GeomNodes[k][0] = r * cos(a);
        GeomNodes[k][1] = r * sin(a);
        GeomNodes[k][2] = 1.0/2.0 - 0.005;

        GeomNorms[k][0] =  0.0;
        GeomNorms[k][1] =  0.0;
        GeomNorms[k][2] = -1.0;
    }

/*
    for( i=0; i<4; i++,j+=3 )
    {
        GeomElem[i][0]  = j +  0;
        GeomElem[i][1]  = j +  3;
        GeomElem[i][2]  = j + 16;
        GeomElem[i][3]  = j + 13;
        GeomElem[i][4]  = j +  1;
        GeomElem[i][5]  = j +  2;
        GeomElem[i][6]  = j + 29;
        GeomElem[i][7]  = j + 42;
        GeomElem[i][8]  = j + 15;
        GeomElem[i][9]  = j + 14;
        GeomElem[i][10] = j + 39;
        GeomElem[i][11] = j + 26;
        GeomElem[i][12] = j + 27;
        GeomElem[i][13] = j + 28;
        GeomElem[i][14] = j + 41;
        GeomElem[i][15] = j + 40;
    }
*/

    n = 0;
{
    int d,s=0;

    for( d=0; d<5; d++, n+=(RN+1)*(CN+1) )
    for( j=0; j<CN/3; j++ )
    for( i=0; i<RN/3; i++,s++ )
    {
        GeomElem[s][0]  = n + 3*i + 3*(RN+1)*j + 0;
        GeomElem[s][1]  = n + 3*i + 3*(RN+1)*j + 3 ;
        GeomElem[s][2]  = n + 3*i + 3*(RN+1)*j + 3*(RN+1)+3;
        GeomElem[s][3]  = n + 3*i + 3*(RN+1)*j + 3*(RN+1); 
        GeomElem[s][4]  = n + 3*i + 3*(RN+1)*j + 1; 
        GeomElem[s][5]  = n + 3*i + 3*(RN+1)*j + 2; 
        GeomElem[s][6]  = n + 3*i + 3*(RN+1)*j + RN+1+3;
        GeomElem[s][7]  = n + 3*i + 3*(RN+1)*j + 2*(RN+1)+3; 
        GeomElem[s][8]  = n + 3*i + 3*(RN+1)*j + 3*(RN+1)+2; 
        GeomElem[s][9]  = n + 3*i + 3*(RN+1)*j + 3*(RN+1)+1; 
        GeomElem[s][10] = n + 3*i + 3*(RN+1)*j + 2*(RN+1); 
        GeomElem[s][11] = n + 3*i + 3*(RN+1)*j + (RN+1); 
        GeomElem[s][12] = n + 3*i + 3*(RN+1)*j + RN+1+1; 
        GeomElem[s][13] = n + 3*i + 3*(RN+1)*j + RN+1+2; 
        GeomElem[s][14] = n + 3*i + 3*(RN+1)*j + 2*(RN+1)+2; 
        GeomElem[s][15] = n + 3*i + 3*(RN+1)*j + 2*(RN+1)+1; 
    }

    NElements = s;
    fprintf( stderr, "ELEMENTS: %d\n", s );
}

fprintf( stderr, "%d\n", 1 );
   Elements = malloc( NElements*sizeof(Geometry_t) );

fprintf( stderr, "%d\n", 2 );
   for( i=0; i<NElements; i++ )
   {
        XMin = YMin = ZMin =  DBL_MAX;
        XMax = YMax = ZMax = -DBL_MAX;

        Elements[i].GeometryType = GEOMETRY_BICUBIC;
        Elements[i].BiCubic = (BiCubic_t *)calloc( sizeof(BiCubic_t),1 );

        for( j=0; j<16; j++ )
        {
            XMin = MIN( XMin,GeomNodes[GeomElem[i][j]][0] );
            YMin = MIN( YMin,GeomNodes[GeomElem[i][j]][1] );
            ZMin = MIN( ZMin,GeomNodes[GeomElem[i][j]][2] );

            XMax = MAX( XMax,GeomNodes[GeomElem[i][j]][0] );
            YMax = MAX( YMax,GeomNodes[GeomElem[i][j]][1] );
            ZMax = MAX( ZMax,GeomNodes[GeomElem[i][j]][2] );

            for( k=0; k< 16; k++ )
            for( n=0; n < 3; n++ )
            {
                Elements[i].BiCubic->PolyFactors[n][j]   += ShapeFunctionMatrix[k][j]*GeomNodes[GeomElem[i][k]][n];
                Elements[i].BiCubic->PolyFactors[n+3][j] += ShapeFunctionMatrix[k][j]*GeomNorms[GeomElem[i][k]][n];
            }
        }

        for( j=0; j<6; j++ )
        {
            BiCubicMonomialToBezier( Elements[i].BiCubic->PolyFactors[j],Elements[i].BiCubic->BezierFactors[j] );
        }

        Elements[i].BBox.XMin = XMin;
        Elements[i].BBox.XMax = XMax;

        Elements[i].BBox.YMin = YMin;
        Elements[i].BBox.YMax = YMax;

        Elements[i].BBox.ZMin = ZMin;
        Elements[i].BBox.ZMax = ZMax;
   }

fprintf( stderr, "%d\n", 3 );
#if 1
   fprintf( stderr, "SIZE: %ld\n", NElements*sizeof(Geometry_t) );
   Geometry = malloc( NElements*sizeof(Geometry_t) );
   memcpy( Geometry,Elements,NElements*sizeof(Geometry_t) );

fprintf( stderr, "%d\n", 4 );
   for( i=0; i<NElements; i++ )
   {
       Geometry[i].GeometryType = GEOMETRY_BICUBIC;

       Geometry[i].BBox.XMin = Elements[i].BBox.XMin;
       Geometry[i].BBox.XMax = Elements[i].BBox.XMax;
       Geometry[i].BBox.YMin = Elements[i].BBox.YMin;
       Geometry[i].BBox.YMax = Elements[i].BBox.YMax;
       Geometry[i].BBox.ZMin = Elements[i].BBox.ZMin;
       Geometry[i].BBox.ZMax = Elements[i].BBox.ZMax;

       Geometry[i].BiCubic = (BiCubic_t *)calloc( sizeof(BiCubic_t),1 );
       memcpy( Geometry[i].BiCubic,Elements[i].BiCubic,sizeof(BiCubic_t) );
   }
fprintf( stderr, "%d\n", 5 );
#endif

   NGeomElem = NElements;
   fprintf( stderr, "ELEMS: %d %d\n", NGeomElem, NElements );
}

void MakeTestModelLinear()
{
    double a,r,PI=2*acos(0.0),XMin,XMax,YMin,YMax,ZMin,ZMax,*nx,*ny,*nz;

    int i,j,k,n,NN,NE;

    FILE *fp = fopen( "qq.qq", "r" );
    char str[200];

    elm_16node_quad_shape_functions( ShapeFunctionMatrix );
    elm_4node_quad_shape_functions(  ShapeFunctionMatrix4 );

    ShapeFunctionMatrix3[0][0] =  1.0;
    ShapeFunctionMatrix3[0][1] = -1.0;
    ShapeFunctionMatrix3[0][2] = -1.0;

    ShapeFunctionMatrix3[1][0] =  0.0;
    ShapeFunctionMatrix3[1][1] =  1.0;
    ShapeFunctionMatrix3[1][2] =  0.0;

    ShapeFunctionMatrix3[2][0] =  0.0;
    ShapeFunctionMatrix3[2][1] =  0.0;
    ShapeFunctionMatrix3[2][2] =  1.0;

    fgets( str,100,fp ); 
    sscanf( str, "%d %d", &NN, &NE );

    for( j=0; j<NN; j++ )
    {
       fgets( str,100, fp );
       sscanf( str, "%lf %lf %lf", &GeomNodes[j][0],&GeomNodes[j][1],&GeomNodes[j][2] );
    }

    nx = (double *)malloc( NE*sizeof(double) );
    ny = (double *)malloc( NE*sizeof(double) );
    nz = (double *)malloc( NE*sizeof(double) );

    for( j=0; j<NE; j++ )
    {
       fgets( str,100, fp );
       sscanf( str, "%d %d %d %d %d %d %lf %lf %lf ", 
               &i, &Type[j],&GeomElem[j][0],&GeomElem[j][1],
          &GeomElem[j][2],&GeomElem[j][3], &nx[j],&ny[j],&nz[j] );
    }
 
    for( j=0; j<NE; j++ )
    {
       GeomNorms[j][0] = nx[j];
       GeomNorms[j][1] = ny[j];
       GeomNorms[j][2] = nz[j];
    }

    free( nx ); free( ny ); free( nz );

    Elements = (Geometry_t *)calloc( NE,sizeof(Geometry_t) );

    NElements = NE;
    fprintf( stderr, "ELEMENTS: %d\n", NElements );

    for( i=0; i<NElements; i++ )
    {
        if ( Type[i] == 404 )
        {
           Elements[i].GeometryType = GEOMETRY_BILINEAR;
           Elements[i].BiLinear = (BiLinear_t *)calloc( 1,sizeof(BiLinear_t) );

           for( j=0; j<4; j++ )
           {
              for( k=0; k<4; k++ )
              for( n=0; n<3; n++ )
              {
                  Elements[i].BiLinear->PolyFactors[n][j] +=
                      ShapeFunctionMatrix4[k][j]*GeomNodes[GeomElem[i][k]][n];

                  Elements[i].BiLinear->PolyFactors[n+3][j] +=
                      ShapeFunctionMatrix4[k][j]*GeomNorms[i][n];
              }
           }
        } else
        {
           Elements[i].GeometryType = GEOMETRY_TRIANGLE;
           Elements[i].Triangle = (Triangle_t *)calloc(sizeof(Triangle_t),1);

           for( j=0; j<3; j++ )
           {
              for( k=0; k<3; k++ )
              for( n=0; n<3; n++ )
              {
                  Elements[i].Triangle->PolyFactors[n][j]   +=
                      ShapeFunctionMatrix3[k][j]*GeomNodes[GeomElem[i][k]][n];

                  Elements[i].Triangle->PolyFactors[n+3][j] +=
                      ShapeFunctionMatrix3[k][j]*GeomNorms[i][n];
              }
           }
        }
    }


   NGeomElem = NElements;
   fprintf( stderr, "ELEMS: %d %d\n", NGeomElem, NElements );
}

/***************************************************************************************/

#if 0
void MakeTestModelTriangle()
{
    double a,r,PI=2*acos(0.0),XMin,XMax,YMin,YMax,ZMin,ZMax;

    int i,j,k,kk,n,NN,NE;

    FILE *fp = fopen( "qq.qq", "r" );
    char str[200];

    ShapeFunctionMatrix3[0][0] =  1.0;
    ShapeFunctionMatrix3[0][1] = -1.0;
    ShapeFunctionMatrix3[0][2] = -1.0;

    ShapeFunctionMatrix3[1][0] =  0.0;
    ShapeFunctionMatrix3[1][1] =  1.0;
    ShapeFunctionMatrix3[1][2] =  0.0;

    ShapeFunctionMatrix3[2][0] =  0.0;
    ShapeFunctionMatrix3[2][1] =  0.0;
    ShapeFunctionMatrix3[2][2] =  1.0;


    fgets( str,100,fp ); 
    sscanf( str, "%d %d", &NN, &NE );

    for( j=0; j<NN; j++ )
    {
       fgets( str,100, fp );
       sscanf( str, "%lf %lf %lf", &xx[j],&yy[j],&zz[j] );
    }
 
    for( j=0; j<NE; j++ )
    {
       fgets( str,100, fp );
       sscanf( str, "%*d %*d %d %d %d %d %lf %lf %lf ", &e[0][j],&e[1][j],&e[2][j],&e[3][j],&nx[j],&ny[j],&nz[j] );
    }
 
   
#if 1
    i = 0;
    for( j=0; j<NE; j++ )
    {
       for( k=0; k<4; k++,i++ )
       {
          GeomNodes[i][0] = xx[e[k][j]];
          GeomNodes[i][1] = yy[e[k][j]];
          GeomNodes[i][2] = zz[e[k][j]];

          GeomNorms[i][0] = nx[j];
          GeomNorms[i][1] = ny[j];
          GeomNorms[i][2] = nz[j];
       }
    }

    i = 0;
    for( j=0; j<NE; j++ )
      for( k=0; k<4; k++, i++ ) GeomElem[j][k] = i;
#else
    i = 0;
    for( j=0; j<NN; j++ )
    {
       GeomNodes[j][0] = xx[j];
       GeomNodes[j][1] = yy[j];
       GeomNodes[j][2] = zz[j];
    }

    i = 0;
    for( j=0; j<NE; j++ )
      for( k=0; k<4; k++, i++ ) GeomElem[j][k] = e[k][j];
#endif


    NElements = NE;
    fprintf( stderr, "ELEMENTS: %d\n", NElements );

   for( i=0; i<NElements; i++ )
   {
      Elements[2*i].GeometryType = GEOMETRY_TRIANGLE;
      Elements[2*i].Triangle = (Triangle_t *)calloc( sizeof(Triangle_t),1 );

      for( j=0; j<3; j++ )
      {
        for( k=0; k<3; k++ )
        for( n=0; n<3; n++ )
        {
          Elements[2*i].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*GeomNodes[GeomElem[i][k]][n];
#if 1
          Elements[2*i].Triangle->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*GeomNorms[GeomElem[i][k]][n];
#else
if ( n==0 ) {
                Elements[i].BiLinear->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*nx[i];
                Elements[i].BiLinear->PolyFactors[n+4][j] += ShapeFunctionMatrix3[k][j]*ny[i];
                Elements[i].BiLinear->PolyFactors[n+5][j] += ShapeFunctionMatrix3[k][j]*nz[i];
}
#endif
        }
      }

      Elements[2*i+1].GeometryType = GEOMETRY_TRIANGLE;
      Elements[2*i+1].Triangle = (Triangle_t *)calloc( sizeof(Triangle_t),1 );

      for( j=0; j<3; j++ )
      {
        for( k=0; k<3; k++ )
        for( n=0; n<3; n++ )
        {
          if ( k==0 )      kk=0;
          else if ( k==1 ) kk=2;
          else if ( k==2 ) kk=3;

          Elements[2*i+1].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*GeomNodes[GeomElem[i][kk]][n];
#if 1
          Elements[2*i+1].Triangle->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*GeomNorms[GeomElem[i][kk]][n];
#else
if ( n==0 ) {
                Elements[i].BiLinear->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*nx[i];
                Elements[i].BiLinear->PolyFactors[n+4][j] += ShapeFunctionMatrix3[k][j]*ny[i];
                Elements[i].BiLinear->PolyFactors[n+5][j] += ShapeFunctionMatrix3[k][j]*nz[i];
}
#endif
        }
      }  
   }

#if 0
   fprintf( stderr, "SIZE: %d\n", NElements*sizeof(Geometry_t) );
   memcpy( Geometry,Elements,NElements*sizeof(Geometry_t) );

   for( i=0; i<NElements; i++ )
   {
       Geometry[i].GeometryType = GEOMETRY_BICUBIC;

       Geometry[i].BBox.XMin = Elements[i].BBox.XMin;
       Geometry[i].BBox.XMax = Elements[i].BBox.XMax;
       Geometry[i].BBox.YMin = Elements[i].BBox.YMin;
       Geometry[i].BBox.YMax = Elements[i].BBox.YMax;
       Geometry[i].BBox.ZMin = Elements[i].BBox.ZMin;
       Geometry[i].BBox.ZMax = Elements[i].BBox.ZMax;

       Geometry[i].BiLinear = (BiLinear_t *)calloc( sizeof(BiLinear_t),1 );
       memcpy( Geometry[i].BiLinear,Elements[i].BiLinear,sizeof(BiLinear_t) );
   }
#endif

   NElements *= 2;
   NGeomElem = NElements;
   fprintf( stderr, "ELEMS: %d %d\n", NGeomElem, NElements );
}

/***************************************************************************************/
#endif
