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

#define MODULE_MAIN

/******************************************************************************

View factor computation.

Juha Ruokolainen/CSC - 24 Aug 1995

******************************************************************************/

#include <ViewFactors.h>
#include "../../config.h"

#if defined(WIN32) || defined(MINGW32) 
double drand48()
{
    return rand()/(1.0*(1<<15));
}
#endif

extern double ShapeFunctionMatrix3[3][3],ShapeFunctionMatrix4[4][4];


static int MaxLev;
/*******************************************************************************

Compute viewfactor from hierarchy

24 Aug 1995

*******************************************************************************/
static double ComputeViewFactorValue( Geometry_t *Geom,int Level )
{
     double S=0.0, Area=Geom->Area;

     GeometryList_t *Link;

     Link = Geom->Link;
     while( Link )
     {
        S += Area * Link->ViewFactor;
        Link = Link->Next;
     }

     if ( !(Geom->Flags & GEOMETRY_FLAG_LEAF) )
     {
        S += ComputeViewFactorValue( Geom->Left,Level+1 );
        S += ComputeViewFactorValue( Geom->Right,Level+1 );
     } else
       MaxLev = MAX(MaxLev,Level);

     return S;
}

static void PrintMesh( Geometry_t *Geom )
{
   static FILE *nfp=NULL,*efp=NULL;
   static int nodes;

   if ( !Geom->Left )
   {
     if ( !nfp )
     {
        nfp = (FILE *)fopen( "nodes", "w" );
        efp = (FILE *)fopen( "elems", "w" );
     }

     if ( Geom->GeometryType == GEOMETRY_TRIANGLE )
     {
         fprintf( nfp, "%g %g %g\n", 
            TriangleValue(0.0,0.0,Geom->Triangle->PolyFactors[0]),
            TriangleValue(0.0,0.0,Geom->Triangle->PolyFactors[1]),
            TriangleValue(0.0,0.0,Geom->Triangle->PolyFactors[2]) );

         fprintf( nfp, "%g %g %g\n", 
            TriangleValue(1.0,0.0,Geom->Triangle->PolyFactors[0]),
            TriangleValue(1.0,0.0,Geom->Triangle->PolyFactors[1]),
            TriangleValue(1.0,0.0,Geom->Triangle->PolyFactors[2]) );

         fprintf( nfp, "%g %g %g\n", 
            TriangleValue(0.0,1.0,Geom->Triangle->PolyFactors[0]),
            TriangleValue(0.0,1.0,Geom->Triangle->PolyFactors[1]),
            TriangleValue(0.0,1.0,Geom->Triangle->PolyFactors[2]) );

        fprintf( efp, "1 303 %d %d %d\n",nodes,nodes+1,nodes+2 );
        nodes += 3;
     } else {
         fprintf( nfp, "%g %g %g\n", 
            BiLinearValue(0.0,0.0,Geom->BiLinear->PolyFactors[0]),
            BiLinearValue(0.0,0.0,Geom->BiLinear->PolyFactors[1]),
            BiLinearValue(0.0,0.0,Geom->BiLinear->PolyFactors[2]) );

         fprintf( nfp, "%g %g %g\n", 
            BiLinearValue(1.0,0.0,Geom->BiLinear->PolyFactors[0]),
            BiLinearValue(1.0,0.0,Geom->BiLinear->PolyFactors[1]),
            BiLinearValue(1.0,0.0,Geom->BiLinear->PolyFactors[2]) );

         fprintf( nfp, "%g %g %g\n", 
            BiLinearValue(1.0,1.0,Geom->BiLinear->PolyFactors[0]),
            BiLinearValue(1.0,1.0,Geom->BiLinear->PolyFactors[1]),
            BiLinearValue(1.0,1.0,Geom->BiLinear->PolyFactors[2]) );

         fprintf( nfp, "%g %g %g\n", 
            BiLinearValue(0.0,1.0,Geom->BiLinear->PolyFactors[0]),
            BiLinearValue(0.0,1.0,Geom->BiLinear->PolyFactors[1]),
            BiLinearValue(0.0,1.0,Geom->BiLinear->PolyFactors[2]) );

        fprintf( efp, "1 404 %d %d %d %d\n",nodes,nodes+1,nodes+2,nodes+3 );
        nodes += 4;
     }
   } else {
      PrintMesh( Geom->Left );
      PrintMesh( Geom->Right );
   }
}

/*******************************************************************************
*******************************************************************************/
static void FreeLinks( Geometry_t *Geom )
{
     GeometryList_t *Link=Geom->Link, *Link1;

     while( Link )
     {
        Link1 = Link->Next;
        free(Link);
        Link = Link1;

     }
     Geom->Link = NULL;

     if ( Geom->Flags & GEOMETRY_FLAG_LEAF ) 
     {
        Geom->Flags &= ~GEOMETRY_FLAG_LEAF;
        return;
     }

     if ( Geom->Left )  FreeLinks( Geom->Left );
     if ( Geom->Right ) FreeLinks( Geom->Right );
}


/*******************************************************************************
*******************************************************************************/
static void FreeChilds( Geometry_t *Geom )
{
    if ( !Geom ) return;

     FreeChilds( Geom->Left );
     FreeChilds( Geom->Right );

    free( Geom );
}



/*******************************************************************************

Compute viewfactors for elements of the model

24 Aug 1995

*******************************************************************************/
static void IntegrateFromGeometry(int NofRadiators, double *RadiatorCoords, int N,double *Factors)
{
    double T,s,F,Fmin=DBL_MAX,Fmax=-DBL_MAX,Favg=0.0,*RowSums,Fact;
    int i,j,k,l,Imin,Imax,Ns;

    GeometryList_t *Link;

    for( i=0; i<N; i++ )
    {
        Elements[i].Area = (*AreaCompute[Elements[i].GeometryType])(&Elements[i]);
        Elements[i].Flags |= GEOMETRY_FLAG_LEAF;
    }

    MaxLev = 0;
    if (NofRadiators==0) {
#pragma omp parallel
{
      Geometry_t *lel;
      lel = (Geometry_t *)malloc(N*sizeof(Geometry_t));
      memcpy(lel,Elements,N*sizeof(Geometry_t));

      #pragma omp for private(i,j,k,l,Fact) schedule(dynamic,10)
      for( i=0; i<N; i++ )
      {
         for( j=i; j<N; j++ ) Factors[i*N+j] = 0.0;
         if ( lel[i].Area<1.0e-10 ) continue;

         for( j=i+1; j<N; j++ )
         { 
            if ( lel[j].Area<1.0e-10 ) continue;

            FreeLinks( &lel[i] );
            FreeLinks( &lel[j] );

            lel[j].Flags |= GEOMETRY_FLAG_LEAF;
            lel[i].Flags |= GEOMETRY_FLAG_LEAF;

            (*ViewFactorCompute[lel[i].GeometryType])( &lel[i],&lel[j],0,0 );
            Fact = ComputeViewFactorValue( &lel[i],0 );
            Factors[i*N+j] = Fact / lel[i].Area;
            Factors[j*N+i] = Fact / lel[j].Area;
         }

         FreeChilds( lel[i].Left );
         lel[i].Left = NULL;

         FreeChilds( lel[i].Right );
         lel[i].Right = NULL;
      }
      free(lel);
}
    } else {
      Geometry_t *lel;
      lel = (Geometry_t *)malloc(N*sizeof(Geometry_t));
      memcpy(lel,Elements,N*sizeof(Geometry_t));

      for( i=0; i<NofRadiators; i++ )
      {
         for( j=0; j<N; j++ ) Factors[i*N+j] = 0.0;
         for( j=0; j<N; j++ )
         { 
            if ( lel[j].Area<1.0e-10 ) continue;

            FreeLinks( &lel[j] );
            lel[j].Flags |= GEOMETRY_FLAG_LEAF;

            (*RadiatorFactorsCompute[lel[j].GeometryType])( &lel[j],
                RadiatorCoords[i], RadiatorCoords[NofRadiators+i], RadiatorCoords[2*NofRadiators+i], 0 );

            Fact = ComputeViewFactorValue( &lel[j],0 );
            Factors[i*N+j] = Fact; // lel[j].Area;

            FreeChilds( lel[j].Left );
            lel[j].Left = NULL;

            FreeChilds( lel[j].Right );
            lel[j].Right = NULL;
         }
      }

      free(lel);
    }

    k = 0;
    Ns = NofRadiators;
    if ( NofRadiators==0) Ns = N;
    for(i=0; i<Ns; i++ )
    {
         s = 0.0;
         for( j=0; j<N; j++ ) s += Factors[i*N+j];

         if ( s < Fmin )
         {
            Fmin = s;
            Imin = i+1;
         }
         if ( s > Fmax )
         {
            Fmax = s;
            Imax = i+1;
         }
         Favg += s;
    }
   fprintf( stdout, "surfs: %d, min(%d)=%-4.2f, max(%d)=%-4.2f, avg=%-4.2f\n", 
                       N,Imin,Fmin,Imax,Fmax,Favg/N );
}

void MakeRadiatorFactorMatrix(int NofRadiators, double *RadiatorCoords, int N,double *Factors,int NInteg,int NInteg3)
{
    double T[32],S[32];
    long int i,j,k,n;

    n = sqrt( NInteg ) + 0.5;

    switch( n )
    {
      case 1:
        T[0] = 0.0;
        S[0] = 2.0;

        N_Integ = 1;
      break;

      case 2:
        T[0] = -0.577350269189625;
        T[1] =  0.577350269189625;

        S[0] = 1.000000000000000;
        S[1] = 1.000000000000000;

        N_Integ = 2;
      break;

      case 4:
        T[0] = -0.861136311594052;
        T[1] = -0.339981043584856;
        T[2] =  0.339981043584856;
        T[3] =  0.861136311594052;

        S[0] =  0.347854845137454;
        S[1] =  0.652145154862546;
        S[2] =  0.652145154862546;
        S[3] =  0.347854845137454;

        N_Integ = 4;
      break;

      case 7:
        T[0] = -0.949107912342759;
        T[1] = -0.741531185599394;
        T[2] = -0.405845151377397;
        T[3] =  0.000000000000000;
        T[4] =  0.405845151377397;
        T[5] =  0.741531185599394;
        T[6] =  0.949107912342759;

        S[0] =  0.129484966168870;
        S[1] =  0.279705391489277;
        S[2] =  0.381830050505119;
        S[3] =  0.417959183673469;
        S[4] =  0.381830050505119;
        S[5] =  0.279705391489277;
        S[6] =  0.129484966168870;

        N_Integ = 7;
      break;
     }

     k = 0;
     for( i=0; i<N_Integ; i++ )
     {
        for( j=0; j<N_Integ; j++,k++ )
        {
           /* FOR -1,1 */
           U_Integ[k] = T[i];
           V_Integ[k] = T[j];
           S_Integ[k] = S[i]*S[j];

           /* FOR 0-1 */
           U_Integ[k] = 0.5*(T[i]+1.0);
           V_Integ[k] = 0.5*(T[j]+1.0);
           S_Integ[k] = S[i]*S[j]/4.0;
       }
       U_Integ1d[i] = 0.5*(T[i]+1.0);
       S_Integ1d[i] = 0.5*S[i];
    }
    N_Integ1d = N_Integ;
    N_Integ = k;


    switch( NInteg3 )
    {
      case 1:
         N_Integ3 = 1;

        U_Integ3[0] = 1.0/3.0;
        V_Integ3[0] = 1.0/3.0;
        S_Integ3[0] = 1.0/2.0;
      break;

      case 3:
        N_Integ3 = 3;

        U_Integ3[0] = 1.0/6.0;
        U_Integ3[1] = 2.0/3.0;
        U_Integ3[2] = 1.0/6.0;

        V_Integ3[0] = 1.0/6.0;
        V_Integ3[1] = 1.0/6.0;
        V_Integ3[2] = 2.0/3.0;

        S_Integ3[0] = 1.0/6.0;
        S_Integ3[1] = 1.0/6.0;
        S_Integ3[2] = 1.0/6.0;
      break;

      case 6:
        N_Integ3 = 6;

        U_Integ3[0] = 0.091576213509771;
        U_Integ3[1] = 0.816847572980459;
        U_Integ3[2] = 0.091576213509771;
        U_Integ3[3] = 0.445948490915965;
        U_Integ3[4] = 0.108103018168070;
        U_Integ3[5] = 0.445948490915965;

        V_Integ3[0] = 0.091576213509771;
        V_Integ3[1] = 0.091576213509771;
        V_Integ3[2] = 0.816847572980459;
        V_Integ3[3] = 0.445948490915965;
        V_Integ3[4] = 0.445948490915965;
        V_Integ3[5] = 0.108103018168070;

        S_Integ3[0] = 0.109951743655322 / 2;
        S_Integ3[1] = 0.109951743655322 / 2;
        S_Integ3[2] = 0.109951743655322 / 2;
        S_Integ3[3] = 0.223381589678011 / 2;
        S_Integ3[4] = 0.223381589678011 / 2;
        S_Integ3[5] = 0.223381589678011 / 2;
      break;
    }

    IntegrateFromGeometry( NofRadiators, RadiatorCoords, N,Factors );
}

void InitGeometryTypes()
{
    InitRayTracer( RayEPS );

    IntegrateDiffToArea[GEOMETRY_LINE]        = LinearIntegrateDiffToArea;
    IntegrateDiffToArea[GEOMETRY_TRIANGLE]    = TriangleIntegrateDiffToArea;
    IntegrateDiffToArea[GEOMETRY_BILINEAR]    = BiLinearIntegrateDiffToArea;
    IntegrateDiffToArea[GEOMETRY_BICUBIC]     = BiCubicIntegrateDiffToArea;
    IntegrateDiffToArea[GEOMETRY_BIQUADRATIC] = BiQuadraticIntegrateDiffToArea;

    Subdivide[GEOMETRY_LINE]                  = LinearSubdivide;
    Subdivide[GEOMETRY_TRIANGLE]              = TriangleSubdivide;
    Subdivide[GEOMETRY_BILINEAR]              = BiLinearSubdivide;
    Subdivide[GEOMETRY_BICUBIC]               = BiCubicSubdivide;
    Subdivide[GEOMETRY_BIQUADRATIC]           = BiQuadraticSubdivide;

    AreaCompute[GEOMETRY_LINE]                = LinearArea;
    AreaCompute[GEOMETRY_TRIANGLE]            = TriangleArea;
    AreaCompute[GEOMETRY_BILINEAR]            = BiLinearArea;
    AreaCompute[GEOMETRY_BICUBIC]             = BiCubicArea;
    AreaCompute[GEOMETRY_BIQUADRATIC]         = BiQuadraticArea;

    ViewFactorCompute[GEOMETRY_LINE]          = LinearComputeViewFactors;
    ViewFactorCompute[GEOMETRY_TRIANGLE]      = TriangleComputeViewFactors;
    ViewFactorCompute[GEOMETRY_BILINEAR]      = BiLinearComputeViewFactors;
    ViewFactorCompute[GEOMETRY_BICUBIC]       = BiCubicComputeViewFactors;
    ViewFactorCompute[GEOMETRY_BIQUADRATIC]   = BiQuadraticComputeViewFactors;

    RadiatorFactorsCompute[GEOMETRY_LINE]      = LinearComputeRadiatorFactors;
    RadiatorFactorsCompute[GEOMETRY_TRIANGLE]  = TriangleComputeRadiatorFactors;
    RadiatorFactorsCompute[GEOMETRY_BILINEAR]  = BiLinearComputeRadiatorFactors;
/*
    ViewFactorCompute[GEOMETRY_BICUBIC]       = BiCubicComputeViewFactors;
    ViewFactorCompute[GEOMETRY_BIQUADRATIC]   = BiQuadraticComputeViewFactors;
*/
}
void MakeViewFactorMatrix(int N,double *Factors,int NInteg,int NInteg3)
{
    double T[32],S[32];
    long int i,j,k,n;

    n = sqrt( NInteg ) + 0.5;

    switch( n )
    {
      case 1:
        T[0] = 0.0;
        S[0] = 2.0;

        N_Integ = 1;
      break;

      case 2:
        T[0] = -0.577350269189625;
        T[1] =  0.577350269189625;

        S[0] = 1.000000000000000;
        S[1] = 1.000000000000000;

        N_Integ = 2;
      break;

      case 4:
        T[0] = -0.861136311594052;
        T[1] = -0.339981043584856;
        T[2] =  0.339981043584856;
        T[3] =  0.861136311594052;

        S[0] =  0.347854845137454;
        S[1] =  0.652145154862546;
        S[2] =  0.652145154862546;
        S[3] =  0.347854845137454;

        N_Integ = 4;
      break;

      case 7:
        T[0] = -0.949107912342759;
        T[1] = -0.741531185599394;
        T[2] = -0.405845151377397;
        T[3] =  0.000000000000000;
        T[4] =  0.405845151377397;
        T[5] =  0.741531185599394;
        T[6] =  0.949107912342759;

        S[0] =  0.129484966168870;
        S[1] =  0.279705391489277;
        S[2] =  0.381830050505119;
        S[3] =  0.417959183673469;
        S[4] =  0.381830050505119;
        S[5] =  0.279705391489277;
        S[6] =  0.129484966168870;

        N_Integ = 7;
      break;
     }

     k = 0;
     for( i=0; i<N_Integ; i++ )
     {
        for( j=0; j<N_Integ; j++,k++ )
        {
           /* FOR -1,1 */
           U_Integ[k] = T[i];
           V_Integ[k] = T[j];
           S_Integ[k] = S[i]*S[j];

           /* FOR 0-1 */
           U_Integ[k] = 0.5*(T[i]+1.0);
           V_Integ[k] = 0.5*(T[j]+1.0);
           S_Integ[k] = S[i]*S[j]/4.0;
       }
       U_Integ1d[i] = 0.5*(T[i]+1.0);
       S_Integ1d[i] = 0.5*S[i];
    }
    N_Integ1d = N_Integ;
    N_Integ = k;


    switch( NInteg3 )
    {
      case 1:
         N_Integ3 = 1;

        U_Integ3[0] = 1.0/3.0;
        V_Integ3[0] = 1.0/3.0;
        S_Integ3[0] = 1.0/2.0;
      break;

      case 3:
        N_Integ3 = 3;

        U_Integ3[0] = 1.0/6.0;
        U_Integ3[1] = 2.0/3.0;
        U_Integ3[2] = 1.0/6.0;

        V_Integ3[0] = 1.0/6.0;
        V_Integ3[1] = 1.0/6.0;
        V_Integ3[2] = 2.0/3.0;

        S_Integ3[0] = 1.0/6.0;
        S_Integ3[1] = 1.0/6.0;
        S_Integ3[2] = 1.0/6.0;
      break;

      case 6:
        N_Integ3 = 6;

        U_Integ3[0] = 0.091576213509771;
        U_Integ3[1] = 0.816847572980459;
        U_Integ3[2] = 0.091576213509771;
        U_Integ3[3] = 0.445948490915965;
        U_Integ3[4] = 0.108103018168070;
        U_Integ3[5] = 0.445948490915965;

        V_Integ3[0] = 0.091576213509771;
        V_Integ3[1] = 0.091576213509771;
        V_Integ3[2] = 0.816847572980459;
        V_Integ3[3] = 0.445948490915965;
        V_Integ3[4] = 0.445948490915965;
        V_Integ3[5] = 0.108103018168070;

        S_Integ3[0] = 0.109951743655322 / 2;
        S_Integ3[1] = 0.109951743655322 / 2;
        S_Integ3[2] = 0.109951743655322 / 2;
        S_Integ3[3] = 0.223381589678011 / 2;
        S_Integ3[4] = 0.223381589678011 / 2;
        S_Integ3[5] = 0.223381589678011 / 2;
      break;
    }

    IntegrateFromGeometry( 0, NULL, N,Factors );
}


void radiatorfactors3d
  ( int *EL_N,  int *EL_Topo, int *EL_Type, double *EL_Coord, double *EL_Normals,
    int *RT_N0, int *RT_Topo0, int *RT_Type, double *RT_Coord, double *RT_Normals,
    int *NofRadiators, double *RadiatorCoords, double *Factors, double *Feps,
      double *Aeps, double *Reps, int *Nr, int *NInteg,int *NInteg3, int  *Combine )
{
   int i,j,k,l,n,NOFRayElements;
   int RT_N=0, *RT_Topo=NULL;

   AreaEPS   = *Aeps; 
   RayEPS    = *Reps; 
   FactorEPS = *Feps; 
   Nrays     = *Nr;

   ShapeFunctionMatrix2[0][0] =  1.0;
   ShapeFunctionMatrix2[0][1] = -1.0;

   ShapeFunctionMatrix2[1][0] =  0.0;
   ShapeFunctionMatrix2[1][1] =  1.0;

   ShapeFunctionMatrix3[0][0] =  1.0;
   ShapeFunctionMatrix3[0][1] = -1.0;
   ShapeFunctionMatrix3[0][2] = -1.0;

   ShapeFunctionMatrix3[1][0] =  0.0;
   ShapeFunctionMatrix3[1][1] =  1.0;
   ShapeFunctionMatrix3[1][2] =  0.0;

   ShapeFunctionMatrix3[2][0] =  0.0;
   ShapeFunctionMatrix3[2][1] =  0.0;
   ShapeFunctionMatrix3[2][2] =  1.0;

   elm_4node_quad_shape_functions(  ShapeFunctionMatrix4 );

   Elements = (Geometry_t *)calloc( *EL_N,sizeof(Geometry_t) );
   for( i=0; i<*EL_N; i++ )
   {
     switch( EL_Type[i] ) {
     case 202:
        Elements[i].GeometryType = GEOMETRY_LINE;
        Elements[i].Linear = (Linear_t *)calloc( sizeof(Linear_t),1 );

        for( j=0; j<2; j++ )
        {
           for( k=0; k<2; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*EL_Topo[2*i+k]+n;
              Elements[i].Linear->PolyFactors[n][j]   += ShapeFunctionMatrix2[k][j]*EL_Coord[l];
              Elements[i].Linear->PolyFactors[n+3][j] += ShapeFunctionMatrix2[k][j]*EL_Normals[3*i+n];
           }
        }
     break;
     case 404:
        Elements[i].GeometryType = GEOMETRY_BILINEAR;
        Elements[i].BiLinear = (BiLinear_t *)calloc( sizeof(BiLinear_t),1 );

        for( j=0; j<4; j++ )
        {
           for( k=0; k<4; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*EL_Topo[4*i+k]+n;
              Elements[i].BiLinear->PolyFactors[n][j]   += ShapeFunctionMatrix4[k][j]*EL_Coord[l];
              Elements[i].BiLinear->PolyFactors[n+3][j] += ShapeFunctionMatrix4[k][j]*EL_Normals[3*i+n];
           }
        }
     break;
     case 303:
        Elements[i].GeometryType = GEOMETRY_TRIANGLE;
        Elements[i].Triangle = (Triangle_t *)calloc( sizeof(Triangle_t),1 );

        for( j=0; j<3; j++ )
        {
           for( k=0; k<3; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*EL_Topo[3*i+k] + n;
              Elements[i].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*EL_Coord[l];
              Elements[i].Triangle->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*EL_Normals[3*i+n];
           }
        }
     break;
     }
   }

   RT_N = *RT_N0;
   RT_Topo = RT_Topo0;

   if ( RT_N == 0 && EL_Type[0] == 202 && *Combine )
   {
     int maxind = 0, maxnodehits=0, tablesize;
     int *nodehits,*nodetable,i,j;
     int ind0,ind1,ind2;
     double x0, y0,  x1, y1, x2, y2, dx1, dx2,
       dy1, dy2, dp1, dp2, ds1,ds2, eps=1e-16;

     printf("Combining original boundary elements for shading\n");

     for (i=0; i<2**EL_N; i++) 
       if(maxind < EL_Topo[i]) maxind = EL_Topo[i];

     nodehits = (int*) malloc((maxind+1)*sizeof(int));
     for(i=0;i<=maxind;i++)
       nodehits[i] = 0;
     for(i=0; i<2**EL_N; i++)  
       nodehits[EL_Topo[i]]++;

     maxnodehits = 0;
     for(i=0; i<=maxind; i++) 
       if(nodehits[i] > maxnodehits) maxnodehits = nodehits[i];
    
     tablesize = (maxind+1)*maxnodehits;
     nodetable = (int*) malloc(tablesize*sizeof(int));
     for(i=0; i< tablesize; i++) 
       nodetable[i] = 0;

     for(i=0;i<=maxind;i++)
       nodehits[i] = 0;

     for(i=0; i<*EL_N; i++) {
       ind1 = EL_Topo[2*i+1];
       ind2 = EL_Topo[2*i+0];
       nodetable[maxnodehits*ind1 + nodehits[ind1]] = i;
       nodetable[maxnodehits*ind2 + nodehits[ind2]] = i;
       nodehits[ind1] += 1;
       nodehits[ind2] += 1;
     }
 
     RT_Topo = (int*) malloc(2**EL_N*sizeof(int));
     for(i=0;i<2**EL_N;i++)
       RT_Topo[i] = EL_Topo[i];
 
     for (i=0; i<=maxind; i++) {
       int elem1,elem2;
       ind0 = i;
      
       if( nodehits[ind0] != 2) continue;

       elem1 = nodetable[maxnodehits*ind0+0];
       if( RT_Topo[2*elem1+1] == ind0 ) 
 	ind1 = RT_Topo[2*elem1];
       else 
 	ind1 = RT_Topo[2*elem1+1];

       elem2 = nodetable[maxnodehits*ind0+1];
       if( RT_Topo[2*elem2+1] == ind0 ) 
 	ind2 = RT_Topo[2*elem2];
       else 
 	ind2 = RT_Topo[2*elem2+1];

       x0 = EL_Coord[3*ind0];
       x1 = EL_Coord[3*ind1];
       x2 = EL_Coord[3*ind2];
      
       y0 = EL_Coord[3*ind0+1];
       y1 = EL_Coord[3*ind1+1];
       y2 = EL_Coord[3*ind2+1];

       dx1 = x1-x0;
       dx2 = x2-x0;
       dy1 = y1-y0;
       dy2 = y2-y0;
      
       dp1 = dx1*dx2+dy1*dy2;
       ds1 = sqrt(dx1*dx1+dy1*dy1);
       ds2 = sqrt(dx2*dx2+dy2*dy2);
      
       dp1 /= (ds1*ds2);

       /* Boundary elements must be aligned */
       if( dp1 > eps - 1. ) continue;

       /* Make the 1st element bigger  */
       if( RT_Topo[2*elem1] == ind0 ) 
 	 RT_Topo[2*elem1] = ind2;
       else 
	 RT_Topo[2*elem1+1] = ind2;
      
       /* Destroy the 2nd element */
       RT_Topo[2*elem2] = 0;
       RT_Topo[2*elem2+1] = 0;

       /* Update the node information */
       nodehits[ind0] = 0;
       if( nodetable[maxnodehits*ind2] == elem2) 
 	nodetable[maxnodehits*ind2] = elem1;
       else 
 	nodetable[maxnodehits*ind2+1] = elem1;
     }

     /* Free, not needed anymore */
     free((char*)(nodetable));

     j = 0;
     for (i=0; i<*EL_N; i++) {
       if(RT_Topo[2*i+1] || RT_Topo[2*i+0]) {
          RT_Topo[2*j+1] = RT_Topo[2*i+1];
          RT_Topo[2*j+0] = RT_Topo[2*i+0];
          j++;
       }
     }
     RT_N = j;
     printf("The combined set includes %d line segments (vs. %d)\n",RT_N,*EL_N);
   }

   /*
    * check if different geometry elements given for shadowing ...
    */
   if ( RT_N > 0 ) {
     RTElements = (Geometry_t *)calloc( RT_N,sizeof(Geometry_t) );
     for( i=0; i<RT_N; i++ )
     {
       switch( RT_Type[i] ) {
       case 202:
          RTElements[i].GeometryType = GEOMETRY_LINE;
          RTElements[i].Linear = (Linear_t *)calloc( sizeof(Linear_t),1 );

          for( j=0; j<2; j++ )
          {
             for( k=0; k<2; k++ )
             for( n=0; n<3; n++ )
             {
                l = 3*RT_Topo[2*i+k]+n;
                RTElements[i].Linear->PolyFactors[n][j]   += ShapeFunctionMatrix2[k][j]*RT_Coord[l];
              }
          }
       break;
       case 404:
          RTElements[i].GeometryType = GEOMETRY_BILINEAR;
          RTElements[i].BiLinear = (BiLinear_t *)calloc( sizeof(BiLinear_t),1 );

          for( j=0; j<4; j++ )
          {
             for( k=0; k<4; k++ )
             for( n=0; n<3; n++ )
             {
                l = 3*RT_Topo[4*i+k]+n;
                RTElements[i].BiLinear->PolyFactors[n][j]   += ShapeFunctionMatrix4[k][j]*RT_Coord[l];
             }
          }
       break;
       case 303:
          RTElements[i].GeometryType = GEOMETRY_TRIANGLE;
          RTElements[i].Triangle = (Triangle_t *)calloc( sizeof(Triangle_t),1 );

          for( j=0; j<3; j++ )
          {
             for( k=0; k<3; k++ )
             for( n=0; n<3; n++ )
             {
                l = 3*RT_Topo[3*i+k] + n;
                RTElements[i].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*RT_Coord[l];
             }
          }
       break;
       }
     }
     NOFRayElements = RT_N;
   } else {
     NOFRayElements = *EL_N;
     RTElements = Elements;
   }

   InitGeometryTypes();
   InitVolumeBounds( 2, NOFRayElements, RTElements );
   MakeRadiatorFactorMatrix( *NofRadiators, RadiatorCoords, *EL_N,Factors,*NInteg,*NInteg3 );


}

void viewfactors3d
  ( int *EL_N,  int *EL_Topo, int *EL_Type, double *EL_Coord, double *EL_Normals,
    int *RT_N0, int *RT_Topo0, int *RT_Type, double *RT_Coord, double *RT_Normals,
    double *Factors, double *Feps, double *Aeps, double *Reps, int *Nr,
    int *NInteg,int *NInteg3, int  *Combine )
{
   int i,j,k,l,n,NOFRayElements;
   int RT_N=0, *RT_Topo=NULL;

   AreaEPS   = *Aeps; 
   RayEPS    = *Reps; 
   FactorEPS = *Feps; 
   Nrays     = *Nr;

   ShapeFunctionMatrix2[0][0] =  1.0;
   ShapeFunctionMatrix2[0][1] = -1.0;

   ShapeFunctionMatrix2[1][0] =  0.0;
   ShapeFunctionMatrix2[1][1] =  1.0;

   ShapeFunctionMatrix3[0][0] =  1.0;
   ShapeFunctionMatrix3[0][1] = -1.0;
   ShapeFunctionMatrix3[0][2] = -1.0;

   ShapeFunctionMatrix3[1][0] =  0.0;
   ShapeFunctionMatrix3[1][1] =  1.0;
   ShapeFunctionMatrix3[1][2] =  0.0;

   ShapeFunctionMatrix3[2][0] =  0.0;
   ShapeFunctionMatrix3[2][1] =  0.0;
   ShapeFunctionMatrix3[2][2] =  1.0;

   elm_4node_quad_shape_functions(  ShapeFunctionMatrix4 );

   Elements = (Geometry_t *)calloc( *EL_N,sizeof(Geometry_t) );
   for( i=0; i<*EL_N; i++ )
   {
     switch( EL_Type[i] ) {
     case 202:
        Elements[i].GeometryType = GEOMETRY_LINE;
        Elements[i].Linear = (Linear_t *)calloc( sizeof(Linear_t),1 );

        for( j=0; j<2; j++ )
        {
           for( k=0; k<2; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*EL_Topo[2*i+k]+n;
              Elements[i].Linear->PolyFactors[n][j]   += ShapeFunctionMatrix2[k][j]*EL_Coord[l];
              Elements[i].Linear->PolyFactors[n+3][j] += ShapeFunctionMatrix2[k][j]*EL_Normals[3*i+n];
           }
        }
     break;
     case 404:
        Elements[i].GeometryType = GEOMETRY_BILINEAR;
        Elements[i].BiLinear = (BiLinear_t *)calloc( sizeof(BiLinear_t),1 );

        for( j=0; j<4; j++ )
        {
           for( k=0; k<4; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*EL_Topo[4*i+k]+n;
              Elements[i].BiLinear->PolyFactors[n][j]   += ShapeFunctionMatrix4[k][j]*EL_Coord[l];
              Elements[i].BiLinear->PolyFactors[n+3][j] += ShapeFunctionMatrix4[k][j]*EL_Normals[3*i+n];
           }
        }
     break;
     case 303:
        Elements[i].GeometryType = GEOMETRY_TRIANGLE;
        Elements[i].Triangle = (Triangle_t *)calloc( sizeof(Triangle_t),1 );

        for( j=0; j<3; j++ )
        {
           for( k=0; k<3; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*EL_Topo[3*i+k] + n;
              Elements[i].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*EL_Coord[l];
              Elements[i].Triangle->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*EL_Normals[3*i+n];
           }
        }
     break;
     }
   }

   RT_N = *RT_N0;
   RT_Topo = RT_Topo0;

   if ( RT_N == 0 && EL_Type[0] == 202 && *Combine )
   {
     int maxind = 0, maxnodehits=0, tablesize;
     int *nodehits,*nodetable,i,j;
     int ind0,ind1,ind2;
     double x0, y0,  x1, y1, x2, y2, dx1, dx2,
       dy1, dy2, dp1, dp2, ds1,ds2, eps=1e-16;

     printf("Combining original boundary elements for shading\n");

     for (i=0; i<2**EL_N; i++) 
       if(maxind < EL_Topo[i]) maxind = EL_Topo[i];

     nodehits = (int*) malloc((maxind+1)*sizeof(int));
     for(i=0;i<=maxind;i++)
       nodehits[i] = 0;
     for(i=0; i<2**EL_N; i++)  
       nodehits[EL_Topo[i]]++;

     maxnodehits = 0;
     for(i=0; i<=maxind; i++) 
       if(nodehits[i] > maxnodehits) maxnodehits = nodehits[i];
    
     tablesize = (maxind+1)*maxnodehits;
     nodetable = (int*) malloc(tablesize*sizeof(int));
     for(i=0; i< tablesize; i++) 
       nodetable[i] = 0;

     for(i=0;i<=maxind;i++)
       nodehits[i] = 0;

     for(i=0; i<*EL_N; i++) {
       ind1 = EL_Topo[2*i+1];
       ind2 = EL_Topo[2*i+0];
       nodetable[maxnodehits*ind1 + nodehits[ind1]] = i;
       nodetable[maxnodehits*ind2 + nodehits[ind2]] = i;
       nodehits[ind1] += 1;
       nodehits[ind2] += 1;
     }
 
     RT_Topo = (int*) malloc(2**EL_N*sizeof(int));
     for(i=0;i<2**EL_N;i++)
       RT_Topo[i] = EL_Topo[i];
 
     for (i=0; i<=maxind; i++) {
       int elem1,elem2;
       ind0 = i;
      
       if( nodehits[ind0] != 2) continue;

       elem1 = nodetable[maxnodehits*ind0+0];
       if( RT_Topo[2*elem1+1] == ind0 ) 
 	ind1 = RT_Topo[2*elem1];
       else 
 	ind1 = RT_Topo[2*elem1+1];

       elem2 = nodetable[maxnodehits*ind0+1];
       if( RT_Topo[2*elem2+1] == ind0 ) 
 	ind2 = RT_Topo[2*elem2];
       else 
 	ind2 = RT_Topo[2*elem2+1];

       x0 = EL_Coord[3*ind0];
       x1 = EL_Coord[3*ind1];
       x2 = EL_Coord[3*ind2];
      
       y0 = EL_Coord[3*ind0+1];
       y1 = EL_Coord[3*ind1+1];
       y2 = EL_Coord[3*ind2+1];

       dx1 = x1-x0;
       dx2 = x2-x0;
       dy1 = y1-y0;
       dy2 = y2-y0;
      
       dp1 = dx1*dx2+dy1*dy2;
       ds1 = sqrt(dx1*dx1+dy1*dy1);
       ds2 = sqrt(dx2*dx2+dy2*dy2);
      
       dp1 /= (ds1*ds2);

       /* Boundary elements must be aligned */
       if( dp1 > eps - 1. ) continue;

       /* Make the 1st element bigger  */
       if( RT_Topo[2*elem1] == ind0 ) 
 	 RT_Topo[2*elem1] = ind2;
       else 
	 RT_Topo[2*elem1+1] = ind2;
      
       /* Destroy the 2nd element */
       RT_Topo[2*elem2] = 0;
       RT_Topo[2*elem2+1] = 0;

       /* Update the node information */
       nodehits[ind0] = 0;
       if( nodetable[maxnodehits*ind2] == elem2) 
 	nodetable[maxnodehits*ind2] = elem1;
       else 
 	nodetable[maxnodehits*ind2+1] = elem1;
     }

     /* Free, not needed anymore */
     free((char*)(nodetable));

     j = 0;
     for (i=0; i<*EL_N; i++) {
       if(RT_Topo[2*i+1] || RT_Topo[2*i+0]) {
          RT_Topo[2*j+1] = RT_Topo[2*i+1];
          RT_Topo[2*j+0] = RT_Topo[2*i+0];
          j++;
       }
     }
     RT_N = j;
     printf("The combined set includes %d line segments (vs. %d)\n",RT_N,*EL_N);
   }

   /*
    * check if different geometry elements given for shadowing ...
    */
   if ( RT_N > 0 ) {
     RTElements = (Geometry_t *)calloc( RT_N,sizeof(Geometry_t) );
     for( i=0; i<RT_N; i++ )
     {
       switch( RT_Type[i] ) {
       case 202:
          RTElements[i].GeometryType = GEOMETRY_LINE;
          RTElements[i].Linear = (Linear_t *)calloc( sizeof(Linear_t),1 );

          for( j=0; j<2; j++ )
          {
             for( k=0; k<2; k++ )
             for( n=0; n<3; n++ )
             {
                l = 3*RT_Topo[2*i+k]+n;
                RTElements[i].Linear->PolyFactors[n][j]   += ShapeFunctionMatrix2[k][j]*RT_Coord[l];
              }
          }
       break;
       case 404:
          RTElements[i].GeometryType = GEOMETRY_BILINEAR;
          RTElements[i].BiLinear = (BiLinear_t *)calloc( sizeof(BiLinear_t),1 );

          for( j=0; j<4; j++ )
          {
             for( k=0; k<4; k++ )
             for( n=0; n<3; n++ )
             {
                l = 3*RT_Topo[4*i+k]+n;
                RTElements[i].BiLinear->PolyFactors[n][j]   += ShapeFunctionMatrix4[k][j]*RT_Coord[l];
             }
          }
       break;
       case 303:
          RTElements[i].GeometryType = GEOMETRY_TRIANGLE;
          RTElements[i].Triangle = (Triangle_t *)calloc( sizeof(Triangle_t),1 );

          for( j=0; j<3; j++ )
          {
             for( k=0; k<3; k++ )
             for( n=0; n<3; n++ )
             {
                l = 3*RT_Topo[3*i+k] + n;
                RTElements[i].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*RT_Coord[l];
             }
          }
       break;
       }
     }
     NOFRayElements = RT_N;
   } else {
     NOFRayElements = *EL_N;
     RTElements = Elements;
   }

   InitGeometryTypes();
   InitVolumeBounds( 2, NOFRayElements, RTElements );
   MakeViewFactorMatrix( *EL_N,Factors,*NInteg,*NInteg3 );
}

