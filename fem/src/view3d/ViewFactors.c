/*****************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
/******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 02 Jun 1997
 *
 *****************************************************************************/

#define MODULE_MAIN

/******************************************************************************

View factor computation.

Juha Ruokolainen/CSC - 24 Aug 1995

******************************************************************************/

#include <ViewFactors.h>
#include <Ipoints.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include "../../config.h"


/* (somewhat modified) copilot code.... */

#include <stdint.h>

typedef struct {
  uint64_t state;
  uint64_t inc;
} pcg32_rng_t;

static pcg32_rng_t **rbuf = NULL;
static int rbufn = 0;     /* allocated length of rbuf */
static int MPIRank = 0;   /* set by viewfactors3d before parallel region */

static inline uint32_t pcg32_random(pcg32_rng_t *rng)
{
    uint64_t old = rng->state;

    rng->state = old * 6364136223846793005ULL + (rng->inc | 1);
    uint32_t x = ((old >> 18u) ^ old) >> 27u;
    uint32_t r = old >> 59u;
    return (x >> r) | (x << ((-r) & 31));
}

/* Mix a work-item key into a well-separated RNG seed.  splitmix64 is the
 * usual companion seeder for pcg/xoshiro: consecutive keys give unrelated
 * streams, which is what we need since the keys here are consecutive
 * element indices. */
static inline uint64_t splitmix64( uint64_t x )
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/* Reseed the calling thread's stream from a work-item key.  Binding the
 * stream to the work item instead of to the thread is what makes the result
 * independent of the thread count and of the dynamic schedule: with a
 * per-thread stream, which rays a given element pair gets depends on which
 * thread happened to pick that pair up and how far that thread's stream had
 * already advanced. */
void vrand_seed( uint64_t key )
{
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    rbuf[tid]->state = splitmix64( key + 0x853c49e6748fea9bULL );
    rbuf[tid]->inc   = splitmix64( key + 0xda3e39cb94b95bdbULL );
}

/* Key for the element pair (a,b), symmetric so that the same physical pair
 * draws the same rays whichever of the two rows drives the integration. */
uint64_t vrand_pair_key( int a, int b, int n )
{
    int lo = (a<b) ? a : b, hi = (a<b) ? b : a;
    return (uint64_t)lo * (uint64_t)n + (uint64_t)hi;
}

inline double vrand()
{
#ifdef _OPENMP
    int tid  = omp_get_thread_num();
#else
    int tid  = 0;
#endif
    return pcg32_random(rbuf[tid]) / (double)UINT32_MAX;
}

void vrand_init()
{
#ifdef _OPENMP
   int tid = omp_get_thread_num(), tidn = omp_get_num_threads();
#else
   int tid = 0, tidn = 1;
#endif

/* The whole initialization is serialized: it runs once per thread per parallel
 * region, so the cost is irrelevant.  The first thread of a team to get here
 * sizes rbuf for the entire team, so the remaining threads never resize it --
 * no realloc can then race with a vrand() call from a thread that has already
 * finished its own init. */
#pragma omp critical
{
   int i;

   /* rbuf used to be sized once, from the team size of whichever parallel
    * region got here first, and never resized.  vrand_init is also called
    * from the serial radiator path (team size 1), so a later, wider region
    * would index past the end of the table.  Grow it instead. */
   if ( tidn > rbufn || tid >= rbufn ) {
     int newn = ( tidn > tid+1 ) ? tidn : tid+1;
     rbuf = realloc( rbuf, sizeof(pcg32_rng_t *) * newn );
     for( i=rbufn; i<newn; i++ ) rbuf[i] = NULL;
     rbufn = newn;
   }

   /* reuse the slot on a repeat call: reseeding is what matters, and
    * mallocing a fresh one each time just leaked the previous stream */
   if ( !rbuf[tid] ) rbuf[tid] = malloc(sizeof(pcg32_rng_t));

   /* Include MPI rank in seed so different ranks generate independent streams */
   rbuf[tid]->state = 0x853c49e6748fea9bULL + (uint64_t)MPIRank * 128 + tid;
   rbuf[tid]->inc   = 0xda3e39cb94b95bdbULL + (((uint64_t)MPIRank * 128 + tid) << 1);
}
}
/* end copilot code */

extern double ShapeFunctionMatrix2[2][2], ShapeFunctionMatrix3[3][3],ShapeFunctionMatrix4[4][4];


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
/*
 * iStart : global index of first source row assigned to this MPI rank (0 in serial)
 * nLocal : number of source rows assigned to this rank (= N in serial)
 * Factors: sized nLocal*N; row (i - iStart) stored at Factors[(i-iStart)*N + j]
 *
 * Cross-rank symmetry is not exploited: each rank computes its own rows fully
 * (all j != i) without writing to rows owned by other ranks.  Within-rank
 * symmetry (both i and j in [iStart, iStart+nLocal)) is preserved to halve
 * within-block computation.
 */
static void IntegrateFromGeometry(int NofRadiators, double *RadiatorCoords, int LineFlag,
                                   int N, double *Factors, int iStart, int nLocal)
{
    double T,s,F,Fmin=DBL_MAX,Fmax=-DBL_MAX,Favg=0.0,*RowSums,Fact,rx,ry,rz,nx,ny,nz,ct,realtime();
    int i,j,k,l,Imin,Imax,Ns;

    GeometryList_t *Link;


    ct = realtime();
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

      vrand_init();

      #pragma omp for private(i,j,k,l,Fact) schedule(dynamic,10)
      for( i=iStart; i<iStart+nLocal; i++ )
      {
         int li = i - iStart;   /* local (rank-relative) row index */

         /* Zero only entries not pre-filled by within-block symmetry.
          * Entries j in [iStart, i) are written when row j is processed
          * (upper-triangle pass, symmetric fill) — do not zero them here.
          *   j < iStart  : before our block, computed in lower-triangle loop
          *   j >= i      : standard upper triangle, not yet computed        */
         for( j=0;      j<iStart; j++ ) Factors[li*N+j] = 0.0;
         for( j=i;      j<N;      j++ ) Factors[li*N+j] = 0.0;
         if ( lel[i].Area<1.0e-10 ) continue;

         /* upper triangle: j > i — store own entry, also fill symmetric
          * entry when j is within this rank's local block               */
         for( j=i+1; j<N; j++ )
         {
            if ( lel[j].Area<1.0e-10 ) continue;

            FreeLinks( &lel[i] );
            FreeLinks( &lel[j] );

            lel[j].Flags |= GEOMETRY_FLAG_LEAF;
            lel[i].Flags |= GEOMETRY_FLAG_LEAF;

            vrand_seed( vrand_pair_key(i,j,N) );
            (*ViewFactorCompute[lel[i].GeometryType])( &lel[i],&lel[j],0,0 );
            Fact = ComputeViewFactorValue( &lel[i],0 );
            Factors[li*N+j] = Fact / lel[i].Area;

            /* fill F_ji only when row j is also owned by this rank */
            if ( j >= iStart && j < iStart+nLocal )
               Factors[(j-iStart)*N+i] = Fact / lel[j].Area;
         }

         /* lower triangle: j < i — only needed when j is NOT in our block
          * (those entries were already filled by symmetry above)         */
         for( j=0; j<i; j++ )
         {
            if ( j >= iStart && j < iStart+nLocal ) continue; /* done via symmetry */
            if ( lel[j].Area<1.0e-10 ) continue;

            FreeLinks( &lel[i] );
            FreeLinks( &lel[j] );

            lel[j].Flags |= GEOMETRY_FLAG_LEAF;
            lel[i].Flags |= GEOMETRY_FLAG_LEAF;

            vrand_seed( vrand_pair_key(i,j,N) );
            (*ViewFactorCompute[lel[i].GeometryType])( &lel[i],&lel[j],0,0 );
            Fact = ComputeViewFactorValue( &lel[i],0 );
            Factors[li*N+j] = Fact / lel[i].Area;
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

      vrand_init();

      for( i=0; i<NofRadiators; i++ )
      {
	 nx = ny = nz = 0;

	 rx = RadiatorCoords[i];
	 ry = RadiatorCoords[NofRadiators+i];
	 rz = RadiatorCoords[2*NofRadiators+i];

	 if ( LineFlag ) {
	   nx = RadiatorCoords[3*NofRadiators+i] - rx;
	   ny = RadiatorCoords[4*NofRadiators+i] - ry;
	   nz = RadiatorCoords[5*NofRadiators+i] - rz;
         }
         	 
         for( j=0; j<N; j++ ) Factors[i*N+j] = 0.0;
         for( j=0; j<N; j++ )
         { 
            if ( lel[j].Area<1.0e-10 ) continue;

            FreeLinks( &lel[j] );
            lel[j].Flags |= GEOMETRY_FLAG_LEAF;

            (*RadiatorFactorsCompute[lel[j].GeometryType])( &lel[j],LineFlag,rx,ry,rz,nx,ny,nz,0);

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
    if ( NofRadiators==0 ) Ns = nLocal;   /* only local rows */
    for(i=0; i<Ns; i++ )
    {
         s = 0.0;
         for( j=0; j<N; j++ ) s += Factors[i*N+j];   /* Factors row i is local row i */

         if ( s < Fmin )
         {
            Fmin = s;
            Imin = iStart + i + 1;   /* report global row number */
         }
         if ( s > Fmax )
         {
            Fmax = s;
            Imax = iStart + i + 1;
         }
         Favg += s;
    }
   fprintf( stdout, "rank %d: surfs: %d/%d, min(%d)=%-4.2f, max(%d)=%-4.2f, avg=%-4.2f, cput=%-4.2f\n",
                       MPIRank, nLocal, N, Imin,Fmin,Imax,Fmax,Favg/Ns, realtime()-ct );
}


/*
 * Tell the user whether the closed form inner integral actually got used.
 * A large miss count means patches are warped or carry non geometric normals,
 * and those pairs silently kept the old quadrature.
 */
void ReportClosedForm()
{
   long hits,miss,tot;

   if ( !ClosedFormInteg ) return;

   ContourCountSum( &hits,&miss );
   tot = hits + miss;
   if ( tot == 0 ) return;

   fprintf( stdout, "rank %d: closed form inner integral: %ld/%ld patch pairs "
                    "(%.1f%%), %ld fell back to quadrature\n",
                     MPIRank, hits, tot, 100.0*hits/tot, miss );
}


/*
 * How much visibility work would the clipping path actually have to do?
 * Bucket 0 is the prize: pairs with no candidate blocker at all, which need
 * no shadow handling of any kind and can take the closed form directly.
 */
void ReportShaftCull()
{
   static const char *Label[SC_NBUCKET] =
     { "none", "1-2", "3-4", "5-8", "9-16", ">16", "overflow",
       "MISSED BLOCKER (cull bug)", "no shaft (face away)" };
   long b[SC_NBUCKET], tot = 0;
   int i;

   if ( !ShaftStats ) return;

   ShaftCountSum( b );
   for( i=0; i<SC_NBUCKET; i++ ) if ( i != 7 ) tot += b[i];
   if ( tot == 0 ) return;

   fprintf( stdout, "rank %d: shaft cull over %ld resolved patch pairs:\n", MPIRank, tot );
   for( i=0; i<SC_NBUCKET; i++ )
      if ( b[i] )
         fprintf( stdout, "rank %d:   candidate blockers %-8s %10ld  (%5.1f%%)\n",
                    MPIRank, Label[i], b[i], 100.0*b[i]/tot );
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
}


void InitPolyFactors( int N, Geometry_t *Elements, int *Type, int *Topo, double *Coord, double *Normals, double *Data, int *Perm )
{
   int i,j,k,l,n;

   for( i=0; i<N; i++ )
   {
     switch( Type[i] ) {
     case 101:
	 j = Perm[i];
	 Elements[i].GeometryType = GEOMETRY_CIRCLE;
         Elements[i].Circle = (Circle_t *)calloc( sizeof(Circle_t),1 );
	 Elements[i].Circle->RMin = Data[8*j+0];
	 Elements[i].Circle->RMax = Data[8*j+1];
	 Elements[i].Circle->CenterPoint.x = Data[8*j+2];
	 Elements[i].Circle->CenterPoint.y = Data[8*j+3];
	 Elements[i].Circle->CenterPoint.z = Data[8*j+4];
	 GetMatrixToRotateVectorToZAxis(Data[8*i+5], Data[8*j+6], Data[8*j+7],
			 Elements[i].Circle->RotationMatrix, &Elements[i].Circle->IdentMatrix );
     break;
     case 202:
        Elements[i].GeometryType = GEOMETRY_LINE;
        Elements[i].Linear = (Linear_t *)calloc( sizeof(Linear_t),1 );

        for( j=0; j<2; j++ )
        {
           for( k=0; k<2; k++ )
           for( n=0; n<3; n++ )
           {
              l = 3*Topo[2*i+k]+n;
              Elements[i].Linear->PolyFactors[n][j]   += ShapeFunctionMatrix2[k][j]*Coord[l];
              if(Normals)Elements[i].Linear->PolyFactors[n+3][j] += ShapeFunctionMatrix2[k][j]*Normals[3*i+n];
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
              l = 3*Topo[4*i+k]+n;
              Elements[i].BiLinear->PolyFactors[n][j]   += ShapeFunctionMatrix4[k][j]*Coord[l];
              if(Normals)Elements[i].BiLinear->PolyFactors[n+3][j] += ShapeFunctionMatrix4[k][j]*Normals[3*i+n];
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
              l = 3*Topo[4*i+k]+n;
              Elements[i].Triangle->PolyFactors[n][j]   += ShapeFunctionMatrix3[k][j]*Coord[l];
              if(Normals)Elements[i].Triangle->PolyFactors[n+3][j] += ShapeFunctionMatrix3[k][j]*Normals[3*i+n];
           }
        }
     break;
     }
   }
}

void InitShapeFunctions()
{
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
}


void Combine2DRaytraceElements( int N, int *Topo, int *RT_N, int *RT_Topo, double *Coord )
{
   int maxind = 0, maxnodehits=0, tablesize;
   int *nodehits,*nodetable,i,j;
   int ind0,ind1,ind2;
   double x0, y0,  x1, y1, x2, y2, dx1, dx2,
     dy1, dy2, dp1, dp2, ds1,ds2, eps=1e-16;

   printf("Combining original boundary elements for shading\n");

   maxind=0; maxnodehits=0;
   for (i=0; i<2*N; i++) 
     if(maxind < Topo[i]) maxind = Topo[i];

   nodehits = (int*) malloc((maxind+1)*sizeof(int));
   for(i=0;i<=maxind;i++)
     nodehits[i] = 0;
   for(i=0; i<2*N; i++)  
     nodehits[Topo[i]]++;

   maxnodehits = 0;
   for(i=0; i<=maxind; i++) 
     if(nodehits[i] > maxnodehits) maxnodehits = nodehits[i];
    
   tablesize = (maxind+1)*maxnodehits;
   nodetable = (int*) malloc(tablesize*sizeof(int));
   for(i=0; i< tablesize; i++) 
     nodetable[i] = 0;

   for(i=0;i<=maxind;i++)
     nodehits[i] = 0;

   for(i=0; i<N; i++) {
     ind1 = Topo[2*i+1];
     ind2 = Topo[2*i+0];
     nodetable[maxnodehits*ind1 + nodehits[ind1]] = i;
     nodetable[maxnodehits*ind2 + nodehits[ind2]] = i;
     nodehits[ind1] += 1;
     nodehits[ind2] += 1;
   }
 
   for(i=0;i<2*N;i++)
     RT_Topo[i] = Topo[i];
 
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

     x0 = Coord[3*ind0];
     x1 = Coord[3*ind1];
     x2 = Coord[3*ind2];
      
     y0 = Coord[3*ind0+1];
     y1 = Coord[3*ind1+1];
     y2 = Coord[3*ind2+1];

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
   free((char*)(nodehits));

   j = 0;
   for (i=0; i<N; i++) {
     if(RT_Topo[2*i+1] || RT_Topo[2*i+0]) {
        RT_Topo[2*j+1] = RT_Topo[2*i+1];
        RT_Topo[2*j+0] = RT_Topo[2*i+0];
        j++;
     }
   }
   *RT_N = j;
   printf("The combined set includes %d line segments (vs. %d)\n",*RT_N,N);
}

void InitStuff( int N, int *Topo, int *Type, double *Coord, double *Normals, int RT_N0,
    int *RT_Topo0, double *RT_Data, int *RT_Perm, int *RT_Type, double *RT_Coord,
      double Feps, double Aeps, double Reps, int Nr, int NInteg2, int NInteg3, int NInteg4, int Combine,
        int ClosedForm, int ShaftStat, int Clip, int RayCull )
{
   int i,j,k,l,n,NOFRayElements;
   int RT_N=0, *RT_Topo=NULL;

   AreaEPS   = Aeps; 
   RayEPS    = Reps; 
   FactorEPS = Feps; 
   Nrays     = Nr;

   ClosedFormInteg = ClosedForm;
   ContourCountInit();

   ShaftStats = ShaftStat;
   ShaftCountInit();

   ClipShadows   = Clip;
   ShaftRayCull  = RayCull;      /* -1: decide below, once the shadow mesh
                                    size is known */

   InitShapeFunctions();

   Elements = (Geometry_t *)calloc(N,sizeof(Geometry_t));
   InitPolyFactors(N,Elements,Type,Topo,Coord,Normals, NULL, NULL);


   RT_N = RT_N0;
   RT_Topo = RT_Topo0;
   if ( RT_N == 0 && Combine )
   {
     if ( Type[0] == 202 ) {
       RT_Topo = (int *)malloc(2*N*sizeof(int));
       Combine2DRaytraceElements(N,Topo,&RT_N,RT_Topo,Coord);
     }
   }

   /*
    * check if different geometry elements given for shadowing ...
    */
   if ( RT_N > 0 ) {
     RTElements = (Geometry_t *)calloc( RT_N,sizeof(Geometry_t) );
     InitPolyFactors(RT_N,RTElements,RT_Type,RT_Topo,RT_Coord,NULL, RT_Data, RT_Perm);
     NOFRayElements = RT_N;
   } else {
     NOFRayElements = N;
     RTElements = Elements;
   }

   FillIPointArrays(NInteg2,NInteg3,NInteg4);
   InitGeometryTypes();
   InitVolumeBounds(2,NOFRayElements,RTElements);
   ShaftInitBoxes(NOFRayElements,RTElements);

   /*
    * Should the rays of a pair be culled to that pair's shaft candidates?
    *
    * The cull replaces Nrays tree traversals by one shaft build per pair, so
    * it pays when the traversals it saves cost more than the shaft does.  A
    * traversal is not O(1): the tree is built over a surface, so a ray meets
    * on the order of sqrt(NOFRayElements) cells on its way through, while
    * the shaft costs the same whatever the mesh.  Hence the product below
    * rather than a ray count alone -- a ray count alone gets it backwards on
    * a small shadow mesh, where the traversal is nearly free.
    *
    * Calibrated on the two ends of that range at one thread, timing the view
    * factor stage with the cull forced on and off:
    *
    *   shadow mesh   rays   Nrays*sqrt(N)   cull on / cull off
    *     6 (box_in_box)   8       20            0.62x   (a real loss)
    *     6               40       98            1.07x
    *   536 (radiation3d)  1       23            0.78x   (a real loss)
    *   536               4        93            1.16x
    *   536               8       185            1.98x
    *   536             100      2315            7.77x
    *
    * Both meshes break even near 70 despite being 89x apart in size, which
    * is what the sqrt scaling buys.  Near the threshold it is a wash either
    * way; well past it the win grows without bound, so a late switch costs
    * little and an early one costs 40%.
    */
   if ( ShaftRayCull < 0 )
      ShaftRayCull = ( Nrays*sqrt((double)NOFRayElements) >= 70.0 );
}

/* Fortran callable interface routines */

void radiatorfactors3d
  ( int *N,  int *Topo, int *Type, double *Coord, double *Normals, int *RT_N0, int *RT_Topo0, double *RT_Data,
      int *RT_Perm, int *RT_Type, double *RT_Coord, int *NofRadiators, double *RadiatorCoords, int *LineFlag,
        double *Factors, double *Feps, double *Aeps, double *Reps, int *Nr, int *NInteg2, int *NInteg3,int *NInteg4, int  *Combine,
          int *ClosedForm, int *ShaftStat, int *Clip, int *RayCull )
{
   /* Radiator factors: no MPI row decomposition (radiators typically few) */
   InitStuff( *N, Topo, Type, Coord, Normals, *RT_N0, RT_Topo0, RT_Data, RT_Perm, RT_Type, RT_Coord,
       *Feps, *Aeps, *Reps, *Nr, *NInteg2, *NInteg3, *NInteg4, *Combine, *ClosedForm, *ShaftStat, *Clip, *RayCull );

   *RayCull = ShaftRayCull;     /* tell the caller how an automatic (-1) went */

   IntegrateFromGeometry(*NofRadiators,RadiatorCoords,*LineFlag,*N,Factors,0,*NofRadiators);

   ReportClosedForm();
   ReportShaftCull();
}

/*
 * MPI-aware entry point.
 * iStart : 0-based global index of first source row on this rank
 * nLocal : number of source rows on this rank
 * mpiRank: MPI rank (for RNG seeding; 0 in serial)
 * Factors: caller allocates nLocal*N doubles
 */
void viewfactors3d
  ( int *N,  int *Topo, int *Type, double *Coord, double *Normals, int *RT_N0, int *RT_Topo0,
       double *RT_Data, int *RT_Perm, int *RT_Type, double *RT_Coord, double *Factors, double *Feps, double *Aeps,
          double *Reps, int *Nr, int *NInteg2,int *NInteg3, int *NInteg4, int *Combine,
          int *iStart, int *nLocal, int *mpiRank, int *ClosedForm, int *ShaftStat, int *Clip, int *RayCull )
{
   MPIRank = *mpiRank;

   InitStuff( *N, Topo, Type, Coord, Normals, *RT_N0, RT_Topo0, RT_Data, RT_Perm, RT_Type, RT_Coord,
            *Feps, *Aeps, *Reps, *Nr, *NInteg2, *NInteg3, *NInteg4, *Combine, *ClosedForm, *ShaftStat, *Clip, *RayCull );

   *RayCull = ShaftRayCull;     /* tell the caller how an automatic (-1) went */

   IntegrateFromGeometry(0,NULL,0,*N,Factors,*iStart,*nLocal);

   ReportClosedForm();
   ReportShaftCull();
}
