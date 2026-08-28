/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library (in file ../LGPL-2.1); if not, write 
 *  to the Free Software Foundation, Inc., 51 Franklin Street, 
 *  Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/
/******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *****************************************************************************/

#include <ViewFactors.h>

/*******************************************************************************

Shaft culling: given a pair of patches, find the shadow mesh elements that
could possibly block the view between them.

Every ray from patch A to patch B stays inside the convex hull of A and B.
Three conservative tests are used, cheapest first:

  1. the combined bounding box of A and B,
  2. the supporting planes of the convex hull of A and B.  Each candidate is
     spanned by three of the corners -- the two patch planes themselves, and
     the lateral faces given by an edge of one patch and a vertex of the
     other -- and is kept only when every vertex of both patches lies on one
     side of it.  That one sided test is what makes it a hull face.

Note that a patch plane is *not* automatically a face of the hull: when the
patches straddle each other's planes, neither is supporting, and using one
anyway would cull real blockers.  The same one sided test decides all of
them, so this case needs no special handling -- it just yields fewer planes.

Both tests are exact half spaces of the hull, so what survives is the true
shaft up to the looseness of the box in test 1.

Surviving elements then face one more, exact, test of their own: a surface
element can only block a ray from A to B if A and B lie on opposite sides of
its plane.  In a closed enclosure that is what removes the side walls, which
a convex shaft across a cavity otherwise always contains.  A survivor of all
of it is a candidate, never a certainty.

The result is the candidate list the clipping path will project and subtract.
On its own it already separates the pairs that need no visibility work at all
(empty list -> the view is unobstructed, use the closed form directly) from
the ones that do.

Elements coplanar with either patch are dropped: they lie in the boundary of
a half space and can only graze, never block.  That is what removes a patch's
own shadow mesh element, and its flat neighbours, from its own candidate list.

Juha Ruokolainen/CSC

*******************************************************************************/

#define SC_EPS 1.0e-12

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Bucket counters, padded per thread like the contour ones so that counting
 * cannot put a shared cache line in the element pair loop.
 * Bucket 0 counts pairs with no candidate blocker at all.
 */
void ShaftCountInit()
{
    int n;
#ifdef _OPENMP
    n = omp_get_max_threads();
#else
    n = 1;
#endif
    if ( n < 1 ) n = 1;
    ShaftNThreads = n;

    if ( ShaftCount ) free( ShaftCount );
    ShaftCount = (long *)calloc( n*SC_NBUCKET*CF_STRIDE, sizeof(long) );
}

void ShaftCountAdd( int n )
{
    int tid = 0, b;

    if ( !ShaftCount ) return;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if ( tid < 0 || tid >= ShaftNThreads ) return;

    if      ( n == -3) b = 8;          /* no shaft: patches face away        */
    else if ( n == -2) b = 7;          /* cull said clear, a ray was blocked */
    else if ( n <  0 ) b = 6;          /* overflowed the candidate list */
    else if ( n == 0 ) b = 0;
    else if ( n <= 2 ) b = 1;
    else if ( n <= 4 ) b = 2;
    else if ( n <= 8 ) b = 3;
    else if ( n <= 16) b = 4;
    else               b = 5;

    ShaftCount[ (tid*SC_NBUCKET + b)*CF_STRIDE ]++;
}

void ShaftCountSum( long *Buckets )
{
    int i,b;

    for( b=0; b<SC_NBUCKET; b++ ) Buckets[b] = 0;
    if ( !ShaftCount ) return;

    for( i=0; i<ShaftNThreads; i++ )
      for( b=0; b<SC_NBUCKET; b++ )
        Buckets[b] += ShaftCount[ (i*SC_NBUCKET + b)*CF_STRIDE ];
}

/*
 * Per element bounding boxes of the shadow mesh.  Built once, read only
 * afterwards, so every thread shares one copy.
 */
void ShaftInitBoxes( int N, Geometry_t *El )
{
    double U[] = { 0.0,1.0,0.0,1.0 }, V[] = { 0.0,0.0,1.0,1.0 };
    double x,y,z,R;
    int i,j,nc;

    if ( RTElementBBox ) free( RTElementBBox );
    RTElementBBox = (BBox_t *)malloc( N*sizeof(BBox_t) );
    RTElementNof  = N;

    for( i=0; i<N; i++ )
    {
        BBox_t *B = &RTElementBBox[i];

        B->XMin = B->YMin = B->ZMin =  1.0e20;
        B->XMax = B->YMax = B->ZMax = -1.0e20;

        switch( El[i].GeometryType )
        {
        case GEOMETRY_CIRCLE:
            R = El[i].Circle->RMax;
            B->XMin = El[i].Circle->CenterPoint.x - R;
            B->XMax = El[i].Circle->CenterPoint.x + R;
            B->YMin = El[i].Circle->CenterPoint.y - R;
            B->YMax = El[i].Circle->CenterPoint.y + R;
            B->ZMin = B->ZMax = El[i].Circle->CenterPoint.z;
            continue;
        case GEOMETRY_LINE:     nc = 2; break;
        case GEOMETRY_TRIANGLE: nc = 3; break;
        default:                nc = 4; break;
        }

        for( j=0; j<nc; j++ )
        {
            x = FunctionValue( &El[i],U[j],V[j],0 );
            y = FunctionValue( &El[i],U[j],V[j],1 );
            z = FunctionValue( &El[i],U[j],V[j],2 );

            B->XMin = MIN(x,B->XMin); B->XMax = MAX(x,B->XMax);
            B->YMin = MIN(y,B->YMin); B->YMax = MAX(y,B->YMax);
            B->ZMin = MIN(z,B->ZMin); B->ZMax = MAX(z,B->ZMax);
        }
    }
}

/* TRUE when the box lies wholly in the negative half space { N.x < D }. */
static int BoxBehind( BBox_t *B, double *N, double D )
{
    double x = (N[0] > 0.0) ? B->XMax : B->XMin;
    double y = (N[1] > 0.0) ? B->YMax : B->YMin;
    double z = (N[2] > 0.0) ? B->ZMax : B->ZMin;

    return ( N[0]*x + N[1]*y + N[2]*z < D );
}

static int BoxDisjoint( BBox_t *a, BBox_t *b )
{
    return a->XMax < b->XMin || a->XMin > b->XMax ||
           a->YMax < b->YMin || a->YMin > b->YMax ||
           a->ZMax < b->ZMin || a->ZMin > b->ZMax;
}

/* TRUE when the box passes every shaft test, i.e. it might contain a blocker. */
static int BoxInShaft( Shaft_t *S, BBox_t *B )
{
    int i;

    if ( BoxDisjoint(B,&S->BBox) )   return FALSE;

    for( i=0; i<S->NP; i++ )
       if ( BoxBehind(B,S->N[i],S->D[i]) ) return FALSE;

    return TRUE;
}

/*
 * Try the plane spanned by the edge (P,Q) and the point R as a face of the
 * convex hull of the two patches: keep it when every vertex of both patches
 * is on one side of it.
 */
static void TryHullPlane( Shaft_t *S, ContourTarget_t *A, ContourTarget_t *B,
                            double *P, double *Q, double *R )
{
    double n[3],d,u[3],v[3],L,t;
    int i,k,pos = 0,neg = 0;

    if ( S->NP >= SH_MAXPLANE ) return;

    for( k=0; k<3; k++ ) { u[k] = Q[k]-P[k]; v[k] = R[k]-P[k]; }

    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];

    L = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
    if ( L < SC_EPS ) return;                    /* degenerate triple */
    for( k=0; k<3; k++ ) n[k] /= L;

    d = n[0]*P[0] + n[1]*P[1] + n[2]*P[2];

    for( i=0; i<A->NV+B->NV; i++ )
    {
        double *x = (i < A->NV) ? A->V[i] : B->V[i-A->NV];

        t = n[0]*x[0] + n[1]*x[1] + n[2]*x[2] - d;
        if ( t >  S->Tol ) pos++;
        if ( t < -S->Tol ) neg++;
    }

    if ( pos && neg ) return;                    /* not a hull face */
    if ( !pos && !neg ) return;                  /* everything coplanar */

    if ( neg ) { for( k=0; k<3; k++ ) n[k] = -n[k]; d = -d; }

    for( k=0; k<3; k++ ) S->N[S->NP][k] = n[k];
    S->D[S->NP] = d - S->Tol;
    S->NP++;
}

/*
 * Build the shaft of the pair (A,B).  Returns FALSE when the two patches do
 * not see each other at all, in which case there is nothing to compute.
 */
int ShaftInit( Shaft_t *S, ContourTarget_t *A, ContourTarget_t *B )
{
    double CA[3],CB[3],tol;
    int i,k;

    /* combined bounding box and the two centroids */
    S->BBox.XMin = S->BBox.YMin = S->BBox.ZMin =  1.0e20;
    S->BBox.XMax = S->BBox.YMax = S->BBox.ZMax = -1.0e20;

    for( k=0; k<3; k++ ) CA[k] = CB[k] = 0.0;

    for( i=0; i<A->NV; i++ )
    {
        S->BBox.XMin = MIN(A->V[i][0],S->BBox.XMin);
        S->BBox.XMax = MAX(A->V[i][0],S->BBox.XMax);
        S->BBox.YMin = MIN(A->V[i][1],S->BBox.YMin);
        S->BBox.YMax = MAX(A->V[i][1],S->BBox.YMax);
        S->BBox.ZMin = MIN(A->V[i][2],S->BBox.ZMin);
        S->BBox.ZMax = MAX(A->V[i][2],S->BBox.ZMax);
        for( k=0; k<3; k++ ) CA[k] += A->V[i][k]/A->NV;
    }
    for( i=0; i<B->NV; i++ )
    {
        S->BBox.XMin = MIN(B->V[i][0],S->BBox.XMin);
        S->BBox.XMax = MAX(B->V[i][0],S->BBox.XMax);
        S->BBox.YMin = MIN(B->V[i][1],S->BBox.YMin);
        S->BBox.YMax = MAX(B->V[i][1],S->BBox.YMax);
        S->BBox.ZMin = MIN(B->V[i][2],S->BBox.ZMin);
        S->BBox.ZMax = MAX(B->V[i][2],S->BBox.ZMax);
        for( k=0; k<3; k++ ) CB[k] += B->V[i][k]/B->NV;
    }

    /* a blocker has to be strictly inside, so back the planes off by a
       tolerance scaled to the shaft itself */
    tol = MAX( S->BBox.XMax-S->BBox.XMin,
          MAX( S->BBox.YMax-S->BBox.YMin, S->BBox.ZMax-S->BBox.ZMin ) );
    S->Tol = 1.0e-8*tol + SC_EPS;

    S->A = *A;
    S->B = *B;
    S->NP = 0;

    /*
     * Every supporting plane of the hull, by one uniform test: the two patch
     * planes, then every edge of one patch against every vertex of the other.
     */
    TryHullPlane( S,A,B, A->V[0],A->V[1],A->V[2] );
    TryHullPlane( S,A,B, B->V[0],B->V[1],B->V[2] );

    for( i=0; i<A->NV; i++ )
      for( k=0; k<B->NV; k++ )
        TryHullPlane( S,A,B, A->V[i],A->V[(i+1)%A->NV],B->V[k] );

    for( i=0; i<B->NV; i++ )
      for( k=0; k<A->NV; k++ )
        TryHullPlane( S,A,B, B->V[i],B->V[(i+1)%B->NV],A->V[k] );

    return TRUE;
}

/*
 * Exact rejection of a single shadow element.
 *
 * A segment from A to B crosses the element's plane only if A and B are on
 * opposite sides of it.  If every vertex of both patches is on one side, no
 * ray between them can meet the element and it cannot block, whatever its
 * extent.  A patch's own element, and everything coplanar with either patch,
 * falls out of this as the degenerate case where one side is the plane
 * itself: grazing does not block.
 */
static int CannotBlock( Shaft_t *S, int j )
{
    ContourTarget_t T;
    double t;
    int i,k,pos = 0,neg = 0;

    if ( !ContourPrepare( &RTElements[j],&T ) ) return FALSE;   /* keep it */

    for( i=0; i<S->A.NV+S->B.NV; i++ )
    {
        double *x = (i < S->A.NV) ? S->A.V[i] : S->B.V[i-S->A.NV];

        t = 0.0;
        for( k=0; k<3; k++ ) t += T.N[k]*(x[k]-T.V[0][k]);

        if ( t >  S->Tol ) pos++;
        if ( t < -S->Tol ) neg++;
    }

    return !( pos && neg );
}

static int Walk( Shaft_t *S, VolumeBounds_t *Volume, int *List, int n, int MaxN )
{
    int i,j;

    if ( !Volume ) return n;
    if ( !BoxInShaft(S,&Volume->BBox) ) return n;

    if ( Volume->Left )
    {
        n = Walk( S,Volume->Left, List,n,MaxN );
        if ( n < 0 ) return n;                   /* overflow: stop, do not
                                                    feed -1 back as a count */
        return Walk( S,Volume->Right,List,n,MaxN );
    }

    for( i=0; i<Volume->n; i++ )
    {
        j = Volume->Elements[i];

        if ( j < 0 || j >= RTElementNof )    continue;
        if ( !BoxInShaft(S,&RTElementBBox[j]) ) continue;
        if ( CannotBlock(S,j) )              continue;

        /* the tree can reach the same element from more than one leaf */
        for( j=0; j<n; j++ ) if ( List[j] == Volume->Elements[i] ) break;
        if ( j < n ) continue;

        if ( n >= MaxN ) return -1;          /* caller must fall back */
        List[n++] = Volume->Elements[i];
    }
    return n;
}

/*
 * Collect the candidate blockers of the shaft into List.  Returns the count,
 * or -1 when there are more than MaxN of them.
 */
int ShaftCandidates( Shaft_t *S, int *List, int MaxN )
{
    if ( !VolumeBounds || !RTElementBBox ) return 0;
    return Walk( S,VolumeBounds,List,0,MaxN );
}



