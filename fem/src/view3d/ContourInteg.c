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

Closed form evaluation of the inner ("to area") view factor integral

    I = Int_P  cosA * cosB / R^2  dA

over a *planar* polygon P, replacing the Gauss quadrature used by the
*IntegrateDiffToArea() routines.  This is Lambert's contour formula

    I = 1/2 Sum_edges  beta_i * ( nF . unit(A_i x A_i+1) ),

    A_i = P_i - F,   beta_i = angle( A_i, A_i+1 ),

which is exact for a planar polygon and costs one evaluation instead of a
quadrature loop.  The gain is in the near field: with the patches close
together the integrand is nearly singular and the quadrature needs many
points, while the contour formula does not care.

The semantics of the numerical routines are reproduced exactly:

  - the source is one sided, cosA > 0 only: the polygon is clipped to the
    half space in front of dA, which is what dropping the cosA < 0 Gauss
    points does in the limit;

  - the target is one sided, cosB > 0 only: for a planar polygon with a
    constant normal cosB has a constant sign over the whole patch, namely
    that of (F-P_0).NT, so a single test replaces the per point one.

Not applicable, and the caller must fall back to quadrature, when

  - the patch is not planar to within CI_PLANAR_TOL (warped bilinear quad),
  - the supplied patch normal is not the geometric one (CI_NORMAL_TOL),

because in both cases the quadrature is integrating a different surface,
or a different integrand, than the contour formula would.

Juha Ruokolainen/CSC

*******************************************************************************/

/* must exceed the largest polygon handed in, plus the vertex a half space
 * clip can add: ContourPoly is called on clipping fragments, which are
 * CL_MAXV long.  Truncating instead would close the contour through the
 * wrong vertices and return a value with no relation to the polygon. */
#define CI_MAXV        (CL_MAXV+8)
#define CI_PLANAR_TOL  1.0e-6      /* out of plane, relative to patch size   */
#define CI_NORMAL_TOL  1.0e-6      /* 1-|Ng.NT|, geometric vs given normal   */
#define CI_EPS         1.0e-30

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Diagnostics.  Counted once per element pair in ContourPrepare, never in the
 * per integration point path, and written to a padded per thread slot so the
 * element pair loop never contends on a shared cache line.
 */
void ContourCountInit()
{
#ifdef _OPENMP
    ClosedFormNThreads = omp_get_max_threads();
#else
    ClosedFormNThreads = 1;
#endif
    if ( ClosedFormNThreads < 1 ) ClosedFormNThreads = 1;

    if ( ClosedFormCount ) free( ClosedFormCount );
    ClosedFormCount = (long *)calloc( 2*ClosedFormNThreads*CF_STRIDE, sizeof(long) );
}

static void ContourCount( int Used )
{
    int tid = 0;

    if ( !ClosedFormInteg || !ClosedFormCount ) return;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if ( tid < 0 || tid >= ClosedFormNThreads ) return;

    ClosedFormCount[ (2*tid + (Used?0:1))*CF_STRIDE ]++;
}

void ContourCountSum( long *Hits, long *Miss )
{
    int i;

    *Hits = *Miss = 0;
    if ( !ClosedFormCount ) return;

    for( i=0; i<ClosedFormNThreads; i++ )
    {
        *Hits += ClosedFormCount[ (2*i+0)*CF_STRIDE ];
        *Miss += ClosedFormCount[ (2*i+1)*CF_STRIDE ];
    }
}

/*
 * Clip a polygon by the half space { p : (p-F).N >= 0 }.  Sutherland-Hodgman;
 * the input is convex so the output is a single convex loop.  Returns the
 * number of output vertices, or -1 if it would not fit -- truncating instead
 * would close the contour through the wrong vertices and hand back a number
 * with no relation to the polygon.
 */
static int ClipToFrontHalf( double (*V)[3], int n,
      double FX, double FY, double FZ, double NX, double NY, double NZ,
        double C[CI_MAXV][3] )
{
    double dc,dn,t;
    int i,j,k,m = 0;

    j  = n-1;
    dc = (V[j][0]-FX)*NX + (V[j][1]-FY)*NY + (V[j][2]-FZ)*NZ;

    for( i=0; i<n; i++ )
    {
        dn = (V[i][0]-FX)*NX + (V[i][1]-FY)*NY + (V[i][2]-FZ)*NZ;

        if ( dn >= 0.0 )
        {
            if ( dc < 0.0 )
            {
                if ( m >= CI_MAXV ) return -1;
                t = dc/(dc-dn);
                for( k=0; k<3; k++ ) C[m][k] = V[j][k] + t*(V[i][k]-V[j][k]);
                m++;
            }
            if ( m >= CI_MAXV ) return -1;
            for( k=0; k<3; k++ ) C[m][k] = V[i][k];
            m++;
        }
        else if ( dc >= 0.0 )
        {
            if ( m >= CI_MAXV ) return -1;
            t = dc/(dc-dn);
            for( k=0; k<3; k++ ) C[m][k] = V[j][k] + t*(V[i][k]-V[j][k]);
            m++;
        }
        dc = dn;
        j  = i;
    }
    return m;
}

/*
 * Prepare a target patch: corners, unit geometric normal, and the checks that
 * decide whether the closed form may be used at all.  Called once per element
 * pair; the result is valid for every source integration point on that pair.
 *
 * Returns FALSE when the caller must use the quadrature instead:
 *   - geometry type without a polygon boundary,
 *   - patch not planar to CI_PLANAR_TOL (a warped bilinear quad),
 *   - given patch normal is not the geometric one (CI_NORMAL_TOL), because
 *     then the quadrature is integrating a different integrand than we would.
 */
int ContourPrepare( Geometry_t *GB, ContourTarget_t *T )
{
    static const double QU[4] = { 0.0, 1.0, 1.0, 0.0 };
    static const double QV[4] = { 0.0, 0.0, 1.0, 1.0 };
    static const double TU[3] = { 0.0, 1.0, 0.0 };
    static const double TV[3] = { 0.0, 0.0, 1.0 };

    double *BX,*BY,*BZ,*NBX,*NBY,*NBZ;
    double GX,GY,GZ,NTX,NTY,NTZ,A2,D,Size;
    int i,j;

    switch( GB->GeometryType )
    {
    case GEOMETRY_TRIANGLE:
        BX  = GB->Triangle->PolyFactors[0];
        BY  = GB->Triangle->PolyFactors[1];
        BZ  = GB->Triangle->PolyFactors[2];
        NBX = GB->Triangle->PolyFactors[3];
        NBY = GB->Triangle->PolyFactors[4];
        NBZ = GB->Triangle->PolyFactors[5];

        T->NV = 3;
        for( i=0; i<3; i++ )
        {
            T->V[i][0] = TriangleValue(TU[i],TV[i],BX);
            T->V[i][1] = TriangleValue(TU[i],TV[i],BY);
            T->V[i][2] = TriangleValue(TU[i],TV[i],BZ);
        }
        NTX = TriangleValue(1.0/3.0,1.0/3.0,NBX);
        NTY = TriangleValue(1.0/3.0,1.0/3.0,NBY);
        NTZ = TriangleValue(1.0/3.0,1.0/3.0,NBZ);
        break;

    case GEOMETRY_BILINEAR:
        BX  = GB->BiLinear->PolyFactors[0];
        BY  = GB->BiLinear->PolyFactors[1];
        BZ  = GB->BiLinear->PolyFactors[2];
        NBX = GB->BiLinear->PolyFactors[3];
        NBY = GB->BiLinear->PolyFactors[4];
        NBZ = GB->BiLinear->PolyFactors[5];

        T->NV = 4;
        for( i=0; i<4; i++ )
        {
            T->V[i][0] = BiLinearValue(QU[i],QV[i],BX);
            T->V[i][1] = BiLinearValue(QU[i],QV[i],BY);
            T->V[i][2] = BiLinearValue(QU[i],QV[i],BZ);
        }
        /* the patch normal is constant: one normal per element */
        NTX = BiLinearValue(0.5,0.5,NBX);
        NTY = BiLinearValue(0.5,0.5,NBY);
        NTZ = BiLinearValue(0.5,0.5,NBZ);
        break;

    default:
        return FALSE;
    }

    /* no source normal available (radiator path uses cosA == 1) */
    D = NTX*NTX + NTY*NTY + NTZ*NTZ;
    if ( D <= 0.0 ) { ContourCount(FALSE); return FALSE; }
    D = 1.0/sqrt(D);
    NTX *= D; NTY *= D; NTZ *= D;

    /* Newell normal; its length is twice the polygon area */
    GX = GY = GZ = 0.0;
    for( i=0; i<T->NV; i++ )
    {
        j = (i+1) % T->NV;
        GX += T->V[i][1]*T->V[j][2] - T->V[i][2]*T->V[j][1];
        GY += T->V[i][2]*T->V[j][0] - T->V[i][0]*T->V[j][2];
        GZ += T->V[i][0]*T->V[j][1] - T->V[i][1]*T->V[j][0];
    }
    A2 = sqrt( GX*GX + GY*GY + GZ*GZ );
    if ( A2 < CI_EPS ) { ContourCount(FALSE); return FALSE; }

    GX /= A2; GY /= A2; GZ /= A2;
    Size = sqrt(A2);

    /* planarity: a triangle always is, so only quads pay for this */
    for( i=3; i<T->NV; i++ )
    {
        D = (T->V[i][0]-T->V[0][0])*GX
          + (T->V[i][1]-T->V[0][1])*GY
          + (T->V[i][2]-T->V[0][2])*GZ;
        if ( ABS(D) > CI_PLANAR_TOL*Size ) { ContourCount(FALSE); return FALSE; }
    }

    D = GX*NTX + GY*NTY + GZ*NTZ;
    if ( ABS(D) < 1.0 - CI_NORMAL_TOL ) { ContourCount(FALSE); return FALSE; }

    /* orient the stored normal the way the caller's normal points, so the
       one sided target test in ContourEvaluate is a single dot product */
    if ( D < 0.0 ) { GX = -GX; GY = -GY; GZ = -GZ; }

    T->N[0] = GX; T->N[1] = GY; T->N[2] = GZ;

    ContourCount(TRUE);
    return TRUE;
}

/*
 * Evaluate  Int_P cosA*cosB/R^2 dA  for an arbitrary planar convex polygon
 * lying in the plane of a prepared patch, seen from (FX,FY,FZ) with *unit*
 * normal (NFX,NFY,NFZ).  N is the patch normal, oriented as ContourPrepare
 * leaves it.  This is the form the clipping path needs: it is called once
 * per visible fragment, all of which share the target patch's plane.
 */
int ContourPoly( double (*V)[3], int nv, double *N,
      double FX, double FY, double FZ, double NFX, double NFY, double NFZ,
        double *Result )
{
    double C[CI_MAXV][3];
    double (*P)[3];
    double AX,AY,AZ,BX,BY,BZ,RX,RY,RZ,S,D,beta,F;
    int i,j,n;

    *Result = 0.0;
    if ( nv < 3 ) return TRUE;

    /*
     * cosB has a constant sign over a planar patch: the sign of (F-P0).N.
     * A back facing target has every Gauss point dropped, so the answer is 0.
     */
    if ( (FX-V[0][0])*N[0] + (FY-V[0][1])*N[1] + (FZ-V[0][2])*N[2] <= 0.0 )
        return TRUE;

    /*
     * The source is one sided too.  The common case is that the whole polygon
     * lies in front of dA, and then no clipping and no copying is needed.
     */
    P = V;
    n = nv;

    S = 0.0;
    for( i=0; i<nv; i++ )
    {
        D = (V[i][0]-FX)*NFX + (V[i][1]-FY)*NFY + (V[i][2]-FZ)*NFZ;
        if ( D < S ) S = D;
    }
    if ( S < 0.0 )
    {
        n = ClipToFrontHalf( V,nv, FX,FY,FZ, NFX,NFY,NFZ, C );
        if ( n < 0 ) return FALSE;              /* caller falls back */
        if ( n < 3 ) return TRUE;
        P = C;
    }

    /*
     * Lambert's contour formula,
     *     Int cosA cosB / R^2 dA = 1/2 Sum beta_i ( nF . unit(A_i x A_i+1) ).
     */
    F = 0.0;
    for( i=0; i<n; i++ )
    {
        j = (i+1) % n;

        AX = P[i][0]-FX; AY = P[i][1]-FY; AZ = P[i][2]-FZ;
        BX = P[j][0]-FX; BY = P[j][1]-FY; BZ = P[j][2]-FZ;

        RX = AY*BZ - AZ*BY;
        RY = AZ*BX - AX*BZ;
        RZ = AX*BY - AY*BX;

        S = sqrt( RX*RX + RY*RY + RZ*RZ );
        if ( S < CI_EPS ) continue;          /* colinear with F: no solid angle */

        beta = atan2( S, AX*BX + AY*BY + AZ*BZ );
        F += beta*( NFX*RX + NFY*RY + NFZ*RZ )/S;
    }

    /* winding of the polygon is arbitrary; the integral is positive by
       construction once the one sided tests above have passed */
    *Result = 0.5*ABS(F);

    return TRUE;
}

/* The whole prepared patch, unoccluded. */
int ContourEvaluate( ContourTarget_t *T,
      double FX, double FY, double FZ, double NFX, double NFY, double NFZ,
        double *Result )
{
    return ContourPoly( T->V,T->NV,T->N, FX,FY,FZ, NFX,NFY,NFZ, Result );
}
