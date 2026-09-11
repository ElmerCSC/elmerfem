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

Visibility by clipping, in place of ray casting.

Occlusion is not a property of a patch pair: the blocked part of the target
is different for every point of the source, which is what a penumbra is.  So
the visible part is resolved per source integration point x:

    F(A->B) = 1/A_A  Int_A [ Int_{visible part of B seen from x} ... dA ] dA

The inner integral is the exact contour formula of ContourInteg.c, and the
visible part is found by removing shadows:

    pieces = { B }
    for each candidate occluder o:
        w      = shadow wedge of o seen from x
        pieces = union over p in pieces of ( p \ w )
    result = sum over pieces of ContourPoly( piece )

The shadow of a convex occluder is the intersection of half spaces -- one per
occluder edge, spanned by that edge and x, plus the occluder's own plane to
keep only what is beyond it -- and a convex region minus an intersection of
half spaces is a handful of convex pieces.  So no general polygon boolean is
needed and no hole is ever represented: the contour formula is additive over
the fragments, so a ring of pieces around a hole gives the right answer
without anyone naming the hole.  Sliver fragments are harmless for the same
reason, they contribute nothing and disappear, which is the main robustness
argument for doing it this way rather than with a boolean library.

Note that the shadow is *not* built by projecting the occluder onto the plane
of B.  That form needs a division by the distance between the occluder and
the plane through x parallel to B, which blows up for an occluder touching
that plane -- and in a real mesh the source patch's own neighbours sit exactly
there, so their shadows come out infinite and wipe the target out.  The wedge
of half spaces has no such division and no near plane special case.

Compared with the ray cast, which scales the *unoccluded* kernel by the
fraction of rays that got through, this weights each surviving piece by its
own cos cos / r^2.  The two differ by a few percent whenever the blocked part
of the target has a different kernel from the visible part, which is exactly
what a blocker sitting in the middle of the shaft does.

Juha Ruokolainen/CSC

*******************************************************************************/

#define CL_EPS     1.0e-30
#define CL_AREPS   1.0e-14
#define CL_MAXW    (CL_MAXV+1)

/*
 * A shadow: the region { p : N_i.p >= D_i for every i }.
 */
typedef struct
{
    int    n;
    double N[CL_MAXW][3];
    double D[CL_MAXW];
} Wedge_t;

/*
 * Clip a convex polygon by the half space { p : side*(N.p - D) >= 0 }.
 * Sutherland-Hodgman, so the result stays convex.  Returns -1 rather than
 * truncating when the result would not fit: a piece gains a vertex per half
 * space it meets, so with many occluders it can outgrow the buffer, and
 * silently dropping vertices makes the fragment too big, which shows up as
 * shadow that never got subtracted.
 */
static int ClipHalf( ClipPoly_t *In, double *Nrm, double D, double side,
                       ClipPoly_t *Out )
{
    double dc,dn,t;
    int i,j,k,m = 0;

    if ( In->n < 3 ) { Out->n = 0; return 0; }

    j  = In->n-1;
    dc = side*( Nrm[0]*In->V[j][0] + Nrm[1]*In->V[j][1] + Nrm[2]*In->V[j][2] - D );

    for( i=0; i<In->n; i++ )
    {
        dn = side*( Nrm[0]*In->V[i][0] + Nrm[1]*In->V[i][1] + Nrm[2]*In->V[i][2] - D );

        if ( dn >= 0.0 )
        {
            if ( dc < 0.0 )
            {
                if ( m >= CL_MAXV ) { Out->n = 0; return -1; }
                t = dc/(dc-dn);
                for( k=0; k<3; k++ )
                   Out->V[m][k] = In->V[j][k] + t*(In->V[i][k]-In->V[j][k]);
                m++;
            }
            if ( m >= CL_MAXV ) { Out->n = 0; return -1; }
            for( k=0; k<3; k++ ) Out->V[m][k] = In->V[i][k];
            m++;
        }
        else if ( dc >= 0.0 )
        {
            if ( m >= CL_MAXV ) { Out->n = 0; return -1; }
            t = dc/(dc-dn);
            for( k=0; k<3; k++ )
               Out->V[m][k] = In->V[j][k] + t*(In->V[i][k]-In->V[j][k]);
            m++;
        }
        dc = dn;
        j  = i;
    }
    Out->n = m;
    return m;
}

static double PolyArea( ClipPoly_t *P )
{
    double G[3];
    int i,j;

    if ( P->n < 3 ) return 0.0;

    G[0] = G[1] = G[2] = 0.0;
    for( i=0; i<P->n; i++ )
    {
        j = (i+1) % P->n;
        G[0] += P->V[i][1]*P->V[j][2] - P->V[i][2]*P->V[j][1];
        G[1] += P->V[i][2]*P->V[j][0] - P->V[i][0]*P->V[j][2];
        G[2] += P->V[i][0]*P->V[j][1] - P->V[i][1]*P->V[j][0];
    }
    return 0.5*sqrt( G[0]*G[0] + G[1]*G[1] + G[2]*G[2] );
}

/*
 * The shadow wedge of a convex occluder seen from x: one half space per
 * occluder edge, spanned by that edge and x, plus the occluder's own plane
 * oriented away from x so that only what lies beyond the occluder is in
 * shadow.  Returns 0 when the occluder is edge on to x and shadows nothing.
 */
static int ShadowWedge( ContourTarget_t *Occ, double *X, Wedge_t *W )
{
    double a[3],b[3],n[3],c[3],s,L,fx;
    int i,j,k;

    W->n = 0;

    /* the occluder's own plane, oriented away from x */
    fx = Occ->N[0]*(X[0]-Occ->V[0][0]) + Occ->N[1]*(X[1]-Occ->V[0][1])
       + Occ->N[2]*(X[2]-Occ->V[0][2]);
    if ( ABS(fx) < CL_EPS ) return 0;             /* x in the occluder plane */
    s = (fx > 0.0) ? -1.0 : 1.0;

    for( k=0; k<3; k++ ) W->N[0][k] = s*Occ->N[k];
    W->D[0] = s*( Occ->N[0]*Occ->V[0][0] + Occ->N[1]*Occ->V[0][1]
                + Occ->N[2]*Occ->V[0][2] );
    W->n = 1;

    /* centroid, to fix which side of each edge plane is the shadowed one */
    for( k=0; k<3; k++ ) c[k] = 0.0;
    for( i=0; i<Occ->NV; i++ )
      for( k=0; k<3; k++ ) c[k] += Occ->V[i][k]/Occ->NV;

    for( i=0; i<Occ->NV; i++ )
    {
        j = (i+1) % Occ->NV;

        for( k=0; k<3; k++ ) { a[k] = Occ->V[i][k]-X[k]; b[k] = Occ->V[j][k]-X[k]; }

        n[0] = a[1]*b[2] - a[2]*b[1];
        n[1] = a[2]*b[0] - a[0]*b[2];
        n[2] = a[0]*b[1] - a[1]*b[0];

        L = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
        if ( L < CL_EPS ) return 0;               /* x colinear with an edge */
        for( k=0; k<3; k++ ) n[k] /= L;

        s = n[0]*(c[0]-X[0]) + n[1]*(c[1]-X[1]) + n[2]*(c[2]-X[2]);
        if ( ABS(s) < CL_EPS ) return 0;          /* degenerate occluder */
        if ( s < 0.0 ) for( k=0; k<3; k++ ) n[k] = -n[k];

        if ( W->n >= CL_MAXW ) return 0;
        for( k=0; k<3; k++ ) W->N[W->n][k] = n[k];
        W->D[W->n] = n[0]*X[0] + n[1]*X[1] + n[2]*X[2];
        W->n++;
    }
    return W->n;
}

/*
 * P minus the wedge, as disjoint convex pieces.  Appends to Out and returns
 * the new count, or -1 when Out would overflow.  The leftover, P inside the
 * wedge, is the shadowed part and is dropped.
 */
static int ConvexDifference( ClipPoly_t *P, Wedge_t *W,
                               ClipPoly_t *Out, int nout, int maxout )
{
    ClipPoly_t rest = *P, keep, drop;
    double d,dmin,dmax,tol,scale;
    int i,j;

    for( i=0; i<W->n; i++ )
    {
        /*
         * Classify first.  Clipping blindly is wrong when the polygon lies
         * *in* the plane: every distance is then zero, both half spaces keep
         * everything, and the piece is emitted as visible while also staying
         * behind to be peeled again -- so it gets duplicated once per such
         * plane.  An occluder coplanar with the target does exactly that.
         */
        dmin =  1.0e30; dmax = -1.0e30; scale = 0.0;
        for( j=0; j<rest.n; j++ )
        {
            d = W->N[i][0]*rest.V[j][0] + W->N[i][1]*rest.V[j][1]
              + W->N[i][2]*rest.V[j][2];
            if ( ABS(d) > scale ) scale = ABS(d);
            d -= W->D[i];
            if ( d < dmin ) dmin = d;
            if ( d > dmax ) dmax = d;
        }
        tol = 1.0e-12*( scale + ABS(W->D[i]) ) + CL_EPS;

        /*
         * The plane contains the polygon, so it separates nothing: every
         * point is in the closed half space, and both clips would keep
         * everything -- which would emit the piece as visible *and* keep it
         * for further peeling, duplicating it once per such plane.  Treat it
         * as inside and move on.  The case where that would be wrong, an
         * occluder coplanar with the target, is removed by ClipVisible
         * before the wedge is ever built.
         */
        if ( dmax-dmin <= tol ) continue;

        /* wholly on the shadowed side: nothing to emit, keep peeling */
        if ( dmin >= -tol ) continue;

        /* wholly outside the shadow: the piece survives entire */
        if ( dmax <= tol )
        {
            if ( nout >= maxout ) return -1;
            Out[nout++] = rest;
            return nout;
        }

        /* genuinely cut by this plane */
        if ( ClipHalf(&rest,W->N[i],W->D[i],-1.0,&drop) < 0 ) return -1;
        if ( drop.n >= 3 && PolyArea(&drop) > CL_AREPS )
        {
            if ( nout >= maxout ) return -1;
            Out[nout++] = drop;
        }

        if ( ClipHalf(&rest,W->N[i],W->D[i], 1.0,&keep) < 0 ) return -1;
        rest = keep;
        if ( rest.n < 3 ) break;
    }
    return nout;
}

/*
 * Numerical counterpart of ContourPoly: the same integral,
 *
 *     Int cosA cosB / R^2 dA
 *
 * over one clipped fragment, by Gauss quadrature instead of Lambert's contour
 * formula.  A fragment is always convex -- ConvexDifference only ever emits
 * convex pieces -- so a fan from vertex 0 triangulates it validly, and the
 * existing triangle rule integrates each fan triangle.  The rule is the one
 * "Viewfactor Triangle Integration Points" selects.
 *
 * This is the fourth corner of the 2x2, and it exists as a cross check rather
 * than as a path anyone should prefer: it shares the clipped geometry with the
 * closed form exactly, so a disagreement between the two is in the inner
 * integration and cannot be in the shadowing.
 *
 * Note it does not clip the fragment to the front half space of dA the way
 * ContourPoly must.  Quadrature does not need it -- a point behind dA simply
 * has cosA <= 0 and drops out, exactly as in the unclipped Gauss path -- but
 * that also means the two differ at grazing incidence by the quadrature's own
 * error there, which is the largest part of what this comparison measures.
 */
static int PolyIntegNumeric( double (*V)[3], int nv, double *N,
      double FX, double FY, double FZ, double NFX, double NFY, double NFZ,
        double *Result )
{
    double F = 0.0;
    int i,j,k;

    *Result = 0.0;
    if ( nv < 3 ) return TRUE;

    /* cosB has a constant sign over a planar fragment, as in ContourPoly */
    if ( (FX-V[0][0])*N[0] + (FY-V[0][1])*N[1] + (FZ-V[0][2])*N[2] <= 0.0 )
        return TRUE;

    for( i=1; i<nv-1; i++ )
    {
        double E1[3],E2[3],CR[3],A2;

        for( k=0; k<3; k++ )
        {
            E1[k] = V[i][k]   - V[0][k];
            E2[k] = V[i+1][k] - V[0][k];
        }

        CR[0] = E1[1]*E2[2] - E1[2]*E2[1];
        CR[1] = E1[2]*E2[0] - E1[0]*E2[2];
        CR[2] = E1[0]*E2[1] - E1[1]*E2[0];

        /*
         * The magnitude of the cross product is twice the fan triangle's
         * area, which is the Jacobian of the map to the reference triangle.
         * The rule's weights sum to 1/2, so this is the same normalisation
         * the geometry routines use as 2*Area.
         */
        A2 = sqrt( CR[0]*CR[0] + CR[1]*CR[1] + CR[2]*CR[2] );
        if ( A2 <= 0.0 ) continue;

        for( j=0; j<N_Integ3; j++ )
        {
            double D[3],R,R2,cosA,cosB,U = U_Integ3[j],W = V_Integ3[j];

            for( k=0; k<3; k++ ) D[k] = V[0][k] + U*E1[k] + W*E2[k];
            D[0] -= FX; D[1] -= FY; D[2] -= FZ;

            R2 = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
            if ( R2 <= 0.0 ) continue;
            R = sqrt(R2);

            cosA = ( D[0]*NFX + D[1]*NFY + D[2]*NFZ ) / R;
            if ( cosA < 1.0e-8 ) continue;

            cosB = -( D[0]*N[0] + D[1]*N[1] + D[2]*N[2] ) / R;
            if ( cosB < 1.0e-8 ) continue;

            F += A2*cosA*cosB*S_Integ3[j] / R2;
        }
    }

    *Result = F;
    return TRUE;
}

/*
 * The kernel: view factor integrand from dA at X with unit normal NF to the
 * visible part of the prepared target, with the candidate occluders removed.
 *
 * Returns FALSE when the fragment count runs away, and the caller must fall
 * back to ray casting for the pair.
 */
int ClipVisible( ContourTarget_t *Tgt, int *Cand, int nc,
                   double *X, double *NF, double *Result )
{
    ClipPoly_t A[CL_MAXPIECE],B[CL_MAXPIECE];
    ContourTarget_t Occ;
    Wedge_t W;
    double F,f;
    int na,nb,i,k,c;

    *Result = 0.0;

    A[0].n = Tgt->NV;
    for( i=0; i<Tgt->NV; i++ )
      for( k=0; k<3; k++ ) A[0].V[i][k] = Tgt->V[i][k];
    na = 1;

    for( c=0; c<nc; c++ )
    {
        double dmin = 1.0e30, dmax = -1.0e30, d, scale = 0.0, tol;

        if ( !ContourPrepare( &RTElements[Cand[c]],&Occ ) ) return FALSE;

        /*
         * An occluder lying in the plane of the target cannot shadow it: it
         * has no extent in the one direction that would matter.  Left in, it
         * would make the wedge degenerate and the target would be judged by
         * the edge planes alone, which can swallow it whole.
         */
        for( i=0; i<Occ.NV; i++ )
        {
            d = 0.0;
            for( k=0; k<3; k++ ) d += Tgt->N[k]*(Occ.V[i][k]-Tgt->V[0][k]);
            if ( ABS(Occ.V[i][0])+ABS(Occ.V[i][1])+ABS(Occ.V[i][2]) > scale )
               scale = ABS(Occ.V[i][0])+ABS(Occ.V[i][1])+ABS(Occ.V[i][2]);
            if ( d < dmin ) dmin = d;
            if ( d > dmax ) dmax = d;
        }
        tol = 1.0e-10*scale + CL_EPS;
        if ( dmax-dmin <= tol && ABS(dmin) <= tol ) continue;

        if ( !ShadowWedge( &Occ,X,&W ) ) continue;

        nb = 0;
        for( i=0; i<na; i++ )
        {
            nb = ConvexDifference( &A[i],&W,B,nb,CL_MAXPIECE );
            if ( nb < 0 ) return FALSE;
        }
        for( i=0; i<nb; i++ ) A[i] = B[i];
        na = nb;

        if ( na == 0 ) return TRUE;                   /* fully shadowed */
    }

    /*
     * Integrate what survived.  Both branches see the identical fragments,
     * which is the whole point of keeping the numerical one: it separates an
     * error in the inner integral from an error in the shadowing.
     */
    F = 0.0;
    for( i=0; i<na; i++ )
    {
        if ( ClosedFormInteg )
        {
            if ( !ContourPoly( A[i].V,A[i].n,Tgt->N, X[0],X[1],X[2],
                                 NF[0],NF[1],NF[2], &f ) ) return FALSE;
        } else
        {
            if ( !PolyIntegNumeric( A[i].V,A[i].n,Tgt->N, X[0],X[1],X[2],
                                      NF[0],NF[1],NF[2], &f ) ) return FALSE;
        }
        F += f;
    }

    *Result = F;
    return TRUE;
}

