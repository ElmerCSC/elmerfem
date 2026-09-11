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

#include <ViewFactors.h>

static int same;
/*******************************************************************************

Utilities for bilinear elements.

Juha Ruokolainen/CSC - 23 Aug 1995

*******************************************************************************/

extern double ShapeFunctionMatrix4[4][4];

/*******************************************************************************

Subdivide a bilinear element to two parts in first of the parameters.

23 Aug 1995

*******************************************************************************/
void BiLinearSubdivideHalfU(double *I,double *L,double *R)
{
    double t,x[4];
    int i,j;

    x[0] = BiLinearValue(0.0,0.0,I);
    x[1] = BiLinearValue(0.5,0.0,I);
    x[2] = BiLinearValue(0.5,1.0,I);
    x[3] = BiLinearValue(0.0,1.0,I);

    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ ) L[i] += ShapeFunctionMatrix4[j][i]*x[j];

    x[0] = BiLinearValue(0.5,0.0,I);
    x[1] = BiLinearValue(1.0,0.0,I);
    x[2] = BiLinearValue(1.0,1.0,I);
    x[3] = BiLinearValue(0.5,1.0,I);

    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ ) R[i] += ShapeFunctionMatrix4[j][i]*x[j];

    return;
}

/*******************************************************************************

Subdivide a bilinear element to two parts in second of the parameters.

23 Aug 1995

*******************************************************************************/
void BiLinearSubdivideHalfV(double *I,double *L,double *R)
{
    double t,x[4];
    int i,j;

    x[0] = BiLinearValue(0.0,0.0,I);
    x[1] = BiLinearValue(1.0,0.0,I);
    x[2] = BiLinearValue(1.0,0.5,I);
    x[3] = BiLinearValue(0.0,0.5,I);

    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ ) L[i] += ShapeFunctionMatrix4[j][i]*x[j];

    x[0] = BiLinearValue(0.0,0.5,I);
    x[1] = BiLinearValue(1.0,0.5,I);
    x[2] = BiLinearValue(1.0,1.0,I);
    x[3] = BiLinearValue(0.0,1.0,I);

    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ ) R[i] += ShapeFunctionMatrix4[j][i]*x[j];

    return;
}

/*******************************************************************************

Subdivide a bilinear element to two parts in longer of the parameters.

23 Aug 1995

*******************************************************************************/
void BiLinearSubdivide( Geometry_t *Geometry, int SubLev,int Where )
{
     BiLinear_t *LeftLinear, *RightLinear;
     BBox_t *BBox;
     double ULength,VLength;

     int j;

#if 0
{
    volatile static int l=0;
    fprintf( stderr, "add new node (total: %d,level: %d,",2*l,SubLev );
    if ( Where ) fprintf( stderr, " from: Integrate\n" );
    else fprintf( stderr, "from: raytrace\n" );
    l++;
}
#endif

     Geometry->Left = (Geometry_t *)calloc(sizeof(Geometry_t),1);
     LeftLinear = Geometry->Left->BiLinear = (BiLinear_t *)calloc(sizeof(BiLinear_t),1);

     Geometry->Right = (Geometry_t *)calloc(sizeof(Geometry_t),1);
     RightLinear = Geometry->Right->BiLinear = (BiLinear_t *)calloc(sizeof(BiLinear_t),1);

     Geometry->Left->GeometryType = Geometry->Right->GeometryType = GEOMETRY_BILINEAR;

     ULength = BiLinearLength(Geometry,1);
     VLength = BiLinearLength(Geometry,0);

     if ( ULength>VLength )
     {
         for( j=0; j<6; j++ )
         {
             BiLinearSubdivideHalfU( Geometry->BiLinear->PolyFactors[j],
              LeftLinear->PolyFactors[j],RightLinear->PolyFactors[j] );
         }
     } else
     {
         for( j=0; j<6; j++ )
         {
             BiLinearSubdivideHalfV( Geometry->BiLinear->PolyFactors[j],
              LeftLinear->PolyFactors[j],RightLinear->PolyFactors[j] );
         }
     }
}

/*******************************************************************************

Compute element of (iso)line for (bi)linear polynomial

23 Aug 1995

*******************************************************************************/
double BiLinearEofL(double U,double V,double *X,double *Y,double *Z,int UnotV)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = BiLinearPartialU(U,V,X);
    dXdV = BiLinearPartialV(U,V,X);

    dYdU = BiLinearPartialU(U,V,Y);
    dYdV = BiLinearPartialV(U,V,Y);

    dZdU = BiLinearPartialU(U,V,Z);
    dZdV = BiLinearPartialV(U,V,Z);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    if ( UnotV )
        return sqrt(Auu);
    else
        return sqrt(Avv);
}

/*******************************************************************************

Compute element of area for a bilinear surface element.

23 Aug 1995

*******************************************************************************/
double BiLinearEofA(double U,double V,double *X,double *Y,double *Z)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = BiLinearPartialU(U,V,X);
    dXdV = BiLinearPartialV(U,V,X);

    dYdU = BiLinearPartialU(U,V,Y);
    dYdV = BiLinearPartialV(U,V,Y);

    dZdU = BiLinearPartialU(U,V,Z);
    dZdV = BiLinearPartialV(U,V,Z);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    detA = Auu*Avv - Auv*Auv;

    return sqrt(detA);
}

/*******************************************************************************

Compute isoline length of a bilinear element in given parameter.

23 Aug 1995

*******************************************************************************/
double BiLinearLength( Geometry_t *Geometry,int UnotV )
{
    double *X,*Y,*Z,EofL,Length;

    int i,j=0;

    X = Geometry->BiLinear->PolyFactors[0];
    Y = Geometry->BiLinear->PolyFactors[1];
    Z = Geometry->BiLinear->PolyFactors[2];

    Length = 0.0;
    for( j=0; j<2; j++ )
    for( i=0; i<N_Integ1d; i++ )
    {
        if ( UnotV )
            EofL = BiLinearEofL(U_Integ1d[i],(double)j,X,Y,Z,1);
        else
            EofL = BiLinearEofL((double)j,U_Integ1d[i],X,Y,Z,0);

        Length += S_Integ1d[i]*EofL;
    }

    return Length;
}

/*******************************************************************************

Compute area of a bilinear element.

23 Aug 1995

*******************************************************************************/
double BiLinearArea( Geometry_t *Geometry )
{
    double *X,*Y,*Z,EofA,Area;

    int i;

    X = Geometry->BiLinear->PolyFactors[0];
    Y = Geometry->BiLinear->PolyFactors[1];
    Z = Geometry->BiLinear->PolyFactors[2];

    Area = 0.0;
    for( i=0; i<N_Integ; i++ )
    {
        EofA = BiLinearEofA(U_Integ[i],V_Integ[i],X,Y,Z);
        Area += S_Integ[i]*EofA;
    }

    return Area;
}

/*******************************************************************************

Compute differential area to area viewfactor for bilinear surface elements by
direct numerical integration.

24 Aug 1995

*******************************************************************************/
double BiLinearIntegrateDiffToArea( Geometry_t *GB, Cylinder_t *Cyl, 
    double FX,double FY,double FZ, double NFX,double NFY,double NFZ )
{
    double DX,DY,DZ,NTX,NTY,NTZ,U,V;
    double F,R,Rs,cosA,cosB,EA;

    double *BX  = GB->BiLinear->PolyFactors[0];
    double *BY  = GB->BiLinear->PolyFactors[1];
    double *BZ  = GB->BiLinear->PolyFactors[2];

    double *NBX = GB->BiLinear->PolyFactors[3];
    double *NBY = GB->BiLinear->PolyFactors[4];
    double *NBZ = GB->BiLinear->PolyFactors[5];

    int i,j;

    R = NFX*NFX + NFY*NFY + NFZ*NFZ;
    if ( !Cyl ) {
      if ( R != 0 ) {
         R = 1.0/sqrt(R);
         NFX *= R;
         NFY *= R;
         NFZ *= R;
      }
    }
    Rs = R;

    /*
     * Closed form inner integral for a planar patch.  This is the one off
     * path (the subdivision estimates); the integration loops below prepare
     * the target once and call ContourEvaluate directly.
     */
    if ( ClosedFormInteg && !Cyl && Rs != 0.0 )
    {
        ContourTarget_t CT;
        double CF;

        if ( ContourPrepare(GB,&CT) &&
               ContourEvaluate(&CT,FX,FY,FZ,NFX,NFY,NFZ,&CF) ) return CF;
    }

    F  = 0.0;
    cosA = 1;
    for( i=0; i<N_Integ; i++ )
    {
       U = U_Integ[i];
       V = V_Integ[i];
        
       NTX = BiLinearValue(U,V,NBX);
       NTY = BiLinearValue(U,V,NBY);
       NTZ = BiLinearValue(U,V,NBZ);
       R = NTX*NTX + NTY*NTY + NTZ*NTZ;

       if ( ABS(1-R)>1.0E-08 )
       {
           R = 1.0/sqrt(R);
           NTX *= R;
           NTY *= R;
           NTZ *= R;
       }

       DX  = BiLinearValue(U,V,BX) - FX;
       DY  = BiLinearValue(U,V,BY) - FY;
       DZ  = BiLinearValue(U,V,BZ) - FZ;
       R = sqrt(DX*DX + DY*DY + DZ*DZ);

       if ( Cyl ) {
         CylinderNormal(FX,FY,FZ,DX,DY,DZ,Cyl,&NFX,&NFY,&NFZ );
       }


       if ( Rs != 0) {
         cosA = (DX*NFX + DY*NFY + DZ*NFZ) / R;
         if ( cosA < 1.0E-8 ) continue;
       }

       cosB = (-DX*NTX - DY*NTY - DZ*NTZ) / R;
       if ( cosB < 1.0E-8 ) continue;


       EA = BiLinearEofA( U,V,BX,BY,BZ );
       F += EA*cosA*cosB*S_Integ[i] / (R*R);
    }

    return F;
}

/*******************************************************************************

Compute area to area viewfactor for bilinear surface elements by subdivision
and direct numerical integration when the differential viewfactors match given
magnitude criterion or areas of the elements are small enough. Blocking of the
view between the elements is resolved by ray tracing.

24 Aug 1995

*******************************************************************************/
void BiLinearComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int LevelA,int LevelB)
{
    double R,FX,FY,FZ,DX,DY,DZ,U,V,W,Hit;
    double F,Fa,Fb,EA,PI=2*acos(0.0);

    double *X  = GA->BiLinear->PolyFactors[0];
    double *Y  = GA->BiLinear->PolyFactors[1];
    double *Z  = GA->BiLinear->PolyFactors[2];

    double *NX = GA->BiLinear->PolyFactors[3];
    double *NY = GA->BiLinear->PolyFactors[4];
    double *NZ = GA->BiLinear->PolyFactors[5];

    int i,j;

    if ( LevelA & 1 ) 
    {
        Fa = 0; Fb = 1;
        goto subdivide;
    }

    if ( (LevelB & 1) && (GB->GeometryType != GEOMETRY_TRIANGLE) ) 
    {
        Fa = 1; Fb = 0;
        goto subdivide;
    }

    U = V = 0.5;
    FX = BiLinearValue( U,V,X );
    FY = BiLinearValue( U,V,Y );
    FZ = BiLinearValue( U,V,Z );

    DX = BiLinearValue( U,V,NX );
    DY = BiLinearValue( U,V,NY );
    DZ = BiLinearValue( U,V,NZ );

    Fa = Fb = (*IntegrateDiffToArea[GB->GeometryType])( GB,NULL,FX,FY,FZ,DX,DY,DZ );

    if ( GA != GB ) 
    {
       U = V = 0.5;
       if ( GB->GeometryType == GEOMETRY_TRIANGLE ) U = V = 1.0/3.0;

       FX = FunctionValue( GB,U,V,0 );
       FY = FunctionValue( GB,U,V,1 );
       FZ = FunctionValue( GB,U,V,2 );

       DX = FunctionValue( GB,U,V,3 );
       DY = FunctionValue( GB,U,V,4 );
       DZ = FunctionValue( GB,U,V,5 );

       Fb = BiLinearIntegrateDiffToArea( GA,NULL,FX,FY,FZ,DX,DY,DZ );
    }

    if ( Fa < 1.0E-10 && Fb < 1.0E-10 ) return;

    if ( (Fa<FactorEPS || GB->Area<AreaEPS) && (Fb<FactorEPS || GA->Area<AreaEPS) )
    {
        GeometryList_t *Link;

        /*
         * The shaft of this pair: the candidate blockers.  Built once and
         * used by everything below -- the statistics, the clipping, and the
         * ray casting, whose rays can only ever meet these same candidates.
         */
        int Clipped = FALSE, HaveShaft = FALSE, nc = -1;
        int Cand[SC_MAXCAND];
        ContourTarget_t TA,TB;
        Shaft_t Sh;

        if ( (ShaftRayCull || ShaftStats || ClipShadows) && GA != GB )
        {
            if ( ContourPrepare(GA,&TA) && ContourPrepare(GB,&TB) &&
                    ShaftInit(&Sh,&TA,&TB) )
            {
                nc = ShaftCandidates( &Sh,Cand,SC_MAXCAND );
                HaveShaft = ( nc >= 0 );
            }
        }

        if ( ShaftStats ) ShaftCountAdd( HaveShaft ? nc : -3 );

        /*
         * Visibility by clipping.  Falls through to the ray casting below
         * whenever it does not apply: a geometry without a polygon boundary,
         * a shaft with too many candidate blockers, or a fragment count that
         * runs away.
         */
        if ( ClipShadows && HaveShaft )
        {
            {
                /*
                 * Accept or subdivide on the same terms as the ray path: an
                 * unobstructed pair outright, a partly blocked one only when
                 * the factors are small.  The outer integral still has a kink
                 * across the penumbra edge, so that criterion still earns its
                 * keep even though the inner one is now exact.
                 */
                if ( nc > 0 && !((Fa<FactorEPS/2 || GB->Area < AreaEPS/2) &&
                                 (Fb<FactorEPS/2 || GA->Area < AreaEPS/2)) )
                    goto subdivide;

                F = 0.0;
                Clipped = TRUE;

                for( i=0; i<N_Integ; i++ )
                {
                    double XP[3],NF[3],DR,DF;

                    U = U_Integ[i];
                    V = V_Integ[i];

                    XP[0] = BiLinearValue(U,V,X);
                    XP[1] = BiLinearValue(U,V,Y);
                    XP[2] = BiLinearValue(U,V,Z);

                    NF[0] = BiLinearValue(U,V,NX);
                    NF[1] = BiLinearValue(U,V,NY);
                    NF[2] = BiLinearValue(U,V,NZ);

                    DR = NF[0]*NF[0] + NF[1]*NF[1] + NF[2]*NF[2];
                    if ( DR <= 0.0 ) { Clipped = FALSE; break; }
                    DR = 1.0/sqrt(DR);
                    NF[0] *= DR; NF[1] *= DR; NF[2] *= DR;

                    if ( !ClipVisible(&TB,Cand,nc,XP,NF,&DF) )
                       { Clipped = FALSE; break; }

                    F += S_Integ[i]*BiLinearEofA(U,V,X,Y,Z)*DF;
                }

                if ( Clipped )
                {
                    if ( F <= 0.0 ) return;         /* fully shadowed */
                    F /= PI;
                    goto store;
                }
            }
        }

        Hit = Nrays;
        for( i=0; i<Nrays; i++ )
        {
            U = vrand();
	    V = vrand();
            FX = BiLinearValue(U,V,X);
            FY = BiLinearValue(U,V,Y);
            FZ = BiLinearValue(U,V,Z);

	    W = U;
	    U=1-V; V=1-W;
            if ( GB->GeometryType == GEOMETRY_TRIANGLE )
                while( U+V>1 ) { U=1-U; V=1-V; }

            DX = FunctionValue(GB,U,V,0) - FX;
            DY = FunctionValue(GB,U,V,1) - FY;
            DZ = FunctionValue(GB,U,V,2) - FZ;

            /* the rays of this pair can only meet its own candidates */
            if ( !ShaftRayCull || !HaveShaft )
                Hit -= RayHitGeometry( FX,FY,FZ,DX,DY,DZ );
            else if ( nc > 0 )
                Hit -= RayHitCandidates( Cand,nc,FX,FY,FZ,DX,DY,DZ );
        }

        /* the cull claimed nothing could block: the rays must agree */
        if ( ShaftStats && HaveShaft && nc == 0 && Hit < Nrays ) ShaftCountAdd( -2 );

        if ( Hit == 0 ) return;

        if ( Hit == Nrays || ((Fa<FactorEPS/2 || GB->Area < AreaEPS/2) &&
                 (Fb<FactorEPS/2 || GA->Area < AreaEPS/2) ))
        {
            /* the target is the same for every source point, so prepare
               the closed form once here instead of inside the loop */
            ContourTarget_t CT;
            int UseCF = ClosedFormInteg && ContourPrepare( GB,&CT );

            F = 0.0;
            for( i=0; i<N_Integ; i++ )
            {
                double DF, DR;

                U = U_Integ[i];
                V = V_Integ[i];

                FX = BiLinearValue(U,V,X);
                FY = BiLinearValue(U,V,Y);
                FZ = BiLinearValue(U,V,Z);

                DX = BiLinearValue(U,V,NX);
                DY = BiLinearValue(U,V,NY);
                DZ = BiLinearValue(U,V,NZ);

                if ( UseCF )
                {
                    DR = DX*DX + DY*DY + DZ*DZ;
                    if ( DR > 0.0 )
                    {
                        DR = 1.0/sqrt(DR);
                        ContourEvaluate( &CT,FX,FY,FZ,DX*DR,DY*DR,DZ*DR,&DF );
                    } else
                        DF = 0.0;
                } else
                    DF = (*IntegrateDiffToArea[GB->GeometryType])( GB,NULL,FX,FY,FZ,DX,DY,DZ );

                EA = BiLinearEofA(U,V,X,Y,Z);
                F += S_Integ[i]*EA*DF;
            }

            F = Hit*F / (PI*Nrays);

store:
            Fb = F / GB->Area;
            Fa = F / GA->Area;

            Link = (GeometryList_t *)calloc(sizeof(GeometryList_t),1);
            Link->Next = GA->Link;
            GA->Link = Link;

            Link->Entry = GB;
            Link->ViewFactor = Fa;

            if ( GA != GB )
            {
                Link = (GeometryList_t *)calloc(sizeof(GeometryList_t),1);
                Link->Next = GB->Link;
                GB->Link = Link;

                Link->Entry = GA;
                Link->ViewFactor = Fb;
            }

            return;
        }
    }

subdivide:

    if ( GA == GB )
    {
        if ( !GB->Left ) BiLinearSubdivide( GB, LevelB,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;

            GB->Left->Flags   |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags  |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area )
            {
                GB->Left->Area  = BiLinearArea(GB->Left);
                GB->Right->Area = BiLinearArea(GB->Right);
            }
        }

        BiLinearComputeViewFactors( GB->Left, GB->Left, LevelA+1,LevelB+1 );
        BiLinearComputeViewFactors( GB->Right,GB->Right,LevelA+1,LevelB+1 );
        BiLinearComputeViewFactors( GB->Left, GB->Right,LevelA+1,LevelB+1 );
    }
    else if ( Fa > Fb )
    {
        if ( !GB->Left ) (*Subdivide[GB->GeometryType])( GB,LevelB,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {    
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;

            GB->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area)
            {
                GB->Left->Area  = (*AreaCompute[GB->GeometryType])(GB->Left);
                GB->Right->Area = (*AreaCompute[GB->GeometryType])(GB->Right);
            }
        }

        BiLinearComputeViewFactors( GA,GB->Left,LevelA,LevelB+1 );
        BiLinearComputeViewFactors( GA,GB->Right,LevelA,LevelB+1 );
    } else
    {
        if ( !GA->Left ) BiLinearSubdivide( GA, LevelA,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = BiLinearArea(GA->Left);
                GA->Right->Area = BiLinearArea(GA->Right);
            }
        }

        BiLinearComputeViewFactors( GA->Left,GB,LevelA+1,LevelB );
        BiLinearComputeViewFactors( GA->Right,GB,LevelA+1,LevelB );
    }
}



/*******************************************************************************

Compute  differential view factor for linear surface elements by subdivision
and direct numerical integration when the differential viewfactors match given
magnitude criterion or areas of the elements are small enough. Blocking of the
view between the elements is resolved by ray tracing.

24 Aug 1995

*******************************************************************************/
void
BiLinearComputeRadiatorFactors (Geometry_t * GA, int LineFlag, double dx, double dy,
    double dz, double nx, double ny, double nz, int LevelA)
{
  double R, FX, FY, FZ, GX, GY, GZ, U, V, Hit;
  double F, Fa, Fb, EA, PI = 2 * acos (0.0), EPS=1e-12;

  double *X = GA->BiLinear->PolyFactors[0];
  double *Y = GA->BiLinear->PolyFactors[1];
  double *Z = GA->BiLinear->PolyFactors[2];

  double *aX = GA->BiLinear->PolyFactors[3];
  double *aY = GA->BiLinear->PolyFactors[4];
  double *aZ = GA->BiLinear->PolyFactors[5];

  Cylinder_t *Cyl = NULL, CylS;

  int i, j, Ident;

  if (LevelA & 1)
    {
      Fa = 0;
      Fb = 0;
      goto subdivide;
    }

    R = nx*nx + ny*ny +nz*nz;
    if (LineFlag &&  R != 0) {
      R = sqrt(R);
      Cyl = &CylS;
      Cyl->Radius = R/25;
      Cyl->CenterPoint.x = (2*dx+nx)/2;
      Cyl->CenterPoint.y = (2*dy+ny)/2;
      Cyl->CenterPoint.z = (2*dz+nz)/2;
      GetMatrixToRotateVectorToZAxis(nx/R,ny/R,nz/R,Cyl->RotationMatrix,&Ident);

      Fa = 0;
      for( i=0; i<N_Integ1d; i++ )
      {
	 U = U_Integ1d[i];
         Fa += S_Integ1d[i]*BiLinearIntegrateDiffToArea(GA,Cyl,dx+U*nx,dy+U*ny,dz+U*nz,nx,ny,nz);
      }
      Fb = Fa;
    } else {
      Fa = Fb = BiLinearIntegrateDiffToArea( GA,Cyl,dx,dy,dz,nx,ny,nz );
    }

    if ( Fa < 1.0e-10 && Fb < 1.0e-10 ) return;

    if ( Fa<FactorEPS || GA->Area<AreaEPS )
    {
       GeometryList_t *Link;

       Hit = Nrays;
       for( i=0; i<Nrays; i++ )
       {
          U = vrand();
          V = vrand();

          FX = BiLinearValue(U,V,X);
          FY = BiLinearValue(U,V,Y);
          FZ = BiLinearValue(U,V,Z);

          GX = dx - FX;
          GY = dy - FY;
          GZ = dz - FZ;
	  if ( Cyl ) {
            U = 1-V;
	    GX += U*nx;
	    GY += U*ny;
	    GZ += U*nz;
	  }

           if ( RayHitGeometry( FX,FY,FZ,GX,GY,GZ ) ) Hit-=1.0;
        }

        if ( Hit == 0 ) return;

        if ( Hit == Nrays || Fa<FactorEPS/2 || GA->Area < AreaEPS/2 )
        {
            F = Hit*Fa/(1.0*Nrays);
            Fa = Fb = F / GA->Area / (4*PI);

            Link = (GeometryList_t *)calloc(sizeof(GeometryList_t),1);
            Link->Next = GA->Link;
            GA->Link = Link;

            Link->Entry = GA;
            Link->ViewFactor = Fa;

            return;
        }
    }

subdivide:

        if ( !GA->Left ) BiLinearSubdivide( GA, LevelA,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = BiLinearArea(GA->Left);
                GA->Right->Area = BiLinearArea(GA->Right);
            }
        }

        BiLinearComputeRadiatorFactors( GA->Left,LineFlag, dx,dy,dz,nx,ny,nz,LevelA+1 );
        BiLinearComputeRadiatorFactors( GA->Right,LineFlag, dx,dy,dz,nx,ny,nz,LevelA+1 );
}

