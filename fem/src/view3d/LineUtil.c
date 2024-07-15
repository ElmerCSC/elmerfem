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

static int same;
/*******************************************************************************

Utilities for linear elements.

Juha Ruokolainen/CSC - 01 Sep 2008

*******************************************************************************/

extern double ShapeFunctionMatrix2[2][2];

/*******************************************************************************

Subdivide a linear element to two parts.

23 Aug 1995

*******************************************************************************/
void LinearSubdivideHalfU(double *I,double *L,double *R)
{
    double t,x[2];
    int i,j;

    x[0] = LinearValue(0.0,I);
    x[1] = LinearValue(0.5,I);

    for( i=0; i<2; i++ )
    for( j=0; j<2; j++ ) L[i] += ShapeFunctionMatrix2[j][i]*x[j];

    x[0] = LinearValue(0.5,I);
    x[1] = LinearValue(1.0,I);

    for( i=0; i<2; i++ )
    for( j=0; j<2; j++ ) R[i] += ShapeFunctionMatrix2[j][i]*x[j];

    return;
}


/*******************************************************************************

Subdivide a linear element to two parts

23 Aug 1995

*******************************************************************************/
void LinearSubdivide( Geometry_t *Geometry, int SubLev,int Where )
{
     Linear_t *LeftLinear, *RightLinear;
     BBox_t *BBox;
     double ULength,VLength;

     int j;

     Geometry->Left = (Geometry_t *)calloc(sizeof(Geometry_t),1);
     LeftLinear = Geometry->Left->Linear = (Linear_t *)calloc(sizeof(Linear_t),1);

     Geometry->Right = (Geometry_t *)calloc(sizeof(Geometry_t),1);
     RightLinear = Geometry->Right->Linear = (Linear_t *)calloc(sizeof(Linear_t),1);

     Geometry->Left->GeometryType = Geometry->Right->GeometryType = GEOMETRY_LINE;

     for( j=0; j<6; j++ )
     {
       LinearSubdivideHalfU( Geometry->Linear->PolyFactors[j],
          LeftLinear->PolyFactors[j],RightLinear->PolyFactors[j] );
     }
}

/*******************************************************************************

Compute element of (iso)line for linear polynomial

23 Aug 1995

*******************************************************************************/
double LinearEofL(double U,double V,double *X,double *Y,double *Z,int UnotV)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = LinearPartialU(U,X);
    dYdU = LinearPartialU(U,Y);
    return sqrt(dXdU*dXdU+dYdU*dYdU);
}

/*******************************************************************************

Compute element of area for a linear surface element.

23 Aug 1995

*******************************************************************************/
double LinearEofA(double U,double V,double *X,double *Y,double *Z)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = LinearPartialU(U,X);
    dYdU = LinearPartialU(U,Y);
    detA = dXdU*dXdU + dYdU*dYdU;

    return sqrt(detA);
}

/*******************************************************************************

Compute isoline length of a linear element in given parameter.

23 Aug 1995

*******************************************************************************/
double LinearLength( Geometry_t *Geometry,int UnotV )
{
    double *X,*Y,*Z,EofL,Length;

    int i,j=0;

    X = Geometry->Linear->PolyFactors[0];
    Y = Geometry->Linear->PolyFactors[1];
    Z = Geometry->Linear->PolyFactors[2];

    Length = 0.0;
    for( j=0; j<2; j++ )
    for( i=0; i<N_Integ1d; i++ )
    {
        if ( UnotV )
            EofL = LinearEofL(U_Integ1d[i],(double)j,X,Y,Z,1);
        else
            EofL = LinearEofL((double)j,U_Integ1d[i],X,Y,Z,0);

        Length += S_Integ1d[i]*EofL;
    }

    return Length;
}

/*******************************************************************************

Compute area of a linear element.

23 Aug 1995

*******************************************************************************/
double LinearArea( Geometry_t *Geometry )
{
    double *X,*Y,*Z,EofA,Area;

    int i;

    X = Geometry->Linear->PolyFactors[0];
    Y = Geometry->Linear->PolyFactors[1];
    Z = Geometry->Linear->PolyFactors[2];

    Area = 0.0;
    for( i=0; i<N_Integ1d; i++ )
    {
        EofA = LinearEofA(U_Integ1d[i],0.0,X,Y,Z);
        Area += S_Integ1d[i]*EofA;
    }

    return Area;
}


/*******************************************************************************

Compute differential area to area viewfactor for linear surface elements by
direct numerical integration.

24 Aug 1995

*******************************************************************************/
double LinearIntegrateDiffToArea( Geometry_t *GB,
    double FX,double FY,double FZ, double NFX,double NFY,double NFZ )
{
    double DX,DY,DZ,NTX,NTY,NTZ,U,V;
    double F,R,Rs,cosA,cosB,EA,d;

    double *BX  = GB->Linear->PolyFactors[0];
    double *BY  = GB->Linear->PolyFactors[1];
    double *BZ  = GB->Linear->PolyFactors[2];

    double *NBX = GB->Linear->PolyFactors[3];
    double *NBY = GB->Linear->PolyFactors[4];
    double *NBZ = GB->Linear->PolyFactors[5];

    int i,j;

    Rs = NFX*NFX + NFY*NFY;
    if ( Rs != 0 && ABS(1-Rs)>1.0e-9 )
    {
       R = 1.0/sqrt(Rs);
       NFX *= Rs;
       NFY *= Rs;
    }

    F  = 0.0;
    cosA = 1.0;

    for( i=0; i<N_Integ1d; i++ )
    {
       U = U_Integ1d[i];
       V = 0.0;
        
       DX = LinearValue(U,BX) - FX;
       DY = LinearValue(U,BY) - FY;

       R = sqrt(DX*DX + DY*DY);

       if ( Rs>0 ) {
         cosA = (DX*NFX + DY*NFY) / R;
         if ( cosA < 1.0e-9 ) continue;
       }

       NTX = LinearValue(U,NBX);
       NTY = LinearValue(U,NBY);


       d = NTX*NTX + NTY*NTY;
       if ( ABS(1-R)>1.0e-9 ) {
          d = 1.0 / sqrt(d);
          NTX *= d; NTY *= d;
       }

       cosB = (-DX*NTX - DY*NTY)/R;
       if ( cosB < 1.0e-9 ) continue;

       EA = LinearEofA( U,V,BX,BY,BZ );
       F += EA*S_Integ1d[i]*cosA*cosB / (2*R);
    }

    return F;
}


double LinearViewFactor( Geometry_t *GA, Geometry_t *GB )
{
  double *AX  = GA->Linear->PolyFactors[0];
  double *AY  = GA->Linear->PolyFactors[1];
  double *NAX = GA->Linear->PolyFactors[3];
  double *NAY = GA->Linear->PolyFactors[4];

  double *BX  = GB->Linear->PolyFactors[0];
  double *BY  = GB->Linear->PolyFactors[1];
  double *NBX = GB->Linear->PolyFactors[3];
  double *NBY = GB->Linear->PolyFactors[4];

  double fact,d, NTX, NTY, NFX, NFY, cosA, cosB, DX, DY, FX, FY, eps=1e-8;

  FX = LinearValue(0.5,AX);
  FY = LinearValue(0.5,AY);

  NFX = LinearValue(0.5,NAX);
  NFY = LinearValue(0.5,NAY);
  d = NFX*NFX + NFY*NFY;
  if ( ABS(1-d) > eps ) {
    d = sqrt(d); NTX/=d; NTY/=d; 
  }

  NTX = LinearValue(0.5,NBX);
  NTY = LinearValue(0.5,NBY);
  d = NTX*NTX + NTY*NTY;
  if ( ABS(1-d) > eps ) {
    d = sqrt(d); NTX/=d; NTY/=d; 
  }

  DX = LinearValue(0.5,BX)-FX;
  DY = LinearValue(0.5,BY)-FY;
  d = DX*DX + DY*DY;
  if ( d != 0.0 ) { d=1./sqrt(d); DX*=d; DY*=d; }

  cosA =  NFX*DX + NFY*DY;
  cosB = -NTX*DX - NTY*DY;
  if ( cosA < eps || cosB < eps ) return 0.0;

  DX = LinearValue(0.0,BX)-LinearValue(0.0,AX);
  DY = LinearValue(0.0,BY)-LinearValue(0.0,AY);
  fact = sqrt(DX*DX + DY*DY);

  DX = LinearValue(1.0,BX)-LinearValue(1.0,AX);
  DY = LinearValue(1.0,BY)-LinearValue(1.0,AY);
  fact += sqrt(DX*DX + DY*DY);

  DX = LinearValue(0.0,BX)-LinearValue(1.0,AX);
  DY = LinearValue(0.0,BY)-LinearValue(1.0,AY);
  fact -= sqrt(DX*DX + DY*DY);

  DX = LinearValue(1.0,BX)-LinearValue(0.0,AX);
  DY = LinearValue(1.0,BY)-LinearValue(0.0,AY);
  fact -= sqrt(DX*DX + DY*DY);
  fact /= 2;

  return ABS(fact);

}




/*******************************************************************************

Compute area to area viewfactor for linear surface elements by subdivision
and direct numerical integration when the differential viewfactors match given
magnitude criterion or areas of the elements are small enough. Blocking of the
view between the elements is resolved by ray tracing.

24 Aug 1995

*******************************************************************************/
void LinearComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int LevelA,int LevelB)
{
    double R,FX,FY,FZ,DX,DY,DZ,U,V,Hit;
    double F,Fa,Fb,EA,PI=2*acos(0.0);

    double *X  = GA->Linear->PolyFactors[0];
    double *Y  = GA->Linear->PolyFactors[1];
    double *Z  = GA->Linear->PolyFactors[2];

    double *NX = GA->Linear->PolyFactors[3];
    double *NY = GA->Linear->PolyFactors[4];
    double *NZ = GA->Linear->PolyFactors[5];

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

    U = 0.5; V=0.0;
    FX = LinearValue( U,X );
    FY = LinearValue( U,Y );
    FZ = 0.0;

    DX = LinearValue( U,NX );
    DY = LinearValue( U,NY );
    DZ = 0.0;

    Fa = Fb = (*IntegrateDiffToArea[GB->GeometryType])( GB,FX,FY,FZ,DX,DY,DZ );

    if ( GA != GB ) 
    {
       U = 0.5; V=0.0;
       if ( GB->GeometryType == GEOMETRY_TRIANGLE ) U = V = 1.0/3.0;

       FX = FunctionValue( GB,U,V,0 );
       FY = FunctionValue( GB,U,V,1 );
       FZ = 0.0;

       DX = FunctionValue( GB,U,V,3 );
       DY = FunctionValue( GB,U,V,4 );
       DZ = 0.0;

       Fb = LinearIntegrateDiffToArea( GA,FX,FY,FZ,DX,DY,DZ );
    }

    if ( Fa < 1.0e-10 && Fb < 1.0e-10 ) return;

    if ( (Fa<FactorEPS || GB->Area<AreaEPS) && (Fb<FactorEPS || GA->Area<AreaEPS) )
    {
       GeometryList_t *Link;

       Hit = Nrays;
       for( i=0; i<Nrays; i++ )
       {
          U = drand48(); V = drand48();

          FX = LinearValue(U,X);
          FY = LinearValue(U,Y);
          FZ = 0.0;

           U = drand48(); V = drand48();
           if ( GB->GeometryType == GEOMETRY_TRIANGLE )
               while( U+V>1 ) { U=drand48(); V=drand48(); }

           DX = FunctionValue(GB,U,V,0)-FX;
           DY = FunctionValue(GB,U,V,1)-FY;
           DZ = 0.0;

           if ( RayHitGeometry( FX,FY,FZ,DX,DY,DZ ) ) Hit-=1.0;
        }

        if ( Hit == 0 ) return;

        if ( Hit == Nrays || ( ((Fa<FactorEPS/2 || GB->Area < AreaEPS/2) &&
                 (Fb<FactorEPS/2 || GA->Area < AreaEPS/2) ) ))
        {
            F = 0.0;
#if 1
            if ( GA != GB ) F = LinearViewFactor(GA,GB);
#else
            if ( GA != GB )
              for( i=0; i<N_Integ1d; i++ )
              {
                  U = U_Integ1d[i]; V = 0.0;
  
                  FX = LinearValue(U,X);
                  FY = LinearValue(U,Y);
                  FZ = 0.0;
  
                  DX = LinearValue(U,NX);
                  DY = LinearValue(U,NY);
                  DZ = 0.0;
  
                  EA = LinearEofA(U,V,X,Y,Z);
                  F += S_Integ1d[i] * EA *
                      (*IntegrateDiffToArea[GB->GeometryType])(GB,FX,FY,FZ,DX,DY,DZ);
              }
#endif

            F = Hit*F/(1.0*Nrays);
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
        if ( !GB->Left ) LinearSubdivide( GB, LevelB,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;

            GB->Left->Flags   |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags  |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area )
            {
                GB->Left->Area  = LinearArea(GB->Left);
                GB->Right->Area = LinearArea(GB->Right);
            }
        }

        LinearComputeViewFactors( GB->Left, GB->Left, LevelA+1,LevelB+1 );
        LinearComputeViewFactors( GB->Right,GB->Right,LevelA+1,LevelB+1 );
        LinearComputeViewFactors( GB->Left, GB->Right,LevelA+1,LevelB+1 );
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

        LinearComputeViewFactors( GA,GB->Left,LevelA,LevelB+1 );
        LinearComputeViewFactors( GA,GB->Right,LevelA,LevelB+1 );
    } else
    {
        if ( !GA->Left ) LinearSubdivide( GA, LevelA,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = LinearArea(GA->Left);
                GA->Right->Area = LinearArea(GA->Right);
            }
        }

        LinearComputeViewFactors( GA->Left,GB,LevelA+1,LevelB );
        LinearComputeViewFactors( GA->Right,GB,LevelA+1,LevelB );
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
LinearComputeRadiatorFactors (Geometry_t * GA, double dx, double dy,
         double dz, int LevelA)
{
  double R, FX, FY, FZ, GX, GY, GZ, U, V, Hit;
  double F, Fa, Fb, EA, PI = 2 * acos (0.0);

  double *X = GA->Linear->PolyFactors[0];
  double *Y = GA->Linear->PolyFactors[1];
  double *Z = GA->Linear->PolyFactors[2];

  int i, j;

  if (LevelA & 1)
    {
      Fa = 0;
      Fb = 0;
        goto subdivide;
    }

    Fa = Fb = LinearIntegrateDiffToArea( GA,dx,dy,dz,0.0,0.0,0.0);

    if ( Fa < 1.0e-10 && Fb < 1.0e-10 ) return;

    if ( Fa<FactorEPS || GA->Area<AreaEPS )
    {
       GeometryList_t *Link;

       Hit = Nrays;
       for( i=0; i<Nrays; i++ )
       {
          U = drand48();

          FX = LinearValue(U,X);
          FY = LinearValue(U,Y);
          FZ = 0.0;

          GX = dx - FX;
          GY = dy - FY;
          GZ = 0.0;

           if ( RayHitGeometry( FX,FY,FZ,GX,GY,GZ ) ) Hit-=1.0;
        }

        if ( Hit == 0 ) return;

        if ( Hit == Nrays || Fa<FactorEPS/2 || GA->Area < AreaEPS/2 )
        {
            F = Hit*Fa/(1.0*Nrays);
            Fa = Fb = F / GA->Area / PI;

            Link = (GeometryList_t *)calloc(sizeof(GeometryList_t),1);
            Link->Next = GA->Link;
            GA->Link = Link;

            Link->Entry = GA;
            Link->ViewFactor = Fa;

            return;
        }
    }

subdivide:

        if ( !GA->Left ) LinearSubdivide( GA, LevelA,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = LinearArea(GA->Left);
                GA->Right->Area = LinearArea(GA->Right);
            }
        }

        LinearComputeRadiatorFactors( GA->Left,dx,dy,dz,LevelA+1 );
        LinearComputeRadiatorFactors( GA->Right,dx,dy,dz,LevelA+1 );
}
