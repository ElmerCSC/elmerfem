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

/*******************************************************************************

Utilities for triangle elements.

Juha Ruokolainen/CSC - 23 Aug 1995

*******************************************************************************/

extern double ShapeFunctionMatrix3[3][3];


/*******************************************************************************

Subdivide a triangle element to two parts in longer of the parameters.

23 Aug 1995

*******************************************************************************/
void TriangleSubdivide( Geometry_t *Geometry, int SubLev,int Where )
{
     Triangle_t *LeftTriangle, *RightTriangle;
     double X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,XN,YN,ZN,L1,L2,L3;

     int i,j;

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
     LeftTriangle = Geometry->Left->Triangle = (Triangle_t *)calloc(sizeof(Triangle_t),1);

     Geometry->Right = (Geometry_t *)calloc(sizeof(Geometry_t),1);
     RightTriangle = Geometry->Right->Triangle = (Triangle_t *)calloc(sizeof(Triangle_t),1);

     Geometry->Left->GeometryType = Geometry->Right->GeometryType = GEOMETRY_TRIANGLE;

     X0 = TriangleValue( 0.0,0.0,Geometry->Triangle->PolyFactors[0] );
     Y0 = TriangleValue( 0.0,0.0,Geometry->Triangle->PolyFactors[1] );
     Z0 = TriangleValue( 0.0,0.0,Geometry->Triangle->PolyFactors[2] );

     X1 = TriangleValue( 1.0,0.0,Geometry->Triangle->PolyFactors[0] );
     Y1 = TriangleValue( 1.0,0.0,Geometry->Triangle->PolyFactors[1] );
     Z1 = TriangleValue( 1.0,0.0,Geometry->Triangle->PolyFactors[2] );

     X2 = TriangleValue( 0.0,1.0,Geometry->Triangle->PolyFactors[0] );
     Y2 = TriangleValue( 0.0,1.0,Geometry->Triangle->PolyFactors[1] );
     Z2 = TriangleValue( 0.0,1.0,Geometry->Triangle->PolyFactors[2] );

     L1 = (X1-X0)*(X1-X0) + (Y1-Y0)*(Y1-Y0) + (Z1-Z0)*(Z1-Z0);
     L2 = (X2-X0)*(X2-X0) + (Y2-Y0)*(Y2-Y0) + (Z2-Z0)*(Z2-Z0);
     L3 = (X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1) + (Z2-Z1)*(Z2-Z1);

     if ( L1 > L2  && L1 > L3 )
     {
        XN = TriangleValue( 0.5,0.0,Geometry->Triangle->PolyFactors[0] );
        YN = TriangleValue( 0.5,0.0,Geometry->Triangle->PolyFactors[1] );
        ZN = TriangleValue( 0.5,0.0,Geometry->Triangle->PolyFactors[2] );

        LeftTriangle->PolyFactors[0][0]  = X0;
        LeftTriangle->PolyFactors[1][0]  = Y0;
        LeftTriangle->PolyFactors[2][0]  = Z0;

        LeftTriangle->PolyFactors[0][1]  = XN - X0;
        LeftTriangle->PolyFactors[1][1]  = YN - Y0;
        LeftTriangle->PolyFactors[2][1]  = ZN - Z0;

        LeftTriangle->PolyFactors[0][2]  = X2 - X0;
        LeftTriangle->PolyFactors[1][2]  = Y2 - Y0;
        LeftTriangle->PolyFactors[2][2]  = Z2 - Z0;

        RightTriangle->PolyFactors[0][0] = XN;
        RightTriangle->PolyFactors[1][0] = YN;
        RightTriangle->PolyFactors[2][0] = ZN;

        RightTriangle->PolyFactors[0][1] = X1 - XN;
        RightTriangle->PolyFactors[1][1] = Y1 - YN;
        RightTriangle->PolyFactors[2][1] = Z1 - ZN;

        RightTriangle->PolyFactors[0][2] = X2 - XN;
        RightTriangle->PolyFactors[1][2] = Y2 - YN;
        RightTriangle->PolyFactors[2][2] = Z2 - ZN;
     } else if ( L2 > L3 )
     {
        XN = TriangleValue( 0.0,0.5,Geometry->Triangle->PolyFactors[0] );
        YN = TriangleValue( 0.0,0.5,Geometry->Triangle->PolyFactors[1] );
        ZN = TriangleValue( 0.0,0.5,Geometry->Triangle->PolyFactors[2] );

        LeftTriangle->PolyFactors[0][0]  = X0;
        LeftTriangle->PolyFactors[1][0]  = Y0;
        LeftTriangle->PolyFactors[2][0]  = Z0;

        LeftTriangle->PolyFactors[0][1]  = X1 - X0;
        LeftTriangle->PolyFactors[1][1]  = Y1 - Y0;
        LeftTriangle->PolyFactors[2][1]  = Z1 - Z0;
 
        LeftTriangle->PolyFactors[0][2]  = XN - X0;
        LeftTriangle->PolyFactors[1][2]  = YN - Y0;
        LeftTriangle->PolyFactors[2][2]  = ZN - Z0;

        RightTriangle->PolyFactors[0][0] = XN;
        RightTriangle->PolyFactors[1][0] = YN;
        RightTriangle->PolyFactors[2][0] = ZN;

        RightTriangle->PolyFactors[0][1] = X1 - XN;
        RightTriangle->PolyFactors[1][1] = Y1 - YN;
        RightTriangle->PolyFactors[2][1] = Z1 - ZN;

        RightTriangle->PolyFactors[0][2] = X2 - XN;
        RightTriangle->PolyFactors[1][2] = Y2 - YN;
        RightTriangle->PolyFactors[2][2] = Z2 - ZN;
     } else 
     {
        XN = TriangleValue( 0.5,0.5,Geometry->Triangle->PolyFactors[0] );
        YN = TriangleValue( 0.5,0.5,Geometry->Triangle->PolyFactors[1] );
        ZN = TriangleValue( 0.5,0.5,Geometry->Triangle->PolyFactors[2] );

        LeftTriangle->PolyFactors[0][0]  = X0;
        LeftTriangle->PolyFactors[1][0]  = Y0;
        LeftTriangle->PolyFactors[2][0]  = Z0;

        LeftTriangle->PolyFactors[0][1]  = X1 - X0;
        LeftTriangle->PolyFactors[1][1]  = Y1 - Y0;
        LeftTriangle->PolyFactors[2][1]  = Z1 - Z0;

        LeftTriangle->PolyFactors[0][2]  = XN - X0;
        LeftTriangle->PolyFactors[1][2]  = YN - Y0;
        LeftTriangle->PolyFactors[2][2]  = ZN - Z0;

        RightTriangle->PolyFactors[0][0] = X0;
        RightTriangle->PolyFactors[1][0] = Y0;
        RightTriangle->PolyFactors[2][0] = Z0;

        RightTriangle->PolyFactors[0][1] = XN - X0;
        RightTriangle->PolyFactors[1][1] = YN - Y0;
        RightTriangle->PolyFactors[2][1] = ZN - Z0;

        RightTriangle->PolyFactors[0][2] = X2 - X0;
        RightTriangle->PolyFactors[1][2] = Y2 - Y0;
        RightTriangle->PolyFactors[2][2] = Z2 - Z0;
     } 

     for( i=0; i<3; i++ )
     for( j=0; j<3; j++ )
     {
        LeftTriangle->PolyFactors[i+3][j] = RightTriangle->PolyFactors[i+3][j] =
                 Geometry->Triangle->PolyFactors[i+3][j];
     }
}

/*******************************************************************************

Compute element of (iso)line for (bi)linear polynomial

23 Aug 1995

*******************************************************************************/
double TriangleEofL(double U,double V,double *X,double *Y,double *Z,int UnotV)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = TrianglePartialU(U,V,X);
    dXdV = TrianglePartialV(U,V,X);

    dYdU = TrianglePartialU(U,V,Y);
    dYdV = TrianglePartialV(U,V,Y);

    dZdU = TrianglePartialU(U,V,Z);
    dZdV = TrianglePartialV(U,V,Z);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    if ( UnotV )
        return sqrt(Auu);
    else
        return sqrt(Avv);
}

/*******************************************************************************

Compute element of area for a triangle surface element.

23 Aug 1995

*******************************************************************************/
double TriangleEofA(double U,double V,double *X,double *Y,double *Z)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = TrianglePartialU(U,V,X);
    dXdV = TrianglePartialV(U,V,X);

    dYdU = TrianglePartialU(U,V,Y);
    dYdV = TrianglePartialV(U,V,Y);

    dZdU = TrianglePartialU(U,V,Z);
    dZdV = TrianglePartialV(U,V,Z);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    detA = Auu*Avv - Auv*Auv;

    return sqrt(detA);
}

/*******************************************************************************

Compute isoline length of a triangle element in given parameter.

23 Aug 1995

*******************************************************************************/
double TriangleLength( Geometry_t *Geometry,int UnotV )
{
    double *X,*Y,*Z,EofL,Length;

    int i,j=0;

    X = Geometry->Triangle->PolyFactors[0];
    Y = Geometry->Triangle->PolyFactors[1];
    Z = Geometry->Triangle->PolyFactors[2];

    Length = 0.0;
    for( j=0; j<2; j++ )
    for( i=0; i<N_Integ1d; i++ )
    {
        if ( UnotV )
          EofL = TriangleEofL(U_Integ1d[i],(double)j,X,Y,Z,1);
        else
          EofL = TriangleEofL((double)j,U_Integ1d[i],X,Y,Z,0);

        Length += S_Integ1d[i]*EofL;
    }

    return Length;
}

/*******************************************************************************

Compute area of a triangle element.

23 Aug 1995

*******************************************************************************/
double TriangleArea( Geometry_t *Geometry )
{
    double *X,*Y,*Z,EofA,Area;

    int i;

    X = Geometry->Triangle->PolyFactors[0];
    Y = Geometry->Triangle->PolyFactors[1];
    Z = Geometry->Triangle->PolyFactors[2];

    Area = 0.0;
    for( i=0; i<N_Integ3; i++ )
    {
        EofA = TriangleEofA(U_Integ3[i],V_Integ3[i],X,Y,Z);
        Area += S_Integ3[i]*EofA;
    }

    return Area;
}

/*******************************************************************************

Compute differential area to area viewfactor for triangle surface elements by
direct numerical integration.

24 Aug 1995

*******************************************************************************/
double TriangleIntegrateDiffToArea( Geometry_t *GB,
   double FX,double FY,double FZ,double NFX,double NFY,double NFZ)
{
    double DX,DY,DZ,NTX,NTY,NTZ,U,V;
    double F,R,Rs,cosA,cosB;

    double *BX  = GB->Triangle->PolyFactors[0];
    double *BY  = GB->Triangle->PolyFactors[1];
    double *BZ  = GB->Triangle->PolyFactors[2];

    double *NBX = GB->Triangle->PolyFactors[3];
    double *NBY = GB->Triangle->PolyFactors[4];
    double *NBZ = GB->Triangle->PolyFactors[5];

    int i,j;


    Rs = NFX*NFX + NFY*NFY + NFZ*NFZ;
    if ( Rs != 0 && ABS(1-R) > 1.0E-8 )
    {
        R = 1.0/sqrt(Rs);
        NFX *= Rs;
        NFY *= Rs;
        NFZ *= Rs;
    }

    U = V = 1.0 / 3.0;
    NTX = TriangleValue(U,V,NBX);
    NTY = TriangleValue(U,V,NBY);
    NTZ = TriangleValue(U,V,NBZ);

    R = NTX*NTX + NTY*NTY + NTZ*NTZ;
    if ( ABS(1-R) > 1.0E-8 )
    {
        R = 1.0/sqrt(R);
        NTX *= R;
        NTY *= R;
        NTZ *= R;
    }

    F  = 0.0;
    cosA = 1;
    for( i=0; i<N_Integ3; i++ )
    {
       U = U_Integ3[i];
       V = V_Integ3[i];
        
       DX  = TriangleValue(U,V,BX) - FX;
       DY  = TriangleValue(U,V,BY) - FY;
       DZ  = TriangleValue(U,V,BZ) - FZ;
       R = sqrt(DX*DX + DY*DY + DZ*DZ);

       if ( Rs != 0 ) {
         cosA = (DX*NFX + DY*NFY + DZ*NFZ) / R;
         if ( cosA < 1.0E-8 ) continue;
       }

       cosB = (-DX*NTX - DY*NTY - DZ*NTZ) / R;
       if ( cosB < 1.0E-8 ) continue;

       F += 2*GB->Area*cosA*cosB*S_Integ3[i] / (R*R);
    }
    return F;
}

/*******************************************************************************

Compute area to area viewfactor for triangle surface elements by subdivision
and direct numerical integration when the differential viewfactors match given
magnitude criterion or areas of the elements are small enough. Blocking of the
view between the elements is resolved by ray tracing.

24 Aug 1995

*******************************************************************************/
void TriangleComputeViewFactors(Geometry_t *GA,Geometry_t *GB,
                   int LevelA,int LevelB)
{
    double F,Fa,Fb,EA,PI=2*acos(0.0);
    double FX,FY,FZ,DX,DY,DZ,U,V,Hit;

    double *AX = GA->Triangle->PolyFactors[0];
    double *AY = GA->Triangle->PolyFactors[1];
    double *AZ = GA->Triangle->PolyFactors[2];

    double *NX = GA->Triangle->PolyFactors[3];
    double *NY = GA->Triangle->PolyFactors[4];
    double *NZ = GA->Triangle->PolyFactors[5];

    int i,j;

    if ( GA == GB ) return; /* Linear triangles can't see each other... */

    U = V  = 1.0/3.0;
    FX = TriangleValue( U,V,AX );
    FY = TriangleValue( U,V,AY );
    FZ = TriangleValue( U,V,AZ );

    DX = TriangleValue( U,V,NX );
    DY = TriangleValue( U,V,NY );
    DZ = TriangleValue( U,V,NZ );
    Fa = Fb = (*IntegrateDiffToArea[GB->GeometryType])(GB,FX,FY,FZ,DX,DY,DZ);

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

       Fb = TriangleIntegrateDiffToArea( GA,FX,FY,FZ,DX,DY,DZ );
    }

    if ( Fa < 1.0E-10 && Fb < 1.0E-10 ) return;

    if ( (Fa<FactorEPS || GB->Area<AreaEPS) && (Fb<FactorEPS || GA->Area<AreaEPS) )
    {
        GeometryList_t *Link;

        Hit = Nrays;
        for( i=0; i<Nrays; i++ )
        {
            U = drand48(); V=drand48();
            while( U+V>1.0 ) { U = drand48(); V = drand48(); }

            FX = TriangleValue(U,V,AX);
            FY = TriangleValue(U,V,AY);
            FZ = TriangleValue(U,V,AZ);

            U = drand48(); V=drand48();
            if ( GB->GeometryType == GEOMETRY_TRIANGLE )
              while( U+V>1.0 ) { U = drand48(); V = drand48(); }

            DX = FunctionValue(GB,U,V,0) - FX;
            DY = FunctionValue(GB,U,V,1) - FY;
            DZ = FunctionValue(GB,U,V,2) - FZ;

            Hit -= RayHitGeometry( FX,FY,FZ,DX,DY,DZ );
        }

        if ( Hit == 0 ) return;

        if ( Hit == Nrays || ((Fa<FactorEPS/2 || GB->Area < AreaEPS/2) &&
                 (Fb<FactorEPS/2 || GA->Area < AreaEPS/2) ))
        {
            Hit /= Nrays;

            F = 0.0;
            for( i=0; i<N_Integ3; i++ )
            {
                U = U_Integ3[i];
                V = V_Integ3[i];

                FX = TriangleValue( U,V,AX );
                FY = TriangleValue( U,V,AY );
                FZ = TriangleValue( U,V,AZ );

                DX = TriangleValue( U,V,NX );
                DY = TriangleValue( U,V,NY );
                DZ = TriangleValue( U,V,NZ );


                EA = 2*GA->Area;           
                F += S_Integ3[i]*EA*(*IntegrateDiffToArea[GB->GeometryType])(GB,FX,FY,FZ,DX,DY,DZ);
            }

            F = Hit*F / PI;
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
        if ( !GB->Left ) TriangleSubdivide( GB,LevelB,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;
            GB->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area )
            {
                GB->Left->Area  = TriangleArea(GB->Left);
                GB->Right->Area = TriangleArea(GB->Right);
            }
        }

        TriangleComputeViewFactors( GB->Left, GB->Left, LevelB+1,LevelB+1 );
        TriangleComputeViewFactors( GB->Right,GB->Right,LevelB+1,LevelB+1 );
        TriangleComputeViewFactors( GB->Left, GB->Right,LevelB+1,LevelB+1 );
    }
    else if ( Fa > Fb )
    {
        if ( !GB->Left ) (*Subdivide[GB->GeometryType])( GB, LevelB,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {    
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;
            GB->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area )
            {
                GB->Left->Area  = (*AreaCompute[GB->GeometryType])(GB->Left);
                GB->Right->Area = (*AreaCompute[GB->GeometryType])(GB->Right);
            }
        }

        TriangleComputeViewFactors( GA,GB->Left,LevelA,LevelB+1 );
        TriangleComputeViewFactors( GA,GB->Right,LevelA,LevelB+1 );
    } else
    {
        if ( !GA->Left ) TriangleSubdivide( GA, LevelA,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = TriangleArea(GA->Left);
                GA->Right->Area = TriangleArea(GA->Right);
            }
        }

        TriangleComputeViewFactors( GA->Left,GB,LevelA+1,LevelB );
        TriangleComputeViewFactors( GA->Right,GB,LevelA+1,LevelB );
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
TriangleComputeRadiatorFactors (Geometry_t * GA, double dx, double dy,
			      double dz, int LevelA)
{
  double R, FX, FY, FZ, GX, GY, GZ, U, V, Hit;
  double F, Fa, Fb, EA, PI = 2 * acos (0.0);

  double *X = GA->Triangle->PolyFactors[0];
  double *Y = GA->Triangle->PolyFactors[1];
  double *Z = GA->Triangle->PolyFactors[2];

  int i, j;

  if (LevelA & 1)
    {
      Fa = 0;
      Fb = 0;
      goto subdivide;
    }

    Fa = Fb = TriangleIntegrateDiffToArea( GA,dx,dy,dz,0.0,0.0,0.0);
    if ( Fa < 1.0e-10 ) return;

    if ( Fa<FactorEPS || GA->Area<AreaEPS )
    {
       GeometryList_t *Link;

       Hit = Nrays;
       for( i=0; i<Nrays; i++ )
       {
          U = drand48();
          V = drand48();
          while( U+V>1.0 ) { U = drand48(); V = drand48(); }

          FX = TriangleValue(U,V,X);
          FY = TriangleValue(U,V,Y);
          FZ = TriangleValue(U,V,Z);

          GX = dx - FX;
          GY = dy - FY;
          GZ = dz - FZ;

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

        if ( !GA->Left ) TriangleSubdivide( GA, LevelA,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = TriangleArea(GA->Left);
                GA->Right->Area = TriangleArea(GA->Right);
            }
        }

        TriangleComputeRadiatorFactors( GA->Left,dx,dy,dz,LevelA+1 );
        TriangleComputeRadiatorFactors( GA->Right,dx,dy,dz,LevelA+1 );
}
