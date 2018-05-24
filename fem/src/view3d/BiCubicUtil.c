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

/*******************************************************************************

Utilities for bicubic elements.

Juha Ruokolainen/CSC - 23 Aug 1995

*******************************************************************************/

#include <ViewFactors.h>

/*******************************************************************************

Convert monomial basis bicubic polynomial to bezier (bernstein basis) form.

23 Aug 1995

*******************************************************************************/
void BiCubicMonomialToBezier(double *MonomialFactors,double *BezierFactors)
{
     static double CMatrix[4][4] =
     {
         { 1.0,   0.0,     0.0,   0.0 },
         { 1.0, 1.0/3.0,   0.0  , 0.0 },
         { 1.0, 2.0/3.0, 1.0/3.0, 0.0 },
         { 1.0,   1.0,     1.0,   1.0 }
     };

     double s,A[4][4];

     int i,j,k;

/*
     (inv(a')*(inv(a')*r1)')'
     (a*(a*x)')'
*/

     for( i=0; i<4; i++ )
     for( j=0; j<4; j++ )
     {
         s = 0.0;
         for( k=0; k<4; k++ ) s += CMatrix[i][k]*MonomialFactors[4*k+j]; 
         A[i][j] = s;
     }

     for( i=0; i<4; i++ )
     for( j=0; j<4; j++ )
     {
         s = 0.0;
         for( k=0; k<4; k++ ) s += CMatrix[i][k]*A[j][k];
         BezierFactors[4*j+i] = s;
     }
}

/*******************************************************************************

Convert bezier (bernstein basis) form polynomial to monomial form.

23 Aug 1995

*******************************************************************************/
void BiCubicBezierToMonomial(double *MonomialFactors,double *BezierFactors)
{
     static double CMatrix[4][4] =
     {
         { 1.0, -3.0,  3.0, -1.0 },
         { 0.0,  3.0, -6.0,  3.0 },
         { 0.0,  0.0,  3.0, -3.0 },
         { 0.0,  0.0,  0.0,  1.0 }
     };

     double s,A[4][4];

     int i,j,k,n;

     for( i=0; i<4; i++ )
     for( j=0; j<4; j++ )
     {
         s = 0.0;
         for( k=0; k<4; k++ ) s += BezierFactors[4*i+k]*CMatrix[k][j]; 
         A[i][j] = s;
     }

     n = 0;
     for( i=0; i<4; i++ )
     for( j=0; j<4; j++,n++ )
     {
         s = 0.0;
         for( k=0; k<4; k++ ) s += CMatrix[k][i]*A[k][j];
         MonomialFactors[n] = s;
     }
}

/*******************************************************************************

Subdivide bicubic polynomial to two parts in first of the parameters.

23 Aug 1995

*******************************************************************************/
void BiCubicBezierSubdivideHalfU(double *I,double *L,double *R)
{
    double t;
    int j;

    for( j=0; j<4; j++,I+=4,R+=4,L+=4 )
    {
        t = 0.5*(I[1]+I[2]);

        L[0] = I[0];       
        L[1] = 0.5*(I[0]+I[1]);
        L[2] = 0.5*(L[1]+t);

        R[3] = I[3];
        R[2] = 0.5*(I[2]+I[3]);
        R[1] = 0.5*(t+R[2]);

        L[3] = R[0] = 0.5*(L[2]+R[1]);
    }
}

/*******************************************************************************

Subdivide bicubic polynomial to two parts in second of the parameters.

23 Aug 1995

*******************************************************************************/
void BiCubicBezierSubdivideHalfV(double *I,double *L,double *R)
{
    double t;
    int j;

    for( j=0; j<4; j++,I++,R++,L++ )
    {
        t = 0.5*(I[4]+I[8]);

        L[0] = I[0];       
        L[4] = 0.5*(I[0]+I[4]);
        L[8] = 0.5*(L[4]+t);

        R[12] = I[12];
        R[8] = 0.5*(I[8]+I[12]);
        R[4] = 0.5*(t+R[8]);

        L[12] = R[0] = 0.5*(L[8]+R[4]);
    }
}

/*******************************************************************************

Test bicubic polynomial for planarity by testing how well the bezier control
points fulfill the plane equation.

23 Aug 1995

*******************************************************************************/
int BiCubicIsAPlane( double *X,double *Y, double *Z )
{
    double Ax,Ay,Az,Cx,Cy,Cz,Dx,Dy,Dz,D,Nx,Ny,Nz,R;
    int i;

    Ax = X[0];
    Ay = Y[0];
    Az = Z[0];

    for( i=15; i>0; i-- )
    {
        Cx = X[i];
        Cy = Y[i];
        Cz = Z[i];
        if ( ABS(Ax-Cx)>1.0E-8||ABS(Ay-Cy)>1.0E-8||ABS(Az-Cz)>1.0E-8 ) break;
    }

    Dx = Cx - Ax;
    Dy = Cy - Ay;
    Dz = Cz - Az;

    for( ;i>0; i-- )
    {
        Cx = X[i] - Ax;
        Cy = Y[i] - Ay;
        Cz = Z[i] - Az;

        Nx =  Dy*Cz - Cy*Dz;
        Ny = -Dx*Cz - Cx*Dz;
        Nz =  Dx*Cy - Cx*Dy;
        if ( Nx*Nx + Ny*Ny + Nz*Nz > 1.0E-8 ) break;
    }

    D = -Ax*Nx - Ay*Ny - Az*Nz;

    R = 0.0;
    for( i=0; i<16; i++ ) R += ABS(Nx*X[i] + Ny*Y[i] + Nz*Z[i] + D);

    return R<1.0E-4;
}


/*******************************************************************************

Subdive a bicubic polynomial surface in longer of the parameters. Also compute
a bounding volume for the surface.

23 Aug 1995

*******************************************************************************/
void BiCubicSubdivide( Geometry_t *Geometry, int SubLev,int Where )
{
     BiCubic_t *LeftCubic, *RightCubic;
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
     LeftCubic = Geometry->Left->BiCubic = (BiCubic_t *)malloc(sizeof(BiCubic_t));

     Geometry->Right = (Geometry_t *)calloc(sizeof(Geometry_t),1);
     RightCubic = Geometry->Right->BiCubic = (BiCubic_t *)malloc(sizeof(BiCubic_t));

     Geometry->Left->GeometryType = Geometry->Right->GeometryType = GEOMETRY_BICUBIC;

     ULength = BiCubicLength(Geometry,1);
     VLength = BiCubicLength(Geometry,0);

     if ( ULength>VLength )
     {
         for( j=0; j<6; j++ )
         {
             BiCubicBezierSubdivideHalfU( Geometry->BiCubic->BezierFactors[j],
                 LeftCubic->BezierFactors[j],RightCubic->BezierFactors[j] );
         }
     } else
     {
         for( j=0; j<6; j++ )
         {
             BiCubicBezierSubdivideHalfV( Geometry->BiCubic->BezierFactors[j],
                LeftCubic->BezierFactors[j],RightCubic->BezierFactors[j] );
         }
     }

     for( j=0; j<6; j++ )
     {
         BiCubicBezierToMonomial(  LeftCubic->PolyFactors[j],LeftCubic->BezierFactors[j] );
         BiCubicBezierToMonomial( RightCubic->PolyFactors[j],RightCubic->BezierFactors[j] );
     }

#if 0
     if ( Where ) return;
#endif

     BBox = &Geometry->Left->BBox;
     BBox->XMin = BBox->XMax = LeftCubic->BezierFactors[0][0];
     BBox->YMin = BBox->YMax = LeftCubic->BezierFactors[1][0];
     BBox->ZMin = BBox->ZMax = LeftCubic->BezierFactors[2][0];

     for( j=0; j<16; j++ )
     {
         BBox->XMin = MIN(BBox->XMin,LeftCubic->BezierFactors[0][j]);
         BBox->XMax = MAX(BBox->XMax,LeftCubic->BezierFactors[0][j]);

         BBox->YMin = MIN(BBox->YMin,LeftCubic->BezierFactors[1][j]);
         BBox->YMax = MAX(BBox->YMax,LeftCubic->BezierFactors[1][j]);

         BBox->ZMin = MIN(BBox->ZMin,LeftCubic->BezierFactors[2][j]);
         BBox->ZMax = MAX(BBox->ZMax,LeftCubic->BezierFactors[2][j]);
     }

     BBox->XMin = MAX(BBox->XMin,Geometry->BBox.XMin);
     BBox->YMin = MAX(BBox->YMin,Geometry->BBox.YMin);
     BBox->ZMin = MAX(BBox->ZMin,Geometry->BBox.ZMin);
     BBox->XMax = MIN(BBox->XMax,Geometry->BBox.XMax);
     BBox->YMax = MIN(BBox->YMax,Geometry->BBox.YMax);
     BBox->ZMax = MIN(BBox->ZMax,Geometry->BBox.ZMax);

     BBox = &Geometry->Right->BBox;
     BBox->XMin = BBox->XMax = RightCubic->BezierFactors[0][0];
     BBox->YMin = BBox->YMax = RightCubic->BezierFactors[1][0];
     BBox->ZMin = BBox->ZMax = RightCubic->BezierFactors[2][0];

     for( j=0; j<16; j++ )
     {
         BBox->XMin = MIN(BBox->XMin,RightCubic->BezierFactors[0][j]);
         BBox->XMax = MAX(BBox->XMax,RightCubic->BezierFactors[0][j]);
         BBox->YMin = MIN(BBox->YMin,RightCubic->BezierFactors[1][j]);
         BBox->YMax = MAX(BBox->YMax,RightCubic->BezierFactors[1][j]);
         BBox->ZMin = MIN(BBox->ZMin,RightCubic->BezierFactors[2][j]);
         BBox->ZMax = MAX(BBox->ZMax,RightCubic->BezierFactors[2][j]);
     }

     BBox->XMin = MAX(BBox->XMin,Geometry->BBox.XMin);
     BBox->YMin = MAX(BBox->YMin,Geometry->BBox.YMin);
     BBox->ZMin = MAX(BBox->ZMin,Geometry->BBox.ZMin);
     BBox->XMax = MIN(BBox->XMax,Geometry->BBox.XMax);
     BBox->YMax = MIN(BBox->YMax,Geometry->BBox.YMax);
     BBox->ZMax = MIN(BBox->ZMax,Geometry->BBox.ZMax);

     if ( BiCubicIsAPlane( Geometry->Left->BiCubic->BezierFactors[0],
                           Geometry->Left->BiCubic->BezierFactors[1],
                           Geometry->Left->BiCubic->BezierFactors[2] ) )
         Geometry->Left->Flags |= GEOMETRY_FLAG_PLANE;

     if ( BiCubicIsAPlane( Geometry->Right->BiCubic->BezierFactors[0],
                           Geometry->Right->BiCubic->BezierFactors[1],
                           Geometry->Right->BiCubic->BezierFactors[2] ) )
         Geometry->Right->Flags |= GEOMETRY_FLAG_PLANE;
}
 
/*******************************************************************************

Compute element of (iso)line for a (bi)cubic polynomial.

23 Aug 1995

*******************************************************************************/
double BiCubicEofL(double U,double V,double *X,double *Y,double *Z,int UnotV)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = BiCubicPartialU(U,V,X);
    dXdV = BiCubicPartialV(U,V,X);

    dYdU = BiCubicPartialU(U,V,Y);
    dYdV = BiCubicPartialV(U,V,Y);

    dZdU = BiCubicPartialU(U,V,Z);
    dZdV = BiCubicPartialV(U,V,Z);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    if ( UnotV )
        return sqrt(Auu);
    else
        return sqrt(Avv);
}

/*******************************************************************************

Compute element of area for a bicubic polynomial.

23 Aug 1995

*******************************************************************************/
double BiCubicEofA(double U,double V,double *X,double *Y,double *Z)
{
    double dXdU,dYdU,dZdU,dXdV,dYdV,dZdV,Auu,Auv,Avv,detA;

    int i;

    dXdU = BiCubicPartialU(U,V,X);
    dXdV = BiCubicPartialV(U,V,X);

    dYdU = BiCubicPartialU(U,V,Y);
    dYdV = BiCubicPartialV(U,V,Y);

    dZdU = BiCubicPartialU(U,V,Z);
    dZdV = BiCubicPartialV(U,V,Z);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    detA = Auu*Avv - Auv*Auv;

    return sqrt(detA);
}

/*******************************************************************************

Compute isoline length of a bicubic polynomial in given parameter.

23 Aug 1995

*******************************************************************************/
double BiCubicLength( Geometry_t *Geometry,int UnotV )
{
    double *X,*Y,*Z,EofL,Length;

    int i,j=0;

    X = Geometry->BiCubic->PolyFactors[0];
    Y = Geometry->BiCubic->PolyFactors[1];
    Z = Geometry->BiCubic->PolyFactors[2];

    Length = 0.0;
    for( j=0; j<4; j++ )
    for( i=0; i<N_Integ1d; i++ )
    {
        if ( UnotV )
            EofL = BiCubicEofL(U_Integ1d[i],j/3.0,X,Y,Z,1);
        else
            EofL = BiCubicEofL(j/3.0,U_Integ1d[i],X,Y,Z,0);

        Length += S_Integ1d[i]*EofL;
    }

    return Length;
}

/*******************************************************************************

Compute area of a bicubic polynomial surface.

23 Aug 1995

*******************************************************************************/
double BiCubicArea( Geometry_t *Geometry )
{
    double *X,*Y,*Z,EofA,Area;

    int i;

    X = Geometry->BiCubic->PolyFactors[0];
    Y = Geometry->BiCubic->PolyFactors[1];
    Z = Geometry->BiCubic->PolyFactors[2];

    Area = 0.0;
    for( i=0; i<N_Integ; i++ )
    {
        EofA = BiCubicEofA(U_Integ[i],V_Integ[i],X,Y,Z);
        Area += S_Integ[i]*EofA;
    }

    return Area;
}

/*******************************************************************************

Compute differential area to area viewfactor for bicubic surface elements by
direct numerical integration.

24 Aug 1995

*******************************************************************************/
double BiCubicIntegrateDiffToArea( Geometry_t *GB,
  double FX,double FY,double FZ,double NFX,double NFY,double NFZ)
{
    double F,R,cosA,cosB,EA,EAF,EAT,PI=2*acos(0.0);
    double DX,DY,DZ,NTX,NTY,NTZ,U,V;

    double *BX  = GB->BiCubic->PolyFactors[0];
    double *BY  = GB->BiCubic->PolyFactors[1];
    double *BZ  = GB->BiCubic->PolyFactors[2];

    double *NBX = GB->BiCubic->PolyFactors[3];
    double *NBY = GB->BiCubic->PolyFactors[4];
    double *NBZ = GB->BiCubic->PolyFactors[5];

    int i,j;


    F  = 0.0;
    for( i=0; i<N_Integ; i++ )
    {
        U = U_Integ[i];
        V = V_Integ[i];
        
        DX  = BiCubicValue(U,V,BX) - FX;
        DY  = BiCubicValue(U,V,BY) - FY;
        DZ  = BiCubicValue(U,V,BZ) - FZ;

        cosA =  DX*NFX + DY*NFY + DZ*NFZ;
        if ( cosA <= 1.0E-9 ) continue;

        NTX = BiCubicValue(U,V,NBX);
        NTY = BiCubicValue(U,V,NBY);
        NTZ = BiCubicValue(U,V,NBZ);

        R = NTX*NTX + NTY*NTY + NTZ*NTZ;
        if ( ABS(1-R) > 1.0E-7 )
        {
            R = 1.0/sqrt(R);
            NTX *= R;
            NTY *= R;
            NTZ *= R;
        }

        cosB = -DX*NTX - DY*NTY - DZ*NTZ;
        if ( cosB <= 1.0E-9 ) continue;

        R = DX*DX + DY*DY + DZ*DZ;

        EA = BiCubicEofA(U,V,BX,BY,BZ);
        F += EA*cosA*cosB*S_Integ[i]/(R*R);
    }

    return F;
}

/*******************************************************************************

Compute area to area viewfactor for bicubic surface elements by subdivision
and direct numerical numerical integration when the differential viewfactors
match given magnitude criterion or areas of the elements are small enough.
Blocking of the view between the elements is resolved by ray traceing.

24 Aug 1995

*******************************************************************************/
void BiCubicComputeViewFactors(Geometry_t *GA,Geometry_t *GB,int LevelA,int LevelB )
{
    double FX,FY,FZ,DX,DY,DZ,U,V,Hit;
    double F,Fa,Fb,EA,PI=2*acos(0.0);

    double *AX = GA->BiCubic->PolyFactors[0];
    double *AY = GA->BiCubic->PolyFactors[1];
    double *AZ = GA->BiCubic->PolyFactors[2];

    double *BX = GB->BiCubic->PolyFactors[0];
    double *BY = GB->BiCubic->PolyFactors[1];
    double *BZ = GB->BiCubic->PolyFactors[2];

    int i,j;

    GA -> Flags |= GEOMETRY_FLAG_USED_INTEG;
    GB -> Flags |= GEOMETRY_FLAG_USED_INTEG;

#ifdef TODO
    Fa = Fb = BiCubicIntegrateDiffToArea( GA,0.5,0.5,GB );
    if ( GA != GB )
    {
        Fb = BiCubicIntegrateDiffToArea( GB,0.5,0.5,GA );
    }
#endif
    if ( Fa < 1.0E-10 && Fb < 1.0E-10 ) return;

    if ( (Fa<FactorEPS || GB->Area<AreaEPS) && (Fb<FactorEPS || GA->Area<AreaEPS) )
    {
        GeometryList_t *Link;

        Hit = N_Integ;
        for( i=0; i<N_Integ; i++ )
        {
            U = drand48();
            V = drand48();

            FX = BiCubicValue( U,V,AX );
            FY = BiCubicValue( U,V,AY );
            FZ = BiCubicValue( U,V,AZ );

            U = drand48();
            V = drand48();

            DX = BiCubicValue( U,V,BX ) - FX;
            DY = BiCubicValue( U,V,BY ) - FY;
            DZ = BiCubicValue( U,V,BZ ) - FZ;

            Hit -= RayHitGeometry( FX,FY,FZ,DX,DY,DZ );
        }

        if ( Hit == 0 ) return;

        if ( Hit == N_Integ || (Fa<FactorEPS/2 && Fb<FactorEPS/2) )
        {
            Hit /= (double)N_Integ;

            F = 0.0;
            for( i=0; i<N_Integ; i++ )
            {
                U = U_Integ[i];
                V = V_Integ[i];
                EA = BiCubicEofA( U,V,AX,AY,AZ );
#ifdef TODO
                F += S_Integ[i]*EA*BiCubicIntegrateDiffToArea( GA,U,V,GB );
#endif
            }

            F = Hit*F / PI;
            Fb = F / GB->Area;
            Fa = F / GA->Area;

            Link = (GeometryList_t *)calloc( sizeof(GeometryList_t),1 );
            Link->Next = GA->Link;
            GA->Link = Link;

            Link->Entry = GB;
            Link->ViewFactor = Fa;

            if ( GA != GB )
            {
                Link = (GeometryList_t *)calloc( sizeof(GeometryList_t),1 );
                Link->Next = GB->Link;
                GB->Link = Link;

                Link->Entry = GA;
                Link->ViewFactor = Fb;
            }

            return;
        }
    }

#if 0
    if ( GA == GB )
    {
        if ( !GB->Left ) BiCubicSubdivide( GB, Level,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;
            GB->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area )
            {
                GB->Left->Area  = BiCubicArea( GB->Left );
                GB->Right->Area = BiCubicArea( GB->Right );
            }
        }

        BiCubicComputeViewFactors( GB->Left, GB->Left, Level+1 );
        BiCubicComputeViewFactors( GB->Right,GB->Right,Level+1 );
        BiCubicComputeViewFactors( GB->Left, GB->Right,Level+1 );
    }
    else if ( Fa > Fb )
    {
        if ( !GB->Left ) BiCubicSubdivide( GB, Level,1 );

        if ( GB->Flags & GEOMETRY_FLAG_LEAF )
        {    
            GB->Flags &= ~GEOMETRY_FLAG_LEAF;
            GB->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GB->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GB->Left->Area )
            {
                GB->Left->Area  = BiCubicArea( GB->Left );
                GB->Right->Area = BiCubicArea( GB->Right );
            }
        }

        BiCubicComputeViewFactors( GA,GB->Left,Level+1 );
        BiCubicComputeViewFactors( GA,GB->Right,Level+1 );
    } else
    {
        if ( !GA->Left ) BiCubicSubdivide( GA, Level,1 );

        if ( GA->Flags & GEOMETRY_FLAG_LEAF )
        {
            GA->Flags &= ~GEOMETRY_FLAG_LEAF;
            GA->Left->Flags  |= GEOMETRY_FLAG_LEAF;
            GA->Right->Flags |= GEOMETRY_FLAG_LEAF;

            if ( !GA->Left->Area )
            {
                GA->Left->Area  = BiCubicArea(GA->Left);
                GA->Right->Area = BiCubicArea(GA->Right);
            }
        }

        BiCubicComputeViewFactors( GA->Left,GB,Level+1 );
        BiCubicComputeViewFactors( GA->Right,GB,Level+1 );
    }
#endif
}
