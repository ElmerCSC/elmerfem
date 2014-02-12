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

Ray intersection tests for various geometric objects.

Juha Ruokolainen/CSC - 23 Aug 1995

*******************************************************************************/

#include <ViewFactors.h>

static double REPS = 1.0E-4;
#define MAX_LEVEL 16

/*******************************************************************************

Solve for volume box/ray (btw. points (Fx,Fy,Fz)->(Tx,Ty,Tz ), end points not
included) intersection.  Return value is whether there is a hit or not. If 
either of points is inside the box hit is given.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitBBox( BBox_t *BBox,double Fx,double Fy,double Fz,double Tx,double Ty,double Tz)
{
     double s,XMin,XMax,TMin,TMax;

     Tx += Fx;
     Ty += Fy;
     Tz += Fz;

     if (
         BBox->XMin>Fx && BBox->XMin>Tx || BBox->XMax<Fx && BBox->XMax<Tx ||
         BBox->YMin>Fy && BBox->YMin>Ty || BBox->YMax<Fy && BBox->YMax<Ty ||
         BBox->ZMin>Fz && BBox->ZMin>Tz || BBox->ZMax<Fz && BBox->ZMax<Tz 
      ) return FALSE;

     if (
           Fx>=BBox->XMin && Fx<=BBox->XMax && Fy>=BBox->YMin && 
           Fy<=BBox->YMax && Fz>=BBox->ZMin && Fz<=BBox->ZMax
       ) return TRUE;

     if (
           Tx>=BBox->XMin && Tx<=BBox->XMax && Ty>=BBox->YMin && 
           Ty<=BBox->YMax && Tz>=BBox->ZMin && Tz<=BBox->ZMax
       ) return TRUE;

     Tx -= Fx;
     Ty -= Fy;
     Tz -= Fz;

     if ( ABS(Tx)<1.0e-12 && (Fx<BBox->XMin || Fx>BBox->XMax) ) return FALSE;
     if ( ABS(Ty)<1.0e-12 && (Fy<BBox->YMin || Fy>BBox->YMax) ) return FALSE;
     if ( ABS(Tz)<1.0e-12 && (Fz<BBox->ZMin || Fz>BBox->ZMax) ) return FALSE;

     XMin = -DBL_MAX;
     XMax =  DBL_MAX;

     if ( ABS(Tx)>=1.0e-12 )
     {
         Tx = 1.0 / Tx;
         TMax = (BBox->XMax-Fx)*Tx;
         TMin = (BBox->XMin-Fx)*Tx;

         if ( TMin>TMax ) { s = TMin; TMin = TMax; TMax = s; }

         XMax = TMax;
         XMin = TMin;
     }

     if ( ABS(Ty)>=1.0e-12 )
     {
         Ty = 1.0 / Ty;
         TMax = (BBox->YMax - Fy)*Ty;
         TMin = (BBox->YMin - Fy)*Ty;

         if ( TMin>TMax ) { s = TMin; TMin = TMax; TMax = s; }

         if ( TMin>XMin ) XMin = TMin;
         if ( TMax<XMax ) XMax = TMax;
 
         if ( XMin>XMax ) return FALSE;
     }

     if ( ABS(Tz)>=1.0e-12 )
     {
         Tz = 1.0/Tz;
         TMin = (BBox->ZMin-Fz)*Tz;
         TMax = (BBox->ZMax-Fz)*Tz;

         if ( TMin>TMax ) { s = TMin; TMin = TMax; TMax = s; }

         if ( TMin>XMin ) XMin = TMin;
         if ( TMax<XMax ) XMax = TMax;
     }

     return XMin<=XMax && XMin>=0.0 && XMin<=1.0;
}

/*******************************************************************************

Find one intersection of a bicubic element and ray segment (FX,FY,FZ)->(FX+DX,
FY+DY,FZ+DZ) by Newton iteration, ray end points not included. The ray is
represented as intersection of two planes here to get rid of one of the
equations...

LAST Modified: 23 Aug 1995

*******************************************************************************/
int SolveRayBiCubicNewton(
                              double RayPlanes[2][4],
           double FX,double FY,double FZ,double DX,double DY,double DZ,double L,
                            double *X,double *Y,double *Z
                         )
{
    double Xdu,Xdv,Ydu,Ydv,Zdu,Zdv,VX,VY,VZ,J0,J1,J2,J3,EPS=REPS;
    double A0,B0,C0,D0,A1,B1,C1,D1,det,F0,F1,NORM,T,U=0.5,V=0.5;

    int iter,ITMAX=30;

    A0 = RayPlanes[0][0];
    B0 = RayPlanes[0][1];
    C0 = RayPlanes[0][2];
    D0 = RayPlanes[0][3];

    A1 = RayPlanes[1][0];
    B1 = RayPlanes[1][1];
    C1 = RayPlanes[1][2];
    D1 = RayPlanes[1][3];

    for( iter=0; iter<ITMAX; iter++ )
    {
        Xdu = BiCubicPartialU(U,V,X);
        Xdv = BiCubicPartialV(U,V,X);

        Ydu = BiCubicPartialU(U,V,Y);
        Ydv = BiCubicPartialV(U,V,Y);

        Zdu = BiCubicPartialU(U,V,Z);
        Zdv = BiCubicPartialV(U,V,Z);

        J3 =  (A0*Xdu + B0*Ydu + C0*Zdu);
        J1 = -(A0*Xdv + B0*Ydv + C0*Zdv);
        J2 = -(A1*Xdu + B1*Ydu + C1*Zdu);
        J0 =  (A1*Xdv + B1*Ydv + C1*Zdv);

        det = J0*J3 - J1*J2;
        if ( ABS(det)<1.0E-12 ) return FALSE;
        det = 1.0/det; 

        VX = BiCubicValue(U,V,X);
        VY = BiCubicValue(U,V,Y);
        VZ = BiCubicValue(U,V,Z);

        F0 = A0*VX + B0*VY + C0*VZ + D0;
        F1 = A1*VX + B1*VY + C1*VZ + D1;

        NORM = F0*F0 + F1*F1;
        if ( NORM<1.0E-20 ) break;

        U -= det*(J0*F0 + J1*F1);
        V -= det*(J2*F0 + J3*F1);

        if ( U<-0.5 || U>1.5 || V<-0.5 || V>1.5 ) return FALSE;
    }

    if ( iter >= ITMAX )
    {
       fprintf( stderr, "NEWTON ITMAX EXCEEDED: %g %g %g\n", U,V,NORM );
       return -1;
    }

    if ( U<0.0 || U>1.0 || V<0.0 || V>1.0 ) return FALSE;

    if ( ABS(DX)>ABS(DY) && ABS(DX)>ABS(DZ) ) 
    {
        T = (VX-FX)/DX;
    } else if ( ABS(DY)>ABS(DZ) )
    {
        T = (VY-FY)/DY;
    } else
    {
        T = (VZ-FZ)/DZ;
    }

    return T>EPS && T<1-EPS;
}


/*******************************************************************************

Find if a bicubic element and the ray segment (FX,FY,FZ)->(FX+DX,FY+DY,FZ+DZ)
intersect within the element boundaries, ray end points not included, by
subdivision. When the element is considered to be planar enough, apply Newton
iteration. The subdivision side effect is to generate a bounding box
hierarchy for the element, which then can be used to guide the root finding
to different parts of the subdivided element for subsequent rays.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int SolveRayBiCubic(
                   Geometry_t *Geometry,double RayPlanes[2][4],
  double FX,double FY,double FZ,double DX,double DY,double DZ,double L,int SubLev
                   )
{
    int T;

    if ( Geometry->Flags & GEOMETRY_FLAG_PLANE ) 
    {
        T = SolveRayBiCubicNewton(   RayPlanes,FX,FY,FZ,DX,DY,DZ,L,
                                    Geometry->BiCubic->PolyFactors[0],
                                    Geometry->BiCubic->PolyFactors[1],
                                    Geometry->BiCubic->PolyFactors[2] );

         if ( T<0 )
         {
             Geometry->Flags &= ~GEOMETRY_FLAG_PLANE;
         }
         else return T;
    }

    if ( !Geometry->Left )
    {
        if (
              BiCubicIsAPlane( Geometry->BiCubic->BezierFactors[0],
                               Geometry->BiCubic->BezierFactors[1],
                               Geometry->BiCubic->BezierFactors[2] )
           )
        {
            Geometry->Flags |= GEOMETRY_FLAG_PLANE;

            T = SolveRayBiCubicNewton(  RayPlanes,FX,FY,FZ,DX,DY,DZ,L,
                                       Geometry->BiCubic->PolyFactors[0],
                                       Geometry->BiCubic->PolyFactors[1],
                                       Geometry->BiCubic->PolyFactors[2] );

            if ( T<0 )
            {
                 Geometry->Flags &= ~GEOMETRY_FLAG_PLANE;
            }
            else return T;
        }

        BiCubicSubdivide( Geometry,SubLev,0 );
    }

    if ( RayHitBBox( &Geometry->Left->BBox,FX,FY,FZ,DX,DY,DZ ) )
    {
        if ( SolveRayBiCubic( Geometry->Left,RayPlanes,FX,FY,FZ,DX,DY,DZ,L,SubLev+1 ) ) return TRUE;
    }

    if ( RayHitBBox( &Geometry->Right->BBox,FX,FY,FZ,DX,DY,DZ ) )
    {
        if ( SolveRayBiCubic(Geometry->Right,RayPlanes,FX,FY,FZ,DX,DY,DZ,L,SubLev+1) ) return TRUE;
    }

    return FALSE;
}

/******************************************************************************

Find if a bicubic element and the ray segment (FX,FY,FZ)->(FX+DX,FY+DY,FZ+DZ)
intersect within the element boundaries, ray end points not included, by
subdivision. When the element is considered to be planar enough, apply Newton
iteration. The subdivision side effect is to generate a bounding volume
hierarchy of the element, which then can be used to guide the root finding to
different parts of the subdivided element.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitBiCubic
  (
   Geometry_t *Geometry,double FX,double FY,double FZ,double DX,double DY,double DZ,double L
  )
{
    double A0,B0,C0,D0,A1,B1,C1,D1,RayPlanes[2][4];

    if ( !RayHitBBox(&Geometry->BBox,FX,FY,FZ,DX,DY,DZ) ) return FALSE;

    if ( ABS(DX)>ABS(DY) && ABS(DX)>ABS(DZ) ) 
    {
        B0 = 1;
        C0 = 2;
        A0 = (-C0*DZ - B0*DY)/DX;
        B1 = 2;
        C1 = 3;
        A1 = (-B1*DY - C1*DZ)/DX;
    } else if ( ABS(DY)>ABS(DX) && ABS(DY)>ABS(DZ) )
    {
        A0 = 1;
        C0 = 2;
        B0 = (-C0*DZ - A0*DX)/DY;
        A1 = 3;
        C1 = 1;
        B1 = (-A1*DX - C1*DZ)/DY;
    } else 
    {
        A0 = 1;
        B0 = 2;
        C0 = (-B0*DY - A0*DX)/DZ;
        A1 = 3;
        B1 = 1;
        C1 = (-A1*DX - B1*DY)/DZ;
    }

    D0 = -A0*FX - B0*FY - C0*FZ;
    D1 = -A1*FX - B1*FY - C1*FZ;

    RayPlanes[0][0] = A0;
    RayPlanes[0][1] = B0;
    RayPlanes[0][2] = C0;
    RayPlanes[0][3] = D0;

    RayPlanes[1][0] = A1;
    RayPlanes[1][1] = B1;
    RayPlanes[1][2] = C1;
    RayPlanes[1][3] = D1;

    return SolveRayBiCubic( Geometry,RayPlanes,FX,FY,FZ,DX,DY,DZ,L,0 );
}

/*******************************************************************************

Solve for bilinear element and ray segment intersection. Return value is if
there is a hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitBiLinear(
      Geometry_t *Geometry,double FX,double FY,double FZ,
          double DX,double DY,double DZ,double L
       )
{
    double A0,B0,C0,D0,A1,B1,C1,D1,A,B,C,D,E,AC4,U,V,T,*X,*Y,*Z;
    double EPS=REPS,TF,TD,*TX;

    X = Geometry->BiLinear->PolyFactors[0];
    Y = Geometry->BiLinear->PolyFactors[1];
    Z = Geometry->BiLinear->PolyFactors[2];

    A0 = -DY*X[1] + DX*Y[1];
    B0 = -DY*X[2] + DX*Y[2];
    C0 = -DY*X[3] + DX*Y[3];
    D0 = -DY*(FX-X[0]) + DX*(FY-Y[0]);

    A1 = -DY*Z[1] + DZ*Y[1];
    B1 = -DY*Z[2] + DZ*Y[2];
    C1 = -DY*Z[3] + DZ*Y[3];
    D1 = -DY*(FZ-Z[0]) + DZ*(FY-Y[0]);

    if ( ABS(DX)>ABS(DY) && ABS(DX)>ABS(DZ) )
    {
       TD = DX; TF = FX; TX=X;
    } else if ( ABS(DY)>ABS(DZ) )
    {
       TD = DY; TF = FY; TX=Y;
    } else
    {
       TD = DZ; TF = FZ; TX=Z;
    }

    if ( ABS(C0)<1.0e-12 && ABS(C1)<1.0e-12 )
    {
       A = A0*B1 - A1*B0;

       if ( ABS(A)<1.0e-12 ) return FALSE;
       A = 1.0 / A;

       U =  A * (B1*D0 - B0*D1);
       if ( U>=0.0 && U<=1.0 )
       {
          V = A * (A0*D1 - A1*D0);
          if ( V>=0.0 && V<=1.0 )
          {
             T = (BiLinearValue(U,V,TX)-TF) / TD;
             return T>EPS && T<1-EPS;
          }
       }
       return FALSE;
    }

    if ( ABS(C0)<1.0e-12 && ABS(A0)<1.0e-12 )
    {
       V = D0 / B0;
       if ( V>=0.0 && V<=1.0 )
       {
          U = (D1 - B1*V) / (A1 + C1*V);
          if ( U>=0.0 && U<=1.0 )
          {
             T = (BiLinearValue(U,V,TX)-TF) / TD;
             return T>EPS && T<1-EPS;
          }
       }
       return FALSE;
    }

    if ( ABS(C1) < 1.0e-12 && ABS(A1) < 1.0e-12 )
    {
       V = D1 / B1;
       if ( V>=0.0 && V<=1.0 )
       {
          U = (D0 - B0*V) / (A0 + C0*V);
          if ( U>=0.0 && U<=1.0 )
          {
             T = (BiLinearValue(U,V,TX)-TF) / TD;
             return T>EPS && T<1-EPS;
          }
       }
       return FALSE;
    }

    if ( ABS(C0)<1.0e-12 && ABS(B0)<1.0e-12 )
    {
       U = D0 / A0;
       if ( U>=0.0 && U<=1.0 )
       {
          V = (D1 - A1*U) / (B1 + C1*U);
          if ( V>=0.0 && V<=1.0 )
          {
             T = (BiLinearValue(U,V,TX)-TF) / TD;
             return T>EPS && T<1-EPS;
           }
       }
       return FALSE;
    }

    if ( ABS(C1)<1.0e-12 && ABS(B1)<1.0e-12 )
    {
       U = D1 / A1;
       if ( U>=0.0 && U<=1.0 )
       {
          V = (D0 - A0*U) / (B0 + C0*U);
          if ( V>=0.0 && V<=1.0 )
          {
             T = (BiLinearValue(U,V,TX)-TF) / TD;
             return T>EPS && T<1-EPS;
          }
       }
       return FALSE;
    }

    A = C0*B1 - B0*C1;
    B = D0*C1 - B0*A1 + A0*B1 - C0*D1;
    C = D0*A1 - A0*D1;

    if ( ABS(A)<1.0E-12 )
    {
        if ( ABS(B) < 1.0e-12 ) return FALSE;

        V = -C/B;
        if ( V<0.0 || V>1.0 ) return FALSE;

        if ( ABS(A0+C0*V) < 1.0e-12 )
          U = (D1 - B1*V) / (A1 + C1*V);
        else
          U = (D0 - B0*V) / (A0 + C0*V);

        if ( U<0.0 || U>1.0 ) return FALSE;
  
        T = (BiLinearValue(U,V,TX)-TF) / TD;
        return T>EPS && T<1-EPS;
    }

    if ( (D = B*B - 4*A*C)<0.0 ) return FALSE;

    D = sqrt(D);

    if ( B > 0 ) {
      V = -2*C/(B+D);
    } else {
      V = (-B+D)/(2*A);
    }
    if ( V>=0.0 && V<=1.0 )
    {
        if ( ABS(A0+C0*V) < 1.0e-12 )
          U = (D1 - B1*V) / (A1 + C1*V);
        else
          U = (D0 - B0*V) / (A0 + C0*V);

        if ( U>=0.0 && U<=1.0 )
        {
            T = (BiLinearValue(U,V,TX)-TF) / TD;
            if ( T>EPS && T<1-EPS ) return TRUE;
        }
    }

    if ( B > 0 ) {
      V = -(B+D)/(2*A);
    } else {
      V = 2*C/(-B+D);
    }
    if ( V>=0.0 && V<=1.0 )
    {
        if ( ABS(A0+C0*V) < 1.0e-12 )
          U = (D1 - B1*V) / (A1 + C1*V);
        else
          U = (D0 - B0*V) / (A0 + C0*V);

        if ( U>=0.0 && U<=1.0 )
        {
           T = (BiLinearValue(U,V,TX)-TF) / TD;
           if ( T>EPS && T<1-EPS ) return TRUE;
        }
    }

    return FALSE;
}
/*******************************************************************************

Solve for linear line element and ray segment intersection. Return value is if
there is a hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitLine(
      Geometry_t *Geometry,double FX,double FY,double FZ,
          double DX,double DY,double DZ,double L
       )
{
    double A11,A12,A21,A22,U,V,T,*X,*Y,*Z,detA;
    double EPS=REPS,TF,TD,*TX;

    X = Geometry->Linear->PolyFactors[0];
    Y = Geometry->Linear->PolyFactors[1];

    A11 = DX;
    A21 = DY;
    A12 = -X[1];
    A22 = -Y[1];
    detA = A11*A22 - A12*A21;

    if ( ABS(detA) < 1e-9 ) return FALSE;

    U = ( A22 * (X[0]-FX) - A12*(Y[0]-FY)) / detA;
    if ( U<EPS || U>1-EPS ) return FALSE;

    U = (-A21 * (X[0]-FX) + A11*(Y[0]-FY)) / detA;
    if ( U<EPS || U>1-EPS ) return FALSE;

    return TRUE;
}

/*******************************************************************************

Solve for triangle element and ray segment intersection. Return value is if
there is a hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitTriangle(
      Geometry_t *Geometry,double FX,double FY,double FZ,double DX,double DY,double DZ,double L
       )
{
    double A11,A12,A13,A21,A22,A23,A31,A32,A33,U,V,T,*X,*Y,*Z;
    double B11,B12,B13,B21,B22,B23,B31,B32,B33,PX,PY,PZ,detA;
    double EPS=REPS;

/*
 * if ( !RayHitBBox(&Geometry->BBox,FX,FY,FZ,DX,DY,DZ) ) return FALSE;
 */

    X = Geometry->Triangle->PolyFactors[0];
    Y = Geometry->Triangle->PolyFactors[1];
    Z = Geometry->Triangle->PolyFactors[2];

    A11 = X[1];
    A12 = X[2];
    A13 = -DX;

    A21 = Y[1];
    A22 = Y[2];
    A23 = -DY;

    A31 = Z[1];
    A32 = Z[2];
    A33 = -DZ;

    detA = A11 * ( A22*A33 - A32*A23 ) +
           A12 * ( A31*A23 - A21*A33 ) +
           A13 * ( A21*A32 - A31*A22 );

    if ( ABS(detA)<1.0e-9 ) return FALSE;

    detA = 1.0 / detA;

    PX = FX - X[0];
    PY = FY - Y[0];
    PZ = FZ - Z[0];

    B31 = A21*A32 - A22*A31;
    B32 = A12*A31 - A11*A32;
    B33 = A11*A22 - A21*A12;

    T = detA * ( B31*PX + B32*PY + B33*PZ );
    if ( T<EPS || T>1.0-EPS ) return FALSE;

    B11 = A22*A33 - A23*A32;
    B12 = A32*A13 - A12*A33;
    B13 = A12*A23 - A22*A13;

    U = detA * ( B11*PX + B12*PY + B13*PZ );
    if ( U<0.0 || U>1.0 ) return FALSE;

    B21 = A23*A31 - A21*A33;
    B22 = A11*A33 - A31*A13;
    B23 = A21*A13 - A11*A23;

    V = detA * ( B21*PX + B22*PY + B23*PZ );
    if ( V<0.0 || V>1.0 ) return FALSE;

    return (U+V) <= 1.0;
}

/*******************************************************************************

Solve for sphere and ray segment intersection. Return value is if there is a
hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitSphere(
       Geometry_t *Geometry,double FX,double FY,double FZ,
           double DX,double DY,double DZ,double L
          )
{
     Sphere_t *Sph = Geometry->Sphere;

     double A,B,C,D,E,AC4,T,R=Sph->Radius,EPS=REPS;

     FX -= Sph->CenterPoint.x;
     FY -= Sph->CenterPoint.y;
     FZ -= Sph->CenterPoint.z;

     if (
            FX>R && FX+DX>R || FX<-R && FX+DX<-R ||
            FY>R && FY+DY>R || FY<-R && FY+DY<-R ||
            FZ>R && FZ+DZ>R || FZ<-R && FZ+DZ<-R
         ) return FALSE;

     A = DX*DX + DY*DY + DZ*DZ;
     if ( A<1.0E-12 ) return FALSE;

     B = 2*(FX*DX + FY*DY + FZ*DZ);
     C = FX*FX + FY*FY + FZ*FZ - R*R;

     if ( B>=0 && C>=0 ) return FALSE;

     if ( (D = B*B - 4*A*C)<0.0 ) return FALSE;

     D = sqrt(D);

     if ( B > 0 ) {
       T = -2*C/(B+D);
     } else {
       T = (-B+D)/(2*A);
     }
     if ( T>EPS && T<1-EPS ) return TRUE;

     if  ( B > 0 ) {
       T = -(B+D)/(2*A);
     } else {
       T = 2*C/(-B+D);
     }
     return T>EPS && T<1-EPS;
}

/*******************************************************************************

Solve for cylinder, in canonical coordinates, and ray segment intersection.
Return value is if there is a hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitCylinder(
      Geometry_t *Geometry,double FX,double FY,double FZ,
             double DX,double DY,double DZ,double L
        )
{
     Cylinder_t *Cyl = Geometry->Cylinder;

     double A,B,C,D,E,T,X,Y,Z,R=Cyl->Radius,CL=Cyl->Length/2,EPS=REPS;

     FX -= Cyl->CenterPoint.x;
     FY -= Cyl->CenterPoint.y;
     FZ -= Cyl->CenterPoint.z;

     if ( !Cyl->IdentMatrix )
     {
         RotateVector(&FX,&FY,&FZ,Cyl->RotationMatrix);
         RotateVector(&DX,&DY,&DZ,Cyl->RotationMatrix);
     }

     if (
            FX>R && FX+DX>R || FX<-R && FX+DX<-R ||
            FY>R && FY+DY>R || FY<-R && FY+DY<-R ||
            FZ<-CL && FZ+DZ<-CL || FZ>CL && FZ+DZ>CL
        ) return FALSE;

     A = DX*DX + DY*DY;
     B = 2*(FX*DX + FY*DY);
     C = FX*FX + FY*FY - R*R;

     if ( B>=0 && C>=0 ) return FALSE;

     if ( A<1.0E-12 )
     {
         T = -C/B;
         if ( T>EPS && T<1-EPS )
         {
             Z = FZ + T*DZ;
             return Z>=-CL && Z<=CL;
         } else return FALSE;
     }

     if ( (D = (B*B - 4*A*C))<0.0 ) return FALSE;

     D = sqrt(D);

     if ( B > 0 ) {
       T = -2*C/(B+D);
     } else {
       T = (-B+D)/(2*A);
     }
     if ( T>EPS && T<1-EPS )
     {
         Z = FZ + T*DZ;
         if ( Z>=-CL && Z<=CL ) return TRUE;
     }

     if ( B > 0 ) {
       T = (-B+D)/(2*A);
     } else {
       T = 2*C/(-B+D);
     }
     if ( T>EPS && T<1-EPS )
     {
         Z = FZ + T*DZ;
         if ( Z>=-CL && Z<=CL ) return TRUE;
     }

     return FALSE;
}

/*******************************************************************************

Solve for planar circle and ray segment intersection, return value is if there
is hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitCircle(
      Geometry_t *Geometry,double FX,double FY,double FZ,double DX,double DY,double DZ,double L
          )
{
     Circle_t *Cir = Geometry->Circle;

     double T,X,Y,Z,R,R0=Cir->RMin,R1=Cir->RMax,EPS = REPS;

     FX -= Cir->CenterPoint.x;
     FY -= Cir->CenterPoint.y;
     FZ -= Cir->CenterPoint.z;

     if ( !Cir->IdentMatrix )
     {
         RotateVector(&FX,&FY,&FZ,Cir->RotationMatrix);
         RotateVector(&DX,&DY,&DZ,Cir->RotationMatrix);
     }

     if ( ABS(DZ)<1.0E-12 ) return FALSE;

     T = -FZ/DZ;
     if ( T<EPS || T>1-EPS ) return FALSE;

     X = FX + T*DX;
     Y = FY + T*DY;
     R = X*X + Y*Y;

     return (R>=R0*R0) && (R<=R1*R1);
}


/*******************************************************************************

Solve for rotational quadratic and ray segment intersection, return value is
if there is hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitRotationalQuadric(
      Geometry_t *Geometry,double FX,double FY,double FZ,
            double DX,double DY,double DZ,double L
          )
{
     RotQuadric_t *RQd = Geometry->RotQuadric;

     double A,B,C,D,E,T,X,Y,Z,R,R0,R1,EPS = REPS;

     FX -= RQd->CenterPoint.x;
     FY -= RQd->CenterPoint.y;
     FZ -= RQd->CenterPoint.z;

     if ( !RQd->IdentMatrix )
     {
         RotateVector(&FX,&FY,&FZ,RQd->RotationMatrix);
         RotateVector(&DX,&DY,&DZ,RQd->RotationMatrix);
     }

     if (
            FX>RQd->RMax && FX+DX>RQd->RMax || FX<-RQd->RMax && FX+DX<-RQd->RMax ||
            FY>RQd->RMax && FY+DY>RQd->RMax || FY<-RQd->RMax && FY+DY<-RQd->RMax ||
            FZ<RQd->ZMin && FZ+DZ<RQd->ZMin || FZ> RQd->ZMax && FZ+DZ> RQd->ZMax
        ) return FALSE;

     if ( ABS(RQd->ZMax - RQd->ZMin)>1.0E-10 )
     {
         R0  = RQd->RMax - RQd->RMin;
         R0 /= RQd->ZMax - RQd->ZMin;
         R1  = RQd->RMin - R0*RQd->ZMin;
     } else
     {
         if ( ABS(DZ)<1.0E-12 ) return FALSE;

         T = -FZ/DZ;
         if ( T<EPS || T>1-EPS ) return FALSE;

         X = FX + T*DX;
         Y = FY + T*DY;
         R = X*X + Y*Y;

         return (R>=RQd->RMin*RQd->RMin) && (R<=RQd->RMax*RQd->RMax);
     }

     A = DX*DX + DY*DY - R0*R0*DZ*DZ;
     B = 2*(FX*DX + FY*DY - (R0*FZ + R1)*R0*DZ);
     C = FX*FX + FY*FY - (R0*FZ + 2*R1)*R0*FZ - R1*R1;

     if ( A<1.0E-12 )
     {
         T = -C/B;
         if ( T>EPS && T<1-EPS )
         {
             Z = FZ + T*DZ;
             R = R0*Z + R1;

             return Z>=RQd->ZMin && Z<=RQd->ZMax && R>=RQd->RMin && R<=RQd->RMax;
         } else return FALSE;
     }

     if ( (D = (B*B - 4*A*C))<0.0 ) return FALSE;

     D = sqrt(D);

     if ( B > 0 ) {
       T = -2*C/(B+D);
     } else {
       T = (-B+D)/(2*A);
     }
     if ( T>EPS && T<1-EPS )
     {
         Z = FZ + T*DZ;
         R = R0*Z + R1;

         if (Z>=RQd->ZMin && Z<=RQd->ZMax && R>=RQd->RMin && R<=RQd->RMax) return TRUE;
     }

     if ( B > 0 ) {
       T = -(B+D)/(2*A);
     } else {
       T = 2*C/(-B+D);
     }
     if ( T>EPS && T<1-EPS )
     {
         Z = FZ + T*DZ;
         R = R0*Z + R1;

         return Z>=RQd->ZMin && Z<=RQd->ZMax && R>=RQd->RMin && R<=RQd->RMax;
     }

     return FALSE;
}

/*******************************************************************************

Solve for rotational quadratic and ray segment intersection, return value is
if there is hit or not.

LAST Modified: 23 Aug 1995

*******************************************************************************/
int RayHitRotationalQuadric1(
    Geometry_t *Geometry,double FX,double FY,double FZ,
          double DX,double DY,double DZ,double L
        )
{
     RotQuadric_t *RQd = Geometry->RotQuadric;

     double Az=RQd->REqn[2],Bz=RQd->REqn[1],Cz=RQd->REqn[0],EPS = REPS;
     double A,B,C,D,E,T,X,Y,Z,R,R0=RQd->RMin*RQd->RMin,R1=RQd->RMax*RQd->RMax;

     FX -= RQd->CenterPoint.x;
     FY -= RQd->CenterPoint.y;
     FZ -= RQd->CenterPoint.z;

     if ( !RQd->IdentMatrix )
     {
         RotateVector(&FX,&FY,&FZ,RQd->RotationMatrix);
         RotateVector(&DX,&DY,&DZ,RQd->RotationMatrix);
     }

     if (
            FX>RQd->RMax && FX+DX>RQd->RMax || FX<-RQd->RMax && FX+DX<-RQd->RMax ||
            FY>RQd->RMax && FY+DY>RQd->RMax || FY<-RQd->RMax && FY+DY<-RQd->RMax ||
            FZ<RQd->ZMin && FZ+DZ<RQd->ZMin || FZ> RQd->ZMax && FZ+DZ> RQd->ZMax
        ) return FALSE;

     if ( ABS(Az)<1.0E-12 && ABS(Bz)<1.0E-12 && ABS(Cz)<1.0E-12 )
     {
         /*  This might be a planar circle */
         if ( ABS(DZ)<1.0E-12 ) return FALSE;

         T = -FZ/DZ;
         if ( T<EPS || T>1-EPS ) return FALSE;

         X = FX + T*DX;
         Y = FY + T*DY;
         R = X*X + Y*Y;

         return R>=R0 && R<=R1;
     }

     A = DX*DX + DY*DY - Az*DZ*DZ;
     B = 2*(FX*DX + FY*DY) - (2*Az*FZ + Bz)*DZ;
     C = FX*FX + FY*FY - (Az*FZ + 2*Bz)*FZ - Cz;

     if ( A<1.0E-12 )
     {
         T = -C/B;
         if ( T>EPS && T<1-EPS )
         {
             Z = FZ + T*DZ;
             R = (Az*Z + Bz)*Z + Cz;

             return Z>=RQd->ZMin && Z<=RQd->ZMax && R>=R0 && R<=R1;
         } else return FALSE;
     }

     if ( (D = (B*B - 4*A*C))<0.0 ) return FALSE;

     D = sqrt(D);

     if ( B > 0  ) {
       T = -2*C/(B+D);
     } else {
       T = (-B+D)/(2*A);
     }
     if ( T>EPS && T<1-EPS )
     {
         Z = FZ + T*DZ;
         R = (Az*Z + Bz)*Z + Cz;

         if (Z>=RQd->ZMin && Z<=RQd->ZMax && R>=R0 && R<=R1) return TRUE;
     }

     if ( B > 0 ) {
       T = -(B+D)/(2*A);
     } else {
       T = 2*C/(-B+D);
     }
     if ( T>EPS && T<1-EPS )
     {
         Z = FZ + T*DZ;
         R = (Az*Z + Bz)*Z + Cz;

         return Z>=RQd->ZMin && Z<=RQd->ZMax && R>=R0 && R<=R1;
     }

     return FALSE;
}


/*******************************************************************************

Solve for ray segment and hierarchy of bounding volumes of elements, return
value is if there is hit or not.

LAST Modified: 2 Sep 1997

*******************************************************************************/

int RayHitVolumeBounds( VolumeBounds_t *Volume,double FX,double FY,double FZ,
           double DX,double DY,double DZ,double L)
{
  int i,j,n;

  if ( !RayHitBBox( &Volume->BBox,FX,FY,FZ,DX,DY,DZ ) ) return FALSE;

  if ( Volume->Left )
  {
    if ( RayHitVolumeBounds( Volume->Left,FX,FY,FZ,DX,DY,DZ,L ) ) return TRUE;
    return RayHitVolumeBounds( Volume->Right,FX,FY,FZ,DX,DY,DZ,L );
  }

  for( i=0; i<Volume->n; i++ )
  {
     j = Volume->Elements[i];
     n = RTElements[j].GeometryType;
     if ( (*RayHit[n] )( &RTElements[j],FX,FY,FZ,DX,DY,DZ,L) ) return TRUE;
  }

  return FALSE;
}

/*******************************************************************************

Solve for ray segment and geometry model intersection, return value is if
there is hit or not.

LAST Modified: 2 Sep 1997

*******************************************************************************/
int RayHitGeometry( double FX,double FY,double FZ, double DX,double DY,double DZ )
{
    double L;
    L = sqrt(DX*DX + DY*DY + DZ*DZ);
    return  RayHitVolumeBounds( VolumeBounds,FX,FY,FZ,DX,DY,DZ,L );

}

void VolumeBBox( VolumeBounds_t *Volume,Geometry_t *RTElements )
{
    double xMin,yMin,zMin,xMax,yMax,zMax,x,y,z;
    int i,j,k,N, NC;

    double U[] = {0.0,1.0,0.0,1.0}, V[] = {0.0,0.0,1.0,1.0};
 
    xMin = yMin = zMin =  1.0e20;
    xMax = yMax = zMax = -1.0e20;

    N = Volume->n;
    for( i=0; i<N; i++ )
    {
         k = Volume->Elements[i];
         switch(RTElements[k].GeometryType)
         {
            case GEOMETRY_LINE:
              NC = 2; break;
            case GEOMETRY_TRIANGLE:
              NC = 3; break;
            default:
              NC = 4; break;
        }
        for( j=0; j<NC; j++ )
          {
            x = FunctionValue( &RTElements[k], U[j],V[j], 0);
            y = FunctionValue( &RTElements[k], U[j],V[j], 1);
            z = FunctionValue( &RTElements[k], U[j],V[j], 2);

            xMin = MIN( x,xMin );
            yMin = MIN( y,yMin );
            zMin = MIN( z,zMin );

            xMax = MAX( x,xMax );
            yMax = MAX( y,yMax );
            zMax = MAX( z,zMax );
          }
    }

    Volume->BBox.XMin = xMin;
    Volume->BBox.XMax = xMax;

    Volume->BBox.YMin = yMin;
    Volume->BBox.YMax = yMax;

    Volume->BBox.ZMin = zMin;
    Volume->BBox.ZMax = zMax;
}



void VolumeDivide( VolumeBounds_t *Volume,int NBounds,Geometry_t *RT_Elements,int Level )
{
    double L1 = Volume->BBox.XMax - Volume->BBox.XMin;
    double L2 = Volume->BBox.YMax - Volume->BBox.YMin;
    double L3 = Volume->BBox.ZMax - Volume->BBox.ZMin;

    int NC;
    Geometry_t *Geom;

    double U[] = { 0.0,1.0,0.0,1.0 }, V[] = { 0.0,0.0,1.0,1.0 }, x,y,z;

    VolumeBounds_t *LeftVolume,*RightVolume;

    int i,j,k,n,N,left,right;

    static int count;
    double VolV,VolL,VolR;

    LeftVolume  = Volume->Left  = (VolumeBounds_t *)calloc(1,sizeof(VolumeBounds_t));
    RightVolume = Volume->Right = (VolumeBounds_t *)calloc(1,sizeof(VolumeBounds_t));

    LeftVolume->BBox  = Volume->BBox;
    RightVolume->BBox = Volume->BBox;

    if ( L1 > L2 &&  L1 > L3 )
    {
       LeftVolume->BBox.XMax = RightVolume->BBox.XMin = 
           (Volume->BBox.XMax - Volume->BBox.XMin)/2 + Volume->BBox.XMin;
    } else if ( L2 > L3 )
    {
       LeftVolume->BBox.YMax = RightVolume->BBox.YMin = 
           (Volume->BBox.YMax - Volume->BBox.YMin)/2 + Volume->BBox.YMin;
    } else 
    {
       LeftVolume->BBox.ZMax = RightVolume->BBox.ZMin = 
           (Volume->BBox.ZMax - Volume->BBox.ZMin)/2 + Volume->BBox.ZMin;
    }

    N = Volume->n;
    for( i=0; i<N; i++ )
    {
        k = Volume->Elements[i];
        left = right = FALSE;

        switch(RTElements[k].GeometryType)
        {
          case GEOMETRY_LINE:
              NC = 2; break;
          case GEOMETRY_TRIANGLE:
              NC = 3; break;
          default:
              NC = 4; break;
        }

        for( j=0; j<NC; j++ )
        {
           x = FunctionValue( &RTElements[k], U[j],V[j], 0);
           y = FunctionValue( &RTElements[k], U[j],V[j], 1);
           z = FunctionValue( &RTElements[k], U[j],V[j], 2);

           if ( (x >= LeftVolume->BBox.XMin) && (x <= LeftVolume->BBox.XMax) )
           if ( (y >= LeftVolume->BBox.YMin) && (y <= LeftVolume->BBox.YMax) )
           if ( (z >= LeftVolume->BBox.ZMin) && (z <= LeftVolume->BBox.ZMax) )  left = TRUE;

           if ( (x >= RightVolume->BBox.XMin) && (x <= RightVolume->BBox.XMax) )
           if ( (y >= RightVolume->BBox.YMin) && (y <= RightVolume->BBox.YMax) )
           if ( (z >= RightVolume->BBox.ZMin) && (z <= RightVolume->BBox.ZMax) ) right = TRUE;
        }
        if ( left )  LeftVolume->n++; else if ( right ) RightVolume->n++;
    }

    LeftVolume->Elements  = (int *)calloc( LeftVolume->n,sizeof(int) );
    LeftVolume->n = 0;

    RightVolume->Elements = (int *)calloc( RightVolume->n,sizeof(int) );
    RightVolume->n = 0;

    for( i=0; i<N; i++ )
    {
        k = Volume->Elements[i];
        left = right = FALSE;

        switch(RTElements[k].GeometryType)
        {
          case GEOMETRY_LINE:
            NC = 2; break;
          case GEOMETRY_TRIANGLE:
            NC = 3; break;
          default:
            NC = 4; break;
        }
        for( j=0; j<NC; j++ )
        {
           x = FunctionValue( &RTElements[k], U[j],V[j], 0);
           y = FunctionValue( &RTElements[k], U[j],V[j], 1);
           z = FunctionValue( &RTElements[k], U[j],V[j], 2);

           if ( (x >= LeftVolume->BBox.XMin) && (x <= LeftVolume->BBox.XMax) )
           if ( (y >= LeftVolume->BBox.YMin) && (y <= LeftVolume->BBox.YMax) )
           if ( (z >= LeftVolume->BBox.ZMin) && (z <= LeftVolume->BBox.ZMax) )  left = TRUE;

           if ( (x >= RightVolume->BBox.XMin) && (x <= RightVolume->BBox.XMax) )
           if ( (y >= RightVolume->BBox.YMin) && (y <= RightVolume->BBox.YMax) )
           if ( (z >= RightVolume->BBox.ZMin) && (z <= RightVolume->BBox.ZMax) ) right = TRUE;
        }
        if ( left )  LeftVolume->Elements[LeftVolume->n++] = k;
        else if ( right ) RightVolume->Elements[RightVolume->n++] = k;
    }

    VolumeBBox( LeftVolume,  RTElements );
    VolumeBBox( RightVolume, RTElements );

    if ( L1 > L2 && L1 > L3 )
    {
        VolV = L1;
        VolL = LeftVolume->BBox.XMax  - LeftVolume->BBox.XMin;
        VolR = RightVolume->BBox.XMax - RightVolume->BBox.XMin;
    } else if ( L2 > L3 ) {
        VolV = L2;
        VolL = LeftVolume->BBox.YMax  - LeftVolume->BBox.YMin;
        VolR = RightVolume->BBox.YMax - RightVolume->BBox.YMin;
    } else {
        VolV = L3;
        VolL = LeftVolume->BBox.ZMax  - LeftVolume->BBox.ZMin;
        VolR = RightVolume->BBox.ZMax - RightVolume->BBox.ZMin;
    }

    count += 2;
    if ( VolL-VolV > -1.0e-12 || VolR-VolV > -1.0e-12 )
    {
        free( LeftVolume->Elements );
        free( LeftVolume );

        free( RightVolume->Elements );
        free( RightVolume);

        Volume->Left = Volume->Right = NULL;
        count -= 2;

        return;
    }

    if ( LeftVolume->n >NBounds && Level<MAX_LEVEL )
      VolumeDivide( LeftVolume,NBounds,RTElements, Level+1 );

    if ( RightVolume->n>NBounds && Level<MAX_LEVEL )
      VolumeDivide( RightVolume,NBounds,RTElements,Level+1 );

    LeftVolume->BBox.XMin = LeftVolume->BBox.XMin - 
        0.001*(LeftVolume->BBox.XMax-LeftVolume->BBox.XMin);

    LeftVolume->BBox.XMax = LeftVolume->BBox.XMax +
        0.001*(LeftVolume->BBox.XMax-LeftVolume->BBox.XMin);

    LeftVolume->BBox.YMin = LeftVolume->BBox.YMin - 
        0.001*(LeftVolume->BBox.YMax-LeftVolume->BBox.YMin);

    LeftVolume->BBox.YMax = LeftVolume->BBox.YMax +
        0.001*(LeftVolume->BBox.YMax-LeftVolume->BBox.YMin);

    LeftVolume->BBox.ZMin = LeftVolume->BBox.ZMin - 
        0.001*(LeftVolume->BBox.ZMax-LeftVolume->BBox.ZMin);

    LeftVolume->BBox.ZMax = LeftVolume->BBox.ZMax + 
        0.001*(LeftVolume->BBox.ZMax-LeftVolume->BBox.ZMin);


    RightVolume->BBox.XMin = RightVolume->BBox.XMin - 
        0.001*(RightVolume->BBox.XMax-RightVolume->BBox.XMin);

    RightVolume->BBox.XMax = RightVolume->BBox.XMax + 
        0.001*(RightVolume->BBox.XMax-RightVolume->BBox.XMin);

    RightVolume->BBox.YMin = RightVolume->BBox.YMin -
        0.001*(RightVolume->BBox.YMax-RightVolume->BBox.YMin);

    RightVolume->BBox.YMax = RightVolume->BBox.YMax +
        0.001*(RightVolume->BBox.YMax-RightVolume->BBox.YMin);

    RightVolume->BBox.ZMin = RightVolume->BBox.ZMin -
        0.001*(RightVolume->BBox.ZMax-RightVolume->BBox.ZMin);

    RightVolume->BBox.ZMax = RightVolume->BBox.ZMax + 
        0.001*(RightVolume->BBox.ZMax-RightVolume->BBox.ZMin);
}


void InitVolumeBounds( int NBounds,int N,Geometry_t *RTElements )
{
    int i;

    VolumeBounds = (VolumeBounds_t *)calloc( 1,sizeof(VolumeBounds_t) );
    VolumeBounds->n = N;

    VolumeBounds->Elements = (int *)malloc( N*sizeof(int) );
    for( i=0;i<N; i++ ) VolumeBounds->Elements[i] = i;

    VolumeBBox( VolumeBounds,RTElements );
    VolumeDivide( VolumeBounds,NBounds,RTElements,0 );
}

void InitRayTracer(double eps)
{
     REPS = eps;

     RayHit[GEOMETRY_LINE]       = RayHitLine;
     RayHit[GEOMETRY_TRIANGLE]   = RayHitTriangle;
     RayHit[GEOMETRY_BILINEAR]   = RayHitBiLinear;
     RayHit[GEOMETRY_BICUBIC]    = RayHitBiCubic;
     RayHit[GEOMETRY_CIRCLE]     = RayHitCircle;
     RayHit[GEOMETRY_CYLINDER]   = RayHitCylinder;
     RayHit[GEOMETRY_ROTQUADRIC] = RayHitRotationalQuadric1;
}
