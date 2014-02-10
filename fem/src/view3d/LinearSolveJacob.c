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

static void LinearSolveGather( Geometry_t *Geom,double UpdateValue,double ViewFactorRowSum )
{
     GeometryList_t *Link = Geom->Link,*Diag=NULL;
     double Area = Geom->Area;

     while( Link )
     {
         ViewFactorRowSum += Link->ViewFactor;

         if ( Geom != Link->Entry )
         {
             UpdateValue += 0.85*Link->ViewFactor*Link->Entry->B;
         }
         else Diag = Link;

         Link = Link->Next;
     }

     if ( Geom->Flags & GEOMETRY_FLAG_LEAF )
     {
         Geom->Value = Geom->E + UpdateValue/ViewFactorRowSum;
         if ( Diag ) Geom->Value /= (1.0-0.85*Diag->ViewFactor/ViewFactorRowSum);

         return;
     }
    
     LinearSolveGather( Geom->Left,UpdateValue,ViewFactorRowSum );
     LinearSolveGather( Geom->Right,UpdateValue,ViewFactorRowSum );
}

static double LinearSolveUpdate( Geometry_t *Geom,double *Residual )
{
    double B = 0.0,Area = Geom->Area;

    if ( Geom->Flags & GEOMETRY_FLAG_LEAF )
    {
        *Residual += ABS(Geom->B - Geom->Value);
        Geom->B = Geom->Value;

        return Geom->B*Area;
    }

    B += LinearSolveUpdate( Geom->Left,Residual );
    B += LinearSolveUpdate( Geom->Right,Residual );

    Geom->B = B / Area;

    return B;
}

void IntegrateFromGeometry()
{
    int i,j,k,CE,n=0;
    double Residual,sum=0.0,*ViewFactorRowScale,*R,*P,*AP;

    FILE *fp;

    for( i=0; i<NGeomElem; i++ )
    {
        Geometry[i].Area = BiCubicArea( &Geometry[i] );
        Geometry[i].Flags |= GEOMETRY_FLAG_LEAF;
    }

    for( i=0; i<NGeomElem; i++ )
    for( j=i; j<NGeomElem; j++ ) ComputeViewFactors( &Geometry[i],&Geometry[j],0,0,&CE );

    for( i=8; i<12; i++ ) PutValue( &Geometry[i],1.0 );

    for( i=0; i<NGeomElem; i++ ) NumerateLeaves( &Geometry[i] );

    R  = (double *)calloc( NMAX,sizeof(double) );
    P  = (double *)calloc( NMAX,sizeof(double) );
    AP = (double *)calloc( NMAX,sizeof(double) );
    ViewFactorRowScale = (double *)calloc( NMAX,sizeof(double) );

    for( i=0; i<NGeomElem; i++ ) ComputeViewFactorRowSums( &Geometry[i],0.0,ViewFactorRowScale );

#ifdef OUTPUT_EQU
    fprintf( stderr, "NMAX: %d\n", NMAX );

    A = (double *)calloc( NMAX*NMAX,sizeof(double) );
    Y = (double *)calloc( NMAX,sizeof(double) );
    S = (double *)calloc( NMAX,sizeof(double) );

    for( i=0; i<NGeomElem; i++ ) MakeMatrix( &Geometry[i] );

    for( i=0; i<NMAX; i++ ) fprintf( stderr, "ROWSUM1: %g\n", R[i] );

    for( i=0; i<NMAX; i++ )
    {
        sum = 0.0;
        for( j=0; j<NMAX; j++ ) sum += A[i*NMAX+j];
        fprintf( stderr, "ROWSUM2: %g\n", sum );
    }

    fprintf( stdout, "%d\n", NMAX );

    for( i=0; i<NMAX; i++ ) fprintf( stdout, "%g\n",S[i] );
    for( i=0; i<NMAX; i++ ) fprintf( stdout, "%g\n",Y[i] );

    for( i=0; i<NMAX; i++ )
    {
        for( j=0; j<NMAX; j++ )
        {
            if ( i==j )
                A[i*NMAX+j] = S[i]*(1-0.85*A[i*NMAX+j]/R[i]);
            else
                A[i*NMAX+j] = -0.85*A[i*NMAX+j]*S[i]/R[i];

            fprintf( stdout, "%g\n", A[i*NMAX+j] );
        }

        Y[i] *= S[i];
        S[i] =  0.0;
    }
#endif

#ifndef JACOB
{
double T,second(),Alpha,Beta,RNorm,R1Norm;
T = second();

#ifndef OUTPUT_EQU
    for( i=0; i<NGeomElem; i++ ) LinearSolveResidual( &Geometry[i],0.0,ViewFactorRowScale,R );
    for( i=0; i<NGeomElem; i++ ) LinearSolveUpdateB( &Geometry[i],0.0 );
    for( i=0; i<NGeomElem; i++ ) LinearSolveUpdateRP( &Geometry[i],0.0,0.0,R,AP );
#else
    for( i=0; i<NMAX; i++ ) { R[i] = Y[i]; P[i] = Y[i]; }
#endif

    for( j=0; j<200; j++ )
    {
#ifndef OUTPUT_EQU
        for( i=0; i<NGeomElem; i++ ) LinearSolveGather( &Geometry[i],0.0,ViewFactorRowScale,P,AP );
#else

        for( i=0; i<NMAX; i++ )
        {
             AP[i] = 0.0;
             for( k=0; k<NMAX; k++ ) AP[i] += A[i*NMAX+k]*P[k];
        }
#endif

        RNorm = R1Norm = Alpha = Beta = 0.0;

        for( i=0; i<NMAX; i++ ) RNorm += R[i]*R[i];
        for( i=0; i<NMAX; i++ ) Alpha += P[i]*AP[i];
        Alpha = RNorm/Alpha;

        for( i=0; i<NMAX; i++ ) R1Norm += (R[i]-Alpha*AP[i])*(R[i]-Alpha*AP[i]);
        Beta  = R1Norm/RNorm;

        fprintf( stderr, "ITERATION: %d, RES: %g,%g,%g,%g\n", j,R1Norm,Alpha,Beta,RNorm );

#ifdef OUTPUT_EQU
        for( i=0; i<NMAX; i++ ) S[i] += Alpha*P[i];
        for( i=0; i<NMAX; i++ ) R[i] -= Alpha*AP[i];
        for( i=0; i<NMAX; i++ ) P[i] = R[i] + Beta*P[i];
#else
        for( i=0; i<NGeomElem; i++ ) LinearSolveUpdateB( &Geometry[i],Alpha );
        for( i=0; i<NGeomElem; i++ ) LinearSolveUpdateRP( &Geometry[i],Alpha,Beta,R,AP );
#endif

        if ( R1Norm < 1.0E-11 || j==199 )
        {
            sprintf( str, "b%d.p", j );
            fp = fopen( str, "w" );
 
            n = 0;
            for( i=0; i<NGeomElem; i++ ) PrintGeometry(&Geometry[i],&n);

            fprintf( fp, "%d %d 1 1\n", n,n/4 );
            for( i=0; i<n; i++ )fprintf( fp, "%g %g %g\n",XX[i],YY[i],ZZ[i] );
            for( i=0; i<n/4; i++ )
            {
                fprintf( fp, "1 404 %d %d %d %d\n",4*i,4*i+1,4*i+2,4*i+3 );
            }

            for( i=0; i<n; i++ ) fprintf( fp, "%g\n", LL[i] );
            fclose( fp );

            break;
        }
    }
fprintf( stderr, "LINSOLVE TIME: %g\n", second()-T );
}
#else
{
double T,second();
T = second();
    for( j=0; j<200; j++ )
    {
        for( i=0; i<NGeomElem; i++ ) LinearSolveGather( &Geometry[i], 0.0, 0.0 );

        Residual = 0.0;
        for( i=0; i<NGeomElem; i++ ) LinearSolveUpdate( &Geometry[i],&Residual );
        fprintf( stderr, "ITERATION: %d, RES: %g\n", j,Residual );

        if ( Residual < 1.0E-11 || j==199  )
        {
            sprintf( str, "b%d.p", j );
            fp = fopen( str, "w" );
 
            n = 0;
            for( i=0; i<NGeomElem; i++ ) PrintGeometry(&Geometry[i],&n);

            fprintf( fp, "%d %d 1 1\n", n,n/4 );
            for( i=0; i<n; i++ )fprintf( fp, "%g %g %g\n",XX[i],YY[i],ZZ[i] );
            for( i=0; i<n/4; i++ )
            {
                fprintf( fp, "1 404 %d %d %d %d\n",4*i,4*i+1,4*i+2,4*i+3 );
            }

            for( i=0; i<n; i++ ) fprintf( fp, "%g\n", LL[i] );
            fclose( fp );

            break;
        }
    }
fprintf( stderr, "LINSOLVE TIME: %g\n", second()-T );
}
#endif

