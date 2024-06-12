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

static double LinearSolveGather( Geometry_t *Geom,double UpdateValue,double *ViewFactorRowScale,double *Residual )
{
    GeometryList_t *Link=Geom->Link,*Diag=NULL;
    double Area = Geom->Area, B=0.0;

    while( Link )
    {
        if ( Geom != Link->Entry )
        {
            UpdateValue += Link->ViewFactor*Link->Entry->M*Link->Entry->B;
        }
        else Diag = Link;

        Link = Link->Next;
    }

    if ( Geom->Flags & GEOMETRY_FLAG_LEAF )
    {
        double E=Geom->E,M=1.0;
        int N = Geom->N;

        B = E + UpdateValue*ViewFactorRowScale[N];

        if ( Diag ) M -= Diag->ViewFactor*Geom->M*ViewFactorRowScale[N];

        *Residual += ABS(B-M*Geom->B);
        B /= M;

        Geom->B = B;
        return Geom->B*Area;
    }
    
    B += LinearSolveGather( Geom->Left,UpdateValue,ViewFactorRowScale,Residual );
    B += LinearSolveGather( Geom->Right,UpdateValue,ViewFactorRowScale,Residual );

    Geom->B = B / Area;

    return B;
}

void LinearSolveGaussSeidel( Geometry_t *Geometry,int NGeom,double *ViewFactorRowScale )
{
    int i,iter,MAX_GAUSS_SEIDEL_ITER=200;
    double Residual,FirstResidual;

double T,second();
T = second();

    FirstResidual = 0.0;
    for( i=0; i<NGeom; i++ )
    {
       LinearSolveGather( &Geometry[i],0.0,ViewFactorRowScale,&FirstResidual );
    }

    for( iter=0; iter<MAX_GAUSS_SEIDEL_ITER; iter++ )
    {
       Residual = 0.0;
       for( i=0; i<NGeom; i++ )
       {
          LinearSolveGather( &Geometry[i],0.0,ViewFactorRowScale,&Residual );
       }
       if ( Residual < 1.0E-08 ) break;
fprintf( stderr, "LINSOLVE TIME: %g,RES (first,latest): %g, %g, ITER: %d\n", second()-T,FirstResidual,Residual,iter );
    }
fprintf( stderr, "LINSOLVE TIME: %g,RES (first,last): %g, %g, ITER: %d\n", second()-T,FirstResidual,Residual,iter );
}
