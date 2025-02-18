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

/*******************************************************************************
 *
 * Main module for element model descriptions & utility routines.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 12 Apr 1996
 *
 * Modification history:
 *
 * 12 Apr 1996, 
 * Juha R
 *
 *
 ******************************************************************************/

#define MODULE_ELEMENT_MATH

#include "../elmerpost.h"
#include "elements.h"


/*******************************************************************************
 *
 *     Name:        elm_element_gradient_2D
 *
 *     Purpose:     Compute elementwise gradient vector of given quantity for
 *                  surface elements.
 *
 *     Parameters:
 *
 *         Input:   (element_t *) structure describing the element
 *                  (double *)    node values of the quantity.
 *
 *         Output:  Gradient vector (double *)GX
 *                                  (double *)GY
 *                                  (double *)GZ at nodes.
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_gradient_2D( element_model_t *model, element_t *elm, double *F,
       double *GX, double *GY, double *GZ, int NPoints, double *NU, double *NV )
{
    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    int i,n;

    double u,v,a,detA,Auv,Avu,Auu,Avv,dXdU,dYdU,dXdV,dYdV,dZdU,dZdV,dFdU,dFdV;

    element_type_t *elmt = elm->ElementType;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    for( i=0; i<NPoints; i++ )
    {
        u = NU[i];
        v = NV[i];

        dXdU = (*dNdU)(X,u,v);
        dYdU = (*dNdU)(Y,u,v);
        dZdU = (*dNdU)(Z,u,v);
        dFdU = (*dNdU)(F,u,v);

        dXdV = (*dNdV)(X,u,v);
        dYdV = (*dNdV)(Y,u,v);
        dZdV = (*dNdV)(Z,u,v);
        dFdV = (*dNdV)(F,u,v);

#if 0
        Avv = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
        Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
        Auu = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

        detA = Auu*Avv - Auv*Auv;
        if ( ABS(detA) < AEPS ) return FALSE;

        Auv = -Auv / detA;           /*   ij  */
        Auu =  Auu / detA;           /*  a    */
        Avv =  Avv / detA;

        a = dFdU;
        dFdU = Auu*a + Auv*dFdV;          /* raise index of the surface */
        dFdV = Auv*a + Avv*dFdV;          /*   vector (@f/@u,@f/@v)     */

        /* transform to global cartesian frame */
        GX[i] = dXdU*dFdU + dXdV*dFdV;
        GY[i] = dYdU*dFdU + dYdV*dFdV;
        GZ[i] = dZdU*dFdU + dZdV*dFdV;       
#else
        Auu = dYdV;
        Auv = dYdU;
        Avu = dXdV;
        Avv = dXdU;
        detA = Auu*Avv - Auv*Avu;

        Auu =  Auu / detA;
        Auv = -Auv / detA;
        Avu = -Avu / detA;
        Avv =  Avv / detA;

        GX[i] = Auu*dFdU + Auv*dFdV;
        GY[i] = Avu*dFdU + Avv*dFdV;
        GZ[i] = 0.0;
#endif
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_element_gradient_3D
 *
 *     Purpose:     Compute elementwise gradient vector of given quantity
 *                  for volume elements.
 *
 *     Parameters:
 *
 *         Input:   (element_t *) structure describing the element
 *                  (double *)    node values of the quantity.
 *
 *         Output:  Gradient vector (double *)GX
 *                                  (double *)GY
 *                                  (double *)GZ at nodes.
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_gradient_3D( element_model_t *model, element_t *elm, double *F,
   double *GX, double *GY, double *GZ, int NPoints, double *NU, double *NV, double *NW )
{
    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    double u,v,w,a,detJ;

    double Jux,Juy,Juz,Jvx,Jvy,Jvz,Jwx,Jwy,Jwz;
    double Kxu,Kxv,Kxw,Kyu,Kyv,Kyw,Kzu,Kzv,Kzw;
    double dFdU, dFdV, dFdW;

    int i,n;

    element_type_t *elmt = elm->ElementType;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;
    double (*dNdW)() = elmt->PartialW;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    for( i=0; i<NPoints; i++ )
    {
        u = NU[i];
        v = NV[i];
        w = NW[i];

        Jux  = (*dNdU)( X,u,v,w );
        Juy  = (*dNdU)( Y,u,v,w );
        Juz  = (*dNdU)( Z,u,v,w );
        dFdU = (*dNdU)( F,u,v,w );

        Jvx  = (*dNdV)( X,u,v,w );
        Jvy  = (*dNdV)( Y,u,v,w );
        Jvz  = (*dNdV)( Z,u,v,w );
        dFdV = (*dNdV)( F,u,v,w );

        Jwx  = (*dNdW)( X,u,v,w );
        Jwy  = (*dNdW)( Y,u,v,w );
        Jwz  = (*dNdW)( Z,u,v,w );
        dFdW = (*dNdW)( F,u,v,w );

        detJ  = Jux*(Jvy*Jwz - Jvz*Jwy);
        detJ += Juy*(Jvz*Jwx - Jvx*Jwz);
        detJ += Juz*(Jvx*Jwy - Jvy*Jwx);
        if ( ABS(detJ) < AEPS ) return FALSE;
        detJ = 1 / detJ;

        Kxu = (Jvy*Jwz - Jvz*Jwy)*detJ;
        Kyu = (Jvz*Jwx - Jvx*Jwz)*detJ;
        Kzu = (Jvx*Jwy - Jvy*Jwx)*detJ;

        Kxv = (Jwy*Juz - Juy*Jwz)*detJ;
        Kyv = (Jux*Jwz - Jwx*Juz)*detJ;
        Kzv = (Juy*Jwx - Jux*Jwy)*detJ;

        Kxw = (Juy*Jvz - Jvy*Juz)*detJ;
        Kyw = (Jvx*Juz - Jux*Jvz)*detJ;
        Kzw = (Jux*Jvy - Jvx*Juy)*detJ;

        /* transform to global cartesian frame */
        GX[i] = Kxu*dFdU + Kxv*dFdV + Kxw*dFdW; 
        GY[i] = Kyu*dFdU + Kyv*dFdV + Kyw*dFdW;
        GZ[i] = Kzu*dFdU + Kzv*dFdV + Kzw*dFdW;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_gradient
 *
 *     Purpose:     Compute gradient vector of given quantity.
 *
 *     Parameters:
 *
 *         Input:
 *
 *         Output:
 *
 *   Return value:  void
 *
 ******************************************************************************/
VARIABLE *elm_gradient( VARIABLE *A )
{
     static double dFdX[ELM_MAX_ELEMENT_NODES];
     static double dFdY[ELM_MAX_ELEMENT_NODES];
     static double dFdZ[ELM_MAX_ELEMENT_NODES];
     static double F[ELM_MAX_ELEMENT_NODES];

     element_type_t *elmt;
     element_t *elm;

     element_model_t *model = CurrentObject->ElementModel;  /*** TODO: FIX ***/

     unsigned char *References = (unsigned char *)malloc( model->NofNodes*sizeof(unsigned char) );
 
     double *af = MATR( A );
 
     int i,j,k,n,status,dim, cols;

     VARIABLE *res;

     double *rx, *ry, *rz;

     cols = NCOL(A) / model->NofNodes;
     if ( NROW(A) != 1 )
     {
         free( References );
         error( "I can only compute the gradient of a scalar variable\n" );
     }
 
     res = (VARIABLE *)var_temp_new( TYPE_DOUBLE, 3, cols*model->NofNodes );
     rx  = &M(res,0,0);
     ry  = &M(res,1,0);
     rz  = &M(res,2,0);

     dim = 1;
     for( k=0; k<model->NofElements; k++ )
      if ( model->Elements[k].ElementType->ElementCode>=500 ) {
        dim=3;
        break;
      } else if ( model->Elements[k].ElementType->ElementCode>=300 ) 
        dim=2;

     for( k=0; k<cols; k++ )
     {
         for( j=0; j<model->NofNodes; j++ )
         {
             n = k*model->NofNodes + j;
             rx[n] = 0.0;
             ry[n] = 0.0;
             rz[n] = 0.0;
             References[j] = 0;
         }
  
         for( i=0; i<model->NofElements; i++ )
         {
             elm  = &model->Elements[i];
             elmt = elm->ElementType;
 
             for( j=0; j<elmt->NumberOfNodes; j++ )
             {
                 n = k*model->NofNodes + elm->Topology[j];
                 F[j] = af[n];
                 dFdX[j] = af[n];
                 dFdY[j] = af[n];
                 dFdZ[j] = af[n];
             }
             if ( dim==3 && elmt->PartialW )
                 status = elm_element_gradient_3D(
                     model, elm,F,dFdX,dFdY,dFdZ,elmt->NumberOfNodes,elmt->NodeU,elmt->NodeV,elmt->NodeW
                   );
             else if ( dim==2 && elmt->PartialV )
                 status = elm_element_gradient_2D(
                     model, elm,F,dFdX,dFdY,dFdZ,elmt->NumberOfNodes,elmt->NodeU,elmt->NodeV 
                   ); 
             else continue;
 
             if ( !status ) 
             {
                 fprintf( stderr, "ELEMENT [%d]: Jacobian is singular.\n",i );
                 continue;
             }
 
             for( j=0; j<elmt->NumberOfNodes; j++ )
             {
                 n = elm->Topology[j];
                 References[n]++;

                 n += k*model->NofNodes;
                 rx[n] += dFdX[j];
                 ry[n] += dFdY[j];
                 rz[n] += dFdZ[j];
             }
         }

         for( j=0; j<model->NofNodes; j++ )
         {
             if ( References[j] != 0 )
             {
                 n = k*model->NofNodes + j;
                 rx[n] /= References[j];
                 ry[n] /= References[j];
                 rz[n] /= References[j];
             } 
         } 
     }

     free( References );
     return res;
}

/*******************************************************************************
 *
 *     Name:        elm_element_rotor_2D
 *
 *     Purpose:     Compute elementwise rotor of given quantity for a surface
 *                  element.
 *
 *     Parameters:
 *
 *         Input:   (element_t *)elm structure describing the element
 *                  (double *)FX      node values of vector x component.
 *                  (double *)FY      node values of vector x component.
 *                  (double *)FZ      node values of vector x component.
 *
 *         Output:  Rotor (double *)R
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_rotor_2D
     (
       element_model_t *model, element_t *elm, double *FX, double *FY, double *FZ, double *R
     )
{
    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    static double Fu[ELM_MAX_ELEMENT_NODES];
    static double Fv[ELM_MAX_ELEMENT_NODES];

    double u,v,w,a,detA;

    double Auu,Auv,Avu,Avv;
    double dXdU,dYdU,dZdU;
    double dXdV,dYdV,dZdV;

    double dFudU, dFudV;
    double dFvdU, dFvdV;
    double dFwdU, dFwdV;

    int i,n;

    element_type_t *elmt = elm->ElementType;

    double *NodeU = elmt->NodeU;
    double *NodeV = elmt->NodeV;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];

        dXdU = (*dNdU)( X,u,v );
        dYdU = (*dNdU)( Y,u,v );
        dZdU = (*dNdU)( Z,u,v );

        dXdV = (*dNdV)( X,u,v );
        dYdV = (*dNdV)( Y,u,v );
        dZdV = (*dNdV)( Z,u,v );

        Fu[i] = FX[i]*dXdU + FY[i]*dYdU + FZ[i]*dZdU;
        Fv[i] = FX[i]*dXdV + FY[i]*dYdV + FZ[i]*dZdV;
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];

        dXdU  = (*dNdU)( X,u,v );
        dYdU  = (*dNdU)( Y,u,v );
        dZdU  = (*dNdU)( Z,u,v );

        dXdV  = (*dNdV)( X,u,v );
        dYdV  = (*dNdV)( Y,u,v );
        dZdV  = (*dNdV)( Z,u,v );

        dFvdU = (*dNdU)( Fv,u,v );
        dFudV = (*dNdV)( Fu,u,v );

        Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU;
        Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV;
        Avu = dXdV*dXdU + dYdV*dYdU + dZdV*dZdU;
        Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

        detA = Auu*Avv - Auv*Avu;

        if ( ABS(detA) < AEPS ) return FALSE;

        R[i] = (dFvdU - dFudV) / sqrt( detA );
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_rotor_2D
 *
 *     Purpose:     Compute rotor of given quantity for surface elements.
 *
 *     Parameters:
 *
 *         Input:
 *
 *         Output:
 *
 *   Return value:  void
 *
 ******************************************************************************/
VARIABLE *elm_rotor_2D( VARIABLE *A )
{
    static double FX[ELM_MAX_ELEMENT_NODES];
    static double FY[ELM_MAX_ELEMENT_NODES];
    static double FZ[ELM_MAX_ELEMENT_NODES];

    static double R[ELM_MAX_ELEMENT_NODES];

    element_type_t *elmt;
    element_t *elm;

    element_model_t *model = CurrentObject->ElementModel;  /*** TODO: FIX ***/

    unsigned char *References = (unsigned char *)malloc( model->NofNodes*sizeof(unsigned char) );

    int i,j,k,n,cols;

    double *ax = &M(A,0,0);
    double *ay = &M(A,1,0);
    double *az = &M(A,2,0);

    VARIABLE *res;
    double *rf;

    cols = NCOL(A) / model->NofNodes;
    if ( NROW(A) != 3 )
    {
        free( References );
        error( "curl 2D: vector variable needed as input.\n" );
    }

    res = var_temp_new( TYPE_DOUBLE, 1, cols*model->NofNodes );
    rf = MATR(res);
 
    for( k=0; k<cols; k++ )
    {
        for( j=0; j<model->NofNodes; j++ )
        {
            rf[k*model->NofNodes+j] = 0.0;
            References[j] = 0;
        }
 
        for( i=0; i<model->NofElements; i++ )
        {
            elm  = &model->Elements[i];
            elmt = elm->ElementType;

            if ( elmt->ElementCode < 300 || elmt->ElementCode > 500 ) continue;

            for( j=0; j<elmt->NumberOfNodes; j++ )
            {
                n = k*model->NofNodes + elm->Topology[j];
                FX[j] = ax[n];
                FY[j] = ay[n];
                FZ[j] = az[n];
            }

            if ( !elm_element_rotor_2D( model,elm,FX,FY,FZ,R ) ) 
            {
                fprintf( stderr, "ELEMENT [%d]: Jacobian is singular.\n",i );
                continue;
            }

            for( j=0; j<elmt->NumberOfNodes; j++ )
            {
                n = elm->Topology[j];
                References[n]++;

                rf[k*model->NofNodes+n] += R[j];
            }
        }

        for( j=0; j<model->NofNodes; j++ )
        {
            if ( References[j] != 0 )
            {
                n = k*model->NofNodes + j;
                rf[n] /= References[j];
            } 
        } 
    }

    free( References );
    return res;
}

/*******************************************************************************
 *
 *     Name:        elm_element_rotor_3D
 *
 *     Purpose:     Compute elementwise rotor of given quantity for a volume
 *                  element.
 *
 *     Parameters:
 *
 *         Input:   (element_t *)elm structure describing the element
 *                  (double *)FX node values of vector x component.
 *                  (double *)FY node values of vector y component.
 *                  (double *)FZ node values of vector z component.
 *
 *         Output:  Rotor vector at nodes (double *)RX
 *                                        (double *)RY
 *                                        (double *)RZ
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_rotor_3D( element_model_t *model, element_t *elm,
                            double *FX, double *FY, double *FZ,
                            double *RX, double *RY, double *RZ )
{
    static double Fu[ELM_MAX_ELEMENT_NODES];
    static double Fv[ELM_MAX_ELEMENT_NODES];
    static double Fw[ELM_MAX_ELEMENT_NODES];

    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    double u,v,w,a,detA;

    double Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww;
    double dXdU,dYdU,dZdU;
    double dXdV,dYdV,dZdV;
    double dXdW,dYdW,dZdW;

    double dFudU, dFudV, dFudW;
    double dFvdU, dFvdV, dFvdW;
    double dFwdU, dFwdV, dFwdW;

    int i,n;

    element_type_t *elmt = elm->ElementType;

    double *NodeU = elmt->NodeU;
    double *NodeV = elmt->NodeV;
    double *NodeW = elmt->NodeW;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;
    double (*dNdW)() = elmt->PartialW;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];
        w = NodeW[i];

        dXdU = (*dNdU)( X,u,v,w );
        dYdU = (*dNdU)( Y,u,v,w );
        dZdU = (*dNdU)( Z,u,v,w );

        dXdV = (*dNdV)( X,u,v,w );
        dYdV = (*dNdV)( Y,u,v,w );
        dZdV = (*dNdV)( Z,u,v,w );

        dXdW = (*dNdW)( X,u,v,w );
        dYdW = (*dNdW)( Y,u,v,w );
        dZdW = (*dNdW)( Z,u,v,w );

        Fu[i] = FX[i]*dXdU + FY[i]*dYdU + FZ[i]*dZdU;
        Fv[i] = FX[i]*dXdV + FY[i]*dYdV + FZ[i]*dZdV;
        Fw[i] = FX[i]*dXdW + FY[i]*dYdW + FZ[i]*dZdW;
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];
        w = NodeW[i];

        dXdU  = (*dNdU)( X,u,v,w );
        dYdU  = (*dNdU)( Y,u,v,w );
        dZdU  = (*dNdU)( Z,u,v,w );

        dFudU = (*dNdU)( Fu,u,v,w );
        dFvdU = (*dNdU)( Fv,u,v,w );
        dFwdU = (*dNdU)( Fw,u,v,w );

        dXdV  = (*dNdV)( X,u,v,w );
        dYdV  = (*dNdV)( Y,u,v,w );
        dZdV  = (*dNdV)( Z,u,v,w );

        dFudV = (*dNdV)( Fu,u,v,w );
        dFvdV = (*dNdV)( Fv,u,v,w );
        dFwdV = (*dNdV)( Fw,u,v,w );

        dXdW  = (*dNdW)( X,u,v,w );
        dYdW  = (*dNdW)( Y,u,v,w );
        dZdW  = (*dNdW)( Z,u,v,w );

        dFudW = (*dNdW)( Fu,u,v,w );
        dFvdW = (*dNdW)( Fv,u,v,w );
        dFwdW = (*dNdW)( Fw,u,v,w );

        Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU;
        Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV;
        Auw = dXdU*dXdW + dYdU*dYdW + dZdU*dZdW;
        Avu = dXdV*dXdU + dYdV*dYdU + dZdV*dZdU;
        Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;
        Avw = dXdV*dXdW + dYdV*dYdW + dZdV*dZdW;
        Awu = dXdW*dXdU + dYdW*dYdU + dZdW*dZdU;
        Awv = dXdW*dXdV + dYdW*dYdV + dZdW*dZdV;
        Aww = dXdW*dXdW + dYdW*dYdW + dZdW*dZdW;

        detA  = Auu*(Avv*Aww - Avw*Awv);
        detA += Auv*(Avw*Awu - Avu*Aww);
        detA += Auw*(Avu*Awv - Avv*Awu);

        if ( ABS(detA) < AEPS ) return FALSE;

        u = (dFwdV - dFvdW) / sqrt( detA );
        v = (dFudW - dFwdU) / sqrt( detA );
        w = (dFvdU - dFudV) / sqrt( detA );

        RX[i] = u*dXdU + v*dXdV + w*dXdW;
        RY[i] = u*dYdU + v*dYdV + w*dYdW;
        RZ[i] = u*dZdU + v*dZdV + w*dZdW;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_rotor_3D
 *
 *     Purpose:     Compute rotor of given quantity for volume elements.
 *
 *     Parameters:
 *
 *         Input:
 *
 *         Output:
 *
 *   Return value:   void
 *
 ******************************************************************************/
VARIABLE *elm_rotor_3D( VARIABLE *A )
{
    static double RotorX[ELM_MAX_ELEMENT_NODES];
    static double RotorY[ELM_MAX_ELEMENT_NODES];
    static double RotorZ[ELM_MAX_ELEMENT_NODES];

    static double FX[ELM_MAX_ELEMENT_NODES];
    static double FY[ELM_MAX_ELEMENT_NODES];
    static double FZ[ELM_MAX_ELEMENT_NODES];

    element_type_t *elmt;
    element_t *elm;
    element_model_t *model = CurrentObject->ElementModel;  /*** TODO: FIX ***/

    unsigned char *References = (unsigned char *)malloc( model->NofNodes*sizeof(unsigned char) );

    int i,j,k,n, cols;

    double *ax = &M(A,0,0);
    double *ay = &M(A,1,0);
    double *az = &M(A,2,0);

    VARIABLE *res;

    double *rx, *ry, *rz;

    cols = NCOL(A) / model->NofNodes;
    if ( NROW(A) != 3 )
    {
        free( References );
        error( "curl: vector variable needed as input.\n" );
    }

    res = var_temp_new( TYPE_DOUBLE, 3, cols*model->NofNodes );
    rx = &M(res,0,0);
    ry = &M(res,1,0);
    rz = &M(res,2,0);

    for( k=0; k<cols; k++ )
    {
        for( j=0; j<model->NofNodes; j++ )
        {
            n = k*model->NofNodes  + j;
            rx[n] = 0.0;
            ry[n] = 0.0;
            rz[n] = 0.0;
            References[j] = 0;
        }
 
        for( i=0; i<model->NofElements; i++ )
        {
            elm  = &model->Elements[i];
            elmt = elm->ElementType;

            if ( elmt->ElementCode <500 ) continue;

            for( j=0; j<elmt->NumberOfNodes; j++ )
            {
                n = k*model->NofNodes + elm->Topology[j];
                FX[j] = ax[n];
                FY[j] = ay[n];
                FZ[j] = az[n];
            }

            if ( !elm_element_rotor_3D( model,elm,FX,FY,FZ,RotorX,RotorY,RotorZ ) ) 
            {
                fprintf( stderr, "ELEMENT [%d]: Jacobian is singular.\n",i );
                continue;
            }

            for( j=0; j<elmt->NumberOfNodes; j++ )
            {
                n = elm->Topology[j];
                References[n]++;

                n += k*model->NofNodes;
                rx[n] += RotorX[j];
                ry[n] += RotorY[j];
                rz[n] += RotorZ[j];
            }
        }

        for( j=0; j<model->NofNodes; j++ )
        {
            if ( References[j] != 0 )
            {
                n = k*model->NofNodes + j;
                rx[n] /= References[j];
                ry[n] /= References[j];
                rz[n] /= References[j];
            }
        }
    }

    free( References );
    return res;
}

/*******************************************************************************
 *
 *     Name:        elm_element_divergence_2D
 *
 *     Purpose:     Compute elementwise divergence of given quantity for a surface
 *                  element.
 *
 *     Parameters:
 *
 *         Input:   (element_t *)elm structure describing the element
 *                  (double *)FX      node values of vector x component.
 *                  (double *)FY      node values of vector x component.
 *                  (double *)FZ      node values of vector x component.
 *
 *         Output:  Rotor (double *)D
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_divergence_2D
    (
       element_model_t *model, element_t *elm, double *FX, double *FY, double *FZ, double *D
    )
{
    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    static double Fu[ELM_MAX_ELEMENT_NODES];
    static double Fv[ELM_MAX_ELEMENT_NODES];

    static double s,detA[ELM_MAX_ELEMENT_NODES];

    double u,v,w,a;

    double Auu,Auv,Avu,Avv;
    double Buu,Buv,Bvu,Bvv;

    double dXdU,dYdU,dZdU;
    double dXdV,dYdV,dZdV;

    double dFudU, dFudV;
    double dFvdU, dFvdV;
    double dFwdU, dFwdV;

    int i,n;

    element_type_t *elmt = elm->ElementType;

    double *NodeU = elmt->NodeU;
    double *NodeV = elmt->NodeV;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];

        dXdU = (*dNdU)( X,u,v );
        dYdU = (*dNdU)( Y,u,v );
        dZdU = (*dNdU)( Z,u,v );

        dXdV = (*dNdV)( X,u,v );
        dYdV = (*dNdV)( Y,u,v );
        dZdV = (*dNdV)( Z,u,v );

        Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU;
        Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV;
        Avu = dXdV*dXdU + dYdV*dYdU + dZdV*dZdU;
        Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

        s = Auu*Avv - Auv*Avu;
        if ( ABS(s) < AEPS ) return FALSE;

        s = 1 / s;

        Buu =  Avv * s;
        Buv = -Auv * s;
        Bvu = -Avu * s;
        Bvv =  Auu * s;

        u = FX[i]*dXdU + FY[i]*dYdU + FZ[i]*dZdU;
        v = FX[i]*dXdV + FY[i]*dYdV + FZ[i]*dZdV;

        s = detA[i] = sqrt( s );
        Fu[i] = (Buu*u + Buv*v) / s;
        Fv[i] = (Bvu*u + Bvv*v) / s;
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];

        dFudU = (*dNdU)( Fu,u,v );
        dFvdV = (*dNdV)( Fv,u,v );

        D[i] = (dFudU + dFvdV) * detA[i];
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_element_divergence_3D
 *
 *     Purpose:     Compute elementwise divergence of given quantity for a volume
 *                  element.
 *
 *     Parameters:
 *
 *         Input:   (element_t *)elm structure describing the element
 *                  (double *)FX      node values of vector x component.
 *                  (double *)FY      node values of vector x component.
 *                  (double *)FZ      node values of vector x component.
 *
 *         Output:  Rotor (double *)D
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_divergence_3D
    (
       element_model_t *model, element_t *elm, double *FX, double *FY, double *FZ, double *D
    )
{
    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    static double Fu[ELM_MAX_ELEMENT_NODES];
    static double Fv[ELM_MAX_ELEMENT_NODES];
    static double Fw[ELM_MAX_ELEMENT_NODES];
    static double s,detA[ELM_MAX_ELEMENT_NODES];

    double u,v,w,a;

    double Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww;
    double Buu,Buv,Buw,Bvu,Bvv,Bvw,Bwu,Bwv,Bww;

    double dXdU,dXdV,dXdW;
    double dYdU,dYdV,dYdW;
    double dZdU,dZdV,dZdW;

    double dFudU, dFvdV, dFwdW;

    int i,n;

    element_type_t *elmt = elm->ElementType;

    double *NodeU = elmt->NodeU;
    double *NodeV = elmt->NodeV;
    double *NodeW = elmt->NodeW;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;
    double (*dNdW)() = elmt->PartialW;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];
        w = NodeW[i];

        dXdU = (*dNdU)( X,u,v,w );
        dYdU = (*dNdU)( Y,u,v,w );
        dZdU = (*dNdU)( Z,u,v,w );

        dXdV = (*dNdV)( X,u,v,w );
        dYdV = (*dNdV)( Y,u,v,w );
        dZdV = (*dNdV)( Z,u,v,w );

        dXdW = (*dNdW)( X,u,v,w );
        dYdW = (*dNdW)( Y,u,v,w );
        dZdW = (*dNdW)( Z,u,v,w );

        Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU;
        Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV;
        Auw = dXdU*dXdW + dYdU*dYdW + dZdU*dZdW;

        Avu = dXdV*dXdU + dYdV*dYdU + dZdV*dZdU;
        Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;
        Avw = dXdV*dXdW + dYdV*dYdW + dZdV*dZdW;

        Awu = dXdW*dXdU + dYdW*dYdU + dZdW*dZdU;
        Awv = dXdW*dXdV + dYdW*dYdV + dZdW*dZdV;
        Aww = dXdW*dXdW + dYdW*dYdW + dZdW*dZdW;

        s  = Auu*(Avv*Aww - Avw*Awv);
        s += Auv*(Avw*Awu - Avu*Aww);
        s += Auw*(Avu*Awv - Avv*Awu);
        if ( ABS(s) < AEPS ) return FALSE;

        s = 1 / s;

        Buu = (Avv*Aww - Avw*Awv) * s;
        Bvu = (Avw*Awu - Avu*Aww) * s;
        Bwu = (Avu*Awv - Avv*Awu) * s;

        Buv = (Awv*Auw - Auv*Aww) * s;
        Bvv = (Auu*Aww - Awu*Auw) * s;
        Bwv = (Auv*Awu - Auu*Awv) * s;

        Buw = (Auv*Avw - Avv*Auw) * s;
        Bvw = (Avu*Auw - Auu*Avw) * s;
        Bww = (Auu*Avv - Avu*Auv) * s;

        u = FX[i]*dXdU + FY[i]*dYdU + FZ[i]*dZdU;
        v = FX[i]*dXdV + FY[i]*dYdV + FZ[i]*dZdV;
        w = FX[i]*dXdW + FY[i]*dYdW + FZ[i]*dZdV;

        s = detA[i] = sqrt( s );
        Fu[i] = (Buu*u + Buv*v + Buw*w) / s;
        Fv[i] = (Bvu*u + Bvv*v + Bvw*w) / s;
        Fw[i] = (Bwu*u + Bwv*v + Bww*w) / s;
    }

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        u = NodeU[i];
        v = NodeV[i];
        w = NodeV[i];

        D[i]  = (*dNdU)( Fu,u,v,w );
        D[i] += (*dNdV)( Fv,u,v,w );
        D[i] += (*dNdW)( Fw,u,v,w );
        D[i] *= detA[i];
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_divergence
 *
 *     Purpose:     Compute divergence of given quantity.
 *
 *     Parameters:
 *
 *         Input:
 *
 *         Output:
 *
 *   Return value:  void
 *
 ******************************************************************************/
VARIABLE *elm_divergence( VARIABLE *A )
{
    static double FX[ELM_MAX_ELEMENT_NODES];
    static double FY[ELM_MAX_ELEMENT_NODES];
    static double FZ[ELM_MAX_ELEMENT_NODES];

    static double D[ELM_MAX_ELEMENT_NODES];

    element_type_t *elmt;
    element_t *elm;
    element_model_t *model = CurrentObject->ElementModel;  /*** TODO: FIX ***/

    unsigned char *References = (unsigned char *)malloc( model->NofNodes*sizeof(unsigned char) );

    int i,j,k,n,status,dim,cols;

    double *ax = &M(A,0,0);
    double *ay = &M(A,1,0);
    double *az = &M(A,2,0);

    VARIABLE *res;
    double *rf;

     cols = NCOL(A) / model->NofNodes;
    if ( NROW(A) != 3 )
    {
        free( References );
        error( "div: vector variable needed as input.\n" );
    }

    dim = 1;
    for( k=0; k<model->NofElements; k++ )
     if ( model->Elements[k].ElementType->ElementCode>=500 ) {
       dim=3;
       break;
     } else if ( model->Elements[k].ElementType->ElementCode>=300 ) 
       dim=2;

    res = var_temp_new( TYPE_DOUBLE, 1, model->NofNodes*cols );
    rf = MATR(res);

    for( k=0; k<cols; k++ )
    {
        for( j=0; j<model->NofNodes; j++ )
        {
            rf[k*model->NofNodes+j] = 0.0;
            References[j] = 0;
        }
 
        for( i=0; i<model->NofElements; i++ )
        {
            elm  = &model->Elements[i];
            elmt = elm->ElementType;

            for( j=0; j<elmt->NumberOfNodes; j++ )
            {
                n = k*model->NofNodes + elm->Topology[j];
                FX[j] = ax[n];
                FY[j] = ay[n];
                FZ[j] = az[n];
            }

             if ( dim==3 && elmt->PartialW )
                 status = elm_element_divergence_3D( model,elm,FX,FY,FZ,D );
             else if ( dim==2 && elmt->PartialV )
                 status = elm_element_divergence_2D( model,elm,FX,FY,FZ,D );
             else continue;

            if ( !status ) 
            {
                fprintf( stderr, "ELEMENT [%d]: Jacobian is singular.\n",i );
                continue;
            }

            for( j=0; j<elmt->NumberOfNodes; j++ )
            {
                n =  elm->Topology[j];
                References[n]++;
                rf[k*model->NofNodes + n] += D[j];
            }
        }

        for( j=0; j<model->NofNodes; j++ )
        {
            if ( References[j] != 0 )
            {
                n = k*model->NofNodes + j;
                rf[n] /= References[j];
            } 
        } 
    }

    free( References );
    return res;
}


/*******************************************************************************
 *
 *     Name:        elm_element_ddx_2D
 *
 *     Purpose:     Compute elementwise matrix of second partial derivates
 *                  for surface elements at given point u,v.
 *
 *     Parameters:
 *
 *         Input:   (element_t *) structure describing the element
 *                  (double *) F nodel values of the quantity
 *                  (double *) u,v point at which to evaluate
 *
 *         Output:   3x3 matrix of partial derivates
 *
 *   Return value:  FALSE if element is degenerated, TRUE otherwise.
 *
 ******************************************************************************/
int elm_element_ddx_2D(
                         element_model_t *model, element_t *elm, double *F,
                                 double u,double v,double *Values
                      )
{
    static double X[ELM_MAX_ELEMENT_NODES];
    static double Y[ELM_MAX_ELEMENT_NODES];
    static double Z[ELM_MAX_ELEMENT_NODES];

    static double partials[4];

    int i,n;

    double ddfddx,ddfddy,ddfddz,ddfdxdy,ddfdxdz,ddfdydz;
    double ddfddu,ddfddv,ddfdudv,dfdu,dfdv, Cddfddu,Cddfddv,Cddfdudv;

    double ddxddu,ddxddv,ddxdudv,dxdu,dxdv;
    double ddyddu,ddyddv,ddydudv,dydu,dydv;
    double ddzddu,ddzddv,ddzdudv,dzdu,dzdv;

    double a,G_uu,G_vv,G_uv,H_uu,H_uv,H_vu,H_vv, detG;
    double C1_uu_u,C1_uu_v,C1_uv_u,C1_uv_v,C1_vv_u,C1_vv_v;
    double C2_uu_u,C2_uu_v,C2_uv_u,C2_uv_v,C2_vv_u,C2_vv_v;

    element_type_t *elmt = elm->ElementType;

    double (*dNdU)() = elmt->PartialU;
    double (*dNdV)() = elmt->PartialV;

    double (*ddN)()  = elmt->SecondPartials;

    for( i=0; i<elmt->NumberOfNodes; i++ )
    {
        n    = elm->Topology[i];
        X[i] = model->NodeArray[n];
        Y[i] = model->NodeArray[model->NofNodes+n];
        Z[i] = model->NodeArray[2*model->NofNodes+n];
    }

    dxdu = (*dNdU)( X,u,v );
    dydu = (*dNdU)( Y,u,v );
    dzdu = (*dNdU)( Z,u,v );

    dxdv = (*dNdV)( X,u,v );
    dydv = (*dNdV)( Y,u,v );
    dzdv = (*dNdV)( Z,u,v );

    (*ddN)( X,u,v,partials );
    ddxddu  = partials[0];
    ddxdudv = partials[1];
    ddxddv  = partials[3];

    (*ddN)( Y,u,v,partials );
    ddyddu  = partials[0];
    ddydudv = partials[1];
    ddyddv  = partials[3];

    (*ddN)( Z,u,v,partials );
    ddzddu  = partials[0];
    ddzdudv = partials[1];
    ddzddv  = partials[3];

    /*
     * covariant metric tensor
     *
     * TODO: multiply by global metric if not identity
     */
    G_uu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    G_vv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    G_uv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    detG = G_uu*G_vv - G_uv*G_uv;
    if ( ABS(detG)<AEPS ) return FALSE;

    /*
     * contravariant metric tensor
     */
    H_uu =  G_vv / detG;
    H_uv = -G_uv / detG;
    H_vu = -G_uv / detG;
    H_vv =  G_uu / detG;

    /*
     * christoffel symbols of the second kind
     *
     * TODO: this is for cartesian global coordinates only...
     */
    C2_uu_u = ddxddu  * dxdu  + ddyddu  * dydu + ddzddu  * dzdu;
    C2_uu_v = ddxddu  * dxdv  + ddyddu  * dydv + ddzddu  * dzdv;

    C2_uv_u = ddxdudv * dxdu  + ddydudv * dydu + ddzdudv * dzdu;
    C2_uv_v = ddxdudv * dxdv  + ddydudv * dydv + ddzdudv * dzdv;

    C2_vv_u = ddxddv  * dxdu  + ddyddv  * dydu + ddzddv  * dzdu;
    C2_vv_v = ddxddv  * dxdv  + ddyddv  * dydv + ddzddv  * dzdv;

    /*
     * christoffel symbols of the first kind
     */
    C1_uu_u = H_uu * C2_uu_u + H_uv * C2_uu_v;
    C1_uu_v = H_vu * C2_uu_u + H_vv * C2_uu_v;

    C1_uv_u = H_uu * C2_uv_u + H_uv * C2_uv_v;
    C1_uv_v = H_vu * C2_uv_u + H_vv * C2_uv_v;

    C1_vv_u = H_uu * C2_vv_u + H_uv * C2_vv_v;
    C1_vv_v = H_vu * C2_vv_u + H_vv * C2_vv_v;

    /*
     * now we are ready to compute second covariant derivates of f
     */
    dfdu = (*dNdU)( F,u,v );
    dfdv = (*dNdU)( F,u,v );

    (*ddN)( F,u,v,partials );
    ddfddu  = partials[0];
    ddfdudv = partials[1];
    ddfddv  = partials[3];

    ddfddu  -= C1_uu_u * dfdu + C1_uu_v * dfdv;
    ddfddv  -= C1_vv_u * dfdu + C1_vv_v * dfdv;
    ddfdudv -= C1_uv_u * dfdu + C1_uv_v * dfdv;

    /*
     * convert to contravariant form
     */
    Cddfddu  = H_uu * (H_uu * ddfddu + H_uv * ddfdudv) + H_uv * (H_uu * ddfdudv + H_uv * ddfddv);
    Cddfddv  = H_vu * (H_vu * ddfddu + H_vv * ddfdudv) + H_vv * (H_vu * ddfdudv + H_vv * ddfddv);
    Cddfdudv = H_uu * (H_vu * ddfddu + H_vv * ddfdudv) + H_uv * (H_vu * ddfdudv + H_vv * ddfddv);

    /*
     * transform to global coordinates 
     *
     * TODO: lower indexes if not cartesian x,y,z
     */
    ddfddx  = dxdu * dxdu * Cddfddu + dxdv * dxdv * Cddfddv + 2 * dxdu * dxdv * Cddfdudv;
    ddfddy  = dydu * dydu * Cddfddu + dydv * dydv * Cddfddv + 2 * dydu * dydv * Cddfdudv;
    ddfddz  = dzdu * dzdu * Cddfddu + dzdv * dzdv * Cddfddv + 2 * dzdu * dzdv * Cddfdudv;

    ddfdxdy = dxdu * dydu * Cddfddu + (dxdu * dydv + dxdv * dydu) * Cddfdudv + dxdv * dydv * Cddfddv;
    ddfdxdz = dxdu * dzdu * Cddfddu + (dxdu * dzdv + dxdv * dzdu) * Cddfdudv + dxdv * dzdv * Cddfddv;
    ddfdydz = dydu * dzdu * Cddfddu + (dydu * dzdv + dydv * dzdu) * Cddfdudv + dydv * dzdv * Cddfddv;

    Values[0] =  ddfddx;
    Values[1] = ddfdxdy;
    Values[2] = ddfdxdz;
    Values[3] = ddfdxdy;
    Values[4] =  ddfddy;
    Values[5] = ddfdydz;
    Values[6] = ddfdxdz;
    Values[7] = ddfdydz;
    Values[8] =  ddfddz;

    return TRUE;
}
