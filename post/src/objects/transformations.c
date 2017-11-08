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
 * Object transformations
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
 *                       Date: 1 Oct 1995
 *
 ******************************************************************************/

#include "../elmerpost.h"

/*******************************************************************************
 *
 *     Name:          obj_init_matrix( matrix_t )
 *
 *     Purpose:       Set given transformation matrix to unit matrix.
 *                    Internal only.
 *
 *     Parameters:
 *
 *         Input:     none
 *
 *         Output:    (matrix_t) matrix to set
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_init_matrix( matrix_t matrix )
{
    memset( matrix,0, sizeof(matrix_t) );

    matrix[0][0] = matrix[1][1] =
    matrix[2][2] = matrix[3][3] = 1.0;
}

/*******************************************************************************
 *
 *     Name:          obj_init_transform( transform_t * )
 *
 *     Purpose:       Initialize given transform_t structure.
 *                    Internal only.
 *
 *     Parameters:
 *
 *         Input:     none
 *
 *         Output:    (transform_t *)
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_init_transform( transform_t *transform )
{
    transform->RotX = transform->RotY = transform->RotZ = 0.0;
    transform->TrnX = transform->TrnY = transform->TrnZ = 0.0;
    transform->SclX = transform->SclY = transform->SclZ = 1.0;

    transform->TransformPriority = trn_pri_str;
    transform->RotationPriority  = rot_pri_xyz;

    obj_init_matrix( transform->Matrix );

    obj_init_matrix( transform->RotMatrix );
    obj_init_matrix( transform->SclMatrix );
    obj_init_matrix( transform->TrnMatrix );
}

/*******************************************************************************
 *
 *     Name:          obj_get_matrix( matrix_t, object_t * )
 *
 *     Purpose:       Return objects transformation matrix given
 *                    (object_t *) structure.
 *
 *     Parameters:
 *
 *         Input:     (object_t *)
 *
 *         Output:    (matrix_t)
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_get_matrix(matrix_t matrix,object_t *object)
{
    transform_t *transform = &object->Transform;

    memcpy( matrix,transform->Matrix,sizeof(matrix_t) );
}

/*******************************************************************************
 *
 *     Name:          obj_get_matrix_transpose( matrix_t, object_t * )
 *
 *     Purpose:       Return transpose of objects transformation matrix given
 *                    (object_t *) structure.
 *
 *     Parameters:
 *
 *         Input:     (object_t *)
 *
 *         Output:    (matrix_t)
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_get_matrix_transpose( matrix_t matrix,object_t *object )
{
    transform_t *transform = &object->Transform;
    int i,j;

    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ ) matrix[j][i] = transform->Matrix[i][j];
}

/*******************************************************************************
 *
 *     Name:          obj_copy_matrix( matrix_t A, matrix_t B)
 *
 *     Purpose:       Copy matrix A to matrix B
 *
 *     Parameters:
 *
 *         Input:     (matrix_t)
 *
 *         Output:    (matrix_t)
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_copy_matrix( matrix_t A,matrix_t B )
{
    memcpy( A,B,sizeof(matrix_t) );
}

/*******************************************************************************
 *
 *     Name:          obj_mult_matrix( matrix_t A, matrix_t B)
 *                    obj_mult_matrix_left( matrix_t A, matrix_t B)
 *
 *     Purpose:       return A = A*B. Internal only.
 *
 *     Parameters:
 *
 *         Input:     (matrix_t)A, (matrix_t)B
 *
 *         Output:    (matrix_t)A
 *
 *   Return value:    void
 *
 ******************************************************************************/
#define obj_mult_matrix( A, B )  obj_mult_matrix_left( A,B )
static void obj_mult_matrix_left( matrix_t A,matrix_t B ) 
{
    matrix_t R;
    double s;
    int i,j,k;
  
    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ )
    {
        s = 0.0;
        for( k=0; k<4; k++ ) s += A[i][k]*B[k][j];
        R[i][j] = s;
    }

    obj_copy_matrix( A,R );
}

/*******************************************************************************
 *
 *     Name:          obj_mult_matrix_right( matrix_t A, matrix_t B)
 *
 *     Purpose:       return A = B*A. Internal only.
 *
 *     Parameters:
 *
 *         Input:     (matrix_t)A, (matrix_t)B
 *
 *         Output:    (matrix_t)A
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_mult_matrix_right( matrix_t A,matrix_t B ) 
{
    matrix_t R;
    double s;
    int i,j,k;
  
    for( i=0; i<4; i++ )
    for( j=0; j<4; j++ )
    {
        s = 0.0;
        for( k=0; k<4; k++ ) s += B[i][k]*A[k][j];
        R[i][j] = s;
    }

    obj_copy_matrix( A,R );
}

/*******************************************************************************
 *
 *     Name:          obj_translate_matrix( matrix_t,double,double,double )
 *
 *     Purpose:       return matrix_t corresponding to translations x,y,z
 *                    Internal only.
 *
 *     Parameters:
 *
 *         Input:    (double,double,double) translations x,y,and z
 *
 *         Output:   (matrix_t) resulting transformation matrix
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_translate_matrix( matrix_t M, double x, double y, double z )
{
     int i;

     obj_init_matrix( M );

     M[0][3] = x;
     M[1][3] = y;
     M[2][3] = z;
/*
 *   for( i=0; i<4; i++ )
 *   {
 *       M[0][i] += x*M[3][i];
 *       M[1][i] += y*M[3][i];
 *       M[2][i] += z*M[3][i];
 *   }
 */
}

/*******************************************************************************
 *
 *     Name:          obj_scale_matrix( matrix_t,double,double,double )
 *
 *     Purpose:       return matrix_t corresponding to scalings x,y,z
 *                    Internal only.
 *
 *     Parameters:
 *
 *         Input:    (double,double,double) scaleings x,y,and z
 *
 *         Output:   (matrix_t) resulting transformation matrix
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_scale_matrix( matrix_t M,double x, double y, double z )
{
     int i;

     obj_init_matrix( M );

     M[0][0] = x;
     M[1][1] = y;
     M[2][2] = z;
/*
 *   for( i=0; i<4; i++ )
 *   {
 *       M[0][i] *= x;
 *       M[1][i] *= y;
 *       M[2][i] *= z;
 *   }
 */
}

/*******************************************************************************
 *
 *     Name:          obj_internal_rotate_matrix( matrix_t,double,double,int )
 *
 *     Purpose:       return matrix_t corresponding to rotation about axis
 *                    given, by amount given in sin(a), cos(a).
 *                    Internal only.
 *
 *     Parameters:
 *
 *         Input:    (double) sine of the angle to rotate
 *                   (double) cosine of the angle to rotate
 *                   (int)   axis about which to rotate:
 *                       x: axis = 0, y: axis = 1, z: axis = 2
 *
 *         Output:   (matrix_t) resulting transformation matrix
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_internal_rotate_matrix( matrix_t M,double s, double c, int axis )
{
     double t;
     int i;

     switch( axis )
     {
         case 0: for( i=0; i<4; i++ )
                 {
                     t = M[1][i];
                     M[1][i] = c*t - s*M[2][i];
                     M[2][i] = c*M[2][i] + s*t;
                 }
         break;
 
         case 1: for( i=0; i<4; i++ )
                 {
                     t = M[0][i];
                     M[0][i] = c*t + s*M[2][i];
                     M[2][i] = c*M[2][i] - s*t;
                 }
         break;

         case 2: for( i=0; i<4; i++ )
                 {
                     t = M[0][i];
                     M[0][i] = c*t - s*M[1][i];
                     M[1][i] = c*M[1][i] + s*t;
                 }
         break;
     }
}

/*******************************************************************************
 *
 *     Name:          obj_set_rotation_priority( object_t *,rot_pri_t )
 *
 *     Purpose:       set rotation priority for an object given
 *
 *     Parameters:
 *
 *         Input:    (rot_pri_t) rotation priority, one of:
 *
 *                    rot_pri_xyz
 *                    rot_pri_xzy
 *                    rot_pri_yxz
 *                    rot_pri_yzx
 *                    rot_pri_zxy
 *                    rot_pri_zyx
 *
 *         Output:   (object_t *) is modified
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_set_rotation_priority( object_t *object, rot_pri_t priority)
{
    object->Transform.RotationPriority = priority;
}

/*******************************************************************************
 *
 *     Name:          obj_set_transform_priority( object_t *, trn_pri_t )
 *
 *     Purpose:       set transformation priority for an object given
 *
 *     Parameters:
 *
 *         Input:    (trn_pri_t) transformation priority, one of:
 *
 *                    trn_pri_trs
 *                    trn_pri_tsr
 *                    trn_pri_rts
 *                    trn_pri_rst
 *                    trn_pri_str
 *                    trn_pri_srt
 *
 *         Output:   (object_t *) is modified
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_set_transform_priority( object_t *object, trn_pri_t priority)
{
    object->Transform.TransformPriority = priority;
}

/*******************************************************************************
 *
 *     Name:          obj_rotate_matrix( matrix_t,rot_pri_t,double,double,double )
 *
 *     Purpose:       return matrix_t corresponding to rotation about
 *                    three axes by angles given. Internal only.
 *
 *     Parameters:
 *
 *         Input:    (rot_pri_t) rotation priority, one of:
 *
 *                    rot_pri_xyz
 *                    rot_pri_xzy
 *                    rot_pri_yxz
 *                    rot_pri_yzx
 *                    rot_pri_zxy
 *                    rot_pri_zyx
 *
 *                    (double)  angle about 'x' - axis
 *                    (double)  angle about 'y' - axis
 *                    (double)  angle about 'z' - axis
 *
 *         Output:   (matrix_t) resulting transformation matrix
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_rotate_matrix( matrix_t M,rot_pri_t rot_pri,double x,double y,double z )
{
     int i,j,axis,a_ord[3];
     double a;
     matrix_t N,P;

     obj_init_matrix( M );

     switch( rot_pri )
     {
         case rot_pri_xyz: a_ord[0]=0; a_ord[1]=1; a_ord[2]=2; break;
         case rot_pri_xzy: a_ord[0]=0; a_ord[1]=2; a_ord[2]=1; break;
         case rot_pri_yxz: a_ord[0]=1; a_ord[1]=0; a_ord[2]=2; break;
         case rot_pri_yzx: a_ord[0]=1; a_ord[1]=2; a_ord[2]=0; break;
         case rot_pri_zxy: a_ord[0]=2; a_ord[1]=0; a_ord[2]=1; break;
         case rot_pri_zyx: a_ord[0]=2; a_ord[1]=1; a_ord[2]=0; break;
     }

     for( i=0; i<3; i++ )
     {
         axis = a_ord[i];
         switch(axis)
         {
             case 0: a = x; break;
             case 1: a = y; break;
             case 2: a = z; break;
         }
         a *= PiDiv180;

         obj_init_matrix( N );
         obj_internal_rotate_matrix( N,sin(a),cos(a),axis );

         obj_mult_matrix_right( M,N );
     }
}

/*******************************************************************************
 *
 *     Name:          obj_rotate_mult_matrix( transform_t *,double *,double *,double *)
 *
 *     Purpose:       return matrix_t corresponding to rotation  about
 *                    three axes by angles given multiplied by objects
 *                    transformation matrix. Internal only.
 *
 *     Parameters:
 *
 *         Input:    (transform_t *)
 *
 *                    (double *) angle about 'x' - axis
 *                    (double *) angle about 'y' - axis
 *                    (double *) angle about 'z' - axis
 *
 *         Output:   (transform_t *)->RotMatrix is modified as are the
 *                   angles.
 *
 *   Return value:    void
 *
 ******************************************************************************/
#define ROUND_PI(A) (180*(int)(((A)>0?(A)+90.0:(A)-90.0)/180.0))

static void obj_rotate_mult_matrix(transform_t *transform,double *x,double *y,double *z )
{
   matrix_t M;

   int i,j;
   double a;

    for( i=0; i<3; i++ )
    {
        switch(i)
        {
            case 0: a = *x; break;
            case 1: a = *y; break;
            case 2: a = *z; break;
        }
        a *= PiDiv180;

        obj_init_matrix( M );
        obj_internal_rotate_matrix( M,sin(a),cos(a),i );

        if ( transform->RotationPriority == rot_pri_local )
        {
            obj_mult_matrix_left( transform->RotMatrix,M );
        } else if ( transform->RotationPriority == rot_pri_parent )
        {
            obj_mult_matrix_right( transform->RotMatrix,M );
        }
    }

    obj_copy_matrix( M,transform->RotMatrix );

    *x =  atan2( M[2][1],M[2][2] ) / PiDiv180;
    *y = -asin(  M[2][0] )/ PiDiv180;
    *z =  atan2( M[1][0],M[0][0] ) / PiDiv180;

    *x += ROUND_PI( transform->RotX - *x );
    *y += ROUND_PI( transform->RotY - *y );
    *z += ROUND_PI( transform->RotZ - *z );
}

/*******************************************************************************
 *
 *     Name:          obj_compute_matrix( transform_t * )
 *
 *     Purpose:       Return total transformation matrix_t given rotations,
 *                    scaleing, and translations. Internal only.
 *
 *     Parameters:
 *
 *         Input:    (transform_t *)
 *
 *         Output:   (transform_t *)->Matrix is modified
 *
 *   Return value:    void
 *
 ******************************************************************************/
static void obj_compute_matrix( transform_t *transform )
{
    transform_list_t *child;

    matrix_t M;

    obj_init_matrix( M );

    switch( transform->TransformPriority )
    {
       case trn_pri_trs: obj_mult_matrix( M,transform->TrnMatrix );
                         obj_mult_matrix( M,transform->RotMatrix );
                         obj_mult_matrix( M,transform->SclMatrix );
       break;

       case trn_pri_tsr: obj_mult_matrix( M,transform->TrnMatrix );
                         obj_mult_matrix( M,transform->SclMatrix );
                         obj_mult_matrix( M,transform->RotMatrix );
       break;

       case trn_pri_rts: obj_mult_matrix( M,transform->RotMatrix );
                         obj_mult_matrix( M,transform->TrnMatrix );
                         obj_mult_matrix( M,transform->SclMatrix );
       break;

       case trn_pri_rst: obj_mult_matrix( M,transform->RotMatrix );
                         obj_mult_matrix( M,transform->SclMatrix );
                         obj_mult_matrix( M,transform->TrnMatrix );
       break;

       case trn_pri_str: obj_mult_matrix( M,transform->SclMatrix );
                         obj_mult_matrix( M,transform->TrnMatrix );
                         obj_mult_matrix( M,transform->RotMatrix );
       break;

       case trn_pri_srt: obj_mult_matrix( M,transform->SclMatrix );
                         obj_mult_matrix( M,transform->RotMatrix );
                         obj_mult_matrix( M,transform->TrnMatrix );
       break;
    }

    if ( transform->Parent && transform->Parent != transform )
    {
        memcpy( transform->Matrix, transform->Parent->Matrix, sizeof(matrix_t) );
        obj_mult_matrix( transform->Matrix,M ); 
    } else
    {
        obj_copy_matrix( transform->Matrix,M );
    }

    for( child = transform->Children; child != NULL; child = child->Next )
    {
        obj_compute_matrix( child->Entry );
    }
}


/*******************************************************************************
 *
 *     Name:          obj_set_parent( object_t *obj, object_t *parent )
 *
 *     Purpose:       Set transformation parent of an object 
 *
 *     Parameters:
 *
 *         Input:    both structures are read
 *
 *         Output:   both structures are modified
 *
 *   Return value:    malloc() success
 *
 ******************************************************************************/
int obj_set_parent( object_t *object,object_t *parent )
{
    transform_t *objtrans = &object->Transform;
    transform_t *partrans = &parent->Transform;

    transform_list_t *child;

    if ( object == parent ) return TRUE;

    objtrans->Parent = partrans;
    
    if ( !(child = (transform_list_t *)calloc(1,sizeof(transform_list_t)) ) )
    {
        fprintf( stderr, "obj_set_parent: FATAL: can't allocate a few bytes of memory.\n" );
        return FALSE;
    }

    child->Entry = objtrans;
    child->Next  = partrans->Children;

    partrans->Children = child;
    
    obj_compute_matrix( objtrans );

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:          obj_rotate( object_t *,double,double,double,int,int )
 *
 *     Purpose:       Modifify object transformation matrix given rotations
 *
 *     Parameters:
 *
 *         Input:    (object_t *)
 *                   (double,double,double) input rotations
 *                   (int) flag if rotations should be relative (TRUE)
 *                          or absolute (FALSE)
 *
 *         Output:   (object->Transform_t *)->XXXMatrix are modified
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_rotate( object_t *object,double x,double y,double z,int which,int relative )
{
    double rx,ry,rz;

    transform_t *transform = &object->Transform;

    switch( which )
    {
        case 'x': if ( transform->RotationPriority < rot_pri_local )
                  {
                      if ( relative ) x += transform->RotX;
                      y  = transform->RotY;
                      z  = transform->RotZ;
                  } else
                  { 
                      if ( !relative )
                      {
                          x -= transform->RotX;
                      }
                      y = z = 0.0;
                  }
        break;

        case 'y': if ( transform->RotationPriority < rot_pri_local )
                  {
                      if ( relative ) y = x + transform->RotY; else y = x;
                      x  = transform->RotX;
                      z  = transform->RotZ;
                  } else
                  { 
                      if ( relative )
                      {
                          y = x;
                      } else
                      {
                          y = x - transform->RotY;
                      }
                      x = z = 0.0;
                  }
        break;

        case 'z': if ( transform->RotationPriority < rot_pri_local )
                  {
                      if ( relative ) z = x + transform->RotZ; else z = x;
                      x  = transform->RotX;
                      y  = transform->RotY;
                  } else
                  { 
                      if ( relative )
                      {
                          z = x;
                      } else
                      {
                          z = x - transform->RotZ;
                      }
                      x = y = 0.0;
                  }
        break;

        case 'a': if ( transform->RotationPriority < rot_pri_local )
                  {
                      if ( relative )
                      {
                          x += transform->RotX;
                          y += transform->RotY;
                          z += transform->RotZ;
                      }
                  } else
                  {
                      if ( !relative ) obj_init_matrix( transform->RotMatrix );
                  }
        break;
    }

    if ( transform->RotationPriority >= rot_pri_local )
    {
        obj_rotate_mult_matrix( transform,&x,&y,&z );
    } else
    {
        obj_rotate_matrix( transform->RotMatrix,transform->RotationPriority,x,y,z );
    }

    transform->RotX = x;
    transform->RotY = y;
    transform->RotZ = z;

    obj_compute_matrix( transform );
}

/*******************************************************************************
 *
 *     Name:          obj_scale( object_t *,double,double,double,int,int )
 *
 *     Purpose:       Modifify object transformation matrix given scaleings
 *
 *     Parameters:
 *
 *         Input:    (object_t *)
 *                   (double,double,double) input scaleings
 *                   (int) flag if scaleings should be relative (TRUE)
 *                          or absolute (FALSE)
 *
 *         Output:   (object->Transform_t *)->XXXMatrix are modified
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_scale( object_t *object,double x,double y,double z,int which, int relative )
{
    transform_t *transform = &object->Transform;
    double s;

    if ( relative )
    {
        switch( which )
        {
            case 's': s = 1+x;
                      x = s*transform->SclX;
                      y = s*transform->SclY;
                      z = s*transform->SclZ;
            break;

            case 'x': x = (1+x)*transform->SclX; break;
            case 'y': y = (1+x)*transform->SclY; break;
            case 'z': z = (1+x)*transform->SclZ; break;

            case 'a': x = (1+x)*transform->SclX;
                      y = (1+y)*transform->SclY;
                      z = (1+z)*transform->SclZ;
            break;
        }
    }

    if ( x > 0 )
    {
        x = MAX(x, 1.0e-6);
    } else {
        x = MIN(x,-1.0e-6);
    }

    if ( y > 0 )
    {
        y = MAX(y, 1.0e-6);
    } else {
        y = MIN(y,-1.0e-6);
    }

    if ( z > 0 )
    {
        z = MAX(z, 1.0e-6);
    } else {
        z = MIN(z,-1.0e-6);
    }

    switch( which )
    {
        case 'a': transform->SclX = x;
                  transform->SclY = y;
                  transform->SclZ = z;
        break;

        case 'x': transform->SclX = x; break;
        case 'y': transform->SclY = y; break;
        case 'z': transform->SclZ = z; break;
    }

    obj_scale_matrix(
                                      transform->SclMatrix,
                    transform->SclX, transform->SclY, transform->SclZ
                );

    obj_compute_matrix( transform );
}

/*******************************************************************************
 *
 *     Name:          obj_translate( object_t *,double,double,double,int,int  )
 *
 *     Purpose:       Modifify object transformation matrix given translations
 *
 *     Parameters:
 *
 *         Input:    (object_t *)
 *                   (double,double,double) input translations
 *                   (int) flag if translations should be relative (TRUE)
 *                          or absolute (FALSE)
 *
 *         Output:   (object->Transform_t *)->XXXMatrix are modified
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_translate( object_t *object,double x,double y,double z,int which,int relative )
{
    transform_t *transform = &object->Transform;

    switch( which )
    {
        case 'x': transform->TrnX = (relative?transform->TrnX:0.0)+x; break;
        case 'y': transform->TrnY = (relative?transform->TrnY:0.0)+x; break;
        case 'z': transform->TrnZ = (relative?transform->TrnZ:0.0)+x; break;
        case 'a': if ( relative )
                  {
                      transform->TrnX += x;
                      transform->TrnY += y;
                      transform->TrnZ += z;
                  } else
                  {
                      transform->TrnX  = x;
                      transform->TrnY  = y;
                      transform->TrnZ  = z;
                  }
        break;
    }

    obj_translate_matrix(
                                transform->TrnMatrix,
                    transform->TrnX, transform->TrnY, transform->TrnZ
                );

    obj_compute_matrix( transform );
}

/*******************************************************************************
 *
 *     Name:          obj_set_matrix( object_t * )
 *
 *     Purpose:       Tell graphics module about transformation matrix of
 *                    a given object
 *
 *     Parameters:
 *
 *         Input:     (object_t *)
 *
 *         Output:    none
 *
 *   Return value:    void
 *
 ******************************************************************************/
void obj_set_matrix( object_t *object )
{
    transform_t *transform = &object->Transform;
    matrix_t matrix;

    obj_get_matrix_transpose( matrix,object );

    gra_mult_matrix( matrix );
}
