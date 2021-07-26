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
 *     The character of the routines in this file.
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
 *                       Date: 1 May 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

#include "../elmerpost.h"


typedef struct particle_s
{
    scalar_t *VectorData[3];
    scalar_t *ColorData;

    int NofParticles;
    scalar_t *ParticleData[5];

    int SaveNofParticles;

    double OutDT,MaxDT,Tolerance;

    material_t *Material;
    colormap_t *ColorMap;

    arrow_style_t ArrowStyle;

    double *X,*Y,*Z,*C,*I;

    int Advance;
    particle_integ_policy_t IntegPolicy;
    particle_integ_method_t IntegMethod;

    particle_style_t Style;

    int LineQuality;
    double LineWidth;
    line_style_t LineStyle;
} particle_t;

extern int CurrentTimeStep;

static double vis_particles_interpolate
  (
    particle_t *particles, element_model_t *model, int n, double *velo, double t, double x, double y,double z
  )
{
    static double NX[ELM_MAX_ELEMENT_NODES];
    static double NY[ELM_MAX_ELEMENT_NODES];
    static double NZ[ELM_MAX_ELEMENT_NODES];
    static double NF[ELM_MAX_ELEMENT_NODES];

    double u,v,w;

    double *XC = &model->NodeArray[0];
    double *YC = &model->NodeArray[model->NofNodes];
    double *ZC = &model->NodeArray[model->NofNodes*2];

    int i,j,k,k0,k1,*T;

    element_type_t *elmt;
    element_t *elements = model->Elements;

    if ( particles->ParticleData[4]->f[n] < 0 ) return 0.0;

    k0 = particles->ParticleData[4]->f[n];
    k1 = particles->ParticleData[4]->f[n] - 1;

    t /= particles->OutDT;

    for( j=0; j<model->NofElements; j++ )
    {
        if ( j & 1 )
        {
            if ( k0 >= model->NofElements ) i = k1--; else i = k0++;
        } else {
            if ( k1  < 0  ) i = k0++; else i = k1--;
        }

        T = elements[i].Topology;
        elmt = elements[i].ElementType;

        if ( !elmt->PointInside )
        {
            fprintf(
               stderr, "WARINING: Particles not implemented for element type: [%s]. Element ingnored.\n",
                            elmt->ElementName
                );
            continue;
        }

        for( k=0; k<elmt->NumberOfNodes; k++ )
        {
            NX[k] = XC[T[k]];
            NY[k] = YC[T[k]];
            NZ[k] = ZC[T[k]];
        }

        if ( !(*elmt->PointInside)( NX,NY,NZ,x,y,z,&u,&v,&w ) ) continue;

        if ( CurrentTimeStep == model->NofTimesteps-1 )
        {
            for( k=0; k<elmt->NumberOfNodes; k++ )
            {
                NF[k] = velo[CurrentTimeStep*model->NofNodes + T[k]];
            }
        }
        else
        {
            for( k=0; k<elmt->NumberOfNodes; k++ )
            {
                NF[k]  = (1-t)*velo[CurrentTimeStep*model->NofNodes + T[k]];
                NF[k] += t*velo[(CurrentTimeStep+1)*model->NofNodes + T[k]];
            }
        }

        particles->ParticleData[4]->f[n] = i;
        return (*elmt->FunctionValue)( NF,u,v,w );
    }

    particles->ParticleData[4]->f[n] = -1;
    return 0.0;
}

static void vis_RungeKutta(
              particle_t *particles,
              element_model_t *model,
              int n,
              double *x,
              double *y,
              double *z,
              double *c,
              double t,
              double dt,
              double *f_x0,
              double *f_y0,
              double *f_z0,
              double *f_c0 )
{
    double k1,k2,k3,k4,dx=*x,dy=*y,dz=*z;

    k1 = k2 = k3 = k4 = 0;

    k1 = dt * vis_particles_interpolate( particles, model, n, f_x0, t, dx, dy, dz );
    k2 = dt * vis_particles_interpolate( particles, model, n, f_x0, t + dt / 2.0, dx + k1 / 2.0, dy, dz );
    k3 = dt * vis_particles_interpolate( particles, model, n, f_x0, t + dt / 2.0, dx + k2 / 2.0, dy, dz );
    k4 = dt * vis_particles_interpolate( particles, model, n, f_x0, t + dt, dx + k3, dy, dz );
    *x += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;

    k1 = dt * vis_particles_interpolate( particles, model, n, f_y0, t, dx, dy, dz );
    k2 = dt * vis_particles_interpolate( particles, model, n, f_y0, t + dt / 2.0, dx, dy + k1 / 2.0, dz );
    k3 = dt * vis_particles_interpolate( particles, model, n, f_y0, t + dt / 2.0, dx, dy + k2 / 2.0, dz );
    k4 = dt * vis_particles_interpolate( particles, model, n, f_y0, t + dt, dx, dy + k3, dz );
    *y += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;

    k1 = dt * vis_particles_interpolate( particles, model, n, f_z0, t, dx, dy, dz );
    k2 = dt * vis_particles_interpolate( particles, model, n, f_z0, t + dt / 2.0, dx, dy, dz + k1 / 2.0 );
    k3 = dt * vis_particles_interpolate( particles, model, n, f_z0, t + dt / 2.0, dx, dy, dz + k2 / 2.0 );
    k4 = dt * vis_particles_interpolate( particles, model, n, f_z0, t + dt, dx, dy, dz + k3 );

    *z += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;

    if ( f_c0 )
    {
        *c = vis_particles_interpolate( particles, model, n, f_c0, t, *x, *y, *z );
    }
}

static void vis_Euler(
              particle_t *particles,
              element_model_t *model,
              int n,
              double *x,
              double *y,
              double *z,
              double *c,
              double t,
              double dt,
              double *f_x0,
              double *f_y0,
              double *f_z0,
              double *f_c0 )
{
    double dx=*x,dy=*y,dz=*z;

    *x += dt * vis_particles_interpolate( particles, model, n, f_x0, t, dx, dy, dz );
    *y += dt * vis_particles_interpolate( particles, model, n, f_y0, t, dx, dy, dz );
    *z += dt * vis_particles_interpolate( particles, model, n, f_z0, t, dx, dy, dz );

    if ( f_c0 )
    {
        *c = vis_particles_interpolate( particles, model, n, f_c0, t, *x, *y, *z );
    }
}

extern double xmin,xmax,ymin,ymax,zmin,zmax;

int vis_particle( geometry_t *geometry, element_model_t *model, particle_t *particles, double gtime )
{
    double *f_x0 = NULL, *f_y0 = NULL, *f_z0 = NULL, *f_c0 = NULL;
    double *PartX, *PartY, *PartZ, *PartC, *PartI;

    double x, y, z, c, CAdd = 0.0, CScl = 1.0;

    double LX,LY,LZ, LC;
    double r = 0.05, *C;

    double x0,y0,x1,y1,z0,z1;
    double k1,k2,k3,k4;

    double t, dt, ddt;

    double CEPS = particles->Tolerance;

    int i, j, k, n, quick;

    double s = MAX(MAX(xmax-xmin,ymax-ymin),zmax-zmin);

    float pnt[3], vec[3];

    if ( !GlobalOptions.StereoMode )
      if ( particles->Material->Diffuse[3]  < 1.0 )
      {
          if ( GlobalPass != 0 ) return TRUE;
      } else if ( GlobalPass == 0 )
      {
          return TRUE;
      }

    if ( particles->NofParticles < 1 ) return TRUE;
    if ( !particles->VectorData   || !particles->VectorData[0]->f )   return TRUE;
    if ( !particles->ParticleData || !particles->ParticleData[0]->f ) return TRUE;

    f_x0 = particles->VectorData[0]->f;
    f_y0 = particles->VectorData[1]->f;
    f_z0 = particles->VectorData[2]->f;

    gra_set_material( particles->Material );

    if ( particles->ColorData && particles->ColorData->f )
    {
        f_c0 = particles->ColorData->f;

        CAdd = particles->ColorData->min;
        if ( particles->ColorData->max - particles->ColorData->min != 0.0 )
            CScl = 1.0 / ( particles->ColorData->max - particles->ColorData->min );
        else
            CScl = 1.0;

        gra_set_colormap( particles->ColorMap );
    } else gra_set_colormap( NULL );

    PartX = particles->ParticleData[0]->f;
    PartY = particles->ParticleData[1]->f;
    PartZ = particles->ParticleData[2]->f;
    PartC = particles->ParticleData[3]->f;
    PartI = particles->ParticleData[4]->f;

    n = particles->NofParticles;
    if ( particles->SaveNofParticles != n )
    {
        {
          VARIABLE *var = var_new( "_particle_last", TYPE_DOUBLE,5,n );

          particles->X = &M(var,0,0);
          particles->Y = &M(var,1,0);
          particles->Z = &M(var,2,0);
          particles->C = &M(var,3,0);
          particles->I = &M(var,4,0);
        }

        for( i=0; i<n; i++ )
        {
            particles->X[i] = PartX[i]; 
            particles->Y[i] = PartY[i];
            particles->Z[i] = PartZ[i];
            particles->C[i] = PartC[i];
            particles->I[i] = PartI[i];
        }
        particles->SaveNofParticles = n;
    }

    quick  = (particles->LineStyle == line_style_line);
    quick |= epMouseDown && epMouseDownTakesTooLong;

    if ( quick )
    {
        if ( particles->Style == particle_style_vector )
        {
           gra_beg_lines();
        } else
        {
           gra_polygon_mode( GRA_LINE );
        }
        gra_sphere_quality( 1 );

        if ( !(epMouseDown && epMouseDownTakesTooLong) )
          gra_line_width( particles->LineWidth );
        else
          gra_line_width( 1.0 );
    }

    if ( !quick && (particles->LineStyle == line_style_cylinder) )
    {
        gra_sphere_quality( particles->LineQuality );
    }

    t = 0.0;
    for( i=0; i<n; i++ )
    {
        if ( particles->ParticleData[4]->f[i]<0 ) continue;
        
        x = PartX[i];
        y = PartY[i];
        z = PartZ[i];
        c = PartC[i];

        if ( particles->Advance && !epMouseDown )
        {
            particles->X[i] = x;
            particles->Y[i] = y;
            particles->Z[i] = z;
            particles->C[i] = c;
            particles->I[i] = PartI[i];

            t = 0.0;
            while( t < particles->OutDT-CEPS )
            {
                dt = MIN( particles->OutDT-t,particles->MaxDT );
                while( TRUE )
                {
                    if ( dt < 1.0E-9 ) break;

                    x0 = x1 = x;
                    y0 = y1 = y;
                    z0 = z1 = z;

                    switch( particles->IntegMethod )
                    {
                        case particle_integ_euler:
                           vis_Euler( particles, model, i, &x0, &y0, &z0, &c, t, dt, f_x0, f_y0, f_z0, f_c0 );
                           if ( particles->IntegPolicy==particle_policy_adaptive )
                           {
                               vis_Euler( particles, model, i, &x1, &y1, &z1, &c, t,  dt / 2, f_x0, f_y0, f_z0, f_c0 );
                               vis_Euler( particles, model, i, &x1, &y1, &z1, &c, t + dt / 2, dt / 2, f_x0, f_y0, f_z0, f_c0 );
                           }
                        break;

                        case particle_integ_runge_kutta:
                           vis_RungeKutta( particles, model, i, &x0, &y0, &z0, &c, t,  dt, f_x0, f_y0, f_z0, f_c0 );
                           if ( particles->IntegPolicy==particle_policy_adaptive )
                           {
                               vis_RungeKutta( particles, model, i, &x1, &y1, &z1, &c, t, dt / 2, f_x0, f_y0, f_z0, f_c0 );
                               vis_RungeKutta( particles, model, i, &x1, &y1, &z1, &c, t + dt / 2, dt / 2, f_x0, f_y0, f_z0, f_c0 );
                           }
                        break;
                    }

                    if ( particles->IntegPolicy == particle_policy_fixed )
                    {
                        x1 = x0;
                        y1 = y0;
                        z1 = z0;
                        break;
                    }

                    if ( ABS(x0-x1) < CEPS && ABS(y0-y1) < CEPS && ABS(z0-z1) < CEPS ) break;

                    dt /= 2;

                    if ( BreakLoop ) break;
                }
                if ( BreakLoop ) break;

                t += dt;
                x = x1;
                y = y1;
                z = z1;

                if ( dt < 1.0E-9 )
                {
                    fprintf( stderr, "For Your record: stepping didn't succeed.\n" );
                    break;
                }
                if ( BreakLoop ) break;
            }
            if ( BreakLoop ) break;

            PartX[i] = x;
            PartY[i] = y;
            PartZ[i] = z;
            PartC[i] = c;
        }
        if ( BreakLoop ) break;

        if ( particles->Style == particle_style_vector )
        {
            LX = particles->X[i];
            LY = particles->Y[i];
            LZ = particles->Z[i];
            LC = particles->C[i];

            LX = ( 2.0*(LX-xmin) - (xmax-xmin) ) / s;
            LY = ( 2.0*(LY-ymin) - (ymax-ymin) ) / s;
            LZ = ( 2.0*(LZ-zmin) - (zmax-zmin) ) / s;

            x = ( 2.0*(x-xmin) - (xmax-xmin) ) / s;
            y = ( 2.0*(y-ymin) - (ymax-ymin) ) / s;
            z = ( 2.0*(z-zmin) - (zmax-zmin) ) / s;

#if 0
if ( i==185 )
{
extern double viewx,viewy,viewz,tox,toy,toz,upx,upy,upz;
viewx=LX;
viewy=LY;
viewz=LZ;
tox=x;
toy=y;
toz=z;
upx=0;
if ( ABS(LX-x)>ABS(LY-y) || ABS(LZ-z)>ABS(LY-y) ) { upy=1; upz=0; }
else { upy=0; upz=1; }
continue;
}
if ( (i%10) ) continue;
#endif


            pnt[0] = LX;
            pnt[1] = LY;
            pnt[2] = LZ;

            vec[0] = x - LX;
            vec[1] = y - LY;
            vec[2] = z - LZ;

            r = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
            gra_arrow( pnt,vec,CScl*(c-CAdd),particles->LineStyle,
                   particles->ArrowStyle,0.1*particles->LineWidth*r );
        } else {
            x = ( 2.0*(x-xmin) - (xmax-xmin) ) / s;
            y = ( 2.0*(y-ymin) - (ymax-ymin) ) / s;
            z = ( 2.0*(z-zmin) - (zmax-zmin) ) / s;
            c = CScl*(c-CAdd);
            gra_sphere( x,y,z,c,particles->LineWidth*r*c );
        }

        if ( epMouseDown && ( i & 30 ) )
        {
            if ( RealTime() - gtime > TooLong2 )
                if ( ++epMouseDownTakesTooLong > 3 )
                {
                    if ( quick )
                        if ( particles->Style == particle_style_vector )
                        {
                             gra_end_lines( );
                        } else
                        {
                             gra_polygon_mode( GRA_FILL );
                        }
                    return FALSE;
                } else gtime = RealTime();
        }
    }

    particles->Advance = FALSE;

    if ( quick )
        if ( particles->Style == particle_style_vector )
        {
             gra_end_lines( );
        } else
        {
             gra_polygon_mode( GRA_FILL );
        }

    return TRUE;
}


/*******************************************************************************
 *
 *     Name:        vis_particle_alloc
 *
 *     Purpose:     allocate memory for particle_t structure
 *
 *     Parameters: 
 *
 *         Input:   none
 *
 *         Output:  none
 *   
 *   Return value:  pointer to allocated memory
 *
 ******************************************************************************/
static particle_t *vis_particle_alloc()
{
     particle_t *particle = (particle_t *)calloc(sizeof(particle_t),1);

     if ( !particle )
     {
         fprintf( stderr, "vis_particle_alloc: FATAL: can't alloc a few bytes of memory\n" );
     }

     return particle;
}

/*******************************************************************************
 *
 *     Name:        vis_particle_delete
 *
 *     Purpose:     free memory associated with particle_t structure
 *
 *     Parameters: 
 *
 *         Input:   (particle_t *) pointer to structure
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
static void vis_particle_delete(particle_t *particle)
{
    if ( particle ) free( particle );
}

/*******************************************************************************
 *
 *     Name:        vis_initialize_particle_visual
 *
 *     Purpose:     Register particle visual type
 *
 *     Parameters: 
 *
 *         Input:   none
 *
 *         Output:  none
 *   
 *   Return value:  vis_add_visual_type (malloc success probably)...
 *
 ******************************************************************************/
int vis_initialize_particle_visual()
{
    static char *visual_name = "Particles";

    visual_type_t VisualDef; 

    static particle_t particle;

    static visual_param_t ParticleParams[] =
    {
        { "VectorData1",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "VectorData2",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "VectorData3",   "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ParticleDataX", "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ParticleDataY", "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ParticleDataZ", "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ParticleDataC", "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ParticleDataI", "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "ColorData",     "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, NULL },
        { "Material",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultMaterial },
        { "ColorMap",      "%s", 0, VIS_VISUAL_PARAM_POINTER, 0, 0.0, &DefaultColorMap },
        { "LineQuality",   "%d", 0, VIS_VISUAL_PARAM_INT,     1, 0.0, NULL },
        { "LineStyle",     "%d", 0, VIS_VISUAL_PARAM_INT,     line_style_line, 0.0, NULL },
        { "LineWidth",     "%lf", 0, VIS_VISUAL_PARAM_FLOAT, 0, 1.0, NULL },
        { "ArrowStyle",    "%d", 0, VIS_VISUAL_PARAM_INT,     arrow_style_stick, 0.0, NULL },
        { "OutDT",         "%lf", 0, VIS_VISUAL_PARAM_FLOAT, 0, 1.0E-1, NULL },
        { "MaxDT",         "%lf", 0, VIS_VISUAL_PARAM_FLOAT, 0, 1.0E-3, NULL },
        { "Tolerance",     "%lf", 0, VIS_VISUAL_PARAM_FLOAT, 0, 1.0E-5, NULL },
        { "NofParticles",  "%d", 0, VIS_VISUAL_PARAM_INT, 0, 0.0, NULL },
        { "Advance",       "%d", 0, VIS_VISUAL_PARAM_INT, TRUE, 0.0, NULL },
        { "IntegMethod",   "%d", 0, VIS_VISUAL_PARAM_INT, particle_integ_runge_kutta, 0.0, NULL },
        { "IntegPolicy",   "%d", 0, VIS_VISUAL_PARAM_INT, particle_policy_adaptive, 0.0, NULL },
        { "Style",   "%d", 0, VIS_VISUAL_PARAM_INT, particle_style_vector, 0.0, NULL },
        { NULL, NULL, 0, 0, 0,0.0, NULL }
    };

    int n = 0;

    ParticleParams[n++].Offset = (char *)&particle.VectorData[0]   - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.VectorData[1]   - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.VectorData[2]   - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ParticleData[0] - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ParticleData[1] - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ParticleData[2] - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ParticleData[3] - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ParticleData[4] - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ColorData       - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.Material        - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ColorMap        - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.LineQuality     - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.LineStyle       - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.LineWidth       - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.ArrowStyle      - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.OutDT           - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.MaxDT           - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.Tolerance       - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.NofParticles    - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.Advance         - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.IntegMethod     - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.IntegPolicy     - (char *)&particle;
    ParticleParams[n++].Offset = (char *)&particle.Style           - (char *)&particle;

    VisualDef.VisualName    = visual_name;
    VisualDef.RealizeVisual = (int   (*)()) vis_particle;
    VisualDef.AllocParams   = (void *(*)()) vis_particle_alloc;
    VisualDef.DeleteParams  = (void  (*)()) vis_particle_delete;
    VisualDef.VisualParams  = ParticleParams;

    return vis_add_visual_type( &VisualDef );
}
