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
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax:  +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 27 Sep 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

#define MODULE_MAIN
#include "elmerpost.h"
#include "../config.h"

#include <tcl.h>
#include <tk.h>

static char str[512];

static Tcl_DString TCLCommand;
Tcl_Interp *TCLInterp,*TCLIntInterp;

#define MAX_COLORMAP_ENTRIES 512

double tx,ty,tz,ax,ay,az,sx=1.,sy=1.,sz=1.,br=0.0,bg=0.0,bb=0.0;
double XMin,XMax,YMin,YMax,ZMin,ZMax;

void opengl_draw(),Mouse();

static int UpdatePending = FALSE;

#ifndef MINGW32
#include <sys/time.h>
#endif

#define MAX_OBJECTS 10

static int ShowVectors,ShowContours,ShowColorMesh=1, ShowSpheres,
ShowMeshLines,ShowIsosurfaces,ShowParticles,ShowColorScale;

static int MeshStyle,MeshEdgeStyle=1,MeshLineStyle,MeshQuality = 1;
static double MeshRadius = 1.0;
static int MeshColorSetMinMax,MeshNodeNumbers=0;
static scalar_t MeshColor[MAX_OBJECTS];

static int ContourLines=5,ContourLineStyle,ContourQuality = 1,ContourColorSetMinMax;
static double ContourRadius = 1.0;
static scalar_t ContourColor[MAX_OBJECTS],ContourContour[MAX_OBJECTS];

static int VectorLineStyle,VectorQuality=1;
static double VectorRadius=1.0,VectorLengthScale=1.0;

static scalar_t VectorColor,VectorLength,VectorThreshold;
static scalar_t VectorArrow[4];
static double VectorCeiling,VectorFloor;

static int IsosurfaceStyle,IsosurfaceContours=1,IsosurfaceLineStyle,IsosurfaceQuality = 1;
static int IsosurfaceRecompute=TRUE, IsosurfaceContourSetMinMax;
static double IsosurfaceRadius = 1.0;
static scalar_t IsosurfaceContour,IsosurfaceColor,IsosurfaceNormal[4];
static int IsosurfaceColorSetMinMax;

static int ParticleLineStyle,ParticleQuality=1,ParticleNofParticles,ParticleArrowStyle;
static int ParticleAdvance=1;
static int ParticleIntegMethod=1, ParticleIntegPolicy=1,ParticleStyle;
static double ParticleRadius = 1.0, ParticleOutDT  = 1.0E-1,
              ParticleMaxDT = 1.0E-3, ParticleTolerance = 1.0E-5;
static scalar_t ParticleParticle[5],ParticleVelocity[4],ParticleColor;

static double SphereRadiusScale=1.0,SphereFloor,SphereCeiling;
static int SphereQuality=1;
static scalar_t SphereColor,SphereRadius,SphereThreshold;

static int ColorScaleEntries,ColorScaleDecimals,ColorScaleStyle,ColorScaleFontColor;
static scalar_t ColorScaleColor;
static double ColorScaleX,ColorScaleY,ColorScaleThickness,ColorScaleLength,ColorScaleFontSize = 17.0;
static int ColorScaleColorSetMinMax;

static int NumberOfScalarVariables, NumberOfVectorVariables, NumberOfParticleVariables;
int CurrentTimeStep, GraphicsClearOn = 1, KeepScale = 0, NormalUpdate = 1;

double xmin,xmax,ymin,ymax,zmin,zmax;

// Optional FTGL-text rendering functionality:
//=============================================
#if defined(HAVE_FTGL_NEW) || defined(HAVE_FTGL_OLD)
int FtFont(ClientData,Tcl_Interp*,int,char**);
int FtText(ClientData,Tcl_Interp*,int,char**);
#endif

static void fpe_sig( int sig ) { signal( SIGFPE, fpe_sig ); }

static void int_sig( int sig )
{
    fprintf( stdout, "\n^C\n" );

    BreakLoop = TRUE;
    signal( SIGINT, int_sig );
}

void DrawItSomeTimeWhenIdle()
{
    if ( !UpdatePending )
    {
        UpdatePending = TRUE;
        Tk_DoWhenIdle( (Tk_IdleProc *)opengl_draw,NULL );
    }
}

static void epSwapBuffers()
{
    auxSwapBuffers();
    UpdatePending = FALSE;
}

object_t *obj_find( object_t *, char * );
object_t *obj_add_object( object_t *, char * );

static int SetObject( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   object_t *obj;

   obj = (object_t *)obj_find( &VisualObject,argv[1] );
   if ( !obj ) obj = (object_t *)obj_add_object( &VisualObject, argv[1] );

   Tcl_UnlinkVar( TCLInterp, "MeshColorMin" );
   Tcl_UnlinkVar( TCLInterp, "MeshColorMax" );

   Tcl_LinkVar( TCLInterp,"MeshColorMin",(char *)&MeshColor[obj->Id].min, TCL_LINK_DOUBLE );
   Tcl_LinkVar( TCLInterp,"MeshColorMax",(char *)&MeshColor[obj->Id].max, TCL_LINK_DOUBLE );

   Tcl_UnlinkVar( TCLInterp, "ContourColorMin" );
   Tcl_UnlinkVar( TCLInterp, "ContourColorMax" );

   Tcl_LinkVar( TCLInterp,"ContourColorMin",(char *)&ContourColor[obj->Id].min, TCL_LINK_DOUBLE );
   Tcl_LinkVar( TCLInterp,"ContourColorMax",(char *)&ContourColor[obj->Id].max, TCL_LINK_DOUBLE );

   CurrentObject = obj;
   Tcl_Eval( TCLInterp, "math who" );

   return TCL_OK;
}

static int SetParentObject( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   object_t *obj;

   obj = (object_t *)obj_find( &VisualObject,argv[1] );
   if ( !obj )
   {
      sprintf( interp->result, "parent: no such object [%s]", argv[1] );
      return TCL_ERROR;
   }

   obj_set_parent( CurrentObject, obj);
   return TCL_OK;
}

static int GetInterpolate( ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int r0,g0,b0,r1,g1,b1,n0,n1,i;
    double r,g,b,dr,dg,db;

    static char result[7*MAX_COLORMAP_ENTRIES],color[32];

    sscanf( argv[1], "%d", &n0 );
    sscanf( argv[2], "%d", &n1 );

    sscanf( argv[3], "#%02x%02x%02x", &r0,&g0,&b0 );
    sscanf( argv[4], "#%02x%02x%02x", &r1,&g1,&b1 );

    if ( n0 < n1 )
    {
        dr = (r1-r0)/(double)(n1-n0);
        dg = (g1-g0)/(double)(n1-n0);
        db = (b1-b0)/(double)(n1-n0);
    }

    r = r0;
    g = g0;
    b = b0;

    result[0] = '\0';
    for( i=n0; i<n1; i++ )
    {
        sprintf( color, "#%02x%02x%02x ",(int)r,(int)g,(int)b );
        strcat( result, color );

        r += dr;
        g += dg;
        b += db;
    }

    Tcl_SetResult(interp,result,TCL_STATIC);

    return TCL_OK;
}

static int GetColorMap(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int first = TRUE,rr, gg, bb, i, n;
    char *str;
    colormap_t *cmp = &DefaultColorMap;

    if ( strncmp( argv[1], "-vectors", 8 ) == 0 ) {
       if ( ArrowColorMap != &DefaultColorMap ) {
           free( ArrowColorMap->Values );
           free( ArrowColorMap );
       }
       cmp = ArrowColorMap = (colormap_t *)malloc( sizeof(colormap_t) );
       cmp->Values = NULL;
       argv++;
    }
    else if ( strncmp( argv[1], "-mesh", 5 ) == 0 ) {
       if ( MeshColorMap != &DefaultColorMap ) {
           free( MeshColorMap->Values );
           free( MeshColorMap );
       }
       cmp = MeshColorMap = (colormap_t *)malloc( sizeof(colormap_t) );
       cmp->Values = NULL;
       argv++;
    }
    else if ( strncmp( argv[1], "-contour", 8 ) == 0 ) {
       if ( ContourColorMap != &DefaultColorMap ) {
           free( ContourColorMap->Values );
           free( ContourColorMap );
       }
       cmp = ContourColorMap = (colormap_t *)malloc( sizeof(colormap_t) );
       cmp->Values = NULL;
       argv++;
    }
    else if ( strncmp( argv[1], "-isosurface", 11 ) == 0 ) {
       if ( IsoSurfaceColorMap != &DefaultColorMap ) {
           free( IsoSurfaceColorMap->Values );
           free( IsoSurfaceColorMap );
       }
       cmp = IsoSurfaceColorMap = (colormap_t *)malloc( sizeof(colormap_t) );
       cmp->Values = NULL;
       argv++;
    }
    else if ( strncmp( argv[1], "-sphere", 7 ) == 0 ) {
       if ( SphereColorMap != &DefaultColorMap ) {
           free( SphereColorMap->Values );
           free( SphereColorMap );
       }
       argv++;
       cmp = SphereColorMap = (colormap_t *)malloc( sizeof(colormap_t) );
       cmp->Values = NULL;
    }
    else if ( strncmp( argv[1], "-particles", 9 ) == 0 ) {
       if ( ParticleColorMap != &DefaultColorMap ) {
           free( ParticleColorMap->Values );
           free( ParticleColorMap );
       }
       argv++;
       cmp = ParticleColorMap = (colormap_t *)malloc( sizeof(colormap_t) );
       cmp->Values = NULL;
    } else if ( *argv[1] == '-' ) {
          argv++;
    }

    sscanf( argv[1], "%d", &n );
    if ( cmp->Values ) free( cmp->Values );

    cmp->Values = (rgb_t *)calloc(n,sizeof(rgb_t));

    if ( !cmp->Values ) return TCL_ERROR;

    str = argv[2];
    for( i=0; i<n; i++ )
    {
        sscanf( &str[7*i], "#%02x%02x%02x", &rr,&gg,&bb );
        cmp->Values[i].r = rr;
        cmp->Values[i].g = gg;
        cmp->Values[i].b = bb;
    }

    cmp->Changed = TRUE;
    cmp->NumberOfEntries = n;

    if ( first )
    {
        def_map.Values = (rgb_t *)malloc( n*sizeof(rgb_t) );
        memcpy( def_map.Values,DefaultColorMap.Values,n*sizeof(rgb_t));
        def_map.Changed = TRUE;
        def_map.NumberOfEntries = n;
    }

    UpdateObject();
    DrawItSomeTimeWhenIdle();

    return TCL_OK;
}

static int StopProcessing(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    int_sig( 0 );

    return TCL_OK;
}

static int RecomputeNormals(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{

    geo_vertex_normals( CurrentObject->Geometry, 50.0 );
    return TCL_OK;
}

static int ClipPlane( ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    int i,id;

    if ( argc != 6 ) {
        sprintf( interp->result, "clip: Wrong number of arguments.\n" );
        return TCL_ERROR;
    }

    sscanf( argv[1], "%d",  &id );
    if ( id < 0 || id > 5 ) {
        sprintf( interp->result, "clip: Plane argument should be from a int from 0..5.\n" );
        return TCL_ERROR;
    }
    CurrentObject->ClipPlane[id] = id;
    sscanf( argv[2], "%lf", &CurrentObject->ClipEquation[id][0] );
    sscanf( argv[3], "%lf", &CurrentObject->ClipEquation[id][1] );
    sscanf( argv[4], "%lf", &CurrentObject->ClipEquation[id][2] );
    sscanf( argv[5], "%lf", &CurrentObject->ClipEquation[id][3] );

    return TCL_OK;
}

static int UpdateBackColor(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    argv++; argc--;
    if ( strcmp( *argv, "color") != 0 ) return TCL_ERROR;

    argv++; argc--;
    if ( argc <= 0 ) return TCL_ERROR;

    argv++; argc--;
    if ( argc ) sscanf( *argv, "%lf", &br );

    argv++; argc--;
    if ( argc ) sscanf( *argv, "%lf", &bg );

    argv++; argc--;
    if ( argc ) sscanf( *argv, "%lf", &bb );

    br /= 100.0;
    bg /= 100.0;
    bb /= 100.0;

    DrawItSomeTimeWhenIdle();

    return TCL_OK;
}

static int UpdateColor(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int diffuse  = FALSE;

    static double r=1.,g=1.,b=1.,alpha=1.0,shininess=0.0,spec_r,spec_g,spec_b;

    double rr,gg,bb;

    material_t *mat = &DefaultMaterial;

    if ( strncmp( argv[1], "-vectors", 8 ) == 0 ) {
       if ( ArrowMaterial != &DefaultMaterial ) {
           free( ArrowMaterial );
       }
       mat = ArrowMaterial = (material_t *)malloc( sizeof(material_t) );
       argv++;
       argc--;
    }
    else if ( strncmp( argv[1], "-mesh", 5 ) == 0 ) {
       if ( MeshMaterial != &DefaultMaterial ) {
           free( MeshMaterial );
       }
       mat = MeshMaterial = (material_t *)malloc( sizeof(material_t) );
       argv++;
       argc--;
    }
    else if ( strncmp( argv[1], "-contour", 8 ) == 0 ) {
       if ( ContourMaterial != &DefaultMaterial ) {
           free( ContourMaterial );
       }
       mat = ContourMaterial = (material_t *)malloc( sizeof(material_t) );
       argv++;
       argc--;
    }
    else if ( strncmp( argv[1], "-isosurface", 11 ) == 0 ) {
       if ( IsoSurfaceMaterial != &DefaultMaterial ) {
           free( IsoSurfaceMaterial );
       }
       mat = IsoSurfaceMaterial = (material_t *)malloc( sizeof(material_t) );
       argv++;
       argc--;
    }
    else if ( strncmp( argv[1], "-sphere", 7 ) == 0 ) {
       if ( SphereMaterial != &DefaultMaterial ) {
           free( SphereMaterial );
       }
       argv++;
       argc--;
       mat = SphereMaterial = (material_t *)malloc( sizeof(material_t) );
    }
    else if ( strncmp( argv[1], "-particles", 9 ) == 0 ) {
       if ( ParticleMaterial != &DefaultMaterial ) {
           free( ParticleMaterial );
       }
       argv++;
       argc--;
       mat = ParticleMaterial = (material_t *)malloc( sizeof(material_t) );
    } else if ( strncmp( argv[1],  "-default",8 ) ) {
       argv++;
       argc--;
    }

    while( argc )
    {
        argv++; argc--;
        if ( argc )
        {
            if ( strcmp( *argv, "opacity" ) == 0 )
            {
                argv++; argc--;
                if ( argc ) { sscanf( *argv, "%lf", &alpha ); alpha /= 100; } else return TCL_ERROR;
            } else  if ( strcmp( *argv, "color" ) == 0 )
            {
                argv++; argc--;
                if ( argc <= 0 ) return TCL_ERROR;

                diffuse = FALSE;
                if ( strcmp( *argv,"diffuse" ) == 0 ) diffuse = TRUE;

                argv++; argc--;
                if ( argc > 0 ) sscanf( *argv, "%lf", &rr ); else return TCL_ERROR;

                argv++; argc--;
                if ( argc > 0 ) sscanf( *argv, "%lf", &gg );

                argv++; argc--;
                if ( argc > 0 ) sscanf( *argv, "%lf", &bb );

                if ( diffuse )
                {
                    r = rr / 100.0;
                    g = gg / 100.0;
                    b = bb / 100.0;
                } else
                {
                    spec_r = rr / 100.0;
                    spec_g = gg / 100.0;
                    spec_b = bb / 100.0;
                }
            } else if ( strcmp( *argv, "shininess" ) == 0 )
            {
                argv++; argc--;
                if ( argc ) { sscanf( *argv, "%lf", &shininess ); }  else return TCL_ERROR;
            }
        }
    }

    mat->Changed  = TRUE;

    mat->Shininess = shininess;

    mat->Diffuse[0] = r;
    mat->Diffuse[1] = g;
    mat->Diffuse[2] = b;
    mat->Diffuse[3] = alpha;

    mat->Specular[0] = spec_r;
    mat->Specular[1] = spec_g;
    mat->Specular[2] = spec_b;
    mat->Specular[3] = alpha;
    mat->Changed  = TRUE;


    UpdateObject();
    DrawItSomeTimeWhenIdle();

    return TCL_OK;
}

static int UpdateEdgeColor(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int diffuse  = FALSE;

    static double r=1.,g=1.,b=1.,alpha=1.0,shininess=0.0,spec_r,spec_g,spec_b;

    double rr,gg,bb;

    while( argc )
    {
        argv++; argc--;
        if ( argc )
        {
            if ( strcmp( *argv, "opacity" ) == 0 )
            {
                argv++; argc--;
                if ( argc ) { sscanf( *argv, "%lf", &alpha ); alpha /= 100; } else return TCL_ERROR;
            } else  if ( strcmp( *argv, "color" ) == 0 )
            {
                argv++; argc--;
                if ( argc <= 0 ) return TCL_ERROR;

                diffuse = FALSE;
                if ( strcmp( *argv,"diffuse" ) == 0 ) diffuse = TRUE;

                argv++; argc--;
                if ( argc > 0 ) sscanf( *argv, "%lf", &rr ); else return TCL_ERROR;

                argv++; argc--;
                if ( argc > 0 ) sscanf( *argv, "%lf", &gg );

                argv++; argc--;
                if ( argc > 0 ) sscanf( *argv, "%lf", &bb );

                if ( diffuse )
                {
                    r = rr / 100.0;
                    g = gg / 100.0;
                    b = bb / 100.0;
                } else
                {
                    spec_r = rr / 100.0;
                    spec_g = gg / 100.0;
                    spec_b = bb / 100.0;
                }
            } else if ( strcmp( *argv, "shininess" ) == 0 )
            {
                argv++; argc--;
                if ( argc ) { sscanf( *argv, "%lf", &shininess ); }  else return TCL_ERROR;
            }
        }
    }

    DefaultEdgeMaterial.Changed  = TRUE;

    DefaultEdgeMaterial.Shininess = shininess;

    DefaultEdgeMaterial.Diffuse[0] = r;
    DefaultEdgeMaterial.Diffuse[1] = g;
    DefaultEdgeMaterial.Diffuse[2] = b;
    DefaultEdgeMaterial.Diffuse[3] = alpha;

    DefaultEdgeMaterial.Specular[0] = spec_r;
    DefaultEdgeMaterial.Specular[1] = spec_g;
    DefaultEdgeMaterial.Specular[2] = spec_b;
    DefaultEdgeMaterial.Specular[3] = alpha;

    DrawItSomeTimeWhenIdle();

    return TCL_OK;
}


static void GetScalarVariable( scalar_t *Func, char *name, element_model_t *model, int TimeAll, int SetMinMax )
{
    char str[256];

    VARIABLE *var;

    double *T;

    int len = strlen(name)-1, index = 0, n, t, i;

    strcpy(str,name);

    if ( str[len-1] == '_' )
    {
        index = str[len] - 'x';
        str[len-1] = 0;
    }

    Func->f = NULL;

    if ( !(var = (VARIABLE *)lst_find( VARIABLES, str ) ) ) return;

    n = NCOL(var) % model->NofNodes;
    t = NCOL(var) / model->NofNodes - 1;

    if ( n || t < 0 ) return;

    T = &M(var,index,0);
    if ( TimeAll )
      Func->f = T;
    else
      Func->f = &T[MIN(CurrentTimeStep,t)*model->NofNodes];

    if ( !SetMinMax ) return;

    Func->min =  1.0e20;
    Func->max = -1.0e20;
    for( i=0; i<NCOL(var); i++ )
    {
        Func->min = MIN( Func->min, T[i] );
        Func->max = MAX( Func->max, T[i] );
    }

   if ( Func->name ) free( Func->name );
   Func->name = (char *)malloc( len+2 );
   strcpy( Func->name, name );
}

static void GetVectorVariable( scalar_t Func[4], char *name, element_model_t *model, int TimeAll )
{
    VARIABLE *var;

    double *T;

    int i,n,t,index;
    Func[0].f = Func[1].f =
        Func[2].f = Func[3].f = NULL;

    if ( !(var = (VARIABLE *)lst_find( VARIABLES, name ) ) ) return;

    n = NCOL(var) % model->NofNodes;
    t = NCOL(var) / model->NofNodes - 1;

    if ( n || t < 0 ) return;

    if ( TimeAll )
    {
        Func[1].f = &M(var,0,0);
        Func[2].f = &M(var,1,0);
        Func[3].f = &M(var,2,0);
    } else
    {
        T = &M(var,0,0);
        Func[1].f = &T[MIN(CurrentTimeStep,t)*model->NofNodes];

        T = &M(var,1,0);
        Func[2].f = &T[MIN(CurrentTimeStep,t)*model->NofNodes];

        T = &M(var,2,0);
        Func[3].f = &T[MIN(CurrentTimeStep,t)*model->NofNodes];
    }

    for( n=1; n<4; n++ )
    {
        Func[n].min =  1.0e20;
        Func[n].max = -1.0e20;
        for( i=0; i<NCOL(var); i++ )
        {
            Func[n].min = MIN( Func[n].min, M(var,n-1,i) );
            Func[n].max = MAX( Func[n].max, M(var,n-1,i) );
        }
    }
}

static void GetParticleVariable( scalar_t Func[5], char *name )
{
    VARIABLE *var;

    double *T;

    int i,n,index;

    var = (VARIABLE *)lst_find( VARIABLES, name );

    Func[0].f =
        Func[1].f =
        Func[2].f =
        Func[3].f =
        Func[4].f = NULL;

    if ( !var || NROW(var) != 5 ) return;

    Func[0].f = &M(var,0,0);
    Func[1].f = &M(var,1,0);
    Func[2].f = &M(var,2,0);
    Func[3].f = &M(var,3,0);
    Func[4].f = &M(var,4,0);

    ParticleNofParticles = NCOL(var);
}


static int UpdateVariable( ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    element_model_t *model = CurrentObject->ElementModel;

    if ( !model ) return TCL_OK;

    if ( strcmp( argv[1], "VectorColor") == 0 )

        GetScalarVariable( &VectorColor, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, TRUE );

    else if ( strcmp( argv[1], "VectorLength" ) == 0 )

        GetScalarVariable( &VectorLength, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, TRUE );

    else if ( strcmp( argv[1], "VectorThreshold" ) == 0 )

        GetScalarVariable( &VectorThreshold, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, TRUE );

    else if ( strcmp( argv[1], "MeshColor" ) == 0 )

        GetScalarVariable( &MeshColor[CurrentObject->Id],
            (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model,0,!MeshColorSetMinMax );

   else if ( strcmp( argv[1], "ColorScaleColor" ) == 0 )

        GetScalarVariable( &ColorScaleColor, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, !ColorScaleColorSetMinMax );

    else if ( strcmp( argv[1], "ContourColor" ) == 0 )

        GetScalarVariable( &ContourColor[CurrentObject->Id],
             (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, !ContourColorSetMinMax );

    else if ( strcmp( argv[1], "ContourContour" ) == 0 )

        GetScalarVariable( &ContourContour[0], (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, TRUE );

    else if ( strcmp( argv[1], "IsosurfaceColor" ) == 0 )

        GetScalarVariable( &IsosurfaceColor, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, !IsosurfaceColorSetMinMax );

    else if ( strcmp( argv[1], "IsosurfaceContour" ) == 0 )

        GetScalarVariable( &IsosurfaceContour, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, !IsosurfaceContourSetMinMax );

    else if ( strcmp( argv[1], "ParticleColor" ) == 0 )

        GetScalarVariable( &ParticleColor, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, TRUE );

    else if ( strcmp( argv[1], "SphereThreshold" ) == 0 )

        GetScalarVariable( &SphereThreshold, (char *)Tcl_GetVar( TCLInterp,argv[1],TCL_GLOBAL_ONLY ),model, 0, TRUE );

    return TCL_OK;
}


static int TimeStep(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    double s;
    static int n,i,j;

    element_model_t *ElementModel = CurrentObject->ElementModel;

    if ( argc != 2 )
    {
        sprintf( interp->result, "c_TimeStep: wrong number of parameters.\n" );
        return TCL_ERROR;
    }

    if ( !ElementModel ) return TCL_OK;

    n = 0;
    sscanf( argv[1], "%d", &n );

    if ( n < 0  || n >= ElementModel->NofTimesteps )
    {
        sprintf( interp->result, "c_TimeStep: Invalid timestep.\n" );
        return TCL_ERROR;
    }

    CurrentTimeStep = n;

    {
        static vertex_t *vertex;

        vertex = CurrentObject->Geometry->Vertices;

        if ( !KeepScale ) {
           xmin = ymin = zmin =  1.0e20;
           xmax = ymax = zmax = -1.0e20;

           for( i=0; i<CurrentObject->ElementModel->NofNodes; i++ )
           {
               xmin = MIN( xmin,ElementModel->NodeArray[i] );
               ymin = MIN( ymin,ElementModel->NodeArray[ElementModel->NofNodes+i] );
               zmin = MIN( zmin,ElementModel->NodeArray[2*ElementModel->NofNodes+i] );

               xmax = MAX( xmax,ElementModel->NodeArray[i] );
               ymax = MAX( ymax,ElementModel->NodeArray[ElementModel->NofNodes+i] );
               zmax = MAX( zmax,ElementModel->NodeArray[2*ElementModel->NofNodes+i] );
           }

           s = MAX(MAX(xmax-xmin,ymax-ymin),zmax-zmin);

           CurrentObject->Geometry->Scale = s;
           CurrentObject->Geometry->MinMax[0].x[0] = xmin;
           CurrentObject->Geometry->MinMax[0].x[1] = ymin;
           CurrentObject->Geometry->MinMax[0].x[2] = zmin;

           CurrentObject->Geometry->MinMax[1].x[0] = xmax;
           CurrentObject->Geometry->MinMax[1].x[1] = ymax;
           CurrentObject->Geometry->MinMax[1].x[2] = zmax;
        } else {
           s = CurrentObject->Geometry->Scale;
           xmin = CurrentObject->Geometry->MinMax[0].x[0];
           ymin = CurrentObject->Geometry->MinMax[0].x[1];
           zmin = CurrentObject->Geometry->MinMax[0].x[2];

           xmax = CurrentObject->Geometry->MinMax[1].x[0];
           ymax = CurrentObject->Geometry->MinMax[1].x[1];
           zmax = CurrentObject->Geometry->MinMax[1].x[2];
        }

        for( i=0; i<ElementModel->NofNodes; i++,vertex++ )
        {
            vertex->x[0] = (2.0 * (ElementModel->NodeArray[i]-xmin)-(xmax-xmin))/s;
            vertex->x[1] = (2.0 * (ElementModel->NodeArray[ElementModel->NofNodes+i]-ymin)-(ymax-ymin)) / s;
            vertex->x[2] = (2.0 * (ElementModel->NodeArray[2*ElementModel->NofNodes+i]-zmin)-(zmax-zmin)) / s;
        }
    }


    /*
    * DrawItSomeTimeWhenIdle();
        * opengl_draw();
    */

    return TCL_OK;
}

static void cam_set_current( camera_t *camera )
{
    static char str[200];

    Tcl_SetVar2( TCLInterp, "camera","name", camera->Name, TCL_GLOBAL_ONLY );

    Tcl_SetVar2( TCLInterp, "camera","display", "on", TCL_GLOBAL_ONLY );

    if ( camera->ProjectionType == camera_proj_ortho )
        Tcl_SetVar2( TCLInterp, "camera","projection", "ortho",TCL_GLOBAL_ONLY );
    else
        Tcl_SetVar2( TCLInterp, "camera","projection", "perspective",TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->LookFromX );
    Tcl_SetVar2( TCLInterp, "camera","from_x", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->LookFromY );
    Tcl_SetVar2( TCLInterp, "camera","from_y", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->LookFromZ );
    Tcl_SetVar2( TCLInterp, "camera","from_z", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->LookAtX );
    Tcl_SetVar2( TCLInterp, "camera","to_x", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->LookAtY );
    Tcl_SetVar2( TCLInterp, "camera", "to_y", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->LookAtZ );
    Tcl_SetVar2( TCLInterp, "camera","to_z", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->UpX );
    Tcl_SetVar2( TCLInterp, "camera","up_x", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->UpY );
    Tcl_SetVar2( TCLInterp, "camera","up_y", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->UpZ );
    Tcl_SetVar2( TCLInterp, "camera","up_z", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->FieldAngle );
    Tcl_SetVar2( TCLInterp, "camera","field_angle", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->ClipNear );
    Tcl_SetVar2( TCLInterp, "camera","near", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->ClipFar );
    Tcl_SetVar2( TCLInterp, "camera","far", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->ViewportLowX );
    Tcl_SetVar2( TCLInterp, "camera","view_low_x", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->ViewportLowY );
    Tcl_SetVar2( TCLInterp, "camera","view_low_y", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->ViewportHighX );
    Tcl_SetVar2( TCLInterp, "camera","view_high_x", str, TCL_GLOBAL_ONLY );

    sprintf( str, "%g", camera->ViewportHighY );
    Tcl_SetVar2( TCLInterp, "camera","view_high_y", str, TCL_GLOBAL_ONLY );
}

/*******************************************************************************
*
*     Name:        cam_load_cameras
*
*     Purpose:     Initialize camerass from default cam setting file
*
*     Parameters:
*
*         Input:   none
*
*         Output:  none
*
*   Return value:  New list of (camera_t *)cameras
*
*******************************************************************************/
static camera_t *cam_load_cameras( camera_t *camera,char *FileName )
{
    static char name[200],str[200],variable[64],value[64];
    FILE *fp = NULL;

    int i;

    float x,y,z,w;

    camera_t *first_cam = camera;

    if ( !FileName || !(fp = fopen( FileName,"r") ) )
    {
        sprintf( name, "%s/lib/cameras/default", getenv("ELMER_POST_HOME") );

        if ( !(fp = fopen( name, "r" ) ) )
        {
            camera = (camera_t *)cam_add_camera( NULL, "camera1" );
            if ( !first_cam ) first_cam = camera;

            cam_set_viewport( camera,0.0,0.0,1.0,1.0 );
            cam_set_look_to( camera,0.0,0.0,0.0,FALSE );
            cam_set_up_vector( camera,0.0,1.0,0.0 );
            cam_set_look_from( camera,0.0,0.0,5.0,FALSE );
            cam_set_clip( camera,1.0,20.0 );
            cam_set_onoff( camera, TRUE );
        }
     }


    if ( fp ) {
        while( fgets( str,100, fp ) )
        {
            if ( sscanf( str, "camera %s", name ) == 1 )
            {
                camera = (camera_t *)cam_add_camera( camera, name );
                if ( !first_cam ) first_cam = camera;
                cam_set_onoff( camera, TRUE );
            } else if ( sscanf( str, "from %f %f %f", &x,&y,&z ) == 3 )
            {
                cam_set_look_from( camera,x,y,z,FALSE );

            } else if ( sscanf( str, "to %f %f %f", &x,&y,&z ) == 3)
            {
                cam_set_look_to( camera,x,y,z,FALSE );

            } else if ( sscanf( str, "viewport %f %f %f %f", &x,&y,&z,&w ) == 4 )
            {
                cam_set_viewport( camera,x,y,z,w );

            } else if ( sscanf( str, "up %f %f %f", &x, &y, &z ) == 3 )
            {
                cam_set_up_vector( camera, x,y,z );

            } else if ( sscanf( str, "projection %s", name ) == 1 )
            {
                if ( strcmp( name,"perspective" ) == 0 )
                    cam_set_projection( camera, camera_proj_perspective );
                else
                    cam_set_projection( camera, camera_proj_ortho );
            } else if ( sscanf( str, "field angle %f", &x ) == 1 )
            {
                cam_set_field_angle( camera, x );
            } else if ( sscanf( str, "clip %f %f", &x, &y ) == 2 )
            {
                cam_set_clip( camera,x,y );
            }
        }
    }

    if ( camera ) cam_set_current( camera );
    
    if ( fp ) 
    {
      fclose( fp );
    }

    i = 0;
    for( camera=first_cam; camera != NULL; camera=camera->Next,i++ )
    {
        sprintf( str, "%d", i );
        Tcl_SetVar2( TCLInterp, "camera_names",str, camera->Name, TCL_GLOBAL_ONLY );
    }
    sprintf( str, "%d", i );
    Tcl_SetVar( TCLInterp, "camera_n", str, TCL_GLOBAL_ONLY );

    return first_cam;
}


static int CurrentCamera(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    camera_t *cam;
    int i;

    if ( argc != 2 )
    {
        sprintf( interp->result, "camera: wrong number of arguments\n" );
        return TCL_ERROR;
    }

    i = 0;
    for( cam=Camera; cam!=NULL; cam=cam->Next,i++ )
    {
        if ( strcmp( cam->Name, argv[1] ) == 0 )
        {
            cam_set_current( cam );
            sprintf( interp->result, "%d", i );
            return TCL_OK;
        }
    }

    return TCL_ERROR;
}


static int SetCamera(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static char str[200];

    float x,y,z,w;

    camera_t *camera = NULL;

    camera = (camera_t *)cam_add_camera( Camera,
           (char *)Tcl_GetVar2( TCLInterp,"camera","name", TCL_GLOBAL_ONLY ) );

    if ( camera ) {
        sscanf( Tcl_GetVar2( TCLInterp, "camera","display", TCL_GLOBAL_ONLY ), "%s", str );
        if ( strcmp( str, "off" ) == 0 )
            cam_set_onoff( camera, FALSE );
        else
            cam_set_onoff( camera, TRUE );

        sscanf( Tcl_GetVar2( TCLInterp, "camera","from_x", TCL_GLOBAL_ONLY ), "%f", &x );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","from_y", TCL_GLOBAL_ONLY ), "%f", &y );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","from_z", TCL_GLOBAL_ONLY ), "%f", &z );
        cam_set_look_from( camera,x,y,z,FALSE );

        sscanf( Tcl_GetVar2( TCLInterp, "camera","to_x", TCL_GLOBAL_ONLY  ), "%f", &x );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","to_y", TCL_GLOBAL_ONLY  ), "%f", &y );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","to_z", TCL_GLOBAL_ONLY  ), "%f", &z );
        cam_set_look_to( camera,x,y,z,FALSE );

        sscanf( Tcl_GetVar2( TCLInterp, "camera","up_x", TCL_GLOBAL_ONLY  ), "%f", &x );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","up_y", TCL_GLOBAL_ONLY  ), "%f", &y );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","up_z", TCL_GLOBAL_ONLY  ), "%f", &z );
        cam_set_up_vector( camera,x,y,z );

        sscanf( Tcl_GetVar2( TCLInterp, "camera","view_low_x", TCL_GLOBAL_ONLY  ),  "%f", &x );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","view_low_y", TCL_GLOBAL_ONLY  ),  "%f", &y );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","view_high_x", TCL_GLOBAL_ONLY  ), "%f", &z );
        sscanf( Tcl_GetVar2( TCLInterp, "camera","view_high_y", TCL_GLOBAL_ONLY  ), "%f", &w );
        cam_set_viewport( camera,x,y,z,w );

        sscanf( Tcl_GetVar2( TCLInterp,"camera","projection", TCL_GLOBAL_ONLY  ), "%s", str );
        if ( strcmp( str,"perspective" ) == 0 )
            cam_set_projection( camera, camera_proj_perspective );
        else
            cam_set_projection( camera, camera_proj_ortho );

        sscanf( Tcl_GetVar2( TCLInterp,"camera","field_angle", TCL_GLOBAL_ONLY ), "%f", &x );
        cam_set_field_angle( camera, x );

        sscanf( Tcl_GetVar2( TCLInterp,"camera","near", TCL_GLOBAL_ONLY ), "%f", &x );
        sscanf( Tcl_GetVar2( TCLInterp,"camera","far", TCL_GLOBAL_ONLY ),  "%f", &y );
        cam_set_clip( camera, x,y );
    }

    DrawItSomeTimeWhenIdle();

    return TCL_OK;
}


static int LoadCamera(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    if ( argc != 2 )
    {
        fprintf( stderr, "camera: wrong number of arguments\n" );
        return TCL_ERROR;
    }

    cam_delete_list( Camera );
    Camera = (camera_t *)cam_load_cameras( NULL,argv[1] );

    DrawItSomeTimeWhenIdle();

    return TCL_OK;
}

static int GroupDisplay(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    int i,gid,groupid;
    group_t *group;

    if ( argc > 1 )
    {
       if ( strcmp( argv[2],"on" ) == 0 )
       {
        for( group=CurrentObject->ElementModel->Groups; group!=NULL; group=group->Next )
                     if ( strcmp( argv[1],group->Name ) == 0 ) {
                            group->status = 1;
                            break;
                        }
       } else {
        for( group=CurrentObject->ElementModel->Groups; group!=NULL; group=group->Next )
            if ( strcmp( argv[1],group->Name ) == 0 ) {
                            group->status = 0;
                            break;
                        }

       }
    }

    for( i=0; i<CurrentObject->ElementModel->NofElements; i++ )
    {
        CurrentObject->ElementModel->Elements[i].DisplayFlag = FALSE;
        for ( gid=0; gid<MAX_GROUP_IDS; gid++ )
        {
           groupid = 0;
           for( group=CurrentObject->ElementModel->Groups; group!=NULL; group=group->Next,groupid++ )
              if ( groupid==CurrentObject->ElementModel->Elements[i].GroupIds[gid] ) break;

           if ( group ) CurrentObject->ElementModel->Elements[i].DisplayFlag |= group->status;
    }
    }

    return TCL_OK;
}

char *mtc_domath( char * );

static int MathCommand(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    VARIABLE *var;
    LIST *lst;

    static char name[256],value[256],*result;

    int i,j,k,NodesChanged;

    element_model_t *ElementModel = CurrentObject->ElementModel;

    /*
     * do the given math command first ...
     *
     */
    var_reset_status( "nodes" );

    result = (char *)mtc_domath( argv[1] );

    NodesChanged = var_get_status( "nodes" );


    /*
     * then check the variables which are good for then
     * element model...
     *
     */
    if ( result && strncmp( result, "MATC ERROR: ", 12 ) == 0 )
    {
      Tcl_SetResult( interp, result, TCL_STATIC );
      return TCL_ERROR;
    }


    /*
    *   if empty list return.
    */
    if ( listheaders[VARIABLES].next == (LIST *)NULL )
    {
      Tcl_SetResult( interp, result, TCL_STATIC );
      return TCL_OK;
    }

    i = j = k = 0;
    for( lst = listheaders[VARIABLES].next; lst; lst = NEXT(lst))
    {
      var = (VARIABLE *)lst;
      if ( !NAME(var)) continue;

      if ( NROW(var) == 5 )
      {
        sprintf( name,  "%d", k++ );
        sprintf( value, "%s", NAME(var) );
        Tcl_SetVar2( TCLInterp,"ParticleVariableNames", name, value, TCL_GLOBAL_ONLY );
      }

      if ( !ElementModel ) continue;

      if ( NROW(var) == 3  )
      {
        if ( !(NCOL(var) % ElementModel->NofNodes) )
        {
          sprintf( name,  "%d", i++ );
          sprintf( value, "%s_x", NAME(var) );
          Tcl_SetVar2( TCLInterp, "ScalarVariableNames", name, value, TCL_GLOBAL_ONLY );

          sprintf( name,  "%d", i++ );
          sprintf( value, "%s_y", NAME(var) );
          Tcl_SetVar2( TCLInterp, "ScalarVariableNames", name, value, TCL_GLOBAL_ONLY );

          sprintf( name,  "%d", i++ );
          sprintf( value, "%s_z", NAME(var) );
          Tcl_SetVar2( TCLInterp, "ScalarVariableNames", name, value, TCL_GLOBAL_ONLY );

          sprintf( name,  "%d", j++ );
          sprintf( value, "%s", NAME(var) );
          Tcl_SetVar2( TCLInterp, "VectorVariableNames", name, value, TCL_GLOBAL_ONLY );
        }
      }
      else if ( NROW(var) == 1 && !(NCOL(var) % ElementModel->NofNodes) )
      {
        sprintf( name,  "%d", i++ );
        Tcl_SetVar2( TCLInterp, "ScalarVariableNames", name, NAME(var), TCL_GLOBAL_ONLY );
      }
    }


    NumberOfScalarVariables   = i;
    NumberOfVectorVariables   = j;
    NumberOfParticleVariables = k;

  
    /*
     * check if scaling required...
     */
    if ( ElementModel && NodesChanged )
    {
      static vertex_t *vertex;
      static double s;

      vertex = CurrentObject->Geometry->Vertices;

      if ( !KeepScale )
      {
         xmin = ymin = zmin =  1.0e20;
         xmax = ymax = zmax = -1.0e20;

         for( i=0; i<ElementModel->NofNodes; i++ )
         {
           xmin = MIN( xmin,ElementModel->NodeArray[i] );
           ymin = MIN( ymin,ElementModel->NodeArray[ElementModel->NofNodes+i] );
           zmin = MIN( zmin,ElementModel->NodeArray[2*ElementModel->NofNodes+i] );

           xmax = MAX( xmax,ElementModel->NodeArray[i] );
           ymax = MAX( ymax,ElementModel->NodeArray[ElementModel->NofNodes+i] );
           zmax = MAX( zmax,ElementModel->NodeArray[2*ElementModel->NofNodes+i] );
         }
         s = MAX(MAX(xmax-xmin,ymax-ymin),zmax-zmin);

         CurrentObject->Geometry->Scale = s;
         CurrentObject->Geometry->MinMax[0].x[0] = xmin;
         CurrentObject->Geometry->MinMax[0].x[1] = ymin;
         CurrentObject->Geometry->MinMax[0].x[2] = zmin;

         CurrentObject->Geometry->MinMax[1].x[0] = xmax;
         CurrentObject->Geometry->MinMax[1].x[1] = ymax;
         CurrentObject->Geometry->MinMax[1].x[2] = zmax;

      } else {
         s = CurrentObject->Geometry->Scale;
         xmin = CurrentObject->Geometry->MinMax[0].x[0];
         ymin = CurrentObject->Geometry->MinMax[0].x[1];
         zmin = CurrentObject->Geometry->MinMax[0].x[2];

         xmax = CurrentObject->Geometry->MinMax[1].x[0];
         ymax = CurrentObject->Geometry->MinMax[1].x[1];
         zmax = CurrentObject->Geometry->MinMax[1].x[2];
      }

      for( i=0; i<CurrentObject->ElementModel->NofNodes; i++,vertex++ )
      {
        vertex->x[0] = (2.0 * (ElementModel->NodeArray[i] - xmin) - (xmax-xmin)) / s;
        vertex->x[1] = (2.0 * (ElementModel->NodeArray[ElementModel->NofNodes+i] - ymin) - (ymax-ymin)) / s;
        vertex->x[2] = (2.0 * (ElementModel->NodeArray[2*ElementModel->NofNodes+i] - zmin) - (zmax-zmin)) / s;
      }

      if ( NormalUpdate ) geo_vertex_normals( CurrentObject->Geometry,50.0 );
    }


    Tcl_SetResult( interp, result, TCL_STATIC );
    return TCL_OK;
}

#ifndef MINGW32
#include <GL/glx.h>
#else
#include <winuser.h>
#endif
#include <X11/Xlib.h>
#include <X11/Xutil.h>

static unsigned int FontBase;

Display *tkXDisplay()
{
    Display *ptr = NULL;
#ifndef MINGW32
    ptr = auxXDisplay();
#endif
    return ptr;
}

Window tkXWindow()
{
    Window ptr = 0;

#ifdef MINGW32
    ptr = auxGetHWND();
#else
    ptr = auxXWindow();
#endif

    return ptr;
}


static int ActivateGraphicsWindow( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
#ifdef MINGW32
    ShowWindow( tkXWindow(), SW_SHOWNORMAL );
#endif
    return TCL_OK;
}




static GLubyte rasterFont[][13] = {
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18},
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x36, 0x36, 0x36},
    {0x00, 0x00, 0x00, 0x66, 0x66, 0xff, 0x66, 0x66, 0xff, 0x66, 0x66, 0x00, 0x00},
    {0x00, 0x00, 0x18, 0x7e, 0xff, 0x1b, 0x1f, 0x7e, 0xf8, 0xd8, 0xff, 0x7e, 0x18},
    {0x00, 0x00, 0x0e, 0x1b, 0xdb, 0x6e, 0x30, 0x18, 0x0c, 0x76, 0xdb, 0xd8, 0x70},
    {0x00, 0x00, 0x7f, 0xc6, 0xcf, 0xd8, 0x70, 0x70, 0xd8, 0xcc, 0xcc, 0x6c, 0x38},
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x1c, 0x0c, 0x0e},

    {0x00, 0x00, 0x0c, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c},
    {0x00, 0x00, 0x30, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x30},
    {0x00, 0x00, 0x00, 0x00, 0x99, 0x5a, 0x3c, 0xff, 0x3c, 0x5a, 0x99, 0x00, 0x00},
    {0x00, 0x00, 0x00, 0x18, 0x18, 0x18, 0xff, 0xff, 0x18, 0x18, 0x18, 0x00, 0x00},
    {0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06, 0x03, 0x03},
    {0x00, 0x00, 0x3c, 0x66, 0xc3, 0xe3, 0xf3, 0xdb, 0xcf, 0xc7, 0xc3, 0x66, 0x3c},
    {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0xe7, 0x7e},
    {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0x07, 0x03, 0x03, 0xe7, 0x7e},
    {0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0xff, 0xcc, 0x6c, 0x3c, 0x1c, 0x0c},
    {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x03, 0x03, 0xff},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
    {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x03, 0x7f, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
    {0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x1c, 0x1c, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x06, 0x0c, 0x18, 0x30, 0x60, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06},
    {0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x06, 0x0c, 0x18, 0x30, 0x60},
    {0x00, 0x00, 0x18, 0x00, 0x00, 0x18, 0x18, 0x0c, 0x06, 0x03, 0xc3, 0xc3, 0x7e},
    {0x00, 0x00, 0x3f, 0x60, 0xcf, 0xdb, 0xd3, 0xdd, 0xc3, 0x7e, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18},
    {0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
    {0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e},
    {0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06},
    {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3},
    {0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e},
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
    {0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c},
    {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
    {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff},
    {0x00, 0x00, 0x3c, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x3c},
    {0x00, 0x03, 0x03, 0x06, 0x06, 0x0c, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x60, 0x60},
    {0x00, 0x00, 0x3c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x3c},
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18},
    {0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x38, 0x30, 0x70},
    {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0x7f, 0x03, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
    {0x00, 0x00, 0x7e, 0xc3, 0xc0, 0xc0, 0xc0, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03, 0x03},
    {0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x33, 0x1e},
    {0x7e, 0xc3, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0},
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x00, 0x18, 0x00},
    {0x38, 0x6c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x00, 0x00, 0x0c, 0x00},
    {0x00, 0x00, 0xc6, 0xcc, 0xf8, 0xf0, 0xd8, 0xcc, 0xc6, 0xc0, 0xc0, 0xc0, 0xc0},
    {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78},
    {0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xfc, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x7c, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x7c, 0x00, 0x00, 0x00, 0x00},
    {0xc0, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00},
    {0x03, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xfe, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x1c, 0x36, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x00},
    {0x00, 0x00, 0x7e, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xdb, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18, 0x3c, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00},
    {0xc0, 0x60, 0x60, 0x30, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0xff, 0x60, 0x30, 0x18, 0x0c, 0x06, 0xff, 0x00, 0x00, 0x00, 0x00},
    {0x00, 0x00, 0x0f, 0x18, 0x18, 0x18, 0x38, 0xf0, 0x38, 0x18, 0x18, 0x18, 0x0f},
    {0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18},
    {0x00, 0x00, 0xf0, 0x18, 0x18, 0x18, 0x1c, 0x0f, 0x1c, 0x18, 0x18, 0x18, 0xf0},
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x8f, 0xf1, 0x60, 0x00, 0x00, 0x00}
};

void MakeRasterFontDefault()
{
    GLuint i;

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    FontBase = glGenLists(128);
    for (i = 32; i < 127; i++) {
        glNewList( i+FontBase, GL_COMPILE );
        glBitmap( 8, 13, 0.0f, 2.0f, 10.0f, 0.0f, rasterFont[i-32] );
        glEndList();
    }
}


#ifndef MINGW32  

void InitializeXFonts()
{
    static char str[32], here = 0;
    int i,n;
    char **FontNames;

    FontNames = XListFonts( tkXDisplay(), "*", 10000, &n );
    
     for( i=0; i<n; i++ )
     {
       sprintf( str, "%d", i );
       Tcl_SetVar2( TCLInterp, "FontNames", str,FontNames[i], TCL_GLOBAL_ONLY );
     }
     
     sprintf( str, "%d", i );
     Tcl_SetVar( TCLInterp, "NumberOfFonts", str, TCL_GLOBAL_ONLY );
     
     MakeRasterFontDefault();
}


void MakeRasterFont( char *name )
{
    unsigned int first, last;

    XFontStruct *fontInfo;
    static char str[32], here = 0;

    CurrentXFont = fontInfo = XLoadQueryFont( tkXDisplay(),name );
    fprintf( stdout, "Font: [%s] %x\n", name, CurrentXFont );
    if ( fontInfo == NULL )
    {
        fprintf( stderr, "Can not find font: [%s]\n", name );
        return;
    }

    first = fontInfo->min_char_or_byte2;
    last  = fontInfo->max_char_or_byte2;
    ColorScaleFontSize = fontInfo->max_bounds.rbearing - fontInfo->min_bounds.lbearing;

    FontBase = glGenLists( last + 1 );
    glXUseXFont( fontInfo->fid,first,last-first+1,FontBase+first );
}
#endif 

void PrintString( char *s )
{
    double param[4],rgba[4];
    if ( GlobalOptions.OutputPS ) {
       glGetDoublev( GL_CURRENT_COLOR, rgba );
       glGetDoublev( GL_CURRENT_RASTER_POSITION, param );
       OutputPSString( param[0],param[1],ColorScaleFontSize, rgba[0],rgba[1],rgba[2], s );
    } else {
       glDisable( GL_LIGHTING );
       glDisable( GL_TEXTURE_1D );
       glPushAttrib( GL_LIST_BIT );
       glListBase( FontBase );
       glCallLists( strlen(s),GL_UNSIGNED_BYTE,(unsigned char *)s );
       glPopAttrib();
       glEnable( GL_TEXTURE_1D );
       glEnable( GL_LIGHTING );
    }
}

static int ShowString( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
    float x,y,z;

    if ( argc < 5 )
    {
        sprintf( interp->result, "string: Wrong number of arguments.\n" );
        return TCL_ERROR;
    }

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

    glDisable( GL_TEXTURE_1D );
    glDisable( GL_DEPTH_TEST );
    glDrawBuffer( GL_FRONT_AND_BACK );

    glColor3f( 1.0,1.0,1.0 );

    x = atof( argv[1] );
    y = atof( argv[2] );
    z = atof( argv[3] );

    glRasterPos3f( x,y,z );

    // if ( argc >= 6 ) MakeRasterFont( argv[5] );
    PrintString( argv[4] );

    glPopMatrix();

    epSwapBuffers();

    return TCL_OK;
}


static int SetFont( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
  int i;
    if ( argc < 1 )
    {
        sprintf( interp->result, "Set font: Wrong number of arguments.\n" );
        return TCL_ERROR;
    }

#ifdef MINGW32
    MakeRasterFontDefault();
#else
    MakeRasterFont( argv[1] );
#endif

    return TCL_OK;
}


int UpdateObject( ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static double L[1024],I[256];
    int i;

    element_model_t *model;

    char *str;

    model = CurrentObject->ElementModel;
    if ( !model ) return TCL_OK;

    vis_delete_visual( CurrentObject->VisualList );
    CurrentObject->VisualList = NULL;

    if ( ShowMeshLines )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Mesh" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        vis_set_param( VL, "ColorMap",          0,          0.0, MeshColorMap );
        vis_set_param( VL, "Material",  0, 0.0, MeshMaterial );
        vis_set_param( VL, "ColorData",         0,          0.0, NULL );
        vis_set_param( VL, "Style",      mesh_style_line,   0.0, NULL );
        vis_set_param( VL, "LineStyle",  line_style_line,   0.0, NULL );
        vis_set_param( VL, "LineQuality",      1,           0.0, NULL );
        vis_set_param( VL, "EdgeStyle",  edge_style_all,    0.0, NULL );
    }

    if  ( ShowParticles )
    {
        visual_t *VL = CurrentObject->VisualList;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Particles" );

        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &ParticleColor, (char *)Tcl_GetVar( TCLInterp,
              "ParticleColor",TCL_GLOBAL_ONLY ), model, 1, TRUE );
        GetVectorVariable( ParticleVelocity, (char *)Tcl_GetVar( TCLInterp,
              "ParticleVelocity",TCL_GLOBAL_ONLY ),model,  1 );
        GetParticleVariable( ParticleParticle, (char *)Tcl_GetVar( TCLInterp,
              "ParticleParticle",TCL_GLOBAL_ONLY ) );

        vis_set_param( VL, "Style", ParticleStyle, 0.0, NULL );

        vis_set_param( VL, "VectorData1", 0, 0.0, &ParticleVelocity[1] );
        vis_set_param( VL, "VectorData2", 0, 0.0, &ParticleVelocity[2] );
        vis_set_param( VL, "VectorData3", 0, 0.0, &ParticleVelocity[3] );

        vis_set_param( VL, "ParticleDataX", 0, 0.0, &ParticleParticle[0] );
        vis_set_param( VL, "ParticleDataY", 0, 0.0, &ParticleParticle[1] );
        vis_set_param( VL, "ParticleDataZ", 0, 0.0, &ParticleParticle[2] );
        vis_set_param( VL, "ParticleDataC", 0, 0.0, &ParticleParticle[3] );
        vis_set_param( VL, "ParticleDataI", 0, 0.0, &ParticleParticle[4] );

        vis_set_param( VL, "ColorData",   0, 0.0, &ParticleColor );

        vis_set_param( VL, "LineStyle",  ParticleLineStyle, 0.0, NULL );
        vis_set_param( VL, "LineQuality", ParticleQuality, 0.0, NULL );
        vis_set_param( VL, "LineWidth",    0, ParticleRadius, NULL );
        vis_set_param( VL, "ArrowStyle",   ParticleArrowStyle, 0.0, NULL );
        vis_set_param( VL, "OutDT", 0, ParticleOutDT, NULL );
        vis_set_param( VL, "MaxDT", 0, ParticleMaxDT, NULL );
        vis_set_param( VL, "Tolerance", 0, ParticleTolerance, NULL );
        vis_set_param( VL, "NofParticles", ParticleNofParticles, 0.0, NULL );
        vis_set_param( VL, "ColorMap", 0, 0.0, ParticleColorMap );
        vis_set_param( VL, "Material", 0, 0.0, ParticleMaterial );

        vis_set_param( VL, "Advance", ParticleAdvance, 0.0, NULL );
        vis_set_param( VL, "IntegMethod", ParticleIntegMethod, 0.0, NULL );
        vis_set_param( VL, "IntegPolicy", ParticleIntegPolicy, 0.0, NULL );
        ParticleAdvance = FALSE;
    }

    if  ( ShowVectors )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Arrows" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &VectorColor, (char *)Tcl_GetVar( TCLInterp,"VectorColor",TCL_GLOBAL_ONLY ),model, 0, TRUE );
        GetScalarVariable( &VectorLength, (char *)Tcl_GetVar( TCLInterp,"VectorLength",TCL_GLOBAL_ONLY ),model, 0, TRUE );
        GetScalarVariable( &VectorThreshold, (char *)Tcl_GetVar( TCLInterp,"VectorThreshold",TCL_GLOBAL_ONLY ),model, 0, TRUE );
        GetVectorVariable( VectorArrow, (char *)Tcl_GetVar( TCLInterp,"VectorArrow",TCL_GLOBAL_ONLY ),model, 0 );

        vis_set_param( VL, "VectorData0", 0, 0.0, &VectorLength );
        vis_set_param( VL, "VectorData1", 0, 0.0, &VectorArrow[1] );
        vis_set_param( VL, "VectorData2", 0, 0.0, &VectorArrow[2] );
        vis_set_param( VL, "VectorData3", 0, 0.0, &VectorArrow[3] );


        vis_set_param( VL, "ColorData",   0, 0.0, &VectorColor );
        vis_set_param( VL, "LengthData",  0, 0.0, &VectorLength );

        vis_set_param( VL, "ThresholdData",  0, 0.0, &VectorThreshold );
        vis_set_param( VL, "Ceiling",  0, VectorCeiling,  NULL );
        vis_set_param( VL, "Floor",    0, VectorFloor, NULL );

        vis_set_param( VL, "LineStyle",  VectorLineStyle, 0.0, NULL );

        vis_set_param( VL, "LengthScale",  0, VectorLengthScale, NULL );

        vis_set_param( VL, "EqualLength", FALSE, 0.0, NULL );

        vis_set_param( VL, "LineQuality", VectorQuality, 0.0, NULL );
        vis_set_param( VL, "ColorMap", 0, 0.0, ArrowColorMap );
        vis_set_param( VL, "Material", 0, 0.0, ArrowMaterial );
    }

    if  ( ShowContours )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Contour Lines" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &ContourColor[CurrentObject->Id], (char *)Tcl_GetVar( TCLInterp,"ContourColor",
                   TCL_GLOBAL_ONLY ),model, 0, !ContourColorSetMinMax );

        vis_set_param( VL, "ColorData",       0,            0.0, &ContourColor[CurrentObject->Id] );
        vis_set_param( VL, "ContourData",     0,            0.0, &ContourColor[CurrentObject->Id] );
        vis_set_param( VL, "LineWidth",       0,   ContourRadius,     NULL );
        vis_set_param( VL, "LineStyle",   ContourLineStyle, 0.0,     NULL );
        vis_set_param( VL, "LineQuality", ContourQuality,   0.0,     NULL );

        for( i=0; i<ContourLines; i++ )
        {
            static char name[32];

            sprintf( name, "%d", i );
            str = (char *)Tcl_GetVar2( TCLInterp, "ContourValues",name, TCL_GLOBAL_ONLY );
            if ( str ) sscanf( str,  "%lf", &L[i] );
        }

        vis_set_param( VL, "Levels",0,0.0,L );
        vis_set_param( VL, "NofLevels",ContourLines,0.0,NULL );
        vis_set_param( VL, "ColorMap", 0, 0.0, ContourColorMap );
        vis_set_param( VL, "Material", 0, 0.0, ContourMaterial );
    }

    if  ( ShowIsosurfaces )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Isosurfaces" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &IsosurfaceColor, (char *)Tcl_GetVar( TCLInterp,"IsosurfaceColor",
                   TCL_GLOBAL_ONLY ),model, 0, !IsosurfaceColorSetMinMax );

        GetScalarVariable( &IsosurfaceContour, (char *)Tcl_GetVar( TCLInterp,"IsosurfaceContour",
                 TCL_GLOBAL_ONLY ),model, 0, !IsosurfaceContourSetMinMax );

        GetVectorVariable( IsosurfaceNormal, (char *)Tcl_GetVar( TCLInterp,
                 "IsosurfaceNormal",TCL_GLOBAL_ONLY ),model, 0 );
        vis_set_param( VL, "Style",   IsosurfaceStyle+1, 0.0,     NULL );
        vis_set_param( VL, "ColorData",       0,            0.0, &IsosurfaceColor );
        vis_set_param( VL, "ContourData",     0,            0.0, &IsosurfaceContour );
        vis_set_param( VL, "NormalData0",     0,            0.0, &IsosurfaceNormal[1] );
        vis_set_param( VL, "NormalData1",     0,            0.0, &IsosurfaceNormal[2] );
        vis_set_param( VL, "NormalData2",     0,            0.0, &IsosurfaceNormal[3] );
        vis_set_param( VL, "LineWidth",       0,   IsosurfaceRadius,    NULL );
        vis_set_param( VL, "LineStyle",   IsosurfaceLineStyle,  0.0,    NULL );
        vis_set_param( VL, "LineQuality", IsosurfaceQuality,    0.0,    NULL );

        vis_set_param( VL, "Recompute",   IsosurfaceRecompute,  0.0,    NULL );
        IsosurfaceRecompute = FALSE;

        for( i=0; i<IsosurfaceContours; i++ )
        {
            static char name[32];

            sprintf( name, "%d", i );
            sscanf( Tcl_GetVar2( TCLInterp, "IsosurfaceValues",name, TCL_GLOBAL_ONLY ), "%lf", &I[i] );
        }

        vis_set_param( VL, "Levels",0,0.0,I );
        vis_set_param( VL, "NofLevels",IsosurfaceContours,0.0,NULL );
        vis_set_param( VL, "ColorMap",  0, 0.0, IsoSurfaceColorMap );
        vis_set_param( VL, "Material", 0, 0.0,  IsoSurfaceMaterial );
    }

    if ( ShowSpheres )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Spheres" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &SphereColor,(char *)Tcl_GetVar(TCLInterp,"SphereColor",TCL_GLOBAL_ONLY),model, 0, TRUE );
        GetScalarVariable( &SphereRadius,(char *)Tcl_GetVar(TCLInterp,"SphereRadius",TCL_GLOBAL_ONLY),model, 0, TRUE );
        GetScalarVariable( &SphereThreshold,(char *)Tcl_GetVar(TCLInterp,"SphereThreshold",TCL_GLOBAL_ONLY),model, 0, TRUE );

        vis_set_param( VL, "ColorMap",          0,          0.0, SphereColorMap );
        vis_set_param( VL, "Material",          0,          0.0, SphereMaterial );

        vis_set_param( VL, "ColorData",         0,      0.0, &SphereColor );
        vis_set_param( VL, "RadiusData",        0,      0.0, &SphereRadius );
        vis_set_param( VL, "RadiusScale",       0,      SphereRadiusScale, NULL );
        vis_set_param( VL, "ThresholdData",     0,      0.0, &SphereThreshold );
        vis_set_param( VL, "Floor",             0,      SphereFloor, NULL );
        vis_set_param( VL, "Ceiling",           0,      SphereCeiling, NULL );
        vis_set_param( VL, "Quality", SphereQuality,    0.0, NULL );
    }

    if ( ShowColorMesh )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "Mesh" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &MeshColor[CurrentObject->Id], (char *)Tcl_GetVar( TCLInterp,"MeshColor",TCL_GLOBAL_ONLY ),
                                     model, 0, !MeshColorSetMinMax );

        vis_set_param( VL, "ColorData",         0,        0.0, &MeshColor[CurrentObject->Id] );
        vis_set_param( VL, "Style",       MeshStyle+1,    0.0, NULL );
        vis_set_param( VL, "LineStyle",   MeshLineStyle,  0.0, NULL );
        vis_set_param( VL, "LineQuality", MeshQuality,    0.0, NULL );
        vis_set_param( VL, "LineWidth",   0,    MeshRadius, NULL );
        vis_set_param( VL, "EdgeStyle",   MeshEdgeStyle,  0.0, NULL );
        vis_set_param( VL, "NodeNumbers", MeshNodeNumbers,0.0, NULL );
        vis_set_param( VL, "ColorMap",  0, 0.0, MeshColorMap );
        vis_set_param( VL, "Material",  0, 0.0, MeshMaterial );
        vis_set_param( VL, "EdgeMaterial",  0, 0.0, &DefaultEdgeMaterial );
    }

    if ( ShowColorScale )
    {
        visual_t *VL = NULL;

        if ( !VL )
        {
            VL = (visual_t *)vis_new_visual( "ColorScale" );
        } else VL->Next = NULL;
        CurrentObject->VisualList = (visual_t *)vis_link_visual( CurrentObject->VisualList,VL );

        GetScalarVariable( &ColorScaleColor, (char *)Tcl_GetVar( TCLInterp,"ColorScaleColor",TCL_GLOBAL_ONLY ),
                                model, 0, !ColorScaleColorSetMinMax );

        vis_set_param( VL, "ColorData",         0,        0.0, &ColorScaleColor );
        vis_set_param( VL, "ColorMap",  0, 0.0, &DefaultColorMap );
        vis_set_param( VL, "Material",  0, 0.0, &DefaultMaterial );
        vis_set_param( VL, "XPosition", 0, ColorScaleX, NULL );
        vis_set_param( VL, "YPosition", 0, ColorScaleY, NULL );
        vis_set_param( VL, "Length",    0, ColorScaleLength, NULL );
        vis_set_param( VL, "Thickness", 0, ColorScaleThickness, NULL );
        vis_set_param( VL, "Style",     ColorScaleStyle, 0.0, NULL );
        vis_set_param( VL, "Font Color",ColorScaleFontColor, 0.0, NULL );
        vis_set_param( VL, "Font Size", 0, ColorScaleFontSize, NULL );
        vis_set_param( VL, "Entries",   ColorScaleEntries, 0.0, NULL );
        vis_set_param( VL, "Decimals",  ColorScaleDecimals, 0.0, NULL );
    }

    return TCL_OK;
}

void opengl_draw()
{
    Tcl_Eval( TCLInterp, ".buts.play configure -back red -relief sunken; update" );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glClearColor( br,bg,bb,1.0 );
    if ( GraphicsClearOn )
    {
        glDrawBuffer( GL_BACK );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
    } else {
        glDrawBuffer( GL_FRONT_AND_BACK );
    }

#if 0
    if ( !CurrentObject->Geometry || CurrentObject->Geometry->VertexCount <= 0 ) {
        Tcl_Eval( TCLInterp, ".buts.play configure -back green -relief raised; update" );
        epSwapBuffers(); return;
    }
#endif

    if ( epMouseDown && epMouseDownTakesTooLong )
    {
        glDisable( GL_BLEND );
        glDisable( GL_DEPTH_TEST );
    } else
    {
#if 0
        glEnable( GL_BLEND );
#endif
        glEnable( GL_DEPTH_TEST );
    }

    if ( cam_display_list( Camera,&VisualObject ) ) epSwapBuffers();

    Tcl_Eval( TCLInterp, ".buts.play configure -back green -relief raised; update" );
}


void Reshape( GLsizei x,GLsizei y)
{
    GraphicsXSize  = x;
    GraphicsYSize  = y;
    GraphicsAspect = (double)x/(double)y;

    DrawItSomeTimeWhenIdle();
}

void CompRot( matrix_t M ) {
  static double bx, by, bz;

  bx = M[0][0]*ax + M[1][0]*ay + M[2][0]*az;
  by = M[0][1]*ax + M[1][1]*ay + M[2][1]*az;
  bz = M[0][2]*ax + M[1][2]*ay + M[2][2]*az;

  ax = bx;
  ay = by;
  az = bz;
}

void epMouseDownProc(int Xpos, int Ypos)
{
    int x,y,x_root,y_root;
    Window root,child;
    static GLint viewport[4];
    transform_t *transform;
    double scale;

#ifndef MINGW32
    XQueryPointer( tkXDisplay(),tkXWindow(),&root,&child,&x_root,&y_root,&x,&y,&epMouseDown );
    epMouseDown &= (ShiftMask | Button1Mask | Button2Mask);
#else
    auxGetMouseLoc( &x, &y );

    epMouseDown = 0;
    if ( GetAsyncKeyState( VK_LBUTTON ) & 0x8000 ) epMouseDown |= Button1Mask;
    if ( GetAsyncKeyState( VK_MBUTTON ) & 0x8000 ) epMouseDown |= Button2Mask;
    if ( GetAsyncKeyState( VK_RBUTTON ) & 0x8000 ) epMouseDown |= Button2Mask;
    if ( GetAsyncKeyState(  VK_SHIFT  ) & 0x8000 ) epMouseDown |= ShiftMask;
#endif

    epMouseDownTakesTooLong = FALSE;

    while( epMouseDown )
    {
        if ( epMouseDown & Button2Mask )
        {
            ax = ay = az = 0.0;
            sx = sy = sz = 0.0;

            if ( epMouseDown & Button1Mask )
            {
                if ( ABS(Ypos-y) > ABS(Xpos-x) )
                {
                    sx = 0.025*(Ypos-y);
                    sy = 0.025*(Ypos-y);
                    sz = 0.025*(Ypos-y);
                } else
                {
                    sx = 0.025*(x-Xpos);
                    sy = 0.025*(x-Xpos);
                    sz = 0.025*(x-Xpos);
                }

                obj_scale( CurrentObject,sx,sy,sz,'a',TRUE );
            } else
            {
	        // scale = 0.4;
	        glGetIntegerv(GL_VIEWPORT, viewport);
	        transform = &CurrentObject->Transform;
  	        scale=180.0/(double)(viewport[3]+1)/transform->SclZ;

                if ( epMouseDown & ShiftMask )
                {
                    if ( ABS(Ypos-y) > ABS(Xpos-x) )
                        az = scale*(y-Ypos);
                    else
                        az = scale*(Xpos-x);
                } else {
                    if ( ABS(y-Ypos) > ABS(x-Xpos) )
                        ax = scale*(y-Ypos);
                    else
                        ay = scale*(x-Xpos);

		    // Compensate for rotation matrix:
		    //--------------------------------
		    CompRot( transform->RotMatrix );

                }

                obj_rotate( CurrentObject,ax,ay,az,'a',TRUE );
            }
        }
        else if ( epMouseDown & Button1Mask )
        {
            tx = ty = tz = 0;
            if ( epMouseDown & ShiftMask )
            {
                if ( ABS(Ypos-y) > ABS(Xpos-x) )
                    tz = 0.03*(Ypos-y);
                else
                    tz = 0.03*(x-Xpos);
            } else {
	        // scale = 0.01;
  	        glGetIntegerv(GL_VIEWPORT, viewport);
	        transform = &CurrentObject->Transform;
	        scale = 2.0/(double)(viewport[3]+1)/transform->SclZ;

                tx = scale*(x-Xpos);
                ty = scale*(Ypos-y);
            }
            obj_translate( CurrentObject,tx,ty,tz,'a',TRUE );
        }

        Xpos = x;
        Ypos = y;

        opengl_draw();

#ifndef MINGW32
        XQueryPointer( tkXDisplay(),tkXWindow(),&root,&child,&x_root,&y_root,&x,&y,&epMouseDown );
	epMouseDown &= (ShiftMask | Button1Mask | Button2Mask);
#else
        epMouseDown = 0;
        if ( GetAsyncKeyState( VK_LBUTTON ) & 0x8000 ) epMouseDown |= Button1Mask;
        if ( GetAsyncKeyState( VK_MBUTTON ) & 0x8000 ) epMouseDown |= Button2Mask;
        if ( GetAsyncKeyState( VK_RBUTTON ) & 0x8000 ) epMouseDown |= Button2Mask;
        if ( GetAsyncKeyState(  VK_SHIFT  ) & 0x8000 ) epMouseDown |= ShiftMask;

        auxGetMouseLoc( &x, &y );
#endif
    }

    epMouseDown = FALSE;
    if ( epMouseDownTakesTooLong ) { opengl_draw(); }
}

#ifndef MINGW32
void Mouse( AUX_EVENTREC *event )
{
    int MouseXPosition =  event->data[AUX_MOUSEX];
    int MouseYPosition =  event->data[AUX_MOUSEY];
    int MouseStatus    =  event->data[AUX_MOUSESTATUS];

    if ( event->event == AUX_MOUSEDOWN )
    {
        if ( MouseXPosition >= 0 && MouseXPosition < GraphicsXSize &&
            MouseYPosition >= 0 && MouseYPosition < GraphicsYSize )
        {
            epMouseDownProc(  MouseXPosition, MouseYPosition );
        }
    }
}
#else
void Mouse( )
{
    int MouseXPosition;
    int MouseYPosition;


    if ( GetFocus() != tkXWindow() ) return;


    auxGetMouseLoc( &MouseXPosition, &MouseYPosition );

    if ( MouseXPosition >= 0 && MouseXPosition < GraphicsXSize &&
         MouseYPosition >= 0 && MouseYPosition < GraphicsYSize )
    {
        epMouseDownProc(  MouseXPosition, MouseYPosition );
    }
}
#endif

static int UpdateDisplay(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    int x,y,x_root,y_root;
    Window root,child;


#ifndef MINGW32
    XQueryPointer( tkXDisplay(),tkXWindow(),&root,&child,&x_root,&y_root,&x,&y,&epMouseDown );
    epMouseDown &= (ShiftMask | Button1Mask | Button2Mask);
#else
    auxGetMouseLoc( &x, &y );

    epMouseDown = 0;
    if ( GetAsyncKeyState( VK_LBUTTON ) & 0x8000 ) epMouseDown |= Button1Mask;
    if ( GetAsyncKeyState( VK_MBUTTON ) & 0x8000 ) epMouseDown |= Button2Mask;
    if ( GetAsyncKeyState( VK_RBUTTON ) & 0x8000 ) epMouseDown |= Button2Mask;
    if ( GetAsyncKeyState(  VK_SHIFT  ) & 0x8000 ) epMouseDown |= ShiftMask;
#endif

    if ( epMouseDown )  epMouseDownProc( x,y );

#ifdef USE_TK
    tkExec( 0 );
#endif
    return TCL_OK;
}


void TestTkEvent()
{
    BreakLoop = FALSE;

    while( Tcl_DoOneEvent(TCL_DONT_WAIT) );
#ifdef MINGW32
    if ( GetKeyState( VK_LBUTTON ) & 0x8000 ) Mouse();
    if ( GetKeyState( VK_MBUTTON ) & 0x8000 ) Mouse();
    if ( GetKeyState( VK_RBUTTON ) & 0x8000 ) Mouse();
#endif
    Tk_Sleep(100);
}

static int WindowSize( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   unsigned int width, height;
   int dx, dy, ox, oy, nx, ny, viewp[4];


   if ( argc < 3 ) {
      sprintf( interp->result, "Usage: winsize width height" );
      return TCL_ERROR;
   }

   width  = atoi( *++argv );
   height = atoi( *++argv );

#ifdef MINGW32
   SetWindowPos( auxGetHWND(), HWND_TOP, 0, 0, width, height, SWP_NOMOVE | SWP_NOACTIVATE );

   // Measure the viewport size and make corrections winsize if necessary:
   //---------------------------------------------------------------------
   opengl_draw();
   glGetIntegerv( GL_VIEWPORT, viewp );

   ox = viewp[0];
   oy = viewp[1];
   nx = viewp[2]+1;
   ny = viewp[3]+1;

   dx = width - nx;
   dy = height -ny;

   width += dx;
   height += dy;

   SetWindowPos( auxGetHWND(), HWND_TOP, 0, 0, width, height, SWP_NOMOVE | SWP_NOACTIVATE );

#else
   XResizeWindow( tkXDisplay(), tkXWindow(), width, height );
#endif

   return TCL_OK;
}

static int WindowPosition( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   if ( argc < 3 ) {
      sprintf( interp->result, "Usage: winpos xpos ypos" );
      return TCL_ERROR;
   }

   int ox = atoi( *++argv );
   int oy = atoi( *++argv );

#ifdef MINGW32
   SetWindowPos( auxGetHWND(), HWND_TOP, ox, oy, 0, 0, SWP_NOSIZE | SWP_NOACTIVATE );
#else
   XMoveWindow( (Display *)tkXDisplay(), tkXWindow(), ox, oy );
#endif

   return TCL_OK;
}


static int MPlayer( ClientData cl, Tcl_Interp *interp, int argc, char **argv )
{
#if defined(MINGW32)

  return TCL_OK;

#else

  if( argc < 2 ) {
    sprintf( interp->result, "Usage: mplayer filename");
    return TCL_ERROR;
  }

  // File name:
  char *fileName = *++argv;

  // Get Window Id:
  GLXDrawable drawable = glXGetCurrentDrawable();
  int winId = (int)drawable;

  // Player command:
  char playCmd[1024];
  sprintf( playCmd, "mplayer -wid %d %s", winId, fileName );

  // Call mplayer:
  system( playCmd );

  return TCL_OK;

#endif
}


int main(int argc,char **argv)
{
  static char init[1024],initcommands[1024],tmp[1024],ephome[512];
  int i,size[4];

  /* For MinGW */
  static char szAppPath[512] = "";
  static char szAppDirectory[512] = "";
  char *exeName;
    
    if((argc > 1) && (!strcmp(argv[1], "-v"))) {
      fprintf(stdout, "ElmerPost v.5.4\n");
      return 0;
    }
    
    if ( getenv("ELMER_POST_HOME") == NULL )
    {
      /* use default installation directory just if nothing is set */
#if defined(MINGW32)
      GetModuleFileName(NULL, szAppPath, 512);
      exeName = strrchr(szAppPath, '\\');
      i = (int)(exeName-szAppPath);
      if(i < 0) i = 0;
      if(i > 512) i = 512;
      strncpy(szAppDirectory, szAppPath, i);

      _snprintf(ephome, 512,
		"ELMER_POST_HOME=%s\\..\\share\\elmerpost",
		szAppDirectory);

      printf("%s\n", ephome);
#else
      snprintf( ephome, 512, "ELMER_POST_HOME=%s", ELMER_POST_HOME );
#endif

      putenv( ephome );
    }


    Tcl_FindExecutable( *argv++ );
    TCLInterp = Tcl_CreateInterp();

    auxInitDisplayMode( AUX_DOUBLE | AUX_RGB | AUX_DEPTH | AUX_DIRECT );
    auxInitPosition( 0, 0, 500, 500 );
    auxInitWindow( "ELMER POST GRAPHICS" );

    strcpy( initcommands, "" );

    while ( *argv ) {

      if ( strcmp( *argv, "-id" ) == 0 ) {
        if ( *++argv )
        {
          Tcl_SetVar( TCLInterp, "elmerpost_id", *argv++,
                          TCL_GLOBAL_ONLY );
        }

      } else if ( strcmp( *argv, "-file" ) == 0 ) {
        if ( *++argv )
        {
          strcat( initcommands, "readfile " );

          /* Quote filename if not already quoted */
          if ( (*argv)[strlen(*argv) - 1] != '\"' )
             strcat( initcommands, "\"" );

          strcat( initcommands, *argv );

          if ( (*argv)[strlen(*argv) - 1] != '\"' )
             strcat( initcommands, "\"" );

          strcat( initcommands, " " );
          ++argv;
        }

      } else if ( 1020 - strlen(initcommands) > strlen( *argv ) ) {
        strcat( initcommands, *argv++ );
        strcat( initcommands, " " );

      } else break;
    }

    Tcl_SetVar( TCLInterp, "argv", "-name \"ELMER POST PROCESSING\"",
                        TCL_GLOBAL_ONLY );

    Tcl_Init( TCLInterp );
    Tk_Init( TCLInterp );

    Tcl_StaticPackage(TCLInterp, "Tk", Tk_Init, Tk_SafeInit);

    Tcl_CreateCommand( TCLInterp, "GetInterpolate",  (Tcl_CmdProc *)GetInterpolate,
                 (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "GetColorMap",     (Tcl_CmdProc *)GetColorMap,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "StopProcessing",  (Tcl_CmdProc *)StopProcessing,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "UpdateBackColor", (Tcl_CmdProc *)UpdateBackColor,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "UpdateColor",     (Tcl_CmdProc *)UpdateColor,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "UpdateEdgeColor", (Tcl_CmdProc *)UpdateEdgeColor,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "UpdateVariable",  (Tcl_CmdProc *)UpdateVariable,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "c_TimeStep",      (Tcl_CmdProc *)TimeStep,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "c_LoadCamera",    (Tcl_CmdProc *)LoadCamera,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "c_CurrentCamera", (Tcl_CmdProc *)CurrentCamera,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "c_SetCamera",     (Tcl_CmdProc *)SetCamera,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "UpdateDisplay",   (Tcl_CmdProc *)UpdateDisplay,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "c_MathCommand",   (Tcl_CmdProc *)MathCommand,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "UpdateObject",    (Tcl_CmdProc *)UpdateObject,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "group",          (Tcl_CmdProc *)GroupDisplay,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "setfont",        (Tcl_CmdProc *)SetFont,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "object", (Tcl_CmdProc *)SetObject,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "parent", (Tcl_CmdProc *)SetParentObject,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "normals", (Tcl_CmdProc *)RecomputeNormals,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "clip", (Tcl_CmdProc *)ClipPlane,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "ActivateGraphicsWindow", (Tcl_CmdProc *)ActivateGraphicsWindow,
                  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "winsize", (Tcl_CmdProc *)WindowSize,
		       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
    
    Tcl_CreateCommand( TCLInterp, "winpos", (Tcl_CmdProc *)WindowPosition, 
		       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

#if defined(HAVE_FTGL_NEW) || defined(HAVE_FTGL_OLD)
    Tcl_CreateCommand( TCLInterp, "fttext", (Tcl_CmdProc *)FtText, 
		       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    Tcl_CreateCommand( TCLInterp, "ftfont", (Tcl_CmdProc *)FtFont, 
		       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
#endif

    Tcl_CreateCommand( TCLInterp, "mplayer", (Tcl_CmdProc *)MPlayer, 
		       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

    CurrentObject = &VisualObject;
    CurrentObject->Name = strcpy( malloc(strlen("default")+1), "default" );

    obj_object_initialize( CurrentObject );

    Tcl_LinkVar( TCLInterp, "BreakLoop", (char *)&BreakLoop, TCL_LINK_INT );

    if ( !(Camera = (camera_t *)cam_load_cameras( Camera,NULL ) ) )
    {
        fprintf( stderr, "ElmerPost: Can't initialize default camera. "
            "Something is definitely wrong here...\n" );
        exit( 0 );
    }

    if ( !vis_initialize_visual_types() )
    {
        fprintf( stderr, "ElmerPost: Can't initialize visual types. "
            "Something is definitely wrong here...\n" );
        exit( 0 );
    }

    if ( !elm_initialize_element_types() )
    {
        fprintf( stderr, "ElmerPost: Can't initialize element types. "
            "Something is definitely wrong here...\n" );
        exit( 0 );
    }

    mtc_init( NULL, stdout, stderr );

    Tcl_LinkVar( TCLInterp, "NumberOfTimesteps", (char *)&ElementModel.NofTimesteps, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "NumberOfScalarVariables", (char *)&NumberOfScalarVariables, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "NumberOfVectorVariables", (char *)&NumberOfVectorVariables, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "NumberOfParticleVariables", (char *)&NumberOfParticleVariables, TCL_LINK_INT );

    strcpy( tmp , "DisplayStyle(Vectors)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowVectors, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(Contours)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowContours, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(MeshLines)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowMeshLines, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(ColorMesh)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowColorMesh, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(Spheres)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowSpheres, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(Isosurfaces)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowIsosurfaces, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(Particles)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowParticles, TCL_LINK_INT );
    strcpy( tmp , "DisplayStyle(ColorScale)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&ShowColorScale, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "MeshStyle", (char *)&MeshStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "MeshLineStyle", (char *)&MeshLineStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "MeshEdgeStyle", (char *)&MeshEdgeStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "MeshQuality", (char *)&MeshQuality, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "MeshRadius", (char *)&MeshRadius, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "MeshNodeNumbers", (char *)&MeshNodeNumbers, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "MeshColorMin", (char *)&MeshColor[0].min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "MeshColorMax", (char *)&MeshColor[0].max, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "MeshColorSetMinMax", (char *)&MeshColorSetMinMax, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "ContourLines", (char *)&ContourLines, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ContourLineStyle", (char *)&ContourLineStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ContourQuality", (char *)&ContourQuality, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ContourRadius", (char *)&ContourRadius, TCL_LINK_DOUBLE );

    Tcl_LinkVar( TCLInterp, "ContourColorMin", (char *)&ContourColor[0].min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ContourColorMax", (char *)&ContourColor[0].max, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ContourColorSetMinMax", (char *)&ContourColorSetMinMax, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "VectorLineStyle", (char *)&VectorLineStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "VectorQuality", (char *)&VectorQuality, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "VectorRadius", (char *)&VectorRadius, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "VectorLengthScale", (char *)&VectorLengthScale, TCL_LINK_DOUBLE );

    Tcl_LinkVar( TCLInterp, "VectorCeiling", (char *)&VectorCeiling, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "VectorFloor", (char *)&VectorFloor, TCL_LINK_DOUBLE );

    Tcl_LinkVar( TCLInterp, "VectorThresholdMin", (char *)&VectorThreshold.min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "VectorThresholdMax", (char *)&VectorThreshold.max, TCL_LINK_DOUBLE );

    Tcl_LinkVar( TCLInterp, "IsosurfaceStyle", (char *)&IsosurfaceStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "IsosurfaceContours", (char *)&IsosurfaceContours, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "IsosurfaceLineStyle", (char *)&IsosurfaceLineStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "IsosurfaceQuality", (char *)&IsosurfaceQuality, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "IsosurfaceRadius", (char *)&IsosurfaceRadius, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "IsosurfaceRecompute", (char *)&IsosurfaceRecompute, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "IsosurfaceColorMin", (char *)&IsosurfaceColor.min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "IsosurfaceColorMax", (char *)&IsosurfaceColor.max, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "IsosurfaceColorSetMinMax", (char *)&IsosurfaceColorSetMinMax, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "IsosurfaceContourSetMinMax", (char *)&IsosurfaceContourSetMinMax, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "IsosurfaceContourMin", (char *)&IsosurfaceContour.min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "IsosurfaceContourMax", (char *)&IsosurfaceContour.max, TCL_LINK_DOUBLE );

    Tcl_LinkVar( TCLInterp, "ParticleStyle", (char *)&ParticleStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleLineStyle", (char *)&ParticleLineStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleArrowStyle", (char *)&ParticleArrowStyle, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleQuality", (char *)&ParticleQuality, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleRadius", (char *)&ParticleRadius, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ParticleOutDT", (char *)&ParticleOutDT, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ParticleMaxDT", (char *)&ParticleMaxDT, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ParticleTolerance", (char *)&ParticleTolerance, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ParticleNofParticles", (char *)&ParticleNofParticles, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleAdvance", (char *)&ParticleAdvance, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleIntegMethod", (char *)&ParticleIntegMethod, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ParticleIntegPolicy", (char *)&ParticleIntegPolicy, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "SphereQuality", (char *)&SphereQuality, TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "SphereRadiusScale", (char *)&SphereRadiusScale, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "SphereFloor", (char *)&SphereFloor, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "SphereCeiling", (char *)&SphereCeiling, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "SphereThresholdMin", (char *)&SphereThreshold.min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "SphereThresholdMax", (char *)&SphereThreshold.max, TCL_LINK_DOUBLE );

    Tcl_LinkVar( TCLInterp, "ColorScaleStyle",     (char *)&ColorScaleStyle,     TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ColorScaleEntries",   (char *)&ColorScaleEntries,   TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ColorScaleDecimals",   (char *)&ColorScaleDecimals,   TCL_LINK_INT );
    Tcl_LinkVar( TCLInterp, "ColorScaleX",         (char *)&ColorScaleX,         TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ColorScaleY",         (char *)&ColorScaleY,         TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ColorScaleThickness", (char *)&ColorScaleThickness, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ColorScaleLength",    (char *)&ColorScaleLength,    TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ColorScaleFontColor", (char *)&ColorScaleFontColor, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "ColorScaleColorMin", (char *)&ColorScaleColor.min, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ColorScaleColorMax", (char *)&ColorScaleColor.max, TCL_LINK_DOUBLE );
    Tcl_LinkVar( TCLInterp, "ColorScaleColorSetMinMax", (char *)&ColorScaleColorSetMinMax, TCL_LINK_INT );


    Tcl_LinkVar( TCLInterp, "GraphicsClearOn", (char *)&GraphicsClearOn, TCL_LINK_INT );


    Tcl_LinkVar( TCLInterp, "KeepScale",    (char *)&KeepScale, TCL_LINK_INT );

    Tcl_LinkVar( TCLInterp, "NormalUpdate", (char *)&NormalUpdate, TCL_LINK_INT );

    strcpy( tmp , "GlobalOptions(SurfaceSides)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.SurfaceSides, TCL_LINK_INT );

    strcpy( tmp , "GlobalOptions(VolumeSides)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.VolumeSides, TCL_LINK_INT );

    strcpy( tmp , "GlobalOptions(VolumeEdges)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.VolumeEdges, TCL_LINK_INT );

    strcpy( tmp , "GlobalOptions(StereoMode)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.StereoMode, TCL_LINK_INT );
    strcpy( tmp , "GlobalOptions(StereoTran)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.StereoTran, TCL_LINK_DOUBLE );
    strcpy( tmp , "GlobalOptions(StereoRot)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.StereoRot, TCL_LINK_DOUBLE );

    strcpy( tmp , "GlobalOptions(OutputPS)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.OutputPS, TCL_LINK_INT );
    strcpy( tmp , "GlobalOptions(FitToPagePS)" );
    Tcl_LinkVar( TCLInterp, tmp, (char *)&GlobalOptions.FitToPagePS, TCL_LINK_INT );

    Misc_Init( TCLInterp );
    Transforms_Init( TCLInterp );
    Readfile_Init( TCLInterp );
    Matctcl_Init( TCLInterp );


    signal( SIGFPE, SIG_IGN );
    signal( SIGINT, int_sig );


    *init = '\0';
    if ( getenv("ELMER_POST_HOME") )
    {
        strncat( init,getenv("ELMER_POST_HOME"),511);
        strncat( init,"/",511 );
    }
    strncat( init,"tcl/init.tcl",511 );
    fprintf( stdout, "Initialization File: [%s]\n", init );
    fflush(stdout);
    Tcl_EvalFile( TCLInterp,init );

    while( Tk_DoOneEvent(TCL_DONT_WAIT) );

#ifdef MINGW32
    auxReshapeFunc( (AUXRESHAPEPROC)Reshape );
    auxExposeFunc( (AUXEXPOSEPROC)Reshape );
    auxIdleFunc( (AUXIDLEPROC)TestTkEvent );
#else
    auxReshapeFunc( Reshape );
    auxExposeFunc( Reshape );
    auxIdleFunc( TestTkEvent );
    auxMouseFunc( AUX_LEFTBUTTON,   AUX_MOUSEDOWN, Mouse );
    auxMouseFunc( AUX_MIDDLEBUTTON, AUX_MOUSEDOWN, Mouse );
#endif


    gra_init();

#ifndef MINGW32
    InitializeXFonts();
#else
    MakeRasterFontDefault();
#endif

    {
      Tcl_DString dstring;
      char *buf;

      buf = Tcl_ExternalToUtfDString( NULL, initcommands,strlen(initcommands),&dstring);
      Tcl_Eval( TCLInterp, buf );
      Tcl_DStringFree( &dstring );
    }

    if ( getenv("ELMER_POST_INIT") )
    {
        *init = '\0';
        strncat( init,getenv("ELMER_POST_INIT"),511);
        fprintf( stdout, "User initialization file: [%s]\n", init );
        fflush(stdout);
        Tcl_EvalFile( TCLInterp,init );
        while( Tk_DoOneEvent(TCL_DONT_WAIT) );
    }
    
#ifdef MINGW32
    auxMainLoop( (AUXMAINPROC)DrawItSomeTimeWhenIdle );
#else
    auxMainLoop(TestTkEvent);
#endif
}
