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
 *     File reading + MATC model setting command
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
 *                       Date: 27 Sep 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/


#include "../elmerpost.h"
#include <tcl.h>

extern double XMin,XMax,YMin,YMax,ZMin,ZMax;
extern double xmin,xmax,ymin,ymax,zmin,zmax;
extern int CurrentTimeStep,KeepScale;
static int E[50],code,NV,NE,NF,NT;

#define BUFFER_SIZE 8192

static int epReadFile( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
    int i,j,k,n,t = 0,total = 0,NamesGiven;

    FILE *fp;

    static vertex_t vertex;
    element_type_t *EL;

    static char *str,name[512],*ptr;

    double s,*NodeArray,*Velo,*Vabs,*Temp,*Pres,fdummy;

    double *Vector[1000], *Scalar[1000], *Times,*tvar;

    VARIABLE *Var1;

    group_t *grp;

    struct {
       int type;
       char name[256];
    } variable[64];

    int groupid,gid,StartTime=1,EndTime=1,IncTime=1,ToRead,last,gotit;
    static char groupname[128];

    group_t *group;
    extern Tcl_Interp *TCLInterp;

    if ( argc < 2 ) return TCL_ERROR;

    fp = fopen( argv[1], "r" ); 
    if ( !fp )
    {
       sprintf( interp->result, "ReadModel: can't open file [%s]\n",argv[1] );
       return TCL_ERROR;
    }

    if ( argc > 2 ) StartTime=atoi(argv[2]);
    if ( argc > 3 ) EndTime=atoi(argv[3]);
    if ( argc > 4 ) IncTime=atoi(argv[4]);

    str = (char *)malloc( BUFFER_SIZE*sizeof(char) );
    if ( !str ) {
      fprintf( stderr, "ERROR: ElmerPost: memory allocation error.\n" );
      exit(0);
    }

    NV = NE = NT = 0;
    fgets( str, BUFFER_SIZE-1, fp );
    sscanf( str, "%d %d %d %d", &NV,&NE,&NF,&NT );

    if ( NV <= 0 || NE <=0 )
    {
        Tcl_SetResult( interp, "Bad element model file.\n",TCL_STATIC );
        return TCL_ERROR;
    }

    ptr = str;
    for( i=0; i<4; i++ )
    {
        while( *ptr &&  isspace(*ptr) ) ptr++;
        while( *ptr && !isspace(*ptr) ) ptr++;
    }
    while( *ptr &&  isspace(*ptr) ) ptr++;

    i = 0;
    while ( *ptr )
    {
       if ( sscanf( ptr, "vector:%s", name ) == 1 )
       {
          variable[i].type = 1;
#if 0
          variable[i].name = (char *)malloc( strlen(name)+1 );
#endif
          strcpy( variable[i].name, name );

          i++;
       } else if ( sscanf( ptr, "scalar:%s", name ) == 1 )
       {
          variable[i].type = 2;
#if 0
          variable[i].name = (char *)malloc( strlen(name)+1 );
#endif
          strcpy( variable[i].name, name );

          i++;
       }

       while( *ptr && !isspace(*ptr) ) ptr++;
       while( *ptr &&  isspace(*ptr) ) ptr++;

       while( *ptr && !isspace(*ptr) ) ptr++;
       while( *ptr &&  isspace(*ptr) ) ptr++;
    }
    NamesGiven = i;

   if ( !CurrentObject->ElementModel )
     CurrentObject->ElementModel =  (element_model_t *)calloc( 1,sizeof(element_model_t) );

   if ( CurrentObject == &VisualObject )
     CurrentObject->ElementModel->NodeArray = NodeArray =
       MATR( var_new( "nodes", TYPE_DOUBLE, 3, NV ) );
   else
   {
      sprintf( str, "nodes_%s", CurrentObject->Name );

      CurrentObject->ElementModel->NodeArray = NodeArray =
       MATR( var_new( str, TYPE_DOUBLE, 3, NV ) );
   }

    Tcl_LinkVar( TCLInterp, "NumberOfTimesteps", (char *)&CurrentObject->ElementModel->NofTimesteps, TCL_LINK_INT );

    xmin = ymin = zmin =  DBL_MAX;
    xmax = ymax = zmax = -DBL_MAX;
    for( i=0; i<NV; i++ )
    {
         fgets(str,BUFFER_SIZE-1,fp);
		 if ( *str == '#' ) { i--; continue; }
		 
         sscanf( str, "%lf %lf %lf", &NodeArray[i],&NodeArray[NV+i],&NodeArray[2*NV+i] );

         xmin  = MIN( xmin, NodeArray[i] );
         ymin  = MIN( ymin, NodeArray[NV+i] );
         zmin  = MIN( zmin, NodeArray[2*NV+i] );

         xmax  = MAX( xmax, NodeArray[i] );
         ymax  = MAX( ymax, NodeArray[NV+i] );
         zmax  = MAX( zmax, NodeArray[2*NV+i] );
    }

    if (CurrentObject->ElementModel->Elements )
    {
        for( i=0; i<CurrentObject->ElementModel->NofElements; i++ )
            if (CurrentObject->ElementModel->Elements[i].Topology )
            {
                free( CurrentObject->ElementModel->Elements[i].Topology );
            }
        free(CurrentObject->ElementModel->Elements );
    }
   CurrentObject->ElementModel->Elements = Elements = (element_t *)calloc( NE,sizeof(element_t) );

    geo_free_groups(CurrentObject->ElementModel->Groups );
   CurrentObject->ElementModel->Groups = NULL;

    for( i=0; i<NE; i++ )
    {
        fgets(str,BUFFER_SIZE-1,fp);

        if ( *str == '#' )
        {
           if ( strncmp( str, "#group ", 7 ) == 0 ) 
           {
               sscanf( str, "#group %s", groupname );
               groupid = geo_group_id( CurrentObject->ElementModel,groupname,1 );
           } else if ( strncmp( str, "#endgroup ",10 ) == 0 )
           {
               sscanf( str, "#endgroup %s", groupname );
               groupid = geo_group_id( CurrentObject->ElementModel,groupname,0 );
           }
           i--;
           continue;
        }

        for( gid=0; gid<MAX_GROUP_IDS; gid++ )CurrentObject->ElementModel->Elements[i].GroupIds[gid] = -1;

        n = sscanf( str, "%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                groupname, &code, &E[0],&E[1],&E[2],&E[3],&E[4],&E[5],&E[6],&E[7],&E[8],&E[9],
                  &E[10],&E[11],&E[12],&E[13],&E[14],&E[15],&E[16],&E[17],&E[18],&E[19],
                       &E[20],&E[21],&E[22],&E[23],&E[24],&E[25],&E[26] );

        n = n - 2;

        while( n  < (code - 100 * (code / 100)) )
        {
          fgets( str,BUFFER_SIZE-1,fp );
          n = n + sscanf( str, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                            &E[n],&E[n+1],&E[n+2],&E[n+3],&E[n+4],&E[n+5],&E[n+6],&E[n+7],&E[n+8],&E[n+9],
                            &E[n+10],&E[n+11],&E[n+12],&E[n+13],&E[n+14],&E[n+15],&E[n+16],&E[n+17],&E[n+18],
                            &E[n+19],&E[n+20],&E[n+21],&E[n+22],&E[n+23],&E[n+24],&E[n+25],&E[n+26] );
        }


        groupid = geo_group_id( CurrentObject->ElementModel,groupname,1 );

        for( EL=ElementDefs.ElementTypes; EL!=NULL; EL=EL->Next )
            if ( code == EL->ElementCode )
            {
                Elements[total].ElementType = EL;
                Elements[total].DisplayFlag = TRUE;

                for( gid=0; gid<MAX_GROUP_IDS; gid++ ) if ( Elements[total].GroupIds[gid] < 0 ) break;

                groupid = 0;
                for( grp=CurrentObject->ElementModel->Groups; grp != NULL; grp=grp->Next, groupid++ )
                   if ( grp->Open ) if ( gid < MAX_GROUP_IDS ) Elements[total].GroupIds[gid++]  = groupid;

                Elements[total].Topology = (int *)malloc( EL->NumberOfNodes*sizeof(int) );

                if ( !Elements[total].Topology )
                {
                    Tcl_SetResult( interp,"FATAL: can't alloc element connection tables.\n",TCL_STATIC );
                    return TCL_ERROR;
                }

                for( j=0; j<EL->NumberOfNodes; j++ ) Elements[total].Topology[j] = E[j];

                total++;

                break;
            }

        groupid = geo_group_id( CurrentObject->ElementModel,groupname,0 );

        if ( EL == NULL )
        {
          if ( code != 101 )
            fprintf( stderr,"Unknown element type: [%d]. Skipping Element.\n", code );
        }
    }

   CurrentObject->ElementModel->NofElements  = NE = total;
   CurrentObject->ElementModel->NofNodes     = NV;
   CurrentObject->ElementModel->NofTimesteps = (EndTime-StartTime+IncTime) / IncTime;
   ToRead = CurrentObject->ElementModel->NofTimesteps;

    last = 0;
    t = 0;
    if ( NF>0 )
    {
       Times = malloc( 3*ToRead*sizeof(double) );
       for( i=0; i<ToRead; i++ )
       {
           Times[0*ToRead+i] = i+1;
           Times[1*ToRead+i] = i+1;
           Times[2*ToRead+i] = i+1;
       }

       if ( !NamesGiven ) {
         if ( NF>=3 )
         {
            Vabs = MATR( var_new( "vabs", TYPE_DOUBLE, 1, ToRead*NV ) );
            Velo = MATR( var_new( "velo", TYPE_DOUBLE, 3, ToRead*NV ) );
         }
         Pres = MATR( var_new( "pres", TYPE_DOUBLE, 1, ToRead*NV ) );
         Temp = MATR( var_new( "temp", TYPE_DOUBLE, 1, ToRead*NV ) );

         if ( NF>=3 )
         {
   	    t = 0;
            for( i=0; i<EndTime; i++ )
            {
               if ( feof(fp) ) goto exit_loop;

               if ( i<StartTime-1 || ((i-StartTime+1)%IncTime) ) {
		   for( k=0; k<NV; k++ ) {
 		      fgets( str,BUFFER_SIZE-1,fp );
                      if ( feof(fp) ) goto exit_loop;
		      if ( *str == '#' ) k--;
		   }
  	       } else {
		   for( k=0; k<NV; k++ )
		   {
	 	      fgets( str,BUFFER_SIZE-1,fp );
                      if ( feof( fp ) ) goto exit_loop;
		      if ( *str == '#' ) { k--; continue; }
		          sscanf( str, "%lf %lf %lf %lf %lf",&Velo[t*NV+k], &Velo[NV*(t+ToRead)+k],
                               &Velo[NV*(t+2*ToRead)+k], &Pres[t*NV+k],&Temp[t*NV+k] );
		   }
		   t++;
 	       }
            }

            for( i=0; i<NV*ToRead; i++ )
            {
               Vabs[i] = sqrt( Velo[i]*Velo[i]+Velo[NV*ToRead+i]*Velo[NV*ToRead+i]+
				   Velo[2*NV*ToRead+i]*Velo[2*NV*ToRead+i] );
            }
         } else {
	     t = 0; last = 0;
             for( i=0; i<EndTime; i++ )
             {
                if ( feof(fp) ) goto exit_loop;
                 if ( i < StartTime-1 || ((i-StartTime+1)%IncTime) )
		 {
		    for( k=0; k<NV; k++ ) {
	  	       fgets( str,BUFFER_SIZE-1,fp );
                       if ( feof( fp ) ) goto exit_loop;
		       if ( *str == '#' ) k--;
		    }
		 } else {
		    for( k=0; k<NV; k++ )
		    {
		       fgets( str,BUFFER_SIZE-1,fp );
                       if ( feof( fp ) ) goto exit_loop;
		       if ( *str == '#' ) { k--; continue; }
		       sscanf( str, "%lf %lf", &Pres[t*NV+k],&Temp[t*NV+k]  );
		    }
	            t++;
		 }
             }
	 }
       } else {
          for( i=0; i<NamesGiven; i++ ) {
             if ( variable[i].type == 1 ) {
                Vector[i] = MATR( var_new( variable[i].name, TYPE_DOUBLE, 3, ToRead*NV ) );
                strcpy( str, variable[i].name ); strcat( str, "_abs" );
                Scalar[i] = MATR( var_new( str, TYPE_DOUBLE, 1, ToRead*NV ) );
             } else {
                Scalar[i] = MATR( var_new( variable[i].name, TYPE_DOUBLE, 1, ToRead*NV ) );
             }
          }

	  t = 0;
	  for( i=0; i<EndTime; i++ )
          {
#if 0
             if ( feof(fp) ) goto exit_loop;
             if ( i<StartTime-1 || ((i-StartTime+1)%IncTime) )
	     {
                 for( j=0; j<NV*NF; j++ )
		 {
                   if ( feof(fp) ) goto exit_loop;
		    while ( 1 != fscanf( fp, "%lf", &fdummy ) ) 
		    {
			 fgets( str, BUFFER_SIZE-1, fp );
                          if ( feof( fp ) ) goto exit_loop;
                         if ( strncmp( str, "#time ", 6 ) == 0 )
                         {
                             double t1,t2,t3;

                             sscanf( str, "#time %lf %lf %lf", &t1,&t2,&t3 );
                             Times[0*ToRead+t] = t1;
                             Times[1*ToRead+t] = t2;
                             Times[2*ToRead+t] = t3;
                         }
		    }
		 }
             } else
#endif
             {
		for( k=0; k<NV; k++ )
		{
                   if ( feof( fp ) ) goto exit_loop;
		   for( j=0; j<NamesGiven; j++ ) 
		   {
                     if ( feof(fp) ) goto exit_loop;

		      if ( variable[j].type == 1 ) 
                      {
                          if ( feof(fp) ) goto exit_loop;
                          while ( 1 != fscanf( fp, "%lf",&Vector[j][t*NV+k] ) ) {
			      fgets( str,BUFFER_SIZE-1,fp );
                              if ( feof( fp ) ) goto exit_loop;
                              if ( strncmp( str, "#time ", 6 ) == 0 )
                              {
                                  double t1,t2,t3;

                                  sscanf( str, "#time %lf %lf %lf", &t1,&t2,&t3 );
                                  Times[0*ToRead+t] = t1;
                                  Times[1*ToRead+t] = t2;
                                  Times[2*ToRead+t] = t3;
                              }
			  }

                          if ( feof(fp) ) goto exit_loop;
                          while ( 1 != fscanf( fp, "%lf",&Vector[j][(t+ToRead)*NV+k] ) ) {
		   	      fgets( str,BUFFER_SIZE-1,fp );
                              if ( feof(fp) ) goto exit_loop;
                              if ( strncmp( str, "#time ", 6 ) == 0 )
                              {
                                  double t1,t2,t3;

                                  sscanf( str, "#time %lf %lf %lf", &t1,&t2,&t3 );
                                  Times[0*ToRead+t] = t1;
                                  Times[1*ToRead+t] = t2;
                                  Times[2*ToRead+t] = t3;
                              }
			  }

                          if ( feof(fp) ) goto exit_loop;
                          while ( 1 != fscanf( fp, "%lf",&Vector[j][(t+2*ToRead)*NV+k] ) ) {
			      fgets( str,BUFFER_SIZE-1,fp );
                              if ( feof(fp) ) goto exit_loop;
                              if ( strncmp( str, "#time ", 6 ) == 0 )
                              {
                                  double t1,t2,t3;

                                  sscanf( str, "#time %lf %lf %lf", &t1,&t2,&t3 );
                                  Times[0*ToRead+t] = t1;
                                  Times[1*ToRead+t] = t2;
                                  Times[2*ToRead+t] = t3;
                              }
			  }
   
                          Scalar[j][t*NV+k] = sqrt( 
                               Vector[j][t*NV+k]*Vector[j][t*NV+k]+
                               Vector[j][NV*(t+ToRead)+k]*Vector[j][NV*(t+ToRead)+k] +
                               Vector[j][NV*(t+2*ToRead)+k]*Vector[j][NV*(t+2*ToRead)+k] );

                      } else {
                          if ( feof(fp) ) goto exit_loop;
			  while ( 1 != fscanf( fp, "%lf",&Scalar[j][t*NV+k] ) )
                          {
			      fgets( str,BUFFER_SIZE-1,fp );
                              if ( feof(fp) ) goto exit_loop;
                              if ( strncmp( str, "#time ", 6 ) == 0 )
                              {
                                  double t1,t2,t3;

                                  sscanf( str, "#time %lf %lf %lf", &t1,&t2,&t3 );
                                  Times[0*ToRead+t] = t1;
                                  Times[1*ToRead+t] = t2;
                                  Times[2*ToRead+t] = t3;
                              }
			  }
                      }
		   }
		}
                last  = 0;
                if ( i>=StartTime-1 && !((i-StartTime+1)%IncTime) ) { last=1; t++; }
                if ( t >= ToRead ) goto exit_loop;
             }
          }
       }
    }
exit_loop:

    if ( t == 0 && ToRead > 0 && i > 0 ) t++;

    fclose( fp );

    if ( t > 0 ) {
      CurrentObject->ElementModel->NofTimesteps = t;
      if ( CurrentObject == &VisualObject )
        tvar = MATR( var_new( "times", TYPE_DOUBLE, 3, t ) );
      else
      {
        sprintf( str, "times_%s", CurrentObject->Name );
        tvar = MATR( var_new( str, TYPE_DOUBLE, 3, t ) );
      }

      for( i=0; i<t; i++ ) {
        tvar[0*t+i] = Times[0*ToRead+i];
        tvar[1*t+i] = Times[1*ToRead+i];
        tvar[2*t+i] = Times[2*ToRead+i];
      }

      free( Times );
    }

   if ( t != ToRead ) {
     fprintf( stderr,"WARNING: ElmerPost: Not enough data for all timesteps"
              " requested. REQUEST: %d, GOT: %d\n", ToRead, t );

      for( i=0; i<NamesGiven; i++ )
      {
         if ( t > 0 ) {
           if ( variable[i].type == 1 )
           {
             Var1 = lst_find( VARIABLES, variable[i].name );
#if 0
             NCOL(Var1) = t*NV;
#endif
 
             strcpy( str, variable[i].name ); strcat( str, "_abs" );
             Var1 = lst_find( VARIABLES, str );
#if 0
             NCOL(Var1) = t*NV;
#endif
           } else {
             Var1 = lst_find( VARIABLES, variable[i].name );
#if 0
             NCOL(Var1) = t*NV;
#endif
           }
         } else {
           if ( variable[i].type == 1 )
           {
             var_delete( variable[i].name );
            strcpy( str, variable[i].name ); strcat( str, "_abs" );
             var_delete( str  );
           } else {
             var_delete( variable[i].name );
           }
         }
      }
   }

    groupid = 0;
    for( group =CurrentObject->ElementModel->Groups; group!=NULL; group=group->Next,groupid++ )
    {
       sprintf( name,  "Groups(%d)", groupid );
       Tcl_SetVar( TCLInterp, name, group->Name, TCL_GLOBAL_ONLY );
    }

    sprintf( name, "NumberOfGroups" );
    sprintf( str,  "%d", groupid );
    Tcl_SetVar( TCLInterp, name, str, TCL_GLOBAL_ONLY );

    if ( CurrentObject->Geometry )
    {
       geo_free_edge_tables( CurrentObject->Geometry );
       geo_free_vertex_face_tables( CurrentObject->Geometry );
    } else {
       CurrentObject->Geometry = (geometry_t *)calloc( 1,sizeof(geometry_t) );
    }

    CurrentObject->Geometry->VertexCount = 0;

    if ( !KeepScale ) {
       s = MAX( MAX( xmax-xmin, ymax-ymin ), zmax-zmin );
    } else {
       s = CurrentObject->Geometry->Scale;
       xmin = CurrentObject->Geometry->MinMax[0].x[0];
       ymin = CurrentObject->Geometry->MinMax[0].x[1];
       zmin = CurrentObject->Geometry->MinMax[0].x[2];

       xmax = CurrentObject->Geometry->MinMax[1].x[0];
       ymax = CurrentObject->Geometry->MinMax[1].x[1];
       zmax = CurrentObject->Geometry->MinMax[1].x[2];
    }
    XMin = YMin = ZMin =  DBL_MAX;
    XMax = YMax = ZMax = -DBL_MAX;

    for( i=0; i<NV; i++ )
    {
        vertex.x[0] = ( 2.0 * ( NodeArray[i] - xmin) - (xmax - xmin)) / s;
        vertex.x[1] = ( 2.0 * ( NodeArray[NV+i] - ymin) - (ymax - ymin)) / s;
        vertex.x[2] = ( 2.0 * ( NodeArray[2*NV+i] - zmin) - (zmax - zmin)) / s;

        XMin = MIN( XMin,vertex.x[0] );
        YMin = MIN( YMin,vertex.x[1] );
        ZMin = MIN( ZMin,vertex.x[2] );

        XMax = MAX( XMax,vertex.x[0] );
        YMax = MAX( YMax,vertex.x[1] );
        ZMax = MAX( ZMax,vertex.x[2] );

        vertex.ElementModelNode = TRUE;
        geo_add_vertex( CurrentObject->Geometry, &vertex );
    }

    CurrentObject->Geometry->Scale = s;
    CurrentObject->Geometry->MinMax[0].x[0] = xmin;
    CurrentObject->Geometry->MinMax[0].x[1] = ymin;
    CurrentObject->Geometry->MinMax[0].x[2] = zmin;

    CurrentObject->Geometry->MinMax[1].x[0] = xmax;
    CurrentObject->Geometry->MinMax[1].x[1] = ymax;
    CurrentObject->Geometry->MinMax[1].x[2] = zmax;

    CurrentObject->Geometry->Edges = (edge_t *)calloc( CurrentObject->Geometry->VertexCount, sizeof(edge_t) );

    if ( !CurrentObject->Geometry->Edges )
    {
        fprintf( stderr, "Can't alloc edge array.\n" );
        exit( 0 );
    }

    CurrentObject->Geometry->TriangleCount = 0;

    for( i=0; i<NE; i++ ) 
    {
        EL = CurrentObject->ElementModel->Elements[i].ElementType;
        (*EL->Triangulate)( CurrentObject->Geometry,(element_t *)&Elements[i],(element_t *)&Elements[i] );
    }

    geo_vertex_normals( CurrentObject->Geometry,50.0 );

    CurrentTimeStep = 0;

    UpdateObject( 0,NULL,0,NULL );
    DrawItSomeTimeWhenIdle();

    free( str );

    return TCL_OK;
}

#if 0
static VARIABLE *SetModel( VARIABLE *ptr )
{
    static vertex_t vertex;
    element_type_t *EL;

    double s,*NodeArray,*Velo,*Vabs,*Temp,*Pres;

    int i,j,k,total;

    double *Topology  = MATR(NEXT(ptr));
    double *Type      = MATR(NEXT(NEXT(ptr)));

    NV = NCOL(ptr);

    NE = NROW(NEXT(ptr));

    if ( NEXT(NEXT(NEXT(ptr))) )
        NT = M(NEXT(NEXT(NEXT(ptr))),0,0);
    else
       NT = 1;

   CurrentObject->ElementModel->NodeArray = NodeArray = MATR( ptr );

    xmin = ymin = zmin =  DBL_MAX;
    xmax = ymax = zmax = -DBL_MAX;
    for( i=0; i<NV; i++ )
    {
         xmin  = MIN( xmin, NodeArray[i] );
         ymin  = MIN( ymin, NodeArray[NV+i] );
         zmin  = MIN( zmin, NodeArray[2*NV+i] );

         xmax  = MAX( xmax, NodeArray[i] );
         ymax  = MAX( ymax, NodeArray[NV+i] );
         zmax  = MAX( zmax, NodeArray[2*NV+i] );
    }

    if (CurrentObject->ElementModel.Elements )
    {
        for( i=0; i<ElementModel.NofElements; i++ )
            if (CurrentObject->ElementModel.Elements[i].Topology )
            {
                free( CurrentObject->ElementModel.Elements[i].Topology );
            }
        free(CurrentObject->ElementModel.Elements );
    }
   CurrentObject->ElementModel.Elements = Elements = (element_t *)calloc( NE,sizeof(element_t) );

    total = 0;
    for( i=0; i<NE; i++ )
    {
        for( EL=ElementDefs.ElementTypes; EL != NULL; EL=EL->Next )
            if ( Type[i] == EL->ElementCode )
            {
                Elements[total].ElementType = EL;

                Elements[total].Topology = (int *)malloc( EL->NumberOfNodes*sizeof(int) );

                if ( !Elements[total].Topology )
                {
                    error( "FATAL: can't alloc element connection tables.\n" );
                }

                for( j=0; j<EL->NumberOfNodes; j++ ) Elements[total].Topology[j] = Topology[i*NCOL(NEXT(ptr))+j];
                total++;

                break;
            }

        if ( EL == NULL )
        {
            fprintf( stderr, "Unknown element type: [%d]. Skipping Element.\n", Type[i] );
        }
    }

   CurrentObject->ElementModel.NofElements   = NE = total;
   CurrentObject->ElementModel.NofNodes      = NV;
   CurrentObject->ElementModel.NofTimesteps  = NT;

    geo_free_edge_tables( &Geometry );
    geo_free_vertex_face_tables( &Geometry );

    Geometry.VertexCount = 0;

    s = MAX( MAX( xmax-xmin, ymax-ymin ), zmax-zmin );
    XMin = YMin = ZMin =  DBL_MAX;
    XMax = YMax = ZMax = -DBL_MAX;

    for( i=0; i<NV; i++ )
    {
        vertex.x[0] = ( 2.0 * ( NodeArray[i] - xmin) - (xmax - xmin)) / s;
        vertex.x[1] = ( 2.0 * ( NodeArray[NV+i] - ymin) - (ymax - ymin)) / s;
        vertex.x[2] = ( 2.0 * ( NodeArray[2*NV+i] - zmin) - (zmax - zmin)) / s;

        XMin = MIN( XMin,vertex.x[0] );
        YMin = MIN( YMin,vertex.x[1] );
        ZMin = MIN( ZMin,vertex.x[2] );

        XMax = MAX( XMax,vertex.x[0] );
        YMax = MAX( YMax,vertex.x[1] );
        ZMax = MAX( ZMax,vertex.x[2] );

        vertex.ElementModelNode = TRUE;
        geo_add_vertex( &Geometry, &vertex );
    }

    Geometry.MinMax[0].x[0] = xmin;
    Geometry.MinMax[0].x[1] = ymin;
    Geometry.MinMax[0].x[2] = zmin;

    Geometry.MinMax[1].x[0] = xmax;
    Geometry.MinMax[1].x[1] = ymax;
    Geometry.MinMax[1].x[2] = zmax;

    Geometry.Edges = (edge_t *)calloc( Geometry.VertexCount, sizeof(edge_t) );

    if ( !Geometry.Edges )
    {
        error( "SetModel: FATAL: Can't alloc edge array.\n" );
    }

    Geometry.TriangleCount = 0;

    for( i=0; i<NE; i++ ) 
    {
        EL = Elements[i].ElementType;
        (*EL->Triangulate)( &Geometry,(element_t *)&Elements[i],(element_t *)&Elements[i] );
    }

    geo_vertex_normals( &Geometry,50.0 );

    CurrentTimeStep = 0;

    fprintf( stderr, "TRIANGLES: %d %d, TS: %d\n", Geometry.TriangleCount,Geometry.MaxTriangleCount,NT );

    UpdateObject( 0,NULL,0,NULL );
    DrawItSomeTimeWhenIdle();

    return (VARIABLE *)NULL;
}
#endif

int Readfile_Init( Tcl_Interp *interp )
{
    Tcl_CreateCommand( interp,"cReadFile",epReadFile,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
#if 0
    com_init
        (
             "model", FALSE, FALSE,SetModel,3,4,
         "Usage: model(node-array,topology,element-type,[timesteps])\nSet current element model."
        );
#endif

    return TCL_OK;
}
