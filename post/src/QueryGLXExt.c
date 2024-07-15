#include <stdio.h>
#include <stdlib.h>
#include <GL/glx.h>

int main( int argc, char **argv )
{
   Display *display;
   int err, event, server_major, server_minor, client_major, client_minor;

   if ( argc > 1 )
     display = XOpenDisplay( argv[1] );
   else
     display = XOpenDisplay( "" );

   if ( !display ) 
   {
     fprintf( stderr, "FATAL: Can't connect to X Server\n" );
     fprintf( stdout, "fatal\n" );
     exit(0);
   }

   if ( glXQueryExtension( display, &err, &event ) )
   {
     sscanf( glXQueryServerString( display, 0, GLX_VERSION ),
          "%d.%d", &server_major, &server_minor );

     sscanf( glXGetClientString( display, GLX_VERSION ), 
          "%d.%d", &client_major, &client_minor );

     if ( server_major != client_major || server_minor != client_minor )
     {
        fprintf( stdout, "glx_mismatch\n" );
     } else
        fprintf( stdout, "success\n" );
   } else
     fprintf( stdout, "no_glx_found\n" );

   XCloseDisplay( display );
   return 1;
}
