#!/bin/sh -f

if test "$ELMER_LIB" = ""; then
  if  test "$ELMER_HOME" = ""; then
     LIBDIR=@ELMER_SOLVER_HOME@/../../@ELMER_INSTALL_LIB_DIR@
     INCLUDE=@ELMER_SOLVER_HOME@/include
  else
     LIBDIR=$ELMER_HOME/@ELMER_INSTALL_LIB_DIR@
     INCLUDE=$ELMER_HOME/share/elmersolver/include
  fi
else
  LIBDIR=$ELMER_LIB
  INCLUDE=$ELMER_LIB/../include
fi

if test "$ELMER_Fortran_COMPILER" = ""; then
  FC=@CMAKE_Fortran_COMPILER@
else
  FC=$ELMER_Fortran_COMPILER
fi


if test @HAVE_ELMERICE@ = "TRUE"; then
    ELMERICE_LIB=$LIBDIR/../../share/elmersolver/lib
    LIBELMERICE="-Xlinker -rpath -Xlinker $ELMERICE_LIB $ELMERICE_LIB/ElmerIceSolvers@SHL_EXTENSION@ $ELMERICE_LIB/ElmerIceUSF@SHL_EXTENSION@"
    printf "with elmerice\n"
else
    LIBELMERICE=""
    printf "no elmerice\n"
fi

if test @HAVE_MMG@ = "TRUE"; then
    MMGLIBDIR="-L@MMG_LIBDIR@"
    MMGINCLUDE="-I@MMG_INCLUDE_DIR@"
    printf "with MMG\n"
else
    MMGLIBDIR=""
    MMGINCLUDE=""
fi

if test @HAVE_PARMMG@ = "TRUE"; then
    PARMMGLIBDIR="-L@PARMMG_LIBDIR@"
    PARMMGINCLUDE="-I@PARMMG_INCLUDE_DIR@"
    printf "with ParMMG\n"
    if test "$MMGLIBDIR" = "$PARMMGLIBDIR"; then
	PARMMGLIBDIR=""
	printf "MMG and ParMMG share the same lib dir\n"
    fi    
    if test "$MMGINCLUDE" = "$PARMMGINCLUDE"; then
	PARMMGINCLUDE=""
	printf "MMG and ParMMG share the same include dir\n"
    fi    
else
    PARMMGLIBDIR=""
    PARMMGINCLUDE=""
fi

cmd="$FC $* @CMAKE_Fortran_FLAGS@ @ELMER_F90FLAGS@ @CMAKE_SHARED_LIBRARY_Fortran_FLAGS@ @CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS@ -I$INCLUDE -L$LIBDIR $LIBELMERICE $MMGINCLUDE $MMGLIBDIR $PARMMGINCLUDE $PARMMGLIBDIR -shared -lelmersolver "
printf "%s " $cmd
printf "\n"
$FC $* @CMAKE_Fortran_FLAGS@ @ELMER_F90FLAGS@ @CMAKE_SHARED_LIBRARY_Fortran_FLAGS@ @CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS@ -I$INCLUDE -L$LIBDIR $LIBELMERICE $MMGINCLUDE $MMGLIBDIR $PARMMGINCLUDE $PARMMGLIBDIR -shared -lelmersolver 
    
# $FC @CMAKE_Fortran_FLAGS@ @ELMERF90_INCLUDE_DIRS@ $*
#@USE_WINDOWS_COMPILER_FALSE@$FC @FCFLAGS@ @SH_LINKING_TO_FLAGS@ @INCLUDE_MODULE_FLAG@$INCLUDE @B64FLAGS@ @SH_LDFLAGS@ @EXTRA_LIBS@ $*
#@USE_WINDOWS_COMPILER_TRUE@$FC @FCFLAGS@ @INCLUDE_MODULE_FLAG@$INCLUDE @B64FLAGS@ @SH_LDFLAGS@ @EXTRA_LIBS@ $* -L$LIBDIR -lelmersolver
