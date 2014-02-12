#!/bin/bash
#
# Juha Vierinen 2005
#

#
# Modules to fetch and build (in this order)
#
if test "$ELMER_MODULES" = ""; then
    modules="umfpack matc meshgen2d elmergrid mathlibs eio hutiter fem"
else
    modules=$ELMER_MODULES
fi

#
# if anything resembling a help flag is given, we print help 
#
if test "`echo $* |grep "\-h"`" != ""; then 
    printf "\n"
    printf " -----------------------\n"
    printf " Elmer build script v0.1\n"
    printf " -----------------------\n"
    printf "\n"
    printf "Fetches elmer from cvs into: ./build-HOST-DATE\n"
    printf "and compiles it into directory under: /tmp/username/build-HOST-DATE\n"
    printf "\n"
    printf "DATE is replaced with the current date and HOST with hostname.\n"
    printf "It is possible to do several builds at the same time,\n"
    printf "\n"
    printf "If ELMER_USER is unset, whoami (in the future anonymous) will be used to\n"
    printf "fetch the source.\n"
    printf "\n"
    exit
fi

#
# Try to use ELMER_USER as cvs user, otherwise revert to whoami
#
if test "$ELMER_USER" != ""; then
    export CVSROOT=$ELMER_USER@ibmsc.csc.fi:/fs/proj1/elmer/elmer-opensource    
else
if test "$ELMER_CVSROOT" != ""; then
    export CVSROOT=$ELMER_CVSROOT
else
    export CVSROOT=$USER@ibmsc.csc.fi:/fs/proj1/elmer/elmer-opensource
fi
fi

printf "Fetching source from CVS using CVSROOT=%s\n" $CVSROOT
export CVS_RSH="ssh"

#
# uset ELMER_HOME and remove it from ld_library_path to avoid confusion
#
if test "$ELMER_HOME" != ""; then
    printf "Variable ELMER_HOME, will be cleared for the duration of the build.\n"
    printf "Also, \$ELMER_HOME/lib will be removed from LD_LIBRARY_PATH.\n"
    OLD_ELMER_HOME=$ELMER_HOME
    unset ELMER_HOME
    OLD_LD_LIBRARY_PATH=$LD_LIBRRARY_PATH
    LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | sed -e "s,$ELMER_HOME/lib[/]*,,g"`
    export LD_LIBRARY_PATH

    printf "LD_LIBRARY_PATH is now: %s\n" $LD_LIBRARY_PATH
    printf "Original values will be restored after the compilation.\n"
fi

#
# Unset ELMER_LIB to avoid confusion
#
if test "$ELMER_LIB" != ""; then
    OLD_ELMER_LIB=$ELMER_LIB
    unset ELMER_LIB
fi

#
# If no processor number set, use only 1 (useful on clusters)
#
if test "$NPROCS" = ""; then
    NPROCS=1
fi
printf "Using %s processors for compiling\n" $NPROCS

#
# Create temp names
#
datestr=`date '+%Y-%m-%d-%H-%M-%S'`
tmpdir=build.`hostname`.$datestr
mkdir -p $tmpdir

topdir=`pwd`

TESTPREFIX="${topdir}/${tmpdir}/dist"

printf "Using %s as build name\n" $TESTPREFIX

cd $tmpdir

#
# Get modules from cvs, configure ; make ; make install into temp directory
#
cvs co $modules
for m in $modules; do
    cd $m
    if test "$USE_OWN_MATHLIBS" = yes; then
	./configure --prefix=$TESTPREFIX $CONFFLAGS --with-blas=$TESTPREFIX/lib/libblas.a --with-lapack=$TESTPREFIX/lib/liblapack.a
    else
	./configure --prefix=$TESTPREFIX $CONFFLAGS 
    fi
    if test "`which gmake`" != ""; then
	gmake -j$NPROCS
    else
	make -j$NPROCS
    fi
    make install
    make dist
    if test "$m" = fem; then
    if test "$ELMER_CHECK" != "no"; then
	make check
    fi
    fi
    cd ..
done

#
# Copy source distribution files to dist dir
#
cp `find . -name *.tar.gz` dist

#
# Create binary distribution
#
distArch=`uname -m`
cd dist; tar cvf - bin include lib share |gzip > "elmer-bin-${distArch}.tar.gz"
cd ..

cd $topdir
