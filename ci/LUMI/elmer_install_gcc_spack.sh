#!/bin/bash -l
umask 0002
echo "########################################################"
echo "# This is the automatic build of Elmer on LUMI "
echo "# "
echo "# process launched:" $(date)
echo "#######################################################"
# set compiler suite
export COMPILER="gcc/11.2.0"
export MPI="cray-mpich/8.1.23"
ml purge
module load PrgEnv-gnu/8.3.3
module load $COMPILER $MPI cray-libsci/22.12.1.1
ml use /appl/local/csc/soft/eng/elmer/spack/23.03/0.19.2/modules/tcl/linux-sles15-zen2
ml use /appl/lumi/spack/23.03/0.19.2/share/spack/modules/linux-sles15-zen2
ml load netcdf-c/4.9.0-gcc-znb netcdf-fortran/4.6.0-gcc-waz42 hdf5/1.12.2-gcc-5ut mumps/5.5.1-gcc-ypwio scotch/7.0.1-gcc-pk2f5 hypre/2.26.0-gcc-kmjvj boost/1.80.0-gcc-423 mmg/5.6.0-gcc-z2ip3
BRANCH="devel"

# set sources, script-dir and inquire build and installation directory (set these!)
BASEDIR="$PWD"
ELMERSRC="${BASEDIR}/elmerfem"
SCRIPTDIR="${BASEDIR}"
cd ${ELMERSRC}
# updating version 
echo "checking for updates"
git checkout $BRANCH
#git fetch 
#git pull
git status -uno


ELMERDEP="/appl/local/csc/soft/eng/elmer/spack/elmerdependencies/"
XIOS_DIR=$ELMERDEP

VERSION=$(git log -1 --pretty=format:%h)
TIMESTAMP=$(date +"%m-%d-%y")
ELMER_REV="Elmer_${BRANCH}_${VERSION}_${TIMESTAMP}"
IDIR=${ELMER_REV}
IPATH="/appl/local/csc/soft/eng/elmer/spack/elmer/${COMPILER/\//-}/${MPI/\//-}/${IDIR}"
export ELMER_HOME=${IPATH}
BUILDDIR="${BASEDIR}/spack_${ELMER_REV}"

echo "#######################################################"
echo "# installation into " ${IPATH}
echo "#######################################################"
# set up modules
echo "building with following modules:"




echo "-------------------------------------"
# create new build-dir
echo "creating"  $BUILDDIR ":"
echo "-------------------------------------"
if [[ ! -e $BUILDDIR ]]; then
    mkdir $BUILDDIR

  # configure
  cd ${BUILDDIR}
  pwd

  echo "configuring:"
  CMAKE=cmake

  TOOLCHAIN=${SCRIPTDIR}/Elmer-linux-gcc11.cmake
  # all main definitions are in this file, so no Precache needed
  # PRECACHE=${SCRIPTDIR}/Elmer-linux-precache.cmake

  echo "-------------------------------------"
  echo "Building Elmer from source " ${ELMERSRC}
  echo "within build directory " ${BUILDDIR}
  #echo "using following toolchain file " ${TOOLCHAIN}
  echo "installation into " ${IDIR}
  echo "Elmer version is:" ${VERSION}
  echo "using pre-cache file"  ${PRECACHE}
  echo "-------------------------------------"
  LOGFILE="${BASEDIR}/Logs/installation_${BRANCH}_${TIMESTAMP}.log"
  echo "-------------------------------------" > ${LOGFILE}
  echo "Building Elmer from source " ${ELMERSRC} >> ${LOGFILE}
  echo "within build directory " ${BUILDDIR} >> ${LOGFILE}
  echo "using following toolchain file " ${TOOLCHAIN} >> ${LOGFILE}
  echo "installation into " ${IDIR} >> ${LOGFILE}
  echo "Elmer version is:" ${VERSION} >> ${LOGFILE}
  echo "using pre-cache file"  ${PRECACHE} >> ${LOGFILE}
  echo "-------------------------------------" >> ${LOGFILE}

  if ! $CMAKE $ELMERSRC -DCMAKE_TOOLCHAIN_FILE=$TOOLCHAIN -DCMAKE_INSTALL_PREFIX=$IPATH -Wno-dev \
       -DELMER_SOLVER_HOME=${IPATH} \
       -DWITH_MPI:BOOL=TRUE \
       -DWITH_LUA:BOOL=TRUE \
       -DWITH_OpenMP:BOOL=TRUE \
       -DWITH_Zoltan:BOOL=TRUE \
       -DBLAS_LIBRARIES="-L ${LIBSCI_BASE_DIR}/gnu/9.1/x86_64/lib -lsci_gnu" \
       -DLAPACK_LIBRARIES="-L ${LIBSCI_BASE_DIR}/gnu/9.1/x86_64/lib -lsci_gnu" \
       -DHDF5_INCLUDE_DIR="${HDF5_INSTALL_ROOT}/include" \
       -DHDF5_LIBRARY="${HDF5_INSTALL_ROOT}/lib/libhdf.so" \
       -DWITH_Mumps:BOOL=TRUE \
       -DMUMPS_ROOT="$MUMPS_ROOT" \
       -DSCALAPACK_LIBRARIES="-L${LIBSCI_BASE_DIR}/gnu/9.1/x86_64/lib -lsci_gnu_mpi_mp -lsci_gnu" \
       -DWITH_Hypre:BOOL=TRUE \
       -DHYPREROOT="$HYPRE_INSTALL_ROOT" \
       -DWITH_Trilinos:BOOL=FALSE \
       -DWITH_NETCDF:BOOL=TRUE \
       -DNETCDF_INCLUDE_DIR="${NETCDF_C_INSTALL_ROOT}/include;${NETCDF_FORTRAN_INSTALL_ROOT}/include" \
       -DNETCDF_LIBRARY="${NETCDF_C_INSTALL_ROOT}/lib/libnetcdf.so" \
       -DNETCDFF_LIBRARY="${NETCDF_FORTRAN_INSTALL_ROOT}/lib/libnetcdff.so" \
       -DWITH_MMG:BOOL=TRUE \
       -DMMG_root=="$MMG_ROOT" \
       -DWITH_XIOS:BOOL=TRUE \
       -DXIOS_INCLUDE_DIR=${XIOS_DIR}/inc \
       -DXIOS_ROOT=${XIOS_DIR}  \
       -DXIOS_LIBRARY="${XIOS_DIR}/lib/libxios.so" \
       -DWITH_ELMERICE:BOOL=TRUE |& tee -a ${LOGFILE};then
      echo "Configuration unsuccessful, exiting..."
      exit 1
  fi 
  echo "done configuration"
  make -j 4 && make install
elif [[ ! -d $dir ]]; then
    echo $BUILDDIR "exists - updating" 1>&2
    cd ${BUILDDIR}
    pwd
fi

