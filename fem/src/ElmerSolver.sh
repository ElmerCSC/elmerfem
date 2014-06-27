#!/bin/sh

if [ -z ${ELMER_HOME} ]; then
  PREFIX=
else
  PREFIX="${ELMER_HOME}/bin/"
fi

exec ${PREFIX}ElmerSolver_mpi $*
