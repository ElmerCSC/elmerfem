#!/bin/csh -f
# set basename = %1   
# set directory = %2  
# set ProblemDirectory = %3

cd $2
setenv LD_LIBRARY_PATH $3
ln -s $3/libcxa.so.3 .
$3/gid2elmer < $1.dat
$3/findparents
/bin/mv mesh.boundary.corrected mesh.boundary
/bin/rm libcxa.so.3
