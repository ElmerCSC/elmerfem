#!/bin/sh -f 
#
# List of files not to be included in open source release.
#
# This will be tricky, because they also have to be removed from the Makefile 
# 
SOLVERS="PhaseChange StatMagSolve Acoustics"
DIRS="tests/StaticMag tests/phase tests/acoustics"

echo ${distdir}
cd ${distdir}

echo `pwd`
for s in $SOLVERS; do
        cd src 
#	echo "grep -v $s Makefile.in > Makefile.tmp ; mv Makefile.tmp Makefile.in"
	grep -v $s Makefile.in > Makefile.tmp ; mv Makefile.tmp Makefile.in
#	echo "grep -v $s Makefile.am > Makefile.tmp ; mv Makefile.tmp Makefile.am"
	grep -v $s Makefile.am > Makefile.tmp ; mv Makefile.tmp Makefile.am
	cd ..
	rm src/$s.*
done

for d in $DIRS; do
	rm -Rf $d
done
