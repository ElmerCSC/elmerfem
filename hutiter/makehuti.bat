
echo off
rem
rem A kludge to compile HUTIter with gmake and bash
rem

set OldPath=%Path%
set Path=..\util;%Path%

cd src
gmake -f Makefile.win libhuti
cd ..

gmake -f Makefile.win install

set Path=%OldPath%
set OldPath=
