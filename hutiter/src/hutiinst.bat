
echo off
set OldPath=%Path%
set Path=..\util;%Path%

gmake -f Makefile.win install

set Path=%OldPath%

echo on

