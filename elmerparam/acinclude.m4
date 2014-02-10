AC_DEFUN([ee_FC_MODFILE],
[AC_CACHE_CHECK(for *.mod-file naming convention, modfilename,
[
mkdir conftestdir
cd conftestdir
cat > foo.f90 << STOP
module modulename
    integer :: bar
end module modulename
STOP
${FC-fc} -c foo.f90
modfilename=""
for m in modulename.mod MODULENAME.mod; do
    if test -r $m ; then
        modfilename=$m
        break;
    fi
done
cd ..
rm -rf conftestdir
if test $modsearch_opt = "" ; then
    AC_MSG_WARN(Can't find *.mod-file naming conventions for ${FC-fc})
fi
])
])


AC_DEFUN([ee_FC_MODSEARCH_OPT],
[AC_CACHE_CHECK(for option needed by ${FC-fc} to find non-local *.mod-files,
modsearch_opt,
[
mkdir conftestdir
cd conftestdir
mkdir prefix
mkdir prefix/include
mkdir src
cd src
cat > foo.f90 << STOP
module foo
    integer :: bar
end module foo
STOP
${FC-fc} -c foo.f90
mv *.mod ../prefix/include
cat > foouser.f90 << STOP
program foouser
    use foo
end program foouser
STOP
modsearch_opt=""
for opt in -I -M -B; do
    if ${FC-fc} ${opt}../prefix/include foouser.f90 > /dev/null 2>&1 ; then
        modsearch_opt=$opt
        break;
    fi
done
cd ../..
rm -rf conftestdir
if test $modsearch_opt = "" ; then
    AC_MSG_WARN(Can't find *.mod-file search option for ${FC-fc})
fi
])
AC_SUBST(modsearch_opt)
])

AC_DEFUN([ee_MEXEXT], [
AC_CACHE_CHECK(for suffix of Matlab/MEX files, mex_ext, [
mkdir conftestdir
cd conftestdir
cat > test.c << STOP
#include "mex.h"

void mexFunction(int nlhs, mxArray **plhs[], int nrhs, const mxArray **prhs[])
{
        return;
}
STOP
export MATLAB_HOME="/v/solaris9/appl/math/matlab/matlab7"
mex test.c
mex_ext=`ls test.mex* | sed 's/test\.//' | sed 's/\*//'`
cd ..
rm -rf conftestdir
])
AC_SUBST(mex_ext)
])


AC_DEFUN([ee_IARGC_DECL], [
AC_CACHE_CHECK(for iargc() declaration needed in Fortran code,
iargc_decl,
[
mkdir conftestdir
cd conftestdir
cat > test.f90 << STOP
program test
    implicit none
    print *, iargc()
end program test
STOP
if ${FC-fc} test.f90 > /dev/null 2>&1; then
    iargc_decl="!none"
else
    cat > test.f90 << STOP
program test
    implicit none
    integer, external :: iargc
    print *, iargc()
end program test
STOP
    if ${FC-fc} test.f90 > /dev/null 2>&1; then
        iargc_decl="integer, external :: iargc"
    else
        AC_MSG_WARN("Fortran function iargc doesn't seem to be supported")
    fi
fi
cd ..
rm -rf conftestdir
AC_SUBST(iargc_decl)
])
])
