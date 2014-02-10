SUBDIRS=src tests 

EXTRA_DIST=autogen.sh acx_elmer.m4 acx_blas.m4 acx_lapack.m4 acx_mpi.m4 GPL-2 LICENSES

#
# Copy the tests directory 
#
dist-hook:
	cp -r tests $(distdir)
	cp -r examples $(distdir)
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name CVS`
	rm -rf `find $(distdir)/ -name \#*`
	rm -rf `find $(distdir)/ -name \*~`
	rm -rf `find $(distdir)/ -name \*.so`
#	if test `ls nodist.sh` != ""; then     \
#	  sh nodist.sh;                        \
#	fi

CLEANFILES = a.out
