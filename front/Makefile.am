SUBDIRS=src include 

#
# Copy the tests directory 
#
dist-hook:
	cp -r src/lib/ $(distdir)/src/
	cp -r src/tcl/ $(distdir)/src/
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name CVS`
	rm -rf `find $(distdir)/ -name \#*`
	rm -rf `find $(distdir)/ -name \*~`
	rm -rf `find $(distdir)/ -name \*.so`


