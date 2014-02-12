all:

install:
	R CMD build elmerparam
	R CMD INSTALL -l /usr/local/lib/R elmerparam

distdir: clean
	cp -r ./* $(distdir)/

clean:
	rm -f elmerparam/src/*.{o,so}

check:
