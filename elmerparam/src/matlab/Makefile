all: elmer_param.

elmer_param.: elmer_param.c
	mex elmer_param.c -I.. -L.. -lelmerparam -lm

install: elmer_param.
	cp elmer_param. /usr/local/lib

distdir: clean
	cp -r ./* $(distdir)/

clean:
	rm -f elmer_param.

check:
