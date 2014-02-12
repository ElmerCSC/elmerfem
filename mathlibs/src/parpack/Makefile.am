## Makefile.am for elmer's default libarpack.
COMMON_SRCS =  \
	pcgetv0.f\
	pclarnv.f\
	pcnaitr.f\
	pcnapps.f\
	pcnaup2.f\
	pcnaupd.f\
	pcneigh.f\
	pcneupd.f\
	pcngets.f\
	pdgetv0.f\
	pdlamch.f\
	pdlarnv.f\
	pdnaitr.f\
	pdnapps.f\
	pdnaup2.f\
	pdnaupd.f\
	pdneigh.f\
	pdneupd.f\
	pdngets.f\
	pdnorm2.f\
	pdsaitr.f\
	pdsapps.f\
	pdsaup2.f\
	pdsaupd.f\
	pdseigt.f\
	pdseupd.f\
	pdsgets.f\
	pdznorm2.f\
	pscnorm2.f\
	psgetv0.f\
	pslamch.f\
	pslarnv.f\
	psnaitr.f\
	psnapps.f\
	psnaup2.f\
	psnaupd.f\
	psneigh.f\
	psneupd.f\
	psngets.f\
	psnorm2.f\
	pssaitr.f\
	pssapps.f\
	pssaup2.f\
	pssaupd.f\
	psseigt.f\
	psseupd.f\
	pssgets.f\
	pzgetv0.f\
	pzlarnv.f\
	pznaitr.f\
	pznapps.f\
	pznaup2.f\
	pznaupd.f\
	pzneigh.f\
	pzneupd.f\
	pzngets.f\
	debug.h \
	pcmout.f \
	pcvout.f \
	pdmout.f \
	pdvout.f \
	pivout.f \
	psmout.f \
	psvout.f \
	pzmout.f \
	pzvout.f \
	stat.h 

EXTRA_DIST=mpif_stub.h

# disable parpack, if it MPI not found

if USE_MPI

lib_LIBRARIES = libparpack.a
libparpack_a_SOURCES = mpif $(COMMON_SRCS)
$(COMMON_SRCS): mpif

endif

if USE_LOCAL_MPIF_H
mpif:
	cp mpif_stub.h mpif.h
	cp mpif_stub.h mpif
else
mpif:
	cp mpif_stub.h mpif
endif


clean-local:
	$(RM) $(ALL) mpif mpif.h
