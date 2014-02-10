
SET(parpack_SRCS pcgetv0.f pclarnv.f pcnaitr.f pcnapps.f pcnaup2.f
	pcnaupd.f pcneigh.f pcneupd.f pcngets.f pdgetv0.f pdlamch.f
	pdlarnv.f pdnaitr.f pdnapps.f pdnaup2.f pdnaupd.f pdneigh.f
	pdneupd.f pdngets.f pdnorm2.f pdsaitr.f pdsapps.f pdsaup2.f
	pdsaupd.f pdseigt.f pdseupd.f pdsgets.f pdznorm2.f pscnorm2.f
	psgetv0.f pslamch.f pslarnv.f psnaitr.f psnapps.f psnaup2.f
	psnaupd.f psneigh.f psneupd.f psngets.f psnorm2.f pssaitr.f
	pssapps.f pssaup2.f pssaupd.f psseigt.f psseupd.f pssgets.f
	pzgetv0.f pzlarnv.f pznaitr.f pznapps.f pznaup2.f pznaupd.f
	pzneigh.f pzneupd.f pzngets.f debug.h pcmout.f pcvout.f
	pdmout.f pdvout.f pivout.f psmout.f psvout.f pzmout.f pzvout.f
	stat.h)

ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
LINK_DIRECTORIES(${MPI_Fortran_LIBRARIES})

ADD_LIBRARY(parpack SHARED ${parpack_SRCS})

ADD_DEPENDENCIES(parpack arpack)

INSTALL(TARGETS parpack DESTINATION "lib")
