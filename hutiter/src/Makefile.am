## Makefile.am for hutiter/src

INCLUDES = \
	-I$(top_builddir)/include       \
	-I$(top_srcdir)/include         

COMMON_SRCS =  \
	huti_aux_C.F90 \
	huti_aux_D.F90 \
	huti_aux_S.F90 \
	huti_aux_Z.F90 \
	huti_bicgstab_2_C.F90 \
	huti_bicgstab_2_D.F90 \
	huti_bicgstab_2_S.F90 \
	huti_bicgstab_2_Z.F90 \
	huti_bicgstab_C.F90 \
	huti_bicgstab_D.F90 \
	huti_bicgstab_S.F90 \
	huti_bicgstab_Z.F90 \
	huti_cg_C.F90 \
	huti_cg_D.F90 \
	huti_cgs_C.F90 \
	huti_cgs_D.F90 \
	huti_cg_S.F90 \
	huti_cgs_S.F90 \
	huti_cgs_Z.F90 \
	huti_cg_Z.F90 \
	huti_fdefs.h \
	huti_gmres_C.F90 \
	huti_gmres_D.F90 \
	huti_gmres_S.F90 \
	huti_gmres_Z.F90 \
	huti_intdefs.h \
	huti_qmr_C.F90 \
	huti_qmr_D.F90 \
	huti_qmr_S.F90 \
	huti_qmr_Z.F90 \
	huti_tfqmr_C.F90 \
	huti_tfqmr_D.F90 \
	huti_tfqmr_S.F90 \
	huti_defs.h \
	huti_main.c \
	huti_sfe_c_C.c \
	huti_sfe_c_D.c \
	huti_sfe_c_S.c \
	huti_sfe_c_Z.c \
	huti_tfqmr_Z.F90




lib_LIBRARIES = libhuti.a
libhuti_a_SOURCES = $(COMMON_SRCS)

#################################################################################
# If sources are missing generate them with the preprocessor.
# this is the section to watch out for, when compiling on more exotic (non-linux) 
# platforms.
#
#
HUTI_CPP=$(CPP)
HUTI_PPFLAGS = -traditional -P
# and now for some sed generated makefile...
huti_aux_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_aux.src huti_aux_C.F90
huti_aux_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_aux.src huti_aux_D.F90
huti_aux_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_aux.src huti_aux_S.F90
huti_aux_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_aux.src huti_aux_Z.F90
huti_bicgstab_2_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_bicgstab_2.src huti_bicgstab_2_C.F90
huti_bicgstab_2_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_bicgstab_2.src huti_bicgstab_2_D.F90
huti_bicgstab_2_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_bicgstab_2.src huti_bicgstab_2_S.F90
huti_bicgstab_2_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_bicgstab_2.src huti_bicgstab_2_Z.F90
huti_bicgstab_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_bicgstab.src huti_bicgstab_C.F90
huti_bicgstab_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_bicgstab.src huti_bicgstab_D.F90
huti_bicgstab_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_bicgstab.src huti_bicgstab_S.F90
huti_bicgstab_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_bicgstab.src huti_bicgstab_Z.F90
huti_cg_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_cg.src huti_cg_C.F90
huti_cg_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_cg.src huti_cg_D.F90
huti_cgs_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_cgs.src huti_cgs_C.F90
huti_cgs_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_cgs.src huti_cgs_D.F90
huti_cg_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_cg.src huti_cg_S.F90
huti_cgs_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_cgs.src huti_cgs_S.F90
huti_cgs_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_cgs.src huti_cgs_Z.F90
huti_cg_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_cg.src huti_cg_Z.F90
huti_gmres_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_gmres.src huti_gmres_C.F90
huti_gmres_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_gmres.src huti_gmres_D.F90
huti_gmres_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_gmres.src huti_gmres_S.F90
huti_gmres_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_gmres.src huti_gmres_Z.F90
huti_qmr_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_qmr.src huti_qmr_C.F90
huti_qmr_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_qmr.src huti_qmr_D.F90
huti_qmr_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_qmr.src huti_qmr_S.F90
huti_qmr_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_qmr.src huti_qmr_Z.F90
huti_tfqmr_C.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_tfqmr.src huti_tfqmr_C.F90
huti_tfqmr_D.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_tfqmr.src huti_tfqmr_D.F90
huti_tfqmr_S.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_tfqmr.src huti_tfqmr_S.F90
huti_tfqmr_Z.F90:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_tfqmr.src huti_tfqmr_Z.F90
huti_sfe_c_C.c:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DC_PRE huti_sfe_c.src huti_sfe_c_C.c
huti_sfe_c_D.c:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DD_PRE huti_sfe_c.src huti_sfe_c_D.c
huti_sfe_c_S.c:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DS_PRE huti_sfe_c.src huti_sfe_c_S.c
huti_sfe_c_Z.c:
	$(HUTI_CPP) $(HUTI_PPFLAGS) -DZ_PRE huti_sfe_c.src huti_sfe_c_Z.c

##########################################################
# The old for loop (which screwed up the dependencies, but was nicer anyway)
#
# $(PP_SRCS): preprocess
#PRECS	=	S D C Z
#precdf90=huti_cg huti_tfqmr huti_cgs huti_qmr huti_bicgstab \
#	 huti_gmres huti_bicgstab_2 huti_aux
#precdc=huti_sfe_c
#
#
#preprocess:
#	@echo
#	@echo 'Generating source files from *.src files, WARNING, only known to work on Linux :)'
#	@echo
#
# First make the precision dependent Fortran90 files
#
#	@if [ "$(precdf90)X" != "X" ]; then \
#	  for p in $(PRECS); do \
#	    for pd in $(precdf90); do \
#		$(RM) $${pd}_$${p}.F90; \
#		echo $(HUTI_CPP) $(HUTI_PPFLAGS) -D$${p}_PRE $${pd}.src $${pd}_$${p}.F90; \
#		$(HUTI_CPP) $(HUTI_PPFLAGS) -D$${p}_PRE $${pd}.src $${pd}_$${p}.F90; \
#	    done; \
#	  done; \
#	fi 
#
# Third make the precision dependent C files
#
#	@if [ "$(precdc)X" != "X" ]; then \
#	  for p in $(PRECS); do \
#	    for pd in $(precdc); do \
#		$(RM) $${pd}_$${p}.c; \
#		$(HUTI_CPP)  $(HUTI_PPFLAGS) -D$${p}_PRE $${pd}.src $${pd}_$${p}.c; \
#	    done; \
#	  done; \
#	fi 


# if LOCAL_BLAS_AND_LAPACK
#libhuti_a_LIBADD = \
#        ../blas/libblas.a \
#        ../lapack/liblapack.a
#else 
#if LOCAL_BLAS
#libhuti_a_LIBADD = \
#	../blas/libblas.a
#else 
#if LOCAL_LAPACK
#libhuti_a_LIBADD = \
#	../lapack/liblapack.a
#
#endif
#endif
#endif

EXTRA_DIST=\
	huti_aux.src \
	huti_bicgstab.src \
	huti_bicgstab_2.src \
	huti_cg.src \
	huti_cgs.src \
	huti_gmres.src \
	huti_qmr.src \
	huti_sfe_c.src \
	huti_tfqmr.src 


# stupid hack around shitty visual fortran CPP
if USING_MS_C_COMPILER
.F90.obj:
	rm -f tmp.f90
	$(CPP) $< > tmp.f90
	$(PPFCCOMPILE) -c -o $@ tmp.f90
endif
