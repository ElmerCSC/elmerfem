SET(LIBSRC 
	huti_aux_C.F90 
	huti_aux_D.F90 
	huti_aux_S.F90 
	huti_aux_Z.F90 
	huti_bicgstab_2_C.F90 
	huti_bicgstab_2_D.F90 
	huti_bicgstab_2_S.F90 
	huti_bicgstab_2_Z.F90 
	huti_bicgstab_C.F90 
	huti_bicgstab_D.F90 
	huti_bicgstab_S.F90 
	huti_bicgstab_Z.F90 
	huti_cg_C.F90 
	huti_cg_D.F90 
	huti_cgs_C.F90 
	huti_cgs_D.F90 
	huti_cg_S.F90 
	huti_cgs_S.F90 
	huti_cgs_Z.F90 
	huti_cg_Z.F90 
	huti_fdefs.h 
	huti_gmres_C.F90 
	huti_gmres_D.F90 
	huti_gmres_S.F90 
	huti_gmres_Z.F90 
	huti_intdefs.h 
	huti_qmr_C.F90 
	huti_qmr_D.F90 
	huti_qmr_S.F90 
	huti_qmr_Z.F90 
	huti_tfqmr_C.F90 
	huti_tfqmr_D.F90 
	huti_tfqmr_S.F90 
	huti_defs.h 
	huti_main.c 
	huti_sfe_c_C.c 
	huti_sfe_c_D.c 
	huti_sfe_c_S.c 
	huti_sfe_c_Z.c 
	huti_tfqmr_Z.F90)


ADD_LIBRARY(huti SHARED ${LIBSRC})

INSTALL(TARGETS huti LIBRARY DESTINATION "lib")
