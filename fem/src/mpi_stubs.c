#include "../config.h"

#ifndef HAVE_MPI_STUBS
void mpi_file_open_(){};
void mpi_file_sync_(){};
void mpi_file_close_(){};
void MPI_File_sync(){};
void MPI_Comm_dup(){};
void MPI_Type_extent(){};
void MPI_Initialized(){};
void MPI_File_write_at_all(){};
void MPI_File_open(){};
void MPI_File_close(){};
void MPI_Type_hindexed(){};
void MPI_Type_contiguous(){};
void MPI_File_read_at(){};
void MPI_Type_hvector(){};
void MPI_Error_string(){};
void MPI_Info_dup(){};
void MPI_Type_size(){};
void MPI_File_read_at_all(){};
void MPI_File_set_size(){};
void MPI_File_set_view(){};
void MPI_File_write_at(){};
void MPI_Get_elements(){};
void MPI_Info_free(){};
void MPI_Comm_free(){};

void MPI_Comm_test_inter() {}

void MPI_Comm_f2c() {}


void STDCALLBULL FC_FUNC_(mpi_init,MPI_INIT) 
     (int *p) { *p = 0; }
void STDCALLBULL FC_FUNC_(mpi_comm_size,MPI_COMM_SIZE) 
     (int *a, int *b, int *c) { *b = 1; *c = 0;}
void STDCALLBULL FC_FUNC_(mpi_comm_rank,MPI_COMM_RANK) 
     (int *a, int *b, int *c) { *b = 0; *c = 0;}
void STDCALLBULL FC_FUNC_(mpi_recv,MPI_RECV) 
     (int *a,int *b,int *c,int *d,int *e,int *f,int *g,int *h) {}
void STDCALLBULL FC_FUNC_(mpi_send,MPI_SEND)
     (int *a,int *b,int *c,int *d,int *e,int *f,int *g) {}
void STDCALLBULL FC_FUNC_(mpi_bcast,MPI_BCAST) () {}
void STDCALLBULL FC_FUNC_(mpi_barrier,MPI_BARRIER)
     (int *a,int *b) {}
void STDCALLBULL FC_FUNC_(mpi_finalize,MPI_FINALIZE)
     (int *a) {}
void STDCALLBULL FC_FUNC_(mpi_conversion_fn_null,MPI_CONVERSION_FN_NULL) () {}
void STDCALLBULL FC_FUNC_(mpi_dup_fn,MPI_DUP_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_null_copy_fn,MPI_NULL_COPY_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_buffer_detach,MPI_BUFFER_DETACH) ( void *buf, int *i, int *ierr ) {}
void STDCALLBULL FC_FUNC_(mpi_bsend,MPI_BSEND) (void *a, void *b, void *c, void *d, void *e, void *f, void *g ) {}
void STDCALLBULL FC_FUNC_(mpi_null_delete_fn,MPI_NULL_DELETE_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_buffer_attach,MPI_BUFFER_ATTACH) ( void *buf, int *i, int *ierr ) {}
void STDCALLBULL FC_FUNC_(mpi_allreduce,MPI_ALLREDUCE) () {}
void STDCALLBULL FC_FUNC_(mpi_wtime,MPI_WTIME) () {}
void STDCALLBULL FC_FUNC_(mpi_irecv,MPI_IRECV) () {}
void STDCALLBULL FC_FUNC_(mpi_isend,MPI_ISEND) () {}
void STDCALLBULL FC_FUNC_(mpi_test,MPI_TEST) () {}
void STDCALLBULL FC_FUNC_(mpi_cancel,MPI_CANCEL) () {}
void STDCALLBULL FC_FUNC_(mpi_ibsend,MPI_IBSEND) () {}

void STDCALLBULL FC_FUNC_(mpi_allgather,MPI_ALLGATHER) () {}

void STDCALLBULL FC_FUNC_(mpi_wtick,MPI_WTICK) () {}
void STDCALLBULL FC_FUNC_(pmpi_wtime,PMPI_WTIME) () {}
void STDCALLBULL FC_FUNC_(pmpi_wtick,PMPI_WTICK) () {}

void STDCALLBULL FC_FUNC_(mpi_type_null_copy_fn,MPI_TYPE_NULL_COPY_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_comm_dup_fn,MPI_COMM_DUP_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_type_null_delete_fn,MPI_TYPE_NULL_DELETE_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_win_dup_fn,MPI_WIN_DUP_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_win_null_copy_fn,MPI_WIL_NULL_COPY_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_type_dup_fn,MPI_TYPE_DUP_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_win_null_delete_fn,MPI_WIN_NULL_DELETE_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_comm_null_delete_fn,MPI_COMM_NULL_DELETE_FN) () {}
void STDCALLBULL FC_FUNC_(mpi_comm_null_copy_fn,MPI_COMM_NULL_DELETE_FN) () {}

void STDCALLBULL FC_FUNC_(mpi_waitall,MPI_WAITALL) () {}
void STDCALLBULL FC_FUNC_(mpi_waitany,MPI_WAITANY) () {}
void STDCALLBULL FC_FUNC_(mpi_wait,MPI_WAIT) () {}
void STDCALLBULL FC_FUNC_(mpi_win_start,MPI_WIN_START) () {}
void STDCALLBULL FC_FUNC_(mpi_win_test,MPI_WIN_TEST) () {}
void STDCALLBULL FC_FUNC_(mpi_win_complete,MPI_WIN_COMPLETE) () {}
void STDCALLBULL FC_FUNC_(mpi_testsome,MPI_TESTSOME) () {}
void STDCALLBULL FC_FUNC_(mpi_startall,MPI_STARTALL) () {}
void STDCALLBULL FC_FUNC_(mpi_start,MPI_START) () {}
void STDCALLBULL FC_FUNC_(mpi_testall,MPI_TESTALL) () {}
void STDCALLBULL FC_FUNC_(mpi_testany,MPI_TESTANY) () {}
void STDCALLBULL FC_FUNC_(mpi_win_post,MPI_WIN_POST) () {}
void STDCALLBULL FC_FUNC_(mpi_win_wait,MPI_WIN_WAIT) () {}
void STDCALLBULL FC_FUNC_(mpi_waitsome,MPI_WAITSOME) () {}

void STDCALLBULL FC_FUNC_(mpi_comm_create,MPI_COMM_CREATE) () {}
void STDCALLBULL FC_FUNC_(mpi_comm_free,MPI_COMM_FREE) () {}
void STDCALLBULL FC_FUNC_(mpi_comm_group,MPI_COMM_GROUP) () {}
void STDCALLBULL FC_FUNC_(mpi_group_incl,MPI_GROUP_INCL) () {}

void STDCALLBULL FC_FUNC(dmumps,DMUMPS) () {}

void STDCALLBULL FC_FUNC_(mpi_statuses_ignore,MPI_STATUSES_IGNORE) () {}

void STDCALLBULL MPI_Testall() {}
void STDCALLBULL MPI_Iprobe() {}
void STDCALLBULL MPI_Test() {}
void STDCALLBULL MPI_Type_free() {}
void STDCALLBULL MPI_Startall() {}
void STDCALLBULL MPI_Gather() {}
void STDCALLBULL MPI_Alltoall() {}
void STDCALLBULL MPI_Address() {}
void STDCALLBULL MPI_Allreduce() {}
void STDCALLBULL MPI_Comm_size() {}
void STDCALLBULL MPI_Recv() {}
void STDCALLBULL MPI_Isend() {}
void STDCALLBULL MPI_Allgather() {}
void STDCALLBULL MPI_Abort() {}
void STDCALLBULL MPI_Waitall() {}
void STDCALLBULL MPI_Wait() {}
void STDCALLBULL MPI_Comm_rank() {}
void STDCALLBULL MPI_Type_struct() {}
void STDCALLBULL MPI_Barrier() {}
void STDCALLBULL MPI_Allgatherv() {}
void STDCALLBULL MPI_Irecv() {}
void STDCALLBULL MPI_Bcast() {}
void STDCALLBULL MPI_Reduce() {}
void STDCALLBULL MPI_Request_free() {}
void STDCALLBULL MPI_Scatter() {}
void STDCALLBULL MPI_Type_commit() {}
void STDCALLBULL MPI_Recv_init() {}
void STDCALLBULL MPI_Wtime() {}
void STDCALLBULL MPI_Send_init() {}
void MPI_Probe() {}
void MPI_Send() {}
void MPI_Get_count() {}
void MPI_Type_vector() {}
void MPI_Irsend() {}

void STDCALLBULL FC_FUNC_(mpi_reduce_scatter,MPI_REDUCE_SCATTER)() {}
void STDCALLBULL FC_FUNC_(mpi_gather,MPI_GATHER)(){}
void STDCALLBULL FC_FUNC_(mpi_iprobe,MPI_IPROBE)(){}
void STDCALLBULL FC_FUNC_(mpi_reduce,MPI_REDUCE)(){}
void STDCALLBULL FC_FUNC_(mpi_op_create,MPI_OP_CREATE)(){}
void STDCALLBULL FC_FUNC_(mpi_request_free,MPI_REQUEST_FREE)(){}
void STDCALLBULL FC_FUNC_(mpi_unpack,MPI_UNPACK)(){}
void STDCALLBULL FC_FUNC_(mpi_pack_size,MPI_PACK_SIZE)(){}
void STDCALLBULL FC_FUNC_(mpi_ssend,MPI_SSEND)(){}
void STDCALLBULL FC_FUNC_(mpi_op_free,MPI_OP_FREE)(){}
void STDCALLBULL FC_FUNC_(mpi_alltoall,MPI_ALLTOALL)(){}
void STDCALLBULL FC_FUNC_(mpi_get_processor_name,MPI_GET_PROCESSOR_NAME)(){}
void STDCALLBULL FC_FUNC_(mpi_comm_dup,MPI_COMM_DUP)(){}
void STDCALLBULL FC_FUNC_(mpi_abort,MPI_ABORT)(){}
void STDCALLBULL FC_FUNC_(mpi_pack,MPI_PACK)(){}
void STDCALLBULL FC_FUNC_(mpi_get_count,MPI_GET_COUNT)(){}
void STDCALLBULL FC_FUNC_(mpi_comm_split,MPI_COMM_SPLIT)(){}
void STDCALLBULL FC_FUNC_(mpi_probe,MPI_PROBE)(){}
void STDCALLBULL FC_FUNC_(mpi_gatherv,MPI_GATHERV)(){}
void STDCALLBULL FC_FUNC_(mpi_initialized,MPI_INITIALIZED)(){}
void MPI_Group_free(){}
void MPI_Error_class(){}
void MPI_Attr_get() {}
void MPI_Comm_create(){}
void MPI_Op_free(){}
void MPI_Rsend(){}
void MPI_Sendrecv(){}
void MPI_Op_create(){}
void MPI_Comm_group(){}
void MPI_Comm_split(){}
void MPI_Pack_size(){}
void MPI_Group_incl(){}
void MPI_Pack(){}
void MPI_Init(){}

void STDCALLBULL FC_FUNC_(mpi_init_thread,MPI_INIT_THREAD)(int *r, int *p, int *ierr) {*p = *r;}

#endif

void hpmp_flinteroperate() {}
void hpmp_max() {}
void hpmp_f_mpi_bottom() {}
void hpmp_min() {}
void hpmp_sum() {}
void hpmp_int() {}
void hpmp_double() {}
void hpmp_comm_world() {}

void hpmp_long() {}
void hpmp_lor() {}

/* parpack */
void STDCALLBULL FC_FUNC(pdseupd,PDSEUPD) ( void *a, void *b, void *c, void *d, void *e, void *f, void *g, void *a8, void *b9, void *c10, void *d11, void *e12, void *f13, void *g14, void *a15, void *b16, void *c17, void *d18, void *e19, void *a20, void *b21, void *c22, void *d23, void *e24, void *f25, void *g26  ) {}
void STDCALLBULL FC_FUNC(pdneupd,PDNEUPD) ( void *a, void *b, void *c, void *d, void *e, void *f, void *g, void *a8, void *b9, void *c10, void *d11, void *e12, void *f13, void *g14, void *a15, void *b16, void *c17, void *d18, void *e19, void *a20, void *b21, void *c22, void *d23, void *e24, void *f25, void *g26, void *g27, void *g28,void *g29 ) {}
void STDCALLBULL FC_FUNC(pzneupd,PZNEUPD) ( void *a, void *b, void *c, void *d, void *e, void *f, void *g, void *a8, void *b9, void *c10, void *d11, void *e12, void *f13, void *g14, void *a15, void *b16, void *c17, void *d18, void *e19 ) {}

void STDCALLBULL FC_FUNC(pdsaupd,PDSAUPD) ( void *a, void *b, void *c, void *d, void *e, void *f, void *g, void *a8, void *b9, void *c10, void *d11, void *e12, void *f13, void *g14, void *a15, void *b16, void *c17, void *d18, void *e19 ) {}
void STDCALLBULL FC_FUNC(pdnaupd,PDNAUPD) ( void *a, void *b, void *c, void *d, void *e, void *f, void *g, void *a8, void *b9, void *c10, void *d11, void *e12, void *f13, void *g14, void *a15, void *b16, void *c17, void *d18, void *e19 ) {}
void STDCALLBULL FC_FUNC(pznaupd,PZNAUPD) ( void *a, void *b, void *c, void *d, void *e, void *f, void *g, void *a8, void *b9, void *c10, void *d11, void *e12, void *f13, void *g14, void *a15, void *b16, void *c17, void *d18, void *e19 ) {}



