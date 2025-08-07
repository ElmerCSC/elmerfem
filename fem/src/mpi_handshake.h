/**
 * dummy call of Common MPI communicator splitting
 */
void mpi_handshake_dummy(MPI_Comm comm);

/**
 * Common MPI communicator splitting
 *
 * @param[in]  group_names names of process groups
 * @param[out] group_comms MPI communicators containing all processes
 *                         that provided the same group names
 * @param[in]  n           number of group names
 * @param[in]  comm        MPI communicator used for the splitting
 * @remark this call is collective for all processes in comm (if no splitting
 *         is required by the local processes \ref mpi_handshake_dummy can be
 *         called instead)
 */
void mpi_handshake(
  char const ** group_names, MPI_Comm * group_comms, size_t n, MPI_Comm comm);
