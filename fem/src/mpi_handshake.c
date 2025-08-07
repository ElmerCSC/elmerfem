/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://gitlab.dkrz.de/dkrz-sw/mpi-handshake
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

//taken from http://beige.ucs.indiana.edu/I590/node85.html
static void mh_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}

static void mh_abort(
  const char * error_string, MPI_Comm comm, const char * file, int line) {

  int rank;
  MPI_Comm_rank(comm, &rank);

  fprintf(stderr, "%3d: %s\n", rank, error_string);
  fprintf(stderr, "Aborting in file %s, line %i ...\n", file, line );
  MPI_Abort(comm, EXIT_FAILURE);
}

#define mh_mpi_call(call, comm)               \
  do {                                          \
    int error_code = (call);                    \
    if (error_code != MPI_SUCCESS)              \
      mh_mpi_error(error_code, comm);          \
  } while(0)

#define mh_assert(exp, msg, comm) \
  {if(!((exp))) mh_abort(((msg)), ((comm)), __FILE__, __LINE__);}

void mpi_handshake(
  char const ** group_names, MPI_Comm * group_comms, size_t n, MPI_Comm comm) {

  // check whether comm is an intercomm
  int is_intercomm;
  mh_mpi_call(MPI_Comm_test_inter(comm, &is_intercomm), comm);
  mh_assert(
    !is_intercomm,
    "ERROR(mpi_handshake): inter-communicators are not supported", comm);

  // STEP 1: Version exchange
  enum {MPI_HANDSHAKE_VERSION = 1};
  int mh_version = MPI_HANDSHAKE_VERSION;
  mh_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &mh_version, 1, MPI_INT, MPI_MIN, comm), comm);
  mh_assert(
    mh_version == MPI_HANDSHAKE_VERSION,
      "ERROR(mpi_handshake): version mismatch."
      "This implementation of the MPI handshake only supports version 1", comm);

  // generate communicators
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  for (size_t i = 0; i < n; ++i) group_comms[i] = MPI_COMM_NULL;

  while(1){
    // STEP 2: determine broadcasting rank
    size_t group_idx = SIZE_MAX;
    for (size_t i = 0; (i < n) && (group_idx == SIZE_MAX); ++i)
      if (group_comms[i] == MPI_COMM_NULL) group_idx = i;
    int broadcasting_rank = group_idx != SIZE_MAX ? rank : size;
    mh_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, &broadcasting_rank, 1, MPI_INT, MPI_MIN, comm), comm);
    mh_assert(broadcasting_rank >= 0 && broadcasting_rank <= size,
      "ERROR(mpi_handshake): "
      "broadcasting rank cannot be negativ or greater than communicator size.",
              comm);
    if (broadcasting_rank == size) break;

    // STEP 3: broadcast group name
    int group_name_buffer_size  = 0;
    if(broadcasting_rank == rank){
      size_t len = strlen(group_names[group_idx]);
      mh_assert(len <= INT_MAX,
        "ERROR(yac_mpi_handshake): group name is too long", comm);
      group_name_buffer_size = (int)len;
    }
    MPI_Bcast(&group_name_buffer_size, 1, MPI_INT, broadcasting_rank, comm);
    char * group_name_buffer =
      malloc((group_name_buffer_size+1) * sizeof(*group_name_buffer));
    mh_assert(
      group_name_buffer,
      "ERROR(mpi_handshake): failed to allocate group name buffer", comm);
    if (broadcasting_rank == rank)
      strcpy(group_name_buffer, group_names[group_idx]);
    mh_mpi_call(
      MPI_Bcast(
        group_name_buffer, group_name_buffer_size, MPI_CHAR, broadcasting_rank, comm), comm);
    group_name_buffer[group_name_buffer_size] = '\0';

    // STEP 4: split communicator
    group_idx = SIZE_MAX;
    for (size_t i = 0; (i < n) && (group_idx == SIZE_MAX); ++i)
      if (!strcmp(group_name_buffer, group_names[i])){
        mh_assert(group_comms[i] == MPI_COMM_NULL,
          "ERROR(yac_mpi_handshake): "
          "Group communicator for a group that was already created "
          "was broadcasted again.", comm);
        group_idx = i;
      }
    free(group_name_buffer);
    MPI_Comm group_comm;
    mh_mpi_call(
      MPI_Comm_split(
        comm, (group_idx != SIZE_MAX)?0:MPI_UNDEFINED, rank, &group_comm),
      comm);
    if (group_idx != SIZE_MAX) {
      group_comms[group_idx] = group_comm;
    }
  }
}

void mpi_handshake_c2f(int n, char const ** group_names,
                       MPI_Fint * group_comms, MPI_Fint comm)
{
  MPI_Comm comm_c = MPI_Comm_f2c(comm);
  MPI_Comm * group_comms_c = malloc(n*sizeof(*group_comms_c));
  mpi_handshake(group_names, group_comms_c, n, comm_c);
  for(int i = 0; i<n; ++i)
    group_comms[i] = MPI_Comm_c2f(group_comms_c[i]);
  free(group_comms_c);
}

void mpi_handshake_dummy(MPI_Comm comm) {mpi_handshake(NULL, NULL, 0, comm);}
