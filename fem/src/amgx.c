/* Copyright (c) 2011-2017, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#define MAX_MSG_LEN 4096

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "cuda_runtime.h"

/* CUDA error macro */
#define CUDA_SAFE_CALL(call) do {                                 \
  cudaError_t err = call;                                         \
  if(cudaSuccess != err) {                                        \
    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
            __FILE__, __LINE__, cudaGetErrorString( err) );       \
    exit(EXIT_FAILURE);                                           \
  } } while (0)

/* standard or dynamically load library */
#ifdef AMGX_DYNAMIC_LOADING
#include "amgx_capi.h"
#else
#include "amgx_c.h"
#endif

/* print error message and exit */
void errAndExit(const char *err)
{
    printf("%s\n", err);
    fflush(stdout);
    exit(1);
}

/* print callback (could be customized) */
void print_callback(const char *msg, int length)
{
    printf("%s", msg);
}

/* print usage and exit */
void printUsageAndExit()
{
    printf("%s", "Usage: ./amgx_capi [-mode [hDDI | hDFI | hFFI | dDDI | dDFI | dFFI]] [-m file] [-c config_file] [-amg \"variable1=value1 ... variable3=value3\"]\n");
    printf("%s", "     -mode:   select the solver mode\n");
    printf("%s", "     -m file: read matrix stored in the file\n");
    printf("%s", "     -c:      set the amg solver options from the config file\n");
    printf("%s", "     -amg:    set the amg solver options from the command line\n");
    exit(0);
}

/* parse parameters */
int findParamIndex(const char **argv, int argc, const char *parm)
{
    int count = 0;
    int index = -1;

    for (int i = 0; i < argc; i++)
    {
        if (strncmp(argv[i], parm, 100) == 0)
        {
            index = i;
            count++;
        }
    }

    if (count == 0 || count == 1)
    {
        return index;
    }
    else
    {
        printf("Error, parameter %s has been specified more than once, exiting\n", parm);
        exit(1);
    }

    return -1;
}

/* reade geometry (advanced input parmeters) */
void readGeometry( const char *fname, double **geo_x, double **geo_y, double **geo_z, int *dim, int *numrows)
{
    printf("Reading geometry from file: '%s'\n", fname);
    FILE *fin = fopen(fname, "r");

    if (!fin)
    {
        printf("Error opening file '%s'\n", fname);
        exit(1);
    }

    int n, dimension;

    if (2 != fscanf(fin, "%d %d\n", &n, &dimension))
    {
        errAndExit("Bad format\n");
    }

    *geo_x = (double *)malloc(n * sizeof(double));
    *geo_y = (double *)malloc(n * sizeof(double));

    if (dimension == 3)
    {
        *geo_y = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i ++)
            if (3 != fscanf(fin, "%lf %lf %lf\n", *geo_x + i, *geo_y + i, *geo_z + i))
            {
                errAndExit("Bad format\n");
            }
    }
    else if (dimension == 2)
    {
        for (int i = 0; i < n; i ++)
            if ( 2 != fscanf(fin, "%lf %lf\n", *geo_x + i, *geo_y + i))
            {
                errAndExit("Bad format\n");
            }
    }

    *dim = dimension;
    *numrows = n;
}

/* reade coloring (advanced input parmeters) */
void readColoring( const char *fname, int **row_coloring, int *colored_rows, int *num_colors)
{
    printf("Reading coloring from file: '%s'\n", fname);
    FILE *fin = fopen(fname, "r");

    if (!fin)
    {
        printf("Error opening file '%s'\n", fname);
        exit(1);
    }

    int n, colors_num;

    if (2 != fscanf(fin, "%d %d\n", &n, &colors_num))
    {
        errAndExit("Bad format\n");
    }

    *row_coloring = (int *)malloc(n * sizeof(int));

    for (int i = 0; i < n; i ++)
        if ( 1 != fscanf(fin, "%d\n", *row_coloring + i))
        {
            errAndExit("Bad format\n");
        }

    *colored_rows = n;
    *num_colors = colors_num;
}


int main(int argc, const char **argv)
{
    //parameter parsing
    int pidx = 0;
    int pidy = 0;
    //versions
    int major, minor;
    char *ver, *date, *time;
    //input geometry
    double *gx = NULL;
    double *gy = NULL;
    double *gz = NULL;
    //input coloring
    int dim = 0;
    int numrows = 0;
    int num_colors = 0;
    int colored_rows = 0;
    int *row_coloring = NULL;
    //input matrix and rhs/solution
    int n = 0;
    int bsize_x = 0;
    int bsize_y = 0;
    int sol_size = 0;
    int sol_bsize = 0;
    //library handles
    AMGX_Mode mode;
    AMGX_config_handle cfg;
    AMGX_resources_handle rsrc;
    AMGX_matrix_handle A;
    AMGX_vector_handle b, x;
    AMGX_solver_handle solver;
    //status handling
    AMGX_SOLVE_STATUS status;

    /* check arguments */
    if (argc == 1)
    {
        printUsageAndExit();
    }

    /* load the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    void *lib_handle = NULL;
    //open the library
#ifdef _WIN32
    lib_handle = amgx_libopen("amgxsh.dll");
#else
    lib_handle = amgx_libopen("libamgxsh.so");
#endif

    if (lib_handle == NULL)
    {
        errAndExit("ERROR: can not load the library");
    }

    //load all the routines
    if (amgx_liblink_all(lib_handle) == 0)
    {
        amgx_libclose(lib_handle);
        errAndExit("ERROR: corrupted library loaded\n");
    }

#endif
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());

    /* get api and build info */
    if ((pidx = findParamIndex(argv, argc, "--version")) != -1)
    {
        AMGX_get_api_version(&major, &minor);
        printf("amgx api version: %d.%d\n", major, minor);
        AMGX_get_build_info_strings(&ver, &date, &time);
        printf("amgx build version: %s\nBuild date and time: %s %s\n", ver, date, time);
        AMGX_SAFE_CALL(AMGX_finalize_plugins());
        AMGX_SAFE_CALL(AMGX_finalize());
        /* close the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
        amgx_libclose(lib_handle);
#endif
        exit(0);
    }

    /* get mode */
    if ((pidx = findParamIndex(argv, argc, "-mode")) != -1)
    {
        if (strncmp(argv[pidx + 1], "hDDI", 100) == 0)
        {
            mode = AMGX_mode_hDDI;
        }
        else if (strncmp(argv[pidx + 1], "hDFI", 100) == 0)
        {
            mode = AMGX_mode_hDFI;
        }
        else if (strncmp(argv[pidx + 1], "hFFI", 100) == 0)
        {
            mode = AMGX_mode_hFFI;
        }
        else if (strncmp(argv[pidx + 1], "dDDI", 100) == 0)
        {
            mode = AMGX_mode_dDDI;
        }
        else if (strncmp(argv[pidx + 1], "dDFI", 100) == 0)
        {
            mode = AMGX_mode_dDFI;
        }
        else if (strncmp(argv[pidx + 1], "dFFI", 100) == 0)
        {
            mode = AMGX_mode_dFFI;
        }
        else if (strncmp(argv[pidx + 1], "hCCI", 100) == 0)
        {
            mode = AMGX_mode_hZZI;
        }
        else if (strncmp(argv[pidx + 1], "hZCI", 100) == 0)
        {
            mode = AMGX_mode_hZCI;
        }
        else if (strncmp(argv[pidx + 1], "hZZI", 100) == 0)
        {
            mode = AMGX_mode_hZZI;
        }
        else if (strncmp(argv[pidx + 1], "dCCI", 100) == 0)
        {
            mode = AMGX_mode_dCCI;
        }
        else if (strncmp(argv[pidx + 1], "dZCI", 100) == 0)
        {
            mode = AMGX_mode_dZCI;
        }
        else if (strncmp(argv[pidx + 1], "dZZI", 100) == 0)
        {
            mode = AMGX_mode_dZZI;
        }
        else
        {
            errAndExit("ERROR: invalid mode");
        }
    }
    else
    {
        printf("Warning: No mode specified, using dDDI by default.\n");
        mode = AMGX_mode_dDDI;
    }

    /* create config */
    pidx = findParamIndex(argv, argc, "-amg");
    pidy = findParamIndex(argv, argc, "-c");

    if ((pidx != -1) && (pidy != -1))
    {
        printf("%s\n", argv[pidx + 1]);
        AMGX_SAFE_CALL(AMGX_config_create_from_file_and_string(&cfg, argv[pidy + 1], argv[pidx + 1]));
    }
    else if (pidy != -1)
    {
        AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, argv[pidy + 1]));
    }
    else if (pidx != -1)
    {
        printf("%s\n", argv[pidx + 1]);
        AMGX_SAFE_CALL(AMGX_config_create(&cfg, argv[pidx + 1]));
    }
    else
    {
        errAndExit("ERROR: no config was specified");
    }

    /* example of how to handle errors */
    //char msg[MAX_MSG_LEN];
    //AMGX_RC err_code = AMGX_resources_create_simple(NULL, cfg);
    //AMGX_SAFE_CALL(AMGX_get_error_string(err_code, msg, MAX_MSG_LEN));
    //printf("ERROR: %s\n",msg);
    /* switch on internal error handling (no need to use AMGX_SAFE_CALL after this point) */
    //AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));
    /* create resources, matrix, vector and solver */
    AMGX_resources_create_simple(&rsrc, cfg);
    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_vector_create(&x, rsrc, mode);
    AMGX_vector_create(&b, rsrc, mode);
    AMGX_solver_create(&solver, rsrc, mode, cfg);

    /* read the input system: matrix [and rhs & solution]
       Please refer to AMGX_read_system description in the AMGX_Reference.pdf
       manual for details on how to specify the rhs and the solution inside
       the input file. If these are not specified than rhs=[1,...,1]^T and
       (initial guess) sol=[0,...,0]^T. */
    if ((pidx = findParamIndex(argv, argc, "-m")) != -1)
    {
        AMGX_read_system(A, b, x, argv[pidx + 1]);
        AMGX_matrix_get_size(A, &n, &bsize_x, &bsize_y);
        AMGX_vector_get_size(x, &sol_size, &sol_bsize);

        if (sol_size == 0 || sol_bsize == 0)
        {
            AMGX_vector_set_zero(x, n, bsize_x);
        }
    }
    else
    {
        errAndExit("ERROR: no linear system was specified");
    }

// example of getting initial residual norm
    /*
      {
        void* t_norm = calloc(bsize_x, sizeof(double));
        AMGX_solver_calculate_residual_norm(solver, A, b, x, t_norm);
        printf("Initial norm: ");
        for (int i = 0; i < bsize_x; i++)
        {
          printf("%f ", ((double*)t_norm)[i]);
        }
        printf("\n");
        free(t_norm);
      }
    */

    /* read the input geometry */
    if ((pidx = findParamIndex(argv, argc, "-geo")) != -1)
    {
        readGeometry(argv[pidx + 1], &gx, &gy, &gz, &dim, &numrows);

        if (dim != 3) { gz = NULL; }

        if (dim != 0)
        {
            AMGX_matrix_attach_geometry(A, gx, gy, gz, numrows);
        }

        if (gx) { free(gx); }

        if (gy) { free(gy); }

        if (gz) { free(gz); }
    }

    /* read the input coloring */
    if ((pidx = findParamIndex(argv, argc, "-color")) != -1)
    {
        readColoring(argv[pidx + 1], &row_coloring, &colored_rows, &num_colors);

        if (num_colors != 0)
        {
            AMGX_matrix_attach_coloring( A, row_coloring, colored_rows, num_colors);
        }

        if (row_coloring) { free(row_coloring); }
    }

    /* solver setup */
    AMGX_solver_setup(solver, A);
    /* solver solve */
    AMGX_solver_solve(solver, b, x);
    /* example of how to change parameters between non-linear iterations */
    //AMGX_config_add_parameters(&cfg, "config_version=2, default:tolerance=1e-12");
    //AMGX_solver_solve(solver, b, x);
    AMGX_solver_get_status(solver, &status);
    /* example of how to print the residual history */
    //int nit;
    //double res;
    //AMGX_solver_get_iterations_number(solver, &nit);
    //for (int i=0; i<nit; i++) {
    //  printf("residual from iteration %d=", i);
    //  for (int j=0; j<bsize_y; j++) {
    //    AMGX_solver_get_iteration_residual(solver, i, j, &res);
    //    printf("%f ", (float)(res));
    //  }
    //  printf("\n");
    //}
    /* example of how to write the linear system to the output */
    //AMGX_write_system(A, b, x, "output.system.mtx");
    /* destroy resources, matrix, vector and solver */
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(x);
    AMGX_vector_destroy(b);
    AMGX_matrix_destroy(A);
    AMGX_resources_destroy(rsrc);
    /* destroy config (need to use AMGX_SAFE_CALL after this point) */
    AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
    /* shutdown and exit */
    AMGX_SAFE_CALL(AMGX_finalize_plugins());
    AMGX_SAFE_CALL(AMGX_finalize());
    /* close the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    amgx_libclose(lib_handle);
#endif
    //CUDA_SAFE_CALL(cudaDeviceReset());
    return status;
}
