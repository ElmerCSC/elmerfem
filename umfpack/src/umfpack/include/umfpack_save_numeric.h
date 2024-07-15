/* ========================================================================== */
/* === umfpack_save_numeric ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  CISE Dept,   */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_save_numeric
(
    void *Numeric,
    char *filename
) ;

long umfpack_dl_save_numeric
(
    void *Numeric,
    char *filename
) ;

int umfpack_zi_save_numeric
(
    void *Numeric,
    char *filename
) ;

long umfpack_zl_save_numeric
(
    void *Numeric,
    char *filename
) ;

/*
double int Syntax:

    #include "umfpack.h"
    int status ;
    char *filename ;
    void *Numeric ;
    status = umfpack_di_save_numeric (Numeric, filename) ;

double long Syntax:

    #include "umfpack.h"
    long status ;
    char *filename ;
    void *Numeric ;
    status = umfpack_dl_save_numeric (Numeric, filename) ;

complex int Syntax:

    #include "umfpack.h"
    int status ;
    char *filename ;
    void *Numeric ;
    status = umfpack_zi_save_numeric (Numeric, filename) ;

complex long Syntax:

    #include "umfpack.h"
    long status ;
    char *filename ;
    void *Numeric ;
    status = umfpack_zl_save_numeric (Numeric, filename) ;

Purpose:

    Saves a Numeric object to a file, which can later be read by
    umfpack_*_load_numeric.  The Numeric object is not modified.

Returns:

    UMFPACK_OK if successful.
    UMFPACK_ERROR_invalid_Numeric_object if Numeric is not valid.
    UMFPACK_ERROR_file_IO if an I/O error occurred.

Arguments:

    void *Numeric ;	    Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric or loaded by umfpack_*_load_numeric.

    char *filename ;	    Input argument, not modified.

	A string that contains the filename to which the Numeric
	object is written.
*/
