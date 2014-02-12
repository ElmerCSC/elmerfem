/*

  ElmerParam - A simple system for parametrized computing
 
  Copyright (C) 2006  CSC - IT Center for Science Ltd.

  Authors: Erik Edelmann <Erik.Edelmann@csc.fi>
           Peter Råback <Peter.Raback@csc.fi>
  Address: CSC - IT Center for Science Ltd.
           Keilaranta 14
           02101 Espoo, Finland
            
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program (in file elmerparam/COPYING); if not, write to
  the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
  Boston, MA 02110-1301, USA.

 */
#ifndef _DYNARRAY_H_
#define _DYNARRAY_H_

/* This module defines arrays of ints (daint_t) and reals (dareal_t) that can
 * grow dynamically when new elements are assigned.  Note that it's the types
 * daint_t and dareal_t and the accompanying functions di_* and dr_* the caller
 * is supposed to use -- dynarray_t and the dynarray_* functions are for
 * internal use.
 *
 * The correct way of using da*_t array is like this:
 *
 * daint_t *x = NULL;  // Shall begin it's life NULLified.
 * ...
 * x = di_set(x, i, k);  // Assign the value of di_set() back to x (important!)
 * ...
 * printf("%i\n", di_get(x, i));  // No suprises here.
 * ...
 * da_kill(x); // Using just free() is not enough!
 */

#define DYNARRAY_ALEN 100

typedef union {
    int i;
    double r;
} da_numeric_t;

typedef struct dynarray_t {
    size_t n; /* Number of elements in the (rest of) the dynarray_t.  */
    da_numeric_t a[DYNARRAY_ALEN];
    struct dynarray_t *next;
} dynarray_t;


typedef struct {
    dynarray_t da;
} daint_t;

typedef struct {
    dynarray_t da;
} dareal_t;


dynarray_t *dynarray_set(dynarray_t *da, int i, da_numeric_t val);
dynarray_t *dynarray_set_from_matc(dynarray_t *da, char type, const char *var);
da_numeric_t dynarray_get(dynarray_t *da, int i);
void dynarray_kill(dynarray_t *da);


/* Assign 'val' to the i:th element of '*di'. 'di' shall be either NULL, or
 * returned from a previous call to this function. */

static inline daint_t *di_set(daint_t *di, int i, int val) 
{
    da_numeric_t u;
    u.i = val;
    return (daint_t *)dynarray_set((dynarray_t *)di, i, u);
}


/* Assign 'val' to the i:th element of '*dr'. 'dr' shall be either NULL, or
 * returned from a previous call to this function. */

static inline dareal_t *dr_set(dareal_t *dr, int i, double val)
{
    da_numeric_t u;
    u.r = val;
    return (dareal_t *)dynarray_set((dynarray_t *)dr, i, u);
}


/* Set 'dr' to the MATC variable 'var'.  */

static inline dareal_t *dr_set_from_matc(dareal_t *dr, const char *var)
{
    return (dareal_t *)dynarray_set_from_matc((dynarray_t *)dr, 'r', var);
}

/* Set 'dr' to the MATC variable 'var'.  */

static inline daint_t *di_set_from_matc(daint_t *di, const char *var)
{
    return (daint_t *)dynarray_set_from_matc((dynarray_t *)di, 'i', var);
}

/* Get the value of i:th element of '*di'.  If it hasn't been assigned (using
 * di_set(), it's value is undefined (not necessarily 0!)  */

static inline int di_get(daint_t *di, int i)
{
    return dynarray_get((dynarray_t *)di, i).i;
}


/* Get the value of i:th element of '*dr'.   If it hasn't been assigned (using
 * dr_set(), it's value is undefined (not necessarily 0.0!)  */

static inline double dr_get(dareal_t *dr, int i)
{
    return dynarray_get((dynarray_t *)dr, i).r;
}


/* Free all space used by a daint_t or dareal_t (works for both types). */

static inline void da_kill(void *d)
{
    dynarray_kill((dynarray_t *)d);
}


/* Return number of elements in 'd' (works for both daint_t and dareal_t).  */

static inline size_t da_n(void *d)
{
    return d ? ((dynarray_t *)d)->n : 0;
}

#endif /* _DYNARRAY_H_ */
