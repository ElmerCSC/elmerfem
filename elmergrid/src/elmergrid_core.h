/* elmergrid_core.h
 *
 * Aggregate header for the ElmerGrid core routines that are shared between the
 * ElmerGrid command line tool and the ElmerGrid plugin built into ElmerGUI.
 *
 * The core is written in C and is compiled as C.  ElmerGUI is C++, so the
 * declarations are wrapped in an extern "C" block here rather than in each of
 * the individual headers, which carry no include guards of their own.
 *
 * C++ callers should include this header instead of the individual egXXX.h
 * headers.  C callers may continue to include the individual headers directly.
 */

#ifndef ELMERGRID_CORE_H
#define ELMERGRID_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egmesh.h"
#include "egnative.h"
#include "egconvert.h"

#ifdef __cplusplus
}
#endif

#endif /* ELMERGRID_CORE_H */
