#ifndef FILE_CSG
#define FILE_CSG

/* *************************************************************************/
/* File:   geoml.hpp                                                        */
/* Author: Joachim Schoeberl                                               */
/* Date:   21. Jun. 98                                                     */
/* *************************************************************************/

#include <myadt.hpp>
#include <gprim.hpp>
#include <meshing.hpp>

#include <geometry2d.hpp>

namespace netgen
{
#include "surface.hpp"
#include "solid.hpp"
#include "identify.hpp"
#include "singularref.hpp"
#include "csgeom.hpp"
#include "csgparser.hpp"

#ifndef SMALLLIB
#define _INCLUDE_MORE
#endif
//#ifdef LINUX
#define _INCLUDE_MORE
//#endif

#ifdef _INCLUDE_MORE
#include "triapprox.hpp"

#include "algprim.hpp"
#include "brick.hpp"
#include "spline3d.hpp"
#include "manifold.hpp"
#include "curve2d.hpp"
#include "explicitcurve2d.hpp"
#include "gencyl.hpp"
#include "polyhedra.hpp"
#include "extrusion.hpp"
#include "revolution.hpp"
#include "specpoin.hpp"
#include "edgeflw.hpp"
#include "meshsurf.hpp"
#endif
}

#endif
