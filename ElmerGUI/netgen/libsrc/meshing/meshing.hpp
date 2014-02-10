#ifndef FILE_MESHING
#define FILE_MESHING



#define CURVEDELEMS_NEW


#include "../include/myadt.hpp"
#include "../include/gprim.hpp"
#include "../include/linalg.hpp"
#include "../include/opti.hpp"




namespace netgen
{

  extern int printmessage_importance;

  class CSGeometry;
  
  
#include "msghandler.hpp"

#include "meshtype.hpp"
#include "localh.hpp"
#include "meshclass.hpp"
#include "global.hpp"


#include "meshtool.hpp"
#include "ruler2.hpp"
#include "adfront2.hpp"
#include "meshing2.hpp"
#include "improve2.hpp"


#include "geomsearch.hpp"
#include "adfront3.hpp"
#include "ruler3.hpp"

#ifndef SMALLLIB
#define _INCLUDE_MORE
#endif
#ifdef LINUX
#define _INCLUDE_MORE
#endif

#ifdef _INCLUDE_MORE
#include "meshing3.hpp"
#include "improve3.hpp"
#endif
#include "findip.hpp"
#include "findip2.hpp"

#include "topology.hpp"

#ifdef CURVEDELEMS_NEW
#include "curvedelems_new.hpp"
#else
#include "curvedelems.hpp"
#endif
#include "clusters.hpp"

#ifdef _INCLUDE_MORE
#include "meshfunc.hpp"
#endif
#include "bisect.hpp"
#include "hprefinement.hpp"
#include "boundarylayer.hpp"
#include "specials.hpp"

#include "validate.hpp"

#ifdef PARALLEL
#include "../parallel/paralleltop.hpp"
// #include "../parallel/parallelmesh.hpp"
#endif



}

#endif
