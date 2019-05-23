/* femtypes.h */
/* Defines the types used in the FEM model. */

/* Definiotins used in allocating space for the structures. */
#define DIM 2               /* dimension of the space */
#define MAXDOFS 20          /* maximum number of variables, e.g. T,P */ 
#define MAXCELLS 100        /* maximum number of subcells in given direction */
#define MAXBOUNDARIES 50    /* maximum number of boundaries for BCs */
#define MAXMATERIALS  50    /* maximum index of materials */
#define MAXCASES    12      /* maximum number of coexisting cases */ 
#define MAXFILESIZE 600      /* maximum filenamesize for i/o files */
#define MAXLINESIZE 200     /* maximum length of line to be read */
#define MAXNAMESIZE 30      /* maximum size of the variablename */
#define MAXPARAMS 30        /* maximum number of parameters */
#define MAXVARS 20          /* maximum number of variables at the sides */
#define MAXNODESD2 27       /* maximum number of 2D nodes */ 
#define MAXNODESD1 9        /* maximum number of 1D nodes */
#define MAXMAPPINGS 10      /* maximum number of geometry mappings */
#define MAXCONNECTIONS 100  /* maximum number of connections in dual graph */
#define MAXBCS 1000         /* maximum number of BCs in naming */
#define MAXBODIES 100       /* maximum number of bodies in naming */
#define MAXPARTITIONS 512   /* maximum number of partitions */
#define MAXFORMATS 15

#define CONPLAIN 0
#define CONDISCONT 1
#define CONPERIODIC 2
#define CONCONSTRAINT 3

/* Struture GridType includes the subcell structure of the 
   geometry and the meshing information. The elements may be 
   directly derived from this structures but it takes some 
   time and is not easy to comprehend. Therefore structures 
   CellType and FemType are derived from this data. The special 
   subcell structure is, however, utilized in some mapping 
   subroutines that in general cases would be much more difficult 
   (and expensive) to perform.
   */
struct GridType {
  int dimension,
    triangles,
    layeredbc, 
    partitions,
    coordsystem,   /* 2D cartesian or axisymmetric? */ 
    layered,
    autoratio,     /* set the scale in x and y automatically? */
    minxelems,     /* minimum number of elements */
    minyelems,
    minzelems,
    totxelems,     /* total number of elements */
    totyelems,
    totzelems,
    elemorder,     
    elemmidpoints, 
    wantedelems,
    limitdxverify,
    firstmaterial, /* first material to be included in mesh */
    lastmaterial,  /* last material to be included in mesh */
    nocells,       /* number of subcells */
    xcells,        /* number of subcells in x-direction */
    ycells,
    zcells,
    noelements,    /* number of elements in the mesh */
    noknots,       /* number of knots in the mesh */
    nonodes,       /* number of nodes in one element */
    numbering,     /* numbering scheme */
    maxwidth,      /* maxwidth of the band matrix */
    noboundaries;  /* number of boundaries for BCs */
  int xlinear[MAXCELLS+1],    /* linearity flag within the subcells */
    ylinear[MAXCELLS+1],
    zlinear[MAXCELLS+1],
    xelems[MAXCELLS+1],       /* number of elements within subcells */
    yelems[MAXCELLS+1],
    zelems[MAXCELLS+1],
    zfirstmaterial[MAXCELLS+1], 
    zlastmaterial[MAXCELLS+1],  
    zmaterial[MAXCELLS+1],
    boundint[MAXBOUNDARIES],  /* internal material for boundary */
    boundext[MAXBOUNDARIES],  /* external material for boundary */
    boundsolid[MAXBOUNDARIES],/* which of these is the solid? */
    boundtype[MAXBOUNDARIES]; /* type of the boundary */
  int structure[MAXCELLS+2][MAXCELLS+2], /* material structure of subcells */
    numbered[MAXCELLS+2][MAXCELLS+2];    /* numbering order of the subcells */
  Real dx0,    /* global mesh scale in x-direction */
    dy0,
    dz0,
    limitdx, 
    triangleangle,
    xyratio, /* ratio between dx0 and dy0 */
    xzratio;
  Real rotateradius1,rotateradius2,rotateimprove;
  int rotate,rotateblocks,rotatecurve,rotatecartesian,mappings,
    reduceordermatmin,reduceordermatmax;
  Real curverad,curveangle,curvezet,polarradius;
  Real x[MAXCELLS+1],     /* vertical lines in the goemetry */
    y[MAXCELLS+1],        /* horizontal lines in the geometry */
    z[MAXCELLS+1],
    xexpand[MAXCELLS+1],  /* local expand ratio in the subcells */
    yexpand[MAXCELLS+1],
    zexpand[MAXCELLS+1],
    xratios[MAXCELLS+1],  /* relative mesh scale ratios in subcells */
    yratios[MAXCELLS+1],
    zratios[MAXCELLS+1],
    dx[MAXCELLS+1],       /* local mesh scale in the subcells */
    dy[MAXCELLS+1],
    dz[MAXCELLS+1],
    xdens[MAXCELLS+1],    /* local density of the mesh in the subcells */
    ydens[MAXCELLS+1],
    zdens[MAXCELLS+1];
  int mappingtype[MAXMAPPINGS],
    mappingline[MAXMAPPINGS],
    mappingpoints[MAXMAPPINGS];
  Real mappinglimits[2*MAXMAPPINGS],
    *mappingparams[MAXMAPPINGS];
};

/* The elements are numbered in the program without allocating 
   space for the knot numbers. Only a limited number of information 
   for each subcell is saved to structure CellType. Specific subroutines 
   are then used to calculate element or knot information using this 
   information. Cell is one macroscopic building block that may be 
   devided to M x N elements. It may even consist of one element. */
struct CellType {
  int nonodes,  /* number of nodes within an element */
    dimension,  /* 1D or 2D */
    numbering,  /* numbering scheme */
    xelem,      /* number of elements in the subcell */
    yelem,   
    levelwidth, /* width in knot numbering */ 
    left1st,    /* first index in the first line */
    left2nd,    /* first index in the second line */
    leftlast,   /* first index in the last line */
    levelwidthcenter,
    leftcenter, /* first index for 8 and 9-node elements */
    left2center,/* first index in the second line of 12- and 16-node elements */ 
    elem1st,    /* index of the lower left element */
    elemwidth,  /* width in element numbering */
    xlinear,    /* linearity flag */
    ylinear,
    material,   /* material flag */
    xind, yind; /* Indexes of the cell */
  int boundary[8], /* material indeces of neighbouring cells */
    neighbour[8];  /* number of neighbouring cells */
  Real xwidth,  /* size of the subcell */
    ywidth,
    xratio,     /* ratio of elements in the subcell */
    yratio,
    dx1,        /* local mesh scale */
    dy1;
  Real xcorner[4], /* coordinates of the subcell corners */
    ycorner[4];
};


/* This type includes all the element information needed for a 
   FEM model: the element topology, node coordinates, node indexing 
   and all the degrees of freedom. */
struct FemType {
  int created,     /* is the structure created? */
    noknots,       /* number of knots */
    noelements,    /* number of elements */
    coordsystem,   /* coordsystem flag */
    nocells,       /* number of subcells */
    maxnodes,      /* maximum number of nodes */
    dim,           /* dimension of space */
    variables,     /* number of variables */
    *dualgraph[MAXCONNECTIONS],  
    dualmaxconnections,
    indexwidth,
    dualexists,

    *partitiontable[MAXCONNECTIONS],  
    maxpartitiontable,
    partitiontableexists, 

    *invtopo[MAXCONNECTIONS],
    maxinvtopo,
    invtopoexists,
    timesteps,     /* number of timesteps */
    periodicexist, /* does the periodic vector exist? */
    *periodic,     /* peridic ordering vector, if needed */
    connectexist,  /* does the connection vector exist? */
    *connect,      /* connections between nodes, if needed */
    partitionexist,/* does the partitioning exist? */
    nopartitions,  /* number of partitions */
    *elempart,     /* which partition owns the element */
    *nodepart,     /* which partition owns the node */
    *elementtypes, /* types of elements using Elmer convention */
    *material,     /* material for each element */
    **topology,    /* element topology */
    bodynamesexist,
    boundarynamesexist;
  int edofs[MAXDOFS],   /* number of dofs in each node */
    alldofs[MAXDOFS];   /* total number of variables */
  Real minsize,maxsize;
  Real *x,  /* in axisymmetric case r */ 
      *y,   /* in axisymmetric case z */
      *z,   /* in cylindrical case theta */       
      *times;
  Real *dofs[MAXDOFS];  /* degrees of freedom in the mesh */
  char dofname[MAXDOFS][MAXNAMESIZE]; 
  char bodyname[MAXBODIES][MAXNAMESIZE];
  char boundaryname[MAXBCS][MAXNAMESIZE];
  int noboundaries,              /* number of boundaries */
      boundint[MAXBOUNDARIES],   /* internal material in the boundary */
      boundext[MAXBOUNDARIES],   /* external material in the boundary */
      boundsolid[MAXBOUNDARIES], /* which one is solid? */
      boundtype[MAXBOUNDARIES];  /* type of the boundary */
};

/* The boundaries between different materials or domains
   are saved into this structure. It is used for setting
   the boundary conditions. In physics it is typical that
   the BCs are more complicated than the equations in the 
   bulk and therefore the stucture must be such that it 
   enables the use of a wide variety of BCs. */
struct BoundaryType {
  int created,       /* is boundary created? */
    nosides,         /* sides on the boundary */
    maxsidenodes,  /* number of sidenodes on the element */
    fixedpoints,     /* number of fixed points allowed */
    coordsystem,     /* coordinate system flag */
    vfcreated,       /* are view factors created */
    gfcreated,       /* are Gephart factors created */
    maparea,         /* mappings of the area */
    mapvf,           /* mappings of the view factors */ 
    open,            /* is the closure partially open? */
    echain,          /* does the chain exist? */
    ediscont,        /* does the discontinuous boundary exist */
    chainsize;       /* size of the chain */ 
  int *parent,       /* primary parents of the sides */
    *parent2,        /* secondary parents of the sides */
    *material,       /* material of the sides */
    *side,           /* side in the primary parent element */
    *side2,          /* side in the secondary parent element */
    *chain,          /* indices in the chain representation */
    *types,
    *discont,        /* type of discontinuous and periodic BCs */
    *normal,         /* direction of the normal */
    *elementtypes,   /* side element types if needed */
    **topology,       /* topology if needed */
    points[MAXVARS], /* how many points for each side? */
    evars[MAXVARS];  /* does the variables exist? */
  Real totalarea,       /* total area of the side */
    areasexist,
    *areas,             /* side areas */
    **vf,               /* view factors */
    **gf,               /* Gephart factors */
    *vars[MAXVARS];     /* variables on the sides */
  char varname[MAXVARS][MAXNAMESIZE]; /* variable name */
};


#define MAXSIDEBULK 10
struct ElmergridType {

  int dim,
    silent,
    center,
    scale,      /* scale the geometry */
    order,      /* reorder the nodes */
    merge,      /* merge mesges */
    translate,  /* translate the mesh */
    rotate,     /* rotate the mesh */
    clone[3],   /* clone the mesh the number of given times */
    mirror[3],  /* mirror the mash around the given axis */
    canter, 
    decimals,   /* save the mesh with number of decimals */
    layers,     /* create boundary layers */
    layerbounds[MAXBOUNDARIES], 
    layernumber[MAXBOUNDARIES], 
    layermove,  /* map the created layer to the original geometry */
    metis,      /* number of Metis partitions */
    partopt,    /* free parameter for optimization */
    partitions, /* number of simple geometric partitions */
    partdim[3],
    inmethod,   /* method in which mesh is read in to ElmerGrid */
    outmethod,  /* method in which the mesh is written by ElmerGrid */
    sidemap[3*MAXBOUNDARIES],
    sidemappings,
    bulkmap[3*MAXMATERIALS],
    bulkmappings,
    boundorder, 
    bulkorder, 
    boundbounds,
    boundbound[3*MAXBOUNDARIES],
    bulkbounds,
    bulkbound[3*MAXBOUNDARIES], 
    mirrorbc,
    layerparents[MAXBOUNDARIES],
    sidebulk[MAXSIDEBULK],
    triangles,
    polar,
    usenames, 
    isoparam,
    cylinder,
    info,
    unitemeshes,
    reduce,
    removelowdim,
    removeunused,
    increase,
    reducemat1,
    reducemat2,
    findsides,
    saveboundaries,
    nodes3d,
    elements3d,
    periodic, 
    periodicdim[3],
    discont,
    discontbounds[MAXBOUNDARIES],
    connect,
    connectbounds[MAXBOUNDARIES],
    partorder,
    partitionhalo, /* create halo for the partitioning */
    partitionindirect, /* should one create indirect connections between nodes */
    nofilesin,
    saveinterval[3],
    elementsredone,
    pelemmap[4*MAXMATERIALS],pelems,
    belemmap[4*MAXMATERIALS], belems,
    advancedelem[7*MAXMATERIALS], advancedmat,
    bcoffset;

  Real cscale[3], 
    corder[3],
    cmerge,
    ctranslate[3],
    crotate[3],
    clonesize[3],
    layerratios[MAXBOUNDARIES], 
    layerthickness[MAXBOUNDARIES],
    layereps, 
    triangleangle, 
    partcorder[3],
    polarradius,
    relh;

  char filesin[MAXCASES][MAXFILESIZE],
    filesout[MAXCASES][MAXFILESIZE], 
    mapfile[MAXFILESIZE];
};

