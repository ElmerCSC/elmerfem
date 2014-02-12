
  // HP_PYRAMID
  int refpyramid_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refpyramid_newelstypes[] =
    {
      HP_PYRAMID,
      HP_NONE,
    };
  int refpyramid_newels[][8] =
    {
      { 1, 2, 3, 4, 5 }
    };
  HPRef_Struct refpyramid =
    {
      HP_PYRAMID,
      refpyramid_splitedges, 
      0, 0,
      refpyramid_newelstypes, 
      refpyramid_newels
    };


// singular point 1      
  // HP_PYRAMID_0E_1V
  int refpyramid_0e_1v_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refpyramid_0e_1v_newelstypes[] =
    {
      HP_TET_0E_1V,
      HP_TET,
      HP_NONE,
    };
  int refpyramid_0e_1v_newels[][8] =
    {
      { 1, 2, 4, 5 },
      { 2, 3, 4, 5 },
    };
  HPRef_Struct refpyramid_0e_1v =
    {
      HP_PYRAMID,
      refpyramid_0e_1v_splitedges, 
      0, 0,
      refpyramid_0e_1v_newelstypes, 
      refpyramid_0e_1v_newels
    };


// singular edges 1-2 1-4 singular point 1 
  // HP_PYRAMID_EDGES
  int refpyramid_edges_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refpyramid_edges_newelstypes[] =
    {
      HP_TET_1E_1VA,
      HP_TET_1E_1VA,
      HP_NONE,
    };
  int refpyramid_edges_newels[][8] =
    {
      { 1, 2, 3, 5 },
      { 1, 4, 5, 3 },
    };
  HPRef_Struct refpyramid_edges =
    {
      HP_PYRAMID,
      refpyramid_edges_splitedges, 
      0, 0,
      refpyramid_edges_newelstypes, 
      refpyramid_edges_newels
    };



// singular face 1-2-5 singular point 5
  // HP_PYRAMID_1FB_0E_1VA
  int refpyramid_1fb_0e_1va_splitedges[][3] =
    {
      { 1, 4, 6 },
      { 2, 3, 7 },
      { 5, 1, 8 },
      { 5, 2, 9 },
      { 5, 3, 10 },
      { 5, 4, 11 },
      { 0, 0, 0 },
    };

  HPREF_ELEMENT_TYPE refpyramid_1fb_0e_1va_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PYRAMID_1FB_0E_1VA,
      HP_PRISM,
      HP_NONE,
    };
  int refpyramid_1fb_0e_1va_newels[][8] =
    {
      { 1, 8, 9, 2, 6, 11, 10, 7 },
      { 8, 9, 10, 11, 5 },
      { 3, 7, 10, 4, 6, 11 }
    };
  HPRef_Struct refpyramid_1fb_0e_1va =
    {
      HP_PYRAMID,
      refpyramid_1fb_0e_1va_splitedges, 
      0, 0,
      refpyramid_1fb_0e_1va_newelstypes, 
      refpyramid_1fb_0e_1va_newels
    };




