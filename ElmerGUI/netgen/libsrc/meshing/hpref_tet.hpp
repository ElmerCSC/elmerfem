
 

// HP_TET
int reftet_splitedges[][3] =
{
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_newelstypes[] =
{
  HP_TET,
  HP_NONE,
};
int reftet_newels[][8] =
{
  { 1, 2, 3, 4 },
};
HPRef_Struct reftet =
{
  HP_TET,
  reftet_splitedges, 
  0, 0,
  reftet_newelstypes, 
  reftet_newels
};



/* *********** Tet - Refinement - 0 edges *************** */

// HP_TET_0E_1V
int reftet_0e_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_0e_1v_newelstypes[] =
{
  HP_TET_0E_1V,
  HP_PRISM,
  HP_NONE,
};
int reftet_0e_1v_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 3, 4 }
};
HPRef_Struct reftet_0e_1v =
{
  HP_TET,
  reftet_0e_1v_splitedges, 
  0, 0,
  reftet_0e_1v_newelstypes, 
  reftet_0e_1v_newels
};



// HP_TET_0E_2V
int reftet_0e_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_0e_2v_newelstypes[] =
{
  HP_TET_0E_1V,
  HP_TET_0E_1V,
  HP_PRISM,
  HP_PRISM,
  HP_NONE,
};
int reftet_0e_2v_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 2, 10, 9, 8 },
  { 5, 6, 7, 8, 9, 10 },
  { 4, 10, 7, 3, 9, 6 },
};
HPRef_Struct reftet_0e_2v =
{
  HP_TET,
  reftet_0e_2v_splitedges, 
  0, 0,
  reftet_0e_2v_newelstypes, 
  reftet_0e_2v_newels
};





// HP_TET_0E_3V
int reftet_0e_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_0e_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 2, 3, 1, 15 },
    { 3, 1, 2, 16 },
    { 0, 0, 0, 0 },
  };
HPREF_ELEMENT_TYPE reftet_0e_3v_newelstypes[] =
{
  HP_PYRAMID_0E_1V,
  HP_PYRAMID_0E_1V,
  HP_PYRAMID_0E_1V,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_TET,
  HP_NONE,
};
int reftet_0e_3v_newels[][8] =
{
  { 1, 5, 14, 6, 7 },
  { 2, 9, 15, 8, 10 },
  { 3, 11, 16, 12, 13 },
  { 5, 14, 7, 8, 15, 10 },
  { 9, 15, 10, 12, 16, 13 },
  { 6, 7, 14, 11, 13, 16 },
  { 14, 15, 16, 7, 10, 13 },
  { 7, 10, 13, 4 }
};
HPRef_Struct reftet_0e_3v =
{
  HP_TET,
  reftet_0e_3v_splitedges, 
  reftet_0e_3v_splitfaces, 
  0,
  reftet_0e_3v_newelstypes, 
  reftet_0e_3v_newels
};





// HP_TET_0E_4V
int reftet_0e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_0e_4v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 1, 2, 4, 18 },
    { 1, 3, 4, 19 },

    { 2, 1, 3, 20 },
    { 2, 1, 4, 21 },
    { 2, 3, 4, 22 },

    { 3, 1, 2, 23 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },

    { 4, 1, 2, 26 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 },
  };
int reftet_0e_4v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 29 },
    { 2, 3, 4, 1, 30 },
    { 3, 4, 1, 2, 31 },
    { 4, 1, 2, 3, 32 },
    { 0 },
  };
HPREF_ELEMENT_TYPE reftet_0e_4v_newelstypes[] =
{
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_TET,
  HP_NONE,
};
int reftet_0e_4v_newels[][8] =
{
  { 1, 5, 17, 6, 7, 18, 29, 19 },
  { 2, 9, 20, 8, 10, 22, 30, 21 },
  { 3, 11, 23, 12, 13, 24, 31, 25 },
  { 4, 15, 26, 14, 16, 28, 32, 27 },
  { 5, 17, 18, 8, 20, 21 },
  { 18, 17, 29, 21, 20, 30 },
  { 6, 19, 17,  11, 24, 23 },
  { 17, 19, 29,  23, 24, 31 },
  { 7, 18, 19, 14, 26, 27 },
  { 19, 18, 29, 27, 26, 32 },
  { 9, 20, 22, 12, 23, 25 },
  { 22, 20, 30, 25, 23, 31 },
  { 10, 22, 21, 15, 28, 26 },
  { 21, 22, 30, 26, 28, 32 },
  { 13, 24, 25, 16, 27, 28 },
  { 25, 24, 31, 28, 27, 32 },
  { 17, 20, 23, 29, 30, 31 },
  { 18, 26, 21, 29, 32, 30 },
  { 19, 24, 27, 29, 31, 32 },
  { 22, 28, 25, 30, 32, 31 },
  { 29, 30, 31, 32 },
};
HPRef_Struct reftet_0e_4v =
{
  HP_TET,
  reftet_0e_4v_splitedges, 
  reftet_0e_4v_splitfaces, 
  reftet_0e_4v_splitelements, 
  reftet_0e_4v_newelstypes, 
  reftet_0e_4v_newels
};

















/* *********** Tet - Refinement - 1 edge *************** */



// HP_TET_1E_0V
int reftet_1e_0v_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_0v_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1e_0v_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 7, 3, 5, 8, 4, 6 }
};
HPRef_Struct reftet_1e_0v =
{
  HP_TET,
  reftet_1e_0v_splitedges, 
  0, 0,
  reftet_1e_0v_newelstypes, 
  reftet_1e_0v_newels
};





// HP_TET_1E_1VA
int reftet_1e_1va_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 1, 2, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_1va_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1e_1va_newels[][8] =
{
  { 1, 9, 5, 6 },
  { 9, 5, 6, 2, 7, 8 },
  { 7, 3, 5, 8, 4, 6 }
};
HPRef_Struct reftet_1e_1va =
{
  HP_TET,
  reftet_1e_1va_splitedges, 
  0, 0,
  reftet_1e_1va_newelstypes, 
  reftet_1e_1va_newels
};






// HP_TET_1E_1VB
int reftet_1e_1vb_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 4, 1, 9 },
  { 4, 2, 10 },
  { 4, 3, 11 },
  { 0, 0, 0 }
};
int reftet_1e_1vb_splitelements[][5] =
{
  { 4, 1, 2, 3, 12 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_1vb_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_1vb_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 4, 11, 10, 9 },
  { 7, 8, 10, 11, 12 },
  { 3, 7, 11, 12 },
  { 5, 11, 9, 6, 12 },
  { 5, 3, 11, 12 },
  { 6, 9, 10, 8, 12 },
  { 5, 7, 3, 12 },
  { 5, 6, 8, 7, 12 },
  { 9, 11, 10, 12 }
};
HPRef_Struct reftet_1e_1vb =
{
  HP_TET,
  reftet_1e_1vb_splitedges, 
  0,
  reftet_1e_1vb_splitelements, 
  reftet_1e_1vb_newelstypes, 
  reftet_1e_1vb_newels
};








// HP_TET_1E_2VA
int reftet_1e_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_2va_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1e_2va_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 5, 6, 7, 8, 9, 10 },
  { 4, 10, 7, 3, 9, 6 },
};
HPRef_Struct reftet_1e_2va =
{
  HP_TET,
  reftet_1e_2va_splitedges, 
  0, 0,
  reftet_1e_2va_newelstypes, 
  reftet_1e_2va_newels
};







// HP_TET_1E_2VB
int reftet_1e_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 1, 10 },
  { 3, 2, 11 },
  { 3, 4, 12 },
  { 0, 0, 0 }
};
int reftet_1e_2vb_splitelements[][5] =
{
  { 3, 4, 1, 2, 13 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_2vb_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_2vb_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 8, 9 },
  { 3, 10, 11, 12 },

  { 8, 9, 12, 11, 13 },
  { 4, 12, 9, 13 },
  { 6, 10, 12, 7, 13 },
  { 4, 7, 12, 13 },
  { 6, 8, 11, 10, 13 },
  { 4, 9, 7, 13 },
  { 6, 7, 9, 8, 13 },
  { 10, 11, 12, 13 },
};
HPRef_Struct reftet_1e_2vb =
{
  HP_TET,
  reftet_1e_2vb_splitedges, 
  0,
  reftet_1e_2vb_splitelements, 
  reftet_1e_2vb_newelstypes, 
  reftet_1e_2vb_newels
};






// HP_TET_1E_2VC
int reftet_1e_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 4, 1, 10 },
  { 4, 2, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};
int reftet_1e_2vc_splitelements[][5] =
{
  { 4, 1, 2, 3, 13 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_2vc_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_2vc_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 8, 9 },
  { 4, 11, 10, 12 },
  { 8, 9, 11, 12, 13 },
  { 3, 8, 12, 13 },
  { 7, 6, 12, 10, 13 },
  { 3, 12, 6, 13 },
  { 9, 7, 10, 11, 13 },
  { 3, 6, 8, 13 },
  { 6, 7, 9, 8, 13 },
  { 10, 12, 11, 13 }
};
HPRef_Struct reftet_1e_2vc =
{
  HP_TET,
  reftet_1e_2vc_splitedges, 
  0,
  reftet_1e_2vc_splitelements, 
  reftet_1e_2vc_newelstypes, 
  reftet_1e_2vc_newels
};








/*

// HP_TET_1E_2VD
int reftet_1e_2vd_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 3, 1, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 1, 12 },
  { 4, 2, 13 },
  { 4, 3, 14 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_2vd_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_TET_0E_1V,
  HP_PRISM,
  HP_HEX,
  HP_NONE,
};
int reftet_1e_2vd_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 4, 13, 12, 14 },
  { 3, 10, 11, 9 },
  { 14, 13, 12, 11, 10, 9 },
  { 6, 12, 13, 8, 5, 9, 10, 7 },
};
HPRef_Struct reftet_1e_2vd =
{
  HP_TET,
  reftet_1e_2vd_splitedges, 
  0, 0,
  reftet_1e_2vd_newelstypes, 
  reftet_1e_2vd_newels
};

*/




//  HP_TET_1E_2VD,  // 1 v on edge
int reftet_1e_2vd_splitedges[][3] =
{
  // { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  // { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_2vd_splitfaces[][4] =
  {
    { 1, 3, 4, 19 },
    { 2, 3, 4, 22 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_1e_2vd_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_TET_0E_1V,
    HP_TET_0E_1V,
    HP_PRISM,
    HP_HEX,
    HP_PYRAMID,
    HP_HEX,
    HP_PYRAMID,
    HP_PRISM,
    HP_PRISM,
    HP_NONE,
  };
int reftet_1e_2vd_newels[][8] =
{
  { 1, 6, 7, 2, 9, 10 },
  { 3, 11, 12, 13 },
  { 4, 16, 15, 14 },
  { 7, 6, 19, 10, 9, 22 },
  { 7, 19, 27, 14, 10, 22, 28, 15 },
  { 14, 15, 28, 27, 16 },
  { 9, 6, 19, 22, 12, 11, 24, 25 },
  { 12, 11, 24, 25, 13 },
  { 19, 24, 27, 22, 25, 28 },
  { 16, 28, 27, 13, 25, 24 }
};
HPRef_Struct reftet_1e_2vd =
{
  HP_TET,
  reftet_1e_2vd_splitedges, 
  reftet_1e_2vd_splitfaces, 
  0,
  reftet_1e_2vd_newelstypes, 
  reftet_1e_2vd_newels
};















// HP_TET_1E_3VA
int reftet_1e_3va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_1e_3va_splitelements[][5] =
{
  { 1, 2, 3, 4, 14 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_3va_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  HP_TET_0E_1V,

  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_3va_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 11, 12, 13 },

  { 6, 7, 10, 9, 14 },
  { 4, 10, 7, 14 },
  { 9, 10, 13, 12, 14 },
  { 4, 13, 10, 14 },
  { 6, 11, 13, 7, 14 },
  { 4, 7, 13, 14 },
  { 6, 11, 12, 9, 14 },
  { 11, 13, 12, 14 },
};

HPRef_Struct reftet_1e_3va =
{
  HP_TET,
  reftet_1e_3va_splitedges, 
  0,
  reftet_1e_3va_splitelements, 
  reftet_1e_3va_newelstypes, 
  reftet_1e_3va_newels
};






















//  HP_TET_1E_3VB,  // 1 v on edge
int reftet_1e_3vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  // { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_3vb_splitfaces[][4] =
  {
    { 1, 3, 4, 19 },
    { 2, 3, 4, 22 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_1e_3vb_newelstypes[] =
  {
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_TET_0E_1V,
    HP_TET_0E_1V,
    HP_PRISM,
    HP_HEX,
    HP_PYRAMID,
    HP_HEX,
    HP_PYRAMID,
    HP_PRISM,
    HP_PRISM,
    HP_NONE,
  };
int reftet_1e_3vb_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 9, 10 },
  { 3, 11, 12, 13 },
  { 4, 16, 15, 14 },
  { 7, 6, 19, 10, 9, 22 },
  { 7, 19, 27, 14, 10, 22, 28, 15 },
  { 14, 15, 28, 27, 16 },
  { 9, 6, 19, 22, 12, 11, 24, 25 },
  { 12, 11, 24, 25, 13 },
  { 19, 24, 27, 22, 25, 28 },
  { 16, 28, 27, 13, 25, 24 }
};
HPRef_Struct reftet_1e_3vb =
{
  HP_TET,
  reftet_1e_3vb_splitedges, 
  reftet_1e_3vb_splitfaces, 
  0,
  reftet_1e_3vb_newelstypes, 
  reftet_1e_3vb_newels
};






/*
// HP_TET_1E_4V
int reftet_1e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_4v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 1, 2, 4, 18 },
    { 1, 3, 4, 19 },

    { 2, 1, 3, 20 },
    { 2, 1, 4, 21 },
    { 2, 3, 4, 22 },

    { 3, 1, 2, 23 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },

    { 4, 1, 2, 26 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 },
  };
int reftet_1e_4v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 29 },
    { 2, 3, 4, 1, 30 },
    { 3, 4, 1, 2, 31 },
    { 4, 1, 2, 3, 32 },
    { 0 },
  };
HPREF_ELEMENT_TYPE reftet_1e_4v_newelstypes[] =
{
  HP_HEX_1E_1V,
  HP_HEX_1E_1V,
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_PRISM_SINGEDGE, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_TET,
  HP_NONE,
};
int reftet_1e_4v_newels[][8] =
{
  { 1, 5, 17, 6, 7, 18, 29, 19 },
  //  { 2, 9, 20, 8, 10, 22, 30, 21 },
  { 2, 8, 21, 10, 9, 20, 30, 22 },
  { 3, 11, 23, 12, 13, 24, 31, 25 },
  { 4, 15, 26, 14, 16, 28, 32, 27 },
  { 5, 17, 18, 8, 20, 21 },
  { 18, 17, 29, 21, 20, 30 },
  { 6, 19, 17,  11, 24, 23 },
  { 17, 19, 29,  23, 24, 31 },
  { 7, 18, 19, 14, 26, 27 },
  { 19, 18, 29, 27, 26, 32 },
  { 9, 20, 22, 12, 23, 25 },
  { 22, 20, 30, 25, 23, 31 },
  { 10, 22, 21, 15, 28, 26 },
  { 21, 22, 30, 26, 28, 32 },
  { 13, 24, 25, 16, 27, 28 },
  { 25, 24, 31, 28, 27, 32 },
  { 17, 20, 23, 29, 30, 31 },
  { 18, 26, 21, 29, 32, 30 },
  { 19, 24, 27, 29, 31, 32 },
  { 22, 28, 25, 30, 32, 31 },

  { 29, 30, 31, 32 },
};
HPRef_Struct reftet_1e_4v =
{
  HP_TET,
  reftet_1e_4v_splitedges, 
  reftet_1e_4v_splitfaces, 
  reftet_1e_4v_splitelements, 
  reftet_1e_4v_newelstypes, 
  reftet_1e_4v_newels
};
*/




// HP_TET_1E_4V
int reftet_1e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_4v_splitfaces[][4] =
  {
    { 1, 3, 4, 17 },
    { 2, 3, 4, 18 },

    { 3, 1, 4, 19 },
    { 3, 2, 4, 20 },

    { 4, 1, 3, 21 },
    { 4, 2, 3, 22 },
    { 0, 0, 0, 0 },
  };

HPREF_ELEMENT_TYPE reftet_1e_4v_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  //  HP_TET_1E_1VA,
  //  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_HEX, 
  HP_HEX, 
  HP_PRISM,
  HP_PRISM,

  HP_PYRAMID,
  HP_TET_0E_1V,

  HP_PYRAMID,
  HP_TET_0E_1V,

  HP_NONE,
};

int reftet_1e_4v_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },

  { 5, 6, 7, 8, 9, 10 },
  { 7, 6, 17, 10, 9, 18 },

  { 7, 10, 18, 17, 14, 15, 22, 21 },
  { 9, 6, 17, 18, 12, 11, 19, 20 },

  { 17, 19, 21, 18, 20, 22 },
  { 16, 22, 21, 13, 20, 19 },

  { 14, 15, 22, 21, 16 },
  { 4, 14, 16, 15 },
  { 12, 11, 19, 20, 13 },
  { 3, 11, 12, 13 },



  { 1, 5, 17, 6, 7, 18, 29, 19 },
  //  { 2, 9, 20, 8, 10, 22, 30, 21 },
  { 2, 8, 21, 10, 9, 20, 30, 22 },
  { 3, 11, 23, 12, 13, 24, 31, 25 },
  { 4, 15, 26, 14, 16, 28, 32, 27 },
  { 5, 17, 18, 8, 20, 21 },
  { 18, 17, 29, 21, 20, 30 },
  { 6, 19, 17,  11, 24, 23 },
  { 17, 19, 29,  23, 24, 31 },
  { 7, 18, 19, 14, 26, 27 },
  { 19, 18, 29, 27, 26, 32 },
  { 9, 20, 22, 12, 23, 25 },
  { 22, 20, 30, 25, 23, 31 },
  { 10, 22, 21, 15, 28, 26 },
  { 21, 22, 30, 26, 28, 32 },
  { 13, 24, 25, 16, 27, 28 },
  { 25, 24, 31, 28, 27, 32 },
  { 17, 20, 23, 29, 30, 31 },
  { 18, 26, 21, 29, 32, 30 },
  { 19, 24, 27, 29, 31, 32 },
  { 22, 28, 25, 30, 32, 31 },

  { 29, 30, 31, 32 },
};
HPRef_Struct reftet_1e_4v =
{
  HP_TET,
  reftet_1e_4v_splitedges, 
  reftet_1e_4v_splitfaces, 
  0, 
  reftet_1e_4v_newelstypes, 
  reftet_1e_4v_newels
};













//  HP_TET_2EA_0V,  // 2 edges connected
int reftet_2ea_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_0v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_0v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 3, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_0v =
{
  HP_TET,
  reftet_2ea_0v_splitedges, 
  reftet_2ea_0v_splitfaces, 
  0,
  reftet_2ea_0v_newelstypes, 
  reftet_2ea_0v_newels
};






//  HP_TET_2EA_1VA,  // 2 edges connected
int reftet_2ea_1va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_1va_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_1va_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_1va_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 5, 17, 7, 8, 9, 10 },
  { 2, 8, 10, 9 },
  { 6, 7, 17, 3, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_1va =
{
  HP_TET,
  reftet_2ea_1va_splitedges, 
  reftet_2ea_1va_splitfaces, 
  0,
  reftet_2ea_1va_newelstypes, 
  reftet_2ea_1va_newels
};








//  HP_TET_2EA_1VB, 
int reftet_2ea_1vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_1vb_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_1vb_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_1vb_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 3, 11, 12, 13 },
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_1vb =
{
  HP_TET,
  reftet_2ea_1vb_splitedges, 
  reftet_2ea_1vb_splitfaces, 
  0,
  reftet_2ea_1vb_newelstypes, 
  reftet_2ea_1vb_newels
};






//  HP_TET_2EA_1VC,  // 2 edges connected
int reftet_2ea_1vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  //  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_1vc_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_1vc_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_1vc_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    //    HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_1vc_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  // { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 3, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_1vc =
{
  HP_TET,
  reftet_2ea_1vc_splitedges, 
  reftet_2ea_1vc_splitfaces, 
  reftet_2ea_1vc_splitelements, 
  reftet_2ea_1vc_newelstypes, 
  reftet_2ea_1vc_newels
};











 
//  HP_TET_2EA_2VA, 
int reftet_2ea_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_2va_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_2va_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_2va_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 3, 11, 12, 13 },
  { 2, 8, 10, 9 },
  { 5, 17, 7, 8, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_2va =
{
  HP_TET,
  reftet_2ea_2va_splitedges, 
  reftet_2ea_2va_splitfaces, 
  0,
  reftet_2ea_2va_newelstypes, 
  reftet_2ea_2va_newels
};











//  HP_TET_2EA_2VB,  // 2 edges connected
int reftet_2ea_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  //  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_2vb_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_2vb_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_2vb_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //  HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_2vb_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 10, 9 },
  //  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 8, 9, 10 },
  { 6, 7, 17, 3, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_2vb =
{
  HP_TET,
  reftet_2ea_2vb_splitedges, 
  reftet_2ea_2vb_splitfaces, 
  reftet_2ea_2vb_splitelements, 
  reftet_2ea_2vb_newelstypes, 
  reftet_2ea_2vb_newels
};







 


//  HP_TET_2EA_2VC,  // 2 edges connected
int reftet_2ea_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  //  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_2vc_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_2vc_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_2vc_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_2vc_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  //  { 2, 8, 10, 9 },
  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_2vc =
{
  HP_TET,
  reftet_2ea_2vc_splitedges, 
  reftet_2ea_2vc_splitfaces, 
  reftet_2ea_2vc_splitelements, 
  reftet_2ea_2vc_newelstypes, 
  reftet_2ea_2vc_newels
};








//  HP_TET_2EA_3V,  // 2 edges connected
int reftet_2ea_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_3v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_3v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_3v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 8, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_3v =
{
  HP_TET,
  reftet_2ea_3v_splitedges, 
  reftet_2ea_3v_splitfaces, 
  reftet_2ea_3v_splitelements, 
  reftet_2ea_3v_newelstypes, 
  reftet_2ea_3v_newels
};







//  HP_TET_2EB_0V,  // 2 opposite edges
int reftet_2eb_0v_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 3, 1, 9 },
  { 3, 2, 10 },
  { 4, 1, 11 },
  { 4, 2, 12 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_0v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_0v_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 3, 9, 10, 4, 11, 12 },
  { 6, 11, 12, 8, 5, 9, 10, 7 },
};
HPRef_Struct reftet_2eb_0v =
{
  HP_TET,
  reftet_2eb_0v_splitedges, 
  0, 0,
  reftet_2eb_0v_newelstypes, 
  reftet_2eb_0v_newels
};


//  HP_TET_2EB_1V,    // V1


//  HP_TET_2EB_1V,  // 2 opposite edges, V1
int reftet_2eb_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_1v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_1v_newels[][8] =
{
  { 5, 6, 7, 2, 9, 10 },
  { 4, 15, 14, 3, 12, 11 },
  { 1, 5, 6, 7 },
  //  { 2, 8, 10, 9 },
  //  { 3, 13, 11, 12 },
  //  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_1v =
{
  HP_TET,
  reftet_2eb_1v_splitedges, 
  0, 0,
  reftet_2eb_1v_newelstypes, 
  reftet_2eb_1v_newels
};



//  HP_TET_2EB_2VA,  // 2 opposite edges, V1,2
int reftet_2eb_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_2va_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_2va_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 4, 15, 14, 3, 12, 11 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  //  { 3, 13, 11, 12 },
  //  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_2va =
{
  HP_TET,
  reftet_2eb_2va_splitedges, 
  0, 0,
  reftet_2eb_2va_newelstypes, 
  reftet_2eb_2va_newels
};


//  HP_TET_2EB_2VB,   // V1,3
int reftet_2eb_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_2vb_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_2vb_newels[][8] =
{
  { 5, 6, 7, 2, 9, 10 },
  { 4, 15, 14, 13, 12, 11 },
  { 1, 5, 6, 7 },
  // { 2, 8, 10, 9 },
  { 3, 13, 11, 12 },
  // { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_2vb =
{
  HP_TET,
  reftet_2eb_2vb_splitedges, 
  0, 0,
  reftet_2eb_2vb_newelstypes, 
  reftet_2eb_2vb_newels
};




//  HP_TET_2EB_2VC,   // V1,4
int reftet_2eb_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_2vc_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_2vc_newels[][8] =
{
  { 5, 6, 7, 2, 9, 10 },
  { 16, 15, 14, 3, 12, 11 },
  { 1, 5, 6, 7 },
  // { 2, 8, 10, 9 },
  // { 3, 13, 11, 12 },
  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_2vc =
{
  HP_TET,
  reftet_2eb_2vc_splitedges, 
  0, 0,
  reftet_2eb_2vc_newelstypes, 
  reftet_2eb_2vc_newels
};






//  HP_TET_2EB_3V,    // V1,2,3
int reftet_2eb_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_3v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_3v_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 4, 15, 14, 13, 12, 11 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 13, 11, 12 },
  // { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_3v =
{
  HP_TET,
  reftet_2eb_3v_splitedges, 
  0, 0,
  reftet_2eb_3v_newelstypes, 
  reftet_2eb_3v_newels
};






//  HP_TET_2EB_4V,  // 2 opposite edges
int reftet_2eb_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_4v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_4v_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 16, 15, 14, 13, 12, 11 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 13, 11, 12 },
  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_4v =
{
  HP_TET,
  reftet_2eb_4v_splitedges, 
  0, 0,
  reftet_2eb_4v_newelstypes, 
  reftet_2eb_4v_newels
};

















//  HP_TET_3EA_0V,  
int reftet_3ea_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 0, 0, 0 }
};
int reftet_3ea_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_0v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_0v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_HEX_1E_0V,
    HP_HEX_1E_0V,
    HP_HEX_1E_0V,
    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_0v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },
  { 5, 2, 8, 14, 15, 9, 17, 20 },
  { 3, 6, 14, 10, 11, 16, 20, 18 },
  { 7, 4, 12, 15, 16, 13, 19, 20 },
  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_0v =
{
  HP_TET,
  reftet_3ea_0v_splitedges, 
  reftet_3ea_0v_splitfaces, 
  reftet_3ea_0v_splitelements, 
  reftet_3ea_0v_newelstypes, 
  reftet_3ea_0v_newels
};










//  HP_TET_3EA_1V,  
int reftet_3ea_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 2, 1, 21 },
  { 3, 1, 22 },
  { 4, 1, 23 },
  { 0, 0, 0 }
};
int reftet_3ea_1v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_1v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_1v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,

    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_1v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },

  { 2, 21, 9, 8 },
  { 5, 14, 15, 21, 8, 9 },
  { 15, 14, 20, 9, 8, 17 },
  //  { 3, 22, 10, 11 },
  //  { 6, 16, 14, 22, 11, 10 },
  { 6, 16, 14, 3, 11, 10 },
  { 14, 16, 20, 10, 11, 18 },
  //  { 4, 23, 13, 12 },
  //  { 7, 15, 16, 23, 12, 13 },
  { 7, 15, 16, 4, 12, 13 },
  { 16, 15, 20, 13, 12, 19 },

  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_1v =
{
  HP_TET,
  reftet_3ea_1v_splitedges, 
  reftet_3ea_1v_splitfaces, 
  reftet_3ea_1v_splitelements, 
  reftet_3ea_1v_newelstypes, 
  reftet_3ea_1v_newels
};










//  HP_TET_3EA_2V,  
int reftet_3ea_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 2, 1, 21 },
  { 3, 1, 22 },
  { 4, 1, 23 },
  { 0, 0, 0 }
};
int reftet_3ea_2v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_2v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_2v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,

    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_2v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },

  { 2, 21, 9, 8 },
  { 5, 14, 15, 21, 8, 9 },
  { 15, 14, 20, 9, 8, 17 },
  { 3, 22, 10, 11 },
  { 6, 16, 14, 22, 11, 10 },
  { 14, 16, 20, 10, 11, 18 },
  //  { 4, 23, 13, 12 },
  { 7, 15, 16, 4, 12, 13 },
  { 16, 15, 20, 13, 12, 19 },

  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_2v =
{
  HP_TET,
  reftet_3ea_2v_splitedges, 
  reftet_3ea_2v_splitfaces, 
  reftet_3ea_2v_splitelements, 
  reftet_3ea_2v_newelstypes, 
  reftet_3ea_2v_newels
};








//  HP_TET_3EA_3V,  
int reftet_3ea_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 2, 1, 21 },
  { 3, 1, 22 },
  { 4, 1, 23 },
  { 0, 0, 0 }
};
int reftet_3ea_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_3v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_3v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,

    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_3v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },

  { 2, 21, 9, 8 },
  { 5, 14, 15, 21, 8, 9 },
  { 15, 14, 20, 9, 8, 17 },
  { 3, 22, 10, 11 },
  { 6, 16, 14, 22, 11, 10 },
  { 14, 16, 20, 10, 11, 18 },
  { 4, 23, 13, 12 },
  { 7, 15, 16, 23, 12, 13 },
  { 16, 15, 20, 13, 12, 19 },

  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_3v =
{
  HP_TET,
  reftet_3ea_3v_splitedges, 
  reftet_3ea_3v_splitfaces, 
  reftet_3ea_3v_splitelements, 
  reftet_3ea_3v_newelstypes, 
  reftet_3ea_3v_newels
};







//  HP_TET_3EV_0V,  
int reftet_3eb_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  //  { 3, 2, 12 },
  { 3, 4, 13 },
  //  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3eb_0v_splitfaces[][4] =
  {
    { 1, 2, 4, 17 },
    { 2, 1, 3, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3eb_0v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3eb_0v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3eb_0v_newels[][8] =
{
  { 1, 7, 17, 5, 6 },
  { 2, 9, 18, 8, 10 },
  //  { 3, 12, 13, 11 },
  //  { 4, 14, 16, 15 },
  { 5, 6, 17, 8, 18, 10 },
  { 7, 17, 6, 4, 15, 16 },
  { 9, 18, 10, 3, 11, 13 },
  
  { 10, 15, 16, 13, 20 },
  { 6, 11, 13, 16, 20 },
  { 10, 17, 15, 20 },
  { 6, 18, 11, 20 },
  { 6, 17, 10, 18, 20 },
  { 6, 16, 15, 17, 20 },
  { 18, 10, 13, 11, 20 },
};
HPRef_Struct reftet_3eb_0v =
{
  HP_TET,
  reftet_3eb_0v_splitedges, 
  reftet_3eb_0v_splitfaces, 
  reftet_3eb_0v_splitelements, 
  reftet_3eb_0v_newelstypes, 
  reftet_3eb_0v_newels
};









//  HP_TET_3EV_1V,  
int reftet_3eb_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  //  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3eb_1v_splitfaces[][4] =
  {
    { 1, 2, 4, 17 },
    { 2, 1, 3, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3eb_1v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3eb_1v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3eb_1v_newels[][8] =
{
  { 1, 7, 17, 5, 6 },
  { 2, 9, 18, 8, 10 },
  { 3, 12, 13, 11 },
  //  { 4, 14, 16, 15 },
  { 5, 6, 17, 8, 18, 10 },
  { 7, 17, 6, 4, 15, 16 },
  { 9, 18, 10, 12, 11, 13 },
  
  { 10, 15, 16, 13, 20 },
  { 6, 11, 13, 16, 20 },
  { 10, 17, 15, 20 },
  { 6, 18, 11, 20 },
  { 6, 17, 10, 18, 20 },
  { 6, 16, 15, 17, 20 },
  { 18, 10, 13, 11, 20 },
};
HPRef_Struct reftet_3eb_1v =
{
  HP_TET,
  reftet_3eb_1v_splitedges, 
  reftet_3eb_1v_splitfaces, 
  reftet_3eb_1v_splitelements, 
  reftet_3eb_1v_newelstypes, 
  reftet_3eb_1v_newels
};








//  HP_TET_3EV_2V,  
int reftet_3eb_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3eb_2v_splitfaces[][4] =
  {
    { 1, 2, 4, 17 },
    { 2, 1, 3, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3eb_2v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3eb_2v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3eb_2v_newels[][8] =
{
  { 1, 7, 17, 5, 6 },
  { 2, 9, 18, 8, 10 },
  { 3, 12, 13, 11 },
  { 4, 14, 16, 15 },
  { 5, 6, 17, 8, 18, 10 },
  { 7, 17, 6, 14, 15, 16 },
  { 9, 18, 10, 12, 11, 13 },
  
  { 10, 15, 16, 13, 20 },
  { 6, 11, 13, 16, 20 },
  { 10, 17, 15, 20 },
  { 6, 18, 11, 20 },
  { 6, 17, 10, 18, 20 },
  { 6, 16, 15, 17, 20 },
  { 18, 10, 13, 11, 20 },
};
HPRef_Struct reftet_3eb_2v =
{
  HP_TET,
  reftet_3eb_2v_splitedges, 
  reftet_3eb_2v_splitfaces, 
  reftet_3eb_2v_splitelements, 
  reftet_3eb_2v_newelstypes, 
  reftet_3eb_2v_newels
};













//  HP_TET_3EC_0V,  
int reftet_3ec_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  //  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  //  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3ec_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 1, 4, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3ec_0v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ec_0v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3ec_0v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 18, 10, 9 },
  //  { 3, 11, 12, 13 },
  //  { 4, 15, 14, 16 },
  { 5, 17, 7, 8, 9, 18 },
  { 6, 7, 17, 3, 13, 12 },
  { 10, 9, 18, 4, 16, 14 },
  
  { 9, 16, 13, 12, 20 },
  { 7, 13, 16, 14, 20 },
  { 7, 14, 18, 20 },
  { 9, 12, 17, 20 },
  { 17, 7, 18, 9, 20 },
  { 7, 17, 12, 13, 20 },
  { 9, 18, 14, 16, 20 },
};
HPRef_Struct reftet_3ec_0v =
{
  HP_TET,
  reftet_3ec_0v_splitedges, 
  reftet_3ec_0v_splitfaces, 
  reftet_3ec_0v_splitelements, 
  reftet_3ec_0v_newelstypes, 
  reftet_3ec_0v_newels
};






 


//  HP_TET_3EC_1V,  
int reftet_3ec_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  // { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3ec_1v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 1, 4, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3ec_1v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ec_1v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3ec_1v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 18, 10, 9 },
  { 3, 11, 12, 13 },
  //  { 4, 15, 14, 16 },
  { 5, 17, 7, 8, 9, 18 },
  { 6, 7, 17, 11, 13, 12 },
  { 10, 9, 18, 4, 16, 14 },
  
  { 9, 16, 13, 12, 20 },
  { 7, 13, 16, 14, 20 },
  { 7, 14, 18, 20 },
  { 9, 12, 17, 20 },
  { 17, 7, 18, 9, 20 },
  { 7, 17, 12, 13, 20 },
  { 9, 18, 14, 16, 20 },
};
HPRef_Struct reftet_3ec_1v =
{
  HP_TET,
  reftet_3ec_1v_splitedges, 
  reftet_3ec_1v_splitfaces, 
  reftet_3ec_1v_splitelements, 
  reftet_3ec_1v_newelstypes, 
  reftet_3ec_1v_newels
};








//  HP_TET_3EC_2V,  
int reftet_3ec_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3ec_2v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 1, 4, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3ec_2v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ec_2v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3ec_2v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 18, 10, 9 },
  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 },
  { 5, 17, 7, 8, 9, 18 },
  { 6, 7, 17, 11, 13, 12 },
  { 10, 9, 18, 15, 16, 14 },
  
  { 9, 16, 13, 12, 20 },
  { 7, 13, 16, 14, 20 },
  { 7, 14, 18, 20 },
  { 9, 12, 17, 20 },
  { 17, 7, 18, 9, 20 },
  { 7, 17, 12, 13, 20 },
  { 9, 18, 14, 16, 20 },
};
HPRef_Struct reftet_3ec_2v =
{
  HP_TET,
  reftet_3ec_2v_splitedges, 
  reftet_3ec_2v_splitfaces, 
  reftet_3ec_2v_splitelements, 
  reftet_3ec_2v_newelstypes, 
  reftet_3ec_2v_newels
};










/* ************************ 1 singular face ******************** */


// HP_TET_1F_0E_0V
int reftet_1f_0e_0v_splitedges[][3] =
{
  { 2, 1, 5 },
  { 3, 1, 6 },
  { 4, 1, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_0v_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_TET,
  HP_NONE,
};
int reftet_1f_0e_0v_newels[][8] =
{
  { 3, 2, 4, 6, 5, 7 },
  { 5, 7, 6, 1 }
};
HPRef_Struct reftet_1f_0e_0v =
{
  HP_TET,
  reftet_1f_0e_0v_splitedges, 
  0, 0,
  reftet_1f_0e_0v_newelstypes, 
  reftet_1f_0e_0v_newels
};





// HP_TET_1F_0E_1VA    ... singular vertex in face
int reftet_1f_0e_1va_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 3, 6 },
  { 2, 4, 7 },
  { 3, 1, 8 },
  { 4, 1, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_1va_newelstypes[] =
{
  HP_HEX_1F_0E_0V,
  HP_TET_1F_0E_1VA,
  HP_TET,
  HP_NONE,
};
int reftet_1f_0e_1va_newels[][8] =
{
  { 3, 6, 7, 4, 8, 5, 5, 9 },
  { 5, 2, 6, 7 },
  { 5, 9, 8, 1 },
};
HPRef_Struct reftet_1f_0e_1va =
{
  HP_TET,
  reftet_1f_0e_1va_splitedges, 
  0, 0,
  reftet_1f_0e_1va_newelstypes, 
  reftet_1f_0e_1va_newels
};





// HP_TET_1F_0E_1VB    ... singular vertex not in face
int reftet_1f_0e_1vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 3, 1, 9 },
  { 4, 1, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_1vb_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_PRISM,
  HP_TET_0E_1V,
  HP_NONE,
};
int reftet_1f_0e_1vb_newels[][8] =
{
  { 2, 4, 3, 8, 10, 9 },
  { 8, 10, 9, 5, 7, 6 }, 
  { 1, 5, 6, 7 },
};
HPRef_Struct reftet_1f_0e_1vb =
{
  HP_TET,
  reftet_1f_0e_1vb_splitedges, 
  0, 0,
  reftet_1f_0e_1vb_newelstypes, 
  reftet_1f_0e_1vb_newels
};








// HP_TET_1F_1EA_0V  ... sing edge is 1..2
int reftet_1f_1ea_0v_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 1, 10 },
  { 4, 1, 11 },
  { 0, 0, 0 }
};

int reftet_1f_1ea_0v_splitfaces[][4] =
  {
    { 2, 1, 3, 12 },
    { 2, 1, 4, 13 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_1f_1ea_0v_newelstypes[] =
{
  HP_HEX_1F_0E_0V,
  //  HP_PRISM,
  HP_PYRAMID_1FB_0E_1VA,
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1f_1ea_0v_newels[][8] =
{
  { 3, 8, 9, 4, 10, 12, 13, 11 },
  // { 2, 9, 8, 7, 13, 12 },
  { 8, 9, 13, 12, 2 },
  { 2, 7, 13, 12 },
  { 7, 13, 12, 1, 6, 5 },
  { 6, 11, 13, 5, 10, 12 }
};
HPRef_Struct reftet_1f_1ea_0v =
{
  HP_TET,
  reftet_1f_1ea_0v_splitedges, 
  reftet_1f_1ea_0v_splitfaces, 
  0, 
  reftet_1f_1ea_0v_newelstypes, 
  reftet_1f_1ea_0v_newels
};








// HP_TET_1F_1EB_0V     singular edge in face, edge is 2-3
int reftet_1f_1eb_0v_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 4, 6 },
  { 3, 1, 7 },
  { 3, 4, 8 },
  { 4, 1, 9 },
  { 0, 0, 0 }
};


HPREF_ELEMENT_TYPE reftet_1f_1eb_0v_newelstypes[] =
{
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FA_0E_0V,
  HP_TET,
  HP_NONE,
};
int reftet_1f_1eb_0v_newels[][8] =
{
  // { 2, 5, 6, 3, 7, 8 },
  { 3, 8, 7, 2, 6, 5 },
  { 6, 4, 8, 5, 9, 7 },
  { 5, 9, 7, 1}
};
HPRef_Struct reftet_1f_1eb_0v =
{
  HP_TET,
  reftet_1f_1eb_0v_splitedges, 
  0, 0, 
  reftet_1f_1eb_0v_newelstypes, 
  reftet_1f_1eb_0v_newels
};










/* ************************ 2 singular faces ******************** */


// HP_TET_2F_0E_0V
int reftet_2f_0e_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 2, 1, 6 },
  { 3, 1, 7 },
  { 3, 2, 8 },
  { 4, 1, 9 },
  { 4, 2, 10 },
  { 0, 0, 0 }
};

int reftet_2f_0e_0v_splitfaces[][4] =
  {
    { 3, 1, 2, 11 },
    { 4, 1, 2, 12 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_2f_0e_0v_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_TET,
  HP_NONE,
};
int reftet_2f_0e_0v_newels[][8] =
{
  { 2, 10, 8, 6, 12, 11 },
  { 1, 7, 9, 5, 11, 12 },
  //   { 3, 11, 8, 4, 12, 10 },
  { 4, 10, 12, 3, 8, 11 }, 
  { 3, 7, 11, 4, 9, 12 },
  { 5, 6, 11, 12 }
};
HPRef_Struct reftet_2f_0e_0v =
{
  HP_TET,
  reftet_2f_0e_0v_splitedges, 
  reftet_2f_0e_0v_splitfaces, 
  0, 
  reftet_2f_0e_0v_newelstypes, 
  reftet_2f_0e_0v_newels
};

