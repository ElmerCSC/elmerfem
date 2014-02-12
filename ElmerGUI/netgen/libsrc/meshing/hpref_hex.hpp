// SZ 

// HP_HEX  ... no refinement
int refhex_splitedges[][3] =
  {
      { 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE refhex_newelstypes[] =
  {
    HP_HEX,
    HP_NONE,
  };
int refhex_newels[][8] =
  {
    { 1, 2, 3, 4, 5, 6, 7, 8 }
  };
HPRef_Struct refhex =
  {
    HP_HEX,
    refhex_splitedges, 
    0, 0,
    refhex_newelstypes, 
    refhex_newels
  };

// HP_HEX_1F  ... face (1 - 4 - 3 -2) singular 
int refhex_1f_0e_0v_splitedges[][3] =
  {
    { 1, 5, 9 },
    { 2, 6, 10 },
    { 3, 7, 11 },
    { 4, 8, 12 }, 
    { 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE refhex_1f_0e_0v_newelstypes[] =
  {
    HP_HEX,
    HP_HEX_1F_0E_0V,
    HP_NONE,
  };
int  refhex_1f_0e_0v_newels[][8] =
  {
    { 9, 10, 11, 12, 5, 6, 7, 8 }, 
    { 1, 2, 3, 4, 9, 10, 11, 12}  
 }; 
HPRef_Struct refhex_1f_0e_0v =
  {
    HP_HEX,
    refhex_1f_0e_0v_splitedges, 
    0, 0,
    refhex_1f_0e_0v_newelstypes, 
    refhex_1f_0e_0v_newels
  };



// HP_HEX_1FA_1FB  ... face (1 - 4 - 3 -2) and face (1-2-6-5) singular 
int refhex_1fa_1fb_0e_0v_splitedges[][3] =
  {
    { 1, 5, 9 },
    { 2, 6, 10 },
    { 3, 7, 11 },
    { 4, 8, 12 },
    { 1, 4, 13 },
    { 2, 3, 14 },  
    { 6, 7, 15 }, 
    { 5, 8, 16 }, 
    { 0, 0, 0 }
  };

int refhex_1fa_1fb_0e_0v_splitfaces[][4] =
  {
    { 2, 3, 6, 17 },
    { 1, 4, 5, 18 },
    { 0, 0, 0, 0 },
  };
HPREF_ELEMENT_TYPE refhex_1fa_1fb_0e_0v_newelstypes[] =
  {
    HP_HEX,
    HP_HEX_1F_0E_0V,
    HP_HEX_1F_0E_0V, 
    HP_HEX_1FA_1FB_0E_0V, 
    HP_NONE,
  };
int  refhex_1fa_1fb_0e_0v_newels[][8] =
  {
    {18, 17, 11, 12, 16, 15, 7, 8}, 
    {13, 14, 3, 4, 18, 17, 11, 12},
    { 5, 6, 10, 9, 16, 15, 17, 18}, 
    { 1, 2, 14, 13, 9, 10, 17, 18} 
  }; 
HPRef_Struct refhex_1fa_1fb_0e_0v =
  {
    HP_HEX,
    refhex_1fa_1fb_0e_0v_splitedges, 
    refhex_1fa_1fb_0e_0v_splitfaces, 0,
    refhex_1fa_1fb_0e_0v_newelstypes, 
    refhex_1fa_1fb_0e_0v_newels
  };



// Refine Dummies 
  // HP_HEX_0E_1V
  int refhex_0e_1v_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refhex_0e_1v_newelstypes[] =
    {
      HP_TET_0E_1V,
      HP_TET,
      HP_TET,
      HP_TET,
      HP_TET,
      HP_TET,
      HP_NONE,
    };
  int refhex_0e_1v_newels[][8] =
    {
      { 1, 5, 2, 4 },
      { 7, 3, 6, 8 },
      { 2, 8, 5, 6 },
      { 2, 8, 6, 3 },
      { 2, 8, 3, 4 },
      { 2, 8, 4, 5 },
    };
  HPRef_Struct refhex_0e_1v =
    {
      HP_HEX,
      refhex_0e_1v_splitedges, 
      0, 0,
      refhex_0e_1v_newelstypes, 
      refhex_0e_1v_newels
    };



// Refine Dummies 
  // HP_HEX_1E_1V
  int refhex_1e_1v_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refhex_1e_1v_newelstypes[] =
    {
      HP_TET_1E_1VA,
      HP_TET,
      HP_TET_0E_1V,
      HP_TET_0E_1V,
      HP_TET_0E_1V,
      HP_TET_0E_1V,
      HP_NONE,
    };
  int refhex_1e_1v_newels[][8] =
    {
      // { 1, 5, 2, 4 }, 
      { 1, 2, 4, 5 },
      { 7, 3, 6, 8 },
      { 2, 8, 5, 6 },
      { 2, 8, 6, 3 },
      { 2, 8, 3, 4 },
      { 2, 8, 4, 5 },
    };
  HPRef_Struct refhex_1e_1v =
    {
      HP_HEX,
      refhex_1e_1v_splitedges, 
      0, 0,
      refhex_1e_1v_newelstypes, 
      refhex_1e_1v_newels
    };


// Refine Dummies 
  // HP_HEX_3E_0V
  int refhex_3e_0v_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refhex_3e_0v_newelstypes[] =
    {
      HP_TET_1E_1VA,
      HP_TET_1E_1VA,
      HP_TET_1E_1VA,
      HP_TET_0E_1V,
      HP_TET,
      HP_NONE,
    };
  int refhex_3e_0v_newels[][8] =
    {
      { 1, 2, 3, 6 },
      { 1, 4, 8, 3 },
      { 1, 5, 6, 8 },
      { 1, 6, 3, 8 },
      { 3, 8, 6, 7 },
    };
  HPRef_Struct refhex_3e_0v =
    {
      HP_HEX,
      refhex_3e_0v_splitedges, 
      0, 0,
      refhex_3e_0v_newelstypes, 
      refhex_3e_0v_newels
    };



// Refine Dummies 
  // HP_HEX_1E_0V 
  int refhex_1e_0v_splitedges[][3] =
    {
      { 0, 0, 0 }
    };

  HPREF_ELEMENT_TYPE refhex_1e_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE,         // HP_PRISM_SINGEDGE_H1,
      HP_PRISM,
      HP_NONE,
    };
  int refhex_1e_0v_newels[][8] =
    {
      { 1, 4, 5, 2, 3, 6 },
      { 5, 4, 8, 6, 3, 7 },
    };
  HPRef_Struct refhex_1e_0v =
    {
      HP_HEX,
      refhex_1e_0v_splitedges, 
      0, 0,
      refhex_1e_0v_newelstypes, 
      refhex_1e_0v_newels
    };


