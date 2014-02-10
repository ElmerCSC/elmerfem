
  // HP_PRISM  ... no refinement
  int refprism_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_newelstypes[] =
    {
      HP_PRISM,
      HP_NONE,
    };
  int refprism_newels[][8] =
    {
      { 1, 2, 3, 4, 5, 6 }
    };
  HPRef_Struct refprism =
    {
      HP_PRISM,
      refprism_splitedges, 
      0, 0,
      refprism_newelstypes, 
      refprism_newels
    };



  // HP_PRISM_SINGEDGE  ... vertical edge 1-4 is singular
  int refprism_singedge_splitedges[][3] =
    {
      { 1, 2, 7 },
      { 1, 3, 8 },
      { 4, 5, 9 },
      { 4, 6, 10 },
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_singedge_newelstypes[] =
    {
      HP_PRISM_SINGEDGE,
      HP_HEX,
      HP_NONE,
    };
  int refprism_singedge_newels[][8] =
    {
      { 1, 7, 8, 4, 9, 10 },
      { 3, 8, 7, 2, 6, 10, 9, 5 }
    };
  HPRef_Struct refprism_singedge =
    {
      HP_PRISM,
      refprism_singedge_splitedges, 
      0, 0,
      refprism_singedge_newelstypes, 
      refprism_singedge_newels
    };






  // HP_PRISM_SINGEDGE_V12  vertical edges 1-4 and 2-5 are singular 
  int refprism_singedge_v12_splitedges[][3] =
    {
      { 1, 2, 7 },
      { 1, 3, 8 },
      { 2, 1, 9 },
      { 2, 3, 10 },
      { 4, 5, 11 },
      { 4, 6, 12 },
      { 5, 4, 13 },
      { 5, 6, 14},
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_singedge_v12_newelstypes[] =
    {
      HP_HEX,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_singedge_v12_newels[][8] =
    {
      { 7, 9, 10, 8, 11, 13, 14, 12 },
      { 1, 7, 8, 4, 11, 12 },
      { 2, 10, 9, 5, 14, 13 },
      { 3, 8, 10, 6, 12, 14 },
    };
  HPRef_Struct refprism_singedge_v12 =
    {
      HP_PRISM,
      refprism_singedge_v12_splitedges, 
      0, 0,
      refprism_singedge_v12_newelstypes, 
      refprism_singedge_v12_newels
    };






  // HP_PRISM_SINGEDGE_H12
  int refprism_singedge_h12_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 1, 8 },
      { 2, 3, 9 },
      { 3, 1, 10 },

      { 4, 6, 12 },
      { 5, 4, 13 },
      { 5, 6, 14 },
      { 6, 4, 15 },

      { 0, 0, 0 }
    };

  int refprism_singedge_h12_splitfaces[][4] =
    {
      { 2, 1, 3, 11 },
      { 5, 4, 6, 16 },
      { 0, 0, 0, 0 },
    };

  HPREF_ELEMENT_TYPE refprism_singedge_h12_newelstypes[] =
    {
      HP_HEX,
      HP_HEX,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_singedge_h12_newels[][8] =
    {
      { 1, 8, 11, 7, 4, 13, 16, 12 },
      { 9, 3, 10, 11, 14, 6, 15, 16 },
      { 7, 11, 10, 12, 16, 15 },
      { 2, 9, 11, 5, 14, 16 },
      { 8, 2, 11, 13, 5, 16 }
    };
  HPRef_Struct refprism_singedge_h12 =
    {
      HP_PRISM,
      refprism_singedge_h12_splitedges, 
      refprism_singedge_h12_splitfaces, 
      0,
      refprism_singedge_h12_newelstypes, 
      refprism_singedge_h12_newels
    };






  // HP_PRISM_SINGEDGE_H1
  int refprism_singedge_h1_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 4, 6, 9 },
      { 5, 6, 10 },
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_singedge_h1_newelstypes[] =
    {
      HP_HEX,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_singedge_h1_newels[][8] =
    {
      { 1, 2, 8, 7, 4, 5, 10, 9 },
      { 3, 7, 8, 6, 9, 10 }
    };
  HPRef_Struct refprism_singedge_h1 =
    {
      HP_PRISM,
      refprism_singedge_h1_splitedges, 
      0, 0,
      refprism_singedge_h1_newelstypes, 
      refprism_singedge_h1_newels
    };



//  HP_PRISM_1FA_0E_0V
  int refprism_1fa_0e_0v_splitedges[][3] =
    {
      { 1, 4, 16 },
      { 2, 5, 17 },
      { 3, 6, 18 },
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fa_0e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_PRISM_1FA_0E_0V,
      HP_NONE,
    };
  int refprism_1fa_0e_0v_newels[][8] =
    {
      { 16, 17, 18, 4, 5, 6 },      
      { 1, 2, 3, 16, 17, 18 }
    };
  HPRef_Struct refprism_1fa_0e_0v =
    {
      HP_PRISM,
      refprism_1fa_0e_0v_splitedges, 
      0, 0,
      refprism_1fa_0e_0v_newelstypes, 
      refprism_1fa_0e_0v_newels
    };

//  HP_PRISM_1FA_1E_0V
  int refprism_1fa_1e_0v_splitedges[][3] =
    {
      { 1, 4, 16 },
      { 2, 5, 17 },
      { 3, 6, 18 },
      { 1, 2, 7},
      { 1, 3, 12},
      { 4, 6, 45},
      { 4, 5, 40},
      { 0, 0, 0 }
    };
  int refprism_1fa_1e_0v_splitfaces[][4] =
    {
      {1,2,4,19},
      {1,3,4,24},
      {0,0,0,0}
    }; 

  HPREF_ELEMENT_TYPE refprism_1fa_1e_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE,
      HP_HEX,
      HP_PRISM_1FA_1E_0V,
      HP_HEX_1F_0E_0V,
      HP_NONE,
    };
  int refprism_1fa_1e_0v_newels[][8] =
    {
      { 16, 19, 24, 4, 40, 45 }, 
      { 24, 19,  17, 18, 45 , 40, 5, 6 }, 
      { 1, 7 , 12 , 16, 19, 24 }, 
      { 7, 2, 3, 12,  19, 17, 18, 24 }
    };
  HPRef_Struct refprism_1fa_1e_0v =
    {
      HP_PRISM,
      refprism_1fa_1e_0v_splitedges, 
      refprism_1fa_1e_0v_splitfaces, 
      0,
      refprism_1fa_1e_0v_newelstypes, 
      refprism_1fa_1e_0v_newels
    };

//  HP_PRISM_2FA_1E_0V
  int refprism_2fa_1e_0v_splitedges[][3] =
    {
      { 1, 4, 16 },
      { 2, 5, 17 },
      { 3, 6, 18 },
      { 1, 2, 7}, 
      { 1, 3, 12},
      { 4, 6, 45},
      { 4, 5, 40},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 0, 0, 0 }
    };
  int refprism_2fa_1e_0v_splitfaces[][4] =
    {
      {1,2,4,19},
      {1,3,4,24},
      {4,1,5,31},
      {4,1,6,36},
      {0,0,0,0}
    }; 

  HPREF_ELEMENT_TYPE refprism_2fa_1e_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE,
      HP_HEX,
      HP_PRISM_1FA_1E_0V,
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_NONE,
    };
  int refprism_2fa_1e_0v_newels[][8] =
    {
      { 16, 19, 24, 28, 31, 36 }, 
      { 24, 19,  17, 18, 36, 31, 29, 30 }, 
      { 1, 7 , 12 , 16, 19, 24 }, 
      { 12, 7, 2, 3, 24, 19, 17, 18 },
      { 4, 45, 40, 28, 36, 31 }, 
      { 40, 45, 6, 5, 31, 36, 30, 29,}
    };
  HPRef_Struct refprism_2fa_1e_0v =
    {
      HP_PRISM,
      refprism_2fa_1e_0v_splitedges, 
      refprism_2fa_1e_0v_splitfaces,
      0,
      refprism_2fa_1e_0v_newelstypes, 
      refprism_2fa_1e_0v_newels
    };

//  HP_PRISM_1FB_0E_0V   ... quad face 1-2-4-5
  int refprism_1fb_0e_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 4, 6, 9 },
      { 5, 6, 10 },
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fb_0e_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_1fb_0e_0v_newels[][8] =
    {
      { 1, 4, 5, 2, 7, 9, 10, 8  },
      { 7, 8, 3, 9, 10, 6 }
    };
  HPRef_Struct refprism_1fb_0e_0v =
    {
      HP_PRISM,
      refprism_1fb_0e_0v_splitedges, 

      0, 0,
      refprism_1fb_0e_0v_newelstypes, 
      refprism_1fb_0e_0v_newels
    };


//  HP_PRISM_1FB_1EA_0V   ... quad face 1-2-4-5
  int refprism_1fb_1ea_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 4, 6, 9 },
      { 5, 6, 10 },
      { 1, 2, 11 },
      { 4, 5, 12 },
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fb_1ea_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_1fb_1ea_0v_newels[][8] =
    {
      { 11, 12, 5, 2, 7, 9, 10, 8  },
      { 1, 11, 7, 4, 12, 9 },
      { 7, 8, 3, 9, 10, 6 }
    };
  HPRef_Struct refprism_1fb_1ea_0v =
    {
      HP_PRISM,
      refprism_1fb_1ea_0v_splitedges, 
      0, 0,
      refprism_1fb_1ea_0v_newelstypes, 
      refprism_1fb_1ea_0v_newels
    };

//  HP_PRISM_1FB_1EC_0V   ... quad face 1-2-4-5 with singular edge 3-6 
  int refprism_1fb_1ec_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {3,2,10},
      {3,1,11},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE, 
      HP_HEX, 
      HP_HEX_1F_0E_0V,
      HP_NONE,
    };
  int refprism_1fb_1ec_0v_newels[][8] =
    {
      { 3, 11, 10, 6, 44, 43},
      { 12, 9, 10, 11, 45, 42, 43, 44}, 
      { 4, 5, 2, 1, 45, 42, 9, 12 } 
    };
  HPRef_Struct refprism_1fb_1ec_0v =
    {
      HP_PRISM,
      refprism_1fb_1ec_0v_splitedges, 
      0, 0,
      refprism_1fb_1ec_0v_newelstypes, 
      refprism_1fb_1ec_0v_newels
    };

//  HP_PRISM_1FA_1FB_1EC_0V   ... bot-trig face, quad face 1-2-4-5 with singular edge 3-6 
  int refprism_1fa_1fb_1ec_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {3,2,10},
      {3,1,11},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      { 0, 0, 0 }
    };

 int refprism_1fa_1fb_1ec_0v_splitfaces[][4] =
   {
     {2,3,5,21},
     {3,2,6,22},
     {3,1,6,23},
     {1,3,4,24},
     {0,0,0,0}
   }; 
  HPREF_ELEMENT_TYPE refprism_1fa_1fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE, 
      HP_HEX, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_NONE,
    };
  int refprism_1fa_1fb_1ec_0v_newels[][8] =
    {
      { 18, 23, 22, 6, 44, 43},
      { 24, 21, 22, 23, 45, 42, 43, 44}, 
      { 4, 5, 17, 16, 45, 42, 21, 24}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 12, 9, 10, 11, 24, 21, 22, 23},
      { 1, 2, 9, 12, 16, 17, 21, 24}
    };
  HPRef_Struct refprism_1fa_1fb_1ec_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_1ec_0v_splitedges,
      refprism_1fa_1fb_1ec_0v_splitfaces, 0,
      refprism_1fa_1fb_1ec_0v_newelstypes, 
      refprism_1fa_1fb_1ec_0v_newels
    };


//  HP_PRISM_1FA_1FB_2EB_0V  
  int refprism_1fa_1fb_2eb_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {3,2,10},
      {3,1,11},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      { 4, 5, 40},
      { 4, 6, 45},
      { 1, 2, 7},
      { 0, 0, 0 }
    };

 int refprism_1fa_1fb_2eb_0v_splitfaces[][4] =
   {
     {2,3,5,21},
     {3,2,6,22},
     {3,1,6,23},
     {1,3,4,24},
     {1,2,4,19},
     {0,0,0,0}
   }; 
  HPREF_ELEMENT_TYPE refprism_1fa_1fb_2eb_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE, 
      HP_HEX, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_NONE,
    };
  int refprism_1fa_1fb_2eb_0v_newels[][8] =
    {
      { 18, 23, 22, 6, 44, 43},
      { 24, 21, 22, 23, 45, 42, 43, 44}, 
      { 40, 5, 17, 19, 45, 42, 21, 24}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 12, 9, 10, 11, 24, 21, 22, 23},
      { 7, 2, 9, 12, 19, 17, 21, 24},
      {16,19,24,4,40,45},
      {1,7,12,16,19,24}
    };
  HPRef_Struct refprism_1fa_1fb_2eb_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_2eb_0v_splitedges, 
      refprism_1fa_1fb_2eb_0v_splitfaces, 0,
      refprism_1fa_1fb_2eb_0v_newelstypes, 
      refprism_1fa_1fb_2eb_0v_newels
    };

 //  HP_PRISM_1FA_1FB_2EC_0V 
  int refprism_1fa_1fb_2ec_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {3,2,10},
      {3,1,11},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,4,41},
      {2,1,8},
      { 0, 0, 0 }
    };

 int refprism_1fa_1fb_2ec_0v_splitfaces[][4] =
   {
     {2,3,5,21},
     {3,2,6,22},
     {3,1,6,23},
     {1,3,4,24},
     {2,1,5,20},
     {0,0,0,0}
   }; 
  HPREF_ELEMENT_TYPE refprism_1fa_1fb_2ec_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE, 
      HP_HEX, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_NONE,
    };
  int refprism_1fa_1fb_2ec_0v_newels[][8] =
    {
      { 18, 23, 22, 6, 44, 43},
      { 24, 21, 22, 23, 45, 42, 43, 44}, 
      { 4, 41, 20, 16, 45, 42, 21, 24}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 12, 9, 10, 11, 24, 21, 22, 23},
      { 1, 8, 9, 12, 16, 20, 21, 24}, 
      {8,2,9,20,17,21}, 
      {5,41,42,17,20,21}
    };
  HPRef_Struct refprism_1fa_1fb_2ec_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_2ec_0v_splitedges, 
      refprism_1fa_1fb_2ec_0v_splitfaces,
      0,
      refprism_1fa_1fb_2ec_0v_newelstypes, 
      refprism_1fa_1fb_2ec_0v_newels
    };







//  HP_PRISM_2FA_1FB_1EC_0V   ... trig faces, quad face 1-2-4-5 with singular edge 3-6 
  int refprism_2fa_1fb_1ec_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {3,2,10},
      {3,1,11},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 0, 0, 0 }
    };

 int refprism_2fa_1fb_1ec_0v_splitfaces[][4] =
   {
     {2,3,5,21},
     {3,2,6,22},
     {3,1,6,23},
     {1,3,4,24},
     {5,2,6,33},
     {6,5,3,34},
     {6,4,3,35},
     {4,1,6,36},
     {0,0,0,0}
   }; 
  HPREF_ELEMENT_TYPE refprism_2fa_1fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE, 
      HP_HEX, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_NONE,
    };
  int refprism_2fa_1fb_1ec_0v_newels[][8] =
    {
      { 18, 23, 22, 30, 35, 34},
      { 24, 21, 22, 23, 36, 33, 34, 35}, 
      { 28, 29, 17, 16, 36, 33, 21, 24}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 12, 9, 10, 11, 24, 21, 22, 23},
      { 1, 2, 9, 12, 16, 17, 21, 24},
      { 6, 43, 44, 30, 34, 35},
      { 44, 43, 42, 45, 35, 34, 33, 36}, 
      { 5, 4, 45, 42, 29, 28, 36, 33 }, 
    };
  HPRef_Struct refprism_2fa_1fb_1ec_0v =
    {
      HP_PRISM,
      refprism_2fa_1fb_1ec_0v_splitedges, 
      refprism_2fa_1fb_1ec_0v_splitfaces,
      0,
      refprism_2fa_1fb_1ec_0v_newelstypes, 
      refprism_2fa_1fb_1ec_0v_newels
    };

//  HP_PRISM_2FA_1FB_2EB_0V  
  int refprism_2fa_1fb_2eb_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {3,2,10},
      {3,1,11},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      {4,5,40},
      {1,2,7},
      { 0, 0, 0 }
    };

 int refprism_2fa_1fb_2eb_0v_splitfaces[][4] =
   {
     {2,3,5,21},
     {3,2,6,22},
     {3,1,6,23},
     {1,3,4,24},
     {5,6,2,33},
     {6,5,3,34},
     {6,4,3,35},
     {4,1,6,36},
     {4,1,5,31},
     {1,2,4,19},
     {0,0,0,0}
   }; 
  HPREF_ELEMENT_TYPE refprism_2fa_1fb_2eb_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE, 
      HP_HEX, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FA_1E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_NONE,
    };
  int refprism_2fa_1fb_2eb_0v_newels[][8] =
    {
      { 18, 23, 22, 30, 35, 34},
      { 24, 21, 22, 23, 36, 33, 34, 35}, 
      { 31, 29, 17, 19, 36, 33, 21, 24}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 12, 9, 10, 11, 24, 21, 22, 23},
      { 7, 2, 9, 12, 19, 17, 21, 24},
      { 6, 43, 44, 30, 34, 35},
      { 44, 43, 42, 45, 35, 34, 33, 36}, 
      { 5, 40, 45, 42, 29, 31, 36, 33 }, 
      { 1, 7, 12, 16, 19, 24 }, 
      { 16, 19, 24, 28, 31, 36 }, 
      { 40, 4, 45, 31, 28, 36 }, 
    };
  HPRef_Struct refprism_2fa_1fb_2eb_0v =
    {
      HP_PRISM,
      refprism_2fa_1fb_2eb_0v_splitedges, 
      refprism_2fa_1fb_2eb_0v_splitfaces, 0,
      refprism_2fa_1fb_2eb_0v_newelstypes, 
      refprism_2fa_1fb_2eb_0v_newels
    };

//  HP_PRISM_1FB_2EA_0V   ... quad face 1-2-4-5 with singular edges 1-4, 2-5
  int refprism_1fb_2ea_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 1, 2, 9 }, 
      { 2, 1, 10 }, 
      { 4, 6, 11 }, 
      { 5, 6, 12 }, 
      { 4, 5, 13 }, 
      { 5, 4, 14 }, 
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fb_2ea_0v_newelstypes[] =
    {
      HP_PRISM, 
      HP_PRISM_1FB_1EA_0V,
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V, 
      HP_NONE,
    };
  int refprism_1fb_2ea_0v_newels[][8] =
    {
      { 7, 8, 3, 11, 12, 6 }, 
      { 1, 9, 7, 4, 13, 11 }, 
      { 13, 14, 10, 9, 11, 12, 8, 7 }, 
      { 5, 14, 12, 2, 10, 8 }, 
    };
  HPRef_Struct refprism_1fb_2ea_0v =
    {
      HP_PRISM,
      refprism_1fb_2ea_0v_splitedges, 
      0, 0,
      refprism_1fb_2ea_0v_newelstypes, 
      refprism_1fb_2ea_0v_newels
    };

//  HP_PRISM_1FB_2EB_0V   ... quad face 1-2-4-5 with singular edges 1-4, 3-6 
  int refprism_1fb_2eb_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 4, 5, 40},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
HPREF_ELEMENT_TYPE refprism_1fb_2eb_0v_newelstypes[] =
    {
      HP_PRISM_SINGEDGE,
      HP_HEX, 
      HP_PRISM_1FB_1EA_0V, 
      HP_HEX_1F_0E_0V, 
      HP_NONE,
    };
  int refprism_1fb_2eb_0v_newels[][8] =
    {
      { 3, 11, 10, 6, 44, 43 },  
      { 12, 9, 10, 11, 45, 42, 43, 44}, 
      { 1, 7, 12, 4, 40, 45}, 
      { 40, 5, 2, 7, 45, 42, 9, 12}
    };
  HPRef_Struct refprism_1fb_2eb_0v =
    {
      HP_PRISM,
      refprism_1fb_2eb_0v_splitedges, 
      0, 0,
      refprism_1fb_2eb_0v_newelstypes, 
      refprism_1fb_2eb_0v_newels
    };

//  HP_PRISM_1FB_3E_0V   ... quad face 1-2-4-5 with singular edges 1-4, 3-6
  int refprism_1fb_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fb_3e_0v_newelstypes[] =
    { 
      HP_PRISM_SINGEDGE,
      HP_HEX, 
      HP_PRISM_1FB_1EA_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V, 
      HP_NONE,
    };
  int refprism_1fb_3e_0v_newels[][8] =
    {
      { 3, 11, 10, 6, 44, 43 }, 
      { 12, 9, 10, 11, 45, 42, 43, 44}, 
      { 1, 7, 12, 4, 40, 45 }, 
      { 40, 41, 8, 7, 45, 42, 9, 12}, 
      { 5, 41, 42, 2, 8, 9}, 
    };
  HPRef_Struct refprism_1fb_3e_0v =
    {
      HP_PRISM,
      refprism_1fb_3e_0v_splitedges, 
      0, 0,
      refprism_1fb_3e_0v_newelstypes, 
      refprism_1fb_3e_0v_newels
    };



//  HP_PRISM_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3
  int refprism_2fb_0e_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 1, 2, 9 }, 
      { 3, 2, 10 }, 
      { 4, 6, 11 },
      { 5, 6, 12 },
      { 4, 5, 13 },
      { 6, 5, 14 },
      { 0, 0, 0 }
    };
 int refprism_2fb_0e_0v_splitfaces[][4] =
    {
      { 1, 2, 3, 15 },
      { 4, 5, 6, 16 },
      { 0, 0, 0, 0 },
    };
  HPREF_ELEMENT_TYPE refprism_2fb_0e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_NONE,
    };
  int refprism_2fb_0e_0v_newels[][8] =
    {
      { 15, 8, 10, 16, 12, 14 }, 
      { 13, 5, 2, 9, 16, 12, 8, 15}, 
      { 11, 7, 3, 6, 16, 15, 10, 14 }, 
      { 1, 9, 15, 4, 13, 16 }, 
      { 4, 11, 16, 1,7, 15 }
    };
  HPRef_Struct refprism_2fb_0e_0v =
    {
      HP_PRISM,
      refprism_2fb_0e_0v_splitedges, 
      refprism_2fb_0e_0v_splitfaces,
      0,
      refprism_2fb_0e_0v_newelstypes, 
      refprism_2fb_0e_0v_newels
    };

//  HP_PRISM_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3 and sing edge 3-6
  int refprism_2fb_1ec_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 1, 2, 9 }, 
      { 3, 2, 10 }, 
      { 4, 6, 11 },
      { 5, 6, 12 },
      { 4, 5, 13 },
      { 6, 5, 14 },
      { 3, 1, 17},
      { 6, 4, 18}, 
      { 0, 0, 0 }
    };
 int refprism_2fb_1ec_0v_splitfaces[][4] =
    {
      { 1, 2, 3, 15 },
      { 4, 5, 6, 16 },
      { 0, 0, 0, 0 },
    };
  HPREF_ELEMENT_TYPE refprism_2fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_NONE,
    };
  int refprism_2fb_1ec_0v_newels[][8] =
    {
      { 15, 8, 10, 16, 12, 14 }, 
      { 13, 5, 2, 9, 16, 12, 8, 15}, 
      { 11, 7, 17, 18, 16, 15, 10, 14 }, 
      { 1, 9, 15, 4, 13, 16 }, 
      { 4, 11, 16, 1,7, 15 }, 
      { 3, 17, 10, 6, 18, 14 } 
    };
  HPRef_Struct refprism_2fb_1ec_0v =
    {
      HP_PRISM,
      refprism_2fb_1ec_0v_splitedges, 
      refprism_2fb_1ec_0v_splitfaces,
      0,
      refprism_2fb_1ec_0v_newelstypes, 
      refprism_2fb_1ec_0v_newels
    };



//  HP_PRISM_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3 and 3 sing edges
  int refprism_2fb_3e_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 1, 2, 9 }, 
      { 3, 2, 10 }, 
      { 4, 6, 11 },
      { 5, 6, 12 },
      { 4, 5, 13 },
      { 6, 5, 14 },
      { 3, 1, 17},
      { 6, 4, 18}, 
      { 2, 1, 19},
      { 5, 4, 20}, 
      { 0, 0, 0 }
    };
 int refprism_2fb_3e_0v_splitfaces[][4] =
    {
      { 1, 2, 3, 15 },
      { 4, 5, 6, 16 },
      { 0, 0, 0, 0 },
    };
  HPREF_ELEMENT_TYPE refprism_2fb_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_NONE,
    };
  int refprism_2fb_3e_0v_newels[][8] =
    {
      { 15, 8, 10, 16, 12, 14 }, 
      { 13, 20, 19, 9, 16, 12, 8, 15}, 
      { 11, 7, 17, 18, 16, 15, 10, 14 }, 
      { 1, 9, 15, 4, 13, 16 }, 
      { 4, 11, 16, 1,7, 15 }, 
      { 3, 17, 10, 6, 18, 14 }, 
      { 5, 20, 12, 2, 19, 8 }
    };
  HPRef_Struct refprism_2fb_3e_0v =
    {
      HP_PRISM,
      refprism_2fb_3e_0v_splitedges, 
      refprism_2fb_3e_0v_splitfaces, 0,
      refprism_2fb_3e_0v_newelstypes, 
      refprism_2fb_3e_0v_newels
    };



//  HP_PRISM_1FA_1FB_0E_0V   ... quad face 1-2-4-5 and trig face 1-2-3
  int refprism_1fa_1fb_0e_0v_splitedges[][3] = 
    {
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {2,3,9},
      {1,3,12},
      {5,6,42},
      {4,6,45},
      {0,0,0}
    };
  int refprism_1fa_1fb_0e_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      { 0, 0, 0, 0 }
    };

HPREF_ELEMENT_TYPE refprism_1fa_1fb_0e_0v_newelstypes[] =
    {
      HP_PRISM, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_NONE,
    };
  int refprism_1fa_1fb_0e_0v_newels[][8] =
    {
      { 24, 21, 18, 45, 42, 6 }, 
      { 4, 5, 17, 16, 45, 42, 21, 24 },
      { 12, 9, 3, 24, 21, 18 }, 
      { 1, 2, 9, 12, 16, 17, 21, 24 } 
    };
  HPRef_Struct refprism_1fa_1fb_0e_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_0e_0v_splitedges, 

      refprism_1fa_1fb_0e_0v_splitfaces, 0,
      refprism_1fa_1fb_0e_0v_newelstypes, 
      refprism_1fa_1fb_0e_0v_newels
    };

/*
//  HP_PRISM_1FA_1FB_1EC_0V   ... quad face 1-2-4-5 and trig face 1-2-3
int refprism_1fa_1fb_1ec_0v_splitedges[][3] =
    {
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {2,3,9},
      {1,3,12},
      {5,6,42},
      {4,6,45},
      {6,5,43},
      {6,4,44},
      {3,2,10},
      {3,1,11},
      {0,0,0}
    };
  int refprism_1fa_1fb_1ec_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      { 0, 0, 0, 0 }
    };

  HPREF_ELEMENT_TYPE refprism_1fa_1fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_SINGEDGE,
      HP_PRISM_1FA_1E_0V, 
      HP_PRISM_
      HP_NONE,
    };
  int refprism_1fa_1fb_0e_0v_newels[][8] =
    {
      { 24, 21, 18, 45, 42, 6 }, 
      { 4, 5, 17, 16, 45, 42, 21, 24 },
      { 12, 9, 3, 24, 21, 18 }, 
      { 1, 2, 9, 12, 16, 17, 21, 24 } 
    };
  HPRef_Struct refprism_1fa_1fb_0e_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_1ec_0v_splitedges, 

      refprism_1fa_1fb_1ec_0v_splitfaces, 0,
      refprism_1fa_1fb_1ec_0v_newelstypes, 
      refprism_1fa_1fb_1ec_0v_newels
    };


*/




//  HP_PRISM_2FA_1FB_0E_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
  int refprism_2fa_1fb_0e_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {4,1,28},
      {5,2,29},
      {6,3,30},
      {0,0,0}
      
    };
  int refprism_2fa_1fb_0e_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {5,6,2,33},
      {4,1,6,36},
      {0,0,0,0}
    };

  HPREF_ELEMENT_TYPE refprism_2fa_1fb_0e_0v_newelstypes[] =
    {  
      HP_HEX_1F_0E_0V,
      HP_PRISM,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_NONE,
    };
  int refprism_2fa_1fb_0e_0v_newels[][8] =
    {
      {28,29,17,16,36,33,21,24}, 
      {24,21,18, 36, 33, 30}, 
      {12,9,3,24,21,18},
      {1,2,9,12,16,17,21,24}, 
      {6,42,45,30,33,36},
      {4,5,29,28,45,42,33,36}
    };
  HPRef_Struct refprism_2fa_1fb_0e_0v =
    {
      HP_PRISM,
      refprism_2fa_1fb_0e_0v_splitedges, 

      refprism_2fa_1fb_0e_0v_splitfaces, 0,
      refprism_2fa_1fb_0e_0v_newelstypes, 
      refprism_2fa_1fb_0e_0v_newels
    };


//  HP_PRISM_1FA_1FB_1EA_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
  int refprism_1fa_1fb_1ea_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {4,5,40},
      {1,2,7},
      {0,0,0}, 
    };
  int refprism_1fa_1fb_1ea_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {1,2,4,19},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_1fa_1fb_1ea_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PRISM,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_NONE
    };
  int refprism_1fa_1fb_1ea_0v_newels[][8] =
    {
      {40,5,17,19,45,42,21,24}, 
      {24,21,18,45,42,6},
      {12,9,3,24,21,18},
      {7,2,9,12,19,17,21,24}, 
      {16,19,24,4,40,45},
      {1,7,12,16,19,24}
      
    };
  HPRef_Struct refprism_1fa_1fb_1ea_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_1ea_0v_splitedges, 
      refprism_1fa_1fb_1ea_0v_splitfaces, 0,
      refprism_1fa_1fb_1ea_0v_newelstypes, 
      refprism_1fa_1fb_1ea_0v_newels
    };

//  HP_PRISM_2FA_1FB_1EA_0V   
  int refprism_2fa_1fb_1ea_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {4,1,28},
      {5,2,29},
      {6,3,30},
      {4,5,40},
      {1,2,7},
      {0,0,0}, 
    };
  int refprism_2fa_1fb_1ea_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {1,2,4,19},
      {4,1,6,36},
      {4,1,5,31},
      {5,6,2,33},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_2fa_1fb_1ea_0v_newelstypes[] =
    {
      HP_PRISM, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_NONE
    };
  int refprism_2fa_1fb_1ea_0v_newels[][8] =
    {
      { 18, 24, 21, 30, 36, 33}, 
      { 31, 29, 17, 19, 36, 33, 21, 24}, 
      { 16,19, 24, 28, 31, 36 }, 
      { 3, 12, 9, 18, 24, 21 }, 
      { 7, 2, 9, 12, 19, 17, 21, 24},  
      { 1, 7, 12, 16, 19, 24 }, 
      { 6, 42, 45, 30, 33, 36 }, 
      { 40, 5, 29, 31, 45, 42, 33, 36 },
      { 40, 4, 45, 31, 28, 36} 
    };
  HPRef_Struct refprism_2fa_1fb_1ea_0v =
    {
      HP_PRISM,
      refprism_2fa_1fb_1ea_0v_splitedges, 
      refprism_2fa_1fb_1ea_0v_splitfaces, 0,
      refprism_2fa_1fb_1ea_0v_newelstypes, 
      refprism_2fa_1fb_1ea_0v_newels
    };


//  HP_PRISM_2FA_1FB_2EA_0V   
  int refprism_2fa_1fb_2ea_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {4,1,28},
      {5,2,29},
      {6,3,30},
      {4,5,40},
      {1,2,7},
      { 5, 4, 41},
      { 2, 1, 8},
      {0,0,0}, 
    };
  int refprism_2fa_1fb_2ea_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {1,2,4,19},
      {4,1,6,36},
      {4,1,5,31},
      {5,6,2,33},
      {5,4,2,32},
      {2,1,5,20},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_2fa_1fb_2ea_0v_newelstypes[] =
    {
      HP_PRISM, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_NONE
    };
  int refprism_2fa_1fb_2ea_0v_newels[][8] =
    {
      { 18, 24, 21, 30, 36, 33}, 
      { 31, 32, 20, 19, 36, 33, 21, 24}, 
      { 16,19, 24, 28, 31, 36 }, 
      { 3, 12, 9, 18, 24, 21 }, 
      {7,8,9,12,19,20,21,24},
      { 1, 7, 12, 16, 19, 24 }, 
      { 6, 42, 45, 30, 33, 36 }, 
      { 40, 41, 32, 31, 45, 42, 33, 36}, 
      { 40, 4, 45, 31, 28, 36}, 
      { 8, 2, 9, 20, 17, 21 },  
      { 29, 32, 33, 17, 20, 21 }, 
      { 5, 41, 42, 29, 32, 33 }, 
    };
  HPRef_Struct refprism_2fa_1fb_2ea_0v =
    {
      HP_PRISM,
      refprism_2fa_1fb_2ea_0v_splitedges, 
      refprism_2fa_1fb_2ea_0v_splitfaces, 0,
      refprism_2fa_1fb_2ea_0v_newelstypes, 
      refprism_2fa_1fb_2ea_0v_newels
    };

//  HP_PRISM_2FA_1FB_3E_0V   
  int refprism_2fa_1fb_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      {0,0,0}, 
    };
  int refprism_2fa_1fb_3e_0v_splitfaces[][4] = 
    {
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,1,5,31},
      {5,4,2,32},
      {5,6,2,33},
      {6,5,3,34},
      {6,4,3,35},
      {4,1,6,36},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_2fa_1fb_3e_0v_newelstypes[] =
    {
      HP_HEX,
      HP_PRISM_SINGEDGE,
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FB_1EA_0V, 

      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_HEX_1FA_1FB_0E_0V,
      
      HP_NONE
    };
  int refprism_2fa_1fb_3e_0v_newels[][8] =
    {
      {24, 21, 22, 23, 36, 33, 34, 35},
      {18, 23, 22, 30, 35, 34}, 
      { 31, 32, 20, 19, 36, 33, 21, 24}, 
      { 16,19, 24, 28, 31, 36 }, 
      { 29, 32, 33, 17, 20, 21},
      
      
      { 12, 9,10,11, 24, 21, 22, 23 }, 
      { 3, 11, 10, 18,23,22}, 
      { 1, 7, 12 , 16, 19, 24}, 
      { 8,2,9, 20, 17,21}, 
      { 7,8,9,12,19, 20, 21, 24}, 
      
      { 44, 43, 42, 45, 35, 34, 33, 36}, 
      { 6, 43, 44, 30, 34, 35}, 
      { 40, 4, 45, 31,28, 36}, 
      { 5, 41,42, 29, 32, 33},  
      { 40, 41, 32, 31, 45, 42, 33, 36}, 
    };
  HPRef_Struct refprism_2fa_1fb_3e_0v =
    {
      HP_PRISM,
      refprism_2fa_1fb_3e_0v_splitedges, 

      refprism_2fa_1fb_3e_0v_splitfaces, 0,
      refprism_2fa_1fb_3e_0v_newelstypes, 
      refprism_2fa_1fb_3e_0v_newels
    };




//  HP_PRISM_1FA_1FB_1EB_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
  int refprism_1fa_1fb_1eb_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {5,4,41},
      {2,1,8},
      {0,0,0}, 
    };
  int refprism_1fa_1fb_1eb_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {2,1,5,20},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_1fa_1fb_1eb_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PRISM,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V ,
      HP_NONE
    };
  int refprism_1fa_1fb_1eb_0v_newels[][8] =
    {
      {4,41,20,16,45,42,21,24}, 
      {24,21,18,45,42,6},
      {12,9,3,24,21,18},
      {1,8,9,12,16,20,21,24},
      {5,41,42,17,20,21},
      {8,2,9,20,17,21}
    };
  HPRef_Struct refprism_1fa_1fb_1eb_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_1eb_0v_splitedges, 

       refprism_1fa_1fb_1eb_0v_splitfaces, 0,
      refprism_1fa_1fb_1eb_0v_newelstypes, 
      refprism_1fa_1fb_1eb_0v_newels
    };


//  HP_PRISM_1FA_1FB_2EA_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
  int refprism_1fa_1fb_2ea_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {5,4,41},
      {2,1,8},
      {4,5,40},
      {1,2,7},
      {0,0,0},

    };
  int refprism_1fa_1fb_2ea_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {2,1,5,20},
      {1,2,4,19},
      {0,0,0,0},
    };

  HPREF_ELEMENT_TYPE refprism_1fa_1fb_2ea_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PRISM,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V ,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_NONE
    };
  int refprism_1fa_1fb_2ea_0v_newels[][8] =
    {
      {40,41,20,19,45,42,21,24}, 
      {24,21,18,45,42,6},
      {12,9,3,24,21,18},
      {7,8,9,12,19,20,21,24},
      {5,41,42,17,20,21},
      {8,2,9,20,17,21},
      {16,19,24,4,40,45},
      {1,7,12,16,19,24}
    };
  HPRef_Struct refprism_1fa_1fb_2ea_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_2ea_0v_splitedges, 

      refprism_1fa_1fb_2ea_0v_splitfaces, 0,
      refprism_1fa_1fb_2ea_0v_newelstypes, 
      refprism_1fa_1fb_2ea_0v_newels
    };


//  HP_PRISM_1FA_1FB_3E_0V   
  int refprism_1fa_1fb_3e_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {5,4,41},
      {2,1,8},
      {4,5,40},
      {1,2,7},
      { 3, 2, 10},
      { 3, 1, 11},
      { 6, 5, 43},
      { 6, 4, 44},
      {0,0,0},

    };
  int refprism_1fa_1fb_3e_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {2,1,5,20},
      {1,2,4,19},
      {3,2,6,22},
      {3,1,6,23},
      {0,0,0,0},
    };

  HPREF_ELEMENT_TYPE refprism_1fa_1fb_3e_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_HEX,
      HP_PRISM_SINGEDGE,
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V ,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_NONE
    };
  int refprism_1fa_1fb_3e_0v_newels[][8] =
    {
      {40,41,20,19,45,42,21,24}, 
      {24, 21, 22, 23, 45, 42, 43, 44},
      {18, 23, 22, 6, 44, 43}, 
      {12, 9, 10, 11, 24, 21, 22, 23}, 
      {3, 11, 10, 18, 23, 22}, 
      {7,8,9,12,19,20,21,24},
      {5,41,42,17,20,21},
      {8,2,9,20,17,21},
      {16,19,24,4,40,45},
      {1,7,12,16,19,24}
    };
  HPRef_Struct refprism_1fa_1fb_3e_0v =
    {
      HP_PRISM,
      refprism_1fa_1fb_3e_0v_splitedges, 

      refprism_1fa_1fb_3e_0v_splitfaces, 0,
      refprism_1fa_1fb_3e_0v_newelstypes, 
      refprism_1fa_1fb_3e_0v_newels
    };








//  HP_PRISM_2FA_0E_0V  singular trig faces
  int refprism_2fa_0e_0v_splitedges[][3] =
    {
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {4,1,28},
      {5,2,29},
      {6,3,30},
      {0,0,0}
    };
  
HPREF_ELEMENT_TYPE refprism_2fa_0e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_NONE
    };
  int refprism_2fa_0e_0v_newels[][8] =
    {
      {16,17,18,28,29,30},
      {1,2,3,16,17,18},
      {4,6,5,28,30,29}, 
    };

HPRef_Struct refprism_2fa_0e_0v = 

    {
      HP_PRISM,
      refprism_2fa_0e_0v_splitedges, 
      0, 0,
      refprism_2fa_0e_0v_newelstypes, 
      refprism_2fa_0e_0v_newels
    };





//  HP_PRISM_1FA_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3
int refprism_1fa_2fb_0e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 3, 9},
      { 3, 2, 10},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 6, 42},
      { 6, 5, 43},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
int refprism_1fa_2fb_0e_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,3,5,21},
      {3,2,6,22},
      {1,3,4,24},
      {4,5,6,46},
      { 0, 0, 0, 0 }
    };
int refprism_1fa_2fb_0e_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {0,0,0,0,0} 
  };
  
HPREF_ELEMENT_TYPE refprism_1fa_2fb_0e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_NONE,
    };
  int refprism_1fa_2fb_0e_0v_newels[][8] =
    {
      { 25, 21, 22, 46, 42, 43 }, 
      { 40, 5, 17, 19, 46, 42, 21, 25 }, 
      { 24, 18, 6, 45, 25, 22, 43, 46}, 
      { 16, 19, 25, 4, 40, 46 }, 
      { 4, 45, 46, 16, 24, 25 }, 
      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 2, 9, 13, 19, 17, 21, 25 }, 
      { 3, 12, 13, 10, 18, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 }, 
      { 12, 1, 13, 24, 16, 25 }
      
    };
  HPRef_Struct refprism_1fa_2fb_0e_0v =
    {
      HP_PRISM,
      refprism_1fa_2fb_0e_0v_splitedges, 
      refprism_1fa_2fb_0e_0v_splitfaces, 
      refprism_1fa_2fb_0e_0v_splitelement, 
      refprism_1fa_2fb_0e_0v_newelstypes, 
      refprism_1fa_2fb_0e_0v_newels
    };

//  HP_PRISM_1FA_2FB_1EC    ... quad face 1-2-4-5 and quad face 1-4-6-3
int refprism_1fa_2fb_1ec_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
int refprism_1fa_2fb_1ec_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,5,6,46},
      { 0, 0, 0, 0 }
    };
int refprism_1fa_2fb_1ec_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {0,0,0,0,0} 
  };
  
HPREF_ELEMENT_TYPE refprism_1fa_2fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V, 
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      
      HP_NONE,
    };
  int refprism_1fa_2fb_1ec_0v_newels[][8] =
    {
      { 25, 21, 22, 46, 42, 43 }, 
      { 40, 5, 17, 19, 46, 42, 21, 25 }, 
      { 24, 23, 44, 45, 25, 22, 43, 46}, 
      { 16, 19, 25, 4, 40, 46 }, 
      { 4, 45, 46, 16, 24, 25 }, 
      { 18, 23, 22, 6, 44, 43}, 


      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 2, 9, 13, 19, 17, 21, 25 }, 
      { 11, 12, 13, 10, 23, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 }, 
      { 12, 1, 13, 24, 16, 25 }, 
      { 3, 11, 10, 18, 23, 22},
      
    };
  HPRef_Struct refprism_1fa_2fb_1ec_0v =
    {
      HP_PRISM,
      refprism_1fa_2fb_1ec_0v_splitedges, 
      refprism_1fa_2fb_1ec_0v_splitfaces, 
      refprism_1fa_2fb_1ec_0v_splitelement, 
      refprism_1fa_2fb_1ec_0v_newelstypes, 
      refprism_1fa_2fb_1ec_0v_newels
    };


//  HP_PRISM_1FA_2FB_3E    ... quad face 1-2-4-5 and quad face 1-4-6-3
int refprism_1fa_2fb_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
int refprism_1fa_2fb_3e_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,5,6,46},
      { 0, 0, 0, 0 }
    };
int refprism_1fa_2fb_3e_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {0,0,0,0,0} 
  };
  
HPREF_ELEMENT_TYPE refprism_1fa_2fb_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FB_1EA_0V,
      
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V,
      
      HP_NONE,
    };
  int refprism_1fa_2fb_3e_0v_newels[][8] =
    {
      { 25, 21, 22, 46, 42, 43 }, 
      { 40, 41, 20, 19, 46, 42, 21, 25 }, 
      { 24, 23, 44, 45, 25, 22, 43, 46}, 
      { 16, 19, 25, 4, 40, 46 }, 
      { 4, 45, 46, 16, 24, 25 }, 
      { 18, 23, 22, 6, 44, 43}, 
      { 5, 41, 42, 17, 20, 21}, 
      

      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 8, 9, 13, 19, 20, 21, 25 }, 
      { 11, 12, 13, 10, 23, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 }, 
      
      { 12, 1, 13, 24, 16, 25 }, 
      { 3, 11, 10, 18, 23, 22},
      { 8, 2, 9, 20, 17, 21}, 
      
    };
  HPRef_Struct refprism_1fa_2fb_3e_0v =
    {
      HP_PRISM,
      refprism_1fa_2fb_3e_0v_splitedges, 
      refprism_1fa_2fb_3e_0v_splitfaces, 
      refprism_1fa_2fb_3e_0v_splitelement, 
      refprism_1fa_2fb_3e_0v_newelstypes, 
      refprism_1fa_2fb_3e_0v_newels
    };









//  HP_PRISM_1FA_2FB_1eb    ... quad face 1-2-4-5 and quad face 1-4-6-3
int refprism_1fa_2fb_1eb_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
int refprism_1fa_2fb_1eb_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {1,3,4,24},
      {4,5,6,46},
      { 0, 0, 0, 0 }
    };
int refprism_1fa_2fb_1eb_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {0,0,0,0,0} 
  };


HPREF_ELEMENT_TYPE refprism_1fa_2fb_1eb_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V, 
      HP_PRISM_1FB_1EA_0V,
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EB_0V,
      
      HP_NONE,
    };

  int refprism_1fa_2fb_1eb_0v_newels[][8] =
    {
      { 25, 21, 22, 46, 42, 43 }, 
      { 40, 41, 20, 19, 46, 42, 21, 25 }, 
      { 24, 18, 6, 45, 25, 22, 43, 46}, 
      { 16, 19, 25, 4, 40, 46 },
      { 4, 45, 46, 16, 24, 25 }, 
      { 5, 41, 42, 17, 20, 21 },


      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 8, 9, 13, 19, 20, 21, 25 }, 
      { 3, 12, 13, 10, 18, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 },  
      { 12, 1, 13, 24, 16, 25 }, 
      { 8, 2, 9, 20, 17, 21}, 
      
    };
  HPRef_Struct refprism_1fa_2fb_1eb_0v =
    {
      HP_PRISM,
      refprism_1fa_2fb_1eb_0v_splitedges, 
      refprism_1fa_2fb_1eb_0v_splitfaces, 
      refprism_1fa_2fb_1eb_0v_splitelement, 
      refprism_1fa_2fb_1eb_0v_newelstypes, 
      refprism_1fa_2fb_1eb_0v_newels
    };






//  HP_PRISM_2FA_2FB 
int refprism_2fa_2fb_0e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 3, 9},
      { 3, 2, 10},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 6, 42},
      { 6, 5, 43},
      { 4, 6, 45},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 0, 0, 0 }
    };
int refprism_2fa_2fb_0e_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,3,5,21},
      {3,2,6,22},
      {1,3,4,24},
      {4,5,6,46},  
      {4,1,5,31},
      {5,6,2,33},
      {6,5,3,34},
      {4,1,6,36},
      { 0, 0, 0, 0 }
    };
int refprism_2fa_2fb_0e_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {4,1,6,5,37},
    {0,0,0,0,0} 
  };
  
HPREF_ELEMENT_TYPE refprism_2fa_2fb_0e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      
      HP_NONE,
    };
  int refprism_2fa_2fb_0e_0v_newels[][8] =
    {
      { 25, 21, 22, 37, 33, 34}, 
      { 31, 29, 17, 19, 37, 33, 21, 25}, 
      { 36, 24, 18, 30, 37, 25, 22, 34}, 
      { 16, 19, 25, 28, 31, 37}, 
      { 28, 36, 37, 16, 24, 25},
      
      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 2, 9, 13, 19, 17, 21, 25 }, 
      { 3, 12, 13, 10, 18, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 }, 
      { 12, 1, 13, 24, 16, 25 }, 

      {  46, 43, 42 ,37, 34, 33},
      { 40, 5, 29, 31, 46, 42, 33, 37 }, 
      { 6, 45, 36, 30, 43, 46, 37, 34 }, 
      { 40, 4, 46, 31, 28, 37 }, 
      { 4, 45, 46, 28, 36, 37},  
      
    };
  HPRef_Struct refprism_2fa_2fb_0e_0v =
    {
      HP_PRISM,
      refprism_2fa_2fb_0e_0v_splitedges, 
      refprism_2fa_2fb_0e_0v_splitfaces, 
      refprism_2fa_2fb_0e_0v_splitelement, 
      refprism_2fa_2fb_0e_0v_newelstypes, 
      refprism_2fa_2fb_0e_0v_newels
    };


//  HP_PRISM_2FA_2FB_1EC 
int refprism_2fa_2fb_1ec_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 4, 5, 40},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
int refprism_2fa_2fb_1ec_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,5,6,46},  
      {4,1,5,31},
      {5,6,2,33},
      {6,5,3,34},
      {6,4,3,35},
      {4,1,6,36},
      { 0, 0, 0, 0 }
    };
int refprism_2fa_2fb_1ec_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {4,1,6,5,37},
    {0,0,0,0,0} 
  };
  
HPREF_ELEMENT_TYPE refprism_2fa_2fb_1ec_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V, 
      
      HP_NONE,
    };
  int refprism_2fa_2fb_1ec_0v_newels[][8] =
    {
      { 25, 21, 22, 37, 33, 34}, 
      { 31, 29, 17, 19, 37, 33, 21, 25}, 
      { 36, 24, 23, 35, 37, 25, 22, 34}, 
      { 16, 19, 25, 28, 31, 37}, 
      { 28, 36, 37, 16, 24, 25},
      { 18, 23, 22, 30, 35, 34}, 
            
      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 2, 9, 13, 19, 17, 21, 25 }, 
      { 11, 12, 13, 10, 23, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 }, 
      { 12, 1, 13, 24, 16, 25 }, 
      { 3, 11, 10, 18, 23, 22 }, 

      { 46, 43, 42 ,37, 34, 33},
      { 40, 5, 29, 31, 46, 42, 33, 37 }, 
      { 44, 45, 36, 35, 43, 46, 37, 34 }, 
      { 40, 4, 46, 31, 28, 37 }, 
      { 4, 45, 46, 28, 36, 37},  
      { 44, 6, 43, 35, 30, 34}, 
      
    };
  HPRef_Struct refprism_2fa_2fb_1ec_0v =
    {
      HP_PRISM,
      refprism_2fa_2fb_1ec_0v_splitedges, 
      refprism_2fa_2fb_1ec_0v_splitfaces, 
      refprism_2fa_2fb_1ec_0v_splitelement, 
      refprism_2fa_2fb_1ec_0v_newelstypes, 
      refprism_2fa_2fb_1ec_0v_newels
    };



//  HP_PRISM_2FA_2FB_3E 
int refprism_2fa_2fb_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0 }
    };
int refprism_2fa_2fb_3e_0v_splitfaces[][4] =
    {
      {1,2,3,13},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,5,6,46},  
      {4,1,5,31},
      {5,4,2,32},
      {5,6,2,33},
      {6,5,3,34},
      {6,4,3,35},
      {4,1,6,36},
      { 0, 0, 0, 0 }
    };
int refprism_2fa_2fb_3e_0v_splitelement[][5] = 
  {
    {1,2,3,4,25}, 
    {4,1,6,5,37},
    {0,0,0,0,0} 
  };
  
HPREF_ELEMENT_TYPE refprism_2fa_2fb_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V, 
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EB_0V,
      HP_PRISM_1FA_1FB_1EA_0V,
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      
      HP_NONE,
    };
  int refprism_2fa_2fb_3e_0v_newels[][8] =
    {
      { 25, 21, 22, 37, 33, 34}, 
      { 31, 32, 20, 19, 37, 33, 21, 25}, 
      { 36, 24, 23, 35, 37, 25, 22, 34}, 
      { 16, 19, 25, 28, 31, 37}, 
      { 28, 36, 37, 16, 24, 25},
      { 18, 23, 22, 30, 35, 34}, 
      { 29, 32, 33, 17, 20, 21}, 
            
      { 13, 9, 10, 25, 21, 22 }, 
      { 7, 8, 9, 13, 19, 20, 21, 25 }, 
      { 11, 12, 13, 10, 23, 24, 25, 22 }, 
      { 1, 7, 13, 16, 19, 25 }, 
      { 12, 1, 13, 24, 16, 25 }, 
      { 3, 11, 10, 18, 23, 22 }, 
      { 8, 2, 9, 20, 17, 21 }, 

      { 46, 43, 42 ,37, 34, 33},
      { 40, 41, 32, 31, 46, 42, 33, 37 }, 
      { 44, 45, 36, 35, 43, 46, 37, 34 }, 
      { 40, 4, 46, 31, 28, 37 }, 
      { 4, 45, 46, 28, 36, 37},  
      { 44, 6, 43, 35, 30, 34},
      { 5, 41, 42, 29, 32, 33}, 
      
    };
  HPRef_Struct refprism_2fa_2fb_3e_0v =
    {
      HP_PRISM,
      refprism_2fa_2fb_3e_0v_splitedges, 
      refprism_2fa_2fb_3e_0v_splitfaces, 
      refprism_2fa_2fb_3e_0v_splitelement, 
      refprism_2fa_2fb_3e_0v_newelstypes, 
      refprism_2fa_2fb_3e_0v_newels
    };




//  HP_PRISM_1FA_2E_0V  
  int refprism_1fa_2e_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {5,4,41},
      {2,1,8},
      {4,5,40},
      {1,2,7},
      {0,0,0},

    };
  int refprism_1fa_2e_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {2,1,5,20},
      {1,2,4,19},
      {0,0,0,0},
    };

  HPREF_ELEMENT_TYPE refprism_1fa_2e_0v_newelstypes[] =
    {
      HP_HEX,
      HP_PRISM,
      HP_PRISM_1FA_0E_0V, 
      HP_HEX_1F_0E_0V,
      HP_PRISM_SINGEDGE, 
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_SINGEDGE,
      HP_PRISM_1FA_1E_0V,
      HP_NONE
    };
  int refprism_1fa_2e_0v_newels[][8] =
    {
      {40,41,20,19,45,42,21,24}, 
      {24,21,18,45,42,6},
      {12,9,3,24,21,18},
      {9, 12, 7, 8, 21, 24, 19, 20}, 
      { 17, 21, 20, 5, 42, 41},
      {2, 9, 8, 17, 21, 20},
      {16,19,24,4,40,45},
      {1,7,12,16,19,24}
    };
  HPRef_Struct refprism_1fa_2e_0v =
    {
      HP_PRISM,
      refprism_1fa_2e_0v_splitedges, 

      refprism_1fa_2e_0v_splitfaces, 0,
      refprism_1fa_2e_0v_newelstypes, 
      refprism_1fa_2e_0v_newels
    };
    
//  HP_PRISM_2FA_2E_0V   
  int refprism_2fa_2e_0v_splitedges[][3] =
    {
      {2,3,9},
      {1,3,12},
      {1,4,16}, 
      {2,5,17},
      {3,6,18},
      {5,6,42},
      {4,6,45},
      {4,1,28},
      {5,2,29},
      {6,3,30},
      {4,5,40},
      {1,2,7},
      { 5, 4, 41},
      { 2, 1, 8},
      {0,0,0}, 
    };
  int refprism_2fa_2e_0v_splitfaces[][4] = 
    {
      {2,3,5,21},
      {1,3,4,24},
      {1,2,4,19},
      {4,1,6,36},
      {4,1,5,31},
      {5,6,2,33},
      {5,4,2,32},
      {2,1,5,20},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_2fa_2e_0v_newelstypes[] =
    {
      HP_PRISM, 
      HP_HEX, 
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE, 
      
      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_1E_0V, 
      HP_PRISM_1FA_1E_0V, 

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_1E_0V, 
      HP_PRISM_1FA_1E_0V, 
      HP_NONE,
      
    };
  int refprism_2fa_2e_0v_newels[][8] =
    {
      { 24, 21, 18, 36, 33, 30}, 
      { 19, 20, 21, 24, 31, 32, 33, 36}, 
      { 16, 19, 24, 28, 31, 36}, 
      { 17, 21, 20, 29, 33, 32}, 
      
      { 12, 9, 3, 24, 21, 18}, 
      { 7, 8, 9, 12, 19, 20, 21, 24}, 
      { 1, 7, 12, 16, 19, 24},
      { 2, 9, 8, 17, 21, 20}, 
      
      { 45, 6, 42, 36, 30, 33}, 
      { 40, 45, 42, 41, 31, 36, 33, 32}, 
      { 4, 45, 40, 28, 36, 31 }, 
      { 5, 41, 42, 29, 32, 33 },
    };
  HPRef_Struct refprism_2fa_2e_0v =
    {
      HP_PRISM,
      refprism_2fa_2e_0v_splitedges, 
      refprism_2fa_2e_0v_splitfaces, 0,
      refprism_2fa_2e_0v_newelstypes, 
      refprism_2fa_2e_0v_newels
    };


 
//  HP_PRISM_3E_0V   
  int refprism_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0}, 
    };
  int refprism_3e_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX,
      HP_HEX,
      HP_HEX,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,
      HP_NONE
    };
  int refprism_3e_0v_newels[][8] =
    {
      { 13, 14, 15, 46, 47, 48}, 
      { 7, 8, 14, 13, 40, 41, 47, 46}, 
      { 15, 14, 9, 10, 48, 47, 42, 43}, 
      { 12, 13, 15, 11, 45, 46, 48, 44}, 
      { 14, 8, 9, 47, 41, 42 }, 
      { 11, 15, 10, 44, 48, 43 }, 
      { 7, 13, 12, 40, 46, 45}, 
      { 1, 7, 12, 4, 40, 45}, 
      { 2, 9, 8, 5, 42, 41 }, 
      { 3, 11, 10, 6, 44, 43 }
    };
  HPRef_Struct refprism_3e_0v =
    {
      HP_PRISM,
      refprism_3e_0v_splitedges, 
      refprism_3e_0v_splitfaces, 0,
      refprism_3e_0v_newelstypes, 
      refprism_3e_0v_newels
    };


//  HP_PRISM_3E_0V   
int refprism_1fa_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      
      { 0, 0, 0}, 
    };
int refprism_1fa_3e_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48}, 
      {0,0,0,0}, 
    };

int refprism_1fa_3e_0v_splitelements[][5] = 
  {
      {1,2,3,4,25},
      {2,1,3,5,26},
      {3,1,2,6,27}, 
      {0,0,0,0,0},
  };

  HPREF_ELEMENT_TYPE refprism_1fa_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX,
      HP_HEX,
      HP_HEX,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_NONE
    };
int refprism_1fa_3e_0v_newels[][8] =
    {
      { 25, 26, 27, 46, 47, 48}, 
      { 19, 20, 26, 25, 40, 41, 47, 46},  
      { 27, 26, 21, 22, 48, 47, 42, 43}, 
      { 23, 24, 25, 27, 44, 45, 46, 48}, 
      { 19, 25, 24, 40, 46, 45}, 
      { 26, 20, 21, 47, 41, 42},
      { 23, 27, 22, 44, 48, 43}, 
      { 16, 19, 24, 4, 40, 45}, 
      { 17, 21, 20, 5, 42, 41}, 
      { 18, 23, 22, 6, 44, 43}, 

      { 13, 14, 15, 25, 26, 27}, 
      { 7, 8, 14, 13, 19, 20, 26, 25},
      { 15, 14, 9, 10, 27, 26, 21, 22}, 
      { 12, 13, 15, 11, 24, 25, 27, 23}, 
      { 14, 8, 9, 26, 20, 21}, 
      { 11, 15, 10, 23, 27, 22}, 
      { 7, 13 , 12, 19, 25, 24}, 
      { 2, 9, 8, 17, 21, 20}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 1, 7, 12, 16, 19, 24}, 
    };
  HPRef_Struct refprism_1fa_3e_0v =
    {
      HP_PRISM,
      refprism_1fa_3e_0v_splitedges, 
      refprism_1fa_3e_0v_splitfaces, 
      refprism_1fa_3e_0v_splitelements, 
      refprism_1fa_3e_0v_newelstypes, 
      refprism_1fa_3e_0v_newels
    };



//  HP_PRISM_2FA_3E_0V   
int refprism_2fa_3e_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0}, 
    };
int refprism_2fa_3e_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,1,5,31},
      {5,4,2,32},
      {5,6,2,33},
      {6,5,3,34},
      {6,4,3,35},
      {4,1,6,36},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48}, 
      {0,0,0,0}, 
    };

int refprism_2fa_3e_0v_splitelements[][5] = 
  {
      {1,2,3,4,25},
      {2,1,3,5,26},
      {3,1,2,6,27}, 
      {4,1,6,5,37},
      {5,2,4,6,38},
      {6,4,5,3,39}, 
      {0,0,0,0,0},
  };

  HPREF_ELEMENT_TYPE refprism_2fa_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX,
      HP_HEX,
      HP_HEX,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,

      HP_NONE
    };

  int refprism_2fa_3e_0v_newels[][8] =
    {
      { 25, 26, 27, 37, 38, 39}, 
      { 19, 20, 26, 25, 31, 32, 38, 37},  
      { 27, 26, 21, 22, 39, 38, 33, 34}, 
      { 23, 24, 25, 27, 35, 36, 37, 39}, 
      { 19, 25, 24, 31, 37, 36}, 
      { 26, 20, 21, 38, 32, 33},
      { 23, 27, 22, 35, 39, 34}, 
      { 16, 19, 24, 28, 31, 36}, 
      { 17, 21, 20, 29, 33, 32}, 
      { 18, 23, 22, 30, 35, 34}, 

      { 13, 14, 15, 25, 26, 27}, 
      { 7, 8, 14, 13, 19, 20, 26, 25},
      { 15, 14, 9, 10, 27, 26, 21, 22}, 
      { 12, 13, 15, 11, 24, 25, 27, 23}, 
      { 14, 8, 9, 26, 20, 21}, 
      { 11, 15, 10, 23, 27, 22}, 
      { 7, 13 , 12, 19, 25, 24}, 
      { 2, 9, 8, 17, 21, 20}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 1, 7, 12, 16, 19, 24}, 

      { 48, 47, 46, 39, 38, 37 }, 
      { 48, 43, 42, 47, 39, 34, 33, 38}, 
      { 45, 44, 48, 46, 36, 35, 39, 37},
      { 46, 47, 41, 40, 37, 38, 32, 31}, 
      { 47, 42, 41, 38, 33, 32}, 
      { 45, 46, 40, 36, 37, 31}, 
      { 44, 43, 48, 35, 34, 39},
      { 6, 43, 44, 30, 34, 35}, 
      { 5, 41, 42, 29, 32, 33}, 
      { 4, 45, 40, 28, 36, 31},
    };

HPRef_Struct refprism_2fa_3e_0v =
  {
    HP_PRISM,
    refprism_2fa_3e_0v_splitedges, 
    refprism_2fa_3e_0v_splitfaces, 
    refprism_2fa_3e_0v_splitelements, 
    refprism_2fa_3e_0v_newelstypes, 
    refprism_2fa_3e_0v_newels
  };



//  HP_PRISM_3FB_0V   
  int refprism_3fb_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0}, 
    };
  int refprism_3fb_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48},
      {0,0,0,0}, 
    };

  HPREF_ELEMENT_TYPE refprism_3fb_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_NONE
    };
  int refprism_3fb_0v_newels[][8] =
    {
      { 13, 14, 15, 46, 47, 48}, 
      { 8, 7, 40, 41, 14,13, 46, 47 }, 
      { 10, 9, 42, 43, 15, 14, 47, 48 }, 
      { 44, 45, 12, 11, 48, 46, 13, 15}, 
      { 1, 7, 13, 4, 40, 46 }, 
      { 4, 45, 46, 1, 12, 13}, 
      { 2, 9, 14, 5, 42, 47 }, 
      { 5, 41, 47, 2, 8, 14 }, 
      { 3, 11, 15, 6, 44, 48}, 
      { 6, 43, 48, 3, 10, 15},

    };
  HPRef_Struct refprism_3fb_0v =
    {
      HP_PRISM,
      refprism_3fb_0v_splitedges, 
      refprism_3fb_0v_splitfaces, 0,
      refprism_3fb_0v_newelstypes, 
      refprism_3fb_0v_newels
    };


//  HP_PRISM_3FB_0V   
int refprism_1fa_3fb_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0}, 
    };
int refprism_1fa_3fb_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48}, 
      {0,0,0,0}, 
    };

int refprism_1fa_3fb_0v_splitelements[][5] = 
  {
      {1,2,3,4,25},
      {2,1,3,5,26},
      {3,1,2,6,27}, 
      {0,0,0,0,0},
  };

  HPREF_ELEMENT_TYPE refprism_1fa_3fb_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      
      HP_NONE
    };
  int refprism_1fa_3fb_0v_newels[][8] =
    {
      { 25, 26, 27, 46, 47, 48}, 
      { 19, 40, 41, 20, 25, 46, 47, 26}, 
      { 22, 21, 42, 43, 27, 26, 47, 48}, 
      { 24, 23, 44, 45, 25, 27, 48, 46},
      
      { 16, 19, 25, 4, 40, 46 }, 
      { 4, 45, 46, 16, 24, 25 }, 
      { 17, 21, 26, 5, 42, 47 }, 
      { 5, 41, 47, 17, 20, 26}, 
      { 18, 23, 27, 6, 44, 48}, 
      { 6, 43, 48, 18, 22, 27},

      { 13, 14, 15, 25, 26, 27}, 
      { 7, 8, 14, 13, 19, 20, 26, 25}, 
      { 9, 10, 15, 14, 21, 22, 27, 26},
      { 11, 12, 13, 15, 23, 24, 25, 27},

      { 2, 9, 14, 17, 21, 26}, 
      { 8, 2, 14, 20, 17, 26}, 
      { 1, 7, 13, 16, 19, 25}, 
      { 12, 1, 13, 24, 16, 25 }, 
      { 3, 11, 15, 18, 23, 27 },
      { 10, 3, 15, 22, 18, 27}, 
      
      };
  HPRef_Struct refprism_1fa_3fb_0v =
    {
      HP_PRISM,
      refprism_1fa_3fb_0v_splitedges, 
      refprism_1fa_3fb_0v_splitfaces, 
      refprism_1fa_3fb_0v_splitelements, 
      refprism_1fa_3fb_0v_newelstypes, 
      refprism_1fa_3fb_0v_newels
    };
     


//  HP_PRISM_2FA_3E_0V   
int refprism_2fa_3fb_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0}, 
    };
int refprism_2fa_3fb_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,1,5,31},
      {5,4,2,32},
      {5,6,2,33},
      {6,5,3,34},
      {6,4,3,35},
      {4,1,6,36},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48}, 
      {0,0,0,0}, 
    };

int refprism_2fa_3fb_0v_splitelements[][5] = 
  {
      {1,2,3,4,25},
      {2,1,3,5,26},
      {3,1,2,6,27}, 
      {4,1,6,5,37},
      {5,2,4,6,38},
      {6,4,5,3,39}, 
      {0,0,0,0,0},
  };

  HPREF_ELEMENT_TYPE refprism_2fa_3fb_0v_newelstypes[] =
    {

      HP_PRISM,
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1FA_1FB_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V, 
      HP_HEX_1FA_1FB_0E_0V,
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 
      HP_PRISM_1FA_1FB_1EA_0V, 
      HP_PRISM_1FA_1FB_1EB_0V, 

      HP_NONE
    };
  int refprism_2fa_3fb_0v_newels[][8] =
    {
      { 25, 26, 27, 37, 38, 39}, 
      { 19, 31, 32, 20, 25, 37, 38, 26}, 
      { 33, 34, 22, 21, 38, 39, 27, 26}, 
      { 35, 36, 24, 23, 39, 37, 25, 27}, 

      { 16, 19, 25, 28, 31, 37}, 
      { 28, 36, 37, 16, 24, 25 }, 
      { 17, 21, 26, 29, 33, 38 }, 
      { 29, 32, 38, 17, 20, 26}, 
      { 18, 23, 27, 30, 35, 39}, 
      { 30, 34, 39, 18, 22, 27},

 
      { 13, 14, 15, 25, 26, 27}, 
      { 7, 8, 14, 13, 19, 20, 26, 25}, 
      { 9, 10, 15, 14, 21, 22, 27, 26},
      { 11, 12, 13, 15, 23, 24, 25, 27},

      { 2, 9, 14, 17, 21, 26}, 
      { 8, 2, 14, 20, 17, 26}, 
      { 1, 7, 13, 16, 19, 25}, 
      { 12, 1, 13, 24, 16, 25 }, 
      { 3, 11, 15, 18, 23, 27 },
      { 10, 3, 15, 22, 18, 27}, 

      
      { 48, 47, 46, 39, 38, 37 }, 
      { 44, 45, 36, 35, 48, 46, 37, 39}, 
      { 40, 41, 32, 31, 46, 47, 38, 37}, 
      { 42, 43, 34, 33, 47, 48, 39, 38}, 
      
      { 6, 43, 48, 30, 34, 39}, 
      { 44, 6, 48, 35, 30, 39}, 
      { 4, 45, 46, 28, 36, 37}, 
      { 40, 4, 46, 31, 28, 37}, 
      { 5, 41, 47, 29, 32, 38}, 
      { 42, 5, 47, 33, 29, 38},
    };

HPRef_Struct refprism_2fa_3fb_0v =
  {
    HP_PRISM,
    refprism_2fa_3fb_0v_splitedges, 
    refprism_2fa_3fb_0v_splitfaces, 
    refprism_2fa_3fb_0v_splitelements, 
    refprism_2fa_3fb_0v_newelstypes, 
    refprism_2fa_3fb_0v_newels
  };


/* 


//  HP_PRISM_3E_4EH
int refprism_3e_4eh_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0},

    };
int refprism_3e_4eh_splitfaces[][4] = 
    {
      {3,1,2,15},
      {6,4,5,48}, 
      {0,0,0,0}, 
    };

HPREF_ELEMENT_TYPE refprism_2fa_3fb_0v_newelstypes[] =
  {
    HP_PRISM, 
    HP_HEX_2EH_0V,
    HP_HEX_2EH_0V,
    HP_TET_2E,
    HP_TET_2E,
    HP_PRISM_1E_2EH_0V, 
    HP_PRISM_1E_2EH_0V, 
    HP_NONE
    };
  int refprism_2fa_3fb_0v_newels[][8] =
    {
      {15, 7, 8, 48, 40, 41 }, 
      
    };

HPRef_Struct refprism_2fa_3fb_0v =
  {
    HP_PRISM,
    refprism_2fa_3fb_0v_splitedges, 
    refprism_2fa_3fb_0v_splitfaces, 
    refprism_2fa_3fb_0v_splitelements, 
    refprism_2fa_3fb_0v_newelstypes, 
    refprism_2fa_3fb_0v_newels
  };
*/ 

/*
//  HP_PRISM_2FA_3E_0V   
int refprism_3e_4_0v_splitedges[][3] =
    {
      { 1, 2, 7},
      { 2, 1, 8},
      { 2, 3, 9},
      { 3, 2, 10},
      { 3, 1, 11},
      { 1, 3, 12},
      { 1, 4, 16}, 
      { 2, 5, 17},
      { 3, 6, 18},
      { 4, 1, 28},
      { 5, 2, 29},
      { 6, 3, 30},
      { 4, 5, 40},
      { 5, 4, 41},
      { 5, 6, 42},
      { 6, 5, 43},
      { 6, 4, 44},
      { 4, 6, 45},
      { 0, 0, 0}, 
    };
int refprism_2fa_3e_0v_splitfaces[][4] = 
    {
      {1,2,3,13},
      {2,3,1,14},
      {3,1,2,15},
      {1,2,4,19},
      {2,1,5,20},
      {2,3,5,21},
      {3,2,6,22},
      {3,1,6,23},
      {1,3,4,24},
      {4,1,5,31},
      {5,4,2,32},
      {5,6,2,33},
      {6,5,3,34},
      {6,4,3,35},
      {4,1,6,36},
      {4,5,6,46},
      {5,4,6,47},
      {6,4,5,48}, 
      {0,0,0,0}, 
    };

int refprism_2fa_3e_0v_splitelements[][5] = 
  {
      {1,2,3,4,25},
      {2,1,3,5,26},
      {3,1,2,6,27}, 
      {4,1,6,5,37},
      {5,2,4,6,38},
      {6,4,5,3,39}, 
      {0,0,0,0,0},
  };

  HPREF_ELEMENT_TYPE refprism_2fa_3e_0v_newelstypes[] =
    {
      HP_PRISM,
      HP_HEX,
      HP_HEX,
      HP_HEX,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,
      HP_PRISM_SINGEDGE,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,

      HP_PRISM_1FA_0E_0V,
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_HEX_1F_0E_0V, 
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_0E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,
      HP_PRISM_1FA_1E_0V,

      HP_NONE
    };

  int refprism_2fa_3e_0v_newels[][8] =
    {
      { 25, 26, 27, 37, 38, 39}, 
      { 19, 20, 26, 25, 31, 32, 38, 37},  
      { 27, 26, 21, 22, 39, 38, 33, 34}, 
      { 23, 24, 25, 27, 35, 36, 37, 39}, 
      { 19, 25, 24, 31, 37, 36}, 
      { 26, 20, 21, 38, 32, 33},
      { 23, 27, 22, 35, 39, 34}, 
      { 16, 19, 24, 28, 31, 36}, 
      { 17, 21, 20, 29, 33, 32}, 
      { 18, 23, 22, 30, 35, 34}, 

      { 13, 14, 15, 25, 26, 27}, 
      { 7, 8, 14, 13, 19, 20, 26, 25},
      { 15, 14, 9, 10, 27, 26, 21, 22}, 
      { 12, 13, 15, 11, 24, 25, 27, 23}, 
      { 14, 8, 9, 26, 20, 21}, 
      { 11, 15, 10, 23, 27, 22}, 
      { 7, 13 , 12, 19, 25, 24}, 
      { 2, 9, 8, 17, 21, 20}, 
      { 3, 11, 10, 18, 23, 22}, 
      { 1, 7, 12, 16, 19, 24}, 

      { 48, 47, 46, 39, 38, 37 }, 
      { 48, 43, 42, 47, 39, 34, 33, 38}, 
      { 45, 44, 48, 46, 36, 35, 39, 37},
      { 46, 47, 41, 40, 37, 38, 32, 31}, 
      { 47, 42, 41, 38, 33, 32}, 
      { 45, 46, 40, 36, 37, 31}, 
      { 44, 43, 48, 35, 34, 39},
      { 6, 43, 44, 30, 34, 35}, 
      { 5, 41, 42, 29, 32, 33}, 
      { 4, 45, 40, 28, 36, 31},
    };

HPRef_Struct refprism_2fa_3e_0v =
  {
    HP_PRISM,
    refprism_2fa_3e_0v_splitedges, 
    refprism_2fa_3e_0v_splitfaces, 
    refprism_2fa_3e_0v_splitelements, 
    refprism_2fa_3e_0v_newelstypes, 
    refprism_2fa_3e_0v_newels
  };

*/
/*

//  HP_PRISM_1FB_1EB_0V   ... quad face 1-2-4-5
  int refprism_1fb_1eb_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 3, 8 },
      { 4, 6, 9 },
      { 5, 6, 10 },
      { 2, 1, 11 },
      { 5, 4, 12 },
      { 0, 0, 0 }
    };
  HPREF_ELEMENT_TYPE refprism_1fb_1eb_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FB_1EB_0V,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_1fb_1eb_0v_newels[][8] =
    {
      { 1, 4, 12, 11, 7, 9, 10, 8  },
      { 11, 2, 8, 12, 5, 10 },
      { 7, 8, 3, 9, 10, 6 }
    };
  HPRef_Struct refprism_1fb_1eb_0v =
    {
      HP_PRISM,
      refprism_1fb_1eb_0v_splitedges, 
      0, 0,
      refprism_1fb_1eb_0v_newelstypes, 
      refprism_1fb_1eb_0v_newels
    };










  // HP_PRISM_2F_0E_0V
  int refprism_2f_0e_0v_splitedges[][3] =
    {
      { 1, 3, 7 },
      { 2, 1, 8 },
      { 2, 3, 9 },
      { 3, 1, 10 },

      { 4, 6, 12 },
      { 5, 4, 13 },
      { 5, 6, 14 },
      { 6, 4, 15 },

      { 0, 0, 0 }
    };

  int refprism_2f_0e_0v_splitfaces[][4] =
    {
      { 2, 1, 3, 11 },
      { 5, 4, 6, 16 },
      { 0, 0, 0, 0 },
    };

  HPREF_ELEMENT_TYPE refprism_2f_0e_0v_newelstypes[] =
    {
      HP_HEX_1F_0E_0V,
      HP_HEX_1F_0E_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM_1FB_1EA_0V,
      HP_PRISM,
      HP_NONE,
    };
  int refprism_2f_0e_0v_newels[][8] =
    {
      //{ 1, 8, 11, 7, 4, 13, 16, 12 },
      // { 9, 3, 10, 11, 14, 6, 15, 16 },
      { 1, 4, 13, 8, 7, 12, 16, 11 },
      { 9, 14, 6, 3, 11, 16, 15, 10 },
      { 2, 9, 11, 5, 14, 16 },
      // { 8, 2, 11, 13, 5, 16 },
      { 5, 13, 16, 2, 8, 11 },
      { 7, 11, 10, 12, 16, 15 }
    };
  HPRef_Struct refprism_2f_0e_0v =
    {
      HP_PRISM,
      refprism_2f_0e_0v_splitedges, 
      refprism_2f_0e_0v_splitfaces, 
      0,
      refprism_2f_0e_0v_newelstypes, 
      refprism_2f_0e_0v_newels
    };

*/
