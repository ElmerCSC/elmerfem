  // HP_SEGM
  int refsegm_splitedges[][3] =
    {
      { 0, 0, 0 }
    };
  
  HPREF_ELEMENT_TYPE refsegm_newelstypes[] =
    {
      HP_SEGM,
      HP_NONE,
    };
  int refsegm_newels[][8] =
    {
      { 1, 2 },
    };
  HPRef_Struct refsegm =
    {
      HP_SEGM, 
      refsegm_splitedges, 
      0, 0, 
      refsegm_newelstypes, 
      refsegm_newels
    };

  // HP_SEGM_SINGCORNERL = 2,
  int refsegm_scl_splitedges[][3] =
    {
      { 1, 2, 3 }, 
      { 0, 0, 0 }
    };
  
  HPREF_ELEMENT_TYPE refsegm_scl_newelstypes[] = 
    {
      HP_SEGM_SINGCORNERL,
      HP_SEGM,
      HP_NONE,
    };
  
  int refsegm_scl_newels[][8] =
    {
      { 1, 3 },
      { 3, 2 },
      { 0, 0 },
    };
  HPRef_Struct refsegm_scl =
    {
      HP_SEGM,
      refsegm_scl_splitedges,
      0, 0,
      refsegm_scl_newelstypes,
      refsegm_scl_newels
    };



  // HP_SEGM_SINGCORNERR
  int refsegm_scr_splitedges[][3] =
    {
      { 2, 1, 3 },
      { 0, 0, 0 }
    };

  HPREF_ELEMENT_TYPE refsegm_scr_newelstypes[] =
    {
      HP_SEGM,
      HP_SEGM_SINGCORNERR,
      HP_NONE,
    };
  int refsegm_scr_newels[][8] =
    {
      { 1, 3 },
      { 3, 2 },
      { 0, 0 },
    };
  HPRef_Struct refsegm_scr =
    {
      HP_SEGM,
      refsegm_scr_splitedges,
      0, 0,
      refsegm_scr_newelstypes,
      refsegm_scr_newels
    };






  // HP_SEGM_SINGCORNERS = 3,
  int refsegm_sc2_splitedges[][3] =
    {
      { 1, 2, 3 },
      { 2, 1, 4 },
      { 0, 0, 0 }
    };

  HPREF_ELEMENT_TYPE refsegm_sc2_newelstypes[] =
    {
      HP_SEGM_SINGCORNERL,
      HP_SEGM_SINGCORNERR,
      HP_SEGM,
      HP_NONE,
    };
  int refsegm_sc2_newels[][8] =
    {
      { 1, 3 },
      { 4, 2 },
      { 3, 4 },
      { 0, 0 },
    };
  HPRef_Struct refsegm_sc2 =
    {  
      HP_SEGM,
      refsegm_sc2_splitedges,
      0, 0, 
      refsegm_sc2_newelstypes,
      refsegm_sc2_newels
    };




