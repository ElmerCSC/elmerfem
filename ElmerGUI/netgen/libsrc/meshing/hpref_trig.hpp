

// HP_TRIG
int reftrig_splitedges[][3] =
{
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_newelstypes[] =
{
  HP_TRIG,
  HP_NONE,
};
int reftrig_newels[][8] =
{
  { 1, 2, 3 },
};
HPRef_Struct reftrig =
{
  HP_TRIG, 
  reftrig_splitedges, 
  0, 0, 
  reftrig_newelstypes, 
  reftrig_newels
};



// HP_TRIG_SINGCORNER
int reftrig_singcorner_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singcorner_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_NONE,
};
int reftrig_singcorner_newels[][8] =
{
  { 1, 4, 5 },
  { 2, 3, 5, 4 },
};
HPRef_Struct reftrig_singcorner =
{
  HP_TRIG,
  reftrig_singcorner_splitedges, 
  0, 0,
  reftrig_singcorner_newelstypes, 
  reftrig_singcorner_newels
};


/*
// HP_TRIG_SINGCORNER, trigs only
int reftrig_singcorner_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singcorner_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_TRIG,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singcorner_newels[][8] =
{
  { 1, 4, 5 },
  { 4, 2, 5 },
  { 5, 2, 3 },
};
HPRef_Struct reftrig_singcorner =
{
  HP_TRIG,
  reftrig_singcorner_splitedges, 
  0, 0,
  reftrig_singcorner_newelstypes, 
  reftrig_singcorner_newels
};
*/





// HP_TRIG_SINGCORNER12
int reftrig_singcorner12_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singcorner12_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singcorner12_newels[][8] =
{
  { 1, 4, 5 },
  { 2, 7, 6 },
  { 4, 6, 7, 5 },
  { 5, 7, 3 },
};
HPRef_Struct reftrig_singcorner12 =
{
  HP_TRIG,
  reftrig_singcorner12_splitedges, 
  0, 0,
  reftrig_singcorner12_newelstypes, 
  reftrig_singcorner12_newels
};




// HP_TRIG_SINGCORNER123_2D
int reftrig_singcorner123_2D_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singcorner123_2D_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int reftrig_singcorner123_2D_newels[][8] =
{
  { 1, 4, 5 },
  { 2, 7, 6 },
  { 3, 8, 9 },
  { 4, 6, 8, 5 },
  { 6, 7, 9, 8 },
};
HPRef_Struct reftrig_singcorner123_2D =
{
  HP_TRIG,
  reftrig_singcorner123_2D_splitedges, 
  0, 0,
  reftrig_singcorner123_2D_newelstypes, 
  reftrig_singcorner123_2D_newels
};






// HP_TRIG_SINGCORNER123
int reftrig_singcorner123_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 0, 0, 0 }
};

int reftrig_singcorner123_splitfaces[][4] =
{
  { 1, 2, 3, 10 },
  { 2, 3, 1, 11 },
  { 3, 1, 2, 12 },
  { 0, 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singcorner123_newelstypes[] =
{
  HP_DUMMY_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  //  HP_TRIG_SINGCORNER,
  //  HP_TRIG,
  //  HP_TRIG_SINGCORNER,
  //  HP_TRIG,
  //  HP_TRIG_SINGCORNER,
  //  HP_TRIG,
  HP_QUAD,
  HP_QUAD,
  HP_QUAD,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singcorner123_newels[][8] =
{
  { 1, 4, 10, 5 },
  { 2, 7, 11, 6 },
  { 3, 8, 12, 9 },
  //  { 1, 4, 5 },
  //  { 5, 4, 10 },
  //  { 2, 7, 6 },
  //  { 6, 7, 11 },
  //  { 3, 8, 9 },
  //  { 8, 12, 9 },
  { 4, 6, 11, 10 },
  { 7, 9, 12, 11 },
  { 8, 5, 10, 12 },
  { 10, 11, 12 },
};
HPRef_Struct reftrig_singcorner123 =
{
  HP_TRIG,
  reftrig_singcorner123_splitedges, 
  reftrig_singcorner123_splitfaces, 
  0, 
  reftrig_singcorner123_newelstypes, 
  reftrig_singcorner123_newels
};

// HP_TRIG_SINGEDGE
int reftrig_singedge_splitedges[][3] =
{
  { 2, 3, 4 },
  { 1, 3, 5 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedge_newelstypes[] =
{
  HP_TRIG,
  HP_QUAD_SINGEDGE,
  HP_NONE,
};
int reftrig_singedge_newels[][8] =
{
  { 4, 3, 5 },
  { 1, 2, 4, 5 },
};
HPRef_Struct reftrig_singedge =
{
  HP_TRIG,
  reftrig_singedge_splitedges, 
  0, 0,
  reftrig_singedge_newelstypes, 
  reftrig_singedge_newels
};






// HP_TRIG_SINGEDGECORNER1
int reftrig_singedgecorner1_splitedges[][3] =
{
  { 1, 2, 6 },
  { 1, 3, 5 },
  { 2, 3, 4 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner1_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedgecorner1_newels[][8] =
{
  { 1, 6, 5 },
  { 6, 2, 4, 5 },
  { 5, 4, 3 },
};
HPRef_Struct reftrig_singedgecorner1 =
{
  HP_TRIG,
  reftrig_singedgecorner1_splitedges, 
  0, 0, 
  reftrig_singedgecorner1_newelstypes, 
  reftrig_singedgecorner1_newels
};








// HP_TRIG_SINGEDGECORNER2
int reftrig_singedgecorner2_splitedges[][3] =
{
  { 2, 1, 6 },
  { 1, 3, 5 },
  { 2, 3, 4 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner2_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedgecorner2_newels[][8] =
{
  { 6, 2, 4},
  { 1, 6, 4, 5 },
  { 5, 4, 3 },
};
HPRef_Struct reftrig_singedgecorner2 =
{
  HP_TRIG,
  reftrig_singedgecorner2_splitedges, 
  0, 0,
  reftrig_singedgecorner2_newelstypes, 
  reftrig_singedgecorner2_newels
};




// HP_TRIG_SINGEDGECORNER12
int reftrig_singedgecorner12_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner12_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedgecorner12_newels[][8] =
{
  { 1, 4, 5 },
  { 6, 2, 7 },
  { 4, 6, 7, 5 },
  { 5, 7, 3 },
};
HPRef_Struct reftrig_singedgecorner12 =
{
  HP_TRIG,
  reftrig_singedgecorner12_splitedges, 
  0, 0,
  reftrig_singedgecorner12_newelstypes, 
  reftrig_singedgecorner12_newels
};







// HP_TRIG_SINGEDGECORNER3
int reftrig_singedgecorner3_splitedges[][3] =
{
  { 1, 3, 4 },
  { 3, 1, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner3_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int reftrig_singedgecorner3_newels[][8] =
{
  { 1, 2, 6, 4 },
  { 4, 6, 7, 5 },
  { 3, 5, 7 },
};
HPRef_Struct reftrig_singedgecorner3 =
{
  HP_TRIG,
  reftrig_singedgecorner3_splitedges, 
  0, 0,
  reftrig_singedgecorner3_newelstypes, 
  reftrig_singedgecorner3_newels
};




// HP_TRIG_SINGEDGECORNER13
int reftrig_singedgecorner13_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 3, 6 },
  { 3, 1, 7 },
  { 3, 2, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner13_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int reftrig_singedgecorner13_newels[][8] =
{
  { 1, 4, 5 },
  { 4, 2, 6, 5 },
  { 5, 6, 8, 7 },
  { 3, 7, 8 },
};
HPRef_Struct reftrig_singedgecorner13 =
{
  HP_TRIG,
  reftrig_singedgecorner13_splitedges, 
  0, 0,
  reftrig_singedgecorner13_newelstypes, 
  reftrig_singedgecorner13_newels
};





// HP_TRIG_SINGEDGECORNER23
int reftrig_singedgecorner23_splitedges[][3] =
{
  { 1, 3, 4 },
  { 2, 1, 5 },
  { 2, 3, 6 },
  { 3, 1, 7 },
  { 3, 2, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner23_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int reftrig_singedgecorner23_newels[][8] =
{
  { 5, 2, 6 },
  { 1, 5, 6, 4 },
  { 4, 6, 8, 7 },
  { 3, 7, 8 },
};
HPRef_Struct reftrig_singedgecorner23 =
{
  HP_TRIG,
  reftrig_singedgecorner23_splitedges, 
  0, 0,
  reftrig_singedgecorner23_newelstypes, 
  reftrig_singedgecorner23_newels
};



// HP_TRIG_SINGEDGECORNER123
int reftrig_singedgecorner123_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftrig_singedgecorner123_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int reftrig_singedgecorner123_newels[][8] =
{
  { 1, 4, 5 },
  { 6, 2, 7 },
  { 4, 6, 7, 5 },
  { 5, 7, 9, 8 },
  { 3, 8, 9 },
};
HPRef_Struct reftrig_singedgecorner123 =
{
  HP_TRIG,
  reftrig_singedgecorner123_splitedges, 
  0, 0,
  reftrig_singedgecorner123_newelstypes, 
  reftrig_singedgecorner123_newels
};

// HP_TRIG_SINGEDGES
int reftrig_singedges_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 0, 0, 0 }
};
int reftrig_singedges_splitfaces[][4] =
{
  { 1, 2, 3, 8 },
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftrig_singedges_newelstypes[] =
{
  //  HP_QUAD_2E,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedges_newels[][8] =
{
  // { 1, 4, 8, 5 },
  { 1, 4, 8 },
  { 5, 1, 8 },
  { 4, 2, 6, 8 },
  { 3, 5, 8, 7 },
  { 6, 7, 8 },
};
HPRef_Struct reftrig_singedges =
{
  HP_TRIG,
  reftrig_singedges_splitedges, 
  reftrig_singedges_splitfaces, 
  0,
  reftrig_singedges_newelstypes, 
  reftrig_singedges_newels
};








// HP_TRIG_SINGEDGES2
int reftrig_singedges2_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 3, 2, 8 },
  { 0, 0, 0 }
};
int reftrig_singedges2_splitfaces[][4] =
{
  { 1, 2, 3, 9 },
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftrig_singedges2_newelstypes[] =
{
  //  HP_QUAD_2E,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedges2_newels[][8] =
{
  //  { 1, 4, 9, 5 },
  { 1, 4, 9 },
  { 5, 1, 9 },
  { 4, 6, 7, 9 },
  { 3, 5, 9, 8 },
  { 6, 2, 7 },
  { 7, 8, 9 },
};
HPRef_Struct reftrig_singedges2 =
{
  HP_TRIG,
  reftrig_singedges2_splitedges, 
  reftrig_singedges2_splitfaces, 
  0,
  reftrig_singedges2_newelstypes, 
  reftrig_singedges2_newels
};




// HP_TRIG_SINGEDGES3
int reftrig_singedges3_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 3, 6 },
  { 3, 1, 7 },
  { 3, 2, 8 },
  { 0, 0, 0 }
};
int reftrig_singedges3_splitfaces[][4] =
{
  { 1, 2, 3, 9 },
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftrig_singedges3_newelstypes[] =
{
  //  HP_QUAD_2E,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedges3_newels[][8] =
{
  //  { 1, 4, 9, 5 },
  { 1, 4, 9 },
  { 5, 1, 9 },
  { 4, 2, 6, 9 },
  { 7, 5, 9, 8 },
  { 3, 7, 8 },
  { 6, 8, 9 },
};
HPRef_Struct reftrig_singedges3 =
{
  HP_TRIG,
  reftrig_singedges3_splitedges, 
  reftrig_singedges3_splitfaces, 
  0,
  reftrig_singedges3_newelstypes, 
  reftrig_singedges3_newels
};






// HP_TRIG_SINGEDGES23
int reftrig_singedges23_splitedges[][3] =
{
  { 1, 2, 4 },
  { 1, 3, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 0, 0, 0 }
};
int reftrig_singedges23_splitfaces[][4] =
{
  { 1, 2, 3, 10 },
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftrig_singedges23_newelstypes[] =
{
  //  HP_QUAD_2E,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,

  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG,
  HP_NONE,
};
int reftrig_singedges23_newels[][8] =
{
  //  { 1, 4, 10, 5 },
  { 1 , 4, 10 },
  { 5, 1, 10 },
  { 4, 6, 7, 10 },
  { 8, 5, 10, 9 },
  { 6, 2, 7 },
  { 3, 8, 9 },
  { 7, 9, 10 },
};
HPRef_Struct reftrig_singedges23 =
{
  HP_TRIG,
  reftrig_singedges23_splitedges, 
  reftrig_singedges23_splitfaces, 
  0,
  reftrig_singedges23_newelstypes, 
  reftrig_singedges23_newels
};


// HP_TRIG_3SINGEDGES
int reftrig_3singedges_splitedges[][3] =
{
  { 1, 2, 4 },
  { 2, 1, 5 }, 
  { 2, 3, 6 }, 
  { 3, 2, 7 }, 
  { 3, 1, 8 }, 
  { 1, 3, 9 }, 
  { 0, 0, 0 }
};
int reftrig_3singedges_splitfaces[][4] =
{
  { 1, 2, 3, 10 },
  { 2, 3, 1, 11 },
  { 3, 1, 2, 12 }, 
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftrig_3singedges_newelstypes[] =
{
  HP_TRIG, 
  HP_QUAD_SINGEDGE, 
  HP_QUAD_SINGEDGE, 
  HP_QUAD_SINGEDGE, 
  HP_TRIG_SINGEDGECORNER1, 
  HP_TRIG_SINGEDGECORNER2, 
  HP_TRIG_SINGEDGECORNER1, 
  HP_TRIG_SINGEDGECORNER2, 
  HP_TRIG_SINGEDGECORNER1, 
  HP_TRIG_SINGEDGECORNER2, 
  HP_NONE, 
};
int reftrig_3singedges_newels[][8] =
{
  { 10, 11, 12 }, 
  { 4, 5, 11, 10 }, 
  { 6, 7, 12, 11 },
  { 8, 9, 10, 12 }, 
  { 1, 4, 10 }, 
  { 9, 1, 10 }, 
  { 2, 6, 11 }, 
  { 5, 2, 11 }, 
  { 3, 8, 12 }, 
  { 7, 3, 12 }, 
};
HPRef_Struct reftrig_3singedges =
{
  HP_TRIG,
  reftrig_3singedges_splitedges, 
  reftrig_3singedges_splitfaces, 
  0,
  reftrig_3singedges_newelstypes, 
  reftrig_3singedges_newels
};
