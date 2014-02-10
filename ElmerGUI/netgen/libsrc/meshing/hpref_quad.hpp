// HP_QUAD
int refquad_splitedges[][3] =
{
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_newelstypes[] =
{
  HP_QUAD,
  HP_NONE,
};
int refquad_newels[][8] =
{
  { 1, 2, 3, 4 },
};
HPRef_Struct refquad =
{
  HP_QUAD,
  refquad_splitedges, 
  0, 0,
  refquad_newelstypes, 
  refquad_newels
};







// HP_QUAD_SINGCORNER
int refquad_singcorner_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_singcorner_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_TRIG,
  HP_NONE,
};
int refquad_singcorner_newels[][8] =
{
  { 1, 5, 6 },
  { 2, 4, 6, 5 },
  { 2, 3, 4 },
};
HPRef_Struct refquad_singcorner =
{
  HP_QUAD,
  refquad_singcorner_splitedges, 
  0, 0,
  refquad_singcorner_newelstypes, 
  refquad_singcorner_newels
};





// HP_DUMMY_QUAD_SINGCORNER
int refdummyquad_singcorner_splitedges[][3] =
{
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refdummyquad_singcorner_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_TRIG,
  HP_NONE,
};
int refdummyquad_singcorner_newels[][8] =
{
  { 1, 2, 4 },
  { 4, 2, 3 },
};
HPRef_Struct refdummyquad_singcorner =
{
  HP_QUAD,
  refdummyquad_singcorner_splitedges, 
  0, 0,
  refdummyquad_singcorner_newelstypes, 
  refdummyquad_singcorner_newels
};







// HP_QUAD_SINGEDGE
int refquad_singedge_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_singedge_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_singedge_newels[][8] =
{
  { 1, 2, 6, 5 },
  { 5, 6, 3, 4 },
};
HPRef_Struct refquad_singedge =
{
  HP_QUAD,
  refquad_singedge_splitedges, 
  0, 0,
  refquad_singedge_newelstypes, 
  refquad_singedge_newels
};






// HP_QUAD_0E_2VA
int refquad_0e_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_0e_2va_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int refquad_0e_2va_newels[][8] =
{
  { 1, 5, 6 },
  { 2, 8, 7 },
  { 5, 7, 8, 6 },
  { 6, 8, 3, 4 },
};
HPRef_Struct refquad_0e_2va =
{
  HP_QUAD,
  refquad_0e_2va_splitedges, 
  0, 0,
  refquad_0e_2va_newelstypes, 
  refquad_0e_2va_newels
};



// HP_QUAD_0E_2VB
int refquad_0e_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 3, 4, 7 },
  { 3, 2, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_0e_2vb_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int refquad_0e_2vb_newels[][8] =
{
  { 1, 5, 6 },
  { 3, 7, 8 },
  { 5, 2, 4, 6 },
  { 2, 8, 7, 4 },
};
HPRef_Struct refquad_0e_2vb =
{
  HP_QUAD,
  refquad_0e_2vb_splitedges, 
  0, 0,
  refquad_0e_2vb_newelstypes, 
  refquad_0e_2vb_newels
};




// HP_QUAD_0E_3V
int refquad_0e_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 0, 0, 0 }
};

int refquad_0e_3v_splitfaces[][4] =
{
  { 2, 3, 1, 14 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_0e_3v_newelstypes[] =
{
  HP_TRIG_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int refquad_0e_3v_newels[][8] =
{
  { 1, 5, 6 },
  { 2, 8, 14, 7 },
  { 3, 10, 9 },
  { 5, 7, 14, 6 },
  { 8, 9, 10, 14 },
  { 6, 14, 10, 4 },
};
HPRef_Struct refquad_0e_3v =
{
  HP_QUAD,
  refquad_0e_3v_splitedges, 
  refquad_0e_3v_splitfaces, 
  0,
  refquad_0e_3v_newelstypes, 
  refquad_0e_3v_newels
};




// HP_QUAD_0E_4V
int refquad_0e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 4, 1, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int refquad_0e_4v_splitfaces[][4] =
{
  { 1, 2, 4, 13 },
  { 2, 3, 1, 14 },
  { 3, 4, 2, 15 },
  { 4, 1, 3, 16 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_0e_4v_newelstypes[] =
{
  HP_DUMMY_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,

  HP_QUAD,
  HP_QUAD,
  HP_QUAD,
  HP_QUAD,

  HP_QUAD,
  HP_NONE,
};
int refquad_0e_4v_newels[][8] =
{
  { 1, 5, 13, 6 },
  { 2, 8, 14, 7 },
  { 3, 10, 15, 9 },
  { 4, 11, 16, 12 },
  { 5, 7, 14, 13 },
  { 8, 9, 15, 14 },
  { 10, 12, 16, 15 },
  { 11, 6, 13, 16 },
  { 13, 14, 15, 16 }
};
HPRef_Struct refquad_0e_4v =
{
  HP_QUAD,
  refquad_0e_4v_splitedges, 
  refquad_0e_4v_splitfaces, 
  0,
  refquad_0e_4v_newelstypes, 
  refquad_0e_4v_newels
};








// HP_QUAD_1E_1VA
int refquad_1e_1va_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_1va_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGEDGECORNER1,
  HP_NONE,
};
int refquad_1e_1va_newels[][8] =
{
  { 7, 2, 6, 5 },
  { 5, 6, 3, 4 },
  { 1, 7, 5 },
};
HPRef_Struct refquad_1e_1va =
{
  HP_QUAD,
  refquad_1e_1va_splitedges, 
  0, 0,
  refquad_1e_1va_newelstypes, 
  refquad_1e_1va_newels
};




// HP_QUAD_1E_1VB
int refquad_1e_1vb_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 2, 1, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_1vb_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGEDGECORNER2,
  HP_NONE,
};
int refquad_1e_1vb_newels[][8] =
{
  { 1, 7, 6, 5 },
  { 5, 6, 3, 4 },
  { 7, 2, 6 },
};
HPRef_Struct refquad_1e_1vb =
{
  HP_QUAD,
  refquad_1e_1vb_splitedges, 
  0, 0,
  refquad_1e_1vb_newelstypes, 
  refquad_1e_1vb_newels
};



// HP_QUAD_1E_1VC
int refquad_1e_1vc_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 3, 4, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_1vc_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_1vc_newels[][8] =
{
  { 1, 2, 6, 5 },
  { 5, 6, 4 },
  { 4, 6, 7, 8 },
  { 3, 8, 7 }
};
HPRef_Struct refquad_1e_1vc =
{
  HP_QUAD,
  refquad_1e_1vc_splitedges, 
  0, 0,
  refquad_1e_1vc_newelstypes, 
  refquad_1e_1vc_newels
};



// HP_QUAD_1E_1VD
int refquad_1e_1vd_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 4, 1, 7 },
  { 4, 3, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_1vd_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_1vd_newels[][8] =
{
  { 1, 2, 6, 5 },
  { 5, 6, 3 },
  { 5, 3, 8, 7 },
  { 4, 7, 8 }
};
HPRef_Struct refquad_1e_1vd =
{
  HP_QUAD,
  refquad_1e_1vd_splitedges, 
  0, 0,
  refquad_1e_1vd_newelstypes, 
  refquad_1e_1vd_newels
};







// HP_QUAD_1E_2VA
int refquad_1e_2va_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 2, 1, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_2va_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_NONE,
};
int refquad_1e_2va_newels[][8] =
{
  { 7, 8, 6, 5 },
  { 5, 6, 3, 4 },
  { 1, 7, 5 },
  { 8, 2, 6 }
};
HPRef_Struct refquad_1e_2va =
{
  HP_QUAD,
  refquad_1e_2va_splitedges, 
  0, 0,
  refquad_1e_2va_newelstypes, 
  refquad_1e_2va_newels
};




// HP_QUAD_1E_2VB
int refquad_1e_2vb_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 3, 2, 8 },
  { 3, 4, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_2vb_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_2vb_newels[][8] =
{
  { 7, 2, 6, 5 },
  { 1, 7, 5 },
  { 5, 6, 4 },
  { 4, 6, 8, 9 },
  { 3, 9, 8 }
};
HPRef_Struct refquad_1e_2vb =
{
  HP_QUAD,
  refquad_1e_2vb_splitedges, 
  0, 0,
  refquad_1e_2vb_newelstypes, 
  refquad_1e_2vb_newels
};




// HP_QUAD_1E_2VC
int refquad_1e_2vc_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 4, 1, 8 },
  { 4, 3, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_2vc_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_2vc_newels[][8] =
{
  { 7, 2, 6, 5 },
  { 1, 7, 5 },
  { 5, 6, 3 },
  { 5, 3, 9, 8 },
  { 4, 8, 9 }
};
HPRef_Struct refquad_1e_2vc =
{
  HP_QUAD,
  refquad_1e_2vc_splitedges, 
  0, 0,
  refquad_1e_2vc_newelstypes, 
  refquad_1e_2vc_newels
};




// HP_QUAD_1E_2VD
int refquad_1e_2vd_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 2, 1, 7 },
  { 3, 2, 8 },
  { 3, 4, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_2vd_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_2vd_newels[][8] =
{
  { 1, 7, 6, 5 },
  { 7, 2, 6 },
  { 5, 6, 4 },
  { 4, 6, 8, 9 },
  { 3, 9, 8 }
};
HPRef_Struct refquad_1e_2vd =
{
  HP_QUAD,
  refquad_1e_2vd_splitedges, 
  0, 0,
  refquad_1e_2vd_newelstypes, 
  refquad_1e_2vd_newels
};





// HP_QUAD_1E_2VE
int refquad_1e_2ve_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 2, 1, 7 },
  { 4, 1, 8 },
  { 4, 3, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_2ve_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_2ve_newels[][8] =
{
  { 1, 7, 6, 5 },
  { 7, 2, 6 },
  { 5, 6, 3 },
  { 5, 3, 9, 8 },
  { 4, 8, 9 }
};
HPRef_Struct refquad_1e_2ve =
{
  HP_QUAD,
  refquad_1e_2ve_splitedges, 
  0, 0,
  refquad_1e_2ve_newelstypes, 
  refquad_1e_2ve_newels
};






// HP_QUAD_1E_2VF
int refquad_1e_2vf_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 4, 1, 7 },
  { 4, 3, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_2vf_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_1e_2vf_newels[][8] =
{
  { 1, 2, 6, 5 },
  { 5, 6, 9, 7 },
  { 7, 9, 10, 8 },
  { 4, 7, 8 },
  { 3, 10, 9 },
};
HPRef_Struct refquad_1e_2vf =
{
  HP_QUAD,
  refquad_1e_2vf_splitedges, 
  0, 0,
  refquad_1e_2vf_newelstypes, 
  refquad_1e_2vf_newels
};





// HP_QUAD_1E_3VA
int refquad_1e_3va_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 2, 1, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_3va_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG,
  HP_NONE,
};
int refquad_1e_3va_newels[][8] =
{
  { 1, 7, 5 },
  { 8, 2, 6 },
  { 3, 10, 9 },
  { 7, 8, 6, 5 },
  { 4, 6, 9, 10 },
  { 5, 6, 4 }
};
HPRef_Struct refquad_1e_3va =
{
  HP_QUAD,
  refquad_1e_3va_splitedges, 
  0, 0,
  refquad_1e_3va_newelstypes, 
  refquad_1e_3va_newels
};





// HP_QUAD_1E_3VB
int refquad_1e_3vb_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 2, 1, 8 },
  { 4, 1, 9 },
  { 4, 3, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_3vb_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG,
  HP_NONE,
};
int refquad_1e_3vb_newels[][8] =
{
  { 1, 7, 5 },
  { 8, 2, 6 },
  { 4, 9, 10 },
  { 7, 8, 6, 5 },
  { 5, 3, 10, 9 },
  { 5, 6, 3 }
};
HPRef_Struct refquad_1e_3vb =
{
  HP_QUAD,
  refquad_1e_3vb_splitedges, 
  0, 0,
  refquad_1e_3vb_newelstypes, 
  refquad_1e_3vb_newels
};





// HP_QUAD_1E_3VC
int refquad_1e_3vc_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 3, 2, 8 },
  { 3, 4, 9 },
  { 4, 3, 10 },
  { 4, 1, 11 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_3vc_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int refquad_1e_3vc_newels[][8] =
{
  { 1, 7, 5 },
  { 3, 9, 8 },
  { 4, 11, 10 },
  { 7, 2, 6, 5 },
  { 5, 6, 8, 11 },
  { 11, 8, 9, 10 }
};
HPRef_Struct refquad_1e_3vc =
{
  HP_QUAD,
  refquad_1e_3vc_splitedges, 
  0, 0,
  refquad_1e_3vc_newelstypes, 
  refquad_1e_3vc_newels
};




// HP_QUAD_1E_3VD
int refquad_1e_3vd_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 2, 1, 7 },
  { 3, 2, 8 },
  { 3, 4, 9 },
  { 4, 3, 10 },
  { 4, 1, 11 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_3vd_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int refquad_1e_3vd_newels[][8] =
{
  { 7, 2, 6 },
  { 3, 9, 8 },
  { 4, 11, 10 },
  { 1, 7, 6, 5 },
  { 5, 6, 8, 11 },
  { 11, 8, 9, 10 }
};
HPRef_Struct refquad_1e_3vd =
{
  HP_QUAD,
  refquad_1e_3vd_splitedges, 
  0, 0,
  refquad_1e_3vd_newelstypes, 
  refquad_1e_3vd_newels
};






// HP_QUAD_1E_4V
int refquad_1e_4v_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 2, 1, 8 },
  { 4, 1, 9 },
  { 3, 2, 10 },
  { 4, 3, 11 },
  { 3, 4, 12 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_1e_4v_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_QUAD,
  HP_NONE,
};
int refquad_1e_4v_newels[][8] =
{
  { 1, 7, 5 },
  { 8, 2, 6 },
  { 3, 12, 10 },
  { 4, 9, 11 },
  { 7, 8, 6, 5 },
  { 5, 6, 10, 9 },
  { 9, 10, 12, 11 }
};
HPRef_Struct refquad_1e_4v =
{
  HP_QUAD,
  refquad_1e_4v_splitedges, 
  0, 0,
  refquad_1e_4v_newelstypes, 
  refquad_1e_4v_newels
};

////////////////////////////////////////////////////////////////////////////////

// HP_QUAD_2E
int refquad_2e_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 3, 8 },
  { 0, 0, 0 }
};
int refquad_2e_splitfaces[][4] =
{
  { 1, 2, 4, 9 },
  { 0, 0, 0, 0 },
};


/* 
   HPREF_ELEMENT_TYPE refquad_2e_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2e_newels[][8] =
{
  { 1, 5, 9 },
  { 6, 1, 9 },
  { 5, 2, 7, 9 },
  { 4, 6, 9, 8 },
  { 9, 7, 3, 8 },
};
*/ 

// SZ refine to 4 quads 
HPREF_ELEMENT_TYPE refquad_2e_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2e_newels[][8] =
{
  { 1, 5, 9, 6 },
  { 5, 2, 7, 9 },
  { 4, 6, 9, 8 },
  { 9, 7, 3, 8 },
};

HPRef_Struct refquad_2e =
{
  HP_QUAD,
  refquad_2e_splitedges, 
  refquad_2e_splitfaces, 
  0,
  refquad_2e_newelstypes, 
  refquad_2e_newels
};


// HP_QUAD_2E_1VA
int refquad_2e_1va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 3, 8 },
  { 2, 1, 10 },
  { 0, 0, 0 }
};
int refquad_2e_1va_splitfaces[][4] =
{
  { 1, 2, 4, 9 },
  { 0, 0, 0, 0 },
};

/* 
HPREF_ELEMENT_TYPE refquad_2e_1va_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGEDGECORNER2,
  HP_NONE,
};
int refquad_2e_1va_newels[][8] =
{
  { 1, 5, 9 },
  { 6, 1, 9 },
  { 5, 10, 7, 9 },
  { 4, 6, 9, 8 },
  { 9, 7, 3, 8 },
  { 10, 2, 7 },
};
*/ 
// SZ Quad_2e refinement 
HPREF_ELEMENT_TYPE refquad_2e_1va_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGEDGECORNER2,
  HP_NONE,
};
int refquad_2e_1va_newels[][8] =
{
  { 1, 5, 9, 6 },
  { 5, 10, 7, 9 },
  { 4, 6, 9, 8 },
  { 9, 7, 3, 8 },
  { 10, 2, 7 },
};

HPRef_Struct refquad_2e_1va =
{
  HP_QUAD,
  refquad_2e_1va_splitedges, 
  refquad_2e_1va_splitfaces, 
  0,
  refquad_2e_1va_newelstypes, 
  refquad_2e_1va_newels
};



// HP_QUAD_2E_1VB
int refquad_2e_1vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 3, 8 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 0, 0, 0 }
};
int refquad_2e_1vb_splitfaces[][4] =
{
  { 1, 2, 4, 9 },
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2e_1vb_newelstypes[] =
{
  // HP_TRIG_SINGEDGECORNER1,
  // HP_TRIG_SINGEDGECORNER2,
  // SZ QUAD_2E 
  HP_QUAD_2E,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_NONE,
};
int refquad_2e_1vb_newels[][8] =
{
  //{ 1, 5, 9 },
  //{ 6, 1, 9 },
  { 1, 5, 9, 6 },
  { 5, 2, 7, 9 },
  { 4, 6, 9, 8 },
  { 7, 8, 9 },
  { 8, 7, 10, 11 },
  { 3, 11, 10 }
};
HPRef_Struct refquad_2e_1vb =
{
  HP_QUAD,
  refquad_2e_1vb_splitedges, 
  refquad_2e_1vb_splitfaces, 
  0,
  refquad_2e_1vb_newelstypes, 
  refquad_2e_1vb_newels
}
;

// HP_QUAD_2E_1VC
int refquad_2e_1vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 1, 8 },
  { 4, 3, 9 },
  { 0, 0, 0 }
};
int refquad_2e_1vc_splitfaces[][4] =
{
  { 1, 2, 4, 10 },
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2e_1vc_newelstypes[] =
{
  //  HP_TRIG_SINGEDGECORNER1,
  // HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_2E, 
  HP_TRIG_SINGEDGECORNER1,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2e_1vc_newels[][8] =
{
  //{ 1, 5, 10 },
  //{ 6, 1, 10 },
  { 1, 5, 10, 6}, 
  { 4, 8, 9 },
  { 5, 2, 7, 10 },
  { 8, 6, 10, 9 },
  { 10, 7, 3, 9 },
};
HPRef_Struct refquad_2e_1vc =
{
  HP_QUAD,
  refquad_2e_1vc_splitedges, 
  refquad_2e_1vc_splitfaces, 
  0,
  refquad_2e_1vc_newelstypes, 
  refquad_2e_1vc_newels
};

// HP_QUAD_2E_2VA
int refquad_2e_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 3, 8 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 2, 1, 12 },
  { 0, 0, 0 }
};
int refquad_2e_2va_splitfaces[][4] =
{
  { 1, 2, 4, 9 },
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2e_2va_newelstypes[] =
{
  //HP_TRIG_SINGEDGECORNER1,
  //HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_2E,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGEDGECORNER2,
  HP_NONE,
};
int refquad_2e_2va_newels[][8] =
{
  // { 1, 5, 9 },
  // { 6, 1, 9 },
  { 1, 5, 9, 6 }, 
  { 5, 12, 7, 9 },
  { 4, 6, 9, 8 },
  { 7, 8, 9 },
  { 8, 7, 10, 11 },
  { 3, 11, 10 },
  { 12, 2, 7 }
};
HPRef_Struct refquad_2e_2va =
{
  HP_QUAD,
  refquad_2e_2va_splitedges, 
  refquad_2e_2va_splitfaces, 
  0,
  refquad_2e_2va_newelstypes, 
  refquad_2e_2va_newels
};






// HP_QUAD_2E_2VB
int refquad_2e_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 4, 1, 9 },
  { 4, 3, 10 },
  { 0, 0, 0 }
};
int refquad_2e_2vb_splitfaces[][4] =
{
  { 1, 2, 4, 11 },
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2e_2vb_newelstypes[] =
{
  // HP_TRIG_SINGEDGECORNER1,
  // HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_2E, 
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2e_2vb_newels[][8] =
{
  //{ 1, 5, 11 },
  //{ 6, 1, 11 },
  { 1, 5, 11, 6 }, 
  { 4, 9, 10 },
  { 7, 2, 8 },
  { 5, 7, 8, 11 },
  { 9, 6, 11, 10 },
  { 3, 10, 11, 8 },
};
HPRef_Struct refquad_2e_2vb =
{
  HP_QUAD,
  refquad_2e_2vb_splitedges, 
  refquad_2e_2vb_splitfaces, 
  0,
  refquad_2e_2vb_newelstypes, 
  refquad_2e_2vb_newels
};

// HP_QUAD_2E_2VC
int refquad_2e_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 3, 8 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 1, 12 },
  { 0, 0, 0 }
};
int refquad_2e_2vc_splitfaces[][4] =
{
  { 1, 2, 4, 9 },
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2e_2vc_newelstypes[] =
{
  // HP_TRIG_SINGEDGECORNER1,
  // HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_2E, 
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGEDGECORNER1, //SZ (vorher: SINGEDGECORNER2) 
  HP_NONE,
};
int refquad_2e_2vc_newels[][8] =
{
  { 1, 5, 9 },
  { 6, 1, 9 },
  { 5, 2, 7, 9 },
  { 12, 6, 9, 8 },
  { 7, 8, 9 },
  { 8, 7, 10, 11 },
  { 3, 11, 10 },
  { 4, 12, 8 }
};
HPRef_Struct refquad_2e_2vc =
{
  HP_QUAD,
  refquad_2e_2vc_splitedges, 
  refquad_2e_2vc_splitfaces, 
  0,
  refquad_2e_2vc_newelstypes, 
  refquad_2e_2vc_newels
};

// HP_QUAD_2E_3V  
int refquad_2e_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 4, 3, 8 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 2, 1, 12 },
  { 4, 1, 13 },
  { 0, 0, 0 }
};
int refquad_2e_3v_splitfaces[][4] =
{
  { 1, 2, 4, 9 },
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2e_3v_newelstypes[] =
{
  // HP_TRIG_SINGEDGECORNER1,
  // HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_2E, 
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG,
  HP_QUAD,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER1,
  HP_NONE,
};
int refquad_2e_3v_newels[][8] =
{
  //{ 1, 5, 9 },
  //{ 6, 1, 9 },
  { 1, 5, 9, 6 }, 
  { 5, 12, 7, 9 },
  { 13, 6, 9, 8 },
  { 7, 8, 9 },
  { 8, 7, 10, 11 },
  { 3, 11, 10 },
  { 12, 2, 7 },
  { 4, 13, 8 }
};
HPRef_Struct refquad_2e_3v =
{
  HP_QUAD,
  refquad_2e_3v_splitedges, 
  refquad_2e_3v_splitfaces, 
  0,
  refquad_2e_3v_newelstypes, 
  refquad_2e_3v_newels
};

// HP_QUAD_2EB_0V
int refquad_2eb_0v_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 0, 0, 0 }
};
int refquad_2eb_0v_splitfaces[][4] =
{
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2eb_0v_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_0v_newels[][8] =
{
  { 1, 2, 6, 5 },
  { 3, 4, 8, 7 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_0v =
{
  HP_QUAD,
  refquad_2eb_0v_splitedges, 
  refquad_2eb_0v_splitfaces, 
  0,
  refquad_2eb_0v_newelstypes, 
  refquad_2eb_0v_newels
};


// HP_QUAD_2EB_1VA
int refquad_2eb_1va_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 1, 2, 9 },
  { 0, 0, 0 }
};
int refquad_2eb_1va_splitfaces[][4] =
{
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2eb_1va_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_1va_newels[][8] =
{
  { 9, 2, 6, 5 },
  { 3, 4, 8, 7 },
  { 1, 9, 5 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_1va =
{
  HP_QUAD,
  refquad_2eb_1va_splitedges, 
  refquad_2eb_1va_splitfaces, 
  0,
  refquad_2eb_1va_newelstypes, 
  refquad_2eb_1va_newels
};

// HP_QUAD_2EB_1VB
int refquad_2eb_1vb_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 2, 1, 9 },
  { 0, 0, 0 }
};
int refquad_2eb_1vb_splitfaces[][4] =
{
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2eb_1vb_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_1vb_newels[][8] =
{
  { 1, 9, 6, 5 },
  { 3, 4, 8, 7 },
  { 9, 2, 6 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_1vb =
{
  HP_QUAD,
  refquad_2eb_1vb_splitedges, 
  refquad_2eb_1vb_splitfaces, 
  0,
  refquad_2eb_1vb_newelstypes, 
  refquad_2eb_1vb_newels
};

// HP_QUAD_2EB_2VA
int refquad_2eb_2va_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 1, 2, 9 },
  { 2, 1, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_2eb_2va_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_2va_newels[][8] =
{
  { 9, 10, 6, 5 },
  { 3, 4, 8, 7 },
  { 1, 9, 5 },
  { 10, 2, 6 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_2va =
{
  HP_QUAD,
  refquad_2eb_2va_splitedges, 
  0, 0,
  refquad_2eb_2va_newelstypes, 
  refquad_2eb_2va_newels
};



// HP_QUAD_2EB_2VB
int refquad_2eb_2vb_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 1, 2, 9 },
  { 3, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_2eb_2vb_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER1,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_2vb_newels[][8] =
{
  { 9, 2, 6, 5 },
  { 10, 4, 8, 7 },
  { 1, 9, 5 },
  { 3, 10, 7 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_2vb =
{
  HP_QUAD,
  refquad_2eb_2vb_splitedges, 
  0, 0,
  refquad_2eb_2vb_newelstypes, 
  refquad_2eb_2vb_newels
};



// HP_QUAD_2EB_2VC
int refquad_2eb_2vc_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 1, 2, 9 },
  { 4, 3, 10 },
  { 0, 0, 0 }
};
int refquad_2eb_2vc_splitfaces[][4] =
{
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2eb_2vc_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_2vc_newels[][8] =
{
  { 9, 2, 6, 5 },
  { 3, 10, 8, 7 },
  { 1, 9, 5 },
  { 10, 4, 8 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_2vc =
{
  HP_QUAD,
  refquad_2eb_2vc_splitedges, 
  refquad_2eb_2vc_splitfaces, 
  0,
  refquad_2eb_2vc_newelstypes, 
  refquad_2eb_2vc_newels
};


// HP_QUAD_2EB_2VD
int refquad_2eb_2vd_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 2, 1, 9 },
  { 4, 3, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_2eb_2vd_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_2vd_newels[][8] =
{
  { 1, 9, 6, 5 },
  { 3, 10, 8, 7 },
  { 9, 2, 6 },
  { 10, 4, 8 },
  { 5, 6, 7, 8 }
};
HPRef_Struct refquad_2eb_2vd =
{
  HP_QUAD,
  refquad_2eb_2vd_splitedges, 
  0, 0,
  refquad_2eb_2vd_newelstypes, 
  refquad_2eb_2vd_newels
};


// HP_QUAD_2EB_3VA
int refquad_2eb_3va_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 2, 1, 8 },
  { 3, 2, 9 },
  { 4, 1, 10 },
  { 3, 4, 11 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_2eb_3va_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER1,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_3va_newels[][8] =
{
  { 1, 7, 5 },
  { 8, 2, 6 },
  { 3, 11, 9},
  { 7, 8, 6, 5 },
  { 11, 4, 10, 9 },
  { 5, 6, 9, 10 }
};
HPRef_Struct refquad_2eb_3va =
{
  HP_QUAD,
  refquad_2eb_3va_splitedges, 
  0, 0,
  refquad_2eb_3va_newelstypes, 
  refquad_2eb_3va_newels
};


// HP_QUAD_2EB_3VB
int refquad_2eb_3vb_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 1, 2, 7 },
  { 2, 1, 8 },
  { 3, 2, 9 },
  { 4, 1, 10 },
  { 4, 3, 11 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refquad_2eb_3vb_newelstypes[] =
{
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER2,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_NONE,
};
int refquad_2eb_3vb_newels[][8] =
{
  { 1, 7, 5 },
  { 8, 2, 6 },
  { 11, 4, 10 },
  { 7, 8, 6, 5 },
  { 3, 11, 10, 9 },
  { 5, 6, 9, 10 }
};
HPRef_Struct refquad_2eb_3vb =
{
  HP_QUAD,
  refquad_2eb_3vb_splitedges, 
  0, 0,
  refquad_2eb_3vb_newelstypes, 
  refquad_2eb_3vb_newels
};


// HP_QUAD_2EB_4V
int refquad_2eb_4v_splitedges[][3] =
{
  { 1, 4, 5 },
  { 2, 3, 6 },
  { 3, 2, 7 },
  { 4, 1, 8 },
  { 1, 2, 9 },
  { 2, 1, 10 },
  { 3, 4, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};
int refquad_2eb_4v_splitfaces[][4] =
{
  { 0, 0, 0, 0 },
};
HPREF_ELEMENT_TYPE refquad_2eb_4v_newelstypes[] =
{
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_NONE,
};
int refquad_2eb_4v_newels[][8] =
{
  { 9, 10, 6, 5 },
  { 11, 12, 8, 7 },
  { 5, 6, 7, 8 },
  { 1, 9, 5 },
  { 10, 2, 6 },
  { 3, 11, 7 },
  { 12, 4, 8 },
};
HPRef_Struct refquad_2eb_4v =
{
  HP_QUAD,
  refquad_2eb_4v_splitedges, 
  refquad_2eb_4v_splitfaces, 
  0,
  refquad_2eb_4v_newelstypes, 
  refquad_2eb_4v_newels
};



// HP_QUAD_3E
int refquad_3e_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 4, 10 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int refquad_3e_splitfaces[][4] =
{
  { 1, 2, 4, 13 },
  { 2, 3, 1, 14 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_3e_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_2E,
//   HP_TRIG_SINGEDGECORNER1,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER1,

  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,

  HP_QUAD,
  HP_NONE,
};
int refquad_3e_newels[][8] =
{
//   { 1, 5, 13 },
//   { 6, 1, 13 },
//   { 7, 2, 14 },
//   { 2, 8, 14 },
  { 1, 5, 13, 6 },
  { 2, 8, 14, 7 },
  { 5, 7, 14, 13 },
  { 8, 3, 10, 14 },
  { 4, 6, 13, 12 },
  { 13, 14, 10, 12 }
};
HPRef_Struct refquad_3e =
{
  HP_QUAD,
  refquad_3e_splitedges, 
  refquad_3e_splitfaces, 
  0,
  refquad_3e_newelstypes, 
  refquad_3e_newels
};







// HP_QUAD_3E_3VA
int refquad_3e_3va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 4, 10 },
  { 3, 2, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int refquad_3e_3va_splitfaces[][4] =
{
  { 1, 2, 4, 13 },
  { 2, 3, 1, 14 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_3e_3va_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_2E,

//   HP_TRIG_SINGEDGECORNER1,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,

  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,

  HP_QUAD,
  HP_NONE,
};
int refquad_3e_3va_newels[][8] =
{
//   { 1, 5, 13 },
//   { 6, 1, 13 },
//   { 7, 2, 14 },
//   { 2, 8, 14 },
  { 1, 5, 13, 6 },
  { 2, 8, 14, 7 },
  { 11, 3, 10 },
  { 5, 7, 14, 13 },
  { 8, 11, 10, 14 },
  { 4, 6, 13, 12 },
  { 13, 14, 10, 12 }
};
HPRef_Struct refquad_3e_3va =
{
  HP_QUAD,
  refquad_3e_3va_splitedges, 
  refquad_3e_3va_splitfaces, 
  0,
  refquad_3e_3va_newelstypes, 
  refquad_3e_3va_newels
};

// HP_QUAD_3E_3VB
int refquad_3e_3vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 4, 10 },
  { 4, 1, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int refquad_3e_3vb_splitfaces[][4] =
{
  { 1, 2, 4, 13 },
  { 2, 3, 1, 14 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_3e_3vb_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_2E,

//   HP_TRIG_SINGEDGECORNER1,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER1,

  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,

  HP_QUAD,
  HP_NONE,
};
int refquad_3e_3vb_newels[][8] =
{
//   { 1, 5, 13 },
//   { 6, 1, 13 },
//   { 7, 2, 14 },
//   { 2, 8, 14 },
  { 1, 5, 13, 6 },
  { 2, 8, 14, 7 },
  { 4, 11, 12 },
  { 5, 7, 14, 13 },
  { 8, 3, 10, 14 },
  { 11, 6, 13, 12 },
  { 13, 14, 10, 12 }
};
HPRef_Struct refquad_3e_3vb =
{
  HP_QUAD,
  refquad_3e_3vb_splitedges, 
  refquad_3e_3vb_splitfaces, 
  0,
  refquad_3e_3vb_newelstypes, 
  refquad_3e_3vb_newels
};









// HP_QUAD_3E_4V
int refquad_3e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 4, 10 },
  { 3, 2, 11 },
  { 4, 3, 12 },
  { 4, 1, 15 },
  { 0, 0, 0 }
};

int refquad_3e_4v_splitfaces[][4] =
{
  { 1, 2, 4, 13 },
  { 2, 3, 1, 14 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_3e_4v_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_2E,

//   HP_TRIG_SINGEDGECORNER1,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER2,
//   HP_TRIG_SINGEDGECORNER1,
  HP_TRIG_SINGEDGECORNER2,
  HP_TRIG_SINGEDGECORNER1,

  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,

  HP_QUAD,
  HP_NONE,
};
int refquad_3e_4v_newels[][8] =
{
//   { 1, 5, 13 },
//   { 6, 1, 13 },
//   { 7, 2, 14 },
//   { 2, 8, 14 },
  { 1, 5, 13, 6 },
  { 2, 8, 14, 7 },
  { 11, 3, 10 },
  { 4, 15, 12 },
  { 5, 7, 14, 13 },
  { 8, 11, 10, 14 },
  { 15, 6, 13, 12 },
  { 13, 14, 10, 12 }
};
HPRef_Struct refquad_3e_4v =
{
  HP_QUAD,
  refquad_3e_4v_splitedges, 
  refquad_3e_4v_splitfaces, 
  0,
  refquad_3e_4v_newelstypes, 
  refquad_3e_4v_newels
};









// HP_QUAD_4E
int refquad_4e_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 4, 1, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int refquad_4e_splitfaces[][4] =
{
  { 1, 2, 4, 13 },
  { 2, 3, 1, 14 },
  { 3, 4, 2, 15 },
  { 4, 1, 3, 16 },
  { 0, 0, 0, 0 },
};

HPREF_ELEMENT_TYPE refquad_4e_newelstypes[] =
{
  HP_QUAD_2E,
  HP_QUAD_2E,
  HP_QUAD_2E,
  HP_QUAD_2E,

  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,
  HP_QUAD_SINGEDGE,

  HP_QUAD,
  HP_NONE,
};
int refquad_4e_newels[][8] =
{
  { 1, 5, 13, 6 },
  { 2, 8, 14, 7 },
  { 3, 10, 15, 9 },
  { 4, 11, 16, 12 },
  { 5, 7, 14, 13 },
  { 8, 9, 15, 14 },
  { 10, 12, 16, 15 },
  { 11, 6, 13, 16 },
  { 13, 14, 15, 16 }
};
HPRef_Struct refquad_4e =
{
  HP_QUAD,
  refquad_4e_splitedges, 
  refquad_4e_splitfaces, 
  0,
  refquad_4e_newelstypes, 
  refquad_4e_newels
};
