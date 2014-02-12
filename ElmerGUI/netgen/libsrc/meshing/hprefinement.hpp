#ifndef FILE_HPREFINEMENT
#define FILE_HPREFINEMENT

/**************************************************************************/
/* File:   hprefinement.hh                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Oct. 2000                                                  */
/**************************************************************************/

/*
  HP Refinement
*/




enum HPREF_ELEMENT_TYPE {
  HP_NONE=0,

  HP_SEGM = 1,
  HP_SEGM_SINGCORNERL,
  HP_SEGM_SINGCORNERR,
  HP_SEGM_SINGCORNERS,

  HP_TRIG = 10,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER12,
  HP_TRIG_SINGCORNER123,
  HP_TRIG_SINGCORNER123_2D,   // not rotational symmetric
  HP_TRIG_SINGEDGE = 20,
  HP_TRIG_SINGEDGECORNER1,   // E = 100, V = 100
  HP_TRIG_SINGEDGECORNER2,   // E = 100, V = 010
  HP_TRIG_SINGEDGECORNER12,  // E = 100, V = 110
  HP_TRIG_SINGEDGECORNER3,
  HP_TRIG_SINGEDGECORNER13,
  HP_TRIG_SINGEDGECORNER23,
  HP_TRIG_SINGEDGECORNER123,
  HP_TRIG_SINGEDGES = 30,
  HP_TRIG_SINGEDGES2,
  HP_TRIG_SINGEDGES3,
  HP_TRIG_SINGEDGES23,
  HP_TRIG_3SINGEDGES = 40,

  HP_QUAD = 50,
  HP_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD_0E_2VA,  // V = 1100
  HP_QUAD_0E_2VB,  // V = 1010
  HP_QUAD_0E_3V,
  HP_QUAD_0E_4V,

  // one edge: marked edge is always edge from vertex 1 to vertex 2 (E = 1000)
  HP_QUAD_1E_1VA,  // vertex on beginning of edge: V = 1000
  HP_QUAD_1E_1VB,  // vertex on end of edge: V = 0100
  HP_QUAD_1E_1VC,  // V = 0010
  HP_QUAD_1E_1VD,  // V = 0001

  HP_QUAD_1E_2VA,  // V = 1100
  HP_QUAD_1E_2VB,  // V = 1010
  HP_QUAD_1E_2VC,  // V = 1001
  HP_QUAD_1E_2VD,  // V = 0110
  HP_QUAD_1E_2VE,  // V = 0101
  HP_QUAD_1E_2VF,  // V = 0011

  HP_QUAD_1E_3VA,  // V = 1110
  HP_QUAD_1E_3VB,  // V = 1101
  HP_QUAD_1E_3VC,  // V = 1011
  HP_QUAD_1E_3VD,  // V = 0111

  HP_QUAD_1E_4V,   // V = 1111


  HP_QUAD_2E,      // E = 1001, V = 1000
  HP_QUAD_2E_1VA,  // E = 1001, V = 1100
  HP_QUAD_2E_1VB,  // E = 1001, V = 1010
  HP_QUAD_2E_1VC,  // E = 1001, V = 1001
  HP_QUAD_2E_2VA,  // E = 1001, V = 1110
  HP_QUAD_2E_2VB,  // E = 1001, V = 1101
  HP_QUAD_2E_2VC,  // E = 1001, V = 1011
  HP_QUAD_2E_3V,   // E = 1001, V = 1111

  HP_QUAD_2EB_0V,   // E = 1010, V = 0000
  HP_QUAD_2EB_1VA,  // E = 1010, V = 1000
  HP_QUAD_2EB_1VB,  // E = 1010, V = 0100
  HP_QUAD_2EB_2VA,  // E = 1010, V = 1100
  HP_QUAD_2EB_2VB,  // E = 1010, V = 1010
  HP_QUAD_2EB_2VC,  // E = 1010, V = 1001
  HP_QUAD_2EB_2VD,  // E = 1010, V = 0101
  HP_QUAD_2EB_3VA,  // E = 1010, V = 1110
  HP_QUAD_2EB_3VB,  // E = 1010, V = 1101

  HP_QUAD_2EB_4V,


  HP_QUAD_3E,      // E = 1101, V = 1100
  HP_QUAD_3E_3VA,  // E = 1101, V = 1110
  HP_QUAD_3E_3VB,  // E = 1101, V = 1101
  HP_QUAD_3E_4V,   // E = 1101, V = 1111

  HP_QUAD_4E,


  HP_TET = 100,     // no singular vertex/edge
  HP_TET_0E_1V,     // V1
  HP_TET_0E_2V,     // V1,2
  HP_TET_0E_3V,     // V1,2,3  
  HP_TET_0E_4V,     // V1,2,3,4
  HP_TET_1E_0V = 200,   // E1-2
  HP_TET_1E_1VA,    // V1
  HP_TET_1E_1VB,    // V3
  HP_TET_1E_2VA,    // V1,2
  HP_TET_1E_2VB,    // V1,3
  HP_TET_1E_2VC,    // V1,4
  HP_TET_1E_2VD,    // V3,4
  HP_TET_1E_3VA,    // V1,2,3
  HP_TET_1E_3VB,    // V1,3,4
  HP_TET_1E_4V,     // V1,2,3,4


  // 2 connected edges, additonally marked Vs
  HP_TET_2EA_0V = 220,    // E1-2, E1-3
  HP_TET_2EA_1VA,   // V2
  HP_TET_2EA_1VB,   // V3
  HP_TET_2EA_1VC,   // V4
  HP_TET_2EA_2VA,   // V2,3
  HP_TET_2EA_2VB,   // V2,4
  HP_TET_2EA_2VC,   // V3,4
  HP_TET_2EA_3V,    // V2,3,4

  // 2 opposite edges
  HP_TET_2EB_0V = 230,    // E1-2, E3-4
  HP_TET_2EB_1V,    // V1
  HP_TET_2EB_2VA,   // V1,2
  HP_TET_2EB_2VB,   // V1,3
  HP_TET_2EB_2VC,   // V1,4
  HP_TET_2EB_3V,    // V1,2,3
  HP_TET_2EB_4V,    // V1,2,3,4

  HP_TET_3EA_0V = 400,  // E1-2, E1-3, E1-4, 3 edges connected
  HP_TET_3EA_1V,        // V2
  HP_TET_3EA_2V,        // V2,3
  HP_TET_3EA_3V,        // V2,3,4

  HP_TET_3EB_0V = 420,  // E1-2, E1-4, E2-3  3 edges chain
  HP_TET_3EB_1V,        // 
  HP_TET_3EB_2V,        // 
  HP_TET_3EC_0V = 430,  // 3 edges chain, alter
  HP_TET_3EC_1V,        // 3 edges chain, alter
  HP_TET_3EC_2V,        // 3 edges chain, alter


  HP_TET_1F_0E_0V = 500,  // 1 singular face
  HP_TET_1F_0E_1VA,       // 1 sing vertex in face (V2)
  HP_TET_1F_0E_1VB,       // 1 sing vertex not in face (V1)
  HP_TET_1F_1EA_0V,       // 1 sing edge not in face
  HP_TET_1F_1EB_0V,       // 1 sing edge in face
  HP_TET_2F_0E_0V = 600,  // 2 singular faces

  HP_PRISM = 1000,
  HP_PRISM_SINGEDGE,
  HP_PRISM_SINGEDGE_V12,
  HP_PRISM_SINGEDGE_H1,
  HP_PRISM_SINGEDGE_H12,

  HP_PRISM_1FA_0E_0V,     // 1 singular trig face
  HP_PRISM_2FA_0E_0V,     // 2 singular trig faces
  HP_PRISM_1FB_0E_0V,     // 1 singular quad face  1-2-4-5

  HP_PRISM_1FB_1EA_0V,     // 1 singular quad face, edge is 1-2
  HP_PRISM_1FA_1E_0V, 
  HP_PRISM_2FA_1E_0V, 
  HP_PRISM_1FA_1FB_0E_0V, 
  HP_PRISM_2FA_1FB_0E_0V,
  HP_PRISM_1FA_1FB_1EA_0V, 
  HP_PRISM_1FA_1FB_1EB_0V, 
  HP_PRISM_2FA_1FB_1EA_0V,
  HP_PRISM_1FB_1EC_0V, 
  HP_PRISM_1FA_1FB_1EC_0V, 
  HP_PRISM_2FA_1FB_1EC_0V,
  HP_PRISM_1FB_2EA_0V, 
  HP_PRISM_1FA_1FB_2EA_0V, 
  HP_PRISM_2FA_1FB_2EA_0V,
  HP_PRISM_1FB_2EB_0V,
  HP_PRISM_1FA_1FB_2EB_0V,  
  HP_PRISM_1FA_1FB_2EC_0V, 
  HP_PRISM_2FA_1FB_2EB_0V, 
  HP_PRISM_1FB_3E_0V, 
  HP_PRISM_1FA_1FB_3E_0V, 
  HP_PRISM_2FA_1FB_3E_0V, 
  HP_PRISM_2FB_0E_0V, 
  HP_PRISM_1FA_2FB_0E_0V, 
  HP_PRISM_2FA_2FB_0E_0V,
  HP_PRISM_2FB_1EC_0V, 
  HP_PRISM_1FA_2FB_1EC_0V,
  HP_PRISM_1FA_2FB_1EB_0V,
  HP_PRISM_2FA_2FB_1EC_0V,
  HP_PRISM_2FB_3E_0V, 
  HP_PRISM_1FA_2FB_3E_0V, 
  HP_PRISM_2FA_2FB_3E_0V, 
  HP_PRISM_1FA_2E_0V, 
  HP_PRISM_2FA_2E_0V,
  HP_PRISM_3E_0V, 
  HP_PRISM_1FA_3E_0V, 
  HP_PRISM_2FA_3E_0V,  
  HP_PRISM_3FB_0V, 
  HP_PRISM_1FA_3FB_0V, 
  HP_PRISM_2FA_3FB_0V,  
  HP_PRISM_3E_4EH,
  
 

  /*  HP_PRISM_1FB_1EA_0V,     // 1 singular quad face, edge is 1-4
  HP_PRISM_1FB_1EB_0V,     // 1 singular quad face, edge is 2-5
  HP_PRISM_2F_0E_0V,      // 2 singular quad faces
  */

  HP_PYRAMID = 2000,
  HP_PYRAMID_0E_1V,
  HP_PYRAMID_EDGES,
  HP_PYRAMID_1FB_0E_1VA,  // 1 trig face, top vertex

  HP_HEX = 3000,
  HP_HEX_0E_1V,
  HP_HEX_1E_1V,
  HP_HEX_1E_0V,
  HP_HEX_3E_0V,
  HP_HEX_1F_0E_0V,
  HP_HEX_1FA_1FB_0E_0V, 
};



struct HPRef_Struct {
  HPREF_ELEMENT_TYPE geom;
  int (*splitedges)[3];
  int (*splitfaces)[4];
  int (*splitelements)[5];
  HPREF_ELEMENT_TYPE * neweltypes;
  int (*newels)[8];
};




class HPRefElement
{
private:
  void Reset(void);

public:
  HPRefElement (); 
  HPRefElement(Element & el);
  HPRefElement(Element2d & el);
  HPRefElement(Segment & el);	
  HPRefElement(HPRefElement & el);

  void SetType( HPREF_ELEMENT_TYPE t);
  // HPRefElement(HPRefElement & el, HPREF_ELEMENT_TYPE t); 
	       
  /* HPRefElement(HPRefElement & el, HPREF_ELEMENT_TYPE t)
  { 
    type = t; 
    HPRef_Struct * hprs = Get_HPRef_Struct(t);
    for (int i=0; i<np ; i++) 
      {
	pnums[i] = el[i];
	for(int l=0; l<np; l++) param[i][l] = el.param[i][l]; 
      }
    switch(hprs->geom)
      {
      case HP_SEGM: np=2; sing_edge_left=0; sing_edge_right=0; break; 
      case HP_QUAD: np=4; break; 
      case HP_TRIG: np=3; break; 
      case HP_HEX: np=8; break; 
      case HP_PRISM: np=6; break;
      case HP_TET: np=4; break; 
      case HP_PYRAMID: np=5; break; 
      }
    index = el.index; 
    levelx = el.levelx; 
    levely = el.levely; 
    levelz = el.levelz; 
    type = el.type; 
    coarse_elnr = el.coarse_elnr;
    singedge_left = el.singedge_left; 
    singedge_right = el.singedge_left; 
    } */ 
  
  HPREF_ELEMENT_TYPE type;
  PointIndex pnums[8];
  double param[8][3];
  int index;
  int levelx;
  int levely;
  int levelz;
  int np; 
  int coarse_elnr;
  int domin, domout; // he: needed for segment!! in 3d there should be surf1, surf2!!
  // int coarse_hpelnr; 
  PointIndex & operator[](int i) { return(pnums[i]);}
  PointIndex & PNumMod(int i) { return pnums[(i-1) % np]; };
  PointIndex & PNum(int i) {return pnums[(i-1)]; };
  int GetIndex () const { return index; }; 
  double singedge_left, singedge_right; 
  

  //  EdgePointGeomInfo epgeominfo[2];
  
};



extern void HPRefinement (Mesh & mesh, Refinement * ref, int levels, 
			  double fac1=0.125, bool setorders=true, bool ref_level = false);


#endif

