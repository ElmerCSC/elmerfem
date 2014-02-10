#include <mystdlib.h>
#include "meshing.hpp"
#include "hprefinement.hpp" 

namespace netgen
{

#include "hpref_segm.hpp"
#include "hpref_trig.hpp"  
#include "hpref_quad.hpp"
#include "hpref_tet.hpp"  
#include "hpref_prism.hpp"
#include "hpref_hex.hpp" 
#include "hpref_pyramid.hpp" 
#include "classifyhpel.hpp"  


  void HPRefElement :: Reset(void)
  {
    np = 8; 
    for (int i = 0; i < 8; i++)
      {
	pnums[i] = -1;
	param[i][0] = param[i][1] = param[i][2] = 0;
        domin=-1; domout=-1; // he:
      }
  } 

  HPRefElement :: HPRefElement () 
  {
    Reset();
  }

  HPRefElement :: HPRefElement(Element & el)
  { 
    //Reset();
    np = el.GetNV(); 
    for (int i=0; i<np ; i++) 
      pnums[i] = el[i]; 
    
    index = el.GetIndex(); 
    const Point3d * points = 
      MeshTopology :: GetVertices (el.GetType());
    for(int i=0;i<np;i++)
      for(int l=0;l<3;l++) 
	param[i][l] = points[i].X(l+1); 
    type = HP_NONE; 
    domin=-1; domout=-1; // he: needed for segments
  }

  
  HPRefElement :: HPRefElement(Element2d & el)
  { 
    //Reset();
    np = el.GetNV(); 
    
    for (int i=0; i<np ; i++) 
      pnums[i] = el[i]; 
    
    index = el.GetIndex(); 
    const Point3d * points = 
      MeshTopology :: GetVertices (el.GetType());
    for(int i=0;i<np;i++)
      for(int l=0;l<3;l++) 
	param[i][l] = points[i].X(l+1); 
    type = HP_NONE; 
    domin=-1; domout=-1; // he: needed for segments
  }

  HPRefElement :: HPRefElement(Segment & el)
  { 
    //Reset();
    np = 2; 
    for (int i=0; i<np ; i++) 
      pnums[i] = el[i];
    const Point3d * points = 
      MeshTopology :: GetVertices (SEGMENT); 
    for(int i=0;i<np;i++)
      for(int l=0;l<3;l++) 
        param[i][l] = points[i].X(l+1); 

    /*
    for (int i=0; i<np; i++)
    {
      param[i][0] = i;   
      param[i][1] = -1; param[i][2] = -1;
    }
    */

    singedge_left = el.singedge_left; 
    singedge_right = el.singedge_right; 
    type = HP_NONE; 
    // he: needed for orientation!
    domin = el.domin;
    domout = el.domout;
  }
  
  HPRefElement :: HPRefElement(HPRefElement & el)
  {
    //Reset();
    np = el.np; 
    for (int i=0; i<np ; i++) 
      {
	pnums[i] = el[i];
	for(int l=0; l<3; l++) param[i][l] = el.param[i][l]; 
      }
    index = el.index; 
    levelx = el.levelx; 
    levely = el.levely; 
    levelz = el.levelz; 
    type = el.type; 
    coarse_elnr = el.coarse_elnr; 
    singedge_left = el.singedge_left; 
    singedge_right = el.singedge_right; 
    domin = el.domin; // he: needed for segments
    domout=el.domout;
          
  }

  void HPRefElement :: SetType( HPREF_ELEMENT_TYPE t) 
  {
    type = t; 
    switch(type)
      {
      case HP_SEGM: np=2; break; 
      case HP_TRIG: np=3; break; 
      case HP_QUAD: np=4; break; 
      case HP_TET: np=4; break; 
      case HP_PRISM: np=6; break; 
      case HP_PYRAMID: np=5; break; 
      case HP_HEX: np=8; break;      
      } 
    for(int k=0;k<8;k++)
      {
	pnums[k]=0;
	for(int l=0;l<3;l++) param[k][l]=0.;
      }
  }
  

  HPRef_Struct * Get_HPRef_Struct (HPREF_ELEMENT_TYPE type)
  {
    HPRef_Struct * hps = NULL;

    switch (type)
      {
      case HP_SEGM:
	hps = &refsegm; break;
      case HP_SEGM_SINGCORNERL:
	hps = &refsegm_scl; break;
      case HP_SEGM_SINGCORNERR:
	hps = &refsegm_scr; break;
      case HP_SEGM_SINGCORNERS:
	hps = &refsegm_sc2; break;
	
      case HP_TRIG:
	hps = &reftrig; break;
      case HP_TRIG_SINGCORNER:
	hps = &reftrig_singcorner; break;
      case HP_TRIG_SINGCORNER12:
	hps = &reftrig_singcorner12; break; 
      case HP_TRIG_SINGCORNER123:
	hps = &reftrig_singcorner123; break;
      case HP_TRIG_SINGCORNER123_2D:
	hps = &reftrig_singcorner123_2D; break;
      case HP_TRIG_SINGEDGE:
	hps = &reftrig_singedge; break;
      case HP_TRIG_SINGEDGECORNER1:
	hps = &reftrig_singedgecorner1; break;
      case HP_TRIG_SINGEDGECORNER2:
	hps = &reftrig_singedgecorner2; break;
      case HP_TRIG_SINGEDGECORNER12:
	hps = &reftrig_singedgecorner12; break;
      case HP_TRIG_SINGEDGECORNER3:
	hps = &reftrig_singedgecorner3; break;
      case HP_TRIG_SINGEDGECORNER13:
	hps = &reftrig_singedgecorner13; break;
      case HP_TRIG_SINGEDGECORNER23:
	hps = &reftrig_singedgecorner23; break;
      case HP_TRIG_SINGEDGECORNER123:
	hps = &reftrig_singedgecorner123; break;
      case HP_TRIG_SINGEDGES:
	hps = &reftrig_singedges; break; 
      case HP_TRIG_SINGEDGES2: 
	hps = &reftrig_singedges2; break;
      case HP_TRIG_SINGEDGES3:
	hps = &reftrig_singedges3; break;
      case HP_TRIG_SINGEDGES23:  
	hps = &reftrig_singedges23; break;
      case HP_TRIG_3SINGEDGES:
	hps = &reftrig_3singedges; break;
 
 
      case HP_QUAD:
	hps = &refquad; break;
      case HP_DUMMY_QUAD_SINGCORNER:
	hps = &refdummyquad_singcorner; break;
      case HP_QUAD_SINGCORNER:
	hps = &refquad_singcorner; break;
      case HP_QUAD_SINGEDGE:
	hps = &refquad_singedge; break;

      case HP_QUAD_0E_2VA:
	hps = &refquad_0e_2va; break;
      case HP_QUAD_0E_2VB:
	hps = &refquad_0e_2vb; break;

      case HP_QUAD_0E_3V:
	hps = &refquad_0e_3v; break;
      case HP_QUAD_0E_4V:
	hps = &refquad_0e_4v; break;

      case HP_QUAD_1E_1VA:
	hps = &refquad_1e_1va; break;
      case HP_QUAD_1E_1VB:
	hps = &refquad_1e_1vb; break;
      case HP_QUAD_1E_1VC:
	hps = &refquad_1e_1vc; break;
      case HP_QUAD_1E_1VD:
	hps = &refquad_1e_1vd; break;

      case HP_QUAD_1E_2VA:
	hps = &refquad_1e_2va; break;
      case HP_QUAD_1E_2VB:
	hps = &refquad_1e_2vb; break;
      case HP_QUAD_1E_2VC:
	hps = &refquad_1e_2vc; break;
      case HP_QUAD_1E_2VD:
	hps = &refquad_1e_2vd; break;
      case HP_QUAD_1E_2VE:
	hps = &refquad_1e_2ve; break;
      case HP_QUAD_1E_2VF:
	hps = &refquad_1e_2vf; break;

      case HP_QUAD_1E_3VA:
	hps = &refquad_1e_3va; break;
      case HP_QUAD_1E_3VB:
	hps = &refquad_1e_3vb; break;
      case HP_QUAD_1E_3VC:
	hps = &refquad_1e_3vc; break;
      case HP_QUAD_1E_3VD:
	hps = &refquad_1e_3vd; break;
      case HP_QUAD_1E_4V:
	hps = &refquad_1e_4v; break;


      case HP_QUAD_2E:
	hps = &refquad_2e; break;
      case HP_QUAD_2E_1VA:
	hps = &refquad_2e_1va; break;
      case HP_QUAD_2E_1VB:
	hps = &refquad_2e_1vb; break;
      case HP_QUAD_2E_1VC:
	hps = &refquad_2e_1vc; break;
      case HP_QUAD_2E_2VA:
	hps = &refquad_2e_2va; break;
      case HP_QUAD_2E_2VB:
	hps = &refquad_2e_2vb; break;
      case HP_QUAD_2E_2VC:
	hps = &refquad_2e_2vc; break;
      case HP_QUAD_2E_3V:
	hps = &refquad_2e_3v; break;

      case HP_QUAD_2EB_0V:
	hps = &refquad_2eb_0v; break;

      case HP_QUAD_2EB_1VA:
	hps = &refquad_2eb_1va; break;
      case HP_QUAD_2EB_1VB:
	hps = &refquad_2eb_1vb; break;


      case HP_QUAD_2EB_2VA:
	hps = &refquad_2eb_2va; break;
      case HP_QUAD_2EB_2VB:
	hps = &refquad_2eb_2vb; break;
      case HP_QUAD_2EB_2VC:
	hps = &refquad_2eb_2vc; break;
      case HP_QUAD_2EB_2VD:
	hps = &refquad_2eb_2vd; break;

      case HP_QUAD_2EB_3VA:
	hps = &refquad_2eb_3va; break;
      case HP_QUAD_2EB_3VB:
	hps = &refquad_2eb_3vb; break;

      case HP_QUAD_2EB_4V:
	hps = &refquad_2eb_4v; break;

      case HP_QUAD_3E:
	hps = &refquad_3e; break;
      case HP_QUAD_3E_3VA:
	hps = &refquad_3e_3va; break;
      case HP_QUAD_3E_3VB:
	hps = &refquad_3e_3vb; break;
      case HP_QUAD_3E_4V:
	hps = &refquad_3e_4v; break;


      case HP_QUAD_4E:
	hps = &refquad_4e; break;


      case HP_TET:
	hps = &reftet; break;
      case HP_TET_0E_1V:
	hps = &reftet_0e_1v; break;
      case HP_TET_0E_2V:
	hps = &reftet_0e_2v; break;
      case HP_TET_0E_3V:
	hps = &reftet_0e_3v; break;
      case HP_TET_0E_4V:
	hps = &reftet_0e_4v; break;

      case HP_TET_1E_0V:      
	hps = &reftet_1e_0v; break;
      case HP_TET_1E_1VA:
	hps = &reftet_1e_1va; break;
      case HP_TET_1E_1VB:
	hps = &reftet_1e_1vb; break;

      case HP_TET_1E_2VA:
	hps = &reftet_1e_2va; break;
      case HP_TET_1E_2VB:
	hps = &reftet_1e_2vb; break;
      case HP_TET_1E_2VC:
	hps = &reftet_1e_2vc; break;
      case HP_TET_1E_2VD:
	hps = &reftet_1e_2vd; break;

      case HP_TET_1E_3VA:
	hps = &reftet_1e_3va; break;
      case HP_TET_1E_3VB:
	hps = &reftet_1e_3vb; break;
      case HP_TET_1E_4V:
	hps = &reftet_1e_4v; break;

      case HP_TET_2EA_0V:
	hps = &reftet_2ea_0v; break;
      case HP_TET_2EA_1VB:
	hps = &reftet_2ea_1vb; break;
      case HP_TET_2EA_1VC:
	hps = &reftet_2ea_1vc; break;
      case HP_TET_2EA_1VA:
	hps = &reftet_2ea_1va; break;
      case HP_TET_2EA_2VA:
	hps = &reftet_2ea_2va; break;
      case HP_TET_2EA_2VB:
	hps = &reftet_2ea_2vb; break;
      case HP_TET_2EA_2VC:
	hps = &reftet_2ea_2vc; break;
      case HP_TET_2EA_3V:
	hps = &reftet_2ea_3v; break;

      case HP_TET_2EB_0V:
	hps = &reftet_2eb_0v; break;
      case HP_TET_2EB_1V:
	hps = &reftet_2eb_1v; break;
      case HP_TET_2EB_2VA:
	hps = &reftet_2eb_2va; break;
      case HP_TET_2EB_2VB:
	hps = &reftet_2eb_2vb; break;
      case HP_TET_2EB_2VC:
	hps = &reftet_2eb_2vc; break;
      case HP_TET_2EB_3V:
	hps = &reftet_2eb_3v; break;
      case HP_TET_2EB_4V:
	hps = &reftet_2eb_4v; break;


      case HP_TET_3EA_0V:
	hps = &reftet_3ea_0v; break;
      case HP_TET_3EA_1V:
	hps = &reftet_3ea_1v; break;
      case HP_TET_3EA_2V:
	hps = &reftet_3ea_2v; break;
      case HP_TET_3EA_3V:
	hps = &reftet_3ea_3v; break;

      case HP_TET_3EB_0V:
	hps = &reftet_3eb_0v; break;
      case HP_TET_3EB_1V:
	hps = &reftet_3eb_1v; break;
      case HP_TET_3EB_2V:
	hps = &reftet_3eb_2v; break;
      case HP_TET_3EC_0V:
	hps = &reftet_3ec_0v; break;
      case HP_TET_3EC_1V:
	hps = &reftet_3ec_1v; break;
      case HP_TET_3EC_2V:
	hps = &reftet_3ec_2v; break;


      case HP_TET_1F_0E_0V:
	hps = &reftet_1f_0e_0v; break;
      case HP_TET_1F_0E_1VA:
	hps = &reftet_1f_0e_1va; break;
      case HP_TET_1F_0E_1VB:
	hps = &reftet_1f_0e_1vb; break;
      case HP_TET_1F_1EA_0V:
	hps = &reftet_1f_1ea_0v; break;
      case HP_TET_1F_1EB_0V:
	hps = &reftet_1f_1eb_0v; break;
      case HP_TET_2F_0E_0V:
	hps = &reftet_2f_0e_0v; break;


      case HP_PRISM:
	hps = &refprism; break;
      case HP_PRISM_SINGEDGE:
	hps = &refprism_singedge; break;
	//      case HP_PRISM_SINGEDGE_H1:
	//	hps = &refprism_singedge_h1; break;
	// case HP_PRISM_SINGEDGE_H12:
	//	hps = &refprism_singedge_h12; break;
      case HP_PRISM_SINGEDGE_V12:
	hps = &refprism_singedge_v12; break;
	

      case HP_PRISM_1FA_0E_0V:
	hps = &refprism_1fa_0e_0v; break;
      case HP_PRISM_2FA_0E_0V:
	hps = &refprism_2fa_0e_0v; break;
      case HP_PRISM_1FB_0E_0V:
	hps = &refprism_1fb_0e_0v; break;
      case HP_PRISM_1FB_1EA_0V: 
	hps = &refprism_1fb_1ea_0v; break;
 
      case HP_PRISM_1FA_1E_0V:
	hps = &refprism_1fa_1e_0v; break;
      case HP_PRISM_2FA_1E_0V:
	hps = &refprism_2fa_1e_0v; break; 
      case HP_PRISM_1FA_1FB_0E_0V: 
	hps = &refprism_1fa_1fb_0e_0v; break;
      case HP_PRISM_2FA_1FB_0E_0V: 
	hps = &refprism_2fa_1fb_0e_0v; break; 
      case HP_PRISM_1FA_1FB_1EA_0V: 
	hps = &refprism_1fa_1fb_1ea_0v; break; 
      case HP_PRISM_1FA_1FB_1EB_0V: 
	hps = &refprism_1fa_1fb_1eb_0v; break; 
      case HP_PRISM_2FA_1FB_1EA_0V:  
	hps = &refprism_2fa_1fb_1ea_0v; break; 
      case HP_PRISM_1FB_1EC_0V: 
	hps = &refprism_1fb_1ec_0v; break; 
      case HP_PRISM_1FA_1FB_1EC_0V: 
	hps = &refprism_1fa_1fb_1ec_0v; break; 
      case HP_PRISM_2FA_1FB_1EC_0V: 
	hps = &refprism_2fa_1fb_1ec_0v; break;  
      case HP_PRISM_1FB_2EA_0V: 
	hps = &refprism_1fb_2ea_0v; break; 
      case HP_PRISM_1FA_1FB_2EA_0V:   
	hps = &refprism_1fa_1fb_2ea_0v; break; 
      case HP_PRISM_2FA_1FB_2EA_0V:  
	hps = &refprism_2fa_1fb_2ea_0v; break;   
      case HP_PRISM_1FB_2EB_0V: 
	hps = &refprism_1fb_2eb_0v; break; 
      case HP_PRISM_1FA_1FB_2EB_0V:   
	hps = &refprism_1fa_1fb_2eb_0v; break;  
      case HP_PRISM_1FA_1FB_2EC_0V: 
	hps = &refprism_1fa_1fb_2ec_0v; break; 
      case HP_PRISM_2FA_1FB_2EB_0V: 
	hps = &refprism_2fa_1fb_2eb_0v; break;  
      case HP_PRISM_1FB_3E_0V: 
	hps = &refprism_1fb_3e_0v; break; 
      case HP_PRISM_1FA_1FB_3E_0V:  
	hps = &refprism_1fa_1fb_3e_0v; break;  
      case HP_PRISM_2FA_1FB_3E_0V:  
        hps = &refprism_2fa_1fb_3e_0v; break; 
      case HP_PRISM_2FB_0E_0V: 
	hps = &refprism_2fb_0e_0v; break;  
      case HP_PRISM_1FA_2FB_0E_0V: 
	hps = &refprism_1fa_2fb_0e_0v; break; 
      case HP_PRISM_2FA_2FB_0E_0V:    
        hps = &refprism_2fa_2fb_0e_0v; break;  
      case HP_PRISM_2FB_1EC_0V:  
	hps = &refprism_2fb_1ec_0v; break;  
      case HP_PRISM_1FA_2FB_1EC_0V: 
        hps = &refprism_1fa_2fb_1ec_0v; break; 
      case HP_PRISM_2FA_2FB_1EC_0V: 
	hps = &refprism_2fa_2fb_1ec_0v; break; 
      case HP_PRISM_1FA_2FB_1EB_0V: 
	hps = &refprism_1fa_2fb_1eb_0v; break; 
      case HP_PRISM_2FB_3E_0V:    
	hps = &refprism_2fb_3e_0v; break;
      case HP_PRISM_1FA_2FB_3E_0V: 
	hps = &refprism_1fa_2fb_3e_0v; break;
      case HP_PRISM_2FA_2FB_3E_0V:  
	hps = &refprism_2fa_2fb_3e_0v; break;
      case HP_PRISM_1FA_2E_0V:   
	hps = &refprism_1fa_2e_0v; break; 
      case HP_PRISM_2FA_2E_0V:  
	hps = &refprism_2fa_2e_0v; break;   
      case HP_PRISM_3E_0V:  
	hps = &refprism_3e_0v; break;
      case HP_PRISM_1FA_3E_0V:   
	hps = &refprism_1fa_3e_0v; break; 
      case HP_PRISM_2FA_3E_0V:  
	hps = &refprism_2fa_3e_0v; break;   
      case HP_PRISM_3FB_0V:  
	hps = &refprism_3fb_0v; break;
      case HP_PRISM_1FA_3FB_0V:   
	hps = &refprism_1fa_3fb_0v; break; 
      case HP_PRISM_2FA_3FB_0V:  
	hps = &refprism_2fa_3fb_0v; break;   
	//  case HP_PRISM_3E_4EH:
	//  hps = &refprism_3e_4eh; break;   
	
	
	/*case HP_PRISM_1FB_1EB_0V:
	hps = &refprism_1fb_1eb_0v; break;
      case HP_PRISM_2F_0E_0V:
	hps = &refprism_2f_0e_0v; break;
	*/
	
	
      case HP_PYRAMID:
	hps = &refpyramid; break;
      case HP_PYRAMID_0E_1V:
	hps = &refpyramid_0e_1v; break;
      case HP_PYRAMID_EDGES:
	hps = &refpyramid_edges; break;
      case HP_PYRAMID_1FB_0E_1VA:
	hps = &refpyramid_1fb_0e_1va; break;

	
      case HP_HEX:
	hps = &refhex; break;
      case HP_HEX_0E_1V:
	hps = &refhex_0e_1v; break;
      case HP_HEX_1E_1V:
	hps = &refhex_1e_1v; break;
      case HP_HEX_1E_0V:
	hps = &refhex_1e_0v; break;
      case HP_HEX_3E_0V:
	hps = &refhex_3e_0v; break;

      case HP_HEX_1F_0E_0V:
	hps = &refhex_1f_0e_0v; break;
      case HP_HEX_1FA_1FB_0E_0V: 
	hps = &refhex_1fa_1fb_0e_0v; break; 
      }

    /*
    if (type != HP_TET_1E_4V && type != HP_TET_1E_2VD)
      {
	if (hps->geom == HP_TET)
	  hps = &reftet;
	if (hps->geom == HP_TRIG)
	  hps = &reftrig;
      }
    */

    if (!hps)
      {
	cout << "Attention hps : hp-refinement not implemented for case " << type << endl;
	PrintSysError ("hp-refinement not implemented for case ", type);
      }

    return hps;
  }

  bool CheckSingularities(Mesh & mesh, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoiclt_dom, 
		       BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
			INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint, int & levels, int & act_ref); 

  bool ClassifyHPElements (Mesh & mesh, ARRAY<HPRefElement> & elements, int & act_ref, int & levels);
  
  
  void  InitHPElements(Mesh & mesh, ARRAY<HPRefElement> & elements) 
  { 
    for(ElementIndex i=0;i<mesh.GetNE();i++) 
      {
	HPRefElement hpel(mesh[i]); 
	hpel.coarse_elnr=i; 
	
	switch (mesh[i].GetType()) 
	  { 
	  case PRISM:
	    hpel.type = HP_PRISM; 
	    break; 
	  case HEX:
	    hpel.type = HP_HEX; 
	    break; 
	  case TET: 
	    hpel.type = HP_TET; 
	    break; 
	  case PYRAMID: 
	    hpel.type = HP_PYRAMID; 
	    break; 
	  } 
	elements.Append(hpel); 
      }
	    
    for(SurfaceElementIndex i=0;i<mesh.GetNSE();i++)
      {
	HPRefElement hpel(mesh.SurfaceElement(i));
	hpel.coarse_elnr = i; 
	switch(mesh.SurfaceElement(i).GetType())
	  { 
	  case TRIG: 
	    hpel.type = HP_TRIG;
	    break; 
	  case QUAD: 
	    hpel.type = HP_QUAD; 
	    break; 
	  } 
	elements.Append(hpel);
      } 
        
    for(int i=1;i<=mesh.GetNSeg();i++) 
      { 
	Segment & seg = mesh.LineSegment(i); 
	HPRefElement hpel(seg); 
	hpel.coarse_elnr = i-1; 
	hpel.type = HP_SEGM; 
	hpel.index = seg.edgenr + 10000*seg.si; 
	if(seg.edgenr >= 10000)
	  {
	    throw NgException("assumption that seg.edgenr < 10000 is wrong");
	  }
	elements.Append(hpel); 

      }
  }

 
 
  /* *******************************  DoRefinement *************************************** */
  void DoRefinement (Mesh & mesh, ARRAY<HPRefElement> & elements,
		     Refinement * ref, double fac1) 
  {
    elements.SetAllocSize (5 * elements.Size());
    INDEX_2_HASHTABLE<int> newpts(elements.Size()+1);
    INDEX_3_HASHTABLE<int> newfacepts(elements.Size()+1);

    // prepare new points  
    
    fac1 = max(0.001,min(0.33,fac1));
    cout << " in HP-REFINEMENT with fac1 " << fac1 << endl; 
    *testout << " in HP-REFINEMENT with fac1 " << fac1 <<  endl; 
   

    int oldelsize = elements.Size();
       
    for (int i = 0; i < oldelsize; i++)
      {
	HPRefElement & el = elements[i]; 
	HPRef_Struct * hprs = Get_HPRef_Struct (el.type);
		
	if (!hprs) 
	  {
	    cout << "Refinementstruct not defined for element " << el.type << endl;
	    continue;
	  }

	int j = 0;
	while (hprs->splitedges[j][0])
	  {
	    INDEX_2 i2(el.pnums[hprs->splitedges[j][0]-1],
		       el.pnums[hprs->splitedges[j][1]-1]);
	    if (!newpts.Used (i2))
	      {
		Point<3> np; 
		for( int l=0;l<3;l++)
		  np(l) = (1-fac1)*mesh.Point(i2.I1())(l) 
		    + fac1 * mesh.Point(i2.I2())(l); 
	
		int npi = mesh.AddPoint (np);
		newpts.Set (i2, npi);
	      }
	    j++;
	  }
	
	j = 0;
	if (hprs->splitfaces)
	  while (hprs->splitfaces[j][0])
	    {
	      INDEX_3 i3(el.pnums[hprs->splitfaces[j][0]-1],
			 el.pnums[hprs->splitfaces[j][1]-1],
			 el.pnums[hprs->splitfaces[j][2]-1]);

	      if (i3.I2() > i3.I3()) Swap (i3.I2(), i3.I3());
	      
	      if (!newfacepts.Used (i3))
		{
		  Point<3> np; 
		  	for( int l=0;l<3;l++)
			  np(l) = (1-2*fac1)*mesh.Point(i3.I1())(l) 
			    + fac1*mesh.Point(i3.I2())(l)  + fac1*mesh.Point(i3.I3())(l);  
		  int npi = mesh.AddPoint (np);
		  newfacepts.Set (i3, npi);
		}
	      j++;
	    }
      }
     
    for (int i = 0; i < oldelsize; i++)
      {
	HPRefElement el = elements[i];
	HPRef_Struct * hprs = Get_HPRef_Struct (el.type);
	int newlevel = el.levelx + 1;

	int oldnp(0);
	switch (hprs->geom)
	  {
	  case HP_SEGM: oldnp = 2; break;
	  case HP_TRIG: oldnp = 3; break;
	  case HP_QUAD: oldnp = 4; break;
	  case HP_TET: oldnp = 4; break;
	  case HP_PYRAMID: oldnp = 5; break;
	  case HP_PRISM: oldnp = 6; break;
	  case HP_HEX: oldnp = 8; break;
	  }


	if (el.type == HP_SEGM ||
	    el.type == HP_TRIG ||
	    el.type == HP_QUAD ||
	    el.type == HP_TET ||
	    el.type == HP_PRISM ||
	    el.type == HP_HEX || 
	    el.type == HP_PYRAMID)
	  newlevel = el.levelx;

	if (!hprs) continue;

	int newpnums[64];
	double newparam[64][3];

	int j;
	for (j = 0; j < oldnp; j++)
	  {
	    newpnums[j] = el.pnums[j];
	    for (int l = 0; l < 3; l++)
	      newparam[j][l] = el.param[j][l];
	  }

	// split edges, incl. transferring curvature
	j = 0;
	while (hprs->splitedges[j][0])
	  {
	    INDEX_2 i2(el.pnums[hprs->splitedges[j][0]-1],
		       el.pnums[hprs->splitedges[j][1]-1]);

	    int npi = newpts.Get(i2);
	    newpnums[hprs->splitedges[j][2]-1] = npi;

	    for (int l = 0; l < 3; l++)
	      newparam[hprs->splitedges[j][2]-1][l] =
		(1-fac1) * el.param[hprs->splitedges[j][0]-1][l] + 
		fac1 * el.param[hprs->splitedges[j][1]-1][l];
	      
	    j++;
	  }

	// split faces
	j = 0;
	if (hprs->splitfaces)
	  while (hprs->splitfaces[j][0])
	    {
	      INDEX_3 i3(el.pnums[hprs->splitfaces[j][0]-1],
			 el.pnums[hprs->splitfaces[j][1]-1],
			 el.pnums[hprs->splitfaces[j][2]-1]);
	      if (i3.I2() > i3.I3())
		Swap (i3.I2(), i3.I3());
	      int npi = newfacepts.Get(i3);
	      newpnums[hprs->splitfaces[j][3]-1] = npi;
	    

	      for (int l = 0; l < 3; l++)
		newparam[hprs->splitfaces[j][3]-1][l] =
		  (1-2*fac1) * el.param[hprs->splitfaces[j][0]-1][l] + 
		  fac1 * el.param[hprs->splitfaces[j][1]-1][l] + 
		  fac1 * el.param[hprs->splitfaces[j][2]-1][l];
	      j++;
	    }
	// split elements
	j = 0;
	if (hprs->splitelements)
	  while (hprs->splitelements[j][0])
	    {
	      //int pi1 = el.pnums[hprs->splitelements[j][0]-1];
	      Point<3> np; 
	      	for( int l=0;l<3;l++)
		  np(l) = (1-3*fac1)* mesh.Point(el.pnums[hprs->splitelements[j][0]-1])(l) 
		    + fac1* mesh.Point(el.pnums[hprs->splitelements[j][1]-1])(l)
		    + fac1* mesh.Point(el.pnums[hprs->splitelements[j][2]-1])(l)
		    + fac1* mesh.Point(el.pnums[hprs->splitelements[j][3]-1])(l); 
	      
	      int npi = mesh.AddPoint (np);
	      
	      newpnums[hprs->splitelements[j][4]-1] = npi;
	      
	  	    
	      for (int l = 0; l  < 3; l++)
		newparam[hprs->splitelements[j][4]-1][l] =
		  (1-3*fac1) * el.param[hprs->splitelements[j][0]-1][l] + 
		  fac1 * el.param[hprs->splitelements[j][1]-1][l] + 
		  fac1 * el.param[hprs->splitelements[j][2]-1][l] + 
		  fac1 * el.param[hprs->splitelements[j][3]-1][l];

	      j++;
	    }
 
	j = 0;

	/*
	*testout << " newpnums = ";
	for (int hi = 0; hi < 64; hi++)
	  *testout << newpnums[hi] << " ";
	*testout << endl;
	*/

	while (hprs->neweltypes[j])
	  {
	    HPRef_Struct * hprsnew = Get_HPRef_Struct (hprs->neweltypes[j]);
	    HPRefElement newel(el);

	    newel.type = hprs->neweltypes[j]; 
            
	    // newel.index = elements[i].index;
	    // newel.coarse_elnr = elements[i].coarse_elnr;
	    newel.levelx = newel.levely = newel.levelz = newlevel;
            switch(hprsnew->geom) 
	      {
	      case HP_SEGM: newel.np=2; break; 
	      case HP_QUAD: newel.np=4; break; 
	      case HP_TRIG: newel.np=3; break; 
	      case HP_HEX: newel.np=8; break; 
	      case HP_PRISM: newel.np=6; break;
	      case HP_TET: newel.np=4; break; 
	      case HP_PYRAMID: newel.np=5; break; 
	      }

	    for (int k = 0; k < newel.np; k++)
	      newel.pnums[k] = newpnums[hprs->newels[j][k]-1];
	    
	    /*
	    *testout  << " newel pnums " ; 
	    for (int k = 0; k < newel.np; k++)  
	      *testout  << newel.pnums[k] << "\t"; 
	    *testout << endl; 
	    */

	    for (int k = 0; k < newel.np; k++)  
	      { 
		for (int l = 0; l < 3; l++)
		  { 
		    newel.param[k][l] = newparam[hprs->newels[j][k]-1][l];
		    //    *testout << newel.param[k][l] << " \t ";
		  } 
		// *testout << endl; 
	      } 
	    
	    if (j == 0) 
	      elements[i] = newel; // overwrite old element
	    else        
	      elements.Append (newel);
	    j++;
	  }
      } 
  }






  /* ************************** DoRefineDummies ******************************** */

  void DoRefineDummies (Mesh & mesh, ARRAY<HPRefElement> & elements,
			Refinement * ref)
  {
    int oldelsize = elements.Size();

    for (int i = 0; i < oldelsize; i++)
      {
	HPRefElement el = elements[i];

	HPRef_Struct * hprs = Get_HPRef_Struct (el.type);
	if (!hprs) continue;

	if (el.type != HP_DUMMY_QUAD_SINGCORNER &&
	    el.type != HP_PYRAMID_EDGES &&
	    el.type != HP_PYRAMID_0E_1V &&
	    el.type != HP_HEX_0E_1V &&
	    el.type != HP_HEX_1E_1V &&
	    el.type != HP_HEX_1E_0V &&
	    el.type != HP_HEX_3E_0V
	    ) continue;

	int newlevel = el.levelx;

	int newpnums[8];
	int j;
	for (j = 0; j < 8; j++)
	  newpnums[j] = el.pnums[j];

	double newparam[8][3];
	for (j = 0; j < 8; j++)
	  for (int k = 0; k < 3; k++)
	    newparam[j][k] = el.param[j][k];

	j = 0;
	while (hprs->neweltypes[j])
	  {
	    HPRef_Struct * hprsnew = Get_HPRef_Struct (hprs->neweltypes[j]);
	    HPRefElement newel(el);
	    switch(hprsnew->geom)
	      {
	      case HP_SEGM: newel.np=2; break; 
	      case HP_QUAD: newel.np=4; break; 
	      case HP_TRIG: newel.np=3; break; 
	      case HP_HEX: newel.np=8; break; 
	      case HP_PRISM: newel.np=6; break;
	      case HP_TET: newel.np=4; break; 
	      case HP_PYRAMID: newel.np=5; break; 
	      }
	    newel.type = hprs->neweltypes[j];
	    for (int k = 0; k < 8; k++)
	      newel.pnums[k] = newpnums[hprs->newels[j][k]-1];
	    newel.index = el.index;
	    newel.coarse_elnr = el.coarse_elnr;
	    newel.levelx = newel.levely = newel.levelz = newlevel;

	    for (int k = 0; k < 8; k++)
	      for (int l = 0; l < 3; l++)
		newel.param[k][l] = newparam[hprs->newels[j][k]-1][l];
		
	    if (j == 0)
	      elements[i] = newel;
	    else
	      elements.Append (newel);
	    j++;
	  }
      }
  }







  void SubdivideDegeneratedHexes (Mesh & mesh, ARRAY<HPRefElement> & elements, double fac1)
  {
    int oldne = elements.Size();
    for (int i = 0; i < oldne; i++)
      if (Get_HPRef_Struct (elements[i].type)->geom == HP_HEX)
	{
	  bool common = 0;
	  for (int j = 0; j < 8; j++)
	    for (int k = 0; k < j; k++)
	      if (elements[i].pnums[j] == elements[i].pnums[k])
		common = 1;
	  if (common)
	    {

               
	      cout << " Degenerate Hex found " << endl; 
              *testout << " Degenerate Hex found " << endl; 
	      HPRefElement el = elements[i];
	      HPRefElement newel = el;

	      Point<3> center(0,0,0);
	      double newparam[3] = { 0, 0, 0 };

	      for (int j = 0; j < 8; j++)
		{
		  
		  
		  center += 0.125 * Vec<3>(mesh[el.pnums[j]]); 
		  // 0.125 originates form 8 points not from fac1;
                  
		  for (int l = 0; l < 3; l++)
		    newparam[l] += 0.125 * el.param[j][l];
                  
		}

	      int npi = mesh.AddPoint (center);

	      const ELEMENT_FACE * faces = MeshTopology::GetFaces (HEX);

	      for (int j = 0; j < 6; j++)  
		{
		  ARRAY<int> pts;
		  for (int k = 0; k < 4; k++)
		    {
		      bool same = 0;
		      for (int l = 0; l < pts.Size(); l++)
			if (el.pnums[pts[l]] == el.pnums[faces[j][k]-1])
			  same = 1;
		      if (!same)
			pts.Append (faces[j][k]-1);

		    }
		  
		  
		  if (pts.Size() == 3) // TrigFace -> TET 
		    {
		      
		      for (int k = 0; k < 3; k++)
			{
			  newel.pnums[k] = el.pnums[pts[2-k]];
			  for (int l = 0; l < 3; l++)
			    newel.param[k][l] = el.param[pts[2-k]][l];
			}
		      newel.pnums[3] = npi;
		      for (int l = 0; l < 3; l++)
			newel.param[3][l] = newparam[l];

		      newel.type = HP_TET;
		      newel.np = 4; 
		    }
		  else
		    {
		      for (int k = 0; k < 4; k++)
			{
			  newel.pnums[k] = el.pnums[pts[3-k]];
			  for (int l = 0; l < 3; l++)
			    newel.param[k][l] = el.param[pts[3-k]][l];
			}

		      newel.pnums[4] = npi;
		      for (int l = 0; l < 3; l++)
			newel.param[4][l] = newparam[l];

		      newel.type = HP_PYRAMID;
		      newel.np = 5;
		    }
		  
		  if (j == 0)
		    elements[i] = newel;
		  else
		    elements.Append (newel); 

		  
		}

	      /*     const ELEMENT_EDGE * edges = MeshTopology::GetEdges (HEX);
	       
		for(int k=0;k<12;k++) 
		  { 
		    int e[2];  
		    for(int l=0;l<2;l++) e[l] = edges[k][l]-1; 
		    if(el.PNum(e[0]+1)!=el.PNum(e[1]+1)) 
		      { 
			newel.SetType(HP_SEGM);
			for(int l=0;l<2;l++) 
			  { 
			    newel.pnums[0] = el.PNum(e[l]+1); 
			    newel.pnums[1] = npi; 
			    for(int j=0;j<3;j++) 
			      {
				//	newel.param[0][j] = el.param[e[l]][j]; 
				//	newel.param[1][j] = newparam[j]; 
			      } 
			    
			    elements.Append(newel);
			  }
			newel.SetType(HP_TRIG);
			newel.pnums[0] = el.PNum(e[0]+1); 			
			newel.pnums[1] = el.PNum(e[1]+1); 			
			newel.pnums[2] = npi; 
			
			*testout << "DEGHEX TRIG :: newpnums " << newel.pnums[0] << "\t"  << newel.pnums[1] << "\t"  << newel.pnums[2] << endl;  
	cout << "DEGHEX TRIG :: newpnums " << newel.pnums[0] << "\t"  << newel.pnums[1] << "\t"  << newel.pnums[2] << endl;  
			for(int j=0;j<3;j++) 
			  {
			    // newel.param[0][j] = el.param[e[0]][j]; 
			    //   newel.param[1][j] = el.param[e[1]][j]; 
			    //   newel.param[2][j] = newparam[j]; 
			  } 
			
			elements.Append(newel);
		      }
			
		      }*/
	    }
	}
  }


  void CalcStatistics (ARRAY<HPRefElement> & elements)
  {
    return;
#ifdef ABC    
    int i, p;
    int nsegm = 0, ntrig = 0, nquad = 0;
    int nhex = 0, nprism = 0, npyramid = 0, ntet = 0;
    int maxlevel = 0;

    for (i = 1; i <= elements.Size(); i++)
      {
	const HPRefElement & el = elements.Get(i);
	maxlevel = max2 (el.level, maxlevel);
	switch (Get_HPRef_Struct (el.type)->geom)
	  {
	  case HP_SEGM:

	    {
	      nsegm++;
	      break;
	    }
	  case HP_TRIG:
	    {
	      ntrig ++;
	      break;
	    }
	  case HP_QUAD:
	    {
	      nquad++;
	      break;
	    }
	  case HP_TET:
	    {
	      ntet++;
	      break;
	    }

	  case HP_PRISM:
	    {
	      nprism++;
	      break;
	    }

	  case HP_PYRAMID:
	    {
	      npyramid++;
	      break;
	    }

	  case HP_HEX:
	    {	
	      nhex++;
	      break;
	    }

	  default:
	    {
	      cerr << "statistics error, unknown element type" << endl;
	    }
	  }
      }

    cout << "level = " << maxlevel << endl;
    cout << "nsegm = " << nsegm << endl;
    cout << "ntrig = " << ntrig << ", nquad = " << nquad << endl;
    cout << "ntet = " << ntet << ", npyr = " << npyramid
	 << ", nprism = " << nprism << ", nhex = " << nhex << endl;

    return;

    double memcost = 0, cpucost = 0;
    for (p = 1; p <= 20; p++)
      {
	memcost = (ntet + nprism + nhex) * pow (static_cast<double>(p), 6.0);
	cpucost = (ntet + nprism + nhex) * pow (static_cast<double>(p), 9.0);
	cout << "costs for p = " << p << ": mem = " << memcost << ", cpu = " << cpucost << endl;
      }

    double memcosttet = 0;
    double memcostprism = 0;
    double memcosthex = 0;
    double memcostsctet = 0; 
    double memcostscprism = 0;
    double memcostschex = 0;
    double cpucosttet = 0;
    double cpucostprism = 0;
    double cpucosthex = 0;

    for (i = 1; i <= elements.Size(); i++)
      {
	const HPRefElement & el = elements.Get(i);
	switch (el.type)
	  {
	  case HP_TET:
	  case HP_TET_0E_1V:
	  case HP_TET_1E_0V:
	  case HP_TET_1E_1VA:
	    {
	      int p1 = maxlevel - el.level + 1;
	      (*testout) << "p1 = " << p1 << ", P1^6 = " << pow (static_cast<double>(p1), 6.0)
			 << " (p1-3)^6 = " << pow ( static_cast<double>(max2(p1-3, 0)), 6.0) 
			 << " p1^3 = " << pow ( static_cast<double>(p1), 3.0) 
			 << " (p1-3)^3 = " << pow ( static_cast<double>(p1-3), 3.0) 
			 << " [p1^3-(p1-3)^3]^2 = " << sqr (pow (static_cast<double>(p1),3.0) - pow ( static_cast<double>(p1-3), 3.0))
			 << endl;

	      p1 /= 2 +1;
	      memcosttet += pow (static_cast<double>(p1), 6.0);
	      memcostsctet += pow (static_cast<double>(p1), 6.0) - pow ( static_cast<double>(max2(p1-3, 1)), 6.0);
	      cpucosttet += pow (static_cast<double>(p1), 9.0);
	      break;
	    }
	  case HP_PRISM:
	  case HP_PRISM_SINGEDGE:
	    {
	      int p1 = maxlevel - el.level + 1;
	      p1 /= 2 +1;
	      memcostprism += pow (static_cast<double>(p1), 6.0);
	      memcostscprism += pow (static_cast<double>(p1), 6.0) - pow ( static_cast<double>(max2(p1-3, 1)), 6.0);
	      cpucostprism += pow (static_cast<double>(p1), 9.0);
	      break;
	    }
	  case HP_HEX:
	    {	
	      int p1 = maxlevel - el.level + 1;
	      int p2 = maxlevel;
	      p1 /= 2 +1;
	      p2 /= 2 +1;
	      memcosthex += pow (static_cast<double>(p1), 4.0) * pow (static_cast<double>(p2), 2.0);
	      memcostschex += pow (static_cast<double>(p1), 6.0) - pow ( static_cast<double>(max2(p1-2, 0)), 6.0);
	      cpucosthex += pow (static_cast<double>(p1), 6.0) * pow (static_cast<double>(p2), 3.0);
	      break;
	    }
	  default:
	    ;
	  }
      }
    cout << "TET: hp-memcost = " << memcosttet 
	 << ", scmemcost = " << memcostsctet
	 << ", cpucost = " << cpucosttet
	 << endl;
    cout << "PRI: hp-memcost = " << memcostprism
	 << ", scmemcost = " << memcostscprism
	 << ", cpucost = " << cpucostprism << endl;
    cout << "HEX: hp-memcost = " << memcosthex
	 << ", scmemcost = " << memcostschex
	 << ", cpucost = " << cpucosthex << endl;
#endif
  }



  void ReorderPoints (Mesh & mesh, ARRAY<HPRefElement> & hpelements)
  {
    ARRAY<int, 1> map (mesh.GetNP());
    
    for (int i = 1; i <= mesh.GetNP(); i++)
      map[i] = i;

    int nwrong(0), nright(0);
    for (int k = 0; k < 5; k++)
      {
        nwrong = nright = 0;
        for (int i = 0; i < hpelements.Size(); i++)
          {
            const HPRefElement & hpel = hpelements[i];
            
            if (Get_HPRef_Struct (hpel.type) -> geom == HP_PRISM)
              {
                int minbot = 0, mintop = 0;
                for (int j = 0; j < 3; j++)
                  {
                    if (map[hpel.pnums[j]] < map[hpel.pnums[minbot]]) minbot = j;
                    if (map[hpel.pnums[j+3]] < map[hpel.pnums[mintop+3]]) mintop = j;
                  }
                if (minbot != mintop) 
                  nwrong++;
                else
                  nright++;
                
                if (minbot != mintop)
                  {
                    if (map[hpel.pnums[minbot]] < map[hpel.pnums[mintop+3]])
                      swap (map[hpel.pnums[3+minbot]], map[hpel.pnums[3+mintop]]);
                    else
                      swap (map[hpel.pnums[minbot]], map[hpel.pnums[mintop]]);
                  }
              }
          }
        // cout << nwrong << " wrong prisms, " << nright << " right prisms" << endl;
      }

    cout << nwrong << " wrong prisms, " << nright << " right prisms" << endl;


    ARRAY<MeshPoint, 1> hpts(mesh.GetNP());

    for (int i = 1; i <= mesh.GetNP(); i++)
      hpts[map[i]] = mesh.Point(i);

    for (int i = 1; i <= mesh.GetNP(); i++)
      mesh.Point(i) = hpts[i];

    for (int i = 0; i < hpelements.Size(); i++)
      {
        HPRefElement & hpel = hpelements[i];
        for (int j = 0; j < hpel.np; j++)
          hpel.pnums[j] = map[hpel.pnums[j]];
      }
  }



  /* ***************************** HPRefinement ********************************** */

  void HPRefinement (Mesh & mesh, Refinement * ref, int levels, double fac1, bool setorders, bool reflevels)
  {
    PrintMessage (1, "HP Refinement called, levels = ", levels);

 
    NgLock mem_lock (mem_mutex,1);

    mesh.coarsemesh = new Mesh; 
    *mesh.coarsemesh = mesh;
    
#ifdef CURVEDELEMS_NEW
    const_cast<CurvedElements&> (mesh.coarsemesh->GetCurvedElements() ).
      BuildCurvedElements (ref, mesh.GetCurvedElements().GetOrder());
#endif


    delete mesh.hpelements;
    mesh.hpelements = new ARRAY<HPRefElement>;
        
    ARRAY<HPRefElement> & hpelements = *mesh.hpelements; 
        
    InitHPElements(mesh,hpelements); 
    
    ARRAY<int> nplevel;
    nplevel.Append (mesh.GetNP());
    
    int act_ref=1;
    bool sing = ClassifyHPElements(mesh,hpelements, act_ref, levels); 

    sing = true; // iterate at least once
    while(sing) 
      {
	cout << " Start new hp-refinement: step " <<  act_ref  << endl; 
		
	DoRefinement (mesh, hpelements, ref, fac1); 
	DoRefineDummies (mesh, hpelements, ref);
	
	nplevel.Append (mesh.GetNP());
	CalcStatistics (hpelements);
	
	SubdivideDegeneratedHexes (mesh, hpelements,fac1);

        ReorderPoints (mesh, hpelements);

	mesh.ClearSegments();
	mesh.ClearSurfaceElements();
  	mesh.ClearVolumeElements();

	for (int i = 0; i < hpelements.Size(); i++)
	  {
	    HPRefElement & hpel = hpelements[i];
	    if (Get_HPRef_Struct (hpel.type))
	      switch (Get_HPRef_Struct (hpel.type) -> geom)
		{
		case HP_SEGM:
		  {
		    Segment seg;
		    seg.p1 = hpel.pnums[0];
		    seg.p2 = hpel.pnums[1];
		    // NOTE: only for less than 10000 elements (HACK) !!!
		    seg.edgenr = hpel.index % 10000;
		    seg.si     = hpel.index / 10000;

                    /*
                    seg.epgeominfo[0].dist = hpel.param[0][0]; // he: war hpel.param[0][0]
                    seg.epgeominfo[1].dist = hpel.param[1][0]; // he: war hpel.param[1][0]
                    */
                    
                    const Segment & coarseseg = mesh.coarsemesh->LineSegment(hpel.coarse_elnr+1);
                    double d1 = coarseseg.epgeominfo[0].dist;
                    double d2 = coarseseg.epgeominfo[1].dist;

                    // seg.epgeominfo[0].dist = hpel.param[0][0]; // he: war hpel.param[0][0]
                    // seg.epgeominfo[1].dist = hpel.param[1][0]; // he: war hpel.param[1][0]

                    seg.epgeominfo[0].dist = d1 + hpel.param[0][0] * (d2-d1); // JS, June 08
                    seg.epgeominfo[1].dist = d1 + hpel.param[1][0] * (d2-d1); 


		    seg.epgeominfo[0].edgenr = seg.edgenr;
		    seg.epgeominfo[1].edgenr = seg.edgenr;
                    seg.domin = hpel.domin; seg.domout=hpel.domout; // he: needed for segments!
		    seg.hp_elnr = i;
		    seg.singedge_left = hpel.singedge_left; 
		    seg.singedge_right = hpel.singedge_right; 
		    mesh.AddSegment (seg); 
		    break;
		  }
		  
		case HP_TRIG: 
		case HP_QUAD: 
		  { 
		    Element2d el(hpel.np); 
		    for(int j=0;j<hpel.np;j++) 
		      el.PNum(j+1) = hpel.pnums[j]; 
		    el.hp_elnr = i; 
		    el.SetIndex(hpel.index);
		    if(setorders)
		      el.SetOrder(act_ref+1,act_ref+1,0); 
		    mesh.AddSurfaceElement(el);
		    break; 
		  } 
		case HP_HEX:
		case HP_TET:
		case HP_PRISM:
		case HP_PYRAMID:
		  { 
		    Element el(hpel.np); 
		    for(int j=0;j<hpel.np;j++) 
		      el.PNum(j+1) = hpel.pnums[j]; 
		    el.SetIndex(hpel.index); 
		    el.hp_elnr = i; 
		    if(setorders)
		      el.SetOrder(act_ref+1,act_ref+1,act_ref+1);
		    mesh.AddVolumeElement(el); 
		    break;
		  } 
		      
		default:
		  PrintSysError ("hpref, backconversion failed for element ", 
				 int(Get_HPRef_Struct (hpel.type) -> geom));
		}
	  }
	cout << " Start with Update Topology " << endl; 
	mesh.UpdateTopology();
	cout << " Mesh Update Topology done " << endl; 

	act_ref++; 
	
	sing = ClassifyHPElements(mesh,hpelements, act_ref, levels); 
      }

    cout << " HP-Refinement done with " << --act_ref << " refinement steps." << endl; 

    if(act_ref>=1)
      { 
	for(ElementIndex i=0;i<mesh.GetNE(); i++) 
	  { 
	    Element el = mesh[i] ;
	    HPRefElement & hpel = hpelements[mesh[i].hp_elnr];
	    const ELEMENT_EDGE * edges = MeshTopology::GetEdges (mesh[i].GetType());
	    double dist[3] = {0,0,0}; 
	    int ord_dir[3] = {0,0,0}; 
	    int edge_dir[12] = {0,0,0,0,0,0,0,0,0,0,0,0}; 
	    int ned = 4; 
	    
	    switch (mesh[i].GetType())
	      {
	      case TET: 
		/* cout << " TET " ; 
		for(int k=0;k<4;k++) cout << el[k] << "\t" ; 
		cout << endl; */ 
		break; 
	      case PRISM:
		/* cout << " PRISM " ; 
		for(int k=0;k<6;k++) cout << el[k] << "\t" ; 
		cout << endl;  */ 
		for(int l=6;l<9;l++) edge_dir[l] = 2; 
		ord_dir[2] = 2; 
		ned = 9; 
		break; 
	      case HEX: 
		/* cout << " HEX " ; 
		for(int k=0;k<8;k++) cout << el[k] << "\t" ; 
		cout << endl; */
		for(int l=8;l<12; l++) edge_dir[l] = 2; 
		edge_dir[2] = edge_dir[3] = edge_dir[6] = edge_dir[7] = 1;
		ord_dir[1] = 1; 
		ord_dir[2] = 2; 
		ned = 12; 
		break;  
	      case PYRAMID: 
		/*	cout << " PYRAMID " ; 
		for(int k=0;k<5;k++) cout << el[k] << "\t" ; 
		cout << endl; */ 
		for(int l=4;l<8;l++) edge_dir[l] = 2; 
		edge_dir[2] = edge_dir[3] = 1; 
		ord_dir[1] = 1; 
		ord_dir[2] = 2; 
		ned = 8;  
		break; 
	      }
	
	    for (int j=0;j<ned;j++) 
	      { 
			
		Vec<3> v(hpel.param[edges[j][0]-1][0]-hpel.param[edges[j][1]-1][0],
			    hpel.param[edges[j][0]-1][1]-hpel.param[edges[j][1]-1][1],
			    hpel.param[edges[j][0]-1][2]-hpel.param[edges[j][1]-1][2]);
		dist[edge_dir[j]] = max(v.Length(),dist[edge_dir[j]]);
	      }
	    
	    int refi[3];  
	    for(int j=0;j<3;j++) 
	      refi[j] = int(max(double(floor(log(dist[ord_dir[j]]/sqrt(2.))/log(fac1))),0.)); 	
	    
	    // cout << " ref " << refi[0] << "\t" << refi[1] << "\t" << refi[2] << endl; 
	    // cout << " order " << act_ref +1 - refi[0] << "\t" << act_ref +1 - refi[1] << "\t" << act_ref +1 - refi[2] << endl; 
	   	      
	    if(setorders)
	      mesh[i].SetOrder(act_ref+1-refi[0],act_ref+1-refi[1],act_ref+1-refi[2]); 
	  }
	for(SurfaceElementIndex i=0;i<mesh.GetNSE(); i++) 
	  { 
	    Element2d el = mesh[i] ;
	    HPRefElement & hpel = hpelements[mesh[i].hp_elnr];
	    const ELEMENT_EDGE * edges = MeshTopology::GetEdges (mesh[i].GetType());
	    double dist[3] = {0,0,0}; 
	    int ord_dir[3] = {0,0,0}; 
	    int  edge_dir[4] = {0,0,0,0} ; 
	    int ned = 3; 
	   
	    if(mesh[i].GetType() == QUAD)
	      {
		/*	cout << " QUAD " ; 
		for(int k=0;k<4;k++) cout << el[k] << "\t" ; 
		cout << endl; 	*/ 
 
		edge_dir[2] = edge_dir[3] = 1; 
		ord_dir[1] = 1; 
		ned = 4; 
	      }
	    /*  else 
	      { 
		cout << " TRIG " ; 
		for(int k=0;k<3;k++) cout << el[k] << "\t" ; 
		cout << endl; 
		} */ 
	    
	    for (int j=0;j<ned;j++) 
	      { 
		Vec<3> v(hpel.param[edges[j][0]-1][0]-hpel.param[edges[j][1]-1][0],
			    hpel.param[edges[j][0]-1][1]-hpel.param[edges[j][1]-1][1],
			    hpel.param[edges[j][0]-1][2]-hpel.param[edges[j][1]-1][2]);
		dist[edge_dir[j]] = max(v.Length(),dist[edge_dir[j]]);
	      }
	    
	    int refi[3]; 
	    for(int j=0;j<3;j++) 
	      refi[j] = int(max(double(floor(log(dist[ord_dir[j]]/sqrt(2.))/log(fac1))),0.)); 	
	    
	    if(setorders)
	      mesh[i].SetOrder(act_ref+1-refi[0],act_ref+1-refi[1],act_ref+1-refi[2]); 

	      // cout << " ref " << refi[0] << "\t" << refi[1] << endl; 
	      // cout << " order " << act_ref +1 - refi[0] << "\t" << act_ref +1 - refi[1] << endl; 
	  }
      }
  }

bool CheckSingularities(Mesh & mesh, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		       BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
			INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint, int & levels, int & act_ref)
{ 
  bool sing=0; 
  if (mesh.GetDimension() == 3)
      {
	/*
	// check, if point has as least 3 different surfs:

	ARRAY<INDEX_3, PointIndex::BASE> surfonpoint(mesh.GetNP());
  	surfonpoint = INDEX_3(0,0,0);

	for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	  {
	    const Element2d & el = mesh[sei];
	    int ind = el.GetIndex();
	    for (int j = 0; j < el.GetNP(); j++)
	      {
		INDEX_3 & i3 = surfonpoint[el[j]];
		if (ind != i3.I1() && ind != i3.I2() && ind != i3.I3())
		  {
		    i3.I1() = i3.I2();
		    i3.I2() = i3.I3();
		    i3.I3() = ind;
		  }
	      }
	  }
	for (int i = 1; i <= mesh.GetNP(); i++)
	  if (surfonpoint.Get(i).I1())
	    cornerpoint.Set(i);
	*/
	cornerpoint.Clear();
	
	for (int i = 1; i <= mesh.GetNP(); i++)
	  {
	    if (mesh.Point(i).Singularity() * levels >= act_ref)
	      {
		cornerpoint.Set(i);
		sing = 1; 
	      } 
	  }
	cout << endl; 

	for (int i = 1; i <= mesh.GetNSeg(); i++)
	  if (mesh.LineSegment(i).singedge_left * levels >= act_ref)
	    {
	      INDEX_2 i2 (mesh.LineSegment(i).p1, 
			  mesh.LineSegment(i).p2);

	      /*
		// before
	      edges.Set (i2, 1);
	      i2.Sort();   
	      INDEX_2 i2s(i2.I2(), i2.I1());
	      edges.Set (i2s, 1);
	      */

	      edges.Set (i2, 1);
	      INDEX_2 i2s(i2.I2(), i2.I1());
	      edges.Set (i2s, 1);


	      edgepoint.Set (i2.I1());
	      edgepoint.Set (i2.I2());
	      sing = 1; 
	    }

	// if 2 adjacent edges of an element are singular, the 
	// commen point must be a singular point
	for (int i = 1; i <= mesh.GetNE(); i++)
	  {
	    const Element & el = mesh.VolumeElement(i);
	    const ELEMENT_EDGE * eledges = MeshTopology::GetEdges (el.GetType());
	    int nedges = MeshTopology::GetNEdges (el.GetType());
	    for (int j = 0; j < nedges; j++)
	      for (int k = 0; k < nedges; k++)
		if (j != k)
		  {
		    INDEX_2 ej(el.PNum(eledges[j][0]), el.PNum(eledges[j][1]));
		    ej.Sort();
		    INDEX_2 ek(el.PNum(eledges[k][0]), el.PNum(eledges[k][1]));
		    ek.Sort();
		    if (edges.Used(ej) && edges.Used(ek))
		      {
			if (ej.I1() == ek.I1()) cornerpoint.Set (ek.I1());
			if (ej.I1() == ek.I2()) cornerpoint.Set (ek.I2());
			if (ej.I2() == ek.I1()) cornerpoint.Set (ek.I1());
			if (ej.I2() == ek.I2()) cornerpoint.Set (ek.I2());
		      }
		  }
	  }

	edgepoint.Or (cornerpoint);
	(*testout) << "cornerpoint = " << endl << cornerpoint << endl;
	(*testout) << "edgepoint = " << endl << edgepoint << endl;

	facepoint = 0;
	for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	  {
	    const Element2d & el = mesh[sei];
	    const FaceDescriptor & fd = mesh.GetFaceDescriptor (el.GetIndex());
	  
	    int domnr = 0;
	    if (fd.domin_singular * levels < act_ref && fd.domout_singular * levels < act_ref) 
	      { domnr=0;  continue;}
	    
	    if (fd.domin_singular * levels >= act_ref) 
	      {
		domnr = fd.DomainIn();
		sing = 1;
	      }
	    if (fd.domout_singular * levels >= act_ref)
	      {
		domnr = fd.DomainOut();
		sing = 1; 
	      } 
	    if (fd.domin_singular * levels >= act_ref 
		&& fd.domout_singular * levels >= act_ref) 
	      {
		domnr = -1;
		sing = 1;
	      } 
  
	    INDEX_3 i3;
	    if (el.GetNP() == 3) 
	      i3 = INDEX_3::Sort (el[0], el[1], el[2]);
	    else
	      {
		INDEX_4 i4 (el[0], el[1], el[2], el[3]);
		i4.Sort();
		i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3());
	      }
	    faces.Set (i3, domnr);
	
	    for (int j = 0; j < el.GetNP(); j++)
	      {
		face_edges.Set (INDEX_2::Sort (el[j], el[(j+1)%el.GetNP()]), domnr);
	
		surf_edges.Set (INDEX_2::Sort (el[j], el[(j+1)%el.GetNP()]), fd.SurfNr()+1);
		
		facepoint[el[j]] = domnr;
	      }
	   
	  }
	(*testout) << "singular faces = " << faces << endl;
	(*testout) << "singular faces_edges = " << face_edges << endl;
      }
    else
      {
	// 2D case

	// check, if point has as least 3 different surfs:
	ARRAY<INDEX_3, PointIndex::BASE> surfonpoint(mesh.GetNP());

	for (int i = 1; i <= mesh.GetNP(); i++)
	  surfonpoint.Elem(i) = INDEX_3(0,0,0);
      
	for (int i = 1; i <= mesh.GetNSeg(); i++)
	  {
	    const Segment & seg = mesh.LineSegment(i);
	    int ind = seg.edgenr;

	   
		if (seg.singedge_left * levels >= act_ref)
		  {
		    INDEX_2 i2 (mesh.LineSegment(i).p1, 
				mesh.LineSegment(i).p2);
		    edges.Set(i2,1); 
		    edgepoint.Set(i2.I1());
		    edgepoint.Set(i2.I2());
		    *testout << " singleft " << endl;  
		    *testout << " mesh.LineSegment(i).domout " << mesh.LineSegment(i).domout << endl;      
		    *testout << " mesh.LineSegment(i).domin " << mesh.LineSegment(i).domin << endl;      
		    edgepoint_dom.Set (INDEX_2(mesh.LineSegment(i).domin, i2.I1()), 1);
		    edgepoint_dom.Set (INDEX_2(mesh.LineSegment(i).domin, i2.I2()), 1);
		    sing = 1; 
		    
		  }
		
		  if (seg.singedge_right * levels >= act_ref)
		    {
		      INDEX_2 i2 (mesh.LineSegment(i).p2, 
				  mesh.LineSegment(i).p1);  
		      edges.Set (i2, 1);
		      edgepoint.Set(i2.I1());
		      edgepoint.Set(i2.I2());

		      *testout << " singright " << endl;  
		      *testout << " mesh.LineSegment(i).domout " << mesh.LineSegment(i).domout << endl;      
		      *testout << " mesh.LineSegment(i).domin " << mesh.LineSegment(i).domin << endl;      
		      
		      edgepoint_dom.Set (INDEX_2(mesh.LineSegment(i).domout, i2.I1()), 1);
		      edgepoint_dom.Set (INDEX_2(mesh.LineSegment(i).domout, i2.I2()), 1);
		      sing = 1;
		    }
	
	    // (*testout) << "seg = " << ind << ", " << seg.p1 << "-" << seg.p2 << endl;


	    if (seg.singedge_left * levels >= act_ref
		|| seg.singedge_right* levels >= act_ref)
	      {
		for (int j = 0; j < 2; j++)
		  {
		    int pi = (j == 0) ? seg.p1 : seg.p2;
		    INDEX_3 & i3 = surfonpoint.Elem(pi);
		    if (ind != i3.I1() &&
			ind != i3.I2())
		      {
			i3.I1() = i3.I2();
			i3.I2() = ind;
		      }
		  }
	      }
	  }


	for (int i = 1; i <= mesh.GetNP(); i++)
	  {
	    // mark points for refinement that are in corners between two anisotropic edges 
	    if (surfonpoint.Get(i).I1())
	      {
		cornerpoint.Set(i);
		edgepoint.Set(i);
	      }
	
	    // mark points for refinement that are explicity specified in input file
	    if (mesh.Point(i).Singularity()*levels >= act_ref)
	      {
		cornerpoint.Set(i);
		edgepoint.Set(i);
		sing =  1; 
	      }
	  }

	edgepoint.Or (cornerpoint);

	(*testout) << "2d sing edges: " << endl << edges << endl;
	(*testout) << "2d cornerpoints: " << endl << cornerpoint << endl
		   << "2d edgepoints: " << endl << edgepoint << endl;
	
	facepoint = 0;
      }

    if (!sing)
      {
	cout << "PrepareElements no more to do for actual refinement " << act_ref << endl; 
	return(sing);
      } 
    return(sing); 
}



  bool ClassifyHPElements (Mesh & mesh, ARRAY<HPRefElement> & elements, int & act_ref, int & levels)
  {
    
    INDEX_2_HASHTABLE<int> edges(mesh.GetNSeg()+1);
    BitArray edgepoint(mesh.GetNP());
    INDEX_2_HASHTABLE<int> edgepoint_dom(mesh.GetNSeg()+1);

    edgepoint.Clear();
    BitArray cornerpoint(mesh.GetNP());
    cornerpoint.Clear();

    // value = nr > 0 ... refine elements in domain nr
    // value = -1   ..... refine elements in any domain
    INDEX_3_HASHTABLE<int> faces(mesh.GetNSE()+1);
    INDEX_2_HASHTABLE<int> face_edges(mesh.GetNSE()+1);
    INDEX_2_HASHTABLE<int> surf_edges(mesh.GetNSE()+1);
    ARRAY<int, PointIndex::BASE> facepoint(mesh.GetNP());

    bool sing = CheckSingularities(mesh, edges, edgepoint_dom, 
			      cornerpoint, edgepoint, faces, face_edges, 
			      surf_edges, facepoint, levels, act_ref); 
  
    if(sing==0) return(sing); 

    int cnt_undef = 0, cnt_nonimplement = 0;
    ARRAY<int> misses(10000);
    misses = 0;

    (*testout) << "edgepoint_dom = " << endl << edgepoint_dom << endl;

    for( int i = 0; i<elements.Size(); i++) 
      {
	// *testout << "classify element " << i << endl;

	HPRefElement & hpel = elements[i]; 
	HPRef_Struct * hprs = Get_HPRef_Struct (hpel.type);
	HPRefElement old_el = elements[i]; 
	int dd=3; 


	if(act_ref !=1 && (hpel.type == HP_HEX || hpel.type == HP_PRISM || hpel.type == HP_TET 
			   || hpel.type == HP_PYRAMID || hpel.type == HP_QUAD || hpel.type == HP_TRIG || hpel.type == HP_SEGM)) 
	  continue; 
	
	sing = 1; 
	switch (hprs->geom)
	  {
	  case HP_TET:
	    {
	      hpel.type = ClassifyTet(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces,face_edges, surf_edges, facepoint); 
	      break;
	    }
	  case HP_PRISM:
	    {
	      hpel.type = ClassifyPrism(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces,
					face_edges, surf_edges, facepoint); 	    	    
	 
	    
	      break;
	    }
	  case HP_HEX:
	    { 
	      hpel.type = hpel.type = ClassifyHex(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces,
						  face_edges, surf_edges, facepoint); 	    	    
	      break; 
	    } 
	  case HP_TRIG: 
	    { 
	      int dim = mesh.GetDimension(); 
	      const FaceDescriptor & fd = mesh.GetFaceDescriptor (hpel.GetIndex());	
	      
	      hpel.type = ClassifyTrig(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, 
				       faces, face_edges, surf_edges, facepoint, dim, fd);    
	     
	      dd = 2; 
	      break; 
	    } 
	  case HP_QUAD: 
	    { 
	      int dim = mesh.GetDimension(); 
	      const FaceDescriptor & fd = mesh.GetFaceDescriptor (hpel.GetIndex());	
	      hpel.type = ClassifyQuad(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, 
				  faces, face_edges, surf_edges, facepoint, dim, fd);    

	      dd = 2; 
	      break; 
	    }
	  case HP_SEGM: 
	    {
	      hpel.type = ClassifySegm(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, 
				       faces, face_edges, surf_edges, facepoint);    
	      dd = 1; 
	      break; 
	    }
	  case HP_PYRAMID: 
	    {
	      hpel.type = ClassifyPyramid(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces,
						  face_edges, surf_edges, facepoint); 	    	    
	      
	      cout << " ** Pyramid classified  " << hpel.type << endl; 
	      break; 
	    }
	  default:
	    {
	      cout << "illegal element type for hp-prepare elements " << hpel.type << endl;
	      throw NgException ("hprefinement.cpp: don't know how to set parameters");
	    }
	  }
	    
	if(hpel.type == HP_NONE) 
	  cnt_undef++; 

	//else 
	//cout << "elem " << i << " classified type " << hpel.type << endl; 

	
	
	if (!Get_HPRef_Struct (hpel.type)) 
	  {
	    (*testout) << "hp-element-type " << hpel.type << " not implemented   " << endl;
	    (*testout) << " elType " << hprs->geom << endl; 
 (cout) << " elType " << hprs->geom << endl;        
	    cnt_nonimplement++;
	    misses[hpel.type]++;
	  }
	
  
	for(int j=0; j<hpel.np; j++)
	  {
	    for( int k=0; k<hpel.np; k++) 
	      if(hpel[j] == old_el.pnums[k]) 
		{ 
		  for(int l=0;l<dd;l++) 
		    hpel.param[j][l] = old_el.param[k][l];
		  break;
		}
	  } 

      }
    
    
    cout << "undefined elements update classification: " << cnt_undef << endl;
    cout << "non-implemented in update classification: " << cnt_nonimplement << endl;

    for (int i = 0; i < misses.Size(); i++)
      if (misses[i])
	cout << " in update classification missing case " << i << " occured " << misses[i] << " times" << endl;

    return(sing); 
  }
}
  
