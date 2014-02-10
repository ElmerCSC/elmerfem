HPREF_ELEMENT_TYPE ClassifyTet(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
		 INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint)
{
  int ep1(0), ep2(0), ep3(0), ep4(0), cp1(0), cp2(0), cp3(0), cp4(0), fp1, fp2, fp3, fp4;
  int isedge1(0), isedge2(0), isedge3(0), isedge4(0), isedge5(0), isedge6(0);
  int isfedge1, isfedge2, isfedge3, isfedge4, isfedge5, isfedge6;
  int isface1(0), isface2(0), isface3(0), isface4(0);

  HPREF_ELEMENT_TYPE type = HP_NONE; 
  

  int debug = 0;
  for (int j = 0;j < 4; j++)
    {
      if (el.pnums[j] == 444) debug++;
      if (el.pnums[j] == 115) debug++;
      if (el.pnums[j] == 382) debug++;
      if (el.pnums[j] == 281) debug++;
    }
  if (debug < 4) debug = 0;
  


  for (int j = 0; j < 4; j++)
    for (int k = 0; k < 4; k++)
      {
	if (j == k) continue;
	if (type) break;
	
	int pi3 = 0;
	while (pi3 == j || pi3 == k) pi3++;
	int pi4 = 6 - j - k - pi3;
	
	// preserve orientation
	int sort[4];
	sort[0] = j; sort[1] = k; sort[2] = pi3; sort[3] = pi4;
	int cnt = 0;
	for (int jj = 0; jj < 4; jj++)
	  for (int kk = 0; kk < 3; kk++)
	    if (sort[kk] > sort[kk+1])
	      {
		cnt++;
		Swap (sort[kk], sort[kk+1]); 
	      }
	if (cnt % 2 == 1) Swap (pi3, pi4);
	
	ep1 = edgepoint.Test (el.pnums[j]);
	ep2 = edgepoint.Test (el.pnums[k]);
	ep3 = edgepoint.Test (el.pnums[pi3]);
	ep4 = edgepoint.Test (el.pnums[pi4]);
	
	cp1 = cornerpoint.Test (el.pnums[j]);
	cp2 = cornerpoint.Test (el.pnums[k]);
	cp3 = cornerpoint.Test (el.pnums[pi3]);
	cp4 = cornerpoint.Test (el.pnums[pi4]);
	
	isedge1 = edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[k]));
	isedge2 = edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[pi3]));
	isedge3 = edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[pi4]));
	isedge4 = edges.Used (INDEX_2::Sort (el.pnums[k], el.pnums[pi3]));
	isedge5 = edges.Used (INDEX_2::Sort (el.pnums[k], el.pnums[pi4]));
	isedge6 = edges.Used (INDEX_2::Sort (el.pnums[pi3], el.pnums[pi4]));
	
	if (debug)
	  {
	    cout << "debug" << endl;
	    *testout  << "debug" << endl;
	    *testout << "ep = " << ep1 << ep2 << ep3 << ep4 << endl;
	    *testout << "cp = " << cp1 << cp2 << cp3 << cp4 << endl;
	    *testout << "edge = " << isedge1 << isedge2 << isedge3 << isedge4 << isedge5 << isedge6 << endl;
	  }


	isface1 = isface2 = isface3 = isface4 = 0;
	for (int l = 0; l < 4; l++)
	  {
	    INDEX_3 i3(0,0,0);
	    switch (l)
	      {
		case 0: i3.I1() = el.pnums[k]; i3.I1() = el.pnums[pi3]; i3.I1() = el.pnums[pi4]; break;
		case 1: i3.I1() = el.pnums[j]; i3.I1() = el.pnums[pi3]; i3.I1() = el.pnums[pi4]; break;
		case 2: i3.I1() = el.pnums[j]; i3.I1() = el.pnums[k]; i3.I1() = el.pnums[pi4]; break;
		case 3: i3.I1() = el.pnums[j]; i3.I1() = el.pnums[k]; i3.I1() = el.pnums[pi3]; break;
	      }
	    i3.Sort();
	    if (faces.Used (i3))
	      {
		int domnr = faces.Get(i3);
		if (domnr == -1 || domnr == el.GetIndex())
		  {
		    switch (l)
		      {
		      case 0: isface1 = 1; break;
		      case 1: isface2 = 1; break;
		      case 2: isface3 = 1; break;
		      case 3: isface4 = 1; break;
		      }
		  }
	      }
	  }
	/*
	  isface1 = faces.Used (INDEX_3::Sort (el.pnums[k], el.pnums[pi3], el.pnums[pi4]));
	  isface2 = faces.Used (INDEX_3::Sort (el.pnums[j], el.pnums[pi3], el.pnums[pi4]));
	  isface3 = faces.Used (INDEX_3::Sort (el.pnums[j], el.pnums[k], el.pnums[pi4]));
	  isface4 = faces.Used (INDEX_3::Sort (el.pnums[j], el.pnums[k], el.pnums[pi3]));
	*/
	
	isfedge1 = isfedge2 = isfedge3 = isfedge4 = isfedge5 = isfedge6 = 0;
	for (int l = 0; l < 6; l++)
	  {
	    INDEX_2 i2(0,0);
	    switch (l)
	      {
		case 0: i2.I1() = el.pnums[j]; i2.I2() = el[k]; break;
		case 1: i2.I1() = el.pnums[j]; i2.I2() = el.pnums[pi3]; break;
		case 2: i2.I1() = el.pnums[j]; i2.I2() = el.pnums[pi4]; break;
		case 3: i2.I1() = el.pnums[k]; i2.I2() = el.pnums[pi3]; break;
		case 4: i2.I1() = el.pnums[k]; i2.I2() = el.pnums[pi4]; break;
		case 5: i2.I1() = el.pnums[pi3]; i2.I2() = el.pnums[pi4]; break;
	      }
	    i2.Sort();
	    if (face_edges.Used (i2))
	      {
		int domnr = face_edges.Get(i2);
		if (domnr == -1 || domnr == el.GetIndex())
		  {
		    switch (l)
		      {
		      case 0: isfedge1 = 1; break;
		      case 1: isfedge2 = 1; break;
		      case 2: isfedge3 = 1; break;
		      case 3: isfedge4 = 1; break;
		      case 4: isfedge5 = 1; break;
		      case 5: isfedge6 = 1; break;
		      }
		  }
	      }
	  }
	/*
	  isfedge1 = face_edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[k]));
	  isfedge2 = face_edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[pi3]));
	  isfedge3 = face_edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[pi4]));
	  isfedge4 = face_edges.Used (INDEX_2::Sort (el.pnums[k], el.pnums[pi3]));
	  isfedge5 = face_edges.Used (INDEX_2::Sort (el.pnums[k], el.pnums[pi4]));
	  isfedge6 = face_edges.Used (INDEX_2::Sort (el.pnums[pi3], el.pnums[pi4]));
	*/
	
	fp1 = fp2 = fp3 = fp4 = 0;
	for (int l = 0; l < 4; l++)
	  {
	    int pti(0);
	    switch (l)
	      {
	      case 0: pti = el.pnums[j]; break;
	      case 1: pti = el.pnums[k]; break;
	      case 2: pti = el.pnums[pi3]; break;
	      case 3: pti = el.pnums[pi4]; break;
	      }
	    int domnr = facepoint[pti];
	    if (domnr == -1 || domnr == el.GetIndex())
	      {
		switch (l)
		  {
		  case 0: fp1 = 1; break;
		  case 1: fp2 = 1; break;
		  case 2: fp3 = 1; break;
		  case 3: fp4 = 1; break;
		  }
	      }
	  }
	
	/*
	  fp1 = facepoint[el.pnums[j]] != 0;
	  fp2 = facepoint[el.pnums[k]] != 0;
	  fp3 = facepoint[el.pnums[pi3]] != 0;
	  fp4 = facepoint[el.pnums[pi4]] != 0;
	*/
	
	
	switch (isface1+isface2+isface3+isface4)
	  {
	  case 0:
	    {
	      isedge1 |= isfedge1;
	      isedge2 |= isfedge2;
	      isedge3 |= isfedge3;
	      isedge4 |= isfedge4;
	      isedge5 |= isfedge5;
	      isedge6 |= isfedge6;
	      
	      ep1 |= fp1;
	      ep2 |= fp2;
	      ep3 |= fp3;
	      ep4 |= fp4;
	      
	      switch (isedge1+isedge2+isedge3+isedge4+isedge5+isedge6)
		{
		case 0:
		  {		
		    if (!ep1 && !ep2 && !ep3 && !ep4)
		      type = HP_TET;
				
		    if (ep1 && !ep2 && !ep3 && !ep4)
		      type = HP_TET_0E_1V;
		    
		    if (ep1 && ep2 && !ep3 && !ep4)
		      type = HP_TET_0E_2V;
		    
		    if (ep1 && ep2 && ep3 && !ep4)
		      type = HP_TET_0E_3V;
		    
		    if (ep1 && ep2 && ep3 && ep4)
		      type = HP_TET_0E_4V;
		    
		    break;
		  }
		  
		case 1:
		  {
		    if (!isedge1) break;
		    
		    if (!cp1 && !cp2 && !ep3 && !ep4)
		      type = HP_TET_1E_0V;
		    
		    if (cp1 && !cp2 && !ep3 && !ep4)
		      type = HP_TET_1E_1VA;
		    
		    if (!cp1 && !cp2 && !ep3 && ep4)
		      type = HP_TET_1E_1VB;
		    
		    if (cp1 && cp2 && !ep3 && !ep4)
		      type = HP_TET_1E_2VA;
		    
		    if (cp1 && !cp2 && ep3 && !ep4)
		      type = HP_TET_1E_2VB;
		    
		    if (cp1 && !cp2 && !ep3 && ep4)
		      type = HP_TET_1E_2VC;
		    
		    if (!cp1 && !cp2 && ep3 && ep4)
		      type = HP_TET_1E_2VD;
		    
		    if (cp1 && cp2 && ep3 && !ep4)
		      type = HP_TET_1E_3VA;
		    
		    if (cp1 && !cp2 && ep3 && ep4)
		      type = HP_TET_1E_3VB;
		    
		    if (cp1 && cp2 && ep3 && ep4)
		      type = HP_TET_1E_4V;
		    
		    break;
		  }
		case 2:
		  {
		    if (isedge1 && isedge2)
		      {
			if (!cp2 && !cp3 && !ep4)
			  type = HP_TET_2EA_0V;
			
			if (cp2 && !cp3 && !ep4)
			  type = HP_TET_2EA_1VA;
			if (!cp2 && cp3 && !ep4)
			  type = HP_TET_2EA_1VB;
			
			if (!cp2 && !cp3 && ep4)
			  type = HP_TET_2EA_1VC;
			
			if (cp2 && cp3 && !ep4)
			  type = HP_TET_2EA_2VA;
			if (cp2 && !cp3 && ep4)
			  type = HP_TET_2EA_2VB;
			if (!cp2 && cp3 && ep4)
			  type = HP_TET_2EA_2VC;
			
			if (cp2 && cp3 && ep4)
			  type = HP_TET_2EA_3V;
		      }
		    if (isedge1 && isedge6)
		      {
			if (!cp1 && !cp2 && !cp3 && !cp4)
			  type = HP_TET_2EB_0V;
			if (cp1 && !cp2 && !cp3 && !cp4)
			  type = HP_TET_2EB_1V;
			if (cp1 && cp2 && !cp3 && !cp4)
			  type = HP_TET_2EB_2VA;
			if (cp1 && !cp2 && cp3 && !cp4)
			  type = HP_TET_2EB_2VB;
			if (cp1 && !cp2 && !cp3 && cp4)
			  type = HP_TET_2EB_2VC;
			if (cp1 && cp2 && cp3 && !cp4)
			  type = HP_TET_2EB_3V;
			if (cp1 && cp2 && cp3 && cp4)
			  type = HP_TET_2EB_4V;
		      }
		    break;
		  }
		case 3:
		  {
		    if (isedge1 && isedge2 && isedge3)
		      {
			if (!cp2 && !cp3 && !cp4)
			  type = HP_TET_3EA_0V;
			if (cp2 && !cp3 && !cp4)
			  type = HP_TET_3EA_1V;
			if (cp2 && cp3 && !cp4)
			  type = HP_TET_3EA_2V;
			if (cp2 && cp3 && cp4)
			  type = HP_TET_3EA_3V;
		      }
		    if (isedge1 && isedge3 && isedge4)
		      {
			if (!cp3 && !cp4)
			  type = HP_TET_3EB_0V;
			if (cp3 && !cp4)
				      type = HP_TET_3EB_1V;
			if (cp3 && cp4)
			  type = HP_TET_3EB_2V;
		      }
		    if (isedge1 && isedge2 && isedge5)
		      {
			if (!cp3 && !cp4)
			  type = HP_TET_3EC_0V;
			if (cp3 && !cp4)
			  type = HP_TET_3EC_1V;
			if (cp3 && cp4)
			  type = HP_TET_3EC_2V;
		      }
		    break;
		  }
		}
	      break;
	    }
	    
	    
	    
	  case 1:  // one singular face
	    {
	      if (!isface1) break;
	      
	      switch (isfedge1+isfedge2+isfedge3+isedge4+isedge5+isedge6)
		{
		case 0:
		  {
		    if (!fp1 && !ep2 && !ep3 && !ep4)
		      type = HP_TET_1F_0E_0V;
		    if (fp1 && !ep2 && !ep3 && !ep4)
		      type = HP_TET_1F_0E_1VB;
		    if (!fp1 && ep2 && !ep3 & !ep4)
		      type = HP_TET_1F_0E_1VA;
		    break;
		  }
		case 1:
		  {
		    if (isfedge1)
		      {
			if (!ep1 && !ep3 && !ep4)
			  type = HP_TET_1F_1EA_0V;
		      }
		    if (isedge4) // V1-V3
		      {
			if (!ep1 && !cp2 && !cp3 && !ep4)
			  type = HP_TET_1F_1EB_0V;
		      }
		    break;
		  }
		}
	      break;
	    }
	    
	    
	  case 2:  // two singular faces
	    {
	      if (!isface1 || !isface2) break;
	      
	      switch (isfedge1+isedge2+isedge3+isedge4+isedge5)
		{
		case 0:
		  {
		    if (!ep1 && !ep2 && !cp3 && !cp4)
		      type = HP_TET_2F_0E_0V;
		    break;
		  }
		}
	      break;
	    }
	    
	    
	  }
	
	if (type != HP_NONE)
	  {
	    int pnums[4]; 
	    pnums[0] = el.pnums[j];
	    pnums[1] = el.pnums[k];
	    pnums[2] = el.pnums[pi3];
	    pnums[3] = el.pnums[pi4];
	    for(k=0;k<4;k++) el.pnums[k] = pnums[k]; 
	    break;
	  }
      }
  
  
  if (debug) cout << "type = " << type << endl;

  if (type == HP_NONE)
    {
      //     cnt_undef++;
      (*testout) << "undefined element" << endl
		 << "cp = " << cp1 << cp2 << cp3 << cp4 << endl
		 << "ep = " << ep1 << ep2 << ep3 << ep4 << endl
		 << "isedge = " << isedge1 << isedge2 << isedge3 
		 << isedge4 << isedge5 << isedge6 << endl
		 << "isface = " << isface1 << isface2 << isface3 << isface4 << endl;
      cout << "undefined element !!! " << endl;

      
    }
  return(type); 
}



HPREF_ELEMENT_TYPE ClassifyPrism(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
		 INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint)
{

  HPREF_ELEMENT_TYPE type = HP_NONE;
  
  int p[6];
  for(int m=1;m<=6;m++)
    {
      int point_sing[6]={0,0,0,0,0,0}; 
      int face_sing[5]={0,0,0,0,0};
      int edge_sing[9]={0,0,0,0,0,0,0,0,0}; 
      
      if(m<4)
	{ 
	  p[0]= m; p[1]=m%3+1; p[2]=(m%3+1)%3+1;
	  for(int l=3;l<6;l++) p[l]=p[l-3]+3;  
	}
      else
	{
	  p[0] = m; p[1]=(m%3+1)%3+4; p[2]=m%3+4;
	  for(int l=3;l<6;l++) p[l]=p[l-3]-3; 
	}
      
      for(int j=0;j<6;j++) 
	{ 
	  if(cornerpoint.Test(el.PNum(p[j])))  { point_sing[p[j]-1]=3;}
	  else if(edgepoint.Test(el.PNum(p[j]))) point_sing[p[j]-1]=2;
	  else if (facepoint[el.PNum(p[j])] == -1 || facepoint[el.PNum(p[j])] == el.GetIndex())
	    point_sing[p[j]-1] = 1;  
	}
      
      const ELEMENT_EDGE * eledges = MeshTopology::GetEdges (PRISM);
      for(int k=0;k<9;k++)
	{
	  INDEX_2 i2 = INDEX_2 :: Sort(el.PNum(p[eledges[k][0]-1]),el.PNum(p[eledges[k][1]-1])); 
	  if (edges.Used(i2)) edge_sing[k] = 2;
	  else edge_sing[k] = face_edges.Used(i2);
	}
      
      const ELEMENT_FACE * elfaces  = MeshTopology::GetFaces (PRISM);
      for (int k=0;k<5;k++)
	{
	  INDEX_3 i3; 
	  
	  if(k<2) 
	    i3 = INDEX_3::Sort(el.pnums[p[elfaces[k][0]-1]-1], el.pnums[p[elfaces[k][1]-1]-1], 
			       el.pnums[p[elfaces[k][2]-1]-1]); 
	  else 
	    { 
	      INDEX_4  i4 = INDEX_4(el.pnums[p[elfaces[k][0]-1]-1], el.pnums[p[elfaces[k][1]-1]-1], el.pnums[p[elfaces[k][2]-1]-1],el.pnums[p[elfaces[k][3]-1]-1]); 
	      i4.Sort();
	      i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3()); 
	    }
	  
	  if (faces.Used (i3))
	    {
	      int domnr = faces.Get(i3); 
	      if (domnr == -1 || domnr == el.GetIndex())
		face_sing[k] = 1; 
	      
	    } 
	} 
      if (face_sing[1] > face_sing[0]) {m=m+2; continue;}  
      
      
      //int cp = 0;  
      
      int qfsing = face_sing[2] + face_sing[3] + face_sing[4];
      int tfsing = face_sing[0] + face_sing[1]; 
      int evsing = edge_sing[6] + edge_sing[7] + edge_sing[8];
      int ehsing = edge_sing[0] + edge_sing[1] + edge_sing[2] + edge_sing[3] + edge_sing[4] + edge_sing[5];
      
      if (qfsing + tfsing + evsing + ehsing == 0)  
	{ type = HP_PRISM;  break;}
      
      HPREF_ELEMENT_TYPE types[] = {HP_NONE,HP_NONE,HP_NONE};   
      
      int fb = (1-face_sing[4])* face_sing[3] * (face_sing[2] + face_sing[3]) + 3*face_sing[4]*face_sing[3]*face_sing[2];  
      int sve[3] = {edge_sing[7] , edge_sing[8], edge_sing[6]}; 
      
            
      if(fb!=qfsing) continue; 
      
      
      switch(fb)
	{ 
	case 0: 
	  if (evsing == 0 && ehsing==3*tfsing) 
	    {
	      types[0] = HP_PRISM; 
	      types[1] = HP_PRISM_1FA_0E_0V;   
	      types[2] = HP_PRISM_2FA_0E_0V; 
	    } 
	  if(evsing > 0 &&  sve[0] == evsing) // 1 vertical edge 1-4 
	    { 
	      types[0] = HP_PRISM_SINGEDGE;
	      types[1] = HP_PRISM_1FA_1E_0V;
	      types[2] = HP_PRISM_2FA_1E_0V;   
	    }
	  
	  if(sve[0] > 0 && sve[1] > 0 && sve[2] == 0)
	    {
	      types[0] = HP_PRISM_SINGEDGE_V12;
	      types[1] = HP_PRISM_1FA_2E_0V; 
	      types[2] = HP_PRISM_2FA_2E_0V; 
	    }
	  if(sve[0] > 0 && sve[1] > 0 && sve[2] > 0) 
	    {
	      types[0] = HP_PRISM_3E_0V;
	      types[1] = HP_PRISM_1FA_3E_0V;
	      types[2] = HP_PRISM_2FA_3E_0V;
	      
	      if ( edge_sing[0] > 1 && edge_sing[2] > 1 &&  
		   edge_sing[4] > 1 && edge_sing[5] > 1 && tfsing==0)
		types[0] = HP_PRISM_3E_4EH; 
	    }
	  
	  break;
	case 1:
	  if(sve[0] <= 1 && sve[1] <= 1)  
	    if(sve[2]==0)
	      { 
		types[0] = HP_PRISM_1FB_0E_0V;
		types[1] = HP_PRISM_1FA_1FB_0E_0V;
		types[2] = HP_PRISM_2FA_1FB_0E_0V;  
	      }
	    else
	      { 
		types[0] = HP_PRISM_1FB_1EC_0V;
		types[1] = HP_PRISM_1FA_1FB_1EC_0V;
		types[2] = HP_PRISM_2FA_1FB_1EC_0V; 
	      }
	  
	  if(sve[0] > 1 && sve[2] >= 1 && sve[1] <= 1)
	    { 
	      types[0] = HP_PRISM_1FB_2EB_0V;  
	      types[1] = HP_PRISM_1FA_1FB_2EB_0V;
	      types[2] = HP_PRISM_2FA_1FB_2EB_0V; 
	    }
	  
	  if(sve[0] > 1 && sve[1] <= 1 && sve[2] == 0) // ea && !eb  
	    {
	      types[0] = HP_PRISM_1FB_1EA_0V;
	      types[1] = HP_PRISM_1FA_1FB_1EA_0V;
	      types[2] = HP_PRISM_2FA_1FB_1EA_0V; 
	    } 
	  
	  if(sve[0] <= 1 && sve[1] > 1 && sve[2] == 0)
	    types[1] = HP_PRISM_1FA_1FB_1EB_0V; 
	  
	  if(sve[0] > 1 && sve[1]>1) 
	    if(sve[2] == 0)  // ea && eb 
	      {
		types[0] = HP_PRISM_1FB_2EA_0V;
		types[1] = HP_PRISM_1FA_1FB_2EA_0V;
		types[2] = HP_PRISM_2FA_1FB_2EA_0V; 
	      }
	  if(sve[0] <= 1 && sve[1] > 1 && sve[2] >0)
	    types[1] = HP_PRISM_1FA_1FB_2EC_0V; 
	  
	  if(sve[0] > 1 && sve[1] > 1 && sve[2] >= 1) //sve[2] can also be a face-edge  
	    {
	      types[0] = HP_PRISM_1FB_3E_0V;  
	      types[1] = HP_PRISM_1FA_1FB_3E_0V; 
	      types[2] = HP_PRISM_2FA_1FB_3E_0V; 
	    } 
	  
	  break;  
	  
	case 2:
	  if(sve[0] <= 1) 
	    cout << " **** WARNING: Edge between to different singular faces should be marked singular " << endl; 
		      
	  if(sve[1] <= 1)   
	    if(sve[2] <=1) 
	      { 
		types[0] = HP_PRISM_2FB_0E_0V; 
		types[1] = HP_PRISM_1FA_2FB_0E_0V;
		types[2] = HP_PRISM_2FA_2FB_0E_0V;
	      }
	    else
	      { 
		types[0] = HP_PRISM_2FB_1EC_0V; 
		types[1] = HP_PRISM_1FA_2FB_1EC_0V; 
		types[2] = HP_PRISM_2FA_2FB_1EC_0V;   
	      }
	  else
	    if(sve[2] <= 1) 
	      types[1] = HP_PRISM_1FA_2FB_1EB_0V; 
	    else
	      { 
		types[0] = HP_PRISM_2FB_3E_0V; 
		types[1] = HP_PRISM_1FA_2FB_3E_0V; 
		types[2] = HP_PRISM_2FA_2FB_3E_0V; 
	      }
	  
	  break;
	  
	case 3: 
	  types[0] = HP_PRISM_3FB_0V; 
	  types[1] = HP_PRISM_1FA_3FB_0V; 
	  types[2] = HP_PRISM_2FA_3FB_0V; 
	  break;
	}
      type = types[tfsing];
      
         
      if(type != HP_NONE)  
	break;
    }
	 
  /*
   *testout << " Prism with pnums " << endl; 
   for(int j=0;j<6;j++) *testout << el.pnums[j] << "\t"; 
   *testout << endl; 
   */
  
  if(type != HP_NONE) 
    {
      int pnums[6]; 
      for(int j=0;j<6;j++) pnums[j] = el.PNum (p[j]);
      for(int k=0;k<6;k++) el.pnums[k] = pnums[k]; 
    }

  /* *testout << " Classified Prism with pnums " << endl; 
  for(int j=0;j<6;j++) *testout << el.pnums[j] << "\t"; 
  *testout << endl; 
  */ 
  return(type); 
}


// #ifdef SABINE 
HPREF_ELEMENT_TYPE ClassifyTrig(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int>  & face_edges, 
				INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint, int dim, const FaceDescriptor & fd)

{
  HPREF_ELEMENT_TYPE type = HP_NONE; 
  
  int pnums[3]; 
  int p[3];   
  
  INDEX_3 i3 (el.pnums[0], el.pnums[1], el.pnums[2]);
  i3.Sort();
  bool sing_face = faces.Used (i3);
  
  // *testout << " facepoint " << facepoint << endl;  
      

  // Try all rotations of the trig 
  for (int j=0;j<3;j++) 
    {
      int point_sing[3] = {0,0,0}; 
      int edge_sing[3] = {0,0,0}; 
      // *testout << " actual rotation of trig points " ;  
      for(int m=0;m<3;m++) 
	{ 
	  p[m] = (j+m)%3 +1; // local vertex number
	  pnums[m] = el.PNum(p[m]); // global vertex number 
	  // *testout << pnums[m] << " \t "; 
	}
      // *testout << endl ; 
      
      if(dim == 3) 
	{
	  // face point 
	  for(int k=0;k<3;k++)
	    if(!sing_face)
	      { 
		//	*testout << " fp [" << k << "] = " << facepoint[pnums[k]] << endl;   
		//	*testout << " fd.DomainIn()" <<  fd.DomainIn() << endl; 
		//	*testout  << " fd.DomainOut()" <<  fd.DomainOut() << endl; 
		if( facepoint[pnums[k]]  && (facepoint[pnums[k]] ==-1 || 
					     facepoint[pnums[k]] == fd.DomainIn() ||   facepoint[pnums[k]] == fd.DomainOut()))
		  point_sing[p[k]-1] = 1; 
	      } 
	  // if point is on face_edge in next step sing = 2 

	  /*	  *testout << " pointsing NACH FACEPOints ... FALLS EDGEPOINT UMSETZEN" ; 
	  for (int k=0;k<3;k++) *testout << "\t" << point_sing[p[k]-1] ;
	  *testout << endl; */
	}
      
      const ELEMENT_EDGE * eledges = MeshTopology::GetEdges(TRIG); 
      
      if(dim==3)
	{
	  for(int k=0;k<3;k++) 
	    { 
	      int ep1=p[eledges[k][0]-1];  
	      int ep2=p[eledges[k][1]-1];  
	      INDEX_2 i2(el.PNum(ep1),el.PNum(ep2)); 
	      
	      if(edges.Used(i2)) 
		{
		  
		  edge_sing[k]=2;
		  point_sing[ep1-1] = 2; 
		  point_sing[ep2-1] = 2; 
		}
	      else // face_edge? 
		{	  
		  i2.Sort();  
		  if(surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr()+1)  // edge not face_edge acc. to surface in which trig lies
		    if(face_edges.Get(i2)==-1 ||face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut() )
		      { 
			edge_sing[k]=1;
		      } 
		    else
		      { 
			point_sing[ep1-1] = 0; // set to edge_point 
			point_sing[ep2-1] = 0; // set to edge_point
		      } 
		}
	      
	      /*  *testout << " pointsing NACH edges UND FACEEDGES UMSETZEN ... " ; 
	      for (int k=0;k<3;k++) *testout << "\t" << point_sing[p[k]-1] ; 
	      *testout << endl;          
	      */
	    }
	}
      /*
      *testout << " dim " << dim << endl; 
      *testout << " edgepoint_dom " << edgepoint_dom << endl; 
      */
      if(dim==2)
	{
	  for(int k=0;k<3;k++) 
	    { 
	      int ep1=p[eledges[k][0]-1];  
	      int ep2=p[eledges[k][1]-1];  
	     
	      INDEX_2 i2(el.PNum(ep1),el.PNum(ep2));  
	     
	      if(edges.Used(i2)) 
		{
		  
		  if(edgepoint_dom.Used(INDEX_2(fd.SurfNr(),pnums[ep1-1])) || 
		     edgepoint_dom.Used(INDEX_2(-1,pnums[ep1-1])) || 
		     edgepoint_dom.Used(INDEX_2(fd.SurfNr(),pnums[ep2-1])) || 
		     edgepoint_dom.Used(INDEX_2(-1,pnums[ep2-1]))) 
		    {
		      edge_sing[k]=2;
		      point_sing[ep1-1] = 2;
		      point_sing[ep2-1] = 2; 
		    }
		}
	     
	    }
	}

     
	 
      for(int k=0;k<3;k++) 
	if(edgepoint.Test(pnums[k])) //edgepoint, but not member of sing_edge on trig -> cp 
	  {
	    INDEX_2 i2a=INDEX_2::Sort(el.PNum(p[k]), el.PNum(p[(k+1)%3])); 
	    INDEX_2 i2b=INDEX_2::Sort(el.PNum(p[k]), el.PNum(p[(k+2)%3])); 
	    
	    if(!edges.Used(i2a) && !edges.Used(i2b)) 
	      point_sing[p[k]-1] = 3; 	
	  } 
      
      for(int k=0;k<3;k++) 
	if(cornerpoint.Test(el.PNum(p[k]))) 
	  point_sing[p[k]-1] = 3; 
      

     if(edge_sing[0] + edge_sing[1] + edge_sing[2] == 0) 
       { 
	 int ps = point_sing[0] + point_sing[1] + point_sing[2]; 
	 
	 if(ps==0) 
	   type = HP_TRIG; 
	 else if(point_sing[p[0]-1]  && !point_sing[p[1]-1] && !point_sing[p[2]-1])
	   type = HP_TRIG_SINGCORNER;
	 else if(point_sing[p[0]-1] && point_sing[p[1]-1] && !point_sing[p[2]-1]) 
	   type = HP_TRIG_SINGCORNER12; 
	 else if(point_sing[p[0]-1] && point_sing[p[1]-1] && point_sing[p[2]-1]) 
	   { 
	     if(dim==2) type = HP_TRIG_SINGCORNER123_2D; 
	     else type = HP_TRIG_SINGCORNER123; 
	   } 
       } 
     else
       if (edge_sing[2] && !edge_sing[0] && !edge_sing[1]) //E[2]=(1,2) 
      { 
	int code = 0; 
	if(point_sing[p[0]-1] > edge_sing[2]) code+=1; 
	if(point_sing[p[1]-1] > edge_sing[2]) code+=2; 
	if(point_sing[p[2]-1]) code+=4; 
	
	HPREF_ELEMENT_TYPE types[] =
	  {
	    HP_TRIG_SINGEDGE, 
	    HP_TRIG_SINGEDGECORNER1, 
	    HP_TRIG_SINGEDGECORNER2,
	    HP_TRIG_SINGEDGECORNER12, 
	    HP_TRIG_SINGEDGECORNER3, 
	    HP_TRIG_SINGEDGECORNER13, 
	    HP_TRIG_SINGEDGECORNER23, 
	    HP_TRIG_SINGEDGECORNER123, 
	  };
	type = types[code]; 
	
      }  // E[0] = [0,2], E[1] =[1,2], E[2] = [0,1]
    else 
      if(edge_sing[2] && !edge_sing[1] && edge_sing[0])
	{
	  if(point_sing[p[2]-1] <= edge_sing[0] ) 
	    { 
	      if(point_sing[p[1]-1]<= edge_sing[2]) type = HP_TRIG_SINGEDGES; 
	      else type = HP_TRIG_SINGEDGES2; 
	    } 
	  else 
	    {
	      if(point_sing[p[1]-1]<= edge_sing[2]) 
		type = HP_TRIG_SINGEDGES3; 
	      else type = HP_TRIG_SINGEDGES23; 
	    }
	}
      else if (edge_sing[2] && edge_sing[1] && edge_sing[0])
	type = HP_TRIG_3SINGEDGES; 
     
     //  cout << " run for " <<  j << " gives type " << type << endl; 
     //*testout << " run for " <<  j << " gives type " << type << endl; 
     if(type!=HP_NONE) break;
    }
    
  for(int k=0;k<3;k++) el[k] = pnums[k]; 
  /*if(type != HP_NONE) 
   {
     
      cout << " TRIG with pnums " << pnums[0] << "\t"  << 
	pnums[1] << "\t"  << pnums[2] << endl; 
      cout << " type "  << type << endl; 
      }
  */
  return(type);
}
#ifdef HPREF_OLD 
HPREF_ELEMENT_TYPE ClassifyTrig(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
				BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
				INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint, int dim, const FaceDescriptor & fd)
{
  HPREF_ELEMENT_TYPE type = HP_NONE; 
  
  int pnums[3]; 
	      
  INDEX_3 i3 (el.pnums[0], el.pnums[1], el.pnums[2]);
  i3.Sort();
  bool sing_face = faces.Used (i3);
   
  
  for (int j = 1; j <= 3; j++)
    {
      int ep1 = edgepoint.Test (el.PNumMod (j));
      int ep2 = edgepoint.Test (el.PNumMod (j+1));
      int ep3 = edgepoint.Test (el.PNumMod (j+2));
      
      if (dim == 2)
	{
	  // JS, Dec 11
	  ep1 = edgepoint_dom.Used (INDEX_2 (fd.SurfNr(), el.PNumMod(j))) ||
	    edgepoint_dom.Used (INDEX_2 (-1, el.PNumMod(j)));
	  ep2 = edgepoint_dom.Used (INDEX_2 (fd.SurfNr(), el.PNumMod(j+1))) ||
	    edgepoint_dom.Used (INDEX_2 (-1, el.PNumMod(j+1)));
	  ep3 = edgepoint_dom.Used (INDEX_2 (fd.SurfNr(), el.PNumMod(j+2))) ||
	    edgepoint_dom.Used (INDEX_2 (-1, el.PNumMod(j+2)));
	  /*
	  ep1 = edgepoint_dom.Used (INDEX_2 (el.index, el.PNumMod(j)));
	  ep2 = edgepoint_dom.Used (INDEX_2 (el.index, el.PNumMod(j+1)));
	  ep3 = edgepoint_dom.Used (INDEX_2 (el.index, el.PNumMod(j+2)));
	  */
	  // ep3 = edgepoint_dom.Used (INDEX_2 (mesh.SurfaceElement(i).GetIndex(), el.PNumMod(j+2)));
	}
      
      
      
      int cp1 = cornerpoint.Test (el.PNumMod (j));
      int cp2 = cornerpoint.Test (el.PNumMod (j+1));
      int cp3 = cornerpoint.Test (el.PNumMod (j+2));
      
      ep1 |= cp1;
      ep2 |= cp2;
      ep3 |= cp3;
      
      
      // (*testout) << "cp = " << cp1 << cp2 << cp3 << ", ep = " << ep1 << ep2 << ep3 << endl;

      int p[3] = { el.PNumMod (j), el.PNumMod (j+1), el.PNumMod (j+2)};
      if(ep1)
	{
	  INDEX_2 i2a=INDEX_2::Sort(p[0], p[1]); 
	  INDEX_2 i2b=INDEX_2::Sort(p[0], p[2]); 
	  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    cp1 = 1; 
	}
      if(ep2)
	{
	  INDEX_2 i2a=INDEX_2::Sort(p[1], p[0]); 
	  INDEX_2 i2b=INDEX_2::Sort(p[1], p[2]); 
	  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    cp2 = 1; 
	}
      if(ep3)
	{
	  INDEX_2 i2a=INDEX_2::Sort(p[2], p[0]); 
	  INDEX_2 i2b=INDEX_2::Sort(p[2], p[1]); 
	  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    cp3= 1; 
	}		      
      
      
      int isedge1=0, isedge2=0, isedge3=0; 
      if(dim == 3 )
	{
	  INDEX_2 i2;
	  i2 = INDEX_2(el.PNumMod (j), el.PNumMod (j+1));
	  isedge1 = edges.Used (i2);
	  i2.Sort();
	  if(surf_edges.Used(i2) &&  surf_edges.Get(i2)   != fd.SurfNr()+1 && 
	     (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
	    {
	      isedge1=1;
	      ep1 = 1; ep2=1;
	    }
	  
	  i2 = INDEX_2(el.PNumMod (j+1), el.PNumMod (j+2));
	  isedge2 = edges.Used (i2);
	  i2.Sort();
	  if(surf_edges.Used(i2) &&  surf_edges.Get(i2)   != fd.SurfNr()+1 &&
	     (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
	    {
	      isedge2=1;
	      ep2 = 1; ep3=1;
	    }
	  i2 = INDEX_2(el.PNumMod (j+2), el.PNumMod (j+3));
	  isedge3 = edges.Used (i2);
	  i2.Sort();
	  if(surf_edges.Used(i2) &&  surf_edges.Get(i2)   != fd.SurfNr()+1 && 
	     (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
	    {
	      isedge3=1;
	      ep1 = 1; ep3=1;
	    }
	  
	  // cout << " isedge " << isedge1 << " \t " << isedge2 << " \t " << isedge3 << endl;  
	
	  if (!sing_face)
	{
	  /*
	    if (!isedge1)  { cp1 |= ep1; cp2 |= ep2; }
	    if (!isedge2)  { cp2 |= ep2; cp3 |= ep3; }
	    if (!isedge3)  { cp3 |= ep3; cp1 |= ep1; }
	  */
	  ep1 |= facepoint [el.PNumMod(j)] != 0;
	  ep2 |= facepoint [el.PNumMod(j+1)] != 0;
	  ep3 |= facepoint [el.PNumMod(j+2)] != 0;
	  
	  
	  isedge1 |= face_edges.Used (INDEX_2::Sort (el.PNumMod(j), el.PNumMod(j+1)));
	  isedge2 |= face_edges.Used (INDEX_2::Sort (el.PNumMod(j+1), el.PNumMod(j+2)));
	  isedge3 |= face_edges.Used (INDEX_2::Sort (el.PNumMod(j+2), el.PNumMod(j+3)));
	}
	}
      
      if(dim ==2) 
	{ 
	  INDEX_2 i2;
	  i2 = INDEX_2(el.PNumMod (j), el.PNumMod (j+1));
	  i2.Sort();
	  isedge1 = edges.Used (i2);
	  if(isedge1)
	    {
	      ep1 = 1; ep2=1;
	    }
	  
	  i2 = INDEX_2(el.PNumMod (j+1), el.PNumMod (j+2));
	  i2.Sort();
	  isedge2 = edges.Used (i2);
	  if(isedge2)
	    {
	      ep2 = 1; ep3=1;
	    }
	  i2 = INDEX_2(el.PNumMod (j+2), el.PNumMod (j+3));
	  i2.Sort();
	  isedge3 = edges.Used (i2);
	  if(isedge3)
	    {
	      ep1 = 1; ep3=1;
	    }
	  
	  
	}
      
		  
		  /*
		  cout << " used " << face_edges.Used (INDEX_2::Sort (el.PNumMod(j), el.PNumMod(j+1))) << endl; 

		  cout << " isedge " << isedge1 << " \t " << isedge2 << " \t " << isedge3 << endl; 
		  cout << " ep " << ep1 << "\t" << ep2 << " \t " << ep3 << endl; 
		  cout << " cp " << cp1 << "\t" << cp2 << " \t " << cp3 << endl; 
		  */
		  

      
      if (isedge1 + isedge2 + isedge3 == 0)
	{
	  if (!ep1 && !ep2 && !ep3)
	    type = HP_TRIG;
	  
	  if (ep1 && !ep2 && !ep3)
	    type = HP_TRIG_SINGCORNER;
	  
	  if (ep1 && ep2 && !ep3)
	    type = HP_TRIG_SINGCORNER12;
	  
	  if (ep1 && ep2 && ep3)
	    {
	      if (dim == 2)
			    type = HP_TRIG_SINGCORNER123_2D;
	      else
		type = HP_TRIG_SINGCORNER123;
	    }
	  
	  if (type != HP_NONE)
	    {
	      pnums[0] = el.PNumMod (j);
	      pnums[1] = el.PNumMod (j+1);
	      pnums[2] = el.PNumMod (j+2);
	      break;
	    }
	}
      
      if (isedge1 && !isedge2 && !isedge3)
	{
	  int code = 0;
	  if (cp1) code += 1;
	  if (cp2) code += 2;
	  if (ep3) code += 4;
	  
	  HPREF_ELEMENT_TYPE types[] =
	    {
	      HP_TRIG_SINGEDGE, 
	      HP_TRIG_SINGEDGECORNER1, 
	      HP_TRIG_SINGEDGECORNER2,
	      HP_TRIG_SINGEDGECORNER12, 
	      HP_TRIG_SINGEDGECORNER3, 
	      HP_TRIG_SINGEDGECORNER13, 
	      HP_TRIG_SINGEDGECORNER23, 
	      HP_TRIG_SINGEDGECORNER123, 
	    };
	  type = types[code];
	  pnums[0] = el.PNumMod (j);
	  pnums[1] = el.PNumMod (j+1);
	  pnums[2] = el.PNumMod (j+2);
	  break;
	}
      
      
      if (isedge1 && !isedge2 && isedge3)
	{
	  if (!cp3)
	    {
	      if (!cp2) type = HP_TRIG_SINGEDGES;
	      else      type = HP_TRIG_SINGEDGES2;
	    }
	  else
	    { 
	      if (!cp2) type = HP_TRIG_SINGEDGES3;
	      else      type = HP_TRIG_SINGEDGES23;
	    }
	  
	  pnums[0] = el.PNumMod (j);
	  pnums[1] = el.PNumMod (j+1);
	  pnums[2] = el.PNumMod (j+2);
	  break;
	}
       
      if (isedge1 && isedge2 && isedge3)
	{
	  type = HP_TRIG_3SINGEDGES;
	  pnums[0] = el.PNumMod (j);
	  pnums[1] = el.PNumMod (j+1);
	  pnums[2] = el.PNumMod (j+2);
	  break;
	}
    }
  
  for(int k=0;k<3;k++) el[k] = pnums[k]; 
  /*if(type != HP_NONE) 
   {
     
      cout << " TRIG with pnums " << pnums[0] << "\t"  << 
	pnums[1] << "\t"  << pnums[2] << endl; 
      cout << " type "  << type << endl; 
      }
  */
  return(type);
}
#endif
HPREF_ELEMENT_TYPE ClassifyQuad(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int>  & face_edges, 
				INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint, int dim, const FaceDescriptor & fd)
{
  HPREF_ELEMENT_TYPE type = HP_NONE; 
  
  int ep1(-1), ep2(-1), ep3(-1), ep4(-1), cp1(-1), cp2(-1), cp3(-1), cp4(-1);
	      int isedge1, isedge2, isedge3, isedge4;

	      for (int j = 1; j <= 4; j++)
		{
		  ep1 = edgepoint.Test (el.PNumMod (j));
		  ep2 = edgepoint.Test (el.PNumMod (j+1));
		  ep3 = edgepoint.Test (el.PNumMod (j+2));
		  ep4 = edgepoint.Test (el.PNumMod (j+3));

		  if (dim == 2)
		    {
		      ep1 = edgepoint_dom.Used (INDEX_2 (el.GetIndex(), el.PNumMod(j)));
		      ep2 = edgepoint_dom.Used (INDEX_2 (el.GetIndex(), el.PNumMod(j+1)));
		      ep3 = edgepoint_dom.Used (INDEX_2 (el.GetIndex(), el.PNumMod(j+2)));
		      ep4 = edgepoint_dom.Used (INDEX_2 (el.GetIndex(), el.PNumMod(j+3)));
		    }

		  cp1 = cornerpoint.Test (el.PNumMod (j));
		  cp2 = cornerpoint.Test (el.PNumMod (j+1));
		  cp3 = cornerpoint.Test (el.PNumMod (j+2));
		  cp4 = cornerpoint.Test (el.PNumMod (j+3));

		  ep1 |= cp1;
		  ep2 |= cp2;
		  ep3 |= cp3;
		  ep4 |= cp4;
		
		  int p[4] = { el.PNumMod (j), el.PNumMod (j+1), el.PNumMod (j+2), el.PNumMod(j+4)};
		  //int epp[4] = { ep1, ep2, ep3, ep4}; 
		  int cpp[4] = { cp1, cp2, cp3, cp4};
		  for(int k=0;k<0;k++)
		    {
		      INDEX_2 i2a=INDEX_2::Sort(p[k], p[(k+1)%4]); 
		      INDEX_2 i2b=INDEX_2::Sort(p[k], p[(k-1)%4]); 
		      if(!edges.Used(i2a) && !edges.Used(i2b)) 
			cpp[k] = 1; 
		    }
		  cp1= cpp[0]; cp2=cpp[1]; cp3=cpp[2]; cp4=cpp[3];
		  

		  if(dim ==3) 
		    { 
		      INDEX_2 i2;
		      i2 = INDEX_2(el.PNumMod (j), el.PNumMod (j+1));
		      // i2.Sort();
		      isedge1 = edges.Used (i2);
		      i2.Sort();
		      if(surf_edges.Used(i2)  &&  surf_edges.Get(i2)   != fd.SurfNr()+1 && 
			 (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
			{
			  isedge1=1;
			  ep1 = 1; ep2=1;
			}
		      i2 = INDEX_2(el.PNumMod (j+1), el.PNumMod (j+2));
		      // i2.Sort();
		      isedge2 = edges.Used (i2);
		      i2.Sort();
		      if(surf_edges.Used(i2)  &&  surf_edges.Get(i2)   != fd.SurfNr()+1 && 
			 (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
			{  
			  isedge2=1;
			  ep2=1; ep3=1;
			}
		      i2 = INDEX_2(el.PNumMod (j+2), el.PNumMod (j+3));
		      // i2.Sort();
		      isedge3 = edges.Used (i2); 
		      i2.Sort();
		      if(surf_edges.Used(i2)   &&  surf_edges.Get(i2)   != fd.SurfNr()+1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
			{
			  isedge3=1;
			  ep3=1; ep4=1;
			}
		      i2 = INDEX_2(el.PNumMod (j+3), el.PNumMod (j+4));
		      // i2.Sort();
		      isedge4 = edges.Used (i2);
		      i2.Sort();
		      if(surf_edges.Used(i2)  &&  surf_edges.Get(i2)   != fd.SurfNr()+1 && 
			 (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()) ) 
			{ 
			  isedge4=1;
			  ep4=1; ep1=1;
			} 
		    

//MH***********************************************************************************************************
		      if(ep1)
			if(edgepoint.Test(p[0]))
			{
			  INDEX_2 i2a=INDEX_2::Sort(p[0], p[1]); 
	  		  INDEX_2 i2b=INDEX_2::Sort(p[0], p[3]); 
	  		  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    		    cp1 = 1; 
			}
		      if(ep2)
			if(edgepoint.Test(p[1]))
			{
			  INDEX_2 i2a=INDEX_2::Sort(p[0], p[1]); 
	  		  INDEX_2 i2b=INDEX_2::Sort(p[1], p[2]); 
	  		  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    		    cp2 = 1; 
			}
		      if(ep3)
			if(edgepoint.Test(p[2]))
			{
			  INDEX_2 i2a=INDEX_2::Sort(p[2], p[1]); 
	  		  INDEX_2 i2b=INDEX_2::Sort(p[3], p[2]); 
	  		  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    		    cp3 = 1; 
			}
		      if(ep4)
			if(edgepoint.Test(p[3]))
			{
			  INDEX_2 i2a=INDEX_2::Sort(p[0], p[3]); 
	  		  INDEX_2 i2b=INDEX_2::Sort(p[3], p[2]); 
	  		  if(!edges.Used(i2a) && !edges.Used(i2b)) 
	    		    cp4 = 1; 
			}
//MH*****************************************************************************************************************************
		}
		  else
		    { 
		      INDEX_2 i2;
		      i2 = INDEX_2(el.PNumMod (j), el.PNumMod (j+1));
		      i2.Sort();
		      isedge1 = edges.Used (i2);
		      if(isedge1)
			{
			  ep1 = 1; ep2=1;
			}
		      i2 = INDEX_2(el.PNumMod (j+1), el.PNumMod (j+2));
		      i2.Sort();
		      isedge2 = edges.Used (i2);
		      if(isedge2)
			{  
			  ep2=1; ep3=1;
			}
		      i2 = INDEX_2(el.PNumMod (j+2), el.PNumMod (j+3));
		      i2.Sort();
		      isedge3 = edges.Used (i2); 
		      
		      if(isedge3)
			{
			  ep3=1; ep4=1;
			}
		      i2 = INDEX_2(el.PNumMod (j+3), el.PNumMod (j+4));
		      i2.Sort();
		      isedge4 = edges.Used (i2);
		      if(isedge4)
			{ 
			  ep4=1; ep1=1;
			} 
		    }

		  int sumcp = cp1 + cp2 + cp3 + cp4;
		  int sumep = ep1 + ep2 + ep3 + ep4;
		  int sumedge = isedge1 + isedge2 + isedge3 + isedge4;

		  switch (sumedge)
		    {
		    case 0:
		      {
			switch (sumep)
			  {
			  case 0: 
			    type = HP_QUAD; 
			    break;
			  case 1: 
			    if (ep1) type = HP_QUAD_SINGCORNER;
			    break; 
			  case 2:
			    {
			      if (ep1 && ep2) type = HP_QUAD_0E_2VA;
			      if (ep1 && ep3) type = HP_QUAD_0E_2VB;
			      break;
			    }
			  case 3: 
			    if (!ep4) type = HP_QUAD_0E_3V; 
			    break; 
			  case 4: 
			    type = HP_QUAD_0E_4V; 
			    break; 
			  }
			break;
		      }
		    case 1:
		      {
			if (isedge1)
			  {
			    switch (cp1+cp2+ep3+ep4)
			      {
			      case 0: 
				type = HP_QUAD_SINGEDGE; 
				break;
			      case 1:
				{
				  if (cp1) type = HP_QUAD_1E_1VA;
				  if (cp2) type = HP_QUAD_1E_1VB;
				  if (ep3) type = HP_QUAD_1E_1VC;
				  if (ep4) type = HP_QUAD_1E_1VD; 
				  break; 
				}
			      case 2:
				{
				  if (cp1 && cp2) type = HP_QUAD_1E_2VA; 
				  if (cp1 && ep3) type = HP_QUAD_1E_2VB; 
				  if (cp1 && ep4) type = HP_QUAD_1E_2VC; 
				  if (cp2 && ep3) type = HP_QUAD_1E_2VD; 
				  if (cp2 && ep4) type = HP_QUAD_1E_2VE; 
				  if (ep3 && ep4) type = HP_QUAD_1E_2VF; 
				  break; 
				}
			      case 3:
				{
				  if (cp1 && cp2 && ep3) type = HP_QUAD_1E_3VA;
				  if (cp1 && cp2 && ep4) type = HP_QUAD_1E_3VB;
				  if (cp1 && ep3 && ep4) type = HP_QUAD_1E_3VC;
				  if (cp2 && ep3 && ep4) type = HP_QUAD_1E_3VD;
				  break;
				}
			      case 4:
				{
				  type = HP_QUAD_1E_4V; 
				  break;
				}
			      }
			  }
			break;
		      }
		    case 2:
		      {
			if (isedge1 && isedge4)
			  {
			    if (!cp2 && !ep3 && !cp4)
			      type = HP_QUAD_2E;
			  
			    if (cp2 && !ep3 && !cp4)
			      type = HP_QUAD_2E_1VA;
			    if (!cp2 && ep3 && !cp4)
			      type = HP_QUAD_2E_1VB;
			    if (!cp2 && !ep3 && cp4)
			      type = HP_QUAD_2E_1VC;

			    if (cp2 && ep3 && !cp4)
			      type = HP_QUAD_2E_2VA;
			    if (cp2 && !ep3 && cp4)
			      type = HP_QUAD_2E_2VB;
			    if (!cp2 && ep3 && cp4)
			      type = HP_QUAD_2E_2VC;

			    if (cp2 && ep3 && cp4)
			      type = HP_QUAD_2E_3V;
			  }
			if (isedge1 && isedge3)
			  {
			    switch (sumcp)
			      {
			      case 0: 
				type = HP_QUAD_2EB_0V; break;
			      case 1:
				{
				  if (cp1) type = HP_QUAD_2EB_1VA; 
				  if (cp2) type = HP_QUAD_2EB_1VB; 
				  break;
				}
			      case 2:
				{
				  if (cp1 && cp2) { type = HP_QUAD_2EB_2VA; }
				  if (cp1 && cp3) { type = HP_QUAD_2EB_2VB; }
				  if (cp1 && cp4) { type = HP_QUAD_2EB_2VC; }
				  if (cp2 && cp4) { type = HP_QUAD_2EB_2VD; }
				  break;
				}
			      case 3:
				{
				  if (cp1 && cp2 && cp3) { type = HP_QUAD_2EB_3VA; }
				  if (cp1 && cp2 && cp4) { type = HP_QUAD_2EB_3VB; }
				  break;
				}
			      case 4:
				{
				  type = HP_QUAD_2EB_4V; break;
				}
			      }
			  }
			break;
		      }

		    case 3:
		      {
			if (isedge1 && isedge2 && isedge4)
			  {
			    if (!cp3 && !cp4) type = HP_QUAD_3E;
			    if (cp3 && !cp4) type = HP_QUAD_3E_3VA;
			    if (!cp3 && cp4) type = HP_QUAD_3E_3VB;
			    if (cp3 && cp4) type = HP_QUAD_3E_4V;
			  }
			break;
		      }

		    case 4:
		      {
			type = HP_QUAD_4E;
			break;
		      }
		    }

		  if (type != HP_NONE)
		    {
		      int pnums[4]; 
		      pnums[0] = el.PNumMod (j); 
		      pnums[1] = el.PNumMod (j+1);
		      pnums[2] = el.PNumMod (j+2); 
		      pnums[3] = el.PNumMod (j+3);
		      for (int k=0;k<4;k++) el[k] = pnums[k]; 
	
		      /*  cout << " QUAD with pnums " << pnums[0] << "\t"  << 
			pnums[1] << "\t"  << pnums[2] << "\t"  << pnums[3] 
			<< endl << " of type " << type << endl; */
		   		      
		      break;
		    }
		}
	      if (type == HP_NONE)
		{
		  (*testout) << "undefined element" << endl
			     << "cp = " << cp1 << cp2 << cp3 << cp4 << endl
			     << "ep = " << ep1 << ep2 << ep3 << ep4 << endl
			     << "isedge = " << isedge1 << isedge2 << isedge3 
			     << isedge4 << endl;
		}
	    

	      return(type);  
}	    


HPREF_ELEMENT_TYPE ClassifyHex(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
		 INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint)
{
  HPREF_ELEMENT_TYPE type = HP_NONE;
  
  // implementation only for HP_HEX_1F_0E_0V
  //                         HP_HEX_1FA_1FB_0E_0V
  //                         HP_HEX 
  // up to now other cases are refine dummies 
  
  // indices of bot,top-faces combinations
  int index[6][2] = {{0,1},{1,0},{2,4},{4,2},{3,5},{5,3}}; 
  int p[8]; 
  const ELEMENT_FACE * elfaces  = MeshTopology::GetFaces (HEX);
  const ELEMENT_EDGE * eledges = MeshTopology::GetEdges (HEX);
  
  for(int m=0;m<6 && type == HP_NONE;m++) 
    for(int j=0;j<4 && type == HP_NONE;j++) 
      { 
	int point_sing[8]={0,0,0,0,0,0,0,0}; 
	int face_sing[6] = {0,0,0,0,0,0};
	int edge_sing[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	int spoint=0, sface=0, sedge=0; 
	for(int l=0;l<4;l++) 
	  {
	    p[l] = elfaces[index[m][0]][(4-j-l)%4]; 
	    p[l+4] = elfaces[index[m][1]][(j+l)%4];
	  }
	
	for(int l=0;l<8;l++) 
	  if(cornerpoint.Test(el.PNum(p[l])))  
	    { 
	      point_sing[p[l]-1]=3;
	      spoint++; 
	    }
	  else if(edgepoint.Test(el.PNum(p[l]))) point_sing[p[l]-1]=2;
	  else if (facepoint[el.PNum(p[l])] == -1 || facepoint[el.PNum(p[l])] == el.GetIndex())
	    point_sing[p[l]-1] = 1;   
	
	for(int k=0;k<12;k++)
	{
	  INDEX_2 i2 = INDEX_2 :: Sort(el.PNum(p[eledges[k][0]-1]),el.PNum(p[eledges[k][1]-1])); 
	  if (edges.Used(i2)) 
	    { 
	      edge_sing[k] = 2;
	      sedge++; 
	    }
	  else edge_sing[k] = face_edges.Used(i2);
	}
	
	for (int k=0;k<6;k++)
	{
	  INDEX_3 i3; 
	  
	
	  INDEX_4  i4 = INDEX_4(el.pnums[p[elfaces[k][0]-1]-1], el.pnums[p[elfaces[k][1]-1]-1], el.pnums[p[elfaces[k][2]-1]-1],el.pnums[p[elfaces[k][3]-1]-1]); 
	  i4.Sort();
	  i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3()); 
	  
	  if (faces.Used (i3))
	    {
	      
	      int domnr = faces.Get(i3); 
	      if (domnr == -1 || domnr == el.GetIndex())
		{
		  face_sing[k] = 1;
		  sface++;
		}
	      
	    } 
	} 
	
	if(!sface && !sedge && !spoint) type = HP_HEX; 
	if(!sedge && !spoint) 
	  {
	    if(face_sing[0] && face_sing[2] && sface==2) 
	      type = HP_HEX_1FA_1FB_0E_0V; 
	    if (face_sing[0] && sface==1)  
	      type = HP_HEX_1F_0E_0V; 
	  }
	
	el.type=type; 

	if(type != HP_NONE) 
	  {
	    int pnums[8]; 
	    for(int l=0;l<8;l++) pnums[l] = el[p[l]-1];
	    for(int l=0;l<8;l++) el[l] = pnums[l];
	    /* cout << " HEX with pnums " << pnums[0] << "\t"  << 
	      pnums[1] << "\t"  << pnums[2] << "\t"  << pnums[3] << "\t"  << 
	      pnums[4] << "\t"  <<  pnums[5] << endl << " of type " << type << endl; */
	    break; 
	  }
      }
  
  return (type); 

}

HPREF_ELEMENT_TYPE ClassifySegm(HPRefElement & hpel, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
		 INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint)
{
  
  int cp1 = cornerpoint.Test (hpel[0]);
  int cp2 = cornerpoint.Test (hpel[1]);
  
  INDEX_2 i2;
  i2 = INDEX_2(hpel[0], hpel[1]);
  i2.Sort();
  if (!edges.Used (i2))
    {
      cp1 = edgepoint.Test (hpel[0]);
      cp2 = edgepoint.Test (hpel[1]);
    }
  
  if(!edges.Used(i2) && !face_edges.Used(i2))
    {
      if(facepoint[hpel[0]]!=0) cp1=1; 
      if(facepoint[hpel[1]]!=0) cp2=1; 
    }
  
  if(edges.Used(i2) && !face_edges.Used(i2))
    {
      if(facepoint[hpel[0]]) cp1 = 1; 
	    if(facepoint[hpel[1]]) cp2 = 1; 
    }
  
  if (!cp1 && !cp2)
    hpel.type = HP_SEGM;
  else if (cp1 && !cp2)
    hpel.type = HP_SEGM_SINGCORNERL;
  else if (!cp1 && cp2)
    hpel.type = HP_SEGM_SINGCORNERR;
  else
    hpel.type = HP_SEGM_SINGCORNERS;
  
  // cout << " SEGM found with " << hpel[0] << " \t" << hpel[1] << endl << " of type " << hpel.type << endl; 
  return(hpel.type) ; 
} 


HPREF_ELEMENT_TYPE ClassifyPyramid(HPRefElement & el, INDEX_2_HASHTABLE<int> & edges, INDEX_2_HASHTABLE<int> & edgepoint_dom, 
		 BitArray & cornerpoint, BitArray & edgepoint, INDEX_3_HASHTABLE<int> & faces, INDEX_2_HASHTABLE<int> & face_edges, 
		 INDEX_2_HASHTABLE<int> & surf_edges, ARRAY<int, PointIndex::BASE> & facepoint)
{
  HPREF_ELEMENT_TYPE type = HP_NONE;
  
  // implementation only for HP_PYRAMID
  //                         HP_PYRAMID_0E_1V
  //                         HP_PYRAMID_EDGES
  //                         HP_PYRAMID_1FB_0E_1VA
  // up to now other cases are refine dummies 
  
  // indices of bot,top-faces combinations
  // int index[6][2] = {{0,1},{1,0},{2,4},{4,2},{3,5},{5,3}}; 

  const ELEMENT_FACE * elfaces  = MeshTopology::GetFaces (PYRAMID);
  const ELEMENT_EDGE * eledges = MeshTopology::GetEdges (PYRAMID);
  
  int point_sing[5]={0,0,0,0,0}; 
  int face_sing[5] = {0,0,0,0,0};
  int edge_sing[8] = {0,0,0,0,0,0,0,0};
  
  int spoint=0, sedge=0, sface=0; 
   
  for(int m=0;m<4 && type == HP_NONE;m++) 
    {
      int p[5] = {m%4, m%4+1, m%4+2, m%4+3, 4}; 

      for(int l=0;l<5;l++) 
	{
	  if(cornerpoint.Test(el.pnums[p[l]]))  
	    point_sing[l]=3;
	  
	     else if(edgepoint.Test(el.pnums[p[l]]))
	    point_sing[l]=2;
	  
	  else if (facepoint[el.pnums[p[l]]] == -1 || facepoint[el.pnums[p[l]]] == el.GetIndex())
	    point_sing[l] = 1;   
	  
	  spoint += point_sing[l]; 
	}
      
      for(int k=0;k<8;k++)
	{
	  INDEX_2 i2 = INDEX_2 :: Sort(el.pnums[p[eledges[k][0]-1]],
				       el.pnums[p[eledges[k][1]-1]]); 
	  if (edges.Used(i2)) 
	    edge_sing[k] = 2;
	  else 
	    edge_sing[k] = face_edges.Used(i2);
	  
	  sedge += edge_sing[k]; 
	}
  
      for (int k=0;k<5;k++)
	{
	  INDEX_3 i3; 
	  INDEX_4  i4 = INDEX_4(el.pnums[p[elfaces[k][0]-1]], el.pnums[p[elfaces[k][1]-1]], el.pnums[p[elfaces[k][2]-1]],
				el.pnums[p[elfaces[k][3]-1]]); 
	  i4.Sort();
	  i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3()); 
	  
	  if (faces.Used (i3))
	    {
	      
	      int domnr = faces.Get(i3); 
	      if (domnr == -1 || domnr == el.GetIndex())
		face_sing[k] = 1;
	    } 
	  sface +=face_sing[k]; 
	} 
  
      if(!sface && !spoint && !sedge) return(HP_PYRAMID); 
      
      if(!sface && !sedge && point_sing[p[0]] == spoint) 
	type = HP_PYRAMID_0E_1V; 
      
      if(!sface && edge_sing[0] + edge_sing[2] == sedge && 
	 spoint == point_sing[0] + point_sing[1] + point_sing[3]) 
	type = HP_PYRAMID_EDGES; 
      
      if(sface && sface == face_sing[0] && spoint == point_sing[4] + 2)
	type = HP_PYRAMID_1FB_0E_1VA; 
      
      
      if(type != HP_NONE) 
	{ 
	  int pnums[8]; 
	  for(int l=0;l<5;l++) pnums[l] = el[p[l]];
	  for(int l=0;l<5;l++) el[l] = pnums[l];
	  el.type=type; 
	  break; 
	} 
    }
  
  return (type); 
  
}
