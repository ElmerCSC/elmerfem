
#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <acisgeom.hpp>
#include <meshing.hpp>

namespace netgen
{

#include "writeuser.hpp"
  
  
  void WriteTETFormat (const Mesh & mesh,
		       const string & filename)//, const string& problemType )
  {
    string problemType = "";
    if(!mesh.PureTetMesh())
      throw NgException("Can only export pure tet mesh in this format");

    cout << "starting .tet export to file " << filename << endl;


    ARRAY<int> point_ids,edge_ids,face_ids;
    ARRAY<int> elnum(mesh.GetNE());
    elnum = -1;

    
    ARRAY<int> userdata_int;
    ARRAY<double> userdata_double;
    ARRAY<int> ports;

    ARRAY<int> uid_to_group_3D, uid_to_group_2D, uid_to_group_1D, uid_to_group_0D;

    int pos_int = 0;
    int pos_double = 0;
    
    bool haveuserdata = 
      (mesh.GetUserData("TETmesh:double",userdata_double) &&
       mesh.GetUserData("TETmesh:int",userdata_int) && 
       mesh.GetUserData("TETmesh:ports",ports) &&
       mesh.GetUserData("TETmesh:point_id",point_ids,PointIndex::BASE) &&
       mesh.GetUserData("TETmesh:uid_to_group_3D",uid_to_group_3D) &&
       mesh.GetUserData("TETmesh:uid_to_group_2D",uid_to_group_2D) &&
       mesh.GetUserData("TETmesh:uid_to_group_1D",uid_to_group_1D) &&
       mesh.GetUserData("TETmesh:uid_to_group_0D",uid_to_group_0D));


    int version,subversion;

    if(haveuserdata)
      {
	version = int(userdata_double[0]);
	subversion = int(10*(userdata_double[0] - version));
	pos_double++;
      }
    else
      {
	version = 2;
	subversion = 0;
      }

    
    if(version >= 2)
      {
	// test if ids are disjunct, if not version 2.0 not possible
	int maxbc(-1),mindomain(-1);
	
	for(ElementIndex i=0; i<mesh.GetNE(); i++)
	  if(i==0 || mesh[i].GetIndex() < mindomain)
	    mindomain = mesh[i].GetIndex();
	for(int i=1; i<=mesh.GetNFD(); i++)
	  if(i==1 || mesh.GetFaceDescriptor(i).BCProperty() > maxbc)
	    maxbc = mesh.GetFaceDescriptor(i).BCProperty();
	
	if(maxbc >= mindomain)
	  {
	    cout << "WARNING: writing version " << version << "." << subversion << " tetfile not possible, ";
	    version = 1; subversion = 1;
	    cout << "using version " << version << "." << subversion << endl;
	  }
      }



    int startsize = point_ids.Size();
    point_ids.SetSize(mesh.GetNP()+1);
    for(int i=startsize; i<point_ids.Size(); i++)
      point_ids[i] = -1;


    for(int i=0; i<PointIndex::BASE; i++)
      point_ids[i] = -1;


    INDEX_2_CLOSED_HASHTABLE<int> edgenumbers(6*mesh.GetNE()+3*mesh.GetNSE());;
    INDEX_3_CLOSED_HASHTABLE<int> facenumbers(4*mesh.GetNE()+mesh.GetNSE());

    ARRAY<INDEX_2> edge2node;
    ARRAY<INDEX_3> face2edge;
    ARRAY<INDEX_4> element2face;

    int numelems(0),numfaces(0),numedges(0),numnodes(0);

    for(SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
      {
	const Segment & seg = mesh[si];
	INDEX_2 i2(seg.p1,seg.p2);
	i2.Sort();
	if(edgenumbers.Used(i2))
	  continue;

	numedges++;
	edgenumbers.Set(i2,numedges);
	edge2node.Append(i2);

	edge_ids.Append(seg.edgenr);

	if(point_ids[seg.p1] == -1)
	  point_ids[seg.p1] = (version >= 2) ? seg.edgenr : 0;
	if(point_ids[seg.p2] == -1)
	  point_ids[seg.p2] = (version >= 2) ? seg.edgenr : 0;
      }

    for(SurfaceElementIndex si = 0; si < mesh.GetNSE(); si++)
      {
	if(mesh[si].IsDeleted())
	  continue;

	const Element2d & elem = mesh[si];

	numfaces++;
	INDEX_3 i3(elem[0], elem[1], elem[2]);

	int min = i3[0];
	int minpos = 0;
	for(int j=1; j<3; j++)
	  if(i3[j] < min)
	    {
	      min = i3[j]; minpos = j;
	    }
	if(minpos == 1)
	  {
	    int aux = i3[0]; i3[0] = i3[1]; i3[1] = i3[2]; i3[2] = aux;
	  }
	else if(minpos == 2)
	  {
	    int aux = i3[0]; i3[0] = i3[2]; i3[2] = i3[1]; i3[1] = aux;
	  }
	facenumbers.Set(i3,numfaces);

	int bc = mesh.GetFaceDescriptor(elem.GetIndex()).BCProperty();
	face_ids.Append(bc);

	for(int j=0; j<3; j++)
	  if(point_ids[elem[j]] == -1)
	    point_ids[elem[j]] = (version >= 2) ? bc : 0;

	INDEX_2 i2a,i2b;
	INDEX_3 f_to_n;
	for(int j=0; j<3; j++)
	  {
	    i2a = INDEX_2(i3[j],i3[(j+1)%3]);
	    i2b[0] = i2a[1]; i2b[1] = i2a[0];
	    if(edgenumbers.Used(i2a))
	      f_to_n[j] = edgenumbers.Get(i2a);
	    else if(edgenumbers.Used(i2b))
	      f_to_n[j] = -edgenumbers.Get(i2b);
	    else
	      {
		numedges++;
		edgenumbers.Set(i2a,numedges);
		edge2node.Append(i2a);
		f_to_n[j] = numedges;
		if(version >= 2)
		  edge_ids.Append(bc);
		else
		  edge_ids.Append(0);
	      }
	  }
	face2edge.Append(f_to_n);
      }
    
    for(ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      {
	const Element & el = mesh[ei];

	if(el.IsDeleted())
	  continue;

	numelems++;
	elnum[ei] = numelems;

	static int tetfaces[4][3] =
	  { { 0, 2, 1 },
	    { 0, 1, 3 },
	    { 1, 2, 3 },
	    { 2, 0, 3 } };
	
	for(int j=0; j<4; j++)
	  if(point_ids[el[j]] == -1)
	    point_ids[el[j]] = (version >= 2) ? el.GetIndex() : 0;

	INDEX_4 e_to_f;

	for(int i = 0; i < 4; i++)
	  {
	    INDEX_3 i3a(el[tetfaces[i][0]],el[tetfaces[i][1]],el[tetfaces[i][2]]);
	    
	    int min = i3a[0];
	    int minpos = 0;
	    for(int j=1; j<3; j++)
	      if(i3a[j] < min)
		{
		  min = i3a[j]; minpos = j;
		}
	    if(minpos == 1)
	      {
		int aux = i3a[0]; i3a[0] = i3a[1]; i3a[1] = i3a[2]; i3a[2] = aux;
	      }
	    else if(minpos == 2)
	      {
		int aux = i3a[0]; i3a[0] = i3a[2]; i3a[2] = i3a[1]; i3a[1] = aux;
	      }
	    INDEX_3 i3b(i3a[0],i3a[2],i3a[1]);
	    

	    if(facenumbers.Used(i3a))
	      e_to_f[i] = facenumbers.Get(i3a);
	    else if(facenumbers.Used(i3b))
	      e_to_f[i] = -facenumbers.Get(i3b);
	    else
	      {
		numfaces++;
		facenumbers.Set(i3a,numfaces);
		e_to_f[i] = numfaces;
		if(version >= 2)
		  face_ids.Append(el.GetIndex());
		else
		  face_ids.Append(0);

		INDEX_2 i2a,i2b;
		INDEX_3 f_to_n;
		for(int j=0; j<3; j++)
		  {
		    i2a = INDEX_2(i3a[j],i3a[(j+1)%3]);
		    i2b[0] = i2a[1]; i2b[1] = i2a[0];
		    if(edgenumbers.Used(i2a))
		      f_to_n[j] = edgenumbers.Get(i2a);
		    else if(edgenumbers.Used(i2b))
		      f_to_n[j] = -edgenumbers.Get(i2b);
		    else
		      {
			numedges++;
			edgenumbers.Set(i2a,numedges);
			edge2node.Append(i2a);
			f_to_n[j] = numedges;
			if(version >= 2)
			  edge_ids.Append(el.GetIndex());
			else
			  edge_ids.Append(0);
		      }
		  }
		face2edge.Append(f_to_n);	  
	      }
	  }
	element2face.Append(e_to_f);
      }




    ofstream outfile(filename.c_str());

    outfile.precision(16);

    int unitcode;
    double tolerance;
    double dS1,dS2, alphaDeg;
    double x3D,y3D,z3D;
    int modelverts(0), modeledges(0), modelfaces(0), modelcells(0);

    int numObj0D,numObj1D,numObj2D,numObj3D;
    int numports = ports.Size();

    ARRAY<int> nodenum(point_ids.Size()+1);

    nodenum = -1;
	    


    numnodes = 0;
    for(int i=0; i<point_ids.Size(); i++)
      {
	if(point_ids[i] != -1)
	  {
	    numnodes++;
	    nodenum[i] = numnodes;
	  }
      }


    if(haveuserdata)
      {
	unitcode = userdata_int[pos_int];
	pos_int++;

	tolerance = userdata_double[pos_double];
	pos_double++;

	dS1 = userdata_double[pos_double];
	pos_double++;
	dS2 = userdata_double[pos_double];
	pos_double++;
	alphaDeg = userdata_double[pos_double];
	pos_double++;

	x3D = userdata_double[pos_double];
	pos_double++;
	y3D = userdata_double[pos_double];
	pos_double++;
	z3D = userdata_double[pos_double];
	pos_double++;

	if(version == 2)
	  {
	    modelverts = userdata_int[pos_int];
	    pos_int++;
	    modeledges = userdata_int[pos_int];
	    pos_int++;
	    modelfaces = userdata_int[pos_int];
	    pos_int++;
	    modelcells = userdata_int[pos_int];
	    pos_int++;
	  }

	numObj3D = userdata_int[pos_int];
	pos_int++;
	numObj2D = userdata_int[pos_int];
	pos_int++;
	numObj1D = userdata_int[pos_int];
	pos_int++;
	numObj0D = userdata_int[pos_int];
	pos_int++;
      }
    else
      {
	unitcode = 3;

	tolerance = 1e-5;

	dS1 = dS2 = alphaDeg = 0;

	x3D = y3D = z3D = 0;

	modelverts = modeledges = modelfaces = modelcells = 0;
	
	numObj3D = numObj2D = numObj1D = numObj0D = 0;
      }

    string uidpid;
    if(version == 1)
      uidpid = "PID";
    else if (version == 2)
      uidpid = "UID";
    

    ARRAY< ARRAY<int,PointIndex::BASE>* > idmaps;
    for(int i=1; i<=mesh.GetIdentifications().GetMaxNr(); i++)
      {
	if(mesh.GetIdentifications().GetType(i) == Identifications::PERIODIC)
	  {
	    idmaps.Append(new ARRAY<int,PointIndex::BASE>);
	    mesh.GetIdentifications().GetMap(i,*idmaps.Last(),true);
	  }
      }

    ARRAY<int> id_num,id_type;
    ARRAY< ARRAY<int> *> id_groups;


	// sst 2008-03-12: Write problem class...
	{
		std::string block;
		block  = "// CST Tetrahedral ";
		block += !problemType.empty() ? problemType : "High Frequency";
		block += " Mesh, Version no.:\n";
		
		size_t size = block.size()-3;
		block += "// ";
		block.append( size, '^' );
		block += "\n";

		outfile
			<< block
			<< version << "." << subversion << "\n\n";
	}

	outfile 
	    << "// User Units Code (1=CM 2=MM 3=M 4=MIC 5=NM 6=FT 7=IN 8=MIL):\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << unitcode << "\n\n"					\
	    << "// Geometric coord \"zero\" tolerance threshold:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << tolerance << "\n\n"				  \
	    << "// Periodic UnitCell dS1 , dS2 , alphaDeg:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << dS1 << " " << dS2 << " " << alphaDeg <<"\n\n"	\
	    << "// Periodic UnitCell origin in global coords (x3D,y3D,z3D):\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << x3D << " " << y3D << " " << z3D << "\n" << endl;

    if(version == 2)
      {
	outfile << "// Model entity count: Vertices, Edges, Faces, Cells:\n" \
		<< "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
		<< modelverts << " " << modeledges << " " << modelfaces << " " << modelcells << endl << endl;
      }


    outfile << "// Topological mesh-entity counts (#elements,#faces,#edges,#nodes):\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    outfile << numelems << " " 
	    << numfaces << " "
	    << numedges << " " 
	    << numnodes << endl << endl;

    outfile << "// NodeID, X, Y, Z, Type (0=Reg 1=PMaster 2=PSlave 3=CPMaster 4=CPSlave), "<< uidpid <<":\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
       


    id_num.SetSize(mesh.GetNP()+1);
    id_type.SetSize(mesh.GetNP()+1);
    id_num = 0;
    id_type = 0;

    int n2,n4,n8;
    n2 = n4 = n8 = 0;

 
    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
      {
	if(id_num[i] != 0)
	  continue;

	if(nodenum[i] == -1)
	  continue;

	ARRAY<int> group;
	group.Append(i);
	for(int j=0; j<idmaps.Size(); j++)
	  {
	    startsize = group.Size();
	    for(int k=0; k<startsize; k++)
	      {
		int id = (*idmaps[j])[group[k]];
		if(id != 0 && !group.Contains(id) && nodenum[id] != -1)
		  {
		    group.Append(id);
		    id_num[id] = j+1+id_num[group[k]];
		  }
	      }
	  }
	if(group.Size() > 1)
	  {
	    id_groups.Append(new ARRAY<int>(group));
	    if(group.Size() == 2)
	      {
		id_type[i] = 1;
		id_type[group[1]] = 2;
		n2++;
	      }
	    else if(group.Size() == 4)
	      {
		id_type[i] = 3;
		for(int j=1; j<group.Size(); j++)
		  id_type[group[j]] = 4;
		n4++;
	      }
	    else if(group.Size() == 8)
	      {
		id_type[i] = 5;
		for(int j=1; j<group.Size(); j++)
		  id_type[group[j]] = 6;
		n8++;
	      }
	    else
	      cerr << "ERROR: Identification group size = " << group.Size() << endl;
	  }
	
      }


    for(PointIndex i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
      {
	if(nodenum[i] == -1)
	  continue;
	outfile << nodenum[i] << " "
		<< mesh[i](0) << " "
		<< mesh[i](1) << " "
		<< mesh[i](2) << " " << id_type[i] << " ";
	if(i-PointIndex::BASE < point_ids.Size())
	  outfile << point_ids[i];
	else
	  outfile << "0";
	outfile << "\n";
      }
    outfile << endl;

    outfile << "\n// Number of Periodic Master Nodes:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << n2 << "\n"			       \
	    << "\n" \
	    << "// MasterNodeID, SlaveNodeID, TranslCode (1=dS1 2=dS2 3=dS1+dS2):\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<id_groups.Size(); i++)
      {
	if(id_groups[i]->Size() != 2)
	  continue;

	for(int j=0; j<id_groups[i]->Size(); j++)
	  outfile << nodenum[(*id_groups[i])[j]] << " ";
	for(int j=1; j<id_groups[i]->Size(); j++)
	  outfile << id_num[(*id_groups[i])[j]] << " ";
	outfile << "\n";

	delete id_groups[i];
	id_groups[i] = NULL;
      }
    outfile << endl;
	
	
    outfile << "// Number of Corner Periodic Master Nodes:\n"	      \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << n4 << "\n"				      \
	    << "\n" \
	    << "// MasterNodeID, 3-SlaveNodeID's, 3-TranslCodes (1=dS1 2=dS2 3=dS1+dS2):\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";


    for(int i=0; i<id_groups.Size(); i++)
      {
	if(!id_groups[i] || id_groups[i]->Size() != 4)
	  continue;

	for(int j=0; j<id_groups[i]->Size(); j++)
	  outfile << nodenum[(*id_groups[i])[j]] << " ";
	for(int j=1; j<id_groups[i]->Size(); j++)
	  {
	    outfile << id_num[(*id_groups[i])[j]] << " ";
	  }
	outfile << "\n";

	delete id_groups[i];
	id_groups[i] = NULL;
      }
    outfile << endl;


    outfile << "// Number of Cubic Periodic Master Nodes:\n"	     \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << n8 << "\n"				     \
	    << "\n" \
	    << "// MasterNodeID, 7-SlaveNodeID's, TranslCodes:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<id_groups.Size(); i++)
      {
	if(!id_groups[i] || id_groups[i]->Size() != 8)
	  continue;

	for(int j=0; j<id_groups[i]->Size(); j++)
	  outfile << nodenum[(*id_groups[i])[j]] << " ";
	for(int j=1; j<id_groups[i]->Size(); j++)
	  outfile << id_num[(*id_groups[i])[j]] << " ";
	outfile << "\n";

	delete id_groups[i];
	id_groups[i] = NULL;
      }
    outfile << endl;

    


    outfile << "// EdgeID, NodeID0, NodeID1, Type (0=Reg 1=PMaster 2=PSlave 3=CPMaster 4=CPSlave), "<<uidpid<<":\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";

    
      
    ARRAY< ARRAY<int>* > vertex_to_edge(mesh.GetNP()+1);
    for(int i=0; i<=mesh.GetNP(); i++)
      vertex_to_edge[i] = new ARRAY<int>;

    ARRAY< ARRAY<int,PointIndex::BASE>* > idmaps_edge(idmaps.Size());
    for(int i=0; i<idmaps_edge.Size(); i++)
      {
	idmaps_edge[i] = new ARRAY<int,PointIndex::BASE>(numedges);
	(*idmaps_edge[i]) = 0;
      }

    ARRAY<int> possible;
    for(int i=0; i<edge2node.Size(); i++)
      {
	const INDEX_2 & v = edge2node[i];
	for(int j=0; j<idmaps.Size(); j++)
	  {
	    INDEX_2 vid((*idmaps[j])[v[0]], (*idmaps[j])[v[1]]);
	    if(vid[0] != 0 && vid[0] != v[0] && vid[1] != 0 && vid[1] != v[1])
	      {
		Intersection(*vertex_to_edge[vid[0]],*vertex_to_edge[vid[1]],possible);
		if(possible.Size() == 1)
		  {
		    (*idmaps_edge[j])[possible[0]] = i+1;
		    (*idmaps_edge[j])[i+1] = possible[0];
		  }
		else if(possible.Size() > 0)
		  {
		    cerr << "ERROR: too many possible edge identifications" << endl;
		    (*testout) << "ERROR: too many possible edge identifications" << endl
			       << "*vertex_to_edge["<<vid[0]<<"] " << *vertex_to_edge[vid[0]] << endl
			       << "*vertex_to_edge["<<vid[1]<<"] " << *vertex_to_edge[vid[1]] << endl
			       << "possible " << possible << endl;
		  }
	      }
	  }
	vertex_to_edge[v[0]]->Append(i+1);
	vertex_to_edge[v[1]]->Append(i+1);
      }


    for(int i=0; i<vertex_to_edge.Size(); i++)
      delete vertex_to_edge[i];


    id_groups.SetSize(0);
    id_num.SetSize(numedges+1);
    id_num = 0;
    id_type.SetSize(numedges+1);
    id_type = 0;

    n2 = n4 = n8 = 0;

    for(int i=1; i<=edge2node.Size(); i++)
      {
	if(id_num[i] != 0)
	  continue;


	ARRAY<int> group;
	group.Append(i);
	for(int j=0; j<idmaps_edge.Size(); j++)
	  {
	    startsize = group.Size();
	    for(int k=0; k<startsize; k++)
	      {
		int id = (*idmaps_edge[j])[group[k]];
		if(id != 0 && !group.Contains(id))
		  {
		    group.Append(id);
		    id_num[id] = j+1+id_num[group[k]];
		  }
	      }
	  }
	if(group.Size() > 1)
	  {
	    id_num[i] = 1;
	    id_groups.Append(new ARRAY<int>(group));
	    if(group.Size() == 2)
	      {
		id_type[i] = 1;
		id_type[group[1]] = 2;
		n2++;
	      }
	    else if(group.Size() == 4)
	      {
		id_type[i] = 3;
		for(int j=1; j<group.Size(); j++)
		  id_type[group[j]] = 4;
		n4++;
	      }
	    else
	      {
		cerr << "ERROR: edge identification group size = " << group.Size() << endl;
		(*testout) << "edge group " << group << endl;
		for(int j=0; j<idmaps_edge.Size(); j++)
		  {
		    (*testout) << "edge id map " << j << endl << *idmaps_edge[j] << endl;
		  }
	      }
	  }
      }



    for(int i=1; i<=edge2node.Size(); i++)
      {
	if(id_num[i] != 0)
	  continue;


	ARRAY<int> group;
	group.Append(i);
	for(int j=0; j<idmaps_edge.Size(); j++)
	  {
	    startsize = group.Size();
	    for(int k=0; k<startsize; k++)
	      {
		int id = (*idmaps_edge[j])[group[k]];
		if(id != 0 && !group.Contains(id))
		  {
		    group.Append(id);
		    id_num[id] = j+1+id_num[group[k]];
		  }
	      }
	  }
	if(group.Size() > 1)
	  {
	    id_num[i] = 1;
	    id_groups.Append(new ARRAY<int>(group));
	    if(group.Size() == 2)
	      {
		id_type[i] = 1;
		id_type[group[1]] = 2;
		n2++;
	      }
	    else if(group.Size() == 4)
	      {
		id_type[i] = 3;
		for(int j=1; j<group.Size(); j++)
		  id_type[group[j]] = 4;
		n4++;
	      }
	    else
	      {
		cerr << "ERROR: edge identification group size = " << group.Size() << endl;
		(*testout) << "edge group " << group << endl;
		for(int j=0; j<idmaps_edge.Size(); j++)
		  {
		    (*testout) << "edge id map " << j << endl << *idmaps_edge[j] << endl;
		  }
	      }
	  }
	
      }

    
    for(int i=0; i<edge2node.Size(); i++)
      outfile << i+1 << " " << nodenum[edge2node[i][0]] << " " << nodenum[edge2node[i][1]] 
	      << " " << id_type[i+1] << " " << edge_ids[i] << "\n";

    outfile << endl;

    

    outfile << "// Number of Periodic Master Edges:\n"\
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"\
	    << n2 << "\n"			      \
	    << "\n"\
	    << "// MasterEdgeID, SlaveEdgeID, TranslCode (1=dS1 2=dS2 3=dS1+dS2):\n"\
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<id_groups.Size(); i++)
      {
	if(id_groups[i]->Size() != 2)
	  continue;

	for(int j=0; j<id_groups[i]->Size(); j++)
	  outfile << (*id_groups[i])[j] << " ";
	for(int j=1; j<id_groups[i]->Size(); j++)
	  outfile << id_num[(*id_groups[i])[j]] << " ";
	outfile << "\n";

	delete id_groups[i];
	id_groups[i] = NULL;
      }
    outfile << endl;

    outfile << "// Number of Corner Periodic Master Edges:\n"		\
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"\
	    << n4 << "\n"				     \
	    << "\n"\
	    << "// MasterEdgeID, 3 SlaveEdgeID's, 3 TranslCode (1=dS1 2=dS2 3=dS1+dS2):\n"\
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<id_groups.Size(); i++)
      {
	if(!id_groups[i] || id_groups[i]->Size() != 4)
	  continue;

	for(int j=0; j<id_groups[i]->Size(); j++)
	  outfile << (*id_groups[i])[j] << " ";
	for(int j=1; j<id_groups[i]->Size(); j++)
	  outfile << id_num[(*id_groups[i])[j]] << " ";
	outfile << "\n";

	delete id_groups[i];
	id_groups[i] = NULL;
      }
    outfile << endl;


    outfile << "// FaceID, EdgeID0, EdgeID1, EdgeID2, FaceType (0=Reg 1=PMaster 2=PSlave), "<<uidpid<<":\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";

    
    
    ARRAY< ARRAY<int>* > edge_to_face(numedges+1);
    for(int i=0; i<edge_to_face.Size(); i++)
      edge_to_face[i] = new ARRAY<int>;

    
    for(int i=0; i<idmaps.Size(); i++)
      {
	idmaps[i]->SetSize(numfaces);
	(*idmaps[i]) = 0;
      }

    
    for(int i=0; i<face2edge.Size(); i++)
      {
	for(int j=0; j<idmaps_edge.Size(); j++)
	  {
	    int e1id,e2id,e3id;
	    e1id = (*idmaps_edge[j])[abs(face2edge[i][0])];
	    e2id = (*idmaps_edge[j])[abs(face2edge[i][1])];
	    e3id = (*idmaps_edge[j])[abs(face2edge[i][2])];
	    if(e1id != 0 && e1id != abs(face2edge[i][0]) &&
	       e2id != 0 && e2id != abs(face2edge[i][1]) &&
	       e3id != 0 && e3id != abs(face2edge[i][2]))
	      {
		Intersection(*edge_to_face[e1id],*edge_to_face[e2id],*edge_to_face[e3id],possible);
		if(possible.Size() == 1)
		  {
		    (*idmaps[j])[possible[0]] = i+1;
		    (*idmaps[j])[i+1] = possible[0];
		  }
		else if(possible.Size() > 0)
		  cerr << "ERROR: too many possible face identifications" << endl;
	      }
	  }

	edge_to_face[abs(face2edge[i][0])]->Append(i+1);
	edge_to_face[abs(face2edge[i][1])]->Append(i+1);
	edge_to_face[abs(face2edge[i][2])]->Append(i+1);
      }

    for(int i=0; i<edge_to_face.Size(); i++)
      delete edge_to_face[i];


    for(int i=0; i<idmaps_edge.Size(); i++)
      delete idmaps_edge[i];

    
    id_groups.SetSize(0);
    id_num.SetSize(numfaces+1);
    id_num = 0;

    n2 = n4 = n8 = 0;

    for(int i=1; i<=numfaces; i++)
      {
	if(id_num[i] != 0)
	  continue;

	ARRAY<int> group;
	group.Append(i);
	for(int j=0; j<idmaps.Size(); j++)
	  {
	    startsize = group.Size();
	    for(int k=0; k<startsize; k++)
	      {
		int id = (*idmaps[j])[group[k]];
		if(id != 0 && !group.Contains(id))
		  {
		    group.Append(id);
		    id_num[id] = j+1+id_num[group[k]];
		  }
	      }
	  }
	if(group.Size() > 1)
	  {
	    id_num[i] = -1;
	    id_groups.Append(new ARRAY<int>(group));
	    if(group.Size() == 2)
	      n2++;
	    else
	      cerr << "ERROR: face identification group size = " << group.Size() << endl;
	  }
	
      }


    for(int i=0; i<idmaps.Size(); i++)
      delete idmaps[i];




    for(int i=0; i<face2edge.Size(); i++)
      {	
	outfile << i+1 << " ";
	for(int j=0; j<3; j++)
	  outfile << face2edge[i][j] << " ";

	if(id_num[i+1] == 0)
	  outfile << 0;
	else if(id_num[i+1] == -1)
	  outfile << 1;
	else
	  outfile << 2;

	outfile << " " << face_ids[i] <<"\n";
      }
    outfile << endl;


    outfile << "// Number of Periodic Master Faces:\n"\
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"\
	    << n2 << "\n"			      \
	    << "\n"\
	    << "// MasterFaceID, SlaveFaceID, TranslCode (1=dS1 2=dS2):\n"\
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<id_groups.Size(); i++)
      {
	if(id_groups[i]->Size() != 2)
	  continue;

	for(int j=0; j<id_groups[i]->Size(); j++)
	  outfile << (*id_groups[i])[j] << " ";
	for(int j=1; j<id_groups[i]->Size(); j++)
	  outfile << id_num[(*id_groups[i])[j]] << " ";
	outfile << "\n";

	delete id_groups[i];
      }
    outfile << endl;

    


    outfile << "// ElemID, FaceID0, FaceID1, FaceID2, FaceID3, "<<uidpid<<":\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";

    for(ElementIndex i=0; i<mesh.GetNE(); i++)
      {
	if(elnum[i] >= 0)
	  {
	    outfile << elnum[i] << " ";
	    for(int j=0; j<4; j++)
	      outfile << element2face[elnum[i]-1][j] << " ";

	    outfile << mesh[i].GetIndex() << "\n";
	  }
      }
    outfile << endl;

    outfile << "// ElemID, NodeID0, NodeID1, NodeID2, NodeID3:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";

    
    for(ElementIndex i=0; i<mesh.GetNE(); i++)
      {
	if(elnum[i] >= 0)
	  outfile << elnum[i] << " "
		  << nodenum[mesh[i][1]] << " " << nodenum[mesh[i][0]] << " " << nodenum[mesh[i][2]] << " " << nodenum[mesh[i][3]] << "\n";
      }
    outfile << endl;
    

    

    outfile << "// Physical Object counts (#Obj3D,#Obj2D,#Obj1D,#Obj0D):\n"
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"
	    << " "<< numObj3D << " " << numObj2D << " " << numObj1D << " " << numObj0D << "\n" \
	    << "\n" \
	    << "// Number of Ports (Ports are a subset of Object2D list):\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" \
	    << numports << "\n"						\
	    << endl;


    ARRAY< ARRAY<int> * > groups;

    int maxg = -1;
    for(int i = 0; i<uid_to_group_3D.Size(); i++)
      if(uid_to_group_3D[i] > maxg)
	maxg = uid_to_group_3D[i];
    for(int i = 0; i<uid_to_group_2D.Size(); i++)
      if(uid_to_group_2D[i] > maxg)
	maxg = uid_to_group_2D[i];
    for(int i = 0; i<uid_to_group_1D.Size(); i++)
      if(uid_to_group_1D[i] > maxg)
	maxg = uid_to_group_1D[i];
    for(int i = 0; i<uid_to_group_0D.Size(); i++)
      if(uid_to_group_0D[i] > maxg)
	maxg = uid_to_group_0D[i];

    groups.SetSize(maxg+1);
    for(int i=0; i<groups.Size(); i++)
      groups[i] = new ARRAY<int>;

    for(ElementIndex i=0; i<mesh.GetNE(); i++)
      if(uid_to_group_3D[mesh[i].GetIndex()] >= 0)
	groups[uid_to_group_3D[mesh[i].GetIndex()]]->Append(i+1);
      
    


    outfile << "// Object3D GroupID, #Elems <immediately followed by> ElemID List:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<numObj3D; i++)
      {
	outfile << i << " " << groups[i]->Size() << "\n";
	for(int j=0; j<groups[i]->Size(); j++)
	  outfile << (*groups[i])[j] << "\n";
      }

    for(int i=0; i<groups.Size(); i++)
      groups[i]->SetSize(0);

    for(int i=0; i<face_ids.Size(); i++)
      if(uid_to_group_2D[face_ids[i]] >= 0)
	groups[uid_to_group_2D[face_ids[i]]]->Append(i+1);
      

    outfile << "// Object2D GroupID, #Faces <immediately followed by> FaceID List:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<numObj2D; i++)
      {
	outfile << i << " " << groups[i]->Size() << "\n";
	for(int j=0; j<groups[i]->Size(); j++)
	  {
	    outfile << (*groups[i])[j];
	    if(ports.Contains(face_ids[(*groups[i])[j]-1]))
	      outfile << " P";
	    outfile << "\n";
	  }
      }
    outfile << endl;

    
    for(int i=0; i<groups.Size(); i++)
      groups[i]->SetSize(0);

    for(int i=0; i<edge_ids.Size(); i++)
      if(uid_to_group_1D[edge_ids[i]] >= 0)
	groups[uid_to_group_1D[edge_ids[i]]]->Append(i+1);



    outfile << "// Object1D GroupID, #Edges <immediately followed by> EdgeID List:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<numObj1D; i++)
      {
	outfile << i << " " << groups[i]->Size() << "\n";
	for(int j=0; j<groups[i]->Size(); j++)
	  outfile << (*groups[i])[j] << "\n";
      }
    outfile << endl;

    
    for(int i=0; i<groups.Size(); i++)
      groups[i]->SetSize(0);
    for(PointIndex i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
      {
	if(i-PointIndex::BASE < point_ids.Size())
	  {
	    if(uid_to_group_0D[point_ids[i]] >= 0)
	      groups[uid_to_group_0D[point_ids[i]]]->Append(i+1-PointIndex::BASE);
	  }
	else
	  groups[uid_to_group_0D[0]]->Append(i+1-PointIndex::BASE);
      }


    outfile << "// Object0D GroupID, #Nodes <immediately followed by> NodeID List:\n" \
	    << "// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    for(int i=0; i<numObj0D; i++)
      {
	outfile << i << " " << groups[i]->Size() << "\n";
	for(int j=0; j<groups[i]->Size(); j++)
	  outfile << (*groups[i])[j] << "\n";
      }
    outfile << endl;

    for(int i=0; i<groups.Size(); i++)
      delete groups[i];


    outfile.close();

    cout << ".tet export done" << endl;
  }
}
