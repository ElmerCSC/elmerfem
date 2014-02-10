#include <mystdlib.h>

#include "meshing.hpp"
#include <opti.hpp>

namespace netgen
{

  class ImprovementRule
  {
  public:
    ARRAY<Element2d> oldels;
    ARRAY<Element2d> newels;
    ARRAY<INDEX_2> deledges;
    ARRAY<int> incelsonnode;
    ARRAY<int> reused;
    int bonus;
    int onp;
  };


  void MeshOptimize2d :: GenericImprove (Mesh & mesh)
  {
    if (!faceindex)
      {
	if (writestatus)
	  PrintMessage (3, "Generic Improve");

	for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	  GenericImprove (mesh);
      
	faceindex = 0;
      }

    // int j, k, l, ri;
    int np = mesh.GetNP();
    int ne = mesh.GetNSE();
    //    SurfaceElementIndex sei;

    
//     for (SurfaceElementIndex sei = 0; sei < ne; sei++)
//       {
// 	const Element2d & el = mesh[sei];
// 	(*testout) << "element " << sei << ": " <<flush;
// 	for(int j=0; j<el.GetNP(); j++)
// 	  (*testout) << el[j] << " " << flush;
// 	(*testout) << "IsDeleted() " << el.IsDeleted()<< endl;
//       }

    bool ok;
    int olddef, newdef;

    ARRAY<ImprovementRule*> rules;
    ARRAY<SurfaceElementIndex> elmap;
    ARRAY<int> elrot;
    ARRAY<PointIndex> pmap;
    ARRAY<PointGeomInfo> pgi;

    int surfnr = mesh.GetFaceDescriptor (faceindex).SurfNr();
  
    ImprovementRule * r1;

    // 2 triangles to quad
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3));
    r1->oldels.Append (Element2d (3, 2, 4));
    r1->newels.Append (Element2d (1, 2, 4, 3));
    r1->deledges.Append (INDEX_2 (2,3));
    r1->onp = 4;
    r1->bonus = 2;
    rules.Append (r1);

    // 2 quad to 1 quad
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (4, 3, 2, 5));
    r1->newels.Append (Element2d (1, 2, 5, 4));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->deledges.Append (INDEX_2 (3, 4));
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);

    // swap quads
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (3, 2, 5, 6));
    r1->newels.Append (Element2d (1, 6, 3, 4));
    r1->newels.Append (Element2d (1, 2, 5, 6));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->onp = 6;
    r1->bonus = 0;
    rules.Append (r1);

    // three quads to 2
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 5, 6, 3));
    r1->oldels.Append (Element2d (3, 6, 7, 4));
    r1->newels.Append (Element2d (1, 2, 5, 4));
    r1->newels.Append (Element2d (4, 5, 6, 7));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->deledges.Append (INDEX_2 (3, 4));
    r1->deledges.Append (INDEX_2 (3, 6));
    r1->onp = 7;
    r1->bonus = -1;
    rules.Append (r1);

    // quad + 2 connected trigs to quad
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 5, 3));
    r1->oldels.Append (Element2d (3, 5, 4));
    r1->newels.Append (Element2d (1, 2, 5, 4));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->deledges.Append (INDEX_2 (3, 4));
    r1->deledges.Append (INDEX_2 (3, 5));
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);

    // quad + 2 non-connected trigs to quad (a and b)
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 6, 3));
    r1->oldels.Append (Element2d (1, 4, 5));
    r1->newels.Append (Element2d (1, 3, 4, 5));
    r1->newels.Append (Element2d (1, 2, 6, 3));
    r1->deledges.Append (INDEX_2 (1, 4));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->onp = 6;
    r1->bonus = 0;
    rules.Append (r1);

    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 6, 3));
    r1->oldels.Append (Element2d (1, 4, 5));
    r1->newels.Append (Element2d (1, 2, 4, 5));
    r1->newels.Append (Element2d (4, 2, 6, 3));
    r1->deledges.Append (INDEX_2 (1, 4));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->onp = 6;
    r1->bonus = 0;
    rules.Append (r1);

    // two quad + trig -> one quad + trig
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 5, 6, 3));
    r1->oldels.Append (Element2d (4, 3, 6));
    r1->newels.Append (Element2d (1, 2, 6, 4));
    r1->newels.Append (Element2d (2, 5, 6));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->deledges.Append (INDEX_2 (3, 4));
    r1->deledges.Append (INDEX_2 (3, 6));
    r1->onp = 6;
    r1->bonus = -1;
    rules.Append (r1);

    // swap quad + trig (a and b)
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 5, 3));
    r1->newels.Append (Element2d (2, 5, 3, 4));
    r1->newels.Append (Element2d (1, 2, 4));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);

    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (2, 5, 3));
    r1->newels.Append (Element2d (1, 2, 5, 3));
    r1->newels.Append (Element2d (1, 3, 4));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);


    // 2 quads to quad + 2 trigs
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (1, 2, 3, 4));
    r1->oldels.Append (Element2d (3, 2, 5, 6));
    r1->newels.Append (Element2d (1, 5, 6, 4));
    r1->newels.Append (Element2d (1, 2, 5));
    r1->newels.Append (Element2d (4, 6, 3));
    r1->deledges.Append (INDEX_2 (2, 3));
    r1->onp = 6;
    r1->bonus = 0;
    //    rules.Append (r1);




    ARRAY<int> mapped(rules.Size());
    ARRAY<int> used(rules.Size());
    used = 0;
    mapped = 0;

  

    for (int ri = 0; ri < rules.Size(); ri++)
      {
	ImprovementRule & rule = *rules[ri];
	rule.incelsonnode.SetSize (rule.onp);
	rule.reused.SetSize (rule.onp);

	for (int j = 1; j <= rule.onp; j++)
	  {
	    rule.incelsonnode.Elem(j) = 0;
	    rule.reused.Elem(j) = 0;
	  }

	for (int j = 1; j <= rule.oldels.Size(); j++)
	  {
	    const Element2d & el = rule.oldels.Elem(j);
	    for (int k = 1; k <= el.GetNP(); k++)
	      rule.incelsonnode.Elem(el.PNum(k))--;
	  }

	for (int j = 1; j <= rule.newels.Size(); j++)
	  {
	    const Element2d & el = rule.newels.Elem(j);
	    for (int k = 1; k <= el.GetNP(); k++)
	      {
		rule.incelsonnode.Elem(el.PNum(k))++;
		rule.reused.Elem(el.PNum(k)) = 1;
	      }
	  }
      }



  
    TABLE<int,PointIndex::BASE> elonnode(np);
    ARRAY<int,PointIndex::BASE> nelonnode(np);
    TABLE<SurfaceElementIndex> nbels(ne);

    nelonnode = -4;

    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      {
	const Element2d & el = mesh[sei];

	if (el.GetIndex() == faceindex && !el.IsDeleted())
	  {
	    for (int j = 0; j < el.GetNP(); j++)
	      elonnode.Add (el[j], sei);
	  }
	if(!el.IsDeleted())
	  {
	    for (int j = 0; j < el.GetNP(); j++)
	      nelonnode[el[j]]++;
	  }
      }

    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      {
	const Element2d & el = mesh[sei];
	if (el.GetIndex() == faceindex && !el.IsDeleted())
	  {
	    for (int j = 0; j < el.GetNP(); j++)
	      {
		for (int k = 0; k < elonnode[el[j]].Size(); k++)
		  {
		    int nbel = elonnode[el[j]] [k];
		    bool inuse = false;
		    for (int l = 0; l < nbels[sei].Size(); l++)
		      if (nbels[sei][l] == nbel)
			inuse = true;
		    if (!inuse)
		      nbels.Add (sei, nbel);
		  }
	      }
	  }
      }


    for (int ri = 0; ri < rules.Size(); ri++)
      {
	const ImprovementRule & rule = *rules[ri];
      
	elmap.SetSize (rule.oldels.Size());
	elrot.SetSize (rule.oldels.Size());
	pmap.SetSize (rule.onp);
	pgi.SetSize (rule.onp);


	for (SurfaceElementIndex sei = 0; sei < ne; sei++)
	  {
	    if (multithread.terminate)
	      break;
	    if (mesh[sei].IsDeleted()) continue;

	    elmap[0] = sei;
	    FlatArray<SurfaceElementIndex> neighbours = nbels[sei];
	    
	    for (elrot[0] = 0; elrot[0] < mesh[sei].GetNP(); elrot[0]++)
	      {
		const Element2d & el0 = mesh[sei];
		const Element2d & rel0 = rule.oldels[0];

		if (el0.GetIndex() != faceindex) continue;
		if (el0.IsDeleted()) continue;
		if (el0.GetNP() != rel0.GetNP()) continue;


		pmap = -1;
 
		for (int k = 0; k < el0.GetNP(); k++)
		  {
		    pmap.Elem(rel0[k]) = el0.PNumMod(k+elrot[0]+1);
		    pgi.Elem(rel0[k]) = el0.GeomInfoPiMod(k+elrot[0]+1);
		  }
		
		ok = 1;
		for (int i = 1; i < elmap.Size(); i++)
		  {
		    // try to find a mapping for reference-element i

		    const Element2d & rel = rule.oldels[i];
		    bool possible = 0;

		    for (elmap[i] = 0; elmap[i] < neighbours.Size(); elmap[i]++)
		      {
			const Element2d & el = mesh[neighbours[elmap[i]]];
			if (el.IsDeleted()) continue;
			if (el.GetNP() != rel.GetNP()) continue;

			for (elrot[i] = 0; elrot[i] < rel.GetNP(); elrot[i]++)
			  {
			    possible = 1;

			    for (int k = 0; k < rel.GetNP(); k++)
			      if (pmap.Elem(rel[k]) != -1 &&
				  pmap.Elem(rel[k]) != el.PNumMod (k+elrot[i]+1))
				possible = 0;

			    if (possible) 
			      {
				for (int k = 0; k < el.GetNP(); k++)
				  {
				    pmap.Elem(rel[k]) = el.PNumMod(k+elrot[i]+1);
				    pgi.Elem(rel[k]) = el.GeomInfoPiMod(k+elrot[i]+1);
				  }
				break;
			      }
			  }
			if (possible) break;
		      }

		    if (!possible) 
		      {
			ok = 0;
			break;
		      }

		    elmap[i] = neighbours[elmap[i]];
		  }

		for(int i=0; ok && i<rule.deledges.Size(); i++)
		  {
		    ok = !mesh.IsSegment(pmap.Elem(rule.deledges[i].I1()),
					 pmap.Elem(rule.deledges[i].I2()));
		  }
								    
								    
		
		
		if (!ok) continue;

		mapped[ri]++;

		olddef = 0;
		for (int j = 1; j <= pmap.Size(); j++)
		  olddef += sqr (nelonnode[pmap.Get(j)]);
		olddef += rule.bonus;

		newdef = 0;
		for (int j = 1; j <= pmap.Size(); j++)
		  if (rule.reused.Get(j))
		    newdef += sqr (nelonnode[pmap.Get(j)] + 
				   rule.incelsonnode.Get(j));

		if (newdef > olddef)
		  continue;

		// calc metric badness
		double bad1 = 0, bad2 = 0;
		Vec<3> n;

		SelectSurfaceOfPoint (mesh.Point(pmap.Get(1)), pgi.Get(1));
		GetNormalVector (surfnr, mesh.Point(pmap.Get(1)), pgi.Elem(1), n);
		  
		for (int j = 1; j <= rule.oldels.Size(); j++)
		  bad1 += mesh.SurfaceElement(elmap.Get(j)).CalcJacobianBadness (mesh.Points(), n);
		  
		// check new element:
		for (int j = 1; j <= rule.newels.Size(); j++)
		  {
		    const Element2d & rnel = rule.newels.Get(j);
		    Element2d nel(rnel.GetNP());
		    for (int k = 1; k <= rnel.GetNP(); k++)
		      nel.PNum(k) = pmap.Get(rnel.PNum(k));

		    bad2 += nel.CalcJacobianBadness (mesh.Points(), n);
		  }

		if (bad2 > 1e3) continue;

		if (newdef == olddef && bad2 > bad1) continue;
		  

		// generate new element:
		for (int j = 1; j <= rule.newels.Size(); j++)
		  {
		    const Element2d & rnel = rule.newels.Get(j);
		    Element2d nel(rnel.GetNP());
		    nel.SetIndex (faceindex);
		    for (int k = 1; k <= rnel.GetNP(); k++)
		      {
			nel.PNum(k) = pmap.Get(rnel.PNum(k));
			nel.GeomInfoPi(k) = pgi.Get(rnel.PNum(k));
		      }
		      
		    mesh.AddSurfaceElement(nel);
		  }
		  
		for (int j = 0; j < rule.oldels.Size(); j++)
		  mesh.DeleteSurfaceElement ( elmap[j] );

		for (int j = 1; j <= pmap.Size(); j++)
		  nelonnode[pmap.Get(j)] += rule.incelsonnode.Get(j);

		used[ri]++;
	      }
	  }
      }

    mesh.Compress();

    for (int ri = 0; ri < rules.Size(); ri++)
      {
	PrintMessage (5, "rule ", ri+1, " ",
		      mapped[ri], "/", used[ri], " mapped/used");
      }
  }




}
