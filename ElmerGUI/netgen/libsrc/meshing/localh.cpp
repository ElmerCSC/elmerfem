#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{

GradingBox :: GradingBox (const double * ax1, const double * ax2)
{
  h2 = 0.5 * (ax2[0] - ax1[0]);
  for (int i = 0; i <= 2; i++)
    {
      /*
      x1[i] = ax1[i];
      x2[i] = ax2[i];
      */
      xmid[i] = 0.5 * (ax1[i] + ax2[i]);
    }

  /*
  (*testout) << "new box: " << xmid[0] << "-" << xmid[1] << "-" << xmid[2] 
	     << " h = " << (x2[0] - x1[0]) << endl;
  */

  for (int i = 0; i < 8; i++)
    childs[i] = NULL;
  father = NULL;

  flags.cutboundary = 0;
  flags.isinner = 0;
  flags.oldcell = 0;
  flags.pinner = 0;

  //  hopt = x2[0] - x1[0];
  hopt = 2 * h2;
}



BlockAllocator GradingBox :: ball(sizeof (GradingBox));

void * GradingBox :: operator new(size_t)
{
  return ball.Alloc();
}

void GradingBox :: operator delete (void * p)
{
  ball.Free (p);
}







void GradingBox :: DeleteChilds()
{
  int i;
  for (i = 0; i < 8; i++)
    if (childs[i])
      {
	childs[i]->DeleteChilds();
	delete childs[i];
	childs[i] = NULL;
      }
}


LocalH :: LocalH (const Point3d & pmin, const Point3d & pmax, double agrading)
{
  double x1[3], x2[3];
  double hmax;
  int i;

  boundingbox = Box3d (pmin, pmax);
  grading = agrading;

  // a small enlargement, non-regular points 
  double val = 0.0879;
  for (i = 1; i <= 3; i++)
    {

      x1[i-1] = (1 + val * i) * pmin.X(i) - val * i * pmax.X(i);
      x2[i-1] = 1.1 * pmax.X(i) - 0.1 * pmin.X(i);
    }

  hmax = x2[0] - x1[0];
  for (i = 1; i <= 2; i++)
    if (x2[i] - x1[i] > hmax)
      hmax = x2[i] - x1[i];

  for (i = 0; i <= 2; i++)
    x2[i] = x1[i] + hmax;

  root = new GradingBox (x1, x2);
  boxes.Append (root);
}

LocalH :: ~LocalH ()
{
  root->DeleteChilds();
  delete root;
}

void LocalH :: Delete ()
{
  root->DeleteChilds();
}

void LocalH :: SetH (const Point3d & p, double h)
{
  /*
  (*testout) << "Set h at " << p << " to " << h << endl;
  if (h < 1e-8)
    {
      cout << "do not set h to " << h << endl;
      return;
    }
  */

  if (fabs (p.X() - root->xmid[0]) > root->h2 ||
      fabs (p.Y() - root->xmid[1]) > root->h2 ||
      fabs (p.Z() - root->xmid[2]) > root->h2)
    return;

  /*      
  if (p.X() < root->x1[0] || p.X() > root->x2[0] ||
      p.Y() < root->x1[1] || p.Y() > root->x2[1] ||
      p.Z() < root->x1[2] || p.Z() > root->x2[2])
    return;
  */


  if (GetH(p) <= 1.2 * h) return;


  GradingBox * box = root;
  GradingBox * nbox = root;
  GradingBox * ngb;
  int childnr;
  double x1[3], x2[3];

  while (nbox)
    {
      box = nbox;
      childnr = 0;
      if (p.X() > box->xmid[0]) childnr += 1;
      if (p.Y() > box->xmid[1]) childnr += 2;
      if (p.Z() > box->xmid[2]) childnr += 4;
      nbox = box->childs[childnr];
    };


  while (2 * box->h2 > h)
    {
      childnr = 0;
      if (p.X() > box->xmid[0]) childnr += 1;
      if (p.Y() > box->xmid[1]) childnr += 2;
      if (p.Z() > box->xmid[2]) childnr += 4;

      double h2 = box->h2;
      if (childnr & 1)
	{
	  x1[0] = box->xmid[0];
	  x2[0] = x1[0]+h2;   // box->x2[0];
	}
      else
	{
	  x2[0] = box->xmid[0];
	  x1[0] = x2[0]-h2;   // box->x1[0];
	}

      if (childnr & 2)
	{
	  x1[1] = box->xmid[1];
	  x2[1] = x1[1]+h2;   // box->x2[1];
	}
      else
	{
	  x2[1] = box->xmid[1];
	  x1[1] = x2[1]-h2;   // box->x1[1];
	}

      if (childnr & 4)
	{
	  x1[2] = box->xmid[2];
	  x2[2] = x1[2]+h2;  // box->x2[2];
	}
      else
	{
	  x2[2] = box->xmid[2];
	  x1[2] = x2[2]-h2;  // box->x1[2];
	}

      ngb = new GradingBox (x1, x2);
      box->childs[childnr] = ngb;
      ngb->father = box;

      boxes.Append (ngb);
      box = box->childs[childnr];
    }

  box->hopt = h;


  double hbox = 2 * box->h2;  // box->x2[0] - box->x1[0];
  double hnp = h + grading * hbox;

  Point3d np;
  int i;
  for (i = 1; i <= 3; i++)
    {
      np = p;
      np.X(i) = p.X(i) + hbox;
      SetH (np, hnp);

      np.X(i) = p.X(i) - hbox;
      SetH (np, hnp);
    }
  /*
  Point3d np;
  int i1, i2, i3;
  for (i1 = -1; i1 <= 1; i1++)
    for (i2 = -1; i2 <= 1; i2++)
      for (i3 = -1; i3 <= 1; i3++)
	{
	  np.X() = p.X() + hbox * i1;
	  np.Y() = p.Y() + hbox * i2;
	  np.Z() = p.Z() + hbox * i3;

	  SetH (np, hnp);
	}
  */
}



double LocalH :: GetH (const Point3d & x) const
{
  const GradingBox * box = root;
  const GradingBox * nbox;
  int childnr;

  while (1)
    {
      childnr = 0;
      if (x.X() > box->xmid[0]) childnr += 1;
      if (x.Y() > box->xmid[1]) childnr += 2;
      if (x.Z() > box->xmid[2]) childnr += 4;
      nbox = box->childs[childnr];
      if (nbox)
	box = nbox;
      else
	{
	  //	  (*testout) << "diam = " << (box->x2[0] - box->x1[0])
	  //		     << " h = " << box->hopt << endl;
	  return box->hopt;
	}
    }
}


/// minimal h in box (pmin, pmax)
double LocalH :: GetMinH (const Point3d & pmin, const Point3d & pmax) const
{ 
  Point3d pmin2, pmax2;
  for (int j = 1; j <= 3; j++)
    if (pmin.X(j) < pmax.X(j))
      { pmin2.X(j) = pmin.X(j); pmax2.X(j) = pmax.X(j); }
    else
      { pmin2.X(j) = pmax.X(j); pmax2.X(j) = pmin.X(j); }

  return GetMinHRec (pmin2, pmax2, root); 
}


double LocalH :: GetMinHRec (const Point3d & pmin, const Point3d & pmax,
			     const GradingBox * box) const
{
  double h2 = box->h2;
  if (pmax.X() < box->xmid[0]-h2 || pmin.X() > box->xmid[0]+h2 ||
      pmax.Y() < box->xmid[1]-h2 || pmin.Y() > box->xmid[1]+h2 ||
      pmax.Z() < box->xmid[2]-h2 || pmin.Z() > box->xmid[2]+h2)
    return 1e8;
  /*
  if (pmax.X() < box->x1[0] || pmin.X() > box->x2[0] ||
      pmax.Y() < box->x1[1] || pmin.Y() > box->x2[1] ||
      pmax.Z() < box->x1[2] || pmin.Z() > box->x2[2])
    return 1e8;
  */

      
  double hmin = 2 * box->h2; // box->x2[0] - box->x1[0];
  int i;
  
  for (i = 0; i <= 7; i++)
    {
      if (box->childs[i])
	{
	  double hi = GetMinHRec (pmin, pmax, box->childs[i]);
	  if (hi < hmin)
	    hmin = hi;
	}	  
    }

  return hmin;
}


void LocalH :: CutBoundaryRec (const Point3d & pmin, const Point3d & pmax,
			       GradingBox * box)
{
  double h2 = box->h2;
  if (pmax.X() < box->xmid[0]-h2 || pmin.X() > box->xmid[0]+h2 ||
      pmax.Y() < box->xmid[1]-h2 || pmin.Y() > box->xmid[1]+h2 ||
      pmax.Z() < box->xmid[2]-h2 || pmin.Z() > box->xmid[2]+h2)
    return;
  /*
  if (pmax.X() < box->x1[0] || pmin.X() > box->x2[0] ||
      pmax.Y() < box->x1[1] || pmin.Y() > box->x2[1] ||
      pmax.Z() < box->x1[2] || pmin.Z() > box->x2[2])
    return;
  */

  box->flags.cutboundary = 1;
  for (int i = 0; i < 8; i++)
    if (box->childs[i])
      CutBoundaryRec (pmin, pmax, box->childs[i]);
}




void LocalH :: FindInnerBoxes ( // int (*sameside)(const Point3d & p1, const Point3d & p2),
			       AdFront3 * adfront,
			       int (*testinner)(const Point3d & p1))
{
  int i;

  int nf = adfront->GetNF();

  for (i = 0; i < boxes.Size(); i++)
    boxes[i] -> flags.isinner = 0;

  root->flags.isinner = 0;

  Point3d rpmid(root->xmid[0], root->xmid[1], root->xmid[2]);
  Vec3d rv(root->h2, root->h2, root->h2);
  Point3d rx2 = rpmid + rv;
  Point3d rx1 = rpmid - rv;


  root->flags.pinner = !adfront->SameSide (rpmid, rx2);
  
  if (testinner)
    (*testout) << "inner = " << root->flags.pinner << " =?= " 
	       << testinner(Point3d(root->xmid[0], root->xmid[1], root->xmid[2])) << endl;
  
  ARRAY<int> faceinds(nf);
  ARRAY<Box3d> faceboxes(nf);

  for (i = 1; i <= nf; i++)
    {
      faceinds.Elem(i) = i;
      adfront->GetFaceBoundingBox(i, faceboxes.Elem(i));
    }
  
  for (i = 0; i < 8; i++)
    FindInnerBoxesRec2 (root->childs[i], adfront, faceboxes, faceinds, nf);
}


void LocalH :: 
FindInnerBoxesRec2 (GradingBox * box,
		    class AdFront3 * adfront, 
		    ARRAY<Box3d> & faceboxes,
		    ARRAY<int> & faceinds, int nfinbox)
{
  if (!box) return;
  
  int i, j;
  
  GradingBox * father = box -> father;
  
  Point3d c(box->xmid[0], box->xmid[1], box->xmid[2]);
  Vec3d v(box->h2, box->h2, box->h2);
  Box3d boxc(c-v, c+v);

  Point3d fc(father->xmid[0], father->xmid[1], father->xmid[2]);
  Vec3d fv(father->h2, father->h2, father->h2);
  Box3d fboxc(fc-fv, fc+fv);

  Box3d boxcfc(c,fc);


  static ARRAY<int> faceused;
  static ARRAY<int> faceused2;
  static ARRAY<int> facenotused;

  faceused.SetSize(0);
  facenotused.SetSize(0);
  faceused2.SetSize(0);

  for (j = 1; j <= nfinbox; j++)
    {
      //      adfront->GetFaceBoundingBox (faceinds.Get(j), facebox);
      const Box3d & facebox = faceboxes.Get(faceinds.Get(j));
  
      if (boxc.Intersect (facebox))
	faceused.Append(faceinds.Get(j));
      else
	facenotused.Append(faceinds.Get(j));

      if (boxcfc.Intersect (facebox))
	faceused2.Append (faceinds.Get(j));
    }
  
  for (j = 1; j <= faceused.Size(); j++)
    faceinds.Elem(j) = faceused.Get(j);
  for (j = 1; j <= facenotused.Size(); j++)
    faceinds.Elem(j+faceused.Size()) = facenotused.Get(j);

  
  if (!father->flags.cutboundary)
    {
      box->flags.isinner = father->flags.isinner;
      box->flags.pinner = father->flags.pinner;
    }
  else
    {
      Point3d cf(father->xmid[0], father->xmid[1], father->xmid[2]);
      
      if (father->flags.isinner)
	box->flags.pinner = 1;
      else
	{
	  if (adfront->SameSide (c, cf, &faceused2))
	    box->flags.pinner = father->flags.pinner;
	  else
	    box->flags.pinner = 1 - father->flags.pinner;
	}
      
      if (box->flags.cutboundary)
	box->flags.isinner = 0;
      else
	box->flags.isinner = box->flags.pinner;
    }

  int nf = faceused.Size();
  for (i = 0; i < 8; i++)
    FindInnerBoxesRec2 (box->childs[i], adfront, faceboxes, faceinds, nf);
}












/*
void LocalH :: FindInnerBoxes ( // int (*sameside)(const Point3d & p1, const Point3d & p2),
			       AdFront3 * adfront,
			       int (*testinner)(const Point3d & p1))
{
  int i;
  for (i = 1; i <= boxes.Size(); i++)
    boxes.Elem(i)->flags.isinner = 0;

  root->flags.isinner = 0;

  Point3d rpmid(root->xmid[0], root->xmid[1], root->xmid[2]);
  Point3d rx2 = rpmid + Vec3d (root->h2, root->h2, root->h2);

  root->flags.pinner = !adfront->SameSide (rpmid, rx2);
  
  if (testinner)
    (*testout) << "inner = " << root->flags.pinner << " =?= " 
	       << testinner(Point3d(root->xmid[0], root->xmid[1], root->xmid[2])) << endl;
  

  for (i = 2; i <= boxes.Size(); i++)
    {
      GradingBox * box = boxes.Elem(i);
      GradingBox * father = box -> father;

      Point3d c(box->xmid[0], box->xmid[1], box->xmid[2]);
      Vec3d v(box->h2, box->h2, box->h2);
      Point3d x1 = c-v;
      Point3d x2 = c+v;


      if (!father->flags.cutboundary)
	{
	  box->flags.isinner = father->flags.isinner;
	  box->flags.pinner = father->flags.pinner;
	}
      else
	{
	  Point3d cf(father->xmid[0], father->xmid[1], father->xmid[2]);

	  if (father->flags.isinner)
	    box->flags.pinner = 1;
	  else
	    {
	      if (adfront->SameSide (c, cf))
		box->flags.pinner = father->flags.pinner;
	      else
		box->flags.pinner = 1 - father->flags.pinner;
	    }

	  if (box->flags.cutboundary)
	    box->flags.isinner = 0;
	  else
	    box->flags.isinner = box->flags.pinner;
	}
    }
  //  FindInnerBoxesRec (inner, root);
}
*/


void LocalH :: FindInnerBoxesRec ( int (*inner)(const Point3d & p),
				   GradingBox * box)
{
  int i;
  if (box->flags.cutboundary)
    {
      for (i = 0; i < 8; i++)
	if (box->childs[i])
	  FindInnerBoxesRec (inner, box->childs[i]);
    }
  else
    {
      if (inner (Point3d (box->xmid[0], box->xmid[1], box->xmid[2])))
	SetInnerBoxesRec (box);
    }
}


void LocalH :: SetInnerBoxesRec (GradingBox * box)
{
  box->flags.isinner = 1;
  for (int i = 0; i < 8; i++)
    if (box->childs[i])
      ClearFlagsRec (box->childs[i]);
}

void LocalH :: ClearFlagsRec (GradingBox * box)
{
  box->flags.cutboundary = 0;
  box->flags.isinner = 0;
  for (int i = 0; i < 8; i++)
    if (box->childs[i])
      ClearFlagsRec (box->childs[i]);
}


void LocalH :: WidenRefinement ()
{
  int nb = boxes.Size(); 
  int i;
  //  (*testout) << "old boxes: " << nb << endl;
  for (i = 1; i <= nb; i++)
    {
      GradingBox * box = boxes.Get(i);
      //      double h = box->x2[0] - box->x1[0];
      double h = box->hopt;
      Point3d c(box->xmid[0], box->xmid[1], box->xmid[2]);
      //      (*testout) << " i = " << i 
      //		 << " c = " << c << " h = " << h << endl;

      for (int i1 = -1; i1 <= 1; i1++)
	for (int i2 = -1; i2 <= 1; i2++)
	  for (int i3 = -1; i3 <= 1; i3++)
	    SetH (Point3d (c.X() + i1 * h, 
			   c.Y() + i2 * h,
			   c.Z() + i3 * h), 1.001 * h);     
    }
}

void LocalH :: GetInnerPoints (ARRAY<Point3d> & points)
{
  int i, nb = boxes.Size(); 

  for (i = 1; i <= nb; i++)
    {
      GradingBox * box = boxes.Get(i);
      /*
      if (box->flags.pinner)
	points.Append (box->randomip);
      */
      //      if (box->flags.pinner)
      if (box->flags.isinner)
	{
	  Point3d c(box->xmid[0], box->xmid[1], box->xmid[2]);
	  points.Append (c);
	  /*
	  cout << "add point " << c << "; h = " << box->hopt
	       << "; max-min = " << (box->x2[0]-box->x1[0]) << endl;
	  */
	}
    }
}



void LocalH :: GetOuterPoints (ARRAY<Point3d> & points)
{
  int i, nb = boxes.Size(); 

  for (i = 1; i <= nb; i++)
    {
      GradingBox * box = boxes.Get(i);
      if (!box->flags.isinner &&
	  !box->flags.cutboundary)
	{
	  Point3d c(box->xmid[0], box->xmid[1], box->xmid[2]);
	  points.Append (c);
	}
    }
}



void LocalH :: Convexify ()
{
  ConvexifyRec (root);
}

void LocalH :: ConvexifyRec (GradingBox * box)
{
  Point3d center(box->xmid[0], box->xmid[1], box->xmid[2]);
  Point3d hp;

  double size = 2 * box->h2; // box->x2[0] - box->x1[0];
  double dx = 0.6 * size;

  double maxh = box->hopt;
  int i;

  
  
  for (i = 1; i <= 6; i++)
    {
      hp = center;
      switch (i)
	{
	case 1: hp.X() += dx; break;
	case 2: hp.X() -= dx; break;
	case 3: hp.Y() += dx; break;
	case 4: hp.Y() -= dx; break;
	case 5: hp.Z() += dx; break;
	case 6: hp.Z() -= dx; break;
	}
      
      double hh = GetH (hp);
      if (hh > maxh) maxh = hh;
    }

  if (maxh < 0.95 * box->hopt)
    SetH (center, maxh);

  for (i = 0; i < 8; i++)
    if (box->childs[i])
      ConvexifyRec (box->childs[i]);  
}

void LocalH :: PrintMemInfo (ostream & ost) const
{
  ost << "LocalH: " << boxes.Size() << " boxes of " << sizeof(GradingBox)
      << " bytes = " << boxes.Size()*sizeof(GradingBox) << " bytes" << endl;
}
}
