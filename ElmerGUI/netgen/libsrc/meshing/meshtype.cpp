#include <mystdlib.h>

#include "meshing.hpp"  

namespace netgen
{


  ostream & operator<<(ostream  & s, const MeshPoint & pt)
  {
    s << Point<3> (pt);
    return s;
  }

  /*
  MultiPointGeomInfo :: MultiPointGeomInfo()
  {
    cnt = 0;
  }
  */

  int MultiPointGeomInfo :: 
  AddPointGeomInfo (const PointGeomInfo & gi)
  {
    for (int k = 0; k < cnt; k++)
      if (mgi[k].trignum == gi.trignum)
	return 0;
  
    if (cnt < MULTIPOINTGEOMINFO_MAX)
      {
	mgi[cnt] = gi;
	cnt++;
	return 0;
      }

    throw NgException ("Please report error: MPGI Size too small\n");
  }
  

  /*
  void MultiPointGeomInfo :: 
  Init ()
  {
    cnt = 0;
  }

  void MultiPointGeomInfo :: 
  DeleteAll ()
  {
    cnt = 0;
  }
  */



  Segment :: Segment() 
  {
    p1 = -1;
    p2 = -1; 
    edgenr = -1;

    singedge_left = 0.;
    singedge_right = 0.;
    seginfo = 0;

    si = -1;

    domin = -1;
    domout = -1;
    tlosurf = -1;

    surfnr1 = -1;
    surfnr2 = -1;
    pmid = -1;
    meshdocval = 0;
    /*
    geominfo[0].trignum=-1; 
    geominfo[1].trignum=-1; 

    epgeominfo[0].edgenr = 1;
    epgeominfo[0].dist = 0;
    epgeominfo[1].edgenr = 1;
    epgeominfo[1].dist = 0;
    */

    bcname = 0;
  }    

  Segment::Segment (const Segment & other)
    : p1(other.p1),
      p2(other.p2),
      edgenr(other.edgenr),
      singedge_left(other.singedge_left),
      singedge_right(other.singedge_right),
      seginfo(other.seginfo),
      si(other.si),
      domin(other.domin),
      domout(other.domout),
      tlosurf(other.tlosurf),
      geominfo(),
      surfnr1(other.surfnr1),
      surfnr2(other.surfnr2),
      epgeominfo(),
      pmid(other.pmid),
      meshdocval(other.meshdocval),
      hp_elnr(other.hp_elnr)
  {
    geominfo[0] = other.geominfo[0];
    geominfo[1] = other.geominfo[1];
    epgeominfo[0] = other.epgeominfo[0];
    epgeominfo[1] = other.epgeominfo[1];
    bcname = other.bcname;
  }

  Segment& Segment::operator=(const Segment & other)
  {
    if (&other != this)
      {
	p1 = other.p1;
	p2 = other.p2;
	edgenr = other.edgenr;
	singedge_left = other.singedge_left;
	singedge_right = other.singedge_right;
	seginfo = other.seginfo;
	si = other.si;
	domin = other.domin;
	domout = other.domout;
	tlosurf = other.tlosurf;
	geominfo[0] = other.geominfo[0];
	geominfo[1] = other.geominfo[1];
	surfnr1 = other.surfnr1;
	surfnr2 = other.surfnr2;
	epgeominfo[0] = other.epgeominfo[0];
	epgeominfo[1] = other.epgeominfo[1];
	pmid = other.pmid;
	meshdocval = other.meshdocval;
	hp_elnr = other.hp_elnr;
	bcname = other.bcname;
      }
    
    return *this;
  }


  ostream & operator<<(ostream  & s, const Segment & seg)
  {
    s << seg.p1 << "(gi=" << seg.geominfo[0].trignum << ") - "
      << seg.p2 << "(gi=" << seg.geominfo[1].trignum << ")"
      << " domin = " << seg.domin << ", domout = " << seg.domout 
      << " si = " << seg.si << ", edgenr = " << seg.edgenr;
    return s;
  }


  Element2d :: Element2d ()
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }
    np = 3;
    index = 0;
    badel = 0;
    deleted = 0;
    typ = TRIG;
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
#ifdef PARALLEL
    isghost = 0;
#endif
  } 


  Element2d :: Element2d (int anp)
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }
    np = anp;
    index = 0;
    badel = 0;
    deleted = 0;
    switch (np)
      {
      case 3: typ = TRIG; break;
      case 4: typ = QUAD; break;
      case 6: typ = TRIG6; break;
      case 8: typ = QUAD8; break;
      }
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
#ifdef PARALLEL
    isghost = 0;
#endif
  } 

  Element2d :: Element2d (ELEMENT_TYPE atyp)
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }

    SetType (atyp);

    index = 0;
    badel = 0;
    deleted = 0;
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
#ifdef PARALLEL
  isghost = 0;
#endif

  } 



  Element2d :: Element2d (int pi1, int pi2, int pi3)
{
  pnum[0] = pi1;
  pnum[1] = pi2;
  pnum[2] = pi3;
  np = 3;
  typ = TRIG;
  pnum[3] = 0;
  pnum[4] = 0;
  pnum[5] = 0;
  
  for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
    geominfo[i].trignum = 0;
  index = 0;
  badel = 0;
  refflag = 1;
  strongrefflag = false;
  deleted = 0;
  orderx = ordery = 1;

#ifdef PARALLEL
  isghost = 0;
#endif

}

Element2d :: Element2d (int pi1, int pi2, int pi3, int pi4)
{
  pnum[0] = pi1;
  pnum[1] = pi2;
  pnum[2] = pi3;
  pnum[3] = pi4;
  np = 4;
  typ = QUAD;

  pnum[4] = 0;
  pnum[5] = 0;
  
  for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
    geominfo[i].trignum = 0;
  index = 0;
  badel = 0;
  refflag = 1;
  strongrefflag = false;
  deleted = 0;
  orderx = ordery = 1;

#ifdef PARALLEL
  isghost = 0;
#endif
}


/*
void Element2d :: SetType (ELEMENT_TYPE atyp)
{
  typ = atyp;
  switch (typ)
    {
    case TRIG: np = 3; break;
    case QUAD: np = 4; break;
    case TRIG6: np = 6; break;
    case QUAD6: np = 6; break;
    default:
      PrintSysError ("Element2d::SetType, illegal type ", typ);
    }
}
*/


void Element2d :: GetBox (const T_POINTS & points, Box3d & box) const
{
  box.SetPoint (points.Get(pnum[0]));
  for (unsigned i = 1; i < np; i++)
    box.AddPoint (points.Get(pnum[i]));
}

bool Element2d :: operator==(const Element2d & el2) const
{
  bool retval = (el2.GetNP() == np);
  for(int i= 0; retval && i<np; i++)
    retval = (el2[i] == (*this)[i]);

  return retval;
}


void Element2d :: Invert2()
{
  switch (typ)
    {
    case TRIG:
      {
	Swap (pnum[1], pnum[2]);
	break;
      }
    case QUAD:
      {
	Swap (pnum[0], pnum[3]);
	Swap (pnum[1], pnum[2]);
	break;
      }
    default:
      {
	cerr << "Element2d::Invert2, illegal element type " << int(typ) << endl;
      }
    }
}

int Element2d::HasFace(const Element2d& el) const
{
  //nur für tets!!! hannes
  for (int i = 1; i <= 3; i++)
    {
      if (PNumMod(i)   == el[0] && 
	  PNumMod(i+1) == el[1] && 
	  PNumMod(i+2) == el[2])
	{
	  return 1;
	}
    }
  return 0;
}

void Element2d :: NormalizeNumbering2 ()
{
  if (GetNP() == 3)
    {
      if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
	return;
      else
	{
	  if (PNum(2) < PNum(3))
	    {
	      PointIndex pi1 = PNum(2);
	      PNum(2) = PNum(3);
	      PNum(3) = PNum(1);
	      PNum(1) = pi1;
	    }
	  else
	    {
	      PointIndex pi1 = PNum(3);
	      PNum(3) = PNum(2);
	      PNum(2) = PNum(1);
	      PNum(1) = pi1;
	    }
	}
    }
  else
    {
      int mini = 1;
      for (int i = 2; i <= GetNP(); i++)
	if (PNum(i) < PNum(mini)) mini = i;
      
      Element2d hel = (*this);
      for (int i = 1; i <= GetNP(); i++)
	PNum(i) = hel.PNumMod (i+mini-1);
    }
}




ARRAY<IntegrationPointData*> ipdtrig;
ARRAY<IntegrationPointData*> ipdquad;


int Element2d :: GetNIP () const
{
  int nip;
  switch (np)
    {
    case 3: nip = 1; break;
    case 4: nip = 4; break;
    default: nip = 0; break;
    }
  return nip;
}

void Element2d :: 
GetIntegrationPoint (int ip, Point2d & p, double & weight) const
{
  static double eltriqp[1][3] =
  {
    { 1.0/3.0, 1.0/3.0, 0.5 }
  };

  static double elquadqp[4][3] =
  { 
    { 0, 0, 0.25 },
    { 0, 1, 0.25 },
    { 1, 0, 0.25 },
    { 1, 1, 0.25 }
  };
  
  double * pp = 0;
  switch (typ)
    {
    case TRIG: pp = &eltriqp[0][0]; break;
    case QUAD: pp = &elquadqp[ip-1][0]; break;
    default:
      PrintSysError ("Element2d::GetIntegrationPoint, illegal type ", typ);
    }

  p.X() = pp[0];
  p.Y() = pp[1];
  weight = pp[2];
}

void Element2d :: 
GetTransformation (int ip, const ARRAY<Point2d> & points,
		   DenseMatrix & trans) const
{
  int np = GetNP();
  static DenseMatrix pmat(2, np), dshape(2, np);
  pmat.SetSize (2, np);
  dshape.SetSize (2, np);

  Point2d p;
  double w;

  GetPointMatrix (points, pmat);
  GetIntegrationPoint (ip, p, w);
  GetDShape (p, dshape);
  
  CalcABt (pmat, dshape, trans);

  /*
  (*testout) << "p = " << p  << endl
	     << "pmat = " << pmat << endl
	     << "dshape = " << dshape << endl
	     << "tans = " << trans << endl;
  */
}

void Element2d :: 
GetTransformation (int ip, class DenseMatrix & pmat,
		   class DenseMatrix & trans) const
{
  int np = GetNP();

#ifdef DEBUG
  if (pmat.Width() != np || pmat.Height() != 2)
    {
      (*testout) << "GetTransofrmation: pmat doesn't fit" << endl;
      return;
    }
#endif

  ComputeIntegrationPointData ();
  DenseMatrix * dshapep;
  switch (typ)
    {
    case TRIG: dshapep = &ipdtrig.Get(ip)->dshape; break;
    case QUAD: dshapep = &ipdquad.Get(ip)->dshape; break;
    default:
      PrintSysError ("Element2d::GetTransformation, illegal type ", typ);
    }
  
  CalcABt (pmat, *dshapep, trans);
}



void Element2d :: GetShape (const Point2d & p, Vector & shape) const
{
  if (shape.Size() != GetNP())
    {
      cerr << "Element::GetShape: Length not fitting" << endl;
      return;
    }

  switch (typ)
    {
    case TRIG:
      shape.Elem(1) = 1 - p.X() - p.Y();
      shape.Elem(2) = p.X();
      shape.Elem(3) = p.Y();
      break;
    case QUAD:
      shape.Elem(1) = (1-p.X()) * (1-p.Y());
      shape.Elem(2) = p.X() * (1-p.Y());
      shape.Elem(3) = p.X() * p.Y();
      shape.Elem(4) = (1-p.X()) * p.Y();
      break;
    default:
      PrintSysError ("Element2d::GetShape, illegal type ", typ);
    }
}



void Element2d :: GetShapeNew (const Point<2> & p, FlatVector & shape) const
{
  switch (typ)
    {
    case TRIG:
      {
	shape(0) = p(0);
	shape(1) = p(1);
	shape(2) = 1-p(0)-p(1);
	break;
      }

    case QUAD:
      {
	shape(0) = (1-p(0))*(1-p(1));
	shape(1) =    p(0) *(1-p(1));
	shape(2) =    p(0) *   p(1) ;
	shape(3) = (1-p(0))*   p(1) ;
	break;
      }
    }
}









void Element2d :: 
GetDShape (const Point2d & p, DenseMatrix & dshape) const
{
#ifdef DEBUG
  if (dshape.Height() != 2 || dshape.Width() != np)
    {
      PrintSysError ("Element::DShape: Sizes don't fit");
      return;
    }
#endif

  switch (typ)
    {
    case TRIG:
      dshape.Elem(1, 1) = -1;
      dshape.Elem(1, 2) = 1;
      dshape.Elem(1, 3) = 0;
      dshape.Elem(2, 1) = -1;
      dshape.Elem(2, 2) = 0;
      dshape.Elem(2, 3) = 1;
      break;
    case QUAD:
      dshape.Elem(1, 1) = -(1-p.Y());
      dshape.Elem(1, 2) = (1-p.Y());
      dshape.Elem(1, 3) = p.Y();
      dshape.Elem(1, 4) = -p.Y();
      dshape.Elem(2, 1) = -(1-p.X());
      dshape.Elem(2, 2) = -p.X();
      dshape.Elem(2, 3) = p.X();
      dshape.Elem(2, 4) = (1-p.X());
      break;

    default:
      PrintSysError ("Element2d::GetDShape, illegal type ", typ);
    }
}




void Element2d :: 
GetDShapeNew (const Point<2> & p, MatrixFixWidth<2> & dshape) const
{
  switch (typ)
    {
    case TRIG:
      {
	dshape = 0;
	dshape(0,0) = 1;
	dshape(1,1) = 1;
	dshape(2,0) = -1;
	dshape(2,1) = -1;
	break;
      }
    case QUAD:
      {
	dshape(0,0) = -(1-p(1));
	dshape(0,1) = -(1-p(0));

	dshape(1,0) =  (1-p(1));
	dshape(1,1) =  -p(0);

	dshape(2,0) = p(1);
	dshape(2,1) = p(0);

	dshape(3,0) = -p(1);
	dshape(3,1) = (1-p(0));
	break;
      }
    }
}





void Element2d :: 
GetPointMatrix (const ARRAY<Point2d> & points,
		DenseMatrix & pmat) const
{
  int np = GetNP();

#ifdef DEBUG
  if (pmat.Width() != np || pmat.Height() != 2)
    {
      cerr << "Element::GetPointMatrix: sizes don't fit" << endl;
      return;
    }
#endif
  
  for (int i = 1; i <= np; i++)
    {
      const Point2d & p = points.Get(PNum(i));
      pmat.Elem(1, i) = p.X();
      pmat.Elem(2, i) = p.Y();
    }
}





double Element2d :: CalcJacobianBadness (const ARRAY<Point2d> & points) const
{
  int i, j;
  int nip = GetNIP();
  static DenseMatrix trans(2,2);
  static DenseMatrix pmat;
  
  pmat.SetSize (2, GetNP());
  GetPointMatrix (points, pmat);

  double err = 0;
  for (i = 1; i <= nip; i++)
    {
      GetTransformation (i, pmat, trans);

      // Frobenius norm
      double frob = 0;
      for (j = 1; j <= 4; j++)
	frob += sqr (trans.Get(j));
      frob = sqrt (frob);
      frob /= 2;

      double det = trans.Det();

      if (det <= 0)
	err += 1e12;
      else
	err += frob * frob / det;
    }

  err /= nip;
  return err;
}



static const int qip_table[4][4] =
  { { 0, 1, 0, 3 },
    { 0, 1, 1, 2 },
    { 3, 2, 0, 3 },
    { 3, 2, 1, 2 }
  };

double Element2d :: 
CalcJacobianBadnessDirDeriv (const ARRAY<Point2d> & points,
			     int pi, Vec2d & dir, double & dd) const
{
  if (typ == QUAD)
    {
      Mat<2,2> trans, dtrans;
      Mat<2,4> vmat, pmat;
      
      for (int j = 0; j < 4; j++)
	{
	  const Point2d & p = points.Get( (*this)[j] );
	  pmat(0, j) = p.X();
	  pmat(1, j) = p.Y();
	}

      vmat = 0.0;
      vmat(0, pi-1) = dir.X();
      vmat(1, pi-1) = dir.Y();
      
      double err = 0;
      dd = 0;

      for (int i = 0; i < 4; i++)
	{
	  int ix1 = qip_table[i][0];
	  int ix2 = qip_table[i][1];
	  int iy1 = qip_table[i][2];
	  int iy2 = qip_table[i][3];
	      
	  trans(0,0) = pmat(0, ix2) - pmat(0,ix1);
	  trans(1,0) = pmat(1, ix2) - pmat(1,ix1);
	  trans(0,1) = pmat(0, iy2) - pmat(0,iy1);
	  trans(1,1) = pmat(1, iy2) - pmat(1,iy1);

	  double det = trans(0,0)*trans(1,1)-trans(1,0)*trans(0,1);

	  if (det <= 0)
	    {
	      dd = 0;
	      return 1e12;
	    }
	  
	  dtrans(0,0) = vmat(0, ix2) - vmat(0,ix1);
	  dtrans(1,0) = vmat(1, ix2) - vmat(1,ix1);
	  dtrans(0,1) = vmat(0, iy2) - vmat(0,iy1);
	  dtrans(1,1) = vmat(1, iy2) - vmat(1,iy1);


	  // Frobenius norm
	  double frob = 0;
	  for (int j = 0; j < 4; j++) 
	    frob += sqr (trans(j));
	  frob = sqrt (frob);
	  
	  double dfrob = 0;
	  for (int j = 0; j < 4; j++)
	    dfrob += trans(j) * dtrans(j);
	  dfrob = dfrob / frob;
	  
	  frob /= 2;      
	  dfrob /= 2;
	  
	  
	  // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
	  double ddet 
	    = dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0)
	    + trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);
	  
	  err += frob * frob / det;
	  dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
	}
      
      err /= 4;
      dd /= 4;
      return err;
    }

  int nip = GetNIP();
  static DenseMatrix trans(2,2), dtrans(2,2);
  static DenseMatrix pmat, vmat;
  
  pmat.SetSize (2, GetNP());
  vmat.SetSize (2, GetNP());

  GetPointMatrix (points, pmat);
  
  vmat = 0.0;
  vmat.Elem(1, pi) = dir.X();
  vmat.Elem(2, pi) = dir.Y();


  double err = 0;
  dd = 0;

  for (int i = 1; i <= nip; i++)
    {
      GetTransformation (i, pmat, trans);
      GetTransformation (i, vmat, dtrans);

      // Frobenius norm
      double frob = 0;
      for (int j = 1; j <= 4; j++)
	frob += sqr (trans.Get(j));
      frob = sqrt (frob);
      
      double dfrob = 0;
      for (int j = 1; j <= 4; j++)
	dfrob += trans.Get(j) * dtrans.Get(j);
      dfrob = dfrob / frob;
      
      frob /= 2;      
      dfrob /= 2;
      
      double det = trans(0,0)*trans(1,1)-trans(1,0)*trans(0,1);

      // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
      double ddet 
	= dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0)
	+ trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);

      if (det <= 0)
	err += 1e12;
      else
	{
	  err += frob * frob / det;
	  dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
	}
    }

  err /= nip;
  dd /= nip;
  return err;
}



double Element2d :: 
CalcJacobianBadness (const T_POINTS & points, const Vec<3> & n) const
{
  int i, j;
  int nip = GetNIP();
  static DenseMatrix trans(2,2);
  static DenseMatrix pmat;
  
  pmat.SetSize (2, GetNP());

  Vec<3> t1, t2;
  t1 = n.GetNormal();
  t2 = Cross (n, t1);

  for (i = 1; i <= GetNP(); i++)
    {
      Point3d p = points.Get(PNum(i));
      pmat.Elem(1, i) = p.X() * t1(0) + p.Y() * t1(1) + p.Z() * t1(2);
      pmat.Elem(2, i) = p.X() * t2(0) + p.Y() * t2(1) + p.Z() * t2(2);
    }

  double err = 0;
  for (i = 1; i <= nip; i++)
    {
      GetTransformation (i, pmat, trans);

      // Frobenius norm
      double frob = 0;
      for (j = 1; j <= 4; j++)
	frob += sqr (trans.Get(j));
      frob = sqrt (frob);
      frob /= 2;

      double det = trans.Det();
      if (det <= 0)
	err += 1e12;
      else
	err += frob * frob / det;
    }

  err /= nip;
  return err;
}



void Element2d :: ComputeIntegrationPointData () const
{
  switch (np)
    {
    case 3: if (ipdtrig.Size()) return; break;
    case 4: if (ipdquad.Size()) return; break;
    }

  for (int i = 1; i <= GetNIP(); i++)
    {
      IntegrationPointData * ipd = new IntegrationPointData;
      Point2d hp;
      GetIntegrationPoint (i, hp, ipd->weight);
      ipd->p(0) = hp.X();
      ipd->p(1) = hp.Y();
      ipd->p(2) = 0;

      ipd->shape.SetSize(GetNP());
      ipd->dshape.SetSize(2, GetNP());

      GetShape (hp, ipd->shape);
      GetDShape (hp, ipd->dshape);

      switch (np)
	{
	case 3: ipdtrig.Append (ipd); break;
	case 4: ipdquad.Append (ipd); break;
	}
    }
}










ostream & operator<<(ostream  & s, const Element2d & el)
{
  s << "np = " << el.GetNP();
  for (int j = 1; j <= el.GetNP(); j++)
    s << " " << el.PNum(j);
  return s;
}


ostream & operator<<(ostream  & s, const Element & el)
{
  s << "np = " << el.GetNP();
  for (int j = 0; j < el.GetNP(); j++)
    s << " " << int(el[j]);
  return s;
}


Element :: Element ()
{
  typ = TET;
  np = 4;
  for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
    pnum[i] = 0;
  index = 0;
  flags.marked = 1;
  flags.badel = 0;
  flags.reverse = 0;
  flags.illegal = 0;
  flags.illegal_valid = 0;
  flags.badness_valid = 0;
  flags.refflag = 1;
  flags.strongrefflag = false;
  flags.deleted = 0;
  flags.fixed = 0;
  orderx = ordery = orderz = 1;

#ifdef PARALLEL
  partitionNumber = -1;
  isghost = 0;
#endif

}


Element :: Element (int anp)
{
  np = anp;
  int i;
  for (i = 0; i < ELEMENT_MAXPOINTS; i++)
    pnum[i] = 0;
  index = 0;
  flags.marked = 1;
  flags.badel = 0;
  flags.reverse = 0;
  flags.illegal = 0;
  flags.illegal_valid = 0;
  flags.badness_valid = 0;
  flags.refflag = 1;
  flags.strongrefflag = false;
  flags.deleted = 0;
  flags.fixed = 0;

  switch (np)
    {
    case 4: typ = TET; break;
    case 5: typ = PYRAMID; break;
    case 6: typ = PRISM; break;
    case 8: typ = HEX; break;
    case 10: typ = TET10; break;
    default: cerr << "Element::Element: unknown element with " << np << " points" << endl;
    }
  orderx = ordery = orderz = 1;

#ifdef PARALLEL
  isghost = 0;
#endif
}

void Element :: SetOrder (const int aorder) 
  { 
    orderx = aorder; 
    ordery = aorder; 
    orderz = aorder;
  }


void Element :: SetOrder (const int ox, const int oy, const int oz) 
{ 
  orderx = ox; 
  ordery = oy;
  orderz = oz; 
}


Element :: Element (ELEMENT_TYPE type)
{
  SetType (type);

  int i;
  for (i = 0; i < ELEMENT_MAXPOINTS; i++)
    pnum[i] = 0;
  index = 0;
  flags.marked = 1;
  flags.badel = 0;
  flags.reverse = 0;
  flags.illegal = 0;
  flags.illegal_valid = 0;
  flags.badness_valid = 0;
  flags.refflag = 1;
  flags.strongrefflag = false;
  flags.deleted = 0;
  flags.fixed = 0;
  orderx = ordery = orderz = 1;

#ifdef PARALLEL
  isghost = 0;
#endif
}





Element & Element :: operator= (const Element & el2)
{
  typ = el2.typ;
  np = el2.np;
  for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
    pnum[i] = el2.pnum[i];
  index = el2.index;
  flags = el2.flags;
  orderx = el2.orderx;
  ordery = el2.ordery;
  orderz = el2.orderz;
  hp_elnr = el2.hp_elnr;
  flags = el2.flags;
  return *this;
}



void Element :: SetNP (int anp)
{
  np = anp; 
  switch (np)
    {
    case 4: typ = TET; break;
    case 5: typ = PYRAMID; break;
    case 6: typ = PRISM; break;
    case 8: typ = HEX; break;
    case 10: typ = TET10; break;
      // 
    default: break;
      cerr << "Element::SetNP unknown element with " << np << " points" << endl;
    }
}



void Element :: SetType (ELEMENT_TYPE atyp)
{
  typ = atyp;
  switch (atyp)
    {
    case TET: np = 4; break;
    case PYRAMID: np = 5; break;
    case PRISM: np = 6; break;
    case HEX: np = 8; break;
    case TET10: np = 10; break;
    case PRISM12: np = 12; break;
    }
}



void Element :: Invert()
{
  switch (GetNP())
    {
    case 4:
      {
	Swap (PNum(3), PNum(4));
	break;
      }
    case 5:
      {
	Swap (PNum(1), PNum(4));
	Swap (PNum(2), PNum(3));
	break;
      }
    case 6:
      {
	Swap (PNum(1), PNum(4));
	Swap (PNum(2), PNum(5));
	Swap (PNum(3), PNum(6));
	break;
      }
    }
}


void Element :: Print (ostream & ost) const
{
  ost << np << " Points: ";
  for (int i = 1; i <= np; i++)
    ost << pnum[i-1] << " " << endl;
}

void Element :: GetBox (const T_POINTS & points, Box3d & box) const
{
  box.SetPoint (points.Get(PNum(1)));
  box.AddPoint (points.Get(PNum(2)));
  box.AddPoint (points.Get(PNum(3)));
  box.AddPoint (points.Get(PNum(4)));
}

double Element :: Volume (const T_POINTS & points) const
{
  Vec<3> v1 = points.Get(PNum(2)) - points.Get(PNum(1));
  Vec<3> v2 = points.Get(PNum(3)) - points.Get(PNum(1));
  Vec<3> v3 = points.Get(PNum(4)) - points.Get(PNum(1)); 
  
  return -(Cross (v1, v2) * v3) / 6;	 
}  


void Element :: GetFace2 (int i, Element2d & face) const
{
  static const int tetfaces[][5] = 
  { { 3, 2, 3, 4, 0 },
    { 3, 3, 1, 4, 0 },
    { 3, 1, 2, 4, 0 },
    { 3, 2, 1, 3, 0 } };

  static const int pyramidfaces[][5] =
  { { 4, 1, 4, 3, 2 },
    { 3, 1, 2, 5, 0 },
    { 3, 2, 3, 5, 0 },
    { 3, 3, 4, 5, 0 },
    { 3, 4, 1, 5, 0 } };

  static const int prismfaces[][5] =
  {
    { 3, 1, 3, 2, 0 },
    { 3, 4, 5, 6, 0 },
    { 4, 1, 2, 5, 4 },
    { 4, 2, 3, 6, 5 },
    { 4, 3, 1, 4, 6 }
  };

  switch (np)
    {
    case 4: // tet
    case 10: // tet
      {
	face.SetType(TRIG);
	for (int j = 1; j <= 3; j++)
	  face.PNum(j) = PNum(tetfaces[i-1][j]);
	break;
      }
    case 5: // pyramid
      {
	// face.SetNP(pyramidfaces[i-1][0]);
	face.SetType ( (i == 1) ? QUAD : TRIG);
	for (int j = 1; j <= face.GetNP(); j++)
	  face.PNum(j) = PNum(pyramidfaces[i-1][j]);
	break;
      }
    case 6: // prism
      {
	//	face.SetNP(prismfaces[i-1][0]);
	face.SetType ( (i >= 3) ? QUAD : TRIG);
	for (int j = 1; j <= face.GetNP(); j++)
	  face.PNum(j) = PNum(prismfaces[i-1][j]);
	break;
      }
    }
}



void Element :: GetTets (ARRAY<Element> & locels) const
{
  GetTetsLocal (locels);
  int i, j;
  for (i = 1; i <= locels.Size(); i++)
    for (j = 1; j <= 4; j++)
      locels.Elem(i).PNum(j) = PNum ( locels.Elem(i).PNum(j) );
}

void Element :: GetTetsLocal (ARRAY<Element> & locels) const
{
  int i, j;
  locels.SetSize(0);
  switch (GetType())
    {
    case TET:
      {
	int linels[1][4] = 
	{ { 1, 2, 3, 4 },
	};
	for (i = 0; i < 1; i++)
	  {
	    Element tet(4);
	    for (j = 1; j <= 4; j++)
	      tet.PNum(j) = linels[i][j-1];
	    locels.Append (tet);
	  }
	break;
      }
    case TET10:
      {
	int linels[8][4] = 
	{ { 1, 5, 6, 7 },
	  { 5, 2, 8, 9 },
	  { 6, 8, 3, 10 },
	  { 7, 9, 10, 4 },
	  { 5, 6, 7, 9 },
	  { 5, 6, 9, 8 },
	  { 6, 7, 9, 10 },
	  { 6, 8, 10, 9 } };
	for (i = 0; i < 8; i++)
	  {
	    Element tet(4);
	    for (j = 1; j <= 4; j++)
	      tet.PNum(j) = linels[i][j-1];
	    locels.Append (tet);
	  }
	break;
      }
    case PYRAMID:
      {
	int linels[2][4] = 
	{ { 1, 2, 3, 5 },
	  { 1, 3, 4, 5 } };
	for (i = 0; i < 2; i++)
	  {
	    Element tet(4);
	    for (j = 1; j <= 4; j++)
	      tet.PNum(j) = linels[i][j-1];
	    locels.Append (tet);
	  }
	break;
      }
    case PRISM:
    case PRISM12:
      {
	int linels[3][4] = 
	{ { 1, 2, 3, 4 },
	  { 4, 2, 3, 5 },
	  { 6, 5, 4, 3 }
	};
	for (i = 0; i < 3; i++)
	  {
	    Element tet(4);
	    for (j = 0; j < 4; j++)
	      tet[j] = linels[i][j];
	    locels.Append (tet);
	  }
	break;
      }
    case HEX:
      {
	int linels[6][4] = 
	{ { 1, 7, 2, 3 },
	  { 1, 7, 3, 4 },
	  { 1, 7, 4, 8 },
	  { 1, 7, 8, 5 },
	  { 1, 7, 5, 6 },
	  { 1, 7, 6, 2 }
	};
	for (i = 0; i < 6; i++)
	  {
	    Element tet(4);
	    for (j = 0; j < 4; j++)
	      tet[j] = linels[i][j];
	    locels.Append (tet);
	  }
	break;
      }
    default:
      {
	cerr << "GetTetsLocal not implemented for el with " << GetNP() << " nodes" << endl;
      }
    }
}

bool Element :: operator==(const Element & el2) const
{
  bool retval = (el2.GetNP() == np);
  for(int i= 0; retval && i<np; i++)
    retval = (el2[i] == (*this)[i]);

  return retval;
}


#ifdef OLD
void Element :: GetNodesLocal (ARRAY<Point3d> & points) const
{
  const static double tetpoints[4][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 }};
  
  const static double prismpoints[6][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 1, 0, 1 },
      { 0, 1, 1 } };
  
  const static double pyramidpoints[6][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 } };
  
  const static double tet10points[10][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 0.5, 0, 0 },
      { 0, 0.5, 0 },
      { 0, 0, 0.5 },
      { 0.5, 0.5, 0 },
      { 0.5, 0, 0.5 },
      { 0, 0.5, 0.5 } };

  const static double hexpoints[8][3] =
    { 
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 1, 0, 1 },
      { 1, 1, 1 },
      { 0, 1, 1 }
    };
  
  int np, i;
  const double (*pp)[3];
  switch (GetType())
    {
    case TET:
      {
	np = 4;
	pp = tetpoints;
	break;
      }
    case PRISM:
    case PRISM12:
      {
	np = 6;
	pp = prismpoints;
	break;
      }
    case TET10:
      {
	np = 10;
	pp = tet10points;
	break;
      }
    case PYRAMID:
      {
	np = 5;
	pp = pyramidpoints;
	break;
      }
    case HEX:
      {
	np = 8;
	pp = hexpoints;
	break;
      }
    default:
      {
	cout << "GetNodesLocal not impelemented for element " << GetType() << endl;
	np = 0;
      }
    }
  
  points.SetSize(0);
  for (i = 0; i < np; i++)
    points.Append (Point3d (pp[i][0], pp[i][1], pp[i][2]));
}
#endif






void Element :: GetNodesLocalNew (ARRAY<Point<3> > & points) const
{
  const static double tetpoints[4][3] =
    {      
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 0, 0, 0 }
    };
  
  const static double prismpoints[6][3] =
    {
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 0, 0, 1 }
    };
  
  const static double pyramidpoints[6][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 } };
  
  const static double tet10points[10][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 0.5, 0, 0 },
      { 0, 0.5, 0 },
      { 0, 0, 0.5 },
      { 0.5, 0.5, 0 },
      { 0.5, 0, 0.5 },
      { 0, 0.5, 0.5 } };

  const static double hexpoints[8][3] =
    { 
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 1, 0, 1 },
      { 1, 1, 1 },
      { 0, 1, 1 }
    };
  

  
  int np, i;
  const double (*pp)[3];
  switch (GetType())
    {
    case TET:
      {
	np = 4;
	pp = tetpoints;
	break;
      }
    case PRISM:
    case PRISM12:
      {
	np = 6;
	pp = prismpoints;
	break;
      }
    case TET10:
      {
	np = 10;
	pp = tet10points;
	break;
      }
    case PYRAMID:
      {
	np = 5;
	pp = pyramidpoints;
	break;
      }
    case HEX:
      {
	np = 8;
	pp = hexpoints;
	break;
      }
    default:
      {
	cout << "GetNodesLocal not impelemented for element " << GetType() << endl;
	np = 0;
      }
    }
  
  points.SetSize(0);
  for (i = 0; i < np; i++)
    points.Append (Point<3> (pp[i][0], pp[i][1], pp[i][2]));
}

















void Element :: GetSurfaceTriangles (ARRAY<Element2d> & surftrigs) const
{
  static int tet4trigs[][3] = 
  { { 2, 3, 4 },
    { 3, 1, 4 },
    { 1, 2, 4 },
    { 2, 1, 3 } };

  static int tet10trigs[][3] = 
  { { 2, 8, 9 }, { 3, 10, 8}, { 4, 9, 10 }, { 9, 8, 10 },
    { 3, 6, 10 }, { 1, 7, 6 }, { 4, 10, 7 }, { 6, 7, 10 },
    { 1, 5, 7 }, { 2, 9, 5 }, { 4, 7, 9 }, { 5, 9, 7 },
    { 1, 6, 5 }, { 2, 5, 8 }, { 3, 8, 6 }, { 5, 6, 8 }
  };

  static int pyramidtrigs[][3] =
  {
    { 1, 3, 2 },
    { 1, 4, 3 },
    { 1, 2, 5 },
    { 2, 3, 5 },
    { 3, 4, 5 },
    { 4, 1, 5 }
  };

  static int prismtrigs[][3] =
    {
      { 1, 3, 2 },
      { 4, 5, 6 },
      { 1, 2, 4 },
      { 4, 2, 5 },
      { 2, 3, 5 },
      { 5, 3, 6 },
      { 3, 1, 6 },
      { 6, 1, 4 }
    };
  
  static int hextrigs[][3] = 
    {
      { 1, 3, 2 },
      { 1, 4, 3 }, 
      { 5, 6, 7 },
      { 5, 7, 8 },
      { 1, 2, 6 },
      { 1, 6, 5 },
      { 2, 3, 7 },
      { 2, 7, 6 },
      { 3, 4, 8 },
      { 3, 8, 7 },
      { 4, 1, 8 },
      { 1, 5, 8 }
    };

  int j;

  int nf;
  int (*fp)[3];

  switch (GetType())
    {
    case TET:
      {
	nf = 4;
	fp = tet4trigs;
	break;
      }
    case PYRAMID:
      {
	nf = 6;
	fp = pyramidtrigs;
	break;
      }
    case PRISM:
    case PRISM12:
      {
	nf = 8;
	fp = prismtrigs;
	break;
      }
    case TET10:
      {
	nf = 16;
	fp = tet10trigs;
	break;
      }
    case HEX:
      {
	nf = 12;
	fp = hextrigs;
	break;
      }
    default:
      {
	nf = 0;
	fp = NULL;
      }
    }

  
  surftrigs.SetSize (nf);
  for (j = 0; j < nf; j++)
    {
      surftrigs.Elem(j+1) = Element2d(TRIG);
      surftrigs.Elem(j+1).PNum(1) = fp[j][0];
      surftrigs.Elem(j+1).PNum(2) = fp[j][1];
      surftrigs.Elem(j+1).PNum(3) = fp[j][2];
    }
}





ARRAY< AutoPtr < IntegrationPointData > > ipdtet;
ARRAY< AutoPtr < IntegrationPointData > > ipdtet10;



int Element :: GetNIP () const
{
  int nip;
  switch (typ)
    {
    case TET: nip = 1; break;
    case TET10: nip = 8; break;
    default: nip = 0; break;
    }
  return nip;
}

void Element :: 
GetIntegrationPoint (int ip, Point<3> & p, double & weight) const
{
  static double eltetqp[1][4] =
  {
    { 0.25, 0.25, 0.25, 1.0/6.0 }
  };

  static double eltet10qp[8][4] =
  {
    { 0.585410196624969, 0.138196601125011, 0.138196601125011, 1.0/24.0 },
    { 0.138196601125011, 0.585410196624969, 0.138196601125011, 1.0/24.0 },
    { 0.138196601125011, 0.138196601125011, 0.585410196624969, 1.0/24.0 },
    { 0.138196601125011, 0.138196601125011, 0.138196601125011, 1.0/24.0 },
    { 1, 0, 0, 1 },
    { 0, 1, 0, 1 },
    { 0, 0, 1, 1 },
    { 0, 0, 0, 1 },
  };
  
  double * pp;
  switch (typ)
    {
    case TET: pp = &eltetqp[0][0]; break;
    case TET10: pp = &eltet10qp[ip-1][0]; break;
    }

  p(0) = pp[0];
  p(1) = pp[1];
  p(2) = pp[2];
  weight = pp[3];
}

void Element :: 
GetTransformation (int ip, const T_POINTS & points,
		   DenseMatrix & trans) const
{
  int np = GetNP();
  static DenseMatrix pmat(3, np), dshape(3, np);
  pmat.SetSize (3, np);
  dshape.SetSize (3, np);

  Point<3> p;
  double w;

  GetPointMatrix (points, pmat);
  GetIntegrationPoint (ip, p, w);
  GetDShape (p, dshape);
  
  CalcABt (pmat, dshape, trans);

  /*
  (*testout) << "p = " << p  << endl
	     << "pmat = " << pmat << endl
	     << "dshape = " << dshape << endl
	     << "tans = " << trans << endl;
  */
}

void Element :: 
GetTransformation (int ip, class DenseMatrix & pmat,
		   class DenseMatrix & trans) const
{
  int np = GetNP();

  if (pmat.Width() != np || pmat.Height() != 3)
    {
      (*testout) << "GetTransofrmation: pmat doesn't fit" << endl;
      return;
    }

  ComputeIntegrationPointData ();
  DenseMatrix * dshapep;
  switch (GetType())
    {
    case TET: dshapep = &ipdtet.Get(ip)->dshape; break;
    case TET10: dshapep = &ipdtet10.Get(ip)->dshape; break;
    }
  
  CalcABt (pmat, *dshapep, trans);
}



void Element :: GetShape (const Point<3> & hp, Vector & shape) const
{
  Point3d p = hp;

  if (shape.Size() != GetNP())
    {
      cerr << "Element::GetShape: Length not fitting" << endl;
      return;
    }

  switch (typ)
    {
    case TET:
      {
	shape.Elem(1) = 1 - p.X() - p.Y() - p.Z(); 
	shape.Elem(2) = p.X();
	shape.Elem(3) = p.Y();
	shape.Elem(4) = p.Z();
	break;
      }
    case TET10:
      {
	double lam1 = 1 - p.X() - p.Y() - p.Z();
	double lam2 = p.X();
	double lam3 = p.Y();
	double lam4 = p.Z();
	
	shape.Elem(5) = 4 * lam1 * lam2;
	shape.Elem(6) = 4 * lam1 * lam3;
	shape.Elem(7) = 4 * lam1 * lam4;
	shape.Elem(8) = 4 * lam2 * lam3;
	shape.Elem(9) = 4 * lam2 * lam4;
	shape.Elem(10) = 4 * lam3 * lam4;
	
	shape.Elem(1) = lam1 - 
	  0.5 * (shape.Elem(5) + shape.Elem(6) + shape.Elem(7));
	shape.Elem(2) = lam2 - 
	  0.5 * (shape.Elem(5) + shape.Elem(8) + shape.Elem(9));
	shape.Elem(3) = lam3 - 
	  0.5 * (shape.Elem(6) + shape.Elem(8) + shape.Elem(10));
	shape.Elem(4) = lam4 - 
	  0.5 * (shape.Elem(7) + shape.Elem(9) + shape.Elem(10));
	break;
      }

    case PRISM:
      {
	Point<3> hp = p; 
	shape(0) = hp(0) * (1-hp(2));
	shape(1) = hp(1) * (1-hp(2));
	shape(2) = (1-hp(0)-hp(1)) * (1-hp(2));
	shape(3) = hp(0) * hp(2);
	shape(4) = hp(1) * hp(2);
	shape(5) = (1-hp(0)-hp(1)) * hp(2);
	break;
      }
    case HEX:
      {
	Point<3> hp = p; 
	shape(0) = (1-hp(0))*(1-hp(1))*(1-hp(2));
	shape(1) = (  hp(0))*(1-hp(1))*(1-hp(2));
	shape(2) = (  hp(0))*(  hp(1))*(1-hp(2));
	shape(3) = (1-hp(0))*(  hp(1))*(1-hp(2));
	shape(4) = (1-hp(0))*(1-hp(1))*(  hp(2));
	shape(5) = (  hp(0))*(1-hp(1))*(  hp(2));
	shape(6) = (  hp(0))*(  hp(1))*(  hp(2));
	shape(7) = (1-hp(0))*(  hp(1))*(  hp(2));
	break;
      }
    }
}



void Element :: GetShapeNew (const Point<3> & p, FlatVector & shape) const
{
  /*
  if (shape.Size() < GetNP())
    {
      cerr << "Element::GetShape: Length not fitting" << endl;
      return;
    }
  */

  switch (typ)
    {
    case TET:
      {
	shape(0) = p(0);
	shape(1) = p(1);
	shape(2) = p(2);
	shape(3) = 1-p(0)-p(1)-p(2);
	break;
      }

    case PYRAMID:
      {
	double noz = 1-p(2);
	if (noz == 0.0) noz = 1e-10;

	double xi  = p(0) / noz;
	double eta = p(1) / noz;
	shape(0) = (1-xi)*(1-eta) * (noz);
	shape(1) = (  xi)*(1-eta) * (noz);
	shape(2) = (  xi)*(  eta) * (noz);
	shape(3) = (1-xi)*(  eta) * (noz);
	shape(4) = p(2);
	break;
      }

    case PRISM:
      {
	shape(0) = p(0) * (1-p(2));
	shape(1) = p(1) * (1-p(2));
	shape(2) = (1-p(0)-p(1)) * (1-p(2));
	shape(3) = p(0) * p(2);
	shape(4) = p(1) * p(2);
	shape(5) = (1-p(0)-p(1)) * p(2);
	break;
      }
    case HEX:
      {
	shape(0) = (1-p(0))*(1-p(1))*(1-p(2));
	shape(1) = (  p(0))*(1-p(1))*(1-p(2));
	shape(2) = (  p(0))*(  p(1))*(1-p(2));
	shape(3) = (1-p(0))*(  p(1))*(1-p(2));
	shape(4) = (1-p(0))*(1-p(1))*(  p(2));
	shape(5) = (  p(0))*(1-p(1))*(  p(2));
	shape(6) = (  p(0))*(  p(1))*(  p(2));
	shape(7) = (1-p(0))*(  p(1))*(  p(2));
	break;
      }
    }
}




void Element :: 
GetDShape (const Point<3> & hp, DenseMatrix & dshape) const
{
  Point3d p = hp;

  int np = GetNP();
  if (dshape.Height() != 3 || dshape.Width() != np)
    {
      cerr << "Element::DShape: Sizes don't fit" << endl;
      return;
    }

  int i, j;
  double eps = 1e-6;
  Vector shaper(np), shapel(np);

  for (i = 1; i <= 3; i++)
    {
      Point3d pr(p), pl(p);
      pr.X(i) += eps;
      pl.X(i) -= eps;
      
      GetShape (pr, shaper);
      GetShape (pl, shapel);
      for (j = 1; j <= np; j++)
	dshape.Elem(i, j) = (shaper.Get(j) - shapel.Get(j)) / (2 * eps);
    }
}



void Element :: 
GetDShapeNew (const Point<3> & p, MatrixFixWidth<3> & dshape) const
{
  switch (typ)
    {
    case TET:
      {
	dshape = 0;
	dshape(0,0) = 1;
	dshape(1,1) = 1;
	dshape(2,2) = 1;
	dshape(3,0) = -1;
	dshape(3,1) = -1;
	dshape(3,2) = -1;
	break;
      }
    case PRISM:
      {
	dshape = 0;
	dshape(0,0) = 1-p(2);
	dshape(0,2) = -p(0);
	dshape(1,1) = 1-p(2);
	dshape(1,2) = -p(1);
	dshape(2,0) = -(1-p(2));
	dshape(2,1) = -(1-p(2));
	dshape(2,2) = -(1-p(0)-p(1));

	dshape(3,0) = p(2);
	dshape(3,2) = p(0);
	dshape(4,1) = p(2);
	dshape(4,2) = p(1);
	dshape(5,0) = -p(2);
	dshape(5,1) = -p(2);
	dshape(5,2) = 1-p(0)-p(1);
	break;
      }

    default:
      {
	int np = GetNP();
	double eps = 1e-6;
	Vector shaper(np), shapel(np);
	
	for (int i = 1; i <= 3; i++)
	  {
	    Point3d pr(p), pl(p);
	    pr.X(i) += eps;
	    pl.X(i) -= eps;
	    
	    GetShapeNew (pr, shaper);
	    GetShapeNew (pl, shapel);
	    for (int j = 1; j <= np; j++)
	      dshape.Elem(j, i) = (shaper.Get(j) - shapel.Get(j)) / (2 * eps);
	  }
      }
    }
}

void Element :: 
GetPointMatrix (const T_POINTS & points,
		DenseMatrix & pmat) const
{
  int np = GetNP();
  /*
  if (pmat.Width() != np || pmat.Height() != 3)
    {
      cerr << "Element::GetPointMatrix: sizes don't fit" << endl;
      return;
    }
  */
  for (int i = 1; i <= np; i++)
    {
      const Point3d & p = points.Get(PNum(i));
      pmat.Elem(1, i) = p.X();
      pmat.Elem(2, i) = p.Y();
      pmat.Elem(3, i) = p.Z();
    }
}






double Element :: CalcJacobianBadness (const T_POINTS & points) const
{
  int i, j;
  int nip = GetNIP();
  static DenseMatrix trans(3,3);
  static DenseMatrix pmat;
  
  pmat.SetSize (3, GetNP());
  GetPointMatrix (points, pmat);

  double err = 0;
  for (i = 1; i <= nip; i++)
    {
      GetTransformation (i, pmat, trans);

      // Frobenius norm
      double frob = 0;
      for (j = 1; j <= 9; j++)
	frob += sqr (trans.Get(j));
      frob = sqrt (frob);
      frob /= 3;

      double det = -trans.Det();
      
      if (det <= 0)
	err += 1e12;
      else
	err += frob * frob * frob / det;
    }

  err /= nip;
  return err;
}

double Element :: 
CalcJacobianBadnessDirDeriv (const T_POINTS & points,
			     int pi, Vec<3> & dir, double & dd) const
{
  int i, j, k, l;
  int nip = GetNIP();
  static DenseMatrix trans(3,3), dtrans(3,3), hmat(3,3);
  static DenseMatrix pmat, vmat;
  
  pmat.SetSize (3, GetNP());
  vmat.SetSize (3, GetNP());

  GetPointMatrix (points, pmat);
  
  for (i = 1; i <= np; i++)
    for (j = 1; j <= 3; j++)
      vmat.Elem(j, i) = 0;
  for (j = 1; j <= 3; j++)
    vmat.Elem(j, pi) = dir(j-1);



  double err = 0;
  dd = 0;

  for (i = 1; i <= nip; i++)
    {
      GetTransformation (i, pmat, trans);
      GetTransformation (i, vmat, dtrans);


      // Frobenius norm
      double frob = 0;
      for (j = 1; j <= 9; j++)
	frob += sqr (trans.Get(j));
      frob = sqrt (frob);
      
      double dfrob = 0;
      for (j = 1; j <= 9; j++)
	dfrob += trans.Get(j) * dtrans.Get(j);
      dfrob = dfrob / frob;
      
      frob /= 3;      
      dfrob /= 3;

      
      double det = trans.Det();
      double ddet = 0;
      
      for (j = 1; j <= 3; j++)
	{
	  hmat = trans;
	  for (k = 1; k <= 3; k++)
	    hmat.Elem(k, j) = dtrans.Get(k, j);
	  ddet += hmat.Det();
	}


      det *= -1;
      ddet *= -1;

      
      if (det <= 0)
	err += 1e12;
      else
	{
	  err += frob * frob * frob / det;
	  dd += (3 * frob * frob * dfrob * det - frob * frob * frob * ddet) / (det * det);
	}
    }

  err /= nip;
  dd /= nip;
  return err;
}

double Element :: 
CalcJacobianBadnessGradient (const T_POINTS & points,
			     int pi, Vec<3> & grad) const
{
  int i, j, k, l;
  int nip = GetNIP();
  static DenseMatrix trans(3,3), dtrans(3,3), hmat(3,3);
  static DenseMatrix pmat, vmat;
  
  pmat.SetSize (3, GetNP());
  vmat.SetSize (3, GetNP());

  GetPointMatrix (points, pmat);
  
  for (i = 1; i <= np; i++)
    for (j = 1; j <= 3; j++)
      vmat.Elem(j, i) = 0;
  for (j = 1; j <= 3; j++)
    vmat.Elem(j, pi) = 1.;


  double err = 0;

  double dfrob[3];

  grad = 0;

  for (i = 1; i <= nip; i++)
    {
      GetTransformation (i, pmat, trans);
      GetTransformation (i, vmat, dtrans);
 
      // Frobenius norm
      double frob = 0;
      for (j = 1; j <= 9; j++)
	frob += sqr (trans.Get(j));
      frob = sqrt (frob);

      for(k = 0; k<3; k++)
	{
	  dfrob[k] = 0;
	  for (j = 1; j <= 3; j++)
	    dfrob[k] += trans.Get(k+1,j) * dtrans.Get(k+1,j);
	  dfrob[k] = dfrob[k] / (3.*frob);
	}

      frob /= 3;      

      double det = trans.Det();
      double ddet[3]; // = 0;
      
      for(k=1; k<=3; k++)
	{
	  int km1 = (k > 1) ? (k-1) : 3;
	  int kp1 = (k < 3) ? (k+1) : 1;
	  ddet[k-1] = 0;
	  for(j=1; j<=3; j++)
	    {
	      int jm1 = (j > 1) ? (j-1) : 3;
	      int jp1 = (j < 3) ? (j+1) : 1;
	      
	      ddet[k-1] += (-1.)* dtrans.Get(k,j) * ( trans.Get(km1,jm1)*trans.Get(kp1,jp1) - 
						     trans.Get(km1,jp1)*trans.Get(kp1,jm1) );
	    }
	}

      
      det *= -1;
      
      if (det <= 0)
	err += 1e12;
      else
	{
	  err += frob * frob * frob / det;
	  double fac = (frob * frob)/(det * det);
	  for(j=0; j<3; j++)
	    grad(j) += fac * (3 * dfrob[j] * det - frob * ddet[j]);
	}
    }

  err /= nip;
  grad *= 1./nip;
  return err;
}






void Element :: ComputeIntegrationPointData () const
{
  switch (GetType())
    {
    case TET: if (ipdtet.Size()) return; break;
    case TET10: if (ipdtet10.Size()) return; break;
    default:
      PrintSysError ("Element::ComputeIntegrationPoint, illegal type ", int(typ));
    }

  switch (GetType())
    {
    case TET: ipdtet.SetSize(GetNIP()); break;
    case TET10: ipdtet10.SetSize(GetNIP()); break;
    }


  for (int i = 1; i <= GetNIP(); i++)
    {
      IntegrationPointData * ipd = new IntegrationPointData;
      GetIntegrationPoint (i, ipd->p, ipd->weight);
      ipd->shape.SetSize(GetNP());
      ipd->dshape.SetSize(3, GetNP());

      GetShape (ipd->p, ipd->shape);
      GetDShape (ipd->p, ipd->dshape);

      switch (GetType())
	{
	case TET: ipdtet.Elem(i).Reset(ipd); break;
	case TET10: ipdtet10.Elem(i).Reset(ipd); break;
	default:
	  PrintSysError ("Element::ComputeIntegrationPoint(2), illegal type ", int(typ));
	}
    }
}








FaceDescriptor ::  FaceDescriptor()
{ 
  surfnr = domin = domout  = bcprop = 0; 
  domin_singular = domout_singular = 0.;
  tlosurf = -1; 
  bcname = 0;
  firstelement = -1;
}

FaceDescriptor ::  FaceDescriptor(const FaceDescriptor& other)
  : surfnr(other.surfnr), domin(other.domin), domout(other.domout),
    tlosurf(other.tlosurf), bcprop(other.bcprop), bcname(other.bcname),
    domin_singular(other.domin_singular), domout_singular(other.domout_singular)
{ 
  firstelement = -1;
}

FaceDescriptor :: 
FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi)
{ 
  surfnr = surfnri; 
  domin = domini; 
  domout = domouti;
  tlosurf = tlosurfi; 
  bcprop = surfnri;
  domin_singular = domout_singular = 0.;
  bcname = 0;
  firstelement = -1;
}

FaceDescriptor :: FaceDescriptor(const Segment & seg)
{ 
  surfnr = seg.si; 
  domin = seg.domin+1;
  domout = seg.domout+1;
  tlosurf = seg.tlosurf+1;
  bcprop = 0;
  domin_singular = domout_singular = 0.;
  bcname = 0;
  firstelement = -1;
}

int FaceDescriptor ::  SegmentFits (const Segment & seg)
{
  return
    surfnr == seg.si &&
    domin == seg.domin+1 &&
    domout == seg.domout+1  &&
    tlosurf == seg.tlosurf+1;
}


string FaceDescriptor :: GetBCName () const
{
  if ( bcname )
    return *bcname;
  else 
    return "default";
  
}

/*
void FaceDescriptor :: SetBCName (string * bcn)
{
  bcname = bcn;
}
*/


ostream & operator<<(ostream  & s, const FaceDescriptor & fd)
{
  s << "surfnr = " << fd.surfnr 
    << ", domin = " << fd.domin
    << ", domout = " << fd.domout
    << ", tlosurf = " << fd.tlosurf
    << ", bcprop = " << fd.bcprop
    << ", domin_sing = " << fd.domin_singular
    << ", domout_sing = " << fd.domout_singular;
  return s;
}






Identifications :: Identifications (Mesh & amesh)
  : mesh(amesh)
{
  identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
  identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
  maxidentnr = 0;
}

Identifications :: ~Identifications ()
{
  delete identifiedpoints;
  delete identifiedpoints_nr;
}

void Identifications :: Delete ()
{
  delete identifiedpoints;
  identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
  delete identifiedpoints_nr;
  identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
  maxidentnr = 0;
}

void Identifications :: Add (PointIndex pi1, PointIndex pi2, int identnr)
{
  //  (*testout) << "Identification::Add, pi1 = " << pi1 << ", pi2 = " << pi2 << ", identnr = " << identnr << endl;
  INDEX_2 pair (pi1, pi2);
  identifiedpoints->Set (pair, identnr);

  INDEX_3 tripl (pi1, pi2, identnr);
  identifiedpoints_nr->Set (tripl, 1);

  if (identnr > maxidentnr) maxidentnr = identnr;

  if (identnr+1 > idpoints_table.Size())
    idpoints_table.ChangeSize (identnr+1);
  idpoints_table.Add (identnr, pair);
  
  //  timestamp = NextTimeStamp();
}

int Identifications :: Get (PointIndex pi1, PointIndex pi2) const
{
  INDEX_2 pair(pi1, pi2);
  if (identifiedpoints->Used (pair))
    return identifiedpoints->Get(pair);
  else
    return 0;
}

bool Identifications :: Get (PointIndex pi1, PointIndex pi2, int nr) const
{
  INDEX_3 tripl(pi1, pi2, nr);
  if (identifiedpoints_nr->Used (tripl))
    return 1;
  else
    return 0;
}



int Identifications :: GetSymmetric (PointIndex pi1, PointIndex pi2) const
{
  INDEX_2 pair(pi1, pi2);
  if (identifiedpoints->Used (pair))
    return identifiedpoints->Get(pair);

  pair = INDEX_2 (pi2, pi1);
  if (identifiedpoints->Used (pair))
    return identifiedpoints->Get(pair);

  return 0;
}


void Identifications :: GetMap (int identnr, ARRAY<int,PointIndex::BASE> & identmap, bool symmetric) const
{
  identmap.SetSize (mesh.GetNP());
  identmap = 0;

  if (identnr)
    for (int i = 0; i < idpoints_table[identnr].Size(); i++)
      {
	INDEX_2 pair = idpoints_table[identnr][i];
	identmap[pair.I1()] = pair.I2();
	if(symmetric)
	  identmap[pair.I2()] = pair.I1();
      }

  else
    {
      cout << "getmap, identnr = " << identnr << endl;

      for (int i = 1; i <= identifiedpoints_nr->GetNBags(); i++)
	for (int j = 1; j <= identifiedpoints_nr->GetBagSize(i); j++)
	  {
	    INDEX_3 i3;
	    int dummy;
	    identifiedpoints_nr->GetData (i, j, i3, dummy);
	    
	    if (i3.I3() == identnr || !identnr)
	      {
		identmap.Elem(i3.I1()) = i3.I2();
		if(symmetric)
		  identmap.Elem(i3.I2()) = i3.I1();
	      }
	  }  
    }

}


void Identifications :: GetPairs (int identnr, 
				  ARRAY<INDEX_2> & identpairs) const
{
  identpairs.SetSize(0);
  
  if (identnr == 0)
    for (int i = 1; i <= identifiedpoints->GetNBags(); i++)
      for (int j = 1; j <= identifiedpoints->GetBagSize(i); j++)
	{
	  INDEX_2 i2;
	  int nr;
	  identifiedpoints->GetData (i, j, i2, nr);
	  identpairs.Append (i2);
	}  
  else
    for (int i = 1; i <= identifiedpoints_nr->GetNBags(); i++)
      for (int j = 1; j <= identifiedpoints_nr->GetBagSize(i); j++)
	{
	  INDEX_3 i3;
	  int dummy;
	  identifiedpoints_nr->GetData (i, j, i3 , dummy);
	  
	  if (i3.I3() == identnr)
	    identpairs.Append (INDEX_2(i3.I1(), i3.I2()));
	}  
}


void Identifications :: SetMaxPointNr (int maxpnum)
{
  for (int i = 1; i <= identifiedpoints->GetNBags(); i++)
    for (int j = 1; j <= identifiedpoints->GetBagSize(i); j++)
      {
	INDEX_2 i2;
	int nr;
	identifiedpoints->GetData (i, j, i2, nr);
	
	if (i2.I1() > maxpnum || i2.I2() > maxpnum)
	  {
	    i2.I1() = i2.I2() = -1;
	    identifiedpoints->SetData (i, j, i2, -1);	    
	  }
      }
}


void Identifications :: Print (ostream & ost) const
{
  ost << "Identifications:" << endl;
  ost << "pairs: " << endl << *identifiedpoints << endl;
  ost << "pairs and nr: " << endl << *identifiedpoints_nr << endl;
  ost << "table: " << endl << idpoints_table << endl;
}


MeshingParameters :: MeshingParameters ()
{
  optimize3d = "cmdmustm";
  //optimize3d = "cmdmstm";
  optsteps3d = 3;
  optimize2d = "smsmsmSmSmSm";
  optsteps2d = 3;
  opterrpow = 2;
  blockfill = 1;
  filldist = 0.1;
  safety = 5;
  relinnersafety = 3;
  uselocalh = 1;
  grading = 0.3;
  delaunay = 1;
  maxh = 1e10;
  minh = 0;
  meshsizefilename = NULL;
  startinsurface = 0;
  checkoverlap = 1;
  checkoverlappingboundary = 1;
  checkchartboundary = 1;
  curvaturesafety = 2;
  segmentsperedge = 1;
  parthread = 0;

  elsizeweight = 0.2;
  giveuptol2d = 200;
  giveuptol = 10;
  maxoutersteps = 10;
  starshapeclass = 5;
  baseelnp = 0;
  sloppy = 1;

  badellimit = 175;
  check_impossible = 0;
  secondorder = 0;
}

void MeshingParameters :: Print (ostream & ost) const
{
  ost << "Meshing parameters: " << endl
      << "optimize3d = " << optimize3d << endl
      << "optsteps3d = " << optsteps3d << endl
      << " optimize2d = " <<  optimize2d << endl
      << " optsteps2d = " <<  optsteps2d << endl
      << " opterrpow = " <<  opterrpow << endl
      << " blockfill = " <<  blockfill << endl
      << " filldist = " <<  filldist << endl
      << " safety = " <<  safety << endl
      << " relinnersafety = " <<  relinnersafety << endl
      << " uselocalh = " <<  uselocalh << endl
      << " grading = " <<  grading << endl
      << " delaunay = " <<  delaunay << endl
      << " maxh = " <<  maxh << endl;
  if(meshsizefilename)
    ost << " meshsizefilename = " <<  meshsizefilename << endl;
  else
    ost << " meshsizefilename = NULL" << endl;
  ost << " startinsurface = " <<  startinsurface << endl
      << " checkoverlap = " <<  checkoverlap << endl
      << " checkchartboundary = " <<  checkchartboundary << endl
      << " curvaturesafety = " <<  curvaturesafety << endl
      << " segmentsperedge = " <<  segmentsperedge << endl
      << " parthread = " <<  parthread << endl
      << " elsizeweight = " <<  elsizeweight << endl
      << " giveuptol2d = " <<  giveuptol2d << endl
      << " giveuptol = " <<  giveuptol << endl
      << " maxoutersteps = " <<  maxoutersteps << endl
      << " starshapeclass = " <<  starshapeclass << endl
      << " baseelnp        = " <<  baseelnp        << endl
      << " sloppy = " <<  sloppy << endl
      << " badellimit = " <<  badellimit << endl
      << " secondorder = " <<  secondorder << endl
      << " elementorder = " <<  elementorder << endl
      << " quad = " <<  quad << endl
      << " inverttets = " <<  inverttets << endl
      << " inverttrigs = " <<  inverttrigs << endl;
}

void MeshingParameters :: CopyFrom(const MeshingParameters & other)
{
  //strcpy(optimize3d,other.optimize3d); 
  optimize3d = other.optimize3d;
  optsteps3d = other.optsteps3d;
  //strcpy(optimize2d,other.optimize2d); 
  optimize2d = other.optimize2d;
  optsteps2d = other.optsteps2d;
  opterrpow = other.opterrpow;
  blockfill = other.blockfill;
  filldist = other.filldist;
  safety = other.safety;
  relinnersafety = other.relinnersafety;
  uselocalh = other.uselocalh;
  grading = other.grading;
  delaunay = other.delaunay;
  maxh = other.maxh;
  //strcpy(const_cast<char*>(meshsizefilename), other.meshsizefilename);
  //const_cast<char*>(meshsizefilename) = other.meshsizefilename; //???
  startinsurface = other.startinsurface;
  checkoverlap = other.checkoverlap;
  checkoverlappingboundary = other.checkoverlappingboundary;
  checkchartboundary = other.checkchartboundary;
  curvaturesafety = other.curvaturesafety;
  segmentsperedge = other.segmentsperedge;
  parthread = other.parthread;
  elsizeweight = other.elsizeweight;
  giveuptol2d = other.giveuptol2d;
  giveuptol = other.giveuptol;
  maxoutersteps = other.maxoutersteps;
  starshapeclass = other.starshapeclass;
  baseelnp = other.baseelnp;       
  sloppy = other.sloppy;
  badellimit = other.badellimit;
  secondorder = other.secondorder;
  elementorder = other.elementorder;
  quad = other.quad;
  inverttets = other.inverttets;
  inverttrigs = other.inverttrigs;
}


DebugParameters :: DebugParameters ()
{
  slowchecks = 0;
  haltsuccess = 0;
  haltnosuccess = 0;
  haltlargequalclass = 0;
  haltsegment = 0;
  haltsegmentp1 = 0;
  haltsegmentp2 = 0;
};



}

