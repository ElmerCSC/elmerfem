#ifndef FILE_CSGEOM
#define FILE_CSGEOM

/**************************************************************************/
/* File:   csgeom.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Nov. 97                                                    */
/**************************************************************************/

/**
  Constructive Solid Geometry
*/


class TriangleApproximation;
class TATriangle;


/**
   A top level object is an entity to be meshed.
   I can be either a solid, or one surface patch of a solid.
 */
class TopLevelObject
{
  Solid * solid;
  Surface * surface;

  double red, blue, green;
  bool visible, transp;
  double maxh;
  string material;
  int layer;
  int bc;     // for surface patches, only
  string bcname;

public:
  TopLevelObject (Solid * asolid,
		  Surface * asurface = NULL);

  const Solid * GetSolid() const { return solid; }
  Solid * GetSolid() { return solid; }

  const Surface * GetSurface () const { return surface; }
  Surface  * GetSurface () { return surface; }

  void GetData (ostream & ost);
  void SetData (istream & ist);

  void SetMaxH (double amaxh) { maxh = amaxh; } 
  double GetMaxH () const { return maxh; }

  void SetRGB (double ared, double agreen, double ablue)
  {
    red = ared;
    green = agreen;
    blue = ablue;
  }

  double GetRed () const { return red; }
  double GetGreen () const { return green; }
  double GetBlue () const { return blue; }

  void SetTransparent (bool atransp) 
  { transp = atransp; }
  bool GetTransparent () const { return transp; }

  void SetVisible (bool avisible)
  { visible = avisible; }
  bool GetVisible () const { return visible; }

  const string GetMaterial () const { return material; }
  void SetMaterial (const string & mat) { material = mat; }

  int GetLayer () const { return layer; }
  void SetLayer (int alayer) { layer = alayer; }

  void SetBCProp (int abc) { bc = abc; }
  int GetBCProp () const { return bc; }

  void SetBCName (string abc) { bcname = abc; }
  const string GetBCName () const { return bcname; }
};


/**
   CSGeometry has the whole geometric information
 */
class CSGeometry
{
private:
  /// all surfaces
  SYMBOLTABLE<Surface*> surfaces;

public:
  /// primitive of surface
  ARRAY<const Primitive*> surf2prim;

private:
  ARRAY<Surface*> delete_them;

  /// all named solids
  SYMBOLTABLE<Solid*> solids;

  /// all 2d splinecurves
  SYMBOLTABLE< SplineGeometry<2>* > splinecurves2d;
  /// all 3d splinecurves
  SYMBOLTABLE< SplineGeometry<3>* > splinecurves3d;

  /// all top level objects: solids and surfaces
  ARRAY<TopLevelObject*> toplevelobjects;

  /// additional points specified by user
  ARRAY<Point<3> > userpoints;
  ARRAY<double> userpoints_ref_factor;

  mutable ARRAY<Point<3> > identpoints;

  /// triangular approximation of top level objects
  ARRAY<TriangleApproximation*> triapprox;

  /// increment, if geometry is changed
  static int changeval;
  
  /// bounding box of geometry
  Box<3> boundingbox;

  /// bounding box, if not set by input file
  static Box<3> default_boundingbox;

  /// identic surfaces are stored by pair of indizes, val = inverse
  INDEX_2_HASHTABLE<int> identicsurfaces;
  ARRAY<int> isidenticto;
  /// identification of boundaries (periodic, thin domains, ...)

  double ideps;


  /// filename of inputfile
  string filename;


public:
  CSGeometry ();
  CSGeometry (const string & afilename);
  ~CSGeometry ();

  void Clean ();

  void Save (ostream & ost);
  void Load (istream & ist);

  void SaveSurfaces (ostream & out);
  void LoadSurfaces (istream & in);

  int GetChangeVal() { return changeval; }
  void Change() { changeval++; }

  void AddSurface (Surface * surf);
  void AddSurface (char * name, Surface * surf);
  void AddSurfaces (Primitive * prim);

  int GetNSurf () const { return surfaces.Size(); }
  const Surface * GetSurface (const char * name) const;
  const Surface * GetSurface (int i) const
  { return surfaces[i]; }

  void SetSolid (const char * name, Solid * sol);
  const Solid * GetSolid (const char * name) const;
  const Solid * GetSolid (const string & name) const;
  int GetNSolids () const { return solids.Size(); }
  const Solid * GetSolid (int i) const { return solids[i]; }
  const SYMBOLTABLE<Solid*> & GetSolids () const { return solids; }


  void SetSplineCurve (const char * name, SplineGeometry<2> * spl);
  void SetSplineCurve (const char * name, SplineGeometry<3> * spl);
  const SplineGeometry<2> * GetSplineCurve2d (const string & name) const;
  const SplineGeometry<3> * GetSplineCurve3d (const string & name) const;
    

  void SetFlags (const char * solidname, const Flags & flags);


  int GetNTopLevelObjects () const
  { return toplevelobjects.Size(); }
  int SetTopLevelObject (Solid * sol, Surface * surf = NULL);
  void GetTopLevelObject (int nr, Solid *& sol, Surface *& surf)
  {
    sol = toplevelobjects[nr]->GetSolid();
    surf = toplevelobjects[nr]->GetSurface();
  }
  void GetTopLevelObject (int nr, const Solid *& sol, const Surface *& surf) const
  {
    sol = toplevelobjects[nr]->GetSolid();
    surf = toplevelobjects[nr]->GetSurface();
  }

  TopLevelObject * GetTopLevelObject (const Solid * sol, const Surface * surf = NULL);
  TopLevelObject * GetTopLevelObject (int nr)
  { return toplevelobjects[nr]; }
  const TopLevelObject * GetTopLevelObject (int nr) const
  { return toplevelobjects[nr]; }
  void RemoveTopLevelObject (Solid * sol, Surface * surf = NULL); 


  void AddUserPoint (const Point<3> & p, double ref_factor = 0)
  { userpoints.Append (p); userpoints_ref_factor.Append (ref_factor); }
  int GetNUserPoints () const
  { return userpoints.Size(); }
  const Point<3> & GetUserPoint (int nr) const
  { return userpoints[nr]; }
  double GetUserPointRefFactor (int nr) const
  { return userpoints_ref_factor[nr]; }
  
  void AddIdentPoint (const Point<3> & p) const
  { identpoints.Append(p);}
  int GetNIdentPoints (void) const
  { return identpoints.Size();}
  const Point<3> & GetIdentPoint(int nr) const
  { return identpoints[nr]; }
  void DeleteIdentPoints(void) const
  { identpoints.DeleteAll();}


  // quick implementations:
  ARRAY<SingularFace*> singfaces;
  ARRAY<SingularEdge*> singedges;
  ARRAY<SingularPoint*> singpoints;
  ARRAY<Identification*> identifications;

  int GetNIdentifications (void) const { return identifications.Size(); }
  void AddIdentification (Identification * ident);


  ///
  void CalcTriangleApproximation(const Box<3> & boundingbox,
				 double detail, double facets);

  ///
  void FindIdenticSurfaces (double eps);
  ///
  void GetSurfaceIndices (const Solid * sol, 
			  const BoxSphere<3> & box, 
			  ARRAY<int> & locsurf) const;
  ///
  void GetIndependentSurfaceIndices (const Solid * sol, 
				     const BoxSphere<3> & box, 
				     ARRAY<int> & locsurf) const;
  ///
  void GetIndependentSurfaceIndices (const Solid * sol, 
				     const Point<3> & p, Vec<3> & v,
				     ARRAY<int> & locsurf) const;
  ///
  void GetIndependentSurfaceIndices (ARRAY<int> & locsurf) const;

  ///
  int GetSurfaceClassRepresentant (int si) const
    { return isidenticto[si]; }

  ///
  const TriangleApproximation * GetTriApprox (int msnr)
  {
    if (msnr < triapprox.Size())
      return triapprox[msnr];
    return 0;
  }
  

  void IterateAllSolids (SolidIterator & it, bool only_once = false);

  void RefineTriangleApprox (Solid * locsol, 
			     int surfind,
			     const BoxSphere<3> & box, 
			     double detail,
			     const TATriangle & tria, 
			     TriangleApproximation & tams,
			     IndexSet & iset);

  const Box<3> & BoundingBox () const { return boundingbox; }

  void SetBoundingBox (const Box<3> & abox)
  {
    boundingbox = abox;
  }


  static void SetDefaultBoundingBox (const Box<3> & abox)
  {
    default_boundingbox = abox;
  }

  double MaxSize () const;

  void SetIdEps(double eps){ideps = eps;}
  double GetIdEps(void) const {return ideps;}

  class BCModification {
  public:
    int si;
    int tlonr;
    int bcnr;
    string * bcname;
  };

  ARRAY<BCModification> bcmodifications;

};
#endif

