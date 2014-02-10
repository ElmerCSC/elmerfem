#ifndef LOCALH
#define LOCALH

/**************************************************************************/
/* File:   localh.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   29. Jan. 97                                                    */
/**************************************************************************/




/// box for grading
class GradingBox
{
  /*
  /// xmin
  float x1[3];
  /// xmax
  float x2[3];
  */
  /// xmid
  float xmid[3];
  /// half edgelength
  float h2;
  ///
  GradingBox * childs[8];
  ///
  GradingBox * father;
  ///
  double hopt;
  ///
  struct 
  {
    unsigned int cutboundary:1;
    unsigned int isinner:1;
    unsigned int oldcell:1;
    unsigned int pinner:1;
  } flags;
public:
  ///
  GradingBox (const double * ax1, const double * ax2);
  ///
  void DeleteChilds();
  ///
  friend class LocalH;


  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};



/**
   Control of 3D mesh grading
 */
class LocalH 
{
  ///
  GradingBox * root;
  ///
  double grading;
  ///
  ARRAY<GradingBox*> boxes;
  ///
  Box3d boundingbox;
public:
  ///
  LocalH (const Point3d & pmin, const Point3d & pmax, double grading);
  ///
  ~LocalH();
  ///
  void Delete();
  ///
  void SetGrading (double agrading) { grading = agrading; }
  ///
  void SetH (const Point3d & x, double h);
  ///
  double GetH (const Point3d & x) const;
  /// minimal h in box (pmin, pmax)
  double GetMinH (const Point3d & pmin, const Point3d & pmax) const;

  /// mark boxes intersecting with boundary-box
  void CutBoundary (const Point3d & pmin, const Point3d & pmax)
    { CutBoundaryRec (pmin, pmax, root); }

  /// find inner boxes
  void FindInnerBoxes ( // int (*sameside)(const Point3d & p1, const Point3d & p2),
		       class AdFront3 * adfront,
		       int (*testinner)(const Point3d & p1));

  /// clears all flags 
  void ClearFlags ()
    { ClearFlagsRec(root); }

  /// widen refinement zone
  void WidenRefinement ();

  /// get points in inner elements
  void GetInnerPoints (ARRAY<Point3d> & points);

  /// get points in outer closure
  void GetOuterPoints (ARRAY<Point3d> & points);

  ///
  void Convexify ();
  ///
  int GetNBoxes () { return boxes.Size(); } 
  const Box3d & GetBoundingBox () const
  { return boundingbox; }
  ///
  void PrintMemInfo (ostream & ost) const;
private:
  /// 
  double GetMinHRec (const Point3d & pmin, const Point3d & pmax,
		     const GradingBox * box) const;
  ///
  void CutBoundaryRec (const Point3d & pmin, const Point3d & pmax,
		       GradingBox * box);

  ///
  void FindInnerBoxesRec ( int (*inner)(const Point3d & p),
			   GradingBox * box);

  ///
  void FindInnerBoxesRec2 (GradingBox * box,
			   class AdFront3 * adfront,
			   ARRAY<Box3d> & faceboxes,
			   ARRAY<int> & finds, int nfinbox);


  ///
  void SetInnerBoxesRec (GradingBox * box);

  ///
  void ClearFlagsRec (GradingBox * box);
  
  ///
  void ConvexifyRec (GradingBox * box);
};


#endif
