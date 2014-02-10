#ifndef FILE_STLLINE
#define FILE_STLLINE


/**************************************************************************/
/* File:   stlline.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

class STLGeometry;
class STLTopology;

class STLEdge
{
public:
  int pts[2];
  int trigs[2]; //left and right trig

  STLEdge (const int * apts) {pts[0] = apts[0]; pts[1] = apts[1];}
  STLEdge (int v1, int v2) {pts[0] = v1; pts[1] = v2;}
  STLEdge () {pts[0]=0;pts[1]=0;}
  int PNum(int i) const {return pts[(i-1)];}

  int LeftTrig() const {return trigs[0];}
  int RightTrig() const {return trigs[1];}
  void SetLeftTrig(int i) {trigs[0] = i;}
  void SetRightTrig(int i) {trigs[1] = i;}
};

enum STL_ED_STATUS { ED_EXCLUDED, ED_CONFIRMED, ED_CANDIDATE, ED_UNDEFINED };
                       

/*

class STLEdgeData
{
public:
  //  float angle;
  int p1;
  int p2;
  int lt; //left trig
  int rt; //right trig
  //  int status;

  STLTopology * top;  // pointer to stl topology
  int topedgenr;  // number of corresponding topology edge

  STLEdgeData() {}; 
  STLEdgeData(float anglei, int p1i, int p2i, int lti, int rti) 
{
//     angle = anglei; 
p1 = p1i; p2 = p2i;
      lt = lti; rt = rti;
    }

  int GetStatus () const;
  void SetStatus (int stat);

  void SetExcluded() { SetStatus (ED_EXCLUDED); }
  void SetConfirmed() { SetStatus (ED_CONFIRMED); }
  void SetCandidate() { SetStatus (ED_CANDIDATE); }
  void SetUndefined() { SetStatus (ED_UNDEFINED); }

  int Excluded() const {return GetStatus() == ED_EXCLUDED;}
  int Confirmed() const {return GetStatus() == ED_CONFIRMED;}
  int Candidate() const {return GetStatus() == ED_CANDIDATE;}
  int Undefined() const {return GetStatus() == ED_UNDEFINED;}
  int ConfCand() const {return GetStatus() == ED_CONFIRMED || GetStatus() == ED_CANDIDATE;}

  float CosAngle() const; 

  void Write(ofstream& of) const;
  void Read(ifstream& ifs);
};

class STLEdgeDataList
{
private:
  INDEX_2_HASHTABLE<int> hashtab;
  ARRAY<STLEdgeData> edgedata;
  TABLE<int> edgesperpoint;
  
public:

  STLEdgeDataList():edgedata(),hashtab(1),edgesperpoint() {};
  const STLEdgeDataList& operator=(const STLEdgeDataList& edl); 
  void SetSize(int size) 
    {
      edgedata.SetSize(size);
      hashtab.SetSize(size);
      edgesperpoint.SetSize(size);
    }
  void Clear() {SetSize(0);}
  int Size() const {return edgedata.Size();}
  const STLEdgeData& Get(int i) const {return edgedata.Get(i);}
  STLEdgeData& Elem(int i) {return edgedata.Elem(i);}
  void Add(const STLEdgeData& ed, int i);

  int GetNEPP(int pn) const 
    {
      return edgesperpoint.EntrySize(pn);
    };
  int GetEdgePP(int pn, int vi) const
    {
      return edgesperpoint.Get(pn,vi);
    };
  void AddEdgePP(int pn, int vn) {edgesperpoint.Add(pn,vn);};

  void ResetAll();
  void ResetCandidates();
  void ConfirmCandidates();
  int GetEdgeNum(int np1, int np2) const;

  int GetNConfEdges() const;

  void Write(ofstream& of) const;
  void Read(ifstream& ifs);

  void BuildLineWithEdge(int ep1, int ep2, ARRAY<twoint>& line);

  int GetNEPPStat(int p, int status) const;
  int GetNConfCandEPP(int p) const;
};
*/
















//a line defined by several points (polyline)
class STLLine
{
private:
  const STLGeometry * geometry;
  ARRAY<int> pts;
  ARRAY<int> lefttrigs;
  ARRAY<int> righttrigs;
  ARRAY<double> dists;
  int split;

public:
  STLLine(const STLGeometry * ageometry);
  void AddPoint(int i) {pts.Append(i);}
  int PNum(int i) const {return pts.Get(i);}
  int NP() const {return pts.Size();}
  int GetNS() const;
  void GetSeg(int nr, int& p1, int& p2) const;
  double GetSegLen(const ARRAY<Point<3> >& ap, int nr) const;
  int GetLeftTrig(int nr) const;
  int GetRightTrig(int nr) const;
  double GetDist(int nr) const { return dists.Get(nr);};
  void GetBoundingBox (const ARRAY<Point<3> > & ap, Box<3> & box) const;

  void AddLeftTrig(int nr) {lefttrigs.Append(nr);}
  void AddRightTrig(int nr) {righttrigs.Append(nr);}
  void AddDist (double dist) {dists.Append(dist); }
  int StartP() const {return pts.Get(1);}
  int EndP() const {return pts.Get(pts.Size());}
    
  double GetLength(const ARRAY<Point<3> >& ap) const;

  //suche punkt in entfernung (in linienkoordinaten) dist
  //in index ist letzter punkt VOR dist (d.h. max pts.Size()-1)
  Point<3> GetPointInDist(const ARRAY<Point<3> >& ap, double dist, int& index) const;

  //return a meshed polyline
  STLLine* Mesh(const ARRAY<Point<3> >& ap, 
		ARRAY<Point3d>& mp, double ghi,
		class Mesh& mesh) const;

  void DoSplit() {split = 1;}
  int ShouldSplit() const {return split;}
};

#endif
