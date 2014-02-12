#ifndef FILE_STLGEOM
#define FILE_STLGEOM

/**************************************************************************/
/* File:   stlgeom.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   26. Jul. 99                                                    */
/**************************************************************************/

/**
   STL Geometry


   Terminology:
   
   Point ... coordinates of STL triangles
   Triangle  (short Trig)  STL triangle
   TopEdge .... edge in topology, boundary of STL triangles (many)
   Edge .... Edges which will occur in the mesh (confirmed edges, less)
*/


#include <gprim.hpp>
#include <meshing.hpp>



namespace netgen
{
  extern int IsInArray(int n, const ARRAY<int>& ia);
  extern int AddIfNotExists(ARRAY<int>& list, int x);


#include "stltopology.hpp"
#include "stltool.hpp"
#include "stlline.hpp"







  class STLEdgeDataList
  {
    ARRAY<int> storedstatus;
    STLTopology & geom;
  public:
  
    STLEdgeDataList(STLTopology & ageom);
    ~STLEdgeDataList();

    void Store ();
    void Restore ();

    void SetSize(int /* size */) { };
    void Clear() { };
    int Size() const { return geom.GetNTE(); }
    const STLTopEdge & Get(int i) const { return geom.GetTopEdge(i); }
    STLTopEdge & Elem(int i) { return geom.GetTopEdge(i); }

    int GetNEPP(int pn) const {return geom.NTopEdgesPerPoint(pn); }
    int GetEdgePP(int pn, int vi) const {return geom.TopEdgePerPoint(pn, vi);};

    //void AddEdgePP(int pn, int vn) { } ;

    void ResetAll();
    void ChangeStatus(int status1, int status2);

    int GetEdgeNum(int np1, int np2) const
    { return geom.GetTopEdgeNum (np1, np2); }

    int GetNConfEdges() const;

    void Write(ofstream& of) const;
    void Read(ifstream& ifs);

    void BuildLineWithEdge(int ep1, int ep2, ARRAY<twoint>& line);
    void BuildClusterWithEdge(int ep1, int ep2, ARRAY<twoint>& line);

    int GetNEPPStat(int p, int status) const;
    int GetNConfCandEPP(int p) const;
  };






  class STLGeometry : public STLTopology
  {
    // edges to be meshed:
    ARRAY<STLEdge> edges;
    //edges per point
    TABLE<int> edgesperpoint;

    // line: a connection of edges
    ARRAY<STLLine*> lines;
    ARRAY<int> lineendpoints; //per geometrypoint, 1 = is endpoint; 0 = no endpoint,

    ARRAY<Vec3d> normals; //normals belong to points!

    ARRAY<twoint> externaledges;

    int undoexternaledges;
    ARRAY<twoint> storedexternaledges;

    STLEdgeDataList * edgedata;
    //  STLEdgeDataList edgedata_store;
    int calcedgedataanglesnew;

    int edgedatastored;



    int facecnt; 
    //meshpoint is only set, if an edge is at this point!!!

    ARRAY<int> vicinity; //is one, if a triangle belongs to vicinity (eg. of selecttrig)
    ARRAY<int> markedtrigs; //is one, if a triangle belongs to marked triangles (calcdirtystrigs)
    ARRAY<Point3d> markedsegs; //every pointpair is a segment!!!  
    ARRAY<twoint> selectedmultiedge;


    //spiralpoints:
    ARRAY<int> spiralpoints;
    //
    ARRAY<STLChart*> atlas;
    //marks all already charted trigs with chartnumber
    ARRAY<int> chartmark; 
    //outerchartspertrig, ascending sorted
    TABLE<int> outerchartspertrig;


    //for meshing and project:
    ARRAY<int> meshcharttrigs; //per trig: 1=belong to chart, 0 not
    int meshchart;

    ARRAY<int> ha_points;  // help array, np long, filled with 0 


    // sharp geometric edges not declared as edges
    // (not considered for spiral check)
    INDEX_2_HASHTABLE<int> * smoothedges;


    //transformation:
    Vec<3> meshtrignv;
    Vec<3> ex, ey, ez;
    Point<3> p1;

  public:
    int edgesfound;
    int surfacemeshed;
    int surfaceoptimized;
    int volumemeshed;

    int trigsconverted; //when STLTriangles exist -> 1

    //for selecting nodes
    //int selecttrig, nodeofseltrig;

    //only for testing;
    ARRAY<STLLine*> meshlines;
    ARRAY<Point3d> meshpoints;

  public:
    STLGeometry();
    virtual ~STLGeometry();


    void Clear();



    void STLInfo(double* data);
    //stldoctor:
    void SmoothNormals();
    void MarkNonSmoothNormals();

    void CalcEdgeData();
    void CalcEdgeDataAngles();

    const STLEdgeDataList& EdgeDataList() const {return *edgedata;}

    void UndoEdgeChange();
    void StoreEdgeData();
    void RestoreEdgeData();

    //void ClearSelectedMultiEdge() {selectedmultiedge.SetSize(0);}
    //void AddSelectedMultiEdge(twoint ep) {selectedmultiedge.Append(ep);}
    //int SelectedMultiEdgeSize() {return selectedmultiedge.Size();}
    const ARRAY<twoint>& SelectedMultiEdge() {return selectedmultiedge;}
    twoint GetNearestSelectedDefinedEdge();
    void BuildSelectedMultiEdge(twoint ep);
    void BuildSelectedEdge(twoint ep);
    void BuildSelectedCluster(twoint ep);

    void ImportEdges();
    void AddEdges(const ARRAY<Point<3> >& eps);
    void ExportEdges();
    void LoadEdgeData(const char* file);
    void SaveEdgeData(const char* file);
    //  void SetEdgeAtSelected(int mode);
  

    void STLDoctorConfirmEdge();
    void STLDoctorCandidateEdge();
    void STLDoctorExcludeEdge();
    void STLDoctorUndefinedEdge();

    void STLDoctorSetAllUndefinedEdges();
    void STLDoctorEraseCandidateEdges();
    void STLDoctorConfirmCandidateEdges();
    void STLDoctorConfirmedToCandidateEdges();

    void STLDoctorDirtyEdgesToCandidates();
    void STLDoctorLongLinesToCandidates();

    void UndoExternalEdges();
    void StoreExternalEdges();
    void RestoreExternalEdges();

    void ImportExternalEdges(const char * filename);  // Flame edges, JS
    //  void LoadExternalEdges();

    void BuildExternalEdgesFromEdges();
    void SaveExternalEdges();
    void AddExternalEdgeAtSelected();
    void AddClosedLinesToExternalEdges();
    void AddLongLinesToExternalEdges();
    void AddAllNotSingleLinesToExternalEdges();
    void STLDoctorBuildEdges();
    void AddExternalEdgesFromGeomLine();
    void DeleteDirtyExternalEdges();
    void DeleteExternalEdgeAtSelected();
    void DeleteExternalEdgeInVicinity();
    void AddExternalEdge(int p1, int p2);
    void DeleteExternalEdge(int p1, int p2);
    int IsExternalEdge(int p1, int p2);
    int NOExternalEdges() const {return externaledges.Size();}
    twoint GetExternalEdge(int i) const {return externaledges.Get(i);}

    void DestroyDirtyTrigs();
    void CalcNormalsFromGeometry();
    void MoveSelectedPointToMiddle();
    void NeighbourAnglesOfSelectedTrig();
    void PrintSelectInfo();
    void ShowSelectedTrigChartnum();
    void ShowSelectedTrigCoords();
    void SmoothGeometry ();


    void LoadMarkedTrigs();
    void SaveMarkedTrigs();
    void ClearMarkedSegs() {markedsegs.SetSize(0);}
    void AddMarkedSeg(const Point<3> & ap1, const Point<3> & ap2) 
    {
      markedsegs.Append(ap1);markedsegs.Append(ap2);
    }

    void GetMarkedSeg(int i, Point<3> & ap1, Point<3> & ap2) 
    {
      ap1=markedsegs.Get(i*2-1); 
      ap2=markedsegs.Get(i*2);
    }
    int GetNMarkedSegs() {return markedsegs.Size()/2;}
    void CalcVicinity(int starttrig);
    void GetVicinity(int starttrig, int size, ARRAY<int>& vic);

    int Vicinity(int trig) const;

    void InitMarkedTrigs();
    void MarkDirtyTrigs();
    void SmoothDirtyTrigs();
    void GeomSmoothRevertedTrigs();
    void MarkRevertedTrigs();
    double CalcTrigBadness(int i);
    int IsMarkedTrig(int trig) const;
    void SetMarkedTrig(int trig, int num);
    void MarkTopErrorTrigs ();

    //Selected triangle
    void SetSelectTrig(int trig);
    int GetSelectTrig() const;
    void SetNodeOfSelTrig(int n);
    int GetNodeOfSelTrig() const;


    int AddNormal(const Vec3d& n) {return normals.Append(n);}
    const Vec3d & GetNormal(int nr) const {return normals.Get(nr);}
    void SetNormal(int nr, const Vec3d& n) {normals.Elem(nr) = n;}

    int AddEdge(const STLEdge& v) {return edges.Append(v);}
    int AddEdge(int p1, int p2);

    STLEdge GetEdge(int nr) {return edges.Get(nr);}
    int GetNE() {return edges.Size();}

    double Area();

    double GetAngle(int t1, int t2);
    double GetGeomAngle(int t1, int t2);
    //if triangles t1 and t2 touch, return 1 and in p1, p2 the touching points
    //int TrigsTouch(int t1, int t2, int& p1, int& p2);


  
    ///

    ///ReadTriangle->STLTriangle, initialise some important variables, always after load!!!
    virtual void InitSTLGeometry (const ARRAY<STLReadTriangle> & readtrigs);
    virtual void TopologyChanged(); //do some things, if topology changed!
    int CheckGeometryOverlapping();

    //get NO edges per point
    int GetEPPSize() const {return edgesperpoint.Size();};
    int GetNEPP(int pn) 
    {
      if (edgesperpoint.Size() == 0) {BuildEdgesPerPoint();}
      return edgesperpoint.EntrySize(pn);
    };
    int GetEdgePP(int pn, int vi)
    {
      if (edgesperpoint.Size() == 0) {BuildEdgesPerPoint();}
      return edgesperpoint.Get(pn,vi);
    };
    void AddEdgePP(int pn, int vn) {edgesperpoint.Add1(pn,vn);};
    //von 2 punkten ermitteln, ob sie eine Kante sind
    int IsEdge(int p1, int p2);
    int IsEdgeNum(int p1, int p2);

    ///Build EdgeSegments
    void ClearEdges();
    void BuildEdges();
    void BuildEdgesPerPoint();
    void UseExternalEdges();


    void FindEdgesFromAngles();
    void CalcFaceNums();
    int GetNOBodys();
    int GetNOFaces() {return facecnt;}
    void LinkEdges();

    void AddConeAndSpiralEdges();
    void AddFaceEdges(); //each face should have at least one starting edge (outherwise it won't be meshed)

    void GetDirtyChartTrigs(int chartnum, STLChart& chart, const ARRAY<int>& outercharttrigs, 
			    ARRAY<int>& chartpointchecked, ARRAY<int>& dirtytrigs);

    void ClearSpiralPoints();
    void SetSpiralPoint(int pn) {spiralpoints.Elem(pn) = 1;};
    int GetSpiralPoint(int pn) const {return spiralpoints.Get(pn);};

    void GetSortedTrianglesAroundPoint(int p, int starttrig, ARRAY<int>& trigs);

    // smooth edges: sharp geometric edges not declared as edges
    void BuildSmoothEdges ();
    int IsSmoothEdge (int pi1, int pi2) const;


    //make charts with regions of a max. angle
    void MakeAtlas(class Mesh & mesh);

    //outerchartspertrig, sorted!
    int GetOCPTSize() const {return outerchartspertrig.Size();};
    int GetNOCPT(int tn) const {return outerchartspertrig.EntrySize(tn);};
    int GetOCPT(int tn, int vi) const {return outerchartspertrig.Get(tn,vi);};
    void SetOCPT(int tn, int vi, int ocn) {outerchartspertrig.Set(tn,vi,ocn);};
    void AddOCPT(int tn, int ocn) {outerchartspertrig.Add1(tn, ocn);};
    int TrigIsInOC(int tn, int ocn) const;
 
    //get chart number of a trig or 0 if unmarked
    int GetChartNr(int i) const;
    int GetMarker(int i) const 
    { return chartmark.Get(i); }
    void SetMarker(int nr, int m);
    int GetNOCharts() const;
    //get a chart from atlas
    const STLChart& GetChart(int nr) const;
    STLChart& GetChart(int nr) {return *(atlas.Get(nr));};
    int AtlasMade() const;
  
    void GetInnerChartLimes(ARRAY<twoint>& limes, int chartnum);

    //FOR MESHING
    int GetMeshChartNr () { return meshchart; }
    void GetMeshChartBoundary (ARRAY<Point2d > & points,
			       ARRAY<Point3d > & points3d,
			       ARRAY<INDEX_2> & lines, double h);


    Point<3> PointBetween(const Point<3> & p1, int t1, const Point<3> & p2, int t2);

    //select triangles in meshcharttrigs of actual (defined by trig) whole chart
    void PrepareSurfaceMeshing();
    //
    void DefineTangentialPlane(const Point<3> & ap1, const Point<3> & ap2, int trig);
    //
    void SelectChartOfTriangle (int trignum);
    //
    void SelectChartOfPoint (const Point<3> & p);
    //
    const Vec<3> & GetChartNormalVector () const { return meshtrignv; }

    // list of trigs
    void ToPlane (const Point<3> & locpoint, int * trigs, Point<2> & plainpoint, 
		  double h, int& zone, int checkchart);
    //return 0, wenn alles OK, 1 sonst
    int FromPlane (const Point<2> & plainpoint, Point<3> & locpoint, double h);
  
    //get nearest point in actual chart and return any triangle where it lies on
    int ProjectNearest(Point<3> & p3d) const;
    //project point with normal nv from last define tangential plane

    int LastTrig() const;
    int Project(Point<3> & p3d) const;
    int ProjectOnWholeSurface (Point<3> & p3d) const;

    int GetNLines() const {return lines.Size();}
    int AddLine(STLLine* line) {return lines.Append(line);}
    STLLine* GetLine(int nr) const {return lines.Get(nr);}
    int GetLineP(int lnr, int pnr) const {return lines.Get(lnr)->PNum(pnr);}
    int GetLineNP(int nr) const {return lines.Get(nr)->NP();}

    void SetLineEndPoint(int pn);
    int IsLineEndPoint(int pn);
    int LineEndPointsSet() const {return lineendpoints.Size() == GetNP();}
    void ClearLineEndPoints();

    void RestrictLocalH(class Mesh & mesh, double gh);
    void RestrictLocalHCurv(class Mesh & mesh, double gh);
    void RestrictHChartDistOneChart(int chartnum, ARRAY<int>& acttrigs, class Mesh & mesh, 
				    double gh, double fact, double minh);

    friend class MeshingSTLSurface;
  };
 

#include "meshstlsurface.hpp"


  extern int STLMeshingDummy (STLGeometry* stlgeometry, Mesh*& mesh,
			      int perfstepsstart, int perfstepsend, char* optstring);


}
#endif
