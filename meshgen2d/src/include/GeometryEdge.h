#if !defined( BL_GEOMETRYEDGE_H )
#define BL_GEOMETRYEDGE_H

#include <fstream>
#include <vector>
#include <map>

#include "GeometryNode.h"
#include "MeshNode.h"
#include "BoundaryElement.h"

enum edge_type { VIRTUAL = 0, OUTER, INNER };

class BGMesh;

typedef std::map< int, Node * > NodeMap;
typedef std::map< int, Node * >::iterator NodeMapIt;

class GeometryEdge
{
public:
	friend std::ostream& operator<< (std::ostream& o, const GeometryEdge& A);
	
	GeometryEdge( int nSegments = 0 );
	GeometryEdge( const int t, int nSegments = 0 );
	~GeometryEdge() { }
	
	void addNode(GeometryNode* a) { dots.push_back(a); }

	void addBGMesh(BGMesh *m) { bgMeshes.push_back(m); }
	
	virtual int size() { return nodes.size(); }
	
	virtual void discretize( NodeMap& allNodes );
	
	void exportNodes(std::vector<Node*>& strip, int direction);
	
	virtual void exportGeometryNodes(
		std::vector<GeometryNode*>& strip, const int direction)
	{
		if( direction > 0 ) 
			std::copy( dots.begin(), dots.end(), std::back_inserter( strip ) );
		else 
			std::copy( dots.rbegin(), dots.rend(), std::back_inserter( strip ) );
	}
	
	void elements(std::vector< BoundaryElement * >& strip, int direction);
	
	GeometryNode* base() { return dots.front(); }
	GeometryNode* end() { return dots.back(); }
	
	void setBoundaryTag( const int t ) { boundaryTag = t; }
	void makeOuter() { type = OUTER; }
	void makeInner() { type = INNER; }
	void makeVirtual() { type = VIRTUAL; }
	
	bool isVirtual() { return type == VIRTUAL; }
	bool isOuter() { return type == OUTER; }

	bool isConstant() { return segments > 0; }
	void setSegments(int s) { segments = s; }

	void midNodes( Node*& a, Node*& b, int dir );
	
	int tag;
	
	void getGridPoints(double ox, double oy, double cellsize, std::vector<int> &x, std::vector<int> &y, std::vector<double> &delta);
protected:
	void discretizeConstantSegment( int nSeg, NodeMap& allNodes, GeometryNode *from, GeometryNode *to );
	void discretizeGradedSegment( NodeMap& allNodes, GeometryNode *from, GeometryNode *to );

	double interpolate(double x, double y);
	
	std::vector< GeometryNode* > dots;
	std::vector<Node*> nodes;
	
	int segments;
	int boundaryTag;
	edge_type type;
	
	std::vector< BoundaryElement * > bels;

	std::vector< BGMesh * > bgMeshes;
};

std::ostream& operator<< (std::ostream& o, const GeometryEdge& A);

#endif /* BL_GEOMETRYEDGE_H */
