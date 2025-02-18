#ifndef SSVORONOISEGMENT_H
#define SSVORONOISEGMENT_H

#include "Connect.h"
#include "PQ.h"

class SSVoronoiSegment : public Connect
{
public:
	SSVoronoiSegment( const int t ) : Connect(t) { }
	SSVoronoiSegment( const int t, BGMesh *bgMesh ) :
	Connect(t, bgMesh) { }

	void makeExplicitSeed( const int t1, const int t2, const int t3 )
	{
		explicitSeed = true;
		seedTags[0] = t1;
		seedTags[1] = t2;
		seedTags[2] = t3;
	}

	void makeImplicitSeed( GeometryEdge *ed, int dir )
	{
		explicitSeed = false;
		baseEdge = ed;
		baseDirection = dir;
	}

protected:
	pq actives;

	bool explicitSeed;
	int seedTags[3];
	MeshNode *seed[3];
	int baseDirection;
	GeometryEdge *baseEdge;
};

#endif
