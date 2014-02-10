#if !defined( MESH_BGVERTEX_H )
#define MESH_BGVERTEX_H

#include "Vertex.h"
#include<vector>
#include "GeometryNode.h"

class BGVertex : public Vertex
{
public:
	BGVertex() : Vertex() { }
	void initInterpolation( );
	double interpolate( const double ix, const double iy );
	
	double coeff[3];
};

#endif /* MESH_BGVERTEX_H */
