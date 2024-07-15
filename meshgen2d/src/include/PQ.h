#if !defined( MESH_PQ_H )
#define MESH_PQ_H

#include <list>
#include "Vertex.h"


#include <vector>
class pq
{
public:
	pq(int chunk = 1024);
	int size();
	Vertex *first();
	
	void insert( Vertex *vtx );
	void removeRelevants( std::list< Vertex* >& vl );
	void remove( Vertex *vtx );
	
	void debug();
	
private:
	void upheap( const int p );
	void downheap( const int p );
	void pq_reserve(int sz);
	
	Vertex** store;
	int last;
	int pq_size;
};

#endif /* MESH_PQ_H */
