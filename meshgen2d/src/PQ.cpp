#include "PQ.h"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string.h>

pq::
pq(int chunk)
{ last = -1; pq_size = 0; pq_reserve(chunk); }

int pq::
size()
{
	return last + 1;
}

Vertex * pq::
first()
{
	if( last < 0 ) return (Vertex *) 0;
	
	Vertex *v = store[0];
	
	if( last > 0 ) remove( v );
	else
	{
		v->setHeap( -1 );
		last = -1;
		store[0] = (Vertex *) 0;
	}
	
	return v;
}

void pq::
remove( Vertex *vtx )
{
	int h = vtx->atHeap();
	
	if( h == last )
	{
		vtx->setHeap( -1 );
		store[last] = (Vertex *) 0;
		last = last - 1;
		return;
	}
	
	vtx->setHeap( -1 );
	
	store[h] = store[last];
	store[h]->setHeap( h );
	
	store[last] = (Vertex *) 0;
	last = last - 1;
	
	if (h > 0 && store[(h - 1) / 2]->orderingValue() < store[h]->orderingValue())
		upheap( h );
	else
		downheap( h );
}

void pq::
insert( Vertex *vtx )
{
	last = last + 1;
	
	if(last == pq_size)
	{
		pq_reserve(2 * pq_size);
	}
	
	store[last] = vtx;
	store[last]->setHeap( last );
	
	upheap( last );
}

void pq::
removeRelevants( std::list< Vertex* >& vl )
{
	std::list< Vertex* >::iterator it;
	for( it = vl.begin(); it != vl.end(); ++it)
	{
		if( (*it)->isAtHeap() )
		{
			remove( *it );
		}
	}	
}

void pq::
debug()
{
	int i;
	std::cout << "Size: " << last + 1 << std::endl;
	/*
	for( i = 0; i <= last; ++i )
	{
	std::cout << store[i]->atHeap() << ' ' << store[i]->orderingValue() << std::endl;
	}
	*/
}

void pq::
upheap( const int p )
{
	int k = p;
	Vertex *v = store[p];
	
	while( k != 0 && store[ ( k - 1 ) / 2 ]->orderingValue() <= v->orderingValue() )
	{
		store[ k ] = store[ ( k - 1 ) / 2 ];
		store[ k ]->setHeap( k );
		
		k = ( k - 1 ) / 2;
	}
	store[ k ] = v;
	store[ k ]->setHeap( k );
}

void pq::
downheap( const int p )
{
	int k = p;
	Vertex *v = store[p];
	int N = last;
	
	while( ( k + 1 ) <= ( N + 1 ) / 2 )
	{
		int j = k + k + 1;
		if( j < N )
		{
			if( store[ j ]->orderingValue() <= store[ j + 1 ]->orderingValue() )
			{
				j = j + 1;
			}
		}
		if( store[ j ]->orderingValue() <= v->orderingValue() )
		{
			break;
		}
		store[ k ] = store[ j ];
		store[ k ]->setHeap( k );
		k = j;
	}
	store[ k ] = v;
	store[ k ]->setHeap( k );	
}

void pq::
pq_reserve(int sz)
{
	int osz = pq_size;
	pq_size = sz;
	
	if(osz == 0)
	{
		store = new Vertex*[pq_size];
	}
	else
	{
		Vertex **hold = store;
		store = new Vertex*[pq_size];
		memcpy((void *)store, (const void *)hold, osz * sizeof(Vertex *));
		delete [] hold;
	}
}

