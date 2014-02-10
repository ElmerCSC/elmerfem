#include "TriangleNELayer.h"
#include "TriangleElement.h"

#define LL grid[j*m+i+1]
#define LR grid[(j+1)*m+i+1]
#define UR grid[(j+1)*m+i]
#define UL grid[j*m+i]

void TriangleNELayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	makeGrid( allNodes, m, n ); 
	
	std::vector<BoundaryElement*> bels[4];
	for (i = 0; i < 4; i++)
		edges[i]->elements(bels[i], directions[i]);
	
	for( j = 0; j < (n - 1); ++j)
		for( i = 0; i < (m - 1); ++i)
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
		}
		
	delete [] grid;
}

void TriangleNWLayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	makeGrid( allNodes, m, n ); 
	
	std::vector<BoundaryElement*> bels[4];
	for (i = 0; i < 4; i++)
		edges[i]->elements(bels[i], directions[i]);
	
	for( j = 0; j < (n - 1); ++j)
		for( i = 0; i < (m - 1); ++i)
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}

	delete [] grid;
}

void TriangleUJNELayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	makeGrid( allNodes, m, n ); 
	
	std::vector<BoundaryElement*> bels[4];
	for (i = 0; i < 4; i++)
		edges[i]->elements(bels[i], directions[i]);
	
	for( j = 0; j < (n - 1); j += 2 )
	{
		for( i = m-2; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
		}
		for( i = m-3; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}
	}
	
	for( j = 1; j < (n - 1); j += 2 )
	{
		for( i = m-2; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}
		for( i = m-3; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
		}
	}

	delete [] grid;
}

void TriangleUJNWLayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	makeGrid( allNodes, m, n ); 
	
	std::vector<BoundaryElement*> bels[4];
	for (i = 0; i < 4; i++)
		edges[i]->elements(bels[i], directions[i]);
	
	for( j = 0; j < (n - 1); j += 2 )
	{
		for( i = m-2; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}
		for( i = m-3; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
		}
	}
	
	for( j = 1; j < (n - 1); j += 2 )
	{
		for( i = m-2; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
		}
		for( i = m-3; i >= 0; i -= 2 )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}
	}

	delete [] grid;
}

void TriangleFBNELayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	makeGrid( allNodes, m, n ); 
	
	std::vector<BoundaryElement*> bels[4];
	for (i = 0; i < 4; i++)
		edges[i]->elements(bels[i], directions[i]);
	
	for( j = 0; j < (n - 1); j += 2 )
		for( i = m-2; i >= 0; --i )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
		}

	for( j = 1; j < (n - 1); j += 2 )
		for( i = m-2; i >= 0; --i )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}

	delete [] grid;
}

void TriangleFBNWLayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	makeGrid( allNodes, m, n ); 
	
	std::vector<BoundaryElement*> bels[4];
	for (i = 0; i < 4; i++)
		edges[i]->elements(bels[i], directions[i]);
	
	for( j = 0; j < (n - 1); j += 2 )
		for( i = m-2; i >= 0; --i )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UL);
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == 0)
				bels[3][i]->setLeft(t->elementId());
			
			t = new TriangleElement( LR, UR, UL);
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
		}

	for( j = 1; j < (n - 1); j += 2 )
		for( i = m-2; i >= 0; --i )
		{
			TriangleElement *t;
			
			t = new TriangleElement( LL, LR, UR );
			allElements.push_back( t );
			
			if (i == m - 2)
				bels[0][j]->setLeft(t->elementId());
			if (j == n - 2)
				bels[1][m - 2 - i]->setLeft(t->elementId());
			
			t = new TriangleElement( LL, UR, UL );
			allElements.push_back( t );
			
			if (i == 0)
				bels[2][n - 2 - j]->setLeft(t->elementId());
		}

	delete [] grid;
}
