#if !defined( BL_ELEMENT_H )
#define BL_ELEMENT_H

#include <vector>
#include <fstream>
#include <set>

#include "../../config.h"
#if defined(_NO_STD_MINMAX) || defined(WIN32)
	#include "minmaxpatch.h"
#endif

#include "MeshNode.h"

class Element
{
public:
	friend std::ostream& operator<< (std::ostream& o, const Element& A);
	
	Element();
	Element( const int t ) { elementTag = t; }
	
	~Element() { };
	
	void newTag();
	int elementId() const { return elementTag; }
	
	void setBody( const int bd ) { body = bd; }
	int partOfBody() const { return body; }
	bool isSameBody( const int bd ) const { return body == bd; }
	
	int size() const { return sz; }
	
	Node* nodeAt( const int t) { return nodes[t]; }
	
	bool isBoundaryConnector( const std::set< std::pair< int, int > > &links ) const;

	void upgrade(std::vector<Node *> &newNodes)
	{
		Node **n = new Node*[sz+newNodes.size()];

		int i;
		for (i = 0; i < sz; i++)
			n[i] = nodes[i];
		for (i = 0; i < newNodes.size(); i++)
			n[sz+i] = newNodes[i];

		delete[] nodes;
		nodes = n;

		sz += newNodes.size();
	}

	virtual int elementType() const = 0;
	
protected:
	int elementTag;
	int sz;
	int body;
	Node** nodes;
};

std::ostream& operator<< (std::ostream& o, const Element& A);

#endif /* BL_ELEMENT_H */
