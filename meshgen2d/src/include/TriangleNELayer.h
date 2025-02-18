#if !defined( BL_TRIANGLENELAYER_H )
#define BL_TRIANGLENELAYER_H

#include "QuadLayer.h"

class Node;

class TriangleNELayer : public QuadLayer
{
public:
	TriangleNELayer(const int t) : QuadLayer(t) { }
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
};

class TriangleNWLayer : public QuadLayer
{
public:
	TriangleNWLayer(const int t) : QuadLayer(t) { }
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
};

class TriangleUJNELayer : public QuadLayer
{
public:
	TriangleUJNELayer(const int t) : QuadLayer(t) { }
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
};

class TriangleUJNWLayer : public QuadLayer
{
public:
	TriangleUJNWLayer(const int t) : QuadLayer(t) { }
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
};

class TriangleFBNELayer : public QuadLayer
{
public:
	TriangleFBNELayer(const int t) : QuadLayer(t) { }
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
};

class TriangleFBNWLayer : public QuadLayer
{
public:
	TriangleFBNWLayer(const int t) : QuadLayer(t) { }
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
};

#endif /* BL_TRIANGLENELAYER_H */
