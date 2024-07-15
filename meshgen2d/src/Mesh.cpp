#include "../config.h"
#include <stdio.h>
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <math.h>
#include "Mesh.h"
#include "GeometryNode.h"
#include "GeometryEdge.h"
#include "Connect.h"
#include "BoundaryLayer.h"
#include "VoronoiVertex.h"
#include "VoronoiSegment.h"
#include "SSMFVoronoiSegment.h"
#include "SSSFVoronoiSegment.h"
#include "Border.h"
#include "Body.h"
#include "QuadLayer.h"
#include "TriangleNELayer.h"

#include "Node.h"
#include "Element.h"

#include <vector>

Mesh::Mesh()
{
}

void Mesh::
convertEdgeFormat( MeshParser& parser )
{
	int i, len;
	
	std::map<int, std::list<EdgeToken*> > nodeEdges;
	std::map<int, std::list<LayerToken*> > nodeLayers;
	std::map<int, int> edgeBody;

	std::vector<Layer *> layers;
	
	len = parser.nodes.size();
	for( i = 0; i < len; ++i )
	{
		NodeToken *nt = parser.nodes[i];
		GeometryNode *nd = new GeometryNode( nt->tag, nt->x, nt->y);
		nd->boundaryTag = nt->boundaryTag;
		// 	nd->setDelta( nt->delta );
		geometryNodes[ nd->tag ] = nd;
	}
	
	len = parser.edges.size();
	for( i = 0; i < len; ++i )
	{
		EdgeToken *et = parser.edges[i];
		
		GeometryEdge *ed = new GeometryEdge( et->tag, et->segCount );
		ed->setBoundaryTag(et->boundaryTag);

		int sz = et->nodes.size();
		for(int j = 0; j < sz; ++j)
		{
			nodeEdges[et->nodes[j]].push_back(et);
			
			GeometryNodeMap::iterator it = geometryNodes.find(et->nodes[j]);
			if (it == geometryNodes.end())
				blm_error("Nonexisting node in edge", et->tag);

			ed->addNode( it->second );
		}
		
		geometryEdges[ ed->tag ] = ed;
	}
	
	len = parser.bodies.size();
	for( i = 0; i < len; ++i )
	{
		BodyToken *bt = parser.bodies[i];
		
		Body *body = new Body( bt->tag, bt->parabolic );
		int layerLen = bt->layers.size();
		for( int m = 0; m < layerLen; ++m )
		{
			LayerToken *ly = bt->layers[m];
			
			if (ly->delta <= 0.0)
				ly->delta = bt->delta;

			BGMesh *bgmesh = NULL;
			Layer *layer;
			
			// Here comes the selection depending on type
			if( ly->type == QUAD_GRID )
				layer = new QuadLayer( bt->tag );
			else if( ly->type == TRIANGLE_NE_GRID )
				layer = new TriangleNELayer( bt->tag );
			else if( ly->type == TRIANGLE_NW_GRID )
				layer = new TriangleNWLayer( bt->tag );
			else if( ly->type == TRIANGLE_UJ_NE_GRID )
				layer = new TriangleUJNELayer( bt->tag );
			else if( ly->type == TRIANGLE_UJ_NW_GRID )
				layer = new TriangleUJNWLayer( bt->tag );
			else if( ly->type == TRIANGLE_FB_NE_GRID )
				layer = new TriangleFBNELayer( bt->tag );
			else if( ly->type == TRIANGLE_FB_NW_GRID )
				layer = new TriangleFBNWLayer( bt->tag );
			else
			{
				if( ly->type == CONNECT )
				{
					bgmesh = new BGTriangleMesh;
					layer = new Connect( bt->tag, bgmesh );
				}
				else if ( ly->type == BOUNDARY_MESH )
				{
					bgmesh = new BGTriangleMesh;
					layer = new BoundaryLayer( bt->tag, bgmesh );
				}
				else
				{
					if ( ly->bg->type == "Grid" )
						bgmesh = new BGGridMesh(2.0);
					else if ( ly->bg->type == "DenseGrid" )
						bgmesh = new BGGridMesh(1.0);
					else if ( ly->bg->type == "SparseGrid" )
						bgmesh = new BGGridMesh(3.0);
					else if ( ly->bg->type == "Delaunay" || ly->bg->type == "Explicit" )
						bgmesh = new BGTriangleMesh;
					else
						blm_error("Unknown background mesh type:", ly->bg->type.c_str());
					
					if( ly->bg->type == "Explicit" )
					{
						std::vector< GeometryNode * > bgnodes;
						std::vector< GeometryEdge * > dummy;

						for(int ind = 0; ind < ly->bg->nodes.size(); ++ind)
						{
							NodeToken nt = ly->bg->nodes[ind];
							GeometryNode *nd = new GeometryNode( 0, nt.x, nt.y);
							nd->setDelta( nt.delta * parser.scale );
							bgnodes.push_back( nd );
						}

						bgmesh->initialize(bgnodes, dummy);
					}
					
					if( ly->type == VORONOI_VERTEX )
						layer = new VoronoiVertex( bt->tag, bgmesh );
					else if( ly->type == VORONOI_SEGMENT )
						layer = new VoronoiSegment( bt->tag, bgmesh );
					else if( ly->type == SSMF_VORONOI_SEGMENT )
						layer = new SSMFVoronoiSegment( bt->tag, bgmesh );
					else if( ly->type == SSSF_VORONOI_SEGMENT )
						layer = new SSSFVoronoiSegment( bt->tag, bgmesh );
					else
						blm_error("Unknown layer type code!!!!");
				}

				// Insert the fixed nodes to the layer
				for (int i = 0; i < ly->fixedNodes.size(); i++)
				{
					NodeToken *nt = ly->fixedNodes[i];
					GeometryNode *nd = new GeometryNode( nt->tag, nt->x, nt->y);
					
					if (nt->delta <= 0.0)
						nt->delta = ly->delta > 0.0 ? ly->delta : parser.delta;
					
					std::cout << "Node " << nt->tag << ": " << nt->delta * parser.scale << std::endl;
					nd->setDelta( nt->delta * parser.scale );
					
					layer->addFixedNode( nd );
				}
			}

			int seedDirection = 1;

			Border *bor = new Border( ly->loops.size() );
			int loopLen = ly->loops.size();
			for( int j = 0; j < loopLen; ++j )
			{
				int edgeLen = ly->loops[j]->edges.size();

				for( int k = 0; k < edgeLen; k++ )
				{
					int edgeTag;
					if (ly->loops[j]->direction > 0)
						edgeTag = ly->loops[j]->edges[k];
					else
						edgeTag = -ly->loops[j]->edges[edgeLen-1-k];

					int dir = edgeTag > 0 ? 1 : -1;
					edgeTag = abs( edgeTag );

					GeometryEdgeMap::iterator it = geometryEdges.find(edgeTag);
					if (it == geometryEdges.end())
						blm_error("Nonexisting edge in layer", ly->tag);
					GeometryEdge *ed = it->second;

					if (ly->gridh > 0 && ly->gridv > 0)
						ed->setSegments(k & 1 ? ly->gridv : ly->gridh);
					
					std::map<int,int>::iterator mi;
					if ((mi = edgeBody.find(edgeTag)) == edgeBody.end())
						edgeBody[edgeTag] = bt->tag;
					else if (mi->second == bt->tag)
						ed->makeVirtual();
					else
						ed->makeInner();
					
					// Add background mesh to edge
					if (bgmesh != NULL)
						ed->addBGMesh(bgmesh);

					// Push layer to node layer list, last node from next edge
					std::vector<GeometryNode*> nds;
					ed->exportGeometryNodes(nds, dir);
					for (int n = 0; n < nds.size() - 1; n++)
						nodeLayers[nds[n]->tag].push_back(ly);
					
					bor->addLoopEdge( j, dir, ed );

					if (ly->seed && ly->seed->edge == edgeTag)
						seedDirection = dir;
				}
			}
			
			// Set the seed
			if (ly->type == SSMF_VORONOI_SEGMENT || ly->type == SSSF_VORONOI_SEGMENT)
			{
				SSVoronoiSegment *vs = static_cast<SSVoronoiSegment *>( layer );
				if( ly->seed->type == "Explicit" )
				{
					vs->makeExplicitSeed( ly->seed->nodes[0], 
						ly->seed->nodes[1], ly->seed->nodes[2] );
				}
				else if( ly->seed->type == "Implicit" )
				{
					vs->makeImplicitSeed( geometryEdges[ ly->seed->edge ], seedDirection );
				}
				else
					blm_error("Unknown seed type:", ly->seed->type.c_str());
			}
			
			layer->setBounds( bor );
			body->addLayer( layer );
			layers.push_back(layer);
		}
		
		bodies[ body->tag ] = body;
	}
	
	len = parser.nodes.size();
	for (i = 0; i < len; i++)
	{
		NodeToken *nt = parser.nodes[i];
		
		if (nt->delta <= 0.0)
		{
			double mean = 1.0;
			int n = 0;
			
			std::list<EdgeToken*> &edges = nodeEdges[nt->tag];
			std::list<EdgeToken*>::iterator it;
			for (it = edges.begin(); it != edges.end(); it++)
			{
				if ((*it)->delta > 0.0)
				{
					mean *= (*it)->delta;
					n++;
				}
			}
			
			if (n == 0)
			{
				std::list<LayerToken*> &layers = nodeLayers[nt->tag];
				std::list<LayerToken*>::iterator it;
				for (it = layers.begin(); it != layers.end(); it++)
				{
					if ((*it)->delta > 0.0)
					{
						mean *= (*it)->delta;
						n++;
					}
				}
			}
			
			if (n > 0)
				nt->delta = pow(mean, 1.0 / n);
			else
				nt->delta = parser.delta;
		}
		
		std::cout << "Node " << nt->tag << ": " << nt->delta * parser.scale << std::endl;
		geometryNodes[nt->tag]->setDelta(nt->delta * parser.scale);
	}
	
	for (std::vector<Layer*>::iterator it = layers.begin(); it != layers.end(); it++)
	{
		(*it)->initialize();
	}

	std::cout << "Conversion completed." << std::endl;
}

void Mesh::
discretize()
{
	GeometryNodeMapIt nit;
	for( nit = geometryNodes.begin(); nit != geometryNodes.end(); ++nit )
	{
		MeshNode *mnd = new MeshNode( *((*nit).second) );
		fixedNodes[ mnd->tag ] = mnd;
		mnd->boundarynode = true;
	}
	
	GeometryEdgeMapIt eit;
	for( eit = geometryEdges.begin(); eit != geometryEdges.end(); ++eit )
	{
		GeometryEdge *ed = (*eit).second;
		if (ed->isConstant()) ed->discretize( fixedNodes );
	}
	
	for( eit = geometryEdges.begin(); eit != geometryEdges.end(); ++eit )
	{
		GeometryEdge *ed = (*eit).second;
		if (!ed->isConstant()) ed->discretize( fixedNodes );
	}

	BodyMapIt bit;
	for( bit = bodies.begin(); bit != bodies.end(); ++bit )
	{
		Body *bd = (*bit).second;
		bd->discretize( fixedNodes, meshNodes, elements );
		std::cout << "Body " << bd->tag << " completed!" << std::endl;
	}
	
	for( eit = geometryEdges.begin(); eit != geometryEdges.end(); ++eit )
	{
		GeometryEdge *ed = (*eit).second;
		ed->elements(boundaryElements, 0);
	}
	
	createMiddleNodes();

	std::cout << "Nodes: " << meshNodes.size() << std::endl;
	std::cout << "Elements: " << elements.size() << std::endl;
}

void Mesh::
createMiddleNodes()
{
	std::map<std::pair<int, int>, Node *> edgeNodeMap;

	std::list<Element *>::iterator i;
	for (i = elements.begin(); i != elements.end(); i++)
	{
		Element *e = *i;
		if (bodies[e->partOfBody()]->isParabolic())
		{
			int size = e->size();
			std::vector<Node *> newNodes;
			for (int i = 0; i < size; i++)
			{
				Node *node;

				Node *node0 = e->nodeAt(i), *node1 = e->nodeAt((i+1)%size);
				std::pair<int, int> p(std::min(node0->tag, node1->tag), std::max(node0->tag, node1->tag));
				std::map<std::pair<int, int>, Node *>::iterator ni;
				if ((ni = edgeNodeMap.find(p)) == edgeNodeMap.end())
				{
					double x = (node0->x + node1->x) / 2;
					double y = (node0->y + node1->y) / 2;
					node = new MeshNode(x, y);
					edgeNodeMap[p] = node;
					meshNodes[node->tag] = node;
				}
				else
				{
					node = ni->second;
				}

				newNodes.push_back(node);
			}

#if 0
			if (e->elementType() > 400)
			{
				double x = (e->nodeAt(0)->x + e->nodeAt(1)->x + e->nodeAt(2)->x + e->nodeAt(3)->x) / 4.0;
				double y = (e->nodeAt(0)->y + e->nodeAt(1)->y + e->nodeAt(2)->y + e->nodeAt(3)->y) / 4.0;
				Node *node = new MeshNode(x, y);
				meshNodes[node->tag] = node;
				newNodes.push_back(node);
			}
#endif

			e->upgrade(newNodes);
		}
	}

	std::vector<BoundaryElement *>::iterator bi;
	for (bi = boundaryElements.begin(); bi != boundaryElements.end(); bi++)
	{
		BoundaryElement *e = *bi;
		int t0 = e->from()->tag, t1 = e->to()->tag;
		std::pair<int, int> p(std::min(t0, t1), std::max(t0, t1));
		std::map<std::pair<int, int>, Node *>::iterator ni;
		if ((ni = edgeNodeMap.find(p)) != edgeNodeMap.end())
			e->addMiddleNode(ni->second);
	}
}

void	Mesh::output(std::ofstream& o)
{
	NodeMapIt noit;
	std::list< Element * >::iterator elit;
	GeometryEdgeMapIt git;
	
	int index = 1;
	
	o << meshNodes.size() << ' ' << elements.size() << std::endl;
	for( noit = meshNodes.begin(); noit != meshNodes.end(); ++noit)
	{
		Node *nd = (*noit).second; 
		nd->tag = index++;
		o << *nd;
	}
	
	for( elit = elements.begin(); elit != elements.end(); ++elit)
	{
		o << **elit;
	}
}

void Mesh::outputHeader(std::ofstream& o)
{
	int ngnds = 0;
	GeometryNodeMap::iterator it;
	for (it = geometryNodes.begin(); it != geometryNodes.end(); it++)
		if (it->second->boundaryTag > 0)
			ngnds++;

	o << meshNodes.size() << ' ' << elements.size() << ' ' << boundaryElements.size() + ngnds << '\n';
	
	std::map<int, int> n;

	std::list< Element* >::const_iterator e;
	for (e = elements.begin(); e != elements.end(); e++)
		n[(*e)->elementType()]++;
	
	std::vector< BoundaryElement* >::const_iterator be;
	for (be = boundaryElements.begin(); be != boundaryElements.end(); be++)
		n[(*be)->middle() != NULL ? 203 : 202]++;
	
	n[101] += ngnds;

	o << n.size() << '\n';

	std::map<int, int>::iterator i;
	for (i = n.begin(); i != n.end(); i++)
		o << i->first << ' ' << i->second << '\n';
}

void Mesh::outputNodes(std::ofstream &o)
{
	int tag = 1;
	NodeMapIt noit;

//      o.setf(ios::scientific);
        o.precision(16);

	for( noit = meshNodes.begin(); noit != meshNodes.end(); ++noit)
	{
		Node *nd = (*noit).second; 
		nd->tag = tag++;
		o << *nd;
	}
}

void Mesh::outputElements(std::ofstream &o)
{
	std::list< Element * >::iterator elit;
	
	for( elit = elements.begin(); elit != elements.end(); ++elit)
	{
		o << **elit;
	}
}

void Mesh::outputBoundary(std::ofstream &o)
{
	std::vector< BoundaryElement* >::iterator bel;
	for (bel = boundaryElements.begin(); bel != boundaryElements.end(); bel++)
	{
		o << **bel;
	}
	
	int id = boundaryElements.size() + 1;
	GeometryNodeMapIt n;
	for (n = geometryNodes.begin(); n != geometryNodes.end(); n++)
	{
		if (n->second->boundaryTag > 0)
			o << id++ << ' ' << n->second->boundaryTag << " -1 -1 101 " << n->first << '\n';
	}
}
