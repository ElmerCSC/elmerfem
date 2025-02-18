#include "MeshParser.h"
#include "MGError.h"
#include <string>
#include <iostream>
#include <ctype.h>

MeshParser::MeshParser( const char *path )
{
	s.open( path );
	if (s.fail())
		blm_error( "Error opening file:", path );

	bgmeshes = NULL;
}

MeshParser::MeshParser( const char *path, std::map<int, BGMeshToken *> *ebgs )
{
	s.open( path );
	if (s.fail())
		blm_error( "Error opening file:", path );

	bgmeshes = ebgs;
}

MeshParser::~MeshParser()
{
	s.close();
}

void MeshParser::readInputFile()
{
	if( !readHeader() )
	{
		blm_error( "Error reading header." );
	}
	
	if( !readNodes() )
	{
		blm_error( "Error reading nodes." );
	}
	
	if( !readEdges() )
	{
		blm_error( "Error reading edges." );
	}
	
	if( !readBodies() )
	{
		blm_error( "Error reading bodies." );
	}

	std::cout << "Parse OK" << std::endl;
}

int MeshParser::
readHeader()
{
	if( !readWord( "Geometry2D:" ) )  return 0;
	if( !readWord( "H:" ) ) return 0;
	if( !readDouble( delta ) ) return 0;
	if( !readWord( "MeshScalingFactor:" ) ) return 0;
	if( !readDouble( scale ) ) return 0;
	if( !readWord( "Nodes:" ) )  return 0;
	if( !readInteger( nodeCount ) ) return 0;
	if( !readWord( "Edges:" ) )  return 0;
	if( !readInteger( edgeCount ) ) return 0;
	if( !readWord( "Bodies:" ) )	return 0;
	if( !readInteger( bodyCount ) ) return 0;
	
	nodes.reserve( nodeCount );
	edges.reserve( edgeCount );
	bodies.reserve( bodyCount );
	
	return 1;
}

int MeshParser::
readNodes()
{
	int i;
	for( i = 0; i < nodeCount; ++i )
	{
		NodeToken *nd = new NodeToken;
		if( !readWord( "NodeId:" ) )	return 0;
		if( !readInteger( nd->tag ) ) return 0;
		if( !readInteger( nd->boundaryTag ) ) return 0;
		
		while (ignoreComment());
		char ch = s.peek();
		if (ch == 'H')
		{
			if (!readWord( "H:" ) ) return 0;
			if (!readDouble( nd->delta ) ) return 0;
		}
		else if (ch == 'R')
		{
			if (!readWord( "R:" ) ) return 0;
			if (!readDouble( nd->delta ) ) return 0;
			nd->delta *= delta;
		}
		else
		{
			nd->delta = -1.0;
		}
		
		if( !readDouble( nd->x ) ) return 0;
		if( !readDouble( nd->y ) ) return 0;
		nodes.push_back( nd );		
	}
	return 1;
}

int MeshParser::
readEdges()
{
	int i;
	int edt;
	int len;
	
	for( i = 0; i < edgeCount; ++i )
	{
		EdgeToken *ed = new EdgeToken;
		if( !readWord( "EdgeId:" ) ) return 0;
		if( !readInteger( ed->tag ) ) return 0;
		if( !readInteger( ed->boundaryTag ) ) return 0;
		
		while (ignoreComment());
		char ch = s.peek();
		if (ch == 'H')
		{
			if (!readWord("H:")) return 0;
			if (!readDouble(ed->delta)) return 0;
			ed->segCount = 0;
		}
		else if (ch == 'R')
		{
			if (!readWord("R:")) return 0;
			if (!readDouble(ed->delta)) return 0;
			ed->delta *= delta;
			ed->segCount = 0;
		}
		else if (ch == 'N')
		{
			if (!readWord("N:")) return 0;
			if (!readInteger(ed->segCount)) return 0;
			ed->delta = -1.0;
		}
		else
		{
			ed->delta = -1.0;
			ed->segCount = 0;
		}
		
		if( !readInteger( len ) ) return 0;
		for(int j = 0; j < len; ++j)
		{
			if( !readInteger( edt ) ) return 0;
			ed->nodes.push_back( edt );
		}
		
		edges.push_back( ed );
	}
	return 1;
}

int MeshParser::
readBodies()
{
	int i;
	for( i = 0; i < bodyCount; ++i )
	{
		int layerCount;
		currentBody = new BodyToken;
		if( !readWord( "BodyId:" ) )	return 0;
		if( !readInteger( currentBody->tag ) ) return 0;
		
		while (ignoreComment());
		char ch = s.peek();
		if (ch == 'H')
		{
			if (!readWord( "H:" ) ) return 0;
			if (!readDouble( currentBody->delta ) ) return 0;
		}
		else if (ch == 'R')
		{
			if (!readWord( "R:" ) ) return 0;
			if (!readDouble( currentBody->delta ) ) return 0;
			currentBody->delta *= delta;
		}
		else
		{
			currentBody->delta = -1.0;
		}
		
		if( !readWord( "ElementOrder:" ) ) return 0;
		std::string order;
		if( !readString( order ) ) return 0;
		if (order == "Parabolic")
			currentBody->parabolic = true;
		else if (order == "Linear")
			currentBody->parabolic = false;
		else
			return 0;
		
		if( !readWord( "Layers:" ) )	return 0;
		if( !readInteger( layerCount ) ) return 0;
		if( !readLayers( layerCount ) ) return 0;
		bodies.push_back( currentBody );
	}
	return 1;
}

int MeshParser::
readLayers( const int layerCount )
{
	int i;
	for( i = 0; i < layerCount; ++i )
	{
		int fixedNodeCount;
		int loopCount;
		currentLayer = new LayerToken;
		if( !readWord( "LayerId:" ) )  return 0;
		if( !readInteger( currentLayer->tag ) ) return 0;
		
		while (ignoreComment());
		char ch = s.peek();
		if (ch == 'H')
		{
			if (!readWord( "H:" ) ) return 0;
			if (!readDouble( currentLayer->delta ) ) return 0;
		}
		else if (ch == 'R')
		{
			if (!readWord( "R:" ) ) return 0;
			if (!readDouble( currentLayer->delta ) ) return 0;
			currentLayer->delta *= delta;
		}
		else
		{
			currentLayer->delta = -1.0;
		}
		
		if( !readWord( "LayerType:" ) )	return 0;
		std::string type;
		if( !readString( type ) )	return 0;
		if (type == "Connect")
			currentLayer->type = CONNECT;
		else if (type == "BoundaryMesh")
			currentLayer->type = BOUNDARY_MESH;
		else if (type == "VoronoiVertex")
			currentLayer->type = VORONOI_VERTEX;
		else if (type == "MovingFront")
			currentLayer->type = VORONOI_SEGMENT;
		else if (type == "SSSFMovingFront")
			currentLayer->type = SSSF_VORONOI_SEGMENT;
		else if (type == "SSMFMovingFront")
			currentLayer->type = SSMF_VORONOI_SEGMENT;
		else if (type == "QuadGrid")
			currentLayer->type = QUAD_GRID;
		else if (type == "TriangleNEGrid")
			currentLayer->type = TRIANGLE_NE_GRID;
		else if (type == "TriangleNWGrid")
			currentLayer->type = TRIANGLE_NW_GRID;
		else if (type == "TriangleUJNEGrid")
			currentLayer->type = TRIANGLE_UJ_NE_GRID;
		else if (type == "TriangleUJNWGrid")
			currentLayer->type = TRIANGLE_UJ_NW_GRID;
		else if (type == "TriangleFBNEGrid")
			currentLayer->type = TRIANGLE_FB_NE_GRID;
		else if (type == "TriangleFBNWGrid")
			currentLayer->type = TRIANGLE_FB_NW_GRID;
		else
			return 0;
		
		if( currentLayer->type == QUAD_GRID ||
			currentLayer->type == TRIANGLE_NE_GRID ||
			currentLayer->type == TRIANGLE_NW_GRID ||
			currentLayer->type == TRIANGLE_UJ_NE_GRID ||
			currentLayer->type == TRIANGLE_UJ_NW_GRID ||
			currentLayer->type == TRIANGLE_FB_NE_GRID ||
			currentLayer->type == TRIANGLE_FB_NW_GRID)
		{
			std::string keyword;
			if (!readString(keyword)) return 0;
			if (keyword == "GridSize:")
			{
				if (!readInteger(currentLayer->gridh)) return 0;
				if (!readInteger(currentLayer->gridv)) return 0;
				if (!readString(keyword)) return 0;
			}
			
			if (keyword != "Loops:")
				return 0;

			if (!readInteger(loopCount)) return 0;
			if (loopCount != 1) return 0;
			if (!readLoops(loopCount)) return 0;
		}
		else
		{
			if( currentLayer->type == VORONOI_VERTEX ||
				currentLayer->type == VORONOI_SEGMENT ||
				currentLayer->type == SSMF_VORONOI_SEGMENT ||
				currentLayer->type == SSSF_VORONOI_SEGMENT)
			{
				if( !readWord( "FixedNodes:" ) )  return 0;
				if( !readInteger( fixedNodeCount ) ) return 0;
				if( !readFixedNodes( fixedNodeCount ) ) return 0;
				if( !readWord( "BGMesh:" ) )	return 0;
				if( !readBGMesh() ) return 0;
			}
			
			if( currentLayer->type == SSMF_VORONOI_SEGMENT ||
				currentLayer->type == SSSF_VORONOI_SEGMENT)
			{
				if( !readWord( "Seed:" ) )  return 0;
				currentSeed = new SeedToken;
				if( !readSeed() ) return 0;
				currentLayer->seed = currentSeed;
			}
			
			if( !readWord( "Loops:" ) )  return 0;
			if( !readInteger( loopCount ) ) return 0;
			if( !readLoops( loopCount ) ) return 0;
		}

		currentBody->layers.push_back( currentLayer );
	}
	return 1;
}

int MeshParser::
readSeed()
{
	if( !readString( currentSeed->type ) )  return 0;
	if( currentSeed->type == "Explicit" )
	{
		int tag;
		if( !readWord( "Nodes:" ) )  return 0;
		if( !readInteger( tag ) ) return 0;
		currentSeed->nodes.push_back( tag );
		if( !readInteger( tag ) ) return 0;
		currentSeed->nodes.push_back( tag );
		if( !readInteger( tag ) ) return 0;
		currentSeed->nodes.push_back( tag );
	}
	else if( currentSeed->type == "Implicit" )
	{
		if( !readWord( "Edge:" ) )  return 0;
		if( !readInteger( currentSeed->edge ) ) return 0;
	}
	else return 0;
	return 1;
}

int MeshParser::
readLoops( const int loopCount )
{
	int i;
	for( i = 0; i < loopCount; ++i )
	{
		int edCount;
		LoopToken *loop = new LoopToken;
		if( !readWord( "LoopId:" ) )	return 0;
		if( !readInteger( loop->tag ) ) return 0;
		if( !readWord( "Direction:" ) ) return 0;
		if( !readInteger( loop->direction ) ||
			loop->direction != 1 && loop->direction != -1 ) return 0;
		if( !readWord( "Edges:" ) )  return 0;
		if( !readInteger( edCount ) ) return 0;
		for( int j = 0; j < edCount; ++j )
		{
			int ed;
			if( !readInteger( ed ) ) return 0;
			loop->edges.push_back( ed );
		}
		currentLayer->loops.push_back( loop );
	}
	return 1;
}

int MeshParser::
readFixedNodes( const int fixedNodeCount )
{
	int i;
	for( i = 0; i < fixedNodeCount; ++i )
	{
		NodeToken *nd = new NodeToken;
		if( !readWord( "NodeId:" ) )	return 0;
		if( !readInteger( nd->tag ) ) return 0;
		if( !readInteger( nd->boundaryTag ) ) return 0;
		
		while (ignoreComment());
		char ch = s.peek();
		if (ch == 'H')
		{
			if (!readWord( "H:" ) ) return 0;
			if (!readDouble( nd->delta ) ) return 0;
		}
		else if (ch == 'R')
		{
			if (!readWord( "R:" ) ) return 0;
			if (!readDouble( nd->delta ) ) return 0;
			nd->delta *= delta;
		}
		else
		{
			nd->delta = -1.0;
		}
		
		if( !readDouble( nd->x ) ) return 0;
		if( !readDouble( nd->y ) ) return 0;
		
		currentLayer->fixedNodes.push_back( nd ); 	
	}
	return 1;
}

int MeshParser::
readBGMesh()
{
	currentLayer->bg = NULL;
	if (bgmeshes != NULL)
	{
		std::map<int, BGMeshToken*>::iterator mesh = bgmeshes->find(currentLayer->tag);
		if (mesh != bgmeshes->end())
			currentLayer->bg = mesh->second;
		else if ((mesh = bgmeshes->find(-1)) != bgmeshes->end())
			currentLayer->bg = mesh->second;
	}

	BGMeshToken *bg = new BGMeshToken;
	if( !readString( bg->type ) ) return 0;
	if( bg->type == "Explicit" )
	{
		if (!readExplicitBGMesh(bg))
			return 0;
	}
	else if ( bg->type == "External" )
	{
		std::string filename;
		if (!readString( filename ) ) return 0;

		if (currentLayer->bg == NULL)
		{
			MeshParser bgmeshparser(filename.c_str());
			if (!bgmeshparser.readExplicitBGMesh(bg))
				return 0;
		}
		bg->type = "Explicit";
	}

	if (currentLayer->bg == NULL)
		currentLayer->bg = bg;
	else
		delete bg;

	return 1;
}

int MeshParser::
readExplicitBGMesh(BGMeshToken *bg)
{
	int len;
	if( !readInteger( len ) ) return 0;
	
	for(int i = 0; i < len; ++i)
	{
		NodeToken nd;
		if( !readDouble( nd.x ) ) return 0;
		if( !readDouble( nd.y ) ) return 0;
		if( !readDouble( nd.delta ) ) return 0;
		bg->nodes.push_back( nd );
	}

        return 1;
}

int MeshParser::
readBoundary()
{
	int i, len, ed;
	boundary = new BoundaryToken;
	if( !readWord( "Boundaries:" ) ) return 0;
	if( !readWord( "OuterBoundaries:" ) ) return 0;
	if( !readInteger( len ) ) return 0;
	for( i = 0; i < len; ++i )
	{
		if( !readInteger( ed ) ) return 0;
		boundary->outerBoundaries.push_back( ed );
	}
	if( !readWord( "InnerBoundaries:" ) ) return 0;
	if( !readInteger( len ) ) return 0;
	for( i = 0; i < len; ++i )
	{
		if( !readInteger( ed ) ) return 0;
		boundary->innerBoundaries.push_back( ed );
	}
	return 1;
}

// Read a whitespace-delimited word, with possible quoting for whitespace and comment characters
int MeshParser::
readString( std::string& ref )
{
	ref.erase();

	bool head = true, inquote = false, incomment = false;

	while (!s.eof())
	{
		int ch = s.get();

		if (incomment)
		{
			if (ch != '\n')
				continue;
			else
				incomment = false;
		}
		
		// We are not in a comment
		if (ch == '\"')
		{
			inquote = !inquote;
		}
		else if (!inquote && isspace(ch))
		{
			if (head)
				continue;
			else
				break;
		}
		else if (!inquote && (ch == '!' || ch == '#'))
		{
			incomment = true;
		}
		else
		{
			head = false;
			ref.insert(ref.end(), static_cast<char>(ch));
		}
	}

	return !s.fail();
}

int MeshParser::
readInteger( int& ref )
{
	while (ignoreComment());
	s >> ref;
	return !s.fail();
}

int MeshParser::
readDouble( double& ref )
{
	while (ignoreComment());
	s >> ref;
	return !s.fail();
}

int MeshParser::
readWord( const char *ref )
{
	while (ignoreComment());
	std::string word; 
	s >> word;
	if( word != ref ) return 0;
	return !s.fail();
}

int MeshParser::
ignoreComment()
{
	int ch = 0;
	while (!s.eof() && isspace(ch = s.peek()))
		s.ignore(1);

	if( ch == '!' || ch == '#' )
	{
		s.ignore(0x7fffffff, '\n');
		return 1;
	}
	return 0;
}
