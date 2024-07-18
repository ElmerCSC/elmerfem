/*
*/

// Stdlib
#include "../config.h"
#include <stdlib.h>
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <stdio.h>

#include <string.h>

#include <iostream>
#include <fstream>

/*
*	BoundaryLayerMesh specific includes
*/
#include "MeshParser.h"
#include "Mesh.h"
#include "MGError.h"

/*
*	Simple error&exit routine.
*/

template<class T>
void blm_error(const char* msg, const T &opt);

int main(int argc, char **argv)
{
	const char *usage = "Usage: Mesh2D [--bgmesh=filename] [--bgcontrol=filename] <input file> [mesh directory]";

	// Read command line parameters
	char *modelfile = NULL;
	char *meshdir = NULL;
	std::map<int, BGMeshToken*> externalBGMeshes;

	for (int i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "--", 2) != 0)
		{
			if (modelfile == NULL) // First parameter without "--" is the input file
				modelfile = argv[i];
			else if (meshdir == NULL) // Second parameter without "--" is the output directory
				meshdir = argv[i];
			else
			{
				std::cerr << usage << std::endl;
				return 1;
			}
		}
		else if (strncmp(argv[i], "--bgcontrol=", 12) == 0)
		{
			std::ifstream file(argv[i] + 12);
			while (!file.fail() && !file.eof())
			{
				int id;
				std::string filename;
				file >> id;
				file >> filename;
				
				if (!file.fail())
				{
					MeshParser bgmeshparser(filename.c_str());
					BGMeshToken *token = new BGMeshToken;
					token->type = "Explicit";
					if (!bgmeshparser.readExplicitBGMesh(token))
						blm_error("Could not read background mesh from file:", filename);
					externalBGMeshes[id] = token;
				}
			}
		}
		else if (strncmp(argv[i], "--bgmesh=", 9) == 0)
		{
			BGMeshToken *token = new BGMeshToken;
			token->type = "Explicit";

			MeshParser bgmeshparser(argv[i] + 9);
			if (!bgmeshparser.readExplicitBGMesh(token))
				blm_error("Could not read background mesh from file:", argv[i] + 9);

			externalBGMeshes[-1] = token;
		}
		else
		{
			std::cerr << usage << std::endl;
			return 1;
		}
	}

	if (modelfile == NULL)
	{
		std::cerr << usage << std::endl;
		return 1;
	}

	std::cout << "Mesh input file: " << modelfile << std::endl;
	std::cout << "Mesh output directory: ";
	if (meshdir != NULL)
		std::cout << meshdir;
	else
	{
		char buf[1024];
		if (getcwd(buf, 1024) != NULL)
			std::cout << buf;
	}
	std::cout << std::endl;
	
	Mesh mesh;
	MeshParser parser( modelfile, &externalBGMeshes );
	parser.readInputFile();
	
	if (meshdir && chdir(meshdir) != 0)
		blm_error("Could not access directory:", meshdir);
	
	mesh.convertEdgeFormat( parser );
	mesh.discretize();
	
	std::ofstream header("mesh.header", std::ios::out | std::ios::trunc);
	if(header.fail())
		blm_error("Could not open header file.");
	
	mesh.outputHeader( header );
	header.close();
	
	std::ofstream nodes("mesh.nodes", std::ios::out | std::ios::trunc);
	if (nodes.fail())
		blm_error("Could not open nodes file.");
	
	mesh.outputNodes( nodes );
	nodes.close();
	
	std::ofstream elements("mesh.elements", std::ios::out | std::ios::trunc);
	if (elements.fail())
		blm_error("Could not open elements file.");
	
	mesh.outputElements( elements );
	elements.close();
	
	std::ofstream boundary("mesh.boundary", std::ios::out | std::ios::trunc);
	if (boundary.fail())
		blm_error("Could not open boundary file.");
	
	mesh.outputBoundary( boundary );
	boundary.close();
	
	std::cout << "*** ALL DONE" << std::endl;
	
	return 0;
}
