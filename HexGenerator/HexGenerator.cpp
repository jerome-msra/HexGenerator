// HexGenerator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

// algorithms
#include "SurfaceCylinderMapper.h"
#include "TetCylinderHarmonicMapper.h"
#include "TetToHexGenerator.h"
#include "OptimizeHex.h"
#include "OptimizeTet.h"
#include "MergeHex.h"

// meshes
#include "GeneralQuadMesh.h"
#include "HexGenerateHMesh.h"
#include "HexGenerateTMesh.h"
#include "TetSurfaceMesh.h"
#include "TetGeneralTMesh.h"
#include "HarmonicMapTMesh.h"
#include "OptimizeHexMesh.h"
#include "OptimizeTetMesh.h"
#include "MergeHexMesh.h"

#define DEFAULT_LAYERS 40
#define DEFAULT_XQUAD 7
#define DEFAULT_YQUAD 7


using namespace MeshLib;

void print_help();
void parse_io();

void test(int argc, char * argv[])
{
	std::vector<std::string> inputNames;
	std::vector<CGQMesh*> quadMeshes;
	std::vector<HMeshLib::CHGHMesh*> hexMeshes;
	std::vector<std::vector<int>> pairList;
	for (int i = 0; i < argc; i++)
	{
		std::string para(argv[i]);
		if (para == "-input")
		{
			std::stringstream strStream;
			strStream << argv[i + 1];
			unsigned int numMeshes;
			strStream >> numMeshes;
			for (size_t h = 0; h < numMeshes; h++)
			{
				CGQMesh * mesh = new CGQMesh();
				std::cout << "Reading Quad Mesh " << argv[i + h + 2] << "..." << std::endl;
				mesh->read_qm(argv[i + h + 2]);
				quadMeshes.push_back(mesh);
			}
		}
		else if (para == "-hex")
		{
			std::stringstream strStream;
			strStream << argv[i + 1];
			unsigned int numMeshes;
			strStream >> numMeshes;
			for (size_t h = 0; h < numMeshes; h++)
			{
				HMeshLib::CHGHMesh * mesh = new HMeshLib::CHGHMesh();
				inputNames.push_back(argv[i + h + 2]);
				std::cout << "Reading Quad Mesh " << argv[i + h + 2] << "..." << std::endl;
				mesh->_load_hm(argv[i + h + 2]);
				hexMeshes.push_back(mesh);
			}
		}
		else if (para == "-pair")
		{
			std::ifstream input;
			input.open(argv[i + 1]);
			int first, second;
			while (input.good())
			{
				input >> first >> second;
				std::vector<int> tmp;
				tmp.push_back(first);
				tmp.push_back(second);
				pairList.push_back(tmp);
			}
		}
	}

	for (size_t i = 0; i < quadMeshes.size(); i++)
	{
		std::map<size_t, size_t> mymap;
		CGQMesh * mesh = quadMeshes[i];
		int v = 0;
		for (CGQMesh::MeshVertexIterator vIter(mesh); !vIter.end(); vIter++, v++)
		{
			CGQMesh::CVertex * pV = *vIter;
			mymap[v] = pV->id();
		}

		for (size_t p = 0; p < pairList.size(); p++)
		{
			pairList[p][i] = (int)mymap[pairList[p][i]];
		}
	}

	for (size_t p = 0; p < pairList.size(); p++)
	{

		CGQMesh::CVertex * pV1 = quadMeshes[0]->idVertex(pairList[p][0]);
		CGQMesh::CVertex * pV2 = quadMeshes[1]->idVertex(pairList[p][1]);

		HMeshLib::CHGHMesh::CVertex * pHV1 = hexMeshes[0]->findVertex(pV1->point());
		HMeshLib::CHGHMesh::CVertex * pHV2 = hexMeshes[1]->findVertex(pV2->point());

		CPoint p1 = pHV1->position();
		CPoint p2 = pHV2->position();

		pHV1->position() = (p1 + p2) / 2.0;
		pHV2->position() = (p1 + p2) / 2.0;
	}

	for (size_t t = 0; t < inputNames.size(); t++)
	{
		size_t idx = inputNames[t].find_last_of(".");
		inputNames[t] = inputNames[t].substr(0, idx);
		inputNames[t] += "_merge.hm";
	}

	hexMeshes[0]->_write_hm(inputNames[0].c_str());
	hexMeshes[1]->_write_hm(inputNames[1].c_str());
};

void parse_io(int argc, char * argv[], std::string & inputName, std::string & outputName)
{

	for (int c = 1; c < argc; c++)
	{
		if (strcmp(argv[c], "-input") == 0)
		{
			inputName = argv[c + 1];
		}
		else if (strcmp(argv[c], "-output") == 0)
		{
			outputName = argv[c + 1];
		}
	}

	if (inputName.size() == 0 || outputName.size() == 0)
	{
		std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
		print_help();
		exit(0);
	}

	return;
};

void generate_hex(int argc, char * argv[])
{
	std::string inputName, outputName;
	parse_io(argc, argv, inputName, outputName);

	int numLayers = 0;
	int xQuad = 0;
	int yQuad = 0;
	bool zero = false;

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-layer") == 0)
		{
			numLayers = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-xquad") == 0)
		{
			xQuad = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-yquad") == 0)
		{
			yQuad = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-zero") == 0)
		{
			zero = true;
		}
	}

	std::cout << "[GenHex]Reading Tetrahedron Mesh " << inputName << "..." << std::endl;
	TMeshLib::CHGTMesh * tetMesh = new TMeshLib::CHGTMesh();
	tetMesh->_load_t(inputName.c_str());

	numLayers = (numLayers == 0) ? DEFAULT_LAYERS : numLayers;
	xQuad = (xQuad == 0) ? DEFAULT_XQUAD : xQuad;
	yQuad = (yQuad == 0) ? DEFAULT_YQUAD : yQuad;

	std::string outputzero = (zero) ? "true" : "false";

	std::cout << std::endl;
	std::cout << "==============================================================" << std::endl;
	std::cout << "[GenHex] # of X quad: " << xQuad << std::endl;
	std::cout << "[GenHex] # of Y quad: " << yQuad << std::endl;
	std::cout << "[GenHex] # of layers: " << numLayers << std::endl;
	std::cout << "[GenHex] Rotate to Zero: " << outputzero << std::endl;
	std::cout << "==============================================================" << std::endl;

	HMeshLib::CTetHexGenerator<TMeshLib::CHGTMesh, HMeshLib::CHGHMesh> * hexGenerator = new HMeshLib::CTetHexGenerator<TMeshLib::CHGTMesh, HMeshLib::CHGHMesh>(tetMesh);
	
	// set properties
	hexGenerator->_zero() = zero;
	hexGenerator->_numLayer() = numLayers;
	hexGenerator->_xQuad() = xQuad;
	hexGenerator->_yQuad() = yQuad;

	hexGenerator->_genHex();

	hexGenerator->_writeHex(outputName.c_str());
};

void merge_hex(int argc, char * argv[])
{
	std::string outputName;

	std::vector<HMeshLib::CMHMesh *> meshes;

	HMeshLib::MergeHex<HMeshLib::CMHMesh> * hexmerger = new HMeshLib::MergeHex<HMeshLib::CMHMesh>();

	for (int i = 0; i < argc; i++)
	{
		std::string para(argv[i]);
		if (para == "-input")
		{
			std::stringstream strStream;
			strStream << argv[i + 1];
			unsigned int numMeshes;
			strStream >> numMeshes;
			for (size_t h = 0; h < numMeshes; h++)
			{
				HMeshLib::CMHMesh * mesh = new HMeshLib::CMHMesh();
				std::cout << "Reading Hex Mesh " << argv[i + h + 2] << "..." << std::endl;
				mesh->_load_hm(argv[i + h + 2]);
				mesh->labelBoundaryVertices();
				mesh->labelBaseZeros();
				mesh->_write_hm(("mergeHexTest_" + std::to_string(h) + ".hm").c_str());
				std::cout << "loaded" << std::endl;
				hexmerger->_addMesh(mesh);
			}
		}
		else if (para == "-output")
		{
			outputName = argv[i + 1];
		}
	}

	hexmerger->_merge();

	return;
};

//-mapcylinder -input data\c2_pert\c2_0_pert_group.t -surface data\c2_pert\c2_0_pert.cylinder.m -output data\c2_pert\c2_0_pert.cylinder.t
void map_cylinder(int argc, char * argv[])
{
	std::string inputName, outputName;
	parse_io(argc, argv, inputName, outputName);

	std::cout << "[map_cylinder]Reading Tetrahedron Mesh " << inputName << "..." << std::endl;
	TMeshLib::CHTMesh * tetMesh = new TMeshLib::CHTMesh();
	tetMesh->_load_t(inputName.c_str());

	CTSMesh * surface = new CTSMesh();
	std::string surfaceName;
	for (int i = 0; i < argc - 1; i++)
	{
		if (strcmp(argv[i], "-surface") == 0)
		{
			surfaceName = argv[i + 1];
		}
	}

	if (surfaceName.length() == 0)
	{
		std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
		print_help();
		return;
	}

	std::cout << "[map_cylinder]Reading Surface Mesh " << surfaceName << "..." << std::endl;
	surface->read_m(surfaceName.c_str());

	TMeshLib::CTetCylinderHarmonicMapper<CTSMesh, TMeshLib::CHTMesh> * mapper = 
		new TMeshLib::CTetCylinderHarmonicMapper<CTSMesh, TMeshLib::CHTMesh>(surface, tetMesh);

	std::cout << "[map_cylinder]Mapping Tetrahedron Mesh ..." << std::endl;
	//mapper->_iterativeMap(1, 1e-6);
	mapper->_map();

	std::cout << "[map_cylinder]Writing Tetrahedron Mesh " << outputName << " ..." << std::endl;
	mapper->_write_t(outputName.c_str());

};
//-mapsurface -input data\c2_pert\c2_0_pert_group.m -tet data\c2_pert\c2_0_pert_group.t -domain data\c2_pert\c2_0_pert.uv.m -height 2 -output data\c2_pert\c2_0_pert.cylinder.m
void map_surface_cylinder(int argc, char * argv[])
{
	std::string inputName, outputName;
	parse_io(argc, argv, inputName, outputName);

	std::cout << "[map_surface]Reading Tetrahedron Surface Mesh " << inputName << "..." << std::endl;

	CTSMesh * surfaceMesh = new CTSMesh();
	surfaceMesh->read_m(inputName.c_str());

	CTSMesh * domainMesh = new CTSMesh();

	TMeshLib::CHTMesh * harmonicTetMesh = new TMeshLib::CHTMesh();
	std::string tetName;
	std::string domainName;
	double height = 1.0;
	MAP_CYLINDER cylinderOption = MAP_CYLINDER::MAP_V;
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-tet") == 0)
		{
			tetName = argv[i + 1];
		}
		else if (strcmp(argv[i], "-domain") == 0)
		{
			domainName = argv[i + 1];
		}
		else if (strcmp(argv[i], "-height") == 0)
		{
			height = atof(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-u") == 0)
		{
			cylinderOption = MAP_CYLINDER::MAP_U;
		}
		else if (strcmp(argv[i], "-v") == 0)
		{
			cylinderOption = MAP_CYLINDER::MAP_V;
		}
	}

	if (tetName.length() == 0 || domainName.length() == 0)
	{
		std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
		print_help();
		return;
	}
	
	std::cout << "[map_surface]Reading Tetrahedron Mesh " << tetName << "..." << std::endl;
	harmonicTetMesh->_load_t(tetName.c_str());
	std::cout << "[map_surface]Reading Domain Mesh " << domainName << "..." << std::endl;
	domainMesh->read_m(domainName.c_str());

	CSurfaceCylinderMapper<MeshLib::CTSMesh, TMeshLib::CHTMesh> * mapper = new CSurfaceCylinderMapper<MeshLib::CTSMesh, TMeshLib::CHTMesh>(surfaceMesh);

	mapper->_cylinderOption() = cylinderOption;
	std::cout << "[map_surface]Connecting With Tetrahedron Mesh..." << std::endl;
	mapper->_connectTet(harmonicTetMesh);

	std::cout << height << std::endl;
	std::cout << "[map_surface]Mapping..." << std::endl;
	mapper->_map(domainMesh, height);

	std::cout << "[map_surface]Writing Cylinder Mesh " << outputName << "..." << std::endl;
	mapper->_write_m(outputName.c_str());

	return;
};

void optimize_tet(int argc, char * argv[])
{
	std::string inputName, outputName;
	parse_io(argc, argv, inputName, outputName);

	TMeshLib::OptimizeTet<TMeshLib::COTMesh> * tetOptimizer = new TMeshLib::OptimizeTet<TMeshLib::COTMesh>();

	TMeshLib::COTMesh * tmesh = new TMeshLib::COTMesh();

	std::cout << "[optimize_tet]Reading Tetrahedron Mesh " << inputName << "..." << std::endl;
	tmesh->_load_t(inputName.c_str());

	tetOptimizer->_addTetMesh(tmesh);

	bool smoothSurface = false, smoothInterior = false;

	int iteration = 1;

	for (int i = 0; i < argc; i++)
	{
		std::string para(argv[i]);
		if (para == "-interior")
		{
			smoothInterior = true;
		}
		else if (para == "-surface")
		{
			smoothSurface = true;
		}
		else if (para == "-iteration")
		{
			iteration = atoi(argv[i + 1]);
		}
	}

	if (smoothInterior)
	{
		tetOptimizer->_smoothInterior(iteration);
	}
	else if (smoothSurface)
	{
		tetOptimizer->_smoothSurface(iteration);
	}

	tetOptimizer->_writeTet(outputName.c_str());
}

//-optihex - input 3 data\eight_40\eight_1_40.hm data\eight_40\eight_2_40.hm data\eight_40\eight_3.hm - output data\eight_40\output
void optimize_hex(int argc, char * argv[])
{

	HMeshLib::OptimizeHex<HMeshLib::COHMesh> * hexOptimizer = new HMeshLib::OptimizeHex<HMeshLib::COHMesh>();

	std::string outputName;

	bool smoothSurface = false, smoothInterior = false;

	int iteration = 1;

	for (int i = 0; i < argc; i++)
	{
		std::string para(argv[i]);
		if (para == "-input")
		{
			std::stringstream strStream;
			strStream << argv[i + 1];
			unsigned int numMeshes;
			strStream >> numMeshes;
			for (size_t h = 0; h < numMeshes; h++)
			{
				HMeshLib::COHMesh * mesh = new HMeshLib::COHMesh();
				std::cout << "Reading Hex Mesh " << argv[i + h + 2] << "..." << std::endl;
				mesh->_load_hm(argv[i + h + 2]);
				hexOptimizer->_addHexMesh(mesh);
			}
		}
		else if (para == "-output")
		{
			outputName = argv[i + 1];
		}
		else if (para == "-interior")
		{
			smoothInterior = true;
		}
		else if (para == "-surface")
		{
			smoothSurface = true;
		}
		else if (para == "-iteration")
		{
			iteration = atoi(argv[i + 1]);
		}
	}

	if (smoothSurface)
	{
		hexOptimizer->_smoothSurface(iteration);
	}
	else if (smoothInterior)
	{
		hexOptimizer->_smoothInterior(iteration);
	}

	hexOptimizer->_writeHex(outputName.c_str());

	return;
};

void print_help()
{
	std::cout << std::endl << "Usage: HexGenerator.exe -[Function Parameters] -input [Input File Name] -output [Output File Name]" << std::endl << std::endl;

	std::cout << "         -h : Show usage" << std::endl << std::endl;

	std::cout << "		   -mapcylinder	: Use harmonic map to map the tetrahedron mesh to a standard cylinder" << std::endl;
	std::cout << "						- Input: input.t tetrahedron mesh" << std::endl;
	std::cout << "						- surface: surface.m the surface mesh of tetrahedron mesh that has been mapped into a cylinder" << std::endl;
	std::cout << "						- Output: output.t A starndard tetrahedron cylinder" << std::endl;

	std::cout << "		   -mapsurface	: map the surface of the tetrahedron mesh to a standard cylinder" << std::endl;
	std::cout << "						- Input: input.m the surface of the tetrahedron mesh" << std::endl;
	std::cout << "						- tet: input.t the tetrahedron mesh used to find the correspondence between the surface mesh and the tet mesh" << std::endl;
	std::cout << "						- uv: input.uv.m the mesh with parameter domain, used to parameterize the cylinder" << std::endl;
	std::cout << "						- Output: output.m A starndard surface cylinder" << std::endl;

	std::cout << "		   -genhex		: With given tetrahedron mesh that has been mapped to standard cylinder, generate the hex mesh with observe the original geometry of cylinder" << std::endl;
	std::cout << "						- Input: input.t tetrahedron mesh mapped onto standard cylinder." << std::endl;
	std::cout << "						- Output: output.hm The hex mesh observing the geometry of original cylinder." << std::endl;

	std::cout << "		   -optihex		: Optimize the hex mesh i.e. smoothing" << std::endl;
	std::cout << "						- Input: input.hm the input hex mesh." << std::endl;
	std::cout << "						- Output: output.hm the output optimized hex mesh." << std::endl;

	return;
};

int main(int argc, char * argv[])
{
	if (argc < 2)
	{
		std::cout << std::endl << "No Pamameters! Please use '-h' to show usage." << std::endl;
		print_help();
		return 0;
	}

	if (strcmp(argv[1], "-?") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0)
	{
		print_help();
	}

	if (strcmp(argv[1], "-mapcylinder") == 0)
	{
		if (argc < 5)
		{
			std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
			print_help();
			return 0;
		}

		map_cylinder(argc, argv);
		return 0;
	}

	if (strcmp(argv[1], "-mapsurface") == 0)
	{
		if (argc < 5)
		{
			std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
			print_help();
			return 0;
		}

		map_surface_cylinder(argc, argv);
		return 0;
	}

	if (strcmp(argv[1], "-genhex") == 0)
	{
		if (argc < 6)
		{
			std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
			print_help();
			return 0;
		}

		generate_hex(argc, argv);
		return 0;
	}

	if (strcmp(argv[1], "-mergehex") == 0)
	{
		if (argc < 6)
		{
			std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
			print_help();
			return 0;
		}

		merge_hex(argc, argv);
		return 0;
	}

	if (strcmp(argv[1], "-optihex") == 0)
	{
		if (argc < 5)
		{
			std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
			print_help();
			return 0;
		}

		optimize_hex(argc, argv);
		return 0;
	}

	if (strcmp(argv[1], "-optitet") == 0)
	{
		if (argc < 5)
		{
			std::cout << "Input/Output Not Enough. Use '-h' to show usage." << std::endl;
			print_help();
			return 0;
		}

		optimize_tet(argc, argv);
		return 0;
	}

	if (strcmp(argv[1], "-test") == 0)
	{
		test(argc, argv);
		return 0;
	}

	return 0;
}

