// SPH3DCPPCuda.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
#include "math.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "SPH3DCPPCuda.h"
#include "SPH3DCuda.h"
#include "readNumberOfParticles.h"
#include "readFile.h"
#include "readFileImportDensity.h"
#include <string>




int main(int argc, char* argv[])
{

	//declare up front; in the heap
	std::string inputFile;
	std::string outputDir;
	std::string temp;
	std::size_t offsetOfSlash;
	size_t maxPathLength = 300;  //limit the size jsut to be safe, 300 is arbitrary

	//check the length and populate
	if (argc > 1) {
		inputFile = argv[1];
		if (inputFile.length() > maxPathLength) {
			printf("Input file name must be <= %u characters.\n", maxPathLength);
			printf("...................................Exiting\n");
			return 0;
		}

	}

	if (argc > 2) {
		outputDir = argv[2];
		if (outputDir.length() > maxPathLength) {
			printf("Output dir name must be <= %u characters.\n", maxPathLength);
			printf("...................................Exiting\n");
			return 0;
		}

	}


	//get input and output directories from command line
	switch (argc)
	{
	case 1: //no parameters
		printf("No input file");
		return 0;
	case 2: //1 parameter make outputDir same as directory where input exists
		//find the last "\"
		//keep up to that point
		//append "dataFiles\"
		offsetOfSlash = inputFile.find_last_of("\\");
		if (offsetOfSlash == std::string::npos) {
			outputDir.append(".\\");
		}
		else {
			outputDir.assign(inputFile, 0, offsetOfSlash + 1);
		}
		outputDir.append("dataFiles\\");  //append with dataFiles

		temp.append("mkdir ");
		temp.append(outputDir);
		system(temp.c_str());  //takes a const char*

		std::cout << inputFile << "\n";
		std::cout << outputDir << "\n";

		break;

	case 3: //2 parameters - use both parameters
		std::cout << inputFile << "\n";
		std::cout << outputDir << "\n";

		break;

	default:
		printf("Too many input arguments");
		return 0;

	}

	//determine CUDA capability
	int CUDAEnable;
	cudaGetDeviceCount(&CUDAEnable);
	std::cout << CUDAEnable;
	if (CUDAEnable >= 1) { CUDAEnable = 1; };
	std::cout << "\nHi,is this working?";

	//leave it as a std::vector to read in
	struct paramsType* params = new paramsType;         //parameters of the simulation
	params->importDensity = 0;
	std::vector<struct kinematicsFunctionStructure>* kinematicsFunction = new std::vector<struct kinematicsFunctionStructure>;             //vector of kinematic functions
	(*kinematicsFunction).reserve(10);  //max number of kinematics functions
	//There is an intermittant run time error associated with the kinematicsFunction.push_back operation
	//The above seems to alleviate the error

	int* nP = new int[2];
	double* xP = 0;
	int output = readNumberOfParticles(params, nP, inputFile);

	struct particleStructure* particles = new particleStructure;

	//need to create - relax particles

	int totalNumberOfParticles = nP[0] + nP[1];

	/*
	//make particle arrays
	//convert to flat arrays; just an array of particles; do not distinguish between free and constrained

	double* pX     = new double[totalNumberOfParticles];  //particle position x
	double* pY     = new double[totalNumberOfParticles];  //particle position y
	double* pVx    = new double[totalNumberOfParticles];  //particle velocity x
	double* pVy    = new double[totalNumberOfParticles];  //particle velocity y
	double* pFx    = new double[totalNumberOfParticles];  //particle force x
	double* pFy    = new double[totalNumberOfParticles];  //particle force y
	double* pR     = new double[totalNumberOfParticles];  //smoothing length
	double* pM     = new double[totalNumberOfParticles];  //mass
	double* pRho   = new double[totalNumberOfParticles];  //density
	int* pColor    = new int[totalNumberOfParticles];     //color
	*/

	//allocate memory for arrays stored in particles
	particles->x = new double[totalNumberOfParticles];
	particles->y = new double[totalNumberOfParticles];
	particles->z = new double[totalNumberOfParticles]; //3D
	particles->vx = new double[totalNumberOfParticles];
	particles->vy = new double[totalNumberOfParticles];
	particles->vz = new double[totalNumberOfParticles]; //3D
	particles->fx = new double[totalNumberOfParticles];
	particles->fy = new double[totalNumberOfParticles];
	particles->fz = new double[totalNumberOfParticles]; //3D
	particles->radius = new double[totalNumberOfParticles];
	particles->mass = new double[totalNumberOfParticles];
	particles->density = new double[totalNumberOfParticles];
	particles->color = new int[totalNumberOfParticles];
	particles->ocx = new double[totalNumberOfParticles];
	particles->ocy = new double[totalNumberOfParticles];
	particles->ocz = new double[totalNumberOfParticles]; //3D
	particles->previousX = new double[totalNumberOfParticles];
	particles->previousY = new double[totalNumberOfParticles];
	particles->previousZ = new double[totalNumberOfParticles]; //3D

	//initialize some of the memory
	for (int ind1 = 0; ind1 < totalNumberOfParticles; ind1++) {
		particles->vx[ind1] = 0;
		particles->vy[ind1] = 0;
		particles->vz[ind1] = 0;                               //3D
		particles->density[ind1] = 1000;
	};

	if (params->importDensity) {
		output = readFileImportDensity(nP, particles, params, kinematicsFunction, inputFile);
	}
	else {
		output = readFile(nP, particles, params, kinematicsFunction, inputFile);
	}


	//initialize the original constrained particle location
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		particles->ocx[ind1] = particles->x[ind1 + nP[0]];
		particles->ocy[ind1] = particles->y[ind1 + nP[0]];
		particles->ocz[ind1] = particles->z[ind1 + nP[0]];         //3D
	};


	//compute some params
	params->nTotal = nP[0] + nP[1];
	//std::cout << '\n' << params->nTotal;
	params->nFree = nP[0];
	params->nConstrained = nP[1];

	//parameters governing Tait's equation
	//following reference: Weakly compressible SPH for free surface flows, M Becker, M Teschner - Proceedings of the 2007 ACM SIGGRAPH/
	double eta = 0.01;                  //parameter governing compressibility
	double cs = params->vf / sqrt(eta);  //sound speed
//	params->gamma = 7;
//	params->B = params->rRef*cs*cs/params->gamma;

	std::cout << "\nIs it working till here?";


	if (CUDAEnable == 1) {
		SPH3DCuda(particles, params, kinematicsFunction, outputDir);
	}
	else {
		printf("No CUDA card found");
	};

	std::cout << "\nTill here?";

	return 0;
}
