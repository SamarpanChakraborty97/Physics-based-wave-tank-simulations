#include <stdio.h>
#include "stdafx.h"
#include <iostream>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include "SPH2DCPPCuda.h"

int exportSystem(struct particleStructure* particles, unsigned int dataFileNumber, struct paramsType* params, std::string outputDir) {

	int nTotal = (*params).nTotal;
	//std::cout << '\n' << nTotal;

	FILE* pFile;
	// const char path[] = "../../dataFiles/"; defined as a variable in header
	const char base[] = "dataFile";
	const char suffix[] = ".txt";
	char filename[350];

	sprintf(filename, "%s%s%05d%s", outputDir.c_str(), base, dataFileNumber, suffix);

	pFile = fopen(filename, "w");

	//write all particles
	for (int ind1 = 0; ind1 < nTotal; ind1++) {

		fprintf(pFile, "%16.16f ", particles->x[ind1]);
		fprintf(pFile, "%16.16f ", particles->y[ind1]);
		//fprintf(pFile, "%16.16f ", particles->vx[ind1]);
		//fprintf(pFile, "%16.16f ", particles->vy[ind1]);
		fprintf(pFile, "%16.16f ", particles->density[ind1]);
		//fprintf(pFile, "%16.16f ", particles->pressure[ind1]);
		//fprintf(pFile, "%16.16f ", particles->rhoGradX[ind1]);
		//fprintf(pFile, "%16.16f ", particles->L1[ind1]);
		if (ind1 < nTotal - 1) {
			fprintf(pFile, "%16.16f\n", (double)particles->color[ind1]);
		}
		else {
			fprintf(pFile, "%16.16f", (double)particles->color[ind1]);
		}

	};



	fclose(pFile); //fclose deletes the pointer
	//std::cout << "\nIs it working here?";

	return 0;
}