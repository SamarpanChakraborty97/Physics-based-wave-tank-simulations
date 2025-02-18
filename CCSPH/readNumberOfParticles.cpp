#include "stdafx.h"
#include <iostream>
#include <vector>
#include "SPH2DCPPCuda.h"
#include "readFile.h"
#include <stdio.h>


int readNumberOfParticles(paramsType* params, int* nP, std::string inputFile) {




	float temp;
	int tempInt;
	int nQuantities;  //4 quantities for regular input, 5 quantities for version002


	FILE* pFile;

	pFile = fopen(inputFile.c_str(), "r");
	//ifstream infile(inputFile,ios::in | ios::nocreate);
	std::cout << "\nHas the problem begun?";

	fscanf(pFile, "%d", &tempInt);//read 1st line
	if (tempInt == -2) {            //indicates version 2
		params->importDensity = 1;
		nQuantities = 3;
		fscanf(pFile, "%d", &nP[0]); //store # free particles
	}
	else {
		*nP = tempInt;  //store the line read as the # free particles
		nQuantities = 2;
	}
	//endif
	std::cout << "\nHere perhaps?";

	//skip over the quantities associated with free particles
	//these are x,y, mass, smoothing length - 4 items
	for (int ind1 = 0; ind1 < nQuantities * nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
	}
	std::cout << '\n' << nP[0];
	for (int ind1 = 0; ind1 < nP[0]; ind1++) { //skip over color
		fscanf(pFile, "%d", &tempInt);
	}

	//skip over 7 parameters of free particles
	for (int ind1 = 0; ind1 < 9; ind1++) {
		fscanf(pFile, "%f", &temp);
	}

	//read number of constrained particles
	fscanf(pFile, "%d", &nP[1]); //# constrained particles
	std::cout << '\n' << nP[1];

	/*
	//skip over the quantities associated with constrained  particles
	//these are x,y, mass, smoothing length - 4 items
	for (int ind1 = 0; ind1 < nQuantities * nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
	}
	//std::cout << '\n' << nP[1];
	for (int ind1 = 0; ind1 < nP[1]; ind1++) { //skip over color
		fscanf(pFile, "%d", &tempInt);
	}
	//read number of measured particles
	fscanf(pFile, "%d", &nP[2]); //# constrained particles
	std::cout << '\n' << nP[2];
	*/

	//that's all the information needed at this point
	fclose(pFile);
	//std::cout << '\n'<<nP;


	return 0;
}