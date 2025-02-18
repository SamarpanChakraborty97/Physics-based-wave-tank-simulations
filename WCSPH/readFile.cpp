#include "stdafx.h"
#include <iostream>
#include <vector>
#include "SPH2DCPPCuda.h"
#include "readFile.h"



int readFile(int* nP, particleStructure* particles, paramsType* params, std::vector<struct kinematicsFunctionStructure>* kinematicsFunction, std::string inputFile) {

	float temp;
	int iTemp;

	FILE* pFile;


	pFile = fopen(inputFile.c_str(), "r");

	//read free particles
	fscanf(pFile, "%d", &nP[0]);  //re-read number of free particles
	//read x values
	//write dummy values for the other values

	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->x[ind1] = temp;
		//std::cout << '\n' << particles->y[ind1]
	};
	std::cout << "\nIs it still working?";
	//read y values
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->y[ind1] = temp;
		//std::cout << '\n' << particles->y[ind1];
	};

	/*****	NOT READING IN VELOCITY IN THIS VERSION*****
	//read vx values
	for (int ind1=0;ind1<nP[0];ind1++){
		fscanf_s(pFile,"%f",&temp);
		particles->vx[ind1] = temp;
	};

	//read vy values
	for (int ind1=0;ind1<nP[0];ind1++){
		fscanf_s(pFile,"%f",&temp);
		particles->vy[ind1] = temp;
	};
	*/

	std::cout << "\nIs it failing here?";

	//read mass
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->mass[ind1] = temp;
	};

	//read radius
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->radius[ind1] = temp;
	};

	//read color
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%d", &iTemp);
		particles->color[ind1] = iTemp;
	};

	//free particle parameters
	fscanf(pFile, "%f", &temp);
	params->vf = temp;  //maximum fluid velocity, used to determine sound speed; Tait's equation
	std::cout << '\n' << params->vf;

	fscanf(pFile, "%f", &temp);
	params->nu = temp;  //viscosity

	fscanf(pFile, "%f", &temp);
	params->delta = temp;  //delta-SPH parameter

	fscanf(pFile, "%f", &temp);
	params->viscoAlpha = temp;  //viscoAlpha

	fscanf(pFile, "%f", &temp);
	params->viscoBeta = temp;  //visco Beta

	fscanf(pFile, "%f", &temp);
	params->epsilon = temp;  //epsilon - XSPH epsilon

	fscanf(pFile, "%f", &temp);
	params->rRef = temp; //reference density



	//load constrained particles
	fscanf(pFile, "%d", &nP[1]);


	//read x-position
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->x[nP[0] + ind1] = temp;
	};

	//read y-position
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->y[nP[0] + ind1] = temp;
	};

	//read mass
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->mass[nP[0] + ind1] = temp;
	};
	std::cout << '\n' << temp;

	//read smoothing length
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->radius[nP[0] + ind1] = temp;
	};

	//read color
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%d", &iTemp);
		particles->color[nP[0] + ind1] = iTemp;
	};

	/*
	//load measured particles
	fscanf(pFile, "%d", &nP[2]);


	//read x-position
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->x[nP[0] + nP[1] + ind1] = temp;
	};

	//read y-position
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->y[nP[0] + nP[1] + ind1] = temp;
	};

	//read mass
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->mass[nP[0] + nP[1] + ind1] = temp;
	};

	//read smoothing length
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->radius[nP[0] + nP[1] +ind1] = temp;
	};

	//read color
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%d", &iTemp);
		particles->color[nP[0] + nP[1] + ind1] = iTemp;
	};
	*/

	//read kinematic functions: used as current = original+function
	fscanf(pFile, "%d", &iTemp);
	std::cout << '\n' << iTemp;
	struct kinematicsFunctionStructure kinematicsTemp;
	if (iTemp == 0) {
		kinematicsTemp.xFunctionType = 10;
		kinematicsTemp.yFunctionType = 10;
		kinematicsTemp.r1 = 1;
		kinematicsTemp.r2 = 2;
		(*kinematicsFunction).push_back(kinematicsTemp);
	}

	int ind1 = 0;
	while (iTemp != 0) {
		kinematicsTemp.xFunctionType = iTemp;

		if (iTemp == 0) //no more functions
			break;

		else if (iTemp == 1)     //unitStep - A*unitStep(t-B);   A; B;
		{
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

		}

		else if (iTemp == 2) { //linear   - A*t;               A
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;
		}

		else if (iTemp == 3) { //sinusoidal - A*sin(2*pi*t*B+C); A; B; C
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3x = temp;

			//fscanf(pFile, "%f", &temp);
			//kinematicsTemp.p4x = temp;

			//fscanf(pFile, "%f", &temp);
			//kinematicsTemp.p5x = temp;
		}

		else if (iTemp == 4) {//rotation - xpt, ypt, theta
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3x = temp;
		}
		else if (iTemp == 5) { // sloshing motion
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p4x = temp;
		}
		else if (iTemp == 6) {  //peregrine soliton
			fscanf(pFile, "%f", &temp); //a0
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp); //wavelength
			kinematicsTemp.p2x = temp;

			//fscanf(pFile, "%f", &temp); 
			//kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp); //lead Time
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp); //target distance
			kinematicsTemp.p4x = temp;

			fscanf(pFile, "%f", &temp); //proportionality
			kinematicsTemp.p5x = temp;

		}
		else if (iTemp == 7) {  //rotational peregrine soliton
			fscanf(pFile, "%f", &temp); //a0
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp); //l
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp); //leadTime
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp); //targetDistance
			kinematicsTemp.p4x = temp;

			fscanf(pFile, "%f", &temp); //centerOfRotationX
			kinematicsTemp.p5x = temp;

			fscanf(pFile, "%f", &temp); //centerOfRotationY
			kinematicsTemp.p6x = temp;

			fscanf(pFile, "%f", &temp); //amplitude Factor
			kinematicsTemp.p7x = temp;

			fscanf(pFile, "%f", &temp); //bias angle
			kinematicsTemp.p8x = temp;


		}

		else if (iTemp == 10)     //Do nothing
		{ //no parameters
		}


		//y-function//
		fscanf(pFile, "%d", &iTemp);   //scan again for y function
		kinematicsTemp.yFunctionType = iTemp;
		if (iTemp == 1)     //unitStep - A*unitStep(t-B);   A; B;
		{
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;

		}

		else if (iTemp == 2) { //linear   - A*t;               A
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;
		}

		else if (iTemp == 3) { //sinusoidal - A*sin(2*pi*t*B+C); A; B; C
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3y = temp;
		}
		else if (iTemp == 4) {//rotation - xpt, ypt, theta
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3y = temp;
		}
		else if (iTemp == 5) { //sloshing motion
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;
		}
		else if (iTemp == 6) {  //peregrine soliton
			fscanf(pFile, "%f", &temp); //a0
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp); //l
			kinematicsTemp.p2y = temp;

			fscanf(pFile, "%f", &temp); //leadTime
			kinematicsTemp.p3y = temp;

			fscanf(pFile, "%f", &temp); //targetDistance
			kinematicsTemp.p4y = temp;

		}



		else if (iTemp == 10)     //Do nothing
		{ //no parameters
		}

		fscanf(pFile, "%d", &iTemp);   //range 1
		kinematicsTemp.r1 = iTemp - 1;  //convert from 1's based index to 0's based index

		fscanf(pFile, "%d", &iTemp);   //range 2
		kinematicsTemp.r2 = iTemp - 1;

		//store the result
		(*kinematicsFunction).push_back(kinematicsTemp);

		fscanf(pFile, "%d", &iTemp);
		ind1++; //increment the counter
	};




	//computational parameters
	fscanf(pFile, "%f", &temp);
	params->gravity = temp;

	fscanf(pFile, "%d", &iTemp);
	params->nTime = iTemp;

	fscanf(pFile, "%f", &temp);
	params->dt = temp;

	fscanf(pFile, "%d", &iTemp);
	params->storageStride = iTemp;

	fclose(pFile);




	return 0;
}