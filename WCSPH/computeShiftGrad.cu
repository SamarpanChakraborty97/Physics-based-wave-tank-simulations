#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <iostream>
#include "smoothingKernels.cuh"
#include <stdio.h>

__device__ void shiftGradInCell2(int2 neighboor, int index, double posX, double posY, double C, struct particleStructure* pparticles, struct paramsType* pparams, double* shiftGrad);

__global__ void computeShiftGrad(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nFree) return;

	//read particle data - host particle
	double posXi = pparticles->sortedX[index];  //these are sorted, I is the receiver
	double posYi = pparticles->sortedY[index];
	double Ci = pparticles->sortedShift[index];

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);

	int2 gridPos = { tempX,tempY }; // grid position of host particle

	//examine neighbooring cells
	double shiftGrad[2] = { 0, 0 };	//shiftingCoefficient calculation

	for (int y = -1; y <= 1; y++) {
		int currentY = gridPos.y + y;
		if ((currentY > -1) && (currentY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int currentX = gridPos.x + x;
				if ((currentX > -1) && (currentX < (*pparams).nCellsX)) {
					int2 neighboor = { currentX,currentY };  //2D index in grid
					shiftGradInCell2(neighboor, index, posXi, posYi, Ci, pparticles, pparams, shiftGrad);
				}
			}
		}
	}

	//pparticles->sortedShiftGradX = shiftGrad[0];
	//pparticles->sortedShiftGradY = shiftGrad[1];

	int originalIndex = pparticles->gridParticleIndex[index];
	pparticles->shiftGradX[originalIndex] = shiftGrad[0];
	pparticles->shiftGradY[originalIndex] = shiftGrad[1];


	return;
}


// loop over the particles in the host cell and surrounding cells; compute density
//__device__ double densityInCell(int2 neighboor,int index,double posX,double posY,double2* dPosSorted,double2* massRadius,int* cellStart,int* cellEnd, struct paramsType* pparams) {
__device__ void shiftGradInCell2(int2 neighboor, int index, double posXi, double posYi, double Ci, struct particleStructure* pparticles, struct paramsType* pparams, double* shiftGrad) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;
	//double constantSpikyImprovedD = pparams->spikyImprovedD;
	double constWendlandD = pparams->constwendlandD;

	int startIndex = pparticles->cellStart[hash];

	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];

		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			//remember to include self density
			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			double rhoJ = pparticles->sortedRho[ind1]; // density of the neighbouring particle
			double Cj = pparticles->sortedShift[ind1]; //shifrting coefficient of neighbouring particle
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dx = posXi - posXj;
			double dy = posYi - posYj;
			double rSq = dx * dx + dy * dy;
			//double diffSq = (*pparams).h2 - rSq;
			if ((rSq <= 4*(*pparams).h2) && (rSq > 0)) {
				double r = sqrt(rSq);
				double rOh = r / pparams->h;
				double normalizedGradientInfluence = (1 / r) * wendlandD(constWendlandD, rOh);
				double shiftGradient = (Cj - Ci) * (m2 / rhoJ) * normalizedGradientInfluence;
				shiftGrad[0] += shiftGradient * dx;
				shiftGrad[1] += shiftGradient * dy;
				//debug
				int index = blockIdx.x * blockDim.x + threadIdx.x;
				int origIndex = pparticles->gridParticleIndex[index];
				if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
					int senderIdx = pparticles->gridParticleIndex[ind1];
					printf("Ini D; it=%u, P# %u -> w/ %u at %f and %f \n", pparams->ind1, pparams->DEBUGpNum, senderIdx, posXj, posYj);
				}


			}; //end checking closeness
		};  //end the for loop
	};//end the if statement

	return;
}