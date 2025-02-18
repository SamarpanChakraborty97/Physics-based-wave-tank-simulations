#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <iostream>
#include "smoothingKernels.cuh"
#include <stdio.h>

__device__ void mlsInCell2(int2 neighboor, int index, double posX, double posY, struct particleStructure* pparticles, struct paramsType* pparams, double* L);

__global__ void computeMLS(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//read particle data - host particle
	double posXi = pparticles->sortedX[index];  //these are sorted, I is the receiver
	double posYi = pparticles->sortedY[index];

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);

	int2 gridPos = { tempX,tempY }; // grid position of host particle

	//examine neighbooring cells
	double L[4] = { 0, 0, 0, 0 };	//shiftingCoefficient calculation

	for (int y = -1; y <= 1; y++) {
		int currentY = gridPos.y + y;
		if ((currentY > -1) && (currentY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int currentX = gridPos.x + x;
				if ((currentX > -1) && (currentX < (*pparams).nCellsX)) {
					int2 neighboor = { currentX,currentY };  //2D index in grid
					mlsInCell2(neighboor, index, posXi, posYi, pparticles, pparams, L);
				}
			}
		}
	}

	//pparticles->sortedShiftGradX = shiftGrad[0];
	//pparticles->sortedShiftGradY = shiftGrad[1];

	//int originalIndex = pparticles->gridParticleIndex[index];
	double det = (L[0] * L[3]) - (L[1] * L[2]);
	//printf("%d\n",det);
	//printf("%d\n", L[0]);
	//pparticles->L1[index] = L[3] / det;
	if (det == 0) {
		pparticles->L1[index] = 0;
		pparticles->L2[index] = 0;
		pparticles->L3[index] = 0;
		pparticles->L4[index] = 0;
	}
	else {
		pparticles->L1[index] = L[3] / det;
		pparticles->L2[index] = -L[1] / det;
		pparticles->L3[index] = -L[2] / det;
		pparticles->L4[index] = L[0] / det;
	}
	return;
}


// loop over the particles in the host cell and surrounding cells; compute density
//__device__ double densityInCell(int2 neighboor,int index,double posX,double posY,double2* dPosSorted,double2* massRadius,int* cellStart,int* cellEnd, struct paramsType* pparams) {
__device__ void mlsInCell2(int2 neighboor, int index, double posXi, double posYi, struct particleStructure* pparticles, struct paramsType* pparams, double* L) {

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
			//double Cj = pparticles->sortedShift[ind1]; //shifrting coefficient of neighbouring particle
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
				//double shiftGradient = (Cj - Ci) * (m2 / rhoJ) * normalizedGradientInfluence;
				L[0] += (m2 / rhoJ) * (-dx * normalizedGradientInfluence * dx);
				L[1] += (m2 / rhoJ) * (-dx * normalizedGradientInfluence * dy);
				L[2] += (m2 / rhoJ) * (-dy * normalizedGradientInfluence * dx);
				L[3] += (m2 / rhoJ) * (-dy * normalizedGradientInfluence * dy);
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