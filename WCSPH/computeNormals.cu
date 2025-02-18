#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <stdio.h>
#include "smoothingKernels.cuh"

__device__ void forcesInCell2(int2 neighboor, int index, double posX, double posY,  struct particleStructure* pparticles, struct paramsType* params, double* stateRates);

__global__ void computeNormals(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//	if (colorSorted[index]==0) return;  //its a boundary particle
	double posXi = pparticles->sortedX[index];
	double posYi = pparticles->sortedY[index];

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);
	int2 gridPos = { tempX,tempY };

	double stateRates[2] = { 0,0 };

	for (int y = -1; y <= 1; y++) {
		int newY = gridPos.y + y;
		if ((newY > -1) && (newY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int newX = gridPos.x + x;
				if ((newX > -1) && (newX < (*pparams).nCellsX)) {
					int2 neighboor = { newX,newY };  //2D index in grid
					forcesInCell2(neighboor, index, posXi, posYi, pparticles, pparams, stateRates);
					//forcesTemp = forcesInCell2(neighboor,index,posXi,posYi,velXi,velYi,rhoi,pressurei,pparticles,pparams);
					//forces.x += forcesTemp.x;
					//forces.y += forcesTemp.y;
				}
			}
		}

	}
	// revised - no longer requires copmute dRhoDt and SPHinfluence
	//int originalIndex = pparticles->gridParticleIndex[index];
	pparticles->sortedNormalX[index] = stateRates[0];
	pparticles->sortedNormalY[index] = stateRates[1];

	return;
}

__device__ void forcesInCell2(int2 neighboor, int index, double posXi, double posYi, struct particleStructure* pparticles, struct paramsType* pparams, double* stateRates) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;

	//required parameters
	double constWendlandD = pparams->constwendlandD;

	int startIndex = pparticles->cellStart[hash];
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];
		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double dx = (posXi - posXj);
			double dy = (posYi - posYj);
			double rSq = dx * dx + dy * dy;
			if ((rSq <= 4 * (*pparams).h2) && (rSq > 0)) {  //if they are close enough, proceede
				double h = pparams->h;
				double dist = sqrt(rSq);  //expensive but necessary
				double rOh = dist / h;
				double mj = pparticles->mass[0];
				double rhoj = pparticles->sortedRho[ind1]; //rho of sender
				double normalizedGradientInfluence = (1 / dist) * wendlandD(constWendlandD, rOh);
				stateRates[0] += h * (mj / rhoj) * normalizedGradientInfluence * dx;
				stateRates[1] += h * (mj / rhoj) * normalizedGradientInfluence * dy;
			};
		};
	};
	return;
};
				
	