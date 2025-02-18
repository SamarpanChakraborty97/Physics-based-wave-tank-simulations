#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <stdio.h>
#include "smoothingKernels.cuh"

//__device__ double poly6(double cd, double r);
//__device__ double spikyImprovedD(double cd, double r);
//__device__ double computePressure(double rhoi, double rhoRef, double tenVMaxSq);

__device__ void computeGrads(int2 neighboor, double posX, double posY, struct paramsType* params, struct particleStructure* particles, double* changes);

__global__ void computeKernelCorrectionTerms(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	double det_threshold = 0.6;
	double rhoRef = pparams->rRef;
	double cSquared = pparams->cSound * pparams->cSound;

	//read primary particle data - this is sorted data
	double posXi = pparticles->sortedX[index];
	double posYi = pparticles->sortedY[index];

	int origIndex = pparticles->gridParticleIndex[index];

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);
	int2 gridPos = { tempX,tempY };

	//examine neighbooring cells
	double changes[3] = { 0,0,0 };  //pointer to array of {fx,fy,XSPHx,XSPHy}

	for (int y = -1; y <= 1; y++) {
		int newY = gridPos.y + y;
		if ((newY > -1) && (newY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int newX = gridPos.x + x;
				if ((newX > -1) && (newX < (*pparams).nCellsX)) {
					int2 neighboor = { newX,newY };  //2D index in grid
					computeGrads(neighboor, posXi, posYi, pparams, pparticles, changes);

				}
			}
		}

	}


	double det = changes[0] * changes[2] - changes[1] * changes[1];

	pparticles->det_values[index] = det;

	if (det < det_threshold) {
		pparticles->sortedA11[index] = 1;
		pparticles->sortedA12[index] = 0;
		pparticles->sortedA22[index] = 1;
	}
	else {
		pparticles->sortedA11[index] = -changes[0];
		pparticles->sortedA12[index] = -changes[1];
		pparticles->sortedA22[index] = -changes[2];
	}

	//if (index == 2500) {
	//	printf("A11 for the particle with index %d is %f\n", index, pparticles->sortedA11[index]);
	//}


	//printf("%d\n",pparticles->sorteddRhodt[index]);

	//debug
#if (0)
	int orig = pparticles->gridParticleIndex[index];
	if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
		printf("AP; it=%u, P# %u; fx %f, fy %f, ", pparams->ind1, origIndex, pparticles->fx[originalIndex], pparticles->fy[originalIndex]);
	}
#endif


	return;
}
/*
__location__(global) void computePressure(particleStructure* pdParticles, paramsType* dParams)
{
	return __location__(global) void();
}
*/

// loop over the particles in the host cell and surrounding cells; compute density
//__device__ double2 forcesInCell(int2 neighboor,int index,double posX,double posY,double velX, double velY, double rho,double2* posSorted,double2* velSorted,int* colorSorted, double2* dXSPHVelTemp,double2* massRadius,double* pRhoSorted,int* cellStart,int* cellEnd,const paramsType* params) {
__device__ void computeGrads(int2 neighboor, double posXi, double posYi, struct paramsType* pparams, struct particleStructure* pparticles, double* changes) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;

	double constWendlandD = pparams->constwendlandD;

	int startIndex = pparticles->cellStart[hash];
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];
		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {

			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double dx = (posXj - posXi);
			double dy = (posYj - posYi);
			double rSq = dx * dx + dy * dy;

			if ((rSq < 4 * pparams->h2) && (rSq > 0)) {
				double r = sqrt(rSq);
				double rOh = r / pparams->h;
				double Grad = (1 / r) * wendlandD(pparams->constwendlandD, rOh);

				double rhoj = pparticles->sortedRho[ind1];

				changes[0] += (dx * dx) * Grad * pparams->mass / rhoj;
				changes[1] += (dx * dy) * Grad * pparams->mass / rhoj;
				changes[2] += (dy * dy) * Grad * pparams->mass / rhoj;

			};
		};
	};

	return;
}