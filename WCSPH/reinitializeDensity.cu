#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"
#include "smoothingKernels.cuh"


__device__ double2 computeComponents(int2 neighboor, int index, double posX, double posY, struct particleStructure* pparticles, struct paramsType* pparams);


__global__ void reinitializeDensity(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//read particle data - host particle
	double posXi = pparticles->sortedX[index];  //these are sorted
	double posYi = pparticles->sortedY[index];

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);

	int2 gridPos = { tempX,tempY }; // grid position of host particle

	//examine neighbooring cells
	double2 numDenom = { 0,0 };	//need 
	double2 temp = { 0,0 };
	for (int y = -1; y <= 1; y++) {
		int currentY = gridPos.y + y;
		if ((currentY > -1) && (currentY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int currentX = gridPos.x + x;
				if ((currentX > -1) && (currentX < (*pparams).nCellsX)) {
					int2 neighboor = { currentX,currentY };  //2D index in grid
					temp = computeComponents(neighboor, index, posXi, posYi, pparticles, pparams);
					numDenom.x += temp.x;
					numDenom.y += temp.y;
				}
			}
		}
	}




	//if particles exceede the boundaries they will have denom = 0 and rhoi = 0;
	//In this case, they are no longer an important part of the calculation, so 
	//set their rhoi = density of a single particle

	//rhoi = m2*(*pparams).constDensity*(1-rOhSq)*(1-rOhSq)*(1-rOhSq);
	// = m2*(*pparams).constDensity*(1-0)*(1-0)*(1-0);
	// = m2*(*pparams).constDensity;

	//debug
	int origIndex = pparticles->gridParticleIndex[index];
	if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
		printf("reini rho; it=%u, P# %u; num %f, denom %f \n", pparams->ind1, origIndex, numDenom.x, numDenom.y);
	}

	//double filteredDensity = numerator/denom;
	double filteredDensity = numDenom.x / numDenom.y;

	if (numDenom.y == 0)
	{
		//maybe keeping it the same is the answer
		filteredDensity = 2 * pparticles->mass[0] * (*pparams).constDensity;
		printf("particle exceeded boundary\n");
	}


	if (filteredDensity == 0) {
		printf("something is wrong\n");
	}


	//store the density in a temporary array
	pparticles->sortedRhoFiltered[index] = filteredDensity;  //sortedRhoFiltered no longer exists

	//write new filtered density back to original unsorted position
	//int originalIndex = pparticles->gridParticleIndex[index];
	//pparticles->unsortedRhoFiltered[originalIndex] = rhoi;


}


// loop over the particles in the host cell and surrounding cells; compute density
//__device__ double densityInCell(int2 neighboor,int index,double posX,double posY,double2* dPosSorted,double2* massRadius,int* cellStart,int* cellEnd, struct paramsType* pparams) {
__device__ double2 computeComponents(int2 neighboor, int index, double posXi, double posYi, struct particleStructure* pparticles, struct paramsType* pparams) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;

	int startIndex = pparticles->cellStart[hash];
	double2 numDenom = { 0,0 };
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];

		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			//no reason to include self in drhodt
			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles

			//The kernel is W = 
			double dxij = (posXi - posXj);
			double dyij = (posYi - posYj);

			double rSq = dxij * dxij + dyij * dyij;
			//double diffSq = 4*(*pparams).h2 - rSq;
			if ((rSq <= 4*(*pparams).h2) && (rSq >= 0))
			{
				double rhoj = pparticles->sortedRho[ind1];
				double dist = sqrt(rSq);
				double rOh = dist / pparams->h;

				double kernelInfluence = wendland(pparams->constwendland, rOh);
				numDenom.x += kernelInfluence * m2;        //numerator
				numDenom.y += kernelInfluence * m2 / rhoj;   //denominator




			};  //end the for loop
		};//end the if statement
	};

	return numDenom;
}