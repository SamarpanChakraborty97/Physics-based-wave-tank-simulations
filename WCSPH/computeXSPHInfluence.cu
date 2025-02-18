#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <stdio.h>

__device__ double2 XSPHVelocityInCell(int2 neighboor, int index, double posX, double posY, double velxi, double velyi, double rhoi, struct particleStructure* pparticles, struct paramsType* pparams);

__global__ void computeXSPHInfluence(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//this is accessing the sorted velocities; but the unsorted veolicties have been updated with the new accelerations
	//it should access the unsorted velocities.


	double epsilon = (*pparams).epsilon;  //parameter for XSPH influence

	//read particle data - host particle
	double posX = pparticles->sortedX[index];  //these are sorted
	double posY = pparticles->sortedY[index];

	double velxi = pparticles->sortedVx[index];  //these are sorted
	double velyi = pparticles->sortedVy[index];

	double rhoi = pparticles->sortedRho[index];

	//get address in grid
	int tempX = floor((posX - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posY - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);

	int2 gridPos = { tempX,tempY }; // grid position of host particle

	//examine neighbooring cells
	double2 XSPHVel = { 0,0 };	//need XSPH
	double2 temp = { 0,0 };
	for (int y = -1; y <= 1; y++) {
		int currentY = gridPos.y + y;
		if ((currentY > -1) && (currentY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int currentX = gridPos.x + x;
				if ((currentX > -1) && (currentX < (*pparams).nCellsX)) {
					int2 neighboor = { currentX,currentY };  //2D index in grid
					temp = XSPHVelocityInCell(neighboor, index, posX, posY, velxi, velyi, rhoi, pparticles, pparams);
					XSPHVel.x += temp.x;
					XSPHVel.y += temp.y;
				}
			}
		}
	}

	//write XSPH density back to original unsorted position
	int originalIndex = pparticles->gridParticleIndex[index];

	pparticles->XSPHVelX[originalIndex] = epsilon * XSPHVel.x;
	pparticles->XSPHVelY[originalIndex] = epsilon * XSPHVel.y;

	//DEBUG
	int origIndex = pparticles->gridParticleIndex[index];
	if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
		printf("XSPHx %f, XSPHY %f \n", pparticles->XSPHVelX[origIndex], pparticles->XSPHVelY[origIndex]);
	}



	return;
}


// loop over the particles in the host cell and surrounding cells; compute density
//__device__ double densityInCell(int2 neighboor,int index,double posX,double posY,double2* dPosSorted,double2* massRadius,int* cellStart,int* cellEnd, struct paramsType* pparams) {
__device__ double2 XSPHVelocityInCell(int2 neighboor, int index, double posXi, double posYi, double velxi, double velyi, double rhoi, struct particleStructure* pparticles, struct paramsType* pparams) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;

	int startIndex = pparticles->cellStart[hash];
	double2 XSPHvel = { 0,0 };
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];

		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {

			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dxji = (posXj - posXi);
			double dyji = (posYj - posYi);
			double rSq = dxji * dxji + dyji * dyji;
			double diffSq = (*pparams).h2 - rSq;
			if (diffSq > 0) {
				//don't include constrained particles
			//	int originalIndex = pparticles->gridParticleIndex[ind1];
			//  if (originalIndex<pparams->nFree) {
				if (1) {

					double h = pparams->h;
					double rhoj = pparticles->sortedRho[ind1];
					double vxj = pparticles->sortedVx[ind1];
					double vyj = pparticles->sortedVy[ind1];

					double dvxji = vxj - velxi;
					double dvyji = vyj - velyi;

					double rhoBarij = (rhoi + rhoj) / 2;
					double rOh = sqrt(rSq) / h;
					double rOhSq = rOh * rOh;
					double mutualInfluence = m2 / rhoBarij * pparams->constDensity * (1 - rOhSq) * (1 - rOhSq) * (1 - rOhSq);

					XSPHvel.x += mutualInfluence * dvxji;
					XSPHvel.y += mutualInfluence * dvyji;

				}; //end excluding constrained particles
			}; //end checking closeness
		};  //end the for loop
	};//end the if statement - populated cells

	return XSPHvel;
}