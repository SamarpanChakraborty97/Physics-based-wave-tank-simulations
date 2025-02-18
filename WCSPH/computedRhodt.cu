#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "smoothingKernels.cuh"
#include "stdio.h"

__device__ double dRhoInCell(int2 neighboor, int index, double posX, double posY, double velxi, double velyi, struct particleStructure* pparticles, struct paramsType* pparams);

__device__ void NeighGrad(int2 neighbour, int ind1, double posXj, double posYj, double velxj, double velyj, struct particleStructure* pparticles, struct paramsType* pparams, double* NeighbourGrad);

__global__ void computedRhodt(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//read particle data - host particle
	double posXi = pparticles->sortedX[index];  //these are sorted
	double posYi = pparticles->sortedY[index];
	double velxi = pparticles->sortedVx[index];
	double velyi = pparticles->sortedVy[index];

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);

	int2 gridPos = { tempX,tempY }; // grid position of host particle

	//examine neighbooring cells
	double dRho = 0;	//need density
	for (int y = -1; y <= 1; y++) {
		int currentY = gridPos.y + y;
		if ((currentY > -1) && (currentY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int currentX = gridPos.x + x;
				if ((currentX > -1) && (currentX < (*pparams).nCellsX)) {
					int2 neighboor = { currentX,currentY };  //2D index in grid
					dRho += dRhoInCell(neighboor, index, posXi, posYi, velxi, velyi, pparticles, pparams);
				}
			}
		}
	}

#if 0 //limit the denisty
	if (newDensity > 6.0 * (*pparams).rRef)
	{
		newDensity = 6.0 * (*pparams).rRef;
		//printf("corrected density");
	}
#endif

	//write drhodt to sorted position
	pparticles->sorteddRhodt[index] = dRho;



	return;
}


// loop over the particles in the host cell and surrounding cells; compute density
//__device__ double densityInCell(int2 neighboor,int index,double posX,double posY,double2* dPosSorted,double2* massRadius,int* cellStart,int* cellEnd, struct paramsType* pparams) {
__device__ double dRhoInCell(int2 neighboor, int index, double posXi, double posYi, double velxi, double velyi, struct particleStructure* pparticles, struct paramsType* pparams) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;

	double constwendlandD = pparams->constwendlandD;
	double h = pparams->h;
	double rhoi = pparticles->sortedRho[index];

	int startIndex = pparticles->cellStart[hash];
	double drdt = 0;
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];

		double L_inv[4] = { 0, 0, 0, 0 };
		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			//no reason to include self in drhodt
			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dxij = (posXi - posXj);
			double dyij = (posYi - posYj);

			double rSq = dxij * dxij + dyij * dyij;
			if ((rSq < 4 * (*pparams).h2) && (rSq > 0)) {
				double dist = sqrt(rSq);
				double rOh = dist / h;
				double normalizedGradInfluence = (1 / dist) * wendlandD(constwendlandD, rOh);
				double normX = normalizedGradInfluence * dxij;
				double normY = normalizedGradInfluence * dyij;
				L_inv[0] += -dxij * normX;
				L_inv[1] += -dxij * normY;
				L_inv[2] += -dyij * normX;
				L_inv[3] += -dyij * normY;
			}
		}

		//now we have to calculate the inverse of the above matrix or tensor
		double L[4] = { 0, 0, 0, 0 };
		double det = (L_inv[0] * L_inv[3]) - (L_inv[1] * L_inv[2]);
		L[0] = L_inv[3] / det;
		L[1] = -L_inv[1] / det;
		L[2] = -L_inv[2] / det;
		L[3] = L_inv[0] / det;

		//we have to use that inverse in the calculation of density gradient renormalization done below
		double rho_grad[2] = { 0,0 };
		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			//no reason to include self in drhodt
			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dxij = (posXi - posXj);
			double dyij = (posYi - posYj);

			double rSq = dxij * dxij + dyij * dyij;
			if ((rSq < 4 * (*pparams).h2) && (rSq > 0)) {
				double dist = sqrt(rSq);
				double rOh = dist / h;
				double normalizedGradInfluence = (1 / dist) * wendlandD(constwendlandD, rOh);
				double normX = normalizedGradInfluence * dxij;
				double normY = normalizedGradInfluence * dyij;
				double rhoj = pparticles->sortedRho[ind1];

				rho_grad[0] += (rhoj - rhoi) * (L[0] * normX + L[1] * normY);
				rho_grad[1] += (rhoj - rhoi) * (L[2] * normX + L[3] * normY);
			}
		}

		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			//no reason to include self in drhodt
			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dxij = (posXi - posXj);
			double dyij = (posYi - posYj);

			double rSq = dxij * dxij + dyij * dyij;
			if ((rSq < 4 * (*pparams).h2) && (rSq > 0)) {

				double vxj = pparticles->sortedVx[ind1];
				double vyj = pparticles->sortedVy[ind1];

				double dvxij = velxi - vxj;
				double dvyij = velyi - vyj;

				double dist = sqrt(rSq);
				double rOh = dist / pparams->h;
				double normalizedGradInfluence = (1 / dist) * wendlandD(constwendlandD, rOh);
				double normX = normalizedGradInfluence * dxij;
				double normY = normalizedGradInfluence * dyij;
				double rhoj = pparticles->sortedRho[ind1];

				int tempXj = floor((posXj - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
				int tempYj = floor((posYj - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);

				int2 gridPosj = { tempXj,tempYj }; // grid position of host particle

				double NeighbourGrad[2] = { 0,0 };
				for (int y = -1; y <= 1; y++) {
					int currentY = gridPosj.y + y;
					if ((currentY > -1) && (currentY < (*pparams).nCellsY)) {

						for (int x = -1; x <= 1; x++) {
							int currentX = gridPosj.x + x;
							if ((currentX > -1) && (currentX < (*pparams).nCellsX)) {
								int2 neighbour = { currentX,currentY };  //2D index in grid
								NeighGrad(neighbour, ind1, posXj, posYj, vxj, vyj, pparticles, pparams, NeighbourGrad);
							}
						}
					}
				}
				double term1 = rhoi * (m2 / rhoj) * (dvxij * normX + dvyij * normY);
				double term2 = (NeighbourGrad[0] + rho_grad[0]) * (-dxij) + (NeighbourGrad[1] + rho_grad[1]) * (-dyij);
				double term3 = ((-dxij * normX) + (-dyij * normY)) / rSq;
				drdt += term1 + 2 * pparams->delta * h * 10 * pparams->vf * (m2 / rhoj) * (rhoj - rhoi - 0.5 * term2) * term3;


			}; //end checking closeness				


		};  //end the for loop
	};//end checking for populated cells

	return drdt;
}

__device__ void NeighGrad(int2 neighbour, int ind1, double posXj, double posYj, double velxj, double velyj, struct particleStructure* pparticles, struct paramsType* pparams, double* NeighbourGrad) {
	int hash2 = neighbour.y * (*pparams).nCellsX + neighbour.x;

	double constwendlandD = pparams->constwendlandD;
	double h = pparams->h;
	double rhoj = pparticles->sortedRho[ind1];

	int startIndex = pparticles->cellStart[hash2];
	double drdt = 0;
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash2];

		double L_inv[4] = { 0, 0, 0, 0 };
		for (int ind2 = startIndex; ind2 < endIndex; ind2++) {
			//no reason to include self in drhodt
			double posXk = pparticles->sortedX[ind2];  //get position of sending particles
			double posYk = pparticles->sortedY[ind2];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dxjk = (posXj - posXk);
			double dyjk = (posYj - posYk);

			double rSq = dxjk * dxjk + dyjk * dyjk;
			if ((rSq < 4 * (*pparams).h2) && (rSq > 0)) {
				double dist = sqrt(rSq);
				double rOh = dist / h;
				double normalizedGradInfluence = (1 / dist) * wendlandD(constwendlandD, rOh);
				double normX = normalizedGradInfluence * dxjk;
				double normY = normalizedGradInfluence * dyjk;
				L_inv[0] += -dxjk * normX;
				L_inv[1] += -dxjk * normY;
				L_inv[2] += -dyjk * normX;
				L_inv[3] += -dyjk * normY;
			}
		}

		//now we have to calculate the inverse of the above matrix or tensor
		double L[4] = { 0, 0, 0, 0 };
		double det = (L_inv[0] * L_inv[3]) - (L_inv[1] * L_inv[2]);
		L[0] = L_inv[3] / det;
		L[1] = -L_inv[1] / det;
		L[2] = -L_inv[2] / det;
		L[3] = L_inv[0] / det;

		//we have to use that inverse in the calculation of density gradient renormalization done below
		//double rho_grad[2] = { 0,0 };
		for (int ind2 = startIndex; ind2 < endIndex; ind2++) {
			//no reason to include self in drhodt
			double posXk = pparticles->sortedX[ind2];  //get position of sending particles
			double posYk = pparticles->sortedY[ind2];
			double m2 = pparticles->mass[0];  //mass; right now these are identical  for all particles
			//compute density;  We use Monaghan's formulation with Muller's skPoly6 smoothing kernel normalized to 2D
			//The kernel is W = 
			double dxjk = (posXj - posXk);
			double dyjk = (posYj - posYk);

			double rSq = dxjk * dxjk + dyjk * dyjk;
			if ((rSq < 4 * (*pparams).h2) && (rSq > 0)) {
				double dist = sqrt(rSq);
				double rOh = dist / h;
				double normalizedGradInfluence = (1 / dist) * wendlandD(constwendlandD, rOh);
				double normX = normalizedGradInfluence * dxjk;
				double normY = normalizedGradInfluence * dyjk;
				double rhok = pparticles->sortedRho[ind2];

				NeighbourGrad[0] += (rhok - rhoj) * (L[0] * normX + L[1] * normY);
				NeighbourGrad[1] += (rhok - rhoj) * (L[2] * normX + L[3] * normY);
			}
		}
	}
}