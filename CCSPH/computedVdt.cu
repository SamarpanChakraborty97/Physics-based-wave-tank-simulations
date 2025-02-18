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

__device__ void forcesInCell2(int2 neighboor, int index, double posX, double posY, double velX, double velY, double rho, double pressurei, double a11, double a12, double a22, struct particleStructure* pparticles, struct paramsType* params, double* stateRates);
//__device__ void forcesInCell2(int2 neighboor, int index, double posX, double posY, double velX, double velY, double rho, double pressurei, struct particleStructure* pparticles, struct paramsType* params, double* stateRates);


__global__ void computedVdt(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//	if (colorSorted[index]==0) return;  //its a boundary particle
	double rhoRef = pparams->rRef;
	double cSquared = pparams->cSound * pparams->cSound;

	//read primary particle data - this is sorted data
	double posXi = pparticles->sortedX[index];
	double posYi = pparticles->sortedY[index];
	double velXi = pparticles->sortedVx[index];
	double velYi = pparticles->sortedVy[index];
	double rhoi = pparticles->sortedRho[index];
	//double rhoxi = pparticles->rhoGradX[index];
	//double rhoyi = pparticles->rhoGradY[index];
	double pressurei = computePressure(rhoi, rhoRef, cSquared);

	double a11I = pparticles->sortedA11[index];
	double a12I = pparticles->sortedA12[index];
	double a22I = pparticles->sortedA22[index];

	//debug
	int origIndex = pparticles->gridParticleIndex[index];
	if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
		printf("AP; it=%u, P# %u press. = %f, rhoi = %f, xi = %f, yi = %f \n", pparams->ind1, origIndex, pressurei, rhoi, posXi, posYi);
	}

	//get address in grid
	int tempX = floor((posXi - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
	int tempY = floor((posYi - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);
	int2 gridPos = { tempX,tempY };

	//examine neighbooring cells
	double stateRates[5] = { 0,0,0,0,0 };  //pointer to array of {fx,fy,XSPHx,XSPHy}

	for (int y = -1; y <= 1; y++) {
		int newY = gridPos.y + y;
		if ((newY > -1) && (newY < (*pparams).nCellsY)) {

			for (int x = -1; x <= 1; x++) {
				int newX = gridPos.x + x;
				if ((newX > -1) && (newX < (*pparams).nCellsX)) {
					int2 neighboor = { newX,newY };  //2D index in grid
					forcesInCell2(neighboor, index, posXi, posYi, velXi, velYi, rhoi, pressurei, a11I, a12I, a22I, pparticles, pparams, stateRates);
					//forcesInCell2(neighboor, index, posXi, posYi, velXi, velYi, rhoi, pressurei, pparticles, pparams, stateRates);
					//forcesTemp = forcesInCell2(neighboor,index,posXi,posYi,velXi,velYi,rhoi,pressurei,pparticles,pparams);
					//forces.x += forcesTemp.x;
					//forces.y += forcesTemp.y;
				}
			}
		}

	}

	/*
	//these are actually accelerations according to A. Vorobyev thesis
	//acceleration due to gravity is added at this point
	//write new forces back to original unsorted position
	int originalIndex = pparticles->gridParticleIndex[index];
	pparticles->fx[originalIndex] = forces.x;  //{x,y} //check sign
	pparticles->fy[originalIndex] = forces.y;  //{x,y}
	*/


	// revised - no longer requires copmute dRhoDt and SPHinfluence
	int originalIndex = pparticles->gridParticleIndex[index];
	pparticles->fx[originalIndex] = stateRates[0];
	pparticles->fy[originalIndex] = stateRates[1];

	pparticles->XSPHVelX[originalIndex] = stateRates[2];
	pparticles->XSPHVelY[originalIndex] = stateRates[3];

	pparticles->sorteddRhodt[index] = stateRates[4];
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
__device__ void forcesInCell2(int2 neighboor, int index, double posXi, double posYi, double velXi, double velYi, double rhoi,  double pressurei, double a11I, double a12I, double a22I, struct particleStructure* pparticles, struct paramsType* pparams, double* stateRates) {
//__device__ void forcesInCell2(int2 neighboor, int index, double posXi, double posYi, double velXi, double velYi, double rhoi, double pressurei, struct particleStructure* pparticles, struct paramsType* pparams, double* stateRates) {

	//compute 1D hash value
	int hash = neighboor.y * (*pparams).nCellsX + neighboor.x;

	//required parameters
	double rhoRef = pparams->rRef;
	double cSquared = pparams->cSound * pparams->cSound;
	//double constantSpikyImprovedD = pparams->spikyImprovedD;
	//double constPoly6 = pparams->constDensity;
	double constWendland = pparams->constwendland;
	double constWendlandD = pparams->constwendlandD;

	int startIndex = pparticles->cellStart[hash];
	//double2 forces = {0,0};
	if (startIndex != 0xffffffff) {
		int endIndex = pparticles->cellEnd[hash];
		for (int ind1 = startIndex; ind1 < endIndex; ind1++) {
			//remember to exclude self-force; stay within desired domain

			/*
			//DEBUG
			if (ind1==index) {
				printf("index == ind1, x position is %f  %f\n", pparticles->sortedX[ind1], posXi);
			}
			*/


			//if (ind1 != index) {

				//many SPH references cite pAB = pA - pB; where A is the primary particle.
			//following this convention we have
			//This gives a vector pointing from particle B to particle A

			double posXj = pparticles->sortedX[ind1];  //get position of sending particles
			double posYj = pparticles->sortedY[ind1];
			double dx = (posXi - posXj);
			double dy = (posYi - posYj);
			double rSq = dx * dx + dy * dy;

			if ((rSq <= 4*(*pparams).h2) && (rSq > 0)) {  //if they are close enough, proceede

				double h = pparams->h;
				double dist = sqrt(rSq);  //expensive but necessary
				double dvxij = velXi - pparticles->sortedVx[ind1];  //velocity of sending particle
				double dvyij = velYi - pparticles->sortedVy[ind1];
				double rOh = dist / h;
				double mj = pparams->mass;
				double rhoj = pparticles->sortedRho[ind1]; //rho of sender
				//double rhoxj = pparticles->rhoGradX[ind1];
				//double rhoyj = pparticles->rhoGradY[ind1];
				double pressurej = computePressure(rhoj, rhoRef, cSquared);
				//double Cij = sqrt(tenVMaxSq);

				double normalizedGradientInfluence = (1 / dist) * wendlandD(constWendlandD, rOh);
				/*

				double vMorTiInner = mj * (pparams->nu) * (rhoi + rhoj) / (rhoi * rhoj) * normalizedGradientInfluence;
				double vMorTiX = vMorTiInner * dvxij;
				double vMorTiY = vMorTiInner * dvyij;
				double dirVel = dvxij * dx + dvyij * dy;
				double rhoBarij = (rhoi + rhoj) / 2;      //used in XSPH as well

				double vTix = vMorTiX;  //it gets vTi regardless
				double vTiy = vMorTiY;

				if (dirVel < 0)   //it may get additional terms
				{
					//				double muij        = h*dirVel/(rSq+0.01*h*h);
					double muij = h  * dirVel / rSq ;

					double addedViscosity1 = mj * pparams->viscoBeta * muij * muij / rhoBarij * normalizedGradientInfluence;
					double addedViscosity2 = mj * (-pparams->viscoAlpha) * Cij* muij / rhoBarij * normalizedGradientInfluence;
					double addedViscosity = addedViscosity1 + addedViscosity2;
					vTix += addedViscosity * dx ;  //already has 1/|rij|
					vTiy += addedViscosity * dy;
				}

				double sharedTerm = mj * (pressurei / (rhoi * rhoi) + pressurej / (rhoj * rhoj)) * normalizedGradientInfluence;
				//double sharedTerm = mj * (pressurei / (rhoi * rhoi) + pressurej / (rhoj * rhoj)) * spikyImprovedD(constantSpikyImprovedD, rOh);

				double term1X = sharedTerm * dx;
				double term1Y = sharedTerm * dy;


				//these are actually accelerations
				stateRates[0] += -term1X ;  //fx
				stateRates[1] += -term1Y ;  //fy

				*/

				double vMorTiInner = mj * (pparams->nu) * (rhoi + rhoj) / (rhoi * rhoj) * normalizedGradientInfluence;
				double vTiX = vMorTiInner * dvxij;
				double vTiY = vMorTiInner * dvyij;

				double a11J = pparticles->sortedA11[ind1];
				double a12J = pparticles->sortedA12[ind1];
				double a22J = pparticles->sortedA22[ind1];

				double A11 = 0.5 * (a11I + a11J);
				double A12 = 0.5 * (a12I + a12J);
				double A22 = 0.5 * (a22I + a22J);

				double det = A11 * A22 - A12 * A12;
				double B11 = (1 / det) * A22;
				double B12 = -(1 / det) * A12;
				double B22 = (1 / det) * A11;

				double dirVel = dvxij * dx + dvyij * dy;
				double rhoBarij = (rhoi + rhoj) / 2;      //used in XSPH as well
				double Cij = pparams->cSound;  //mean speed of sound

				double sharedTerm = (1/rhoi) * (mj/rhoj) * (pressurei + pressurej) * normalizedGradientInfluence;
				//double sharedTerm = mj * ((pressurei / (rhoi * rhoi)) + (pressurej / (rhoj * rhoj))) * normalizedGradientInfluence;

				//double term1X = sharedTerm * dx;
				//double term1Y = sharedTerm * dy;
				double term1X = sharedTerm * (B11 * dx + B12 * dy);
				double term1Y = sharedTerm * (B12 * dx + B22 * dy);


				if (dirVel < 0)   //it may get additional terms
				{
					//				double muij        = h*dirVel/(rSq+0.01*h*h);
					double muij = h * dirVel / rSq;

					double addedViscosity1 = mj * pparams->viscoBeta * muij * muij / rhoBarij * normalizedGradientInfluence;
					//double addedViscosity2 = pparams->rRef * (mj/rhoj) * (-pparams->viscoAlpha) * Cij* muij / rhoBarij * normalizedGradientInfluence;
					double addedViscosity2 = mj * (-pparams->viscoAlpha) * Cij * muij / rhoBarij * normalizedGradientInfluence;
					double addedViscosity = addedViscosity1 + addedViscosity2;
					//term1X += addedViscosity * dx;  //already has 1/|rij| 
					//term1Y += addedViscosity * dy;
					term1X += addedViscosity * (B11 * dx + B12 * dy);  //already has 1/|rij| 
					term1Y += addedViscosity * (B12 * dx + B22 * dy);
				}

				//these are actually accelerations
				//stateRates[0] += -term1X; //fx
				//stateRates[1] += -term1Y;  //fy

				stateRates[0] += -term1X + vTiX; //fx
				stateRates[1] += -term1Y + vTiY;  //fy

				//stateRates[0] += -term1X * (B11 * dx + B12 * dy) + vTiX; //fx
				//stateRates[1] += -term1Y * (B12 * dx + B22 * dy) + vTiY;  //fy

				//		forces.x += -term1X+vTix;
				//		forces.y += -term1Y+vTiy;


				//XSPH
				double mutualInfluence = pparams->epsilon * mj / rhoBarij * wendland(constWendland, rOh);
				stateRates[2] += -mutualInfluence * dvxij;  //XSPHx; dxvij = -dvxji; dvxji is called for in the definition
				stateRates[3] += -mutualInfluence * dvyij;  //XSPHy; dxvij = -dvxji; dvxji is called for in the definition

				//DRho/dt
				double term1 = mj * normalizedGradientInfluence * dirVel;
				//double term2 = (rhoj - rhoi) - 0.5 * (((rhoxi + rhoxj) * (-dx)) + ((rhoyi + rhoyj) * (-dy)));
				double term2 = (rhoj - rhoi);
				//double term2 = rho_grad[0] * (-dxij) + rho_grad[1] * (-dyij);
				double term3 = 2 * term2 * (mj / rhoj) * (((-dx * normalizedGradientInfluence * dx) + (-dy * normalizedGradientInfluence * dy)) / rSq);

				stateRates[4] += term1 + pparams->delta * h * pparams->cSound * term3;
				//stateRates[4] += mj * normalizedGradientInfluence * dirVel;
				//stateRates[4] += mj * spikyImprovedD(constantSpikyImprovedD, rOh) * (dvxij + dvyij);



				//debug
				int origIndex = pparticles->gridParticleIndex[index];
				if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
					int sendingIndex = pparticles->gridParticleIndex[ind1];
					printf("AP; it=%u, P# %u -> w/ %u, rhoj %f \n", pparams->ind1, origIndex, sendingIndex, rhoj);
				}



				/*
					int originalIndex = pparticles->gridParticleIndex[index];
					if (originalIndex==0)
					printf("host %d vtiy %10.10f term1Y %10.10f \n",originalIndex,vTiy, term1Y);
					*/
					/*
										int originalIndex = pparticles->gridParticleIndex[index];
										int originalInd1 = pparticles->gridParticleIndex[ind1];
										printf("host %d is close to %d %f %f %f \n",originalIndex, originalInd1,pressurei, dist,forces.y);
					*/

			}; //end excluding self
		}; //end looping over the cell
	};  //if start index is not empty

//return forces;
	return;
}