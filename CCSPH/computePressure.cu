#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <stdio.h>
#include "smoothingKernels.cuh"

__global__ void computePressure(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	//	if (colorSorted[index]==0) return;  //its a boundary particle
	double rhoRef = pparams->rRef;
	double cSquared = pparams->cSound * pparams->cSound;

	//read primary particle data - this is sorted data
	//double posXi = pparticles->sortedX[index];
	//double posYi = pparticles->sortedY[index];
	//double velXi = pparticles->sortedVx[index];
	//double velYi = pparticles->sortedVy[index];
	double rhoi = pparticles->sortedRho[index];
	double pressurei = computePressure(rhoi, rhoRef, cSquared);

	pparticles->sortedPressure[index] = pressurei;
	int originalIndex = pparticles->gridParticleIndex[index];
	pparticles->pressure[originalIndex] = pressurei;   //this one is used
}