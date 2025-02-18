#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <iostream>
#include <stdio.h>

__global__ void initializeDensity2(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	double density = 1000;
	pparticles->sortedRho[index] = density;   //this one is not used

	//	printf("sorted Rho of a particle %u is %f\n",index, density);
	int originalIndex = pparticles->gridParticleIndex[index];
	pparticles->density[originalIndex] = density;   //this one is used

	pparticles->sorteddRhodt[index] = 0;

	return;
}