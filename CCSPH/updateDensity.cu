#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"
#include "smoothingKernels.cuh"

__global__ void updateDensity(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	double rhoRef = pparams->rRef;
	double cSquared = pparams->cSound * pparams->cSound;

	if (index < (*pparams).nTotal) {  //operate over all particles

		double dt = pparams->dt;

		//store the updated sorted density
		pparticles->sortedRho[index] += pparticles->sorteddRhodt[index] * dt;
		pparticles->sortedPressure[index] = computePressure(pparticles->sortedRho[index], rhoRef, cSquared);
		//printf("%d \n", pparticles->sortedRho[index]);

#if 0  //limit the density 
		if (> ) {
		}

#endif

	}

	return;
}