#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"

__global__ void copySortedRhoToDensity(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= (*pparams).nTotal) return;

	int origIndex = pparticles->gridParticleIndex[index];
	pparticles->density[origIndex] = pparticles->sortedRho[index];
	pparticles->pressure[origIndex] = pparticles->sortedPressure[index];

	//store density in an "unsorted form" the same order as x, y, vx, vy
	//at the begining of each step, new hash will be computed, and x, y, vx, vy density will be placed in their sorted arrays
	//therefore an array of unsorted denisty must exist




	return;
}