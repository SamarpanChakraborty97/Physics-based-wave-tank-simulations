#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"
#include "smoothingKernels.cuh"

__global__ void moveFilteredDensityToDensity(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	double rhoRef = pparams->rRef;
	double tenVMaxSq = pparams->tenVMaxSq;
	if (index >= (*pparams).nTotal) return;

	//put filtered density back into "working" density 
	pparticles->sortedRho[index] = pparticles->sortedRhoFiltered[index];
	pparticles->sortedPressure[index] = computePressure(pparticles->sortedRho[index], rhoRef, tenVMaxSq);


	//	pparticles->density[index] = pparticles->unsortedRhoFiltered[index];


	return;
}