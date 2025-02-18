#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"

__global__ void DEBUGdisplayMoveDensity(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < (*pparams).nTotal) {  //operate over all particles

	//DEBUG
		int origIndex = pparticles->gridParticleIndex[index];
		if ((origIndex == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
			printf("RHO filt; it=%u, P# %u; rhoi %f\n", pparams->ind1, origIndex, pparticles->sortedRho[index]);
		}


	}
	return;
}