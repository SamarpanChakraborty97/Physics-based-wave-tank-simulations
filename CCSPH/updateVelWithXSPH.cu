#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"


__global__ void updateVelWithXSPH(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < (*pparams).nFree) {  //operate over free particles


	//DEBUG
		if ((index == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
			printf("V before XSPH; it=%u, P# %u; vx %f, vy %f \n", pparams->ind1, index, pparticles->vx[index], pparticles->vy[index]);
		}

		//correct the velcoity with the XSPH correction
		//pparticles->vx[index] += pparticles->XSPHVelX[index];
		//pparticles->vy[index] += pparticles->XSPHVelY[index];

		pparticles->vxH[index] += pparticles->XSPHVelX[index];
		pparticles->vyH[index] += pparticles->XSPHVelY[index];

		//DEBUG
		if ((index == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
			printf("V+XSPH; it=%u, P# %u; vx %f, vy %f \n", pparams->ind1, index, pparticles->vx[index], pparticles->vy[index]);
		}



#if 0  //limit the velocity
		if (> ) {
		}

#endif

	}


	return;
}