#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"

__global__ void updatePositionFreeParticles(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < (*pparams).nFree)
	{  //only operate over free particles

		double dt = pparams->dt;
		double x0 = pparams->relaxStart;
		double x1 = pparams->relaxEnd;
		double beta = pparams->betaDis;

		//int sortedIndex = pparticles->gridParticleIndex[index];
		//1st order Euler
#if 0
		pparticles->x[index] += pparticles->vx[index] * pparams->dt;
		pparticles->y[index] += pparticles->vy[index] * pparams->dt;
#endif


		//Leapfrog
#if 1
		//pparticles->sortedX[index] += pparticles->vxH[sortedIndex] * dt;
		//pparticles->sortedY[index] += pparticles->vyH[sortedIndex] * dt;

		pparticles->x[index] += pparticles->vxH[index] * dt;
		pparticles->y[index] += pparticles->vyH[index] * dt;

		if ((pparticles->x[index] >= x0) && (pparticles->x[index] <= x1)) {
			double f = 1 - dt * beta * ((pparticles->x[index] - x0) / (x1 - x0)) * ((pparticles->x[index] - x0) / (x1 - x0));
			pparticles->vxH[index] = pparticles->vxH[index] * f;
			pparticles->vyH[index] = pparticles->vyH[index] * f;
		}

		//pparticles->x[index] += pparticles->vx[index] * dt;
		//pparticles->y[index] += pparticles->vy[index] * dt;
#endif
		/*
		if ((pparticles->x[index] > relS) && (pparticles->x[index] < relE))
		{
			double dampV = b * dt * ((pparticles->x[index] - relS) / (relE - pparticles->x[index])) * ((pparticles->x[index] - relS) / (relE - pparticles->x[index]));
			pparticles->vx[index] = pparticles->vx[index] * (1 - dampV);
			pparticles->vy[index] = pparticles->vy[index] * (1 - dampV);
		}
		*/

		//debug
		if ((index == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
			printf("POS; it=%u, P# %u; px %f, py %f \n", pparams->ind1, index, pparticles->x[index], pparticles->y[index]);
		}



	} //end looping over free

	return;
}