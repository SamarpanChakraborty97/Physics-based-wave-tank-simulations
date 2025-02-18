#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "stdio.h"

__global__ void updateShift(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < (*pparams).nFree)
	{  //only operate over free particles

		double dt = pparams->dt;
		//double b = pparams->beta;
		//double relS = pparams->relaxStart;
		//double relE = pparams->relaxEnd;

		double Afst = 1.5;
		double Afsm = 2;
		double A = 2;
		double velX = pparticles->vx[index];
		double velY = pparticles->vy[index];
		double velMag = sqrt(velX * velX + velY * velY);

		if (pparticles->posDiv[index] - Afst < 0) {
			double Afsc = (pparticles->posDiv[index] - Afst) / (Afsm - Afst);
			double drX = -Afsc * A * pparams->h * velMag * dt * pparticles->shiftGradX[index];
			double drY = -Afsc * A * pparams->h * velMag * dt * pparticles->shiftGradY[index];
			pparticles->x[index] += pparticles->vxH[index] * dt + drX;
			pparticles->y[index] += pparticles->vyH[index] * dt + drY;
		}
		else if (pparticles->posDiv[index] - Afst == 0) {
			double drX = -A * pparams->h * velMag * dt * pparticles->shiftGradX[index];
			double drY = -A * pparams->h * velMag * dt * pparticles->shiftGradY[index];
			pparticles->x[index] += pparticles->vxH[index] * dt + drX;
			pparticles->y[index] += pparticles->vyH[index] * dt + drY;
		}
		else {
			pparticles->x[index] += pparticles->vxH[index] * dt;
			pparticles->y[index] += pparticles->vyH[index] * dt;
		}
	}

	return;
}