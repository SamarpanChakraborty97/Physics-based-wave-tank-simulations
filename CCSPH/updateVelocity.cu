#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <iostream>
#include <stdio.h>

__global__ void updateVelocity(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < (*pparams).nFree) {  //only operate over free particles
		//the accelerations are stored; simply add gravity to the y-dir
		//and incorporate the XSPH terms

	//F = m a;
	//a = F/m

//Euler 1st order
#if 0
		double accelX = pparticles->fx[index]; //fx actually stores an acceleration; no need to divide by mass
		double accelY = pparticles->fy[index] + pparams->gravity; //need to add gravity

		double dt = pparams->dt;
		double vNewx = pparticles->vx[index] + accelX * dt;  //XSPH is incorporated previously
		double vNewy = pparticles->vy[index] + accelY * dt;  //
		//store the updated velocity
		pparticles->vx[index] = vNewx;
		pparticles->vy[index] = vNewy;
#endif

		//Leapfrog
#if 1
			//UPDATE THE VELOCITY WITH XSPH VELOCITY
		double dt = pparams->dt;
		if (pparams->ind1 == 0) {
			//pparticles->vxH[index] = pparticles->vx[index] + pparticles->fx[index] * dt / 2;
			//pparticles->vyH[index] = pparticles->vy[index] + (pparticles->fy[index] + pparams->gravity) * dt / 2;

			pparticles->vxH[index] += pparticles->fx[index] * dt / 2;
			pparticles->vyH[index] += (pparticles->fy[index] + pparams->gravity) * dt / 2;
			pparticles->vx[index] = pparticles->vxH[index] + pparticles->fx[index] * dt / 2;
			pparticles->vy[index] = pparticles->vyH[index] + (pparticles->fy[index] + pparams->gravity) * dt / 2;
			//	printf("tStep==0\n");
			//	printf("vy of particle 1 %f\n",pparticles->vy[1]);
		}
		else {
			pparticles->vxH[index] += pparticles->fx[index] * dt;
			pparticles->vyH[index] += (pparticles->fy[index] + pparams->gravity) * dt;
			pparticles->vx[index] = pparticles->vxH[index] + pparticles->fx[index] * dt / 2;
			pparticles->vy[index] = pparticles->vyH[index] + (pparticles->fy[index] + pparams->gravity) * dt / 2;
		}
#endif








		//debug
		if ((index == pparams->DEBUGpNum) && (pparams->DEBUGinfo == 1)) {
			printf("V+update; it=%u, P# %u; vx %f, vy %f \n", pparams->ind1, pparams->DEBUGpNum, pparticles->vx[index], pparticles->vy[index]);
		}


	}
	return;
}