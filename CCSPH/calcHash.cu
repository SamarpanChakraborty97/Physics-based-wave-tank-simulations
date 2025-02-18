#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include <stdio.h>


__global__ void calcHash(struct particleStructure* pparticles, struct paramsType* pparams) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < (pparams->nTotal)) {

		//compute grid position & clamp
		//remember cells go from [0,nCellsX-1] to give nCellsX # of cells
		int tempX = floor((pparticles->x[index] - (*pparams).globalOriginX) * (*pparams).cellSizeRecip);
		if (tempX >= (*pparams).nCellsX) { 
			tempX = (*pparams).nCellsX - 1; 
			printf("x-exceeded"); 
			
		}  //exceeded original x domain
		if (tempX < 0) { tempX = 0; 
			printf("x-negative"); 
		}  //lower than original x domain

		int tempY = floor((pparticles->y[index] - (*pparams).globalOriginY) * (*pparams).cellSizeRecip);
		if (tempY >= (*pparams).nCellsY) { 
			tempY = (*pparams).nCellsY - 1; 
			printf("y-exceeded"); 
		}  //exceeded original y domain
		if (tempY < 0) { 
			tempY = 0; 
			printf("y-negative"); 
		}  //lower than original y domain

		//compute 1D hash value
		int hash = tempY * (*pparams).nCellsX + tempX;


		//won't happen
		if (hash < 0) { printf("hash<0 particle # %d\n", index); }

		if (hash >= (*pparams).nCellsTotal) { printf("hash>=max particle # %d\n", index); }  //this won't happen

		//store 
		pparticles->gridParticleHash[index] = hash;
		pparticles->gridParticleIndex[index] = index; //

		//dGridParticleIndex is just an array of consecutive indices right now, but it will be sorted later
	}



	return;
}