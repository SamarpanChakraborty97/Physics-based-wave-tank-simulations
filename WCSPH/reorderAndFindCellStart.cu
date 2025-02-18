#include <cuda.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "thrust/device_vector.h"


__global__ void reorderAndFindCellStart(struct particleStructure* pparticles, struct paramsType* pparams)
{
	extern __shared__ int sharedHash[];  //the size of this array is passed as the third arguement; it defaults to the blockDim.x; which won't work for us
	int index = blockIdx.x * blockDim.x + threadIdx.x;  //each thread treats a particle

	int hash;
	//# particle not multiple block size
	if (index < (*pparams).nTotal) {
		hash = pparticles->gridParticleHash[index]; //work on valid data; get 1D hash value

		sharedHash[threadIdx.x + 1] = hash;

		if (index > 0 && threadIdx.x == 0) {  //the first one in the block - get it from the previous block
			sharedHash[0] = pparticles->gridParticleHash[index - 1];
		}
	}

	__syncthreads();

	if (index < (*pparams).nTotal) {

		if (index == 0 || hash != sharedHash[threadIdx.x])  //only one thread/cell will satisfy this
		{
			pparticles->cellStart[hash] = index;  //store as a function of hash value
			if (index > 0) { pparticles->cellEnd[sharedHash[threadIdx.x]] = index; };  //rem sharedHash blockSize+1
		}

		if (index == (*pparams).nTotal - 1)
		{
			pparticles->cellEnd[hash] = index + 1;
		};

		int sortedIndex = pparticles->gridParticleIndex[index];  //this was sorted by thrust
		//double dt = pparams->dt;
		//double b = pparams->beta;
		//double relS = pparams->relaxStart;
		//double relE = pparams->relaxEnd;

		//cellStart and cellEnd contain the sorted particle indices which span i.e. begin and end 
		//each hash cell

		//not all hash cells will contain particles, so they are assigned a value of 0xffffffff
		//all particles will be associated with a cell.


		//sort arrays
		//sorted arrays are accessible by cellStart and cellEnd
		//for instance, assuming cell 5 is populated with a few particles
		//cellStart[5] & cellEnd[5] provide an indices which span those particles
		//the particles must be sorted so consecutive particles from cellStart[5] to cellEnd[5] are contained
		//within the cell


		//all individual properties of particles must be sorted
		//there may be a faster way to sort these arrays
		pparticles->sortedX[index] = pparticles->x[sortedIndex];
		pparticles->sortedY[index] = pparticles->y[sortedIndex];
		pparticles->sortedVx[index] = pparticles->vx[sortedIndex];
		pparticles->sortedVy[index] = pparticles->vy[sortedIndex];
		pparticles->sortedRho[index] = pparticles->density[sortedIndex];
		//pparticles->sortedPressure[index] = pparticles->pressure[sortedIndex];

		/*
		if ((pparticles->sortedX[index] > relS) && (pparticles->sortedX[index] < relE))
		{
			double dampV = b * dt * ((pparticles->sortedX[index] - relS) / (relE - pparticles->sortedX[index])) * ((pparticles->sortedX[index] - relS) / (relE - pparticles->sortedX[index]));
			pparticles->sortedVx[index] = pparticles->sortedVx[index] * (1 - dampV);
			pparticles->sortedVy[index] = pparticles->sortedVy[index] * (1 - dampV);
		}
		*/
		//if each particle had a unique smoothing length and mass, they would have to be sorted as sorted arrays as well.
		//pM,  pR
	};

	//printf("%d \n", pparticles->sortedX[index]);
//	if (index < 1){
//		printf("cell start[0] %d\n",dCellStart[0]);
//		printf("cell end[0] %d\n",dCellEnd[0]);
//	};




	return;
}