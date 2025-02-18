//the cuda particles example, included in the SDK, was used as a reference
//the algorithms presented here were developed for this particular implmentation


#include "stdafx.h"
#include <stdio.h>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "exportSystem.h"
#include "moveDataToGPU.h"
#include <omp.h>
//#include <iostream>


//note: had to modify C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v6.5\bin\nvcc.profile to get open mp working
//added "/openmp" with quotes in INCLUDE statement


#include "thrust/sort.h"
#include "thrust/device_vector.h"

#include <time.h>
#include "SPH2DCUDA.h"

#define threadsPerBlock 128

#define timeComponents (0)
#define DEBUGING (0)

__global__ void moveConstrainedParticlesD(struct particleStructure* pdParticles, struct paramsType* dParams, kinematicsFunctionStructure* dKinematicsFunction);
__global__ void moveConstrainedParticlesD2(struct particleStructure* pdParticles, struct paramsType* dParams, kinematicsFunctionStructure* dKinematicsFunction);
__global__ void calcHash(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void reorderAndFindCellStart(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void initializeDensity(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void initializeDensity2(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computedVdt(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computePressure(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void updateVelocity(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computeXSPHInfluence(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void updatePositionFreeParticles(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computeShift(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computeShiftGrad(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computeDiv(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void updateShift(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void computedRhodt(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void reinitializeDensity(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void moveFilteredDensityToDensity(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void updateDensity(struct particleStructure* pdParticles, struct paramsType* dParams);
__global__ void copySortedRhoToDensity(struct particleStructure* pparticles, struct paramsType* pparams);
__global__ void updateVelWithXSPH(struct particleStructure* pdParticles, struct paramsType* pdParams);

//DEBUG
__global__ void DEBUGdisplayDensity(struct particleStructure* pdParticles, struct paramsType* pdParams);
__global__ void DEBUGdisplaydRhodt(struct particleStructure* pdParticles, struct paramsType* pdParams);
__global__ void DEBUGdisplayMoveDensity(struct particleStructure* pdParticles, struct paramsType* pdParams);


void SPH2DCuda(particleStructure* particles, paramsType* params, std::vector<kinematicsFunctionStructure>* kinematicsFunction, std::string outputDir)
{


#if DEBUGING
#define DEBUG(...) \
	do{\
	__VA_ARGS__;\
	} while(0)
#else
#define DEBUG(...)\
	do {;} while(0)
#endif


#if timeComponents
	cudaEvent_t start;
	cudaEvent_t stop;
	float cudaTimeTemp;
	const int nProcedures = 21;
	int nChars = 30;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float cudaTime[nProcedures];
	char* procedureNames[nProcedures];
	for (int ind1 = 0; ind1 < nProcedures; ind1++) {
		cudaTime[ind1] = 0;  //initialize cumulative timers
		procedureNames[ind1] = (char*)malloc(sizeof(char) * nChars);
		procedureNames[ind1] = ".";
	}
	procedureNames[0] = "calcHash                     ";
	procedureNames[1] = "thrust::sort_by_key..........";
	procedureNames[2] = "reorderAndFindCellStart......";
	procedureNames[3] = "computedVdt                  ";
	procedureNames[4] = "updateVelWithXSPH............";
	procedureNames[5] = "updateVelocity               ";
	procedureNames[6] = "updateDensity................";
	procedureNames[7] = "reinitializeDensity          ";
	procedureNames[8] = "moveFilteredDensityToDensity.";
	procedureNames[9] = "copySortedRhoToDensity       ";
	procedureNames[10] = "updatePositionFreeParticles..";
	procedureNames[11] = "moveConstrainedParticlesD    ";
	procedureNames[12] = "cudaMemcpy...................";
	procedureNames[13] = "exportSystem***OMP***        ";
#define TIME(processNumber,...) \
	do {\
	cudaEventRecord(start,0); \
	__VA_ARGS__; \
	cudaEventRecord(stop,0); \
	cudaEventSynchronize(stop);	\
	cudaEventElapsedTime(&cudaTimeTemp, start, stop); \
	cudaTime[processNumber] += cudaTimeTemp; \
		} while(0)
#else
#define TIME(processNumber,...) __VA_ARGS__;
#endif

	//heartbeat
	cudaEvent_t tStart, tStop;
	cudaEventCreate(&tStart);
	cudaEventCreate(&tStop);
	cudaEventRecord(tStart);
	float elapsedTime;
	std::cout << '\n' << omp_get_max_threads();

	omp_set_num_threads(2);  //one thread writes the files, one thread issues GPU commands, more threads could only cause problems


	int smemSize = (threadsPerBlock + 1) * sizeof(int);
	int nTotal = params->nTotal;  //number of particles
	std::cout << '\n' << params->nTotal;
	int numberBlocks = (nTotal % threadsPerBlock != 0) ? (nTotal / threadsPerBlock + 1) : (nTotal / threadsPerBlock); //if condt ? outTrue : outFalse

	//int nFree = params->nFree;//number of free particles;
	int numberBlocksFree = (params->nFree % threadsPerBlock != 0) ? (nTotal / threadsPerBlock + 1) : (nTotal / threadsPerBlock); //if condt ? outTrue : outFalse

	int nC = params->nConstrained;  //number of constrained particles
	int numberBlocksConstrained = (nC % threadsPerBlock != 0) ? (nC / threadsPerBlock + 1) : (nC / threadsPerBlock); //if condt ? outTrue : outFalse

	//int nM = params->nMeasured; //number of constrained particles
	//int numberBlocksMeasured = (nM % threadsPerBlock != 0) ? (nM / threadsPerBlock + 1) : (nM / threadsPerBlock); //if condt ? outTrue : outFalse

	int nT = params->nTime; //number of time steps
	int storageStride = params->storageStride;
	int nStrides = nT / storageStride; //integer division implies floor();
	//printf("%d", nStrides);


	//define the # forces acting on the free particles 
	//in standard simulation, boundary particles administer forces as free particles
	//in Leonard-Jones type forces from boundary particles are handeled differently
	params->nFreeParticleForceDomain = params->nTotal;


	//--determine bin size--
	//bin size in x and y dimensions
	//its the maximum of the smoothing radii
	double binSize = 2 * particles->radius[0];

	//store the bin size and ...
	params->cellSizeRecip = 1 / binSize;

	//--determine domain limits--
	//domain limits are taken as the maximum and minimums in x & y directions from the INITIAL data set
	//if the user would like to specify larger domain limits, they can input 
	//constrained particles at specified locations.

	double2 globalMin = { particles->x[0],particles->y[0] };
	double2 globalMax = { particles->x[0],particles->y[0] };

	for (int ind1 = 0; ind1 < nTotal; ind1++) {
		if (particles->x[ind1] > globalMax.x) { globalMax.x = particles->x[ind1]; };
		if (particles->x[ind1] < globalMin.x) { globalMin.x = particles->x[ind1]; };
		if (particles->y[ind1] > globalMax.y) { globalMax.y = particles->y[ind1]; };
		if (particles->y[ind1] < globalMin.y) { globalMin.y = particles->y[ind1]; };
	}
	params->globalOriginX = globalMin.x;
	params->globalOriginY = globalMin.y;

	double2 delta = { globalMax.x - globalMin.x,globalMax.y - globalMin.y };
	//params->nCellsY     = fmod(delta.y,binSize) !=0 ? (int) (delta.y/binSize+1) : (int) (delta.y/binSize); //if condt ? outTrue : outFalse
	//params->nCellsX     = fmod(delta.x,binSize) !=0 ? (int) (delta.x/binSize+1) : (int) (delta.x/binSize); //if condt ? outTrue : outFalse
	params->nCellsY = (int)(delta.y / binSize) + 1; //cast to int is like floor, if its on the upper bin edge, then an additional bin survives
	params->nCellsX = (int)(delta.x / binSize) + 1;
	//std::cout << '\n' << params->nCellsTotal;
	params->nCellsTotal = (params->nCellsX) * (params->nCellsY);
	std::cout << '\n' << params->nCellsTotal;
	params->nFunctions = (*kinematicsFunction).size(); //store the number of kinematics functions
	//need the +1 becuase a particle might be on the max bin edge


	//precompute some constants used in the kernels, derived in Mathematica, and in Bindell.  Muller in 2D
	//the kernels as used are:
	//density   - 4/(pi h^8)*(h^2-r^2)^3
	//pressure  - -30/(h^5*pi)*(h-r)^2
	//viscosity - 40/(h^5*pi)*(h-r)


	double h = particles->radius[0];  //make particles 2x large here and above; use 2x smoothing length; revert in viscosity muij
	double h2 = h * h;
	double h8 = h2 * h2 * h2 * h2;
	double const1 = 4.0 / (3.141592654 * h2);  //poly6 2D              - verified
	double constPressure = (30.0 / (3.1415926535 * h2 * h2 * h));
	double constViscosity = 40.0 / (3.1415926535 * h2 * h2 * h);
	double quartic = 1 / (7 * 3.141592654 * h2);      //Liu normalized to h=1;
	double quarticD = -60 / (7 * 3.141592654 * h2 * h);  //Liu derivative
	double spikyImprovedD = -12 / (3.141592654 * h * h * h);
	double constWend = 7 / (4 * 3.1415926535 * h2);
	double constWendD = -35 / (4 * 3.1415926535 * h2 * h);

	params->h = h;
	params->h2 = h2;
	params->h8 = h8;
	params->constwendland = constWend;
	params->constwendlandD = constWendD;
	params->tenVMaxSq = (15 * params->vf) * (15 * params->vf);
	params->ind1 = 0;
	params->quartic = quartic;
	params->quarticD = quarticD;
	params->spikyImprovedD = spikyImprovedD;
	params->constwendland = constWend;
	params->constwendlandD = constWendD;

	//DEBUG related
	params->DEBUGinfo = 0;  //display debug info
	params->DEBUGpNum = 7171;  //debug particle #

	//make device pointers
	struct particleStructure* pdParticles = 0;  //device particles; 
	struct particleStructure* pdParticlesHostMirror = new struct particleStructure;  //resides on host, contains device pointers
	struct paramsType* pdParams = 0;  //device parameters
	struct kinematicsFunctionStructure* pdKinematicsFunction = 0;  //device kinematics

		//allocate memory and transfer data to GPU

	moveDataToGPU(particles, params, kinematicsFunction, &pdParticles, &pdParams, &pdKinematicsFunction, &pdParticlesHostMirror);
	//first, get constrained particles to their proper location
	//and reset the device time

	//moveConstrainedParticlesD2 increments the GPU copy of ind1
	//DEBUG - Nothing has been sorted, so the first time is no problem
	moveConstrainedParticlesD << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams, pdKinematicsFunction);
	cudaMemcpy(pdParams, params, sizeof(paramsType), cudaMemcpyHostToDevice); //reset the time (ind1)

	//initialize density and pressure
	calcHash << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams);  //compute hash value
	thrust::sort_by_key(thrust::device_ptr<int>(pdParticlesHostMirror->gridParticleHash), thrust::device_ptr<int>(pdParticlesHostMirror->gridParticleHash + (*params).nTotal), thrust::device_ptr<int>(pdParticlesHostMirror->gridParticleIndex));  //sort
	cudaMemset(pdParticlesHostMirror->cellStart, 0xffffffff, (*params).nCellsTotal * sizeof(int));
	cudaDeviceSynchronize();
	reorderAndFindCellStart << <numberBlocks, threadsPerBlock, smemSize >> > (pdParticles, pdParams);
	if (params->importDensity == 0) {  //need to initialize density
		//initializeDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams);     //compute
		initializeDensity2 << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams);     //compute
		printf("\nNot importing density; initializing density\n");
	}
	else { printf("\nImporting density\n"); }

	computePressure << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams);     //compute pressure
	//cudaMemcpy(particles->density,pdParticlesHostMirror->sortedRho,nTotal*sizeof(double),cudaMemcpyDeviceToHost);  // no need
	//transfer memory back to host in anticipation of writing 0th file

	//enter loop
	unsigned int iteration = 0;  //global iteration #
	while (iteration < nT)
	{
		unsigned int snapshotIteration = iteration / storageStride;  //copy the current iteration to use in the exported file name
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////TRANSFER CURRENT DATA TO HOST IN ANTICIPATION OF FILE EXPORT//////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
		cudaDeviceSynchronize();
		TIME(10, cudaMemcpy(particles->x, pdParticlesHostMirror->x, nTotal * sizeof(double), cudaMemcpyDeviceToHost););
		TIME(10, cudaMemcpy(particles->y, pdParticlesHostMirror->y, nTotal * sizeof(double), cudaMemcpyDeviceToHost););
		//TIME(10, cudaMemcpy(particles->vx, pdParticlesHostMirror->vx, nTotal * sizeof(double), cudaMemcpyDeviceToHost););
		//TIME(10, cudaMemcpy(particles->vy, pdParticlesHostMirror->vy, nTotal * sizeof(double), cudaMemcpyDeviceToHost););
		TIME(10, cudaMemcpy(particles->density, pdParticlesHostMirror->density, nTotal * sizeof(double), cudaMemcpyDeviceToHost););
		//TIME(10, cudaMemcpy(particles->pressure, pdParticlesHostMirror->pressure, nTotal * sizeof(double), cudaMemcpyDeviceToHost););
		//std::cout <<'\n' <<  params->nMeasured;

#if timeComponents																		//
		cudaEventRecord(stop, 0);												//
		cudaEventSynchronize(stop);												//
#endif	

#pragma omp parallel sections
		{
#pragma omp section
			{
				/////////////////////////////////////////////////////////////////////////////////////////////////////					
								//EXPORT DATA from previous transfer while the next iterations are being computed
								//no cuda calls in this section
				/////////////////////////////////////////////////////////////////////////////////////////////////////
#if timeComponents
			//cudaEventRecord(start[13],0); 
				double fileTimeStart = omp_get_wtime();
#endif
				int	output = exportSystem(particles, snapshotIteration, params, outputDir);
#if timeComponents
				//			cudaEventRecord(stop{13],0);
				//			cudaEventSynchronize(stop[13]);
				double fileTimeStop = omp_get_wtime();
				cudaTime[13] = (float)(fileTimeStop - fileTimeStart) * 1000; //file writing time in ms
				//printf("OMP thread number %u\n",omp_get_thread_num());
#endif
			}
#pragma omp section
			{
				/////////////////////////////////////////////////////////////////////////////////////////////////////
				//COMPUTE NEXT STRIDES WORTH OF ITERATIONS
				/////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ind2 = 0; ind2 < storageStride; ind2++)
				{
					TIME(0, calcHash << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););  //compute hash value
					TIME(1, thrust::sort_by_key(thrust::device_ptr<int>(pdParticlesHostMirror->gridParticleHash), thrust::device_ptr<int>(pdParticlesHostMirror->gridParticleHash + (*params).nTotal), thrust::device_ptr<int>(pdParticlesHostMirror->gridParticleIndex)););  //sort
					cudaDeviceSynchronize();
					//dGridParticleHash is now sorted into ascending order
					//dGridParticleIndex is now sorted based on the rearrangement of dGridParticleHash


					/*  Debug
					cudaDeviceSynchronize();
					printf("value of pdParticlesHostMirror->gridParticleHash%p\n",pdParticlesHostMirror->gridParticleHash);

					int* tempInt   = new int[nTotal];

					cudaMemcpy(tempInt,pdParticlesHostMirror->gridParticleHash,nTotal*sizeof(int),cudaMemcpyDeviceToHost);
					cudaMemcpy(particles->x,pdParticlesHostMirror->x,nTotal*sizeof(double),cudaMemcpyDeviceToHost);
					cudaMemcpy(particles->y,pdParticlesHostMirror->y,nTotal*sizeof(double),cudaMemcpyDeviceToHost);

					int cellMax = tempInt[0];
					for (int ind3 = 0;ind3<nTotal;ind3++)
					{if (tempInt[ind3]>cellMax) {cellMax = tempInt[ind3];}
					};
					printf("%d",cellMax);
					double k2 = 2;
					delete tempInt;
					printf("%d\n",ind2);
					*/

					//set all cells to empty
					cudaMemset(pdParticlesHostMirror->cellStart, 0xffffffff, (*params).nCellsTotal * sizeof(int));
					cudaDeviceSynchronize();
					//reorder data and find cell start
					TIME(2, reorderAndFindCellStart << <numberBlocks, threadsPerBlock, smemSize >> > (pdParticles, pdParams););
					/////////////////COMPUTE DV/DT////////////////////////////////////////////////////////////
					TIME(3, computedVdt << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););		//
					//TIME(4, computedRhodt << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					//////////////////////////////////////////////////////////////////////////////////////////


					////////////////MOVE PARTICLES, UPDATE VELOCITY, UPDATE DENISTY////////////////////
					TIME(4, updateVelWithXSPH << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					TIME(5, updateVelocity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					TIME(6, updateDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					//DEBUG(DEBUGdisplayDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););					   //
					/*
					if (iteration % 20 == 0)  //this is tied to "iteration" not the subiterator counting iterations per stride
					{
						TIME(7, reinitializeDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
						TIME(8, moveFilteredDensityToDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
						DEBUG(DEBUGdisplayMoveDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););				   //
					}
					*/
					TIME(7, copySortedRhoToDensity << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					TIME(8, updatePositionFreeParticles << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););  //have to re-bin in order to accurately reinitialize
					//TIME(11, computeShift << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					//TIME(12, computeShiftGrad << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					//TIME(12, computeDiv << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					//TIME(13, updateShift << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams););
					//move the boundary
					//we only need to cover the constrained particles
					if (numberBlocksConstrained > 0) {
						TIME(9, moveConstrainedParticlesD << <numberBlocks, threadsPerBlock >> > (pdParticles, pdParams, pdKinematicsFunction););
					}  //close constrained particles

		////////////////////////////////////////////////////////////////////////////////////////
		//HEART BEAT
		////////////////////////////////////////////////////////////////////////////////////////
					if (ind2 == (storageStride - 1))
					{
						cudaEventRecord(tStop);
						cudaEventSynchronize(tStop);
						cudaEventElapsedTime(&elapsedTime, tStart, tStop);
						printf("%5.0d steps took %f ms; %f ms/step;OMP thread # % u;\n", storageStride, elapsedTime, elapsedTime / storageStride, omp_get_thread_num());
						//OMP thread # % u
						//omp_get_thread_num()
						cudaEventRecord(tStart);

					}  //status update

					iteration++;  //increment global iteration count
					cudaDeviceSynchronize();


				} // close iterating over a stride
			} // close omp section of iterating over a stride

//////////////////////////////////////////////////////////////////////////////////////////////////////////


		} //close OMP SECTIONS; IMPLICIT BARRIER FOR OMP THREADS


#if timeComponents

		double componentSum = 0;
		//compute elapsed times
		for (int ind3 = 0; ind3 < 13; ind3++) { //exclude file export
			componentSum += cudaTime[ind3];
		}

		printf("\n");
		for (int ind3 = 0; ind3 < 19; ind3++) {
			printf("Time for %s is %f ms, %f, %%\n", procedureNames[ind3], cudaTime[ind3], cudaTime[ind3] / componentSum * 100);
			cudaTime[ind3] = 0;
		}
		printf("componentSum                           is %f\n", componentSum);

#endif



	}  //check while loop to see if more iterations necessary

			//free the GPU memory
			//it is probably freed upon exit



	/*
			cudaFree(dPos);
			cudaFree(dmassRadius);
			cudaFree(dpRho);
			cudaFree(dVel);
			cudaFree(dVelHalf);
			cudaFree(dForce);
			cudaFree(dpColor);
			cudaFree(dpColorSorted);
			cudaFree(dParams);

			cudaFree(dSortedPos);
			cudaFree(dSortedVel);

			cudaFree(dsortedpRho);
			cudaFree(dXSPHVel);
			cudaFree(dGridParticleHash);
			cudaFree(dGridParticleIndex);
			cudaFree(dCellStart);
			cudaFree(dCellEnd);
			cudaFree(dInd1);

			*/

			//end time loop

	return;
}