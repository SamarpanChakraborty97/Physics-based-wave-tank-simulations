#ifndef SPH2DCPPCUDA_H
#define SPH2DCPPCUDA_H

// header file for SPH2DCPPCUDA


struct kinematicsFunctionStructure {
	int xFunctionType;
	double p1x;
	double p2x;
	double p3x;
	double p4x;
	double p5x;
	double p6x;
	double p7x;
	double p8x;
	int yFunctionType;
	double p1y;
	double p2y;
	double p3y;
	double p4y;
	int r1;
	int r2;
};

struct paramsType {
	int nTime;            //computational parameters UI
	double dt;         // UI
	bool importDensity;  //set flag to import density
	int storageStride; // UI
	int nTotal;        // CO
	int nConstrained;  // UI
	int nFree;         // UI
	int nFunctions;   //# functions that move constrained particles
	double gravity;   //physical parameters UI
	double nu;        //viscosity UI
	double rRef;      // reference density UI
	double cSound;
	int nCellsX; //binning parameters CO
	int nCellsY;
	int nCellsTotal; // CO
	double globalOriginX; //CO
	double globalOriginY; //CO
	double cellSizeRecip; //CO
	double mass;
	double h; //kernel parameters UI
	double h2; //CO
	int ind1;  //current iteration
	int nFreeParticleForceDomain; //# particles acted on by SPH forces (as opposed to Leonard Jones type forces)
	double viscoAlpha;
	double viscoBeta;
	double constDensity;
	double spikyImprovedD;
	double constwendland;
	double constwendlandD;
	double delta;
	double epsilon;
	double beta;
	double relaxStart;
	double relaxEnd;
	double betaDis;
	int DEBUGinfo;  //display debug info
	int DEBUGpNum;  //debug particle #
};

struct particleStructure {
	double* x;
	double* y;
	double* vx;
	double* vy;
	double* fx;
	double* fy;
	double* density;
	double* ocx;
	double* ocy;
	double* previousX;
	double* previousY;
	int* gridParticleHash;
	int* gridParticleIndex;
	int* cellStart;
	int* cellEnd;
	double* sortedX;
	double* sortedY;
	double* sortedVx;
	double* sortedVy;
	double* sortedRho;
	double* sortedA11;
	double* sortedA12;
	double* sortedA22;
	double* XSPHVelX;
	double* XSPHVelY;
	double* vxH;
	double* vyH;
	int* color;
	double* sortedPressure;
	double* pressure;
	double* sortedRhoFiltered;
	double* sorteddRhodt;
	double* timeSteps;
	double* det_values;
};


#endif
