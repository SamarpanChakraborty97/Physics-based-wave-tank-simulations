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

//UI - user input
//CO  - computed once
struct paramsType {
	int nTime;            //computational parameters UI
	double dt;         // UI
	bool importDensity;  //set flag to import density
	int storageStride; // UI
	int nTotal;        // CO
	int nConstrained;  // UI
	int nFree;         // UI
	int nMeasured;     // UI
	int nFunctions;   //# functions that move constrained particles
	double gravity;   //physical parameters UI
	double vf;        // max velocity of fluid - used in  Tait's eq UI
	double nu;        //viscosity UI
	double rRef;      // reference density UI
	double tenVMaxSq; // used in Tait's equation B = (10vMax)^2/(rho*gamma)
	int nCellsX; //binning parameters CO
	int nCellsY;
	//int nCellsZ;//                 3D
	int nCellsTotal; // CO
	double globalOriginX; //CO
	double globalOriginY; //CO
	//double globalOrginZ; //CO
	double cellSizeRecip; //CO
	double h; //kernel parameters UI
	double h2; //CO
	double h8; //CO
	int ind1;  //current iteration
	int nFreeParticleForceDomain; //# particles acted on by SPH forces (as opposed to Leonard Jones type forces)
	double quartic;
	double quarticD;
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
	int DEBUGinfo;  //display debug info
	int DEBUGpNum;  //debug particle #
};

struct particleStructure {
	double* x;
	double* y;
	//double* z;
	double* vx;
	double* vy;
	double* fx;
	double* fy;
	double* radius;
	double* mass;
	double* density;
	double* L1;
	double* L2;
	double* L3;
	double* L4;
	double* rhoGradX;
	double* rhoGradY;
	double* posDiv;
	double* shiftGradX;
	double* shiftGradY;
	double* ocx;
	double* ocy;
	//double* omx;
	//double* omy;
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
	double* sortedShift;
	double* XSPHVelX;
	double* XSPHVelY;
	double* vxH;
	double* vyH;
	int* color;
	double* sortedPressure;
	double* pressure;
	double* sortedRhoFiltered;
	double* sorteddRhodt;
	//double* phaseValues1;
	//double* phaseValues2;
	//double* phaseValues3;
	//double* phaseValues4;
	//double* phaseValues5;
	//double* phaseValues6;
	double* Ai;
	double* Bi;
	double* fi;
	double* ki;
};


#endif
