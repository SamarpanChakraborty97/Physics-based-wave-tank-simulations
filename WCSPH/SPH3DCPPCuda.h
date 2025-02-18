#ifndef SPH3DCPPCUDA_H
#define SPH3DCPPCUDA_H

// header file for SPH3DCPPCUDA


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
	int zFunctionType; //This is done to incorporate the dam-break in 3-D.
	double p1z;        //parameter1 for Dam-Break in 3-D       
	double p2z;        //parameter2 for Dam-Break in 3-D  
	int r1;
	int r2;
};

//UI - user input
//CO  - computed once
struct paramsType {
	int nTime;         //computational parameters UI
	double dt;         // UI
	bool importDensity;  //set flag to import density
	int storageStride; // UI
	int nTotal;        // CO
	int nConstrained;  // UI
	int nFree;         // UI
	int nFunctions;    //# functions that move constrained particles
	double gravity;    //physical parameters UI
	double vf;         // max velocity of fluid - used in  Tait's eq UI
	double nu;         //viscosity UI
	double rRef;       // reference density UI
	double tenVMaxSq;  // used in Tait's equation B = (10vMax)^2/(rho*gamma)
	int nCellsX;       //binning parameters CO
	int nCellsY;       //                   CO
	int nCellsZ;                                                                 //3D
	int nCellsTotal; // CO
	double globalOriginX; //CO
	double globalOriginY; //CO
	double globalOriginZ;                                                        //3D
	double cellSizeRecip; //CO
	double h; //kernel parameters UI
	double h2; //CO
	double h3; //CO
	//double constDensity; //CO
	//double constViscosity; //CO
	//double constPressure; //CO
	int ind1;  //current iteration
	int nFreeParticleForceDomain; //# particles acted on by SPH forces (as opposed to Leonard Jones type forces)
	//double quartic;
	//double quarticD;
	double viscoAlpha;
	double viscoBeta;
	double delta;
	//double spikyImprovedD;
	double constwendland;
	double constwendlandD;
	double epsilon;
	int DEBUGinfo;  //display debug info
	int DEBUGpNum;  //debug particle #
};

struct particleStructure {
	double* x;
	double* y;
	double* z; // 3D
	double* vx;
	double* vy;
	double* vz; //3D
	double* fx;
	double* fy;
	double* fz; //3D
	double* radius;
	double* mass;
	double* density;
	double* ocx;
	double* ocy;
	double* ocz; //3D
	double* previousX;
	double* previousY;
	double* previousZ; //3D
	int* gridParticleHash;
	int* gridParticleIndex;
	int* cellStart;
	int* cellEnd;
	double* sortedX;
	double* sortedY;
	double* sortedZ; //3D
	double* sortedVx;
	double* sortedVy;
	double* sortedVz; //3D
	double* sortedRho;
	double* XSPHVelX;
	double* XSPHVelY;
	double* XSPHVelZ; //3D
	double* vxH;
	double* vyH;
	double* vzH; //3D
	int* color;
	double* sortedPressure;
	double* pressure;
	double* sortedRhoFiltered;
	double* sorteddRhodt;
};


#endif


