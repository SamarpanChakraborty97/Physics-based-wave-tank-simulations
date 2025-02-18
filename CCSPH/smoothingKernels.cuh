#ifndef SMOOTHINGKERNELS_H
#define SMOOTHINGKERNELS_H

inline __device__ double poly6(double cd, double r);
inline __device__ double spikyImprovedD(double cd, double r);
inline __device__ double computePressure(double rhoi, double rhoRef, double tenVMaxSq);

inline __device__ double wendland(double cd, double rOh) {
	//kernel is only valid for [0 <= r <= h]
	//cd -  kernel coefficient 
	//rOh - r/h 
	double rOhSq = (1 - (rOh/2)) * (1 - (rOh/2)); //  (r^2/h^2)
	double W = cd * (1 + 2*rOh) * rOhSq * rOhSq;  //kernel influence
	return W;
}
inline __device__ double wendlandD(double cd, double rOh) {
	//kernel is only valid for [0 <= r <= h]
	//cd -  kernel coefficient
	//rOh - r/h
	double rOhCu = (1 - (rOh / 2)) * (1 - (rOh / 2)) * (1 - (rOh / 2));
	double W = cd * (rOh * rOhCu);  //kernel influence
	return W;
}

inline __device__ double computePressure(double rho, double rhoRef, double tenVMaxSq) {
	//Equation of State
	double temp = rho / rhoRef;
	double tempSq = temp * temp;
	double temp7 = tempSq * tempSq * tempSq * temp;
	double B = tenVMaxSq * rhoRef / 7;
	double pres = B * (temp7 - 1);

	return pres;
}
/*
inline __device__ double computePressure(double rho, double rhoRef, double tenVMaxSq) {
	//Equation of State
	//double temp = rho / rhoRef;
	//double tempSq = temp * temp;
	//double temp7 = tempSq * tempSq * tempSq * temp;
	//double B = tenVMaxSq * rhoRef / 7;
	double pres = tenVMaxSq * (rho - rhoRef);

	return pres;
}
*/


#endif