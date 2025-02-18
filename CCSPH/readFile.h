#pragma once
#include "SPH2DCPPCuda.h"
#include <vector>
#include <string>

#ifndef READFILE_H
#define READFILE_H

int readFile(int* nP, particleStructure* particles, paramsType* params, std::vector<struct kinematicsFunctionStructure>* kinematicsFunction, std::string inputFile);

#endif