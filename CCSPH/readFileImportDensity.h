#pragma once
#include "SPH2DCPPCuda.h"
#include <vector>
#include <string>
#ifndef READFILEIMPORTDENSITY_H
#define READFILEIMPORTDENSITY_H

int readFileImportDensity(int* nP, particleStructure* particles, paramsType* params, std::vector<struct kinematicsFunctionStructure>* kinematicsFunction, std::string inputFile);

#endif