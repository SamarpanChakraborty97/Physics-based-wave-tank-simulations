#pragma once
#include "SPH2DCPPCuda.h"
#include <vector>
#include <string>

#ifndef SPH2DCUDA_H
#define SPH2DCUDA_H

void SPH2DCuda(particleStructure * particles, paramsType * params, std::vector<kinematicsFunctionStructure> * kinematicsFunction, std::string outputDir);

#endif