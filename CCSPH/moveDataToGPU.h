#pragma once
#include "SPH2DCPPCuda.h"
#include <vector>

#ifndef MOVEDATATOGPU_H
#define MOVEDATATOGPU_H

int moveDataToGPU(particleStructure* particles, paramsType* params, std::vector<kinematicsFunctionStructure>* kinematicsFunction, particleStructure** ppdParticles, struct paramsType** ppdParams, struct kinematicsFunctionStructure** ppdKinematicsFunction, particleStructure** pdParticlesHostMirror);

#endif