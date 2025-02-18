#ifndef SPH3DCUDA_H
#define SPH3DCUDA_H

void SPH3DCuda(particleStructure* particles, paramsType* params, std::vector<kinematicsFunctionStructure>* kinematicsFunction, std::string outputDir);

#endif
