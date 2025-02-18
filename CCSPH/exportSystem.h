#pragma once
#include <string>
#ifndef EXPORTSYSTEM_H
#define EXPORTSYSTEM_H

int exportSystem(struct particleStructure* particles, unsigned int dataFileNumber, struct paramsType* params, std::string outputDir);

#endif