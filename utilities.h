#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <vector>
#include "bamhandle.h"

void calculateFootprintCoverage(std::vector<int> &fc, BamHandle *handle, const std::string &qName, int qStart, int qEnd);

double calculateAverageFootprintCoverage(const std::vector<int> &fc, int qStart, int qEnd);

#endif