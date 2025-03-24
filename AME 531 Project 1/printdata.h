#pragma once
#ifndef PRINTDATA_H
#define PRINTDATA_H

#include <vector>
#include <string> // Include the string header

void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, double cfl, double flow_rate, double error, double ctime_steps, double shear_stress, const std::string& analysisType, double dt, double dy, double avg_timestep);

#endif // PRINTDATA_H