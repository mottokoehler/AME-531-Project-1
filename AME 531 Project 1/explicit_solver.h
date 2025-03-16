#pragma once
#ifndef EXPLICIT_SOLVER_H
#define EXPLICIT_SOLVER_H

#include <vector>

void explicitSolver(std::vector<double>& u, double dt, double dy, double nu, double dpdx, double rho);

#endif // EXPLICIT_SOLVER_H
