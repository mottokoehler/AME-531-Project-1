#pragma once
#ifndef IMPLICIT_SOLVER_H
#define IMPLICIT_SOLVER_H

#include <vector>

void implicitSolver(std::vector<double>& u, double dt, double dy, double nu, double dpdx, double rho);

#endif // IMPLICIT_SOLVER_H
