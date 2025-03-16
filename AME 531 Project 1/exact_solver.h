#ifndef EXACT_SOLVER_H
#define EXACT_SOLVER_H

#include <vector>

void exactSolver(std::vector<double>& u_exact, double dpdx, double rho, double nu, double Ly, int Ny);

#endif // EXACT_SOLVER_H