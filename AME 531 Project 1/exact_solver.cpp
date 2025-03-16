#include "exact_solver.h"

void exactSolver(std::vector<double>& u_exact, double dpdx, double rho, double nu, double Ly, int Ny) {
    double dy = Ly / (Ny - 1);
    for (int j = 0; j < Ny; ++j) {
        double y = j * dy;
        u_exact[j] = (dpdx / (2 * nu * rho)) * (y * y - Ly * y);
    }
}