#include "explicit_solver.h"

void explicitSolver(std::vector<double>& u, double dt, double dy, double nu, double dpdx, double rho) {
    std::vector<double> u_new(u.size(), 0.0);

    for (int j = 1; j < u.size() - 1; ++j) {
        u_new[j] = u[j] + dt * (nu * (u[j + 1] - 2 * u[j] + u[j - 1]) / (dy * dy) - dpdx / rho);
    }

    u = u_new;
}
