#include "runge_kutta_solver.h"

void rungeKuttaSolver(std::vector<double>& u, double dt, double dy, double nu, double dpdx, double rho) {
    std::vector<double> k1(u.size(), 0.0);
    std::vector<double> k2(u.size(), 0.0);
    std::vector<double> k3(u.size(), 0.0);
    std::vector<double> k4(u.size(), 0.0);
    std::vector<double> u_temp(u.size(), 0.0);

    // Calculate k1
    for (int j = 1; j < u.size() - 1; ++j) {
        k1[j] = dt * (nu * (u[j + 1] - 2 * u[j] + u[j - 1]) / (dy * dy) - dpdx / rho);
    }

    // Calculate k2
    for (int j = 1; j < u.size() - 1; ++j) {
        u_temp[j] = u[j] + 0.5 * k1[j];
    }
    for (int j = 1; j < u.size() - 1; ++j) {
        k2[j] = dt * (nu * (u_temp[j + 1] - 2 * u_temp[j] + u_temp[j - 1]) / (dy * dy) - dpdx / rho);
    }

    // Calculate k3
    for (int j = 1; j < u.size() - 1; ++j) {
        u_temp[j] = u[j] + 0.5 * k2[j];
    }
    for (int j = 1; j < u.size() - 1; ++j) {
        k3[j] = dt * (nu * (u_temp[j + 1] - 2 * u_temp[j] + u_temp[j - 1]) / (dy * dy) - dpdx / rho);
    }

    // Calculate k4
    for (int j = 1; j < u.size() - 1; ++j) {
        u_temp[j] = u[j] + k3[j];
    }
    for (int j = 1; j < u.size() - 1; ++j) {
        k4[j] = dt * (nu * (u_temp[j + 1] - 2 * u_temp[j] + u_temp[j - 1]) / (dy * dy) - dpdx / rho);
    }

    // Update u
    for (int j = 1; j < u.size() - 1; ++j) {
        u[j] = u[j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6.0;
    }
}
