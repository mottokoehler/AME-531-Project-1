#include "implicit_solver.h"
#include <vector>
#include <Eigen/Dense> // Include Eigen library for matrix operations

using namespace Eigen;

void implicitSolver(std::vector<double>& u, double dt, double dy, double nu, double dpdx, double rho) {
    int Ny = u.size();
    MatrixXd A = MatrixXd::Zero(Ny, Ny);
    VectorXd b = VectorXd::Zero(Ny);
    VectorXd u_vec = VectorXd::Zero(Ny);

    // Fill the matrix A and vector b
    for (int j = 1; j < Ny - 1; ++j) {
        A(j, j - 1) = -nu * dt / (dy * dy);
        A(j, j) = 1 + 2 * nu * dt / (dy * dy);
        A(j, j + 1) = -nu * dt / (dy * dy);
        b(j) = u[j] - dt * dpdx / rho;
    }

    // Apply boundary conditions
    A(0, 0) = 1;
    A(Ny - 1, Ny - 1) = 1;
    b(0) = 0;
    b(Ny - 1) = 0;

    // Solve the linear system
    u_vec = A.colPivHouseholderQr().solve(b);

    // Update the velocity vector
    for (int j = 0; j < Ny; ++j) {
        u[j] = u_vec(j);
    }
}