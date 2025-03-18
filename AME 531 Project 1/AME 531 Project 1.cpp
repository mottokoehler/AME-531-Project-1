#include <iostream>
#include <vector>
#include <cmath>
#include "explicit_solver.h"
#include "implicit_solver.h"
#include "exact_solver.h"
#include "runge_kutta_solver.h" // Include Runge-Kutta solver
#include "matplotlibcpp.h" // Include matplotlib-cpp

namespace plt = matplotlibcpp;

// User-supplied variables
double dpdx; // Change in pressure over x-direction
double rho;  // Density
double nu;   // Kinematic viscosity
double dt;   // Time step
double endTime; // End time

// Define constants
const int Ny = 21;       // Number of grid points in y-direction
const double Ly = 2.0;   // Length of the domain in y-direction

// Define grid and variables
std::vector<double> u(Ny, 0.0);
std::vector<double> u_exact(Ny, 0.0); // Exact solution

// Function prototypes
void initialize();
void applyBoundaryConditions();
void solveNavierStokes(int solverType, double dt, int numSteps);
void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact);
void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact);

int main() {
    int solverType;
    std::cout << "Enter dp/dx: ";
    std::cin >> dpdx;
    std::cout << "Enter density (rho): ";
    std::cin >> rho;
    std::cout << "Enter kinematic viscosity (nu): ";
    std::cin >> nu;
    std::cout << "Enter time step (dt): ";
    std::cin >> dt;
    std::cout << "Enter end time: ";
    std::cin >> endTime;
    std::cout << "Select solver type (0 for explicit, 1 for implicit, 2 for Runge-Kutta): ";
    std::cin >> solverType;

    int numSteps = static_cast<int>(endTime / dt);

    initialize();
    solveNavierStokes(solverType, dt, numSteps);

    return 0;
}

void initialize() {
    // Initialize the grid and variables
    for (int j = 0; j < Ny; ++j) {
        u[j] = 0.0;
        u_exact[j] = 0.0;
    }
}

void applyBoundaryConditions() {
    // Apply no-slip boundary conditions
    u[0] = 0.0;
    u[Ny - 1] = 0.0;
}

void solveNavierStokes(int solverType, double dt, int numSteps) {
    double dy = Ly / (Ny - 1);
    for (int t = 0; t < numSteps; ++t) {
        if (solverType == 1) {
            implicitSolver(u, dt, dy, nu, dpdx, rho);
        }
        else if (solverType == 2) {
            rungeKuttaSolver(u, dt, dy, nu, dpdx, rho);
        }
        else {
            explicitSolver(u, dt, dy, nu, dpdx, rho);
        }
        applyBoundaryConditions();
        exactSolver(u_exact, dpdx, rho, nu, Ly, Ny); // Calculate exact solution
        printVelocities(u, u_exact); // Print velocities and exact solution after each time step
        //plotVelocities(u, u_exact); // Plot velocities
    }
    plotVelocities(u, u_exact); // Plot velocities and exact solution after the final time step
}

void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact) {
    std::cout << "Velocities: ";
    for (const auto& vel : u) {
        std::cout << vel << " ";
    }
    std::cout << std::endl;

    std::cout << "Exact Velocities: ";
    for (const auto& vel : u_exact) {
        std::cout << vel << " ";
    }
    std::cout << std::endl;
}

void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact) {
    std::vector<double> y(Ny);
    double dy = Ly / (Ny - 1);
    for (int j = 0; j < Ny; ++j) {
        y[j] = j * dy;
    }

    plt::named_plot("Numerical", u, y, "r-");
    plt::named_plot("Exact", u_exact, y, "b--");
    plt::xlabel("y");
    plt::ylabel("u");
    plt::legend();
    plt::show();
}
