#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> // Include algorithm for std::max_element
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
double dy; // Grid spacing in y-direction
double y_bottom; // Location of the bottom wall
double y_top; // Location of the top wall
double printInterval; // Time interval for printing velocity values

// Define grid and variables
std::vector<double> u;
std::vector<double> u_exact; // Exact solution

// Function prototypes
void initialize(int Ny);
void applyBoundaryConditions(int Ny);
void solveNavierStokes(int solverType, double dt, int numSteps, int Ny);
void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, double cfl);
void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, int Ny);

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
    std::cout << "Enter grid spacing (dy): ";
    std::cin >> dy;
    std::cout << "Enter location of the bottom wall (y_bottom): ";
    std::cin >> y_bottom;
    std::cout << "Enter location of the top wall (y_top): ";
    std::cin >> y_top;
    std::cout << "Enter time interval for printing velocity values: ";
    std::cin >> printInterval;
    std::cout << "Select solver type (0 for explicit, 1 for implicit, 2 for Runge-Kutta): ";
    std::cin >> solverType;

    // Validate printInterval
    if (printInterval < dt || static_cast<int>(printInterval / dt) * dt != printInterval) {
        std::cerr << "Invalid print interval. It must be greater than or equal to the time step and a multiple of the time step." << std::endl;
        return 1;
    }

    int Ny = static_cast<int>((y_top - y_bottom) / dy) + 1;
    int numSteps = static_cast<int>(endTime / dt);

    u.resize(Ny, 0.0);
    u_exact.resize(Ny, 0.0);

    initialize(Ny);
    solveNavierStokes(solverType, dt, numSteps, Ny);

    return 0;
}

void initialize(int Ny) {
    // Initialize the grid and variables
    for (int j = 0; j < Ny; ++j) {
        u[j] = 0.0;
        u_exact[j] = 0.0;
    }
}

void applyBoundaryConditions(int Ny) {
    // Apply no-slip boundary conditions
    u[0] = 0.0;
    u[Ny - 1] = 0.0;
}

void solveNavierStokes(int solverType, double dt, int numSteps, int Ny) {
    int printStepInterval = static_cast<int>(printInterval / dt);
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
        applyBoundaryConditions(Ny);
        exactSolver(u_exact, dpdx, rho, nu, y_top - y_bottom, Ny); // Calculate exact solution

        // Calculate CFL number if explicit solver is used
        double cfl = 0.0;
        if (solverType == 0) {
            double u_max = *std::max_element(u.begin(), u.end());
            cfl = u_max * dt / dy;
        }

        if (t % printStepInterval == 0) {
            printVelocities(u, u_exact, cfl); // Print velocities and exact solution at specified intervals
        }
        //plotVelocities(u, u_exact, Ny); // Plot velocities
    }
    plotVelocities(u, u_exact, Ny); // Plot velocities and exact solution after the final time step
}

void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, double cfl) {
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

    std::cout << "CFL Number: " << cfl << std::endl;
}

void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, int Ny) {
    std::vector<double> y(Ny);
    for (int j = 0; j < Ny; ++j) {
        y[j] = y_bottom + j * dy;
    }

    plt::named_plot("Numerical", u, y, "r-");
    plt::named_plot("Exact", u_exact, y, "b--");
    plt::xlabel("y");
    plt::ylabel("u");
    plt::legend();
    plt::show();
}
