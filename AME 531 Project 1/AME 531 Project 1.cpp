#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> // Include algorithm for std::max_element
#include <fstream>
#include <sstream>
#include <unordered_map>
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
std::vector<double> cfl_values; // CFL values
std::vector<double> time_steps; // Time steps
std::vector<double> u_max_values; // Maximum velocity values
std::vector<double> flow_rate_values; // Flow rate values

// Function prototypes
void initialize(int Ny);
void applyBoundaryConditions(int Ny);
void solveNavierStokes(int solverType, double dt, int numSteps, int Ny);
void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, double cfl);
void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, int Ny, const std::vector<double>& cfl_values, const std::vector<double>& time_steps, const std::vector<double>& u_max_values, const std::vector<double>& flow_rate_values);
std::unordered_map<std::string, double> readVariablesFromFile(const std::string& filename);

int main() {
    int solverType;
    std::unordered_map<std::string, double> variables = readVariablesFromFile("input.k");

    if (variables.find("dpdx") != variables.end()) {
        dpdx = variables["dpdx"];
    }
    else {
        std::cout << "Enter dp/dx: ";
        std::cin >> dpdx;
    }

    if (variables.find("rho") != variables.end()) {
        rho = variables["rho"];
    }
    else {
        std::cout << "Enter density (rho): ";
        std::cin >> rho;
    }

    if (variables.find("nu") != variables.end()) {
        nu = variables["nu"];
    }
    else {
        std::cout << "Enter kinematic viscosity (nu): ";
        std::cin >> nu;
    }

    if (variables.find("dt") != variables.end()) {
        dt = variables["dt"];
    }
    else {
        std::cout << "Enter time step (dt): ";
        std::cin >> dt;
    }

    if (variables.find("endTime") != variables.end()) {
        endTime = variables["endTime"];
    }
    else {
        std::cout << "Enter end time: ";
        std::cin >> endTime;
    }

    if (variables.find("dy") != variables.end()) {
        dy = variables["dy"];
    }
    else {
        std::cout << "Enter grid spacing (dy): ";
        std::cin >> dy;
    }

    if (variables.find("y_bottom") != variables.end()) {
        y_bottom = variables["y_bottom"];
    }
    else {
        std::cout << "Enter location of the bottom wall (y_bottom): ";
        std::cin >> y_bottom;
    }

    if (variables.find("y_top") != variables.end()) {
        y_top = variables["y_top"];
    }
    else {
        std::cout << "Enter location of the top wall (y_top): ";
        std::cin >> y_top;
    }

    if (variables.find("printInterval") != variables.end()) {
        printInterval = variables["printInterval"];
    }
    else {
        std::cout << "Enter time interval for printing velocity values: ";
        std::cin >> printInterval;
    }

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
    plt::ion(); // Enable interactive mode
    solveNavierStokes(solverType, dt, numSteps, Ny);
    //plt::ioff(); // Disable interactive mode

    // Wait for user input before terminating
    std::cout << "Press any key to exit..." << std::endl;
    std::cin.get();
    std::cin.get();

    return 0;
}

std::unordered_map<std::string, double> readVariablesFromFile(const std::string& filename) {
    std::unordered_map<std::string, double> variables;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string key;
            double value;
            if (iss >> key >> value) {
                variables[key] = value;
            }
        }
        file.close();
    }

    return variables;
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
        double u_max = *std::max_element(u.begin(), u.end());
        if (solverType == 0) {
            cfl = u_max * dt / dy;
            cfl_values.push_back(cfl);
            time_steps.push_back(t * dt);
        }

        // Store u_max values
        u_max_values.push_back(u_max);

        // Calculate and store flow rate
        double flow_rate = 0.0;
        for (int j = 0; j < Ny; ++j) {
            flow_rate += u[j] * dy;
        }
        flow_rate_values.push_back(flow_rate);

        if (t % printStepInterval == 0) {
            printVelocities(u, u_exact, cfl); // Print velocities and exact solution at specified intervals
            plotVelocities(u, u_exact, Ny, cfl_values, time_steps, u_max_values, flow_rate_values);
            plt::draw(); // Update the plot
            plt::pause(0.01); // Pause for a short time to create animation effect
            plt::clf(); // Clear the current figure
        }
    }
    plotVelocities(u, u_exact, Ny, cfl_values, time_steps, u_max_values, flow_rate_values); // Plot velocities and exact solution after the final time step
    plt::show(); // Keep the plot window open
    // Wait for user input before terminating
    std::cout << "N O R M A L  T E R M I N A T I O N" << std::endl;
    std::cin.get();
    std::cin.get();
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

void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, int Ny, const std::vector<double>& cfl_values, const std::vector<double>& time_steps, const std::vector<double>& u_max_values, const std::vector<double>& flow_rate_values) {
    std::vector<double> y(Ny);
    for (int j = 0; j < Ny; ++j) {
        y[j] = y_bottom + j * dy;
    }

    plt::subplot(4, 1, 1);
    plt::named_plot("Numerical", u, y, "r-");
    plt::named_plot("Exact", u_exact, y, "b--");
    plt::xlabel("u");
    plt::ylabel("y");
    plt::legend();

    plt::subplot(4, 1, 2);
    plt::semilogy(time_steps, cfl_values, "g-");
    plt::xlabel("Time");
    plt::ylabel("CFL");

    plt::subplot(4, 1, 3);
    plt::plot(time_steps, u_max_values, "b-");
    plt::xlabel("Time");
    plt::ylabel("u_max");

    plt::subplot(4, 1, 4);
    plt::plot(time_steps, flow_rate_values, "m-");
    plt::xlabel("Time");
    plt::ylabel("Flow Rate");

    plt::show();
}
