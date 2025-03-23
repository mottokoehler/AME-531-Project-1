#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>                // Include algorithm for std::max_element
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include "explicit_solver.h"        // Include explicit solver
#include "implicit_solver.h"        // Include implicit solver
#include "exact_solver.h"           // Include exact solver
#include "runge_kutta_solver.h"     // Include Runge-Kutta solver
#include "matplotlibcpp.h"          // Include matplotlib-cpp for plotting
#include <cstdlib>                  // for std::exit()
#include "printdata.h"
//#include "plotdata.h"

namespace plt = matplotlibcpp;

// User-supplied variables
double dpdx;                        // Change in pressure over x-direction
double rho;                         // Density
double nu;                          // Kinematic viscosity
double dt;                          // Time step
double endTime;                     // End time
double dy;                          // Grid spacing in y-direction
double y_bottom;                    // Location of the bottom wall
double y_top;                       // Location of the top wall
double printInterval;               // Time interval for printing velocity values
double ctime_steps;                  // Current Time Step

// Define grid and variables
std::vector<double> u;
std::vector<double> u_exact;        // Exact solution
std::vector<double> cfl_values;     // CFL values
std::vector<double> time_steps;     // Time steps
std::vector<double> u_max_values;   // Maximum velocity values
std::vector<double> flow_rate_values;       // Flow rate values
std::vector<double> shear_stress_values;    // Shear stress values
double error;           // Error Calculation

bool cflWarningDisplayed = false; // Flag to track if the CFL warning has been displayed

// Function prototypes
void initialize(int Ny);
void applyBoundaryConditions(int Ny);
void solveNavierStokes(int solverType, double dt, int numSteps, int Ny, const std::string& analysisType);
//void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, double cfl, double flow_rateconst, double error, double ctime_steps, double shear_stress);
void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, int Ny, const std::vector<double>& cfl_values, const std::vector<double>& time_steps, const std::vector<double>& u_max_values, const std::vector<double>& flow_rate_values, const std::vector<double>& shear_stress_values);
std::unordered_map<std::string, double> readVariablesFromFile(const std::string& filename);

int main() {
    // Print the current working directory
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;
    std::cout << "Place input.k file in current working directory" << std::endl;
    std::cout << "Sample format for input.k:" << std::endl;
    std::cout << "dpdx -1  //input.k has variables to be read in and their values" << std::endl;
    std::cout << "................................................................" << std::endl;

    int solverType;
    std::string analysisType;

    // Read in variables from a file. If the variables do not exist, prompt the user to enter them.
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

    // Set analysisType based on solverType
    if (solverType == 0) {
        analysisType = "explicit";
    }
    else if (solverType == 1) {
        analysisType = "implicit";
    }
    else if (solverType == 2) {
        analysisType = "RK";
    }

    // Validate printInterval
    if (printInterval < dt || static_cast<int>(printInterval / dt) * dt != printInterval) {
        std::cerr << "Invalid print interval. It must be greater than or equal to the time step and a multiple of the time step." << std::endl;
        return 1;
    }

    // Calculate the number of grid points and time steps
    int Ny = static_cast<int>((y_top - y_bottom) / dy) + 1;
    int numSteps = static_cast<int>(1 + endTime / dt);

    // Initialize the grid and variables
    u.resize(Ny, 0.0); // Initialize the velocity vector
    u_exact.resize(Ny, 0.0); // Initialize the exact solution vector

    initialize(Ny);
    plt::figure_size(1200, 300);
    solveNavierStokes(solverType, dt, numSteps, Ny, analysisType); // Solve the Navier-Stokes equations

    plt::show();

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

void solveNavierStokes(int solverType, double dt, int numSteps, int Ny, const std::string& analysisType) {
    int printStepInterval = static_cast<int>(printInterval / dt);
    for (int t = 1; t < numSteps; ++t) {

        if (solverType == 1) {
            implicitSolver(u, dt, dy, nu, dpdx, rho);               // Call Implicit
        }
        else if (solverType == 2) {
            rungeKuttaSolver(u, dt, dy, nu, dpdx, rho);             // Call 4th Order RK
        }
        else {
            explicitSolver(u, dt, dy, nu, dpdx, rho);               // Call Explicit
        }
        applyBoundaryConditions(Ny);
        exactSolver(u_exact, dpdx, rho, nu, y_top - y_bottom, Ny); // Calculate exact solution

        // Calculate CFL number if explicit solver is used
        double cfl = 0.0;
        double u_max = *std::max_element(u.begin(), u.end());
        double u_max_exact = *std::max_element(u_exact.begin(), u_exact.end());
        double error = 100 * (u_max_exact - u_max) / u_max_exact;
        cfl = u_max * dt / dy;

        cfl_values.push_back(cfl);
        time_steps.push_back(t * dt);
        ctime_steps = ctime_steps + dt;

        // Store u_max values
        u_max_values.push_back(u_max);

        // Calculate and store flow rate
        double flow_rate = 0.0;
        for (int j = 0; j < Ny; ++j) {
            flow_rate += u[j] * dy;
        }
        flow_rate_values.push_back(flow_rate);

        // Calculate and store shear stress at the wall
        double shear_stress = nu * (u[1] - u[0]) / dy;
        shear_stress_values.push_back(shear_stress);

        //Issue Warning if CFL exceeds 1, option to about or continue
        if (cfl > 1 && !cflWarningDisplayed) {
            std::cout << "CFL exceeds safe limit of 1, solution likely unstable" << std::endl;
            std::cout << "Do you want to exit the program? (y/n): ";
            char choice;
            std::cin >> choice;
            if (choice == 'y' || choice == 'Y') {
                printVelocities(u, u_exact, cfl, flow_rate, error, ctime_steps, shear_stress, analysisType, dt, dy);
                exit(0);
            }
            else {
                cflWarningDisplayed = true;
            }
		}
        if (t % printStepInterval == 0) {
            printVelocities(u, u_exact, cfl, flow_rate, error, ctime_steps, shear_stress, analysisType, dt, dy); // Print velocities, exact solution, CFL, and flow rate at specified intervals
            plotVelocities(u, u_exact, Ny, cfl_values, time_steps, u_max_values, flow_rate_values, shear_stress_values);
            plt::draw();            // Update the plot
            plt::pause(0.01);       // Pause for a short time to create animation effect
            plt::clf();             // Clear the current figure
        }
    }
    printVelocities(u, u_exact, cfl_values.back(), flow_rate_values.back(), error, ctime_steps, shear_stress_values.back(), analysisType, dt, dy); // Print final velocities and exact solution
    plotVelocities(u, u_exact, Ny, cfl_values, time_steps, u_max_values, flow_rate_values, shear_stress_values); // Plot velocities and exact solution after the final time step

}

void plotVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, int Ny, const std::vector<double>& cfl_values, const std::vector<double>& time_steps, const std::vector<double>& u_max_values, const std::vector<double>& flow_rate_values, const std::vector<double>& shear_stress_values) {
    std::vector<double> y(Ny);
    //plt::figure_size(800, 300);
    for (int j = 0; j < Ny; ++j) {
        y[j] = y_bottom + j * dy;
    }

    
    plt::subplot(1, 4, 1);
    plt::named_plot("Numerical", u, y, "r-");
    plt::named_plot("Exact", u_exact, y, "b--");
    plt::xlabel("Velocity, u");
    plt::ylabel("Length, y");
    plt::legend();

    plt::subplot(1, 4, 2);
    plt::semilogy(time_steps, cfl_values, "g-");
    plt::xlabel("Time");
    plt::ylabel("CFL");

    //plt::subplot(2, 2, 3);
    //plt::plot(time_steps, u_max_values, "b-");
    //plt::xlabel("Time");
    //plt::ylabel("u_max");

    plt::subplot(1, 4, 3);
    plt::plot(time_steps, flow_rate_values, "m-");
    plt::xlabel("Time");
    plt::ylabel("Flow Rate, Q");

    plt::subplot(1, 4, 4);
    plt::plot(time_steps, shear_stress_values, "c-");
    plt::xlabel("Time");
    plt::ylabel("Shear Stress");

	plt::tight_layout();

    //plt::show();
}

