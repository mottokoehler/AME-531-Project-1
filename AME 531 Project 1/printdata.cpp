#include "printdata.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> // For std::max_element

void printVelocities(const std::vector<double>& u, const std::vector<double>& u_exact, double cfl, double flow_rate, double error, double ctime_steps, double shear_stress, const std::string& analysisType, double dt, double dy, double avg_timestep) {
    // Construct the filenames based on analysis type, dt, and dy
    std::ostringstream velocityFilename;
    velocityFilename << "results_velocity_" << analysisType << "_dt" << dt << "_dy" << dy << ".out";

    std::ostringstream cflFilename;
    cflFilename << "results_cfl_" << analysisType << "_dt" << dt << "_dy" << dy << ".out";

    std::ostringstream flowRateFilename;
    flowRateFilename << "results_flow_rate_" << analysisType << "_dt" << dt << "_dy" << dy << ".out";

    std::ostringstream shearStressFilename;
    shearStressFilename << "results_shear_stress_" << analysisType << "_dt" << dt << "_dy" << dy << ".out";

    std::ostringstream summaryFilename;
    summaryFilename << "results_summary_" << analysisType << "_dt" << dt << "_dy" << dy << ".out";

    std::ofstream outf{ velocityFilename.str(), std::ios::trunc };  // Velocity out file
    outf << "Time = ";
    outf << ctime_steps;
    outf << "\n";
    for (const auto& vel : u) {
        outf << vel;
        outf << ",";
    }
    outf << "\n";

    std::ofstream outf2{ cflFilename.str(), std::ios::trunc };
    outf2 << "Time = ";
    outf2 << ctime_steps;
    outf2 << "\n";
    outf2 << cfl;
    outf2 << "\n";

    std::ofstream outf3{ flowRateFilename.str(), std::ios::trunc };
    outf3 << "Time = ";
    outf3 << ctime_steps;
    outf3 << "\n";
    outf3 << flow_rate;
    outf3 << "\n";

    std::ofstream outf4{ shearStressFilename.str(), std::ios::trunc };
    outf4 << "Time = ";
    outf4 << ctime_steps;
    outf4 << "\n";
    outf4 << shear_stress;
    outf4 << "\n";

    // Find the maximum velocity
    double max_velocity = *std::max_element(u.begin(), u.end());

    // Write summary data to the summary file
    std::ofstream outfSummary{ summaryFilename.str(), std::ios::trunc };
    outfSummary << "Time = " << ctime_steps << "\n";
    outfSummary << "Max Velocity = " << max_velocity << "\n";
    outfSummary << "CFL = " << cfl << "\n";
    outfSummary << "Flow Rate = " << flow_rate << "\n";
    outfSummary << "Shear Stress = " << shear_stress << "\n";
    outfSummary << "dy = " << dy << "\n";
    outfSummary << "dt = " << dt << "\n";
    outfSummary << "Average Timestep Length = " << avg_timestep << " seconds\n";

    // Print values to screen //
    std::cout << "Time Step: " << ctime_steps << std::endl;
    std::cout << "Velocities: ";  // Print Velocities from Numerical
    for (const auto& vel : u) {
        std::cout << vel << " ";
    }
    std::cout << std::endl;

    std::cout << "Exact Velocities: ";  // Print Velocities from Exact
    for (const auto& vel : u_exact) {
        std::cout << vel << " ";
    }
    std::cout << std::endl;

    std::cout << "CFL Number: " << cfl << std::endl;  // Print CFL (should add something to stop this for implicit)
    std::cout << "Flow Rate: " << flow_rate << std::endl;  // Print Flow Rate
    std::cout << "% Error for u_max = " << error << std::endl;  // Print % Error
    std::cout << "Average Timestep Length: " << avg_timestep << " seconds" << std::endl;
}
