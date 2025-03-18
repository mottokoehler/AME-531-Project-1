#pragma once
#ifndef RUNGE_KUTTA_SOLVER_H
#define RUNGE_KUTTA_SOLVER_H

#include <vector>

void rungeKuttaSolver(std::vector<double>& u, double dt, double dy, double nu, double dpdx, double rho);

#endif // RUNGE_KUTTA_SOLVER_H
