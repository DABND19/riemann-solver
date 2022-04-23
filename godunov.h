#pragma once

#include "riemann_problem.h"

static const double CFL = 0.2;

double U1(const GasParameters* flow);
double U2(const GasParameters* flow);
double U3(const GasParameters* flow);

double M(const GasParameters* flow, double knot_velocity);
double J(const GasParameters* flow, double knot_velocity);
double E(const GasParameters* flow, double knot_velocity);

double F1(const GasParameters* flow);
double F2(const GasParameters* flow);
double F3(const GasParameters* flow);

GasParameters from_conservative(double u1, double u2, double u3);

double calculus_schema(double u, double f_left, double f_right,
                       double x_left, double x_right, double dt);
