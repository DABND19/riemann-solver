#pragma once

#include "gas_flow.h"
#include "riemann_problem.h"

static const double CFL = 0.2;

double U1(const GasFlow* flow);
double U2(const GasFlow* flow);
double U3(const GasFlow* flow);

double M(const GasFlow* flow, double knot_velocity);
double J(const GasFlow* flow, double knot_velocity);
double E(const GasFlow* flow, double knot_velocity);

double F1(const GasFlow* flow);
double F2(const GasFlow* flow);
double F3(const GasFlow* flow);

GasFlow from_conservative(double u1, double u2, double u3);

double flat_schema(double u, double f_left, double f_right,
                   double x_left, double x_right, double dt);

double b1(const GasFlow* flow);
double b2(const GasFlow* flow);
double b3(const GasFlow* flow);
double spherical_schema(double u, 
                       double f_left, double f_right,
                       double b, double B,
                       double x_left, double x_right,
                       double dt);
