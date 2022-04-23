#include "godunov.h"

#include <gsl/gsl_math.h>


double U1(const GasParameters* flow) {
  return flow->density;
}

double U2(const GasParameters* flow) {
  return flow->density * flow->velocity;
}

double U3(const GasParameters* flow) {
  return flow->density * total_energy(flow);
}

double M(const GasParameters* flow, double knot_velocity) {
  return flow->density * (flow->velocity - knot_velocity);
}

double J(const GasParameters* flow, double knot_velocity) {
  return flow->pressure + flow->velocity * M(flow, knot_velocity);
}

double E(const GasParameters* flow, double knot_velocity) {
  return total_energy(flow) * M(flow, knot_velocity) + flow->pressure * flow->velocity;
}

double F1(const GasParameters* flow) {
  return M(flow, 0.);
}

double F2(const GasParameters* flow) {
  return J(flow, 0.);
}

double F3(const GasParameters* flow) {
  return E(flow, 0.);
}

GasParameters from_conservative(double u1, double u2, double u3) {
  double density = u1;
  double velocity = u2 / density;
  double pressure = (GAMMA - 1.) * (u3 - 0.5 * density * gsl_pow_2(velocity));
  return (GasParameters){pressure, density, velocity};
}

double calculus_schema(double u, double f_left, double f_right,
                       double x_left, double x_right, double dt) {
  return u - dt / (x_right - x_left) * (f_right - f_left);
}
