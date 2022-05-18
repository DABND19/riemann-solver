#include "base.hpp"

#include <gsl/gsl_math.h>

#include <tuple>

ConservativeVariable to_conservative(const GasFlow& flow) noexcept {
  double r = flow.getDensity();
  double u = flow.getVelocity();
  double e = flow.getTotalEnergy();
  return std::make_tuple(r, r * u, r * e);
}

ConservativeVariable to_flux(const GasFlow& flow) noexcept {
  double p = flow.getPressure();
  double r = flow.getDensity();
  double u = flow.getVelocity();
  double e = flow.getTotalEnergy();
  return std::make_tuple(r * u, p + r * gsl_pow_2(u), r * e * u + p * u);
}

GasFlow from_conservative(const ConservativeVariable& u) noexcept {
  auto [u1, u2, u3] = u;
  double density = u1;
  double velocity = u2 / u1;
  double pressure =
      (GasFlow::GAMMA - 1.) * (u3 - 0.5 * density * gsl_pow_2(velocity));
  return GasFlow(pressure, density, velocity);
}

ConservativeVariable operator+(const ConservativeVariable& left,
                               const ConservativeVariable& right) {
  auto [left_1, left_2, left_3] = left;
  auto [right_1, right_2, right_3] = right;
  return std::make_tuple(left_1 + right_1, left_2 + right_2, left_3 + right_3);
}

ConservativeVariable operator-(const ConservativeVariable& left,
                               const ConservativeVariable& right) {
  return left + (-right);
}

ConservativeVariable operator*(double alpha,
                               const ConservativeVariable& value) {
  auto [value_1, value_2, value_3] = value;
  return std::make_tuple(alpha * value_1, alpha * value_2, alpha * value_3);
}

ConservativeVariable operator*(const ConservativeVariable& value,
                               double alpha) {
  return alpha * value;
}

ConservativeVariable operator-(const ConservativeVariable& value) {
  return -1. * value;
}

ConservativeVariable operator/(const ConservativeVariable& value, double alpha) {
  return value * (1. / alpha);
}

std::pair<GasFlow, GasFlow> RiemannSolver::getState() const noexcept {
  return std::make_pair(this->left, this->right);
}

void RiemannSolver::setState(const GasFlow& left, const GasFlow& right) {
  this->left = left;
  this->right = right;
}
