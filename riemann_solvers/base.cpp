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
  double pressure = (GasFlow::GAMMA - 1.) * (u3 - 0.5 * density * gsl_pow_2(velocity));
  return GasFlow(pressure, density, velocity);
}

std::pair<GasFlow, GasFlow> RiemannSolver::getState() const noexcept {
  return std::make_pair(this->left, this->right);
}

void RiemannSolver::setState(const GasFlow& left, const GasFlow& right) {
  this->left = left;
  this->right = right;
}
