#include "hll.hpp"

#include <tuple>

std::pair<double, double> HllSolver::getWaveSpeed() const noexcept {
  double u_l = this->left.getVelocity();
  double c_l = this->left.getSoundSpeed();

  double u_r = this->right.getVelocity();
  double c_r = this->right.getSoundSpeed();

  return std::make_pair(gsl_min(u_l, u_r) - gsl_max(c_l, c_r),
                        gsl_max(u_l, u_r) + gsl_max(c_l, c_r));
}

ConservativeVariable HllSolver::getFlux() const noexcept {
  double s = 0.;
  auto [S_l, S_r] = this->getWaveSpeed();

  auto u_l = to_conservative(this->left);
  auto F_l = to_flux(this->left);

  auto u_r = to_conservative(this->right);
  auto F_r = to_flux(this->right);

  if (s <= S_l) {
    return F_l;
  }

  if (s >= S_r) {
    return F_r;
  }

  return (S_r * F_l - S_l * F_r + S_l * S_r * (u_r - u_l)) / (S_r - S_l);
}
