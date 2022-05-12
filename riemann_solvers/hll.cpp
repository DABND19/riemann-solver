#include "hll.hpp"

#include <tuple>

double HllSolver::F_c(double u_l, double u_r, double F_l, double F_r,
                      double S_l, double S_r) noexcept {
  return (S_r * F_l - S_l * F_r + S_l * S_r * (u_r - u_l)) / (S_r - S_l);
}

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

  auto [u1_l, u2_l, u3_l] = u_l;
  auto [F1_l, F2_l, F3_l] = F_l;

  auto [u1_r, u2_r, u3_r] = u_r;
  auto [F1_r, F2_r, F3_r] = F_r;

  double F1_c = HllSolver::F_c(u1_l, u1_r, F1_l, F1_r, S_l, S_r);
  double F2_c = HllSolver::F_c(u2_l, u2_r, F2_l, F2_r, S_l, S_r);
  double F3_c = HllSolver::F_c(u3_l, u3_r, F3_l, F3_r, S_l, S_r);
  return std::make_tuple(F1_c, F2_c, F3_c);
}
