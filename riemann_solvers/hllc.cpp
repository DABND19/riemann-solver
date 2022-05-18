#include "hllc.hpp"

ConservativeVariable HllcSolver::U_c_k(const GasFlow& flow, double S_k,
                                       double S_c_) noexcept {
  double r_k = flow.getDensity();
  double p_k = flow.getPressure();
  double u_k = flow.getVelocity();
  double E_k = flow.getTotalEnergy();

  return r_k * (S_k - u_k) / (S_k - S_c_) *
         std::make_tuple(
             1., S_c_, E_k + (S_c_ - u_k) * (S_c_ + p_k / (r_k * (S_k - u_k))));
}

std::pair<double, double> HllcSolver::getWaveSpeed() const noexcept {
  double u_l = this->left.getVelocity();
  double c_l = this->left.getSoundSpeed();

  double u_r = this->right.getVelocity();
  double c_r = this->right.getSoundSpeed();

  return std::make_pair(gsl_min(u_l, u_r) - gsl_max(c_l, c_r),
                        gsl_max(u_l, u_r) + gsl_max(c_l, c_r));
}

double HllcSolver::getContactVelocity() const noexcept {
  double p_l = this->left.getPressure();
  double r_l = this->left.getDensity();
  double u_l = this->left.getVelocity();

  double p_r = this->right.getPressure();
  double r_r = this->right.getDensity();
  double u_r = this->right.getVelocity();

  auto [S_l, S_r] = this->getWaveSpeed();

  double nominator =
      p_r - p_l + r_l * u_l * (S_l - u_l) - r_r * u_r * (S_r - u_r);
  double denominator = r_l * (S_l - u_l) - r_r * (S_r - u_r);
  return nominator / denominator;
}

ConservativeVariable HllcSolver::getFlux() const noexcept {
  auto F_l = to_flux(this->left);
  auto F_r = to_flux(this->right);

  auto [S_l, S_r] = this->getWaveSpeed();
  auto S_c = this->getContactVelocity();

  if (0 <= S_l) {
    return F_l;
  }

  if (0 >= S_r) {
    return F_r;
  }

  auto U_l = to_conservative(this->left);
  auto U_r = to_conservative(this->right);
  auto U_c_l = this->U_c_k(this->left, S_l, S_c);
  auto U_c_r = this->U_c_k(this->right, S_r, S_c);

  if (0 < S_c) {
    return F_l + S_l * (U_c_l - U_l);
  } else {
    return F_r + S_r * (U_c_r - U_r);
  }
}
