#include "hllc.hpp"

double HllcSolver::j_k(double r, double u, double S) noexcept {
  return r * (u - S);
}

double HllcSolver::S_c(double j_l, double j_r, double u_l, double u_r,
                       double p_l, double p_r) noexcept {
  double nominator = j_r * u_r - j_l * u_l + p_r - p_l;
  double denominator = j_r - j_l;
  return nominator / denominator;
}

double HllcSolver::p_c(double j_l, double j_r, double u_l, double u_r,
                       double p_l, double p_r) noexcept {
  double nominator = j_r * p_l - j_l * p_r - j_l * j_r * (u_r - u_l);
  double denominator = j_r - j_l;
  return nominator / denominator;
}

double HllcSolver::r_c_k(double r, double u, double S, double S_c) noexcept {
  return r * (S - u) / (S - S_c);
}

double HllcSolver::e_c_k(double e, double p, double p_c, double r,
                         double r_c) noexcept {
  return e - 0.5 * (p + p_c) * (1. / r_c - 1. / r);
}

std::pair<double, double> HllcSolver::getWaveSpeed() const noexcept {
  double u_l = this->left.getVelocity();
  double c_l = this->left.getSoundSpeed();

  double u_r = this->right.getVelocity();
  double c_r = this->right.getSoundSpeed();

  return std::make_pair(gsl_min(u_l, u_r) - gsl_max(c_l, c_r),
                        gsl_max(u_l, u_r) + gsl_max(c_l, c_r));
}

ConservativeVariable HllcSolver::getFlux() const noexcept {
  double p_l = this->left.getPressure();
  double r_l = this->left.getDensity();
  double u_l = this->left.getVelocity();
  double e_l = this->right.getInternalEnergy();

  double p_r = this->right.getPressure();
  double r_r = this->right.getDensity();
  double u_r = this->right.getVelocity();
  double e_r = this->right.getInternalEnergy();

  double s = 0.;
  auto [S_l, S_r] = this->getWaveSpeed();

  if (s <= S_l) {
    return to_flux(this->left);
  }

  if (s >= S_r) {
    return to_flux(this->right);
  }

  double j_l = this->j_k(r_l, u_l, S_l);
  double j_r = this->j_k(r_r, u_r, S_r);
  double S_c = this->S_c(j_l, j_r, u_l, u_r, p_l, p_r);
  double p_c = this->p_c(j_l, j_r, u_l, u_r, p_l, p_r);

  double r_c = s >= S_c ? this->r_c_k(r_r, u_r, S_r, S_c)
                        : this->r_c_k(r_l, u_l, S_l, S_c);
  double e_c = s >= S_c ? this->e_c_k(e_r, p_r, p_c, r_r, r_c)
                        : this->e_c_k(e_l, p_l, p_c, r_l, r_c);
  return std::make_tuple(r_c * S_c, p_c + r_c * gsl_pow_2(S_c),
                         r_c * (e_c + 0.5 * gsl_pow_2(S_c)) * S_c + p_c * S_c);
}
