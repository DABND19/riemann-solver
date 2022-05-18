#include "exact.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <stdlib.h>

#include <cmath>
#include <sstream>
#include <stdexcept>

const double EPSILON = 1e-6;

ExactSolver::ExactSolver() {
  this->bisection_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  this->newton_solver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
}

ExactSolver::~ExactSolver() {
  gsl_root_fsolver_free(this->bisection_solver);
  gsl_root_fdfsolver_free(this->newton_solver);
}

void ExactSolver::setState(const GasFlow& left, const GasFlow& right) {
  this->left = left;
  this->right = right;
  this->calculateContactPressure();
}

double ExactSolver::Q_k(double p_c, const GasFlow& flow) noexcept {
  double p_k = flow.getPressure();
  double r_k = flow.getDensity();
  double A_k = 2. / ((GasFlow::GAMMA + 1.) * r_k);
  double B_k = (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.) * p_k;
  return std::sqrt((p_c + B_k) / A_k);
}

double ExactSolver::f_k(double p_c, const GasFlow& flow) noexcept {
  double p = p_c;
  double p_k = flow.getPressure();
  if (p > p_k) {
    return (p - p_k) / Q_k(p, flow);
  } else {
    double a_k = flow.getSoundSpeed();
    double tmp =
        std::pow(p / p_k, 0.5 * (GasFlow::GAMMA - 1.) / GasFlow::GAMMA);
    return 2. * a_k / (GasFlow::GAMMA - 1.) * (tmp - 1.);
  }
}

double ExactSolver::U_c(double p_c, const GasFlow& left,
                        const GasFlow& right) noexcept {
  double u_l = left.getVelocity();
  double u_r = right.getVelocity();
  double f_l = f_k(p_c, left);
  double f_r = f_k(p_c, right);
  return 0.5 * (u_l + u_r) + 0.5 * (f_r - f_l);
}

double ExactSolver::df_k(double p_c, const GasFlow& flow) noexcept {
  double p = p_c;
  double p_k = flow.getPressure();
  double r_k = flow.getDensity();
  if (p > p_k) {
    double B_k = (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.) * p_k;
    return (1. - 0.5 * (p - p_k) / (B_k + p)) / Q_k(p, flow);
  } else {
    double a_k = flow.getSoundSpeed();
    return 1. / (r_k * a_k) *
           std::pow(p / p_k, -0.5 * (GasFlow::GAMMA + 1.) / GasFlow::GAMMA);
  }
}

void ExactSolver::fdf_k(double p_c, const GasFlow& flow, double* f_,
                        double* df_) noexcept {
  *f_ = f_k(p_c, flow);
  *df_ = df_k(p_c, flow);
}

double ExactSolver::f(double p_c, const GasFlow& left,
                      const GasFlow& right) noexcept {
  double du = right.getVelocity() - left.getVelocity();
  double f_l = f_k(p_c, left);
  double f_r = f_k(p_c, right);
  return f_l + f_r + du;
}

double ExactSolver::df(double p_c, const GasFlow& left,
                       const GasFlow& right) noexcept {
  double df_l = df_k(p_c, left);
  double df_r = df_k(p_c, right);
  return df_l + df_r;
}

void ExactSolver::fdf(double p_c, const GasFlow& left, const GasFlow& right,
                      double* f_, double* df_) {
  double du = right.getVelocity() - left.getVelocity();
  double f_l, f_r, df_l, df_r;
  fdf_k(p_c, left, &f_l, &df_l);
  fdf_k(p_c, right, &f_r, &df_r);
  *f_ = f_l + f_r + du;
  *df_ = df_l + df_r;
}

struct gsl_f_params {
  const GasFlow& left;
  const GasFlow& right;
};

double ExactSolver::gsl_f(double x, void* parameters) noexcept {
  struct gsl_f_params* params = (struct gsl_f_params*)parameters;
  return f(x, params->left, params->right);
}

double ExactSolver::gsl_df(double x, void* parameters) noexcept {
  struct gsl_f_params* params = (struct gsl_f_params*)parameters;
  return df(x, params->left, params->right);
}

void ExactSolver::gsl_fdf(double x, void* parameters, double* y,
                          double* dy) noexcept {
  struct gsl_f_params* params = (struct gsl_f_params*)parameters;
  fdf(x, params->left, params->right, y, dy);
}

double ExactSolver::P_TR(const GasFlow& left, const GasFlow& right) noexcept {
  double a_l = left.getSoundSpeed();
  double a_r = right.getSoundSpeed();
  double u_l = left.getVelocity();
  double u_r = right.getVelocity();
  double p_l = left.getPressure();
  double p_r = right.getPressure();
  double G = 0.5 * (GasFlow::GAMMA - 1.) / GasFlow::GAMMA;
  double nominator = (a_l + a_r - 0.5 * (GasFlow::GAMMA - 1.) * (u_r - u_l));
  double denominator = (a_l / std::pow(p_l, G) + a_r / std::pow(p_r, G));
  return std::pow(nominator / denominator, 1. / G);
}

const size_t MAX_ITER = 1000;

void ExactSolver::calculateContactPressure() {
  double p_min = gsl_min(this->left.getPressure(), this->right.getPressure());
  double p_max = gsl_max(this->left.getPressure(), this->right.getPressure());
  double f_min = f(p_min, this->left, this->right);
  double f_max = f(p_max, this->left, this->right);
  double f_vacuum = f(0., this->left, this->right);

  int error = GSL_SUCCESS;

  if (gsl_fcmp(f_vacuum, 0., EPSILON) >= 0) {
    this->p_contact = 0.;
    return;
  }

  if (!gsl_fcmp(f_max, 0., EPSILON)) {
    this->p_contact = p_max;
    return;
  }

  if (!gsl_fcmp(f_min, 0., EPSILON)) {
    this->p_contact = p_min;
    return;
  }

  if (gsl_fcmp(f_min, 0., EPSILON) > 0 && gsl_fcmp(f_max, 0., EPSILON) > 0) {
    this->p_contact = P_TR(this->left, this->right);
  }

  double p_0 = f_max >= 0 ? p_min : p_max;
  double p = p_0;
  struct gsl_f_params equation_params = {left, right};
  gsl_function_fdf fdf_equation;
  fdf_equation.f = &gsl_f;
  fdf_equation.df = &gsl_df;
  fdf_equation.fdf = &gsl_fdf;
  fdf_equation.params = &equation_params;
  gsl_root_fdfsolver_set(this->newton_solver, &fdf_equation, p_0);

  error = GSL_CONTINUE;
  for (size_t i = 0; i < MAX_ITER && error == GSL_CONTINUE; ++i) {
    error = gsl_root_fdfsolver_iterate(this->newton_solver);
    if (error != GSL_SUCCESS) {
      break;
    }

    p_0 = p;
    p = gsl_root_fdfsolver_root(this->newton_solver);

    error = gsl_root_test_delta(p, p_0, 0., EPSILON);
  }

  if (error && f_max < 0) {
    gsl_function f_equation;
    f_equation.function = &gsl_f;
    f_equation.params = &equation_params;

    double p_l = p_min;
    double p_r = p_max;

    gsl_root_fsolver_set(this->bisection_solver, &f_equation, p_l, p_r);

    error = GSL_CONTINUE;
    for (size_t i = 0; i < MAX_ITER && error == GSL_CONTINUE; ++i) {
      error = gsl_root_fsolver_iterate(this->bisection_solver);
      if (error != GSL_SUCCESS) {
        break;
      }

      p = gsl_root_fsolver_root(this->bisection_solver);
      p_l = gsl_root_fsolver_x_lower(this->bisection_solver);
      p_r = gsl_root_fsolver_x_upper(this->bisection_solver);

      error = gsl_root_test_interval(p_l, p_r, 0., EPSILON);
    }
  }

  if (error) {
    std::stringstream error_message;
    error_message << "Failed to calculate pressure: \n";
    error_message << "Left: " << this->left.getPressure() << '\t'
                  << this->left.getDensity() << '\t' << this->left.getVelocity()
                  << '\n';
    error_message << "Right: " << this->right.getPressure() << '\t'
                  << this->right.getDensity() << '\t'
                  << this->right.getVelocity() << '\n';
    throw std::runtime_error(error_message.str());
  }
  this->p_contact = p;
}

double ExactSolver::S_shock_left(double p_c, const GasFlow& left) noexcept {
  double r_l = left.getDensity();
  double u_l = left.getVelocity();
  double Q_l = Q_k(p_c, left);
  return u_l - Q_l / r_l;
}

GasFlow ExactSolver::left_shock_wave(double s) const noexcept {
  double r_l = this->left.getDensity();
  double p_l = this->left.getPressure();
  double p_c = this->p_contact;

  double S_l = S_shock_left(p_c, this->left);

  if (s <= S_l) {
    return this->left;
  } else {
    double p_ = p_c / p_l;
    double G = (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.);
    double r_c = r_l * (p_ + G) / (G * p_ + 1.);

    double u_c = U_c(p_c, this->left, this->right);

    return {p_c, r_c, u_c};
  }
}

double ExactSolver::S_rarefraction_left(double p_c,
                                        const GasFlow& left) noexcept {
  double a_l = left.getSoundSpeed();
  double u_l = left.getVelocity();
  return u_l - a_l;
}

GasFlow ExactSolver::left_rarefraction_wave(double s) const noexcept {
  double p_l = this->left.getPressure();
  double r_l = this->left.getDensity();
  double u_l = this->left.getVelocity();
  double a_l = left.getSoundSpeed();
  double p_c = this->p_contact;
  double u_c = U_c(p_c, this->left, this->right);

  double p_ = p_c / p_l;

  double a_c = a_l * std::pow(p_, 0.5 * (GasFlow::GAMMA - 1.) / GasFlow::GAMMA);
  double S_l = S_rarefraction_left(p_c, this->left);
  double S_c = u_c - a_c;

  if (s <= S_l) {
    return this->left;
  }

  double r_c = r_l * std::pow(p_, 1. / GasFlow::GAMMA);
  if (s >= S_c) {
    return {p_c, r_c, u_c};
  }

  double G1 = 2. / (GasFlow::GAMMA + 1.);
  double G2 = (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.);
  double r =
      r_l * std::pow(G1 + G2 / a_l * (u_l - s), 2. / (GasFlow::GAMMA - 1.));
  double u = G1 * (a_l + 0.5 * (GasFlow::GAMMA - 1.) * u_l + s);
  double p = p_l * std::pow(G1 + G2 / a_l * (u_l - s),
                            2. * GasFlow::GAMMA / (GasFlow::GAMMA - 1.));
  return {p, r, u};
}

double ExactSolver::S_shock_right(double p_c, const GasFlow& right) noexcept {
  double r_r = right.getDensity();
  double u_r = right.getVelocity();
  double Q_r = Q_k(p_c, right);
  return u_r + Q_r / r_r;
}

GasFlow ExactSolver::right_shock_wave(double s) const noexcept {
  double r_r = this->right.getDensity();
  double p_r = this->right.getPressure();
  double p_c = this->p_contact;

  double S_r = S_shock_right(p_c, this->right);

  if (s >= S_r) {
    return this->right;
  } else {
    double p_ = p_c / p_r;
    double G = (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.);
    double r_c = r_r * (p_ + G) / (G * p_ + 1.);

    double u_c = U_c(p_c, this->left, this->right);

    return {p_c, r_c, u_c};
  }
}

double ExactSolver::S_rarefraction_right(double p_c,
                                         const GasFlow& right) noexcept {
  double a_r = right.getSoundSpeed();
  double u_r = right.getVelocity();
  return u_r + a_r;
}

GasFlow ExactSolver::right_rarefraction_wave(double s) const noexcept {
  double p_r = this->right.getPressure();
  double r_r = this->right.getDensity();
  double u_r = this->right.getVelocity();
  double a_r = this->right.getSoundSpeed();
  double p_c = this->p_contact;
  double u_c = U_c(p_c, this->left, this->right);

  double p_ = p_c / p_r;

  double a_c = a_r * std::pow(p_, 0.5 * (GasFlow::GAMMA - 1.) / GasFlow::GAMMA);
  double S_r = S_rarefraction_right(p_c, this->right);
  double S_c = u_c + a_c;

  if (s >= S_r) {
    return this->right;
  }

  double r_c = r_r * std::pow(p_, 1. / GasFlow::GAMMA);
  if (s <= S_c) {
    return {p_c, r_c, u_c};
  }

  double G1 = 2. / (GasFlow::GAMMA + 1.);
  double G2 = (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.);
  double r =
      r_r * std::pow(G1 - G2 / a_r * (u_r - s), 2. / (GasFlow::GAMMA - 1.));
  double u = G1 * (-a_r + 0.5 * (GasFlow::GAMMA - 1.) * u_r + s);
  double p = p_r * std::pow(G1 - G2 / a_r * (u_r - s),
                            2. * GasFlow::GAMMA / (GasFlow::GAMMA - 1.));
  return {p, r, u};
}

std::pair<double, double> ExactSolver::getWaveSpeed() const noexcept {
  double left_wave, right_wave;
  if (this->p_contact > this->left.getPressure()) {
    left_wave = this->S_shock_left(this->p_contact, this->left);
  } else {
    left_wave = this->S_rarefraction_left(this->p_contact, this->left);
  }
  if (this->p_contact > this->right.getPressure()) {
    right_wave = this->S_shock_right(this->p_contact, this->right);
  } else {
    right_wave = this->S_rarefraction_right(this->p_contact, this->right);
  }
  return std::make_pair(left_wave, right_wave);
}

ConservativeVariable ExactSolver::getFlux() const noexcept {
  double s = 0.;
  auto S_c = this->U_c(this->p_contact, this->left, this->right);

  if (s < S_c) {
    if (this->p_contact > this->left.getPressure()) {
      return to_flux(this->left_shock_wave(s));
    } else {
      return to_flux(this->left_rarefraction_wave(s));
    }
  } else {
    if (this->p_contact > this->right.getPressure()) {
      return to_flux(this->right_shock_wave(s));
    } else {
      return to_flux(this->right_rarefraction_wave(s));
    }
  }
}
