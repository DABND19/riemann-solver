#pragma once

#include <gsl/gsl_roots.h>

#include <memory>
#include <tuple>

#include "base.hpp"

class ExactSolver : public RiemannSolver {
  gsl_root_fsolver* bisection_solver;
  gsl_root_fdfsolver* newton_solver;

  static double Q_k(double p_c, const GasFlow& flow) noexcept;
  static double f_k(double p_c, const GasFlow& flow) noexcept;
  static double U_c(double p_c, const GasFlow& left, const GasFlow& right) noexcept;
  static double df_k(double p_c, const GasFlow& flow) noexcept;
  static void fdf_k(double p_c, const GasFlow& flow, double* f_, double* df_) noexcept;
  static double f(double p_c, const GasFlow& left, const GasFlow& right) noexcept;
  static double df(double p_c, const GasFlow& left, const GasFlow& right) noexcept;
  static void fdf(double p_c, const GasFlow& left, const GasFlow& right, double* f_, double* df_);
  static double gsl_f(double x, void* parameters) noexcept;
  static double gsl_df(double x, void* parameters) noexcept;
  static void gsl_fdf(double x, void* parameters, double* y, double* dy) noexcept;
  static double P_TR(const GasFlow& left, const GasFlow& right) noexcept;

  static double S_shock_left(double p_c, const GasFlow& left) noexcept;
  GasFlow left_shock_wave(double s) const noexcept;
  static double S_rarefraction_left(double p_c, const GasFlow& left) noexcept;
  GasFlow left_rarefraction_wave(double s) const noexcept;

  static double S_shock_right(double p_c, const GasFlow& right) noexcept;
  GasFlow right_shock_wave(double s) const noexcept;
  static double S_rarefraction_right(double p_c, const GasFlow& right) noexcept;
  GasFlow right_rarefraction_wave(double s) const noexcept;

  double p_contact;
  void calculateContactPressure();

 public:
  ExactSolver();
  ~ExactSolver();

  void setState(const GasFlow& left, const GasFlow& right) override;
  std::pair<double, double> getWaveSpeed() const noexcept override;
  ConservativeVariable getFlux() const noexcept override;
};
