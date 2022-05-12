#pragma once

#include <tuple>

#include "base.hpp"

class HllcSolver : public RiemannSolver {
  static double j_k(double r, double u, double S) noexcept;
  static double S_c(double j_l, double j_r, double u_l, double u_r, double p_l,
                    double p_r) noexcept;
  static double p_c(double j_l, double j_r, double u_l, double u_r, double p_l,
                    double p_r) noexcept;
  static double r_c_k(double r, double u, double S, double S_c) noexcept;
  static double e_c_k(double e, double p, double p_c, double r,
                      double r_c) noexcept;

 public:
  HllcSolver() = default;

  std::pair<double, double> getWaveSpeed() const noexcept override;
  ConservativeVariable getFlux() const noexcept override;
};
