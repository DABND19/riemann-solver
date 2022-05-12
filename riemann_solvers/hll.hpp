#pragma once

#include <tuple>

#include "base.hpp"

class HllSolver : public RiemannSolver {
  static double F_c(double u_l, double u_r, 
                    double F_l, double F_r, 
                    double S_l, double S_r) noexcept;

 public:
  HllSolver() = default;

  std::pair<double, double> getWaveSpeed() const noexcept override;
  ConservativeVariable getFlux() const noexcept override;
};
