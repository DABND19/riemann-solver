#pragma once

#include <tuple>

#include "base.hpp"

class HllcSolver : public RiemannSolver {
  static ConservativeVariable U_c_k(const GasFlow& flow, double S_k, double S_c_) noexcept;
  double getContactVelocity() const noexcept;

 public:
  HllcSolver() = default;

  std::pair<double, double> getWaveSpeed() const noexcept override;
  ConservativeVariable getFlux() const noexcept override;
};
