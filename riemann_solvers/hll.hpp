#pragma once

#include <tuple>

#include "base.hpp"

class HllSolver : public RiemannSolver {
 public:
  HllSolver() = default;

  std::pair<double, double> getWaveSpeed() const noexcept override;
  ConservativeVariable getFlux() const noexcept override;
};
