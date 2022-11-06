#pragma once

#include "base_godunov.hpp"

class RiemannProblemGodunovSolver : public BaseGodunovSolver {
 protected:
  ConservativeVariable differenceSchema(size_t i) const noexcept override;

  ConservativeVariable getLeftBoundFlux() const noexcept override;
  ConservativeVariable getRightBoundFlux() const noexcept override;
};
