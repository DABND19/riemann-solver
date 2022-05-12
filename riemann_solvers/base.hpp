#pragma once

#include <tuple>

#include "../gas_flow.hpp"

using ConservativeVariable = std::tuple<double, double, double>;

ConservativeVariable to_conservative(const GasFlow& flow) noexcept;
ConservativeVariable to_flux(const GasFlow& flow) noexcept;
GasFlow from_conservative(const ConservativeVariable& u) noexcept;

class RiemannSolver {
 protected:
  GasFlow left, right;

 public:
  RiemannSolver() = default;

  std::pair<GasFlow, GasFlow> getState() const noexcept;
  virtual void setState(const GasFlow& left, const GasFlow& right);

  virtual std::pair<double, double> getWaveSpeed() const noexcept = 0;
  virtual ConservativeVariable getFlux() const noexcept = 0;
};
