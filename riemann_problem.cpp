#pragma once

#include "riemann_problem.hpp"

ConservativeVariable RiemannProblemGodunovSolver::differenceSchema(const size_t i) const noexcept {
  auto F_left = this->getLeftFlux(i);
  auto F_right = this->getRightFlux(i);
  auto dx = this->getCellSize(i);
  auto dt = this->getTimeStep();
  auto u = this->getU()[i];
  return u - dt / dx * (F_right - F_left);
}

ConservativeVariable RiemannProblemGodunovSolver::getLeftBoundFlux() const noexcept {
  auto flow = from_conservative(this->getU().front());
  return to_flux(flow);
}

ConservativeVariable RiemannProblemGodunovSolver::getRightBoundFlux() const noexcept {
  auto flow = from_conservative(this->getU().back());
  return to_flux(flow);
}
