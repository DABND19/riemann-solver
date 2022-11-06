#include "base_godunov.hpp"

#include <gsl/gsl_math.h>

#include <cmath>
#include <stdexcept>

const double BaseGodunovSolver::CFL = 0.2;

void BaseGodunovSolver::setSolver(
    std::shared_ptr<RiemannSolver> solver) noexcept {
  this->solver = solver;
}

void BaseGodunovSolver::setInitialConditions(
    const std::vector<double>& knots,
    const std::vector<GasFlow>& flow) noexcept {
  size_t dimension = knots.size();

  this->x = knots;
  this->F = std::vector<ConservativeVariable>(dimension);

  this->u = std::vector<ConservativeVariable>(dimension - 1);
  for (size_t i = 0; i < u.size(); ++i) {
    u[i] = to_conservative(flow[i]);
  }
}

const std::vector<double>& BaseGodunovSolver::getX() const noexcept {
  return this->x;
}

const std::vector<ConservativeVariable>& BaseGodunovSolver::getU()
    const noexcept {
  return this->u;
}

const std::vector<ConservativeVariable>& BaseGodunovSolver::getF()
    const noexcept {
  return this->F;
}

double BaseGodunovSolver::getCurrentTime() const noexcept {
  return this->t_current;
}

double BaseGodunovSolver::getTimeStep() const noexcept { return this->dt; }

ConservativeVariable BaseGodunovSolver::getLeftFlux(
    size_t i_cell) const noexcept {
  return this->F[i_cell];
}

ConservativeVariable BaseGodunovSolver::getRightFlux(
    size_t i_cell) const noexcept {
  return this->F[i_cell + 1];
}

double BaseGodunovSolver::getCellSize(size_t i_cell) const noexcept {
  return this->x[i_cell + 1] - this->x[i_cell];
}

GasFlow BaseGodunovSolver::getLeftCellFlow(size_t i_knot) const {
  if (i_knot == 0) {
    throw std::out_of_range("Left cell flow for the first flux is undefined.");
  }
  return from_conservative(this->getU()[i_knot - 1]);
}

GasFlow BaseGodunovSolver::getRightCellFlow(size_t i_knot) const {
  if (i_knot == this->getF().size() - 1) {
    throw std::out_of_range("Right cell flow for the last flux is undefined.");
  }
  return from_conservative(this->getU()[i_knot]);
}

void BaseGodunovSolver::calculateFluxes() {
  this->F.front() = this->getLeftBoundFlux();
  this->F.back() = this->getRightBoundFlux();

  this->dt = INFINITY;
  for (size_t i = 1; i < this->F.size() - 1; ++i) {
    auto left_cell = this->getLeftCellFlow(i);
    auto right_cell = this->getRightCellFlow(i);

    this->solver->setState(left_cell, right_cell);
    this->F[i] = this->solver->getFlux();

    auto [left_wave, right_wave] = this->solver->getWaveSpeed();
    if (gsl_fcmp(left_wave, 0., 1e-16)) {
      double dx = x[i] - x[i - 1];
      this->dt = gsl_min(this->dt, this->CFL * dx / std::fabs(left_wave));
    }
    if (gsl_fcmp(right_wave, 0., 1e-16)) {
      double dx = x[i + 1] - x[i];
      this->dt = gsl_min(this->dt, this->CFL * dx / std::fabs(right_wave));
    }
  }
}

void BaseGodunovSolver::run(double t_end) {
  for (; this->t_current < t_end; this->t_current += this->dt) {
    this->calculateFluxes();
    this->dt = gsl_min(this->dt, t_end - this->t_current);

    for (size_t i = 0; i < this->u.size(); ++i) {
      u[i] = this->differenceSchema(i);
    }
  }
}
