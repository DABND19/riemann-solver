#pragma once

#include <memory>
#include <vector>

#include "riemann_solvers/base.hpp"

class BaseGodunovSolver {
  static const double CFL;

  double t_current = 0.;
  double dt = 0.;

  std::vector<double> x;
  std::vector<ConservativeVariable> u;
  std::vector<ConservativeVariable> F;

  std::shared_ptr<RiemannSolver> solver;

  void calculateFluxes();

 protected:
  ConservativeVariable getLeftFlux(size_t i_cell) const noexcept;
  ConservativeVariable getRightFlux(size_t i_cell) const noexcept;
  double getCellSize(size_t i_cell) const noexcept;

  virtual ConservativeVariable differenceSchema(const size_t i) const noexcept = 0;

  virtual ConservativeVariable getLeftBoundFlux() const noexcept = 0;
  virtual ConservativeVariable getRightBoundFlux() const noexcept = 0;

  virtual GasFlow getLeftCellFlow(size_t i_knot) const;
  virtual GasFlow getRightCellFlow(size_t i_knot) const;

 public:
  BaseGodunovSolver() = default;

  const std::vector<double>& getX() const noexcept;
  const std::vector<ConservativeVariable>& getU() const noexcept;
  const std::vector<ConservativeVariable>& getF() const noexcept;

  double getCurrentTime() const noexcept;
  double getTimeStep() const noexcept;

  void setSolver(std::shared_ptr<RiemannSolver> solver) noexcept;
  void setInitialConditions(const std::vector<double>& knots,
                            const std::vector<GasFlow>& flow) noexcept;

  void run(double t_end);
};
