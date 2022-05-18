#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "riemann_solvers/base.hpp"

using BoundaryCondition =
    std::function<ConservativeVariable(const ConservativeVariable&, double)>;
using DifferenceSchema = std::function<ConservativeVariable(
    const ConservativeVariable&, const ConservativeVariable&,
    const ConservativeVariable&, double, double, double)>;

ConservativeVariable cartesian_schema(const ConservativeVariable& u,
                                      const ConservativeVariable& F_left,
                                      const ConservativeVariable& F_right,
                                      double x_left, double x_right,
                                      double dt) noexcept;
ConservativeVariable spherical_schema(const ConservativeVariable& u,
                                      const ConservativeVariable& F_left,
                                      const ConservativeVariable& F_right,
                                      double x_left, double x_right,
                                      double dt) noexcept;
ConservativeVariable soft_bundary_condition(const ConservativeVariable& u, double t) noexcept;

class GodunovMethod {
  static const double CFL;

  double t_current;

  std::vector<double> x;
  std::vector<ConservativeVariable> u;
  std::vector<ConservativeVariable> F;

  DifferenceSchema difference_schema;

  BoundaryCondition left_boundary_condition;
  BoundaryCondition right_boundary_condition;

  std::shared_ptr<RiemannSolver> solver;

  void calculateFluxes(double& dt);

 public:
  GodunovMethod(DifferenceSchema difference_schema) noexcept;

  void setLeftBoundaryCondition(BoundaryCondition condition) noexcept;
  void setRightBoundaryCondition(BoundaryCondition condition) noexcept;

  void setSolver(std::shared_ptr<RiemannSolver> solver) noexcept;

  void setInitialConditions(const std::vector<double>& knots_coordinates,
                            const std::vector<GasFlow>& initial_flow) noexcept;

  void run(double t_end);

  std::tuple<std::vector<double>, std::vector<GasFlow>> getSolution()
      const noexcept;
  double getCurrentTime() const noexcept;
};
