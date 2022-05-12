#include "godunov.hpp"

#include <cmath>

const double GodunovMethod::CFL = 0.2;

ConservativeVariable cartesian_schema(const ConservativeVariable& u,
                                      const ConservativeVariable& F_left,
                                      const ConservativeVariable& F_right,
                                      double x_left, double x_right,
                                      double dt) noexcept {
  double dx = x_right - x_left;
  auto [u1, u2, u3] = u;
  auto [F1_left, F2_left, F3_left] = F_left;
  auto [F1_right, F2_right, F3_right] = F_right;

  u1 = u1 - dt / dx * (F1_right - F1_left);
  u2 = u2 - dt / dx * (F2_right - F2_left);
  u3 = u3 - dt / dx * (F3_right - F3_left);
  return std::make_tuple(u1, u2, u3);
}

ConservativeVariable spherical_schema(const ConservativeVariable& u,
                                      const ConservativeVariable& F_left,
                                      const ConservativeVariable& F_right,
                                      double x_left, double x_right,
                                      double dt) noexcept {
  double dx = x_right - x_left;
  double x = 0.5 * (x_left + x_right);

  auto flow = from_conservative(u);
  auto [Q1, Q2, Q3] = to_flux(flow);
  double q1 = 0.;
  double q2 = flow.getPressure();
  double q3 = 0.;

  auto [u1, u2, u3] = u;
  auto [F1_left, F2_left, F3_left] = F_left;
  auto [F1_right, F2_right, F3_right] = F_right;

  u1 = u1 - dt / dx * (F1_right - F1_left) + 2. / x * dt * (q1 - Q1);
  u2 = u2 - dt / dx * (F2_right - F2_left) + 2. / x * dt * (q2 - Q2);
  u3 = u3 - dt / dx * (F3_right - F3_left) + 2. / x * dt * (q3 - Q3);
  return std::make_tuple(u1, u2, u3);
}

ConservativeVariable soft_bundary_condition(const ConservativeVariable& u, double t) noexcept {
  return u;
}

GodunovMethod::GodunovMethod(DifferenceSchema difference_schema) noexcept
    : difference_schema(difference_schema) {}

void GodunovMethod::setLeftBoundaryCondition(
    BoundaryCondition condition) noexcept {
  this->left_boundary_condition = condition;
}

void GodunovMethod::setRightBoundaryCondition(
    BoundaryCondition condition) noexcept {
  this->right_boundary_condition = condition;
}

void GodunovMethod::setSolver(std::shared_ptr<RiemannSolver> solver) noexcept {
  this->solver = solver;
}

void GodunovMethod::setInitialConditions(
    const std::vector<double>& knots_coordinates,
    const std::vector<GasFlow>& initial_flow) noexcept {
  this->t_current = 0.;

  this->x = knots_coordinates;
  this->F = std::vector<ConservativeVariable>(x.size());

  this->u = std::vector<ConservativeVariable>(initial_flow.size());
  for (size_t i = 0; i < this->u.size(); ++i) {
    this->u[i] = to_conservative(initial_flow[i]);
  }
}

double GodunovMethod::getCurrentTime() const noexcept {
  return this->t_current;
}

void GodunovMethod::calculateFluxes(double& dt) {
  for (size_t i = 0; i < this->F.size(); ++i) {
    auto left_cell = i > 0 
        ? u[i - 1]
        : this->left_boundary_condition(u.front(), this->getCurrentTime());
    auto right_cell = i < this->u.size()
        ? u[i]
        : this->right_boundary_condition(u.back(), this->getCurrentTime());

    this->solver->setState(from_conservative(left_cell), from_conservative(right_cell));
    this->F[i] = solver->getFlux();

    auto [left_wave, right_wave] = solver->getWaveSpeed();
    if (i > 0 && gsl_fcmp(left_wave, 0., 1e-16)) {
      double dx = x[i] - x[i - 1];
      dt = gsl_min(dt, this->CFL * dx / std::fabs(left_wave));
    }
    if (i < this->F.size() - 1 && gsl_fcmp(right_wave, 0., 1e-16)) {
      double dx = x[i + 1] - x[i];
      dt = gsl_min(dt, this->CFL * dx / std::fabs(right_wave));
    }
  }
}

void GodunovMethod::run(double t_end) {
  while (this->t_current < t_end) {
    double dt = t_end - this->t_current;
    this->calculateFluxes(dt);

    for (size_t i = 0; i < this->u.size(); ++i) {
      this->u[i] = this->difference_schema(this->u[i], 
                                           this->F[i], this->F[i + 1], 
                                           this->x[i], this->x[i + 1], dt);
    }

    this->t_current += dt;
  }
}

std::tuple<std::vector<double>, std::vector<GasFlow>> GodunovMethod::getSolution() 
    const noexcept {
  std::vector<GasFlow> flow(this->u.size());
  for (size_t i = 0; i < flow.size(); ++i) {
    flow[i] = from_conservative(this->u[i]);
  }
  return std::make_pair(this->x, std::move(flow));
}