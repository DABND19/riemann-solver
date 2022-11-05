#include "godunov.hpp"

#include <cmath>

const double GodunovMethod::CFL = 0.2;

ConservativeVariable cartesian_schema(const ConservativeVariable& u,
                                      const ConservativeVariable& F_left,
                                      const ConservativeVariable& F_right,
                                      double x_left, double x_right,
                                      double dt) noexcept {
  double dx = x_right - x_left;
  return u - dt / dx * (F_right - F_left);
}

ConservativeVariable spherical_schema(const ConservativeVariable& u,
                                      const ConservativeVariable& F_left,
                                      const ConservativeVariable& F_right,
                                      double x_left, double x_right,
                                      double dt) noexcept {
  double dx = x_right - x_left;
  double x = 0.5 * (x_left + x_right);

  auto flow = from_conservative(u);
  auto Q = to_flux(flow);
  auto q = std::make_tuple(0., flow.getPressure(), 0.);

  return u - dt / dx * (F_right - F_left) + 2. / x * dt * (q - Q);
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

GasFlow left_approximation(const GasFlow& flow, double r_left, double r_right) {
  if (flow.getMachNumber() < 1) {
    return flow;
  }
  double r_mid = 0.5 * (r_left + r_right);
  return GasFlow(
    flow.getPressure() * std::pow(r_mid / r_right, 2. * GasFlow::GAMMA), 
    flow.getDensity() * gsl_pow_2(r_mid / r_right), 
    flow.getVelocity());
}

GasFlow right_approximation(const GasFlow& flow, double r_left, double r_right) {
  if (flow.getMachNumber() < 1) {
    return flow;
  }
  double r_mid = 0.5 * (r_left + r_right);
  return GasFlow(
    flow.getPressure() * std::pow(r_mid / r_left, 2. * GasFlow::GAMMA), 
    flow.getDensity() * gsl_pow_2(r_mid / r_left), 
    flow.getVelocity());
}

void GodunovMethod::calculateFluxes(double& dt) {
  for (size_t i = 0; i < this->F.size(); ++i) {
    auto left_cell = i > 0 
        ? u[i - 1]
        : this->left_boundary_condition(u.front(), this->getCurrentTime());
    auto right_cell = i < this->u.size()
        ? u[i]
        : this->right_boundary_condition(u.back(), this->getCurrentTime());

    auto left_flow = i > 0
        ? left_approximation(from_conservative(left_cell), this->x[i], this->x[i + 1])
        : left_approximation(from_conservative(left_cell), 
                             this->x[i] - 0.5 * (this->x[i + i] - this->x[i]), this->x[i]);
    auto right_flow = i < this->u.size()
        ? right_approximation(from_conservative(right_cell), this->x[i + 1], 
                              this->x[i + 1] + 0.5 * (this->x[i + 1] - this->x[i]))
        : from_conservative(right_cell);

    this->solver->setState(left_flow, right_flow);
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
