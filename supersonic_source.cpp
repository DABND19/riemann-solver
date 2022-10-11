#include <gsl/gsl_math.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

#include "gas_flow.hpp"
#include "godunov.hpp"
#include "riemann_solvers/exact.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "You must specify 1 required argument: time" << std::endl;
    return 1;
  }

  const double T = std::atof(argv[1]);
  const double X_LEFT = 0.05;
  const double X_RIGHT = 3.;
  const size_t CELLS_COUNT = 10000;
  const double dx = (X_RIGHT - X_LEFT) / CELLS_COUNT;

  auto godunov = std::make_unique<GodunovMethod>(spherical_schema);

  auto exact_solver = std::make_shared<ExactSolver>();
  godunov->setSolver(std::static_pointer_cast<RiemannSolver>(exact_solver));

  auto knots = gen_mesh(X_LEFT, X_RIGHT, CELLS_COUNT + 1);
  auto initial = gen_initial_state(knots, [X_LEFT](double x) {
    double density = gsl_pow_2(X_LEFT / x);
    double velocity = 1;
    double pressure = 0.006 * std::pow(X_LEFT / x, 2 * GasFlow::GAMMA);
    return GasFlow(pressure, density, velocity);
  });
  godunov->setInitialConditions(knots, initial);

  godunov->setLeftBoundaryCondition([X_LEFT, dx](const ConservativeVariable&, double t) {
    double x = X_LEFT + dx;
    double density = gsl_pow_2(X_LEFT / x);
    double velocity = 1;
    double pressure = 0.006 * std::pow(X_LEFT / x, 2 * GasFlow::GAMMA);
    return to_conservative(GasFlow(pressure, density, velocity));
  });
  godunov->setRightBoundaryCondition(soft_bundary_condition);

  godunov->run(T);

  auto [x, flow] = godunov->getSolution();
  for (size_t i = 0; i < flow.size(); ++i) {
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 +
                                   1)
              << std::scientific << 0.5 * (x[i] + x[i + 1]) << '\t'
              << flow[i].getPressure() << '\t' << flow[i].getDensity() << '\t'
              << flow[i].getVelocity() << '\t' << flow[i].getMachNumber()
              << '\n';
  }
  return 0;
}
