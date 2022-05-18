#include <gsl/gsl_math.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "gas_flow.hpp"
#include "godunov.hpp"
#include "riemann_solvers/hll.hpp"
#include "riemann_solvers/hllc.hpp"
#include "riemann_solvers/exact.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "You must specify 1 positional argument" << std::endl;
  }

  const GasFlow LEFT = {10.3333, 3.857143, 2.629369};
  const double T = std::atof(argv[1]);

  const size_t CELLS_COUNT = 10000;
  const double X_LEFT = 0, X_RIGHT = 1.;
  std::vector<double> x_knots(CELLS_COUNT + 1);
  for (size_t i = 0; i < CELLS_COUNT + 1; ++i) {
    double dx = (X_RIGHT - X_LEFT) / CELLS_COUNT;
    x_knots[i] = X_LEFT + dx * i;
  }
  std::vector<GasFlow> initial_flow(CELLS_COUNT);
  for (size_t i = 0; i < CELLS_COUNT; ++i) {
    double x = 0.5 * (x_knots[i + 1] + x_knots[i]);
    if (x < 0.125) {
      initial_flow[i] = LEFT;
    } else {
      initial_flow[i] = {1., 1. + 0.2 * std::sin(2. * M_PI * 8. * x), 0.};
    }
  }

  auto godunov = std::make_unique<GodunovMethod>(cartesian_schema);
  godunov->setInitialConditions(x_knots, initial_flow);

  auto left_boundary_condition = [LEFT](const ConservativeVariable&, double t) {
    return to_conservative(LEFT);
  };
  godunov->setLeftBoundaryCondition(left_boundary_condition);
  godunov->setRightBoundaryCondition(soft_bundary_condition);

  auto riemann_solver = std::shared_ptr<RiemannSolver>(new ExactSolver());
  godunov->setSolver(riemann_solver);

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
