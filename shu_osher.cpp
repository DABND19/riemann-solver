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
#include "utils.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "You must specify 1 positional argument" << std::endl;
    return 1;
  }

  const GasFlow LEFT = {10.3333, 3.857143, 2.629369};
  const double T = std::atof(argv[1]);

  const size_t CELLS_COUNT = 500;
  const double X_LEFT = 0, X_RIGHT = 1.;

  auto godunov = std::make_unique<GodunovMethod>(cartesian_schema);

  auto knots = gen_mesh(X_LEFT, X_RIGHT, CELLS_COUNT + 1);
  auto initial = gen_initial_state(knots, [LEFT](double x) {
    if (x < 0.125) {
      return LEFT;
    } else {
      return GasFlow(1., 1. + 0.2 * std::sin(16 * M_PI * x), 0.);
    }
  });
  godunov->setInitialConditions(knots, initial);

  auto left_boundary_condition = [LEFT](const ConservativeVariable&, double t) {
    return to_conservative(LEFT);
  };
  godunov->setLeftBoundaryCondition(left_boundary_condition);
  godunov->setRightBoundaryCondition(soft_bundary_condition);

  auto exact_solver = std::make_shared<ExactSolver>();
  auto hll_solver = std::make_shared<HllSolver>();
  auto hllc_solver = std::make_shared<HllcSolver>();
  godunov->setSolver(std::static_pointer_cast<RiemannSolver>(hll_solver));

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
