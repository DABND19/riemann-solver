#include <gsl/gsl_math.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "gas_flow.hpp"
#include "godunov.hpp"
#include "riemann_solvers/hllc.hpp"

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "You must specify 4 positional arguments" << std::endl;
  }

  double T = std::atof(argv[1]);
  const auto LEFT =
      GasFlow(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));

  const size_t CELLS_COUNT = 300;
  const double X_LEFT = 0.05, X_RIGHT = 3.;
  std::vector<double> x_knots(CELLS_COUNT + 1);
  for (size_t i = 0; i < CELLS_COUNT + 1; ++i) {
    double dx = (X_RIGHT - X_LEFT) / CELLS_COUNT;
    x_knots[i] = X_LEFT + dx * i;
  }
  std::vector<GasFlow> initial_flow(CELLS_COUNT);
  for (size_t i = 0; i < CELLS_COUNT; ++i) {
    double x = 0.5 * (x_knots[i] + x_knots[i + 1]);
    double density = LEFT.getDensity() * gsl_pow_2(X_LEFT / x);
    double velocity = LEFT.getVelocity();
    double pressure =
        LEFT.getPressure() * std::pow(X_LEFT / x, 2. * GasFlow::GAMMA);
    initial_flow[i] = GasFlow(pressure, density, velocity);
  }

  auto godunov = std::make_unique<GodunovMethod>(spherical_schema);
  godunov->setInitialConditions(x_knots, initial_flow);

  auto left_boundary_condition = [LEFT](const ConservativeVariable&, double t) {
    double pressure = LEFT.getPressure();
    double density = LEFT.getDensity();
    double velocity = LEFT.getVelocity();
    return to_conservative({pressure, density, velocity});
  };
  godunov->setLeftBoundaryCondition(left_boundary_condition);
  godunov->setRightBoundaryCondition(soft_bundary_condition);

  auto hllc_solver = std::shared_ptr<RiemannSolver>(new HllcSolver());
  godunov->setSolver(hllc_solver);

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
