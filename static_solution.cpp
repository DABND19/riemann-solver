#include <iomanip>
#include <iostream>
#include <memory>

#include "godunov.hpp"
#include "riemann_solvers/exact.hpp"
#include "stellar_wind.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "You must pass 3 required argumets: time, pressure, Mach number" << std::endl;
    return 1;
  }

  const double T = std::atof(argv[1]);
  const double K = std::atof(argv[2]);
  const double M_e = std::atof(argv[3]);
  const double X_LEFT = R_e(K);
  const double X_RIGHT = 3;
  const size_t CELLS_COUNT = 100;
  const size_t KNOTS_COUNT = CELLS_COUNT + 1;

  auto godunov = std::make_unique<GodunovMethod>(spherical_schema);

  auto exact_solver = std::make_shared<ExactSolver>();
  godunov->setSolver(std::static_pointer_cast<RiemannSolver>(exact_solver));

  auto knots = gen_mesh(X_LEFT, X_RIGHT, KNOTS_COUNT);
  auto initial = gen_initial_state(knots, [K, M_e](double x) {
    if (x < 1.) {
      return supersonic_flow(R_e(K), M_e)(x);
    } else {
      return subsonic_flow(K)(x);
    }
  });
  godunov->setInitialConditions(knots, initial);

  double dx = (X_RIGHT - X_LEFT) / CELLS_COUNT;
  godunov->setLeftBoundaryCondition([X_LEFT, M_e, dx](const ConservativeVariable&, double) {
    auto flow = supersonic_flow(X_LEFT, M_e)(X_LEFT + 0.3 * dx);
    return to_conservative(flow);
  });
  godunov->setRightBoundaryCondition([dx, K, X_RIGHT](const ConservativeVariable& u, double) {
    auto source = subsonic_flow(K)(X_RIGHT);
    // auto flow = from_conservative(u);
    return to_conservative(source);
  });

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