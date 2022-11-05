#include <iomanip>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <vector>

#include "gas_flow.hpp"
#include "godunov.hpp"
#include "riemann_solvers/exact.hpp"
#include "riemann_solvers/hll.hpp"
#include "riemann_solvers/hllc.hpp"
#include "stellar_wind.hpp"
#include "utils.hpp"

nlohmann::json dump_solution(double t, const std::vector<double>& knots,
                             const std::vector<GasFlow>& solution) {
  std::vector<double> x(solution.size());
  std::vector<double> density(solution.size());
  std::vector<double> pressure(solution.size());
  for (size_t i = 0; i < solution.size(); ++i) {
    x[i] = 0.5 * (knots[i] + knots[i + 1]);
    density[i] = solution[i].getDensity();
    pressure[i] = solution[i].getPressure();
  }
  return nlohmann::json::object(
      {{"t", t}, {"x", x}, {"density", density}, {"pressure", pressure}});
}

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "You must pass 5 required args: (T, M_e, K, St, A)"
              << std::endl;
    std::exit(1);
  }

  const double T = std::atof(argv[1]);
  const double M_e = std::atof(argv[2]);
  const double K = std::atof(argv[3]);
  const double St = std::atof(argv[4]);
  const double A = std::atof(argv[5]);

  auto exact_solver = std::make_shared<ExactSolver>();
  auto hll_solver = std::make_shared<HllSolver>();
  auto hllc_solver = std::make_shared<HllcSolver>();

  auto godunov = std::make_unique<GodunovMethod>(spherical_schema);
  godunov->setSolver(std::static_pointer_cast<RiemannSolver>(exact_solver));

  const double X_LEFT = R_e(K);
  const double X_RIGHT = 3.;
  const size_t KNOTS_COUNT = 1000 + 1;
  auto knots = gen_mesh(X_LEFT, X_RIGHT, KNOTS_COUNT);
  auto initial = gen_initial_state(knots, [K, M_e](double x) {
    if (x < 1.) {
      return supersonic_flow(R_e(K), M_e)(x);
    } else {
      return subsonic_flow(K)(x);
    }
  });
  godunov->setInitialConditions(knots, initial);

  double dx = (X_RIGHT - X_LEFT) / (KNOTS_COUNT - 1);
  // godunov->setLeftBoundaryCondition([dx, X_LEFT, M_e](const
  // ConservativeVariable&, double) {
  //   return to_conservative(supersonic_flow(X_LEFT, M_e)(X_LEFT + dx));
  // });
  godunov->setRightBoundaryCondition(
      [dx, K, X_RIGHT](const ConservativeVariable& u, double) {
        auto source = subsonic_flow(K)(X_RIGHT);
        auto flow = from_conservative(u);
        double pressure = K;
        double density = flow.getDensity();
        double velocity = flow.getVelocity();
        return to_conservative({pressure, density, velocity});
      });
  godunov->setLeftBoundaryCondition(
      [dx, X_LEFT, M_e, A, St](const ConservativeVariable&, double t) {
        auto flow = supersonic_flow(X_LEFT, M_e)(X_LEFT + dx);
        double p = flow.getPressure();
        double r = flow.getDensity() * (1. + A * std::sin(t / St));
        double u = flow.getVelocity();
        return to_conservative({p, r, u});
      });

  auto data = nlohmann::json::array();
  const size_t TIMES_COUNT = 1000;
  for (size_t i = 0; i < TIMES_COUNT; ++i) {
    double dt = T / TIMES_COUNT;
    godunov->run(dt * i);
    auto [x, sol] = godunov->getSolution();
    data.push_back(dump_solution(godunov->getCurrentTime(), x, sol));
  }

  std::cout << std::scientific << data << std::endl;

  return 0;
}
