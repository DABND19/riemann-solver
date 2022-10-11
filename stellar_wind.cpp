#include "stellar_wind.hpp"

#include <gsl/gsl_math.h>

#include <cmath>

double R_e(double K) noexcept {
  return std::sqrt(K * 2. * (GasFlow::GAMMA + 1.) / (GasFlow::GAMMA + 3.));
}

std::function<GasFlow(double)> supersonic_flow(double r_e,
                                               double M_e) noexcept {
  return [r_e, M_e](double x) {
    double r = 1. * gsl_pow_2(r_e / x);
    double u = 1.;
    double p = 1. / (GasFlow::GAMMA * gsl_pow_2(M_e)) *
               std::pow(r_e / x, 2. * GasFlow::GAMMA);
    return GasFlow(p, r, u);
  };
}

std::function<GasFlow(double)> subsonic_flow(double K) noexcept {
  return [K](double x) {
    double r_e = R_e(K);
    double r = (GasFlow::GAMMA + 1.) / (GasFlow::GAMMA - 1.) * gsl_pow_2(r_e);
    double u =
        (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.) * gsl_pow_2(1. / x);
    double p = K - 0.5 * (GasFlow::GAMMA - 1.) / (GasFlow::GAMMA + 1.) *
                       gsl_pow_2(r_e) / gsl_pow_4(x);
    return GasFlow(p, r, u);
  };
}
