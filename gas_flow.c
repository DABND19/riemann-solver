#include "gas_flow.h"

#include <gsl/gsl_math.h>
#include <math.h>

double sound_speed(const GasFlow* flow) {
  return sqrt(GAMMA * flow->pressure / flow->density);
}

double internal_energy(const GasFlow* flow) {
  return flow->pressure / ((GAMMA - 1.) * flow->density);
}

double total_energy(const GasFlow* flow) {
  return internal_energy(flow) + 0.5 * gsl_pow_2(flow->velocity);
}
