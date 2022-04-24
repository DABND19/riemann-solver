#pragma once

static const double GAMMA = 5. / 3.;

typedef struct __gas_flow {
  double pressure;
  double density;
  double velocity;
} GasFlow;

double sound_speed(const GasFlow* flow);
double internal_energy(const GasFlow* flow);
double total_energy(const GasFlow* flow);
