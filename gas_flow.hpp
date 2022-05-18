#pragma once

#include <gsl/gsl_math.h>

#include <cmath>

class GasFlow {
  double pressure;
  double density;
  double velocity;
  double sound_speed;

 public:
  static const double GAMMA;

  GasFlow() = default;
  GasFlow(double pressure, double density, double velocity) noexcept;

  double getInternalEnergy() const noexcept;
  double getTotalEnergy() const noexcept;
  double getMachNumber() const noexcept;

  double getPressure() const noexcept { return this->pressure; }
  double getDensity() const noexcept { return this->density; }
  double getVelocity() const noexcept { return this->velocity; }
  double getSoundSpeed() const noexcept { return this->sound_speed; }
};
