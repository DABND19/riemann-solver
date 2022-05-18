#include "gas_flow.hpp"

const double GasFlow::GAMMA = 1.4;

GasFlow::GasFlow(double pressure, double density, double velocity) noexcept
    : pressure(pressure),
      density(density),
      velocity(velocity),
      sound_speed(sqrt(GAMMA * pressure / density)) {}

double GasFlow::getInternalEnergy() const noexcept {
  return this->pressure / ((this->GAMMA - 1.) * this->density);
}

double GasFlow::getTotalEnergy() const noexcept {
  return this->getInternalEnergy() + 0.5 * gsl_pow_2(this->velocity);
}

double GasFlow::getMachNumber() const noexcept {
  return fabs(this->velocity) / this->sound_speed;
}
