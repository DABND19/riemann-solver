#pragma once

#include <functional>

#include "gas_flow.hpp"

double R_e(double K) noexcept;
std::function<GasFlow(double)> supersonic_flow(double r_e,
                                               double M_e) noexcept;
std::function<GasFlow(double)> subsonic_flow(double K) noexcept;
