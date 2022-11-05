#pragma once

#include <functional>
#include <vector>

#include "gas_flow.hpp"

std::vector<double> gen_mesh(double x_left, double x_right,
                             size_t knots_count) noexcept;
std::vector<GasFlow> gen_initial_state(
    const std::vector<double> knots,
    std::function<GasFlow(double)> generator) noexcept;
