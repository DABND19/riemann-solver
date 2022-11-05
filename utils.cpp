#include "utils.hpp"

std::vector<double> gen_mesh(double x_left, double x_right,
                             size_t knots_count) noexcept {
  double dx = (x_right - x_left) / (knots_count - 1);
  std::vector<double> mesh(knots_count);
  for (size_t i = 0; i < knots_count; ++i) {
    mesh[i] = x_left + dx * i;
  }
  return mesh;
}

std::vector<GasFlow> gen_initial_state(
    const std::vector<double> knots,
    std::function<GasFlow(double)> generator) noexcept {
  std::vector<GasFlow> initial_state(knots.size() - 1);
  for (size_t i = 0; i < initial_state.size(); ++i) {
    double x = 0.5 * (knots[i] + knots[i + 1]);
    initial_state[i] = generator(x);
  }
  return initial_state;
}
