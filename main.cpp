#include <iostream>
#include <vector>

#include <nlohmann/json.hpp>

#include "shared.hpp"

using json = nlohmann::json;

int main(int argc, char** argv) {
  std::vector<double> x(10 + 1);
  linspace(x, 0., 1.);
  std::cout << json(GasFlow{1, 2, 3}) << std::endl;
  return 0;
}
