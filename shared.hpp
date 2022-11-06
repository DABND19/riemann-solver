#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <nlohmann/json.hpp>

#include "gas_flow.hpp"

using json = nlohmann::json;

void to_json(json& payload, const GasFlow& flow) {
  payload = json{
    {"pressure", flow.getPressure()},
    {"density", flow.getDensity()},
    {"velocity", flow.getVelocity()},
    {"Mach", flow.getMachNumber()},
  };
}

void linspace(std::vector<double>& v, double left, double right) {
  double step = (right - left) / (v.size() - 1);
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = left + step * i;
  }
}
