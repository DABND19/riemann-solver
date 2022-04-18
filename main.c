#include <stdio.h>
#include <stdlib.h>

#include "riemann_problem.h"

void print_gas_parameters(const GasParameters* flow) {
  printf("pressure=%f, density=%f, velocity=%f", flow->pressure, flow->density,
         flow->velocity);
}

int main(int argc, char** argv) {
  if (argc < 7) {
    printf("Please specify left and right parameters");
    return 1;
  }

  GasParameters left = {atof(argv[1]), atof(argv[2]), atof(argv[3])};
  printf("Left parameters: ");
  print_gas_parameters(&left);
  printf("\n");
  GasParameters right = {atof(argv[4]), atof(argv[5]), atof(argv[6])};
  printf("Right parameters: ");
  print_gas_parameters(&right);
  printf("\n");
  int status;
  RiemannProblemSolution solution =
      solve_riemann_problem(&left, &right, &status);
  printf("Contact pressure: %f\n", solution.contact_pressure);
  printf("Status code: %d\n", status);

  double D_left = left_wave_velocity(&solution);
  double D_right = right_wave_velocity(&solution);
  double U_contact = contact_surface_velocity(&solution);
  printf("D_left=%f U_contact=%f D_right=%f\n", D_left, U_contact, D_right);
  double step = (D_right - D_left + 2) / 100.;
  for (float curve = D_left - 1; curve < D_right + 1; curve += step) {
    GasParameters params = solution_on_curve(&solution, curve);
    printf("curve=%f ", curve);
    print_gas_parameters(&params);
    printf("\n");
  }
  return 0;
}
