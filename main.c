#include <stdio.h>
#include <stdlib.h>

#include "riemann_problem.h"

int main(int argc, char** argv) {
  GasParameters left = {atof(argv[1]), atof(argv[2]), atof(argv[3])};
  GasParameters right = {atof(argv[4]), atof(argv[5]), atof(argv[6])};
  int status;
  RiemannProblemSolution solution =
      solve_riemann_problem(&left, &right, &status);
  // printf("Contact pressure: %f\n", solution.contact_pressure);
  // printf("Status code: %d\n", status);

  printf("Exact solution %f %f %f %f %f %f\n", 
         left.pressure, left.density, left.velocity,
         right.pressure, right.density, right.velocity);

  double D_left = left_wave_velocity(&solution);
  double D_right = right_wave_velocity(&solution);
  // double U_contact = contact_surface_velocity(&solution);
  double step = (D_right - D_left + 2) / 1000;
  for (float curve = D_left - 1; curve < D_right + 1; curve += step) {
    GasParameters flow = solution_on_curve(&solution, curve);
    printf("%f\t%f\n", curve, flow.density);
  }
  return 0;
}
