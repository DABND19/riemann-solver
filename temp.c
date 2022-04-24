#include <gsl/gsl_math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "riemann_problem.h"
#include "godunov.h"

bool is_flow_equals(const GasFlow* left, const GasFlow* right) {
  int equals_velocity = !gsl_fcmp(left->velocity, right->velocity, EPSILON);
  int equals_density = !gsl_fcmp(left->density, right->density, EPSILON);
  int equals_pressure = !gsl_fcmp(left->pressure, right->pressure, EPSILON);
  return equals_velocity & equals_density & equals_pressure;
}

int main(int argc, char** argv) {
  if (argc != 8) {
    fprintf(stderr, "Required 7 parameters.\n");
    return 1;
  }

  const GasFlow LEFT = {atof(argv[1]), atof(argv[2]), atof(argv[3])};
  const GasFlow RIGHT = {atof(argv[4]), atof(argv[5]), atof(argv[6])};
  double T = atof(argv[7]);

  printf("Godunov %f %f %f %f %f %f %f\n", T,
         LEFT.pressure, LEFT.density, LEFT.velocity,
         RIGHT.pressure, RIGHT.density, RIGHT.velocity);

  const double X_LEFT = -1;
  const double X_RIGHT = 1;
  const double X_DIAPH = 0;
  const size_t CELLS_NUM = 1000;

  const double dx = (X_RIGHT - X_LEFT) / CELLS_NUM;

  double* u1 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* u2 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* u3 = (double*)malloc(CELLS_NUM * sizeof(double));

  for (size_t i = 0; i < CELLS_NUM; ++i) {
    double x = X_LEFT + i * dx + 0.5 * dx;
    GasFlow flow = x < X_DIAPH ? LEFT : RIGHT;
    u1[i] = U1(&flow);
    u2[i] = U2(&flow);
    u3[i] = U3(&flow);
  }

  const size_t KNOTS_NUM = CELLS_NUM + 1;
  double* f1 = (double*)malloc(KNOTS_NUM * sizeof(double));
  double* f2 = (double*)malloc(KNOTS_NUM * sizeof(double));
  double* f3 = (double*)malloc(KNOTS_NUM * sizeof(double));

  double t = 0;
  while (t < T) {
    GasFlow left_boundary_flow, right_boundary_flow;

    left_boundary_flow = from_conservative(u1[0], u2[0], u3[0]);
    f1[0] = F1(&left_boundary_flow);
    f2[0] = F2(&left_boundary_flow);
    f3[0] = F3(&left_boundary_flow);

    right_boundary_flow = from_conservative(u1[CELLS_NUM - 1], u2[CELLS_NUM - 1], u3[CELLS_NUM - 1]);
    f1[KNOTS_NUM - 1] = F1(&right_boundary_flow);
    f2[KNOTS_NUM - 1] = F2(&right_boundary_flow);
    f3[KNOTS_NUM - 1] = F3(&right_boundary_flow);

    double max_velocity = 0;
    for (size_t i = 0; i < CELLS_NUM - 1; ++i) {
      GasFlow left_flow = from_conservative(u1[i], u2[i], u3[i]);
      GasFlow right_flow = from_conservative(u1[i + 1], u2[i + 1], u3[i + 1]);

      GasFlow knot_flow;
      double knot_velocity = 0;
      if (is_flow_equals(&left_flow, &right_flow)) {
        f1[i + 1] = F1(&left_flow);
        f2[i + 1] = F2(&left_flow);
        f3[i + 1] = F3(&left_flow);
      } else {
        int error;
        RiemannProblemSolution exact_solution = solve_riemann_problem(&left_flow, &right_flow, &error);
        if (error) {
          fprintf(stderr, "Riemann problem solution error.\n");
          break;
        }

        knot_flow = solution_on_curve(&exact_solution, knot_velocity);
        f1[i + 1] = F1(&knot_flow);
        f2[i + 1] = F2(&knot_flow);
        f3[i + 1] = F3(&knot_flow);

        max_velocity = gsl_max(max_velocity, -riemann_problem_left_wave_velocity(&exact_solution));
        max_velocity = gsl_max(max_velocity, riemann_problem_right_wave_velocity(&exact_solution));
      }
    }

    double dt = T - t;
    if (gsl_fcmp(max_velocity, 0., EPSILON) > 0) {
      dt = gsl_min(dt, CFL * dx / max_velocity);
    }

    for (size_t i = 0; i < CELLS_NUM; ++i) {
      double x_left = X_LEFT + i * dx;
      double x_right = x_left + dx;
      u1[i] = calculus_schema(u1[i], f1[i], f1[i + 1], x_left, x_right, dt);
      u2[i] = calculus_schema(u2[i], f2[i], f2[i + 1], x_left, x_right, dt);
      u3[i] = calculus_schema(u3[i], f3[i], f3[i + 1], x_left, x_right, dt);
    }

    t += dt;
  }

  int error;
  RiemannProblemSolution exact = solve_riemann_problem(&LEFT, &RIGHT, &error);
  for (size_t i = 0; i < CELLS_NUM; ++i) {
    double x_left = X_LEFT + i * dx;
    double x_right = x_left + dx;
    GasFlow flow = from_conservative(u1[i], u2[i], u3[i]);
    GasFlow flow_ = solution_on_curve(&exact, 0.5 * (x_left + x_right) / T);
    printf("%f\t%f\t%f\n", 0.5 * (x_left + x_right), flow.density, flow_.density);
  }

  free(u1);
  free(u2);
  free(u3);
  free(f1);
  free(f2);
  free(f3);
  return 0;
}
