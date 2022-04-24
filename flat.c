#include <gsl/gsl_math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "riemann_problem.h"
#include "godunov.h"

static const double EPSILON = 1e-5;

int main(int argc, char** argv) {
  if (argc != 8) {
    fprintf(stderr, "Required 7 parameters.\n");
    return 1;
  }

  const GasFlow LEFT = {atof(argv[1]), atof(argv[2]), atof(argv[3])};
  const GasFlow RIGHT = {atof(argv[4]), atof(argv[5]), atof(argv[6])};
  double T = atof(argv[7]);

  // printf("Godunov %f %f %f %f %f %f %f\n", T,
  //        LEFT.pressure, LEFT.density, LEFT.velocity,
  //        RIGHT.pressure, RIGHT.density, RIGHT.velocity);

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
  RiemannSolver* riemann_solver = riemann_solver_alloc();
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

    double max_velocity = 0.;
    for (size_t i = 0; i < CELLS_NUM - 1; ++i) {
      GasFlow left_flow = from_conservative(u1[i], u2[i], u3[i]);
      GasFlow right_flow = from_conservative(u1[i + 1], u2[i + 1], u3[i + 1]);

      int error = riemann_solver_set(riemann_solver, &left_flow, &right_flow);
      if (error) {
        fprintf(stderr, "Riemann solver error:\n");
        fprintf(stderr, "Left: %f\t%f\t%f\n", 
                left_flow.pressure, left_flow.density, left_flow.velocity);
        fprintf(stderr, "Right: %f\t%f\t%f\n", 
                right_flow.pressure, right_flow.density, right_flow.velocity);
        return 1;
      }

      double knot_velocity = 0.;
      GasFlow knot_flow = riemann_problem_solution(riemann_solver, knot_velocity);
      f1[i + 1] = F1(&knot_flow);
      f2[i + 1] = F2(&knot_flow);
      f3[i + 1] = F3(&knot_flow);

      double S_l = riemann_problem_left_wave_velocity(riemann_solver);
      double S_r = riemann_problem_right_wave_velocity(riemann_solver);
      max_velocity = gsl_max(max_velocity, fabs(S_l));
      max_velocity = gsl_max(max_velocity, fabs(S_r));
    }

    double dt = T - t;
    if (gsl_fcmp(max_velocity, 0., EPSILON) > 0) {
      dt = gsl_min(dt, CFL * 0.5 * dx / max_velocity);
    }

    for (size_t i = 0; i < CELLS_NUM; ++i) {
      double x_left = X_LEFT + dx * i;
      double x_right = X_LEFT + dx * (i + 1);
      u1[i] = flat_schema(u1[i], f1[i], f1[i + 1], x_left, x_right, dt);
      u2[i] = flat_schema(u2[i], f2[i], f2[i + 1], x_left, x_right, dt);
      u3[i] = flat_schema(u3[i], f3[i], f3[i + 1], x_left, x_right, dt);
    }

    t += dt;
  }

  int error = riemann_solver_set(riemann_solver, &LEFT, &RIGHT);
  if (error) {
    fprintf(stderr, "Riemann solver error:\n");
    fprintf(stderr, "Left: %f\t%f\t%f\n", 
            LEFT.pressure, LEFT.density, LEFT.velocity);
    fprintf(stderr, "Right: %f\t%f\t%f\n", 
            RIGHT.pressure, RIGHT.density, RIGHT.velocity);
    return 1;
  }
  for (size_t i = 0; i < CELLS_NUM; ++i) {
    double x_left = X_LEFT + i * dx;
    double x_right = X_LEFT + (i + 1) * dx;
    double x = 0.5 * (x_left + x_right);
    GasFlow solution = from_conservative(u1[i], u2[i], u3[i]);
    GasFlow exact_solution = riemann_problem_solution(riemann_solver, x / T);
    printf("%f\t%f\t%f\n", x, solution.density, exact_solution.density);
  }

  riemann_solver_free(riemann_solver);
  free(u1);
  free(u2);
  free(u3);
  free(f1);
  free(f2);
  free(f3);
  return 0;
}
