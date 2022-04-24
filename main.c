#include <gsl/gsl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "godunov.h"
#include "riemann_problem.h"

static const TOL = 1e-5;

int main(int argc, char* argv[]) {
  const double T = 1.;
  const double X_LEFT = 0.;
  const double X_RIGHT = 1.;

  const size_t CELLS_NUM = 1000;
  const double dx = (X_RIGHT - X_LEFT) / CELLS_NUM;

  double* u1 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* u2 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* u3 = (double*)malloc(CELLS_NUM * sizeof(double));

  for (size_t i = 0; i < CELLS_NUM; ++i) {
    // Initial conditions...
  }

  const size_t KNOTS_NUM = CELLS_NUM + 1;
  double* f1 = (double*)malloc(KNOTS_NUM * sizeof(double));
  double* f2 = (double*)malloc(KNOTS_NUM * sizeof(double));
  double* f3 = (double*)malloc(KNOTS_NUM * sizeof(double));

  double* q1 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* q2 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* q3 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* Q1 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* Q2 = (double*)malloc(CELLS_NUM * sizeof(double));
  double* Q3 = (double*)malloc(CELLS_NUM * sizeof(double));

  RiemannSolver* solver = riemann_solver_alloc();

  double t = 0.;
  while (t < T) {
    double dt = T - t;

    for (size_t i = 0; i < CELLS_NUM; ++i) {
      GasFlow cell_flow = from_conservative(u1[i], u2[i], u3[i]);

      Q1[i] = F1(&cell_flow);
      Q2[i] = F2(&cell_flow);
      Q3[i] = F3(&cell_flow);

      q1[i] = b1(&cell_flow);
      q2[i] = b2(&cell_flow);
      q3[i] = b3(&cell_flow);
    }

    // TODO: Add left boundary condition
    f1[0] = 0.;
    f2[0] = 0.;
    f3[0] = 0.;

    GasFlow right_bound = from_conservative(u1[CELLS_NUM - 1], 
                                            u2[CELLS_NUM - 1],
                                            u3[CELLS_NUM - 1]);
    f1[KNOTS_NUM - 1] = F1(&right_bound);
    f2[KNOTS_NUM - 1] = F2(&right_bound);
    f3[KNOTS_NUM - 1] = F3(&right_bound);

    // Calculate fluxes between cells
    double max_velocity = 0.;
    for (size_t i = 1; i < KNOTS_NUM - 1; ++i) {
      GasFlow left_cell = from_conservative(u1[i - 1], u2[i - 1], u3[i - 1]);
      GasFlow right_cell = from_conservative(u1[i], u2[i], u3[i]);

      int error = riemann_solver_set(solver, &left_cell, &right_cell);
      if (error) {
        fprintf(stderr, "Riemann solver error:\n");
        fprintf(stderr, "Left: %f\t%f\t%f\n", 
                left_cell.pressure, left_cell.density, left_cell.velocity);
        fprintf(stderr, "Right: %f\t%f\t%f\n", 
                right_cell.pressure, right_cell.density, right_cell.velocity);
        return 1;
      }

      double knot_speed = 0.;
      GasFlow knot_flow = riemann_problem_solution(solver, knot_speed);
      f1[i] = F1(&knot_flow);
      f2[i] = F2(&knot_flow);
      f3[i] = F3(&knot_flow);

      double left_wave_speed = riemann_problem_left_wave_velocity(solver);
      double right_wave_speed = riemann_problem_right_wave_velocity(solver);
      max_velocity = gsl_max(max_velocity, fabs(left_wave_speed));
      max_velocity = gsl_max(max_velocity, fabs(right_wave_speed));
    }

    if (gsl_fcmp(max_velocity, 0., TOL) > 0) {
      dt = gsl_min(dt, CFL * 0.5 * dx / max_velocity);
    }

    for (size_t i = 0; i < CELLS_NUM; ++i) {
      double x_left = X_LEFT + dx * i;
      double x_right = X_LEFT + dx * (i + 1);

      u1[i] = spherical_schema(u1[i], 
                               f1[i], f1[i + 1], 
                               q1[i], Q1[i], 
                               x_left, x_right,
                               dt);
      u2[i] = spherical_schema(u2[i], 
                               f2[i], f2[i + 1], 
                               q2[i], Q2[i], 
                               x_left, x_right,
                               dt);
      u3[i] = spherical_schema(u3[i], 
                               f3[i], f3[i + 1], 
                               q3[i], Q3[i], 
                               x_left, x_right,
                               dt);
    }

    t += dt;
  }

  riemann_solver_free(solver);

  free(q1);
  free(q2);
  free(q3);
  free(Q1);
  free(Q2);
  free(Q3);

  free(f1);
  free(f2);
  free(f3);

  free(u1);
  free(u2);
  free(u3);

  return 0;
}
