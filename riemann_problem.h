#pragma once

static const double EPSILON = 1e-6;
static const double GAMMA = 5. / 3.;
static const int MAX_ITERATIONS_COUNT = 1000;

typedef struct __gas_parameters {
  double pressure;
  double density;
  double velocity;
} GasParameters;

double speed_of_sound(const GasParameters* flow);

typedef struct __riemann_problem_solution {
  const GasParameters left;
  const GasParameters right;
  const double contact_pressure;
} RiemannProblemSolution;

RiemannProblemSolution solve_riemann_problem(const GasParameters* left,
                                             const GasParameters* right,
                                             int* error_code);
double contact_surface_velocity(const RiemannProblemSolution* solution);
double left_wave_velocity(const RiemannProblemSolution* solution);
double right_wave_velocity(const RiemannProblemSolution* solution);
GasParameters solution_on_curve(const RiemannProblemSolution* solution,
                                double curve);
