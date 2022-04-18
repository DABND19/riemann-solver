#include "riemann_problem.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>

double speed_of_sound(const GasParameters* gas) {
  return sqrt(GAMMA * gas->pressure / gas->density);
}

const double G1 = (GAMMA - 1) / (2 * GAMMA);
const double G2 = (GAMMA + 1) / (2 * GAMMA);
const double G3 = (GAMMA - 1) / 2;
const double G4 = (GAMMA + 1) / 2;
const double G5 = (GAMMA - 1) / (GAMMA + 1);

double f(double contact_pressure, const GasParameters* flow) {
  double P = contact_pressure;
  double p = flow->pressure;
  double c = speed_of_sound(flow);
  double pi = P / p;
  double r = flow->density;

  if (P >= p) {
    double nominator = P - p;
    double denominator = r * c * sqrt(G2 * pi + G1);
    return nominator / denominator;
  } else {
    return 1 / G3 * c * (pow(pi, G1) - 1);
  }
}

double df(double contact_pressure, const GasParameters* flow) {
  double P = contact_pressure;
  double p = flow->pressure;
  double c = speed_of_sound(flow);
  double pi = P / p;
  double r = flow->density;

  if (P >= p) {
    double nominator = (GAMMA + 1) * pi + (3 * GAMMA - 1);
    double denominator = 4 * GAMMA * r * c * sqrt(gsl_pow_3(G2 * pi + G1));
    return nominator / denominator;
  } else {
    return 1. / (GAMMA * P) * c * pow(pi, G1);
  }
}

double contact_pressure_equation(double contact_pressure,
                                 const GasParameters* left,
                                 const GasParameters* right) {
  double f_left = f(contact_pressure, left);
  double f_right = f(contact_pressure, right);
  return f_left + f_right;
}

double derivative_of_contact_pressure_equation(double contact_pressure,
                                               const GasParameters* left,
                                               const GasParameters* right) {
  double df_left = df(contact_pressure, left);
  double df_right = df(contact_pressure, right);
  return df_left + df_right;
}

struct EquationParameters {
  const GasParameters* left;
  const GasParameters* right;
};

double F(double x, void* parameters) {
  struct EquationParameters* parsed_parameters =
      (struct EquationParameters*)parameters;
  double delta_U =
      parsed_parameters->left->velocity - parsed_parameters->right->velocity;
  return contact_pressure_equation(x, parsed_parameters->left,
                                   parsed_parameters->right) -
         delta_U;
}

double dF(double x, void* parameters) {
  struct EquationParameters* parsed_parameters =
      (struct EquationParameters*)parameters;
  return derivative_of_contact_pressure_equation(x, parsed_parameters->left,
                                                 parsed_parameters->right);
}

void FdF(double x, void* parameters, double* y, double* dy) {
  *y = F(x, parameters);
  *dy = dF(x, parameters);
}

double mass_velocity(double contact_pressure, const GasParameters* flow) {
  double P = contact_pressure;
  double p = flow->pressure;
  double r = flow->density;
  double c = speed_of_sound(flow);
  if (P >= p) {
    return sqrt(r * (G4 * P + G3 * p));
  } else {
    double pi = P / p;
    return G1 * r * c * (1 - pi) / (1 - pow(pi, G1));
  }
}

double contact_surface_velocity(const RiemannProblemSolution* solution) {
  double a1 = mass_velocity(solution->contact_pressure, &(solution->left));
  double a2 = mass_velocity(solution->contact_pressure, &(solution->right));
  double u1 = solution->left.velocity;
  double u2 = solution->right.velocity;
  double p1 = solution->left.pressure;
  double p2 = solution->right.pressure;

  return (a1 * u1 + a2 * u2 + p1 - p2) / (a1 + a2);
}

RiemannProblemSolution solve_riemann_problem(const GasParameters* left,
                                             const GasParameters* right,
                                             int* error_code) {
  double max_pressure = gsl_max(left->pressure, right->pressure);
  double min_pressure = gsl_min(left->pressure, right->pressure);
  double U_shock = contact_pressure_equation(max_pressure, left, right);
  double U_rarefaction = contact_pressure_equation(min_pressure, left, right);
  double U_vacuum = contact_pressure_equation(0, left, right);
  double delta_U = left->velocity - right->velocity;
  double contact_pressure;
  *error_code = GSL_SUCCESS;

  struct EquationParameters equation_parameters = {left, right};
  if (delta_U >= U_shock) {
    gsl_function_fdf equation;
    equation.f = &F;
    equation.df = &dF;
    equation.fdf = &FdF;
    equation.params = &equation_parameters;

    double contact_pressure_start = max_pressure;
    contact_pressure = contact_pressure_start;

    gsl_root_fdfsolver* solver =
        gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
    gsl_root_fdfsolver_set(solver, &equation, contact_pressure_start);

    *error_code = GSL_CONTINUE;
    for (int iteration = 0;
         iteration < MAX_ITERATIONS_COUNT && *error_code == GSL_CONTINUE;
         ++iteration) {
      *error_code = gsl_root_fdfsolver_iterate(solver);
      if (*error_code != GSL_SUCCESS) {
        break;
      }

      contact_pressure_start = contact_pressure;
      contact_pressure = gsl_root_fdfsolver_root(solver);

      *error_code = gsl_root_test_delta(contact_pressure,
                                        contact_pressure_start, 0, EPSILON);
    }
    gsl_root_fdfsolver_free(solver);
  } else if (delta_U > U_rarefaction) {
    gsl_function equation;
    equation.function = &F;
    equation.params = &equation_parameters;

    double contact_pressure_left = min_pressure;
    double contact_pressure_right = max_pressure;

    gsl_root_fsolver* solver =
        gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
    gsl_root_fsolver_set(solver, &equation, contact_pressure_left,
                         contact_pressure_right);

    *error_code = GSL_CONTINUE;
    for (int iteration = 0;
         iteration < MAX_ITERATIONS_COUNT && *error_code != GSL_SUCCESS;
         ++iteration) {
      *error_code = gsl_root_fsolver_iterate(solver);
      if (*error_code != GSL_SUCCESS) {
        break;
      }

      contact_pressure = gsl_root_fsolver_root(solver);
      contact_pressure_left = gsl_root_fsolver_x_lower(solver);
      contact_pressure_right = gsl_root_fsolver_x_upper(solver);

      *error_code = gsl_root_test_interval(contact_pressure_left,
                                           contact_pressure_right, 0, EPSILON);
    }
    gsl_root_fsolver_free(solver);
  } else if (delta_U > U_vacuum) {
    contact_pressure =
        min_pressure *
        pow((delta_U - U_vacuum) / (U_rarefaction - U_vacuum), 1 / G1);
  } else {
    contact_pressure = 0;
  }

  return (RiemannProblemSolution){*left, *right, contact_pressure};
}

const int LEFT_DIRECTION = -1;
const int RIGHT_DIRECTION = 1;

double left_wave_velocity(const RiemannProblemSolution* solution) {
  double P = solution->contact_pressure;
  double u1 = solution->left.velocity;
  double p1 = solution->left.pressure;
  double r1 = solution->left.density;
  double a1 = mass_velocity(P, &(solution->left));
  if (P >= p1) {
    return u1 - a1 / r1;
  } else {
    double c1 = speed_of_sound(&(solution->left));
    return u1 - c1;
  }
}

double right_wave_velocity(const RiemannProblemSolution* solution) {
  double P = solution->contact_pressure;
  double u2 = solution->right.velocity;
  double p2 = solution->right.pressure;
  double r2 = solution->right.density;
  double a2 = mass_velocity(P, &(solution->right));
  if (P >= p2) {
    return u2 + a2 / r2;
  } else {
    double c2 = speed_of_sound(&(solution->right));
    return u2 + c2;
  }
}

double left_rarefaction_trailing_edge_valocity(
    const RiemannProblemSolution* solution) {
  double c1 = speed_of_sound(&(solution->left));
  double u1 = solution->left.velocity;
  double U = contact_surface_velocity(solution);
  double c1_star = c1 + G3 * (u1 - U);
  return U - c1_star;
}

double right_rarefaction_trailing_edge_valocity(
    const RiemannProblemSolution* solution) {
  double c2 = speed_of_sound(&(solution->right));
  double u2 = solution->right.velocity;
  double U = contact_surface_velocity(solution);
  double c2_star = c2 - G3 * (u2 - U);
  return U + c2_star;
}

double left_contact_density(const RiemannProblemSolution* solution) {
  double a1 = mass_velocity(solution->contact_pressure, &(solution->left));
  double r1 = solution->left.density;
  double u1 = solution->left.velocity;
  double p1 = solution->left.pressure;
  double U = contact_surface_velocity(solution);
  double P = solution->contact_pressure;
  if (P >= p1) {
    return r1 * a1 / (a1 - r1 * (u1 - U));
  } else {
    double c1 = speed_of_sound(&(solution->left));
    double c1_star = c1 + G3 * (u1 - U);
    return GAMMA * P / (gsl_pow_2(c1_star));
  }
}

double right_contact_density(const RiemannProblemSolution* solution) {
  double a2 = mass_velocity(solution->contact_pressure, &(solution->right));
  double r2 = solution->right.density;
  double u2 = solution->right.velocity;
  double p2 = solution->right.pressure;
  double U = contact_surface_velocity(solution);
  double P = solution->contact_pressure;
  if (P >= p2) {
    return r2 * a2 / (a2 - r2 * (u2 - U));
  } else {
    double c2 = speed_of_sound(&(solution->right));
    double c2_star = c2 + G3 * (u2 - U);
    return GAMMA * P / (gsl_pow_2(c2_star));
  }
}

GasParameters rarefaction_wave_solution(const RiemannProblemSolution* solution,
                                        double curve, int direction) {
  GasParameters flow;
  if (direction == LEFT_DIRECTION) {
    flow = solution->left;
  } else {
    flow = solution->right;
  }

  double P = solution->contact_pressure;
  double c = speed_of_sound(&flow);
  double p = flow.pressure;
  double u = flow.velocity;
  double D = u + direction * c;

  if (direction * curve >= direction * D) {
    return flow;
  }

  double U = contact_surface_velocity(solution);
  double c_star = c + direction * G3 * (U - u);
  double D_star = U + direction * c_star;
  double R = GAMMA * P / gsl_pow_2(c_star);

  if (direction * curve <= direction * D_star) {
    return (GasParameters){P, R, U};
  }

  c_star = 1 / G4 * c + direction * G5 * (curve - u);
  U = c_star - direction * curve;
  P = p * pow(c_star / c, 1 / G1);
  R = GAMMA * P / gsl_pow_2(c_star);
  return (GasParameters){P, R, U};
}

GasParameters shock_wave_solution(const RiemannProblemSolution* solution,
                                  double curve, int direction) {
  GasParameters flow;
  if (direction == LEFT_DIRECTION) {
    flow = solution->left;
  } else {
    flow = solution->right;
  }

  double P = solution->contact_pressure;
  double a = mass_velocity(P, &flow);
  double r = flow.density;
  double u = flow.velocity;
  double D = u + direction * a / r;

  if (direction * curve >= direction * D) {
    return flow;
  }

  double U = contact_surface_velocity(solution);
  double R = r * a / (a + direction * (u - U));
  return (GasParameters){P, R, U};
}

GasParameters solution_on_curve(const RiemannProblemSolution* solution,
                                double curve) {
  double U = contact_surface_velocity(solution);
  double P = solution->contact_pressure;
  int direction = curve >= U ? RIGHT_DIRECTION : LEFT_DIRECTION;
  GasParameters flow = curve >= U ? solution->right : solution->left;
  if (P >= flow.pressure) {
    return shock_wave_solution(solution, curve, direction);
  } else {
    return rarefaction_wave_solution(solution, curve, direction);
  }
}
