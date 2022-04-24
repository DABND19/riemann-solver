#include "riemann_problem.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>
#include <stdlib.h>

static const double EPSILON = 1e-6;

typedef struct __riemann_solver {
  gsl_root_fdfsolver* newton_solver;
  gsl_root_fsolver* bisection_solver;
  void* state;
} RiemannSolver;

struct __riemann_solver_state {
  GasFlow left;
  GasFlow right;
  double p_contact;
};

double Q_k(double p_c, const GasFlow* flow) {
  double p_k = flow->pressure;
  double r_k = flow->density;
  double A_k = 2. / ((GAMMA + 1.) * r_k);
  double B_k = (GAMMA - 1.) / (GAMMA + 1.) * p_k;
  return sqrt((p_c + B_k) / A_k);
}

double f_k(double p_c, const GasFlow* flow) {
  double p = p_c;
  double p_k = flow->pressure;
  if (p > p_k) {
    return (p - p_k) / Q_k(p, flow);
  } else {
    double a_k = sound_speed(flow);
    double tmp = pow(p / p_k, 0.5 * (GAMMA - 1.) / GAMMA);
    return 2. * a_k / (GAMMA - 1.) * (tmp - 1.);
  }
}

double U_c(double p_c, const GasFlow* left, const GasFlow* right) {
  double u_l = left->velocity;
  double u_r = right->velocity;
  double f_l = f_k(p_c, left);
  double f_r = f_k(p_c, right);
  return 0.5 * (u_l + u_r) + 0.5 * (f_r - f_l);
}

double df_k(double p_c, const GasFlow* flow) {
  double p = p_c;
  double p_k = flow->pressure;
  double r_k = flow->density;
  if (p > p_k) {
    double B_k = (GAMMA - 1.) / (GAMMA + 1.) * p_k;
    return (1. - 0.5 * (p - p_k) / (B_k + p)) / Q_k(p, flow);
  } else {
    double a_k = sound_speed(flow);
    return 1. / (r_k * a_k) * pow(p / p_k, - 0.5 * (GAMMA + 1.) / GAMMA);
  }
}

void fdf_k(double p_c, const GasFlow* flow, double* f_, double* df_) {
  *f_ = f_k(p_c, flow);
  *df_ = df_k(p_c, flow);
}

double f(double p_c, const GasFlow* left, const GasFlow* right) {
  double du = right->velocity - left->velocity;
  double f_l = f_k(p_c, left);
  double f_r = f_k(p_c, right);
  return f_l + f_r + du;
}

double df(double p_c, const GasFlow* left, const GasFlow* right) {
  double df_l = df_k(p_c, left);
  double df_r = df_k(p_c, right);
  return df_l + df_r;
}

void fdf(double p_c, const GasFlow* left, const GasFlow* right, double* f_, double* df_) {
  double du = right->velocity - left->velocity;
  double f_l, f_r, df_l, df_r;
  fdf_k(p_c, left, &f_l, &df_l);
  fdf_k(p_c, right, &f_r, &df_r);
  *f_ = f_l + f_r + du;
  *df_ = df_l + df_r;
}

struct gsl_f_params {
  const GasFlow* left;
  const GasFlow* right;
};

double gsl_f(double x, void* parameters) {
  struct gsl_f_params* params = (struct gsl_f_params*)parameters;
  return f(x, params->left, params->right);
}

double gsl_df(double x, void* parameters) {
  struct gsl_f_params* params = (struct gsl_f_params*)parameters;
  return df(x, params->left, params->right);
}

void gsl_fdf(double x, void* parameters, double* y, double* dy) {
  struct gsl_f_params* params = (struct gsl_f_params*)parameters;
  fdf(x, params->left, params->right, y, dy);
}

double P_TR(const GasFlow* left, const GasFlow* right) {
  double a_l = sound_speed(left);
  double a_r = sound_speed(right);
  double u_l = left->velocity;
  double u_r = right->velocity;
  double p_l = left->pressure;
  double p_r = right->pressure;
  double G = 0.5 * (GAMMA - 1.) / GAMMA;
  double nominator = (a_l + a_r - 0.5 * (GAMMA - 1.) * (u_r - u_l));
  double denominator = (a_l / pow(p_l, G) + a_r / pow(p_r, G));
  return pow(nominator / denominator, 1. / G);
}

const size_t MAX_ITER = 1000;

double P_c(const GasFlow* left, const GasFlow* right,
           gsl_root_fdfsolver* newton, gsl_root_fsolver* bisection,
           int* error) {
  double p_min = gsl_min(left->pressure, right->pressure);
  double p_max = gsl_max(left->pressure, right->pressure);
  double f_min = f(p_min, left, right);
  double f_max = f(p_max, left, right);
  double f_vacuum = f(0., left, right);

  *error = GSL_SUCCESS;

  if (gsl_fcmp(f_vacuum, 0., EPSILON) >= 0) {
    return 0.;
  }

  if (!gsl_fcmp(f_max, 0., EPSILON)) {
    return p_max;
  }

  if (!gsl_fcmp(f_min, 0., EPSILON)) {
    return p_min;
  }

  if (gsl_fcmp(f_min, 0., EPSILON) > 0 && gsl_fcmp(f_max, 0., EPSILON) > 0) {
    return P_TR(left, right);
  }

  double p_0 = f_max >= 0 ? p_min : p_max;
  double p = p_0;
  struct gsl_f_params equation_params = {left, right};
  gsl_function_fdf fdf_equation;
  fdf_equation.f = &gsl_f;
  fdf_equation.df = &gsl_df;
  fdf_equation.fdf = &gsl_fdf;
  fdf_equation.params = &equation_params;
  gsl_root_fdfsolver_set(newton, &fdf_equation, p_0);

  *error = GSL_CONTINUE;
  for (size_t i = 0; i < MAX_ITER && *error == GSL_CONTINUE; ++i) {
    *error = gsl_root_fdfsolver_iterate(newton);
    if (*error != GSL_SUCCESS) {
      break;
    }

    p_0 = p;
    p = gsl_root_fdfsolver_root(newton);

    *error = gsl_root_test_delta(p, p_0, 0., EPSILON);
  }

  if (*error && f_max < 0) {
    gsl_function f_equation;
    f_equation.function = &gsl_f;
    f_equation.params = &equation_params;

    double p_l = p_min;
    double p_r = p_max;

    gsl_root_fsolver_set(bisection, &f_equation, p_l, p_r);

    *error = GSL_CONTINUE;
    for (size_t i = 0; i < MAX_ITER && *error == GSL_CONTINUE; ++i) {
      *error = gsl_root_fsolver_iterate(bisection);
      if (*error != GSL_SUCCESS) {
        break;
      }

      p = gsl_root_fsolver_root(bisection);
      p_l = gsl_root_fsolver_x_lower(bisection);
      p_r = gsl_root_fsolver_x_upper(bisection);

      *error = gsl_root_test_interval(p_l, p_r, 0., EPSILON);
    }
  }

  return p;
}

RiemannSolver* riemann_solver_alloc() {
  RiemannSolver* solver = (RiemannSolver*)malloc(sizeof(RiemannSolver));
  solver->newton_solver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
  solver->bisection_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  solver->state = malloc(sizeof(struct __riemann_solver_state));
  return solver;
}

void riemann_solver_free(RiemannSolver* solver) {
  gsl_root_fdfsolver_free(solver->newton_solver);
  gsl_root_fsolver_free(solver->bisection_solver);
  free(solver->state);
  free(solver);
}

int riemann_solver_set(RiemannSolver* solver, const GasFlow* left, const GasFlow* right) {
  int error = GSL_SUCCESS;
  struct __riemann_solver_state* state = (struct __riemann_solver_state*)solver->state;

  state->left = *left;
  state->right = *right;

  double p_c = P_c(left, right, solver->newton_solver, solver->bisection_solver, &error);
  state->p_contact = p_c;

  return error;
}

double S_shock_left(double p_c, const GasFlow* left) {
  double r_l = left->density;
  double u_l = left->velocity;
  double Q_l = Q_k(p_c, left);
  return u_l - Q_l / r_l;
}

GasFlow left_shock_wave(const struct __riemann_solver_state* state, double s) {
  double r_l = state->left.density;
  double p_l = state->left.pressure;
  double p_c = state->p_contact;

  double S_l = S_shock_left(p_c, &(state->left));

  if (s <= S_l) {
    return state->left;
  } else {
    double p_ = p_c / p_l;
    double G = (GAMMA - 1.) / (GAMMA + 1.);
    double r_c = r_l * (p_ + G) / (G * p_ + 1.);

    double u_c = U_c(p_c, &(state->left), &(state->right));

    return (GasFlow){p_c, r_c, u_c};
  }
}

double S_rarefraction_left(double p_c, const GasFlow* left) {
  double a_l = sound_speed(left);
  double u_l = left->velocity;
  return u_l - a_l;
}

GasFlow left_rarefraction_wave(const struct __riemann_solver_state* state, double s) {
  double p_l = state->left.pressure;
  double r_l = state->left.density;
  double u_l = state->left.velocity;
  double a_l = sound_speed(&(state->left));
  double p_c = state->p_contact;
  double u_c = U_c(p_c, &(state->left), &(state->right));

  double p_ = p_c / p_l;

  double a_c = a_l * pow(p_, 0.5 * (GAMMA - 1.) / GAMMA);
  double S_l = S_rarefraction_left(p_c, &(state->left));
  double S_c = u_c - a_c;

  if (s <= S_l) {
    return state->left;
  }

  double r_c = r_l * pow(p_, 1. / GAMMA);
  if (s >= S_c) {
    return (GasFlow){p_c, r_c, u_c};
  }

  double G1 = 2. / (GAMMA + 1.);
  double G2 = (GAMMA - 1.) / (GAMMA + 1.);
  double r = r_l * pow(G1 + G2 / a_l * (u_l - s), 2. / (GAMMA - 1.));
  double u = G1 * (a_l + 0.5 * (GAMMA - 1.) * u_l + s);
  double p = p_l * pow(G1 + G2 / a_l * (u_l - s), 2. * GAMMA / (GAMMA - 1.));
  return (GasFlow){p, r, u};
}

double S_shock_right(double p_c, const GasFlow* right) {
  double r_r = right->density;
  double u_r = right->velocity;
  double Q_r = Q_k(p_c, right);
  return u_r + Q_r / r_r;
}

GasFlow right_shock_wave(const struct __riemann_solver_state* state, double s) {
  double r_r = state->right.density;
  double p_r = state->right.pressure;
  double p_c = state->p_contact;

  double S_r = S_shock_right(p_c, &(state->right));

  if (s >= S_r) {
    return state->right;
  } else {
    double p_ = p_c / p_r;
    double G = (GAMMA - 1.) / (GAMMA + 1.);
    double r_c = r_r * (p_ + G) / (G * p_ + 1.);

    double u_c = U_c(p_c, &(state->left), &(state->right));

    return (GasFlow){p_c, r_c, u_c};
  }
}

double S_rarefraction_right(double p_c, const GasFlow* right) {
  double a_r = sound_speed(right);
  double u_r = right->velocity;
  return u_r + a_r;
}

GasFlow right_rarefraction_wave(const struct __riemann_solver_state* state, double s) {
  double p_r = state->right.pressure;
  double r_r = state->right.density;
  double u_r = state->right.velocity;
  double a_r = sound_speed(&(state->right));
  double p_c = state->p_contact;
  double u_c = U_c(p_c, &(state->left), &(state->right));

  double p_ = p_c / p_r;

  double a_c = a_r * pow(p_, 0.5 * (GAMMA - 1.) / GAMMA);
  double S_r = S_rarefraction_right(p_c, &(state->right));
  double S_c = u_c + a_c;

  if (s >= S_r) {
    return state->right;
  }

  double r_c = r_r * pow(p_, 1. / GAMMA);
  if (s <= S_c) {
    return (GasFlow){p_c, r_c, u_c};
  }

  double G1 = 2. / (GAMMA + 1.);
  double G2 = (GAMMA - 1.) / (GAMMA + 1.);
  double r = r_r * pow(G1 - G2 / a_r * (u_r - s), 2. / (GAMMA - 1.));
  double u = G1 * (-a_r + 0.5 * (GAMMA - 1.) * u_r + s);
  double p = p_r * pow(G1 - G2 / a_r * (u_r - s), 2. * GAMMA / (GAMMA - 1.));
  return (GasFlow){p, r, u};
} 

double riemann_problem_contact_surface_velocity(const RiemannSolver* solver) {
  struct __riemann_solver_state* state = (struct __riemann_solver_state*)solver->state;
  return U_c(state->p_contact, &(state->left), &(state->right));
}

double riemann_problem_left_wave_velocity(const RiemannSolver* solver) {
  struct __riemann_solver_state* state = (struct __riemann_solver_state*)solver->state;
  if (state->p_contact > state->left.pressure) {
    return S_shock_left(state->p_contact, &(state->left));
  } else {
    return S_rarefraction_left(state->p_contact, &(state->left));
  }
}

double riemann_problem_right_wave_velocity(const RiemannSolver* solver) {
  struct __riemann_solver_state* state = (struct __riemann_solver_state*)solver->state;
  if (state->p_contact > state->right.pressure) {
    return S_shock_right(state->p_contact, &(state->right));
  } else {
    return S_rarefraction_right(state->p_contact, &(state->right));
  }
}

GasFlow riemann_problem_solution(const RiemannSolver* solver, double s) {
  struct __riemann_solver_state* state = (struct __riemann_solver_state*)solver->state;

  double p_c = state->p_contact;
  double u_c = U_c(p_c, &(state->left), &(state->right));

  if (gsl_fcmp(s, u_c, EPSILON) < 0) {
    if (p_c > state->left.pressure) {
      return left_shock_wave(state, s);
    } else {
      return left_rarefraction_wave(state, s);
    }
  } else {
    if (p_c > state->right.pressure) {
      return right_shock_wave(state, s);
    } else {
      return right_rarefraction_wave(state, s);
    }
  }
}
