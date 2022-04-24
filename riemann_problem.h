#pragma once

#include "gas_flow.h"

typedef struct __riemann_solver RiemannSolver;

RiemannSolver* riemann_solver_alloc();
void riemann_solver_free(RiemannSolver* solver);
int riemann_solver_set(RiemannSolver* solver, const GasFlow* left, const GasFlow* right);
double riemann_problem_contact_surface_velocity(const RiemannSolver* solution);
double riemann_problem_left_wave_velocity(const RiemannSolver* solution);
double riemann_problem_right_wave_velocity(const RiemannSolver* solution);
GasFlow riemann_problem_solution(const RiemannSolver* solver, double s);
