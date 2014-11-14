/*
 * Copyright 2014 Christopher Nilsen
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <math.h>
#include "utils.h"
#include "basis.h"
#include "linalg.h"
#include "tripalg.h"
#include "time.h"
#include "poisson.h"
#include "solver.h"
#include "convection.h"
#include "stokes.h"
#include "helmholtz.h"
#include "global.h"
#include "navier_stokes.h"

#define ABSTOL 1.0e-10
#define RELTOL 1.0e-10

static void integrate_convection_operator(struct element **sel,
		struct trivec **u_tilde);
static double residual_tolerance(struct element **sel, double ****rhs,
		double abs_tol, double rel_tol);
static void estimate_pressure(struct element **sel);
static void compute_pressure_rhs(struct element **sel, struct trivec **u_tilde,
		double ****rhs);
static void solve_pressure_poisson(struct element **sel,
		struct trivec **u_tilde);
static void estimate_velocity(struct element **sel);
static void compute_velocity_rhs(struct element **sel, struct trivec **u_tilde,
		double ****rhs[3]);
static void solve_velocity_helmholtz(struct element **sel,
		struct trivec **u_tilde);

const static int deriv[3] = { X_DERIV, Y_DERIV, Z_DERIV };

static void integrate_convection_operator(struct element **sel,
		struct trivec **u_tilde)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		time_shift(sel[i]->conv, sel[i]->old->conv,
				sel[i]->old->old->conv, sel[i]->ntot);
		convection_operator(sel[i], sel[i]->u, sel[i]->u, sel[i]->conv);
	}
	integrate_convection_se3(sel, u_tilde);
	for (int i = 0; i < (*sel)->nel; ++i) {
		scale_trivec(u_tilde[i], sel[i]->ntot, -1.0);
		add_trivec(sel[i]->f, u_tilde[i], sel[i]->ntot,
				sel[i]->params->dt);
		stiffly_stable(sel[i]->u, sel[i]->old->u, sel[i]->old->old->u,
				u_tilde[i], sel[i]->ntot, ADD_NEW, IMPLICIT);
	}
}

static void lagrangian_convection_integration(struct element **sel,
		struct trivec **u_tilde)
{
	const double dt = (*sel)->params->dt;
	const double dtconv = dt / 2.0;
	double t[4] = { (*sel)->params->t - 2.0 * dt, (*sel)->params->t - dt,
			(*sel)->params->t, (*sel)->params->t + dt };

	for (int i = 0; i < (*sel)->nel; ++i) {
		time_shift(sel[i]->conv, sel[i]->old->conv,
				sel[i]->old->old->conv, sel[i]->ntot);
		convection_operator(sel[i], sel[i]->u, sel[i]->u, sel[i]->conv);
	}

	integrate_convection_rk4(sel, u_tilde, t[0], t[1], dtconv, 1.0 / 3.0);
	integrate_convection_rk4(sel, u_tilde, t[1], t[2], dtconv, -3.0 / 2.0);
	integrate_convection_rk4(sel, u_tilde, t[2], t[3], dtconv, 3.0);

	for (int i = 0; i < (*sel)->nel; ++i)
		add_trivec(sel[i]->f, u_tilde[i], sel[i]->ntot,
				sel[i]->params->dt);
}

static double residual_tolerance(struct element **sel, double ****rhs,
		double abs_tol, double rel_tol)
{
	const double tmp = sqrt(global_dot_product(sel, rhs, rhs));
	return dmax(1.0e-10, dmin(abs_tol, tmp * rel_tol));
}

static void estimate_pressure(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		double  ***tmp = new_3d_array(sel[i]->n, sel[i]->n, sel[i]->n);
		single_stiffly_stable(sel[i]->p, sel[i]->old->p,
				sel[i]->old->old->p, tmp, sel[i]->ntot, 1.0,
				EXPLICIT3);
		single_time_shift(sel[i]->p, sel[i]->old->p,
				sel[i]->old->old->p, sel[i]->ntot);
		copy_array(**tmp, **sel[i]->p, sel[i]->ntot);
		free_3d_array(tmp);
	}
}

static void compute_pressure_rhs(struct element **sel, struct trivec **u_tilde,
		double ****rhs)
{
	const int nel = (*sel)->nel;
	const int n = (*sel)->basis->n;
	double ****tmp[3] = { new_4d_array(n, n, n, nel), new_4d_array(n, n, n,
			nel), new_4d_array(n, n, n, nel) };

	for (int i = 0; i < (*sel)->nel; ++i) {
		struct trivec *dpd = new_trivec(n);
		estimate_pressure_gradient(sel[i], dpd);
		for (int j = 0; j < 3; ++j)
			copy_array(***dpd->ind[j], **tmp[j][i], n * n * n);
		free_trivec(dpd);
	}
	for (int j = 0; j < 3; ++j) {
		direct_stiffness_summation(sel, tmp[j]);
		global_inverse_mass_matrix_operator(sel, tmp[j], tmp[j]);
	}
	for (int i = 0; i < (*sel)->nel; ++i) {
		discrete_divergence_operator(sel[i], u_tilde[i], rhs[i]);
		scale_array(**rhs[i], sel[i]->ntot, -sel[i]->params->dens
				/ sel[i]->params->dt);

		for (int j = 0; j < 3; ++j)
			add_surface_integral(sel[i], tmp[j][i], rhs[i],
					deriv[j], ~sel[i]->bc);
	}
	direct_stiffness_summation(sel, rhs);
	for (int i = 0; i < 3; ++i)
		free_4d_array(tmp[i]);
}

static void solve_pressure_poisson(struct element **sel,
		struct trivec **u_tilde)
{
	const int n = (*sel)->basis->n;
	double ****rhs = new_4d_array(n, n, n, (*sel)->nel);

	compute_pressure_rhs(sel, u_tilde, rhs);
	const double eps = residual_tolerance(sel, rhs, ABSTOL, RELTOL);
	estimate_pressure(sel);
	solve_poisson_cg(sel, rhs, PRESSURE, eps);
	free_4d_array(rhs);
}

static void estimate_velocity(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		struct trivec *tmp = new_trivec(sel[i]->n);
		stiffly_stable(sel[i]->u, sel[i]->old->u, sel[i]->old->old->u,
				tmp, sel[i]->ntot, 1.0, EXPLICIT3);
		time_shift(sel[i]->u, sel[i]->old->u, sel[i]->old->old->u,
				sel[i]->ntot);
		copy_trivec(tmp, sel[i]->u, sel[i]->ntot);
		for (int j = 0; j < 3; ++j)
			set_dirichlet_boundaries(sel[i], *sel[i]->u->ind[j],
					sel[i]->bc, sel[i]->bc_func[j]);
		free_trivec(tmp);
	}
}

static void compute_velocity_rhs(struct element **sel, struct trivec **u_tilde,
		double ****rhs[3])
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		scale_trivec(u_tilde[i], sel[i]->ntot, 1.0 / (sel[i]->params->dt
					* sel[i]->params->visc));
		struct trivec *tmp = new_trivec(sel[i]->n);
		discrete_gradient_operator(sel[i], sel[i]->p, tmp);

		for (int j = 0; j < 3; ++j) {
			copy_array(***u_tilde[i]->ind[j], **rhs[j][i],
					(*sel)->basis->ntot);
			test_function_inner_product(sel[i], rhs[j][i],
					rhs[j][i]);
			add_surface_integral(sel[i], sel[i]->p, *tmp->ind[j],
					deriv[j], sel[i]->bc);
			add_array(***tmp->ind[j], **rhs[j][i], sel[i]->ntot,
					-1.0 / (sel[i]->params->visc
					* sel[i]->params->dens));
		}
		free_trivec(tmp);
	}
	for (int j = 0; j < 3; ++j)
		direct_stiffness_summation(sel, rhs[j]);
}

static void solve_velocity_helmholtz(struct element **sel,
		struct trivec **u_tilde)
{
	const int nel = (*sel)->nel;
	const int n = (*sel)->basis->n;
	double ****rhs[3] = { new_4d_array(n, n, n, nel), new_4d_array(n, n, n,
			nel), new_4d_array(n, n, n, nel) };
	compute_velocity_rhs(sel, u_tilde, rhs);

	const int velocity[3] = { X_VELOCITY, Y_VELOCITY, Z_VELOCITY };
	const double eps[3] = { residual_tolerance(sel, rhs[0], ABSTOL, RELTOL),
		residual_tolerance(sel, rhs[1], ABSTOL, RELTOL),
		residual_tolerance(sel, rhs[2], ABSTOL, RELTOL) };

	estimate_velocity(sel);
	for (int i = 0; i < 3; ++i) {
		solve_helmholtz_cg(sel, rhs[i], velocity[i], eps[i]);
		free_4d_array(rhs[i]);
	}
}

void integrate_navier_stokes(struct element **sel)
{
	struct trivec **u_tilde = new_ptr_array((*sel)->nel, (*sel)->n,
			(void *(*) ()) new_trivec);
	//integrate_convection_operator(sel, u_tilde);
	lagrangian_convection_integration(sel, u_tilde);
	solve_pressure_poisson(sel, u_tilde);
	solve_velocity_helmholtz(sel, u_tilde);

	free_ptr_array(u_tilde, (*sel)->nel, free_trivec);
	(*sel)->params->t += (*sel)->params->dt;
}
