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
#include "structs.h"
#include "basis.h"
#include "linalg.h"
#include "tripalg.h"
#include "convection.h"
#include "poisson.h"
#include "global.h"
#include "time.h"

void time_shift(struct trivec *u, struct trivec *old, struct trivec *old_old,
		int n)
{
	for (int i = 0; i < 3; ++i)
		single_time_shift(*u->ind[i], *old->ind[i], *old_old->ind[i], n);
}

void single_time_shift(double ***data, double ***old, double ***old_old, int n)
{
	copy_array(**old, **old_old, n);
	copy_array(**data, **old, n);
}

void stiffly_stable(struct trivec *u1, struct trivec *u2, struct trivec *u3,
		struct trivec *u, size_t ntot, double c, int imp_exp)
{
	for (int i = 0; i < 3; ++i)
		single_stiffly_stable(*u1->ind[i], *u2->ind[i], *u3->ind[i],
				*u->ind[i], ntot, c, imp_exp);
}

void single_stiffly_stable(double ***u1, double ***u2, double ***u3,
		double ***u, size_t ntot, double c, int imp_exp)
{
	double a1, a2, a3;

	if (imp_exp == IMPLICIT) {
		a1 = 3.0 * c;
		a2 = -(3.0 / 2.0) * c;
		a3 = (1.0 / 3.0) * c;
	} else if (imp_exp == EXPLICIT2) {
		a1 = 2.0 * c;
		a2 = -1.0 * c;
		a3 = 0.0 * c;
	} else if (imp_exp == EXPLICIT3) {
		a1 = 3.0 * c;
		a2 = -3.0 * c;
		a3 = 1.0 * c;
	}

	add_array(**u1, **u, ntot, a1);
	add_array(**u2, **u, ntot, a2);
	add_array(**u3, **u, ntot, a3);
}

void integrate_convection_se3(struct element **sel, struct trivec **u)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		stiffly_stable(sel[i]->conv, sel[i]->old->conv, sel[i]->old->old
				->conv, u[i], sel[i]->ntot, sel[i]->params->dt,
				EXPLICIT3);
	trivec_inverse_mass_matrix(sel, u);
}

void integrate_convection_rk4(struct element **sel, struct trivec **u,
		double t0, double t, double dt, double beta)
{
	const int n = (*sel)->basis->n;
	const int nel = (*sel)->nel;
	const double a[4] = { 0.0, 0.5, 0.5, 1.0 };
	const double b[4] = { 1.0, 2.0, 2.0, 1.0 };
	void *(*newtri) () = (void *(*) ()) new_trivec;

	struct trivec **f[4] = { new_ptr_array(nel, n, newtri),
			new_ptr_array(nel, n, newtri), new_ptr_array(nel, n,
			newtri),  new_ptr_array(nel, n, newtri) };
	struct trivec **y[4] = { new_ptr_array(nel, n, newtri),
			new_ptr_array(nel, n, newtri), new_ptr_array(nel, n,
			newtri),  new_ptr_array(nel, n, newtri) };
	struct trivec **utmp = new_ptr_array(nel, n, newtri);

	add_interpolated_velocity(sel, u, t0, beta);
	int k = 0;
	while(t0 + (double) k++ * dt < t - 0.5 * dt) {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < (*sel)->nel; ++j) {
				copy_trivec(u[j], y[i][j], n * n * n);
				if (i > 0)
					add_trivec(f[i - 1][j], y[i][j], n * n
							* n, -dt * a[i]);
				scale_trivec(utmp[j], n * n * n, 0.0);
			}
			add_interpolated_velocity(sel, utmp, t0
					+ ((double) (k - 1) + a[i]) * dt, 1.0);
			for (int j = 0; j < (*sel)->nel; ++j)
				convection_operator(sel[j], utmp[j], y[i][j],
						f[i][j]);
			trivec_inverse_mass_matrix(sel, f[i]);
		}
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < (*sel)->nel; ++j)
				add_trivec(f[i][j], u[j], n * n * n,
						-dt * b[i] / 6.0);
	}

	free_ptr_array(utmp, nel, free_trivec);
	for (int i = 0; i < 4; ++i) {
		free_ptr_array(f[i], nel, free_trivec);
		free_ptr_array(y[i], nel, free_trivec);
	}
}

void add_interpolated_velocity(struct element **sel, struct trivec **u,
		double t, double alpha)
{
	const double t0 = (*sel)->params->t;
	const double dt0 = (*sel)->params->dt;
	const double dt[3] = { t - t0 + 2.0 * dt0, t - t0 + dt0, t - t0 };
	const double a = 0.5 * dt[1] * dt[2] / pow(dt0, 2);
	const double b = -dt[0] * dt[2] / pow(dt0, 2);
	const double c = 0.5 * dt[0] * dt[1] / pow(dt0, 2);
	const double eps = 1.0e-6;
	const int ntot = (*sel)->basis->ntot;

	for (int i = 0; i < (*sel)->nel; ++i) {
		if (fabs(a) > eps)
			add_trivec(sel[i]->old->old->u, u[i], ntot, a * alpha);
		if (fabs(b) > eps)
			add_trivec(sel[i]->old->u, u[i], ntot, b * alpha);
		if (fabs(c) > eps)
			add_trivec(sel[i]->u, u[i], ntot, c * alpha);
	}
}

void trivec_inverse_mass_matrix(struct element **sel, struct trivec **u)
{
	const int n = (*sel)->basis->n;
	double ****tmp = new_4d_array(n, n, n, (*sel)->nel);
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < (*sel)->nel; ++i)
			copy_array(***u[i]->ind[j], **tmp[i], n * n * n);
		direct_stiffness_summation(sel, tmp);
		global_inverse_mass_matrix_operator(sel, tmp, tmp);
		for (int i = 0; i < (*sel)->nel; ++i)
			copy_array(**tmp[i], ***u[i]->ind[j], n * n * n);
	}
	free_4d_array(tmp);
}
