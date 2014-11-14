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

#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "structs.h"
#include "basis.h"
#include "tripalg.h"
#include "geometry.h"
#include "poisson.h"
#include "stokes.h"
#include "global.h"
#include "solver.h"
#include "precon.h"
#include "test.h"

#define A 0.78539816
#define D 1.57079633

static double ethier_steinman_u(double x, double y, double z, double t);
static double ethier_steinman_v(double x, double y, double z, double t);
static double ethier_steinman_w(double x, double y, double z, double t);
static double ethier_steinman_p(double x, double y, double z, double t);
static double ethier_steinman_dudt(double x, double y, double z, double t);
static double ethier_steinman_dvdt(double x, double y, double z, double t);
static double ethier_steinman_dwdt(double x, double y, double z, double t);

double zero_func(double x, double y, double z, double t)
{
	return 0.0;
}

double one_func(double x, double y, double z, double t)
{
	return 1.0;
}

static double ethier_steinman_u(double x, double y, double z, double t)
{
	return -A * (exp(A * x) * sin(A * y + D * z) + exp(A * z)
			* cos(A * x + D * y)) * exp(-D * D * t);
}

static double ethier_steinman_v(double x, double y, double z, double t)
{
	return -A * (exp(A * y) * sin(A * z + D * x) + exp(A * x)
			* cos(A * y + D * z)) * exp(-D * D * t);
}

static double ethier_steinman_w(double x, double y, double z, double t)
{
	return -A * (exp(A * z) * sin(A * x + D * y) + exp(A * y)
			* cos(A * z + D * x)) * exp(-D * D * t);
}

static double ethier_steinman_p(double x, double y, double z, double t)
{
	const double tmp1 = exp(2.0 * A * x) + exp(2.0 * A * y)
			+ exp(2.0 * A * z);
	const double tmp2 = 2.0 * sin(A * x + D * y) * cos(A * z + D * x)
			* exp(A * (y + z));
	const double tmp3 = 2.0 * sin(A * y + D * z) * cos(A * x + D * y)
			* exp(A * (z + x));
	const double tmp4 = 2.0 * sin(A * z + D * x) * cos(A * y + D * z)
			* exp(A * (x + y));
	return -A * A * (tmp1 + tmp2 + tmp3 + tmp4) * exp(-2.0 * D * D * t)
			/ 2.0;
}

static double ethier_steinman_dudt(double x, double y, double z, double t)
{
	return -D * D * ethier_steinman_u(x, y, z, t);
}

static double ethier_steinman_dvdt(double x, double y, double z, double t)
{
	return -D * D * ethier_steinman_v(x, y, z, t);
}

static double ethier_steinman_dwdt(double x, double y, double z, double t)
{
	return -D * D * ethier_steinman_w(x, y, z, t);
}

double max_error(struct element **sel, int var, double (*exact_solution) ())
{
	double ***u;
	double err = 0.0;
	double tmp;

	for (int i = 0; i < (*sel)->nel; ++i) {
		const double * const x = **sel[i]->geom->grid->x;
		const double * const y = **sel[i]->geom->grid->y;
		const double * const z = **sel[i]->geom->grid->z;
		u = variable_array(sel[i], var);
		for (int j = 0; j < (*sel)->ntot; ++j) {
			tmp = fabs(*(**u + j) - exact_solution(*(x + j),
					*(y + j), *(z + j), (*sel)->params->t));
			if (tmp > err)
				err = tmp;
		}
	}
	err = global_max(sel, err);

	return err;
}

double maxmin_val(struct element **sel, int var, double maxmin)
{
	double err = 0.0;
	double tmp;

	for (int i = 0; i < (*sel)->nel; ++i) {
		double ***u = variable_array(sel[i], var);
		for (int j = 0; j < (*sel)->ntot; ++j) {
			tmp = *(**u + j);
			if (maxmin * tmp > err)
				err = maxmin * tmp;
		}
	}
	err = global_max(sel, err);

	return maxmin * err;
}

double max_divergence(struct element **sel)
{
	double maxdiv = 0.0;
	const int n = (*sel)->basis->n;
	double ****div = new_4d_array(n, n, n, (*sel)->nel);
	double tmp;

	for (int i = 0; i < (*sel)->nel; ++i)
		discrete_divergence_operator(sel[i], sel[i]->u, div[i]);
	direct_stiffness_summation(sel, div);
	global_inverse_mass_matrix_operator(sel, div, div);

	for (int i = 0; i < (*sel)->nel; ++i) {
		for (int j = 0; j < (*sel)->ntot; ++j) {
			tmp = fabs(*(**div[i] + j));
			if (tmp > maxdiv)
				maxdiv = tmp;
		}
	}
	maxdiv = global_max(sel, maxdiv);
	free_4d_array(div);

	return maxdiv;
}

void set_exact_pressure(struct element **sel, double (*exact_solution) (),
		double ****p)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		const double * const x = **sel[i]->geom->grid->x;
		const double * const y = **sel[i]->geom->grid->y;
		const double * const z = **sel[i]->geom->grid->z;
		for (int j = 0; j < (*sel)->ntot; ++j)
			*(**p[i] + j) = exact_solution(*(x + j), *(y + j),
					*(z + j), (*sel)->params->t);
	}
}

double pressure_error(struct element **sel, double (*exact_solution) ())
{
	const int n = (*sel)->basis->n;
	double ****p = new_4d_array(n, n, n, (*sel)->nel);

	set_exact_pressure(sel, exact_solution, p);
	subtract_nullspace(sel, p, PRESSURE);

	double tmp;
	double max;
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < (*sel)->ntot; ++j) {
			tmp = fabs(*(**p[i] + j) - *(**sel[i]->p + j));
			if (tmp > max)
				max = tmp;
		}
	max = global_max(sel, max);
	free_4d_array(p);
	
	return max;
}

void initialise_test(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		sel[i]->bc_func[0] = ethier_steinman_u;
		sel[i]->bc_func[1] = ethier_steinman_v;
		sel[i]->bc_func[2] = ethier_steinman_w;
		sel[i]->bc_dfunc[0] = ethier_steinman_dudt;
		sel[i]->bc_dfunc[1] = ethier_steinman_dvdt;
		sel[i]->bc_dfunc[2] = ethier_steinman_dwdt;
	}
	initialise_from_function(sel, ethier_steinman_u, ethier_steinman_v,
			ethier_steinman_w);
}
