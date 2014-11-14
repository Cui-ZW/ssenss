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

#include "utils.h"
#include "linalg.h"
#include "tripalg.h"
#include "basis.h"
#include "deriv.h"
#include "time.h"
#include "convection.h"
#include "stokes.h"
#include "test.h"
#include "dataio.h"
#include "geometry.h"
#include "precon.h"
#include "global.h"
#include "poisson.h"

void discrete_laplace_operator(struct element *sel, double ***u,
		double ***lap_u)
{
	struct trivec *tmp = new_trivec(sel->basis->n);
	struct trivec *dud = new_trivec(sel->basis->n);

	discrete_x_differentiation(sel->basis, u, dud->x, ADD_NEW,
			REPLACE_OLD);
	discrete_y_differentiation(sel->basis, u, dud->y, ADD_NEW,
			REPLACE_OLD);
	discrete_z_differentiation(sel->basis, u, dud->z, ADD_NEW,
			REPLACE_OLD);

	trimat_trivec_product(sel->geom->geom_mat, dud, sel->basis->n, tmp);

	for (int i = 0; i < 3; ++i)
		test_function_inner_product(sel, *tmp->ind[i], *tmp->ind[i]);

	transpose_divergence(sel->basis, tmp, lap_u);

	free_trivec(tmp);
	free_trivec(dud);
}

void test_function_inner_product(struct element *sel, double ***u, double ***vu)
{
	double *tmp = new_1d_array(sel->basis->ntot);
	multiply(**u, **sel->basis->weights_array, tmp, sel->basis->ntot,
			ADD_NEW, REPLACE_OLD);
	multiply(tmp, **sel->geom->jac_det, **vu, sel->basis->ntot, ADD_NEW,
			REPLACE_OLD);
	free_1d_array(tmp);
}

void inverse_mass_matrix_operator(struct element *sel, double ***u,
		double ***vu)
{
	double *tmp = new_1d_array(sel->basis->ntot);
	divide(**u, **sel->basis->weights_array, tmp, sel->basis->ntot,
			ADD_NEW, REPLACE_OLD);
	divide(tmp, **sel->geom->jac_det, **vu, sel->basis->ntot, ADD_NEW,
			REPLACE_OLD);
	free_1d_array(tmp);
}

void global_inverse_mass_matrix_operator(struct element **sel, double ****u,
		double ****mu)
{
	double *tmp = new_1d_array((*sel)->basis->ntot * (*sel)->nel);
	multiply(***u, **(*sel)->precon->inverse_mass, tmp, (*sel)->basis->ntot
			* (*sel)->nel, 1.0, 0.0);
	copy_array(tmp, ***mu, (*sel)->basis->ntot * (*sel)->nel);
	free_1d_array(tmp);
}

void mask_dirichlet_boundaries(struct element **sel, double ****u, int var)
{
	if (var == PRESSURE)
		return;
	const int n = (*sel)->basis->n;
	double ****tmp = new_4d_array(n, n, n, (*sel)->nel);
	multiply(***u, **(*sel)->precon->dirichlet_mask[(*sel)->precon->level],
			***tmp, n * n * n * (*sel)->nel, 1.0, 0.0);
	copy_array(***tmp, ***u, n * n * n * (*sel)->nel);
	free_4d_array(tmp);
}

void set_dirichlet_boundaries(struct element *sel, double ***u, int bc,
		double (*bc_func) ())
{
	const int n = sel->basis->n;
	double ***x = sel->geom->grid->x;
	double ***y = sel->geom->grid->y;
	double ***z = sel->geom->grid->z;
	double t = sel->params->t + sel->params->dt;

	if (bc & DIRICHLET_W)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) 
				u[j][i][0] = bc_func(x[j][i][0], y[j][i][0],
						z[j][i][0], t);
	if (bc & DIRICHLET_E)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) 
				u[j][i][n - 1] = bc_func(x[j][i][n - 1],
						y[j][i][n - 1], z[j][i][n - 1],
						t);
	if (bc & DIRICHLET_S)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) 
				u[j][0][i] = bc_func(x[j][0][i], y[j][0][i],
						z[j][0][i], t);
	if (bc & DIRICHLET_N)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) 
				u[j][n - 1][i] = bc_func(x[j][n - 1][i],
						y[j][n - 1][i], z[j][n - 1][i],
						t);
	if (bc & DIRICHLET_B)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) 
				u[0][j][i] = bc_func(x[0][j][i], y[0][j][i],
						z[0][j][i], t);
	if (bc & DIRICHLET_T)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) 
				u[n - 1][j][i] = bc_func(x[n - 1][j][i],
						y[n - 1][j][i], z[n - 1][j][i],
						t);
}

void add_surface_integral(struct element *sel, double ***g, double ***f,
		int deriv, int bc)
{
	const int m = sel->n - 1;
	struct surface *norm;

	if (deriv == X_DERIV)
		norm = sel->geom->surf_norm_x;
	else if (deriv == Y_DERIV)
		norm = sel->geom->surf_norm_y;
	else if (deriv == Z_DERIV)
		norm = sel->geom->surf_norm_z;

	if (!(bc & DIRICHLET_W))
		for (int j = 0; j < sel->n; ++j)
			for (int i = 0; i < sel->n; ++i)
				f[j][i][0] += sel->basis->weights[i] *
						sel->basis->weights[j] *
						norm->w[j][i] * g[j][i][0];
	if (!(bc & DIRICHLET_E))
		for (int j = 0; j < sel->n; ++j)
			for (int i = 0; i < sel->n; ++i)
				f[j][i][m] += sel->basis->weights[i] *
						sel->basis->weights[j] *
						norm->e[j][i] * g[j][i][m];
	if (!(bc & DIRICHLET_S))
		for (int j = 0; j < sel->n; ++j)
			for (int i = 0; i < sel->n; ++i)
				f[j][0][i] += sel->basis->weights[i] *
						sel->basis->weights[j] *
						norm->s[j][i] * g[j][0][i];
	if (!(bc & DIRICHLET_N))
		for (int j = 0; j < sel->n; ++j)
			for (int i = 0; i < sel->n; ++i)
				f[j][m][i] += sel->basis->weights[i] *
						sel->basis->weights[j] *
						norm->n[j][i] * g[j][m][i];
	if (!(bc & DIRICHLET_B))
		for (int j = 0; j < sel->n; ++j)
			for (int i = 0; i < sel->n; ++i)
				f[0][j][i] += sel->basis->weights[i] *
						sel->basis->weights[j] *
						norm->b[j][i] * g[0][j][i];
	if (!(bc & DIRICHLET_T))
		for (int j = 0; j < sel->n; ++j)
			for (int i = 0; i < sel->n; ++i)
				f[m][j][i] += sel->basis->weights[i] *
						sel->basis->weights[j] *
						norm->t[j][i] * g[m][j][i];
}

double ***variable_array(struct element *sel, int var)
{
	switch (var) {
		case PRESSURE:
			return sel->p;
		case X_VELOCITY:
			return sel->u->x;
		case Y_VELOCITY:
			return sel->u->y;
		case Z_VELOCITY:
			return sel->u->z;
		default:
			return NULL;
	}
}

double ***source_array(struct element *sel, int var)
{
	switch (var) {
		case X_SOURCE:
			return sel->f->x;
		case Y_SOURCE:
			return sel->f->y;
		case Z_SOURCE:
			return sel->f->z;
		default:
			return NULL;
	}
}

void set_lifted_velocity(struct element **sel, double t, double (*ud) (),
		double (*vd) (), double (*wd) ())
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		const double * const x = **sel[i]->geom->grid->x;
		const double * const y = **sel[i]->geom->grid->y;
		const double * const z = **sel[i]->geom->grid->z;
		for (int j = 0; j < (*sel)->ntot; ++j) {
			*(**sel[i]->ud->x + j) = ud(*(x + j), *(y + j),
					*(z + j), t);
			*(**sel[i]->ud->y + j) = vd(*(x + j), *(y + j),
					*(z + j), t);
			*(**sel[i]->ud->z + j) = wd(*(x + j), *(y + j),
					*(z + j), t);
		}
	}

}

void initialise_from_function(struct element **sel, double (*ud) (),
		double (*vd) (), double (*wd) ())
{
	const double dt = (*sel)->params->dt;
	double t = (*sel)->params->t - 2.0 * dt;

	for (int j = 0; j < 3; ++j) {
		set_lifted_velocity(sel, t, ud, vd, wd);
		for (int i = 0; i < (*sel)->nel; ++i) {
			time_shift(sel[i]->u, sel[i]->old->u,
					sel[i]->old->old->u, sel[i]->ntot);
			copy_trivec(sel[i]->ud, sel[i]->u, (*sel)->ntot);

			if (j >= 2)
				continue;

			time_shift(sel[i]->conv, sel[i]->old->conv,
					sel[i]->old->old->conv, sel[i]->ntot);
			time_shift(sel[i]->rot, sel[i]->old->rot,
					sel[i]->old->old->rot, sel[i]->ntot);

			convection_operator(sel[i], sel[i]->u, sel[i]->u,
					sel[i]->conv);
			double_curl_operator(sel[i], sel[i]->u, sel[i]->rot);
		}
		t += dt;
	}
	set_lifted_velocity(sel, t, zero_func, zero_func, zero_func);
}

void initialise_from_file(struct element **sel, char *filename)
{
	const double dt = (*sel)->params->dt;
	double t = (*sel)->params->t - 2.0 * dt;
	const int n = (*sel)->basis->n;
	double ****tmp = new_4d_array(n, n, n, (*sel)->nel);

	for (int j = 0; j < 3; ++j) {
		read_velocity(sel, filename);
		for (int k = 0; k < 3; ++k) {
			for (int i = 0; i < (*sel)->nel; ++i)
				copy_array(***sel[i]->u->ind[k], **tmp[i],
						n * n * n);
			direct_stiffness_averaging(sel, tmp);
			for (int i = 0; i < (*sel)->nel; ++i) {
				set_dirichlet_boundaries(sel[i], tmp[i], sel[i]
						->bc, sel[i]->bc_func[k]);
				copy_array(**tmp[i], ***sel[i]->u->ind[k],
						n * n * n);
			}
		}
		for (int i = 0; i < (*sel)->nel; ++i) {
			time_shift(sel[i]->u, sel[i]->old->u,
					sel[i]->old->old->u, sel[i]->ntot);
			if (j >= 2)
				continue;

			time_shift(sel[i]->conv, sel[i]->old->conv,
					sel[i]->old->old->conv, sel[i]->ntot);
			time_shift(sel[i]->rot, sel[i]->old->rot,
					sel[i]->old->old->rot, sel[i]->ntot);

			convection_operator(sel[i], sel[i]->u, sel[i]->u,
					sel[i]->conv);
			double_curl_operator(sel[i], sel[i]->u, sel[i]->rot);
		}
		t += (*sel)->params->dt;
	}
	free_4d_array(tmp);
}

void set_source_term(struct element **sel, int var, double (*func) (),
		double alpha)
{
	double ***f;

	for (int i = 0; i < (*sel)->nel; ++i) {
		const double * const x = **sel[i]->geom->grid->x;
		const double * const y = **sel[i]->geom->grid->y;
		const double * const z = **sel[i]->geom->grid->z;
		f = source_array(sel[i], var);
		for (int j = 0; j < (*sel)->ntot; ++j)
			*(**f + j) = func(*(x + j), *(y + j), *(z + j),
					sel[i]->params->t + sel[i]->params->dt);
		scale_array(**f, sel[i]->ntot, alpha);
	}
}
