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
#include "basis.h"
#include "geometry.h"
#include "linalg.h"
#include "tripalg.h"
#include "deriv.h"
#include "poisson.h"
#include "helmholtz.h"
#include "global.h"
#include "interpolation.h"
#include "solver.h"
#include "init.h"
#include "test.h"
#include "precon.h"

static void initialise_multigrid_precon(struct element **sel);
static void compute_multigrid_sendcounts(struct element **sel);
static double ***precon_array(struct element *sel, int var, int level);
static double multigrid_recurse(struct element **sel, double ****f,
		double ****x, int var, int level, double eps,
		void (*opeval) ());
static void switch_multigrid_level(struct element **sel, int level);
static void multigrid_restriction(struct element **sel, double ****x,
		double ****y, int level);
static void multigrid_prolongation(struct element **sel, double ****x,
		double ****y, int level);
static void compute_residual(struct element **sel, double ****x, double ****f,
		double ****r, int var, void (*opeval) ());
static void diagonal_preconditioner(struct element **sel, double ****x,
		double ****px, int var, void (*opeval) (), double omega,
		double alpha);
static void reconnect_elements(struct element **sel);
static void initialise_diagonal_precon(struct element **sel, int var);
static void initialise_stiffness_summation_scaling(struct element **sel);
static void initialise_dirichlet_mask(struct element **sel);
static void initialise_inverse_mass_matrix(struct element **sel);

struct preconditioner *new_preconditioner(int n, int nlevels)
{
	struct preconditioner *precon = ec_malloc(sizeof(*precon));
	precon->level = 0;
	precon->nmg = new_int_array(nlevels);
	for (int i = 0; i < nlevels; ++i) {
		precon->nmg[i] = max(2, (n - 1) / int_pow(2, i) + 1);
		if (precon->nmg[i] == 2) {
			nlevels = i + 1;
			break;
		}
	}
	precon->nlevels = nlevels;
	precon->geom = new_ptr_array(nlevels, 0, NULL);
	precon->basis = new_ptr_array(2 * nlevels - 1, 0, NULL);
	precon->poisson_diag = new_ptr_array(nlevels, 0, NULL);
	precon->helmholtz_diag = new_ptr_array(nlevels, 0, NULL);
	precon->stiffsum_scale = new_ptr_array(nlevels, 0, NULL);
	precon->dirichlet_mask = new_ptr_array(nlevels, 0, NULL);
	precon->sendcounts = new_ptr_array(nlevels, 0, NULL);
	precon->senddispls = new_ptr_array(nlevels, 0, NULL);
	precon->poisson_eigen_max = new_1d_array(nlevels);
	precon->helmholtz_eigen_max = new_1d_array(nlevels);

	return precon;
}

void free_preconditioner(struct preconditioner *precon, int id)
{
	if (id == 0) {
		free_3d_array(precon->inverse_mass);
		for (int i = 0; i < precon->nlevels; ++i) {
			free_3d_array(precon->poisson_diag[i]);
			free_3d_array(precon->helmholtz_diag[i]);
			free_3d_array(precon->stiffsum_scale[i]);
			free_3d_array(precon->dirichlet_mask[i]);
		}
	}
	for (int i = 1; i < precon->nlevels; ++i) {
		if (precon->geom[i])
			free_geometry(precon->geom[i]);
		if (precon->sendcounts[i])
			free_int_array(precon->sendcounts[i]);
		if (precon->senddispls[i])
			free_int_array(precon->senddispls[i]);
	}
	for (int i = 1; i < (2 * precon->nlevels - 1); ++i)
		if (precon->basis[i])
			free_basis(precon->basis[i]);
	free_ptr_array(precon->basis, 0, NULL);
	free_ptr_array(precon->geom, 0, NULL);
	free_ptr_array(precon->poisson_diag, 0, NULL);
	free_ptr_array(precon->helmholtz_diag, 0, NULL);
	free_ptr_array(precon->stiffsum_scale, 0, NULL);
	free_ptr_array(precon->dirichlet_mask, 0, NULL);
	free_ptr_array(precon->sendcounts, 0, NULL);
	free_ptr_array(precon->senddispls, 0, NULL);
	free_int_array(precon->nmg);
	free_1d_array(precon->poisson_eigen_max);
	free_1d_array(precon->helmholtz_eigen_max);
	free(precon);
}

static void initialise_multigrid_precon(struct element **sel)
{
	const int nlevels = (*sel)->precon->nlevels;
	for (int i = 0; i < (*sel)->nel; ++i) {
		sel[i]->precon->basis[0] = sel[i]->basis;
		sel[i]->precon->geom[0] = sel[i]->geom;
		for (int j = 1; j < 2 * nlevels - 1; ++j) {
			const int n = sel[i]->precon->nmg[j / 2];
			sel[i]->precon->basis[j] = new_gll_basis(n);
			compute_svv_derivative_matrix(sel[i]->precon->basis[j],
					sel[i]->params->svvn,
					sel[i]->params->svva);
		}
		for (int j = 2; j < 2 * nlevels - 1; j += 2)
			create_interpolation_matrices(sel[i]->precon->basis[j - 1],
					sel[i]->precon->basis[j]);
		for (int j = 1; j < nlevels; ++j) {
			struct basis *basis1 = sel[i]->precon->basis[j * 2];
			struct basis *basis2 = sel[i]->precon->basis[j * 2 - 1];
			sel[i]->precon->geom[j] = new_geometry(basis1,
					(*sel)->nel, (*sel)->geom->neltot);
			copy_geometry_connectivity(sel[i]->geom, sel[i]->precon
					->geom[j]);
			struct trivec *grid = sel[i]->precon->geom[j]->grid;
			interpolate_to_grid(basis1, grid->x, basis2,
					sel[i]->precon->geom[j - 1]->grid->x);
			interpolate_to_grid(basis1, grid->y, basis2,
					sel[i]->precon->geom[j - 1]->grid->y);
			interpolate_to_grid(basis1, grid->z, basis2,
					sel[i]->precon->geom[j - 1]->grid->z);
			compute_geometry_transformations(sel[i]->precon
					->geom[j], basis1);
		}
	}
	const int n = (*sel)->precon->nmg[0];
	double ****tmp = new_4d_array(n, n, n, (*sel)->nel);
	for (int i = 0; i < (*sel)->nel; ++i)
		sel[i]->precon->inverse_mass = tmp[i];
	free(tmp);
	for (int j = 0; j < nlevels; ++j) {
		const int n = sel[0]->precon->nmg[j];
		double ****tmp1 = new_4d_array(n, n, n, (*sel)->nel);
		double ****tmp2 = new_4d_array(n, n, n, (*sel)->nel);
		double ****tmp3 = new_4d_array(n, n, n, (*sel)->nel);
		double ****tmp4 = new_4d_array(n, n, n, (*sel)->nel);
		for (int i = 0; i < (*sel)->nel; ++i) {
			sel[i]->precon->poisson_diag[j] = tmp1[i];
			sel[i]->precon->helmholtz_diag[j] = tmp2[i];
			sel[i]->precon->stiffsum_scale[j] = tmp3[i];
			sel[i]->precon->dirichlet_mask[j] = tmp4[i];
		}
		free(tmp1);
		free(tmp2);
		free(tmp3);
		free(tmp4);
	}
	for (int j = 1; j < nlevels; ++j) {
		(*sel)->precon->sendcounts[j] = new_int_array((*sel)->params
				->nproc);
		(*sel)->precon->senddispls[j] = new_int_array((*sel)->params
				->nproc);
	}
}

static void compute_multigrid_sendcounts(struct element **sel)
{
	(*sel)->precon->sendcounts[0] = (*sel)->comm->sendcounts;
	(*sel)->precon->senddispls[0] = (*sel)->comm->senddispls;
	for (int j = 1; j < (*sel)->precon->nlevels; ++j) {
		(*sel)->comm->sendcounts = (*sel)->precon->sendcounts[j];
		(*sel)->comm->senddispls = (*sel)->precon->senddispls[j];
		(*sel)->basis = (*sel)->precon->basis[(j + 1) / 2 + 1];
		compute_buffer_sizes(sel);
	}
	(*sel)->comm->sendcounts = (*sel)->precon->sendcounts[0];
	(*sel)->comm->senddispls = (*sel)->precon->senddispls[0];
	(*sel)->basis = (*sel)->precon->basis[0];
}

void initialise_preconditioners(struct element **sel)
{
	initialise_multigrid_precon(sel);
	compute_multigrid_sendcounts(sel);
	initialise_diagonal_precon(sel, PRESSURE);
	initialise_diagonal_precon(sel, X_VELOCITY);
	initialise_stiffness_summation_scaling(sel);
	initialise_dirichlet_mask(sel);
	initialise_inverse_mass_matrix(sel);
}

static double ***precon_array(struct element *sel, int var, int level)
{
	switch (var) {
		case PRESSURE:
			return sel->precon->poisson_diag[level];
		case X_VELOCITY:
			return sel->precon->helmholtz_diag[level];
		case Y_VELOCITY:
			return sel->precon->helmholtz_diag[level];
		case Z_VELOCITY:
			return sel->precon->helmholtz_diag[level];
		default:
			return NULL;
	}
}

void multigrid_preconditioner(struct element **sel, double ****x, double ****px,
		int var, double eps, void (*opeval) ())
{
	jacobi_preconditioner(sel, x, px, var, eps, opeval);
	for (int i = 0; i < 1; ++i) {
		double err = multigrid_recurse(sel, x, px, var, 0, eps, opeval);
		if (err < eps)
			break;
	}
}

static double multigrid_recurse(struct element **sel, double ****f,
		double ****x, int var, int level, double eps,
		void (*opeval) ())
{
	if (level == (*sel)->precon->nlevels - 1)
		return conjugate_gradient(sel, f, x, var, opeval,
				jacobi_preconditioner, eps, 500);

	const int m = 15;
	const int n1 = (*sel)->basis->n;
	const int n2 = (*sel)->precon->nmg[level + 1];
	double ****r = new_4d_array(n1, n1, n1, (*sel)->nel);
	double ****fn = new_4d_array(n2, n2, n2, (*sel)->nel);
	double ****e = new_4d_array(n2, n2, n2, (*sel)->nel);

	if (level == 0)
		compute_residual(sel, x, f, r, var, opeval);
	else
		copy_array(***f, ***r, n1 * n1 * n1 * (*sel)->nel);

	switch_multigrid_level(sel, level + 1);
	multigrid_restriction(sel, r, fn, level + 1);
	multigrid_recurse(sel, fn, e, var, level + 1, eps, opeval);
	switch_multigrid_level(sel, level);
	multigrid_prolongation(sel, e, r, level);
	add_array(***r, ***x, n1 * n1 * n1 * (*sel)->nel, 1.0);

	double err = conjugate_gradient(sel, f, x, var, opeval,
			jacobi_preconditioner, eps, m);
	free_4d_array(r);
	free_4d_array(fn);
	free_4d_array(e);
	return err;
}

static void switch_multigrid_level(struct element **sel, int level)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		sel[i]->precon->level = level;
		sel[i]->basis = sel[i]->precon->basis[2 * level];
		sel[i]->geom = sel[i]->precon->geom[level];
	}
	reconnect_elements(sel);
}

static void multigrid_restriction(struct element **sel, double ****x,
		double ****y, int level)
{
	struct basis *xbasis = (*sel)->precon->basis[level * 2 - 1];
	struct basis *ybasis = (*sel)->precon->basis[level * 2];
	for (int i = 0; i < (*sel)->nel; ++i)
		integrate_to_grid(ybasis, y[i], xbasis, x[i]);
	direct_stiffness_averaging(sel, y);
}

static void multigrid_prolongation(struct element **sel, double ****x,
		double ****y, int level)
{
	struct basis *xbasis = (*sel)->precon->basis[level * 2 + 2];
	struct basis *ybasis = (*sel)->precon->basis[level * 2 + 1];
	for (int i = 0; i < (*sel)->nel; ++i)
		interpolate_to_grid(ybasis, y[i], xbasis, x[i]);
}

static void compute_residual(struct element **sel, double ****x, double ****f,
	double ****r, int var, void (*opeval) ())
{
	const int n = (*sel)->basis->ntot;
	global_opeval(sel, x, r, var, opeval);
	scale_array(***r, n * (*sel)->nel, -1.0);
	add_array(***f, ***r, n * (*sel)->nel, 1.0);
}

void jacobi_preconditioner(struct element **sel, double ****x, double ****px,
		int var, double eps, void (*opeval) ())
{
	diagonal_preconditioner(sel, x, px, var, opeval, 1.0, 0.0);
}

static void diagonal_preconditioner(struct element **sel, double ****x,
		double ****px, int var, void (*opeval) (), double omega,
		double alpha)
{
	const size_t n = (*sel)->basis->ntot * (*sel)->nel;
	const int k = (*sel)->precon->level;
	multiply(***x, **precon_array(*sel, var, k), ***px, n, omega, alpha);
}

static void reconnect_elements(struct element **sel)
{
	const int level = (*sel)->precon->level;
	(*sel)->comm->sendcounts = (*sel)->precon->sendcounts[level];
	(*sel)->comm->senddispls = (*sel)->precon->senddispls[level];
	connect_global_elements(sel);
	connect_global_edges(sel);
	connect_global_corners(sel);
}

static void initialise_diagonal_precon(struct element **sel, int var)
{
	const double lambda = (11.0 / 6.0) * (double) (var != PRESSURE)
			/ ((*sel)->params->dt * (*sel)->params->visc);
	double ****p = new_ptr_array((*sel)->nel, 0, NULL);

	for (int j = 0; j < (*sel)->precon->nlevels; ++j) {
		switch_multigrid_level(sel, j);
		const int n = (*sel)->basis->n;
		double ****a = new_4d_array(n, n, n, (*sel)->nel);
		double ****b = new_4d_array(n, n, n, (*sel)->nel);
		for (int i = 0; i < (*sel)->nel; ++i) {
			p[i] = precon_array(sel[i], var, j);
			set_helmholtz_diagonal(sel[i], a[i], lambda);
		}
		direct_stiffness_summation(sel, a);
		ones(***b, (*sel)->basis->ntot * (*sel)->nel);
		divide(***b, ***a, ***p, n * n * n * (*sel)->nel, 1.0, 0.0);
		free_4d_array(a);
		free_4d_array(b);
	}
	free_ptr_array(p, 0, NULL);
	switch_multigrid_level(sel, 0);
}

static void initialise_stiffness_summation_scaling(struct element **sel)
{
	double ****x = new_ptr_array((*sel)->nel, 0, NULL);
	for (int j = 0; j < (*sel)->precon->nlevels; ++j) {
		switch_multigrid_level(sel, j);
		const int n = (*sel)->basis->n;
		double ****a = new_4d_array(n, n, n, (*sel)->nel);
		double ****b = new_4d_array(n, n, n, (*sel)->nel);
		for (int i = 0; i < (*sel)->nel; ++i)
			x[i] = sel[i]->precon->stiffsum_scale[j];

		ones(***a, (*sel)->basis->ntot * (*sel)->nel);
		ones(***b, (*sel)->basis->ntot * (*sel)->nel);
		direct_stiffness_summation(sel, a);
		divide(***b, ***a, ***x, n * n * n * (*sel)->nel, 1.0, 0.0);
		free_4d_array(a);
		free_4d_array(b);
	}
	free_ptr_array(x, 0, NULL);
	switch_multigrid_level(sel, 0);
}

static void initialise_dirichlet_mask(struct element **sel)
{
	for (int j = 0; j < (*sel)->precon->nlevels; ++j) {
		switch_multigrid_level(sel, j);
		ones(**(*sel)->precon->dirichlet_mask[j],
				(*sel)->basis->ntot * (*sel)->nel);
		for (int i = 0; i < (*sel)->nel; ++i)
			set_dirichlet_boundaries(sel[i], sel[i]->precon
					->dirichlet_mask[j], sel[i]->bc,
					zero_func);
	}
	switch_multigrid_level(sel, 0);
}

static void initialise_inverse_mass_matrix(struct element **sel)
{
	const int n = (*sel)->basis->n;
	double ****x = new_4d_array(n, n, n, (*sel)->nel);
	double ****y = new_4d_array(n, n, n, (*sel)->nel);

	ones(***x, n * n * n * (*sel)->nel);
	ones(***y, n * n * n * (*sel)->nel);
	for (int i = 0; i < (*sel)->nel; ++i)
		test_function_inner_product(sel[i], y[i], y[i]);
	direct_stiffness_summation(sel, y);
	divide(***x, ***y, **(*sel)->precon->inverse_mass,
			n * n * n * (*sel)->nel, 1.0, 0.0);
		
	free_4d_array(x);
	free_4d_array(y);
}
