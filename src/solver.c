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

#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "basis.h"
#include "linalg.h"
#include "poisson.h"
#include "helmholtz.h"
#include "global.h"
#include "precon.h"
#include "test.h"
#include "solver.h"

#define BC sel[i]->bc * (var != PRESSURE)

double conjugate_gradient(struct element **sel, double ****f, double ****u,
		const int var, void (*opeval) (), void (*precon) (),
		const double eps, int max_iter)
{
	const int n = (*sel)->basis->n;
	const int nel = (*sel)->nel;
	const size_t ntot = n * n * n;

	double ****r = new_4d_array(n, n, n, nel);
	double ****z = new_4d_array(n, n, n, nel);
	double ****p = new_4d_array(n, n, n, nel);
	double ****ap = new_4d_array(n, n, n, nel);

	mask_dirichlet_boundaries(sel, f, var);
	global_opeval(sel, u, r, var, opeval);
	scale_array(***r, ntot * nel, -1.0);
	add_array(***f, ***r, ntot * nel, 1.0);
	subtract_nullspace(sel, r, var);

	precon(sel, r, z, var, eps, opeval);
	subtract_nullspace(sel, z, var);
	copy_array(***z, ***p, ntot * nel);
	double rdz = global_dot_product(sel, r, z);
	size_t j = 0;
	double err;
	do {
		global_opeval(sel, p, ap, var, opeval);
		const double alpha = rdz / global_dot_product(sel, p, ap);
		add_array(***p, ***u, nel * ntot, alpha);
		add_array(***ap, ***r, nel * ntot, -alpha);

		precon(sel, r, z, var, eps, opeval);
		subtract_nullspace(sel, z, var);

		const double rdz_old = rdz;
		rdz = global_dot_product(sel, r, z);
		const double beta = rdz / rdz_old;
		scale_array(***p, nel * ntot, beta);
		add_array(***z, ***p, nel * ntot, ADD_NEW);
		mask_dirichlet_boundaries(sel, p, var);
	} while ((err = sqrt(global_dot_product(sel, r, r))) > eps
			&& ++j < max_iter);

	free_4d_array(p);
	free_4d_array(ap);
	free_4d_array(r);
	free_4d_array(z);
	return err;
}

double flexible_conjugate_gradient(struct element **sel, double ****f,
		double ****u, const int var, void (*opeval) (),
		void (*precon) (), const double eps, int max_iter)
{
	const int n = (*sel)->basis->n;
	const int nel = (*sel)->nel;
	const size_t ntot = n * n * n;
	int mmax = 1;
	if (precon == multigrid_preconditioner)
		mmax = 10;

	double ****r = new_4d_array(n, n, n, nel);
	double ****z = new_4d_array(n, n, n, nel);
	double *****p = new_ptr_array(mmax + 1, 0, NULL);
	double *****ap = new_ptr_array(mmax, 0, NULL);

	for (int i = 0; i < mmax; ++i)
		ap[i] = new_4d_array(n, n, n, nel);
	for (int i = 0; i < mmax + 1; ++i)
		p[i] = new_4d_array(n, n, n, nel);

	mask_dirichlet_boundaries(sel, f, var);
	global_opeval(sel, u, r, var, opeval);
	scale_array(***r, ntot * nel, -1.0);
	add_array(***f, ***r, ntot * nel, 1.0);
	subtract_nullspace(sel, r, var);

	size_t j = 0;
	double *pap = new_1d_array(mmax);
	double *beta = new_1d_array(mmax);
	double err;
	do {
		int m = min(j, max(1, j % (mmax + 1)));
		precon(sel, r, z, var, eps, opeval);
		subtract_nullspace(sel, z, var);

		for (int k = 0; k < m; ++k)
			beta[k] = global_dot_product(sel, z, ap[k]) / pap[k];
		for (int k = 0; k < m; ++k)
			copy_array(***p[m - k - 1], ***p[m - k], ntot * nel);

		copy_array(***z, ***p[0], ntot * nel);
		for (int k = 1; k < m + 1; ++k)
			add_array(***p[k], ***p[0], ntot * nel, -beta[k - 1]);
		mask_dirichlet_boundaries(sel, p[0], var);

		m = min(j + 1, max(1, (j + 1) % (mmax + 1)));
		for (int k = 0; k < m - 1; ++k)
			copy_array(***ap[m - k - 2], ***ap[m - k - 1],
					nel * ntot);

		global_opeval(sel, p[0], ap[0], var, opeval);

		for (int k = 0; k < m - 1; ++k)
			pap[m - k - 1] = pap[m - k - 2];
		pap[0] = global_dot_product(sel, p[0], ap[0]);
		const double alpha = global_dot_product(sel, p[0], r) / pap[0];

		add_array(***p[0], ***u, ntot * nel, alpha);
		add_array(***ap[0], ***r, ntot * nel, -alpha);
	} while ((err = sqrt(global_dot_product(sel, r, r))) > eps
			&& ++j < max_iter);
	free_1d_array(pap);
	free_1d_array(beta);

	for (int i = 0; i < mmax; ++i)
		free_4d_array(ap[i]);
	for (int i = 0; i < mmax + 1; ++i)
		free_4d_array(p[i]);
	free_ptr_array(p, mmax + 1, NULL);
	free_ptr_array(ap, mmax, NULL);
	free_4d_array(r);
	free_4d_array(z);
	return err;
}

void solve_poisson_cg(struct element **sel, double ****f, int var, double eps)
{
	const int n = (*sel)->basis->n;
	double ****u = new_4d_array(n, n, n, (*sel)->nel);

	for (int i = 0; i < (*sel)->nel; ++i)
		copy_array(**variable_array(sel[i], var), **u[i],
				(*sel)->basis->ntot);
	flexible_conjugate_gradient(sel, f, u, var, discrete_laplace_operator,
			multigrid_preconditioner, eps, 5000);
	for (int i = 0; i < (*sel)->nel; ++i)
		copy_array(**u[i], **variable_array(sel[i], var),
				(*sel)->basis->ntot);
	free_4d_array(u);
}

void solve_helmholtz_cg(struct element **sel, double ****f, int var, double eps)
{
	const int n = (*sel)->basis->n;
	double ****u = new_4d_array(n, n, n, (*sel)->nel);

	for (int i = 0; i < (*sel)->nel; ++i)
		copy_array(**variable_array(sel[i], var), **u[i],
				(*sel)->basis->ntot);
	conjugate_gradient(sel, f, u, var, discrete_helmholtz_operator,
			jacobi_preconditioner, eps, 5000);
	for (int i = 0; i < (*sel)->nel; ++i)
		copy_array(**u[i], **variable_array(sel[i], var),
				(*sel)->basis->ntot);
	free_4d_array(u);
}

void global_opeval(struct element **sel, double ****p, double ****ap, int var,
		void (*opeval) ())
{
	for (int i = 0; i < (*sel)->nel; ++i)
		opeval(sel[i], p[i], ap[i]);
	mask_dirichlet_boundaries(sel, ap, var);
	direct_stiffness_summation(sel, ap);
}

void subtract_nullspace(struct element **sel, double ****u, int var)
{
	if (var != PRESSURE)
		return;
	const int n = (*sel)->basis->n;

	double ****null_vec = new_4d_array(n, n, n, (*sel)->nel);
	ones(***null_vec, n * n * n * (*sel)->nel);
	const double nunique = global_dot_product(sel, null_vec, null_vec);
	const double null_comp = global_dot_product(sel, u, null_vec);
	add_array(***null_vec, ***u, n * n * n * (*sel)->nel,
			-null_comp / nunique);
	free_4d_array(null_vec);
}
