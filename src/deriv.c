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
#include "linalg.h"
#include "tripalg.h"
#include "basis.h"
#include "deriv.h"

static double svv_func(int n, int m, int i);

void discrete_x_differentiation(struct basis *disc_basis, double ***u,
		double ***du, double alpha, double beta)
{
	const int n = disc_basis->n;
	matrix_matrix_product(disc_basis->deriv_matrix, *u, *du, n, n, n * n,
			alpha, beta);
}

void discrete_y_differentiation(struct basis *disc_basis, double ***u,
		double ***du, double alpha, double beta)
{
	const int n = disc_basis->n;
	for (int k = 0; k < n; ++k)
		matrix_matrix_product(*(u + k), disc_basis->trans_deriv_matrix,
				*(du + k), n, n, n, alpha, beta);
}

void discrete_z_differentiation(struct basis *disc_basis, double ***u,
		double ***du, double alpha, double beta)
{
	const int n = disc_basis->n;
	matrix_matrix_product(*u, disc_basis->trans_deriv_matrix, *du, n * n,
			n, n, alpha, beta);
}

void transpose_divergence(struct basis *disc_basis, struct trivec *u,
		double ***du)
{
	const int n = disc_basis->n;

	matrix_matrix_product(disc_basis->trans_deriv_matrix, *u->x, *du, n, n,
			n * n, ADD_NEW, REPLACE_OLD);
	for (int k = 0; k < n; ++k)
		matrix_matrix_product(*(u->y + k), disc_basis->deriv_matrix,
				*(du + k), n, n, n, ADD_NEW, ADD_OLD);
	matrix_matrix_product(*u->z, disc_basis->deriv_matrix, *du, n * n, n, n,
			ADD_NEW, ADD_OLD);
}

void compute_jacobian(struct basis *disc_basis, struct trivec *grid,
		struct trimat *jac)
{
	for (int i = 0; i < 3; ++i) {
		discrete_x_differentiation(disc_basis, *grid->ind[i],
				*jac->ind[i][0], ADD_NEW, REPLACE_OLD);
		discrete_y_differentiation(disc_basis, *grid->ind[i],
				*jac->ind[i][1], ADD_NEW, REPLACE_OLD);
		discrete_z_differentiation(disc_basis, *grid->ind[i],
				*jac->ind[i][2], ADD_NEW, REPLACE_OLD);
	}
}

void compute_svv_derivative_matrix(struct basis *basis, int m, double c)
{
	const double *xi = basis->points;
	const double *rho = basis->weights;
	int n = basis->n;
	double **leg = legendre_transform_matrix(xi, rho, n);
	double **inv_leg = inverse_legendre_transform_matrix(xi, n);
	double **d = lagrange_derivative_matrix(xi, n);
	double **tmp = new_2d_array(n, n);

	matrix_matrix_product(leg, d, tmp, n, n, n, ADD_NEW, REPLACE_OLD);

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			tmp[j][i] *= sqrt(1.0 + c * svv_func(n - 1, m, i));

	matrix_matrix_product(inv_leg, tmp, basis->svv_deriv_matrix, n, n, n,
			ADD_NEW, REPLACE_OLD);

	free_2d_array(tmp);
	free_2d_array(d);
	free_2d_array(leg);
	free_2d_array(inv_leg);
}

static double svv_func(int n, int m, int i)
{
	if (i <= m)
		return 0.0;
	else
		return exp(-pow(((double) (n - i)) / ((double) (m - i)), 2));
}
