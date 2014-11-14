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
#include "basis.h"
#include "interpolation.h"

static double lagrange(double xi, double *x, int j, int n);
static int first_interpolation_point(double *x, double xp, int n, int m);
static void allocate_interpolation_matrices(struct basis *basis1,
		struct basis *basis2);
static void compute_interpolation_matrix(struct basis *to_basis,
		struct basis *from_basis);

static double lagrange(double xi, double *x, int j, int n)
{
	double li = 1.0;
	for (int i = 0; i < j; ++i)
		li *= (xi - x[i]) / (x[j] - x[i]);
	for (int i = j + 1; i < n; ++i)
		li *= (xi - x[i]) / (x[j] - x[i]);
	return li;
}

double lagrangian_interpolation(double *x, double ***f, double xp, double yp,
		double zp, int n)
{
	double f_int = 0.0;
	double *lx = new_1d_array(n);
	double *ly = new_1d_array(n);
	double *lz = new_1d_array(n);
	double **tmp = new_2d_array(n, n);

	for (int i = 0; i < n; ++i) {
		lx[i] = lagrange(xp, x, i, n);
		ly[i] = lagrange(yp, x, i, n);
		lz[i] = lagrange(zp, x, i, n);
	}

	matrix_vector_product(*f, lz, *tmp, n * n, n, ADD_NEW, REPLACE_OLD);
	matrix_vector_product(tmp, ly, lz, n, n, ADD_NEW, REPLACE_OLD);
	f_int = dot_product(lx, lz, n);

	free_1d_array(lx);
	free_1d_array(ly);
	free_1d_array(lz);
	free_2d_array(tmp);
	return f_int;
}

double polynomial_interpolation(double *x, double ***f, double xp, double yp,
		double zp, int n)
{
	const int m = 5;
	if (n <= 5)
		return lagrangian_interpolation(x, f, xp, yp, zp, n);

	double tmp = 0.0;
	int i1 = first_interpolation_point(x, xp, n, m);
	int j1 = first_interpolation_point(x, yp, n, m);
	int k1 = first_interpolation_point(x, zp, n, m);
	double *lx = new_1d_array(m);
	double *ly = new_1d_array(m);
	double *lz = new_1d_array(m);

	for (int i = 0; i < m; ++i) {
		lx[i] = lagrange(xp, x + i1, i, m);
		ly[i] = lagrange(yp, x + j1, i, m);
		lz[i] = lagrange(zp, x + k1, i, m);
	}
	for (int k = 0; k < m; ++k)
		for (int j = 0; j < m; ++j)
			for (int i = 0; i < m; ++i)
				tmp += f[k1 + k][j1 + j][i1 + i]
					* lx[i] * ly[j] * lz[k];
	free_1d_array(lx);
	free_1d_array(ly);
	free_1d_array(lz);
	return tmp;
}

static int first_interpolation_point(double *x, double xp, int n, int m)
{
	for (int i = 0; i < n - m + 1; ++i)
		if (x[i] > xp)
			return max(0, i - m / 2);
	return n - m;
}

static void allocate_interpolation_matrices(struct basis *basis1,
		struct basis *basis2)
{
	basis1->interp_matrix = new_2d_array(basis1->n, basis2->n);
	basis1->trans_interp_matrix = new_2d_array(basis2->n, basis1->n);
	basis2->interp_matrix = new_2d_array(basis2->n, basis1->n);
	basis2->trans_interp_matrix = new_2d_array(basis1->n, basis2->n);
}

static void compute_interpolation_matrix(struct basis *to_basis,
		struct basis *from_basis)
{
	for (int p = 0; p < from_basis->n; ++p)
		for (int i = 0; i < to_basis->n; ++i)
			to_basis->interp_matrix[p][i] = lagrange(
					to_basis->points[i],
					from_basis->points, p, from_basis->n);
}

void create_interpolation_matrices(struct basis *basis1, struct basis *basis2)
{
	allocate_interpolation_matrices(basis1, basis2);
	compute_interpolation_matrix(basis1, basis2);
	compute_interpolation_matrix(basis2, basis1);
	transpose(basis1->interp_matrix, basis1->trans_interp_matrix,
			basis2->n, basis1->n);
	transpose(basis2->interp_matrix, basis2->trans_interp_matrix,
			basis1->n, basis2->n);
}

void interpolate_to_grid(struct basis *to_basis, double ***int_u,
		struct basis *from_basis, double ***u)
{
	const int m = to_basis->n;
	const int n = from_basis->n;
	double ***tmp1 = new_3d_array(m, n, n);
	double ***tmp2 = new_3d_array(m, m, n);

	matrix_matrix_product(to_basis->interp_matrix, *u, *tmp1, m, n, n * n,
			1.0, 0.0);
	for (int k = 0; k < n; ++k)
		matrix_matrix_product(*(tmp1 + k),
				to_basis->trans_interp_matrix, *(tmp2 + k), m,
				n, m, 1.0, 0.0);
	matrix_matrix_product(*tmp2, to_basis->trans_interp_matrix, *int_u,
			m * m, n, m, 1.0, 0.0);

	free_3d_array(tmp1);
	free_3d_array(tmp2);
}

void integrate_to_grid(struct basis *to_basis, double ***int_u,
		struct basis *from_basis, double ***u)
{
	const int m = to_basis->n;
	const int n = from_basis->n;
	double ***tmp1 = new_3d_array(m, n, n);
	double ***tmp2 = new_3d_array(m, m, n);

	matrix_matrix_product(from_basis->trans_interp_matrix, *u, *tmp1, m, n,
			n * n, 1.0, 0.0);
	for (int k = 0; k < n; ++k)
		matrix_matrix_product(*(tmp1 + k),
				from_basis->interp_matrix, *(tmp2 + k), m, n, m,
				1.0, 0.0);
	matrix_matrix_product(*tmp2, from_basis->interp_matrix, *int_u, m * m,
			n, m, 1.0, 0.0);

	free_3d_array(tmp1);
	free_3d_array(tmp2);
}
