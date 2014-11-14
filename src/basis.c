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
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "structs.h"
#include "linalg.h"
#include "basis.h"

static double legendre(double x, int n);
static double legendre_derivative(double x, int n);
static double *gl_gridpoints(int n);

struct basis *new_gll_basis(int n)
{
	struct basis *gll_basis = ec_malloc(sizeof(*gll_basis));

	gll_basis->n = n;
	gll_basis->ntot = n * n * n;
	gll_basis->points = gll_gridpoints(n);
	gll_basis->weights = gll_weights(gll_basis->points, n);
	gll_basis->deriv_matrix = lagrange_derivative_matrix(
			gll_basis->points, n);
	gll_basis->trans_deriv_matrix = new_2d_array(n, n);
	transpose(gll_basis->deriv_matrix, gll_basis->trans_deriv_matrix, n, n);
	gll_basis->weights_array = gll_weights_array(gll_basis->weights, n);
	gll_basis->interp_matrix = NULL;
	gll_basis->trans_interp_matrix = NULL;
	gll_basis->q = NULL;
	gll_basis->qt = NULL;
	gll_basis->lambda = NULL;
	gll_basis->svv_deriv_matrix = new_2d_array(n, n);

	return gll_basis;
}

void free_basis(struct basis *discrete_basis)
{
	free_1d_array(discrete_basis->points);
	free_1d_array(discrete_basis->weights);
	free_2d_array(discrete_basis->deriv_matrix);
	free_2d_array(discrete_basis->trans_deriv_matrix);
	free_3d_array(discrete_basis->weights_array);
	free_2d_array(discrete_basis->svv_deriv_matrix);
	if (discrete_basis->interp_matrix)
		free_2d_array(discrete_basis->interp_matrix);
	if (discrete_basis->trans_interp_matrix)
		free_2d_array(discrete_basis->trans_interp_matrix);
	if (discrete_basis->q)
		free_2d_array(discrete_basis->q);
	if (discrete_basis->qt)
		free_2d_array(discrete_basis->qt);
	if (discrete_basis->lambda)
		free_1d_array(discrete_basis->lambda);
	free(discrete_basis);
}

static double legendre(double x, int n)
{
	double ln[3] = { 1.0, x, 0.0 };

	if (n <= 1)
		return ln[n];
	for (int i = 1; i < n; ++i) {
		ln[2] = ((2.0 * (double) i + 1.0) / ((double) i + 1.0))
				* x * ln[1] - (((double) i) /
				((double) i + 1.0)) * ln[0];
		ln[0] = ln[1];
		ln[1] = ln[2];
	}
	return ln[2];
}

static double legendre_derivative(double x, int n)
{
	double dln;

	if (pow(x, 2) == 1.0)
		dln = pow(x, n-1) * 0.5 * ((double) n) * ((double) (n+1));
	else
		dln = (x * legendre(x, n) - legendre(x, n - 1)) * ((double) n)
				/ (pow(x, 2) - 1.0);
	return dln;
}

static double *gl_gridpoints(int n)
{
	double **a = new_2d_array(n, n);
	double *xi = new_1d_array(n);

	a[1][0] = 1.0;
	if (n > 2)
		for (int i = 2; i < n; ++i) {
			a[i - 2][i - 1] = ((double) (i - 1)) / 
					((double) (2 * i - 1));
			a[i][i - 1] = ((double) i) / ((double) (2 * i - 1));
		}
	a[n-2][n-1] = ((double) (n - 1)) / ((double) (2 * n - 1));
	
	compute_eigenvalues(a, xi, n);
	sort_array(xi, n, SORT_INC);
	free_2d_array(a);
	return xi;
}

double *gll_gridpoints(int n)
{
	const double eps = 1.e-12;
	double *xi = new_1d_array(n);

	xi[0] = -1.0;
	xi[n-1] = 1.0;
	if (n < 3)
		return xi;
	double *gl = gl_gridpoints(n - 1);

	for (int i = 2; i < n; ++i) {
		xi[i - 1] = 0.5 * (gl[i - 2] + gl[i - 1]);
		double tmp = 0.0;

		while (fabs(xi[i - 1] - tmp) > eps) {
			tmp = xi[i - 1];
			xi[i - 1] = tmp + ((1.0 - pow(tmp, 2)) *
					legendre_derivative(tmp, n - 1)) /
					(((double) (n - 1)) * ((double) n) *
					legendre(tmp, n - 1));
		}
	}
	free_1d_array(gl);
	return xi;
}

double *gll_weights(const double *xi, int n)
{
	double *rho = new_1d_array(n);

	for (int i = 1; i <= n; ++i)
		rho[i - 1] = 2.0 / (((double) (n - 1)) * ((double) n) *
				pow(legendre(xi[i - 1], n - 1), 2));
	return rho;
}

double **lagrange_derivative_matrix(const double *xi, int n)
{
	double **d = new_2d_array(n, n);

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i) {
			if (i == j)
				d[j][i] = 0.0;
			else
				d[j][i] = legendre(xi[i], n - 1) /
						(legendre(xi[j], n - 1) *
						(xi[i] - xi[j]));
		}
	d[0][0] = -((double) n) * ((double) (n-1)) / 4.0;
	d[n - 1][n - 1] = -d[0][0];

	return d;
}

double ***gll_weights_array(const double *rho, int n)
{
	double ***v = new_3d_array(n, n, n);

	for (int k = 0; k < n; ++k)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i)
				v[k][j][i] = rho[i] * rho[j] * rho[k];
	return v;
}

double **legendre_transform_matrix(const double *xi, const double *rho, int n)
{
	double **a = new_2d_array(n, n);

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			a[j][i] = rho[j] * legendre(xi[j], i)
				* ((double) i + 0.5);
	return a;
}

double **inverse_legendre_transform_matrix(const double *xi, int n)
{
	double **a = new_2d_array(n, n);

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			a[j][i] = legendre(xi[i], j);
	return a;
}

double **stiffness_matrix_1d(const double *rho, double **d, int n)
{
	double **a = new_2d_array(n, n);

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			for (int k = 0; k < n; ++k)
				a[j][i] += rho[k] * d[i][k] * d[j][k];
	return a;
}

double **mass_matrix_1d(const double *rho, int n)
{
	double **b = new_2d_array(n, n);
	for (int i = 0; i < n; ++i)
		b[i][i] = rho[i];
	return b;
}
