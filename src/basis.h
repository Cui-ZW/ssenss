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

#ifndef BASIS_H
#define BASIS_H

struct basis {
	int n;
	int ntot;
	double *points;
	double *weights;
	double **deriv_matrix;
	double **svv_deriv_matrix;
	double **trans_deriv_matrix;
	double ***weights_array;
	double **interp_matrix;
	double **trans_interp_matrix;
	double **q;
	double **qt;
	double *lambda;
};

struct basis *new_gll_basis(int n);
void free_basis(struct basis *discrete_basis);
double *gll_gridpoints(int n);
double *gll_weights(const double *xi, int n);
double **lagrange_derivative_matrix(const double *xi, int n);
double ***gll_weights_array(const double *rho, int n);
double **legendre_transform_matrix(const double *xi, const double *rho, int n);
double **inverse_legendre_transform_matrix(const double *xi, int n);
double **stiffness_matrix_1d(const double *rho, double **d, int n);
double **mass_matrix_1d(const double *rho, int n);

#endif
