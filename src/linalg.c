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
#include "utils.h"
#include "structs.h"
#include "lapack.h"
#include "blas.h"
#include "linalg.h"

void matrix_vector_product(double **a, double *x, double *y, int m, int n,
		double alpha, double beta)
{
	char trans = 'N';
	int incx = 1;
	int incy = 1;
	dgemv_(&trans, &m, &n, &alpha, *a, &m, x, &incx, &beta, y, &incy);
}

void matrix_matrix_product(double **a, double **b, double **c, int n1, int n2,
		int n3, double alpha, double beta)
{
	char transa = 'N';
	char transb = 'N';

	dgemm_(&transa, &transb, &n1, &n3, &n2, &alpha, *a, &n1, *b, &n2,
			&beta, *c, &n1);
}

void add_array(double *x, double *y, int n, double alpha)
{
	int incx = 1;
	int incy = 1;

	daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void add_inc_array(double *x, double *y, int n, double alpha, int incx,
		int incy)
{
	daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void copy_array(double *x, double *y, int n)
{
	int incx = 1;
	int incy = 1;

	dcopy_(&n, x, &incx, y, &incy);
}

void copy_int_array(int *x, int *y, int n)
{
	for (int i = 0; i < n; ++i)
		y[i] = x[i];
}

void copy_inc_array(double *x, double *y, int n, int incx, int incy)
{
	dcopy_(&n, x, &incx, y, &incy);
}

void scale_array(double *x, int n, double alpha)
{
	int incx = 1;

	dscal_(&n, &alpha, x, &incx);
}

double dot_product(double *x, double *y, int n)
{
	int incx = 1;
	int incy = 1;

	double prod = ddot_(&n, x, &incx, y, &incy);
	return prod;
}

void compute_eigenvalues(double **a, double *lambda, int n)
{
	char jobvl = 'N';
	char jobvr = 'N';
	int lda = n;
	int ldvl = 1;
	int ldvr = 1;
	double *wi = new_1d_array(n);
	double *vl = NULL;
	double *vr = NULL;
	int lwork = 3*n;
	double *work = new_1d_array(lwork);
	int info;

	dgeev_(&jobvl, &jobvr, &n, *a, &lda, lambda, wi, vl, &ldvl, vr, &ldvr,
			work, &lwork, &info);
	free_1d_array(wi);
	free_1d_array(work);
}

void generalised_eigenvalue_problem(double **a, double **b, int n,
		double *lambda)
{
	int itype = 1;
	char jobz = 'V';
	char uplo = 'U';
	int lwork = 3 * n;
	double *work = new_1d_array(lwork);
	int info;

	dsygv_(&itype, &jobz, &uplo, &n, *a, &n, *b, &n, lambda, work, &lwork,
			&info);
	free_1d_array(work);
}

void sort_array(double *x, int n, char inc_dec)
{
	int info;
	dlasrt_(&inc_dec, &n, x, &info);
}

void transpose(double **a, double **at, int n1, int n2)
{
	for (int j = 0; j < n2; ++j)
		for (int i = 0; i < n1; ++i)
			at[j][i] = a[i][j];
}

void multiply(const double * restrict x, const double * restrict y,
		double * restrict z, size_t n, double alpha, double beta)
{
	for (int i = 0; i < n; ++i) {
		const double a = x[i];
		const double b = y[i];
		const double c = z[i];
		const double d = alpha * a * b + beta * c;
		z[i] = d;
	}
}

void divide(const double * restrict x, const double * restrict y,
		double * restrict z, size_t n, double alpha, double beta)
{
	for (int i = 0; i < n; ++i) {
		const double a = x[i];
		const double b = y[i];
		const double c = z[i];
		const double d = alpha * a / b + beta * c;
		z[i] = d;
	}
}
