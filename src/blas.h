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

#ifndef BLAS_H
#define BLAS_H

void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
		int *lda, double *x, int *incx, double *beta, double *y,
		int *incy);
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
		double *a, int *lda, double *b, int *ldb, double *beta,
		double *c, int *ldc);
void daxpy_(int *n, double *a, double *x, int *incx, double *y, int *incy);
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
void dscal_(int *n, double *a, double *x, int *incx);
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void dsbmv_(char *uplo, int *n, int *k, double *alpha, double *a, int *lda,
		double *x, int *incx, double *beta, double *y, int *incy);

#endif
