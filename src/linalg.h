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

#ifndef LINALG_H
#define LINALG_H

#include <stdio.h>

#define SORT_INC 'I'
#define SORT_DEC 'D'
#define ADD_NEW 1.0
#define SUBTRACT_NEW -1.0
#define ADD_OLD 1.0
#define SUBTRACT_OLD -1.0
#define REPLACE_OLD 0.0

void matrix_vector_product(double **a, double *x, double *y, int m, int n,
		double alpha, double beta);
void matrix_matrix_product(double **a, double **b, double **c, int n1, int n2,
		int n3, double alpha, double beta);
void add_array(double *x, double *y, int n, double alpha);
void add_inc_array(double *x, double *y, int n, double alpha, int incx,
		int incy);
void copy_array(double *x, double *y, int n);
void copy_inc_array(double *x, double *y, int n, int incx, int incy);
void copy_int_array(int *x, int *y, int n);
void scale_array(double *x, int n, double alpha);
double dot_product(double *x, double *y, int n);
void compute_eigenvalues(double **a, double *lambda, int n);
void generalised_eigenvalue_problem(double **a, double **b, int n,
		double *lambda);
void sort_array(double *x, int n, char inc_dec);
void transpose(double **a, double **at, int n1, int n2);
void multiply(const double * restrict x, const double * restrict y,
		double * restrict z, size_t n, double alpha, double beta);
void divide(const double * restrict x, const double *restrict y,
		double * restrict z, size_t n, double alpha, double beta);

#endif
