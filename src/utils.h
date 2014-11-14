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

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

#define PI 3.1415926535

double *new_1d_array(size_t n);
double **new_2d_array(size_t n1, size_t n2);
double ***new_3d_array(size_t n1, size_t n2, size_t n3);
double ****new_4d_array(size_t n1, size_t n2, size_t n3, size_t n4);
void *new_ptr_array(size_t nptr, size_t n, void *(*new_data) ());
int *new_int_array(size_t n);
void free_1d_array(double *array);
void free_2d_array(double **array);
void free_3d_array(double ***array);
void free_4d_array(double ****array);
void free_ptr_array(void *ptr_array, size_t nptr, void (*free_data) ());
void free_int_array(int *array);
void *ec_malloc(size_t size);
void *ec_realloc(void *old_ptr, size_t size);
void fatal_error(char *message);
void zeros(double *array, size_t n);
void ones(double *array, size_t n);
int int_pow(int x, int p);
int min(int a, int b);
int max(int a, int b);
double dmin(double a, double b);
double dmax(double a, double b);
double random_double(void);
int id_offset(int proc, int neltot, int nproc);
int local_nel(int proc, int neltot, int nproc);
int id_of_elnum(int elnum, int proc, int neltot, int nproc);
int elnum_of_id(int id, int neltot, int nproc);
int proc_of_id(int id, int neltot, int nproc);
#endif
