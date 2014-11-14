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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "structs.h"
#include "final.h"

double *new_1d_array(size_t n)
{
	double *array = ec_malloc(n * sizeof(*array));
	zeros(array, n);

	return array;
}

double **new_2d_array(size_t n1, size_t n2)
{
	double **array;
	array = ec_malloc(n2 * sizeof(*array));
	*array = ec_malloc(n1 * n2 * sizeof(**array));
	
	for (size_t i = 1; i < n2; ++i)
		array[i] = array[i - 1] + n1;
	zeros(*array, n1 * n2);

	return array;
}

double ***new_3d_array(size_t n1, size_t n2, size_t n3)
{
	double ***array;
	array = ec_malloc(n3 * sizeof(*array));
	*array = ec_malloc(n2 * n3 * sizeof(**array));
	**array = ec_malloc(n1 * n2 * n3 * sizeof(***array));

	for (size_t i = 1; i < n3; ++i)
		array[i] = array[i - 1] + n2;
	for (size_t i = 1; i < n2 * n3; ++i)
		array[0][i] = array[0][i - 1] + n1;
	zeros(**array, n1 * n2 * n3);

	return array;
}

double ****new_4d_array(size_t n1, size_t n2, size_t n3, size_t n4)
{
	double ****array;
	array = ec_malloc(n4 * sizeof(*array));
	*array = ec_malloc(n3 * n4 * sizeof(**array));
	**array = ec_malloc(n2 * n3 * n4 * sizeof(***array));
	***array = ec_malloc(n1 * n2 * n3 * n4 * sizeof(****array));

	for (size_t i = 1; i < n4; ++i)
		array[i] = array[i - 1] + n3;
	for (size_t i = 1; i < n3 * n4; ++i)
		array[0][i] = array[0][i - 1] + n2;
	for (size_t i = 1; i < n2 * n3 * n4; ++i)
		array[0][0][i] = array[0][0][i - 1] + n1;
	zeros(***array, n1 * n2 * n3 * n4);

	return array;
}

void *new_ptr_array(size_t nptr, size_t n, void *(*new_data) ())
{
	void **ptr_array = ec_malloc(nptr * sizeof(ptr_array));
	for (size_t i = 0; i < nptr; ++i)
		if (new_data == NULL)
			ptr_array[i] = NULL;
		else if (new_data == (void *(*) ()) new_3d_array)
			ptr_array[i] = new_data(n, n, n);
		else
			ptr_array[i] = new_data(n);
	return ptr_array;
}

int *new_int_array(size_t n)
{
	int *array = ec_malloc(n * sizeof(*array));

	for (size_t i = 0; i < n; ++i)
		array[i] = 0;
	return array;
}

void free_1d_array(double *array)
{
	free(array);
}

void free_2d_array(double **array)
{
	free(*array);
	free(array);
}

void free_3d_array(double ***array)
{
	free(**array);
	free(*array);
	free(array);
}

void free_4d_array(double ****array)
{
	free(***array);
	free(**array);
	free(*array);
	free(array);
}

void free_ptr_array(void *ptr_array, size_t nptr, void (*free_data) ())
{
	for (size_t i = 0; i < nptr; ++i)
		if (free_data == NULL)
			break;
		else
			free_data(((void **) ptr_array)[i]);
	free(ptr_array);
}

void free_int_array(int *array)
{
	free(array);
}

void *ec_malloc(size_t size)
{
	void *ptr = malloc(size);
	if (ptr == NULL)
		fatal_error("in ec_malloc() on memory allocation");
	return ptr;
}

void *ec_realloc(void *old_ptr, size_t size)
{
	void *ptr = realloc(old_ptr, size);
	if (ptr == NULL)
		fatal_error("in ec_realloc() on memory allocation");
	return ptr;
}

void fatal_error(char *message)
{
	char error_message[100];
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	strcpy(error_message, "[!!] Fatal Error ");
	strncat(error_message, message, 82);
	if (rank == 0)
		fprintf(stderr, "%s\n", error_message);
	stop_mpi();	
	exit(0);
}

void zeros(double *array, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		array[i] = 0.0;
}

void ones(double *array, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		array[i] = 1.0;
}

int int_pow(int x, int p)
{
	if (p == 0)
		return 1;
	if (p == 1)
		return x;
	if (p < 0)
		return 0;

	const int tmp = int_pow(x, p / 2);
	if (p % 2 == 0)
		return tmp * tmp;
	else
		return x * tmp * tmp;
}

int min(int a, int b)
{
	if (a < b)
		return a;
	else
		return b;
}

int max(int a, int b)
{
	if (a > b)
		return a;
	else
		return b;
}

double dmin(double a, double b)
{
	if (a < b)
		return a;
	else
		return b;
}

double dmax(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
}

double random_double(void)
{
	return (double) rand() / (double) RAND_MAX;
}

int id_offset(int proc, int neltot, int nproc)
{
	return proc * (neltot / nproc) + min(proc, neltot % nproc);
}

int local_nel(int proc, int neltot, int nproc)
{
	return neltot / nproc + (proc < (neltot % nproc));
}

int id_of_elnum(int elnum, int proc, int neltot, int nproc)
{
	if (elnum < 0 || proc < 0)
		return elnum;
	return id_offset(proc, neltot, nproc) + elnum;
}

int elnum_of_id(int id, int neltot, int nproc)
{
	if (id < 0)
		return id;
	return id - id_offset(proc_of_id(id, neltot, nproc), neltot, nproc);
}

int proc_of_id(int id, int neltot, int nproc)
{
	for (int i = nproc; i--;)
		if (id >= id_offset(i, neltot, nproc))
			return i;
	return id;
}
