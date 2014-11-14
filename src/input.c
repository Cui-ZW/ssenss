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
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "input.h"
#include "utils.h"
#include "poisson.h"
#include "parse.h"
#include "init.h"

#define N 500

static int read_parameters(struct parameters *params);
static void distribute_parameters(struct parameters *params);
static int read_block(FILE *fp, char *name, char *data, int n1, int n2);
static char *read_field(char *data, char *name, char *val);

void get_parameters(struct parameters *params)
{
	int file_error;
	if (params->rank == 0)
		file_error = read_parameters(params);
	int mpi_error = MPI_Bcast(&file_error, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (mpi_error)
		fatal_error("in get_parameters for MPI_Bcast");
	if (file_error)
		fatal_error("in read_parameters: unable to read setup file");
	distribute_parameters(params);
}

static int read_parameters(struct parameters *params)
{
	char *data, *name, *val;
	char *buff = malloc(3 * N);
	data = buff;
	name = data + N;
	val = name + N;

	FILE *fp = fopen("setup", "r");
	
	if (fp == NULL)
		return 1;

	while (read_block(fp, name, data, N, N)) {
		if (!strcmp(name, "parameters")) {
			while ((data = read_field(data, name, val)))
				set_parameter_value(params, name, val);
			break;
		}
		data = buff;
	}
	free(buff);
	fclose(fp);
	return 0;
}

static void distribute_parameters(struct parameters *params)
{
	const int rank = params->rank;
	int error = MPI_Bcast(params, sizeof(*params), MPI_BYTE, 0,
			MPI_COMM_WORLD);
	if (error)
		fatal_error("in distribute_parameters for MPI_Bcast");
	params->rank = rank;
}

static int read_block(FILE *fp, char *name, char *data, int n1, int n2)
{
	int c;
	int i = 0, j = 0;
	int read_name = 0, read_data = 0;

	while ((c = getc(fp)) != EOF && i < n1 && j < n2) {
		if (!isspace(c)) {
			if (read_name)
				name[i++] = c;
			else if (read_data)
				data[j++] = c;
		}
		if (c == '#') {
			read_name = 1;
		} else if (c == '{') {
			read_name = 0;
			read_data = 1;
			name[i - 1] = '\0';
		} else if (c == '}') {
			data[j - 1] = '\0';
			return j + i - 2;
		}
	}
	return 0;
}

static char *read_field(char *buff, char *name, char *val)
{
	int read_name = 1, read_val = 0;

	--buff;
	while (*++buff != '\0') {
		if (read_val)
			*val++ = *buff;
		else if (read_name && *buff != ';')
			*name++ = *buff;
		if (*buff == ':' && read_name) {
			read_name = 0;
			read_val = 1;
			*(name - 1) = '\0';
		} else if (*buff == ';' && read_val) {
			*(val - 1) = '\0';
			return buff;
		}
	}
	return NULL;
}
