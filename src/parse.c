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
#include <string.h>
#include "utils.h"
#include "structs.h"
#include "poisson.h"
#include "parse.h"

static double read_double(const char *val);
static int read_int(const char *val);

void parse_command_arguments(int argc, char **argv, struct parameters *params)
{
	for (int i = 1; i < argc; i += 1) {
		if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--n"))
			params->n = read_int(argv[++i]);
		if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--dens"))
			params->dens = read_double(argv[++i]);
		if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--visc"))
			params->visc = read_double(argv[++i]);
		if (!strcmp(argv[i], "-dt") || !strcmp(argv[i], "--dt"))
			params->dt = read_double(argv[++i]);
		if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--tmax"))
			params->tmax = read_double(argv[++i]);
		if (!strcmp(argv[i], "-rt") || !strcmp(argv[i], "--run-test"))
			params->test = 1;
	}
}

void set_parameter_value(struct parameters *params, char *name, char *val)
{
	if (!strcmp(name, "dt"))
		params->dt = read_double(val);
	if (!strcmp(name, "tmax"))
		params->tmax = read_double(val);
	if (!strcmp(name, "visc"))
		params->visc = read_double(val);
	if (!strcmp(name, "dens"))
		params->dens = read_double(val);
	if (!strcmp(name, "gx"))
		params->gx = read_double(val);
	if (!strcmp(name, "gy"))
		params->gy = read_double(val);
	if (!strcmp(name, "gz"))
		params->gz = read_double(val);
	if (!strcmp(name, "n"))
		params->n = read_int(val);
	if (!strcmp(name, "write_every"))
		params->write_every = read_int(val);
	if (!strcmp(name, "datdir"))
		strcpy(params->datdir, val);
	if (!strcmp(name, "gridfile"))
		strcpy(params->gridfile, val);
	if (!strcmp(name, "initfile"))
		strcpy(params->initfile, val);
}

static double read_double(const char *val)
{
	double tmp;
	sscanf(val, "%lf", &tmp);
	return tmp;
}

static int read_int(const char *val)
{
	int tmp;
	sscanf(val, "%d", &tmp);
	return tmp;
}
