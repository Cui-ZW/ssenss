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

#ifndef PRECON_H
#define PRECON_H

#include "structs.h"

#define NORMAL 1
#define TRANSPOSE 2

struct preconditioner {
	int nlevels;
	int level;
	int *nmg;
	struct basis **basis;
	struct geometry **geom;
	double ****helmholtz_diag;
	double ****poisson_diag;
	double ****stiffsum_scale;
	double ****dirichlet_mask;
	double ***inverse_mass;
	int **sendcounts;
	int **senddispls;
	double *poisson_eigen_max;
	double *helmholtz_eigen_max;
};

struct preconditioner *new_preconditioner(int n, int nlevels);
void free_preconditioner(struct preconditioner *precon, int id);
void initialise_preconditioners(struct element **sel);
void multigrid_preconditioner(struct element **sel, double ****x, double ****px,
		int var, double eps, void (*opeval) ());
void jacobi_preconditioner(struct element **sel, double ****x, double ****px,
		int var, double eps, void (*opeval) ());
void fdm_preconditioner(struct element **sel, double ****x, double ****px,
		int var, double eps, void (*opeval) ());
#endif
