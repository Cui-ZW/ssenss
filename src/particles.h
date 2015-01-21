/*
 * Copyright 2015 Christopher Nilsen
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

#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdio.h>
#include <mpi.h>

struct particle {
	size_t id;
	int type;
	double x;
	double y;
	double z;
	double u;
	double v;
	double w;
	double t;
	int el;
	int old_el;
	int proc;
	int old_proc;
	int status;
	double tau;
};

struct part_params {
	size_t n;
	int ntypes;
	double tau_min;
	double tau_max;
};

struct plist {
	struct particle *list;
	size_t n;
	size_t ntot;
	size_t nbuff;
	int ntypes;
	double *tau;
	int nproc;
	MPI_Comm mpi_comm;
};

struct part_params *new_particle_parameters(void);
void free_particle_parameters(struct part_params *pparams);
struct plist *new_particle_list(struct part_params *pparams,
		struct element **sel);
void free_particle_list(struct plist *p);
size_t local_particle_count(struct element **sel, size_t ntot);
void get_particle_timescales(double tau_min, double tau_max, int n,
		double *tau);
void initialise_particles(struct plist *p, struct element **sel, size_t nel);
void initialise_particle_position(struct particle *p, struct element *sel);
void integrate_particle_equations(struct plist *p, struct element **sel);
double physical_particle_velocity(struct particle *p, struct element *sel,
		int comp);
double parametric_particle_velocity(struct particle *p, struct element *sel,
		int comp);
double physical_particle_coordinate(struct particle *p, struct element *sel,
		int comp);
int normal_first(const void *a, const void *b);
int increasing_proc(const void *a, const void *b);
double max_particle_velocity(struct plist *p, struct element **sel);
#endif
