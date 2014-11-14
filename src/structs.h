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

#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdio.h>
#include <mpi.h>

struct surface {
	double **w;
	double **e;
	double **s;
	double **n;
	double **b;
	double **t;
	double ***ind[6];
};

struct edges {
	double *sw;
	double *se;
	double *nw;
	double *ne;
	double *bw;
	double *be;
	double *tw;
	double *te;
	double *bs;
	double *bn;
	double *ts;
	double *tn;
	double **ind[12];
};

struct parameters {
	int n;
	int nconv;
	int nel;
	int neltot;
	int svvn;
	int rank;
	int nproc;
	int test;
	int write_every;
	double dt;
	double visc;
	double dens;
	double gx;
	double gy;
	double gz;
	double t;
	double tmax;
	double svva;
	double wtime[3];
	char datdir[100];
	char gridfile[100];
	char initfile[100];
};

struct communicator {
	MPI_Comm mpi_comm;
	MPI_Info mpi_info;
	int rank;
	int nproc;
	int *sendcounts;
	int *senddispls;
	int buffsize;
};

struct element {
	int n;
	int nsurf;
	size_t ntot;
	int nel;
	double *sendbuff;
	double *recvbuff;
	struct parameters *params;
	struct communicator *comm;
	struct trivec *u;
	struct trivec *ud;
	struct trivec *ur;
	double ***p;
	struct geometry *geom;
	struct geometry *conv_geom;
	struct basis *basis;
	struct basis *conv_basis;
	struct trivec *conv;
	struct trivec *rot;
	struct trivec *f;
	struct preconditioner *precon;
	int bc;
	double (*bc_func[3]) ();
	double (*bc_dfunc[3]) ();
	struct surface *bound_val;
	struct surface *neigh_val;
	struct edges *bound_edge;
	struct edges *neigh_edge;
	double *bound_corn[8];
	double *neigh_corn[8];
	struct old_element *old;
};

struct old_element {
	struct trivec *u;
	struct trivec *ur;
	struct trivec *conv;
	struct trivec *rot;
	double ***p;
	struct old_element *old;
};

struct element *new_element(struct geometry *geom, struct geometry *conv_geom,
		struct basis *disc_basis, struct basis *conv_basis,
		struct parameters *params, struct communicator *comm);
void free_element(struct element *sel);
struct element **new_elements(struct geometry **geom,
		struct geometry **conv_geom, struct basis *disc_basis,
		struct basis *conv_basis, struct parameters *params,
		struct communicator *comm);
void free_elements(struct element **sel);
struct surface *new_surface(int n);
void free_surface(struct surface *surf);
struct edges *new_edges(int n);
void free_edges(struct edges *edg);
struct parameters *new_parameters(int rank, int nproc);
void free_parameters(struct parameters *params);
struct communicator *new_communicator(int rank, int nproc, MPI_Comm mpi_comm,
		MPI_Info mpi_info);
void free_communicator(struct communicator *comm);

#endif
