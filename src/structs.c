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

#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "basis.h"
#include "tripalg.h"
#include "geometry.h"
#include "test.h"
#include "precon.h"
#include "structs.h"

static struct old_element *new_old_element(struct element *sel, int *count);
static void free_old_element(struct old_element *old);

struct element *new_element(struct geometry *geom, struct geometry *conv_geom,
		struct basis *disc_basis, struct basis *conv_basis,
		struct parameters *params, struct communicator *comm)
{
	struct element *sel = ec_malloc(sizeof(*sel));
	int *count = ec_malloc(sizeof(*count));
	*count = 1;
	const int n = params->n;

	sel->n = n;
	sel->nsurf = n * n;
	sel->ntot = ((size_t) n) * ((size_t) n) * ((size_t) n);
	sel->params = params;
	sel->nel = geom->nel;
	sel->bc = geom->bc;

	sel->comm = comm;
	sel->geom = geom;
	sel->conv_geom = conv_geom;
	sel->basis = disc_basis;
	sel->conv_basis = conv_basis;
	sel->u = new_trivec(n);
	sel->p = new_3d_array(n, n, n);
	sel->conv = new_trivec(n);
	sel->rot = new_trivec(n);
	sel->f = new_trivec(n);
	sel->precon = new_preconditioner(n, 10);
	sel->ud = new_trivec(n);
	sel->ur = new_trivec(n);
	sel->bound_val = new_surface(0);
	sel->neigh_val = new_surface(0);
	sel->bound_edge = new_edges(0);
	sel->neigh_edge = new_edges(0);

	for (int i = 0; i < 3; ++i) {
		sel->bc_func[i] = zero_func;
		sel->bc_dfunc[i] = zero_func;
	}
	for (int i = 0; i < 8; ++i) {
		sel->bound_corn[i] = NULL;
		sel->neigh_corn[i] = NULL;
	}

	sel->old = new_old_element(sel, count);

	free(count);

	return sel;
}

void free_element(struct element *sel)
{
	free_trivec(sel->u);
	free_trivec(sel->ud);
	free_trivec(sel->ur);
	free_trivec(sel->conv);
	free_trivec(sel->rot);
	free_trivec(sel->f);
	free_3d_array(sel->p);
	free_preconditioner(sel->precon, sel->geom->id);
	free(sel->bound_val);
	free(sel->neigh_val);
	free(sel->bound_edge);
	free(sel->neigh_edge);

	if (sel->old)
		free_old_element(sel->old);
	free(sel);
}

struct element **new_elements(struct geometry **geom,
		struct geometry **conv_geom, struct basis *disc_basis,
		struct basis *conv_basis, struct parameters *params,
		struct communicator *comm)
{
	struct element **sel = new_ptr_array((*geom)->nel, 0, NULL);
	for (int i = 0; i < (*geom)->nel; ++i)
		sel[i] = new_element(geom[i], conv_geom[i], disc_basis,
				conv_basis, params, comm);
	return sel;
}

void free_elements(struct element **sel)
{
	const int nel = (*sel)->nel;
	for (int i = 0; i < nel; ++i)
		free_element(sel[i]);
	free_ptr_array(sel, 0, NULL);
}

static struct old_element *new_old_element(struct element *sel, int *count)
{
	if (*count >= 3)
		return NULL;
	struct old_element *old = ec_malloc(sizeof(*old));
	(*count)++;

	old->u = new_trivec(sel->n);
	old->ur = new_trivec(sel->n);
	old->p = new_3d_array(sel->n, sel->n, sel->n);
	old->conv = new_trivec(sel->n);
	old->rot = new_trivec(sel->n);
	old->old = new_old_element(sel, count);

	return old;
}

static void free_old_element(struct old_element *old)
{
	free_trivec(old->u);
	free_trivec(old->ur);
	free_trivec(old->conv);
	free_trivec(old->rot);
	free_3d_array(old->p);

	if (old->old)
		free_old_element(old->old);
	free(old);
}

struct surface *new_surface(int n)
{
	struct surface *surf = ec_malloc(sizeof(*surf));

	surf->ind[0] = &surf->w;
	surf->ind[1] = &surf->e;
	surf->ind[2] = &surf->s;
	surf->ind[3] = &surf->n;
	surf->ind[4] = &surf->b;
	surf->ind[5] = &surf->t;

	if (n == 0)
		for (int i = 0; i < 6; ++i)
			*surf->ind[i] = NULL;
	else
		for (int i = 0; i < 6; ++i)
			*surf->ind[i] = new_2d_array(n, n);

	return surf;
}

void free_surface(struct surface *surf)
{
	for (int i = 0; i < 6; ++i)
		free_2d_array(*surf->ind[i]);
	free(surf);
}

struct edges *new_edges(int n)
{
	struct edges *edg = ec_malloc(sizeof(*edg));

	edg->ind[0] = &edg->sw;
	edg->ind[1] = &edg->se;
	edg->ind[2] = &edg->nw;
	edg->ind[3] = &edg->ne;
	edg->ind[4] = &edg->bw;
	edg->ind[5] = &edg->be;
	edg->ind[6] = &edg->tw;
	edg->ind[7] = &edg->te;
	edg->ind[8] = &edg->bs;
	edg->ind[9] = &edg->bn;
	edg->ind[10] = &edg->ts;
	edg->ind[11] = &edg->tn;

	if (n == 0)
		for (int i = 0; i < 12; ++i)
			*edg->ind[i] = NULL;
	else
		for (int i = 0; i < 12; ++i)
			*edg->ind[i] = new_1d_array(n);

	return edg;
}

void free_edges(struct edges *edg)
{
	for (int i = 0; i < 12; ++i)
		free_1d_array(*edg->ind[i]);
	free(edg);
}

struct parameters *new_parameters(int rank, int nproc)
{
	struct parameters *params = ec_malloc(sizeof(*params));

	params->rank = rank;
	params->nproc = nproc;
	params->dt = 1.0;
	params->visc = 1.0;
	params->dens = 1.0;
	params->t = 0.0;
	params->tmax = 1.0;
	params->gx = 0.0;
	params->gy = 0.0;
	params->gz = 0.0;
	params->n = 2;
	params->nconv = 2;
	params->nel = 0;
	params->neltot = 0;
	params->test = 0;
	params->svvn = 2;
	params->svva = 0.0;
	params->wtime[0] = 0.0;
	params->wtime[1] = 0.0;
	params->wtime[2] = 0.0;
	params->write_every = 1000;
	strcpy(params->gridfile, "dat/grid.h5");
	strcpy(params->initfile, "dat/init.h5");
	strcpy(params->datdir, "dat");

	return params;
}

void free_parameters(struct parameters *params)
{
	free(params);
}

struct communicator *new_communicator(int rank, int nproc, MPI_Comm mpi_comm,
		MPI_Info mpi_info)
{
	struct communicator *comm = ec_malloc(sizeof(*comm));
	comm->rank = rank;
	comm->nproc = nproc;
	comm->mpi_comm = mpi_comm;
	comm->mpi_info = mpi_info;
	return comm;
}

void free_communicator(struct communicator *comm)
{
	free(comm);
}
