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

#include <mpi.h>
#include "structs.h"
#include "basis.h"
#include "linalg.h"
#include "geometry.h"
#include "utils.h"
#include "precon.h"
#include "global.h"

static void populate_surface(double ***array, struct surface *surf, int n);
static void add_surface(double ***array, struct surface *surf, int n);
static void populate_edges(double ***array, struct edges *edg, int n);
static void add_edges(double ***array, struct edges *edg, int n);
static void populate_corners(double ***array, double **corners, int n);
static void add_corners(double ***array, double **corners, int n);
static void exchange_boundaries(struct element **sel);
static void sum_edges(struct element **sel);
static void sum_corners(struct element **sel);

static void populate_surface(double ***array, struct surface *surf, int n)
{
	if (surf->w)
		copy_inc_array(**array, *surf->w, n * n, n, 1);
	if (surf->e)
		copy_inc_array(**array + n - 1, *surf->e, n * n, n, 1);
	if (surf->s)
		for (int j = 0; j < n; ++j)
			copy_inc_array(array[j][0], *surf->s + n * j, n, 1, 1);
	if (surf->n)
		for (int j = 0; j < n; ++j)
			copy_inc_array(array[j][n - 1],
					*surf->n + n * j, n, 1, 1);
	if (surf->b)
		copy_inc_array(*array[0], *surf->b, n * n, 1, 1);
	if (surf->t)
		copy_inc_array(*array[n - 1], *surf->t, n * n, 1, 1);
}

static void add_surface(double ***array, struct surface *surf, int n)
{
	if (surf->w)
		add_inc_array(*surf->w, **array, n * n, 1.0, 1, n);
	if (surf->e)
		add_inc_array(*surf->e, **array + n - 1, n * n, 1.0, 1, n);
	if (surf->s)
		for (int j = 0; j < n; ++j)
			add_inc_array(*surf->s + n * j,
					array[j][0], n, 1.0, 1, 1);
	if (surf->n)
		for (int j = 0; j < n; ++j)
			add_inc_array(*surf->n + n * j,
					array[j][n - 1], n, 1.0, 1, 1);
	if (surf->b)
		add_inc_array(*surf->b, *array[0], n * n, 1.0, 1, 1);
	if (surf->t)
		add_inc_array(*surf->t, *array[n - 1], n * n, 1.0, 1, 1);
}

static void populate_edges(double ***array, struct edges *edg, int n)
{
	if (edg->sw)
		for (int i = 0; i < n; ++i)
			edg->sw[i] = array[i][0][0];
	if (edg->se)
		for (int i = 0; i < n; ++i)
			edg->se[i] = array[i][0][n - 1];
	if (edg->nw)
		for (int i = 0; i < n; ++i)
			edg->nw[i] = array[i][n - 1][0];
	if (edg->ne)
		for (int i = 0; i < n; ++i)
			edg->ne[i] = array[i][n - 1][n - 1];
	if (edg->bw)
		for (int i = 0; i < n; ++i)
			edg->bw[i] = array[0][i][0];
	if (edg->be)
		for (int i = 0; i < n; ++i)
			edg->be[i] = array[0][i][n - 1];
	if (edg->tw)
		for (int i = 0; i < n; ++i)
			edg->tw[i] = array[n - 1][i][0];
	if (edg->te)
		for (int i = 0; i < n; ++i)
			edg->te[i] = array[n - 1][i][n - 1];
	if (edg->bs)
		for (int i = 0; i < n; ++i)
			edg->bs[i] = array[0][0][i];
	if (edg->bn)
		for (int i = 0; i < n; ++i)
			edg->bn[i] = array[0][n - 1][i];
	if (edg->ts)
		for (int i = 0; i < n; ++i)
			edg->ts[i] = array[n - 1][0][i];
	if (edg->tn)
		for (int i = 0; i < n; ++i)
			edg->tn[i] = array[n - 1][n - 1][i];
}

static void add_edges(double ***array, struct edges *edg, int n)
{
	if (edg->sw)
		for (int i = 0; i < n; ++i)
			array[i][0][0] += edg->sw[i];
	if (edg->se)
		for (int i = 0; i < n; ++i)
			array[i][0][n - 1] += edg->se[i];
	if (edg->nw)
		for (int i = 0; i < n; ++i)
			array[i][n - 1][0] += edg->nw[i];
	if (edg->ne)
		for (int i = 0; i < n; ++i)
			array[i][n - 1][n - 1] += edg->ne[i];
	if (edg->bw)
		for (int i = 0; i < n; ++i)
			array[0][i][0] += edg->bw[i];
	if (edg->be)
		for (int i = 0; i < n; ++i)
			array[0][i][n - 1] += edg->be[i];
	if (edg->tw)
		for (int i = 0; i < n; ++i)
			array[n - 1][i][0] += edg->tw[i];
	if (edg->te)
		for (int i = 0; i < n; ++i)
			array[n - 1][i][n - 1] += edg->te[i];
	if (edg->bs)
		for (int i = 0; i < n; ++i)
			array[0][0][i] += edg->bs[i];
	if (edg->bn)
		for (int i = 0; i < n; ++i)
			array[0][n - 1][i] += edg->bn[i];
	if (edg->ts)
		for (int i = 0; i < n; ++i)
			array[n - 1][0][i] += edg->ts[i];
	if (edg->tn)
		for (int i = 0; i < n; ++i)
			array[n - 1][n - 1][i] += edg->tn[i];
}

static void populate_corners(double ***array, double **corners, int n)
{
	for (int i = 0; i < 8; ++i) {
		const int xid = (i % 2) * (n - 1);
		const int yid = ((i / 2) % 2) * (n - 1);
		const int zid = (i / 4) * (n - 1);
		if (corners[i])
			*corners[i] = array[zid][yid][xid];
	}
}

static void add_corners(double ***array, double **corners, int n)
{
	for (int i = 0; i < 8; ++i) {
		const int xid = (i % 2) * (n - 1);
		const int yid = ((i / 2) % 2) * (n - 1);
		const int zid = (i / 4) * (n - 1);
		if (corners[i])
			array[zid][yid][xid] += *corners[i];
	}
}

static void exchange_boundaries(struct element **sel)
{
	int *sendcounts = (*sel)->comm->sendcounts;
	int *senddispls = (*sel)->comm->senddispls;
	int error = MPI_Alltoallv((*sel)->sendbuff, sendcounts, senddispls,
			MPI_DOUBLE, (*sel)->recvbuff, sendcounts, senddispls,
			MPI_DOUBLE, (*sel)->comm->mpi_comm);
	if (error)
		fatal_error("in exchange_boundaries for MPI_Alltoallv");
}

static void sum_edges(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 12; ++j)
			for (int k = 1; k < sel[i]->geom->edge_size[j]; ++k)
				add_array(*sel[i]->neigh_edge->ind[j] + k,
						*sel[i]->neigh_edge->ind[j],
						sel[i]->nsurf, ADD_NEW);
}

static void sum_corners(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 8; ++j)
			for (int k = 1; k < sel[i]->geom->corner_size[j]; ++k)
				*sel[i]->neigh_corn[j] +=
						*(sel[i]->neigh_corn[j] + k);
}

double global_sum(struct element **sel, double u)
{
	double sum;
	int error = MPI_Allreduce(&u, &sum, 1, MPI_DOUBLE, MPI_SUM,
			(*sel)->comm->mpi_comm);
	if (error)
		fatal_error("in global_sum for MPI_Allreduce");

	return sum;
}

int global_int_sum(MPI_Comm mpi_comm, int n)
{
	int ntot;
	int error = MPI_Allreduce(&n, &ntot, 1, MPI_INT, MPI_SUM, mpi_comm);
	if (error)
		fatal_error("in global_int_sum for MPI_Allreduce");

	return ntot;
}

double global_max(struct element **sel, double u)
{
	double max;
	int error = MPI_Allreduce(&u, &max, 1, MPI_DOUBLE, MPI_MAX,
			(*sel)->comm->mpi_comm);
	if (error)
		fatal_error("in global_sum for MPI_Allreduce");

	return max;
}

void direct_stiffness_summation(struct element **sel, double ****u)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		populate_surface(u[i], sel[i]->bound_val, sel[i]->basis->n);
		populate_edges(u[i], sel[i]->bound_edge, sel[i]->basis->n);
		populate_corners(u[i], sel[i]->bound_corn, sel[i]->basis->n);
	}
	exchange_boundaries(sel);
/*	sum_edges(sel);
	sum_corners(sel);
*/
	for (int i = 0; i < (*sel)->nel; ++i) {
		add_surface(u[i], sel[i]->neigh_val, sel[i]->basis->n);
		add_edges(u[i], sel[i]->neigh_edge, sel[i]->basis->n);
		add_corners(u[i], sel[i]->neigh_corn, sel[i]->basis->n);
	}
}

void direct_stiffness_averaging(struct element **sel, double ****u)
{
	const int n = (*sel)->basis->n;
	double ****tmp = new_4d_array(n, n, n, (*sel)->nel);
	copy_array(***u, ***tmp, n * n * n * (*sel)->nel);
	direct_stiffness_summation(sel, tmp);
	multiply(***tmp, **(*sel)->precon->stiffsum_scale[(*sel)->precon
			->level], ***u, n * n * n * (*sel)->nel, 1.0, 0.0);
	free_4d_array(tmp);
}

double global_dot_product(struct element **sel, double ****x, double ****y)
{
	const int n = (*sel)->basis->n;
	double ****tmp = new_4d_array(n, n, n, (*sel)->nel);

	multiply(***x, **(*sel)->precon->stiffsum_scale[(*sel)->precon->level],
			***tmp, n * n * n * (*sel)->nel, 1.0, 0.0);
	double prod = dot_product(***tmp, ***y, n * n * n * (*sel)->nel);
	double globprod = global_sum(sel, prod);
	free_4d_array(tmp);

	return globprod;
}
