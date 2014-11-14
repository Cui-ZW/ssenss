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
#include "utils.h"
#include "structs.h"
#include "basis.h"
#include "geometry.h"
#include "init.h"

static void get_displs(int *send, int *recv, int *displs, int n, int sproc,
		int rproc);
static void connect_buffer(double **surf, double *buff, int n, int displ);
static int surfind(int surf[6], int el_id, int first);
static void allocate_boundary_arrays(struct element **sel);
static void allocate_edge_arrays(struct element **sel);
static void allocate_corner_arrays(struct element **sel);
static void allocate_buffers(struct element **sel);
static void connect_local_elements(struct element **sel);
static void connect_local_edges(struct element **sel);
static void connect_local_corners(struct element **sel);

void connect_elements(struct element **sel)
{
	allocate_buffers(sel);
	allocate_boundary_arrays(sel);
	allocate_edge_arrays(sel);
	allocate_corner_arrays(sel);
	connect_local_elements(sel);
	connect_local_edges(sel);
	connect_local_corners(sel);
	connect_global_elements(sel);
	connect_global_edges(sel);
	connect_global_corners(sel);
}

static void allocate_buffers(struct element **sel)
{
	(*sel)->comm->sendcounts = new_int_array((*sel)->params->nproc);
	(*sel)->comm->senddispls = new_int_array((*sel)->params->nproc);

	size_t buffsize = compute_buffer_sizes(sel);

	(*sel)->sendbuff = new_1d_array(buffsize);
	(*sel)->recvbuff = new_1d_array(buffsize);

	for (int i = 1; i < (*sel)->nel; ++i) {
		sel[i]->sendbuff = (*sel)->sendbuff;
		sel[i]->recvbuff = (*sel)->recvbuff;
	}
}

static void connect_local_elements(struct element **sel)
{
	int si;
	int elid;
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 6; ++j) {
			si = -1;
			elid = sel[i]->geom->neigh_el[j];
			if (sel[i]->geom->neigh_proc[j]
					== sel[i]->params->rank) {
				while (si < 0 || sel[elid]->geom->neigh_proc[si]
						!= sel[i]->params->rank)
					si = surfind(sel[elid]->geom->neigh_el,
							i, si + 1);
				*sel[i]->neigh_val->ind[j] = *sel[elid]
						->bound_val->ind[si];
			}
		}
}

static void connect_local_edges(struct element **sel)
{
	int ei;
	int elid;
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 12; ++j) {
			ei = -1;
			elid = sel[i]->geom->neigh_edge_el[j];
			if (sel[i]->geom->neigh_edge_proc[j]
					== sel[i]->params->rank) {
				while (ei < 0 || sel[elid]->geom
						->neigh_edge_proc[ei] 
						!= sel[i]->params->rank)
					ei = edgeind(sel[elid]->geom
							->neigh_edge_el, i,
							ei + 1);
				*sel[i]->neigh_edge->ind[j] = *sel[elid]
						->bound_edge->ind[ei];
			}
		}
}

static void connect_local_corners(struct element **sel)
{
	int ci;
	int elid;
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 8; ++j) {
			ci = -1;
			elid = sel[i]->geom->neigh_corn_el[j];
			if (sel[i]->geom->neigh_corn_proc[j]
					== sel[i]->params->rank) {
				while (ci < 0 || sel[elid]->geom
						->neigh_corn_proc[ci]
						!= sel[i]->params->rank)
					ci = cornerind(sel[elid]->geom
							->neigh_corn_el, i,
							ci + 1);
				sel[i]->neigh_corn[j]
						= sel[elid]->bound_corn[ci];
			}
		}
}

void connect_global_elements(struct element **sel)
{
	int nsurf;
	const int n = (*sel)->basis->n;
	const int m = n * n;
	int *send, *recv, *displs, *surf;
	int k, si;
	int *senddispls = (*sel)->comm->senddispls;

	for (int p = 0; p < (*sel)->params->nproc; ++p) {
		if (p == (*sel)->params->rank)
			continue;
		nsurf = 0;
		for (int i = 0; i < (*sel)->nel; ++i)
			for (int j = 0; j < 6; ++j)
				if (sel[i]->geom->neigh_proc[j] == p)
					++nsurf;
		send = new_int_array(nsurf);
		recv = new_int_array(nsurf);
		surf = new_int_array(nsurf);
		displs = new_int_array(nsurf);

		k = 0;
		for (int i = 0; i < (*sel)->nel; ++i) {
			si = -1;
			while ((si = surfind(sel[i]->geom->neigh_proc, p, si + 1))
					>= 0) {
				send[k] = i;
				recv[k] = sel[i]->geom->neigh_el[si];
				surf[k++] = si;
			}
		}

		get_displs(send, recv, displs, nsurf, (*sel)->params->rank, p);

		for (int i = 0; i < nsurf; ++i) {
			connect_buffer(*sel[send[i]]->neigh_val->ind[surf[i]],
					(*sel)->recvbuff, n, senddispls[p]
					+ m * displs[i]);
			connect_buffer(*sel[send[i]]->bound_val->ind[surf[i]],
					(*sel)->sendbuff, n, senddispls[p]
					+ m * displs[i]);
		}
		free_int_array(send);
		free_int_array(recv);
		free_int_array(surf);
		free_int_array(displs);
	}

}

/* The following two functions are not general, assumes edge/corner_size = 1 */
void connect_global_edges(struct element **sel)
{
	int nsurf, nedge;
	const int n = (*sel)->basis->n;
	const int m = n * n;
	int k, ei;
	int *send, *recv, *edge, *displs;
	int *senddispls = (*sel)->comm->senddispls;

	for (int p = 0; p < (*sel)->params->nproc; ++p) {
		if (p == (*sel)->params->rank)
			continue;
		nedge = 0;
		nsurf = 0;
		for (int i = 0; i < (*sel)->nel; ++i) {
			for (int j = 0; j < 12; ++j)
				if (sel[i]->geom->neigh_edge_proc[j] == p)
					++nedge;
			for (int j = 0; j < 6; ++j)
				if (sel[i]->geom->neigh_proc[j] == p)
					++nsurf;
		}
		send = new_int_array(nedge);
		recv = new_int_array(nedge);
		edge = new_int_array(nedge);
		displs = new_int_array(nedge);

		k = 0;
		for (int i = 0; i < (*sel)->nel; ++i) {
			ei = -1;
			while ((ei = edgeind(sel[i]->geom->neigh_edge_proc, p,
					ei + 1)) >= 0) {
				send[k] = i;
				recv[k] = sel[i]->geom->neigh_edge_el[ei];
				edge[k++] = ei;
			}
		}

		get_displs(send, recv, displs, nedge, (*sel)->params->rank, p);

		for (int i = 0; i < nedge; ++i) {
			*sel[send[i]]->neigh_edge->ind[edge[i]] =
					(*sel)->recvbuff + senddispls[p]
					+ m * nsurf + n * displs[i];
			*sel[send[i]]->bound_edge->ind[edge[i]] =
					(*sel)->sendbuff + senddispls[p]
					+ m * nsurf + n * displs[i];
		}
		free_int_array(send);
		free_int_array(recv);
		free_int_array(edge);
		free_int_array(displs);
	}
}

void connect_global_corners(struct element **sel)
{
	int nsurf, nedge, ncorn;
	const int n = (*sel)->basis->n;
	const int m = n * n;
	int k, ci;
	int *send, *recv, *corner, *displs;
	int *senddispls = (*sel)->comm->senddispls;

	for (int p = 0; p < (*sel)->params->nproc; ++p) {
		if (p == (*sel)->params->rank)
			continue;
		nedge = 0;
		nsurf = 0;
		ncorn = 0;
		for (int i = 0; i < (*sel)->nel; ++i) {
			for (int j = 0; j < 8; ++j)
				if (sel[i]->geom->neigh_corn_proc[j] == p)
					++ncorn;
			for (int j = 0; j < 12; ++j)
				if (sel[i]->geom->neigh_edge_proc[j] == p)
					++nedge;
			for (int j = 0; j < 6; ++j)
				if (sel[i]->geom->neigh_proc[j] == p)
					++nsurf;
		}
		send = new_int_array(ncorn);
		recv = new_int_array(ncorn);
		corner = new_int_array(ncorn);
		displs = new_int_array(ncorn);

		k = 0;
		for (int i = 0; i < (*sel)->nel; ++i) {
			ci = -1;
			while ((ci = cornerind(sel[i]->geom->neigh_corn_proc, p,
					ci + 1)) >= 0) {
				send[k] = i;
				recv[k] = sel[i]->geom->neigh_corn_el[ci];
				corner[k++] = ci;
			}
		}

		get_displs(send, recv, displs, ncorn, (*sel)->params->rank, p);

		for (int i = 0; i < ncorn; ++i) {
			sel[send[i]]->neigh_corn[corner[i]] =
					(*sel)->recvbuff + senddispls[p]
					+ m * nsurf + n * nedge + displs[i];
			sel[send[i]]->bound_corn[corner[i]] =
					(*sel)->sendbuff + senddispls[p]
					+ m * nsurf + n * nedge + displs[i];
		}
		free_int_array(send);
		free_int_array(recv);
		free_int_array(corner);
		free_int_array(displs);
	}
}

static int surfind(int surf[6], int el_id, int first)
{
	for (int i = first; i < 6; ++i)
		if (surf[i] == el_id)
			return i;
	return -1;
}

int edgeind(int edge[12], int el_id, int first)
{
	for (int i = first; i < 12; ++i)
		if (edge[i] == el_id)
			return i;
	return -1;
}

int cornerind(int corner[8], int el_id, int first)
{
	for (int i = first; i < 8; ++i)
		if (corner[i] == el_id)
			return i;
	return -1;
}

/* This function assumes that only one surface/edge/corner is communicated
 * between the two elements on the two processors. Shift-size is quite
 * arbitrary.
 */
static void get_displs(int *send, int *recv, int *displs, int n, int sproc,
		int rproc)
{
	int sshift, rshift;
	long *unique = ec_malloc(n * sizeof(*unique));
	int k;

	if (sproc < rproc) {
		sshift = 16;
		rshift = 0;
	} else {
		sshift = 0;
		rshift = 16;
	}

	for (int j = 0; j < n; ++j)
		unique[j] = (send[j] << sshift) + (recv[j] << rshift);
	for (int j = 0; j < n; ++j) {
		k = 0;
		for (int i = 0; i < n; ++i)
		       if (unique[i] < unique[j])
			       ++k;
		displs[j] = k;
	}
	free(unique);
}

static void connect_buffer(double **surf, double *buff, int n, int displ)
{
	surf[0] = buff + displ;
	for (int i = 1; i < n; ++i)
		surf[i] = surf[i - 1] + n;
}

size_t compute_buffer_sizes(struct element **sel)
{
	const int n = (*sel)->basis->n;
	const int m = n * n;

	size_t buffsize = 0;
	int *sendcounts = (*sel)->comm->sendcounts;
	int *senddispls = (*sel)->comm->senddispls;
	for (int j = 0; j < (*sel)->nel; ++j) {
		/* Add surfaces */
		for (int i = 0; i < 6; ++i) {
			if (sel[j]->geom->neigh_proc[i] < 0)
				continue;
			if (sel[j]->geom->neigh_proc[i]
					!= sel[j]->params->rank) {
				sendcounts[sel[j]->geom->neigh_proc[i]] += m;
				buffsize += m;
			}
		}
		/* Add edges */
		for (int i = 0; i < 12; ++i) {
			if (sel[j]->geom->neigh_edge_proc[i] < 0)
				continue;
			if (sel[j]->geom->neigh_edge_proc[i] !=
					sel[j]->params->rank) {
				sendcounts[sel[j]->geom->neigh_edge_proc[i]] +=
						n * sel[j]->geom->edge_size[i];
				buffsize += n * sel[j]->geom->edge_size[i];
			}
		}
		/* Add corners */
		for (int i = 0; i < 8; ++i) {
			if (sel[j]->geom->neigh_corn_proc[i] < 0)
				continue;
			if (sel[j]->geom->neigh_corn_proc[i] !=
					sel[j]->params->rank) {
				sendcounts[sel[j]->geom->neigh_corn_proc[i]]
						+= sel[j]->geom->corner_size[i];
				buffsize += sel[j]->geom->corner_size[i];
			}
		}
	}
	for (int i = 1; i < (*sel)->params->nproc; ++i)
		senddispls[i] = senddispls[i - 1] + sendcounts[i - 1];
	return buffsize;
}

static void allocate_boundary_arrays(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 6; ++j)
			if (sel[i]->geom->neigh_proc[j]
					== sel[i]->params->rank) {
				*sel[i]->bound_val->ind[j] = new_2d_array(
						sel[i]->n, sel[i]->n);
			} else if (sel[i]->geom->neigh_proc[j] >= 0) {
				*sel[i]->bound_val->ind[j] = new_ptr_array(
						sel[i]->n, 0, NULL);
				*sel[i]->neigh_val->ind[j] = new_ptr_array(
						sel[i]->n, 0, NULL);
			}
}

static void allocate_edge_arrays(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 12; ++j)
			if (sel[i]->geom->neigh_edge_proc[j]
					== sel[i]->params->rank)
				*sel[i]->bound_edge->ind[j] = new_1d_array(
						sel[i]->n * sel[i]
						->geom->edge_size[j]);
}

static void allocate_corner_arrays(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 8; ++j)
			if (sel[i]->geom->neigh_corn_proc[j]
					== sel[i]->params->rank)
				sel[i]->bound_corn[j] = new_1d_array(sel[i]
						->geom->corner_size[j]);
}
