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
#include <mpi.h>
#include "utils.h"
#include "structs.h"
#include "geometry.h"
#include "final.h"

static void free_buffers(struct element **sel);
static void free_boundary_arrays(struct element **sel);
static void free_edge_arrays(struct element **sel);
static void free_corner_arrays(struct element **sel);

void disconnect_elements(struct element **sel)
{
	free_buffers(sel);
	free_boundary_arrays(sel);
	free_edge_arrays(sel);
	free_corner_arrays(sel);
}

static void free_buffers(struct element **sel)
{
	free_int_array((*sel)->comm->sendcounts);
	free_int_array((*sel)->comm->senddispls);
	free_1d_array((*sel)->sendbuff);
	free_1d_array((*sel)->recvbuff);
}

static void free_boundary_arrays(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 6; ++j)
			if (sel[i]->geom->neigh_proc[j]
					== sel[i]->params->rank) {
				free_2d_array(*sel[i]->bound_val->ind[j]);
			} else if (sel[i]->geom->neigh_proc[j] >= 0) {
				free_ptr_array(*sel[i]->bound_val->ind[j], 0,
						NULL);
				free_ptr_array(*sel[i]->neigh_val->ind[j], 0,
						NULL);
			}
}

static void free_edge_arrays(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 12; ++j)
			if (sel[i]->geom->neigh_edge_proc[j]
					== sel[i]->params->rank)
				free_1d_array(*sel[i]->bound_edge->ind[j]);
}

static void free_corner_arrays(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i)
		for (int j = 0; j < 8; ++j)
			if (sel[i]->geom->neigh_corn_proc[j]
					== sel[i]->params->rank)
				free_1d_array(sel[i]->bound_corn[j]);
}

void stop_mpi(void)
{
	MPI_Finalize();
}
