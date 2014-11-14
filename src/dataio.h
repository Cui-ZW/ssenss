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

#ifndef DATAIO_H
#define DATAIO_H

#include <mpi.h>
#include <hdf5.h>
#include "structs.h"
#include "geometry.h"

enum { READ_GRID = 01, READ_VELOCITY = 02, READ_PARTICLES = 04,
		WRITE_GRID = 010, WRITE_VELOCITY = 020,	WRITE_PARTICLES = 040 };

void write_log(struct element **sel);
void write_to_file(struct element **sel, const char *filename, int write_data);
void read_velocity(struct element **sel, const char *filename);
void get_grid_size(struct parameters *params, struct communicator *comm,
		const char *filename);
void read_grid(struct geometry **geom, struct communicator *comm,
		struct basis *basis, const char *filename);

#endif
