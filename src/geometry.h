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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "structs.h"

struct geometry {
	int id;
	int n;
	int nel;
	int neltot;
	int bc;
	struct trivec *grid;
	struct trimat *geom_mat;
	double ***jac_det;
	struct surface *surf_norm_x;
	struct surface *surf_norm_y;
	struct surface *surf_norm_z;
	struct trimat *jac;
	struct trimat *tinv_jac;
	int neigh_proc[6];
	int neigh_el[6];
	int neigh_rot[6];
	int edge_size[12];
	int neigh_edge_el[12];
	int neigh_edge_proc[12];
	int corner_size[8];
	int neigh_corn_el[8];
	int neigh_corn_proc[8];
};

struct geometry *new_geometry(struct basis *disc_basis, int nel, int neltot);
void free_geometry(struct geometry *geom);
struct geometry **new_geometries(struct basis *disc_basis, int nel, int neltot);
void free_geometries(struct geometry **geom);
void copy_geometry_connectivity(struct geometry *geom1, struct geometry *geom2);
void compute_geometry_transformations(struct geometry *geom,
		struct basis *basis);
void cartesian_grid(double *xi, double xa, double xb, double ya, double yb,
		double za, double zb, int n, struct trivec *grid);
void compute_geometry_transformation(struct trimat *jac, double ***jacd,
		int n, struct trimat *g);
void compute_surface_normal(struct trimat *jac, int n,
		struct surface *nx, struct surface *ny, struct surface *nz);
void compute_reference_domain_velocity(struct element **sel);
double element_volume(struct element *sel);
void get_neighbour_ids(int *neigh_el, int *neigh_proc, int *neigh_id,
		int neltot, int nproc, int m);
void get_neighbour_elnums_procs(int *neigh_el, int *neigh_proc, int *neigh_id,
		int neltot, int nproc, int m);
void set_boundary_conditions(struct geometry *geom);

#endif
