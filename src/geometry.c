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
#include <stdlib.h>
#include "utils.h"
#include "linalg.h"
#include "tripalg.h"
#include "basis.h"
#include "time.h"
#include "deriv.h"
#include "poisson.h"
#include "geometry.h"

struct geometry *new_geometry(struct basis *disc_basis, int nel, int neltot)
{
	struct geometry *geom = ec_malloc(sizeof(*geom));

	const int n = disc_basis->n;
	geom->n = n;
	geom->nel = nel;
	geom->neltot = neltot;
	geom->geom_mat = new_sym_trimat(n);
	geom->jac_det = new_3d_array(n, n, n);
	geom->grid = new_trivec(n);
	geom->bc = 0;

	geom->tinv_jac = new_trimat(n);
	geom->jac = new_trimat(n);

	geom->surf_norm_x = new_surface(n);
	geom->surf_norm_y = new_surface(n);
	geom->surf_norm_z = new_surface(n);

	return geom;
}

void free_geometry(struct geometry *geom)
{
	free_trivec(geom->grid);
	free_sym_trimat(geom->geom_mat);
	free_trimat(geom->tinv_jac);
	free_trimat(geom->jac);
	free_3d_array(geom->jac_det);
	free_surface(geom->surf_norm_x);
	free_surface(geom->surf_norm_y);
	free_surface(geom->surf_norm_z);
	free(geom);
}

struct geometry **new_geometries(struct basis *disc_basis, int nel, int neltot)
{
	struct geometry **geom = new_ptr_array(nel, 0, NULL);
	for (int i = 0; i < nel; ++i)
		geom[i] = new_geometry(disc_basis, nel, neltot);
	return geom;
}

void free_geometries(struct geometry **geom)
{
	const int nel = (*geom)->nel;
	for (int i = 0; i < nel; ++i)
		free_geometry(geom[i]);
	free_ptr_array(geom, 0, NULL);
}

void copy_geometry_connectivity(struct geometry *geom1, struct geometry *geom2)
{
	geom2->id = geom1->id;
	for (int i = 0; i < 6; ++i) {
		geom2->neigh_proc[i] = geom1->neigh_proc[i];
		geom2->neigh_el[i] = geom1->neigh_el[i];
		geom2->neigh_rot[i] = geom1->neigh_rot[i];
	}
	for (int i = 0; i < 12; ++i) {
		geom2->neigh_edge_proc[i] = geom1->neigh_edge_proc[i];
		geom2->neigh_edge_el[i] = geom1->neigh_edge_el[i];
		geom2->edge_size[i] = geom1->edge_size[i];
	}
	for (int i = 0; i < 8; ++i) {
		geom2->neigh_corn_proc[i] = geom1->neigh_corn_proc[i];
		geom2->neigh_corn_el[i] = geom1->neigh_corn_el[i];
		geom2->corner_size[i] = geom1->corner_size[i];
	}
}

void compute_geometry_transformations(struct geometry *geom,
		struct basis *basis)
{
	compute_jacobian(basis, geom->grid, geom->jac);
	compute_trimat_determinant(geom->jac, basis->n, geom->jac_det);
	invert_trimat(geom->jac, geom->jac_det, basis->n, geom->tinv_jac);
	compute_geometry_transformation(geom->tinv_jac, geom->jac_det, basis->n,
			geom->geom_mat);
	compute_surface_normal(geom->jac, basis->n, geom->surf_norm_x,
			geom->surf_norm_y, geom->surf_norm_z);
	transpose_trimat(geom->tinv_jac);
}

void cartesian_grid(double *xi, double xa, double xb, double ya, double yb,
		double za, double zb, int n, struct trivec *grid)
{
	for (int k = 0; k < n; ++k)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i) {
				grid->x[k][j][i] = xa + 0.5 * (1.0 + xi[i])
						* (xb - xa);
				grid->y[k][j][i] = ya + 0.5 * (1.0 + xi[j])
						* (yb - ya);
				grid->z[k][j][i] = za + 0.5 * (1.0 + xi[k])
						* (zb - za);
			}
}

void compute_geometry_transformation(struct trimat *jaci, double ***jacd, int n,
		struct trimat *g)
{
	for (int i = 0; i < 3; ++i)
		trimat_trivec_product(jaci, jaci->row[i], n, g->row[i]);
}

void compute_surface_normal(struct trimat *jac, int n, struct surface *nx,
		struct surface *ny, struct surface *nz)
{
	const int m = n - 1;

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i) {
			nx->b[j][i] = -(jac->y.x[0][j][i] * jac->z.y[0][j][i] -
					jac->z.x[0][j][i] * jac->y.y[0][j][i]);
			ny->b[j][i] = (jac->x.x[0][j][i] * jac->z.y[0][j][i] -
					jac->z.x[0][j][i] * jac->x.y[0][j][i]);
			nz->b[j][i] = -(jac->x.x[0][j][i] * jac->y.y[0][j][i] -
					jac->y.x[0][j][i] * jac->x.y[0][j][i]);

			nx->t[j][i] = (jac->y.x[m][j][i] * jac->z.y[m][j][i] -
					jac->z.x[m][j][i] * jac->y.y[m][j][i]);
			ny->t[j][i] = -(jac->x.x[m][j][i] * jac->z.y[m][j][i] -
					jac->z.x[m][j][i] * jac->x.y[m][j][i]);
			nz->t[j][i] = (jac->x.x[m][j][i] * jac->y.y[m][j][i] -
					jac->y.x[m][j][i] * jac->x.y[m][j][i]);

			nx->s[j][i] = (jac->y.x[j][0][i] * jac->z.z[j][0][i] -
					jac->z.x[j][0][i] * jac->y.z[j][0][i]);
			ny->s[j][i] = -(jac->x.x[j][0][i] * jac->z.z[j][0][i] -
					jac->z.x[j][0][i] * jac->x.z[j][0][i]);
			nz->s[j][i] = (jac->x.x[j][0][i] * jac->y.z[j][0][i] -
					jac->y.x[j][0][i] * jac->x.z[j][0][i]);

			nx->n[j][i] = -(jac->y.x[j][m][i] * jac->z.z[j][m][i] -
					jac->z.x[j][m][i] * jac->y.z[j][m][i]);
			ny->n[j][i] = (jac->x.x[j][m][i] * jac->z.z[j][m][i] -
					jac->z.x[j][m][i] * jac->x.z[j][m][i]);
			nz->n[j][i] = -(jac->x.x[j][m][i] * jac->y.z[j][m][i] -
					jac->y.x[j][m][i] * jac->x.z[j][m][i]);

			nx->w[j][i] = -(jac->y.y[j][i][0] * jac->z.z[j][i][0] -
					jac->z.y[j][i][0] * jac->y.z[j][i][0]);
			ny->w[j][i] = (jac->x.y[j][i][0] * jac->z.z[j][i][0] -
					jac->z.y[j][i][0] * jac->x.z[j][i][0]);
			nz->w[j][i] = -(jac->x.y[j][i][0] * jac->y.z[j][i][0] -
					jac->y.y[j][i][0] * jac->x.z[j][i][0]);

			nx->e[j][i] = (jac->y.y[j][i][m] * jac->z.z[j][i][m] -
					jac->z.y[j][i][m] * jac->y.z[j][i][m]);
			ny->e[j][i] = -(jac->x.y[j][i][m] * jac->z.z[j][i][m] -
					jac->z.y[j][i][m] * jac->x.z[j][i][m]);
			nz->e[j][i] = (jac->x.y[j][i][m] * jac->y.z[j][i][m] -
					jac->y.y[j][i][m] * jac->x.z[j][i][m]);
		}
}

void compute_reference_domain_velocity(struct element **sel)
{
	for (int i = 0; i < (*sel)->nel; ++i) {
		time_shift(sel[i]->ur, sel[i]->old->ur, sel[i]->old->old->ur,
				sel[i]->ntot);
		transpose_trimat(sel[i]->geom->tinv_jac);
		trimat_trivec_product(sel[i]->geom->tinv_jac, sel[i]->u,
				sel[i]->n, sel[i]->ur);
		transpose_trimat(sel[i]->geom->tinv_jac);
	}
}

double element_volume(struct element *sel)
{
	const int n = sel->basis->n;
	double ***a = new_3d_array(n, n, n);
	double ***b = new_3d_array(n, n, n);
	ones(**a, n * n * n);
	test_function_inner_product(sel, a, b);
	double vol = dot_product(**a, **b, n * n * n);
	free_3d_array(a);
	free_3d_array(b);
	return vol;
}

void get_neighbour_ids(int *neigh_el, int *neigh_proc, int *neigh_id,
		int neltot, int nproc, int m)
{
	for (int i = 0; i < m; ++i)
		neigh_id[i] = id_of_elnum(neigh_el[i], neigh_proc[i], neltot,
				nproc);
}

void get_neighbour_elnums_procs(int *neigh_el, int *neigh_proc, int *neigh_id,
		int neltot, int nproc, int m)
{
	for (int i = 0; i < m; ++i) {
		neigh_el[i] = elnum_of_id(neigh_id[i], neltot, nproc);
		neigh_proc[i] = proc_of_id(neigh_id[i], neltot, nproc);
	}
}

void set_boundary_conditions(struct geometry *geom)
{
	int diric_bc[6] = {DIRICHLET_W, DIRICHLET_E, DIRICHLET_S, DIRICHLET_N,
			DIRICHLET_B, DIRICHLET_T};
	for (int i = 0; i < 6; ++i)
		if (geom->neigh_el[i] == -1)
			geom->bc = geom->bc | diric_bc[i];
}
