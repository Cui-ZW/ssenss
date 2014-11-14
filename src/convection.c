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

#include "utils.h"
#include "structs.h"
#include "basis.h"
#include "linalg.h"
#include "tripalg.h"
#include "geometry.h"
#include "deriv.h"
#include "time.h"
#include "poisson.h"
#include "interpolation.h"
#include "convection.h"

void convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv)
{
	antialiased_convection_operator(sel, u, v, udv);
}

void antialiased_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv)
{
	const int n = sel->conv_basis->n;
	struct trivec *u_int = new_trivec(n);
	struct trivec *v_int = new_trivec(n);
	struct trivec *udv_int = new_trivec(n);

	struct basis * const tmp_basis = sel->basis;
	struct geometry * const tmp_geometry = sel->geom;
	sel->basis = sel->conv_basis;
	sel->geom = sel->conv_geom;

	for (int i = 0; i < 3; ++i) {
		interpolate_to_grid(sel->conv_basis, *u_int->ind[i], tmp_basis,
				*u->ind[i]);
		interpolate_to_grid(sel->conv_basis, *v_int->ind[i], tmp_basis,
				*v->ind[i]);
	}

	skew_symmetric_convection_operator(sel, u_int, v_int, udv_int);

	for (int i = 0; i < 3; ++i)
		integrate_to_grid(tmp_basis, *udv->ind[i], sel->conv_basis,
				*udv_int->ind[i]);

	sel->basis = tmp_basis;
	sel->geom = tmp_geometry;

	free_trivec(u_int);
	free_trivec(v_int);
	free_trivec(udv_int);
}

void convective_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv)
{
	const int n = sel->basis->n;
	const int ntot = sel->basis->ntot;
	struct trimat *tmp = new_trimat(n);
	struct trimat *dv = new_trimat(n);
	const double old[3] = { REPLACE_OLD, ADD_OLD, ADD_OLD };

	compute_jacobian(sel->basis, v, tmp);

	for (int j = 0; j < 3; ++j) {
		trimat_trivec_product(sel->geom->tinv_jac, tmp->row[j], n,
				dv->row[j]);
		for (int i = 0; i < 3; ++i)
			multiply(***u->ind[i], ***dv->ind[j][i], ***udv->ind[j],
					ntot, ADD_NEW, old[i]);
		test_function_inner_product(sel, *udv->ind[j], *udv->ind[j]);
	}

	free_trimat(tmp);
	free_trimat(dv);
}

void conservative_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv)
{
	const int n = sel->basis->n;
	const int ntot = sel->basis->ntot;
	struct trivec *tmp = new_trivec(n);
	struct trivec *duvd = new_trivec(n);
	double ***uv = new_3d_array(n, n, n);

	for (int i = 0; i < 3; ++i) {
		multiply(**u->x, ***v->ind[i], **uv, ntot, ADD_NEW,
				REPLACE_OLD);
		discrete_x_differentiation(sel->basis, uv, tmp->x, ADD_NEW,
				REPLACE_OLD);
		multiply(**u->y, ***v->ind[i], **uv, ntot, ADD_NEW,
				REPLACE_OLD);
		discrete_y_differentiation(sel->basis, uv, tmp->y, ADD_NEW,
				REPLACE_OLD);
		multiply(**u->z, ***v->ind[i], **uv, ntot, ADD_NEW,
				REPLACE_OLD);
		discrete_z_differentiation(sel->basis, uv, tmp->z, ADD_NEW,
				REPLACE_OLD);

		trimat_trivec_product(sel->geom->tinv_jac, tmp, n, duvd);
		copy_array(**duvd->x, ***udv->ind[i], ntot);
		add_array(**duvd->y, ***udv->ind[i], ntot, ADD_NEW);
		add_array(**duvd->z, ***udv->ind[i], ntot, ADD_NEW);
		test_function_inner_product(sel, *udv->ind[i], *udv->ind[i]);
	}
	free_trivec(tmp);
	free_trivec(duvd);
	free_3d_array(uv);
}

void skew_symmetric_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv)
{
	convective_convection_operator(sel, u, v, udv);
	struct trivec *tmp = new_trivec(sel->basis->n);
	conservative_convection_operator(sel, u, v, tmp);

	scale_trivec(udv, sel->basis->ntot, 0.5);
	add_trivec(tmp, udv, sel->basis->ntot, 0.5);

	free_trivec(tmp);
}
