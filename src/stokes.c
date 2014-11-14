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
#include "linalg.h"
#include "tripalg.h"
#include "basis.h"
#include "deriv.h"
#include "poisson.h"
#include "geometry.h"
#include "convection.h"
#include "time.h"

void discrete_gradient_operator(struct element *sel, double ***p,
		struct trivec *grad_p)
{
	struct trivec *tmp1 = new_trivec(sel->n);
	double ***tmp2 = new_3d_array(sel->n, sel->n, sel->n);

	test_function_inner_product(sel, p, tmp2);
	
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < 3; ++i)
			multiply(***sel->geom->tinv_jac->ind[j][i], **tmp2,
					***tmp1->ind[i], sel->ntot, ADD_NEW,
					REPLACE_OLD);
		transpose_divergence(sel->basis, tmp1, *grad_p->ind[j]);
	}

	scale_trivec(grad_p, sel->ntot, -1.0);

	free_trivec(tmp1);
	free_3d_array(tmp2);
}

void discrete_divergence_operator(struct element *sel, struct trivec *u,
		double ***div_u)
{
	struct trimat *tmp = new_trimat(sel->n);
	const double old[2] = { REPLACE_OLD, ADD_OLD };

	compute_jacobian(sel->basis, u, tmp);

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			multiply(***tmp->ind[i][j], ***sel->geom->tinv_jac
					->ind[i][j], **div_u, sel->ntot,
					ADD_NEW, old[i + j != 0]);

	test_function_inner_product(sel, div_u, div_u);

	free_trimat(tmp);
}

void discrete_curl_operator(struct element *sel, struct trivec *u,
		struct trivec *omega)
{
	struct trimat *tmp = new_trimat(sel->n);
	const double old[3] = { REPLACE_OLD, ADD_OLD, ADD_OLD };
	const double new[3] = { ADD_NEW, SUBTRACT_NEW, ADD_NEW };

	compute_jacobian(sel->basis, u, tmp);

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			const int k = (5 - j) / 2;
			const int l = (2 - j) / 2;
			multiply(***tmp->ind[k][i], ***sel->geom->tinv_jac
					->ind[l][i], ***omega->ind[j],
					sel->ntot, new[j], old[i]);
			multiply(***tmp->ind[l][i], ***sel->geom->tinv_jac
					->ind[k][i], ***omega->ind[j],
					sel->ntot, -new[j], ADD_OLD);
		}
	free_trimat(tmp);
}

void double_curl_operator(struct element *sel, struct trivec *u,
		struct trivec *nxnxu)
{
	struct trivec *tmp = new_trivec(sel->basis->n);
	discrete_curl_operator(sel, u, tmp);
	discrete_curl_operator(sel, tmp, nxnxu);

	for (int i = 0; i < 3; ++i)
		test_function_inner_product(sel, *nxnxu->ind[i],
				*nxnxu->ind[i]);
	free_trivec(tmp);
}

void estimate_pressure_gradient(struct element *sel, struct trivec *dpd)
{
	const int n = sel->basis->n;
	double ***tmp = new_3d_array(n, n, n);
	copy_trivec(sel->f, dpd, sel->ntot);

	for (int i = 0; i < 3; ++i) {
		set_dirichlet_boundaries(sel, tmp, sel->bc,
				sel->bc_dfunc[i]);
		add_array(**tmp, ***dpd->ind[i], n * n * n, SUBTRACT_NEW);
		test_function_inner_product(sel, *dpd->ind[i], *dpd->ind[i]);
	}
	free_3d_array(tmp);

	time_shift(sel->rot, sel->old->rot, sel->old->old->rot, sel->ntot);
	double_curl_operator(sel, sel->u, sel->rot);

	stiffly_stable(sel->conv, sel->old->conv, sel->old->old->conv, dpd,
			sel->ntot, SUBTRACT_NEW, EXPLICIT3);
	stiffly_stable(sel->rot, sel->old->rot, sel->old->old->rot, dpd,
			sel->ntot, -sel->params->visc, EXPLICIT3);

	scale_trivec(dpd, sel->ntot, sel->params->dens);
}
