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

#include <math.h>
#include "utils.h"
#include "structs.h"
#include "linalg.h"
#include "tripalg.h"
#include "basis.h"
#include "geometry.h"
#include "poisson.h"
#include "global.h"
#include "helmholtz.h"

static double helmholtz_diagonal(struct element *sel, int i, int j, int k,
		double lambda);

void discrete_helmholtz_operator(struct element *sel, double ***u,
		double ***helm_u)
{
	double ***tmp = new_3d_array(sel->basis->n, sel->basis->n,
			sel->basis->n);
	const double lambda = (11.0 / 6.0) / (sel->params->dt
			* sel->params->visc);

	double **tmp_deriv = sel->basis->deriv_matrix;
	sel->basis->deriv_matrix = sel->basis->svv_deriv_matrix;
	discrete_laplace_operator(sel, u, helm_u);
	sel->basis->deriv_matrix = tmp_deriv;

	test_function_inner_product(sel, u, tmp);
	add_array(**tmp, **helm_u, sel->basis->ntot, lambda);
	free_3d_array(tmp);
}

static double helmholtz_diagonal(struct element *sel, int i, int j, int k,
		double lambda)
{
	double **d = sel->basis->deriv_matrix;
	double ***g11 = sel->geom->geom_mat->x.x;
	double ***g22 = sel->geom->geom_mat->y.y;
	double ***g33 = sel->geom->geom_mat->z.z;
	double ***rho = sel->basis->weights_array;
	
	double tmp = 0.0;
	for (int p = 0; p < sel->basis->n; ++p) {
		tmp += rho[k][j][p] * pow(d[i][p], 2.0) * g11[k][j][p];
		tmp += rho[k][p][i] * pow(d[j][p], 2.0) * g22[k][p][i];
		tmp += rho[p][j][i] * pow(d[k][p], 2.0) * g33[p][j][i];
	}
	tmp += lambda * rho[k][j][i] * sel->geom->jac_det[k][j][i];
	return tmp;
}

void set_helmholtz_diagonal(struct element *sel, double ***x, double lambda)
{
	for (int k = 0; k < sel->basis->n; ++k)
		for (int j = 0; j < sel->basis->n; ++j)
			for (int i = 0; i < sel->basis->n; ++i)
				x[k][j][i] = helmholtz_diagonal(sel, i, j, k,
						lambda);
}
