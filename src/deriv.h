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

#ifndef DERIV_H
#define DERIV_H

#include "structs.h"

void discrete_x_differentiation(struct basis *disc_basis, double ***u,
		double ***du, double alpha, double beta);
void discrete_y_differentiation(struct basis *disc_basis, double ***u,
		double ***du, double alpha, double beta);
void discrete_z_differentiation(struct basis *disc_basis, double ***u,
		double ***du, double alpha, double beta);
void transpose_divergence(struct basis *disc_basis, struct trivec *u,
		double ***du);
void compute_jacobian(struct basis *disc_basis, struct trivec *grid,
		struct trimat *jac);
void compute_svv_derivative_matrix(struct basis *basis, int m, double c);

#endif
