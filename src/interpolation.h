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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "structs.h"

double lagrangian_interpolation(double *x, double ***f, double xp, double yp,
		double zp, int n);
double polynomial_interpolation(double *x, double ***f, double xp, double yp,
		double zp, int n);
void create_interpolation_matrices(struct basis *basis1, struct basis *basis2);
void interpolate_to_grid(struct basis *to_basis, double ***int_u,
		struct basis *from_basis, double ***u);
void integrate_to_grid(struct basis *to_basis, double ***int_u,
		struct basis *from_basis, double ***u);
#endif
