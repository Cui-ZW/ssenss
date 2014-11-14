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

#ifndef TIME_H
#define TIME_H

#include "structs.h"

enum { IMPLICIT, EXPLICIT2, EXPLICIT3 };

void time_shift(struct trivec *u, struct trivec *old, struct trivec *old_old,
		int n);
void single_time_shift(double ***data, double ***old, double ***old_old, int n);
void stiffly_stable(struct trivec *u1, struct trivec *u2, struct trivec *u3,
		struct trivec  *u, size_t ntot, double c, int expl);
void single_stiffly_stable(double ***u1, double ***u2, double ***u3,
		double ***u, size_t ntot, double c, int imp_exp);
void integrate_convection_se3(struct element **sel, struct trivec **u);
void integrate_convection_rk4(struct element **sel, struct trivec **u,
		double t0, double t, double dt, double beta);
void add_interpolated_velocity(struct element **sel, struct trivec **u,
		double t, double alpha);
void trivec_inverse_mass_matrix(struct element **sel, struct trivec **u);

#endif
