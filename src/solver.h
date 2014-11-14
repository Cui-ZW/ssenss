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

#ifndef SOLVER_H
#define SOLVER_H

#include "structs.h"

double conjugate_gradient(struct element **sel, double ****f, double ****u,
		int var, void (*opeval) (), void (*precon) (), double eps,
		int max_iter);
double flexible_conjugate_gradient(struct element **sel, double ****f,
		double ****u, int var, void (*opeval) (), void (*precon) (),
		double eps, int max_iter);
void solve_poisson_cg(struct element **sel, double ****f, int var, double eps);
void solve_helmholtz_cg(struct element **sel, double ****f, int var,
		double eps);
void global_opeval(struct element **sel, double ****p, double ****ap, int var,
		void (*opeval) ());
void subtract_nullspace(struct element **sel, double ****u, int var);

#endif
