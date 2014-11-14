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

#ifndef POISSON_H
#define POISSON_H

#include "structs.h"

enum { DIRICHLET_W = 01, DIRICHLET_E = 02, DIRICHLET_S = 04, DIRICHLET_N = 010,
		DIRICHLET_B = 020, DIRICHLET_T = 040 };
enum { X_DERIV, Y_DERIV, Z_DERIV };
enum { PRESSURE, X_VELOCITY, Y_VELOCITY, Z_VELOCITY };
enum { X_SOURCE, Y_SOURCE, Z_SOURCE };

void discrete_laplace_operator(struct element *sel, double ***u,
		double ***lap_u);
void test_function_inner_product(struct element *sel, double ***u,
		double ***vu);
void inverse_mass_matrix_operator(struct element *sel, double ***u,
		double ***vu);
void global_inverse_mass_matrix_operator(struct element **sel, double ****u,
		double ****mu);
void mask_dirichlet_boundaries(struct element **sel, double ****u, int var);
void set_dirichlet_boundaries(struct element *sel, double ***u, int bc,
		double (*bc_func) ());
void add_surface_integral(struct element *sel, double ***g, double ***f,
		int deriv, int bc);
double ***variable_array(struct element *sel, int var);
double ***source_array(struct element *sel, int var);
void set_lifted_velocity(struct element **sel, double t, double (*ud) (),
		double (*vd) (), double (*wd) ());
void initialise_from_function(struct element **sel, double (*ud) (),
		double (*vd) (), double (*wd) ());
void initialise_from_file(struct element **sel, char *filename);
void set_source_term(struct element **sel, int var, double (*f) (),
		double alpha);

#endif
