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

#ifndef TEST_H
#define TEST_H

#include "structs.h"

double zero_func(double x, double y, double z, double t);
double one_func(double x, double y, double z, double t);
double max_error(struct element **sel, int var, double (*exact_solution) ());
double maxmin_val(struct element **sel, int var, double maxmin);
double max_divergence(struct element **sel);
void set_exact_pressure(struct element **sel, double (*exact_solution) (),
		double ****p);
double pressure_error(struct element **sel, double (*exact_solution) ());
void initialise_test(struct element **sel);

#endif
