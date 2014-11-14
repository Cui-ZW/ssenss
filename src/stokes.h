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

#ifndef STOKES_H
#define STOKES_H

#include "structs.h"

void discrete_gradient_operator(struct element *sel, double ***p,
		struct trivec *grad_p);
void discrete_divergence_operator(struct element *sel, struct trivec *u,
		double ***div_u);
void discrete_curl_operator(struct element *sel, struct trivec *u,
		struct trivec *omega);
void double_curl_operator(struct element *sel, struct trivec *u,
		struct trivec *nxnxu);
void estimate_pressure_gradient(struct element *sel, struct trivec *dpd);

#endif
