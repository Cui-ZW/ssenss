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

#ifndef CONVECTION_H
#define CONVECTION_H

void convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv);
void antialiased_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv);
void convective_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv);
void conservative_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv);
void skew_symmetric_convection_operator(struct element *sel, struct trivec *u,
		struct trivec *v, struct trivec *udv);
#endif
