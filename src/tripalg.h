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

#ifndef TRIPALG_H
#define TRIPALG_H

struct trivec {
	double ***x;
	double ***y;
	double ***z;
	double ****ind[3];
};

struct trimat {
	struct trivec x;
	struct trivec y;
	struct trivec z;
	double ****ind[3][3];
	struct trivec *row[3];
};

struct trivec *new_trivec(int n);
void free_trivec(struct trivec *trip);
struct trimat *new_trimat(int n);
void free_trimat(struct trimat *ttrip);
struct trimat *new_sym_trimat(int n);
void free_sym_trimat(struct trimat *ttrip);
void add_trivec(struct trivec *u, struct trivec *v, int n, double alpha);
void scale_trivec(struct trivec *u, int n, double alpha);
void copy_trivec(struct trivec *u, struct trivec *v, int n);
void compute_trimat_determinant(struct trimat *a, int n, double ***ad);
void transpose_trimat(struct trimat *a);
void invert_trimat(struct trimat *a, double ***ad, int n, struct trimat *ai);
void trimat_trivec_product(struct trimat *a, struct trivec *b, int n,
		struct trivec *c);

#endif
