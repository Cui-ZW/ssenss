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

#include <stdlib.h>
#include "structs.h"
#include "utils.h"
#include "linalg.h"
#include "tripalg.h"

struct trivec *new_trivec(int n)
{
	struct trivec *trip = ec_malloc(sizeof(*trip));

	trip->x = new_3d_array(n, n, n);
	trip->y = new_3d_array(n, n, n);
	trip->z = new_3d_array(n, n, n);
	trip->ind[0] = &trip->x;
	trip->ind[1] = &trip->y;
	trip->ind[2] = &trip->z;

	return trip;
}

void free_trivec(struct trivec *trip)
{
	free_3d_array(trip->x);
	free_3d_array(trip->y);
	free_3d_array(trip->z);
	free(trip);
}

struct trimat *new_trimat(int n)
{
	struct trimat *ttrip = ec_malloc(sizeof(*ttrip));

	ttrip->ind[0][0] = &ttrip->x.x;
	ttrip->ind[0][1] = &ttrip->x.y;
	ttrip->ind[0][2] = &ttrip->x.z;
	ttrip->ind[1][0] = &ttrip->y.x;
	ttrip->ind[1][1] = &ttrip->y.y;
	ttrip->ind[1][2] = &ttrip->y.z;
	ttrip->ind[2][0] = &ttrip->z.x;
	ttrip->ind[2][1] = &ttrip->z.y;
	ttrip->ind[2][2] = &ttrip->z.z;

	ttrip->x.ind[0] = &ttrip->x.x;
	ttrip->x.ind[1] = &ttrip->x.y;
	ttrip->x.ind[2] = &ttrip->x.z;
	ttrip->y.ind[0] = &ttrip->y.x;
	ttrip->y.ind[1] = &ttrip->y.y;
	ttrip->y.ind[2] = &ttrip->y.z;
	ttrip->z.ind[0] = &ttrip->z.x;
	ttrip->z.ind[1] = &ttrip->z.y;
	ttrip->z.ind[2] = &ttrip->z.z;

	ttrip->row[0] = &ttrip->x;
	ttrip->row[1] = &ttrip->y;
	ttrip->row[2] = &ttrip->z;

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			*ttrip->ind[i][j] = new_3d_array(n, n, n);

	return ttrip;
}

void free_trimat(struct trimat *ttrip)
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			free_3d_array(*ttrip->ind[i][j]);

	free(ttrip);
}

struct trimat *new_sym_trimat(int n)
{
	struct trimat *ttrip = ec_malloc(sizeof(*ttrip));

	ttrip->ind[0][0] = &ttrip->x.x;
	ttrip->ind[0][1] = &ttrip->x.y;
	ttrip->ind[0][2] = &ttrip->x.z;
	ttrip->ind[1][0] = &ttrip->y.x;
	ttrip->ind[1][1] = &ttrip->y.y;
	ttrip->ind[1][2] = &ttrip->y.z;
	ttrip->ind[2][0] = &ttrip->z.x;
	ttrip->ind[2][1] = &ttrip->z.y;
	ttrip->ind[2][2] = &ttrip->z.z;

	ttrip->x.ind[0] = &ttrip->x.x;
	ttrip->x.ind[1] = &ttrip->x.y;
	ttrip->x.ind[2] = &ttrip->x.z;
	ttrip->y.ind[0] = &ttrip->y.x;
	ttrip->y.ind[1] = &ttrip->y.y;
	ttrip->y.ind[2] = &ttrip->y.z;
	ttrip->z.ind[0] = &ttrip->z.x;
	ttrip->z.ind[1] = &ttrip->z.y;
	ttrip->z.ind[2] = &ttrip->z.z;

	ttrip->row[0] = &ttrip->x;
	ttrip->row[1] = &ttrip->y;
	ttrip->row[2] = &ttrip->z;

	ttrip->x.x = new_3d_array(n, n, n);
	ttrip->x.y = new_3d_array(n, n, n);
	ttrip->x.z = new_3d_array(n, n, n);

	ttrip->y.x = ttrip->x.y;
	ttrip->y.y = new_3d_array(n, n, n);
	ttrip->y.z = new_3d_array(n, n, n);

	ttrip->z.x = ttrip->x.z;
	ttrip->z.y = ttrip->y.z;
	ttrip->z.z = new_3d_array(n, n, n);

	return ttrip;
}

void free_sym_trimat(struct trimat *ttrip)
{
	free_3d_array(ttrip->x.x);
	free_3d_array(ttrip->x.y);
	free_3d_array(ttrip->x.z);

	free_3d_array(ttrip->y.y);
	free_3d_array(ttrip->y.z);

	free_3d_array(ttrip->z.z);
	free(ttrip);
}

void add_trivec(struct trivec *u, struct trivec *v, int n, double alpha)
{
	for (int i = 0; i < 3; ++i)
		add_array(***u->ind[i], ***v->ind[i], n, alpha);
}

void scale_trivec(struct trivec *u, int n, double alpha)
{
	for (int i = 0; i < 3; ++i)
		scale_array(***u->ind[i], n, alpha);
}

void copy_trivec(struct trivec *u, struct trivec *v, int n)
{
	for (int i = 0; i < 3; ++i)
		copy_array(***u->ind[i], ***v->ind[i], n);
}

void compute_trimat_determinant(struct trimat *a, int n, double ***ad)
{
	double ***tmp = new_3d_array(n, n, n);
	const size_t m = ((size_t) n) * ((size_t) n) * ((size_t) n);
	const double new[3] = { ADD_NEW, SUBTRACT_NEW, ADD_NEW };
	const double old[3] = { REPLACE_OLD, ADD_OLD, ADD_OLD };

	for (int i = 0; i < 3; ++i) {
		const int j = (2 - i) / 2;
		const int k = 3 - i - j;
		multiply(***a->ind[1][j], ***a->ind[2][k], **tmp, m, ADD_NEW,
				REPLACE_OLD);
		multiply(***a->ind[1][k], ***a->ind[2][j], **tmp, m,
				SUBTRACT_NEW, ADD_OLD);
		multiply(***a->ind[0][i], **tmp, **ad, m, new[i], old[i]);
	}
	free_3d_array(tmp);
}

void transpose_trimat(struct trimat *a)
{
	double ***ptr;

	ptr = a->x.y;
	a->x.y = a->y.x;
	a->y.x = ptr;

	ptr = a->x.z;
	a->x.z = a->z.x;
	a->z.x = ptr;

	ptr = a->y.z;
	a->y.z = a->z.y;
	a->z.y = ptr;
}

void invert_trimat(struct trimat *a, double ***ad, int n, struct trimat *ai)
{
	const size_t m = ((size_t) n) * ((size_t) n) * ((size_t) n);
	double ***tmp = new_3d_array(n, n, n);
	const int k[4][3] = { { 1, 0, 0 }, { 2, 2, 1 }, { 2, 0, 1 },
			{ 1, 2, 0 } };
	const double new[3] = { ADD_NEW, SUBTRACT_NEW, ADD_NEW };
	

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			multiply(***a->ind[k[0][j]][k[3][i]],
					***a->ind[k[1][j]][k[2][i]], **tmp, m,
					ADD_NEW, REPLACE_OLD);
			multiply(***a->ind[k[0][j]][k[2][i]],
					***a->ind[k[1][j]][k[3][i]], **tmp, m,
					SUBTRACT_NEW, ADD_OLD);
			divide(**tmp, **ad, ***ai->ind[i][j], m, new[j],
					REPLACE_OLD);
		}
	free_3d_array(tmp);
}

void trimat_trivec_product(struct trimat *a, struct trivec *b, int n,
		struct trivec *c)
{
	const size_t m = ((size_t) n) * ((size_t) n) * ((size_t) n);
	const double old[3] = { REPLACE_OLD, ADD_OLD, ADD_OLD };

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			multiply(***a->ind[i][j], ***b->ind[j], ***c->ind[i], m,
					ADD_NEW, old[j]);
}
