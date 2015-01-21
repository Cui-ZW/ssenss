/*
 * Copyright 2015 Christopher Nilsen
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
#include <math.h>
#include <mpi.h>
#include "utils.h"
#include "tripalg.h"
#include "basis.h"
#include "interpolation.h"
#include "geometry.h"
#include "poisson.h"
#include "global.h"
#include "particles.h"

#define BUFFSIZE(N) (((N) * 12) / 10)

static void get_elemental_particle_count(struct element **sel, size_t np,
		size_t *npel);
static int crossed_boundary(struct particle *p, double y[6]);
static void compute_particle_equation_rhs(struct particle *p,
		struct element *sel, double y[6], double t, double f[6]);
static void integrate_particle_rk4(struct particle *p, struct element *sel,
		double dt, int flag);
static void integrate_to_boundary(struct particle *p, struct element *sel,
		double dt);
static double timestep_to_boundary(double y[6], double f[6]);
static int unfinished_particles(struct plist *p);
static void find_new_particle_element(struct particle *p, struct element *sel);
static void find_new_particle_location(struct particle *p, struct element *sel);
static int closest_face(struct particle *p);
static void rotate_face_coordinates(double x[2], int rot);
static void relocate_particles(struct plist *p);
static void get_particle_comm_sizes(struct plist *p, int *sendcounts,
		int *recvcounts, int *senddispls, int *recvdispls);
static void prepare_particle_buffers(struct plist *p, int *recvcounts,
		struct particle **sendbuff, struct particle **recvbuff);
static void remove_old_particles(struct plist *p);
static double fluid_point_velocity(struct element *sel, double y[6], double t,
		int comp);
static void sort_particles(struct plist *p, int (*compare_func) ());
static void sort_particles_for_sending(struct plist *p);

enum { P_NORMAL, P_COLLISION, P_OUTSIDE, P_NEW_PROC, P_BUFFER };
enum { RK4_NORMAL, RK4_FORCE };

struct part_params *new_particle_parameters(void)
{
	struct part_params *pparams = ec_malloc(sizeof(*pparams));
	pparams->n = 100;
	pparams->ntypes = 10;
	pparams->tau_min = 0.1;
	pparams->tau_max = 1.0;

	return pparams;
}

void free_particle_parameters(struct part_params *pparams)
{
	free(pparams);
}

struct plist *new_particle_list(struct part_params *pparams,
		struct element **sel)
{
	struct plist *p = ec_malloc(sizeof(*p));
	p->n = local_particle_count(sel, pparams->n);
	p->ntot = global_int_sum((*sel)->comm->mpi_comm, p->n);
	p->nbuff = BUFFSIZE(p->n);
	p->mpi_comm = (*sel)->comm->mpi_comm;
	p->nproc = (*sel)->params->nproc;
	p->list = ec_malloc(sizeof(*p->list) * p->nbuff);

	for (size_t i = 0; i < p->n; ++i) {
		p->list[i].status = P_NORMAL;
		p->list[i].proc = (*sel)->comm->rank;
		p->list[i].old_proc = (*sel)->comm->rank;
	}
	for (size_t i = p->n; i < p->nbuff; ++i)
		p->list[i].status = P_BUFFER;

	p->ntypes = pparams->ntypes;
	p->tau = new_1d_array(p->ntypes);
	get_particle_timescales(pparams->tau_min, pparams->tau_max, p->ntypes,
			p->tau);
	return p;
}

void free_particle_list(struct plist *p)
{
	free(p->list);
	free_1d_array(p->tau);
	free(p);
}

size_t local_particle_count(struct element **sel, size_t ntot)
{
	double vol = 0.0;
	for (int i = 0; i < (*sel)->nel; ++i)
		vol += element_volume(sel[i]);
	double totvol = global_sum(sel, vol);
	size_t n = (size_t) ((double) ntot * vol / totvol);
	return n;
}

static void get_elemental_particle_count(struct element **sel, size_t np,
		size_t *npel)
{
	double *vol = new_1d_array((*sel)->nel);
	double totvol = 0.0;
	for (int i = 0; i < (*sel)->nel; ++i) {
		vol[i] = element_volume(sel[i]);
		totvol += vol[i];
	}
	size_t nsum = 0;
	for (int i = 0; i < (*sel)->nel - 1; ++i) {
		npel[i] = (size_t) ((double) np * vol[i] / totvol);
		nsum += npel[i];
	}
	npel[(*sel)->nel - 1] = np - nsum;
	free_1d_array(vol);
}

void get_particle_timescales(double tau_min, double tau_max, int n, double *tau)
{
	const double a = 0.5 * (log(tau_max) + log(tau_min));
	const double b = 0.5 * (log(tau_max) - log(tau_min));
	double *xi = gll_gridpoints(n);

	for (int i = 0; i < n; ++i)
		tau[i] = exp(a + xi[i] * b);
	free_1d_array(xi);
}

void initialise_particles(struct plist *p, struct element **sel, size_t nel)
{
	srand(((unsigned int) MPI_Wtime()) * (1 + (*sel)->params->rank));
	size_t *npel = ec_malloc((*sel)->nel * sizeof(*npel));
	get_elemental_particle_count(sel, p->n, npel);
	size_t i = 0;
	for (int j = 0; j < (*sel)->nel; ++j)
		for (int k = 0; k < npel[j]; ++k) {
			p->list[i].id = p->n * (*sel)->params->rank + i;
			p->list[i].el = j;
			p->list[i].old_el = p->list[i].el;
			initialise_particle_position(&p->list[i], sel[j]);
			p->list[i].t = 0.0;
			p->list[i].type = k % p->ntypes;
			p->list[i].tau = p->tau[p->list[i].type];
			++i;
		}
	free(npel);
}

void initialise_particle_position(struct particle *p, struct element *sel)
{
	double tmp[3] = { 2.0 * random_double() - 1.0, 2.0 * random_double()
			- 1.0, 2.0 * random_double() - 1.0 };
	p->x = tmp[0];
	p->y = tmp[1];
	p->z = tmp[2];
	p->u = fluid_point_velocity(sel, tmp, sel->params->t, 0);
	p->v = fluid_point_velocity(sel, tmp, sel->params->t, 1);
	p->w = fluid_point_velocity(sel, tmp, sel->params->t, 2);
}

void integrate_particle_equations(struct plist *p, struct element **sel)
{
	for (size_t i = 0; i < p->n; ++i)
		p->list[i].status = P_OUTSIDE;
	int j = 0;
	do {
		for (size_t i = 0; i < p->n; ++i) {
			if (p->list[i].status == P_NORMAL)
				continue;
			p->list[i].status = P_NORMAL;
			if (j > 10) {
				initialise_particle_position(&p->list[i],
						sel[p->list[i].el]);
				p->list[i].t = (*sel)->params->t;
				printf("Particle %d reinitialised\n",
						p->list[i].id);
				continue;
			}
			integrate_particle_rk4(&p->list[i], sel[p->list[i].el],
					(*sel)->params->t - p->list[i].t,
					RK4_NORMAL);
			if (p->list[i].status == P_OUTSIDE)
				integrate_to_boundary(&p->list[i],
						sel[p->list[i].el],
						2.0 * (*sel)->params->dt);
			if (p->list[i].status == P_OUTSIDE)
				find_new_particle_element(&p->list[i],
						sel[p->list[i].el]);
		}
		relocate_particles(p);
		for (size_t i = 0; i < p->n; ++i)
			if (p->list[i].status == P_OUTSIDE)
				find_new_particle_location(&p->list[i],
						sel[p->list[i].el]);
		++j;
	} while (unfinished_particles(p));
}

static int crossed_boundary(struct particle *p, double y[6])
{
	for (int i = 0; i < 3; ++i)
		if (fabs(y[i]) > 1.0)
			return 1;
	return 0;
}

static void compute_particle_equation_rhs(struct particle *p,
		struct element *sel, double y[6], double t, double f[6])
{
	const double a = 1.0 / p->tau;
	for (int i = 0; i < 3; ++i) {
		f[i] = y[i + 3];
		f[i + 3] = a * (fluid_point_velocity(sel, y, t, i) - y[i + 3]);
	}
}

static void integrate_particle_rk4(struct particle *p, struct element *sel,
		double dt, int flag)
{
	double t;
	const double a[4] = { 0.0, 0.5, 0.5, 1.0 };
	double y[6];
	double y0[6] = { p->x, p->y, p->z, p->u, p->v, p->w };
	double *yn[6] = { &p->x, &p->y, &p->z, &p->u, &p->v, &p->w };
	double f[4][6];
	
	for (int j = 0; j < 4; ++j) {
		t = p->t + a[j] * dt;
		for (int i = 0; i < 6; ++i)
			y[i] = y0[i] + a[j] * dt * f[max(j - 1, 0)][i];
		if (j > 0)
			if (crossed_boundary(p, y) && flag == RK4_NORMAL) {
				p->status = P_OUTSIDE;
				return;
			}
		compute_particle_equation_rhs(p, sel, y, t, f[j]);
	}
	for (int i = 0; i < 6; ++i)
		y[i] = y0[i] +  dt * (f[0][i] / 6.0 + f[1][i] / 3.0
				+ f[2][i] / 3.0 + f[3][i] / 6.0);
	if (crossed_boundary(p, y)) {
		p->status = P_OUTSIDE;
		return;
	}
	for (int i = 0; i < 6; ++i)
		*yn[i] = y[i];
	p->t = sel->params->t;
}

static void integrate_to_boundary(struct particle *p, struct element *sel,
		double dt)
{
	double t;
	const double a[4] = { 0.0, 0.5, 0.5, 1.0 };
	double y[6];
	double y0[6] = { p->x, p->y, p->z, p->u, p->v, p->w };
	double *yn[6] = { &p->x, &p->y, &p->z, &p->u, &p->v, &p->w };
	double f[5][6];
	double dtau = 0.0;

	for (int j = 0; j < 4; ++j) {
		t = p->t + a[j] * dtau;
		for (int i = 0; i < 6; ++i)
			y[i] = y0[i] + a[j] * dtau * f[max(j - 1, 0)][i];
		compute_particle_equation_rhs(p, sel, y, t, f[j]);
		dtau = timestep_to_boundary(y0, f[j]);
	}
	for (int i = 0; i < 6; ++i)
		f[4][i] = (f[0][i] / 6.0 + f[1][i] / 3.0 + f[2][i] / 3.0
				+ f[3][i] / 6.0);
	dtau = timestep_to_boundary(y0, f[4]);

	if (dtau < dt) {
		p->t += dtau;
		for (int i = 0; i < 6; ++i)
			*yn[i] += dtau * f[4][i];
	} else {
		p->status = P_NORMAL;
		integrate_particle_rk4(p, sel, dt, RK4_FORCE);
		if (p->status != P_NORMAL) {
			p->status = P_NORMAL;
			initialise_particle_position(p, sel);
			p->t = sel->params->t;
			printf("Particle %d reinitialised\n", p->id);
		}
	}
}

static double timestep_to_boundary(double y[6], double f[6])
{
	double dt[3];
	for (int i = 0; i < 3; ++i)
		dt[i] = (copysign(1.0, f[i]) - y[i]) / f[i];
	return dmin(dmin(dt[0], dt[1]), dt[2]);
}

static int unfinished_particles(struct plist *p)
{
	int tmp[2] = { 0, 0 };
	for (size_t i = 0; i < p->n; ++i)
		if (p->list[i].status != P_NORMAL) {
			tmp[0] = 1;
			break;
		}
	int error = MPI_Allreduce(tmp, tmp + 1, 1, MPI_INT, MPI_SUM,
			p->mpi_comm);
	if (error)
		fatal_error("in relocate_particles for MPI_Allreduce");
	return tmp[1];
}

static void find_new_particle_element(struct particle *p, struct element *sel)
{
	double *u[3] = { &p->u, &p->v, &p->w };
	int j = closest_face(p);

	if (sel->geom->neigh_el[j] >= 0) {
		p->el = sel->geom->neigh_el[j];
		p->proc = sel->geom->neigh_proc[j];
		for (int i = 0; i < 3; ++i)
			*u[i] = physical_particle_velocity(p, sel, i);
		if (p->proc != sel->params->rank)
			p->status = P_NEW_PROC;
	} else {
		p->status = P_COLLISION;
		*u[j / 2] *= -1.0;
	}
}

static void find_new_particle_location(struct particle *p, struct element *sel)
{
	double y0[6] = { p->x, p->y, p->z, p->u, p->v, p->w };
	double *yn[6] = { &p->x, &p->y, &p->z, &p->u, &p->v, &p->w };
	double x[2];
	int k = closest_face(p);

	int j = 0;
	for (int i = 0; i < 3; ++i)
		if (i != k / 2)
			x[j++] = y0[i];

	for (int i = 0; i < 6; ++i)
		if (sel->geom->neigh_el[i] == p->old_el &&
				sel->geom->neigh_proc[i] == p->old_proc)
			k = i;
	*yn[k / 2] = pow(-1.0, k + 1);
	rotate_face_coordinates(x, sel->geom->neigh_rot[k]);

	j = 0;
	for (int i = 0; i < 3; ++i) {
		if (i != k / 2)
			*yn[i] = x[j++];
		*yn[i + 3] = parametric_particle_velocity(p, sel, i);
	}
	p->old_el = p->el;
	p->old_proc = p->proc;
}

static int closest_face(struct particle *p)
{
	double x[3] = { p->x, p->y, p->z };
	double hmin = 2.0;
	int k;
	for (int i = 0; i < 6; ++i) {
		const double h = fabs(x[i / 2] - pow(-1.0, i + 1));
		if (h < hmin) {
			hmin = h;
			k = i;
		}
	}
	return k;
}

static void rotate_face_coordinates(double x[2], int rot)
{
	for (int i = 0; i < rot; ++i) {
		double tmp = x[2];
		x[2] = x[1];
		x[1] = -tmp;
	}
}

static void relocate_particles(struct plist *p)
{
	sort_particles(p, normal_first);
	int comm[2] = { p->list[p->n - 1].status == P_NEW_PROC, 0 };
	int error = MPI_Allreduce(comm, comm + 1, 1, MPI_INT, MPI_SUM,
			p->mpi_comm);
	if (error)
		fatal_error("in relocate_particles for MPI_Allreduce");
	if (!comm[1])
		return;

	const size_t nold = p->n;
	int *sendcounts = new_int_array(p->nproc);
	int *recvcounts = new_int_array(p->nproc);
	int *senddispls = new_int_array(p->nproc);
	int *recvdispls = new_int_array(p->nproc);
	struct particle *sendbuff;
	struct particle *recvbuff;

	get_particle_comm_sizes(p, sendcounts, recvcounts, senddispls,
			recvdispls);
	prepare_particle_buffers(p, recvcounts, &sendbuff, &recvbuff);
	error = MPI_Alltoallv(sendbuff, sendcounts, senddispls, MPI_BYTE,
			recvbuff, recvcounts, recvdispls, MPI_BYTE,
			p->mpi_comm);
	if (error)
		fatal_error("in relocate_particles for MPI_Alltoallv");
	for (size_t i = nold; i < p->n; ++i)
		p->list[i].status = P_OUTSIDE;
	remove_old_particles(p);

	free_int_array(sendcounts);
	free_int_array(recvcounts);
	free_int_array(senddispls);
	free_int_array(recvdispls);
}

static void get_particle_comm_sizes(struct plist *p, int *sendcounts,
		int *recvcounts, int *senddispls, int *recvdispls)
{
	for (size_t i = 0; i < p->n; ++i)
		if (p->list[i].status == P_NEW_PROC)
			sendcounts[p->list[i].proc] += sizeof(*p->list);

	int error = MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT,
			p->mpi_comm);
	if (error)
		fatal_error("in get_particle_comm_sizes for MPI_Alltoall");

	for (int i = 1; i < p->nproc; ++i) {
		senddispls[i] = senddispls[i - 1] + sendcounts[i - 1];
		recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
	}
}

static void prepare_particle_buffers(struct plist *p, int *recvcounts,
		struct particle **sendbuff, struct particle **recvbuff)
{
	sort_particles_for_sending(p);
	int nnew = p->n;
	for (int i = 0; i < p->nproc; ++i)
		nnew += recvcounts[i] / sizeof(*p->list);
	if (nnew > p->nbuff) {
		p->nbuff = BUFFSIZE(nnew);
		p->list = ec_realloc(p->list, sizeof(*p->list) * p->nbuff); 
	}
	for (size_t i = p->n; i < p->nbuff; ++i)
		p->list[i].status = P_BUFFER;
	size_t k = -1;
	while (++k < p->n)
		if (p->list[k].status == P_NEW_PROC)
			break;
	*sendbuff = p->list + k;
	*recvbuff = p->list + p->n;
	p->n = nnew;
}

static void remove_old_particles(struct plist *p)
{
	size_t nnew = p->n;
	for (size_t i = 0; i < p->n; ++i)
		if (p->list[i].status == P_NEW_PROC) {
			p->list[i].status = P_BUFFER;
			--nnew;
		}
	sort_particles(p, normal_first);
	p->n = nnew;
}

static double fluid_point_velocity(struct element *sel, double y[6], double t,
		int comp)
{
	const double eps = 1.0e-6;
	double u[3] = { 0.0, 0.0, 0.0 };
	const double dt[3] = { t - sel->params->t + 2.0 * sel->params->dt, t -
			sel->params->t + sel->params->dt, t - sel->params->t };
	const double a = 0.5 * dt[1] * dt[2] / pow(sel->params->dt, 2);
	const double b = -dt[0] * dt[2] / pow(sel->params->dt, 2);
	const double c = 0.5 * dt[0] * dt[1] / pow(sel->params->dt, 2);

	if (fabs(a) > eps)
		u[0] = polynomial_interpolation(sel->basis->points, *sel->old
				->old->ur->ind[comp], y[0], y[1], y[2], sel->n);
	if (fabs(b) > eps)
		u[1] = polynomial_interpolation(sel->basis->points, *sel->old
				->ur->ind[comp], y[0], y[1], y[2], sel->n);
	if (fabs(c) > eps)
		u[2] = polynomial_interpolation(sel->basis->points, *sel
				->ur->ind[comp], y[0], y[1], y[2], sel->n);

	return a * u[0] + b * u[1] + c * u[2];
}

double physical_particle_velocity(struct particle *p, struct element *sel,
		int comp)
{
	double dx[3];

	for (int i = 0; i < 3; ++i)
		dx[i] = polynomial_interpolation(sel->basis->points,
				*sel->geom->jac->ind[comp][i], p->x, p->y, p->z,
				sel->n);
	return (p->u * dx[0] + p->v * dx[1] + p->w * dx[2]);
}

double parametric_particle_velocity(struct particle *p, struct element *sel,
		int comp)
{
	double dx[3];

	for (int i = 0; i < 3; ++i)
		dx[i] = polynomial_interpolation(sel->basis->points,
				*sel->geom->tinv_jac->ind[i][comp], p->x, p->y,
				p->z, sel->n);
	return (p->u * dx[0] + p->v * dx[1] + p->w * dx[2]);
}

double physical_particle_coordinate(struct particle *p, struct element *sel,
		int comp)
{
	return polynomial_interpolation(sel->basis->points,
			*sel->geom->grid->ind[comp], p->x, p->y, p->z, sel->n);
}

static void sort_particles(struct plist *p, int (*compare_func) ())
{
	qsort(p->list, p->n, sizeof(*p->list), compare_func);
}

static void sort_particles_for_sending(struct plist *p)
{
	qsort(p->list, p->n, sizeof(*p->list), normal_first);
	size_t m = p->n - 1;
	for (size_t i = 0; i < p->n; ++i)
		if (p->list[i].status == P_NEW_PROC) {
			m = i;
			break;
		}
	qsort(p->list + m, p->n - m, sizeof(*p->list), increasing_proc);
}

int normal_first(const void *a, const void *b)
{
	const size_t a_val = (* (struct particle *) a).status;
	const size_t b_val = (* (struct particle *) b).status;
	return a_val - b_val;
}

int increasing_proc(const void *a, const void *b)
{
	const size_t a_val = (* (struct particle *) a).proc;
	const size_t b_val = (* (struct particle *) b).proc;
	return a_val - b_val;
}

double max_particle_velocity(struct plist *p, struct element **sel)
{
	double umax = 0.0;
	const double t = (*sel)->params->t;
	for (size_t i = 0; i < p->n; ++i) {
		double u = physical_particle_velocity(&p->list[i],
				sel[p->list[i].el], 0);
		double v = physical_particle_velocity(&p->list[i],
				sel[p->list[i].el], 1);
		double w = physical_particle_velocity(&p->list[i],
				sel[p->list[i].el], 2);
		double unorm = sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));
		if (unorm > umax)
			umax = unorm;
	}
	return global_max(sel, umax);
}
