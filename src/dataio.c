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
#include <stdio.h>
#include <mpi.h>
#include <hdf5.h>
#include "structs.h"
#include "utils.h"
#include "linalg.h"
#include "tripalg.h"
#include "basis.h"
#include "interpolation.h"
#include "geometry.h"
#include "poisson.h"
#include "test.h"
#include "dataio.h"

#define INC_PTR(A, B) ((void *) ((char *) (A) + (B)))

/* Velocity buffer definition
 */
#define UTSIZE sizeof(double)
#define VEL_INC(A, B) (vel_const_inc[A] + vel_var_inc[A] * (B))
#define VEL_SIZE(A) (3 * UTSIZE * (A))
enum { U, V, W };
const static int vel_const_inc[3] = { [U] = 0, [V] = 0, [W] = 0 };
const static int vel_var_inc[3] = { [U] = 0, [V] = UTSIZE, [W] = 2 * UTSIZE };

/* Grid buffer definition
 */
#define XTSIZE sizeof(double)
#define NTSIZE sizeof(int)
#define GRID_INC(A, B) (grid_const_inc[A] + grid_var_inc[A] * (B))
#define GRID_SIZE(A) (3 * XTSIZE * (A) + 52 * NTSIZE)

enum { NEIGHBOURS, NEIGH_ROT, EDGE_SIZE, EDGE_NEIGHBOURS, CORNER_SIZE,
		CORNER_NEIGHBOURS, X, Y, Z };
const static int grid_const_inc[9] = { [NEIGHBOURS] = 0,
		[NEIGH_ROT] = 6 * NTSIZE, [EDGE_SIZE] = 12 * NTSIZE,
		[EDGE_NEIGHBOURS] = 24 * NTSIZE, [CORNER_SIZE] = 36 * NTSIZE,
		[CORNER_NEIGHBOURS] = 44 * NTSIZE, [X] = 52 * NTSIZE,
		[Y] = 52 * NTSIZE, [Z] = 52 * NTSIZE };
const static int grid_var_inc[9] = { [NEIGHBOURS] = 0,
		[NEIGH_ROT] = 0, [EDGE_SIZE] = 0, [EDGE_NEIGHBOURS] = 0,
		[CORNER_SIZE] = 0, [CORNER_NEIGHBOURS] = 0, [X] = 0,
		[Y] = XTSIZE, [Z] = 2 * XTSIZE };

struct hdf5_data {
	hid_t type;
	void *data;
};

struct file_params {
	int nel;
	int n;
};

static hid_t open_hdf5_file(const char *filename, MPI_Comm comm, MPI_Info info,
		int rw);
static herr_t close_hdf5_file(hid_t file_id);
static struct hdf5_data new_grid_data(int n, int nel);
static struct hdf5_data new_velocity_data(int n, int nel);
static void free_hdf5_data(struct hdf5_data h5data);
static void copy_velocity_to_struct(struct element **sel,
		struct hdf5_data h5data);
static void copy_grid_to_struct(struct element **sel, struct hdf5_data h5data);
static void copy_velocity_from_struct(struct element **sel,
		struct hdf5_data h5data, struct file_params fparams);
static void copy_grid_from_struct(struct geometry **geom,
		struct communicator *comm, struct basis *basis,
		struct hdf5_data h5data, struct file_params fparams);
static herr_t read_write_parameters(hid_t file_id, struct file_params *fparams,
		int rw);
static herr_t read_write_element_array(hid_t file_id, struct hdf5_data h5data,
		hsize_t *dimsf, hsize_t *count, hsize_t *offset, int rw);

void write_log(struct element **sel)
{
	const double div = max_divergence(sel);
	const double umax = maxmin_val(sel, X_VELOCITY, 1.0);
	const double vmax = maxmin_val(sel, Y_VELOCITY, 1.0);
	const double wmax = maxmin_val(sel, Z_VELOCITY, 1.0);
	const double pmax = maxmin_val(sel, PRESSURE, 1.0);
	const double dt1 = (*sel)->params->wtime[1] - (*sel)->params->wtime[0];
	const double dt2 = (*sel)->params->wtime[2] - (*sel)->params->wtime[1];
	static int iter = 0;
	char hfmt[100] = "%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t%-15s\n";
	char vfmt[100] = "%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n";

	if ((*sel)->comm->rank == 0) {
		if ((iter++ % 20) == 0)
			printf(hfmt, "Time", "Max div", "Umax", "Vmax", "Wmax",
					"Pmax", "Wtime 1", "Wtime 2");
		printf(vfmt, (*sel)->params->t, div, umax, vmax, wmax, pmax,
				dt1, dt2);
	}
}

void write_to_file(struct element **sel, const char *filename, int write_data)
{
	int rank = (*sel)->params->rank;
	int nel = (*sel)->params->nel;
	int neltot = (*sel)->params->neltot;
	int nproc = (*sel)->params->nproc;

	hsize_t dimsf[1] = {neltot};
	hsize_t count[1] = {nel};
	hsize_t offset[1] = {id_offset(rank, neltot, nproc)};

	hid_t file_id = open_hdf5_file(filename, (*sel)->comm->mpi_comm,
			(*sel)->comm->mpi_info, write_data);

	struct hdf5_data h5data;
	if (write_data == WRITE_GRID) {
		h5data = new_grid_data((*sel)->params->n, nel);
		copy_grid_to_struct(sel, h5data);
	} else if (write_data == WRITE_VELOCITY) {
		h5data = new_velocity_data((*sel)->params->n, nel);
		copy_velocity_to_struct(sel, h5data);
	}

	struct file_params fparams;
	fparams.n = (*sel)->params->n;
	fparams.nel = (*sel)->params->neltot;

	herr_t status = read_write_parameters(file_id, &fparams, write_data);
	status = read_write_element_array(file_id, h5data, dimsf, count,
			offset, write_data);

	free_hdf5_data(h5data);
	close_hdf5_file(file_id);
}

void read_velocity(struct element **sel, const char *filename)
{
	int rank = (*sel)->params->rank;
	int nel = (*sel)->params->nel;
	int neltot = (*sel)->params->neltot;
	int nproc = (*sel)->params->nproc;

	hsize_t dimsf[1] = {neltot};
	hsize_t count[1] = {nel};
	hsize_t offset[1] = {id_offset(rank, neltot, nproc)};

	hid_t file_id = open_hdf5_file(filename, (*sel)->comm->mpi_comm,
			(*sel)->comm->mpi_info, READ_VELOCITY);

	struct file_params fparams;
	herr_t status = read_write_parameters(file_id, &fparams, READ_VELOCITY);
	if (fparams.nel != (*sel)->params->neltot)
		fatal_error("in read_velocity: neltot in file is wrong");

	struct hdf5_data h5data;
	h5data = new_velocity_data(fparams.n, fparams.nel);

	status = read_write_element_array(file_id, h5data, dimsf, count,
			offset, READ_VELOCITY);
	copy_velocity_from_struct(sel, h5data, fparams);

	free_hdf5_data(h5data);
	close_hdf5_file(file_id);
}

void get_grid_size(struct parameters *params, struct communicator *comm,
		const char *filename)
{
	hid_t file_id = open_hdf5_file(filename, comm->mpi_comm,
			comm->mpi_info, READ_GRID);

	struct file_params fparams;
	herr_t status = read_write_parameters(file_id, &fparams, READ_GRID);
	close_hdf5_file(file_id);

	if (fparams.nel < comm->nproc)
		fatal_error("nproc > neltot");
	int nel = local_nel(comm->rank, fparams.nel, comm->nproc);
	params->neltot = fparams.nel;
	params->nel = nel;
}

void read_grid(struct geometry **geom, struct communicator *comm,
		struct basis *basis, const char *filename)
{
	int rank = comm->rank;
	int nel = (*geom)->nel;
	int neltot = (*geom)->neltot;
	int nproc = comm->nproc;

	hsize_t dimsf[1] = {neltot};
	hsize_t count[1] = {nel};
	hsize_t offset[1] = {id_offset(rank, neltot, nproc)};

	hid_t file_id = open_hdf5_file(filename, comm->mpi_comm,
			comm->mpi_info, READ_GRID);

	struct file_params fparams;
	herr_t status = read_write_parameters(file_id, &fparams, READ_GRID);

	struct hdf5_data h5data;
	h5data = new_grid_data(fparams.n, fparams.nel);

	status = read_write_element_array(file_id, h5data, dimsf, count,
			offset, READ_GRID);
	copy_grid_from_struct(geom, comm, basis, h5data, fparams);

	free_hdf5_data(h5data);
	close_hdf5_file(file_id);

	for (int i = 0; i < nel; ++i)
		compute_geometry_transformations(geom[i], basis);
}

static hid_t open_hdf5_file(const char *filename, MPI_Comm comm, MPI_Info info,
		int rw)
{
	H5E_auto2_t old_func;
	void *old_client_data;
	H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, comm, info);

	hid_t file_id;
	if (rw & (READ_GRID | READ_VELOCITY | READ_PARTICLES))
		file_id = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
	else if (rw & (WRITE_GRID | WRITE_VELOCITY | WRITE_PARTICLES))
		file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT,
				plist_id);
	H5Pclose(plist_id);

	if (file_id < 0) {
		char tmp[100];
		sprintf(tmp, "in open_hdf5_file: unable to open file '%s'",
				filename);
		fatal_error(tmp);
	}
	H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
	return file_id;
}

static herr_t close_hdf5_file(hid_t file_id)
{
	return H5Fclose(file_id);
}

static struct hdf5_data new_grid_data(int n, int nel)
{
	size_t ntot = (size_t) n * (size_t) n * (size_t) n;
	hsize_t array_dim[] = {ntot};
	void *tmp = ec_malloc(GRID_SIZE(ntot) * nel);
	hsize_t intdim6[] = {6};
	hsize_t intdim8[] = {8};
	hsize_t intdim12[] = {12};

	hid_t array_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, array_dim);
	hid_t intarr6 = H5Tarray_create(H5T_NATIVE_INT, 1, intdim6);
	hid_t intarr8 = H5Tarray_create(H5T_NATIVE_INT, 1, intdim8);
	hid_t intarr12 = H5Tarray_create(H5T_NATIVE_INT, 1, intdim12);
	hid_t h5struct = H5Tcreate(H5T_COMPOUND, GRID_SIZE(ntot));

	H5Tinsert(h5struct, "Neighbours", GRID_INC(NEIGHBOURS, ntot), intarr6);
	H5Tinsert(h5struct, "Neighbour rotation", GRID_INC(NEIGH_ROT, ntot),
			intarr6);
	H5Tinsert(h5struct, "Corner size", GRID_INC(CORNER_SIZE, ntot),
			intarr8);
	H5Tinsert(h5struct, "Corner neighbours", GRID_INC(CORNER_NEIGHBOURS,
			ntot), intarr8);
	H5Tinsert(h5struct, "Edge size", GRID_INC(EDGE_SIZE, ntot), intarr12);
	H5Tinsert(h5struct, "Edge neighbours", GRID_INC(EDGE_NEIGHBOURS, ntot),
			intarr12);
	H5Tinsert(h5struct, "X", GRID_INC(X, ntot), array_tid);
	H5Tinsert(h5struct, "Y", GRID_INC(Y, ntot), array_tid);
	H5Tinsert(h5struct, "Z", GRID_INC(Z, ntot), array_tid);

	struct hdf5_data h5data;
	h5data.type = h5struct;
	h5data.data = tmp;

	H5Tclose(array_tid);
	H5Tclose(intarr6);
	H5Tclose(intarr8);
	H5Tclose(intarr12);

	return h5data;
}

static struct hdf5_data new_velocity_data(int n, int nel)
{
	size_t ntot = (size_t) n * (size_t) n * (size_t) n;
	hsize_t array_dim[] = {ntot};
	void *tmp = ec_malloc(VEL_SIZE(ntot) * nel);

	hid_t array_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, array_dim);
	hid_t h5struct = H5Tcreate(H5T_COMPOUND, VEL_SIZE(ntot));
	H5Tinsert(h5struct, "U", VEL_INC(U, ntot), array_tid);
	H5Tinsert(h5struct, "V", VEL_INC(V, ntot), array_tid);
	H5Tinsert(h5struct, "W", VEL_INC(W, ntot), array_tid);

	struct hdf5_data h5data;
	h5data.type = h5struct;
	h5data.data = tmp;

	H5Tclose(array_tid);

	return h5data;
}

static void free_hdf5_data(struct hdf5_data h5data)
{
	H5Tclose(h5data.type);
	free(h5data.data);
}

static void copy_velocity_to_struct(struct element **sel,
		struct hdf5_data h5data)
{
	size_t ntot = (*sel)->ntot;
	void *u = h5data.data;

	for (int i = 0; i < (*sel)->nel; ++i) {
		copy_array(**sel[i]->u->x, INC_PTR(u, VEL_INC(U, ntot)
				+ i * VEL_SIZE(ntot)), ntot);
		copy_array(**sel[i]->u->y, INC_PTR(u, VEL_INC(V, ntot)
				+ i * VEL_SIZE(ntot)), ntot);
		copy_array(**sel[i]->u->z, INC_PTR(u, VEL_INC(W, ntot)
				+ i * VEL_SIZE(ntot)), ntot);
	}
}

static void copy_grid_to_struct(struct element **sel, struct hdf5_data h5data)
{
	void *x = h5data.data;
	size_t ntot = (*sel)->ntot;
	int neltot = (*sel)->params->neltot;
	int nproc = (*sel)->params->nproc;

	for (int i = 0; i < (*sel)->nel; ++i) {
		const size_t j = i * GRID_SIZE(ntot);
		copy_int_array(sel[i]->geom->neigh_rot, INC_PTR(x,
				GRID_INC(NEIGH_ROT, ntot) + j), 6);
		copy_int_array(sel[i]->geom->edge_size, INC_PTR(x,
				GRID_INC(EDGE_SIZE, ntot) + j), 12);
		copy_int_array(sel[i]->geom->corner_size, INC_PTR(x,
				GRID_INC(CORNER_SIZE, ntot) + j), 8);
		get_neighbour_ids(sel[i]->geom->neigh_el, sel[i]->geom
				->neigh_proc, INC_PTR(x, GRID_INC(NEIGHBOURS,
				ntot) + j), neltot, nproc, 6);
		get_neighbour_ids(sel[i]->geom->neigh_edge_el, sel[i]->geom
				->neigh_edge_proc, INC_PTR(x,
				GRID_INC(EDGE_NEIGHBOURS, ntot) + j), neltot,
				nproc, 12);
		get_neighbour_ids(sel[i]->geom->neigh_corn_el, sel[i]->geom
				->neigh_corn_proc, INC_PTR(x,
				GRID_INC(CORNER_NEIGHBOURS, ntot) + j), neltot,
				nproc, 8);
		copy_array(**sel[i]->geom->grid->x, INC_PTR(x, GRID_INC(X, ntot)
				+ j), (*sel)->ntot);
		copy_array(**sel[i]->geom->grid->y, INC_PTR(x, GRID_INC(Y, ntot)
				+ j), (*sel)->ntot);
		copy_array(**sel[i]->geom->grid->z, INC_PTR(x, GRID_INC(Z, ntot)
				+ j), (*sel)->ntot);
	}
}

static void copy_velocity_from_struct(struct element **sel,
		struct hdf5_data h5data, struct file_params fparams)
{
	size_t ntot = (size_t) fparams.n * (size_t) fparams.n
			* (size_t) fparams.n;
	void *u = h5data.data;
	struct basis *fbasis = new_gll_basis(fparams.n);
	struct basis *sbasis = new_gll_basis((*sel)->params->n);
	create_interpolation_matrices(fbasis, sbasis);
	double ***tmp = new_3d_array(fparams.n, fparams.n, fparams.n);

	for (int i = 0; i < (*sel)->nel; ++i) {
		const size_t j = i * VEL_SIZE(ntot);
		copy_array(INC_PTR(u, VEL_INC(U, ntot) + j), **tmp, ntot);
		interpolate_to_grid(sbasis, sel[i]->u->x, fbasis, tmp);
		copy_array(INC_PTR(u, VEL_INC(V, ntot) + j), **tmp, ntot);
		interpolate_to_grid(sbasis, sel[i]->u->y, fbasis, tmp);
		copy_array(INC_PTR(u, VEL_INC(W, ntot) + j), **tmp, ntot);
		interpolate_to_grid(sbasis, sel[i]->u->z, fbasis, tmp);
	}
	free_3d_array(tmp);
	free_basis(fbasis);
	free_basis(sbasis);
}

static void copy_grid_from_struct(struct geometry **geom,
		struct communicator *comm, struct basis *basis,
		struct hdf5_data h5data, struct file_params fparams)
{
	size_t ntot = (size_t) fparams.n * (size_t) fparams.n
			* (size_t) fparams.n;
	void *x = h5data.data;
	struct basis *fbasis = new_gll_basis(fparams.n);
	struct basis *sbasis = new_gll_basis(basis->n);
	create_interpolation_matrices(fbasis, sbasis);
	double ***tmp = new_3d_array(fparams.n, fparams.n, fparams.n);
	int neltot = (*geom)->neltot;
	int nproc = comm->nproc;

	for (int i = 0; i < (*geom)->nel; ++i) {
		const size_t j = i * GRID_SIZE(ntot);
		geom[i]->id = i;
		copy_int_array(INC_PTR(x, GRID_INC(NEIGH_ROT, ntot) + j),
				geom[i]->neigh_rot, 6);
		copy_int_array(INC_PTR(x, GRID_INC(EDGE_SIZE, ntot) + j),
				geom[i]->edge_size, 12);
		copy_int_array(INC_PTR(x, GRID_INC(CORNER_SIZE, ntot) + j),
				geom[i]->corner_size, 8);
		get_neighbour_elnums_procs(geom[i]->neigh_el, geom[i]
				->neigh_proc, INC_PTR(x, GRID_INC(NEIGHBOURS,
				ntot) + j), neltot, nproc, 6);
		get_neighbour_elnums_procs(geom[i]->neigh_edge_el, geom[i]
				->neigh_edge_proc, INC_PTR(x,
				GRID_INC(EDGE_NEIGHBOURS, ntot) + j), neltot,
				nproc, 12);
		get_neighbour_elnums_procs(geom[i]->neigh_corn_el, geom[i]
				->neigh_corn_proc, INC_PTR(x,
				GRID_INC(CORNER_NEIGHBOURS, ntot) + j), neltot,
				nproc, 8);
		set_boundary_conditions(geom[i]);
		copy_array(INC_PTR(x, GRID_INC(X, ntot) + j), **tmp, ntot);
		interpolate_to_grid(sbasis, geom[i]->grid->x, fbasis, tmp);
		copy_array(INC_PTR(x, GRID_INC(Y, ntot) + j), **tmp, ntot);
		interpolate_to_grid(sbasis, geom[i]->grid->y, fbasis, tmp);
		copy_array(INC_PTR(x, GRID_INC(Z, ntot) + j), **tmp, ntot);
		interpolate_to_grid(sbasis, geom[i]->grid->z, fbasis, tmp);
	}
	free_3d_array(tmp);
	free_basis(fbasis);
	free_basis(sbasis);
}

static herr_t read_write_parameters(hid_t file_id, struct file_params *fparams,
		int rw)
{
	hsize_t dimsf[1] = {1};
	hid_t param_type = H5Tcreate(H5T_COMPOUND, sizeof(*fparams));
	H5Tinsert(param_type, "n", HOFFSET(struct file_params, n),
			H5T_NATIVE_INT);
	H5Tinsert(param_type, "nel", HOFFSET(struct file_params, nel),
			H5T_NATIVE_INT);

	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	herr_t status;
	if (rw & (READ_GRID | READ_VELOCITY | READ_PARTICLES)) {
		hid_t dset_id = H5Dopen(file_id, "Parameters", H5P_DEFAULT);
		status = H5Dread(dset_id, param_type, H5S_ALL, H5S_ALL,
				plist_id, fparams);
		H5Dclose(dset_id);
	} else if (rw & (WRITE_GRID | WRITE_VELOCITY | WRITE_PARTICLES)) {
		hid_t filespace = H5Screate_simple(1, dimsf, NULL);
		hid_t dset_id = H5Dcreate(file_id, "Parameters", param_type,
				filespace, H5P_DEFAULT, H5P_DEFAULT,
				H5P_DEFAULT);
		status = H5Dwrite(dset_id, param_type, H5S_ALL, H5S_ALL,
				plist_id, fparams);
		H5Dclose(dset_id);
		H5Sclose(filespace);
	}

	H5Tclose(param_type);
	H5Pclose(plist_id);
	return status;
}

static herr_t read_write_element_array(hid_t file_id, struct hdf5_data h5data,
		hsize_t *dimsf, hsize_t *count, hsize_t *offset, int rw)
{
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	
	herr_t status;
	if (rw & (READ_GRID | READ_VELOCITY | READ_PARTICLES)) {
		hid_t dset_id = H5Dopen(file_id, "Element array", H5P_DEFAULT);
		hid_t memspace = H5Screate_simple(1, count, NULL);
		hid_t filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
				count, NULL);

		status = H5Dread(dset_id, h5data.type, memspace, filespace,
				plist_id, h5data.data);
		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Sclose(memspace);
	} else if (rw & (WRITE_GRID | WRITE_VELOCITY | WRITE_PARTICLES)) {
		hid_t filespace = H5Screate_simple(1, dimsf, NULL);
		hid_t dset_id = H5Dcreate(file_id, "Element array", h5data.type,
				filespace, H5P_DEFAULT,
				H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(filespace);

		hid_t memspace = H5Screate_simple(1, count, NULL);
		filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
				count, NULL);

		status = H5Dwrite(dset_id, h5data.type, memspace, filespace,
				plist_id, h5data.data);
		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Sclose(memspace);
	}

	H5Pclose(plist_id);
	return status;
}
