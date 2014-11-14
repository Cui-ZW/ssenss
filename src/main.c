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

#include <mpi.h>
#include "structs.h"
#include "basis.h"
#include "poisson.h"
#include "navier_stokes.h"
#include "init.h"
#include "final.h"
#include "input.h"
#include "dataio.h"
#include "test.h"
#include "parse.h"
#include "interpolation.h"
#include "geometry.h"
#include "deriv.h"
#include "precon.h"

int main(int argc, char **argv)
{
	/* Initialise MPI and communicator
	 */
	int rank, nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	struct communicator *comm = new_communicator(rank, nproc,
			MPI_COMM_WORLD, MPI_INFO_NULL);

	/* Initialise and read parameters
	 */
	struct parameters *params = new_parameters(rank, nproc);
	get_parameters(params);
	parse_command_arguments(argc, argv, params);
	params->nconv = (3 * params->n) / 2;
	params->svvn = 2 * params->n / 3;
	params->svva = 0.0;
	get_grid_size(params, comm, params->gridfile);

	/* Initialise basis
	 */
	struct basis *gll_basis = new_gll_basis(params->n);
	struct basis *conv_basis = new_gll_basis(params->nconv);
	create_interpolation_matrices(gll_basis, conv_basis);
	compute_svv_derivative_matrix(gll_basis, params->svvn, params->svva);

	/* Initialise and read geometry
	 */
	struct geometry **geom = new_geometries(gll_basis, params->nel,
			params->neltot);
	struct geometry **conv_geom = new_geometries(conv_basis, params->nel,
			params->neltot);
	read_grid(geom, comm, gll_basis, params->gridfile);
	read_grid(conv_geom, comm, conv_basis, params->gridfile);

	/* Initialise elements
	 */
	struct element **sel = new_elements(geom, conv_geom, gll_basis,
			conv_basis, params, comm);
	connect_elements(sel);
	initialise_preconditioners(sel);

	initialise_from_file(sel, params->initfile);	
	if (params->test)
		initialise_test(sel);

	/* Write initial condition to file
	 */
	char filename[200];
	sprintf(filename, "%s/vel_%08d.h5", params->datdir, 0);
	write_to_file(sel, filename, WRITE_VELOCITY);

	/* Run simulation until tmax
	 */
	int iter = 0;
	while (params->t < params->tmax - 0.5 * params->dt) {	
		++iter;
		params->wtime[0] = MPI_Wtime();
		set_source_term(sel, X_SOURCE, one_func, params->gx);
		set_source_term(sel, Y_SOURCE, one_func, params->gy);
		set_source_term(sel, Z_SOURCE, one_func, params->gz);
		integrate_navier_stokes(sel);
		params->wtime[1] = MPI_Wtime();
		params->wtime[2] = MPI_Wtime();
		write_log(sel);
		if ((iter % params->write_every) == 0) {
			sprintf(filename, "%s/vel_%08d.h5", params->datdir,
					iter);
			write_to_file(sel, filename, WRITE_VELOCITY);
		}
	}

	/* Free allocated resources and stop mpi
	 */
	free_basis(gll_basis);
	free_basis(conv_basis);
	disconnect_elements(sel);
	free_elements(sel);
	free_communicator(comm);
	free_geometries(geom);
	free_geometries(conv_geom);
	free_parameters(params);

	stop_mpi();
	return 0;
}
