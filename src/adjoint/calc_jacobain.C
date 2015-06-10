/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *******************************************************************
 */
//Jan 30, 2015
//haghakha
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define KEY0   8
#define KEY1   3
#define ITER   5
#define EFFELL 4
#define J      1
#define JACIND 0

#define DEBUG

void reset_resflag(ResFlag resflag[5]);

void calc_jacobian(DualMesh* dualmesh, MatProps* matprops_ptr, TimeProps* timeprops_ptr,
    MapNames *mapname_ptr, double const increment) {
	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int neigh_flag;

	int iter = timeprops_ptr->iter;
	double tiny = GEOFLOW_TINY;

	cout << "computing jacobian for time iteration " << iter << endl;

	//here are some dummy values that we need for calc_edge_state
	int order_flag = 1;
	double outflow = 0;

	//this array holds ResFlag for element itself and its neighbors
	ResFlag resflag[5];
	reset_resflag(resflag);

	dualmesh->calc_flux();

	int Ny = dualmesh->get_Ny();
	int Nx = dualmesh->get_Nx();

	for (int yind = 0; yind < Ny; ++yind)
		for (int xind = 0; xind < Nx; ++xind) {

			if (yind == 0 || xind == 0 || yind == Ny - 1 || xind == Nx - 1)
				// for boundary cells, Jacobian is zero because they are independent from neighbours
				(dualmesh->get_dualcell(yind, xind))->set_jacobian();
			else {
				DualCell* cell = dualmesh->get_dualcell(yind, xind);

				double state_vars[NUM_STATE_VARS], gravity[NUM_STATE_VARS];
				double d_gravity[DIMENSION], curvature[DIMENSION];
				double d_state_vars_old[NUM_STATE_VARS * DIMENSION];
				double *prev_state_vars = cell->get_prev_state_vars();
				//Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
				double bedfrict = matprops_ptr->bedfrict[1];
				double dx[2] = { dualmesh->get_dx(), dualmesh->get_dy() };
				double kactxy[DIMENSION];
				double orgSrcSgn[2];

				for (int ind = 0; ind < NUM_STATE_VARS; ++ind) {

					state_vars[ind] = *(cell->get_state_vars() + ind);
					gravity[ind] = *(cell->get_gravity() + ind);

				}

				for (int ind = 0; ind < DIMENSION; ++ind) {

					d_gravity[ind] = *(cell->get_d_gravity() + ind);
					curvature[ind] = *(cell->get_curvature() + ind);

				}

				for (int ind = 0; ind < NUM_STATE_VARS * DIMENSION; ++ind)
					d_state_vars_old[ind] = *(cell->get_d_state_vars() + ind);

				cell->calc_slopes(dualmesh);

				double fluxold[4][NUM_STATE_VARS];

				record_flux(dualmesh, cell, fluxold);

				if (timeprops_ptr->iter < 51)
					matprops_ptr->frict_tiny = 0.1;
				else
					matprops_ptr->frict_tiny = 0.000000001;

				orgSourceSgn(cell, matprops_ptr->frict_tiny, orgSrcSgn);

				double dt = timeprops_ptr->dt.at(iter - 1); //at final time step we do not need the computation of adjoint and we always compute it for the previouse time so we need iter.
				double dtdx = dt / dx[0];
				double dtdy = dt / dx[1];

				double org_residual[3];

				//here we compute the residuals
				residual(org_residual, state_vars, prev_state_vars, fluxold[0], //4
				    fluxold[1], fluxold[2], fluxold[3], dtdx, dtdy, dt, d_state_vars_old, //7
				    (d_state_vars_old + NUM_STATE_VARS), curvature, //2
				    matprops_ptr->intfrict, //1
				    bedfrict, gravity, d_gravity, cell->get_kact(), //4
				    matprops_ptr->frict_tiny, orgSrcSgn, 0. /*here we set increment0.*/, //3
				    matprops_ptr->epsilon, 1/*here we set srcflag=1*/); //2

				for (int effelement = 0; effelement < 5; effelement++) {

					double void_res[NUM_STATE_VARS] = { 0., 0., 0. };

					if (effelement == 0 && prev_state_vars[0] == 0.)

						//this is a void element so the residual vector does not change by changing it's values
						cell->set_jacobianMat_zero(effelement);

					else if (effelement != 0 && void_neigh_cell(dualmesh, cell, effelement))

						cell->set_jacobianMat_zero(effelement);

					else {

						unsigned key_xp[2] = { yind, xind + 1 }, key_xm[2] = { yind, xind - 1 }, key_yp[2] = {
						    yind + 1, xind }, key_ym[2] = { yind - 1, xind };

						for (int j = 0; j < NUM_STATE_VARS; j++) {

							double vec_res[NUM_STATE_VARS];
							double total_res[NUM_STATE_VARS] = { 0., 0., 0. };

							int scheme = 0;
							for (; scheme < 2; scheme++) {

								int updateflux, srcflag;
								reset_resflag(resflag);

								// here we modify increment to one time compute forward and one time compute backward difference if it is necessary
								double signe = pow(-1., scheme);
								double incr = signe * increment;
								increment_state(dualmesh, cell, incr, effelement, j, &updateflux, &srcflag,
								    resflag);
								int aaa = 0, bbb = 1;
								if (yind == KEY0 && xind == KEY1 && effelement == EFFELL && iter == ITER && j == J)
									aaa = bbb;

								calc_flux_slope_kact(dualmesh, cell, matprops_ptr, effelement, updateflux, resflag);

								double flux[4][NUM_STATE_VARS];

								for (int ind = 0; ind < NUM_STATE_VARS; ++ind) {

									flux[0][ind] = *(cell->get_xflux() + ind);
									flux[1][ind] = *(cell->get_yflux() + ind);

									flux[2][ind] =
									    *((dualmesh->get_dualcell(key_xm[0], key_xm[1]))->get_xflux() + ind);
									flux[3][ind] =
									    *((dualmesh->get_dualcell(key_ym[0], key_ym[1]))->get_yflux() + ind);

								}

#ifdef DEBUG
								int dbgflag = 0, printflag = 0;
								if (yind == KEY0 && xind == KEY1 && effelement == EFFELL && iter == ITER && j == J)

									flux_debug(cell, fluxold, flux, effelement, j, iter, dt, dx);
#endif
								double *d_state_vars = cell->get_d_state_vars();

								//here we compute the residuals
								residual(vec_res, state_vars, prev_state_vars, flux[0], //4
								    flux[1], flux[2], flux[3], dtdx, dtdy, dt, d_state_vars, //7
								    (d_state_vars + NUM_STATE_VARS), curvature, //2
								    matprops_ptr->intfrict, //1
								    bedfrict, gravity, d_gravity, cell->get_kact(), //4
								    matprops_ptr->frict_tiny, orgSrcSgn, incr, //3
								    matprops_ptr->epsilon, srcflag); //2

								if (scheme == 0) //forward difference
									for (int ind = 0; ind < NUM_STATE_VARS; ind++)
										vec_res[ind] = vec_res[ind]-org_residual[ind];
								else
									for (int ind = 0; ind < NUM_STATE_VARS; ind++)
										vec_res[ind] = org_residual[ind]- vec_res[ind];

								//we have to return everything back
								restore(dualmesh, cell, matprops_ptr, effelement, j, increment, fluxold,
								    d_state_vars_old);

								for (int ind = 0; ind < NUM_STATE_VARS; ind++)
									total_res[ind] += signe * vec_res[ind];

//											if (scheme == 0 && fabs(vec_res[0] / incr) < 5.
//													&& fabs(vec_res[1] / incr) < 5.
//													&& fabs(vec_res[2] / incr) < 5.)
//												//this means that forward difference is enough, and we do not need to compute central difference
								break;
							}

							double jacincr = increment; //= scheme == 0 ? increment : 2. * increment;

							cell->set_jacobian(effelement, total_res, j,
							// following term is necessary to consider the scheme that whether it is forward difference or central difference
							    jacincr); //sets the propper components of the Jacobian for this element

						}
					}
				}
			}
		}
}

int jac_mat_index(int effelement, int xp) {

	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;
	// effelement - 1 shows which neighbor we are processing
	if (effelement == 0)
		return 0;
	else if (effelement - 1 == xp)
		return 1;
	else if (effelement - 1 == yp)
		return 2;
	else if (effelement - 1 == xm)
		return 3;
	else if (effelement - 1 == ym)
		return 4;
	else
		exit(1);
}

//this function returns 1 if the neighbor element is void
int void_neigh_cell(DualMesh* dualmesh, DualCell* cell, int effelement) {

	int a, b;
	set_ab(&a, &b, effelement);

	int i = *(cell->get_key()), j = *(cell->get_key() + 1);

	DualCell* neigh_cell = dualmesh->get_dualcell(i + a, j + b); //basically we are checking all neighbor elements, and start from xp neighbor

	assert(neigh_cell);

	if (*(neigh_cell->get_prev_state_vars()) == 0.)
		return 1;

	return 0;
}

//this function returns 1 if the neighbor element is void
int void_neigh_elem(HashTable* El_Table, Element* Curr_El, int effelement) {

	Element* neigh_elem = (Element*) (El_Table->lookup(
	    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));

	assert(neigh_elem);

	if (*(neigh_elem->get_prev_state_vars()) == 0.)
		return 1;

	return 0;
}

void restore(DualMesh* dualmesh, DualCell* cell, MatProps* matprops_ptr, int effelement, int j,
    double increment, double fluxold[4][NUM_STATE_VARS],
    double d_state_vars_old[NUM_STATE_VARS * DIMENSION]) {

	double *prev_state_vars = cell->get_prev_state_vars();

	unsigned* key = cell->get_key();

	unsigned key_xp[2] = { key[0], key[1] + 1 }, key_xm[2] = { key[0], key[1] - 1 }, key_yp[2] = {
	    key[0] + 1, key[1] }, key_ym[2] = { key[0] - 1, key[1] };

	if (effelement == 0) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		prev_state_vars[j] -= increment; //changing the state varibales at the element itself

		for (int ind = 0; ind < NUM_STATE_VARS; ++ind) {

			*(cell->get_xflux() + ind) = fluxold[0][ind];
			*(cell->get_yflux() + ind) = fluxold[1][ind];
			*((dualmesh->get_dualcell(key_xm[0], key_xm[1]))->get_xflux() + ind) = fluxold[2][ind];
			*((dualmesh->get_dualcell(key_ym[0], key_ym[1]))->get_yflux() + ind) = fluxold[3][ind];

		}

	} else {

		int a, b;
		set_ab(&a, &b, effelement);

		DualCell* neigh_cell = dualmesh->get_dualcell(key[0] + a, key[1] + b); //basically we are checking all neighbor elements, and start from xp neighbor

		*(neigh_cell->get_prev_state_vars() + j) -= increment;

		if (effelement == 1)

			for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
				*(cell->get_xflux() + ind) = fluxold[0][ind];

		else if (effelement == 2)

			for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
				*(cell->get_yflux() + ind) = fluxold[1][ind];

		else if (effelement == 3) {

			DualCell* cell_xm = dualmesh->get_dualcell(key_xm[0], key_xm[1]);
			for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
				*(cell_xm->get_xflux() + ind) = fluxold[2][ind];

		} else {

			DualCell* cell_ym = dualmesh->get_dualcell(key_ym[0], key_ym[1]);
			for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
				*(cell_ym->get_yflux() + ind) = fluxold[3][ind];
		}

	}
	cell->calc_slopes(dualmesh);
	return;
}

void calc_flux_slope_kact(DualMesh* dualmesh, DualCell* cell, MatProps* matprops_ptr,
    int effelement, int updateflux, ResFlag resflag[5]) {

	cell->calc_slopes(dualmesh);
	unsigned* key = cell->get_key();

	int lgft[2] = { 0, 0 };

//	if (srcflag && effelement == 0) {
//
//		gmfggetcoef_(cell->get_prev_state_vars(), d_state_vars, (d_state_vars + NUM_STATE_VARS), dx,
//		    &(matprops_ptr->bedfrict[0]), &(matprops_ptr->intfrict), &kactxy[0], &kactxy[1], &tiny,
//		    &(matprops_ptr->epsilon));
//
//		Curr_El->put_kactxy(kactxy);
//		Curr_El->calc_stop_crit(matprops_ptr);
//	}

	if (effelement == 0 && updateflux) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		lgft[0] = resflag[0].lgft;
		lgft[1] = resflag[1].lgft;
		//change of the state_vars causes the all around fluxes change, this update xp ,yp
		cell->update_flux(dualmesh, lgft);
		// earlier in this file, we made sure that this element is not a boundary element
		// so here we do not require to check the neighbor elements to make sure they are not located on the boundary

		DualCell* cell_xm;
		cell_xm = dualmesh->get_dualcell(key[0], key[1] - 1);
		assert(cell_xm);

		lgft[0] = resflag[3].lgft;
		lgft[1] = resflag[0].lgft;

		//this update the flux on share edge with xm
		cell_xm->update_flux_x(dualmesh, lgft);

		DualCell* cell_ym = dualmesh->get_dualcell(key[0] - 1, key[1]);
		assert(cell_ym);

		lgft[0] = resflag[4].lgft;
		lgft[1] = resflag[0].lgft;

		//this update the flux on share edge with ym
		cell_ym->update_flux_y(dualmesh, lgft);

	} else if (effelement != 0 && updateflux) {

		if (effelement == 1) {

			//this is for xp cell

			lgft[0] = resflag[0].lgft;
			lgft[1] = resflag[1].lgft;

			//if we change the state variables in xp or yp, just the flux at this element has to be updated
			cell->update_flux_x(dualmesh, lgft);

		} else if (effelement == 2) {

			//this is for yp cell

			lgft[0] = resflag[0].lgft;
			lgft[1] = resflag[2].lgft;

			//if we change the state variables in xp or yp, just the flux at this element has to be updated
			cell->update_flux_y(dualmesh, lgft);

		} else if (effelement == 3) {

			//this is for xm cell

			DualCell* cell_xm = dualmesh->get_dualcell(key[0], key[1] - 1);
			assert(cell_xm);

			lgft[0] = resflag[3].lgft;
			lgft[1] = resflag[0].lgft;

			//this update the flux on share edge with xm
			cell_xm->update_flux_x(dualmesh, lgft);

		} else {

			DualCell* cell_ym = dualmesh->get_dualcell(key[0] - 1, key[1]);
			assert(cell_ym);

			lgft[0] = resflag[4].lgft;
			lgft[1] = resflag[0].lgft;

			//this update the flux on share edge with ym
			cell_ym->update_flux_y(dualmesh, lgft);
		}

	}

}

void increment_state(DualMesh* dualmesh, DualCell* cell, double increment, int effelement, int j,
    int* updateflux, int* srcflag, ResFlag resflag[5]) {

	*updateflux = 1;
	*srcflag = 1;
	double *prev_state_vars = cell->get_prev_state_vars();

	if (effelement == 0) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		if (j == 0 && prev_state_vars[j] < GEOFLOW_TINY) {
			*updateflux = 0;
			*srcflag = 0;
			resflag[0].lgft = 1;
		}

		prev_state_vars[j] += increment; //changing the state varibales at the element itself

	} else {

		int a, b;
		set_ab(&a, &b, effelement);

		int indy = *(cell->get_key()), indx = *(cell->get_key() + 1);

		DualCell* neigh_cell = dualmesh->get_dualcell(indy + a, indx + b); //basically we are checking all neighbor elements, and start from xp neighbor

		assert(neigh_cell);

		if (j == 0 && *(neigh_cell->get_prev_state_vars() + j) < GEOFLOW_TINY) {
			*updateflux = 0;
			resflag[effelement].lgft = 1;
		}

		*(neigh_cell->get_prev_state_vars() + j) += increment;

	}
	return;
}

void reset_resflag(ResFlag resflag[5]) {

	for (int i = 0; i < 5; i++) {
		resflag[i].callflag = 1;
		resflag[i].lgft = 0;
	}

	return;
}
