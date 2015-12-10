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

#define KEY0   1612867174
#define KEY1   1717986918
#define ITER   249
#define EFFELL 0
#define J      0
#define JACIND 0

#define DEBUG

void reset_resflag(ResFlag resflag[EFF_ELL]);

int check_restore(double* prev_state_vars, double* d_state_vars, double* prev_state_vars_old,
    double* d_state_vars_old, double flux[][NUM_STATE_VARS], double fluxold[][NUM_STATE_VARS]);

void set_fluxes_hsens(Element* Curr_El, const Mat3x3*& jac_flux_n_x, const Mat3x3*& jac_flux_p_x,
    const Mat3x3*& jac_flux_n_y, const Mat3x3*& jac_flux_p_y, double* dh_sens, const int& effelement);

const Mat3x3 ZERO_MATRIX;

void calc_jacobian(MeshCTX* meshctx, PropCTX* propctx, PertElemInfo* eleminfo) {

	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	int neigh_flag;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;

	int iter = timeprops_ptr->iter;
	double tiny = GEOFLOW_TINY;

#ifdef DEBUGFILE
	ofstream myfile;
	myfile.open("debug.txt",ios::app);
#endif

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					int boundary = 0;

					//this part handles if the Curr_El is a boundary element
					//cout<<"I am running do not worry"<<endl;
					for (int neighnum = 0; neighnum < 4; neighnum++)
						if (*(Curr_El->get_neigh_proc() + neighnum) == INIT) {
							boundary = 1;
							// this for the element that are on the boundary
							Curr_El->set_jacobian();
							break;
						}
					if (!boundary) {

						double *prev_state_vars = Curr_El->get_prev_state_vars();
						double *state_vars = Curr_El->get_state_vars();
						double *d_state_vars = Curr_El->get_d_state_vars();
						double *gravity = Curr_El->get_gravity();
						double *d_gravity = Curr_El->get_d_gravity();
						double *curvature = Curr_El->get_curvature();
						Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
						double bedfrict = Curr_El->get_effect_bedfrict();
						double *dx = Curr_El->get_dx();

						if (timeprops_ptr->iter < 51)
							matprops_ptr->frict_tiny = 0.1;
						else
							matprops_ptr->frict_tiny = 0.000000001;

						double flux[4][NUM_STATE_VARS];
						get_flux(El_Table, NodeTable, Curr_El->pass_key(), matprops_ptr, myid, flux);

						double dt = timeprops_ptr->dt.at(iter - 1);	//at final time step we do not need the computation of adjoint and we always compute it for the previouse time so we need iter.
						double dtdx = dt / dx[0];
						double dtdy = dt / dx[1];

						int stop[DIMENSION] = { 0, 0 };
						double orgSrcSgn[DIMENSION] = { 0., 0. };

						update_states(state_vars, prev_state_vars, flux[0], flux[1], flux[2], flux[3], dtdx,
						    dtdy, dt, d_state_vars, (d_state_vars + NUM_STATE_VARS), curvature,
						    matprops_ptr->intfrict, bedfrict, gravity, d_gravity, *(Curr_El->get_kactxy()),
						    matprops_ptr->frict_tiny, stop, orgSrcSgn);

						Vec_Mat<9>& jacobian = Curr_El->get_jacobian();

						for (int side = 0; side < 4; side++)
							if (*(Curr_El->get_neigh_proc() + side) == INIT) // this is a boundary!
								for (int ind = 0; ind < NUM_STATE_VARS; ind++)
									state_vars[ind] = 0.;

						for (int effelement = 0; effelement < EFF_ELL; effelement++) { //0 for the element itself, and the rest id for neighbour elements

#ifdef DEBUG
							char filename[] = "jacobian";

							if (/*Curr_El->get_ithelem() == 4*/ *(Curr_El->pass_key()) == KEY0
							 && *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER
							 && effelement == EFFELL)
								Curr_El->write_elem_info(NodeTable, filename, timeprops_ptr->iter, dt);
#endif

							if ((effelement == 0 && prev_state_vars[0] == 0.) || //this is a void element so the residual vector does not change by changing it's values
							    (effelement > 4 && *(Curr_El->get_neigh_proc() + (effelement - 1)) == -2) || //one neighbor in this side
							    (effelement != 0 && void_neigh_elem(El_Table, Curr_El, effelement)) || //this is a void neighbor element so the residual of the curr_el does not depend on this neighbor
							    (effelement > 0 && *(Curr_El->get_neigh_proc() + (effelement - 1)) == INIT))

								Curr_El->set_jacobianMat_zero(effelement);

							else {
								const Mat3x3 *jac_flux_n_x, *jac_flux_p_x, *jac_flux_n_y, *jac_flux_p_y;
								double dh_sens[2];

								set_fluxes_hsens(Curr_El, jac_flux_n_x, jac_flux_p_x, jac_flux_n_y, jac_flux_p_y,
								    dh_sens, effelement);

								calc_jacobian_elem(jacobian(effelement), *jac_flux_n_x, *jac_flux_p_x,
								    *jac_flux_n_y, *jac_flux_p_y, prev_state_vars, d_state_vars,
								    (d_state_vars + NUM_STATE_VARS), curvature, gravity, d_gravity, dh_sens,
								    matprops_ptr->intfrict, bedfrict, *(Curr_El->get_kactxy()), effelement, dtdx,
								    dtdy, dt, stop, orgSrcSgn);

							}
						}
					}
				}
				currentPtr = currentPtr->next;
			}
		}
	}

#ifdef DEBUGFILE
	myfile.close();
#endif

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

int check_restore(double* prev_state_vars, double* d_state_vars, double* prev_state_vars_old,
    double* d_state_vars_old, double flux[][NUM_STATE_VARS], double fluxold[][NUM_STATE_VARS]) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		if (prev_state_vars[i] != prev_state_vars_old[i])
			return 1;

	for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
		if (d_state_vars[i] != d_state_vars_old[i])
			return 1;

	for (int j = 0; j < 4; ++j)
		for (int i = 0; i < NUM_STATE_VARS; ++i)
			if (flux[j][i] != fluxold[j][i])
				return 1;

	return 0;
}

void set_fluxes_hsens(Element* Curr_El, const Mat3x3 *&jac_flux_n_x, const Mat3x3 *&jac_flux_p_x,
    const Mat3x3 *&jac_flux_n_y, const Mat3x3 *&jac_flux_p_y, double* dh_sens, const int& effelement) {

	// (side:x=0,y=1)(direction:neg=0,pos=1)
	FluxJac& flux_jac = (Curr_El->get_flx_jac_cont());
	Matrix<double, 2, 5>& h_slope_sens = Curr_El->get_hslope_sens();

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	if (effelement == 0) {
		int index = 0;

		jac_flux_n_x = &flux_jac(0, 0, index);
		jac_flux_p_x = &flux_jac(0, 1, index);
		jac_flux_n_y = &flux_jac(1, 0, index);
		jac_flux_p_y = &flux_jac(1, 1, index);

		dh_sens[0] = h_slope_sens(0, 0);
		dh_sens[1] = h_slope_sens(1, 0);

	} else if (effelement <= 4) {

		int index = 1;

		if ((effelement - 1) % 4 == xp) {

			dh_sens[0] = h_slope_sens(0, 2);
			dh_sens[1] = 0.;

			jac_flux_n_x = &ZERO_MATRIX;
			jac_flux_p_x = &flux_jac(0, 1, index);
			jac_flux_n_y = &ZERO_MATRIX;
			jac_flux_p_y = &ZERO_MATRIX;

		} else if ((effelement - 1) % 4 == xm) {

			dh_sens[0] = h_slope_sens(0, 1);
			dh_sens[1] = 0.;

			jac_flux_n_x = &flux_jac(0, 0, index);
			jac_flux_p_x = &ZERO_MATRIX;
			jac_flux_n_y = &ZERO_MATRIX;
			jac_flux_p_y = &ZERO_MATRIX;

		} else if ((effelement - 1) % 4 == yp) {

			dh_sens[0] = 0.;
			dh_sens[1] = h_slope_sens(1, 2);

			jac_flux_n_x = &ZERO_MATRIX;
			jac_flux_p_x = &ZERO_MATRIX;
			jac_flux_n_y = &ZERO_MATRIX;
			jac_flux_p_y = &flux_jac(1, 1, index);

		} else {

			dh_sens[0] = 0.;
			dh_sens[1] = h_slope_sens(1, 1);

			jac_flux_n_x = &ZERO_MATRIX;
			jac_flux_p_x = &ZERO_MATRIX;
			jac_flux_n_y = &flux_jac(1, 0, index);
			jac_flux_p_y = &ZERO_MATRIX;
		}

	} else {

		int index = 2;

		if ((effelement - 1) % 4 == xp) {

			dh_sens[0] = h_slope_sens(0, 4);
			dh_sens[1] = 0.;

			jac_flux_n_x = &ZERO_MATRIX;
			jac_flux_p_x = &flux_jac(0, 1, index);
			jac_flux_n_y = &ZERO_MATRIX;
			jac_flux_p_y = &ZERO_MATRIX;

		} else if ((effelement - 1) % 4 == xm) {

			dh_sens[0] = h_slope_sens(0, 3);
			dh_sens[1] = 0.;

			jac_flux_n_x = &flux_jac(0, 0, index);
			jac_flux_p_x = &ZERO_MATRIX;
			jac_flux_n_y = &ZERO_MATRIX;
			jac_flux_p_y = &ZERO_MATRIX;

		} else if ((effelement - 1) % 4 == yp) {

			dh_sens[0] = 0.;
			dh_sens[1] = h_slope_sens(1, 4);

			jac_flux_n_x = &ZERO_MATRIX;
			jac_flux_p_x = &ZERO_MATRIX;
			jac_flux_n_y = &ZERO_MATRIX;
			jac_flux_p_y = &flux_jac(1, 1, index);

		} else {

			dh_sens[0] = 0.;
			dh_sens[1] = h_slope_sens(1, 3);

			jac_flux_n_x = &ZERO_MATRIX;
			jac_flux_p_x = &ZERO_MATRIX;
			jac_flux_n_y = &flux_jac(1, 0, index);
			jac_flux_p_y = &ZERO_MATRIX;
		}

	}

}

