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

#define KEY0   3788876458
#define KEY1   2863311530
#define ITER   5

void error_compute(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	int myid = propctx->myid, numprocs = propctx->numproc;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	ErrorElem* Curr_El = NULL;

	int iter = propctx->timeprops->iter;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (ErrorElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					double *state_vars = Curr_El->get_bilin_state();
					double *prev_state_vars = Curr_El->get_bilin_prev_state();
					double *gravity = Curr_El->get_gravity();
					double *d_gravity = Curr_El->get_d_gravity();
					double *curvature = Curr_El->get_curvature();
					Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
					double bedfrict = Curr_El->get_effect_bedfrict();
					double *dx = Curr_El->get_dx();
					double orgSrcSgn[4];
					double *vec_res = Curr_El->get_residual();
					double *el_error = Curr_El->get_el_error();
					double *bilin_adj = Curr_El->get_bilin_adj();
					double *adjoint = Curr_El->get_adjoint();
					double *correction = Curr_El->get_correction();
					double dt = timeprops_ptr->dt.at(iter - 1); // if we have n iter size of dt vector is n-1
					double dtdx = dt / dx[0];
					double dtdy = dt / dx[1];
					double tiny = GEOFLOW_TINY;
					double *d_state_vars = Curr_El->get_d_state_vars();

					if (timeprops_ptr->iter < 52)
						matprops_ptr->frict_tiny = 0.1;
					else
						matprops_ptr->frict_tiny = 0.000000001;

					double flux[4][NUM_STATE_VARS];
					get_flux(El_Table, NodeTable, Curr_El->pass_key(), matprops_ptr, myid, flux);

					int stop[2];
					double tmp[] = { 0., 0., 0. };

					update_states(tmp, prev_state_vars, flux[0], flux[1], flux[2], flux[3], dtdx, dtdy, dt,
					    d_state_vars, (d_state_vars + NUM_STATE_VARS), curvature, matprops_ptr->intfrict,
					    bedfrict, gravity, d_gravity, *(Curr_El->get_kactxy()), matprops_ptr->frict_tiny,
					    stop, orgSrcSgn);

					el_error[0] = el_error[1] = 0.;
					*correction = 0.;

					for (int j = 0; j < NUM_STATE_VARS; j++) {

						vec_res[j] = tmp[j] - state_vars[j];

						el_error[1] += fabs(vec_res[j] * (adjoint[j] - bilin_adj[j]));

						*correction += vec_res[j] * adjoint[j];

//						double r_adj = 0., r_res = 0., r_prev = 0., r_state = 0.;
//
//						if (adjoint[j] != 0)
//							r_adj = fabs((adjoint[j] - bilin_adj[j]) / adjoint[j]);
//						if (state_vars[j] != 0)
//							r_res = fabs(vec_res[j] / state_vars[j]);
//						if (*(Curr_El->get_state_vars() + j) != 0.)
//							r_state = fabs(
//							    (state_vars[j] - *(Curr_El->get_state_vars() + j))
//							        / *(Curr_El->get_state_vars() + j));
//
//						if (*(Curr_El->get_prev_state_vars() + j) != 0.)
//							r_prev = fabs(
//							    (prev_state_vars[j] - *(Curr_El->get_prev_state_vars() + j))
//							        / *(Curr_El->get_prev_state_vars() + j));
//
//						if ((r_adj > .1 || r_state > .1 || r_res > .1 || r_prev > .1) && el_error[1] > 1e-7)
//							printf(
//							    "ratio adj %-8f ratio res %-8f ratio prev  %-8f  ratio bilin %-8f error %-8e  state %-1d key: %-u, %-u \n",
//							    r_adj, r_res, r_prev, r_state, el_error[1], j, *(Curr_El->pass_key()),
//							    *(Curr_El->pass_key() + 1));

					}

				}
				currentPtr = currentPtr->next;
			}
		}

//	meshplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., 2);
#ifdef DEBUG
	if (checkElement(El_Table))
	exit(22);
//		(timeprops_ptr->iter)++;
//		tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr,
//				dummyv_star, adjflag);
	cout << "number of elements -7   " << num_nonzero_elem(El_Table, -7) << endl
	<< "number of elements -6   " << num_nonzero_elem(El_Table, -6) << endl
	<< "number of elements  0   " << num_nonzero_elem(El_Table, 0) << endl
	<< "number of elements  1   " << num_nonzero_elem(El_Table, 1) << endl
	<< "number of elements  2   " << num_nonzero_elem(El_Table, 2) << endl
	<< "number of elements  3   " << num_nonzero_elem(El_Table, 3) << endl
	<< "number of elements  4   " << num_nonzero_elem(El_Table, 4) << endl
	<< "number of elements  5   " << num_nonzero_elem(El_Table, 5) << endl;
//		getchar();
	int nonz = num_nonzero_elem(El_Table);

	cout << "number of elements after refinement  " << nonz << endl;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {

					int index = Curr_El->get_sol_rec_ind();
					dbgvec[index] += 1;
					pass[index] = 1;

				}
			}
		}
	}

	for (int i = 0; i < nonz1; i++) {
		if (dbgvec[i] != 4) {
			cout << "these elements have problem:   " << i << endl
			<< "the value is  " << dbgvec[i] << endl;
			exit(EXIT_FAILURE);

		} else if (pass[i] != 1) {
			cout << "this index has not been passed  " << i << endl;
			exit(EXIT_FAILURE);

		}
	}
#endif

}

