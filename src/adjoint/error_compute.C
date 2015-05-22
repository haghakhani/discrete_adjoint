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

void error_compute(HashTable* El_Table, HashTable* NodeTable,
		TimeProps* timeprops_ptr, MatProps* matprops_ptr, int iter, int myid,
		int numprocs) {

	setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr,
			timeprops_ptr);

	int order_flag = 1;	//this is dummy here
	double outflow[1];	//this is dummy here
	ResFlag resflag;
	resflag.callflag = 1;
	resflag.lgft = 0;
	calc_edge_states(El_Table, NodeTable, matprops_ptr, timeprops_ptr, myid,
			&order_flag, outflow, resflag);

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;

	if (iter != 0) {
		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() == NEWSON) {

						int dbgflag;
						if (*(Curr_El->pass_key()) == KEY0
								&& *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER)
							dbgflag = 1;

						double *state_vars = Curr_El->get_state_vars();
						double *prev_state_vars = Curr_El->get_prev_state_vars();
						double *gravity = Curr_El->get_gravity();
						double *d_gravity = Curr_El->get_d_gravity();
						double *curvature = Curr_El->get_curvature();
						Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
						double bedfrict = Curr_El->get_effect_bedfrict();
						double velocity[DIMENSION];
						double *dx = Curr_El->get_dx();
						double kactxy[DIMENSION];
						double orgSrcSgn[2], vec_res[3];

						Curr_El->get_slopes_prev(El_Table, NodeTable, matprops_ptr->gamma);
						double *d_state_vars = Curr_El->get_d_state_vars();

						if (timeprops_ptr->iter < 50)
							matprops_ptr->frict_tiny = 0.1;
						else
							matprops_ptr->frict_tiny = 0.000000001;

						orgSourceSgn(Curr_El, matprops_ptr->frict_tiny, orgSrcSgn);

						double dt = timeprops_ptr->dt.at(iter - 1); // if we have n iter size of dt vector is n-1
						double dtdx = dt / dx[0];
						double dtdy = dt / dx[1];

						double *el_error = Curr_El->get_el_error();

						double *constAdj = Curr_El->get_const_adj();
						double *correction = Curr_El->get_correction();

						int xp = Curr_El->get_positive_x_side(); //finding the direction of element
						int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

						const int state_num = NUM_STATE_VARS - 2;

						Node* nxp = (Node*) NodeTable->lookup(
								Curr_El->getNode() + (xp + 4) * 2);

						Node* nyp = (Node*) NodeTable->lookup(
								Curr_El->getNode() + (yp + 4) * 2);

						Node* nxm = (Node*) NodeTable->lookup(
								Curr_El->getNode() + (xm + 4) * 2);

						Node* nym = (Node*) NodeTable->lookup(
								Curr_El->getNode() + (ym + 4) * 2);

						double fluxxp[state_num], fluxyp[state_num]; //we just need to compute jacobian for h,u,v not for cont. adjoint so we don't store the fluxes for the adjoint
						double fluxxm[state_num], fluxym[state_num];

						for (int ivar = 0; ivar < state_num; ivar++)
							fluxxp[ivar] = nxp->flux[ivar];

						for (int ivar = 0; ivar < state_num; ivar++)
							fluxyp[ivar] = nyp->flux[ivar];

						for (int ivar = 0; ivar < state_num; ivar++)
							fluxxm[ivar] = nxm->flux[ivar];

						for (int ivar = 0; ivar < state_num; ivar++)
							fluxym[ivar] = nym->flux[ivar];

						if (*(Curr_El->pass_key()) == KEY0
								&& *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER)
							dbgflag = 1;

						residual(vec_res, state_vars, prev_state_vars, fluxxp, fluxyp,
								fluxxm, fluxym, dtdx, dtdy, dt, d_state_vars,
								(d_state_vars + NUM_STATE_VARS), curvature,
								matprops_ptr->intfrict, bedfrict, gravity, d_gravity,
								Curr_El->get_kactxy(), matprops_ptr->frict_tiny, orgSrcSgn,
								0./*here increment is zero*/, matprops_ptr->epsilon);

						state_vars[1] = vec_res[0];
						state_vars[4] = vec_res[1];
						state_vars[5] = vec_res[2];

						el_error[1] = 0.0;
						*correction = 0.0;
						for (int j = 0; j < 3; j++) {
							el_error[1] += vec_res[j]
									* (state_vars[NUM_STATE_VARS + j] - constAdj[j]);

							*correction += vec_res[j] * state_vars[NUM_STATE_VARS + j];
						}
						//if (el_error[1]!=0)
						//	cout<<"it should print blue"<<endl ;

					}
					currentPtr = currentPtr->next;
				}
			}
	}

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
	return;
}
