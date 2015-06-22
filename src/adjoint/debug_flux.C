/*
 * debug_flux.C
 *
 *  Created on: May 15, 2015
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

void record_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key,
		MatProps* matprops_ptr, int effelement, int myid, double* fluxxpold,
		double* fluxypold, double* fluxxmold, double* fluxymold) {

	double dummydt = 0., outflow = 0.;
	int order_flag = 1;
	ResFlag resflag;
	resflag.callflag = 1;
	resflag.lgft = 0;

	Element* Curr_El = (Element*) (El_Table->lookup(key));

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	Element* neigh_elem;

	if (effelement == 0) {

		Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt,
				&order_flag, &outflow, resflag, resflag); //chenge of the state_vars causes the all around fluxes change, this update xp ,yp

		if ((*(Curr_El->get_neigh_proc() + xm)) != INIT) { //we have to make sure that there exit an element in xm side

			Element* elem_xm = (Element*) (El_Table->lookup(
					Curr_El->get_neighbors() + xm * KEYLENGTH));
			assert(elem_xm);
			elem_xm->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid,
					dummydt, &order_flag, &outflow, resflag, resflag); //this update the flux on share edge with xm

		}

		if ((*(Curr_El->get_neigh_proc() + ym)) != INIT) { //we have to make sure that there exit an element in ym side

			Element* elem_ym = (Element*) (El_Table->lookup(
					Curr_El->get_neighbors() + ym * KEYLENGTH));
			assert(elem_ym);
			elem_ym->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid,
					dummydt, &order_flag, &outflow, resflag, resflag); //this update the flux on share edge with ym
		}

	} else {
		if ((*(Curr_El->get_neigh_proc() + (effelement - 1))) != INIT) {

			neigh_elem = (Element*) (El_Table->lookup(
					Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
			assert(neigh_elem);

			if ((effelement - 1) == xp || (effelement - 1) == yp)
				Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid,
						dummydt, &order_flag, &outflow, resflag,resflag); //if we change the state variables in xp or yp, just the flux at this element has to be updated
			else
				neigh_elem->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid,
						dummydt, &order_flag, &outflow, resflag, resflag); //otherwise the flux at the corresponding element has to be updated

		}
	}

	Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

	Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

	Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

	Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
		fluxxpold[ivar] = nxp->flux[ivar];
		fluxypold[ivar] = nyp->flux[ivar];
		fluxxmold[ivar] = nxm->flux[ivar];
		fluxymold[ivar] = nym->flux[ivar];
	}
}

void flux_debug(Element* Curr_El, double* fluxxpold, double* fluxxmold, double* fluxypold,
    double* fluxymold, double* fluxxp, double* fluxxm, double* fluxyp, double* fluxym,
    int effelement, int j, int iter, double dt) {

	double diff1[NUM_STATE_VARS], diff2[NUM_STATE_VARS], diff3[NUM_STATE_VARS], diff4[NUM_STATE_VARS];
	double fluxxold[NUM_STATE_VARS], fluxxnew[NUM_STATE_VARS], fluxyold[NUM_STATE_VARS], fluxynew[NUM_STATE_VARS];
	double abs_fluxx_diff[NUM_STATE_VARS], abs_fluxy_diff[NUM_STATE_VARS];

	double *state_vars = Curr_El->get_state_vars();
	double *prev_state_vars = Curr_El->get_prev_state_vars();
	double *d_state_vars = Curr_El->get_d_state_vars();
	double *dx = Curr_El->get_dx();
	double dtdx = dt / dx[0];
	double dtdy = dt / dx[1];

	FILE *fp;
	fp = fopen("debugfile", "a");
	fprintf(fp, "this is for jacobian corresponding to state vars %d \n", j);
	fprintf(fp,
	    "In corrector time step %d with dt=%f dtdx=%f dtdy=%f kactx=%f , kacty=%f , x=%6f, y=%6f \n state vars are: \n",
	    iter, dt, dtdx, dtdy, *(Curr_El->get_kactxy()), *(Curr_El->get_kactxy() + 1),
	    *(Curr_El->get_coord()), *(Curr_El->get_coord() + 1));

	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", state_vars[state]);
	fprintf(fp, "\n prev state vars are: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", prev_state_vars[state]);
	fprintf(fp, "\n xp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxxp[state]);
	fprintf(fp, "\n xm: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxxm[state]);
	fprintf(fp, "\n yp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxyp[state]);
	fprintf(fp, "\n ym: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxym[state]);
	fprintf(fp, "\n dUdx: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state]);
	fprintf(fp, "\n dUdy: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state + NUM_STATE_VARS]);
	fprintf(fp, "\n");

	fclose(fp);

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		diff1[ivar] = dabs(fluxxpold[ivar] - fluxxp[ivar]);
		diff2[ivar] = dabs(fluxxmold[ivar] - fluxxm[ivar]);
		diff3[ivar] = dabs(fluxypold[ivar] - fluxyp[ivar]);
		diff4[ivar] = dabs(fluxymold[ivar] - fluxym[ivar]);
	}

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		fluxxold[ivar] = fluxxpold[ivar] - fluxxmold[ivar];
		fluxxnew[ivar] = fluxxp[ivar] - fluxxm[ivar];
		fluxyold[ivar] = fluxypold[ivar] - fluxymold[ivar];
		fluxynew[ivar] = fluxyp[ivar] - fluxym[ivar];
	}

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		abs_fluxx_diff[ivar] = dabs(fluxxold[ivar] - fluxxnew[ivar]);
		if (abs_fluxx_diff[ivar] > 0. && ivar != 1)
			cout << setw(10) << setprecision(8) << "change effects the flux_x for var  " << ivar
			    << " eff_elem   " << effelement << "  j  " << j << "  value  " << abs_fluxx_diff[ivar]
			    << endl;

		abs_fluxy_diff[ivar] = dabs(fluxxold[ivar] - fluxxnew[ivar]);
		if (abs_fluxy_diff[ivar] > 0. && ivar != 1)
			cout << "change effects the flux_y for var  " << ivar << " eff_elem   " << effelement
			    << "  j  " << j << "  value  " << abs_fluxy_diff[ivar] << endl;
	}

	return;
}

