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
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: get_coef_and_eigen.C,v 1.4 2004/08/11 15:58:46 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define DO_EROSION
//#define DEBUG

#define KEY0   2198474436
#define KEY1   3964585199
#define ITER   11
#define EFFELL 1
#define J      0
#define JACIND 0
#define X      1.97853625e-01
#define Y      2.30663500e-01

#include "../header/hpfem.h"

void correct(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
    FluxProps *fluxprops, TimeProps *timeprops, Element *EmTemp) {
//	Element *EmTemp = (Element *) EmTemp_in;
	double *dx = EmTemp->get_dx();
	double dtdx = dt / dx[0];
	double dtdy = dt / dx[1];
	double kactxy[DIMENSION];

	int xp = EmTemp->get_positive_x_side();
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	int ivar, i;
	double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
	double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];

	Node* nxp = (Node*) NodeTable->lookup(EmTemp->getNode() + (xp + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxp[ivar] = nxp->flux[ivar];

	Node* nyp = (Node*) NodeTable->lookup(EmTemp->getNode() + (yp + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxyp[ivar] = nyp->flux[ivar];

	Node* nxm = (Node*) NodeTable->lookup(EmTemp->getNode() + (xm + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxm[ivar] = nxm->flux[ivar];

	Node* nym = (Node*) NodeTable->lookup(EmTemp->getNode() + (ym + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxym[ivar] = nym->flux[ivar];

	double *state_vars = EmTemp->get_state_vars();
	double *prev_state_vars = EmTemp->get_prev_state_vars();
	double *d_state_vars = EmTemp->get_d_state_vars();
	double *gravity = EmTemp->get_gravity();
	double *d_gravity = EmTemp->get_d_gravity();
	double *curvature = EmTemp->get_curvature();
	double *tan_bedfrict = EmTemp->get_tanbedfrict();

	double Vsolid[DIMENSION];

	if (state_vars[0] > GEOFLOW_TINY) {
		for (i = 0; i < DIMENSION; i++)
			kactxy[i] = *(EmTemp->get_kactxy() + i);

		double inv_h = 1. / state_vars[0];

		// fluid velocities
		Vsolid[0] = state_vars[1] * inv_h;
		Vsolid[1] = state_vars[2] * inv_h;

	} else {
		for (i = 0; i < DIMENSION; i++) {
			kactxy[i] = matprops_ptr->epsilon;
			Vsolid[i] = 0.;
		}
	}

	int debuging, ggg = 0;
	if (/*EmTemp->get_ithelem() == 8251*/ EmTemp->pass_key()[0] == KEY0 && EmTemp->pass_key()[1] == KEY1
	&& timeprops->iter ==51)
//	if (dabs(*(EmTemp->get_coord()) - X) < INCREMENT&& dabs(*(EmTemp->get_coord()+1) - Y)<INCREMENT
//	&& timeprops->iter == ITER)
		debuging = ggg = 1;

	int stop[2];
	double orgSrcSgn[4];

	update_states(state_vars, prev_state_vars, //2
	    fluxxp, fluxyp, fluxxm, fluxym, dtdx, //5
	    dtdy, dt, d_state_vars, (d_state_vars + NUM_STATE_VARS), //4
	    curvature, (matprops_ptr->intfrict), tan_bedfrict, gravity, //4
	    d_gravity, kactxy[0], matprops_ptr->frict_tiny, stop, orgSrcSgn);

#ifdef DEBUG
	if (*(EmTemp->pass_key()) == KEY0 && *(EmTemp->pass_key() + 1) == KEY1 && timeprops->iter == ITER) {
		EmTemp->write_elem_info(NodeTable, filename, timeprops->iter, dt);
	}
#endif

	bool print_vars = false;
	for (i = 0; i < NUM_STATE_VARS; i++)
		if (state_vars[i] > 1e3) // || dabs(state_vars[5])>5 )
			print_vars = false;

	if (print_vars) {
		double tempU[NUM_STATE_VARS];
		for (i = 0; i < NUM_STATE_VARS; i++)
			tempU[i] = prev_state_vars[i] - dtdx * (fluxxp[i] - fluxxm[i])
			    - dtdy * (fluxyp[i] - fluxym[i]);
		printf("ElemKey: %d   ElemKey2: %d\n ", *EmTemp->pass_key(), *(EmTemp->pass_key() + 1));
		printf("Kactxy = %10.5e, %10.5e\n", kactxy[0], kactxy[1]);
		printf("TAN BedFrict: %10.5e: IntFrict: %10.5e\n", *tan_bedfrict, matprops_ptr->intfrict);

		printf("state_vars: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", state_vars[i]);
		printf("\n");

		printf("prev_state_vars: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", prev_state_vars[i]);
		printf("\n");

		printf("Ustore: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", tempU[i]);
		printf("\n");

		printf("fluxes: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("fluxxp:%10.5e, fluxxm:%10.5e, fluxyp:%10.5e, fluxym:%10.5e \n ", fluxxp[i], fluxxm[i],
			    fluxyp[i], fluxym[i]);

		exit(1);
	}

	return;
}
