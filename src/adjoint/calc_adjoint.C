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
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define KEY0 3938123776
#define KEY1 0
#define ITER 141

struct Func_CTX {

	double dx;
	double dy;
	TimeProps* timeprops;
	int iter;
	int adjiter;
};

void calc_adjoint(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;
	int iter = propctx->timeprops->iter;

	double aa, bb = .1;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER)
						aa = bb;

					calc_adjoint_elem(meshctx, propctx, Curr_El);
				}
				currentPtr = currentPtr->next;
			}
		}
	}
	return;
}

void calc_adjoint_elem(MeshCTX* meshctx, PropCTX* propctx, Element *Curr_El) {

	Func_CTX ctx;
	ctx.dx = *(Curr_El->get_dx());
	ctx.dy = *(Curr_El->get_dx() + 1);
	ctx.timeprops = propctx->timeprops;
	ctx.iter = propctx->timeprops->iter;
	ctx.adjiter = propctx->timeprops->adjiter;

	double* adjoint = Curr_El->get_adjoint();

	Curr_El->calc_func_sens((const void *) &ctx);

	HashTable* El_Table = meshctx->el_table;

	if (propctx->timeprops->adjiter == 0) {

		for (int i = 0; i < NUM_STATE_VARS; ++i)
			adjoint[i] = *(Curr_El->get_func_sens() + i);

	} else {

		Element *neigh_elem;
		double* adjoint_pointer;
		double adjcontr[NUM_STATE_VARS] = { 0.0, 0.0, 0.0 };
		double*** jacobianmat;

		for (int effelement = 0; effelement < EFF_ELL; effelement++) { //0 for the element itself, and the rest id for neighbour elements

			if (effelement == 0) {		    //this part of code

				adjoint_pointer = (Curr_El->get_prev_adjoint());

				jacobianmat = Curr_El->get_jacobian();

				for (int k = 0; k < NUM_STATE_VARS; ++k)
					for (int l = 0; l < NUM_STATE_VARS; ++l)
						adjcontr[k] += adjoint_pointer[l] * jacobianmat[0][k][l];

			} else {

				neigh_elem = Curr_El->get_side_neighbor(El_Table, effelement - 1);//basically we are checking all neighbor elements, and start from xp neighbor
				if (neigh_elem) {

					adjoint_pointer = neigh_elem->get_prev_adjoint();
					jacobianmat = neigh_elem->get_jacobian();

					int jacind;

					switch (effelement) {
						case 1:	//in xp neighbor I have to read jacobian of xm, because position of curr_el is in xm side of that neighbor
							jacind = 3;
							break;
						case 2:		    //for yp return ym
							jacind = 4;
							break;
						case 3:		    //for xm return xp
							jacind = 1;
							break;
						case 4:		    //for ym return yp
							jacind = 2;
							break;
						default:
							cout << "invalid neighbor position" << endl;
					}

					for (int k = 0; k < NUM_STATE_VARS; ++k)
						for (int l = 0; l < NUM_STATE_VARS; ++l)
							adjcontr[k] += adjoint_pointer[l] * jacobianmat[jacind][k][l];

				}
			}
		}

		for (int j = 0; j < NUM_STATE_VARS; j++)
			adjoint[j] = *(Curr_El->get_func_sens() + j) - adjcontr[j];
	}

	for (int i = 0; i < NUM_STATE_VARS; i++)
		if (isnan(adjoint[i]) || isinf(adjoint[i]))
			cout << "it is incorrect  " << endl;

	return;
}
void Element::calc_func_sens(const void * ctx) {

	Func_CTX* contx = (Func_CTX *) ctx;

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		func_sens[i] = 0.;

	if (contx->adjiter == 0) {

		double dt_n = contx->timeprops->dt.at(contx->iter - 1);
		func_sens[0] = .5 * state_vars[0] * dt_n * contx->dx * contx->dy;

	} else if (contx->iter == 1) {

		double dt_n = contx->timeprops->dt.at(contx->iter - 1);
		func_sens[0] = .5 * prev_state_vars[0] * dt_n * contx->dx * contx->dy;

	} else {

		double dt_n = contx->timeprops->dt.at(contx->iter - 1);
		double dt_p = contx->timeprops->dt.at(contx->iter - 2);
		func_sens[0] = .5 * prev_state_vars[0] * (dt_n + dt_p) * contx->dx * contx->dy;

	}
}

