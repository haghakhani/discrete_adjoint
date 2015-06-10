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

#define KEY0 3916612844
#define KEY1 1321528399
#define ITER 1

struct Func_CTX {

	double dx;
	double dy;
	TimeProps* timeprops;
	int iter;
	int adjiter;
};

void calc_adjoint(DualMesh* dualmesh, TimeProps* timeprops_ptr, int iter, int adjiter, int myid) {

	int Ny = dualmesh->get_Ny();
	int Nx = dualmesh->get_Nx();

	Func_CTX ctx;
	ctx.dx = dualmesh->get_dx();
	ctx.dy = dualmesh->get_dx();
	ctx.timeprops = timeprops_ptr;
	ctx.iter = iter;
	ctx.adjiter = adjiter;

	for (int i = 0; i < Ny; ++i)
		for (int j = 0; j < Nx; ++j) {

			DualCell* dualcell = dualmesh->get_dualcell(i, j);
			dualcell->calc_func_sens((const void *) &ctx);

			if (adjiter == 0) {

				double* adjoint = dualcell->get_curr_adjoint();

				for (int k = 0; k < NUM_STATE_VARS; ++k)
					adjoint[k] = *(dualcell->get_funcsens() + k);

			} else {

				double* adjoint = dualcell->get_curr_adjoint();

				DualCell *neigh_cell;
				double* adjoint_pointer;
				double adjcontr[NUM_STATE_VARS] = { 0., 0., 0. };
				double*** jacobianmat;

				for (int effelement = 0; effelement < 5; effelement++) { //0 for the element itself, and the rest id for neighbour elements

					if (effelement == 0) {		    //this part of code

						adjoint_pointer = dualcell->get_prev_adjoint();

						jacobianmat = dualcell->get_jacobian();

						for (int k = 0; k < NUM_STATE_VARS; ++k)
							for (int l = 0; l < NUM_STATE_VARS; ++l)
								adjcontr[k] += adjoint_pointer[l] * jacobianmat[0][k][l];

					} else {

						int a, b;
						set_ab(&a, &b, effelement);

						neigh_cell = dualmesh->get_dualcell(i + a, j + b);//basically we are checking all neighbor elements, and start from xp neighbor

						if (neigh_cell) {

							adjoint_pointer = neigh_cell->get_prev_adjoint();
							jacobianmat = neigh_cell->get_jacobian();

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

				for (int k = 0; k < NUM_STATE_VARS; k++)
					adjoint[k] = *(dualcell->get_funcsens() + k) - adjcontr[k];

			}
		}

	return;
}

void set_ab(int* a, int* b, int effelement) {

	//b is for index in x, and a is index for y
	*a = *b = 0;
	switch (effelement) {
		case 1:
			*b = 1;
			break;
		case 2:
			*a = 1;
			break;
		case 3:
			*b = -1;
			break;
		case 4:
			*a = -1;
			break;
		default:
			cerr << "This effective element is invalis" << endl;
	}
}

void DualCell::calc_func_sens(const void * ctx) {

	Func_CTX* contx = (Func_CTX *) ctx;

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

