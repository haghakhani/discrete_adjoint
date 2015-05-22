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

void calc_adjoint_elem(HashTable* El_Table, vector<Jacobian*>* solHyst,
		Element *Curr_El, int iter, int adjiter, int myid);

void calc_adjoint(HashTable* El_Table, vector<Jacobian*>* solHyst, int iter,
		int adjiter, int myid) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;

	double aa, bb = .1;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->pass_key()) == KEY0
							&& *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER)
						aa = bb;

					calc_adjoint_elem(El_Table, solHyst, Curr_El, iter, adjiter, myid);
				}
				currentPtr = currentPtr->next;
			}
		}
	}
	return;
}

void calc_adjoint_elem(HashTable* El_Table, vector<Jacobian*>* solHyst,
		Element *Curr_El, int iter, int adjiter, int myid) {

	double* adjoint = (Curr_El->get_state_vars() + 6);
	Jacobian *jacobian, *neighjac;
	jacobian = solHyst->at(Curr_El->get_sol_rec_ind());

	if (adjiter == 0) {

		for (int i = 0; i < 3; ++i)
			adjoint[i] = *(jacobian->get_funcsens(iter) + i);

	} else {

		Element *neigh_elem;
		double* adjoint_pointer;
		double adjcontr[3] = { 0.0, 0.0, 0.0 };
		double*** jacobianmat;

		for (int effelement = 0; effelement < 5; effelement++) { //0 for the element itself, and the rest id for neighbour elements

			if (effelement == 0) {		    //this part of code

				adjoint_pointer = (Curr_El->get_prev_state_vars() + 6);

				jacobianmat = jacobian->get_jacobian();

				for (int k = 0; k < 3; ++k)
					for (int l = 0; l < 3; ++l)
						adjcontr[k] += adjoint_pointer[l] * jacobianmat[0][k][l];

			} else {

				neigh_elem = Curr_El->get_side_neighbor(El_Table, effelement - 1);//basically we are checking all neighbor elements, and start from xp neighbor
				if (neigh_elem) {

					adjoint_pointer = (neigh_elem->get_prev_state_vars() + 6);
					neighjac = solHyst->at(neigh_elem->get_sol_rec_ind());
					jacobianmat = neighjac->get_jacobian();

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

					for (int k = 0; k < 3; ++k)
						for (int l = 0; l < 3; ++l)
							adjcontr[k] += adjoint_pointer[l] * jacobianmat[jacind][k][l];

				}
			}
		}

		for (int j = 0; j < 3; j++)
			adjoint[j] = *(jacobian->get_funcsens(iter-1) + j) - adjcontr[j];
	}

	for (int i = 0; i < 3; i++)
		if (isnan(adjoint[i]) || isinf(adjoint[i]))
			cout << "it is incorrect  " << endl;

	return;
}
