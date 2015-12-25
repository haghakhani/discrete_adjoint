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

#define DEBUGFILE1

#define KEY0   3781669179
#define KEY1   330382100
#define ITER   11

struct Func_CTX {

	double dx;
	double dy;
	TimeProps* timeprops;
	int iter;
	int adjiter;

};

void calc_func_sens(HashTable *El_Table);
int get_jacind(int effelement);

void calc_adjoint(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;
	int iter = propctx->timeprops->iter;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					int boundary = 0;
					//this part handles if the Curr_El is a boundary element
					for (int neighnum = 0; neighnum < 4; neighnum++)
						if (*(Curr_El->get_neigh_proc() + neighnum) == INIT) {
							boundary = 1;
							for (int j = 0; j < NUM_STATE_VARS; ++j)
								*(Curr_El->get_adjoint() + j) = 0.;

							break;
						}
					if (!boundary)
						calc_adjoint_elem(meshctx, propctx, Curr_El);

				}
				currentPtr = currentPtr->next;
			}
		}
	}

}

void calc_adjoint_elem(MeshCTX* meshctx, PropCTX* propctx, Element *Curr_El) {

	Func_CTX ctx;
	ctx.dx = *(Curr_El->get_dx());
	ctx.dy = *(Curr_El->get_dx() + 1);
	ctx.timeprops = propctx->timeprops;
	ctx.iter = propctx->timeprops->iter;
	ctx.adjiter = propctx->timeprops->adjiter;

	double* adjoint = Curr_El->get_adjoint();
	double adjcontr[NUM_STATE_VARS] = { 0., 0., 0. };

//	Curr_El->calc_func_sens((const void *) &ctx);

	HashTable* El_Table = meshctx->el_table;

#ifdef DEBUGFILE
	ofstream myfile;
	myfile.open("adjelem.txt", ios::app);
	ofstream myfile;
	myfile.open("adjdebug.txt", ios::app);

	myfile << "Elem Key[0]= " << *(Curr_El->pass_key()) << "  Key[1]= " << *(Curr_El->pass_key())
	<< " iter= " << propctx->timeprops->iter << endl;

#endif

	if (propctx->timeprops->adjiter == 0) {

		Curr_El->calc_func_sens((void*) propctx);

		for (int i = 0; i < NUM_STATE_VARS; ++i)
			adjoint[i] = *(Curr_El->get_func_sens() + i);

	} else {

		for (int effelement = 0; effelement < EFF_ELL; effelement++) { //0 for the element itself, and the rest id for neighbour elements

			if (effelement == 0) {		    //this part of code

				double* adjoint_prev = (Curr_El->get_prev_adjoint());

				Vec_Mat<9>& jacobianmat = Curr_El->get_jacobian();

				for (int k = 0; k < NUM_STATE_VARS; ++k)
					for (int l = 0; l < NUM_STATE_VARS; ++l)
						adjcontr[k] += adjoint_prev[l] * jacobianmat(effelement, l, k);

			} else if (effelement <= 4
			    || (effelement > 4 && *(Curr_El->get_neigh_proc() + (effelement - 1)) >= 0)) {

				//basically we are checking all neighbor elements, and start from xp neighbor
				Element * neigh_elem = (Element*) (El_Table->lookup(
				    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));

				if (neigh_elem) {

					double* adjoint_prev = neigh_elem->get_prev_adjoint();
					Vec_Mat<9>& jacobianmat = neigh_elem->get_jacobian();

					int jacind = neigh_elem->which_neighbor(Curr_El->pass_key());
					// because we have to consider the element itself which is in jacind=0
					jacind++;

					for (int k = 0; k < NUM_STATE_VARS; ++k)
						for (int l = 0; l < NUM_STATE_VARS; ++l)
							adjcontr[k] += adjoint_prev[l] * jacobianmat(jacind, l, k);

				}
			}
		}

		for (int j = 0; j < NUM_STATE_VARS; j++)
			adjoint[j] = /**(Curr_El->get_func_sens() + j)*/-adjcontr[j];
	}

	for (int i = 0; i < NUM_STATE_VARS; i++)
		if (isnan(adjoint[i]) || isinf(adjoint[i]))
			cout << "it is incorrect  " << endl;

#ifdef DEBUGFILE
	ofstream adjdebug;
	adjdebug.open("adjdebug.txt", ios::app);

	adjdebug << "Elem Key[0]= " << *(Curr_El->pass_key()) << "  Key[1]= "
	<< *(Curr_El->pass_key() + 1) << " iter= " << propctx->timeprops->iter << " elem pos x= "
	<< *(Curr_El->get_coord()) << " y=" << *(Curr_El->get_coord() + 1) << endl;
	for (int k = 0; k < NUM_STATE_VARS; ++k)
	adjdebug << " adjoint[" << k << "]= " << adjoint[k];
	adjdebug << "\n";

	myfile.close();
#endif
}
void Element::calc_func_sens(const void * ctx) {

	Func_CTX* contx = (Func_CTX *) ctx;

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		func_sens[i] = 0.;

	if (prev_state_vars[0] > 0.)
		func_sens[0] = dx[0] * dx[1];

//	if (contx->adjiter == 0) {
//
//		double dt_n = contx->timeprops->dt.at(contx->iter - 1);
//		func_sens[0] = .5 * state_vars[0] * dt_n * contx->dx * contx->dy;
//
//	} else if (contx->iter == 1) {
//
//		double dt_n = contx->timeprops->dt.at(contx->iter - 1);
//		func_sens[0] = .5 * prev_state_vars[0] * dt_n * contx->dx * contx->dy;
//
//	} else {
//
//		double dt_n = contx->timeprops->dt.at(contx->iter - 1);
//		double dt_p = contx->timeprops->dt.at(contx->iter - 2);
//		func_sens[0] = .5 * prev_state_vars[0] * (dt_n + dt_p) * contx->dx * contx->dy;
//
//	}
}

void calc_func_sens(HashTable *El_Table) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;
	double max = 0.;

//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//
//				if (Curr_El->get_adapted_flag() > 0 )
//					Curr_El->calc_func_sens((void *) &max);
//
//				currentPtr = currentPtr->next;
//			}
//		}

//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//
//				if (Curr_El->get_adapted_flag() > 0 && *(Curr_El->get_state_vars()) > max)
//					max = *(Curr_El->get_state_vars());
//
//				currentPtr = currentPtr->next;
//			}
//		}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					if (*(Curr_El->get_state_vars()) > 0.) {
						*(Curr_El->get_func_sens()) = -*(Curr_El->get_dx()) + *(Curr_El->get_dx() + 1);
						*(Curr_El->get_func_sens() + 1) = *(Curr_El->get_func_sens() + 2) = 0.;
					} else {
						*(Curr_El->get_func_sens()) = *(Curr_El->get_func_sens() + 1) =
						    *(Curr_El->get_func_sens() + 2) = 0.;
					}

				}
				currentPtr = currentPtr->next;
			}
		}

}

int get_jacind(int effelement) {
	int jacind;

	switch (effelement) {
		case 1:	//in xp neighbor I have to read jacobian of xm, because position of curr_el is in xm side of that neighbor
			jacind = 3;
			break;
		case 2:	//for yp return ym
			jacind = 4;
			break;
		case 3:	//for xm return xp
			jacind = 1;
			break;
		case 4:	//for ym return yp
			jacind = 2;
			break;
		default:
			cout << "invalid neighbor position" << endl;
	}
	return jacind;
}

