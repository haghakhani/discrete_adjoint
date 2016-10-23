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

#define KEY0   3937086542
#define KEY1   3303820997
#define ITER   31

struct Func_CTX {
	MeshCTX* meshctx;
	PropCTX* propctx;

};

void calc_func_sens(MeshCTX* meshctx, PropCTX* propctx);
int get_jacind(int effelement);
void sens_on_boundary(MeshCTX* meshctx, PropCTX* propctx, DualElem* eff_el, int side);
void calc_func_sens_hmax(MeshCTX* meshctx, PropCTX* propctx);

void calc_adjoint(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	DualElem* Curr_El = NULL;
	int iter = propctx->timeprops->iter;
	int aa = 0, bb = 1;

//	if (propctx->timeprops->adjiter == 0)
//		calc_func_sens_hmax(meshctx, propctx);

	calc_func_sens(meshctx, propctx);

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);

//				if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER)
				if (Curr_El->get_ithelem() == 338)
					bb = aa;
//				if (fabs(*(Curr_El->get_coord()) - 161.3) < .04
//				    && fabs(*(Curr_El->get_coord() + 1) - 539.58) < .04 && (iter == 545 || iter == 546))
//					bb = aa;
				if (Curr_El->get_adapted_flag() > 0) {
					calc_adjoint_elem(meshctx, propctx, Curr_El);

				}
				currentPtr = currentPtr->next;
			}
		}
	}

	move_dual_data(meshctx, propctx);

	apply_adjoint_bc(meshctx, propctx);

}

void calc_adjoint_elem(MeshCTX* meshctx, PropCTX* propctx, DualElem *Curr_El) {

	double* adjoint = Curr_El->get_adjoint();
	double adjcontr[NUM_STATE_VARS] = { 0., 0., 0. };

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

//		Curr_El->calc_func_sens((void*) propctx);

		for (int i = 0; i < NUM_STATE_VARS; ++i)
			adjoint[i] = -*(Curr_El->get_func_sens() + i);

	} else {

		for (int effelement = 0; effelement < EFF_ELL; effelement++) { //0 for the element itself, and the rest id for neighbour elements

			if (effelement == 0) {
				bool boundary = false;
				for (int bound = 0; bound < 4; ++bound)
					if (Curr_El->get_neigh_proc()[bound] == INIT)
						boundary = true;

				if (!boundary) {

					double* adjoint_prev = Curr_El->get_prev_adjoint();
					double* adjoint_pre3 = Curr_El->get_pre3_adjoint();

					Vec_Mat<9>& jacobianmat = Curr_El->get_jacobian();

					double coef = 0.25;
					if (propctx->timeprops->adjiter < 3)
						for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
							adjoint_pre3[ind] = adjoint_prev[ind];

					double dry = 1.;
					if (Curr_El->get_prev_state_vars()[0] == 0.)
						dry = 0.;

					for (int k = 0; k < NUM_STATE_VARS; ++k)
						for (int l = 0; l < NUM_STATE_VARS; ++l)
							if (k == l)
								adjcontr[k] += -0.75 * adjoint_prev[l]
								    * (1. * dry - 2. * jacobianmat(effelement, l, k))
								    - coef * adjoint_pre3[l] * dry;
							else
								adjcontr[k] += 1.5 * adjoint_prev[l] * jacobianmat(effelement, l, k);
				}

			} else if (effelement <= 4
			    || (effelement > 4 && *(Curr_El->get_neigh_proc() + (effelement - 1)) > -2)) {

				//basically we are checking all neighbor elements, and start from xp neighbor
				DualElem * neigh_elem = (DualElem*) (El_Table->lookup(
				    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));

				if (neigh_elem) {

					double* adjoint_prev = neigh_elem->get_prev_adjoint();
					Vec_Mat<9>& jacobianmat = neigh_elem->get_jacobian();

					int jacind = neigh_elem->which_neighbor(Curr_El->pass_key());
					// because we have to consider the element itself which is in jacind=0
					jacind++;

					for (int k = 0; k < NUM_STATE_VARS; ++k)
						for (int l = 0; l < NUM_STATE_VARS; ++l)
							adjcontr[k] += 1.5 * adjoint_prev[l] * jacobianmat(jacind, l, k);

				}
			}
		}

		for (int j = 0; j < NUM_STATE_VARS; j++)
			adjoint[j] = -*(Curr_El->get_func_sens() + j) - adjcontr[j];
	}

//	int cc = 1, bb = 0;
//	for (int i = 0; i < NUM_STATE_VARS; i++)
//		if (fabs(adjoint[i]) < 1.e-16)
//			adjoint[i]=0.;

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

void apply_adjoint_bc(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	DualElem* Curr_El = NULL;
	int iter = propctx->timeprops->iter;
	int myid = propctx->myid;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					for (int neighbor = 0; neighbor < 8; ++neighbor)
						if (Curr_El->get_neigh_proc()[neighbor] >= 0) {

							DualElem* neig_elem = (DualElem*) El_Table->lookup(
							    Curr_El->get_neighbors() + KEYLENGTH * neighbor);
							for (int side = 0; side < 4; ++side)
								if (neig_elem->get_neigh_proc()[side] == INIT && (neighbor % 4) == side) {

									int diff_gen = neig_elem->get_gen() - Curr_El->get_gen();

									switch (diff_gen) {
										case -1:
											// Curr_El with another cell is neighbor of a cell adjacent to interface
											//so just half of that cell will effect its adjoint
											for (int comp = 0; comp < NUM_STATE_VARS; ++comp)
												Curr_El->get_adjoint()[comp] += .5 * neig_elem->get_adjoint()[comp];
											break;
										case 0:
											for (int comp = 0; comp < NUM_STATE_VARS; ++comp)
												Curr_El->get_adjoint()[comp] += neig_elem->get_adjoint()[comp];
											break;
										case 1:
											//there has to be another neighbor cell of Curr_El which
											//its effect will be added later or has been added earlier
											for (int comp = 0; comp < NUM_STATE_VARS; ++comp)
												Curr_El->get_adjoint()[comp] += .5 * neig_elem->get_adjoint()[comp];
											break;
										default:
											exit(1);
											break;
									}
								}
						}
				currentPtr = currentPtr->next;
			}
		}
}

void DualElem::calc_func_sens(const void * ctx) {

	Func_CTX* contx = (Func_CTX *) ctx;
	MeshCTX* meshctx = contx->meshctx;
	PropCTX* propctx = contx->propctx;

	if (propctx->timeprops->adjiter == 0)
		if (state_vars[0] > 0.)
			func_sens[0] += dx[0] * dx[1];

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

void calc_func_sens_hmax(MeshCTX* meshctx, PropCTX* propctx) {
	HashTable* El_Table = meshctx->el_table;
	TimeProps* timeprops = propctx->timeprops;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;

	double hmax = 0.;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0 && hmax < *(Curr_El->get_state_vars())) {
					hmax = *(Curr_El->get_state_vars());
				}
				currentPtr = currentPtr->next;
			}
		}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0 && hmax == *(Curr_El->get_state_vars())) {
					*(Curr_El->get_func_sens()) = 1.;
				}

				for (int i = 0; i < NUM_STATE_VARS; i++)
					if (isnan(*(Curr_El->get_func_sens()+i)) || isinf(*(Curr_El->get_func_sens() + i)))
						cout << "it is incorrect  " << endl;

				currentPtr = currentPtr->next;
			}
		}

}

void calc_func_sens(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;
	TimeProps* timeprops = propctx->timeprops;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	int numproc = propctx->numproc;
	int myid = propctx->myid;
	int iter = timeprops->iter;
	double dt = timeprops->dt.at(iter - 1);

	// this for the elements that one of their sides is boundary
	// and the opposite side is in another processor
	// in following vector we pack the keys and side and
	// send it to other processors
	vector<unsigned> complicate_elements;
	vector<int> complicate_elements_side;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					for (int j = 0; j < NUM_STATE_VARS; ++j)
						*(Curr_El->get_func_sens() + j) = 0.;

				}
				currentPtr = currentPtr->next;
			}
		}

	//the following section is for taking care of the elements that their neighbors are adjacent to a boundary
	//and so the flux coming out from that boundary is depend on this cell
	//there is also a complicated case that the cell is a neighbor of a cell to another processor which that cell
	//is adjacent to the boundary
	if (propctx->timeprops->adjiter != 0)
		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					DualElem* Curr_El = (DualElem*) (currentPtr->value);

					if (Curr_El->get_adapted_flag() > 0) {

						int aa = 1, bb = 0;
						if (*(Curr_El->pass_key()) == 2868635045&& *(Curr_El->pass_key() + 1) == 2779096474
						&& propctx->timeprops->iter == ITER)
							bb = aa;

						int* neigh_proc = Curr_El->get_neigh_proc();
						int xp = Curr_El->get_positive_x_side(); //finding the direction of element
						int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

						for (int j = 0; j < 4; j++)
							if (neigh_proc[j] == INIT && j == ym) { // this is a boundary!

								double* dx = Curr_El->get_dx();
								double* func_sens = Curr_El->get_func_sens();

//								double flux[4][NUM_STATE_VARS];
//								get_flux(El_Table, NodeTable, Curr_El, flux);
//
//								if (j == xp && flux[0][0] != 0.) {
//									assert(flux[0][0] == flux[2][0]);
//									func_sens[1] += -dt * dx[1];
//								} else if (j == xm && flux[2][0] != 0.) {
//									assert(flux[0][0] == flux[2][0]);
//									func_sens[1] += dt * dx[1];
//								} else if (j == yp && flux[1][0] != 0.) {
//									assert(flux[1][0] == flux[3][0]);
//									func_sens[2] += -dt * dx[0];
//								} else if (j == ym && flux[3][0] != 0.) {
//									assert(flux[1][0] == flux[3][0]);
//									func_sens[2] += dt * dx[0];
//								}

								unsigned* neigh_key = Curr_El->get_neighbors();

								int opos_dirc = (j + 2) % 4;

								DualElem* eff_el = (DualElem*) El_Table->lookup(neigh_key + opos_dirc * KEYLENGTH);

								if (eff_el->get_myprocess() != myid) {
									complicate_elements.push_back(eff_el->pass_key()[0]);
									complicate_elements.push_back(eff_el->pass_key()[1]);
									complicate_elements_side.push_back(j);
								} else
									sens_on_boundary(meshctx, propctx, eff_el, j);

								// if the adjacent cell to the boundary has two neighbors on the opposite side
								if (neigh_proc[opos_dirc + 4] != -2) {
									DualElem* eff_el = (DualElem*) El_Table->lookup(
									    neigh_key + (opos_dirc + 4) * KEYLENGTH);

									if (eff_el->get_myprocess() != myid) {
										complicate_elements.push_back(eff_el->pass_key()[0]);
										complicate_elements.push_back(eff_el->pass_key()[1]);
										complicate_elements_side.push_back(j);
									} else
										sens_on_boundary(meshctx, propctx, eff_el, j);

								}

							}
					}
					currentPtr = currentPtr->next;
				}
			}

	int sizes[numproc];
	int size = complicate_elements.size();
	MPI_Allgather(&size, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

	int sum = 0;
	int mystart[numproc];
	for (int i = 0; i < numproc; ++i) {
		mystart[i] = sum;
		sum += sizes[i];
	}

	if (sum) {

		unsigned* check_elem = new unsigned[sum];

		MPI_Allgatherv(&complicate_elements[0], sizes[myid], MPI_UNSIGNED, check_elem, sizes, mystart,
		MPI_UNSIGNED, MPI_COMM_WORLD);

		for (int i = 0; i < numproc; ++i) {
			mystart[i] = mystart[i] / 2;
			sizes[i] = sizes[i] / 2;
		}

		int* check_side = new int[sum / 2];

		MPI_Allgatherv(&complicate_elements_side[0], sizes[myid], MPI_INT, check_side, sizes, mystart,
		MPI_INT, MPI_COMM_WORLD);

		for (int i = 0; i < sum / 2; i = i + 2) {
			unsigned key[] = { check_elem[i], check_elem[i + 1] };
			DualElem* eff_el = (DualElem*) El_Table->lookup(key);
			if (eff_el)
				sens_on_boundary(meshctx, propctx, eff_el, check_side[i]);
		}

		delete[] check_elem;
		delete[] check_side;

	}

//	Func_CTX contx;
//	contx.meshctx = meshctx;
//	contx.propctx = propctx;
//
//	if (timeprops->adjiter == 0)
//		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//			if (*(buck + i)) {
//				currentPtr = *(buck + i);
//				while (currentPtr) {
//					DualElem* Curr_El = (DualElem*) (currentPtr->value);
//
//					if (Curr_El->get_adapted_flag() > 0)
//						Curr_El->calc_func_sens(&contx);
//
//					currentPtr = currentPtr->next;
//				}
//			}

}

//void calc_func_sens(MeshCTX* meshctx, PropCTX* propctx) {
//
//	HashTable* El_Table = meshctx->el_table;
//	HashTable* NodeTable = meshctx->nd_table;
//	TimeProps* timeprops = propctx->timeprops;
//
//	HashEntryPtr currentPtr;
//	int numproc = propctx->numproc;
//	int myid = propctx->myid;
//	int iter = timeprops->iter;
//	double dt = timeprops->dt.at(iter - 1);
//
//	if (propctx->timeprops->adjiter == 0) {
//		HashEntryPtr* buck = El_Table->getbucketptr();
//		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//			if (*(buck + i)) {
//				currentPtr = *(buck + i);
//				while (currentPtr) {
//					DualElem* Curr_El = (DualElem*) (currentPtr->value);
//
//					if (Curr_El->get_adapted_flag() > 0)
//						for (int j = 0; j < NUM_STATE_VARS; ++j)
//							Curr_El->get_func_sens()[j] = 0.;
//
//					currentPtr = currentPtr->next;
//				}
//			}
//	} else {
//
//		DISCHARGE* discharge = propctx->discharge;
//		HashEntryPtr* buck = El_Table->getbucketptr();
//		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//			if (*(buck + i)) {
//				currentPtr = *(buck + i);
//				while (currentPtr) {
//					DualElem* Curr_El = (DualElem*) (currentPtr->value);
//
//					if (Curr_El->get_adapted_flag() > 0) {
//
//						unsigned *nodes = Curr_El->getNode();
//						double nodescoord[9][2], *coord;
//						Node* node;
//
//						for (int inode = 0; inode < 8; inode++) {
//							node = (Node*) NodeTable->lookup(nodes + 2 * inode);
//							coord = node->get_coord();
//
//							nodescoord[inode][0] = coord[0];
//							nodescoord[inode][1] = coord[1];
//
//						}
//						nodescoord[8][0] = *(Curr_El->get_coord());
//						nodescoord[8][1] = *(Curr_El->get_coord() + 1);
//
//						discharge->discharge_sens(nodescoord, dt, Curr_El->get_func_sens());
//
//					}
//					currentPtr = currentPtr->next;
//				}
//			}
//	}
//}

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

void sens_on_boundary(MeshCTX* meshctx, PropCTX* propctx, DualElem* eff_el, int side) {

	double* dx = eff_el->get_dx();

	TimeProps* timeprops = propctx->timeprops;
	int iter = timeprops->iter;
	double dt = timeprops->dt.at(iter - 1);	//at final time step we do not need the computation of adjoint and we always compute it for the previouse time so we need iter.
	FluxJac& flux_jac = eff_el->get_flx_jac_cont();
	double* func_sens = eff_el->get_func_sens();

	const int xp = eff_el->get_positive_x_side(); //finding the direction of element
	const int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	if (side == xm) {
		for (int k = 0; k < NUM_STATE_VARS; ++k)
			func_sens[k] += -1.5 * dt * dx[1] * flux_jac(0, 0, 0)(0, k);
	} else if (side == xp) {
		for (int k = 0; k < NUM_STATE_VARS; ++k)
			func_sens[k] += 1.5 * dt * dx[1] * flux_jac(0, 1, 0)(0, k);
	} else if (side == ym) {
		for (int k = 0; k < NUM_STATE_VARS; ++k)
			func_sens[k] += -1.5 * dt * dx[0] * flux_jac(1, 0, 0)(0, k);
	} else {
		for (int k = 0; k < NUM_STATE_VARS; ++k)
			func_sens[k] += 1.5 * dt * dx[0] * flux_jac(1, 1, 0)(0, k);
	}

}

double discharge = 0.;

void update_discharge(HashTable* El_Table, HashTable* NodeTable, DISCHARGE* dsge, double dt) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	double local_discharge = 0.;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Element* Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					for (int bound = 0; bound < 4; ++bound)
						if (Curr_El->get_neigh_proc()[bound] == INIT) {

							int xp = Curr_El->get_positive_x_side(); //finding the direction of element
							int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;
							if (bound == ym) {
								double flux[4][NUM_STATE_VARS];
								get_flux(El_Table, NodeTable, Curr_El, flux);

//								if (flux[3][0] < 0.)
								local_discharge += flux[3][0] * Curr_El->get_dx()[0] * dt;
							}
						}

//					unsigned *nodes = Curr_El->getNode();
//					double nodescoord[9][2], *coord;
//					Node* node;
//
//					for (int inode = 0; inode < 8; inode++) {
//						node = (Node*) NodeTable->lookup(nodes + 2 * inode);
//						coord = node->get_coord();
//
//						nodescoord[inode][0] = coord[0];
//						nodescoord[inode][1] = coord[1];
//
//					}
//					nodescoord[8][0] = *(Curr_El->get_coord());
//					nodescoord[8][1] = *(Curr_El->get_coord() + 1);
//
//					discharge->update(nodescoord, Curr_El->get_prev_state_vars(), dt);

				}

				currentPtr = currentPtr->next;
			}
		}

	double g_local_discharge = 0.;
	MPI_Reduce(&local_discharge, &g_local_discharge, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	discharge += g_local_discharge;

}

