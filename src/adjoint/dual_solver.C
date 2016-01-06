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
 * $Id: calc_jacobian.C 164 2013-02-27 15:27:22Z haghakhani $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"
#include <map>

#define DEBUG1

#define KEY0   3920807148
#define KEY1   1321528399
#define ITER   10
#define J      0

double initial_dot(HashTable* El_Table);

void copy_jacobian(HashTable* El_Table, map<int, Vec_Mat<9>>& jac_code);

void compare_jacobians(map<int, Vec_Mat<9>>& jac_code, map<int, Vec_Mat<9>>& jac_diff);

void clean_jacobian(HashTable* El_Table);

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx, PertElemInfo* eleminfo) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	const int rescomp = 1;
	const double increment = INCREMENT;
	const int maxiter = timeprops_ptr->iter;

	reset_adaption_flag(El_Table);

	double functional = 0.0, dt;

	cout << "computing ADJOINT time step " << maxiter << endl;

	set_ithm(El_Table);

	plot_ithm(El_Table);

	calc_adjoint(meshctx, propctx);

	int tecflag = 2;

//	uinform_refine(meshctx, propctx);

//	error_compute(meshctx, propctx, maxiter);

//	dual_unrefine(meshctx, propctx);

	meshplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., tecflag);

	setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr, timeprops_ptr);

	tecflag = 1;
	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		cout << "computing ADJOINT time step " << iter - 1 << endl;
		dt = timeprops_ptr->dt.at(iter - 1);
		timeprops_ptr->adjiter++;

		setup_dual_flow(solrec, meshctx, propctx);

//		cout << "test of adjoint: " << simple_test(El_Table) << endl;

		timeprops_ptr->adjoint_time(iter - 1);

//		compute_functional(El_Table, &functional, timeprops_ptr);

//		eleminfo->update_dual_func(functional);

		calc_jacobian(meshctx, propctx, eleminfo);

//		cout<<" max jac is: "<<max_jac<<endl;

		calc_adjoint(meshctx, propctx);

//		if (iter - 1 == 1)
		cout << "test of adjoint: " << simple_test(El_Table, timeprops_ptr, matprops_ptr) << endl;

//		map<int, Vec_Mat<9>> jac_code;

//		copy_jacobian(El_Table, jac_code);

//		calc_jacobian_old(meshctx, propctx);

//		map<int, Vec_Mat<9>> jac_diff;

//		copy_jacobian(El_Table, jac_diff);

//		compare_jacobians(jac_code, jac_diff);

//		clean_jacobian(El_Table);

//		print_Elem_Table(El_Table, NodeTable, timeprops_ptr->iter, 1);

//		check_state_vars_with_record(El_Table, solrec, iter);

//		print_jacobian(El_Table, iter);

//		if (eleminfo->iter == iter - 1)
//			fill_pertelem_info(El_Table, eleminfo);

//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
//sensitivity w.r.t to parameters

//		uinform_refine(meshctx, propctx);

//		error_compute(meshctx, propctx, iter);

// in dual weighted error estimation if solver performs n step, we'll have n+1
// solution and n+1 adjoint solution, but we'll have just n residual and as a
// result n error estimate. The point is that at initial step (0'th step),
// we know the solution from initial condition  so the error of 0th step is zero,
// and we have to compute the error for other time steps.

//		dual_unrefine(meshctx, propctx);
//		set_ithm(El_Table);
//
//		plot_ithm(El_Table);
//		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)
		meshplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., tecflag);

	}

	return;
}

int num_nonzero_elem(HashTable *El_Table) {
	int num = 0;			//myid
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

int num_nonzero_elem(HashTable *El_Table, int type) {
	int num = 0;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() == type)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

void reverse_states(HashTable* El_Table, HashTable* solrec, int iter, ElemPtrList* refinelist,
    ElemPtrList* unrefinelist) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

#ifdef DEBUG
	if (checkElement(El_Table))
	exit(22);
	for (int i = 0; i < nonz1; i++) {
		dbgvec[i] = 0;
		pass[i] = 0;
	}
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
		if (dbgvec[i] != 1) {
			cout << "these elements have problem:   " << i << endl
			<< "the value is  " << dbgvec[i] << endl;
			exit(EXIT_FAILURE);

		} else if (pass[i] != 1) {
			cout << "this index has not been passed  " << i << endl;
			exit(EXIT_FAILURE);

		}
	}

	delete[] dbgvec;

	cout << "number of elements after unrefinement"
	<< num_nonzero_elem(El_Table) << endl;
//		getchar();
	if (checkElement(El_Table))
	exit(22);
#endif

	int reg = 0, ref = 0, unref = 0;
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					Curr_El->rev_state_vars(solrec, El_Table, iter, &reg, &unref, &ref, refinelist,
					    unrefinelist);

				}
				currentPtr = currentPtr->next;
			}
		}
	}

//	cout << "reg= " << reg << "  ref=" << ref << "  unref= " << unref << endl;
}

void print_jacobian(HashTable* El_Table, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->print_jacobian(iter);

				currentPtr = currentPtr->next;
			}
		}
	}

}

void compute_functional(HashTable* El_Table, double* functional, TimeProps* timeprops_ptr) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	double const *dx;
	double const *state_vars;
	double const *prev_state_vars;
	double dt;

	dt = timeprops_ptr->dt.at(timeprops_ptr->iter - 1);

//	printf("iter=%4d  dt=%8f \n", timeprops_ptr->iter, dt);

//we do not have make it zero here, because we want to compute the integration over the time and space
//*functional = 0.0;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					dx = Curr_El->get_coord();
					state_vars = Curr_El->get_state_vars();
					prev_state_vars = Curr_El->get_prev_state_vars();

					// we used trapezoidal integration on time
					// flag is for time step 0

					*functional += .5
					    * (state_vars[0] * state_vars[0] + prev_state_vars[0] * prev_state_vars[0]) * dx[0]
					    * dx[1] * dt;

				}
				currentPtr = currentPtr->next;
			}
		}

}

void orgSourceSgn(Element* Curr_El, double frictiny, double* orgSgn) {

	double* d_state_vars_x = Curr_El->get_d_state_vars();
	double* d_state_vars_y = d_state_vars_x + NUM_STATE_VARS;
	double* prev_state_vars = Curr_El->get_prev_state_vars();
	double h_inv;
	double tmp = 0.0;
	double velocity[2];
	for (int i = 0; i < 2; i++)
		orgSgn[i] = 0.0;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		velocity[0] = prev_state_vars[1] / prev_state_vars[0];
		velocity[1] = prev_state_vars[2] / prev_state_vars[0];

	} else {
		for (int k = 0; k < DIMENSION; k++)
			velocity[k] = 0.;
	}

	if (prev_state_vars[0] > 0.0)
		h_inv = 1. / prev_state_vars[0];

	tmp = h_inv * (d_state_vars_y[1] - velocity[0] * d_state_vars_y[0]);
	orgSgn[0] = tiny_sgn(tmp, frictiny);

	tmp = h_inv * (d_state_vars_x[2] - velocity[1] * d_state_vars_x[0]);
	orgSgn[1] = tiny_sgn(tmp, frictiny);

}

void refinement_report(HashTable* El_Table) {

	int newbuffer = 0, buffer = 0, newson = 0, newfather = 0, norecadapt = 0, tobedeleted = 0,
	    oldfather = 0, oldson = 0;

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				switch (Curr_El->get_adapted_flag()) {
					case 5:
						newbuffer++;
						break;
					case 4:
						buffer++;
						break;
					case 3:
						newson++;
						break;
					case 2:
						newfather++;
						break;
					case 1:
						norecadapt++;
						break;
					case 0:
						tobedeleted++;
						break;
					case -6:
						oldfather++;
						break;
					case -7:
						oldson++;
						break;
					default:
						cout << "this case is irregular \n";
				}
				currentPtr = currentPtr->next;
			}
		}

	cout << " new buffer: " << newbuffer << "\n buffer:     " << buffer << "\n newson:     " << newson
	    << "\n newfather:  " << newfather << "\n norecadapt: " << norecadapt << "\n tobedeleted:"
	    << tobedeleted << "\n oldfather:  " << oldfather << "\n oldson:     " << oldson << "\n";

}

void dual_refine_unrefine(MeshCTX* meshctx, PropCTX* propctx, ElemPtrList* refinelist,
    ElemPtrList* unrefinelist) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	ElemPtrList NewFatherList, OtherProcUpdate;

	int rescomp = 0;

//	unsigned keyy[2] = { 3410598297, 2576980374 };

	delete_unused_elements_nodes(El_Table, NodeTable, myid);

	double target = 0.1;

//	cout << "1 \n";
//	refinement_report(El_Table);

	if (refinelist->get_num_elem()) {
		// start refinement
		for (int i = 0; i < refinelist->get_num_elem(); ++i) {
			refine(refinelist->get(i), El_Table, NodeTable, matprops_ptr, rescomp);
			(refinelist->get(i))->put_adapted_flag(OLDFATHER);
			(refinelist->get(i))->put_refined_flag(1);
		}

//		cout << "2 \n";
//		refinement_report(El_Table);

		refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) refinelist, timeprops_ptr);

//		cout << "3 \n";
//		refinement_report(El_Table);

		move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);

//		cout << "4 \n";
//		refinement_report(El_Table);

		int refdel = 0;
		int hash_size = El_Table->get_no_of_buckets();
		for (int i = 0; i < hash_size; i++) {
			HashEntryPtr entryp = *(El_Table->getbucketptr() + i);
			while (entryp) {
				Element* EmTemp = (Element*) (entryp->value);
				entryp = entryp->next;

				if (EmTemp->get_adapted_flag() == TOBEDELETED) {
					El_Table->remove(EmTemp->pass_key());
					refdel++;

				}
			}
		}
//		cout << "5 \n";
//		refinement_report(El_Table);
//		cout << "number of deleted elem after ref " << refdel << "  number of ref list  "
//		    << refinelist->get_num_elem() << endl;
	}

	if (unrefinelist->get_num_elem()) {

		Element* brothers[4];

//		cout << "6 \n";
//		refinement_report(El_Table);

		int unrefined = 0;

		do {

			NewFatherList.trashlist();

			for (int i = 0; i < unrefinelist->get_num_elem(); ++i) {
				Element* Curr_El = (unrefinelist->get(i));
				if ((Curr_El->get_which_son() == 0) && (Curr_El->get_adapted_flag() != OLDSON))

					Curr_El->find_brothers(El_Table, NodeTable, target, myid, matprops_ptr, &NewFatherList,
					    &OtherProcUpdate, rescomp);

			}

			unrefined += NewFatherList.get_num_elem();

//			cout << "7 \n";
//			refinement_report(El_Table);

			unrefine_neigh_update(El_Table, NodeTable, myid, (void*) &NewFatherList);

//			cout << "8 \n";
//			refinement_report(El_Table);

			unrefine_interp_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &OtherProcUpdate);

			for (int k = 0; k < NewFatherList.get_num_elem(); k++)
				delete_oldsons(El_Table, NodeTable, myid, NewFatherList.get(k));

			reset_adaption_flag(El_Table);
//			reset_newson_adaption_flag(El_Table);

//			cout<<"this is the counter  "<<counter++<<endl;

		} while (unrefined != (unrefinelist->get_num_elem() / 4));

//		cout << "9 \n";
//		refinement_report(El_Table);
//
//		cout << "number of created elem after unref " << unrefined << "  number of unref list  "
//		    << unrefinelist->get_num_elem() << endl;

		move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);
	}

//	calc_d_gravity(El_Table);

	refinelist->trashlist();
	unrefinelist->trashlist();
	reset_adaption_flag(El_Table);

}

void reset_adaption_flag(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
					Curr_El->put_refined_flag(0);
				}

				currentPtr = currentPtr->next;
			}
		}

}

void reset_newson_adaption_flag(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() == NEWSON) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
					Curr_El->put_refined_flag(0);
				}

				currentPtr = currentPtr->next;
			}
		}

}

void setup_dual_flow(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;

	ElemPtrList refinelist, unrefinelist;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	if (solrec->get_first_solution() >= iter) {
		solrec->free_all_available_sol();
		solrec->load_new_set_of_solution();
	}

	if (iter % 5 == 4 && propctx->adapt_flag != 0) {

		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);

					if (Curr_El->get_adapted_flag() > 0)
						Curr_El->check_refine_unrefine(solrec, El_Table, iter, &refinelist, &unrefinelist);

					currentPtr = currentPtr->next;
				}
			}

		dual_refine_unrefine(meshctx, propctx, &refinelist, &unrefinelist);

//			setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr, timeprops_ptr);

		calc_d_gravity(El_Table);

	}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->update_state(solrec, El_Table, iter);

				currentPtr = currentPtr->next;
			}
		}

//	allocJacoMat(El_Table);

	int num_node_buckets = NodeTable->get_no_of_buckets();
	buck = NodeTable->getbucketptr();
	for (int i = 0; i < num_node_buckets; i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Node* Curr_Node = (Node*) (currentPtr->value);
				Curr_Node->zero_flux();

				currentPtr = currentPtr->next;
			}
		}

	// this function computes fluxes based on prev_state_vars (we need for dual problem),
	// and jacobian of fluxes and store the in elements
	calc_flux(meshctx, propctx);

	//this function computes slopes based on prev_state_vars and dh/dh_e where h_e is pile height in neighbor element
	// we need this term to compute jacobian of elements
	slopes(El_Table, NodeTable, matprops_ptr, 1);

}

void check_state_vars_with_record(HashTable* El_Table, HashTable* solrec, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					if (Curr_El->check_state(solrec, El_Table, iter))
						cout << "the idea did not work \n";

				currentPtr = currentPtr->next;
			}
		}

}

void print_Elem_Table(HashTable* El_Table, HashTable* NodeTable, int iter, int place) {

	ofstream myfile;
	char filename[50];
	sprintf(filename, "El_Table_%d_%08d", place, iter);

	myfile.open(filename, ios::app);

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					int print = 1;
					for (int j = 0; j < 4; j++)
						if (*(Curr_El->get_neigh_proc() + j) == INIT) {
							print = 0;
							break;
						}
					if (print) {

						int xp = Curr_El->get_positive_x_side(); //finding the direction of element
						int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

						Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

						Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

						Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

						Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

						double flux[4][NUM_STATE_VARS];

						for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
							flux[0][ivar] = *(nxp->get_flux() + ivar);
							flux[1][ivar] = *(nyp->get_flux() + ivar);
							flux[2][ivar] = *(nxm->get_flux() + ivar);
							flux[3][ivar] = *(nym->get_flux() + ivar);
						}

						double flux_diff[DIMENSION][NUM_STATE_VARS];

						for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
							flux_diff[0][ivar] = flux[0][ivar] - flux[2][ivar];
							flux_diff[1][ivar] = flux[1][ivar] - flux[3][ivar];
						}

						myfile << "key: ";
						myfile << *(Curr_El->pass_key()) << " " << *(Curr_El->pass_key() + 1) << " ";

						myfile << "neighbors: ";
						for (int i = 0; i < 8 * 2; ++i)
							myfile << *(Curr_El->get_neighbors() + i) << " ";

						myfile << "neighbors_gen: ";
						for (int i = 0; i < 8; ++i)
							myfile << *(Curr_El->get_neigh_gen() + i) << " ";

						myfile << "neighbors_proc: ";
						for (int i = 0; i < 8; ++i)
							myfile << *(Curr_El->get_neigh_proc() + i) << " ";

						myfile << "state_vars: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << *(Curr_El->get_state_vars() + i) << " ";

						myfile << "prev_state_vars: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << *(Curr_El->get_prev_state_vars() + i) << " ";

						myfile << "d_state_vars: ";
						for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
							myfile << *(Curr_El->get_d_state_vars() + i) << " ";

//					myfile << "flux_diff_x: ";
//					for (int i = 0; i < NUM_STATE_VARS; ++i)
//						myfile << flux_diff[0][i] << " ";
//
//					myfile << "flux_diffy: ";
//					for (int i = 0; i < NUM_STATE_VARS; ++i)
//						myfile << flux_diff[1][i] << " ";
//
						myfile << "flux_xp: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[0][i] << " ";

						myfile << "flux_yp: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[1][i] << " ";

						myfile << "flux_xm: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[2][i] << " ";

						myfile << "flux_ym: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[3][i] << " ";

						myfile << "elem_id: " << Curr_El->get_ithelem() << " ";

//				myfile <<"gravity: ";
//				for (int i = 0; i < NUM_STATE_VARS; ++i)
//					myfile << *(Curr_El->get_gravity() + i) << " ";
//
//				myfile <<"d_gravity: ";
//				for (int i = 0; i < DIMENSION; ++i)
//					myfile << *(Curr_El->get_d_gravity() + i) << " ";

						myfile << endl;
					}
				}

				currentPtr = currentPtr->next;
			}
		}

	myfile.close();

}

bool must_write(MemUse* memuse_ptr) {

	struct sysinfo memInfo;
	sysinfo(&memInfo);
	unsigned long totalPhysMem = memInfo.totalram;
	unsigned long freeram = memInfo.freeram;
	unsigned long current_physMemUsed = memInfo.totalram - memInfo.freeram;
	long long last_time_step_use = current_physMemUsed - memuse_ptr->usedmem;

//	if (memuse_ptr->usedmem > 0) {
//		memuse_ptr->usedmem = current_physMemUsed;
//
//		double ratio = ((long double) last_time_step_use) / ((long double) memInfo.freeram);
	double ratio1 = ((long double) current_physMemUsed) / ((long double) totalPhysMem);
	double ratio2 = ((long double) freeram) / ((long double) totalPhysMem);
//		cout << "ratio of last time use over free mem " << ratio << "  and ratio od used mem" << endl;
	printf(" ratio of used mem over total mem is %4f \n", ratio1);

	if (ratio2 < .05)
		return true;
//	} else
//		memuse_ptr->usedmem = current_physMemUsed;

	return false;

}

void dual_unrefine(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	vector<RefUnref> unref_list;

	ElemPtrList NewFatherList, OtherProcUpdate;

	HashEntryPtr currentPtr;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Element *Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() == NEWSON && Curr_El->get_which_son() == 0) {
					RefUnref sample(Curr_El->pass_key());
					unref_list.push_back(sample);
				}

				currentPtr = currentPtr->next;
			}
		}

	int unrefined = 0;
	int rescomp = 1;
	double target = 0.1;

	int to_be_unrefined = unref_list.size();

	while (unrefined != to_be_unrefined) {

		NewFatherList.trashlist();
		OtherProcUpdate.trashlist();

		for (int i = 0; i < to_be_unrefined; ++i) {
			Element* Curr_El = (Element*) El_Table->lookup(unref_list[i].key);
			if (Curr_El && Curr_El->get_adapted_flag() != OLDSON)

				Curr_El->find_brothers(El_Table, NodeTable, target, myid, matprops_ptr, &NewFatherList,
				    &OtherProcUpdate, rescomp);

		}

		unrefined += NewFatherList.get_num_elem();
//		cout << unrefined << endl;

//			cout << "7 \n";
//			refinement_report(El_Table);

		unrefine_neigh_update(El_Table, NodeTable, myid, (void*) &NewFatherList);

//			cout << "8 \n";
//			refinement_report(El_Table);

		unrefine_interp_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &OtherProcUpdate);

		for (int k = 0; k < NewFatherList.get_num_elem(); k++)
			delete_oldsons(El_Table, NodeTable, myid, NewFatherList.get(k));

//		refinement_report(El_Table);

		reset_adaption_flag(El_Table);
//			reset_newson_adaption_flag(El_Table);

//			cout<<"this is the counter  "<<counter++<<endl;

	}

}

double simple_test(HashTable* El_Table, TimeProps* timeprops, MatProps* matprops_ptr) {

	double dot = 0.;
	vector<int> wrong_elem;
	vector<double> wrong_value, wrong_value1;

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	cout << "results:  \n";

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					double* state_vars_prev = Curr_El->get_prev_state_vars();
					double* adjoint = Curr_El->get_prev_adjoint();
					double* gravity = Curr_El->get_gravity();
					double* curve = Curr_El->get_curvature();
					double vel[2], h_inv, orgSrcSgn[2];
					double unitvx, unitvy;
					double* d_state_vars_x = Curr_El->get_d_state_vars();
					double* d_state_vars_y = (Curr_El->get_d_state_vars() + 3);
					double fric_tiny = matprops_ptr->frict_tiny;
					double kact = *(Curr_El->get_kactxy());
					double* d_grav = Curr_El->get_d_gravity();

					if (state_vars_prev[0] > GEOFLOW_TINY) {
						h_inv = 1. / state_vars_prev[0];
						vel[0] = state_vars_prev[1] / state_vars_prev[0];
						vel[1] = state_vars_prev[2] / state_vars_prev[0];

					} else {
						h_inv = 0.;
						vel[0] = 0.;
						vel[1] = 0.;
						unitvx = unitvy = 0.;

					}

					double speed = sqrt(vel[0] * vel[0] + vel[1] * vel[1]);

					if (speed > 0.) {

						unitvx = vel[0] / speed;
						unitvy = vel[1] / speed;
					}

					double test1 = 0.;

					if (max(gravity[2] * state_vars_prev[0] + vel[0] * state_vars_prev[1] * curve[0], 0.))
						test1 = adjoint[1] * unitvx
						    * (state_vars_prev[0] * gravity[2] + state_vars_prev[1] * vel[0] * curve[0]);

					if (max(gravity[2] * state_vars_prev[0] + vel[1] * state_vars_prev[2] * curve[1], 0.))
						test1 += adjoint[2] * unitvy
						    * (state_vars_prev[0] * gravity[2] + state_vars_prev[2] * vel[1] * curve[1]);

					dot += test1;

					double tmp = h_inv * (d_state_vars_y[1] - vel[0] * d_state_vars_y[0]);
					orgSrcSgn[0] = tiny_sgn(tmp, fric_tiny);

					tmp = h_inv * (d_state_vars_x[2] - vel[1] * d_state_vars_x[0]);
					orgSrcSgn[1] = tiny_sgn(tmp, fric_tiny);

					double test2 = state_vars_prev[0] * kact
					    * (adjoint[1] * orgSrcSgn[0]
					        * (gravity[2] * d_state_vars_y[0] + d_grav[1] * state_vars_prev[0])
					        + adjoint[2] * orgSrcSgn[1]
					            * (gravity[2] * d_state_vars_x[0] + d_grav[0] * state_vars_prev[0]));

					// this is for simple test
//					double test1, test2 = 0.;
//					test1 = adjoint[1] + adjoint[2] ;

					if (fabs(test1) > 1e-16 || fabs(test2) > 1e-16) {
						wrong_elem.push_back(Curr_El->get_ithelem());
						wrong_value.push_back(test1);
						wrong_value1.push_back(test2);
					}
					cout << test1 << " , " << test2 << endl;

				}

				currentPtr = currentPtr->next;
			}
		}

	ofstream f("wrong_elem.txt", ios::app);
	f << "time step: " << timeprops->iter << " size of vector: " << wrong_elem.size() << endl;
	for (int i = 0; i < wrong_elem.size(); ++i)
		f << wrong_elem[i] << " , " << wrong_value[i] << " , " << wrong_value1[i] << '\n';

	return dot;
}

double initial_dot(HashTable* El_Table) {

	double dot = 0.;

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0 && *(Curr_El->get_state_vars()) > dot)
					dot = *(Curr_El->get_state_vars());

				currentPtr = currentPtr->next;
			}
		}

	return dot;
}

void set_ithm(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	int count = 0;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					Curr_El->put_ithelem(count);
					count++;
				}
			}
		}
}

void plot_ithm(HashTable* El_Table) {
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					*(Curr_El->get_residual()) = Curr_El->get_ithelem();
				}
			}
		}
}

void copy_jacobian(HashTable* El_Table, map<int, Vec_Mat<9>>& jac_map) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					jac_map[Curr_El->get_ithelem()] = Curr_El->get_jacobian();
				}
			}
		}

}

void compare_jacobians(map<int, Vec_Mat<9>>& jac_code, map<int, Vec_Mat<9>>& jac_diff) {

	for (map<int, Vec_Mat<9>>::iterator it = jac_code.begin(); it != jac_code.end(); ++it) {
		Vec_Mat<9>& jacdiff = jac_diff[it->first];
		Vec_Mat<9>& jaccode = jac_code[it->first];

		for (int i = 0; i < EFF_ELL; ++i)
			for (int j = 0; j < NUM_STATE_VARS; ++j)
				for (int k = 0; k < NUM_STATE_VARS; ++k)
					if (fabs(jacdiff(i, j, k) - jaccode(i, j, k)) > 1e-11 && jacdiff(i, j, k) != 0.) {
						double denum;
						if (jacdiff(i, j, k) != 0.)
							denum = dabs(jacdiff(i, j, k));
						cout << "in element  " << it->first << " in indices: " << i << " , " << j << " , " << k
						    << "  there is a difference of  "
						    << 100 * dabs(jacdiff(i, j, k) - jaccode(i, j, k)) / denum << "  jacobian diff= "
						    << jacdiff(i, j, k) << "  jacobian code=" << jaccode(i, j, k) << endl;
					}

	}

}

void clean_jacobian(HashTable* El_Table) {
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					Vec_Mat<9>& jacobian = Curr_El->get_jacobian();
					for (int i = 0; i < EFF_ELL; ++i)
						jacobian(i) = ZERO_MATRIX;
				}
			}
		}

}
