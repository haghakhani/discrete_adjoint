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
#include <algorithm>
#include <map>

double max_err1 = 0., max_err2 = 0.;
unsigned key1_1, key2_1, iter_1, key1_2, key2_2, iter_2;

#define DEBUG1
//#define Error

#define KEY0   3920807148
#define KEY1   1321528399
#define ITER   10
#define J      0

void copy_jacobian(HashTable* El_Table, map<int, Vec_Mat<9>>& jac_code);

void compare_jacobians(map<int, Vec_Mat<9>>& jac_code, map<int, Vec_Mat<9>>& jac_diff);

bool is_old_son(Element* elem) {
	return (elem->get_adapted_flag() == OLDSON);
}
;

void check_kackxy(HashTable* Err_El_Tab);

template<typename T1, typename T2>
void copy_hashtables_objects(HashTable* El_Table, HashTable* cp_El_Table) {

	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				T1 *Curr = (T1*) (currentPtr->value);
				T2 *cp = new T2(Curr);
				cp_El_Table->add(cp->pass_key(), cp);
				currentPtr = currentPtr->next;

			}
		}
}

template<typename T1>
void delete_hashtables_objects(HashTable* El_Table) {
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				T1 *Curr = (T1*) (currentPtr->value);
				delete Curr;
				currentPtr = currentPtr->next;

			}
		}

	delete El_Table;
}

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	const int maxiter = timeprops_ptr->iter;

	cout << endl << "******************************\n" << "*          start             *\n"
	    << "*     ADJOINT SOLUTION       *\n" << "*                            *\n"
	    << "*                            *\n" << "******************************\n"
			    "Generating grids for computing the Adjoint and the Error ....\n";

	HashTable *Dual_El_Tab = new HashTable(El_Table);

	MeshCTX dual_meshctx;
	dual_meshctx.el_table = Dual_El_Tab;
	dual_meshctx.nd_table = NodeTable;

	copy_hashtables_objects<Element, DualElem>(El_Table, Dual_El_Tab);

	cout << "The Adjoint grid has been generated ....\n";

	calc_adjoint(&dual_meshctx, propctx);

	dualplotter(Dual_El_Tab, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., 2);

#ifdef Error

	HashTable *Err_El_Tab = new HashTable(El_Table);
	HashTable *Err_Nod_Tab = new HashTable(NodeTable);

	copy_hashtables_objects<Element, ErrorElem>(El_Table, Err_El_Tab);
	copy_hashtables_objects<Node, Node>(NodeTable, Err_Nod_Tab);

	MeshCTX error_meshctx;
	error_meshctx.el_table = Err_El_Tab;
	error_meshctx.nd_table = Err_Nod_Tab;

#endif

	delete_hashtables_objects<Element>(El_Table);

	cout << "The Error grid has been generated ....\n";

	cout << "computing ADJOINT time step " << maxiter << endl;

#ifdef Error

	// this function refines and do constant reconstruction
	uinform_refine(&error_meshctx, propctx);

	make_dual_err_link(Dual_El_Tab, Err_El_Tab);

	send_from_dual_to_error(Dual_El_Tab, Err_El_Tab);

	bilinear_interp(Err_El_Tab);

	update_bilinear_error_grid(&error_meshctx, propctx);

	error_compute(&error_meshctx, propctx);

	errorplotter(Err_El_Tab, Err_Nod_Tab, matprops_ptr, timeprops_ptr, mapname_ptr, 0.);
#else
	MeshCTX error_meshctx;
#endif

//	set_ithm(El_Table);
//
//	plot_ithm(El_Table);

//	cout << "elements number of original grid " << num_nonzero_elem(El_Table) << endl;
//	cout << "elements number of refined grid " << num_nonzero_elem(cp_El_Table) << endl;

//	print_Elem_Table(El_Table, NodeTable, timeprops_ptr->iter, 0);
//
//	print_Elem_Table(cp_El_Table, cp_NodeTable, timeprops_ptr->iter, 1);

//	int tecflag = 2;

//	set_ithm(El_Table);
//	plot_ithm(El_Table);

//	reset_adaption_flag(Dual_El_Tab);

//	refinement_report(El_Table);
//
//	refinement_report(cp_El_Table);

//	this function reconstruct bilinear interpolation
//	set_ithm(cp_El_Table);
//	plot_ithm(cp_El_Table);
//	cout<<num_nonzero_elem(El_Table)<<endl;
//	cout<<num_nonzero_elem(cp_El_Table)<<endl;
//
//	timeprops_ptr->iter--;
//	meshplotter(cp_El_Table, cp_NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., tecflag);

//	setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr, timeprops_ptr);

	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		cout << "computing ADJOINT time step " << iter - 1 << endl;
		timeprops_ptr->adjiter++;

		setup_dual_flow(solrec, &dual_meshctx, &error_meshctx, propctx);

//		cout << "elements number of original grid " << num_nonzero_elem(El_Table) << endl;
//		cout << "elements number of refined grid " << num_nonzero_elem(cp_El_Table) << endl;

		timeprops_ptr->adjoint_time(iter - 1);

//		compute_functional(El_Table, &functional, timeprops_ptr);

//		eleminfo->update_dual_func(functional);

		calc_jacobian(&dual_meshctx, propctx);

//		cout<<" max jac is: "<<max_jac<<endl;

		calc_adjoint(&dual_meshctx, propctx);

//		cout << "test of adjoint: " << simple_test(Dual_El_Tab, timeprops_ptr, matprops_ptr) << endl;

#ifdef Error

		send_from_dual_to_error(Dual_El_Tab, Err_El_Tab);

		bilinear_interp(Err_El_Tab);

		update_bilinear_error_grid(&error_meshctx, propctx);

		error_compute(&error_meshctx, propctx);
		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)

		errorplotter(Err_El_Tab, Err_Nod_Tab, matprops_ptr, timeprops_ptr, mapname_ptr, 0.);
#endif

//		if (iter - 1 == 1)
//		cout << "test of adjoint: " << simple_test(El_Table, timeprops_ptr, matprops_ptr) << endl;

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

//		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)
		dualplotter(Dual_El_Tab, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., 1);

	}

	delete_hashtables_objects<DualElem>(Dual_El_Tab);
	delete_hashtables_objects<Node>(NodeTable);

#ifdef Error
	delete_hashtables_objects<DualElem>(Err_El_Tab);
	delete_hashtables_objects<Node>(Err_Nod_Tab);
#endif

	delete_hashtables_objects<Jacobian>(solrec);

//	if (fabs(max_err1) > fabs(max_err2))
//		cout << "max error occurred in test 1" << max_err1 << "  at iter " << iter_1 << " key is "
//		    << key1_1 << " , " << key2_1 << endl;
//	else
//		cout << "max error occurred in test 2" << max_err2 << "  at iter " << iter_2 << " key is "
//		    << key1_2 << " , " << key2_2 << endl;

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

void print_jacobian(HashTable* El_Table, int iter) {

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
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

template<typename T>
void dual_refine_unrefine(MeshCTX* meshctx, PropCTX* propctx, ElemPtrList<T>* refinelist,
    ElemPtrList<T>* unrefinelist) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	ElemPtrList<Element> NewFatherList, OtherProcUpdate;

	vector<T*> first_son;
	vector<pair<unsigned, unsigned> > new_father;

//	unsigned keyy[2] = { 3410598297, 2576980374 };

	reset_adaption_flag(El_Table);

	delete_unused_elements_nodes(El_Table, NodeTable, myid);

	double target = 0.1;

//	cout << "1 \n";
//	refinement_report(El_Table);

	if (refinelist->get_num_elem()) {
		// start refinement
		for (int i = 0; i < refinelist->get_num_elem(); ++i) {
			refine(refinelist->get(i), El_Table, NodeTable, matprops_ptr, 1);
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

		for (int i = 0; i < unrefinelist->get_num_elem(); ++i)
			if (unrefinelist->get(i)->get_which_son() == 0) {
				first_son.push_back(unrefinelist->get(i));

				//we need the the new fathers for future
				//this is not required but is good for checking purposes
//				new_father.push_back(
//				    make_pair(*(unrefinelist->get(i)->getfather()),
//				        *(unrefinelist->get(i)->getfather() + 1)));
			}

		//		cout << "6 \n";
		//		refinement_report(El_Table);

		while (first_son.size() > 0) {

			for (int i = 0; i < first_son.size(); ++i)
				first_son[i]->template find_brothers<T>(El_Table, NodeTable, target, myid, matprops_ptr,
				    &NewFatherList, &OtherProcUpdate, 1);

			unrefine_neigh_update(El_Table, NodeTable, myid, (void*) &NewFatherList);

			unrefine_interp_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &OtherProcUpdate);

			first_son.erase(remove_if(first_son.begin(), first_son.end(), is_old_son), first_son.end());

			for (int k = 0; k < NewFatherList.get_num_elem(); k++)
				delete_oldsons(El_Table, NodeTable, myid, NewFatherList.get(k));

			reset_newfather_adaption_flag(El_Table);

			NewFatherList.trashlist();

		}

		move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);
	}

//	calc_d_gravity(El_Table);
//	cout << "8 \n";
	//this is not required but good for checks
//	set_new_fathers(El_Table, new_father);
//	refinement_report(El_Table);

	refinelist->trashlist();
	unrefinelist->trashlist();
//	reset_adaption_flag(El_Table);

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

void reset_newfather_adaption_flag(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() == NEWFATHER) {
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

void setup_dual_flow(SolRec* solrec, MeshCTX* dual_meshctx, MeshCTX* err_meshctx,
    PropCTX* propctx) {

	HashTable* El_Table = dual_meshctx->el_table;
	HashTable* NodeTable = dual_meshctx->nd_table;
#ifdef Error
	HashTable* cp_El_Table = err_meshctx->el_table;
	HashTable* cp_NodeTable = err_meshctx->nd_table;
#endif
	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;

	if (solrec->get_first_solution() >= iter - 1) {
//		solrec->free_all_available_sol();
		solrec->load_new_set_of_solution();
	}

	HashEntryPtr *buck;
#ifdef Error
	if (iter % 5 == 0 && propctx->adapt_flag != 0) {

		ElemPtrList<ErrorElem> refinelist, unrefinelist;

		buck = cp_El_Table->getbucketptr();

		for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				ErrorElem *Curr_El = (ErrorElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
				Curr_El->error_check_refine_unrefine(solrec, cp_El_Table, iter, &refinelist,
						&unrefinelist);

				currentPtr = currentPtr->next;
			}
		}
//		cout<<"in refined table"<<endl;
//		cout<<"has to be refined "<<refinelist->get_num_elem()<<endl;
//		cout<<"has to be unrefined "<<unrefinelist->get_num_elem()<<endl;

		dual_refine_unrefine<ErrorElem>(err_meshctx, propctx, &refinelist, &unrefinelist);

		calc_d_gravity(cp_El_Table);

	}

	if (iter != 1)
	update_error_grid(solrec, err_meshctx, propctx);
#endif
	if (iter % 5 == 4 && propctx->adapt_flag != 0) {

		ElemPtrList<DualElem> refinelist, unrefinelist;

		buck = El_Table->getbucketptr();
		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				HashEntryPtr currentPtr = *(buck + i);
				while (currentPtr) {
					DualElem *Curr_El = (DualElem*) (currentPtr->value);

					if (Curr_El->get_adapted_flag() > 0)
						Curr_El->dual_check_refine_unrefine(solrec, El_Table, iter, &refinelist, &unrefinelist);

					currentPtr = currentPtr->next;
				}
			}

//		cout<<"in original table"<<endl;
//		cout<<"has to be refined "<<refinelist->get_num_elem()<<endl;
//		cout<<"has to be unrefined "<<unrefinelist->get_num_elem()<<endl;

		dual_refine_unrefine<DualElem>(dual_meshctx, propctx, &refinelist, &unrefinelist);

//		set_ithm(El_Table);
//
//		plot_ithm(El_Table);

		calc_d_gravity(El_Table);

	}

	// inside updating state_vars we also delete the solution that we used
	// since we no longer need it
	update_dual_grid(solrec, dual_meshctx, propctx);

	clear_empty_jacobians(solrec, iter);

}

void check_state_vars_with_record(HashTable* El_Table, SolRec* solrec, int iter) {

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

	if (ratio1 < .05)
		return true;
//	} else
//		memuse_ptr->usedmem = current_physMemUsed;

	return false;

}

double simple_test(HashTable* El_Table, TimeProps* timeprops, MatProps* matprops_ptr) {

	double dot = 0.;
	vector<pair<unsigned, unsigned>> wrong_elem;
	vector<double> wrong_value, wrong_value1;

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	cout << "results:  \n";

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					double* state_vars_prev = Curr_El->get_prev_state_vars();
					double* adjoint = Curr_El->get_adjoint();
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
					orgSrcSgn[2] = tiny_sgn(vel[0], fric_tiny);

					tmp = h_inv * (d_state_vars_x[2] - vel[1] * d_state_vars_x[0]);
					orgSrcSgn[1] = tiny_sgn(tmp, fric_tiny);
					orgSrcSgn[3] = tiny_sgn(vel[1], fric_tiny);

					double test2 = state_vars_prev[0] * kact
					    * (adjoint[1] * orgSrcSgn[0]
					        * (gravity[2] * d_state_vars_y[0] + d_grav[1] * state_vars_prev[0])
					        + adjoint[2] * orgSrcSgn[1]
					            * (gravity[2] * d_state_vars_x[0] + d_grav[0] * state_vars_prev[0]));

					// this is for simple test
//					double test1, test2 = 0.;
//					test1 = adjoint[1] + adjoint[2] ;

					if (fabs(test1) > 1e-16 || fabs(test2) > 1e-16) {
						wrong_elem.push_back(make_pair(*(Curr_El->pass_key()), *(Curr_El->pass_key() + 1)));
						wrong_value.push_back(test1);
						wrong_value1.push_back(test2);
//						cout << test1 << " , " << test2 << endl;

						if (fabs(test1) > max_err1) {
							max_err1 = fabs(test1);
							key1_1 = *(Curr_El->pass_key());
							key2_1 = *(Curr_El->pass_key() + 1);
							iter_1 = timeprops->iter;
						} else if (fabs(test2) > max_err2) {
							max_err2 = fabs(test2);
							key1_2 = *(Curr_El->pass_key());
							key2_2 = *(Curr_El->pass_key() + 1);
							iter_2 = timeprops->iter;
						}

					}

				}

				currentPtr = currentPtr->next;
			}
		}

	ofstream f("wrong_elem.txt", ios::app);
	f << "time step: " << timeprops->iter << " size of vector: " << wrong_elem.size() << endl;
	for (int i = 0; i < wrong_elem.size(); ++i)
		f << setw(10) << wrong_elem[i].first << " , " << wrong_elem[i].second << " , " << wrong_value[i]
		    << " , " << wrong_value1[i] << '\n';

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
	ErrorElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (ErrorElem*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					*(Curr_El->get_residual()) = Curr_El->get_ithelem();
				}
			}
		}
}

void copy_jacobian(HashTable* El_Table, map<int, Vec_Mat<9>>& jac_map) {

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
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
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					Vec_Mat<9>& jacobian = Curr_El->get_jacobian();
					for (int i = 0; i < EFF_ELL; ++i)
						jacobian(i) = ZERO_MATRIX;
				}
			}
		}

}

void copy_hashtables(HashTable* El_Table, HashTable* NodeTable, HashTable* cp_El_Table,
    HashTable* cp_NodeTable) {

	Node *Curr_node;
	HashEntryPtr *buck = NodeTable->getbucketptr();

	for (int i = 0; i < NodeTable->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_node = (Node*) (currentPtr->value);
				Node *cp_node = new Node(Curr_node);
				cp_NodeTable->add(cp_node->pass_key(), cp_node);
				currentPtr = currentPtr->next;

			}
		}

	Element *Curr_El;
	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				DualElem *cp_elem = new DualElem(Curr_El);
				cp_El_Table->add(cp_elem->pass_key(), cp_elem);
				currentPtr = currentPtr->next;

			}
		}

}

void update_error_grid(SolRec* solrec, MeshCTX* cp_meshctx, PropCTX* propctx) {

	HashTable* cp_El_Table = cp_meshctx->el_table;
	HashTable* cp_NodeTable = cp_meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;

	HashEntryPtr currentPtr;
	ErrorElem *Curr_El;

	HashEntryPtr *buck = cp_El_Table->getbucketptr();
	for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (ErrorElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->error_update_state(solrec, iter);

				currentPtr = currentPtr->next;
			}
		}

}

void set_new_fathers(HashTable* El_Table, vector<pair<unsigned, unsigned> >& new_father) {

	HashEntryPtr currentPtr;
	Element *Curr_El;

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() == NOTRECADAPTED) {

					for (int j = 0; j < new_father.size(); ++j) {
						unsigned key[] = { new_father[j].first, new_father[j].second };
						if (compare_key(Curr_El->pass_key(), key))
							Curr_El->put_adapted_flag(NEWFATHER);
					}
				}
				currentPtr = currentPtr->next;
			}
		}
}

void check_kackxy(HashTable* El_Tab) {

	HashEntryPtr currentPtr;
	Element *Curr_El;

	int aa = 1, bb = 0;

	HashEntryPtr *buck = El_Tab->getbucketptr();
	for (int i = 0; i < El_Tab->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					if (*(Curr_El->get_kactxy()) < 0 || *(Curr_El->get_kactxy() + 1) < 0)
						aa = bb;

				currentPtr = currentPtr->next;
			}
		}

}

void make_dual_err_link(HashTable *Dual_El_Tab, HashTable *Err_El_Tab) {

	HashEntryPtr *buck = Err_El_Tab->getbucketptr();
	HashEntryPtr currentPtr;

	for (int i = 0; i < Err_El_Tab->get_no_of_buckets(); ++i)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				ErrorElem *son[4];

				son[0] = (ErrorElem*) (currentPtr->value);

				son[0]->calc_which_son();

				if (son[0]->get_adapted_flag() > 0 && son[0]->get_which_son() == 0) {

					unsigned* father_key = son[0]->getfather();

					DualElem* father = (DualElem*) Dual_El_Tab->lookup(father_key);

					vector<ErrorElem*>& mysons = father->get_son_addresses();

					for (int j = 1; j < 4; ++j)
						son[j] = (ErrorElem*) Err_El_Tab->lookup(son[0]->get_brothers() + j * KEYLENGTH);

					for (int j = 0; j < 4; ++j) {
						mysons.push_back(son[j]);
						son[j]->put_father_address(father);
					}
				}
				currentPtr = currentPtr->next;
			}
		}
}

void send_from_dual_to_error(HashTable *Dual_El_Tab, HashTable *Err_El_Tab) {

	HashEntryPtr *buck = Dual_El_Tab->getbucketptr();
	HashEntryPtr currentPtr;

	for (int i = 0; i < Dual_El_Tab->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				DualElem* dual_el = (DualElem*) (currentPtr->value);

				if (dual_el->get_adapted_flag() > 0) {

					vector<ErrorElem*> error_el = dual_el->get_son_addresses();

					for (int j = 0; j < error_el.size(); ++j)
						for (int k = 0; k < NUM_STATE_VARS; ++k) {

							*(error_el[j]->get_prev_state_vars() + k) = *(dual_el->get_prev_state_vars() + k);
							*(error_el[j]->get_state_vars() + k) = *(dual_el->get_state_vars() + k);
							*(error_el[j]->get_adjoint() + k) = *(dual_el->get_adjoint() + k);

						}
				}
				currentPtr = currentPtr->next;
			}

		}
}

void update_dual_grid(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->update_state(solrec, El_Table, iter);

				currentPtr = currentPtr->next;
			}
		}

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

void update_bilinear_error_grid(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;
	int yek = 1;
	double outflow = 0.;
	double tiny = GEOFLOW_TINY;

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				ErrorElem * Curr_El = (ErrorElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					gmfggetcoef_(Curr_El->get_bilin_prev_state(), Curr_El->get_d_state_vars(),
					    (Curr_El->get_d_state_vars() + NUM_STATE_VARS), Curr_El->get_dx(),
					    &(matprops_ptr->bedfrict[Curr_El->get_material()]), &(matprops_ptr->intfrict),
					    (Curr_El->get_kactxy()), (Curr_El->get_kactxy() + 1), &tiny,
					    &(matprops_ptr->epsilon));

				currentPtr = currentPtr->next;
			}
		}

	calc_edge_states<ErrorElem>(El_Table, NodeTable, matprops_ptr, timeprops_ptr, myid, &yek,
	    &outflow);

	slopes(El_Table, NodeTable, matprops_ptr, 2);

}

void clear_empty_jacobians(SolRec* solrec, int iter) {

	HashEntryPtr *buck = solrec->getbucketptr();
	for (int i = 0; i < solrec->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Jacobian * jacobian = (Jacobian*) (currentPtr->value);

				currentPtr = currentPtr->next;

				jacobian->erase_solution(iter - 1);

				if (jacobian->is_container_empty()) {
					solrec->remove(jacobian->get_key());
					delete jacobian;
				}

			}
		}

}
