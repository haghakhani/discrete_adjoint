/*
 * dual_refine_unrefine.C
 *
 *  Created on: Mar 3, 2016
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"
#include <algorithm>

bool is_old_son(Element* elem) {
	return (elem->get_adapted_flag() == OLDSON);
}
;

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
//	cout << " dual has been called \n";

	reset_adaption_flag(El_Table);

//	delete_unused_elements_nodes(El_Table, NodeTable, myid);

	double target = 0.1;

//	cout << "1 \n";
//	refinement_report(El_Table);

//	if (refinelist->get_num_elem()) {
	// start refinement
	for (int i = 0; i < refinelist->get_num_elem(); ++i) {
		refine(refinelist->get(i), El_Table, NodeTable, matprops_ptr);
		(refinelist->get(i))->put_adapted_flag(OLDFATHER);
		(refinelist->get(i))->put_refined_flag(1);
	}

//		cout << "2 \n";
//		refinement_report(El_Table);

	refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) refinelist, timeprops_ptr);

//		cout << "3 \n";
//		refinement_report(El_Table, myid);

	move_dual_data(meshctx, propctx);

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
				delete EmTemp;
				refdel++;
			}
		}
	}
//		cout << "5 \n";
//		refinement_report(El_Table);
//		cout << "number of deleted elem after ref " << refdel << "  number of ref list  "
//		    << refinelist->get_num_elem() << endl;
//	}

	MPI_Barrier(MPI_COMM_WORLD);

//	AssertMeshErrorFree(El_Table, NodeTable, numprocs, myid, -2.0);

//	if (unrefinelist->get_num_elem()) {

	for (int i = 0; i < unrefinelist->get_num_elem(); ++i)
		if (unrefinelist->get(i)->get_which_son() == 0) {
			first_son.push_back(unrefinelist->get(i));

			//we need the the new fathers for future
			//this is not required but is good for checking purposes
			new_father.push_back(
			    make_pair(*(unrefinelist->get(i)->getfather()),
			        *(unrefinelist->get(i)->getfather() + 1)));
		}

	//		cout << "6 \n";
	//		refinement_report(El_Table);
	int global_f_son_size = 0, size = first_son.size();
	MPI_Allreduce(&size, &global_f_son_size, 1, MPI_INT, MPI_SUM,
	MPI_COMM_WORLD);

	int counter = 0;

	while (global_f_son_size > 0) {

		for (int i = 0; i < first_son.size(); ++i)
//			first_son[i]->template find_brothers<T>(El_Table, NodeTable, target,
//					myid, matprops_ptr, &NewFatherList, &OtherProcUpdate, 1);
			first_son[i]->template dual_find_brothers<T>(El_Table, NodeTable, target, myid, matprops_ptr,
			    &NewFatherList, &OtherProcUpdate);

		unrefine_neigh_update(El_Table, NodeTable, myid, (void*) &NewFatherList);

		unrefine_interp_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &OtherProcUpdate);

		first_son.erase(remove_if(first_son.begin(), first_son.end(), is_old_son), first_son.end());

		for (int k = 0; k < NewFatherList.get_num_elem(); k++)
			delete_oldsons(El_Table, NodeTable, myid, NewFatherList.get(k));

		reset_newfather_adaption_flag(El_Table);

		move_dual_data(meshctx, propctx);

		NewFatherList.trashlist();
		OtherProcUpdate.trashlist();
		size = first_son.size();
		MPI_Allreduce(&size, &global_f_son_size, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
		counter++;
//		if (counter > 5)
//			cout << "myid " << myid << " size " << first_son.size() << " counter " << counter << endl;

	}
//	}

//	calc_d_gravity(El_Table);
//	cout << "8 \n";
//this is not required but good for checks
	set_new_fathers(El_Table, new_father);
//	refinement_report(El_Table, myid);
//	refine_flag_report(El_Table, myid);

	refinelist->trashlist();
	unrefinelist->trashlist();
//	reset_adaption_flag(El_Table);

}

template<>
void dual_refine_unrefine<ErrorElem>(MeshCTX* meshctx, PropCTX* propctx,
    ElemPtrList<ErrorElem>* refinelist, ElemPtrList<ErrorElem>* unrefinelist) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	ElemPtrList<Element> NewFatherList, OtherProcUpdate;

	vector<ErrorElem*> first_son;
	vector<pair<unsigned, unsigned> > new_father;

	reset_adaption_flag(El_Table);

//	delete_unused_elements_nodes(El_Table, NodeTable, myid);

	double target = 0.1;

//	cout << "1 \n";
//	refinement_report(El_Table);

//	if (refinelist->get_num_elem()) {
	// start refinement
	for (int i = 0; i < refinelist->get_num_elem(); ++i) {
		refine(refinelist->get(i), El_Table, NodeTable, matprops_ptr);
		(refinelist->get(i))->put_adapted_flag(OLDFATHER);
		(refinelist->get(i))->put_refined_flag(1);
	}

//		cout << "2 \n";
//	refinement_report(El_Table, myid);

	refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) refinelist, timeprops_ptr);

//		cout << "3 \n";
//		refinement_report(El_Table, myid);

	move_err_data(meshctx, propctx);

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
				delete EmTemp;
				refdel++;
			}
		}
	}
//		cout << "5 \n";
//		refinement_report(El_Table);
//		cout << "number of deleted elem after ref " << refdel << "  number of ref list  "
//		    << refinelist->get_num_elem() << endl;
//	}

	MPI_Barrier(MPI_COMM_WORLD);

//	AssertMeshErrorFree(El_Table, NodeTable, numprocs, myid, -2.0);

//	if (unrefinelist->get_num_elem()) {

	for (int i = 0; i < unrefinelist->get_num_elem(); ++i)
		if (unrefinelist->get(i)->get_which_son() == 0) {
			first_son.push_back(unrefinelist->get(i));

			//we need the the new fathers for future
			//this is not required but is good for checking purposes
			new_father.push_back(
			    make_pair(*(unrefinelist->get(i)->getfather()),
			        *(unrefinelist->get(i)->getfather() + 1)));
		}

	//		cout << "6 \n";
	//		refinement_report(El_Table);
	int global_f_son_size = 0, size = first_son.size();
	MPI_Allreduce(&size, &global_f_son_size, 1, MPI_INT, MPI_SUM,
	MPI_COMM_WORLD);

	int counter = 0;

	while (global_f_son_size > 0) {

		for (int i = 0; i < first_son.size(); ++i)
//			first_son[i]->template find_brothers<T>(El_Table, NodeTable, target,
//					myid, matprops_ptr, &NewFatherList, &OtherProcUpdate, 1);
			first_son[i]->dual_find_brothers<ErrorElem>(El_Table, NodeTable, target, myid, matprops_ptr,
			    &NewFatherList, &OtherProcUpdate);

		unrefine_neigh_update(El_Table, NodeTable, myid, (void*) &NewFatherList);

		unrefine_interp_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &OtherProcUpdate);

		first_son.erase(remove_if(first_son.begin(), first_son.end(), is_old_son), first_son.end());

		for (int k = 0; k < NewFatherList.get_num_elem(); k++)
			delete_oldsons(El_Table, NodeTable, myid, NewFatherList.get(k));

		reset_newfather_adaption_flag(El_Table);

		move_err_data(meshctx, propctx);

		NewFatherList.trashlist();
		OtherProcUpdate.trashlist();
		size = first_son.size();
		MPI_Allreduce(&size, &global_f_son_size, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
		counter++;
//		if (counter > 5)
//			cout << "myid " << myid << " size " << first_son.size() << " counter " << counter << endl;

	}
//	}

//	calc_d_gravity(El_Table);
//	cout << "8 \n";
//this is not required but good for checks
	set_new_fathers(El_Table, new_father);
//	refinement_report(El_Table, myid);
//	refine_flag_report(El_Table, myid);

	refinelist->trashlist();
	unrefinelist->trashlist();
//	reset_adaption_flag(El_Table);
//	refinement_report(El_Table, myid);

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

static void run_from_snapshot(MeshCTX* dual_meshctx, PropCTX* propctx, SolRec* solrec, vector<Snapshot>& snapshot_vec){

	propctx->runcond = (run_mode) (propctx->runcond & ~RECORD);
	const int numprocs = propctx->numproc;
	const int myid = propctx->myid;
	TimeProps* timeprops = propctx->timeprops;
	MatProps* matprops = propctx->matprops;
	const int vec_size = snapshot_vec.size() - 1;
	int index , j = 0;
	for (unsigned i = vec_size; i >= 0; --i){
		if (((timeprops->iter - snapshot_vec[i].get_iter()) >= solrec->get_range())|| snapshot_vec[i].get_iter()==0){
			index = i;
			break;
		}
		else
			snapshot_vec.erase(snapshot_vec.end() - j);
		j++;
	}

	// reconstructing hash tables
	HashTable *newel_table = new HashTable(snapshot_vec[index].get_elem_table_minimal());
	HashTable *newnd_table = new HashTable(snapshot_vec[index].get_node_table_minimal());

	// filling elem table and node table new information
	const vector<Node_minimal>* node_vector = snapshot_vec[index].get_node_vector();
	const int num_node = node_vector->size();
	for (int i = 0; i < num_node; ++i){
		Node* newnode = new Node(&((*node_vector)[i]));
		newnd_table->add(newnode->pass_key(), newnode);
	}

	const vector<Elem_minimal>* elem_vector = snapshot_vec[index].get_elem_vector();
	const int num_elem = elem_vector->size();
	for (int i = 0; i < num_elem; ++i){
		Element* newelem = new Element(&((*elem_vector)[i]), matprops);
		newel_table->add(newelem->pass_key(), newelem);
	}

//	setup_geoflow(newel_table, newnd_table, myid, numprocs, matprops, timeprops);
	setup_from_snapshot(newel_table, newnd_table, myid, numprocs, matprops, timeprops);

	move_data(numprocs, myid, newel_table, newnd_table, timeprops,matprops);

	MeshCTX newmeshctx;
	newmeshctx.nd_table = newnd_table;
	newmeshctx.el_table = newel_table;

	/* the time steps iteration is always one more than
	 * maximum time step at the end
	 */
	/*FIXME the max iter in forward solve is always +1 that it should*/
	const int final_iter = timeprops->iter-1;
	snapshot_vec[index].adjust_timeprops(timeprops,final_iter);

	solrec->free_all_available_sol();

	forward_solve(newmeshctx, *propctx, solrec);

	delete_hashtables_objects<Element>(newel_table);

	delete_hashtables_objects<Node>(newnd_table);

	snapshot_vec.erase(snapshot_vec.begin()+index);
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

	 vector<Snapshot> *snapshot_vec= dual_meshctx->snapshot_vec;

	read_solution.start();
	if (solrec->get_first_solution() && (solrec->get_first_solution() >= iter - 1)){
		run_from_snapshot(dual_meshctx, propctx, solrec, *snapshot_vec);
		//because the hash tabeles change in run_from_snapshot
//		El_Table = dual_meshctx->el_table;
//		NodeTable = dual_meshctx->nd_table;
	}
	read_solution.stop();

	//timing inside the function
	if (timeprops_ptr->ifrepartition() && propctx->adapt_flag && numprocs > 1) {
		dual_err_repartition(solrec, dual_meshctx, err_meshctx, propctx);
	}

	ElemPtrList<DualElem> refinelist, unrefinelist;
	ElemPtrList<ErrorElem> err_refinelist, err_unrefinelist;
	HashEntryPtr *buck;
	dual_adapt.start();
	if (timeprops_ptr->ifrefine() && propctx->adapt_flag
	    && !(timeprops_ptr->ifrepartition() && numprocs > 1)) {

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
#ifdef Error
		error_adapt.start();
		make_refine_unrefine_list_from_father(dual_meshctx, err_meshctx, &refinelist, &unrefinelist,
		    &err_refinelist, &err_unrefinelist);
		error_adapt.stop();
#endif

		dual_refine_unrefine<DualElem>(dual_meshctx, propctx, &refinelist, &unrefinelist);

		calc_d_gravity(El_Table);
	}

	update_dual_grid(solrec, dual_meshctx, propctx);
// inside updating state_vars we also delete the solution that we used
// since we no longer need it
	MPI_Barrier(MPI_COMM_WORLD);

#ifdef Error
	error_adapt.start();
	if (timeprops_ptr->ifrefine() && propctx->adapt_flag
	    && !(timeprops_ptr->ifrepartition() && numprocs > 1)) {

		dual_refine_unrefine<ErrorElem>(err_meshctx, propctx, &err_refinelist, &err_unrefinelist);

		correct_dual_err_link(err_meshctx, dual_meshctx);

		calc_d_gravity(cp_El_Table);
	}
	error_adapt.stop();
#endif

	clear_empty_jacobians(solrec, iter);
	dual_adapt.stop();
}
