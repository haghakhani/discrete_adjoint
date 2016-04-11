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

	reset_adaption_flag(El_Table);

//	delete_unused_elements_nodes(El_Table, NodeTable, myid);

	double target = 0.1;

//	cout << "1 \n";
//	refinement_report(El_Table);

//	if (refinelist->get_num_elem()) {
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
		MPI_Allreduce(&size, &global_f_son_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		while (global_f_son_size > 0) {

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
			OtherProcUpdate.trashlist();
			size = first_son.size();
			MPI_Allreduce(&size, &global_f_son_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		}

		move_dual_data(meshctx, propctx);

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
		solrec->load_new_set_of_solution(myid);
	}

	HashEntryPtr *buck;
#ifdef Error
	if (iter % 5 == 0 && propctx->adapt_flag) {

		ElemPtrList<ErrorElem> refinelist, unrefinelist;

		buck = cp_El_Table->getbucketptr();

		for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				ErrorElem *Curr_El = (ErrorElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
				Curr_El->error_check_refine_unrefine(solrec, cp_El_Table, iter,
						&refinelist, &unrefinelist);

				currentPtr = currentPtr->next;
			}
		}
//		cout<<"in refined table"<<endl;
//		cout<<"has to be refined "<<refinelist->get_num_elem()<<endl;
//		cout<<"has to be unrefined "<<unrefinelist->get_num_elem()<<endl;

		dual_refine_unrefine<ErrorElem>(err_meshctx, propctx, &refinelist,
				&unrefinelist);

		calc_d_gravity(cp_El_Table);

	}

	if (iter != 1)
	update_error_grid(solrec, err_meshctx, propctx);
#endif

	if (timeprops_ptr->ifrepartition() && propctx->adapt_flag && numprocs>1) {
		dual_repartition(solrec, dual_meshctx, propctx);
		//		cout<<"in original table"<<endl;
//		cout<<"has to be refined "<<refinelist->get_num_elem()<<endl;
//		cout<<"has to be unrefined "<<unrefinelist->get_num_elem()<<endl;

	}

	if (timeprops_ptr->ifrefine() && propctx->adapt_flag) {

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

		dual_refine_unrefine<DualElem>(dual_meshctx, propctx, &refinelist, &unrefinelist);

		calc_d_gravity(El_Table);

	}

// inside updating state_vars we also delete the solution that we used
// since we no longer need it
	MPI_Barrier(MPI_COMM_WORLD);
	update_dual_grid(solrec, dual_meshctx, propctx);

	clear_empty_jacobians(solrec, iter);

}
