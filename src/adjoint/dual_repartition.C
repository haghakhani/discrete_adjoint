/*
 * dual_repartition.C
 *
 *  Created on: Mar 3, 2016
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"
#include "../header/exvar.h"
#include <algorithm>

void dual_repartition(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = timeprops_ptr->iter;
	vector<TRANSKEY> trans_keys_vec;
	vector<int> trans_keys_status;

	//first of all each procs search it hashtable for missing elements
	ElemPtrList<DualElem> refinelist, unrefinelist, repart_list;

	HashEntryPtr * buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem *Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->dual_check_refine_unrefine_repartition(solrec, El_Table, iter, &refinelist,
					    &unrefinelist, trans_keys_vec, trans_keys_status, &repart_list);

				currentPtr = currentPtr->next;
			}
		}

	MPI_Barrier(MPI_COMM_WORLD);
	// now we have to look in other procs to find the missing elements

	// make a connection pairs between available processors
	//following lines initialize the connection loop
	int missing = trans_keys_vec.size(), receive_size;
	int receive_from = 0, send_to = 0, count = 0, found;

	vector<TRANSKEY> keys_to_check_vec;
	vector<int> my_keys_status, other_keys_status;

	MPI_Status status;
	MPI_Request trans_keys_req, trans_keys_send, key_status_rec, key_status_send, elem_status_rec,
	    elem_status_send;

	do {
		count++;
		set_send_receive_proc(count, myid, numprocs, receive_from, send_to);

		MPI_Send(&missing, 1, MPI_INT, send_to, count, MPI_COMM_WORLD);
		MPI_Recv(&receive_size, 1, MPI_INT, receive_from, count, MPI_COMM_WORLD, &status);

		if (missing > 0)
			MPI_Send(&trans_keys_vec[0], missing, TRANSKEYS, send_to, count, MPI_COMM_WORLD);

		if (receive_size > 0) {
			keys_to_check_vec.resize(receive_size);
			other_keys_status.resize(receive_size);

			MPI_Recv(&keys_to_check_vec[0], receive_size, TRANSKEYS, receive_from, count,
			MPI_COMM_WORLD, &status);
		}
		// first we open a window to show we are waiting to receive the data from the other proc
		MPI_Irecv(&my_keys_status, missing, MPI_INT, send_to, count, MPI_COMM_WORLD, &key_status_rec);

		found = 0;
		check_received_keys(solrec, keys_to_check_vec, other_keys_status, iter, found);

		MPI_Isend(&other_keys_status, receive_size, MPI_INT, receive_from, count, MPI_COMM_WORLD,
		    &key_status_send);

		DualElemPack* receive_array;

		if (found) {
			receive_array = new DualElemPack[found];
			MPI_Irecv(&receive_array, found, DUALELEMTYPE, receive_from, count, MPI_COMM_WORLD,
			    &elem_status_rec);
		}

		//now we have to ask the source proc to send us the elements that belong to us
		//first each proc asks its connection how many of elements have been found
		//note that we send and receive from different processors so it's better to
		//have non-blocking Communications
		int IfSentRecvd;
		do {
			MPI_Test(&key_status_rec, &IfSentRecvd, &status);
			if (IfSentRecvd) {
				int to_be_sent = 0;
				for (int i = 0; i < my_keys_status.size(); ++i)
					if (my_keys_status[i] > 0)
						to_be_sent++;

				if (to_be_sent) {

					DualElemPack* send_array = new DualElemPack[to_be_sent];

					int component = 0;
					for (int i = 0; i < my_keys_status.size(); ++i)
						if (my_keys_status[i] > 0) {
							repart_list.get(i)->Pack_element((send_array + component), NodeTable, send_to);
							component++;
							if (my_keys_status[i] == 1 || my_keys_status[i] == 2 || my_keys_status[i] == 12) {

								El_Table->remove(repart_list.get(i)->pass_key());
								delete repart_list.get(i);

							} else
								cerr << " it is not good, but my sibling is far from me\n";
						}

					MPI_Isend(send_array, to_be_sent, DUALELEMTYPE, send_to, count, MPI_COMM_WORLD,
					    &elem_status_send);
				}
			}
		} while (IfSentRecvd != 1);

		if (found)
			do {
				MPI_Test(&elem_status_rec, &IfSentRecvd, &status);

				if (IfSentRecvd)
					for (int i = 0; i < found; ++i) {
						DualElem* elm = (DualElem*) El_Table->lookup(receive_array[i].key);
						if (elm == NULL) { // this elm doesn't exist on this proc
							DualElem* new_elm = new DualElem((receive_array + i), NodeTable, myid);
							if ((new_elm->get_adapted_flag() < 0) && (new_elm->get_adapted_flag() >= -BUFFER))
								new_elm->put_myprocess(myid);
							El_Table->add(new_elm->pass_key(), new_elm);

						} else {

							cerr << "this should not happen hear \n";
							//this elm is already on this proc, rather than delete old copy
							//and allocate space for a new one, save time by only copying the
							//new element data to the old element.
//								elm->update((receive_array+i), NodeTable, myid);
//								if ((elm->get_adapted_flag() < 0) && (elm->get_adapted_flag() >= -BUFFER))
//									elm->put_myprocess(myid);
						}
					}
			} while (IfSentRecvd != 1);

		MPI_Barrier(MPI_COMM_WORLD);
		move_dual_data(meshctx, propctx);

		my_keys_status.clear();
		other_keys_status.clear();

	} while (count < numprocs - 1);

	dual_refine_unrefine<DualElem>(meshctx, propctx, &refinelist, &unrefinelist);

//		set_ithm(El_Table);
//
//		plot_ithm(El_Table);

	calc_d_gravity(El_Table);
}

void set_send_receive_proc(int count, int myid, int numprocs, int& receive_from, int& send_to) {

	int odd = count % 2;

	if (odd) {
		send_to = (myid + (int) ceil(count / 2.)) % numprocs;
		receive_from = (myid - (int) ceil(count / 2.) + numprocs) % numprocs;

	} else {
		receive_from = (myid + (int) ceil(count / 2.) + numprocs) % numprocs;
		send_to = (myid - (int) ceil(count / 2.) + numprocs) % numprocs;

	}

}

void SetTransPack(HashTable* El_Table, Element* elem, TRANSKEY* transelem) {

	for (int i = 0; i < KEYLENGTH; ++i) {
		transelem->key[i] = *(elem->pass_key() + i);
		transelem->key[i + KEYLENGTH] = *(elem->getfather() + i);
	}

	elem->gen_my_sons_key(El_Table, &(transelem->key[KEYLENGTH * KEYLENGTH]));
}

void check_received_keys(SolRec* solrec, vector<TRANSKEY>& keys_to_check_vec,
    vector<int>& keys_status, int iter, int &found) {

	for (int i = 0; i < keys_to_check_vec.size(); ++i) {
		Solution * prev_sol = solrec->lookup((keys_to_check_vec[i].key), iter - 1);

		//initialize to zero
		keys_status[i] = 0;

		if (prev_sol) {
			//means element founded
			keys_status[i] = 1;
			found++;
			continue;
		}

		//checks if the father exist
		prev_sol = solrec->lookup(&(keys_to_check_vec[i].key[KEYLENGTH]), iter - 1);

		if (prev_sol) {
			//means element founded and it has been unrefined
			found++;
			keys_status[i] = 2;
			continue;
		}

		int flag = 0;

		for (int j = 0; j < 4; ++j) {

			prev_sol = solrec->lookup(&(keys_to_check_vec[i].key[(2 + j) * KEYLENGTH]), iter - 1);
			if (prev_sol) {
				// for each son we add three points
				keys_status[i] += 3;
				if (!flag)
					flag = 1;
			}
		}

		if (flag)
			found++;
	}
}
