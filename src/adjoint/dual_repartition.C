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
#include <float.h>

double doubleKeyRange;

void set_proc_state(HashTable* El_Table);

void adjust_range(HashTable* El_Table, ElemPtrList<DualElem>& refList,
    ElemPtrList<DualElem>& unRefList, double* myKeyRange);

void dual_repartition(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = timeprops_ptr->iter;

	doubleKeyRange = *(El_Table->get_doublekeyrange() + 1);

	//get the first and last key for this proc
	double myKeyRange[] = { DBL_MAX, -1. };

//	find_my_key_range(solrec, myKeyRange, iter);

//	cout << "myid " << myid << " my key range " << myKeyRange[0] << " , " << myKeyRange[1] << endl;

//	double *allKeyRange = new double[2 * numprocs];
//
//	MPI_Allgather(myKeyRange, 2, MPI_DOUBLE, allKeyRange, 2, MPI_DOUBLE, MPI_COMM_WORLD);

//	for (int i = 0; i < 2 * numprocs; ++i)
//		cout << allKeyRange[i] << " , ";

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

				if (Curr_El->get_adapted_flag() > 0) {

					Curr_El->dual_check_refine_unrefine_repartition(solrec, El_Table, iter, &refinelist,
					    &unrefinelist, trans_keys_vec, trans_keys_status, &repart_list, myKeyRange);

				} else {
					//during repartitioning we delete the GHOST elements as we do in forward run
					El_Table->remove(Curr_El->pass_key());
					delete Curr_El;
				}

				currentPtr = currentPtr->next;
			}
		}

	delete_extra_nodes(El_Table, NodeTable);

	// make a connection pairs between available processors
	//following lines initialize the connection loop
	int send = trans_keys_vec.size(), receive_size;
	int receive_from = 0, send_to = 0, count = 0, found = 0;

	vector<TRANSKEY> keys_to_check_vec;
	vector<int> my_keys_status(send), other_keys_status;

	MPI_Status status;
	MPI_Request* s_request = new MPI_Request[2];
	MPI_Request* r_request = new MPI_Request[2];

	int IfSentRecvd;
	DualElemPack *receive_array, *send_array;

	do {
		count++;
		set_send_receive_proc(count, myid, numprocs, receive_from, send_to);

		MPI_Send(&send, 1, MPI_INT, send_to, count, MPI_COMM_WORLD);
		MPI_Recv(&receive_size, 1, MPI_INT, receive_from, count, MPI_COMM_WORLD, &status);

		if (send > 0) {
			MPI_Send(&trans_keys_vec[0], send, TRANSKEYS, send_to, count, MPI_COMM_WORLD);
			// first we open a window to show we are waiting to receive the data from the other proc
			MPI_Irecv(&my_keys_status[0], send, MPI_INT, send_to, count, MPI_COMM_WORLD, &(r_request[0]));
		}

		if (receive_size > 0) {

			keys_to_check_vec.resize(receive_size);
			other_keys_status.resize(receive_size);

			MPI_Recv(&keys_to_check_vec[0], receive_size, TRANSKEYS, receive_from, count, MPI_COMM_WORLD,
			    &status);

			found = 0;
			check_received_keys(solrec, keys_to_check_vec, other_keys_status, iter, found);

			update_my_key_range_receiving_elements(keys_to_check_vec, other_keys_status, myKeyRange);

			MPI_Isend(&other_keys_status[0], receive_size, MPI_INT, receive_from, count, MPI_COMM_WORLD,
			    &(s_request[0]));

			if (found) {
				receive_array = new DualElemPack[found];
				MPI_Irecv(receive_array, found, DUALELEMTYPE, receive_from, count, MPI_COMM_WORLD,
				    &(r_request[1]));
			}
		}

		if (send > 0) {
			//now we have to ask the source proc to send us the elements that belong to us
			//first each proc asks its connection how many of elements have been found
			//note that we send and receive from different processors so it's better to
			//have non-blocking Communications

			do {
				MPI_Test(&(r_request[0]), &IfSentRecvd, &status);
				if (IfSentRecvd) {
					int to_be_sent = 0;
					for (int i = 0; i < my_keys_status.size(); ++i)
						if (my_keys_status[i] > 0)
							to_be_sent += 1;

					if (to_be_sent) {

						send_array = new DualElemPack[to_be_sent];

						int component = 0;
						for (int i = 0; i < my_keys_status.size(); ++i)
							if (my_keys_status[i] > 0) {
								repart_list.get(i)->Pack_element((send_array + component), NodeTable, send_to);
								component += 1;
								if (my_keys_status[i] == 1 || my_keys_status[i] == 2 || my_keys_status[i] == 12) {

									El_Table->remove(repart_list.get(i)->pass_key());
									delete repart_list.get(i);

								} else
									cerr << " it is not good, but my sibling is far from me since my status is:"
									    << my_keys_status[i] << " \n";
							}

						assert(component == to_be_sent);

						MPI_Isend(send_array, to_be_sent, DUALELEMTYPE, send_to, count, MPI_COMM_WORLD,
						    &(s_request[1]));

						delete_extra_nodes(El_Table, NodeTable);
					}
//					cout << to_be_sent << endl;
				}
			} while (IfSentRecvd != 1);
		}

		if (found) {
//			cout << found << endl;
			do {
				MPI_Test(&(r_request[1]), &IfSentRecvd, &status);
				if (IfSentRecvd) {
					found = 0;
					for (int i = 0; i < other_keys_status.size(); ++i) {
						if (other_keys_status[i] > 0) {

							DualElem* elm = (DualElem*) El_Table->lookup(receive_array[found].key);
							assert(elm == NULL);

							elm = new DualElem((receive_array + i), NodeTable, myid);

							El_Table->add(elm->pass_key(), elm);

							if (other_keys_status[i] == 1) {

								double doublekey = *(elm->pass_key()) * doubleKeyRange + *(elm->pass_key() + 1);

								if (doublekey < myKeyRange[0])
									myKeyRange[0] = doublekey;

								if (doublekey > myKeyRange[1])
									myKeyRange[1] = doublekey;

							} else if (other_keys_status[i] == 2) {
								unrefinelist.add(elm);

								double doublekey = *(elm->getfather()) * doubleKeyRange + *(elm->getfather() + 1);

								if (doublekey < myKeyRange[0])
									myKeyRange[0] = doublekey;

								if (doublekey > myKeyRange[1])
									myKeyRange[1] = doublekey;

							} else if (other_keys_status[i] == 3 || other_keys_status[i] == 6
							    || other_keys_status[i] == 9 || other_keys_status[i] == 12) {

								refinelist.add(elm);

								unsigned son_key[4][2];
								elm->gen_my_sons_key(El_Table, son_key);

								for (int i = 0; i < 4; ++i) {
									Solution* prev_sol = solrec->lookup(son_key[i], iter - 1);

									if (prev_sol) {
										double doublekey = son_key[i][0] * doubleKeyRange + son_key[i][1];

										if (doublekey < myKeyRange[0])
											myKeyRange[0] = doublekey;

										if (doublekey > myKeyRange[1])
											myKeyRange[1] = doublekey;
									}
								}

							} else
								cerr << "element status is not correct, and repartitioning fails \n";

							found++;
						}
					}
					delete[] receive_array;
				}

			} while (IfSentRecvd != 1);
		}

		if (send)
			do {

				MPI_Test(&(s_request[1]), &IfSentRecvd, &status);
				if (IfSentRecvd)
					delete[] send_array;

			} while (IfSentRecvd != 1);

		MPI_Barrier(MPI_COMM_WORLD);

		my_keys_status.clear();
		other_keys_status.clear();

		// or when there is no missing element anymore
	} while (count < numprocs - 1);

	delete_extra_nodes(El_Table, NodeTable);

	// we need this because we still have not refined and unrefined
//	adjust_range(El_Table, refinelist, unrefinelist, myKeyRange);
	double *allKeyRange = new double[2 * numprocs];

	MPI_Allgather(myKeyRange, 2, MPI_DOUBLE, allKeyRange, 2, MPI_DOUBLE, MPI_COMM_WORLD);

	allKeyRange[0] = -1;

//	for (int i = 0; i < 2 * numprocs; ++i)
//		cout << allKeyRange[i] << " , ";

	update_neighbor_proc(propctx, El_Table, allKeyRange);

	move_dual_data(meshctx, propctx);

	dual_refine_unrefine<DualElem>(meshctx, propctx, &refinelist, &unrefinelist);


	delete[] s_request;
	delete[] r_request;
	delete[] allKeyRange;

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

void SetTransPack(TRANSKEY& transelem, unsigned* key, unsigned* father_key, unsigned son_key[][2]) {

	for (int i = 0; i < KEYLENGTH; ++i) {

		transelem.key[i] = key[i];
		transelem.key[i + KEYLENGTH] = father_key[i];

		for (int j = 0; j < 4; ++j)
			transelem.key[i + KEYLENGTH * (KEYLENGTH + j)] = son_key[j][i];
	}

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

void delete_extra_nodes(HashTable* El_Table, HashTable* NodeTable) {

	int NodeTable_num_buck = NodeTable->get_no_of_buckets();
	HashEntryPtr *NodeTable_bucket0 = NodeTable->getbucketptr();
	HashEntryPtr NodeTable_entry_ptr;
	int inodebucket, inode;
	Node *NdTemp;

//zero the number of elems each node is associated with
	for (inodebucket = 0; inodebucket < NodeTable_num_buck; inodebucket++) {
		NodeTable_entry_ptr = *(NodeTable_bucket0 + inodebucket);

		while (NodeTable_entry_ptr) {

			NdTemp = (Node*) (NodeTable_entry_ptr->value);
			NodeTable_entry_ptr = NodeTable_entry_ptr->next;
			assert(NdTemp);
			NdTemp->put_num_assoc_elem(0);
		}
	}

	int num_buck = El_Table->get_no_of_buckets();
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element *EmTemp;

	for (int ibuck = 0; ibuck < num_buck; ibuck++) {
		currentPtr = *(buck + ibuck);

		while (currentPtr) {
			EmTemp = (Element*) (currentPtr->value);
			currentPtr = currentPtr->next;
			assert(EmTemp);
			NdTemp = (Node*) NodeTable->lookup(EmTemp->pass_key());
			assert(NdTemp);
			NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem() + 1);

			for (inode = 0; inode < 8; inode++) {
				NdTemp = (Node*) NodeTable->lookup(EmTemp->getNode() + inode * KEYLENGTH);
				assert(NdTemp);
				NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem() + 1);
			}
		}
	}

	//if a node is not associated with any elements delete it
	for (inodebucket = 0; inodebucket < NodeTable_num_buck; inodebucket++) {
		NodeTable_entry_ptr = *(NodeTable_bucket0 + inodebucket);

		while (NodeTable_entry_ptr) {

			NdTemp = (Node*) (NodeTable_entry_ptr->value);
			NodeTable_entry_ptr = NodeTable_entry_ptr->next;
			assert(NdTemp);

			if (NdTemp->get_num_assoc_elem() == 0) {
				NodeTable->remove(NdTemp->pass_key());
				delete NdTemp;
			}
		}
	}
}

void update_neighbor_proc(PropCTX* propctx, HashTable* El_Table, double * allKeyRange) {

	int myid = propctx->myid, numprocs = propctx->numproc;
//	int num_elem=0;

	HashEntryPtr * buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem *Curr_El = (DualElem*) (currentPtr->value);

				assert(Curr_El);
				assert(Curr_El->get_adapted_flag()>=NOTRECADAPTED);
//				num_elem++;

				Curr_El->put_myprocess(myid);

				for (int ineigh = 0; ineigh < 8; ineigh++)
					if (*(Curr_El->get_neigh_proc() + ineigh) >= 0) {

						//make a double precision version of the neighbor's key
						double doublekey = *(Curr_El->get_neighbors() + ineigh * KEYLENGTH) * doubleKeyRange
						    + *(Curr_El->get_neighbors() + ineigh * KEYLENGTH + 1);

						//check which processor the neighbor's key belongs to
						for (int iproc = 0; iproc < numprocs; iproc++)
							if ((allKeyRange[2 * iproc] <= doublekey)
							    && (allKeyRange[2 * iproc + 1] >= doublekey)) {
								Curr_El->put_neigh_proc(ineigh, iproc);
								break;
							}
					}

				currentPtr = currentPtr->next;
			}
		}
}

void update_my_key_range_leaving_elements(vector<TRANSKEY>& trans_keys_vec,
    vector<int>& trans_keys_status, vector<double>& myDoubleKeys) {

	assert(trans_keys_status.size() == trans_keys_vec.size());
	vector<double> repDoubleKey;

	for (int i = 0; i < trans_keys_status.size(); ++i)
//		if (trans_keys_status[i] == 1 || trans_keys_status[i] == 2 || trans_keys_status[i] == 12) {
		repDoubleKey.push_back(trans_keys_vec[i].key[0] * doubleKeyRange + trans_keys_vec[i].key[1]);

}

void update_my_key_range_receiving_elements(vector<TRANSKEY>& keys_to_check_vec,
    vector<int>& other_keys_status, double* keyRange) {

	assert(other_keys_status.size() == keys_to_check_vec.size());

	for (int i = 0; i < other_keys_status.size(); ++i)

		if (other_keys_status[i] > 0) {

			double doubleKey = keys_to_check_vec[i].key[0] * doubleKeyRange + keys_to_check_vec[i].key[1];

			if (doubleKey < keyRange[0])
				keyRange[0] = doubleKey;

			if (doubleKey > keyRange[1])
				keyRange[1] = doubleKey;

		}
}

void find_my_key_range(SolRec* solrec, double* myKeyRange, int iter) {

	HashEntryPtr * buck = solrec->getbucketptr();
	for (int i = 0; i < solrec->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Jacobian *jacobian = (Jacobian*) (currentPtr->value);

				Solution* solution = jacobian->get_solution(iter - 1);
				if (solution) {
					unsigned* key = jacobian->get_key();
					double doubleKey = key[0] * doubleKeyRange + key[1];
					if (doubleKey > myKeyRange[1])
						myKeyRange[1] = doubleKey;
					if (doubleKey < myKeyRange[0])
						myKeyRange[0] = doubleKey;
				}

				currentPtr = currentPtr->next;
			}
		}
}

void adjust_range(HashTable* El_Table, ElemPtrList<DualElem>& refList,
    ElemPtrList<DualElem>& unRefList, double* myKeyRange) {

	for (int i = 0; i < refList.get_num_elem(); ++i) {

		unsigned* key = refList.get(i)->pass_key();
		double doubleKey = key[0] * doubleKeyRange + key[1];

		if (doubleKey < myKeyRange[0]) {
			myKeyRange[0] = doubleKey;
			cout << "key min updated\n";
		}

		if (doubleKey > myKeyRange[1]) {
			myKeyRange[1] = doubleKey;
			cout << "key max updated\n";
		}
	}

	for (int i = 0; i < unRefList.get_num_elem(); ++i) {

		unsigned* key = unRefList.get(i)->pass_key();
		double doubleKey = key[0] * doubleKeyRange + key[1];

		if (doubleKey < myKeyRange[0]) {
			myKeyRange[0] = doubleKey;
			cout << "key min updated\n";
		}

		if (doubleKey > myKeyRange[1]) {
			myKeyRange[1] = doubleKey;
			cout << "key max updated\n";
		}
	}
//	if (status == 2) {
//		unsigned* father_key = elm->getfather();
//		double doubleKey = father_key[0] * doubleKeyRange + father_key[1];
//
//		if (doubleKey < myKeyRange[0]) {
//			myKeyRange[0] = doubleKey;
//			cout << "key min updated\n";
//		}
//
//		if (doubleKey > myKeyRange[1]) {
//			myKeyRange[1] = doubleKey;
//			cout << "key max updated\n";
//		}
//	}
//
//	if (status == 3 || status == 4 || status == 8 || status == 12) {
//
//		unsigned son_key[4][2];
//		elm->gen_my_sons_key(El_Table, son_key);
//		double doubleKey[4];
//
//		for (int i = 0; i < 4; ++i) {
//			doubleKey[i] = son_key[i][0] * doubleKeyRange + son_key[i][1];
//
//			if (doubleKey[i] < myKeyRange[0]) {
//				myKeyRange[0] = doubleKey[i];
//				cout << "key min updated\n";
//			}
//
//			if (doubleKey[i] > myKeyRange[1]) {
//				myKeyRange[1] = doubleKey[i];
//				cout << "key max updated\n";
//			}
//		}
//	}
}

void set_proc_state(HashTable* El_Table) {

	HashEntryPtr * buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem *Curr_El = (DualElem*) (currentPtr->value);

				*(Curr_El->get_prev_state_vars()) = 0.;
				*(Curr_El->get_prev_state_vars() + 1) = 0.;

				if (Curr_El->get_adapted_flag() > 0) {
					*(Curr_El->get_prev_state_vars()) = Curr_El->get_myprocess();

					int *neigh_proc = Curr_El->get_neigh_proc();

					for (int i = 0; i < 8; ++i)
						if (neigh_proc[i] >= 0 && neigh_proc[i] != Curr_El->get_myprocess()) {
							*(Curr_El->get_prev_state_vars() + 1) = GHOST;
							break;
						}
				}
				currentPtr = currentPtr->next;
			}
		}
}

