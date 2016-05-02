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
#include <algorithm>

double doubleKeyRange;

void inspect_element(HashTable* El_Table, unsigned* key);

void dual_repartition(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = timeprops_ptr->iter;

	doubleKeyRange = *(El_Table->get_doublekeyrange() + 1);

	unsigned dbg_key[2] = { 1, 2081 };

	//get the first and last key for this proc
	double myKeyRange[] = { DBL_MAX, -1. };

	vector<TRANSKEY> trans_keys_vec;
	vector<int> trans_keys_status;
	vector<DualElem*> repart_list;

	//first of all each procs search it hashtable for missing elements
	ElemPtrList<DualElem> refinelist, unrefinelist;

	HashEntryPtr * buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem *Curr_El = (DualElem*) (currentPtr->value);

				// because we may delete the it we first send it to next
				currentPtr = currentPtr->next;

				if (Curr_El->get_adapted_flag() > 0) {

					Curr_El->dual_check_refine_unrefine_repartition(solrec, El_Table, iter, &refinelist,
					    &unrefinelist, trans_keys_vec, trans_keys_status, repart_list, myKeyRange);

				} else {
					//during repartitioning we delete the GHOST elements as we do in forward run
					El_Table->remove(Curr_El->pass_key());
					delete Curr_El;
				}
			}
		}

	delete_extra_nodes(El_Table, NodeTable);

//	AssertMeshErrorFree(El_Table, NodeTable, numprocs, myid, 2.0);

//	inspect_element(El_Table, dbg_key);

	int count = 0, remaining;

	MPI_Status status;
	MPI_Request* s_request = new MPI_Request[2];
	MPI_Request* r_request = new MPI_Request[2];

	int IfSentRecvd;
	DualElemPack *receive_array, *send_array;

	do {
		count++;

		// make a connection pairs between available processors
		//following lines initialize the connection loop
		int send = trans_keys_vec.size(), receive_size;
		int receive_from = 0, send_to = 0, found = 0;
		vector<TRANSKEY> keys_to_check_vec;
		vector<int> my_keys_status(send), other_keys_status;

		set_send_receive_proc(count, myid, numprocs, receive_from, send_to);

		MPI_Send(&send, 1, MPI_INT, send_to, count, MPI_COMM_WORLD);
		MPI_Recv(&receive_size, 1, MPI_INT, receive_from, count, MPI_COMM_WORLD, &status);

		if (send > 0) {
			MPI_Send(&trans_keys_vec[0], send, TRANSKEYS, send_to, count, MPI_COMM_WORLD);
			// first we open a window to show we are waiting to receive the data from the other proc
			MPI_Irecv(&my_keys_status[0], send, MPI_INT, send_to, count, MPI_COMM_WORLD, &(r_request[0]));
		}

//		cout << "proc " << myid << " send " << send << " receive_size " << receive_size << " send to "
//		    << send_to << " receive_from " << receive_from << endl;

		if (receive_size > 0) {

			keys_to_check_vec.resize(receive_size);
			other_keys_status.resize(receive_size);

			MPI_Recv(&keys_to_check_vec[0], receive_size, TRANSKEYS, receive_from, count, MPI_COMM_WORLD,
			    &status);

			found = 0;
			check_received_keys(solrec, keys_to_check_vec, other_keys_status, iter, found);

//			update_my_key_range_receiving_elements(keys_to_check_vec, other_keys_status, myKeyRange);

			MPI_Isend(&other_keys_status[0], receive_size, MPI_INT, receive_from, count, MPI_COMM_WORLD,
			    &(s_request[0]));

			if (found) {
				receive_array = new DualElemPack[found];
				MPI_Irecv(receive_array, found, DUALELEMTYPE, receive_from, count, MPI_COMM_WORLD,
				    &(r_request[1]));
			}
		}

//		cout << "receive_size passed and " << found << " found \n";
		int to_be_sent = 0;

		if (send > 0) {
			//now we have to ask the source proc to send us the elements that belong to us
			//first each proc asks its connection how many of elements have been found
			//note that we send and receive from different processors so it's better to
			//have non-blocking Communications

			do {
				MPI_Test(&(r_request[0]), &IfSentRecvd, &status);
				if (IfSentRecvd) {
					to_be_sent = 0;
					for (int i = 0; i < my_keys_status.size(); ++i)
						if (my_keys_status[i] > 0)
							to_be_sent += 1;

					if (to_be_sent) {

						send_array = new DualElemPack[to_be_sent];

						int component = 0;
						for (int i = 0; i < my_keys_status.size(); ++i)
							if (my_keys_status[i] > 0) {
								repart_list[i]->Pack_element((send_array + component), NodeTable, send_to);
								component += 1;
//								cout << "status[" << i << "] " << my_keys_status[i] << endl;
								if (my_keys_status[i] == 1 || my_keys_status[i] == 2 || my_keys_status[i] == 12) {

									El_Table->remove(repart_list[i]->pass_key());
									delete repart_list[i];

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

//		cout << "send passed \n";

		if (found) {
//			cout << found << endl;
			do {
				MPI_Test(&(r_request[1]), &IfSentRecvd, &status);
				if (IfSentRecvd) {
					int component = 0;
					for (int i = 0; i < other_keys_status.size(); ++i) {
						if (other_keys_status[i] > 0) {

							DualElem* elm = (DualElem*) El_Table->lookup(receive_array[component].key);
							assert(elm == NULL);
//							cout << component << endl;

							elm = new DualElem((receive_array + component), NodeTable, myid);

							El_Table->add(elm->pass_key(), elm);

							double doublekey = *(elm->pass_key()) * doubleKeyRange + *(elm->pass_key() + 1);

							if (doublekey < myKeyRange[0])
								myKeyRange[0] = doublekey;

							if (doublekey > myKeyRange[1])
								myKeyRange[1] = doublekey;

							if (other_keys_status[i] == 2)

								unrefinelist.add(elm);

							else if (other_keys_status[i] == 3 || other_keys_status[i] == 6
							    || other_keys_status[i] == 9 || other_keys_status[i] == 12)

								refinelist.add(elm);
							else if (other_keys_status[i] != 1)
								cerr << "element status is not correct, and repartitioning fails \n";

							component++;
						}
					}
					delete[] receive_array;
				}

			} while (IfSentRecvd != 1);
		}

//		cout << " found passed " << endl;

		if (to_be_sent) {
			vector<TRANSKEY> cp_trans_keys_vec;
			vector<int> cp_trans_keys_status;
			vector<DualElem*> cp_repart_list;

			assert(trans_keys_vec.size() == my_keys_status.size());
			assert(trans_keys_status.size() == my_keys_status.size());
			assert(repart_list.size() == my_keys_status.size());
			// we must clean those elements that we sent from the list
			// and update status of those elements that we have found their sons in other procs
			for (int i = 0; i < my_keys_status.size(); ++i)
				// the following condition means those elements that we have not found them yet or just some of their sons have been found
				if (my_keys_status[i] % 3 == 0 && my_keys_status[i] != 12) {
					cp_repart_list.push_back(repart_list[i]);
					cp_trans_keys_status.push_back(my_keys_status[i] + trans_keys_status[i]);
					cp_trans_keys_vec.push_back(trans_keys_vec[i].key);
				}

			trans_keys_vec.clear();
			trans_keys_status.clear();
			repart_list.clear();

			trans_keys_vec = cp_trans_keys_vec;
			trans_keys_status = cp_trans_keys_status;
			repart_list = cp_repart_list;

			do {

				MPI_Test(&(s_request[1]), &IfSentRecvd, &status);
				if (IfSentRecvd)
					delete[] send_array;

			} while (IfSentRecvd != 1);

		}

//		cout << " second send passed \n";

		int nsize = trans_keys_vec.size();

		MPI_Allreduce(&nsize, &remaining, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

//		cout << "nsize " << nsize << " remaining  " << remaining << endl;

		// or when there is no missing element anymore
	} while (count < numprocs - 1 && remaining > 0);

//	cout << " out of while loop \n";

	delete_extra_nodes(El_Table, NodeTable);

//	cout << " extra node deleted \n";

	// we need this because we still have not refined and unrefined
	//	adjust_range(El_Table, refinelist, unrefinelist, myKeyRange);
	double *allKeyRange = new double[2 * numprocs];

	MPI_Allgather(myKeyRange, 2, MPI_DOUBLE, allKeyRange, 2, MPI_DOUBLE, MPI_COMM_WORLD);

//	cout << "after all gather my range is " << fixed << myKeyRange[0] << " , " << myKeyRange[1]
//	    << endl;

	allKeyRange[0] = -1.0;
	allKeyRange[2 * numprocs - 1] = DBL_MAX;

//	for (int i = 0; i < 2 * numprocs; ++i)
//		cout << allKeyRange[i] << " , ";

	update_neighbor_proc(propctx, El_Table, allKeyRange);

//	cout << " after update neighb proc \n";

	move_dual_data(meshctx, propctx);

//	cout << "after move data \n";

	dual_refine_unrefine<DualElem>(meshctx, propctx, &refinelist, &unrefinelist);

//	AssertMeshErrorFree(El_Table, NodeTable, numprocs, myid, 1.0);

//	cout << "after ref and unref \n";

	delete[] s_request;
	delete[] r_request;
	delete[] allKeyRange;

	calc_d_gravity(El_Table);
}

void set_send_receive_proc(int count, int myid, int numprocs, int& receive_from, int& send_to) {

	int odd = count % 2;

	if (odd) {
		send_to = (myid - (int) ceil(count / 2.) + numprocs) % numprocs;
		receive_from = (myid + (int) ceil(count / 2.)) % numprocs;

	} else {
		receive_from = (myid - (int) ceil(count / 2.) + numprocs) % numprocs;
		send_to = (myid + (int) ceil(count / 2.) + numprocs) % numprocs;

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

	//by the following we will adjust the overlap if it exists
	for (int i = 0; i < numprocs - 1; ++i)
		if (allKeyRange[2 * i + 1] > allKeyRange[2 * (i + 1)]) {
			double med = allKeyRange[2 * i + 1];
			allKeyRange[2 * i + 1] = allKeyRange[2 * (i + 1)];
			allKeyRange[2 * (i + 1)] = med;
		}

	vector<Element*> wrong_elem;

	HashEntryPtr * buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Element *Curr_El = (Element*) (currentPtr->value);

				assert(Curr_El);
				assert(Curr_El->get_adapted_flag()>=NOTRECADAPTED);

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

				double doublekey2 = *(Curr_El->pass_key()) * doubleKeyRange + *(Curr_El->pass_key() + 1);

				if (!(doublekey2 <= allKeyRange[myid * 2 + 1] && doublekey2 >= allKeyRange[myid * 2]))
					wrong_elem.push_back(Curr_El);

				currentPtr = currentPtr->next;
			}
		}

	vector<unsigned> lost_neighb;
	vector<int> lost_neighb_proc;

	for (int i = 0; i < wrong_elem.size(); ++i) {
		unsigned *key = wrong_elem[i]->pass_key();
		for (int neigh = 0; neigh < 8; neigh++) {
			if (*(wrong_elem[i]->get_neigh_proc() + neigh) >= 0) {
				unsigned *neigh_key = (wrong_elem[i]->get_neighbors() + neigh * KEYLENGTH);
				Element *neigh_elem = (Element*) El_Table->lookup(neigh_key);

				if (neigh_elem) {
					//luckily the neighbor also lives here
					int side = neigh_elem->which_neighbor(key);
					assert(side < 8 && side > -1);
					// we will update its information explicitly
					neigh_elem->put_neigh_proc(side, myid);
				} else {
					//must be updated by communication
					lost_neighb.push_back(neigh_key[0]);
					lost_neighb.push_back(neigh_key[1]);
					lost_neighb.push_back(key[0]);
					lost_neighb.push_back(key[1]);
					lost_neighb_proc.push_back(myid);
				}
			}
		}
	}

//	cout<<"size vec "<<lost_neighb_proc.size()<<endl;

	int sizes[numprocs];
	int size = lost_neighb_proc.size();
	MPI_Allgather(&size, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

	int sum = 0;
	int mystart[numprocs];
	for (int i = 0; i < numprocs; ++i) {
		mystart[i] = sum;
		sum += sizes[i];
	}

	int *check_lost_proc = new int[sum];

	MPI_Allgatherv(&lost_neighb_proc[0], sizes[myid], MPI_INT, check_lost_proc, sizes, mystart,
	MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < numprocs; ++i) {
		mystart[i] = 4 * mystart[i];
		sizes[i] = 4 * sizes[i];
	}

	unsigned *check_lost = new unsigned[4 * sum];
	MPI_Allgatherv(&lost_neighb[0], sizes[myid], MPI_UNSIGNED, check_lost, sizes, mystart,
	MPI_UNSIGNED, MPI_COMM_WORLD);
	unsigned *check = new unsigned[sum];
	unsigned *all_check = new unsigned[sum];

	for (int i = 0; i < sum; ++i)
		all_check[i] = check[i] = 0;

	int start = mystart[myid] / 4, end = (mystart[myid] + size) / 4;

	for (int i = 0; i < sum; ++i) {
		if (i < start || i >= end) {
			unsigned neigh_key[] = { check_lost[4 * i], check_lost[4 * i + 1] };
			Element *neigh_elem = (Element*) El_Table->lookup(neigh_key);
			if (neigh_elem) {
				// we found the element
				unsigned key[] = { check_lost[4 * i + 2], check_lost[4 * i + 3] };
				int side = neigh_elem->which_neighbor(key);
				assert(side < 8 && side > -1);
				neigh_elem->put_neigh_proc(side, check_lost_proc[i]);
				check[i] = 1;
				assert(check_lost_proc[i] >= 0 && check_lost_proc[i] < numprocs);
			}
		}
	}

	MPI_Allreduce(check, all_check, sum, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	for (int i = 0; i < sum; ++i)
		assert(all_check[i] == 1);

	delete[] check_lost;
	delete[] check_lost_proc;
	delete[] check;
	delete[] all_check;
	MPI_Barrier(MPI_COMM_WORLD);

}

void inspect_element(HashTable* El_Table, unsigned* key) {

	HashEntryPtr * buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Element *Curr_El = (Element*) (currentPtr->value);
				if (compare_key(Curr_El->pass_key(), key)) {
					cout << "this is something you have to avoid \n";
					exit(1);
				}

				currentPtr = currentPtr->next;
			}
		}
}
