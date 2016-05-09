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
//Jan 30, 2015
//haghakha
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define KEY0   3777862041
#define KEY1   2576980374
//#define DEBUG

void uinform_refine(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	int myid = propctx->myid, numprocs = propctx->numproc;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	ErrorElem* Curr_El = NULL;

	//for debugging perpose
	unsigned key[2] = { KEY0, KEY1 };
	double max = 0;

#ifdef DEBUG
	double dummyv_star = 0.0;
	int adjflag = 1;
	tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0.,adjflag);

	int nonz1 = num_nonzero_elem(El_Table);

	cout << "number of elements before refinement  " << nonz1 << endl;

	int *dbgvec = new int[nonz1];
	int *pass = new int[nonz1];
	for (int i = 0; i < nonz1; i++) {
		dbgvec[i] = 0;
		pass[i] = 0;
	}

#endif

//	htflush(El_Table, NodeTable, 1);
	move_err_data(meshctx, propctx);
	ElemPtrList<ErrorElem> RefinedList(num_nonzero_elem(El_Table));

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (ErrorElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() >= NOTRECADAPTED) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
					RefinedList.add(Curr_El);
				}
				currentPtr = currentPtr->next;
			}
		}
	}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (ErrorElem*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() == NOTRECADAPTED) {
					refine(Curr_El, El_Table, NodeTable, matprops_ptr);
				}

			}
		}
	}

//	cout << "here is the problem   "<<checkElement(El_Table, &max, key) << endl;

#ifdef DEBUG

	cout << "number of elements -7   " << num_nonzero_elem(El_Table, -7) << endl
	<< "number of elements -6   " << num_nonzero_elem(El_Table, -6) << endl
	<< "number of elements  0   " << num_nonzero_elem(El_Table, 0) << endl
	<< "number of elements  1   " << num_nonzero_elem(El_Table, 1) << endl
	<< "number of elements  2   " << num_nonzero_elem(El_Table, 2) << endl
	<< "number of elements  3   " << num_nonzero_elem(El_Table, 3) << endl
	<< "number of elements  4   " << num_nonzero_elem(El_Table, 4) << endl
	<< "number of elements  5   " << num_nonzero_elem(El_Table, 5) << endl;
#endif

	//this function delete old father elements
	uniform_refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &RefinedList,
	    timeprops_ptr);

//	cout << "here is the problem   "<<checkElement(El_Table, &max, key) << endl;

	move_err_data(meshctx, propctx);

	adjust_node_info(meshctx, propctx);

//	AssertMeshErrorFree(El_Table, NodeTable, numprocs, myid, 0);

	reset_adaption_flag(El_Table);
	delete_extra_nodes(El_Table, NodeTable);

//	double outflow=0.;
//	int order_flag=1;
//
//	calc_edge_states(El_Table, NodeTable, matprops_ptr, timeprops_ptr, myid, &order_flag, &outflow);
//
//	slopes(El_Table, NodeTable, matprops_ptr, 1);

}

void uniform_refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int nump, int myid,
    void* RL, TimeProps* timeprops_ptr) {

//	Element* EmFather;
//	Element* EmSon[4];
//	Element* EmNeighNew[4];
//	Element* EmNeighOld[2];
	ElemPtrList<Element>* RefinedList = (ElemPtrList<Element>*) RL;
//	Node* NdTemp;
//	int ifather, iside, ineigh, ineighp4, isonA, isonB;
//	int ineighme, ineighmep4, ineighson, ikey, inewcase, inode;
//
//	int iproc, ineighm4;
//
//	/*************************************************************/
//	/* now do the on processor updates while I'm waiting to      */
//	/* receive neighbor update information from other processors */
//	/*************************************************************/
//
//	for (ifather = RefinedList->get_inewstart(); ifather < RefinedList->get_num_elem(); ifather++) {
//
//		EmFather = RefinedList->get(ifather); //Hello I'm the OLDFATHER
//		assert(EmFather); //Help I've been abducted call the FBI!!!
//		assert(EmFather->adapted==OLDFATHER); //sanity check
//
//		NdTemp = (Node*) NodeTable->lookup(EmFather->key);
//		assert(NdTemp);
//		NdTemp->info = CORNER;
//
//		//These are my sons, I'm going to introduce them to my neighbors
//		for (isonA = 0; isonA < 4; isonA++) {
//			EmSon[isonA] = (Element*) El_Table->lookup(EmFather->son[isonA]);
//			assert(EmSon[isonA]); //MY son has been abducted call the FBI!!!
//
//			NdTemp = (Node*) NodeTable->lookup(EmSon[isonA]->key);
//			assert(NdTemp);
//			NdTemp->info = BUBBLE;
//
//			NdTemp = (Node*) NodeTable->lookup(EmSon[isonA]->node_key[(isonA + 1) % 4 + 4]);
//			assert(NdTemp);
//			NdTemp->info = SIDE;
//		}
//
//		//visit my neighbors on each side
//		for (iside = 0; iside < 4; iside++) {
//
//			ineigh = iside;
//			ineighp4 = ineigh + 4;
//			isonA = ineigh;
//			isonB = (ineighp4 + 1) % 4;
//
//			if (EmFather->neigh_proc[ineigh] == -1) {
//				//handle map boundary special
//				for (ikey = 0; ikey < KEYLENGTH; ikey++) {
//					EmSon[isonA]->neighbor[ineigh][ikey] = EmSon[isonA]->neighbor[ineighp4][ikey] =
//					    EmSon[isonB]->neighbor[ineigh][ikey] = EmSon[isonB]->neighbor[ineighp4][ikey] = 0;
//				}
//				EmSon[isonA]->neigh_gen[ineigh] = EmSon[isonA]->neigh_gen[ineighp4] =
//				    EmSon[isonB]->neigh_gen[ineigh] = EmSon[isonB]->neigh_gen[ineighp4] = 0;
//
//				EmSon[isonA]->neigh_proc[ineigh] = EmSon[isonB]->neigh_proc[ineigh] = -1;
//
//				EmSon[isonA]->neigh_proc[ineighp4] = EmSon[isonB]->neigh_proc[ineighp4] = -2;
//				//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//			} else if ((EmFather->neigh_proc[ineigh] == myid)
//			    && ((EmFather->neigh_proc[ineighp4] == myid) || (EmFather->neigh_proc[ineighp4] == -2))) {
//				//case where one neighbor on this side is on my proc while the other
//				//is on another proc has already been handled up above, when packing
//				//the information to send to the other proc.
//
//				//knock knock, Hello Neighbors
//				EmNeighOld[0] = (Element*) El_Table->lookup(EmFather->neighbor[ineigh]);
//				assert(EmNeighOld[0]);
//
//				EmNeighOld[1] = (Element*) El_Table->lookup(EmFather->neighbor[ineighp4]);
//				assert(EmNeighOld[1]);
//				EmNeighNew[0] = EmNeighNew[1] = EmNeighNew[2] = EmNeighNew[3] = NULL;
//
//				for (ineighme = 0; ineighme < 8; ineighme++) {
//					if (compare_key(EmFather->key, EmNeighOld[1]->neighbor[ineighme]))
//						break;
//
//				}
//				if (!(ineighme < 8)) {
//					printf(
//					    "FUBAR 0 detected in refine_neigh_update\nEmFather={%10u,%10u}\nEmNeighOld[0]={%10u,%10u} ineigh=%d isonA=%d isonB=%d\n",
//					    EmFather->key[0], EmFather->key[1], EmNeighOld[0]->key[0], EmNeighOld[0]->key[1],
//					    ineigh, isonA, isonB);
//					printf("aborting!!");
//
//					assert(ineighme < 8);
//				}
//
//				//There are 5 cases I need to worry about about
//				//A: my old neighbor is one generation older than me, only
//				//   possible if we've both been refined
//				//B: my old neighbor is of my generation and hasn't been refined
//				//C: my old neighbor is of my generation and has been refined
//				//D: my old neighbor is one generation younger than me and hasn't
//				//   been refined
//				//E: my old neighbor is one generation younger than me and has
//				//   been refined
//				//
//				//I'm going to compress this into one of 3 cases I need to handle
//				//0: case A: same as case E but from the other side, only need
//				//   to do once so don't do anything this time
//				//1: my new neighbor is my generation
//				//2: my new neighbor is one generation younger than me
//				//3: my new neighbor is two generations younger than me
//
//				//this switch loop compresses cases
//				switch (EmNeighOld[0]->generation - EmFather->generation) {
//					case -1:
//						//this is a case A
//						inewcase = 0;
//						assert(EmNeighOld[0]->adapted==OLDFATHER);
//						ineighme = ineighmep4 = -1; //for sanity check
//						break;
//					case 0:
//						assert(ineighme < 4);
//						ineighmep4 = ineighme + 4;
//
//						if (EmNeighOld[0]->adapted == OLDFATHER) {
//							//this is a case C
//							inewcase = 2;
//
//							EmNeighNew[0] = (Element*) El_Table->lookup(EmNeighOld[0]->son[(ineighme + 1) % 4]);
//							assert(EmNeighNew[0]);
//
//							EmNeighNew[1] = (Element*) El_Table->lookup(EmNeighOld[0]->son[ineighme]);
//							assert(EmNeighNew[1]);
//							EmNeighNew[2] = EmNeighNew[3] = NULL;
//						} else {
//							//this is a case B
//							if (EmNeighOld[0]->adapted > TOBEDELETED) {
//								inewcase = 1;
//								EmNeighNew[0] = EmNeighOld[0];
//								EmNeighNew[1] = EmNeighNew[2] = EmNeighNew[3] = NULL;
//							} else
//								inewcase = 0;
//
//						}
//						break;
//					case 1:
//						assert(ineighme < 4);
//						ineighmep4 = ineighme + 4;
//
//						if ((EmNeighOld[0]->adapted == OLDFATHER) || (EmNeighOld[1]->adapted == OLDFATHER)) {
//							//this is a case E
//							inewcase = 3;
//
//							if (EmNeighOld[0]->adapted == OLDFATHER) {
//								EmNeighNew[0] = (Element*) El_Table->lookup(EmNeighOld[0]->son[(ineighme + 1) % 4]);
//								assert(EmNeighNew[0]);
//
//								EmNeighNew[1] = (Element*) El_Table->lookup(EmNeighOld[0]->son[ineighme]);
//								assert(EmNeighNew[1]);
//							} else
//								EmNeighNew[1] = EmNeighNew[0] = EmNeighOld[0];
//
//							if (EmNeighOld[1]->adapted == OLDFATHER) {
//								EmNeighNew[2] = (Element*) El_Table->lookup(EmNeighOld[1]->son[(ineighme + 1) % 4]);
//								assert(EmNeighNew[2]);
//
//								EmNeighNew[3] = (Element*) El_Table->lookup(EmNeighOld[1]->son[ineighme]);
//								assert(EmNeighNew[3]);
//							} else
//								EmNeighNew[3] = EmNeighNew[2] = EmNeighOld[1];
//						} else {
//							//this is a case D
//							inewcase = 2;
//
//							EmNeighNew[0] = EmNeighOld[0];
//							EmNeighNew[1] = EmNeighOld[1];
//							EmNeighNew[2] = EmNeighNew[3] = NULL;
//						}
//						break;
//					default:
//						inewcase = -1;
//
//						printf("FUBAR 1 detected in refine_neigh_update! aborting.\n");
//						assert(0);
//						break;
//				} //switch based on difference in generation between me and my old neighbor, this is used to reduce the number of cases from 5 to 3 (based on new neighbor generation)
//
//				//sanity check
//				assert((ineigh >= 0) && (ineigh < 4));
//				assert(ineighp4 == ineigh + 4);
//				if (inewcase) {
//					assert((ineighme >= 0) && (ineighme < 4));
//					assert(ineighmep4 == ineighme + 4);
//				}
//
//				//now only deal with the new cases, and yes I know that I
//				//am resetting neighbor information in ghost cells but
//				//not neighbor information of the original cells on other
//				//processors, I'm going to fix that in a minute
//				switch (inewcase) {
//					case 0:
//						//case A
//						break;
//					case 1:
//						//case B
//						//new neighbor generation is my (the OLDFATHER) generation
//						for (ikey = 0; ikey < KEYLENGTH; ikey++) {
//							EmNeighNew[0]->neighbor[ineighme][ikey] = EmSon[isonB]->key[ikey];
//							EmNeighNew[0]->neighbor[ineighmep4][ikey] = EmSon[isonA]->key[ikey];
//
//							EmSon[isonA]->neighbor[ineigh][ikey] = EmSon[isonA]->neighbor[ineighp4][ikey] =
//							    EmSon[isonB]->neighbor[ineigh][ikey] = EmSon[isonB]->neighbor[ineighp4][ikey] =
//							        EmNeighNew[0]->key[ikey];
//						}
//
//						EmNeighNew[0]->neigh_gen[ineighme] = EmNeighNew[0]->neigh_gen[ineighmep4] =
//						    EmSon[isonA]->generation;
//
//						EmSon[isonA]->neigh_gen[ineigh] = EmSon[isonA]->neigh_gen[ineighp4] =
//						    EmSon[isonB]->neigh_gen[ineigh] = EmSon[isonB]->neigh_gen[ineighp4] =
//						        EmNeighNew[0]->generation;
//
//						EmSon[isonA]->neigh_proc[ineighp4] = EmSon[isonB]->neigh_proc[ineighp4] = -2;
//
//						EmNeighNew[0]->neigh_proc[ineighme] = EmNeighNew[0]->neigh_proc[ineighmep4] =
//						    EmFather->myprocess;
//
//						//update the nodes on this side
//						//The new difference in generation tells me the OLDFATHER's
//						//2 corner nodes on this side are actually CORNER's and not
//						//S_C_CON's
//
//						inode = ineigh;
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(CORNER);
//
//						inode = ineighp4;
//						NdTemp = (Node*) NodeTable->lookup(EmSon[isonA]->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(S_S_CON);
//
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(S_C_CON);
//
//						NdTemp = (Node*) NodeTable->lookup(EmSon[isonB]->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(S_S_CON);
//
//						inode = (ineigh + 1) % 4;
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(CORNER);
//
//						//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//
//						break;
//					case 2:
//						//cases C & D
//						//new neighbor generation is my son's generation
//
//						for (ikey = 0; ikey < KEYLENGTH; ikey++) {
//
//							EmNeighNew[0]->neighbor[ineighme][ikey] = EmNeighNew[0]->neighbor[ineighmep4][ikey] =
//							    EmSon[isonA]->key[ikey];
//							EmNeighNew[1]->neighbor[ineighme][ikey] = EmNeighNew[1]->neighbor[ineighmep4][ikey] =
//							    EmSon[isonB]->key[ikey];
//
//							EmSon[isonA]->neighbor[ineigh][ikey] = EmSon[isonA]->neighbor[ineighp4][ikey] =
//							    EmNeighNew[0]->key[ikey];
//							EmSon[isonB]->neighbor[ineigh][ikey] = EmSon[isonB]->neighbor[ineighp4][ikey] =
//							    EmNeighNew[1]->key[ikey];
//						}
//
//						EmNeighNew[0]->neigh_gen[ineighme] = EmNeighNew[0]->neigh_gen[ineighmep4] =
//						    EmNeighNew[1]->neigh_gen[ineighme] = EmNeighNew[1]->neigh_gen[ineighmep4] =
//						        EmSon[isonA]->generation;
//
//						EmNeighNew[0]->neigh_proc[ineighmep4] = EmNeighNew[1]->neigh_proc[ineighmep4] =
//						    EmSon[isonA]->neigh_proc[ineighp4] = EmSon[isonB]->neigh_proc[ineighp4] = -2;
//
//						EmSon[isonA]->neigh_gen[ineigh] = EmSon[isonA]->neigh_gen[ineighp4] =
//						    EmSon[isonB]->neigh_gen[ineigh] = EmSon[isonB]->neigh_gen[ineighp4] =
//						        EmNeighNew[0]->generation;
//
//						EmSon[isonA]->neigh_proc[ineigh] = EmNeighNew[0]->myprocess;
//
//						EmSon[isonB]->neigh_proc[ineigh] = EmNeighNew[1]->myprocess;
//
//						EmNeighNew[0]->neigh_proc[ineighme] = EmNeighNew[1]->neigh_proc[ineighme] =
//						    EmFather->myprocess;
//
//						//update the nodes on this side
//						//don't update my corner nodes because they could be S_C_CON's
//						//if they should be S_C_CON's and I reset them to CORNERs I
//						//will no longer conserve mass/volume in a dramatically
//						//observable fashion
//
//						if (EmSon[isonA]->neigh_gen[(ineigh + 3) % 4] == EmSon[isonA]->generation) {
//							//neighbor before (tested here) and after this (the ineigh)
//							//corner (i.e. the ineigh neighbor) are the same generation
//							//as me, therefor this (the ineigh) node is a CORNER and not
//							//an S_C_CON node
//							inode = ineigh;
//							NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//							assert(NdTemp);
//							NdTemp->putinfo(CORNER);
//						}
//
//						inode = ineigh + 4;
//						NdTemp = (Node*) NodeTable->lookup(EmSon[isonA]->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(SIDE);
//
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(CORNER);
//
//						NdTemp = (Node*) NodeTable->lookup(EmSon[isonB]->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(SIDE);
//
//						if (EmSon[isonB]->neigh_gen[(ineigh + 1) % 4] == EmSon[isonB]->generation) {
//							//neighbor before (i.e. the ineigh neighbor) and after
//							//(tested here) this (the (ineigh+1)%4) corner are the
//							//the same generation as me, therefore this (the
//							//(ineigh+1)%4) node is a CORNER and not an S_C_CON node
//							inode = (ineigh + 1) % 4;
//							NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//							assert(NdTemp);
//							NdTemp->putinfo(CORNER);
//						}
//
//						break;
//					case 3:
//						//case E
//
//						//update the nodes on this side
//
//						inode = ineigh; //father corner node
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(CORNER);
//
//						inode = ineighp4; //father edge node
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(CORNER);
//
//						inode = (ineigh + 1) % 4; //father corner node
//						NdTemp = (Node*) NodeTable->lookup(EmFather->node_key[inode]);
//						assert(NdTemp);
//						NdTemp->putinfo(CORNER);
//
//						for (ikey = 0; ikey < KEYLENGTH; ikey++) {
//
//							EmNeighNew[0]->neighbor[ineighme][ikey] = EmNeighNew[0]->neighbor[ineighmep4][ikey] =
//							    EmNeighNew[1]->neighbor[ineighme][ikey] =
//							        EmNeighNew[1]->neighbor[ineighmep4][ikey] = EmSon[isonA]->key[ikey];
//
//							EmNeighNew[2]->neighbor[ineighme][ikey] = EmNeighNew[2]->neighbor[ineighmep4][ikey] =
//							    EmNeighNew[3]->neighbor[ineighme][ikey] =
//							        EmNeighNew[3]->neighbor[ineighmep4][ikey] = EmSon[isonB]->key[ikey];
//
//							EmSon[isonA]->neighbor[ineigh][ikey] = EmNeighNew[0]->key[ikey];
//							EmSon[isonA]->neighbor[ineighp4][ikey] = EmNeighNew[1]->key[ikey];
//
//							EmSon[isonB]->neighbor[ineigh][ikey] = EmNeighNew[2]->key[ikey];
//							EmSon[isonB]->neighbor[ineighp4][ikey] = EmNeighNew[3]->key[ikey];
//						}
//
//						//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//
//						EmNeighNew[0]->neigh_gen[ineighme] = EmNeighNew[0]->neigh_gen[ineighmep4] =
//						    EmNeighNew[1]->neigh_gen[ineighme] = EmNeighNew[1]->neigh_gen[ineighmep4] =
//						        EmNeighNew[2]->neigh_gen[ineighme] = EmNeighNew[2]->neigh_gen[ineighmep4] =
//						            EmNeighNew[3]->neigh_gen[ineighme] = EmNeighNew[3]->neigh_gen[ineighmep4] =
//						                EmSon[isonA]->generation;
//
//						EmSon[isonA]->neigh_gen[ineigh] = EmSon[isonA]->neigh_gen[ineighp4] =
//						    EmNeighNew[0]->generation;
//
//						EmSon[isonB]->neigh_gen[ineigh] = EmSon[isonB]->neigh_gen[ineighp4] =
//						    EmNeighNew[2]->generation;
//
//						EmNeighNew[0]->neigh_proc[ineighmep4] = EmNeighNew[1]->neigh_proc[ineighmep4] =
//						    EmNeighNew[2]->neigh_proc[ineighmep4] = EmNeighNew[3]->neigh_proc[ineighmep4] = -2;
//
//						//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//
//						EmSon[isonA]->neigh_proc[ineigh] = EmFather->neigh_proc[ineigh];
//
//						//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//
//						EmSon[isonB]->neigh_proc[ineigh] = EmFather->neigh_proc[ineighp4];
//
//						//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//
//						inode = ineighp4; //sonA edge node
//						NdTemp = (Node*) NodeTable->lookup(EmSon[isonA]->node_key[inode]);
//						assert(NdTemp);
//						if (compare_key(EmSon[isonA]->neighbor[ineigh], EmSon[isonA]->neighbor[ineighp4])) {
//							EmSon[isonA]->neigh_proc[ineighp4] = -2;
//							NdTemp->putinfo(SIDE);
//						} else {
//							EmSon[isonA]->neigh_proc[ineighp4] = EmSon[isonA]->neigh_proc[ineigh];
//
//							NdTemp->putinfo(S_C_CON);
//
//							inode = ineighmep4;
//
//							NdTemp = (Node*) NodeTable->lookup(EmNeighNew[0]->node_key[inode]);
//							assert(NdTemp);
//							NdTemp->putinfo(S_S_CON);
//
//							NdTemp = (Node*) NodeTable->lookup(EmNeighNew[1]->node_key[inode]);
//							assert(NdTemp);
//							NdTemp->putinfo(S_S_CON);
//						}
//
//						//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
//						inode = ineighp4; //sonB edge node
//						NdTemp = (Node*) NodeTable->lookup(EmSon[isonB]->node_key[inode]);
//						assert(NdTemp);
//						if (compare_key(EmSon[isonB]->neighbor[ineigh], EmSon[isonB]->neighbor[ineighp4])) {
//							EmSon[isonB]->neigh_proc[ineighp4] = -2;
//							NdTemp->putinfo(SIDE);
//						} else {
//							EmSon[isonB]->neigh_proc[ineighp4] = EmSon[isonB]->neigh_proc[ineigh];
//
//							NdTemp->putinfo(S_C_CON);
//
//							inode = ineighmep4;
//
//							NdTemp = (Node*) NodeTable->lookup(EmNeighNew[2]->node_key[inode]);
//							assert(NdTemp);
//							NdTemp->putinfo(S_S_CON);
//
//							NdTemp = (Node*) NodeTable->lookup(EmNeighNew[3]->node_key[inode]);
//							assert(NdTemp);
//							NdTemp->putinfo(S_S_CON);
//						}
//
//						break;
//					default:
//						printf("FUBAR 2 detected in refine_neigh_update! aborting.\n");
//						assert(0);
//
//						break;
//				} //switch(inewcase), case based on generation of my new neighbor
//
//			} //else: not a map boundary
//
//		} //iside loop
//
//	} //ifather loop

	/*************************************************************/
	/* The interprocessor update information should be here by   */
	/* now or at the very least I won't have to weight very long */
	/* receive neighbor update information from other processors */
	/*************************************************************/

//	vector<unsigned> neighb;
//	vector<int> casevec;
	vector<ErrorElem*> old_neigh_father;
	Element* EmFather;

	for (int ifather = 0; ifather < RefinedList->get_num_elem(); ifather++) {
		EmFather = RefinedList->get(ifather);
		for (int i = 0; i < 8; ++i)

			if (EmFather->neigh_proc[i] == INIT) {
				assert(i < 4);

				int son1_ind = i, son2_ind = (i + 1) % 4;
				Element* son1 = (Element*) El_Table->lookup(EmFather->son[son1_ind]);
				Element* son2 = (Element*) El_Table->lookup(EmFather->son[son2_ind]);
				int neigh1 = i, neigh2 = i + 4;

				for (int key = 0; key < KEYLENGTH; ++key) {
					son1->neighbor[neigh1][key] = 0;
					son1->neighbor[neigh2][key] = 0;
					son2->neighbor[neigh1][key] = 0;
					son2->neighbor[neigh2][key] = 0;
				}
				son1->neigh_gen[neigh1] = son1->neigh_gen[neigh2] = son2->neigh_gen[neigh1] =
				    son2->neigh_gen[neigh2] = 0;
				son1->neigh_proc[neigh1] = son2->neigh_proc[neigh1] = INIT;
				son1->neigh_proc[neigh2] = son2->neigh_proc[neigh2] = -2;
				assert(son1->generation == EmFather->generation + 1);
				assert(son2->generation == EmFather->generation + 1);

//				Node* nodeTemp = (Node*) NodeTable->lookup(EmFather->key);
//				assert(nodeTemp);
//				nodeTemp->info = CORNER;
//
//				nodeTemp = (Node*) NodeTable->lookup(son1->key);
//				assert(nodeTemp);
//				nodeTemp->info = BUBBLE;
//
//				nodeTemp = (Node*) NodeTable->lookup(son2->key);
//				assert(nodeTemp);
//				nodeTemp->info = BUBBLE;
//
//				nodeTemp = (Node*) NodeTable->lookup(EmFather->node_key[i + 4]);
//				assert(nodeTemp);
//				nodeTemp->info = CORNER;
//
//				nodeTemp = (Node*) NodeTable->lookup(son1->node_key[i + 4]);
//				assert(nodeTemp);
//				nodeTemp->info = SIDE;
//
//				nodeTemp = (Node*) NodeTable->lookup(son2->node_key[i + 4]);
//				assert(nodeTemp);
//				nodeTemp->info = SIDE;

			} else if (EmFather->neigh_proc[i] >= 0 /*&& EmFather->neigh_proc[i] != myid*/) {
				// this old neighbor should be deleted at the end
				ErrorElem* old_neigh = (ErrorElem*) El_Table->lookup(EmFather->neighbor[i]);
				unsigned neigh_sonkeys[4][2];
				old_neigh->gen_my_sons_key(El_Table, neigh_sonkeys);
				int difgen = EmFather->generation - old_neigh->generation;
				switch (difgen) {
					case -1: {
						//this means that I have 8 neighbors
						int son = ((i < 4) ? i : (i + 1) % 4);
						Element* son_elem = (Element*) El_Table->lookup(EmFather->son[son]);
						assert(son_elem);
						int neigh1 = i % 4;
						int neigh2 = i % 4 + 4;
						for (int key = 0; key < KEYLENGTH; ++key) {
							son_elem->neighbor[neigh1][key] = neigh_sonkeys[(neigh1 + 3) % 4][key];
							son_elem->neighbor[neigh2][key] = neigh_sonkeys[(neigh2 + 2) % 4][key];
						}

						son_elem->neigh_proc[neigh2] = son_elem->neigh_proc[neigh1] = old_neigh->myprocess;
						assert(son_elem->generation == EmFather->generation + 1);
						assert(son_elem->generation == EmFather->generation + 1);
						son_elem->neigh_gen[neigh1] = old_neigh->generation + 1;
						son_elem->neigh_gen[neigh2] = old_neigh->generation + 1;

//						for (int key = 0; key < KEYLENGTH; ++key)
//							neighb.push_back(son_elem->key[key]);
//
//						for (int key = 0; key < KEYLENGTH; ++key)
//							neighb.push_back(neigh_sonkeys[(neigh1 + 3) % 4][key]);
//						casevec.push_back(-1);
//
//						for (int key = 0; key < KEYLENGTH; ++key)
//							neighb.push_back(son_elem->key[key]);
//
//						for (int key = 0; key < KEYLENGTH; ++key)
//							neighb.push_back(neigh_sonkeys[(neigh2 + 2) % 4][key]);
//						casevec.push_back(-1);

//						Node* nodeTemp = (Node*) NodeTable->lookup(EmFather->key);
//						assert(nodeTemp);
//						nodeTemp->info = CORNER;
//
//						nodeTemp = (Node*) NodeTable->lookup(EmFather->node_key[(i % 4) + 4]);
//						assert(nodeTemp);
//						nodeTemp->info = CORNER;
//
//						nodeTemp = (Node*) NodeTable->lookup(son_elem->key);
//						assert(nodeTemp);
//						nodeTemp->info = BUBBLE;
//
//						nodeTemp = (Node*) NodeTable->lookup(son_elem->node_key[(i % 4) + 4]);
//						assert(nodeTemp);
//						nodeTemp->info = S_C_CON;
						break;
					}
					case 0: {
						// here we set both neighbors
						if (i < 4) {
							int son1_ind = i, son2_ind = (i + 1) % 4;
							Element* son1 = (Element*) El_Table->lookup(EmFather->son[son1_ind]);
							Element* son2 = (Element*) El_Table->lookup(EmFather->son[son2_ind]);
							int neigh1 = i, neigh2 = i + 4;

							for (int key = 0; key < KEYLENGTH; ++key) {
								son1->neighbor[neigh1][key] = neigh_sonkeys[(son1_ind + 3) % 4][key];
								son1->neighbor[neigh2][key] = neigh_sonkeys[(son1_ind + 3) % 4][key];
								son2->neighbor[neigh1][key] = neigh_sonkeys[(son2_ind + 1) % 4][key];
								son2->neighbor[neigh2][key] = neigh_sonkeys[(son2_ind + 1) % 4][key];
							}

							assert(son1->generation == EmFather->generation + 1);
							assert(son2->generation == EmFather->generation + 1);
							son1->neigh_gen[neigh1] = old_neigh->generation + 1;
							son1->neigh_gen[neigh2] = old_neigh->generation + 1;
							son2->neigh_gen[neigh1] = old_neigh->generation + 1;
							son2->neigh_gen[neigh2] = old_neigh->generation + 1;

							son1->neigh_proc[neigh1] = old_neigh->myprocess;
							son1->neigh_proc[neigh2] = -2;
							son2->neigh_proc[neigh1] = old_neigh->myprocess;
							son2->neigh_proc[neigh2] = -2;

//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son1->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[(son1_ind + 3) % 4][key]);
//							casevec.push_back(0);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son1->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[(son1_ind + 3) % 4][key]);
//							casevec.push_back(0);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son2->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[(son2_ind + 1) % 4][key]);
//							casevec.push_back(0);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son2->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[(son2_ind + 1) % 4][key]);
//							casevec.push_back(0);

//							Node* nodeTemp = (Node*) NodeTable->lookup(EmFather->key);
//							assert(nodeTemp);
//							nodeTemp->info = CORNER;
//
//							///?????????????????????????????
//							nodeTemp = (Node*) NodeTable->lookup(EmFather->node_key[i + 4]);
//							assert(nodeTemp);
//							nodeTemp->info = CORNER;
//
//							nodeTemp = (Node*) NodeTable->lookup(son1->key);
//							assert(nodeTemp);
//							nodeTemp->info = BUBBLE;
//
//							nodeTemp = (Node*) NodeTable->lookup(son2->key);
//							assert(nodeTemp);
//							nodeTemp->info = BUBBLE;
//
//							///?????????????????????????????
//							nodeTemp = (Node*) NodeTable->lookup(son1->node_key[i + 4]);
//							assert(nodeTemp);
//							nodeTemp->info = SIDE;
//
//							///?????????????????????????????
//							nodeTemp = (Node*) NodeTable->lookup(son2->node_key[i + 4]);
//							assert(nodeTemp);
//							nodeTemp->info = SIDE;
						}
						break;
					}
					case 1: {
						if (i < 4) {
							assert(EmFather->neigh_proc[i + 4] == -2);
							int son1_ind = i, son2_ind = (i + 1) % 4;
							Element* son1 = (Element*) El_Table->lookup(EmFather->son[son1_ind]);
							Element* son2 = (Element*) El_Table->lookup(EmFather->son[son2_ind]);
							assert(son1 && son2);
							int neigh1 = i, neigh2 = i + 4;

							int w_neigh = old_neigh->which_neighbor(EmFather->key);
							assert(old_neigh->neigh_proc[(w_neigh + 4) % 4] != -2);

							int son_neigh = (w_neigh < 4) ? w_neigh : (w_neigh + 1) % 4;

							for (int key = 0; key < KEYLENGTH; ++key) {
								son1->neighbor[neigh1][key] = neigh_sonkeys[son_neigh][key];
								son1->neighbor[neigh2][key] = neigh_sonkeys[son_neigh][key];
								son2->neighbor[neigh1][key] = neigh_sonkeys[son_neigh][key];
								son2->neighbor[neigh2][key] = neigh_sonkeys[son_neigh][key];
							}

							assert(son1->generation == EmFather->generation + 1);
							assert(son2->generation == EmFather->generation + 1);
							son1->neigh_gen[neigh1] = old_neigh->generation + 1;
							son1->neigh_gen[neigh2] = old_neigh->generation + 1;
							son2->neigh_gen[neigh1] = old_neigh->generation + 1;
							son2->neigh_gen[neigh2] = old_neigh->generation + 1;

							son1->neigh_proc[neigh1] = old_neigh->myprocess;
							son1->neigh_proc[neigh2] = -2;
							son2->neigh_proc[neigh1] = old_neigh->myprocess;
							son2->neigh_proc[neigh2] = -2;

//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son1->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[son_neigh][key]);
//							casevec.push_back(1);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son1->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[son_neigh][key]);
//							casevec.push_back(1);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son2->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[son_neigh][key]);
//							casevec.push_back(1);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(son2->key[key]);
//
//							for (int key = 0; key < KEYLENGTH; ++key)
//								neighb.push_back(neigh_sonkeys[son_neigh][key]);
//							casevec.push_back(1);

//							Node* nodeTemp = (Node*) NodeTable->lookup(EmFather->key);
//							assert(nodeTemp);
//							nodeTemp->info = CORNER;
//
//							nodeTemp = (Node*) NodeTable->lookup(EmFather->node_key[i + 4]);
//							assert(nodeTemp);
//							nodeTemp->info = S_C_CON;
//
//							nodeTemp = (Node*) NodeTable->lookup(son1->key);
//							assert(nodeTemp);
//							nodeTemp->info = BUBBLE;
//
//							nodeTemp = (Node*) NodeTable->lookup(son2->key);
//							assert(nodeTemp);
//							nodeTemp->info = BUBBLE;
//
//							nodeTemp = (Node*) NodeTable->lookup(son1->node_key[i + 4]);
//							assert(nodeTemp);
//							nodeTemp->info = S_S_CON;
//
//							nodeTemp = (Node*) NodeTable->lookup(son2->node_key[i + 4]);
//							assert(nodeTemp);
//							nodeTemp->info = S_S_CON;
//
//							nodeTemp = (Node*) NodeTable->lookup(son1->node_key[i]);
//							assert(nodeTemp);
//							nodeTemp->info = CORNER;

						}
						break;
					}
					default:
						assert(0);

				}
			}
	}

//	int recv, sentid, recvid;
//	int mysend = neighb.size();
//	if (myid == 0)
//		sentid = recvid = 1;
//	else
//		sentid = recvid = 0;
//	MPI_Status status;
//
//	MPI_Send(&mysend, 1, MPI_INT, sentid, 1000, MPI_COMM_WORLD);
//	MPI_Recv(&recv, 1, MPI_INT, recvid, 1000, MPI_COMM_WORLD, &status);
//
//	unsigned* check = new unsigned[recv];
//	int* check_case = new int[recv / 4];
//
//	MPI_Send(&neighb[0], mysend, MPI_UNSIGNED, sentid, 1000, MPI_COMM_WORLD);
//	MPI_Recv(check, recv, MPI_UNSIGNED, recvid, 1000, MPI_COMM_WORLD, &status);
//
//	MPI_Send(&casevec[0], mysend / 4, MPI_INT, sentid, 1000, MPI_COMM_WORLD);
//	MPI_Recv(check_case, recv / 4, MPI_INT, recvid, 1000, MPI_COMM_WORLD, &status);
//
//	for (int i = 0; i < recv; i = i + 4) {
//		unsigned key[] = { check[i], check[i + 1] };
//		unsigned key_neighb[] = { check[i + 2], check[i + 3] };
//		Element* myelm = (Element*) El_Table->lookup(key_neighb);
//
//		int neigh = myelm->which_neighbor(key);
//		cout << check_case[i / 4] << endl;
//
//		assert(myelm->neigh_proc[neigh] != myid && myelm->neigh_proc[neigh] >= 0);
//	}

//	for (int i = 0; i < old_neigh_father.size(); ++i) {
//		ErrorElem* oldneigh = old_neigh_father[i]; //Hello I'm the OLDFATHER
//		assert(oldneigh);
//		El_Table->remove(oldneigh->key);
//		delete oldneigh;
//	}
	int count = 0;
	for (int ifather = 0; ifather < RefinedList->get_num_elem(); ifather++) {
		EmFather = RefinedList->get(ifather); //Hello I'm the OLDFATHER
		assert(EmFather); //Help I've been abducted call the FBI!!!
		assert(EmFather->adapted==OLDFATHER); //sanity check
		assert(EmFather->refined == 1);
		EmFather->adapted = TOBEDELETED; //I've lived a good life, it's my time to die

		for (int i = 0; i < 8; ++i) {
			if (EmFather->neigh_proc[i] != myid && EmFather->neigh_proc[i] >= 0) {
				Element* oldneigh = (Element*) El_Table->lookup(EmFather->neighbor[i]);
				if (oldneigh) {
					El_Table->remove(oldneigh->key, 1, stdout, myid, 15);
					delete oldneigh;
				}
			}
		}
		count++;

		El_Table->remove(EmFather->key, 1, stdout, myid, 15);
		delete EmFather;
	}

	RefinedList->trashlist();

}

void adjust_node_info(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* Node_Table = meshctx->nd_table;
	int myid = propctx->myid, side;

	int NodeTable_num_buck = Node_Table->get_no_of_buckets();
	HashEntryPtr *NodeTable_bucket0 = Node_Table->getbucketptr();
	HashEntryPtr NodeTable_entry_ptr;
	int inodebucket;

//zero the number of elems each node is associated with
	for (inodebucket = 0; inodebucket < NodeTable_num_buck; inodebucket++) {
		NodeTable_entry_ptr = *(NodeTable_bucket0 + inodebucket);

		while (NodeTable_entry_ptr) {

			Node *NdTemp = (Node*) (NodeTable_entry_ptr->value);
			NodeTable_entry_ptr = NodeTable_entry_ptr->next;
			assert(NdTemp);
			NdTemp->info = -1;
		}
	}

	HashEntryPtr currentPtr;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Element *Curr_El = (Element*) (currentPtr->value);
				Node* node = (Node*) Node_Table->lookup(Curr_El->key);
				assert(node);
				node->info = BUBBLE;

				for (int side = 0; side < 4; ++side) {
					node = (Node*) Node_Table->lookup(Curr_El->node_key[side]);
					assert(node);
					node->info = CORNER;
				}

				for (int side = 4; side < 8; ++side) {
					node = (Node*) Node_Table->lookup(Curr_El->node_key[side]);
					assert(node);
					node->info = SIDE;
				}

				currentPtr = currentPtr->next;
			}
		}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Element *Curr_El = (Element*) (currentPtr->value);
				Node* node;

				for (int side = 0; side < 4; ++side) {
					if (Curr_El->generation > Curr_El->neigh_gen[side] && Curr_El->neigh_proc[side] != INIT) {
						node = (Node*) Node_Table->lookup(Curr_El->node_key[side + 4]);
						assert(node);
						node->info = S_S_CON;
					}
					if (Curr_El->generation < Curr_El->neigh_gen[side] && Curr_El->neigh_proc[side] != INIT) {
						node = (Node*) Node_Table->lookup(Curr_El->node_key[side + 4]);
						assert(node);
						node->info = S_C_CON;
					}
				}

				currentPtr = currentPtr->next;
			}
		}

//	HashEntryPtr currentPtr;
//	HashEntryPtr *buck = El_Table->getbucketptr();
//
//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Element *Curr_El = (Element*) (currentPtr->value);
//				Node* node = (Node*) Node_Table->lookup(Curr_El->key);
//				assert(node);
//				node->info == BUBBLE;
//				for (int side = 0; side < 4; ++side) {
//					if (Curr_El->neigh_proc[side] == INIT) {
//						node = (Node*) Node_Table->lookup(Curr_El->node_key[side]);
//						assert(node);
//						if (node->info != S_C_CON)
//							node->info = CORNER;
//
//						node = (Node*) Node_Table->lookup(Curr_El->node_key[(side + 1) % 4]);
//						assert(node);
//						if (node->info != S_C_CON)
//							node->info = CORNER;
//
//						node = (Node*) Node_Table->lookup(Curr_El->node_key[side + 4]);
//						assert(node);
//						node->info = SIDE;
//
//						// the following condition means that it is ok if element itself is not
//						// in this processor but the side that we are studding must be in this processor
//					} else if (Curr_El->myprocess == myid || Curr_El->neigh_proc[side] == myid) {
//						Element *neigh_elem = (Element*) El_Table->lookup(Curr_El->neighbor[side]);
//						assert(neigh_elem);
//						switch (Curr_El->generation - neigh_elem->generation) {
//							case -1:
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[side + 4]);
//								assert(node);
//								node->info = S_C_CON;
//
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[side]);
//								assert(node);
//								if (node->info != S_C_CON)
//									node->info = CORNER;
//
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[(side + 1) % 4]);
//								assert(node);
//								if (node->info != S_C_CON)
//									node->info = CORNER;
//
//								break;
//							case 0:
//
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[side + 4]);
//								assert(node);
//								node->info = SIDE;
//
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[side]);
//								assert(node);
//								if (node->info != S_C_CON)
//									node->info = CORNER;
//
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[(side + 1) % 4]);
//								assert(node);
//								if (node->info != S_C_CON)
//									node->info = CORNER;
//
//								assert(Curr_El->neigh_proc[side + 4] == -2);
//								assert(
//								    Curr_El->neigh_gen[side + 4] == Curr_El->neigh_gen[side]
//								        && Curr_El->neigh_gen[side + 4] == Curr_El->generation);
//
//								break;
//							case 1:
//
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[side + 4]);
//								assert(node);
//								node->info = S_S_CON;
//
//								assert(Curr_El->neigh_proc[side + 4] == -2);
//
////								node = (Node*) Node_Table->lookup(Curr_El->node_key[side ]);
////								assert(node);
////								node->info == CORNER;
////
////								node = (Node*) Node_Table->lookup(Curr_El->node_key[(side+1)%4 ]);
////								assert(node);
////								node->info == CORNER;
//
//								break;
//							default:
//								assert(0);
//								break;
//
//						}
//
//					}
//
//				}
//
//				currentPtr = currentPtr->next;
//			}
//		}

//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Element *Curr_El = (Element*) (currentPtr->value);
//				if (Curr_El->get_adapted_flag() > 0) {
//					for (int j = 0; j < 8; ++j)
//						if (Curr_El->neigh_proc[j] >= 0 && Curr_El->neigh_proc[j] != myid) {
//							Element *neigh_elem = (Element*) El_Table->lookup(Curr_El->neighbor[j]);
//							assert(neigh_elem);
//							assert(neigh_elem->get_adapted_flag() <= 0);
//							Node* node;
//							for (int k = 0; k < 4; ++k) {
//								node = (Node*) Node_Table->lookup(Curr_El->node_key[k]);
//								node->info == CORNER;
//								node = (Node*) Node_Table->lookup(neigh_elem->node_key[k]);
//								node->info == CORNER;
//							}
//							switch (Curr_El->generation - Curr_El->neigh_gen[j]) {
//								case -1:
//									node = (Node*) Node_Table->lookup(Curr_El->node_key[(j % 4) + 4]);
//									assert(node);
//									node->info = S_C_CON;
//
//									node = (Node*) Node_Table->lookup(neigh_elem->node_key[(j + 2) % 4 + 4]);
//									assert(node);
//									node->info = S_S_CON;
//									break;
//								case 0:
//									node = (Node*) Node_Table->lookup(Curr_El->node_key[(j % 4) + 4]);
//									assert(node);
//									node->info = SIDE;
//
//									assert(
//									    compare_key(Curr_El->node_key[j % 4 + 4],
//									        neigh_elem->node_key[(j + 2) % 4 + 4]));
//
//									break;
//								case 1:
//									node = (Node*) Node_Table->lookup(Curr_El->node_key[(j % 4) + 4]);
//									assert(node);
//									node->info = S_S_CON;
//
//									node = (Node*) Node_Table->lookup(neigh_elem->node_key[(j + 2) % 4 + 4]);
//									assert(node);
//									node->info = S_C_CON;
//
//									break;
//								default:
//									assert(0);
//									break;
//
//							}
//
//						}
//
//				}
//
//				currentPtr = currentPtr->next;
//			}
//		}

}
