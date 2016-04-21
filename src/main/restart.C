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
 * $Id: restart.C,v 1.4 2004/08/11 15:54:22 kdalbey Exp $ 
 */

#include "../header/hpfem.h"
//#define DEBUGSAVESPLIT
//#define DEBUGSAVEHEADER
#define NUM_CHAR_IN_SAVE_HEADER 16384 //equiv to 4096 integers of space available for save header

void save_forward(const MeshCTX& meshctx, const PropCTX& propctx, SolRec *solrec) {

	TimeProps* timeprops = propctx.timeprops;
	// we do not need to write them we can re-make them when we restart
//	MatProps* matprops = propctx.matprops;
//	MapNames* mapnames = propctx.mapnames;
	OutLine* outline = propctx.outline;
	int numprocs = propctx.numproc;
	int myid = propctx.myid;
	int adapt_flag = propctx.adapt_flag;

	HashTable* El_Table = meshctx.el_table;
	HashTable* NodeTable = meshctx.nd_table;
	char filename[50];

	move_data(numprocs, myid, El_Table, NodeTable, timeprops);

	sprintf(filename, "restart_%04d", myid);
	gzFile myfile = gzopen(filename, "wb");

	timeprops->wrtie_to_file(myfile);

	NodeTable->write_table(myfile);
	int numnode = table_members(NodeTable);
	gzwrite(myfile, (void*) &(numnode), sizeof(int));
	cout << "num node " << numnode << endl;

	int count = 0;
	HashEntryPtr * buck = NodeTable->getbucketptr();
	for (int i = 0; i < NodeTable->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Node *Curr_Node = (Node*) (currentPtr->value);

				Curr_Node->write_node(myfile);
				count++;
				currentPtr = currentPtr->next;
			}
		}

	assert(count == numnode);

	El_Table->write_table(myfile);

	int numelem = table_members(El_Table);
	gzwrite(myfile, (void*) &(numelem), sizeof(int));
	cout << "num elem " << numelem << endl;

	count = 0;
	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Element *Curr_El = (Element*) (currentPtr->value);

				Curr_El->write_elem(myfile);
				count++;
				currentPtr = currentPtr->next;
			}
		}

	assert(count == numelem);

	solrec->write_table(myfile);

	gzclose(myfile);
}

int loadrun(int myid, int numprocs, HashTable** NodeTable, HashTable** ElemTable, SolRec** solrec,
    MatProps* matprops_ptr, TimeProps* timeprops_ptr) {

	char filename[50];
	sprintf(filename, "restart_%04d", myid);
	//printf("filename=\"%s\"\n",filename);
	gzFile myfile = gzopen(filename, "rb");

	if (myfile == NULL)
		return (0);

	timeprops_ptr->read_from_file(myfile);

	//recreate the node hashtable
	*NodeTable = new HashTable(myfile);

	int numnode;
	gzread(myfile, (void*) &(numnode), sizeof(int));

	//read in all the nodes
	for (int inode = 0; inode < numnode; inode++) {
		Node* NodeP = new Node(myfile, matprops_ptr);
		(*NodeTable)->add(NodeP->pass_key(), NodeP);
	}

	//recreate the element hashtable
	*ElemTable = new HashTable(myfile);

	int numelem;
	gzread(myfile, (void*) &(numelem), sizeof(int));

	//read in all the elements
	int maxgen = 0;
	Element* ElemP;
	for (int ielem = 0; ielem < numelem; ielem++) {
		ElemP = new Element(myfile, *NodeTable, matprops_ptr, myid);
		(*ElemTable)->add(ElemP->pass_key(), ElemP);
		if (ElemP->get_gen() > maxgen)
			maxgen = ElemP->get_gen();
	}

	double dx = *(ElemP->get_dx() + 0), dy = *(ElemP->get_dx() + 1);
	if (dx < dy)
		dx = dy;

	REFINE_LEVEL = ElemP->get_gen()
	    + ceil(
	        log(dx * (matprops_ptr->number_of_cells_across_axis) / (matprops_ptr->smallest_axis))
	            / log(2.0));

	if (maxgen > REFINE_LEVEL)
		REFINE_LEVEL = maxgen;

	maxgen = REFINE_LEVEL;
	MPI_Allreduce(&maxgen, &REFINE_LEVEL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	*solrec = new SolRec(myfile);

	gzclose(myfile);

	//calc_d_gravity
	int no_of_elm_buckets = (*ElemTable)->get_no_of_buckets();
	for (int ibucket = 0; ibucket < no_of_elm_buckets; ibucket++) {
		HashEntryPtr entryp = *((*ElemTable)->getbucketptr() + ibucket);
		while (entryp) {
			Element* EmTemp = (Element*) entryp->value;
			assert(EmTemp);

			if (EmTemp->get_adapted_flag() > 0) {
				EmTemp->calc_d_gravity(*ElemTable);
				EmTemp->find_opposite_brother(*ElemTable);
			}
			entryp = entryp->next;
		}
	}

	move_data(numprocs, myid, *ElemTable, *NodeTable, timeprops_ptr);

	return (1);
}

