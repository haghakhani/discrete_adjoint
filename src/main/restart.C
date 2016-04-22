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
#include <set>

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
	gzwrite(myfile, &(numnode), sizeof(int));
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
	gzwrite(myfile, &(numelem), sizeof(int));
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

	// now writing outline

	//output maximum flow depth a.k.a. flow outline
	OutLine outline2;
	double dxy[2];
	dxy[0] = outline->dx;
	dxy[1] = outline->dy;
	outline2.init2(dxy, outline->xminmax, outline->yminmax);
	int NxNyout = outline->Nx * outline->Ny;
	MPI_Allreduce(*(outline->pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE,
	MPI_SUM, MPI_COMM_WORLD);

	// we write every single file to when we want to load can read from the file concurrently
	sprintf(filename, "outline_%04d", myid);
	myfile = gzopen(filename, "wb");
	int nx = outline->Nx, ny = outline->Ny;
	gzwrite(myfile, dxy, sizeof(double) * 2);
	gzwrite(myfile, outline->xminmax, sizeof(double) * 2);
	gzwrite(myfile, outline->yminmax, sizeof(double) * 2);
	for (int i = 0; i < ny; ++i)
		gzwrite(myfile, outline2.pileheight[i], sizeof(double) * nx);
	gzclose(myfile);

	MPI_Barrier(MPI_COMM_WORLD);
}

int loadrun(int myid, int numprocs, HashTable** NodeTable, HashTable** ElemTable, SolRec** solrec,
    MatProps* matprops, TimeProps* timeprops, OutLine* outline) {

	char filename[50];
	sprintf(filename, "restart_%04d", myid);
	//printf("filename=\"%s\"\n",filename);
	gzFile myfile = gzopen(filename, "rb");

	if (myfile == NULL)
		return (0);

	timeprops->read_from_file(myfile);

	//recreate the node hashtable
	*NodeTable = new HashTable(myfile);

	int numnode;
	gzread(myfile, (void*) &(numnode), sizeof(int));

	//read in all the nodes
	for (int inode = 0; inode < numnode; inode++) {
		Node* NodeP = new Node(myfile, matprops);
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
		ElemP = new Element(myfile, *NodeTable, matprops, myid);
		(*ElemTable)->add(ElemP->pass_key(), ElemP);
		if (ElemP->get_gen() > maxgen)
			maxgen = ElemP->get_gen();
	}

	double dx = *(ElemP->get_dx() + 0), dy = *(ElemP->get_dx() + 1);
	if (dx < dy)
		dx = dy;

	REFINE_LEVEL = ElemP->get_gen()
	    + ceil(
	        log(dx * (matprops->number_of_cells_across_axis) / (matprops->smallest_axis)) / log(2.0));

	if (maxgen > REFINE_LEVEL)
		REFINE_LEVEL = maxgen;

	maxgen = REFINE_LEVEL;
	MPI_Allreduce(&maxgen, &REFINE_LEVEL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	*solrec = new SolRec(myfile);

	gzclose(myfile);

	move_data(numprocs, myid, *ElemTable, *NodeTable, timeprops);
//	//calc_d_gravity
//	int no_of_elm_buckets = (*ElemTable)->get_no_of_buckets();
//	for (int ibucket = 0; ibucket < no_of_elm_buckets; ibucket++) {
//		HashEntryPtr entryp = *((*ElemTable)->getbucketptr() + ibucket);
//		while (entryp) {
//			Element* EmTemp = (Element*) entryp->value;
//			assert(EmTemp);
//
//			if (EmTemp->get_adapted_flag() > 0) {
////				EmTemp->calc_which_son();
////				EmTemp->find_positive_x_side(*NodeTable);
////				EmTemp->calculate_dx(*NodeTable);
////				EmTemp->calc_topo_data(matprops);
////				EmTemp->calc_gravity_vector(matprops);
//				EmTemp->calc_d_gravity(*ElemTable);
//				EmTemp->find_opposite_brother(*ElemTable);
//			}
//			entryp = entryp->next;
//		}
//	}
//
//	move_data(numprocs, myid, *ElemTable, *NodeTable, timeprops);

	// now reading outline

	sprintf(filename, "outline_%04d", myid);
	myfile = gzopen(filename, "rb");
	int nx = outline->Nx, ny = outline->Ny;
	double dxy[2], xminmax[2], yminmax[2];
	gzread(myfile, dxy, sizeof(double) * 2);
	gzread(myfile, xminmax, sizeof(double) * 2);
	gzread(myfile, yminmax, sizeof(double) * 2);
	outline->init2(dxy, xminmax, yminmax);
	for (int i = 0; i < ny; ++i)
		gzread(myfile, outline->pileheight[i], sizeof(double) * nx);
	gzclose(myfile);

	MPI_Barrier(MPI_COMM_WORLD);

	return (1);
}

class Data {
private:
	unsigned* key;
	double* state;

public:

	Data(Element* elem) {
		key = elem->pass_key();
		state = elem->get_state_vars();
	}

	unsigned* get_key() const {
		return key;
	}

	double* get_state() const {
		return state;
	}

	bool operator<(const Data& rdata) const {
		if (key[0] < rdata.key[0] || (key[0] == rdata.key[0] && key[1] < rdata.key[1]))
			return true;

		return false;
	}
	;
};

void write_alldata_ordered(HashTable* El_Table, int myid) {

	set<Data> mydata;

	int no_of_elm_buckets = El_Table->get_no_of_buckets();
	for (int ibucket = 0; ibucket < no_of_elm_buckets; ibucket++) {
		HashEntryPtr entryp = *(El_Table->getbucketptr() + ibucket);
		while (entryp) {
			Element* EmTemp = (Element*) entryp->value;
			if (EmTemp->get_adapted_flag() > 0)
				mydata.insert(Data(EmTemp));

			entryp = entryp->next;
		}
	}

	char filename[50];
	sprintf(filename, "ordered_%04d", myid);
//	gzFile myfile = gzopen(filename, "wb");
	FILE *fp = fopen(filename, "w");

	set<Data>::iterator it;
	for (it = mydata.begin(); it != mydata.end(); ++it) {
		fprintf(fp, "%u %u %16.10f %16.10f %16.10f \n", it->get_key()[0], it->get_key()[1],
		    it->get_state()[0], it->get_state()[1], it->get_state()[2]);
//		gzwrite(myfile, (it->get_key()), sizeof(unsigned) * 2);
//		gzwrite(myfile, (it->get_state()), sizeof(double) * 3);
	}

	fclose(fp);

//	gzclose (myfile);

}

