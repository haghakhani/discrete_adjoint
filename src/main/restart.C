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

	unsigned status = 1;

	gzwrite(myfile, &(status), sizeof(unsigned));
	gzwrite(myfile, &(min_dx[0]), sizeof(double));
	gzwrite(myfile, &(min_dx[1]), sizeof(double));
	gzwrite(myfile, &(min_gen), sizeof(double));

	timeprops->wrtie_to_file(myfile);

	NodeTable->write_table(myfile);
	int numnode = table_members(NodeTable);
	gzwrite(myfile, &(numnode), sizeof(int));
//	cout << "num node " << numnode << endl;

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

//	assert(count == numnode);

	El_Table->write_table(myfile);

	int numelem = table_members(El_Table);
	gzwrite(myfile, &(numelem), sizeof(int));
//	cout << "num elem " << numelem << endl;

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

//	assert(count == numelem);

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

void save_dual(const MeshCTX* meshctx, const MeshCTX* err_meshctx, const PropCTX* propctx,
    SolRec *solrec) {

	TimeProps* timeprops = propctx->timeprops;

	int numprocs = propctx->numproc;
	int myid = propctx->myid;
	int adapt_flag = propctx->adapt_flag;

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	HashTable* Err_El_Table = err_meshctx->el_table;
	HashTable* Err_NodeTable = err_meshctx->nd_table;

	char filename[50];

	sprintf(filename, "restart_%04d", myid);
	gzFile myfile = gzopen(filename, "wb");

	unsigned status = 2;

	gzwrite(myfile, &(status), sizeof(unsigned));
	gzwrite(myfile, &(min_dx[0]), sizeof(double));
	gzwrite(myfile, &(min_dx[1]), sizeof(double));
	gzwrite(myfile, &(min_gen), sizeof(double));

	gzwrite(myfile, &(FUNC_VAR[0]), 2 * sizeof(double));

	timeprops->wrtie_to_file(myfile);

	//Writing error node table
	NodeTable->write_table(myfile);
	int numnode = table_members(NodeTable);
	gzwrite(myfile, &(numnode), sizeof(int));
//	cout << "num node " << numnode << endl;

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
//	unsigned check = 1111;
//	gzwrite(myfile, &(check), sizeof(unsigned));

#ifdef Error
	//Writing error node table
	Err_NodeTable->write_table(myfile);
	numnode = table_members(Err_NodeTable);
	gzwrite(myfile, &(numnode), sizeof(int));
	//	cout << "num node " << numnode << endl;
	count = 0;
	buck = Err_NodeTable->getbucketptr();
	for (int i = 0; i < Err_NodeTable->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Node *Curr_Node = (Node*) (currentPtr->value);

				Curr_Node->write_node(myfile);
				count++;
				currentPtr = currentPtr->next;
			}
		}
//	check = 2222;
//	gzwrite(myfile, &(check), sizeof(unsigned));

//	assert(count == numnode);
#endif

	//write element table
	El_Table->write_table(myfile);

	int numelem = table_members(El_Table);
	gzwrite(myfile, &(numelem), sizeof(int));
//	cout << "num elem " << numelem << endl;

	count = 0;
	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem *Curr_El = (DualElem*) (currentPtr->value);

				Curr_El->write_elem(myfile);
				count++;
				currentPtr = currentPtr->next;
			}
		}
//	check = 3333;
//	gzwrite(myfile, &(check), sizeof(unsigned));

#ifdef Error
	//write error element table
	Err_El_Table->write_table(myfile);

	numelem = table_members(Err_El_Table);
	gzwrite(myfile, &(numelem), sizeof(int));
//	cout << "num elem " << numelem << endl;

	count = 0;
	buck = Err_El_Table->getbucketptr();
	for (int i = 0; i < Err_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				ErrorElem *Curr_El = (ErrorElem*) (currentPtr->value);

				Curr_El->write_elem(myfile);
				count++;
				currentPtr = currentPtr->next;
			}
		}

//	check = 4444;
//	gzwrite(myfile, &(check), sizeof(unsigned));
//	assert(count == numelem);
#endif

	solrec->write_table(myfile);

//	check = 5555;
//	gzwrite(myfile, &(check), sizeof(unsigned));

	gzclose(myfile);

	MPI_Barrier(MPI_COMM_WORLD);
}

int loadrun(int myid, int numprocs, HashTable** NodeTable, HashTable** ElemTable,
    HashTable** Err_NodeTable, HashTable** Err_ElemTable, SolRec** solrec, MatProps* matprops,
    TimeProps* timeprops, OutLine* outline) {

	char filename[50];
	sprintf(filename, "restart_%04d", myid);
	//printf("filename=\"%s\"\n",filename);
	gzFile myfile = gzopen(filename, "rb");

	if (myfile == NULL)
		return (0);

	// this is a flag that shows the place of restart. status=1 shows that the restart will take place in forward run and status=2
	unsigned status = 0;

	gzread(myfile, (void*) &(status), sizeof(unsigned));

	if (status == 1) {

		gzread(myfile, (void*) &(min_dx[0]), sizeof(double));
		gzread(myfile, (void*) &(min_dx[1]), sizeof(double));
		gzread(myfile, (void*) &(min_gen), sizeof(double));

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
		        log(dx * (matprops->number_of_cells_across_axis) / (matprops->smallest_axis))
		            / log(2.0));

		if (maxgen > REFINE_LEVEL)
			REFINE_LEVEL = maxgen;

		maxgen = REFINE_LEVEL;
		MPI_Allreduce(&maxgen, &REFINE_LEVEL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		*solrec = new SolRec(myfile);

		gzclose(myfile);

		move_data(numprocs, myid, *ElemTable, *NodeTable, timeprops);

		char filename[50];
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

	} else if (status == 2) {

		gzread(myfile, (void*) &(min_dx[0]), sizeof(double));
		gzread(myfile, (void*) &(min_dx[1]), sizeof(double));
		gzread(myfile, (void*) &(min_gen), sizeof(double));

		gzread(myfile, (void*) &(FUNC_VAR[0]), 2 * sizeof(double));

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

//		unsigned check;
//		gzread(myfile, &(check), sizeof(unsigned));
//		assert(check == 1111);

#ifdef Error
		//recreate error node hashtable
		*Err_NodeTable = new HashTable(myfile);

		gzread(myfile, (void*) &(numnode), sizeof(int));

		//read in all the nodes
		for (int inode = 0; inode < numnode; inode++) {
			Node* NodeP = new Node(myfile, matprops);
			(*Err_NodeTable)->add(NodeP->pass_key(), NodeP);
		}

//		gzread(myfile, &(check), sizeof(unsigned));
//		assert(check == 2222);

#endif

		//recreate the element hashtable
		*ElemTable = new HashTable(myfile);

		int numelem;
		gzread(myfile, (void*) &(numelem), sizeof(int));

		//read in all the elements
		int maxgen = 0;
		DualElem* ElemP;
		for (int ielem = 0; ielem < numelem; ielem++) {
			ElemP = new DualElem(myfile, *NodeTable, matprops, myid);

			(*ElemTable)->add(ElemP->pass_key(), ElemP);
			if (ElemP->get_gen() > maxgen)
				maxgen = ElemP->get_gen();
		}

//		gzread(myfile, &(check), sizeof(unsigned));
//		assert(check == 3333);

		double dx = *(ElemP->get_dx() + 0), dy = *(ElemP->get_dx() + 1);
		if (dx < dy)
			dx = dy;

		REFINE_LEVEL = ElemP->get_gen()
		    + ceil(
		        log(dx * (matprops->number_of_cells_across_axis) / (matprops->smallest_axis))
		            / log(2.0));

		if (maxgen > REFINE_LEVEL)
			REFINE_LEVEL = maxgen;

		maxgen = REFINE_LEVEL;
		MPI_Allreduce(&maxgen, &REFINE_LEVEL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		PropCTX propctx;
		propctx.myid = myid;
		propctx.numproc = numprocs;
		propctx.outline = outline;
		propctx.matprops = matprops;
		propctx.timeprops = timeprops;

		MeshCTX dualmesh;
		dualmesh.el_table = *ElemTable;
		dualmesh.nd_table = *NodeTable;
		move_dual_data(&dualmesh, &propctx);

#ifdef Error
		//recreate error element hashtable
		*Err_ElemTable = new HashTable(myfile);
		int err_numelem;

		gzread(myfile, (void*) &(err_numelem), sizeof(int));

		//read in all the elements
		ErrorElem* ElemPe;
		for (int ielem = 0; ielem < err_numelem; ielem++) {
			ElemPe = new ErrorElem(myfile, *Err_NodeTable, matprops, myid);
			(*Err_ElemTable)->add(ElemPe->pass_key(), ElemPe);
		}

//		gzread(myfile, &(check), sizeof(unsigned));
//		assert(check == 4444);

		MeshCTX errmesh;
		errmesh.el_table = *Err_ElemTable;
		errmesh.nd_table = *Err_NodeTable;

		move_err_data(&errmesh, &propctx);

		HashEntryPtr *buck = (*Err_ElemTable)->getbucketptr();
		for (int i = 0; i < (*Err_ElemTable)->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				HashEntryPtr currentPtr = *(buck + i);
				while (currentPtr) {
					ErrorElem *Curr_El = (ErrorElem*) (currentPtr->value);

					set_link(Curr_El, *ElemTable, *Err_ElemTable);
					currentPtr = currentPtr->next;
				}
			}

#endif

		*solrec = new SolRec(myfile, timeprops->iter, myid);

//		gzread(myfile, &(check), sizeof(unsigned));
//		assert(check == 5555);

		gzclose(myfile);

	}

	return (status);
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

class DualData: public Data {
private:
	double *adjoint;

public:
	DualData(DualElem* elem) :
			Data(elem) {
		adjoint = elem->get_adjoint();
	}

	double* get_adjoint() const {
		return adjoint;
	}
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

void write_alldualdata_ordered(HashTable* El_Table, int myid) {

	set<DualData> mydata;

	int no_of_elm_buckets = El_Table->get_no_of_buckets();
	for (int ibucket = 0; ibucket < no_of_elm_buckets; ibucket++) {
		HashEntryPtr entryp = *(El_Table->getbucketptr() + ibucket);
		while (entryp) {
			DualElem* EmTemp = (DualElem*) entryp->value;
			if (EmTemp->get_adapted_flag() > 0)
				mydata.insert(DualData(EmTemp));

			entryp = entryp->next;
		}
	}

	char filename[50];
	sprintf(filename, "ordered_%04d", myid);
//	gzFile myfile = gzopen(filename, "wb");
	FILE *fp = fopen(filename, "w");

	set<DualData>::iterator it;
	for (it = mydata.begin(); it != mydata.end(); ++it) {
		fprintf(fp, "%u %u %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f\n", it->get_key()[0],
		    it->get_key()[1], it->get_state()[0], it->get_state()[1], it->get_state()[2],
		    it->get_adjoint()[0], it->get_adjoint()[1], it->get_adjoint()[2]);
//		gzwrite(myfile, (it->get_key()), sizeof(unsigned) * 2);
//		gzwrite(myfile, (it->get_state()), sizeof(double) * 3);
	}

	fclose(fp);

//	gzclose (myfile);

}

