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

#define KEY0   3788876458
#define KEY1   2863311530
#define ITER   5

void uinform_refine(HashTable* El_Table, HashTable* NodeTable,
		TimeProps* timeprops_ptr, MatProps* matprops_ptr, int numprocs, int myid) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;
	int rescomp = 1;

	//for debugging perpose
	unsigned key[2] = { KEY0, KEY1 };
	double max=0;

#ifdef DEBUG
	double dummyv_star = 0.0;
	int adjflag = 1;
	tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr,
			dummyv_star, adjflag);
	int nonz1 = num_nonzero_elem(El_Table);

	cout << "number of elements before refinement  " << nonz1 << endl;

	int *dbgvec = new int[nonz1];
	int *pass = new int[nonz1];
	for (int i = 0; i < nonz1; i++) {
		dbgvec[i] = 0;
		pass[i] = 0;
	}
	if (checkElement(El_Table))
	exit(23);
#endif

//	if (checkElement(El_Table, &max, key))
//		cout << "here is the problem" << endl;

	htflush(El_Table, NodeTable, 1);
	move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() >= NOTRECADAPTED) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
				}
				currentPtr = currentPtr->next;
			}
		}
	}

	ElemPtrList RefinedList(num_nonzero_elem(El_Table));

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() == NOTRECADAPTED) {

					refinewrapper(El_Table, NodeTable, matprops_ptr, &RefinedList,
							Curr_El, rescomp);
				}

			}
		}
	}
//	if (checkElement(El_Table,&max, key))
//		cout << "here is the problem" << endl;

#ifdef DEBUG
	if (checkElement(El_Table))
	exit(24);
	cout << "number of elements -7   " << num_nonzero_elem(El_Table, -7) << endl
	<< "number of elements -6   " << num_nonzero_elem(El_Table, -6) << endl
	<< "number of elements  0   " << num_nonzero_elem(El_Table, 0) << endl
	<< "number of elements  1   " << num_nonzero_elem(El_Table, 1) << endl
	<< "number of elements  2   " << num_nonzero_elem(El_Table, 2) << endl
	<< "number of elements  3   " << num_nonzero_elem(El_Table, 3) << endl
	<< "number of elements  4   " << num_nonzero_elem(El_Table, 4) << endl
	<< "number of elements  5   " << num_nonzero_elem(El_Table, 5) << endl;
#endif

	bilinear_interp(El_Table);	//this function reconstruct linear interpolation

//	if (checkElement(El_Table, &max, key))
//		cout << "here is the problem" << endl;

	refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &RefinedList,
			timeprops_ptr);	//this function delete old father elements

//	if (checkElement(El_Table, &max, key))
//		cout << "here is the problem" << endl;
	RefinedList.trashlist();

	move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);

//	if (checkElement(El_Table, &max, key))
//		cout << "here is the problem" << endl;
	return;
}
