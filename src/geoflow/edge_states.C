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
 * $Id: edge_states.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#define KEY0   2218506450
#define KEY1   2027059408
#define ITER   1

/*! calc_edge_states() cycles through the element Hashtable (listing of all 
 *  elements) and for each element (that has not been refined this iteration 
 *  and is not a ghost_element) calls Element member function 
 *  Element::calc_edge_states() (which calculates the Riemann fluxes across 
 *  the element's boundaries), and adds local boundary-element outflow to 
 *  GIS map's cummulative outflow (defined as the mass flow off of the
 *  GIS map).  Also, the elements are checked for multiple pile-height values 
 */

template<class T>
void calc_edge_states(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    TimeProps* timeprops_ptr, int myid, int* order_flag, double *outflow) {

	vector<T*> x_elem_list, y_elem_list;

	//-------------------go through all the elements of the subdomain and
	//-------------------find the edge states

	HashEntryPtr* buck = El_Table->getbucketptr();
	double localoutflow;
	*outflow = 0.0;
	/* mdj 2007-04 */
	double localoutflow_sum = 0.0;
	HashEntryPtr currentPtr;
	T* Curr_El;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (T*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					int aa=0,bb=1;
					if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1)
						bb=aa;

					Curr_El->calc_edge_states<T>(El_Table, NodeTable, &x_elem_list, &y_elem_list,
					    matprops_ptr, myid, timeprops_ptr->dtime, order_flag, &localoutflow);
					localoutflow_sum += localoutflow;
				}
				currentPtr = currentPtr->next;
			}
		}

	for (int i = 0; i < x_elem_list.size(); ++i)
		x_elem_list[i]->boundary_flux(El_Table, NodeTable, myid, 0);

	for (int i = 0; i < y_elem_list.size(); ++i)
		y_elem_list[i]->boundary_flux(El_Table, NodeTable, myid, 1);

	*outflow = localoutflow_sum;
	return;
}

void calc_flux(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	DualElem* Curr_El;

	vector<DualElem*> x_elem_list, y_elem_list;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					//if this element doesn't belong on this processor don't involve

//					if (Curr_El->get_ithelem() == 8255
//					//*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1
//					    /*&& timeprops_ptr->iter == ITER*/) {
//						int ddd, aa = 0;
//						int gg = ddd;
//					}

					Curr_El->calc_fluxes(El_Table, NodeTable, &x_elem_list, &y_elem_list, myid);

				}
				currentPtr = currentPtr->next;
			}
		}

	for (int i = 0; i < x_elem_list.size(); ++i)
		x_elem_list[i]->boundary_flux(El_Table, NodeTable, myid, 0);

	for (int i = 0; i < y_elem_list.size(); ++i)
		y_elem_list[i]->boundary_flux(El_Table, NodeTable, myid, 1);

	return;
}

void usefull_link() {

	HashTable* El_Table;
	MatProps* matprops_ptr;
	TimeProps* timeprops_ptr;
	int myid;
	double yek;

	calc_edge_states<Element>(El_Table, El_Table, matprops_ptr, timeprops_ptr, myid, &myid, &yek);
	calc_edge_states<ErrorElem>(El_Table, El_Table, matprops_ptr, timeprops_ptr, myid, &myid, &yek);

}

