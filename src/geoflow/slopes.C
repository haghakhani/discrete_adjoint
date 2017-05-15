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
 * $Id: slopes.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/geoflow.h"

void slopes(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, run_mode dualcall) {
	int i;
	//-------------------go through all the elements of the subdomain------------------------
	//-------------------and   --------------------------

	HashEntryPtr* buck = El_Table->getbucketptr();
	/* mdj 2007-02 */
	HashEntryPtr currentPtr;
	Element* Curr_El;
	if (dualcall & FORWARD) {
		for (i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) { //if this element does not belong on this processor don't involve!!!
						Curr_El->get_slopes(El_Table, NodeTable, matprops_ptr->gamma);
					}
					currentPtr = currentPtr->next;
				}
			}
	} else if (dualcall & ADJOINT) {
		for (i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (DualElem*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) { //if this element does not belong on this processor don't involve!!!
//								if (Curr_El->get_ithelem() == 4
//								/**(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1
//								 /*&& timeprops_ptr->iter == ITER*/) {
//									int ddd, aa = 0;
//									int gg = ddd;
//								}
						Curr_El->get_slopes_prev(El_Table, NodeTable, matprops_ptr->gamma);
					}
					currentPtr = currentPtr->next;
				}
			}
	} else if (dualcall & ERROR) {

		for (i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (ErrorElem*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) { //if this element does not belong on this processor don't involve!!!

						Curr_El->get_slopes_prev(El_Table, NodeTable, matprops_ptr->gamma);
					}
					currentPtr = currentPtr->next;
				}
			}

	}

	return;
}
