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
 * $Id: delete_tab.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

template<typename T1, typename T2>
void copy_hashtables_objects(HashTable* El_Table, HashTable* cp_El_Table) {

	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				T1 *Curr = (T1*) (currentPtr->value);
				T2 *cp = new T2(Curr);
				cp_El_Table->add(cp->pass_key(), cp);
				currentPtr = currentPtr->next;

			}
		}
}

template void copy_hashtables_objects<Element, DualElem>(HashTable* El_Table, HashTable* cp_El_Table);

template<typename T1>
void delete_hashtables_objects(HashTable* El_Table) {
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				T1 *Curr = (T1*) (currentPtr->value);
				delete Curr;
				currentPtr = currentPtr->next;

			}
		}

	delete El_Table;
}

template void delete_hashtables_objects<DualElem>(HashTable* El_Table);

template void delete_hashtables_objects<Node>(HashTable* El_Table);

template void delete_hashtables_objects<Jacobian>(HashTable* El_Table);

template void delete_hashtables_objects<Element>(HashTable* El_Table);

template void delete_hashtables_objects<ErrorElem>(HashTable* El_Table);

void delete_data(SolRec* solrec, MeshCTX* meshctx, MeshCTX* error_meshctx,
		PropCTX* propctx) {

	HashTable* elem_table = meshctx->el_table;
	HashTable* node_table = meshctx->nd_table;

	if (elem_table)
		delete_hashtables_objects<DualElem>(elem_table);

	if (node_table)
		delete_hashtables_objects<Node>(node_table);

	if (solrec)
		delete_hashtables_objects<Jacobian>(solrec);

	delete propctx->discharge;
	delete propctx->fluxprops;
	delete propctx->mapnames;
	delete propctx->matprops;
	delete propctx->pileprops;
	delete propctx->statprops;
	delete propctx->timeprops;

#ifdef Error
	error.start();
	HashTable* err_elem_table = error_meshctx->el_table;
	HashTable* err_node_table = error_meshctx->nd_table;
	if (err_elem_table)
		delete_hashtables_objects<ErrorElem>(err_elem_table);

	if (err_node_table)
		delete_hashtables_objects<Node>(err_node_table);
	error.stop();
#endif
}

