/*
 * dual_grid_routines.C
 *
 *  Created on: Mar 3, 2016
 *      Author: haghakha
 */


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"
#include "../header/exvar.h"

void update_dual_grid(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

//order of calling these functions are very important
//first we update prev_state_vatrs
//then compute k_act and fluxes

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;

//	move_dual_data(meshctx, propctx);

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->update_state(solrec, El_Table, iter);

				currentPtr = currentPtr->next;
			}
		}

	move_dual_data(meshctx, propctx);

//	int num_node_buckets = NodeTable->get_no_of_buckets();
//	buck = NodeTable->getbucketptr();
//	for (int i = 0; i < num_node_buckets; i++)
//		if (*(buck + i)) {
//			HashEntryPtr currentPtr = *(buck + i);
//			while (currentPtr) {
//				Node* Curr_Node = (Node*) (currentPtr->value);
//				Curr_Node->zero_flux();
//
//				currentPtr = currentPtr->next;
//			}
//		}

//this function computes slopes based on prev_state_vars and dh/dh_e where h_e is pile height in neighbor element
// we need this term to compute jacobian of elements
	slopes(El_Table, NodeTable, matprops_ptr, ADJOINT);

//	move_dual_data(meshctx, propctx);

	double tiny = GEOFLOW_TINY;
	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					gmfggetcoef_(Curr_El->get_prev_state_vars(), Curr_El->get_d_state_vars(),
					    (Curr_El->get_d_state_vars() + NUM_STATE_VARS), Curr_El->get_dx(),
					    &(matprops_ptr->bedfrict[Curr_El->get_material()]), &(matprops_ptr->intfrict),
					    (Curr_El->get_kactxy()), (Curr_El->get_kactxy() + 1), &tiny,
					    &(matprops_ptr->epsilon));

				currentPtr = currentPtr->next;
			}
		}

	// this is very necessary for ghost element jacobian computation
	move_dual_data(meshctx, propctx);

// this function computes fluxes based on prev_state_vars (we need for dual problem),
// and jacobian of fluxes and store the in elements
	calc_flux(meshctx, propctx);

}

void clear_empty_jacobians(SolRec* solrec, int iter) {

	HashEntryPtr *buck = solrec->getbucketptr();
	for (int i = 0; i < solrec->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Jacobian * jacobian = (Jacobian*) (currentPtr->value);

				currentPtr = currentPtr->next;

				jacobian->erase_solution(iter - 1);

				if (jacobian->is_container_empty()) {
					solrec->remove(jacobian->get_key());
					delete jacobian;
				}

			}
		}
}


