/*
 * error_grid_routines.C
 *
 *  Created on: Mar 3, 2016
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

void update_bilinear_error_grid(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;
	int yek = 1;
	double outflow = 0.;
	double tiny = GEOFLOW_TINY;

	slopes(El_Table, NodeTable, matprops_ptr, 2);

	move_err_data(meshctx, propctx);

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				ErrorElem * Curr_El = (ErrorElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					gmfggetcoef_(Curr_El->get_bilin_prev_state(), Curr_El->get_d_state_vars(),
					    (Curr_El->get_d_state_vars() + NUM_STATE_VARS), Curr_El->get_dx(),
					    &(matprops_ptr->bedfrict[Curr_El->get_material()]), &(matprops_ptr->intfrict),
					    (Curr_El->get_kactxy()), (Curr_El->get_kactxy() + 1), &tiny,
					    &(matprops_ptr->epsilon));

				currentPtr = currentPtr->next;
			}
		}

	move_err_data(meshctx, propctx);

	calc_edge_states<ErrorElem>(El_Table, NodeTable, matprops_ptr, timeprops_ptr, myid, &yek,
	    &outflow);

}

void update_error_grid(SolRec* solrec, MeshCTX* cp_meshctx, PropCTX* propctx) {

	HashTable* cp_El_Table = cp_meshctx->el_table;
	HashTable* cp_NodeTable = cp_meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;
	int iter = propctx->timeprops->iter;

	HashEntryPtr currentPtr;
	ErrorElem *Curr_El;

	HashEntryPtr *buck = cp_El_Table->getbucketptr();
	for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (ErrorElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->error_update_state(solrec, iter);

				currentPtr = currentPtr->next;
			}
		}
}

void set_new_fathers(HashTable* El_Table, vector<pair<unsigned, unsigned> >& new_father) {

	HashEntryPtr currentPtr;
	Element *Curr_El;

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() == NOTRECADAPTED) {

					for (int j = 0; j < new_father.size(); ++j) {
						unsigned key[] = { new_father[j].first, new_father[j].second };
						if (compare_key(Curr_El->pass_key(), key))
							Curr_El->put_adapted_flag(NEWFATHER);
					}
				}
				currentPtr = currentPtr->next;
			}
		}
}

void make_dual_err_link(HashTable *Dual_El_Tab, HashTable *Err_El_Tab) {

	HashEntryPtr *buck = Err_El_Tab->getbucketptr();
	HashEntryPtr currentPtr;

	for (int i = 0; i < Err_El_Tab->get_no_of_buckets(); ++i)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				ErrorElem *son[4];

				son[0] = (ErrorElem*) (currentPtr->value);

				son[0]->calc_which_son();

				if (son[0]->get_adapted_flag() > 0 && son[0]->get_which_son() == 0) {

					unsigned* father_key = son[0]->getfather();

					DualElem* father = (DualElem*) Dual_El_Tab->lookup(father_key);

					vector<ErrorElem*>& mysons = father->get_son_addresses();

					for (int j = 1; j < 4; ++j)
						son[j] = (ErrorElem*) Err_El_Tab->lookup(son[0]->get_brothers() + j * KEYLENGTH);

					for (int j = 0; j < 4; ++j) {
						mysons.push_back(son[j]);
						son[j]->put_father_address(father);
					}
				}
				currentPtr = currentPtr->next;
			}
		}
}

void send_from_dual_to_error(HashTable *Dual_El_Tab, HashTable *Err_El_Tab, int last) {

	HashEntryPtr *buck = Dual_El_Tab->getbucketptr();
	HashEntryPtr currentPtr;

	for (int i = 0; i < Dual_El_Tab->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				DualElem* dual_el = (DualElem*) (currentPtr->value);
//				double max_x = 0., min_x = 1.e8, max_y = 0., min_y = 1.e8;

				if (dual_el->get_adapted_flag() > 0) {

					vector<ErrorElem*>& error_el = dual_el->get_son_addresses();

					for (int j = 0; j < error_el.size(); ++j) {
						ErrorElem* err_El = error_el[j];

						for (int k = 0; k < NUM_STATE_VARS; ++k) {
//							if (!last)
							*(err_El->get_prev_state_vars() + k) = *(dual_el->get_prev_state_vars() + k);
							*(err_El->get_state_vars() + k) = *(dual_el->get_state_vars() + k);
							*(err_El->get_adjoint() + k) = *(dual_el->get_adjoint() + k);

						}

//						if (*(err_El->get_coord()) + .5 * *(err_El->get_dx()) > max_x)
//							max_x = *(err_El->get_coord()) + .5 * *(err_El->get_dx());
//
//						if (*(err_El->get_coord()) - .5 * *(err_El->get_dx()) < min_x)
//							min_x = *(err_El->get_coord()) - .5 * *(err_El->get_dx());
//
//						if (*(err_El->get_coord() + 1) + .5 * *(err_El->get_dx() + 1) > max_y)
//							max_y = *(err_El->get_coord() + 1) + .5 * *(err_El->get_dx() + 1);
//
//						if (*(err_El->get_coord() + 1) - .5 * *(err_El->get_dx() + 1) < min_y)
//							min_y = *(err_El->get_coord() + 1) - .5 * *(err_El->get_dx() + 1);

					}

//					assert(fabs(max_x - (*(dual_el->get_coord()) + .5 * *(dual_el->get_dx()))) < 1e-8);
//					assert(fabs(min_x - (*(dual_el->get_coord()) - .5 * *(dual_el->get_dx()))) < 1e-8);
//					assert(
//					    fabs(max_y - (*(dual_el->get_coord() + 1) + .5 * *(dual_el->get_dx() + 1))) < 1e-8);
//					assert(
//					    fabs(min_y - (*(dual_el->get_coord() + 1) - .5 * *(dual_el->get_dx() + 1))) < 1e-8);

				}
				currentPtr = currentPtr->next;
			}

		}
}

