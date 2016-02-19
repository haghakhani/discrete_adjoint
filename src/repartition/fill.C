#include "../header/hpfem.h"

// for BSFC repartitioning scheme
void Pack_element(void *sendel_in, ElemPack* elem, HashTable* HT_Node_Ptr, int destination_proc) {
	int j, i = 0;
	Element *sendel = (Element *) sendel_in;
	if ((sendel->key[0] == 0) && (sendel->key[1] == 0)) {
		printf("Element key={0,0} is being packed\n");
		i++;
	}

	Node* node;

	elem->myprocess = destination_proc;
	elem->generation = sendel->generation;
	elem->opposite_brother_flag = sendel->opposite_brother_flag;
	elem->material = sendel->material;

	for (i = 0; i < 8; i++) {
		elem->neigh_proc[i] = sendel->neigh_proc[i];
		elem->neigh_gen[i] = sendel->neigh_gen[i];
	}

	elem->refined = sendel->refined;
	elem->adapted = sendel->adapted;
	elem->which_son = sendel->which_son;
	elem->new_old = sendel->new_old;

	for (i = 0; i < KEYLENGTH; i++)
		elem->key[i] = sendel->key[i];

	for (i = 0; i < 4; i++)
		for (j = 0; j < KEYLENGTH; j++)
			elem->brothers[i][j] = sendel->brothers[i][j];

	for (i = 0; i < 8; i++)
		for (j = 0; j < KEYLENGTH; j++) {
			elem->node_key[i][j] = sendel->node_key[i][j];
			elem->neighbor[i][j] = sendel->neighbor[i][j];
			if (i < 4)
				elem->son[i][j] = sendel->son[i][j];
		}
	for (i = 0; i < EQUATIONS; i++) {
		elem->el_error[i] = sendel->el_error[i];

	}

	//and the node info:
	for (i = 0; i < 8; i++) {
		node = (Node*) HT_Node_Ptr->lookup(elem->node_key[i]);
		assert(node);
		elem->n_info[i] = node->info;
		for (j = 0; j < 2; j++)
			elem->n_coord[i][j] = node->coord[j];
		elem->node_elevation[i] = node->elevation;
	}

	node = (Node*) HT_Node_Ptr->lookup(elem->key);
	assert(node);
	elem->n_info[8] = node->info;
	for (j = 0; j < 2; j++)
		elem->n_coord[8][j] = node->coord[j];
	elem->node_elevation[8] = node->elevation;

	//geoflow stuff
	elem->positive_x_side = sendel->positive_x_side;
	elem->elevation = sendel->elevation;
	for (i = 0; i < DIMENSION; i++) {
		elem->dx[i] = sendel->dx[i];
		elem->eigenvxymax[i] = sendel->eigenvxymax[i];
		elem->kactxy[i] = sendel->kactxy[i];
		elem->zeta[i] = sendel->zeta[i];
		elem->curvature[i] = sendel->curvature[i];
		elem->d_gravity[i] = sendel->d_gravity[i];
	}
	for (i = 0; i < NUM_STATE_VARS; i++) {
		elem->state_vars[i] = sendel->state_vars[i];
		elem->prev_state_vars[i] = sendel->prev_state_vars[i];
		elem->Influx[i] = sendel->Influx[i];
	}
	for (i = 0; i < 3; i++)
		elem->gravity[i] = sendel->gravity[i];

	for (i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
		elem->d_state_vars[i] = sendel->d_state_vars[i];

	elem->shortspeed = sendel->shortspeed;
	elem->lb_weight = sendel->lb_weight;
	elem->elm_loc[0] = sendel->elm_loc[0];
	elem->elm_loc[1] = sendel->elm_loc[1];

	elem->iwetnode = sendel->iwetnode;
	elem->Awet = sendel->Awet;
	elem->Swet = sendel->Swet;
	elem->drypoint[0] = sendel->drypoint[0];
	elem->drypoint[1] = sendel->drypoint[1];
}

//void Pack_dual_element(void *sendel_in, ElemPack* elem, HashTable* HT_Node_Ptr,
//		int destination_proc) {
//	int j, i = 0;
//	DualElem *sendel = (DualElem *) sendel_in;
//
//	Node* node;
//
//	elem->myprocess = destination_proc;
//	elem->generation = sendel->get_gen();
//	elem->opposite_brother_flag = sendel->get_opposite_brother_flag();
//	elem->material = sendel->get_material();
//
//	for (i = 0; i < 8; i++) {
//		elem->neigh_proc[i] = *(sendel->get_neigh_proc() + i);
//		elem->neigh_gen[i] = *(sendel->get_neigh_gen() + i);
//	}
//
//	elem->refined = sendel->get_refined_flag();
//	elem->adapted = sendel->get_adapted_flag();
//	elem->which_son = sendel->get_which_son();
//	elem->new_old = sendel->get_new_old();
//
//	for (i = 0; i < KEYLENGTH; i++)
//		elem->key[i] = *(sendel->pass_key()+i);
//
//	for (i = 0; i < 4; i++)
//		for (j = 0; j < KEYLENGTH; j++)
//			elem->brothers[i][j] = *(sendel->get_brothers() + i * KEYLENGTH + j);
//
//	for (i = 0; i < 8; i++)
//		for (j = 0; j < KEYLENGTH; j++) {
//			elem->node_key[i][j] = *(sendel->getNode() + i * KEYLENGTH + j);
//			elem->neighbor[i][j] = *(sendel->get_neighbors() + i * KEYLENGTH + j);
//			if (i < 4)
//				elem->son[i][j] = *(sendel->getson() + i * KEYLENGTH + j);
//
//		}
//	for (i = 0; i < EQUATIONS; i++) {
//		elem->el_error[i] = *(sendel->get_el_err()+i);
//
//	}
//
//	//and the node info:
//	for (i = 0; i < 8; i++) {
//		node = (Node*) HT_Node_Ptr->lookup(elem->node_key[i]);
//		assert(node);
//		elem->n_info[i] = node->info;
//		for (j = 0; j < 2; j++)
//			elem->n_coord[i][j] = node->coord[j];
//		elem->node_elevation[i] = node->elevation;
//	}
//
//	node = (Node*) HT_Node_Ptr->lookup(elem->key);
//	assert(node);
//	elem->n_info[8] = node->info;
//	for (j = 0; j < 2; j++)
//		elem->n_coord[8][j] = node->coord[j];
//	elem->node_elevation[8] = node->elevation;
//
//	//geoflow stuff
//	elem->positive_x_side = sendel->positive_x_side;
//	elem->elevation = sendel->elevation;
//	for (i = 0; i < DIMENSION; i++) {
//		elem->dx[i] = sendel->dx[i];
//		elem->eigenvxymax[i] = sendel->eigenvxymax[i];
//		elem->kactxy[i] = sendel->kactxy[i];
//		elem->zeta[i] = sendel->zeta[i];
//		elem->curvature[i] = sendel->curvature[i];
//		elem->d_gravity[i] = sendel->d_gravity[i];
//	}
//	for (i = 0; i < NUM_STATE_VARS; i++) {
//		elem->state_vars[i] = sendel->state_vars[i];
//		elem->prev_state_vars[i] = sendel->prev_state_vars[i];
//		elem->Influx[i] = sendel->Influx[i];
//	}
//	for (i = 0; i < 3; i++)
//		elem->gravity[i] = sendel->gravity[i];
//
//	for (i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
//		elem->d_state_vars[i] = sendel->d_state_vars[i];
//
//	elem->shortspeed = sendel->shortspeed;
//	elem->lb_weight = sendel->lb_weight;
//	elem->elm_loc[0] = sendel->elm_loc[0];
//	elem->elm_loc[1] = sendel->elm_loc[1];
//
//	elem->iwetnode = sendel->iwetnode;
//	elem->Awet = sendel->Awet;
//	elem->Swet = sendel->Swet;
//	elem->drypoint[0] = sendel->drypoint[0];
//	elem->drypoint[1] = sendel->drypoint[1];
//}

void Pack_neighbor(int target_proc, ELinkPtr* EL_head, int* counter, NePtr* packed_neighbor_info) {

	//creates a packing for the changed neighbor info
	*counter = 0;
	int counter2 = 0;

	ELinkPtr EL_temp = *EL_head;
	NeighborPack* pack_try;

	if (*EL_head) {
		while (EL_temp) {
			if (EL_temp->target_proc == target_proc)
				(*counter)++;
			EL_temp = EL_temp->next;

		}

		if (*counter) {
			pack_try = new NeighborPack[*counter];
			EL_temp = *EL_head;

			while (counter2 < *counter && EL_temp) {

				if (EL_temp->target_proc == target_proc) {
					pack_try[counter2].target_proc = target_proc;
					pack_try[counter2].new_proc = EL_temp->new_proc;
					for (int i = 0; i < KEYLENGTH; i++) {
						pack_try[counter2].elkey[i] = EL_temp->elkey[i];
						pack_try[counter2].targetkey[i] = EL_temp->targetkey[i];
					}
					counter2++;
				}

				EL_temp = EL_temp->next;
			}
		}
		*packed_neighbor_info = pack_try;
	}
}
