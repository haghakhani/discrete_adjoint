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
 * $Id: element2.C 143 2007-06-25 17:58:08Z dkumar $ 
 */
//#define DEBUG_SAVE_ELEM
//#define SHORTSPEED
#define DISABLE_DRY_FLUX_ZEROING   //this disables the zeroing of fluxes through dry sides facet of thin layer control

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

//stopping criteria define statements have been moved to geoflow.h

//#define PRINT_GIS_ERRORS

/*  original element   */
Element::Element(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[], int mat,
    int* elm_loc_in, double pile_height, int myid, unsigned* opposite_brother) {

	adapted = NOTRECADAPTED;

	for (int i = 0; i < NUM_STATE_VARS; i++)
		prev_state_vars[i] = 0.;

	for (int i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
		d_state_vars[i] = 0.;

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	for (int i = 0; i < 4; i++)
		son[i][0] = son[i][1] = brothers[i][0] = brothers[i][1] = NULL;
	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	new_old = OLD;
	generation = 0; //--first generation
	material = mat;
	for (int i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (int i = 0; i < KEYLENGTH; i++) {
		father[i] = NULL;
		key[i] = nodekeys[8][i]; //--using bubble key to represent the element
	}

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < KEYLENGTH; j++)
			node_key[i][j] = nodekeys[i][j];

	for (int i = 0; i < 4; i++) {
		neigh_proc[i] = n_pro[i];
		neigh_proc[i + 4] = -2; //-- -2 means regular element
		if (neigh_proc[i] != -1)
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = neigh[i][j];
		else
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = NULL;
	}

	for (int i = 0; i < 8; i++)
		neigh_gen[i] = 0;

	refined = 0;

	myprocess = myid;
	elm_loc[0] = elm_loc_in[0];
	elm_loc[1] = elm_loc_in[1];
	calc_which_son();
	for (int i = 0; i < KEYLENGTH; i++) {
		brothers[which_son][i] = key[i];
		brothers[(which_son + 2) % 4][i] = opposite_brother[i];
	}

	switch (which_son) {
		case 0:
			for (int i = 0; i < KEYLENGTH; i++) {
				brothers[1][i] = neighbor[1][i];
				brothers[3][i] = neighbor[2][i];
			}
			break;
		case 1:
			for (int i = 0; i < KEYLENGTH; i++) {
				brothers[0][i] = neighbor[3][i];
				brothers[2][i] = neighbor[2][i];
			}
			break;
		case 2:
			for (int i = 0; i < KEYLENGTH; i++) {
				brothers[1][i] = neighbor[0][i];
				brothers[3][i] = neighbor[3][i];
			}
			break;
		case 3:
			for (int i = 0; i < KEYLENGTH; i++) {
				brothers[0][i] = neighbor[0][i];
				brothers[2][i] = neighbor[1][i];
			}
			break;
	}
	opposite_brother_flag = 1;

	new_old = OLD;
	state_vars[0] = pile_height;
	state_vars[2] = state_vars[1] = 0.;
	if (state_vars[0] != 0)
		state_vars[1] = state_vars[2] = 0.0001;

	prev_state_vars[0] = pile_height;
	prev_state_vars[1] = prev_state_vars[2] = 0.;
	if (prev_state_vars[0] != 0)
		prev_state_vars[1] = prev_state_vars[2] = 0.0001;

	for (int i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
		d_state_vars[i] = 0;
	for (int i = 0; i < NUM_STATE_VARS; i++)
		Influx[i] = 0.0;

	// initialize kactxy
	kactxy[0] = kactxy[1] = 0.;
	tan_bed_frict = 0.;
}

//used for refinement
Element::Element(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[], int gen,
    int elm_loc_in[], int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
    HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr) {

	adapted = NEWSON;

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		prev_state_vars[i] = 0.;
		Influx[i] = 0.;
	}

	for (int i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
		d_state_vars[i] = 0.;

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;
	for (int i = 0; i < 4; i++)
		son[i][0] = son[i][1] = brothers[i][0] = brothers[i][1] = NULL;
	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	myprocess = myid;
	generation = gen; //--first generation
	opposite_brother_flag = 1;
	material = mat;
	for (int i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (int i = 0; i < KEYLENGTH; i++) {
		father[i] = fthTemp->key[i];
		key[i] = nodekeys[8][i]; //--using buble key to represent the element
	}

	elm_loc[0] = elm_loc_in[0];
	elm_loc[1] = elm_loc_in[1];

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < KEYLENGTH; j++)
			node_key[i][j] = nodekeys[i][j];

	for (int i = 0; i < 4; i++) {
		neigh_proc[i] = n_pro[i];
		neigh_proc[i + 4] = -2; //-- -2 means regular element
		if (neigh_proc[i] != -1)
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = neigh[i][j];
		else
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = NULL;

		neigh_gen[i] = neigh_gen[i + 4] = gen_neigh[i];
	}

	refined = 0;

	new_old = NEW;
	//geoflow stuff
	dx[0] = .5 * fthTemp->dx[0];  //assume constant refinement
	dx[1] = .5 * fthTemp->dx[1];

	double myfractionoffather=1.0;

	find_positive_x_side(NodeTable);
	calc_topo_data(matprops_ptr);
	calc_gravity_vector(matprops_ptr);
	calc_d_gravity(El_Table);

	coord[0] = coord_in[0];
	coord[1] = coord_in[1];

	calc_which_son();

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		state_vars[i] = fthTemp->state_vars[i] * myfractionoffather;
		prev_state_vars[i] = fthTemp->prev_state_vars[i] * myfractionoffather;
	}
}
/*********************************
 making a father element from its sons
 *****************************************/
Element::Element(Element* sons[], HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr) {
	adapted = NEWFATHER;

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		prev_state_vars[i] = 0.;
		Influx[i] = 0.;
		Influx[i] = 0.;
		d_state_vars[i] = d_state_vars[NUM_STATE_VARS + i] = 0.;
	}

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	int i, j, ikey, ison, isonneigh, ineigh;

	for (ikey = 0; ikey < KEYLENGTH; ikey++)
		key[ikey] = *(sons[2]->getNode() + ikey);

	for (ison = 0; ison < 4; ison++) {
		sons[ison]->put_adapted_flag(OLDSON);
		for (ikey = 0; ikey < KEYLENGTH; ikey++) {
			son[ison][ikey] = *(sons[ison]->pass_key() + ikey);
			sons[ison]->put_father(key);
		}
	}

	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	new_old = NEW;
	unsigned* son_nodes[4];
	opposite_brother_flag = 0;
	for (i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (ison = 0; ison < 4; ison++)
		son_nodes[ison] = sons[ison]->getNode();

	for (ikey = 0; ikey < KEYLENGTH; ikey++) {
		father[ikey] = NULL;
		node_key[0][ikey] = son_nodes[0][ikey];
		node_key[1][ikey] = son_nodes[1][KEYLENGTH + ikey];
		node_key[2][ikey] = son_nodes[2][2 * KEYLENGTH + ikey];
		node_key[3][ikey] = son_nodes[3][3 * KEYLENGTH + ikey];
		node_key[4][ikey] = son_nodes[0][KEYLENGTH + ikey];
		node_key[5][ikey] = son_nodes[1][2 * KEYLENGTH + ikey];
		node_key[6][ikey] = son_nodes[2][3 * KEYLENGTH + ikey];
		node_key[7][ikey] = son_nodes[3][ikey];

		elm_loc[ikey] = (*(sons[0]->get_elm_loc() + ikey)) / 2;
	}
	myprocess = sons[0]->get_myprocess();
	generation = sons[0]->get_gen() - 1;
	material = sons[0]->get_material();

	calc_which_son();

	refined = 1; // not an active element yet!!!

	// neighbor information
	for (ison = 0; ison < 4; ison++) {
		isonneigh = ison;
		ineigh = isonneigh;
		neigh_gen[ineigh] = *(sons[ison]->get_neigh_gen() + isonneigh);
		for (ikey = 0; ikey < KEYLENGTH; ikey++)
			neighbor[ineigh][ikey] = *(sons[ison]->get_neighbors() + isonneigh * KEYLENGTH + ikey);
		neigh_proc[ineigh] = *(sons[ison]->get_neigh_proc() + isonneigh);

		isonneigh = (ison + 3) % 4;
		ineigh = isonneigh + 4;
		neigh_gen[ineigh] = *(sons[ison]->get_neigh_gen() + isonneigh);
		for (ikey = 0; ikey < KEYLENGTH; ikey++)
			neighbor[ineigh][ikey] = *(sons[ison]->get_neighbors() + isonneigh * KEYLENGTH + ikey);
		if ((*(sons[ison]->get_neigh_gen() + isonneigh) == generation)
		    || (*(sons[ison]->get_neigh_proc() + isonneigh) == -1))
			neigh_proc[ineigh] = -2;
		else
			neigh_proc[ineigh] = *(sons[ison]->get_neigh_proc() + isonneigh);
	}

	/* brother information -- requires that at least one of this
	 element's neighboring brothers is on this process in
	 order to get information on the brother that is not a neighbor */
	Element* EmTemp;
	switch (which_son) {
		case 0:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[0][i] = key[i];
			if (neigh_proc[1] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = 0;
			} else if (neigh_gen[1] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = neighbor[1][i];
			} else if (neigh_gen[1] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[1]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = bro_key[i];
			} else if (neigh_gen[1] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[2] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = 0;
			} else if (neigh_gen[2] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = neighbor[2][i];
			} else if (neigh_gen[2] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[2]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = bro_key[i];
			} else if (neigh_gen[2] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 1:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[1][i] = key[i];
			if (neigh_proc[3] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = 0;
			} else if (neigh_gen[3] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = neighbor[3][i];
			} else if (neigh_gen[3] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[3]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = bro_key[i];
			} else if (neigh_gen[3] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[2] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = 0;
			} else if (neigh_gen[2] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = neighbor[2][i];
			} else if (neigh_gen[2] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[2]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = bro_key[i];
			} else if (neigh_gen[2] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 2:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[2][i] = key[i];
			if (neigh_proc[0] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = 0;
			} else if (neigh_gen[0] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = neighbor[0][i];
			} else if (neigh_gen[0] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[0]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = bro_key[i];
			} else if (neigh_gen[0] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[3] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = 0;
			} else if (neigh_gen[3] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = neighbor[3][i];
			} else if (neigh_gen[3] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[3]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = bro_key[i];
			} else if (neigh_gen[3] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 3:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[3][i] = key[i];
			if (neigh_proc[0] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = 0;
			} else if (neigh_gen[0] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = neighbor[0][i];
			} else if (neigh_gen[0] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[0]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = bro_key[i];
			} else if (neigh_gen[0] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[1] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = 0;
			} else if (neigh_gen[1] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = neighbor[1][i];
			} else if (neigh_gen[1] == generation + 1) {
				EmTemp = (Element*) El_Table->lookup(neighbor[1]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = bro_key[i];
			} else if (neigh_gen[1] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
	}

	find_positive_x_side(NodeTable);  //also inserts the coordinates
	calculate_dx(NodeTable);
	find_opposite_brother(El_Table);
	calc_topo_data(matprops_ptr);
	calc_gravity_vector(matprops_ptr);
	calc_d_gravity(El_Table);

	for (i = 0; i < NUM_STATE_VARS; i++) {
		state_vars[i] = 0.;
		prev_state_vars[i] = 0.;
		for (j = 0; j < 4; j++) {
			state_vars[i] += *(sons[j]->get_state_vars() + i) * 0.25;
			prev_state_vars[i] += *(sons[j]->get_prev_state_vars() + i) * 0.25;
		}
	}
}

unsigned* Element::getfather() {
	switch (which_son) {
		case 0:
			return node_key[2];
			break;
		case 1:
			return node_key[3];
			break;
		case 2:
			return node_key[0];
			break;
		case 3:
			return node_key[1];
			break;
	}
	printf("my key is %u %u in getfather on proc %d\n", key[0], key[1], myprocess);
	assert(0); // 0 <= which_son <= 3 !!!
}

int Element::which_neighbor(unsigned* FindNeigh) {
	int i;
	for (i = 0; i < 8; i++)
		if (compare_key(neighbor[i], FindNeigh) && (neigh_proc[i] >= 0))
			return i;

	assert(i < 8);

	return i;
}

void Element::change_neighbor(unsigned* newneighbs, int which_side, int proc, int reg) {
	int j;
	switch (reg) {
		case 1:
			j = 0;
			break;
		case 3:
			assert(which_side < 4);
			for (j = 0; j < KEYLENGTH; j++) {
				neighbor[which_side][j] = *(newneighbs + j);
				neighbor[which_side + 4][j] = *(newneighbs + KEYLENGTH + j);
			}
			neigh_proc[which_side + 4] = proc; //assuming no element movement
			neigh_gen[which_side] = neigh_gen[which_side + 4] = neigh_gen[which_side] + 1;
			break;
		case 4:
			j = 0;
			break;
		case 2:
			j = 0;
			break;
		case 5:
			for (j = 0; j < KEYLENGTH; j++)
				neighbor[which_side][j] = *(newneighbs + j);
			neigh_gen[which_side] = neigh_gen[which_side] + 1;
			break;

		case 6:
			for (j = 0; j < KEYLENGTH; j++)
				neighbor[which_side][j] = neighbor[which_side + 4][j] = *(newneighbs + j);
			neigh_gen[which_side] = neigh_gen[which_side + 4] = neigh_gen[which_side] + 1;
			break;

			/*Andrew's section called from update_interproc*/
		case 10: //the refined element and old neighbor have the same gen.
			assert(which_side < 4);
			for (j = 0; j < KEYLENGTH; j++) {
				neighbor[which_side][j] = *(newneighbs + j);
				neighbor[which_side + 4][j] = *(newneighbs + KEYLENGTH + j);
			}
			neigh_proc[which_side + 4] = proc;

			neigh_gen[which_side] = neigh_gen[which_side + 4] = neigh_gen[which_side] + 1;
			break;

		case 11:
			for (j = 0; j < KEYLENGTH; j++)
				neighbor[which_side][j] = neighbor[which_side + 4][j] = *(newneighbs + j);
			neigh_gen[which_side] = neigh_gen[which_side + 4] = neigh_gen[which_side] + 1;
			break;

	}
}

/* routine also puts in the coords of the center node in the elm */
void Element::find_positive_x_side(HashTable* nodetable) {
	int i, j, side;
	double xmax;
	Node* nodeptr;

	nodeptr = (Node*) (nodetable->lookup(key));
	xmax = nodeptr->coord[0];
	coord[0] = xmax;
	coord[1] = nodeptr->coord[1];

	for (i = 4; i < 8; i++) {
		nodeptr = (Node*) (nodetable->lookup(node_key[i]));
		double xcoord = nodeptr->coord[0];
		if (xcoord > xmax) {
			xmax = xcoord;
			side = i - 4;
		}
	}

	positive_x_side = side;

	return;
}

void Element::get_slopes_prev(HashTable* El_Table, HashTable* NodeTable, double gamma) {

	int j = 0, bc = 0;
	/* check to see if this is a boundary */
	while (j < 4 && bc == 0) {
		if (neigh_proc[j] == INIT)
			bc = 1;
		j++;
	}
	if (bc == 1) {
		for (j = 0; j < NUM_STATE_VARS * DIMENSION; j++)
			d_state_vars[j] = 0;
		return;
	}

	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}
	/* x direction */
	Element *ep = (Element*) (El_Table->lookup(&neighbor[xp][0]));
	Element *em = (Element*) (El_Table->lookup(&neighbor[xm][0]));
	Element *ep2 = NULL;
	Element *em2 = NULL;
	//check if element has 2 neighbors on either side
	Node* ndtemp = (Node*) NodeTable->lookup(&node_key[xp + 4][0]);
	if (ndtemp->info == S_C_CON) {
		ep2 = (Element*) (El_Table->lookup(&neighbor[xp + 4][0]));
		assert(neigh_proc[xp + 4] >= 0 && ep2);
	}
	ndtemp = (Node*) NodeTable->lookup(&node_key[xm + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (Element*) (El_Table->lookup(&neighbor[xm + 4][0]));
		assert(neigh_proc[xm + 4] >= 0 && em2);
	}

	double dp, dm, dc, dxp, dxm;
	dxp = ep->coord[0] - coord[0];
	dxm = coord[0] - em->coord[0];
	for (j = 0; j < NUM_STATE_VARS; j++) {
		dp = (ep->prev_state_vars[j] - prev_state_vars[j]) / dxp;
		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->prev_state_vars[j] - prev_state_vars[j]) / dxp);
		dm = (prev_state_vars[j] - em->prev_state_vars[j]) / dxm;
		if (em2 != NULL)
			dm = .5 * (dm + (prev_state_vars[j] - em2->prev_state_vars[j]) / dxm);

		dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
		//do slope limiting
		d_state_vars[j] = .5 * (c_sgn(dp) + c_sgn(dm))
		    * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
	}

	/* y direction */
	ep = (Element*) (El_Table->lookup(&neighbor[yp][0]));
	em = (Element*) (El_Table->lookup(&neighbor[ym][0]));
	ep2 = NULL;
	em2 = NULL;
	//check if element has 2 neighbors on either side
	ndtemp = (Node*) NodeTable->lookup(&node_key[yp + 4][0]);
	if (ndtemp->info == S_C_CON) {
		ep2 = (Element*) (El_Table->lookup(&neighbor[yp + 4][0]));
		assert(neigh_proc[yp + 4] >= 0 && ep2);
	}
	ndtemp = (Node*) NodeTable->lookup(&node_key[ym + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (Element*) (El_Table->lookup(&neighbor[ym + 4][0]));
		if (!(neigh_proc[ym + 4] >= 0 && em2)) {
			printf("ym=%d neigh_proc[ym+4]=%d em2=%d\n", ym, neigh_proc[ym + 4], em2);
		}
		assert(neigh_proc[ym + 4] >= 0 && em2);
	}

	dxp = ep->coord[1] - coord[1];
	dxm = coord[1] - em->coord[1];
	for (j = 0; j < NUM_STATE_VARS; j++) {
		dp = (ep->prev_state_vars[j] - prev_state_vars[j]) / dxp;
		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->prev_state_vars[j] - prev_state_vars[j]) / dxp);
		dm = (prev_state_vars[j] - em->prev_state_vars[j]) / dxm;
		if (em2 != NULL)
			dm = .5 * (dm + (prev_state_vars[j] - em2->prev_state_vars[j]) / dxm);

		dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
		//do slope limiting
		d_state_vars[j + NUM_STATE_VARS] = .5 * (c_sgn(dp) + c_sgn(dm))
		    * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
	}

	return;
}

void Element::get_slopes(HashTable* El_Table, HashTable* NodeTable, double gamma) {
	int j = 0, bc = 0;
	/* check to see if this is a boundary */
	while (j < 4 && bc == 0) {
		if (neigh_proc[j] == INIT)
			bc = 1;
		j++;
	}
	if (bc == 1) {
		for (j = 0; j < NUM_STATE_VARS * DIMENSION; j++)
			d_state_vars[j] = 0;
		return;
	}

	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}
	/* x direction */
	Element *ep = (Element*) (El_Table->lookup(&neighbor[xp][0]));
	Element *em = (Element*) (El_Table->lookup(&neighbor[xm][0]));
	Element *ep2 = NULL;
	Element *em2 = NULL;
//check if element has 2 neighbors on either side
	Node* ndtemp = (Node*) NodeTable->lookup(&node_key[xp + 4][0]);
	if (ndtemp->info == S_C_CON) {
		ep2 = (Element*) (El_Table->lookup(&neighbor[xp + 4][0]));
		assert(neigh_proc[xp + 4] >= 0 && ep2);
	}
	ndtemp = (Node*) NodeTable->lookup(&node_key[xm + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (Element*) (El_Table->lookup(&neighbor[xm + 4][0]));
		assert(neigh_proc[xm + 4] >= 0 && em2);
	}

	double dp, dm, dc, dxp, dxm;
	dxp = ep->coord[0] - coord[0];
	dxm = coord[0] - em->coord[0];
	for (j = 0; j < NUM_STATE_VARS; j++) {
		dp = (ep->state_vars[j] - state_vars[j]) / dxp;
		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->state_vars[j] - state_vars[j]) / dxp);
		dm = (state_vars[j] - em->state_vars[j]) / dxm;
		if (em2 != NULL)
			dm = .5 * (dm + (state_vars[j] - em2->state_vars[j]) / dxm);

		dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
		//do slope limiting
		d_state_vars[j] = .5 * (c_sgn(dp) + c_sgn(dm))
		    * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
	}

	/* y direction */
	ep = (Element*) (El_Table->lookup(&neighbor[yp][0]));
	em = (Element*) (El_Table->lookup(&neighbor[ym][0]));
	ep2 = NULL;
	em2 = NULL;
//check if element has 2 neighbors on either side
	ndtemp = (Node*) NodeTable->lookup(&node_key[yp + 4][0]);
	if (ndtemp->info == S_C_CON) {
		ep2 = (Element*) (El_Table->lookup(&neighbor[yp + 4][0]));
		assert(neigh_proc[yp + 4] >= 0 && ep2);
	}
	ndtemp = (Node*) NodeTable->lookup(&node_key[ym + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (Element*) (El_Table->lookup(&neighbor[ym + 4][0]));
		if (!(neigh_proc[ym + 4] >= 0 && em2)) {
			printf("ym=%d neigh_proc[ym+4]=%d em2=%d\n", ym, neigh_proc[ym + 4], em2);
		}
		assert(neigh_proc[ym + 4] >= 0 && em2);
	}

	dxp = ep->coord[1] - coord[1];
	dxm = coord[1] - em->coord[1];
	for (j = 0; j < NUM_STATE_VARS; j++) {
		dp = (ep->state_vars[j] - state_vars[j]) / dxp;
		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->state_vars[j] - state_vars[j]) / dxp);
		dm = (state_vars[j] - em->state_vars[j]) / dxm;
		if (em2 != NULL)
			dm = .5 * (dm + (state_vars[j] - em2->state_vars[j]) / dxm);

		dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
		//do slope limiting
		d_state_vars[j + NUM_STATE_VARS] = .5 * (c_sgn(dp) + c_sgn(dm))
		    * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
	}

	return;
}

void Element::set_mins(HashTable* NodeTable) {
	int i, j;
	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}

	Node *np, *nm;

	np = (Node*) NodeTable->lookup(node_key[xp + 4]);
	nm = (Node*) NodeTable->lookup(node_key[xm + 4]);

	dx[0] = (np->coord[0] - nm->coord[0]) /*(zeta[0]*zeta[0]+1)*/;

	if (dx[0] == 0)
		cout<<"ERROR in dx[0] \n";
//		printf("np %p,nm %p,dx, np_coord, nm_coord %e %e\n", np, nm, np->coord[0], nm->coord[0]);

	np = (Node*) NodeTable->lookup(node_key[yp + 4]);
	nm = (Node*) NodeTable->lookup(node_key[ym + 4]);

	dx[1] = (np->coord[1] - nm->coord[1]) /*(zeta[1]*zeta[1]+1)*/;

	min_gen=generation;
	min_dx[0]=dx[0];
	min_dx[1]=dx[1];

	if (dx[1] == 0)
		cout<<"ERROR in dx[1] \n";
//		printf("dy, np_coord, nm_coord %e %e\n", np->coord[1], nm->coord[1]);

	return;
}


void Element::calculate_dx(HashTable* NodeTable) {
	int i, j;
	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}

	Node *np, *nm;

	np = (Node*) NodeTable->lookup(node_key[xp + 4]);
	nm = (Node*) NodeTable->lookup(node_key[xm + 4]);

	dx[0] = (np->coord[0] - nm->coord[0]) /*(zeta[0]*zeta[0]+1)*/;
//	double dif_gen = generation - min_gen;
//	dx[0] = min_dx[0] * pow(.5, dif_gen);
//	dx[1] = min_dx[1] * pow(.5, dif_gen);

	if (dx[0] == 0)
		cout<<"ERROR in dx[0] \n";
//		printf("np %p,nm %p,dx, np_coord, nm_coord %e %e\n", np, nm, np->coord[0], nm->coord[0]);

	np = (Node*) NodeTable->lookup(node_key[yp + 4]);
	nm = (Node*) NodeTable->lookup(node_key[ym + 4]);

	dx[1] = (np->coord[1] - nm->coord[1]) /*(zeta[1]*zeta[1]+1)*/;

	if (dx[1] == 0)
		cout<<"ERROR in dx[1] \n";
//		printf("dy, np_coord, nm_coord %e %e\n", np->coord[1], nm->coord[1]);

	return;
}

void Element::insert_coord(HashTable* NodeTable) {

	Node* node = (Node*) (NodeTable->lookup(key));
	int i;
	for (i = 0; i < DIMENSION; i++)
		coord[i] = node->coord[i];

	return;
}

double max(double x, double y) {
	if (x > y)
		return (x);

	return (y);
}

double min(double x, double y) {
	if (x > y)
		return (y);

	return (x);
}

//x direction flux in current cell
void Element::xdirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]) {
	int i, j;
	double a, Vel;

	if (state_vars[0] < GEOFLOW_TINY) {

		for (i = 0; i < 3; i++)
			for (j = 0; j < NUM_STATE_VARS; j++)
				hfv[i][j] = 0.;

	} else {
		//state variables
		for (i = 0; i < NUM_STATE_VARS; i++)
			hfv[0][i] = state_vars[i] + d_state_vars[i] * dz;

//		if ((0.0 < Awet) && (Awet < 1.0))
//			for (i = 0; i < NUM_STATE_VARS; i++)
//				hfv[0][i] *= wetnessfactor;

		// Solid-phase velocity in x-dir
		Vel = hfv[0][1] / hfv[0][0];

		// sound-speed : a^2 = k_ap*h*g(3)
		a = sqrt(kactxy[0] * hfv[0][0] * gravity[2]);

		//fluxes
		hfv[1][0] = hfv[0][1];
		hfv[1][1] = hfv[0][1] * Vel + 0.5 * a * a * hfv[0][0];
		hfv[1][2] = hfv[0][2] * Vel;

		//wave speeds
		hfv[2][0] = Vel - a;
		hfv[2][1] = Vel;
		hfv[2][2] = Vel + a;
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			hrfv[i][j] = hfv[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hrfv[i][j]))
				cout << "flux is NAN" << endl;

}

//x direction flux in current cell
void Element::dual_xdirflux(Mat3x3& hfv, Mat3x3& flux_jac, Mat3x3& s_jac) {

	if ((prev_state_vars[0] < GEOFLOW_TINY)) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < NUM_STATE_VARS; j++) {
				hfv(i, j) = 0.;
				flux_jac(i, j) = 0.;
				s_jac(i, j) = 0.;
			}
	} else {

		for (int i = 0; i < NUM_STATE_VARS; i++)
			hfv(0, i) = prev_state_vars[i];

		double h_inv = 1. / hfv(0, 0);

		// Solid-phase velocity in x-dir
		double Vel = hfv(0, 1) * h_inv;

		// sound-speed : a^2 = k_ap*h*g(3)
		double a = sqrt(kactxy[0] * hfv(0, 0) * gravity[2]);

		//fluxes
		hfv(1, 0) = hfv(0, 1);
		hfv(1, 1) = hfv(0, 1) * Vel + 0.5 * a * a * hfv(0, 0);
		hfv(1, 2) = hfv(0, 2) * Vel;

		//wave speeds
		hfv(2, 0) = Vel - a;
		hfv(2, 1) = Vel;
		hfv(2, 2) = Vel + a;

		// jacobian of flux
		flux_jac(0, 0) = 0.;
		flux_jac(0, 1) = 1.;
		flux_jac(0, 2) = 0.;

		flux_jac(1, 0) = a * a - Vel * Vel;
		flux_jac(1, 1) = 2 * Vel;
		flux_jac(1, 2) = 0.;

		flux_jac(2, 0) = -Vel * hfv(0, 2) * h_inv;
		flux_jac(2, 1) = hfv(0, 2) * h_inv;
		flux_jac(2, 2) = Vel;

		//jacobian of speed
		s_jac(0, 0) = -.5 * a * h_inv - Vel * h_inv;
		s_jac(0, 1) = h_inv;
		s_jac(0, 2) = 0.;

		s_jac(1, 0) = -Vel * h_inv;
		s_jac(1, 1) = h_inv;
		s_jac(1, 2) = 0.;

		s_jac(2, 0) = .5 * a * h_inv - Vel * h_inv;
		s_jac(2, 1) = h_inv;
		s_jac(2, 2) = 0.;

	}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hfv(i, j)))
				cout << "flux is NAN" << endl;

}

//y direction flux in current cell
void Element::ydirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]) {
	int i, j;
	double Vel, a;

//the "update flux" values (hfv) are the fluxes used to update the solution,
// they may or may not be "reset" from their standard values based on whether
// or not the stopping criteria is triggering a change intended to cause the flow to stop.
	if (state_vars[0] < GEOFLOW_TINY) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < NUM_STATE_VARS; j++)
				hfv[i][j] = 0.0; //state variables

	} else {
		//state variables
		for (i = 0; i < NUM_STATE_VARS; i++)
			hfv[0][i] = state_vars[i] + d_state_vars[NUM_STATE_VARS + i] * dz;

//		if ((0.0 < Awet) && (Awet < 1.0))
//			for (i = 0; i < NUM_STATE_VARS; i++)
//				hfv[0][i] *= wetnessfactor;

		// a = speed of sound through the medium
		a = sqrt(kactxy[1] * hfv[0][0] * gravity[2]);

		// Solid-phase velocity in y-dir
		Vel = hfv[0][2] / hfv[0][0];

		//fluxes
		hfv[1][0] = hfv[0][2];
		hfv[1][1] = hfv[0][1] * Vel;
		hfv[1][2] = hfv[0][2] * Vel + 0.5 * a * a * hfv[0][0];

		//wave speeds
		hfv[2][0] = Vel - a;
		hfv[2][1] = Vel;
		hfv[2][2] = Vel + a;
	}

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			hrfv[i][j] = hfv[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hrfv[i][j]))
				cout << "flux is NAN" << endl;
}

//y direction flux in current cell
void Element::dual_ydirflux(Mat3x3& hfv, Mat3x3& flux_jac, Mat3x3& s_jac) {

	if ((prev_state_vars[0] < GEOFLOW_TINY)) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < NUM_STATE_VARS; j++) {
				hfv(i, j) = 0.;
				flux_jac(i, j) = 0.;
				s_jac(i, j) = 0.;
			}
	} else {		//state variables

		for (int i = 0; i < NUM_STATE_VARS; i++)
			hfv(0, i) = prev_state_vars[i];

		// a = speed of sound through the medium
		double a = sqrt(kactxy[1] * hfv(0, 0) * gravity[2]);

		double h_inv = 1. / hfv(0, 0);

		// Solid-phase velocity in y-dir
		double Vel = hfv(0, 2) * h_inv;

		//fluxes
		hfv(1, 0) = hfv(0, 2);
		hfv(1, 1) = hfv(0, 1) * Vel;
		hfv(1, 2) = hfv(0, 2) * Vel + 0.5 * a * a * hfv(0, 0);

		//wave speeds
		hfv(2, 0) = Vel - a;
		hfv(2, 1) = Vel;
		hfv(2, 2) = Vel + a;

		// jacobian of flux
		flux_jac(0, 0) = 0.;
		flux_jac(0, 1) = 0.;
		flux_jac(0, 2) = 1.;

		flux_jac(1, 0) = -Vel * hfv(0, 1) * h_inv;
		flux_jac(1, 1) = Vel;
		flux_jac(1, 2) = hfv(0, 1) * h_inv;

		flux_jac(2, 0) = a * a - Vel * Vel;
		flux_jac(2, 1) = 0.;
		flux_jac(2, 2) = 2 * Vel;

		//jacobian of speed
		s_jac(0, 0) = -.5 * a * h_inv - Vel * h_inv;
		s_jac(0, 1) = 0.;
		s_jac(0, 2) = h_inv;

		s_jac(1, 0) = -Vel * h_inv;
		s_jac(1, 1) = 0.;
		s_jac(1, 2) = h_inv;

		s_jac(2, 0) = .5 * a * h_inv - Vel * h_inv;
		s_jac(2, 1) = 0.;
		s_jac(2, 2) = h_inv;

	}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hfv(i, j)))
				cout << "flux is NAN" << endl;

}

//note z is not "z" but either x or y
void Element::zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    int order_flag, int dir, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS],
    Element *EmNeigh, double dt) {
	double dz = 0.0;
//	int ineigh = which_neighbor(EmNeigh->pass_key());
//
//	if (!((-1 < ineigh) && (ineigh < 8))) {
//		printf("zdirflux: ineigh=%d, dir=%d\n", ineigh, dir);
//		printf("this element******************************\n");
//		ElemBackgroundCheck2(El_Table, NodeTable, this, stdout);
//		fflush(stdout);
//		printf("Neigh element******************************\n");
//		ElemBackgroundCheck2(El_Table, NodeTable, EmNeigh, stdout);
//		fflush(stdout);
//		exit(ineigh);
//	}
//
//	double wetnessfactor = calc_elem_edge_wetness_factor(ineigh, dt);

	if (order_flag == 2)
		dz = (1.0 + dir % 2 - dir) * 0.5 * dx[dir % 2]; //+ or - 1/2 dx or dy

	if (dir % 2 == 0)
		xdirflux(matprops_ptr, dz, 1., hfv, hrfv);
	else if (dir % 2 == 1)
		ydirflux(matprops_ptr, dz, 1., hfv, hrfv);
	else {
		printf("zdirflux: direction %d not known\n", dir);
		exit(1);
	}
	return;
}

//note z is not "z" but either x or y
void Element::dual_zdirflux(int dir, Mat3x3& hfv, Mat3x3& flux_jac, Mat3x3& s_jac) {

	if (dir % 2 == 0)
		dual_xdirflux(hfv, flux_jac, s_jac);
	else if (dir % 2 == 1)
		dual_ydirflux(hfv, flux_jac, s_jac);
	else {
		printf("zdirflux: direction %d not known\n", dir);
		exit(1);
	}
	return;
}

void interflux_x(double hfvl[3][NUM_STATE_VARS], double hfvr[3][NUM_STATE_VARS],
    double interFlux[3], double k_ap_l, double k_ap_r, double g_l, double g_r) {
// this function computes the intermediate coeffiecient vector
// to impliment the approximate Riemann solver for non_conservative form of hyperbolic equation
// we need to a conservative path with some properties that are explained in the litrature
// we selected sai = left_state + S * (right_state - left_state)

//hfv: h=state variable, f=flux, v=wave speeds
//l="left" (the minus side), r="right" (the plus side)

//
	/*
	 ////    [ 0  k_ap*g*h-u^2  -u*v]        [0  a1 a2]
	 ////A = [ 1      2*u         v ]    A = [a3 a4 a5]
	 ////    [ 0       0          u ]        [0  0  a6]

	 */

	int i;
	double h_l, h_r, vel_x_l, vel_x_r, vel_y_l, vel_y_r;

	if (hfvl[0][0] < GEOFLOW_TINY && hfvr[0][0] > GEOFLOW_TINY) {
		h_l = hfvl[0][0];
		h_r = hfvr[0][0];
		vel_x_l = 0;
		vel_x_r = hfvr[0][2] / h_r;
		vel_y_l = 0;
		vel_y_r = hfvr[0][3] / h_r;
	} else if (hfvl[0][0] > GEOFLOW_TINY && hfvr[0][0] < GEOFLOW_TINY) {
		h_l = hfvl[0][0];
		h_r = hfvr[0][0];
		vel_x_l = hfvl[0][2] / h_l;
		vel_x_r = 0;
		vel_y_l = hfvl[0][3] / h_l;
		vel_y_r = 0;
	} else if (hfvl[0][0] < GEOFLOW_TINY && hfvr[0][0] < GEOFLOW_TINY) {
		h_l = hfvl[0][0];
		h_r = hfvr[0][0];
		vel_x_l = 0;
		vel_x_r = 0;
		vel_y_l = 0;
		vel_y_r = 0;
	} else {
		h_l = hfvl[0][0];
		h_r = hfvr[0][0];
		vel_x_l = hfvl[0][2] / h_l;
		vel_x_r = hfvr[0][2] / h_r;
		vel_y_l = hfvl[0][3] / h_l;
		vel_y_r = hfvr[0][3] / h_r;
	}

	interFlux[0] = 0.5
	    * (k_ap_r * g_r * h_r - vel_x_r * vel_x_r + k_ap_l * g_l * h_l - vel_x_l * vel_x_l)
	    * (hfvr[0][4] - hfvl[0][4])
	    - 0.5 * (vel_x_r * vel_y_r + vel_x_l * vel_y_l) * (hfvr[0][5] - hfvl[0][5]);

	interFlux[1] = (hfvr[0][1] - hfvl[0][1])
	    + 0.5 * 2 * (vel_x_r + vel_x_l) * (hfvr[0][4] - hfvl[0][4])
	    + 0.5 * (vel_y_r + vel_y_l) * (hfvr[0][5] - hfvl[0][5]);

	interFlux[2] = 0.5 * (vel_x_r + vel_x_l) * (hfvr[0][5] - hfvl[0][5]);

//  for(i=0;i<3;i++)
//    if(interFlux[i]>1e3)
//      printf("this is big in x %d \n", i);

	return;
}

void interflux_y(double hfvl[3][NUM_STATE_VARS], double hfvr[3][NUM_STATE_VARS],
    double interFlux[3], double k_ap_l, double k_ap_r, double g_l, double g_r) {
// this function computes the intermediate coef_matrix based on the Roe linearization matrix
// to impliment the approximate Riemann solver for non_conservative form of hyperbolic equation
// we need to a conservative path with some properties that are explained in the litrature
// we selected sai = left_state + S * (right_state - left_state)

//hfv: h=state variable, f=flux, v=wave speeds
//l="left" (the minus side), r="right" (the plus side)

//
	/*
	 ////    [ 0  -u*v  k_ap*g*h-v^2]        [0  b1 b2]
	 ////B = [ 0    v        0      ]    B = [0  b3  0]
	 ////    [ 1    u       2*v     ]        [b4 b5 b6]

	 */

	int i;
	double h_l, h_r, vel_x_l, vel_x_r, vel_y_l, vel_y_r;

	if (hfvl[0][0] < GEOFLOW_TINY && hfvr[0][0] > GEOFLOW_TINY) {
		h_l = 0;  //hfvl[0][0];
		h_r = hfvr[0][0];
		vel_x_l = 0;
		vel_x_r = hfvr[0][2] / h_r;
		vel_y_l = 0;
		vel_y_r = hfvr[0][3] / h_r;
	} else if (hfvl[0][0] > GEOFLOW_TINY && hfvr[0][0] < GEOFLOW_TINY) {
		h_l = hfvl[0][0];
		h_r = 0;  //hfvr[0][0];
		vel_x_l = hfvl[0][2] / h_l;
		vel_x_r = 0;
		vel_y_l = hfvl[0][3] / h_l;
		vel_y_r = 0;
	} else if (hfvl[0][0] < GEOFLOW_TINY && hfvr[0][0] < GEOFLOW_TINY) {
		h_l = 0;  //hfvl[0][0];
		h_r = 0;  //hfvr[0][0];
		vel_x_l = 0;
		vel_x_r = 0;
		vel_y_l = 0;
		vel_y_r = 0;
	} else {
		h_l = hfvl[0][0];
		h_r = hfvr[0][0];
		vel_x_l = hfvl[0][2] / h_l;
		vel_x_r = hfvr[0][2] / h_r;
		vel_y_l = hfvl[0][3] / h_l;
		vel_y_r = hfvr[0][3] / h_r;
	}

	interFlux[0] = -0.5 * (vel_x_r * vel_y_r + vel_x_l * vel_y_l) * (hfvr[0][4] - hfvl[0][4])
	    + 0.5 * (k_ap_r * g_r * h_r - vel_y_r * vel_y_r + k_ap_l * g_l * h_l - vel_y_l * vel_y_l)
	        * (hfvr[0][5] - hfvl[0][5]);

	interFlux[1] = 0.5 * (vel_y_r + vel_y_l) * (hfvr[0][4] - hfvl[0][4]);

	interFlux[2] = (hfvr[0][1] - hfvl[0][1]) + 0.5 * (vel_x_r + vel_x_l) * (hfvr[0][4] - hfvl[0][4])
	    + 0.5 * 2 * (vel_y_r + vel_y_l) * (hfvr[0][5] - hfvl[0][5]);

//  for(i=0;i<3;i++)
//    if(interFlux[i]>1e3)
//      printf("this is big in y %d \n",i);

	return;
}

//need move this to step.C
void riemannflux(double hfvl[3][NUM_STATE_VARS], double hfvr[3][NUM_STATE_VARS],
    double flux[NUM_STATE_VARS]) {
//hfv: h=state variable, f=flux, v=wave speeds
//l="left" (the minus side), r="right" (the plus side)
	int ivar, i;

	if ((hfvl[0][0] == 0.) && (hfvr[0][0] == 0.))
		for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			flux[ivar] = 0.;
	else {

		double sl, sr;
		if (hfvl[0][0] == 0.) {
			sl = min(0, 2.0 * hfvr[2][0] - hfvr[2][1]);
			sr = max(0, 2.0 * hfvr[2][2] - hfvr[2][1]);
		} else if (hfvr[0][0] == 0.) {
			sl = min(0, 2.0 * hfvl[2][0] - hfvl[2][1]);
			sr = max(0, 2.0 * hfvl[2][2] - hfvl[2][1]);
		} else {
			sl = min(0, min(hfvl[2][0], hfvr[2][0]));
			sr = max(0, max(hfvl[2][2], hfvr[2][2]));
		}
		if (sl >= 0.0)
			for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				flux[ivar] = hfvl[1][ivar];

		else if (sr <= 0.0)
			for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				flux[ivar] = hfvr[1][ivar];

		else
			for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				flux[ivar] = (sr * hfvl[1][ivar] - sl * hfvr[1][ivar]
				    + sl * sr * (hfvr[0][ivar] - hfvl[0][ivar])) / (sr - sl);
	}

}
//x direction flux in current cell
void Element::xdirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], ResFlag resflag) {
	int i, j;
	double a, Vel;

	if ((state_vars[0] < GEOFLOW_TINY && resflag.callflag == 0)
	    || (prev_state_vars[0] < GEOFLOW_TINY && resflag.callflag == 1)
	    || (resflag.lgft && resflag.callflag)) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < NUM_STATE_VARS; j++)
				hfv[i][j] = 0.0;
	} else {
		if (resflag.callflag == 0)
			//state variables
			for (i = 0; i < NUM_STATE_VARS; i++)
				hfv[0][i] = state_vars[i] + d_state_vars[i] * dz;
		else
			for (i = 0; i < NUM_STATE_VARS; i++)
				hfv[0][i] = prev_state_vars[i];

//		if ((0.0 < Awet) && (Awet < 1.0))
//			for (i = 0; i < NUM_STATE_VARS; i++)
//				hfv[0][i] *= wetnessfactor;

		// Solid-phase velocity in x-dir
		Vel = hfv[0][1] / hfv[0][0];

		// sound-speed : a^2 = k_ap*h*g(3)
		a = sqrt(kactxy[0] * hfv[0][0] * gravity[2]);

		//fluxes
		hfv[1][0] = hfv[0][1];
		hfv[1][1] = hfv[0][1] * Vel + 0.5 * a * a * hfv[0][0];
		hfv[1][2] = hfv[0][2] * Vel;

		//wave speeds
		hfv[2][0] = Vel - a;
		hfv[2][1] = Vel;
		hfv[2][2] = Vel + a;

	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			hrfv[i][j] = hfv[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hrfv[i][j]))
				cout << "flux is NAN" << endl;

	return;
}

//y direction flux in current cell
void Element::ydirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], ResFlag resflag) {
	int i, j;
	double Vel, a;

//the "update flux" values (hfv) are the fluxes used to update the solution,
// they may or may not be "reset" from their standard values based on whether
// or not the stopping criteria is triggering a change intended to cause the flow to stop.
	if ((state_vars[0] < GEOFLOW_TINY && resflag.callflag == 0)
	    || (prev_state_vars[0] < GEOFLOW_TINY && resflag.callflag == 1)
	    || (resflag.lgft && resflag.callflag)) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < NUM_STATE_VARS; j++)
				hfv[i][j] = 0.0; //state variables
	} else {
		//state variables
		if (resflag.callflag == 0)
			for (i = 0; i < NUM_STATE_VARS; i++)
				hfv[0][i] = state_vars[i] + d_state_vars[NUM_STATE_VARS + i] * dz;
		else
			for (i = 0; i < NUM_STATE_VARS; i++)
				hfv[0][i] = prev_state_vars[i];

//		if ((0.0 < Awet) && (Awet < 1.0))
//			for (i = 0; i < NUM_STATE_VARS; i++)
//				hfv[0][i] *= wetnessfactor;

		// a = speed of sound through the medium
		a = sqrt(kactxy[1] * hfv[0][0] * gravity[2]);

		// Solid-phase velocity in y-dir
		Vel = hfv[0][2] / hfv[0][0];

		//fluxes
		hfv[1][0] = hfv[0][2];
		hfv[1][1] = hfv[0][1] * Vel;
		hfv[1][2] = hfv[0][2] * Vel + 0.5 * a * a * hfv[0][0];

		//wave speeds
		hfv[2][0] = Vel - a;
		hfv[2][1] = Vel;
		hfv[2][2] = Vel + a;

	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			hrfv[i][j] = hfv[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hrfv[i][j]))
				cout << "flux is NAN" << endl;

	return;
}

//note z is not "z" but either x or y
void Element::zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    int order_flag, int dir, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS],
    Element *EmNeigh, double dt, ResFlag resflag) {
	double dz = 0.0;
//	int ineigh = which_neighbor(EmNeigh->pass_key());
//
//	if (!((-1 < ineigh) && (ineigh < 8))) {
//		printf("zdirflux: ineigh=%d, dir=%d\n", ineigh, dir);
//		printf("this element******************************\n");
//		ElemBackgroundCheck2(El_Table, NodeTable, this, stdout);
//		fflush(stdout);
//		printf("Neigh element******************************\n");
//		ElemBackgroundCheck2(El_Table, NodeTable, EmNeigh, stdout);
//		fflush(stdout);
//		exit(ineigh);
//	}
//
//	double wetnessfactor = calc_elem_edge_wetness_factor(ineigh, dt);

	if (order_flag == 2)
		dz = (1.0 + dir % 2 - dir) * 0.5 * dx[dir % 2]; //+ or - 1/2 dx or dy

	if (dir % 2 == 0)
		xdirflux(matprops_ptr, dz, 1., hfv, hrfv, resflag);
	else if (dir % 2 == 1)
		ydirflux(matprops_ptr, dz, 1., hfv, hrfv, resflag);
	else {
		printf("zdirflux: direction %d not known\n", dir);
		exit(1);
	}
	return;
}

void Element::calc_edge_states(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    int myid, double dt, int* order_flag, double *outflow, ResFlag lresflag, ResFlag rresflag) {
	Node *np, *np1, *np2, *nm, *nm1, *nm2;
	Element *elm1, *elm2;
	int side, zp, zm;
	int zp2, zm2; //positive_z_side_2 minus_z_side_2
	int zelmpos = -100, zelmpos_2 = -100;
	int ivar;
	const int SolFlux = 0, RefineFlux = 1; //indicate for which purpose flux is being calculated 0 is for solution and 1 is for refinement

//neutral resflag
	ResFlag nresflag;
	nresflag.callflag = rresflag.callflag;
	nresflag.lgft = 0;

	double hfv[3][NUM_STATE_VARS], hfv1[3][NUM_STATE_VARS], hfv2[3][NUM_STATE_VARS]; //update flux
	double hrfv[3][NUM_STATE_VARS], hrfv1[3][NUM_STATE_VARS], hrfv2[3][NUM_STATE_VARS]; //refinement flux

//ghost elements don't have nodes so you have to make temp storage for flux
	double ghostflux[NUM_STATE_VARS]; //, (*fluxptr)[NUM_STATE_VARS];
	*outflow = 0.0;

//  if (key[0]==3978454630 && key[1]==1717986917)
//    cout<<"this element is being checked"<<endl;

	for (side = 0; side < 2; side++) {
		zp = (positive_x_side + side) % 4;
		zm = (zp + 2) % 4;
		np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);

		if (neigh_proc[zp] == -1) {
			nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
			*outflow += (nm->flux[0]) * dx[!side];

			//outflow boundary conditions
			for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
				np->flux[ivar] = nm->flux[ivar];
				np->refinementflux[ivar] = nm->refinementflux[ivar];
			}
		} else if (neigh_proc[zp] != myid) {
			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
			elm1 = (Element*) El_Table->lookup(&neighbor[zp][0]);
			assert(elm1);

			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt, lresflag);
			elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this,
			    dt, rresflag);

			riemannflux(hfv, hfv1, np->flux);
			riemannflux(hrfv, hrfv1, np->refinementflux);

			elm2 = (Element*) El_Table->lookup(&neighbor[zp + 4][0]);
			assert(elm2);
			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt, lresflag);
			elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this,
			    dt, rresflag);

			//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
			if (neigh_proc[zp + 4] == myid) {
				zm2 = elm2->which_neighbor(pass_key()) % 4;
				nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zm2 + 4][0]);

				riemannflux(hfv, hfv2, nm2->flux);
				riemannflux(hrfv, hrfv2, nm2->refinementflux);

				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
					np->flux[ivar] = 0.5 * (np->flux[ivar] + nm2->flux[ivar]);
					np->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + nm2->refinementflux[ivar]);
				}
			} else {
				riemannflux(hfv, hfv2, ghostflux);
				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					np->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

				riemannflux(hrfv, hrfv2, ghostflux);
				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					np->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + ghostflux[ivar]);
			}
		} else {

			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
			elm1 = (Element*) El_Table->lookup(&neighbor[zp][0]);
			assert(elm1);

			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt, lresflag);
			elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this,
			    dt, rresflag);

			riemannflux(hfv, hfv1, np->flux);
			riemannflux(hrfv, hrfv1, np->refinementflux);

			/* CASE I
			 ------------------- -------------------
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |      this       np|         elm1      |
			 |        h          |         hp        |
			 |    kactxy_gz      |    kactxy_gz_n    |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 ------------------- -------------------
			 */

			/*
			 Case II
			 --------- ----------------------------
			 |         |         |                  |
			 |positive_z_side--->|<----zelmpos      |
			 |         | this  np|                  |
			 |         |   h     |                  |
			 |         |         |                  |
			 |---------|---------|nm1   elm1        |
			 |         |         |      hp          |
			 |         |         |                  |
			 |         | elm2 np2|                  |
			 |         |         |                  |
			 positive_z_side_2-->|                  |
			 --------- ----------------------------
			 */

			if (np->info == S_S_CON) {
				nm1 = NULL;
				np2 = NULL;

				zelmpos = elm1->which_neighbor(pass_key());
				assert(zelmpos > -1);
				nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos % 4 + 4][0]);

				elm2 = (Element*) El_Table->lookup(&elm1->neighbor[(zelmpos + 4) % 8][0]);
				assert(elm2);

				elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, elm2, dt,
				    rresflag);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, elm1, dt,
				    nresflag);

				if (*(elm1->get_neigh_proc() + (zelmpos + 4) % 8) == myid) {
					zp2 = elm2->which_neighbor(elm1->pass_key()) % 4;
					np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);
					riemannflux(hfv2, hfv1, np2->flux);

					riemannflux(hrfv2, hrfv1, np2->refinementflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = 0.5 * (np->flux[ivar] + np2->flux[ivar]);
						nm1->refinementflux[ivar] = 0.5
						    * (np->refinementflux[ivar] + np2->refinementflux[ivar]);
					}
				} else {

					riemannflux(hfv2, hfv1, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm1->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

					riemannflux(hrfv2, hrfv1, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm1->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + ghostflux[ivar]);
				}
			}

			/*  Case III
			 ------------------- --------- ---------
			 |                   |         |         |
			 |positive_z_side--->|<----zelmpos_2     |
			 |                   |nm2 elm2 |         |
			 |                   |     hp2 |         |
			 |                   |         |         |
			 |       this        |---------|---------
			 |        h          |         |         |
			 |                   |         |         |
			 |                   |nm1 elm1 |         |
			 |                   |     hp1 |         |
			 |                   |<----zelmpos       |
			 ------------------- --------- ---------
			 */

			else if (np->info == S_C_CON) {

				nm1 = NULL;
				nm2 = NULL;

				zelmpos = elm1->which_neighbor(pass_key()) % 4;
				nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos + 4][0]);

				elm2 = (Element*) (El_Table->lookup(&neighbor[zp + 4][0]));
				assert(elm2);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt,
				    lresflag);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this,
				    dt, rresflag);

				if (neigh_proc[zp + 4] == myid) {
					zelmpos_2 = elm2->which_neighbor(pass_key()) % 4;
					nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zelmpos_2 + 4][0]);
					riemannflux(hfv, hfv2, nm2->flux);
					riemannflux(hrfv, hrfv2, nm2->refinementflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = np->flux[ivar];
						np->flux[ivar] = 0.5 * (nm1->flux[ivar] + nm2->flux[ivar]);

						nm1->refinementflux[ivar] = np->refinementflux[ivar];
						np->refinementflux[ivar] = 0.5
						    * (nm1->refinementflux[ivar] + nm2->refinementflux[ivar]);
					}
				} else {
					riemannflux(hfv, hfv2, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = np->flux[ivar];
						np->flux[ivar] = 0.5 * (nm1->flux[ivar] + ghostflux[ivar]);
					}

					riemannflux(hrfv, hrfv2, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->refinementflux[ivar] = np->refinementflux[ivar];
						np->refinementflux[ivar] = 0.5 * (nm1->refinementflux[ivar] + ghostflux[ivar]);
					}
				}

			}

		}

		if (neigh_proc[zm] != myid) {

			if (neigh_proc[zm] == -1) {
				np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
				nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
				*outflow -= (np->flux[0]) * dx[!side];
				//outflow boundary conditions
				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
					nm->flux[ivar] = np->flux[ivar];
					nm->refinementflux[ivar] = np->refinementflux[ivar];
				}
			} else {
				/* if an interface is on the x-minus or y-minus side, need to
				 calculate those edgestates in this element */
				// x-minus side
				/*
				 interface
				 |
				 |
				 left cells-GHOST CELLS   |
				 (no need to        |
				 calculate         |
				 fluxes )         |
				 v
				 Case I
				 ------------------- -------------------
				 |                   |                   |
				 |                   |<--zm              |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 |       elm1        |nm      this       |
				 |        h          |         hp        |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 ------------------- -------------------
				 Case II
				 --------- -----------------------------
				 |         |         |                   |
				 |         |         |<--zm(z-minus side)|
				 |         |  elm1   |                   |
				 |         |   h     |                   |
				 |         |         |                   |
				 |---------|---------|nm    this         |
				 |         |         |       hp          |
				 |         |         |                   |
				 |         | elm2    |                   |
				 |         |   h2    |                   |
				 |         |         |                   |
				 --------- -----------------------------
				 */

				nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
				elm1 = (Element*) El_Table->lookup(&neighbor[zm][0]);
				assert(elm1);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm1, dt,
				    rresflag);
				elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, this, dt,
				    lresflag);
				riemannflux(hfv1, hfv, nm->flux);

				riemannflux(hrfv1, hrfv, nm->refinementflux);

				elm2 = (Element*) El_Table->lookup(&neighbor[zm + 4][0]);
				assert(elm2);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm2, dt,
				    rresflag);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, this, dt,
				    nresflag);

				//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
				if (neigh_proc[zm + 4] == myid) {
					zp2 = elm2->which_neighbor(pass_key()) % 4;
					np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);

					riemannflux(hfv2, hfv, np2->flux);

					riemannflux(hrfv2, hrfv, np2->refinementflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm->flux[ivar] = 0.5 * (nm->flux[ivar] + np2->flux[ivar]);
						nm->refinementflux[ivar] = 0.5 * (nm->refinementflux[ivar] + np2->refinementflux[ivar]);
					}
				} else {
					riemannflux(hfv2, hfv, ghostflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm->flux[ivar] = 0.5 * (nm->flux[ivar] + ghostflux[ivar]);

					riemannflux(hrfv2, hrfv, ghostflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm->refinementflux[ivar] = 0.5 * (nm->refinementflux[ivar] + ghostflux[ivar]);
				}

			}
		}

	}

}

template<class T>
void Element::calc_edge_states(HashTable* El_Table, HashTable* NodeTable, vector<T*>* x_elem_list,
    vector<T*>* y_elem_list, MatProps* matprops_ptr, int myid, double dt, int* order_flag,
    double *outflow) {
	Node *np, *np1, *np2, *nm, *nm1, *nm2;
	T *elm1, *elm2;
	int side, zp, zm;
	int zp2, zm2; //positive_z_side_2 minus_z_side_2
	int zelmpos = -100, zelmpos_2 = -100;
	int ivar;

	double hfv[3][NUM_STATE_VARS], hfv1[3][NUM_STATE_VARS], hfv2[3][NUM_STATE_VARS]; //update flux
	double hrfv[3][NUM_STATE_VARS], hrfv1[3][NUM_STATE_VARS], hrfv2[3][NUM_STATE_VARS]; //refinement flux

//ghost elements don't have nodes so you have to make temp storage for flux
	double ghostflux[NUM_STATE_VARS]; //, (*fluxptr)[NUM_STATE_VARS];
	*outflow = 0.0;

//  if (key[0]==3978454630 && key[1]==1717986917)
//    cout<<"this element is being checked"<<endl;

	for (side = 0; side < 2; side++) {

		vector<T*>* elem_list;
		if (side)
			elem_list = y_elem_list;
		else
			elem_list = x_elem_list;

		zp = (positive_x_side + side) % 4;
		zm = (zp + 2) % 4;
		np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);

		if (neigh_proc[zp] == -1) {
			elem_list->push_back((T*) this);
//			nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
//			*outflow += (nm->flux[0]) * dx[!side];
//
//			//outflow boundary conditions
//			for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
//				np->flux[ivar] = nm->flux[ivar];
//				np->refinementflux[ivar] = nm->refinementflux[ivar];
//			}
		} else if (neigh_proc[zp] != myid) {
			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
			elm1 = (T*) El_Table->lookup(&neighbor[zp][0]);
			assert(elm1);

			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt);
			elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this,
			    dt);

			riemannflux(hfv, hfv1, np->flux);
			riemannflux(hrfv, hrfv1, np->refinementflux);

			elm2 = (T*) El_Table->lookup(&neighbor[zp + 4][0]);
			assert(elm2);
			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt);
			elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this,
			    dt);

			//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
			if (neigh_proc[zp + 4] == myid) {
				zm2 = elm2->which_neighbor(pass_key()) % 4;
				nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zm2 + 4][0]);

				riemannflux(hfv, hfv2, nm2->flux);
				riemannflux(hrfv, hrfv2, nm2->refinementflux);

				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
					np->flux[ivar] = 0.5 * (np->flux[ivar] + nm2->flux[ivar]);
					np->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + nm2->refinementflux[ivar]);
				}
			} else {
				riemannflux(hfv, hfv2, ghostflux);
				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					np->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

				riemannflux(hrfv, hrfv2, ghostflux);
				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					np->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + ghostflux[ivar]);
			}
		} else {

			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
			elm1 = (T*) El_Table->lookup(&neighbor[zp][0]);
			assert(elm1);

			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt);
			elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this,
			    dt);

			riemannflux(hfv, hfv1, np->flux);
			riemannflux(hrfv, hrfv1, np->refinementflux);

			/* CASE I
			 ------------------- -------------------
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |      this       np|         elm1      |
			 |        h          |         hp        |
			 |    kactxy_gz      |    kactxy_gz_n    |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 ------------------- -------------------
			 */

			/*
			 Case II
			 --------- ----------------------------
			 |         |         |                  |
			 |positive_z_side--->|<----zelmpos      |
			 |         | this  np|                  |
			 |         |   h     |                  |
			 |         |         |                  |
			 |---------|---------|nm1   elm1        |
			 |         |         |      hp          |
			 |         |         |                  |
			 |         | elm2 np2|                  |
			 |         |         |                  |
			 positive_z_side_2-->|                  |
			 --------- ----------------------------
			 */

			if (np->info == S_S_CON) {
				nm1 = NULL;
				np2 = NULL;

				zelmpos = elm1->which_neighbor(pass_key());
				assert(zelmpos > -1);
				nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos % 4 + 4][0]);

				elm2 = (T*) El_Table->lookup(&elm1->neighbor[(zelmpos + 4) % 8][0]);
				assert(elm2);

				elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, elm2, dt);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, elm1, dt);

				if (*(elm1->get_neigh_proc() + (zelmpos + 4) % 8) == myid) {
					zp2 = elm2->which_neighbor(elm1->pass_key()) % 4;
					np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);
					riemannflux(hfv2, hfv1, np2->flux);

					riemannflux(hrfv2, hrfv1, np2->refinementflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = 0.5 * (np->flux[ivar] + np2->flux[ivar]);
						nm1->refinementflux[ivar] = 0.5
						    * (np->refinementflux[ivar] + np2->refinementflux[ivar]);
					}
				} else {

					riemannflux(hfv2, hfv1, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm1->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

					riemannflux(hrfv2, hrfv1, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm1->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + ghostflux[ivar]);
				}
			}

			/*  Case III
			 ------------------- --------- ---------
			 |                   |         |         |
			 |positive_z_side--->|<----zelmpos_2     |
			 |                   |nm2 elm2 |         |
			 |                   |     hp2 |         |
			 |                   |         |         |
			 |       this        |---------|---------
			 |        h          |         |         |
			 |                   |         |         |
			 |                   |nm1 elm1 |         |
			 |                   |     hp1 |         |
			 |                   |<----zelmpos       |
			 ------------------- --------- ---------
			 */

			else if (np->info == S_C_CON) {

				nm1 = NULL;
				nm2 = NULL;

				zelmpos = elm1->which_neighbor(pass_key()) % 4;
				nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos + 4][0]);

				elm2 = (T*) (El_Table->lookup(&neighbor[zp + 4][0]));
				assert(elm2);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this,
				    dt);

				if (neigh_proc[zp + 4] == myid) {
					zelmpos_2 = elm2->which_neighbor(pass_key()) % 4;
					nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zelmpos_2 + 4][0]);
					riemannflux(hfv, hfv2, nm2->flux);
					riemannflux(hrfv, hrfv2, nm2->refinementflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = np->flux[ivar];
						np->flux[ivar] = 0.5 * (nm1->flux[ivar] + nm2->flux[ivar]);

						nm1->refinementflux[ivar] = np->refinementflux[ivar];
						np->refinementflux[ivar] = 0.5
						    * (nm1->refinementflux[ivar] + nm2->refinementflux[ivar]);
					}
				} else {
					riemannflux(hfv, hfv2, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = np->flux[ivar];
						np->flux[ivar] = 0.5 * (nm1->flux[ivar] + ghostflux[ivar]);
					}

					riemannflux(hrfv, hrfv2, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->refinementflux[ivar] = np->refinementflux[ivar];
						np->refinementflux[ivar] = 0.5 * (nm1->refinementflux[ivar] + ghostflux[ivar]);
					}
				}

			}

		}

		if (neigh_proc[zm] != myid) {

			if (neigh_proc[zm] == -1) {

				elem_list->push_back((T*) this);
//				np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
//				nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
//				*outflow -= (np->flux[0]) * dx[!side];
//				//outflow boundary conditions
//				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
//					nm->flux[ivar] = np->flux[ivar];
//					nm->refinementflux[ivar] = np->refinementflux[ivar];
//				}
			} else {
				/* if an interface is on the x-minus or y-minus side, need to
				 calculate those edgestates in this element */
				// x-minus side
				/*
				 interface
				 |
				 |
				 left cells-GHOST CELLS   |
				 (no need to        |
				 calculate         |
				 fluxes )         |
				 v
				 Case I
				 ------------------- -------------------
				 |                   |                   |
				 |                   |<--zm              |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 |       elm1        |nm      this       |
				 |        h          |         hp        |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 ------------------- -------------------
				 Case II
				 --------- -----------------------------
				 |         |         |                   |
				 |         |         |<--zm(z-minus side)|
				 |         |  elm1   |                   |
				 |         |   h     |                   |
				 |         |         |                   |
				 |---------|---------|nm    this         |
				 |         |         |       hp          |
				 |         |         |                   |
				 |         | elm2    |                   |
				 |         |   h2    |                   |
				 |         |         |                   |
				 --------- -----------------------------
				 */

				nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
				elm1 = (T*) El_Table->lookup(&neighbor[zm][0]);
				assert(elm1);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm1, dt);
				elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, this, dt);
				riemannflux(hfv1, hfv, nm->flux);

				riemannflux(hrfv1, hrfv, nm->refinementflux);

				elm2 = (T*) El_Table->lookup(&neighbor[zm + 4][0]);
				assert(elm2);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm2, dt);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, this, dt);

				//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
				if (neigh_proc[zm + 4] == myid) {
					zp2 = elm2->which_neighbor(pass_key()) % 4;
					np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);

					riemannflux(hfv2, hfv, np2->flux);

					riemannflux(hrfv2, hrfv, np2->refinementflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm->flux[ivar] = 0.5 * (nm->flux[ivar] + np2->flux[ivar]);
						nm->refinementflux[ivar] = 0.5 * (nm->refinementflux[ivar] + np2->refinementflux[ivar]);
					}
				} else {
					riemannflux(hfv2, hfv, ghostflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm->flux[ivar] = 0.5 * (nm->flux[ivar] + ghostflux[ivar]);

					riemannflux(hrfv2, hrfv, ghostflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm->refinementflux[ivar] = 0.5 * (nm->refinementflux[ivar] + ghostflux[ivar]);
				}
			}
		}
	}
}


template<>
void Element::calc_edge_states<ErrorElem>(HashTable* El_Table, HashTable* NodeTable, vector<ErrorElem*>* x_elem_list,
    vector<ErrorElem*>* y_elem_list, MatProps* matprops_ptr, int myid, double dt, int* order_flag,
    double *outflow) {
	Node *np, *np1, *np2, *nm, *nm1, *nm2;
	ErrorElem *elm1, *elm2;
	int side, zp, zm;
	int zp2, zm2; //positive_z_side_2 minus_z_side_2
	int zelmpos = -100, zelmpos_2 = -100;
	int ivar;

	double hfv[3][NUM_STATE_VARS], hfv1[3][NUM_STATE_VARS], hfv2[3][NUM_STATE_VARS]; //update flux
	double hrfv[3][NUM_STATE_VARS], hrfv1[3][NUM_STATE_VARS], hrfv2[3][NUM_STATE_VARS]; //refinement flux

//ghost elements don't have nodes so you have to make temp storage for flux
	double ghostflux[NUM_STATE_VARS]; //, (*fluxptr)[NUM_STATE_VARS];
	*outflow = 0.0;

//  if (key[0]==3978454630 && key[1]==1717986917)
//    cout<<"this element is being checked"<<endl;

	for (side = 0; side < 2; side++) {

		vector<ErrorElem*>* elem_list;
		if (side)
			elem_list = y_elem_list;
		else
			elem_list = x_elem_list;

		zp = (positive_x_side + side) % 4;
		zm = (zp + 2) % 4;
		np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);

		if (neigh_proc[zp] == -1) {
			elem_list->push_back((ErrorElem*) this);
//			nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
//			*outflow += (nm->flux[0]) * dx[!side];
//
//			//outflow boundary conditions
//			for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
//				np->flux[ivar] = nm->flux[ivar];
//				np->refinementflux[ivar] = nm->refinementflux[ivar];
//			}
		} else if (neigh_proc[zp] != myid) {
			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
			elm1 = (ErrorElem*) El_Table->lookup(&neighbor[zp][0]);
			assert(elm1);

			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt);
			elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this,
			    dt);

			riemannflux(hfv, hfv1, np->flux);

			elm2 = (ErrorElem*) El_Table->lookup(&neighbor[zp + 4][0]);
			assert(elm2);
			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt);
			elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this,
			    dt);

			//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
			if (neigh_proc[zp + 4] == myid) {
				zm2 = elm2->which_neighbor(pass_key()) % 4;
				nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zm2 + 4][0]);

				riemannflux(hfv, hfv2, nm2->flux);

				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
					np->flux[ivar] = 0.5 * (np->flux[ivar] + nm2->flux[ivar]);
				}
			} else {
				riemannflux(hfv, hfv2, ghostflux);
				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					np->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);
			}
		} else {

			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
			elm1 = (ErrorElem*) El_Table->lookup(&neighbor[zp][0]);
			assert(elm1);

			zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt);
			elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this,
			    dt);

			riemannflux(hfv, hfv1, np->flux);

			/* CASE I
			 ------------------- -------------------
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |      this       np|         elm1      |
			 |        h          |         hp        |
			 |    kactxy_gz      |    kactxy_gz_n    |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 ------------------- -------------------
			 */

			/*
			 Case II
			 --------- ----------------------------
			 |         |         |                  |
			 |positive_z_side--->|<----zelmpos      |
			 |         | this  np|                  |
			 |         |   h     |                  |
			 |         |         |                  |
			 |---------|---------|nm1   elm1        |
			 |         |         |      hp          |
			 |         |         |                  |
			 |         | elm2 np2|                  |
			 |         |         |                  |
			 positive_z_side_2-->|                  |
			 --------- ----------------------------
			 */

			if (np->info == S_S_CON) {
				nm1 = NULL;
				np2 = NULL;

				zelmpos = elm1->which_neighbor(pass_key());
				assert(zelmpos > -1);
				nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos % 4 + 4][0]);

				elm2 = (ErrorElem*) El_Table->lookup(&elm1->neighbor[(zelmpos + 4) % 8][0]);
				assert(elm2);

				elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, elm2, dt);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, elm1, dt);

				if (*(elm1->get_neigh_proc() + (zelmpos + 4) % 8) == myid) {
					zp2 = elm2->which_neighbor(elm1->pass_key()) % 4;
					np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);
					riemannflux(hfv2, hfv1, np2->flux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = 0.5 * (np->flux[ivar] + np2->flux[ivar]);
					}
				} else {

					riemannflux(hfv2, hfv1, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm1->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

				}
			}

			/*  Case III
			 ------------------- --------- ---------
			 |                   |         |         |
			 |positive_z_side--->|<----zelmpos_2     |
			 |                   |nm2 elm2 |         |
			 |                   |     hp2 |         |
			 |                   |         |         |
			 |       this        |---------|---------
			 |        h          |         |         |
			 |                   |         |         |
			 |                   |nm1 elm1 |         |
			 |                   |     hp1 |         |
			 |                   |<----zelmpos       |
			 ------------------- --------- ---------
			 */

			else if (np->info == S_C_CON) {

				nm1 = NULL;
				nm2 = NULL;

				zelmpos = elm1->which_neighbor(pass_key()) % 4;
				nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos + 4][0]);

				elm2 = (ErrorElem*) (El_Table->lookup(&neighbor[zp + 4][0]));
				assert(elm2);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this,
				    dt);

				if (neigh_proc[zp + 4] == myid) {
					zelmpos_2 = elm2->which_neighbor(pass_key()) % 4;
					nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zelmpos_2 + 4][0]);
					riemannflux(hfv, hfv2, nm2->flux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = np->flux[ivar];
						np->flux[ivar] = 0.5 * (nm1->flux[ivar] + nm2->flux[ivar]);

					}
				} else {
					riemannflux(hfv, hfv2, ghostflux);
					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm1->flux[ivar] = np->flux[ivar];
						np->flux[ivar] = 0.5 * (nm1->flux[ivar] + ghostflux[ivar]);
					}

					riemannflux(hrfv, hrfv2, ghostflux);
				}

			}

		}

		if (neigh_proc[zm] != myid) {

			if (neigh_proc[zm] == -1) {

				elem_list->push_back((ErrorElem*) this);
//				np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
//				nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
//				*outflow -= (np->flux[0]) * dx[!side];
//				//outflow boundary conditions
//				for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
//					nm->flux[ivar] = np->flux[ivar];
//					nm->refinementflux[ivar] = np->refinementflux[ivar];
//				}
			} else {
				/* if an interface is on the x-minus or y-minus side, need to
				 calculate those edgestates in this element */
				// x-minus side
				/*
				 interface
				 |
				 |
				 left cells-GHOST CELLS   |
				 (no need to        |
				 calculate         |
				 fluxes )         |
				 v
				 Case I
				 ------------------- -------------------
				 |                   |                   |
				 |                   |<--zm              |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 |       elm1        |nm      this       |
				 |        h          |         hp        |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 |                   |                   |
				 ------------------- -------------------
				 Case II
				 --------- -----------------------------
				 |         |         |                   |
				 |         |         |<--zm(z-minus side)|
				 |         |  elm1   |                   |
				 |         |   h     |                   |
				 |         |         |                   |
				 |---------|---------|nm    this         |
				 |         |         |       hp          |
				 |         |         |                   |
				 |         | elm2    |                   |
				 |         |   h2    |                   |
				 |         |         |                   |
				 --------- -----------------------------
				 */

				nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
				elm1 = (ErrorElem*) El_Table->lookup(&neighbor[zm][0]);
				assert(elm1);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm1, dt);
				elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, this, dt);
				riemannflux(hfv1, hfv, nm->flux);

				elm2 = (ErrorElem*) El_Table->lookup(&neighbor[zm + 4][0]);
				assert(elm2);

				zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm2, dt);
				elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, this, dt);

				//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
				if (neigh_proc[zm + 4] == myid) {
					zp2 = elm2->which_neighbor(pass_key()) % 4;
					np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);

					riemannflux(hfv2, hfv, np2->flux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
						nm->flux[ivar] = 0.5 * (nm->flux[ivar] + np2->flux[ivar]);
					}
				} else {
					riemannflux(hfv2, hfv, ghostflux);

					for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
						nm->flux[ivar] = 0.5 * (nm->flux[ivar] + ghostflux[ivar]);

				}
			}
		}
	}
}


void Element::boundary_flux(HashTable* El_Table, HashTable* NodeTable, const int myid,
    const int side) {

	int zp = (positive_x_side + side) % 4;
	int zm = (zp + 2) % 4;
	Node* np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
	Node* nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);

	if (neigh_proc[zp] == -1) {

		// note the this may causes a bug if there is 2 level of refinement
		// for three layers of the cells adjacent to the boundary

		//in the boundary is in the right side of the cell

		//outflow boundary conditions
		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			np->flux[ivar] = nm->flux[ivar];

		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			np->refinementflux[ivar] = nm->refinementflux[ivar];

	}
	if (neigh_proc[zm] == -1) {

		//outflow boundary conditions
		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			nm->flux[ivar] = np->flux[ivar];

		//in above line we do not compute the flux on negative side and just
		// set it to positive flux which means outlet=inlet
		//by default fluxes_jac has been initialized to zero, but to be on the safe side

		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			nm->refinementflux[ivar] = np->refinementflux[ivar];

	}
}

void Element::eval_velocity(double xoffset, double yoffset, double Vel[]) {
	int i;
	double temp_state_vars[NUM_STATE_VARS];
	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		temp_state_vars[ivar] = state_vars[ivar] + d_state_vars[ivar] * xoffset + //distfromcenter[0]+
		    d_state_vars[NUM_STATE_VARS + ivar] * yoffset; //distfromcenter[1];

	for (i = 0; i < 4; i++)
		Vel[i] = 0;

	if (temp_state_vars[0] > GEOFLOW_TINY) {
		Vel[0] = temp_state_vars[1] / temp_state_vars[0];
		Vel[1] = temp_state_vars[2] / temp_state_vars[0];
		Vel[2] = temp_state_vars[1] / temp_state_vars[0];
		Vel[3] = temp_state_vars[2] / temp_state_vars[0];
	}
	return;
}

void Element::calc_gravity_vector(MatProps* matprops_ptr) {
	double max_slope = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1]);
	double max_angle = atan(max_slope);

	double down_slope_gravity = 9.8 * sin(max_angle);
	if (dabs(down_slope_gravity) > GEOFLOW_TINY) {
		gravity[0] = -down_slope_gravity * zeta[0] / max_slope;
		gravity[1] = -down_slope_gravity * zeta[1] / max_slope;
		gravity[2] = 9.8 * cos(max_angle);
	} else {
		gravity[0] = 0;
		gravity[1] = 0;
		gravity[2] = 9.8;
	}

	for (int i = 0; i < 3; i++)
		gravity[i] = gravity[i] / matprops_ptr->GRAVITY_SCALE;

	return;
}

int Element::determine_refinement(double target) {
	int flag = 0, i;

	if (state_vars[0] > target)
		flag = 1;
	return flag;
}

void Element::calc_d_gravity(HashTable* El_Table) {

	unsigned keyy[2] = { 3941335040, 0 };
	int aa = 0, bb = 1;
	if (key[0] == keyy[0] && key[1] == keyy[1])
		bb = aa;

	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}
	/* x direction */
	Element* ep = (Element*) (El_Table->lookup(&neighbor[xp][0]));
	Element* em = (Element*) (El_Table->lookup(&neighbor[xm][0]));
	int j;
	if (ep != NULL && em != NULL) {
		double dp, dm, dxp, dxm;
		dxp = ep->coord[0] - coord[0];
		dp = (ep->gravity[2] - gravity[2]) / dxp;
		dxm = coord[0] - em->coord[0];
		dm = (gravity[2] - em->gravity[2]) / dxm;

		d_gravity[0] = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
	} else if (em != NULL) {
		double dm, dxm;
		dxm = coord[0] - em->coord[0];
		d_gravity[0] = (gravity[2] - em->gravity[2]) / dxm;
	} else if (ep != NULL) {
		double dp, dxp;
		dxp = ep->coord[0] - coord[0];
		d_gravity[0] = (ep->gravity[2] - gravity[2]) / dxp;
	} else
//no neighbors on either side -- assume that the ground is flat
		d_gravity[0] = 0;

	/* y direction */
	ep = (Element*) (El_Table->lookup(&neighbor[yp][0]));
	em = (Element*) (El_Table->lookup(&neighbor[ym][0]));
	if (ep != NULL && em != NULL) {
		double dp, dm, dxp, dxm;
		dxp = ep->coord[1] - coord[1];
		dp = (ep->gravity[2] - gravity[2]) / dxp;
		dxm = coord[1] - em->coord[1];
		dm = (gravity[2] - em->gravity[2]) / dxm;

		d_gravity[1] = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
	} else if (em != NULL) {
		double dm, dxm;
		dxm = coord[1] - em->coord[1];
		d_gravity[1] = (gravity[2] - em->gravity[2]) / dxm;
	} else if (ep != NULL) {
		double dp, dxp;
		dxp = ep->coord[1] - coord[1];
		d_gravity[1] = (ep->gravity[2] - gravity[2]) / dxp;
	} else
//no neighbors on either side -- assume that the ground is flat
		d_gravity[1] = 0;

	return;
}

double* Element::get_zeta() {
	return zeta;
}

void Element::calc_topo_data(MatProps* matprops_ptr) {

	double resolution = (dx[0]/*/(zeta[0]*zeta[0]+1)*/+ dx[1]/*/(zeta[1]*zeta[1]+1)*/)
	    * (matprops_ptr->LENGTH_SCALE) / 2.0;  // element "size"
	double xcoord = coord[0] * (matprops_ptr->LENGTH_SCALE);
	double ycoord = coord[1] * (matprops_ptr->LENGTH_SCALE);
//double eldif = elevation;
	int i = Get_elevation(resolution, xcoord, ycoord, &elevation);
#ifdef PRINT_GIS_ERRORS
	if(i != 0) {
		printf("Error in Get_elevation(error code %d)\n", i);
		exit(1);
	}
#endif
	elevation = elevation / matprops_ptr->LENGTH_SCALE;
//eldif=(elevation-eldif)*matprops_ptr->LENGTH_SCALE;
//if(fabs(eldif)>1.0) printf("calc_topo_data() after-before=%g\n",eldif);
	i = Get_slope(resolution, xcoord, ycoord, zeta, (zeta + 1));
#ifdef PRINT_GIS_ERRORS
	if(i != 0) {
		printf("Error in Get_slope(error code %d)\n", i);
		exit(1);
	}
#endif
	i = Get_curvature(resolution, xcoord, ycoord, curvature, (curvature + 1));
#ifdef PRINT_GIS_ERRORS
	if(i != 0) {
		printf("Error in Get_curvature(error code %d)\n", i);
		exit(1);
	}
#endif
	curvature[0] = curvature[0] * (matprops_ptr->LENGTH_SCALE);
	curvature[1] = curvature[1] * (matprops_ptr->LENGTH_SCALE);

	if (matprops_ptr->material_count == 1) //only one material so don't need map
		material = 0; //GIS material id tag/index starts from 1
	else
//more than one material so need to get material from map
		Get_raster_id(resolution, xcoord, ycoord, &material);

	tan_bed_frict = matprops_ptr->tanbedfrict[material];

//flat plane!!!
	/*  elevation = 0;
	 zeta[0] = 0;
	 zeta[1] = 0;
	 curvature[0] = 0;
	 curvature[1] = 0; */

	return;
}

void Element::calc_flux_balance(HashTable* NodeTable) {
	int i, j;
	double flux[3] = { 0, 0, 0 };
	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}
	Node *nd_xp, *nd_xn, *nd_yp, *nd_yn;
	nd_xp = (Node*) NodeTable->lookup(node_key[xp + 4]);
	nd_xn = (Node*) NodeTable->lookup(node_key[xm + 4]);
	nd_yp = (Node*) NodeTable->lookup(node_key[yp + 4]);
	nd_yn = (Node*) NodeTable->lookup(node_key[ym + 4]);
	for (j = 0; j < 3; j++)
		flux[j] = dabs(nd_xp->refinementflux[j] - nd_xn->refinementflux[j])
		    + dabs(nd_yp->refinementflux[j] - nd_yn->refinementflux[j]);

	el_error[0] = 0;
	for (j = 0; j < NUM_STATE_VARS; j++)
		el_error[0] += flux[j];

	el_error[0] = 2. * el_error[0] * el_error[0] / (dx[0] + dx[1]) + WEIGHT_ADJUSTER; //W_A is so that elements with pile height = 0 have some weight.

	return;
}

// load-balancing stuff ...
void Element::put_lb_key(unsigned* in_key) {
	int i;
	for (i = 0; i < KEYLENGTH; i++)
		lb_key[i] = in_key[i];
	return;
}

void Element::copy_key_to_lb_key() {
	int i;
	for (i = 0; i < KEYLENGTH; i++)
		lb_key[i] = key[i];
	return;
}

void Element::put_coord(double* coord_in) {
	int i;
	for (i = 0; i < KEYLENGTH; i++)
		coord[i] = coord_in[i];
	return;
}

/*
 using elm_loc, which_son is calculated
 */
void Element::calc_which_son() {
	if (elm_loc[0] % 2 == 0) {
		if (elm_loc[1] % 2 == 0)
			which_son = 0;
		else
			which_son = 3;
	} else {
		if (elm_loc[1] % 2 == 0)
			which_son = 1;
		else
			which_son = 2;
	}

}

//should be full proof way to get key of opposite brother,
//as long as you know your own coord, dx, and which_son
void Element::find_opposite_brother(HashTable* El_Table) {

	if (opposite_brother_flag == 1)
		return;

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		brothers[(which_son + 2) % 4][ikey] = 0;
	unsigned nullkey[2] = { 0, 0 };
	if (!(compare_key(brothers[(which_son + 1) % 4], nullkey)
	    && compare_key(brothers[(which_son + 3) % 4], nullkey))) {
//use space filling curve to compute the key of opposite
//brother from it's bubble node coordinates
		double bro_norm_coord[2];
		unsigned nkey = KEYLENGTH;

		if ((which_son == 0) || (which_son == 3))
			bro_norm_coord[0] = El_Table->get_invdxrange()
			    * (coord[0] + dx[0] - *(El_Table->get_Xrange() + 0));
		else
			bro_norm_coord[0] = El_Table->get_invdxrange()
			    * (coord[0] - dx[0] - *(El_Table->get_Xrange() + 0));

		if ((which_son == 0) || (which_son == 1))
			bro_norm_coord[1] = El_Table->get_invdyrange()
			    * (coord[1] + dx[1] - *(El_Table->get_Yrange() + 0));
		else
			bro_norm_coord[1] = El_Table->get_invdyrange()
			    * (coord[1] - dx[1] - *(El_Table->get_Yrange() + 0));

		fhsfc2d_(bro_norm_coord, &nkey, brothers[(which_son + 2) % 4]);

		opposite_brother_flag = 1;
	}

	return;

}

int Element::if_pile_boundary(HashTable *ElemTable, double contour_height) {

	int ineigh;
	Element* ElemNeigh;

	assert(state_vars[0] >= 0.0);

	if (state_vars[0] >= contour_height) {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0) //don't check outside map boundary or duplicate neighbor
			    {
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				if (ElemNeigh == NULL) {
					printf(
					    "ElemNeigh==NULL ineigh=%d\n mykey   ={%u,%u} myprocess =%d generation=%d refined=%d adapted=%d\n",
					    ineigh, key[0], key[1], myprocess, generation, refined, adapted);
					printf(" neighbor={%u,%u} neigh_proc=%d neigh_gen =%d\n\n", neighbor[ineigh][0],
					    neighbor[ineigh][1], neigh_proc[ineigh], neigh_gen[ineigh]);
					fflush(stdout);
				}
				assert(ElemNeigh);
				if (*(ElemNeigh->get_state_vars() + 0) < contour_height)
					return (2); //inside of pileheight contour line
			}
	} else {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0) //don't check outside map boundary or duplicate neighbor
			    {
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				if (ElemNeigh == NULL) {
					printf(
					    "ElemNeigh==NULL\n mykey   ={%u,%u} myprocess =%d generation=%d refined=%d adapted=%d\n",
					    key[0], key[1], myprocess, generation, refined, adapted);
					printf(" neighbor={%u,%u} neigh_proc=%d neigh_gen =%d\n ineigh=%d\n", neighbor[ineigh][0],
					    neighbor[ineigh][1], neigh_proc[ineigh], neigh_gen[ineigh], ineigh);
					fflush(stdout);
				}
				assert(ElemNeigh);
				assert(*(ElemNeigh->get_state_vars() + 0) >= 0.0);
				if (*(ElemNeigh->get_state_vars() + 0) >= contour_height)
					return (1); //outside of pileheight contour line
			}
	}

	return (0); //not on pileheight contour line
}

int Element::if_source_boundary(HashTable *ElemTable) {

	int ineigh;
	Element* ElemNeigh;

	if (!(Influx[0] >= 0.0)) {
		printf("if_source_boundary() Influx[0]=%g\n", Influx[0]);
		fflush(stdout);
	}
	assert(Influx[0] >= 0.0); //currently mass sinks are not allowed

	if (Influx[0] > 0.0) {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0) { //don't check outside map boundary or duplicate neighbor
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				if (ElemNeigh == NULL) {
					printf(
					    "ElemNeigh==NULL\n mykey   ={%u,%u} myprocess =%d generation=%d refined=%d adapted=%d\n",
					    key[0], key[1], myprocess, generation, refined, adapted);
					printf(" neighbor={%u,%u} neigh_proc=%d neigh_gen =%d\n\n", neighbor[ineigh][0],
					    neighbor[ineigh][1], neigh_proc[ineigh], neigh_gen[ineigh]);
					fflush(stdout);
				}
				assert(ElemNeigh);
				if (*(ElemNeigh->get_influx() + 0) <= 0.0)
					return (2); //inside of line bounding area with a mass source
			}
//else if(neigh_proc[ineigh%4]==-1) return(2); //mass source on boundary of domain
	}

	else if (Influx[0] == 0.0) {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0.0) { //don't check outside map boundary or duplicate neighbor
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				if (ElemNeigh == NULL) {
					printf(
					    "ElemNeigh==NULL\n mykey   ={%u,%u} myprocess =%d generation=%d refined=%d adapted=%d\n",
					    key[0], key[1], myprocess, generation, refined, adapted);
					printf(" neighbor={%u,%u} neigh_proc=%d neigh_gen =%d\n\n", neighbor[ineigh][0],
					    neighbor[ineigh][1], neigh_proc[ineigh], neigh_gen[ineigh]);
					fflush(stdout);
				}
				assert(ElemNeigh);
				assert(*(ElemNeigh->get_influx() + 0) >= 0.0);
				if (*(ElemNeigh->get_influx() + 0) != 0.0)
					return (1); //outside of line bounding area with a mass source/sink
			}
	} else if (Influx[0] < 0.0) {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0.0) { //don't check outside map boundary or duplicate neighbor
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				if (ElemNeigh == NULL) {
					printf(
					    "ElemNeigh==NULL\n mykey   ={%u,%u} myprocess =%d generation=%d refined=%d adapted=%d\n",
					    key[0], key[1], myprocess, generation, refined, adapted);
					printf(" neighbor={%u,%u} neigh_proc=%d neigh_gen =%d\n\n", neighbor[ineigh][0],
					    neighbor[ineigh][1], neigh_proc[ineigh], neigh_gen[ineigh]);
					fflush(stdout);
				}
				assert(ElemNeigh);
				if (*(ElemNeigh->get_influx() + 0) >= 0.0)
					return (-1); //inside of line bounding area with a mass sink
			}
//else if(neigh_proc[ineigh%4]==-1) return(-1); //mass sink on boundary of domain
	}

	return (0); //not on line bounding area with mass source/sink
}

int Element::if_first_buffer_boundary(HashTable *ElemTable, double contour_height) {

	int ineigh;
	Element* ElemNeigh;
	int iffirstbuffer = 0;

	if (adapted <= 0)
		return (adapted - 1);

	assert(state_vars[0] >= 0.0);
	assert(Influx[0] >= 0.0);
	if ((state_vars[0] < contour_height) && (Influx[0] == 0.0)) {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0) { //don't check outside map boundary or duplicate neighbor
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				assert(ElemNeigh);
				if ((*(ElemNeigh->get_state_vars() + 0) >= contour_height)
				    || (*(ElemNeigh->get_influx() + 0) > 0.0)) {
					iffirstbuffer = 1;
					break;
				}
			}
	} else {
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0) //don't check outside map boundary or duplicate neighbor
			    {
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				assert(ElemNeigh);
				if ((*(ElemNeigh->get_state_vars() + 0) < contour_height)
				    && (*(ElemNeigh->get_influx() + 0) == 0.0)) {
					iffirstbuffer = 1;
					break;
				}
			}
	}

	if (iffirstbuffer) {
		if ((adapted >= NEWSON) || (generation == REFINE_LEVEL))
			return (2); //is a member of the buffer but doesn't need to be refined
		else
			return (1); //needs to be refined and some of its sons will be members
	}

	return (0);
}

int Element::if_next_buffer_boundary(HashTable *ElemTable, HashTable *NodeTable,
    double contour_height) {

	int ineigh;
	Element* ElemNeigh;
	int ifnextbuffer;
	ifnextbuffer = 0;
	if (adapted <= 0)
//GHOST element or element that should be deleted soon
		return (adapted - 1);

	if ((adapted != BUFFER) && //this element is not in the buffer
	    ((Influx[0] == 0.0))) //&& //this element is OUTSIDE the buffer layer "circle"
		for (ineigh = 0; ineigh < 8; ineigh++)
			if (neigh_proc[ineigh] >= 0) //don't check outside map boundary or duplicate neighbor
			    {
				ElemNeigh = (Element*) ElemTable->lookup(neighbor[ineigh]);
				if (!ElemNeigh) {
					printf("Elem={%10u,%10u} missing neighbor ineigh=%d {%10u,%10u}\n", key[0], key[1],
					    ineigh, neighbor[ineigh][0], neighbor[ineigh][1]);
					ElemBackgroundCheck(ElemTable, NodeTable, key, stdout);
					assert(ElemNeigh);
				}

				if ((abs(ElemNeigh->get_adapted_flag()) == BUFFER)
				    && (state_vars[0] <= *(ElemNeigh->get_state_vars()))) { //this element is next to a member of the old buffer layer
					ifnextbuffer = 1; //which means this element is a member of the next outer boundary of the buffer layer
					break;
				}
			}

	if (ifnextbuffer == 1) {
		if ((adapted >= NEWSON) || (generation == REFINE_LEVEL))
			return (2); //is a member of the buffer but doesn't need to be refined
		else
			return (1); //needs to be refined and some of its sons will be members
	}

	return (0);
}

template<class T>
T* Element::get_side_neighbor(HashTable *El_Table, int side) {

//       __________yp__________
//       |  |   |  |   |  |   |
//       |__|___|__|___|__|___|
//       |  |   |  |   |  |   |
//       |__| 9_|5_|1__| 8|___|
//       |  |   |      |  |   |
//	  xm |__|_2_|      |_4|___|   xp
//       |  |   |      |  |   |
//       |__|_6_|__ ___|_0|___|
//       |  |   |  |   |  |   |
//       |__|_10|_3|7__|11|___|
//       |  |   |  |   |  |   |
//       |__|___|__|___|__|___|
//	              ym
//  unsigned keyy[2]={3819350698,2863311530};
//	if (keyy[0]==key[0] && keyy[1]==key[1])
//		cout<<"why?"<<endl;

	T *neigh = NULL, *reserve;
	int xp = positive_x_side;
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	switch (side) {

		case 0: {
			if (neigh_proc[xp] != INIT)
				neigh = (T*) (El_Table->lookup(&neighbor[xp][0]));
		}
			break;

		case 1: {
			if (neigh_proc[yp] != INIT)
				neigh = (T*) (El_Table->lookup(&neighbor[yp][0]));
		}
			break;

		case 2: {
			if (neigh_proc[xm] != INIT)
				neigh = (T*) (El_Table->lookup(&neighbor[xm][0]));
		}
			break;

		case 3: {
			if (neigh_proc[ym] != INIT)
				neigh = (T*) (El_Table->lookup(&neighbor[ym][0]));
		}
			break;

		case 4:
			if (neigh_proc[xp + 4] < 0)
				neigh = (T*) (El_Table->lookup(&neighbor[xp][0]));
			else
				neigh = (T*) (El_Table->lookup(&neighbor[xp + 4][0]));
			break;

		case 5:

			if (neigh_proc[yp] < 0)
				neigh = (T*) (El_Table->lookup(&neighbor[yp][0]));
			else
				neigh = (T*) (El_Table->lookup(&neighbor[yp + 4][0]));
			break;

		case 6:

			if (neigh_proc[xm] < 0)
				neigh = (T*) (El_Table->lookup(&neighbor[xm][0]));
			else
				neigh = (T*) (El_Table->lookup(&neighbor[xm + 4][0]));
			break;

		case 7:

			if (neigh_proc[ym] < INIT)
				neigh = (T*) (El_Table->lookup(&neighbor[ym][0]));
			else
				neigh = (T*) (El_Table->lookup(&neighbor[ym + 4][0]));
			break;

		case 8:
			reserve = get_side_neighbor<T>(El_Table, 1);
			if (reserve)
				neigh = reserve->get_side_neighbor<T>(El_Table, 0);
			break;

		case 9:
			reserve = get_side_neighbor<T>(El_Table, 2);
			if (reserve)
				neigh = reserve->get_side_neighbor<T>(El_Table, 1);
			break;
		case 10:
			reserve = get_side_neighbor<T>(El_Table, 3);
			if (reserve)
				neigh = reserve->get_side_neighbor<T>(El_Table, 2);
			break;
		case 11:
			reserve = get_side_neighbor<T>(El_Table, 0);
			if (reserve)
				neigh = reserve->get_side_neighbor<T>(El_Table, 3);
			break;

		default:
			cout << "you entered an incorrect number for side" << endl;

	}

	return neigh;

}

void Element::for_link_temp1() {
	HashTable *El_Table;
	Element* elem;
	ErrorElem* elem1;
	vector<Element*> gg;
	vector<ErrorElem*> ff;
	MatProps* matprops_ptr;
	int aa;
	double bb;
//	DualElem* dualelem;
//	ErrorElem* errelem;

	elem->get_side_neighbor<Element>(El_Table, 1);
	elem->get_side_neighbor<DualElem>(El_Table, 1);
	elem->get_side_neighbor<ErrorElem>(El_Table, 1);

	elem->calc_edge_states<Element>(El_Table, El_Table, &gg, &gg, matprops_ptr, aa, bb, &aa, &bb);
	elem1->calc_edge_states<ErrorElem>(El_Table, El_Table, &ff, &ff, matprops_ptr, aa, bb, &aa, &bb);

}

void Element::gen_my_sons_key(HashTable* El_Table, unsigned son_key[4][KEYLENGTH]) {

	/*	---------
	 | 3 | 2 |
	 |   |   |
	 ---------
	 | 0 | 1 |
	 |   |   |
	 ---------*/

	static double XRange[KEYLENGTH], YRange[KEYLENGTH];
	double Xson[4] = { coord[0] - 0.25 * dx[0], coord[0] + 0.25 * dx[0], coord[0] + 0.25 * dx[0],
	    coord[0] - 0.25 * dx[0] };
	double Yson[4] = { coord[1] - 0.25 * dx[1], coord[1] - 0.25 * dx[1], coord[1] + 0.25 * dx[1],
	    coord[1] + 0.25 * dx[1] };
	double norm_coord[4][KEYLENGTH];
	unsigned nkey = 2;

	for (int i = 0; i < KEYLENGTH; ++i) {
		XRange[i] = *(El_Table->get_Xrange() + i);
		YRange[i] = *(El_Table->get_Yrange() + i);
	}

	for (int i = 0; i < 4; ++i) {
		norm_coord[i][0] = (Xson[i] - XRange[0]) / (XRange[1] - XRange[0]);
		norm_coord[i][1] = (Yson[i] - YRange[0]) / (YRange[1] - YRange[0]);
	}

	for (int i = 0; i < 4; ++i)
		fhsfc2d_(norm_coord[i], &nkey, son_key[i]);
}

void Element::gen_my_sons_key(HashTable* El_Table, unsigned* son_key) {

	/*	---------
	 | 3 | 2 |
	 |   |   |
	 ---------
	 | 0 | 1 |
	 |   |   |
	 ---------*/

	static double XRange[KEYLENGTH], YRange[KEYLENGTH];
	double Xson[4] = { coord[0] - 0.25 * dx[0], coord[0] + 0.25 * dx[0], coord[0] + 0.25 * dx[0],
	    coord[0] - 0.25 * dx[0] };
	double Yson[4] = { coord[1] - 0.25 * dx[1], coord[1] - 0.25 * dx[1], coord[1] + 0.25 * dx[1],
	    coord[1] + 0.25 * dx[1] };
	double norm_coord[4][KEYLENGTH];
	unsigned nkey = 2;

	for (int i = 0; i < KEYLENGTH; ++i) {
		XRange[i] = *(El_Table->get_Xrange() + i);
		YRange[i] = *(El_Table->get_Yrange() + i);
	}

	for (int i = 0; i < 4; ++i) {
		norm_coord[i][0] = (Xson[i] - XRange[0]) / (XRange[1] - XRange[0]);
		norm_coord[i][1] = (Yson[i] - YRange[0]) / (YRange[1] - YRange[0]);
		fhsfc2d_(norm_coord[i], &nkey, (son_key + i * KEYLENGTH));
	}

}

int Element::check_state(SolRec* solrec, HashTable* El_Table, int iter) {

	Solution* curr_sol = solrec->lookup(key, iter);

	if (curr_sol)
		for (int i = 0; i < NUM_STATE_VARS; i++)
			if (dabs(state_vars[i] - *(curr_sol->get_solution() + i)) > 1e-10)
				return 1;

	return 0;

}

void Element::write_elem_info(HashTable* NodeTable, char* filename, int iter, double dt) {

	double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
	double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];
	int ivar;
	int xp = positive_x_side;
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	Node* nxp = (Node*) NodeTable->lookup(&node_key[xp + 4][0]);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxp[ivar] = nxp->flux[ivar];

	Node* nyp = (Node*) NodeTable->lookup(&node_key[yp + 4][0]);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxyp[ivar] = nyp->flux[ivar];

	Node* nxm = (Node*) NodeTable->lookup(&node_key[xm + 4][0]);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxm[ivar] = nxm->flux[ivar];

	Node* nym = (Node*) NodeTable->lookup(&node_key[ym + 4][0]);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxym[ivar] = nym->flux[ivar];

	FILE *fp;
	fp = fopen("debugfile", "a");

	fprintf(fp,
	    "In %s time step %d with dt=%f dtdx=%f dtdy=%f kactx=%f , kacty=%f , x=%6f, y=%6f \n state vars are: \n",
	    filename, iter, dt, dt * dx[0], dt * dx[1], kactxy[0], kactxy[1], coord[0], coord[1]);

	for (int i = 0; i < NUM_STATE_VARS; i++)
		fprintf(fp, "%10e ", state_vars[i]);
	fprintf(fp, "\n prev state vars are: \n");
	for (int i = 0; i < NUM_STATE_VARS; i++)
		fprintf(fp, "%10e ", prev_state_vars[i]);
	fprintf(fp, "\n xp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxxp[state]);
	fprintf(fp, "\n xm: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxxm[state]);
	fprintf(fp, "\n yp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxyp[state]);
	fprintf(fp, "\n ym: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxym[state]);
	fprintf(fp, "\n dUdx: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state]);
	fprintf(fp, "\n dUdy: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state + NUM_STATE_VARS]);
	fprintf(fp, "\n");

	fclose(fp);

}

void Element::write_elem(gzFile& myfile) {

	gzwrite(myfile, &(generation), sizeof(int));
	gzwrite(myfile, (key), sizeof(unsigned) * 2);

	for (int i = 0; i < 8; ++i) {
		gzwrite(myfile, (node_key[i]), sizeof(unsigned) * 2);
		gzwrite(myfile, (neighbor[i]), sizeof(unsigned) * 2);
	}

	for (int i = 0; i < 4; ++i) {
		gzwrite(myfile, (son[i]), sizeof(unsigned) * 2);
		gzwrite(myfile, (brothers[i]), sizeof(unsigned) * 2);
	}

	gzwrite(myfile, (neigh_proc), sizeof(int) * 8);
	gzwrite(myfile, (neigh_gen), sizeof(int) * 8);
	gzwrite(myfile, &(refined), sizeof(int));
	gzwrite(myfile, &(adapted), sizeof(int));

	gzwrite(myfile, (coord), sizeof(double) * 2);
	gzwrite(myfile, (elm_loc), sizeof(int) * 2);
	gzwrite(myfile, (state_vars), sizeof(double) * 3);
	assert(state_vars[0] >= 0.);
	gzwrite(myfile, (prev_state_vars), sizeof(double) * 3);

	gzwrite(myfile, &(which_son), sizeof(int));
	gzwrite(myfile, &(positive_x_side), sizeof(int));
	gzwrite(myfile, dx, sizeof(double) * 2);
	gzwrite(myfile, curvature, sizeof(double) * 2);
	gzwrite(myfile, &elevation, sizeof(double));
	gzwrite(myfile, gravity, sizeof(double) * 3);
	gzwrite(myfile, d_gravity, sizeof(double) * 2);
	gzwrite(myfile, &opposite_brother_flag, sizeof(int));

	gzwrite(myfile, kactxy, sizeof(double) * 2);

	gzwrite(myfile, &material, sizeof(int));
}

void Element::save_elem(FILE* fp, FILE *fptxt) {
}

Element::Element(gzFile& myfile, HashTable* NodeTable, MatProps* matprops_ptr, int myid) {

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	for (int i = 0; i < NUM_STATE_VARS; i++)
		Influx[i] = 0.0;
	myprocess = myid;

	gzread(myfile, &(generation), sizeof(int));
	gzread(myfile, (key), sizeof(unsigned) * 2);

	for (int i = 0; i < 8; ++i) {
		gzread(myfile, (node_key[i]), sizeof(unsigned) * 2);
		gzread(myfile, (neighbor[i]), sizeof(unsigned) * 2);
	}

	for (int i = 0; i < 4; ++i) {
		gzread(myfile, (son[i]), sizeof(unsigned) * 2);
		gzread(myfile, (brothers[i]), sizeof(unsigned) * 2);
	}

	gzread(myfile, (neigh_proc), sizeof(int) * 8);
	gzread(myfile, (neigh_gen), sizeof(int) * 8);
	gzread(myfile, &(refined), sizeof(int));
	gzread(myfile, &(adapted), sizeof(int));

	gzread(myfile, (coord), sizeof(double) * 2);
	gzread(myfile, (elm_loc), sizeof(int) * 2);
	gzread(myfile, (state_vars), sizeof(double) * 3);
	assert(state_vars[0] >= 0.);
	gzread(myfile, (prev_state_vars), sizeof(double) * 3);

	//extra
	gzread(myfile, &(which_son), sizeof(int));
	gzread(myfile, &(positive_x_side), sizeof(int));
	gzread(myfile, dx, sizeof(double) * 2);
	gzread(myfile, curvature, sizeof(double) * 2);
	gzread(myfile, &elevation, sizeof(double));
	gzread(myfile, gravity, sizeof(double) * 3);
	gzread(myfile, d_gravity, sizeof(double) * 2);
	gzread(myfile, &opposite_brother_flag, sizeof(int));

	gzread(myfile, kactxy, sizeof(double) * 2);

	//super extra
	gzread(myfile, &material, sizeof(int));
	tan_bed_frict = matprops_ptr->tanbedfrict[material];

//	if (adapted > 0) {
//		calc_which_son();
//		find_positive_x_side(NodeTable);
//		calculate_dx(NodeTable);
//		calc_topo_data(matprops_ptr);
//		calc_gravity_vector(matprops_ptr);
//	}
}
