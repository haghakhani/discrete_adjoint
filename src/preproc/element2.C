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
 * $Id: element2.C 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include "element.h"
#include "useful_lib.h"

Element::Element() {
	int i;
	for (i = 0; i < 4; i++)
		neighbor[i] = NULL;
	opposite_brother = NULL;
}

void Element::setparameters(int di, Node* nodes[], int mat, int* elm_loc_in) {
	int i;
	elementid = di;
	for (i = 0; i < 8; i++)
		element_nodes[i] = nodes[i];

	for (i = 0; i < 4; i++) {
		neighbor[i] = 0;
	}

	material = mat;
	elm_loc[0] = elm_loc_in[0];
	elm_loc[1] = elm_loc_in[1];
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

void Element::create_m_node(double* max, double *min) {

	Node* newPtr = new Node();

	for (int j = 0; j < 2; j++) {
		newPtr->node_coord[j] = 0;
		for (int i = 0; i < 4; i++)
			newPtr->node_coord[j] += 0.25 * element_nodes[i]->node_coord[j];
	}
	newPtr->written = 0;

	element_nodes[8] = newPtr;
	element_nodes[8]->nodeid = -elementid;

	if (element_nodes[8]->node_coord[0] > *max)
		*max = element_nodes[8]->node_coord[0];
	if (element_nodes[8]->node_coord[1] > *(max + 1))
		*(max + 1) = element_nodes[8]->node_coord[1];
	if (element_nodes[8]->node_coord[0] < *min)
		*min = element_nodes[8]->node_coord[0];
	if (element_nodes[8]->node_coord[1] < *(min + 1))
		*(min + 1) = element_nodes[8]->node_coord[1];

	return;
}

void Element::determine_neighbors(int array_loc, Element* el) {
	int i, neigh_loc;

	for (i = 4; i < 8; i++) {
		neigh_loc = element_nodes[i]->get_element_array_loc(array_loc);
		if (neigh_loc >= 0)
			neighbor[i - 4] = (el + neigh_loc);
	}
}

//*************************finding the opposite brother***********************
void Element::determine_opposite_brother() {
	if (which_son != 3) {
		if (neighbor[which_son + 1] != NULL) {
			if (which_son != 2)
				opposite_brother = neighbor[which_son + 1]->get_neighbors(which_son + 2);
			else
				opposite_brother = neighbor[which_son + 1]->get_neighbors(0);
		}
	} else {
		if (neighbor[0] != NULL)
			opposite_brother = neighbor[0]->get_neighbors(1);
	}

}

Element* Element::get_neighbors(int side) {
	return neighbor[side];
}

void Element::myproc(int np, int i, int ec) {
	myprocess = i / (ec / np);

	if (i >= (ec / np) * np)
		myprocess = np - 1;
}

void Element::write_element_data_bin(FILE *fp) {
	int i;
	for (i = 0; i < 9; i++) {
		fwriteU(fp, element_nodes[i]->key[0]);
		fwriteU(fp, element_nodes[i]->key[1]);
	}

	for (i = 0; i < 4; i++)
		if (neighbor[i] != 0) {
			fwriteI(fp, neighbor[i]->myprocess);
			fwriteU(fp, neighbor[i]->element_nodes[8]->key[0]);
			fwriteU(fp, neighbor[i]->element_nodes[8]->key[1]);
		} else
			fwriteI(fp, -1);

	fwriteI(fp, elm_loc[0]);
	fwriteI(fp, elm_loc[1]);

	if (opposite_brother == NULL) {
		fwriteI(fp, 0);
		fwriteI(fp, 0);
	} else {
		fwriteU(fp, *(opposite_brother->pass_key()));
		fwriteU(fp, *(opposite_brother->pass_key() + 1));
	}

	fwriteI(fp, material);
}

void Element::reset_written_flag() {

	for (int i = 0; i < 8; i++)
		if (element_nodes[i]->written != 0)
			element_nodes[i]->written = 0;

}

unsigned* Element::pass_key() {
	unsigned* mykey = element_nodes[8]->get_key();

	return mykey;
}
