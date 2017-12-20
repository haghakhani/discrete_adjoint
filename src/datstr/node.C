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
 * $Id: node.C 156 2007-06-27 20:33:08Z dkumar $ 
 */
//#define DEBUG_SAVE_NODE
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef NODE_C
#define NODE_C
#include "../header/node.h"
#include "../gisapi/GisApi.h"
#include "../header/properties.h"
#include "../header/dualmesh.h"

//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR
#include <mpi.h>
#include <assert.h>

Node::Node() {
	printf("creating a node without setting its values\n");
}

Node::Node(unsigned* keyi, double* coordi, MatProps* matprops_ptr) {
	int i;
	id = 0; /* keith added this so save_node wouldn't write an uninitialized
	 variable and valgrind wouldn't flag an error.  id is used in
	 ../repartition/BSFC_update_and_send_elements.C */

	info = INIT;

	for (i = 0; i < DIMENSION; i++)
		coord[i] = *(coordi + i);

	for (i = 0; i < KEYLENGTH; i++)
		key[i] = *(keyi + i);

	zero_flux();
	// find the max resolution of the GIS info and then get the elevation at this node
	double resolution = 0;
	double xcoord = coord[0] * (matprops_ptr->LENGTH_SCALE);
	double ycoord = coord[1] * (matprops_ptr->LENGTH_SCALE);
	i = Get_max_resolution(&resolution);
	if (i != 0) {
		printf("error in Get_max_resolution\n");
		exit(1);
	}
	i = Get_elevation(resolution, xcoord, ycoord, &elevation);
	if (i != 0) {
		printf("error in Get_elevation(%d) r=%g (x,y)=(%g,%g) e=%g\n", i, resolution, xcoord, ycoord,
		    elevation);
		exit(1);
	}
	elevation = elevation / matprops_ptr->LENGTH_SCALE;
	/*  if((unsigned) 1548032885 == key[0])
	 printf("created the missing node...\n"); */

	/*
	 if((coord[0]==0)||(coord[1]==0)){
	 printf("node={%u,%u} coord=(%20g,%20g)\n",key[0],key[1],coord[0],coord[1]);
	 assert(coord[0]*coord[1]);
	 }
	 */
}

Node::Node(unsigned* keyi, double* coordi, int inf, MatProps* matprops_ptr)  //for refined
    {
	int i;
	id = 0; /* keith added this so save_node wouldn't write an uninitialized
	 variable and valgrind wouldn't flag an error.  id is used in
	 ../repartition/BSFC_update_and_send_elements.C */

	info = INIT;

	for (i = 0; i < DIMENSION; i++)
		coord[i] = *(coordi + i);

	for (i = 0; i < KEYLENGTH; i++)
		key[i] = *(keyi + i);

	info = inf;
	//geoflow stuff
	zero_flux();
	// find the max resolution of the GIS info and then get the elevation at this node
	double resolution = 0;
	i = Get_max_resolution(&resolution);
	if (i != 0) {
		printf("error in Get_max_resolution\n");
		exit(1);
	}
	double xcoord = coord[0] * (matprops_ptr->LENGTH_SCALE);
	double ycoord = coord[1] * (matprops_ptr->LENGTH_SCALE);
	i = Get_elevation(resolution, xcoord, ycoord, &elevation);
	if (i != 0) {
		printf("error in Get_elevation\n");
		exit(1);
	}
	elevation = elevation / matprops_ptr->LENGTH_SCALE;
	/*  if((unsigned) 1548032885 == key[0])
	 printf("created the missing node111111...\n");*/

	/*
	 if((coord[0]==0)||(coord[1]==0)){
	 printf("node={%u,%u} coord=(%20g,%20g)\n",key[0],key[1],coord[0],coord[1]);
	 assert(coord[0]*coord[1]);
	 }
	 */
	return;
}

Node::Node(unsigned* keyi, double* coordi, int inf, double elev, int yada) {
	int i;
	id = 0; /* keith added this so save_node wouldn't write an uninitialized
	 variable and valgrind wouldn't flag an error.  id is used in
	 ../repartition/BSFC_update_and_send_elements.C */

	info = INIT;

	for (i = 0; i < DIMENSION; i++)
		coord[i] = *(coordi + i);

	for (i = 0; i < KEYLENGTH; i++)
		key[i] = *(keyi + i);

	info = inf;
	//geoflow stuff
	zero_flux();
	elevation = elev;
	/*
	 if((coord[0]==0)||(coord[1]==0)){
	 printf("inode=%d node={%u,%u} coord=(%20g,%20g)\n",yada,key[0],key[1],coord[0],coord[1]);
	 assert(coord[0]*coord[1]);
	 }
	 */
	return;
}

Node::Node(Node* node) {

	id = node->id;

	num_assoc_elem = node->num_assoc_elem;

	info = node->info;

	elevation = node->elevation;

	for (int i = 0; i < DIMENSION; ++i) {

		coord[i] = node->coord[i];

		key[i] = node->key[i];

	}

	for (int i = 0; i < NUM_STATE_VARS; ++i) {

		flux[i] = node->flux[i];

		refinementflux[i] = node->refinementflux[i];
	}
}

Node::Node(gzFile& myfile, MatProps* matprops_ptr) {

	gzread(myfile, &(id), sizeof(int));
	gzread(myfile, &(info), sizeof(int));
	gzread(myfile, (key), sizeof(unsigned) * 2);
	gzread(myfile, (coord), sizeof(double) * 2);
	gzread(myfile, &(elevation), sizeof(double));

//	double resolution = 0;
//	int i = Get_max_resolution(&resolution);
//	if (i != 0) {
//		printf("error in Get_max_resolution\n");
//		exit(1);
//	}
//	double xcoord = coord[0] * (matprops_ptr->LENGTH_SCALE);
//	double ycoord = coord[1] * (matprops_ptr->LENGTH_SCALE);
//	i = Get_elevation(resolution, xcoord, ycoord, &elevation);
//	if (i != 0) {
//		printf("error in Get_elevation\n");
//		exit(1);
//	}
//	elevation = elevation / matprops_ptr->LENGTH_SCALE;

	zero_flux();
	num_assoc_elem = 0;
}

Node::Node(const Node_minimal* node_minimal){

	id = node_minimal->id;
	info = node_minimal->info;
	for (int i=0;i<DIMENSION;++i){
		key[i]= node_minimal->key[i];
		coord[i] = node_minimal->coord[i];
	}
	elevation = node_minimal->elevation;
	zero_flux();
	num_assoc_elem = 0;
}

void Node::set_parameters(int inf) {
	info = inf;
}
void Node::putinfo(int in) {
	info = in;
}

void Node::zero_flux() {
	int i;
	for (i = 0; i < NUM_STATE_VARS; i++)
		flux[i] = refinementflux[i] = 0.0;
}

void Node::set_elevation(MatProps* matprops_ptr) {
	double resolution = 0;
	double xcoord = coord[0] * (matprops_ptr->LENGTH_SCALE);
	double ycoord = coord[1] * (matprops_ptr->LENGTH_SCALE);
	int i = Get_max_resolution(&resolution);
	if (i != 0) {
		printf("error in Get_max_resolution\n");
		exit(1);
	}
	i = Get_elevation(resolution, xcoord, ycoord, &elevation);
	if (i != 0) {
		printf("error in Get_elevation\n");
		exit(1);
	}
	elevation = elevation / matprops_ptr->LENGTH_SCALE;

}

void Node::write_node(gzFile myfile) {

	gzwrite(myfile, &(id), sizeof(int));
	gzwrite(myfile, &(info), sizeof(int));
	gzwrite(myfile, (key), sizeof(unsigned) * 2);
	gzwrite(myfile, (coord), sizeof(double) * 2);
	gzwrite(myfile, &(elevation), sizeof(double));

}

int Node::get_info(){
	return info;
}
#endif

