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
 * $Id: node.h 152 2007-06-27 20:29:54Z dkumar $ 
 */

#ifndef NODE_H
#define NODE_H

#include "properties.h"
#include "constant.h"
#include "hashtab.h"
#include "struct.h"
//#include "jacobian.h"
class SolRec;
class Jacobian;
class Element;
class DualElem;
class Node_minimal;
class Node {

	friend class Element;

	friend class DualElem;

	friend class ErrorElem;

	friend void correct(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
	    FluxProps *fluxprops, TimeProps *timeprops, Element *EmTemp);

	friend void calc_jacobian(MeshCTX* meshctx, PropCTX* propctx);

	friend void get_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key,
	    MatProps* matprops_ptr, int myid, double fluxold[4][NUM_STATE_VARS]);

	friend void AssertMeshErrorFree(HashTable *El_Table, HashTable* NodeTable, int numprocs, int myid,
	    double loc);

	friend void ElemBackgroundCheck(HashTable* El_Table, HashTable* NodeTable, unsigned *debugkey,
	    FILE *fp);

	friend void ElemBackgroundCheck2(HashTable* El_Table, HashTable* NodeTable, void *EmDebug,
	    FILE *fp);

	friend void NodeBackgroundCheck(HashTable* El_Table, HashTable* NodeTable, unsigned *debugkey,
	    FILE *fp);

	friend void delete_oldsons(HashTable* El_Table, HashTable* NodeTable, int myid, void *EmFather);

	friend void refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int numprocs, int myid,
	    void* RefinedList, TimeProps* timeprops_ptr);

	friend void restore(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El, int effelement,
	    int j, double increment, double fluxold[4][NUM_STATE_VARS],
	    double d_state_vars_old[DIMENSION * NUM_STATE_VARS]);

	friend void record_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key,
	    MatProps* matprops_ptr, int myid, double fluxold[4][NUM_STATE_VARS]);

	friend int checkElement(HashTable *El_Table, HashTable *NodeTable, double *max, unsigned *key);

	/*
	 friend void unrefine_interp_neigh_update(HashTable* El_Table,
	 HashTable* NodeTable, int nump,
	 int myid, int NumOtherProcUpdate,
	 Element **OtherProcUpdate);
	 */
	friend void unrefine_interp_neigh_update(HashTable* El_Table, HashTable* NodeTable, int nump,
	    int myid, void* OtherProcUpdate);

	//friend void Pack_element(Element* sendel, ElemPack** elemptr, HashTable* HT_Node_Ptr, int);

	friend void Pack_element(void *sendel, ElemPack* elem, HashTable* HT_Node_Ptr, int);

	friend void destroy_element(void *r_element, HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr);

	friend void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
	    double* e_error, MatProps *matprops);

	friend void uniform_refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int nump, int myid, void* RL,
	    TimeProps* timeprops_ptr);

	friend void adjust_node_info(MeshCTX* meshctx, PropCTX* propctx);

public:
	//! this is the constructor that creates a node when the initial grid is read in
	Node(unsigned *keyi, double *coordi, MatProps *matprops_ptr);

	//! this is the constructor that creates bubble and edge nodes for son Elements when the father Element is refined
	Node(unsigned *keyi, double *coordi, int inf, MatProps *matprops_ptr);/*for refined*/

	//! this is the node constructor that is called in construct_el() in update_element_info.C
	Node(unsigned* keyi, double* coordi, int inf, double elev, int yada);

	//! this is the constructor that recreates/restores a node that was saved in a restart file.
	Node(gzFile& myfile, MatProps* matprops_ptr);

	Node(const Node_minimal* node_minimal);

	//! constructor that creates a node without setting any of its values
	Node();

	//! copy constructor
	Node(Node* node);

	~Node() {
	}
	;

	//! this function writes all of one Node's data necessary for restart to a file in a single fwrite statement
	void write_node(gzFile myfile);

	//! this function is legacy afeapi, it is extraneous for the finite difference/volume version of titan, however it appears once in htflush.C
	void putdof(int lower, int up);

	//! this function is legacy afeapi, it is extraneous for the finite difference/volume version of titan
	int* getdof();

	//! this function is legacy afeapi
	void putglnum(int);

	//! this function is legacy afeapi
	void putsol(double* s);

	//! this function is legacy afeapi
	int getglnum();

	//! this function returns the node type, the options are listed in constant.h and include: NODEINIT, CORNER, BUBBLE, SIDE, CONSTRAINED, S_C_CON, S_S_CON, ASSIGNED,and UNASSIGNED.
	int getinfo();

	//! this function sets the node type, the options are listed in constant.h and include: NODEINIT, CORNER, BUBBLE, SIDE, CONSTRAINED, S_C_CON, S_S_CON, ASSIGNED,and UNASSIGNED.
	void putinfo(int in);

	//! this function returns the node key, a key is a single number that is 2 unsigned variables long and is used to access the pointer to a Node or Element through the HashTable
	unsigned* pass_key();

	//! this function returns the global x and y coordinates of the node, in the finite difference version of Titan this is not always reliable, use the coordinates of the element instead.  It is reliable in the Discontinuous Galerkin version of Titan however.
	double* get_coord();

	//! this is legacy afeapi and is not used
	double* getsol();

	//! this is legacy afeapi and is not used
	int get_order();

	//! this is legacy afeapi and is not used
	void put_order(int);

	//! this is legacy afeapi and is not used
	void increase_order();

	//! this function sets the node information and order, node order is legacy afeapi but node information is currently used, this function is called in update_element_info.C, another distict function with a similar name refined_neighbor::set_parameters also existis and is used in updatenei.C, these should not be confused
	void set_parameters(int inf);

	//! this is legacy afeapi and is not used
	int get_reconstructed();

	//! this is legacy afeapi and is not used at all except once in htflush.C
	void put_reconstructed(int);

	//! this is legacy afeapi and is not used at all except once in htflush.C
	int get_sol_deleted();

	//! this is legacy afeapi and is not used at all except once in htflush.C
	void put_sol_deleted(int flag);

	//! this function sets the id of a node, it is used in repartitioning,
	void put_id(int id_in);

	//! this function returns the id of a node, it is used in repartitioning,
	int get_id();

	int get_info();

	//! this function returns the vector of fluxes stored in an edge node between elements
	double* get_flux();

	//! this function zeros the flux used during refinement, this is only distinct from the regular flux if flux velocity is being zero'd because the experimental stopping criteria says it should be.  This feature is disabled by default.  Keith implemented it.
	void zero_flux();

	//! this function returns the elevation of a node, in the finite difference version of Titan this is not always reliable, use the elevation of the element instead.  It is reliable in the Discontinuous Galerkin version of Titan however.
	double get_elevation();

	//! this function sets the elevation of a node
	void set_elevation(MatProps* matprops_ptr);

	//! this function stores the number of elements associated with this node
	void put_num_assoc_elem(int numin);

	//! this function returns the number of elements associated with this node
	int get_num_assoc_elem();

protected:
	//! used in delete_unused_nodes_and_elements() function
	int id;

	//! the number of associated elements, it is used in extraneous node
	//deletion and debugging function AssertMeshErrorFree()
	int num_assoc_elem;

	//! says what type of node this is see the comments of Node::get_info()
	int info;

	//! the global x and y coordinates of the node
	double coord[DIMENSION];

	//! this is the node key, a key is a single number that is 2 unsigned variables long and is used to access the pointer to a Node or Element through the HashTable
	unsigned key[KEYLENGTH];

	//! this elevation should currently be the GIS elevation at the finest "scale"
	double elevation;

	//! these are the so called "regular fluxes" that is the ones that are used to update the elements, assume that element normal is parallel to either the x or y axis, Keith is the one who introduced a distinction between regular and refinement fluxes for use with the stopping criteria, this distinction is disabled by default
	double flux[NUM_STATE_VARS];

	//! the "refinement flux" is necessary when using the stopping criteria to reset the "regular" fluxes to what they would be if velocity was zero in the cell(s) involved.  The refinement flux is what the flux would have been if it had not been reset, they are needed since refinement is based on fluxes (and also pileheight gradient but that's not relevant here) Keith is the one who introduced a distinction between regular and refinement fluxes for use with the stopping criteria, this distinction is disabled by default.
	double refinementflux[NUM_STATE_VARS];
};

inline int Node::getinfo() {
	return info;
}
;

inline double* Node::get_coord() {
	return coord;
}
;

inline void Node::put_id(int id_in) {
	id = id_in;
}
;

inline int Node::get_id() {
	return id;
}
;

inline double* Node::get_flux() {
	return flux;
}
;

inline double Node::get_elevation() {
	return elevation;
}
;

inline unsigned* Node::pass_key() {
	return key;
}

inline void Node::put_num_assoc_elem(int numin) {
	num_assoc_elem = numin;
}
;

inline int Node::get_num_assoc_elem() {
	return num_assoc_elem;
}
;

#endif
