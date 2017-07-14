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
 * $Id: element2.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifndef ELEMENT_H
#define ELEMENT_H
#include "boundary.h"
#include "hashtab.h"
#include "node.h"
#include "struct.h"

#include <fstream>
#include <iostream>

using namespace std;

template<typename T>
class ElemPtrList;
#include "../header/matrix.h"

//#define USE_FATHER

//! The Element class is a data structure designed to hold all the information need for an h (cell edge length) p (polynomial order) adaptive finite element.  Titan doesn't use p adaptation because it is a finite difference/volume code, hence many of the members are legacy from afeapi (adaptive finite element application programmers interface) which serves as the core of titan.  There is a seperate Discontinuous Galerkin Method (finite elements + finite volumes) version of titan and the polynomial information is not legacy there.  However in this version of Titan elements function simply as finite volume cells.

class Element {

	friend class HashTable;

	friend class DualElem;

	friend class ErrorElem;

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

	friend void unrefine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int myid,
	    void* NewFatherList);

	friend void unrefine_interp_neigh_update(HashTable* El_Table, HashTable* NodeTable, int nump,
	    int myid, void* OtherProcUpdate);

	friend void BSFC_combine_elements(int side, Element *EmTemp, HashTable *HT_Elem_Ptr,
	    HashTable *HT_Node_Ptr, int destination_proc);

	friend void Pack_element(void *sendel, ElemPack* elem, HashTable* HT_Node_Ptr,
	    int destination_proc);

	friend void destroy_element(void *r_element, HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
	    int target_pro, ELinkPtr* EL_head);

	friend void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
	    int myid, double* e_error, MatProps *matprops);

	friend void construct_el(Element* newelement, ElemPack* elem2, HashTable* HT_Node_Ptr, int myid,
	    double* e_error, MatProps *matprops);

	friend void uniform_refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int nump, int myid, void* RL,
	    TimeProps* timeprops_ptr);

	friend void adjust_node_info(MeshCTX* meshctx, PropCTX* propctx);

public:

	//! default constructor, does nothing except set stoppedflags=2, this should never be used
	Element() {
		father[0] = father[1] = 0; //initialize the father key to zero
		for (int i = 0; i < NUM_STATE_VARS; i++)
			state_vars[i] = -1;
		for (int i = 0; i < NUM_STATE_VARS; i++)
			Influx[i] = 0.;

		adapted = TOBEDELETED;
		refined = 1;
	};

	//! constructor that creates an original element when funky is read in
	Element(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[], int mat,
	    int *elm_loc_in, double pile_height, int myid, unsigned *opposite_brother);

	//! constructor that creates a son element from its father during refinement
	Element(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[], int gen,
	    int elm_loc_in[], int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
	    HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr);

	//! constructor that creates a father element from its four sons during unrefinement
	Element(Element *sons[], HashTable *NodeTable, HashTable *El_Table, MatProps *matprops_ptr);

	//! constructor that creates/restores a saved element during restart
	Element(FILE* fp, HashTable* NodeTable, MatProps* matprops_ptr, int myid);

	Element(gzFile& myfile, HashTable* NodeTable, MatProps* matprops_ptr, int myid);

	//! destructor that does nothing except delete boundary condition pointer
	virtual ~Element() {
	}
	;

	//! this member function saves a single element to a file with a single fwrite call, this allows the element to be recreated/restored upon restart of a simulation
	void save_elem(FILE* fp, FILE* fptxt); //for restart

	//! returns address of element (same as bubble node, node 8 out of 0->8) hashtable key
	unsigned* pass_key();

	//! returns the integer material flag for this element, needed for use of a material map which allows bedfriction to vary with physical position
	int get_material();

	//! returns the address of the first of 8 (nodes 0-7) node keys in an array, the node keys are used to access the nodes through the node hashtable
	unsigned* getNode();

	//! set the generation (number of times it's been refined -8<=gen<=+3) of this "element"/cell
	void put_gen(int);

	//! store the keys for the four son "elements" in the father element, used temporarily during refinement
	void putson(unsigned*);

	//! when a father element is refined into 4 son elements, the 4 son elements are "brothers" (they can be recombined into the father), this function stores the keys of all four brothers in one of them, it should be called 4 times one for each brother
	void putbrothers(unsigned*);

	//! this function returns the keys of an element's 4 brothers (an element is considered to be it's own brother) this is used during unrefinement to combine 4 brothers into their father element
	unsigned* get_brothers();

	//! this function stores the key "n" of neighbor "i" in the array of the 8 keys of the neighbor keys
	void putneighbor(unsigned *n, int i);

	//! this function stores the processor id "proc" of neighbor "i" in the 8 element array of neighbor processors, use this function instead of putassoc.
	void put_neigh_proc(int i, int proc);

	//! find and return what the key of this element's father element would be, very simple since the bubble node has the same key as the element, so all this function does is find which of its corner nodes will be the father element's bubble node, which it knows since it knows which_son it is.
	unsigned* getfather();

	//! store the father's key in the "father" variable, the "father's" key is zero until an element has been unrefined (and has not yet been deleted) it is only used in unrefinement. The getfather() member function computes the father key from "which_son" and it's nodes and is totally unrelated to the father variable.
	void put_father(unsigned fatherin[KEYLENGTH]);

	//! return the element keys of this element's 4 sons, used during refinement
	unsigned* getson();

	//! returns the element's error vector
	double* get_el_error();

	//! returns the array of keys for this element's 8 neighbors
	unsigned* get_neighbors();

	//! returns the array of processor ids for this element's 8 neighbors
	int* get_neigh_proc();

	//! returns this elements generation, that is how many times it's been refined -8<=generation<=+3, negative means courser than original mesh
	int get_gen();

	//! compare the FindNeigh key against the keys of this element's 8 neighbors to determine which if any neighbor FindNeigh is
	int which_neighbor(unsigned *FindNeigh);

	//! refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag().  refined can be permanently set to GHOST (defined in constant.h) or zero or temporarily set to 1 (with in the refinement and unrefinement routines), Keith believes it's not being unset (set from 1 to 0) when it should be after the refinement is done.  Keith believes the problem is located within H_adapt() or a function called from within it, recurse down.
	int get_refined_flag();

	//! magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell.  This allowed Keith to implement one time only immunity to unrefinement for recently refined (NEWSON) elements, which allowed him to protect a refined layer of buffer cells around piles.  Keith has partially replaced refined, get_refined_flag() and put_refined_flag() with adapted, get_adapted_flag() and put_adapted_flag(), but has left the if statements in the code responsible for refinement and unrefinement untouched because he encountered a bug, that he narrowed to within H_adapt() or a function called from within H_adapt(), recurse down, but has not pinpointed.  Keith believes the bug is related to the refined flag being inappropriately set to 1, or not unset to zero when it should be.
	int get_adapted_flag();

	//! call this function after this element's neighbor(s) have been refined, proc is processor id for neighbor[which_side+4]
	void change_neighbor(unsigned *newneighbs, int which_side, int proc, int reg);

	//! set this element's refined flag to i, can set it to normal (hasn't just been refined and isn't a ghost cell), "temporarily" set to "refined" (has just been refined so don't refine again), or say that it's a GHOST cell, see constant.h, (which means you don't update it, instead you get new values from the processor that owns it and you don't refine it.) refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag().
	void put_refined_flag(int i);

	//! refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag(). The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell. These values are defined in constant.h.  The NEWSON value has allowed Keith to provide one time only immunity from unrefinement to recently refined elements, after which the "adapted" flag is resent to NOTRECADAPTED.
	void put_adapted_flag(int new_adapted_status);

	//! this function returns an array holding the generation of all 8 of this element's neighbors
	int* get_neigh_gen();

	//! this function sets the ith neighbor's generation to "gen"
	void put_neigh_gen(int i, int gen);

	//! this function sets the which_son flag when a father element is refined into its 4 sons, the which_son flag tells the portion of the father element that this element is physically located in
	void put_which_son(int);

	//! returns the which_son flag, which tells the portion of the father element that this element is physically located in
	int get_which_son();

	//! this function calculates the which_son flag for the original (or restored in case of a restart) element.  It also calculates which son of the grandfather element the father is durring unrefinement.
	void calc_which_son();

	//! this function sets the new or old flag, it is initialized in htflush.C and reset during repartitioning (repartition_BSFC.C and BSFC_update_and_send_elements.C)
	void put_new_old(int i);

	//! this function returns the vlaue of the new_old flag which is used during mesh adaptation and repartitioning
	int get_new_old();

	//! this function is called during repartitioning when one of an element's neighbors is sent to another processor
	void change_neighbor_process(int which, int newp);

	//! the function returns the vector of element "error", element error is used to say when a function should be refined
	double* get_el_err();

	//! this function returns the Load Balancing weight of an element which is used in repartitioning
	double get_lb_weight();

	//! this function stores an element's load balancing weight, which is used during repartitioning
	void put_lb_weight(double dd_in);

	//! this function returns the load balancing key, which is used during repartitioning
	unsigned* get_lb_key();

	//! this function sets the load balancing key, which is used during repartitioning
	void put_lb_key(unsigned* in_key);

	//! this function copies the elmenent key to the load balancing key
	void copy_key_to_lb_key();

	//! this function sets the process(or) id of an element, it says which processor owns the element.
	void put_myprocess(int in_proc);

	//! this function returns the process(or) id of an element, it says which processor owns the element
	int get_myprocess();

	//! this function returns the opposite_brother_flag, I (Keith) am not entirely sure what this flag is for, but I know that it is used in repartioning, see BSFC_combine_elements, I think it says if an element has an opposite brother, that is, can it be combined with it's brothers to form their father
	int get_opposite_brother_flag();

	//! this function computes searches for an element's brother, i.e. the brother (son of the same father) that is located diagonally from it, to get the brother information requires that atleast one of this element's neighboring brothers is on this process in order to get information onthe brother that is not a neighbor
	void find_opposite_brother(HashTable*);

	/* geoflow functions */

	//! this function initializes pileheight, momentums and shortspeed (also known as the L'Hosptial speed see calc_shortspeed for an explanation),this function is called in init_piles.C
	void put_height_mom(double pile_height, double vfract, double xmom, double ymom);

	//! this function assigns a specified value to the pileheight and zeros to the momentums and shortspeed
	void put_height(double pile_height);

	//! this function returns the vector of state variables
	double* get_state_vars();

	//! this function returns the vector of adjoint
	double* get_adjoint();

	//! this function returns the vector of prev_adjoint
	double* get_prev_adjoint();

	// ! updates prev_adjoint before computing adjoint
	void update_adjoint();

	//! this function returns the vector of x and y derivatives of state variables, all the x derivatives come first as a group followed by the y derivatives as a group
	double* get_d_state_vars();

	//! this function returns the x and y slopes of the terrain elevation
	double* get_zeta();

	//! this function returns the length of an element in the x and y directions
	double* get_dx();

	//! this function computes which side of the element is facing the positive x direction
	void find_positive_x_side(HashTable*);

	//! this function returns which side of the element is facing the positive x direction
	int get_positive_x_side();

	//! this function computes the x and y derivatives of the state variables
	void get_slopes(HashTable*, HashTable*, double);
	virtual void get_slopes_prev(HashTable*, HashTable*, double);

	//! this function returns a vector containing the previous state variables, previous mean beginning of timestep before the finite difference predictor halfstep
	double* get_prev_state_vars();

	//! updates prev_states variables to current states, for first order-calculations
	void update_prev_state_vars();

	//! this function calculates the lengths of the element in the (global) x and y directions
	void calculate_dx(HashTable* NodeTable);

	void set_mins(HashTable* NodeTable);

	//! this function assigns the element's coordinates to be its bubble node's coordinates
	void insert_coord(HashTable* NodeTable);

	//! this function, based on the dir flag, chooses between calling xdirflux and ydirflux, which respectively, calculate either the x or y direction analytical cell center fluxes (or the fluxes at the the boundary if 2nd order flux option is checked on the gui). Keith wrote this.
	void zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int order_flag,
	    int dir, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], Element *EmNeigh,
	    double dt);

	void zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int order_flag,
	    int dir, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], Element *EmNeigh,
	    double dt, ResFlag resflag);

	//! this function calculates the analytical cell center (or cell boundary if 2nd order flux flag is checked on the gui) x direction fluxes. Keith wrote this
	void xdirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
	    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], ResFlag resflag);

	//! this function calculates the analytical cell center (or cell boundary if 2nd order flux flag is checked on the gui) y direction fluxes. Keith wrote this
	void ydirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
	    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], ResFlag resflag);

	//! this function calculates the analytical cell center (or cell boundary if 2nd order flux flag is checked on the gui) x direction fluxes. Keith wrote this
	virtual void xdirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
	    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]);

	//! this function calculates the analytical cell center (or cell boundary if 2nd order flux flag is checked on the gui) y direction fluxes. Keith wrote this
	virtual void ydirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
	    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]);

	void calc_flux(HashTable* El_Table, HashTable* NodeTable, vector<Element*>* elem_list, int myid,
	    int side);

	void calc_fluxes(HashTable* El_Table, HashTable* NodeTable, vector<Element*>* x_elem_list,
	    vector<Element*>* y_elem_list, int myid);

	virtual void boundary_flux(HashTable* El_Table, HashTable* NodeTable, int myid, int side);

	void dual_zdirflux(int dir, Mat3x3& hfv, Mat3x3& flux_jac, Mat3x3& s_jac);
	void dual_ydirflux(Mat3x3& hfv, Mat3x3& flux_jac, Mat3x3& s_jac);
	void dual_xdirflux(Mat3x3& hfv, Mat3x3& flux_jac, Mat3x3& s_jac);

	//! this function (indirectly) calculates the fluxes that will be used to perform the finite volume
	//corrector step and stores them in element edge nodes, indirectly because it calls other functions
	//to calculate the analytical fluxes and then calls another function to compute the riemann fluxes
	//from the analytical fluxes. Talk to me (Keith) before you modify this, as I am fairly certain that
	//it is now completely bug free and parts of it can be slightly confusing.
	template<class T>
	void calc_edge_states(HashTable* El_Table, HashTable* NodeTable, vector<T*>* x_elem_list,
	    vector<T*>* y_elem_list, MatProps* matprops_ptr, int myid, double dt, int* order_flag,
	    double *outflow);

	/*! this function computes the velocity */
	void eval_velocity(double xoffset, double yoffset, double Vel[]);

	//! this function returns the already calculated value(s) of k active passive, which comes from using th Coulomb friction model of granular flows (this is problem specific to titan and thus does not appear in the standard afeapi code)
	double* get_kactxy();

	//! returns the already computed gravity vector in local coordinates, the local z direction is normal to the terrain surface and the projection of the local x and y components into the horizontal plane are aligned with global x (UTM E) and y (UTM N) directions.
	double* get_gravity();

	//! this function is titan legacy code it is defined in Element2.C but is not called anywhere
	int determine_refinement(double);

	//! this function returns the precomputed elevation
	double get_elevation();

	//! this function returns the precomputed derivatives of the z component of gravity, this is a purely terrain geometry dependant derivative, that is little diffent than curvature
	double* get_d_gravity();

	//! this function returns the precomputed local terrain curvature.  Curvature itself is the inverse of radius of curvature.  The exact value of curvature  is the spatial second derivative of the normal coordinate of the surface along directions tangent to the surface at that point (local x and y).  However I believe that Laercio Namikawa implemented it approximately, i.e. as the global x and y second derivatives of terrain elevation.
	double* get_curvature();

	//! this function is called in element_weight.C, it is used in computing the load balancing weight
	void calc_flux_balance(HashTable *NodeTable);

	//! this function calculates topographic data, it calls GIS commands to compute elevation, slopes, and curvatures from the GIS map and scales them appropriately
	void calc_topo_data(MatProps *matprops_ptr);

	//! this function calculates the (global) x and y derivatives of the local z component of gravity as an approximation of the local derivatives, it wouldn't be that difficult to correct incorporating the terrain slopes in the calculation it is calculated in the creation of a father element, after mesh refinement and, during a restart.
	void calc_d_gravity(HashTable *El_Table);

	//! this function calculates the gravity vector in local coordinates
	void calc_gravity_vector(MatProps *matprops_ptr);

	//! this function is defined in unrefine.C, it is also called in that file, it finds this element's brothers
	template<class T>
	int find_brothers(HashTable* El_Table, HashTable* NodeTable, double target, int myid,
	    MatProps* matprops_ptr, void* NewFatherList, void* OtherProcUpdate);

	template<class T>
	int dual_find_brothers(HashTable* El_Table, HashTable* NodeTable, double target, int myid,
	    MatProps* matprops_ptr, void* NewFatherList, void* OtherProcUpdate);

	void for_link_temp();

	void dual_find_brothers(HashTable* El_Table, HashTable* NodeTable, int myid,
	    MatProps* matprops_ptr, void* NewFatherList, void* OtherProcUpdate, int rescomp);

	//! this function is defined in unrefine.C, it is also called in that file and no where else, it prevents refinement when one or more of the brothers does not belong to this processor
	int check_unrefinement(HashTable *El_Table, double target);
	int check_dual_unrefinement(HashTable* El_Table, double target);

	//! this function updates this elements neighbor info when one of its neighbors has been unrefined
	void change_neigh_info(unsigned *fth_key, unsigned *ng_key, int neworder, int ng_gen,
	    int fth_proc);

	//! this function returns the elm_loc variable, which is used in unrefinement beyond the initial coarse grid
	int* get_elm_loc();

	//! this function sets the elm_loc variable, which is used in unrefinement beyond the initial coarse grid
	void put_elm_loc(int* int_in);

	//! this function returns the precomputed and scaled coordinates of this element (which would be the same as its bubble node's coordinates)
	double* get_coord();

	//! this function stores the coordinates of this element (which would be the same as its bubble node's coordinates)
	void put_coord(double* coord_in);

	//! interface to change value of earth-pressure cofficients
	void put_kactxy(double kap[]) {
		kactxy[0] = kap[0];
		kactxy[1] = kap[1];
	}

	double* get_tanbedfrict();

	//! this function zeros the extrusion (out of the ground) fluxes in this element
	void zero_influx();

	//! this function returns the stored value of the extrusion (out of the ground) fluxes in this element
	double *get_influx();

	//! this function calculates the extrusion (out of the ground) fluxes for this elements
	void calc_flux(HashTable *NodeTable, FluxProps *fluxprops, TimeProps *timeprops);

	//! this function returns 2 if this element contains pileheight>=contour_height and has a neighbor who contains pileheight<contour_height.  It returns 1 if this element contains pileheight<contour_height and has a neighbor who contains pileheight>=contour_height.  It returns 0 otherwise. The intended use if if(EmTemp->if_pile_boundary(ElemTable,contour_height)) but I (Keith) added the distinction bewteen 1 and 2 to allow future developers to distinguish between the inside and outside of a pileheight contour line, as this functionality could be useful in the future.
	int if_pile_boundary(HashTable *ElemTable, double contour_height);

	//! this function returns 2 if this element has Influx[0]>0 and has a neighbor who has Influx[0]<=0.  It returns 1 if this element has Influx[0]==0 and has a neighbor who has Influx[0]!=0.  It returns -1 if this element has Influx[0]<0 and a neighbor with Influx[0]>=0. It returns 0 otherwise. Influx[0] is a pileheight per unit time source term.  Currently Influx[0] is restricted to be non-negative (a source or no source with sinks not allowed), but I (Keith) have added the extra functionality because it may be useful at a future date. The intended use if if(EmTemp->if_source_boundary(ElemTable)), but the distinction between 1 and 2 allows futuredevelopers to distinguish between the strictly inside and strictly outside of an area with a flux source term.
	int if_source_boundary(HashTable *ElemTable);

	//! the buffer layer is a layer of refined cells on the outside of the pile, i.e. ((pileheight<contour_height)&&(Influx[0]==0)) and adjacent to the pile.  It is "N" elements wide, and the "N" element width is increased one element at a time.  This function returns 2 if this element a member of the innermost boundary of the buffer and does not need to be adapted.  It returns 1 if this elment needs to be refined and some of its sons will be members of the innermost boundary of the buffer layer
	int if_first_buffer_boundary(HashTable *ElemTable, double contour_height);

	//! the buffer layer is a layer of refined cells on the outside of the pile, i.e. ((pileheight<contour_height)&&(Influx[0]==0)) and adjacent to the pile.  It is "N" elements wide, and the "N" element width is increased one element at a time.  This function returns 2 if this element a member of the boundary of the buffer that is one element wider than the current buffer and does not need to be adapted.  It returns 1 if this elment needs to be refined and some of its sons will be in the next buffer boundary
	int if_next_buffer_boundary(HashTable *ElemTable, HashTable *NodeTable, double contour_height);

	//! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
	int get_ithelem() const;

	//! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
	void put_ithelem(int i);

	//! sgn of double
	double sgn(double a) {
		return (a < 0 ? -1. : 1.);
	}

	void for_link_temp1();

	//! this function returns the neighbor of element considering its position, sides:xp=0, yp=1, xm=2, ym=3, xpyp=4, xpym=5, xmym=6, xpym=7
	template<class T>
	T* get_side_neighbor(HashTable *El_Table, int side);

	void gen_my_sons_key(HashTable* El_Table, unsigned son_key[4][KEYLENGTH]);

	void gen_my_sons_key(HashTable* El_Table, unsigned* son_key);

	int check_state(SolRec* solrec, HashTable* El_Table, int iter);

	virtual void write_elem_info(HashTable* NodeTable, char* filename, int iter, double dt);

	void calc_edge_states(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int myid,
	    double dt, int* order_flag, double *outflow, ResFlag lresflag, ResFlag rresflag);

	void write_elem(gzFile& myfile);

protected:
	//! myprocess is id of the process(or) that owns this element
	int myprocess;

	//! generation is how many times this element has been refined, currently -8<=generation<=3, a negative generation number means it has been unrefined beyond the orignal coarse mesh, a positive generation number means it has been refined (is smaller than the original element size)
	int generation;

	//! opposite_brother_flag indicate if we have the correct key for the non-neighbor brother (0:= don't have info, 1:= have info)
	int opposite_brother_flag;

	//! the material flag indicates which material should be used to set this element's bed friction, this is for when a GIS material map, specifying different materials in different spatial regions of the map, the GIS material map is a non standard grass map format that Laercio Namikawa developed, it's stored in the "cats" folder under a grass mapset directory
	int material;/*! ! ! THE MAT. FLAG ! ! !*/

	//! this is the load-balancing weight
	double lb_weight;

	//! this is the key for load-balancing, if there is no constrained node, it is the element key, otherwise it is a construct of the element "bunch", keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
	unsigned lb_key[KEYLENGTH];

	//! this is the element key, which has the same value as the key of the element's bubble node, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
	unsigned key[KEYLENGTH];

	//! this array holds the first 8 (0->7) of this element's nodes' keys, the n9th (8 out of 0->8) node is the bubble node it's key is not stored separately since it has the same key as the element, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
	unsigned node_key[8][KEYLENGTH];

	//! this array holds the keys of this element's 8 neighbors (2 neigbors to a side if the neighbor is more refined than this element, otherwise the two neighbor keys for that side are identical in value), having 8 neighbors is an outcome of the 1 irregularity refinement rule, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
	unsigned neighbor[8][KEYLENGTH];

	//! the key of the father it is assigned in the refine() and unrefine_elements() functions
	unsigned father[KEYLENGTH];

	//! this array holds the keys of this element's 4 sons, it is only used temporarily in the refinement process before the father (this element) is deleted, there's was an old comment "garantee ccw" associated with this variable, I don't know what it means, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
	unsigned son[4][KEYLENGTH];

	//! this array holds the process(or) id of this element's 8 neighbors, there can be 8 neighbors because of the 1 irregularity rule.  neigh_proc[4:7] != -2 only if it has 2 neighbors on that side, a value of -1 for neigh_proc means that this edge is a boundary of the computational domain.
	int neigh_proc[8];

	//! neigh_gen is an array that holds the "generation" (how refined it is) of this element's 8 neighbors, there can-be/are 2 neighbors to a side because of the 1 irregularity rule
	int neigh_gen[8];

	//! this holds the "error" in the element's solution, which is useful in determining refinement, this may actually be afeapi legacy
	double el_error[EQUATIONS];

	//! refined is a flag that usually has the value 0, but will be 1 if the element has been refined this iteration (used to enforce the 1 irregularity rule), or have the value "GHOST" if it is a ghost cell, refined and ghost cells are not updated, see constant.h for the value of GHOST
	int refined;

	//! The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell. This allowed Keith to implement one time only immunity to unrefinement for recently refined (NEWSON) elements, which allowed him to protect a refined layer of buffer cells around piles.  Keith has partially replaced refined, get_refined_flag() and put_refined_flag() with adapted, get_adapted_flag() and put_adapted_flag(), but has left the if statements in the code responsible for refinement and unrefinement untouched because he encountered a bug, that he narrowed to within H_adapt() or a function called from within H_adapt(), recurse down, but has not pinpointed.  Keith believes the bug is related to the refined flag being inappropriately set to 1, or not unset to zero when it should be.
	int adapted;

	//! which_son holds the value of which son this element is, which of the 4 squares that makes up the father elements square.
	int which_son;

	//! the new_old flag is used in mesh adaptation and repartitioning
	int new_old;

	//! this array holds the keys of this element's 4 brothers (an element is considered to be it's own brother), this information is used during mesh unrefinement (combining the 4 brothers to make their father), keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
	unsigned brothers[4][KEYLENGTH];

	//! coord holds the coordinates of the elements cell center, these are the same as the coordinates of the element's bubble node's
	double coord[DIMENSION];

	//! elm_loc is used in unrefining beyond the original coarse mesh
	int elm_loc[2];

	/* variables for hyperbolic geoflow problem */

	//! state_vars is an array that holds the current state variables: h, hVx, and hVy
	double state_vars[NUM_STATE_VARS];

	//! these are the values of the state variables from before the predictor step
	double prev_state_vars[NUM_STATE_VARS];

	//! these are the spatial (x and y) derivatives of the state variables: (dh/dx, dhVx/dx, dhVy/dx, dh/dy, dhVx/dy, dhVy/dy)
	double d_state_vars[NUM_STATE_VARS * DIMENSION];

	//! length of the element in the global x and y directions: dx and dy
	double dx[DIMENSION];

	//! for structured grid, tells which side is the positive x direction
	int positive_x_side;

	//! k active/passive in the x and y directions, k active/passive is part of the coulomb friction model for Granular Flows
	double kactxy[DIMENSION];

	//this is tan bed frict angle but it changes w.r.t stop criteria
	double tan_bed_frict;

	//! terrain elevation at this elements center/bubble node
	double elevation;

	//! terrain slope in the global x and y directions
	double zeta[DIMENSION];

	//! Curvature itself is the inverse of radius of curvature.  The exact value of curvature is the spatial second derivative of the normal coordinate of the surface along directions tangent to the surface at that point (local x and y).  However I believe that Laercio Namikawa implemented it approximately, i.e. as the global x and y second derivatives of terrain elevation.
	double curvature[DIMENSION];

	//! the gravity vector in local x,y,z coordinates (z is normal to the terrain surface, the projections of the x and y local directions onto a horizontal plane are aligned with the global x and y directions)
	double gravity[3];

	//! the spatial (x and y) derivatives of the local z component of the gravity vector
	double d_gravity[DIMENSION];

	//! extrusion flux rate for this timestep for this element, used when having material flow out of the ground, a volume per unit area influx rate source term
	double Influx[NUM_STATE_VARS];

	//! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
	int ithelem;
};

inline int Element::get_ithelem() const {
	return ithelem;
}
;

inline void Element::put_ithelem(int i) {
	ithelem = i;
}
;

inline unsigned* Element::pass_key() {
	return key;
}
;

inline int Element::get_material() {
	return material;
}
;

inline unsigned* Element::get_brothers() {
	return &brothers[0][0];
}
;

inline double Element::get_lb_weight() {
	return lb_weight;
}
;

inline void Element::put_lb_weight(double dd_in) {
	lb_weight = dd_in;
}
;

inline unsigned* Element::get_lb_key() {
	return lb_key;
}
;

inline void Element::put_myprocess(int in_proc) {
	myprocess = in_proc;
}
;

inline int Element::get_myprocess() {
	return myprocess;
}
;

inline int Element::get_opposite_brother_flag() {
	return opposite_brother_flag;
}
;

inline void Element::put_height_mom(double pile_height, double volf, double xmom, double ymom) {
	prev_state_vars[0] = state_vars[0] = pile_height;
	prev_state_vars[1] = state_vars[1] = xmom;
	prev_state_vars[2] = state_vars[2] = ymom;
}

inline void Element::put_height(double pile_height) {
	put_height_mom(pile_height, 1., 0., 0.);
	return;
}

inline double* Element::get_state_vars() {
	return state_vars;
}
;

inline double* Element::get_d_state_vars() {
	return d_state_vars;
}
;

inline double* Element::get_dx() {
	return dx;
}
;

inline int Element::get_positive_x_side() {
	return positive_x_side;
}
;

inline double* Element::get_prev_state_vars() {
	return prev_state_vars;
}
;

inline void Element::update_prev_state_vars() {
	for (int i = 0; i < NUM_STATE_VARS; i++)
		prev_state_vars[i] = state_vars[i];
}

inline double* Element::get_kactxy() {
	return kactxy;
}
;

inline double* Element::get_tanbedfrict(){
	return &tan_bed_frict;
}
;

inline double* Element::get_gravity() {
	return gravity;
}
;

inline double Element::get_elevation() {
	return elevation;
}
;

inline double* Element::get_d_gravity() {
	return d_gravity;
}
;

inline double* Element::get_curvature() {
	return curvature;
}
;

inline int* Element::get_elm_loc() {
	return elm_loc;
}
;

inline void Element::put_elm_loc(int* int_in) {
	elm_loc[0] = int_in[0];
	elm_loc[1] = int_in[1];
}
;

inline double* Element::get_coord() {
	return coord;
}
;

inline unsigned* Element::getNode() {
	return &(node_key[0][0]);
}

inline void Element::put_gen(int g) {
	generation = g;
}

inline void Element::put_neigh_proc(int i, int proc) {
	neigh_proc[i] = proc;
}

inline unsigned* Element::getson() {
	return &(son[0][0]);
}

inline void Element::putson(unsigned* s) {

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < KEYLENGTH; j++)
			son[i][j] = *(s + i * KEYLENGTH + j);

	refined = 1;
	adapted = OLDFATHER;
}

inline void Element::putbrothers(unsigned* s) {

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < KEYLENGTH; j++)
			brothers[i][j] = *(s + i * KEYLENGTH + j);

	return;
}

inline void Element::putneighbor(unsigned* n, int i) {
	int j;
	for (j = 0; j < KEYLENGTH; j++)
		neighbor[i][j] = *(n + j);
}

inline double* Element::get_el_error() {
	return el_error;
}

inline unsigned* Element::get_neighbors() {
	return &neighbor[0][0];
}

inline int* Element::get_neigh_proc() {
	return neigh_proc;
}

inline int Element::get_gen() {
	return generation;
}

inline int Element::get_refined_flag() {
	return refined;
}

inline void Element::put_refined_flag(int i) {
	refined = i;
}

inline int Element::get_adapted_flag() {
	return adapted;
}

inline void Element::put_adapted_flag(int new_adapted_status) {
	adapted = new_adapted_status;
}

inline int* Element::get_neigh_gen() {
	return neigh_gen;
}

inline void Element::put_neigh_gen(int i, int gen) {
	neigh_gen[i] = gen;
}

inline void Element::put_which_son(int i) {
	which_son = i;
}

inline int Element::get_which_son() {
	return which_son;
}

inline int Element::get_new_old() {
	return new_old;
}

inline void Element::put_new_old(int i) {
	new_old = i;
}

inline void Element::change_neighbor_process(int which, int newp) {
	neigh_proc[which] = newp;
}

inline void Element::zero_influx() {
	for (int i = 0; i < NUM_STATE_VARS; i++)
		Influx[i] = 0.;
}

inline double* Element::get_influx() {
	return Influx;
}

inline void Element::put_father(unsigned *fatherin) {
	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = fatherin[ikey];
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

//! The ElemPtrList class is basically just a "smart array" of pointers to Elements, by smart I mean it keeps track of its size and number of Elements in the list and expands/reallocates itself whenever you add an element ptr to the list when you've run out of space, it also keeps a record of the index of the first "new" element pointer you've added in the current series, which is useful for the intended purpose... ElemList was designed for use in refinement and unrefinement to replace fixed sized arrays (length=297200) of pointers to Elements.  The reason for this upgrade was it was causing valgrind to issue all kinds of warnings about the "client switching stacks" and "invalid write/read of size blah blah blah" because the stacksize was too large.  My 20061121 rewrite of hadapt.C and unrefine.C to make them "fast" caused this problem because I added a second (large) fixed sized array to both of them so I could reduce the number of hashtable scans by only revisiting the "new" additions to the array of pointers of Elements. --Keith wrote this on 20061124, i.e. the day after Thanksgiving, and I'm very thankful for having the inspiration to figure out the cause of valgrid warning
template<typename T>
class ElemPtrList {
public:

	/*
	 ElemPtrList();
	 ElemPtrList(int initial_size);
	 ~ElemPtrList();
	 void add(Element* EmTemp);
	 */

//! this constructor allocates space for an array of the default initial size (1024 Element pointers), the size_increment equal the initial size
	ElemPtrList() {
		init(1024);
		return;
	}
	;

//! this constructor allocates space for user specified initial-size, the size_increment equals the initial size.
	ElemPtrList(int initial_size) {
		if (initial_size == 0)
			initial_size = 1024;
		init(initial_size);
		return;
	}
	;

//! the destructor frees the list space so the programmer never has to worry about it
	~ElemPtrList() {
		//printf("list_space=%d, num_elem=%d, inewstart=%d\n",list_space,num_elem,inewstart);
		free(list);
		return;
	}

//! add an element pointer to the list, it will increase the size of the list by the size_increment if necessary, the size_increment is the initial size.
	void add(T* EmTemp);

//! returns the ith Element pointer stored in the list
	T* get(int i);

//! returns the key of the ith Element whose pointer is stored in the list
	unsigned* get_key(int i);

//! returns the number of elements in the list.
	int get_num_elem();

//! marks the "starting position" of a new "series" of Element pointers stored
	void set_inewstart(int inewstart_in);

//! returns the "starting position" of the new "series" of Element pointers stored in the list
	int get_inewstart();

//! zeros the list (does not deallocate space)
	void trashlist();

private:
//! actually creates the list, is called by the constructors
//void      init(int initial_size);
	void init(int initial_size) {
		list_space = size_increment = initial_size;
		num_elem = inewstart = 0;
		list = (T **) malloc(list_space * sizeof(T*));
		for (int i = 0; i < list_space; i++)
			list[i] = NULL;
	}
	;

//! number of elements whose pointers are stored in the list
	int num_elem;

//! when the list runs out of space increase it by this much
	int size_increment;

//! the current size of the list (the ammount of memory allocated not number of nonzero entries)
	int list_space;

//! a convienient way to mark the start of a new series in the list
	int inewstart;

//! the "array" holding the list of element pointers, space is allocated by the constructor, increased automatically whenever needed, and freed automatically by the destructor
	T** list;
};

template<typename T>
inline T* ElemPtrList<T>::get(int i) {
	return (((i >= 0) && (i < num_elem)) ? list[i] : NULL);
}
;

template<typename T>
inline unsigned* ElemPtrList<T>::get_key(int i) {
	return (((i >= 0) && (i < num_elem)) ? list[i]->pass_key() : NULL);
}
;

template<typename T>
inline int ElemPtrList<T>::get_inewstart() {
	return inewstart;
}
;

template<typename T>
inline void ElemPtrList<T>::set_inewstart(int inewstart_in) {
	inewstart = inewstart_in;
	return;
}
;

template<typename T>
inline int ElemPtrList<T>::get_num_elem() {
	return num_elem;
}
;

template<typename T>
inline void ElemPtrList<T>::trashlist() {
	for (int i = 0; i < num_elem; i++)
		list[i] = NULL;
	num_elem = inewstart = 0;
	return;
}
;

template<typename T>
inline void ElemPtrList<T>::add(T* EmTemp) {
	if (num_elem == list_space - 1) {
		list_space += size_increment;
		list = (T **) realloc(list, list_space * sizeof(T *));
	}

	list[num_elem] = EmTemp;
	num_elem++;
	return;
}
;

#endif
