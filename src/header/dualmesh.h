/*
 * dualmesh.h
 *
 *  Created on: May 25, 2015
 *      Author: haghakha
 */

#ifndef DUALMESH_H
#define DUALMESH_H

#include "hashtab.h"
#include "element2.h"
#include "jacobian.h"
#include <hdf5.h>
#include "GMFG_hdfapi.h"
#include "geoflow.h"

class ErrorElem;
class Node;

class Node_minimal{

	friend class Node;

private:
	int id, info;
	double coord[2],elevation;
	unsigned key[2];

public:
	Node_minimal(Node *node);

};

class Elem_minimal{

	friend class Element;

private:
	int myprocess, generation;
	int neigh_proc[8],neigh_gen[8];
	int refined,adapted;
	int elm_loc[2],material;
	int opposite_brother_flag;

	unsigned key[DIMENSION];
	unsigned node_key[8][DIMENSION];
	unsigned neighbor[8][DIMENSION];
	unsigned son[4][DIMENSION];
	unsigned brothers[4][DIMENSION];

	double coord[DIMENSION];
	double state_vars[NUM_STATE_VARS];
	double prev_state_vars[NUM_STATE_VARS];

public:
	Elem_minimal(Element *elem);

};

class Table_minimal{
	friend class HashTable;

private:
	unsigned MinKey[2];
	unsigned MaxKey[2];
	double doublekeyrange[2];
	double hashconstant;
	double Xrange[2];
	double Yrange[2];
	double invdxrange, invdyrange;
	int NBUCKETS, PRIME;

public:
	Table_minimal(HashTable* table);

};

class Snapshot{
private:

	int iter;
	int num_nodes, num_elems;
	double time;

	vector<Node_minimal> node_vec;
	vector<Elem_minimal> elem_vec;
	Table_minimal node_tab, elem_tab;

public:
	Snapshot(const MeshCTX& meshctx, const PropCTX& propctx);
	vector<Node_minimal>* get_node_vector(){
		return &node_vec;
	};
	vector<Elem_minimal>* get_elem_vector(){
		return &elem_vec;
	};
	Table_minimal* get_node_table_minimal(){
		return &node_tab;
	};
	Table_minimal* get_elem_table_minimal(){
		return &elem_tab;
	};
	int get_iter(){
		return iter;
	};
	double get_time(){
		return time;
	};

	void adjust_timeprops(TimeProps* timeprops,const int final_iter){
		double inv_t_sc = 1./timeprops->TIME_SCALE;
		timeprops->iter=iter;
		timeprops->time=time;
		timeprops->isave = (int) (time / (timeprops->timesave * inv_t_sc));
		timeprops->ioutput = (int) (time / (timeprops->timeoutput * inv_t_sc));
		timeprops->ndnextsave = ((timeprops->isave + 1) * timeprops->timesave) * inv_t_sc;
		timeprops->ndnextoutput = ((timeprops->ioutput + 1) * timeprops->timeoutput) * inv_t_sc;;
		timeprops->maxiter = final_iter;
	}

	~Snapshot()	{
	};
};

class SolRec: public HashTable {

private:

	// this integer shows the first time step that its solution is available in SolRec
	int first_solution_time_step;

	// this integer shows the last time step that its solution is available in SolRec
	int last_solution_time_step;

	const int range;

public:

	// constructor
	SolRec(double *doublekeyrangein, int size, int prime, double XR[], double YR[]);

	SolRec(gzFile& myfile);

	// this function records the solution of last time step
	void record_solution(MeshCTX* meshctx, PropCTX* propctx);

	// this function writes the recorded solution to files, and for each time step separately
	void wrtie_sol_to_disk(int myid);

	void wrtie_sol_to_disk_hdf5(int myid);

	// this function reads the recorded solution from files, and for each time step separately and store it in SolRec
	void read_sol_from_disk(int myid, int iter);

	void read_sol_from_disk_hdf5(int myid, int iter);

	// this functions deletes the SolRec and deallocate the allocated memory
	void delete_jacobians_after_writes();

	// this function changes first_solution_time_step to the next time step
	void update_first_sol_time(int iter);

	// this function changes last_solution_time_step to the new time step
	void update_last_sol_time(int iter);

	int get_first_solution();

	int get_last_solution();

	int get_range(){return range;}

	void load_new_set_of_solution(int myid);

	void free_all_available_sol();

	int data_range();

	int write_sol();

	int read_sol();

	using HashTable::lookup;

	Solution* lookup(unsigned* key, int iter, int erase_map=0);

	void write_table(gzFile& myfile,run_mode = NORMAL);

    virtual	~SolRec();

};

class DualElem: public Element {

public:

	DualElem(Element* element);

	//used for refinement
	DualElem(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[], int gen,
	    int elm_loc_in[], int gen_neigh[], int mat, DualElem *fthTemp, double *coord_in,
	    HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr);

	//used for unrefinement
	DualElem(DualElem* sons[], HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr);

	DualElem(DualElemPack* elem2, HashTable* HT_Node_Ptr, MatProps* matprops_ptr, int myid);

	DualElem(gzFile& myfile, HashTable* NodeTable, MatProps* matprops_ptr, int myid);

	void update(DualElemPack* elem2, HashTable* HT_Node_Ptr, int myid);

	void Pack_jacobian(JacPack* jac, int neigh);

	void put_jacpack(JacPack* jac);

	void Pack_element(DualElemPack* elem, HashTable* HT_Node_Ptr, int destination_proc);

	void get_slopes_prev(HashTable* El_Table, HashTable* NodeTable, double gamma);

	void calc_flux(HashTable* El_Table, HashTable* NodeTable, vector<DualElem*>* elem_list, int myid,
	    int side);

	void calc_fluxes(HashTable* El_Table, HashTable* NodeTable, vector<DualElem*>* x_elem_list,
	    vector<DualElem*>* y_elem_list, int myid);

	void boundary_flux(HashTable* El_Table, HashTable* NodeTable, const int myid, const int side);

	void set_jacobian(int neigh_num, double elemjacob[NUM_STATE_VARS], int state_vars_num,
	    const double incr);

	void dual_check_refine_unrefine(SolRec* solrec, HashTable* El_Table, int iter,
	    ElemPtrList<DualElem>* refinelist, ElemPtrList<DualElem>* unrefinelist);

	void update_state(SolRec* solrec, HashTable* El_Table, int iter);

	void set_jacobian();

	double* get_prev_adjoint();

	double* get_adjoint();

	void update_adjoint();

	double* get_func_sens();

	FluxJac& get_flx_jac_cont();

	Matrix<double, 2, 5>& get_hslope_sens();

	Vec_Mat<9>& get_jacobian();

	void print_jacobian(int iter);

	void set_jacobianMat_zero(int jacmatind);

	void calc_func_sens(const void * ctx);

	ErrorElem** get_son_addresses();

	void dual_check_refine_unrefine_repartition(SolRec* solrec, HashTable* El_Table, int iter,
	    ElemPtrList<DualElem>* refinelist, ElemPtrList<DualElem>* unrefinelist,
	    vector<TRANSKEY>& trans_keys_vec, vector<int>& trans_keys_status,
	    vector<DualElem*>& repart_list, double* allKeyRange);

	void write_elem(gzFile& myfile);

private:

	//! adjoint vector
	double adjoint[NUM_STATE_VARS];

	//! prev. time step adjoint vector
	double prev_adjoint[NUM_STATE_VARS];

	//! double jacobianMat [8][3][3],self[3][3],neigh1[3][3],neigh2[3][3],neigh3[3][3],neigh4[3][3]
	Vec_Mat<9> jacobianMat;

	double func_sens[NUM_STATE_VARS];

	// this array holds the jacobian of fluxes,and its indexing is [side:x=0,y=1][direction:neg=0,pos=1]
	//Matrix<Matrix> FluxJac (2,2);
	FluxJac flx_jac_cont;

	//2: 0->for x direction, 1-> for y direction
	//5: 0-> element itself, 1->first neighbor in left(down), 2->first neighbor in right(up),
	//3->second neighbor in left(down), 4->second neighbor in right(up),
	Matrix<double, 2, 5> hslope_sens;

	// this is to keep the address of his four sons in Error grid
	// keeping these addresses let us to avoid calling lookup function
	// for many times
	// actually we just need to keep the address for 4 sons and 16 is just for the cases
	// that the four element in error grid has been refined but still the element in dual grid
	//
	ErrorElem* mysons[4];
};

inline double* DualElem::get_adjoint() {
	return adjoint;
}
;

inline double* DualElem::get_prev_adjoint() {
	return prev_adjoint;
}
;

inline void DualElem::update_adjoint() {
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		prev_adjoint[i] = adjoint[i];
}
;

inline double* DualElem::get_func_sens() {
	return func_sens;
}

inline FluxJac& DualElem::get_flx_jac_cont() {

	return flx_jac_cont;
}

inline Matrix<double, 2, 5>& DualElem::get_hslope_sens() {
	return hslope_sens;
}

inline ErrorElem** DualElem::get_son_addresses() {
	return mysons;
}
;

class ErrorElem: public Element {

public:

	ErrorElem(Element* element);

	//used for refinement
	ErrorElem(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[], int gen,
	    int elm_loc_in[], int gen_neigh[], int mat, ErrorElem *fthTemp, double *coord_in,
	    HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr);

	//used for unrefinement
	ErrorElem(ErrorElem* sons[], HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr);

	ErrorElem(ErrElemPack* elem2, HashTable* HT_Node_Ptr, MatProps *matprops_ptr, int myid);

	ErrorElem(gzFile& myfile, HashTable* NodeTable, MatProps* matprops_ptr, int myid);

	void update(ErrElemPack* elem2, HashTable* HT_Node_Ptr, int myid);

	void get_slopes_prev(HashTable* El_Table, HashTable* NodeTable, double gamma);

	void error_update_state(SolRec* solrec, int iter);

	void error_check_refine_unrefine(SolRec* solrec, HashTable* El_Table, int iter,
	    ElemPtrList<ErrorElem>* refinelist, ElemPtrList<ErrorElem>* unrefinelist);

	double* get_bilin_adj();

	double* get_correction();

	double* get_residual();

	double* get_adjoint();

	double* get_bilin_state();

	double* get_bilin_prev_state();

	DualElem* get_father_address();

	void put_father_address(DualElem* fth);

	void xdirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
	    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]);

	void ydirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
	    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]);

	void Pack_element(ErrElemPack* elem, HashTable* HT_Node_Ptr, int destination_proc);

	void write_elem(gzFile& myfile);

private:

	//! adjoint vector
	double adjoint[NUM_STATE_VARS];

	//! this array is for bilinear reconstruction of adjoints
	double bilin_adj[NUM_STATE_VARS];

	//! this array is for bilinear reconstruction of state_vars
	double bilin_state[NUM_STATE_VARS];

	//! this array is for bilinear reconstruction of state_vars
	double bilin_prev_state[NUM_STATE_VARS];

	//! residual vector from error compute
	double residual[NUM_STATE_VARS];

	//! this term holds the correction term computed from inner product of linear construction of residual and linear construction of adjoint
	double correction;

	// this is to keep the address of his father in Dual grid
	// keeping this address let us to avoid calling lookup function
	// for many times
	DualElem* myfather;

};

inline double* ErrorElem::get_bilin_adj() {
	return bilin_adj;
}
;

inline double* ErrorElem::get_bilin_state() {
	return bilin_state;
}
;

inline double* ErrorElem::get_bilin_prev_state() {
	return bilin_prev_state;
}
;

inline double* ErrorElem::get_correction() {
	return &correction;
}
;

inline double* ErrorElem::get_residual() {
	return residual;
}
;

inline double* ErrorElem::get_adjoint() {
	return adjoint;
}
;

inline DualElem* ErrorElem::get_father_address() {
	return myfather;
}
;

inline void ErrorElem::put_father_address(DualElem* fth) {
	myfather = fth;
	return;
}
;

#endif /* DUALMESH_H */
