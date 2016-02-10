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
 * $Id: geoflow.h 152 2007-06-27 20:29:54Z dkumar $ 
 */

#ifndef __GEOFLOW
#define __GEOFLOW

/* geoflow header file */
#define WEIGHT_ADJUSTER 1

#define NUM_FREEFALLS_2_STOP 2 //stopping criteria parameter
//#define STOPCRIT_CHANGE_FLUX
//#define STOPCRIT_CHANGE_BED
//#define STOPCRIT_CHANGE_SOURCE
//#define DO_EROSION

//#define REFINE_LEVEL 3
extern int REFINE_LEVEL; //make REFINE_LEVEL a global variable that can be changed set it in  Read_grid() (datread.C) or loadrun() (restart.C)
//(mdj)2007-04-11 #define MIN_GENERATION -1 //minimum refinement level
#define MIN_GENERATION -3 //minimum refinement level

//! non member C++ function that wraps the fortran correct_() function
void correct(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
    FluxProps *fluxprops, TimeProps *timeprops, Element *EmTemp, double *forceint, double *forcebed,
    double *eroded, double *deposited);

//! this function is legacy, the prototype exists but the function is not defined
void checknodesol(HashTable*);

//! This function assigns a global_weight to the collection of elements based on the sum of their element_weight
double element_weight(HashTable* El_Table, HashTable*, int myid, int nump);

//! This function calculates the vast majority of statistics used for output, including most of what appears in output_summary.######, the friction body forces however are not calculated in here, Keith wrote this to replace calc_volume()
void calc_stats(HashTable* El_Table, HashTable* NodeTable, int myid, MatProps* matprops,
    TimeProps* timeprops, StatProps* statprops, DISCHARGE* discharge, double d_time);

//! calc_volume() has been replaced by calc_stats(), calc_volume() is out of date legacy code, the function is still defined in step.C but it is not called.
void calc_volume(HashTable* El_Table, int myid, MatProps* matprops_ptr, TimeProps* timeprops_ptr,
    double d_time, double* v_star, double* nz_star);

//! get_max_momentum() is legacy, it has been replaced by calc_stats()
double get_max_momentum(HashTable* El_Table, MatProps* matprops_ptr);

//! this function prints a warning message at end of the simulation to say if the flow is still moving and thus should be run longer before using the data to make decisions
void sim_end_warning(HashTable* El_Table, MatProps* matprops_ptr, TimeProps* timeprops_ptr,
    double v_star);

//! this function outputs final stats for one run in a collection of stochastic/probabilistic runs
void out_final_stats(TimeProps* timeprops_ptr, StatProps* statprops_ptr);

//! this function loops through the nodes zeroing the fluxes, then loops through the elements and finds the positive x direction of the element, calculates element size, calculates local terrain elevation, slopes, and curvatures, and calculates the gravity vector in local coordinates.
void setup_geoflow(HashTable* El_Table, HashTable* NodeTable, int myid, int nump,
    MatProps* matprops_ptr, TimeProps *timeprops_ptr);

void calc_d_gravity(HashTable* El_Table);

//! this function calculates the spatial derivatives of the state variables
void slopes(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int dualcall);

//! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model) calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum allowable timestep for this iteration.
double get_coef_and_eigen(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    FluxProps* fluxprops_ptrs, TimeProps* timeprops_ptr, int ghost_flag);

//! this function transfers information during events such as ghost element data exchange and repartitioning
void move_data(int nump, int myid, HashTable* El_Table, HashTable* NodeTable,
    TimeProps* timeprops_ptr);

//! this function deletes the current ghost elements
void delete_ghost_elms(HashTable* El_Table, int myid);

//! This function loops through all the non-ghost current elements and calls the Element member function Element::calc_edge_states() which calculates the Riemann fluxes between elements and stores the Riemann fluxes in the edge nodes.
template<class T>
void calc_edge_states(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    TimeProps* timeprops_ptr, int myid, int* order_flag, double *outflow);

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx);

void restore(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El, MatProps* matprops_ptr,
    int effelement, int j, int myid, double increment);

void restore(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El, MatProps* matprops_ptr,
    int effelement, int j, int myid, double increment, double* d_state_vars_old);

void record_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key, MatProps* matprops_ptr,
    int myid, double fluxold[4][NUM_STATE_VARS]);

void increment_state(HashTable* El_Table, Element* Curr_El, double increment, int effelement, int j,
    int* updateflux, int* srcflag, ResFlag resflag[EFF_ELL]);

void calc_flux_slope_kact(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El,
    MatProps* matprops_ptr, int myid, int effelement, int updateflux, int srcflag,
    ResFlag resflag[5]);

void zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int order_flag,
    int dir, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], Element *EmNeigh,
    double dt, ResFlag resflag);

void calc_jacobian(MeshCTX* meshctx, PropCTX* propctx);

void calc_jacobian_old(MeshCTX* meshctx, PropCTX* propctx);

void calc_jacobian_elem(Mat3x3& jacobian, const Mat3x3& jac_flux_n_x, const Mat3x3& jac_flux_p_x,
    const Mat3x3& jac_flux_n_y, const Mat3x3& jac_flux_p_y, double* prev_state_vars,
    double* d_state_vars_x, double* d_state_vars_y, double *curvature, double* gravity,
    double* d_gravity, double* dh_sens, double int_fric, double bedfrict, double kact, int effelem,
    double dtdx, double dtdy, double dt, int* stop, double* OrgSgn);

void error_compute(MeshCTX* meshctx, PropCTX* propctx);

void update_error_grid(MeshCTX* meshctx, PropCTX* propctx);

void update_bilin_error_grid(MeshCTX* meshctx, PropCTX* propctx);

void dual_unrefine(MeshCTX* meshctx, PropCTX* propctx);

void calc_adjoint_elem(MeshCTX* meshctx, PropCTX* propctx, DualElem *Curr_El);

void calc_adjoint(MeshCTX* meshctx, PropCTX* propctx);

void uinform_refine(MeshCTX* meshctx, PropCTX* propctx);

void bilinear_interp(HashTable* El_Table);

int num_nonzero_elem(HashTable *El_Table);

void make_dual_err_link(HashTable *Dual_El_Tab, HashTable *Err_El_Tab);

void send_from_dual_to_error(HashTable *Dual_El_Tab, HashTable *Err_El_Tab);

void update_dual_err_link(HashTable *Dual_El_Tab, HashTable *Err_El_Tab);

void adjoint_init(HashTable* BT_Elem_Ptr, HashTable* BT_Node_Ptr);

void compute_funcsens(Element* element, double dt, double* func_sens);

void compute_functional(HashTable* El_Table, double* functional, TimeProps* timeprops_ptr);

void residual(double* residual, double *state_vars, double *prev_state_vars, //3
    double *fluxxp, double *fluxyp, double *fluxxm, double *fluxym, double dtdx, //5
    double dtdy, double dt, double *d_state_vars_x, double *d_state_vars_y, //4
    double *curvature, double intfrictang, double bedfrictin, double *gravity, //4
    double *dgdx, double kactxyelem, double fric_tiny, double* orgSrcSgn, //4
    double increment, double epsilon, int* check_stop_crit, int srcflag = 1, int org_res_flag = 1); //5

void residual(double* residual, double *state_vars, double *prev_state_vars, double *fluxxp,
    double *fluxyp, double *fluxxm, double *fluxym, double dtdx, double dtdy, double dt,
    double *d_state_vars_x, double *d_state_vars_y, double *curvature, double intfrictang,
    double bedfrict, double *gravity, double *dgdx, double kactxyelem, double fric_tiny, int* stop,
    double* orgSrcSgn);

void update_states(double *state_vars, double *prev_state_vars, //2
    double *fluxxp, double *fluxyp, double *fluxxm, double *fluxym, double dtdx, //5
    double dtdy, double dt, double *d_state_vars_x, double *d_state_vars_y, //4
    double *curvature, double intfrictang, double bedfrict, double *gravity, //4
    double *dgdx, double kactxyelem, double fric_tiny, int* stop, double* orgSrcSgn); //5

void save_solution(HashTable* El_Table, HashTable* NodeTable, HashTable* solHystPtr,
    TimeProps* timeprops_ptr);

void allocJacoMat(HashTable *El_Table);

//this function is such that returns 0 for element itself, 1 for xp element, 2 yp element, 3 xm element, 4 ym element
int jac_mat_index(int effelement, int xp);

void orgSourceSgn(Element* Curr_El, double frictiny, double* orgSgn);

int checkElement(HashTable *El_Table, HashTable* NodeTable, double *max, unsigned *key);

int num_nonzero_elem(HashTable *El_Table, int type);

void bilinear_interp(HashTable* El_Table);

void bilinear_interp_elem(ErrorElem *elem11, ErrorElem *elem21, ErrorElem *elem12,
    ErrorElem *elem22, ErrorElem *Curr_El);

int void_neigh_elem(HashTable* El_Table, Element* Curr_El, int effelement);

void get_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key, MatProps* matprops_ptr,
    int myid, double fluxold[4][NUM_STATE_VARS]);

void flux_debug(Element* Curr_El, double* fluxxpold, double* fluxxmold, double* fluxypold,
    double* fluxymold, double* fluxxp, double* fluxxm, double* fluxyp, double* fluxym,
    int effelement, int j, int iter, double dt);

void print_jacobian(HashTable* El_Table, int iter);

void calc_flux(MeshCTX* meshctx, PropCTX* propctx);

void refinement_report(HashTable* El_Table);

template<typename T>
void dual_refine_unrefine(MeshCTX* meshctx, PropCTX* propctx, ElemPtrList<T>* refinelist,
    ElemPtrList<T>* unrefinelist);

void force_unrefine(HashTable* El_Table, HashTable* NodeTable, Element* Curr_el, int myid,
    MatProps* matprops_ptr, void* NewFatherList, void* OtherProcUpdate, int rescomp);

void update_error_grid(SolRec* solrec, MeshCTX* cp_meshctx, PropCTX* propctx);

void update_bilinear_error_grid(MeshCTX* meshctx, PropCTX* propctx);

void update_dual_grid(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx);

void clear_empty_jacobians(SolRec* solrec, int iter);

void correct_neighbor_information(Element* newfth, Element* newfth2, Element* bros[4],
    Element* bros_2[4], int neigh);

void set_new_fathers(HashTable* El_Table, vector<pair<unsigned, unsigned> >& new_father);

void reset_adaption_flag(HashTable* El_Table);

void reset_newson_adaption_flag(HashTable* El_Table);

void reset_newfather_adaption_flag(HashTable* El_Table);

void setup_dual_flow(SolRec* solrec, MeshCTX* meshctx, MeshCTX* cp_meshctx, PropCTX* propctx);

double initial_dot(HashTable* El_Table);

void clean_jacobian(HashTable* El_Table);

void copy_hashtables(HashTable* El_Table, HashTable* NodeTable, HashTable* cp_El_Table,
    HashTable* cp_NodeTable);

void check_state_vars_with_record(HashTable* El_Table, SolRec* solrec, int iter);

void print_Elem_Table(HashTable* El_Table, HashTable* NodeTable, int iter, int place);

bool must_write(MemUse* memuse_ptr);

extern double max_jac;

double simple_test(HashTable* El_Table, TimeProps* timeprops, MatProps* matprops_ptr);

void plot_ithm(HashTable* El_Table);

void set_ithm(HashTable* El_Table);

extern Mat3x3 ZERO_MATRIX;

void usefull_link();

//===========function that are used for the test mode========================
void perturbU(HashTable* El_Table, PertElemInfo* pelinf, int iter);

void find_test_elem(HashTable* El_Table, PertElemInfo** pelinf, int iter);

void find_adjoint_sol(HashTable* El_Table, vector<Jacobian*>* solHyst, PertElemInfo* eleminfo,
    TimeProps* timeprops_ptr);

int check_elem_exist(HashTable *El_Table, unsigned *key);

void fill_pertelem_info(HashTable* El_Table, PertElemInfo* eleminfo);

double tiny_sgn(double num, double tiny);

//! c++ sgn function 
inline double c_sgn(double zz) {
	double sgn;
	if (zz > GEOFLOW_TINY)
		sgn = 1;
	else if (zz < -GEOFLOW_TINY)
		sgn = -1;
	else
		sgn = 0;

	return sgn;
}

//! c++ dmin1 function
inline double c_dmin1(double d1, double d2) {

	if (d1 > d2)
		d1 = d2;

	return d1;
}

//! another c++ dmin1 function
inline double c_dmin1(double d1, double d2, double d3) {

	if (d1 > d2)
		d1 = d2;
	if (d1 > d3)
		d1 = d3;

	return d1;
}

//! a c++ dmax1 function
inline double c_dmax1(double d1, double d2) {

	if (d1 < d2)
		d1 = d2;

	return d1;
}

//! another c++ dmax1 function
inline double c_dmax1(double d1, double d2, double d3) {

	if (d1 < d2)
		d1 = d2;
	if (d1 < d3)
		d1 = d3;

	return d1;
}

// a c++ dabs function 
inline double dabs(double dd) {
	if (dd < 0.)
		dd = -dd;

	return dd;
}

/* fortran calls */
#ifdef SUNOS 
//! the actual calculation of k active passive is done by a fortran call this should be ripped out and rewritten as a C++ Element member function
extern "C" void gmfggetcoef_(double*, double*, double*, double*, double*,
		double*, double*, double*, double*, double*);

//! the actual calculation of wave speeds (eigen vectors of the flux jacoboians) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void eigen_(double *Uvec, double *eigenvxmax, double *eigenvymax,
		double *evalue, double *tiny, double *kactxy,
		double *gravity, double *Vs, double *eps);

//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void predict_(double *Uvec, double *dUdx, double *dUdy,
		double *Uprev, double *tiny, double *kactxy,
		double *dt2, double *g, double *curv,
		double *bedfrictang, double *intfrictang,
		double *dgdx, double *frict_tiny, int *order_flag,
		double *VxVy, int *if_stopped, double *fluxcoef);

//! the actual corrector timestep update 
extern "C" void correct_(double *Uvec, double *Uprev, double *fluxxp,
		double *fluxyp, double *fluxxm, double *fluxym,
		double *tiny, double *dtdx, double *dtdy, double *dt,
		double *dUdx, double *dUdy, double *xslope,
		double *yslope, double *curv, double *intfrictang,
		double *bedfrictang, double *g, double *d_gravity,double *kactxy,
		double *frict_tiny, double *forceint, double *forcebed, double *dragf,
		int *do_erosion, double *eroded,
		double *Vsolid, double *terminal_vel, double *eps,
		int *if_stopped, double *fluxcoef);
#endif
#ifdef IBMSP
extern "C" void gmfggetcoef(double*, double*, double*, double*, double*,
		double*, double*, double*, double*, double*);
extern "C" void eigen(double *Uvec, double *eigenvxmax, double *eigenvymax,
		double *evalue, double *tiny, double *kactxy,
		double *gravity, double *VxVy);

#endif

#endif
