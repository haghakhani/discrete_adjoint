/*
 * adjoint.h
 *
 *  Created on: Mar 3, 2016
 *      Author: haghakha
 */

#ifndef SRC_HEADER_ADJOINT_H_
#define SRC_HEADER_ADJOINT_H_

#include <vector>
#define Error

//! this function transfers information during events such as ghost element data exchange and repartitioning
void move_dual_data(MeshCTX* meshctx, PropCTX* propctx);

void move_err_data(MeshCTX* meshctx, PropCTX* propctx);

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

void comminucate_jacobians(MeshCTX* meshctx, PropCTX* propctx);

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

int table_members(HashTable *NodeTable);

void make_dual_err_link(HashTable *Dual_El_Tab, HashTable *Err_El_Tab);

void send_from_dual_to_error(HashTable *Dual_El_Tab, HashTable *Err_El_Tab, int last);

void set_link(ErrorElem* son0, HashTable* Dual_Table, HashTable* Err_Table);

void check_link(MeshCTX* err_meshctx, MeshCTX* dual_meshctx);

void correct_dual_err_link(MeshCTX* err_meshctx, MeshCTX* dual_meshctx);

void correct_dual_err_link(MeshCTX* err_meshctx, MeshCTX* dual_meshctx,
    vector<pair<unsigned, unsigned> >& imported_elem);

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

void check_elem_size(HashTable *El_Table);

void myround(double *num);

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

void refinement_report(HashTable* El_Table, int myid);

void refine_flag_report(HashTable* El_Table, int myid);

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

bool must_write(MemUse* memuse_ptr, int myid);

extern double max_jac;

double simple_test(HashTable* El_Table, TimeProps* timeprops, MatProps* matprops_ptr);

void plot_ithm(HashTable* El_Table);

void set_ithm(HashTable* El_Table);

extern Mat3x3 ZERO_MATRIX;

extern double min_gen, min_dx[2];

void set_mins(HashTable *El_Table);

void compute_dx(HashTable *El_Table);

void usefull_link();

void set_send_receive_proc(int count, int myid, int numprocs, int& receive_from, int& send_to);

void free_mpi_types();

void SetTransPack(TRANSKEY& transelem, unsigned* key, unsigned* father_key, unsigned son_key[][2]);

void check_received_keys(SolRec* solrec, vector<TRANSKEY>& keys_to_check_vec,
    vector<int>& keys_status, int iter, int &found);

void dual_repartition(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx);

void dual_err_repartition(SolRec* solrec, MeshCTX* dual_meshctx, MeshCTX* err_meshctx,
    PropCTX* propctx);

void delete_extra_nodes(HashTable* El_Table, HashTable* NodeTable);

void update_neighbor_proc(PropCTX* propctx, HashTable* El_Table, double * allKeyRange);

void save_forward(const MeshCTX& meshctx, const PropCTX& propctx, SolRec *solrec);

extern Timer dual, dual_vis, jacobian, adjoint_sol, dual_repart, dual_adapt, read_solution,
    dual_init, dual_neigh_update;

extern Timer error, error_init, error_repart, error_adapt, bilin_interp, error_comp, read_dual,
    update_error, error_vis, error_neigh_update;

extern Timer primal, stept, adaption, visualization, write_solution, repartition_f,
    initialization_f, total;

void print_timings(int myid);

void make_refine_unrefine_list_from_father(MeshCTX* dual_meshctx, MeshCTX* err_meshctx,
    ElemPtrList<DualElem> *refinelist, ElemPtrList<DualElem> *unrefinelist,
    ElemPtrList<ErrorElem> *err_refinelist, ElemPtrList<ErrorElem> *err_unrefinelist);

void uniform_refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int nump, int myid,
    void* RL, TimeProps* timeprops_ptr);

void adjust_node_info(MeshCTX* meshctx, PropCTX* propctx);

unsigned* makekey(unsigned k1, unsigned k2);

void check_the_list(vector<ErrorElem*> imported_elem, HashTable* El_Table);

//===========function that are used for the test mode========================
void perturbU(HashTable* El_Table, PertElemInfo* pelinf, int iter);

void find_test_elem(HashTable* El_Table, PertElemInfo** pelinf, int iter);

void find_adjoint_sol(HashTable* El_Table, vector<Jacobian*>* solHyst, PertElemInfo* eleminfo,
    TimeProps* timeprops_ptr);

int check_elem_exist(HashTable *El_Table, unsigned *key);

void fill_pertelem_info(HashTable* El_Table, PertElemInfo* eleminfo);

#endif /* SRC_HEADER_ADJOINT_H_ */
