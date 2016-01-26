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
 * $Id: jacobian.h  2014-02-5 10:09:10 haghakha $ 
 */

#ifndef JACOBIAN_H
#define JACOBIAN_H

//class DualMesh;

#include <unordered_map>

class Solution {

public:

	Solution(double* curr_sol, double kactxy);

	double* get_solution();

	double get_kact();

	~Solution();

protected:

	double states[NUM_STATE_VARS]; //to save the solution
	double kact; //to save kact

};

inline double* Solution::get_solution() {
	return states;
}

inline double Solution::get_kact() {
	return kact;
}

class Jacobian {
	//friend functions and classes

public:

	//constructors
	Jacobian(unsigned* key, double* position);

	Jacobian(unsigned* key);

	double* get_position();

	unsigned* get_key();

	void put_solution(Solution* sol, int iter);

	Solution* get_solution(int iter);

	void erase_solution(int iter);

	bool is_container_empty();

	void clear_container();

//	unsigned get_create_time(){
//		return create_time;
//	}

	//destructor
	virtual ~Jacobian();

	//members
protected:

	unordered_map<int, Solution*> solContainer;
	unsigned key[DIMENSION];
//	unsigned create_time;
//	double position[DIMENSION];

};

//inline double* Jacobian::get_position() {
//	return position;
//}

inline void Jacobian::put_solution(Solution* solution, int iter) {

	//because insert is faster than emplace
	solContainer[iter] = solution;

	return;
}

inline Solution* Jacobian::get_solution(int iter) {

	return solContainer[iter];
}

inline unsigned* Jacobian::get_key() {
	return key;
}

inline void Jacobian::erase_solution(int iter) {
	Solution* sol = get_solution(iter);
	delete sol;
	// this function should call destructor of solution, so there is no need to call them explicitly
	solContainer.erase(iter);
}

inline void Jacobian::clear_container() {

	solContainer.clear();
}

inline bool Jacobian::is_container_empty() {

	return solContainer.empty();

}

//class DualCell: public Jacobian {
//
//public:
//
//	DualCell(unsigned* key, double* position);
//
//	void allocMem();
//
//	double* get_state_vars();
//
//	double* get_prev_state_vars();
//
//	double* get_curr_adjoint();
//
//	double* get_prev_adjoint();
//
//	double* get_curvature();
//
//	double* get_gravity();
//
//	double* get_d_state_vars();
//
//	double get_elevation();
//
//	double* get_d_gravity();
//
//	double* get_xflux();
//
//	double* get_yflux();
//
//	double* get_zeta();
//
//	double get_kact();
//
//	void rev_state_vars(int iter);
//
//	void xdir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft);
//
//	void ydir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft);
//
//	void zdir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft, int side);
//
//	void update_flux(DualMesh* dualmesh, int* lgft, int side);
//
//	void update_flux_x(DualMesh* dualmesh, int* lgft);
//
//	void update_flux_y(DualMesh* dualmesh, int* lgft);
//
//	void update_flux(DualMesh* dualmesh, int* lgft);
//
//	void calc_gravity(MatProps* matprops_ptr);
//
//	void calc_d_gravity(DualMesh* dualmesh);
//
//	void calc_slopes(DualMesh* dualmesh);
//
//	void calc_topo_data(DualMesh* dualmesh, MatProps* matprops_ptr);
//
//	double* get_funcsens();
//
//	void calc_func_sens(const void * ctx);
//
//	void print_cell_info(int iter);
//
//	void print_cell_neighb_info(DualMesh*dualmesh, int iter);
//
//	//destructor
//	~DualCell();
//
//private:
//
//	double* curr_adjoint;
//	double* prev_adjoint;
//	double* func_sens;
//	double* state_vars;
//	double* prev_state_vars;
//	double* flux;
//	double* gravity;
//	double* d_gravity;
//	double* d_state_vars;
//	double* zeta;
//	double* curvature;
//
//	double kact;
//	double elevation;
//
//};
//
//void set_ab(int* a, int* b, int effelement);

#endif
