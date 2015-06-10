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

class DualMesh;

class Solution {

public:

	Solution(double* curr_sol, double kactxy);

	double* get_solution(void);

	double get_kact(void);

	~Solution();

protected:

	double states[NUM_STATE_VARS]; //to save the solution
	double kact; //to save kact

};

class Jacobian {
	//friend functions and classes

public:

	//constructors
	Jacobian(unsigned* key, double* position);

	void set_jacobian(int neigh_num, double elemjacob[3], int state_vars_num, const double incr);

	// this function sets the jacobian for a boundary element
	void set_jacobian();

	void print_jacobian(int iter);

	double*** get_jacobian(void);

	double* get_solution(void);

	void new_jacobianMat(void);

	void put_solution(Solution* sol);

	virtual void rev_state_vars(void* element, int iter);

	void set_jacobianMat_zero(int jacmatind);

	void add_state_func_sens(double* func_sens_prev, int iter);

	void del_jacobianMat();

	double* get_position();

	unsigned* get_key();

	//destructor
	virtual ~Jacobian();

	//members
protected:

	vector<Solution*> solvector;
	unsigned key[DIMENSION];
	double position[DIMENSION];
	double ***jacobianMat; //double jacobianMat [5][3][3],self[3][3],neigh1[3][3],neigh2[3][3],neigh3[3][3],neigh4[3][3]

};

class DualCell: public Jacobian {

public:

	DualCell(unsigned* key, double* position);

	void allocMem();

	double* get_state_vars();

	double* get_prev_state_vars();

	double* get_curr_adjoint();

	double* get_prev_adjoint();

	double* get_curvature();

	double* get_gravity();

	double* get_d_state_vars();

	double get_elevation();

	double* get_d_gravity();

	double* get_xflux();

	double* get_yflux();

	double* get_zeta();

	double get_kact();

	void rev_state_vars(int iter);

	void xdir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft);

	void ydir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft);

	void zdir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft, int side);

	void update_flux(DualMesh* dualmesh, int* lgft, int side);

	void update_flux_x(DualMesh* dualmesh, int* lgft);

	void update_flux_y(DualMesh* dualmesh, int* lgft);

	void update_flux(DualMesh* dualmesh, int* lgft);

	void calc_gravity(MatProps* matprops_ptr);

	void calc_d_gravity(DualMesh* dualmesh);

	void calc_slopes(DualMesh* dualmesh);

	void calc_topo_data(DualMesh* dualmesh, MatProps* matprops_ptr);

	double* get_funcsens();

	void calc_func_sens(const void * ctx);

	//destructor
	~DualCell();

private:

	double* curr_adjoint;
	double* prev_adjoint;
	double* func_sens;
	double* state_vars;
	double* prev_state_vars;
	double* flux;
	double* gravity;
	double* d_gravity;
	double* d_state_vars;
	double* zeta;
	double* curvature;

	double kact;
	double elevation;

};

#endif
