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
#include <math.h>
#include <vector>

class Solution {

public:

	Solution(double* curr_sol, double* kactxy, const double* funcsensitivity);
	double* get_solution(void);
	double* get_kact(void);
	double* get_funcsens(void);
	~Solution();

protected:

	double funcsens[3]; //this variable keeps the value of sensitivity at each time step for this element.
	double states[4]; //to save the solution
	double kact[2]; //to save kact

};

class Jacobian {
	//friend functions and classes

public:
	//constructors
	Jacobian(const int myid, unsigned* key, double* position);

	void set_jacobian(int neigh_num, double elemjacob[3], int state_vars_num,
			const double incr);
	// this function sets the jacobian for a boundary element
	void set_jacobian();
	void print_jacobian(int iter);
	double*** get_jacobian(void);
	double* get_solution(void);
	double* get_kact(void);
	void new_jacobianMat(void);
	double* get_funcsens(int iter);
	void put_solution(Solution* sol);
	void rev_state_vars(void* element, int iter);
	void set_jacobianMat_zero(int jacmatind);
	void add_state_func_sens(double* func_sens_prev, int iter);

	//destructor
	~Jacobian();

	//members
private:

	vector<Solution*> solvector;
	unsigned key[2];
	double position[2];
	int myid;
	double ***jacobianMat; //double jacobianMat [5][3][3],self[3][3],neigh1[3][3],neigh2[3][3],neigh3[3][3],neigh4[3][3]

};
#endif
