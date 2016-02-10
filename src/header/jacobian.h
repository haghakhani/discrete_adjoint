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

#include <map>

class Solution {

public:

	Solution(double* curr_sol, double kactxy);

	double* get_solution();

	double get_kact();

	~Solution();

	static Solution solution_zero;

protected:

	Solution();

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

	void clear_container(int iter);

//	unsigned get_create_time(){
//		return create_time;
//	}

	//destructor
	virtual ~Jacobian();

	//members
protected:

	map<int, Solution*> solContainer;
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
	// this function should call destructor of solution, so there is no need to call them explicitly
	solContainer.erase(iter);
}

inline void Jacobian::clear_container() {

	solContainer.clear();
}

inline bool Jacobian::is_container_empty() {

	return solContainer.empty();

}

#endif
