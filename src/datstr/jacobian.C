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
 * $Id: jacobian.C  2014-04-10 10:33:10 haghakha $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

Solution::Solution(double* curr_sol, double kactxy) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		states[i] = curr_sol[i];

	kact = kactxy;

}

double* Solution::get_solution() {
	return states;
}

double Solution::get_kact() {
	return kact;
}

Solution::~Solution() {

}

Jacobian::Jacobian(unsigned* key, double* position) {

	for (int i = 0; i < KEYLENGTH; ++i) {
		Jacobian::key[i] = key[i];
		Jacobian::position[i] = position[i];
	}

}

Jacobian::Jacobian(unsigned* key) {

	for (int i = 0; i < KEYLENGTH; ++i)
		Jacobian::key[i] = key[i];

}

double* Jacobian::get_position() {
	return position;
}

void Jacobian::put_solution(Solution* solution, int iter) {

	//because insert is faster than emplace
	solContainer[iter] = solution;

	return;
}

Solution* Jacobian::get_solution(int iter) {

	return solContainer[iter];
}

unsigned* Jacobian::get_key() {
	return key;
}

void Jacobian::erase_solution(int iter) {
Solution* sol=get_solution(iter);
delete sol;
	// this function should call destructor of solution, so there is no need to call them explicitly
	solContainer.erase(iter);
}

void Jacobian::clear_container() {

	solContainer.clear();
}

bool Jacobian::is_container_empty() {

	return solContainer.empty();

}

Jacobian::~Jacobian() {

	solContainer.clear();

}
