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

Solution::Solution(double* curr_sol) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		states[i] = curr_sol[i];
}

Solution::Solution() {
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		states[i] = 0.;
}

Solution::~Solution() {

}

Jacobian::Jacobian(unsigned* key, double* position) {
	for (int i = 0; i < KEYLENGTH; ++i) {
		Jacobian::key[i] = key[i];
	}
}

Jacobian::Jacobian(unsigned* key) {
	for (int i = 0; i < KEYLENGTH; ++i)
		Jacobian::key[i] = key[i];
}

void Jacobian::clear_container(int iter) {
	for (map<int, Solution*>::iterator it = solContainer.begin(); it != solContainer.end(); ++it)
		if (it->first >= iter - 1) {
			delete it->second;
			solContainer.erase(it);
		}
}

Jacobian::~Jacobian() {
	for (map<int, Solution*>::iterator it = solContainer.begin(); it != solContainer.end(); ++it){
		if (it->second)
			delete it->second;
		solContainer.erase(it);
	}
}
