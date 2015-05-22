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

//#include "../header/jacobian.h"
#include "../header/hpfem.h"
//#include <math.h>

Solution::Solution(double* curr_sol, double* kactxy,
		const double* funcsensitivity) {

	for (int i = 0; i < 4; i++)
		states[i] = curr_sol[i];
	for (int i = 0; i < 2; i++)
		kact[i] = kactxy[i];

	for (int i = 0; i < 3; i++)
		funcsens[i] = funcsensitivity[i];
}
double* Solution::get_solution() {
	return states;
}

double* Solution::get_kact() {
	return kact;
}

double* Solution::get_funcsens() {
	return (funcsens);
}

Jacobian::Jacobian(const int myid, unsigned* key, double* position) {
	for (int i = 0; i < 2; ++i) {
		Jacobian::key[i] = key[i];
		Jacobian::position[i] = position[i];
	}

	Jacobian::myid = myid;
	jacobianMat = NULL;
}

void Jacobian::new_jacobianMat() //in forward run we just save the solution and in backward run we compute the jacobian
{
	int i, j, k;
	jacobianMat = new double**[5];
	for (i = 0; i < 5; i++) {
		jacobianMat[i] = new double*[3];
		for (j = 0; j < 3; j++)
			jacobianMat[i][j] = new double[3];
	}

	for (i = 0; i < 5; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				jacobianMat[i][j][k] = 0.0;

	return;
}

Jacobian::~Jacobian() {
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 3; ++j)
			delete[] jacobianMat[i][j];

		delete[] jacobianMat[i];
	}
	delete[] jacobianMat;
}

void Jacobian::set_jacobian(int neigh_num, double elemjacob[3],
		int state_vars_num, const double incr) {

	int i, j;

	if (state_vars_num < 1) //since state_vars=1 is for first component of adjoint
		i = state_vars_num;
	else
		i = state_vars_num - 1;

	for (j = 0; j < 3; j++)
		jacobianMat[neigh_num][i][j] = elemjacob[j] / incr;

	return;
}

void Jacobian::set_jacobian() {

	for (int i = 0; i < 5; i++)
		for (int j = 3; j < 3; j++)
			for (int k = 0; k < 3; k++)
				jacobianMat[i][j][k] = 0.0;

	return;
}

double*** Jacobian::get_jacobian() {
	return jacobianMat;
}

void Jacobian::print_jacobian(int iter) {

	cout << "iter:  " << iter << '\n';
	cout << "key1:  " << key[0] << "  key2:  " << key[1] << '\n';
	cout << "X:  " << position[0] << "  Y:  " << position[1] << '\n';
	cout << "Jacobian: " << '\n';
	//cout << "self"<<"
	for (int i = 0; i < 5; i++) {
		cout << "Matrix=  " << i << "," << '\n';

		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				cout << scientific << setw(10) << setprecision(8)
						<< jacobianMat[i][j][k] << "  ";
				if (dabs(jacobianMat[i][j][k]) > 10.)
					cout << "Jedi begir mano" << endl;
			}
			cout << '\n';
		}
	}
	return;
}

void Jacobian::put_solution(Solution* sol) {
	solvector.push_back(sol);

	return;
}

void Jacobian::rev_state_vars(void* elementin, int iter) {

	double *state_vars, *prev_state_vars, *kactxy;
	Element* element;
	element = (Element*) elementin;

	state_vars = element->get_state_vars();
	prev_state_vars = element->get_prev_state_vars();
	kactxy = element->get_kactxy();

	for (int i = 0; i < NUM_STATE_VARS - 2; i++)
		state_vars[i] = *((solvector.at(iter))->get_solution() + i);

	for (int i = 0; i < NUM_STATE_VARS - 2; i++)
		prev_state_vars[i] = *((solvector.at(iter - 1))->get_solution() + i);

	for (int i = 0; i < 2; i++)
		kactxy[i] = *((solvector.at(iter))->get_kact() + i);

	for (int i = 0; i < 3; ++i)
		prev_state_vars[6 + i] = state_vars[6 + i];

	return;
}

double* Jacobian::get_funcsens(int iter) {
	return (solvector.at(iter)->get_funcsens());
}

void Jacobian::set_jacobianMat_zero(int jacmatind) {

	for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
			jacobianMat[jacmatind][j][k] = 0.0;

	return;
}

void Jacobian::add_state_func_sens(double* func_sens_prev, int iter) {

	for (int ind = 0; ind < 3; ind++)
		*(get_funcsens(iter) + ind) += func_sens_prev[ind];

	return;
}

