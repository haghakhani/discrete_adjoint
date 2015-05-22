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
 * $Id: calc_jacobian.C 164 2013-02-27 15:27:22Z haghakhani $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define DEBUG1

#define KEY0    3777862041
#define KEY1    2576980374
#define EFFELL  0
#define ITER    187
#define J       0

void dual_solver(HashTable* El_Table, HashTable* NodeTable,
		vector<Jacobian*>* solHyst, MatProps* matprops_ptr,
		TimeProps* timeprops_ptr, MapNames *mapname_ptr, PertElemInfo* eleminfo) {

	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	const int rescomp = 1;
	const double increment = INCREMENT;
	const int maxiter = timeprops_ptr->iter;

//	// here we do this because iter in timeprops is such that it is one iter more than
//	// actual iteration at the end of forward run, so we have to correct that.
//	timeprops_ptr->iter = timeprops_ptr->maxiter;

	double functional = 0.0, dt;

	int adjiter = 0;

	int hrs, mins;
	double secs;

	allocJacoMat(*solHyst); //this function allocates memory to store Jacobian matrices
	//int unsigned key[2] = { KEY0, KEY1 };

	calc_adjoint(El_Table, solHyst, maxiter, adjiter, myid);

	uinform_refine(El_Table, NodeTable, timeprops_ptr, matprops_ptr, numprocs,
			myid);

	error_compute(El_Table, NodeTable, timeprops_ptr, matprops_ptr, maxiter, myid,
			numprocs);

	double UNREFINE_TARGET = .01;	//dummy value is not used in the function
	unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr,
			matprops_ptr, rescomp);

	int tecflag = 2;
	tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr,
			functional, tecflag);

	tecflag = 1;

	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		dt = timeprops_ptr->dt.at(iter - 1);
		adjiter++;

		// we need this even for  iter = maxiter because after refine and unrefine
		// the state variables are not same as forward run
		reverse_states(El_Table, solHyst, iter);

		timeprops_ptr->adjoint_time(iter - 1);

		setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr,
				timeprops_ptr);

		compute_functional(El_Table, &functional, timeprops_ptr);
		eleminfo->update_dual_func(functional);

		calc_jacobian(El_Table, NodeTable, solHyst, matprops_ptr, timeprops_ptr,
				mapname_ptr, increment);

//		print_jacobian(El_Table, solHyst, iter);

		calc_adjoint(El_Table, solHyst, iter, adjiter, myid);

		if (eleminfo->iter == iter - 1)
			fill_pertelem_info(El_Table, solHyst, eleminfo);

		//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
		//sensitivity w.r.t to parameters

		uinform_refine(El_Table, NodeTable, timeprops_ptr, matprops_ptr, numprocs,
				myid);

		error_compute(El_Table, NodeTable, timeprops_ptr, matprops_ptr, iter, myid,
				numprocs);

		// in dual weighted error estimation if solver performs n step, we'll have n+1
		// solution and n+1 adjoint solution, but we'll have just n residual and as a
		// result n error estimate. The point is that at initial step (0'th step),
		// we know the solution from initial condition  so the error of 0th step is zero,
		// and we have to compute the error for other time steps.

		double UNREFINE_TARGET = .01;	//dummy value is not used in the function
		unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs,
				timeprops_ptr, matprops_ptr, rescomp);

		if (/*adjiter*/timeprops_ptr->ifadjoint_out() /*|| adjiter == 1*/)
			tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr,
					functional, tecflag);

	}

	return;
}

int num_nonzero_elem(HashTable *El_Table) {
	int num = 0;			//myid
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

void initSolRec(HashTable* El_Table, HashTable* NodeTable,
		vector<Jacobian*> *solHyst, TimeProps* timeprops_ptr, int myid) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El;
	Jacobian* jacobian;
	double functionalsens[3] = { 0., 0., 0. };
	int num = 0;

	solHyst->reserve(num_nonzero_elem(El_Table));

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) { // this part allocate memory and initialize jacobian matrices inside the corresponding Jacobian
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					jacobian = new Jacobian(myid, Curr_El->pass_key(),
							Curr_El->get_coord());
					// at time step 0 we do not compute functional sensitivity,
					// we compute the contribution of this time step n functional sensitivity on time step 1
//					compute_funcsens(Curr_El, timeprops_ptr, functionalsens);
					Solution *solution = new Solution(Curr_El->get_state_vars(),
							Curr_El->get_kactxy(), functionalsens);
					Curr_El->put_sol_rec_ind(num);
					jacobian->put_solution(solution);
					solHyst->push_back(jacobian);
					num++;

				}
				currentPtr = currentPtr->next;
			}
		}
	}
	return;
}

void allocJacoMat(vector<Jacobian*> solHyst) {

	vector<Jacobian*>::iterator it;
	for (it = solHyst.begin(); it != solHyst.end(); ++it)
		(*it)->new_jacobianMat();

	return;
}

double tiny_sgn(double num, double tiny) {
	if (dabs(num) < tiny)
		return 0.;
	else if (num > tiny)
		return 1.;
	else
		return -1.;
}

void orgSourceSgn(Element* Curr_El, double frictiny, double* orgSgn) {

	double* d_state_vars_x = Curr_El->get_d_state_vars();
	double* d_state_vars_y = d_state_vars_x + NUM_STATE_VARS;
	double* prev_state_vars = Curr_El->get_prev_state_vars();
	double h_inv;
	double tmp = 0.0;
	double velocity[2];
	for (int i = 0; i < 2; i++)
		orgSgn[i] = 0.0;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		velocity[0] = prev_state_vars[2] / prev_state_vars[0];
		velocity[1] = prev_state_vars[3] / prev_state_vars[0];

	} else {
		for (int k = 0; k < DIMENSION; k++)
			velocity[k] = 0.;
	}

	if (prev_state_vars[0] > 0.0)
		h_inv = 1. / prev_state_vars[0];

	tmp = h_inv * (d_state_vars_y[2] - velocity[0] * d_state_vars_y[0]);
	orgSgn[0] = tiny_sgn(tmp, frictiny);

	tmp = h_inv * (d_state_vars_x[3] - velocity[1] * d_state_vars_x[0]);
	orgSgn[1] = tiny_sgn(tmp, frictiny);

	return;
}


int num_nonzero_elem(HashTable *El_Table, int type) {
	int num = 0;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() == type)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

void reverse_states(HashTable* El_Table, vector<Jacobian*>* solHyst, int iter) {
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

#ifdef DEBUG
	if (checkElement(El_Table))
	exit(22);
	for (int i = 0; i < nonz1; i++) {
		dbgvec[i] = 0;
		pass[i] = 0;
	}
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {

					int index = Curr_El->get_sol_rec_ind();
					dbgvec[index] += 1;
					pass[index] = 1;

				}
			}
		}
	}
	for (int i = 0; i < nonz1; i++) {
		if (dbgvec[i] != 1) {
			cout << "these elements have problem:   " << i << endl
			<< "the value is  " << dbgvec[i] << endl;
			exit(EXIT_FAILURE);

		} else if (pass[i] != 1) {
			cout << "this index has not been passed  " << i << endl;
			exit(EXIT_FAILURE);

		}
	}

	delete[] dbgvec;

	cout << "number of elements after unrefinement"
	<< num_nonzero_elem(El_Table) << endl;
//		getchar();
	if (checkElement(El_Table))
	exit(22);
#endif

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					Jacobian *jacobian = solHyst->at(Curr_El->get_sol_rec_ind());

					if (iter != 0)
						jacobian->rev_state_vars(Curr_El, iter);

				}
				currentPtr = currentPtr->next;
			}
		}
	}

	return;
}

void print_jacobian(HashTable* El_Table, vector<Jacobian*>* solHyst, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					Jacobian *jacobian = solHyst->at(Curr_El->get_sol_rec_ind());
					jacobian->print_jacobian(iter);
				}
				currentPtr = currentPtr->next;
			}
		}
	}
	return;
}

void compute_functional(HashTable* El_Table, double* functional,
		TimeProps* timeprops_ptr) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	double const *dx;
	double const *state_vars;
	double const *prev_state_vars;
	double dt;

	dt = timeprops_ptr->dt.at(timeprops_ptr->iter - 1);

	printf("iter=%4d  dt=%8f \n", timeprops_ptr->iter, dt);

	//we do not have make it zero here, because we want to compute the integration over the time and space
	//*functional = 0.0;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					dx = Curr_El->get_coord();
					state_vars = Curr_El->get_state_vars();
					prev_state_vars = Curr_El->get_prev_state_vars();

					// we used trapezoidal integration on time
					// flag is for time step 0

					*functional += .5
							* (state_vars[0] * state_vars[0]
									+ prev_state_vars[0] * prev_state_vars[0]) * dx[0] * dx[1]
							* dt;

				}
				currentPtr = currentPtr->next;
			}
		}

	cout << "functional is: " << *functional << endl;

	return;
}
