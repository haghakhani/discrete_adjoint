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

void dual_solver(DualMesh* dualmesh, MatProps* matprops_ptr, TimeProps* timeprops_ptr,
    MapNames *mapname_ptr, PertElemInfo* eleminfo) {

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

	dualmesh->allocCellMem(); //this function allocates memory for all dual cells

	reverse_states(dualmesh, maxiter);

	dualmesh->initialize_dual_flow(matprops_ptr);

	calc_adjoint(dualmesh, timeprops_ptr, timeprops_ptr->iter, adjiter, myid);

//	uinform_refine(El_Table, NodeTable, timeprops_ptr, matprops_ptr, numprocs, myid);
//
//	error_compute(El_Table, NodeTable, timeprops_ptr, matprops_ptr, maxiter, myid, numprocs);
//
//	double UNREFINE_TARGET = .01;	//dummy value is not used in the function
//	unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr, matprops_ptr,
//	    rescomp);

	int plotflag = 2;
	dualplot(dualmesh, matprops_ptr, timeprops_ptr, mapname_ptr, functional, plotflag);
//	tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, functional, tecflag);

	plotflag = 1;

	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		dt = timeprops_ptr->dt.at(iter - 1);
		adjiter++;

		// we need this even for  iter = maxiter because after refine and unrefine
		// the state variables are not same as forward run
		reverse_states(dualmesh, iter);

		timeprops_ptr->adjoint_time(iter - 1);

//		dualmesh->initialize_dual_flow(matprops_ptr);
		dualmesh->calc_flux();
		dualmesh->calc_slopes();

		compute_functional(dualmesh, &functional, timeprops_ptr);

		eleminfo->update_dual_func(functional);

		calc_jacobian(dualmesh, matprops_ptr, timeprops_ptr, mapname_ptr, increment);

//		print_jacobian(dualmesh, iter);

		calc_adjoint(dualmesh, timeprops_ptr, iter, adjiter, myid);

		if (eleminfo->iter == iter - 1)
			fill_pertelem_info(dualmesh, eleminfo);

		//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
		//sensitivity w.r.t to parameters

//		uinform_refine(El_Table, NodeTable, timeprops_ptr, matprops_ptr, numprocs, myid);
//
//		error_compute(El_Table, NodeTable, timeprops_ptr, matprops_ptr, iter, myid, numprocs);
//
//		// in dual weighted error estimation if solver performs n step, we'll have n+1
//		// solution and n+1 adjoint solution, but we'll have just n residual and as a
//		// result n error estimate. The point is that at initial step (0'th step),
//		// we know the solution from initial condition  so the error of 0th step is zero,
//		// and we have to compute the error for other time steps.
//
//		double UNREFINE_TARGET = .01;	//dummy value is not used in the function
//		unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr, matprops_ptr,
//		    rescomp);

		if (adjiter/*timeprops_ptr->ifadjoint_out() /*|| adjiter == 1*/)
			dualplot(dualmesh, matprops_ptr, timeprops_ptr, mapname_ptr, functional, plotflag);

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

void initSolRec(HashTable* El_Table, HashTable* NodeTable, DualMesh *dualmesh, double dt,
    int myid) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El;
	int num = 0;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) { // this part allocate memory and initialize jacobian matrices inside the corresponding Jacobian
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					Solution *solution = new Solution(Curr_El->get_state_vars(), *(Curr_El->get_kactxy()));

					dualmesh->update_sol(Curr_El, solution);

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

void orgSourceSgn(Element* cell, double frictiny, double* orgSgn) {

}

void orgSourceSgn(DualCell* cell, double frictiny, double* orgSgn) {

	double* d_state_vars_x = cell->get_d_state_vars();
	double* d_state_vars_y = d_state_vars_x + NUM_STATE_VARS;
	double* prev_state_vars = cell->get_prev_state_vars();
	double h_inv;
	double tmp = 0.0;
	double velocity[2];

	for (int i = 0; i < DIMENSION; i++)
		orgSgn[i] = 0.0;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		velocity[0] = prev_state_vars[1] / prev_state_vars[0];
		velocity[1] = prev_state_vars[2] / prev_state_vars[0];

	} else {
		for (int k = 0; k < DIMENSION; k++)
			velocity[k] = 0.;
	}

	if (prev_state_vars[0] > 0.0)
		h_inv = 1. / prev_state_vars[0];

	tmp = h_inv * (d_state_vars_y[1] - velocity[0] * d_state_vars_y[0]);
	orgSgn[0] = tiny_sgn(tmp, frictiny);

	tmp = h_inv * (d_state_vars_x[2] - velocity[1] * d_state_vars_x[0]);
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
void reverse_states(DualMesh* dualmesh, int iter) {

	int Ny = dualmesh->get_Ny();
	int Nx = dualmesh->get_Nx();

	for (int i = 0; i < Ny; ++i)
		for (int j = 0; j < Nx; ++j) {
			DualCell* dualcell = dualmesh->get_dualcell(i, j);
			dualcell->rev_state_vars(iter);
		}
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

void print_jacobian(DualMesh*dualmesh, int iter) {

	int Ny = dualmesh->get_Ny();
	int Nx = dualmesh->get_Nx();
	DualCell* dualcell;

	for (int i = 0; i < Ny; ++i)
		for (int j = 0; j < Nx; ++j) {
			dualcell = dualmesh->get_dualcell(i, j);
			dualcell->print_jacobian(iter);

		}
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

void compute_functional(DualMesh* dualmesh, double* functional, TimeProps* timeprops_ptr) {

	double const *state_vars;
	double const *prev_state_vars;
	double dt;

	dt = timeprops_ptr->dt.at(timeprops_ptr->iter - 1);

//	printf("iter=%4d  dt=%8f \n", timeprops_ptr->iter, dt);

	int Ny = dualmesh->get_Ny();
	int Nx = dualmesh->get_Nx();
	double dx = dualmesh->get_dx();
	double dy = dualmesh->get_dy();

	for (int i = 0; i < Ny; ++i)
		for (int j = 0; j < Nx; ++j) {

			state_vars = (dualmesh->get_dualcell(i, j))->get_state_vars();
			prev_state_vars = (dualmesh->get_dualcell(i, j))->get_prev_state_vars();

			*functional += .5 * (state_vars[0] * state_vars[0] + prev_state_vars[0] * prev_state_vars[0])
			    * dx * dy * dt;

		}

//	cout << "functional is: " << *functional << endl;

}
void compute_functional(HashTable* El_Table, double* functional, TimeProps* timeprops_ptr) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	double const *dx;
	double const *state_vars;
	double const *prev_state_vars;
	double dt;

	dt = timeprops_ptr->dt.at(timeprops_ptr->iter - 1);

//	printf("iter=%4d  dt=%8f \n", timeprops_ptr->iter, dt);

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
					    * (state_vars[0] * state_vars[0] + prev_state_vars[0] * prev_state_vars[0]) * dx[0]
					    * dx[1] * dt;

				}
				currentPtr = currentPtr->next;
			}
		}

//	cout << "functional is: " << *functional << endl;

	return;
}
