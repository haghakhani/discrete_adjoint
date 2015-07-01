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

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx, PertElemInfo* eleminfo) {

	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	const int rescomp = 1;
	const double increment = INCREMENT;
	const int maxiter = timeprops_ptr->iter;

	// dummy variables
	int order_flag = 1;
	double outflow = 0;

	allocJacoMat(El_Table);

	double functional = 0.0, dt;

	calc_adjoint(meshctx, propctx);

	uinform_refine(meshctx, propctx, numprocs, myid);

	error_compute(meshctx, propctx, maxiter, myid, numprocs);

	double UNREFINE_TARGET = .01;	//dummy value is not used in the function
	unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr, matprops_ptr,
	    rescomp);

	int tecflag = 2;
	tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, functional, tecflag);

	tecflag = 1;

	for (int iter = maxiter; iter > 0; --iter) {

		allocJacoMat(El_Table);

		timeprops_ptr->iter = iter;
		dt = timeprops_ptr->dt.at(iter - 1);
		timeprops_ptr->adjiter++;

		// we need this even for  iter = maxiter because after refine and unrefine
		// the state variables are not same as forward run
		reverse_states(El_Table, solrec, iter);

		timeprops_ptr->adjoint_time(iter - 1);

		//this array holds ResFlag for element itself and its neighbors
		ResFlag resflag;
		resflag.callflag = 1;
		resflag.lgft = 0;

		calc_d_gravity(El_Table);

		calc_edge_states(El_Table, NodeTable, matprops_ptr, timeprops_ptr, myid, &order_flag, &outflow,
		    resflag);
		slopes(El_Table, NodeTable, matprops_ptr, 1);

		compute_functional(El_Table, &functional, timeprops_ptr);

		eleminfo->update_dual_func(functional);

		calc_jacobian(meshctx, propctx, eleminfo, INCREMENT);

//		print_jacobian(El_Table, iter);

		calc_adjoint(meshctx, propctx);

		if (eleminfo->iter == iter - 1)
			fill_pertelem_info(El_Table, eleminfo);

		//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
		//sensitivity w.r.t to parameters

		uinform_refine(meshctx, propctx, numprocs, myid);

		error_compute(meshctx, propctx, iter, myid, numprocs);

		// in dual weighted error estimation if solver performs n step, we'll have n+1
		// solution and n+1 adjoint solution, but we'll have just n residual and as a
		// result n error estimate. The point is that at initial step (0'th step),
		// we know the solution from initial condition  so the error of 0th step is zero,
		// and we have to compute the error for other time steps.

		double UNREFINE_TARGET = .01;	//dummy value is not used in the function
		unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr, matprops_ptr,
		    rescomp);

		if (timeprops_ptr->adjiter/*timeprops_ptr->ifadjoint_out() /*|| adjiter == 1*/)
			tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, functional,
			    tecflag);

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

//void initSolRec(HashTable* El_Table, HashTable* NodeTable, DualMesh *dualmesh, double dt,
//    int myid) {
//
//	HashEntryPtr* buck = El_Table->getbucketptr();
//	HashEntryPtr currentPtr;
//	Element* Curr_El;
//	int num = 0;
//
//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) { // this part allocate memory and initialize jacobian matrices inside the corresponding Jacobian
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//				if (Curr_El->get_adapted_flag() > 0) {
//
//					Solution *solution = new Solution(Curr_El->get_state_vars(), *(Curr_El->get_kactxy()));
//
//					dualmesh->update_sol(Curr_El, solution);
//
//				}
//				currentPtr = currentPtr->next;
//			}
//		}
//	}
//	return;
//}

void record_solution(MeshCTX* meshctx, PropCTX* propctx, SolRec* solrec) {

	HashTable* El_Table = meshctx->el_table;

	TimeProps* timeptr = propctx->timeprops;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El;

	if (timeptr->iter == 0) {

		for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) {
						Jacobian *jacobian = new Jacobian(Curr_El->pass_key(), Curr_El->get_coord());

						if (*(Curr_El->get_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_state_vars(),
							    *(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter);

						} else
							jacobian->put_solution(solrec->get_zero_solution(), timeptr->iter);

						solrec->add(jacobian->get_key(), jacobian);

					}
					currentPtr = currentPtr->next;
				}
			}
		}

	} else {

		for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) {
						Jacobian *jacobian = (Jacobian *) solrec->lookup(Curr_El->pass_key());

						if (!jacobian) {
							jacobian = new Jacobian(Curr_El->pass_key(), Curr_El->get_coord());
							solrec->add(jacobian->get_key(), jacobian);

						}
						if (*(Curr_El->get_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_state_vars(),
							    *(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter);

						} else
							jacobian->put_solution(solrec->get_zero_solution(), timeptr->iter);

					}
					currentPtr = currentPtr->next;
				}
			}
		}
	}
}

//void initSolRec(HashTable* El_Table, HashTable* NodeTable, DualMesh *dualmesh, double dt,
//    int myid) {
//
//	HashEntryPtr* buck = El_Table->getbucketptr();
//	HashEntryPtr currentPtr;
//	Element* Curr_El;
//	int num = 0;
//
//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) { // this part allocate memory and initialize jacobian matrices inside the corresponding Jacobian
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//				if (Curr_El->get_adapted_flag() > 0) {
//
//					Solution *solution = new Solution(Curr_El->get_state_vars(), *(Curr_El->get_kactxy()));
//
//					dualmesh->update_sol(Curr_El, solution);
//
//				}
//				currentPtr = currentPtr->next;
//			}
//		}
//	}
//	return;
//}

void allocJacoMat(HashTable *El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->new_jacobianMat();
				currentPtr = currentPtr->next;
			}
		}

}

double tiny_sgn(double num, double tiny) {
	if (dabs(num) < tiny)
		return 0.;
	else if (num > tiny)
		return 1.;
	else
		return -1.;
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

void reverse_states(HashTable* El_Table, HashTable* solrec, int iter) {

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
					Jacobian *jacobian = (Jacobian *) solrec->lookup(Curr_El->pass_key());
					Curr_El->rev_state_vars(jacobian, iter);

				}
				currentPtr = currentPtr->next;
			}
		}
	}

}

void print_jacobian(HashTable* El_Table, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->print_jacobian(iter);

				currentPtr = currentPtr->next;
			}
		}
	}

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

}
