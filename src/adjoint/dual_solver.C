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
#include "../header/exvar.h"
#if HAVE_HDF5
#include "../header/GMFG_hdfapi.h"
#endif

#define DEBUG1

double FUNC_VAR[] = { 0., 0. };

template<typename T1, typename T2>
void copy_hashtables_objects(HashTable* El_Table, HashTable* cp_El_Table) {

	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				T1 *Curr = (T1*) (currentPtr->value);
				T2 *cp = new T2(Curr);
				cp_El_Table->add(cp->pass_key(), cp);
				currentPtr = currentPtr->next;

			}
		}
}

template<typename T1>
void delete_hashtables_objects(HashTable* El_Table) {
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				T1 *Curr = (T1*) (currentPtr->value);
				delete Curr;
				currentPtr = currentPtr->next;

			}
		}

	delete El_Table;
}

//Timers for dual, dual_init, dual_adapt and dual_repart have been set such that their timings include the error
//part as well, so to compute the timing for dual itself for these timing the error part has to be considered

Timer dual("dual"), dual_vis("dual visualization"), jacobian("jacobian"), adjoint_sol(
    "adjoint solver"), dual_repart("dual repartitioning"), dual_adapt("dual adaption"),
    read_solution("read solution"), dual_init("dual initialization"), dual_neigh_update(
        "dual repart. neighbor update");

#ifdef Error
Timer error("error"), error_init("error initialization"), error_repart("error repartitioning"),
    error_adapt("error adaption"), bilin_interp("bilinear interpolation"), error_comp(
        "error computation"), read_dual("read from dual"), update_error("updating error grid"),
    error_vis("error visualization"), error_neigh_update("error repart. neighbor update");
#endif

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	timeprops_ptr->adjust_save_time();

	const int maxiter = timeprops_ptr->iter;
	timeprops_ptr->update_savetime();

	dual_init.start();

	if (myid == 0)

		cout << endl << "******************************\n" << "*          start             *\n"
		    << "*     ADJOINT SOLUTION       *\n" << "*                            *\n"
		    << "*                            *\n" << "******************************\n"
				    "Generating grids for computing the Adjoint and the Error ....\n";

	HashTable *Dual_El_Tab = new HashTable(El_Table);

	MeshCTX dual_meshctx;
	dual_meshctx.el_table = Dual_El_Tab;
	dual_meshctx.nd_table = NodeTable;

	copy_hashtables_objects<Element, DualElem>(El_Table, Dual_El_Tab);

	if (myid == 0)
		cout << "The Adjoint grid has been generated ....\n";

	calc_adjoint(&dual_meshctx, propctx);

//	compute_functional_variation(&dual_meshctx, propctx);

	dual_vis.start();
	write_dual_xdmf(Dual_El_Tab, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_NEW, 2);
	dual_vis.stop();

#ifdef Error
	error.start();
	error_init.start();

	HashTable *Err_El_Tab = new HashTable(El_Table);
	HashTable *Err_Nod_Tab = new HashTable(NodeTable);

	copy_hashtables_objects<Element, ErrorElem>(El_Table, Err_El_Tab);
	copy_hashtables_objects<Node, Node>(NodeTable, Err_Nod_Tab);

	MeshCTX error_meshctx;
	error_meshctx.el_table = Err_El_Tab;
	error_meshctx.nd_table = Err_Nod_Tab;
	if (myid == 0)
		cout << "The Error grid has been generated ....\n";

	error_init.stop();
	error.stop();
#endif

	delete_hashtables_objects<Element>(El_Table);

	if (myid == 0)
		cout << "computing ADJOINT time step " << maxiter << endl;

	dual_init.stop();
#ifdef Error
	error.start();
	// this function refines and do constant reconstruction
	error_init.start();
	uinform_refine(&error_meshctx, propctx);

	make_dual_err_link(Dual_El_Tab, Err_El_Tab);
	error_init.stop();

	read_dual.start();
	send_from_dual_to_error(Dual_El_Tab, Err_El_Tab, 1);

	move_err_data(&error_meshctx, propctx);
	read_dual.stop();

	bilin_interp.start();
	bilinear_interp(Err_El_Tab);

	move_err_data(&error_meshctx, propctx);
	bilin_interp.stop();

	update_error.start();
	update_bilinear_error_grid(&error_meshctx, propctx);
	update_error.stop();

	error_comp.start();
	error_compute(&error_meshctx, propctx);
	error_comp.stop();

	error_vis.start();
	write_err_xdmf(Err_El_Tab, Err_Nod_Tab, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_NEW, 2);
	error_vis.stop();
	error.stop();
#else
	MeshCTX error_meshctx;
#endif

	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		if (myid == 0)
			cout << "computing ADJOINT time step " << iter - 1 << endl;
		timeprops_ptr->adjiter++;

		setup_dual_flow(solrec, &dual_meshctx, &error_meshctx, propctx);

		timeprops_ptr->adjoint_time(iter - 1);

		jacobian.start();
		calc_jacobian(&dual_meshctx, propctx);

		comminucate_jacobians(&dual_meshctx, propctx);
		jacobian.stop();

		adjoint_sol.start();
		calc_adjoint(&dual_meshctx, propctx);
		adjoint_sol.stop();
//		cout << "test of adjoint: " << simple_test(Dual_El_Tab, timeprops_ptr, matprops_ptr) << endl;

		compute_param_sens(&dual_meshctx, propctx);
		compute_functional_variation(&dual_meshctx, propctx);

#ifdef Error
		error.start();

		read_dual.start();
		send_from_dual_to_error(Dual_El_Tab, Err_El_Tab, 0);

		move_err_data(&error_meshctx, propctx);
		read_dual.stop();

		bilin_interp.start();
		bilinear_interp(Err_El_Tab);

		move_err_data(&error_meshctx, propctx);
		bilin_interp.stop();

		update_error.start();
		update_bilinear_error_grid(&error_meshctx, propctx);
		update_error.stop();

		error_comp.start();
		error_compute(&error_meshctx, propctx);
		error_comp.stop();

		error_vis.start();
		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/) {
			write_err_xdmf(Err_El_Tab, Err_Nod_Tab, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD,
			    1);
			print_func_var(propctx);
		}
		error_vis.stop();
//		}
		error.stop();
#else
		dual_vis.start();
//		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)
		write_dual_xdmf(Dual_El_Tab, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD,
				1);
		dual_vis.stop();
#endif
		if (timeprops_ptr->ifsave_adj()) {
			move_dual_data(&dual_meshctx, propctx);
			move_err_data(&error_meshctx, propctx);
			save_dual(&dual_meshctx, &error_meshctx, propctx, solrec);
		}
//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
//sensitivity w.r.t to parameters

// in dual weighted error estimation if solver performs n step, we'll have n+1
// solution and n+1 adjoint solution, but we'll have just n residual and as a
// result n error estimate. The point is that at initial step (0'th step),
// we know the solution from initial condition  so the error of 0th step is zero,
// and we have to compute the error for other time steps.
	}

	delete_hashtables_objects<DualElem>(Dual_El_Tab);
	delete_hashtables_objects<Node>(NodeTable);
	close_xdmf_files(myid);

#ifdef Error
	error.start();
	delete_hashtables_objects<ErrorElem>(Err_El_Tab);
	delete_hashtables_objects<Node>(Err_Nod_Tab);
	error.stop();
#endif

	delete_hashtables_objects<Jacobian>(solrec);
}

void dual_solver(SolRec* solrec, MeshCTX* dual_meshctx, MeshCTX* error_meshctx, PropCTX* propctx) {

	HashTable* Dual_El_Tab = dual_meshctx->el_table;
	HashTable* NodeTable = dual_meshctx->nd_table;

	HashTable* Err_El_Tab = error_meshctx->el_table;
	HashTable* Err_Nod_Tab = error_meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	const int maxiter = timeprops_ptr->iter - 1;

	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		if (myid == 0)
			cout << "computing ADJOINT time step " << iter - 1 << endl;
		timeprops_ptr->adjiter++;

		setup_dual_flow(solrec, dual_meshctx, error_meshctx, propctx);

		timeprops_ptr->adjoint_time(iter - 1);

		jacobian.start();
		calc_jacobian(dual_meshctx, propctx);

		comminucate_jacobians(dual_meshctx, propctx);
		jacobian.stop();

		adjoint_sol.start();
		calc_adjoint(dual_meshctx, propctx);
		adjoint_sol.stop();
//		cout << "test of adjoint: " << simple_test(Dual_El_Tab, timeprops_ptr, matprops_ptr) << endl;

		compute_param_sens(dual_meshctx, propctx);
		compute_functional_variation(dual_meshctx, propctx);

#ifdef Error
		error.start();

		read_dual.start();
		send_from_dual_to_error(Dual_El_Tab, Err_El_Tab, 0);

		move_err_data(error_meshctx, propctx);
		read_dual.stop();

		bilin_interp.start();
		bilinear_interp(Err_El_Tab);

		move_err_data(error_meshctx, propctx);
		bilin_interp.stop();

		update_error.start();
		update_bilinear_error_grid(error_meshctx, propctx);
		update_error.stop();

		error_comp.start();
		error_compute(error_meshctx, propctx);
		error_comp.stop();

		error_vis.start();
		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/) {
			write_err_xdmf(Err_El_Tab, Err_Nod_Tab, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD,
			    1);
			print_func_var(propctx);
		}
		error_vis.stop();
//		}
		error.stop();
#else
		dual_vis.start();
//		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)
		write_dual_xdmf(Dual_El_Tab, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD,
				1);
		dual_vis.stop();
#endif
		if (timeprops_ptr->ifsave_adj()) {
			move_dual_data(dual_meshctx, propctx);
			move_err_data(error_meshctx, propctx);
			save_dual(dual_meshctx, error_meshctx, propctx, solrec);
		}
//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
//sensitivity w.r.t to parameters

// in dual weighted error estimation if solver performs n step, we'll have n+1
// solution and n+1 adjoint solution, but we'll have just n residual and as a
// result n error estimate. The point is that at initial step (0'th step),
// we know the solution from initial condition  so the error of 0th step is zero,
// and we have to compute the error for other time steps.
	}

	delete_hashtables_objects<DualElem>(Dual_El_Tab);
	delete_hashtables_objects<Node>(NodeTable);
	close_xdmf_files(myid);

#ifdef Error
	error.start();
	delete_hashtables_objects<ErrorElem>(Err_El_Tab);
	delete_hashtables_objects<Node>(Err_Nod_Tab);
	error.stop();
#endif

	delete_hashtables_objects<Jacobian>(solrec);

}

bool must_write(MemUse* memuse_ptr, int myid) {

	struct sysinfo memInfo;
	sysinfo(&memInfo);
	unsigned long totalPhysMem = memInfo.totalram;
	unsigned long freeram = memInfo.freeram;
	unsigned long current_physMemUsed = memInfo.totalram - memInfo.freeram;
	long long last_time_step_use = current_physMemUsed - memuse_ptr->usedmem;

//	if (memuse_ptr->usedmem > 0) {
//		memuse_ptr->usedmem = current_physMemUsed;
//
//		double ratio = ((long double) last_time_step_use) / ((long double) memInfo.freeram);
	double ratio1 = ((long double) current_physMemUsed) / ((long double) totalPhysMem);
	double ratio2 = ((long double) freeram) / ((long double) totalPhysMem);

	double glob_ratio = 0.;

	MPI_Allreduce(&ratio1, &glob_ratio, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if (myid == 0)
//		cout << "ratio of last time use over free mem " << ratio << "  and ratio od used mem" << endl;
		printf(" ratio of used mem over total mem is %4f \n", ratio1);

	if (glob_ratio < .05)
		return true;
//	} else
//		memuse_ptr->usedmem = current_physMemUsed;

	return false;

}

void print_timings(int myid) {

	if (myid == 0) {
		cout << "\n=========== TIMING of PRIMAL PROBLOM =========\n";
		initialization_f.print();
		adaption.print();
		repartition_f.print();
		stept.print();
		write_solution.print();
		visualization.print();
		cout << "\n============ TIMING of DUAL PROBLOM ==========\n";
		dual_init.print();
		dual_adapt.print();
		dual_repart.print();
		dual_neigh_update.print();
		jacobian.print();
		adjoint_sol.print();
		read_solution.print();
		dual_vis.print();
#ifdef Error
		cout << "\n=========== TIMING of ERROR PROBLOM ==========\n";
		error_init.print();
		error_adapt.print();
		error_repart.print();
		error_neigh_update.print();
		update_error.print();
		read_dual.print();
		bilin_interp.print();
		error_comp.print();
		error_vis.print();
#endif
		cout << "\n===============  TOTAL TIMING  ===============\n";
		primal.print();
		dual.print();
#ifdef Error
		error.print();
#endif
		total.print();
	}
}

void compute_functional_variation(MeshCTX* dual_meshctx, PropCTX* propctx) {

	HashTable* El_Table = dual_meshctx->el_table;
	HashTable* NodeTable = dual_meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	HashEntryPtr currentPtr;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					double* adjoint = Curr_El->get_adjoint();
					double* phi_sens = Curr_El->get_phi_sens();
					double* pint_sens = Curr_El->get_pint_sens();

					// the minus sign comes from the adjoint equation
					FUNC_VAR[0] += -(adjoint[1] * phi_sens[1] + adjoint[2] * phi_sens[2]);

					FUNC_VAR[1] += -(adjoint[1] * pint_sens[1] + adjoint[2] * pint_sens[2]);

//					int cc = 0, bb = 1;
//					for (int j = 0; j < 2; ++j)
//						if (isnan(FUNC_VAR[j]))
//							bb = cc;

				}
				currentPtr = currentPtr->next;
			}
		}

	double global_funcvar[2];

	MPI_Allreduce(&FUNC_VAR[0], &global_funcvar[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (numprocs > 1)
		for (int i = 0; i < 2; ++i)
			FUNC_VAR[i] = global_funcvar[i];
}

void print_func_var(PropCTX* propctx) {

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

//	double hscale = matprops_ptr->HEIGHT_SCALE;
//	double lscale = matprops_ptr->LENGTH_SCALE;
//	double tscale = timeprops_ptr->TIME_SCALE;
//	double gsacel = matprops_ptr->GRAVITY_SCALE;
//
//	double velocity_scale = sqrt(lscale * gsacel); // scaling factor for the velocities
//	double momentum_scale = hscale * velocity_scale; // scaling factor for the momentums

	char filename[50] = "func_var.out";

	FILE *file = fopen(filename, "a");

	double functional_scale = (matprops_ptr->LENGTH_SCALE) * (matprops_ptr->LENGTH_SCALE)
	    * (matprops_ptr->HEIGHT_SCALE);

	fprintf(file, "%d %8.8f %8.8f %8.8f\n", timeprops_ptr->iter,
	    timeprops_ptr->time * timeprops_ptr->TIME_SCALE, FUNC_VAR[0] * functional_scale,
	    FUNC_VAR[1] * functional_scale);

	fclose(file);
}

