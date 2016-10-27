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

		set_ithm(Dual_El_Tab);

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

//		write_jacobian_to_compute_eigen(&dual_meshctx, propctx);

		adjoint_sol.start();
		calc_adjoint(&dual_meshctx, propctx);
		adjoint_sol.stop();
//		cout << "test of adjoint: " << simple_test(Dual_El_Tab, timeprops_ptr, matprops_ptr) << endl;

//		wrtie_El_Table_ordered(&dual_meshctx, propctx, "DUAL");

//		compute_param_sens(&dual_meshctx, propctx);
//		compute_functional_variation(&dual_meshctx, propctx);

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
//			print_func_var(propctx);
		}
		error_vis.stop();
//		}
		error.stop();
#else
		dual_vis.start();
		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/) {
			write_dual_xdmf(Dual_El_Tab, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD,
			    1);
//		print_func_var(propctx);
		}
		dual_vis.stop();
#endif

		if (timeprops_ptr->ifsave_adj()) {
			move_dual_data(&dual_meshctx, propctx);

#ifdef Error
			move_err_data(&error_meshctx, propctx);
#endif

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

//	compute_init_location_variation(&dual_meshctx, propctx);
	compute_init_volume_variation(&dual_meshctx, propctx);

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

//		compute_param_sens(dual_meshctx, propctx);
//		compute_functional_variation(dual_meshctx, propctx);

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
//			print_func_var(propctx);
		}
		error_vis.stop();
//		}
		error.stop();
#else
		dual_vis.start();
		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/) {
//		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)
			write_dual_xdmf(Dual_El_Tab, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD,
			    1);
//			print_func_var(propctx);
		}
		dual_vis.stop();
#endif

		if (timeprops_ptr->ifsave_adj()) {
			move_dual_data(dual_meshctx, propctx);

#ifdef Error
			move_err_data(error_meshctx, propctx);
#endif
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

//	compute_init_location_variation(dual_meshctx, propctx);
	compute_init_volume_variation(dual_meshctx, propctx);

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

		//writing sameinformation to  a file
		ofstream fp;
		fp.open("timing.txt", ofstream::out);
		fp << "\n=========== TIMING of PRIMAL PROBLOM =========\n";
		initialization_f.write(fp);
		adaption.write(fp);
		repartition_f.write(fp);
		stept.write(fp);
		write_solution.write(fp);
		visualization.write(fp);
		fp << "\n============ TIMING of DUAL PROBLOM ==========\n";
		dual_init.write(fp);
		dual_adapt.write(fp);
		dual_repart.write(fp);
		dual_neigh_update.write(fp);
		jacobian.write(fp);
		adjoint_sol.write(fp);
		read_solution.write(fp);
		dual_vis.write(fp);
#ifdef Error
		fp << "\n=========== TIMING of ERROR PROBLOM ==========\n";
		error_init.write(fp);
		error_adapt.write(fp);
		error_repart.write(fp);
		error_neigh_update.write(fp);
		update_error.write(fp);
		read_dual.write(fp);
		bilin_interp.write(fp);
		error_comp.write(fp);
		error_vis.write(fp);
#endif
		fp << "\n===============  TOTAL TIMING  ===============\n";
		primal.write(fp);
		dual.write(fp);
#ifdef Error
		error.write(fp);
#endif
		total.write(fp);

		fp.close();

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
	int elem = 0;
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		matprops_ptr->sensitivity[i] = 0.;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					// note that at time step n adjoint_prev is for time step n+1 and adjoint is at time step n
					double* adjoint = Curr_El->get_prev_adjoint();
					double* phi_sens = Curr_El->get_phi_sens();
					double* pint_sens = Curr_El->get_pint_sens();

					// the minus sign comes from the adjoint equation
					matprops_ptr->sensitivity[0] = 0.;
					matprops_ptr->sensitivity[1] += adjoint[1] * phi_sens[1] + adjoint[2] * phi_sens[2];
					matprops_ptr->sensitivity[2] += adjoint[1] * pint_sens[1] + adjoint[2] * pint_sens[2];

//					int cc = 0, bb = 1;
////					for (int j = 0; j < 2; ++j)
//					if (isnan(matprops_ptr->sensitivity[0]) || isinf(matprops_ptr->sensitivity[0])
//					    || fabs(matprops_ptr->sensitivity[0]) > 1e6)
//						cout << "func_var " << matprops_ptr->sensitivity[0] << " adjoint [1] " << adjoint[1]
//						    << " adjoint[2] " << adjoint[2] << " phi_sens[1] " << phi_sens[1]
//						    << "  phi_sens[2] " << phi_sens[2] << "  elem " << elem << "iter "
//						    << timeprops_ptr->iter << endl;
//					elem++;

				}
				currentPtr = currentPtr->next;
			}
		}

}

void compute_init_location_variation(MeshCTX* dual_meshctx, PropCTX* propctx) {

	HashTable* El_Table = dual_meshctx->el_table;
	HashTable* NodeTable = dual_meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	HashEntryPtr currentPtr;
	HashEntryPtr *buck;

	HashTable* new_hashtab = new HashTable(El_Table);

	buck = El_Table->getbucketptr();
	//first we save everything
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					double* adjoint = Curr_El->get_adjoint();
					double* state_vars = Curr_El->get_state_vars();
					unsigned* key = Curr_El->pass_key();

					Container* contain = new Container(Curr_El);
					new_hashtab->add(key, contain);
				}
				currentPtr = currentPtr->next;
			}
		}

	//then we perturb initial location in x_dir
	PileProps* pileprops = propctx->pileprops;
	pileprops->xCen[0] += INCREMENT;

	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					elliptical_pile_height(NodeTable, Curr_El, matprops_ptr, pileprops);

				currentPtr = currentPtr->next;
			}
		}

	move_dual_data(dual_meshctx, propctx);

	calc_d_gravity(El_Table);

	slopes(El_Table, NodeTable, matprops_ptr, 1);

	double tiny = GEOFLOW_TINY;
	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					gmfggetcoef_(Curr_El->get_prev_state_vars(), Curr_El->get_d_state_vars(),
					    (Curr_El->get_d_state_vars() + NUM_STATE_VARS), Curr_El->get_dx(),
					    &(matprops_ptr->bedfrict[Curr_El->get_material()]), &(matprops_ptr->intfrict),
					    (Curr_El->get_kactxy()), (Curr_El->get_kactxy() + 1), &tiny,
					    &(matprops_ptr->epsilon));

				currentPtr = currentPtr->next;
			}
		}

	move_dual_data(dual_meshctx, propctx);

	calc_flux(dual_meshctx, propctx);

	move_dual_data(dual_meshctx, propctx);

	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					double *prev_state_vars = Curr_El->get_prev_state_vars();
					double *state_vars = Curr_El->get_state_vars();
					double *d_state_vars = Curr_El->get_d_state_vars();
					double *gravity = Curr_El->get_gravity();
					double *d_gravity = Curr_El->get_d_gravity();
					double *curvature = Curr_El->get_curvature();
					Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
					double bedfrict = Curr_El->get_effect_bedfrict();
					double *dx = Curr_El->get_dx();

					double flux[4][NUM_STATE_VARS];
					get_flux(El_Table, NodeTable, Curr_El, flux);
					int iter = timeprops_ptr->iter;
					double dt = timeprops_ptr->dt.at(iter - 1);
					double dtdx = dt / dx[0];
					double dtdy = dt / dx[1];

					int stop[DIMENSION] = { 0, 0 };
					double orgSrcSgn[4] = { 0., 0., 0., 0. };

					update_states(state_vars, prev_state_vars, flux[0], flux[1], flux[2], flux[3], dtdx, dtdy,
					    dt, d_state_vars, (d_state_vars + NUM_STATE_VARS), curvature, matprops_ptr->intfrict,
					    bedfrict, gravity, d_gravity, *(Curr_El->get_kactxy()), matprops_ptr->frict_tiny,
					    stop, orgSrcSgn);

				}

				currentPtr = currentPtr->next;
			}
		}

//	 now we can compute the sensitivities w.r.t x perturbation

	buck = new_hashtab->getbucketptr();
	for (int i = 0; i < new_hashtab->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {

				Container * container = (Container*) (currentPtr->value);

				DualElem* Curr_El = (DualElem*) El_Table->lookup(container->pass_key());

				double sens[] = { 0., 0., 0. };

				for (int i = 0; i < NUM_STATE_VARS; ++i) {
					sens[i] = (Curr_El->get_state_vars()[i] - container->get_state()[i]) / INCREMENT;
					matprops_ptr->sensitivity[1] += sens[i] * container->get_adjoint()[i];
				}
				currentPtr = currentPtr->next;
			}
		}

	//now we can compute the sens w.r.t y perturb

	pileprops->xCen[0] -= INCREMENT;
	pileprops->yCen[0] += INCREMENT;

	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					elliptical_pile_height(NodeTable, Curr_El, matprops_ptr, pileprops);

				currentPtr = currentPtr->next;
			}
		}

	move_dual_data(dual_meshctx, propctx);

	calc_d_gravity(El_Table);

	slopes(El_Table, NodeTable, matprops_ptr, 1);

	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					gmfggetcoef_(Curr_El->get_prev_state_vars(), Curr_El->get_d_state_vars(),
					    (Curr_El->get_d_state_vars() + NUM_STATE_VARS), Curr_El->get_dx(),
					    &(matprops_ptr->bedfrict[Curr_El->get_material()]), &(matprops_ptr->intfrict),
					    (Curr_El->get_kactxy()), (Curr_El->get_kactxy() + 1), &tiny,
					    &(matprops_ptr->epsilon));

				currentPtr = currentPtr->next;
			}
		}

	move_dual_data(dual_meshctx, propctx);

	calc_flux(dual_meshctx, propctx);

	move_dual_data(dual_meshctx, propctx);

	buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem * Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					double *prev_state_vars = Curr_El->get_prev_state_vars();
					double *state_vars = Curr_El->get_state_vars();
					double *d_state_vars = Curr_El->get_d_state_vars();
					double *gravity = Curr_El->get_gravity();
					double *d_gravity = Curr_El->get_d_gravity();
					double *curvature = Curr_El->get_curvature();
					Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
					double bedfrict = Curr_El->get_effect_bedfrict();
					double *dx = Curr_El->get_dx();

					double flux[4][NUM_STATE_VARS];
					get_flux(El_Table, NodeTable, Curr_El, flux);
					int iter = timeprops_ptr->iter;
					double dt = timeprops_ptr->dt.at(iter - 1);
					double dtdx = dt / dx[0];
					double dtdy = dt / dx[1];

					int stop[DIMENSION] = { 0, 0 };
					double orgSrcSgn[4] = { 0., 0., 0., 0. };

					update_states(state_vars, prev_state_vars, flux[0], flux[1], flux[2], flux[3], dtdx, dtdy,
					    dt, d_state_vars, (d_state_vars + NUM_STATE_VARS), curvature, matprops_ptr->intfrict,
					    bedfrict, gravity, d_gravity, *(Curr_El->get_kactxy()), matprops_ptr->frict_tiny,
					    stop, orgSrcSgn);

				}

				currentPtr = currentPtr->next;
			}
		}

	//	 now we can compute the sensitivities w.r.t y perturbation

	buck = new_hashtab->getbucketptr();
	for (int i = 0; i < new_hashtab->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {

				Container * container = (Container*) (currentPtr->value);

				DualElem* Curr_El = (DualElem*) El_Table->lookup(container->pass_key());

				double sens[] = { 0., 0., 0. };

				for (int i = 0; i < NUM_STATE_VARS; ++i) {
					sens[i] = (Curr_El->get_state_vars()[i] - container->get_state()[i]) / INCREMENT;
					matprops_ptr->sensitivity[2] += sens[i] * container->get_adjoint()[i];
				}
				currentPtr = currentPtr->next;
			}
		}

	print_func_var(propctx);

	delete_hashtables_objects<Container>(new_hashtab);
}

void print_func_var(PropCTX* propctx) {

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	double hscale = matprops_ptr->HEIGHT_SCALE;
	double lscale = matprops_ptr->LENGTH_SCALE;
//	double tscale = timeprops_ptr->TIME_SCALE;
//	double gsacel = matprops_ptr->GRAVITY_SCALE;
//
//	double velocity_scale = sqrt(lscale * gsacel); // scaling factor for the velocities
//	double momentum_scale = hscale * velocity_scale; // scaling factor for the momentums

	double global_sens[] = { 0., 0., 0. };

	MPI_Reduce(matprops_ptr->sensitivity, global_sens, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myid == 0) {

		FILE *file = fopen("func_var.out", "a");

		double functional_scale = (matprops_ptr->LENGTH_SCALE) * (matprops_ptr->LENGTH_SCALE)
		    * (matprops_ptr->HEIGHT_SCALE);

		fprintf(file, "%d %f %e %e %e\n", timeprops_ptr->iter,
		    timeprops_ptr->time * timeprops_ptr->TIME_SCALE, global_sens[0] * functional_scale,
		    global_sens[1] * functional_scale, global_sens[2] * functional_scale);

		fclose(file);
	}
}

void compute_init_volume_variation(MeshCTX* dual_meshctx, PropCTX* propctx) {

	HashTable* El_Table = dual_meshctx->el_table;
	HashTable* NodeTable = dual_meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	PileProps* pileprops = propctx->pileprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	int iter = timeprops_ptr->iter;
	double dt = timeprops_ptr->dt.at(iter - 1);

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashTable* new_hashtab = new HashTable(El_Table);

	double sens_x0 = 0., sens_y0 = 0., sens_majrad = 0., sens_minrad = 0.;
	double majrad = pileprops->majorrad[0], minrad = pileprops->minorrad[0],
	    xcen = pileprops->xCen[0], ycen = pileprops->yCen[0];

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

//					Container* contain = new Container(Curr_El);
//					new_hashtab->add(Curr_El->pass_key(), contain);
//					double* sensitivity = contain->get_sens();
					double sensitivity = 0.;
					double* coord = Curr_El->get_coord();

					for (int effelement = 0; effelement < EFF_ELL; ++effelement) {

						if (effelement == 0) {
							double* adjoint = Curr_El->get_prev_adjoint();

							Vec_Mat<9>& jacobianmat = Curr_El->get_jacobian();

							for (int l = 0; l < NUM_STATE_VARS; ++l)
								sensitivity += adjoint[l] * jacobianmat(effelement, l, 0);

						} else if (effelement <= 4
						    || (effelement > 4 && *(Curr_El->get_neigh_proc() + (effelement - 1)) > -2)) {

							//basically we are checking all neighbor elements, and start from xp neighbor
							DualElem * neigh_elem = (DualElem*) (El_Table->lookup(
							    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));

							if (neigh_elem) {

								double* adjoint = neigh_elem->get_prev_adjoint();
								Vec_Mat<9>& jacobianmat = neigh_elem->get_jacobian();

								int jacind = neigh_elem->which_neighbor(Curr_El->pass_key());
								// because we have to consider the element itself which is in jacind=0
								jacind++;

								for (int l = 0; l < NUM_STATE_VARS; ++l)
									sensitivity += adjoint[l] * jacobianmat(jacind, l, 0);

							}
						}
					}

					sens_x0 += sensitivity * 2. * (coord[0] - xcen) / (majrad * majrad);
					sens_y0 += sensitivity * 2. * (coord[1] - ycen) / (minrad * minrad);
					sens_majrad += sensitivity * 2. * (coord[0] - xcen) * (coord[0] - xcen)
					    / (majrad * majrad * majrad);
					sens_minrad += sensitivity * 2. * (coord[1] - ycen) * (coord[1] - ycen)
					    / (minrad * minrad * minrad);

				}
				currentPtr = currentPtr->next;
			}
		}
	}
	double g_sens_x0 = 0., g_sens_y0 = 0., g_sens_majrad = 0., g_sens_minrad = 0.;

	MPI_Reduce(&sens_x0, &g_sens_x0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sens_y0, &g_sens_y0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sens_majrad, &g_sens_majrad, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sens_minrad, &g_sens_minrad, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	double hscale = matprops_ptr->HEIGHT_SCALE;
	double lscale = matprops_ptr->LENGTH_SCALE;
	double functional_scale = lscale * lscale * hscale;

	if (myid == 0) {
		char filename[50];
		sprintf(filename, "func_var");

		FILE *file = fopen(filename, "a");
		fprintf(file, "discharge,sens_x0,sens_y0,sens_majrad,sens_minrad\n");

		fprintf(file, "%e,%e,%e,%e,%e\n", fabs(discharge) * functional_scale,
		    g_sens_x0 * functional_scale / lscale, g_sens_y0 * functional_scale / lscale,
		    g_sens_majrad * functional_scale / lscale, g_sens_minrad * functional_scale / lscale);

//	buck = new_hashtab->getbucketptr();
//	for (int i = 0; i < new_hashtab->get_no_of_buckets(); i++)
//		if (*(buck + i)) {
//			HashEntryPtr currentPtr = *(buck + i);
//			while (currentPtr) {
//
//				Container * container = (Container*) (currentPtr->value);
//
//				fprintf(file, "%f,%f,%f,%e\n", container->get_coord()[0] * lscale,
//				    container->get_coord()[1] * lscale, container->get_state()[0] * hscale,
//				    *(container->get_sens()) * lscale * lscale);
//
//				currentPtr = currentPtr->next;
//			}
//		}
		fclose(file);
	}

	delete_hashtables_objects<Container>(new_hashtab);
}

void write_jacobian_to_compute_eigen(MeshCTX* dual_meshctx, PropCTX* propctx) {

	HashTable* El_Table = dual_meshctx->el_table;

	set_ithm(El_Table);

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashTable* new_hashtab = new HashTable(El_Table);

	char filename[20];

	sprintf(filename, "Jacobian%d.csv", propctx->timeprops->iter);

	FILE* file = fopen(filename, "w");

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				DualElem* Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					for (int effel = 0; effel < EFF_ELL; ++effel) {

						Vec_Mat<9>& jacobianmat = Curr_El->get_jacobian();
						int row = Curr_El->get_ithelem();
						int col;
						if (effel == 0) {
							col = row;

							for (int sub_row = 0; sub_row < NUM_STATE_VARS; ++sub_row)
								for (int sub_col = 0; sub_col < NUM_STATE_VARS; ++sub_col)
									if (jacobianmat(effel, sub_row, sub_col) != 0.)
										fprintf(file, "%d,%d,%f\n", row * NUM_STATE_VARS + sub_row,
										    col * NUM_STATE_VARS + sub_col, jacobianmat(effel, sub_row, sub_col));

						} else if (effel <= 4
						    || (effel > 4 && *(Curr_El->get_neigh_proc() + (effel - 1)) > -2)) {

							DualElem * neigh_elem = (DualElem*) (El_Table->lookup(
							    Curr_El->get_neighbors() + (effel - 1) * KEYLENGTH));

							if (neigh_elem) {
								col = neigh_elem->get_ithelem();
								for (int sub_row = 0; sub_row < NUM_STATE_VARS; ++sub_row)
									for (int sub_col = 0; sub_col < NUM_STATE_VARS; ++sub_col)
										if (jacobianmat(effel, sub_row, sub_col) != 0.)
											fprintf(file, "%d,%d,%f\n", row * NUM_STATE_VARS + sub_row,
											    col * NUM_STATE_VARS + sub_col, jacobianmat(effel, sub_row, sub_col));

							}
						}
					}
				}
				currentPtr = currentPtr->next;
			}
		}
	}

	fclose(file);
}
