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
//#define Error

#define KEY0   3920807148
#define KEY1   1321528399
#define ITER   10
#define J      0

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

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	const int maxiter = timeprops_ptr->iter;

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

//	dualplotter(Dual_El_Tab, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 2);
	write_dual_xdmf(El_Table, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_NEW,2);
//	write_xdmf(El_Table, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_NEW);

#ifdef Error

	HashTable *Err_El_Tab = new HashTable(El_Table);
	HashTable *Err_Nod_Tab = new HashTable(NodeTable);

	copy_hashtables_objects<Element, ErrorElem>(El_Table, Err_El_Tab);
	copy_hashtables_objects<Node, Node>(NodeTable, Err_Nod_Tab);

	MeshCTX error_meshctx;
	error_meshctx.el_table = Err_El_Tab;
	error_meshctx.nd_table = Err_Nod_Tab;

#endif

	delete_hashtables_objects<Element>(El_Table);

	if (myid == 0) {

		cout << "The Error grid has been generated ....\n";

		cout << "computing ADJOINT time step " << maxiter << endl;
	}

#ifdef Error

	// this function refines and do constant reconstruction
	uinform_refine(&error_meshctx, propctx);

	make_dual_err_link(Dual_El_Tab, Err_El_Tab);

	send_from_dual_to_error(Dual_El_Tab, Err_El_Tab, 1);

	bilinear_interp(Err_El_Tab);

	update_bilinear_error_grid(&error_meshctx, propctx);

	error_compute(&error_meshctx, propctx);

	errorplotter(Err_El_Tab, Err_Nod_Tab, matprops_ptr, timeprops_ptr,
			mapname_ptr, maxiter);

#else
	MeshCTX error_meshctx;
#endif

//	set_ithm(El_Table);
//
//	plot_ithm(El_Table);

//	cout << "elements number of original grid " << num_nonzero_elem(El_Table) << endl;
//	cout << "elements number of refined grid " << num_nonzero_elem(cp_El_Table) << endl;

//	print_Elem_Table(El_Table, NodeTable, timeprops_ptr->iter, 0);
//
//	print_Elem_Table(cp_El_Table, cp_NodeTable, timeprops_ptr->iter, 1);

//	int tecflag = 2;

//	set_ithm(El_Table);
//	plot_ithm(El_Table);

//	reset_adaption_flag(Dual_El_Tab);

//	refinement_report(El_Table);
//
//	refinement_report(cp_El_Table);

//	this function reconstruct bilinear interpolation
//	set_ithm(cp_El_Table);
//	plot_ithm(cp_El_Table);
//	cout<<num_nonzero_elem(El_Table)<<endl;
//	cout<<num_nonzero_elem(cp_El_Table)<<endl;
//
//	timeprops_ptr->iter--;
//	meshplotter(cp_El_Table, cp_NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., tecflag);

//	setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr, timeprops_ptr);

//	refinement_report(Dual_El_Tab, myid);
//	refine_flag_report(Dual_El_Tab, myid);

	for (int iter = maxiter; iter > 0; --iter) {

		set_ithm(Dual_El_Tab);

		timeprops_ptr->iter = iter;
		if (myid == 0)
			cout << "computing ADJOINT time step " << iter - 1 << endl;
		timeprops_ptr->adjiter++;

		setup_dual_flow(solrec, &dual_meshctx, &error_meshctx, propctx);

//		cout << "elements number of original grid " << num_nonzero_elem(El_Table) << endl;
//		cout << "elements number of refined grid " << num_nonzero_elem(cp_El_Table) << endl;

		timeprops_ptr->adjoint_time(iter - 1);

//		compute_functional(El_Table, &functional, timeprops_ptr);

//		eleminfo->update_dual_func(functional);

		calc_jacobian(&dual_meshctx, propctx);

		comminucate_jacobians(&dual_meshctx, propctx);

//		cout<<" max jac is: "<<max_jac<<endl;

		calc_adjoint(&dual_meshctx, propctx);

//		cout << "test of adjoint: " << simple_test(Dual_El_Tab, timeprops_ptr, matprops_ptr) << endl;

#ifdef Error
		if (iter > 1) {

			send_from_dual_to_error(Dual_El_Tab, Err_El_Tab, 0);

			bilinear_interp(Err_El_Tab);

			update_bilinear_error_grid(&error_meshctx, propctx);

			error_compute(&error_meshctx, propctx);
			if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)

			errorplotter(Err_El_Tab, Err_Nod_Tab, matprops_ptr, timeprops_ptr,
					mapname_ptr, iter - 1);
		}
#endif

//		if (iter - 1 == 1)
//		cout << "test of adjoint: " << simple_test(El_Table, timeprops_ptr, matprops_ptr) << endl;

//		map<int, Vec_Mat<9>> jac_code;

//		copy_jacobian(El_Table, jac_code);

//		calc_jacobian_old(meshctx, propctx);

//		map<int, Vec_Mat<9>> jac_diff;

//		copy_jacobian(El_Table, jac_diff);

//		compare_jacobians(jac_code, jac_diff);

//		clean_jacobian(El_Table);

//		print_Elem_Table(El_Table, NodeTable, timeprops_ptr->iter, 1);

//		check_state_vars_with_record(El_Table, solrec, iter);

//		print_jacobian(El_Table, iter);

//		if (eleminfo->iter == iter - 1)
//			fill_pertelem_info(El_Table, eleminfo);

//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
//sensitivity w.r.t to parameters

//		uinform_refine(meshctx, propctx);

//		error_compute(meshctx, propctx, iter);

// in dual weighted error estimation if solver performs n step, we'll have n+1
// solution and n+1 adjoint solution, but we'll have just n residual and as a
// result n error estimate. The point is that at initial step (0'th step),
// we know the solution from initial condition  so the error of 0th step is zero,
// and we have to compute the error for other time steps.

//		dual_unrefine(meshctx, propctx);

		if (/*timeprops_ptr->adjiter*/timeprops_ptr->ifadjoint_out()/*|| adjiter == 1*/)
//			dualplotter(Dual_El_Tab, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 1);
			write_dual_xdmf(El_Table, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD, 1);
//			write_xdmf(El_Table, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_OLD);
	}
//	write_dual_xdmf(El_Table, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_CLOSE, 1);
//	write_xdmf(El_Table, NodeTable, timeprops_ptr, matprops_ptr, mapname_ptr, XDMF_CLOSE);

	delete_hashtables_objects<DualElem>(Dual_El_Tab);
	delete_hashtables_objects<Node>(NodeTable);

#ifdef Error
	delete_hashtables_objects<DualElem>(Err_El_Tab);
	delete_hashtables_objects<Node>(Err_Nod_Tab);
#endif

	delete_hashtables_objects<Jacobian>(solrec);

//	if (fabs(max_err1) > fabs(max_err2))
//		cout << "max error occurred in test 1" << max_err1 << "  at iter " << iter_1 << " key is "
//		    << key1_1 << " , " << key2_1 << endl;
//	else
//		cout << "max error occurred in test 2" << max_err2 << "  at iter " << iter_2 << " key is "
//		    << key1_2 << " , " << key2_2 << endl;

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

