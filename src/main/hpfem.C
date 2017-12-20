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
 * $Id: hpfem.C 211 2009-06-16 20:02:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#if HAVE_HDF5
#include "../header/GMFG_hdfapi.h"
#endif

int REFINE_LEVEL = 3;

Mat3x3 ZERO_MATRIX;

double min_gen = 10000., min_dx[] = { 10000., 10000. };

Timer primal("primal"), stept("step"), adaption("forward adaption"), visualization(
    "forward visualization"), write_solution("writing solution"), repartition_f(
    "forward repartition"), initialization_f("forward initialization"), total("total");

int main(int argc, char *argv[]) {

	total.start();
	primal.start();
	initialization_f.start();
	MPI_Init(&argc, &argv);

	HashTable *Node_Table, *Err_Node_Table;
	HashTable *El_Table, *Err_El_Table;

	SolRec* solrec;
	MemUse memuse;
	memuse.usedmem = 0;

	//-- MPI
	int myid, numprocs;
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Get_processor_name(processor_name, &namelen);

	/* create new MPI datastructures for class objects */
	MPI_New_Datatype();

	/* read original data from serial preprocessing
	 code and then initialize element
	 stiffness routines info */
	int material_count = 0;
	double epsilon = 1., intfrictang = 1, *bedfrictang = NULL, gamma = 1;
	double frict_tiny = .1, mu = 1.0e-03, rho = 1600;
	char **matnames = NULL;

	StatProps statprops;
	MatProps matprops(material_count, matnames, intfrictang, bedfrictang, mu, rho, epsilon, gamma,
	    frict_tiny, 1.0, 1.0, 1.0);
	TimeProps timeprops;
	timeprops.starttime = time(NULL);
	timeprops.REFINE = 5;

	MapNames mapnames;
	PileProps pileprops;
	FluxProps fluxprops;
	OutLine outline;
	DISCHARGE discharge;

	int adaptflag;

	vector<Snapshot> snapshot_vec;

	int viz_flag = 0, order_flag; //savefileflag will be flipped so first savefile will end in 0
	int srctype;

	Read_data(myid, &matprops, &pileprops, &statprops, &timeprops, &fluxprops, &adaptflag, &viz_flag,
	    &order_flag, &mapnames, &discharge, &outline, &srctype);

	run_mode runcond = loadrun(myid, numprocs, &Node_Table, &El_Table, &Err_Node_Table, &Err_El_Table,
	    &solrec, &matprops, &timeprops, &outline, argv[1], argv[2]);

	if (runcond & NORMAL) {
		Read_grid(myid, numprocs, &Node_Table, &El_Table, &matprops, &timeprops, &outline, &solrec);

		setup_geoflow(El_Table, Node_Table, myid, numprocs, &matprops, &timeprops);

		move_data(numprocs, myid, El_Table, Node_Table, &timeprops,&matprops);

		AssertMeshErrorFree(El_Table, Node_Table, numprocs, myid, -1.0);

		//initialize pile height and if appropriate perform initial adaptation
		init_piles(El_Table, Node_Table, myid, numprocs, adaptflag, &matprops, &timeprops, &mapnames,
		    &pileprops, &fluxprops, &statprops);

		setup_geoflow(El_Table, Node_Table, myid, numprocs, &matprops, &timeprops);
	}

	MeshCTX meshctx;
	meshctx.el_table = El_Table;
	meshctx.nd_table = Node_Table;
	meshctx.snapshot_vec = &snapshot_vec;

	PropCTX propctx;
	propctx.timeprops = &timeprops;
	propctx.matprops = &matprops;
	propctx.mapnames = &mapnames;
	propctx.outline = &outline;
	propctx.discharge = &discharge;
	propctx.statprops = &statprops;
	propctx.fluxprops = &fluxprops;
	propctx.numproc = numprocs;
	propctx.myid = myid;
	propctx.adapt_flag = adaptflag;
	propctx.runcond = (run_mode) (runcond | RECORD);

	//we need to record the initial configuration
	meshctx.snapshot_vec->push_back(Snapshot(meshctx,propctx));

	if ((runcond & NORMAL) || (runcond & FORWARD)) {
		forward_solve(meshctx, propctx, solrec);
	}

	dual.start();

	MeshCTX error_meshctx;
	error_meshctx.el_table = Err_El_Table;
	error_meshctx.nd_table = Err_Node_Table;

	timeprops.adjiter=timeprops.maxiter;

	dual_solver(solrec, &meshctx, &error_meshctx, &propctx, runcond);

	dual.stop();

	delete_data(solrec, &meshctx, &error_meshctx, &propctx);
	free_mpi_types();

	total.stop();

	print_timings(myid);

	MPI_Finalize();
	return (0);
}

