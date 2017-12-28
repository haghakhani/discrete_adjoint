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
 * $Id: hpfem.C 211 2009-06-16 20:02:10Z haghakha $
 */

#include "../header/hpfem.h"

void forward_solve(MeshCTX &meshctx, PropCTX &propctx, SolRec *solrec) {

	TimeProps* timeprops = propctx.timeprops;
	MapNames* mapname = propctx.mapnames;
	MatProps* matprops = propctx.matprops;
	StatProps* statprops= propctx.statprops;
	DISCHARGE* discharge = propctx.discharge;
	FluxProps* fluxprops = propctx.fluxprops;
	OutLine* outline = propctx.outline;
	PileProps* pileprops = propctx.pileprops;
	int myid = propctx.myid;
	int numprocs = propctx.numproc;
	const int adaptflag = propctx.adapt_flag;

	HashTable *El_Table = meshctx.el_table;
	HashTable *Node_Table = meshctx.nd_table;

	const int rescomp = 0;
	int adjflag = 0;
	int order_flag = 1;
	int savefileflag = 1;

	if (myid == 0) {
		for (int imat = 0; imat < matprops->material_count; imat++)
		printf("bed friction angle for \"%s\" is %g\n", matprops->matnames[imat],
				matprops->bedfrict[imat] * 180.0 / PI);

		printf("internal friction angle is %g, epsilon is %g \n",
				matprops->intfrict * 180.0 / PI, matprops->epsilon);
		printf("REFINE_LEVEL=%d\n", REFINE_LEVEL);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	calc_stats(El_Table, Node_Table, myid, matprops, timeprops, statprops, discharge, 0.0);

	if (timeprops->verbose)
	output_discharge(matprops, timeprops, discharge, myid);

	move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops);

	int OUTPUT=0;
	if (OUTPUT)
		write_xdmf(El_Table,Node_Table,timeprops,matprops,mapname,XDMF_NEW);

	initialization_f.stop();

	/*
	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 Time Stepping Loop

	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 */

	while (!(timeprops->ifend())) //(timeprops.ifend(0.5*statprops.vmean)) && !ifstop)
	{
		/*
		 *  mesh adaption routines
		 */
		double TARGET = .05;
		//double UNREFINE_TARGET = .005;
		double UNREFINE_TARGET = .01;
		int h_count = 0;
		if (timeprops->iter < 50)
		matprops->frict_tiny = 0.1;
		else
		matprops->frict_tiny = 0.000000001;

		if ((propctx.adapt_flag != 0) && timeprops->ifrefine()) {

			adaption.start();

			H_adapt(El_Table, Node_Table, h_count, TARGET, matprops, fluxprops, timeprops, 5);

			move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops);

			unrefine(El_Table, Node_Table, UNREFINE_TARGET, myid, numprocs, timeprops, matprops,
					rescomp);

			move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work
			adaption.stop();

			if ((numprocs > 1) && timeprops->ifrepartition()) {
				repartition_f.start();

				repartition2(El_Table, Node_Table, timeprops, matprops);

				move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work
				repartition_f.stop();
			}

			calc_d_gravity(El_Table);
		}

		stept.start();

		step(El_Table, Node_Table, myid, numprocs, matprops, timeprops, pileprops, fluxprops,
				statprops, &order_flag, outline, discharge, adaptflag);

		stept.stop();

//		char filename[50];
//		sprintf(filename,"forward_%d_%d",timeprops->iter,myid);
//
//		write_alldata_ordered(El_Table, filename);
//
//		write_solution.start();

		solrec->record_solution(&meshctx, &propctx);

//		if (solrec->write_sol()/* || must_write(&memuse, myid)*/) {
//			solrec->wrtie_sol_to_disk(myid);
//
//			solrec->delete_jacobians_after_writes();
//		}
		write_solution.stop();

		/*
		 * output results to file
		 */

		visualization.start();
		if (timeprops->ifoutput() && OUTPUT) {
			move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops);

			output_discharge(matprops, timeprops, discharge, myid);

			if (myid == 0 && timeprops->verbose)
			output_summary(timeprops, statprops, savefileflag);

			write_xdmf(El_Table,Node_Table, timeprops, matprops, mapname ,XDMF_OLD);

		}
		visualization.stop();

		if (timeprops->ifsave() && (propctx.runcond & RECORD)) {
			move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops);
//			meshctx.snapshot_vec->push_back(Snapshot(meshctx,propctx));
			solrec->update_first_sol_time(timeprops->iter);
			save_forward(meshctx, propctx, solrec);
//			solrec->wrtie_sol_to_disk(myid);
			solrec->delete_jacobians_after_writes();
			if (propctx.runcond & RECORD)
				meshctx.snapshot_vec->push_back(Snapshot(meshctx,propctx));
		}
	}

	move_data(numprocs, myid, El_Table, Node_Table, timeprops, matprops);

	output_discharge(matprops, timeprops, discharge, myid);
	MPI_Barrier(MPI_COMM_WORLD);

	if (timeprops->verbose) {
		if (myid == 0)
		output_summary(timeprops, statprops, savefileflag);

		MPI_Barrier(MPI_COMM_WORLD);
		// write out ending warning, maybe flow hasn't finished moving
		sim_end_warning(El_Table, matprops, timeprops, statprops->vstar);
		MPI_Barrier(MPI_COMM_WORLD);

		//write out the final pile statistics (and run time)
		if (myid == 0)
		out_final_stats(timeprops, statprops);

		MPI_Barrier(MPI_COMM_WORLD);

		//write out stochastic simulation statistics
		//if(statprops.lhs.runid>=0)
		if (myid == 0)
		output_stoch_stats(matprops, statprops);
		MPI_Barrier(MPI_COMM_WORLD);

		//output maximum flow depth a.k.a. flow outline

		OutLine outline2;
		double dxy[2];
		dxy[0] = outline->dx;
		dxy[1] = outline->dy;
		outline2.init2(dxy, outline->xminmax, outline->yminmax);
		int NxNyout = outline->Nx * outline->Ny;
		MPI_Reduce(*(outline->pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE,
				MPI_SUM, 0, MPI_COMM_WORLD);
		if (myid == 0)
		outline2.output(matprops, statprops);

		// we deallocate these to make more space in memory
		// adjoint solution
		outline->dealloc();
		outline2.dealloc();
		MPI_Barrier(MPI_COMM_WORLD);
	}

	primal.stop();
}
