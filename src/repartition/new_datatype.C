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
 * $Id: new_datatype.C 135 2007-06-07 20:15:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "../header/refined_neighbor_info.h"
#include "./repartition_BSFC.h"

MPI_Datatype ELEMTYPE;
MPI_Datatype DUALELEMTYPE;
MPI_Datatype JACTYPE;
MPI_Datatype TRANSKEYS;
MPI_Datatype NEIGHTYPE;
MPI_Datatype REFINED_INFO;
MPI_Datatype ENRICHED_INFO;
MPI_Datatype NSOLTYPE;
MPI_Datatype LB_VERT_TYPE;

struct Pointer_Holder {
	ElemPack* elem;
	DualElemPack* dualelem;
	JacPack* jacelem;
	NeighborPack* neigh;
	refined_neighbor_pack* fine;
	Neigh_Sol_Pack* neigh_sol;
	BSFC_VERTEX* sfc_vert_ptr;
} pholder;

void MPI_New_Datatype() {
	/* Create new MPI datatype: ElemPack type definition in struct.h
	 so structures of ElemPack and NeighborPack can be sent and received */
	/*  MPI_Datatype ELEMTYPE;
MPI_Datatype DUALELEMTYPE;
MPI_Datatype JACTYPE;
MPI_Datatype TRANSKEYS;
MPI_Datatype NEIGHTYPE;
MPI_Datatype REFINED_INFO;
MPI_Datatype ENRICHED_INFO;
MPI_Datatype NSOLTYPE;
MPI_Datatype LB_VERT_TYPE;*/

	int blockcounts[3] = { 37, 25 * KEYLENGTH, 66 };
	MPI_Datatype types[3];
	MPI_Aint displs[3];
	int d;
	ElemPack* elem = new ElemPack;
	pholder.elem = elem;

	MPI_Address(&(elem->myprocess), &displs[0]);
	MPI_Address(&(elem->key[0]), &displs[1]);
	MPI_Address(&(elem->elevation), &displs[2]);

	types[0] = MPI_INT;
	types[1] = MPI_UNSIGNED;
	types[2] = MPI_DOUBLE;

	for (d = 2; d >= 0; d--)
		displs[d] -= displs[0];

	MPI_Type_struct(3, blockcounts, displs, types, &ELEMTYPE);
	MPI_Type_commit(&ELEMTYPE);

	//create the 2nd new d_type

	int blockcounts_d[3] = { 37, 25 * KEYLENGTH, 72 };
	MPI_Datatype types_d[3];
	MPI_Aint displs_d[3];
	int dd;
	DualElemPack* dualelem = new DualElemPack;
	pholder.dualelem = dualelem;

	MPI_Address(&(dualelem->myprocess), &displs_d[0]);
	MPI_Address(&(dualelem->key[0]), &displs_d[1]);
	MPI_Address(&(dualelem->elevation), &displs_d[2]);

	types_d[0] = MPI_INT;
	types_d[1] = MPI_UNSIGNED;
	types_d[2] = MPI_DOUBLE;

	for (dd = 2; dd >= 0; dd--)
		displs_d[dd] -= displs_d[0];

	MPI_Type_struct(3, blockcounts_d, displs_d, types_d, &DUALELEMTYPE);
	MPI_Type_commit(&DUALELEMTYPE);

	//create the 3nd new d_type

	int blockcounts_j[2] = { 2 * KEYLENGTH, 9 };
	MPI_Datatype types_j[2];
	MPI_Aint displs_j[2];
	JacPack* jacelem = new JacPack;
	pholder.jacelem = jacelem;

	MPI_Address(&(jacelem->key_send[0]), &displs_j[0]);
	MPI_Address(&(jacelem->mat[0]), &displs_j[1]);

	types_j[0] = MPI_UNSIGNED;
	types_j[1] = MPI_DOUBLE;

	for (dd = 1; dd >= 0; dd--)
		displs_j[dd] -= displs_j[0];

	MPI_Type_struct(2, blockcounts_j, displs_j, types_j, &JACTYPE);
	MPI_Type_commit(&JACTYPE);

	//create the 4th new d_type

	MPI_Type_contiguous(12, MPI_UNSIGNED, &TRANSKEYS);
	MPI_Type_commit(&TRANSKEYS);

	//create the 5th new d_type

	int blockcounts2[2] = { 2, 2 * KEYLENGTH };
	MPI_Datatype types2[2];
	MPI_Aint displs2[2];

	NeighborPack* neigh = new NeighborPack;
	pholder.neigh = neigh;

	MPI_Address(&(neigh->target_proc), &displs2[0]);
	MPI_Address(&(neigh->elkey), &displs2[1]);

	types2[0] = MPI_INT;
	types2[1] = MPI_UNSIGNED;

	for (d = 1; d >= 0; d--)
		displs2[d] -= displs2[0];

	MPI_Type_struct(2, blockcounts2, displs2, types2, &NEIGHTYPE);
	MPI_Type_commit(&NEIGHTYPE);

	//create the 6th new d_type

	int blockcounts3[2] = { 1, 4 * KEYLENGTH };
	MPI_Datatype types3[2] = { MPI_INT, MPI_UNSIGNED };
	MPI_Aint displs3[2] = { 0, 0 };

	refined_neighbor_pack* fine = new refined_neighbor_pack;
	pholder.fine = fine;

	MPI_Address(&(fine->orig_gen), &displs3[0]);
	MPI_Address(&(fine->target_element), &displs3[1]);

	for (d = 1; d >= 0; d--)
		displs3[d] -= displs3[0];

	MPI_Type_struct(2, blockcounts3, displs3, types3, &REFINED_INFO);
	MPI_Type_commit(&REFINED_INFO);

	//create the 7th new d_type
	// for getting the neighbor solution in the new error estimator, when the neighbor is in diff subdomain

	int blockcounts5[3] = { 6, KEYLENGTH, 260 };
	MPI_Datatype types5[3];
	MPI_Aint displs5[3];

	Neigh_Sol_Pack* neigh_sol = new Neigh_Sol_Pack;
	pholder.neigh_sol = neigh_sol;

	MPI_Address(&(neigh_sol->nside), &displs5[0]);
	MPI_Address((neigh_sol->key), &displs5[1]);
	MPI_Address((neigh_sol->solu), &displs5[2]);

	types5[0] = MPI_INT;
	types5[1] = MPI_UNSIGNED;
	types5[2] = MPI_DOUBLE;

	for (d = 2; d >= 0; d--)
		displs5[d] -= displs5[0];

	MPI_Type_struct(3, blockcounts5, displs5, types5, &NSOLTYPE);
	MPI_Type_commit(&NSOLTYPE);

	//delete neigh_sol;  //added by acbauer 4/3/02 -- may be a bug?????

	int blockcounts6[3] = { 3, KEYLENGTH + 1, 1 };
	MPI_Datatype types6[3] = { MPI_INT, MPI_UNSIGNED, MPI_FLOAT };
	MPI_Aint displs6[3];

	BSFC_VERTEX* sfc_vert_ptr = new BSFC_VERTEX;
	pholder.sfc_vert_ptr = sfc_vert_ptr;

	MPI_Address(&(sfc_vert_ptr->destination_proc), &displs6[0]);
	MPI_Address(&(sfc_vert_ptr->sfc_key[0]), &displs6[1]);
	MPI_Address(&(sfc_vert_ptr->lb_weight), &displs6[2]);

	for (d = 2; d >= 0; d--)
		displs6[d] -= displs6[0];

	MPI_Type_struct(3, blockcounts6, displs6, types6, &LB_VERT_TYPE);
	MPI_Type_commit(&LB_VERT_TYPE);

	//New data types are created at this point

}

void free_mpi_types() {

	MPI_Type_free(&ELEMTYPE);
	MPI_Type_free(&DUALELEMTYPE);
	MPI_Type_free(&JACTYPE);
	MPI_Type_free(&NEIGHTYPE);
	MPI_Type_free(&REFINED_INFO);
	MPI_Type_free(&NSOLTYPE);
	MPI_Type_free(&LB_VERT_TYPE);
	MPI_Type_free(&TRANSKEYS);


	delete (pholder.dualelem);
	delete (pholder.elem);
	delete (pholder.fine);
	delete (pholder.jacelem);
	delete (pholder.neigh);
	delete (pholder.neigh_sol);
	delete (pholder.sfc_vert_ptr);

}
