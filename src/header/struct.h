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
 * $Id: struct.h 129 2007-06-07 19:54:28Z dkumar $ 
 */

#ifndef STRUCT_H
#define STRUCT_H

//! ElemPack is a smaller (memory spacewise) version of Element that can be sent from one processor to another via MPI calls
struct ElemPack {
	//see ../repartition/new_datatype.C blockcounts[3]={37,25*KEYLENGTH,66}
	int myprocess;                                              //  1
	int generation;                                             //  2
	int material;/*flag added by andrew*/                       //  3
	int neigh_proc[8];                                          // 11
	int neigh_gen[8];                                           // 19
	int refined;                                                // 20
	int adapted;                                                // 21
	int which_son;                                              // 22
	int new_old;                                                // 23
	int positive_x_side;                                        // 24
	int elm_loc[2];                                             // 26
	int opposite_brother_flag;                                  // 27
	int iwetnode;                                               // 28
	int n_info[9];                                              // 37

	unsigned key[KEYLENGTH];/*contains the 9th node key*/   //  1
	unsigned node_key[8][KEYLENGTH];                        //  9
	unsigned neighbor[8][KEYLENGTH];                        // 17
	unsigned son[4][KEYLENGTH];                             // 21
	unsigned brothers[4][KEYLENGTH];                        // 25 * KEYLENGTH

	/*NODE DEF*/
	double elevation;                                   //  1
	double n_coord[9][2];                               // 19 
	double el_error[EQUATIONS];                         // 21
	double state_vars[NUM_STATE_VARS];                  // 24
	double prev_state_vars[NUM_STATE_VARS];             // 27
	double d_state_vars[NUM_STATE_VARS * DIMENSION];    // 33
	double shortspeed;                                  // 34
	double dx[DIMENSION];                               // 36
	double eigenvxymax[DIMENSION];                      // 38
	double kactxy[DIMENSION];                           // 40
	double zeta[DIMENSION];                             // 42
	double curvature[DIMENSION];                        // 44
	double gravity[3];                                  // 47
	double d_gravity[DIMENSION];                        // 49
	double lb_weight;                                   // 50
	double node_elevation[9];                           // 59
	double Influx[NUM_STATE_VARS];                      // 62
	double Awet;                                        // 63
	double Swet;                                        // 64
	double drypoint[DIMENSION];                         // 66

};

// this structures is exactly same as ElemPack but has some more data
struct DualElemPack {
	//see ../repartition/new_datatype.C blockcounts[3]={37,25*KEYLENGTH,72}
	int myprocess;                                              //  1
	int generation;                                             //  2
	int material;/*flag added by andrew*/                       //  3
	int neigh_proc[8];                                          // 11
	int neigh_gen[8];                                           // 19
	int refined;                                                // 20
	int adapted;                                                // 21
	int which_son;                                              // 22
	int new_old;                                                // 23
	int positive_x_side;                                        // 24
	int elm_loc[2];                                             // 26
	int opposite_brother_flag;                                  // 27
	int iwetnode;                                               // 28
	int n_info[9];                                              // 37

	unsigned key[KEYLENGTH];/*contains the 9th node key*/   //  1
	unsigned node_key[8][KEYLENGTH];                        //  9
	unsigned neighbor[8][KEYLENGTH];                        // 17
	unsigned son[4][KEYLENGTH];                             // 21
	unsigned brothers[4][KEYLENGTH];                        // 25 * KEYLENGTH

	/*NODE DEF*/
	double elevation;                                   //  1
	double n_coord[9][2];                               // 19
	double el_error[EQUATIONS];                         // 21
	double state_vars[NUM_STATE_VARS];                  // 24
	double prev_state_vars[NUM_STATE_VARS];             // 27
	double adjoint[NUM_STATE_VARS];                     // 30
	double prev_adjoint[NUM_STATE_VARS];                // 33
	double d_state_vars[NUM_STATE_VARS * DIMENSION];    // 39
	double shortspeed;                                  // 40
	double dx[DIMENSION];                               // 42
	double eigenvxymax[DIMENSION];                      // 44
	double kactxy[DIMENSION];                           // 46
	double zeta[DIMENSION];                             // 48
	double curvature[DIMENSION];                        // 50
	double gravity[NUM_STATE_VARS];                     // 53
	double d_gravity[DIMENSION];                        // 55
	double lb_weight;                                   // 56
	double node_elevation[9];                           // 65
	double Influx[NUM_STATE_VARS];                      // 68
	double Awet;                                        // 69
	double Swet;                                        // 70
	double drypoint[DIMENSION];                         // 72

};

// this structures is exactly same as ElemPack but has some more data
struct ErrElemPack {
	//see ../repartition/new_datatype.C blockcounts[3]={37,25*KEYLENGTH,78}
	int myprocess;                                              //  1
	int generation;                                             //  2
	int material;/*flag added by andrew*/                       //  3
	int neigh_proc[8];                                          // 11
	int neigh_gen[8];                                           // 19
	int refined;                                                // 20
	int adapted;                                                // 21
	int which_son;                                              // 22
	int new_old;                                                // 23
	int positive_x_side;                                        // 24
	int elm_loc[2];                                             // 26
	int opposite_brother_flag;                                  // 27
	int iwetnode;                                               // 28
	int n_info[9];                                              // 37

	unsigned key[KEYLENGTH];/*contains the 9th node key*/   //  1
	unsigned node_key[8][KEYLENGTH];                        //  9
	unsigned neighbor[8][KEYLENGTH];                        // 17
	unsigned son[4][KEYLENGTH];                             // 21
	unsigned brothers[4][KEYLENGTH];                        // 25 * KEYLENGTH

	/*NODE DEF*/
	double elevation;                                   //  1
	double n_coord[9][2];                               // 19
	double el_error[EQUATIONS];                         // 21
	double state_vars[NUM_STATE_VARS];                  // 24
	double prev_state_vars[NUM_STATE_VARS];             // 27
	double adjoint[NUM_STATE_VARS];                     // 30
	double d_state_vars[NUM_STATE_VARS * DIMENSION];    // 36
	double shortspeed;                                  // 37
	double dx[DIMENSION];                               // 39
	double eigenvxymax[DIMENSION];                      // 41
	double kactxy[DIMENSION];                           // 43
	double zeta[DIMENSION];                             // 45
	double curvature[DIMENSION];                        // 47
	double gravity[NUM_STATE_VARS];                     // 50
	double d_gravity[DIMENSION];                        // 52
	double lb_weight;                                   // 53
	double node_elevation[9];                           // 62
	double Influx[NUM_STATE_VARS];                      // 65
	double Awet;                                        // 66
	double Swet;                                        // 67
	double drypoint[DIMENSION];                         // 69
	double bilin_adj[NUM_STATE_VARS];                   // 72
	double bilin_state[NUM_STATE_VARS];                 // 75
	double bilin_prev_state[NUM_STATE_VARS];            // 78
};

//see ../repartition/new_datatype.C blockcounts[2]={2*KEYLENGTH,9}
struct JacPack {
	unsigned key_send[KEYLENGTH];
	unsigned key_neigh[KEYLENGTH];
	double mat[9];
};

//see ../repartition/new_datatype.C blockcounts={6*KEYLENGTH}
struct TRANSKEY {

	unsigned key[12];

	TRANSKEY() {
	}
	;

	TRANSKEY(unsigned* keys) {
		for (unsigned i = 0; i < 12; ++i)
			key[i] = keys[i];
	}

	TRANSKEY(const TRANSKEY& transkey) {
		for (unsigned i = 0; i < 12; ++i)
			key[i] = transkey.key[i];
	}
};
//                             \|||/    
//                             (o o)   
//---------Elementlink------oo0-(_)-0oo----------STARTS HERE-------

struct ElementLink {
	int target_proc;
	int new_proc;
	unsigned elkey[KEYLENGTH];
	unsigned targetkey[KEYLENGTH];

	ElementLink* pre;
	ElementLink* next;

	ElementLink(unsigned* keyi, unsigned* key2, int tp, int np) {
		target_proc = tp;
		new_proc = np;
		int i;
		for (i = 0; i < KEYLENGTH; i++) {
			elkey[i] = *(keyi + i);
			targetkey[i] = *(key2 + i);
		}
		next = NULL;
		pre = NULL;
	}
	ElementLink() {
		next = NULL;
		pre = NULL;
	}
	~ElementLink() {
		if (next)
			next->pre = pre;
		if (pre)
			pre->next = next;
	}

};

typedef ElementLink* ELinkPtr;

//                              \|||/    
//                              (o o)   
//--------MPI datatype-------oo0-(_)-0oo----------STARTS HERE-------

struct NeighborPack {

	int target_proc;
	int new_proc;
	unsigned elkey[KEYLENGTH];
	unsigned targetkey[KEYLENGTH];

};

struct Neigh_Sol_Pack {

	int nside;
	int norder[5];
	unsigned key[KEYLENGTH];
	double solu[2][121];
	double Xnod[18];

};

typedef NeighborPack* NePtr;

#endif
