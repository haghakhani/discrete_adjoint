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
 * $Id: constant.h 129 2007-06-07 19:54:28Z dkumar $ 
 */

#ifndef CONSTANTS
#define CONSTANTS

const int KEYLENGTH = 2;
const int MAX_PROCS = 2056;

const int PRIME = 2017;        //used for creating hash key
const int DIMENSION = 2;
const int EQUATIONS = 2; /* be careful EQUATIONS and NUM_STATE_VARS
 are completely different things and are
 not interchangeable, EQUATIONS appears
 to not be used except for the el_error
 and el_solution variables (members of
 element) I believe they are only for
 FEM which we aren't really using in the
 Finite difference/volume version of
 titan.  Just leave it at 2 and otherwise
 ignore it. */

const double PI = 3.1415926;

const int CORNER = 2;//node is in corner of cells
const int BUBBLE = 6;//node is in center of cell
const int SIDE = 4;//node is in side of the cell and neighbor cell is in same generation, or when the cell is adjucent to the boundary
const int S_C_CON = 7;//node is corner of one cell and side of another cell, neighbor cells are not in same generation
const int S_S_CON = 5;//node is in side of both neighbor cells, neighbor cells are not in same generation

const int NEW = 1;
const int OLD = 0;

const int INIT = -1;
const int ONE_SIDE_NEIGH = -2;
const int MUST_BE_CORRECTED = -3;

const float LOAD_BALANCE_TOLERANCE = 1.0001;

extern void fhsfc3d(double*, unsigned*, unsigned*);
extern void fhsfc2d_(double*, unsigned*, unsigned*);

/* geoflow data */
const int NUM_STATE_VARS = 3;
const int NUM_ADJ_EQS = 3;
const int EFF_ELL = 9;

const size_t STRING_SIZE = 30;

const double GEOFLOW_TINY = 1.e-04;
const double GEOFLOW_SHORT = 3.0e-03;

const int GHOST = -9876; //"refined" GHOST CELL FLAG

//adjoint constants
const double INCREMENT = 1.49e-08;

//The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell.
#define NEWBUFFER      5  //NEWBUFFER elements have at least one current buffer element as a neighbor, a temporary marking needed for building the buffer layer, to be remarked as BUFFER elements
#define BUFFER         4  //Do not refine or unrefine elements in the buffer layer, has to be a separate flag from NEWSON because of triggered refinement
#define NEWSON         3  //Do not refine or unrefine new son elements
#define NEWFATHER      2  //Do not unrefine new father elements 
#define NOTRECADAPTED  1  //you can do anything you want to this element
#define TOBEDELETED    0  //This element has been marked for deletion, do ot involve
#define OLDFATHER     -6  //This is a temporary marking that says neighbors need to be updated with NEWSON information and then this element should be remarked as TOBEDELETED
#define OLDSON        -7  //The plan is make this analogous to OLDFATHER to allow unrefinement NEXT TO but not across interprocessor boundaries, this has not been implemented yet (date of this comment: 20061026)

typedef enum{
	FORWARD=0x01,
	ADJOINT=0x02,
	ERROR=0x04,
	RESTART=0x08,
	NORMAL=0x10,
	RECORD=0x20
}run_mode;

#endif
