/*
 * dualmesh.h
 *
 *  Created on: May 25, 2015
 *      Author: haghakha
 */

#ifndef DUALMESH_H
#define DUALMESH_H

#include "element2.h"
#include "jacobian.h"

class SolRec: public HashTable {

private:

	// this is a solution class for void region to avoid memmory allocation for each of the void region
	Solution* solution_zero;

	// this integer shows the first time step that its solution is vailable in SolRec
	int first_solution_time_step;

	// this integer shows the last time step that its solution is vailable in SolRec
	int last_solution_time_step;

	int readflag;

	int writeflag;

	const int range;

public:

	// constructor
	SolRec(double *doublekeyrangein, int size, int prime, double XR[], double YR[], int ifrestart);

	// this function returns he pointer of solution_zero
	Solution* get_zero_solution();

	// this function records the solution of last time step
	void record_solution(MeshCTX* meshctx, PropCTX* propctx);

	// this function writes the recorded solution to files, and for each time step separately
	void wrtie_sol_to_disk();

	// this function reads the recorded solution from files, and for each time step separately and store it in SolRec
	void read_sol_from_disk(int iter);

	// this functions deletes the SolRec and deallocate the allocated memory
	void delete_empty_jacobians();

	// this function changes first_solution_time_step to the next time step
	void update_first_sol_time(int iter);

	// this function changes last_solution_time_step to the new time step
	void update_last_sol_time(int iter);

	int get_first_solution();

	int get_last_solution();

	void load_new_set_of_solution();

	void free_all_available_sol();

	int data_range();

	int write_sol();

	int read_sol();

	using HashTable::lookup;

	Solution* lookup(unsigned* key,int iter);

};

#endif /* DUALMESH_H */
