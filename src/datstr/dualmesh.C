/*
 * dualmesh.C
 *
 *  Created on: May 25, 2015
 *      Author: haghakha
 */

#include "../header/hpfem.h"

SolRec::SolRec(double *doublekeyrangein, int size, int prime, double XR[], double YR[],
		int ifrestart) :
		HashTable(doublekeyrangein, size, prime, XR, YR, ifrestart) {

	double zero_sol[3] = { 0., 0., 0. };
	double kact_z = 0.;
	solution_zero = new Solution(zero_sol, kact_z);
	first_solution_time_step = 0;
	last_solution_time_step = 0;
	readflag = 0;
	writeflag = 0;

}

Solution* SolRec::get_zero_solution() {
	return solution_zero;
}

void SolRec::record_solution(MeshCTX* meshctx, PropCTX* propctx) {

	HashTable* El_Table = meshctx->el_table;

	TimeProps* timeptr = propctx->timeprops;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					int aa = 0, bb = 1;
					unsigned keyy[2] = { 4041453883, 330382100 };
					if (*(Curr_El->pass_key()) == keyy[0] && *(Curr_El->pass_key() + 1) == keyy[1]
							&& timeptr->iter == 140)
						bb = aa;

					Jacobian *jacobian = (Jacobian *) lookup(Curr_El->pass_key());
					if (jacobian) {
						if (*(Curr_El->get_prev_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
									*(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter - 1);

						} else
							jacobian->put_solution(get_zero_solution(), timeptr->iter - 1);

					} else {
						jacobian = new Jacobian(Curr_El->pass_key());
						add(jacobian->get_key(), jacobian);

						if (*(Curr_El->get_prev_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
									*(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter - 1);

						} else
							jacobian->put_solution(get_zero_solution(), timeptr->iter - 1);

					}
				}
				currentPtr = currentPtr->next;
			}
		}

	last_solution_time_step = timeptr->iter;

}

void SolRec::wrtie_sol_to_disk() {

	FILE *myfile;
	char filename[50];
	HashEntryPtr currentPtr;
	double *solution, kact;

	for (int step = first_solution_time_step; step < last_solution_time_step; ++step) {
		int count = 0;

		sprintf(filename, "solution_%08d", step);
		myfile = fopen(filename, "w");

		for (int i = 0; i < NBUCKETS; ++i)
			if (*(bucket + i)) {

				currentPtr = *(bucket + i);
				while (currentPtr) {

					Jacobian* jacobian = (Jacobian*) (currentPtr->value);
					Solution* sol = jacobian->get_solution(step);
					count++;

//					unsigned key[2]={541694361,2576980377};
//					if (key[0]==*(jacobian->get_key()) && key[1]==*(jacobian->get_key()+1))
//						cout<<"problem found \n";

//					becasue some of the elements are created and deleted in refinement and unrefinement
					if (sol) {

						solution = sol->get_solution();
						kact = sol->get_kact();

						fprintf(myfile, "%u %u %.8e %.8e %.8e %.8e\n", *(jacobian->get_key()),
								*(jacobian->get_key() + 1), solution[0], solution[1], solution[2], sol->get_kact());

						if (sol != solution_zero)
							jacobian->erase_solution(step);
					}

					currentPtr = currentPtr->next;

				}
			}

		fclose(myfile);
		cout << "number of written elem  " << count << endl;
	}
	first_solution_time_step = last_solution_time_step;

}

void SolRec::delete_empty_jacobians() {

	HashEntryPtr currentPtr;

	for (int i = 0; i < NBUCKETS; ++i)
		if (*(bucket + i)) {

			currentPtr = *(bucket + i);
			while (currentPtr) {

				Jacobian* jacobian = (Jacobian*) (currentPtr->value);
				currentPtr = currentPtr->next;
//				if (jacobian->is_container_empty()) {

//					unsigned key[2] = { 3780766798, 3303820997 };
//					if (key[0] == *(jacobian->get_key()) && key[1] == *(jacobian->get_key() + 1))
//						cout << "problem found \n";
				this->remove(jacobian->get_key());
				delete jacobian;

//				}

			}
		}
}

void SolRec::update_first_sol_time(int iter) {

	first_solution_time_step = iter;
}

void SolRec::update_last_sol_time(int iter) {

	last_solution_time_step = iter;
}

void SolRec::read_sol_from_disk(int iter) {

	FILE *myfile;
	char filename[50];
	double state_vars[NUM_STATE_VARS] = { 0., 0., 0. };
	double kact = 0.;
	unsigned key[DIMENSION] = { 0, 0 };
	Solution * solution;
	int dbg, count = 0;

	sprintf(filename, "solution_%08d", iter);
	myfile = fopen(filename, "r");

	while (!feof(myfile)) {

		count++;

		dbg = fscanf(myfile, "%u %u %le %le %le %le\n", key, key + 1, state_vars, state_vars + 1,
				state_vars + 2, &kact);

		Jacobian *jacobian = (Jacobian *) lookup(key);

		if (state_vars[0] > 0.)
			solution = new Solution(state_vars, kact);
		else
			solution = solution_zero;

		if (jacobian)

			jacobian->put_solution(solution, iter);

		else {

			jacobian = new Jacobian(key);
			add(key, jacobian);

			jacobian->put_solution(solution, iter);

		}

	}

	fclose(myfile);
	cout << "number of readed elem  " << count << endl;
}

int SolRec::get_first_solution() {
	return first_solution_time_step;
}

int SolRec::get_last_solution() {
	return last_solution_time_step;
}

void SolRec::free_all_available_sol() {

	HashEntryPtr currentPtr;

	for (int i = 0; i < NBUCKETS; ++i)
		if (*(bucket + i)) {

			currentPtr = *(bucket + i);
			while (currentPtr) {

				Jacobian* jacobian = (Jacobian*) (currentPtr->value);

				jacobian->clear_container();

				currentPtr = currentPtr->next;

			}
		}

}

int SolRec::data_range() {
	return last_solution_time_step - first_solution_time_step;
}

int SolRec::write_sol() {
	if (data_range() > range)
		return 1;

	return 0;
}

int SolRec::read_sol() {
	if (data_range() < range)
		return 1;

	return 0;
}

void SolRec::load_new_set_of_solution() {

	struct sysinfo memInfo;
	double ratio = 1.;
	unsigned long totalPhysMem, freeram;
	last_solution_time_step = first_solution_time_step - 1;

	while ((data_range() < 2 || read_sol()) && first_solution_time_step) {

		read_sol_from_disk(first_solution_time_step - 1);
		first_solution_time_step--;
		sysinfo(&memInfo);
		totalPhysMem = memInfo.totalram;
		freeram = memInfo.freeram;
		ratio = (((long double) freeram) / ((long double) totalPhysMem));

	}

}
