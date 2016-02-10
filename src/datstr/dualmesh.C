/*
 * dualmesh.C
 *
 *  Created on: May 25, 2015
 *      Author: haghakha
 */

#include "../header/hpfem.h"

SolRec::SolRec(double *doublekeyrangein, int size, int prime, double XR[], double YR[],
    int ifrestart) :
		HashTable(doublekeyrangein, size, prime, XR, YR, ifrestart), range(50) {

	first_solution_time_step = 0;
	last_solution_time_step = 0;
	readflag = 0;
	writeflag = 0;

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

//					int aa = 0, bb = 1;
//					unsigned keyy[2] = { 4041453883, 330382100 };
//					if (*(Curr_El->pass_key()) == keyy[0] && *(Curr_El->pass_key() + 1) == keyy[1]
//							&& timeptr->iter == 140)
//						bb = aa;

					Jacobian *jacobian = (Jacobian *) lookup(Curr_El->pass_key());
					if (jacobian) {
						if (*(Curr_El->get_prev_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
							    *(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter - 1);

						} else
							jacobian->put_solution(&(Solution::solution_zero), timeptr->iter - 1);

					} else {
						jacobian = new Jacobian(Curr_El->pass_key());
						add(jacobian->get_key(), jacobian);

						if (*(Curr_El->get_prev_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
							    *(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter - 1);

						} else
							jacobian->put_solution(&(Solution::solution_zero), timeptr->iter - 1);

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

						fflush(myfile);
						//fsync(fileno(myfile));

						if (sol != &(Solution::solution_zero))
							delete sol;
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

void SolRec::delete_jacobians_after_writes() {

	HashEntryPtr currentPtr;

	for (int i = 0; i < NBUCKETS; ++i)
		if (*(bucket + i)) {

			currentPtr = *(bucket + i);
			while (currentPtr) {

				Jacobian* jacobian = (Jacobian*) (currentPtr->value);
				currentPtr = currentPtr->next;
				this->remove(jacobian->get_key());
				delete jacobian;

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
			solution = &(Solution::solution_zero);

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

				currentPtr = currentPtr->next;

				jacobian->clear_container(first_solution_time_step);

				if (jacobian->is_container_empty()) {
					this->remove(jacobian->get_key());
					delete jacobian;
				}

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

Solution* SolRec::lookup(unsigned* key, int iter) {

	Jacobian *jacobian = (Jacobian *) lookup(key);
	if (jacobian)
		return jacobian->get_solution(iter);
	else
		return NULL;

}

//==============================================================

DualElem::DualElem(Element* element) {

	myprocess = element->myprocess;

	generation = element->generation;

	opposite_brother_flag = element->opposite_brother_flag;

	material = element->material;

	lb_weight = element->lb_weight;

	refined = element->refined;

	adapted = element->adapted;

	which_son = element->which_son;

	new_old = element->new_old;

	shortspeed = element->shortspeed;

	positive_x_side = element->positive_x_side;

	elevation = element->elevation;

	stoppedflags = element->stoppedflags;

	effect_bedfrict = element->effect_bedfrict;

	effect_tanbedfrict = element->effect_tanbedfrict;

	counted = element->counted;

	ithelem = element->ithelem;

	iwetnode = element->iwetnode;

	Awet = element->Awet;

	Swet = element->Swet;

	for (int i = 0; i < NUM_STATE_VARS; ++i) {

		state_vars[i] = element->state_vars[i];

		prev_state_vars[i] = element->prev_state_vars[i];

		gravity[i] = element->gravity[i];

		Influx[i] = element->Influx[i];

		d_state_vars[i] = element->d_state_vars[i];
		d_state_vars[i + NUM_STATE_VARS] = element->d_state_vars[i + NUM_STATE_VARS];

	}

	for (int i = 0; i < DIMENSION; ++i) {

		dx[i] = element->dx[i];

		eigenvxymax[i] = element->eigenvxymax[i];

		kactxy[i] = element->kactxy[i];

		zeta[i] = element->zeta[i];

		curvature[i] = element->curvature[i];

		d_gravity[i] = element->d_gravity[i];

		effect_kactxy[i] = element->effect_kactxy[i];

		drypoint[i] = element->drypoint[i];

		coord[i] = element->coord[i];

		elm_loc[i] = element->elm_loc[i];

		lb_key[i] = element->lb_key[i];

		key[i] = element->key[i];

		father[i] = element->father[i];

		el_error[i] = element->el_error[i];
	}

	for (int i = 0; i < 8; ++i) {

		neigh_gen[i] = element->neigh_gen[i];

		neigh_proc[i] = element->neigh_proc[i];

		for (int j = 0; j < KEYLENGTH; ++j) {

			node_key[i][j] = element->node_key[i][j];

			neighbor[i][j] = element->neighbor[i][j];
		}
	}

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 2; ++j) {
			son[i][j] = element->son[i][j];
			brothers[i][j] = element->brothers[i][j];
		}

	for (int i = 0; i < NUM_STATE_VARS; ++i) {

		adjoint[i] = 0.;

		prev_adjoint[i] = 0.;

		func_sens[i] = 0.;
	}
}

//used for refinement
DualElem::DualElem(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[],
    int gen, int elm_loc_in[], int gen_neigh[], int mat, DualElem *fthTemp, double *coord_in,
    HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather,
    double Awetfather, double *drypoint_in, int SETLINK) {
	counted = 0; //for debugging only

	adapted = NEWSON;

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		prev_state_vars[i] = 0.;
		Influx[i] = 0.;
		adjoint[i] = 0.;
		prev_adjoint[i] = 0.;
		Influx[i] = 0.;
	}

	for (int i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
		d_state_vars[i] = 0.;

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	for (int i = 0; i < NUM_STATE_VARS; i++)
		func_sens[i] = 0.;

	for (int i = 0; i < 4; i++)
		son[i][0] = son[i][1] = brothers[i][0] = brothers[i][1] = NULL;
	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	myprocess = myid;
	generation = gen; //--first generation
	opposite_brother_flag = 1;
	material = mat;
	for (int i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (int i = 0; i < KEYLENGTH; i++) {
		father[i] = fthTemp->key[i];
		key[i] = nodekeys[8][i]; //--using buble key to represent the element
	}

//	int aa = 0, bb = 1;
//	unsigned keyy[2] = { 3410598297, 2576980374 };
//	if (key[0] == keyy[0] && key[1] == keyy[1])
//		bb = aa;

	elm_loc[0] = elm_loc_in[0];
	elm_loc[1] = elm_loc_in[1];

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < KEYLENGTH; j++)
			node_key[i][j] = nodekeys[i][j];

	for (int i = 0; i < 4; i++) {
		neigh_proc[i] = n_pro[i];
		neigh_proc[i + 4] = -2; //-- -2 means regular element
		if (neigh_proc[i] != -1)
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = neigh[i][j];
		else
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = NULL;

		neigh_gen[i] = neigh_gen[i + 4] = gen_neigh[i];
	}

	refined = 0;

	new_old = NEW;
	//geoflow stuff
	dx[0] = .5 * fthTemp->dx[0];  //assume constant refinement
	dx[1] = .5 * fthTemp->dx[1];

	iwetnode = iwetnodefather;
	drypoint[0] = drypoint_in[0];
	drypoint[1] = drypoint_in[1];

	double myfractionoffather;
	if ((Awetfather == 0.0) || (Awetfather == 1.0)) {
		Awet = Awetfather;
		myfractionoffather = 1.0;
	} else {
		Awet = convect_dryline(dx, 0.0); //dx is a dummy stand in for convection speed... value doesn't matter because it's being multiplied by a timestep of zero
		myfractionoffather = Awet / Awetfather;
	}
	Swet = 1.0;

	double dxx = coord_in[0] - fthTemp->coord[0];
	double dyy = coord_in[1] - fthTemp->coord[1];

	if (state_vars[0] < 0.)
		state_vars[0] = 0.;

	find_positive_x_side(NodeTable);
	calc_topo_data(matprops_ptr);
	calc_gravity_vector(matprops_ptr);
	calc_d_gravity(El_Table);

	coord[0] = coord_in[0];
	coord[1] = coord_in[1];

	calc_which_son();

	stoppedflags = fthTemp->stoppedflags;

	kactxy[0] = fthTemp->kactxy[0];
	kactxy[1] = fthTemp->kactxy[1];
	effect_kactxy[0] = fthTemp->effect_kactxy[0];
	effect_kactxy[1] = fthTemp->effect_kactxy[1];
	for (int i = 0; i < NUM_STATE_VARS; i++) {
		state_vars[i] = fthTemp->state_vars[i] * myfractionoffather;
		prev_state_vars[i] = fthTemp->prev_state_vars[i] * myfractionoffather;
		adjoint[i] = 0.25 * fthTemp->adjoint[i];
		prev_adjoint[i] = 0.25 * fthTemp->prev_adjoint[i];
		Influx[i] = 0.;

	}

	if (SETLINK) {

		// correcting the links
		vector<ErrorElem*>& errel_vec = fthTemp->get_son_addresses();

		for (int j = 0; j < errel_vec.size(); ++j) {

//			assert(errel_vec[j]->get_father_address() == fthTemp);

			unsigned* father_key = errel_vec[j]->getfather();

			if (compare_key(father_key, key)) {
				errel_vec[j]->put_father_address(this);
				son_address.push_back(errel_vec[j]);

			}
		}
	}

}

/*********************************
 making a father element from its sons
 *****************************************/
DualElem::DualElem(DualElem* sons[], HashTable* NodeTable, HashTable* El_Table,
    MatProps* matprops_ptr, int SETLINK) {
	counted = 0; //for debugging only

	adapted = NEWFATHER;

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		prev_state_vars[i] = 0.;
		Influx[i] = 0.;
		adjoint[i] = 0.;
		prev_adjoint[i] = 0.;
		Influx[i] = 0.;
		d_state_vars[i] = d_state_vars[NUM_STATE_VARS + i] = 0.;
		func_sens[i] = 0.;
	}

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	int i, j, ikey, ison, isonneigh, ineigh;

	for (ikey = 0; ikey < KEYLENGTH; ikey++)
		key[ikey] = *(sons[2]->getNode() + ikey);

	for (ison = 0; ison < 4; ison++) {
		sons[ison]->put_adapted_flag(OLDSON);
		for (ikey = 0; ikey < KEYLENGTH; ikey++) {
			son[ison][ikey] = *(sons[ison]->pass_key() + ikey);
			sons[ison]->put_father(key);
		}
	}

	//for(i=0;i<4;i++) son[i][0]=son[i][1]=brothers[i][0]=brothers[i][1]=NULL;
	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	new_old = NEW;
	unsigned* son_nodes[4];
	opposite_brother_flag = 0;
	stoppedflags = 2;
	for (i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (ison = 0; ison < 4; ison++) {
		son_nodes[ison] = sons[ison]->getNode();
		if (sons[ison]->stoppedflags < stoppedflags)
			stoppedflags = sons[ison]->stoppedflags;
	}

	for (ikey = 0; ikey < KEYLENGTH; ikey++) {
		father[ikey] = NULL;
		node_key[0][ikey] = son_nodes[0][ikey];
		node_key[1][ikey] = son_nodes[1][KEYLENGTH + ikey];
		node_key[2][ikey] = son_nodes[2][2 * KEYLENGTH + ikey];
		node_key[3][ikey] = son_nodes[3][3 * KEYLENGTH + ikey];
		node_key[4][ikey] = son_nodes[0][KEYLENGTH + ikey];
		node_key[5][ikey] = son_nodes[1][2 * KEYLENGTH + ikey];
		node_key[6][ikey] = son_nodes[2][3 * KEYLENGTH + ikey];
		node_key[7][ikey] = son_nodes[3][ikey];

		elm_loc[ikey] = (*(sons[0]->get_elm_loc() + ikey)) / 2;
	}
	myprocess = sons[0]->get_myprocess();
	generation = sons[0]->get_gen() - 1;
	material = sons[0]->get_material();

	calc_which_son();

	refined = 1; // not an active element yet!!!

	// neighbor information
	for (ison = 0; ison < 4; ison++) {
		isonneigh = ison;
		ineigh = isonneigh;
		neigh_gen[ineigh] = *(sons[ison]->get_neigh_gen() + isonneigh);
		for (ikey = 0; ikey < KEYLENGTH; ikey++)
			neighbor[ineigh][ikey] = *(sons[ison]->get_neighbors() + isonneigh * KEYLENGTH + ikey);
		neigh_proc[ineigh] = *(sons[ison]->get_neigh_proc() + isonneigh);

		isonneigh = (ison + 3) % 4;
		ineigh = isonneigh + 4;
		neigh_gen[ineigh] = *(sons[ison]->get_neigh_gen() + isonneigh);
		for (ikey = 0; ikey < KEYLENGTH; ikey++)
			neighbor[ineigh][ikey] = *(sons[ison]->get_neighbors() + isonneigh * KEYLENGTH + ikey);
		if ((*(sons[ison]->get_neigh_gen() + isonneigh) == generation)
		    || (*(sons[ison]->get_neigh_proc() + isonneigh) == -1))
			neigh_proc[ineigh] = -2;
		else
			neigh_proc[ineigh] = *(sons[ison]->get_neigh_proc() + isonneigh);
	}

	/* brother information -- requires that at least one of this
	 element's neighboring brothers is on this process in
	 order to get information on the brother that is not a neighbor */
	DualElem* EmTemp;
	switch (which_son) {
		case 0:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[0][i] = key[i];
			if (neigh_proc[1] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = 0;
			} else if (neigh_gen[1] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = neighbor[1][i];
			} else if (neigh_gen[1] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[1]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = bro_key[i];
			} else if (neigh_gen[1] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[2] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = 0;
			} else if (neigh_gen[2] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = neighbor[2][i];
			} else if (neigh_gen[2] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[2]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = bro_key[i];
			} else if (neigh_gen[2] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 1:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[1][i] = key[i];
			if (neigh_proc[3] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = 0;
			} else if (neigh_gen[3] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = neighbor[3][i];
			} else if (neigh_gen[3] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[3]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = bro_key[i];
			} else if (neigh_gen[3] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[2] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = 0;
			} else if (neigh_gen[2] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = neighbor[2][i];
			} else if (neigh_gen[2] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[2]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = bro_key[i];
			} else if (neigh_gen[2] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 2:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[2][i] = key[i];
			if (neigh_proc[0] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = 0;
			} else if (neigh_gen[0] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = neighbor[0][i];
			} else if (neigh_gen[0] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[0]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = bro_key[i];
			} else if (neigh_gen[0] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[3] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = 0;
			} else if (neigh_gen[3] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = neighbor[3][i];
			} else if (neigh_gen[3] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[3]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = bro_key[i];
			} else if (neigh_gen[3] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 3:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[3][i] = key[i];
			if (neigh_proc[0] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = 0;
			} else if (neigh_gen[0] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = neighbor[0][i];
			} else if (neigh_gen[0] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[0]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = bro_key[i];
			} else if (neigh_gen[0] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[1] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = 0;
			} else if (neigh_gen[1] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = neighbor[1][i];
			} else if (neigh_gen[1] == generation + 1) {
				EmTemp = (DualElem*) El_Table->lookup(neighbor[1]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = bro_key[i];
			} else if (neigh_gen[1] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
	}

	find_positive_x_side(NodeTable);  //also inserts the coordinates
	calculate_dx(NodeTable);
	find_opposite_brother(El_Table);
	calc_topo_data(matprops_ptr);
	calc_gravity_vector(matprops_ptr);
	calc_d_gravity(El_Table);

	for (i = 0; i < NUM_STATE_VARS; i++) {
		state_vars[i] = 0.;
		prev_state_vars[i] = 0.;
		adjoint[i] = 0.;
		prev_adjoint[i] = 0.;

		for (j = 0; j < 4; j++) {
			state_vars[i] += *(sons[j]->get_state_vars() + i) * 0.25;
			prev_state_vars[i] += *(sons[j]->get_prev_state_vars() + i) * 0.25;
			adjoint[i] += *(sons[j]->get_adjoint() + i);
			prev_adjoint[i] += *(sons[j]->get_prev_adjoint() + i);

		}
	}

	Awet = 0.0;
	for (int ison = 0; ison < 4; ison++)
		Awet += sons[ison]->get_Awet();
	Awet *= 0.25;

	//uninitialized flag values... will fix shortly
	drypoint[0] = drypoint[1] = 0.0;
	iwetnode = 8;
	Swet = 1.0;

	//calculate the shortspeed
	shortspeed = 0.0;
	for (j = 0; j < 4; j++)
		shortspeed += *(sons[j]->get_state_vars()) * sons[j]->get_shortspeed();
	if (state_vars[0] > 0.0)
		shortspeed /= (4.0 * state_vars[0]);

	kactxy[0] = effect_kactxy[0] = 0.0;
	kactxy[1] = effect_kactxy[1] = 0.0;

	for (j = 0; j < 4; j++) {
//		for (i = 0; i < EQUATIONS; i++) {
//			kactxy[i] += *(sons[j]->get_kactxy() + i) * 0.25;
//			effect_kactxy[i] += *(sons[j]->get_effect_kactxy() + i) * 0.25;
//			el_error[i] += *(sons[j]->get_el_error() + i) * 0.25;
//		}

	}
	if (SETLINK) {
		// correcting the links
		for (int i = 0; i < 4; ++i) {
			vector<ErrorElem*>& errel_vec = sons[i]->get_son_addresses();

			for (int j = 0; j < errel_vec.size(); ++j) {

				assert(errel_vec[j]->get_father_address() == sons[i]);
				errel_vec[j]->put_father_address(this);
				son_address.push_back(errel_vec[j]);

			}
		}
	}

}

void DualElem::get_slopes_prev(HashTable* El_Table, HashTable* NodeTable, double gamma) {
	int j = 0, bc = 0;
	/* check to see if this is a boundary */
	while (j < 4 && bc == 0) {
		if (neigh_proc[j] == INIT)
			bc = 1;
		j++;
	}
	if (bc == 1) {
		for (j = 0; j < NUM_STATE_VARS * DIMENSION; j++)
			d_state_vars[j] = 0;

		for (int i = 0; i < DIMENSION; ++i)
			for (int k = 0; k < 5; ++k)
				hslope_sens(i, k) = 0.;

		return;
	}

	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}

	/* x direction */
	DualElem *ep = (DualElem*) (El_Table->lookup(&neighbor[xp][0]));
	DualElem *em = (DualElem*) (El_Table->lookup(&neighbor[xm][0]));
	DualElem *ep2 = NULL;
	DualElem *em2 = NULL;

//check if element has 2 neighbors on either side
	Node* ndtemp = (Node*) NodeTable->lookup(&node_key[xp + 4][0]);

	if (ndtemp->info == S_C_CON) {
		ep2 = (DualElem*) (El_Table->lookup(&neighbor[xp + 4][0]));
		assert(neigh_proc[xp + 4] >= 0 && ep2);
	}

	ndtemp = (Node*) NodeTable->lookup(&node_key[xm + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (DualElem*) (El_Table->lookup(&neighbor[xm + 4][0]));
		assert(neigh_proc[xm + 4] >= 0 && em2);
	}

	double dp, dm, dc, dxp, dxm;
	dxp = ep->coord[0] - coord[0];
	dxm = coord[0] - em->coord[0];

	double inv_dxp = 1 / dxp, inv_dxm = 1 / dxm, inv_sdxpdxm = 1 / (dxm + dxp);

	double min_slopes;

	for (j = 0; j < NUM_STATE_VARS; j++) {

		dp = (ep->prev_state_vars[j] - prev_state_vars[j]) * inv_dxp;

		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->prev_state_vars[j] - prev_state_vars[j]) * inv_dxp);

		dm = (prev_state_vars[j] - em->prev_state_vars[j]) * inv_dxm;

		if (em2 != NULL)
			dm = .5 * (dm + (prev_state_vars[j] - em2->prev_state_vars[j]) * inv_dxm);

		dc = (dp * dxm + dm * dxp) * inv_sdxpdxm;  // weighted average

		min_slopes = c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));

		//do slope limiting
		d_state_vars[j] = .5 * (c_sgn(dp) + c_sgn(dm)) * min_slopes;

		if (j == 0)
			// if dp and dm have different signes
			if (.5 * (c_sgn(dp) + c_sgn(dm)) == 0)
				for (int i = 0; i < 5; ++i)
					hslope_sens(0, i) = 0.;

			else if (min_slopes == gamma * dabs(dp)) {

				hslope_sens(0, 0) = -inv_dxp;
				hslope_sens(0, 1) = 0.;
				hslope_sens(0, 2) = inv_dxp;
				hslope_sens(0, 3) = 0.;
				hslope_sens(0, 4) = 0.;

				if (ep2 != NULL)
					hslope_sens(0, 2) = hslope_sens(0, 4) = .5 * inv_dxp;

			} else if (min_slopes == gamma * dabs(dm)) {

				hslope_sens(0, 0) = inv_dxm;
				hslope_sens(0, 1) = -inv_dxm;
				hslope_sens(0, 2) = 0.;
				hslope_sens(0, 3) = 0.;
				hslope_sens(0, 4) = 0.;

				if (em2 != NULL)
					hslope_sens(0, 1) = hslope_sens(0, 3) = -.5 * inv_dxm;

			} else {

				// under any circumstances hslope_sens[0][0] is this
				hslope_sens(0, 0) = inv_dxm - inv_dxp;
				hslope_sens(0, 1) = -inv_dxm;
				hslope_sens(0, 2) = inv_dxp;
				hslope_sens(0, 3) = 0.;
				hslope_sens(0, 4) = 0.;

				if (ep2 != NULL)
					hslope_sens(0, 2) = hslope_sens(0, 4) = .5 * inv_dxp;

				if (em2 != NULL)
					hslope_sens(0, 1) = hslope_sens(0, 3) = -.5 * inv_dxm;

				hslope_sens(0, 1) *= dxp * inv_sdxpdxm;
				hslope_sens(0, 2) *= dxm * inv_sdxpdxm;
				hslope_sens(0, 3) *= dxp * inv_sdxpdxm;
				hslope_sens(0, 4) *= dxm * inv_sdxpdxm;

			}
	}

	/* y direction */
	ep = (DualElem*) (El_Table->lookup(&neighbor[yp][0]));
	em = (DualElem*) (El_Table->lookup(&neighbor[ym][0]));
	ep2 = NULL;
	em2 = NULL;
//check if element has 2 neighbors on either side
	ndtemp = (Node*) NodeTable->lookup(&node_key[yp + 4][0]);
	if (ndtemp->info == S_C_CON) {
		ep2 = (DualElem*) (El_Table->lookup(&neighbor[yp + 4][0]));
		assert(neigh_proc[yp + 4] >= 0 && ep2);
	}
	ndtemp = (Node*) NodeTable->lookup(&node_key[ym + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (DualElem*) (El_Table->lookup(&neighbor[ym + 4][0]));
		if (!(neigh_proc[ym + 4] >= 0 && em2)) {
			printf("ym=%d neigh_proc[ym+4]=%d em2=%d\n", ym, neigh_proc[ym + 4], em2);
		}
		assert(neigh_proc[ym + 4] >= 0 && em2);
	}

	dxp = ep->coord[1] - coord[1];
	dxm = coord[1] - em->coord[1];

	inv_dxp = 1 / dxp;
	inv_dxm = 1 / dxm;
	inv_sdxpdxm = 1 / (dxm + dxp);

	for (j = 0; j < NUM_STATE_VARS; j++) {
		dp = (ep->prev_state_vars[j] - prev_state_vars[j]) * inv_dxp;
		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->prev_state_vars[j] - prev_state_vars[j]) * inv_dxp);
		dm = (prev_state_vars[j] - em->prev_state_vars[j]) * inv_dxm;
		if (em2 != NULL)
			dm = .5 * (dm + (prev_state_vars[j] - em2->prev_state_vars[j]) * inv_dxm);

		dc = (dp * dxm + dm * dxp) * inv_sdxpdxm;  // weighted average

		min_slopes = c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
		//do slope limiting
		d_state_vars[j + NUM_STATE_VARS] = .5 * (c_sgn(dp) + c_sgn(dm)) * min_slopes;

		if (j == 0)
			// if dp and dm have different signes
			if (.5 * (c_sgn(dp) + c_sgn(dm)) == 0)
				for (int i = 0; i < 5; ++i)
					hslope_sens(1, i) = 0.;

			else if (min_slopes == gamma * dabs(dp)) {

				hslope_sens(1, 0) = -inv_dxp;
				hslope_sens(1, 1) = 0.;
				hslope_sens(1, 2) = inv_dxp;
				hslope_sens(1, 3) = 0.;
				hslope_sens(1, 4) = 0.;

				if (ep2 != NULL)
					hslope_sens(1, 2) = hslope_sens(1, 4) = .5 * inv_dxp;

			} else if (min_slopes == gamma * dabs(dm)) {

				hslope_sens(1, 0) = inv_dxm;
				hslope_sens(1, 1) = -inv_dxm;
				hslope_sens(1, 2) = 0.;
				hslope_sens(1, 3) = 0.;
				hslope_sens(1, 4) = 0.;

				if (em2 != NULL)
					hslope_sens(1, 1) = hslope_sens(1, 3) = -.5 * inv_dxm;

			} else {

				// under any circumstances hslope_sens[0][0] is this
				hslope_sens(1, 0) = inv_dxm - inv_dxp;
				hslope_sens(1, 1) = -inv_dxm;
				hslope_sens(1, 2) = inv_dxp;
				hslope_sens(1, 3) = 0.;
				hslope_sens(1, 4) = 0.;

				if (ep2 != NULL)
					hslope_sens(1, 2) = hslope_sens(1, 4) = .5 * inv_dxp;

				if (em2 != NULL)
					hslope_sens(1, 1) = hslope_sens(1, 3) = -.5 * inv_dxm;

				hslope_sens(1, 1) *= dxp * inv_sdxpdxm;
				hslope_sens(1, 2) *= dxm * inv_sdxpdxm;
				hslope_sens(1, 3) *= dxp * inv_sdxpdxm;
				hslope_sens(1, 4) *= dxm * inv_sdxpdxm;

			}
	}
}

void dual_riemannflux(Mat3x3& hfvl, Mat3x3& hfvr, double flux[NUM_STATE_VARS], Mat3x3& flux_jac_l,
    Mat3x3& flux_jac_r, Mat3x3& s_jac_l, Mat3x3& s_jac_r, Mat3x3& jac_l, Mat3x3& jac_r) {
//hfv: h=state variable, f=flux, v=wave speeds
//l="left" (the minus side), r="right" (the plus side)

	if ((hfvl(0, 0) == 0.) && (hfvr(0, 0) == 0.)) {
		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			flux[ivar] = 0.;

		for (int i = 0; i < NUM_STATE_VARS; ++i)
			for (int j = 0; j < NUM_STATE_VARS; ++j)
				jac_r(i, j) = jac_l(i, j) = 0.;

	} else {

		double sl, sr;
		double d_sl_l[3], d_sl_r[3], d_sr_l[3], d_sr_r[3];

		if (hfvl(0, 0) == 0.) {

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				d_sr_l[i] = d_sl_l[i] = 0.;

			sl = fmin(0., 2. * hfvr(2, 0) - hfvr(2, 1));

			if (sl == 0) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sl_r[i] = 0.;

			} else {

				// now we have to implement chain rule
				// \frac{\partial s_l}{\partial h}=
				// 2.0*\frac{\partial hfvr[2][0]}{\partial h_r}-\frac{\partial hfvr[2][1]}{\partial h_r}
				// and similarly for other terms
				d_sl_r[0] = 2 * s_jac_r(0, 0) - s_jac_r(1, 0);
				d_sl_r[1] = 2 * s_jac_r(0, 1) - s_jac_r(1, 1);
				d_sl_r[2] = 2 * s_jac_r(0, 2) - s_jac_r(1, 2);

			}

			sr = fmax(0., 2. * hfvr(2, 2) - hfvr(2, 1));

			if (sr == 0) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sr_r[i] = 0.;

			} else {

				// now we have to implement chain rule
				d_sr_r[0] = 2 * s_jac_r(2, 0) - s_jac_r(1, 0);
				d_sr_r[1] = 2 * s_jac_r(2, 1) - s_jac_r(1, 1);
				d_sr_r[2] = 2 * s_jac_r(2, 2) - s_jac_r(1, 2);

			}

		} else if (hfvr(0, 0) == 0.) {

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				d_sr_r[i] = d_sl_r[i] = 0.;

			sl = fmin(0., 2. * hfvl(2, 0) - hfvl(2, 1));

			if (sl == 0) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sl_l[i] = 0.;

			} else {

				// now we have to implement chain rule
				d_sl_l[0] = 2 * s_jac_l(0, 0) - s_jac_l(1, 0);
				d_sl_l[1] = 2 * s_jac_l(0, 1) - s_jac_l(1, 1);
				d_sl_l[2] = 2 * s_jac_l(0, 2) - s_jac_l(1, 2);

			}

			sr = fmax(0., 2. * hfvl(2, 2) - hfvl(2, 1));

			if (sr == 0) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sr_l[i] = 0.;

			} else {

				// now we have to implement chain rule
				d_sr_l[0] = 2 * s_jac_l(2, 0) - s_jac_l(1, 0);
				d_sr_l[1] = 2 * s_jac_l(2, 1) - s_jac_l(1, 1);
				d_sr_l[2] = 2 * s_jac_l(2, 2) - s_jac_l(1, 2);

			}

		} else {
			sl = fmin(0., fmin(hfvl(2, 0), hfvr(2, 0)));

			if (sl == 0) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sl_l[i] = d_sl_r[i] = 0.;

			} else if (sl == hfvl(2, 0)) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sl_r[i] = 0.;

				d_sl_l[0] = s_jac_l(0, 0);
				d_sl_l[1] = s_jac_l(0, 1);
				d_sl_l[2] = s_jac_l(0, 2);

			} else {

				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sl_l[i] = 0.;

				d_sl_r[0] = s_jac_r(0, 0);
				d_sl_r[1] = s_jac_r(0, 1);
				d_sl_r[2] = s_jac_r(0, 2);

			}
			sr = fmax(0., fmax(hfvl(2, 2), hfvr(2, 2)));

			if (sr == 0) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sr_l[i] = d_sr_r[i] = 0.;

			} else if (sr == hfvl(2, 2)) {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sr_r[i] = 0.;

				d_sr_l[0] = s_jac_l(2, 0);
				d_sr_l[1] = s_jac_l(2, 1);
				d_sr_l[2] = s_jac_l(2, 2);

			} else {
				for (int i = 0; i < NUM_STATE_VARS; ++i)
					d_sr_l[i] = 0.;

				d_sr_r[0] = s_jac_r(2, 0);
				d_sr_r[1] = s_jac_r(2, 1);
				d_sr_r[2] = s_jac_r(2, 2);

			}

		}

		if (sl >= 0.) {
			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				flux[ivar] = hfvl(1, ivar);

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				for (int j = 0; j < NUM_STATE_VARS; ++j) {
					jac_r(i, j) = 0.;
					jac_l(i, j) = flux_jac_l(i, j);
				}

		} else if (sr <= 0.) {
			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				flux[ivar] = hfvr(1, ivar);

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				for (int j = 0; j < NUM_STATE_VARS; ++j) {
					jac_r(i, j) = flux_jac_r(i, j);
					jac_l(i, j) = 0.;
				}

		} else {

			double inv_speed_diff = 1. / (sr - sl);
			double v_diff[NUM_STATE_VARS];
			double d_flux_d_sr[NUM_STATE_VARS], d_flux_d_sl[NUM_STATE_VARS];

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				v_diff[i] = hfvr(0, i) - hfvl(0, i);

			// in this formulation flux is a function of 5 variables: sr, sl, Fr, Fl, v_diff
			// then we find the derivative by chain rule
			for (int i = 0; i < NUM_STATE_VARS; ++i)
				flux[i] = (sr * hfvl(1, i) - sl * hfvr(1, i) + sl * sr * v_diff[i]) * inv_speed_diff;

			for (int i = 0; i < NUM_STATE_VARS; i++) {
				d_flux_d_sr[i] = (hfvl(1, i) + sl * v_diff[i] - flux[i]) * inv_speed_diff;
				d_flux_d_sl[i] = (-hfvr(1, i) + sr * v_diff[i] + flux[i]) * inv_speed_diff;
			}

			Mat3x3 jac_v_diff_r, jac_v_diff_l;

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				for (int j = 0; j < NUM_STATE_VARS; ++j) {

					if (i == j) {
						jac_v_diff_r(i, j) = sr * sl * inv_speed_diff;
						jac_v_diff_l(i, j) = -sr * sl * inv_speed_diff;
					} else

						jac_v_diff_r(i, j) = jac_v_diff_l(i, j) = 0.;
				}

			//---------------------------------------
			double d_flux_d_Fr = -sl * inv_speed_diff;
			double d_flux_d_Fl = sr * inv_speed_diff;
			// very interestingly I noticed that in jacobian computation I have following terms:
			// d_flux_d_Fr*jac_flux_r[i][j]*d_flux_d_Fl*jac_flux_r[i][j]
			// which can simplify to:
			// jac_flux_r[i][j]

			for (int i = 0; i < NUM_STATE_VARS; ++i)
				for (int j = 0; j < NUM_STATE_VARS; ++j) {

					jac_r(i, j) = d_flux_d_sr[i] * d_sr_r[j] + d_flux_d_sl[i] * d_sl_r[j]
					    + d_flux_d_Fr * flux_jac_r(i, j) + jac_v_diff_r(i, j);

					jac_l(i, j) = d_flux_d_sr[i] * d_sr_l[j] + d_flux_d_sl[i] * d_sl_l[j]
					    + d_flux_d_Fl * flux_jac_l(i, j) + jac_v_diff_l(i, j);
				}

		}
	}

}

void DualElem::calc_fluxes(HashTable* El_Table, HashTable* NodeTable,
    vector<DualElem*>* x_elem_list, vector<DualElem*>* y_elem_list, int myid) {

	calc_flux(El_Table, NodeTable, x_elem_list, myid, 0);
	calc_flux(El_Table, NodeTable, y_elem_list, myid, 1);

}

void DualElem::calc_flux(HashTable* El_Table, HashTable* NodeTable, vector<DualElem*>* elem_list,
    int myid, int side) {

	Node *np, *np1, *np2, *nm, *nm1, *nm2;
	DualElem *elm1, *elm2;
	int zp, zm;
	int zp2, zm2; //positive_z_side_2 minus_z_side_2
	int zelmpos = -100, zelmpos_2 = -100;

	Mat3x3 hfv, hfv1, hfv2;
	Mat3x3 flux_jac, flux_jac1, flux_jac2;
	Mat3x3 s_jac, s_jac1, s_jac2;
	Mat3x3 jac, jac_neigh1, jac_neigh2, jac_res; // the latest one is for times we need to evaluate riemannflux for the current element twice

//ghost elements don't have nodes so you have to make temp storage for flux
	double ghostflux[NUM_STATE_VARS];

//  if (key[0]==3978454630 && key[1]==1717986917)
//    cout<<"this element is being checked"<<endl;

	zp = (positive_x_side + side) % 4;
	zm = (zp + 2) % 4;
	np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);

	if (neigh_proc[zp] == -1) {

		elem_list->push_back(this);

	} else if (neigh_proc[zp] != myid) {

		np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
		elm1 = (DualElem*) El_Table->lookup(&neighbor[zp][0]);
		assert(elm1);

		dual_zdirflux(side, hfv, flux_jac, s_jac);
		elm1->dual_zdirflux(side + 2, hfv1, flux_jac1, s_jac1);

		dual_riemannflux(hfv, hfv1, np->flux, flux_jac, flux_jac1, s_jac, s_jac1, jac, jac_neigh1);

		//now we can store jacobians in elements
		flx_jac_cont.set(side, 1, 0, jac);
		flx_jac_cont.set(side, 1, 1, jac_neigh1);
		flx_jac_cont.set(side, 1, 2, ZERO_MATRIX);

		// here the element elm1 must be a ghost element
		(elm1->get_flx_jac_cont()).set(side, 0, 0, jac_neigh1);
		(elm1->get_flx_jac_cont()).set(side, 0, 1, jac);
		(elm1->get_flx_jac_cont()).set(side, 0, 2, ZERO_MATRIX);

		// why element elm2 which must be a ghost element always exists?????(Hossein asking)
		elm2 = (DualElem*) El_Table->lookup(&neighbor[zp + 4][0]);
		assert(elm2);
		dual_zdirflux(side, hfv, flux_jac, s_jac);
		elm2->dual_zdirflux(side + 2, hfv2, flux_jac2, s_jac2);

		//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
		if (neigh_proc[zp + 4] == myid) {

			zm2 = elm2->which_neighbor(pass_key()) % 4;
			nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zm2 + 4][0]);

			dual_riemannflux(hfv, hfv2, nm2->flux, flux_jac, flux_jac2, s_jac, s_jac2, jac_res,
			    jac_neigh2);

			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				np->flux[ivar] = 0.5 * (np->flux[ivar] + nm2->flux[ivar]);

			flx_jac_cont.set(side, 1, 0, jac, jac_res);
			flx_jac_cont.set(side, 1, 1, jac_neigh1, ZERO_MATRIX);
			flx_jac_cont.set(side, 1, 2, jac_neigh2, ZERO_MATRIX);

			// here the element elm1 must be a ghost element
			(elm2->get_flx_jac_cont()).set(side, 0, 0, jac_neigh2);
			(elm2->get_flx_jac_cont()).set(side, 0, 1, jac_res);
			(elm2->get_flx_jac_cont()).set(side, 0, 2, ZERO_MATRIX);

		} else {

			dual_riemannflux(hfv, hfv2, ghostflux, flux_jac, flux_jac2, s_jac, s_jac2, jac_res,
			    jac_neigh2);
			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
				np->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

			flx_jac_cont.set(side, 1, 0, jac, jac_res);
			flx_jac_cont.set(side, 1, 1, jac_neigh1, ZERO_MATRIX);
			flx_jac_cont.set(side, 1, 2, jac_neigh2);

			// here the element elm1 must be a ghost element
			(elm2->get_flx_jac_cont()).set(side, 0, 0, jac_neigh2);
			(elm2->get_flx_jac_cont()).set(side, 0, 1, jac_res);
			(elm2->get_flx_jac_cont()).set(side, 0, 2, ZERO_MATRIX);

		}

	} else {

		np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
		elm1 = (DualElem*) El_Table->lookup(&neighbor[zp][0]);
		assert(elm1);

		dual_zdirflux(side, hfv, flux_jac, s_jac);
		elm1->dual_zdirflux(side + 2, hfv1, flux_jac1, s_jac1);

		dual_riemannflux(hfv, hfv1, np->flux, flux_jac, flux_jac1, s_jac, s_jac1, jac, jac_neigh1);

		//now we can store jacobians in elements
		flx_jac_cont.set(side, 1, 0, jac);
		flx_jac_cont.set(side, 1, 1, jac_neigh1);
		flx_jac_cont.set(side, 1, 2, ZERO_MATRIX);

		// here the element elm1 must be a ghost element
		(elm1->get_flx_jac_cont()).set(side, 0, 0, jac_neigh1);
		(elm1->get_flx_jac_cont()).set(side, 0, 1, jac);
		(elm1->get_flx_jac_cont()).set(side, 0, 2, ZERO_MATRIX);

		/* CASE I
		 ------------------- -------------------
		 |                   |                   |
		 |                   |                   |
		 |                   |                   |
		 |                   |                   |
		 |                   |                   |
		 |      this       np|         elm1      |
		 |        h          |         hp        |
		 |    kactxy_gz      |    kactxy_gz_n    |
		 |                   |                   |
		 |                   |                   |
		 |                   |                   |
		 ------------------- -------------------
		 */

		/*
		 Case II
		 --------- ----------------------------
		 |         |         |                  |
		 |positive_z_side--->|<----zelmpos      |
		 |         | this  np|                  |
		 |         |   h     |                  |
		 |         |         |                  |
		 |---------|---------|nm1   elm1        |
		 |         |         |      hp          |
		 |         |         |                  |
		 |         | elm2 np2|                  |
		 |         |         |                  |
		 positive_z_side_2-->|                  |
		 --------- ----------------------------
		 */

		if (np->info == S_S_CON) {
			nm1 = NULL;
			np2 = NULL;

			zelmpos = elm1->which_neighbor(pass_key());
			assert(zelmpos > -1);
			nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos % 4 + 4][0]);

			elm2 = (DualElem*) El_Table->lookup(&elm1->neighbor[(zelmpos + 4) % 8][0]);
			assert(elm2);

			elm1->dual_zdirflux(side, hfv1, flux_jac1, s_jac1);
			elm2->dual_zdirflux(side, hfv2, flux_jac2, s_jac2);

			if (*(elm1->get_neigh_proc() + (zelmpos + 4) % 8) == myid) {

				zp2 = elm2->which_neighbor(elm1->pass_key()) % 4;
				np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);

				dual_riemannflux(hfv2, hfv1, np2->flux, flux_jac2, flux_jac1, s_jac2, s_jac1, jac_neigh2,
				    jac_res);

				for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					nm1->flux[ivar] = 0.5 * (np->flux[ivar] + np2->flux[ivar]);

			} else {

				dual_riemannflux(hfv2, hfv1, ghostflux, flux_jac2, flux_jac1, s_jac2, s_jac1, jac_neigh2,
				    jac_res);
				for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					nm1->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);

			}

			(elm1->get_flx_jac_cont()).set(side, 0, 0, jac_neigh1, jac_res);

			if (elm1->which_neighbor(key) < 4) {
				(elm1->get_flx_jac_cont()).set(side, 0, 1, jac, ZERO_MATRIX);
				(elm1->get_flx_jac_cont()).set(side, 0, 2, jac_neigh2, ZERO_MATRIX);
			} else {
				(elm1->get_flx_jac_cont()).set(side, 0, 1, jac_neigh2, ZERO_MATRIX);
				(elm1->get_flx_jac_cont()).set(side, 0, 2, jac, ZERO_MATRIX);
			}
		}

		/*  Case III
		 ------------------- --------- ---------
		 |                   |         |         |
		 |positive_z_side--->|<----zelmpos_2     |
		 |                   |nm2 elm2 |         |
		 |                   |     hp2 |         |
		 |                   |         |         |
		 |       this        |---------|---------
		 |        h          |         |         |
		 |                   |         |         |
		 |                   |nm1 elm1 |         |
		 |                   |     hp1 |         |
		 |                   |<----zelmpos       |
		 ------------------- --------- ---------
		 */

		else if (np->info == S_C_CON) {

			nm1 = NULL;
			nm2 = NULL;

			zelmpos = elm1->which_neighbor(pass_key()) % 4;
			nm1 = (Node*) NodeTable->lookup(&elm1->node_key[zelmpos + 4][0]);

			elm2 = (DualElem*) (El_Table->lookup(&neighbor[zp + 4][0]));
			assert(elm2);

			dual_zdirflux(side, hfv, flux_jac, s_jac);
			elm2->dual_zdirflux(side + 2, hfv2, flux_jac2, s_jac2);

			if (neigh_proc[zp + 4] == myid) {

				zelmpos_2 = elm2->which_neighbor(pass_key()) % 4;
				nm2 = (Node*) NodeTable->lookup(&elm2->node_key[zelmpos_2 + 4][0]);

				dual_riemannflux(hfv, hfv2, nm2->flux, flux_jac, flux_jac2, s_jac, s_jac2, jac_res,
				    jac_neigh2);

				for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
					nm1->flux[ivar] = np->flux[ivar];
					np->flux[ivar] = 0.5 * (nm1->flux[ivar] + nm2->flux[ivar]);

				}
			} else {
				dual_riemannflux(hfv, hfv2, ghostflux, flux_jac, flux_jac2, s_jac, s_jac2, jac_res,
				    jac_neigh2);
				for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
					nm1->flux[ivar] = np->flux[ivar];
					np->flux[ivar] = 0.5 * (nm1->flux[ivar] + ghostflux[ivar]);
				}

			}

			//now we can store jacobians in elements
			flx_jac_cont.set(side, 1, 0, jac, jac_res);
			flx_jac_cont.set(side, 1, 1, jac_neigh1, ZERO_MATRIX); // to have the effect of this element
			flx_jac_cont.set(side, 1, 2, jac_neigh2, ZERO_MATRIX); // to have the effect of this element

			// here the element elm1 must be a ghost element
			(elm1->get_flx_jac_cont()).set(side, 0, 0, jac_neigh1);
			(elm1->get_flx_jac_cont()).set(side, 0, 1, jac);
			(elm1->get_flx_jac_cont()).set(side, 0, 2, ZERO_MATRIX);

			// here the element elm1 must be a ghost element
			(elm2->get_flx_jac_cont()).set(side, 0, 0, jac_neigh2);
			(elm2->get_flx_jac_cont()).set(side, 0, 1, jac_res);
			(elm2->get_flx_jac_cont()).set(side, 0, 2, ZERO_MATRIX);

		}

	}

	if (neigh_proc[zm] != myid) {

		if (neigh_proc[zm] == -1) {

			elem_list->push_back(this);
			// we have to remove them in another function

//			np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
//			nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
//
//			//outflow boundary conditions
//			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
//				nm->flux[ivar] = np->flux[ivar];
//
//			//in above line we do not compute the flux on negative side and just
//			// set it to positive flux which means outlet=inlet
//			//by default fluxes_jac has been initialized to zero, but to be on the safe side
//			for (int i = 0; i < NUM_STATE_VARS; ++i)
//				flx_jac_cont.set(side, 0, i, ZERO_MATRIX);

		} else {
			/* if an interface is on the x-minus or y-minus side, need to
			 calculate those edgestates in this element */
			// x-minus side
			/*
			 interface
			 |
			 |
			 left cells-GHOST CELLS   |
			 (no need to        |
			 calculate         |
			 fluxes )         |
			 v
			 Case I
			 ------------------- -------------------
			 |                   |                   |
			 |                   |<--zm              |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |       elm1        |nm      this       |
			 |        h          |         hp        |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 |                   |                   |
			 ------------------- -------------------
			 Case II
			 --------- -----------------------------
			 |         |         |                   |
			 |         |         |<--zm(z-minus side)|
			 |         |  elm1   |                   |
			 |         |   h     |                   |
			 |         |         |                   |
			 |---------|---------|nm    this         |
			 |         |         |       hp          |
			 |         |         |                   |
			 |         | elm2    |                   |
			 |         |   h2    |                   |
			 |         |         |                   |
			 --------- -----------------------------
			 */

			nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);
			elm1 = (DualElem*) El_Table->lookup(&neighbor[zm][0]);
			assert(elm1);

			dual_zdirflux(side + 2, hfv, flux_jac, s_jac);
			elm1->dual_zdirflux(side, hfv1, flux_jac1, s_jac1);
			dual_riemannflux(hfv1, hfv, nm->flux, flux_jac1, flux_jac, s_jac1, s_jac, jac_neigh1, jac);

			//now we can store jacobians in elements
			flx_jac_cont.set(side, 0, 0, jac);
			flx_jac_cont.set(side, 0, 1, jac_neigh1);
			flx_jac_cont.set(side, 0, 2, ZERO_MATRIX);

			//*** in this section elem1 & elem2 are ghost element, and we do not need to update them here
			// here the element elm1 must be a ghost element
//			elm1->set_fluxjac(side, 1, 0, jac_neigh1);
//			elm1->set_fluxjac(side, 1, 1, jac);
//			elm1->set_fluxjac(side, 1, 2, jac_zero);

			elm2 = (DualElem*) El_Table->lookup(&neighbor[zm + 4][0]);
			assert(elm2);

			dual_zdirflux(side + 2, hfv, flux_jac, s_jac);
			elm2->dual_zdirflux(side, hfv2, flux_jac2, s_jac2);

			//note a rectangular domain ensures that neigh_proc[zm+4]!=-1
			if (neigh_proc[zm + 4] == myid) {
				zp2 = elm2->which_neighbor(pass_key()) % 4;
				np2 = (Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);

				dual_riemannflux(hfv2, hfv, np2->flux, flux_jac2, flux_jac, s_jac2, s_jac, jac_neigh2,
				    jac_res);

				for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					nm->flux[ivar] = 0.5 * (nm->flux[ivar] + np2->flux[ivar]);

				//now we can store jacobians in elements
				flx_jac_cont.set(side, 0, 0, jac, jac_res);
				flx_jac_cont.set(side, 0, 1, jac_neigh1);
				flx_jac_cont.set(side, 0, 2, jac_neigh2);

			} else {
				dual_riemannflux(hfv2, hfv, ghostflux, flux_jac2, flux_jac, s_jac2, s_jac, jac_neigh2,
				    jac_res);

				for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
					nm->flux[ivar] = 0.5 * (nm->flux[ivar] + ghostflux[ivar]);

				//now we can store jacobians in elements
				flx_jac_cont.set(side, 0, 0, jac, jac_res);
				flx_jac_cont.set(side, 0, 1, jac_neigh1);
				flx_jac_cont.set(side, 0, 2, jac_neigh2);

			}

		}
	}
}

void DualElem::boundary_flux(HashTable* El_Table, HashTable* NodeTable, const int myid,
    const int side) {

	int zp = (positive_x_side + side) % 4;
	int zm = (zp + 2) % 4;
	Node* np = (Node*) NodeTable->lookup(&node_key[zp + 4][0]);
	Node* nm = (Node*) NodeTable->lookup(&node_key[zm + 4][0]);

	if (neigh_proc[zp] == -1) {

		// note the this may causes a bug if there is 2 level of refinement
		// for three layers of the cells adjacent to the boundary

		//in the boundary is in the right side of the cell

		//outflow boundary conditions
		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			np->flux[ivar] = nm->flux[ivar];

		//by default fluxes_jac has been initialized to zero, but to be on the safe side
		for (int i = 0; i < NUM_STATE_VARS; ++i)
			flx_jac_cont.set(side, 1, i, flx_jac_cont(side, 0, i));

	}
	if (neigh_proc[zm] == -1) {

		//outflow boundary conditions
		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
			nm->flux[ivar] = np->flux[ivar];

		//in above line we do not compute the flux on negative side and just
		// set it to positive flux which means outlet=inlet
		//by default fluxes_jac has been initialized to zero, but to be on the safe side

		//by default fluxes_jac has been initialized to zero, but to be on the safe side
		for (int i = 0; i < NUM_STATE_VARS; ++i)
			flx_jac_cont.set(side, 0, i, flx_jac_cont(side, 1, i));

	}
}

void DualElem::set_jacobian(int neigh_num, double elemjacob[NUM_STATE_VARS], int state_vars_num,
    const double incr) {

	for (int j = 0; j < NUM_STATE_VARS; ++j)
		jacobianMat(neigh_num)(j, state_vars_num) = (elemjacob[j] / incr);

	return;
}

void DualElem::set_jacobian() {

	Mat3x3 zero_mat;
	for (int i = 0; i < EFF_ELL; i++)
		jacobianMat(i) = zero_mat;

}

Vec_Mat<9>& DualElem::get_jacobian() {
	return jacobianMat;
}

void DualElem::print_jacobian(int iter) {

	cout << "iter:  " << iter << '\n';
	cout << "key1:  " << key[0] << "  key2:  " << key[1] << '\n';
	cout << "X:  " << coord[0] << "  Y:  " << coord[1] << '\n';
	cout << "Jacobian: " << '\n';
//cout << "self"<<"
	for (int i = 0; i < EFF_ELL; i++) {
		cout << "Matrix=  " << i << "," << '\n';

		for (int j = 0; j < NUM_STATE_VARS; j++) {
			for (int k = 0; k < NUM_STATE_VARS; k++) {
				cout << scientific << setw(10) << setprecision(8) << jacobianMat(i, j, k) << "  ";
				if (dabs(jacobianMat(i, j, k)) > 10.)
					cout << "Jedi begir mano" << endl;
			}
			cout << '\n';
		}
	}
}

void DualElem::update_state(SolRec* solrec, HashTable* El_Table, int iter) {

//	int aa = 0, bb = 1;
//	unsigned keyy[2] = { 3781669179, 330382100 };
//	if (key[0] == keyy[0] && key[1] == keyy[1])
//		bb = aa;

	Solution* prev_sol = solrec->lookup(key, iter - 1);

	for (int i = 0; i < NUM_STATE_VARS; i++)
		prev_state_vars[i] = *(prev_sol->get_solution() + i);

	kactxy[0] = prev_sol->get_kact();
	kactxy[1] = kactxy[0];

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		prev_adjoint[i] = adjoint[i];

	if (prev_sol != &(Solution::solution_zero))
		delete prev_sol;

}

void DualElem::set_jacobianMat_zero(int jacmatind) {

	jacobianMat(jacmatind) = ZERO_MATRIX;

}

void DualElem::dual_check_refine_unrefine(SolRec* solrec, HashTable* El_Table, int iter,
    ElemPtrList<DualElem>* refinelist, ElemPtrList<DualElem>* unrefinelist) {

//	int aa = 0, bb = 1;
//	unsigned keyy[2] = { 3410598297, 2576980374 };
//	if (key[0] == keyy[0] && key[1] == keyy[1])
//		bb = aa;

	Solution* prev_sol = solrec->lookup(key, iter - 1);

	if (!prev_sol) {
// first we check to see whether the element has been refined, so we have to read from its father
		prev_sol = solrec->lookup(getfather(), iter - 1);

		if (prev_sol)
			unrefinelist->add(this);

		else
			// the only remaining case is that it has been unrefined so, we have to read from its sons
			refinelist->add(this);
	}
}

//==========================================================================

ErrorElem::ErrorElem(Element* element) {

	myprocess = element->myprocess;

	generation = element->generation;

	opposite_brother_flag = element->opposite_brother_flag;

	material = element->material;

	lb_weight = element->lb_weight;

	refined = element->refined;

	adapted = element->adapted;

	which_son = element->which_son;

	new_old = element->new_old;

	shortspeed = element->shortspeed;

	positive_x_side = element->positive_x_side;

	elevation = element->elevation;

	stoppedflags = element->stoppedflags;

	effect_bedfrict = element->effect_bedfrict;

	effect_tanbedfrict = element->effect_tanbedfrict;

	counted = element->counted;

	ithelem = element->ithelem;

	iwetnode = element->iwetnode;

	Awet = element->Awet;

	Swet = element->Swet;

	for (int i = 0; i < NUM_STATE_VARS; ++i) {

		state_vars[i] = element->state_vars[i];

		prev_state_vars[i] = element->prev_state_vars[i];

		gravity[i] = element->gravity[i];

		Influx[i] = element->Influx[i];

		d_state_vars[i] = element->d_state_vars[i];
		d_state_vars[i + NUM_STATE_VARS] = element->d_state_vars[i + NUM_STATE_VARS];

	}

	for (int i = 0; i < DIMENSION; ++i) {

		dx[i] = element->dx[i];

		eigenvxymax[i] = element->eigenvxymax[i];

		kactxy[i] = element->kactxy[i];

		zeta[i] = element->zeta[i];

		curvature[i] = element->curvature[i];

		d_gravity[i] = element->d_gravity[i];

		effect_kactxy[i] = element->effect_kactxy[i];

		drypoint[i] = element->drypoint[i];

		coord[i] = element->coord[i];

		elm_loc[i] = element->elm_loc[i];

		lb_key[i] = element->lb_key[i];

		key[i] = element->key[i];

		father[i] = element->father[i];

		el_error[i] = element->el_error[i];
	}

	for (int i = 0; i < 8; ++i) {

		neigh_gen[i] = element->neigh_gen[i];

		neigh_proc[i] = element->neigh_proc[i];

		for (int j = 0; j < KEYLENGTH; ++j) {

			node_key[i][j] = element->node_key[i][j];

			neighbor[i][j] = element->neighbor[i][j];
		}
	}

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 2; ++j) {
			son[i][j] = element->son[i][j];
			brothers[i][j] = element->brothers[i][j];
		}

	for (int i = 0; i < NUM_STATE_VARS; ++i) {

		adjoint[i] = 0.;

		bilin_adj[i] = 0.;

		residual[i] = 0.;

		bilin_state[i] = 0.;
	}

	correction = 0.;

	father_address = NULL;

}

//used for refinement
ErrorElem::ErrorElem(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], int n_pro[],
    int gen, int elm_loc_in[], int gen_neigh[], int mat, ErrorElem *fthTemp, double *coord_in,
    HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather,
    double Awetfather, double *drypoint_in, int SETLINK) {

	counted = 0; //for debugging only

	adapted = NEWSON;

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		prev_state_vars[i] = 0.;
		Influx[i] = 0.;
		adjoint[i] = 0.;
		bilin_adj[i] = 0.;
		Influx[i] = 0.;
		residual[i] = 0.;
	}
	for (int i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
		d_state_vars[i] = 0.;

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	for (int i = 0; i < 4; i++)
		son[i][0] = son[i][1] = brothers[i][0] = brothers[i][1] = NULL;
	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	myprocess = myid;
	generation = gen; //--first generation
	opposite_brother_flag = 1;
	material = mat;
	for (int i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (int i = 0; i < KEYLENGTH; i++) {
		father[i] = fthTemp->key[i];
		key[i] = nodekeys[8][i]; //--using buble key to represent the element
	}

	elm_loc[0] = elm_loc_in[0];
	elm_loc[1] = elm_loc_in[1];

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < KEYLENGTH; j++)
			node_key[i][j] = nodekeys[i][j];

	for (int i = 0; i < 4; i++) {
		neigh_proc[i] = n_pro[i];
		neigh_proc[i + 4] = -2; //-- -2 means regular element
		if (neigh_proc[i] != -1)
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = neigh[i][j];
		else
			for (int j = 0; j < KEYLENGTH; j++)
				neighbor[i][j] = neighbor[i + 4][j] = NULL;

		neigh_gen[i] = neigh_gen[i + 4] = gen_neigh[i];
	}

	refined = 0;

	new_old = NEW;
//geoflow stuff
	dx[0] = .5 * fthTemp->dx[0];  //assume constant refinement
	dx[1] = .5 * fthTemp->dx[1];

	iwetnode = iwetnodefather;
	drypoint[0] = drypoint_in[0];
	drypoint[1] = drypoint_in[1];

	double myfractionoffather;
	if ((Awetfather == 0.0) || (Awetfather == 1.0)) {
		Awet = Awetfather;
		myfractionoffather = 1.0;
	} else {
		Awet = convect_dryline(dx, 0.0); //dx is a dummy stand in for convection speed... value doesn't matter because it's being multiplied by a timestep of zero
		myfractionoffather = Awet / Awetfather;
	}
	Swet = 1.0;

	double dxx = coord_in[0] - fthTemp->coord[0];
	double dyy = coord_in[1] - fthTemp->coord[1];

	if (state_vars[0] < 0.)
		state_vars[0] = 0.;

	find_positive_x_side(NodeTable);
	calc_topo_data(matprops_ptr);
	calc_gravity_vector(matprops_ptr);
	calc_d_gravity(El_Table);

	coord[0] = coord_in[0];
	coord[1] = coord_in[1];

	calc_which_son();

	stoppedflags = fthTemp->stoppedflags;

	correction = 0.0;
	kactxy[0] = fthTemp->kactxy[0];
	kactxy[1] = fthTemp->kactxy[1];
	effect_kactxy[0] = fthTemp->effect_kactxy[0];
	effect_kactxy[1] = fthTemp->effect_kactxy[1];
	for (int i = 0; i < NUM_STATE_VARS; i++) {
		state_vars[i] = fthTemp->state_vars[i] * myfractionoffather;
		prev_state_vars[i] = fthTemp->prev_state_vars[i] * myfractionoffather;
		adjoint[i] = 0.25 * fthTemp->adjoint[i];
		bilin_adj[i] = 0.25 * fthTemp->adjoint[i];
		Influx[i] = 0.;
		residual[i] = 0.;
	}

	int aa = 1, bb = 0;
	if (fthTemp->key[0] == 3897323835 && fthTemp->key[1] == 330382099)
		aa = bb;

	if (SETLINK) {

		father_address = fthTemp->father_address;

		// we just do that for the son 0
		int count = 0;
		if (which_son == 0) { // first deleting the old links
			vector<ErrorElem*>& mysons = father_address->get_son_addresses();
			vector<ErrorElem*>::iterator it = mysons.begin();

			while (it != mysons.end())
				if (*it == fthTemp) {
					mysons.erase(it);
					count++;
				} else
					++it;
		}

		if (which_son == 0)
			assert(count == 1);

		// now adding the new link
		(father_address->get_son_addresses()).push_back(this);
	}
}

/*********************************
 making a father element from its sons
 *****************************************/
ErrorElem::ErrorElem(ErrorElem* sons[], HashTable* NodeTable, HashTable* El_Table,
    MatProps* matprops_ptr, int SETLINK) {
	counted = 0; //for debugging only

	adapted = NEWFATHER;

	for (int ikey = 0; ikey < KEYLENGTH; ikey++)
		father[ikey] = brothers[0][ikey] = brothers[1][ikey] = brothers[2][ikey] = brothers[3][ikey] =
		    son[0][ikey] = son[1][ikey] = son[2][ikey] = son[3][ikey] = 0;

	int i, j, ikey, ison, isonneigh, ineigh;

	for (ikey = 0; ikey < KEYLENGTH; ikey++)
		key[ikey] = *(sons[2]->getNode() + ikey);

	for (ison = 0; ison < 4; ison++) {
		sons[ison]->put_adapted_flag(OLDSON);
		for (ikey = 0; ikey < KEYLENGTH; ikey++) {
			son[ison][ikey] = *(sons[ison]->pass_key() + ikey);
			sons[ison]->put_father(key);
		}
	}

	lb_key[0] = lb_key[1] = NULL;
	lb_weight = 1.0;
	new_old = NEW;
	unsigned* son_nodes[4];
	opposite_brother_flag = 0;
	stoppedflags = 2;
	for (i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.0;

	for (ison = 0; ison < 4; ison++) {
		son_nodes[ison] = sons[ison]->getNode();
		if (sons[ison]->stoppedflags < stoppedflags)
			stoppedflags = sons[ison]->stoppedflags;
	}

	for (ikey = 0; ikey < KEYLENGTH; ikey++) {
		father[ikey] = NULL;
		node_key[0][ikey] = son_nodes[0][ikey];
		node_key[1][ikey] = son_nodes[1][KEYLENGTH + ikey];
		node_key[2][ikey] = son_nodes[2][2 * KEYLENGTH + ikey];
		node_key[3][ikey] = son_nodes[3][3 * KEYLENGTH + ikey];
		node_key[4][ikey] = son_nodes[0][KEYLENGTH + ikey];
		node_key[5][ikey] = son_nodes[1][2 * KEYLENGTH + ikey];
		node_key[6][ikey] = son_nodes[2][3 * KEYLENGTH + ikey];
		node_key[7][ikey] = son_nodes[3][ikey];

		elm_loc[ikey] = (*(sons[0]->get_elm_loc() + ikey)) / 2;
	}
	myprocess = sons[0]->get_myprocess();
	generation = sons[0]->get_gen() - 1;
	material = sons[0]->get_material();

	calc_which_son();

	refined = 1; // not an active element yet!!!

	// neighbor information
	for (ison = 0; ison < 4; ison++) {
		isonneigh = ison;
		ineigh = isonneigh;
		neigh_gen[ineigh] = *(sons[ison]->get_neigh_gen() + isonneigh);
		for (ikey = 0; ikey < KEYLENGTH; ikey++)
			neighbor[ineigh][ikey] = *(sons[ison]->get_neighbors() + isonneigh * KEYLENGTH + ikey);
		neigh_proc[ineigh] = *(sons[ison]->get_neigh_proc() + isonneigh);

		isonneigh = (ison + 3) % 4;
		ineigh = isonneigh + 4;
		neigh_gen[ineigh] = *(sons[ison]->get_neigh_gen() + isonneigh);
		for (ikey = 0; ikey < KEYLENGTH; ikey++)
			neighbor[ineigh][ikey] = *(sons[ison]->get_neighbors() + isonneigh * KEYLENGTH + ikey);
		if ((*(sons[ison]->get_neigh_gen() + isonneigh) == generation)
		    || (*(sons[ison]->get_neigh_proc() + isonneigh) == -1))
			neigh_proc[ineigh] = -2;
		else
			neigh_proc[ineigh] = *(sons[ison]->get_neigh_proc() + isonneigh);
	}

	/* brother information -- requires that at least one of this
	 element's neighboring brothers is on this process in
	 order to get information on the brother that is not a neighbor */
	ErrorElem* EmTemp;
	switch (which_son) {
		case 0:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[0][i] = key[i];
			if (neigh_proc[1] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = 0;
			} else if (neigh_gen[1] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = neighbor[1][i];
			} else if (neigh_gen[1] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[1]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = bro_key[i];
			} else if (neigh_gen[1] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[2] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = 0;
			} else if (neigh_gen[2] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = neighbor[2][i];
			} else if (neigh_gen[2] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[2]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = bro_key[i];
			} else if (neigh_gen[2] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 1:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[1][i] = key[i];
			if (neigh_proc[3] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = 0;
			} else if (neigh_gen[3] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = neighbor[3][i];
			} else if (neigh_gen[3] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[3]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = bro_key[i];
			} else if (neigh_gen[3] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[2] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = 0;
			} else if (neigh_gen[2] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = neighbor[2][i];
			} else if (neigh_gen[2] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[2]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = bro_key[i];
			} else if (neigh_gen[2] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 2:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[2][i] = key[i];
			if (neigh_proc[0] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = 0;
			} else if (neigh_gen[0] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = neighbor[0][i];
			} else if (neigh_gen[0] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[0]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = bro_key[i];
			} else if (neigh_gen[0] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[1][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[3] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = 0;
			} else if (neigh_gen[3] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = neighbor[3][i];
			} else if (neigh_gen[3] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[3]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = bro_key[i];
			} else if (neigh_gen[3] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[3][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
		case 3:
			for (i = 0; i < KEYLENGTH; i++)
				brothers[3][i] = key[i];
			if (neigh_proc[0] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = 0;
			} else if (neigh_gen[0] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = neighbor[0][i];
			} else if (neigh_gen[0] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[0]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = bro_key[i];
			} else if (neigh_gen[0] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[0][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			if (neigh_proc[1] == INIT) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = 0;
			} else if (neigh_gen[1] == generation) {
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = neighbor[1][i];
			} else if (neigh_gen[1] == generation + 1) {
				EmTemp = (ErrorElem*) El_Table->lookup(neighbor[1]);
				assert(EmTemp);
				unsigned* bro_key = EmTemp->getfather();
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = bro_key[i];
			} else if (neigh_gen[1] == generation + 2) {
				//this should not happen except in unrefinement in dual rin
				for (i = 0; i < KEYLENGTH; i++)
					brothers[2][i] = MUST_BE_CORRECTED;
			} else
				assert(0);
			break;
	}

	find_positive_x_side(NodeTable);  //also inserts the coordinates
	calculate_dx(NodeTable);
	find_opposite_brother(El_Table);
	calc_topo_data(matprops_ptr);
	calc_gravity_vector(matprops_ptr);
	calc_d_gravity(El_Table);

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		state_vars[i] = 0.;
		prev_state_vars[i] = 0.;
		Influx[i] = 0.;
		adjoint[i] = 0.;
		bilin_adj[i] = 0.;
		Influx[i] = 0.;
		residual[i] = 0.;
		d_state_vars[i] = d_state_vars[NUM_STATE_VARS + i] = 0.;
	}

	Awet = 0.0;
	for (int ison = 0; ison < 4; ison++)
		Awet += sons[ison]->get_Awet();
	Awet *= 0.25;

	//uninitialized flag values... will fix shortly
	drypoint[0] = drypoint[1] = 0.0;
	iwetnode = 8;
	Swet = 1.0;

	//calculate the shortspeed
	shortspeed = 0.0;
	for (j = 0; j < 4; j++)
		shortspeed += *(sons[j]->get_state_vars()) * sons[j]->get_shortspeed();
	if (state_vars[0] > 0.)
		shortspeed /= (4.0 * state_vars[0]);

	kactxy[0] = effect_kactxy[0] = 0.;
	kactxy[1] = effect_kactxy[1] = 0.;
	correction = 0.0;

	for (i = 0; i < EQUATIONS; i++)
		el_error[i] = 0.;

	if (SETLINK) {
		// correcting the link

		for (i = 1; i < 4; i++)
			assert(sons[i]->get_father_address() == sons[0]->get_father_address());

		father_address = sons[0]->get_father_address();

		// first deleting the old links
		vector<ErrorElem*>& mysons = father_address->get_son_addresses();

		int count = 0;
		for (int j = 0; j < 4; ++j) {
			vector<ErrorElem*>::iterator it = mysons.begin();
			while (it != mysons.end())
				if (*it == sons[j]) {
					mysons.erase(it);
					count++;
				} else
					++it;
		}

		assert(count == 4);

		// now adding the new link
		(father_address->get_son_addresses()).push_back(this);
	}
}

void ErrorElem::get_slopes_prev(HashTable* El_Table, HashTable* NodeTable, double gamma) {
	int j = 0, bc = 0;
	/* check to see if this is a boundary */
	while (j < 4 && bc == 0) {
		if (neigh_proc[j] == INIT)
			bc = 1;
		j++;
	}
	if (bc == 1) {
		for (j = 0; j < NUM_STATE_VARS * DIMENSION; j++)
			d_state_vars[j] = 0;

		return;
	}

	int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	xp = positive_x_side;
	switch (positive_x_side) {
		case 0:
			xm = 2;
			yp = 1;
			ym = 3;
			break;
		case 1:
			xm = 3;
			yp = 2;
			ym = 0;
			break;
		case 2:
			xm = 0;
			yp = 3;
			ym = 1;
			break;
		case 3:
			xm = 1;
			yp = 0;
			ym = 2;
			break;
	}

	/* x direction */
	ErrorElem *ep = (ErrorElem*) (El_Table->lookup(&neighbor[xp][0]));
	ErrorElem *em = (ErrorElem*) (El_Table->lookup(&neighbor[xm][0]));
	ErrorElem *ep2 = NULL;
	ErrorElem *em2 = NULL;

//check if element has 2 neighbors on either side
	Node* ndtemp = (Node*) NodeTable->lookup(&node_key[xp + 4][0]);

	if (ndtemp->info == S_C_CON) {
		ep2 = (ErrorElem*) (El_Table->lookup(&neighbor[xp + 4][0]));
		assert(neigh_proc[xp + 4] >= 0 && ep2);
	}

	ndtemp = (Node*) NodeTable->lookup(&node_key[xm + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (ErrorElem*) (El_Table->lookup(&neighbor[xm + 4][0]));
		assert(neigh_proc[xm + 4] >= 0 && em2);
	}

	double dp, dm, dc, dxp, dxm;
	dxp = ep->coord[0] - coord[0];
	dxm = coord[0] - em->coord[0];

	double inv_dxp = 1 / dxp, inv_dxm = 1 / dxm, inv_sdxpdxm = 1 / (dxm + dxp);

	double min_slopes;

	for (j = 0; j < NUM_STATE_VARS; j++) {

		dp = (ep->bilin_prev_state[j] - bilin_prev_state[j]) * inv_dxp;

		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->bilin_prev_state[j] - bilin_prev_state[j]) * inv_dxp);

		dm = (bilin_prev_state[j] - em->bilin_prev_state[j]) * inv_dxm;

		if (em2 != NULL)
			dm = .5 * (dm + (bilin_prev_state[j] - em2->bilin_prev_state[j]) * inv_dxm);

		dc = (dp * dxm + dm * dxp) * inv_sdxpdxm;  // weighted average

		min_slopes = c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));

		//do slope limiting
		d_state_vars[j] = .5 * (c_sgn(dp) + c_sgn(dm)) * min_slopes;
	}

	/* y direction */
	ep = (ErrorElem*) (El_Table->lookup(&neighbor[yp][0]));
	em = (ErrorElem*) (El_Table->lookup(&neighbor[ym][0]));
	ep2 = NULL;
	em2 = NULL;
//check if element has 2 neighbors on either side
	ndtemp = (Node*) NodeTable->lookup(&node_key[yp + 4][0]);
	if (ndtemp->info == S_C_CON) {
		ep2 = (ErrorElem*) (El_Table->lookup(&neighbor[yp + 4][0]));
		assert(neigh_proc[yp + 4] >= 0 && ep2);
	}
	ndtemp = (Node*) NodeTable->lookup(&node_key[ym + 4][0]);
	if (ndtemp->info == S_C_CON) {
		em2 = (ErrorElem*) (El_Table->lookup(&neighbor[ym + 4][0]));
		if (!(neigh_proc[ym + 4] >= 0 && em2)) {
			printf("ym=%d neigh_proc[ym+4]=%d em2=%d\n", ym, neigh_proc[ym + 4], em2);
		}
		assert(neigh_proc[ym + 4] >= 0 && em2);
	}

	dxp = ep->coord[1] - coord[1];
	dxm = coord[1] - em->coord[1];

	inv_dxp = 1 / dxp;
	inv_dxm = 1 / dxm;
	inv_sdxpdxm = 1 / (dxm + dxp);

	for (j = 0; j < NUM_STATE_VARS; j++) {
		dp = (ep->bilin_prev_state[j] - bilin_prev_state[j]) * inv_dxp;
		if (ep2 != NULL)
			dp = .5 * (dp + (ep2->bilin_prev_state[j] - bilin_prev_state[j]) * inv_dxp);
		dm = (bilin_prev_state[j] - em->bilin_prev_state[j]) * inv_dxm;
		if (em2 != NULL)
			dm = .5 * (dm + (prev_state_vars[j] - em2->prev_state_vars[j]) * inv_dxm);

		dc = (dp * dxm + dm * dxp) * inv_sdxpdxm;  // weighted average

		min_slopes = c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
		//do slope limiting
		d_state_vars[j + NUM_STATE_VARS] = .5 * (c_sgn(dp) + c_sgn(dm)) * min_slopes;

	}
}

void ErrorElem::error_update_state(SolRec* solrec, int iter) {

	Solution* prev_sol = solrec->lookup(this->getfather(), iter - 2);

	for (int i = 0; i < NUM_STATE_VARS; i++)
		prev_state_vars[i] = *(prev_sol->get_solution() + i);

	kactxy[0] = prev_sol->get_kact();
	kactxy[1] = kactxy[0];

}

void ErrorElem::error_check_refine_unrefine(SolRec* solrec, HashTable* El_Table, int iter,
    ElemPtrList<ErrorElem>* refinelist, ElemPtrList<ErrorElem>* unrefinelist) {

	Solution* prev_sol = solrec->lookup(getfather(), iter - 2);

	if (!prev_sol) {

		prev_sol = solrec->lookup(key, iter - 2);

		if (prev_sol)
			refinelist->add(this);

		else
// the only remaining case is that it has been unrefined so, we have to read from its sons
			unrefinelist->add(this);
	}
}

//x direction flux in current cell
void ErrorElem::xdirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]) {
	int i, j;
	double a, Vel;

	if (bilin_prev_state[0] < GEOFLOW_TINY) {

		for (i = 0; i < 3; i++)
			for (j = 0; j < NUM_STATE_VARS; j++)
				hfv[i][j] = 0.;

	} else {
		//state variables
		for (i = 0; i < NUM_STATE_VARS; i++)
			hfv[0][i] = bilin_prev_state[i];

		// Solid-phase velocity in x-dir
		Vel = hfv[0][1] / hfv[0][0];

		// sound-speed : a^2 = k_ap*h*g(3)
		a = sqrt(kactxy[0] * hfv[0][0] * gravity[2]);

		//fluxes
		hfv[1][0] = hfv[0][1];
		hfv[1][1] = hfv[0][1] * Vel + 0.5 * a * a * hfv[0][0];
		hfv[1][2] = hfv[0][2] * Vel;

		//wave speeds
		hfv[2][0] = Vel - a;
		hfv[2][1] = Vel;
		hfv[2][2] = Vel + a;
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			hrfv[i][j] = hfv[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hrfv[i][j]))
				cout << "flux is NAN" << endl;

}

void ErrorElem::ydirflux(MatProps* matprops_ptr, double dz, double wetnessfactor,
    double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]) {
	int i, j;
	double Vel, a;

//the "update flux" values (hfv) are the fluxes used to update the solution,
// they may or may not be "reset" from their standard values based on whether
// or not the stopping criteria is triggering a change intended to cause the flow to stop.
	if (bilin_prev_state[0] < GEOFLOW_TINY) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < NUM_STATE_VARS; j++)
				hfv[i][j] = 0.0; //state variables

	} else {
		//state variables
		for (i = 0; i < NUM_STATE_VARS; i++)
			hfv[0][i] = bilin_prev_state[i];

		// a = speed of sound through the medium
		a = sqrt(kactxy[1] * hfv[0][0] * gravity[2]);

		// Solid-phase velocity in y-dir
		Vel = hfv[0][2] / hfv[0][0];

		//fluxes
		hfv[1][0] = hfv[0][2];
		hfv[1][1] = hfv[0][1] * Vel;
		hfv[1][2] = hfv[0][2] * Vel + 0.5 * a * a * hfv[0][0];

		//wave speeds
		hfv[2][0] = Vel - a;
		hfv[2][1] = Vel;
		hfv[2][2] = Vel + a;
	}

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			hrfv[i][j] = hfv[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < NUM_STATE_VARS; j++)
			if (isnan(hrfv[i][j]))
				cout << "flux is NAN" << endl;
}

