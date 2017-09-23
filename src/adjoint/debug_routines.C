/*
 * test_routines.C
 *
 *  Created on: Mar 3, 2016
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"
#include "../header/exvar.h"
#include <algorithm>
#include <map>

void copy_jacobian(HashTable* El_Table, map<int, Vec_Mat<9> >& jac_map) {

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					jac_map[Curr_El->get_ithelem()] = Curr_El->get_jacobian();
				}
			}
		}

}

void compare_jacobians(map<int, Vec_Mat<9> >& jac_code, map<int, Vec_Mat<9> >& jac_diff) {

	for (map<int, Vec_Mat<9> >::iterator it = jac_code.begin(); it != jac_code.end(); ++it) {
		Vec_Mat<9>& jacdiff = jac_diff[it->first];
		Vec_Mat<9>& jaccode = jac_code[it->first];

		for (int i = 0; i < EFF_ELL; ++i)
			for (int j = 0; j < NUM_STATE_VARS; ++j)
				for (int k = 0; k < NUM_STATE_VARS; ++k)
					if (fabs(jacdiff(i, j, k) - jaccode(i, j, k)) > 1e-11 && jacdiff(i, j, k) != 0.) {
						double denum;
						if (jacdiff(i, j, k) != 0.)
							denum = dabs(jacdiff(i, j, k));
						cout << "in element  " << it->first << " in indices: " << i << " , " << j << " , " << k
						    << "  there is a difference of  "
						    << 100 * dabs(jacdiff(i, j, k) - jaccode(i, j, k)) / denum << "  jacobian diff= "
						    << jacdiff(i, j, k) << "  jacobian code=" << jaccode(i, j, k) << endl;
					}

	}

}

void check_kackxy(HashTable* El_Tab) {

	HashEntryPtr currentPtr;
	Element *Curr_El;

	int aa = 1, bb = 0;

	HashEntryPtr *buck = El_Tab->getbucketptr();
	for (int i = 0; i < El_Tab->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					if (*(Curr_El->get_kactxy()) < 0 || *(Curr_El->get_kactxy() + 1) < 0)
						aa = bb;

				currentPtr = currentPtr->next;
			}
		}

}

void clean_jacobian(HashTable* El_Table) {
	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					Vec_Mat<9>& jacobian = Curr_El->get_jacobian();
					for (int i = 0; i < EFF_ELL; ++i)
						jacobian(i) = ZERO_MATRIX;
				}
			}
		}
}

void check_state_vars_with_record(HashTable* El_Table, SolRec* solrec, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					if (Curr_El->check_state(solrec, El_Table, iter))
						cout << "the idea did not work \n";

				currentPtr = currentPtr->next;
			}
		}
}

void print_Elem_Table(HashTable* El_Table, HashTable* NodeTable, int iter, int place) {

	ofstream myfile;
	char filename[50];
	sprintf(filename, "El_Table_%d_%08d", place, iter);

	myfile.open(filename, ios::app);

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					int print = 1;
					for (int j = 0; j < 4; j++)
						if (*(Curr_El->get_neigh_proc() + j) == INIT) {
							print = 0;
							break;
						}
					if (print) {

						int xp = Curr_El->get_positive_x_side(); //finding the direction of element
						int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

						Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

						Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

						Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

						Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

						double flux[4][NUM_STATE_VARS];

						for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
							flux[0][ivar] = *(nxp->get_flux() + ivar);
							flux[1][ivar] = *(nyp->get_flux() + ivar);
							flux[2][ivar] = *(nxm->get_flux() + ivar);
							flux[3][ivar] = *(nym->get_flux() + ivar);
						}

						double flux_diff[DIMENSION][NUM_STATE_VARS];

						for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
							flux_diff[0][ivar] = flux[0][ivar] - flux[2][ivar];
							flux_diff[1][ivar] = flux[1][ivar] - flux[3][ivar];
						}

						myfile << "key: ";
						myfile << *(Curr_El->pass_key()) << " " << *(Curr_El->pass_key() + 1) << " ";

						myfile << "neighbors: ";
						for (int i = 0; i < 8 * 2; ++i)
							myfile << *(Curr_El->get_neighbors() + i) << " ";

						myfile << "neighbors_gen: ";
						for (int i = 0; i < 8; ++i)
							myfile << *(Curr_El->get_neigh_gen() + i) << " ";

						myfile << "neighbors_proc: ";
						for (int i = 0; i < 8; ++i)
							myfile << *(Curr_El->get_neigh_proc() + i) << " ";

						myfile << "state_vars: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << *(Curr_El->get_state_vars() + i) << " ";

						myfile << "prev_state_vars: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << *(Curr_El->get_prev_state_vars() + i) << " ";

						myfile << "d_state_vars: ";
						for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
							myfile << *(Curr_El->get_d_state_vars() + i) << " ";

//					myfile << "flux_diff_x: ";
//					for (int i = 0; i < NUM_STATE_VARS; ++i)
//						myfile << flux_diff[0][i] << " ";
//
//					myfile << "flux_diffy: ";
//					for (int i = 0; i < NUM_STATE_VARS; ++i)
//						myfile << flux_diff[1][i] << " ";
//
						myfile << "flux_xp: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[0][i] << " ";

						myfile << "flux_yp: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[1][i] << " ";

						myfile << "flux_xm: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[2][i] << " ";

						myfile << "flux_ym: ";
						for (int i = 0; i < NUM_STATE_VARS; ++i)
							myfile << flux[3][i] << " ";

						myfile << "elem_id: " << Curr_El->get_ithelem() << " ";

//				myfile <<"gravity: ";
//				for (int i = 0; i < NUM_STATE_VARS; ++i)
//					myfile << *(Curr_El->get_gravity() + i) << " ";
//
//				myfile <<"d_gravity: ";
//				for (int i = 0; i < DIMENSION; ++i)
//					myfile << *(Curr_El->get_d_gravity() + i) << " ";

						myfile << endl;
					}
				}

				currentPtr = currentPtr->next;
			}
		}

	myfile.close();

}

void set_ithm(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	int count = 0;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					Curr_El->put_ithelem(count);
					count++;
				}
			}
		}
}

void plot_ithm(HashTable* El_Table) {
	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {
					Curr_El->put_ithelem(Curr_El->get_ithelem());
				}
			}
		}
}

void refinement_report(HashTable* El_Table, int myid) {

	int newbuffer = 0, buffer = 0, newson = 0, newfather = 0, norecadapt = 0, ghnewbuffer = 0,
	    ghbuffer = 0, ghnewson = 0, ghnewfather = 0, ghnotrecadapt = 0, tobedeleted = 0,
	    oldfather = 0, oldson = 0;

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				switch (Curr_El->get_adapted_flag()) {
					case NEWBUFFER:
						newbuffer++;
						break;
					case BUFFER:
						buffer++;
						break;
					case NEWSON:
						newson++;
						break;
					case NEWFATHER:
						newfather++;
						break;
					case NOTRECADAPTED:
						norecadapt++;
						break;
					case TOBEDELETED:
						tobedeleted++;
						break;
					case -NOTRECADAPTED:
						ghnotrecadapt++;
						break;
					case -NEWFATHER:
						ghnewfather++;
						break;
					case -NEWSON:
						ghnewson++;
						break;
					case -BUFFER:
						ghbuffer++;
						break;
					case -NEWBUFFER:
						ghnewbuffer++;
						break;
					case OLDFATHER:
						oldfather++;
						break;
					case OLDSON:
						oldson++;
						break;
					default:
						cout << "this case is irregular \n";
				}
				currentPtr = currentPtr->next;
			}
		}

	cout << "Proc ID:" << myid << "\n new buffer:       " << newbuffer << "\n buffer:           "
	    << buffer << "\n newson:           " << newson << "\n newfather:        " << newfather
	    << "\n norecadapt:       " << norecadapt << "\n tobedeleted:      " << tobedeleted
	    << "\n ghost new buffer: " << ghnewbuffer << "\n ghost buffer:     " << ghbuffer
	    << "\n ghostnewson:      " << ghnewson << "\n ghost newfather:  " << ghnewfather
	    << "\n ghost norecadapt: " << ghnotrecadapt << "\n oldfather:        " << oldfather
	    << "\n oldson:           " << oldson << "\n";

}

void refine_flag_report(HashTable* El_Table, int myid) {

	int one = 0, zero = 0, ghost = 0;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				switch (Curr_El->get_refined_flag()) {
					case 0:
						zero++;
						break;
					case 1:
						one++;
						break;
					case GHOST:
						ghost++;
						break;
					default:
						cout << "this case is irregular" << Curr_El->get_refined_flag() << "\n";
				}
				currentPtr = currentPtr->next;
			}
		}

	cout << "Proc ID:" << myid << "\n refine=1:    " << one << "\n refine=0:    " << zero
	    << "\n refine=GHOST:" << ghost << "\n";
}

void print_jacobian(HashTable* El_Table, int iter) {

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->print_jacobian(iter);

				currentPtr = currentPtr->next;
			}
		}
	}
}

int num_nonzero_elem(HashTable *El_Table, int type) {
	int num = 0;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() == type)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

int num_nonzero_elem(HashTable *El_Table) {
	int num = 0;			//myid
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

int table_members(HashTable *NodeTable) {

	int num = 0;
	HashEntryPtr currentPtr;
	HashEntryPtr *buck = NodeTable->getbucketptr();

	for (int i = 0; i < NodeTable->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				num++;
				currentPtr = currentPtr->next;
			}
		}
return num;
}

void get_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key, MatProps* matprops_ptr,
    int myid, double flux[4][NUM_STATE_VARS]) {

	Element* Curr_El = (Element*) (El_Table->lookup(key));

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

	Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

	Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

	Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
		flux[0][ivar] = nxp->flux[ivar];
		flux[1][ivar] = nyp->flux[ivar];
		flux[2][ivar] = nxm->flux[ivar];
		flux[3][ivar] = nym->flux[ivar];
	}
}

void flux_debug(Element* Curr_El, double* fluxxpold, double* fluxxmold, double* fluxypold,
    double* fluxymold, double* fluxxp, double* fluxxm, double* fluxyp, double* fluxym,
    int effelement, int j, int iter, double dt) {

	double diff1[NUM_STATE_VARS], diff2[NUM_STATE_VARS], diff3[NUM_STATE_VARS], diff4[NUM_STATE_VARS];
	double fluxxold[NUM_STATE_VARS], fluxxnew[NUM_STATE_VARS], fluxyold[NUM_STATE_VARS],
	    fluxynew[NUM_STATE_VARS];
	double abs_fluxx_diff[NUM_STATE_VARS], abs_fluxy_diff[NUM_STATE_VARS];

	double *state_vars = Curr_El->get_state_vars();
	double *prev_state_vars = Curr_El->get_prev_state_vars();
	double *d_state_vars = Curr_El->get_d_state_vars();
	double *dx = Curr_El->get_dx();
	double dtdx = dt / dx[0];
	double dtdy = dt / dx[1];

	FILE *fp;
	fp = fopen("debugfile", "a");
	fprintf(fp, "this is for jacobian corresponding to state vars %d \n", j);
	fprintf(fp,
	    "In corrector time step %d with dt=%f dtdx=%f dtdy=%f kactx=%f , kacty=%f , x=%6f, y=%6f \n state vars are: \n",
	    iter, dt, dtdx, dtdy, *(Curr_El->get_kactxy()), *(Curr_El->get_kactxy() + 1),
	    *(Curr_El->get_coord()), *(Curr_El->get_coord() + 1));

	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", state_vars[state]);
	fprintf(fp, "\n prev state vars are: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", prev_state_vars[state]);
	fprintf(fp, "\n xp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxxp[state]);
	fprintf(fp, "\n xm: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxxm[state]);
	fprintf(fp, "\n yp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxyp[state]);
	fprintf(fp, "\n ym: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", fluxym[state]);
	fprintf(fp, "\n dUdx: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state]);
	fprintf(fp, "\n dUdy: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state + NUM_STATE_VARS]);
	fprintf(fp, "\n");

	fclose(fp);

	cout << "for element " << Curr_El->get_ithelem() << endl;
	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		diff1[ivar] = (fluxxpold[ivar] - fluxxp[ivar]);
		diff2[ivar] = (fluxxmold[ivar] - fluxxm[ivar]);
		diff3[ivar] = (fluxypold[ivar] - fluxyp[ivar]);
		diff4[ivar] = (fluxymold[ivar] - fluxym[ivar]);
		cout << "jacobian of flux of " << ivar << " for neighbor " << effelement << " for state " << j
		    << "  xp: " << (fluxxpold[ivar] - fluxxp[ivar]) / INCREMENT << " xm: "
		    << (fluxxmold[ivar] - fluxxm[ivar]) / INCREMENT << " yp: "
		    << (fluxypold[ivar] - fluxyp[ivar]) / INCREMENT << " ym: "
		    << (fluxymold[ivar] - fluxym[ivar]) / INCREMENT << endl;
	}

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		fluxxold[ivar] = fluxxpold[ivar] - fluxxmold[ivar];
		fluxxnew[ivar] = fluxxp[ivar] - fluxxm[ivar];
		fluxyold[ivar] = fluxypold[ivar] - fluxymold[ivar];
		fluxynew[ivar] = fluxyp[ivar] - fluxym[ivar];
	}

//	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
//
//		abs_fluxx_diff[ivar] = dabs(fluxxold[ivar] - fluxxnew[ivar]);
//		if (abs_fluxx_diff[ivar] > 0.)
//			cout << setw(10) << setprecision(8) << "change effects the flux_x for var  " << ivar
//			    << " eff_elem   " << effelement << "  j  " << j << "  value  " << abs_fluxx_diff[ivar]
//			    << endl;
//
//		abs_fluxy_diff[ivar] = dabs(fluxyold[ivar] - fluxynew[ivar]);
//		if (abs_fluxy_diff[ivar] > 0.)
//			cout << "change effects the flux_y for var  " << ivar << " eff_elem   " << effelement
//			    << "  j  " << j << "  value  " << abs_fluxy_diff[ivar] << endl;
//	}

	return;
}

void record_flux(HashTable* El_Table, HashTable* NodeTable, unsigned* key, MatProps* matprops_ptr,
    int myid, double fluxold[4][NUM_STATE_VARS]) {

	Element* Curr_El = (Element*) (El_Table->lookup(key));

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

	Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

	Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

	Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
		fluxold[0][ivar] = nxp->flux[ivar];
		fluxold[1][ivar] = nyp->flux[ivar];
		fluxold[2][ivar] = nxm->flux[ivar];
		fluxold[3][ivar] = nym->flux[ivar];
	}
}

unsigned* makekey(unsigned k1,unsigned k2){
	 unsigned *key=new unsigned[2];
	 key[0]=k1;key[1]=k2;
	 return key;
}
