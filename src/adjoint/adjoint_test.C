/*
 * adjoint_test.C
 *
 *  Created on: May 22, 2015
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"
#include <fenv.h>

double max_err1 = 0., max_err2 = 0.;
unsigned key1_1, key2_1, iter_1, key1_2, key2_2, iter_2;

double simple_test(HashTable* El_Table, TimeProps* timeprops, MatProps* matprops_ptr) {

	double dot = 0.;
	vector<pair<unsigned, unsigned> > wrong_elem;
	vector<double> wrong_value, wrong_value1;

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	cout << "results:  \n";

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					double* state_vars_prev = Curr_El->get_prev_state_vars();
					double* adjoint = Curr_El->get_adjoint();
					double* gravity = Curr_El->get_gravity();
					double* curve = Curr_El->get_curvature();
					double vel[2], h_inv, orgSrcSgn[2];
					double unitvx, unitvy;
					double* d_state_vars_x = Curr_El->get_d_state_vars();
					double* d_state_vars_y = (Curr_El->get_d_state_vars() + 3);
					double fric_tiny = matprops_ptr->frict_tiny;
					double kact = *(Curr_El->get_kactxy());
					double* d_grav = Curr_El->get_d_gravity();

					if (state_vars_prev[0] > GEOFLOW_TINY) {
						h_inv = 1. / state_vars_prev[0];
						vel[0] = state_vars_prev[1] / state_vars_prev[0];
						vel[1] = state_vars_prev[2] / state_vars_prev[0];

					} else {
						h_inv = 0.;
						vel[0] = 0.;
						vel[1] = 0.;
						unitvx = unitvy = 0.;

					}

					double speed = sqrt(vel[0] * vel[0] + vel[1] * vel[1]);

					if (speed > 0.) {

						unitvx = vel[0] / speed;
						unitvy = vel[1] / speed;
					}

					double test1 = 0.;

					if (max(gravity[2] * state_vars_prev[0] + vel[0] * state_vars_prev[1] * curve[0], 0.))
						test1 = adjoint[1] * unitvx
						    * (state_vars_prev[0] * gravity[2] + state_vars_prev[1] * vel[0] * curve[0]);

					if (max(gravity[2] * state_vars_prev[0] + vel[1] * state_vars_prev[2] * curve[1], 0.))
						test1 += adjoint[2] * unitvy
						    * (state_vars_prev[0] * gravity[2] + state_vars_prev[2] * vel[1] * curve[1]);

					dot += test1;

					double tmp = h_inv * (d_state_vars_y[1] - vel[0] * d_state_vars_y[0]);
					orgSrcSgn[0] = tiny_sgn(tmp, fric_tiny);
					orgSrcSgn[2] = tiny_sgn(vel[0], fric_tiny);

					tmp = h_inv * (d_state_vars_x[2] - vel[1] * d_state_vars_x[0]);
					orgSrcSgn[1] = tiny_sgn(tmp, fric_tiny);
					orgSrcSgn[3] = tiny_sgn(vel[1], fric_tiny);

					double test2 = state_vars_prev[0] * kact
					    * (adjoint[1] * orgSrcSgn[0]
					        * (gravity[2] * d_state_vars_y[0] + d_grav[1] * state_vars_prev[0])
					        + adjoint[2] * orgSrcSgn[1]
					            * (gravity[2] * d_state_vars_x[0] + d_grav[0] * state_vars_prev[0]));

					// this is for simple test
//					double test1, test2 = 0.;
//					test1 = adjoint[1] + adjoint[2] ;

					if (fabs(test1) > 1e-16 || fabs(test2) > 1e-16) {
						wrong_elem.push_back(make_pair(*(Curr_El->pass_key()), *(Curr_El->pass_key() + 1)));
						wrong_value.push_back(test1);
						wrong_value1.push_back(test2);
//						cout << test1 << " , " << test2 << endl;

						if (fabs(test1) > max_err1) {
							max_err1 = fabs(test1);
							key1_1 = *(Curr_El->pass_key());
							key2_1 = *(Curr_El->pass_key() + 1);
							iter_1 = timeprops->iter;
						} else if (fabs(test2) > max_err2) {
							max_err2 = fabs(test2);
							key1_2 = *(Curr_El->pass_key());
							key2_2 = *(Curr_El->pass_key() + 1);
							iter_2 = timeprops->iter;
						}

					}

				}

				currentPtr = currentPtr->next;
			}
		}

	ofstream f("wrong_elem.txt", ios::app);
	f << "time step: " << timeprops->iter << " size of vector: " << wrong_elem.size() << endl;
	for (int i = 0; i < wrong_elem.size(); ++i)
		f << setw(10) << wrong_elem[i].first << " , " << wrong_elem[i].second << " , " << wrong_value[i]
		    << " , " << wrong_value1[i] << '\n';

	return dot;
}

void perturbU(HashTable* El_Table, PertElemInfo* pelinf, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	const double epsilon = INCREMENT;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (dabs(*(Curr_El->get_coord()) - pelinf->elempos[0]) < epsilon
					    && dabs(*(Curr_El->get_coord() + 1) - pelinf->elempos[1]) < epsilon) {

						// perturbing U
						*(Curr_El->get_state_vars()) += epsilon;

						// to break the for loop
						i = El_Table->get_no_of_buckets();

						// break the while loop
						break;
					}

				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}

void fill_pertelem_info(HashTable* El_Table, PertElemInfo* eleminfo) {

	HashEntryPtr currentPtr;
	DualElem *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	Jacobian *jacobian, *neighjac;
	DualElem *neigh_elem;
	Vec_Mat<9> jacobianmat;
	double *curr_adj_ptr, *prev_adj_ptr;

	const double epsilon = INCREMENT;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (DualElem*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (dabs(*(Curr_El->get_coord()) - eleminfo->elempos[0]) < epsilon
					    && dabs(*(Curr_El->get_coord() + 1) - eleminfo->elempos[1]) < epsilon) {

						eleminfo->h = *(Curr_El->get_state_vars());

						eleminfo->elem_size[0] = *(Curr_El->get_dx());

						eleminfo->elem_size[1] = *(Curr_El->get_dx() + 1);

//						eleminfo->

						jacobianmat = Curr_El->get_jacobian();
						eleminfo->func_sens = *(Curr_El->get_func_sens());

						for (int effelement = 0; effelement < 5; effelement++) { //0 for the element itself, and the rest id for neighbour elements

							if (effelement == 0) {		    //this part of code

								curr_adj_ptr = (Curr_El->get_adjoint());
								prev_adj_ptr = (Curr_El->get_prev_adjoint());

								for (int j = 0; j < 3; ++j) {
									eleminfo->neigh_jac[effelement].curr_adj[j] = curr_adj_ptr[j];
									eleminfo->neigh_jac[effelement].prev_adj[j] = prev_adj_ptr[j];
								}

								jacobianmat = Curr_El->get_jacobian();

								for (int k = 0; k < 3; ++k)
									for (int l = 0; l < 3; ++l)
										eleminfo->neigh_jac[effelement].jacobianMat[k][l] = jacobianmat(effelement, k,
										    l);

							} else {

								neigh_elem = Curr_El->get_side_neighbor<DualElem>(El_Table, effelement - 1);//basically we are checking all neighbor elements, and start from xp neighbor
								if (neigh_elem) {

									curr_adj_ptr = (neigh_elem->get_state_vars() + 6);
									prev_adj_ptr = (neigh_elem->get_prev_state_vars() + 6);

									for (int j = 0; j < 3; ++j) {
										eleminfo->neigh_jac[effelement].curr_adj[j] = curr_adj_ptr[j];
										eleminfo->neigh_jac[effelement].prev_adj[j] = prev_adj_ptr[j];
									}

									jacobianmat = neigh_elem->get_jacobian();

									int jacind;

									switch (effelement) {
										case 1:	//in xp neighbor I have to read jacobian of xm, because position of curr_el is in xm side of that neighbor
											jacind = 3;
											break;
										case 2:		    //for yp return ym
											jacind = 4;
											break;
										case 3:		    //for xm return xp
											jacind = 1;
											break;
										case 4:		    //for ym return yp
											jacind = 2;
											break;
										default:
											cout << "invalid neighbor position" << endl;
									}

									for (int k = 0; k < 3; ++k)
										for (int l = 0; l < 3; ++l)
											eleminfo->neigh_jac[effelement].jacobianMat[k][l] = jacobianmat(effelement, k,
											    l);

								}
							}
						}

						// to break the for loop
						i = El_Table->get_no_of_buckets();

						// break the while loop
						break;
					}
				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}
void find_test_elem(HashTable* El_Table, PertElemInfo** pelinf, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	double maxh = 0.;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->get_state_vars()) > maxh)
						maxh = *(Curr_El->get_state_vars());

				}
				currentPtr = currentPtr->next;
			}
		}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->get_state_vars()) == maxh) {

						double* pos = Curr_El->get_coord();

						(*pelinf) = new PertElemInfo(pos, maxh, iter);

						// to break the for loop
						i = El_Table->get_no_of_buckets();

						// break the while loop
						break;
					}

				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}

int check_elem_exist(HashTable *El_Table, unsigned *key) {
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
					if (*(Curr_El->pass_key()) == key[0] && *(Curr_El->pass_key() + 1) == key[1])
						return (1);

				currentPtr = currentPtr->next;
			}
		}

	return (0);
}

int checkElement(HashTable *El_Table, HashTable *NodeTable, double *max, unsigned *key) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	int gg = 0;

	double fluxold[4][NUM_STATE_VARS];

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->pass_key()) == key[0] && *(Curr_El->pass_key() + 1) == key[1]) {
						int xp = Curr_El->get_positive_x_side(); //finding the direction of element
						int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

						Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

						Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

						Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

						Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

//						for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
//							fluxold[0][ivar] = nxp->flux[ivar];
//							fluxold[1][ivar] = nyp->flux[ivar];
//							fluxold[2][ivar] = nxm->flux[ivar];
//							fluxold[3][ivar] = nym->flux[ivar];
//						}
						cout << "flux yp  " << nyp->flux[0] << endl;
						gg = 1;
					}
				}
				currentPtr = currentPtr->next;
			}
		}

//	cout << " \n flux xp: \n";
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		cout << " " << fluxold[0][i];
//
//	cout << " \n flux yp: \n";
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		cout << " " << fluxold[1][i];
//
//	cout << " \n flux xm: \n";
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		cout << " " << fluxold[2][i];
//
//	cout << " \n flux ym: \n";
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		cout << " " << fluxold[3][i];

	return (gg);
}

void check_elem_size(HashTable *El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

//	double min_dx[] = { 10000, 10000 };
//	int min_gen = 10000;
//
//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//				if (Curr_El->get_adapted_flag() > 0) {
//					if (Curr_El->get_gen() < min_gen) {
//						min_gen = Curr_El->get_gen();
//						min_dx[0] = Curr_El->get_dx()[0];
//						min_dx[1] = Curr_El->get_dx()[1];
//					}
//
//				}
//				currentPtr = currentPtr->next;
//			}
//		}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					int gen = Curr_El->get_gen();
					double dif_gen = gen - min_gen;
					if (!(Curr_El->get_dx()[0] == min_dx[0] * pow(.5, dif_gen)))
						cout << "this element does not pass the first test \n";
					if (!(Curr_El->get_dx()[1] == min_dx[1] * pow(.5, dif_gen)))
						cout << "this element does not pass the second test \n";
				}
				currentPtr = currentPtr->next;
			}
		}

}

void myround(double *num) {
	fesetround(FE_TONEAREST);
	*num = rint(*num * 1.e8) / 1.e8;
}

void compute_dx(HashTable *El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					double dif_gen = Curr_El->get_gen() - min_gen;
					Curr_El->get_dx()[0] = min_dx[0] * pow(.5, dif_gen);
					Curr_El->get_dx()[1] = min_dx[1] * pow(.5, dif_gen);
					if (Curr_El->get_dx()[0] == 0. || Curr_El->get_dx()[1] == 0.)
						cout << "error in dx\n";

				}
				currentPtr = currentPtr->next;
			}
		}
}

