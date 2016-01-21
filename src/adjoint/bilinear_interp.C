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
 *******************************************************************
 */
//Jan 30, 2015
//haghakha
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define KEY0 3777862041
#define KEY1 2576980374
//#define DEBUG
#include <algorithm>
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb,
    int *info);

struct ElLess {
	bool operator()(Element* a, Element* b) {

		if ((*(a->get_coord()) < *(b->get_coord())
		    && dabs(*(a->get_coord()) - *(b->get_coord())) > 1e-12)
		    || dabs(*(a->get_coord()) - *(b->get_coord())) < 1e-12
		        && *(a->get_coord() + 1) < *(b->get_coord() + 1))
			return true;
		return false;
	}
};

void bilinear_interp(HashTable* El_Table, HashTable* cp_El_Table) {

	HashEntryPtr currentPtr;
	Element *father;
	HashEntryPtr *buck = cp_El_Table->getbucketptr();
	int which_son;
	Element* son[4];
	Element* neighbors[12];

	for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				son[0] = (Element*) (currentPtr->value);

				if (son[0]->get_adapted_flag() > 0 && son[0]->get_which_son() == 0) {

					father = (Element*) El_Table->lookup(son[0]->getfather());

					assert(father);

					for (int j = 0; j < 12; ++j)
						neighbors[j] = father->get_side_neighbor(El_Table, j);

					int xp = father->get_positive_x_side();
					int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

					// this condition means that 10 and 6 are same elements
					if (neighbors[6] && neighbors[10]
					    && *(neighbors[6]->get_coord() + 1) == *(neighbors[10]->get_coord() + 1))

						// this condition means that 2 and 6 are same elements
						if (*(father->get_neigh_proc() + xm + 4) < 0)
							bilinear_interp_elem(neighbors[10], neighbors[3], neighbors[9], father, son[0]);
						else
							bilinear_interp_elem(neighbors[10], neighbors[3], neighbors[2], father, son[0]);

					// this condition means that 10 and 3 are same elements
					else if (neighbors[3] && neighbors[10]
					    && *(neighbors[3]->get_coord()) == *(neighbors[10]->get_coord()))
						// this condition means that 3 and 7 are same elements
						if (*(father->get_neigh_proc() + ym + 4) < 0)
							bilinear_interp_elem(neighbors[10], neighbors[11], neighbors[6], father, son[0]);
						else
							bilinear_interp_elem(neighbors[10], neighbors[7], neighbors[6], father, son[0]);

					else
						bilinear_interp_elem(neighbors[10], neighbors[3], neighbors[6], father, son[0]);

					for (int j = 1; j < 4; ++j) {
						son[j] = (Element*) cp_El_Table->lookup(son[0]->get_brothers() + j * KEYLENGTH);

						switch (j) {

							case 1: {

								// this condition means that 0 and 11 are same elements
								if (neighbors[11] && neighbors[0]
								    && *(neighbors[11]->get_coord() + 1) == *(neighbors[0]->get_coord() + 1))

									//0 , 4 are same
									if (*(father->get_neigh_proc() + xp + 4) < 0)
										bilinear_interp_elem(neighbors[7], neighbors[11], father, neighbors[8], son[j]);
									else
										bilinear_interp_elem(neighbors[7], neighbors[11], father, neighbors[4], son[j]);

								// this condition means that 7 and 11 are same elements
								else if (neighbors[7] && neighbors[11]
								    && *(neighbors[7]->get_coord()) == *(neighbors[11]->get_coord()))

									if (*(father->get_neigh_proc() + ym + 4) < 0)
										bilinear_interp_elem(neighbors[10], neighbors[11], father, neighbors[0],
										    son[j]);
									else
										bilinear_interp_elem(neighbors[3], neighbors[11], father, neighbors[0], son[j]);

								else
									bilinear_interp_elem(neighbors[7], neighbors[11], father, neighbors[0], son[j]);

							}
								break;
							case 2: {

								// this condition means that 4 and 8 are same elements
								if (neighbors[4] && neighbors[8]
								    && *(neighbors[4]->get_coord() + 1) == *(neighbors[8]->get_coord() + 1))
									//0 , 4 are same
									if (*(father->get_neigh_proc() + xp + 4) < 0)
										bilinear_interp_elem(father, neighbors[11], neighbors[1], neighbors[8], son[j]);
									else
										bilinear_interp_elem(father, neighbors[0], neighbors[1], neighbors[8], son[j]);

								// this condition means that 8 and 1 are same elements
								else if (neighbors[8] && neighbors[1]
								    && *(neighbors[8]->get_coord()) == *(neighbors[1]->get_coord()))

									if (*(father->get_neigh_proc() + yp + 4) < 0)
										bilinear_interp_elem(father, neighbors[4], neighbors[9], neighbors[8], son[j]);
									else
										bilinear_interp_elem(father, neighbors[4], neighbors[5], neighbors[8], son[j]);
								else
									bilinear_interp_elem(father, neighbors[4], neighbors[1], neighbors[8], son[j]);

							}
								break;
							case 3: {
								// this condition means that 9 and 2 are same elements
								if (neighbors[9] && neighbors[2]
								    && *(neighbors[9]->get_coord() + 1) == *(neighbors[2]->get_coord() + 1))

									if (*(father->get_neigh_proc() + xm + 4) < 0)
										bilinear_interp_elem(neighbors[10], father, neighbors[9], neighbors[5], son[j]);
									else
										bilinear_interp_elem(neighbors[6], father, neighbors[9], neighbors[5], son[j]);

								// this condition means that 10 and 3 are same elements
								else if (neighbors[9] && neighbors[5]
								    && *(neighbors[9]->get_coord()) == *(neighbors[5]->get_coord()))

									if (*(father->get_neigh_proc() + yp + 4) < 0)
										bilinear_interp_elem(neighbors[2], father, neighbors[9], neighbors[8], son[j]);
									else
										bilinear_interp_elem(neighbors[2], father, neighbors[9], neighbors[1], son[j]);
								else
									bilinear_interp_elem(neighbors[2], father, neighbors[9], neighbors[5], son[j]);

							}
								break;
							default:
								cout << "incorrect son please check me" << endl;

						}
					}

				}
				currentPtr = currentPtr->next;
			}
		}
}

double bilinear_interp_value(double x1, double x2, double y1, double y2, double f11, double f21,
    double f12, double f22, double xinterp, double yinterp, int type) {

//	xmin = x1
//	xmax = x2
//	ymin = y1
//	ymax = y2
//
//	f12--------f22
//	|           |
//	|           |
//	|  interp   |
//	|           |
//	|           |
//	f11---------f21

	double interp = 0.0;
	switch (type) {
		case 0:			//for bilinear interpolation
			interp = (f11 * (x2 - xinterp) * (y2 - yinterp) + f21 * (xinterp - x1) * (y2 - yinterp)
			    + f12 * (x2 - xinterp) * (yinterp - y1) + f22 * (xinterp - x1) * (yinterp - y1))
			    / ((x2 - x1) * (y2 - y1));
			break;
		case 1:
			//interpolation in y
			//no other modification on bilinear interpolation is required,
			//since we call the function such that put zero for those elements that do not exist
			//we just removed (x2 - x1) from denominator
			interp = (f11 * (x2 - xinterp) * (y2 - yinterp) + f21 * (xinterp - x1) * (y2 - yinterp)
			    + f12 * (x2 - xinterp) * (yinterp - y1) + f22 * (xinterp - x1) * (yinterp - y1))
			    / (y2 - y1);
			break;
		case 2:
			//interpolation in x
			//no other modification on bilinear interpolation is required,
			//since we call the function such that put zero for those elements that do not exist
			//we just removed (y2 - y1) from denominator
			interp = (f11 * (x2 - xinterp) * (y2 - yinterp) + f21 * (xinterp - x1) * (y2 - yinterp)
			    + f12 * (x2 - xinterp) * (yinterp - y1) + f22 * (xinterp - x1) * (yinterp - y1))
			    / (x2 - x1);
			break;

		default:
			cout << "not a valid type in interp_value function " << endl;
	}
	if (isnan(interp) || isinf(interp))
		cout << "it is so sad that I found you" << endl;
	return (interp);
}

inline double linear_interp(double x1, double x2, double f1, double f2, double x) {

	return (f2 - f1) / (x2 - x1) * (x - x1);
}

double bilinear_interpolation(double x1, double x2, double x3, double x4, double y1, double y2,
    double y3, double y4, double f1, double f2, double f3, double f4, double xin, double yin) {

//  P=a0+a1*x+a2*y+a3*x*y
//	Ab=x
//	[1,x0,y0,x0y0][a0] [phi0]
//	[1,x1,y1,x1y1][a1] [phi1]
//	[1,x2,y2,x2y2][a2]=[phi2]
//	[1,x3,y3,x3y3][a3] [phi3]

	int dim = 4, one = 1, info, ipiv[dim];

	double A[dim * dim], x[] = { x1, x2, x3, x4 }, y[] = { y1, y2, y3, y4 }, coef[] =
	    { f1, f2, f3, f4 };

	for (int i = 0; i < dim; ++i) {
		A[i] = 1;
		A[i + 4] = x[i];
		A[i + 8] = y[i];
		A[i + 12] = x[i] * y[i];
	}

//	cout << "Matrix A and vector phi" << endl;
//	for (int i = 0; i < dim; ++i) {
//		cout << "[ " << A[i] << " , " << A[i + 4] << " , " << A[i + 8] << " , " << A[i + 12] << " ]";
//		cout << "[ " << phi[i] << " ]" << endl;
//	}

	dgesv_(&dim, &one, A, &dim, ipiv, coef, &dim, &info);

	return coef[0] + coef[1] * xin + coef[2] * yin + coef[3] * xin * yin;

//	cout << "Solution" << endl;
//	for (int i = 0; i < dim; ++i)
//		cout << "[ " << phi[i] << " ]" << endl;

}

//inline double bilinear_interpolation(double x1, double x2, double y1, double y2, double f11,
//    double f12, double f21, double f22, double x, double y) {
//
//	//	xmin = x1
//	//	xmax = x2
//	//	ymin = y1
//	//	ymax = y2
//	//
//	//	f12--------f22
//	//	|           |
//	//	|           |
//	//	|  interp   |
//	//	|           |
//	//	|           |
//	//	f11---------f21
//
//	return (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1) * (y2 - y) + f12 * (x2 - x) * (y - y1)
//	    + f22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1));
//}

inline double barycentric_interpolation(double x1, double x2, double x3, double y1, double y2,
    double y3, double f1, double f2, double f3, double x, double y) {

	double inv_det = 1 / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
	double landa1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) * inv_det;
	double landa2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) * inv_det;
	double landa3 = 1 - landa1 - landa2;

	return landa1 * f1 + landa2 * f2 + landa3 * f3;

}

void bilinear_interp_elem(Element *elem11, Element *elem21, Element *elem12, Element *elem22,
    Element *Curr_El) {

	double *state_vars = Curr_El->get_state_vars();
	double *prev_state_vars = Curr_El->get_prev_state_vars();
	double *adjoint = Curr_El->get_adjoint();
	double *prev_adjoint = Curr_El->get_prev_adjoint();
	double *x = Curr_El->get_coord();
	Element* temp[4] = { elem11, elem21, elem12, elem22 };
	vector<Element*> interpolant;
	int add = 0;

	for (int i = 0; i < 4; ++i)
		if (temp[i]) {
			add = 1;
			for (int j = 0; j < interpolant.size(); ++j)
				if (temp[i] == interpolant[j]) {
					add = 0;
					break;
				}
			if (add)
				interpolant.push_back(temp[i]);
		}

	switch (interpolant.size()) {
		case 0:
			cout << "ERROR: NO ELEMENT FOR BILINEAR INTERPOLATION " << endl;
			break;

		case 1:
			//father is in corner so there is no element for interp
			//we do not do any extrapolation and leave as it is,
			//which in refinement constructor should be the value of father element

			break;
		case 2: {

			ElLess customLess;

			sort(interpolant.begin(), interpolant.end(), customLess);

			if (*(interpolant[1]->get_coord()) - *(interpolant[0]->get_coord()) < 1e-12)
				//interpolation in y
				for (int i = 0; i < NUM_STATE_VARS; ++i) {

					state_vars[i] = linear_interp(*(interpolant[0]->get_coord() + 1),
					    *(interpolant[1]->get_coord() + 1), *(interpolant[0]->get_state_vars() + i),
					    *(interpolant[1]->get_state_vars() + i), x[1]);

					prev_state_vars[i] = linear_interp(*(interpolant[0]->get_coord() + 1),
					    *(interpolant[1]->get_coord() + 1), *(interpolant[0]->get_prev_state_vars() + i),
					    *(interpolant[1]->get_prev_state_vars() + i), x[1]);

					adjoint[i] = linear_interp(*(interpolant[0]->get_coord() + 1),
					    *(interpolant[1]->get_coord() + 1), *(interpolant[0]->get_adjoint() + i),
					    *(interpolant[1]->get_adjoint() + i), x[1]);

					prev_adjoint[i] = linear_interp(*(interpolant[0]->get_coord() + 1),
					    *(interpolant[1]->get_coord() + 1), *(interpolant[0]->get_prev_adjoint() + i),
					    *(interpolant[1]->get_prev_adjoint() + i), x[1]);

				}

			else

				for (int i = 0; i < NUM_STATE_VARS; ++i) {

					state_vars[i] = linear_interp(*(interpolant[0]->get_coord()),
					    *(interpolant[1]->get_coord()), *(interpolant[0]->get_state_vars() + i),
					    *(interpolant[1]->get_state_vars() + i), x[0]);

					prev_state_vars[i] = linear_interp(*(interpolant[0]->get_coord()),
					    *(interpolant[1]->get_coord()), *(interpolant[0]->get_prev_state_vars() + i),
					    *(interpolant[1]->get_prev_state_vars() + i), x[0]);

					adjoint[i] = linear_interp(*(interpolant[0]->get_coord()), *(interpolant[1]->get_coord()),
					    *(interpolant[0]->get_adjoint() + i), *(interpolant[1]->get_adjoint() + i), x[0]);

					prev_adjoint[i] = linear_interp(*(interpolant[0]->get_coord()),
					    *(interpolant[1]->get_coord()), *(interpolant[0]->get_prev_adjoint() + i),
					    *(interpolant[1]->get_prev_adjoint() + i), x[0]);
				}

			break;
		}
		case 3: {
			// this comes out of refinement near the boundary, and it's not possible to find all of the neighbors
			// we do “barycentric interpolation” a.k.a. “convex combination”	a.k.a. “affine linear extension

			for (int i = 0; i < NUM_STATE_VARS; ++i) {
				state_vars[i] = barycentric_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[0]->get_coord() + 1), *(interpolant[1]->get_coord() + 1),
				    *(interpolant[2]->get_coord() + 1), *(interpolant[0]->get_state_vars() + i),
				    *(interpolant[1]->get_state_vars() + i), *(interpolant[2]->get_state_vars() + i), x[0],
				    x[1]);

				prev_state_vars[i] = barycentric_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[0]->get_coord() + 1), *(interpolant[1]->get_coord() + 1),
				    *(interpolant[2]->get_coord() + 1), *(interpolant[0]->get_prev_state_vars() + i),
				    *(interpolant[1]->get_prev_state_vars() + i),
				    *(interpolant[2]->get_prev_state_vars() + i), x[0], x[1]);

				adjoint[i] = barycentric_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[0]->get_coord() + 1), *(interpolant[1]->get_coord() + 1),
				    *(interpolant[2]->get_coord() + 1), *(interpolant[0]->get_adjoint() + i),
				    *(interpolant[1]->get_adjoint() + i), *(interpolant[2]->get_adjoint() + i), x[0], x[1]);

				prev_adjoint[i] = barycentric_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[0]->get_coord() + 1), *(interpolant[1]->get_coord() + 1),
				    *(interpolant[2]->get_coord() + 1), *(interpolant[0]->get_prev_adjoint() + i),
				    *(interpolant[1]->get_prev_adjoint() + i), *(interpolant[2]->get_prev_adjoint() + i),
				    x[0], x[1]);
			}
			break;
		}
		case 4: {
			// bilinear interpolation
			for (int i = 0; i < NUM_STATE_VARS; ++i) {

				state_vars[i] = bilinear_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[3]->get_coord()), *(interpolant[0]->get_coord() + 1),
				    *(interpolant[1]->get_coord() + 1), *(interpolant[2]->get_coord() + 1),
				    *(interpolant[3]->get_coord() + 1), *(interpolant[0]->get_state_vars() + i),
				    *(interpolant[1]->get_state_vars() + i), *(interpolant[2]->get_state_vars() + i),
				    *(interpolant[3]->get_state_vars() + i), x[0], x[1]);

				prev_state_vars[i] = bilinear_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[3]->get_coord()), *(interpolant[0]->get_coord() + 1),
				    *(interpolant[1]->get_coord() + 1), *(interpolant[2]->get_coord() + 1),
				    *(interpolant[3]->get_coord() + 1), *(interpolant[0]->get_prev_state_vars() + i),
				    *(interpolant[1]->get_prev_state_vars() + i),
				    *(interpolant[2]->get_prev_state_vars() + i),
				    *(interpolant[3]->get_prev_state_vars() + i), x[0], x[1]);

				adjoint[i] = bilinear_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[3]->get_coord()), *(interpolant[0]->get_coord() + 1),
				    *(interpolant[1]->get_coord() + 1), *(interpolant[2]->get_coord() + 1),
				    *(interpolant[3]->get_coord() + 1), *(interpolant[0]->get_adjoint() + i),
				    *(interpolant[1]->get_adjoint() + i), *(interpolant[2]->get_adjoint() + i),
				    *(interpolant[3]->get_adjoint() + i), x[0], x[1]);

				prev_adjoint[i] = bilinear_interpolation(*(interpolant[0]->get_coord()),
				    *(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
				    *(interpolant[3]->get_coord()), *(interpolant[0]->get_coord() + 1),
				    *(interpolant[1]->get_coord() + 1), *(interpolant[2]->get_coord() + 1),
				    *(interpolant[3]->get_coord() + 1), *(interpolant[0]->get_prev_adjoint() + i),
				    *(interpolant[1]->get_prev_adjoint() + i), *(interpolant[2]->get_prev_adjoint() + i),
				    *(interpolant[3]->get_prev_adjoint() + i), x[0], x[1]);
			}
			break;
		}
		default:
			cout << "ERROR IN BILINEAR INTERPOLATION " << endl;
			break;

	}

//	double *state_vars, *elem_state[4], *prev_state_vars, *elem_prev_state[4], *adjoint,
//	    *elem_adjoint[4], *coord, *elem_coord[4];
//
//	int type = 0;			//this is just a flag that indicates the type of element
//
//	state_vars = Curr_El->get_state_vars();
//	prev_state_vars = Curr_El->get_prev_state_vars();
//	coord = Curr_El->get_coord();
//	adjoint = Curr_El->get_adjoint();
//
//	for (int i = 0; i < 4; ++i) {
//		elem_state[i] = elem[i]->get_state_vars();
//		elem_prev_state[i] = elem[i]->get_prev_state_vars();
//		elem_adjoint[i] = elem[i]->get_adjoint();
//		elem_coord[i] = elem[i]->get_coord();
//	}
//
//	for (int i = 0; i < 4; ++i) {
//		if (elem[i])
//			for (int j = 0; j < NUM_STATE_VARS; j++)
//				if ( isnan(elem_state[i][j]) || isinf(elem_state[i][j]))
//					cout << "input to interpolation is not NAN or INF" << endl;
//	}
//
//
//
//	if (elem11 && elem12 && elem21 && elem22) {
//		//this is an ordinary case for an element inside the domain
//		//type = 0; we initialized type=0
//
//		for (int j = 0; j < NUM_STATE_VARS; j++) {
//
//			state_vars[j] = bilinear_interp_value(elem_coord[0][0], elem21_coord[0], elem21_coord[1],
//			    elem22_coord[1], elem11_state[j], elem21_state[j], elem12_state[j], elem22_state[j],
//			    coord[0], coord[1], type);
//
//			prev_state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], elem21_coord[1],
//			    elem22_coord[1], elem11_prev_state[j], elem21_prev_state[j], elem12_prev_state[j],
//			    elem22_prev_state[j], coord[0], coord[1], type);
//
//			adjoint[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], elem21_coord[1],
//			    elem22_coord[1], elem11_adjoint[j], elem21_adjoint[j], elem12_adjoint[j],
//			    elem22_adjoint[j], coord[0], coord[1], type);
//		}
//	} else if ((!elem11 && !elem12 && elem21 && elem22)	//interpolation only in y
//	|| (elem11 && elem12 && !elem21 && !elem22)	//left or right side of father is boundary
//	    ) {
//		type = 1;
//		if (elem11) {	// in this case elem21 & elem22 do not exist, so we replace their value with zero
//			for (int j = 0; j < NUM_STATE_VARS; j++) {
//				state_vars[j] = bilinear_interp_value(0,
//				    0,	//interpolation is in y, so x position is not important
//				    elem11_coord[1], elem12_coord[1], elem11_state[j], 0, elem12_state[j], 0, coord[0],
//				    coord[1], type);
//
//				prev_state_vars[j] = bilinear_interp_value(0,
//				    0,	//interpolation is in y, so x position is not important
//				    elem11_coord[1], elem12_coord[1], elem11_prev_state[j], 0, elem12_prev_state[j], 0,
//				    coord[0], coord[1], type);
//
//				adjoint[j] = bilinear_interp_value(0,
//				    0,	//interpolation is in y, so x position is not important
//				    elem11_coord[1], elem12_coord[1], elem11_adjoint[j], 0, elem12_adjoint[j], 0, coord[0],
//				    coord[1], type);
//			}
//
//		} else {	// in this case elem11 & elem12 do not exist, so we replace their values with zero
//
//			for (int j = 0; j < NUM_STATE_VARS; j++) {
//				state_vars[j] = bilinear_interp_value(0,
//				    0,	//interpolation is in y, so x position is not important
//				    elem21_coord[1], elem22_coord[1], 0, elem21_state[j], 0, elem22_state[j], coord[0],
//				    coord[1], type);
//
//				prev_state_vars[j] = bilinear_interp_value(0,
//				    0,	//interpolation is in y, so x position is not important
//				    elem21_coord[1], elem22_coord[1], 0, elem21_prev_state[j], 0, elem22_prev_state[j],
//				    coord[0], coord[1], type);
//
//				adjoint[j] = bilinear_interp_value(0,
//				    0,	//interpolation is in y, so x position is not important
//				    elem21_coord[1], elem22_coord[1], 0, elem21_adjoint[j], 0, elem22_adjoint[j], coord[0],
//				    coord[1], type);
//			}
//		}
//	} else if ((!elem11 && elem12 && !elem21 && elem22)	//interpolation only in x
//	|| (elem11 && !elem12 && elem21 && !elem22)	//top or bottom side of father is boundary
//	    ) {
//		type = 2;
//
//		if (elem11) {	// in this case elem12 & elem22 do not exist, so we replace their value with zero
//			for (int j = 0; j < NUM_STATE_VARS; j++) {
//				state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], 0, 0,	//interpolation is in x, so y position is not important
//				    elem11_state[j], elem21_state[j], 0, 0, coord[0], coord[1], type);
//
//				prev_state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], 0, 0,//interpolation is in x, so y position is not important
//				    elem11_prev_state[j], elem21_prev_state[j], 0, 0, coord[0], coord[1], type);
//
//				adjoint[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], 0, 0,//interpolation is in x, so y position is not important
//				    elem11_adjoint[j], elem21_adjoint[j], 0, 0, coord[0], coord[1], type);
//			}
//
//		} else {	// in this case elem11 & elem21 do not exist, so we replace their value with zero
//
//			for (int j = 0; j < NUM_STATE_VARS; j++) {
//				state_vars[j] = bilinear_interp_value(elem12_coord[0], elem22_coord[0], 0, 0,	//interpolation is in x, so y position is not important
//				    0, 0, elem12_state[j], elem22_state[j], coord[0], coord[1], type);
//
//				prev_state_vars[j] = bilinear_interp_value(elem12_coord[0], elem22_coord[0], 0, 0,//interpolation is in x, so y position is not important
//				    0, 0, elem12_prev_state[j], elem22_prev_state[j], coord[0], coord[1], type);
//
//				adjoint[j] = bilinear_interp_value(elem12_coord[0], elem22_coord[0], 0, 0,//interpolation is in x, so y position is not important
//				    0, 0, elem12_adjoint[j], elem22_adjoint[j], coord[0], coord[1], type);
//			}
//		}
//
//	} else if ((elem11 && !elem12 && !elem21 && !elem22)//father is in corner so there is no element for interp
//	|| (!elem11 && elem12 && !elem21 && !elem22)	      //we do not do any extrapolation and leave as it is,
//	    || (!elem11 && !elem12 && elem21 && !elem22)    //which in refinement constructor should be the value of father element
//	    || (!elem11 && !elem12 && !elem21 && elem22)) {
//
//		//do not do anything
//
//	} else {
//		// this comes out of refinement near the boundary, and it's not possible to find all of the neighbors
//		// we do “barycentric interpolation” a.k.a. “convex combination”	a.k.a. “affine linear extension
//
//		*(Curr_El->get_residual() + 1) = 200;
////		cout << "something is wrong in this configuration" << endl;
////		*(Curr_El->get_state_vars()) = 50;
//
//	}

}
