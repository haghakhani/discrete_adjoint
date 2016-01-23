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
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
		double *b, int *ldb, int *info);

struct Elnear {

	double x, y;
	Elnear(double xin, double yin) {
		x = xin;
		y = yin;
	}
	bool operator()(Element* a, Element* b) {

		double *xa = a->get_coord();
		double *xb = b->get_coord();

		if (((xa[0] - x) * (xa[0] - x) + (xa[1] - y) * (xa[1] - y))
				< ((xb[0] - x) * (xb[0] - x) + (xb[1] - y) * (xb[1] - y)))
			return true;
		return false;
	}
};

struct Elnear_x {

	bool operator()(Element* a, const double x) {
		return *(a->get_coord()) < x;
	}

	bool operator()(const double x, Element* a) {
		return x < *(a->get_coord());
	}
};

struct Elnear_y {

	bool operator()(Element* a, const double x) {
		return *(a->get_coord() + 1) < x;
	}

	bool operator()(const double x, Element* a) {
		return x < *(a->get_coord() + 1);
	}
};

void bilinear_interp(HashTable* El_Table, HashTable* cp_El_Table) {

	HashEntryPtr currentPtr;
	HashEntryPtr *buck = cp_El_Table->getbucketptr();

	for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				Element* son[4];

				son[0] = (Element*) (currentPtr->value);

				if (son[0]->get_adapted_flag() > 0 && son[0]->get_which_son() == 0) {

					for (int j = 1; j < 4; ++j)
						son[j] = (Element*) cp_El_Table->lookup(
								son[0]->get_brothers() + j * KEYLENGTH);

					Element * father = (Element*) El_Table->lookup(son[0]->getfather());

					assert(father);

					vector<Element*> neighbors;
					neighbors.push_back(father);

					// following section let to just save valid and not repeated neighbors
					int add = 0;
					for (int j = 0; j < 12; ++j) {
						Element* neigh = father->get_side_neighbor(El_Table, j);
						if (neigh) {
							add = 1;
							for (int k = 0; k < neighbors.size(); ++k)
								if (neigh == neighbors[k]) {
									add = 0;
									break;
								}
							if (add)
								neighbors.push_back(neigh);

						}
					}

					for (int j = 0; j < 4; ++j)
						bilinear_interp_elem(neighbors, son[j]);

				}
				currentPtr = currentPtr->next;
			}
		}
}

double bilinear_interp_value(double x1, double x2, double y1, double y2,
		double f11, double f21, double f12, double f22, double xinterp,
		double yinterp, int type) {

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
		interp = (f11 * (x2 - xinterp) * (y2 - yinterp)
				+ f21 * (xinterp - x1) * (y2 - yinterp)
				+ f12 * (x2 - xinterp) * (yinterp - y1)
				+ f22 * (xinterp - x1) * (yinterp - y1)) / ((x2 - x1) * (y2 - y1));
		break;
	case 1:
		//interpolation in y
		//no other modification on bilinear interpolation is required,
		//since we call the function such that put zero for those elements that do not exist
		//we just removed (x2 - x1) from denominator
		interp = (f11 * (x2 - xinterp) * (y2 - yinterp)
				+ f21 * (xinterp - x1) * (y2 - yinterp)
				+ f12 * (x2 - xinterp) * (yinterp - y1)
				+ f22 * (xinterp - x1) * (yinterp - y1)) / (y2 - y1);
		break;
	case 2:
		//interpolation in x
		//no other modification on bilinear interpolation is required,
		//since we call the function such that put zero for those elements that do not exist
		//we just removed (y2 - y1) from denominator
		interp = (f11 * (x2 - xinterp) * (y2 - yinterp)
				+ f21 * (xinterp - x1) * (y2 - yinterp)
				+ f12 * (x2 - xinterp) * (yinterp - y1)
				+ f22 * (xinterp - x1) * (yinterp - y1)) / (x2 - x1);
		break;

	default:
		cout << "not a valid type in interp_value function " << endl;
	}
//	if (isnan(interp) || isinf(interp))
//		cout << "it is so sad that I found you" << endl;
	return (interp);
}

inline double linear_interp(double x1, double x2, double f1, double f2,
		double x) {

	return (f2 - f1) / (x2 - x1) * (x - x1);
}

double bilinear_interpolation(double x1, double x2, double x3, double x4,
		double y1, double y2, double y3, double y4, double f1, double f2, double f3,
		double f4, double xin, double yin) {

//  P=a0+a1*x+a2*y+a3*x*y
//	Ab=x
//	[1,x0,y0,x0y0][a0] [phi0]
//	[1,x1,y1,x1y1][a1] [phi1]
//	[1,x2,y2,x2y2][a2]=[phi2]
//	[1,x3,y3,x3y3][a3] [phi3]

	int dim = 4, one = 1, info, ipiv[dim];

	double A[dim * dim], x[] = { x1, x2, x3, x4 }, y[] = { y1, y2, y3, y4 },
			coef[] = { f1, f2, f3, f4 };
	if (f1 == 0. && f2 == 0. && f3 == 0. && f4 == 0.)
		return 0.;
	else {

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

		if (info != 0)
			cout << "solution is not correct " << info << endl;

		return coef[0] + coef[1] * xin + coef[2] * yin + coef[3] * xin * yin;
	}

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

inline double barycentric_interpolation(double x1, double x2, double x3,
		double y1, double y2, double y3, double f1, double f2, double f3, double x,
		double y) {

	double inv_det = 1 / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
	double landa1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) * inv_det;
	double landa2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) * inv_det;
	double landa3 = 1 - landa1 - landa2;

	return landa1 * f1 + landa2 * f2 + landa3 * f3;

}

void bilinear_interp_elem(vector<Element*>& interpolant, Element *Curr_El) {

	double *state_vars = Curr_El->get_state_vars();
	double *prev_state_vars = Curr_El->get_prev_state_vars();
	double *adjoint = Curr_El->get_adjoint();
	double *prev_adjoint = Curr_El->get_prev_adjoint();
	double *x = Curr_El->get_coord();

	Elnear_x near_x;
	Elnear_y near_y;

	vector<Element*>::iterator el0 = lower_bound(interpolant.begin(),
			interpolant.end(), *(Curr_El->get_coord()), near_x);
	vector<Element*>::iterator el1 = upper_bound(interpolant.begin(),
			interpolant.end(), *(Curr_El->get_coord()), near_x);

	vector<Element*>::iterator el2 = lower_bound(interpolant.begin(),
			interpolant.end(), *(Curr_El->get_coord() + 1), near_y);
	vector<Element*>::iterator el3 = upper_bound(interpolant.begin(),
			interpolant.end(), *(Curr_El->get_coord() + 1), near_y);

	if (!(*el0) || !(*el1) || !(*el2) || !(*el3)) {
		// we do “barycentric interpolation” a.k.a. “convex combination”	a.k.a. “affine linear extension
		// this should happen just when father is in the corner
//		cout << "ERROR: NO ELEMENT FOR BILINEAR INTERPOLATION " << endl;
		for (int i = 0; i < NUM_STATE_VARS; ++i) {
			state_vars[i] = barycentric_interpolation(*(interpolant[0]->get_coord()),
					*(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
					*(interpolant[0]->get_coord() + 1),
					*(interpolant[1]->get_coord() + 1),
					*(interpolant[2]->get_coord() + 1),
					*(interpolant[0]->get_state_vars() + i),
					*(interpolant[1]->get_state_vars() + i),
					*(interpolant[2]->get_state_vars() + i), x[0], x[1]);

			prev_state_vars[i] = barycentric_interpolation(
					*(interpolant[0]->get_coord()), *(interpolant[1]->get_coord()),
					*(interpolant[2]->get_coord()), *(interpolant[0]->get_coord() + 1),
					*(interpolant[1]->get_coord() + 1),
					*(interpolant[2]->get_coord() + 1),
					*(interpolant[0]->get_prev_state_vars() + i),
					*(interpolant[1]->get_prev_state_vars() + i),
					*(interpolant[2]->get_prev_state_vars() + i), x[0], x[1]);

			adjoint[i] = barycentric_interpolation(*(interpolant[0]->get_coord()),
					*(interpolant[1]->get_coord()), *(interpolant[2]->get_coord()),
					*(interpolant[0]->get_coord() + 1),
					*(interpolant[1]->get_coord() + 1),
					*(interpolant[2]->get_coord() + 1),
					*(interpolant[0]->get_adjoint() + i),
					*(interpolant[1]->get_adjoint() + i),
					*(interpolant[2]->get_adjoint() + i), x[0], x[1]);

			prev_adjoint[i] = barycentric_interpolation(
					*(interpolant[0]->get_coord()), *(interpolant[1]->get_coord()),
					*(interpolant[2]->get_coord()), *(interpolant[0]->get_coord() + 1),
					*(interpolant[1]->get_coord() + 1),
					*(interpolant[2]->get_coord() + 1),
					*(interpolant[0]->get_prev_adjoint() + i),
					*(interpolant[1]->get_prev_adjoint() + i),
					*(interpolant[2]->get_prev_adjoint() + i), x[0], x[1]);
		}

	} else {

		//		cout << interpolant.size() << endl;
		for (int i = 0; i < NUM_STATE_VARS; ++i) {

			state_vars[i] = bilinear_interpolation(*((*el0)->get_coord()),
					*((*el1)->get_coord()), *((*el2)->get_coord()),
					*((*el3)->get_coord()), *((*el0)->get_coord() + 1),
					*((*el1)->get_coord() + 1), *((*el2)->get_coord() + 1),
					*((*el3)->get_coord() + 1), *((*el0)->get_state_vars() + i),
					*((*el1)->get_state_vars() + i), *((*el2)->get_state_vars() + i),
					*((*el3)->get_state_vars() + i), x[0], x[1]);

			prev_state_vars[i] = bilinear_interpolation(*((*el0)->get_coord()),
					*((*el1)->get_coord()), *((*el2)->get_coord()),
					*((*el3)->get_coord()), *((*el0)->get_coord() + 1),
					*((*el1)->get_coord() + 1), *((*el2)->get_coord() + 1),
					*((*el3)->get_coord() + 1), *((*el0)->get_prev_state_vars() + i),
					*((*el1)->get_prev_state_vars() + i),
					*((*el2)->get_prev_state_vars() + i),
					*((*el3)->get_prev_state_vars() + i), x[0], x[1]);

			adjoint[i] = bilinear_interpolation(*((*el0)->get_coord()),
					*((*el1)->get_coord()), *((*el2)->get_coord()),
					*((*el3)->get_coord()), *((*el0)->get_coord() + 1),
					*((*el1)->get_coord() + 1), *((*el2)->get_coord() + 1),
					*((*el3)->get_coord() + 1), *((*el0)->get_adjoint() + i),
					*((*el1)->get_adjoint() + i), *((*el2)->get_adjoint() + i),
					*((*el3)->get_adjoint() + i), x[0], x[1]);

			prev_adjoint[i] = bilinear_interpolation(*((*el0)->get_coord()),
					*((*el1)->get_coord()), *((*el2)->get_coord()),
					*((*el3)->get_coord()), *((*el0)->get_coord() + 1),
					*((*el1)->get_coord() + 1), *((*el2)->get_coord() + 1),
					*((*el3)->get_coord() + 1), *((*el0)->get_prev_adjoint() + i),
					*((*el1)->get_prev_adjoint() + i), *((*el2)->get_prev_adjoint() + i),
					*((*el3)->get_prev_adjoint() + i), x[0], x[1]);
		}

	}
}
