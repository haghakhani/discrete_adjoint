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

#define barycentric_interpolation_wrapeper(interpolant,data,i,xi,yi) barycentric_interpolation(interpolant[0].x,\
	    interpolant[1].x, interpolant[2].x,interpolant[0].y, interpolant[1].y, interpolant[2].y,\
		  interpolant[0].data[i], interpolant[1].data[i], interpolant[2].data[i], xi, yi)

#define bilinear_interpolation_wrapper(interpolant,data,i,xi,yi)  bilinear_interpolation(interpolant[0].x,\
				    interpolant[1].x, interpolant[2].x,interpolant[3].x, interpolant[0].y ,\
				    interpolant[1].y, interpolant[2].y, interpolant[3].y, interpolant[0].data[i],\
					  interpolant[1].data[i], interpolant[2].data[i],	interpolant[3].data[i], xi,yi)

#define linearx_interp_wrapper(interpolant,data,i,xi) linear_interp(interpolant[0].x,\
					    interpolant[1].x, interpolant[0].data[i], interpolant[1].data[i],xi)

#define lineary_interp_wrapper(interpolant,data,i,xi) linear_interp(interpolant[0].y,\
					    interpolant[1].y, interpolant[0].data[i], interpolant[1].data[i],xi)

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

//int three=0,four=0,two=0;

struct Neighbor {

	Neighbor(ErrorElem* input) {

		for (int i = 0; i < NUM_STATE_VARS; ++i) {
			// we need to find the properties of its father

			state_vars[i] = *(input->get_state_vars() + i);
			prev_state_vars[i] = *(input->get_prev_state_vars() + i);
			adjoint[i] = *(input->get_adjoint() + i);
		}
		// we need to find the father coordinate

		int which_son = input->get_which_son();
		double* dx = input->get_dx();
		double* xi = input->get_coord();

		switch (which_son) {
			case 0:
				x = xi[0] + .5 * dx[0];
				y = xi[1] + .5 * dx[1];
				break;
			case 1:
				x = xi[0] - .5 * dx[0];
				y = xi[1] + .5 * dx[1];
				break;
			case 2:
				x = xi[0] - .5 * dx[0];
				y = xi[1] - .5 * dx[1];
				break;
			case 3:
				x = xi[0] + .5 * dx[0];
				y = xi[1] - .5 * dx[1];
				break;
			default:
				cerr << "this is not a valid element \n";
		}

	}
	;

	double x, y;
	double state_vars[3], prev_state_vars[3], adjoint[3];

	bool operator <(const Neighbor& a) const {
		if ((x < a.x && (fabs(x - a.x) > 1e-12)) || ((fabs(x - a.x) < 1e-12) && y < a.y))
			return true;
		return false;
	}

	bool operator ==(const Neighbor& a) const {
		if ((fabs(x - a.x) < 1e-12) && (fabs(y - a.y) < 1e-12))
			return true;
		return false;
	}

};

struct NeighLess {
	bool operator()(const Neighbor& a, const Neighbor& b) const {

		if (((a.x < b.x) && dabs(a.x - b.x) > 1e-12) || (dabs(a.x - b.x) < 1e-12 && (a.y < b.y)))
			return true;
		return false;
	}
};

int counter = 0;

void bilinear_interp(HashTable* cp_El_Table) {

	HashEntryPtr currentPtr;
	HashEntryPtr *buck = cp_El_Table->getbucketptr();

	for (int i = 0; i < cp_El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {

				ErrorElem *son[4];

				son[0] = (ErrorElem*) (currentPtr->value);

				if (son[0]->get_adapted_flag() > 0 && son[0]->get_which_son() == 0) {

					counter++;

					ErrorElem *neighbors[12];

					for (int j = 0; j < 12; ++j)
						neighbors[j] = NULL;

					// orientation for the father and sons are same
					int xp = son[0]->get_positive_x_side();
					int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

					for (int j = 1; j < 4; ++j)
						son[j] = (ErrorElem*) cp_El_Table->lookup(son[0]->get_brothers() + j * KEYLENGTH);

					neighbors[3] = son[0]->get_side_neighbor<ErrorElem>(cp_El_Table, 3);

					neighbors[6] = son[0]->get_side_neighbor<ErrorElem>(cp_El_Table, 6);

					neighbors[10] = son[0]->get_side_neighbor<ErrorElem>(cp_El_Table, 10);

					neighbors[0] = son[1]->get_side_neighbor<ErrorElem>(cp_El_Table, 0);

					neighbors[7] = son[1]->get_side_neighbor<ErrorElem>(cp_El_Table, 7);

					neighbors[11] = son[1]->get_side_neighbor<ErrorElem>(cp_El_Table, 11);

					neighbors[8] = son[2]->get_side_neighbor<ErrorElem>(cp_El_Table, 8);

					neighbors[4] = son[2]->get_side_neighbor<ErrorElem>(cp_El_Table, 4);

					neighbors[1] = son[2]->get_side_neighbor<ErrorElem>(cp_El_Table, 1);

					neighbors[5] = son[3]->get_side_neighbor<ErrorElem>(cp_El_Table, 5);

					neighbors[9] = son[3]->get_side_neighbor<ErrorElem>(cp_El_Table, 9);

					neighbors[2] = son[3]->get_side_neighbor<ErrorElem>(cp_El_Table, 2);

					for (int j = 0; j < 4; ++j) {
						switch (j) {

							// it is very important to point that we actually need the father of these elements
							// but we retrieve the father information from them. Because state, prev_state and adjoint
							// of the sons are same as father (from constant reconstruction), and the coord of
							// the father can be computed by knowing the coord of son and knowing which son it is
							case 0: {
								// this condition means that 10 and 6 are same elements
								if (neighbors[6] && neighbors[10] && neighbors[6] == neighbors[10])

									// this condition means that 2 and 6 are same elements
									if (*(son[j]->get_neigh_proc() + xm + 4) < 0)
										bilinear_interp_elem(neighbors[10], neighbors[3], neighbors[9], son[j], son[j]);
									else
										bilinear_interp_elem(neighbors[10], neighbors[3], neighbors[2], son[j], son[j]);

								// this condition means that 10 and 3 are same elements
								else if (neighbors[3] && neighbors[10] && neighbors[3] == neighbors[10])
									// this condition means that 3 and 7 are same elements
									if (*(son[j]->get_neigh_proc() + ym + 4) < 0)
										bilinear_interp_elem(neighbors[10], neighbors[11], neighbors[6], son[j],
										    son[j]);
									else
										bilinear_interp_elem(neighbors[10], neighbors[7], neighbors[6], son[j], son[j]);

								else
									bilinear_interp_elem(neighbors[10], neighbors[3], neighbors[6], son[j], son[j]);
							}
								break;

							case 1: {

								// this condition means that 0 and 11 are same elements
								if (neighbors[11] && neighbors[0] && neighbors[11] == neighbors[0])

									//0 , 4 are same
									if (*(son[j]->get_neigh_proc() + xp + 4) < 0)
										bilinear_interp_elem(neighbors[7], neighbors[11], son[j], neighbors[8], son[j]);
									else
										bilinear_interp_elem(neighbors[7], neighbors[11], son[j], neighbors[4], son[j]);

								// this condition means that 7 and 11 are same elements
								else if (neighbors[7] && neighbors[11] && neighbors[7] == neighbors[11])

									if (*(son[j]->get_neigh_proc() + ym + 4) < 0)
										bilinear_interp_elem(neighbors[10], neighbors[11], son[j], neighbors[0],
										    son[j]);
									else
										bilinear_interp_elem(neighbors[3], neighbors[11], son[j], neighbors[0], son[j]);

								else
									bilinear_interp_elem(neighbors[7], neighbors[11], son[j], neighbors[0], son[j]);

							}
								break;
							case 2: {

								// this condition means that 4 and 8 are same elements
								if (neighbors[4] && neighbors[8] && neighbors[4] == neighbors[8])
									//0 , 4 are same
									if (*(son[j]->get_neigh_proc() + xp + 4) < 0)
										bilinear_interp_elem(son[j], neighbors[11], neighbors[1], neighbors[8], son[j]);
									else
										bilinear_interp_elem(son[j], neighbors[0], neighbors[1], neighbors[8], son[j]);

								// this condition means that 8 and 1 are same elements
								else if (neighbors[8] && neighbors[1] && neighbors[8] == neighbors[1])

									if (*(son[j]->get_neigh_proc() + yp + 4) < 0)
										bilinear_interp_elem(son[j], neighbors[4], neighbors[9], neighbors[8], son[j]);
									else
										bilinear_interp_elem(son[j], neighbors[4], neighbors[5], neighbors[8], son[j]);
								else
									bilinear_interp_elem(son[j], neighbors[4], neighbors[1], neighbors[8], son[j]);

							}
								break;
							case 3: {
								// this condition means that 9 and 2 are same elements
								if (neighbors[9] && neighbors[2] && neighbors[9] == neighbors[2])

									if (*(son[j]->get_neigh_proc() + xm + 4) < 0)
										bilinear_interp_elem(neighbors[10], son[j], neighbors[9], neighbors[5], son[j]);
									else
										bilinear_interp_elem(neighbors[6], son[j], neighbors[9], neighbors[5], son[j]);

								// this condition means that 10 and 3 are same elements
								else if (neighbors[9] && neighbors[5] && neighbors[9] == neighbors[5])

									if (*(son[j]->get_neigh_proc() + yp + 4) < 0)
										bilinear_interp_elem(neighbors[2], son[j], neighbors[9], neighbors[8], son[j]);
									else
										bilinear_interp_elem(neighbors[2], son[j], neighbors[9], neighbors[1], son[j]);
								else
									bilinear_interp_elem(neighbors[2], son[j], neighbors[9], neighbors[5], son[j]);

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

//	cout<<"counter is "<<counter<<" this is something new "<<num_nonzero_elem(cp_El_Table)<<endl;
}

//	cout<<"four: "<< four<<" three: "<< three<<" two: "<<two<<endl;

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
//	if (isnan(interp) || isinf(interp))
//		cout << "it is so sad that I found you" << endl;
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

inline double barycentric_interpolation(double x1, double x2, double x3, double y1, double y2,
    double y3, double f1, double f2, double f3, double x, double y) {

	double inv_det = 1 / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
	double landa1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) * inv_det;
	double landa2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) * inv_det;
	double landa3 = 1 - landa1 - landa2;

	return landa1 * f1 + landa2 * f2 + landa3 * f3;

}

void bilinear_interp_elem(ErrorElem *elem11, ErrorElem *elem21, ErrorElem *elem12,
    ErrorElem *elem22, ErrorElem *Curr_El) {

	double *bi_state_vars = Curr_El->get_bilin_state();
	double *bi_prev_state_vars = Curr_El->get_bilin_prev_state();
	double *bi_adjoint = Curr_El->get_bilin_adj();

	double *x = Curr_El->get_coord();
	ErrorElem* temp[4] = { elem11, elem21, elem12, elem22 };
	vector<Neighbor> input;
	vector<Neighbor> interpolant;

	for (int i = 0; i < 4; ++i)
		if (temp[i])
			input.push_back(Neighbor(temp[i]));

	int add = 0;
	for (int i = 0; i < input.size(); ++i) {
		add = 1;
		for (int j = 0; j < interpolant.size(); ++j)
			if (input[i] == interpolant[j]) {
				add = 0;
				break;
			}
		if (add)
			interpolant.push_back(input[i]);
	}

	switch (interpolant.size()) {
		case 0:
			cout << "ERROR: NO ELEMENT FOR BILINEAR INTERPOLATION " << endl;
			break;

		case 1:
			//father is in corner so there is no element for interp
			//we do not do any extrapolation and leave as it is,
			//which in refinement constructor should be the value of father element
//			cout << "WARNING: ONE ELEMENT FOR BILINEAR INTERPOLATION " << endl;

			for (int i = 0; i < NUM_STATE_VARS; ++i) {

				bi_state_vars[i] = *(Curr_El->get_state_vars() + i);

				bi_prev_state_vars[i] = *(Curr_El->get_prev_state_vars() + i);

				bi_adjoint[i] = *(Curr_El->get_adjoint() + i);

			}

			break;
		case 2: {
//			two++;

			NeighLess customLess;

			sort(interpolant.begin(), interpolant.end(), customLess);

			int xdir = 0, ydir = 1;

			if (interpolant[1].x - interpolant[0].x < 1e-12)
				//interpolation in y
				for (int i = 0; i < NUM_STATE_VARS; ++i) {

					bi_state_vars[i] = lineary_interp_wrapper(interpolant, state_vars, i, x[1]);

					bi_prev_state_vars[i] = lineary_interp_wrapper(interpolant, prev_state_vars, i, x[1]);

					bi_adjoint[i] = lineary_interp_wrapper(interpolant, adjoint, i, x[1]);

				}

			else

				for (int i = 0; i < NUM_STATE_VARS; ++i) {

					bi_state_vars[i] = linearx_interp_wrapper(interpolant, state_vars, i, x[0]);

					bi_prev_state_vars[i] = linearx_interp_wrapper(interpolant, prev_state_vars, i, x[0]);

					bi_adjoint[i] = linearx_interp_wrapper(interpolant, adjoint, i, x[0]);

				}

			break;
		}
		case 3: {
			// this comes out of refinement near the boundary, and it's not possible to find all of the neighbors
			// we do “barycentric interpolation” a.k.a. “convex combination”	a.k.a. “affine linear extension
//			three++;

			for (int i = 0; i < NUM_STATE_VARS; ++i) {

				bi_state_vars[i] = barycentric_interpolation_wrapeper(interpolant, state_vars, i, x[0],
				    x[1]);

				bi_prev_state_vars[i] = barycentric_interpolation_wrapeper(interpolant, prev_state_vars, i,
				    x[0], x[1]);

				bi_adjoint[i] = barycentric_interpolation_wrapeper(interpolant, adjoint, i, x[0], x[1]);

			}
			break;
		}
		case 4: {
//			four++;
			// bilinear interpolation
			for (int i = 0; i < NUM_STATE_VARS; ++i) {

				bi_state_vars[i] = bilinear_interpolation_wrapper(interpolant, state_vars, i, x[0], x[1]);

				bi_prev_state_vars[i] = bilinear_interpolation_wrapper(interpolant, prev_state_vars, i,
				    x[0], x[1]);

				bi_adjoint[i] = bilinear_interpolation_wrapper(interpolant, adjoint, i, x[0], x[1]);

			}
			break;
		}
		default:
			cout << "ERROR IN BILINEAR INTERPOLATION " << endl;
			break;

	}

	for (int i=0;i<NUM_STATE_VARS;++i){
		assert(!isinf(bi_state_vars[i]) && !isinf(bi_prev_state_vars[i]) && !isinf(bi_adjoint[i]));
		assert(!isnan(bi_state_vars[i]) && !isnan(bi_prev_state_vars[i]) && !isnan(bi_adjoint[i]));
	}
}

