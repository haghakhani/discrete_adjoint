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
 * Author: Hossein Aghakhani
 * Date: Jan 9 2015
 *******************************************************************
 */
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define	DEBUG

void residual(double* residual, double *state_vars, double *prev_state_vars, //3
    double *fluxxp, double *fluxyp, double *fluxxm, double *fluxym, double dtdx, //5
    double dtdy, double dt, double *d_state_vars_x, double *d_state_vars_y, //4
    double *curvature, double intfrictang, double bedfrict, double *gravity, //4
    double *dgdx, double kactxyelem, double fric_tiny, double* orgSrcSgn, //4
    double increment, double epsilon, int* check_stop_crit, int srcflag, int org_res_flag) {

	double velocity[DIMENSION];
	double kactxy[DIMENSION];
	double tmp[NUM_STATE_VARS];
	//double bedfrict;

	if (prev_state_vars[0] > GEOFLOW_TINY) {
		for (int k = 0; k < DIMENSION; k++)
			kactxy[k] = kactxyelem;

		if ((prev_state_vars[1] == 0. && prev_state_vars[2] == increment)
		    || (prev_state_vars[2] == 0. && prev_state_vars[1] == increment)) {
			velocity[0] = 0.;
			velocity[1] = 0.;
		} else {

			// fluid velocities
			velocity[0] = prev_state_vars[1] / prev_state_vars[0];
			velocity[1] = prev_state_vars[2] / prev_state_vars[0];
		}

	} else {
		for (int k = 0; k < DIMENSION; k++) {
			kactxy[k] = epsilon;
			velocity[k] = 0.;
		}
		//bedfrict = bedfrictin;
	}

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		residual[i] = 0.0;
		tmp[i] = 0.;
	}

	for (int i = 0; i < NUM_STATE_VARS; i++)
		tmp[i] = prev_state_vars[i] - dtdx * (fluxxp[i] - fluxxm[i]) - dtdy * (fluxyp[i] - fluxym[i]);

	if (prev_state_vars[0] > GEOFLOW_TINY && srcflag) {

		double unitvx = 0., unitvy = 0., h_inv = 0., speed = 0.;

		speed = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]);

		if (speed > 0.) {
			unitvx = velocity[0] / speed;
			unitvy = velocity[1] / speed;
		}

		//x dir
		double s1 = gravity[0] * prev_state_vars[0];

		double sin_int_fric = sin(intfrictang);
		double s2 = orgSrcSgn[0] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_y[0] + dgdx[1] * prev_state_vars[0]) * sin_int_fric;

		double tan_bed_fric = tan(bedfrict);
		double s3 = unitvx
		    * max(gravity[2] * prev_state_vars[0] + velocity[0] * prev_state_vars[1] * curvature[0],
		        0.0) * tan_bed_fric;
		if (prev_state_vars[1] == increment)
			s3 = 0;

		if (dabs(tmp[1] + dt * s1) > dabs(dt * (s2 + s3)) && !check_stop_crit[0])
			tmp[1] += dt * (s1 - s2 - s3);
		else {
			tmp[1] = 0.;
			if (org_res_flag)
				check_stop_crit[0] = 1;
		}

		//y dir

		s1 = gravity[1] * prev_state_vars[0];

		s2 = orgSrcSgn[1] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_x[0] + dgdx[0] * prev_state_vars[0]) * sin_int_fric;

		s3 = unitvy
		    * max(gravity[2] * prev_state_vars[0] + velocity[1] * prev_state_vars[2] * curvature[1],
		        0.0) * tan_bed_fric;

		if (prev_state_vars[2] == increment)
			s3 = 0;

		if (dabs(tmp[2] + dt * s1) > dabs(dt * (s2 + s3)) && !check_stop_crit[1])
			tmp[2] += dt * (s1 - s2 - s3);
		else {
			tmp[2] = 0.;
			if (org_res_flag)
				check_stop_crit[1] = 1;
		}

	}

	if (org_res_flag)
		for (int i = 0; i < NUM_STATE_VARS; i++)
			state_vars[i] = tmp[i];

	for (int i = 0; i < NUM_STATE_VARS; i++)
		residual[i] = state_vars[i] - tmp[i];

#ifdef DEBUG

//	for (int k = 0; k < 3; k++)
//		if (residual[k] > 1e-5) {
//			cout << "something that has to be checked" << endl << flush;
//			exit(-2);
//		}

	for (int k = 0; k < 3; k++)
		if (isnan(residual[k])) {
			cout << "exit for NAN in residual" << endl << flush;
			exit(-1);
		}

	for (int k = 0; k < 3; k++)
		if (isinf(residual[k])) {
			cout << "exit for Inf in residual" << endl << flush;
			exit(-2);
		}
#endif

	return;
}

double tiny_sgn(double num, double tiny) {
	if (dabs(num) < tiny)
		return 0.;
	else if (num > tiny)
		return 1.;
	else
		return -1.;
}

void update_states(double *state_vars, double *prev_state_vars, //2
    double *fluxxp, double *fluxyp, double *fluxxm, double *fluxym, double dtdx, //5
    double dtdy, double dt, double *d_state_vars_x, double *d_state_vars_y, //4
    double *curvature, double intfrictang, double bedfrict, double *gravity, //4
    double *dgdx, double kactxyelem, double fric_tiny, int* stop, double* orgSrcSgn) { //5

	double velocity[DIMENSION];
	double kactxy[DIMENSION];
	double tmp;

	if (prev_state_vars[0] > GEOFLOW_TINY) {
		for (int k = 0; k < DIMENSION; k++)
			kactxy[k] = kactxyelem;

		// velocities
		velocity[0] = prev_state_vars[1] / prev_state_vars[0];
		velocity[1] = prev_state_vars[2] / prev_state_vars[0];

	} else {

		for (int k = 0; k < DIMENSION; k++) {
			kactxy[k] = 0.;
			velocity[k] = 0.;
		}

	}

	for (int i = 0; i < NUM_STATE_VARS; i++)
		state_vars[i] = prev_state_vars[i] - dtdx * (fluxxp[i] - fluxxm[i])
		    - dtdy * (fluxyp[i] - fluxym[i]);

	if (state_vars[0] < 0.)
		state_vars[0] = 0;

	/// I have added this part
//	double cte=1.,bte=1.;
//
//	state_vars[1] += dt * cte;
//	state_vars[2] += dt * bte;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

//		double unitvx = 0., unitvy = 0., h_inv = 0., speed = 0.;
		double h_inv = 0.;

		h_inv = 1. / prev_state_vars[0];

		tmp = h_inv * (d_state_vars_y[1] - velocity[0] * d_state_vars_y[0]);
		orgSrcSgn[0] = tiny_sgn(tmp, fric_tiny);
		orgSrcSgn[2] = tiny_sgn(velocity[0], fric_tiny);

		tmp = h_inv * (d_state_vars_x[2] - velocity[1] * d_state_vars_x[0]);
		orgSrcSgn[1] = tiny_sgn(tmp, fric_tiny);
		orgSrcSgn[3] = tiny_sgn(velocity[1], fric_tiny);

//		speed = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]);
//
//		if (speed > 0.) {
//			unitvx = velocity[0] / speed;
//			unitvy = velocity[1] / speed;
//		}

		//x dir
		double s1 = gravity[0] * prev_state_vars[0];

		double sin_int_fric = sin(intfrictang);
		double s2 = orgSrcSgn[0] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_y[0] + dgdx[1] * prev_state_vars[0]) * sin_int_fric;

		double tan_bed_fric = tan(bedfrict);
		double s3 = orgSrcSgn[2]
		    * max(gravity[2] * prev_state_vars[0] + velocity[0] * prev_state_vars[1] * curvature[0],
		        0.0) * tan_bed_fric;

		if (s3 == 0. && orgSrcSgn[2])
			stop[0] = 1;

		state_vars[1] += dt * (s1 - s2 - s3);

//		if (dabs(state_vars[1] + dt * s1) > dabs(dt * (s2 + s3)))
//			state_vars[1] += dt * (s1 - s2 - s3);
//		else {
//			state_vars[1] = 0.;
//			stop[0] = 1;
//		}

		//y dir
		s1 = gravity[1] * prev_state_vars[0];

		s2 = orgSrcSgn[1] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_x[0] + dgdx[0] * prev_state_vars[0]) * sin_int_fric;

		s3 = orgSrcSgn[3]
		    * max(gravity[2] * prev_state_vars[0] + velocity[1] * prev_state_vars[2] * curvature[1],
		        0.0) * tan_bed_fric;

		if (s3 == 0. && orgSrcSgn[3])
			stop[1] = 1;

		state_vars[2] += dt * (s1 - s2 - s3);

//		if (dabs(state_vars[2] + dt * s1) > dabs(dt * (s2 + s3)))
//			state_vars[2] += dt * (s1 - s2 - s3);
//		else {
//			state_vars[2] = 0.;
//			stop[1] = 1;
//		}

	}

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		assert(!isnan(state_vars[i]) && !isinf(state_vars[i]));

}

void orgSourceSgn(Element* Curr_El, double frictiny, double* orgSgn) {

	double* d_state_vars_x = Curr_El->get_d_state_vars();
	double* d_state_vars_y = d_state_vars_x + NUM_STATE_VARS;
	double* prev_state_vars = Curr_El->get_prev_state_vars();
	double h_inv;
	double tmp = 0.0;
	double velocity[2];
	for (int i = 0; i < 2; i++)
		orgSgn[i] = 0.0;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		velocity[0] = prev_state_vars[1] / prev_state_vars[0];
		velocity[1] = prev_state_vars[2] / prev_state_vars[0];

	} else {
		for (int k = 0; k < DIMENSION; k++)
			velocity[k] = 0.;
	}

	if (prev_state_vars[0] > 0.0)
		h_inv = 1. / prev_state_vars[0];

	tmp = h_inv * (d_state_vars_y[1] - velocity[0] * d_state_vars_y[0]);
	orgSgn[0] = tiny_sgn(tmp, frictiny);

	tmp = h_inv * (d_state_vars_x[2] - velocity[1] * d_state_vars_x[0]);
	orgSgn[1] = tiny_sgn(tmp, frictiny);

}

void residual(double *state_vars, double *prev_state_vars, double *fluxxp, double *fluxyp,
    double *fluxxm, double *fluxym, double dtdx, double dtdy, double dt, double *d_state_vars_x,
    double *d_state_vars_y, double *curvature, double intfrictang, double bedfrict, double *gravity,
    double *dgdx, double kactxyelem, double fric_tiny, int* stop, double* orgSrcSgn, int iter,
    double *pre3_state, double* adjusted_tan_phi_bed, double* adjusted_sin_phi_int) {

	double coef = 0.;
	if (iter > 2)
		coef = 0.25;

	double res_vec[] = { 0., 0., 0. };
	for (int i = 0; i < NUM_STATE_VARS; i++)
		res_vec[i] = -dtdx * (fluxxp[i] - fluxxm[i]) - dtdy * (fluxyp[i] - fluxym[i]);

	//multi-step 3rd order TVD time scheme p512 Lecture notes in Comp. Phys.
	state_vars[0] = 0.75 * prev_state_vars[0] + 1.5 * res_vec[0] + coef * pre3_state[0];

	if (state_vars[0] < 0.)
		state_vars[0] = 0.;

	double velocity[DIMENSION];
	double kactxy[DIMENSION];
	double tmp;

	stop[0] = stop[1] = 0;

	for (int i = 0; i < 4; ++i)
		orgSrcSgn[i] = 0.;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		for (int k = 0; k < DIMENSION; k++)
			kactxy[k] = kactxyelem;
		// velocities
		velocity[0] = prev_state_vars[1] / prev_state_vars[0];
		velocity[1] = prev_state_vars[2] / prev_state_vars[0];

		double h_inv = 1. / prev_state_vars[0];

		tmp = h_inv * (d_state_vars_y[1] - velocity[0] * d_state_vars_y[0]);
		orgSrcSgn[0] = tiny_sgn(tmp, fric_tiny);
		orgSrcSgn[2] = tiny_sgn(velocity[0], fric_tiny);

		tmp = h_inv * (d_state_vars_x[2] - velocity[1] * d_state_vars_x[0]);
		orgSrcSgn[1] = tiny_sgn(tmp, fric_tiny);
		orgSrcSgn[3] = tiny_sgn(velocity[1], fric_tiny);

		//x dir
		double s1 = gravity[0] * prev_state_vars[0];

		double sin_int_fric = sin(intfrictang);
		double s2 = orgSrcSgn[0] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_y[0] + dgdx[1] * prev_state_vars[0]) * sin_int_fric;

		double tan_bed_fric = tan(bedfrict);
		double s3 = orgSrcSgn[2]
		    * max(gravity[2] * prev_state_vars[0] + velocity[0] * prev_state_vars[1] * curvature[0],
		        0.0) * tan_bed_fric;

		if (s3 == 0. && orgSrcSgn[2])
			stop[0] = 1;

		if (fabs(1.5 * dt * (s2 + s3))
		    > fabs(1.5 * (res_vec[1] + dt * s1) + .75 * prev_state_vars[1] + coef * pre3_state[1])) {
			//friction are big enough to stop the flow
			state_vars[1] = 0.;
			//we adjust the friction terms to handle this situation
			double fract = fabs(
			    1.5 * (res_vec[1] + dt * s1) + .75 * prev_state_vars[1] + coef * pre3_state[1])
			    / fabs(1.5 * dt * (s2 + s3));

			adjusted_tan_phi_bed[0] = fract * tan_bed_fric;
			adjusted_sin_phi_int[0] = fract * sin_int_fric;

		} else {
			res_vec[1] += dt * (s1 - s2 - s3);
			state_vars[1] = 0.75 * prev_state_vars[1] + 1.5 * res_vec[1] + coef * pre3_state[1];
			adjusted_tan_phi_bed[0] = tan_bed_fric;
			adjusted_sin_phi_int[0] = sin_int_fric;
		}

		//y dir
		s1 = gravity[1] * prev_state_vars[0];

		s2 = orgSrcSgn[1] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_x[0] + dgdx[0] * prev_state_vars[0]) * sin_int_fric;

		s3 = orgSrcSgn[3]
		    * max(gravity[2] * prev_state_vars[0] + velocity[1] * prev_state_vars[2] * curvature[1],
		        0.0) * tan_bed_fric;

		if (s3 == 0. && orgSrcSgn[3])
			stop[1] = 1;

		if (fabs(1.5 * dt * (s2 + s3))
		    > fabs(1.5 * (res_vec[2] + dt * s1) + .75 * prev_state_vars[2] + coef * pre3_state[2])) {
			//friction are big enough to stop the flow
			state_vars[2] = 0.;
			//we have to adjust the frictions otherwise we will have problem in jacobian computation may be just one is enough
			double fract = fabs(
			    1.5 * (res_vec[2] + dt * s1) + .75 * prev_state_vars[2] + coef * pre3_state[2])
			    / fabs(1.5 * dt * (s2 + s3));

			adjusted_tan_phi_bed[1] = fract * tan_bed_fric;
			adjusted_sin_phi_int[1] = fract * sin_int_fric;

		} else {
			res_vec[2] += dt * (s1 - s2 - s3);
			state_vars[2] = 0.75 * prev_state_vars[2] + 1.5 * res_vec[2] + coef * pre3_state[2];
			adjusted_tan_phi_bed[1] = tan_bed_fric;
			adjusted_sin_phi_int[1] = sin_int_fric;
		}

	} else {
		state_vars[1] = 0.75 * prev_state_vars[1] + 1.5 * res_vec[1] + coef * pre3_state[1];
		state_vars[2] = 0.75 * prev_state_vars[2] + 1.5 * res_vec[2] + coef * pre3_state[2];
	}

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		assert(!isnan(res_vec[i]) && !isinf(res_vec[i]));

}

