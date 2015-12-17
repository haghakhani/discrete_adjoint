/*
 * compute_jacobian_elem.C
 *
 *  Created on: Aug 31, 2015
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define DEBUG

void calc_jacobian_elem(Mat3x3& jacobian, const Mat3x3& jac_flux_n_x, const Mat3x3& jac_flux_p_x,
    const Mat3x3& jac_flux_n_y, const Mat3x3& jac_flux_p_y, double* prev_state_vars,
    double* d_state_vars_x, double* d_state_vars_y, double *curvature, double* gravity,
    double* d_gravity, double* dh_sens, double int_fric, double bedfrict, double kact, int effelem,
    double dtdx, double dtdy, double dt, int* stop, double* OrgSgn) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		for (int j = 0; j < NUM_STATE_VARS; ++j)
			jacobian(i, j) = +dtdx * (jac_flux_p_x(i, j) - jac_flux_n_x(i, j))
			    + dtdy * (jac_flux_p_y(i, j) - jac_flux_n_y(i, j));

	if (effelem == 0)
		// this part is for element itself
		for (int i = 0; i < NUM_STATE_VARS; ++i)
			jacobian(i, i) -= 1.;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		double sin_int_fric = sin(int_fric);

		if (effelem > 0) {

			if (!stop[0])
				jacobian(1, 0) -= -dt * prev_state_vars[0] * kact * OrgSgn[0] * sin_int_fric * gravity[2]
				    * dh_sens[1];

			if (!stop[1])
				jacobian(2, 0) -= -dt * prev_state_vars[0] * kact * OrgSgn[1] * sin_int_fric * gravity[2]
				    * dh_sens[0];

		} else {

			double velocity[2] = { prev_state_vars[1] / prev_state_vars[0], prev_state_vars[2]
			    / prev_state_vars[0] }, unitvx = 0., unitvy = 0., speed_inv = 0.;

			double vx_sq = velocity[0] * velocity[0], vy_sq = velocity[1] * velocity[1];

			double speed = sqrt(vx_sq + vy_sq);

			double tan_bed_frict = tan(bedfrict);

			if (speed > 0.) {
				speed_inv = 1. / speed;
				unitvx = velocity[0] * speed_inv;
				unitvy = velocity[1] * speed_inv;
			}

			//effect of prev_state_vars on residual vector
			if (!stop[0]) {

				double alpha = sin_int_fric * OrgSgn[0] * kact;

				jacobian(1, 0) -= dt
				    * (gravity[0]
				        - alpha * prev_state_vars[0] * (2 * d_gravity[1] + gravity[2] * dh_sens[1])
				        - alpha * gravity[2] * d_state_vars_y[0]
				        - unitvx * tan_bed_frict * (gravity[2] - vx_sq * curvature[0]));

				if (speed > 0. && unitvx != 0.) {

					jacobian(1, 1) -=
					    dt * speed_inv * tan_bed_frict
					        * ((unitvx * unitvx - 3.) * curvature[0] * vx_sq
					            + (unitvx * unitvx - 1.) * gravity[2]);

					jacobian(1, 2) -= dt * speed_inv * tan_bed_frict * unitvx * unitvy
					    * (gravity[2] + vx_sq * curvature[0]);
				}
			}

			if (!stop[1]) {

				double beta = sin_int_fric * OrgSgn[1] * kact;

				jacobian(2, 0) -= dt
				    * (gravity[1] - beta * prev_state_vars[0] * (2 * d_gravity[0] + gravity[2] * dh_sens[0])
				        - beta * gravity[2] * d_state_vars_x[0]
				        - unitvy * tan_bed_frict * (gravity[2] - vy_sq * curvature[1]));

				if (speed > 0. && unitvy != 0.) {
					jacobian(2, 1) -= dt * speed_inv * tan_bed_frict * unitvx * unitvy
					    * (gravity[2] + vy_sq * curvature[1]);

					jacobian(2, 2) -=
					    dt * speed_inv * tan_bed_frict
					        * ((unitvy * unitvy - 3.) * curvature[1] * vy_sq
					            + (unitvy * unitvy - 1.) * gravity[2]);

				}
			}
		}
	}

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		for (int j = 0; j < NUM_STATE_VARS; ++j)
			if (isnan(jacobian(i, j)) || isinf(jacobian(i, j)))
				cout << "hello, I found you \n";
#ifdef DEBUG

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		for (int j = 0; j < NUM_STATE_VARS; ++j)
			if (dabs(jacobian(i, j)) > 5.) {
				cout << "WARNING for Jacobian \n";
				if (dabs(jacobian(i, j)) > max_jac)
					max_jac = dabs(jacobian(i, j));
				jacobian(i, j)=sign(jacobian(i, j))*5;
			}
#endif
}

