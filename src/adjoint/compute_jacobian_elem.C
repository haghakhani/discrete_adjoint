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

void calc_jacobian_elem(Mat3x3<double>& jacobian, const Mat3x3<double>& jac_flux_n_x,
    const Mat3x3<double>& jac_flux_p_x, const Mat3x3<double>& jac_flux_n_y,
    const Mat3x3<double>& jac_flux_p_y, double* state_vars, double* d_state_vars_x,
    double* d_state_vars_y, double *curvature, double* gravity, double* d_gravity, double* dh_sens,
    double int_fric, double bedfrict, double kact, int effelem, double dtdx, double dtdy, double dt,
    int* stop, double* OrgSgn) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		for (int j = 0; j < NUM_STATE_VARS; ++j)
			jacobian(i, j) = dtdx * (jac_flux_p_x(i, j) - jac_flux_n_x(i, j))
			    + dtdy * (jac_flux_p_y(i, j) - jac_flux_n_y(i, j));

	if (effelem == 0)
		// this part is for element itself
		for (int i = 0; i < NUM_STATE_VARS; ++i)
			jacobian(i, i) -= 1.;

	if (state_vars[0] > GEOFLOW_TINY) {

		double sin_int_fric = sin(int_fric);

		if (effelem > 0) {

			if (!stop[0])
				jacobian(1, 0) -= dt * state_vars[0] * kact * OrgSgn[0] * sin_int_fric * gravity[2]
				    * dh_sens[1];

			if (!stop[1])
				jacobian(2, 0) -= dt * state_vars[0] * kact * OrgSgn[1] * sin_int_fric * gravity[2]
				    * dh_sens[0];

		} else {

			double velocity[2] = { 0., 0. }, unitvx, unitvy;

			double speed = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]);

			double tan_bed_frict = tan(bedfrict);

			if (speed > 0.) {
				unitvx = velocity[0] / speed;
				unitvy = velocity[1] / speed;
			}

			//effect of prev_state_vars on residual vector
			if (!stop[0]) {

				jacobian(1, 0) -= dt
				    * (gravity[0] - gravity[2] * unitvx * tan_bed_frict
				        - state_vars[0] * kact * OrgSgn[0] * sin_int_fric
				            * (gravity[2] * dh_sens[1] + d_gravity[1])
				        - kact * OrgSgn[0] * sin_int_fric
				            * (gravity[2] * d_state_vars_y[0] + state_vars[0] * d_gravity[1]));

				if (speed == 0.) {
					jacobian(1, 1) -= dt * (-unitvx * tan_bed_frict * curvature[0]);

					jacobian(1, 2) -= 0.;
				} else {
					jacobian(1, 1) -= dt
					    * ((tan_bed_frict * (gravity[2] * state_vars[0] + state_vars[1] * curvature[0])
					        * (unitvx * unitvx / speed)) - unitvx * tan_bed_frict * curvature[0]
					        - (gravity[2] * state_vars[0] + state_vars[1] * curvature[0]) * tan_bed_frict
					            / speed);
					jacobian(1, 2) -= dt
					    * (tan_bed_frict * (gravity[2] * state_vars[0] + state_vars[1] * curvature[0])
					        * (unitvx * unitvy / speed));

				}
			}

			if (!stop[1]) {

				jacobian(2, 0) -= dt
				    * (gravity[1] - gravity[2] * unitvy * tan_bed_frict
				        - state_vars[0] * kact * OrgSgn[1] * sin_int_fric
				            * (gravity[2] * dh_sens[0] + d_gravity[0])
				        - kact * OrgSgn[1] * sin_int_fric
				            * (gravity[2] * d_state_vars_x[0] + state_vars[0] * d_gravity[0]));

				if (speed == 0) {
					jacobian(2, 1) -= 0.;

					jacobian(2, 2) -= dt * (-unitvy * tan_bed_frict * curvature[1]);

				} else {

					jacobian(2, 1) -= dt
					    * (tan_bed_frict * (gravity[2] * state_vars[0] + state_vars[2] * curvature[1])
					        * (unitvx * unitvy / speed));

					jacobian(2, 2) -= dt
					    * ((tan_bed_frict * (gravity[2] * state_vars[0] + state_vars[2] * curvature[1])
					        * (unitvy * unitvy / speed)) - unitvy * tan_bed_frict * curvature[1]
					        - (gravity[2] * state_vars[0] + state_vars[2] * curvature[1]) * tan_bed_frict
					            / speed);
				}

			}
		}
	}
}

