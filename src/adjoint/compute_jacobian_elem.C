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

//#define DEBUG

void calc_jacobian_elem(Mat3x3& jacobian, const Mat3x3& jac_flux_n_x, const Mat3x3& jac_flux_p_x,
    const Mat3x3& jac_flux_n_y, const Mat3x3& jac_flux_p_y, double* prev_state_vars,
    double* d_state_vars_x, double* d_state_vars_y, double *curvature, double* gravity,
    double* d_gravity, double* dh_sens, double kact, int effelem, double dtdx, double dtdy,
    double dt, int* stop, double* OrgSgn, double* adjusted_tan_phi_bed,
    double* adjusted_sin_phi_int) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		for (int j = 0; j < NUM_STATE_VARS; ++j)
			jacobian(i, j) = +dtdx * (jac_flux_p_x(i, j) - jac_flux_n_x(i, j))
			    + dtdy * (jac_flux_p_y(i, j) - jac_flux_n_y(i, j));

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		double alpha = adjusted_sin_phi_int[0] * OrgSgn[0] * kact;

		double betta = adjusted_sin_phi_int[1] * OrgSgn[1] * kact;

		if (effelem > 0) {

			if (!stop[0])
				jacobian(1, 0) -= -dt * prev_state_vars[0] * alpha * gravity[2] * dh_sens[1];

			if (!stop[1])
				jacobian(2, 0) -= -dt * prev_state_vars[0] * betta * gravity[2] * dh_sens[0];

		} else {

			double velocity[2] = { prev_state_vars[1] / prev_state_vars[0], prev_state_vars[2]
			    / prev_state_vars[0] };

			double vx_sq = velocity[0] * velocity[0], vy_sq = velocity[1] * velocity[1];

			jacobian(1, 0) -= dt
			    * (gravity[0] - alpha * prev_state_vars[0] * (2 * d_gravity[1] + gravity[2] * dh_sens[1])
			        - alpha * gravity[2] * d_state_vars_y[0]);

			if (OrgSgn[2] != 0. && !stop[0]) {

				jacobian(1, 0) -= -dt * OrgSgn[2] * adjusted_tan_phi_bed[0]
				    * (gravity[2] - vx_sq * curvature[0]);

				jacobian(1, 1) -= -dt * OrgSgn[2] * adjusted_tan_phi_bed[0]
				    * (2 * velocity[0] * curvature[0]);

				jacobian(1, 2) -= 0.;
			}

			jacobian(2, 0) -= dt
			    * (gravity[1] - betta * prev_state_vars[0] * (2 * d_gravity[0] + gravity[2] * dh_sens[0])
			        - betta * gravity[2] * d_state_vars_x[0]);

			if (OrgSgn[3] != 0. && !stop[1] /*&& unitvy != 0.*/) {
				jacobian(2, 0) -= -dt * OrgSgn[3] * adjusted_tan_phi_bed[1]
				 * (gravity[2] - vy_sq * curvature[1]);

				jacobian(2, 1) -= 0.;

				jacobian(2, 2) -= -dt * OrgSgn[3] * adjusted_tan_phi_bed[1]
				 * (2 * velocity[1] * curvature[1]);

			}
		}
	}

#ifdef DEBUG

	for (int i = 0; i < NUM_STATE_VARS; ++i)
	for (int j = 0; j < NUM_STATE_VARS; ++j)
	if (isnan(jacobian(i, j)) || isinf(jacobian(i, j)))
	cout << "hello, I found you \n";

	for (int i = 0; i < NUM_STATE_VARS; ++i)
	for (int j = 0; j < NUM_STATE_VARS; ++j)
	if (dabs(jacobian(i, j)) > 5.) {
		cout << "WARNING for Jacobian \n";
		if (dabs(jacobian(i, j)) > max_jac)
		max_jac = dabs(jacobian(i, j));
		jacobian(i, j) = sign(jacobian(i, j)) * 5;
	}
#endif
}

