/*
 * debug_flux.C
 *
 *  Created on: May 15, 2015
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

void record_flux(DualMesh* dualmesh, DualCell* cell, double fluxold[4][NUM_STATE_VARS]) {

	int yind = *(cell->get_key()), xind = *(cell->get_key() + 1);

	unsigned key_xp[2] = { yind, xind + 1 }, key_xm[2] = { yind, xind - 1 }, key_yp[2] = { yind + 1,
	    xind }, key_ym[2] = { yind - 1, xind };

	for (int ind = 0; ind < NUM_STATE_VARS; ++ind) {

		fluxold[0][ind] = *(cell->get_xflux() + ind);
		fluxold[1][ind] = *(cell->get_yflux() + ind);

		fluxold[2][ind] = *((dualmesh->get_dualcell(key_xm[0], key_xm[1]))->get_xflux() + ind);
		fluxold[3][ind] = *((dualmesh->get_dualcell(key_ym[0], key_ym[1]))->get_yflux() + ind);

	}

}
void flux_debug(DualCell* cell, double fluxold[4][NUM_STATE_VARS], double flux[4][NUM_STATE_VARS],
    int effelement, int j, int iter, double dt, double* dx) {

	double diff1[NUM_STATE_VARS], diff2[NUM_STATE_VARS], diff3[NUM_STATE_VARS], diff4[NUM_STATE_VARS];
	double diff_fluxxold[NUM_STATE_VARS], diff_fluxxnew[NUM_STATE_VARS],
	    diff_fluxyold[NUM_STATE_VARS], diff_fluxynew[NUM_STATE_VARS];
	double abs_fluxx_diff[NUM_STATE_VARS], abs_fluxy_diff[NUM_STATE_VARS];

	const double *state_vars = cell->get_state_vars();
	const double *prev_state_vars = cell->get_prev_state_vars();
	const double *d_state_vars = cell->get_d_state_vars();

	double dtdx = dt / dx[0];
	double dtdy = dt / dx[1];

	FILE *fp;
	fp = fopen("debugfile", "a");
	fprintf(fp, "this is for jacobian corresponding to state vars %d \n", j);
	fprintf(fp,
	    "In jacobian time step %d with dt=%f dtdx=%f dtdy=%f kactx=%f , x=%6f, y=%6f \n state vars are: \n",
	    iter, dt, dtdx, dtdy, cell->get_kact(), *(cell->get_position()),
	    *(cell->get_position() + 1));

	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", state_vars[state]);
	fprintf(fp, "\n prev state vars are: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", prev_state_vars[state]);
	fprintf(fp, "\n xp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", flux[0][state]);
	fprintf(fp, "\n xm: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", flux[2][state]);
	fprintf(fp, "\n yp: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", flux[1][state]);
	fprintf(fp, "\n ym: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", flux[3][state]);
	fprintf(fp, "\n dUdx: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state]);
	fprintf(fp, "\n dUdy: \n");
	for (int state = 0; state < NUM_STATE_VARS; state++)
		fprintf(fp, "%10e ", d_state_vars[state + NUM_STATE_VARS]);
	fprintf(fp, "\n");

	fclose(fp);

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		diff1[ivar] = dabs(fluxold[0][ivar] - flux[0][ivar]);
		diff2[ivar] = dabs(fluxold[2][ivar] - flux[2][ivar]);
		diff3[ivar] = dabs(fluxold[1][ivar] - flux[1][ivar]);
		diff4[ivar] = dabs(fluxold[3][ivar] - flux[3][ivar]);
	}

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		diff_fluxxold[ivar] = fluxold[0][ivar] - fluxold[2][ivar];
		diff_fluxxnew[ivar] = flux[0][ivar] - flux[2][ivar];
		diff_fluxyold[ivar] = fluxold[1][ivar] - fluxold[3][ivar];
		diff_fluxynew[ivar] = flux[1][ivar] - flux[3][ivar];
	}

	for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {

		abs_fluxx_diff[ivar] = dabs(diff_fluxxold[ivar] - diff_fluxxnew[ivar]);

		if (abs_fluxx_diff[ivar] > 0.)
			cout << setw(10) << setprecision(8) << "change effects the flux_x for var  " << ivar
			    << " eff_elem   " << effelement << "  j  " << j << "  value  " << abs_fluxx_diff[ivar]
			    << endl;

		abs_fluxy_diff[ivar] = dabs(diff_fluxyold[ivar] - diff_fluxynew[ivar]);
		if (abs_fluxy_diff[ivar] > 0.)
			cout << "change effects the flux_y for var  " << ivar << " eff_elem   " << effelement
			    << "  j  " << j << "  value  " << abs_fluxy_diff[ivar] << endl;
	}

	return;
}

void flux_debug(Element* Curr_El, double fluxxpold[4], double fluxxmold[4], double fluxypold[4],
    double fluxymold[4], double fluxxp[4], double fluxxm[4], double fluxyp[4], double fluxym[4],
    int effelement, int j, int iter, double dt) {

	int state_num = 4;

	double diff1[state_num], diff2[state_num], diff3[state_num], diff4[state_num];
	double fluxxold[state_num], fluxxnew[state_num], fluxyold[state_num], fluxynew[state_num];
	double abs_fluxx_diff[state_num], abs_fluxy_diff[state_num];

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

	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", state_vars[state]);
	fprintf(fp, "\n prev state vars are: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", prev_state_vars[state]);
	fprintf(fp, "\n xp: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", fluxxp[state]);
	fprintf(fp, "\n xm: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", fluxxm[state]);
	fprintf(fp, "\n yp: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", fluxyp[state]);
	fprintf(fp, "\n ym: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", fluxym[state]);
	fprintf(fp, "\n dUdx: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", d_state_vars[state]);
	fprintf(fp, "\n dUdy: \n");
	for (int state = 0; state < state_num; state++)
		fprintf(fp, "%10e ", d_state_vars[state + NUM_STATE_VARS]);
	fprintf(fp, "\n");

	fclose(fp);

	for (int ivar = 0; ivar < state_num; ivar++) {

		diff1[ivar] = dabs(fluxxpold[ivar] - fluxxp[ivar]);
		diff2[ivar] = dabs(fluxxmold[ivar] - fluxxm[ivar]);
		diff3[ivar] = dabs(fluxypold[ivar] - fluxyp[ivar]);
		diff4[ivar] = dabs(fluxymold[ivar] - fluxym[ivar]);
	}

	for (int ivar = 0; ivar < state_num; ivar++) {

		fluxxold[ivar] = fluxxpold[ivar] - fluxxmold[ivar];
		fluxxnew[ivar] = fluxxp[ivar] - fluxxm[ivar];
		fluxyold[ivar] = fluxypold[ivar] - fluxymold[ivar];
		fluxynew[ivar] = fluxyp[ivar] - fluxym[ivar];
	}

	for (int ivar = 0; ivar < state_num; ivar++) {

		abs_fluxx_diff[ivar] = dabs(fluxxold[ivar] - fluxxnew[ivar]);
		if (abs_fluxx_diff[ivar] > 0. && ivar != 1)
			cout << setw(10) << setprecision(8) << "change effects the flux_x for var  " << ivar
			    << " eff_elem   " << effelement << "  j  " << j << "  value  " << abs_fluxx_diff[ivar]
			    << endl;

		abs_fluxy_diff[ivar] = dabs(fluxxold[ivar] - fluxxnew[ivar]);
		if (abs_fluxy_diff[ivar] > 0. && ivar != 1)
			cout << "change effects the flux_y for var  " << ivar << " eff_elem   " << effelement
			    << "  j  " << j << "  value  " << abs_fluxy_diff[ivar] << endl;
	}

	return;
}

