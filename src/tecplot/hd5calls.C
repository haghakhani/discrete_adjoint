/*
 * =====================================================================================
 *
 *       Filename:  hd5calls.c : If you don't know anything about HDF5
 *                               THEN LEAVE THIS FILE ALONE
 *
 *    Description:  This file contains function calls to hdf5 library
 *                  and provides simpler titan specific interface
 *
 *        Created:  05/17/2007 09:36:27 AM EDT
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
using namespace std;
#include<vector>

#include "GMFG_hdfconstants.h"

hid_t GH5_openfile(const char *filename, char mode) {
	hid_t fid, gid;
	herr_t status;
	// if its a new file
	if (mode == 'n') {
		fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (fid < 0) {
			fprintf(stderr,
			    "Unable to open new hdf file.\
	       If %s already exists  remove it or move it to a different directory\n",
			    filename);
			exit(0);
		}
	}
	// if the file already exists
	else if (mode == 'o') {
		fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (fid < 0)
			fprintf(stderr, "Unable to open %s. Make sure its in the current directory\n", filename);
		return fid;
	} else {
		fprintf(stderr, "GH5 ERROR: Unknown file access option, EXITING\n");
		exit(1);
	}
	//Create base groups for data sets
	gid = H5Gcreate(fid, "/Mesh", 0);
	// this should not happen
	if (gid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Mesh\"\n");
		exit(1);
	}
	status = H5Gclose(gid);
	gid = H5Gcreate(fid, "/Properties", 0);
	if (gid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Properties\"\n");
		exit(1);
	}
	status = H5Gclose(gid);
	return fid;
}

hid_t GH5_open_sol_file(const char *filename, char mode) {
	hid_t fid, gid;
	herr_t status;
	// if its a new file
	if (mode == 'n') {
		fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (fid < 0) {
			fprintf(stderr,
			    "Unable to open new hdf file.\
	       If %s already exists  remove it or move it to a different directory\n",
			    filename);
			exit(0);
		}
	}
	// if the file already exists
	else if (mode == 'o') {
		fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (fid < 0)
			fprintf(stderr, "Unable to open %s. Make sure its in the current directory\n", filename);
		return fid;
	} else {
		fprintf(stderr, "GH5 ERROR: Unknown file access option, EXITING\n");
		exit(1);
	}
	//Create base groups for data sets
	gid = H5Gcreate(fid, "/ZERO_CELLS", 0);
	// this should not happen
	if (gid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Mesh\"\n");
		exit(1);
	}
	status = H5Gclose(gid);
	gid = H5Gcreate(fid, "/NON_ZERO_CELLS", 0);
	if (gid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Properties\"\n");
		exit(1);
	}
	status = H5Gclose(gid);
	return fid;
}

hid_t GH5_open_zip_sol_file(const char *filename, char mode) {
	hid_t fid, gid;
	herr_t status;
	// if its a new file
	if (mode == 'n') {
		hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);

		//The number of chunk slots in the raw data chunk cache for this dataset
		size_t nslots = 521;

		//The total size of the raw data chunk cache for this dataset
		size_t nbytes = 1 * 1024 * 1024;

		//The chunk preemption policy for this dataset
		double w0 = .75;

		H5Pset_cache(fapl, 0, nslots, nbytes, w0);
		fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
		if (fid < 0) {
			fprintf(stderr,
			    "Unable to open new hdf file.\
	       If %s already exists  remove it or move it to a different directory\n",
			    filename);
			exit(0);
		}
	}
	// if the file already exists
	else if (mode == 'o') {
		fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (fid < 0)
			fprintf(stderr, "Unable to open %s. Make sure its in the current directory\n", filename);
		return fid;
	} else {
		fprintf(stderr, "GH5 ERROR: Unknown file access option, EXITING\n");
		exit(1);
	}
	//Create base groups for data sets
	gid = H5Gcreate(fid, "/ZERO_CELLS", 0);
	// this should not happen
	if (gid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Mesh\"\n");
		exit(1);
	}
	status = H5Gclose(gid);
	gid = H5Gcreate(fid, "/NON_ZERO_CELLS", 0);
	if (gid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Properties\"\n");
		exit(1);
	}
	status = H5Gclose(gid);
	return fid;
}

// Create a group within file
hid_t GH5_open_group(hid_t fp, const char *name) {
	hid_t grpid;
	grpid = H5Gopen(fp, name);
	// this should not happen, if life is running usual
	if (grpid < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to open group %s\n", name);
		exit(1);
	}
	return grpid;
}

// General fucntion to create dataset
hid_t GH5_createdataset(hid_t gid, hid_t dataspace, const char *dsetname, unsigned type) {
	hid_t dataset, datatype;
	hsize_t *dimf;
	int i;
	// set datatype
	switch (type) {
		case 1:
			datatype = H5Tcopy(H5T_NATIVE_INT);
			break;
		case 2:
			datatype = H5Tcopy(H5T_NATIVE_FLOAT);
			break;
		case 3:
			datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
			break;
		case 4:
			datatype = H5Tcopy(H5T_NATIVE_CHAR);
			break;
		case 5:
			datatype = H5Tcopy(H5T_NATIVE_UINT);
			break;
		default:
			fprintf(stderr, "GH5 ERROR, Unkown datatype passed to dataset %s\n", dsetname);
			exit(1);
	}
	dataset = H5Dcreate(gid, dsetname, datatype, dataspace, H5P_DEFAULT);
	//Not expecting this error
	if (dataset < 0) {
		fprintf(stderr, "GH5 ERROR: Failed to create dataset %s\n", dsetname);
		exit(1);
	}
	return dataset;
}

// Write mesh
void GH5_write_mesh_data(hid_t fp, int conns, int points, int *conndata, double *ptsdata) {
	hid_t mesh_grp, conn_dataid, xyz_dataid;
	hid_t conn_spc, pts_spc;
	hsize_t conn_dim[2] = { conns, 4 };
	hsize_t pts_dim[2] = { points, 3 };
	herr_t status;

	// create group of Meshdata
	mesh_grp = GH5_open_group(fp, "/Mesh");

	//create data spaces
	conn_spc = H5Screate_simple(2, conn_dim, 0);
	pts_spc = H5Screate_simple(2, pts_dim, 0);

	//create dataset for Connections
	conn_dataid = GH5_createdataset(mesh_grp, conn_spc, "Connections", INT);

	//create dataset for points
	xyz_dataid = GH5_createdataset(mesh_grp, pts_spc, "Points", DOUBLE);

	//write connectivity data
	status = H5Dwrite(conn_dataid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conndata);

	//close conn dataset and dataspcae
	status = H5Dclose(conn_dataid);
	status = H5Sclose(conn_spc);

	//write XYZ data as points
	status = H5Dwrite(xyz_dataid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptsdata);

	//close xyz dataset and dataspace
	status = H5Dclose(xyz_dataid);
	status = H5Sclose(pts_spc);

	//close Mesh Group
	status = H5Gclose(mesh_grp);

	return;
}

void GH5_write_zero_keys(hid_t h5fid, int size_zero, vector<unsigned>& key) {

	hsize_t key_dim[2] = { size_zero, 2 };
//	int ch_size = size_zero > 2056 ? 2056 : size_zero;
//	hsize_t chunk_dims[2] = { 256, 2 };

	// create group of Meshdata
	hid_t zero_cell = GH5_open_group(h5fid, "/ZERO_CELLS");

	//create data spaces
	hid_t key_spc = H5Screate_simple(2, key_dim, 0);

//	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);

	herr_t status;

//status = H5Pset_chunk(plist_id, 2, chunk_dims);

//	status = H5Pset_deflate(plist_id, 9);

	//create dataset for Connections
//	hid_t key_dataid = H5Dcreate2(h5fid, "KEYS", H5T_NATIVE_UINT, key_spc, H5P_DEFAULT, plist_id,
//	H5P_DEFAULT);
	hid_t key_dataid = H5Dcreate2(zero_cell, "KEYS", H5T_NATIVE_UINT, key_spc, H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);
//	hid_t key_dataid = GH5_createdataset(zero_cell, key_spc, "KEYS", UNSINT);

//write connectivity data
	status = H5Dwrite(key_dataid, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &key[0]);

	//close conn dataset and dataspcae
//	status = H5Pclose(plist_id);
	status = H5Dclose(key_dataid);
	status = H5Sclose(key_spc);

	//close Mesh Group
	status = H5Gclose(zero_cell);
}

void GH5_write_non_zero_keys_sol(hid_t h5fid, int size, vector<unsigned>& key,
    vector<double>& state) {

	hsize_t key_dim[2] = { size, 2 };
	hsize_t sol_dim[2] = { size, 3 };
////	int ch_size = size_zero > 2056 ? 2056 : size_zero;
//	hsize_t key_chunk_dims[2] = { 256, 2 };
//	hsize_t sol_chunk_dims[2] = { 256, 3 };

// create group of Meshdata
	hid_t non_zero_cell = GH5_open_group(h5fid, "/NON_ZERO_CELLS");

//	status = H5Pset_chunk (key_plist_id, 2, key_chunk_dims);
//	status = H5Pset_chunk (sol_plist_id, 2, sol_chunk_dims);

//	hid_t key_dataid = GH5_createdataset(zero_cell, key_spc, "KEYS", UNSINT);
//create data spaces
	hid_t key_spc = H5Screate_simple(2, key_dim, 0);
	hid_t sol_spc = H5Screate_simple(2, sol_dim, 0);


//	hid_t key_plist_id = H5Pcreate(H5P_DATASET_CREATE);
//	hid_t sol_plist_id = H5Pcreate(H5P_DATASET_CREATE);
//	hid_t pal = H5Pcreate(H5P_DATASET_ACCESS);

//	//The number of chunk slots in the raw data chunk cache for this dataset
//	size_t nslots = 521;
//
//	//The total size of the raw data chunk cache for this dataset
//	size_t nbytes = 1 * 1024 * 1024;
//
//	//The chunk preemption policy for this dataset
//	double w0 = .75;
//
//	herr_t status;
//
////	status = H5Pset_chunk_cache(key_plist_id, nslots, nbytes, w0);
//	status = H5Pset_chunk_cache(pal, nslots, nbytes, w0);
	herr_t status;
//	status = H5Pset_deflate(key_plist_id, 9);
//	status = H5Pset_deflate(sol_plist_id, 9);

	//create dataset for Connections
//	hid_t key_dataid = H5Dcreate2(h5fid, "KEYS", H5T_NATIVE_UINT, key_spc, H5P_DEFAULT, plist_id,
//	H5P_DEFAULT);
	hid_t key_dataid = H5Dcreate2(non_zero_cell, "KEYS", H5T_NATIVE_UINT, key_spc, H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);
	hid_t sol_dataid = H5Dcreate2(non_zero_cell, "SOLUTIONS", H5T_NATIVE_DOUBLE, sol_spc, H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);

//write data
	status = H5Dwrite(key_dataid, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &key[0]);
	status = H5Dwrite(sol_dataid, H5T_NATIVE_DOUBLE_g, H5S_ALL, H5S_ALL, H5P_DEFAULT, &state[0]);

	//close dataset and dataspcae
//	status = H5Pclose(key_plist_id);
//	status = H5Pclose(sol_plist_id);
//	status = H5Pclose(pal);
	status = H5Dclose(key_dataid);
	status = H5Dclose(sol_dataid);
	status = H5Sclose(key_spc);
	status = H5Sclose(sol_spc);

	//close Mesh Group
	status = H5Gclose(non_zero_cell);
}

//write variables
void GH5_write_state_vars(hid_t fid, int num_elms, double *state_var, const char *var_name) {
	hid_t vars_grp, dataid, dataspc;
	herr_t status;
	hsize_t dims[2] = { num_elms, 1 };
	int i;
	//create group for state_vars
	vars_grp = GH5_open_group(fid, "/Properties");

	// Create dataspace
	dataspc = H5Screate_simple(2, dims, 0);
	//start dataset for var[i]
	dataid = GH5_createdataset(vars_grp, dataspc, var_name, DOUBLE);

	//write var[i] to the dataset
	status = H5Dwrite(dataid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_var);

	//close dataset and dataspace for var[i]
	status = H5Dclose(dataid);
	status = H5Sclose(dataspc);

	//close variables group
	status = H5Gclose(vars_grp);

	return;
}

///*
// * =====================================================================================
// *
// *       Filename:  hd5calls.c : If you don't know anything about HDF5
// *                               THEN LEAVE THIS FILE ALONE
// *
// *    Description:  This file contains function calls to hdf5 library
// *                  and provides simpler titan specific interface
// *
// *        Created:  05/17/2007 09:36:27 AM EDT
// *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
// *
// * This software can be redistributed free of charge.  See COPYING
// * file in the top distribution directory for more details.
// *
// * This software is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// *
// * =====================================================================================
// * $Id:$
// */
//
//#include <hdf5.h>
//#include <stdio.h>
//#include <string.h>
//
//#include "GMFG_hdfconstants.h"
//
//hid_t GH5_openfile(const char *filename, char mode) {
//	hid_t fid, gid;
//	herr_t status;
//	// if its a new file
//	if (mode == 'n') {
//		fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//		if (fid < 0) {
//			fprintf(stderr,
//			    "Unable to open new hdf file.\
//	       If %s already exists  remove it or move it to a different directory\n",
//			    filename);
//			exit(0);
//		}
//	}
//	// if the file already exists
//	else if (mode == 'o') {
//		fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
//		if (fid < 0)
//			fprintf(stderr, "Unable to open %s. Make sure its in the current directory\n", filename);
//		return fid;
//	} else {
//		fprintf(stderr, "GH5 ERROR: Unknown file access option, EXITING\n");
//		exit(1);
//	}
//	//Create base groups for data sets
//	gid = H5Gcreate(fid, "/Mesh", 0);
//	// this should not happen
//	if (gid < 0) {
//		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Mesh\"\n");
//		exit(1);
//	}
//	status = H5Gclose(gid);
//	gid = H5Gcreate(fid, "/Properties", 0);
//	if (gid < 0) {
//		fprintf(stderr, "GH5 ERROR: Failed to create base group \"/Properties\"\n");
//		exit(1);
//	}
//	status = H5Gclose(gid);
//	return fid;
//}
//
//// Create a group within file
//hid_t GH5_open_group(hid_t fp, const char *name) {
//	hid_t grpid;
//	grpid = H5Gopen(fp, name);
//	// this should not happen, if life is running usual
//	if (grpid < 0) {
//		fprintf(stderr, "GH5 ERROR: Failed to open group %s\n", name);
//		exit(1);
//	}
//	return grpid;
//}
//
//// General fucntion to create dataset
//hid_t GH5_createdataset(hid_t gid, hid_t dataspace, const char *dsetname, unsigned type) {
//	hid_t dataset, datatype;
//	hsize_t *dimf;
//	int i;
//	// set datatype
//	switch (type) {
//		case 1:
//			datatype = H5Tcopy(H5T_NATIVE_INT);
//			break;
//		case 2:
//			datatype = H5Tcopy(H5T_NATIVE_FLOAT);
//			break;
//		case 3:
//			datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
//			break;
//		case 4:
//			datatype = H5Tcopy(H5T_NATIVE_CHAR);
//			break;
//		default:
//			fprintf(stderr, "GH5 ERROR, Unkown datatype passed to dataset %s\n", dsetname);
//			exit(1);
//	}
//	dataset = H5Dcreate(gid, dsetname, datatype, dataspace, H5P_DEFAULT);
//	//Not expecting this error
//	if (dataset < 0) {
//		fprintf(stderr, "GH5 ERROR: Failed to create dataset %s\n", dsetname);
//		exit(1);
//	}
//	return dataset;
//}
//
//// Write mesh
//void GH5_write_mesh_data(hid_t fp, int conns, int points, int *conndata, double *ptsdata) {
//	hid_t mesh_grp, conn_dataid, xyz_dataid;
//	hid_t conn_spc, pts_spc;
//	hsize_t conn_dim[2] = { conns, 4 };
//	hsize_t pts_dim[2] = { points, 3 };
//	herr_t status;
//
//	// create group of Meshdata
//	mesh_grp = GH5_open_group(fp, "/Mesh");
//
//	//create data spaces
//	conn_spc = H5Screate_simple(2, conn_dim, 0);
//	pts_spc = H5Screate_simple(2, pts_dim, 0);
//
//	//create dataset for Connections
//	conn_dataid = GH5_createdataset(mesh_grp, conn_spc, "Connections", INT);
//
//	//create dataset for points
//	xyz_dataid = GH5_createdataset(mesh_grp, pts_spc, "Points", DOUBLE);
//
//	//write connectivity data
//	status = H5Dwrite(conn_dataid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conndata);
//
//	//close conn dataset and dataspcae
//	status = H5Dclose(conn_dataid);
//	status = H5Sclose(conn_spc);
//
//	//write XYZ data as points
//	status = H5Dwrite(xyz_dataid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptsdata);
//
//	//close xyz dataset and dataspace
//	status = H5Dclose(xyz_dataid);
//	status = H5Sclose(pts_spc);
//
//	//close Mesh Group
//	status = H5Gclose(mesh_grp);
//
//	return;
//}
//
////write variables
//void GH5_write_state_vars(hid_t fid, int num_elms, double *state_var, const char *var_name) {
//	hid_t vars_grp, dataid, dataspc;
//	herr_t status;
//	hsize_t dims[2] = { num_elms, 1 };
//	int i;
//	//create group for state_vars
//	vars_grp = GH5_open_group(fid, "/Properties");
//
//	// Create dataspace
//	dataspc = H5Screate_simple(2, dims, 0);
//	//start dataset for var[i]
//	dataid = GH5_createdataset(vars_grp, dataspc, var_name, DOUBLE);
//
//	//write var[i] to the dataset
//	status = H5Dwrite(dataid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_var);
//
//	//close dataset and dataspace for var[i]
//	status = H5Dclose(dataid);
//	status = H5Sclose(dataspc);
//
//	//close variables group
//	status = H5Gclose(vars_grp);
//
//	return;
//}
