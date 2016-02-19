/*
 * matrix.C
 *
 *  Created on: Sep 9, 2015
 *      Author: haghakha
 */

#ifndef __MATRIX
#define __MATRIX

#include "../header/matrix.h"
#include <cassert>

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col>::Matrix() {
	mat = new T[Row * Col];
	rows = Row;
	cols = Col;
}

template<typename T, unsigned Row, unsigned Col>
Matrix<T, Row, Col>::~Matrix() {
	delete[] mat;
}

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col>& Matrix<T, Row, Col>::operator=(const Matrix<T, Row, Col>& rhs) {

//	std::assert(rows == rhs.get_rows() && cols == rhs.get_cols());

	for (unsigned i = 0; i < rows * cols; ++i)
		mat[i] = rhs.mat[i];

	return *this;
}

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col> Matrix<T, Row, Col>::operator+(const Matrix<T, Row, Col>& rhs) {

	Matrix<T, Row, Col> result;

	for (unsigned i = 0; i < rows * cols; ++i)
		result.mat[i] = this->mat[i] + rhs.mat[i];

	return result;
}

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col>& Matrix<T, Row, Col>::operator+=(const Matrix<T, Row, Col>& rhs) {

	for (unsigned i = 0; i < rows * cols; ++i)
		this->mat[i] += rhs.mat[i];
	return *this;
}

// Matrix/scalar operations
template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col> Matrix<T, Row, Col>::operator+(const T& rhs) {

	Matrix<T, Row, Col> result;

	for (unsigned i = 0; i < rows * cols; ++i)
		result.mat[i] = this->mat[i] + rhs;

	return result;
}

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col> Matrix<T, Row, Col>::operator-(const T& rhs) {

	Matrix<T, Row, Col> result;

	for (unsigned i = 0; i < rows * cols; ++i)
		result.mat[i] = this->mat[i] - rhs;

	return result;
}

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col> Matrix<T, Row, Col>::operator*(const T& rhs) {

	Matrix<T, Row, Col> result;

	for (unsigned i = 0; i < rows * cols; ++i)
		result.mat[i] = this->mat[i] * rhs;

	return result;
}

template<typename T, unsigned Row, unsigned Col>
inline Matrix<T, Row, Col> Matrix<T, Row, Col>::operator/(const T& rhs) {

	Matrix<T, Row, Col> result;

	for (unsigned i = 0; i < rows * cols; ++i)
		result.mat[i] = this->mat[i] / rhs;

	return result;
}

// Access the individual elements
template<typename T, unsigned Row, unsigned Col>
inline T& Matrix<T, Row, Col>::operator()(const unsigned& row, const unsigned& col) {

	return this->mat[row * cols + col];
}

template<typename T, unsigned Row, unsigned Col>
inline const T& Matrix<T, Row, Col>::operator()(const unsigned& row, const unsigned& col) const {

	return this->mat[row * cols + col];
}

template<unsigned size>
inline Vec_Mat<size>::Vec_Mat() {

	vec_mat = new Mat3x3[size];
}

template<unsigned size>
inline Vec_Mat<size>::~Vec_Mat() {

	delete[] vec_mat;
}

template<unsigned size>
inline Mat3x3& Vec_Mat<size>::operator()(const unsigned& index) {
	return vec_mat[index];

}

template<unsigned size>
inline double Vec_Mat<size>::operator()(const unsigned& index, const unsigned& i,
    const unsigned& j) {
	Mat3x3& A = vec_mat[index];
	return A(i, j);
}

inline void FluxJac::set(int side, int dir, int index, const Mat3x3 jac_matrix) {

	(mat_flux_jac(side, dir))(index) = jac_matrix;

}

inline void FluxJac::set(int side, int direction, int indx, Mat3x3& jac1, Mat3x3& jac2) {
	Mat3x3 dummy;
	dummy = (jac1 + jac2) * .5;
	set(side, direction, indx, dummy);

}

inline Mat3x3& FluxJac::operator()(const unsigned& i, const unsigned& j, const unsigned& index) {
	return (mat_flux_jac(i, j))(index);
}

#endif
