/*
 * matrix.h
 *
 *  Created on: Sep 9, 2015
 *      Author: haghakha
 */

#ifndef SRC_HEADER_MATRIX_H_
#define SRC_HEADER_MATRIX_H_

#include <vector>

template<typename T> class Matrix {

protected:
	std::vector<std::vector<T> > mat;
	unsigned rows;
	unsigned cols;

public:
	Matrix();
	Matrix(unsigned _rows, unsigned _cols, const T& _initial);
	Matrix(unsigned _rows, unsigned _cols);
	Matrix(const Matrix<T>& rhs);
	virtual ~Matrix();

	// Operator overloading, for "standard" mathematical matrix operations
	Matrix<T>& operator=(const Matrix<T>& rhs);

	// Matrix mathematical operations
	Matrix<T> operator+(const Matrix<T>& rhs);
	Matrix<T>& operator+=(const Matrix<T>& rhs);
	Matrix<T> operator-(const Matrix<T>& rhs);
	Matrix<T>& operator-=(const Matrix<T>& rhs);
	Matrix<T> operator*(const Matrix<T>& rhs);
	Matrix<T>& operator*=(const Matrix<T>& rhs);
	Matrix<T> transpose();

	// Matrix/scalar operations
	Matrix<T> operator+(const T& rhs);
	Matrix<T> operator-(const T& rhs);
	Matrix<T> operator*(const T& rhs);
	Matrix<T> operator/(const T& rhs);

	// Matrix/vector operations
	std::vector<T> operator*(const std::vector<T>& rhs);
	std::vector<T> diag_vec();

	// Access the individual elements
	T& operator()(const unsigned& row, const unsigned& col);
	const T& operator()(const unsigned& row, const unsigned& col) const;

	// Access the row and column sizes
	unsigned get_rows() const;
	unsigned get_cols() const;

};

template<typename T>
class Mat3x3: public Matrix<T> {

public:
	Mat3x3() :
			Matrix<T>(3, 3, 0.) {
	}
	;

	using Matrix<T>::operator=;

	~Mat3x3() {
	}
	;
};

//this matrix is a 3x3x3 matrix that keeps the sensitivity of each elements flux
//w.r.t to its two possible neighbors and itself when we are computing its flux
//in specific direction
class Vec_Mat {

private:
	std::vector<Mat3x3<double> > vec_mat;

public:
	Vec_Mat() {
	}
	;

	Vec_Mat(int num_vec) :
			vec_mat(num_vec) {
	}
	;

	void set_size(int num_vec) {
		vec_mat.resize(num_vec);
	}
	;

	Mat3x3<double>& operator()(const unsigned& index);

	double operator()(const unsigned& index,const unsigned& i,const unsigned& j);

	~Vec_Mat() {
	}
	;

};

//this is a 2x2 matrix so has four Flux_Sens_Mat the order of indexing is:
// (side:x=0,y=1)(direction:neg=0,pos=1)
class FluxJac {

private:
	Matrix<Vec_Mat> mat_flux_jac;

public:
	FluxJac();

	void set(int side, int dir, int index, const Mat3x3<double> jac_matrix);
	void set(int side, int direction, int indx, Mat3x3<double>& jac1, Mat3x3<double>& jac2);
	Mat3x3<double>& operator()(const unsigned& index,const unsigned& i,const unsigned& j);

	~FluxJac() {
	}
	;

};

#endif /* SRC_HEADER_MATRIX_H_ */
