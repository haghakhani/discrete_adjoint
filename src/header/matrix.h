/*
 * matrix.h
 *
 *  Created on: Sep 9, 2015
 *      Author: haghakha
 */

#ifndef SRC_HEADER_MATRIX_H_
#define SRC_HEADER_MATRIX_H_

template<typename T, unsigned Row, unsigned Col>
class Matrix {

protected:
	T* mat;
	unsigned rows;
	unsigned cols;

public:
	Matrix() {
		mat = new T[Row * Col];
		rows = Row;
		cols = Col;

		//this line initialize the matrix to zero depending upon the type of T
		for (unsigned i = 0; i < rows * cols; ++i)
			mat[i] = T();

	}
	;

	Matrix(const Matrix& rhs) {

		rows = rhs.get_rows();
		cols = rhs.get_cols();

		mat = new T[rows * cols];

		for (unsigned i = 0; i < rows * cols; ++i)
			mat[i] = rhs.mat[i];
	}
	;

	virtual ~Matrix() {
		delete[] mat;
//		std::cout<<"Now the matrix is getting deleted \n";

	}
	;

	// Operator overloading, for "standard" mathematical matrix operations
	Matrix<T, Row, Col>& operator=(const Matrix<T, Row, Col>& rhs) {

		//	std::assert(rows == rhs.get_rows() && cols == rhs.get_cols());
		if (this != &rhs)
			for (unsigned i = 0; i < rows * cols; ++i)
				mat[i] = rhs.mat[i];

		return *this;
	}
	;

	// Matrix mathematical operations
	Matrix<T, Row, Col> operator+(const Matrix<T, Row, Col>& rhs) {

		Matrix<T, Row, Col> result;

		for (unsigned i = 0; i < rows * cols; ++i)
			result.mat[i] = this->mat[i] + rhs.mat[i];

		return result;
	}
	;

	Matrix<T, Row, Col>& operator+=(const Matrix<T, Row, Col>& rhs) {

		for (unsigned i = 0; i < rows * cols; ++i)
			this->mat[i] += rhs.mat[i];
		return *this;
	}
	;

	// Matrix/scalar operations
	Matrix<T, Row, Col> operator+(const T& rhs) {

		Matrix<T, Row, Col> result;

		for (unsigned i = 0; i < rows * cols; ++i)
			result.mat[i] = this->mat[i] + rhs;

		return result;
	}
	;

	Matrix<T, Row, Col> operator-(const T& rhs) {

		Matrix<T, Row, Col> result;

		for (unsigned i = 0; i < rows * cols; ++i)
			result.mat[i] = this->mat[i] - rhs;

		return result;
	}
	;
	Matrix<T, Row, Col> operator*(const T& rhs) {

		Matrix<T, Row, Col> result;

		for (unsigned i = 0; i < rows * cols; ++i)
			result.mat[i] = this->mat[i] * rhs;

		return result;
	}
	;
	Matrix<T, Row, Col> operator/(const T& rhs) {

		Matrix<T, Row, Col> result;

		for (unsigned i = 0; i < rows * cols; ++i)
			result.mat[i] = this->mat[i] / rhs;

		return result;
	}
	;

	// Access the individual elements
	T& operator()(const unsigned& row, const unsigned& col) {

		return this->mat[row * cols + col];
	}
	;
	const T& operator()(const unsigned& row, const unsigned& col) const {

		return this->mat[row * cols + col];
	}
	;

	// Access the row and column sizes
	unsigned get_rows() const {
		return rows;
	}
	;
	unsigned get_cols() const {
		return cols;
	}
	;

};

typedef Matrix<double, 3, 3> Mat3x3;

//this matrix is a 3x3x3 matrix that keeps the sensitivity of each elements flux
//w.r.t to its two possible neighbors and itself when we are computing its flux
//in specific direction
template<unsigned size>
class Vec_Mat {

private:
	Mat3x3* vec_mat;

public:
	Vec_Mat() {

		vec_mat = new Mat3x3[size];
	}
	;

	Vec_Mat(const Vec_Mat& rhs) {

		vec_mat = new Mat3x3[size];
		for (unsigned i = 0; i < size; ++i)
			vec_mat[i] = rhs.vec_mat[i];

	}
	;

	~Vec_Mat() {

		delete[] vec_mat;
//		std::cout<<"Now the vec_mat is getting deleted \n";
	}
	;

	Mat3x3& operator()(const unsigned& index) {
		return vec_mat[index];

	}
	;

	double operator()(const unsigned& index, const unsigned& i, const unsigned& j) const{
		Mat3x3& A = vec_mat[index];
		return A(i, j);
	}
	;

	double& operator()(const unsigned& index, const unsigned& i, const unsigned& j) {
		Mat3x3& A = vec_mat[index];
		return A(i, j);
	}
	;

	Vec_Mat<size>& operator=(const Vec_Mat<size>& rhs) {
		if (this != &rhs)
			for (unsigned i = 0; i < size; ++i)
				vec_mat[i] = rhs.vec_mat[i];

		return *this;

	}

};

//this is a 2x2 matrix so has four Flux_Sens_Mat the order of indexing is:
// (side:x=0,y=1)(direction:neg=0,pos=1)
class FluxJac {

private:
	Matrix<Vec_Mat<3>, 2, 2> mat_flux_jac;

public:

	void set(int side, int dir, int index, const Mat3x3 jac_matrix) {

		(mat_flux_jac(side, dir))(index) = jac_matrix;

	}
	;
	void set(int side, int direction, int indx, Mat3x3& jac1, Mat3x3& jac2) {
		Mat3x3 dummy;
		dummy = (jac1 + jac2) * .5;
		set(side, direction, indx, dummy);

	}
	;
	Mat3x3& operator()(const unsigned& i, const unsigned& j, const unsigned& index) {
		return (mat_flux_jac(i, j))(index);
	}
	;

	virtual ~FluxJac() {
//		std::cout<<"Now flux_jac is getting deleted \n";
	}
	;

};

#endif /* SRC_HEADER_MATRIX_H_ */
