/*
 * matrix.h
 *
 *  Created on: Apr 9, 2010
 *      Author: schuette
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>

class Matrix {
public:
    /**
     * Default constructor. 
     */
	Matrix();
    /**
     * Copy constructor of class matrix.
     * @param m the matrix to be copied.
     */
	Matrix(const Matrix& m);
    /**
     * construct a matrix of dimension (rows,cols)
     * @param rows number of rows of the matrix.
     * @param cols number of columns of the matrix.
     */
	Matrix(int rows, int cols);
    /**
     * operator= for class matrix.
     * @param rhs the matrix to be assigned.
     * @return reference to ourselves.
     */
	Matrix& operator=(const Matrix& rhs);
    
    /**
     * destructs and resets matrix.
     */
	void erase();
    /**
     * calculates self = m1 * m2;
     * @param m1 operand 1
     * @param m2 operand 2
     */
    void scalarMatMatMult(double &a, const Matrix& m1, const Matrix& m2);
    void matMatMult(const Matrix &m1, const Matrix &m2);
    /**
     * helper function to copy matrix m.
     * @param m matrix to be copied.
     */
    //	void copyMat(const Matrix& m);
    /**
     * helper function to allocate the necessary memory.
     * @param rows number of rows to allocate.
     * @param cols number of columns to allocate.
     */
	void resize(int rows, int cols);
    /**
     * sets an element in the matrix.
     * @param row row of the element.
     * @param col column of the element.
     * @param value value of the element.
     */
	void set(int row, int col, double value);
    /**
     * retrieves element from the matrix.
     * @param row row of the element.
     * @param col column of the element.
     * @return value of the element.
     */
	double get(int row, int col) const;
    /**
     * copies data from matrix smat from (offset_rows,offset_cols) on.
     * @param smat source matrix
     * @param offset_rows row offset
     * @param offset_cols column offset
     * @param pref scalar multiplicator (defaults to 1).
     */
    /**
     * transposes the matrix.
     */
    void transpose();
    /**
     * returns the transposed of the matrix without harming the actual matrix
     */
    Matrix returnTransposed();
    /**
     * returns a block out of the matrix
     */
    Matrix cutBlock(const int fromRow, const int toRow, const int fromCol, const int toCol);
    /**
     * initializes the matrix with zero's.
     */
	void zero();
    /**
     * initialize the matrix with a unity matrix.
     */
    void one();
    
    /**
     * diagonalizes the matrix using lapack function dysev
     * and dstev
     */
    std::vector<double> diag();
    std::vector<double> tdiag();
    
    /**
     * multiplies the matrix by a scalar number
     * @param scalar scalar number
     */
	void multiply(double scalar);
    /**
     * print the matrix to std::cout.
     */
	void print();
    
    /**
     * Destructor for the Matrix class. Releases all reserved memory.
     */
	virtual ~Matrix();
    
    /**
     * returns number of matrix columns.
     * @return number of columns
     */
    int getCols() const;
    /**
     * return number of matrix rows.
     * @return number of rows
     */
    int getRows() const;
    
    /**
     * check whether matrix is square;
     * @return 
     */
    bool isSquare() const;
    
    /**
     * number of elements in the matrix
     * @return 
     */
	int getSize();
    /**
     * standard matrix operations
     */
    Matrix &operator+=(const Matrix &rhs);
    Matrix &operator-=(const Matrix &rhs);
    Matrix &operator*=(const Matrix &rhs);
    
    const Matrix operator+(const Matrix &rhs) const;
    const Matrix operator-(const Matrix &rhs) const;
    const Matrix operator*(const Matrix &rhs) const;
    void d() const;
private:
    //! number of rows in the matrix.
    int m_rows;
    //! number of columns in the matrix.
	int m_cols;
    //! pointer to the data memory.
	double* m_data;
    /**
     * retrieves pointer to data memory.
     * @return pointer 
     */
    double* getData() const;
};

#endif /* MATRIX_H_ */
