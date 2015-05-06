/*
 * Matrix.cpp
 *
 *  Created on: Apr 9, 2010
 *      Author: schuette
 */

#include "matrix.h"
#include "diag.h"
#include <memory.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;

Matrix::Matrix()
{
	m_data = NULL;
	m_rows = 0;
	m_cols = 0;
}

void Matrix::erase()
{
	if(m_data != NULL)
		delete[] m_data;
	m_rows = 0;
	m_cols = 0;
	m_data = NULL;
}

void Matrix::resize(int rows, int cols)
{
    erase();
    if ((rows > 0) && (cols > 0))
    {
        m_data = new double[rows*cols];
        m_rows = rows;
        m_cols = cols;
    }
}

int Matrix::getCols() const
{
	return m_cols;
}

int Matrix::getRows() const
{
	return m_rows;
}

Matrix& Matrix::operator=(const Matrix& rhs)
{
    if (this != &rhs)
    {
        if ((rhs.getCols() > 0) && (rhs.getRows() > 0))
        {
            if((m_rows != rhs.getRows()) || (m_cols != rhs.getCols()))
                resize(rhs.getRows(), rhs.getCols());
            for (int i = 0; i < m_rows; i++)
                for (int j = 0; j < m_cols; j++)
                    set(i, j, rhs.get(i, j));
        }
        else
            erase();
    }
    
	return *this;
}

Matrix::Matrix(const Matrix& m)
{
    m_data = NULL;
    m_rows = 0;
    m_cols = 0;
    resize(m.getRows(), m.getCols());
    for (int i = 0; i < m_rows; i++)
        for (int j = 0; j < m_cols; j++)
            set(i, j, m.get(i, j));
}

bool Matrix::isSquare() const {
    if(getRows() == 0)
        return false;
    if(getRows() != getCols())
        return false;
    return true;
}

Matrix::Matrix(int rows, int cols)
{
    if ((rows > 0) && (cols > 0))
    {
        m_rows = rows;
        m_cols = cols;
        m_data = new double[rows*cols];
    }
    else {
        m_data = NULL;
        m_rows = 0;
        m_cols = 0;
    }
}

void Matrix::transpose()
{
	Matrix result(m_cols, m_rows);
	for(int i = 0; i < m_rows; i++)
	{
		for(int j = 0; j < m_cols; j++)
		{
			result.set(j, i, get(i, j));
		}
	}
	*this = result;
}

Matrix Matrix::returnTransposed()
{
    Matrix result = *this;
    result.transpose();
    
    return result;
}

Matrix Matrix::cutBlock(const int fromRow, const int toRow, const int fromCol, const int toCol)
{
	if ((fromRow >= 0) && (toRow < m_rows) && (fromCol >= 0) && (toCol < m_cols))
	{
		int rows = toRow-fromRow, cols = toCol-fromCol;
		Matrix result(rows+1, cols+1);//std::cout << result.getRows() << std::endl;
		for (int i = 0; i <= rows; i++)
			for (int j = 0; j <= cols; j++)
				result.set(i, j, m_data[fromRow+i + (fromCol+j)*m_cols]);
		return result;
	} else {
		std::cout << "in cutBlock: matrix dimension does not fit" << std::endl;
	}

	return *this;
}

void Matrix::scalarMatMatMult(double &a, const Matrix& m1, const Matrix& m2)
{
    if(m2.getRows() != m1.getCols())
    {
        std::cerr << "Matrix dimensions incompatible" << std::endl;
    }
    if((m_rows != m1.getRows()) || (m_cols != m2.getCols()) )
        resize(m1.getRows(), m2.getCols());
    int w = m1.getRows();
    int wp = m2.getCols();
    int k_max = m1.getCols();
    double scalar = a;
    double zero = 0.0;
    char N = 'N';
    utils::dgemm_(&N, &N, &w, &wp, &k_max, &scalar, m1.m_data, &w, m2.m_data, &k_max, &zero, m_data, &w);
}

void Matrix::matMatMult(const Matrix &m1, const Matrix &m2)
{
	double one = 1;
    scalarMatMatMult(one, m1, m2);
}

void Matrix::print()
{
	for(int i = 0; i < m_rows; i++)
	{
		for(int j = 0; j < m_cols; j++)
            std::cout << std::scientific << std::setprecision(6) << get(i,j) << "\t";
		std::cout << std::endl;
	}
}

void Matrix::multiply(double scalar)
{
    if (m_data != NULL)
    {
        int end = m_rows*m_cols;
        for(int i = 0; i < end; i++)
            m_data[i] *= scalar;
    }
    
}

void Matrix::zero()
{
    int end = m_rows*m_cols;
	for(int i = 0; i < end; i++)
		m_data[i] = 0;
}

void Matrix::one()
{
    if(m_rows != m_cols)
    {
        std::cout << "not a square matrix!" << std::endl;
        return;
    }
    zero();
    for(int i = 0; i < m_rows; i++)
        m_data[i + m_cols * i] = 1.0;
}

void Matrix::set(int row, int col, double value)
{
    if ((m_rows >= row) && (m_cols >= col))
        m_data[(row) + (col) * m_rows] = value;
    else
        std::cout << "Element (" << row << ", " << col << ") does not exist" << std::endl;
}

double Matrix::get(int row, int col) const
{
    if ((m_rows >= row) && (m_cols >= col))
        return m_data[(row) + (col) * m_rows];
    else
        std::cout << "Element (" << row << ", " << col << ") does not exist" << std::endl;
    
    return 0;
}

std::vector<double> Matrix::diag()
{
    int N = getRows();
	std::vector<double> eigenvalues(N);
    char v = 'V', u = 'U';
	int info, lwork = max(1, 3*N-1);
    double *work = NULL;
    if(!isSquare())
        return eigenvalues;
    
    work = new double[lwork];
	utils::dsyev_(&v, &u, &N, getData(), &N, &eigenvalues.at(0), work, &lwork, &info);
    
	if (info != 0)
		cout << "diagonalization routine dsyev_ error: info = " << info << endl;
    
//   (*this).transpose();

	delete[] work;
	return eigenvalues;
}


Matrix::~Matrix() {
	erase();
    
}

int Matrix::getSize()
{
    return m_rows * m_cols;
}

Matrix &Matrix::operator+=(const Matrix &rhs)
{
    if ((m_data == NULL) || (rhs.m_data == NULL))
        std::cout << "Error in +=: Matrix not initialized" << std::endl;
    if (this != &rhs)
    {
        if ((m_cols == rhs.getCols()) && (m_rows == rhs.getRows()))
        {
            for (int i = 0; i < m_rows; i++)
                for (int j = 0; j < m_cols; j++)
                    set(i, j, get(i, j) + rhs.get(i, j));
        }
    } else {
        Matrix help = *this;
        *this = help + help;
    }
    
    return *this;  
}

Matrix &Matrix::operator-=(const Matrix &rhs)
{
    if ((m_data == NULL) || (rhs.m_data == NULL))
        std::cout << "Error in -=: Matrix not initialized" << std::endl;
    if (this != &rhs)
    {
        if ((m_cols == rhs.getCols()) && (m_rows == rhs.getRows()))
        {
            for (int i = 0; i < m_rows; i++)
                for (int j = 0; j < m_cols; j++)
                    set(i, j, get(i, j) - rhs.get(i, j));
        }
    } else {
        zero();
    }
    
    return *this;  
}

Matrix &Matrix::operator*=(const Matrix &rhs)
{
    
    if ((m_data == NULL) || (rhs.m_data == NULL))
        std::cout << "Error in *=: Matrix not initialized" << std::endl;
    Matrix m = *this;
    matMatMult(m, rhs);
    
    return *this;  
}


const Matrix Matrix::operator+(const Matrix &rhs) const
{
    return Matrix(*this) += rhs;
}

const Matrix Matrix::operator-(const Matrix &rhs) const
{
    return Matrix(*this) -= rhs;
}

const Matrix Matrix::operator*(const Matrix &rhs) const
{
    return Matrix(*this) *= rhs;
}

double* Matrix::getData() const
{
	return m_data;
}

void Matrix::d() const
{
	cout << m_rows << " " << m_cols << endl;
}

