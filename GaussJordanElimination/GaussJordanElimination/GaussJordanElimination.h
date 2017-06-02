#pragma once
#include <iostream>
#include <cstdlib>
#include <array>

class GaussJordanElimination {
private:
	const static int MAX_MATRIX_SIZE = 10;
	typedef std::array<std::array<int, MAX_MATRIX_SIZE>, MAX_MATRIX_SIZE> matrix;
	typedef std::array<std::array<double, MAX_MATRIX_SIZE>, MAX_MATRIX_SIZE> iMatrix;
	typedef std::array<int, MAX_MATRIX_SIZE> vector;
	typedef std::array<double, MAX_MATRIX_SIZE> iVector;
	int matrixSize;
	matrix originalMatrix;
	vector resultValue;
	void gauss(matrix &tempMatrix, vector &tempValue);
	void gauss(matrix &tempMatrix, matrix &inverseMatrix);
public:
	GaussJordanElimination();
	void printInverseMatrix();
	void printAnswerValue();
};