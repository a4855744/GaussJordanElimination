#include "GaussJordanElimination.h"
using namespace std;

int gcd(int num1, int num2) {
	if (num2 == 0)
		return num1;
	gcd(num2, num1 % num2);
}

GaussJordanElimination::GaussJordanElimination() {
	cout << "input square matrix size : ";
	cin >> matrixSize;
	if (matrixSize > MAX_MATRIX_SIZE) {
		cout << "Error!" << endl;
		return;
	}
	cout << "input matrix" << endl;
	for (int col = 0; col < matrixSize; col++)
		for (int row = 0; row < matrixSize; row++)
			cin >> originalMatrix[col][row];
	cout << "input linear system result value" << endl;
	for (int col = 0; col < matrixSize; col++)
		cin >> resultValue[col];
}

void GaussJordanElimination::printInverseMatrix() {
	matrix tempMatrix;
	matrix unitMatrix;

	for (int col = 0; col < matrixSize; col++)
		for (int row = 0; row < matrixSize; row++)
			tempMatrix[col][row] = originalMatrix[col][row];
	for (int col = 0; col < matrixSize; col++) // initialize unit matrix
		for (int row = 0; row < matrixSize; row++)
			if (col == row)
				unitMatrix[col][row] = 1;
			else
				unitMatrix[col][row] = 0;

	// check no inverse matrix

	gauss(tempMatrix, unitMatrix);

	for (int pivotPoint = matrixSize - 1; pivotPoint >= 0; pivotPoint--) { // make dialog matrix
		for (int col = 0; col < pivotPoint; col++) {
			if (tempMatrix[pivotPoint][pivotPoint] == 0) {
				cout << "no inverse matrix" << endl;
				return;
			}
			if (tempMatrix[pivotPoint][pivotPoint] != 0 && tempMatrix[col][pivotPoint] != 0) {
				int lcm = tempMatrix[pivotPoint][pivotPoint] * tempMatrix[col][pivotPoint] / gcd(tempMatrix[pivotPoint][pivotPoint], tempMatrix[col][pivotPoint]);
				auto pivotInverseVector = unitMatrix[pivotPoint];
				int mulValue1 = lcm / tempMatrix[pivotPoint][pivotPoint];
				int mulValue2 = lcm / tempMatrix[col][pivotPoint];

				tempMatrix[col][pivotPoint] = 0;
				for (int row = 0; row < matrixSize; row++) {
					tempMatrix[col][row] = tempMatrix[col][row] * mulValue2;
					pivotInverseVector[row] = pivotInverseVector[row] * mulValue1;
					unitMatrix[col][row] = unitMatrix[col][row] * mulValue2;
					unitMatrix[col][row] = unitMatrix[col][row] - pivotInverseVector[row];
				}
			}
		}
	}

	iMatrix inverseMatrix;

	for (int col = 0; col < matrixSize; col++)
		for (int row = 0; row < matrixSize; row++)
			inverseMatrix[col][row] = unitMatrix[col][row];

	for (int col = 0; col < matrixSize; col++) {
		for (int row = 0; row < matrixSize; row++)
			cout << inverseMatrix[col][row] / tempMatrix[col][col] << "\t";
		cout << endl;
	}
}

void GaussJordanElimination::printAnswerValue() {
	matrix tempMatrix = originalMatrix;
	vector tempValue = resultValue, answerValue;

	gauss(tempMatrix, tempValue);

	for (int col = matrixSize - 1; col >= 0; col--) { // backsub
		if (tempMatrix[col][col] == 0) {
			if (tempValue[col] == 0) {
				cout << "many answer" << endl;
				return;
			}
			else {
				cout << "no answer" << endl;
				return;
			}
		}
		answerValue[col] = tempValue[col] / tempMatrix[col][col];
		for (int loopValue = 0; loopValue < col; loopValue++) {
			tempValue[loopValue] = tempValue[loopValue] - answerValue[col] * tempMatrix[loopValue][col];
			tempMatrix[loopValue][col] = 0;
		}
	}

	for (int col = 0; col < matrixSize; col++)
		cout << "x" << col << " = " << answerValue[col] << endl;
}

void GaussJordanElimination::gauss(matrix &tempMatrix, vector &tempValue) {
	for (int pivotPoint = 0; pivotPoint < matrixSize; pivotPoint++) { // gauss
		int pivot = pivotPoint;

		for (int col = pivotPoint; col < matrixSize; col++)
			if (tempMatrix[col][pivotPoint] < 0) {
				for (int row = 0; row < matrixSize; row++)
					tempMatrix[col][row] = -tempMatrix[col][row];
				tempValue[col] = -tempValue[col];
			}

		for (int col = pivotPoint + 1; col < matrixSize; col++) // search pivot
			if ((tempMatrix[pivot][pivotPoint] == 0) || (tempMatrix[col][pivotPoint] != 0 && (tempMatrix[pivot][pivotPoint] > tempMatrix[col][pivotPoint])))
				pivot = col;

		auto tempVector = tempMatrix[pivotPoint]; // swap matrix[pivotPoint] and matrix[pivot]
		tempMatrix[pivotPoint] = tempMatrix[pivot];
		tempMatrix[pivot] = tempVector;

		auto temp = tempValue[pivotPoint]; // swap result[pivotPoint] and result[pivot]
		tempValue[pivotPoint] = tempValue[pivot];
		tempValue[pivot] = temp;

		for (int col = pivotPoint + 1; col < matrixSize; col++) {
			if (tempMatrix[pivotPoint][pivotPoint] != 0 && tempMatrix[col][pivotPoint] != 0) {
				int lcm = tempMatrix[pivotPoint][pivotPoint] * tempMatrix[col][pivotPoint] / gcd(tempMatrix[pivotPoint][pivotPoint], tempMatrix[col][pivotPoint]);
				vector pivotVector = tempMatrix[pivotPoint];
				int pivotResultValue = tempValue[pivotPoint];
				int mulValue1 = lcm / pivotVector[pivotPoint]; // for pivot vector
				int mulValue2 = lcm / tempMatrix[col][pivotPoint]; // for other vector

				for (int row = pivotPoint; row < matrixSize; row++) {
					pivotVector[row] = pivotVector[row] * mulValue1;
					tempMatrix[col][row] = tempMatrix[col][row] * mulValue2;
					tempMatrix[col][row] = tempMatrix[col][row] - pivotVector[row];
				}
				pivotResultValue = pivotResultValue * mulValue1;
				tempValue[col] = tempValue[col] * mulValue2;
				tempValue[col] = tempValue[col] - pivotResultValue;
			}
		}
	}
}

void GaussJordanElimination::gauss(matrix &tempMatrix, matrix &inverseMatrix) {
	for (int pivotPoint = 0; pivotPoint < matrixSize; pivotPoint++) { // gauss
		int pivot = pivotPoint;

		for (int col = pivotPoint; col < matrixSize; col++)
			if (tempMatrix[col][pivotPoint] < 0) {
				for (int row = 0; row < matrixSize; row++) {
					tempMatrix[col][row] = -tempMatrix[col][row];
					inverseMatrix[col][row] = -inverseMatrix[col][row];
				}
			}

		for (int col = pivotPoint + 1; col < matrixSize; col++) // search pivot
			if ((tempMatrix[pivot][pivotPoint] == 0) || (tempMatrix[col][pivotPoint] != 0 && (tempMatrix[pivot][pivotPoint] > tempMatrix[col][pivotPoint])))
				pivot = col;

		auto tempVector = tempMatrix[pivotPoint]; // swap matrix[pivotPoint] and matrix[pivot]
		tempMatrix[pivotPoint] = tempMatrix[pivot];
		tempMatrix[pivot] = tempVector;

		tempVector = inverseMatrix[pivotPoint]; // swap inverseMatrix[pivotPoint] and inverseMatrix[pivot]
		inverseMatrix[pivotPoint] = inverseMatrix[pivot];
		inverseMatrix[pivot] = tempVector;

		for (int col = pivotPoint + 1; col < matrixSize; col++) { // make triangular matrix
			if (tempMatrix[pivotPoint][pivotPoint] != 0 && tempMatrix[col][pivotPoint] != 0) {
				int lcm = tempMatrix[pivotPoint][pivotPoint] * tempMatrix[col][pivotPoint] / gcd(tempMatrix[pivotPoint][pivotPoint], tempMatrix[col][pivotPoint]);
				vector pivotVector = tempMatrix[pivotPoint];
				vector inverseVector = inverseMatrix[pivotPoint];
				int mulValue1 = lcm / pivotVector[pivotPoint]; // for pivot vector
				int mulValue2 = lcm / tempMatrix[col][pivotPoint]; // for other vector

				for (int row = 0; row < matrixSize; row++) {
					pivotVector[row] = pivotVector[row] * mulValue1;
					tempMatrix[col][row] = tempMatrix[col][row] * mulValue2;
					tempMatrix[col][row] = tempMatrix[col][row] - pivotVector[row];
					inverseVector[row] = inverseVector[row] * mulValue1;
					inverseMatrix[col][row] = inverseMatrix[col][row] * mulValue2;
					inverseMatrix[col][row] = inverseMatrix[col][row] - inverseVector[row];
				}
			}
		}
	}
}