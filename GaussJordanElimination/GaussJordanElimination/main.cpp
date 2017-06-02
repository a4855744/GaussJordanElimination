#include "GaussJordanElimination.h"
#include <Windows.h>

int main() {
	GaussJordanElimination &system1 = GaussJordanElimination();

	system1.printInverseMatrix();
	system1.printAnswerValue();

	system("pause");

	return 0;
}