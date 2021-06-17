/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"


// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}		

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	return Out;
}

// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.6f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}


// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{
	int i = 0;
	int k = 0;

	for (i=0; i < _A.rows; i++)
		for (k = 0; k < _A.cols; k++)
			_A.at[i][k] = _val;
}

// Create matrix of all zeros
Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 0);

	return Out;
}

void copyMat(Matrix _A, Matrix Out)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j];
}

Matrix eye(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols;j++)
		{
			if (i == j)
				Out.at[i][j] = 1;
			else 
				Out.at[i][j] = 0;
		}
	return Out;
}

void Rexchange(Matrix _A, int row1, int row2)
{
	Matrix OUT = createMat(1, _A.cols);
	initMat(OUT, 0);

	for (int j = 0; j < _A.cols; j++)
	{
		OUT.at[0][j]   = _A.at[row1][j];
		_A.at[row1][j] = _A.at[row2][j];
		_A.at[row2][j] = OUT.at[0][j];
	}
 }

Matrix Multi(Matrix _A, Matrix _B)
{
	Matrix _OUT = zeros(_A.rows, _B.cols);

	for (int j = 0; j < _B.cols; j++)
		for (int i = 0; i < _A.rows; i++)
			for (int k = 0; k < _B.rows; k++)
				_OUT.at[i][j] += _A.at[i][k] * _B.at[k][j];
	return _OUT;
}

Matrix transpose(Matrix _A)
{
	Matrix out = createMat(_A.cols, _A.rows);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			out.at[j][i] = _A.at[i][j];
	return out;
}
double norm(Matrix _A)
{
	double out = 0;

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			out += _A.at[i][j] * _A.at[i][j];

	return sqrt(out);
}

Matrix Minus(Matrix _A)
{
	Matrix out = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			out.at[i][j] = -1 * _A.at[i][j];
	return out;
}


Matrix transpose1(Matrix _A)
{
	Matrix out = createMat(_A.cols, _A.rows);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			out.at[j][i] = _A.at[i][j];

	return out;
}