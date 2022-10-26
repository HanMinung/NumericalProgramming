/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 10-24-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"
#include "myNP.h"

/*---------------------------------		BASIC TOOL		---------------------------------------*/

// Create Matrix with specified size
Matrix createMat(int _rows, int _cols)
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
void freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix txt2Mat(std::string _filePath, std::string _fileName)
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
void printMat(Matrix _A, const char* _name)
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
void initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {

			_A.at[i][j] = _val;
		}

	}

}

// Create matrix of all zeros
Matrix zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	
	initMat(Out,0);

	return Out;
}

Matrix  ones(int _rows, int _cols) {

	Matrix Out = createMat(_rows, _cols);

	initMat(Out, 1);

	return Out;

}


Matrix transpose(Matrix _A) {
	
	Matrix Out = zeros(_A.cols, _A.rows);

	FOR_LOOP(i, 0, _A.rows, 1) {
		FOR_LOOP(j, 0, _A.cols, 1) {

			Out.at[j][i] = _A.at[i][j];

		}
	}
	return Out;
}


Matrix eye(int _rows, int _cols) {

	Matrix Out = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++) {

		for (int j = 0; j < _cols; j++) {


			Out.at[i][j] = 0;

			if (i == j) {
				Out.at[i][j] = 1;
			}
		}
	}
	return Out;
}


// Copy matrix A to B ( A --> B )
void copyMat(Matrix _A, Matrix _B) {

	if (_A.cols != _B.cols || _A.rows != _B.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error - 'copyMat' function");
		printf("\n*************************************************\n");
		/*return createMat(0, 0);*/

		exit(1);
	}

	FOR_LOOP(i, 0, _A.rows, 1) {
		FOR_LOOP(j, 0, _A.cols, 1) {

			_B.at[i][j] = _A.at[i][j];
		}
	}	

}

Matrix vectExt(Matrix _A, int in) {

	Matrix Out = createMat(_A.rows, 1);

	FOR_LOOP(i,0,_A.rows,1){

		Out.at[i][0] = _A.at[i][in];
	}

	return Out;
}


void vectins(Matrix _v, int in, Matrix _A) {

	if (_v.cols != 1) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'vectins' function");
		printf("\n*************************************************\n");
		exit(1);
	}

	if (_v.rows != _A.rows)
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'vectins' function");
		printf("\n*************************************************\n");
		exit(1);
	}

	for (int i = 0; i < _v.rows; i++)
	{
		_A.at[i][in] = _v.at[i][0];
	}

	//printMat(_A, "_A in vectins");
}


Matrix multiScalar(double scalar, Matrix _A) {

	Matrix Out = zeros(_A.rows, _A.cols);

	copyMat(_A, Out);

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {

			Out.at[i][j] = scalar * _A.at[i][j];

		}
	}

	return Out;
}


/*----------------------------------------------------------------------------------*/

// NxN matrix
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d) {

	if (_A.cols != _b.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR ! Dimension Error : starting gaussElim ");
		printf("\n*************************************************\n");

		exit(1);
	}

	double mult = 0;
	uint8_t row = _A.rows;	// n(rows)

	FOR_LOOP(k, 0, row - 1, 1) {

		FOR_LOOP(i, k + 1, row, 1) {

			mult = _A.at[i][k] / _A.at[k][k];
			_A.at[i][k] = 0;

			FOR_LOOP(j, k+1, row, 1) {

				_A.at[i][j] = _A.at[i][j] - (mult * _A.at[k][j]);
				
			}
			
			// input에 대한 항은 따로 계산
			_b.at[i][0] = _b.at[i][0] - mult * _b.at[k][0];
			
		}

	}

	copyMat(_A, _U);	// A --> U
	copyMat(_b, _d);	// b --> d

}

// NxN 행렬
// Main statement : _L = eye(n,n); 
//					_U = zeros(n, n);
void LUdecomp(Matrix _A, Matrix _L, Matrix _U) {

	if ((_A.cols != _L.cols) || (_A.rows != _L.rows) || (_A.rows != _U.cols) || (_A.cols != _U.rows)) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at : 'LUdecomp' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	uint8_t row = _A.rows;
	double mult = 0;

	copyMat(_A, _U);

	FOR_LOOP(k, 0, row - 1, 1) {

		FOR_LOOP(i, k + 1, row, 1) {

			mult = _U.at[i][k] / _U.at[k][k];
			
			_L.at[i][k] = mult;
			_U.at[i][k] = 0;

			FOR_LOOP(j, k+1, row, 1) {

				_U.at[i][j] = _U.at[i][j] - mult * _U.at[k][j];
			}
		}
		
	}

}


