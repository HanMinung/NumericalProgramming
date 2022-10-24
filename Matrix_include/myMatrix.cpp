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


double ele(Matrix _A, uint8_t i, uint8_t j) {

	return _A.at[i][j];
}


/*----------------------------------------------------------------------------------*/

// NxN matrix
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d) {

	if (_A.cols != _b.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR ! Dimension Error (starting gaussElim) ");
		printf("\n*************************************************\n");

		exit(1);
	}

	double mult = 0;
	uint8_t m = _A.rows;	// n(rows)

	FOR_LOOP(k, 0, m - 1, 1) {

		FOR_LOOP(i, k + 1, m, 1) {

			mult = ele(_A, i, k) / ele(_A, k, k);
			_A.at[i][k] = 0;

			FOR_LOOP(j, k+1, m, 1) {

				_A.at[i][j] = _A.at[i][j] - (mult * _A.at[k][j]);
				
			}
			
			// input에 대한 항은 따로 계산
			_b.at[i][0] = _b.at[i][0] - mult * _b.at[k][0];
			
		}

	}

	copyMat(_A, _U);	// A --> U
	copyMat(_b, _d);	// b --> d

}



