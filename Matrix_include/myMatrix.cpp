/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 11-02-2022
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

// Expect to extract third row : in = 3
Matrix rowExt(Matrix _A, int in) {

	Matrix Out = createMat(1, _A.cols);

	FOR_LOOP(i, 0, _A.cols, 1) {

		Out.at[0][i] = _A.at[in][i];
	}

	return Out;
}

// 실제로 어떤 행을 인자 _v로 받아서 _A의 특정 행에 넣어주는 함수 : 실제 matrix의 행값을 넣어주면 된다.
void rowIns(Matrix _v, int in, Matrix _A) {

	if (_v.rows != 1 || _v.cols != _A.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'rowIns' function");
		printf("\n*************************************************\n");
		exit(1);
	}

	FOR_LOOP(i, 0, _A.cols, 1) {

		_A.at[in][i] = _v.at[0][i];
	}
}

// exchange seq1 th row and seq2 th row
void rowExchange(Matrix _A, int seq1, int seq2) {

	Matrix row_f = rowExt(_A, seq1);
	Matrix row_s = rowExt(_A, seq2);

	rowIns(row_f, seq2, _A);
	rowIns(row_s, seq1, _A);

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

	FOR_LOOP(i,0,_v.rows,1){
	
		_A.at[i][in] = _v.at[i][0];
	}

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


Matrix row_max(Matrix _A, uint8_t k) {

	Matrix Out = zeros(_A.rows-k,1);
	double max = 0;
	double ele = 0;

	FOR_LOOP(i, k, _A.rows, 1) {

		max = _A.at[i][0];

		FOR_LOOP(j, 0, _A.cols, 1) {
			
			ele = fabs(_A.at[i][j]);
			if (ele > max)	max = ele;
		}

		Out.at[i-k][0] = max;
	}

	printf("================\n");
	printf("|  At k = %d    |\n", k);
	printf("================\n\n");
	printMat(Out, "Max elements of each row");

	return Out;

	freeMat(Out);
}


int row_SP(Matrix _A, Matrix _rowMax, uint8_t k) {

	// 각 행의 최대 SP를 저장하기 위한 매트릭스
	Matrix SP = zeros(_A.rows-k, 1);
	double max = 0;
	double SP_val = 0;
	int max_idx = 0;

	FOR_LOOP(i,k,_A.rows,1) {
		
		SP_val = fabs(_A.at[i][k] / _rowMax.at[i-k][0]);
		SP.at[i-k][0] = SP_val;

	}

	printMat(SP, "SPval");
	
	SP_val = SP.at[0][0];

	FOR_LOOP(m, 0, SP.rows, 1) {

		if (SP.at[m][0] > SP_val) {

			SP_val = SP.at[m][0];
			max_idx = m;

		}
	}

	return max_idx ;

	freeMat(SP);
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


void pivoting(Matrix _A, Matrix _b, Matrix P, int k) {

	Matrix Max = zeros(_A.rows - k, 1);
	Matrix SPval = zeros(_A.rows - k, 1);

	Max = row_max(_A, k);
	
	// 여기까지 각 행의 가장 큰 원소가 Max matrix에 저장되어있다.
	


}


void gaussElim_pivoing(Matrix _A, Matrix _b, Matrix _U, Matrix _d, Matrix _P) {

	if (_A.cols != _b.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR ! Dimension Error : starting gaussElim ");
		printf("\n*************************************************\n");

		exit(1);
	}

	//	rowExchange(_A, k, max_idx + k);
	double mult = 0;
	double nRow = 0;

	FOR_LOOP(k, 0, _A.rows-1, 1) {
		
		nRow = row_SP(_A,row_max(_A, k),k);
		
		rowExchange(_A, k, nRow + k);
		rowExchange(_b, k, nRow + k);
		rowExchange(_P, k, nRow + k);

		//printMat(_A,"\nrow exchange completed");

		FOR_LOOP(i, k + 1, _A.rows, 1) {

			mult = _A.at[i][k] / _A.at[k][k];
			_A.at[i][k] = 0;

			FOR_LOOP(j, k + 1, _A.rows, 1) {

				_A.at[i][j] = _A.at[i][j] - (mult * _A.at[k][j]);

			}

			_b.at[i][0] = _b.at[i][0] - mult * _b.at[k][0];

		}

		printMat(_A,"reduced form");
		printf("------------------------------------------------------------------------------------------------");
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


// 값이 이상하게 나옴
void LUdecomp_pivoting(Matrix _A, Matrix _L, Matrix _U, Matrix _P) {

	if ((_A.cols != _L.cols) || (_A.rows != _L.rows) || (_A.rows != _U.cols) || (_A.cols != _U.rows)) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at : 'LUdecomp' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	uint8_t row = _A.rows;
	double mult = 0;
	double nRow = 0;

	copyMat(_A, _U);

	FOR_LOOP(k, 0, row - 1, 1) {

		nRow = row_SP(_U, row_max(_U, k), k);

		rowExchange(_P, k, nRow + k);
		rowExchange(_U, k, nRow + k);

		printMat(_U, "Row exchanged U");

		FOR_LOOP(i, k + 1, row, 1) {

			mult = _U.at[i][k] / _U.at[k][k];

			_L.at[i][k] = mult;
			_U.at[i][k] = 0;

			FOR_LOOP(j, k + 1, row, 1) {

				_U.at[i][j] = _U.at[i][j] - mult * _U.at[k][j];
			}
		}

		printMat(_P, "P");
		printMat(_L, "L");
		printMat(_U, "U");

		printf("------------------------------------------------------------------------------------------------\n\n");
	}

}


void solveLU(Matrix _L, Matrix _U, Matrix _b, Matrix _x) {

	Matrix _y = zeros(_b.rows, 1);

	_y = fwdSub(_L, _b);
	_y = backSub(_U, _y);

	copyMat(_y, _x);

	freeMat(_y);
}


// For square matrix
Matrix invMat(Matrix _A) {

	if (_A.rows != _A.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at : 'invMat' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	uint8_t N = _A.rows;
	
	Matrix _L = eye(N, N);		Matrix _U = zeros(N, N);	Matrix _invA = zeros(N, N);
	Matrix _bk = zeros(N, 1);	Matrix _x = zeros(N, 1);	Matrix _I = eye(N, N);

	LUdecomp(_A, _L, _U);

	FOR_LOOP(k, 0, N, 1) {

		if (_U.at[k][k] == 0) {

			printf("\n*************************************************************");
			printf("\n  ERROR!!: dimension error at : 'Do not exist Inverse matrix'");
			printf("\n*************************************************************\n");
		}

		copyMat(vectExt(_I, k), _bk);

		solveLU(_L, _U, _bk, _x);

		vectins(_x, k, _invA);
				
	}
	
	return _invA;

	freeMat(_L);	freeMat(_U);	freeMat(_invA);
	freeMat(_bk);	freeMat(_x);	freeMat(_I);

}

