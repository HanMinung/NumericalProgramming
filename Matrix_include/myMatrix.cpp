/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 12-12-2022
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
		return zeros(0, 0);
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
		return zeros(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = zeros(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}


// Change array to a matrix form
Matrix	arr2Mat(double* _1Darray, int _rows, int _cols){

	Matrix Output = zeros(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
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


//diagnal을 nxn으로 뽑아내는 function 
Matrix getDiagonal(Matrix _A) {

	Matrix diag = zeros(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++) {
		int j = i;

		diag.at[i][j] = _A.at[i][j];
	}

	return diag;
}



//diagnal을 nx1으로 뽑아내는 함수
Matrix getDiagonal_vec(Matrix _A) {

	Matrix diag = zeros(_A.rows, 1);

	for (int i = 0; i < _A.rows; i++) {
		int j = i;

		diag.at[i][0] = _A.at[i][j];
	}

	return diag;
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

	//printf("================\n");
	//printf("|  At k = %d    |\n", k);
	//printf("================\n\n");
	//printMat(Out, "Max elements of each row");

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

	//printMat(SP, "SPval");
	
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


double lengthVec(Matrix _b) {

	if (_b.cols != 1) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'lengthVec' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	double length = 0;

	for (int i = 0; i < _b.rows; i++) {

		length += pow((_b.at[i][0]), 2);

	}

	return sqrt(length);
}


double	Norm(Matrix _c){
	
	double Out = 0;

	for (int i = 0; i < _c.rows; i++){

		Out += pow(_c.at[i][0],2);
	}

	Out = sqrt(Out);

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

		//printMat(_A,"reduced form");
		//printf("------------------------------------------------------------------------------------------------");
	}

	copyMat(_A, _U);	// A --> U
	copyMat(_b, _d);	// b --> d

}


Matrix gaussJor(Matrix _A, Matrix _b) {

	if (_A.cols != _b.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'gaussJor' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	double pivot = 0;
	int N = _A.rows;
	double mult = 0;

	Matrix Out = zeros(N, N);
	copyMat(_A, Out);

	FOR_LOOP(k,0,N,1){

		pivot = Out.at[k][k];

		FOR_LOOP(i,k,N,1){

			// 한 행에 대해서 pivot으로 나누고 시작함
			Out.at[k][i] = Out.at[k][i] / pivot;
		}

		_b.at[k][0] = _b.at[k][0] / pivot;

		// 모든 행에 대해서 가우스 소거를 진행하나, mult의 분자는 1로 고정된 상태
		FOR_LOOP(i,0,N,1){

			if (i != k) {

				mult = Out.at[i][k];
				Out.at[i][k] = 0;

				FOR_LOOP(j,k+1,N,1){

					Out.at[i][j] = Out.at[i][j] - (mult * Out.at[k][j]);
				}

				_b.at[i][0] = _b.at[i][0] - (mult * _b.at[k][0]);
			}

		}

		//printMat(Out, "gaussjor check matrix B in gaussJor");

		//탈출문 만약 전체 행이 0일때(충분히 작을때) 탈출한다
		double sum = 0;

		for (int i = 0; i < N; i++) {

			if (k + 1 >= N) {

				break;
			}
			sum += sum + fabs(Out.at[k + 1][i]);

		}

		if (sum < 1e-5) {

			//printMat(Out, "after gaussjor U in gaussJor");
			return Out;
		}
	}

	return Out;
}


// NxN 행렬
// Main statement : _L = eye(n,n); 
//					_U = zeros(n, n);
void LUdecomp_orig(Matrix _A, Matrix _L, Matrix _U) {

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


// pivoting
void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P) {

	if ((_A.cols != _L.cols) || (_A.rows != _L.rows) || (_A.rows != _U.cols) || (_A.cols != _U.rows)) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at : 'LUdecomp' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix diag = eye(_A.rows, _A.cols);
	uint8_t row = _A.rows;
	double mult = 0;
	double nRow = 0;

	copyMat(_A, _U);

	FOR_LOOP(k, 0, row - 1, 1) {

		nRow = row_SP(_U, row_max(_U, k), k);

		rowExchange(_P, k, nRow + k);
		rowExchange(_U, k, nRow + k);
		rowExchange(_L, k, nRow + k);

		//printMat(_U, "Row exchanged U");

		FOR_LOOP(i, k + 1, row, 1) {

			mult = _U.at[i][k] / _U.at[k][k];

			_L.at[i][k] = mult;
			_U.at[i][k] = 0;

			FOR_LOOP(j, k + 1, row, 1) {

				_U.at[i][j] = _U.at[i][j] - mult * _U.at[k][j];
			}
		}

		//printMat(_P, "P");
		//printMat(_L, "L");
		//printMat(_U, "U");

		//printf("------------------------------------------------------------------------------------------------\n\n");
	}

	copyMat(addMat(diag, _L),_L);

}


void solveLU_orig(Matrix _L, Matrix _U, Matrix _b, Matrix _x) {

	Matrix _y = zeros(_b.rows, 1);

	_y = fwdSub(_L, _b);
	_y = backSub(_U, _y);

	copyMat(_y, _x);

	freeMat(_y);
}


void solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix _x) {

	Matrix _y = zeros(_b.rows, 1);

	_y = Matproduct(_P, _y);

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

	LUdecomp_orig(_A, _L, _U);

	FOR_LOOP(k, 0, N, 1) {

		if (_U.at[k][k] == 0) {

			printf("\n*************************************************************");
			printf("\n  ERROR!!: dimension error at : 'Do not exist Inverse matrix'");
			printf("\n*************************************************************\n");
		}

		copyMat(vectExt(_I, k), _bk);

		solveLU_orig(_L, _U, _bk, _x);

		vectins(_x, k, _invA);
				
	}

	return _invA;

	freeMat(_L);	freeMat(_U);	freeMat(_invA);
	freeMat(_bk);	freeMat(_x);	freeMat(_I);

}


void QRHousehold(Matrix _A, Matrix _Q, Matrix _R) {

	//---------------------------------------------------------------------------------------
	//	_Q : matrix I
	//	_R : A가 eigen value를 유지한 채로 upper triangular martrix로 reduce되는 대상
	//	declaration lists : orthogonal matrix I
	//						H , e , c , v : to make Householder matrix
	//---------------------------------------------------------------------------------------

	double N = _A.rows;

	//Matrix for c, e, c, H, Q, R
	Matrix I = eye(N, N);		Matrix H = zeros(N, N);		Matrix e = zeros(N, 1);		
	Matrix c = zeros(N, 1);		Matrix v = zeros(N, 1);

	copyMat(_A, _R);
	copyMat(I, _Q);

	FOR_LOOP(k, 0, N-1, 1){

		c = vectExt(_R, k);

		FOR_LOOP(i, 0, k , 1){

			c.at[i][0] = 0;
		}

		if (c.at[k][0] < 0)	e.at[k][0] = 1;
		else				e.at[k][0] = -1;

		c = addMat(c, multiScalar(lengthVec(c), e));
		copyMat(c, v);

		Matrix vt = transpose(v);
		Matrix vvt = Matproduct(v, vt);
		Matrix vtv = Matproduct(vt, v);

		H = subMat(I, multiScalar(2.0 / vtv.at[0][0], vvt));

		// state update
		copyMat(Matproduct(_Q, H), _Q);
		copyMat(Matproduct(H, _R), _R);

	}

	freeMat(I);		freeMat(H);		freeMat(e);
	freeMat(c);		freeMat(v);

}


Matrix eig(Matrix _A) {

	// checking dimension error
	Dim_error(_A,"eig");

	Matrix Out = zeros(_A.rows, _A.cols);
	Matrix Out_eigvalue = zeros(_A.rows, 1);

	copyMat(_A, Out);
	uint16_t n = _A.rows;		int k = 0;
	double sum = 0;  

	do {

		Matrix _Q = zeros(_A.rows, _A.cols);
		Matrix _R = zeros(_A.rows, _A.cols);
		sum = 0;

		QRHousehold(Out, _Q, _R);

		copyMat(Matproduct(_R, _Q), Out);
		//printMat(Out, "Out (after)");
		//printf("\n\n");

		freeMat(_Q);	freeMat(_R);

		FOR_LOOP(k, 0, n-1, 1){
			FOR_LOOP(i,k+1,n,1){

				sum += fabs(Out.at[i][k]);
			}
		}

		k++;

	} while(sum > 1e-7);
	
	Out_eigvalue = getDiagonal_vec(Out);

	//printf("totla loop k = %d\n", k);
	freeMat(Out);

	return Out_eigvalue;

}


Matrix eigVec(Matrix _A)
{

	Matrix eigvect = zeros(_A.rows, _A.cols);	
	Matrix I = eye(_A.rows, _A.cols);	
	Matrix vec0 = zeros(_A.rows, 1);
	uint16_t idx = 0;

	// checking dimension error
	Dim_error(_A, "eigVec");

	Matrix eigVal = eig(_A);

	FOR_LOOP(i,0,_A.cols,1){
	
		Matrix _U = zeros(_A.rows, _A.cols);
		Matrix _B = subMat(_A, multiScalar(eigVal.at[i][0], I));

		//printMat(_A, "A");
		//printMat(_B, "B = A - (lamda)*I");
		//printf("eigen value [%d] : %lf\n\n", i, eigVal.at[i][0]);

		copyMat(gaussJor(_B, vec0), _U);
		//printMat(_U, "U");

		FOR_LOOP(j,0, eigvect.rows,1){
			
			if (fabs(_U.at[j][j]) < 1e-5){

				Matrix change = multiScalar(-1.0, vectExt(_U, j));
				
				change.at[j][0] = 1.0;

				//eignen vector nomalize해주는 과정
				double norm = lengthVec(change);

				//normalize
				copyMat(multiScalar(1 / norm, change), change);

				//printMat(change, "change");
				vectins(change, idx, eigvect);

				idx++;

			}

		}
		//printMat(eigvect, "eigvect");
		
		freeMat(_U);	freeMat(_B);

	}

	return eigvect;
}



//Curve Fitting
double LinRegress(Matrix _x, Matrix _y, double find_x) {

	if ((_x.rows != _y.rows) || (_x.cols != _y.cols)) {
		
		print_str("Dimension error at 'LinRegress' func");

		exit(1);
	}

	double find_y = 0;

	double S_x = 0;
	double S_y = 0;
	double S_xx = 0;
	double S_xy = 0;
	double n = _x.rows;
	double a1, a0 = 0;

	FOR_LOOP(i,0,_x.rows,1){

		S_x += _x.at[i][0];
		S_y += _y.at[i][0];
		S_xx += _x.at[i][0] * _x.at[i][0];
		S_xy += _x.at[i][0] * _y.at[i][0];

	}

	a1 = (n * S_xy - S_x * S_y) / (n * S_xx - S_x * S_x);
	a0 = (S_xx * S_y - S_x * S_xy) / (n * S_xx - S_x * S_x);

	//printf("a1 = %lf \n\na0 = %lf\n\n", a1, a0);

	find_y = a0 + a1 * find_x;

	return find_y;
}


Matrix polyfit(Matrix _x, Matrix _y, int n) {

	if ((_x.cols != _y.cols) || (_x.rows != _y.rows)) {

		print_str("Dimension error at 'polyfit' func");

		exit(1);
	}

	double Sx_temp = 0;
	double Sxy_temp = 0;

	Matrix Sxy = zeros(n + 1, 1);			// A transpose * y
	Matrix Sx = zeros(2 * n + 1, 1);
	Matrix S = zeros(n + 1, n + 1);			// A tranpose * A
	Matrix z = zeros(n + 1, 1);				// least square fit

	// create matrix A
	FOR_LOOP_INC(i, 0, 2 * n, 1) {

		Sx_temp = 0;

		FOR_LOOP(k, 0, _x.rows, 1) {

			Sx_temp += pow(_x.at[k][0], i);
		}

		Sx.at[i][0] = Sx_temp;
	}


	FOR_LOOP_INC(i, 0, n, 1) {
		FOR_LOOP_INC(j, 0, n, 1) {

			S.at[i][j] = Sx.at[2*n - i - j][0];
		}
	}

	FOR_LOOP_INV_INC(j, n, 0, 1) {
		
		Sxy_temp = 0;

		FOR_LOOP(k, 0, _x.rows, 1) {

			Sxy_temp += pow(_x.at[k][0], j) * _y.at[k][0];
		}

		Sxy.at[n - j][0] = Sxy_temp;


	}

	z = Matproduct(invMat(S), Sxy);

	freeMat(Sxy);	freeMat(Sx);	freeMat(S);

	return z;
}

// ----------------------------	 JACOBIAN : Non linear solver	---------------------------------=

Matrix solve_nonlinear(Matrix X, Matrix F(const double x, const double y), Matrix J(const double x, const double y)){

	double epsilon = 10.00;		double tol = 1e-9;

	Matrix L = eye(X.rows, X.rows);
	Matrix U = zeros(X.rows, X.rows);
	Matrix H = zeros(X.rows, 1);

	double x = X.at[0][0];
	double y = X.at[1][0];

	int i = 0;			int Nmax = 10;
	double x_tol = 0;	double y_tol = 0;

	do {

		printf("ITERATION : %d \n", i);
		//printMat(F(x, y), "F");

		// Jacobian matrix LU decomposition 
		LUdecomp_orig(J(x, y), L, U);

		solveLU_orig(L, U, multiScalar(-1.0, F(x, y)), H);
		
		//update
		X = addMat(X, H);

		x = X.at[0][0];   y = X.at[1][0];

		epsilon = Norm(F(x, y));

		i++;

	} while (i<Nmax && epsilon>tol);

	freeMat(H);

	New_line(2);

	return X;

}




