/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 11-02-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myNP.h"
#include "myMatrix.h"

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

Matrix	backSub(Matrix _U, Matrix _d)
{
	if (_U.rows != _d.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'backSub' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = createMat(_d.rows, 1);
	initMat(Out, 0);

	double sum = 0;

	FOR_LOOP_INV_INC(i,_U.rows-1, 0, 1){

		sum = 0;
		
		FOR_LOOP(j,i+1,_U.cols,1){
		
			sum = sum + (_U.at[i][j] * _d.at[j][0]);

		}
		_d.at[i][0] = (_d.at[i][0] - sum) / _U.at[i][i];
	}
	
	copyMat(_d, Out);

	return Out;
}

// 후에 LU 분해할 때 검증해볼것
Matrix  fwdSub(Matrix _L, Matrix _b) {

	if (_L.rows != _b.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'fwdSub' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = createMat(_b.rows, 1);
	Matrix _copyb = zeros(_b.rows, _b.cols);

	copyMat(_b, _copyb);

	double sum = 0;

	FOR_LOOP(i,0,_copyb.rows,1){
		sum = 0;  

		FOR_LOOP(j,0,i,1) {

			sum = sum + _L.at[i][j] * _copyb.at[j][0];
		}

		_copyb.at[i][0] = (_copyb.at[i][0] - sum) / _L.at[i][i];
	}

	copyMat(_copyb, Out);
	freeMat(_copyb);

	return Out;
}


Matrix Matproduct(Matrix _A, Matrix _B) {

	if (_A.cols != _B.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'Matproduct' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = zeros(_A.rows, _B.cols);
	
	FOR_LOOP(k,0,_A.rows,1){

		FOR_LOOP(i,0,_A.cols,1){

			FOR_LOOP(j, 0, _B.cols, 1) {

				Out.at[k][j] += _A.at[k][i] * _B.at[i][j];

			}
		}
	}

	return Out;
}