/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 19-10-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.cpp
----------------------------------------------------------------*/

#include "myNP.h"


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


Matrix	backSub(Matrix _U, Matrix _b)
{
	if (_U.rows != _b.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'backSub' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = createMat(_b.rows, 1);
	initMat(Out, 0);

	double sum = 0;

	FOR_LOOP_INV_INC(i,_U.rows-1, 0, 1){

		sum = 0;
		
		FOR_LOOP(j,i+1,_U.cols,1){
		
			sum = sum + (_U.at[i][j] * _b.at[j][0]);

		}

		_b.at[i][0] = (_b.at[i][0] - sum) / _U.at[i][i];
	}
	
	copyMat(_b, Out);

	return Out;

}


// Matrix forwardSubstitution ¸¸µé±â
