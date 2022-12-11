/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 12-12-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

#define		FOR_LOOP(i,Initial,Final,Interval)			for(int i = Initial ; i < Final ; i = i + Interval)
#define		FOR_LOOP_INC(i,Initial,Final,Interval)		for(int i = Initial ; i <= Final ; i = i + Interval)
#define		FOR_LOOP_INV(i,Initial,Final,Interval)		for(int i = Initial ; i > Final ; i = i - Interval)
#define		FOR_LOOP_INV_INC(i,Initial,Final,Interval)	for(int i = Initial ; i >= Final ; i = i - Interval)
#define		New_line(N)		FOR_LOOP(i,0,N,1){printf("\n");}


typedef struct { 
	double** at;
	int rows, cols;
}Matrix;


//using namespace std;

// Create Matrix with specified size
extern Matrix createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern void	freeMat(Matrix _A);

// Create a matrix from a text file
extern Matrix txt2Mat(std::string _filePath, std::string _fileName);

// Change array to matrix
Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);

extern void	printMat(Matrix _A, const char* _name);

// initialization of Matrix elements
extern void	initMat(Matrix _A, double _val);

// Create matrix of all zeros
extern Matrix zeros(int _rows, int _cols);

// Create matrix of all ones
extern Matrix ones(int _rows, int _cols);

// Create identity 
extern Matrix eye(int _rows, int _cols);

// Create Transpose matrix
extern Matrix transpose(Matrix _A);

// Copy matrix ( A --> B )
extern void copyMat(Matrix _A, Matrix _B);

// Extract specific row of a column
extern Matrix rowExt(Matrix _A, int in);

// Insert specific row to a Matrix
extern void rowIns(Matrix _v, int in, Matrix _A);

// Exchange two rows in Matrix
extern void rowExchange(Matrix _A, int seq1, int seq2);

// Extract specific column of a matrix 
extern Matrix vectExt(Matrix _A, int in);

// Insert specific column to a matrix
extern void vectins(Matrix _v, int in, Matrix _A);

//diagnal을 nxn으로 뽑아내는 function 
extern Matrix getDiagonal(Matrix _A);

//diagnal을 nx1으로 뽑아내는 함수
extern Matrix getDiagonal_vec(Matrix _A);

// Scalar multiplication to a matrix
extern Matrix multiScalar(double scalar, Matrix _A);

// Insert specific column of a matrix
extern void vectins(Matrix _v, int in, Matrix _A);

// Find max elements of each row
extern Matrix row_max(Matrix _A, uint8_t k);

extern int row_SP(Matrix _A, Matrix _rowMax, uint8_t k);

// Norm of Vector
extern double lengthVec(Matrix _b);

extern double Norm(Matrix _c);

/*------------------------------------------------------------------*/

extern void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

extern void gaussElim_pivoing(Matrix _A, Matrix _b, Matrix _U, Matrix _d, Matrix _P);

extern Matrix gaussJor(Matrix _A, Matrix _b);

extern void LUdecomp_orig(Matrix _A, Matrix _L, Matrix _U);

extern void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

extern void solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix _x);

extern void solveLU_orig(Matrix _L, Matrix _U, Matrix _b, Matrix _x);

extern Matrix invMat(Matrix _A);

extern void QRHousehold(Matrix _A, Matrix _Q, Matrix _R);

extern Matrix eig(Matrix _A);

extern Matrix eigVec(Matrix _A);

extern double LinRegress(Matrix _x, Matrix _y, double find_x);

extern Matrix polyfit(Matrix _x, Matrix _y, int n);

/*------------------------------------------------------------------*/

Matrix solve_nonlinear(Matrix X, Matrix F(const double x, const double y), Matrix J(const double x, const double y));

#endif