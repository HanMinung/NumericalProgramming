/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 11-02-2022
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

// Scalar multiplication to a matrix
extern Matrix multiScalar(double scalar, Matrix _A);

// Insert specific column of a matrix
extern void vectins(Matrix _v, int in, Matrix _A);

// Find max elements of each row
extern Matrix row_max(Matrix _A, uint8_t k);

extern int row_SP(Matrix _A, Matrix _rowMax, uint8_t k);

/*------------------------------------------------------------------*/

extern void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

extern void gaussElim_pivoing(Matrix _A, Matrix _b, Matrix _U, Matrix _d, Matrix _P);

extern void LUdecomp(Matrix _A, Matrix _L, Matrix _U);

extern void LUdecomp_pivoting(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

extern void solveLU(Matrix _L, Matrix _U, Matrix _b, Matrix _x);

extern Matrix invMat(Matrix _A);

extern void LUpivot(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

#endif