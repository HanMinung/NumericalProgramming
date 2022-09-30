/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [한민웅]
Created          : 26-03-2018
Modified         : 09-08-2022
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#include "myNP.h"

/*--------------------------------------------------------------*/
/*		            Handle Exceptional Cases 				   	*/
/*--------------------------------------------------------------*/

void ErrorHandling(char* message)
{

	fputs(message, stderr);

	fputc('\n', stderr);

	exit(1);

}

/*--------------------------------------------------------------*/
/*					    Basic Functions 					   	*/
/*--------------------------------------------------------------*/

double	Square_Number(double _x) {

	return _x * _x;
	
}

double	Cubic_Number(double _x) {

	return _x * _x * _x;

}

double	Four_Square_Number(double _x) {

	return _x * _x * _x * _x;

}

double	Factorial(double _x) {

	if (_x <= 1)	
		return 1;
	
	else
		return  _x * Factorial(_x - 1);

	}


/*		sin 함수의 taylor series 근사 : input [rad]		with predefined function	  */					

double	SinTaylor(double _x) {

	int N_max = 20;
	int N = 0;
	double epsilon = 1e-5;
	double S_N = 0, S_N_prev = 0, rel_chg = 0;

	do {
	
		N++;
		S_N_prev = S_N;
		S_N = 0;
		
		for (int k = 0; k < N; k++)	S_N += pow(-1, k) * pow(_x, 2 * k + 1) / Factorial(2 * k + 1);

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);
	
	} while (N < N_max && rel_chg >= epsilon);

	return S_N;

}

/*		 sin 함수의 taylor series 근사 : input [deg]		*/

double	SindTaylor(double _x) {

	return	SinTaylor(_x * PI / 180);
	return	SinTaylor(_x * PI / 180);

}


/*				Without Predefined function					*/	

double	SinTaylor_2(double _x) {

	int N_max = 20;
	int N = 0;
	double epsilon = 1e-5;
	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	
	do {

		N++;
		S_N_prev = S_N;
		S_N = 0;


		for (int k = 0; k < N; k++) {

			int sign_part = -1;

			/*		(-1)^ n		*/
			for (int i = 1; i < N; i++)
				sign_part *= -1;

			double pow_part = 1;

			/*		(-1)^ n		*/
			for (int i = 1; i <= 2 * k + 1; i++)
				pow_part *= _x;

			double fac_part = 1;

			/*		(-1)^ n		*/
			for (int i = 1; i <= 2 * k + 1; i++)
				fac_part *= i;

			S_N += sign_part * pow_part / fac_part;

		}

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);


	return S_N;

}


double Bisection(double func(double x), double a, double b, double tol) {

	double xn = 0;

	double epsilon = 1.0;
	int N_max = 100;
	int i = 0;

	do {

		printf("iteration : %d	|	x(n) = %lf	|	Tolerance : %.6f\n", i, xn, epsilon);
		
		i++;

		xn = (a + b) / 2;
		epsilon = fabs(func(xn));

		if (func(xn) == 0) {

			break;
		}

		if (func(a) * func(xn) < 0)	b = xn;

		if (func(b) * func(xn) < 0)	a = xn;

	} while (i < N_max && epsilon > tol);

		return xn;

}


double NewtonRaphson(double func(double x), double d_func(double x), double x0 , double tol) {
	
	double xn = x0;

	double epsilon = 1;
	int N_max = 100;	
	int i = 0;
	double h = 0;

	do {

		printf("iteration : %d	|	x(n) = %.10lf	|	Tolerance : %.10f\n", i, xn, epsilon);

		i++;
		
		if (d_func(xn) != 0) {

			h = -func(xn) / d_func(xn);

			xn = xn + h;

		}

		else {

			printf("Derivative value equals 0 -- ERROR ! \n\n\n");
			
			break;

		}

		epsilon = fabs(func(xn));

	} while (i < N_max && epsilon > tol);

	return xn;

}


double Secant(double func(double x), double x0, double x1, double tol) {

	double xn = x1;

	double epsilon = 1; //일단 tol보다는 큰 값을 설정해준다
	int Nmax = 1000;
	int i = 0;

	double df_xn = 0;
	double h = 0;

	do {

		printf("iteration : %d	|	x(n) = %lf	|	Tolerance : %.6f\n", i, xn, epsilon);
		i++;

		df_xn = ((func(x1) - func(x0)) / (x1 - x0));

		if (df_xn != 0) {

			h = -func(xn) / df_xn;
			xn = xn + h;
		}

		else {

			printf("Derivative value equals 0 -- ERROR ! \n");

			break;
		}

		epsilon = fabs(func(xn));

		x0 = x1;
		x1 = xn;

	} while (i < Nmax && epsilon > tol);

	return xn;

}


/*--------------------------------------------------------------------------------------------------------*/

void gradient1D_2Points(double x[], double y[], double dydx[], int m) {

	double h = x[1] - x[0];

	if (sizeof(x) == sizeof(y)) {

		printf("\nSAME VECTOR LENGTH ! \n\n");

	}

	else
	{
		printf("WARNING ! Length of two vectors is differenct ! \n");
		return;
	}

	/*		FORWARD VALUE		*/

	dydx[0] = (y[1] - y[0]) / h;

	/*		CENTRAL VALIE		*/

	for (int k = 1; k < m - 1; k++) {

		dydx[k] = (y[k + 1] - y[k - 1]) / (2 * h);

	}

	/*		BACKWARD VALUE		*/

	dydx[m - 1] = (y[m - 1] - y[m - 2]) / h;

}

// 표준함수로 설정 : 제일 많이 쓰는걸로
void gradient1D(double x[], double y[], double dydx[], int m) {

	double h = x[1] - x[0];		// calculate h value

	if (sizeof(x) == sizeof(y)) {

		printf("\nSAME VECTOR LENGTH ! \n\n");

	}

	else
	{
		printf("WARNING ! Length of two vectors is differenct ! \n");
		return;
	}

	// Forward Value : 3 points
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);

	//	CENTRAL VALUE	: Two points central value 
	for (int k = 1; k < m - 1; k++) {

		dydx[k] = (y[k + 1] - y[k - 1]) / (2 * h);

	}

	// Forward value : 3points
	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);

}


void gradient1D_4Points(double x[], double y[], double dydx[], int m) {

	double h = x[1] - x[0];

	if (sizeof(x) == sizeof(y)) {

		printf("\nSAME VECTOR LENGTH ! \n\n");

	}

	else
	{
		printf("WARNING ! Length of two vectors is differenct ! \n");
		return;
	}

	// Forward Value : 2개
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);
	dydx[1] = (-3 * y[1] + 4 * y[2] - y[3]) / (2 * h);

	FOR_LOOP(i, 2, m-2, 1) {

		dydx[i] = (y[i - 2] - 8 * y[i - 1] + 8 * y[i + 1] - y[i + 2]) / (12 * h);

	}

	// Backword Value : 2개
	dydx[m - 2] = (y[m - 4] - 4 * y[m - 3] + 3 * y[m - 2]) / (2 * h);
	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);

}


void gradientFunc(double func(const double x), double x[], double dydx[], int m) {

	double* y = (double*)malloc(sizeof(double) * m);


	// Your algorithm goes here
	if (sizeof(x) != sizeof(y)) {
		printf("WARNING!: Legth of x and y are different\n\n\n");
		return;
	}
	double h = x[1] - x[0]; // x축의 넓이



	//y 안에 임의의 f(x)의 값을 넣어줘야한다.
	for (int i = 0; i < m; i++) {

		y[i] = func(x[i]);

	}


	//forward
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);

	//central
	for (int k = 1; k < m - 1; k++) {

		dydx[k] = (y[k + 1] - y[k - 1]) / (2 * h);

	}
	//end of for

	//backwrad
	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);

	free(y);
}

// gradient2D_Func : 2계 미분 함수
void acceleration(double func(const double x), double x[], double dy2dx2[], int m) {

	double* y = (double*)malloc(sizeof(double) * m);

	if (sizeof(x) == sizeof(y)) {

		printf("\nSame Vector length !\n\n");
	}

	else {
		printf("WARNING!: Legth of x and y are different\n\n\n");

		return;
	}

	double h = x[1] - x[0]; // x축의 넓이


	for (int i = 0; i < m; i++) {

		y[i] = func(x[i]);

	}

	//Forward : Four-point forward differentiation
	dy2dx2[0] = (2*y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h * h);

	//Three - point central differenciation
	for (int k = 1; k < m - 1; k++) {

		dy2dx2[k] = (y[k + 1] - 2 * y[k] + y[k - 1]) / (h * h);

	}

	//Backward : Four-point backward differentiation
	dy2dx2[m - 1] = ( -y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / (h * h);


	free(y);

}

/*		Want to know specific Y value according to X_in		*/
void Function_call(double func(const double x), double X_in) {

	double y = func(X_in);
	//printf("Y according to X_in is : %lf", y);

}


void printVec(double* _vec, int _row) {

	for (int i = 0; i < _row; i++) {

		printf("Vector[%d] = %.5f \n", i, _vec[i]);

	}

	printf("\n\n");
}


/*--------------------------------------------------------------------------------------------------------*/


/*		m개의 index를 적분을 한다 : 0 ~ m-1까지 적분		*/
double Integral_Rect(double x[], double y[], int m) {

	double Integral = 0;

	for (int i = 0; i < m - 1; i++) {

		Integral += y[i] * (x[i + 1] - x[i]);

	}

	return Integral;

}


double Trapz(double x[], double y[], int m) {

	double	I = 0;

	FOR_LOOP(i,0,m-1,1) {

		I += (y[i] + y[i + 1]) * (x[i+1]-x[i]);

	}

	return I / 2;

}

// Distributed & Discrete data : Simpson 13	
double  simpson13(double x[], double y[], int m) {

	double I = 0;
	double Interval = x[1] - x[0];
	
	I = y[0] + y[m - 1] + 4 * y[m-1];

	FOR_LOOP(i, 1, m - 2, 2) {

		I += (4 * y[i]) + (2 * y[i+1]);

	}

	return I * (Interval / 3);

}

// not-evenly distributed data & function : Simpson 13
double Integral(double func(const double x), double a, double b, int N) {

	double I = 0;
	double Interval = (b - a) / N;

	I = func(a) + func(b) + 4 * func(b - Interval);

	FOR_LOOP(i, 1, N - 2, 2) {

		double xi = a + i * Interval;

		I += 4 * func(xi) + 2 * func(xi + Interval);

	}

	return I * (Interval / 3);

}


double	Integral_Simpson_38(double func(const double x), double a, double b, int N) {

	double I = 0;
	double h = (b - a) / N;

	I = func(a) + func(b) + 3 * func(b - 2 * h) + 3 * func(b - h);

	FOR_LOOP(i, 1, N - 4, 3) {

		double xi = a + (i * h);

		I += 3 * func(xi) + 3 * func(xi + h) + 2 * func(xi + 2 * h);

	}

	return h * I * 3 / 8;

}



