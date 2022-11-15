/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : HanMinung
Created          : 26-03-2018
Modified         : 10-12-2022
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/
#define _CRT_SECURE_NO_WARNIGNS

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

void print_str(const char* result) {

	printf("\n****************************************************\n");
	printf("	    Result of % s          \n",result);
	printf("****************************************************\n\n\n");

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


double	CosTaylor(double _x) {

	int N_max = 20;
	int N = 0;
	double epsilon = 1e-5;
	double S_N = 0, S_N_prev = 0, rel_chg = 0;

	do {

		N++;
		S_N_prev = S_N;
		S_N = 0;

		for (int k = 0; k < N; k++)	S_N += pow(-1, k) * pow(_x, 2 * k) / Factorial(2 * k);

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	return S_N;

}

/*		 sin 함수의 taylor series 근사 : input [deg]		*/

double	SindTaylor(double _x) {

	return	SinTaylor(_x * PI / 180);

}


double	CosdTaylor(double _x) {

	return	CosTaylor(_x * PI / 180);

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

double ExpTaylor(double _x){

	int N_max = 20;
	int N = 0;
	double epsilon = 1e-5;
	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	
	do {
	
		N++;
		S_N_prev = S_N;
		S_N = 0;
	
		for (int k = 0; k < N; k++)	S_N += pow(_x, k) / Factorial(k);
	
		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);
	
	} while (N < N_max && rel_chg >= epsilon);
	
	return S_N;

}


double Bisection(double func(double x), double a, double b, double tol) {

	double xn = 0;

	double epsilon = 1.0;
	int N_max = 100;
	int i = 0;

	if (func(a) * func(b) < 0) printf("\nSolution exists between a , b \n\n");
	else printf("\nSolution does not exist between a , b \n\n");

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

// 2point central difference
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

// 3 point forward difference
void gradient1D_3Points(double x[], double y[], double dydx[], int m) {

	double h = x[1] - x[0];		// calculate h value

	if (sizeof(x) == sizeof(y)) {

		printf("\nSAME VECTOR LENGTH ! \n\n");

	}

	else
	{
		printf("WARNING ! Length of two vectors is differenct ! \n");
		return;
	}

	//	CENTRAL VALUE	: 3 points forward difference
	for (int k = 0; k < m - 2; k++) {

		dydx[k] = (-3 * y[k] + 4 * y[k + 1] - y[k + 2]) / (2 * h);

	}

	// 마지막 값들 : 3points backward 
	dydx[m - 2] = (y[m - 4] - 4 * y[m - 3] + 3 * y[m - 2]) / (2 * h);
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
	double h = x[1] - x[0]; // x interval



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


// 2계 미분함수 : x (array) , y (function)
void gradient2D(double func(const double x), double x[], double dy2dx2[], int m) {

	double* y = (double*)malloc(sizeof(double) * m);

	if (sizeof(x) == sizeof(y)) {

		printf("\nSame Vector length !\n\n");
	}

	else {
		printf("WARNING!: Length of x and y is different\n\n\n");

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



// 2계 미분함수 : x,y 모두 array
void gradient2D_array(double x[], double y[],double dy2dx2[], int m) {

	double h = 0;

	if (sizeof(x) != sizeof(y)){ 
		
		printf("Warning ! : Length of x and y is different\n\n\n"); 

		return;
	}

	h = x[1] - x[0];

	//Forward : Four-point forward differentiation
	dy2dx2[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h * h);

	//Three - point central differenciation
	FOR_LOOP(k,1,m-1,1) {

		h = x[k + 1] - x[k];

		dy2dx2[k] = (y[k + 1] - 2 * y[k] + y[k - 1]) / (h * h);

	}

	//Backward : Four-point backward differentiation
	
	h = x[m] - x[m-1];

	dy2dx2[m - 1] = (-y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / (h * h);

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
	double h = x[1] - x[0];

	I = y[0] + y[m-1] + 4 * y[m-2];

	FOR_LOOP(i, 1, m-2, 2) {

		I += (4 * y[i]) + (2 * y[i + 1]);

	}

	return I * (h / 3);

}

// Call function : Simpson 13
double Integral_13(double func(const double x), double a, double b, int N) {

	double I = 0;
	double Interval = (b - a) / N;

	I = func(a) + func(b) + 4 * func(b - Interval);

	FOR_LOOP(i, 1, N - 2, 2) {

		double xi = a + i * Interval;

		I += 4 * func(xi) + 2 * func(xi + Interval);

	}

	return I * (Interval / 3.0);

}


double	Integral_38(double func(const double x), double a, double b, int N) {

	double I = 0;
	double h = (b - a) / N;

	I = func(a) + func(b) + 3 * func(b - 2 * h) + 3 * func(b - h);

	// index : 0 ~ N-1
	FOR_LOOP(i, 1, N - 4, 3) {

		double xi = a + (i * h);

		I += 3 * func(xi) + 3 * func(xi + h) + 2 * func(xi + 2 * h);

	}

	return h * I * (3.0) / (8.0);

}


/*--------------------------------------------------------------------------------------------------------*/

double lagrange_1st(double x, double x0, double x1, double y0, double y1) {

	double f_x = y0 * (x - x1) / (x0 - x1) + y1 * (x - x0) / (x1 - x0);

	return f_x;

}

double lagrange_2nd(double x, double x0, double x1,double x2, double y0, double y1, double y2) {

	double f_x = y0 * ((x-x1) * (x-x2))/((x0-x1)*(x0-x2)) + y1 * ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) + y2 * ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1));

	return f_x;

}


double lagrange_nth(double x[], double y[], double xx,int length) {
	//int length = sizeof(_x)/sizeof(double);

	double sum = 0.0;
	double product = 0.0;

	// 인덱스 수 찾는 과정이 오류가 남
	FOR_LOOP(i, 0, length, 1) {

		double product = y[i];
		FOR_LOOP(j,0, length,1) {

			if (i != j) {

				product = product * ((xx - x[j]) / (x[i] - x[j]));

			}
		}
		sum = sum + product;
	}

	return sum;
}


double Newton_1st(double x, double x0, double x1, double y0, double y1) {

	double f_x = y0 + ((y1 - y0) / (x1 - x0)) * (x - x0);
	
	return f_x;

}


double Newton_2nd(double x, double x0, double x1, double x2, double y0, double y1, double y2) {

	double a0 = y0;
	double a1 = (y1 - y0) / (x1 - x0);
	double fac = (y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0);
	double a2 = (fac) / (x2 - x0);
	
	double f_x = a0 + a1 * (x - x0) + a2 * (x - x0) * (x - x1);

	return f_x;

}


/*--------------------------------------------------------------------------------------------------------*/

// select which method is used to solve ODE
// --- method ---
//	  EU	0
//	  RK2	1
//	  RK3   2
//	  RK4	3
void ode(double func(const double t, const double v), double t0, double tf, double h, double v0, uint8_t method) {
	
	switch (method) {

		case(0) :

			odeEU(func, t0, tf, h, v0);
			break;

		case(1) :

			odeRK2(func, t0, tf, h, v0);
			break;
		
		case(2):

			odeRK3(func, t0, tf, h, v0);
			break;

		case(3):

			odeRK4(func, t0, tf, h, v0);
			break;

	}
}


// 오일러 method
void odeEU(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double slope = 0;
	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1)); //t  N+1 = 101 : 0~100 ==> 101개 
	v = (double*)malloc(sizeof(double) * (N + 1));

	t[0] = t0;
	v[0] = v0;

	FOR_LOOP(i,0,N,1){
		// first slope 
		slope = func(t[i], v[i]);

		// update
		v[i + 1] = v[i] + slope * h;
		t[i + 1] = t[i] + h;

	}

	print_str("odeEU");

	FOR_LOOP(i,0,N+1,1){

		printf("The odeEU v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp; 
	_fp = fopen("odeEU.txt", "wt"); 

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1){

		fprintf(_fp, "%f\n", v[i]);  
	}

	fclose(_fp);

	free(t);	free(v);
}

// 오일러 modified
void odeEM(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1));		//t  N+1 = 101 : 0~100 ==> 101개 
	v = (double*)malloc(sizeof(double) * (N + 1));

	//intial vaule
	t[0] = t0;	
	v[0] = v0;

	// 1st slope & next slope : average --> Euler modified
	double slope_f = 0;
	double slope_b = 0;
	
	FOR_LOOP(i,0,N,1){

		slope_f = func(t[i], v[i]);
		v[i + 1] = v[i] + slope_f * h;

		t[i + 1] = t[i] + h;
		slope_b = func(t[i + 1], v[i + 1]);

		// update y[i+1]
		v[i + 1] = v[i] + (slope_b + slope_f) * h / 2;

	}

	print_str("odeEM");

	FOR_LOOP(i, 0, N+1, 1){

		printf("The odeEM t[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeEM.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);

}


void odeRK2(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1));
	v = (double*)malloc(sizeof(double) * (N + 1));

	t[0] = t0;	v[0] = v0;

	double K_1 = 0;	double K_2 = 0;

	//default
	double alpha = 1;
	double beta = alpha;

	double c2 = 1 / (2 * alpha);
	double c1 = 1 - c2;

	FOR_LOOP(i,0,N,1){

		K_1 = func(t[i], v[i]);
		K_2 = func(t[i] + alpha * h, v[i] + beta * h * K_1);

		v[i + 1] = v[i] + (c1 * K_1 + c2 * K_2) * h;

		t[i + 1] = t[i] + h;
	}

	print_str("RK2 method");

	FOR_LOOP(i,0,N+1,1){

		printf("The RK2 v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeRK2.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1){

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);

}


void odeRK3(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1));
	v = (double*)malloc(sizeof(double) * (N + 1));

	t[0] = t0;	v[0] = v0;

	double K_1 = 0;		
	double K_2 = 0;		
	double K_3 = 0;

	double alpha2 = (1.0 / 2.0);		
	double alpha3 = 1.0;
	
	double beta21 = (1.0 / 2.0);		
	double beta31 = -1.0;		
	double beta32 = 2.0;
	
	double c1 = 1.0 / 6.0;	
	double c2 = 4.0 / 6.0;	
	double c3 = 1.0 / 6.0;


	FOR_LOOP(i, 0, N, 1) {

		K_1 = func(t[i], v[i]);
		K_2 = func(t[i] + alpha2 * h, v[i] + beta21 * K_1 * h);
		K_3 = func(t[i] + alpha3 * h, v[i] + beta31 * K_1 * h + beta32 * K_2 * h);

		v[i + 1] = v[i] + (c1 * K_1 + c2 * K_2 + c3 * K_3) * h;

		t[i + 1] = t[i] + h;
	}

	print_str("RK3 method");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The RK3 v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeRK3.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N+1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);

}


void odeRK4(double func(double t, double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1)); //t  N+1 = 101 : 0~100 ==> 101개 
	v = (double*)malloc(sizeof(double) * (N + 1));

	//intial vaule
	if (t != NULL && v != NULL) {
		t[0] = t0;
		v[0] = v0;
	}
	else
		printf("NULL POINTER");

	double K1 = 0;
	double K2 = 0;
	double K3 = 0;
	double K4 = 0;

	double alpha2 = 1.0 / 2.0;
	double alpha3 = 1.0 / 2.0;
	double alpha4 = 1.0;
	double beta21 = 1.0 / 2.0;
	double beta31 = 0;
	double beta32 = 1.0 / 2.0;
	double beta41 = 0.0;
	double beta42 = 0.0;
	double beta43 = 1.0;

	double c1 = 1.0 / 6.0;
	double c2 = 2.0 / 6.0;
	double c3 = 2.0 / 6.0;
	double c4 = 1.0 / 6.0;

	FOR_LOOP(i, 0, N, 1) {

		K1 = func(t[i], v[i]);
		K2 = func(t[i] + 0.5 * h, v[i] + 0.5 * K1 * h);
		K3 = func(t[i] + 0.5 * h, v[i] + 0.5 * K2 * h);
		K4 = func(t[i] + 0.5 * h, v[i] + 0.5 * K3 * h);

		v[i + 1] = v[i] + (c1 * K1 + c2 * K2 + c3 * K3 + c4 * K4) * h;
		t[i + 1] = t[i] + h;
	}

	print_str("RK4 method");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The RK4 v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeRK4.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N+1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);		free(v);
	
}