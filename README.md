# Numerical Programming

* Author : Han Minung
* Modified : 2022.10.03
* Mechanical and Control Engineering , Handong Global University



## Functions

### myNP.c (Function definition for implementation)

* Used macro
  * for loop iteration ( FOR_LOOP )
  * newline ( New_line )

* Taylor series to approximate Sin function , Cos function , Exponential function

  * Basic idea : taylor series approximation of functions

  <img src="https://user-images.githubusercontent.com/99113269/197355079-e1d2453d-b1e2-4396-bb15-2ad112b3dfa5.jpg" alt="Taylor-Series" style="zoom: 67%;" />

  * Input : [ rad ]
  * Input : [ deg ]

* Non-linear solver
  * Bisection method
  * Newton - Raphson method
  * Secant method

* Differentiation

  * 2 point differentiation

  * 3 point differentiation

  * 4 point differentiation

  * Input element : array  or  user-defined function


* Integration
  * Rectangular method
  * Midpoint method
  * Trapezoidal method
  * Simpson 13 ( 2 order polynomial - integrate )
    * parameter : function, two points , N value
  * Simpson 38 ( 3 order polynomial - integrate )
    * parameter : function, two points , N value
* Interpolation ( parameter : array elements , value to find )
  * Lagrange 1st order polynomial
  * Lagrange 2nd order polynomial
  * Lagrange nth order polynomial
    * order 3 : need 3 array elements to interpolate
    * order 4 : need 4 array elements to interpolate
    * ...
  * Newton 1st order polynomial
  * Newton 2nd order polynomial





------------------

### myMatrix.h (Function definition for implementation)

* **<u>*BAISC STRUCTURE OF MATRIX*</u>**

  ```c
  typedef struct { 
      
  	double** at;
  	int rows, cols;
      
  }Matrix;
  
  - need double pointer to create 2nd order matrix like MxN matrix
  ```

  -------------------------------

  **<u>Basic Function lists</u>**

  

* Matrix createMat( int _rows, int _cols )

  * Function return type : Matrix

  * Code structure

    ```c
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
    ```



* void freeMat ( Matrix _A )

  * Free allocated column memory
  * Free allocated row memory

* void printMat ( Matrix _A , const char* _name )

  * Print matrix elements
  * Print form : name - elements

* void initMat ( Matrix _A , double _val )

  * Fill elements of Matrix _A with 'val'

* Matrix zeros ( int _rows , int _cols )

  * parameter : n ( column ) , n ( row )

  * return value : Matrix 

  * Code structure

    ```c
    Matrix	zeros(int _rows, int _cols)
    {
    	Matrix Out = createMat(_rows, _cols);
    	initMat(Out,0);
        
    	return Out;
    }
    ```

    

* void copyMat ( Matrix _A , Matrix _B )

  * parameter : 2 matrix

  * copy elements of matrix A to those of matrix B

  * Code structure

    ```c
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
    ```

    

-----------------------------------

***<u>Main function lists</u>***

* void gaussElim ( Matrix _A , Matrix _b , Matrix _U , Matrix _d )

  * Solving a system of linear equation

  * parameter : user-defined matrix _A and _b

  * Usually, it is efficient to get matrix by typing elements in txt file and use 'txt2Mat' function.

    * It is essential to check the path of txt file !
    * It is essential to double check the size of matrices !
    * If required sizes of matrices to solve any problem are not satisfied, it is set to output an error.

  * A --> reduces echelon form (A) --> copy elements of A to U ( copyMat )

  * size of Matrices

    * _A : n x n
    * _b : n x 1
    * _U : n x n
    * _d : n x 1

    

    

    

    

    
    
    
    
    
    
    
    
    
    
    
