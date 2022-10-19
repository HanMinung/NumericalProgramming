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

