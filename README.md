# AlgLinComp
Repository for Computational Linear Algebra Codes

## Instructions to Use The Code - P1-Tasks

You need to have the matrix you want to work on inside your directory

That also applies to the vectors which are needed to solve the linear equations

Towards the end of every python script there is a call to a function:

> matrix_A = read_matrix_A('Matriz_A.dat')

You can change that argument to the name of the file holding your matrix

That also applies to the vectors

> vector_B_1 = read_vector_B('Vetor_B_01.dat')

Furthermore, when you call the script:

> python P1Task01.py

You must add the ICOD and the tol, which are the code of the method you are trying to use and the tolerance which will be used to stop the iterative methods. Eg:

> python P1Task01.py 1 0.001

## Legend for ICOD
### P1Task01.py
ICOD = 1 -> LU Decomposition

ICOD = 2 -> Cholesky Decomposition

ICOD = 3 -> Jacobi Iterative Method

ICOD = 4 -> Gauss-Seidel Iterative Method

### P1Task02.py
ICOD = 1 -> Power Method

ICOD = 2 -> Jacobi Method

Note: If you choose a very low tolerance, the algorithm might take longer to converge, especially using the Jacobi Method

## Instructions to Use The Code - P2-Tasks

This one is more simple to execute than P1-Tasks, you only have to execute the python scripts to get the results.

### P2Task01.py

This python code calculates the roots of a system of non linear equations using both the Newtons Method and Broydens method. 

The system is being defined inside the "func" function. 

There is also the definition of the Jacobian of that function, required for the methods.

If you want to calculate the roots for a different system of non linear equations, you must change both the functions "func" and "jacobian".

### P2Task02.py

This python code calculates the integral defined by the limits "a" and "b" using an adaptative integration function and the gauss quadrature method. 

It comparates both solutions and outputs it into the terminal.

It also calculates two different integrals, m0 and m2, having its particularities each.

### P2Task03.py

This python code calculates the EDO in the problem of the falling body into the water

Just by executing it you can find the tests and iterations available
