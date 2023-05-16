# AlgLinComp
Repository for Computational Linear Algebra Codes

# Instructions to Use The Code

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

# Legend for ICOD
## P1Task01.py
ICOD = 1 -> LU Decomposition

ICOD = 2 -> Cholesky Decomposition

ICOD = 3 -> Jacobi Iterative Method

ICOD = 4 -> Gauss-Seidel Iterative Method

## P1Task02.py
ICOD = 1 -> Power Method

ICOD = 2 -> Jacobi Method

Note: If you choose a very low tolerance, the algorithm might take longer to converge, especially using the Jacobi Method