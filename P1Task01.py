import sys

def read_matrix_A(matrix_A_file):
    matrix_A = []
    with open(matrix_A_file, 'r') as f:
        lines = f.readlines()
        size_of_matrix = len(lines)
        for l in lines:
            line = l.split()
            for i in range(size_of_matrix):
                line[i] = float(line[i])
            matrix_A.append(line)
    f.close()
    return matrix_A

def read_vector_B(vector_B_file):
    vector_B = []
    with open(vector_B_file, 'r') as f2:
        lines = f2.readlines()
        for l in lines:
            vector_B.append(float(l))
    f2.close()
    return vector_B

def symmetrical(matrix_A):
    # Check whether matrix is symmetrical or not
    n = len(matrix_A)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if matrix_A[i][j] != matrix_A[j][i]:
                return False
    return True

def LU_Only_Decomposition(matrix_A):
    n = len(matrix_A)
    for i in range(n):
        for j in range(i + 1, n):
            if matrix_A[i][i] == 0:
                sys.exit("Null pivot detected, LU decomposition without pivoting is not possible")
            matrix_A[j][i] = matrix_A[j][i] / matrix_A[i][i]
        for k in range(i + 1, n):
            for l in range(i + 1, n):
                matrix_A[l][k] -= (matrix_A[l][i] * matrix_A[i][k])
    return matrix_A
    
def Cholesky_Only_Decomposition(matrix_A):
    n = len(matrix_A)
    # In order to do cholesky, the matrix needs to be symmetrical, definite and positive
    if not symmetrical(matrix_A):
        sys.exit("Matrix must be definite positive symmetric")
    for i in range(n):
        for j in range(i):
            matrix_A[i][i] -= (matrix_A[i][j] ** 2)
        if matrix_A[i][i] <= 0:
            sys.exit("Matrix must be definite positive symmetric")
        matrix_A[i][i] = matrix_A[i][i] ** (1/2)
        for k in range(i + 1, n):
            matrix_A[k][i] = matrix_A[i][k]
            for l in range(i):
                matrix_A[k][i] -= (matrix_A[i][l] * matrix_A[k][l])
            matrix_A[k][i] = matrix_A[k][i] / matrix_A[i][i]
    return matrix_A

def LU_Decomposition(matrix_A, vector_B):
    n = len(matrix_A)
    # Turning B into Y
    for i in range(1, n):
        for j in range(i):
            vector_B[i] -= (matrix_A[i][j] * vector_B[j])
    # Turning B into X
    vector_B[n - 1] = vector_B[n - 1] / matrix_A[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, n):
            vector_B[i] -= (matrix_A[i][j] * vector_B[j])
        vector_B[i] = vector_B[i] / matrix_A[i][i]
    return vector_B

def Cholesky(matrix_A, vector_B):
    n = len(matrix_A)
    # Turning B into Y
    vector_B[0] = vector_B[0] / matrix_A[0][0]
    for i in range(1, n):
        for j in range(i):
            vector_B[i] -= (matrix_A[i][j] * vector_B[j])
        vector_B[i] = vector_B[i] / matrix_A[i][i]
    # Turning B into X
    vector_B[n - 1] = vector_B[n - 1] / matrix_A[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, n):
            vector_B[i] -= (matrix_A[j][i] * vector_B[j])
        vector_B[i] = vector_B[i] / matrix_A[i][i]
    return vector_B

def diagonal_dominant(matrix_A):
    n = len(matrix_A)
    for i in range(n):
        line_sum = 0
        column_sum = 0
        for j in range(n):
            if i == j:
                continue
            line_sum += abs(matrix_A[i][j])
            column_sum += abs(matrix_A[j][i])
        if abs(matrix_A[i][i]) < column_sum or abs(matrix_A[i][i]) < column_sum:
            return False
    return True

def residue(x0, x1, n):
    d = [0 for _ in range(n)]
    numerator = 0
    for i in range(n):
        d[i] = (x1[i] - x0[i]) ** 2
        numerator += d[i]
    numerator = numerator ** (1/2)

    denominator = 0
    for index, j in enumerate(x1):
        d[index] = j ** 2
        denominator += d[index]
    denominator = denominator ** (1/2)

    return numerator / denominator

def jacobi(matrix_A, vector_B, tol):
    if not diagonal_dominant(matrix_A):
        print(0, 0)
        sys.exit("Matrix must be diagonal dominant")
    n = len(matrix_A)
    # Initialazing x0 and x1
    x0 = [1 for _ in range(n)]
    x1 = [0 for _ in range(n)]
    storage = []  # Storage of errors
    while True:
        for i in range(n):
            x1[i] = vector_B[i] / matrix_A[i][i]
            for j in range(n):
                if j == i:
                    continue
                x1[i] -= (matrix_A[i][j] * x0[j]) / matrix_A[i][i]
        erro = residue(x0, x1, n)
        storage.append(erro)
        if erro < tol:
            num_iter = len(storage)
            return x1, storage, num_iter
        x0 = x1[:] # Copy elements from x1 to x0

def gauss_seidel(matrix_A, vector_B, tol):
    if not diagonal_dominant(matrix_A):
        print(0, 0)
        sys.exit("Matrix must be diagonal dominant")
    n = len(matrix_A)
    x0 = [1 for _ in range(n)]
    x1 = [0 for _ in range(n)]
    storage = []
    while True:
        for i in range(n):
            x1[i] = vector_B[i] / matrix_A[i][i]
            for j in range(i):
                x1[i] -= (matrix_A[i][j] * x1[j]) / matrix_A[i][i]
            for j in range(i + 1, n):
                x1[i] -= (matrix_A[i][j] * x0[j]) / matrix_A[i][i]
        erro = residue(x0, x1, n)
        storage.append(erro)
        if erro < tol:
            num_iter = len(storage)
            return x1, storage, num_iter
        x0 = x1[:]

ICOD = int(sys.argv[1])
tol  = float(sys.argv[2])

matrix_A   = read_matrix_A('Matriz_A.dat')
vector_B_1 = read_vector_B('Vetor_B_01.dat')
vector_B_2 = read_vector_B('Vetor_B_02.dat')
vector_B_3 = read_vector_B('Vetor_B_03.dat')

if ICOD == 1:
    print("You chose LU Decomposition")
    LU_Only_Decomposition(matrix_A)
    print("Solution for Vector B 1 is: ")
    print(LU_Decomposition(matrix_A, vector_B_1))
    print("Solution for Vector B 2 is: ")
    print(LU_Decomposition(matrix_A, vector_B_2))
    print("Solution for Vector B 3 is: ")
    print(LU_Decomposition(matrix_A, vector_B_3))
elif ICOD == 2:
    print("You chose Cholesky Decomposition")
    Cholesky_Only_Decomposition(matrix_A)
    print("Solution for Vector B 1 is: ")
    print(Cholesky(matrix_A, vector_B_1))
    print("Solution for Vector B 2 is: ")
    print(Cholesky(matrix_A, vector_B_2))
    print("Solution for Vector B 3 is: ")
    print(Cholesky(matrix_A, vector_B_3))
elif ICOD == 3:
    print("You chose Jacobi Iteration Method")
    print("Solution for Vector B 1 is: ")
    print(jacobi(matrix_A, vector_B_1, tol))
    print("Solution for Vector B 2 is: ")
    print(jacobi(matrix_A, vector_B_2, tol))
    print("Solution for Vector B 3 is: ")
    print(jacobi(matrix_A, vector_B_3, tol))
elif ICOD == 4:
    print("You chose Gauss Seidel Iteration Method")
    print("Solution for Vector B 1 is: ")
    print(gauss_seidel(matrix_A, vector_B_1, tol))
    print("Solution for Vector B 2 is: ")
    print(gauss_seidel(matrix_A, vector_B_2, tol))
    print("Solution for Vector B 3 is: ")
    print(gauss_seidel(matrix_A, vector_B_3, tol))
else:
    sys.exit("Wrong ICOD, please insert values between 1 and 4")