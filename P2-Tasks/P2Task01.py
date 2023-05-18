import sys

def Jacobian(vector_X):
    c2, c3, c4 = vector_X

def setting_non_linear_func(vector_X, c1, c2):
    c2, c3, c4 = vector_X

def norma(vector_X):
    result = 0


def LU_Decomposition(matrix_A, vector_B):
    n = len(matrix_A)
    # Turning A into L and U
    for i in range(n):
        for j in range(i + 1, n):
            if matrix_A[i][i] == 0:
                sys.exit("Null pivot detected, LU decomposition without pivoting is not possible")
            matrix_A[j][i] = matrix_A[j][i] / matrix_A[i][i]
        for k in range(i + 1, n):
            for l in range(i + 1, n):
                matrix_A[l][k] -= (matrix_A[l][i] * matrix_A[i][k])
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

def Newton_Method(n1, n2, tol=0.00001, num_iter=1000):
    x0 = [1, 0, 0]
    x1 = [1, 0, 0]
    for i in range(num_iter):
        Jacob = Jacobian(x0)
        funct = setting_non_linear_func(x0, n1, n2)
        delta = LU_Decomposition(Jacob, funct)
        for j in range(len(delta)):
            delta[j] = -delta[j]
        for k in range(len(x1)):
            x1[k] = x0[k] + delta[k]
        tol_test = norma(delta) / norma(x1)
        if tol_test < tol:
            return x1
        x0 = x1[:]
    return "Converge acquired"
