import sys
import numpy as np
import sympy as sp

def compute_jacobian(func, variables):
    J = []
    for i in range(len(func)):
        row = []
        for j in range(len(variables)):
            derivative = sp.diff(func[i], variables[j])
            row.append(derivative)
        J.append(row)
    return J

def evaluate_functions(func, values, variables):
    return [expr.subs(zip(variables, values)) for expr in func]

def func(x):
    c2, c3, c4 = x
    return [
        16 * (c2 ** 4) + 16 * (c3 ** 4) + (c4 ** 4) - 16,
        (c2 ** 2) + (c3 ** 2) + (c4 ** 2) - 3,
        (c2 ** 3) - c3 + c4 - 1
    ]

def norma(vector_X):
    result = 0

def Newton_Method(tol=0.00001, num_iter=1000):
    x0 = [1, 0, 0]
    x1 = [1, 0, 0]
    Jacob_result = []
    variables = sp.symbols('c2 c3 c4')
    for i in range(num_iter):
        Jacob = compute_jacobian(func(variables), variables)
        for jacob_funct in Jacob:
            results = evaluate_functions(jacob_funct, x0, variables)
            Jacob_result.append(results)
        funct = evaluate_functions(func(variables), x0, variables)
        print(Jacob_result)
        print(funct)
        delta = LU_Decomposition(Jacob_result, funct)
        for j in range(len(delta)):
            delta[j] = -delta[j]
        for k in range(len(x1)):
            x1[k] = x0[k] + delta[k]
        tol_test = norma(delta) / norma(x1)
        if tol_test < tol:
            return x1
        x0 = x1[:]
    return "Converge acquired"

def LU_Decomposition(matrix_A, vector_B):
    matrix_A = np.array(matrix_A)
    vector_B = np.array(vector_B)
    n = len(matrix_A)
    # Turning A into L and U
    for i in range(n):
        # Pivoting
        if matrix_A[i][i] == 0:
            max_row_index = i + np.argmax(np.abs(matrix_A[i:, i]))
            matrix_A[[i, max_row_index]] = matrix_A[[max_row_index, i]]
            vector_B[[i, max_row_index]] = vector_B[[max_row_index, i]]
        for j in range(i + 1, n):
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

print(Newton_Method())