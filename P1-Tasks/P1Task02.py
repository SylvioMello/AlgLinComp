import math
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

def residue(lambda0, lambda1):
    return abs(lambda1 - lambda0) / abs(lambda1)
 
def Power_Method(matrix_A, tol):
    n = len(matrix_A)
    x0 = [1 for _ in range(n)]
    lambda0, storage = 1, []
    while True:
        x1 = [0 for _ in range(n)]
        for i in range(n):
            for j in range(n):
                x1[i] += matrix_A[i][j] * x0[j]
        lambda1 = max(x1)
        for i in range(n):
            x1[i] = x1[i] / lambda1
        error = residue(lambda0, lambda1)
        storage.append(error)
        if error < tol:
            num_iter = len(storage)
            return lambda1, x1, num_iter
        x0 = x1[:]
        lambda0 = lambda1

def symmetrical(matrix_A):
    n = len(matrix_A)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if matrix_A[i][j] != matrix_A[j][i]:
                return False
    return True

def check_smaller_than_tol(matrix_A, n, tol):
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if abs(matrix_A[i][j]) < tol:
                continue
            else:
                return False
    return True

def Jacobi_Method(matrix_A, tol):
    if not symmetrical(matrix_A):
        sys.exit("Matrix must be symmetrical")
    n = len(matrix_A)
    num_iter = 1
    X0 = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    while True:
        bigger, bigger_position = 0, [0, 0]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if bigger < abs(matrix_A[i][j]):
                    bigger = abs(matrix_A[i][j])
                    bigger_position = [i, j]
        if matrix_A[bigger_position[0]][bigger_position[0]] == matrix_A[bigger_position[1]][bigger_position[1]]:
            phi = math.pi / 4
        else:
            arctg_num = (2 * bigger)
            arctg_deno = matrix_A[bigger_position[0]][bigger_position[0]] - matrix_A[bigger_position[1]][bigger_position[1]]
            phi = (1 / 2) * math.atan(arctg_num / arctg_deno)
        sen = math.sin(phi)
        cos = math.cos(phi)
        for j in range(n):
            x = matrix_A[bigger_position[0]][j]
            y = matrix_A[bigger_position[1]][j]
            matrix_A[bigger_position[0]][j] = cos * x + sen * y
            matrix_A[bigger_position[1]][j] = -sen * x + cos * y
        for i in range(n):
            x = matrix_A[i][bigger_position[0]]
            y = matrix_A[i][bigger_position[1]]
            matrix_A[i][bigger_position[0]] = cos * x + sen * y
            matrix_A[i][bigger_position[1]] = -sen * x + cos * y
            x = X0[i][bigger_position[0]]
            y = X0[i][bigger_position[1]]
            X0[i][bigger_position[0]] = cos * x + sen * y
            X0[i][bigger_position[1]] = -sen * x + cos * y
        if check_smaller_than_tol(matrix_A, n, tol):
            return matrix_A, X0, num_iter
        num_iter += 1

def Jacobi_Determinant(matrix_A):
    n = len(matrix_A)
    det = 1
    for i in range(n):
        for j in range(n):
            if i == j:
                det *= matrix_A[i][j]
    return det

matrix_A   = read_matrix_A('Matriz_A.dat')
ICOD = int(sys.argv[1])
tol  = float(sys.argv[2])

if ICOD == 1:
    print("You chose the Power Method")
    print("This solution method does not allow the calculation of the determinant")
    print("The highest Eignvalue, its Eignvector and the number of iterations are: ")
    print(Power_Method(matrix_A, tol))
elif ICOD == 2:
    print("You chose the Jacobi Method")
    flag_det = str(input("Would you like to calculate the determinant? "))
    print("The Matrix of Eignvalues, the Eignvectors and the numbers of iterations: ")
    print(Jacobi_Method(matrix_A, tol))
    if "yes" in flag_det or "sim" in flag_det:
        print(f"Determinant: {Jacobi_Determinant(matrix_A)}")
else:
    sys.exit("Wrong ICOD, please insert values between 1 and 2")