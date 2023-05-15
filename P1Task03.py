import sys

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

def Lagrange(points, n, x_expected):
    fx, x, y = 0, [], []
    for point in points:
        x.append(point[0])
        y.append(point[1])
    for i in range(n):
        term = y[i]
        for j in range(n):
            if i == j:
                continue
            try:
                term *= (x_expected - x[j]) / (x[i] - x[j])
            except ZeroDivisionError:
                return "Points have the same coordinate X."
        fx += term
    return fx

def Regression(points, n, x):
    A = [[n, 0], [0, 0]]
    C = [0, 0]
    for i in range(n):
        A[0][1] += points[i][0]
        A[1][0] += points[i][0]
        A[1][1] += (points[i][0]) ** 2
        C[0] += points[i][1]
        C[1] += points[i][0] * points[i][1]
    B = LU_Decomposition(A, C)
    if type(B) == str:
        return B
    return B[0] + x * B[1]

coordinates = []
ICOD        = int(input("Type ICOD Value: "))
N           = int(input("Type the number of pairs of points: "))
for i in range(N):
    xi = float(input(f"Type the {i+1}st value for x: "))
    yi = float(input(f"Type the {i+1}st value for y: "))
    coordinate = [xi, yi]
    coordinates.append(coordinate)
x_expected  = float(input("Type the value of x to which you want to calculate y: "))

if ICOD == 1:
    print("You chose Lagrange")
    print(f"Y = {Lagrange(coordinates, N, x_expected)}")
elif ICOD == 2:
    print("You chose Regression")
    print(f"Y = {Regression(coordinates, N, x_expected)}")