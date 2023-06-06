import numpy as np
def metodo_newton(tol=0.0001, n_iter=10000):
    x0 = [1, 0, 0]
    x1 = [1, 0, 0]
    for _ in range(n_iter):
        j = jacobiano(x0)
        f = funcao(x0)
        delta_x = decomposicao_lu(j, f, 3)
        print(delta_x)
        # Multiplicando o vetor delta_x por -1
        for i in range(len(delta_x)):
            delta_x[i] = -delta_x[i]
        for i in range(len(x1)):
            x1[i] = x0[i] + delta_x[i]
        tol1 = norma(delta_x) / norma(x1)
        if tol1 < tol:
            return x1
        x0 = x1[:]
    return "Convergência não obtida. Número máximo de iterações alcançada."


def metodo_broyden(t1, t2, tol=0.0001, n_iter=10000):
    x0 = [1, 0, 0]
    x1 = [0, 0, 0]
    b = jacobiano(x0)
    for _ in range(n_iter):
        j = b
        f0 = funcao(x0, t1, t2)
        delta_x = decomposicao_lu(j, f0, 3)
        for i in range(len(delta_x)):
            delta_x[i] = -delta_x[i]
        for i in range(len(x1)):
            x1[i] = x0[i] + delta_x[i]
        f1 = funcao(x1, t1, t2)
        y = subt_vetores(f1, f0)
        tol1 = norma(delta_x) / norma(x1)
        if tol1 < tol:
            return x1
        # Multiplicação entre matriz B e vetor delta_x
        bAx = mult_matriz_vetor(b, delta_x)
        # Subtração do Y - B * delta_x
        y_bAx = subt_vetores(y, bAx)
        # Achar o valor correspondente da multiplicação de delta_x transposto por delta_x
        AxTAx = 0
        for i in range(len(delta_x)):
            AxTAx += delta_x[i] * delta_x[i]
        # Multiplicar o vetor Y - B * delta_x com delta_x transposto,
        # Além disso, dividimos por AxTAx e também já adicionamos com a matriz B,
        # atualizando a mesma
        for i in range(3):
            for j in range(3):
                b[i][j] += y_bAx[i] * delta_x[j] / AxTAx
        x0 = x1[:]
    return "Convergência não obtida. Número máximo de iterações alcançada."


def funcao(vetor_x):
    """Função que retorna os valores das funções dadas as constantes"""
    x, y, z = vetor_x
    f1 = 16 * (x ** 4) + 16 * (y ** 4) + (z ** 4) - 16
    f2 = (x ** 2) + (y ** 2) + (z ** 2) - 3
    f3 = (x ** 3) - y + z - 1
    return [f1, f2, f3]


def jacobiano(vetor_x):
    """Função que retorna a matriz Jacobiana dada as constantes"""
    x, y, z = vetor_x
    J = []
    j11, j12, j13 = 64 * (x ** 3), 64 * (y ** 3), 4  * (z ** 3)
    line = [j11, j12, j13]
    J.append(line)
    j21, j22, j23 = 2 * x, 2 * y, 2 * z
    line = [j21, j22, j23]
    J.append(line)
    j31, j32, j33 = 3 * (x ** 2), -1, 1
    line = [j31, j32, j33]
    J.append(line)
    return J


def mult_matriz_vetor(matriz, vetor):
    """Dado uma matriz e um vetor, fazemos sua multiplicação"""
    dim = len(matriz)
    resultado = [0 for _ in range(dim)]

    for i in range(dim):
        for j in range(dim):
            resultado[i] += matriz[i][j] * vetor[j]

    return resultado


def norma(vetor):
    """Retorna a norma de um vetor"""
    resultado = 0

    for elem in vetor:
        resultado += elem ** 2

    return resultado ** (1/2)


def add_vetores(v1, v2):
    """Retorna o vetor resultante da soma de dois vetores"""
    dim = len(v1)
    v = []
    for i in range(dim):
        v.append(v1[i] + v2[i])
    return v


def subt_vetores(v1, v2):
    """Retorna o vetor resultante da subtração de dois vetores"""
    dim = len(v1)
    v = []
    for i in range(dim):
        v.append(v1[i] - v2[i])
    return v


def decomposicao_lu(matriz_orig, b_orig, n):
    # Modificação da matriz A em LU, sem pivotamento
    matriz = matriz_orig[:]
    b = b_orig[:]
    for k in range(n):
        for i in range(k + 1, n):
            if matriz[k][k] == 0:
                max_row_index = i + np.argmax(np.abs(matriz[i:, i]))
                matriz[[i, max_row_index]] = matriz[[max_row_index, i]]
                b[[i, max_row_index]] = b[[max_row_index, i]]
            matriz[i][k] = matriz[i][k] / matriz[k][k]
        for j in range(k + 1, n):
            for i in range(k + 1, n):
                matriz[i][j] -= (matriz[i][k] * matriz[k][j])
    # Modificamos o vetor B para se transformar em Y
    for i in range(1, n):
        for j in range(i):
            b[i] -= (matriz[i][j] * b[j])
    # Modificar o vetor B para se transforar em X
    b[n - 1] = b[n - 1] / matriz[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, n):
            b[i] -= (matriz[i][j] * b[j])
        b[i] = b[i] / matriz[i][i]
    return b

print(metodo_newton())