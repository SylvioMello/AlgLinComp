import math

def metodo_potencia(matriz, n, tol):
    # ICOD == 1
    x0 = [1 for _ in range(n)]
    lambda0 = 1
    historico = []
    while True:
        # inicializar x1, ou Y (no slide)
        x1 = [0 for _ in range(n)]
        # Fazer a multiplicação Y = A.X0
        for i in range(n):
            for j in range(n):
                x1[i] += matriz[i][j] * x0[j]
        # Normalizar x1
        lambda1 = max(x1)
        for i in range(n):
            x1[i] = x1[i] / lambda1
        erro = residuo(lambda0, lambda1)
        historico.append(erro)
        if erro < tol:
            num_iter = len(historico)
            return x1, lambda1, num_iter
        x0 = x1[:]
        lambda0 = lambda1

def metodo_jacobi(matriz, n, tol):
    # ICOD == 2
    # Matriz A deve ser simétrica
    if not simetrica(matriz):
        return 0, 0, "Execução parada. Matriz deve ser simétrica."
    num_iter = 1
    X0 = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    while True:
        # Achar o maior elemento de A, fora da diagonal principal
        maior = 0
        pos_maior = [0, 0]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if maior < abs(matriz[i][j]):
                    maior = abs(matriz[i][j])
                    pos_maior = [i, j]
        # Achar valor de phi
        if matriz[pos_maior[0]][pos_maior[0]] == matriz[pos_maior[1]][pos_maior[1]]:
            phi = math.pi / 4
        else:
            arctg_num = (2 * maior)
            arctg_deno = matriz[pos_maior[0]][pos_maior[0]] - matriz[pos_maior[1]][pos_maior[1]]
            phi = (1 / 2) * math.atan(arctg_num / arctg_deno)
        # Seno e cosseno de phi
        sen = math.sin(phi)
        cos = math.cos(phi)
        # Fazendo PT.A
        for j in range(n):
            x = matriz[pos_maior[0]][j]
            y = matriz[pos_maior[1]][j]

            matriz[pos_maior[0]][j] = cos * x + sen * y
            matriz[pos_maior[1]][j] = -sen * x + cos * y
        # Fazendo PTA.P e Xk.P
        for i in range(n):
            x = matriz[i][pos_maior[0]]
            y = matriz[i][pos_maior[1]]

            matriz[i][pos_maior[0]] = cos * x + sen * y
            matriz[i][pos_maior[1]] = -sen * x + cos * y

            x = X0[i][pos_maior[0]]
            y = X0[i][pos_maior[1]]

            X0[i][pos_maior[0]] = cos * x + sen * y
            X0[i][pos_maior[1]] = -sen * x + cos * y
        if elem_fora_diag_zero(matriz, n, tol):
            return matriz, X0, num_iter
        num_iter += 1

# Rotinas acessórias
def residuo(lambda0, lambda1):
    return abs(lambda1 - lambda0) / abs(lambda1)

def simetrica(matriz):
    n = len(matriz)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if matriz[i][j] != matriz[j][i]:
                # Não é simétrica
                return False
    return True

def elem_fora_diag_zero(matriz, n, tol):
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if abs(matriz[i][j]) < tol:
                continue
            else:
                return False
    return True

def determinante_jacobi(matriz, n):
    det = 1
    for i in range(n):
        for j in range(n):
            if i == j:
                det *= matriz[i][j]
    return det

A = [[1, 0.2, 0], [0.2, 1, 0.5], [0, 0.5, 1]]
# print(metodo_potencia(A, 3, 0.0001))
print(metodo_potencia(A, 3, 0.001))