
def f(v):
    return -g - (Kd * A * v * abs(v))



def taylor_series_method(z0, v0, h, num_iterations):
    t = 0.0
    z = z0
    v = v0
    print("Iteração\tt\tz(t)")
    print("-------------------------------")
    for i in range(num_iterations):
        z_next = z + v * h + f(v) * (h**2)/2
        v_next = v + f(v) * h
        t += h
        z = z_next
        v = v_next
        print(f"{i+1}\t\t{t:.3f}\t{z:.3f}")


def runge_kutta_nyström_method(z0, v0, h, num_iterations):
    t = 0.0
    z = z0
    v = v0
    print("Iteração\tt\tz(t)")
    print("-------------------------------")
    for i in range(num_iterations):
        k1z = v * h
        k1v = f(v) * h
        k2z = (v + k1v / 2) * h
        k2v = f(v + k1v / 2) * h
        k3z = (v + k2v / 2) * h
        k3v = f(v + k2v / 2) * h
        k4z = (v + k3v) * h
        k4v = f(v + k3v) * h
        t += h
        z += (k1z + 2 * k2z + 2 * k3z + k4z) / 6
        v += (k1v + 2 * k2v + 2 * k3v + k4v) / 6
        print(f"{i+1}\t\t{t:.4f}\t{z:.4f}")



# Parâmetros da equação
m = 1.0    # Massa
g = 9.81   # Aceleração da gravidade
Kd = 1.0   # Coeficiente de arrasto
A = 1.0    # Área
z0 = 0.0   # posição inicial  -> z(t)
v0 = 0.0   # velocidade inicial -> z'(t)
h = 0.1   # Tamanho do passo de tempo -> delta t
num_iterations = 200  # Número de iterações

# Opção do método
method = input("Escolha o método de solução (T - Taylor, R - Runge-Kutta): ")
if method.upper() == 'T':
    print("Método de Taylor")
    taylor_series_method(z0, v0, h, num_iterations)
elif method.upper() == 'R':
    print("Método de Runge-Kutta")
    runge_kutta_nyström_method(z0, v0, h, num_iterations)
else:
    print("Método inválido")

