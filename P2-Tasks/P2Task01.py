import numpy as np

def norm(vector):
    result = 0
    for elem in vector:
        result += elem ** 2
    return result ** (1/2)

def jacobian(X):
    x, y, z = X
    J = np.zeros((3, 3))
    J[0, 0] = 64 * (x ** 3)
    J[0, 1] = 64 * (y ** 3)
    J[0, 2] = 4 * (z ** 3)

    J[1, 0] = 2 * x
    J[1, 1] = 2 * y
    J[1, 2] = 2 * z

    J[2, 0] = 3 * (x ** 2)
    J[2, 1] = -1
    J[2, 2] = 1
    return J

def function(X):
    x, y, z = X
    F = np.zeros(3)
    F[0] = 16 * (x ** 4) + 16 * (y ** 4) + (z ** 4) - 16
    F[1] = (x ** 2) + (y ** 2) + (z ** 2) - 3
    F[2] = (x ** 3) - y + z - 1
    return F

def newton_method(X0, tol, NITER):
    Xk = X0
    for k in range(1, NITER + 1):
        J = jacobian(Xk)
        F = function(Xk)
        try:
            deltaX = -np.linalg.solve(J, F)
        except np.linalg.LinAlgError:
            deltaX = -np.linalg.lstsq(J, F, rcond=None)[0]
        Xk = Xk + deltaX
        tolk = norm(deltaX) / norm(X0)
        if tolk < tol:
            print("Convergence reached.")
            return Xk
    print("Convergence not reached.")
    return None

def broyden_method(x0, tol=1e-6, max_iter=100):
    x = x0.copy()
    B = jacobian(x0)
    for k in range(max_iter):
        F = function(x)
        try:
            delta_x = np.linalg.solve(B, -F)
        except np.linalg.LinAlgError:
            delta_x = -np.linalg.lstsq(B, -F, rcond=None)[0]
        x_new = x + delta_x
        tolk = np.linalg.norm(delta_x) / np.linalg.norm(x_new)
        if tolk < tol:
            return x_new
        Y = function(x_new) - F
        delta_XT = delta_x[:, np.newaxis]
        B = B + (Y - np.dot(B, delta_XT) * delta_XT) / np.dot(delta_XT, delta_x)
        x = x_new
    print("Convergence not reached")
    return None
