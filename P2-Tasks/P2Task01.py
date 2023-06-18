import numpy as np
def newton_method(initial_vector, tol, max_iter):
    vector = initial_vector
    for _ in range(max_iter):
        J = jacobian(vector)
        F = func(vector)
        delta_x = -np.linalg.solve(J, F)
        vector += delta_x
        tol_k = np.linalg.norm(delta_x) / np.linalg.norm(vector)
        if tol_k < tol:
            print(f"Convergance acquired with {_} iterations")
            return vector
    print("Convergence not reached")
    return None

def broyden_method(initial_vector, initial_jacobian, tol, max_iterations):
    current_vector = np.array(initial_vector)
    current_jacobian = np.array(initial_jacobian)
    for _ in range(max_iterations):
        F = np.array(func(current_vector))
        deltaX = -np.linalg.solve(current_jacobian, F)
        next_vector = current_vector + deltaX
        tolk = np.linalg.norm(deltaX) / np.linalg.norm(next_vector)
        if tolk < tol:
            print(f"Convergance acquired with {_} iterations")
            return next_vector
        Y = func(next_vector) - F
        deltaXT = np.transpose(deltaX)
        current_jacobian = current_jacobian + (np.outer(Y - np.dot(current_jacobian, deltaX), deltaXT)) / np.dot(deltaXT, deltaX)
        current_vector = next_vector
    print("Convergence not reached")
    return None

def func(vector):
    x, y, z = vector
    return [
        16 * (x ** 4) + 16 * (y ** 4) + (z ** 4) - 16,
        (x ** 2) + (y ** 2) + (z ** 2) - 3,
        (x ** 3) - y + z - 1
    ]

def jacobian(vector):
    x, y, z = vector
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

# Set the tolerance, maximum number of iterations and the initial vector
tolerance = 1e-6
max_iterations = 10000
initial_vector = [1, 1, 1]

# Apply Newton's method
solution_newton = newton_method(initial_vector, tolerance, max_iterations)
print("Solution using Newton's method:", solution_newton)

# Set the tolerance, maximum number of iterations and the initial vector
tolerance = 1e-6
max_iterations = 10000
initial_vector = [1, 1, 1]

# Apply Broyden's method
solution_broyden = broyden_method(initial_vector, jacobian(initial_vector), tolerance, max_iterations)
print("Solution using Broyden's method:", solution_broyden)