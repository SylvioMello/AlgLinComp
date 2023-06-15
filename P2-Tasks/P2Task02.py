import sys
import numpy as np

def function_m0(x, Sn1_or_Sn2):
    if Sn1_or_Sn2 == "Sn1":
        Sn = 2
        return Sn / (x**4 - 1.99 * x**2 + 1)
    elif Sn1_or_Sn2 == "Sn2":
        numerator_Sn = 36 * np.pi ** 3
        denominator_Sn = 625 * np.exp((16 * np.pi ** 3) / (625 * x ** 4)) * x ** 5
        Sn = numerator_Sn/denominator_Sn
        return Sn / (x**4 - 1.99 * x**2 + 1)
    else:
        sys.exit("Please provide Sn1 or Sn2")

def function_m2(x, Sn1_or_Sn2):
    if Sn1_or_Sn2 == "Sn1":
        Sn = 2
        return (x ** 2) * (Sn / (x**4 - 1.99 * x**2 + 1))
    elif Sn1_or_Sn2 == "Sn2":
        numerator_Sn = 36 * np.pi ** 3
        denominator_Sn = 625 * np.exp((16 * np.pi ** 3) / (625 * x ** 4)) * x ** 5
        Sn = numerator_Sn/denominator_Sn
        return (x ** 2) * (Sn / (x**4 - 1.99 * x**2 + 1))
    else:
        sys.exit("Please provide Sn1 or Sn2")

def adaptive_integration(a, b, tol, m0_or_m2, Sn1_or_Sn2):
    # Initialize variables
    num_iterations = 0
    num_points = 3
    integral_values = []
    while True:
        num_iterations += 1
        integral = 0.0
        h = (b - a) / num_points
        # Calculate the integral using Simpson's rule
        for i in range(num_points):
            x0 = a + i * h
            x1 = x0 + h
            x2 = x1 + h
            if m0_or_m2 == "m0":
                integral += (h / 6) * (function_m0(x0, Sn1_or_Sn2) + 4 * function_m0(x1, Sn1_or_Sn2) + function_m0(x2, Sn1_or_Sn2))
            elif m0_or_m2 == "m2":
                integral += (h / 6) * (function_m2(x0, Sn1_or_Sn2) + 4 * function_m2(x1, Sn1_or_Sn2) + function_m2(x2, Sn1_or_Sn2))
            else:
                sys.exit("Please use m0 or m2 as the input for the function")
        integral_values.append(integral)
        # Check for convergence
        if num_iterations > 1:
            prev_integral = integral_values[-2]
            error_ratio = np.abs(integral - prev_integral) / np.abs(integral)
            if error_ratio <= tol:
                break
        # Subdivide the intervals
        new_num_points = 2 * num_points - 1
        new_h = (b - a) / (new_num_points - 1)
        midpoints = np.linspace(a, b, new_num_points)
        if m0_or_m2 == "m0":
            f_midpoints = function_m0(midpoints, Sn1_or_Sn2)
        elif m0_or_m2 == "m2":
            f_midpoints = function_m2(midpoints, Sn1_or_Sn2)
        else:
            sys.exit("Please use m0 or m2 as the input for the function")
        integral = 0.0
        for i in range(1, new_num_points - 1, 2):
            integral += (new_h / 6) * (f_midpoints[i - 1] + 4 * f_midpoints[i] + f_midpoints[i + 1])
        num_points = new_num_points
    print("Number of iterations:", num_iterations)
    print("Number of points of integration at the end:", num_points)
    print("Integral values on each step:", integral_values)
    return integral

def gauss_quadrature(a, b, n_points, m0_or_m2, Sn1_or_Sn2):
    if n_points > 11:
        sys.exit("This program only executes up to 11 points of integration")
    table_weights = {
        2:  [1.0, 1.0],
        3:  [0.556, 0.889, 0.556],
        4:  [0.348, 0.652, 0.652, 0.348],
        5:  [0.237, 0.479, 0.569, 0.479, 0.237],
        6:  [0.171, 0.361, 0.468, 0.468, 0.361, 0.171],
        7:  [0.129, 0.280, 0.382, 0.418, 0.382, 0.280, 0.129],
        8:  [0.101, 0.222, 0.314, 0.363, 0.363, 0.314, 0.222, 0.101],
        9:  [0.081, 0.181, 0.261, 0.312, 0.330, 0.312, 0.261, 0.181, 0.081],
        10: [0.067, 0.149, 0.219, 0.269, 0.296, 0.296, 0.269, 0.219, 0.149, 0.067],
        11: [0.056, 0.126, 0.186, 0.233, 0.263, 0.273, 0.263, 0.233, 0.186, 0.126, 0.056]
    }
    table_points = {
        2:  [-0.577, 0.577],
        3:  [-0.775, 0.0, 0.775],
        4:  [-0.861, -0.340, 0.340, 0.861],
        5:  [-0.906, -0.538, 0.0, 0.538, 0.906],
        6:  [-0.932, -0.661, -0.239, 0.239, 0.661, 0.932],
        7:  [-0.949, -0.742, -0.406, 0.0, 0.406, 0.742, 0.949],
        8:  [-0.960, -0.797, -0.526, -0.183, 0.183, 0.526, 0.797, 0.960],
        9:  [-0.968, -0.836, -0.613, -0.324, 0.0, 0.324, 0.613, 0.836, 0.968],
        10: [-0.974, -0.865, -0.679, -0.433, -0.149, 0.149, 0.433, 0.679, 0.865, 0.974],
        11: [-0.978, -0.887, -0.730, -0.519, -0.270, 0.000, 0.270, 0.519, 0.730, 0.887, 0.978]
    }
    result = 0
    weights = table_weights[n_points]
    points = table_points[n_points]
    for i in range(n_points):
        z = points[i]
        xi = (a + b + z * (b - a)) / 2
        if m0_or_m2 == "m0":
            result += function_m0(xi, Sn1_or_Sn2) * weights[i]
        elif m0_or_m2 == "m2":
            result += function_m2(xi, Sn1_or_Sn2) * weights[i]
        else:
            sys.exit("Please use m0 or m2 as the input for the function")
    return round(result * ((b - a) / 2), 4)

# m0 and Sn1
print("\n RESULT FOR M0 AND SN = 2 (SN1) \n")
a          = 0
b          = 10
n_points   = 11
m0_or_m2   = "m0"
tolerance  = 1e-6
Sn1_or_Sn2 = "Sn1"
result_adpative_integration = adaptive_integration(a, b, tolerance, m0_or_m2, Sn1_or_Sn2)
result_gauss_quadrature     = gauss_quadrature(a, b, n_points, m0_or_m2, Sn1_or_Sn2)
print("Final result Adaptative Integration:", result_adpative_integration)
print("Final result Gauss Quadrature:", result_gauss_quadrature)

# m0 and Sn2
print("\n RESULT FOR M0 AND SN = Expression (SN2) \n")
a          = 0.01
b          = 10
n_points   = 11
m0_or_m2   = "m0"
tolerance  = 1e-6
Sn1_or_Sn2 = "Sn2"
result_adpative_integration = adaptive_integration(a, b, tolerance, m0_or_m2, Sn1_or_Sn2)
result_gauss_quadrature     = gauss_quadrature(a, b, 11, m0_or_m2, Sn1_or_Sn2)
print("Final result Adaptative Integration:", result_adpative_integration)
print("Final result Gauss Quadrature:", result_gauss_quadrature)

# m2 and Sn1
print("\n RESULT FOR M2 AND SN = 2 (SN1) \n")
a          = 0
b          = 10
n_points   = 11
m0_or_m2   = "m2"
tolerance  = 1e-6
Sn1_or_Sn2 = "Sn1"
result_adpative_integration = adaptive_integration(a, b, tolerance, m0_or_m2, Sn1_or_Sn2)
result_gauss_quadrature     = gauss_quadrature(a, b, 11, m0_or_m2, Sn1_or_Sn2)
print("Final result Adaptative Integration:", result_adpative_integration)
print("Final result Gauss Quadrature:", result_gauss_quadrature)

# m2 and Sn2
print("\n RESULT FOR M2 AND SN = Expression (SN2) \n")
a          = 0.01
b          = 10
n_points   = 11
m0_or_m2   = "m2"
tolerance  = 1e-6
Sn1_or_Sn2 = "Sn2"
result_adpative_integration = adaptive_integration(a, b, tolerance, m0_or_m2, Sn1_or_Sn2)
result_gauss_quadrature     = gauss_quadrature(a, b, 11, m0_or_m2, Sn1_or_Sn2)
print("Final result Adaptative Integration:", result_adpative_integration)
print("Final result Gauss Quadrature:", result_gauss_quadrature)