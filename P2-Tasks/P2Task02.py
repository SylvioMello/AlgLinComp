import numpy as np

def function(x):
    """Define your function here"""
    return np.sin(x)

def adaptive_integration(a, b, tol):
    # Initialize variables
    num_iterations = 0
    num_points = 3
    intervals = np.array([[a, b]])
    integral_values = []

    while True:
        num_iterations += 1
        integral = 0.0
        new_intervals = []

        for interval in intervals:
            # Divide the interval into two sub-intervals
            sub_interval_points = np.linspace(interval[0], interval[1], num_points)

            # Calculate function values at the midpoints
            midpoints = (sub_interval_points[:-1] + sub_interval_points[1:]) / 2
            f_midpoints = function(midpoints)

            # Calculate the integral in each sub-interval using Simpson's rule
            integral_sub = 0.0
            for i in range(0, num_points - 1, 2):
                h = sub_interval_points[i + 2] - sub_interval_points[i]
                integral_sub += (h / 6) * (function(sub_interval_points[i]) + 4 * f_midpoints[i // 2] + function(sub_interval_points[i + 2]))

            integral += integral_sub

            # Check if the sub-intervals need further subdivision
            if num_iterations > 1:
                prev_integral = integral_values[-1]
                error_ratio = np.abs(integral_sub - prev_integral) / np.abs(prev_integral)

                if error_ratio > tol:
                    new_intervals.append([sub_interval_points[i], sub_interval_points[i + 2]])

        integral_values.append(integral)

        # Check for convergence
        if num_iterations > 1:
            prev_integral = integral_values[-2]
            error_ratio = np.abs(integral - prev_integral) / np.abs(prev_integral)

            if error_ratio <= tol:
                break

        # Update the intervals for the next iteration
        intervals = new_intervals
        num_points *= 2

    print("Number of iterations:", num_iterations)
    print("Number of points of integration:", num_points)
    print("Integral values on each step:", integral_values)

    return integral

# Test the adaptive integration function
a = 0
b = np.pi
tolerance = 1e-6
result = adaptive_integration(a, b, tolerance)
print("Final result:", result)