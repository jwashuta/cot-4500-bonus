import numpy as np

a = np.array([[3, 1, 1], [1, 4, 1], [2, 3, 7]])
b = np.array([1, 3, 0])
tolerance = 1e-6
iterations = 50

# 1
def func(x, x_0, tolerance):
    return (max(abs(x - x_0))) / (max(abs(x_0)) + tolerance)

def gauss_seidel(a, b, tolerance, iterations):
    length = len(b)
    x = np.zeros((length), dtype=np.double)
    k = 1
    while (k <= iterations):
        x_0 = x.copy()
        for i in range(length):
            first_sum = second_sum = 0
            for j in range(i):
                first_sum += (a[i][j] * x[j])
            for j in range(i + 1, length):
                second_sum += (a[i][j] * (x_0[j]))
            x[i] = (1 / a[i][i]) * (-first_sum - second_sum + b[i])
            if (func(x, x_0, tolerance) < tolerance):
                return k
        k += 1
    return k

print(gauss_seidel(a, b, tolerance, iterations), "\n")

# 2
def jacobi(a, b, tolerance, iterations):
    length = len(b)
    x = np.zeros((length), dtype=np.double)
    k = 1
    while (k <= iterations):
        x_0 = x.copy()
        for i in range(length):
            sum = 0
            for j in range(length):
                if j != i:
                    sum += (a[i][j] * x_0[j])
            x[i] = (1 / a[i][i]) * (-sum + b[i])
            if (func(x, x_0, tolerance) < tolerance):
                return k
        k += 1
    return k

print(jacobi(a, b, tolerance, iterations), "\n")

# 3
initial_approx: float = 0.5
tol: float = .000001
sequence: str = "x**3 - (x**2) + 2"

def deriv(value):
    return (3 * value * value) - (2 * value)

def newton_raphson(initial_approx: float, tol: float, sequence: str):
    count = 0
    x = initial_approx
    f = eval(sequence)
    f_prime = deriv(initial_approx)
    approx: float = f / f_prime

    while(abs(approx) >= tol):
        x = initial_approx
        f = eval(sequence)
        f_prime = deriv(initial_approx)
        approx = f / f_prime
        initial_approx -= approx
        count += 1
    return count

print(newton_raphson(initial_approx, tol, sequence), "\n")

# 4
def apply_div_diff(matrix: np.array):
    # use a for loop to go through the matrix
    for i in range(2, len(matrix)):
        for j in range(2, i + 2):
            
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue
            numerator = matrix[i][j - 1] - matrix[i - 1][j - 1]
            denominator = matrix[i][0] - matrix[i - j + 1][0]
            operation = numerator / denominator           
            matrix[i][j] = operation 
    return matrix
  
def hermite_interpolation():
    x_points = [0, 1, 2]
    y_points = [1, 2, 4]
    slopes = [1.06, 1.23, 1.55]
    size = len(x_points) 
    matrix = np.zeros((size * 2, size * 2))
    index = 0
   
    for x in range(0, size * 2, 2):      
        matrix[x][0] = x_points[index]
        matrix[x + 1][0] = x_points[index]
        index += 1
    index = 0
  
    for y in range(0, size * 2, 2):
        matrix[y][1] = y_points[index]
        matrix[y + 1][1] = y_points[index]
        index += 1
    index = 0
   
    for i in range(1, size * 2, 2):
        matrix[i][2] = slopes[index]
        index += 1
    filled_matrix = apply_div_diff(matrix)
    print(filled_matrix, "\n")

hermite_interpolation()

# 5
in_point = 0.5
p_a = 0
p_b = 3
n = 100  

def function_given(t, y):
    return y - (t**3)

def modified_eulers(in_point, p_a, p_b, n):
    h = (p_b - p_a) / n
    t = p_a
    w = in_point   
    for i in range(n):
        w = w + ((h / 2) * (function_given(t, w) + 
                            function_given(t + h, w + (h * function_given(t, w)))))
        t += h    
    return w

print("%.5f" %modified_eulers(in_point, p_a, p_b, n))
