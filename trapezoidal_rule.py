import matplotlib.pyplot as plt
import numpy as np

#define the function f(x)
f = lambda t, x: np.cos(x*t**2)

def trapezoidal_rule(f, lower_bound, upper_bound, x, n):
    """
    Approximates the definite integral of f(x) from a to b using the trapezoidal rule
    with n equally spaced panels.

    Arguments:
    f -- the integrand (function)
    lower_bound -- the lower limit of integration
    upper_bound -- the upper limit of integration
    x -- constant variable
    n -- the number of panels (1/delta_t)

    Returns:
    The approximate value of the definite integral
    """
    n = int(n)
    if n <= 0:
        print("Error: Number of panels must be positive.")
        return None

    # Step size
    delta_t = (upper_bound - lower_bound) / n

    # Initialize sum
    integral = 0

    # Sum over all panels
    for i in range(n):
        t0 = lower_bound + i * delta_t
        t1 = t0 + delta_t
        integral += (f(t0, x) + f(t1, x)) * delta_t / 2

    return integral


def taylor_series_approximation(t, x):
    """
    Approximation of the indefinat integral of cos(x*t**2) using a taylor expanion at x=0

    Arguments:
    t -- function argument value
    x -- constant variable

    Returns:
    The approximated function value of the indefinate integral
    """
    return  - (t**13 * x**6) / 9360 + (t**9 * x**4) / 216 - (t**5 * x**2) / 10 + t


def richardson_rule(f, a, b, x, tol, k_max):
    """
    Error estimation using Richardson extrapolation for a second-order convergence scheme such
    as the trapezoidal rule.
    
    Arguments:
    f -- the integrand (function)
    a -- the lower limit of integration
    b -- the upper limit of integration
    x -- constant variable
    tol -- error tolerance
    k_max -- max number of iterations
    
    Returns:
    trap -- value of f(x)
    n -- n that gives f(x) within the given error tolerance tol
    k -- the number of iterations
    booleans - True/Fales whether the given tolerance was archived or not
    """
    n = 10
    error = np.inf
    k = 0
    while error > tol and k < k_max:
        n *= 2
        trap = trapezoidal_rule(f, a, b, x, n)
        rich_h = (4/3)*(trapezoidal_rule(f, a, b, x, n//2)-trap)
        rich_h_2 = (4/3)*(trapezoidal_rule(f, a, b, x, n//4)-trapezoidal_rule(f, a, b, x, n//2))
        error = np.abs(rich_h-rich_h_2) / (2**2 -1)
        k += 1
    
    if k == k_max:
        return trap, n, k, False
    else:
        return trap, n, k, True


def estimate_f(x, tol):
    """
    Using Richardson error estimation to find the value of f(x) for a given errror tolerance 

    Arguments:
    x -- constant variable
    tol -- error tolerance

    Returns:
    res -- vaule of f(x)
    n -- n that gives f(x) within the given error tolerance tol
    k -- number of iterations
    booleans - True/Fales weather the given tolerance was archived ore not
    """

    a = 0; b = 1; k_max = 20
    res, n, k, success = richardson_rule(f, a, b, x, tol, k_max)
    if success:
        return res, n, k, True
    else:
        return res, n, k, False

#part b)
a = 0; b = 1; n = 1/0.01
x_vals = np.linspace(0,4, num=100)
f_val_trapezoidal = []
f_val_taylor = []
for x in x_vals:
    f_val_trapezoidal.append( trapezoidal_rule(f, a, b, x, n) )
    f_val_taylor.append( taylor_series_approximation(b,x)-taylor_series_approximation(a,x) ) 
plt.plot(x_vals, f_val_trapezoidal, label="Trapezoidal Rule")
plt.plot(x_vals, f_val_taylor, label="Taylor Series Approximation")
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Comparison of Trapezoidal Rule and Taylor Series Approximation")
plt.show()   

#part c)
x = 1; a = 0; b = 1
n_vals = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
trapezoidal_rule_vals = [trapezoidal_rule(f, a, b, x, n) for n in n_vals]
error = np.abs(trapezoidal_rule_vals - trapezoidal_rule_vals[-1])
plt.plot(n_vals, error, 'o-', label="absolut error")
plt.plot(n_vals, ((b-a)/np.array(n_vals))**2, ":", label="step size dt^2")
plt.xscale('log', base=10)
plt.yscale('log', base=10)
plt.xlabel("Number of panels (log scale)")
plt.ylabel("Absolute error (log scale)")
plt.title("Convergence Study for Trapezoidal Rule")
plt.legend()
plt.show()

#part d)
x = 1; tol = 0.000000001
res, n, k, convergance = estimate_f(x, tol)
if convergance == True:
    print(f"f(x) = {res:.8f} for the given tolerance of {tol}. Number of optimal points {n} found within {k} Iterations.")
else:
    print(f"Requested accuracy not met. Try to increase maxiumum number of iteration or change input.")

 
