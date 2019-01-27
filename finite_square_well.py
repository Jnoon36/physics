import numpy as np
from scipy.misc import derivative
from sympy import *
import math
import matplotlib.pyplot as plt


def bisection(f, a, b, tol, N):
    """
    N is the maximum number of iterations, tol is the tolerance (how accurate of a solution you want).
    Generally, as the tolerance becomes smaller, the number of iterations N must grow. a and b are the endpoints of the close interval, while f is your function.
    """
    if f(a) * f(b) > 0:
        return 'No real roots contained within the chosen interval.'
    else:
        n = 1
        while n <= N:
            c = (b + a) / 2
            if abs(f(c)) < tol:
                return c
                break
            elif f(a) * f(c) < 0:
                b = c
                n = n + 1
            else:
                a = c
                n = n + 1
        return 'No root found within chosen tolerance. Max iterations exceeded.'


# There will be two functions to find the roots of; one is the symmetric solution and the other is the anti-symmetric solution.

def symmetric(x):
    """
    x is a measure of energy, and x_0 is a measure of how deep and wide the well is. Larger x_0 corresponds to a deeper and shallower well.
    """
    return np.tan(x) - np.sqrt(((x_0 / x)**2) - 1)


def antisymmetric(x):
    """
    x is a measure of energy, and x_0 is a measure of how deep and wide the well is. Larger x_0 corresponds to a deeper and shallower well.
    """
    return (np.cos(x) / np.sin(x)) - np.sqrt(((x_0 / x)**2) - 1)


# First we can obtain a qualitative picture for what the discrete energy states look like based on where the functions tan(x) and cot(x) both equal sqrt(((x_0/x)**2)-1). We will look at various values of x_0.


def g(x, x_0):
    return np.sqrt(((x_0 / x)**2) - 1)

 # tan has asymptotic bounds at npi/2, where n is odd
domain1 = np.linspace(0, 1.5, 100)
domain2 = np.linspace(1.6, 4.7, 200)
domain3 = np.linspace(4.8, 7.7, 200)

domain4 = np.linspace(0.5, 6, 400)  # domain for g(x, 6)
domain5 = np.linspace(0.2, 4, 400)  # domain for g(x, 4)
domain6 = np.linspace(0.2, 9, 600)  # domain for g(x, 9)

# cot has asymptotic bounds at npi, where n is an integer
domain7 = np.linspace(0.1, 3.1, 200)
domain8 = np.linspace(3.2, 6.2, 200)
domain9 = np.linspace(6.3, 9.2, 200)


# x_0 = 4, symmetric states
plt.figure(1)
plt.plot(domain1, np.tan(domain1), color='blue', label='tan(x)')
plt.plot(domain2, np.tan(domain2), color='blue')
plt.plot(domain5, g(domain5, 4), color='red', label='sqrt((4/x)**2 - 1)')
plt.vlines(np.pi / 2, 0, 10, linestyle='dashed')
plt.vlines(3 * np.pi / 2, 0, 10, linestyle='dashed')
plt.ylim(0, 10)
plt.xlabel('x')
plt.title('Symmetric States (x_0 = 4)\nDiscrete Energies at Intersections')
plt.legend(loc=1)

"""
The intersection of these curves represents the discrete states. With x_0 = 4 in this plot, only one state is allowed.
"""


# x_0 = 6, symmetric states
plt.figure(2)
plt.plot(domain1, np.tan(domain1), color='blue', label='tan(x)')
plt.plot(domain2, np.tan(domain2), color='blue')
plt.plot(domain4, g(domain4, 6), color='red', label='sqrt((6/x)**2 - 1)')
plt.vlines(np.pi / 2, 0, 10, linestyle='dashed')
plt.vlines(3 * np.pi / 2, 0, 10, linestyle='dashed')
plt.ylim(0, 10)
plt.xlabel('x')
plt.title('Symmetric States (x_0 = 6)\nDiscrete Energies at Intersections')
plt.legend()

"""
We can now do examples of an anti-symmetric states.
"""
#x_0 = 6, anti-symmetric
plt.figure(3)
plt.plot(domain7, np.cos(domain7) / np.sin(domain7), color='green', label='cot(x)')
plt.plot(domain8, np.cos(domain8) / np.sin(domain8), color='green')
plt.plot(domain4, g(domain4, 6), color='red', label='sqrt((6/x)**2 - 1)')
plt.vlines(0, 0, 10, linestyle='dashed')
plt.vlines(np.pi, 0, 10, linestyle='dashed')
plt.ylim(0, 10)
plt.xlabel('x')
plt.title('Anti-symmetric States (x_0 = 6)\nDiscrete Energies at Intersections')
plt.legend()


#x_0 = 9, anti-symmetric
plt.figure(4)
plt.plot(domain7, np.cos(domain7) / np.sin(domain7), color='green', label='cot(x)')
plt.plot(domain8, np.cos(domain8) / np.sin(domain8), color='green')
plt.plot(domain9, np.cos(domain9) / np.sin(domain9), color='green')
plt.plot(domain6, g(domain6, 9), color='red', label='sqrt((9/x)**2 - 1)')
plt.vlines(0, 0, 10, linestyle='dashed')
plt.vlines(np.pi, 0, 10, linestyle='dashed')
plt.vlines(2 * np.pi, 0, 10, linestyle='dashed')
plt.ylim(0, 10)
plt.xlabel('x')
plt.title('Anti-symmetric States (x_0 = 9)\nDiscrete Energies at Intersections')
plt.legend()

"""
Now we want to see the relationship between energy (x) and the well depth and width (x_0). This is where our root finding algorithms come in. We will choose a value of x_0 (must be greater than zero, otherwise we don't have bound states) and then use one of the algorithms to find the value of x that best satisfies the symmetric and anti-symmetric functions when they are equal to zero. From the figures that were generated above, we can get a good idea of the interavals that will be used in searching for the zeros (roots).
"""

# symmetric states
symmetric_roots1 = []
symmetric_roots2 = []
symmetric_roots3 = []

# interval 1
for i in np.arange(2, 10, 0.5):
    x_0 = i
    root = bisection(symmetric, 0.25, 1.5, .01, 20)
    symmetric_roots1.append(root)

# interval 2
for i in np.arange(5, 13, 0.5):
    x_0 = i
    root = bisection(symmetric, 3.16, 4.7, .01, 20)
    symmetric_roots2.append(root)

# interval 3
for i in np.arange(8, 17, 0.5):
    x_0 = i
    root = bisection(symmetric, 6.3, 7.8, .01, 20)
    symmetric_roots3.append(root)


# anti-symmetric states
antisymmetric_roots1 = []
antisymmetric_roots2 = []
antisymmetric_roots3 = []

# interval 1
for i in np.arange(5, 13, 0.5):
    x_0 = i
    root = bisection(antisymmetric, 3.16, 4.7, .01, 20)
    antisymmetric_roots1.append(root)

# interval 2
for i in np.arange(8, 17, 0.5):
    x_0 = i
    root = bisection(antisymmetric, 6.3, 7.8, .01, 20)
    antisymmetric_roots2.append(root)

# interval 3
for i in np.arange(11, 20, 0.5):
    x_0 = i
    root = bisection(antisymmetric, 9.5, 10.9, .01, 20)
    antisymmetric_roots3.append(root)


"""
Now it's time to fit a polynomial through the data points using the numpy packages polyfit and poly1d. First we will find the coefficients, then generate a polynomial using those coefficients. Lastly, we we will plot the discrete data points along with the polynomials that fit those points.
"""

# generating the coefficients for the symmetric points, 7th degree polynomial
symmetric_coeff1 = np.polyfit(np.arange(2, 10, 0.5), symmetric_roots1, 7)
symmetric_coeff2 = np.polyfit(np.arange(5, 13, 0.5), symmetric_roots2, 7)
symmetric_coeff3 = np.polyfit(np.arange(8, 17, 0.5), symmetric_roots3, 7)

# generating the coefficients for the anti-symmetric points, 7th degree polynomial
antisymmetric_coeff1 = np.polyfit(np.arange(5, 13, 0.5), antisymmetric_roots1, 7)
antisymmetric_coeff2 = np.polyfit(np.arange(8, 17, 0.5), antisymmetric_roots2, 7)
antisymmetric_coeff3 = np.polyfit(np.arange(11, 20, 0.5), antisymmetric_roots3, 7)

# genrating the polynomial for the symmetric points
symmetric_poly1 = np.poly1d(symmetric_coeff1)
symmetric_poly2 = np.poly1d(symmetric_coeff2)
symmetric_poly3 = np.poly1d(symmetric_coeff3)

# generating the polynomial for the anti-symmetric points
antisymmetric_poly1 = np.poly1d(antisymmetric_coeff1)
antisymmetric_poly2 = np.poly1d(antisymmetric_coeff2)
antisymmetric_poly3 = np.poly1d(antisymmetric_coeff3)


# symmetric plot
plt.figure(5)
plt.plot(np.arange(2, 10, 0.5), symmetric_poly1(np.arange(2, 10, 0.5)), color='blue')
plt.plot(np.arange(5, 13, 0.5), symmetric_poly2(np.arange(5, 13, 0.5)), color='blue')
plt.plot(np.arange(8, 17, 0.5), symmetric_poly3(np.arange(8, 17, 0.5)), color='blue')
plt.scatter(np.arange(2, 10, 0.5), symmetric_roots1, marker='d', label='Symmetric', color='blue')
plt.scatter(np.arange(5, 13, 0.5), symmetric_roots2, marker='d', color='blue')
plt.scatter(np.arange(8, 17, 0.5), symmetric_roots3, marker='d', color='blue')
plt.xlabel('Well Depth/Width (x_0)')
plt.ylabel('Energy (x)')
plt.title('Symmetric Energy States')


# anti-symmetric plot
plt.figure(6)
plt.plot(np.arange(5, 13, 0.5), antisymmetric_poly1(np.arange(5, 13, 0.5)), color='red')
plt.plot(np.arange(8, 17, 0.5), antisymmetric_poly2(np.arange(8, 17, 0.5)), color='red')
plt.plot(np.arange(11, 20, 0.5), antisymmetric_poly3(np.arange(11, 20, 0.5)), color='red')
plt.scatter(np.arange(5, 13, 0.5), antisymmetric_roots1, marker='.', label='Anti-Symmetric', color='red')
plt.scatter(np.arange(8, 17, 0.5), antisymmetric_roots2, marker='.', color='red')
plt.scatter(np.arange(11, 20, 0.5), antisymmetric_roots3, marker='.', color='red')
plt.xlabel('Well Depth/Width (x_0)')
plt.ylabel('Energy (x)')
plt.title('Anti-Symmetric Energy States')
plt.show()


"""
One very import point to notice from these plots: As the well depth/width gets larger, the energy states for the symmetric solution get larger, but the anti-symmetric states get smaller. As the well depth/width gets smaller, the symmetric states get smaller, but the anti-symmetric energy states get larger. This is also evident from figures 1-4; the shapes of the tangent and cotanget functions would have given this fact away even before figures 5 and 6. The sqrt function intersects the tangent function at a larger x value as the x-intercept of the sqrt function moves further away from the origin (i.e., x_0 gets larger). Conversely, the sqrt function intersects the cotangent function at a smaller x value as the x-intercept of the sqrt function moves further from the origin.
"""
