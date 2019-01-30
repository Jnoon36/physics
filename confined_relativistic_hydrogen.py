import numpy as np
from scipy import optimize
from mpmath import *  # mpmath is better at evaluating special functions than scipy
from math import sqrt
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')


rest = 510998.9461  # rest energy of the electron
a = 1 / 137.03599  # fine structure constant
b = 2589.605192  # rest energy/hbar*c
bohr = .0529177  # bohr radius in nm


def f(E):
    """
    R is the confining radius in units of the bohr radius, k = +/-(1, 2, 3, ...), E is the energy in units of mc^2, Z is the amount of nucleons (will be equal to one in this case), L = (1/2)*sqrt((2*k+1)**2), kv is the modified bessel function of the second kind, V is the potential energy term (in units of mc^2), s = sqrt(k**2 - (Za)**2)
    """
    A = 1 + E - V
    B = 1 - E + V
    return ((((1 - sqrt((2 * k + 1)**2) + 2 * k) / (2 * R * b * A)) - (sqrt(B / A)) * (besselk(L - 1, sqrt(A * B) * b * R) / besselk(L, sqrt(A * B) * b * R)) - sqrt((1 - E) / (1 + E))) * (((E * Z * a) / sqrt(1 - E**2)) - s) * hyp1f1(s - (E * Z * a / sqrt(1 - E**2)) + 1, (2 * s) + 1, 2 * b * sqrt(1 - E**2) * R)) + ((((1 - sqrt((2 * k + 1)**2) + 2 * k) / (2 * R * b * A)) - (sqrt(B / A)) * (besselk(L - 1, sqrt(A * B) * b * R) / besselk(L, sqrt(A * B) * b * R)) + sqrt((1 - E) / (1 + E))) * (k - (Z * a / sqrt(1 - E**2))) * hyp1f1(s - (E * Z * a / sqrt(1 - E**2)), (2 * s) + 1, 2 * b * sqrt(1 - E**2) * R))


list1 = []  # 1s1/2 radii, energies (ground state), k = -1
list2 = []  # 2s1/2 radii, energies, k = -1
list3 = []  # 3s1/2 radii, energies, k = -1
list4 = []  # 4s1/2 radii, energies, k = -1
list5 = []  # 2p1/2 radii, energies, k = +1
list6 = []  # 3p1/2 radii, energies, k = +1
list7 = []  # 4p1/2 radii, energies, k = +1
list8 = []  # 2p3/2 radii, energies, k = -2
list9 = []  # 3p3/2 radii, energies, k = -2
list10 = []  # 4p3/2 radii, energies, k = -2

"""
Some states you will see two for loops. This is because if we were to stretch out R futher, the function f no longer has opposite signs at the endpoints. So we need to move the endpoint b in closer to endpoint a (the unconfined binding energy) in order to find more roots and assure asymptotic behavior as R gets further away. In simpler english, moving R further out requires a slight change in the interval we hope to find roots.
"""

# k=-1 states

# 1s1/2
for i in np.arange(1.84, 5.44, 0.3):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999733756789446, 0.9999997337419158)  # interval in which we expect a solution, between ground state binding energy and mc2, in units of mc2
    bind = root * rest - rest  # solving for binding energy, which is negative
    list1.append((i, bind))

# 2s1/2
for i in np.arange(6.25, 12.75, 0.5):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999933424520228, 0.9999997337419158)  # interval in which we expect a solution, between ground state binding energy and mc2, in units of mc2
    bind = root * rest - rest  # solving for binding energy, which is negative
    list2.append((i, bind))

for i in np.arange(13.25, 16.25, 0.5):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999933424520228, 0.9999981485284711)  # interval in which we expect a solution, between ground state binding energy and mc2, in units of mc2
    bind = root * rest - rest  # solving for binding energy, which is negative
    list2.append((i, bind))

# 3s1/2
for i in np.arange(13.4, 22.5, 0.7):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999970410897879, 0.9999997337419158)
    bind = root * rest - rest
    list3.append((i, bind))

for i in np.arange(23, 25.5, 0.5):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999970410897879, 0.9999981485284711)
    bind = root * rest - rest
    list3.append((i, bind))

# 4s1/2
for i in np.arange(23.5, 34.3, 0.9):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999983365914814, 0.9999997337419158)
    bind = root * rest - rest
    list4.append((i, bind))

# k = +1 states

# 2p1/2
for i in np.arange(5.2, 11.7, 0.5):
    Z = 1
    k = 1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999933424520228, 0.9999997337419158)
    bind = root * rest - rest
    list5.append((i, bind))

for i in np.arange(12.2, 15.2, 0.5):
    Z = 1
    k = 1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999933424520228, 0.9999981485284711)
    bind = root * rest - rest
    list5.append((i, bind))

# 3p1/2
for i in np.arange(12.4, 21.5, 0.7):
    Z = 1
    k = 1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999970410897879, 0.9999997337419158)
    bind = root * rest - rest
    list6.append((i, bind))

for i in np.arange(22, 25, 0.5):
    Z = 1
    k = 1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999970410897879, 0.9999981485284711)
    bind = root * rest - rest
    list6.append((i, bind))

# 4p1/2
for i in np.arange(23, 34.7, 0.9):
    Z = 1
    k = 1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999983365914814, 0.9999997337419158)
    bind = root * rest - rest
    list7.append((i, bind))

# k = -2 states

# 2p3/2
for i in np.arange(5.2, 11.7, 0.5):
    Z = 1
    k = -2
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999933424520228, 0.9999997337419158)
    bind = root * rest - rest
    list8.append((i, bind))

for i in np.arange(12.2, 15.2, 0.5):
    Z = 1
    k = -2
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999933424520228, 0.9999981485284711)
    bind = root * rest - rest
    list8.append((i, bind))

# 3p1/2
for i in np.arange(12.5, 21.6, 0.7):
    Z = 1
    k = -2
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999970410897879, 0.9999997337419158)
    bind = root * rest - rest
    list9.append((i, bind))

for i in np.arange(22.1, 25.1, 0.5):
    Z = 1
    k = -2
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999970410897879, 0.9999981485284711)
    bind = root * rest - rest
    list9.append((i, bind))

# 4p3/2
for i in np.arange(23, 34.7, 0.9):
    Z = 1
    k = -2
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = i * bohr
    root = optimize.bisect(f, 0.9999983365914814, 0.9999997337419158)
    bind = root * rest - rest
    list10.append((i, bind))

"""
These domains for R values were a result of trial and error. They mostly coincide with the principle quantum numbers (2s1/2 has the same ranges as the 3p1/2 and 3p3/2), but they can be slightly different (s orbitals have different shapes than p orbitals). Generally speaking, the higher the n value for ns1/2, np1/2, or np3/2, the further out R needs to be.
"""

# Next we need to fit polynomials (we will try for 7th degree) through each state. Before that, we must unzip each list of tuples into a list of separate R values (multiples of the bohr radius) and E (bind) values to be used for fitting.

R1, E1 = zip(*list1)
R2, E2 = zip(*list2)
R3, E3 = zip(*list3)
R4, E4 = zip(*list4)
R5, E5 = zip(*list5)
R6, E6 = zip(*list6)
R7, E7 = zip(*list7)
R8, E8 = zip(*list8)
R9, E9 = zip(*list9)
R10, E10 = zip(*list10)

# Next is deriving the coefficients for the polynomials.

coeff1 = np.polyfit(R1, E1, 7)  # 1s1/2 coefficients
coeff2 = np.polyfit(R2, E2, 7)  # 2s1/2 coefficients
coeff3 = np.polyfit(R3, E3, 7)  # 3s1/2 coefficients
coeff4 = np.polyfit(R4, E4, 7)  # 4s1/2 coefficients
coeff5 = np.polyfit(R5, E5, 7)  # 2p1/2 coefficients
coeff6 = np.polyfit(R6, E6, 7)  # 3p1/2 coefficients
coeff7 = np.polyfit(R7, E7, 7)  # 4p1/2 coefficients
coeff8 = np.polyfit(R8, E8, 7)  # 2p3/2 coefficients
coeff9 = np.polyfit(R9, E9, 7)  # 3p3/2 coefficients
coeff10 = np.polyfit(R10, E10, 7)  # 4p3/2 coefficients

# Then we need to generate the 7th degree polynomials

poly1 = np.poly1d(coeff1)  # 1s1/2 polynomial
poly2 = np.poly1d(coeff2)  # 2s1/2 polynomial
poly3 = np.poly1d(coeff3)  # 3s1/2 polynomial
poly4 = np.poly1d(coeff4)  # 4s1/2 polynomial
poly5 = np.poly1d(coeff5)  # 2p1/2 polynomial
poly6 = np.poly1d(coeff6)  # 3p1/2 polynomial
poly7 = np.poly1d(coeff7)  # 4p1/2 polynomial
poly8 = np.poly1d(coeff8)  # 2p3/2 polynomial
poly9 = np.poly1d(coeff9)  # 3p3/2 polynomial
poly10 = np.poly1d(coeff10)  # 4p3/2 polynomial

# Now we can make plots with the polynomials and the discrete data points

# ground state
plt.figure(1)
plt.plot(np.arange(R1[0], R1[len(R1) - 1], .01), poly1(np.arange(R1[0], R1[len(R1) - 1], .01)), color='blue')
plt.scatter(R1, E1, color='blue', marker='.')
plt.xlabel('Confining Radius (Multiples of Bohr Radius)')
plt.ylabel('Binding Energy (eV)')
plt.title('Ground State Binding Energy vs Confinement Radius')

# k=-1 states (except the ground state)
plt.figure(2)
plt.plot(np.arange(R2[0], R2[len(R2) - 1], .01), poly2(np.arange(R2[0], R2[len(R2) - 1], .01)), color='red', label='2s1/2')  # 2s1/2 polynomial
plt.scatter(R2, E2, marker='.', color='red')  # 2s1/2
plt.plot(np.arange(R3[0], R3[len(R3) - 1], .01), poly3(np.arange(R3[0], R3[len(R3) - 1], .01)), color='green', label='3s1/2')  # 3s1/2 polynomial
plt.scatter(R3, E3, marker='.', color='green')  # 3s1/2 data points
plt.plot(np.arange(R4[0], R4[len(R4) - 1], .01), poly4(np.arange(R4[0], R4[len(R4) - 1], .01)), color='blue', label='4s1/2')  # 4s1/2 polynomial
plt.scatter(R4, E4, marker='.', color='blue')  # 4s1/2 data points
plt.xlabel('Confining Radius (Multiples of Bohr Radius)')
plt.ylabel('Binding Energy (eV)')
plt.title('Binding Energy vs Confinement (ns1/2 States)')
plt.legend()

# k=+1 states
plt.figure(3)
plt.plot(np.arange(R5[0], R5[len(R5) - 1], .01), poly5(np.arange(R5[0], R5[len(R5) - 1], .01)), color='red', label='2p1/2')  # 2p1/2 polynomial
plt.scatter(R5, E5, marker='.', color='red')  # 2p1/2
plt.plot(np.arange(R6[0], R6[len(R6) - 1], .01), poly6(np.arange(R6[0], R6[len(R6) - 1], .01)), color='green', label='3p1/2')  # 3p1/2 polynomial
plt.scatter(R6, E6, marker='.', color='green')  # 3p1/2 data points
plt.plot(np.arange(R7[0], R7[len(R7) - 1], .01), poly7(np.arange(R7[0], R7[len(R7) - 1], .01)), color='blue', label='4p1/2')  # 4p1/2 polynomial
plt.scatter(R7, E7, marker='.', color='blue')  # 4p1/2 data points
plt.xlabel('Confining Radius (Multiples of Bohr Radius)')
plt.ylabel('Binding Energy (eV)')
plt.title('Binding Energy vs Confinement (np1/2 States)')
plt.legend()

# k=-2 states
plt.figure(4)
plt.plot(np.arange(R8[0], R8[len(R8) - 1], .01), poly8(np.arange(R8[0], R8[len(R8) - 1], .01)), color='red', label='2p3/2')  # 2p3/2 polynomial
plt.scatter(R8, E8, marker='.', color='red')  # 2p3/2
plt.plot(np.arange(R9[0], R9[len(R9) - 1], .01), poly9(np.arange(R9[0], R9[len(R9) - 1], .01)), color='green', label='3p3/2')  # 3p3/2 polynomial
plt.scatter(R9, E9, marker='.', color='green')  # 3p3/2 data points
plt.plot(np.arange(R10[0], R10[len(R10) - 1], .01), poly10(np.arange(R10[0], R10[len(R10) - 1], .01)), color='blue', label='4p3/2')  # 4p3/2 polynomial
plt.scatter(R10, E10, marker='.', color='blue')
plt.xlabel('Confining Radius (Multiples of Bohr Radius)')
plt.ylabel('Binding Energy (eV)')
plt.title('Binding Energy vs Confinement (np3/2 States)')
plt.legend()


"""
Figures 2, 3, and 4 look like each individual plot is the same as in the other figures, but they are not. If you plot all 9 polynomials (with the data points), you will see a slight difference in radial value and energy value. This is hard to see from three separate plots, but you can go ahead and plot all 9 on the same graph to check it out yourself.

If you look up the energies for each state when the atom is unconfined (which is what each of the confined states here approach as R goes to infinity), you will notice, for example, that the 2s1/2 and 2p1/2 have the same binding energy, but are slightly lower in energy than the 2p3/2, because j = 3/2. This is a result of the spin-orbit interaction. Another point to note is that even though 2s1/2 and 3s/12 states approach the exact same energy, one approaches it faster than the other because of the shape differences betwen the p and s orbitals.

The last thing I want to plot is the comparison of ground state energies between the relativistic and non-relativistic confined hydrogen atoms with the exact same barrier height (V = 1). These numbers come from reference #4 in my paper. We will use the same R values as they do.
"""

r = [5.77827, 4.87924, 4.08889, 3.45203, 3.15412, 2.51487, 2.04918]
nonrelE = [-.9998, -.9990, -.9960, -.9881, -.9803, -.9426, -.8734]  # in units of the rydberg energy

# To transfer their units of energy into our, we just multiply each by 13.605
nonrelEbind = []

for i in range(0, 7):
    bind = nonrelE[i] * (13.605)
    nonrelEbind.append(bind)

# Now we generate the relativistic (ground state) binding energies using their R values. We just look for energies in the same range of values as we did above for the ground state, but set R equal to their R.
relE = []

for i in range(0, 7):
    Z = 1
    k = -1
    V = 1
    s = sqrt((k**2) - ((Z * a)**2))
    L = (1 / 2) * sqrt((2 * k + 1)**2)
    R = r[i] * bohr
    root = optimize.bisect(f, 0.9999733739682344, 0.9999997337419158)
    bind = root * rest - rest
    relE.append(bind)

# We will once again do the same as we did above, by fitting polynomials through their data and ours.

coeff11 = np.polyfit(r, relE, 6)  # relativistic coefficients
coeff12 = np.polyfit(r, nonrelEbind, 5)  # non-relativistic coefficients

poly11 = np.poly1d(coeff11)  # relativistic polynomial
poly12 = np.poly1d(coeff12)  # non-relativistic polynomial

# Now we plot both curves
plt.figure(5)
plt.plot(np.arange(r[6], r[0], .01), poly12(np.arange(r[6], r[0], .01)), color='green', label='Non-relativistic')
plt.scatter(r, nonrelEbind, marker='.', color='green')
plt.plot(np.arange(r[6], r[0], .01), poly11(np.arange(r[6], r[0], .01)), color='blue', label='Relativistic')
plt.scatter(r, relE, marker='.', color='blue')
plt.xlabel('Confining Radius (Multiples of Bohr Radius)')
plt.ylabel('Binding Energy (eV)')
plt.title('Binding Energy vs Confinement Radius (Hydrogen Ground State)')
plt.legend()
plt.show()
