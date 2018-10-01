import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


def getEff(PressureRatio):
    """
    Get the efficiency of  Bryton Cycle using a given Pressure ratio.
    Max T = T3 = 1000 C
    amb T = T1 = 20 C
    """
    P1 = 1.0  # Pa
    V1 = 1.0  # m^3
    P2 = PressureRatio  # Pa
    T1 = 20.0 + 272.15  # Kelvin
    T3 = 1000.0 + 272.15  # Kelvin
    R = 8.3144598  # J * /(mol * K)
    cp = 7.0 * R / 2.0  # specific heat cap (J * /(mol * K))
    gamma = 1.4  # air constant

    V2 = V1 * np.power(P1 / P2, 1 / gamma)

    T2 = T1 * np.power(V1 / V2, gamma - 1)

    V3 = T3 * V2 / T2
    # = V2 + nR(T3 - T2) / P2  # find an eq

    P3 = P2  # Isobar!
    P4 = P1  # Isobar!

    V4 = V3 * np.power(P3 / P4, 1 / gamma)

    T4 = T3 * np.power(V3 / V4, gamma - 1)  # check it again

    T1 = T4 * V1 / V4

    n = P2 * V2 / (R * T2)

    # calc the areas

    # plt1
    a = V2
    b = V1
    f = lambda x: P1 * np.power(V1 / x, gamma)
    plt1Area, error = integrate.quad(f, a, b)

    # plt2
    a = V3
    b = V4
    f = lambda x: P4 * np.power(V4 / x, gamma)
    plt2Area, error = integrate.quad(f, a, b)

    # plt3
    a = V2
    b = V3
    f = lambda x: P2
    plt3Area, error = integrate.quad(f, a, b)

    # plt4
    a = V1
    b = V4
    f = lambda x: P1
    plt4Area, error = integrate.quad(f, a, b)

    Work = plt2Area + plt3Area - plt1Area - plt4Area
    Q_h = n * cp * (T3 - T2)

    eta1 = Work / Q_h

    return eta1

P1 = 1.0  # Pa
V1 = 1.0  # m^3
P2 = 15.0  # Pa
T1 = 20.0 + 272.15  # Kelvin
T3 = 1000.0 + 272.15  # Kelvin
R = 8.3144598  # J * /(mol * K)
cp = 7.0 * R / 2.0  # specific heat cap (J * /(mol * K))
gamma = 1.4  # air constant

V2 = V1 * np.power(P1 / P2, 1 / gamma)

T2 = T1 * np.power(V1 / V2, gamma - 1)

V3 = T3 * V2 / T2
# = V2 + nR(T3 - T2) / P2  # find an eq

P3 = P2  # Isobar!
P4 = P1  # Isobar!

V4 = V3 * np.power(P3 / P4, 1 / gamma)

T4 = T3 * np.power(V3 / V4, gamma - 1)  # check it again

T1 = T4 * V1 / V4

n = P2 * V2 / (R * T2)  # calculate the mols


# plot time!
plt.figure(figsize=(20, 20))
V = np.linspace(V1, V2)
# the curve on the left P1 -> P2 (ADIABATIC)
plt1 = P1 * np.power(V1 / V, gamma)
plt.plot(V, plt1, 'r--', linewidth=1.5, label="ADIABATIC")

V = np.linspace(V4, V3)  # the curve on the right P3 -> P4 (ADIABATIC)
plt2 = P4 * np.power(V4 / V, gamma)
plt.plot(V, plt1, 'b--', linewidth=1.5, label="ADIABATIC")

V = np.linspace(V2, V3)
plt3 = P2 * np.ones(V.shape)  # the line at the top (ISOBAR)
plt.plot(V, plt3, 'g-', linewidth=1.5, label="ISOBAR")

V = np.linspace(V1, V4)
plt4 = P1 * np.ones(V.shape)  # the line at the bottom (ISOBAR)
plt.plot(V, plt4, 'k-', linewidth=1.5, label="ISOBAR")


plt.xlabel('Volume Ratio', fontsize=18)
plt.ylabel('Pressure Ratio', fontsize=18)
plt.title("Pressure as a function of Volueme for a Brayton Cycle", fontsize=30)
plt.grid('on')
plt.legend(loc=0, prop={'size': 15})
# plt.axis([0, 2, 0, 50])


# calc the areas

plt1
a = V2
b = V1
f = lambda x: P1 * np.power(V1 / x, gamma)
plt1Area, error = integrate.quad(f, a, b)

# plt2
a = V3
b = V4
f = lambda x: P4 * np.power(V4 / x, gamma)
plt2Area, error = integrate.quad(f, a, b)

# plt3
a = V2
b = V3
f = lambda x: P2
plt3Area, error = integrate.quad(f, a, b)

# plt4
a = V1
b = V4
f = lambda x: P1
plt4Area, error = integrate.quad(f, a, b)

Work = plt2Area + plt3Area - plt1Area - plt4Area
Q_h = n * cp * (T3 - T2)

eta1 = 100 * Work / Q_h
eta2 = (1 - np.power(P1 / P2, 1 - 1 / gamma)) * 100
print("Eta1: {} , Eta2: {}".format(eta1, eta2))

print("Area: {}".format(Work))


PressureRatio = np.linspace(1, 100)
T_efficiency = (1 - np.power(1 / PressureRatio, 1 - (1 / gamma)))
E_efficiency = [getEff(x) for x in PressureRatio]


plt.figure(figsize=(20, 20))
plt.plot(PressureRatio, T_efficiency, 'b-', linewidth=1.5, label="Theoretical efficiency")
plt.plot(PressureRatio, E_efficiency, 'ro', linewidth=2, label="Calculated Efficiency")


plt.xlabel('Pressure Ratio', fontsize=18)
plt.ylabel('Efficiency', fontsize=18)
plt.title(
    "Thermal efficiency vs. pressure ratio for a Brayton Cycle ",
    fontsize=30)
plt.grid('on')
plt.legend(loc=0, prop={'size': 15})

plt.show()
