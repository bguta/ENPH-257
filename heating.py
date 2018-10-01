import matplotlib.pyplot as plt
import numpy as np
import math

L = 0.5  # length of rod
inital = 20 + 273.15
Dx = 0.01  # steps of x
N = int(L / Dx)  # the number of steps
Dt = 0.1  # steps in time

k = 45.0  # W / m / K
c = 490.0  # J / kg / K
p = 2712.0  # kg / m^3


r = 0.005  # m radius
Pin = 2  # W  power in
T_amb = 20 + 273.15  # ambient temp
k_c = 11.09  # 1.32 * (70 - 20) ** (1 / 4) / (2 * r)  # W/m^2/K convection constant
epi = 0.01
sigma = 5.67 * 10 ** (-8)  # W/m^2/K (stefan-Boltzmann constant)

C = k / (c * p)  # constant
x = np.linspace(0.0, L, N)

T = [inital] + [20.0 + 273.15] * (N - 1)
T = np.array(T)

t = 0
print(str(x))
q = 50000
plt.figure(figsize=(20, 20))
while True:
    if t > q:
        print("time: " + str(t))
    # print(T[1])
    for j in range(1000):
        for i in range(1, N - 1):
            #k_c = 1.32 * (T[i] - T_amb) ** (1 / 4) / (2 * r) ** (1/4)

            T[i] = T[i] + C * Dt * (T[i - 1] - 2 * T[i] + T[i + 1]) / (Dx ** 2)  # heat transfer (conduction)

            T[i] = T[i] - 2 * Dt * k_c * (T[i] - T_amb) / (c * p * r)  # convection loss

            T[i] = T[i] - 2 * Dt * epi * sigma * (T[i] ** (4) - T_amb ** (4)) / (c * p * r)  # radiation loss

        T[0] = T[0] - C * Dt * (T[0] - T[1]) / (Dx ** 2)  # heat transfer (conduction)loss

        T[0] = T[0] + Pin * Dt / (c * p * math.pi * r ** (2) * Dx)  # % power gain
        #k_c = 1.32 * (T[0] - T_amb) ** (1 / 4) / (2 * r) ** (1/4)

        T[0] = T[0] - Dt * (2 / (c * r * p) + 1 / (c * p * Dx)) * k_c * (T[0] - T_amb)  # convection loss

        T[0] = T[0] - Dt * (2 / (c * r * p) + 1 / (c * p * Dx)) * epi * sigma * (T[0]**(4) - T_amb ** 4)  # radiation loss
        #k_c = 1.32 * (T[-1] - T_amb) ** (1 / 4) / (2 * r) ** (1/4)

        T[-1] = T[-1] + C * Dt * (T[-2] - T[-1]) / (Dx**2)  # heat gain (conduction)
        T[-1] = T[-1] - Dt * (2 / (c * r * p) + 1 / (c * p * Dx)) * k_c * (T[-1] - T_amb)  # convection loss

        T[-1] = T[-1] - Dt * (2 / (c * r * p) + 1 / (c * p * Dx)) * epi * sigma * (T[-1]**(4) - T_amb ** 4)  # radiation loss

        t += Dt

    if t >= 10000:

        plt.title(
            "Steady state Temperature versus position on a Basketball",
            fontsize=30)
        plt.xlabel("Distance from start (m)", fontsize=18)
        plt.ylabel("Temperature (deg C)", fontsize=18)
        plt.grid('on')

        plt.axis([0, L, 0, 100])
        plt.plot(
            x, T - 273.15 * np.ones(T.shape), "b-", label="Left End: {}, Right End: {}".format(T[0] - 273.15, T[-1] - 273.15))
        plt.legend(loc=0, prop={'size': 15})
        plt.pause(0.01)
        plt.clf()
        # input()
