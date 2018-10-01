import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

b_animate = True


R = 8.3144598  # universal gas constant J/mol/K
g = 9.81  # m/s^2
m = 28.9647  # g/mol  mean molar mass of air
T_mean = 288  # mean surface temp (K)
gamma = 1.4
dz = 0.1  # km

# p = m * P / (R * T_mean)


gamma_d = -m * g * (1 - 1 / gamma) / (R)

P_0 = 101 * 10 ** 3  # Pa
T_0 = 15 + 273.15  # K
C = T_0 * np.power(P_0, 1 / gamma - 1)


T = lambda z: T_0 + gamma_d * z
P = lambda z: np.power(C / T(z), 1 / (1 / gamma - 1))
P_z = lambda z: P_0 * np.exp(-m * g * z / (R * T_mean))
P_z2 = lambda z: P_0 * np.exp(-m * g * z / (R * T(z)))


# plt.figure(figsize=(20, 20))
z = np.arange(0, 20, dz)


P_ = [P_0]
T_ = [T_0]
for i in range(200):
    P_.append(P_[i] * (1 - m * g * dz / (R * T_[i])))
    T_.append(T_[i] - dz * m * g * (1 - 1 / gamma) / (R))


dP_dz = lambda z: m * P(z) * g / (R * T(z))
rho = lambda z: m * P(z) / (R * T(z))
drho_dz = lambda z: m * dP_dz(z) / (R * T(z)) - m * P_z2(z) * gamma / (R * T(z)**2)

# plt.figure(figsize=(20, 20))
# plt.title("Pressure as a function of Altitude", fontsize=20)
# plt.xlabel("Altitude (km)", fontsize=18)
# plt.ylabel("Pressure (Pa)", fontsize=18)
# plt.xticks(z)
# plt.yticks(range(0, P_0 + 1, 1000))
# # plt.tick_params([0.0001, 0.0002, ..., 20])
# plt.grid('on')
# plt.plot(z, [P_0 / np.exp(1)] * z.size, 'r--', label="P*T")
# plt.plot(z, P(z), 'k--', label="Exp")
# plt.legend('best', prop={'size': 15})
# plt.axis([5, 10, 30000, 50000])


plt.figure(figsize=(20, 20))
plt.title("Pressure as a function of Altitude", fontsize=20)
plt.xlabel("Altitude (km)", fontsize=18)
plt.ylabel("Pressure (Pa)", fontsize=18)
# plt.xticks(z)
#plt.yticks(range(0, P_0 + 1, 1000))
# plt.tick_params([0.0001, 0.0002, ..., 20])
plt.grid('on')
l1, = plt.plot(z, [P_0 / np.exp(1)] * z.size, 'r--', label="1/e", linewidth=5)
l2, = plt.plot(z, P_[:-1], 'k-', label="Pressure v Altitude", linewidth=10)
plt.legend('best', prop={'size': 15}, handles=[l1, l2])
#plt.axis([5, 10, 30000, 50000])

'''
plt.figure(figsize=(20, 20))

plt.title(
    "Rate of change of Pressure as a function of Altitude",
    fontsize=30)
plt.xlabel("Altitude (km)", fontsize=18)
plt.ylabel("dPressure (Pa/km)", fontsize=18)
plt.grid('on')
plt.plot(z, dP_dz(z), "r--")

plt.figure(figsize=(20, 20))

plt.title(
    "Density as a function of Altitude",
    fontsize=30)
plt.xlabel("Altitude (km)", fontsize=18)
plt.ylabel("Density (g/$m^3$)", fontsize=18)
plt.grid('on')
plt.plot(z, rho(z), "b--")

plt.figure(figsize=(20, 20))

plt.title(
    "Rate of change of Density as a function of Altitude",
    fontsize=30)
plt.xlabel("Altitude (km)", fontsize=18)
plt.ylabel("dDensity (g/$m^3$/km)", fontsize=18)
plt.grid('on')
plt.plot(z, drho_dz(z), "b-")
plt.plot(z, -6.79645849755056 * z + 140, "k-")
'''

plt.show()


# dP_dz = lambda z: -m * P(z) * g / (R * (T_mean + gamma_d))
