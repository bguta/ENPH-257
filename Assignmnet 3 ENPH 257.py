import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integral


vmin = 0
vmax = 3000  # m/s
nv = 5000
v_ = np.linspace(vmin, vmax, nv)
# dv
Temp = 20 + 273.15  # kelvin
m_He = 4 * 1.67 * 10 ** -27  # kg
m_air = 28.97 * 1.67 * 10 ** -27  # kg
k_b = 1.3806503 * 10 ** -23

f_v = lambda v, m, T: np.sqrt((m / (2 * np.pi * k_b * T))**3) * 4 * np.pi * v**2 * np.exp(-m * v**2 / (2 * k_b * T))
f_he = lambda v: np.sqrt((m_He / (2 * np.pi * k_b * Temp))**3) * 4 * np.pi * v**2 * np.exp(-m_He * v**2 / (2 * k_b * Temp))
f_air = lambda v: np.sqrt((m_air / (2 * np.pi * k_b * Temp))**3) * 4 * np.pi * v**2 * np.exp(-m_air * v**2 / (2 * k_b * Temp))

f_he2 = lambda v: f_he(v) * v**2
f_air2 = lambda v: f_air(v) * v**2


plt.figure(figsize=(10, 10))
plt.xlabel("Velocity $m/s$", fontsize=10)
plt.ylabel("Probabliity $s/m$", fontsize=10)
plt.title("Maxwell-Boltzmann Distribution @ T = {} $^oC$".format(Temp - 273.15), fontsize=15)

plt.plot(v_, f_v(v_, m_air * np.ones(v_.shape), Temp * np.ones(v_.shape)), "r--", label="Air")
plt.plot(v_, f_v(v_, m_He * np.ones(v_.shape), Temp * np.ones(v_.shape)), "b--", label="Helium")

plt.legend(loc=0, prop={'size': 15})

Air2, error = integral.quad(f_air2, 0, 10000)
Helium2, error = integral.quad(f_he2, 0, 10000)


Helium_RMS = np.sqrt(3 * k_b * (20 + 273.15) / m_He)
Air_RMS = np.sqrt(3 * k_b * (20 + 273.15) / m_air)
print("sound Velocity of Air {} @ {} deg C".format(np.sqrt(7 / (5 * 3)) * Air_RMS, 20))
print("sound Velocity of Helium {} @ {} deg C".format(np.sqrt(5 / (3 * 3)) * Helium_RMS, 20))

print(Air_RMS)
print(Helium_RMS)

print(np.sqrt(Air2))
print(np.sqrt(Helium2))


plt.show()
