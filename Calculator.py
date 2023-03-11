import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Capacitor
capacitance = 3400  # [=] farads
esr = 0.00013  # [=] ohms
initial_voltage = 2.85  # [=] volts

# Railgun, projectile, leads, construction
w = 0.00635 * 2  # width of the rails [=] meters
h = 0.00635 * 2  # height of rail [=] meters
l = 1  # length of the rails [=] meters
d = 0.00635  # separation of the rails and width of the bar [=] meters

lp = 0.00635  # length of projectile [=] meters
hp = 0.00635  # height of projectile [=] meters

lc = 0.3048 / 3  # length of conductor (both leads added) [=] meters
dc = 0.018288  # diameter or connector wire [=] meters

angle = 0  # angle of launch (from ground) [=] degrees
initial_velocity = 1  # meters per second

# Physical constants (should not need to change)
mu0 = (4 * np.pi) * (10 ** -7)  # magnetic constant [=] newtons per ampere squared
muR = 0.999994  # relative permeability of copper, very close to permeability of free space (mu0) [=] mu / mu0
copper_resistivity = 1.68 * (10 ** -8)  # [=] ohms * meters
aluminum_resistivity = 2.65 * (10 ** -8)  # [=] ohms * meters
cross_c = np.pi * (dc / 2) ** 2  # [=] meters squared
connection_resistance = copper_resistivity * lc / cross_c  # [=] ohms
aluminum_resistance = aluminum_resistivity * lc / cross_c  # [=] ohms
copper_density = 8950 * 1000  # [=] g per cubic meter
friction_coefficient = 0.2  # friction coefficient of copper [=]
# mass = lp * d * hp * copper_density  # [=] grams
mass = 0.69  # 0.25 in. aluminum cube [=] grams
# inductance_gradient = ((4 * mu0 * muR * (d + w)) / h)  # inductance constant (multiplies with position for inductance)
inductance_gradient = 3 * (10 ** -7)  # measured value [=] henris / meters
resistance_gradient = 2 * aluminum_resistivity / (w * h)  # resistance coefficient of rails [=] ohms / position
projectile_resistance = 1.63934426 * aluminum_resistivity * d / (
        hp * lp)  # coefficient for aluminum projectile [=] ohms
weight = mass * 9.81
friction_force = friction_coefficient * weight * np.cos(np.radians(angle))  # [=] newtons
static_resistance = esr + projectile_resistance + connection_resistance  # [=] ohms
# inductance_leads = (mu0 * muR * lc * (d + w)) / h  # [=] volts / dI/dt
inductance_leads = 4.46 * (10 ** -7)  # measured value [=] henris
capacitor_energy = 1 / 2 * capacitance * initial_voltage ** 2  # [=] joules
initial_current_rate = initial_voltage / inductance_leads


def dydt(y, t):
    position, velocity, current, current_rate, voltage = y

    acceleration = (1 / 2 * inductance_gradient * current ** 2 - friction_force) / mass
    d_capacitor_voltage = 1 / capacitance * current
    d_resistance = static_resistance * current_rate
    d_resistance_gradient = resistance_gradient * (current * velocity + position * current_rate)
    d_EMF = inductance_gradient * (current * acceleration + velocity * current_rate)
    # d_inductance = inductance_leads * current_rate_rate
    # d_inductance_gradient = inductance_gradient * (position * current_rate_rate + current_rate * velocity)

    # 'd_inductance' and 'd_inductance_gradient' contain the term 'current_rate_rate' which must be solved for, and
    # they are included in the final equation via the addition of 'inductance1' and multiplication of 'inductance2'
    inductance1 = inductance_gradient * current_rate * velocity
    inductance2 = -1 / (inductance_leads + inductance_gradient * position)

    current_rate_rate = inductance2 * (d_capacitor_voltage + d_resistance + d_resistance_gradient + d_EMF + inductance1)

    if position > l:
        acceleration = 0
        current_rate_rate = 0
        current_rate = 0
        current = 0

    return velocity, acceleration, current_rate, current_rate_rate, -current / capacitance


time = np.linspace(0, 0.22, 101)
y0 = [0.0, initial_velocity, 0, initial_current_rate, initial_voltage]
y1 = odeint(dydt, y0, time)

final_velocity = (y1[:, 1][-1])
projectile_energy = 1 / 2 * mass * (final_velocity - initial_velocity) ** 2
energy_efficiency = projectile_energy / capacitor_energy
print('final_velocity =', round(final_velocity, 4), 'm/s')
print('energy_efficiency =', round(energy_efficiency * 1000, 4), "%")

# plt.plot(time, y1[:, 0], 'b', label='position')
plt.plot(time, y1[:, 1], 'r', label='velocity')
# plt.plot(time, y1[:, 2], 'g', label='current')
# plt.plot(time, y1[:, 3], 'o' label='current_rate')
plt.plot(time, y1[:, 4], 'k', label='voltage')
plt.legend(loc='best')
plt.xlabel('t')
plt.show()
plt.grid()
