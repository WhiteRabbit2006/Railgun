import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Capacitor
capacitance = 3400  # [=] farads
esr = 0.00015  # [=] ohms
initial_voltage = 2.85  # [=] volts

# Railgun, projectile, leads, construction
w = 0.00635 * 2  # width of the rails [=] meters
h = 0.00635 * 2  # height of rail [=] meters
l = 1  # length of the rails [=] meters
d = 0.00635  # separation of the rails and width of the bar [=] meters

lp = 0.00635  # length of projectile [=] meters
hp = 0.00635  # height of projectile [=] meters
material = -1  # 0 for copper, -1 for aluminum, else specify 'mass' and 'projectile_resistance' values below

lc = 0.3048 / 3  # length of conductor (both leads added) [=] meters
dc = 0.018288  # diameter or connector wire [=] meters

angle = 0  # angle of launch (from ground) [=] degrees
initial_velocity = 1  # meters per second

# Physical constants (should not need to change)
mu0 = (4 * np.pi) * (10 ** -7)  # magnetic constant [=] newtons per ampere squared
copper_resistivity = 1.68 * (10 ** -8)  # [=] ohms * meters
aluminum_resistivity = 2.65 * (10 ** -8)  # [=] ohms * meters
cross_c = np.pi * (dc / 2) ** 2  # [=] meters squared
connection_resistance = copper_resistivity * lc / cross_c  # [=] ohms
copper_density = 8950 * 1000  # [=] g per cubic meter
aluminum_density = 2710 * 1000  # [=] g per cubic meter
friction_coefficient = 0.2  # friction coefficient of copper [=]
if material == 0:
    mass = lp * d * hp * copper_density  # [=] grams
    projectile_resistance = copper_resistivity * d / (hp * lp)  # resistance of copper projectile [=] ohms
elif material == -1:
    mass = lp * d * hp * aluminum_density  # [=] grams
    projectile_resistance = aluminum_resistivity * d / (hp * lp)  # resistance aluminum projectile [=] ohms
else:
    mass = "specify here"  # for other material [=] grams
    projectile_resistance = "specify"  # for other material [=] ohms
inductance_gradient = 3 * (10 ** -7)  # measured value [=] henris / meters
inductance_leads = 4.46 * (10 ** -7)  # measured value [=] henris
static_resistance = esr + projectile_resistance + connection_resistance  # [=] ohms
resistance_gradient = 2 * aluminum_resistivity / (w * h)  # resistance coefficient of rails [=] ohms / position
weight = mass * 9.81  # [=] newtons
friction_force = friction_coefficient * weight * np.cos(np.radians(angle))  # [=] newtons
capacitor_energy = 1 / 2 * capacitance * initial_voltage ** 2  # [=] joules
initial_current_rate = initial_voltage / inductance_leads


def dydt(y, t):
    position, velocity, current, current_rate, voltage = y

    if position < l:
        acceleration = (1 / 2 * inductance_gradient * current ** 2 - friction_force) / mass
        current_rate_rate = -1 * (
                1 / capacitance * current + static_resistance * current_rate + resistance_gradient * (
                current * velocity + position * current_rate) + inductance_gradient * current_rate * velocity +
                inductance_gradient * (current * acceleration + velocity * current_rate)) / (
                                    inductance_leads + inductance_gradient * position)
    else:
        acceleration = 0
        current_rate_rate = 0
        current_rate = 0
        current = 0

    return velocity, acceleration, current_rate, current_rate_rate, -current / capacitance


time = np.linspace(0, 0.5, 101)
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
