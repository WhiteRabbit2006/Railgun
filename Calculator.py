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
p_material = -1  # 0 for copper, -1 for aluminum, else specify 'mass' and 'projectile_resistance' values below

lc = 0.3048 / 3  # length of conductor (both leads added) [=] meters
dc = 0.018288  # diameter of connector wire [=] meters

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
friction_coefficient = 1.4  # friction coefficient of copper [=] newtons / newtons (no units)
if p_material == 0:
    mass = lp * d * hp * copper_density  # [=] grams
    projectile_resistance = copper_resistivity * d / (hp * lp)  # resistance of copper projectile [=] ohms
elif p_material == -1:
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
initial_current_rate = initial_voltage / inductance_leads  # derivative of current with respect to time


def closest_value(input_list, input_value):
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    return i


def dydt(y, t):
    position, velocity, current, current_rate, voltage = y

    acceleration = (1 / 2 * inductance_gradient * current ** 2 - friction_force) / mass
    d_voltage = 1 / capacitance * current
    d_resistance = static_resistance * current_rate
    d_resistance_gradient = resistance_gradient * (current * velocity + position * current_rate)
    d_EMF = inductance_gradient * (current * acceleration + velocity * current_rate)
    # d_inductance = inductance_leads * 'current_rate_rate'
    # d_inductance_gradient = inductance_gradient * (position * 'current_rate_rate' + current_rate * velocity)
    # 'd_inductance' and 'd_inductance_gradient' contain the term 'current_rate_rate' which must be solved for, and
    # they are included in the final equation via the addition of 'inductance1' and multiplication of 'inductance2'
    inductance1 = inductance_gradient * current_rate * velocity
    inductance2 = -1 / (inductance_leads + inductance_gradient * position)

    if position < l and velocity > 0:
        current_rate_rate = inductance2 * (inductance1 + d_voltage + d_resistance + d_resistance_gradient + d_EMF)
    else:
        acceleration = 0
        current_rate_rate = 0
        current_rate = 0
        current = 0

    return velocity, acceleration, current_rate, current_rate_rate, -current / capacitance


length = 3
time = np.linspace(0, length, 100000)
y0 = [0.0, initial_velocity, 0, initial_current_rate, initial_voltage]
y1 = odeint(dydt, y0, time)
interval = 1.2 * length * (closest_value(y1[:, 0], l) / len(y1[:, 0]))  # used for graphing, to set correct x-axis window
time = np.linspace(0, interval, 1000)
y1 = odeint(dydt, y0, time)

final_velocity = (y1[:, 1][-1])
projectile_energy = 1 / 2 * mass * (final_velocity - initial_velocity) ** 2
capacitor_energy_used = capacitor_energy - (1 / 2 * capacitance * y1[:, 4][-1] ** 2)  # charge lost by capacitor
energy_efficiency = projectile_energy / capacitor_energy_used * 100  # percentage of energy transferred to projectile
total_energy_efficiency = projectile_energy / capacitor_energy
if final_velocity <= 0:
    energy_efficiency = 0
    projectile_energy = 0
    total_energy_efficiency = 0
print('final_velocity =', round(final_velocity, 4), 'm/s')
print('energy_efficiency =', round(energy_efficiency, 4), "%")
print('total_energy_efficiency =', round(total_energy_efficiency, 4), '%')
print("projectile_energy =", projectile_energy, "joules")

# plt.plot(time, y1[:, 0], 'b', label='position')
plt.plot(time, y1[:, 1], 'r', label='velocity')
# plt.plot(time, y1[:, 2], 'g', label='current')
# plt.plot(time, y1[:, 3], 'o', label='current_rate')
plt.plot(time, y1[:, 4], 'k', label='voltage')
plt.xlim(0, interval)
plt.legend(loc='best')
plt.xlabel('t')
plt.show()
plt.grid()
