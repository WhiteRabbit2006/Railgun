import numpy as np
import scipy.integrate
from scipy.integrate import odeint, quadrature
import matplotlib.pyplot as plt

# Railgun and projectile construction
D = 0.00635  # separation of the rails and width of the bar [=] meters
w = 0.0508  # width of the rails [=] meters
h = 0.003175  # height of rail [=] meters
L = 1  # length of the rails [=] meters
lp = 0.00635  # length of projectile [=] meters
hp = 0.00635  # height of projectile [=] meters
lc = 0.3048  # length of conductor (both leads added) [=] meters
angle = 0  # angle of launch (from ground) [=] degrees
initial_velocity = 1  # meters per second

# Capacitor
capacitance = 1800  # [=] farads
esr = 0.00016  # [=] ohms
initial_voltage = 2.85  # volts

# Physical constants
mu0 = (4 * np.pi) * (10 ** -7)  # magnetic constant [=] newtons per ampere squared
muR = 0.999994  # relative permeability of copper, very close to permeability of free space (mu0) [=] mu / mu0
copper_resistivity = 1.68 * (10 ** -8)  # [=] ohms * meters
copper_density = 8950  # [=] kg per cubic meter
friction_coefficient = 0.2  # friction coefficient of copper [=]
mass = lp * D * hp * copper_density  # [=] kilograms
inductance_gradient = ((4 * mu0 * muR * (lp + w)) / h)  # inductance constant (multiplies with position for inductance)
resistance_gradient = 2 * copper_resistivity / (w * h)  # resistance of rails (multiplies with position for resistance)
projectile_resistance = copper_resistivity * D / (hp * lp)
connection_resistance = 0
weight = mass * 9.81
friction_force = friction_coefficient * weight * np.cos(np.radians(angle))  # [=] newtons
static_resistance = esr + projectile_resistance + connection_resistance  # [=] ohms
inductance_leads = (mu0 * muR * lc * (lp + w)) / h  # [=] volts / dI/dt


def capacitor_voltage(charge):  # [=] volts
    v = charge / capacitance
    return v


def EMF(velocity):  # returns resistance due to EMF [=] ohms
    return inductance_gradient * velocity


def dydt(y, t):
    position, velocity, current, current_rate = y

    if position < L:
        acceleration = (1 / 2 * inductance_gradient * current ** 2 - friction_force) / mass
        current_rate_rate = -1 * (
                -1 / capacitance * current + static_resistance * current_rate + resistance_gradient * (
                current * velocity + position * current_rate) + inductance_gradient * current_rate * velocity +
                inductance_gradient * (current * acceleration + velocity * current_rate)) / (
                                    inductance_leads + inductance_gradient * position)
    else:
        acceleration = 0
        current_rate_rate = 0
        current_rate = 0

    return velocity, acceleration, current_rate, current_rate_rate


time = np.linspace(0, 0.02, 101)

initial_current_rate = initial_voltage / inductance_leads
y0 = [0.0, initial_velocity, 0, initial_current_rate]
y1 = odeint(dydt, y0, time)

# plt.plot(time, y1[:, 0], 'b', label='position')
# plt.plot(time, y1[:, 1], 'r', label='velocity')
plt.plot(time, y1[:, 2], 'g', label='current')
# plt.plot(time, y1[:, 3], 'o' label='current_rate')
plt.legend(loc='best')
plt.xlabel('t')
plt.show()
plt.grid()
