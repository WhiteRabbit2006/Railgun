import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Physical constants
capacitance = 1800  # farads
esr = 0.00016  # ohms
mu0 = (4 * np.pi) * (10 ** -7)  # magnetic constant newtons per ampere squared
copper_resistivity = 1.68 * (10 ** -8)  # ohms * meters
copper_density = 8950  # kg per cubic meter

# Parameters
D = 0.00635  # separation of the rails and width of the bar (meters)
w = 0.00635  # width of the rails (meters)
h = 0.00635  # height of rail
lp = 0.00635  # length of projectile
L = 1  # length of the rails (meters)
mass = lp * w * h * copper_density  # kilogram

initial_voltage = 3.0  # volts
initial_velocity = 1  # meters per second

projectile_resistance = copper_resistivity * D / (D * w)
connection_resistance = 0
e = np.e  # eulers number


def total_resistance(position):
    r_rail = copper_resistivity * (position / (w * h))
    resistance = (2 * r_rail) + esr + projectile_resistance + connection_resistance
    return resistance


def dydt(y, t):
    position, velocity, voltage, current = y

    resistance = total_resistance(position)
    voltage_derivative = -current / capacitance
    current_derivative = voltage_derivative / resistance
    force_rails = ((mu0 * position * (current ** 2)) / np.pi) * (
                (w ** 2) + (position ** 2) - (position * np.sqrt((w ** 2) + (position ** 2)))) / (
                      np.sqrt(w ** 2 + (position ** 2)))
    # force_friction = -velocity * area * coefficient_of_dynamic_friction

    if position < L:
        acceleration = force_rails / mass
    else:
        acceleration = 0

    y = [velocity, acceleration, voltage_derivative, current_derivative]
    return y


time = np.linspace(0, 1, 101)

initial_current = (initial_voltage / total_resistance(0))
y0 = [0.0, initial_velocity, initial_voltage, initial_current]
y1 = odeint(dydt, y0, time)

# plt.plot(time, y1[:, 0], 'b', label='position')
plt.plot(time, y1[:, 1], 'g', label='velocity')
# plt.plot(time, y1[:, 2], 'g', label='voltage')
# plt.plot(time, y1[:, 3], 'g', label='current')
plt.legend(loc='best')
plt.xlabel('t')
plt.show()
plt.grid()
