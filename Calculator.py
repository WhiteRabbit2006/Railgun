import numpy as np
import scipy.integrate
from scipy.integrate import odeint, quadrature
import matplotlib.pyplot as plt

# Physical constants
capacitance = 1800  # [=] farads
esr = 0.00016  # [=] ohms
mu0 = (4 * np.pi) * (10 ** -7)  # magnetic constant [=] newtons per ampere squared
muR = 0.999994  # relative permeability of copper, very close to permeability of free space (mu0) [=] mu / mu0
copper_resistivity = 1.68 * (10 ** -8)  # [=] ohms * meters
copper_density = 8950  # [=] kg per cubic meter
friction_coefficient = 0.2  # friction coefficient of copper [=]

# Parameters
D = 0.00635  # separation of the rails and width of the bar [=] meters
w = 0.0508  # width of the rails [=] meters
h = 0.003175  # height of rail [=] meters
lp = 0.00635  # length of projectile [=] meters
hp = 0.00635  # height of projectile [=] meters
L = 1  # length of the rails [=] meters
lc = 0.3048  # length of conductor (leads) [=] meters
# hc =  # height of projectile [=] meters
mass = lp * D * hp * copper_density  # [=] kilograms
angle = 0  # angle of launch (from ground) [=] degrees

# mass = lp * w * h * copper_density  # kilogram

initial_voltage = 3.0  # volts
initial_velocity = 1  # meters per second

projectile_resistance = copper_resistivity * D / (hp * lp)
connection_resistance = 0
e = np.e  # eulers number


def friction_force():
    return friction_coefficient * mass * -9.81 * np.cosd(angle)

def magnetic_force(current, position):  # given current and position returns force [=] newtons

    def dFdp(pp):  # parameter is position along projectile width (D), returns force for pp [=] newtons / meters
        part_1 = np.sqrt(position**2 + pp**2)
        part_2 = np.sqrt(position**2 + (D-pp)**2)
        result = 1/part_1 + 1/part_2
        return result

    constant = mu0 * current**2 * position / (4 * np.pi)
    result = scipy.integrate.quadrature(dFdp, 0, D)
    return result * constant

def current_derivative(current):


def total_resistance(position):
    cross_section_area = w * h
    r_rail = copper_resistivity * position / cross_section_area  # position = length of rail
    resistance = (2 * r_rail) + esr + projectile_resistance + connection_resistance
    return resistance


def total_inductance_mtu(position):
    inductance_rails = ((4 * mu0 * muR * (lp + w)) / h) * position
    inductance_leads = (mu0 * muR * lc * (lp + w)) / h
    return inductance_rails + inductance_leads


def total_inductance_plos(position):
    big_equation = 1  #TODO: write out equation
    return big_equation


# Total inductance L(position) = L'(2 * position)
# L' = (2 * mu0 * muR * 2 * (lp + w)) / h
# L(position) = ((4 * mu0 * muR * 2 * (lp + w)) / h) * position

# Inductance leads = (mu0 * muR * lc * (lp + w)) / h

# Current derivative = (1 / (


def dFdp(p):
    return (1 / (((L **2) + (p ** 2)) ** (1/2))) + (1 / (((L ** 2) + (D - p) ** 2) ** (1/2)))



def dydt(y, t):
    position, velocity, voltage_capacitor, current = y

    resistance = total_resistance(position)
    voltage_capacitor_derivative = -current / capacitance
    current_derivative = voltage_capacitor_derivative / resistance
    force_friction = (friction_coefficient * mass * -9.81 * np.cos(angle))



    # force_rails = 0.5 * ((2 * mu0 * muR * 2 * (lp + w)) / h) * (current ** 2)

    # force_rails = ((mu0 * position * (current ** 2)) / np.pi) * (
    #             (w ** 2) + (position ** 2) - (position * np.sqrt((w ** 2) + (position ** 2)))) / (
    #                   np.sqrt(w ** 2 + (position ** 2)))

    if position < L:
        acceleration = force_rails + force_friction / mass
    else:
        acceleration = 0

    y = [velocity, acceleration, voltage_capacitor_derivative, current_derivative]
    return y


time = np.linspace(0, 0.2, 101)

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
