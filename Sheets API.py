from __future__ import print_function

import os.path

from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from enum import Enum
from pprint import pprint

# Before using program, download the Credentials file from the google project credentials page for your project, and
# save it to the working directory with name 'credentials.json'

# If modifying these scopes, delete the file token.json.
SCOPES = ['https://www.googleapis.com/auth/spreadsheets']

# The ID and range of the Google sheet.
spreadsheet_id = '13CeT2uwaPntn6KWsTfxRYzwRJTfuDqsKefQ4t2nhP48'
range_name = 'D2:D16'

copper_resistivity = 1.68 * (10 ** -8)  # [=] ohms * meters
aluminum_resistivity = 2.65 * (10 ** -8)  # [=] ohms * meters
copper_density = 8950 * 1000  # [=] g per cubic meter
aluminum_density = 2710 * 1000  # [=] g per cubic meter
inductance_gradient = 3 * (10 ** -7)  # measured value [=] henris / meters
inductance_leads = 4.46 * (10 ** -7)  # measured value [=] henris  # constants


def calculations(vals):  # vals = [cap, esr, vol, w, h, l, d, lp, hp, mat, lc, dc, angle, vel, fric]
    # Capacitor
    capacitance = float(vals[0])  # [=] farads
    esr = float(vals[1]) / 1000  # [=] ohms
    initial_voltage = float(vals[2])  # [=] volts

    # Railgun, projectile, leads, construction
    w = float(vals[3]) / 1000  # width of the rails [=] meters
    h = float(vals[4]) / 1000  # height of rail [=] meters
    l = float(vals[5])  # length of the rails [=] meters
    d = float(vals[6]) / 1000  # separation of the rails and width of the bar [=] meters

    lp = float(vals[7]) / 1000  # length of projectile [=] meters
    hp = float(vals[8]) / 1000  # height of projectile [=] meters
    if vals[9] == "Aluminum":
        p_material = -1
    elif vals[9] == "Copper":
        p_material = 0
    else:
        p_material = -1

    lc = float(vals[10])  # length of conductor (both leads added) [=] meters
    dc = float(vals[11]) / 1000  # diameter or connector wire [=] meters

    angle = float(vals[12])  # angle of launch (from ground) [=] degrees
    initial_velocity = 0.44704 * float(vals[14])  # meters per second

    # Physical constants (should not need to change)
    cross_c = np.pi * (dc / 2) ** 2  # cross-section of lead wire [=] meters squared
    connection_resistance = copper_resistivity * lc / cross_c  # [=] ohms
    friction_coefficient = float(vals[13])  # friction coefficient of copper [=] newtons / newtons (no units)
    if p_material == 0:
        mass = lp * d * hp * copper_density  # [=] grams
        projectile_resistance = copper_resistivity * d / (hp * lp)  # resistance of copper projectile [=] ohms
    elif p_material == -1:
        mass = lp * d * hp * aluminum_density  # [=] grams
        projectile_resistance = aluminum_resistivity * d / (hp * lp)  # resistance aluminum projectile [=] ohms
    else:
        mass = "specify here"  # for other material [=] grams
        projectile_resistance = "specify"  # for other material [=] ohms
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

        current_rate_rate = inductance2 * (inductance1 + d_voltage + d_resistance + d_resistance_gradient + d_EMF)

        if position > l:
            acceleration = 0
            current_rate_rate = 0
            current_rate = 0
            current = 0

        return velocity, acceleration, current_rate, current_rate_rate, -current / capacitance

    length = 3
    time = np.linspace(0, length, 999996)
    y0 = [0.0, initial_velocity, 0, initial_current_rate, initial_voltage]
    y1 = odeint(dydt, y0, time)

    pos = closest_value(y1[:, 0], l) / len(y1[:, 0])  # used for graphing, to set correct x-axis window
    interval = int(1.2 * 999996 * pos)
    final_velocity = (y1[:, 1][-1])
    projectile_energy = 1 / 2 * mass * (final_velocity - initial_velocity) ** 2
    capacitor_energy_used = 1 / 2 * capacitance * (initial_voltage - y1[:, 4][-1]) ** 2  # charge lost by capacitor
    energy_efficiency = projectile_energy / capacitor_energy_used * 100  # percentage of energy transferred to projectile
    if final_velocity <= 0:
        energy_efficiency = 0

    # final_velocity = (y1[:, 1][-1])
    # if final_velocity <= 0:
    #     projectile_energy = 0
    # else:
    #     projectile_energy = 1 / 2 * mass * (final_velocity - initial_velocity) ** 2
    # capacitor_energy_used = 1 / 2 * capacitance * (initial_voltage - y1[:, 4][-1]) ** 2  # charge lost by capacitor
    # energy_efficiency = projectile_energy / capacitor_energy_used * 100  # percentage of energy transferred
    # pos = list(y1[:, 1]).index(final_velocity)  # finds when projectile stops accelerating
    # interval = pos + int(0.1 * pos)  # defines time interval of data, 10% greater than window of acceleration
    # if interval >= (100000 - 3):  # prevents from exceeding google sheets maximum cell count
    #     interval = 999999 - 3

    velocity_data = [[i] for i in list(y1[:, 1])]  # convert each item in list to list within list
    voltage_data = [[i] for i in list(y1[:, 4])]  # convert each item in list to list within list

    return [velocity_data, voltage_data, energy_efficiency, interval]


def main():

    creds = None
    # The file token.json stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json', SCOPES)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('token.json', 'w') as token:
            token.write(creds.to_json())

    try:
        service = build('sheets', 'v4', credentials=creds)

        # Call the Sheets API
        sheet = service.spreadsheets()
        result = sheet.values().get(spreadsheetId=spreadsheet_id,
                                    range=range_name).execute()
        values = result.get('values', [])

        if not values:
            print('No data found.')
            return

        param = [i for [i] in values]
        calc = calculations(param)

        print(sheet.values().clear(spreadsheetId=spreadsheet_id, range='F3:F999999').execute())
        print(sheet.values().update(
            spreadsheetId=spreadsheet_id,
            range='D18:D19',
            valueInputOption="USER_ENTERED",
            body={"values": [[calc[0][calc[3]][0]], [calc[2]]]}
        ).execute())
        print(sheet.values().update(
            spreadsheetId=spreadsheet_id,
            range='F3:F' + str(calc[3] + 3),
            valueInputOption="USER_ENTERED",
            body={"values": calc[0][0:(calc[3])]}
        ).execute())

    except HttpError as err:
        print(err)


if __name__ == '__main__':
    main()

# TODO: implement double run of dydt() in order to determine proper time interval, so more datapoints can be obtained
#  for relevant interval efficiently (run once to determine interval, run again with interval)
