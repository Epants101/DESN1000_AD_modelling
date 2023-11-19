import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import streams

def calculate_dialysis_water_temp(solar_thermal_temp, tank_temperature, vmass_evap, hex_water_bed1, hex_water_bed2, t):
    mdot_t = streams.mdot_tu(t)
    mdot_r = streams.mdot_mp(hex_water_bed1, t)+streams.mdot_lq(hex_water_bed2, t)
    return streams.temp_u(tank_temperature, solar_thermal_temp, vmass_evap, mdot_t, mdot_r)

results = pd.read_csv("./results.csv")
times = results["t"].to_numpy()
solar_thermal_temp = results["temperature solar thermal"].to_numpy()

vcalcuate_temp = np.vectorize(calculate_dialysis_water_temp)
fig = plt.figure(1)
plt.plot(times, vcalcuate_temp(results["temperature solar thermal"].to_numpy(), results["tank temperature"].to_numpy(), results["vapour mass evaporator"].to_numpy(), results["hex water bed one"].to_numpy(), results["hex water bed two"].to_numpy(), times), 'r-')
plt.xlim([0, 10800])
plt.ylim([0, 30])
plt.xlabel("Time (s)")
plt.ylabel("Temperature of Dialysis Input (°C)")
plt.title("The Temperature of the Dialysis Input over Time")
fig.savefig("./dtemp_graph.png")

fig2 = plt.figure(2)
plt.plot(times, solar_thermal_temp, 'r-')
plt.xlim([0,10800])
plt.ylim([0, 100])
plt.xlabel("Time (s)")
plt.ylabel("Temperature of Solar Thermal Panel (°C)")
plt.title("The Temperature of the Solar Thermal Panel over Time")
fig2.savefig("./stemp_graph.png")
