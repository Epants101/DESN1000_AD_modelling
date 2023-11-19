import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.linalg import lstsq
import numpy as np

def least_squares(xdat, ydat):
    xdat_matrix = np.hstack(((xdat[:, np.newaxis])**0, (xdat[:, np.newaxis])**1))
    return lstsq(xdat_matrix, ydat)[0]


results = pd.read_csv("./results.csv")
times = results["t"].to_numpy()
waste_output = results["waste mass"].to_numpy()
remaining_water = results["water in tank"].to_numpy()
water_used = -remaining_water+300
max_places = argrelextrema(water_used, np.greater)
water = np.polynomial.Polynomial(least_squares(times[max_places], water_used[max_places]))
print(water)
print(water(21600))

fig = plt.figure(1)
plt.plot(times, water_used, 'r-', label="Simulation Results")
plt.plot(times, water(times), 'm--', label="Trendline of Maximum Water Usage")
plt.legend(loc="upper left")
plt.xlim([0, 10800])
plt.ylim([0, 225])
plt.xlabel("Time (s)")
plt.ylabel("Clean Water Used (L)")
plt.title("Clean Water Requirements for the Dialysis System")
fig.savefig("water_usage.png")

fig2 = plt.figure(2)
waste = np.polynomial.Polynomial(least_squares(times[np.nonzero(waste_output)], waste_output[np.nonzero(waste_output)]))
plt.plot(times, waste_output, 'r-', label="Simulation Results")
plt.plot(times, waste(times), 'm--', label="Trendline of Waste Production")
plt.legend(loc="upper left")
plt.xlim([0, 10800])
plt.ylim([-10, 90])
plt.yticks(ticks=np.linspace(-10, 90, 11))
plt.xlabel("Time (s)")
plt.ylabel("Waste Ouput (L)")
plt.title("Waste Production from the Dialysis Water Recycling")
fig2.savefig("waste_production.png")
print(waste)
print(waste(21600))