import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd
from functions import *
import system

def save_results(filename, results_t, results_y):
    d = np.vstack((results_t,results_y)).T
    column_names = ["t", "mass combiner", "temperature combiner", "tds combiner", "vapour mass bed one", "silica gel water mass bed one", "hex water bed one", "temperature bed one", "tds hex bed one", "vapour mass bed two", "silica gel water mass bed two", "hex water bed two", 'temperature bed two', 'tds hex bed two', 'temperature solar thermal', 'waste mass', 'water in tank', 'tank temperature', 'vapour mass evaporator', 'temperature condenser']
    df = pd.DataFrame(data=d, columns=column_names)
    df.to_csv(filename, index=False)

#np.seterr(all="ignore")

results = solve_ivp(system.run_system, [0,10800], system.INITIAL_SYSTEM, method='BDF', rtol=1e-10, atol=1e-8)
print(results)

final_result = results.y[:,-1]
end_time = results.t[-1]
print(end_time)
print(final_results)

if results.success:
    save_results("./results.csv", results.t, results.y)
    
