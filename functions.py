import numpy as np

AMBIENT_HUMIDITY = 0.05

def num_diff(fun, x, **kwargs):
    eps = np.finfo(float).eps
    h = max(eps*abs(x), eps**(1/3))
    return (fun(x+h, **kwargs)-fun(x-h, **kwargs))/(2*h)

#smoothed square wave, switches every T units (default on)
def periodic_switch(x, T=1200):
    return (1 + np.sin(np.pi*x/T)/np.sqrt(np.sin(np.pi*x/T)**2+T**(-2)))/2

#logistic function, switches at T (default off)
def single_switch(x, T=0, k=1):
    return 1/(1+np.exp(-k*(x-T)))

#clamp function
def clamp(x, lowerlimit, upperlimit):
    if x < lowerlimit:
        x = lowerlimit
    if x > upperlimit:
        x = upperlimit
    return x

#modified from javascript implementation on wikipedia
def pascal_triangle(a, b):
    result = 1
    for i in range(b):
        result *= (a-i)/(i+1)
    return result

#step function for switching
def smooth_step(x, N=5):
    x = clamp(x, 0, 1)
    result = 0
    for n in range(N+1):
        result += pascal_triangle(-N-1, n)*pascal_triangle(2*N+1, N-n)*(x**(N+n+1))
    return result

#Buck equation (in kPa), temperature in celsius
def vapour_pressure(temperature):
    return 0.61121*np.exp((18.678-temperature/234.5)*(temperature/(257.14+temperature)))

#From https://vortex.plymouth.edu/~stmiller/stmiller_content/Publications/AtmosRH_Equations_Rev.pdf, pressure in kPa
def pressure_to_mass_ratio(pressure):
    eps=0.622
    atmospheric_pressure = 101.325
    mass_ratio = pressure*eps/(atmospheric_pressure-(1-eps)*pressure)
    if mass_ratio < 0:
        return 144
    return smooth_min(mass_ratio, 144)

#From https://docs.vaisala.com/r/M211280EN-D/en-US/GUID-905DDD94-2974-479D-8DCD-33811A9A081B, pressure in kPa, temperature in celcius
def pressure_to_absolute_humidity(pressure, temperature):
    c = 216.679e-3
    abs_hum = (pressure*10) * c / (temperature + 273.15)
    return abs_hum

#pressure in kpa, abs hum in kg/m^3, temp in celcius
def absolute_humidity_to_pressure(absolute_humidity, temperature):
    c = 216.679e-3
    pressure = absolute_humidity*(temperature+273.15)/(10*c)
    return pressure

#absolute humidity in kg/m^3, temperate in degrees celcius
def absolute_humidity_to_mass_ratio(absolute_humidity, temperature):
    return absolute_humidity/air_density(temperature, absolute_humidity_to_pressure(absolute_humidity, temperature))

#temperature in degrees celcius, pressure in kPa
def air_density(temperature, vapour_pressure=0, pressure=101.3):
    Rd = 287.058
    Rv = 461.495
    return 1e3*((pressure-vapour_pressure)/(Rd*(temperature+273.15))+vapour_pressure/(Rv*(temperature+273.15)))

#takes flow rate of input streams provides a quantity for the output stream
def combine_streams(flowrates, quantities):
    if np.sum(flowrates) == 0:
        return 0.0
    else:
        return np.dot(flowrates, quantities)/np.sum(flowrates)

#using a first order taylor series as an approximation, units are kg, m and s. positive indicates that it is flow from 1 to 2.
def ficks_law_mass_flow(abshum_1, abshum_2, A=0.25, L=5e-6, D=2.6e-1):
    MM_water = 18.016e-3
    return D*A*MM_water*(abshum_1-abshum_2)/L

#LSE, from wikipedia
def smooth_max(a, b):
    return np.log(np.exp(a)+np.exp(b))

def smooth_min(a, b):
    return -smooth_max(-a,-b)

def smooth_max2(a, b, n=1e-3):
    return (a+b+np.sqrt((a-b)**2+n))/2

def smooth_min2(a, b, n=1e-3):
    return -smooth_max2(-a, -b, n)

def vapor_mass_in_container(temperature, volume):
    abshum = pressure_to_absolute_humidity(vapour_pressure(temperature), temperature)
    return volume*abshum

def flow_rate_container(mass, scaling_factor=1):
    return smooth_max2(np.exp(mass/scaling_factor), 0, 1e-6)

def smooth_clamp(x, a, b, n=1e-3):
    return smooth_min2(smooth_max2(x, a, n), b, n)