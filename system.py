import streams
import numpy as np
from functions import *

#Array structure
#x0: mass combiner
#x1: temperature combiner
#x2: TDS combiner
#x3: vapour mass bed one
#x4: silica gel water mass bed one
#x5: hex water bed one
#x6: temperature bed one
#x7: tds hex bed one
#x8: vapour mass bed two
#x9: silica gel water mass bed two
#x10: hex water bed two
#x11: temperature bed two
#x12: tds hex bed two
#x13: temperature solar thermal (temperature of evaporator)
#x14: total waste
#x15: mass in tank
#x16: tank temperature
#x17: vapour mass in evaporator
#x18: temperature of condenser
INITIAL_SYSTEM = np.array([0.1, 30, 0, AMBIENT_HUMIDITY*vapor_mass_in_container(30, 0.7), 0, 0.1, 30, 0, AMBIENT_HUMIDITY*vapor_mass_in_container(30, 0.7), 0, 0.1, 30, 0, 30, 0, 300, 30, AMBIENT_HUMIDITY*vapor_mass_in_container(30, 0.05), 30])

def run_system(t, system_array):
    #constants
    thermal_mass_beds = 415.8 #kJ/K
    thermal_mass_solar = 4.4 #kJ/K
    solar_thermal_efficiency = 0.05
    solar_thermal_area = 1 #m^2
    solar_irradiance = 2.22e4 #kJ/m^2 (from BOM)
    daylight_length = 32400 #s (from BOM)
    specific_heat_adsorption = 2700 #kJ/kg
    c_water = 4.186 #kJ/kg.K
    thermal_mass_tank = 46.2 #kJ/K
    overall_heat_transfer = 500
    ambient_temperature = 30
    thermal_mass_condenser = 49.1 #kJ/K
    heat_vap = 2257 #kJ/kg

    #streams
    #dependent on system properties only
    mdot_ab = streams.mdot_ab(t)
    mdot_d = streams.mdot_d(system_array[0])
    mdot_e = streams.mdot_e(system_array[3], system_array[18], t)
    mdot_f = streams.mdot_f(system_array[8], system_array[18], t)
    mdot_j = streams.mdot_j(system_array[17], system_array[3], t)
    mdot_l = streams.mdot_lq(system_array[5], t)
    mdot_m = streams.mdot_mp(system_array[5], t)
    mdot_o = streams.mdot_o(system_array[17], system_array[8], t)
    mdot_p = streams.mdot_mp(system_array[10], t)
    mdot_q = streams.mdot_lq(system_array[10], t)
    mdot_t = mdot_u = streams.mdot_tu(t)
    #dependent on previously calucalated properties
    tds_c = streams.tds_cs(mdot_l, mdot_p, system_array[7], system_array[12])
    temp_c = streams.temp_cs(mdot_l, mdot_p, system_array[6], system_array[11])
    mdot_g = mdot_e + mdot_f
    mdot_h = streams.mdot_h(mdot_d, t)
    mdot_i = streams.mdot_i(mdot_d, t)
    mdot_r = mdot_m + mdot_q
    tds_r = streams.tds_r(mdot_m, mdot_q, system_array[7], system_array[12])
    #dependent on previous set of properties
    mdot_c = streams.mdot_c(mdot_l, mdot_p, tds_c)
    mdot_k = streams.mdot_k(system_array[17], system_array[13], mdot_r, t)
    tds_k = tds_n = streams.tds_kn(tds_r, system_array[17], system_array[13])
    mdot_n = streams.mdot_n(system_array[17], system_array[13], mdot_r, t)
    mdot_s = streams.mdot_s(mdot_l, mdot_p, tds_c)
    temp_u = streams.temp_u(system_array[16], system_array[13], system_array[17], mdot_t, mdot_r)
    #the final one
    temp_ab = streams.temp_ab(temp_u)

    #silica gel flow rates
    mdot_sg1 = streams.mdot_silica(system_array[3], system_array[4], system_array[6])
    mdot_sg2 = streams.mdot_silica(system_array[8], system_array[9], system_array[11])

    #balances
    xdot = np.zeros(19)
    xdot[0] = mdot_ab + mdot_c - mdot_d
    xdot[1] = (mdot_ab*temp_ab+mdot_c*temp_c-mdot_d*system_array[1])/system_array[0]
    xdot[2] = (mdot_ab*streams.TDS_ab+mdot_c*tds_c-mdot_d*system_array[2])/system_array[0]
    xdot[3] = mdot_j - mdot_e - mdot_sg1
    xdot[4] = mdot_sg1
    xdot[5] = mdot_h + mdot_k - mdot_l - mdot_m
    xdot[6] = (mdot_sg1*specific_heat_adsorption+c_water*(mdot_h*system_array[18]+mdot_k*system_array[13]-(mdot_m+mdot_l)*system_array[6]))/(c_water*system_array[5]+thermal_mass_beds)
    xdot[7] = (mdot_h*system_array[2]+mdot_k*tds_k-(mdot_l+mdot_m)*system_array[7])/system_array[5]
    xdot[8] = mdot_o - mdot_f - mdot_sg2
    xdot[9] = mdot_sg2
    xdot[10] = mdot_i + mdot_n - mdot_p - mdot_q
    xdot[11] = (mdot_sg2*specific_heat_adsorption+c_water*(mdot_i*system_array[18]+mdot_n*system_array[13]-(mdot_p+mdot_q)*system_array[11]))/(c_water*system_array[10]+thermal_mass_beds)
    xdot[12] = (mdot_i*system_array[2]+mdot_n*tds_n-(mdot_p+mdot_q)*system_array[12])/system_array[10]
    xdot[13] = (c_water*(mdot_m*system_array[6]+mdot_q*system_array[11]-mdot_r*system_array[13])+solar_thermal_efficiency*solar_thermal_area*solar_irradiance/daylight_length)/thermal_mass_solar
    xdot[14] = mdot_s
    xdot[15] = mdot_g - mdot_t
    xdot[16] = (mdot_g*c_water*system_array[18]-mdot_t*c_water*system_array[16]+(ambient_temperature-system_array[16])*overall_heat_transfer)/(thermal_mass_tank+system_array[15]*c_water)
    xdot[17] = mdot_r - (mdot_k + mdot_n + mdot_j + mdot_o)
    xdot[18] = (c_water*(mdot_d*system_array[1]+mdot_e*system_array[6]+mdot_f*system_array[11]-(mdot_g+mdot_h+mdot_i)*system_array[18])+mdot_g*heat_vap)/thermal_mass_condenser

    #return values
    return np.nan_to_num(xdot)
    
