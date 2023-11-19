from functions import *
import numpy as np

#mass flow rate for the dialysis output
def mdot_ab(t, T=21600, mdot=0.1):
    return mdot*(1-single_switch(t, T))

#input_temp = temp_u
def temp_ab(input_temp):
    return input_temp + 7

TDS_ab = 9000

#from the waste switch to the combiner
def mdot_c(pmdot_l, pmdot_p, tds, limit=108000):
    return (1-single_switch(tds, limit))*(pmdot_l+pmdot_p)

def tds_cs(pmdot_l, pmdot_p, tds_bed1, tds_bed2):
    return combine_streams(np.array([pmdot_l, pmdot_p]), np.array([tds_bed1, tds_bed2]))

def temp_cs(pmdot_l, pmdot_p, temp_bed1, temp_bed2):
    return combine_streams(np.array([pmdot_l, pmdot_p]), np.array([temp_bed1, temp_bed2]))

#from the combiner to the condenser hex
def mdot_d(water_mass):
    return flow_rate_container(water_mass)
#TDS_d = TDS_combiner
#temp_d = temp_combiner

#from bed 1 to the condenser (vapour)
def mdot_e(mvapor_bed, temp_condenser, t, volume_bed=0.7):
    abshum_bed = mvapor_bed/volume_bed
    return (1-periodic_switch(t))*ficks_law_mass_flow(abshum_bed, AMBIENT_HUMIDITY*pressure_to_absolute_humidity(vapour_pressure(temp_condenser), temp_condenser))
#temp_e = temp_bed1

#from bed 2 to the condenser (vapour)
def mdot_f(mvapor_bed, temp_condenser, t, volume_bed=0.7):
    abshum_bed = mvapor_bed/volume_bed
    return periodic_switch(t)*ficks_law_mass_flow(abshum_bed, AMBIENT_HUMIDITY*pressure_to_absolute_humidity(vapour_pressure(temp_condenser), temp_condenser))
#temp_f = temp_bed2

#from condenser to water tank
#mdot_g = mdot_e + mdot_f
#temp_g = temp_condenser

#from condenser to bed 1 heat exchanger
#tds_h = tds_combiner
def mdot_h(pmdot_d, t):
    return pmdot_d*periodic_switch(t)

#temp_hi = temp_combiner

#from condenser to bed 2 heat exchanger
def mdot_i(pmdot_d, t):
    return pmdot_d*(1-periodic_switch(t))

#from the evaporator to bed 1 (vapour)
def mdot_j(mvapor_evap, mvapor_bed, t, volume_bed=0.7, volume_evap=0.05):
    abshum_evap = mvapor_evap/volume_evap
    abshum_bed = mvapor_bed/volume_bed
    return periodic_switch(t)*ficks_law_mass_flow(abshum_evap, abshum_bed)

#for the stream splitting, fraction going to the output
def evaporator_splitting(abshum_evap, temp_evap, clamping=True):
    gamma = 1-1/144*(pressure_to_mass_ratio(vapour_pressure(temp_evap))-absolute_humidity_to_mass_ratio(abshum_evap, temp_evap))
    if clamping:
        return smooth_clamp(gamma, 0, 1)
    return gamma
    
#from the evaporator to bed 1 heat exchanger
def mdot_k(mvapor_evap, temp_evap, pmdot_r, t, volume_evap=0.05):
    abshum_evap = mvapor_evap/volume_evap
    return (1-periodic_switch(t))*evaporator_splitting(abshum_evap, temp_evap)*pmdot_r
#temp_k = temp_solar

def tds_kn(tds_r, mvapor_evap, temp_evap, volume_evap=0.05):
    abshum_evap = mvapor_evap/volume_evap
    return tds_r/evaporator_splitting(abshum_evap, temp_evap)

#from bed 1 heat exchanger to the waste switch
def mdot_lq(mwater_bed, t, hex_capacity=1):
    return (1-periodic_switch(t))*flow_rate_container(mwater_bed, hex_capacity)
#temp_lm = temp_bed1
#tds_lm = tds_bed1

#from bed 1 heat exchamger to solar thermal
def mdot_mp(mwater_bed, t, hex_capacity=1):
    return periodic_switch(t)*flow_rate_container(mwater_bed, hex_capacity)

#from the evaporator to bed 2 heat exchanger
def mdot_n(mvapor_evap, temp_evap, pmdot_r, t, volume_evap=0.05):
    abshum_evap = mvapor_evap/volume_evap
    return periodic_switch(t)*evaporator_splitting(abshum_evap, temp_evap)*pmdot_r
#temp_n = temp_solar

#from the evaporator to bed 2 (vapour)
def mdot_o(mvapor_evap, mvapor_bed, t, volume_bed=0.7, volume_evap=0.05):
    abshum_evap = mvapor_evap/volume_evap
    abshum_bed = mvapor_bed/volume_bed
    return (1-periodic_switch(t))*ficks_law_mass_flow(abshum_evap, abshum_bed)

#from bed 2 heat exchanger to waste switch - p
#from bed 2 heat exchanger to solar thermal - q
#temp_pq = temp_bed2
#tds_pq = tds_bed2

#from solar thermal to evaporator
#mdot_r = mdot_m + mdot_q
def tds_r(pmdot_m, pmdot_q, tds_bed1, tds_bed2):
    return combine_streams(np.array([pmdot_m, pmdot_q]), np.array([tds_bed1, tds_bed2]))
#temp_r = temp_solar

#from waste switch to waste
def mdot_s(pmdot_l, pmdot_p, tds, limit=108000):
    return single_switch(tds, limit)*(pmdot_l+pmdot_p)

#from water tank to evaporator heat sink
def mdot_tu(t, mdot=0.1, treatment_length=21600):
    return mdot*(1-single_switch(t,treatment_length))
#temp_t = temp_tank

#from evaporator to dialysis
def temp_u(temp_tank, temp_evap, mvapor_evap, pmdot_t, pmdot_r, volume_evap=0.05):
    c=4.186
    heat_vap = 2257
    return temp_tank-(1-evaporator_splitting(mvapor_evap/volume_evap,temp_evap))*pmdot_r*heat_vap/(pmdot_t*c)

#flow from the bed to the silica gel
def mdot_silica(mvapor_bed, mwater_silica, temp_bed, volume_bed=0.7, mass_silica=900):
    #literature values from https://www.sciencedirect.com/science/article/pii/S0360544217305753
    Rp = 1.75e-3
    gas_constant = 8.314e3
    f0 = 60
    Dso = 5.9e-4
    Ea = 28.665
    n = 1.68
    energy = 167.74
    max_concentration = 0.36

    #calculated values
    abshum_bed = mvapor_bed/volume_bed
    abshum_sat = pressure_to_absolute_humidity(vapour_pressure(temp_bed), temp_bed)
    temp_kelvin = temp_bed + 273.15
    Ds = Dso*np.exp(-Ea/(gas_constant*temp_kelvin))
    k = f0*Ds/Rp**2
    equilibrium_concentration = max_concentration*np.exp(-(smooth_max2(gas_constant*temp_kelvin*np.log(abshum_sat/abshum_bed)/energy, 0))**n)
    return k*(equilibrium_concentration-(mwater_silica/mass_silica))*mass_silica

