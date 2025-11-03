# SciPy ADM1 input array from Pettigrew (2017) jADM1 and Rosen et al (2006) BSM2 report
# initiate variables (initial values for the reactor state at t0)
import numpy as np



def get_initial_state(mixing_ratio):
    S_su = 0.00001
    S_aa = 0.00001
    S_fa = 0.00001
    S_va = 0.00001
    S_bu = 0.00001
    S_pro = 0.00001
    S_ac = 0.00001
    S_h2 = 10 ** -8
    S_ch4 = 10 ** -5
    S_IC = 0.04
    S_IN = 0.01
    S_I = 0.02

    # Feed1 initial
    X_xc1 = mixing_ratio * (2/(2+5+20+5)) * 0.00001
    X_ch1 = mixing_ratio * (5/(2+5+20+5)) * 0.00001
    X_pr1 = mixing_ratio * (20/(2+5+20+5)) * 0.00001
    X_li1 = mixing_ratio * (5/(2+5+20+5)) * 0.00001

    # Feed2 initial
    X_xc2 = (1-mixing_ratio) * (0/(205+47.76+6.925)) * 0.00001
    X_ch2 = (1-mixing_ratio) * (205/(205+47.76+6.925)) * 0.00001
    X_pr2 = (1-mixing_ratio) * (47.76/(205+47.76+6.925)) * 0.00001
    X_li2 = (1-mixing_ratio) * (6.925/(205+47.76+6.925)) * 0.00001

    variation_factor=5
    X_su = variation_factor*1.41
    X_aa = variation_factor*0.75
    X_fa = variation_factor*0.38
    X_c4 = variation_factor*0.35
    X_pro = variation_factor * 0.197
    X_ac = variation_factor * 0.95
    X_h2 = 0.466
    X_I = 25.6

    S_cation = 0.04000
    S_anion = 0.02

    pH = 7.4655377
    S_H_ion = 6.04e-8
    S_va_ion = 0.0082
    S_bu_ion = 0.0156
    S_pro_ion = 0.0172
    S_ac_ion = 0.11199
    S_hco3_ion = 0.15415
    S_nh3 = 0.0025
    S_nh4_ion = 0.126138
    S_co2 = 0.0093003
    S_gas_h2 = 4.91 * 10 ** -6
    S_gas_ch4 = 1.78
    S_gas_co2 = 0.025

    return {
        'S_su': S_su,
        'S_aa': S_aa,
        'S_fa': S_fa,
        'S_va': S_va,
        'S_bu': S_bu,
        'S_pro': S_pro,
        'S_ac': S_ac,
        'S_h2': S_h2,
        'S_ch4': S_ch4,
        'S_IC': S_IC,
        'S_IN': S_IN,
        'S_I': S_I,
        'X_xc1': X_xc1,
        'X_ch1': X_ch1,
        'X_pr1': X_pr1,
        'X_li1': X_li1,
        'X_xc2': X_xc2,
        'X_ch2': X_ch2,
        'X_pr2': X_pr2,
        'X_li2': X_li2,
        'X_su': X_su,
        'X_aa': X_aa,
        'X_fa': X_fa,
        'X_c4': X_c4,
        'X_pro': X_pro,
        'X_ac': X_ac,
        'X_h2': X_h2,
        'X_I': X_I,
        'S_cation': S_cation,
        'S_anion': S_anion,
        'pH': pH,
        'S_H_ion': S_H_ion,
        'S_va_ion': S_va_ion,
        'S_bu_ion': S_bu_ion,
        'S_pro_ion': S_pro_ion,
        'S_ac_ion': S_ac_ion,
        'S_hco3_ion': S_hco3_ion,
        'S_nh3': S_nh3,
        'S_nh4_ion': S_nh4_ion,
        'S_co2': S_co2,
        'S_gas_h2': S_gas_h2,
        'S_gas_ch4': S_gas_ch4,
        'S_gas_co2': S_gas_co2
    }








def get_SS_initial_state(mixing_ratio):
    S_su = 0.00001
    S_aa = 0.00001
    S_fa = 0.00001
    S_va = 0.00001
    S_bu = 0.00001
    S_pro = 0.00001
    S_ac = 0.00001
    S_h2 = 10 ** -8
    S_ch4 = 10 ** -5
    S_IC = 0.04
    S_IN = 0.01
    S_I = 0.02



    X_xc1 = mixing_ratio * 0
    X_ch1 = mixing_ratio * 53.21
    X_pr1 = mixing_ratio * 33.62*0.2
    X_li1 = mixing_ratio * 26.51
    X_I1 = mixing_ratio *224.85

    variation_factor4 = 1
    X_xc2_total = variation_factor4 * 259.992
    X_ch2_fraction = 0.79
    X_pr2_fraction = 0.184
    X_li2_fraction = 0.026

    X_xc2 = (1 - mixing_ratio) * 0
    X_ch2 = (1 - mixing_ratio) * X_ch2_fraction * X_xc2_total
    X_pr2 = (1 - mixing_ratio) * X_pr2_fraction * X_xc2_total
    X_li2 = (1 - mixing_ratio) * X_li2_fraction * X_xc2_total
    X_I2 = (1 - mixing_ratio) *19.023




    variation_factor=1
    X_su = variation_factor*1.41
    X_aa = variation_factor*0.75
    X_fa = variation_factor*0.38
    X_c4 = variation_factor*0.35
    X_pro = variation_factor * 0.197
    X_ac = variation_factor * 0.95*5
    X_h2 = 0.466
    X_I = X_I2 + X_I1

    S_cation = 0.04000
    S_anion = 0.02

    pH = 7.4655377
    S_H_ion = 6.04e-8
    S_va_ion = 0.0082
    S_bu_ion = 0.0156
    S_pro_ion = 0.0172
    S_ac_ion = 0.11199
    S_hco3_ion = 0.15415
    S_nh3 = 0.0025
    S_nh4_ion = 0.126138
    S_co2 = 0.0093003
    S_gas_h2 = 4.91 * 10 ** -6
    S_gas_ch4 = 1.78
    S_gas_co2 = 0.025

    return {
        'S_su': S_su,
        'S_aa': S_aa,
        'S_fa': S_fa,
        'S_va': S_va,
        'S_bu': S_bu,
        'S_pro': S_pro,
        'S_ac': S_ac,
        'S_h2': S_h2,
        'S_ch4': S_ch4,
        'S_IC': S_IC,
        'S_IN': S_IN,
        'S_I': S_I,
        'X_xc1': X_xc1,
        'X_ch1': X_ch1,
        'X_pr1': X_pr1,
        'X_li1': X_li1,
        'X_xc2': X_xc2,
        'X_ch2': X_ch2,
        'X_pr2': X_pr2,
        'X_li2': X_li2,
        'X_su': X_su,
        'X_aa': X_aa,
        'X_fa': X_fa,
        'X_c4': X_c4,
        'X_pro': X_pro,
        'X_ac': X_ac,
        'X_h2': X_h2,
        'X_I': X_I,
        'S_cation': S_cation,
        'S_anion': S_anion,
        'pH': pH,
        'S_H_ion': S_H_ion,
        'S_va_ion': S_va_ion,
        'S_bu_ion': S_bu_ion,
        'S_pro_ion': S_pro_ion,
        'S_ac_ion': S_ac_ion,
        'S_hco3_ion': S_hco3_ion,
        'S_nh3': S_nh3,
        'S_nh4_ion': S_nh4_ion,
        'S_co2': S_co2,
        'S_gas_h2': S_gas_h2,
        'S_gas_ch4': S_gas_ch4,
        'S_gas_co2': S_gas_co2
    }