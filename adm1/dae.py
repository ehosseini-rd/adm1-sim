import numpy as np
from adm1.params import *


# Pure DAE solver. It DOES NOT read or mutate module-level state.
# Inputs: current scalar state values needed for the algebraic subsystem.
# Returns: (S_H_ion, pH, S_h2, S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_nh3)
# Function for DAE equations adopted from the Rosen et al (2006) BSM2 report bmadm1_report
# Function for DAE equations adopted from the Rosen et al (2006) BSM2 report bmadm1_report

def DAESolve(state, state_input, params):
    for k, v in params.items():
        globals()[k] = v

    S_su_in = state_input[0]
    S_aa_in = state_input[1]
    S_fa_in = state_input[2]
    S_va_in = state_input[3]
    S_bu_in =  state_input[4]
    S_pro_in =  state_input[5]
    S_ac_in =  state_input[6]
    S_h2_in =   state_input[7]
    S_ch4_in = state_input[8]
    S_IC_in = state_input[9]
    S_IN_in =  state_input[10]
    S_I_in = state_input[11]

    X_xc1_in =  state_input[12]
    X_ch1_in = state_input[13]
    X_pr1_in = state_input[14]
    X_li1_in =  state_input[15]
    
    X_xc2_in =  state_input[16]
    X_ch2_in = state_input[17]
    X_pr2_in = state_input[18]
    X_li2_in =  state_input[19]  
        
    X_su_in =  state_input[20]
    X_aa_in =  state_input[21]
    X_fa_in =  state_input[22]
    X_c4_in =  state_input[23]
    X_pro_in =  state_input[24]
    X_ac_in =  state_input[25]
    X_h2_in =  state_input[26]
    X_I_in = state_input[27]
    S_cation_in = state_input[28]
    S_anion_in = state_input[29]
    (S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I,
     X_xc1, X_ch1, X_pr1, X_li1, X_xc2, X_ch2, X_pr2, X_li2,
     X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cation, S_anion,
     S_H_ion, S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_co2, S_nh3, S_nh4_ion,
     S_gas_h2, S_gas_ch4, S_gas_co2) = state
    tol = 1e-12; maxIter = 1000
    shdelta = 1.0; i = 0
    while abs(shdelta) > tol and i < maxIter:
        S_va_ion = K_a_va * S_va / (K_a_va + S_H_ion)
        S_bu_ion = K_a_bu * S_bu / (K_a_bu + S_H_ion)
        S_pro_ion = K_a_pro * S_pro / (K_a_pro + S_H_ion)
        S_ac_ion = K_a_ac * S_ac / (K_a_ac + S_H_ion)
        S_hco3_ion = K_a_co2 * S_IC / (K_a_co2 + S_H_ion)
        S_nh3 = K_a_IN * S_IN / (K_a_IN + S_H_ion)
        shdelta = (S_cation + (S_IN - S_nh3) + S_H_ion - S_hco3_ion - S_ac_ion/64 -
                   S_pro_ion/112 - S_bu_ion/160 - S_va_ion/208 - K_w / S_H_ion - S_anion)
        shgradeq = (1 + K_a_IN * S_IN / (K_a_IN + S_H_ion)**2 + K_a_co2 * S_IC / (K_a_co2 + S_H_ion)**2 +
                    (K_a_ac * S_ac)/(64*(K_a_ac + S_H_ion)**2) + (K_a_pro * S_pro)/(112*(K_a_pro + S_H_ion)**2) +
                    (K_a_bu * S_bu)/(160*(K_a_bu + S_H_ion)**2) + (K_a_va * S_va)/(208*(K_a_va + S_H_ion)**2) +
                    K_w / S_H_ion**2)
        


        S_H_ion -= shdelta / shgradeq
            
        if S_H_ion <= 0: S_H_ion = tol
        i += 1
    pH_value = -np.log10(S_H_ion)
    S_h2delta = 1.0; j = 0
    while abs(S_h2delta) > tol and j < maxIter:
        Rho_5 = k_m_su * (S_su / (K_S_su + S_su + 1e-12)) * X_su
        Rho_6 = k_m_aa * (S_aa / (K_S_aa + S_aa + 1e-12)) * X_aa
        Rho_7 = k_m_fa * (S_fa / (K_S_fa + S_fa + 1e-12)) * X_fa
        Rho_8 = k_m_c4 * (S_va / (K_S_c4 + S_va + 1e-12)) * X_c4 * (S_va / (S_bu + S_va + 1e-6))
        Rho_9 = k_m_c4 * (S_bu / (K_S_c4 + S_bu + 1e-12)) * X_c4 * (S_bu / (S_bu + S_va + 1e-6))
        Rho_10 = k_m_pro * (S_pro / (K_S_pro + S_pro + 1e-12)) * X_pro
        Rho_12 = k_m_h2 * (S_h2 / (K_S_h2 + S_h2 + 1e-12)) * X_h2
        p_gas_h2 = S_gas_h2 * R * T_ad / 16
        Rho_T_8 = k_L_a * (S_h2 - 16 * K_H_h2 * p_gas_h2)
        S_h2delta = (q_ad / V_liq * (S_h2_in - S_h2) + (1 - Y_su) * f_h2_su * Rho_5 + (1 - Y_aa) * f_h2_aa * Rho_6 +
                     (1 - Y_fa) * 0.3 * Rho_7 + (1 - Y_c4) * 0.15 * Rho_8 + (1 - Y_c4) * 0.2 * Rho_9 +
                     (1 - Y_pro) * 0.43 * Rho_10 - Rho_12 - Rho_T_8)
        h = max(1e-9, 1e-6 * S_h2)
        S_h2_fd = S_h2 + h
        Rho_12_fd = k_m_h2 * (S_h2_fd / (K_S_h2 + S_h2_fd + 1e-12)) * X_h2
        S_h2delta_fd = (q_ad / V_liq * (S_h2_in - S_h2_fd) + (1 - Y_su) * f_h2_su * Rho_5 + (1 - Y_aa) * f_h2_aa * Rho_6 +
                        (1 - Y_fa) * 0.3 * Rho_7 + (1 - Y_c4) * 0.15 * Rho_8 + (1 - Y_c4) * 0.2 * Rho_9 +
                        (1 - Y_pro) * 0.43 * Rho_10 - Rho_12_fd - Rho_T_8)
        dS = (S_h2delta_fd - S_h2delta) / h
        if dS == 0: break
        S_h2 -= S_h2delta / dS
        if S_h2 <= 0: S_h2 = tol
        j += 1
    S_va_ion = K_a_va * S_va / (K_a_va + S_H_ion)
    S_bu_ion = K_a_bu * S_bu / (K_a_bu + S_H_ion)
    S_pro_ion = K_a_pro * S_pro / (K_a_pro + S_H_ion)
    S_ac_ion = K_a_ac * S_ac / (K_a_ac + S_H_ion)
    S_hco3_ion = K_a_co2 * S_IC / (K_a_co2 + S_H_ion)
    S_nh3 = K_a_IN * S_IN / (K_a_IN + S_H_ion)
    S_nh4_ion = S_IN - S_nh3
    S_co2 = S_IC - S_hco3_ion
    new_state = list(state)
    new_state[7] = S_h2
    new_state[30] = S_H_ion
    new_state[31] = S_va_ion
    new_state[32] = S_bu_ion
    new_state[33] = S_pro_ion
    new_state[34] = S_ac_ion
    new_state[35] = S_hco3_ion
    new_state[36] = S_co2
    new_state[37] = S_nh3
    new_state[38] = S_nh4_ion
    return new_state, pH_value