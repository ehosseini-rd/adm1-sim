#from adm1.params import *  # imports model parameters and initial states


# Function for calculating the derivatives related to ADM1 system of equations from the Rosen et al (2006) BSM2 report.
# state_zero: current dynamic state vector (length 42)
# state_input: influent / feed state vector (length 30) for this timestep


def ADM1_ODE(t, state_zero, state_input,params):
  for k, v in params.items():
    globals()[k] = v
  global S_nh4_ion, S_co2, p_gas, q_gas, q_ch4
  S_su = state_zero[0]
  S_aa = state_zero[1]
  S_fa = state_zero[2]
  S_va = state_zero[3]
  S_bu = state_zero[4]
  S_pro = state_zero[5]
  S_ac = state_zero[6]
  S_h2 = state_zero[7]
  S_ch4 = state_zero[8]
  S_IC = state_zero[9]
  S_IN = state_zero[10]
  S_I = state_zero[11]
  X_xc1 = state_zero[12]
  X_ch1 = state_zero[13]
  X_pr1 = state_zero[14]
  X_li1 = state_zero[15]
  X_xc2 = state_zero[16]
  X_ch2 = state_zero[17]
  X_pr2 = state_zero[18]
  X_li2 = state_zero[19]
  X_su = state_zero[20]
  X_aa = state_zero[21]
  X_fa = state_zero[22]
  X_c4 = state_zero[23]
  X_pro = state_zero[24]
  X_ac = state_zero[25]
  X_h2 = state_zero[26]
  X_I = state_zero[27]
  S_cation =  state_zero[28]
  S_anion = state_zero[29]
  S_H_ion =  state_zero[30]
  S_va_ion = state_zero[31]
  S_bu_ion = state_zero[32]
  S_pro_ion = state_zero[33]
  S_ac_ion = state_zero[34]
  S_hco3_ion =  state_zero[35]
  S_co2 = state_zero[36]
  S_nh3 = state_zero[37]
  S_nh4_ion =  state_zero[38]
  S_gas_h2 = state_zero[39]
  S_gas_ch4 = state_zero[40]
  S_gas_co2 = state_zero[41]


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

  S_nh4_ion =  (S_IN - S_nh3)

  S_co2 =  (S_IC - S_hco3_ion)

  # Base inhibition factors via shared utility
  from .inhibition import compute_inhibition_factors

  inhib = compute_inhibition_factors(
    S_H_ion=S_H_ion,
    S_IN=S_IN,
    S_h2=S_h2,
    S_nh3=S_nh3,
    params=params,
    disable_inhibition=params.get('disable_inhibition', False),
  )


  I_5 = inhib["I_5"]
  I_6 = inhib["I_6"]
  I_7 = inhib["I_7"]
  I_8 = inhib["I_8"]
  I_9 = inhib["I_9"]
  I_10 = inhib["I_10"]
  I_11 = inhib["I_11"]
  I_12 = inhib["I_12"]
 

  # biochemical process rates from Rosen et al (2006) BSM2 report
    
  ###disintegration

  Rho_1_1 = k_dis1 * X_xc1   # Disintegration of Substrate 1 particulate
  Rho_1_2 = k_dis2 * X_xc2   # Disintegration of Substrate 2 particulate

  ###hydrolysis

  Rho_2_1 = k_hyd_ch1 * X_ch1   # Hydrolysis of Substrate1 carbohydrates
  Rho_3_1 = k_hyd_pr1 * X_pr1   # Hydrolysis of Substrate1 proteins
  Rho_4_1 = k_hyd_li1 * X_li1   # Hydrolysis of Substrate1 lipids

  Rho_2_2 = k_hyd_ch2 * X_ch2   # Hydrolysis of Substrate2 carbohydrates
  Rho_3_2 = k_hyd_pr2 * X_pr2   # Hydrolysis of Substrate2 proteins
  Rho_4_2 = k_hyd_li2 * X_li2   # Hydrolysis of Substrate2 lipids

######################
  Rho_5 =  k_m_su * S_su / (K_S_su + S_su) * X_su * I_5  # Uptake of sugars
  Rho_6 =  (k_m_aa * (S_aa / (K_S_aa + S_aa)) * X_aa * I_6)  # Uptake of amino-acids
  Rho_7 =  (k_m_fa * (S_fa / (K_S_fa + S_fa)) * X_fa * I_7)  # Uptake of LCFA (long-chain fatty acids)
  Rho_8 =  (k_m_c4 * (S_va / (K_S_c4 + S_va )) * X_c4 * (S_va / (S_bu + S_va + 1e-6)) * I_8)  # Uptake of valerate
  Rho_9 =  (k_m_c4 * (S_bu / (K_S_c4 + S_bu )) * X_c4 * (S_bu / (S_bu + S_va + 1e-6)) * I_9)  # Uptake of butyrate
  Rho_10 =  (k_m_pro * (S_pro / (K_S_pro + S_pro)) * X_pro * I_10)  # Uptake of propionate
  Rho_11 =  (k_m_ac * (S_ac / (K_S_ac + S_ac)) * X_ac * I_11)  # Uptake of acetate
  Rho_12 =  (k_m_h2 * (S_h2 / (K_S_h2 + S_h2)) * X_h2 * I_12)  # Uptake of hydrogen
  Rho_13 =  (k_dec_X_su * X_su)  # Decay of X_su
  Rho_14 =  (k_dec_X_aa * X_aa)  # Decay of X_aa
  Rho_15 =  (k_dec_X_fa * X_fa)  # Decay of X_fa
  Rho_16 =  (k_dec_X_c4 * X_c4)  # Decay of X_c4
  Rho_17 =  (k_dec_X_pro * X_pro)  # Decay of X_pro
  Rho_18 =  (k_dec_X_ac * X_ac)  # Decay of X_ac
  Rho_19 =  (k_dec_X_h2 * X_h2)  # Decay of X_h2

  # acid-base rates for the BSM2 ODE implementation from Rosen et al (2006) BSM2 report
  Rho_A_4 =  (k_A_B_va * (S_va_ion * (K_a_va + S_H_ion) - K_a_va * S_va))
  Rho_A_5 =  (k_A_B_bu * (S_bu_ion * (K_a_bu + S_H_ion) - K_a_bu * S_bu))
  Rho_A_6 =  (k_A_B_pro * (S_pro_ion * (K_a_pro + S_H_ion) - K_a_pro * S_pro))
  Rho_A_7 =  (k_A_B_ac * (S_ac_ion * (K_a_ac + S_H_ion) - K_a_ac * S_ac))
  Rho_A_10 =  (k_A_B_co2 * (S_hco3_ion * (K_a_co2 + S_H_ion) - K_a_co2 * S_IC))
  Rho_A_11 =  (k_A_B_IN * (S_nh3 * (K_a_IN + S_H_ion) - K_a_IN * S_IN))

  # gas phase algebraic equations from Rosen et al (2006) BSM2 report
  p_gas_h2 =  (S_gas_h2 * R * T_op / 16)
  p_gas_ch4 =  (S_gas_ch4 * R * T_op / 64)
  p_gas_co2 =  (S_gas_co2 * R * T_op)


  p_gas=  (p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o)
  q_gas =  (k_p * (p_gas- p_atm))
  if q_gas < 0:    q_gas = 0

  q_ch4 = q_gas * (p_gas_ch4/p_gas) # methane flow

  # gas transfer rates from Rosen et al (2006) BSM2 report
  Rho_T_8 =  (k_L_a * (S_h2 - 16 * K_H_h2 * p_gas_h2))
  Rho_T_9 =  (k_L_a * (S_ch4 - 64 * K_H_ch4 * p_gas_ch4))
  Rho_T_10 =  (k_L_a * (S_co2 - K_H_co2 * p_gas_co2))

  ##differential equaitons from Rosen et al (2006) BSM2 report
  # differential equations 1 to 12 (soluble matter)
  diff_S_su = q_ad / V_liq * (S_su_in - S_su) + Rho_2_1 +Rho_2_2 + (1 - f_fa_li) * Rho_4_1 +(1 - f_fa_li) * Rho_4_2 - Rho_5  # eq1

  diff_S_aa = q_ad / V_liq * (S_aa_in - S_aa) + Rho_3_1+ Rho_3_2 - Rho_6  # eq2

  diff_S_fa = q_ad / V_liq * (S_fa_in - S_fa) + (f_fa_li * Rho_4_1)+ (f_fa_li * Rho_4_2) - Rho_7  # eq3

  diff_S_va = q_ad / V_liq * (S_va_in - S_va) + (1 - Y_aa) * f_va_aa * Rho_6 - Rho_8  # eq4

  diff_S_bu = q_ad / V_liq * (S_bu_in - S_bu) + (1 - Y_su) * f_bu_su * Rho_5 + (1 - Y_aa) * f_bu_aa * Rho_6 - Rho_9  # eq5

  diff_S_pro = q_ad / V_liq * (S_pro_in - S_pro) + (1 - Y_su) * f_pro_su * Rho_5 + (1 - Y_aa) * f_pro_aa * Rho_6 + (1 - Y_c4) * 0.54 * Rho_8 - Rho_10  # eq6

  diff_S_ac = q_ad / V_liq * (S_ac_in - S_ac) + (1 - Y_su) * f_ac_su * Rho_5 + (1 - Y_aa) * f_ac_aa * Rho_6 + (1 - Y_fa) * 0.7 * Rho_7 + (1 - Y_c4) * 0.31 * Rho_8 + (1 - Y_c4) * 0.8 * Rho_9 + (1 - Y_pro) * 0.57 * Rho_10 - Rho_11  # eq7

  #diff_S_h2 is defined with DAE paralel equaitons

  diff_S_ch4 = q_ad / V_liq * (S_ch4_in - S_ch4) + (1 - Y_ac) * Rho_11 + (1 - Y_h2) * Rho_12 - Rho_T_9  # eq9


  ## eq10 start##
  s_1_1 =(-1 * C_xc + f_sI_xc1 * C_sI 
                  + f_ch_xc1 * C_ch 
                  + f_pr_xc1 * C_pr 
                  + f_li_xc1 * C_li 
                  + f_xI_xc1 * C_xI)

  s_1_2 =(-1 * (C_xc_ref-C_xc) + f_sI_xc2 * (C_sI_ref-C_sI)
                  + f_ch_xc2 * (C_ch_ref-C_ch)
                  + f_pr_xc2 * (C_pr_ref-C_pr)
                  + f_li_xc2 * (C_li_ref-C_li)
                  + f_xI_xc2 * (C_xI_ref-C_xI))

  s_1= s_1_1+s_1_2

  s_2 =  (-1 * C_ch -1 * (C_ch_ref-C_ch) + C_su)
  s_3 =  (-1 * C_pr -1 * (C_pr_ref-C_pr) + C_aa)
  s_4 =  (-1 * C_li -1 * (C_li_ref-C_li) + (1 - f_fa_li) * C_su + f_fa_li * C_fa)
  s_5 =  (-1 * C_su + (1 - Y_su) * (f_bu_su * C_bu + f_pro_su * C_pro + f_ac_su * C_ac) + Y_su * C_bac)
  s_6 =  (-1 * C_aa + (1 - Y_aa) * (f_va_aa * C_va + f_bu_aa * C_bu + f_pro_aa * C_pro + f_ac_aa * C_ac) + Y_aa * C_bac)
  s_7 =  (-1 * C_fa + (1 - Y_fa) * 0.7 * C_ac + Y_fa * C_bac)
  s_8 =  (-1 * C_va + (1 - Y_c4) * 0.54 * C_pro + (1 - Y_c4) * 0.31 * C_ac + Y_c4 * C_bac)
  s_9 =  (-1 * C_bu + (1 - Y_c4) * 0.8 * C_ac + Y_c4 * C_bac)
  s_10 =  (-1 * C_pro + (1 - Y_pro) * 0.57 * C_ac + Y_pro * C_bac)
  s_11 =  (-1 * C_ac + (1 - Y_ac) * C_ch4 + Y_ac * C_bac)
  s_12 =  ((1 - Y_h2) * C_ch4 + Y_h2 * C_bac)
  s_13 =  (-1 * C_bac + C_xc + (C_xc_ref-C_xc)) 

  Sigma =  (s_1 * (Rho_1_2 + Rho_1_1)
            + s_2 * (Rho_2_1 + Rho_2_2)
            + s_3 * (Rho_3_1 + Rho_3_2)
            + s_4 * (Rho_4_1 + Rho_4_2)
            + s_5 * Rho_5 + s_6 * Rho_6 + s_7 * Rho_7 
            + s_8 * Rho_8 + s_9 * Rho_9 + s_10 * Rho_10 
            + s_11 * Rho_11 + s_12 * Rho_12 + s_13 * (Rho_13 + Rho_14 + Rho_15 + Rho_16 + Rho_17 + Rho_18 + Rho_19))


  diff_S_IC = q_ad / V_liq * (S_IC_in - S_IC) - Sigma - Rho_T_10
  ## eq10 end## 


 
  diff_S_IN = q_ad / V_liq * (S_IN_in - S_IN) + (N_xc - f_xI_xc1 * N_I - f_sI_xc1 * N_I-f_pr_xc1 * N_aa) * Rho_1_1+ (N_xc - f_xI_xc2 * N_I - f_sI_xc2 * N_I-f_pr_xc2 * N_aa) * Rho_1_2 - Y_su * N_bac * Rho_5 + (N_aa - Y_aa * N_bac) * Rho_6 - Y_fa * N_bac * Rho_7 - Y_c4 * N_bac * Rho_8 - Y_c4 * N_bac * Rho_9 - Y_pro * N_bac * Rho_10 - Y_ac * N_bac * Rho_11 - Y_h2 * N_bac * Rho_12 + (N_bac - N_xc) * (Rho_13 + Rho_14 + Rho_15 + Rho_16 + Rho_17 + Rho_18 + Rho_19) # eq11 


  diff_S_I = q_ad / V_liq * (S_I_in - S_I) + f_sI_xc1 * Rho_1_1 + f_sI_xc2 * Rho_1_2  # eq12


  # Differential equations 13 to 24 (particulate matter)
  diff_X_xc1 = q_in1 / (mixing_ratio * V_liq) * (X_xc1_in - X_xc1) - Rho_1_1 + mixing_ratio * (Rho_13 + Rho_14 + Rho_15 + Rho_16 + Rho_17 + Rho_18 + Rho_19)  # eq13 

  diff_X_xc2 = q_in2 / ((1-mixing_ratio) * V_liq) * (X_xc2_in - X_xc2) - Rho_1_2 + (1-mixing_ratio) * (Rho_13 + Rho_14 + Rho_15 + Rho_16 + Rho_17 + Rho_18 + Rho_19)  # eq13 

  diff_X_ch1 = q_in1 / (mixing_ratio * V_liq) * (X_ch1_in - X_ch1) + f_ch_xc1 * Rho_1_1 - Rho_2_1  # eq14

  diff_X_ch2 = q_in2 / ((1-mixing_ratio) * V_liq) * (X_ch2_in - X_ch2) + f_ch_xc2 * Rho_1_2 - Rho_2_2  # eq14

  diff_X_pr1 = q_in1 / (mixing_ratio * V_liq) * (X_pr1_in - X_pr1) + f_pr_xc1 * Rho_1_1 - Rho_3_1    # eq15 

  diff_X_pr2 = q_in2 / ((1-mixing_ratio) * V_liq) * (X_pr2_in - X_pr2) + f_pr_xc2 * Rho_1_2 - Rho_3_2   # eq15 

  diff_X_li1 = q_in1 / (mixing_ratio * V_liq) * (X_li1_in - X_li1) + f_li_xc1 * Rho_1_1 - Rho_4_1    # eq16 

  diff_X_li2 = q_in2 / ((1-mixing_ratio) * V_liq) * (X_li2_in - X_li2) + f_li_xc2 * Rho_1_2 - Rho_4_2   # eq16 
  
  diff_X_su = q_ad / V_liq * (X_su_in - X_su) + Y_su * Rho_5 - Rho_13  # eq17

  diff_X_aa = q_ad / V_liq * (X_aa_in - X_aa) + Y_aa * Rho_6 - Rho_14  # eq18

  diff_X_fa = q_ad / V_liq * (X_fa_in - X_fa) + Y_fa * Rho_7 - Rho_15  # eq19

  diff_X_c4 = q_ad / V_liq * (X_c4_in - X_c4) + Y_c4 * Rho_8 + Y_c4 * Rho_9 - Rho_16  # eq20

  diff_X_pro = q_ad / V_liq * (X_pro_in - X_pro) + Y_pro * Rho_10 - Rho_17  # eq21

  diff_X_ac = q_ad / V_liq * (X_ac_in - X_ac) + Y_ac * Rho_11 - Rho_18  # eq22

  diff_X_h2 = q_ad / V_liq * (X_h2_in - X_h2) + Y_h2 * Rho_12 - Rho_19  # eq23

  diff_X_I = q_ad / V_liq * (X_I_in - X_I) + f_xI_xc1 * Rho_1_1 + f_xI_xc2 * Rho_1_2  # eq24 

  # Differential equations 25 and 26 (cations and anions)
  diff_S_cation = q_ad / V_liq * (S_cation_in - S_cation)  # eq25

  diff_S_anion = q_ad / V_liq * (S_anion_in - S_anion)  # eq26



  diff_S_h2 = 0

    # Differential equations 27 to 32 (ion states, only for ODE implementation)
  diff_S_va_ion = 0  # eq27

  diff_S_bu_ion = 0  # eq28

  diff_S_pro_ion = 0  # eq29

  diff_S_ac_ion = 0  # eq30

  diff_S_hco3_ion = 0  # eq31

  diff_S_nh3 = 0  # eq32

  # Gas phase equations: Differential equations 33 to 35
  diff_S_gas_h2 = (q_gas / V_gas * -1 * S_gas_h2) + (Rho_T_8 * V_liq / V_gas)  # eq33

  diff_S_gas_ch4 = (q_gas / V_gas * -1 * S_gas_ch4) + (Rho_T_9 * V_liq / V_gas)  # eq34

  diff_S_gas_co2 = (q_gas / V_gas * -1 * S_gas_co2) + (Rho_T_10 * V_liq / V_gas)  # eq35

  diff_S_H_ion = diff_S_co2 = diff_S_nh4_ion = 0 #to keep the output same length as input for ADM1_ODE funcion


  return diff_S_su, diff_S_aa, diff_S_fa, diff_S_va, diff_S_bu, diff_S_pro, diff_S_ac, diff_S_h2, diff_S_ch4, diff_S_IC, diff_S_IN, diff_S_I, diff_X_xc1, diff_X_ch1, diff_X_pr1, diff_X_li1, diff_X_xc2, diff_X_ch2, diff_X_pr2, diff_X_li2 ,diff_X_su, diff_X_aa, diff_X_fa, diff_X_c4, diff_X_pro, diff_X_ac, diff_X_h2, diff_X_I, diff_S_cation, diff_S_anion, diff_S_H_ion, diff_S_va_ion,  diff_S_bu_ion, diff_S_pro_ion, diff_S_ac_ion, diff_S_hco3_ion, diff_S_co2,  diff_S_nh3, diff_S_nh4_ion, diff_S_gas_h2, diff_S_gas_ch4, diff_S_gas_co2