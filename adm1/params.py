import numpy as np
from adm1.initial_state import get_initial_state




# ADM1  suggested biochemical parameter sets for different temperature phases (from Batstone et al. 2002 Table 6.2)

MESOPHILIC_HIGH_RATE = {
	"k_dis1": 0.4, #d^-1
    "k_dis2": 0.4, #d^-1
	"k_hyd_ch1": 0.1, #d^-1
	"k_hyd_pr1": 0.1, #d^-1
	"k_hyd_li1": 0.1, #d^-1
    "k_hyd_ch2": 5.22, #d^-1
	"k_hyd_pr2": 1.86, #d^-1
	"k_hyd_li2": 1.24, #d^-1
    #_______________________
	"k_dec_X_su":  0.02, #d^-1
    "k_dec_X_aa":  0.02, #d^-1
    "k_dec_X_fa":  0.02, #d^-1
    "k_dec_X_c4":  0.02, #d^-1
    "k_dec_X_pro":  0.02, #d^-1
    "k_dec_X_ac":  0.02, #d^-1
    "k_dec_X_h2":  0.02, #d^-1
    "K_S_IN":  10 ** -4, #M
    #_______________________
    "pH_UL_aa":  5.5,
    "pH_LL_aa":  4,
    #_______________________
    "k_m_su":  30, #d^-1
    "K_S_su":  0.5, #kgCOD.m^-3
    "Y_su":  0.1, #COD.COD^-1
    #_______________________
    "k_m_aa":  50, #d^-1
    "K_S_aa":  0.3, ##kgCOD.m^-3
    "Y_aa":  0.08, #COD.COD^-1
    #_______________________
    "k_m_fa":  6, #d^-1
    "K_S_fa":  0.4, #kgCOD.m^-3
    "Y_fa":  0.06, #COD.COD^-1
    "K_I_h2_fa":  5 * 10 ** -6, #kgCOD.m^-3
    #_______________________
    "k_m_c4":  20, #d^-1
    "K_S_c4":  0.3, #kgCOD.m^-3
    "Y_c4":  0.04, #COD.COD^-1
    "K_I_h2_c4":  3.5 ** -6, #kgCOD.m^-3
    #_______________________
    "k_m_pro":  13, #d^-1
    "K_S_pro":  0.3, #kgCOD.m^-3
    "Y_pro":  0.04, #COD.COD^-1
    "K_I_h2_pro":  3.5 * 10 ** -6, #kgCOD.m^-3
    #_______________________
    "k_m_ac":  8, #d^-1
    "K_S_ac":  0.15, #kgCOD.m^-3
    "Y_ac":  0.05, #COD.COD^-1
    "pH_UL_ac":  7,
    "pH_LL_ac":  6,
    "K_I_nh3":  0.0018, #M
    #_______________________
    "k_m_h2":  35, #d^-1
    "K_S_h2":  2.5 * 10 ** -5, #kgCOD.m^-3
    "Y_h2":  0.06, #COD.COD^-1
    "pH_UL_h2":  6,
    "pH_LL_h2":  5
}

MESOPHILIC_SOLIDS = {
	"k_dis1": 0.5, #d^-1
    "k_dis2": 0.5,  #d^-1
	"k_hyd_ch1": 0.1, #d^-1
	"k_hyd_pr1": 0.1, #d^-1
	"k_hyd_li1": 0.1, #d^-1
    "k_hyd_ch2": 5.22, #d^-1
	"k_hyd_pr2": 1.86, #d^-1
	"k_hyd_li2": 1.24, #d^-1
    #_______________________
	"k_dec_X_su":  0.02, #d^-1
    "k_dec_X_aa":  0.02, #d^-1
    "k_dec_X_fa":  0.02, #d^-1
    "k_dec_X_c4":  0.02, #d^-1
    "k_dec_X_pro":  0.02, #d^-1
    "k_dec_X_ac":  0.02, #d^-1
    "k_dec_X_h2":  0.02, #d^-1
    "K_S_IN":  10 ** -4, #M
    #_______________________
    "pH_UL_aa":  5.5,
    "pH_LL_aa":  4,
    #_______________________
    "k_m_su":  30, #d^-1
    "K_S_su":  0.5, #kgCOD.m^-3
    "Y_su":  0.1, #COD.COD^-1
    #_______________________
    "k_m_aa":  50, #d^-1
    "K_S_aa":  0.3, ##kgCOD.m^-3
    "Y_aa":  0.08, #COD.COD^-1
    #_______________________
    "k_m_fa":  6, #d^-1
    "K_S_fa":  0.4, #kgCOD.m^-3
    "Y_fa":  0.06, #COD.COD^-1
    "K_I_h2_fa":  5 * 10 ** -6, #kgCOD.m^-3
    #_______________________
    "k_m_c4":  20, #d^-1
    "K_S_c4":  0.3, #kgCOD.m^-3
    "Y_c4":  0.04,
    "K_I_h2_c4":  3.5 ** -6, #kgCOD.m^-3
    #_______________________
    "k_m_pro":  13, #d^-1
    "K_S_pro":  0.1, #kgCOD.m^-3
    "Y_pro":  0.04, #COD.COD^-1
    "K_I_h2_pro":  3.5 * 10 ** -6, #kgCOD.m^-3
    #_______________________x
    "k_m_ac":  8, #d^-1
    "K_S_ac":  0.15, #kgCOD.m^-3
    "Y_ac":  0.05, #COD.COD^-1
    "pH_UL_ac":  7,
    "pH_LL_ac":  6,
    "K_I_nh3":  0.0018, #M
    #_______________________
    "k_m_h2":  35, #d^-1
    "K_S_h2":  7 * 10 ** -6, #kgCOD.m^-3
    "Y_h2":  0.06, #COD.COD^-1
    "pH_UL_h2":  6,
    "pH_LL_h2":  5
}



THERMOPHILIC_SOLIDS = {
	"k_dis1": 1.0 , #d^-1
    "k_dis2": 1.0 , #d^-1
	"k_hyd_ch1": 0.1, #d^-1
	"k_hyd_pr1": 0.1, #d^-1
	"k_hyd_li1": 0.1, #d^-1
    "k_hyd_ch2": 5.22, #d^-1
	"k_hyd_pr2": 1.86, #d^-1
	"k_hyd_li2": 1.24, #d^-1
    #_______________________
	"k_dec_X_su":  0.04, #d^-1
    "k_dec_X_aa":  0.04, #d^-1
    "k_dec_X_fa":  0.04, #d^-1
    "k_dec_X_c4":  0.04, #d^-1
    "k_dec_X_pro":  0.04, #d^-1
    "k_dec_X_ac":  0.04, #d^-1
    "k_dec_X_h2":  0.04, #d^-1
    "K_S_IN":  10 ** -4, #M
    #_______________________
    "pH_UL_aa":  5.5,
    "pH_LL_aa":  4,
    #_______________________
    "k_m_su":  70, #d^-1
    "K_S_su":  1, #kgCOD.m^-3
    "Y_su":  0.1, #COD.COD^-1
    #_______________________
    "k_m_aa":  70, #d^-1
    "K_S_aa":  0.3, ##kgCOD.m^-3
    "Y_aa":  0.08, #COD.COD^-1 
    #_______________________
    "k_m_fa": 10, #d^-1
    "K_S_fa":  0.4, #kgCOD.m^-3
    "Y_fa":  0.06, #COD.COD^-1
    "K_I_h2_fa":  5 * 10 ** -6, #kgCOD.m^-3
    #_______________________
    "k_m_c4":  30, #d^-1
    "K_S_c4":  0.4, #kgCOD.m^-3
    "Y_c4":  0.06, #COD.COD^-1
    "K_I_h2_c4":  3 ** -5, #kgCOD.m^-3
    #_______________________
    "k_m_pro":  20, #d^-1
    "K_S_pro":  0.3, #kgCOD.m^-3
    "Y_pro":  0.05, #COD.COD^-1
    "K_I_h2_pro":  1 * 10 ** -5, #kgCOD.m^-3
    #_______________________
    "k_m_ac":  16, #d^-1
    "K_S_ac":  0.3, #kgCOD.m^-3
    "Y_ac":  0.05, #COD.COD^-1
    "pH_UL_ac":  7,
    "pH_LL_ac":  6,
    "K_I_nh3":  0.011, #M
    #_______________________
    "k_m_h2":  35, #d^-1
    "K_S_h2":  5 * 10 ** -5, #kgCOD.m^-3
    "Y_h2":  0.06, #COD.COD^-1 
    "pH_UL_h2":  6,
    "pH_LL_h2":  5
}

PARAMETER_SETS = {
	"mesophilic_high_rate": MESOPHILIC_HIGH_RATE,
	"mesophilic_solids": MESOPHILIC_SOLIDS,
	"thermophilic_solids": THERMOPHILIC_SOLIDS,
}
#########################

# Choose scenario: 'mesophilic_high_rate', 'mesophilic_solids', or 'thermophilic_solids'
#params = PARAMETER_SETS["mesophilic_solids"]


# Extract parameters from the selected parameter set


def set_global_params_from_dict(param_dict):
  #  globals_ = globals()
    for k, v in param_dict.items():
        globals()[k] = v




def get_adm1_params(T_ad, T_base, parameter_set, mixing_ratio):
        import numpy as np
        R = 0.083145  # bar.M^-1.K^-1

        # Unpack the selected parameter set
        params = parameter_set.copy()

        # Calculated physico-chemical parameters
        params.update({
            "K_w": 10 ** -14.0 * np.exp((55900 / (100 * R)) * (1 / T_base - 1 / T_ad)),
            "K_a_va": 10 ** -4.86,
            "K_a_bu": 10 ** -4.82,
            "K_a_pro": 10 ** -4.88,
            "K_a_ac": 10 ** -4.76,
            "K_a_co2": 10 ** -6.35 * np.exp((7646 / (100 * R)) * (1 / T_base - 1 / T_ad)),
            "K_a_IN": 10 ** -9.25 * np.exp((51965 / (100 * R)) * (1 / T_base - 1 / T_ad)),
            "k_A_B_va": 10 ** 10,
            "k_A_B_bu": 10 ** 10,
            "k_A_B_pro": 10 ** 10,
            "k_A_B_ac": 10 ** 10,
            "k_A_B_co2": 10 ** 10,
            "k_A_B_IN": 10 ** 10,
            "p_atm": 1.013,
            "p_gas_h2o": 0.0313 * np.exp(5290 * (1 / T_base - 1 / T_ad)),
            "k_p": 5 * 10 ** 4,
            "k_L_a": 200.0,
            "K_H_co2": 0.035 * np.exp((-19410 / (100 * R)) * (1 / T_base - 1 / T_ad)),
            "K_H_ch4": 0.0014 * np.exp((-14240 / (100 * R)) * (1 / T_base - 1 / T_ad)),
            "K_H_h2": 7.8 * 10 ** -4 * np.exp(-4180 / (100 * R) * (1 / T_base - 1 / T_ad)),
        })

        # pH inhibition parameters
        params.update({
            "K_pH_aa": 10 ** (-1 * (params["pH_LL_aa"] + params["pH_UL_aa"]) / 2.0),
            "nn_aa": 3.0 / (params["pH_UL_aa"] - params["pH_LL_aa"]),
            "K_pH_ac": 10 ** (-1 * (params["pH_LL_ac"] + params["pH_UL_ac"]) / 2.0),
            "n_ac": 3.0 / (params["pH_UL_ac"] - params["pH_LL_ac"]),
            "K_pH_h2": 10 ** (-1 * (params["pH_LL_h2"] + params["pH_UL_h2"]) / 2.0),
            "n_h2": 3.0 / (params["pH_UL_h2"] - params["pH_LL_h2"]),
        })

        # Composition fractions for particulate COD of each substrate
        params.update({
            "f_sI_xc1": 0.1,
            "f_xI_xc1": 0.2,
            "f_ch_xc1": 0.2,
            "f_pr_xc1": 0.2,
            "f_li_xc1": 0.3,
            "f_sI_xc2": 0.1,
            "f_xI_xc2": 0.2,
            "f_ch_xc2": 0.2,
            "f_pr_xc2": 0.2,
            "f_li_xc2": 0.3,
        })

        # Nitrogen and carbon parameters
        params.update({
            "N_xc": 0.0376 / 14,
            "N_I": 0.06 / 14,
            "N_aa": 0.007,
            "C_xc": mixing_ratio * 0.02786,
            "C_sI": mixing_ratio * 0.03,
            "C_ch": mixing_ratio * 0.0313,
            "C_pr": mixing_ratio * 0.03,
            "C_li": mixing_ratio * 0.022,
            "C_xI": mixing_ratio * 0.03,
            "C_xc_ref": 0.02786,
            "C_sI_ref": 0.03,
            "C_ch_ref": 0.0313,
            "C_pr_ref": 0.03,
            "C_li_ref": 0.022,
            "C_xI_ref": 0.03,
            "C_su": 0.0313,
            "C_aa": 0.03,
            "f_fa_li": 0.95,
            "C_fa": 0.0217,
            "f_h2_su": 0.19,
            "f_bu_su": 0.13,
            "f_pro_su": 0.27,
            "f_ac_su": 0.41,
            "N_bac": 0.08 / 14,
            "C_bu": 0.025,
            "C_pro": 0.0268,
            "C_ac": 0.0313,
            "C_bac": 0.0313,
            "f_h2_aa": 0.06,
            "f_va_aa": 0.23,
            "f_bu_aa": 0.26,
            "f_pro_aa": 0.05,
            "f_ac_aa": 0.40,
            "C_va": 0.024,
            "C_ch4": 0.0156,
            "R":R,
            "T_ad":T_ad,
            "T_op":T_ad,
            "T_base":T_base,
            "parameter_set":parameter_set,
            "mixing_ratio":mixing_ratio
        })

        return params



def get_VSS(influent, q):
    """
    Calculate VSS as the sum of lambda values multiplied by corresponding influent biomass fractions.
    Args:
        influent (dict): Dictionary containing influent values for each biomass fraction.
    Returns:
        float: Total VSS value.
    """
    lambda_dict = {
        'S_su_in': 1.07,
        'S_aa_in': 1.53,
        'S_fa_in': 2.87,
        'S_va_in': 2.04,
        'S_bu_in': 1.82,
        'S_pro_in': 1.51,
        'S_ac_in': 1.07,
        'X_ch1_in': 1.18,
        'X_ch2_in': 1.18,
        'X_pr1_in': 1.53,
        'X_pr2_in': 1.53,
        'X_li1_in': 2.87,
        'X_li2_in': 2.87,
        'X_I_in': 1.18,
        'X_I1_in': 1.18,
        'X_I2_in': 1.18,
        #######

        'X_su': 1.41,
        'X_aa': 1.41,
        'X_fa': 1.41,
        'X_c4': 1.41,
        'X_pro': 1.41,
        'X_ac': 1.41,
        'X_h2': 1.41,
        'S_su': 1.07,
        'S_aa': 1.53,
        'S_fa': 2.87,
        'S_va': 2.04,
        'S_va_ion': 1.98,
        'S_bu': 1.82,
        'S_bu_ion': 1.75,
        'S_pro': 1.51,
        'S_pro_ion': 1.42,
        'S_ac': 1.07,
        'S_ac_ion': 0.95,
        'X_ch1': 1.18,
        'X_ch2': 1.18,
        'X_pr1': 1.53,
        'X_pr2': 1.53,
        'X_li1': 2.87,
        'X_li2': 2.87,
        'X_I': 1.18,
        'S_h2': 7.94,
        'S_ch4': 3.99
    }
    VSS = 0.0
    for key, lam in lambda_dict.items():
        VSS += q*0.001 * influent.get(key, 0)/lam #tonne/day
    return VSS



# Mass fractions for C and N from Table S.1 (per unit mass)
# Only main organic components included; extend as needed
ADM1_MASS_FRACTIONS = {
    'X_ch1_in': (0.4445/1.18, 0.0000),  # Carbohydrates
    'X_pr1_in': (0.5497/1.53, 0.1603/1.53),  # Proteins
    'X_li1_in': (0.7587/2.87, 0.0000/2.87),  # Lipids
    'X_ch2_in': (0.4445/1.18, 0.0000),  # Carbohydrates
    'X_pr2_in': (0.5497/1.53, 0.1603/1.53),  # Proteins
    'X_li2_in': (0.7587/2.87, 0.0000),  # Lipids
    'X_I_in': (0.4445/1.18, 0.0000),  # Inerts
    'S_IN_in': (0.0000, 0.7995*18),  # Inorganic nitrogen

    # Add more as needed
}

def calculate_CN_ratio(influent, mass_fractions=ADM1_MASS_FRACTIONS):
    """
    Calculate the C/N ratio from influent composition and mass fractions.
    influent: dict, keys are component names, values are concentrations (e.g., kg/m3)
    mass_fractions: dict, keys are component names, values are (C_content, N_content) per unit mass
    Returns: dict with C/N ratio, total C, and total N
    """
    total_C = 0.0
    total_N = 0.0
    for comp, conc in influent.items():
        if comp in mass_fractions:
            C_content, N_content = mass_fractions[comp]
            total_C += conc * C_content
            total_N += conc * N_content
    if total_N == 0:
        cn_ratio = float('inf')
    else:
        cn_ratio = total_C / total_N
    return {"C/N ratio": cn_ratio, "total_C": total_C, "total_N": total_N}