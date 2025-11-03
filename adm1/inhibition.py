"""
Shared inhibition factor calculations for ADM1.

This module centralizes the computation of inhibition factors used in
both the ODE/DAE kernel (ode.py/dae.py) and the orchestration layer (coAD.py).

Contract:
- Inputs:
  - S_H_ion, S_IN, S_h2, S_nh3: current state values (floats)
  - params: dict-like with model parameters; supports overrides
  - disable_inhibition: bool, when True returns all inhibitions as 1.0
- Outputs: dict of individual inhibitions and composites (I_5..I_12)

Notes:
- Parameter keys expected in `params` (with defaults if missing):
  K_pH_aa, nn_aa, K_pH_ac, n_ac, K_pH_h2, n_h2, K_S_IN,
  K_I_h2_fa, K_I_h2_c4, K_I_h2_pro, K_I_nh3.
"""

from typing import Dict, Any
import numpy as np



def compute_inhibition_factors(
    *,
    S_H_ion: float,
    S_IN: float,
    S_h2: float,
    S_nh3: float,
    params: Dict[str, Any],
    disable_inhibition: bool = False,
) -> Dict[str, float]:
    """
    Compute individual and composite inhibition factors for ADM1.

    Parameters:
        S_H_ion (float): Hydrogen ion concentration.
        S_IN (float): Inorganic nitrogen concentration.
        S_h2 (float): Hydrogen concentration.
        S_nh3 (float): Ammonia concentration.
        params (Dict[str, Any]): Model parameters, can override defaults.
        disable_inhibition (bool): If True, all inhibitions are set to 1.0.

    Returns:
        Dict[str, float]: Dictionary of inhibition factors (individual and composite).
    """
    
    pH = - np.log10(S_H_ion)


    # Pull parameters with sensible defaults
    K_pH_aa = params.get("K_pH_aa", 1.0)
    nn_aa = params.get("nn_aa", 1.0)
    K_pH_ac = params.get("K_pH_ac", 1.0)
    n_ac = params.get("n_ac", 1.0)
    K_pH_h2 = params.get("K_pH_h2", 1.0)
    n_h2 = params.get("n_h2", 1.0)
    K_S_IN = params.get("K_S_IN", 1e-4)
    K_I_h2_fa = params.get("K_I_h2_fa", 1e-6)
    K_I_h2_c4 = params.get("K_I_h2_c4", 1e-6)
    K_I_h2_pro = params.get("K_I_h2_pro", 1e-6)
    K_I_nh3 = params.get("K_I_nh3", 0.0018)
    
    pH_UL_aa = params.get("pH_UL_aa", 5.5)
    pH_LL_aa = params.get("pH_LL_aa", 4)
    pH_UL_ac = params.get("pH_UL_ac", 7)
    pH_LL_ac = params.get("pH_LL_ac", 6)
    pH_UL_h2 = params.get("pH_UL_h2", 6)
    pH_LL_h2 = params.get("pH_LL_h2", 5)

    # Individual inhibition factors

    if pH < pH_UL_aa:
        I_pH_aa= np.exp(-3 * ((pH - pH_UL_aa) / (pH_UL_aa - pH_LL_aa)) ** 2)
    else:
        I_pH_aa= 1.0

    if pH < pH_UL_ac:
        I_pH_ac= np.exp(-3 * ((pH - pH_UL_ac) / (pH_UL_ac - pH_LL_ac)) ** 2)
    else:
        I_pH_ac= 1.0

    if pH < pH_UL_h2:
        I_pH_h2 =np.exp(-3 * ((pH - pH_UL_h2) / (pH_UL_h2 - pH_LL_h2)) ** 2)
    else:
        I_pH_h2 = 1.0
    


    I_IN_lim =  (1 / (1 + (K_S_IN / S_IN)))
    I_h2_fa =  (1 / (1 + (S_h2 / K_I_h2_fa)))
    I_h2_c4 =  (1 / (1 + (S_h2 / K_I_h2_c4)))
    I_h2_pro =  (1 / (1 + (S_h2 / K_I_h2_pro)))
    I_nh3 =  (1 / (1 + (S_nh3 / K_I_nh3)))


    '''
    I_pH_aa =  ((K_pH_aa ** nn_aa) / (S_H_ion ** nn_aa + K_pH_aa ** nn_aa))
    I_pH_ac =  ((K_pH_ac ** n_ac) / (S_H_ion ** n_ac + K_pH_ac ** n_ac))
    I_pH_h2 =  ((K_pH_h2 ** n_h2) / (S_H_ion ** n_h2 + K_pH_h2 ** n_h2))
    I_IN_lim =  (1 / (1 + (K_S_IN / S_IN)))
    I_h2_fa =  (1 / (1 + (S_h2 / K_I_h2_fa)))
    I_h2_c4 =  (1 / (1 + (S_h2 / K_I_h2_c4)))
    I_h2_pro =  (1 / (1 + (S_h2 / K_I_h2_pro)))
    I_nh3 =  (1 / (1 + (S_nh3 / K_I_nh3)))
    '''

    if disable_inhibition:
        I_pH_aa = 1
        I_pH_ac = 1
        I_pH_h2 = 1
        I_IN_lim = 1
        I_h2_fa = 1
        I_h2_c4 = 1
        I_h2_pro = 1
        I_nh3 = 1

    # Composite inhibitions following model structure
    I_5 = I_pH_aa * I_IN_lim
    I_6 = I_5
    I_7 = I_pH_aa * I_IN_lim * I_h2_fa
    I_8 = I_pH_aa * I_IN_lim * I_h2_c4
    I_9 = I_8
    I_10 = I_pH_aa * I_IN_lim * I_h2_pro
    I_11 = I_pH_ac * I_IN_lim * I_nh3
    I_12 = I_pH_h2 * I_IN_lim

    return {
        "I_pH_aa": I_pH_aa,
        "I_pH_ac": I_pH_ac,
        "I_pH_h2": I_pH_h2,
        "I_IN_lim": I_IN_lim,
        "I_h2_fa": I_h2_fa,
        "I_h2_c4": I_h2_c4,
        "I_h2_pro": I_h2_pro,
        "I_nh3": I_nh3,
        "I_5": I_5,
        "I_6": I_6,
        "I_7": I_7,
        "I_8": I_8,
        "I_9": I_9,
        "I_10": I_10,
        "I_11": I_11,
        "I_12": I_12,
    }
