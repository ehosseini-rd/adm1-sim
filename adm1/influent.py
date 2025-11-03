
import numpy as np
from typing import Optional
from adm1.params import get_VSS 


# ---  Reactor Setup ---
def reactor_setup(
    influent,               # Influent scenario function
    initial_state,          # Initial state scenario function
    q_ad_init,              # Initial influent flow rate [m^3/d]
    density,                # Influent density [tonne/m^3]
    VS_per_TS_PS,              # Volatile solids per total solids [kg VS/kg TS]
    VS_per_TS_SS,              # Volatile solids per total solids [kg VS/kg TS]
    TS_fraction,            # Fraction of water in influent
    mixing_ratio,           # Fraction for feed 1 # feed1 / total (can be overridden per scenario)
    mixing_ratio2,          # Fraction for VS calculation # feed1 / total (can be overridden per scenario)
    OLR,                    # Organic Loading Rate [kg VS/m3/d]
    recycle_ratio,          # Recycle ratio [-]
    Batch_process: bool,    # If True, simulate a batch process (no influent flow)
    VSS: Optional[float],  # Optional manual VSS override [tonne/day]; if None, it will be computed
    V_liq: Optional[float],
):
    """
    Calculate influent and reactor parameters for ADM1.

    Parameters
    - influent: influent scenario function or dict-like used by get_VSS
    - initial_state: initial state scenario function or dict-like used by get_VSS
    - q_ad_init: initial influent flow rate [m^3/d]
    - density: influent density [tonne/m^3]
    - VS_per_TS: volatile solids per total solids [kg VS/kg TS]
    - TS_fraction: fraction of water in influent [-]
    - mixing_ratio: feed split for stream 1 [-]
    - OLR: organic loading rate [kg VS/m3/d]
    - recycle_ratio: recycle ratio [-]
    - VSS: optional manual VSS [tonne/day]; if not provided, computed via get_VSS(influent, q_ad_init)

    Returns a dict with all relevant values.
    """

    VS_per_TS_Mix=mixing_ratio2 * VS_per_TS_PS + (1 - mixing_ratio2) * VS_per_TS_SS

    # Use provided VSS if given; otherwise compute from influent scenario
    VSS = VSS if VSS is not None else get_VSS(influent, q_ad_init)  # tonne/day

    TS=VSS / VS_per_TS_Mix  # Total solids [tonne TS/day]
    TS_fraction_initial=TS/(q_ad_init * density)  # Initial TS fraction [-]
    #water_content=(TS*(1-TS_fraction))/TS_fraction # [tonne/d]
    q_in = (TS_fraction_initial*q_ad_init)/TS_fraction  # New total influent flow [m^3/d]

    # Feed split
    q_in1 = mixing_ratio * q_in
    q_in2 = q_in - q_in1

    #Digestate (liquid effluent) flowrate (m^3/d)
    q_out=q_in/(1-recycle_ratio)

    # Recycle flowrate (m^3/d)
    q_r=recycle_ratio*q_out

    if Batch_process==True:
        V_liq=V_liq
        
        VSS = get_VSS(initial_state,V_liq)
        VS_in=VSS

        TS=VSS / VS_per_TS_Mix  # Total solids [tonne TS/day]
        TS_fraction_initial=TS/(V_liq * density)  # Initial TS fraction [-]
        #water_content=(TS*(1-TS_fraction))/TS_fraction # [tonne/d]
        V_liq = (TS_fraction_initial*V_liq)/TS_fraction  # New total influent flow [m^3/d]
        V_gas = 0.1 * V_liq # [m^3]
        V_ad = V_liq + V_gas # [m^3]

        q_in=10**-10

        q_in1 = mixing_ratio * q_in
        q_in2 = q_in - q_in1
        q_out=q_in/(1-recycle_ratio)

        q_ad=q_out

        HRT = 25  # [days]
    else:
        q_ad = q_out
        # Total VS in: use the (possibly overridden) VSS for the influent part
        VS_in = VSS #+ get_VSS(initial_state, q_ad, mixing_ratio)  # [tonne/d]

        # Recalculate VS fraction after evaporation
        #vs_frac_actual = vs_content / q_ad if q_ad > 0 else 0
        #VS_in = density * q_ad * vs_frac_actual * 0.001  # [tonne/d]

        # Hydraulic Retention Time (HRT) and reactor volumes
        HRT = VS_in*1000 / (q_ad*OLR*(1+recycle_ratio)) if OLR > 0 else 0  # [days]
        V_liq = HRT * q_ad  # [m^3]
        V_gas = 0.1 * V_liq # [m^3]
        V_ad = V_liq + V_gas # [m^3]

    return {
        "TS_fraction_initial": TS_fraction_initial,
        'q_in': q_in,
        'q_in1': q_in1,
        'q_in2': q_in2,
        'q_ad': q_ad,
        'q_out': q_out,
        'q_r': q_r,
        'VS_in': VS_in,
        'HRT': HRT,
        'V_liq': V_liq,
        'V_gas': V_gas,
        'V_ad': V_ad,
        'OLR': OLR,
        'density': density,
        'mixing_ratio': mixing_ratio,
        'TS': TS,
        'VSS': VSS
    }



# --- Recycle mixing helper ---
def mix_influent_with_recycle(
    influent: dict,
    effluent: dict,
    q_in: float,
    q_r: float,
    q_ad: float,
) -> dict:
    """
    Flow-weighted mixing of fresh influent with recycled effluent.

    new_influent = (q_in * influent + q_r * effluent) / q_ad

    Notes:
    - Influent keys are like 'S_su_in', 'X_ch1_in', ...
    - Effluent keys are the base state names without '_in' (e.g., 'S_su', 'X_ch1', ...).
    - For each influent key k, we map to base_key = k without the trailing '_in' and look it up in effluent.
    - If missing, we fall back to the influent value (i.e., no recycle contribution for that key).
    - If q_ad <= 0, returns the original influent to avoid division by zero.
    """
    if q_ad is None or q_ad <= 0:
        return dict(influent)

    mixed: dict = {}
    for k, v_in in influent.items():
        # Map 'S_su_in' -> 'S_su' and so on.
        base_key = k[:-3] if k.endswith('_in') else k
        v_eff = effluent.get(base_key, None) if effluent is not None else None
        if v_eff is None:
            v_eff = v_in  # fallback: no effluent contribution available

        mixed[k] = (q_in * v_in + q_r * v_eff) / q_ad

    return mixed





# --- Influent scenario function ---
def get_influent(mixing_ratio,
                 mixing_ratio2):
    
    total_COD=338.33


    """Base influent concentrations at reference flow q_ref (if provided).

    If q_ref supplied, concentrations are defined at that flow; later you can call
    rescale_influent(...) when actual q_in differs.
    """
    S_su_in = 0.001
    S_aa_in = 0.001
    S_fa_in = 0.001
    S_va_in = 0.001
    S_bu_in = 0.001
    S_pro_in = 0.001
    S_ac_in = 0.001
    S_h2_in = 10 ** -8
    S_ch4_in = 10 ** -5
    S_co2_in = 10 ** -5
    S_IC_in = 0.04
    S_IN_in = 0
    S_I_in = 0.02

    #Primary sludge

    X_xc1_total_PS = total_COD/(2*mixing_ratio2+1/mixing_ratio2-2)
    X_ch1_fraction_PS = 0.2041
    X_pr1_fraction_PS = 0.0045
    X_li1_fraction_PS = 0.0828
    X_I1_in_PS= 1-(X_ch1_fraction_PS + X_pr1_fraction_PS + X_li1_fraction_PS)


    X_xc1_in_PS = mixing_ratio * 0
    X_ch1_in_PS = mixing_ratio * X_ch1_fraction_PS * X_xc1_total_PS
    X_pr1_in_PS = mixing_ratio * X_pr1_fraction_PS * X_xc1_total_PS
    X_li1_in_PS = mixing_ratio * X_li1_fraction_PS * X_xc1_total_PS
    X_I1_in_PS  = mixing_ratio * X_I1_in_PS * X_xc1_total_PS


    #secondary sludge

    X_xc1_total_SS = X_xc1_total_PS*(1/mixing_ratio2-1)
    X_ch1_fraction_SS = 0.0647
    X_pr1_fraction_SS = 0.2822
    X_li1_fraction_SS = 0.0697
    X_I1_in_SS= 1-(X_ch1_fraction_SS + X_pr1_fraction_SS + X_li1_fraction_SS)


    X_xc1_in_SS = mixing_ratio * 0
    X_ch1_in_SS = mixing_ratio * X_ch1_fraction_SS * X_xc1_total_SS
    X_pr1_in_SS = mixing_ratio * X_pr1_fraction_SS * X_xc1_total_SS
    X_li1_in_SS = mixing_ratio * X_li1_fraction_SS * X_xc1_total_SS
    X_I1_in_SS  = mixing_ratio * X_I1_in_SS * X_xc1_total_SS



    # Mixture of primary and secondary sludge

    X_xc1_in = mixing_ratio * (mixing_ratio2*X_xc1_in_PS+(1-mixing_ratio2)*X_xc1_in_SS)
    X_ch1_in = mixing_ratio * (mixing_ratio2*X_ch1_in_PS+(1-mixing_ratio2)*X_ch1_in_SS)
    X_pr1_in = mixing_ratio * (mixing_ratio2*X_pr1_in_PS+(1-mixing_ratio2)*X_pr1_in_SS)
    X_li1_in = mixing_ratio * (mixing_ratio2*X_li1_in_PS+(1-mixing_ratio2)*X_li1_in_SS)
    X_I1_in= mixing_ratio *(mixing_ratio2*X_I1_in_PS+(1-mixing_ratio2)*X_I1_in_SS)




    variation_factor4 = 1
    X_xc2_total = variation_factor4 * 259.992
    X_ch2_fraction = 0.79
    X_pr2_fraction = 0.184
    X_li2_fraction = 0.026

    X_xc2_in = (1 - mixing_ratio) * 0
    X_ch2_in = (1 - mixing_ratio) * X_ch2_fraction * X_xc2_total
    X_pr2_in = (1 - mixing_ratio) * X_pr2_fraction * X_xc2_total
    X_li2_in = (1 - mixing_ratio) * X_li2_fraction * X_xc2_total
    X_I2_in = (1 - mixing_ratio) *19.023

    X_su_in = 0.0
    X_aa_in = 0.0
    X_fa_in = 0.0
    X_c4_in = 0.0
    X_pro_in = 0.0
    X_ac_in = 0.0
    X_h2_in = 0.0

    X_I_in = X_I2_in +X_I1_in
    S_cation_in = 0.04
    S_anion_in = 0.02
    S_nh3_in = 0.0

    return {
        'S_su_in': S_su_in,
        'S_aa_in': S_aa_in,
        'S_fa_in': S_fa_in,
        'S_va_in': S_va_in,
        'S_bu_in': S_bu_in,
        'S_pro_in': S_pro_in,
        'S_ac_in': S_ac_in,
        'S_h2_in': S_h2_in,
        'S_ch4_in': S_ch4_in,
        'S_co2_in': S_co2_in,
        'S_IC_in': S_IC_in,
        'S_IN_in': S_IN_in,
        'S_I_in': S_I_in,
        'X_xc1_in': X_xc1_in,
        'X_ch1_in': X_ch1_in,
        'X_pr1_in': X_pr1_in,
        'X_li1_in': X_li1_in,
        'X_xc2_in': X_xc2_in,
        'X_ch2_in': X_ch2_in,
        'X_pr2_in': X_pr2_in,
        'X_li2_in': X_li2_in,
        'X_su_in': X_su_in,
        'X_aa_in': X_aa_in,
        'X_fa_in': X_fa_in,
        'X_c4_in': X_c4_in,
        'X_pro_in': X_pro_in,
        'X_ac_in': X_ac_in,
        'X_h2_in': X_h2_in,
        'X_I_in': X_I_in,
        'S_cation_in': S_cation_in,
        'S_anion_in': S_anion_in,
        'S_nh3_in': S_nh3_in
        
    }


def rescale_influent(mixing_ratio,influent, q_in, q_ad_init):

    # Feed split
    q_in1 = mixing_ratio * q_in
    q_in2 = q_in - q_in1


    q_ad_init1 = mixing_ratio * q_ad_init
    q_ad_init2 = q_ad_init - q_ad_init1 

    S_su_in = influent['S_su_in']*q_ad_init/q_in
    S_aa_in = influent['S_aa_in']*q_ad_init/q_in
    S_fa_in = influent['S_fa_in']*q_ad_init/q_in
    S_va_in = influent['S_va_in']*q_ad_init/q_in
    S_bu_in = influent['S_bu_in']*q_ad_init/q_in
    S_pro_in = influent['S_pro_in']*q_ad_init/q_in
    S_ac_in = influent['S_ac_in']*q_ad_init/q_in
    S_h2_in = influent['S_h2_in']*q_ad_init/q_in
    S_ch4_in = influent['S_ch4_in']*q_ad_init/q_in
    S_co2_in = influent['S_co2_in']*q_ad_init/q_in
    S_IC_in = influent['S_IC_in']*q_ad_init/q_in
    S_IN_in = influent['S_IN_in']*q_ad_init/q_in
    S_I_in = influent['S_I_in']*q_ad_init/q_in



    X_xc1_in = influent['X_xc1_in']*q_ad_init1/q_in1
    X_ch1_in = influent['X_ch1_in']*q_ad_init1/q_in1
    X_pr1_in = influent['X_pr1_in']*q_ad_init1/q_in1
    X_li1_in = influent['X_li1_in']*q_ad_init1/q_in1


    X_xc2_in = influent['X_xc2_in']*q_ad_init2/q_in2
    X_ch2_in = influent['X_ch2_in']*q_ad_init2/q_in2
    X_pr2_in = influent['X_pr2_in']*q_ad_init2/q_in2
    X_li2_in = influent['X_li2_in']*q_ad_init2/q_in2

    X_su_in = influent['X_su_in']*q_ad_init/q_in
    X_aa_in = influent['X_aa_in']*q_ad_init/q_in
    X_fa_in = influent['X_fa_in']*q_ad_init/q_in
    X_c4_in = influent['X_c4_in']*q_ad_init/q_in
    X_pro_in = influent['X_pro_in']*q_ad_init/q_in
    X_ac_in = influent['X_ac_in']*q_ad_init/q_in
    X_h2_in = influent['X_h2_in']*q_ad_init/q_in
    X_I_in = influent['X_I_in']*q_ad_init/q_in
    S_cation_in = influent['S_cation_in']*q_ad_init/q_in
    S_anion_in = influent['S_anion_in']*q_ad_init/q_in
    S_nh3_in = influent['S_nh3_in']*q_ad_init/q_in

    return {
        'S_su_in': S_su_in,
        'S_aa_in': S_aa_in,
        'S_fa_in': S_fa_in,
        'S_va_in': S_va_in,
        'S_bu_in': S_bu_in,
        'S_pro_in': S_pro_in,
        'S_ac_in': S_ac_in,
        'S_h2_in': S_h2_in,
        'S_ch4_in': S_ch4_in,
        'S_co2_in': S_co2_in,
        'S_IC_in': S_IC_in,
        'S_IN_in': S_IN_in,
        'S_I_in': S_I_in,
        'X_xc1_in': X_xc1_in,
        'X_ch1_in': X_ch1_in,
        'X_pr1_in': X_pr1_in,
        'X_li1_in': X_li1_in,
        'X_xc2_in': X_xc2_in,
        'X_ch2_in': X_ch2_in,
        'X_pr2_in': X_pr2_in,
        'X_li2_in': X_li2_in,
        'X_su_in': X_su_in,
        'X_aa_in': X_aa_in,
        'X_fa_in': X_fa_in,
        'X_c4_in': X_c4_in,
        'X_pro_in': X_pro_in,
        'X_ac_in': X_ac_in,
        'X_h2_in': X_h2_in,
        'X_I_in': X_I_in,
        'S_cation_in': S_cation_in,
        'S_anion_in': S_anion_in,
        'S_nh3_in': S_nh3_in
    }



