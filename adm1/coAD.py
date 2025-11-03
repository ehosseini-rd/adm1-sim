import numpy as np
import pandas as pd
#from adm1.ode import ADM1_ODE
# Import custom plotting functions
from adm1.influent import reactor_setup
from adm1.influent import mix_influent_with_recycle
from adm1.influent import get_influent
from adm1.initial_state import get_initial_state
from adm1.params import PARAMETER_SETS, set_global_params_from_dict
from adm1.params import get_adm1_params
#from adm1.params import *
from adm1.dae import DAESolve  # pure DAE solver
from adm1.solver import simulate
from adm1.params import get_VSS 
from adm1.influent import rescale_influent
from .inhibition import compute_inhibition_factors


def ADM1_coAD(
    q_ad_init,              # Initial influent flow rate [m^3/d]
    density,               # Influent density [tonne/m^3]
    VS_per_TS_PS,            # Volatile solids per  total solids [kg VS/kg TS]     
    VS_per_TS_SS,
    TS_fraction,        # Fraction of TS in influent
    mixing_ratio,       # Fraction for feed 1 / total (can be overridden per scenario)
    mixing_ratio2,      # Fraction for PS / total (can be overridden per scenario)
    OLR,                      # Organic Loading Rate [kg VS/m3/d]
    T_ad, #K
    T_base, #K
    T_op, #k ##T_ad #=35 C
    recycle_ratio,
    influent,              # Optional: pass custom influent dict
    initials,              # Optional: pass custom initial state dict
    VSS,                    # Optional: pass custom VSS value
    days,                   # Number of days to simulate
    timesteps,              #"Day(s)", "15 Minute(s)","Hour(s)"
    V_liq,                    # Optional: pass custom reactor liquid volume [m^3]
    param_overrides,        # Optional: dict of parameter overrides (e.g., k_L_a, k_p, K_H_*)
    disable_inhibition: bool,  # If True, force all inhibition factors to 1
    Batch_process: bool, # If True, simulate a batch process (no influent flow)
):


    if influent is None:
        influent = get_influent(mixing_ratio, mixing_ratio2)
    else:
        influent = influent
    

    ########################################
    # Initial state
    if initials is None:
        initial_state = get_initial_state(mixing_ratio)
    else:
        initial_state = initials

    S_su=initial_state['S_su']
    S_aa=initial_state['S_aa']
    S_fa=initial_state['S_fa']
    S_va=initial_state['S_va']
    S_bu=initial_state['S_bu']
    S_pro=initial_state['S_pro']
    S_ac=initial_state['S_ac']
    S_h2=initial_state['S_h2']
    S_ch4=initial_state['S_ch4']
    S_IC=initial_state['S_IC']
    S_IN=initial_state['S_IN']
    S_I=initial_state['S_I']
    X_xc1=initial_state['X_xc1']
    X_ch1=initial_state['X_ch1']
    X_pr1=initial_state['X_pr1']
    X_li1=initial_state['X_li1']
    X_xc2=initial_state['X_xc2']
    X_ch2=initial_state['X_ch2']
    X_pr2=initial_state['X_pr2']
    X_li2=initial_state['X_li2']
    X_su=initial_state['X_su']
    X_aa=initial_state['X_aa']
    X_fa=initial_state['X_fa']
    X_c4=initial_state['X_c4']
    X_pro=initial_state['X_pro']
    X_ac=initial_state['X_ac']
    X_h2=initial_state['X_h2']
    X_I=initial_state['X_I']
    S_cation=initial_state['S_cation']
    S_anion=initial_state['S_anion']
    pH=initial_state['pH']
    S_H_ion=initial_state['S_H_ion']
    S_va_ion=initial_state['S_va_ion']
    S_bu_ion=initial_state['S_bu_ion']
    S_pro_ion=initial_state['S_pro_ion']
    S_ac_ion=initial_state['S_ac_ion']
    S_hco3_ion=initial_state['S_hco3_ion']
    S_nh3=initial_state['S_nh3']
    S_nh4_ion=initial_state['S_nh4_ion']
    S_co2=initial_state['S_co2']
    S_gas_h2=initial_state['S_gas_h2']
    S_gas_ch4=initial_state['S_gas_ch4']
    S_gas_co2=initial_state['S_gas_co2']

    state_zero = [S_su,
                S_aa,
                S_fa,
                S_va,
                S_bu,
                S_pro,
                S_ac,
                S_h2,
                S_ch4,
                S_IC,
                S_IN,
                S_I,
                X_xc1,
                X_ch1,
                X_pr1,
                X_li1,
                X_xc2,
                X_ch2,
                X_pr2,
                X_li2,
                X_su,
                X_aa,
                X_fa,
                X_c4,
                X_pro,
                X_ac,
                X_h2,
                X_I,
                S_cation,
                S_anion,
                S_H_ion,
                S_va_ion,
                S_bu_ion,
                S_pro_ion,
                S_ac_ion,
                S_hco3_ion,
                S_co2,
                S_nh3,
                S_nh4_ion,
                S_gas_h2,
                S_gas_ch4,
                S_gas_co2]

    
    # --- Set up reactor for this scenario ---
    reactor = reactor_setup(
    influent,               # Influent scenario function
    initial_state,          # Initial state scenario function
    q_ad_init,                # Initial influent flow rate [m^3/d]
    density,                    # Influent density [tonne/m^3]
    VS_per_TS_PS,                # Volatile solids per  total solids [kg VS/kg TS]
    VS_per_TS_SS,    
    TS_fraction,         # Fraction of water in influent
    mixing_ratio,   # Fraction for feed 1 # feed1 / total (can be overridden per scenario)
    mixing_ratio2,
    OLR,                        # Organic Loading Rate [kg VS/m3/d]
    recycle_ratio,                # Recycle ratio [m^3/d]
    Batch_process,                # If True, simulate a batch process (no influent flow)
    VSS,                            # Optional: pass custom VSS value
    V_liq=V_liq,                            # Optional: pass custom reactor liquid volume [m^3]
    )
    # --- Set up reactor for this scenario ---

    q_in = reactor['q_in']
    q_in1 = reactor['q_in1']
    q_in2 = reactor['q_in2']
    q_ad = reactor['q_ad']
    q_out = reactor['q_out']
    q_r = reactor['q_r']
    VS_in = reactor['VS_in']
    HRT = reactor['HRT']
    V_liq = reactor['V_liq']
    V_gas = reactor['V_gas']
    V_ad = reactor['V_ad']
    TS = reactor['TS']
    VSS = reactor['VSS']
    TS_fraction_initial = reactor['TS_fraction_initial']



    print("VS_in:", VS_in)
    print("HRT:", HRT)
    print("V_liq:", V_liq)
    print("V_gas:", V_gas)
    print("V_ad:", V_ad)
    print("TS:", TS)
    print("VSS:", VSS)
    print("TS_fraction_initial:", TS_fraction_initial)

    # --- End  Reactor Setup ---

    # Rescale influent concentrations to actual flow (preserve mass load)
    new_influent = rescale_influent(mixing_ratio,influent, q_in, q_ad_init)


    S_su_in = new_influent['S_su_in']
    S_aa_in = new_influent['S_aa_in']
    S_fa_in = new_influent['S_fa_in']
    S_va_in = new_influent['S_va_in']
    S_bu_in = new_influent['S_bu_in']
    S_pro_in = new_influent['S_pro_in']
    S_ac_in = new_influent['S_ac_in']
    S_h2_in = new_influent['S_h2_in']
    S_ch4_in = new_influent['S_ch4_in']
    S_IC_in = new_influent['S_IC_in']
    S_IN_in = new_influent['S_IN_in']
    S_I_in = new_influent['S_I_in']
    X_xc1_in = new_influent['X_xc1_in']
    X_ch1_in = new_influent['X_ch1_in']
    X_pr1_in = new_influent['X_pr1_in']
    X_li1_in = new_influent['X_li1_in']
    X_xc2_in = new_influent['X_xc2_in']
    X_ch2_in = new_influent['X_ch2_in']
    X_pr2_in = new_influent['X_pr2_in']
    X_li2_in = new_influent['X_li2_in']
    X_su_in = new_influent['X_su_in']
    X_aa_in = new_influent['X_aa_in']
    X_fa_in = new_influent['X_fa_in']
    X_c4_in = new_influent['X_c4_in']
    X_pro_in = new_influent['X_pro_in']
    X_ac_in = new_influent['X_ac_in']
    X_h2_in = new_influent['X_h2_in']
    X_I_in = new_influent['X_I_in']
    S_cation_in = new_influent['S_cation_in']
    S_anion_in = new_influent['S_anion_in']

    state_input = [S_su_in,
                S_aa_in,
                S_fa_in,
                S_va_in,
                S_bu_in,
                S_pro_in,
                S_ac_in,
                S_h2_in,
                S_ch4_in,
                S_IC_in,
                S_IN_in,
                S_I_in,
                X_xc1_in,
                X_ch1_in,
                X_pr1_in,
                X_li1_in,
                X_xc2_in,
                X_ch2_in,
                X_pr2_in,
                X_li2_in,
                X_su_in,
                X_aa_in,
                X_fa_in,
                X_c4_in,
                X_pro_in,
                X_ac_in,
                X_h2_in,
                X_I_in,
                S_cation_in,
                S_anion_in]

    S_co2 =  (S_IC - S_hco3_ion)

    #pH equation

    ######################

    # Physical parameter values

    t0=0
    n=0



    # Choose parameter set based on operating temperature (K) -> (°C)
    # Rule: 20–40°C => MESOPHILIC_SOLIDS, 45–70°C => THERMOPHILIC_SOLIDS, gap (40–45°C) defaults to mesophilic
    T_C_sel = float(T_ad) - 273.15
    if 20.0 <= T_C_sel <= 40.0:
        selected_set = "mesophilic_solids"
    elif 45.0 <= T_C_sel <= 70.0:
        selected_set = "thermophilic_solids"
    else:
        # Default to mesophilic in the transition gap or outside specified bounds
        selected_set = "mesophilic_solids"

    params2 = PARAMETER_SETS[selected_set]
    set_global_params_from_dict(params2)
    params = get_adm1_params(T_ad, T_base, params2, mixing_ratio)

    # Apply any caller-provided overrides after base/temperature-derived params are built
    if isinstance(param_overrides, dict) and param_overrides:
        params.update(param_overrides)

    params.update({
        'q_in': reactor['q_in'],
        'q_in1': reactor['q_in1'],
        'q_in2': reactor['q_in2'],
        'q_ad': reactor['q_ad'],
        'q_out': reactor['q_out'],
        'q_r': reactor['q_r'],
        'VS_in': reactor['VS_in'],
        'HRT': reactor['HRT'],
        'V_liq': reactor['V_liq'],
        'V_gas': reactor['V_gas'],
        'V_ad': reactor['V_ad'],
        'OLR': reactor['OLR'],
        'density': reactor['density'],
        'mixing_ratio': reactor['mixing_ratio'],
        # Toggle for inhibition in ODE path
        'disable_inhibition': disable_inhibition,
    })

    for k, v in params.items():
        globals()[k] = v

    # Local aliases for parameters use 

    # Initiate the cache data frame for storing simulation results
    simulate_results = pd.DataFrame([state_zero])
    columns = ["S_su", "S_aa", "S_fa", "S_va", "S_bu", "S_pro", "S_ac", "S_h2", "S_ch4", "S_IC", "S_IN", "S_I", "X_xc1", "X_ch1", "X_pr1", "X_li1", "X_xc2", "X_ch2", "X_pr2", "X_li2", "X_su", "X_aa", "X_fa", "X_c4", "X_pro", "X_ac", "X_h2", "X_I", "S_cation", "S_anion", "pH", "S_va_ion", "S_bu_ion", "S_pro_ion", "S_ac_ion", "S_hco3_ion", "S_co2", "S_nh3", "S_nh4_ion", "S_gas_h2", "S_gas_ch4", "S_gas_co2"]
    simulate_results.columns = columns

    # Initiate cache data frame for storing gasflow values
    initflow = {'p_gas_h2': [0],'p_gas_ch4': [0],'p_gas_co2': [0],'p_gas': [0],'q_gas': [0], 'q_ch4': [0], 'total_ch4': [0],'p_gas_ch4/p_gas': [0],'p_gas_co2/p_gas': [0], 'ch4_yield':[0], 'co2_yield':[0], 'cumulative_methane_yield':[0], 'h2_yield':[0]}
    gasflow = pd.DataFrame(initflow)
    total_ch4 = 0

    # Initiate cache data frame for storing inhibition values
    init_inhib_val = 1 if disable_inhibition else 0
    initflow = {
        'I_5': [init_inhib_val],
        'I_6': [init_inhib_val],
        'I_7': [init_inhib_val],
        'I_8': [init_inhib_val],
        'I_10': [init_inhib_val],
        'I_12': [init_inhib_val],
        'I_pH_aa': [init_inhib_val],
        'I_pH_ac': [init_inhib_val],
        'I_pH_h2': [init_inhib_val],
        'I_IN_lim': [init_inhib_val],
        'I_h2_fa': [init_inhib_val],
        'I_h2_c4': [init_inhib_val],
        'I_h2_pro': [init_inhib_val],
        'I_nh3': [init_inhib_val],
    }
    inhibition = pd.DataFrame(initflow)

    # Initiate cache data frame for storing ions values
    initflow = {'S_cation': [0],'S_anion': [0],'S_H_ion': [0],'S_va_ion': [0],'S_bu_ion': [0],'S_pro_ion': [0],'S_ac_ion': [0],'S_hco3_ion': [0],'S_nh4_ion': [0]}
    ions = pd.DataFrame(initflow)

    ##############################
    ##time definition

    timeSteps_sets = {
    "Day(s)": days,
    "15 Minute(s)": days*24*4,
    "Hour(s)": days*24,
}
    
   

    selected_timeSteps = timeSteps_sets[timesteps] #every 1 hour 


    t = np.linspace(0, days, selected_timeSteps) #sequence of timesteps as fractions of days



    solvermethod = 'DOP853'
    # solvermethod = 'BDF'

    # --- Recycle stream integration ---
    # Each step, compute fresh+recycle mixed influent using flow-weighted average
    # Also record the mixed influent so users can inspect its evolution over time
    mixed_influent_records = []

    for u in t[1:]:
        n += 1

        # Effluent state MUST use base (no *_in) names. Influent uses *_in names.
        effluent_state = {
            'S_su': S_su, 'S_aa': S_aa, 'S_fa': S_fa, 'S_va': S_va, 'S_bu': S_bu, 'S_pro': S_pro,
            'S_ac': S_ac, 'S_h2': S_h2, 'S_ch4': S_ch4, 'S_IC': S_IC, 'S_IN': S_IN, 'S_I': S_I,
            'X_xc1': X_xc1, 'X_ch1': X_ch1, 'X_pr1': X_pr1, 'X_li1': X_li1, 'X_xc2': X_xc2,
            'X_ch2': X_ch2, 'X_pr2': X_pr2, 'X_li2': X_li2, 'X_su': X_su, 'X_aa': X_aa,
            'X_fa': X_fa, 'X_c4': X_c4, 'X_pro': X_pro, 'X_ac': X_ac, 'X_h2': X_h2,
            'X_I': X_I, 'S_cation': S_cation, 'S_anion': S_anion
        }

        # Flow-weighted mixing using helper ( (q_in * fresh + q_r * effluent)/q_ad )
        mixed_influent = mix_influent_with_recycle(new_influent, effluent_state, q_in, q_r, q_ad)

        # Record history (add time for traceability)
        rec = {'time': u}
        rec.update(mixed_influent)
        mixed_influent_records.append(rec)


        # Update state_input for this time step
        state_input = [mixed_influent['S_su_in'], mixed_influent['S_aa_in'], mixed_influent['S_fa_in'], mixed_influent['S_va_in'], mixed_influent['S_bu_in'], mixed_influent['S_pro_in'], 
                       mixed_influent['S_ac_in'], mixed_influent['S_h2_in'], mixed_influent['S_ch4_in'], mixed_influent['S_IC_in'], mixed_influent['S_IN_in'], mixed_influent['S_I_in'],
                       mixed_influent['X_xc1_in'], mixed_influent['X_ch1_in'], mixed_influent['X_pr1_in'], mixed_influent['X_li1_in'], mixed_influent['X_xc2_in'], 
                       mixed_influent['X_ch2_in'], mixed_influent['X_pr2_in'], mixed_influent['X_li2_in'], mixed_influent['X_su_in'], mixed_influent['X_aa_in'], 
                       mixed_influent['X_fa_in'], mixed_influent['X_c4_in'], mixed_influent['X_pro_in'], mixed_influent['X_ac_in'], mixed_influent['X_h2_in'], 
                       mixed_influent['X_I_in'], mixed_influent['S_cation_in'], mixed_influent['S_anion_in']]

        # Span for next time step
        tstep = [t0, u]

        # ...existing simulation code...
        # After simulation step, update effluent_state with new effluent values (from output)
        # effluent_state = ... (update with output from simulation)

        # Build current state vector (y0)
        current_state = [S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I,
                        X_xc1, X_ch1, X_pr1, X_li1, X_xc2, X_ch2, X_pr2, X_li2,
                        X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cation, S_anion,
                        S_H_ion, S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_co2, S_nh3, S_nh4_ion,
                        S_gas_h2, S_gas_ch4, S_gas_co2]

        # ODE integration
        sim = simulate(tstep, current_state, state_input, solvermethod,params)

        # Unpack solution arrays
        (sim_S_su, sim_S_aa, sim_S_fa, sim_S_va, sim_S_bu, sim_S_pro, sim_S_ac, sim_S_h2, sim_S_ch4, sim_S_IC, sim_S_IN, sim_S_I,
        sim_X_xc1, sim_X_ch1, sim_X_pr1, sim_X_li1, sim_X_xc2, sim_X_ch2, sim_X_pr2, sim_X_li2,
        sim_X_su, sim_X_aa, sim_X_fa, sim_X_c4, sim_X_pro, sim_X_ac, sim_X_h2, sim_X_I, sim_S_cation, sim_S_anion,
        sim_S_H_ion, sim_S_va_ion, sim_S_bu_ion, sim_S_pro_ion, sim_S_ac_ion, sim_S_hco3_ion, sim_S_co2, sim_S_nh3, sim_S_nh4_ion,
        sim_S_gas_h2, sim_S_gas_ch4, sim_S_gas_co2) = sim

        # Take last values
        S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, \
        X_xc1, X_ch1, X_pr1, X_li1, X_xc2, X_ch2, X_pr2, X_li2, \
        X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cation, S_anion, \
        S_H_ion, S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_co2, S_nh3, S_nh4_ion, \
        S_gas_h2, S_gas_ch4, S_gas_co2 = \
            sim_S_su[-1], sim_S_aa[-1], sim_S_fa[-1], sim_S_va[-1], sim_S_bu[-1], sim_S_pro[-1], sim_S_ac[-1], sim_S_h2[-1], sim_S_ch4[-1], sim_S_IC[-1], sim_S_IN[-1], sim_S_I[-1], \
            sim_X_xc1[-1], sim_X_ch1[-1], sim_X_pr1[-1], sim_X_li1[-1], sim_X_xc2[-1], sim_X_ch2[-1], sim_X_pr2[-1], sim_X_li2[-1], \
            sim_X_su[-1], sim_X_aa[-1], sim_X_fa[-1], sim_X_c4[-1], sim_X_pro[-1], sim_X_ac[-1], sim_X_h2[-1], sim_X_I[-1], sim_S_cation[-1], sim_S_anion[-1], \
            sim_S_H_ion[-1], sim_S_va_ion[-1], sim_S_bu_ion[-1], sim_S_pro_ion[-1], sim_S_ac_ion[-1], sim_S_hco3_ion[-1], sim_S_co2[-1], sim_S_nh3[-1], sim_S_nh4_ion[-1], \
            sim_S_gas_h2[-1], sim_S_gas_ch4[-1], sim_S_gas_co2[-1]

        # Algebraic update (pure DAE) - pass state, receive corrected state & pH
        state_for_dae = [S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I,
                        X_xc1, X_ch1, X_pr1, X_li1, X_xc2, X_ch2, X_pr2, X_li2,
                        X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cation, S_anion,
                        S_H_ion, S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_co2, S_nh3, S_nh4_ion,
                        S_gas_h2, S_gas_ch4, S_gas_co2]
        
        
        new_state, pH_value = DAESolve(state_for_dae,state_input,params)


        # Overwrite updated components from new_state (others unchanged)
        S_h2 = new_state[7]
        S_H_ion = new_state[30]
        S_va_ion = new_state[31]
        S_bu_ion = new_state[32]
        S_pro_ion = new_state[33]
        S_ac_ion = new_state[34]
        S_hco3_ion = new_state[35]
        S_co2 = new_state[36]
        S_nh3 = new_state[37]
        S_nh4_ion = new_state[38]
        # pH_value available if needed for direct storage/inhibition calcs

        prevS_H_ion = S_H_ion

         # Base inhibition factors via shared utility


        inhib = compute_inhibition_factors(
            S_H_ion=S_H_ion,
            S_IN=S_IN,
            S_h2=S_h2,
            S_nh3=S_nh3,
            params=params,
            disable_inhibition=disable_inhibition,
        )

        I_pH_aa = inhib['I_pH_aa']
        I_pH_ac = inhib['I_pH_ac']
        I_pH_h2 = inhib['I_pH_h2']
        I_IN_lim = inhib['I_IN_lim']
        I_h2_fa = inhib['I_h2_fa']
        I_h2_c4 = inhib['I_h2_c4']
        I_h2_pro = inhib['I_h2_pro']
        I_nh3 = inhib['I_nh3']

        I_5 = inhib['I_5']
        I_6 = inhib['I_6']
        I_7 = inhib['I_7']
        I_8 = inhib['I_8']
        I_9 = inhib['I_9']
        I_10 = inhib['I_10']
        I_11 = inhib['I_11']
        I_12 = inhib['I_12']

        # Store all computed inhibition factors for consistency and maintainability
        inhibittemp = {'I_5': I_5,'I_6': I_6,'I_7': I_7,'I_8': I_8,'I_9': I_9,'I_10': I_10,'I_11': I_11,'I_12': I_12, 'I_pH_aa': I_pH_aa,'I_pH_ac': I_pH_ac,'I_pH_h2': I_pH_h2, 'I_IN_lim': I_IN_lim,'I_h2_fa': I_h2_fa,'I_h2_c4': I_h2_c4,'I_h2_pro': I_h2_pro,'I_nh3': I_nh3}
        inhibition = pd.concat([inhibition, pd.DataFrame([inhibittemp])], ignore_index=True)

        ################

        S_nh4_ion =  (S_IN - S_nh3)
        S_co2 =  (S_IC - S_hco3_ion)
        #pH = - np.log10(S_H_ion)

        # Algebraic equations 

        p_gas_h2 =  (S_gas_h2 * R * T_op / 16)
        p_gas_ch4 =  (S_gas_ch4 * R * T_op / 64)
        p_gas_co2 =  (S_gas_co2 * R * T_op)
        
        Rho_T_8 =  (k_L_a * (S_h2 - 16 * K_H_h2 * p_gas_h2))
        Rho_T_9 =  (k_L_a * (S_ch4 - 64 * K_H_ch4 * p_gas_ch4))
        Rho_T_10 =  (k_L_a * (S_co2 - K_H_co2 * p_gas_co2))
        
        p_gas=  (p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o)
        q_gas =  (k_p * (p_gas- p_atm))
        
        #q_gas= (R*T_op/(p_atm-p_gas_h2o))*V_liq*(Rho_T_8/16+Rho_T_9/64+Rho_T_10)
        
        if q_gas < 0:    
            q_gas = 0
        
        q_ch4 = q_gas * (p_gas_ch4/p_gas) # methane flow
        q_co2 = q_gas * (p_gas_co2/p_gas) # co2 flow
        q_h2 = q_gas * (p_gas_h2/p_gas) # h2 flow

        if q_ch4 < 0:
            q_ch4 = q_gas * (p_gas_co2/p_gas) # co2 flow
        q_h2 = q_gas * (p_gas_h2/p_gas) # h2 flow

        if q_ch4 < 0:
            q_ch4 = 0

        total_ch4 = total_ch4 + q_ch4 

        flowtemp = {'p_gas_h2': p_gas_h2,'p_gas_ch4': p_gas_ch4,'p_gas_co2': p_gas_co2,'p_gas': p_gas,'q_gas': q_gas, 'q_ch4': q_ch4, 'total_ch4': total_ch4, 'p_gas_ch4/p_gas': p_gas_ch4/p_gas,'p_gas_co2/p_gas': p_gas_co2/p_gas, 'ch4_yield':q_ch4/VS_in,  'co2_yield':q_co2/VS_in, 'cumulative_methane_yield': total_ch4 / VS_in, 'h2_yield':q_h2/VS_in}
        gasflow = pd.concat([gasflow, pd.DataFrame([flowtemp])], ignore_index=True)

        total_ch4 = total_ch4 + q_ch4     
        
        ############
        
        ionstemp = {'S_cation': S_cation,'S_anion': S_anion,'S_H_ion': S_H_ion,'S_va_ion': S_va_ion,'S_bu_ion': S_bu_ion,'S_pro_ion': S_pro_ion,'S_ac_ion': S_ac_ion,'S_hco3_ion': S_hco3_ion,'S_nh4_ion': S_nh4_ion}
        ions = pd.concat([ions, pd.DataFrame([ionstemp])], ignore_index=True)

        ##############

        # Rebuild and append state (store pH later)
        state_zero = [S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, \
                    X_xc1, X_ch1, X_pr1, X_li1, X_xc2, X_ch2, X_pr2, X_li2, \
                    X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cation, S_anion, \
                    S_H_ion, S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_co2, S_nh3, S_nh4_ion, \
                    S_gas_h2, S_gas_ch4, S_gas_co2]

        dfstate_zero = pd.DataFrame([state_zero], columns=columns)
        simulate_results = pd.concat([simulate_results, dfstate_zero], ignore_index=True)
        print(u)
        t0 = u 

    # End of time loop
    ##############################

    p_gas_h2 =  (S_gas_h2 * R * T_op / 16)
    p_gas_ch4 =  (S_gas_ch4 * R * T_op / 64)
    p_gas_co2 =  (S_gas_co2 * R * T_op)
    Rho_T_8 =  (k_L_a * (S_h2 - 16 * K_H_h2 * p_gas_h2))
    Rho_T_9 =  (k_L_a * (S_ch4 - 64 * K_H_ch4 * p_gas_ch4))
    Rho_T_10 =  (k_L_a * (S_co2 - K_H_co2 * p_gas_co2))
                
    p_gas=  (p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o)
    q_gas =  (k_p * (p_gas- p_atm))

    #q_gas= (R*T_op/(p_atm-p_gas_h2o))*V_liq*(Rho_T_8/16+Rho_T_9/64+Rho_T_10)

    if q_gas < 0:    
        q_gas = 0

    q_ch4 = q_gas * (p_gas_ch4/p_gas) # methane flow
    if q_ch4 < 0:
        q_ch4 = 0

    phlogarray = -1 * np.log10(simulate_results['pH'])
    simulate_results['pH'] = phlogarray
            

    VS_out=get_VSS(simulate_results.iloc[-1],q_out)

    VS_reduction=(VS_in-VS_out)*100/VS_in

    return {
        "new_influent": new_influent,
        "VS_reduction": VS_reduction,
        "VS_out": VS_out,
        "simulate_results": simulate_results,
        "gasflow": gasflow,
        "inhibition": inhibition,
        "mixed_influent_history": pd.DataFrame(mixed_influent_records),
        "u": t,
        "S_su": S_su,
        "S_aa": S_aa,
        "S_fa": S_fa,
        "S_va": S_va,
        "S_bu": S_bu,
        "S_pro": S_pro,
        "S_ac": S_ac,
        "S_h2": S_h2,
        "S_ch4": S_ch4,
        "S_IC": S_IC,
        "S_IN": S_IN,
        "S_I": S_I,
        "X_xc1": X_xc1,
        "X_ch1": X_ch1,
        "X_pr1": X_pr1,
        "X_li1": X_li1,
        "X_xc2": X_xc2,
        "X_ch2": X_ch2,
        "X_pr2": X_pr2,
        "X_li2": X_li2,
        "X_su": X_su,
        "X_aa": X_aa,
        "X_fa": X_fa,
        "X_c4": X_c4,
        "X_pro": X_pro,
        "X_ac": X_ac,
        "X_h2": X_h2,
        "X_I": X_I,
        "S_cation": S_cation,
        "S_anion": S_anion,
        "S_H_ion": S_H_ion,
        "S_va_ion": S_va_ion,
        "S_bu_ion": S_bu_ion,
        "S_pro_ion": S_pro_ion,
        "S_ac_ion": S_ac_ion,
        "S_hco3_ion": S_hco3_ion,
        "S_co2": S_co2,
        "S_nh3": S_nh3,
        "S_nh4_ion": S_nh4_ion,
        "S_gas_h2": S_gas_h2,
        "S_gas_ch4": S_gas_ch4,
        "S_gas_co2": S_gas_co2,
        "p_gas_h2": p_gas_h2,
        "p_gas_ch4": p_gas_ch4,
        "p_gas_co2": p_gas_co2,
        "p_gas": p_gas,
        "q_gas": q_gas,
        "q_ch4": q_ch4,
        "cumulative_methane_yield": total_ch4 / VS_in,
        "biomethane_yield": q_ch4 / VS_in,
        "biomethane_yield2": q_ch4 / V_ad,
        "q_in": q_in,
        "q_out": q_out,
        "q_recycle": q_r,
        "HRT": HRT,
        'q_in1': q_in1,
        'q_in2': q_in2,
        'q_ad': q_ad,
        'VS_in': VS_in,
        'V_liq': V_liq,
        'V_gas': V_gas,
        'V_ad': V_ad,
        'TS': TS,
        'VSS': VSS,
        'TS_fraction_initial': TS_fraction_initial
    }