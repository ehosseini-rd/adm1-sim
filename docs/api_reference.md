# API Reference

## Core Modules

### adm1.coAD

Main co-digestion simulation module.

#### ADM1_coAD(mixing_ratio, **kwargs)

**Purpose**: Simulate anaerobic co-digestion with specified mixing ratios.

**Parameters**:
- `mixing_ratio` (float): Ratio of food waste to total feedstock (0-1)
- `simulation_days` (int, optional): Duration in days (default: 100)
- `output_interval` (float, optional): Time step for outputs in days (default: 0.1)
- `reactor_volume` (float, optional): Working volume in m³ (default: 3400)

**Returns**:
- `dict`: Simulation results with timestamps and state variables

**Example**:
```python
from adm1.coAD import ADM1_coAD

# Simulate 50% food waste / 50% PPMS mixture
result = ADM1_coAD(
    mixing_ratio=0.5,
    simulation_days=365,
    output_interval=1.0
)
```

### adm1.solver

Numerical integration methods for the ADM1 system.

#### ODESolver

**Methods**:
- `integrate(t_span, y0, method='DOP853')`: Solve ODE system
- `set_tolerances(rtol, atol)`: Configure solver precision
- `validate_solution()`: Check mass balance and stability

### adm1.params

Parameter management and configuration.

#### get_default_params()

**Returns**: Dictionary of default ADM1 parameters

#### update_params(params, **updates)

**Purpose**: Safely update parameter dictionary with validation.

**Parameters**:
- `params` (dict): Base parameter dictionary
- `**updates`: Parameter updates as keyword arguments

### adm1.inhibition

Inhibition kinetics calculations.

#### pH_inhibition(pH, pH_UL, pH_LL)

**Purpose**: Calculate pH inhibition factor.

**Parameters**:
- `pH` (float): Current pH value
- `pH_UL` (float): Upper pH limit for 50% inhibition  
- `pH_LL` (float): Lower pH limit for 50% inhibition

**Returns**:
- `float`: Inhibition factor (0-1)

## Utility Modules

### plot_utils

Visualization functions for ADM1 results.

#### plot_methane_production(results, **kwargs)

**Purpose**: Create methane production time series plot.

**Parameters**:
- `results` (dict): ADM1 simulation results
- `figsize` (tuple, optional): Figure size (default: (10, 6))
- `save_path` (str, optional): Path to save figure

#### plot_substrate_concentrations(results, substrates=None)

**Purpose**: Plot substrate concentration profiles.

**Parameters**:
- `results` (dict): Simulation results
- `substrates` (list, optional): List of substrates to plot

### scenario_compare_utils

Batch analysis and comparison utilities.

#### compare_scenarios(scenarios, **kwargs)

**Purpose**: Run multiple scenarios and compare results.

**Parameters**:
- `scenarios` (dict): Dictionary of scenario configurations
- `metrics` (list, optional): Metrics to calculate for comparison

**Returns**:
- `pandas.DataFrame`: Comparison results with summary statistics

#### export_comparison(results, format='html', output_dir='exports/')

**Purpose**: Export comparison results to file.

**Parameters**:
- `results` (DataFrame): Comparison results  
- `format` (str): Export format ('html', 'csv', 'markdown')
- `output_dir` (str): Output directory path

## Data Structures

### Simulation Results

Standard result dictionary structure:

```python
results = {
    'time': np.array,           # Time points (days)
    'S_su': np.array,          # Monosaccharides (kg COD/m³)
    'S_aa': np.array,          # Amino acids (kg COD/m³)
    'S_fa': np.array,          # Long-chain fatty acids (kg COD/m³)
    'S_va': np.array,          # Valerate (kg COD/m³)
    'S_bu': np.array,          # Butyrate (kg COD/m³) 
    'S_pro': np.array,         # Propionate (kg COD/m³)
    'S_ac': np.array,          # Acetate (kg COD/m³)
    'S_h2': np.array,          # Hydrogen (kg COD/m³)
    'S_ch4': np.array,         # Methane (kg COD/m³)
    'S_IC': np.array,          # Inorganic carbon (kmol C/m³)
    'S_IN': np.array,          # Inorganic nitrogen (kmol N/m³)
    'pH': np.array,            # pH
    'q_gas': np.array,         # Gas flow rate (m³/d)
    'methane_production': np.array,  # Cumulative methane (m³)
    'mass_balance_error': np.array   # COD balance error (%)
}
```

### Parameter Dictionary

Core parameters structure:

```python
params = {
    # Kinetic parameters
    'k_dis': 0.5,              # Disintegration rate (d⁻¹)
    'k_hyd_ch': 10.0,          # Carbohydrate hydrolysis (d⁻¹)
    'k_hyd_pr': 10.0,          # Protein hydrolysis (d⁻¹)
    'k_hyd_li': 10.0,          # Lipid hydrolysis (d⁻¹)
    
    # Uptake rates
    'k_m_su': 30.0,            # Sugar uptake (d⁻¹)
    'k_m_aa': 50.0,            # Amino acid uptake (d⁻¹)
    'k_m_fa': 6.0,             # LCFA uptake (d⁻¹)
    
    # Inhibition constants
    'K_I_nh3': 0.0018,         # Ammonia inhibition (M)
    'pH_UL_aa': 5.5,           # pH upper limit
    'pH_LL_aa': 4.0,           # pH lower limit
    
    # Stoichiometric coefficients
    'Y_su': 0.1,               # Sugar yield
    'Y_aa': 0.08,              # Amino acid yield
    'Y_fa': 0.06,              # LCFA yield
    
    # Physical constants
    'T_op': 308.15,            # Operating temperature (K)
    'R': 0.083145,             # Gas constant (bar·L/(mol·K))
    'V_liq': 3400,             # Liquid volume (m³)
    'V_gas': 300               # Gas volume (m³)
}
```

## Error Handling

### Common Exceptions

- `ParameterError`: Invalid parameter values or ranges
- `ConvergenceError`: Solver failed to converge
- `MassBalanceError`: Significant COD conservation violation
- `StabilityError`: Solution exhibits unphysical oscillations

### Debugging Tools

#### validate_inputs(params, initial_state)

**Purpose**: Check parameter and initial condition validity.

#### diagnose_convergence_failure(solver_result)

**Purpose**: Analyze failed simulation runs.

#### check_mass_balance(results, tolerance=0.001)

**Purpose**: Verify COD and nitrogen conservation.

## Performance Optimization

### Solver Configuration

```python
# High accuracy for publication results
solver_config = {
    'method': 'DOP853',
    'rtol': 1e-8,
    'atol': 1e-10,
    'max_step': 0.01
}

# Fast screening for parameter studies  
fast_config = {
    'method': 'RK45',
    'rtol': 1e-4,
    'atol': 1e-6,
    'max_step': 0.1
}
```

### Memory Management

- Use `output_interval` to control result size
- Enable `sparse_output` for long simulations
- Implement checkpointing for very long runs

## Integration Examples

### Parameter Sensitivity Analysis

```python
import numpy as np
from adm1.coAD import ADM1_coAD

# Define parameter ranges
k_hyd_range = np.linspace(5, 15, 11)
results = {}

for k_hyd in k_hyd_range:
    result = ADM1_coAD(
        mixing_ratio=0.5,
        k_hyd_ch=k_hyd,
        simulation_days=100
    )
    results[k_hyd] = result['methane_production'][-1]
```

### Custom Feedstock Definition

```python
# Define new feedstock composition
custom_feedstock = {
    'X_ch': 0.35,  # Carbohydrates fraction
    'X_pr': 0.25,  # Proteins fraction  
    'X_li': 0.15,  # Lipids fraction
    'X_I': 0.25,   # Inerts fraction
    'COD_total': 150  # Total COD (kg/m³)
}

# Run simulation with custom feedstock
result = ADM1_coAD(
    mixing_ratio=0.0,  # Pure custom feedstock
    custom_feedstock=custom_feedstock
)
```