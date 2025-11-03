# Quick Start Guide

This guide will get you up and running with the ADM1 simulation platform in 5 minutes.

## Installation

```bash
git clone https://github.com/Erfan76hosseini/adm1-sim.git
cd adm1-sim
pip install -r requirements.txt
```

## Run the Demo

```bash
python demo.py
```

This will showcase all key capabilities and generate a summary report.

## Your First Simulation

```python
from adm1.coAD import ADM1_coAD

# Simulate co-digestion: 60% food waste + 40% PPMS
result = ADM1_coAD(mixing_ratio=0.6, simulation_days=100)

# Check results
print(f"Methane yield: {result['methane_production'][-1]:.1f} m³")
print(f"Average pH: {result['pH'].mean():.2f}")
```

## Compare Scenarios

```python
from scenario_compare_utils import compare_scenarios

scenarios = {
    'high_fw': {'mixing_ratio': 0.8},
    'balanced': {'mixing_ratio': 0.5},
    'low_fw': {'mixing_ratio': 0.2}
}

results = compare_scenarios(scenarios, duration=365)
print(results)
```

## Next Steps

1. Explore `notebooks/scenario_comparison.ipynb`
2. Read `docs/technical_overview.md` for model details
3. Check `docs/api_reference.md` for complete API
4. Modify parameters in `adm1/params.py`

## Key Outputs

- **Methane Production**: Biogas generation profiles
- **pH Dynamics**: Process stability monitoring  
- **Substrate Utilization**: Feedstock conversion efficiency
- **Process Optimization**: Optimal operating conditions