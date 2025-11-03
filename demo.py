#!/usr/bin/env python3
"""
ADM1 Bioprocess Simulation Demo

This script demonstrates the key capabilities of the ADM1 simulation platform
for bioprocess engineering and renewable energy applications.

Author: Erfan Hosseini
Date: November 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Add the adm1 module to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from adm1.coAD import ADM1_coAD
    from adm1.params import get_default_params
    from plot_utils import plot_methane_production, plot_substrate_concentrations
    from scenario_compare_utils import compare_scenarios
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure all required modules are available")
    sys.exit(1)


def demo_basic_simulation():
    """Demonstrate basic ADM1 simulation capabilities."""
    print("=" * 60)
    print("DEMO 1: Basic ADM1 Simulation")
    print("=" * 60)
    
    print("Running ADM1 simulation for co-digestion scenario...")
    print("- Feedstock: 40% Food Waste + 60% PPMS")
    print("- Duration: 100 days")
    print("- Reactor volume: 3400 m³")
    
    # Run simulation
    result = ADM1_coAD(
        mixing_ratio=0.4,
        simulation_days=100,
        output_interval=1.0
    )
    
    # Extract key results
    final_methane = result['methane_production'][-1]
    avg_ph = np.mean(result['pH'][50:])  # Steady-state pH
    max_gas_flow = np.max(result['q_gas'])
    
    print(f"\nResults:")
    print(f"- Total methane production: {final_methane:.1f} m³")
    print(f"- Average steady-state pH: {avg_ph:.2f}")
    print(f"- Peak gas flow rate: {max_gas_flow:.1f} m³/day")
    print(f"- Simulation completed successfully!")
    
    return result


def demo_scenario_comparison():
    """Demonstrate multi-scenario analysis."""
    print("\n" + "=" * 60)
    print("DEMO 2: Scenario Comparison Analysis") 
    print("=" * 60)
    
    print("Comparing different feedstock mixing ratios...")
    
    # Define scenarios
    scenarios = {
        'High Food Waste': {'mixing_ratio': 0.8, 'description': '80% FW + 20% PPMS'},
        'Balanced Mix': {'mixing_ratio': 0.5, 'description': '50% FW + 50% PPMS'},
        'Low Food Waste': {'mixing_ratio': 0.2, 'description': '20% FW + 80% PPMS'},
        'Pure PPMS': {'mixing_ratio': 0.0, 'description': '100% PPMS'}
    }
    
    results = {}
    
    print("\nRunning scenarios:")
    for name, config in scenarios.items():
        print(f"- {name}: {config['description']}")
        
        result = ADM1_coAD(
            mixing_ratio=config['mixing_ratio'],
            simulation_days=100,
            output_interval=2.0
        )
        
        results[name] = {
            'methane_yield': result['methane_production'][-1],
            'avg_ph': np.mean(result['pH'][25:]),
            'max_gas_rate': np.max(result['q_gas']),
            'stability_index': np.std(result['pH'][50:])
        }
    
    # Create comparison table
    df = pd.DataFrame(results).T
    df.index.name = 'Scenario'
    
    print(f"\nComparison Results:")
    print(df.round(2))
    
    # Find optimal scenario
    best_scenario = df['methane_yield'].idxmax()
    print(f"\nOptimal scenario for methane production: {best_scenario}")
    print(f"Methane yield: {df.loc[best_scenario, 'methane_yield']:.1f} m³")
    
    return df


def demo_parameter_sensitivity():
    """Demonstrate parameter sensitivity analysis."""
    print("\n" + "=" * 60)
    print("DEMO 3: Parameter Sensitivity Analysis")
    print("=" * 60)
    
    print("Analyzing sensitivity to temperature variation...")
    
    # Temperature range (15°C to 45°C)
    temp_celsius = np.arange(15, 46, 5)
    temp_kelvin = temp_celsius + 273.15
    
    methane_yields = []
    ph_values = []
    
    print("\nTemperature analysis:")
    for temp_c, temp_k in zip(temp_celsius, temp_kelvin):
        print(f"- Running simulation at {temp_c}°C...")
        
        # Get default parameters and modify temperature
        params = get_default_params()
        params['T_op'] = temp_k
        
        # Adjust kinetic rates for temperature (simplified Arrhenius)
        temp_factor = np.exp(0.05 * (temp_c - 35))  # Simplified temperature correction
        params['k_m_ac'] *= temp_factor
        params['k_m_h2'] *= temp_factor
        
        result = ADM1_coAD(
            mixing_ratio=0.5,
            simulation_days=50,
            **params
        )
        
        methane_yields.append(result['methane_production'][-1])
        ph_values.append(np.mean(result['pH'][25:]))
    
    # Find optimal temperature
    optimal_idx = np.argmax(methane_yields)
    optimal_temp = temp_celsius[optimal_idx]
    
    print(f"\nTemperature Sensitivity Results:")
    print(f"- Temperature range: {temp_celsius[0]}°C to {temp_celsius[-1]}°C")
    print(f"- Optimal temperature: {optimal_temp}°C")
    print(f"- Max methane yield: {methane_yields[optimal_idx]:.1f} m³")
    print(f"- Yield variation: {np.std(methane_yields):.1f} m³ (±{100*np.std(methane_yields)/np.mean(methane_yields):.1f}%)")
    
    return temp_celsius, methane_yields, ph_values


def demo_process_optimization():
    """Demonstrate process optimization capabilities."""
    print("\n" + "=" * 60) 
    print("DEMO 4: Process Optimization")
    print("=" * 60)
    
    print("Optimizing feedstock mixing ratio for maximum methane yield...")
    
    # Create fine grid of mixing ratios
    mixing_ratios = np.linspace(0.0, 1.0, 21)
    yields = []
    stabilities = []
    
    print(f"\nTesting {len(mixing_ratios)} mixing ratios...")
    
    for i, ratio in enumerate(mixing_ratios):
        if i % 5 == 0:  # Progress indicator
            print(f"Progress: {100*i/len(mixing_ratios):.0f}%")
        
        result = ADM1_coAD(
            mixing_ratio=ratio,
            simulation_days=80,
            output_interval=2.0
        )
        
        yields.append(result['methane_production'][-1])
        
        # Calculate stability metric (inverse of pH standard deviation)
        ph_stability = 1.0 / (np.std(result['pH'][20:]) + 0.001)
        stabilities.append(ph_stability)
    
    # Find optimal points
    max_yield_idx = np.argmax(yields)
    max_stability_idx = np.argmax(stabilities)
    
    # Multi-objective optimization (weighted score)
    normalized_yields = np.array(yields) / np.max(yields)
    normalized_stability = np.array(stabilities) / np.max(stabilities)
    composite_score = 0.7 * normalized_yields + 0.3 * normalized_stability
    optimal_idx = np.argmax(composite_score)
    
    print(f"\nOptimization Results:")
    print(f"- Maximum yield: {mixing_ratios[max_yield_idx]:.2f} mixing ratio ({yields[max_yield_idx]:.1f} m³)")
    print(f"- Maximum stability: {mixing_ratios[max_stability_idx]:.2f} mixing ratio")
    print(f"- Optimal balance: {mixing_ratios[optimal_idx]:.2f} mixing ratio")
    print(f"  * Methane yield: {yields[optimal_idx]:.1f} m³")
    print(f"  * Composite score: {composite_score[optimal_idx]:.3f}")
    
    return mixing_ratios, yields, stabilities


def create_summary_report():
    """Create a summary visualization of demo results."""
    print("\n" + "=" * 60)
    print("GENERATING SUMMARY REPORT")
    print("=" * 60)
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('ADM1 Simulation Platform - Demo Results', fontsize=16, fontweight='bold')
    
    # Demo 1: Basic simulation time series
    print("Creating basic simulation plot...")
    result = ADM1_coAD(mixing_ratio=0.4, simulation_days=100)
    ax1.plot(result['time'], result['methane_production'], 'b-', linewidth=2)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Cumulative Methane (m³)')
    ax1.set_title('Demo 1: Methane Production Profile')
    ax1.grid(True, alpha=0.3)
    
    # Demo 2: Scenario comparison
    print("Creating scenario comparison plot...")
    ratios = [0.0, 0.2, 0.5, 0.8]
    scenario_yields = []
    for ratio in ratios:
        result = ADM1_coAD(mixing_ratio=ratio, simulation_days=50)
        scenario_yields.append(result['methane_production'][-1])
    
    ax2.bar(range(len(ratios)), scenario_yields, color=['lightcoral', 'gold', 'lightgreen', 'skyblue'])
    ax2.set_xlabel('Food Waste Fraction')
    ax2.set_ylabel('Total Methane Yield (m³)')
    ax2.set_title('Demo 2: Scenario Comparison')
    ax2.set_xticks(range(len(ratios)))
    ax2.set_xticklabels([f'{r:.1f}' for r in ratios])
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Demo 3: Temperature sensitivity
    print("Creating temperature sensitivity plot...")
    temp_range = np.arange(20, 46, 5)
    temp_yields = []
    for temp in temp_range:
        params = get_default_params()
        params['T_op'] = temp + 273.15
        result = ADM1_coAD(mixing_ratio=0.5, simulation_days=30, **params)
        temp_yields.append(result['methane_production'][-1])
    
    ax3.plot(temp_range, temp_yields, 'ro-', linewidth=2, markersize=8)
    ax3.set_xlabel('Temperature (°C)')
    ax3.set_ylabel('Methane Yield (m³)')
    ax3.set_title('Demo 3: Temperature Sensitivity')
    ax3.grid(True, alpha=0.3)
    
    # Demo 4: Optimization landscape
    print("Creating optimization landscape plot...")
    opt_ratios = np.linspace(0, 1, 11)
    opt_yields = []
    for ratio in opt_ratios:
        result = ADM1_coAD(mixing_ratio=ratio, simulation_days=40)
        opt_yields.append(result['methane_production'][-1])
    
    ax4.plot(opt_ratios, opt_yields, 'g-', linewidth=3)
    ax4.scatter(opt_ratios, opt_yields, color='green', s=50, zorder=5)
    optimal_idx = np.argmax(opt_yields)
    ax4.axvline(opt_ratios[optimal_idx], color='red', linestyle='--', alpha=0.7, label=f'Optimum: {opt_ratios[optimal_idx]:.1f}')
    ax4.set_xlabel('Mixing Ratio (Food Waste Fraction)')
    ax4.set_ylabel('Methane Yield (m³)')
    ax4.set_title('Demo 4: Process Optimization')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the figure
    output_path = 'demo_results_summary.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Summary report saved as: {output_path}")
    
    return fig


def main():
    """Run complete ADM1 simulation platform demo."""
    print("🔬 ADM1 BIOPROCESS SIMULATION PLATFORM DEMO")
    print("=" * 60)
    print("Showcasing advanced bioprocess engineering capabilities")
    print("Author: Erfan Hosseini | November 2025")
    print("=" * 60)
    
    try:
        # Run all demonstrations
        demo_basic_simulation()
        demo_scenario_comparison()
        demo_parameter_sensitivity()
        demo_process_optimization()
        
        # Generate summary report
        create_summary_report()
        
        print("\n" + "=" * 60)
        print("✅ ALL DEMOS COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print("\nKey Capabilities Demonstrated:")
        print("✓ Bioprocess modeling with ADM1")
        print("✓ Multi-scenario comparative analysis")
        print("✓ Parameter sensitivity studies")
        print("✓ Process optimization algorithms")
        print("✓ Advanced data visualization")
        print("✓ Automated reporting capabilities")
        
        print(f"\nThis demo showcases expertise in:")
        print("• Bioprocess Engineering & Biochemical Modeling")
        print("• Numerical Methods & Scientific Computing") 
        print("• Data Science & Advanced Analytics")
        print("• Software Engineering & Clean Code Architecture")
        print("• Renewable Energy & Environmental Applications")
        
        print(f"\n📊 Results exported to: demo_results_summary.png")
        print(f"🔗 Full documentation: docs/technical_overview.md")
        
    except Exception as e:
        print(f"\n❌ Demo failed with error: {e}")
        print("Please check dependencies and module imports")
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)