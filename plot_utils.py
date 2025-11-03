
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Color mapping for substrate and degrader pairs
SUBSTRATE_DEGRADER_COLORS = {
	'S_su': '#1f77b4', 'X_su': '#1f77b4',
	'S_aa': '#ff7f0e', 'X_aa': '#ff7f0e',
	'S_fa': '#2ca02c', 'X_fa': '#2ca02c',
	'S_va': '#d62728', 'X_c4': '#d62728',
	'S_bu': '#9467bd', 'X_c4': '#9467bd',
	'S_pro': '#8c564b', 'X_pro': '#8c564b',
	'S_ac': '#e377c2', 'X_ac': '#e377c2',
	'S_h2': '#7f7f7f', 'X_h2': '#7f7f7f',
}

# Consistent font size
FONT_SIZE = 12


# Biomass and substrate plot
def plot_biomass_and_substrate(
		df,
		mixing_ratio,
		time,
		timesteps):

	fig, axes = plt.subplots(1, 2, figsize=(16, 6))
	# Substrate with full names and custom styles
	axes[0].plot(time, df['S_su'], label='Monosaccharides (S_su)', linestyle='-', linewidth=2)
	axes[0].plot(time, df['S_aa'], label='Amino acids (S_aa)', linestyle='--', linewidth=2)
	axes[0].plot(time, df['S_fa'], label='LCFA (S_fa)', linestyle='-.', linewidth=2)
	axes[0].plot(time, df['S_va'], label='Total valerate (S_va)', linestyle=':', linewidth=2, alpha=0.7)
	axes[0].plot(time, df['S_bu'], label='Total butyrate (S_bu)', linestyle=':', linewidth=2, alpha=0.7)
	axes[0].plot(time, df['S_pro'], label='Total propionate (S_pro)', linestyle=(0, (3, 1, 1, 1)), linewidth=2)
	axes[0].plot(time, df['S_ac'], label='Total acetate (S_ac)', linestyle=(0, (5, 1)), linewidth=2)
	axes[0].set_title('Substrates for mixing_ratio='+ f"{mixing_ratio:.2f}", fontsize=16)
	axes[0].set_xlabel('Time(days)', fontsize=16)
	axes[0].set_ylabel('Concentration (kg COD/m³)', fontsize=16)
	axes[0].legend(loc='upper right', fontsize=10)
	axes[0].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
	axes[0].tick_params(axis='both', labelsize=12)
	# Biomass with full names and custom styles
	axes[1].plot(time, df['X_su'], label='Sugar degraders (X_su)', linestyle='-', linewidth=2)
	axes[1].plot(time, df['X_aa'], label='Amino Acids degraders (X_aa)', linestyle='--', linewidth=2)
	axes[1].plot(time, df['X_fa'], label='LCFA degraders (X_fa)', linestyle='-.', linewidth=2)
	axes[1].plot(time, df['X_c4'], label='Valerate & butyrate degraders (X_c4)', linestyle=':', linewidth=2, alpha=0.7)
	axes[1].plot(time, df['X_pro'], label='Propionate degraders (X_pro)', linestyle=(0, (3, 1, 1, 1)), linewidth=2)
	axes[1].plot(time, df['X_ac'], label='Acetate degraders (X_ac)', linestyle=(0, (5, 1)), linewidth=2)
	axes[1].plot(time, df['X_h2'], label='Hydrogen degraders (X_h2)', linestyle=(0, (1, 5)), linewidth=2)
	axes[1].set_title('Biomass Compounds for mixing_ratio='+ f"{mixing_ratio:.2f}", fontsize=16)
	axes[1].set_xlabel("Time(days)", fontsize=16)
	axes[1].set_ylabel('Concentration (kg COD/m³)', fontsize=16)
	axes[1].legend(loc='upper right', fontsize=10)
	axes[1].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
	axes[1].tick_params(axis='both', labelsize=12)
	plt.tight_layout()
	plt.show()



# Gas and inhibition plot
def plot_gas_and_inhibition(
		gasflow,
		inhibition,
		mixing_ratio,
		time,
		timesteps):
	
	fig, axes = plt.subplots(1, 2, figsize=(16, 6))
	# Gas production yield (side 1)
	if 'ch4_yield' in gasflow:
		axes[0].plot(time, gasflow['ch4_yield'], label='Methane Yield (CH₄)', linewidth=2)
	if 'co2_yield' in gasflow:
		axes[0].plot(time, gasflow['co2_yield'], label='Carbon Dioxide Yield (CO₂)', linewidth=2)
	#if 'h2_yield' in gasflow:
	#	axes[0].plot(gasflow.index, gasflow['h2_yield'], label='Hydrogen Yield (H₂)', linewidth=2)
	axes[0].set_xlabel('Time(days)', fontsize=16)
	axes[0].set_ylabel('Yield (m³/tonne-VS)', fontsize=16)
	axes[0].set_title('Gas Production Yield for mixing_ratio='+ f"{mixing_ratio:.2f}", fontsize=16)
	axes[0].legend(loc='upper right', fontsize=10)
	axes[0].grid(True)
	axes[0].tick_params(axis='both', labelsize=12)
	# Composite inhibition functions (side 2)
	axes[1].plot(time, inhibition['I_pH_h2'], label='I_pH_h2', linewidth=2)
	axes[1].plot(time, inhibition['I_IN_lim'], label='I_IN_lim', linewidth=2)
	axes[1].plot(time, inhibition['I_h2_fa'], label='I_h2_fa', linewidth=2)
	axes[1].plot(time, inhibition['I_h2_c4'], label='I_h2_c4', linewidth=2)
	axes[1].plot(time, inhibition['I_h2_pro'], label='I_h2_pro', linewidth=2)
	axes[1].plot(time, inhibition['I_pH_aa'], label='I_pH_aa', linewidth=2)
	axes[1].plot(time, inhibition['I_pH_ac'], label='I_pH_ac', linewidth=2)
	axes[1].plot(time, inhibition['I_nh3'], label='I_nh3', linewidth=2)
	axes[1].set_title('Composite Inhibition Functions for mixing_ratio='+ f"{mixing_ratio:.2f}", fontsize=16)
	axes[1].set_xlabel('Time(days)', fontsize=16)
	axes[1].set_ylabel('Inhibition Factor (0–1)', fontsize=16)
	axes[1].set_ylim(0, 1.05)
	axes[1].grid(True, linestyle='--', alpha=0.6)
	axes[1].legend(loc='upper right', fontsize=10)
	axes[1].tick_params(axis='both', labelsize=12)
	plt.tight_layout()
	plt.show()

# pH plot
def plot_ph(
		df,
		mixing_ratio,
		time,
		timesteps):
	
	plt.figure(figsize=(8, 6))
	plt.plot(time, df['pH'], label='pH', color='#1f77b4', linewidth=2)
	plt.title('pH for mixing_ratio='+ f"{mixing_ratio:.2f}", fontsize=16)
	plt.xlabel('Time(days)', fontsize=20)
	plt.ylabel('pH', fontsize=20)
	plt.legend(fontsize=12)
	plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.tight_layout()
	plt.show()


# Inorganic nitrogen (S_IN) plot
def plot_s_in(df):
	"""
	Plot soluble inorganic nitrogen (S_IN) over timesteps.

	Parameters
	----------
	df : pandas.DataFrame
	    Time-indexed DataFrame containing an 'S_IN' (or 'S_in') column.
	"""
	plt.figure(figsize=(8, 6))
	col = 'S_IN' if 'S_IN' in df.columns else ('S_in' if 'S_in' in df.columns else None)
	if col is None:
		raise KeyError("DataFrame must contain 'S_IN' or 'S_in' column to plot inorganic nitrogen.")
	plt.plot(df.index, df[col], label=col, color='#17becf', linewidth=2)
	plt.title('S_IN Over Time', fontsize=20)
	plt.xlabel('Time(days)', fontsize=20)
	plt.ylabel('Concentration (kg N/m³)', fontsize=20)
	plt.legend(fontsize=12)
	plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.tight_layout()
	plt.show()





# Gas and inhibition plot
def plot_cumulativegas_and_inhibition(gasflow, inhibition):
	fig, axes = plt.subplots(1, 2, figsize=(16, 6))
	# Gas production yield (side 1)
	if 'cumulative_methane_yield' in gasflow:
		axes[0].plot(gasflow.index, gasflow['cumulative_methane_yield'], label='Methane Yield (CH₄)', linewidth=2)
	#if 'h2_yield' in gasflow:
	#	axes[0].plot(gasflow.index, gasflow['h2_yield'], label='Hydrogen Yield (H₂)', linewidth=2)
	axes[0].set_xlabel('Timesteps (hours)', fontsize=16)
	axes[0].set_ylabel('Yield (m³/tonne-VS)', fontsize=16)
	axes[0].set_title('Gas Production Yield for FW', fontsize=20)
	axes[0].legend(loc='upper right', fontsize=10)
	axes[0].grid(True)
	axes[0].tick_params(axis='both', labelsize=12)
	# Composite inhibition functions (side 2)
	axes[1].plot(inhibition.index, inhibition['I_pH_h2'], label='I_pH_h2', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_IN_lim'], label='I_IN_lim', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_h2_fa'], label='I_h2_fa', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_h2_c4'], label='I_h2_c4', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_h2_pro'], label='I_h2_pro', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_pH_aa'], label='I_pH_aa', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_pH_ac'], label='I_pH_ac', linewidth=2)
	axes[1].plot(inhibition.index, inhibition['I_nh3'], label='I_nh3', linewidth=2)
	axes[1].set_title('Composite Inhibition Functions Over Time', fontsize=16)
	axes[1].set_xlabel('Timesteps (hours)', fontsize=16)
	axes[1].set_ylabel('Inhibition Factor (0–1)', fontsize=16)
	axes[1].set_ylim(0, 1.05)
	axes[1].grid(True, linestyle='--', alpha=0.6)
	axes[1].legend(loc='upper right', fontsize=10)
	axes[1].tick_params(axis='both', labelsize=12)
	plt.tight_layout()
	plt.show()


	

	