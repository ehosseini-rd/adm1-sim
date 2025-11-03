# Notebooks Directory

This directory contains Jupyter notebooks for interactive analysis and demonstration of the ADM1 simulation platform.

## Main Notebooks

### scenario_comparison.ipynb
**Purpose**: Interactive scenario analysis and comparison
- Compare different feedstock mixing ratios
- Analyze process performance metrics
- Generate comparative visualizations
- Export results in multiple formats

### adm1_coad_analysis.ipynb  
**Purpose**: Detailed co-digestion analysis
- In-depth study of food waste + PPMS co-digestion
- Parameter sensitivity analysis
- Process optimization studies
- Advanced visualization techniques

### mass_balance.ipynb
**Purpose**: Model validation and mass balance verification
- COD and nitrogen balance verification
- Model accuracy assessment
- Debugging and troubleshooting tools
- Validation against literature data

### scenario_comparison_batch.ipynb
**Purpose**: Batch processing of multiple scenarios
- High-throughput scenario analysis
- Statistical comparison methods
- Automated report generation
- Performance benchmarking

## Archived Notebooks

The `archived/` subdirectory contains:
- Development versions
- Backup copies  
- Experimental notebooks
- Legacy analysis files

These files are preserved for reference but are not part of the main workflow.

## Usage

1. **Start Jupyter**: `jupyter notebook` or `jupyter lab`
2. **Navigate** to the notebooks directory
3. **Open** any notebook of interest
4. **Run cells** sequentially for complete analysis

## Dependencies

All notebooks require the packages listed in `requirements.txt`. Ensure the ADM1 modules are in your Python path:

```python
import sys
sys.path.append('..')  # Add parent directory to path
from adm1.coAD import ADM1_coAD
```

## Output

Notebooks generate:
- Interactive plots and visualizations
- Exported data files (CSV, Excel)
- Summary reports (HTML, Markdown)
- Figures for publication (PNG, PDF)

Results are typically saved to the `exports/` directory.

## Best Practices

- **Run in sequence**: Execute cells from top to bottom
- **Save frequently**: Use Ctrl+S to save progress
- **Clear outputs**: Before committing, clear outputs to reduce file size
- **Document changes**: Add markdown cells to explain modifications
- **Use descriptive names**: When saving custom versions