# uHAF: Unified Hierarchical Annotation Framework for Single-cell Data

## Introduction

The surge in single-cell transcriptomics has led to vast datasets but also to inconsistencies in cell type annotations across studies. Different researchers may label the same cell type with varying names or levels of detail—for example, "cardiomyocyte," "CM," or "heart muscle cell." These discrepancies hinder the comparability of cell labels, affecting machine learning training, evaluation, and data integration.

While efforts like the Cell Ontology aim to standardize nomenclature, limitations exist, such as lack of organ-specific details and user-friendly mapping tools. To tackle these challenges, we introduce **uHAF**, a unified hierarchical annotation framework providing organ-specific cell type hierarchies. With accompanying Python packages, uHAF allows users to unify cell type annotations to consistent granularity or navigate between different hierarchical levels.

Additionally, we developed **uHAF-GPT**, leveraging large language models to map various cell type aliases to uHAF nodes, facilitating rapid unification of labels from diverse datasets. This makes cell type labels comparable across studies, enhancing data integration, machine learning applications, and downstream analyses. We also offer evaluation metrics based on uHAF, enabling fair and biologically meaningful assessment of cell type annotations.

## Features

- **Hierarchical Cell Type Trees**: Provides organ-specific hierarchical structures of cell types.
- **Python Integration**: Includes Python classes and functions for easy manipulation and mapping of cell types.
- **uHAF-GPT Integration**: Utilizes large language models to map diverse cell type names to standardized uHAF nodes.
- **Customizable Granularity**: Allows users to unify annotations to specific hierarchical levels or trace back to broader categories.
- **Evaluation Metrics**: Offers a set of metrics for assessing cell type annotations based on hierarchical relationships.
- **User-Friendly**: Designed to be accessible for researchers with varying levels of computational expertise.

## Installation

### Prerequisites

- Python 3.6 or higher
- Required Python packages:
   - `pandas`
   - `numpy`
   - `tqdm`
   - `openpyxl` (for reading Excel files)

### Installation Steps

1. **Clone the Repository**

```bash
git clone https://github.com/yourusername/uHAF.git
```

2. **Navigate to the Directory**

```bash
cd uHAF
```

3. **Install Required Packages**

Install the required Python packages using `pip`:

```bash
pip install pandas numpy tqdm openpyxl
```

## Usage

### Loading uHAF
To start using uHAF, you need to load the uHAF Excel file and initialize the `uHAF` class.
```python
import os
import pandas as pd
from yourmodule import uHAF, build_uhaf
# Specify the version of the uHAF Excel file
uhaf_xlsx_version = 'uHAF_v1'
# Build the uHAF object
uhaf = build_uhaf(uhaf_xlsx_version)
```
### Accessing Organ-Specific Hierarchies
You can generate uHAF for specific organs by specifying the sheet names:
```python
target_sheetnames = ['Heart', 'Lung']
uhaf = build_uhaf(uhaf_xlsx_version, target_sheetnames)
```
### Mapping Cell Types
To map custom cell types to the uHAF framework:
```python
sheet_name = 'Heart'
custom_cell_types = ['CM', 'cardiomyocyte', 'heart muscle cell']
# Generate GPT prompts for mapping
prompts = uhaf.generate_uhaf_GPTs_prompts(sheet_name, custom_cell_types)
print(prompts)
```
### Navigating Hierarchical Levels
To unify cell type annotations to a specific hierarchical level:
```python
query_cell_types = ['Ventricular Cardiomyocyte', 'Atrial Cardiomyocyte']
annotation_level = 2  # Specify the desired level
# Get the mapping to the specified level
annotation_level_map = uhaf.set_annotation_level(query_cell_types, sheet_name, annotation_level)
print(annotation_level_map)
```
### Tracing Cell Lineage
To trace the hierarchical path from a specific cell type back to the root:
```python
cell_type_target = 'Ventricular Cardiomyocyte'
trace = uhaf.track_cell_from_uHAF(sheet_name, cell_type_target)
print(" -> ".join(trace))
```
## Example
Here's a complete example demonstrating how to unify custom cell type annotations:
```python
import os
import pandas as pd
from yourmodule import uHAF, build_uhaf
# Build uHAF for the 'Heart' organ
uhaf = build_uhaf('uHAF_v1', target_sheetnames=['Heart'])
# Custom cell types to map
custom_cell_types = ['CM', 'cardiomyocyte', 'heart muscle cell']
# Generate GPT prompts
prompts = uhaf.generate_uhaf_GPTs_prompts('Heart', custom_cell_types)
print(prompts)
# Set annotation level
annotation_level_map = uhaf.set_annotation_level(custom_cell_types, 'Heart', annotation_level=2)
print(annotation_level_map)
```
## Contributing
We welcome contributions from the community. If you have suggestions for improvements or new features,please open an issue or submit a pull request.
## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
## Acknowledgments
We would like to thank all contributors and the wider scientific community for their support andcollaboration.

## Contact

For questions or inquiries, please contact [your.email@example.com](mailto:your.email@example.com).